#!/usr/bin/env python3
"""Synthetic mechanics tests for the S11c-b comparator.

These fixtures do not load or run either measured S11c-b engine.
"""

from __future__ import annotations

import contextlib
import io
import tempfile
import unittest
from pathlib import Path

import sympy as sp
from sympy.core.function import AppliedUndef

import S11c_b_cross_engine_comparator as comparator


class TypedKeyTests(unittest.TestCase):
    def test_fixed_axis_order_and_sector_spelling_fold(self) -> None:
        key = comparator.make_key(
            (
                ("DOF", "THETA"),
                ("SECTOR", "TRANSVERSE_TO_THETA_EW_UL"),
                ("BRANCH", "LAB_HELD"),
            )
        )
        self.assertEqual(
            key,
            (
                ("BRANCH", "LAB_HELD"),
                ("SECTOR", "TRANSVERSE_TO_THICKNESS"),
                ("DOF", "THETA"),
            ),
        )

    def test_duplicate_axis_is_rejected(self) -> None:
        with self.assertRaises(comparator.AxisError):
            comparator.make_key((("SOURCE", "W_BG"), ("SOURCE", "MU_R_BG")))

    def test_direction_and_source_are_typed(self) -> None:
        with self.assertRaises(comparator.AxisError):
            comparator.make_key((("DIRECTION", "0"),))
        with self.assertRaises(comparator.AxisError):
            comparator.make_key((("SOURCE", "BOTH"),))


class LoaderTests(unittest.TestCase):
    def test_wl_payload_is_reassembled_until_the_next_tag(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "wl.out"
            path.write_text(
                "WL_S11CB_FIRST: <|\n\"EXPRESSION\" -> x\n|>\n"
                "WL_S11CB_SECOND: <|\"EXPRESSION\" -> y|>\n",
                encoding="utf-8",
            )
            loaded = comparator.load_wl(path)
        self.assertEqual(set(loaded), {"FIRST", "SECOND"})
        self.assertEqual(comparator.wl_field(loaded["FIRST"], "EXPRESSION"), "x")


class ClosedFoldTests(unittest.TestCase):
    def test_mu_registry_is_branch_specific_and_argument_preserving(self) -> None:
        lab = comparator.make_key((("BRANCH", "LAB_HELD"),))
        material = comparator.make_key((("BRANCH", "MATERIAL_ADVECTED"),))
        lab_value = comparator.parse_wl_value(
            "muThetaOperand[xOne,xTwo,xThree,time]", lab
        )
        material_value = comparator.parse_wl_value(
            "muThetaOperand[xOne,xTwo,xThree,time]", material
        )
        self.assertIsInstance(lab_value, AppliedUndef)
        self.assertEqual(lab_value.func.__name__, "mu_theta_L")
        self.assertEqual(material_value.func.__name__, "mu_theta_M")
        self.assertEqual(lab_value.args, material_value.args)

    def test_mu_heads_are_never_globally_collapsed(self) -> None:
        key = comparator.make_key((("BRANCH", "LAB_HELD"),))
        value = comparator.parse_py_value(
            "Tuple(Symbol('mu_theta_L'), Symbol('mu_theta_M'))",
            key,
            coefficient=False,
        )
        self.assertEqual(value, (sp.Symbol("mu_theta_L"), sp.Symbol("mu_theta_M")))

    def test_s11cb_bookkeepers_are_unit_renamed(self) -> None:
        key = comparator.make_key(())
        value = comparator.parse_wl_value(
            "backgroundOrder + waveOrder + virtualOrder + wave + order", key
        )
        self.assertEqual(value, 5)

    def test_supplied_support_symbols_remain_literal(self) -> None:
        key = comparator.make_key((("BRANCH", "LAB_HELD"),))
        py_value = comparator.parse_py_value(
            "Symbol('t_hold_plus_0_4')", key, coefficient=False
        )
        wl_value = comparator.parse_wl_value(
            "tractionHoldUpperW[xOne,xTwo,xThree]", key
        )
        self.assertEqual(py_value, sp.Symbol("t_hold_plus_0_4"))
        self.assertEqual(py_value, wl_value)

    def test_generated_profile_jet_keeps_all_orders(self) -> None:
        key = comparator.make_key(())
        value = comparator.parse_wl_value("widthProfileJet[1,0,2]", key)
        self.assertEqual(value, sp.Symbol("w1_profile_d1d3d3"))

    def test_energy_control_object_spelling_only(self) -> None:
        py_payload = (
            "Tuple(Tuple(Tuple(Str('ENERGY_BASIS'), Str('LAB_HELD'), "
            "Str('RHO4_CONSTANT'), Str('W_BG'), Integer(1)), "
            "Tuple(Tuple(Str('VALUE'), Symbol('representative_A')))))"
        )
        wl_payload = (
            '<|"ENERGY_BASIS_VARIABLE|LAB_HELD|RHO4_CONSTANT|W_BG|DIRECTION_1" -> '
            '<|"EXPRESSION" -> representativeB|>|>'
        )
        py_cases = comparator.extract_control(
            "PY", "CONTROL_FORM_BASE_OPERAND", py_payload
        )
        wl_cases = comparator.extract_control(
            "WL", "CONTROL_FORM_BASE_OPERAND", wl_payload
        )
        self.assertEqual(py_cases[0].key, wl_cases[0].key)
        self.assertIn("non_unique_energy_quotient", py_cases[0].note or "")


class IntegralTests(unittest.TestCase):
    def test_bound_integral_keeps_binder_and_bounds(self) -> None:
        w = sp.Symbol("w")
        f = sp.Function("f")
        value = comparator.canonical_integrals(sp.Integral(f(w), (w, 0, 2)))
        node = next(
            item
            for item in value.atoms(AppliedUndef)
            if item.func == comparator.BOUND_INTEGRAL
        )
        self.assertEqual(node.args[1:3], (0, 2))


class DivergenceQuotientTests(unittest.TestCase):
    def test_total_derivative_has_zero_spatial_euler_signature(self) -> None:
        theta = sp.Symbol("theta_probe")
        theta_d1 = sp.Symbol("theta_probe_d1")
        width = sp.Symbol("W_bg")
        width_d1 = sp.Symbol("W_bg_d1")
        reduced = comparator.modulo_total_divergence(
            theta_d1 * width + theta * width_d1
        )
        self.assertEqual(reduced.as_dict()["theta_probe"], 0)

    def test_nondivergence_density_is_retained(self) -> None:
        theta = sp.Symbol("theta_probe")
        width = sp.Symbol("W_bg")
        reduced = comparator.modulo_total_divergence(theta * width)
        self.assertEqual(reduced.as_dict()["theta_probe"], width)


class ResidualTests(unittest.TestCase):
    def test_small_rational_leaf_uses_exact_closed_canon(self) -> None:
        x, y = sp.symbols("x y")
        value = comparator.typed_residual(
            (x + y) ** 3 / (x + 1),
            (x - y) ** 3 / (x + 1),
            "SYNTHETIC",
            0.1,
            100,
        )
        expected = sp.cancel(
            sp.together(sp.expand((x + y) ** 3 / (x + 1) - (x - y) ** 3 / (x + 1)))
        )
        self.assertEqual(sp.cancel(sp.together(value - expected)), 0)

    def test_oversized_scalar_is_exact_unevaluated_difference(self) -> None:
        x, y = sp.symbols("x y")
        value = comparator.typed_residual(x, y, "SYNTHETIC", 0.1, 200_001)
        self.assertIsInstance(value, comparator.SymbolicDifference)
        self.assertEqual(value.expression, sp.Add(x, -y, evaluate=False))


class RunSemanticsTests(unittest.TestCase):
    def test_disagreement_is_measurement_exit_zero(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            py_path = root / "py.out"
            wl_path = root / "wl.out"
            py_path.write_text("PY_S11CB_SYNTHETIC: Symbol('left')\n", encoding="utf-8")
            wl_path.write_text("WL_S11CB_SYNTHETIC: right\n", encoding="utf-8")
            stream = io.StringIO()
            with contextlib.redirect_stdout(stream):
                code = comparator.run(
                    ["--py", str(py_path), "--wl", str(wl_path)]
                )
        rendered = stream.getvalue()
        self.assertEqual(code, 0)
        self.assertIn("operand_A =", rendered)
        self.assertIn("operand_B =", rendered)
        self.assertIn("A_minus_B =", rendered)
        self.assertIn("ACCOUNTING SYNTHETIC", rendered)
        self.assertIn("RUN_ACCOUNTING families=1", rendered)


if __name__ == "__main__":
    unittest.main()
