#!/usr/bin/env python3
"""Synthetic mechanics tests for S11c_a_cross_engine_comparator.

These fixtures never load or run either measured S11c-a payload.
"""

from __future__ import annotations

import unittest

import sympy as sp
from sympy.core.function import AppliedUndef

import S11c_a_cross_engine_comparator as comparator


class TypedKeyTests(unittest.TestCase):
    def test_fixed_axis_order_is_independent_of_encounter_order(self) -> None:
        first = comparator.make_key(
            (("DOF", "DELTA_W"), ("BRANCH", "LAB_HELD"), ("FACE", "1"))
        )
        second = comparator.make_key(
            (("FACE", "FACE_PLUS"), ("DOF", "DELTA_W"), ("BRANCH", "LAB_HELD"))
        )
        self.assertEqual(first, second)
        self.assertEqual(
            first,
            (
                ("BRANCH", "LAB_HELD"),
                ("FACE", "PLUS"),
                ("DOF", "DELTA_W"),
            ),
        )

    def test_duplicate_axis_is_rejected(self) -> None:
        with self.assertRaises(comparator.AxisError):
            comparator.make_key((("DOF", "DELTA_W"), ("DOF", "ZETA_C")))

    def test_direction_is_one_two_or_three_not_face_sign(self) -> None:
        with self.assertRaises(comparator.AxisError):
            comparator.make_key((("DIRECTION", "-1"),))


class ClosedFoldTests(unittest.TestCase):
    def test_current_derivative_preserves_all_arguments(self) -> None:
        key = comparator.make_key((("BRANCH", "LAB_HELD"),))
        raw = (
            "Derivative[0, 0, 1, {0, 0}][currentXPerturbation3]"
            "[x1, x2, x3, {normalCoordinate, time}]"
        )
        value = comparator.parse_wl_value(raw, key)
        self.assertIsInstance(value, AppliedUndef)
        self.assertEqual(value.func.__name__, "delta_j_bulk_3_d3")
        self.assertEqual(len(value.args), 4)
        self.assertEqual(value.args[-1], (sp.Symbol("w"), sp.Symbol("time")))
        self.assertNotEqual(value, sp.Symbol("delta_j_bulk_3_d3"))

    def test_current_base_head_preserves_argument_arity(self) -> None:
        key = comparator.make_key((("BRANCH", "LAB_HELD"),))
        value = comparator.parse_wl_value(
            "currentWPerturbation[x1,x2,x3,{normalCoordinate,time}]", key
        )
        self.assertIsInstance(value, AppliedUndef)
        self.assertEqual(value.func.__name__, "delta_j_bulk_4")
        self.assertEqual(len(value.args), 4)

    def test_mu_registry_is_branch_specific_and_arg_preserving(self) -> None:
        lab = comparator.make_key((("BRANCH", "LAB_HELD"),))
        material = comparator.make_key((("BRANCH", "MATERIAL_ADVECTED"),))
        lab_value = comparator.parse_wl_value("muThetaOperand[x1,x2,x3,time]", lab)
        material_value = comparator.parse_wl_value("muThetaOperand[x1,x2,x3,time]", material)
        self.assertEqual(lab_value.func.__name__, "mu_theta_L")
        self.assertEqual(material_value.func.__name__, "mu_theta_M")
        self.assertEqual(lab_value.args, material_value.args)
        self.assertNotEqual(lab_value, sp.Symbol("mu_theta_L"))

    def test_inactive_equal_is_a_relation_not_native_boolean(self) -> None:
        key = comparator.make_key(())
        value = comparator.parse_wl_value("Inactive[Equal][left,right]", key)
        self.assertIsInstance(value, sp.Equality)
        self.assertEqual(value.lhs, sp.Symbol("left"))
        self.assertEqual(value.rhs, sp.Symbol("right"))

    def test_window_and_inactive_integral_bridge(self) -> None:
        key = comparator.make_key(())
        value = comparator.parse_wl_value(
            "Inactive[Integrate][Derivative[1,0][windowFunction][normalCoordinate-a,-normalCoordinate-b],"
            "{normalCoordinate,-Infinity,Infinity}]",
            key,
        )
        integrals = list(value.atoms(AppliedUndef))
        self.assertTrue(any(item.func == comparator.BOUND_INTEGRAL for item in integrals))
        self.assertTrue(any(item.func == comparator.Dwin for item in integrals))
        bound = next(item for item in integrals if item.func == comparator.BOUND_INTEGRAL)
        self.assertEqual(bound.args[1:3], (-sp.oo, sp.oo))


class BoundIntegralTests(unittest.TestCase):
    def bound_nodes(self, value: sp.Basic) -> list[AppliedUndef]:
        return [
            item
            for item in value.atoms(AppliedUndef)
            if item.func == comparator.BOUND_INTEGRAL
        ]

    def test_bounds_are_retained_and_different_bounds_do_not_combine(self) -> None:
        w = sp.Symbol("w")
        f = sp.Function("f")
        value = comparator.canonical_integrals(
            sp.Integral(f(w), (w, 0, 1)) + sp.Integral(f(w), (w, 0, 2))
        )
        bounds = {(item.args[1], item.args[2]) for item in self.bound_nodes(value)}
        self.assertEqual(bounds, {(sp.Integer(0), sp.Integer(1)), (sp.Integer(0), sp.Integer(2))})

    def test_same_bounds_combine_integrands(self) -> None:
        w = sp.Symbol("w")
        f = sp.Function("f")
        g = sp.Function("g")
        value = comparator.canonical_integrals(
            sp.Integral(f(w), (w, -1, 1)) + sp.Integral(g(w), (w, -1, 1))
        )
        nodes = self.bound_nodes(value)
        self.assertEqual(len(nodes), 1)
        self.assertEqual(nodes[0].args[3], f(comparator.BOUND_BINDER) + g(comparator.BOUND_BINDER))

    def test_binder_dependent_factor_stays_inside(self) -> None:
        w = sp.Symbol("w")
        c = sp.Function("c")
        f = sp.Function("f")
        value = comparator.canonical_integrals(sp.Integral(c(w) * f(w), (w, 0, 1)))
        node = self.bound_nodes(value)[0]
        self.assertEqual(
            node.args[3], c(comparator.BOUND_BINDER) * f(comparator.BOUND_BINDER)
        )

    def test_free_factor_is_pulled_out(self) -> None:
        w = sp.Symbol("w")
        c = sp.Symbol("c")
        f = sp.Function("f")
        value = comparator.canonical_integrals(sp.Integral(c * f(w), (w, 0, 1)))
        node = self.bound_nodes(value)[0]
        self.assertEqual(sp.simplify(value / node), c)
        self.assertEqual(node.args[3], f(comparator.BOUND_BINDER))

    def test_substitution_does_not_capture_existing_common_binder(self) -> None:
        w = sp.Symbol("w")
        f = sp.Function("f")
        value = comparator.canonical_integrals(
            sp.Integral(f(w) + comparator.BOUND_BINDER, (w, 0, 1))
        )
        node = self.bound_nodes(value)[0]
        self.assertTrue(node.args[3].has(comparator.CAPTURED_FREE_BINDER))


class ResidualTests(unittest.TestCase):
    def test_budget_expiry_falls_back_to_exact_closed_cas_route(self) -> None:
        x, y = sp.symbols("x y")
        left = (x + y) ** 8 / (x + 1)
        right = (x - y) ** 8 / (x + 1)
        value = comparator.cached_typed_residual(
            left, right, "SYNTHETIC", 1.0e-9
        )
        expected = sp.cancel(sp.together(sp.expand(left - right)))
        self.assertNotEqual(type(value).__name__, "ResidualFailure")
        self.assertEqual(value, expected)


if __name__ == "__main__":
    unittest.main()
