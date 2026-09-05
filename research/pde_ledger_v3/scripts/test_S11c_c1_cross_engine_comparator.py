#!/usr/bin/env python3
"""Synthetic mechanics tests for the S11c-c1 comparator.

These fixtures do not load or run either measured c1 engine.
"""

from __future__ import annotations

import contextlib
import io
import re
import tempfile
import unittest
from pathlib import Path

import sympy as sp
from sympy.core.function import AppliedUndef

import S11c_c1_cross_engine_comparator as comparator


class TypedKeyTests(unittest.TestCase):
    def test_face_and_direction_are_distinct_axes(self) -> None:
        face = comparator.make_key((("FACE", "1"),))
        direction = comparator.make_key((("DIRECTION", "1"),))
        self.assertNotEqual(face, direction)
        self.assertEqual(comparator.make_key((("FACE", "-1"),)), (("FACE", "-1"),))
        with self.assertRaises(comparator.AxisError):
            comparator.make_key((("DIRECTION", "-1"),))
        self.assertEqual(
            comparator.decode_py_key("Tuple(Integer(-1))", ("FACE",)),
            (("FACE", "-1"),),
        )
        with self.assertRaises(comparator.AxisError):
            comparator.decode_py_key("Tuple(Integer(-1))", ("DIRECTION",))

    def test_duplicate_and_unknown_axes_are_rejected(self) -> None:
        with self.assertRaises(comparator.AxisError):
            comparator.make_key((("FACE", "1"), ("FACE", "-1")))
        with self.assertRaises(comparator.AxisError):
            comparator.make_key((("UNTYPED", "x"),))

    def test_wl_prefixes_are_typed_without_cross_axis_folds(self) -> None:
        self.assertEqual(
            comparator.decode_wl_key("FACE_-1", ("FACE",)), (("FACE", "-1"),)
        )
        self.assertEqual(
            comparator.decode_wl_key("DIRECTION_1", ("DIRECTION",)),
            (("DIRECTION", "1"),),
        )
        self.assertEqual(
            comparator.decode_wl_key("PARITY_DELTA_W", ("PARITY",)),
            (("PARITY", "DELTA_W"),),
        )
        with self.assertRaises(comparator.AxisError):
            comparator.decode_wl_key("1", ("FACE",))

    def test_regime_and_parity_literals_do_not_prejoin(self) -> None:
        py_regime = comparator.make_key((("REGIME_OUT", "PROPAGATING"),))
        wl_regime = comparator.make_key((("REGIME_OUT", "OUTPUT_PROPAGATING"),))
        py_parity = comparator.make_key((("PARITY", "THICKNESS"),))
        wl_parity = comparator.make_key((("PARITY", "DELTA_W"),))
        self.assertNotEqual(py_regime, wl_regime)
        self.assertNotEqual(py_parity, wl_parity)

    def test_equal_integer_on_different_axes_does_not_join(self) -> None:
        py = comparator.py_case(
            comparator.make_key((("OBJECT", "DTN"), ("FACE", "1"))), "Symbol('x')"
        )
        wl = comparator.wl_case(
            comparator.make_key((("OBJECT", "DTN"), ("DIRECTION", "1"))), "x"
        )
        stream = io.StringIO()
        with contextlib.redirect_stdout(stream):
            accounting = comparator.compare_family(
                "SYNTHETIC", comparator.FamilyCases([py], [wl]), budget=0.1
            )
        self.assertEqual(accounting.join, 0)
        self.assertEqual(accounting.py_only, 1)
        self.assertEqual(accounting.wl_only, 1)
        self.assertEqual(accounting.axis_set_mismatch, 2)
        self.assertIn("axis_set_mismatch", stream.getvalue())

    def test_closed_axes_and_regime_roles_are_enforced(self) -> None:
        with self.assertRaises(comparator.AxisError):
            comparator.make_key((("SCENARIO", "BOGUS"),))
        with self.assertRaises(comparator.AxisError):
            comparator.make_key((("LIMIT", "BOGUS"),))
        with self.assertRaises(comparator.AxisError):
            comparator.make_key((("REGIME_OUT", "INPUT_PROPAGATING"),))
        with self.assertRaises(comparator.AxisError):
            comparator.make_key((("REGIME_IN", "OUTPUT_PROPAGATING"),))


class LoaderTests(unittest.TestCase):
    def test_wl_multiline_association_is_reassembled(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "wl.out"
            path.write_text(
                'WL_S11CC1_FIRST: <|\n"EXPRESSION" -> x\n|>\n'
                'WL_S11CC1_SECOND: <|"EXPRESSION" -> y|>\n',
                encoding="utf-8",
            )
            loaded = comparator.load_wl(path)
        self.assertEqual(set(loaded), {"FIRST", "SECOND"})
        self.assertEqual(comparator.wl_field(loaded["FIRST"], "EXPRESSION"), "x")

    def test_duplicate_tag_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "py.out"
            path.write_text(
                "PY_S11CC1_SYNTHETIC: Symbol('x')\n"
                "PY_S11CC1_SYNTHETIC: Symbol('y')\n",
                encoding="utf-8",
            )
            with self.assertRaises(comparator.InputError):
                comparator.load_py(path)


class ExtractionCoverageTests(unittest.TestCase):
    def test_dimensions_keeps_ratio_definitions_and_derivations(self) -> None:
        py_payload = (
            "Tuple(Tuple(Str('DTN'), "
            "ImmutableDenseMatrix([[Integer(1)], [Integer(2)], [Integer(3)]])))"
        )
        wl_payload = (
            '<|"OBJECT_DIMENSIONS" -> <|"DTN" -> {1,2,3}|>, '
            '"DERIVATION_EQUATIONS" -> <|"DTN" -> Inactive[Equal][x,x]|>, '
            '"RATIO_DEFINITIONS" -> <|"DTN" -> p/v|>|>'
        )
        cases = comparator.extract_family("DIMENSIONS", py_payload, wl_payload)
        self.assertEqual(len(cases.py), 1)
        self.assertEqual(len(cases.wl), 3)
        self.assertIn(
            "RATIO_DEFINITION", {dict(case.key)["LEAF"] for case in cases.wl}
        )


class ClosedSpellingTests(unittest.TestCase):
    def test_mechanical_map_injectivity_is_checked_first(self) -> None:
        with self.assertRaises(comparator.InputError):
            comparator.checked_mechanical_symbol_map({"alpha_beta", "alphaBeta"})
        inverse = comparator.checked_mechanical_symbol_map({"rho_m", "Lambda_A_0"})
        self.assertEqual(inverse["rhoM"], "rho_m")
        self.assertEqual(inverse["LambdaA0"], "Lambda_A_0")
        with self.assertRaisesRegex(comparator.InputError, "non-injective"):
            comparator.verify_active_spelling_map(
                ("Symbol('alpha_beta') + Symbol('alphaBeta')",)
            )

    def test_mu_theta_representatives_are_never_registry_folded(self) -> None:
        key = comparator.make_key((("BRANCH", "LAB_HELD"),))
        py_names = (
            "s11cc1_mu_theta_lab_held_plus",
            "s11cc1_mu_theta_material_advected_minus",
            "mu_theta_drive",
        )
        for name in py_names:
            self.assertEqual(
                comparator.parse_py_value(f"Symbol('{name}')", key, coefficient=False),
                sp.Symbol(name),
            )
        wl_value = comparator.parse_wl_value("muTheta", key)
        self.assertEqual(wl_value, sp.Symbol("muTheta"))
        difference = comparator.typed_residual(
            sp.Symbol(py_names[0]), wl_value, "SYNTHETIC", 0.1, 20
        )
        self.assertNotEqual(difference, 0)

        applied = comparator.parse_wl_value("muThetaOperand[xOne]", key)
        self.assertIsInstance(applied, AppliedUndef)
        self.assertEqual(applied.func.__name__, "muThetaOperand")
        self.assertEqual(applied.args, (sp.Symbol("xOne"),))

    def test_applied_heads_keep_their_arguments(self) -> None:
        key = comparator.make_key(())
        for raw, head, arity in (
            ("qOut[omega,{kOne,kTwo,kThree}]", "qOut", 2),
            ("w1Profile[xOne/LW,xTwo/LW,xThree/LW]", "w1Profile", 3),
            (
                "rhoBrBgRho4Constant[xOne,xTwo,xThree]",
                "rhoBrBgRho4Constant",
                3,
            ),
        ):
            value = comparator.parse_wl_value(raw, key)
            self.assertIsInstance(value, AppliedUndef)
            self.assertEqual(value.func.__name__, head)
            self.assertEqual(len(value.args), arity)

    def test_inactive_greater_and_fourier_transform_stay_held(self) -> None:
        key = comparator.make_key(())
        greater = comparator.parse_wl_value("Inactive[Greater][x,0]", key)
        transform = comparator.parse_wl_value(
            "Inactive[FourierTransform][f[x],{x},{kOne}]", key
        )
        self.assertIsInstance(greater, AppliedUndef)
        self.assertEqual(greater.func.__name__, "HeldInactiveGreater")
        self.assertIsInstance(transform, AppliedUndef)
        self.assertEqual(transform.func.__name__, "HeldInactiveFourierTransform")

    def test_omega_assumptions_are_not_erased(self) -> None:
        key = comparator.make_key(())
        real_omega = comparator.parse_py_value(
            "Symbol('omega', real=True)", key, coefficient=False
        )
        plain_omega = comparator.parse_py_value(
            "Symbol('omega')", key, coefficient=False
        )
        self.assertTrue(real_omega.is_real)
        self.assertIsNone(plain_omega.is_real)
        self.assertNotEqual(sp.srepr(real_omega), sp.srepr(plain_omega))

    def test_sound_speed_spelling_is_mechanical_and_injective(self) -> None:
        self.assertEqual(comparator.mechanical_lower_camel("c_s0"), "cS0")
        names = comparator.REQUIRED_C1_PY_SYMBOLS
        self.assertEqual(len(names), 12)
        inverse = comparator.checked_mechanical_symbol_map(names)
        self.assertEqual(len(inverse), 12)
        self.assertEqual(inverse, comparator.C1_BARE_SYMBOL)
        active, present = comparator.verify_active_spelling_map(
            (sp.srepr(sp.Tuple(*(sp.Symbol(name) for name in sorted(names)))),)
        )
        self.assertEqual(present, names)
        self.assertEqual(active, inverse)

    def test_multi_range_integral_and_nested_limit_stay_held(self) -> None:
        coordinates = sp.symbols("xOne xTwo xThree")
        raw = (
            "Inactive[Integrate][f[xOne,xTwo,xThree],"
            "{xOne,-Infinity,Infinity},{xTwo,-Infinity,Infinity},"
            "{xThree,-Infinity,Infinity}]"
        )
        integral = sp.Function("HeldInactiveIntegrate")(
            sp.Function("f")(*coordinates),
            *(sp.Tuple(coordinate, -sp.oo, sp.oo) for coordinate in coordinates),
        )
        limit = sp.Function("HeldInactiveLimit")(
            integral, sp.Function("Rule")(sp.Symbol("outwardDistance"), sp.oo)
        )
        for operand, expected in (
            (raw, integral),
            (f"Inactive[Limit][{raw}, outwardDistance -> Infinity]", limit),
        ):
            with self.subTest(operand=operand):
                value = comparator.parse_wl_value(operand, ())
                self.assertIsInstance(value, AppliedUndef)
                self.assertEqual(value, expected)
                self.assertFalse(value.has(sp.Integral, sp.Limit))
                stream = io.StringIO()
                with contextlib.redirect_stdout(stream):
                    accounting = comparator.compare_family(
                        "SYNTHETIC", comparator.FamilyCases(
                            wl=[comparator.wl_case((), operand)]
                        ), budget=0.1,
                    )
                self.assertEqual(accounting.join, 0)
                self.assertEqual(accounting.wl_only, 1)
                self.assertEqual(accounting.parse_failed, 0)
                self.assertIn(f"operand_B = {sp.srepr(expected)}", stream.getvalue())
                self.assertNotIn("<PARSE_FAILED>", stream.getvalue())


class ResidualTests(unittest.TestCase):
    def test_bound_integral_keeps_binder_and_bounds(self) -> None:
        key = comparator.make_key(())
        value = comparator.parse_wl_value(
            "Inactive[Integrate][f[w],{w,0,2}]", key
        )
        nodes = [
            item
            for item in value.atoms(AppliedUndef)
            if item.func == comparator.BOUND_INTEGRAL
        ]
        self.assertEqual(len(nodes), 1)
        self.assertEqual(nodes[0].args[1:3], (0, 2))
        self.assertEqual(
            value,
            comparator.BOUND_INTEGRAL(
                comparator.BOUND_BINDER, 0, 2,
                sp.Function("f")(comparator.BOUND_BINDER),
            ),
        )

    def test_native_boolean_is_not_a_residual_operand(self) -> None:
        value = comparator.typed_residual(True, sp.Integer(1), "SYNTHETIC", 0.1, 2)
        self.assertIsInstance(value, comparator.BooleanNotResidualable)

    def test_nested_boolean_does_not_block_algebraic_sibling(self) -> None:
        x = sp.Symbol("x")
        left = comparator.Association((("BOOL", True), ("ALGEBRAIC", x + 1)))
        right = comparator.Association((("BOOL", False), ("ALGEBRAIC", x)))
        value = comparator.typed_residual(left, right, "SYNTHETIC", 0.1, 20)
        self.assertIsInstance(value, comparator.ResidualAssociation)
        fields = dict(value.entries)
        self.assertIsInstance(fields["BOOL"], comparator.BooleanNotResidualable)
        self.assertEqual(fields["ALGEBRAIC"], 1)

    def test_undecided_and_unsupported_are_not_parse_failures(self) -> None:
        x = sp.Symbol("x")
        undecided = comparator.typed_residual(
            comparator.TextAtom("UNDECIDED"), x, "SYNTHETIC_STATUS", 0.1, 10
        )
        unsupported = comparator.typed_residual({}, {}, "SYNTHETIC", 0.1, 4)
        self.assertIsInstance(undecided, comparator.UndecidedResidual)
        self.assertIsInstance(unsupported, comparator.ResidualFailure)
        case = comparator.build_case(
            "PY", comparator.make_key(()), "raw", lambda: (_ for _ in ()).throw(ValueError("bad"))
        )
        self.assertTrue(comparator.materialize(case))
        self.assertIn("ValueError", case.error or "")

    def test_undecided_object_is_emitted_as_the_residual(self) -> None:
        key = comparator.make_key((("OBJECT", "SYNTHETIC"),))
        py = comparator.ParsedCase(
            key, comparator.TextAtom("UNDECIDED"), "UNDECIDED",
            compare_value=comparator.TextAtom("UNDECIDED"),
        )
        wl = comparator.ParsedCase(
            key, sp.Symbol("x"), "x", compare_value=sp.Symbol("x")
        )
        stream = io.StringIO()
        with contextlib.redirect_stdout(stream):
            accounting = comparator.compare_family(
                "SYNTHETIC_STATUS", comparator.FamilyCases([py], [wl]), budget=0.1
            )
        self.assertEqual(accounting.parse_failed, 0)
        self.assertIn("A_minus_B = UndecidedResidual", stream.getvalue())

    def test_parse_failure_is_emitted_and_accounted_separately(self) -> None:
        key = comparator.make_key((("OBJECT", "SYNTHETIC"),))
        py = comparator.build_case(
            "PY", key, "bad", lambda: (_ for _ in ()).throw(ValueError("bad"))
        )
        wl = comparator.wl_case(key, "x")
        stream = io.StringIO()
        with contextlib.redirect_stdout(stream):
            accounting = comparator.compare_family(
                "SYNTHETIC", comparator.FamilyCases([py], [wl]), budget=0.1
            )
        self.assertEqual(accounting.parse_failed, 1)
        self.assertIn("<PARSE_FAILED ValueError: bad", stream.getvalue())

    def test_oversized_scalar_uses_symbolic_difference(self) -> None:
        x, y = sp.symbols("x y")
        value = comparator.typed_residual(x, y, "SYNTHETIC", 0.1, 200_001)
        self.assertIsInstance(value, comparator.SymbolicDifference)


class RepointAndRawControlTests(unittest.TestCase):
    def test_repointing_a_different_object_moves_the_residual(self) -> None:
        key = comparator.make_key((("OBJECT", "TARGET"),))
        py_objects = {"TARGET": "Symbol('x')", "NEIGHBOUR": "Symbol('y')"}
        wl_objects = {"TARGET": "x"}
        baseline = comparator.typed_residual(
            comparator.parse_py_value(py_objects["TARGET"], key, coefficient=False),
            comparator.parse_wl_value(wl_objects["TARGET"], key),
            "SYNTHETIC", 0.1, 20,
        )
        repointed = comparator.typed_residual(
            comparator.parse_py_value(py_objects["NEIGHBOUR"], key, coefficient=False),
            comparator.parse_wl_value(wl_objects["TARGET"], key),
            "SYNTHETIC", 0.1, 20,
        )
        self.assertEqual(baseline, 0)
        self.assertNotEqual(repointed, 0)
        baseline_cases = comparator.FamilyCases(
            [comparator.py_case(key, py_objects["TARGET"])],
            [comparator.wl_case(key, wl_objects["TARGET"])],
        )
        repointed_cases = comparator.FamilyCases(
            [comparator.py_case(key, py_objects["NEIGHBOUR"])],
            [comparator.wl_case(key, wl_objects["TARGET"])],
        )
        baseline_stream, repointed_stream = io.StringIO(), io.StringIO()
        with contextlib.redirect_stdout(baseline_stream):
            comparator.compare_family("SYNTHETIC", baseline_cases, budget=0.1)
        with contextlib.redirect_stdout(repointed_stream):
            comparator.compare_family("SYNTHETIC", repointed_cases, budget=0.1)
        self.assertIn("A_minus_B = Integer(0)", baseline_stream.getvalue())
        self.assertNotIn("A_minus_B = Integer(0)", repointed_stream.getvalue())

    def test_raw_control_is_refused_when_scalar_leaf_is_reachable(self) -> None:
        self.assertEqual(
            comparator.RAW_CONTROL_WHITELIST,
            {
                ("DTN_OPERATOR", "WHOLE_OBJECT"),
                ("NONINVERTIBILITY_CONDITION", "OPERATOR"),
            },
        )
        for family, leaf in comparator.RAW_CONTROL_WHITELIST:
            case = comparator.raw_control_case(
                "PY", comparator.make_key(()), "Tuple(Symbol('x'))",
                family=family, leaf=leaf,
            )
            self.assertIsInstance(case.value, comparator.RawCAS)
        wl_operator = comparator.raw_control_case(
            "WL", comparator.make_key(()), "OperatorSum[x,y]",
            family="DTN_OPERATOR", leaf="WHOLE_OBJECT",
        )
        self.assertIsInstance(wl_operator.value, comparator.RawCAS)
        self.assertIsNone(wl_operator.error)
        with self.assertRaises(comparator.InputError):
            comparator.raw_control_case(
                "PY", comparator.make_key(()),
                "Tuple(Tuple(Str('OBJECT'), Symbol('x')))",
                family="CONTROL_FORM_BASE_OPERAND", leaf="OBJECT",
            )
        with self.assertRaises(comparator.InputError):
            comparator.raw_control_case(
                "WL", comparator.make_key(()), "x",
                family="NONINVERTIBILITY_CONDITION", leaf="DETERMINANT",
            )

    def test_control_form_descends_to_inner_face_scalars(self) -> None:
        py_payload = (
            "Tuple(Tuple(Tuple(Str('DTN'), Str('LAB_HELD'), "
            "Str('RHO4_CONSTANT'), Integer(1)), "
            "Tuple(Tuple(Integer(1), Symbol('left')), "
            "Tuple(Integer(-1), Symbol('same')))))"
        )
        wl_payload = (
            '<|"DTN|LAB_HELD|RHO4_CONSTANT|DIRECTION_1" -> '
            '<|"OBJECT" -> <|"FACE_1" -> right, "FACE_-1" -> same|>|>|>'
        )
        cases = comparator.extract_family(
            "CONTROL_FORM_BASE_OPERAND", py_payload, wl_payload
        )
        self.assertEqual(len(cases.py), 3)
        self.assertEqual(len(cases.wl), 2)
        for case in (*cases.py, *cases.wl):
            comparator.materialize(case)
            self.assertNotIsInstance(case.compare_value, comparator.RawCAS)
        stream = io.StringIO()
        with contextlib.redirect_stdout(stream):
            accounting = comparator.compare_family(
                "CONTROL_FORM_BASE_OPERAND", cases, budget=0.1
            )
        self.assertEqual(accounting.join, 2)
        self.assertEqual(accounting.py_only, 1)
        self.assertIn("Add(Symbol('left'), Mul(Integer(-1), Symbol('right')))", stream.getvalue())

    def test_control_form_scalar_wl_object_is_not_folded_to_upper_face(self) -> None:
        py_payload = (
            "Tuple(Tuple(Tuple(Str('DTN'), Str('LAB_HELD'), "
            "Str('RHO4_CONSTANT'), Integer(1)), "
            "Tuple(Tuple(Integer(1), Symbol('left')))))"
        )
        wl_payload = (
            '<|"DTN|LAB_HELD|RHO4_CONSTANT|DIRECTION_1" -> '
            '<|"OBJECT" -> right|>|>'
        )
        cases = comparator.extract_family(
            "CONTROL_FORM_BASE_OPERAND", py_payload, wl_payload
        )
        self.assertNotIn(("FACE", "1"), cases.wl[0].key)
        parent = [case for case in cases.py if ("FACE", "1") not in case.key]
        face = [case for case in cases.py if ("FACE", "1") in case.key]
        self.assertEqual(len(parent), 1)
        self.assertEqual(len(face), 1)
        self.assertEqual(parent[0].key, cases.wl[0].key)
        stream = io.StringIO()
        with contextlib.redirect_stdout(stream):
            accounting = comparator.compare_family(
                "CONTROL_FORM_BASE_OPERAND", cases, budget=0.1
            )
        self.assertEqual(accounting.join, 1)
        self.assertEqual(accounting.py_only, 1)
        self.assertIn("Mismatch", stream.getvalue())


class RunSemanticsTests(unittest.TestCase):
    def test_sound_speed_fold_requires_complete_real_py_vocabulary(self) -> None:
        required = comparator.REQUIRED_C1_PY_SYMBOLS
        for names in (required, required - {"c_s0"}, {"c_s0"}, set()):
            with self.subTest(names=sorted(names)), tempfile.TemporaryDirectory() as directory:
                root = Path(directory)
                py_path = root / "py.out"
                wl_path = root / "wl.out"
                vocabulary = sp.srepr(sp.Tuple(*(sp.Symbol(name) for name in sorted(names))))
                py_operand = "c_s0" if "c_s0" in names else "left"
                py_path.write_text(
                    f"PY_S11CC1_SYNTHETIC: Symbol('{py_operand}')\n"
                    f"PY_S11CC1_LOCAL_VOCABULARY: {vocabulary}\n",
                    encoding="utf-8",
                )
                wl_path.write_text("WL_S11CC1_SYNTHETIC: cS0\n", encoding="utf-8")
                stream = io.StringIO()
                with contextlib.redirect_stdout(stream):
                    code = comparator.run(["--py", str(py_path), "--wl", str(wl_path)])
                self.assertEqual(code, 0)
                rendered = stream.getvalue()
                self.assertIn("SPELLING_INJECTIVITY collisions=0", rendered)
                if names == required:
                    self.assertIn("reserved_names=12 mechanical_spellings=12", rendered)
                    self.assertIn("('cS0', 'c_s0')", rendered)
                    self.assertIn("operand_B = Symbol('c_s0')", rendered)
                    self.assertIn("A_minus_B = Integer(0)", rendered)
                else:
                    self.assertIn("active_c1_folds=[]", rendered)
                    self.assertIn("operand_B = Symbol('cS0')", rendered)
                    self.assertNotIn("A_minus_B = Integer(0)", rendered)

    def test_partial_synthetic_vocabulary_does_not_activate_spelling_folds(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            py_path = root / "py.out"
            wl_path = root / "wl.out"
            py_path.write_text(
                "PY_S11CC1_SYNTHETIC: Symbol('rho_m')\n", encoding="utf-8"
            )
            wl_path.write_text(
                "WL_S11CC1_SYNTHETIC: rhoM\n", encoding="utf-8"
            )
            stream = io.StringIO()
            with contextlib.redirect_stdout(stream):
                code = comparator.run(["--py", str(py_path), "--wl", str(wl_path)])
        rendered = stream.getvalue()
        self.assertEqual(code, 0)
        self.assertIn("active_c1_folds=[]", rendered)
        self.assertIn("operand_B = Symbol('rhoM')", rendered)
        self.assertNotIn("A_minus_B = Integer(0)", rendered)

    def test_disagreement_is_measurement_and_exits_zero(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            py_path = root / "py.out"
            wl_path = root / "wl.out"
            py_path.write_text(
                "PY_S11CC1_SYNTHETIC: Symbol('left')\n"
                "PY_S11CC1_LOCAL_NOTE: Str('excluded')\n",
                encoding="utf-8",
            )
            wl_path.write_text(
                "WL_S11CC1_SYNTHETIC: right\n"
                'WL_S11CC1_LOCAL_TAG_NAMES: <|"LOCAL_TAG_NAMES" -> {}|>\n',
                encoding="utf-8",
            )
            stream = io.StringIO()
            with contextlib.redirect_stdout(stream):
                code = comparator.run(["--py", str(py_path), "--wl", str(wl_path)])
        rendered = stream.getvalue()
        self.assertEqual(code, 0)
        self.assertLess(rendered.index("operand_A ="), rendered.index("ACCOUNTING SYNTHETIC"))
        self.assertIn("operand_B =", rendered)
        self.assertIn("A_minus_B =", rendered)
        self.assertIn("extracted_leaves_py=1", rendered)
        self.assertIn("extracted_leaves_wl=1", rendered)
        self.assertIn("RUN_ACCOUNTING families=1 families_with_join=1", rendered)
        self.assertIn("families_with_unpaired=0 parse_failed=0 duplicate_key=0", rendered)
        self.assertIn("MEASUREMENT_SCOPE", rendered)
        self.assertIn(
            "sections_1_to_2,supplied_substrate,mu_theta_operand,supplied_bookkeeping",
            rendered,
        )
        self.assertIn("residual_target=none", rendered)
        self.assertNotIn("ACCOUNTING LOCAL_", rendered)
        self.assertIn("LOCAL_INVENTORY engine=PY", rendered)
        self.assertIn("LOCAL_INVENTORY engine=WL", rendered)
        self.assertIsNone(
            re.search(
                r"^(?:PASS|FAIL|VERDICT|FINAL_STATUS|AGREE|DISAGREE|STATUS)(?:\s|:|=)",
                rendered,
                re.MULTILINE,
            )
        )

    def test_unpaired_family_is_visible_and_still_exits_zero(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            py_path = root / "py.out"
            wl_path = root / "wl.out"
            py_path.write_text(
                "PY_S11CC1_SYNTHETIC: Symbol('left')\n", encoding="utf-8"
            )
            wl_path.write_text(
                'WL_S11CC1_LOCAL_TAG_NAMES: <|"LOCAL_TAG_NAMES" -> {}|>\n',
                encoding="utf-8",
            )
            stream = io.StringIO()
            with contextlib.redirect_stdout(stream):
                code = comparator.run(["--py", str(py_path), "--wl", str(wl_path)])
        self.assertEqual(code, 0)
        rendered = stream.getvalue()
        self.assertIn("ACCOUNTING SYNTHETIC {join=0, py_only=1, wl_only=0", rendered)
        self.assertIn("families_with_join=0 families_with_unpaired=1", rendered)

    def test_malformed_stream_is_operational_error(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            py_path = root / "py.out"
            wl_path = root / "wl.out"
            py_path.write_text("not a tagged row\n", encoding="utf-8")
            wl_path.write_text("WL_S11CC1_SYNTHETIC: x\n", encoding="utf-8")
            with contextlib.redirect_stderr(io.StringIO()):
                code = comparator.run(["--py", str(py_path), "--wl", str(wl_path)])
        self.assertNotEqual(code, 0)


if __name__ == "__main__":
    unittest.main()
