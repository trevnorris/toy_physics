#!/usr/bin/env python3
"""Regression tests for reusable reduction-registry machinery."""

from __future__ import annotations

import subprocess
import sys
import unittest
from pathlib import Path

import sympy as sp

from registry_read import (
    DeclaredValueDefaultWarning,
    EvaluationError,
    Registry,
    RegistryValidationError,
    load_raw_documents,
    load_registry,
)


HERE = Path(__file__).resolve().parent


class ReducibleDimensionRegressionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.registry = load_registry()

    def test_plane_union_line_uses_maximum_component_dimension(self) -> None:
        x, y, z = sp.symbols("x y z", real=True)
        observed = self.registry.constraint_dimension((z * x, z * y), (x, y, z))
        self.assertEqual(observed, 2)

    def test_positive_dim3_union_dim1_uses_maximum_component_dimension(self) -> None:
        x, y, z, w = sp.symbols("x y z w", real=True)
        constraints = (
            (x - y) * (x - 2 * y),
            (x - y) * (z - 3 * y),
            (x - y) * (w - 4 * y),
        )
        observed = self.registry.constraint_dimension(constraints, (x, y, z, w))
        self.assertEqual(observed, 3)
        self.registry.certify_positive_real_dimension(
            constraints,
            (x, y, z, w),
            dimension=observed,
        )

    def test_zero_constraint_certifies_full_ambient_dimension(self) -> None:
        x = sp.Symbol("x", real=True)
        self.registry.certify_positive_real_dimension((0,), (x,), dimension=1)


class AbleToFailProtocolRegressionTests(unittest.TestCase):
    def test_empty_harness_is_error(self) -> None:
        completed = subprocess.run(
            [sys.executable, str(HERE / "able_to_fail.py"), "--demonstrate-empty"],
            check=False,
            capture_output=True,
            text=True,
            timeout=600,
        )
        self.assertEqual(completed.returncode, 1)
        self.assertIn("HARNESS_CASE_COUNT: ERROR expected=5 observed=0", completed.stdout)
        self.assertIn("ABLE_TO_FAIL_HARNESS: ERROR", completed.stdout)
        self.assertNotIn("ABLE_TO_FAIL_HARNESS: PASS", completed.stdout)

    def test_child_cannot_spoof_parent_verdict(self) -> None:
        completed = subprocess.run(
            [sys.executable, str(HERE / "able_to_fail.py"), "--demonstrate-spoof"],
            check=False,
            capture_output=True,
            text=True,
            timeout=600,
        )
        self.assertEqual(completed.returncode, 1)
        self.assertIn(
            "CHILD verdict-spoof stdout: ABLE_TO_FAIL_CHILD_TEXT: PASS",
            completed.stdout,
        )
        self.assertIn("ABLE_TO_FAIL_HARNESS: ERROR", completed.stdout)
        self.assertNotIn("ABLE_TO_FAIL_HARNESS: PASS", completed.stdout)

    def test_no_witness_exception_is_reported_as_error(self) -> None:
        completed = subprocess.run(
            [sys.executable, str(HERE / "able_to_fail.py"), "--demonstrate-crash"],
            check=False,
            capture_output=True,
            text=True,
            timeout=600,
        )
        self.assertEqual(completed.returncode, 1)
        self.assertIn(
            "DEMONSTRATION no-witness-crash: status=ERROR observed_exit=1",
            completed.stdout,
        )
        self.assertIn("ABLE_TO_FAIL_HARNESS: ERROR", completed.stdout)
        self.assertNotIn("ABLE_TO_FAIL_HARNESS: PASS", completed.stdout)


class LiteralConsistencyRegressionTests(unittest.TestCase):
    @staticmethod
    def _documents() -> tuple[dict, dict, dict]:
        return load_raw_documents()

    @staticmethod
    def _r1(relations: dict) -> dict:
        return next(
            row
            for row in relations["relations"]
            if row["relation_id"] == "R1"
        )

    def test_n_eos_has_declared_value(self) -> None:
        registry = load_registry()
        self.assertEqual(registry.quantities["Q.medium.n_eos"].value, 5)

    def test_declared_value_default_is_opt_in_and_reported(self) -> None:
        quantities, relations, schema = self._documents()
        big_k = next(
            row
            for row in quantities["quantities"]
            if row["qid"] == "Q.medium.K"
        )
        big_k["value"] = 1
        registry = Registry.from_documents(quantities, relations, schema)
        inputs = {"rho0": 1, "mass": 5}
        with self.assertRaisesRegex(
            EvaluationError,
            r"missing primitive/residue numeric input: Q\.medium\.K",
        ):
            registry.evaluate_output("h0", inputs)
        with self.assertWarnsRegex(
            DeclaredValueDefaultWarning,
            r"declared value defaults used: Q\.medium\.K=1",
        ):
            observed = registry.evaluate_output(
                "h0",
                inputs,
                allow_declared_defaults=True,
            )
        self.assertEqual(observed, sp.Rational(5, 4))

    def test_r1_coefficient_must_equal_n_eos(self) -> None:
        quantities, relations, schema = self._documents()
        self._r1(relations)["residual"][2][1][1][1] = 3
        with self.assertRaisesRegex(
            RegistryValidationError,
            r"R1: literal .* Q\.medium\.n_eos",
        ):
            Registry.from_documents(quantities, relations, schema)

    def test_r1_exponent_must_equal_n_eos_minus_one(self) -> None:
        quantities, relations, schema = self._documents()
        self._r1(relations)["residual"][2][1][1][3][2] = 3
        with self.assertRaisesRegex(
            RegistryValidationError,
            r"R1: literal .* Q\.medium\.n_eos",
        ):
            Registry.from_documents(quantities, relations, schema)

    def test_changing_n_eos_rechecks_r1_literals(self) -> None:
        quantities, relations, schema = self._documents()
        n_eos = next(
            row
            for row in quantities["quantities"]
            if row["qid"] == "Q.medium.n_eos"
        )
        n_eos["value"] = 3
        with self.assertRaisesRegex(
            RegistryValidationError,
            r"R1: literal .* Q\.medium\.n_eos",
        ):
            Registry.from_documents(quantities, relations, schema)


if __name__ == "__main__":
    unittest.main()
