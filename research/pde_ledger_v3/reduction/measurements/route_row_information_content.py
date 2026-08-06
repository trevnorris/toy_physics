#!/usr/bin/env python3
"""What does the `q2_downstream_route` row report when the routes actually differ?

Thirteen cross-engine rows compare a bare route token -- Mathematica's
`quadraticFormRoute` against SymPy's `M_B` -- and every one reports
NAMING_MISMATCH, because any single symbol can be relabelled onto any other.

This script substitutes a series of replacement payloads for the SymPy side of
one such row and prints the verdict the comparator returns for each, alongside
the verdict for the committed payloads.  It prints the verdicts; it does not say
whether the row is worth keeping.
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parents[2]
HARNESS = ROOT / "reduction" / "engine_output_checks.py"
sys.path.insert(0, str(HARNESS.parent))

spec = importlib.util.spec_from_file_location("harness", HARNESS)
harness = importlib.util.module_from_spec(spec)
sys.modules["harness"] = harness
spec.loader.exec_module(harness)

WL_PAYLOAD = "quadraticFormRoute"

SUBSTITUTIONS = (
    ("the committed SymPy payload", "M_B"),
    ("the same token as Mathematica", "quadraticFormRoute"),
    ("a token naming the OTHER route", "positionSpaceEulerLagrangeRoute"),
    ("a token naming no route at all", "banana"),
    ("a token asserting the routes disagreed", "routesDisagreed"),
    ("a token asserting the step failed", "FAILED"),
)


def parsed(value: object) -> object:
    return harness.ParsedValue(harness.value_kind(value), value)


def main() -> int:
    row = {
        "quantity": "q2_downstream_route",
        "wl": "W",
        "py": "P",
        "cardinality": {"kind": "scalar"},
    }
    config = {
        "symbol_naming": {
            "rule": "registry_snake_case_to_lower_camel",
            "engine_styles": {"wl": "lower_camel", "py": "canonical"},
            "exceptions": [],
        }
    }
    print(f"WL payload held fixed at: {WL_PAYLOAD}")
    print()
    for description, py_payload in SUBSTITUTIONS:
        outputs = {
            "wl": {"W": parsed(sp.Symbol(WL_PAYLOAD))},
            "py": {"P": parsed(sp.Symbol(py_payload))},
        }
        report = harness.check_cross_engine([row], outputs, config=config)
        observed = report.rows[0]
        print(f"py = {py_payload!r}   ({description})")
        print(f"  status               = {observed.status}")
        print(f"  undeclared_spellings = {observed.undeclared_spellings}")
        print(f"  counts as a verdict  = {observed.status in harness.VERDICTS}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
