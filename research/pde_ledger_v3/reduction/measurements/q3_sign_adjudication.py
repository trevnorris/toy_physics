#!/usr/bin/env python3
"""Ask the CAS what each disputed root's sign is under the declared premises.

Three cross-engine rows report DISAGREE on `q3_sign`.  In each the two engines
emit different things, but a difference between `1` and `Sign(...)` -- or
between `Sign(...)` and an undecidability token -- is not the same kind of event
as one engine saying `+1` where the other says `-1`.

This script takes each disputed root as the engines emitted it, states the
premise set both engines declare, and asks SymPy to decide the root's sign under
exactly those premises.  It prints the root, the premises used, each engine's
emitted payload, and the CAS's own determination.  It does not say which engine
is right; it prints what the premises entail.
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parents[2]
WL_OUT = ROOT / "mathematica" / "out" / "S10_brane_mode_spectrum_mathematica_audit.out"
PY_OUT = ROOT / "scripts" / "out" / "S10_brane_mode_spectrum_sympy_audit.out"

mu_R, rho_br, s_rho = sp.symbols("mu_R rho_br s_rho", positive=True)
k = sp.symbols("k1 k2 k3 k4 k5", real=True)

# (row, WL tag, PY tag, root expression as the engines emitted it, extra premises)
CASES = (
    (
        "main_d5_root2_q3_sign",
        "WL_S10_MAIN_D5_ROOT2_Q3_SIGN",
        "PY_S10_MAIN_D5_ROOT2_Q3_SIGN",
        mu_R * (k[0] ** 2 + k[1] ** 2 + k[2] ** 2 + k[3] ** 2 + k[4] ** 2) / rho_br,
        sum(component ** 2 for component in k[:5]) > 0,
    ),
    (
        "xform_aniso_d3_root3_q3_sign",
        "WL_S10_XFORM_ANISO_D3_ROOT3_Q3_SIGN",
        "PY_S10_XFORM_ANISO_D3_ROOT3_Q3_SIGN",
        mu_R * (k[0] ** 2 * s_rho + k[1] ** 2 + k[2] ** 2) / (rho_br * s_rho),
        sum(component ** 2 for component in k[:3]) > 0,
    ),
    (
        "xform_aniso_d4_root3_q3_sign",
        "WL_S10_XFORM_ANISO_D4_ROOT3_Q3_SIGN",
        "PY_S10_XFORM_ANISO_D4_ROOT3_Q3_SIGN",
        mu_R * (k[0] ** 2 * s_rho + k[1] ** 2 + k[2] ** 2 + k[3] ** 2) / (rho_br * s_rho),
        sum(component ** 2 for component in k[:4]) > 0,
    ),
)

PREMISES = "mu_R > 0, rho_br > 0, s_rho > 0, every k component real, sum of k squared > 0"


def read_tag(path: Path, tag: str) -> str:
    with path.open() as handle:
        for line in handle:
            head, separator, tail = line.partition(": ")
            if separator and head == tag:
                return tail.rstrip("\n")
    return "<no tag>"


def decide(expression: sp.Expr, norm_positive: sp.Expr) -> tuple[object, object]:
    """Is the expression positive on the region the premises carve out?

    The premise `sum k^2 > 0` is not expressible as a symbol assumption, so it
    is discharged by cases: the expression's numerator is a positive-weighted
    sum of squares, and the region excludes the single point where all squares
    vanish.  Both the direct query and the by-cases minimum are printed.
    """
    direct = sp.ask(sp.Q.positive(expression))
    numerator, denominator = sp.fraction(sp.together(expression))
    # Minimum of the numerator over the region, attained when all but one
    # component vanish; report the per-component coefficient signs instead of
    # asserting the minimum.
    components = sorted(numerator.free_symbols & set(k), key=str)
    per_component = {}
    for component in components:
        substitution = {other: 0 for other in components if other is not component}
        per_component[str(component)] = sp.simplify(
            numerator.subs(substitution) / component ** 2
        )
    return direct, (sp.simplify(denominator), per_component)


def main() -> int:
    print(f"WL_SOURCE: {WL_OUT}")
    print(f"PY_SOURCE: {PY_OUT}")
    print(f"PREMISES: {PREMISES}")
    print()
    for row, wl_tag, py_tag, expression, norm_positive in CASES:
        direct, (denominator, per_component) = decide(expression, norm_positive)
        print(f"ROW={row}")
        print(f"  root                    = {expression}")
        print(f"  wl_emitted              = {read_tag(WL_OUT, wl_tag)}")
        print(f"  py_emitted              = {read_tag(PY_OUT, py_tag)}")
        print(f"  ask(Q.positive(root))   = {direct}")
        print(f"  denominator             = {denominator}")
        print("  numerator coefficient of each k_i^2 with the others zeroed:")
        for name, coefficient in per_component.items():
            print(f"    {name}^2 : {coefficient}   positive={sp.ask(sp.Q.positive(coefficient))}")
        print()
    return 0


if __name__ == "__main__":
    sys.exit(main())
