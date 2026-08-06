#!/usr/bin/env python3
"""Recover each package's stiffness coefficient from its own emitted Lagrangian.

Each engine emits `Q6_STIFFNESS_COEFFICIENTS` as a reported quantity and its
Lagrangian as a separate object.  This script does not read the reported
quantity.  It takes the emitted Lagrangian at D=3, discards every term carrying
a time derivative, substitutes the independent gradient symbols for the spatial
derivatives, and divides the remainder by that package's emitted §Q7 stiffness
density.

For each package and engine it prints: the stiffness part of the Lagrangian, the
§Q7 density, the quotient, the residual of quotient x density against the
stiffness part, and -- separately -- the coefficient the engine itself reported.
It draws no conclusion from any of them.
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parents[2]
WL_OUT = ROOT / "mathematica" / "out" / "S10_brane_mode_spectrum_mathematica_audit.out"
PY_OUT = ROOT / "scripts" / "out" / "S10_brane_mode_spectrum_sympy_audit.out"

PACKAGES = (
    "MAIN",
    "XFORM_FULLGRAD",
    "XFORM_DIVONLY",
    "XFORM_SIGNFLIP",
    "XFORM_ANISO",
    "XCOEF_SCALE",
)
N = 3

# SymPy emits fields as u{j}(t, x1..xD); Mathematica as
# Derivative[a1..aD, at][u{j}][x1..xD, t].
PY_DERIVATIVE = re.compile(r"Derivative\(u(\d)\(t(?:, x\d)+\), (t|x\d)\)")
WL_DERIVATIVE = re.compile(r"Derivative\[([0-9, ]+)\]\[u(\d)\]\[[^\]]*\]")


def read_tags(path: Path) -> dict[str, str]:
    table: dict[str, str] = {}
    with path.open() as handle:
        for line in handle:
            head, separator, tail = line.partition(": ")
            if separator and re.fullmatch(r"[A-Z0-9_]+", head):
                table[head] = tail.rstrip("\n")
    return table


def py_lagrangian_to_symbols(text: str) -> sp.Expr:
    """t-derivatives become the marker T{j}; x-derivatives become g{i}{j}."""

    def rename(match: re.Match[str]) -> str:
        field, variable = match.group(1), match.group(2)
        if variable == "t":
            return f"T{field}"
        return f"g{variable[1:]}{field}"

    return sp.sympify(PY_DERIVATIVE.sub(rename, text))


def wl_lagrangian_to_symbols(text: str) -> sp.Expr:
    """Mathematica's multi-index is (x1..xD, t) -- see the emitted argument list."""

    def rename(match: re.Match[str]) -> str:
        orders = [int(piece) for piece in match.group(1).split(",")]
        field = match.group(2)
        spatial, time_order = orders[:N], orders[N]
        if time_order == 1 and sum(spatial) == 0:
            return f"T{field}"
        if time_order == 0 and sum(spatial) == 1:
            return f"g{spatial.index(1) + 1}{field}"
        raise ValueError(f"unhandled derivative multi-index {orders}")

    converted = WL_DERIVATIVE.sub(rename, text).replace("^", "**")
    # Mathematica prints the wl spellings of the two coefficients.
    return sp.sympify(converted)


def drop_time_terms(expression: sp.Expr) -> sp.Expr:
    markers = {sp.Symbol(f"T{index}") for index in range(1, N + 1)}
    kept = [
        term
        for term in sp.Add.make_args(sp.expand(expression))
        if not (term.free_symbols & markers)
    ]
    return sp.expand(sp.Add(*kept)) if kept else sp.Integer(0)


def transpose_free_density(text: str) -> sp.Expr:
    return sp.expand(sp.sympify(re.sub(r"\bg(\d)x(\d)\b", r"g\1\2", text.replace("^", "**"))))


def main() -> int:
    wl = read_tags(WL_OUT)
    py = read_tags(PY_OUT)
    print(f"WL_SOURCE: {WL_OUT}")
    print(f"PY_SOURCE: {PY_OUT}")
    print(f"DIMENSION: {N}")
    print()

    for package in PACKAGES:
        for engine, tags, lagrangian_suffix, density_suffix, converter, density_converter in (
            (
                "wl",
                wl,
                "Q1_LAGRANGIAN",
                "Q7_PACKAGE_STIFFNESS_DENSITY",
                wl_lagrangian_to_symbols,
                transpose_free_density,
            ),
            (
                "py",
                py,
                "Q1_LAGRANGIAN_EXPANDED",
                "Q7_STIFFNESS",
                py_lagrangian_to_symbols,
                lambda text: sp.expand(sp.sympify(text)),
            ),
        ):
            prefix = f"{engine.upper()}_S10_{package}_D{N}_"
            lagrangian_tag = prefix + lagrangian_suffix
            density_tag = prefix + density_suffix
            reported_tag = prefix + "Q6_STIFFNESS_COEFFICIENTS"
            if lagrangian_tag not in tags or density_tag not in tags:
                print(f"PACKAGE={package} ENGINE={engine}  MISSING {lagrangian_tag} or {density_tag}")
                continue
            stiffness_part = drop_time_terms(converter(tags[lagrangian_tag]))
            density = density_converter(tags[density_tag])
            quotient = sp.simplify(sp.cancel(stiffness_part / density)) if density != 0 else None
            residual = (
                sp.expand(quotient * density - stiffness_part) if quotient is not None else None
            )
            print(f"PACKAGE={package} ENGINE={engine}")
            print(f"  stiffness_part_of_L = {stiffness_part}")
            print(f"  q7_density          = {density}")
            print(f"  quotient            = {quotient}")
            print(f"  residual            = {residual}")
            print(f"  engine_reported_Q6_STIFFNESS_COEFFICIENTS = {tags.get(reported_tag, '<no tag>')}")
        print()
    return 0


if __name__ == "__main__":
    sys.exit(main())
