#!/usr/bin/env python3
"""How many relabellings of the nine gradient symbols do the Q7 payloads admit?

The declaration in `reduction/checks_S10.yaml` fixes one map from WL's g{r}x{c}
to PY's g{r}{c}.  The payloads cannot check it: an earlier measurement showed the
transposed map reproduces every payload exactly.  This script asks the complete
question -- over all 9! bijections of the nine symbols, how many leave EVERY
emitted D=3 stiffness density invariant?

Each density is read from the committed outputs, expanded, and reduced to a
mapping from degree-2 monomial (an unordered symbol pair) to coefficient.  A
permutation is admitted when it carries every density's coefficient mapping onto
itself.  The script prints the count, the subgroup structure it observes, and
several admitted permutations by name.  It states no conclusion about which
permutation is correct -- that is settled from the engine definitions, not here.
"""

from __future__ import annotations

import itertools
import re
import sys
from collections import defaultdict
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parents[2]
PY_OUT = ROOT / "scripts" / "out" / "S10_brane_mode_spectrum_sympy_audit.out"

# The distinct stiffness densities across the six packages; the emitted
# residuals are differences of these and add no independent constraint.
DENSITY_TAGS = (
    "PY_S10_MAIN_D3_Q7_STIFFNESS",          # curl
    "PY_S10_XFORM_FULLGRAD_D3_Q7_STIFFNESS",
    "PY_S10_XFORM_DIVONLY_D3_Q7_STIFFNESS",
    "PY_S10_MAIN_D3_Q7_CURL_DOT",           # the computed reference
)

SYMBOLS = [f"g{r}{c}" for r in (1, 2, 3) for c in (1, 2, 3)]
INDEX = {name: position for position, name in enumerate(SYMBOLS)}


def read_tags(path: Path) -> dict[str, str]:
    table: dict[str, str] = {}
    with path.open() as handle:
        for line in handle:
            head, separator, tail = line.partition(": ")
            if separator and re.fullmatch(r"[A-Z0-9_]+", head):
                table[head] = tail.rstrip("\n")
    return table


def monomial_map(expression: sp.Expr) -> dict[tuple[int, int], int]:
    """Degree-2 form -> {(i, j) with i <= j: integer coefficient}."""
    result: dict[tuple[int, int], int] = {}
    poly = sp.Poly(sp.expand(expression), *[sp.Symbol(name) for name in SYMBOLS])
    for exponents, coefficient in poly.terms():
        positions: list[int] = []
        for position, power in enumerate(exponents):
            positions.extend([position] * power)
        if len(positions) != 2:
            raise ValueError(f"non-quadratic monomial {exponents}")
        key = (min(positions), max(positions))
        result[key] = int(coefficient)
    return result


def permute(form: dict[tuple[int, int], int], sigma: tuple[int, ...]) -> dict[tuple[int, int], int]:
    out: dict[tuple[int, int], int] = {}
    for (i, j), coefficient in form.items():
        a, b = sigma[i], sigma[j]
        out[(min(a, b), max(a, b))] = coefficient
    return out


def describe(sigma: tuple[int, ...]) -> str:
    return ", ".join(
        f"{SYMBOLS[position]}->{SYMBOLS[image]}"
        for position, image in enumerate(sigma)
        if position != image
    ) or "identity"


def main() -> int:
    tags = read_tags(PY_OUT)
    forms = []
    print(f"SOURCE: {PY_OUT}")
    for tag in DENSITY_TAGS:
        form = monomial_map(sp.sympify(tags[tag]))
        forms.append(form)
        print(f"  {tag}: monomials={len(form)}")
    print()

    admitted: list[tuple[int, ...]] = []
    for sigma in itertools.permutations(range(9)):
        if all(permute(form, sigma) == form for form in forms):
            admitted.append(sigma)

    print(f"TOTAL_PERMUTATIONS_TESTED: {sp.factorial(9)}")
    print(f"ADMITTED_PERMUTATIONS: {len(admitted)}")
    print()

    diagonal = {INDEX[f"g{i}{i}"] for i in (1, 2, 3)}
    preserves_diagonal = sum(
        1 for sigma in admitted if {sigma[position] for position in diagonal} == diagonal
    )
    print(f"ADMITTED_THAT_MAP_DIAGONAL_TO_DIAGONAL: {preserves_diagonal}")

    transpose = tuple(INDEX[f"g{c}{r}"] for r in (1, 2, 3) for c in (1, 2, 3))
    print(f"TRANSPOSE_ADMITTED: {transpose in set(admitted)}")
    diagonal_cycle = list(range(9))
    for i in (1, 2, 3):
        diagonal_cycle[INDEX[f"g{i}{i}"]] = INDEX[f"g{i % 3 + 1}{i % 3 + 1}"]
    print(f"DIAGONAL_3CYCLE_ADMITTED: {tuple(diagonal_cycle) in set(admitted)}")
    print()

    by_size: defaultdict[int, list[str]] = defaultdict(list)
    for sigma in admitted:
        by_size[sum(1 for position, image in enumerate(sigma) if position != image)].append(
            describe(sigma)
        )
    print("ADMITTED PERMUTATIONS BY NUMBER OF MOVED SYMBOLS")
    for size in sorted(by_size):
        print(f"  moved={size} count={len(by_size[size])}")
        for text in by_size[size][:3]:
            print(f"    {text}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
