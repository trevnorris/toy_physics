#!/usr/bin/env python3
"""Ask the Q7 payloads which index map relates the two engines' gradient symbols.

The two engines spell the D=3 gradient symbols differently -- Mathematica emits
g{row}x{column}, SymPy emits g{row}{column}.  Two candidate maps carry one
spelling onto the other:

    IDENTITY   g{r}x{c} -> g{r}{c}
    TRANSPOSE  g{r}x{c} -> g{c}{r}

This script substitutes each candidate into every emitted Mathematica Q7 payload
and prints, per package and per object, both operands and the residual against
the SymPy payload.  It states nothing about which map is correct; it prints what
each map produces.
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

# (object name, Mathematica suffix, SymPy suffix)
OBJECTS = (
    ("package_stiffness_density", "Q7_PACKAGE_STIFFNESS_DENSITY", "Q7_STIFFNESS"),
    ("ordinary_curl_norm", "Q7_ORDINARY_CURL_NORM", "Q7_CURL_DOT"),
    ("residual", "Q7_PACKAGE_STIFFNESS_VS_ORDINARY_CURL_RESIDUAL", "Q7_DIFFERENCE"),
)

WL_SYMBOL = re.compile(r"\bg([1-9])x([1-9])\b")


def read_tags(path: Path) -> dict[str, str]:
    table: dict[str, str] = {}
    with path.open() as handle:
        for line in handle:
            head, separator, tail = line.partition(": ")
            if separator and re.fullmatch(r"[A-Z0-9_]+", head):
                table[head] = tail.rstrip("\n")
    return table


def wl_to_sympy_text(text: str, transpose: bool) -> str:
    def rename(match: re.Match[str]) -> str:
        row, column = match.group(1), match.group(2)
        return f"g{column}{row}" if transpose else f"g{row}{column}"

    return WL_SYMBOL.sub(rename, text.replace("^", "**"))


def main() -> int:
    wl = read_tags(WL_OUT)
    py = read_tags(PY_OUT)

    print(f"WL_SOURCE: {WL_OUT}")
    print(f"PY_SOURCE: {PY_OUT}")
    print()

    missing: list[str] = []
    for package in PACKAGES:
        for name, wl_suffix, py_suffix in OBJECTS:
            wl_tag = f"WL_S10_{package}_D3_{wl_suffix}"
            py_tag = f"PY_S10_{package}_D3_{py_suffix}"
            if wl_tag not in wl or py_tag not in py:
                missing.append(f"{wl_tag if wl_tag not in wl else py_tag}")
                continue
            py_expression = sp.expand(sp.sympify(py[py_tag]))
            for map_name, transpose in (("IDENTITY", False), ("TRANSPOSE", True)):
                mapped = sp.expand(sp.sympify(wl_to_sympy_text(wl[wl_tag], transpose)))
                residual = sp.expand(mapped - py_expression)
                print(f"PACKAGE={package} OBJECT={name} MAP={map_name}")
                print(f"  wl_mapped = {mapped}")
                print(f"  py        = {py_expression}")
                print(f"  residual  = {residual}")
            print()

    print(f"MISSING_TAGS ({len(missing)}): {', '.join(missing) if missing else 'none'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
