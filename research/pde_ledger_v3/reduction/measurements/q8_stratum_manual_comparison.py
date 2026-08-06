#!/usr/bin/env python3
"""Compare the Q8 stratum re-run across engines by hand.

§8's tag grammar has no scope token for a stratum, so the two engines emit
different raw names for these tags -- Mathematica

    WL_S10_<package>_D<n>_STRATUM<s>_ROOT<r>_<suffix>

and SymPy

    PY_S10_<package>_D<n>_Q8_STRATUM<s>_ROOT<r>_<suffix>

-- and `reduction/checks_S10.yaml` declares no cross-engine pair for any of
them.  The harness therefore never compares the stratum mode counts.  This
script pairs them by (package, dimension, stratum, root, suffix) and prints both
operands for every pair it can form, plus every tag on either side that has no
partner.  It emits no verdict.
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
WL_OUT = ROOT / "mathematica" / "out" / "S10_brane_mode_spectrum_mathematica_audit.out"
PY_OUT = ROOT / "scripts" / "out" / "S10_brane_mode_spectrum_sympy_audit.out"

# Scalar count-bearing suffixes only -- the matrices and bases are compared
# nowhere and their representations differ by engine.
SUFFIXES = (
    "N2_RANK",
    "N2_NULLITY",
    "N3_STACKED_RANK",
    "N3_TRANSVERSE_NULLITY",
    "N4_NULLITY_DIFFERENCE",
    "N7_BASIS_COUNT",
    "N7_COUNT_RESIDUAL",
    "N7_BASIS_COUNT_RESIDUAL",
)

WL_KEY = re.compile(
    r"^WL_S10_(?P<package>[A-Z_]+?)_D(?P<dimension>\d+)_STRATUM(?P<stratum>\d+)_"
    r"ROOT(?P<root>\d+)_(?P<suffix>[A-Z0-9_]+)$"
)
PY_KEY = re.compile(
    r"^PY_S10_(?P<package>[A-Z_]+?)_D(?P<dimension>\d+)_Q8_STRATUM(?P<stratum>\d+)_"
    r"ROOT(?P<root>\d+)_(?P<suffix>[A-Z0-9_]+)$"
)


def read_tags(path: Path) -> dict[str, str]:
    table: dict[str, str] = {}
    with path.open() as handle:
        for line in handle:
            head, separator, tail = line.partition(": ")
            if separator and re.fullmatch(r"[A-Z0-9_]+", head):
                table[head] = tail.rstrip("\n")
    return table


def index(tags: dict[str, str], pattern: re.Pattern[str]) -> dict[tuple[str, ...], tuple[str, str]]:
    result: dict[tuple[str, ...], tuple[str, str]] = {}
    for tag, value in tags.items():
        match = pattern.fullmatch(tag)
        if match is None:
            continue
        key = (
            match["package"],
            match["dimension"],
            match["stratum"],
            match["root"],
            match["suffix"],
        )
        result[key] = (tag, value)
    return result


def main() -> int:
    wl = index(read_tags(WL_OUT), WL_KEY)
    py = index(read_tags(PY_OUT), PY_KEY)

    print(f"WL_SOURCE: {WL_OUT}")
    print(f"PY_SOURCE: {PY_OUT}")
    print(f"WL_STRATUM_ROOT_TAGS: {len(wl)}")
    print(f"PY_STRATUM_ROOT_TAGS: {len(py)}")
    print()

    counted = [key for key in sorted(set(wl) | set(py)) if key[4] in SUFFIXES]
    print(f"COUNT_BEARING_KEYS: {len(counted)}")
    for key in counted:
        package, dimension, stratum, root, suffix = key
        wl_entry = wl.get(key)
        py_entry = py.get(key)
        print(f"{package} D{dimension} STRATUM{stratum} ROOT{root} {suffix}")
        print(f"  wl = {wl_entry[1] if wl_entry else '<no tag>'}   [{wl_entry[0] if wl_entry else '-'}]")
        print(f"  py = {py_entry[1] if py_entry else '<no tag>'}   [{py_entry[0] if py_entry else '-'}]")
    print()

    wl_only = sorted(key for key in wl if key not in py)
    py_only = sorted(key for key in py if key not in wl)
    print(f"WL_ONLY_KEYS ({len(wl_only)}): {', '.join('/'.join(k) for k in wl_only) or 'none'}")
    print(f"PY_ONLY_KEYS ({len(py_only)}): {', '.join('/'.join(k) for k in py_only) or 'none'}")

    print()
    print("GENERIC_COUNTERPARTS (same objects away from the stratum, for contrast)")
    wl_tags = read_tags(WL_OUT)
    py_tags = read_tags(PY_OUT)
    for package, dimension in (("XFORM_ANISO", "3"), ("XFORM_ANISO", "4")):
        for root in ("1", "2", "3"):
            for suffix in ("N2_NULLITY", "N3_TRANSVERSE_NULLITY"):
                wl_tag = f"WL_S10_{package}_D{dimension}_ROOT{root}_{suffix}"
                py_tag = f"PY_S10_{package}_D{dimension}_ROOT{root}_{suffix}"
                print(f"{package} D{dimension} GENERIC ROOT{root} {suffix}")
                print(f"  wl = {wl_tags.get(wl_tag, '<no tag>')}")
                print(f"  py = {py_tags.get(py_tag, '<no tag>')}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
