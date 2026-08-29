#!/usr/bin/env python3
"""Reporting tool (NOT a physics computation): read the committed multigrade
instrument stdout and print, per case, the residual's nonzero (a,b) grade
support per leaf path and the remainder leading-grade set.  It only reads the
already-reviewed instrument's emitted objects; it computes no physics."""
import re
import sys
from collections import OrderedDict
import sympy as sp

RUN = sys.argv[1] if len(sys.argv) > 1 else "/home/trevnorris/.s11_build/S11c_b_multigrade_run.out"


def Association(*a, **k):
    entries = k["entries"] if "entries" in k else (a[0] if a else {})
    if isinstance(entries, dict):
        return OrderedDict(entries)
    return OrderedDict(entries)


def TextAtom(*a, **k):
    return str(k["value"]) if "value" in k else (str(a[0]) if a else "")


NS = {name: getattr(sp, name) for name in dir(sp) if not name.startswith("_")}
NS.update({"Association": Association, "TextAtom": TextAtom, "nan": sp.nan, "oo": sp.oo})


def ev(rhs):
    return eval(rhs, {"__builtins__": {}}, NS)


def case_label(assoc):
    fam = assoc["FAMILY"]
    key = assoc["KEY"]
    parts = ",".join(f"{k}={v}" for k, v in key.items())
    return f"{fam} ({parts})"


def is_zero(x):
    # The instrument already emits normalised coefficients, so a zero cell is
    # structurally Integer(0).  Structural test only — no simplify (fast).
    return x == 0


def grade_support(path_assoc):
    """path_assoc: OrderedDict grade-key 'a,b' -> coeff, plus REMAINDER / REMAINDER_LEADING_GRADES."""
    support = []
    leading = None
    for k, v in path_assoc.items():
        if k == "REMAINDER":
            continue
        if k == "REMAINDER_LEADING_GRADES":
            leading = v
            continue
        if re.fullmatch(r"\d+,\d+", str(k)) and not is_zero(v):
            support.append(str(k))
    support.sort(key=lambda s: tuple(int(t) for t in s.split(",")))
    return support, leading


def main():
    lines = open(RUN).read().splitlines()
    blocks = []
    cur = None
    for ln in lines:
        m = re.match(r"^([A-Z_]+) = (.*)$", ln)
        if not m:
            continue
        tag, rhs = m.group(1), m.group(2)
        if tag == "CASE":
            cur = {"CASE": rhs}
            blocks.append(cur)
        elif cur is not None and tag in ("MULTIGRADE_A", "MULTIGRADE_B", "MULTIGRADE_RESIDUAL"):
            cur[tag] = rhs
    print(f"CASES={len(blocks)}")
    for b in blocks:
        assoc = ev(b["CASE"])
        label = case_label(assoc)
        A = ev(b["MULTIGRADE_A"])
        B = ev(b["MULTIGRADE_B"])
        res = ev(b["MULTIGRADE_RESIDUAL"])
        print(f"\n=== {label}")
        for path in res:
            sa, la = grade_support(A[path])
            sb, lb = grade_support(B[path])
            sr, lr = grade_support(res[path])
            print(f"  leaf[{path}]")
            print(f"    A(PY) ={sa}  remA_lead={la}")
            print(f"    B(WL) ={sb}  remB_lead={lb}")
            print(f"    RES   ={sr}  remRES_lead={lr}")


if __name__ == "__main__":
    main()
