#!/usr/bin/env python3
"""Leg script 11: in the COMMITTED exported S11c-b slab_operator payload (read-only copy), list
every leaf (row path) whose free symbols include the physical centre-shift coordinates
zeta_c, zeta_c_t, zeta_c_d1..3, and every leaf containing e_W-type coordinates, per case.
Prints computed path lists only."""
import sys, re
import sympy as sp
from sympy.core.symbol import Str
from sympy.functions.elementary.piecewise import ExprCondPair
sys.setrecursionlimit(100000)
src = open("/tmp/s11cd_clean_review_r4_claude/engine_copy/slab_operator_value_line.txt").read().strip()
code = eval(re.match(r"^\s*'value': _restore\((\".*\")\),\s*$", src, re.S).group(1))
value = eval(code, {'__builtins__': {}, 'Str': Str, 'ExprCondPair': ExprCondPair, **vars(sp),
                    'Equality': lambda l, r: sp.Eq(l, r, evaluate=False)})
def leaves(obj, path, out):
    if isinstance(obj, sp.Tuple) and len(obj) == 2 and isinstance(obj[0], Str):
        leaves(obj[1], path + [str(obj[0])], out); return
    if isinstance(obj, Str): return
    if isinstance(obj, (sp.Tuple, tuple, list)):
        for i, it in enumerate(obj): leaves(it, path + [f"[{i}]"], out)
        return
    if isinstance(obj, sp.MatrixBase):
        for i, it in enumerate(obj): leaves(it, path + [f"<{i}>"], out)
        return
    if isinstance(obj, sp.core.relational.Relational):
        out.append(("/".join(path), obj.lhs - obj.rhs)); return
    if isinstance(obj, sp.Basic) and hasattr(obj, "free_symbols"):
        out.append(("/".join(path), obj))
ZC = re.compile(r"^zeta_c(_t|_d\d)?$")
for case in value:
    out = []; leaves(case[1], [], out)
    hits = sorted({re.sub(r"/\[\d+\]", "/[i]", p) for p, e in out
                   if "DIMENSION" not in p and "GRADE" not in p and any(ZC.match(str(s)) for s in e.free_symbols)})
    print("CASE", case[0], "LEAF_PATHS_WITH_PHYSICAL_ZETA_C", len(hits))
    for h in hits: print("   ", h)
