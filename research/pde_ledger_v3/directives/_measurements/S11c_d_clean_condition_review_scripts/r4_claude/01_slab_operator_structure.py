#!/usr/bin/env python3
"""Leg script 01: load the COMMITTED S11c-b exported slab_operator payload (a copy of
one line of research/pde_ledger_v3/scripts/S11c_b_exports.py, read-only) and print,
for each (anchoring, density) case and each named row/sub-row, the free symbols that
are NOT in directive B's declared formal input domain (B §2.3: X_slab + 20 p/v trace
coordinates) and NOT background/parameter symbols.

Purpose: test whether B §2.3's 20-coordinate trace domain is the complete set of
bulk-trace coordinates the engine's rows actually depend on (R3-1 resolution).
Prints computed sets only.
"""
import sys, re, time
import sympy as sp
from sympy.core.symbol import Str
from sympy.functions.elementary.piecewise import ExprCondPair

sys.setrecursionlimit(100000)
LINE = "/tmp/s11cd_clean_review_r4_claude/engine_copy/slab_operator_value_line.txt"
src = open(LINE).read().strip()
m = re.match(r"^\s*'value': _restore\((\".*\")\),\s*$", src, re.S)
code = eval(m.group(1))  # the string literal
_REL = {
    'Equality': lambda l, r: sp.Eq(l, r, evaluate=False),
    'Unequality': lambda l, r: sp.Ne(l, r, evaluate=False),
}
t0 = time.time()
value = eval(code, {'__builtins__': {}, 'Str': Str, 'ExprCondPair': ExprCondPair, **vars(sp), **_REL})
print("LOAD_SECONDS", round(time.time() - t0, 1))

# Directive B §2.3 declared formal domain (verbatim names).
B_SLAB = ["u_1", "u_2", "u_3", "theta", "e_W", "zeta_c"]
B_TRACE = []
for face in ("plus", "minus"):
    B_TRACE += [f"delta_p_{face}", f"d_w_delta_p_{face}"]
    B_TRACE += [f"delta_v_bulk_{face}_{i}" for i in range(1, 5)]
    B_TRACE += [f"d_w_delta_v_bulk_{face}_{i}" for i in range(1, 5)]
B_DOMAIN = set(B_SLAB + B_TRACE)

# A symbol is a wave/trace coordinate if it is a jet/time-derivative of a declared slab
# field (u_a_*, theta_*, grad_theta_*, e_W_*, zeta_c_*) or a declared trace name.
def is_declared(name):
    if name in B_DOMAIN:
        return True
    for base in ("u_1", "u_2", "u_3", "theta", "e_W", "zeta_c"):
        if name.startswith(base + "_d") or name in (base + "_t", base + "_tt") or name.startswith(base + "_t_d"):
            return True
    if name.startswith("grad_theta_"):
        return True
    return False

BACKGROUND_OR_PARAM = re.compile(
    r"^(B_rho_3|C|G_theta_u|G_W_u|kappa_theta|kappa_theta_W|kappa_W|k_W|Lambda_[AVX]_0|L_W|"
    r"m1_profile.*|w1_profile.*|mu_R|mu_S|mu_W|mu_theta_[LM]|omega|rho_br|rho_m|sigma_W|tau_[AVX]|"
    r"W_0|eta_bg|epsilon_shape|gamma_s11cb_.*|e_W_bg)$")

def undeclared(expr):
    names = sorted(str(s) for s in sp.sympify(expr).free_symbols)
    return [n for n in names if not is_declared(n) and not BACKGROUND_OR_PARAM.match(n)]

def walk(obj, path, out):
    if isinstance(obj, sp.Tuple) and len(obj) == 2 and isinstance(obj[0], Str):
        walk(obj[1], path + [str(obj[0])], out)
        return
    if isinstance(obj, (sp.Tuple, tuple, list)):
        for i, item in enumerate(obj):
            walk(item, path + [f"[{i}]"], out)
        return
    if isinstance(obj, sp.MatrixBase):
        for i, item in enumerate(obj):
            walk(item, path + [f"<{i}>"], out)
        return
    if isinstance(obj, Str):
        return
    try:
        bad = undeclared(obj)
    except Exception as exc:  # non-expression leaf
        return
    out.append(("/".join(path), bad))

for case in value:
    axes, payload = case[0], case[1]
    rows = []
    walk(payload, [], rows)
    print("CASE", axes)
    agg = {}
    for path, bad in rows:
        # collapse index-only path components after the row name for a compact report
        key = re.sub(r"/\[\d+\]", "/[i]", path)
        key = re.sub(r"/<\d+>", "/<i>", key)
        agg.setdefault(key, set()).update(bad)
    for key in sorted(agg):
        print("  ROW", key, "UNDECLARED_INPUT_SYMBOLS", sorted(agg[key]))
