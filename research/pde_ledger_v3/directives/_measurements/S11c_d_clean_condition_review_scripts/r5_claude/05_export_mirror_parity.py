#!/usr/bin/env python3
"""Mirror-parity census of the COMMITTED S11c-b exported operator objects on R1 x class P.

Reads the stored `slab_operator` and `mu_theta_operator` VALUE payloads from
research/pde_ledger_v3/scripts/S11c_b_exports.py (read-only; text extraction + the
exports' own restore namespace), then:
  R1      : every profile jet w1_profile_d.../m1_profile_d... carrying a 2 or 3 index -> 0
  class P : every perturbation/virtual/bulk jet carrying a d3 index -> 0
  mirror  : x3 -> -x3 acts on every symbol by (-1)^([component==3] + #d3)
For every scalar leaf L it prints the counts of EVEN / ODD / MIXED / ZERO leaves, where
MIXED means both (L+ML)/2 and (L-ML)/2 are nonzero (an odd<->even block inside one
output), and WRONG_SIDE counts a leaf whose path component index is 3 (odd output) but
whose content is purely even-nonzero, or vice versa, for the slab rows whose index
semantics are fixed by the engine (U rows, divergence fluxes, face U rows).
FORM control (CTRL_D3): R1 is NOT applied to the direction-3 first/higher profile jets;
they are kept as untransformed background data (not mirrored).  Same census.
Also prints which slab-field families occur in each top-level path (zeta_c presence).
"""
import re, sys, collections, time
import sympy as sp
from sympy.core.symbol import Str
from sympy.functions.elementary.piecewise import ExprCondPair

P = "/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_b_exports.py"
RELS = {
    'Equality': lambda l, r: sp.Eq(l, r, evaluate=False),
    'Unequality': lambda l, r: sp.Ne(l, r, evaluate=False),
    'StrictGreaterThan': lambda l, r: sp.Gt(l, r, evaluate=False),
    'StrictLessThan': lambda l, r: sp.Lt(l, r, evaluate=False),
    'GreaterThan': lambda l, r: sp.Ge(l, r, evaluate=False),
    'LessThan': lambda l, r: sp.Le(l, r, evaluate=False),
}
NS = {'__builtins__': {}, 'Str': Str, 'ExprCondPair': ExprCondPair, **vars(sp), **RELS}

def stored_value(name):
    with open(P, encoding="utf-8") as fh:
        hit = False
        for ln in fh:
            if ln.startswith(f"    '{name}':"):
                hit = True; continue
            if hit and ln.lstrip().startswith("'value': _restore("):
                s = ln.strip()[len("'value': _restore("):]
                s = s[:s.rfind(")")]
                return eval(eval(s), NS)
    raise KeyError(name)

def leaves(obj, path=()):
    if isinstance(obj, sp.Tuple):
        # labelled pair (Str, value)
        if len(obj) == 2 and isinstance(obj[0], Str):
            yield from leaves(obj[1], path + (str(obj[0]),))
            return
        for i, item in enumerate(obj):
            if isinstance(item, Str):
                continue
            yield from leaves(item, path + (i,))
        return
    if isinstance(obj, sp.MatrixBase):
        for i, item in enumerate(obj):
            yield path + (i,), item
        return
    if isinstance(obj, Str):
        return
    if isinstance(obj, sp.Expr):
        yield path, obj

D3 = re.compile(r"d(\d)")
def derivative_indices(name):
    # suffix after the field stem: e.g. u_3_d1d2 -> ['1','2'];  w1_profile_d2d3 -> ['2','3']
    m = re.search(r"_d(\d(?:d\d)*)$", name)
    if not m:
        return []
    return m.group(1).split("d")

def classify_symbol(name):
    """Return (component3:bool, d3count:int, family) for a perturbation-type symbol."""
    idx = derivative_indices(name)
    d3 = sum(1 for i in idx if i == '3')
    stem = re.sub(r"_d\d(?:d\d)*$", "", name)
    comp3 = False
    for pat in (r"^u_3(_t|_tt)?$", r"^delta_v_u_3$", r"^delta_v_bulk_(plus|minus)_3$",
                r"^d_w_delta_v_bulk_(plus|minus)_3$", r"^delta_j_bulk_3$",
                r"^d_w_delta_j_bulk_(plus|minus)_3$", r"^trace_grad_f_3$", r"^d_w_trace_grad_f_3$",
                r"^grad_theta_3$"):
        if re.match(pat, stem) or re.match(pat, name):
            comp3 = True
    # grad_theta_3 / e_W_d3 style first jets: the label is a derivative index, not a component
    if name == "grad_theta_3":
        comp3, d3 = False, 1
    return comp3, d3, stem

PROFILE = re.compile(r"^(w1|m1)_profile_d")
def build_maps(symbols, keep_d3_profile):
    r1p, mirror = {}, {}
    for s in symbols:
        n = s.name
        if PROFILE.match(n):
            idx = derivative_indices(n)
            if any(i in ('2', '3') for i in idx):
                if keep_d3_profile and '3' in idx and '2' not in idx and set(idx) <= {'1','3'}:
                    continue  # untransformed background datum (control): kept, not mirrored
                r1p[s] = 0
            continue
        comp3, d3, stem = classify_symbol(n)
        if d3 > 0:
            r1p[s] = 0          # class P: no x3 dependence of any perturbation/test/bulk field
            continue
        if comp3:
            mirror[s] = -s
    return r1p, mirror

def census(obj, keep_d3_profile, label):
    syms = set()
    for _, leaf in leaves(obj):
        syms |= leaf.free_symbols
    r1p, mirror = build_maps(syms, keep_d3_profile)
    counts = collections.Counter()
    wrong = 0
    mixed_paths = []
    for path, leaf in leaves(obj):
        L = sp.expand(leaf.xreplace(r1p))
        if L == 0:
            counts['ZERO'] += 1; continue
        ML = sp.expand(L.xreplace(mirror))
        E = sp.expand(L + ML); O = sp.expand(L - ML)
        if E != 0 and O != 0:
            counts['MIXED'] += 1; mixed_paths.append(path)
        elif E != 0:
            counts['EVEN'] += 1
            kind = 'EVEN'
        else:
            counts['ODD'] += 1
            kind = 'ODD'
        if E == 0 or O == 0:
            # index-semantics check for rows whose last integer index is a component label
            sp_path = [p for p in path if isinstance(p, str)]
            ints = [p for p in path if isinstance(p, int)]
            if sp_path and sp_path[-1] in ('LOCAL', 'EXPANDED') and 'U_BODY_BALANCE' in sp_path and ints:
                expect = 'ODD' if ints[-1] == 2 else 'EVEN'
                if expect != kind:
                    wrong += 1
            if sp_path and sp_path[-1] == 'DIVERGENCE_FLUX' and 'U_BODY_BALANCE' in sp_path and len(ints) >= 2:
                a, i = ints[-2], ints[-1]
                expect = 'ODD' if ((a == 2) ^ (i == 2)) else 'EVEN'
                if expect != kind:
                    wrong += 1
    print(label, dict(sorted(counts.items())), "WRONG_SIDE", wrong,
          "N_MIRRORED_SYMBOLS", len(mirror), "N_R1P_ZEROED", len(r1p))
    for p in mixed_paths[:6]:
        print("   MIXED_PATH", p)

t0 = time.time()
slab = stored_value("slab_operator")
mu = stored_value("mu_theta_operator")
print("LOADED_SECONDS", round(time.time() - t0, 1))
FAMILIES = {
    'zeta_c': re.compile(r"^(delta_v_)?zeta_c"),
    'u_3': re.compile(r"^u_3"), 'delta_v_u_3': re.compile(r"^delta_v_u_3"),
    'bulk_trace': re.compile(r"^(d_w_)?delta_(p|v_bulk|rho_4D_face|j_bulk)"),
}
for obj_name, obj in (("SLAB_OPERATOR", slab), ("MU_THETA_OPERATOR", mu)):
    for case in obj:
        key = tuple(str(k) for k in case[0])
        payload = case[1]
        value = None
        for item in payload:
            if str(item[0]) == "VALUE":
                value = item[1]
        census(value, False, f"{obj_name} {key} BASELINE")
        census(value, True, f"{obj_name} {key} CTRL_D3")
        if obj_name == "SLAB_OPERATOR":
            for top in value:
                name = str(top[0])
                fs = set()
                for _, lf in leaves(top[1]):
                    fs |= {s.name for s in lf.free_symbols}
                present = {f: any(rx.match(n) for n in fs) for f, rx in FAMILIES.items()}
                print("   PATH_FAMILIES", key, name, present)
print("TOTAL_SECONDS", round(time.time() - t0, 1))
