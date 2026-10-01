#!/usr/bin/env python3
"""Leg script 02: reflection-parity decomposition of the COMMITTED S11c-b exported
slab_operator payload (read-only copy of one line of S11c_b_exports.py), on the
domain R1 (background profiles depend on x_1 only: every profile jet carrying a
derivative label 2 or 3 set to zero) and class P (every wave/trial/test field
independent of x_3: every field jet carrying derivative label 3 set to zero).

The mirror M: x_3 -> -x_3 acts on a symbol by the sign (-1)^(number of '3' component
labels) (derivative-3 labels are already removed by P).  For every leaf expression e
of every row, print the number of terms of its even part (e + M e)/2 and of its odd
part (e - M e)/2.  A leaf with BOTH parts nonzero is a reflection-odd <-> even block.

Control (FORM, re-enters at the data, not at a result): the same decomposition with
R1 lifted for ONE background datum only -- w1_profile_d3 (and m1_profile_d3) kept
live (a direction-3 first jet), everything else identical.  If the decomposition
machinery can see a block, this control must produce mixed leaves.

Prints computed counts only; no conclusion.
"""
import sys, re, time
import sympy as sp
from sympy.core.symbol import Str
from sympy.functions.elementary.piecewise import ExprCondPair

sys.setrecursionlimit(100000)
LINE = "/tmp/s11cd_clean_review_r4_claude/engine_copy/slab_operator_value_line.txt"
src = open(LINE).read().strip()
m = re.match(r"^\s*'value': _restore\((\".*\")\),\s*$", src, re.S)
code = eval(m.group(1))
_REL = {'Equality': lambda l, r: sp.Eq(l, r, evaluate=False),
        'Unequality': lambda l, r: sp.Ne(l, r, evaluate=False)}
value = eval(code, {'__builtins__': {}, 'Str': Str, 'ExprCondPair': ExprCondPair, **vars(sp), **_REL})

ALL = set()
def collect(obj):
    if isinstance(obj, Str):
        return
    if isinstance(obj, (sp.Tuple, tuple, list)):
        for it in obj:
            collect(it)
        return
    if isinstance(obj, sp.MatrixBase):
        for it in obj:
            collect(it)
        return
    if isinstance(obj, sp.Basic):
        ALL.update(obj.free_symbols)
collect(value)

FIELD_BASES = ("u_1", "u_2", "u_3", "theta", "e_W", "zeta_c", "delta_v_u_1", "delta_v_u_2", "delta_v_u_3",
               "delta_j_bulk_1", "delta_j_bulk_2", "delta_j_bulk_3", "delta_j_bulk_4")
PROFILE = re.compile(r"^(w1_profile|m1_profile)_((?:d\d)+)$")
FIELDJET = re.compile(r"^(u_[123]|u_[123]_t|theta|e_W|zeta_c|delta_v_u_[123]|delta_j_bulk_[1234])_((?:d\d)+)$")

def component_sign(name):
    # component-3 labels of vector-valued coordinates
    for pat in (r"^u_3(_|$)", r"^delta_v_u_3(_|$)", r"^delta_v_bulk_(plus|minus)_3$",
                r"^d_w_delta_v_bulk_(plus|minus)_3$", r"^delta_j_bulk_3(_|$)",
                r"^d_w_delta_j_bulk_(plus|minus)_3$", r"^trace_grad_f_3$", r"^d_w_trace_grad_f_3$",
                r"^f_hold_u_3_0$", r"^t_hold_(plus|minus)_0_3$", r"^u_[TL]_3(_|$)"):
        if re.search(pat, name):
            return -1
    return 1

def build_maps(lift_direction3_profile_jet):
    zero = {}
    sign = {}
    for s in ALL:
        n = str(s)
        mp = PROFILE.match(n)
        if mp:
            idx = mp.group(2).replace("d", "")
            if lift_direction3_profile_jet and idx == "3":
                # control: a FIXED background datum along direction 3 that the mirror does
                # NOT transform (the background then violates R1); keep it live, sign +1
                continue
            if "2" in idx or "3" in idx:
                zero[s] = sp.Integer(0)
            continue
        mf = FIELDJET.match(n)
        if mf and "3" in mf.group(2).replace("d", ""):
            zero[s] = sp.Integer(0)
            continue
        if n in ("grad_theta_3",):
            zero[s] = sp.Integer(0)
            continue
        if component_sign(n) == -1:
            sign[s] = -s
    return zero, sign

def leaves(obj, path, out):
    if isinstance(obj, sp.Tuple) and len(obj) == 2 and isinstance(obj[0], Str):
        leaves(obj[1], path + [str(obj[0])], out); return
    if isinstance(obj, Str):
        return
    if isinstance(obj, (sp.Tuple, tuple, list)):
        for i, it in enumerate(obj):
            leaves(it, path + [f"[{i}]"], out)
        return
    if isinstance(obj, sp.MatrixBase):
        for i, it in enumerate(obj):
            leaves(it, path + [f"<{i}>"], out)
        return
    if isinstance(obj, sp.core.relational.Relational):
        out.append(("/".join(path) + "/REL_LHS_MINUS_RHS", obj.lhs - obj.rhs)); return
    if isinstance(obj, (sp.logic.boolalg.Boolean,)) and not isinstance(obj, sp.Expr):
        return
    if isinstance(obj, sp.Basic):
        out.append(("/".join(path), obj))

def nterms(e):
    e = sp.expand(e)
    return 0 if e == 0 else len(sp.Add.make_args(e))

for label, lift in (("BASELINE_R1_P", False), ("CONTROL_LIFT_PROFILE_JET_D3", True)):
    zero, sign = build_maps(lift)
    print("RUN", label, "ZEROED_SYMBOLS", len(zero), "MIRROR_ODD_SYMBOLS", sorted(str(k) for k in sign))
    for case in value:
        axes, payload = case[0], case[1]
        out = []
        leaves(payload, [], out)
        mixed = []
        classes = {}
        n_leaf = 0
        for path, e in out:
            if "DIMENSION_L_T_M" in path or "GRADE" in path:
                continue
            n_leaf += 1
            e0 = e.xreplace(zero)
            em = e0.xreplace(sign)
            ev = nterms((e0 + em) / 2)
            od = nterms((e0 - em) / 2)
            cls = "MIXED" if (ev and od) else ("PURE_EVEN" if ev else ("PURE_ODD" if od else "ZERO"))
            classes.setdefault(cls, []).append(path)
            if ev and od:
                mixed.append((path, ev, od))
        print(" ", label, "CASE", axes, "LEAVES", n_leaf, "MIXED_PARITY_LEAVES", len(mixed))
        for path, ev, od in mixed[:40]:
            print("    MIXED", path, "EVEN_TERMS", ev, "ODD_TERMS", od)
        if len(mixed) > 40:
            print("    MIXED_TRUNCATED_LISTING", len(mixed) - 40)
        for cls in sorted(classes):
            print("    CLASS_COUNT", cls, len(classes[cls]))
        if label == "BASELINE_R1_P":
            for path in classes.get("PURE_ODD", []):
                print("    PURE_ODD_LEAF", path)
