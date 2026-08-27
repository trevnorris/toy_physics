#!/usr/bin/env python3
"""
S11c-a T7 HAND-CODED cross-engine comparison.

This is the reconciliation layer, written out BY HAND. It does NOT auto-classify residuals and does NOT
blanket-collapse applied functions to bare symbols — that blanket collapse deletes real dependence and can
only ever HIDE a genuine disagreement (it hid the in-plane bulk-current w-freeze; see
memory feedback_handcode_comparison_never_blanket_collapse). Instead:

  1. A basic mechanical pass clears the obvious matches: residuals that are literally zero, and the handful of
     PURE VARIABLE RENAMES below that were read and confirmed to be the same variable.
  2. A tiny, explicitly VERIFIED collapse list handles the applied-vs-bare convention where — and only where —
     it was checked that WL never differentiates the object (so its arguments are an inert evaluation point):
       - mu_theta_{L,M}: verified 0 x-jets / 0 Derivative in the whole run.
       - delta_j_bulk_*  : the current jets. WL writes them f(x1,x2,x3,w,time); the w is REAL (§1b) and PY now
         carries it too (Function(name)(w)); the x1,x2,x3,time are the inert convention (verified 0 x-jets).
         So we KEEP w and map only the implicit coordinates away — never strip w.
  3. Everything else is compared AS-IS. Any object the two engines write differently that is not covered above
     is FLAGGED (a candidate finding), never massaged.

Run:  python3 S11c_a_handcoded_comparison.py            # committed transcripts
      python3 S11c_a_handcoded_comparison.py --py <path>  --wl <path>
"""
from __future__ import annotations
import argparse, io, contextlib
from pathlib import Path
import sympy as sp
from sympy import Symbol, Function, Integer, Rational, Mul, Add, Pow, I, Float, Dummy
from sympy.core.function import AppliedUndef
import S11c_a_cross_engine_comparator as C

NS = dict(Symbol=Symbol, Function=Function, Integer=Integer, Rational=Rational, Mul=Mul, Add=Add,
          Pow=Pow, I=I, Float=Float, Dummy=Dummy, Tuple=sp.Tuple, oo=sp.oo, pi=sp.pi)

# ---- HARDCODED, HAND-VERIFIED reconciliation --------------------------------------------------------
WL_TO_PY_RENAME = {            # same variable, different spelling (read + confirmed)
    "capitalX1": "X_1", "capitalX2": "X_2", "capitalX3": "X_3",
    "widthProfile": "W_bg",
    "w1ProfileXJETXd1": "w1_profile_d1", "w1ProfileXJETXd2": "w1_profile_d2", "w1ProfileXJETXd3": "w1_profile_d3",
}
INERT_APPLIED = {"mu_theta_L", "mu_theta_M"}   # WL writes f(x1,x2,x3,time); verified WL never differentiates
                                               # them (0 x-jets / 0 Derivative) -> args inert -> == PY bare.
CURRENT_KEEP_W = "delta_j_bulk"                # current jets: keep the REAL w, map implicit x,t to f(w).
# -----------------------------------------------------------------------------------------------------

W = C.BOUND_BINDER  # the canonical w binder the comparator uses inside bound integrals

def parse(raw: str):
    raw = raw.strip()
    if raw.startswith(("TextAtom", "Mismatch")) or raw == "<MISSING>":
        return None
    v = eval(raw, {"__builtins__": {}}, NS)
    return Add(*v) if isinstance(v, tuple) else v

def reconcile(expr):
    if expr is None or not hasattr(expr, "xreplace"):
        return expr
    subs = {}
    for a in expr.atoms(AppliedUndef):
        nm = a.func.__name__
        if nm in INERT_APPLIED:
            subs[a] = Symbol(nm)
        elif nm in WL_TO_PY_RENAME:
            subs[a] = Function(WL_TO_PY_RENAME[nm])(*a.args)
        elif nm.startswith(CURRENT_KEEP_W) and W in a.free_symbols and a.args != (W,):
            subs[a] = Function(nm)(W)          # keep w, drop the implicit coordinate args
    for s in expr.atoms(Symbol):
        if s.name in WL_TO_PY_RENAME:
            subs[s] = Symbol(WL_TO_PY_RENAME[s.name])
    return expr.xreplace(subs)

def residual_zero(a_raw: str, b_raw: str):
    pa, pb = parse(a_raw), parse(b_raw)
    if pa is None or pb is None:
        return None  # sentinel / coverage — not a comparison
    ra, rb = reconcile(pa), reconcile(pb)
    if isinstance(ra, tuple) or isinstance(rb, tuple):
        return str(ra) == str(rb)
    try:
        return sp.simplify(C.combine_bound_integrals(sp.expand(Add(ra, Mul(-1, rb))))) == 0
    except Exception:
        return False

def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--py", type=Path, default=C.DEFAULT_PY)
    ap.add_argument("--wl", type=Path, default=C.DEFAULT_WL)
    args = ap.parse_args(argv)
    py, wl = C.load_py(args.py), C.load_wl(args.wl)
    print(f"{'FAMILY':32s} VERDICT")
    print("-" * 60)
    for family in C.FAMILIES:
        cases = C.extract_family(family, py.get(family), wl.get(family))
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            C.compare_family(family, cases, leaf_budget_seconds=0.1)
        lines = buf.getvalue().splitlines()
        A = [l[len("operand_A = "):] for l in lines if l.startswith("operand_A = ")]
        B = [l[len("operand_B = "):] for l in lines if l.startswith("operand_B = ")]
        verdicts = [residual_zero(a, b) for a, b in zip(A, B)]
        if not verdicts:
            tag = "EMPTY"
        elif all(v is None for v in verdicts):
            tag = "SKIP  (coverage — one engine only; see ACCOUNTING)"
        elif all(v for v in verdicts if v is not None):
            tag = "MATCH"
        else:
            n = sum(1 for v in verdicts if v is False)
            tag = f"FLAG  ({n}/{len(verdicts)} cases differ — read them)"
        print(f"{family:32s} {tag}")

if __name__ == "__main__":
    import sys
    sys.exit(main())
