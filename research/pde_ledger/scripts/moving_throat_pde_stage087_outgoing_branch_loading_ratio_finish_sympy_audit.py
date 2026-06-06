#!/usr/bin/env python3
"""
SymPy audit for Stage 087.

This stage is a checkpoint-consolidation statement, not a fresh derivation.
The paper card (paper/stages/stage_087.tex) `Purpose` line reads:
"Stage 087 records that the explicit Family-1 support/source side has been
reduced to a single outgoing-branch loading ratio." Its `Inputs` are stages
085 and 086. The cancellation chain that collapses dependence on
s_-, lambda_-, beta_0, mhat_-, Pi_tr, C_mix, Pe_req down to rho_alpha alone
is performed and verified upstream in stages 081-086 (post-renumber):
- scripts/moving_throat_pde_stage081_*_sympy_audit.py / .wl  (a.k.a. former stage 65)
- scripts/moving_throat_pde_stage082_*_sympy_audit.py / .wl  (former stage 69 closure)
- scripts/moving_throat_pde_stage085_quadrupole_demand_cancellation_*  (Pi_tr/C_mix cancellation)
- scripts/moving_throat_pde_stage086_family1_loading_ratio_window_*  (Family-1 window)

This script restates the unblocked one-ratio criterion `zeta_req(rho_alpha; 0) = rho_alpha - 1`
as a downstream-consistency probe and sanity-checks the Family-1 window
literals carried from the Stage-086 notes via their ordering and
constructive-ceiling gap.
"""

from __future__ import annotations
import sympy as sp


def banner(t: str) -> None:
    print("\n" + "="*88)
    print(t)
    print("="*88)


def expect_zero(name: str, expr: sp.Expr) -> None:
    expr = sp.simplify(sp.expand(expr))
    print(f"{name} = {expr}")
    if expr != 0:
        raise AssertionError(f"{name} is not zero")


def expect_close(name: str, value: sp.Expr, target: sp.Expr, tol: sp.Expr) -> None:
    diff = sp.Abs(sp.N(value - target, 40))
    print(f"{name} diff = {diff}")
    if diff > tol:
        raise AssertionError(f"{name} exceeds tolerance {tol}: diff={diff}")


banner("STAGE 087 — FINAL REDUCED FINISH-LINE IN THE LOADING-RATIO VARIABLE")

rho_alpha, eps_blk = sp.symbols("rho_alpha eps_blk", positive=True, real=True)
zeta_req = sp.simplify((rho_alpha - 1) / (1 - eps_blk * (2 - rho_alpha)))
print("zeta_req(rho_alpha,eps_blk) =", zeta_req)
print("d zeta_req / d rho_alpha =", sp.simplify(sp.diff(zeta_req, rho_alpha)))

# Unblocked reduction
expect_zero("unblocked zeta_req", sp.simplify(zeta_req.subs(eps_blk, 0) - (rho_alpha - 1)))

# Family-1 unblocked window values.
rho_suff = sp.Float("3.46622291347846")
rho_fail = sp.Float("3.46752913273870")
rho_max = sp.Float("3.46752922945601")
zeta_suff = sp.N(zeta_req.subs({rho_alpha: rho_suff, eps_blk: 0}), 25)
zeta_fail = sp.N(zeta_req.subs({rho_alpha: rho_fail, eps_blk: 0}), 25)
zeta_max = sp.N(zeta_req.subs({rho_alpha: rho_max, eps_blk: 0}), 25)

print("zeta at success ratio =", zeta_suff)
print("zeta at failure ratio =", zeta_fail)
print("zeta at max ratio     =", zeta_max)

# The Family-1 ratio-window literals below are carried from the notes.
# Their strict ordering and tight constructive-ceiling gap sanity-check
# transcription of the three carried values.
print("rho_suff^(chi) =", rho_suff)
print("rho_fail^(chi) =", rho_fail)
print("rho_max^(F1)   =", rho_max)
rho_suff_before_fail = bool(rho_suff < rho_fail)
rho_fail_before_max = bool(rho_fail < rho_max)
rho_gap = sp.N(rho_max - rho_fail, 40)
rho_gap_ok = bool(rho_gap > 0) and bool(rho_gap < sp.Float("1e-6"))
print("rho_suff < rho_fail =", rho_suff_before_fail)
if not rho_suff_before_fail:
    raise AssertionError("rho_suff is not strictly below rho_fail")
print("rho_fail < rho_max =", rho_fail_before_max)
if not rho_fail_before_max:
    raise AssertionError("rho_fail is not strictly below rho_max")
print("0 < rho_max - rho_fail < 1e-6 =", rho_gap_ok, "; gap =", rho_gap)
if not rho_gap_ok:
    raise AssertionError("rho_max - rho_fail is not in the tight constructive-ceiling gap")

print("\nFINAL LEDGER")
print("The reduced explicit Family-1 theorem is completely equivalent to a one-number test:")
print("  rho_alpha = alpha_req/alpha_mix")
print("with success/failure determined by the Stage-69 ratio window.")
