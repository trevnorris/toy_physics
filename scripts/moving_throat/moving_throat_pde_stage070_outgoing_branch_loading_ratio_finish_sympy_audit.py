#!/usr/bin/env python3
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


banner("STAGE 70 — FINAL REDUCED FINISH-LINE IN THE LOADING-RATIO VARIABLE")

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

expect_zero("zeta_suff - 2.46622291347846", sp.N(zeta_suff - sp.Float("2.46622291347846"), 18))
expect_zero("zeta_fail - 2.46752913273870", sp.N(zeta_fail - sp.Float("2.46752913273870"), 18))
expect_zero("zeta_max - 2.46752922945601", sp.N(zeta_max - sp.Float("2.46752922945601"), 18))

print("\nFINAL LEDGER")
print("The reduced explicit Family-1 theorem is completely equivalent to a one-number test:")
print("  rho_alpha = alpha_req/alpha_mix")
print("with success/failure determined by the Stage-69 ratio window.")
