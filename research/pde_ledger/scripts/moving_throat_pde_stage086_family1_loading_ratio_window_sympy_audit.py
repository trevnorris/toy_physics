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


banner("STAGE 69 — FAMILY-1 PURE LOADING-RATIO WINDOW")

eps = sp.symbols("eps", real=True)
zeta = sp.symbols("zeta", positive=True, real=True)
Q = sp.simplify((1 + (1 - 2 * eps) * zeta) / (1 - eps * zeta))
print("Q(zeta;eps) =", Q)
print("dQ/dzeta =", sp.simplify(sp.diff(Q, zeta)))

# Exact unblocked reduction.
expect_zero("Q(zeta;0) - (1+zeta)", sp.simplify(Q.subs(eps, 0) - (1 + zeta)))

# Numerical Family-1 data carried from Stages 63-64.
zeta_suff_chi = sp.Float("2.46622291347846")
zeta_fail_chi = sp.Float("2.46752913273870")
zeta_suff_J = sp.Float("2.44257571477179")
zeta_max = sp.Float("2.46752922945601")

rho_suff_chi = sp.N(Q.subs({zeta: zeta_suff_chi, eps: 0}), 30)
rho_fail_chi = sp.N(Q.subs({zeta: zeta_fail_chi, eps: 0}), 30)
rho_suff_J = sp.N(Q.subs({zeta: zeta_suff_J, eps: 0}), 30)
rho_max = sp.N(Q.subs({zeta: zeta_max, eps: 0}), 30)

print("rho_suff^(chi)(eps=0) =", rho_suff_chi)
print("rho_fail^(chi)(eps=0) =", rho_fail_chi)
print("rho_suff^(J)(eps=0)   =", rho_suff_J)
print("rho_max^(F1)(eps=0)   =", rho_max)

expect_zero("rho_suff^(chi) - (1+zeta_suff)", sp.N(rho_suff_chi - (1 + zeta_suff_chi), 20))
expect_zero("rho_fail^(chi) - (1+zeta_fail)", sp.N(rho_fail_chi - (1 + zeta_fail_chi), 20))
expect_zero("rho_suff^(J) - (1+zeta_suff_J)", sp.N(rho_suff_J - (1 + zeta_suff_J), 20))
expect_zero("rho_max - (1+zeta_max)", sp.N(rho_max - (1 + zeta_max), 20))

# Exact blocking cap.
eps_cap = sp.N(1 / zeta_max, 25)
print("eps_blk cap = 1/zeta_max =", eps_cap)

Qmax_eps = sp.simplify(Q.subs(zeta, zeta_max))
print("rho_max^(F1)(eps) =", Qmax_eps)
print("d rho_max / d eps =", sp.simplify(sp.diff(Qmax_eps, eps)))

print("\nFINAL LEDGER")
print("rho_alpha <= Q(zeta_suff^(chi);eps_blk)  -> guaranteed success")
print("rho_alpha >= Q(zeta_fail^(chi);eps_blk)  -> guaranteed failure")
print("rho_alpha <  Q(zeta_max^(F1);eps_blk)    -> hard ceiling")
