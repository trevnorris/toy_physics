#!/usr/bin/env python3
"""
5pn_stage65_master_quadrupole_residual.py

Stage 65 audit: master quadrupole residual of the full reduced moving-throat PDE.
"""

from __future__ import annotations

import math
import mpmath as mp
import sympy as sp


def banner(title: str) -> None:
    line = "=" * 88
    print("\n" + line)
    print(title)
    print(line)


def subbanner(title: str) -> None:
    line = "-" * 88
    print("\n" + line)
    print(title)
    print(line)


def expect_zero(name: str, expr) -> None:
    expr_s = sp.simplify(sp.together(sp.expand(expr)))
    print(f"{name} = {expr_s}")
    if expr_s != 0:
        raise AssertionError(f"{name} is not zero")

banner("STAGE 65 — MASTER QUADRUPOLE RESIDUAL")

Pi_tr, C_mix, eps_blk = sp.symbols("Pi_tr C_mix eps_blk", positive=True, real=True)
Xi, eta, kappa = sp.symbols("Xi eta kappa", positive=True, real=True)
Pe_star = sp.Function("Pe_*")(Xi, eta, kappa)
zeta_req = sp.simplify((Pi_tr - C_mix) / (C_mix - eps_blk * (2 * C_mix - Pi_tr)))
print("zeta_req(Pi_tr,C_mix,eps_blk) =")
sp.pprint(zeta_req)

subbanner("65.1 — Exact physical support ratio")
Pe = sp.symbols("Pe", nonnegative=True, real=True)
Omega_Pe = sp.Function("Omega_Pe")(Pe)
yeta = sp.Function("y")(eta)
zeta_phys = Omega_Pe**2 * (kappa + sp.pi**2 / 4) / (kappa + yeta**2)
print("zeta_phys(Pe,eta;kappa) =")
sp.pprint(zeta_phys)

subbanner("65.2 — Master residual")
R_quad = sp.simplify(zeta_req - zeta_phys.subs(Pe, Pe_star))
print("R_quad =")
sp.pprint(R_quad)

subbanner("65.3 — Exact bounded version")
Delta0, Deltainf = sp.symbols("Delta_0 Delta_inf", positive=True, real=True)
zeta_minus = sp.Function("zeta_-")(Xi, eta, kappa)
zeta_plus = sp.Function("zeta_+")(Xi, eta, kappa)
R_minus = sp.Symbol("R_-", real=True)
R_plus = sp.Symbol("R_+", real=True)
print("If Xi Delta_0 <= Pe_* <= Xi Delta_inf and zeta_phys is monotone increasing in Pe, then")
print("  zeta_- <= zeta_phys(Pe_*) <= zeta_+")
print("and therefore")
print("  R_- := zeta_req - zeta_+ <= R_quad <= zeta_req - zeta_- =: R_+.")
print("Hence")
print("  R_+ <= 0  -> guaranteed success,")
print("  R_- > 0   -> guaranteed failure.")

subbanner("65.4 — Direct product thresholds")
zeta_bound = sp.symbols("zeta_bound", positive=True, real=True)
Q = sp.simplify((1 + (1 - 2 * eps_blk) * zeta_bound) / (1 - eps_blk * zeta_bound))
Pi_bound = sp.simplify(C_mix * Q)
print("Pi_bound(C_mix,eps_blk;zeta_bound) =")
sp.pprint(Pi_bound)

subbanner("65.5 — Family-1 specialization")
Upsilon_w, Theta_w = sp.symbols("Upsilon_w Theta_w", positive=True, real=True)
Xi_F1 = sp.simplify(1369 * Upsilon_w)
Xi_F1_theta = sp.simplify(136900 * Theta_w)
print("Xi_F1 =")
sp.pprint(Xi_F1)
print("Xi_F1 =")
sp.pprint(Xi_F1_theta)

banner("STAGE 65 FINAL LEDGER")
print("Stage 65 compresses the whole reduced moving-throat PDE to the single scalar residual")
print("  R_quad = zeta_req - zeta_phys(Pe_*).")
print("Everything else in the reduced program now feeds this one object.")
