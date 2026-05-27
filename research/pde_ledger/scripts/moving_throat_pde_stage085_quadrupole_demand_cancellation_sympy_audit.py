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


banner("STAGE 085 — EXACT CANCELLATION OF OUTGOING-NORMALIZATION FACTORS")

pi = sp.pi
A, beta0, NQ = sp.symbols("A beta0 NQ", positive=True, real=True)
alpha_req, alpha_mix = sp.symbols("alpha_req alpha_mix", positive=True, real=True)
kappa0_sq = sp.Rational(8) / pi**2

R_target = sp.simplify(NQ * A / (beta0 * kappa0_sq))
G_tr = sp.simplify(8 * alpha_req / (pi**2 * A))
M_mix = sp.simplify(8 * alpha_mix / (pi**2 * A))

Pi_tr = sp.simplify(R_target * G_tr)
C_mix = sp.simplify(R_target * M_mix)

print("R_target =", R_target)
print("G_tr =", G_tr)
print("M_mix =", M_mix)
print("Pi_tr =", Pi_tr)
print("C_mix =", C_mix)

expect_zero("Pi_tr - NQ*alpha_req/beta0", Pi_tr - NQ * alpha_req / beta0)
expect_zero("C_mix - NQ*alpha_mix/beta0", C_mix - NQ * alpha_mix / beta0)
expect_zero("Pi_tr/C_mix - alpha_req/alpha_mix", sp.simplify(Pi_tr / C_mix - alpha_req / alpha_mix))

mhat, s_minus, lam_minus = sp.symbols("mhat s_minus lam_minus", positive=True, real=True)
NQ_selected = sp.simplify(mhat**2 * beta0 * s_minus / lam_minus)
Pi_sel = sp.simplify(Pi_tr.subs(NQ, NQ_selected))
C_sel = sp.simplify(C_mix.subs(NQ, NQ_selected))

print("NQ_selected =", NQ_selected)
print("Pi_sel =", Pi_sel)
print("C_sel =", C_sel)

expect_zero("Pi_sel - mhat^2 s_- alpha_req/lambda_-", Pi_sel - mhat**2 * s_minus * alpha_req / lam_minus)
expect_zero("C_sel - mhat^2 s_- alpha_mix/lambda_-", C_sel - mhat**2 * s_minus * alpha_mix / lam_minus)

# Exact selected-demand law in pure loading-ratio form.
eps_blk = sp.symbols("eps_blk", real=True)
rho_alpha = sp.symbols("rho_alpha", positive=True, real=True)
zeta_req = sp.simplify((Pi_tr - C_mix) / (C_mix - eps_blk * (2 * C_mix - Pi_tr)))
zeta_expected = sp.simplify((alpha_req - alpha_mix) / (alpha_mix - eps_blk * (2 * alpha_mix - alpha_req)))
expect_zero("zeta_req product form - loading form", zeta_req - zeta_expected)
expect_zero(
    "zeta_req loading form - rho_alpha form",
    zeta_expected.subs(alpha_req, rho_alpha * alpha_mix) - (rho_alpha - 1) / (1 - eps_blk * (2 - rho_alpha)),
)
expect_zero(
    "unblocked limit",
    sp.simplify(zeta_expected.subs(eps_blk, 0) - (alpha_req / alpha_mix - 1)),
)

print("\nFINAL LEDGER")
print("Pi_tr = (N_Q^(target)/beta0) alpha_req")
print("C_mix = (N_Q^(target)/beta0) alpha_mix")
print("Pi_tr/C_mix = alpha_req/alpha_mix")
print("zeta_req = (alpha_req-alpha_mix)/[alpha_mix - eps_blk(2 alpha_mix-alpha_req)]")
