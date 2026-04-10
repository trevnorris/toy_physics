#!/usr/bin/env python3
"""
5pn_stage98_core_balance_compensation_theorem.py

Stage 98 audit: exact core-balance compensation theorem.
"""

from __future__ import annotations
import sympy as sp

def banner(title: str) -> None:
    line = "=" * 88
    print("\n" + line)
    print(title)
    print(line)

def expect_zero(name: str, expr) -> None:
    expr_s = sp.simplify(sp.factor(sp.together(sp.expand(expr))))
    print(f"{name} = {expr_s}")
    if expr_s != 0:
        raise AssertionError(f"{name} is not zero")

banner("STAGE 98 — EXACT CORE-BALANCE COMPENSATION THEOREM")

Ks, Kq, lam = sp.symbols("K_s K_q lambda", positive=True, real=True)
gs, gq = sp.symbols("g_s g_q", real=True)
kappa0, gamma0, z = sp.symbols("kappa_0 gamma_0 z", real=True)

rc = sp.simplify(lam**2 / (Ks * Kq))
rho_c = sp.simplify(gs**2 / Ks)
sigma_c = sp.simplify((Ks * gq - lam * gs)**2 / (Ks**2 * Kq * (1 + rc)))
kappa_c = sp.simplify(kappa0 / (1 + rc))
gamma_c = sp.simplify(gamma0 / (1 + rc))

balance_eq = sp.simplify(sp.factor(sp.together(rho_c - 4 * sigma_c)))
print("rho_c - 4 sigma_c =", balance_eq)
expect_zero(
    "balance law numerator",
    sp.factor(sp.together(balance_eq)).as_numer_denom()[0]
    - (gs**2 * (Ks * Kq + lam**2) - 4 * (Ks * gq - lam * gs)**2),
)

sol_gq = sp.solve(sp.Eq(gs**2 * (Ks * Kq + lam**2), 4 * (Ks * gq - lam * gs)**2), gq)
print("g_q branches =", sol_gq)
for idx, branch in enumerate(sol_gq, start=1):
    print(f"g_q branch {idx} =", sp.factor(branch))

expect_zero("kappa_0 condition", sp.solve(sp.Eq(kappa_c, sp.Rational(1, 3)), kappa0)[0] - (1 + rc) / 3)
expect_zero("gamma_0 condition", sp.solve(sp.Eq(gamma_c, sp.Rational(1, 9)), gamma0)[0] - (1 + rc) / 9)

Dw_can = 1 - z**2 / 3 - sp.I * z**5 / 9
sigma_star = sp.simplify(gs**2 / (4 * Ks))

delta_core = sp.simplify(
    rho_c - sigma_c / (1 - kappa_c * z**2 - sp.I * gamma_c * z**5)
)

subs_common = {
    kappa0: (1 + rc) / 3,
    gamma0: (1 + rc) / 9,
}
target = sp.simplify(4 * sigma_star - sigma_star / Dw_can)

for idx, branch in enumerate(sol_gq, start=1):
    delta_branch = sp.simplify(delta_core.subs(subs_common).subs(gq, branch))
    expect_zero(f"branch {idx} collapse to compensated outlet", delta_branch - target)

Lambda_out = -3 + z**2 / 3 + z**4 / 9 + sp.I * z**5 / 9
Lambda_tot = sp.series(Lambda_out + target, z, 0, 6).removeO()
Yhat_tot = sp.series(sp.simplify(sp.expand(Lambda_tot).coeff(z, 0) / Lambda_tot), z, 0, 6).removeO()
Yhat_can = 1 + z**2 / 9 + 4 * z**4 / 81 + sp.I * z**5 / 27
expect_zero("normalized response remains canonical", Yhat_tot - Yhat_can)

banner("STAGE 98 FINAL LEDGER")
print("The concrete two-channel core lands on the compensated canonical branch exactly when:")
print("  g_s^2 (K_s K_q + lambda^2) = 4 (K_s g_q - lambda g_s)^2,")
print("  kappa_0 = (1+r_c)/3,   gamma_0 = (1+r_c)/9,")
print("with r_c = lambda^2/(K_s K_q).")
print("On either coupling-balance branch, the eliminated core outlet collapses to")
print("  delta Lambda_core = 4 sigma_* - sigma_*/(1 - z^2/3 - i z^5/9),")
print("sigma_* = g_s^2/(4 K_s),")
print("and the normalized outgoing fingerprint remains exactly canonical.")
