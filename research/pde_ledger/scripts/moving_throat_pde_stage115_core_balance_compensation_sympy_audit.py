
#!/usr/bin/env python3
from __future__ import annotations
import sympy as sp

def banner(title: str) -> None:
    line = "=" * 88
    print("\n" + line)
    print(title)
    print(line)

def expect_zero(name: str, expr: sp.Expr) -> None:
    s = sp.simplify(sp.expand(expr))
    print(f"{name} = {s}")
    if s != 0:
        raise AssertionError(f"{name} is not zero")

banner("STAGE 115 — EXACT CORE-BALANCE COMPENSATION THEOREM")

K_s, K_q, lam = sp.symbols("K_s K_q lam", positive=True, real=True)
g_s, g_q = sp.symbols("g_s g_q", real=True)
kappa0, gamma0, z = sp.symbols("kappa0 gamma0 z", positive=True, real=True)

r_c = lam**2 / (K_s * K_q)
rho_c = g_s**2 / K_s
sigma_c = (K_s * g_q - lam * g_s)**2 / (K_s**2 * K_q * (1 + r_c))

balance_eq = sp.expand(rho_c - 4 * sigma_c)
print("rho_c - 4 sigma_c =")
sp.pprint(balance_eq)

gq_solutions = sp.solve(sp.Eq(balance_eq, 0), g_q)
print("\nExact coupling-balance solutions for g_q:")
for sol in gq_solutions:
    sp.pprint(sp.simplify(sol))

# Pick one branch; both give the same sigma_c.
gq_branch = sp.simplify(gq_solutions[0])

sigma_star = sp.simplify((g_s**2) / (4 * K_s))
expect_zero(
    "sigma_c on balance surface",
    sigma_c.subs(g_q, gq_branch) - sigma_star,
)

kappa0_can = sp.simplify((1 + r_c) / 3)
gamma0_can = sp.simplify((1 + r_c) / 9)

delta_core = sp.simplify(
    rho_c
    - sigma_c / (1 - (kappa0 / (1 + r_c)) * z**2 - sp.I * (gamma0 / (1 + r_c)) * z**5)
).subs({
    g_q: gq_branch,
    kappa0: kappa0_can,
    gamma0: gamma0_can,
})

target_delta = sp.simplify(
    4 * sigma_star - sigma_star / (1 - z**2 / 3 - sp.I * z**5 / 9)
)
expect_zero("exact collapse to compensated branch", delta_core - target_delta)

Lambda_out = -3 + z**2 / 3 + z**4 / 9 + sp.I * z**5 / 9
Lambda_eff = sp.expand(sp.series(Lambda_out + target_delta, z, 0, 6).removeO())
Y_eff = sp.simplify(Lambda_eff.subs(z, 0) / Lambda_eff)
Y_target = 1 + z**2 / 9 + 4 * z**4 / 81 + sp.I * z**5 / 27
expect_zero(
    "normalized outgoing fingerprint preserved",
    sp.expand(sp.series(Y_eff, z, 0, 6).removeO()) - Y_target,
)

print("\nDerived branch data:")
print("sigma_* =", sigma_star)
print("kappa0  =", kappa0_can)
print("gamma0  =", gamma0_can)
