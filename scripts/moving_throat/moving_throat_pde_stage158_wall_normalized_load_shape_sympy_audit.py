#!/usr/bin/env python3
"""
moving_throat_pde_stage158_wall_normalized_load_shape_sympy_audit.py

SymPy-backed audit for Stage 158.

Checks:
1. Exact wall-normalized factorization
      B0 = K chi^2,
      Z0 = K Upsilon,
      N0 = Lambda^2.
2. Exact defect rewrite
      Sigma_B = d ln chi^2,
      Sigma_Z = d ln Upsilon,
      Sigma_N = d ln(Lambda^2/K).
3. Conservative-shape-preserving reduction
      Xi_load = <2 d ln Lambda - dK>_N.
4. Naive common-shape branch gives Xi_load = -dK.
"""

from __future__ import annotations
import sympy as sp


def banner(title: str) -> None:
    line = "=" * 88
    print("\n" + line)
    print(title)
    print(line)


def expect_zero(name: str, expr: sp.Expr) -> None:
    expr = sp.simplify(sp.expand(expr))
    print(f"{name} = {expr}")
    if expr != 0:
        raise AssertionError(f"{name} is not zero")


banner("STAGE 158 — WALL-NORMALIZED LOAD/SHAPE FACTORIZATION")

# ---------------------------------------------------------------------------
# Exact algebraic factorization
# ---------------------------------------------------------------------------

K = sp.symbols("K", positive=True, real=True)
chi = sp.symbols("chi", positive=True, real=True)
varpi = sp.symbols("varpi", positive=True, real=True)
c = sp.sqrt(K) * varpi * chi
B0 = sp.simplify(c**2 / varpi**2)
expect_zero("B0 - K*chi^2", B0 - K * chi**2)

Ou, Ow, R, gu, gw = sp.symbols("Ou Ow R gu gw", positive=True, real=True)
Ou_hat, Ow_hat, R_hat, gu_hat, gw_hat = sp.symbols("Ou_hat Ow_hat R_hat gu_hat gw_hat", positive=True, real=True)

Delta = Ou * Ow - R**2
Q = gu**2 * Ow + 2 * gu * gw * R + gw**2 * Ou
P = Ou * gw + R * gu

subs_hat = {
    Ou: K * Ou_hat,
    Ow: K * Ow_hat,
    R:  K * R_hat,
    gu: K * gu_hat,
    gw: K * gw_hat,
}

expect_zero(
    "Delta - K^2*Delta_hat",
    sp.expand(Delta.subs(subs_hat)) - K**2 * (Ou_hat * Ow_hat - R_hat**2),
)
expect_zero(
    "Q - K^3*Q_hat",
    sp.expand(Q.subs(subs_hat)) - K**3 * (gu_hat**2 * Ow_hat + 2 * gu_hat * gw_hat * R_hat + gw_hat**2 * Ou_hat),
)
expect_zero(
    "P - K^2*P_hat",
    sp.expand(P.subs(subs_hat)) - K**2 * (Ou_hat * gw_hat + R_hat * gu_hat),
)

Upsilon = sp.simplify((Q / (K * Delta)).subs(subs_hat))
Lambda = sp.simplify((P / Delta).subs(subs_hat))
expect_zero("Z0 - K*Upsilon", sp.simplify((Q / Delta).subs(subs_hat)) - K * Upsilon)
expect_zero("N0 - Lambda^2", sp.simplify((P**2 / Delta**2).subs(subs_hat)) - Lambda**2)

# ---------------------------------------------------------------------------
# Differential/logarithmic identities
# ---------------------------------------------------------------------------

banner("Differential defect identities")

eps = sp.symbols("eps", real=True)
kappa = sp.symbols("kappa", real=True)
schi = sp.symbols("schi", real=True)
su, sw, sr, sgu, sgw = sp.symbols("su sw sr sgu sgw", real=True)

K0, chi0, varpi0 = sp.symbols("K0 chi0 varpi0", positive=True, real=True)
Ou0, Ow0, R0, gu0, gw0 = sp.symbols("Ou0 Ow0 R0 gu0 gw0", positive=True, real=True)


def dlog(expr: sp.Expr) -> sp.Expr:
    return sp.simplify(sp.diff(sp.log(sp.simplify(expr)), eps).subs(eps, 0))

subs_eps = {
    K: K0 * sp.exp(kappa * eps),
    chi: chi0 * sp.exp(schi * eps),
    varpi: varpi0,
    Ou_hat: Ou0 * sp.exp(su * eps),
    Ow_hat: Ow0 * sp.exp(sw * eps),
    R_hat:  R0 * sp.exp(sr * eps),
    gu_hat: gu0 * sp.exp(sgu * eps),
    gw_hat: gw0 * sp.exp(sgw * eps),
}

# Sigma_B = 2 dln(c/varpi) - dK = dln(chi^2)
expr_c = c.subs({K: K0 * sp.exp(kappa * eps), chi: chi0 * sp.exp(schi * eps), varpi: varpi0})
Sigma_B_direct = sp.simplify(2 * dlog(expr_c) - kappa)
Sigma_B_shape = sp.simplify(dlog((chi**2).subs(subs_eps)))
expect_zero("Sigma_B - dln(chi^2)", Sigma_B_direct - Sigma_B_shape)

# Sigma_Z = dln(Q/Delta) - dK = dln(Upsilon)
expr_Z = (Q / Delta).subs(subs_hat).subs(subs_eps)
expr_U = Upsilon.subs(subs_eps)
Sigma_Z_direct = sp.simplify(dlog(expr_Z) - kappa)
Sigma_Z_shape = sp.simplify(dlog(expr_U))
expect_zero("Sigma_Z - dln(Upsilon)", Sigma_Z_direct - Sigma_Z_shape)

# Sigma_N = 2 dln(P/Delta) - dK = dln(Lambda^2/K)
expr_ratio = (P / Delta).subs(subs_hat).subs(subs_eps)
expr_L = Lambda.subs(subs_eps)
Sigma_N_direct = sp.simplify(2 * dlog(expr_ratio) - kappa)
Sigma_N_shape = sp.simplify(dlog((Lambda**2 / K).subs(subs_eps)))
expect_zero("Sigma_N - dln(Lambda^2/K)", Sigma_N_direct - Sigma_N_shape)
expect_zero("Sigma_N - (2 dln Lambda - dK)", Sigma_N_direct - (2 * dlog(expr_L) - kappa))

# ---------------------------------------------------------------------------
# Conservative-shape-preserving and common-shape branches
# ---------------------------------------------------------------------------

banner("Conservative-shape-preserving reductions")

# Freeze all normalized shapes.
Sigma_B_cons = sp.simplify(Sigma_B_direct.subs({schi: 0}))
Sigma_Z_cons = sp.simplify(Sigma_Z_direct.subs({su: 0, sw: 0, sr: 0, sgu: 0, sgw: 0}))
Sigma_N_common = sp.simplify(Sigma_N_direct.subs({su: 0, sw: 0, sr: 0, sgu: 0, sgw: 0}))
expect_zero("Conservative-shape branch Sigma_B", Sigma_B_cons)
expect_zero("Conservative-shape branch Sigma_Z", Sigma_Z_cons)
expect_zero("Common-shape branch Sigma_N + dK", Sigma_N_common + kappa)

print("\nConclusions:")
print("  B0 = K chi^2,  Z0 = K Upsilon,  N0 = Lambda^2.")
print("  Therefore")
print("    Sigma_B = d ln chi^2")
print("    Sigma_Z = d ln Upsilon")
print("    Sigma_N = d ln(Lambda^2/K) = 2 d ln Lambda - dK")
print("  If all wall-normalized shapes are frozen, then Sigma_N = -dK.")
print("  So naive common self-similarity does not kill the grouped defect.")
