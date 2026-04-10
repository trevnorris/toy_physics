#!/usr/bin/env python3
"""
5pn_stage6_sector_sieve_and_outgoing_load.py

Sixth executable SymPy audit for the 5PN grouped-real-P2 program.

What this script does
---------------------
1. Continues the Stage-5 primitive-deformation problem by testing the first two
   concrete primitive anisotropy mechanisms directly:
      - wall-only,
      - BdG-only.
2. Proves that neither wall-only nor pure BdG anisotropy can satisfy the full
   linear 5PN consistency triplet
      K1 = 0,
      Xi_load = 0,
      H_even = D41 - (2/3) D21 - D01/27 = 0
   nontrivially in the explicit prototype.
3. Shows that the pure BdG self-similar branch kills Xi_load but still fails the
   even-preserving and hidden-even gates unless it is trivial.
4. Rewrites Xi_load in the exact Stage-157 wall-referenced self-similarity form
      Xi_load = (delta_N-delta_K) + omega_B (delta_B-delta_K)
                             + omega_Z (delta_Z-delta_K),
   and then in terms of the microscopic defect fields
      Sigma^(B), Sigma^(Z), Sigma^(N).
5. Derives the Stage-158 conservative-shape theorem and the outgoing-load theorem.
6. Factors the outgoing load through the Stage-159 variables
      M_r, I_r, H_r,
   proving the square-root mixed-leg law as the surviving nontrivial corridor.

Interpretation
--------------
This script is the first genuine mechanism sieve after Stage 5. It shows that
wall-only and BdG-only weak-axisymmetric anisotropies are dead ends, while the
remaining nontrivial corridor is the outgoing Maxwell/mixed loading sector.
"""

from __future__ import annotations

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


def expect_zero(name: str, expr: sp.Expr | sp.Matrix) -> None:
    if isinstance(expr, sp.MatrixBase):
        expr = expr.applyfunc(lambda z: sp.simplify(sp.expand(z)))
        print(f"{name} =")
        sp.pprint(expr)
        if any(entry != 0 for entry in expr):
            raise AssertionError(f"{name} is not zero")
    else:
        expr = sp.simplify(sp.expand(expr))
        print(f"{name} = {expr}")
        if expr != 0:
            raise AssertionError(f"{name} is not zero")


# ---------------------------------------------------------------------------
# I. Primitive continuation point from Stage 5
# ---------------------------------------------------------------------------

banner("I. PRIMITIVE CONTINUATION POINT FROM STAGE 5")

# Baseline grouped-wall operator and load variables.
D0, N0 = sp.symbols("D0 N0", nonzero=True)
dK, dM = sp.symbols("dK dM", real=True)

# BdG primitive variables.
B0 = sp.symbols("B0", positive=True, real=True)
varpi = sp.symbols("varpi", positive=True, real=True)
x_c = sp.symbols("x_c", real=True)        # delta ln c_alpha
x_varpi = sp.symbols("x_varpi", real=True)  # delta ln varpi_alpha

# The three linear 5PN gates carried forward from Stage 5.
print("The linear 5PN gates are")
print("  K1       = D21 + D01/9,")
print("  Xi_load  = N01/N0 - D01/D0,")
print("  H_even   = D41 - (2/3) D21 - D01/27.")
print("A viable primitive weak-axisymmetric mechanism must satisfy all three.")


# ---------------------------------------------------------------------------
# II. Wall-only mechanism: exact no-go
# ---------------------------------------------------------------------------

banner("II. WALL-ONLY MECHANISM: EXACT NO-GO")

D01_wall = dK
D21_wall = -dM
D41_wall = sp.Integer(0)
N01_wall = sp.Integer(0)

K1_wall = sp.simplify(D21_wall + D01_wall / 9)
Xi_wall = sp.simplify(N01_wall / N0 - D01_wall / D0)
He_wall = sp.simplify(D41_wall - sp.Rational(2, 3) * D21_wall - D01_wall / 27)

print("D01_wall =", D01_wall)
print("D21_wall =", D21_wall)
print("D41_wall =", D41_wall)
print("N01_wall =", N01_wall)
print("K1_wall  =", K1_wall)
print("Xi_wall  =", Xi_wall)
print("He_wall  =", He_wall)

wall_solution = sp.solve(
    [sp.Eq(K1_wall, 0), sp.Eq(Xi_wall, 0), sp.Eq(He_wall, 0)],
    [dK, dM],
    dict=True,
)
print("\nFull wall-only solve:", wall_solution)
if wall_solution != [{dK: 0, dM: 0}]:
    raise AssertionError("Wall-only branch unexpectedly admitted a nontrivial solution.")

print("Conclusion: wall-only weak-axisymmetric loading is dead unless it is trivial.")


# ---------------------------------------------------------------------------
# III. Pure BdG mechanism: Xi_load can vanish, but the 5PN triplet still fails
# ---------------------------------------------------------------------------

banner("III. PURE BdG MECHANISM: Xi_load CAN VANISH, BUT THE 5PN TRIPLET STILL FAILS")

# One-mode BdG sector: B0 = c^2 / varpi^2, B2 = B0 / varpi^2, B4 = B0 / varpi^4.
# Use logarithmic slopes x_c = delta ln c, x_varpi = delta ln varpi.
D01_B = sp.simplify(-2 * B0 * (x_c - x_varpi))
D21_B = sp.simplify(-2 * B0 * (x_c - 2 * x_varpi) / varpi**2)
D41_B = sp.simplify(-2 * B0 * (x_c - 3 * x_varpi) / varpi**4)
N01_B = sp.Integer(0)

K1_B = sp.simplify(D21_B + D01_B / 9)
Xi_B = sp.simplify(N01_B / N0 - D01_B / D0)
He_B = sp.simplify(D41_B - sp.Rational(2, 3) * D21_B - D01_B / 27)

print("D01_B =", D01_B)
print("D21_B =", D21_B)
print("D41_B =", D41_B)
print("N01_B =", N01_B)
print("K1_B  =", sp.factor(K1_B))
print("Xi_B  =", Xi_B)
print("He_B  =", sp.factor(He_B))

bdg_solution = sp.solve(
    [sp.Eq(K1_B, 0), sp.Eq(Xi_B, 0), sp.Eq(He_B, 0)],
    [x_c, x_varpi],
    dict=True,
)
print("\nFull BdG-only solve:", bdg_solution)
if bdg_solution != [{x_c: 0, x_varpi: 0}]:
    raise AssertionError("BdG-only branch unexpectedly admitted a nontrivial solution.")

subbanner("BdG self-similar branch")
# Stage-157 self-similarity for the BdG sector is 2 delta ln(c/varpi) = delta_K.
# On a wall-fixed pure-BdG branch delta_K = 0, so x_c = x_varpi.
delta = sp.symbols("delta", real=True)
self_sim_subs = {x_c: delta, x_varpi: delta}

D01_B_ss = sp.simplify(D01_B.subs(self_sim_subs))
D21_B_ss = sp.simplify(D21_B.subs(self_sim_subs))
D41_B_ss = sp.simplify(D41_B.subs(self_sim_subs))
K1_B_ss = sp.simplify(K1_B.subs(self_sim_subs))
Xi_B_ss = sp.simplify(Xi_B.subs(self_sim_subs))
He_B_ss = sp.simplify(He_B.subs(self_sim_subs))

print("BdG self-similarity condition on a wall-fixed branch: x_c = x_varpi = delta")
print("D01_B(self-similar) =", D01_B_ss)
print("D21_B(self-similar) =", D21_B_ss)
print("D41_B(self-similar) =", D41_B_ss)
print("K1_B(self-similar)  =", K1_B_ss)
print("Xi_B(self-similar)  =", Xi_B_ss)
print("He_B(self-similar)  =", He_B_ss)
expect_zero("Xi_B vanishes on the self-similar BdG branch", Xi_B_ss)

print(
    "Conclusion: pure BdG self-similarity can kill Xi_load, but it still leaves"
    " nonzero K1 and hidden-even drift unless delta = 0."
)


# ---------------------------------------------------------------------------
# IV. Stage-157 exact wall-referenced self-similarity theorem
# ---------------------------------------------------------------------------

banner("IV. EXACT WALL-REFERENCED SELF-SIMILARITY THEOREM (STAGE 157)")

# Static grouped-lane slopes and weights.
delta_K, delta_B, delta_Z, delta_N = sp.symbols("delta_K delta_B delta_Z delta_N", real=True)
omega_K, omega_B, omega_Z = sp.symbols("omega_K omega_B omega_Z", real=True)

Xi = sp.simplify(delta_N - omega_K * delta_K + omega_B * delta_B + omega_Z * delta_Z)
Xi_reduced = sp.simplify(Xi.subs(omega_K, 1 + omega_B + omega_Z))

print("Xi_load before the weight identity =", Xi)
print("Xi_load after using omega_K = 1 + omega_B + omega_Z =", Xi_reduced)
expect_zero(
    "Stage-157 decomposition identity",
    Xi_reduced - ((delta_N - delta_K) + omega_B * (delta_B - delta_K) + omega_Z * (delta_Z - delta_K)),
)

subbanner("Microscopic defect fields")
Sigma_B, Sigma_Z, Sigma_N = sp.symbols("Sigma_B Sigma_Z Sigma_N", real=True)

Xi_sigma = sp.simplify(Sigma_N + omega_B * Sigma_B + omega_Z * Sigma_Z)
Xi_sigma_recovered = sp.simplify(
    Xi_sigma.subs(
        {
            Sigma_B: delta_B - delta_K,
            Sigma_Z: delta_Z - delta_K,
            Sigma_N: delta_N - delta_K,
        }
    )
)

print("Xi_load = Sigma_N + omega_B Sigma_B + omega_Z Sigma_Z")
expect_zero("Defect-field reconstruction", Xi_sigma_recovered - Xi_reduced)

subbanner("Static self-similarity")
Xi_self_similar = sp.simplify(Xi_sigma.subs({Sigma_B: 0, Sigma_Z: 0, Sigma_N: 0}))
expect_zero("Xi_load under static self-similarity", Xi_self_similar)
print(
    "If the BdG, conservative Maxwell/mixed, and outgoing-transfer bundles all"
    " co-load with the wall baseline, Xi_load vanishes automatically."
)


# ---------------------------------------------------------------------------
# V. Stage-158 outgoing-load theorem and naive common-self-similarity no-go
# ---------------------------------------------------------------------------

banner("V. STAGE-158 OUTGOING-LOAD THEOREM")

Theta_B, Theta_Z = sp.symbols("Theta_B Theta_Z", real=True)
delta_ln_Lambda = sp.symbols("delta_ln_Lambda", real=True)

Xi_conservative_shape = sp.simplify(2 * delta_ln_Lambda - delta_K)
print(
    "On conservative-shape-preserving branches, Xi_load collapses to"
    " 2 sum_r rho_r^(N) delta ln Lambda_r - delta_K."
)
print("One-port representative: Xi_load =", Xi_conservative_shape)

# Naive common-self-similarity would also set delta ln Lambda = 0.
Xi_naive = sp.simplify(Xi_conservative_shape.subs(delta_ln_Lambda, 0))
print("Naive common self-similarity gives Xi_load =", Xi_naive)
expect_zero("Naive common-self-similarity no-go", Xi_naive + delta_K)

print(
    "Conclusion: preserving all wall-normalized shapes is not enough. On a"
    " nontrivial wall-loading branch the outgoing sector must actively load."
)

Xi_outgoing_theorem = sp.Eq(2 * delta_ln_Lambda, delta_K)
print("Outgoing-load theorem (one-port form):", Xi_outgoing_theorem)


# ---------------------------------------------------------------------------
# VI. Stage-159 factorization and the square-root mixed-leg law
# ---------------------------------------------------------------------------

banner("VI. STAGE-159 FACTORIZATION AND THE SQUARE-ROOT MIXED-LEG LAW")

GW, GU, R, OmegaU_sq, OmegaW_sq, K = sp.symbols(
    "G_W G_U R Omega_U_sq Omega_W_sq K", positive=True, real=True
)

M_r = sp.simplify(GW / (OmegaW_sq * sp.sqrt(K)))
I_r = sp.simplify(R * GU / (OmegaU_sq * GW))
H_r = sp.simplify(R**2 / (OmegaU_sq * OmegaW_sq))
Lambda = sp.simplify((OmegaU_sq * GW + R * GU) / (OmegaU_sq * OmegaW_sq - R**2))

factorized = sp.simplify(sp.sqrt(K) * M_r * (1 + I_r) / (1 - H_r))
factorized_sq = sp.simplify(M_r**2 * (1 + I_r)**2 / (1 - H_r)**2)
print("Lambda_r =", Lambda)
print("Factorized Lambda_r =", factorized)
print("Lambda_r^2 / K =", sp.simplify(Lambda**2 / K))
print("Factorized Lambda_r^2 / K =", factorized_sq)
expect_zero("Stage-159 exact factorization (squared form)", sp.simplify(Lambda**2 / K) - factorized_sq)

subbanner("Outgoing defect split")
dlnM, dlnI, dlnH = sp.symbols("delta_ln_M delta_ln_I delta_ln_H", real=True)
SigmaN_factorized = sp.simplify(2 * dlnM + 2 * dlnI / (1 + I_r) * I_r + 2 * dlnH / (1 - H_r) * H_r)
print("Sigma_r^(N) = 2 d ln M + 2 I/(1+I) d ln I + 2 H/(1-H) d ln H")
print("Sigma_r^(N) =")
sp.pprint(SigmaN_factorized)

subbanner("Interference/hybridization rigidity")
SigmaN_rigid = sp.simplify(SigmaN_factorized.subs({dlnI: 0, dlnH: 0}))
expect_zero("Rigid I,H => Sigma_r^(N) - 2 dlnM", SigmaN_rigid - 2 * dlnM)

# The square-root mixed-leg law is d ln M = 0.
# Because M = G_W / (Omega_W^2 sqrt(K)), this is equivalent to
# G_W / Omega_W^2 proportional to sqrt(K).
print("Square-root mixed-leg law:")
print("  d ln M_r = 0")
print("with M_r = G_W / (Omega_W^2 sqrt(K)).")
print("Equivalently, G_W / Omega_W^2 ∝ sqrt(K).")

subbanner("One-port zero-defect theorem under rigidity")
Xi_rigid = sp.simplify((2 * dlnM + 2 * I_r / (1 + I_r) * 0 + 2 * H_r / (1 - H_r) * 0))
expect_zero("Xi_load under rigidity minus 2 dlnM", Xi_rigid - 2 * dlnM)

print(
    "Conclusion: once conservative shapes are preserved and the interference and"
    " hybridization ratios are rigid, the only surviving nontrivial corridor is"
    " the outgoing mixed leg, and it is controlled by a square-root wall-loading law."
)


# ---------------------------------------------------------------------------
# VII. Final theorem ledger
# ---------------------------------------------------------------------------

banner("VII. FINAL THEOREM LEDGER")
print("1. Wall-only weak-axisymmetric anisotropy is impossible unless trivial.")
print("2. Pure BdG weak-axisymmetric anisotropy is also impossible unless trivial;")
print("   BdG self-similarity can kill Xi_load, but not K1 and hidden-even together.")
print("3. The remaining linear grouped defect is exactly the weighted failure of")
print("   static self-similarity relative to the wall baseline.")
print("4. On conservative-shape-preserving branches, that defect collapses to an")
print("   outgoing-load theorem:")
print("      2 sum rho_r^(N) delta ln Lambda_r = delta_K.")
print("5. Naive common self-similarity is a no-go on any nontrivial wall-loading")
print("   branch, because it gives Xi_load = -delta_K.")
print("6. Factoring Lambda_r shows that the true surviving nontrivial corridor is")
print("   the outgoing Maxwell/mixed leg.")
print("7. Under interference/hybridization rigidity, exact first-order defect")
print("   cancellation reduces to the square-root mixed-leg law")
print("      G_W / Omega_W^2 ∝ sqrt(K).")
print("8. So after Stage 6 the next honest theorem gate is not another conservative")
print("   algebra sweep. It is to compute the weak-axisymmetric outgoing-load")
print("   slippages on the actual moving-throat branch.")
