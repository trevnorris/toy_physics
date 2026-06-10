#!/usr/bin/env python3
"""
Stage 045 SymPy audit.

Checks:
1. The coherent local D/N support kernel implies g_B g_R = g_W g_S exactly.
2. The mixed and support interference ratios coincide: rho_0 = sigma_0.
3. The common direction factor R_tr has the exact range identities used in the note.
4. The Stage-044 quadratic branch equation collapses to the one-parameter tracking law.
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


banner("STAGE 045 — COHERENT LOCAL D/N KERNEL TRACKING AUDIT")

# ---------------------------------------------------------------------------
# 1. Coherent local support density => exact tracking condition
# ---------------------------------------------------------------------------

banner("1. Exact tracking identity from the coherent local kernel")

lam_W, lam_phi, gamma = sp.symbols("lambda_W lambda_phi gamma", positive=True, real=True)
mu_eta, mu_U, mu_W, mu_phi = sp.symbols("mu_eta mu_U mu_W mu_phi", positive=True, real=True)
K_U = sp.symbols("K_U", positive=True, real=True)
g_U = sp.symbols("g_U", real=True)

# Mass-normalized amplitudes produced by
#   [lambda_W W + lambda_phi phi] [eta - gamma U].
W_sym, phi_sym, eta_sym, U_sym = sp.symbols("W_sym phi_sym eta_sym U_sym")
coupling_density = sp.expand((lam_W*W_sym + lam_phi*phi_sym) * (eta_sym - gamma*U_sym))
# Extract bilinear coefficients from the kernel directly.
c_W_eta   = coupling_density.coeff(W_sym).coeff(eta_sym)
c_W_U     = coupling_density.coeff(W_sym).coeff(U_sym)
c_phi_eta = coupling_density.coeff(phi_sym).coeff(eta_sym)
c_phi_U   = coupling_density.coeff(phi_sym).coeff(U_sym)
g_W_ext = c_W_eta   / sp.sqrt(mu_eta * mu_W)
g_R_ext = -c_W_U    / sp.sqrt(mu_U   * mu_W)
g_B_ext = c_phi_eta / sp.sqrt(mu_eta * mu_phi)
g_S_ext = -c_phi_U  / sp.sqrt(mu_U   * mu_phi)
# Reference (the form the script historically used).
g_W = lam_W / sp.sqrt(mu_eta * mu_W)
g_R = gamma * lam_W / sp.sqrt(mu_U * mu_W)
g_B = lam_phi / sp.sqrt(mu_eta * mu_phi)
g_S = gamma * lam_phi / sp.sqrt(mu_U * mu_phi)
# Cross-check: extracted vs reference (catches sign / coefficient errors in the kernel).
expect_zero("g_W extracted - reference", g_W_ext - g_W)
expect_zero("g_R extracted - reference", g_R_ext - g_R)
expect_zero("g_B extracted - reference", g_B_ext - g_B)
expect_zero("g_S extracted - reference", g_S_ext - g_S)
# Now the kernel-derived identity becomes a non-trivial check.
expect_zero("g_B g_R - g_W g_S", g_B_ext * g_R_ext - g_W_ext * g_S_ext)

rho_0 = sp.simplify(g_R_ext * g_U / (K_U * g_W_ext))
sigma_0 = sp.simplify(g_U * g_S_ext / (K_U * g_B_ext))
print("rho_0   =", rho_0)
print("sigma_0 =", sigma_0)
expect_zero("rho_0 - sigma_0", rho_0 - sigma_0)

chi_0 = sp.symbols("chi_0", positive=True, real=True)
chi_subs = {rho_0: chi_0, sigma_0: chi_0}

# ---------------------------------------------------------------------------
# 2. Common direction factor and exact range identities
# ---------------------------------------------------------------------------

banner("2. Common direction factor and exact range identities")

delta_U = sp.symbols("delta_U", positive=True, real=True)
R_tr = sp.simplify((1 + chi_0 / (1 + delta_U)) / (1 + chi_0))
print("R_tr =", R_tr)

expr1 = sp.simplify(1 - R_tr)
expr2 = sp.simplify(R_tr - 1 / (1 + delta_U))
print("1 - R_tr =", expr1)
print("R_tr - 1/(1+delta_U) =", expr2)

expect_zero(
    "range identity 1",
    expr1 - chi_0 * delta_U / ((1 + chi_0) * (1 + delta_U)),
)
expect_zero(
    "range identity 2",
    expr2 - delta_U / ((1 + chi_0) * (1 + delta_U)),
)

# ---------------------------------------------------------------------------
# 3. Total baseline on the coherent branch
# ---------------------------------------------------------------------------

banner("3. Exact total baseline")

Z_W, Z_phi = sp.symbols("Z_W Z_phi", positive=True, real=True)
eps_eta, eps_W_split, eps_phi_split = sp.symbols(
    "eps_eta eps_W_split eps_phi_split", real=True
)

prefactor = 8 * (1 + chi_0) ** 2 / (sp.pi ** 2 * (1 - eps_eta))
M_mix = sp.simplify(prefactor * Z_W / (1 - eps_W_split))
M_supp = sp.simplify(prefactor * Z_phi / (1 - eps_phi_split))
M_tr = sp.simplify(M_mix + M_supp)
print("M_mix  =", M_mix)
print("M_supp =", M_supp)
print("M_tr   =", M_tr)
# M_mix and M_supp are carried forward from Stages 022 and 026 in symbolic form;
# the substantive verification of the prefactor structure lives in those stages.

# ---------------------------------------------------------------------------
# 4. Stage-044 quadratic branch equation collapses to tracking law
# ---------------------------------------------------------------------------

banner("4. Exact collapse of the Stage-044 branch equation")

xi, delta, lam0 = sp.symbols("xi delta lambda_0", positive=True, real=True)
R_U, R_phi = sp.symbols("R_U R_phi", real=True)
Mmix, Msupp = sp.symbols("Mmix Msupp", real=True)

# Stage-044 continuum-selected branch equation in the reduced notation.
branch_eq = sp.simplify(
    Msupp - (
        xi * (delta + xi) - Mmix * (delta + (1 + lam0 * R_U ** 2) * xi)
    ) / (
        delta + (1 + lam0 * R_phi ** 2) * xi - Mmix * lam0 * (R_U - R_phi) ** 2
    )
)

branch_track = sp.together(sp.simplify(branch_eq.subs(R_phi, R_U)))
num_track = sp.expand(sp.factor(sp.together(branch_track).as_numer_denom()[0]))
den_track = sp.factor(sp.together(branch_track).as_numer_denom()[1])
print("tracking numerator =", num_track)
print("tracking denominator =", den_track)

M_tr_sym = sp.symbols("M_tr", real=True)
collapsed_num = sp.expand(
    xi ** 2 + (delta - M_tr_sym * (1 + lam0 * R_U ** 2)) * xi - delta * M_tr_sym
)
expect_zero(
    "tracking quadratic collapse",
    num_track + collapsed_num.subs(M_tr_sym, Mmix + Msupp),
)

M_tr_req = sp.simplify(sp.solve(sp.Eq(collapsed_num, 0), M_tr_sym)[0])
print("M_tr required on tracking branch =", M_tr_req)
expect_zero(
    "G_tr formula",
    M_tr_req - xi * (delta + xi) / (delta + (1 + lam0 * R_U ** 2) * xi),
)

R_target = sp.symbols("R_target", positive=True, real=True)
lam0_dn = sp.Rational(2, 9)
G_tr_dn = sp.simplify(M_tr_req.subs(lam0, lam0_dn))
G_tr_expected = sp.simplify(9 * xi * (delta + xi) / (9 * delta + (9 + 2 * R_U ** 2) * xi))
expect_zero("G_tr D/N specialization", G_tr_dn - G_tr_expected)

# Import Stage-044 F_cont residual; substitute tracking (R_phi -> R_U)
# plus the D/N value (lambda_0 -> 2/9).
# See: scripts/moving_throat_pde_stage044_continuum_selected_rank2_sympy_audit.py:82-90,140-146
D_cont_stage044 = sp.simplify(
    (delta + xi - Mmix * lam0 * R_U * (R_U - R_phi)) ** 2
    + lam0 * (Mmix * (R_U - R_phi) + R_phi * xi) ** 2
)
F_cont_stage044 = sp.simplify(
    (delta + (1 + lam0 * R_U * R_phi) * xi) ** 2
    * (delta + (1 + lam0 * R_phi) * xi - Mmix * lam0 * (R_U - R_phi) * (R_U - 1)) ** 2
    / ((1 - xi) * D_cont_stage044 ** 2)
)
F_track_stage044 = sp.simplify(F_cont_stage044.subs(R_phi, R_U))
F_track_expected = sp.simplify(
    (delta + (1 + lam0 * R_U ** 2) * xi) ** 2
    * (delta + (1 + lam0 * R_U) * xi) ** 2
    / ((1 - xi) * ((delta + xi) ** 2 + lam0 * R_U ** 2 * xi ** 2) ** 2)
)
expect_zero("Stage-044 tracking F collapse", F_track_stage044 - F_track_expected)

F_tr_from_stage044 = sp.simplify(
    F_cont_stage044.subs([(R_phi, R_U), (lam0, lam0_dn)])
)
F_tr_expected = sp.simplify(
    (9 * delta + (9 + 2 * R_U ** 2) * xi) ** 2
    * (9 * delta + (9 + 2 * R_U) * xi) ** 2
    / (81 * (1 - xi) * (9 * delta ** 2 + 18 * delta * xi + (9 + 2 * R_U ** 2) * xi ** 2) ** 2)
)
expect_zero("F_tr collapse from Stage-044 residual", F_tr_from_stage044 - F_tr_expected)
print("coherent normalization residual =", sp.simplify(R_target - F_tr_from_stage044))

print("\nAll Stage-045 symbolic checks passed.")
