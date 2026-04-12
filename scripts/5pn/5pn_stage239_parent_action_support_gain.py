#!/usr/bin/env python3
from __future__ import annotations

import sympy as sp

from fivepn_stage199_201_common import banner, subbanner, expect_zero

"""
Stage 239 — parent-action projection of the microscopic support/source gain.
"""


def main() -> None:
    banner("STAGE 239 — PARENT-ACTION SUPPORT/SOURCE GAIN")

    rho_star, g_phi, O_sp = sp.symbols("rho_star g_phi O_sigma_phi", positive=True, real=True)
    N_ss, N_pp = sp.symbols("N_sigma_sigma N_phi_phi", positive=True, real=True)
    m, c_s_star = sp.symbols("m c_s_star", positive=True, real=True)
    K_X, T_X, L = sp.symbols("K_X T_X L", positive=True, real=True)
    kappa = sp.symbols("kappa", positive=True, real=True)
    C_sp_sq = sp.symbols("C_sigma_phi_sq", nonnegative=True, real=True)
    chi_sigma_eff = sp.symbols("chi_sigma_eff", positive=True, real=True)
    G_micro = sp.symbols("G_micro", positive=True, real=True)
    Xi_micro = sp.symbols("Xi_micro", positive=True, real=True)

    hprime = sp.simplify(m * c_s_star**2 / rho_star)
    Theta_sigma = sp.simplify(hprime * N_ss)
    chi_sigma_expr = sp.simplify(1 / Theta_sigma)
    Lambda_phi = sp.simplify(g_phi * O_sp)
    G_micro_expr = sp.simplify(chi_sigma_expr * Lambda_phi**2 / K_X)
    Xi_micro_expr = sp.simplify(G_micro_expr * kappa)
    Xi_micro_TX_expr = sp.simplify(Xi_micro_expr.subs(kappa, K_X * L**2 / T_X))

    C_def = sp.simplify(O_sp**2 / (N_ss * N_pp))
    G_micro_coh = sp.simplify((rho_star * g_phi**2 * N_pp / (m * c_s_star**2 * K_X)) * C_sp_sq)
    G_max = sp.simplify(rho_star * g_phi**2 * N_pp / (m * c_s_star**2 * K_X))

    subbanner("I. Exact n=5 compressional stiffness")
    print("h'(rho_*) =")
    sp.pprint(hprime)
    print()
    print("Theta_sigma =")
    sp.pprint(Theta_sigma)
    print()
    print("chi_sigma^(eff) =")
    sp.pprint(chi_sigma_expr)
    expect_zero(
        "exact source-susceptibility identity",
        sp.simplify(chi_sigma_expr - rho_star / (m * c_s_star**2 * N_ss)),
    )

    subbanner("II. Exact parent microscopic gain")
    print("Lambda_phi =")
    sp.pprint(Lambda_phi)
    print()
    print("G_micro =")
    sp.pprint(G_micro_expr)
    expect_zero(
        "exact parent gain identity",
        sp.simplify(G_micro_expr - rho_star * g_phi**2 * O_sp**2 / (m * c_s_star**2 * K_X * N_ss)),
    )

    subbanner("III. Coherence-factor factorization")
    print("C_(sigma,phi)^2 =")
    sp.pprint(C_def)
    print()
    print("G_micro in coherence form =")
    sp.pprint(G_micro_coh)
    expect_zero(
        "coherence-factor identity",
        sp.simplify(G_micro_expr.subs(O_sp**2, C_sp_sq * N_ss * N_pp) - G_micro_coh),
    )
    print()
    print("G_max(g_phi) =")
    sp.pprint(G_max)
    print("(attained only at perfect alignment C_(sigma,phi)^2 = 1)")

    subbanner("IV. Exact fixed-point strength Xi_micro")
    print("Xi_micro = kappa G_micro =")
    sp.pprint(Xi_micro_expr)
    print()
    print("With kappa = K_X L^2 / T_X:")
    sp.pprint(Xi_micro_TX_expr)
    expect_zero(
        "Xi_micro projected identity",
        sp.simplify(Xi_micro_TX_expr - rho_star * g_phi**2 * O_sp**2 * L**2 / (m * c_s_star**2 * T_X * N_ss)),
    )

    banner("STAGE 239 LEDGER")
    print("1. The microscopic support/source gain is no longer a phenomenological placeholder.")
    print("   On the projected parent branch it is")
    print("      G_micro = rho_* g_phi^2 O_(sigma,phi)^2 / [m c_s*^2 K_X N_(sigma,sigma)].")
    print("2. Equivalently,")
    print("      G_micro = [rho_* g_phi^2 N_(phi,phi)/(m c_s*^2 K_X)] C_(sigma,phi)^2,")
    print("   with 0 <= C_(sigma,phi)^2 <= 1 by Cauchy-Schwarz.")
    print("3. The exact best-case gain at fixed g_phi is therefore")
    print("      G_max(g_phi) = rho_* g_phi^2 N_(phi,phi)/(m c_s*^2 K_X).")
    print("4. The operator-selected fixed-point strength is")
    print("      Xi_micro = kappa G_micro = rho_* g_phi^2 O_(sigma,phi)^2 L^2 / [m c_s*^2 T_X N_(sigma,sigma)].")
    print("5. So the remaining support/source theorem gap is now a parent-overlap problem, not an abstract gain problem.")


if __name__ == "__main__":
    main()
