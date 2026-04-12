#!/usr/bin/env python3
"""
Step 11 for the g-2 rebuild from the moving-throat PDE base.

What this script does
---------------------
1. Introduces the microscopic coherent-kernel variables behind Step 10:
       chi_0, delta_U, epsilon_W, epsilon_eta, Z_W/Omega_W^2.
2. Verifies their logarithmic grouped weak-axisymmetric drifts in terms of the
   underlying microscopic kernel drifts
       (gamma_1, c_1, lambda_1, mu_1, tau_1, kappa_U, kappa_eta, kappa_W).
3. Rewrites the Step-10 coherent-branch defect Xi_1 directly in those
   microscopic slippage coordinates
       (Sigma_Z, Sigma_chi, Sigma_epsilon, Sigma_delta, Sigma_eta).
4. Shows that the exact tracking factor depends only on the single combination
       Sigma_tr = (1+chi_0) Sigma_delta + (1+delta_U) Sigma_chi.
5. Defines the genuine nontracking transfer-shape slippage Sigma_nt and proves
   the exact triangular normal form
       Theta_1 = -C_tr Sigma_tr,
       Xi_1    = A_tr Sigma_tr + Sigma_nt,
       R_1+Xi_1 = -(epsilon_eta/(1-epsilon_eta)) Sigma_eta.
6. Inverts the triangular system exactly and records the quartic anomaly gates.

Interpretation
--------------
Step 10 showed that coherent support itself is support-blind at the level of the
first grouped defect. Step 11 sharpens that result: the entire coherent
weak-axisymmetric problem collapses to three branch-adapted scalars,

    Sigma_tr, Sigma_nt, Sigma_eta,

rather than a broad list of microscopic placement drifts. So the next live
question for g-2 is no longer “which raw microscopic variables move?” but
“which of these three exact branch-adapted slippages survives on the actual
moving-throat branch?”
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
        simplified = expr.applyfunc(lambda z: sp.simplify(sp.expand(z)))
        print(f"{name} =")
        sp.pprint(simplified)
        if any(entry != 0 for entry in simplified):
            raise AssertionError(f"{name} is not zero")
    else:
        simplified = sp.simplify(sp.expand(expr))
        print(f"{name} = {simplified}")
        if simplified != 0:
            raise AssertionError(f"{name} is not zero")



def first_log_drift(expr: sp.Expr, s: sp.Symbol, lam: sp.Symbol) -> sp.Expr:
    return sp.simplify(sp.diff(sp.log(expr), s).subs(s, 0) / lam)



def step11_microscopic_slippage_normal_form() -> None:
    banner("STEP 11 — MICROSCOPIC SLIPPAGE NORMAL FORM")

    # ------------------------------------------------------------------
    # I. Microscopic coherent-kernel variables and their log drifts
    # ------------------------------------------------------------------
    subbanner("XI.1 — Microscopic coherent-kernel log coordinates")

    s, lam = sp.symbols("s lambda_A", real=True)
    gamma, c_etaU = sp.symbols("gamma c_etaU", positive=True, real=True)
    lambda_W = sp.symbols("lambda_W", positive=True, real=True)
    K_U, K_eta, K_W = sp.symbols("K_U K_eta K_W", positive=True, real=True)
    mu_W, T_U, L = sp.symbols("mu_W T_U L", positive=True, real=True)
    sigma = sp.symbols("sigma", positive=True, real=True)

    gamma1, c1, lambda1 = sp.symbols("gamma_1 c_1 lambda_1", real=True)
    mu1, tau1 = sp.symbols("mu_1 tau_1", real=True)
    kappa_U, kappa_eta, kappa_W = sp.symbols(
        "kappa_U kappa_eta kappa_W", real=True
    )

    gamma_A = gamma * sp.exp(s * lam * gamma1)
    c_A = c_etaU * sp.exp(s * lam * c1)
    lambda_W_A = lambda_W * sp.exp(s * lam * lambda1)
    K_U_A = K_U * sp.exp(s * lam * kappa_U)
    K_eta_A = K_eta * sp.exp(s * lam * kappa_eta)
    K_W_A = K_W * sp.exp(s * lam * kappa_W)
    mu_W_A = mu_W * sp.exp(s * lam * mu1)
    T_U_A = T_U * sp.exp(s * lam * tau1)

    chi_A = sp.simplify(gamma_A * c_A / K_U_A)
    delta_U_A = sp.simplify(sp.pi**2 * T_U_A / (L**2 * K_U_A))
    epsilon_eta_A = sp.simplify(c_A**2 / (K_U_A * K_eta_A))
    epsilon_W_A = sp.simplify(gamma_A**2 * lambda_W_A**2 * sigma / (K_U_A * K_W_A))
    Zratio_A = sp.simplify(lambda_W_A**2 * mu_W_A / (K_eta_A * K_W_A**2))

    Sigma_chi = first_log_drift(chi_A, s, lam)
    Sigma_delta = first_log_drift(delta_U_A, s, lam)
    Sigma_eta = first_log_drift(epsilon_eta_A, s, lam)
    Sigma_epsilon = first_log_drift(epsilon_W_A, s, lam)
    Sigma_Z = first_log_drift(Zratio_A, s, lam)

    expect_zero("Sigma_chi - (gamma_1 + c_1 - kappa_U)", Sigma_chi - (gamma1 + c1 - kappa_U))
    expect_zero("Sigma_delta - (tau_1 - kappa_U)", Sigma_delta - (tau1 - kappa_U))
    expect_zero("Sigma_eta - (2 c_1 - kappa_U - kappa_eta)", Sigma_eta - (2 * c1 - kappa_U - kappa_eta))
    expect_zero(
        "Sigma_epsilon - (2 gamma_1 + 2 lambda_1 - kappa_U - kappa_W)",
        Sigma_epsilon - (2 * gamma1 + 2 * lambda1 - kappa_U - kappa_W),
    )
    expect_zero(
        "Sigma_Z - (2 lambda_1 + mu_1 - kappa_eta - 2 kappa_W)",
        Sigma_Z - (2 * lambda1 + mu1 - kappa_eta - 2 * kappa_W),
    )

    print("chi_0 =")
    sp.pprint(sp.simplify(gamma * c_etaU / K_U))
    print("delta_U =")
    sp.pprint(sp.simplify(sp.pi**2 * T_U / (L**2 * K_U)))
    print("epsilon_eta =")
    sp.pprint(sp.simplify(c_etaU**2 / (K_U * K_eta)))
    print("epsilon_W =")
    sp.pprint(sp.simplify(gamma**2 * lambda_W**2 * sigma / (K_U * K_W)))
    print("Z_W/Omega_W^2 =")
    sp.pprint(sp.simplify(lambda_W**2 * mu_W / (K_eta * K_W**2)))

    # ------------------------------------------------------------------
    # II. Exact microscopic form of Xi_1
    # ------------------------------------------------------------------
    subbanner("XI.2 — Exact microscopic form of the coherent grouped defect")

    Zratio, chi0 = sp.symbols("Zratio chi_0", positive=True, real=True)
    delta_U = sp.symbols("delta_U", positive=True, real=True)
    epsilon_W, epsilon_eta = sp.symbols("epsilon_W epsilon_eta", positive=True, real=True)

    Sigma_Z_sym, Sigma_chi_sym, Sigma_epsilon_sym, Sigma_delta_sym, Sigma_eta_sym = sp.symbols(
        "Sigma_Z Sigma_chi Sigma_epsilon Sigma_delta Sigma_eta", real=True
    )

    Zratio_t = Zratio * sp.exp(s * lam * Sigma_Z_sym)
    chi_t = chi0 * sp.exp(s * lam * Sigma_chi_sym)
    delta_t = delta_U * sp.exp(s * lam * Sigma_delta_sym)
    epsilon_W_t = epsilon_W * sp.exp(s * lam * Sigma_epsilon_sym)
    epsilon_eta_t = epsilon_eta * sp.exp(s * lam * Sigma_eta_sym)

    epsilon = sp.simplify(epsilon_W * (1 - sp.Rational(2, 11) * delta_U / (1 + delta_U)))
    epsilon_t = sp.simplify(
        epsilon_W_t * (1 - sp.Rational(2, 11) * delta_t / (1 + delta_t))
    )
    T2_t = sp.simplify(Zratio_t * (1 + chi_t) ** 2 / (1 - epsilon_t) ** 2)

    Xi1 = first_log_drift(T2_t, s, lam)
    Xi1_formula = sp.simplify(
        Sigma_Z_sym
        + 2 * chi0 / (1 + chi0) * Sigma_chi_sym
        + 2 * epsilon_W / (1 - epsilon)
        * (
            (11 + 9 * delta_U) / (11 * (1 + delta_U)) * Sigma_epsilon_sym
            - 2 * delta_U / (11 * (1 + delta_U) ** 2) * Sigma_delta_sym
        )
    )

    expect_zero(
        "Xi_1 - exact microscopic slippage law",
        Xi1 - Xi1_formula,
    )

    print("Xi_1 =")
    sp.pprint(Xi1_formula)

    # ------------------------------------------------------------------
    # III. Tracking-factor drift and the tracking combination Sigma_tr
    # ------------------------------------------------------------------
    subbanner("XI.3 — Exact tracking slippage combination")

    Rtr_t = sp.simplify((1 + chi_t / (1 + delta_t)) / (1 + chi_t))
    Theta1 = first_log_drift(Rtr_t, s, lam)

    Sigma_tr = sp.simplify((1 + chi0) * Sigma_delta_sym + (1 + delta_U) * Sigma_chi_sym)
    C_tr = sp.simplify(
        chi0 * delta_U / ((1 + chi0) * (1 + delta_U) * (1 + chi0 + delta_U))
    )

    expect_zero("Theta_1 + C_tr Sigma_tr", Theta1 + C_tr * Sigma_tr)

    print("Sigma_tr =")
    sp.pprint(Sigma_tr)
    print("Theta_1 =")
    sp.pprint(sp.simplify(-C_tr * Sigma_tr))

    # ------------------------------------------------------------------
    # IV. Genuine nontracking transfer-shape slippage and triangular form
    # ------------------------------------------------------------------
    subbanner("XI.4 — Exact nontracking slippage and triangular observable law")

    A_tr = sp.simplify(2 * chi0 / ((1 + chi0) * (1 + delta_U)))
    Sigma_nt = sp.simplify(
        Sigma_Z_sym
        + 2 * epsilon_W / (1 - epsilon) * (11 + 9 * delta_U) / (11 * (1 + delta_U)) * Sigma_epsilon_sym
        - (
            2 * chi0 / (1 + delta_U)
            + 4 * epsilon_W * delta_U / (11 * (1 - epsilon) * (1 + delta_U) ** 2)
        )
        * Sigma_delta_sym
    )

    expect_zero("Xi_1 - (A_tr Sigma_tr + Sigma_nt)", Xi1_formula - (A_tr * Sigma_tr + Sigma_nt))

    Lambda0 = sp.symbols("Lambda_0", positive=True, real=True)
    Rtarget_t = sp.simplify(Lambda0 * (1 - epsilon_eta_t) / T2_t)
    R1 = first_log_drift(Rtarget_t, s, lam)

    expect_zero(
        "R_1 + Xi_1 + epsilon_eta Sigma_eta/(1-epsilon_eta)",
        R1 + Xi1 + epsilon_eta / (1 - epsilon_eta) * Sigma_eta_sym,
    )

    print("Sigma_nt =")
    sp.pprint(Sigma_nt)
    print("Xi_1 =")
    sp.pprint(sp.simplify(A_tr * Sigma_tr + Sigma_nt))
    print("R_1 + Xi_1 =")
    sp.pprint(sp.simplify(-epsilon_eta / (1 - epsilon_eta) * Sigma_eta_sym))

    # ------------------------------------------------------------------
    # V. Exact inverse reconstruction formulas
    # ------------------------------------------------------------------
    subbanner("XI.5 — Exact inverse reconstruction formulas")

    Sigma_tr_recon = sp.simplify(-((1 + chi0) * (1 + delta_U) * (1 + chi0 + delta_U) / (chi0 * delta_U)) * Theta1)
    Sigma_nt_recon = sp.simplify(Xi1 + 2 * (1 + chi0 + delta_U) / delta_U * Theta1)
    Sigma_eta_recon = sp.simplify(-((1 - epsilon_eta) / epsilon_eta) * (R1 + Xi1))

    expect_zero("Sigma_tr - reconstructed Sigma_tr", Sigma_tr - Sigma_tr_recon)
    expect_zero("Sigma_nt - reconstructed Sigma_nt", Sigma_nt - Sigma_nt_recon)
    expect_zero("Sigma_eta - reconstructed Sigma_eta", Sigma_eta_sym - Sigma_eta_recon)

    print("Sigma_tr(Theta_1) =")
    sp.pprint(Sigma_tr_recon)
    print("Sigma_nt(Theta_1, Xi_1) =")
    sp.pprint(Sigma_nt_recon)
    print("Sigma_eta(R_1, Xi_1) =")
    sp.pprint(Sigma_eta_recon)

    # ------------------------------------------------------------------
    # VI. Quartic anomaly gates in the branch-adapted coordinates
    # ------------------------------------------------------------------
    subbanner("XI.6 — Quartic anomaly gates in the branch-adapted coordinates")

    Lambda1 = sp.Float("0.279605891931464")
    tracking_rigid_gate = sp.Eq(Sigma_nt, Lambda1)
    dressing_rigid_gate = sp.Eq(Sigma_eta_sym, -((1 - epsilon_eta) / epsilon_eta) * Lambda1)

    print("Full quartic gate:")
    sp.pprint(sp.Eq(A_tr * Sigma_tr + Sigma_nt, Lambda1))
    print("If exact tracking rigidity holds (Sigma_tr = 0):")
    sp.pprint(tracking_rigid_gate)
    print("If tracking rigidity and selected-branch rigidity (R_1 = 0) both hold:")
    sp.pprint(dressing_rigid_gate)

    # ------------------------------------------------------------------
    # VII. Final theorem ledger
    # ------------------------------------------------------------------
    banner("FINAL STEP-11 LEDGER")
    print("1. The direct microscopic coherent-kernel log drifts are:")
    print("   Sigma_chi     = gamma_1 + c_1 - kappa_U")
    print("   Sigma_delta   = tau_1 - kappa_U")
    print("   Sigma_Z       = 2 lambda_1 + mu_1 - kappa_eta - 2 kappa_W")
    print("   Sigma_epsilon = 2 gamma_1 + 2 lambda_1 - kappa_U - kappa_W")
    print("   Sigma_eta     = 2 c_1 - kappa_U - kappa_eta")
    print()
    print("2. The coherent grouped defect depends only on four direct slippages:")
    print("   Sigma_Z, Sigma_chi, Sigma_epsilon, Sigma_delta.")
    print("   The support lane is absent already at this level.")
    print()
    print("3. The whole coherent weak-axisymmetric problem collapses to three")
    print("   branch-adapted scalars:")
    print("      Sigma_tr = (1+chi_0) Sigma_delta + (1+delta_U) Sigma_chi")
    print("      Sigma_nt = genuine nontracking transfer-shape slippage")
    print("      Sigma_eta = dressing slippage")
    print()
    print("4. Exact triangular normal form:")
    print("   Theta_1 = - C_tr Sigma_tr")
    print("   Xi_1    = A_tr Sigma_tr + Sigma_nt")
    print("   R_1+Xi_1 = - [epsilon_eta/(1-epsilon_eta)] Sigma_eta")
    print()
    print("5. Exact inverse reconstruction:")
    print("   Sigma_tr = -[(1+chi_0)(1+delta_U)(1+chi_0+delta_U)/(chi_0 delta_U)] Theta_1")
    print("   Sigma_nt = Xi_1 + 2(1+chi_0+delta_U) Theta_1 / delta_U")
    print("   Sigma_eta = -[(1-epsilon_eta)/epsilon_eta] (R_1 + Xi_1)")
    print()
    print("6. Quartic anomaly gate:")
    print("   A_tr Sigma_tr + Sigma_nt = Lambda_1 ≈ 0.279605891931464")
    print("   so on an exactly tracking-rigid branch the missing g-2 layer is")
    print("   purely the nontracking slippage Sigma_nt.")
    print()
    print("So the next continuation point is now smaller again:")
    print("  compute whether the actual moving-throat branch preserves or drives")
    print("  the three branch-adapted scalars (Sigma_tr, Sigma_nt, Sigma_eta),")
    print("  rather than chasing a larger raw-placement variable set.")


if __name__ == "__main__":
    step11_microscopic_slippage_normal_form()
