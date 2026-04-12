
#!/usr/bin/env python3
"""
Step 12 for the g-2 rebuild from the moving-throat PDE base.

What this script does
---------------------
1. Starts from the Step-11 branch-adapted slippages
       Sigma_tr, Sigma_nt, Sigma_eta
   on the coherent local D/N tracking branch.
2. Proves that these are exactly the logarithmic drifts of three direct
   microscopic monomials:
       C_tr,*,
       C_nt,*,
       epsilon_eta.
3. Rewrites the full first grouped weak-axisymmetric zero-defect condition as
       delta ln C_tr,* = delta ln C_nt,* = delta ln epsilon_eta = 0.
4. Builds the exact 3x8 monomial-drift matrix M_* from the microscopic drift
   vector
       (lambda_1, c_1, gamma_1, kappa_U, kappa_eta, kappa_W, mu_1, tau_1),
   verifies the closed-form matrix, and checks the useful nonzero minor
       det M_*^(tau_1, kappa_eta, mu_1) = 1 + chi_{0,*}.
5. Solves the resulting compatibility ledger explicitly for
       tau_1, kappa_eta, mu_1
   in terms of the remaining microscopic drifts.
6. Rewrites the quartic anomaly gate directly in the monomial coordinates.

Interpretation
--------------
Step 11 reduced the coherent weak-axisymmetric problem to three exact
branch-adapted scalars. Step 12 sharpens that one step further: the live
coordinates are not abstract slippages anymore, but three direct microscopic
kernel monomials. This is the natural bridge into the exact similarity-orbit /
quotient language of the moving-throat notes.
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


def main() -> None:
    banner("STEP 12 — DIRECT MICROSCOPIC MONOMIALS AND COMPATIBILITY LEDGER")

    # ------------------------------------------------------------------
    # I. Microscopic variables and direct coherent-kernel ratios
    # ------------------------------------------------------------------
    subbanner("XII.1 — Microscopic coherent-kernel variables")

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
    delta_A = sp.simplify(sp.pi**2 * T_U_A / (L**2 * K_U_A))
    eps_eta_A = sp.simplify(c_A**2 / (K_U_A * K_eta_A))
    eps_W_A = sp.simplify(gamma_A**2 * lambda_W_A**2 * sigma / (K_U_A * K_W_A))
    Zratio_A = sp.simplify(lambda_W_A**2 * mu_W_A / (K_eta_A * K_W_A**2))

    Sigma_chi = first_log_drift(chi_A, s, lam)
    Sigma_delta = first_log_drift(delta_A, s, lam)
    Sigma_eta = first_log_drift(eps_eta_A, s, lam)
    Sigma_eps = first_log_drift(eps_W_A, s, lam)
    Sigma_Z = first_log_drift(Zratio_A, s, lam)

    expect_zero(
        "Sigma_chi - (gamma_1 + c_1 - kappa_U)",
        Sigma_chi - (gamma1 + c1 - kappa_U),
    )
    expect_zero(
        "Sigma_delta - (tau_1 - kappa_U)",
        Sigma_delta - (tau1 - kappa_U),
    )
    expect_zero(
        "Sigma_eta - (2 c_1 - kappa_U - kappa_eta)",
        Sigma_eta - (2 * c1 - kappa_U - kappa_eta),
    )
    expect_zero(
        "Sigma_eps - (2 gamma_1 + 2 lambda_1 - kappa_U - kappa_W)",
        Sigma_eps - (2 * gamma1 + 2 * lambda1 - kappa_U - kappa_W),
    )
    expect_zero(
        "Sigma_Z - (2 lambda_1 + mu_1 - kappa_eta - 2 kappa_W)",
        Sigma_Z - (2 * lambda1 + mu1 - kappa_eta - 2 * kappa_W),
    )

    chi0 = sp.simplify(gamma * c_etaU / K_U)
    delta_U = sp.simplify(sp.pi**2 * T_U / (L**2 * K_U))
    eps_eta = sp.simplify(c_etaU**2 / (K_U * K_eta))
    eps_W = sp.simplify(gamma**2 * lambda_W**2 * sigma / (K_U * K_W))
    Zratio = sp.simplify(lambda_W**2 * mu_W / (K_eta * K_W**2))
    eps = sp.simplify(eps_W * (1 - sp.Rational(2, 11) * delta_U / (1 + delta_U)))

    print("chi_0 =")
    sp.pprint(chi0)
    print("delta_U =")
    sp.pprint(delta_U)
    print("epsilon_eta =")
    sp.pprint(eps_eta)
    print("epsilon_W =")
    sp.pprint(eps_W)
    print("Z_W/Omega_W^2 =")
    sp.pprint(Zratio)

    # ------------------------------------------------------------------
    # II. Direct microscopic tracking monomial
    # ------------------------------------------------------------------
    subbanner("XII.2 — Direct microscopic tracking monomial")

    chi_star, delta_star = sp.symbols("chi0_star deltaU_star", positive=True, real=True)

    Ctr_A = sp.simplify(chi_A ** (1 + delta_star) * delta_A ** (1 + chi_star))
    dlog_Ctr = first_log_drift(Ctr_A, s, lam)
    Sigma_tr = sp.simplify((1 + chi_star) * Sigma_delta + (1 + delta_star) * Sigma_chi)

    expect_zero("delta ln C_tr,* - Sigma_tr", dlog_Ctr - Sigma_tr)

    Ctr_explicit = sp.simplify(
        (gamma * c_etaU / K_U) ** (1 + delta_star)
        * (sp.pi**2 * T_U / (L**2 * K_U)) ** (1 + chi_star)
    )

    print("C_tr,* =")
    sp.pprint(Ctr_explicit)
    print("delta ln C_tr,* =")
    sp.pprint(dlog_Ctr)

    # ------------------------------------------------------------------
    # III. Direct microscopic nontracking monomial
    # ------------------------------------------------------------------
    subbanner("XII.3 — Direct microscopic nontracking monomial")

    epsW_star, eps_star = sp.symbols("epsilonW_star epsilon_star", positive=True, real=True)
    E_star = sp.simplify(
        2 * epsW_star / (1 - eps_star) * (11 + 9 * delta_star) / (11 * (1 + delta_star))
    )
    F_star = sp.simplify(
        2 * chi_star / (1 + delta_star)
        + 4 * epsW_star * delta_star / (11 * (1 - eps_star) * (1 + delta_star) ** 2)
    )

    Cnt_A = sp.simplify(Zratio_A * eps_W_A ** E_star * delta_A ** (-F_star))
    dlog_Cnt = first_log_drift(Cnt_A, s, lam)
    Sigma_nt = sp.simplify(Sigma_Z + E_star * Sigma_eps - F_star * Sigma_delta)

    expect_zero("delta ln C_nt,* - Sigma_nt", dlog_Cnt - Sigma_nt)

    Cnt_explicit = sp.simplify(
        (lambda_W**2 * mu_W / (K_eta * K_W**2))
        * (gamma**2 * lambda_W**2 * sigma / (K_U * K_W)) ** E_star
        * (sp.pi**2 * T_U / (L**2 * K_U)) ** (-F_star)
    )

    print("E_* =")
    sp.pprint(E_star)
    print("F_* =")
    sp.pprint(F_star)
    print("C_nt,* =")
    sp.pprint(Cnt_explicit)
    print("delta ln C_nt,* =")
    sp.pprint(dlog_Cnt)

    # ------------------------------------------------------------------
    # IV. Exact zero-defect monomial ledger
    # ------------------------------------------------------------------
    subbanner("XII.4 — Exact zero-defect monomial ledger")

    expect_zero("delta ln epsilon_eta - Sigma_eta", Sigma_eta - first_log_drift(eps_eta_A, s, lam))

    print("Zero-defect theorem:")
    print("  Theta_1 = Xi_1 = R_1 = 0")
    print("  iff")
    print("  delta ln C_tr,* = delta ln C_nt,* = delta ln epsilon_eta = 0")

    # ------------------------------------------------------------------
    # V. Exact monomial-drift matrix M_*
    # ------------------------------------------------------------------
    subbanner("XII.5 — Exact monomial-drift matrix")

    xvec = sp.Matrix([lambda1, c1, gamma1, kappa_U, kappa_eta, kappa_W, mu1, tau1])
    qvec = sp.Matrix([dlog_Ctr, dlog_Cnt, Sigma_eta])

    M = sp.simplify(qvec.jacobian(xvec))

    M_expected = sp.Matrix(
        [
            [
                0,
                1 + delta_star,
                1 + delta_star,
                -(2 + chi_star + delta_star),
                0,
                0,
                0,
                1 + chi_star,
            ],
            [
                2 * (1 + E_star),
                0,
                2 * E_star,
                F_star - E_star,
                -1,
                -(2 + E_star),
                1,
                -F_star,
            ],
            [
                0,
                2,
                0,
                -1,
                -1,
                0,
                0,
                0,
            ],
        ]
    )

    expect_zero("M_* - expected matrix", M - M_expected)

    print("M_* =")
    sp.pprint(M)

    minor = M.extract([0, 1, 2], [7, 4, 6])
    minor_det = sp.simplify(minor.det())
    expect_zero("det M_*^(tau_1,kappa_eta,mu_1) - (1 + chi0_star)", minor_det - (1 + chi_star))

    print("det M_*^(tau_1, kappa_eta, mu_1) =")
    sp.pprint(minor_det)
    print("So rank(M_*) = 3 and dim ker(M_*) = 5 on the physical branch.")

    # ------------------------------------------------------------------
    # VI. Exact microscopic compatibility ledger
    # ------------------------------------------------------------------
    subbanner("XII.6 — Exact microscopic compatibility equations")

    compatibility = [
        sp.Eq(dlog_Ctr, 0),
        sp.Eq(dlog_Cnt, 0),
        sp.Eq(Sigma_eta, 0),
    ]
    sol = sp.solve(compatibility, (tau1, kappa_eta, mu1), dict=True)[0]

    tau_sol = sp.simplify(sol[tau1])
    keta_sol = sp.simplify(sol[kappa_eta])
    mu_sol = sp.simplify(sol[mu1])

    expect_zero("tracking compatibility residual", dlog_Ctr.subs(sol))
    expect_zero("nontracking compatibility residual", dlog_Cnt.subs(sol))
    expect_zero("dressing compatibility residual", Sigma_eta.subs(sol))

    print("tau_1 =")
    sp.pprint(tau_sol)
    print("kappa_eta =")
    sp.pprint(keta_sol)
    print("mu_1 =")
    sp.pprint(mu_sol)

    # ------------------------------------------------------------------
    # VII. Quartic anomaly gate in the monomial coordinates
    # ------------------------------------------------------------------
    subbanner("XII.7 — Quartic anomaly gate in the monomial coordinates")

    A_tr = sp.simplify(2 * chi_star / ((1 + chi_star) * (1 + delta_star)))
    Xi_mono = sp.simplify(A_tr * dlog_Ctr + dlog_Cnt)
    Xi_raw = sp.simplify(A_tr * Sigma_tr + Sigma_nt)

    expect_zero("Xi_1 monomial form - Xi_1 raw form", Xi_mono - Xi_raw)

    Lambda1 = sp.symbols("Lambda_1", real=True)
    print("Xi_1 =")
    sp.pprint(Xi_mono)
    print("Quartic anomaly gate:")
    print("  A_tr delta ln C_tr,* + delta ln C_nt,* = Lambda_1")
    print("with")
    print("  A_tr =")
    sp.pprint(A_tr)
    print("Tracking-rigid specialization:")
    print("  delta ln C_tr,* = 0  =>  delta ln C_nt,* = Lambda_1")

    banner("STEP 12 COMPLETE")


if __name__ == "__main__":
    main()
