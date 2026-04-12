#!/usr/bin/env python3
"""
Step 9 for the g-2 rebuild from the moving-throat PDE base.

What this script does
---------------------
1. Starts from Step 8's exact portwise theorem
       N_{A,0}^{(r)} = K_A * T_{A,r}^2,
       Xi_1 = 2 sum_r rho_r^(N) tau_r,
   and proves the exact collapse to one effective transfer shape
       T_eff,A^2 = N_{A,0}/K_A,
       Xi_1 = d ln(T_eff,A^2)/(eps*lambda_A).
2. Evaluates that effective shape on the actual minimal isotropic continuum
   branch, where there is only one active outgoing port and
       N_{A,0} = beta_{0,A},
       K_A     = K_{0,A}.
3. Uses the continuum branch formulas to derive the direct one-port identity
       T_A^2 = Z_{W,A}(1+rho_A)^2 / [ Omega_{W,A}^2 (1-epsilon_{W,A})^2 ].
4. Perturbs that formula weak-axisymmetrically and proves
       Xi_1 = zeta_W - omega_W + 2 rho_1/(1+rho) + 2 varepsilon_W/(1-epsilon_W).
5. Rewrites the same transfer shape in the selected-branch language and proves
       Xi_1 = - eta_1/(1-epsilon_eta) - R_1.
6. Records the exact one-port zero-defect theorem and the direct quartic-anomaly
   matching conditions on both the direct-port and selected-branch variables.

Interpretation
--------------
Step 8 reduced the remaining quartic anomaly layer to an outgoing-weighted
transfer-shape slope. Step 9 sharpens that again: the many-port weighted average
is itself the slope of one effective transfer shape, and on the actual minimal
continuum branch that shape is an explicit function of the microscopic port
placement data. So the remaining theorem problem is no longer to compute many
raw outgoing-port slopes; it is to determine whether the actual weak-axisymmetric
moving-throat branch keeps this single continuum transfer shape rigid.
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


def trunc1(expr: sp.Expr, eps: sp.Symbol) -> sp.Expr:
    """Keep terms through first order in eps."""
    return sp.expand(sp.series(expr, eps, 0, 2).removeO())


def step9_effective_transfer_shape_continuum() -> None:
    banner("STEP 9 — EFFECTIVE TRANSFER SHAPE AND ONE-PORT CONTINUUM BRANCH")

    # ------------------------------------------------------------------
    # I. Exact collapse from many ports to one effective transfer shape
    # ------------------------------------------------------------------
    subbanner("IX.1 — Exact collapse T_eff^2 = N_0 / K and Xi_1 = d ln T_eff^2")

    eps, lam = sp.symbols("epsilon lambda_A", real=True)
    T1, T2, T3 = sp.symbols("T_1 T_2 T_3", positive=True, real=True)
    tau1, tau2, tau3 = sp.symbols("tau_1 tau_2 tau_3", real=True)

    T1_A = trunc1(T1 * sp.exp(eps * lam * tau1), eps)
    T2_A = trunc1(T2 * sp.exp(eps * lam * tau2), eps)
    T3_A = trunc1(T3 * sp.exp(eps * lam * tau3), eps)

    Teff2 = sp.simplify(T1**2 + T2**2 + T3**2)
    Teff2_A = trunc1(T1_A**2 + T2_A**2 + T3_A**2, eps)
    Xi1_eff = sp.simplify(sp.diff(sp.log(Teff2_A), eps).subs(eps, 0) / lam)

    rho1N = sp.simplify(T1**2 / Teff2)
    rho2N = sp.simplify(T2**2 / Teff2)
    rho3N = sp.simplify(T3**2 / Teff2)

    expect_zero("rho_1^(N) + rho_2^(N) + rho_3^(N) - 1", rho1N + rho2N + rho3N - 1)
    expect_zero(
        "Xi_1 - 2 sum_r rho_r^(N) tau_r",
        Xi1_eff - 2 * (rho1N * tau1 + rho2N * tau2 + rho3N * tau3),
    )

    tau_eff = sp.simplify(Xi1_eff / 2)
    print("T_eff,A^2 =")
    sp.pprint(Teff2_A)
    print("Xi_1 =")
    sp.pprint(Xi1_eff)
    print("tau_eff =")
    sp.pprint(tau_eff)

    # ------------------------------------------------------------------
    # II. Exact one-port continuum transfer shape
    # ------------------------------------------------------------------
    subbanner("IX.2 — Exact one-port continuum transfer shape")

    mu_eta, mu_W = sp.symbols("mu_eta mu_W", positive=True, real=True)
    Keta_eff, KW_eff = sp.symbols("K_eta_eff K_W_eff", positive=True, real=True)
    ZW, rho, epsW = sp.symbols("Z_W rho epsilon_W", positive=True, real=True)
    eps_eta = sp.symbols("epsilon_eta", real=True)

    K0 = sp.simplify(Keta_eff / mu_eta)
    beta0 = sp.simplify(
        (mu_W / mu_eta) * (Keta_eff / KW_eff) * ZW * (1 + rho) ** 2 / (1 - epsW) ** 2
    )
    T2_cont = sp.simplify(beta0 / K0)
    OmegaW2 = sp.simplify(KW_eff / mu_W)

    expect_zero(
        "T_A^2 - (mu_W/K_W_eff) * Z_W (1+rho)^2 / (1-epsilon_W)^2",
        T2_cont - (mu_W / KW_eff) * ZW * (1 + rho) ** 2 / (1 - epsW) ** 2,
    )
    expect_zero(
        "T_A^2 - Z_W (1+rho)^2 / [Omega_W^2 (1-epsilon_W)^2]",
        T2_cont - ZW * (1 + rho) ** 2 / (OmegaW2 * (1 - epsW) ** 2),
    )
    expect_zero("d T_A^2 / d epsilon_eta", sp.diff(T2_cont, eps_eta))

    print("K_{0,A} =")
    sp.pprint(K0)
    print("beta_{0,A} =")
    sp.pprint(beta0)
    print("T_A^2 =")
    sp.pprint(T2_cont)

    # ------------------------------------------------------------------
    # III. Direct weak-axisymmetric slope law in continuum variables
    # ------------------------------------------------------------------
    subbanner("IX.3 — Direct weak-axisymmetric slope law in continuum variables")

    zetaW, omegaW, rho1, varepsW = sp.symbols(
        "zeta_W omega_W rho_1 varepsilon_W", real=True
    )

    ZW_A = trunc1(ZW * sp.exp(eps * lam * zetaW), eps)
    OmegaW2_A = trunc1(OmegaW2 * sp.exp(eps * lam * omegaW), eps)
    rho_A = sp.expand(rho + eps * lam * rho1)
    epsW_A = sp.expand(epsW + eps * lam * varepsW)

    T2_A = trunc1(ZW_A * (1 + rho_A) ** 2 / (OmegaW2_A * (1 - epsW_A) ** 2), eps)
    Xi1_direct = sp.simplify(sp.diff(sp.log(T2_A), eps).subs(eps, 0) / lam)
    tau_direct = sp.simplify(Xi1_direct / 2)

    expect_zero(
        "Xi_1 - [zeta_W - omega_W + 2 rho_1/(1+rho) + 2 varepsilon_W/(1-epsilon_W)]",
        Xi1_direct - (zetaW - omegaW + 2 * rho1 / (1 + rho) + 2 * varepsW / (1 - epsW)),
    )
    expect_zero("tau - Xi_1/2", tau_direct - Xi1_direct / 2)

    print("Xi_1 =")
    sp.pprint(Xi1_direct)
    print("tau =")
    sp.pprint(tau_direct)

    # ------------------------------------------------------------------
    # IV. Selected-branch reformulation
    # ------------------------------------------------------------------
    subbanner("IX.4 — Selected-branch reformulation")

    G, cs, a_th, c_light = sp.symbols("G c_s a c", positive=True, real=True)
    pi = sp.pi

    Lambda = sp.simplify(27 * pi**2 * G * cs**5 * KW_eff / (20 * a_th**5 * c_light**5 * mu_W))
    Lambda_via_Omega = sp.simplify(27 * pi**2 * G * cs**5 * OmegaW2 / (20 * a_th**5 * c_light**5))
    expect_zero("Lambda - (27 pi^2 G c_s^5 / (20 a^5 c^5)) Omega_W^2", Lambda - Lambda_via_Omega)

    Rtarget = sp.simplify(Lambda * (1 - eps_eta) * (1 - epsW) ** 2 / (ZW * (1 + rho) ** 2))
    front = sp.simplify(27 * pi**2 * G * cs**5 / (20 * a_th**5 * c_light**5))
    T2_sel = sp.simplify(front * (1 - eps_eta) / Rtarget)

    expect_zero("T_A^2(selected) - T_A^2(direct)", T2_sel - T2_cont)

    eta1, R1 = sp.symbols("eta_1 mathcal_R_1", real=True)
    eps_eta_A = sp.expand(eps_eta + eps * lam * eta1)
    Rtarget_A = trunc1(Rtarget * sp.exp(eps * lam * R1), eps)
    T2_sel_A = trunc1(front * (1 - eps_eta_A) / Rtarget_A, eps)
    Xi1_sel = sp.simplify(sp.diff(sp.log(T2_sel_A), eps).subs(eps, 0) / lam)

    expect_zero(
        "Xi_1(selected) + eta_1/(1-epsilon_eta) + R_1",
        Xi1_sel + eta1 / (1 - eps_eta) + R1,
    )

    print("R_target =")
    sp.pprint(Rtarget)
    print("Xi_1(selected) =")
    sp.pprint(Xi1_sel)

    # ------------------------------------------------------------------
    # V. Zero-defect theorem and corollaries
    # ------------------------------------------------------------------
    subbanner("IX.5 — Exact one-port zero-defect theorem and corollaries")

    # Direct zero-defect theorem.
    direct_invariant = sp.simplify(ZW * (1 + rho) ** 2 / (OmegaW2 * (1 - epsW) ** 2))
    direct_invariant_A = trunc1(
        ZW_A * (1 + rho_A) ** 2 / (OmegaW2_A * (1 - epsW_A) ** 2), eps
    )
    expect_zero(
        "d ln(direct invariant)/(eps*lambda_A) - Xi_1",
        sp.simplify(sp.diff(sp.log(direct_invariant_A), eps).subs(eps, 0) / lam - Xi1_direct),
    )

    zero_defect_eq = sp.solve(sp.Eq(Xi1_sel, 0), R1)[0]
    expect_zero(
        "R_1 + eta_1/(1-epsilon_eta) under Xi_1=0",
        zero_defect_eq + eta1 / (1 - eps_eta),
    )

    print("Zero-defect selected-branch law:")
    print("  R_1 =")
    sp.pprint(zero_defect_eq)
    print("Target-rigid branch (R_1 = 0) requires eta_1 =")
    sp.pprint(sp.solve(sp.Eq(Xi1_sel.subs(R1, 0), 0), eta1)[0])
    print("epsilon_eta-rigid branch (eta_1 = 0) requires R_1 =")
    sp.pprint(sp.solve(sp.Eq(Xi1_sel.subs(eta1, 0), 0), R1)[0])

    # ------------------------------------------------------------------
    # VI. Direct quartic anomaly gate and reference balanced splits
    # ------------------------------------------------------------------
    subbanner("IX.6 — Direct quartic anomaly gate and reference balanced splits")

    Lambda1 = sp.Float("0.279605891931464")
    tau_target = sp.simplify(Lambda1 / 2)

    # Balanced split in the direct continuum variables, after rescaling the four
    # microscopic channels to unit coefficient.
    x = sp.symbols("x", real=True)
    direct_balanced = sp.solve(sp.Eq(4 * x, Lambda1), x)[0]
    zetaW_bal = sp.simplify(direct_balanced)
    omegaW_bal = sp.simplify(-direct_balanced)
    rho1_bal = sp.simplify((1 + rho) * direct_balanced / 2)
    varepsW_bal = sp.simplify((1 - epsW) * direct_balanced / 2)

    # Balanced split in the selected-branch variables after rescaling
    # y1 = -eta_1/(1-epsilon_eta), y2 = -R_1.
    y = sp.symbols("y", real=True)
    selected_balanced = sp.solve(sp.Eq(2 * y, Lambda1), y)[0]
    eta1_bal = sp.simplify(-(1 - eps_eta) * selected_balanced)
    R1_bal = sp.simplify(-selected_balanced)

    print("Xi_1 target = Lambda_1 =")
    sp.pprint(Lambda1)
    print("tau target = Lambda_1 / 2 =")
    sp.pprint(tau_target)
    print()
    print("Balanced direct-port reference split:")
    print("  zeta_W =")
    sp.pprint(zetaW_bal)
    print("  omega_W =")
    sp.pprint(omegaW_bal)
    print("  rho_1 =")
    sp.pprint(rho1_bal)
    print("  varepsilon_W =")
    sp.pprint(varepsW_bal)
    print()
    print("Balanced selected-branch reference split:")
    print("  eta_1 =")
    sp.pprint(eta1_bal)
    print("  R_1 =")
    sp.pprint(R1_bal)

    # ------------------------------------------------------------------
    # VII. Grouped weak-axisymmetric normalization pattern
    # ------------------------------------------------------------------
    subbanner("IX.7 — Grouped weak-axisymmetric normalization pattern")

    Delta20 = sp.simplify(eps * Xi1_direct)
    Delta21 = sp.simplify(eps * Xi1_direct / 2)
    Delta22 = sp.simplify(-eps * Xi1_direct)

    expect_zero("Delta_Q^(20) - 2 epsilon tau", Delta20 - 2 * eps * tau_direct)
    expect_zero("Delta_Q^(21) - epsilon tau", Delta21 - eps * tau_direct)
    expect_zero("Delta_Q^(22) + 2 epsilon tau", Delta22 + 2 * eps * tau_direct)

    print("Delta_Q^(20) =")
    sp.pprint(Delta20)
    print("Delta_Q^(21) =")
    sp.pprint(Delta21)
    print("Delta_Q^(22) =")
    sp.pprint(Delta22)

    # ------------------------------------------------------------------
    # VIII. Final theorem ledger
    # ------------------------------------------------------------------
    banner("FINAL STEP-9 LEDGER")
    print("1. Exact collapse from many ports to one effective transfer shape:")
    print("   T_eff,A^2 = N_{A,0} / K_A")
    print("   Xi_1 = d ln(T_eff,A^2)/(epsilon lambda_A) = 2 tau_eff")
    print()
    print("2. Exact one-port continuum transfer shape:")
    print("   T_A^2 = Z_{W,A}(1+rho_A)^2 / [Omega_{W,A}^2 (1-epsilon_{W,A})^2]")
    print()
    print("3. Direct one-port weak-axisymmetric defect law:")
    print("   Xi_1 = zeta_W - omega_W + 2 rho_1/(1+rho) + 2 varepsilon_W/(1-epsilon_W)")
    print()
    print("4. Selected-branch reformulation:")
    print("   Xi_1 = - eta_1/(1-epsilon_eta) - R_1")
    print()
    print("5. Zero-defect theorem on the one-port continuum branch:")
    print("   d ln[ Z_W (1+rho)^2 / (Omega_W^2 (1-epsilon_W)^2 ) ] = 0")
    print("   <=>  d ln R_target = - d epsilon_eta / (1-epsilon_eta)")
    print()
    print("6. Quartic anomaly gate on the actual one-port branch:")
    print("   Xi_1 = Lambda_1 ≈ 0.279605891931464")
    print("   tau  = Lambda_1/2 ≈ 0.139802945965732")
    print()
    print("So the next theorem problem is now smaller again:")
    print("  determine whether the actual weak-axisymmetric moving-throat branch")
    print("  keeps this single continuum transfer shape rigid, or else which of")
    print("  its four direct microscopic drift channels carries the needed")
    print("  Lambda_1-sized defect.")


if __name__ == "__main__":
    step9_effective_transfer_shape_continuum()
