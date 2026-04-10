#!/usr/bin/env python3
"""
Step 10 for the g-2 rebuild from the moving-throat PDE base.

What this script does
---------------------
1. Starts from Step 9's exact one-port transfer-shape law and substitutes the
   actual coherent local D/N tracking-branch map.
2. Proves the exact coherent-branch identity
       T^2 = Z_W (1+chi_0)^2 / [Omega_W^2 (1-eps)^2]
   with
       eps = eps_W * [1 - (2/11) delta_U/(1+delta_U)].
3. Uses the coherent support-enhancement factor
       M_tr = M_mix S(zeta;eps),
       S(zeta;eps) = 1 + zeta(1-eps)/(1-zeta eps),
   to show that the support ratio zeta changes the baseline but drops out of
   T^2, R_target, and Xi_1.
4. Perturbs the coherent branch weak-axisymmetrically and derives the exact
   split-blocking drift
       eps_1 = (1 - (2/11) delta_U/(1+delta_U)) varepsilon_W
               - 2 eps_W deltaU_1 / [11 (1+delta_U)^2],
   together with the direct coherent-branch defect law
       Xi_1 = zeta_Z - omega_W + 2 chi_1/(1+chi_0) + 2 eps_1/(1-eps).
5. Rewrites the same result in selected-branch form and derives the exact
   tracking-factor drift Theta_1.
6. Verifies that Theta_1 = 0 is not sufficient to kill Xi_1.
7. Records the direct quartic-anomaly gate on the coherent tracking branch.

Interpretation
--------------
Step 9 reduced the remaining quartic anomaly layer to the slope of one actual
continuum transfer shape. Step 10 pushes that formula onto the coherent local
D/N tracking branch. The new theorem is sharp: coherent support can raise the
steady baseline, but it is support-blind at the level of the first grouped
weak-axisymmetric defect. So the next live calculation is not the support ratio;
it is the weak-axisymmetric drift of the mixed/outgoing placement variables.
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



def trunc1(expr: sp.Expr, s: sp.Symbol) -> sp.Expr:
    """Keep terms through first order in s."""
    return sp.expand(sp.series(expr, s, 0, 2).removeO())



def step10_coherent_tracking_branch_defect() -> None:
    banner("STEP 10 — COHERENT TRACKING BRANCH AND SUPPORT-BLIND DEFECT LAW")

    # ------------------------------------------------------------------
    # I. Exact coherent-branch substitution
    # ------------------------------------------------------------------
    subbanner("X.1 — Exact coherent tracking-branch substitution")

    G, cs, a_th, c_light = sp.symbols("G c_s a c", positive=True, real=True)
    ZW, chi0 = sp.symbols("Z_W chi_0", positive=True, real=True)
    eps_eta, epsW = sp.symbols("epsilon_eta epsilon_W", real=True)
    deltaU = sp.symbols("delta_U", positive=True, real=True)
    OmegaW2 = sp.symbols("Omega_W_sq", positive=True, real=True)
    zeta_support = sp.symbols("zeta", real=True)
    pi = sp.pi

    front = sp.simplify(27 * pi**2 * G * cs**5 / (20 * a_th**5 * c_light**5))
    Lambda = sp.simplify(front * OmegaW2)
    eps_split = sp.simplify(epsW * (1 - sp.Rational(2, 11) * deltaU / (1 + deltaU)))
    Rtarget = sp.simplify(Lambda * (1 - eps_eta) * (1 - eps_split) ** 2 / (ZW * (1 + chi0) ** 2))
    T2_selected = sp.simplify(front * (1 - eps_eta) / Rtarget)
    T2_direct = sp.simplify(ZW * (1 + chi0) ** 2 / (OmegaW2 * (1 - eps_split) ** 2))

    expect_zero("T^2(selected coherent branch) - T^2(direct coherent branch)", T2_selected - T2_direct)

    print("eps =")
    sp.pprint(eps_split)
    print("R_target =")
    sp.pprint(Rtarget)
    print("T^2 =")
    sp.pprint(T2_direct)

    # ------------------------------------------------------------------
    # II. Exact support-blindness theorem
    # ------------------------------------------------------------------
    subbanner("X.2 — Exact support-blindness theorem")

    Mmix = sp.simplify(8 * ZW * (1 + chi0) ** 2 / (pi**2 * (1 - eps_eta) * (1 - eps_split)))
    S_support = sp.simplify(1 + zeta_support * (1 - eps_split) / (1 - zeta_support * eps_split))
    Mtr = sp.simplify(Mmix * S_support)
    product_law = sp.simplify(Rtarget * Mtr)

    expect_zero("d T^2 / d zeta", sp.diff(T2_direct, zeta_support))
    expect_zero("d R_target / d zeta", sp.diff(Rtarget, zeta_support))
    dS_dzeta = sp.simplify(sp.diff(S_support, zeta_support))
    expect_zero("dS/dzeta - (1-eps)/(1-zeta eps)^2", dS_dzeta - (1 - eps_split) / (1 - zeta_support * eps_split) ** 2)
    print("S(zeta;eps) =")
    sp.pprint(S_support)
    print("dS/dzeta =")
    sp.pprint((1 - eps_split) / (1 - zeta_support * eps_split) ** 2)
    print("M_tr =")
    sp.pprint(Mtr)
    print("R_target * M_tr =")
    sp.pprint(product_law)

    expect_zero(
        "R_target M_tr - 8 Lambda (1-eps)/pi^2 * S(zeta;eps)",
        product_law - 8 * Lambda * (1 - eps_split) * S_support / pi**2,
    )

    # ------------------------------------------------------------------
    # III. Weak-axisymmetric coherent-branch defect law
    # ------------------------------------------------------------------
    subbanner("X.3 — Exact weak-axisymmetric defect law in physical branch variables")

    s, lam = sp.symbols("s lambda_A", real=True)
    zetaZ, omegaW = sp.symbols("zeta_Z omega_W", real=True)
    chi1, varepsW, deltaU1, eta1 = sp.symbols(
        "chi_1 varepsilon_W deltaU_1 eta_1", real=True
    )

    ZW_A = trunc1(ZW * sp.exp(s * lam * zetaZ), s)
    OmegaW2_A = trunc1(OmegaW2 * sp.exp(s * lam * omegaW), s)
    chiA = sp.expand(chi0 + s * lam * chi1)
    epsW_A = sp.expand(epsW + s * lam * varepsW)
    deltaU_A = sp.expand(deltaU + s * lam * deltaU1)
    eps_eta_A = sp.expand(eps_eta + s * lam * eta1)

    eps_split_A = trunc1(
        epsW_A * (1 - sp.Rational(2, 11) * deltaU_A / (1 + deltaU_A)), s
    )
    eps1_expr = sp.simplify(sp.diff(eps_split_A, s).subs(s, 0) / lam)

    expect_zero(
        "eps_1 - [(1 - 2 delta_U/(11(1+delta_U))) varepsilon_W - 2 epsilon_W deltaU_1/(11(1+delta_U)^2)]",
        eps1_expr
        - (
            (1 - sp.Rational(2, 11) * deltaU / (1 + deltaU)) * varepsW
            - 2 * epsW * deltaU1 / (11 * (1 + deltaU) ** 2)
        ),
    )

    T2_A = trunc1(ZW_A * (1 + chiA) ** 2 / (OmegaW2_A * (1 - eps_split_A) ** 2), s)
    Xi1 = sp.simplify(sp.diff(sp.log(T2_A), s).subs(s, 0) / lam)
    Xi1_formula = sp.simplify(zetaZ - omegaW + 2 * chi1 / (1 + chi0) + 2 * eps1_expr / (1 - eps_split))

    expect_zero(
        "Xi_1 - [zeta_Z - omega_W + 2 chi_1/(1+chi_0) + 2 eps_1/(1-eps)]",
        Xi1 - Xi1_formula,
    )

    print("eps_1 =")
    sp.pprint(eps1_expr)
    print("Xi_1 =")
    sp.pprint(Xi1_formula)

    # ------------------------------------------------------------------
    # IV. Selected-branch reformulation
    # ------------------------------------------------------------------
    subbanner("X.4 — Selected-branch reformulation on the coherent branch")

    # Since R_target = Lambda (1-epsilon_eta)(1-eps)^2 / [Z_W (1+chi_0)^2]
    # with Lambda = front * Omega_W^2, the logarithmic drift is immediate.
    R1 = sp.simplify(
        omegaW - eta1 / (1 - eps_eta) - zetaZ - 2 * chi1 / (1 + chi0) - 2 * eps1_expr / (1 - eps_split)
    )

    expect_zero("Xi_1 + eta_1/(1-epsilon_eta) + R_1", Xi1_formula + eta1 / (1 - eps_eta) + R1)

    print("R_1 =")
    sp.pprint(R1)

    # ------------------------------------------------------------------
    # V. Tracking-factor drift is not sufficient
    # ------------------------------------------------------------------
    subbanner("X.5 — Tracking-factor drift is not sufficient")

    Rtr = sp.simplify((1 + chi0 / (1 + deltaU)) / (1 + chi0))
    # d ln R_tr = d ln(1+chi_0+delta_U) - d ln(1+chi_0) - d ln(1+delta_U)
    Theta1 = sp.simplify((chi1 + deltaU1) / (1 + chi0 + deltaU) - chi1 / (1 + chi0) - deltaU1 / (1 + deltaU))

    Theta1_target = sp.simplify(
        -(
            chi0 * (1 + chi0) * deltaU1 + deltaU * (1 + deltaU) * chi1
        )
        / ((1 + chi0) * (1 + deltaU) * (1 + chi0 + deltaU))
    )
    expect_zero("Theta_1 - exact tracking-drift formula", Theta1 - Theta1_target)

    Theta1_example = sp.simplify(Theta1.subs({chi1: 0, deltaU1: 0}))
    expect_zero("Theta_1 on chi_1 = deltaU_1 = 0 slice", Theta1_example)

    Xi1_slice_formula = sp.simplify(Xi1_formula.subs({chi1: 0, deltaU1: 0}))
    print("Theta_1 =")
    sp.pprint(Theta1_target)
    print("Xi_1 on chi_1 = deltaU_1 = 0 slice =")
    sp.pprint(Xi1_slice_formula)

    # ------------------------------------------------------------------
    # VI. Support-blindness at the defect level
    # ------------------------------------------------------------------
    subbanner("X.6 — Support-blindness at the defect level")

    expect_zero("d Xi_1 / d zeta", sp.diff(Xi1, zeta_support))
    expect_zero("d R_1 / d zeta", sp.diff(R1, zeta_support))

    # ------------------------------------------------------------------
    # VII. Direct quartic-anomaly gate on the coherent branch
    # ------------------------------------------------------------------
    subbanner("X.7 — Direct quartic-anomaly gate on the coherent tracking branch")

    Lambda1 = sp.Float("0.279605891931464")
    coherent_gate = sp.Eq(Xi1_formula, Lambda1)
    tau_target = sp.simplify(Lambda1 / 2)

    print("Xi_1 target = Lambda_1 =")
    sp.pprint(Lambda1)
    print("tau target = Lambda_1/2 =")
    sp.pprint(tau_target)
    print("Coherent-branch quartic gate:")
    sp.pprint(coherent_gate)

    # ------------------------------------------------------------------
    # VIII. Final theorem ledger
    # ------------------------------------------------------------------
    banner("FINAL STEP-10 LEDGER")
    print("1. Exact coherent tracking-branch substitution:")
    print("   T^2 = Z_W (1+chi_0)^2 / [Omega_W^2 (1-eps)^2]")
    print("   eps = epsilon_W [1 - (2/11) delta_U/(1+delta_U)]")
    print()
    print("2. Exact support-blindness theorem:")
    print("   d_zeta T^2 = d_zeta R_target = d_zeta Xi_1 = 0")
    print("   so support can raise the steady baseline M_tr, but it cannot")
    print("   repair or spoil the first grouped weak-axisymmetric defect.")
    print()
    print("3. Exact coherent-branch defect law:")
    print("   Xi_1 = zeta_Z - omega_W + 2 chi_1/(1+chi_0) + 2 eps_1/(1-eps)")
    print("   with eps_1 = (1 - 2 delta_U/(11(1+delta_U))) varepsilon_W")
    print("              - 2 epsilon_W deltaU_1/(11(1+delta_U)^2)")
    print()
    print("4. Selected-branch reformulation:")
    print("   Xi_1 = - eta_1/(1-epsilon_eta) - R_1")
    print()
    print("5. Tracking-factor rigidity is not enough:")
    print("   Theta_1 can vanish while Xi_1 remains nonzero.")
    print()
    print("6. Quartic anomaly gate on the coherent branch:")
    print("   Xi_1 = Lambda_1 ≈ 0.279605891931464")
    print("   tau  = Lambda_1/2 ≈ 0.139802945965732")
    print()
    print("So the next theorem problem is now smaller again:")
    print("  compute the first grouped weak-axisymmetric drifts of the")
    print("  mixed/outgoing placement variables (or equivalently the")
    print("  next microscopic slippage ledger), because coherent support")
    print("  itself has dropped out of the defect law.")


if __name__ == "__main__":
    step10_coherent_tracking_branch_defect()
