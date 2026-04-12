#!/usr/bin/env python3
"""
Step 8 for the g-2 rebuild from the moving-throat PDE base.

What this script does
---------------------
1. Starts from Step 7's exact actual-port theorem
       Xi_1 = bar{nu}_N - kappa_1,
   where nu_r is the weak-axisymmetric logarithmic slope of the actual static
   outgoing-transfer coefficient N_{0}^{(r)} = P_r^2 / Delta_r^2.
2. Rewrites the actual outgoing-transfer coefficient in wall-normalized form,
   proving the exact factorization
       N_{0}^{(r)} = K * T_r^2,
   with one dimensionless transfer shape T_r.
3. Derives the weak-axisymmetric transport laws for the wall-normalized port
   variables and proves
       nu_r = kappa_1 + 2 tau_r,
   where tau_r is the transfer-shape slope.
4. Collapses the remaining grouped defect to
       Xi_1 = 2 sum_r rho_r^(N) tau_r.
5. Verifies exact equivalence with the earlier slippage language and shows that
   the anomaly-matching condition becomes
       sum_r rho_r^(N) tau_r = Lambda_1 / 2.

Interpretation
--------------
Step 7 said the remaining quartic anomaly layer is the mismatch between the
outgoing-weighted static transfer slope and the wall-baseline slope. Step 8
sharpens that again: each actual outgoing port factors into the wall baseline
K times a dimensionless transfer shape squared. So the whole remaining theorem
problem collapses to the weak-axisymmetric rigidity of one wall-normalized
transfer shape per port.
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


def step8_wall_normalized_transfer_shape() -> None:
    banner("STEP 8 — WALL-NORMALIZED TRANSFER-SHAPE THEOREM")

    # ------------------------------------------------------------------
    # I. Exact wall-normalized factorization of N_0^{(r)}
    # ------------------------------------------------------------------
    subbanner("VIII.1 — Exact factorization N_0^{(r)} = K * T_r^2")

    K = sp.symbols("K", positive=True, real=True)
    OmegaU, OmegaW = sp.symbols("OmegaU OmegaW", positive=True, real=True)
    GWr, GUr, Rr = sp.symbols("GWr GUr Rr", positive=True, real=True)

    P_r = sp.simplify(OmegaU**2 * GWr + Rr * GUr)
    Delta_r = sp.simplify(OmegaU**2 * OmegaW**2 - Rr**2)
    N0_r = sp.simplify(P_r**2 / Delta_r**2)

    GhatW = sp.simplify(GWr / (OmegaW**2 * sp.sqrt(K)))
    GhatU = sp.simplify(GUr / (OmegaU * OmegaW * sp.sqrt(K)))
    Rhat = sp.simplify(Rr / (OmegaU * OmegaW))
    T_r = sp.simplify((GhatW + Rhat * GhatU) / (1 - Rhat**2))

    print("P_r =")
    sp.pprint(P_r)
    print("Delta_r =")
    sp.pprint(Delta_r)
    print("T_r =")
    sp.pprint(T_r)

    expect_zero(
        "P_r - sqrt(K) * OmegaU^2 * OmegaW^2 * (GhatW + Rhat*GhatU)",
        P_r - sp.sqrt(K) * OmegaU**2 * OmegaW**2 * (GhatW + Rhat * GhatU),
    )
    expect_zero(
        "Delta_r - OmegaU^2 * OmegaW^2 * (1 - Rhat^2)",
        Delta_r - OmegaU**2 * OmegaW**2 * (1 - Rhat**2),
    )
    expect_zero(
        "N_0^{(r)} - K * T_r^2",
        N0_r - K * T_r**2,
    )

    # ------------------------------------------------------------------
    # II. Weak-axisymmetric transport of wall-normalized variables
    # ------------------------------------------------------------------
    subbanner("VIII.2 — Weak-axisymmetric slopes of the wall-normalized port variables")

    eps, lam = sp.symbols("epsilon lambda_A", real=True)
    kappa1 = sp.symbols("kappa_1", real=True)
    oUr, oWr = sp.symbols("mathfrak_oU_r mathfrak_oW_r", real=True)
    gWr, gUr = sp.symbols("mathfrak_gW_r mathfrak_gU_r", real=True)
    rfrak = sp.symbols("mathfrak_r_r", real=True)

    # Primitive weak-axisymmetric branch. We only keep first order in eps.
    sqrtK_A = trunc1(sp.sqrt(K) * (1 + eps * lam * kappa1 / 2), eps)
    OmegaU_A = trunc1(OmegaU * (1 + eps * lam * oUr / 2), eps)
    OmegaW_A = trunc1(OmegaW * (1 + eps * lam * oWr / 2), eps)
    GWr_A = trunc1(GWr * (1 + eps * lam * gWr), eps)
    GUr_A = trunc1(GUr * (1 + eps * lam * gUr), eps)
    Rr_A = trunc1(Rr * (1 + eps * lam * rfrak), eps)

    GhatW_A = trunc1(GWr_A / (OmegaW_A**2 * sqrtK_A), eps)
    GhatU_A = trunc1(GUr_A / (OmegaU_A * OmegaW_A * sqrtK_A), eps)
    Rhat_A = trunc1(Rr_A / (OmegaU_A * OmegaW_A), eps)

    wfrak = sp.simplify(sp.diff(sp.log(GhatW_A), eps).subs(eps, 0) / lam)
    ufrak = sp.simplify(sp.diff(sp.log(GhatU_A), eps).subs(eps, 0) / lam)
    cfrak = sp.simplify(sp.diff(sp.log(Rhat_A), eps).subs(eps, 0) / lam)

    expect_zero(
        "mathfrak_w_r - (g_W - o_W - kappa_1/2)",
        wfrak - (gWr - oWr - sp.Rational(1, 2) * kappa1),
    )
    expect_zero(
        "mathfrak_u_r - (g_U - o_U/2 - o_W/2 - kappa_1/2)",
        ufrak - (gUr - sp.Rational(1, 2) * oUr - sp.Rational(1, 2) * oWr - sp.Rational(1, 2) * kappa1),
    )
    expect_zero(
        "mathfrak_c_r - (r - o_U/2 - o_W/2)",
        cfrak - (rfrak - sp.Rational(1, 2) * oUr - sp.Rational(1, 2) * oWr),
    )

    print("mathfrak_w_r =")
    sp.pprint(wfrak)
    print("mathfrak_u_r =")
    sp.pprint(ufrak)
    print("mathfrak_c_r =")
    sp.pprint(cfrak)

    # ------------------------------------------------------------------
    # III. Exact transfer-shape slope and nu_r = kappa_1 + 2 tau_r
    # ------------------------------------------------------------------
    subbanner("VIII.3 — Exact transfer-shape slope and the identity nu_r = kappa_1 + 2 tau_r")

    w, u, c = sp.symbols("mathfrak_w_r mathfrak_u_r mathfrak_c_r", real=True)
    GhatW_A2 = trunc1(GhatW * (1 + eps * lam * w), eps)
    GhatU_A2 = trunc1(GhatU * (1 + eps * lam * u), eps)
    Rhat_A2 = trunc1(Rhat * (1 + eps * lam * c), eps)
    T_A = trunc1((GhatW_A2 + Rhat_A2 * GhatU_A2) / (1 - Rhat_A2**2), eps)
    K_A = trunc1(K * (1 + eps * lam * kappa1), eps)
    N_A0 = trunc1(K_A * T_A**2, eps)

    tau_r = sp.simplify(sp.diff(sp.log(T_A), eps).subs(eps, 0) / lam)
    nu_r = sp.simplify(sp.diff(sp.log(N_A0), eps).subs(eps, 0) / lam)

    ahat = sp.simplify(GhatW / (GhatW + Rhat * GhatU))
    bhat = sp.simplify(Rhat * GhatU / (GhatW + Rhat * GhatU))

    expect_zero("ahat + bhat - 1", ahat + bhat - 1)
    expect_zero(
        "tau_r - [ ahat*w + bhat*(u+c) + 2 Rhat^2/(1-Rhat^2) * c ]",
        tau_r - (ahat * w + bhat * (u + c) + 2 * Rhat**2 / (1 - Rhat**2) * c),
    )
    expect_zero("nu_r - (kappa_1 + 2*tau_r)", nu_r - (kappa1 + 2 * tau_r))

    print("tau_r =")
    sp.pprint(tau_r)
    print("nu_r =")
    sp.pprint(nu_r)

    # ------------------------------------------------------------------
    # IV. Collapse of the remaining grouped defect
    # ------------------------------------------------------------------
    subbanner("VIII.4 — Exact collapse Xi_1 = 2 sum_r rho_r^(N) tau_r")

    N10, N20 = sp.symbols("N_10 N_20", positive=True, real=True)
    tau1, tau2 = sp.symbols("tau_1 tau_2", real=True)
    nu1, nu2 = sp.symbols("nu_1 nu_2", real=True)

    rho1 = sp.simplify(N10 / (N10 + N20))
    rho2 = sp.simplify(N20 / (N10 + N20))
    Xi1 = sp.simplify(rho1 * (nu1 - kappa1) + rho2 * (nu2 - kappa1))
    Xi1_tau = sp.simplify(Xi1.subs({nu1: kappa1 + 2 * tau1, nu2: kappa1 + 2 * tau2}))

    expect_zero("rho_1 + rho_2 - 1", rho1 + rho2 - 1)
    expect_zero(
        "Xi_1 - 2(rho_1 tau_1 + rho_2 tau_2)",
        Xi1_tau - 2 * (rho1 * tau1 + rho2 * tau2),
    )

    print("Xi_1 =")
    sp.pprint(Xi1_tau)

    # ------------------------------------------------------------------
    # V. Exact equivalence to the earlier slippage language
    # ------------------------------------------------------------------
    subbanner("VIII.5 — Exact equivalence to the earlier slippage language")

    Mcal = sp.simplify(GhatW)
    Ical = sp.simplify(Rhat * GhatU / GhatW)
    Hcal = sp.simplify(Rhat**2)

    mfrak, ifrak, hfrak = sp.symbols("mathfrak_m_r mathfrak_i_r mathfrak_h_r", real=True)
    tau_slip = sp.simplify(
        tau_r.subs({w: mfrak, u: ifrak + mfrak - hfrak / 2, c: hfrak / 2})
    )

    expect_zero("Ical - (Rhat*GhatU/GhatW)", Ical - (Rhat * GhatU / GhatW))
    expect_zero("Hcal - Rhat^2", Hcal - Rhat**2)
    expect_zero(
        "tau_r - [ m + I/(1+I) i + H/(1-H) h ]",
        tau_slip - (mfrak + Ical / (1 + Ical) * ifrak + Hcal / (1 - Hcal) * hfrak),
    )

    sigma_r = sp.simplify(2 * tau_r)
    print("sigma_r = 2 tau_r =")
    sp.pprint(sigma_r)

    # ------------------------------------------------------------------
    # VI. Direct anomaly gate in transfer-shape language
    # ------------------------------------------------------------------
    subbanner("VIII.6 — Direct anomaly gate in transfer-shape language")

    Lambda1 = sp.Float("0.279605891931464")
    tau_target = sp.simplify(Lambda1 / 2)
    print("If Xi_1 = Lambda_1, then the transfer-shape target is")
    print("  sum_r rho_r^(N) tau_r = Lambda_1 / 2 =")
    sp.pprint(tau_target)
    print("Dominant-port specialization:")
    print("  tau_* = Lambda_1 / 2 =")
    sp.pprint(tau_target)

    # ------------------------------------------------------------------
    # VII. Reduced theorem ledger
    # ------------------------------------------------------------------
    banner("FINAL STEP-8 LEDGER")
    print("1. Exact wall-normalized factorization:")
    print("   N_0^{(r)} = K * T_r^2")
    print()
    print("2. Exact transfer-shape slope law:")
    print("   nu_r = kappa_1 + 2 tau_r")
    print()
    print("3. Exact collapse of the remaining grouped defect:")
    print("   Xi_1 = 2 sum_r rho_r^(N) tau_r")
    print()
    print("4. Therefore anomaly matching requires:")
    print(f"   sum_r rho_r^(N) tau_r = Lambda_1/2 = {sp.N(tau_target, 16)}")
    print()
    print("5. So the next honest theorem gate is no longer to compute the raw port slopes")
    print("   nu_r. It is to compute the weak-axisymmetric wall-normalized transfer-shape")
    print("   slopes tau_r on the actual moving-throat branch.")


if __name__ == "__main__":
    step8_wall_normalized_transfer_shape()
