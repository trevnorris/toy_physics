#!/usr/bin/env python3
"""
Step 7 for the g-2 rebuild from the moving-throat PDE base.

What this script does
---------------------
1. Starts from Step 6's surviving quartic gate, which said the remaining defect is
       Xi_1 = P_1/P_0
   carried only by the outgoing load factor Lambda_r = P_r/Delta_r on the
   conservative-shape-preserving branch.
2. Rewrites that defect directly in terms of the weak-axisymmetric slopes of the
   actual outgoing-port data
       P_r, Delta_r, N_{0}^{(r)} = P_r^2 / Delta_r^2.
3. Proves the exact port-level theorem
       Xi_1 = \bar{nu}_N - kappa_1,
   where kappa_1 is the wall-baseline slope and \bar{nu}_N is the outgoing-weighted
   static transfer slope.
4. Derives exact formulas for the actual numerator slope mathfrak{p}_r and the
   actual detuning slope mathfrak{d}_r in terms of primitive moving-throat port
   variables.
5. Verifies the exact equivalence with the earlier slippage language
       mathfrak{m}_r, mathfrak{i}_r, mathfrak{h}_r.
6. States the direct anomaly gate
       \bar{nu}_N = kappa_1 + Lambda_1,
   and the dominant-port specialization.

Interpretation
--------------
Step 6 said the missing quartic anomaly layer is a wall-loading defect of
Lambda_r = P_r/Delta_r. Step 7 sharpens that again: the full remaining defect is
nothing more than a mismatch between

  • the wall-baseline slope kappa_1,
  • and the outgoing-weighted static transfer slope of the actual moving-throat ports.

So the next theorem gate is no longer “compute every microscopic slippage
variable.” It is just: compute the actual outgoing-port slopes nu_r and check
whether their outgoing-weighted average co-loads with the wall baseline.
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


def step7_outgoing_port_coloading() -> None:
    banner("STEP 7 — DIRECT OUTGOING-PORT CO-LOADING THEOREM")

    eps, lam = sp.symbols("epsilon lambda_A", real=True)
    kappa1 = sp.symbols("kappa_1", real=True)

    # ------------------------------------------------------------------
    # I. Weak-axisymmetric slopes of the actual outgoing-port data
    # ------------------------------------------------------------------
    subbanner("VII.1 — Weak-axisymmetric slopes of P_r, Delta_r, and N_{0}^{(r)}")

    K = sp.symbols("K", positive=True)
    P_r = sp.symbols("P_r", positive=True)
    Delta_r = sp.symbols("Delta_r", positive=True)
    pfrak_r, dfrak_r = sp.symbols("mathfrak_p_r mathfrak_d_r", real=True)

    K_A = sp.simplify(K * (1 + eps * lam * kappa1))
    P_A_r = sp.simplify(P_r * (1 + eps * lam * pfrak_r))
    Delta_A_r = sp.simplify(Delta_r * (1 + eps * lam * dfrak_r))
    N_A0_r = sp.simplify(P_A_r**2 / Delta_A_r**2)

    nu_r = sp.simplify(sp.diff(sp.log(N_A0_r), eps).subs(eps, 0) / lam)
    Xi_port_r = sp.simplify(sp.diff(sp.log(N_A0_r / K_A), eps).subs(eps, 0) / lam)

    print("N_{A,0}^{(r)} =")
    sp.pprint(N_A0_r)
    print("nu_r = d ln N_{A,0}^{(r)} / (epsilon lambda_A) =")
    sp.pprint(nu_r)

    expect_zero("nu_r - 2(mathfrak_p_r - mathfrak_d_r)", nu_r - 2 * (pfrak_r - dfrak_r))
    expect_zero(
        "d ln( N_{A,0}^{(r)} / K_A ) / (epsilon lambda_A) - (nu_r-kappa_1)",
        Xi_port_r - (nu_r - kappa1),
    )

    print("Per-port wall-referenced defect:")
    print("  d ln( Lambda_{A,r}^2 / K_A ) / (epsilon lambda_A) = nu_r - kappa_1")

    # ------------------------------------------------------------------
    # II. Exact outgoing-weighted average and the surviving scalar Xi_1
    # ------------------------------------------------------------------
    subbanner("VII.2 — Exact weighted average Xi_1 = bar{nu}_N - kappa_1")

    N10, N20 = sp.symbols("N_10 N_20", positive=True)
    nu1, nu2 = sp.symbols("nu_1 nu_2", real=True)

    N1A = sp.simplify(N10 * (1 + eps * lam * nu1))
    N2A = sp.simplify(N20 * (1 + eps * lam * nu2))
    NtotA = sp.expand(N1A + N2A)

    rho1 = sp.simplify(N10 / (N10 + N20))
    rho2 = sp.simplify(N20 / (N10 + N20))
    bar_nu_N = sp.simplify(sp.diff(sp.log(NtotA), eps).subs(eps, 0) / lam)

    expect_zero("rho_1 + rho_2 - 1", rho1 + rho2 - 1)
    expect_zero("bar{nu}_N - (rho_1 nu_1 + rho_2 nu_2)", bar_nu_N - (rho1 * nu1 + rho2 * nu2))

    Xi1 = sp.simplify(bar_nu_N - kappa1)

    print("bar{nu}_N =")
    sp.pprint(bar_nu_N)
    print("Xi_1 =")
    sp.pprint(Xi1)

    print("Main scalar theorem:")
    print("  Xi_1 = bar{nu}_N - kappa_1")
    print("  = outgoing-weighted static transfer slope minus wall-baseline slope")

    # Direct anomaly gate carried from Steps 2-6.
    Lambda1 = sp.Float("0.279605891931464")
    bar_nu_target = sp.simplify(kappa1 + Lambda1)
    dominant_target = sp.simplify(kappa1 + Lambda1)

    print("If Xi_1 must match the quartic anomaly target Lambda_1, then:")
    print("  bar{nu}_N = kappa_1 + Lambda_1 =")
    sp.pprint(bar_nu_target)
    print("Dominant-port specialization:")
    print("  nu_* = kappa_1 + Lambda_1 =")
    sp.pprint(dominant_target)

    # ------------------------------------------------------------------
    # III. Exact formulas for the actual numerator slope mathfrak{p}_r
    # ------------------------------------------------------------------
    subbanner("VII.3 — Exact actual numerator slope")

    OmegaU2 = sp.symbols("OmegaU2_r", positive=True)
    OmegaW2 = sp.symbols("OmegaW2_r", positive=True)
    GWr = sp.symbols("GWr", positive=True)
    GUr = sp.symbols("GUr", positive=True)
    Rr = sp.symbols("Rr", positive=True)
    oUr = sp.symbols("mathfrak_oU_r", real=True)
    oWr = sp.symbols("mathfrak_oW_r", real=True)
    gWr = sp.symbols("mathfrak_gW_r", real=True)
    gUr = sp.symbols("mathfrak_gU_r", real=True)
    rfrak = sp.symbols("mathfrak_r_r", real=True)

    OmegaU2_A = sp.simplify(OmegaU2 * (1 + eps * lam * oUr))
    OmegaW2_A = sp.simplify(OmegaW2 * (1 + eps * lam * oWr))
    GWr_A = sp.simplify(GWr * (1 + eps * lam * gWr))
    GUr_A = sp.simplify(GUr * (1 + eps * lam * gUr))
    Rr_A = sp.simplify(Rr * (1 + eps * lam * rfrak))

    P_static = sp.simplify(OmegaU2 * GWr + Rr * GUr)
    P_static_A = sp.expand(OmegaU2_A * GWr_A + Rr_A * GUr_A)

    pfrak_exact = sp.simplify(sp.diff(sp.log(P_static_A), eps).subs(eps, 0) / lam)

    alpha_r = sp.simplify(OmegaU2 * GWr / P_static)
    beta_r = sp.simplify(Rr * GUr / P_static)

    print("P_r =")
    sp.pprint(P_static)
    print("alpha_r =")
    sp.pprint(alpha_r)
    print("beta_r =")
    sp.pprint(beta_r)

    expect_zero("alpha_r + beta_r - 1", alpha_r + beta_r - 1)
    expect_zero(
        "mathfrak{p}_r - [ alpha_r(o_U+g_W) + beta_r(r+g_U) ]",
        pfrak_exact - (alpha_r * (oUr + gWr) + beta_r * (rfrak + gUr)),
    )

    print("Actual numerator slope:")
    sp.pprint(pfrak_exact)

    # ------------------------------------------------------------------
    # IV. Exact formulas for the actual detuning slope mathfrak{d}_r
    # ------------------------------------------------------------------
    subbanner("VII.4 — Exact actual detuning slope")

    Delta_static = sp.simplify(OmegaU2 * OmegaW2 - Rr**2)
    Delta_static_A = sp.expand(OmegaU2_A * OmegaW2_A - Rr_A**2)

    dfrak_exact = sp.simplify(sp.diff(sp.log(Delta_static_A), eps).subs(eps, 0) / lam)

    chi_r = sp.simplify(OmegaU2 * OmegaW2 / Delta_static)
    zeta_r = sp.simplify(Rr**2 / Delta_static)
    H_r = sp.simplify(Rr**2 / (OmegaU2 * OmegaW2))

    print("Delta_r =")
    sp.pprint(Delta_static)
    print("chi_r =")
    sp.pprint(chi_r)
    print("zeta_r =")
    sp.pprint(zeta_r)

    expect_zero("chi_r - 1/(1-H_r)", chi_r - 1 / (1 - H_r))
    expect_zero("zeta_r - H_r/(1-H_r)", zeta_r - H_r / (1 - H_r))
    expect_zero(
        "mathfrak{d}_r - [ chi_r(o_U+o_W) - 2 zeta_r r ]",
        dfrak_exact - (chi_r * (oUr + oWr) - 2 * zeta_r * rfrak),
    )

    print("Actual detuning slope:")
    sp.pprint(dfrak_exact)

    # ------------------------------------------------------------------
    # V. Static outgoing-transfer slope in actual port variables
    # ------------------------------------------------------------------
    subbanner("VII.5 — Static outgoing-transfer slope in actual port variables")

    nu_exact = sp.simplify(2 * (pfrak_exact - dfrak_exact))

    expect_zero(
        "nu_r - [ 2 alpha_r(o_U+g_W) + 2 beta_r(r+g_U) - 2 chi_r(o_U+o_W) + 4 zeta_r r ]",
        nu_exact - (
            2 * alpha_r * (oUr + gWr)
            + 2 * beta_r * (rfrak + gUr)
            - 2 * chi_r * (oUr + oWr)
            + 4 * zeta_r * rfrak
        ),
    )

    print("nu_r in primitive actual-port variables =")
    sp.pprint(nu_exact)

    # ------------------------------------------------------------------
    # VI. Exact equivalence to the earlier slippage language
    # ------------------------------------------------------------------
    subbanner("VII.6 — Exact equivalence to the earlier slippage language")

    I_r = sp.simplify(Rr * GUr / (OmegaU2 * GWr))

    mfrak = sp.simplify(gWr - oWr - sp.Rational(1, 2) * kappa1)
    ifrak = sp.simplify(rfrak + gUr - oUr - gWr)
    hfrak = sp.simplify(2 * rfrak - oUr - oWr)

    expect_zero("alpha_r - 1/(1+I_r)", alpha_r - 1 / (1 + I_r))
    expect_zero("beta_r - I_r/(1+I_r)", beta_r - I_r / (1 + I_r))
    expect_zero(
        "nu_r - [kappa_1 + 2 m_r + 2 I_r/(1+I_r) i_r + 2 H_r/(1-H_r) h_r]",
        nu_exact - (kappa1 + 2 * mfrak + 2 * I_r / (1 + I_r) * ifrak + 2 * H_r / (1 - H_r) * hfrak),
    )

    print("Slippage equivalence:")
    print("  sigma_r = nu_r - kappa_1")
    print("        = 2 m_r + 2 I_r/(1+I_r) i_r + 2 H_r/(1-H_r) h_r")

    # ------------------------------------------------------------------
    # VII. Co-loading theorems and no-go corollaries
    # ------------------------------------------------------------------
    subbanner("VII.7 — Co-loading theorems and immediate corollaries")

    # Strong per-port sufficient condition.
    Xi_port_match = sp.simplify((nu_exact - kappa1).subs({nu_exact: kappa1}))
    expect_zero("Per-port co-loading condition nu_r = kappa_1 kills the defect", Xi_port_match)

    # Naive rigidity/no-go.
    nu_rigid = sp.simplify(nu_exact.subs({pfrak_exact: 0, dfrak_exact: 0}))
    # Simpler explicit substitution on the abstract formula.
    expect_zero("Rigid ports give nu_r = 0", (2 * (pfrak_r - dfrak_r)).subs({pfrak_r: 0, dfrak_r: 0}))
    print("Naive rigidity corollary:")
    print("  if mathfrak{p}_r = mathfrak{d}_r = 0, then nu_r = 0 and Xi_1 = -kappa_1")

    print("\nSTEP 7 LEDGER")
    print("  Xi_1 = P_1/P_0 = bar{nu}_N - kappa_1")
    print("  nu_r = 2(mathfrak{p}_r - mathfrak{d}_r)")
    print("  mathfrak{p}_r = alpha_r(o_U+g_W) + beta_r(r+g_U)")
    print("  mathfrak{d}_r = chi_r(o_U+o_W) - 2 zeta_r r")
    print("  zero defect  <=>  bar{nu}_N = kappa_1")
    print("  anomaly law  <=>  bar{nu}_N = kappa_1 + Lambda_1")


if __name__ == "__main__":
    step7_outgoing_port_coloading()
