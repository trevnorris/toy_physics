#!/usr/bin/env python3
"""
Step 6 for the g-2 rebuild from the moving-throat PDE base.

What this script does
---------------------
1. Starts from Step 5's exact loading scalar
       Xi_load = N_01/N_0 - D_01/D_0,
   and rewrites the conservative static slope D_01/D_0 in terms of the wall
   baseline, BdG support dressing, and conservative Maxwell/mixed dressing.
2. Verifies that Xi_load is exactly a weighted failure of static
   self-similarity relative to the wall baseline slope delta_K.
3. Derives the microscopic weighted-log formulas for the three bundles:
       BdG support,
       conservative Maxwell/mixed static loading,
       outgoing-transfer loading.
4. Factors the microscopic bundles by the wall baseline into the natural
   wall-normalized variables
       chi_alpha, Upsilon_r, Lambda_r,
   and proves the exact compact formula
       Xi_load = Theta_N + omega_B Theta_B + omega_Z Theta_Z.
5. Specializes to conservative-shape-preserving branches and proves the
   outgoing-load theorem
       Xi_load = sum_r rho_r^(N) (2 d ln Lambda_r - delta_K).
6. Derives two immediate consequences:
       (a) naive common self-similarity gives Xi_load = -delta_K,
       (b) on a uniform outgoing-load branch, matching a nonzero anomaly target
           Lambda_1 requires
               d ln Lambda = (delta_K + Lambda_1)/2,
           i.e. an extra outgoing slippage of Lambda_1/2 above the
           defect-killing self-similar tracking law.

Interpretation
--------------
Step 5 collapsed the quartic g-2 correction to a one-scalar loading mismatch.
Step 6 sharpens that further: on conservative-shape-preserving branches, the
entire remaining defect is not a generic static mismatch of the whole throat
bundle. It is purely an outgoing-load theorem.
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


def step6_static_self_similarity_outgoing_load() -> None:
    banner("STEP 6 — STATIC SELF-SIMILARITY AND THE OUTGOING-LOAD THEOREM")

    # ------------------------------------------------------------------
    # I. Exact static-slope decomposition of the loading scalar
    # ------------------------------------------------------------------
    subbanner("VI.1 — Exact static-slope decomposition of Xi_load")

    K, B0, Z0, N0 = sp.symbols("K B_0 Z_0 N_0", positive=True)
    K1, B1, Z1, N1 = sp.symbols("K_1 B_1 Z_1 N_1", real=True)

    D0 = sp.simplify(K - B0 - Z0)
    D1 = sp.simplify(K1 - B1 - Z1)

    delta_K = sp.simplify(K1 / K)
    delta_B = sp.simplify(B1 / B0)
    delta_Z = sp.simplify(Z1 / Z0)
    delta_N = sp.simplify(N1 / N0)
    delta_D = sp.simplify(D1 / D0)

    omega_K = sp.simplify(K / D0)
    omega_B = sp.simplify(B0 / D0)
    omega_Z = sp.simplify(Z0 / D0)

    Xi_load = sp.simplify(delta_N - delta_D)

    print("D_0 =")
    sp.pprint(D0)
    print("delta_D = D_01 / D_0 =")
    sp.pprint(delta_D)
    print("Xi_load = N_01/N_0 - D_01/D_0 =")
    sp.pprint(Xi_load)

    expect_zero(
        "delta_D - (omega_K*delta_K - omega_B*delta_B - omega_Z*delta_Z)",
        delta_D - (omega_K * delta_K - omega_B * delta_B - omega_Z * delta_Z),
    )
    expect_zero(
        "omega_K - omega_B - omega_Z - 1",
        omega_K - omega_B - omega_Z - 1,
    )
    expect_zero(
        "Xi_load - ((delta_N-delta_K)+omega_B*(delta_B-delta_K)+omega_Z*(delta_Z-delta_K))",
        Xi_load - ((delta_N - delta_K) + omega_B * (delta_B - delta_K) + omega_Z * (delta_Z - delta_K)),
    )

    print("Wall-referenced decomposition:")
    print("  Xi_load = (delta_N-delta_K) + omega_B(delta_B-delta_K) + omega_Z(delta_Z-delta_K)")

    # ------------------------------------------------------------------
    # II. Microscopic weighted-log formulas for the three bundles
    # ------------------------------------------------------------------
    subbanner("VI.2 — Microscopic weighted-log formulas")

    eps = sp.symbols("epsilon", real=True)

    # BdG support bundle: two symbolic modes are enough to expose the exact pattern.
    c1, c2 = sp.symbols("c_1 c_2", positive=True)
    w1, w2 = sp.symbols("varpi_1 varpi_2", positive=True)
    sc1, sc2 = sp.symbols("sigma_c1 sigma_c2", real=True)
    sw1, sw2 = sp.symbols("sigma_w1 sigma_w2", real=True)

    c1A = c1 * (1 + eps * sc1)
    c2A = c2 * (1 + eps * sc2)
    w1A = w1 * (1 + eps * sw1)
    w2A = w2 * (1 + eps * sw2)

    B01_mode = sp.simplify(c1**2 / w1**2)
    B02_mode = sp.simplify(c2**2 / w2**2)
    B0_modes = sp.simplify(B01_mode + B02_mode)
    B0A = sp.simplify(c1A**2 / w1A**2 + c2A**2 / w2A**2)
    delta_B_micro = sp.simplify(sp.diff(B0A, eps).subs(eps, 0) / B0_modes)
    rho_B1 = sp.simplify(B01_mode / B0_modes)
    rho_B2 = sp.simplify(B02_mode / B0_modes)

    print("delta_B (BdG support bundle) =")
    sp.pprint(delta_B_micro)
    expect_zero(
        "delta_B - sum rho_B * 2 dln(c/varpi)",
        delta_B_micro - (rho_B1 * (2 * sc1 - 2 * sw1) + rho_B2 * (2 * sc2 - 2 * sw2)),
    )

    # Conservative Maxwell/mixed static bundle.
    Q1, Q2 = sp.symbols("Q_1 Q_2", positive=True)
    Del1, Del2 = sp.symbols("Delta_1 Delta_2", positive=True)
    q1s, q2s = sp.symbols("sigma_Q1 sigma_Q2", real=True)
    d1s, d2s = sp.symbols("sigma_Delta1 sigma_Delta2", real=True)

    Q1A = Q1 * (1 + eps * q1s)
    Q2A = Q2 * (1 + eps * q2s)
    Del1A = Del1 * (1 + eps * d1s)
    Del2A = Del2 * (1 + eps * d2s)

    Z01_port = sp.simplify(Q1 / Del1)
    Z02_port = sp.simplify(Q2 / Del2)
    Z0_ports = sp.simplify(Z01_port + Z02_port)
    Z0A = sp.simplify(Q1A / Del1A + Q2A / Del2A)
    delta_Z_micro = sp.simplify(sp.diff(Z0A, eps).subs(eps, 0) / Z0_ports)
    rho_Z1 = sp.simplify(Z01_port / Z0_ports)
    rho_Z2 = sp.simplify(Z02_port / Z0_ports)

    print("delta_Z (conservative Maxwell/mixed bundle) =")
    sp.pprint(delta_Z_micro)
    expect_zero(
        "delta_Z - sum rho_Z * dln(Q/Delta)",
        delta_Z_micro - (rho_Z1 * (q1s - d1s) + rho_Z2 * (q2s - d2s)),
    )

    # Outgoing-transfer bundle.
    P1r, P2r = sp.symbols("P_1 P_2", positive=True)
    p1s, p2s = sp.symbols("sigma_P1 sigma_P2", real=True)

    P1A = P1r * (1 + eps * p1s)
    P2A = P2r * (1 + eps * p2s)

    N01_port = sp.simplify(P1r**2 / Del1**2)
    N02_port = sp.simplify(P2r**2 / Del2**2)
    N0_ports = sp.simplify(N01_port + N02_port)
    N0A = sp.simplify(P1A**2 / Del1A**2 + P2A**2 / Del2A**2)
    delta_N_micro = sp.simplify(sp.diff(N0A, eps).subs(eps, 0) / N0_ports)
    rho_N1 = sp.simplify(N01_port / N0_ports)
    rho_N2 = sp.simplify(N02_port / N0_ports)

    print("delta_N (outgoing-transfer bundle) =")
    sp.pprint(delta_N_micro)
    expect_zero(
        "delta_N - sum rho_N * 2 dln(P/Delta)",
        delta_N_micro - (rho_N1 * (2 * p1s - 2 * d1s) + rho_N2 * (2 * p2s - 2 * d2s)),
    )

    # ------------------------------------------------------------------
    # III. Exact wall-referenced self-similarity defect fields
    # ------------------------------------------------------------------
    subbanner("VI.3 — Exact wall-referenced self-similarity defect fields")

    omB, omZ = sp.symbols("omega_B omega_Z", real=True)

    Sigma_B1 = sp.simplify(2 * sc1 - 2 * sw1 - delta_K)
    Sigma_B2 = sp.simplify(2 * sc2 - 2 * sw2 - delta_K)
    Sigma_Z1 = sp.simplify(q1s - d1s - delta_K)
    Sigma_Z2 = sp.simplify(q2s - d2s - delta_K)
    Sigma_N1 = sp.simplify(2 * p1s - 2 * d1s - delta_K)
    Sigma_N2 = sp.simplify(2 * p2s - 2 * d2s - delta_K)

    Xi_micro = sp.simplify((delta_N_micro - delta_K) + omB * (delta_B_micro - delta_K) + omZ * (delta_Z_micro - delta_K))
    Xi_sigma = sp.simplify(
        rho_N1 * Sigma_N1 + rho_N2 * Sigma_N2
        + omB * (rho_B1 * Sigma_B1 + rho_B2 * Sigma_B2)
        + omZ * (rho_Z1 * Sigma_Z1 + rho_Z2 * Sigma_Z2)
    )

    print("Xi_load in wall-referenced microscopic form =")
    sp.pprint(Xi_sigma)
    expect_zero("Xi_micro - Xi_sigma", Xi_micro - Xi_sigma)

    # ------------------------------------------------------------------
    # IV. Factorization by the wall baseline into chi, Upsilon, Lambda
    # ------------------------------------------------------------------
    subbanner("VI.4 — Wall-normalized shape/load variables")

    KA = K * (1 + eps * delta_K)

    chi1 = sp.simplify(c1 / (sp.sqrt(K) * w1))
    chi2 = sp.simplify(c2 / (sp.sqrt(K) * w2))
    chi1A = sp.simplify(c1A / (sp.sqrt(KA) * w1A))
    chi2A = sp.simplify(c2A / (sp.sqrt(KA) * w2A))

    Upsilon1 = sp.simplify(Q1 / (K * Del1))
    Upsilon2 = sp.simplify(Q2 / (K * Del2))
    Upsilon1A = sp.simplify(Q1A / (KA * Del1A))
    Upsilon2A = sp.simplify(Q2A / (KA * Del2A))

    Lambda1 = sp.simplify(P1r / Del1)
    Lambda2 = sp.simplify(P2r / Del2)
    Lambda1A = sp.simplify(P1A / Del1A)
    Lambda2A = sp.simplify(P2A / Del2A)

    expect_zero("B_{0,1} - K chi_1^2", B01_mode - K * chi1**2)
    expect_zero("B_{0,2} - K chi_2^2", B02_mode - K * chi2**2)
    expect_zero("Z_{0,1} - K Upsilon_1", Z01_port - K * Upsilon1)
    expect_zero("Z_{0,2} - K Upsilon_2", Z02_port - K * Upsilon2)
    expect_zero("N_{0,1} - Lambda_1^2", N01_port - Lambda1**2)
    expect_zero("N_{0,2} - Lambda_2^2", N02_port - Lambda2**2)

    dln_chi1_sq = sp.simplify(sp.diff(sp.log(chi1A**2), eps).subs(eps, 0))
    dln_chi2_sq = sp.simplify(sp.diff(sp.log(chi2A**2), eps).subs(eps, 0))
    dln_Upsilon1 = sp.simplify(sp.diff(sp.log(Upsilon1A), eps).subs(eps, 0))
    dln_Upsilon2 = sp.simplify(sp.diff(sp.log(Upsilon2A), eps).subs(eps, 0))
    dln_Lambda1sq_over_K = sp.simplify(sp.diff(sp.log(Lambda1A**2 / KA), eps).subs(eps, 0))
    dln_Lambda2sq_over_K = sp.simplify(sp.diff(sp.log(Lambda2A**2 / KA), eps).subs(eps, 0))

    expect_zero("dln chi_1^2 - Sigma_B1", dln_chi1_sq - Sigma_B1)
    expect_zero("dln chi_2^2 - Sigma_B2", dln_chi2_sq - Sigma_B2)
    expect_zero("dln Upsilon_1 - Sigma_Z1", dln_Upsilon1 - Sigma_Z1)
    expect_zero("dln Upsilon_2 - Sigma_Z2", dln_Upsilon2 - Sigma_Z2)
    expect_zero("dln(Lambda_1^2/K) - Sigma_N1", dln_Lambda1sq_over_K - Sigma_N1)
    expect_zero("dln(Lambda_2^2/K) - Sigma_N2", dln_Lambda2sq_over_K - Sigma_N2)

    Theta_B = sp.simplify(rho_B1 * dln_chi1_sq + rho_B2 * dln_chi2_sq)
    Theta_Z = sp.simplify(rho_Z1 * dln_Upsilon1 + rho_Z2 * dln_Upsilon2)
    Theta_N = sp.simplify(rho_N1 * dln_Lambda1sq_over_K + rho_N2 * dln_Lambda2sq_over_K)

    expect_zero("Xi_micro - (Theta_N + omega_B Theta_B + omega_Z Theta_Z)", Xi_micro - (Theta_N + omB * Theta_B + omZ * Theta_Z))

    print("Compact wall-normalized formula:")
    print("  Xi_load = Theta_N + omega_B Theta_B + omega_Z Theta_Z")

    # ------------------------------------------------------------------
    # V. Conservative-shape theorem and outgoing-load theorem
    # ------------------------------------------------------------------
    subbanner("VI.5 — Conservative-shape theorem and outgoing-load theorem")

    ell1, ell2, ell, Lam1 = sp.symbols("ell_1 ell_2 ell Lambda_1", real=True)
    Xi_conservative_shape = sp.simplify(rho_N1 * (2 * ell1 - delta_K) + rho_N2 * (2 * ell2 - delta_K))
    Xi_conservative_shape_avg = sp.simplify(2 * (rho_N1 * ell1 + rho_N2 * ell2) - delta_K)
    expect_zero("Xi_conservative_shape - (2 weighted outgoing load - delta_K)", Xi_conservative_shape - Xi_conservative_shape_avg)

    Xi_uniform = sp.simplify(Xi_conservative_shape.subs({ell1: ell, ell2: ell}))
    print("On conservative-shape-preserving branches:")
    print("  Xi_load = sum_r rho_r^(N) (2 dln Lambda_r - delta_K)")
    print("Uniform outgoing-load law:")
    print("  Xi_load =")
    sp.pprint(Xi_uniform)

    naive_self_similarity = sp.simplify(Xi_uniform.subs(ell, 0))
    defect_kill_law = sp.simplify(sp.solve(sp.Eq(Xi_uniform, 0), ell)[0])
    anomaly_law = sp.simplify(sp.solve(sp.Eq(Xi_uniform, Lam1), ell)[0])
    extra_slip = sp.simplify(anomaly_law - defect_kill_law)

    print("Naive common self-similarity (dln Lambda = 0) gives:")
    sp.pprint(naive_self_similarity)
    print("Defect-killing wall-loading law (Xi_load = 0):")
    sp.pprint(defect_kill_law)
    print("Anomaly-matching uniform law (Xi_load = Lambda_1):")
    sp.pprint(anomaly_law)
    print("Extra outgoing slippage above the defect-killing law:")
    sp.pprint(extra_slip)

    expect_zero("naive self-similarity + delta_K", naive_self_similarity + delta_K)
    expect_zero("defect-killing law - delta_K/2", defect_kill_law - delta_K / 2)
    expect_zero("anomaly law - (delta_K + Lambda_1)/2", anomaly_law - (delta_K + Lam1) / 2)
    expect_zero("extra slippage - Lambda_1/2", extra_slip - Lam1 / 2)

    # ------------------------------------------------------------------
    # VI. Numerical anomaly benchmark
    # ------------------------------------------------------------------
    subbanner("VI.6 — Numerical benchmark")

    Lam1_num = sp.Float("0.279605891931464", 30)
    extra_num = sp.N(extra_slip.subs(Lam1, Lam1_num), 18)

    print(f"Lambda_1 = {Lam1_num}")
    print(f"Required extra outgoing-load slippage above wall-tracking = {extra_num}")

    # ------------------------------------------------------------------
    # VII. Interpretation
    # ------------------------------------------------------------------
    subbanner("VI.7 — Interpretation")
    print("1. Step 5's loading scalar Xi_load is exactly a wall-referenced")
    print("   self-similarity defect: it measures how the support/transfer bundles")
    print("   fail to co-load with the wall baseline.")
    print("2. After wall-normalization, the conservative bundles become pure shape")
    print("   variables (chi_alpha, Upsilon_r), while the outgoing bundle becomes")
    print("   a wall-loading law for Lambda_r = P_r/Delta_r.")
    print("3. On conservative-shape-preserving branches, the whole remaining linear")
    print("   grouped defect is carried only by the outgoing load factor.")
    print("4. Therefore the quartic g-2 target is not asking for an arbitrary new")
    print("   static bundle deformation. It is asking for a very specific outgoing")
    print("   slippage relative to the wall baseline.")


if __name__ == "__main__":
    step6_static_self_similarity_outgoing_load()
