#!/usr/bin/env python3
"""
same_charge_barrier_audit_stage003_dynamic_mixed_port_kernel_sympy_audit.py

Stage 003 — dynamic mixed-port kernel audit.

What this script does
---------------------
1. Builds the exact frequency-domain one-port wall / U / W bundle with a dynamic
   wall kernel K_B(omega) and an outgoing-port correction Pi(omega) on the mixed
   W lane.
2. Verifies the exact determinant identity

       det(K_dyn) = Delta_Pi * D_Pi,

   where
       Delta_Pi = A(omega) W(omega) - R^2,
       D_Pi     = K_B(omega) - Q_Pi(omega)/Delta_Pi.
3. Verifies the exact inverse-entry formulas generalizing the Stage-002 static
   susceptibilities to the dynamic bundle.
4. Proves the exact collinear-source susceptibility law

       V_mix = -1/2 chi_s(omega) S(omega)^2,

   with chi_s = N_s(omega)/(Delta_Pi D_Pi).
5. Verifies the dynamic product-family decomposition for the primitive reduced
   load families

       S_Q(x) = x^{-3},
       S_Y(x) = exp(-2 kappa x)/x.

   The linear dynamic bundle still produces only
       x^{-6},
       exp(-2 kappa x)/x^4,
       exp(-4 kappa x)/x^2,
   now with complex frequency-dependent coefficients.
6. Proves the exact outgoing-port derivative identity

       d(K_dyn^{-1})/dPi = K_dyn^{-1} E_W K_dyn^{-1},
       dV_mix/dPi        = -1/2 (e_W^T K_dyn^{-1} J)^2.

   So the leading passive/outgoing correction is controlled by the square of one
   exact transfer amplitude.
7. Verifies the conservative/static reduction of the dynamic bundle at
       omega = 0,
       Pi(0) = 0.
"""

from __future__ import annotations

import sympy as sp


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

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
        simplified = expr.applyfunc(lambda z: sp.factor(sp.simplify(z)))
        print(f"{name} =")
        sp.pprint(simplified)
        if any(entry != 0 for entry in simplified):
            raise AssertionError(f"{name} is not zero")
    else:
        simplified = sp.factor(sp.simplify(expr))
        print(f"{name} = {simplified}")
        if simplified != 0:
            raise AssertionError(f"{name} is not zero")


# ---------------------------------------------------------------------------
# Part I. Dynamic one-port bundle and exact determinant
# ---------------------------------------------------------------------------

def verify_dynamic_bundle():
    banner("PART I — DYNAMIC ONE-PORT MIXED BUNDLE")

    KB, GU, GW, A0, W0, R, Pi = sp.symbols(
        "K_B G_U G_W A_0 W_0 R Pi", real=True
    )

    A = A0
    W = sp.simplify(W0 - Pi)
    Delta_Pi = sp.simplify(A * W - R**2)
    Q_Pi = sp.simplify(GU**2 * W + 2 * GU * GW * R + GW**2 * A)
    D_Pi = sp.simplify(KB - Q_Pi / Delta_Pi)

    K_dyn = sp.Matrix(
        [
            [KB, -GU, -GW],
            [-GU, A, -R],
            [-GW, -R, W],
        ]
    )

    print("K_dyn =")
    sp.pprint(K_dyn)
    print("Delta_Pi =")
    sp.pprint(Delta_Pi)
    print("Q_Pi =")
    sp.pprint(Q_Pi)
    print("D_Pi =")
    sp.pprint(D_Pi)

    expect_zero("det(K_dyn) - Delta_Pi*D_Pi", sp.factor(K_dyn.det() - Delta_Pi * D_Pi))

    return KB, GU, GW, A, W, A0, W0, R, Pi, Delta_Pi, Q_Pi, D_Pi, K_dyn


# ---------------------------------------------------------------------------
# Part II. Exact inverse and dynamic susceptibilities
# ---------------------------------------------------------------------------

def verify_inverse_entries(GBundle):
    banner("PART II — EXACT DYNAMIC INVERSE AND SUSCEPTIBILITIES")

    KB, GU, GW, A, W, A0, W0, R, Pi, Delta_Pi, Q_Pi, D_Pi, K_dyn = GBundle

    inv = sp.simplify(K_dyn.inv())
    print("K_dyn^{-1} =")
    sp.pprint(inv)

    P_U = sp.simplify(GU * W + R * GW)
    P = sp.simplify(A * GW + R * GU)

    chi = {
        "qq": sp.simplify(inv[0, 0]),
        "qU": sp.simplify(inv[0, 1]),
        "qW": sp.simplify(inv[0, 2]),
        "UU": sp.simplify(inv[1, 1]),
        "UW": sp.simplify(inv[1, 2]),
        "WW": sp.simplify(inv[2, 2]),
    }

    expect_zero("chi_qq - 1/D_Pi", chi["qq"] - 1 / D_Pi)
    expect_zero("chi_qU - P_U/(Delta_Pi D_Pi)", chi["qU"] - P_U / (Delta_Pi * D_Pi))
    expect_zero("chi_qW - P/(Delta_Pi D_Pi)", chi["qW"] - P / (Delta_Pi * D_Pi))
    expect_zero(
        "chi_UU - (K_B W - G_W^2)/(Delta_Pi D_Pi)",
        chi["UU"] - (KB * W - GW**2) / (Delta_Pi * D_Pi),
    )
    expect_zero(
        "chi_UW - (K_B R + G_U G_W)/(Delta_Pi D_Pi)",
        chi["UW"] - (KB * R + GU * GW) / (Delta_Pi * D_Pi),
    )
    expect_zero(
        "chi_WW - (K_B A - G_U^2)/(Delta_Pi D_Pi)",
        chi["WW"] - (KB * A - GU**2) / (Delta_Pi * D_Pi),
    )

    print("\nDynamic susceptibility entries:")
    for name in ("qq", "qU", "qW", "UU", "UW", "WW"):
        print(f"chi_{name}(omega) =")
        sp.pprint(chi[name])

    return inv, chi, P_U, P


# ---------------------------------------------------------------------------
# Part III. Collinear-source theorem
# ---------------------------------------------------------------------------

def verify_collinear_source_law(GBundle, inv):
    banner("PART III — EXACT COLLINEAR-SOURCE SUSCEPTIBILITY LAW")

    KB, GU, GW, A, W, A0, W0, R, Pi, Delta_Pi, Q_Pi, D_Pi, K_dyn = GBundle

    s_q, s_U, s_W, S = sp.symbols("s_q s_U s_W S", real=True)
    s = sp.Matrix([s_q, s_U, s_W])
    J = S * s
    V_mix = sp.simplify(-sp.Rational(1, 2) * (J.T * inv * J)[0])

    P_U = sp.simplify(GU * W + R * GW)
    P = sp.simplify(A * GW + R * GU)
    N_s = sp.simplify(
        Delta_Pi * s_q**2
        + 2 * P_U * s_q * s_U
        + 2 * P * s_q * s_W
        + (KB * W - GW**2) * s_U**2
        + 2 * (KB * R + GU * GW) * s_U * s_W
        + (KB * A - GU**2) * s_W**2
    )

    print("V_mix(S s) =")
    sp.pprint(V_mix)
    print("N_s(omega) =")
    sp.pprint(N_s)

    expect_zero(
        "V_mix + 1/2 S^2 N_s/(Delta_Pi D_Pi)",
        sp.simplify(V_mix + sp.Rational(1, 2) * S**2 * N_s / (Delta_Pi * D_Pi)),
    )


# ---------------------------------------------------------------------------
# Part IV. Primitive source-family theorem
# ---------------------------------------------------------------------------

def verify_product_families(chi):
    banner("PART IV — DYNAMIC PRODUCT-FAMILY THEOREM")

    x, kappa = sp.symbols("x kappa", positive=True, real=True)
    betaQ, betaU, betaW = sp.symbols("beta_Q beta_U beta_W", real=True)

    SQ = x ** (-3)
    SY = sp.exp(-2 * kappa * x) / x
    J = sp.Matrix([betaQ * SQ, betaU * SY, betaW * SY])

    V_mix = sp.simplify(
        -sp.Rational(1, 2)
        * (
            chi["qq"] * J[0] ** 2
            + 2 * chi["qU"] * J[0] * J[1]
            + 2 * chi["qW"] * J[0] * J[2]
            + chi["UU"] * J[1] ** 2
            + 2 * chi["UW"] * J[1] * J[2]
            + chi["WW"] * J[2] ** 2
        )
    )

    C6 = sp.simplify(chi["qq"] * betaQ**2)
    C4 = sp.simplify(chi["qU"] * betaQ * betaU + chi["qW"] * betaQ * betaW)
    C2 = sp.simplify(
        chi["UU"] * betaU**2 + 2 * chi["UW"] * betaU * betaW + chi["WW"] * betaW**2
    )

    print("V_mix(x,omega) =")
    sp.pprint(V_mix)
    print("C6(omega) =")
    sp.pprint(C6)
    print("C4(omega) =")
    sp.pprint(C4)
    print("C2(omega) =")
    sp.pprint(C2)

    expect_zero(
        "dynamic product-family decomposition",
        sp.simplify(
            V_mix
            + sp.Rational(1, 2)
            * (
                C6 / x**6
                + 2 * C4 * sp.exp(-2 * kappa * x) / x**4
                + C2 * sp.exp(-4 * kappa * x) / x**2
            )
        ),
    )


# ---------------------------------------------------------------------------
# Part V. Outgoing-port derivative identity
# ---------------------------------------------------------------------------

def verify_outgoing_port_identity(GBundle, inv):
    banner("PART V — EXACT OUTGOING-PORT DERIVATIVE IDENTITY")

    KB, GU, GW, A, W, A0, W0, R, Pi, Delta_Pi, Q_Pi, D_Pi, K_dyn = GBundle

    Jq, Ju, Jw = sp.symbols("J_q J_U J_W", real=True)
    J = sp.Matrix([Jq, Ju, Jw])
    eW = sp.Matrix([0, 0, 1])
    EW = eW * eW.T

    dKinv_dPi = sp.diff(inv, Pi)
    expect_zero(
        "d(K_dyn^{-1})/dPi - K_dyn^{-1} E_W K_dyn^{-1}",
        sp.simplify(dKinv_dPi - inv * EW * inv),
    )

    V_mix = sp.simplify(-sp.Rational(1, 2) * (J.T * inv * J)[0])
    dVdPi = sp.simplify(sp.diff(V_mix, Pi))
    T_J = sp.simplify((eW.T * inv * J)[0])

    expect_zero("dV_mix/dPi + (1/2) T_J^2", sp.simplify(dVdPi + sp.Rational(1, 2) * T_J**2))

    print("T_J(omega) = e_W^T K_dyn^{-1} J =")
    sp.pprint(T_J)
    print("\nSo for a small outgoing-port load delta Pi, the exact first correction is")
    print("delta V_mix = -1/2 * deltaPi * T_J(omega)^2 + O(deltaPi^2).")


# ---------------------------------------------------------------------------
# Part VI. Conservative/static reduction check
# ---------------------------------------------------------------------------

def verify_static_reduction() -> None:
    banner("PART VI — CONSERVATIVE/STATIC REDUCTION CHECK")

    K, M, C, varpi, omega = sp.symbols("K M C varpi omega", positive=True, real=True)
    OU, OW, Pi = sp.symbols("Omega_U Omega_W Pi", positive=True, real=True)

    KB_expr = sp.simplify(K - M * omega**2 - C**2 / (varpi**2 - omega**2))
    A_expr = sp.simplify(OU**2 - omega**2)
    W_expr = sp.simplify(OW**2 - omega**2 - Pi)

    K_star = sp.simplify(K - C**2 / varpi**2)

    expect_zero("K_B(0) - K_star", sp.simplify(KB_expr.subs(omega, 0) - K_star))
    expect_zero("A(0) - Omega_U^2", sp.simplify(A_expr.subs(omega, 0) - OU**2))
    expect_zero(
        "W(0)|_{Pi=0} - Omega_W^2",
        sp.simplify(W_expr.subs({omega: 0, Pi: 0}) - OW**2),
    )

    print("So the dynamic bundle reduces exactly to the Stage-002 static bundle when")
    print("omega = 0 and Pi(0) = 0.")


# ---------------------------------------------------------------------------
# Run all parts
# ---------------------------------------------------------------------------

def main() -> None:
    g = verify_dynamic_bundle()
    inv, chi, PU, P = verify_inverse_entries(g)
    verify_collinear_source_law(g, inv)
    verify_product_families(chi)
    verify_outgoing_port_identity(g, inv)
    verify_static_reduction()

    banner("STAGE 003 SUMMARY")
    print("1. The dynamic one-port bundle is governed by the same denominator pair Delta_Pi and D_Pi.")
    print("2. Linear monochromatic driving does not create a new spatial kernel family.")
    print("3. The exact outgoing-port derivative is a perfect square transfer law.")
    print("4. Therefore the first passive/outgoing correction is phase-lag / pumping at linear order, not conservative barrier lowering.")
    print("5. Any surviving linear dynamic corridor must be a resonant dispersive enhancement of the existing short-range families.")


if __name__ == "__main__":
    main()
