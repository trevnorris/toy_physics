#!/usr/bin/env python3
"""
same_charge_barrier_audit_stage002_one_port_mixed_bundle_kernel_sympy_audit.py

Stage 002 — one-port mixed-bundle static kernel audit.

What this script does
---------------------
1. Rebuilds the static one-port wall / mixed-bundle stiffness matrix after the
   stable BdG support mode has been integrated out into the effective wall
   stiffness K_* = K - C^2/varpi^2.
2. Verifies the exact determinant identity

       det(K_red) = Delta * D0,

   where
       Delta = Omega_U^2 Omega_W^2 - R^2,
       D0    = K_* - Q/Delta,
       Q     = G_U^2 Omega_W^2 + 2 G_U G_W R + G_W^2 Omega_U^2.
3. Computes the exact inverse and derives the static susceptibility law

       delta V = -1/2 J^T K_red^{-1} J.

4. Proves the exact general collinear-source formula

       J(r) = S(r) s
       ->
       delta V(r) = -1/2 chi_s S(r)^2,

   with chi_s = N_s / (Delta D0).
5. Records the basis susceptibilities and the identities tying them directly to
   the same outgoing-load factor and prefactor data used in the 5PN / 2.5PN
   normalization chain.
6. Verifies the product-kernel decomposition for the two primitive reduced load
   families used in the same-charge audit:

       S_Q(x) = x^{-3},
       S_Y(x) = exp(-2 kappa x)/x.

   The resulting static mixed correction can only generate the three product
   families
       x^{-6},
       exp(-2 kappa x)/x^4,
       exp(-4 kappa x)/x^2.
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


# ---------------------------------------------------------------------------
# Part I. Static one-port mixed bundle and exact determinant
# ---------------------------------------------------------------------------

def verify_static_bundle() -> tuple[sp.Symbol, sp.Symbol, sp.Symbol, sp.Symbol, sp.Symbol, sp.Symbol, sp.Symbol, sp.Symbol, sp.Expr, sp.Expr, sp.Expr, sp.Matrix]:
    banner("PART I — STATIC ONE-PORT MIXED BUNDLE")

    Kst, GU, GW, OU, OW, R = sp.symbols(
        "K_star G_U G_W Omega_U Omega_W R", positive=True, real=True
    )

    Delta = sp.simplify(OU**2 * OW**2 - R**2)
    Q = sp.simplify(GU**2 * OW**2 + 2 * GU * GW * R + GW**2 * OU**2)
    P = sp.simplify(OU**2 * GW + R * GU)
    PU = sp.simplify(GU * OW**2 + R * GW)
    D0 = sp.simplify(Kst - Q / Delta)

    K_red = sp.Matrix(
        [
            [Kst, -GU, -GW],
            [-GU, OU**2, -R],
            [-GW, -R, OW**2],
        ]
    )

    print("K_red =")
    sp.pprint(K_red)
    print("Delta =")
    sp.pprint(Delta)
    print("Q =")
    sp.pprint(Q)
    print("P =")
    sp.pprint(P)
    print("P_U =")
    sp.pprint(PU)
    print("D0 =")
    sp.pprint(D0)

    expect_zero("det(K_red) - Delta*D0", sp.factor(K_red.det() - Delta * D0))

    return Kst, GU, GW, OU, OW, R, Delta, D0, Q, P, PU, K_red


# ---------------------------------------------------------------------------
# Part II. Exact inverse and static susceptibility kernel
# ---------------------------------------------------------------------------

def verify_inverse_and_susceptibilities(
    Kst: sp.Symbol,
    GU: sp.Symbol,
    GW: sp.Symbol,
    OU: sp.Symbol,
    OW: sp.Symbol,
    R: sp.Symbol,
    Delta: sp.Symbol,
    D0: sp.Expr,
    Q: sp.Expr,
    P: sp.Expr,
    PU: sp.Expr,
    K_red: sp.Matrix,
) -> tuple[sp.Matrix, dict[str, sp.Expr]]:
    banner("PART II — EXACT INVERSE AND STATIC SUSCEPTIBILITIES")

    inv = sp.simplify(K_red.inv())
    print("K_red^{-1} =")
    sp.pprint(inv)

    chi = {
        "qq": sp.simplify(inv[0, 0]),
        "qU": sp.simplify(inv[0, 1]),
        "qW": sp.simplify(inv[0, 2]),
        "UU": sp.simplify(inv[1, 1]),
        "UW": sp.simplify(inv[1, 2]),
        "WW": sp.simplify(inv[2, 2]),
    }

    expect_zero("chi_qq - 1/D0", chi["qq"] - 1 / D0)
    expect_zero("chi_qU - P_U/(Delta D0)", chi["qU"] - PU / (Delta * D0))
    expect_zero("chi_qW - P/(Delta D0)", chi["qW"] - P / (Delta * D0))
    expect_zero(
        "chi_UU - (K_* Omega_W^2 - G_W^2)/(Delta D0)",
        chi["UU"] - (Kst * OW**2 - GW**2) / (Delta * D0),
    )
    expect_zero(
        "chi_UW - (K_* R + G_U G_W)/(Delta D0)",
        chi["UW"] - (Kst * R + GU * GW) / (Delta * D0),
    )
    expect_zero(
        "chi_WW - (K_* Omega_U^2 - G_U^2)/(Delta D0)",
        chi["WW"] - (Kst * OU**2 - GU**2) / (Delta * D0),
    )

    print("\nBasis susceptibilities:")
    for name in ("qq", "qU", "qW", "UU", "UW", "WW"):
        print(f"chi_{name} =")
        sp.pprint(chi[name])

    return inv, chi


# ---------------------------------------------------------------------------
# Part III. General collinear-source formula
# ---------------------------------------------------------------------------

def verify_general_source_formula(
    Kst: sp.Symbol,
    GU: sp.Symbol,
    GW: sp.Symbol,
    OU: sp.Symbol,
    OW: sp.Symbol,
    R: sp.Symbol,
    Delta: sp.Symbol,
    D0: sp.Expr,
    inv: sp.Matrix,
) -> None:
    banner("PART III — GENERAL COLLINEAR-SOURCE LAW")

    s_q, s_U, s_W, S = sp.symbols("s_q s_U s_W S", real=True)
    s = sp.Matrix([s_q, s_U, s_W])
    J = S * s
    delta_V = sp.simplify(-sp.Rational(1, 2) * (J.T * inv * J)[0])

    N_s = sp.simplify(
        Delta * s_q**2
        + 2 * (GU * OW**2 + R * GW) * s_q * s_U
        + 2 * (OU**2 * GW + R * GU) * s_q * s_W
        + (Kst * OW**2 - GW**2) * s_U**2
        + 2 * (Kst * R + GU * GW) * s_U * s_W
        + (Kst * OU**2 - GU**2) * s_W**2
    )

    print("delta V(S s) =")
    sp.pprint(delta_V)
    print("N_s =")
    sp.pprint(N_s)

    expect_zero(
        "delta V + (1/2) S^2 N_s/(Delta D0)",
        sp.simplify(delta_V + sp.Rational(1, 2) * S**2 * N_s / (Delta * D0)),
    )

    print(
        "\nSo for any collinear reduced source J(r)=S(r) s, the exact static correction is"
    )
    print("delta V(r) = -1/2 * chi_s * S(r)^2 with chi_s = N_s/(Delta D0).")


# ---------------------------------------------------------------------------
# Part IV. Exact bridge to outgoing-load / 5PN normalization variables
# ---------------------------------------------------------------------------

def verify_bridge_to_prefactor(
    Delta: sp.Symbol,
    D0: sp.Expr,
    P: sp.Expr,
    chi: dict[str, sp.Expr],
) -> None:
    banner("PART IV — BRIDGE TO OUTGOING-LOAD / 5PN VARIABLES")

    Lambda = sp.simplify(P / Delta)
    N0 = sp.simplify(P**2 / Delta**2)
    P0 = sp.simplify(N0 / D0)

    print("Lambda =")
    sp.pprint(Lambda)
    print("N0 =")
    sp.pprint(N0)
    print("P0 =")
    sp.pprint(P0)

    expect_zero("chi_qW - Lambda/D0", sp.simplify(chi["qW"] - Lambda / D0))
    expect_zero("N0 - Lambda^2", sp.simplify(N0 - Lambda**2))
    expect_zero("P0 - Lambda^2/D0", sp.simplify(P0 - Lambda**2 / D0))
    expect_zero("chi_qW^2 - P0/D0", sp.simplify(chi["qW"] ** 2 - P0 / D0))

    print("\nThe wall-mixed cross susceptibility is controlled by the same outgoing-load factor")
    print("Lambda = P/Delta that feeds the 5PN / 2.5PN prefactor chain.")


# ---------------------------------------------------------------------------
# Part V. Product-kernel theorem for primitive reduced load families
# ---------------------------------------------------------------------------

def verify_product_kernel_theorem(chi: dict[str, sp.Expr]) -> None:
    banner("PART V — PRODUCT-KERNEL THEOREM")

    x, kappa = sp.symbols("x kappa", positive=True, real=True)
    beta_Q, beta_U, beta_W = sp.symbols("beta_Q beta_U beta_W", real=True)

    S_Q = x ** (-3)
    S_Y = sp.exp(-2 * kappa * x) / x

    delta_V = sp.simplify(
        -sp.Rational(1, 2)
        * (
            chi["qq"] * (beta_Q * S_Q) ** 2
            + 2 * chi["qU"] * beta_Q * beta_U * S_Q * S_Y
            + 2 * chi["qW"] * beta_Q * beta_W * S_Q * S_Y
            + chi["UU"] * (beta_U * S_Y) ** 2
            + 2 * chi["UW"] * beta_U * beta_W * S_Y**2
            + chi["WW"] * (beta_W * S_Y) ** 2
        )
    )

    C6 = sp.simplify(chi["qq"] * beta_Q**2)
    C4 = sp.simplify(chi["qU"] * beta_Q * beta_U + chi["qW"] * beta_Q * beta_W)
    C2 = sp.simplify(
        chi["UU"] * beta_U**2 + 2 * chi["UW"] * beta_U * beta_W + chi["WW"] * beta_W**2
    )

    target = sp.simplify(
        -sp.Rational(1, 2)
        * (
            C6 / x**6
            + 2 * C4 * sp.exp(-2 * kappa * x) / x**4
            + C2 * sp.exp(-4 * kappa * x) / x**2
        )
    )

    print("delta V_mix(x) =")
    sp.pprint(delta_V)
    print("C6 =")
    sp.pprint(C6)
    print("C4 =")
    sp.pprint(C4)
    print("C2 =")
    sp.pprint(C2)

    expect_zero("product-kernel decomposition", sp.simplify(delta_V - target))

    print("\nSo the one-port static bundle can only generate the three product families:")
    print("  x^(-6),   exp(-2 kappa x)/x^4,   exp(-4 kappa x)/x^2")
    print("from the primitive reduced source pair S_Q=x^{-3}, S_Y=exp(-2 kappa x)/x.")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    data = verify_static_bundle()
    inv, chi = verify_inverse_and_susceptibilities(*data)
    verify_general_source_formula(*data[:8], inv)
    verify_bridge_to_prefactor(data[6], data[7], data[9], chi)
    verify_product_kernel_theorem(chi)

    banner("STAGE 002 SUMMARY")
    print("1. The first PDE-constrained static mixed correction is an exact quadratic")
    print("   susceptibility kernel of the one-port bundle.")
    print("2. Its denominator is Delta*D0, the same bundle object that feeds the")
    print("   Stage-8/9/18 normalization chain.")
    print("3. For primitive quadrupole and Yukawa reduced drives, the bundle cannot")
    print("   create a new slower-than-source attractive family; it only produces")
    print("   the product kernels x^-6, e^{-2 kappa x}/x^4, e^{-4 kappa x}/x^2.")
    print("4. The wall-mixed cross susceptibility is controlled by the outgoing-load")
    print("   factor Lambda=P/Delta, so static same-charge softening and the 5PN/")
    print("   2.5PN bundle are not independent knobs.")


if __name__ == "__main__":
    main()
