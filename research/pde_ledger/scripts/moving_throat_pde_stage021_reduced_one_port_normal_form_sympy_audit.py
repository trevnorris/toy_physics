#!/usr/bin/env python3
"""
moving_throat_pde_stage4_maxwell_mixed_sympy_audit.py

SymPy audit for the first localized-Maxwell + mixed-sector moving-throat reduction.

Scope
-----
This script backs the next step after the matter-coupled Stage-3 reduction:

  • exact gauge-invariant mixed 4+1 Maxwell combinations E_w and C_a,
  • a reduced wall + brane-like Maxwell mode + mixed A_w/F_{mu w}/J^w-active mode,
  • exact elimination of the conservative internal Maxwell/mixed sector,
  • low-frequency conservative coefficients,
  • outgoing dressing of the mixed block and the induced wall-level odd term,
  • the compact outgoing l=2 PDE fingerprint,
  • a scalar derivative-coupling check showing why the dangerous i*omega scalar
    term need not reappear at wall level.

This is still a controlled reduced-sector calculation, not a full solution of the
coupled GNLS + localized Maxwell + moving-throat PDE.
"""

from __future__ import annotations

import sympy as sp
from sympy.calculus.euler import euler_equations

I = sp.I
sqrt = sp.sqrt


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
# Section I. Exact mixed-field gauge invariants from 4+1 Maxwell
# ---------------------------------------------------------------------------

def mixed_field_gauge_invariance_audit() -> None:
    banner("SECTION I — MIXED 4+1 MAXWELL FIELDS ARE GAUGE INVARIANT")

    t, x, w = sp.symbols("t x w", real=True)
    chi = sp.Function("chi")(t, x, w)
    A0 = sp.Function("A0")(t, x, w)
    Aa = sp.Function("Aa")(t, x, w)
    Aw = sp.Function("Aw")(t, x, w)

    Ew = -sp.diff(Aw, t) - sp.diff(A0, w)
    Ca = sp.diff(Aw, x) - sp.diff(Aa, w)

    # With the scalar-potential convention used in the project summaries,
    # A_0 transforms with a minus time derivative while spatial components
    # transform with plus spatial derivatives.
    A0p = A0 - sp.diff(chi, t)
    Aap = Aa + sp.diff(chi, x)
    Awp = Aw + sp.diff(chi, w)

    Ewp = -sp.diff(Awp, t) - sp.diff(A0p, w)
    Cap = sp.diff(Awp, x) - sp.diff(Aap, w)

    subbanner("I.1 — Mixed electric and mixed spatial fields")
    expect_zero("E_w gauge variation", sp.simplify(Ewp - Ew))
    expect_zero("C_a=F_aw gauge variation", sp.simplify(Cap - Ca))

    print("So the mixed fields E_w = F_w0 and C_a = F_aw are honest gauge-invariant observables.")
    print("This is the exact parent-theory reason the A_w / F_{mu w} sector can carry real reduced dynamics.")


# ---------------------------------------------------------------------------
# Section II. Conservative wall + brane-like Maxwell + mixed mode reduction
# ---------------------------------------------------------------------------

def conservative_maxwell_mixed_reduction_audit() -> dict[str, sp.Expr]:
    banner("SECTION II — CONSERVATIVE WALL + BRANE-MAXWELL + MIXED-MODE REDUCTION")

    t, omega = sp.symbols("t omega", real=True)
    Q = sp.Function("Q")(t)
    A = sp.Function("A")(t)
    W = sp.Function("W")(t)

    M, K = sp.symbols("M K", positive=True, real=True)
    OA, OW = sp.symbols("Omega_A Omega_W", positive=True, real=True)
    R, gA, gW = sp.symbols("R g_A g_W", real=True)

    Lred = (
        sp.Rational(1, 2) * M * sp.diff(Q, t) ** 2
        - sp.Rational(1, 2) * K * Q ** 2
        + sp.Rational(1, 2) * sp.diff(A, t) ** 2
        - sp.Rational(1, 2) * OA**2 * A ** 2
        + sp.Rational(1, 2) * sp.diff(W, t) ** 2
        - sp.Rational(1, 2) * OW**2 * W ** 2
        + R * A * W
        + gA * Q * A
        + gW * Q * W
    )

    EQ_Q = euler_equations(Lred, Q, [t])[0]
    EQ_A = euler_equations(Lred, A, [t])[0]
    EQ_W = euler_equations(Lred, W, [t])[0]

    subbanner("II.1 — Euler–Lagrange equations")
    expect_zero("Q equation", EQ_Q.lhs + M * sp.diff(Q, t, 2) + K * Q - gA * A - gW * W)
    expect_zero("A equation", EQ_A.lhs + sp.diff(A, t, 2) + OA**2 * A - R * W - gA * Q)
    expect_zero("W equation", EQ_W.lhs + sp.diff(W, t, 2) + OW**2 * W - R * A - gW * Q)

    Aker = OA**2 - omega**2
    Wker = OW**2 - omega**2
    Delta = sp.simplify(Aker * Wker - R**2)

    Sigma_cons = sp.simplify((gA**2 * Wker + 2 * gA * gW * R + gW**2 * Aker) / Delta)
    D_cons = sp.simplify(K - M * omega**2 - Sigma_cons)

    subbanner("II.2 — Exact conservative self-energy")
    print("Sigma_EM+mix^cons(omega) =")
    sp.pprint(Sigma_cons)
    print("D_cons(omega) = K - M omega^2 - Sigma_EM+mix^cons(omega) =")
    sp.pprint(D_cons)

    # Direct solve for A and W to verify the Schur-complement form.
    A_sol = sp.simplify((gA * Wker + gW * R) / Delta)
    W_sol = sp.simplify((gW * Aker + gA * R) / Delta)
    expect_zero("A exact solution residual", sp.simplify(Aker * A_sol - R * W_sol - gA))
    expect_zero("W exact solution residual", sp.simplify(Wker * W_sol - R * A_sol - gW))

    subbanner("II.3 — Low-frequency conservative coefficients")
    D0 = sp.symbols("D0", positive=True, real=True)
    S2 = sp.symbols("S2", positive=True, real=True)
    N0 = sp.symbols("N0", real=True)
    G2 = sp.symbols("G2", real=True)
    toy = (N0 - G2 * omega**2) / (D0 - S2 * omega**2 + omega**4)
    toy_series = sp.expand(sp.series(toy, omega, 0, 5).removeO())
    z0_toy = sp.simplify(toy_series.coeff(omega, 0))
    z2_toy = sp.simplify(toy_series.coeff(omega, 2))
    z4_toy = sp.simplify(toy_series.coeff(omega, 4))
    expect_zero("z0 formula", z0_toy - N0 / D0)
    expect_zero("z2 formula", z2_toy - (N0 * S2 - G2 * D0) / D0**2)
    expect_zero("z4 formula", z4_toy - (N0 * (S2**2 - D0) - S2 * G2 * D0) / D0**3)

    subs_dict = {
        D0: OA**2 * OW**2 - R**2,
        S2: OA**2 + OW**2,
        N0: gA**2 * OW**2 + 2 * gA * gW * R + gW**2 * OA**2,
        G2: gA**2 + gW**2,
    }
    Sigma_series = sp.expand(sp.series(Sigma_cons, omega, 0, 5).removeO())
    z0 = sp.simplify(Sigma_series.coeff(omega, 0))
    z2 = sp.simplify(Sigma_series.coeff(omega, 2))
    z4 = sp.simplify(Sigma_series.coeff(omega, 4))
    expect_zero("Sigma z0", z0 - (N0 / D0).subs(subs_dict))
    expect_zero("Sigma z2", z2 - ((N0 * S2 - G2 * D0) / D0**2).subs(subs_dict))
    expect_zero("Sigma z4", z4 - ((N0 * (S2**2 - D0) - S2 * G2 * D0) / D0**3).subs(subs_dict))

    print("z0^(EM+mix) =")
    sp.pprint(z0)
    print("z2^(EM+mix) =")
    sp.pprint(z2)
    print("z4^(EM+mix) =")
    sp.pprint(z4)

    return {
        "Sigma_cons": Sigma_cons,
        "D_cons": D_cons,
        "Aker": Aker,
        "Wker": Wker,
        "Delta": Delta,
    }


# ---------------------------------------------------------------------------
# Section III. Outgoing dressing of the mixed block
# ---------------------------------------------------------------------------

def outgoing_mixed_dressing_audit() -> dict[str, sp.Expr]:
    banner("SECTION III — OUTGOING DRESSING OF THE MIXED BLOCK")

    omega = sp.symbols("omega", real=True)
    OA, OW = sp.symbols("Omega_A Omega_W", positive=True, real=True)
    R, gA, gW = sp.symbols("R g_A g_W", real=True)
    Pi = sp.symbols("Pi", real=True)

    Aker = OA**2 - omega**2
    Wker = OW**2 - omega**2
    Delta = sp.simplify(Aker * Wker - R**2)

    Sigma_cons = sp.simplify((gA**2 * Wker + 2 * gA * gW * R + gW**2 * Aker) / Delta)
    Sigma_full = sp.simplify((gA**2 * (Wker - Pi) + 2 * gA * gW * R + gW**2 * Aker) / (Aker * (Wker - Pi) - R**2))
    Sigma_first = sp.expand(sp.series(Sigma_full, Pi, 0, 2).removeO())
    N_omega = sp.simplify((Sigma_first - Sigma_cons) / Pi)

    subbanner("III.1 — First-order outgoing correction")
    print("Sigma_full(omega) =")
    sp.pprint(Sigma_full)
    print("N(omega) from Sigma_full = Sigma_cons + Pi*N + O(Pi^2):")
    sp.pprint(N_omega)

    expect_zero(
        "N(omega) compact formula",
        N_omega - ((Aker * gW + R * gA) ** 2) / Delta**2,
    )

    N0 = sp.simplify(N_omega.subs(omega, 0))
    expect_zero(
        "N(0) positive-square form",
        N0 - (OA**2 * gW + R * gA) ** 2 / (OA**2 * OW**2 - R**2) ** 2,
    )

    print("N(0) =")
    sp.pprint(N0)
    print("So the mixed-sector outgoing port is transferred to the wall with a manifestly nonnegative coefficient.")

    # Translate to wall-operator convention D = K - M omega^2 - Sigma.
    Gamma_port = sp.symbols("Gamma_port", positive=True, real=True)
    a, c_s = sp.symbols("a c_s", positive=True, real=True)
    Dcorr = sp.simplify(-I * Gamma_port * omega**5 * N0)
    expect_zero(
        "delta D_2^(odd) composed from Section III N(0) closed form and Section IV Gamma5 = a^5/(27 c_s^5)",
        Dcorr.subs(Gamma_port, a**5 / (27 * c_s**5))
        - (
            -I
            * ((OA**2 * gW + R * gA) ** 2 / (OA**2 * OW**2 - R**2) ** 2)
            * a**5
            / (27 * c_s**5)
            * omega**5
        ),
    )
    print("If Pi_out = + i Gamma_port omega^5 + O(omega^7), then")
    print("delta D_wall^(odd) =")
    sp.pprint(Dcorr)
    print("In this wall-operator convention the outgoing branch therefore appears with a negative imaginary sign.")
    print("Equivalently, the normalized response/admittance carries the positive +i omega^5 sign used in the 2.5PN audit.")

    return {
        "N0": N0,
    }


# ---------------------------------------------------------------------------
# Section IV. Compact outgoing l=2 PDE fingerprint
# ---------------------------------------------------------------------------

def outgoing_l2_fingerprint_audit() -> dict[str, sp.Expr]:
    banner("SECTION IV — COMPACT OUTGOING l=2 FINGERPRINT")

    k, a, c_s, omega = sp.symbols("k a c_s omega", positive=True, real=True)
    za = sp.symbols("za", positive=True, real=True)

    j2a = (3 / za**3 - 1 / za) * sp.sin(za) - 3 * sp.cos(za) / za**2
    y2a = -(3 / za**3 - 1 / za) * sp.cos(za) - 3 * sp.sin(za) / za**2
    h2a = sp.simplify(j2a + I * y2a)
    Lambda2 = sp.simplify((k * sp.diff(h2a, za) / h2a).subs(za, k * a))
    Lambda2_series = sp.series(Lambda2, k, 0, 7).removeO()
    Y2 = sp.simplify(sp.series(1 / Lambda2_series, k, 0, 6).removeO())
    Y2_static = sp.simplify(Y2.subs(k, 0))
    Y2_hat = sp.simplify(sp.expand(Y2 / Y2_static))
    Y2_hat_series = sp.series(Y2_hat, k, 0, 6).removeO()
    Y2_hat_omega = sp.expand(Y2_hat_series.subs(k, omega / c_s))

    subbanner("IV.1 — Normalized outgoing l=2 response")
    print("Lambda2(k) =")
    sp.pprint(Lambda2_series)
    print("Y2_hat(omega) =")
    sp.pprint(Y2_hat_omega)

    expect_zero(
        "Y2_hat minimal branch",
        Y2_hat_omega
        - (
            1
            + a**2 * omega**2 / (9 * c_s**2)
            + 4 * a**4 * omega**4 / (81 * c_s**4)
            + I * a**5 * omega**5 / (27 * c_s**5)
        ),
    )

    Gamma5_port = sp.simplify(sp.expand(Y2_hat_omega).coeff(omega, 5) / I)
    expect_zero("Gamma5_port - a^5/(27 c_s^5)", Gamma5_port - a**5 / (27 * c_s**5))

    print("Gamma5_port =")
    sp.pprint(Gamma5_port)
    print("So the compact outgoing l=2 branch starts at + i omega^5 with coefficient a^5/(27 c_s^5).")

    return {
        "Gamma5_port": Gamma5_port,
    }


# ---------------------------------------------------------------------------
# Section V. Scalar derivative-coupling check
# ---------------------------------------------------------------------------

def scalar_derivative_coupling_audit() -> None:
    banner("SECTION V — DERIVATIVE-COUPLED SCALAR MIXED OUTLET STARTS AT i*omega^3")

    omega = sp.symbols("omega", real=True)
    OA, OW = sp.symbols("Omega_A Omega_W", positive=True, real=True)
    R, eta, gamma1 = sp.symbols("R eta gamma1", positive=True, real=True)

    # Breathing/scalar lane with no direct non-derivative wall coupling into the port-active block.
    gA = sp.Integer(0)
    gW = eta * omega
    Aker = OA**2 - omega**2
    Wker = OW**2 - omega**2
    Delta = sp.simplify(Aker * Wker - R**2)
    N_omega = sp.simplify((Aker * gW + R * gA) ** 2 / Delta**2)
    N_series = sp.expand(sp.series(N_omega, omega, 0, 3).removeO())
    print("N_scalar(omega) =")
    sp.pprint(N_series)

    expect_zero(
        "N_scalar leading term",
        N_series - eta**2 * OA**4 * omega**2 / (OA**2 * OW**2 - R**2) ** 2,
    )

    Pi0 = I * gamma1 * omega
    deltaD0 = sp.expand(Pi0 * N_series)
    print("Pi0_out * N_scalar =")
    sp.pprint(deltaD0)
    expect_zero(
        "scalar odd order",
        sp.simplify(deltaD0 - I * gamma1 * eta**2 * OA**4 * omega**3 / (OA**2 * OW**2 - R**2) ** 2),
    )

    print("So a derivative-coupled scalar mixed outlet converts the naive i*omega port law into an i*omega^3 wall correction.")
    print("That matches the scalar rescue pattern isolated in the 2.5PN audit.")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    mixed_field_gauge_invariance_audit()
    conservative_maxwell_mixed_reduction_audit()
    mixed = outgoing_mixed_dressing_audit()
    out = outgoing_l2_fingerprint_audit()
    scalar_derivative_coupling_audit()

    banner("FINAL STAGE-21 LEDGER")
    print("Verified with SymPy:")
    print("  • the mixed 4+1 Maxwell fields E_w = F_w0 and C_a = F_aw are gauge invariant;")
    print("  • a reduced wall + brane-like Maxwell mode + mixed A_w/F_{mu w}/J^w-active mode has the exact conservative self-energy")
    print("    Sigma_EM+mix^cons = [g_A^2 W + 2 g_A g_W R + g_W^2 A] / (A W - R^2);")
    print("  • its low-frequency coefficients are explicit rational functions of the conservative pole data and couplings;")
    print("  • dressing the mixed block by a retarded port Pi_out transfers the odd part to the wall with coefficient")
    print("    N(0) = (Omega_A^2 g_W + R g_A)^2 / (Omega_A^2 Omega_W^2 - R^2)^2 >= 0;")
    print("  • the compact outgoing l=2 branch has the normalized fingerprint")
    print("    Y2_hat = 1 + a^2 omega^2/(9 c_s^2) + 4 a^4 omega^4/(81 c_s^4) + i a^5 omega^5/(27 c_s^5) + ...;")
    print("  • therefore the first wall-level odd quadrupole coefficient is the outgoing-port coefficient")
    print("    a^5/(27 c_s^5) multiplied by the conservative mixed-sector transfer factor N(0);")
    print("  • if the scalar outlet enters only through a derivative mixed coupling, the wall-level scalar odd term starts at i*omega^3 rather than i*omega.")


if __name__ == "__main__":
    main()
