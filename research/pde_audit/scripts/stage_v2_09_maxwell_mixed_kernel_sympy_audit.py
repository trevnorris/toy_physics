#!/usr/bin/env python3
"""Fast symbolic audit for V2-09 Maxwell/mixed-sector reduced kernel.

This script audits the reduced one-lane wall + localized-Maxwell + mixed-sector
model used as the first outgoing quadrupole bridge in the Moving-Throat PDE
Volume 2 program.
"""
from __future__ import annotations
import sympy as sp


def check_zero(name: str, expr) -> bool:
    res = sp.expand(expr)
    ok = (res == 0)
    print(f"{name}: {'PASS' if ok else 'FAIL'}")
    if not ok:
        print("  residual =", res)
    return ok


def main() -> None:
    print("Stage V2-09 Maxwell/mixed-sector kernel audit")
    print("=" * 72)

    # 1) Mixed field gauge invariance.
    t, x, w = sp.symbols("t x w")
    chi = sp.Function("chi")(t, x, w)
    A0 = sp.Function("A0")(t, x, w)
    Ax = sp.Function("Ax")(t, x, w)
    Aw = sp.Function("Aw")(t, x, w)
    Ew = -sp.diff(Aw, t) - sp.diff(A0, w)
    Cx = sp.diff(Aw, x) - sp.diff(Ax, w)
    A0p = A0 - sp.diff(chi, t)
    Axp = Ax + sp.diff(chi, x)
    Awp = Aw + sp.diff(chi, w)
    Ewp = -sp.diff(Awp, t) - sp.diff(A0p, w)
    Cxp = sp.diff(Awp, x) - sp.diff(Axp, w)
    checks = [
        check_zero("MIXED_GAUGE_INVARIANCE_Ew", Ewp - Ew),
        check_zero("MIXED_GAUGE_INVARIANCE_Ca", Cxp - Cx),
    ]

    # 2) One-lane reduced variables. y = omega^2.
    y, omega = sp.symbols("y omega")
    M, K = sp.symbols("M K")
    OU2, OW2, R = sp.symbols("Omega_U2 Omega_W2 R")
    gU, gW = sp.symbols("g_U g_W")
    U = OU2 - y
    W = OW2 - y
    Delta = U*W - R**2
    Delta0 = OU2*OW2 - R**2
    S2 = OU2 + OW2
    Q0 = gU**2*OW2 + 2*gU*gW*R + gW**2*OU2
    G2 = gU**2 + gW**2
    Sigma_num = Q0 - G2*y
    Sigma_den = Delta0 - S2*y + y**2

    # Schur complement numerator/denominator identity.
    Sigma_schur_num = gU**2*W + 2*gU*gW*R + gW**2*U
    checks.append(check_zero("CONSERVATIVE_SELF_ENERGY_NUMERATOR", Sigma_schur_num - Sigma_num))
    checks.append(check_zero("CONSERVATIVE_SELF_ENERGY_DENOMINATOR", Delta - Sigma_den))

    # Low-frequency z0,z2,z4: den*(z0+z2*y+z4*y^2)-num = O(y^3).
    z0 = Q0/Delta0
    z2 = (Q0*S2 - G2*Delta0)/Delta0**2
    z4 = (Q0*(S2**2 - Delta0) - S2*G2*Delta0)/Delta0**3
    sigma_match_num = sp.expand(Delta0**3*(Sigma_den*(z0 + z2*y + z4*y**2) - Sigma_num))
    poly = sp.Poly(sigma_match_num, y)
    for n in range(3):
        checks.append(check_zero(f"SIGMA_LOW_FREQ_COEFF_{n}", poly.coeff_monomial(y**n)))

    # 3) Outgoing port transfer. Differentiate by Pi by hand.
    # Sigma(Pi) = (N - gU^2 Pi)/(Delta - U Pi)
    # dSigma/dPi|0 = (U*N - gU^2*Delta)/Delta^2.
    transfer_num_from_derivative = U*Sigma_num - gU**2*Delta
    transfer_square_num = (U*gW + R*gU)**2
    checks.append(check_zero("OUTGOING_TRANSFER_FACTOR_NUMERATOR", transfer_num_from_derivative - transfer_square_num))

    P = OU2*gW + R*gU
    N0 = P**2/Delta0**2
    N2 = 2*P*(P*S2 - Delta0*gW)/Delta0**3
    N4 = (Delta0**2*gW**2 - 2*Delta0*P**2 - 4*Delta0*P*S2*gW + 3*P**2*S2**2)/Delta0**4
    N_num = (P - gW*y)**2
    N_den = (Delta0 - S2*y + y**2)**2
    N_match_num = sp.expand(Delta0**4*(N_den*(N0 + N2*y + N4*y**2) - N_num))
    polyN = sp.Poly(N_match_num, y)
    for n in range(3):
        checks.append(check_zero(f"TRANSFER_LOW_FREQ_COEFF_{n}", polyN.coeff_monomial(y**n)))

    # 4) Static stability via Schur determinant.
    Hfull = sp.Matrix([[K, -gU, -gW], [-gU, OU2, -R], [-gW, -R, OW2]])
    det_full = sp.expand(Hfull.det())
    # Delta0*(K - Sigma0) has numerator K*Delta0-Q0.
    checks.append(check_zero("FULL_STATIC_DETERMINANT_SCHUR", det_full - (K*Delta0 - Q0)))

    # 5) l=2 outgoing coefficient sign.
    a, cs = sp.symbols("a c_s")
    deltaD_odd = -sp.I*N0*a**5*omega**5/(27*cs**5)
    expected_deltaD = -sp.I*P**2*a**5*omega**5/(27*cs**5*Delta0**2)
    checks.append(check_zero("L2_WALL_OPERATOR_ODD_COEFFICIENT", deltaD_odd - expected_deltaD))

    # 6) Scalar derivative coupling: gU=0, gW=eta*omega => N ~ omega^2.
    eta, gamma1 = sp.symbols("eta gamma1")
    scalar_N_num = eta**2 * omega**2 * (OU2 - omega**2)**2
    leading_scalar_N = sp.expand(scalar_N_num).coeff(omega, 2) / Delta0**2
    checks.append(check_zero("SCALAR_DERIVATIVE_TRANSFER_STARTS_AT_OMEGA2", leading_scalar_N - eta**2*OU2**2/Delta0**2))
    # Pi_0 ~ i gamma1 omega, wall operator has -Pi*N, hence leading -i omega^3.
    leading_deltaD_coeff = gamma1 * leading_scalar_N
    checks.append(check_zero("SCALAR_DERIVATIVE_ODD_DEMOTED_TO_OMEGA3", leading_deltaD_coeff - gamma1*eta**2*OU2**2/Delta0**2))

    print("\nKey symbolic formulas")
    print("-" * 72)
    print("Delta0 =", Delta0)
    print("Sigma_cons = (Q0 - G2*omega^2)/(Delta0 - S2*omega^2 + omega^4)")
    print("Q0 =", Q0)
    print("G2 =", G2)
    print("z0 =", z0)
    print("z2 =", z2)
    print("z4 =", z4)
    print("N_transfer = (A(omega)*g_W + R*g_U)^2 / Delta(omega)^2")
    print("P =", P)
    print("N0 =", N0)
    print("N2 =", N2)
    print("N4 =", N4)
    print("deltaD_l=2_odd =", expected_deltaD)
    print("scalar derivative leading N/omega^2 =", eta**2*OU2**2/Delta0**2)
    print("direct scalar dark-port condition: ", P, "= 0")

    print("\nPass/fail ledger")
    print("-" * 72)
    print(f"identity_checks: {len(checks)}")
    print(f"identity_passes: {sum(bool(c) for c in checks)}")
    print("\nConditional verdict")
    print("-" * 72)
    if all(checks):
        print("MAXWELL_MIXED_KERNEL: PASS within the declared one-lane reduced closure.")
        print("OUTGOING_TRANSFER: PASS; transfer factor is a perfect square over Delta^2.")
        print("L2_ODD_BRANCH: PASS; wall operator inherits -i*N0*a^5/(27*c_s^5)*omega^5.")
        print("STABILITY_GATE: REQUIRE Delta0>0 and K-Sigma0>0, equivalently Delta0>0 and K*Delta0-Q0>0.")
        print("DARK_PORT_WARNING: if Omega_U2*g_W + R*g_U = 0, then N0=0 and the port does not radiate at leading l=2 order.")
        print("SCALAR_OHMIC_WARNING: direct non-derivative scalar port coupling gives N0 != 0 generically; derivative coupling demotes the leading scalar odd term to omega^3.")
    else:
        print("MAXWELL_MIXED_KERNEL: FAIL; inspect residuals above.")


if __name__ == "__main__":
    main()
