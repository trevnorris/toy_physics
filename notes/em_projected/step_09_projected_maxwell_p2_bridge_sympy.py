#!/usr/bin/env python3
"""
Projected Maxwell -> grouped P2 bridge

This script symbolically connects the near-throat projection-first Maxwell
corrections to the grouped real P2 / full-bundle language used by the moving-
throat PDE program.

Main jobs:
1. Start from the isotropic full-bundle coefficients:
      D0 = K - B0 - Z0
      D2 = -(M + B2 + Z2)
      D4 = -(B4 + Z4)
      P0 = N0 / D0
      P2 = (D0 N2 - 2 D2 N0) / D0^2
      P4 = (D0^2 N4 - 2 D0(D2 N2 + D4 N0) + 3 D2^2 N0) / D0^3
2. Add projected-Maxwell near-mouth corrections to the Maxwell/mixed moments:
      Z_n -> Z_n + eps z_n
      N_n -> N_n + eps n_n
   where eps may later be interpreted as the mouth scale ell or, for a
   symmetric interior observer, as sigma^2.
3. Expand the induced corrections to u2, u4, P0, P2, P4, the one-pole defect,
   and the isotropic compatibility equation.
4. Add a weak-axisymmetric Y20 grouped signature:
      (20,21,22) = (1, 1/2, -1)
   and derive the resulting grouped anomalies and prefactor-slope Xi1.

The output is designed to be read as a theorem-style audit notebook.
"""

from __future__ import annotations

import sympy as sp
from sympy.physics.wigner import gaunt


def banner(title: str) -> str:
    return "\n" + "=" * 96 + f"\n{title}\n" + "=" * 96 + "\n"


def line(label: str, expr: sp.Expr | str) -> str:
    return f"{label}\n  {expr}\n"


def assert_zero(label: str, expr: sp.Expr) -> None:
    residue = sp.factor(sp.together(sp.simplify(expr)))
    if residue != 0:
        raise AssertionError(f"{label} failed: {sp.sstr(residue)}")


def assert_nonzero(label: str, expr: sp.Expr) -> None:
    value = sp.factor(sp.together(sp.simplify(expr)))
    if value == 0:
        raise AssertionError(f"{label} unexpectedly vanished")


def real_y20_square_ratio(m: int) -> sp.Expr:
    base = sp.simplify(gaunt(2, 2, 2, 0, 0, 0))
    if m == 0:
        return sp.Integer(1)
    same_sign = sp.simplify(gaunt(2, 2, 2, 0, m, m))
    if same_sign != 0:
        raise AssertionError(f"Real-harmonic same-sign cross term should vanish for m={m}: {same_sign}")
    return sp.simplify((sp.Integer(-1) ** m) * gaunt(2, 2, 2, 0, m, -m) / base)


def main() -> None:
    D0, D2, D4, N0, N2, N4 = sp.symbols("D0 D2 D4 N0 N2 N4", nonzero=True)
    z0, z2, z4, n0, n2, n4 = sp.symbols("z0 z2 z4 n0 n2 n4")
    eps = sp.symbols("eps")

    # Isotropic projected-Maxwell corrections:
    D0p = D0 - eps * z0
    D2p = D2 - eps * z2
    D4p = D4 - eps * z4
    N0p = N0 + eps * n0
    N2p = N2 + eps * n2
    N4p = N4 + eps * n4

    u2p = -D2p / D0p
    u4p = (D2p**2 - D0p * D4p) / D0p**2
    P0p = N0p / D0p
    P2p = (D0p * N2p - 2 * D2p * N0p) / D0p**2
    P4p = (D0p**2 * N4p - 2 * D0p * (D2p * N2p + D4p * N0p) + 3 * D2p**2 * N0p) / D0p**3

    def lin(expr: sp.Expr) -> sp.Expr:
        return sp.expand(sp.series(sp.expand(expr), eps, 0, 2).removeO())

    u2_lin = lin(u2p)
    u4_lin = lin(u4p)
    P0_lin = lin(P0p)
    P2_lin = lin(P2p)
    P4_lin = lin(P4p)

    du2 = sp.simplify((u2_lin + D2 / D0) / eps)
    du4 = sp.simplify((u4_lin - (D2**2 - D0 * D4) / D0**2) / eps)
    dP0 = sp.simplify((P0_lin - N0 / D0) / eps)
    dP2 = sp.simplify((P2_lin - (D0 * N2 - 2 * D2 * N0) / D0**2) / eps)
    dP4 = sp.simplify(
        (P4_lin - (D0**2 * N4 - 2 * D0 * (D2 * N2 + D4 * N0) + 3 * D2**2 * N0) / D0**3) / eps
    )
    assert_zero("delta u2", du2 - (D0*z2 - D2*z0) / D0**2)
    assert_zero("delta P0", dP0 - (D0*n0 + N0*z0) / D0**2)

    B4, Z4, M, B2, Z2, Ptarget = sp.symbols("B4 Z4 M B2 Z2 Ptarget", nonzero=True)
    S = sp.symbols("S", nonzero=True)
    T = sp.symbols("T", nonzero=True)

    # One-pole defect and compatibility surface
    pole_defect = D0 * (B4 + Z4) - 3 * (M + B2 + Z2) ** 2
    pole_defect_p = (D0 - eps * z0) * (B4 + Z4 + eps * z4) - 3 * (M + B2 + Z2 + eps * z2) ** 2
    d_pole = sp.simplify((lin(pole_defect_p) - pole_defect) / eps)
    assert_zero("one-pole first variation", d_pole - (D0*z4 - z0*(B4 + Z4) - 6*z2*(M + B2 + Z2)))

    K, B0, Z0slot = sp.symbols("K B0 Z0slot", nonzero=True)
    K_one_pole_p = sp.solve(
        sp.Eq((K - B0 - (Z0slot + eps * z0)) * (T + eps * z4), 3 * (S + eps * z2) ** 2),
        K,
    )[0]
    K_norm_p = sp.solve(
        sp.Eq((N0 + eps * n0) / (K - B0 - (Z0slot + eps * z0)), Ptarget),
        K,
    )[0]
    K_one_pole = K_one_pole_p.subs(eps, 0)
    K_norm = K_norm_p.subs(eps, 0)
    dK_one_pole = sp.simplify((lin(K_one_pole_p) - K_one_pole) / eps)
    dK_norm = sp.simplify((lin(K_norm_p) - K_norm) / eps)
    compat_direct_p = sp.simplify((N0 + eps * n0) / Ptarget - 3 * (S + eps * z2) ** 2 / (T + eps * z4))
    compat_direct = compat_direct_p.subs(eps, 0)
    d_compat_direct = sp.simplify((lin(compat_direct_p) - compat_direct) / eps)
    d_compat = sp.simplify((lin(K_norm_p - K_one_pole_p) - (K_norm - K_one_pole)) / eps)
    D0target = sp.symbols("D0target", nonzero=True)
    Ptarget_transport_p = sp.simplify((N0 + eps * n0) / D0target)
    K_norm_transport_p = sp.solve(
        sp.Eq((N0 + eps * n0) / (K - B0 - (Z0slot + eps * z0)), Ptarget_transport_p),
        K,
    )[0]
    compat_transport_p = sp.simplify(K_norm_transport_p - K_one_pole_p)
    compat_transport = compat_transport_p.subs(eps, 0)
    d_compat_transport = sp.simplify((lin(compat_transport_p) - compat_transport) / eps)
    assert_zero("compatibility surface after eliminating K", (K_norm_p - K_one_pole_p) - compat_direct_p)
    assert_zero("one-pole K shift", dK_one_pole - (z0 + 6 * S * z2 / T - 3 * S**2 * z4 / T**2))
    assert_zero("normalization K shift", dK_norm - (z0 + n0 / Ptarget))
    assert_zero("compatibility shift from competing K surfaces", d_compat - (dK_norm - dK_one_pole))
    assert_zero("compatibility first variation from eliminated surface", d_compat - d_compat_direct)
    assert_zero("compatibility first variation", d_compat_direct - (n0 / Ptarget - 6 * S * z2 / T + 3 * S**2 * z4 / T**2))
    assert_zero("transported-target normalization K surface", K_norm_transport_p - (B0 + Z0slot + eps * z0 + D0target))
    assert_zero(
        "transported-target compatibility surface",
        compat_transport_p - (D0target - 3 * (S + eps * z2) ** 2 / (T + eps * z4)),
    )
    assert_zero(
        "transported-target compatibility first variation",
        d_compat_transport - (-6 * S * z2 / T + 3 * S**2 * z4 / T**2),
    )
    assert_nonzero(
        "mutated compatibility transport should fail",
        d_compat_direct - (n0 / Ptarget - 6 * S * z2 / T - 3 * S**2 * z4 / T**2),
    )
    assert_nonzero(
        "mutated transported-target compatibility should fail",
        d_compat_transport - (-6 * S * z2 / T - 3 * S**2 * z4 / T**2),
    )

    # Specialize to the constant-prefactor branch
    N2_const = sp.simplify(2 * D2 * N0 / D0)
    N4_const = sp.simplify((2 * D0 * (D2 * N2_const + D4 * N0) - 3 * D2**2 * N0) / D0**2)
    dP2_const = sp.simplify(dP2.subs(N2, N2_const))
    dP4_const = sp.simplify(dP4.subs({N2: N2_const, N4: N4_const}))
    P2_base = (D0 * N2 - 2 * D2 * N0) / D0**2
    P4_base = (D0**2 * N4 - 2 * D0 * (D2 * N2 + D4 * N0) + 3 * D2**2 * N0) / D0**3
    assert_zero("constant-prefactor P2 factorization", P2_base - (N2 - N2_const) / D0)
    assert_zero("constant-prefactor P4 factorization", P4_base.subs(N2, N2_const) - (N4 - N4_const) / D0)

    # Weak-axisymmetric grouped signature
    ea = sp.symbols("ea")
    lam20 = real_y20_square_ratio(0)
    lam21 = real_y20_square_ratio(1)
    lam22 = real_y20_square_ratio(2)
    assert_zero("Y20 overlap lane 20", lam20 - 1)
    assert_zero("Y20 overlap lane 21", lam21 - sp.Rational(1, 2))
    assert_zero("Y20 overlap lane 22", lam22 + 1)

    x0, x1 = sp.symbols("x0 x1")
    x20 = x0 + ea * lam20 * x1
    x21 = x0 + ea * lam21 * x1
    x22 = x0 + ea * lam22 * x1

    xbar = sp.simplify((x20 + 2 * x21 + 2 * x22) / 5)
    ax = sp.simplify((2 * x20 - x21 - x22) / 10)
    bx = sp.simplify((x21 - x22) / 2)
    assert_zero("weak-axisymmetric grouped trace", xbar - x0)
    assert_zero("weak-axisymmetric b=3a", bx - 3 * ax)

    # Static prefactor slope if only the projected-Maxwell static bundle is anisotropic:
    lam = sp.symbols("lam")
    za, na = sp.symbols("za na")  # weak-axisymmetric static slopes in D0 and N0
    P_lane = sp.expand((N0 + ea * lam * na) / (D0 - ea * lam * za))
    P_lane_lin = sp.expand(sp.series(P_lane, ea, 0, 2).removeO())
    Xi1_proj = sp.simplify((sp.diff(P_lane_lin, ea).subs(ea, 0) / lam) / (N0 / D0))
    assert_zero("static Xi1 slope", Xi1_proj - (na/N0 + za/D0))

    # Same for u2 in each lane if only Z0,Z2 carry the anisotropy:
    z2a = sp.symbols("z2a")
    u2_lane = sp.expand(-(D2 - ea * lam * z2a) / (D0 - ea * lam * za))
    u2_lane_lin = sp.expand(sp.series(u2_lane, ea, 0, 2).removeO())
    u2_slope = sp.simplify((sp.diff(u2_lane_lin, ea).subs(ea, 0) / lam))
    assert_zero("u2 projected-Maxwell slope", u2_slope - (D0*z2a - D2*za) / D0**2)

    lines = []
    lines.append(banner("1) Base isotropic full-bundle quantities"))
    lines.append("The moving-throat grouped-bundle backbone uses\n")
    lines.append(line("D0 =", "K - B0 - Z0"))
    lines.append(line("D2 =", "-(M + B2 + Z2)"))
    lines.append(line("D4 =", "-(B4 + Z4)"))
    lines.append(line("P0 =", "N0 / D0"))
    lines.append(line("P2 =", "(D0*N2 - 2*D2*N0)/D0**2"))
    lines.append(line("P4 =", "(D0**2*N4 - 2*D0*(D2*N2 + D4*N0) + 3*D2**2*N0)/D0**3"))
    lines.append(
        "We now perturb only the projected Maxwell / mixed moments,\n"
        "interpreting eps as the projection-first correction scale.\n"
        "For a mouth-anchored kernel, eps ~ ell. For a symmetric interior kernel, eps ~ sigma^2.\n"
    )

    lines.append(banner("2) Exact first-order bridge from projected-Maxwell moments to bundle quantities"))
    lines.append("Projected Maxwell near-throat corrections:\n")
    lines.append(line("D0 ->", "D0 - eps*z0"))
    lines.append(line("D2 ->", "D2 - eps*z2"))
    lines.append(line("D4 ->", "D4 - eps*z4"))
    lines.append(line("N0 ->", "N0 + eps*n0"))
    lines.append(line("N2 ->", "N2 + eps*n2"))
    lines.append(line("N4 ->", "N4 + eps*n4"))
    lines.append("\nInduced normalized-response / prefactor shifts:\n")
    lines.append(line("delta u2 =", du2))
    lines.append(line("delta u4 =", du4))
    lines.append(line("delta P0 =", dP0))
    lines.append(line("delta P2 =", dP2))
    lines.append(line("delta P4 =", dP4))

    lines.append(banner("3) One-pole defect and isotropic compatibility surface"))
    lines.append("One-pole defect:\n")
    lines.append(line("Pole =", "D0*(B4 + Z4) - 3*(M + B2 + Z2)**2"))
    lines.append(line("delta Pole =", d_pole))
    lines.append(
        "So the projected-Maxwell conservative near-throat data enter the one-pole defect only through\n"
        "z0, z2, z4, with the exact linear combination above.\n"
    )
    lines.append("\nCompatibility surface written in the sharp Stage-18 form:\n")
    lines.append(line("Compat =", "N0/Ptarget - 3*S**2/T"))
    lines.append(line("delta Compat =", d_compat))
    lines.append(
        "Important structural result:\n"
        "  z0 drops out of delta Compat after exact elimination of K and D0.\n"
        "So a projected-Maxwell shift of the static conservative denominator D0 changes P0 directly,\n"
        "but does NOT move the isotropic compatibility equation between one-pole and normalization.\n"
    )

    lines.append(banner("4) Constant-prefactor branch: exact linearized transport"))
    lines.append(line("Impose P2 = 0 by setting N2 =", N2_const))
    lines.append(line("Then delta P2 becomes =", dP2_const))
    lines.append(line("Impose also P4 = 0 by setting N4 =", N4_const))
    lines.append(line("Then delta P4 becomes =", dP4_const))
    lines.append(
        "So on the constant-prefactor branch, the projected-Maxwell near-throat data that matter first are\n"
        "  z0, z2, z4, n0, n2, n4,\n"
        "with z0 entering delta P2 and delta P4 only through denominator transport.\n"
    )

    lines.append(banner("5) Weak-axisymmetric grouped signature"))
    lines.append("Use the exact Y20 weak-axisymmetric grouped pattern\n")
    lines.append(line("(lambda20, lambda21, lambda22) =", "(1, 1/2, -1)"))
    lines.append(line("Weighted grouped trace xbar =", xbar))
    lines.append(line("Grouped anomaly a_x =", ax))
    lines.append(line("Grouped anomaly b_x =", bx))
    lines.append(
        "Thus any first weak-axisymmetric projected-Maxwell contribution obeys\n"
        "  xbar = x^(0),   a_x = eps*x^(1)/4,   b_x = 3*eps*x^(1)/4.\n"
        "So the isotropic grouped trace is unchanged at first order, while the anisotropy sits on the line b = 3 a.\n"
    )

    lines.append(banner("6) Direct projected-Maxwell contribution to Xi1 = P1/P0"))
    lines.append(
        "If only the STATIC projected-Maxwell bundle is anisotropic at first order,\n"
        "so that for lane A we have\n"
        "  D0_A = D0 - ea*lambda_A*za,\n"
        "  N0_A = N0 + ea*lambda_A*na,\n"
        "then the prefactor in each lane becomes\n"
    )
    lines.append(line("P_A =", P_lane_lin))
    lines.append(line("Xi1^(proj,static) =", Xi1_proj))
    lines.append(
        "Therefore the projected-Maxwell static near-throat slippage contributes to the weak-axisymmetric bottleneck by\n"
        "  Xi1^(proj,static) = na/N0 + za/D0.\n"
        "Equivalently,\n"
        "  P_A = P0 * (1 + ea*lambda_A*Xi1^(proj,static)).\n"
    )
    lines.append(
        "This is the direct projection-first bridge to the transported prefactor packet:\n"
        "  the grouped trace stays fixed at first order, but Xi1 is shifted linearly.\n"
    )

    lines.append(banner("7) Direct projected-Maxwell contribution to grouped u2 anisotropy"))
    lines.append(
        "If only Z0 and Z2 carry the weak-axisymmetric projected-Maxwell slope,\n"
        "then lanewise\n"
    )
    lines.append(line("u2_A =", u2_lane_lin))
    lines.append(line("u2 slope coefficient =", u2_slope))
    lines.append(
        "So the same projected-Maxwell weak-axisymmetric slippage also feeds the grouped conservative response anisotropy.\n"
    )

    lines.append(banner("8) Reading the result back into the near-throat Maxwell picture"))
    lines.append(
        "The near-throat projected Maxwell work found:\n"
        "  - symmetric interior kernels start at O(sigma^2),\n"
        "  - mouth-anchored kernels start at O(ell).\n"
        "So the bundle formulas above imply:\n"
        "  - interior projection-first corrections feed z_n, n_n only at even order,\n"
        "  - mouth-local projection-first corrections can feed the same bundle data already at first order.\n"
        "\n"
        "That makes the mouth the natural place where projection-first Maxwell can alter\n"
        "  Z2, Z4, N0, N2, N4  (isotropic calibration)\n"
        "and\n"
        "  Xi1 = P1/P0         (weak-axisymmetric prefactor transport).\n"
    )

    print("".join(lines))
    print("STATUS: PASS")


if __name__ == "__main__":
    main()
