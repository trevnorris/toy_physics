#!/usr/bin/env python3
"""Master-note audit for step_08_projected_maxwell_push_bundle_master_notes.md."""
from __future__ import annotations

import sympy as sp
from sympy.physics.wigner import gaunt


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


def grouped_trace_anomaly(x20: sp.Expr, x21: sp.Expr, x22: sp.Expr) -> tuple[sp.Expr, sp.Expr, sp.Expr]:
    xbar = sp.simplify((x20 + 2 * x21 + 2 * x22) / 5)
    ax = sp.simplify((2 * x20 - x21 - x22) / 10)
    bx = sp.simplify((x21 - x22) / 2)
    return xbar, ax, bx


def main() -> None:
    eps = sp.symbols("eps")
    D0, D2, D4, N0, N2, N4 = sp.symbols("D0 D2 D4 N0 N2 N4", nonzero=True)
    z0, z2, z4, n0, n2, n4 = sp.symbols("z0 z2 z4 n0 n2 n4")

    D0p = D0 - eps * z0
    D2p = D2 - eps * z2
    D4p = D4 - eps * z4
    N0p = N0 + eps * n0
    N2p = N2 + eps * n2
    N4p = N4 + eps * n4

    P0p = N0p / D0p
    P2p = (D0p * N2p - 2 * D2p * N0p) / D0p**2
    P4p = (D0p**2 * N4p - 2 * D0p * (D2p * N2p + D4p * N0p) + 3 * D2p**2 * N0p) / D0p**3

    dP0 = sp.diff(P0p, eps).subs(eps, 0)
    dP2 = sp.diff(P2p, eps).subs(eps, 0)
    dP4 = sp.diff(P4p, eps).subs(eps, 0)

    assert_zero("delta P0", dP0 - (n0 / D0 + N0 * z0 / D0**2))
    assert_zero(
        "delta P2",
        dP2 - (
            n2 / D0
            + N2 * z0 / D0**2
            + 2 * N0 * z2 / D0**2
            - 2 * D2 * n0 / D0**2
            - 4 * D2 * N0 * z0 / D0**3
        ),
    )
    assert_zero(
        "delta P4",
        dP4 - (
            n4 / D0
            + N4 * z0 / D0**2
            + 2 * N2 * z2 / D0**2
            - 2 * D2 * n2 / D0**2
            + 2 * N0 * z4 / D0**2
            - 2 * D4 * n0 / D0**2
            - 4 * (D2 * N2 + D4 * N0) * z0 / D0**3
            - 6 * D2 * N0 * z2 / D0**3
            + 3 * D2**2 * n0 / D0**3
            + 9 * D2**2 * N0 * z0 / D0**4
        ),
    )

    K, B0, Z0slot, Ptarget, S, T = sp.symbols("K B0 Z0slot Ptarget S T", nonzero=True)
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
    dK_one_pole = sp.diff(sp.series(K_one_pole_p, eps, 0, 2).removeO(), eps).subs(eps, 0)
    dK_norm = sp.diff(sp.series(K_norm_p, eps, 0, 2).removeO(), eps).subs(eps, 0)
    compat_direct_p = sp.simplify((N0 + eps * n0) / Ptarget - 3 * (S + eps * z2) ** 2 / (T + eps * z4))
    compat_direct = compat_direct_p.subs(eps, 0)
    dcompat_direct = sp.diff(sp.series(compat_direct_p, eps, 0, 2).removeO(), eps).subs(eps, 0)
    compat = K_norm - K_one_pole
    dcompat = sp.simplify(
        (sp.series(K_norm_p - K_one_pole_p, eps, 0, 2).removeO() - compat) / eps
    )
    D0target = sp.symbols("D0target", nonzero=True)
    Ptarget_transport_p = sp.simplify((N0 + eps * n0) / D0target)
    K_norm_transport_p = sp.solve(
        sp.Eq((N0 + eps * n0) / (K - B0 - (Z0slot + eps * z0)), Ptarget_transport_p),
        K,
    )[0]
    compat_transport_p = sp.simplify(K_norm_transport_p - K_one_pole_p)
    compat_transport = compat_transport_p.subs(eps, 0)
    dcompat_transport = sp.diff(sp.series(compat_transport_p, eps, 0, 2).removeO(), eps).subs(eps, 0)
    assert_zero("compatibility surface after eliminating K", (K_norm_p - K_one_pole_p) - compat_direct_p)
    assert_zero("one-pole K shift", dK_one_pole - (z0 + 6 * S * z2 / T - 3 * S**2 * z4 / T**2))
    assert_zero("normalization K shift", dK_norm - (z0 + n0 / Ptarget))
    assert_zero("compatibility first variation from competing K surfaces", dcompat - (dK_norm - dK_one_pole))
    assert_zero("compatibility first variation from eliminated surface", dcompat - dcompat_direct)
    assert_zero("compatibility first variation", dcompat_direct - (n0 / Ptarget - 6 * S * z2 / T + 3 * S**2 * z4 / T**2))
    assert_zero("transported-target normalization K surface", K_norm_transport_p - (B0 + Z0slot + eps * z0 + D0target))
    assert_zero(
        "transported-target compatibility surface",
        compat_transport_p - (D0target - 3 * (S + eps * z2) ** 2 / (T + eps * z4)),
    )
    assert_zero(
        "transported-target compatibility first variation",
        dcompat_transport - (-6 * S * z2 / T + 3 * S**2 * z4 / T**2),
    )
    assert_zero("fixed-target compatibility z0 cancellation", sp.diff(dcompat_direct, z0))
    assert_zero("transported-target compatibility z0 cancellation", sp.diff(dcompat_transport, z0))
    assert_nonzero("normalization K surface still carries z0", sp.diff(dK_norm, z0) - 0)
    assert_nonzero(
        "mutated compatibility transport should fail",
        dcompat_direct - (n0 / Ptarget - 6 * S * z2 / T - 3 * S**2 * z4 / T**2),
    )
    assert_nonzero(
        "mutated transported-target compatibility should fail",
        dcompat_transport - (-6 * S * z2 / T - 3 * S**2 * z4 / T**2),
    )

    x0, x1 = sp.symbols("x0 x1")
    lam20 = real_y20_square_ratio(0)
    lam21 = real_y20_square_ratio(1)
    lam22 = real_y20_square_ratio(2)
    assert_zero("Y20 overlap lane 20", lam20 - 1)
    assert_zero("Y20 overlap lane 21", lam21 - sp.Rational(1, 2))
    assert_zero("Y20 overlap lane 22", lam22 + 1)
    x20 = x0 + eps * lam20 * x1
    x21 = x0 + eps * lam21 * x1
    x22 = x0 + eps * lam22 * x1
    xbar, ax, bx = grouped_trace_anomaly(x20, x21, x22)
    assert_zero("weak-axisymmetric trace", xbar - x0)
    assert_zero("weak-axisymmetric a anomaly", ax - eps * x1 / 4)
    assert_zero("weak-axisymmetric b anomaly", bx - 3 * eps * x1 / 4)
    assert_zero("weak-axisymmetric line b=3a", bx - 3 * ax)

    Q, S2, H, Delta, P, Gw = sp.symbols("Q S2 H Delta P Gw", nonzero=True)
    q1, s1, h1, d1, p1, g1 = sp.symbols("q1 s1 h1 d1 p1 g1")
    D0sym, N0sym = sp.symbols("D0sym N0sym", nonzero=True)

    z0_prim = (Delta * q1 - Q * d1) / Delta**2
    n0_prim = 2 * P * (Delta * p1 - P * d1) / Delta**3
    Xi_static = sp.simplify(n0_prim / N0sym + z0_prim / D0sym)

    assert_zero(
        "primitive static Xi formula",
        Xi_static.subs(N0sym, P**2 / Delta**2)
        - (2*p1/P - 2*d1/Delta + (Delta*q1 - Q*d1)/(D0sym*Delta**2)),
    )

    print("STEP 08 PROJECTED MAXWELL PUSH MASTER AUDIT")
    print("Checked bundle perturbation slots, z0 cancellation from compatibility, grouped signature, and primitive static dependencies.")
    print("Target-transport z0 cancellation guards = PASS")
    print("STATUS: PASS")


if __name__ == "__main__":
    main()
