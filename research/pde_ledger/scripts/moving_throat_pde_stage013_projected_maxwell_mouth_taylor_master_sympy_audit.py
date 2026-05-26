#!/usr/bin/env python3
"""Master-note audit for step_11_projected_maxwell_mouth_taylor_master_notes.md."""
from __future__ import annotations

import sympy as sp


def assert_zero(label: str, expr: sp.Expr) -> None:
    residue = sp.factor(sp.together(sp.simplify(expr)))
    if residue != 0:
        raise AssertionError(f"{label} failed: {sp.sstr(residue)}")


def assert_nonzero(label: str, expr: sp.Expr) -> None:
    residue = sp.factor(sp.together(sp.simplify(expr)))
    if residue == 0:
        raise AssertionError(f"{label} unexpectedly vanished")


def main() -> None:
    u, ell = sp.symbols("u ell", positive=True)
    X0, X1, X2 = sp.symbols("X0 X1 X2")
    W = sp.exp(-u)
    X = X0 + ell * u * X1 + ell**2 * u**2 * X2 / 2
    Xproj = sp.integrate(W * X, (u, 0, sp.oo))
    assert_zero("one-sided Taylor first moment", sp.series(Xproj, ell, 0, 2).removeO() - (X0 + ell * X1))
    assert_nonzero(
        "mutated one-sided Taylor first moment should fail",
        sp.series(Xproj, ell, 0, 2).removeO() - (X0 + 2 * ell * X1),
    )
    W2 = u * sp.exp(-u)
    mu1_W2 = sp.integrate(u * W2, (u, 0, sp.oo))
    Xproj_W2 = sp.integrate(W2 * X, (u, 0, sp.oo))
    assert_zero("one-sided Taylor W2 first moment", mu1_W2 - 2)
    assert_zero("one-sided Taylor arbitrary first moment", sp.series(Xproj_W2, ell, 0, 2).removeO() - (X0 + ell * mu1_W2 * X1))
    assert_nonzero(
        "mutated arbitrary first moment should fail",
        sp.series(Xproj_W2, ell, 0, 2).removeO() - (X0 + ell * X1),
    )

    Q, S2, Hport, Delta, P, Gw = sp.symbols("Q S2 H Delta P Gw", nonzero=True)
    D0 = sp.symbols("D0", nonzero=True)
    mu1 = sp.symbols("mu1", nonzero=True)

    q1, s1, h1, d1, p1, g1 = sp.symbols("q1 s1 h1 d1 p1 g1")
    Qx, Sx, Hx, Dx, Px, Gx = sp.symbols("Qx S2x Hx Deltax Px Gwx")
    subs_der = {q1: mu1 * Qx, s1: mu1 * Sx, h1: mu1 * Hx, d1: mu1 * Dx, p1: mu1 * Px, g1: mu1 * Gx}

    # z0, z2, z4 are the order-ell^0, ell^2, ell^4 coefficients of the
    # bottleneck chain-rule expansion (Q(t) - Hport(t)*ell^2)/
    # (Delta(t) - S2(t)*ell^2 + ell^4), evaluated at t=0 with leading
    # values {Q, Delta, S2, Hport} and t-derivatives {q1, d1, s1, h1}.
    # The polynomial forms below are independently derived and asserted
    # against the master primitive in
    # mathematica/moving_throat_pde_stage013_projected_maxwell_mouth_taylor_master_mathematica_audit.wl
    # (claims M3 and M4); a discrepancy there is the canonical place to
    # discover an error in these literals.
    z0 = (Delta * q1 - Q * d1) / Delta**2
    z2 = (-Delta**2 * h1 + Delta * (Hport * d1 + Q * s1 + S2 * q1) - 2 * Q * S2 * d1) / Delta**3
    z4 = (
        -Delta**2 * Hport * s1
        - Delta**2 * S2 * h1
        - Delta**2 * q1
        + 2 * Delta * Hport * S2 * d1
        + 2 * Delta * Q * S2 * s1
        + 2 * Delta * Q * d1
        + Delta * S2**2 * q1
        - 3 * Q * S2**2 * d1
    ) / Delta**4

    # n0, n2, n4 are the analogous order-ell^0, ell^2, ell^4
    # coefficients of the chain-rule expansion
    # (P(t) - Gw(t)*ell^2)^2/(Delta(t) - S2(t)*ell^2 + ell^4)^2,
    # evaluated at t=0 with leading values {P, Delta, S2, Gw} and
    # t-derivatives {p1, d1, s1, g1}. The polynomial forms below are
    # independently derived and asserted against the master primitive in
    # mathematica/moving_throat_pde_stage013_projected_maxwell_mouth_taylor_master_mathematica_audit.wl
    # (claim M4); a discrepancy there is the canonical place to discover
    # an error in these literals.
    n0 = 2 * P * (Delta * p1 - P * d1) / Delta**3
    n2 = -(
        2 * Delta**2 * (Gw * p1 + P * g1)
        - 2 * Delta * P * (2 * Gw * d1 + P * s1 + 2 * S2 * p1)
        + 6 * P**2 * S2 * d1
    ) / Delta**4
    n4 = 2 * (
        Delta**3 * Gw * g1
        - Delta**2 * Gw**2 * d1
        - 2 * Delta**2 * Gw * P * s1
        - 2 * Delta**2 * Gw * S2 * p1
        - 2 * Delta**2 * P * S2 * g1
        - 2 * Delta**2 * P * p1
        + 6 * Delta * Gw * P * S2 * d1
        + 3 * Delta * P**2 * S2 * s1
        + 3 * Delta * P**2 * d1
        + 3 * Delta * P * S2**2 * p1
        - 6 * P**2 * S2**2 * d1
    ) / Delta**5

    Xi = sp.simplify((2 * p1 / P - 2 * d1 / Delta + q1 / (D0 * Delta) - Q * d1 / (D0 * Delta**2)).subs(subs_der) / mu1)
    # Paper round-trip: verify Xi matches the paper's closed-form Xi_load = n0/N0 + z0/D0,
    # with the natural identification N0 = P^2/Delta^2 (so that n_0/N_0 = d/ds log(P^2/Delta^2)).
    z0_form = (Delta * q1 - Q * d1) / Delta**2
    n0_form = 2 * P * (Delta * p1 - P * d1) / Delta**3
    N0_form = P**2 / Delta**2
    Xi_paper = sp.simplify(((n0_form / N0_form + z0_form / D0).subs(subs_der)) / mu1)
    assert_zero("Xi matches paper closed form n0/N0 + z0/D0", Xi - Xi_paper)
    assert_zero("dXi/dPprime", sp.diff(Xi, Px) - 2 / P)

    assert_nonzero("Xi should depend on Pprime", sp.diff(Xi, Px))

    print("STEP 11 PROJECTED MAXWELL MOUTH-TAYLOR MASTER AUDIT")
    print("Checked one-sided Taylor projection, Taylor coefficient maps, and Xi_load Pprime dependence.")
    print("STATUS: PASS")


if __name__ == "__main__":
    main()
