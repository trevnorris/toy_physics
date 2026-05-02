#!/usr/bin/env python3
"""Master-note audit for step_13_parent_throat_action_master_notes.md."""
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


def boundary_value(expr: sp.Expr, var: sp.Symbol) -> sp.Expr:
    return sp.simplify(sp.limit(expr, var, sp.oo) - sp.limit(expr, var, -sp.oo))


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
    w_ibp = sp.symbols("w_ibp", real=True)
    Aw = sp.Function("A")(w_ibp)
    etaw = sp.Function("eta")(w_ibp)
    quad_boundary_density = -Aw * etaw**2 / 2
    quad_bulk_density = sp.diff(Aw, w_ibp) * etaw**2 / 2
    assert_zero(
        "generic quadratic IBP identity",
        -Aw * etaw * sp.diff(etaw, w_ibp) - (sp.diff(quad_boundary_density, w_ibp) + quad_bulk_density),
    )
    assert_nonzero(
        "mutated quadratic IBP sign should fail",
        -Aw * etaw * sp.diff(etaw, w_ibp) - (sp.diff(-quad_boundary_density, w_ibp) + quad_bulk_density),
    )
    A_concrete = sp.exp(-w_ibp**2)
    eta_concrete = sp.exp(-w_ibp**2 / 2)
    quad_boundary_concrete = boundary_value((-A_concrete * eta_concrete**2 / 2), w_ibp)
    quad_cross_concrete = sp.integrate(
        -A_concrete * eta_concrete * sp.diff(eta_concrete, w_ibp),
        (w_ibp, -sp.oo, sp.oo),
    )
    quad_bulk_concrete = sp.integrate(
        sp.diff(A_concrete, w_ibp) * eta_concrete**2 / 2,
        (w_ibp, -sp.oo, sp.oo),
    )
    assert_zero("concrete quadratic IBP boundary discharge", quad_boundary_concrete)
    assert_zero(
        "concrete quadratic IBP with boundary",
        quad_cross_concrete - (quad_boundary_concrete + quad_bulk_concrete),
    )

    eps = sp.symbols("eps")
    eta, eta_t, eta_w, grad2 = sp.symbols("eta eta_t eta_w grad2")
    R0p = sp.symbols("R0p")
    mu0, Tw0, TO0, U0 = sp.symbols("mu0 Tw0 TO0 U0")
    TwR0, TwRR0, UR0, URR0 = sp.symbols("TwR0 TwRR0 UR0 URR0")
    d_TwR_R0p = sp.symbols("d_TwR_R0p")

    Tw = Tw0 + eps * TwR0 * eta + eps**2 * TwRR0 * eta**2 / 2
    U = U0 + eps * UR0 * eta + eps**2 * URR0 * eta**2 / 2
    Rt = eps * eta_t
    Rw = R0p + eps * eta_w

    L = mu0 * Rt**2 / 2 - Tw * Rw**2 / 2 - TO0 * eps**2 * grad2 / 2 - U
    L2_raw = sp.diff(sp.expand(L), eps, 2).subs(eps, 0) / 2
    K_eta = URR0 - d_TwR_R0p + TwRR0 * R0p**2 / 2
    canonical_L2 = mu0 * eta_t**2 / 2 - Tw0 * eta_w**2 / 2 - TO0 * grad2 / 2 - K_eta * eta**2 / 2

    # The raw cross term must be the one integrated by parts.
    raw_cross_coeff = sp.diff(sp.diff(L2_raw, eta), eta_w)
    assert_zero("raw eta eta_w cross term", raw_cross_coeff + TwR0 * R0p)
    cross_term = -TwR0 * R0p * eta * eta_w
    cross_after_ibp = d_TwR_R0p * eta**2 / 2
    L2_after_ibp_derived = sp.expand(L2_raw - cross_term + cross_after_ibp)
    assert_zero("quadratic wall action after integration by parts", L2_after_ibp_derived - canonical_L2)
    canonical_L2_mutated = mu0 * eta_t**2 / 2 - Tw0 * eta_w**2 / 2 - TO0 * grad2 / 2 - (
        URR0 + d_TwR_R0p + TwRR0 * R0p**2 / 2
    ) * eta**2 / 2
    assert_nonzero("mutated K_eta sign should fail", L2_after_ibp_derived - canonical_L2_mutated)

    wall_w = sp.symbols("wall_w", real=True)
    beta = sp.Function("beta")(wall_w)
    delta_mu = sp.Function("delta_mu")(wall_w)
    delta_Tw = sp.Function("delta_Tw")(wall_w)
    delta_TO = sp.Function("delta_TO")(wall_w)
    delta_Keta = sp.Function("delta_Keta")(wall_w)
    dM_overlap = sp.Integral(delta_mu * beta**2, (wall_w, -sp.oo, sp.oo))
    dK_overlap = sp.Integral(
        delta_Tw * sp.diff(beta, wall_w)**2 + (delta_Keta + 6 * delta_TO) * beta**2,
        (wall_w, -sp.oo, sp.oo),
    )

    dK, dM = sp.symbols("dK dM")
    B01, B21, B41 = sp.symbols("B01 B21 B41")
    Z01, Z21, Z41 = sp.symbols("Z01 Z21 Z41")
    D01_full = dK - B01 - Z01
    D21_full = -(dM + B21 + Z21)
    D41_full = -(B41 + Z41)
    K1_full = sp.expand(D21_full + D01_full / 9)
    H_even_full = sp.expand(D41_full - sp.Rational(2, 3) * D21_full - D01_full / 27)
    wall_only_specialization = {B01: 0, B21: 0, B41: 0, Z01: 0, Z21: 0, Z41: 0}
    K1_wall = sp.expand(K1_full.subs(wall_only_specialization))
    H_even_wall = sp.expand(H_even_full.subs(wall_only_specialization))
    assert_zero("wall-only K1 specialization", K1_wall - (-dM + dK / 9))
    assert_zero("wall-only H_even specialization", H_even_wall - (sp.Rational(2, 3) * dM - dK / 27))
    assert_zero(
        "wall-only K1 from overlap-generated slots",
        K1_wall.subs({dK: dK_overlap, dM: dM_overlap}) - (-dM_overlap + dK_overlap / 9),
    )
    assert_zero(
        "wall-only H_even from overlap-generated slots",
        H_even_wall.subs({dK: dK_overlap, dM: dM_overlap}) - (sp.Rational(2, 3) * dM_overlap - dK_overlap / 27),
    )
    wall_matrix = sp.Matrix(
        [
            [sp.diff(K1_wall, dK), sp.diff(K1_wall, dM)],
            [sp.diff(H_even_wall, dK), sp.diff(H_even_wall, dM)],
        ]
    )
    assert_zero("wall-only even-gate determinant", wall_matrix.det() - sp.Rational(1, 27))
    wall_even_solve = sp.solve([sp.Eq(K1_wall, 0), sp.Eq(H_even_wall, 0)], [dK, dM], dict=True)[0]
    assert_zero("wall-only dK closed form", wall_even_solve[dK])
    assert_zero("wall-only dM closed form", wall_even_solve[dM])

    lam20 = real_y20_square_ratio(0)
    lam21 = real_y20_square_ratio(1)
    lam22 = real_y20_square_ratio(2)
    assert_zero("Y20 overlap lane 20", lam20 - 1)
    assert_zero("Y20 overlap lane 21", lam21 - sp.Rational(1, 2))
    assert_zero("Y20 overlap lane 22", lam22 + 1)

    x0, eps1 = sp.symbols("x0 eps1")
    x20 = x0 + eps1 * lam20
    x21 = x0 + eps1 * lam21
    x22 = x0 + eps1 * lam22
    xbar, ax, bx = grouped_trace_anomaly(x20, x21, x22)
    assert_zero("grouped trace", xbar - x0)
    assert_zero("grouped line b=3a", bx - 3 * ax)

    print("STEP 13 PARENT THROAT ACTION MASTER AUDIT")
    print("Checked promoted-action quadratic limit, concrete boundary discharge, K_eta formula, grouped signature, and wall-only gate obstruction as the zero-support/zero-mixed specialization of the full even gates.")
    print("STATUS: PASS")


if __name__ == "__main__":
    main()
