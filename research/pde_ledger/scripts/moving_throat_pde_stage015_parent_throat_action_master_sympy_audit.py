#!/usr/bin/env python3
"""Master-note audit for step_13_parent_throat_action_master_notes.md."""
from __future__ import annotations

import sympy as sp


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
    nonzero_boundary_probe = boundary_value(sp.atan(w_ibp), w_ibp)
    assert_nonzero("boundary operator detects nonzero endpoint discharge", nonzero_boundary_probe)
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
    # Second concrete probe with asymmetric A profile so the cross and bulk
    # integrals are individually nonzero. Both integrals are even by parity
    # (A_extra is odd, eta is even, products produce even integrands), so
    # the IBP identity asserts a real cancellation, not 0 = 0 + 0.
    A_concrete_asym = w_ibp * sp.exp(-w_ibp**2)
    eta_concrete_asym = sp.exp(-w_ibp**2 / 2)
    quad_boundary_concrete_asym = boundary_value(
        (-A_concrete_asym * eta_concrete_asym**2 / 2), w_ibp
    )
    quad_cross_concrete_asym = sp.integrate(
        -A_concrete_asym * eta_concrete_asym * sp.diff(eta_concrete_asym, w_ibp),
        (w_ibp, -sp.oo, sp.oo),
    )
    quad_bulk_concrete_asym = sp.integrate(
        sp.diff(A_concrete_asym, w_ibp) * eta_concrete_asym**2 / 2,
        (w_ibp, -sp.oo, sp.oo),
    )
    assert_nonzero("asymmetric concrete IBP cross is nontrivial", quad_cross_concrete_asym)
    assert_nonzero("asymmetric concrete IBP bulk is nontrivial", quad_bulk_concrete_asym)
    assert_zero("asymmetric concrete quadratic IBP boundary discharge", quad_boundary_concrete_asym)
    assert_zero(
        "asymmetric concrete quadratic IBP with boundary",
        quad_cross_concrete_asym - (quad_boundary_concrete_asym + quad_bulk_concrete_asym),
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

    print("STEP 13 PARENT THROAT ACTION MASTER AUDIT")
    print("Checked promoted-action quadratic limit, concrete boundary discharge, and K_eta formula.")
    print("Boundary operator nonzero sanity check = PASS")
    print("STATUS: PASS")


if __name__ == "__main__":
    main()
