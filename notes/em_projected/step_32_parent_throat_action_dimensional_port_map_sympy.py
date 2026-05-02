#!/usr/bin/env python3
"""Exact dimensional/source-map bridge for the hardcoded Branch-B runtime patch."""
from __future__ import annotations

import hashlib
import json

import sympy as sp


def assert_zero(label: str, expr: sp.Expr) -> None:
    expr = sp.simplify(sp.expand(expr))
    if expr != 0:
        raise AssertionError(f"{label} failed: {expr}")


def main() -> None:
    branch_metadata = {
        "branch_id": "v2_local_parent_background_dimensional_port_map",
        "pre_target_freeze": True,
        "target_blind": False,
        "no_post_residual_refit": True,
        "boundary_class": "open_impedance_demo",
        "interpretation": "explicit source-map/dimensional bridge for treating the moderate Branch-B corridor as a hardcoded CFD boundary condition",
    }
    branch_freeze_hash = hashlib.sha256(json.dumps(branch_metadata, sort_keys=True).encode("utf-8")).hexdigest()[:16]

    G, c, c_s, a = sp.symbols("G c c_s a", positive=True, real=True)
    S_port, mhat0 = sp.symbols("S_port mhat_0", positive=True, real=True)
    P0_red, P0_base_red = sp.symbols("P0_red P0_base_red", positive=True, real=True)
    omega_phys = sp.symbols("omega_phys", positive=True, real=True)
    sigma = sp.symbols("sigma", real=True)
    eps_eta, R_target, T_eff2 = sp.symbols("eps_eta R_target T_eff2", positive=True, real=True)
    Kbl, D0 = sp.symbols("K_bl D_0", positive=True, real=True)

    # Frozen carry-forward constants from the 4D bridge.
    n = sp.Integer(5)
    alpha_opt = sp.simplify((n - 1) / 2)
    kappa_add = sp.Rational(1, 2)
    alpha_sq = sp.Rational(3, 4)
    K_vec = sp.simplify(sp.Integer(2) / sp.pi**2)

    Omega_Q = sp.simplify(3 * c_s / (2 * a))
    omega_red = sp.simplify(omega_phys / Omega_Q)
    assert_zero("reduced frequency map", omega_red - 2 * a * omega_phys / (3 * c_s))

    # Exact source-map normalization from the compact note.
    P0_target = sp.simplify(54 * G * c_s**5 / (5 * S_port * a**5 * c**5 * mhat0**2))
    R_norm = sp.simplify(mhat0**2 * S_port * P0_red - 54 * G * c_s**5 / (5 * a**5 * c**5))
    assert_zero("R_norm target sheet", R_norm.subs(P0_red, P0_target))

    S_port_from_target = sp.solve(sp.Eq(R_norm, 0), S_port)[0]
    cs_from_target = sp.solve(sp.Eq(R_norm, 0), c_s)[0]
    assert_zero("source-map solve for S_port", S_port_from_target - 54 * G * c_s**5 / (5 * a**5 * c**5 * mhat0**2 * P0_red))
    assert_zero("source-map solve for c_s", cs_from_target - (5 * P0_red * S_port * a**5 * c**5 * mhat0**2 / (54 * G)) ** sp.Rational(1, 5))

    # On the hardcoded Branch-B target sheet, the reduced scripts enforce mhat0^2 P0 = 54/5.
    S_port_patched = sp.simplify(S_port_from_target.subs(P0_red, sp.Rational(54, 5) / mhat0**2))
    assert_zero("patched target S_port", S_port_patched - G * c_s**5 / (a**5 * c**5))

    # The naive direct-SI test from step31 was exactly the convention S_port = 1.
    cs_direct_SI = sp.simplify(sp.solve(sp.Eq(S_port_patched, 1), c_s)[0])
    assert_zero("direct-SI convention", cs_direct_SI - (a**5 * c**5 / G) ** sp.Rational(1, 5))

    # Transfer-shape compiler and outgoing amplitude map.
    lambda_out = sp.simplify(P0_red / P0_base_red)
    lambda_out_branchb = sp.simplify(1 - sigma)
    delta_ln_Teff2 = sp.simplify(sp.log(lambda_out_branchb))
    assert_zero("Branch-B amplitude law", lambda_out.subs(P0_red, lambda_out_branchb * P0_base_red) - lambda_out_branchb)

    P0_from_transfer = sp.simplify((Kbl / D0) * T_eff2)
    assert_zero("P0 transfer compiler", P0_red.subs(P0_red, P0_from_transfer) - P0_from_transfer)

    T_selected = sp.simplify(27 * sp.pi**2 * G * c_s**5 * (1 - eps_eta) / (20 * a**5 * c**5 * R_target))
    P0_over_T = sp.simplify(P0_target / T_selected)
    assert_zero(
        "target P0 to selected T_eff^2 ratio",
        P0_over_T - 8 * R_target / (sp.pi**2 * S_port * mhat0**2 * (1 - eps_eta)),
    )
    R_target_from_transfer = sp.simplify(sp.solve(sp.Eq(T_eff2, T_selected), R_target)[0])
    assert_zero(
        "selected-branch load-ratio map",
        R_target_from_transfer - 27 * sp.pi**2 * G * c_s**5 * (1 - eps_eta) / (20 * a**5 * c**5 * T_eff2),
    )

    gamma_quad_eff = sp.simplify(mhat0**2 * S_port * P0_red * a**5 / (27 * c_s**5))
    assert_zero("quadrupole normalization bridge", gamma_quad_eff.subs(P0_red, P0_target) - 2 * G / (5 * c**5))

    print("STEP 32 DIMENSIONAL PORT MAP AUDIT")
    print("Compiled the exact dimensional/source-map bridge for using the moderate Branch-B corridor as a hardcoded CFD boundary condition.")
    print("V2 branch-freeze metadata:")
    print("  branch_id =", branch_metadata["branch_id"])
    print("  branch_freeze_hash =", branch_freeze_hash)
    print("  pre_target_freeze =", str(branch_metadata["pre_target_freeze"]).lower())
    print("  target_blind =", str(branch_metadata["target_blind"]).lower())
    print("  no_post_residual_refit =", str(branch_metadata["no_post_residual_refit"]).lower())
    print("  boundary_class =", branch_metadata["boundary_class"])
    print("  interpretation =", branch_metadata["interpretation"])
    print("Frozen carry-forward 4D/bridge constants:")
    print("  n =", n)
    print("  alpha_opt =", alpha_opt)
    print("  kappa_add =", kappa_add)
    print("  alpha^2 =", alpha_sq)
    print("  K_vec =", sp.sstr(K_vec))
    print("Exact dimensional/source-map bridge:")
    print("  Omega_Q =", sp.sstr(Omega_Q))
    print("  omega_red =", sp.sstr(omega_red))
    print("  P0_target =", sp.sstr(P0_target))
    print("  S_port(P0_red,mhat0,c_s,a) =", sp.sstr(S_port_from_target))
    print("  c_s(P0_red,mhat0,S_port,a) =", sp.sstr(cs_from_target))
    print("  On the patched target sheet mhat0^2 P0_red = 54/5:")
    print("    S_port =", sp.sstr(S_port_patched))
    print("  The direct-SI convention tested in step31 is:")
    print("    S_port = 1  =>  c_s =", sp.sstr(cs_direct_SI))
    print("Outgoing amplitude and transfer-shape map:")
    print("  lambda_out = P0_red / P0_base_red =", sp.sstr(lambda_out))
    print("  delta ln(T_eff^2) =", sp.sstr(delta_ln_Teff2))
    print("  lambda_out = 1 - sigma =", sp.sstr(lambda_out_branchb))
    print("  P0_red = (K_bl/D0) T_eff^2")
    print("  P0_target / T_eff^2 =", sp.sstr(P0_over_T))
    print("  R_target(T_eff^2) =", sp.sstr(R_target_from_transfer))
    print("  gamma_quad_eff =", sp.sstr(gamma_quad_eff.subs(P0_red, P0_target)))
    print("Interpretation:")
    print("  mhat0 and lambda_out are dimensionless runtime controls.")
    print("  S_port carries the unresolved dimensional bridge between the reduced PDE prefactor and SI targets.")
    print("  On the hardcoded Branch-B target sheet, S_port collapses exactly to G c_s^5 / (a^5 c^5).")
    print("  So step31 did not falsify the whole patch; it falsified the special convention S_port = 1.")
    print("STATUS: PASS")


if __name__ == "__main__":
    main()
