#!/usr/bin/env python3
"""Stage95 branch-B parameter burden for the moderate sigma corridor."""
from __future__ import annotations

import hashlib
import json

import sympy as sp

from step_29_parent_throat_action_moderate_branchb_sector_sympy import export_step29_moderate_branchb_patch


def assert_zero(label: str, expr: sp.Expr) -> None:
    expr = sp.simplify(sp.expand(expr))
    if expr != 0:
        raise AssertionError(f"{label} failed: {expr}")


def assert_close(label: str, actual: float, expected: float, tol: float) -> None:
    if abs(actual - expected) > tol:
        raise AssertionError(f"{label} failed: {actual} vs {expected} (tol={tol})")


def assert_nonzero(label: str, expr: sp.Expr) -> None:
    expr = sp.simplify(sp.expand(expr))
    if expr == 0:
        raise AssertionError(f"{label} unexpectedly vanished")


def main() -> None:
    branch_metadata = {
        "branch_id": "v2_local_parent_background_branchb_parameter_burden",
        "pre_target_freeze": True,
        "target_blind": True,
        "no_post_residual_refit": True,
        "boundary_class": "open_impedance_demo",
        "interpretation": "compare the moderate Stage95 branch-B outlet parameters against the forbidden natural-source Robin burden",
    }
    branch_freeze_hash = hashlib.sha256(json.dumps(branch_metadata, sort_keys=True).encode("utf-8")).hexdigest()[:16]

    sigma, rho, kappa, gamma = sp.symbols("sigma rho kappa gamma", real=True)
    eps = sp.symbols("eps", real=True, nonzero=True)

    # Exact stage95 branch-B canonical-even solve.
    rho_branch_b = 4 * sigma
    kappa_branch_b = sp.Rational(1, 3)
    gamma_canonical = sp.Rational(1, 9)
    chi_branch_b = sp.simplify((1 - 9 * sigma * gamma) / (1 - sigma))

    assert_zero("stage95 branch-B rho law", rho_branch_b - 4 * sigma)
    assert_zero("stage95 branch-B kappa law", kappa_branch_b - sp.Rational(1, 3))
    assert_zero("stage95 branch-B canonical chi_Q", chi_branch_b.subs(gamma, gamma_canonical) - 1)
    assert_nonzero("stage95 branch-B rho law detects mutated coefficient", rho_branch_b - (4 + eps) * sigma)
    assert_nonzero("stage95 branch-B chi_Q detects mutated gamma", chi_branch_b.subs(gamma, gamma_canonical + eps) - 1)

    # Step26 natural-source minimal Robin burdens.
    rho_nat_scale009 = -113614.01985490319
    rho_nat_scale0092 = -137733.1209385474

    step29_patch = export_step29_moderate_branchb_patch()
    point_sigma19 = {
        "label": "moderate same-scale point",
        "sigma": float(step29_patch["same_scale_lambda20"]["sigma"]),
        "lambda_out": float(step29_patch["same_scale_lambda20"]["lambda_out"]),
        "delta_ln_t2": float(step29_patch["same_scale_lambda20"]["delta_ln_t2"]),
    }
    point_sigma49 = {
        "label": "stronger same-scale point",
        "sigma": float(step29_patch["same_scale_lambda50"]["sigma"]),
        "lambda_out": float(step29_patch["same_scale_lambda50"]["lambda_out"]),
        "delta_ln_t2": float(step29_patch["same_scale_lambda50"]["delta_ln_t2"]),
    }

    def branch_b_rho(point: dict[str, float]) -> float:
        return 4.0 * point["sigma"]

    rho_sigma19 = branch_b_rho(point_sigma19)
    rho_sigma49 = branch_b_rho(point_sigma49)
    burden_ratio_19 = abs(rho_nat_scale009) / abs(rho_sigma19)
    burden_ratio_49 = abs(rho_nat_scale0092) / abs(rho_sigma49)
    sigma19_sign_flip_residual = -4.0 * point_sigma19["sigma"] - rho_sigma19
    sigma49_sign_flip_residual = -4.0 * point_sigma49["sigma"] - rho_sigma49

    assert_close("sigma=-19 branch-B rho", rho_sigma19, -76.0, 1e-12)
    assert_close("sigma=-49 branch-B rho", rho_sigma49, -196.0, 1e-12)
    assert_close("sigma=-19 burden reduction factor", burden_ratio_19, 1494.9213138803052, 1e-6)
    assert_close("sigma=-49 burden reduction factor", burden_ratio_49, 702.7200047885071, 1e-6)
    if abs(sigma19_sign_flip_residual) < 1.0 or abs(sigma49_sign_flip_residual) < 1.0:
        raise AssertionError("branch-B rho sign-flip mutation did not move the moderate corridor enough")

    print("STEP 30 BRANCH-B PARAMETER BURDEN AUDIT")
    print("Mapped the moderate branch-B corridor to exact Stage95 outlet parameters and compared that burden to the forbidden natural-source Robin loads.")
    print("V2 branch-freeze metadata:")
    print("  branch_id =", branch_metadata["branch_id"])
    print("  branch_freeze_hash =", branch_freeze_hash)
    print("  pre_target_freeze =", str(branch_metadata["pre_target_freeze"]).lower())
    print("  target_blind =", str(branch_metadata["target_blind"]).lower())
    print("  no_post_residual_refit =", str(branch_metadata["no_post_residual_refit"]).lower())
    print("  boundary_class =", branch_metadata["boundary_class"])
    print("  interpretation =", branch_metadata["interpretation"])
    print("Exact Stage95 branch-B laws:")
    print("  rho =", sp.sstr(rho_branch_b))
    print("  kappa =", sp.sstr(kappa_branch_b))
    print("  canonical gamma =", sp.sstr(gamma_canonical))
    print("  canonical chi_Q =", sp.sstr(chi_branch_b.subs(gamma, gamma_canonical)))
    print("  rho sign-flip mutation residuals =", (sigma19_sign_flip_residual, sigma49_sign_flip_residual))
    print("Moderate corridor points:")
    print(" ", point_sigma19)
    print("    branch-B rho =", rho_sigma19)
    print("    natural-source Robin burden ratio =", burden_ratio_19)
    print(" ", point_sigma49)
    print("    branch-B rho =", rho_sigma49)
    print("    natural-source Robin burden ratio =", burden_ratio_49)
    print("Interpretation:")
    print("  The admissible branch-B corridor is still nontrivial, but its outlet burden is orders of magnitude smaller than the forbidden natural-source Robin loads from step26.")
    print("  The moderate sigma = -19 point corresponds to rho = -76, not rho ~ -1.136e5.")
    print("  The stronger sigma = -49 point corresponds to rho = -196, not rho ~ -1.377e5.")
    print("  So the branch-B escape hatch is not just algebraically distinct from the N_Q story; it is quantitatively much cheaper at the outlet-parameter level.")
    print("STATUS: PASS")


if __name__ == "__main__":
    main()
