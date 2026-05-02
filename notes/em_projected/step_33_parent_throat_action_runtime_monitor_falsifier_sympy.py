#!/usr/bin/env python3
"""Runtime monitor and hard-fail compiler for a Branch-B CFD analog test."""
from __future__ import annotations

import hashlib
import json

import numpy as np
import sympy as sp

from cfd_runtime_monitor_postprocess import analyze_snapshot
from step_34_parent_throat_action_cfd_runtime_postprocessor_sympy import (
    make_periodic_consistency_snapshot,
    make_radial_snapshot,
)


def assert_zero(label: str, expr: sp.Expr) -> None:
    expr = sp.simplify(sp.expand(expr))
    if expr != 0:
        raise AssertionError(f"{label} failed: {expr}")


def assert_nonzero(label: str, expr: sp.Expr) -> None:
    expr = sp.simplify(sp.expand(expr))
    if expr == 0:
        raise AssertionError(f"{label} was unexpectedly zero")


def assert_small(label: str, value: float, tol: float) -> None:
    if not np.isfinite(value) or abs(value) > tol:
        raise AssertionError(f"{label} failed: {value} > {tol}")


def assert_large(label: str, value: float, floor: float) -> None:
    if not np.isfinite(value) or abs(value) < floor:
        raise AssertionError(f"{label} failed: {value} < {floor}")


def main() -> None:
    branch_metadata = {
        "branch_id": "v2_local_parent_background_runtime_monitor_falsifier",
        "pre_target_freeze": True,
        "target_blind": False,
        "no_post_residual_refit": True,
        "boundary_class": "open_impedance_demo",
        "interpretation": "exact CFD monitor identities plus hard falsifiers for the hardcoded Branch-B gravity analog test",
    }
    # Exact brane projection / continuity variables.
    rho_brane, dt_rho, divJ, S_rho = sp.symbols("rho_brane dtrho divJ S_rho", real=True)
    J_dot_grad_rho = sp.symbols("JdotGradRho", real=True)
    div_v = sp.symbols("div_v", real=True)
    lap_phi = sp.symbols("lap_phi", real=True)
    grad_rho_dot_v = sp.symbols("gradRhoDotV", real=True)
    rho0, phi, lap_phi_ext, mu = sp.symbols("rho_0 phi lap_phi_ext mu", positive=True, real=True)
    c_probe, G, M_eff, b = sp.symbols("c_probe G M_eff b", positive=True, real=True)
    N_probe, Phi_eff = sp.symbols("N_probe Phi_eff", real=True)
    Theta_tail, c, c_s = sp.symbols("Theta_tail c c_s", positive=True, real=True)
    r, A = sp.symbols("r A", positive=True, real=True)
    n = sp.Integer(5)

    # Exact projected continuity and divergence identity.
    R_cont = sp.simplify(dt_rho + divJ - S_rho)
    div_v_expected = sp.simplify((S_rho - dt_rho) / rho_brane - J_dot_grad_rho / rho_brane**2)
    R_div = sp.simplify(div_v - div_v_expected)
    assert_zero("continuity residual identity at the exact source", R_cont.subs(divJ, S_rho - dt_rho))
    assert_zero("divergence identity", R_div.subs(div_v, div_v_expected))

    # Exact longitudinal identity and quasi-static Poisson candidate.
    R_pois_exact = sp.simplify(rho_brane * lap_phi - S_rho + dt_rho + grad_rho_dot_v)
    R_pois_lin = sp.simplify(rho0 * lap_phi - S_rho)
    assert_zero("exact longitudinal identity", R_pois_exact.subs(lap_phi, (S_rho - dt_rho - grad_rho_dot_v) / rho_brane))
    assert_zero("linearized quasi-static Poisson identity", R_pois_lin.subs(lap_phi, S_rho / rho0))

    # Massive-scalar / Yukawa exterior diagnostic.
    phi_newton = -A / r
    phi_yukawa = -A * sp.exp(-mu * r) / r
    phi_power = -A / r**2
    phi_log = A * sp.log(r)
    radial_laplacian = lambda expr: sp.simplify(sp.diff(r**2 * sp.diff(expr, r), r) / r**2)
    Q_r_newton = sp.simplify(4 * sp.pi * r**2 * sp.diff(phi_newton, r))
    Q_r_yukawa = sp.simplify(4 * sp.pi * r**2 * sp.diff(phi_yukawa, r))
    Q_r_power = sp.simplify(4 * sp.pi * r**2 * sp.diff(phi_power, r))
    Q_r_log = sp.simplify(4 * sp.pi * r**2 * sp.diff(phi_log, r))
    mu_eff2_newton = sp.simplify(radial_laplacian(phi_newton) / phi_newton)
    mu_eff2_yukawa = sp.simplify(radial_laplacian(phi_yukawa) / phi_yukawa)
    mu_eff2_power = sp.simplify(radial_laplacian(phi_power) / phi_power)
    mu_eff2_log = sp.simplify(radial_laplacian(phi_log) / phi_log)
    assert_zero("Newton flux plateau derivative", sp.diff(Q_r_newton, r))
    assert_zero("Newton exterior masslessness", mu_eff2_newton)
    assert_zero("Yukawa mass signature", mu_eff2_yukawa - mu**2)
    assert_nonzero("1/r^2 impostor has no flux plateau", sp.diff(Q_r_power, r))
    assert_nonzero("log impostor has no flux plateau", sp.diff(Q_r_log, r))
    assert_nonzero("1/r^2 impostor is not massless Newton exterior", mu_eff2_power)
    assert_nonzero("log impostor is not massless Newton exterior", mu_eff2_log)

    # Optical/redshift coefficient from the n=5 bridge.
    alpha_n = sp.simplify((n - 1) / 2)
    Delta_theta = sp.simplify(2 * alpha_n * G * M_eff / (b * c_probe**2))
    Delta_shapiro = sp.simplify(alpha_n * G * M_eff / c_probe**3)
    alpha_fit = sp.simplify(-c_probe**2 * (N_probe - 1) / Phi_eff)
    R_opt = sp.simplify(alpha_fit - alpha_n)
    N_probe_target = sp.simplify(1 - alpha_n * Phi_eff / c_probe**2)
    alpha_wrong = sp.simplify((n - 2) / 2)
    assert_zero("n=5 optical coefficient", alpha_n - 2)
    assert_zero("n=5 deflection coefficient", Delta_theta - 4 * G * M_eff / (b * c_probe**2))
    assert_zero("n=5 Shapiro coefficient", Delta_shapiro - 2 * G * M_eff / c_probe**3)
    assert_zero("optical coefficient fit identity", R_opt.subs(N_probe, N_probe_target))
    assert_nonzero("wrong optical coefficient mutation", alpha_n - alpha_wrong)

    # Tail scalar from the compact packet.
    R_tail = sp.simplify(Theta_tail * (c / c_s) ** 3 - 1)

    periodic_summary = analyze_snapshot(
        make_periodic_consistency_snapshot(),
        c_probe=1.0,
        bins=18,
        tail_fraction=0.35,
        periodic_xyz=True,
    )
    newton_summary = analyze_snapshot(
        make_radial_snapshot(mu=None),
        c_probe=1.0,
        bins=26,
        tail_fraction=0.3,
        periodic_xyz=False,
        center=(0.0, 0.0, 0.0),
    )
    yukawa_summary = analyze_snapshot(
        make_radial_snapshot(mu=1.4),
        c_probe=1.0,
        bins=26,
        tail_fraction=0.3,
        periodic_xyz=False,
        center=(0.0, 0.0, 0.0),
    )
    bad_optics_snapshot = make_radial_snapshot(mu=None)
    bad_optics_snapshot["N_probe"] = 1.0 - 1.4 * bad_optics_snapshot["Phi_eff"]
    bad_optics_summary = analyze_snapshot(
        bad_optics_snapshot,
        c_probe=1.0,
        bins=26,
        tail_fraction=0.3,
        periodic_xyz=False,
        center=(0.0, 0.0, 0.0),
    )

    assert_small("concrete periodic continuity residual", periodic_summary["max_abs_R_cont"], 1e-11)
    assert_small("concrete periodic exact Poisson residual", periodic_summary["max_abs_R_pois_exact"], 1e-11)
    assert_small("concrete periodic alpha_fit error", periodic_summary["alpha_fit_tail_mean"] - 2.0, 5e-3)
    if not (newton_summary["Q_r_tail_cv"] < yukawa_summary["Q_r_tail_cv"]):
        raise AssertionError("concrete Newton exterior should plateau better than Yukawa exterior")
    assert_small("concrete Newton exterior masslessness", newton_summary["mu_eff2_tail_median"], 0.2)
    assert_large("concrete Yukawa exterior massive-scalar signature", yukawa_summary["mu_eff2_tail_median"], 1.0)
    assert_small("concrete Newton optics coefficient", newton_summary["alpha_fit_tail_mean"] - 2.0, 3e-2)
    assert_large("concrete bad-optics alpha mismatch", bad_optics_summary["alpha_fit_tail_mean"] - 2.0, 0.3)

    branch_freeze_payload = {
        "metadata": branch_metadata,
        "symbolic_exports": {
            "R_cont": sp.sstr(R_cont),
            "R_Pois_exact": sp.sstr(R_pois_exact),
            "Q_r_newton": sp.sstr(Q_r_newton),
            "mu_eff2_yukawa": sp.sstr(mu_eff2_yukawa),
            "alpha_n": sp.sstr(alpha_n),
        },
        "numerical_exports": {
            "periodic_max_abs_R_cont": periodic_summary["max_abs_R_cont"],
            "periodic_max_abs_R_pois_exact": periodic_summary["max_abs_R_pois_exact"],
            "newton_Q_r_tail_cv": newton_summary["Q_r_tail_cv"],
            "newton_mu_eff2_tail_median": newton_summary["mu_eff2_tail_median"],
            "yukawa_Q_r_tail_cv": yukawa_summary["Q_r_tail_cv"],
            "yukawa_mu_eff2_tail_median": yukawa_summary["mu_eff2_tail_median"],
            "bad_optics_alpha_fit_tail_mean": bad_optics_summary["alpha_fit_tail_mean"],
        },
    }
    branch_freeze_hash = hashlib.sha256(json.dumps(branch_freeze_payload, sort_keys=True).encode("utf-8")).hexdigest()[:16]

    print("STEP 33 RUNTIME MONITOR / FALSIFIER AUDIT")
    print("Compiled the exact CFD monitor identities and hard-fail signatures for the hardcoded Branch-B gravity analog test.")
    print("V2 branch-freeze metadata:")
    print("  branch_id =", branch_metadata["branch_id"])
    print("  branch_freeze_hash =", branch_freeze_hash)
    print("  pre_target_freeze =", str(branch_metadata["pre_target_freeze"]).lower())
    print("  target_blind =", str(branch_metadata["target_blind"]).lower())
    print("  no_post_residual_refit =", str(branch_metadata["no_post_residual_refit"]).lower())
    print("  boundary_class =", branch_metadata["boundary_class"])
    print("  interpretation =", branch_metadata["interpretation"])
    print("Primary runtime monitors to record from the CFD grid:")
    print("  rho_brane = ∫ W rho dw")
    print("  J_brane   = ∫ W J_xyz dw")
    print("  S_rho     = -[W J_w] + ∫ W' J_w dw")
    print("  v_brane   = J_brane / rho_brane")
    print("  R_ij      = Pi_ij - rho_brane v_i v_j")
    print("  S_Ji      = -[W J_i J_w/rho] + ∫ W' J_i J_w/rho dw")
    print("Exact monitor residuals:")
    print("  R_cont =", sp.sstr(R_cont))
    print("  R_div  =", sp.sstr(R_div))
    print("  R_Pois_exact =", sp.sstr(R_pois_exact))
    print("  R_Pois_lin   =", sp.sstr(R_pois_lin))
    print("Exterior gravity-like field diagnostics:")
    print("  Newton plateau  Q_r =", sp.sstr(Q_r_newton))
    print("  Yukawa flux     Q_r =", sp.sstr(Q_r_yukawa))
    print("  1/r^2 impostor Q_r =", sp.sstr(Q_r_power))
    print("  log impostor   Q_r =", sp.sstr(Q_r_log))
    print("  mu_eff^2(Newton) =", sp.sstr(mu_eff2_newton))
    print("  mu_eff^2(Yukawa) =", sp.sstr(mu_eff2_yukawa))
    print("  mu_eff^2(1/r^2 impostor) =", sp.sstr(mu_eff2_power))
    print("  mu_eff^2(log impostor) =", sp.sstr(mu_eff2_log))
    print("Optics/redshift falsifier:")
    print("  alpha_n =", sp.sstr(alpha_n))
    print("  alpha_fit =", sp.sstr(alpha_fit))
    print("  R_opt =", sp.sstr(R_opt))
    print("  Delta_theta(n=5) =", sp.sstr(Delta_theta))
    print("  Delta_shapiro(n=5 coefficient) =", sp.sstr(Delta_shapiro))
    print("Optional tail monitor:")
    print("  R_tail =", sp.sstr(R_tail))
    print("Concrete-flow monitor exhibit:")
    print("  periodic max |R_cont| =", periodic_summary["max_abs_R_cont"])
    print("  periodic max |R_Pois_exact| =", periodic_summary["max_abs_R_pois_exact"])
    print("  periodic alpha_fit tail =", periodic_summary["alpha_fit_tail_mean"])
    print("  Newton Q_r tail cv =", newton_summary["Q_r_tail_cv"])
    print("  Newton mu_eff^2 tail median =", newton_summary["mu_eff2_tail_median"])
    print("  Yukawa Q_r tail cv =", yukawa_summary["Q_r_tail_cv"])
    print("  Yukawa mu_eff^2 tail median =", yukawa_summary["mu_eff2_tail_median"])
    print("  Bad-optics alpha_fit tail =", bad_optics_summary["alpha_fit_tail_mean"])
    print("Hard fail criteria:")
    print("  1. Continuity / projection failure: R_cont, R_div, or R_Pois_exact do not stay small.")
    print("  2. Massive-scalar failure: Q_r(r) fails to plateau or mu_eff^2(r) tends to a nonzero exterior constant.")
    print("  3. Optical/redshift failure: alpha_fit does not plateau near 2 on the same exterior source that defines Phi_eff.")
    print("  4. Tail failure: if Theta_tail is measured independently, R_tail must stay near zero.")
    print("Interpretation:")
    print("  For CFD, the raw source to monitor is S_rho (the projected leakage term), not the Stage-45 support field phi.")
    print("  The derived gravity-like observable is the brane longitudinal potential phi_3 from the Helmholtz decomposition of v_brane.")
    print("  The decisive hard falsifiers are Yukawa-like screening in the exterior field and the wrong weak-field optical coefficient.")
    print("STATUS: PASS")


if __name__ == "__main__":
    main()
