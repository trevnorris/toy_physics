#!/usr/bin/env python3
"""Self-test the runtime fail-fast classifier on synthetic pass/fail cases."""
from __future__ import annotations

import hashlib
import json
import pathlib
import subprocess
import sys
import tempfile

import numpy as np

from cfd_runtime_failfast import classify_summary, q_tail_cv_max_from_resolution
from cfd_runtime_monitor_postprocess import analyze_snapshot, laplacian3
from step_34_parent_throat_action_cfd_runtime_postprocessor_sympy import (
    make_periodic_consistency_snapshot,
    make_radial_snapshot,
)


def main() -> None:
    branch_metadata = {
        "branch_id": "v2_local_parent_background_failfast_classifier",
        "pre_target_freeze": True,
        "target_blind": False,
        "no_post_residual_refit": True,
        "boundary_class": "open_impedance_demo",
        "interpretation": "synthetic pass/fail self-test for the CFD-side falsification classifier",
    }

    newton_summary = analyze_snapshot(
        make_radial_snapshot(mu=None),
        c_probe=1.0,
        bins=26,
        tail_fraction=0.3,
        periodic_xyz=False,
        center=(0.0, 0.0, 0.0),
    )
    newton_summary["require_projection_source"] = False
    yukawa_summary = analyze_snapshot(
        make_radial_snapshot(mu=1.4),
        c_probe=1.0,
        bins=26,
        tail_fraction=0.3,
        periodic_xyz=False,
        center=(0.0, 0.0, 0.0),
    )
    yukawa_summary["require_projection_source"] = False

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
    bad_optics_summary["require_projection_source"] = False

    projection_break_snapshot = make_periodic_consistency_snapshot()
    x = projection_break_snapshot["x"]
    y = projection_break_snapshot["y"]
    z = projection_break_snapshot["z"]
    spacings = (float(x[1] - x[0]), float(y[1] - y[0]), float(z[1] - z[0]))
    source_like = projection_break_snapshot["rho0"] * laplacian3(projection_break_snapshot["phi3"], spacings, periodic=True)
    projection_break_snapshot["dt_rho"] = 0.3 * source_like
    projection_break_summary = analyze_snapshot(
        projection_break_snapshot,
        c_probe=1.0,
        bins=18,
        tail_fraction=0.35,
        periodic_xyz=True,
    )

    newton_verdict = classify_summary(newton_summary)
    yukawa_verdict = classify_summary(yukawa_summary)
    bad_optics_verdict = classify_summary(bad_optics_summary)
    projection_break_verdict = classify_summary(projection_break_summary)

    missing_optics_summary = dict(newton_summary)
    missing_optics_summary.pop("alpha_fit_tail_mean", None)
    missing_optics_summary.pop("alpha_fit_tail_std", None)
    missing_optics_verdict = classify_summary(missing_optics_summary)

    zero_source_summary = dict(newton_summary)
    zero_source_summary["require_projection_source"] = True
    zero_source_summary["rms_S_rho"] = 0.0
    zero_source_verdict = classify_summary(zero_source_summary)

    default_projection_summary = dict(newton_summary)
    default_projection_summary.pop("require_projection_source", None)
    default_projection_verdict = classify_summary(default_projection_summary)

    q_boundary_summary = dict(newton_summary)
    q_boundary_summary["Q_r_tail_cv"] = 0.05
    q_boundary_verdict = classify_summary(q_boundary_summary)

    q_just_outside_summary = dict(newton_summary)
    q_just_outside_summary["Q_r_tail_cv"] = 0.0501
    q_just_outside_verdict = classify_summary(q_just_outside_summary)

    mu_boundary_summary = dict(newton_summary)
    mu_boundary_summary["mu_eff2_tail_median"] = 0.25
    mu_boundary_verdict = classify_summary(mu_boundary_summary)

    mu_just_outside_summary = dict(newton_summary)
    mu_just_outside_summary["mu_eff2_tail_median"] = 0.2501
    mu_just_outside_verdict = classify_summary(mu_just_outside_summary)

    resolution_q_boundary_summary = dict(newton_summary)
    resolution_q_boundary_summary["tail_n_points"] = 200
    resolution_q_threshold = q_tail_cv_max_from_resolution(resolution_q_boundary_summary)
    resolution_q_boundary_summary["Q_r_tail_cv"] = resolution_q_threshold
    resolution_q_boundary_verdict = classify_summary(
        resolution_q_boundary_summary,
        q_tail_cv_from_resolution=True,
    )

    resolution_q_outside_summary = dict(resolution_q_boundary_summary)
    resolution_q_outside_summary["Q_r_tail_cv"] = resolution_q_threshold + 1.0e-4
    resolution_q_outside_verdict = classify_summary(
        resolution_q_outside_summary,
        q_tail_cv_from_resolution=True,
    )
    with tempfile.TemporaryDirectory() as tmpdir:
        conflict_summary_path = pathlib.Path(tmpdir) / "summary.json"
        conflict_summary_path.write_text(json.dumps(resolution_q_boundary_summary))
        failfast_cli = pathlib.Path(__file__).with_name("cfd_runtime_failfast.py")
        q_threshold_conflict = subprocess.run(
            [
                sys.executable,
                str(failfast_cli),
                str(conflict_summary_path),
                "--q-tail-cv-max",
                "0.05",
                "--q-tail-cv-max-from-resolution",
            ],
            check=False,
            capture_output=True,
            text=True,
        )

    branch_freeze_payload = {
        "metadata": branch_metadata,
        "verdicts": {
            "newton": newton_verdict["status"],
            "yukawa": yukawa_verdict["status"],
            "bad_optics": bad_optics_verdict["status"],
            "projection_break": projection_break_verdict["status"],
            "missing_optics": missing_optics_verdict["status"],
            "zero_source": zero_source_verdict["status"],
            "default_projection_source": default_projection_verdict["status"],
            "q_boundary": q_boundary_verdict["status"],
            "q_just_outside": q_just_outside_verdict["status"],
            "mu_boundary": mu_boundary_verdict["status"],
            "mu_just_outside": mu_just_outside_verdict["status"],
            "resolution_q_boundary": resolution_q_boundary_verdict["status"],
            "resolution_q_outside": resolution_q_outside_verdict["status"],
        },
        "metrics": {
            "newton_Q_r_tail_cv": newton_summary["Q_r_tail_cv"],
            "newton_mu_eff2_tail_median": newton_summary["mu_eff2_tail_median"],
            "yukawa_Q_r_tail_cv": yukawa_summary["Q_r_tail_cv"],
            "yukawa_mu_eff2_tail_median": yukawa_summary["mu_eff2_tail_median"],
            "bad_optics_alpha_fit_tail_mean": bad_optics_summary["alpha_fit_tail_mean"],
            "resolution_q_tail_cv_max": resolution_q_threshold,
        },
    }
    branch_freeze_hash = hashlib.sha256(json.dumps(branch_freeze_payload, sort_keys=True).encode("utf-8")).hexdigest()[:16]

    if newton_verdict["status"] != "PASS":
        raise AssertionError(f"expected PASS for Newton-like exterior, got {newton_verdict}")
    if yukawa_verdict["status"] != "FAIL":
        raise AssertionError(f"expected FAIL for Yukawa exterior, got {yukawa_verdict}")
    if bad_optics_verdict["status"] != "FAIL":
        raise AssertionError(f"expected FAIL for bad optics, got {bad_optics_verdict}")
    if projection_break_verdict["status"] != "FAIL":
        raise AssertionError(f"expected FAIL for projection break, got {projection_break_verdict}")
    if missing_optics_verdict["status"] != "INCOMPLETE":
        raise AssertionError(f"expected INCOMPLETE for missing optics, got {missing_optics_verdict}")
    if zero_source_verdict["status"] != "INCOMPLETE":
        raise AssertionError(f"expected INCOMPLETE for near-zero source scale, got {zero_source_verdict}")
    if not any("source scale is near zero" in warning for warning in zero_source_verdict["warnings"]):
        raise AssertionError(f"near-zero source warning did not fire: {zero_source_verdict}")
    if default_projection_verdict["status"] != "INCOMPLETE":
        raise AssertionError(f"default source-free exterior should be INCOMPLETE, got {default_projection_verdict}")
    if q_boundary_verdict["status"] != "PASS":
        raise AssertionError(f"Q_r threshold boundary should PASS, got {q_boundary_verdict}")
    if q_just_outside_verdict["status"] != "FAIL":
        raise AssertionError(f"Q_r just outside threshold should FAIL, got {q_just_outside_verdict}")
    if mu_boundary_verdict["status"] != "PASS":
        raise AssertionError(f"mu_eff2 threshold boundary should PASS, got {mu_boundary_verdict}")
    if mu_just_outside_verdict["status"] != "FAIL":
        raise AssertionError(f"mu_eff2 just outside threshold should FAIL, got {mu_just_outside_verdict}")
    if resolution_q_boundary_verdict["status"] != "PASS":
        raise AssertionError(f"resolution-derived Q_r threshold boundary should PASS, got {resolution_q_boundary_verdict}")
    if resolution_q_outside_verdict["status"] != "FAIL":
        raise AssertionError(f"resolution-derived Q_r just outside should FAIL, got {resolution_q_outside_verdict}")
    if q_threshold_conflict.returncode == 0 or "mutually exclusive" not in q_threshold_conflict.stderr:
        raise AssertionError(f"expected CLI threshold conflict to raise, got {q_threshold_conflict}")

    print("STEP 35 FAIL-FAST CLASSIFIER SELF-TEST")
    print("Validated the CFD-side falsification classifier on synthetic pass/fail cases.")
    print("V2 branch-freeze metadata:")
    print("  branch_id =", branch_metadata["branch_id"])
    print("  branch_freeze_hash =", branch_freeze_hash)
    print("  pre_target_freeze =", str(branch_metadata["pre_target_freeze"]).lower())
    print("  target_blind =", str(branch_metadata["target_blind"]).lower())
    print("  no_post_residual_refit =", str(branch_metadata["no_post_residual_refit"]).lower())
    print("  boundary_class =", branch_metadata["boundary_class"])
    print("  interpretation =", branch_metadata["interpretation"])
    print("Default fail-fast thresholds:")
    print("  continuity_rel_max = 0.05")
    print("  poisson_rel_max = 0.05")
    print("  q_tail_cv_max = 0.05")
    print("  mu_eff2_tail_abs_max = 0.25")
    print("  alpha_fit_tail_error_max = 0.1")
    print("  alpha_fit_tail_std_max = 0.1")
    print("Resolution-derived threshold:")
    print("  q_tail_cv_max(tail_n_points=200) =", resolution_q_threshold)
    print("Classifier verdicts:")
    print("  Newton-like exterior =", newton_verdict["status"])
    print("  Yukawa exterior =", yukawa_verdict["status"])
    print("  Bad optics exterior =", bad_optics_verdict["status"])
    print("  Projection-broken snapshot =", projection_break_verdict["status"])
    print("  Missing optics snapshot =", missing_optics_verdict["status"])
    print("  Near-zero source snapshot =", zero_source_verdict["status"])
    print("  Default source-free exterior =", default_projection_verdict["status"])
    print("  Q_r threshold-boundary snapshot =", q_boundary_verdict["status"])
    print("  Q_r just-outside-threshold snapshot =", q_just_outside_verdict["status"])
    print("  mu_eff2 threshold-boundary snapshot =", mu_boundary_verdict["status"])
    print("  mu_eff2 just-outside-threshold snapshot =", mu_just_outside_verdict["status"])
    print("  resolution-derived Q_r boundary snapshot =", resolution_q_boundary_verdict["status"])
    print("  resolution-derived Q_r just-outside snapshot =", resolution_q_outside_verdict["status"])
    print("  CLI Q_r threshold conflict guard = PASS")
    print("Representative failure reasons:")
    print("  Yukawa:", "; ".join(yukawa_verdict["failures"]))
    print("  Bad optics:", "; ".join(bad_optics_verdict["failures"]))
    print("  Projection break:", "; ".join(projection_break_verdict["failures"]))
    print("  Missing optics:", "; ".join(missing_optics_verdict["incomplete_reasons"]))
    print("  Near-zero source:", "; ".join(zero_source_verdict["incomplete_reasons"]))
    print("  Default source-free exterior:", "; ".join(default_projection_verdict["incomplete_reasons"]))
    print("  Q_r just outside:", "; ".join(q_just_outside_verdict["failures"]))
    print("  mu_eff2 just outside:", "; ".join(mu_just_outside_verdict["failures"]))
    print("  resolution Q_r just outside:", "; ".join(resolution_q_outside_verdict["failures"]))
    print("CLI:")
    print("  postprocess = python cfd_runtime_monitor_postprocess.py snapshot.npz --output-json summary.json")
    print("  classify    = python cfd_runtime_failfast.py summary.json --output-json verdict.json")
    print("  classify with resolution Q_r threshold = add --q-tail-cv-max-from-resolution")
    print("Interpretation:")
    print("  A real simulation snapshot can now be reduced to a fail-fast verdict with explicit reasons instead of manual note matching.")
    print("STATUS: PASS")


if __name__ == "__main__":
    main()
