#!/usr/bin/env python3
"""Classify runtime-monitor summaries into fast falsification verdicts."""
from __future__ import annotations

import argparse
import json
import math
import pathlib
from typing import Any


DEFAULT_THRESHOLDS = {
    "continuity_rel_max": 5.0e-2,
    "poisson_rel_max": 5.0e-2,
    "q_tail_cv_max": 5.0e-2,
    "mu_eff2_tail_abs_max": 2.5e-1,
    "alpha_fit_tail_error_max": 1.0e-1,
    "alpha_fit_tail_std_max": 1.0e-1,
}


def load_summary(path: pathlib.Path) -> dict[str, Any]:
    return json.loads(path.read_text())


def classify_summary(summary: dict[str, Any], thresholds: dict[str, float] | None = None) -> dict[str, Any]:
    limits = dict(DEFAULT_THRESHOLDS)
    if thresholds:
        limits.update(thresholds)

    failures: list[str] = []
    warnings: list[str] = []
    incomplete_reasons: list[str] = []

    source_scale = float(summary.get("rms_S_rho", 0.0))
    require_projection_source = bool(summary.get("require_projection_source", True))
    continuity_rel = math.nan
    poisson_rel = math.nan
    if source_scale > 1.0e-12:
        continuity_rel = float(summary["rms_R_cont"]) / source_scale
        poisson_rel = float(summary["rms_R_pois_exact"]) / source_scale
        if continuity_rel > limits["continuity_rel_max"]:
            failures.append(f"continuity residual too large: {continuity_rel:.6g}")
        if poisson_rel > limits["poisson_rel_max"]:
            failures.append(f"exact Poisson residual too large: {poisson_rel:.6g}")
    else:
        warnings.append("source scale is near zero; projection residuals were not normalized")
        if require_projection_source:
            incomplete_reasons.append("source scale is near zero; projection residuals are not load-bearing")

    q_tail_cv = float(summary.get("Q_r_tail_cv", math.nan))
    if math.isfinite(q_tail_cv) and q_tail_cv > limits["q_tail_cv_max"]:
        failures.append(f"exterior Q_r plateau is too noisy: {q_tail_cv:.6g}")

    mu_tail = float(summary.get("mu_eff2_tail_median", math.nan))
    if math.isfinite(mu_tail) and abs(mu_tail) > limits["mu_eff2_tail_abs_max"]:
        failures.append(f"effective exterior mass scale is too large: {mu_tail:.6g}")

    alpha_mean = summary.get("alpha_fit_tail_mean")
    alpha_std = summary.get("alpha_fit_tail_std")
    if alpha_mean is None or alpha_std is None:
        warnings.append("optics channel missing; alpha_fit test incomplete")
        incomplete_reasons.append("optics channel missing; alpha_fit test incomplete")
    else:
        alpha_mean = float(alpha_mean)
        alpha_std = float(alpha_std)
        if abs(alpha_mean - 2.0) > limits["alpha_fit_tail_error_max"]:
            failures.append(f"alpha_fit mean is not close to 2: {alpha_mean:.6g}")
        if alpha_std > limits["alpha_fit_tail_std_max"]:
            failures.append(f"alpha_fit tail is not stable: std={alpha_std:.6g}")

    if failures:
        status = "FAIL"
    elif incomplete_reasons:
        status = "INCOMPLETE"
    else:
        status = "PASS"

    return {
        "status": status,
        "thresholds": limits,
        "metrics": {
            "source_scale_rms": source_scale,
            "continuity_rel": continuity_rel,
            "poisson_rel": poisson_rel,
            "Q_r_tail_cv": q_tail_cv,
            "mu_eff2_tail_median": mu_tail,
            "alpha_fit_tail_mean": alpha_mean,
            "alpha_fit_tail_std": alpha_std,
        },
        "failures": failures,
        "warnings": warnings,
        "incomplete_reasons": incomplete_reasons,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("summary_json", help="JSON file produced by cfd_runtime_monitor_postprocess.py")
    parser.add_argument("--output-json", help="write verdict JSON to this path")
    args = parser.parse_args()

    verdict = classify_summary(load_summary(pathlib.Path(args.summary_json)))
    print("CFD FAIL-FAST VERDICT")
    print("  status           =", verdict["status"])
    for key, value in verdict["metrics"].items():
        print(f"  {key:17} = {value}")
    if verdict["failures"]:
        print("  failures:")
        for item in verdict["failures"]:
            print("   -", item)
    if verdict["warnings"]:
        print("  warnings:")
        for item in verdict["warnings"]:
            print("   -", item)
    if verdict["incomplete_reasons"]:
        print("  incomplete reasons:")
        for item in verdict["incomplete_reasons"]:
            print("   -", item)

    if args.output_json:
        pathlib.Path(args.output_json).write_text(json.dumps(verdict, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
