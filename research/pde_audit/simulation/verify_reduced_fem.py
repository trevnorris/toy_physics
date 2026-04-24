#!/usr/bin/env python3
"""Manufactured and structural checks for reduced FEM primitives."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
from typing import Any, Dict, List

import numpy as np

import reduced_fem as reduced


OUTPUT_DIR = Path(__file__).resolve().parent / "output"


def sha256_json(obj: Any) -> str:
    payload = json.dumps(obj, sort_keys=True, separators=(",", ":"), default=float).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def matrix_structure_check() -> Dict[str, Any]:
    s = np.linspace(0.0, 1.0, 41)
    mu = 1.0 + 0.15 * s
    T = 1.0 + 0.20 * s
    V = 0.75 + 0.10 * s * s
    K, M, keep = reduced.assemble_fem_matrices(s, mu, T, V, robin_Y_right=0.05, dirichlet_left=True)
    K_sym = float(np.linalg.norm(K - K.T) / max(np.linalg.norm(K), 1e-30))
    M_sym = float(np.linalg.norm(M - M.T) / max(np.linalg.norm(M), 1e-30))
    K_min = float(np.min(np.linalg.eigvalsh(K)))
    M_min = float(np.min(np.linalg.eigvalsh(M)))
    return {
        "name": "matrix_structure",
        "pass": bool(K_sym < 1e-13 and M_sym < 1e-13 and K_min > 0.0 and M_min > 0.0 and len(keep) == 40),
        "K_symmetry_relative": K_sym,
        "M_symmetry_relative": M_sym,
        "K_min_eigenvalue": K_min,
        "M_min_eigenvalue": M_min,
        "kept_dof_count": len(keep),
    }


def dn_halfwave_convergence_check() -> Dict[str, Any]:
    L = 1.0
    exact_lambda = math.pi * math.pi / (4.0 * L * L)
    grid_points = [41, 81, 161]
    rows: List[Dict[str, float]] = []
    for n in grid_points:
        s = np.linspace(0.0, L, n)
        mu = np.ones_like(s)
        T = np.ones_like(s)
        V = np.zeros_like(s)
        result = reduced.solve_sturm_liouville(
            s=s,
            mu=mu,
            T=T,
            V=V,
            robin_Y_right=0.0,
            mode_index=0,
            label="manufactured_DN_halfwave",
        )
        profile = result["profile"]
        exact_profile = math.sqrt(2.0 / L) * np.sin(0.5 * math.pi * s / L)
        if reduced.trapz(profile * exact_profile, s) < 0.0:
            profile = -profile
        l2_error = math.sqrt(reduced.trapz((profile - exact_profile) ** 2, s))
        rows.append({
            "grid_points": float(n),
            "lambda": float(result["eigenvalue"]),
            "lambda_relative_error": abs(float(result["eigenvalue"]) - exact_lambda) / exact_lambda,
            "profile_l2_error": l2_error,
            "right_derivative_abs": abs(float(result["right_derivative"])),
            "fem_residual_relative": float(result["fem_residual_relative"]),
        })

    lambda_errors = [row["lambda_relative_error"] for row in rows]
    profile_errors = [row["profile_l2_error"] for row in rows]
    right_derivative_errors = [row["right_derivative_abs"] for row in rows]
    monotone_lambda = all(lambda_errors[i + 1] < lambda_errors[i] for i in range(len(lambda_errors) - 1))
    monotone_right_derivative = all(
        right_derivative_errors[i + 1] < right_derivative_errors[i]
        for i in range(len(right_derivative_errors) - 1)
    )
    return {
        "name": "manufactured_DN_halfwave",
        "pass": bool(
            monotone_lambda
            and monotone_right_derivative
            and lambda_errors[-1] < 1e-4
            and max(profile_errors) < 1e-10
        ),
        "exact_lambda": exact_lambda,
        "rows": rows,
        "max_profile_l2_error": max(profile_errors),
        "lambda_error_ratio_41_to_81": lambda_errors[0] / lambda_errors[1],
        "lambda_error_ratio_81_to_161": lambda_errors[1] / lambda_errors[2],
    }


def open_shape_smoke_check() -> Dict[str, Any]:
    s = np.linspace(0.0, 1.85, 81)
    shape = reduced.solve_open_shape_profile(
        s=s,
        a_mouth=1.0,
        R_exit_pref=0.42,
        T_R=0.35,
        K_R=6.0,
        Y_exit=2.5,
    )
    R = shape["R"]
    Rp = shape["R_prime"]
    finite = bool(np.all(np.isfinite(R)) and np.all(np.isfinite(Rp)))
    return {
        "name": "open_shape_profile_smoke",
        "pass": bool(finite and shape["R_min"] > 0.0 and shape["stationary_residual_relative"] < 1e-12),
        "finite": finite,
        "R_mouth": float(shape["R_mouth"]),
        "R_exit": float(shape["R_exit"]),
        "R_min": float(shape["R_min"]),
        "stationary_residual_relative": float(shape["stationary_residual_relative"]),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Verify reduced FEM primitives used by simulation generation")
    parser.add_argument("--output-dir", default=str(OUTPUT_DIR), help="Simulation output directory")
    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    checks = [matrix_structure_check(), dn_halfwave_convergence_check(), open_shape_smoke_check()]
    report = {
        "schema": "pde_audit_simulation_reduced_fem_verification/v1",
        "checks": checks,
        "pass": all(bool(check["pass"]) for check in checks),
        "failed_checks": [check["name"] for check in checks if not check["pass"]],
    }
    report["report_hash"] = sha256_json(report)
    write_json(output_dir / "reduced_fem_verification.json", report)

    lines = [
        "PDE audit reduced FEM verification",
        "=" * 44,
        f"pass: {report['pass']}",
        f"checks_passed: {sum(1 for check in checks if check['pass'])}/{len(checks)}",
        f"failed_checks: {', '.join(report['failed_checks']) if report['failed_checks'] else 'none'}",
        f"report_hash: {report['report_hash']}",
    ]
    for check in checks:
        lines.append(f"{'PASS' if check['pass'] else 'FAIL'}  {check['name']}")
    (output_dir / "reduced_fem_verification.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines))
    return 0 if report["pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
