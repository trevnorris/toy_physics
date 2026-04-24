#!/usr/bin/env python3
"""Pre-target nonlinear solver readiness checks.

This script does not emit V2-22B candidate packets and does not import target
evaluation modules.  It verifies the numerical and freeze-boundary pieces that
must be in place before a nonlinear exporter is allowed to feed the existing
target-evaluation chain.
"""

from __future__ import annotations

import argparse
import ast
import hashlib
import json
import math
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Sequence, Set, Tuple

import numpy as np

import nonlinear_protocol as protocol_lib


SIM_DIR = Path(__file__).resolve().parent
OUTPUT_DIR = SIM_DIR / "output"
BANNED_IMPORT_PREFIXES = ("stage_v2_",)
BANNED_TARGET_KEYS = {
    "dominant_score_component",
    "pass_flags",
    "R_norm",
    "R_P2",
    "R_P4",
    "R_pole",
    "R_tail",
    "residuals",
    "score",
    "score_components",
    "target_packet_pass",
    "target_values",
}
ALLOWED_IMPORTS = {
    "nonlinear_protocol.py": {"__future__", "hashlib", "json", "typing"},
    "verify_nonlinear_solver.py": {
        "__future__",
        "argparse",
        "ast",
        "hashlib",
        "json",
        "math",
        "nonlinear_protocol",
        "numpy",
        "pathlib",
        "typing",
    },
}


def sha256_json(obj: Any) -> str:
    payload = json.dumps(obj, sort_keys=True, separators=(",", ":"), default=float).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def trapz(y: np.ndarray, x: np.ndarray) -> float:
    trapezoid = getattr(np, "trapezoid", np.trapz)
    return float(trapezoid(y, x))


def imported_roots(path: Path) -> Set[str]:
    tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
    roots: Set[str] = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            for alias in node.names:
                roots.add(alias.name.split(".", 1)[0])
        elif isinstance(node, ast.ImportFrom) and node.module:
            roots.add(node.module.split(".", 1)[0] if node.level == 0 else "__future__")
    return roots


def find_forbidden_keys(obj: Any, prefix: str = "") -> List[str]:
    paths: List[str] = []
    if isinstance(obj, dict):
        for key, value in obj.items():
            key_text = str(key)
            path = f"{prefix}.{key_text}" if prefix else key_text
            if key_text in BANNED_TARGET_KEYS:
                paths.append(path)
            paths.extend(find_forbidden_keys(value, path))
    elif isinstance(obj, list):
        for idx, value in enumerate(obj):
            paths.extend(find_forbidden_keys(value, f"{prefix}[{idx}]"))
    return paths


def coeffs(s: np.ndarray, params: Mapping[str, Any], beta_override: Optional[float] = None) -> Dict[str, np.ndarray]:
    L = float(params["L"])
    x = s / L
    beta = float(params["beta"] if beta_override is None else beta_override)
    T = 1.0 + float(params["T_slope"]) * x
    Tp = np.full_like(s, float(params["T_slope"]) / L)
    V = float(params["V_base"]) + float(params["V_quad"]) * x * x
    k = 0.5 * math.pi / L
    u = np.sin(k * s)
    up = k * np.cos(k * s)
    upp = -k * k * np.sin(k * s)
    f = -(Tp * up + T * upp) + V * u + beta * u**3
    g_exit = float(T[-1] * up[-1] + float(params["Y_exit"]) * u[-1])
    return {"T": T, "V": V, "f": f, "u_exact": u, "g_exit": np.array([g_exit]), "beta": np.array([beta])}


def residual_and_jacobian(
    y: np.ndarray,
    params: Mapping[str, Any],
    beta_override: Optional[float] = None,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    n = len(y) + 1
    L = float(params["L"])
    s = np.linspace(0.0, L, n)
    h = s[1] - s[0]
    data = coeffs(s, params, beta_override=beta_override)
    T = data["T"]
    V = data["V"]
    f = data["f"]
    beta = float(data["beta"][0])
    g_exit = float(data["g_exit"][0])
    Y = float(params["Y_exit"])

    u = np.zeros(n, dtype=float)
    u[1:] = y
    T_half = 0.5 * (T[:-1] + T[1:])
    F = np.zeros(n - 1, dtype=float)
    J = np.zeros((n - 1, n - 1), dtype=float)

    for i in range(1, n - 1):
        row = i - 1
        flux_right = T_half[i] * (u[i + 1] - u[i]) / h
        flux_left = T_half[i - 1] * (u[i] - u[i - 1]) / h
        F[row] = -(flux_right - flux_left) / h + V[i] * u[i] + beta * u[i] ** 3 - f[i]
        J[row, i - 1] = (T_half[i] + T_half[i - 1]) / (h * h) + V[i] + 3.0 * beta * u[i] ** 2
        if i - 1 >= 1:
            J[row, i - 2] = -T_half[i - 1] / (h * h)
        if i + 1 <= n - 1:
            J[row, i] = -T_half[i] / (h * h)

    row = n - 2
    F[row] = T[-1] * (3.0 * u[-1] - 4.0 * u[-2] + u[-3]) / (2.0 * h) + Y * u[-1] - g_exit
    J[row, n - 2] = 3.0 * T[-1] / (2.0 * h) + Y
    J[row, n - 3] = -4.0 * T[-1] / (2.0 * h)
    if n > 3:
        J[row, n - 4] = T[-1] / (2.0 * h)
    return F, J, s, data["u_exact"]


def newton_solve(
    params: Mapping[str, Any],
    grid_points: int,
    tol: float,
    max_iterations: int,
    min_step: float,
    beta_override: Optional[float] = None,
    initial: Optional[np.ndarray] = None,
) -> Dict[str, Any]:
    L = float(params["L"])
    s = np.linspace(0.0, L, grid_points)
    exact = coeffs(s, params, beta_override=beta_override)["u_exact"]
    y = exact[1:] * 0.85 if initial is None else np.array(initial, dtype=float).copy()
    history: List[float] = []
    converged = False

    for iteration in range(max_iterations + 1):
        F, J, _, _ = residual_and_jacobian(y, params, beta_override=beta_override)
        norm = float(np.linalg.norm(F, ord=np.inf))
        history.append(norm)
        if norm < tol:
            converged = True
            break
        if iteration == max_iterations:
            break
        delta = np.linalg.solve(J, -F)
        step = 1.0
        accepted = False
        while step >= min_step:
            y_trial = y + step * delta
            F_trial, _, _, _ = residual_and_jacobian(y_trial, params, beta_override=beta_override)
            if float(np.linalg.norm(F_trial, ord=np.inf)) < norm:
                y = y_trial
                accepted = True
                break
            step *= 0.5
        if not accepted:
            y = y + delta

    u = np.zeros(grid_points, dtype=float)
    u[1:] = y
    l2_error = math.sqrt(trapz((u - exact) ** 2, s))
    return {
        "converged": converged,
        "iterations": len(history) - 1,
        "residual_inf": history[-1],
        "residual_history": history,
        "solution": u,
        "s": s,
        "exact": exact,
        "l2_error": l2_error,
    }


def source_import_check() -> Dict[str, Any]:
    issues: List[str] = []
    details: Dict[str, Any] = {}
    for filename, allowed in ALLOWED_IMPORTS.items():
        imports = imported_roots(SIM_DIR / filename)
        banned = sorted(root for root in imports if any(root.startswith(prefix) for prefix in BANNED_IMPORT_PREFIXES))
        unexpected = sorted(imports - allowed)
        details[filename] = {
            "imports": sorted(imports),
            "allowed_imports": sorted(allowed),
            "banned_imports": banned,
            "unexpected_imports": unexpected,
        }
        for root in banned:
            issues.append(f"{filename} imports banned target-evaluation module root {root!r}")
        for root in unexpected:
            issues.append(f"{filename} imports undeclared module root {root!r}")
    return {"name": "source_import_boundary", "pass": not issues, "issues": issues, "details": details}


def protocol_hash_check(protocol: Mapping[str, Any]) -> Dict[str, Any]:
    protocol_without_hash = dict(protocol)
    protocol_without_hash.pop("protocol_hash", None)
    recomputed = sha256_json(protocol_without_hash)
    candidates = list(protocol_lib.iter_candidates(protocol))
    freeze_hashes = [protocol_lib.candidate_freeze_hash(protocol, candidate) for candidate in candidates]
    forbidden_paths = find_forbidden_keys(protocol)
    ok = (
        protocol.get("protocol_hash") == recomputed
        and len(freeze_hashes) == len(set(freeze_hashes))
        and not forbidden_paths
        and protocol.get("target_blind_controls", {}).get("target_residuals_computed") is False
    )
    return {
        "name": "protocol_and_freeze_hashes",
        "pass": bool(ok),
        "stored_protocol_hash": protocol.get("protocol_hash"),
        "recomputed_protocol_hash": recomputed,
        "candidate_count": len(candidates),
        "unique_candidate_freeze_hash_count": len(set(freeze_hashes)),
        "candidate_freeze_hashes": [
            {"candidate": candidate["name"], "freeze_hash": freeze_hash}
            for candidate, freeze_hash in zip(candidates, freeze_hashes)
        ],
        "forbidden_target_key_paths": forbidden_paths,
    }


def jacobian_consistency_check(protocol: Mapping[str, Any]) -> Dict[str, Any]:
    params = next(protocol_lib.iter_candidates(protocol))["parameters"]
    n = int(protocol["jacobian_check"]["grid_points"])
    eps = float(protocol["jacobian_check"]["finite_difference_step"])
    tol = float(protocol["jacobian_check"]["relative_directional_tolerance"])
    s = np.linspace(0.0, float(params["L"]), n)
    exact = coeffs(s, params)["u_exact"]
    y = 0.80 * exact[1:] + 0.03 * np.sin(2.0 * math.pi * s[1:] / float(params["L"]))
    direction = np.cos(1.7 * s[1:]) + 0.25 * np.sin(3.1 * s[1:])
    direction = direction / max(float(np.linalg.norm(direction)), 1e-30)
    F_plus, _, _, _ = residual_and_jacobian(y + eps * direction, params)
    F_minus, _, _, _ = residual_and_jacobian(y - eps * direction, params)
    _, J, _, _ = residual_and_jacobian(y, params)
    fd = (F_plus - F_minus) / (2.0 * eps)
    analytic = J @ direction
    rel = float(np.linalg.norm(fd - analytic) / max(float(np.linalg.norm(analytic)), 1e-30))
    return {
        "name": "jacobian_directional_consistency",
        "pass": bool(rel < tol),
        "relative_error": rel,
        "tolerance": tol,
        "grid_points": n,
        "finite_difference_step": eps,
    }


def manufactured_mesh_check(protocol: Mapping[str, Any]) -> Dict[str, Any]:
    rows: List[Dict[str, Any]] = []
    tol = float(protocol["newton"]["residual_inf_tol"])
    max_iterations = int(protocol["newton"]["max_iterations"])
    min_step = float(protocol["newton"]["line_search_min_step"])
    grid_points = [int(n) for n in protocol["grid_points_for_readiness"]]
    finest_max = float(protocol["mesh_check"]["l2_error_finest_max"])

    for candidate in protocol_lib.iter_candidates(protocol):
        errors: List[float] = []
        candidate_rows: List[Dict[str, Any]] = []
        for n in grid_points:
            result = newton_solve(candidate["parameters"], n, tol, max_iterations, min_step)
            errors.append(float(result["l2_error"]))
            candidate_rows.append({
                "grid_points": n,
                "converged": result["converged"],
                "iterations": result["iterations"],
                "residual_inf": result["residual_inf"],
                "l2_error": result["l2_error"],
            })
        rows.append({
            "candidate": candidate["name"],
            "rows": candidate_rows,
            "monotone_l2_error": all(errors[i + 1] < errors[i] for i in range(len(errors) - 1)),
            "finest_l2_error": errors[-1],
        })

    ok = all(
        item["monotone_l2_error"]
        and item["finest_l2_error"] < finest_max
        and all(row["converged"] and row["residual_inf"] < tol for row in item["rows"])
        for item in rows
    )
    return {
        "name": "manufactured_nonlinear_mesh_convergence",
        "pass": bool(ok),
        "finest_l2_error_max": finest_max,
        "rows": rows,
    }


def continuation_check(protocol: Mapping[str, Any]) -> Dict[str, Any]:
    candidate = next(protocol_lib.iter_candidates(protocol))
    params = candidate["parameters"]
    beta_values = [float(v) for v in protocol["continuation"]["values"]]
    grid_points = 81
    tol = float(protocol["newton"]["residual_inf_tol"])
    max_iterations = int(protocol["newton"]["max_iterations"])
    min_step = float(protocol["newton"]["line_search_min_step"])
    rows: List[Dict[str, Any]] = []
    previous: Optional[np.ndarray] = None

    for beta in beta_values:
        result = newton_solve(
            params,
            grid_points,
            tol,
            max_iterations,
            min_step,
            beta_override=beta,
            initial=previous[1:] if previous is not None else None,
        )
        previous = result["solution"]
        rows.append({
            "beta": beta,
            "converged": result["converged"],
            "iterations": result["iterations"],
            "residual_inf": result["residual_inf"],
            "l2_error": result["l2_error"],
            "solution_min": float(np.min(result["solution"])),
            "solution_max": float(np.max(result["solution"])),
        })

    ok = all(row["converged"] and row["residual_inf"] < tol for row in rows)
    return {"name": "continuation_beta_sequence", "pass": bool(ok), "grid_points": grid_points, "rows": rows}


def main() -> int:
    parser = argparse.ArgumentParser(description="Verify nonlinear PDE/continuation readiness without target evaluation")
    parser.add_argument("--output-dir", default=str(OUTPUT_DIR), help="Simulation output directory")
    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    protocol = protocol_lib.build_protocol()
    write_json(output_dir / "nonlinear_protocol.json", protocol)

    checks = [
        source_import_check(),
        protocol_hash_check(protocol),
        jacobian_consistency_check(protocol),
        manufactured_mesh_check(protocol),
        continuation_check(protocol),
    ]
    report = {
        "schema": "pde_audit_nonlinear_readiness/v1",
        "protocol_path": str(output_dir / "nonlinear_protocol.json"),
        "protocol_hash": protocol["protocol_hash"],
        "packets_emitted": False,
        "target_residuals_computed": False,
        "checks": checks,
        "pass": all(bool(check["pass"]) for check in checks),
        "failed_checks": [check["name"] for check in checks if not check["pass"]],
    }
    report["report_hash"] = sha256_json(report)
    write_json(output_dir / "nonlinear_readiness.json", report)

    lines = [
        "PDE audit nonlinear readiness",
        "=" * 36,
        f"pass: {report['pass']}",
        f"checks_passed: {sum(1 for check in checks if check['pass'])}/{len(checks)}",
        f"packets_emitted: {report['packets_emitted']}",
        f"target_residuals_computed: {report['target_residuals_computed']}",
        f"protocol_hash: {report['protocol_hash']}",
        f"report_hash: {report['report_hash']}",
    ]
    for check in checks:
        lines.append(f"{'PASS' if check['pass'] else 'FAIL'}  {check['name']}")
    (output_dir / "nonlinear_readiness.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines))
    return 0 if report["pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
