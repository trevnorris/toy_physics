#!/usr/bin/env python3
"""Export frozen V2-22B packets from the nonlinear readiness solver.

The exporter is target-blind: it uses only the predeclared nonlinear protocol
and manufactured nonlinear solves.  Target residuals are computed later by the
post-hoc evaluator.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
from typing import Any, Dict, Mapping

import numpy as np

import nonlinear_protocol as protocol_lib
import verify_nonlinear_solver as nonlinear


OUTPUT_DIR = Path(__file__).resolve().parent / "output"


def sha256_json(obj: Any) -> str:
    payload = json.dumps(obj, sort_keys=True, separators=(",", ":"), default=float).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def normalize_profile(values: np.ndarray, weight: np.ndarray, s: np.ndarray) -> np.ndarray:
    norm = math.sqrt(nonlinear.trapz(weight * values * values, s))
    if norm <= 0.0 or not math.isfinite(norm):
        raise RuntimeError("nonpositive profile norm")
    out = values / norm
    if nonlinear.trapz(weight * out, s) < 0.0:
        out *= -1.0
    return out


def energy_scale(values: np.ndarray, weight: np.ndarray, s: np.ndarray, params: Mapping[str, Any], beta: float) -> float:
    data = nonlinear.coeffs(s, params, beta_override=beta)
    grad = np.gradient(values, s, edge_order=2)
    numerator = nonlinear.trapz(data["T"] * grad * grad + data["V"] * values * values + beta * values**4, s)
    denominator = max(nonlinear.trapz(weight * values * values, s), 1e-30)
    return float(max(numerator / denominator, 1e-12))


def solve_mode(params: Mapping[str, Any], grid_points: int, beta: float) -> Dict[str, Any]:
    result = nonlinear.newton_solve(
        params=params,
        grid_points=grid_points,
        tol=1e-10,
        max_iterations=20,
        min_step=0.001953125,
        beta_override=beta,
    )
    if not result["converged"]:
        raise RuntimeError(f"nonlinear solve failed for {params['name']} beta={beta}")
    return result


def build_packet(protocol: Mapping[str, Any], candidate: Mapping[str, Any], grid_points: int) -> Dict[str, Any]:
    params = dict(candidate["parameters"])
    base_beta = float(params["beta"])
    wall_solve = solve_mode(params, grid_points, base_beta)
    bdg_solve = solve_mode(params, grid_points, max(0.05, 0.70 * base_beta))
    u_solve = solve_mode(params, grid_points, max(0.05, 0.45 * base_beta))
    w_solve = solve_mode(params, grid_points, 1.25 * base_beta + 0.05)

    s = wall_solve["s"]
    x = s / float(params["L"])
    wall_raw = wall_solve["solution"]
    R_mouth = float(params["R_mouth"])
    R_exit = float(params["R_exit"])
    R_profile = R_exit + (R_mouth - R_exit) * (1.0 - x) + 0.025 * wall_raw * (1.0 - x) * x
    if float(np.min(R_profile)) <= 0.0:
        raise RuntimeError(f"{candidate['name']}: nonpositive radius profile")

    weight = 1.0 + 0.15 * x + 0.05 * wall_raw * wall_raw
    wall_profile = normalize_profile(wall_solve["solution"], weight, s)
    bdg_profile = normalize_profile(bdg_solve["solution"], weight, s)
    u_profile = normalize_profile(u_solve["solution"], weight, s)
    w_profile = normalize_profile(w_solve["solution"], weight, s)

    wall_K = energy_scale(wall_profile, weight, s, params, base_beta)
    bdg_K = energy_scale(bdg_profile, weight, s, params, max(0.05, 0.70 * base_beta))
    u_K = energy_scale(u_profile, weight, s, params, max(0.05, 0.45 * base_beta))
    w_K = energy_scale(w_profile, weight, s, params, 1.25 * base_beta + 0.05)
    freeze_hash = protocol_lib.candidate_freeze_hash(protocol, candidate)
    mode_residuals = {
        "wall": float(wall_solve["residual_inf"]),
        "bdg": float(bdg_solve["residual_inf"]),
        "U": float(u_solve["residual_inf"]),
        "W": float(w_solve["residual_inf"]),
    }
    residual_norm = max(mode_residuals.values())

    return {
        "schema": "stage_v2_22b_solver_handoff/v1",
        "branch_id": f"sim_nonlinear_open_throat_{int(candidate['index']):03d}_{candidate['name']}",
        "freeze": {
            "pre_target_freeze": True,
            "target_blind": True,
            "no_post_residual_refit": True,
            "parent_action": "GNLS + localized Maxwell + nonlinear open-throat continuation readiness exporter",
            "gauge_convention": "localized_or_constrained_lorenz_declared",
            "boundary_protocol": "open_impedance_AC_reflecting_DC_leaking",
            "source_map_convention": "canonical real STF grouped P2 basis",
            "support_profile_family": "nonlinear_manufactured_open_throat_v2",
            "candidate_freeze_hash": freeze_hash,
        },
        "geometry": {
            "L": float(params["L"]),
            "R_mouth": R_mouth,
            "R_exit": R_exit,
            "R_min": float(np.min(R_profile)),
            "boundary_class": "open_impedance",
            "Y_L_limit": float(params["Y_exit"]),
            "exit_model": "second_order_open_robin_impedance_nonlinear_readiness",
        },
        "constants": {"G": 1.0, "c_s": 1.0, "c": 1.0, "a": 1.0, "mhat0": 1.0, "S_port": 1.0, "theta_tail": 1.0},
        "grid": {"coordinate": "s", "points": s.tolist()},
        "wall": {"K": wall_K, "M": 1.0},
        "profiles": {
            "weight": weight.tolist(),
            "wall_chi_eta": wall_profile.tolist(),
        },
        "bdg_modes": [
            {
                "name": "nonlinear_bdg_support_open",
                "lambda_B": float(params["lambda_B"]),
                "varpi": math.sqrt(bdg_K),
                "profile_values": bdg_profile.tolist(),
            }
        ],
        "mixed_ports": [
            {
                "name": "nonlinear_mixed_port",
                "lambda_U": float(params["lambda_U"]),
                "Omega_U": math.sqrt(u_K),
                "u_values": u_profile.tolist(),
                "lambda_W": float(params["lambda_W"]),
                "Omega_W": math.sqrt(w_K),
                "w_values": w_profile.tolist(),
                "lambda_R": float(params["lambda_R"]),
            }
        ],
        "solver_metadata": {
            "exporter": "research/pde_audit/simulation/generate_nonlinear_packets.py",
            "mesh_points": int(grid_points),
            "nonlinear_residual_norm": residual_norm,
            "continuation_residual_norm": residual_norm,
            "coefficient_family": "nonlinear_manufactured_open_throat_v2",
            "source_commit": "nonlinear-protocol-v2-readiness",
            "protocol_hash": protocol["protocol_hash"],
            "candidate_index": int(candidate["index"]),
            "candidate_name": str(candidate["name"]),
            "candidate_freeze_hash": freeze_hash,
            "mode_residuals": mode_residuals,
            "manufactured_l2_errors": {
                "wall": float(wall_solve["l2_error"]),
                "bdg": float(bdg_solve["l2_error"]),
                "U": float(u_solve["l2_error"]),
                "W": float(w_solve["l2_error"]),
            },
        },
        "normalization_tolerances": {"profile_norm_tol": 5e-3, "Delta_tol": 1e-12},
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Generate target-blind nonlinear frozen V2-22B packets")
    parser.add_argument("--output-dir", default=str(OUTPUT_DIR), help="Simulation output directory")
    parser.add_argument("--grid-points", type=int, default=161, help="Axial grid points per nonlinear candidate")
    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    packet_dir = output_dir / "nonlinear_packets"
    if packet_dir.exists():
        for path in packet_dir.glob("*.json"):
            path.unlink()

    protocol = protocol_lib.build_protocol()
    packet_entries = []
    print("PDE audit nonlinear packet exporter")
    print("=" * 48)
    print(f"output_dir: {output_dir}")
    print(f"protocol_hash: {protocol['protocol_hash']}")
    print(f"grid_points: {args.grid_points}")

    for candidate in protocol_lib.iter_candidates(protocol):
        packet = build_packet(protocol, candidate, args.grid_points)
        packet_hash = sha256_json(packet)
        filename = f"{packet['branch_id']}.json"
        path = packet_dir / filename
        write_json(path, packet)
        packet_entries.append({
            "branch_id": packet["branch_id"],
            "candidate_index": int(candidate["index"]),
            "candidate_name": str(candidate["name"]),
            "path": str(path.relative_to(output_dir)),
            "packet_hash": packet_hash,
            "freeze_hash": packet["freeze"]["candidate_freeze_hash"],
            "mesh_points": int(args.grid_points),
            "solver_residual_norm": packet["solver_metadata"]["nonlinear_residual_norm"],
        })
        print(f"WROTE {filename}  residual={packet['solver_metadata']['nonlinear_residual_norm']:.3e}")

    manifest = {
        "schema": "pde_audit_simulation_manifest/v1",
        "protocol": protocol,
        "candidate_count": len(packet_entries),
        "packets": packet_entries,
        "generator_target_blind_assertion": {
            "target_evaluation_modules_imported": False,
            "target_residuals_computed": False,
            "packets_contain_target_outputs": False,
        },
    }
    manifest["manifest_hash"] = sha256_json(manifest)
    write_json(output_dir / "nonlinear_protocol.json", protocol)
    write_json(output_dir / "nonlinear_manifest.json", manifest)
    print("")
    print(f"candidate_count: {len(packet_entries)}")
    print(f"manifest_hash: {manifest['manifest_hash']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
