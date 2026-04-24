#!/usr/bin/env python3
"""Generate frozen V2-22B packets from a target-blind reduced open-throat sweep.

This generator intentionally does not import the V2-21/V2-22A/V2-22B target
evaluation modules.  It only solves the reduced open-throat profiles and writes
frozen solver packets for later evaluation.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping

import numpy as np

import reduced_fem as reduced


OUTPUT_DIR = Path(__file__).resolve().parent / "output"
PACKET_DIR = OUTPUT_DIR / "packets"


def sha256_json(obj: Any) -> str:
    payload = json.dumps(obj, sort_keys=True, separators=(",", ":"), default=float).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def candidate_protocol(grid_points: int, protocol_name: str) -> Dict[str, Any]:
    """Return the predeclared, target-blind simulation protocol."""
    base_geometries = [
        {
            "name": "baseline_geometry",
            "L": 1.85,
            "a_mouth": 1.0,
            "R_exit_pref": 0.42,
            "T_R": 0.35,
            "K_R": 6.0,
            "Y_exit": 2.5,
        },
        {
            "name": "wide_soft_exit",
            "L": 1.85,
            "a_mouth": 1.0,
            "R_exit_pref": 0.55,
            "T_R": 0.32,
            "K_R": 5.0,
            "Y_exit": 2.0,
        },
        {
            "name": "narrow_stiff_exit",
            "L": 1.85,
            "a_mouth": 1.0,
            "R_exit_pref": 0.32,
            "T_R": 0.40,
            "K_R": 7.0,
            "Y_exit": 3.0,
        },
    ]
    broad_geometries = [
        *base_geometries,
        {
            "name": "long_low_tension",
            "L": 2.15,
            "a_mouth": 1.0,
            "R_exit_pref": 0.48,
            "T_R": 0.25,
            "K_R": 4.5,
            "Y_exit": 1.6,
        },
        {
            "name": "short_high_tension",
            "L": 1.55,
            "a_mouth": 1.0,
            "R_exit_pref": 0.38,
            "T_R": 0.55,
            "K_R": 8.0,
            "Y_exit": 3.5,
        },
        {
            "name": "wide_stiff_exit",
            "L": 1.70,
            "a_mouth": 1.0,
            "R_exit_pref": 0.62,
            "T_R": 0.42,
            "K_R": 7.5,
            "Y_exit": 4.0,
        },
    ]
    base_couplings = [
        {"name": "baseline_couplings", "lambda_B": 0.60, "lambda_U": 0.55, "lambda_W": 0.75, "lambda_R": 0.25},
        {"name": "weak_internal_couplings", "lambda_B": 0.45, "lambda_U": 0.45, "lambda_W": 0.60, "lambda_R": 0.15},
        {"name": "strong_transfer_couplings", "lambda_B": 0.75, "lambda_U": 0.70, "lambda_W": 0.90, "lambda_R": 0.30},
        {"name": "high_mixing_couplings", "lambda_B": 0.60, "lambda_U": 0.60, "lambda_W": 0.80, "lambda_R": 0.45},
    ]
    broad_couplings = [
        *base_couplings,
        {"name": "low_bdg_high_transfer", "lambda_B": 0.35, "lambda_U": 0.85, "lambda_W": 0.95, "lambda_R": 0.35},
        {"name": "high_bdg_low_transfer", "lambda_B": 0.95, "lambda_U": 0.35, "lambda_W": 0.45, "lambda_R": 0.12},
        {"name": "balanced_large_couplings", "lambda_B": 0.85, "lambda_U": 0.85, "lambda_W": 1.05, "lambda_R": 0.55},
        {"name": "near_mixed_limit_safe", "lambda_B": 0.70, "lambda_U": 0.75, "lambda_W": 1.10, "lambda_R": 0.90},
    ]
    operator_variants = [
        {
            "name": "baseline_operator",
            "wall_T_slope": 0.12,
            "wall_K_base": 0.32,
            "wall_K_quad": 0.08,
            "wall_TOmega_num": 0.045,
            "bdg_T_base": 0.92,
            "bdg_T_slope": 0.06,
            "bdg_V_base": 1.65,
            "bdg_V_amp": 0.18,
            "U_T_base": 1.08,
            "U_T_slope": 0.03,
            "U_V_base": 2.05,
            "U_V_slope": 0.10,
            "W_T_base": 0.86,
            "W_T_slope": 0.05,
            "W_V_base": 2.45,
            "W_V_amp": 0.12,
        },
        {
            "name": "soft_internal_modes",
            "wall_T_slope": 0.12,
            "wall_K_base": 0.34,
            "wall_K_quad": 0.07,
            "wall_TOmega_num": 0.050,
            "bdg_T_base": 0.70,
            "bdg_T_slope": 0.04,
            "bdg_V_base": 0.70,
            "bdg_V_amp": 0.08,
            "U_T_base": 0.76,
            "U_T_slope": 0.02,
            "U_V_base": 0.90,
            "U_V_slope": 0.05,
            "W_T_base": 0.70,
            "W_T_slope": 0.03,
            "W_V_base": 1.00,
            "W_V_amp": 0.05,
        },
        {
            "name": "stiff_wall_soft_modes",
            "wall_T_slope": 0.18,
            "wall_K_base": 0.62,
            "wall_K_quad": 0.14,
            "wall_TOmega_num": 0.085,
            "bdg_T_base": 0.72,
            "bdg_T_slope": 0.05,
            "bdg_V_base": 0.65,
            "bdg_V_amp": 0.08,
            "U_T_base": 0.78,
            "U_T_slope": 0.03,
            "U_V_base": 0.85,
            "U_V_slope": 0.04,
            "W_T_base": 0.72,
            "W_T_slope": 0.03,
            "W_V_base": 0.95,
            "W_V_amp": 0.05,
        },
        {
            "name": "split_internal_modes",
            "wall_T_slope": 0.10,
            "wall_K_base": 0.42,
            "wall_K_quad": 0.10,
            "wall_TOmega_num": 0.060,
            "bdg_T_base": 0.62,
            "bdg_T_slope": 0.03,
            "bdg_V_base": 0.55,
            "bdg_V_amp": 0.06,
            "U_T_base": 0.70,
            "U_T_slope": 0.02,
            "U_V_base": 0.80,
            "U_V_slope": 0.04,
            "W_T_base": 1.10,
            "W_T_slope": 0.06,
            "W_V_base": 1.75,
            "W_V_amp": 0.08,
        },
    ]
    if protocol_name == "smoke":
        geometry_variants = base_geometries
        coupling_variants = base_couplings
        selected_operator_variants = [operator_variants[0]]
    elif protocol_name == "broad_v1":
        geometry_variants = broad_geometries
        coupling_variants = broad_couplings
        selected_operator_variants = [operator_variants[0]]
    elif protocol_name == "operator_v1":
        geometry_variants = broad_geometries
        coupling_variants = broad_couplings
        selected_operator_variants = operator_variants
    else:
        raise ValueError(f"unknown protocol {protocol_name!r}")

    return {
        "schema": "pde_audit_reduced_open_throat_sweep_protocol/v1",
        "protocol_name": protocol_name,
        "description": "Target-blind reduced FEM sweep. Candidate grid is fixed before residual evaluation.",
        "grid_points": int(grid_points),
        "constants": {"G": 1.0, "c_s": 1.0, "c": 1.0, "a": 1.0, "mhat0": 1.0, "S_port": 1.0, "theta_tail": 1.0},
        "geometry_variants": geometry_variants,
        "coupling_variants": coupling_variants,
        "operator_variants": selected_operator_variants,
        "target_blind_controls": {
            "generator_imports_target_evaluation_modules": False,
            "no_post_residual_refit": True,
            "candidate_set_declared_before_evaluation": True,
            "evaluation_script": "evaluate_frozen_sweep.py",
        },
    }


def iter_candidates(protocol: Mapping[str, Any]) -> Iterable[Dict[str, Any]]:
    idx = 0
    for geometry in protocol["geometry_variants"]:
        for couplings in protocol["coupling_variants"]:
            for operator in protocol["operator_variants"]:
                idx += 1
                yield {
                    "index": idx,
                    "name": f"{geometry['name']}__{couplings['name']}__{operator['name']}",
                    "geometry": dict(geometry),
                    "couplings": dict(couplings),
                    "operator": dict(operator),
                }


def solve_candidate(protocol: Mapping[str, Any], candidate: Mapping[str, Any]) -> Dict[str, Any]:
    geometry = dict(candidate["geometry"])
    couplings = dict(candidate["couplings"])
    operator = dict(candidate["operator"])
    grid_points = int(protocol["grid_points"])
    L = float(geometry["L"])
    s = np.linspace(0.0, L, grid_points)

    shape = reduced.solve_open_shape_profile(
        s=s,
        a_mouth=float(geometry["a_mouth"]),
        R_exit_pref=float(geometry["R_exit_pref"]),
        T_R=float(geometry["T_R"]),
        K_R=float(geometry["K_R"]),
        Y_exit=float(geometry["Y_exit"]),
    )
    R0 = shape["R"]
    Rp = shape["R_prime"]
    mu_raw = R0**2 * np.sqrt(1.0 + Rp**2)
    mu_s = mu_raw * (L / reduced.trapz(mu_raw, s))
    reduced.assert_finite_array("mu_s", mu_s)

    x = s / L
    T_wall = 1.0 + float(operator["wall_T_slope"]) * x
    K_eta = float(operator["wall_K_base"]) + float(operator["wall_K_quad"]) * x**2
    T_Omega = float(operator["wall_TOmega_num"]) / (R0**2 + 0.08)
    V_wall_l2 = K_eta + 6.0 * T_Omega

    T_bdg = float(operator["bdg_T_base"]) + float(operator["bdg_T_slope"]) * (1.0 - x)
    V_bdg = float(operator["bdg_V_base"]) + float(operator["bdg_V_amp"]) * np.cos(0.5 * np.pi * x) ** 2
    T_U = float(operator["U_T_base"]) + float(operator["U_T_slope"]) * x
    V_U = float(operator["U_V_base"]) + float(operator["U_V_slope"]) * x
    T_W = float(operator["W_T_base"]) + float(operator["W_T_slope"]) * (1.0 - x)
    V_W = float(operator["W_V_base"]) + float(operator["W_V_amp"]) * np.sin(0.5 * np.pi * x) ** 2

    robin_y = 0.0
    wall_mode = reduced.solve_sturm_liouville(s, mu_s, T_wall, V_wall_l2, robin_y, 0, "wall_l2_open")
    bdg_mode = reduced.solve_sturm_liouville(s, mu_s, T_bdg, V_bdg, robin_y, 0, "bdg_support_open")
    u_mode = reduced.solve_sturm_liouville(s, mu_s, T_U, V_U, robin_y, 0, "brane_like_U_open")
    w_mode = reduced.solve_sturm_liouville(s, mu_s, T_W, V_W, robin_y, 0, "mixed_W_open")

    mode_residuals = {
        "shape_stationarity": shape["stationary_residual_relative"],
        "wall_l2": wall_mode["fem_residual_relative"],
        "bdg": bdg_mode["fem_residual_relative"],
        "U": u_mode["fem_residual_relative"],
        "W": w_mode["fem_residual_relative"],
    }
    residual_norm = max(float(v) for v in mode_residuals.values())
    freeze_payload = {
        "protocol_hash": protocol["protocol_hash"],
        "candidate_index": candidate["index"],
        "candidate_name": candidate["name"],
        "geometry": geometry,
        "couplings": couplings,
        "operator": operator,
        "grid_points": grid_points,
        "coefficient_family": "linear_FEM_open_throat_reduced_sweep_v1",
    }
    freeze_hash = sha256_json(freeze_payload)

    packet = {
        "schema": "stage_v2_22b_solver_handoff/v1",
        "branch_id": f"sim_reduced_open_throat_{candidate['index']:03d}_{candidate['name']}",
        "freeze": {
            "pre_target_freeze": True,
            "target_blind": True,
            "no_post_residual_refit": True,
            "parent_action": "GNLS + localized Maxwell + effective S_eta wall closure, reduced FEM simulation surrogate",
            "gauge_convention": "localized_or_constrained_lorenz_declared",
            "boundary_protocol": "open_impedance_AC_reflecting_DC_leaking",
            "source_map_convention": "canonical real STF grouped P2 basis",
            "support_profile_family": "linear_FEM_open_throat_reduced_sweep_v1",
            "candidate_freeze_hash": freeze_hash,
        },
        "geometry": {
            "L": L,
            "R_mouth": float(shape["R_mouth"]),
            "R_exit": float(shape["R_exit"]),
            "R_min": float(shape["R_min"]),
            "R_exit_pref": float(geometry["R_exit_pref"]),
            "boundary_class": "open_impedance",
            "Y_L_limit": 0.0,
            "exit_model": "finite_radius_impedance_opening_reduced_FEM",
        },
        "constants": dict(protocol["constants"]),
        "grid": {"coordinate": "s", "points": s.tolist()},
        "wall": {"K": float(wall_mode["eigenvalue"]), "M": 1.0},
        "profiles": {
            "weight": mu_s.tolist(),
            "wall_chi_eta": wall_mode["profile"].tolist(),
        },
        "bdg_modes": [
            {
                "name": "bdg_support_open",
                "lambda_B": float(couplings["lambda_B"]),
                "varpi": math.sqrt(float(bdg_mode["eigenvalue"])),
                "profile_values": bdg_mode["profile"].tolist(),
            }
        ],
        "mixed_ports": [
            {
                "name": "one_mixed_port",
                "lambda_U": float(couplings["lambda_U"]),
                "Omega_U": math.sqrt(float(u_mode["eigenvalue"])),
                "u_values": u_mode["profile"].tolist(),
                "lambda_W": float(couplings["lambda_W"]),
                "Omega_W": math.sqrt(float(w_mode["eigenvalue"])),
                "w_values": w_mode["profile"].tolist(),
                "lambda_R": float(couplings["lambda_R"]),
            }
        ],
        "solver_metadata": {
            "exporter": "research/pde_audit/simulation/generate_reduced_sweep.py",
            "mesh_points": grid_points,
            "nonlinear_residual_norm": residual_norm,
            "linear_residual_norm": residual_norm,
            "coefficient_family": "linear_FEM_open_throat_reduced_sweep_v1",
            "source_commit": "simulation-protocol-v1",
            "protocol_hash": protocol["protocol_hash"],
            "candidate_index": int(candidate["index"]),
            "candidate_name": str(candidate["name"]),
            "operator_variant": operator["name"],
            "candidate_freeze_hash": freeze_hash,
            "mode_residuals": mode_residuals,
        },
        "normalization_tolerances": {"profile_norm_tol": 5e-3, "Delta_tol": 1e-12},
    }
    return packet


def main() -> int:
    parser = argparse.ArgumentParser(description="Generate target-blind reduced open-throat simulation sweep")
    parser.add_argument("--output-dir", default=str(OUTPUT_DIR), help="Simulation output directory")
    parser.add_argument("--grid-points", type=int, default=81, help="Axial grid points per candidate")
    parser.add_argument("--protocol", choices=["smoke", "broad_v1", "operator_v1"], default="operator_v1", help="Predeclared target-blind candidate protocol")
    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    packet_dir = output_dir / "packets"
    if packet_dir.exists():
        for path in packet_dir.glob("*.json"):
            path.unlink()
    protocol = candidate_protocol(args.grid_points, args.protocol)
    protocol["protocol_hash"] = sha256_json({k: v for k, v in protocol.items() if k != "protocol_hash"})

    packet_entries: List[Dict[str, Any]] = []
    print("PDE audit reduced open-throat simulation generator")
    print("=" * 64)
    print(f"output_dir: {output_dir}")
    print(f"protocol_hash: {protocol['protocol_hash']}")
    print(f"protocol_name: {args.protocol}")
    print(f"grid_points: {args.grid_points}")

    for candidate in iter_candidates(protocol):
        packet = solve_candidate(protocol, candidate)
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
    write_json(output_dir / "protocol.json", protocol)
    write_json(output_dir / "manifest.json", manifest)
    print("")
    print(f"candidate_count: {len(packet_entries)}")
    print(f"manifest_hash: {manifest['manifest_hash']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
