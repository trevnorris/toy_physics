#!/usr/bin/env python3
"""Stage V2-23 — Mesh-convergence audit for the reduced open-throat solver."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping, Sequence

import stage_v2_23_minimal_open_throat_branch_solver as solver


KEYS = [
    "lambda_wall_l2",
    "lambda_bdg",
    "lambda_U",
    "lambda_W",
    "I_eta_phi",
    "I_eta_u",
    "I_eta_w",
    "I_uw",
    "P0",
    "P2",
    "P4",
    "R_pole",
    "R_norm",
]


def relative_delta(a: float, b: float) -> float:
    return abs(float(a) - float(b)) / max(abs(float(a)), abs(float(b)), 1e-30)


def summarize_run(grid_points: int) -> Dict[str, Any]:
    branch = solver.build_branch(grid_points=grid_points)
    coeffs = solver.compute_coefficients(branch)
    gates = solver.evaluate_gates(branch, coeffs)
    return {
        "grid_points": int(grid_points),
        "branch_freeze_hash": branch["branch_freeze_hash"],
        "gates": {
            "open_gate_pass": gates["open_gate_pass"],
            "stability_gate_pass": gates["stability_gate_pass"],
            "outgoing_transfer_gate_pass": gates["outgoing_transfer_gate_pass"],
            "target_packet_pass": gates["target_packet_pass"],
        },
        "mode_residuals": {
            "wall_l2": branch["mode_diagnostics"]["wall"]["fem_residual_relative"],
            "bdg": branch["mode_diagnostics"]["bdg"]["fem_residual_relative"],
            "U": branch["mode_diagnostics"]["U"]["fem_residual_relative"],
            "W": branch["mode_diagnostics"]["W"]["fem_residual_relative"],
            "shape_stationarity": branch["geometry"]["stationary_residual_relative"],
        },
        "values": {
            "lambda_wall_l2": coeffs["wall"]["K"],
            "lambda_bdg": coeffs["bdg"]["varpi2"],
            "lambda_U": coeffs["mixed"]["Omega_U2"],
            "lambda_W": coeffs["mixed"]["Omega_W2"],
            "I_eta_phi": coeffs["overlaps"]["I_eta_phi"],
            "I_eta_u": coeffs["overlaps"]["I_eta_u"],
            "I_eta_w": coeffs["overlaps"]["I_eta_w"],
            "I_uw": coeffs["overlaps"]["I_uw"],
            "P0": coeffs["prefactor"]["P0"],
            "P2": coeffs["prefactor"]["P2"],
            "P4": coeffs["prefactor"]["P4"],
            "R_pole": coeffs["targets"]["R_pole"],
            "R_norm": coeffs["targets"]["R_norm"],
        },
    }


def pairwise_deltas(runs: Sequence[Mapping[str, Any]], keys: Iterable[str]) -> List[Dict[str, Any]]:
    deltas: List[Dict[str, Any]] = []
    for left, right in zip(runs, runs[1:]):
        by_key = {
            key: relative_delta(left["values"][key], right["values"][key])
            for key in keys
        }
        deltas.append({
            "from_grid": left["grid_points"],
            "to_grid": right["grid_points"],
            "max_relative_delta": max(by_key.values()),
            "relative_deltas": by_key,
        })
    return deltas


def convergence_orders(deltas: Sequence[Mapping[str, Any]], keys: Iterable[str]) -> Dict[str, float]:
    if len(deltas) < 2:
        return {}
    coarse = deltas[-2]
    fine = deltas[-1]
    h_coarse = 1.0 / (float(coarse["to_grid"]) - 1.0)
    h_fine = 1.0 / (float(fine["to_grid"]) - 1.0)
    orders: Dict[str, float] = {}
    for key in keys:
        e_coarse = float(coarse["relative_deltas"][key])
        e_fine = float(fine["relative_deltas"][key])
        if e_coarse > 0.0 and e_fine > 0.0 and h_coarse > h_fine:
            orders[key] = math.log(e_coarse / e_fine) / math.log(h_coarse / h_fine)
    return orders


def parse_grids(raw: str) -> List[int]:
    grids = [int(part.strip()) for part in raw.split(",") if part.strip()]
    if len(grids) < 3:
        raise argparse.ArgumentTypeError("at least three comma-separated grid sizes are required")
    if any(n < 25 for n in grids):
        raise argparse.ArgumentTypeError("all grid sizes must be at least 25")
    if any(grids[i + 1] <= grids[i] for i in range(len(grids) - 1)):
        raise argparse.ArgumentTypeError("grid sizes must be strictly increasing")
    return grids


def main() -> int:
    parser = argparse.ArgumentParser(description="Stage V2-23 mesh-convergence audit")
    parser.add_argument("--grids", type=parse_grids, default=parse_grids("91,121,181,241"), help="Comma-separated increasing grid sizes")
    parser.add_argument("--finest-pair-rel-tol", type=float, default=1e-4, help="Maximum relative delta allowed between the two finest grids")
    parser.add_argument("--residual-tol", type=float, default=1e-8, help="Maximum FEM/shape residual allowed on each run")
    parser.add_argument("--out-json", default=None, help="Optional JSON report path")
    args = parser.parse_args()

    formula = solver.sympy_formula_checks()
    runs = [summarize_run(n) for n in args.grids]
    deltas = pairwise_deltas(runs, KEYS)
    finest_pair = deltas[-1]
    observed_orders = convergence_orders(deltas, KEYS)

    all_gate_protocols_pass = all(
        run["gates"]["open_gate_pass"]
        and run["gates"]["stability_gate_pass"]
        and run["gates"]["outgoing_transfer_gate_pass"]
        and not run["gates"]["target_packet_pass"]
        for run in runs
    )
    residual_max = max(max(abs(float(v)) for v in run["mode_residuals"].values()) for run in runs)
    finite_values = all(
        math.isfinite(float(value))
        for run in runs
        for value in list(run["values"].values()) + list(run["mode_residuals"].values())
    )
    finest_pair_pass = finest_pair["max_relative_delta"] <= args.finest_pair_rel_tol
    residual_pass = residual_max <= args.residual_tol
    formula_pass = formula["passed"] == formula["total"]
    freeze_hash_stable = len({run["branch_freeze_hash"] for run in runs}) == 1
    pairwise_delta_decreases = all(
        deltas[i + 1]["max_relative_delta"] <= deltas[i]["max_relative_delta"]
        for i in range(len(deltas) - 1)
    )

    checks = {
        "formula_pass": bool(formula_pass),
        "finite_values": bool(finite_values),
        "all_open_stable_transfer_and_target_failing": bool(all_gate_protocols_pass),
        "freeze_hash_stable_across_grids": bool(freeze_hash_stable),
        "pairwise_max_delta_decreases": bool(pairwise_delta_decreases),
        "finest_pair_converged": bool(finest_pair_pass),
        "residuals_below_tolerance": bool(residual_pass),
    }
    all_pass = all(checks.values())
    report = {
        "schema": "stage_v2_23_mesh_convergence_report/v1",
        "grids": args.grids,
        "keys": KEYS,
        "tolerances": {
            "finest_pair_rel_tol": args.finest_pair_rel_tol,
            "residual_tol": args.residual_tol,
        },
        "formula": formula,
        "runs": runs,
        "pairwise_deltas": deltas,
        "observed_orders_from_last_two_pairs": observed_orders,
        "residual_max": residual_max,
        "checks": checks,
        "pass": bool(all_pass),
    }
    report["report_hash"] = solver.sha256_json(report)

    if args.out_json:
        out_path = Path(args.out_json)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        out_path.write_text(json.dumps(report, indent=2, sort_keys=True), encoding="utf-8")

    print("Stage V2-23: mesh-convergence audit")
    print("=" * 56)
    print(f"formula_sympy_audit: {formula['passed']}/{formula['total']} checks passed")
    print(f"grids: {', '.join(str(n) for n in args.grids)}")
    print(f"finest_pair: {finest_pair['from_grid']} -> {finest_pair['to_grid']}")
    print(f"finest_pair_max_relative_delta: {finest_pair['max_relative_delta']:.3e}")
    print(f"finest_pair_rel_tol: {args.finest_pair_rel_tol:.3e}")
    print(f"residual_max: {residual_max:.3e}")
    print(f"residual_tol: {args.residual_tol:.3e}")
    if observed_orders:
        min_order_key = min(observed_orders, key=observed_orders.get)
        max_order_key = max(observed_orders, key=observed_orders.get)
        print(f"observed_order_min: {observed_orders[min_order_key]:.3g} ({min_order_key})")
        print(f"observed_order_max: {observed_orders[max_order_key]:.3g} ({max_order_key})")
    for name, passed in checks.items():
        print(f"{name}: {passed}")
    print(f"mesh_convergence_pass: {all_pass}")
    return 0 if all_pass else 1


if __name__ == "__main__":
    raise SystemExit(main())
