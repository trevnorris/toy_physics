#!/usr/bin/env python3
"""Pre-target protocol helpers for the nonlinear simulation readiness layer."""

from __future__ import annotations

import hashlib
import json
from typing import Any, Dict, Iterable, Mapping


def stable_hash(obj: Any) -> str:
    payload = json.dumps(obj, sort_keys=True, separators=(",", ":"), default=float).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def build_protocol() -> Dict[str, Any]:
    protocol: Dict[str, Any] = {
        "schema": "pde_audit_nonlinear_protocol/v2",
        "protocol_name": "nonlinear_v2_readiness",
        "description": "Pre-target nonlinear PDE/continuation readiness protocol. This does not evaluate target residuals.",
        "solver_family": "scalar_nonlinear_reaction_diffusion_manufactured_open_throat_readiness",
        "grid_points_for_readiness": [41, 81, 161],
        "newton": {
            "residual_inf_tol": 1e-10,
            "max_iterations": 20,
            "line_search_min_step": 0.001953125,
        },
        "jacobian_check": {
            "grid_points": 51,
            "relative_directional_tolerance": 1e-7,
            "finite_difference_step": 1e-6,
        },
        "mesh_check": {
            "l2_error_finest_max": 2e-4,
            "require_monotone_l2_error": True,
            "require_monotone_residual": False,
        },
        "boundary_conditions": {
            "mouth": "Dirichlet",
            "exit": "open Robin impedance",
            "exit_derivative_order": 2,
            "forbidden": ["hard_cap", "closed_throat_endpoint_substitution", "post_target_refit"],
        },
        "manufactured_solution": {
            "u_exact": "sin(pi*s/(2*L))",
            "operator": "-d/ds(T(s) du/ds) + V(s) u + beta u^3 = f(s)",
            "T": "1 + T_slope*s/L",
            "V": "V_base + V_quad*(s/L)^2",
            "robin_exit": "T(L) u'(L) + Y_exit u(L) = g_exit",
        },
        "candidate_variants": [
            {
                "name": "baseline_open_nonlinear",
                "L": 1.0,
                "R_mouth": 1.0,
                "R_exit": 0.42,
                "T_slope": 0.20,
                "V_base": 0.70,
                "V_quad": 0.15,
                "beta": 0.45,
                "Y_exit": 0.35,
                "lambda_B": 0.42,
                "lambda_U": 0.28,
                "lambda_W": 0.52,
                "lambda_R": 0.35,
            },
            {
                "name": "stiffer_reaction_open",
                "L": 1.0,
                "R_mouth": 1.0,
                "R_exit": 0.35,
                "T_slope": 0.35,
                "V_base": 1.05,
                "V_quad": 0.25,
                "beta": 0.75,
                "Y_exit": 0.50,
                "lambda_B": 0.55,
                "lambda_U": 0.34,
                "lambda_W": 0.64,
                "lambda_R": 0.42,
            },
            {
                "name": "soft_reaction_open",
                "L": 1.0,
                "R_mouth": 1.0,
                "R_exit": 0.50,
                "T_slope": 0.10,
                "V_base": 0.45,
                "V_quad": 0.10,
                "beta": 0.25,
                "Y_exit": 0.20,
                "lambda_B": 0.36,
                "lambda_U": 0.24,
                "lambda_W": 0.45,
                "lambda_R": 0.28,
            },
        ],
        "continuation": {
            "parameter": "beta",
            "values": [0.0, 0.25, 0.5, 0.75],
            "use_previous_solution_as_initial_guess": True,
        },
        "freeze_boundary": {
            "pre_target_freeze": True,
            "target_blind": True,
            "no_post_residual_refit": True,
            "packet_schema": "stage_v2_22b_solver_handoff/v1",
            "packet_output_dir": "simulation/output/packets",
        },
        "target_blind_controls": {
            "forbidden_import_prefixes": ["stage_v2_"],
            "target_evaluation_modules_imported": False,
            "target_residuals_computed": False,
            "packets_emitted_by_readiness_layer": False,
        },
    }
    protocol["protocol_hash"] = stable_hash({k: v for k, v in protocol.items() if k != "protocol_hash"})
    return protocol


def iter_candidates(protocol: Mapping[str, Any]) -> Iterable[Dict[str, Any]]:
    for idx, variant in enumerate(protocol["candidate_variants"], start=1):
        yield {
            "index": idx,
            "name": str(variant["name"]),
            "parameters": dict(variant),
        }


def candidate_freeze_payload(protocol: Mapping[str, Any], candidate: Mapping[str, Any]) -> Dict[str, Any]:
    return {
        "schema": "pde_audit_nonlinear_candidate_freeze_payload/v1",
        "protocol_hash": protocol["protocol_hash"],
        "candidate_index": int(candidate["index"]),
        "candidate_name": str(candidate["name"]),
        "parameters": dict(candidate["parameters"]),
        "newton": dict(protocol["newton"]),
        "mesh_check": dict(protocol["mesh_check"]),
        "boundary_conditions": dict(protocol["boundary_conditions"]),
        "freeze_boundary": dict(protocol["freeze_boundary"]),
    }


def candidate_freeze_hash(protocol: Mapping[str, Any], candidate: Mapping[str, Any]) -> str:
    return stable_hash(candidate_freeze_payload(protocol, candidate))
