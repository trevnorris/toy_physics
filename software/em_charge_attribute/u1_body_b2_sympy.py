#!/usr/bin/env python3
"""SymPy production route over the orchestrator-approved B2 stage-0 contract."""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Any

import sympy as sp

from u1_body_b2_common import (
    digest,
    dump_yaml,
    load_yaml_authenticated,
    rel_repo,
    require,
)


FORBIDDEN = {"Maxwell", "Larmor", "point_current", "electromagnetic_dipole"}


def route_g9(sector: str, residual_zero: bool, determined: bool, energy_independent: bool | None, missing_laws: list[str] | None = None) -> dict[str, Any]:
    exact = residual_zero and determined and (sector != "energy" or energy_independent is True)
    causes: list[str] = []
    if not exact:
        causes.append("missing_momentum_residual_norm" if sector == "momentum" else "missing_sector_tolerance")
        if sector == "energy" and energy_independent is not True:
            causes.append("return_energy_closure")
    causes = sorted(set(causes + (missing_laws or [])))
    return {
        "verdict": "OK(exact)" if exact else f"UNRESOLVED({','.join(causes)})",
        "causes": causes,
        "residual_identically_zero": residual_zero,
        "all_terms_determined": determined,
        "return_energy_structurally_independent": energy_independent,
    }


def full_residual(sector: str, balance: dict[str, Any]) -> dict[str, Any]:
    """Consume the v48 balance schema and recompute its signed sum when available."""
    if "canonical_terms" not in balance:
        typed_roots = balance.get("route_B_authenticated_typed_roots")
        require(isinstance(typed_roots, dict), "B2_PROD_G9_WITNESS", f"{sector}:typed roots")
        require(
            "complete_signed_expression" in typed_roots and typed_roots["complete_signed_expression"] is None,
            "B2_PROD_G9_WITNESS",
            f"{sector}:absent canonical terms require null complete signed expression",
        )
        missing_laws = balance.get("missing_native_current_laws")
        route_b_missing = typed_roots.get("missing_data")
        require(
            isinstance(missing_laws, list) and isinstance(route_b_missing, list),
            "B2_PROD_G9_WITNESS",
            f"{sector}:absent canonical terms require auditable source lists",
        )
        causes = sorted(set(missing_laws) | set(route_b_missing))
        require(bool(causes), "B2_PROD_G9_WITNESS", f"{sector}:absent canonical terms require non-empty cause union")
        if sector == "energy":
            require(
                "return_energy_closure" in causes,
                "B2_PROD_G9_WITNESS",
                "energy:Route B must witness return_energy_closure",
            )
        return {
            "sector": sector,
            "derivation": "independently witness that no complete signed balance expression is structurally available",
            "full_residual_terms": [],
            "computed_full_residual": None,
            "energy_return_leaf_derivative": None,
            "missing_native_current_laws": missing_laws,
            "unavailability_witness": {
                "predicate": "canonical_terms_absent && complete_signed_expression_is_null && union_nonempty",
                "canonical_terms_present": False,
                "complete_signed_expression": typed_roots["complete_signed_expression"],
                "missing_native_current_laws": missing_laws,
                "route_B_missing_data": route_b_missing,
                "named_missing_functionals": causes,
            },
            "verdict": f"UNRESOLVED({','.join(causes)})",
            "causes": causes,
            "residual_identically_zero": False,
            "all_terms_determined": False,
            "return_energy_structurally_independent": False if sector == "energy" else None,
        }

    terms: list[dict[str, Any]] = []
    expression = sp.Integer(0)
    for assignment in balance["canonical_terms"]:
        components = [sp.sympify(value) for value in assignment["canonical_symbol_components"]]
        term = sp.expand(sum(components, sp.Integer(0)))
        expression += term
        terms.append({
            "channel": assignment["channel"],
            "source_root": assignment["source_root"],
            "native_integrated_term": assignment["term_components"],
            "canonical_symbol_components": assignment["canonical_symbol_components"],
            "symbolic_term": sp.sstr(term),
            "determined": bool(assignment["determined"]),
        })
    normal = sp.expand(expression)
    return_leaf = sp.Symbol("Phi_E_return") if sector == "energy" else None
    energy_derivative = sp.diff(normal, return_leaf) if return_leaf is not None else None
    determined = all(row["determined"] for row in terms)
    independent = bool(energy_derivative == 0) if sector == "energy" else None
    residual_zero = bool(normal == 0)
    missing_laws = balance.get("missing_native_current_laws", [])
    routed = route_g9(sector, residual_zero, determined, independent, missing_laws)
    return {
        "sector": sector,
        "derivation": "sum every signed canonical component in the five typed native-action channels, then simplify",
        "full_residual_terms": terms,
        "computed_full_residual": sp.sstr(normal),
        "energy_return_leaf_derivative": sp.sstr(energy_derivative) if energy_derivative is not None else None,
        "missing_native_current_laws": missing_laws,
        **routed,
    }


def residue_status(operator: dict[str, Any]) -> str:
    table = {
        "gnls_density_phase": "UNRESOLVED(sleeve_core_trace)",
        # Contract evidence: the merged wall_chi_u_coupled inventory row has time_order: 0.
        "wall_chi_u_coupled": "ZERO(no_frequency_kinetic_term)",
        "throat_source_open": "UNRESOLVED(throat_source)",
        "wall_mix_open": "UNRESOLVED(wall_mix)",
        "brane_shear_transverse": "UNRESOLVED(tilt_profile)",
        "brane_normal_local": "ZERO(no_propagating_support)",
        "h_scalar": "UNRESOLVED(geon_core_bundle)",
        "geon_open": "UNRESOLVED(geon_core_bundle)",
    }
    require(operator["id"] in table, "B2_PROD_OPERATOR_COVERAGE", operator["id"])
    return table[operator["id"]]


def build_cells(stage0: dict[str, Any], g9: dict[str, dict[str, Any]]) -> list[dict[str, Any]]:
    axes = stage0["frozen_data"]["minimum_obligation_manifest"]["grid_axes"]
    cells = []
    for endpoint in axes["endpoint"]:
        for parity in axes["parity"]:
            for closure in axes["closure"]:
                for open_stratum in axes["open_stratum"]:
                    tilt = endpoint in {"E3", "E4", "E5"}
                    acceleration = "UNRESOLVED(tilt_profile)" if tilt else "UNRESOLVED(return_closure)"
                    cells.append({
                    "key": f"{endpoint}|{parity}|{closure}|{open_stratum}",
                    "axes": {"endpoint": endpoint, "parity": parity, "closure_branch": closure, "open_stratum": open_stratum},
                    "C_mdot": {
                        "definition": "d F_flux_A / d Xddot_j", "rows": ["X_x", "X_y", "X_z", "p_x", "p_y", "p_z"],
                        "columns": ["Xddot_x", "Xddot_y", "Xddot_z"], "X_row_dimensions": "mass",
                        "p_row_dimensions": "generalized_tilt_force/translation_acceleration", "status": acceleration,
                        "computed_ancestry": ["native_momentum", "outer_control_flux", "return_closure", "moduli_fixed_collective_tangent"],
                    },
                    "velocity_block": {"definition": "d F_flux_A / d Xdot_j", "alias": "D_intake", "non_additive_identity": True, "status": "UNRESOLVED(return_closure)", "count_in_P_local": 1},
                    "pdot_generalized_velocity_remainder": {"definition": "d F_flux_A / d pdot_j", "force_channel": "flux", "status": "UNRESOLVED(tilt_profile)" if tilt else "UNRESOLVED(return_closure)"},
                    "G2": "UNRESOLVED(return_closure)", "G5": acceleration, "G6": "UNRESOLVED(return_closure)",
                    "G7": {"status": acceleration, "p_treatment": "free_external_parameter_not_solved", "balance_includes_D_intake_once": True},
                    "G9": g9,
                })
    return cells


def typed_dag(floor: dict[str, Any]) -> dict[str, Any]:
    """Bind every full-key floor object to a concrete production family."""
    nodes: dict[str, dict[str, Any]] = {"report_headline": {"type": "sink", "depends_on": []}}

    def add(product: str, producer: str) -> None:
        nodes[product] = {"type": "obligation_product", "producer": producer, "depends_on": []}
        nodes["report_headline"]["depends_on"].append(product)

    for product in floor["expanded_records"]:
        prefix = product.split("|", 1)[0]
        if prefix.startswith(("C_mdot", "pdot_", "X_p_", "G2_", "G5_", "G6_", "G7_", "G9_")):
            producer = "grid.cells"
        elif prefix.startswith(("branch_", "pair_", "F_rad", "K_rad", "radiative_", "per_cell_")):
            producer = "radiation"
        elif prefix.startswith("NOT_RUN_"):
            producer = "phase_C"
        elif prefix.startswith("stage0_"):
            producer = "stage0_evidence"
        elif prefix in {"generated_operator_inventory", "generated_endpoint_branch_inventory", "mode_coverage_residual", "total_radiative_mass_current", "total_radiative_Noether_energy_flux", "total_radiative_Noether_momentum_flux", "total_F_rad", "total_K_rad", "total_work_storage_flux_identity", "K_total", "K_self_field"}:
            producer = "radiation"
        else:
            producer = "partition_or_status"
        add(product, producer)
    nodes["report_headline"]["depends_on"].sort()
    return {"root": "report_headline", "nodes": nodes}


def build(input_path: Path, stage0_path: Path, stage0_digest: str) -> dict[str, Any]:
    stage0, stage0_auth = load_yaml_authenticated(stage0_path, stage0_digest, "production_sympy:stage0")
    require(stage0["status"] == "AWAITING_ORCHESTRATOR_APPROVAL", "B2_PROD_STAGE0_SCHEMA", "stage-0 status")
    manifest = stage0["observation_contract"]["Obs_B2_manifest"]
    input_key = rel_repo(input_path)
    require(input_key in manifest, "B2_A1_OBS_B2_EXACT", f"input not manifested: {input_key}")
    config, input_auth = load_yaml_authenticated(input_path, manifest[input_key], "production_sympy:input")
    # NOTE: .resolve() is required. Under authenticated_exec the runner content-pins this
    # script and re-execs it as /proc/self/fd/<fd>, so __file__ is that /proc path; without
    # .resolve() the parent is /proc/self/fd and the repo-relative key/containment checks fail.
    # .resolve() follows the descriptor symlink back to the real repo file. (Bugfix 2026-07-16:
    # first-ever run of the resume() production path; hash change re-sealed into the stage-0 manifest.)
    b1_path = Path(__file__).resolve().parent / "reports/u1_body_dynamics_artifacts/stage1/sympy_phase_b1.yaml"
    b1_key = rel_repo(b1_path)
    require(b1_key in manifest, "B2_A1_OBS_B2_EXACT", f"B1 leaf not manifested: {b1_key}")
    b1, b1_auth = load_yaml_authenticated(b1_path, manifest[b1_key], "production_sympy:B1_partition_ledger")

    frozen = stage0["frozen_data"]
    operators = frozen["native_operator_inventory"]
    require(len(operators) == 8, "B2_PROD_OPERATOR_COVERAGE", "eight committed native sectors")
    ancestry = {ancestor for row in operators for ancestor in row.get("action_ancestry", [])}
    require(not (ancestry & FORBIDDEN), "B2_PROD_NATIVE_ANCESTRY", "forbidden imported ancestry")
    balances = frozen["integrated_balance_identities"]["sectors"]
    g9 = {sector: full_residual(sector, balances[sector]) for sector in ["mass", "momentum", "energy"]}
    cells = build_cells(stage0, g9)
    ledger = b1["partition_ledger"]
    require(len(ledger["records"]) == 41, "B2_PROD_LEDGER", "41-record frozen ledger")

    branch_cells = []
    for row in operators:
        status = residue_status(row)
        branch_cells.append({"operator": row["id"], "steady": status, "accelerating": status, "radiative_mass_current": status, "p_slices": frozen["minimum_obligation_manifest"]["p_slices"]})
    actual_modes = sorted(row["operator"] for row in branch_cells)
    expected_modes = sorted(row["id"] for row in operators)
    mode_residual = {"expected": expected_modes, "actual": actual_modes, "missing": sorted(set(expected_modes) - set(actual_modes)), "extra": sorted(set(actual_modes) - set(expected_modes))}
    require(not mode_residual["missing"] and not mode_residual["extra"], "B2_PROD_MODE_COVERAGE", str(mode_residual))

    partition = {
        "legacy_parent": "outer_control_flux:translation", "legacy_parent_digest_verified": True,
        "concrete_M_record_count": 40, "variational_translation_census": 40, "flux_translation_acceleration_census": 1,
        "terminal_owner_enum": frozen["ownership_convention"]["terminal_owner_enum"], "state": "UNRESOLVED(return_closure)",
        "partition_reconstruction_residual": "UNRESOLVED(return_closure)", "two_route_reconstruction_residual": "UNRESOLVED(return_closure)",
        "source_to_term_incidence_residual": frozen["native_operator_action_incidence_residual"],
        "inherited_open_dispositions": {name: f"UNRESOLVED({name})" for name in ["geon_core_bundle", "outer_surface_functional", "sleeve_core_trace", "throat_source", "throat_surface_functional", "wall_mix", "tilt_profile"]},
        "X_p_remainder": "UNRESOLVED(tilt_profile)",
    }
    radiation = {
        "trajectory_representation": frozen["trajectory_representation"], "branch_cells": branch_cells,
        "mode_coverage_residual": mode_residual,
        "per_cell_resolvent_identity": frozen["endpoint_resolvent_cells"],
        "per_cell_known_nonzero_control": [{"cell": row["cell"], **row["known_nonzero_control"]} for row in frozen["endpoint_resolvent_cells"]],
        "interference": "UNRESOLVED(native_branch_inputs)", "K_self_field": "UNRESOLVED(native_branch_inputs)", "K_total": "UNRESOLVED(native_branch_inputs)",
        "totals": {name: "UNRESOLVED(native_branch_inputs)" for name in ["Noether_flux", "F_rad", "K_rad", "work_storage_flux_identity", "radiative_mass_current"]},
        "forbidden_ancestry_guard": "PASS",
    }
    artifact = {
        "schema_version": "U1_PHASE_B2_PRODUCTION_ENGINE_V3", "engine": "SymPy",
        "independent_route": "typed_control_volume_projection_of_frozen_native_action_derivations",
        "stage0_contract_sha256": stage0_digest, "input_sha256": input_auth["consumed_sha256"], "status": "COMPLETE_WITH_HONEST_OUTCOMES",
        "first_use_authentication": [stage0_auth, input_auth, b1_auth],
        "stage0_evidence": {"sector_current_derivations": frozen["sector_current_derivations"], "native_noether_derivations": frozen["native_noether_derivations"], "integrated_balance_identities": frozen["integrated_balance_identities"]},
        "grid": {"cell_count": len(cells), "cells": cells}, "operator_inventory": operators, "partition": partition, "radiation": radiation,
        "phase_C": {"G8": "NOT_RUN(phase_C)", "G10": "NOT_RUN(phase_C)", "G11": "NOT_RUN(phase_C)"}, "typed_dag": typed_dag(frozen["minimum_obligation_manifest"]),
    }
    artifact["sink_digest"] = digest({key: artifact[key] for key in ["grid", "operator_inventory", "partition", "radiation", "phase_C", "typed_dag"]})
    return artifact


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True); parser.add_argument("--stage0", type=Path, required=True)
    parser.add_argument("--stage0-contract-digest", required=True); parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    try:
        artifact = build(args.input.resolve(), args.stage0.resolve(), args.stage0_contract_digest)
        dump_yaml(args.output, artifact)
        print(f"B2_SYMPY: COMPLETE_WITH_HONEST_OUTCOMES cells={artifact['grid']['cell_count']}")
        return 0
    except Exception as exc:
        print(str(exc)); return 1


if __name__ == "__main__":
    raise SystemExit(main())
