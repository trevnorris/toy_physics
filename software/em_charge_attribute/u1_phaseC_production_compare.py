#!/usr/bin/env python3
"""Comparator and report assembler for the U1 Phase-C production engines."""

from __future__ import annotations

import argparse
import copy
import hashlib
import json
import sys
from pathlib import Path
from typing import Any, Callable

import yaml


SCHEMA = "U1_PHASE_C_PRODUCTION_ENGINE_V1"
RESULT_SCHEMA = "U1_PHASE_C_PRODUCTION_RESULTS_V1"
RATIFIED_DIGEST = "83233baabd7f8e27c88d130b911691e76d01d5797da8eeb32c90bbae111ec95a"
MEDIATORS = ("h", "u_T", "u_L", "wall_chi")
ENDPOINTS = ("E1", "E2", "E3", "E4", "E5")
AMBIENTS = ("one_sided_pathA29", "symmetric_postulate")
CLOSURES = ("body_mass_growth", "return_path", "sleeve_exit")
PARITY_TAGS = {"body-only", "ambient-postulate-dependent", "one-sided-asymmetry-map"}
TILT_ENUM = {
    "TILT_LINEAR", "TILT_OTHER", "TILT_ZERO", "TILT_NO_STEADY",
    "TILT_UNSTABLE", "TILT_UNRESOLVED",
}
COUPLING_ENUM = {
    "EXACT_SV", "SV_PLUS_DEPARTURE", "DEPARTURE_ONLY", "NULL",
    "UNRESOLVED", "ILL_POSED",
}
CHANNELS = {"variational", "flux", "constraint/multiplier", "Rayleigh", "radiation"}


class AssertionDeath(RuntimeError):
    def __init__(self, assert_id: str, detail: str):
        super().__init__(detail)
        self.assert_id = assert_id
        self.detail = detail


CHECK_ORDER: list[str] = []


def checked(condition: bool, assert_id: str, detail: str, checks: list[dict[str, str]]) -> None:
    if assert_id not in CHECK_ORDER:
        CHECK_ORDER.append(assert_id)
    if not condition:
        raise AssertionDeath(assert_id, detail)
    checks.append({"assert_id": assert_id, "status": "PASS"})


def load_yaml(path: Path) -> Any:
    if path.suffix == ".json":
        return json.loads(path.read_text(encoding="utf-8"))
    with path.open("rb") as handle:
        return yaml.load(handle, Loader=yaml.CSafeLoader)


def dump_yaml(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        yaml.safe_dump(value, sort_keys=False, allow_unicode=True, width=140),
        encoding="utf-8",
    )


def canonical_bytes(value: Any) -> bytes:
    return json.dumps(
        value, ensure_ascii=False, sort_keys=True, separators=(",", ":")
    ).encode("utf-8")


def digest(value: Any) -> str:
    return hashlib.sha256(canonical_bytes(value)).hexdigest()


def ids(values: list[str] | tuple[str, ...] | set[str]) -> list[str]:
    return sorted(set(values), key=lambda value: (value.lower(), value))


def rows_by(rows: list[dict[str, Any]], key: str) -> dict[str, dict[str, Any]]:
    return {row[key]: row for row in rows}


def enum_valid(record: dict[str, Any], allowed: set[str]) -> bool:
    return record.get("enum") in allowed


def normalize_ablation(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    return sorted(
        [{name: row[name] for name in (
            "term_id", "root_id", "response_nonzero", "support", "root_type"
        )} for row in rows],
        key=lambda row: row["term_id"],
    )


def normalize_force(engine: dict[str, Any]) -> list[dict[str, Any]]:
    out = []
    for endpoint in engine["tilt"]["force_balance_by_endpoint"]:
        out.append({
            "endpoint": endpoint["endpoint"],
            "E1_placement": endpoint["E1_placement"],
            "template": ids(endpoint["constructed_template_terms"]),
            "active": ids(endpoint["active_residual_terms"]),
            "channel_terms": {name: ids(values) for name, values in endpoint["channel_terms"].items()},
            "channel_sum_reconstructs": endpoint["channel_sum_reconstructs_active_residual"],
            "no_double_count": endpoint["no_double_count"],
            "terms": sorted([
                {name: row[name] for name in (
                    "term_id", "source_id", "channel", "support",
                    "formal_expression", "applicability", "dimensions_restored"
                )}
                for row in endpoint["terms"]
            ], key=lambda row: row["term_id"]),
        })
    return sorted(out, key=lambda row: row["endpoint"])


def normalize_tilt_cells(engine: dict[str, Any]) -> list[dict[str, Any]]:
    return sorted([
        {
            "cell_id": row["cell_id"],
            "key": row["key"],
            "formalism_id": row["formalism_id"],
            "availability": row["availability"],
            "physics_status": row["physics_status"],
            "ancestry": ids(row["computed_typed_ancestry_unresolved_slots"]),
            "mapped": ids(row["dependency_map_slots"]),
            "parity": row["parity"],
            "steady_substitution": row["steady_substitution"],
            "susceptibility_status": row["susceptibility_status"],
        }
        for row in engine["tilt"]["cells"]
    ], key=lambda row: row["cell_id"])


def normalize_delta(engine: dict[str, Any]) -> list[dict[str, Any]]:
    return sorted([
        {
            "row": row["row_mediator"], "column": row["column_mediator"],
            "diagonal": row["diagonal"], "functional": row["functional"],
            "availability": row["availability"], "classification": row["classification"],
        }
        for row in engine["coupling_package"]["7.5b"]["deltaO_AB"]
    ], key=lambda row: (row["row"], row["column"]))


def normalize_coupling_cells(engine: dict[str, Any]) -> list[dict[str, Any]]:
    return sorted([
        {
            "cell_id": row["cell_id"], "mediator": row["mediator"], "key": row["key"],
            "formalism_id": row["formalism_id"], "availability": row["availability"],
            "off_shell": row["off_shell_in_p_status"], "physics_status": row["physics_status"],
            "ancestry": ids(row["computed_typed_ancestry_unresolved_slots"]),
            "mapped": ids(row["dependency_map_slots"]), "s_parity": row["s_parity"],
            "O(V)": row["O(V)_classification"], "j_sV": row["j_proportional_sV"],
            "split": row["mass_charge_split_status"], "steady": row["steady_substituted_row"],
        }
        for row in engine["coupling_package"]["7.5d"]["cells"]
    ], key=lambda row: row["cell_id"])


def normalize_ownership(engine: dict[str, Any]) -> list[dict[str, Any]]:
    out = []
    for endpoint in engine["coupling_package"]["channel_ownership_by_endpoint"]:
        out.append({
            "endpoint": endpoint["endpoint"], "E1_placement": endpoint["E1_placement"],
            "channel_terms": {name: ids(values) for name, values in endpoint["channel_terms"].items()},
            "channel_sum_reconstructs": endpoint["channel_sum_reconstructs_active_response"],
            "no_double_count": endpoint["no_double_count"],
            "entries": sorted(endpoint["entries"], key=lambda row: row["entry_id"]),
        })
    return sorted(out, key=lambda row: row["endpoint"])


def normalize_gates(engine: dict[str, Any]) -> dict[str, Any]:
    gates = copy.deepcopy(engine["gates"])
    gates["G8"]["records"] = sorted(gates["G8"]["records"], key=lambda row: row["source_id"])
    return gates


def semantic_view(engine: dict[str, Any]) -> dict[str, Any]:
    formalism = engine["tilt"]["formalism"]
    dimension = engine["dimensional_firewall"]
    return {
        "axes": engine["axes"],
        "availability": {
            "summary": engine["availability_contract"]["ratified_summary"],
            "unresolved": ids(engine["availability_contract"]["unresolved_slot_ids"]),
            "derived": ids(engine["availability_contract"]["derived_slot_ids"]),
            "witness_challenge_references": sorted(
                engine["availability_contract"]["witness_challenge_pairs_consumed_by_reference"],
                key=lambda row: row["slot_id"],
            ),
        },
        "native_action": {
            "S_body": engine["native_action"]["S_body"],
            "term_ids": ids(engine["native_action"]["action_term_ids"]),
            "ablations": normalize_ablation(engine["native_action"]["per_native_term_ablation"]),
        },
        "tilt_formulas": {
            "profile": formalism["profile_family"],
            "residual": formalism["field_equation_residual"],
            "force": formalism["total_force_balance"],
            "statics": formalism["statics"],
            "susceptibility": formalism["susceptibility"],
            "partition": formalism["partition_successor"],
            "enum": formalism["parent_enum"],
        },
        "force": normalize_force(engine),
        "tilt_cells": normalize_tilt_cells(engine),
        "7.5a": engine["coupling_package"]["7.5a"],
        "J": sorted(engine["coupling_package"]["7.5b"]["J_A"], key=lambda row: row["mediator"]),
        "deltaO": normalize_delta(engine),
        "7.5c": sorted(engine["coupling_package"]["7.5c"], key=lambda row: row["endpoint"]),
        "total_response": engine["coupling_package"]["7.5d"]["formal_total_response"],
        "coupling_cells": normalize_coupling_cells(engine),
        "ownership": normalize_ownership(engine),
        "parity_census": sorted(engine["coupling_package"]["parity_census"], key=lambda row: (
            row["endpoint"], row["ambient_branch"], row["field_or_profile"]
        )),
        "mouth": sorted(engine["coupling_package"]["mouth_datum_records"], key=lambda row: row["endpoint"]),
        "gates": normalize_gates(engine),
        "G4": {
            "claimed": ids(engine["G4_known_nonzero_controls"]["claimed_zero_ids"]),
            "covered": ids(engine["G4_known_nonzero_controls"]["covered_zero_ids"]),
            "control_classes": sorted([
                {"control_class": row["control_class"], "covers": ids(row["covers_zero_ids"]),
                 "fixture_nonzero": row["fixture_nonzero"],
                 "dimensions_restored": row["dimensions_restored"]}
                for row in engine["G4_known_nonzero_controls"]["controls"]
            ], key=lambda row: row["control_class"]),
        },
        "reconciliation": engine["reconciliation"],
        "dimension": {
            "classes": sorted([
                {"expression_class": row["expression_class"], "homogeneous": row["homogeneous"]}
                for row in dimension["constructed_expression_classes"]
            ], key=lambda row: row["expression_class"]),
            "all_inline_homogeneous": dimension["all_inline_homogeneous"],
            "cross_expression_consistent": dimension["cross_expression_consistent"],
            "no_back_solved_carrier": dimension["no_back_solved_carrier"],
            "back_solved": dimension["back_solved_free_carriers"],
            "firing_ablation": dimension["firing_ablation"]["heterogeneity_detected"],
        },
        "ancestry": engine["native_ancestry_guard"],
        "symbolic_inventory": engine["symbolic_object_inventory"],
        "A9": {
            "ids": ids(engine["A9_external_verification"]["object_ids"]),
            "map": sorted(engine["A9_external_verification"]["coverage_map"], key=lambda row: row["object_id"]),
            "category_counts": engine["A9_external_verification"]["coverage_category_counts"],
        },
        "headline_entries": engine["headline_entries"],
        "agreement_nonce": engine.get("agreement_nonce", "baseline"),
    }


def apply_mutation(engine: dict[str, Any], mutation: str, side: str) -> None:
    """Plant one comparator-visible defect without changing engine code paths."""
    if mutation == "ASSERT_ENGINE_SCHEMA":
        engine["schema_version"] = "MUTATED"
    elif mutation == "ASSERT_STAGE0_BINDING":
        engine["stage0_binding"]["equal"] = False
    elif mutation == "ASSERT_AXES":
        engine["axes"]["tilt_cell_count"] -= 1
    elif mutation == "ASSERT_STRATA":
        engine["axes"]["axis_strata_exact_set_equal"] = False
    elif mutation == "ASSERT_AVAILABILITY":
        engine["availability_contract"]["ratified_summary"]["UNRESOLVED"] -= 1
    elif mutation == "ASSERT_ACTION_NATIVE_ABLATIONS":
        engine["native_action"]["per_native_term_ablation"][0]["response_nonzero"] = False
    elif mutation == "ASSERT_FORCE_CENSUS":
        engine["tilt"]["force_census_incidence_complete"] = False
    elif mutation == "ASSERT_FORCE_CHANNEL_OWNERSHIP":
        engine["tilt"]["force_balance_by_endpoint"][0]["terms"][0]["channel"] = "not-a-channel"
    elif mutation == "ASSERT_FORCE_NO_DOUBLE_COUNT":
        engine["tilt"]["force_balance_by_endpoint"][0]["no_double_count"] = False
    elif mutation == "ASSERT_ENDPOINT_PLACEMENT":
        engine["tilt"]["force_balance_by_endpoint"][0]["E1_placement"] = "multiplier"
    elif mutation == "ASSERT_TILT_FORMAL_FUNCTIONALS":
        engine["tilt"]["formalism"]["profile_family"]["expression"] = ""
    elif mutation == "ASSERT_TILT_GRID":
        engine["tilt"]["cells"].pop()
    elif mutation == "ASSERT_TILT_DEPENDENCY_MAP":
        engine["tilt"]["cells"][0]["dependency_map_slots"].pop()
    elif mutation == "ASSERT_TILT_STATUS_ENUM":
        engine["tilt"]["cells"][0]["physics_status"]["enum"] = "OK"
    elif mutation == "ASSERT_PARITY_AUTHORITY":
        engine["tilt"]["cells"][0]["parity"]["authority_tag"] = "untagged"
    elif mutation == "ASSERT_SUSCEPTIBILITY_BRANCHES":
        engine["tilt"]["formalism"]["susceptibility"]["branches"] = []
    elif mutation == "ASSERT_COUPLING_75A":
        engine["coupling_package"]["7.5a"]["Sigma_surface"]["expression"] = ""
    elif mutation == "ASSERT_J_MEDIATORS":
        engine["coupling_package"]["7.5b"]["J_A"].pop()
    elif mutation == "ASSERT_DELTAO_ORDERED":
        engine["coupling_package"]["7.5b"]["deltaO_AB"].pop()
    elif mutation == "ASSERT_COUPLING_CENSUS":
        engine["coupling_package"]["coupling_census_incidence_complete"] = False
    elif mutation == "ASSERT_COUPLING_CHANNEL_OWNERSHIP":
        engine["coupling_package"]["channel_ownership_by_endpoint"][0]["exactly_one_channel"] = False
    elif mutation == "ASSERT_COUPLING_GRID":
        engine["coupling_package"]["7.5d"]["cells"].pop()
    elif mutation == "ASSERT_COUPLING_DEPENDENCY_MAP":
        engine["coupling_package"]["7.5d"]["cells"][0]["dependency_map_slots"].pop()
    elif mutation == "ASSERT_TOTAL_COUPLED_RESPONSE":
        engine["coupling_package"]["7.5d"]["formal_total_response"]["full_mixed_kernel_included"] = False
    elif mutation == "ASSERT_PROJECTION_FROZEN":
        engine["coupling_package"]["7.5d"]["formal_total_response"]["fixed_projection_id"] = "post_hoc"
    elif mutation == "ASSERT_MASS_CHARGE_SPLIT":
        del engine["coupling_package"]["mass_charge_split"]["orientation_charge"]
    elif mutation == "ASSERT_COUPLING_CLASSIFICATION":
        engine["coupling_package"]["7.5d"]["cells"][0]["physics_status"]["enum"] = "MAGNETIC"
    elif mutation == "ASSERT_PARITY_CENSUS":
        engine["coupling_package"]["parity_census"].pop()
    elif mutation == "ASSERT_MOUTH_RECORDS":
        engine["coupling_package"]["mouth_datum_records"].pop()
    elif mutation == "ASSERT_NATIVE_ANCESTRY":
        engine["native_ancestry_guard"]["forbidden_nodes"] = ["import:Maxwell_exchange_kernel"]
    elif mutation == "ASSERT_NO_BANNED_TILT":
        engine["symbolic_object_inventory"]["signed_tilt_coordinate_constructed"] = True
    elif mutation == "ASSERT_ONE_BODY_SCOPE":
        engine["symbolic_object_inventory"]["two_body_objects_constructed"] = ["force_sign"]
    elif mutation == "ASSERT_G8":
        engine["gates"]["G8"]["level2_exactly_one"] = False
    elif mutation == "ASSERT_G11":
        engine["gates"]["G11"]["velocity_fixture_ablation_zero"] = False
    elif mutation == "ASSERT_G10":
        engine["gates"]["G10"]["zero_would_be_recorded_honestly"] = False
    elif mutation == "ASSERT_STANDING_GATES":
        engine["gates"]["G5"]["ambient_quarantine_enforced"] = False
    elif mutation == "ASSERT_G4_ZERO_CONTROLS":
        engine["G4_known_nonzero_controls"]["exact_coverage"] = False
    elif mutation == "ASSERT_RECONCILIATION":
        engine["reconciliation"]["exact_set_equal"] = False
    elif mutation == "ASSERT_PARTITION_SUCCESSOR":
        engine["tilt"]["formalism"]["partition_successor"]["upstream_fact_reused_as_closed_owner"] = True
    elif mutation == "ASSERT_DIMENSIONAL_FIREWALL":
        engine["dimensional_firewall"]["all_inline_homogeneous"] = False
    elif mutation == "ASSERT_UNRESOLVED_HONESTY":
        engine["coupling_package"]["7.5d"]["cells"][0]["j_proportional_sV"] = "EXACT_SV"
    elif mutation == "ASSERT_A9_COVERAGE":
        engine["A9_external_verification"]["object_count"] -= 1
    elif mutation == "ASSERT_ENGINE_INDEPENDENCE":
        engine["independent_route"] = "same-route"
    elif mutation == "ASSERT_DUAL_ENGINE_AGREEMENT":
        engine["agreement_nonce"] = f"mutated-{side}"
    elif mutation == "ASSERT_HEADLINES":
        engine["headline_entries"].pop()
    elif mutation in {"ASSERT_SUCCESSOR_ASSEMBLY", "ASSERT_SUMMARY_COMPLETE"}:
        # These defects are planted in their post-comparison assembly paths.
        pass
    else:
        print("MUTATION_NOOP", mutation, file=sys.stderr)


def assert_engine(
    sympy: dict[str, Any],
    wolfram: dict[str, Any],
    projection: dict[str, Any],
) -> tuple[list[dict[str, str]], dict[str, Any]]:
    checks: list[dict[str, str]] = []
    engines = (sympy, wolfram)
    checked(all(row["schema_version"] == SCHEMA for row in engines),
            "ASSERT_ENGINE_SCHEMA", "engine schema mismatch", checks)
    checked(all(
        row["stage0_binding"]["equal"]
        and row["stage0_binding"]["supplied_digest"] == RATIFIED_DIGEST
        for row in engines
    ), "ASSERT_STAGE0_BINDING", "stage-0 digest not bound in both engines", checks)
    checked(all(
        row["axes"]["endpoints"] == list(ENDPOINTS)
        and row["axes"]["ambient_branches"] == list(AMBIENTS)
        and row["axes"]["closure_branches"] == list(CLOSURES)
        and row["axes"]["tilt_cell_count"] == 240
        and row["axes"]["coupling_cell_count"] == 960
        and not row["axes"]["axis_collapse_performed"]
        for row in engines
    ), "ASSERT_AXES", "full production axes/cell counts missing", checks)
    checked(all(
        row["axes"]["axis_strata_exact_set_equal"]
        and len(row["axes"]["GAP_OPEN_FIELD_PROFILE_strata"]) == 8
        and ids(row["axes"]["GAP_OPEN_FIELD_PROFILE_strata"])
        == ids(row["axes"]["generated_active_strata"])
        for row in engines
    ), "ASSERT_STRATA", "OPEN-stratum membership is not exact", checks)
    checked(all(
        row["availability_contract"]["ratified_summary"]["DERIVED"] == 2
        and row["availability_contract"]["ratified_summary"]["UNRESOLVED"] == 56
        and len(row["availability_contract"]["unresolved_slot_ids"]) == 56
        and ids(item["slot_id"] for item in row["availability_contract"]["witness_challenge_pairs_consumed_by_reference"])
        == ids(row["availability_contract"]["unresolved_slot_ids"])
        and all(item["witness_id"] and item["challenge_id"]
                for item in row["availability_contract"]["witness_challenge_pairs_consumed_by_reference"])
        for row in engines
    ), "ASSERT_AVAILABILITY", "ratified slot dispositions changed", checks)
    checked(all(
        row["native_action"]["all_nonvanishing"]
        and len(row["native_action"]["per_native_term_ablation"]) == 15
        and all(item["response_nonzero"] for item in row["native_action"]["per_native_term_ablation"])
        for row in engines
    ), "ASSERT_ACTION_NATIVE_ABLATIONS", "native action ablation did not respond", checks)
    checked(all(
        row["tilt"]["force_census_incidence_complete"]
        and len(row["tilt"]["force_census_expected_terms"]) == 28
        for row in engines
    ), "ASSERT_FORCE_CENSUS", "force census incidence incomplete", checks)
    checked(all(
        all(term["channel"] in CHANNELS for endpoint in row["tilt"]["force_balance_by_endpoint"] for term in endpoint["terms"])
        for row in engines
    ), "ASSERT_FORCE_CHANNEL_OWNERSHIP", "force term lacks exactly-one legal channel", checks)
    checked(all(
        all(endpoint["no_double_count"] and endpoint["channel_sum_reconstructs_active_residual"]
            for endpoint in row["tilt"]["force_balance_by_endpoint"])
        for row in engines
    ), "ASSERT_FORCE_NO_DOUBLE_COUNT", "force residual double-count guard fired", checks)
    checked(all(
        rows_by(row["tilt"]["force_balance_by_endpoint"], "endpoint")["E1"]["E1_placement"]
        == "variational holonomic field boundary condition"
        for row in engines
    ), "ASSERT_ENDPOINT_PLACEMENT", "E1 placement is not the frozen holonomic variational route", checks)
    checked(all(
        all(row["tilt"]["formalism"][name]["expression"] for name in (
            "profile_family", "field_equation_residual", "total_force_balance"
        ))
        for row in engines
    ), "ASSERT_TILT_FORMAL_FUNCTIONALS", "tilt functional/residual is empty", checks)
    checked(all(len(row["tilt"]["cells"]) == 240 for row in engines),
            "ASSERT_TILT_GRID", "tilt grid is not full", checks)
    checked(all(
        all(cell["dependency_exact_set_equal"]
            and ids(cell["dependency_map_slots"]) == ids(cell["computed_typed_ancestry_unresolved_slots"])
            for cell in row["tilt"]["cells"])
        for row in engines
    ), "ASSERT_TILT_DEPENDENCY_MAP", "tilt dependency map is not exact", checks)
    checked(all(
        set(row["tilt"]["formalism"]["parent_enum"]) == TILT_ENUM
        and all(enum_valid(cell["physics_status"], TILT_ENUM) for cell in row["tilt"]["cells"])
        for row in engines
    ), "ASSERT_TILT_STATUS_ENUM", "tilt parent enum/exclusivity violated", checks)
    checked(all(
        all(cell["parity"]["authority_tag"] in PARITY_TAGS for cell in row["tilt"]["cells"])
        for row in engines
    ), "ASSERT_PARITY_AUTHORITY", "parity authority tag missing", checks)
    checked(all(
        row["tilt"]["formalism"]["susceptibility"]["same_residual_linearized"]
        and len(row["tilt"]["formalism"]["susceptibility"]["branches"]) == 4
        and row["tilt"]["formalism"]["susceptibility"]["pole_condition"]
        for row in engines
    ), "ASSERT_SUSCEPTIBILITY_BRANCHES", "susceptibility pole/branch package incomplete", checks)
    checked(all(
        row["coupling_package"]["7.5a"]["domain"]["expression"]
        and row["coupling_package"]["7.5a"]["Sigma_surface"]["expression"]
        and row["coupling_package"]["7.5a"]["partial_Omega_surface"]["expression"]
        for row in engines
    ), "ASSERT_COUPLING_75A", "7.5a domain/surface functional missing", checks)
    checked(all(
        {item["mediator"] for item in row["coupling_package"]["7.5b"]["J_A"]} == set(MEDIATORS)
        and all(item["functional"]["expression"] for item in row["coupling_package"]["7.5b"]["J_A"])
        for row in engines
    ), "ASSERT_J_MEDIATORS", "J_A mediator set is incomplete", checks)
    checked(all(
        {(item["row_mediator"], item["column_mediator"])
         for item in row["coupling_package"]["7.5b"]["deltaO_AB"]}
        == {(left, right) for left in MEDIATORS for right in MEDIATORS}
        for row in engines
    ), "ASSERT_DELTAO_ORDERED", "full ordered deltaO matrix is incomplete", checks)
    checked(all(row["coupling_package"]["coupling_census_incidence_complete"] for row in engines),
            "ASSERT_COUPLING_CENSUS", "coupling census incidence incomplete", checks)
    checked(all(
        all(endpoint["exactly_one_channel"] and endpoint["expected_reachable_exact_set_equal"]
            and endpoint["channel_sum_reconstructs_active_response"] and endpoint["no_double_count"]
            for endpoint in row["coupling_package"]["channel_ownership_by_endpoint"])
        for row in engines
    ), "ASSERT_COUPLING_CHANNEL_OWNERSHIP", "coupling channel ownership/reachability failed", checks)
    checked(all(len(row["coupling_package"]["7.5d"]["cells"]) == 960 for row in engines),
            "ASSERT_COUPLING_GRID", "coupling grid is not full", checks)
    checked(all(
        all(cell["dependency_exact_set_equal"]
            and ids(cell["dependency_map_slots"]) == ids(cell["computed_typed_ancestry_unresolved_slots"])
            for cell in row["coupling_package"]["7.5d"]["cells"])
        for row in engines
    ), "ASSERT_COUPLING_DEPENDENCY_MAP", "coupling dependency map is not exact", checks)
    checked(all(
        row["coupling_package"]["7.5d"]["formal_total_response"]["total_not_source_only"]
        and row["coupling_package"]["7.5d"]["formal_total_response"]["full_mixed_kernel_included"]
        and row["coupling_package"]["7.5d"]["formal_total_response"]["multipole"]
        for row in engines
    ), "ASSERT_TOTAL_COUPLED_RESPONSE", "7.5d is not the total source+mixed-kernel response", checks)
    checked(all(
        row["coupling_package"]["7.5d"]["formal_total_response"]["fixed_projection_id"] == projection["id"]
        and row["coupling_package"]["7.5d"]["formal_total_response"]["projection"] == projection["projection"]
        for row in engines
    ), "ASSERT_PROJECTION_FROZEN", "c_sv/Delta projection differs from ratified artifact", checks)
    checked(all(
        set(row["coupling_package"]["mass_charge_split"]) >= {"mass_drain", "orientation_charge", "total", "V0_orientation_label"}
        for row in engines
    ), "ASSERT_MASS_CHARGE_SPLIT", "mass/charge split incomplete", checks)
    checked(all(
        set(row["coupling_package"]["parent_enum"]) == COUPLING_ENUM
        and all(enum_valid(cell["physics_status"], COUPLING_ENUM)
                for cell in row["coupling_package"]["7.5d"]["cells"])
        for row in engines
    ), "ASSERT_COUPLING_CLASSIFICATION", "coupling enum/exclusivity violated", checks)
    checked(all(
        len(row["coupling_package"]["parity_census"]) == 40
        and all(item["authority_tag"] in PARITY_TAGS for item in row["coupling_package"]["parity_census"])
        for row in engines
    ), "ASSERT_PARITY_CENSUS", "one-body parity census incomplete", checks)
    checked(all(
        {item["endpoint"] for item in row["coupling_package"]["mouth_datum_records"]} == set(ENDPOINTS)
        and all("fixed_source_vs_displacement_or_defect" in item for item in row["coupling_package"]["mouth_datum_records"])
        for row in engines
    ), "ASSERT_MOUTH_RECORDS", "descriptive mouth-datum records incomplete", checks)
    checked(all(
        row["native_ancestry_guard"]["guard_status"] == "PASS"
        and row["native_ancestry_guard"]["runtime_scan_executed"]
        and not row["native_ancestry_guard"]["forbidden_pattern_hits"]
        and not row["native_ancestry_guard"]["forbidden_nodes"]
        and not row["native_ancestry_guard"]["Maxwell_Larmor_Coulomb_ancestry_nodes"]
        for row in engines
    ), "ASSERT_NATIVE_ANCESTRY", "forbidden/Maxwell ancestry entered the graph", checks)
    checked(all(
        not row["symbolic_object_inventory"]["signed_tilt_coordinate_constructed"]
        and row["symbolic_object_inventory"]["banned_signed_tilt_scan_executed"]
        and not row["symbolic_object_inventory"]["banned_signed_tilt_pattern_detected"]
        and row["symbolic_object_inventory"]["collective_coordinates"] == ["X", "p"]
        for row in engines
    ), "ASSERT_NO_BANNED_TILT", "banned signed tilt convention constructed", checks)
    checked(all(
        not row["symbolic_object_inventory"]["two_body_objects_constructed"]
        and not row["symbolic_object_inventory"]["electric_sign_assertions"]
        and not row["symbolic_object_inventory"]["magnetism_sign_assertions"]
        and not row["symbolic_object_inventory"]["current_law_assertions"]
        for row in engines
    ), "ASSERT_ONE_BODY_SCOPE", "two-body sign/current object constructed", checks)
    checked(all(
        row["gates"]["G8"]["level1_source_partition_exact"]
        and row["gates"]["G8"]["level2_exactly_one"]
        and row["gates"]["G8"]["entry_count"] == 20
        and all(item["level2_disposition"] == "entry_witness" for item in row["gates"]["G8"]["records"])
        for row in engines
    ), "ASSERT_G8", "G8 two-level ratified disposition failed", checks)
    checked(all(
        row["gates"]["G11"]["orientation_fixture_ablation_zero"]
        and row["gates"]["G11"]["velocity_fixture_ablation_zero"]
        and row["gates"]["G11"]["static_orientation_fixture_survives"]
        and row["gates"]["G11"]["contamination"].startswith("UNRESOLVED")
        for row in engines
    ), "ASSERT_G11", "G11/converse control failed or contamination hidden", checks)
    checked(all(
        row["gates"]["G10"]["zero_would_be_recorded_honestly"]
        and row["gates"]["G10"]["no_zero_massaging_path"]
        for row in engines
    ), "ASSERT_G10", "G10 process record failed", checks)
    checked(all(
        set(row["gates"]) >= {"G2", "G5", "G6", "G8", "G10", "G11"}
        and row["gates"]["G5"]["ambient_quarantine_enforced"]
        and row["gates"]["G6"]["every_output_endpoint_keyed"]
        for row in engines
    ), "ASSERT_STANDING_GATES", "G2/G5/G6 record/quarantine failed", checks)
    checked(all(
        row["G4_known_nonzero_controls"]["exact_coverage"]
        and ids(row["G4_known_nonzero_controls"]["claimed_zero_ids"])
        == ids(row["G4_known_nonzero_controls"]["covered_zero_ids"])
        and all(control["fixture_nonzero"] for control in row["G4_known_nonzero_controls"]["controls"])
        for row in engines
    ), "ASSERT_G4_ZERO_CONTROLS", "claimed zero lacks a known-nonzero control", checks)
    checked(all(
        row["reconciliation"]["exact_set_equal"]
        and row["reconciliation"]["successor_count"] == 14125
        and row["reconciliation"]["G9_preserved_count"] == 90
        and not row["reconciliation"]["new_witness_minted"]
        and not row["reconciliation"]["upstream_record_modified"]
        for row in engines
    ), "ASSERT_RECONCILIATION", "immutable-overlay successor set is not exact", checks)
    checked(all(
        row["tilt"]["formalism"]["partition_successor"]["computed_owner"]
        == "UNRESOLVED(open_leaf:return_closure)"
        and not row["tilt"]["formalism"]["partition_successor"]["upstream_fact_reused_as_closed_owner"]
        for row in engines
    ), "ASSERT_PARTITION_SUCCESSOR", "partition owner was reused as a closed upstream fact", checks)
    checked(all(
        row["dimensional_firewall"]["all_inline_homogeneous"]
        and row["dimensional_firewall"]["cross_expression_consistent"]
        and row["dimensional_firewall"]["no_back_solved_carrier"]
        and not row["dimensional_firewall"]["back_solved_free_carriers"]
        and row["dimensional_firewall"]["firing_ablation"]["heterogeneity_detected"]
        and row["tilt"]["formalism"]["statics"]["dimensions_restored"]
        and row["tilt"]["formalism"]["susceptibility"]["dimensions_restored"]
        and row["coupling_package"]["7.5d"]["formal_total_response"]["dimensions_restored"]
        and all(item["dimensions_restored"] for item in row["coupling_package"]["7.5c"])
        and all(item["dimensions_restored"] for item in row["G4_known_nonzero_controls"]["controls"])
        for row in engines
    ), "ASSERT_DIMENSIONAL_FIREWALL", "inline dimensional firewall failed", checks)
    checked(all(
        all(cell["physics_status"]["enum"] == "TILT_UNRESOLVED" for cell in row["tilt"]["cells"])
        and all(cell["physics_status"]["enum"] == "UNRESOLVED"
                and cell["j_proportional_sV"] == "NOT_CLASSIFIABLE_UNRESOLVED"
                for cell in row["coupling_package"]["7.5d"]["cells"])
        for row in engines
    ), "ASSERT_UNRESOLVED_HONESTY", "OPEN-descended output promoted to an OK/classification", checks)
    checked(all(
        row["A9_external_verification"]["object_count"] == len(row["A9_external_verification"]["object_ids"])
        and len(row["A9_external_verification"]["object_ids"])
        == len(set(row["A9_external_verification"]["object_ids"]))
        and row["A9_external_verification"]["generated_category_union_exact"]
        and row["A9_external_verification"]["coverage_map_exact"]
        and row["A9_external_verification"]["computed_class_equivalence_all"]
        and ids(item["object_id"] for item in row["A9_external_verification"]["coverage_map"])
        == ids(row["A9_external_verification"]["object_ids"])
        and row["A9_external_verification"]["exactly_one_class_per_object"]
        and row["A9_external_verification"]["status"] == "AWAITING_ORCHESTRATOR_EXTERNAL_VERIFICATION"
        for row in engines
    ), "ASSERT_A9_COVERAGE", "A9 external coverage map is incomplete", checks)
    checked(
        sympy["independent_route"] != wolfram["independent_route"]
        and sympy["independent_route"] != "same-route"
        and wolfram["independent_route"] != "same-route",
        "ASSERT_ENGINE_INDEPENDENCE", "engine derivation routes are not independent", checks,
    )
    sym_view = semantic_view(sympy)
    wol_view = semantic_view(wolfram)
    checked(sym_view == wol_view, "ASSERT_DUAL_ENGINE_AGREEMENT",
            "canonical formulas, grids, dependencies, or verdicts disagree", checks)
    required_headlines = {
        "stage0_binding", "tilt_7.4", "coupling_7.5a", "coupling_7.5b",
        "coupling_7.5c", "coupling_7.5d", "mass_charge_split", "parity_census",
        "mouth_datum", "G2", "G4", "G5", "G6", "G8", "G10", "G11",
        "reconciliation", "dimensional_firewall", "native_ancestry",
        "A9_external_verification",
    }
    checked(all(set(row["headline_entries"]) == required_headlines for row in engines),
            "ASSERT_HEADLINES", "machine headline inventory incomplete", checks)
    return checks, sym_view


def build_successors(
    reconciliation: dict[str, Any],
    tilt_slots: list[str],
    engine_keys: list[str],
    mutation: str | None = None,
) -> dict[str, Any]:
    records = []
    for upstream in reconciliation["records"]:
        key = upstream["successor_key"]
        preserved = upstream["phase_C_stage0_routing"] == "PRESERVED_G9_EXACT_REFERENCE"
        if key.startswith("B1_LEAF|"):
            dependencies = ["tilt:" + key.split("|", 1)[1]]
        else:
            dependencies = tilt_slots
        records.append({
            "successor_key": key,
            "upstream_artifact_digest": upstream["upstream_artifact_digest"],
            "canonical_upstream_id_or_schema_path": upstream["canonical_upstream_id_or_schema_path"],
            "frozen_upstream_disposition_verbatim": upstream["frozen_upstream_disposition_verbatim"],
            "stage0_routing_verbatim": upstream["phase_C_stage0_routing"],
            "phase_C_production_disposition": (
                "PRESERVED_G9_EXACT_REFERENCE"
                if preserved else "TILT_UNRESOLVED(named ratified slots)"
            ),
            "dependency_slots": [] if preserved else dependencies,
            "witness_references": [] if preserved else [f"witness:{slot_id}" for slot_id in dependencies],
            "new_witness_minted": False,
            "upstream_record_modified": False,
            "upstream_derived_value_disagreement": False,
        })
    if mutation == "ASSERT_SUCCESSOR_ASSEMBLY":
        records.pop()
    keys = ids([row["successor_key"] for row in records])
    expected = ids(reconciliation["expected_ids"])
    if keys != expected or keys != ids(engine_keys):
        raise AssertionDeath("ASSERT_SUCCESSOR_ASSEMBLY", "successor assembly changed the exact inventory")
    omitted = set(expected) - set(keys)
    invented = set(keys) - set(expected)
    return {
        "schema_version": "U1_PHASE_C_PRODUCTION_SUCCESSOR_RECORDS_V1",
        "immutable_overlay": True,
        "record_count": len(records),
        "expected_exact_set_equal": keys == expected == ids(engine_keys),
        "duplicate_count": len(records) - len(keys),
        "omitted_count": len(omitted),
        "invented_count": len(invented),
        "G9_preserved_count": sum(
            row["phase_C_production_disposition"] == "PRESERVED_G9_EXACT_REFERENCE"
            for row in records
        ),
        "records": records,
    }


def summary_text(
    results: dict[str, Any], agreement: dict[str, Any], mutation: str | None = None
) -> str:
    gate = results["gates"]
    lines = [
        "# U1 Phase C production summary",
        "",
        f"Run classification: `{results['run_classification']}`. Computational integrity: `{results['computational_integrity']}`.",
        f"Stage-0 binding: `{results['stage0_contract_digest']}`; dual engines: `{agreement['status']}`.",
        "",
        "## Headline entries",
        "",
        "- `stage0_binding`: ratified digest re-asserted; live environment re-assertion is a runner preflight before each evaluation.",
        "- `tilt_7.4`: `TILT_UNRESOLVED` in all 240 cells. The total five-channel residual, endpoint existence branches, stiffness/anchoring dependence, pole equation, and four susceptibility branch classes are emitted; no steady value is fabricated.",
        "- `coupling_7.5a`: native `S_body[Omega_c]` is derived; `Sigma` and `partial_Omega_c` terms are `UNRESOLVED` on named ratified surface slots.",
        "- `coupling_7.5b`: all four `J_A` and all 16 ordered `deltaO_AB` entries, including 12 off-diagonals, are explicit functionals and `UNRESOLVED` on their ratified slots.",
        "- `coupling_7.5c`: E1 is variational/holonomic, E4 is multiplier virtual work, E5 is Rayleigh; endpoint data remain named-slot `UNRESOLVED` where applicable.",
        "- `coupling_7.5d`: all 960 total-response cells include sources plus the full mixed kernel and far-multipole operator; every cell is `UNRESOLVED`, so `j proportional sV` is not classifiable.",
        "- `mass_charge_split`: drain/mass and orientation/charge response functionals are separate; no unresolved contamination is silently set to zero.",
        "- `parity_census`: 40 one-body post-mouth records carry the required three-way authority tags; parity components are honestly profile-blocked.",
        "- `mouth_datum`: five descriptive endpoint records state the held field datum and defer fixed-source versus fixed-displacement/defect character to U2 where required.",
        f"- `G2`: `{gate['G2']['status']}`.",
        "- `G4`: every claimed structural/control zero is exactly covered by a known-nonzero same-pipeline control.",
        f"- `G5`: `{gate['G5']['status']}`; ambient Green/operator quarantine is active.",
        f"- `G6`: `{gate['G6']['status']}` over E1-E5; physics sensitivity is unresolved, not collapsed.",
        f"- `G8`: `{gate['G8']['status']}`; all 20 level-1 entries retain exactly one ratified level-2 `entry_witness` disposition.",
        f"- `G10`: `{gate['G10']['status']}`; any future derived tilt zero has an honest result path.",
        f"- `G11`: `{gate['G11']['status']}`; mass-only and zero-velocity converse controls execute, while physical contamination stays unresolved.",
        f"- `reconciliation`: {results['reconciliation']['successor_count']} immutable-overlay successors, including {results['reconciliation']['G9_preserved_count']} exact G9 references; no upstream record or witness is rewritten.",
        "- `dimensional_firewall`: inline restored-unit classes are homogeneous in both engines, cross-consistent, with no back-solved carrier and a firing dimension defect.",
        "- `native_ancestry`: runtime typed-root graph is green; no Maxwell/Larmor/Coulomb, point-current, or declared `j=sV` ancestry exists.",
        "- `A9_external_verification`: exact coverage inventory emitted; all four external legs remain orchestrator-owned and pending before banking.",
        "",
        "## Scope and conditionality",
        "",
        "Every cell is keyed by E1-E5, both ambient branches, three closure branches, and each of the eight `GAP_OPEN_FIELD_PROFILE` strata. This is a one-body collective-coordinate result: no two-body force/sign, magnetic sign, electric sign, or current law is asserted.",
        "",
        "`UNRESOLVED` is the production physics result wherever a predicate depends on the ratified OPEN leaves. The formal `c_sv/Delta` projection is retained, but it is not evaluated on a nonexistent total response.",
        "",
        "## Armor and external status",
        "",
        "The dual-engine comparator is green. Mutation, process-tree containment, evaluated-code closure, and final runner-terminal evidence are separate machine artifacts produced by the sealed runner. This summary does not claim the four external A9 legs have run.",
    ]
    if mutation == "ASSERT_SUMMARY_COMPLETE":
        lines = [line for line in lines if "`A9_external_verification`" not in line]
    text = "\n".join(lines) + "\n"
    for headline in results["headline_entries"]:
        if f"`{headline}`" not in text:
            raise AssertionDeath("ASSERT_SUMMARY_COMPLETE", f"summary omitted {headline}")
    return text


def compare_and_write(args: argparse.Namespace) -> int:
    sympy = load_yaml(Path(args.sympy).resolve())
    wolfram = load_yaml(Path(args.wolfram).resolve())
    bundle_dir = Path(args.bundle_dir).resolve()
    projection = load_yaml(bundle_dir / "projection_freeze.yaml")
    if args.mutation:
        if args.mutation_side in {"sympy", "both"}:
            apply_mutation(sympy, args.mutation, "sympy")
        if args.mutation_side in {"wolfram", "both"}:
            apply_mutation(wolfram, args.mutation, "wolfram")
    checks, semantic = assert_engine(sympy, wolfram, projection)
    reached_asserts = {row["assert_id"] for row in checks}
    agreement = {
        "schema_version": "U1_PHASE_C_PRODUCTION_ENGINE_AGREEMENT_V1",
        "status": "ENGINE_AGREE",
        "stage0_contract_digest": RATIFIED_DIGEST,
        "sympy_route": sympy["independent_route"],
        "wolfram_route": wolfram["independent_route"],
        "canonical_semantic_digest": digest(semantic),
        "checks": checks,
        "check_count": len(checks),
        "tilt_cell_count": len(sympy["tilt"]["cells"]),
        "coupling_cell_count": len(sympy["coupling_package"]["7.5d"]["cells"]),
        "successor_count": sympy["reconciliation"]["successor_count"],
        "transliteration_fidelity_status": "AWAITING_FRESH_EXTERNAL_AUDIT",
    }
    if args.output:
        dump_yaml(Path(args.output).resolve(), agreement)
    if args.results_output or args.successors_output or args.summary_output:
        reconciliation = load_yaml(bundle_dir / "reconciliation_inventory.yaml")
        results = copy.deepcopy(sympy)
        results["schema_version"] = RESULT_SCHEMA
        results["engine"] = "DUAL_ENGINE_COMPARATOR_ASSERTED"
        results["run_classification"] = (
            "UNSEALED_SELF_TEST_NO_SEALED_RESULT" if args.self_test
            else "PRODUCTION_CANDIDATE_PENDING_RUNNER_ARMOR_AND_EXTERNAL_A9"
        )
        results["stage0_contract_digest"] = RATIFIED_DIGEST
        results["computational_integrity"] = "COMPUTATION_VALID"
        results["dual_engine_agreement"] = {
            "status": "ENGINE_AGREE",
            "canonical_semantic_digest": agreement["canonical_semantic_digest"],
            "check_count": agreement["check_count"],
        }
        results.pop("semantic_sink_digest", None)
        tilt_slots = ids(
            slot for slot in results["availability_contract"]["unresolved_slot_ids"]
            if slot.startswith("tilt:")
        )
        successors = build_successors(
            reconciliation, tilt_slots, results["reconciliation"]["expected_successor_keys"],
            args.mutation,
        )
        reached_asserts.add("ASSERT_SUCCESSOR_ASSEMBLY")
        if args.results_output:
            dump_yaml(Path(args.results_output).resolve(), results)
        if args.successors_output:
            dump_yaml(Path(args.successors_output).resolve(), successors)
        if args.summary_output:
            summary_path = Path(args.summary_output).resolve()
            summary_path.parent.mkdir(parents=True, exist_ok=True)
            summary_path.write_text(
                summary_text(results, agreement, args.mutation), encoding="utf-8"
            )
            reached_asserts.add("ASSERT_SUMMARY_COMPLETE")
    if args.control_assert and args.control_assert not in reached_asserts:
        raise AssertionDeath(
            args.control_assert, "defect-absent control did not reach its target require"
        )
    print(
        f"PHASEC_PRODUCTION_COMPARE_PASS checks={len(checks)} "
        f"tilt=240 coupling=960 successors=14125",
        flush=True,
    )
    return 0


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser()
    result.add_argument("--sympy")
    result.add_argument("--wolfram")
    result.add_argument("--bundle-dir")
    result.add_argument("--output")
    result.add_argument("--results-output")
    result.add_argument("--successors-output")
    result.add_argument("--summary-output")
    result.add_argument("--self-test", action="store_true")
    result.add_argument("--list-asserts", action="store_true")
    result.add_argument("--mutation")
    result.add_argument("--mutation-side", choices=("sympy", "wolfram", "both"), default="both")
    result.add_argument("--control-assert")
    return result


def main() -> int:
    args = parser().parse_args()
    if args.list_asserts:
        # Exercise the baseline to populate the true runtime check order; the
        # caller supplies artifacts so no static assertion list can drift.
        for required in ("sympy", "wolfram", "bundle_dir"):
            if not getattr(args, required):
                raise SystemExit(f"--list-asserts requires --{required.replace('_', '-')}")
        try:
            sympy = load_yaml(Path(args.sympy).resolve())
            wolfram = load_yaml(Path(args.wolfram).resolve())
            projection = load_yaml(Path(args.bundle_dir).resolve() / "projection_freeze.yaml")
            checks, _ = assert_engine(sympy, wolfram, projection)
            print(yaml.safe_dump({"assert_ids": [row["assert_id"] for row in checks]}, sort_keys=False))
            return 0
        except AssertionDeath as death:
            print(f"ASSERTION_FAILED {death.assert_id}: {death.detail}", file=sys.stderr)
            return 1
    for required in ("sympy", "wolfram", "bundle_dir"):
        if not getattr(args, required):
            raise SystemExit(f"missing --{required.replace('_', '-')}")
    try:
        return compare_and_write(args)
    except AssertionDeath as death:
        print(f"ASSERTION_FAILED {death.assert_id}: {death.detail}", file=sys.stderr)
        return 1
    except (KeyError, TypeError, ValueError) as failure:
        print(f"PHASEC_PRODUCTION_COMPARATOR_HARNESS_FAILED {failure}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
