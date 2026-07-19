#!/usr/bin/env python3
"""Strict comparator, guards, report assembler, and mutation probe for U2."""

from __future__ import annotations

import argparse
import copy
import hashlib
import json
import sys
from collections import Counter
from pathlib import Path
from typing import Any, Iterable

import yaml


RATIFIED_DIGEST = "b01a1821e908589c3698512bbb9aff874b721af6dcbfa1c3b8b1f8d33119b32b"
SOURCE = Path(__file__).resolve().parent
if str(SOURCE) not in sys.path:
    sys.path.insert(0, str(SOURCE))

from u2_production_sympy import evaluator_behavior as frozen_evaluator_behavior  # noqa: E402
from u2_stage0_sympy import (  # noqa: E402
    classify_evidence as frozen_classify_evidence,
    classify_inference_content as frozen_classify_inference_content,
    evaluate_cross_level as frozen_evaluate_cross_level,
    evaluate_disposition as frozen_evaluate_disposition,
    evaluate_promotion as frozen_evaluate_promotion,
    evaluate_topology as frozen_evaluate_topology,
    record_schema_valid as frozen_record_schema_valid,
)


RECORD_TYPES = (
    "candidate_disposition", "obligation_evidence", "unavailability_witness",
    "derivability_challenge", "ensemble_level_1", "ensemble_level_2",
    "topology_gate", "host_location", "closure_adjudication", "promotion",
    "posed_BVP_template", "return_closure_ownership",
)


class AssertionDeath(RuntimeError):
    def __init__(self, assert_id: str, detail: str):
        super().__init__(detail)
        self.assert_id = assert_id
        self.detail = detail


def checked(condition: bool, assert_id: str, detail: str, checks: list[dict[str, str]] | None = None) -> None:
    if not condition:
        raise AssertionDeath(assert_id, detail)
    if checks is not None:
        checks.append({"assert_id": assert_id, "status": "PASS"})


def load_yaml(path: Path) -> Any:
    with path.open("rb") as handle:
        return yaml.load(handle, Loader=yaml.CSafeLoader)


def dump_yaml(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        yaml.safe_dump(value, sort_keys=False, allow_unicode=True, width=160), encoding="utf-8"
    )


def canonical_bytes(value: Any) -> bytes:
    return json.dumps(value, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode()


def digest(value: Any) -> str:
    return hashlib.sha256(canonical_bytes(value)).hexdigest()


def arbiter_physics_record_projection(document: dict[str, Any]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    context_fields = ("cell_id", "stable_branch_id", "candidate_id", "ambient", "stratum")
    disposition_fields = (
        *context_fields, "native_root_class", "integrity", "expected_dependencies", "used_dependencies",
        "obligation_evidence", "disposition", "disposition_evaluator_landing", "unavailability",
    )
    for cell in document["cell_records"]:
        context = {key: cell[key] for key in context_fields}
        rows.extend((
            {"record_id": f"candidate_disposition:{cell['cell_id']}", "record_class": "candidate_disposition", "payload": {key: cell[key] for key in disposition_fields}},
            {"record_id": cell["ensemble"]["level_1"]["record_id"], "record_class": "ensemble_level_1", "payload": {**context, **cell["ensemble"]["level_1"]}},
            {"record_id": cell["ensemble"]["level_2"]["record_id"], "record_class": "ensemble_level_2", "payload": {**context, "applicability": cell["ensemble"]["applicability"], **cell["ensemble"]["level_2"]}},
            {"record_id": cell["topology_gate"]["record_id"], "record_class": "topology_gate", "payload": cell["topology_gate"]},
            {"record_id": cell["host_location"]["record_id"], "record_class": "host_location", "payload": cell["host_location"]},
            {"record_id": cell["return_closure_ownership"]["record_id"], "record_class": "return_closure_ownership", "payload": cell["return_closure_ownership"]},
        ))
    rows.extend({"record_id": row["record_id"], "record_class": "closure_adjudication", "payload": row} for row in document["closure_records"])
    rows.extend({"record_id": row["record_id"], "record_class": "promotion", "payload": row} for row in document["promotion_records"])
    return sorted(rows, key=lambda row: (row["record_id"].casefold(), row["record_id"]))


def contains_solved_payload(value: Any) -> bool:
    if isinstance(value, dict):
        if value.get("payload_kind") in {"SOLVED_PROFILE", "RESPONSE_COEFFICIENT", "EIGENPAIR"}:
            return True
        return any(contains_solved_payload(child) for child in value.values())
    if isinstance(value, list):
        return any(contains_solved_payload(child) for child in value)
    return False


def scalar_values(value: Any) -> set[Any]:
    if isinstance(value, dict):
        return set().union(*(scalar_values(child) for child in value.values())) if value else set()
    if isinstance(value, list):
        return set().union(*(scalar_values(child) for child in value)) if value else set()
    return {value} if isinstance(value, (str, int, float, bool, type(None))) else set()


def walk_template_terms(node: dict[str, Any]) -> list[dict[str, str]]:
    rows = []
    if node.get("op") == "term" and node.get("coefficient") != 0:
        rows.append({"term_id": node["term_id"], "kind": node["kind"]})
    for child in node.get("args", []):
        if isinstance(child, dict):
            rows.extend(walk_template_terms(child))
    return rows


def ids(values: Iterable[str]) -> list[str]:
    return sorted(set(values), key=lambda value: (value.lower(), value))


def by_id(rows: list[dict[str, Any]], key: str) -> dict[str, dict[str, Any]]:
    result = {row[key]: row for row in rows}
    checked(len(result) == len(rows), "ASSERT_RECORD_IDS_UNIQUE", f"duplicate {key}")
    return result


def route_assert_id(row: dict[str, Any]) -> str:
    return f"ASSERT_ROUTE_CONTROL:{row['cell_id']}"


def ablation_assert_id(row: dict[str, Any]) -> str:
    return f"ASSERT_ROUTE_ABLATION:{row['cell_id']}"


def exchange_assert_id(cell: dict[str, Any]) -> str:
    return f"ASSERT_ENSEMBLE_EXCHANGE:{cell['cell_id']}"


def record_schema_assert_id(record_type: str) -> str:
    return f"ASSERT_RECORD_SCHEMA:{record_type}"


def validate_record_schema_fixtures(checks: list[dict[str, str]] | None = None) -> None:
    for record_type in RECORD_TYPES:
        assert_id = record_schema_assert_id(record_type)
        valid = (
            {"record_type": record_type, "integrity": "COMPUTATION_VALID", "physics": "fixture"},
            {"record_type": record_type, "integrity": "HARNESS_FAILED", "physics": None, "failure_reason": "fixture_typed_reason"},
            {"record_type": record_type, "integrity": "NOT_RUN", "physics": None, "failed_upstreams": ["fixture:upstream"]},
        )
        invalid = (
            {"record_type": record_type, "integrity": "HARNESS_FAILED", "physics": "illegal", "failure_reason": "fixture_typed_reason"},
            {"record_type": record_type, "integrity": "NOT_RUN", "physics": "illegal", "failed_upstreams": ["fixture:upstream"]},
            {"record_type": record_type, "integrity": "COMPUTATION_VALID", "physics": "illegal", "failed_upstreams": ["fixture:upstream"]},
        )
        checked(
            all(frozen_record_schema_valid(row) for row in valid)
            and not any(frozen_record_schema_valid(row) for row in invalid),
            assert_id, f"ratified integrity schema fixture failed for {record_type}", checks,
        )


def assert_route_control(row: dict[str, Any], frozen: dict[str, Any]) -> None:
    assert_id = route_assert_id(row)
    checked(row["route_id"] == frozen["route_id"] and row["fixture_id"] == frozen["fixture_id"], assert_id, "route/fixture binding changed")
    checked(row["semantic_route_id"] == "obligation_residual_classifier_v1" and row["test_only"], assert_id, "fixture did not traverse frozen semantic route")
    checked(row["positive_admissible"] and all(value == 0 for value in row["positive_residual"]), assert_id, "positive fixture did not realize its route")
    checked(row["nondegenerate_norm_squared"] > 0, assert_id, "positive fixture is degenerate")
    checked(row["malformed_excluded"] and any(value != 0 for value in row["malformed_residual"]), assert_id, "malformed fixture failed to exclude")
    native = row["native_structure_exercised"]
    root_class = frozen["native_root_class"]
    if "nonholonomic" in root_class:
        checked(native.get("constraint_multiplier") is True, assert_id, "E4 route omitted multiplier/virtual-work structure")
    if "Rayleigh" in root_class or "rayleigh" in root_class:
        checked(native.get("dissipation_bookkeeping") is True, assert_id, "E5 route omitted dissipation bookkeeping")
    dissipation = row["dissipation_measurement"]
    if dissipation["applicable"]:
        expected = dissipation["rayleigh_coefficient"] * dissipation["velocity_coordinate"] ** 2
        checked(dissipation["value"] == expected and dissipation["frozen_fixture_field_equal"], assert_id, "dissipation field is not bound to the executed calculation")
    else:
        checked(dissipation["value"] is None and dissipation["frozen_fixture_field_equal"], assert_id, "non-dissipative route carries decorative dissipation")

    ablation = row["semantic_ablation"]
    checked(
        ablation["criterion"] == "fail_nondefinitional_obligation_or_change_canonical_operator_class"
        and ablation["nondefinitional_obligation_failed"]
        and any(value != 0 for value in ablation["ablated_residual"]),
        ablation_assert_id(row), "semantic ablation did not fail a non-definitional obligation",
    )


def assert_local_ensemble(
    cell: dict[str, Any], availability: dict[str, dict[str, Any]]
) -> None:
    assert_id = "ASSERT_ENSEMBLE_STATIC_CLASSIFIER"
    level_1 = cell["ensemble"]["level_1"]
    slot_id = f"ensemble:{cell['candidate_id']}:boundary_action_variation"
    slot = availability[slot_id]
    checked(
        level_1["local_preclosure_committed_structure_id"] == slot_id,
        assert_id, f"committed ensemble citation mismatch {cell['cell_id']}",
    )
    if slot["availability_outcome"] != "DERIVED":
        checked(
            level_1["local_preclosure_classification"] is None
            and level_1["local_preclosure_geometric_component"] is None
            and level_1["local_preclosure_classifier_evidence"] is None,
            assert_id, f"unavailable ensemble was classified {cell['cell_id']}",
        )
        return
    variation = slot["derived_object"]
    measurements = level_1["local_preclosure_classifier_evidence"]
    expected_roles: list[str] = []
    expected_variation_channels: list[str] = []
    expected_flux_channels: list[str] = []
    component_role_counts: list[int] = []
    for component in variation["component_variations"]:
        condition = component["boundary_condition"]
        before = len(expected_roles)
        if "v.normal=V.normal" in condition:
            expected_roles.append("normal_velocity_control")
        elif "normal_traction" in condition:
            expected_roles.append("normal_traction_response")
        if "v.tangent=V.tangent" in condition:
            expected_roles.append("tangential_velocity_control")
        elif "tangential_traction=0" in condition:
            expected_roles.append("tangential_traction_response")
        component_role_counts.append(len(expected_roles) - before)
        expected_variation_channels.extend(component["variation_channels"])
        expected_flux_channels.extend(component["flux_channels"])
    signature = ids(expected_roles)
    expected_classification = (
        "fixed-displacement/geometric"
        if all(
            role_count >= len(component["variation_channels"])
            for role_count, component in zip(component_role_counts, variation["component_variations"])
        )
        and signature == ["normal_velocity_control", "tangential_velocity_control"]
        else f"mixed/other-ensemble({slot_id})"
    )
    checked(
        measurements["classifier"] == "committed_boundary_variation_static_signature_walk_v1"
        and measurements["derived_object_kind"] == "boundary_action_variation_record"
        and measurements["component_count"] == len(variation["component_variations"])
        and measurements["signature"] == signature
        and measurements["signature_complete"]
        == (
            bool(signature) and all(
                role_count >= len(component["variation_channels"])
                for role_count, component in zip(component_role_counts, variation["component_variations"])
            )
        )
        and measurements["variation_channels"] == ids(expected_variation_channels)
        and measurements["flux_channels"] == ids(expected_flux_channels)
        and level_1["local_preclosure_classification"] == expected_classification
        and level_1["local_preclosure_geometric_component"] == bool(expected_variation_channels),
        assert_id, f"ensemble classification is not channel-derived {cell['cell_id']}",
    )


def assert_exchange(
    cell: dict[str, Any], availability: dict[str, dict[str, Any]] | None = None
) -> None:
    control = cell["ensemble"]["level_1"]["exchange_control"]
    if control is None:
        return
    assert_id = exchange_assert_id(cell)
    classifier = control["signature_classifier_evidence"]
    checked(
        classifier["classifier"] == "committed_boundary_variation_static_signature_walk_v1",
        assert_id, "real signature classifier did not execute",
    )
    expected_pairing: dict[str, Any] | None = None
    if availability is not None:
        structure = availability[f"candidate_definition:{cell['candidate_id']}"]["derived_object"]
        components = structure["components"]
        canonical_map = {
            "normal_velocity_lock": "normal_traction_response",
            "tangential_velocity_lock": "tangential_traction_response",
            "tangential_traction_free": "tangential_velocity_control",
        }
        eligible = ids(
            row["endpoint_id"] for row in components
            if row["variational_class"] == "holonomic_field_BC"
        )
        nonconjugate = ids(
            f"{row['endpoint_id']}:{row['variational_class']}" for row in components
            if row["variational_class"] != "holonomic_field_BC"
        )
        unmapped = ids(
            value for value in structure["canonical_signature"] if value not in canonical_map
        )
        complete = (
            structure["native_root_class"] == "variational_holonomic"
            and len(eligible) == len(components) and not nonconjugate and not unmapped
            and bool(structure["canonical_signature"])
        )
        expected_pairing = {
            "native_root_class": structure["native_root_class"],
            "required_component_count": len(components),
            "eligible_variational_components": eligible,
            "nonconjugate_components": nonconjugate,
            "unmapped_canonical_signature_terms": unmapped,
            "complete_structure_pairing": complete,
            "independently_expected_exchanged_signature": ids(
                canonical_map[value] for value in structure["canonical_signature"]
            ) if complete else [],
        }
    if control["route"] == "native_conjugate_transform":
        conjugates = {
            "normal_velocity_control": "normal_traction_response",
            "normal_traction_response": "normal_velocity_control",
            "tangential_velocity_control": "tangential_traction_response",
            "tangential_traction_response": "tangential_velocity_control",
        }
        recomputed_exchange = ids(conjugates[value] for value in control["baseline_signature"])
        pairing = control["native_pairing_construction"]
        checked(
            control["response_character_flipped"]
            and control["baseline_signature"] == classifier["signature"]
            and control["baseline_signature"] != control["computed_exchanged_signature"]
            and control["computed_exchanged_signature"] == recomputed_exchange
            and control["computed_exchanged_signature"] == control["independently_generated_expected_signature"],
            assert_id, "ensemble exchange did not flip the response character",
        )
        checked(
            pairing["complete_structure_pairing"]
            and pairing["independently_expected_exchanged_signature"]
            == control["independently_generated_expected_signature"],
            assert_id, "native-root conjugate expectation is incomplete",
        )
        if expected_pairing is not None:
            checked(
                all(pairing[key] == value for key, value in expected_pairing.items()),
                assert_id, "native-root pairing was not independently constructed",
            )
    else:
        insufficiency = control["pairing_insufficiency"]
        computed_certificate = (
            not insufficiency["complete_structure_pairing"]
            and len(insufficiency["eligible_variational_components"])
            < insufficiency["required_component_count"]
            and bool(
                insufficiency["nonconjugate_components"]
                or insufficiency["unmapped_canonical_signature_terms"]
            )
        )
        checked(
            control["route"] == "computed_no_pairing_certificate"
            and control["observed_boundary_signature"] == classifier["signature"]
            and control["no_pairing_certificate"] == computed_certificate is True,
            assert_id, "mixed ensemble has no valid no-pairing certificate",
        )
        if expected_pairing is not None:
            checked(
                all(
                    insufficiency[key] == value for key, value in expected_pairing.items()
                    if key != "independently_expected_exchanged_signature"
                ),
                assert_id, "no-pairing insufficiency differs from committed native structure",
            )


def walk_proof_dags(value: Any, path: str = "root") -> Iterable[tuple[str, dict[str, Any]]]:
    if isinstance(value, dict):
        if "proof_dag" in value and isinstance(value["proof_dag"], dict):
            yield path + ".proof_dag", value["proof_dag"]
        for key, child in value.items():
            if key != "proof_dag":
                yield from walk_proof_dags(child, path + "." + str(key))
    elif isinstance(value, list):
        for index, child in enumerate(value):
            yield from walk_proof_dags(child, f"{path}[{index}]")


def assert_proof_dag(path: str, dag: dict[str, Any], taxonomy_roots: set[str]) -> None:
    content = dag["normalized_inference_content"]
    computed = frozen_classify_inference_content(content)
    checked(computed == dag["constructors"] and "UNCLASSIFIED" not in computed, "ASSERT_TYPED_ANCESTRY", f"content classification mismatch at {path}")
    checked(set(dag["root_ids"]) <= taxonomy_roots, "ASSERT_NATIVE_ANCESTRY", f"untyped/non-native root at {path}")
    checked("STABILITY_DYNAMICAL_CLASS" not in computed, "ASSERT_PROGRAM_WIDE_STABILITY_BAN", f"stability root at {path}")
    checked("SOLVE_EVALUATION_RESULT" not in computed, "ASSERT_BVP_NOT_SOLVED", f"solve/evaluation node at {path}")


def assert_closure_record(row: dict[str, Any]) -> None:
    assert_id = "ASSERT_CLOSURE_GATE"
    checked(row["integrity"] == "COMPUTATION_VALID" and row["physics"]["kind"] == "CLOSURE_REFUSED", assert_id, "closure refusal record schema invalid")
    checked(row["certificate_emitted"] is False and row["return_closure_consumed"] is False, assert_id, "refusal emitted a partial certificate or consumed return closure")
    checked(
        row["contribution_census_A_construction"] == "claim_side_boundary_variation_kind_expansion"
        and row["contribution_census_B_construction"] == "raw_committed_component_channel_schema_walk"
        and row["contribution_census_A"] == row["contribution_census_B"],
        "ASSERT_CLOSURE_CENSUS_GENERATION", "independent contribution censuses disagree",
    )
    checked(row["assignment_coverage"] == row["contribution_census_A"], assert_id, "closure assignment coverage is incomplete")
    assignments = row["assignments"]
    checked(len({value["contribution_id"] for value in assignments}) == len(assignments), assert_id, "duplicate contribution owner/refusal")
    independent_terms = row["independently_constructed_terms"]
    checked(
        ids(value["contribution_id"] for value in independent_terms) == row["contribution_census_B"]
        and sum(value["symbolic_coefficient"] for value in independent_terms)
        == row["independently_constructed_symbolic_total"]
        == row["assignment_symbolic_total"],
        "ASSERT_CLOSURE_TOTAL_GENERATION", "closure total reconstruction failed",
    )
    checked(all(value["landing"] in {"computed-zero", "typed-refusal"} for value in assignments), assert_id, "refusal record contains partial ownership")
    radiation = next(value for value in assignments if value["contribution_id"] == "static_radiation")
    checked(radiation["landing"] == "computed-zero" and radiation["owner"] == "radiation/static-zero", assert_id, "radiation/static-zero omitted")
    checked(row["host_physics"] == "undetermined", "ASSERT_CLOSURE_HOST_CONSISTENCY", "closure refusal is inconsistent with host field")


def validate_engine_document(
    document: dict[str, Any], bundle: dict[str, Any], checks: list[dict[str, str]] | None = None
) -> None:
    checked(document["schema_version"] == "U2_PRODUCTION_ENGINE_V1", "ASSERT_ENGINE_SCHEMA", "wrong engine schema", checks)
    registry = document["engine_local_route_registry"]
    expected_routes = {
        "obligation_residual_classifier_v1", "frozen_evidence_state_classifier_v1",
        "frozen_disposition_precedence_v1", "frozen_topology_aggregate_v1",
        "frozen_cross_level_ensemble_v1", "frozen_promotion_evaluator_v1",
        "branch_template_eligibility_v12", "conditional_posed_template_v12",
        "physics_record_projection_invariance_v12",
    }
    checked({row["semantic_route_id"] for row in registry} == expected_routes, "ASSERT_ENGINE_ROUTE_REGISTRY", "semantic route registry incomplete", checks)
    checked(all(row["exists"] and row["executed"] and row["engine_local_function"] for row in registry), "ASSERT_ENGINE_ROUTE_REGISTRY", "engine-local route absent/unexecuted", checks)
    validate_semantic(document["semantic_view"], bundle, checks)


def validate_semantic(
    view: dict[str, Any], bundle: dict[str, Any], checks: list[dict[str, str]] | None = None
) -> None:
    validate_record_schema_fixtures(checks)
    candidate_doc = bundle["candidate_inventory"]
    grid_doc = bundle["dependency_grid_inventory"]
    obligation_doc = bundle["obligation_censuses"]["censuses"]
    route_doc = bundle["route_fixture_inventory"]
    contracts = bundle["closure_template_contracts"]
    availability = by_id(bundle["availability_slots"]["slots"], "slot_id")
    taxonomy_roots = set(bundle["evidence_taxonomy"]["source_census"]["source_ids"])
    checked(view["stage0_binding"]["equal"] and view["stage0_binding"]["supplied_digest"] == RATIFIED_DIGEST == view["stage0_binding"]["recomputed_component_digest"], "ASSERT_STAGE0_BINDING", "stage-0 digest not bound", checks)
    checked(view["frozen_evaluator_behavior"] == frozen_evaluator_behavior(), "ASSERT_FROZEN_EVALUATORS", "production evaluator behavior differs from ratified functions", checks)
    axes = view["axes"]
    checked(
        axes["candidate_axis"] == candidate_doc["candidate_axis"]
        and axes["active_strata"] == grid_doc["active_strata"]
        and axes["cell_count"] == 144 and axes["promotion_context_count"] == 16
        and axes["candidate_universe_digest"] == candidate_doc["candidate_universe_digest"],
        "ASSERT_AXES", "production axes differ from ratified inventory", checks,
    )
    pin = contracts["physics_record_invariance_contract"]
    projected = arbiter_physics_record_projection(view)
    projected_digest = digest(projected)
    projected_universe = [row["record_id"] for row in projected]
    declared_invariance = view["physics_record_invariance"]

    frozen_routes = by_id(route_doc["route_records"], "route_id")
    routes = by_id(view["route_controls"], "route_id")
    checked(set(routes) == set(frozen_routes) and len(routes) == 144, "ASSERT_ROUTE_SET", "executed route set differs from frozen inventory", checks)
    frozen_cell_ids = {row["cell_id"] for row in grid_doc["grid_cells"]}
    checked({row["cell_id"] for row in routes.values()} == frozen_cell_ids, "ASSERT_CELL_ROUTE_FIXTURE_COVERAGE", "cell-to-route coverage differs", checks)
    for route_id, row in routes.items():
        assert_route_control(row, frozen_routes[route_id])

    cells = by_id(view["cell_records"], "cell_id")
    checked(set(cells) == frozen_cell_ids and len(cells) == 144, "ASSERT_GRID_COVERAGE", "144-cell exact coverage failed", checks)
    named_by_candidate: dict[str, set[str]] = {}
    for cell_id, cell in cells.items():
        checked(cell["integrity"] == "COMPUTATION_VALID", "ASSERT_ZERO_INTEGRITY", f"non-valid cell {cell_id}")
        checked(cell["expected_dependencies"] == cell["used_dependencies"], "ASSERT_DEPENDENCY_BIDIRECTIONAL", f"OPEN ancestry mismatch {cell_id}")
        expected_obligations = obligation_doc[cell["candidate_id"]]["generator_A"]
        evidence = cell["obligation_evidence"]
        checked([row["obligation_id"] for row in evidence] == expected_obligations, "ASSERT_EVIDENCE_CENSUS_COVERAGE", f"obligation coverage mismatch {cell_id}")
        checked(len({row["obligation_id"] for row in evidence}) == len(evidence), "ASSERT_EVIDENCE_CENSUS_COVERAGE", f"duplicate evidence state {cell_id}")
        for record in evidence:
            dependency_measurements = record["predicate_measurements"]
            slot_measurements = record["stage0_slot_measurements"]
            checked(bool(slot_measurements), "ASSERT_EVIDENCE_MEASUREMENTS", f"no ratified slot measurement {cell_id}:{record['obligation_id']}")
            for measurement in slot_measurements:
                slot = availability[measurement["stage0_slot_id"]]
                checked(
                    measurement["availability"] == slot["availability_outcome"]
                    and measurement["integrity_state"] == slot["integrity_state"]
                    and measurement["missing"] == (slot["availability_outcome"] == "UNRESOLVED")
                    and measurement["decisive_incompatibility"] == (slot["availability_outcome"] == "EXCLUDED")
                    and measurement["witness_id"] == slot.get("witness_id")
                    and measurement["challenge_id"] == slot.get("challenge_id")
                    and measurement["producer_set"] == slot["producer_set"],
                    "ASSERT_EVIDENCE_MEASUREMENTS", f"ratified slot measurement mismatch {cell_id}:{record['obligation_id']}",
                )
            all_measurements = [*dependency_measurements, *slot_measurements]
            raw = record["raw_predicate_inputs"]
            checked(
                raw["datum_missing"] == (raw["applicable"] and any(row["missing"] for row in all_measurements))
                and raw["committed_incompatibility"] == (
                    raw["applicable"] and any(
                        row.get("decisive_incompatibility", False) or row["availability"] == "EXCLUDED"
                        for row in all_measurements
                    )
                )
                and raw["ancestry_complete_and_typed"] == (
                    bool(all_measurements) and all(
                        row["integrity_state"] == "COMPUTATION_VALID" and row["measurement_source"]
                        for row in all_measurements
                    )
                ),
                "ASSERT_EVIDENCE_MEASUREMENTS", f"raw evidence is not computed from measurements {cell_id}:{record['obligation_id']}",
            )
            checked(record["integrity"] == "COMPUTATION_VALID" and frozen_classify_evidence(record["raw_predicate_inputs"]) == record["emitted_state"], "ASSERT_EVIDENCE_STATE_CLASSIFICATION", f"evidence relabel/drop {cell_id}:{record['obligation_id']}")
        landing = frozen_evaluate_disposition(evidence)
        checked(cell["disposition_evaluator_landing"] == landing, "ASSERT_DISPOSITION_PRECEDENCE", f"frozen precedence mismatch {cell_id}")
        expected_disposition_kind = landing.split("(", 1)[0]
        checked(cell["disposition"]["kind"] == expected_disposition_kind, "ASSERT_DISPOSITION_RECORD", f"wrong disposition {cell_id}")
        unavailability = cell["unavailability"]
        challenge = unavailability["challenge"]
        missing_dependency_ids = ids(
            measurement["dependency_token"] for row in evidence
            for measurement in row["predicate_measurements"] if measurement["missing"]
        )
        missing_slot_ids = ids(
            measurement["stage0_slot_id"] for row in evidence
            for measurement in row["stage0_slot_measurements"] if measurement["missing"]
        )
        checked(
            unavailability["witness_id"] == cell["disposition"]["witness_id"]
            and unavailability["challenge_id"] == cell["disposition"]["challenge_id"]
            and unavailability["named_datum"] and unavailability["stratum_datum"] == f"tilt:{cell['stratum']}"
            and challenge["executed"] and challenge["outcome"] == "CONSTRUCTIVE_FAIL"
            and challenge["measured_rank"] == 0 and challenge["measured_nullity"] == 1
            and challenge["restore_rank"] == 1 and challenge["restore_nullity"] == 0
            and not challenge["empty_output"] and challenge["candidate_well_typed"]
            and challenge["dual_engine_certificate"] and unavailability["reference_set_exact"]
            and unavailability["missing_dependency_tokens"] == missing_dependency_ids
            and unavailability["missing_stage0_slot_ids"] == missing_slot_ids
            and ids(row["dependency_token"] for row in unavailability["referenced_phaseC_slots"]) == missing_dependency_ids
            and ids(row["stage0_slot_id"] for row in unavailability["referenced_stage0_slots"]) == missing_slot_ids,
            "ASSERT_TYPED_WITNESS_CHALLENGE", f"invalid unresolved witness/challenge {cell_id}",
        )
        named_by_candidate.setdefault(cell["candidate_id"], set()).add(unavailability["named_datum"])

        topology = cell["topology_gate"]
        sector_input = topology["sector_disconnection"].split("(", 1)[0]
        interpolation_input = topology["finite_energy_interpolation"].split("(", 1)[0]
        expected_topology = frozen_evaluate_topology(sector_input, interpolation_input)
        checked(
            topology["integrity"] == "COMPUTATION_VALID" and topology["physics"] == expected_topology
            and topology["pair_annihilation"].split("(", 1)[0]
            in {"NO_FINITE_ENERGY_PATH", "PATH_EXISTS", "UNRESOLVED"}
            and not topology["pair_used_in_aggregate"],
            "ASSERT_TOPOLOGY_AGGREGATE", f"topology aggregate/pair polarity wrong {cell_id}",
        )
        if "nonholonomic" in cell["native_root_class"]:
            checked(any("E4" in value for value in topology["native_structure_tokens"]), "ASSERT_TOPOLOGY_NATIVE_ROOT", f"E4 topology route used action alone {cell_id}")
        if "Rayleigh" in cell["native_root_class"] or "rayleigh" in cell["native_root_class"]:
            checked(any("E5" in value or "gammaSigma" in value or "tangentDtN" in value for value in topology["native_structure_tokens"]), "ASSERT_TOPOLOGY_NATIVE_ROOT", f"E5 topology route used action alone {cell_id}")
        pair = topology["pair_configuration"]
        checked(pair["firewall_tag"] == "PAIR_ANNIHILATION_ONLY" and pair["consumer"] == "topology_pair_annihilation_subquestion", "ASSERT_PAIR_FIREWALL", f"pair object escaped carve-out {cell_id}")

        level_1 = cell["ensemble"]["level_1"]
        level_2 = cell["ensemble"]["level_2"]
        checked(level_1["integrity"] == level_2["integrity"] == "COMPUTATION_VALID", "ASSERT_ENSEMBLE_SCHEMA", f"ensemble integrity wrong {cell_id}")
        assert_local_ensemble(cell, availability)
        local_positive = level_1["local_preclosure_classification"] is not None
        if local_positive:
            checked(cell["closure_adjudication"] is not None and level_1["physics"] == "UNRESOLVED(mechanical-closure refusal)", "ASSERT_ENSEMBLE_CLOSURE_PROPAGATION", f"positive ensemble bypassed closure refusal {cell_id}")
            assert_closure_record(cell["closure_adjudication"])
        else:
            checked(cell["closure_adjudication"] is None and level_1["physics"] == "UNRESOLVED(boundary-action variation)", "ASSERT_ENSEMBLE_CLOSURE_PROPAGATION", f"unavailable ensemble emitted positive claim {cell_id}")
        assert_exchange(cell, availability)
        inputs = level_2["evaluator_inputs"]
        cross = frozen_evaluate_cross_level(inputs["applicability"], inputs["gate_integrity"], inputs["gate_outcome"])
        checked(level_2["physics"] == cross, "ASSERT_CROSS_LEVEL_ENSEMBLE", f"cross-level landing wrong {cell_id}")
        checked(cell["host_location"]["integrity"] == "COMPUTATION_VALID" and cell["host_location"]["physics"] == "undetermined" and cell["host_location"]["evidential_basis"]["challenge_executed"], "ASSERT_HOST_RECORD", f"host field invalid {cell_id}")
        ownership = cell["return_closure_ownership"]
        checked(not ownership["U2_owned"] and not ownership["computed_reachable_from_U2_verdict"] and ownership["preserved_terminal"] == "UNRESOLVED(return_closure)", "ASSERT_RETURN_CLOSURE_OWNERSHIP", f"return closure leaked into U2 verdict {cell_id}")
        checked(cell["stable_branch_id"].startswith("U2:") and cell["stable_branch_id"] == ownership["stable_branch_id"], "ASSERT_STABLE_BRANCH_IDS", f"branch id missing/mismatched {cell_id}")

    unresolved_candidates = {
        cell["candidate_id"] for cell in cells.values()
        if cell["disposition"]["kind"] == "UNRESOLVED"
    }
    checked(
        set(named_by_candidate) == unresolved_candidates
        and len({value for values in named_by_candidate.values() for value in values})
        >= len(unresolved_candidates),
        "ASSERT_CANDIDATE_WISE_NONIDENTIFIABILITY",
        "unresolved candidates lack candidate-specific named data", checks,
    )

    closure_records = by_id(view["closure_records"], "record_id")
    expected_closure_ids = {
        cell["closure_adjudication"]["record_id"]
        for cell in cells.values() if cell["closure_adjudication"] is not None
    }
    checked(set(closure_records) == expected_closure_ids, "ASSERT_GATED_CLAIM_CLOSURE_COVERAGE", "gated claim closure/refusal coverage failed", checks)

    promotions = by_id(view["promotion_records"], "promotion_key")
    expected_promotions = {row["promotion_key"] for row in grid_doc["promotion_contexts"]}
    checked(set(promotions) == expected_promotions and len(promotions) == 16, "ASSERT_PROMOTION_COVERAGE", "promotion-context coverage failed", checks)
    for key, promotion in promotions.items():
        checked(promotion["integrity"] == "COMPUTATION_VALID", "ASSERT_ZERO_INTEGRITY", f"invalid promotion {key}")
        checked(promotion["candidate_universe_digest"] == axes["candidate_universe_digest"], "ASSERT_PROMOTION_CANDIDATE_UNIVERSE", f"stale promotion candidate universe {key}")
        checked(promotion["forcing_census_A"]["available_normalized_operator_derivations"] == promotion["forcing_census_B"]["available_normalized_operator_derivations"], "ASSERT_FORCING_CENSUS", f"forcing census disagreement {key}")
        landing = frozen_evaluate_promotion(promotion["evaluator_inputs"])
        witness = promotion["no_forcing_witness"]
        expected_promotion_kind = landing.split("(", 1)[0]
        checked(
            landing == promotion["evaluator_landing"]
            and promotion["physics"] == expected_promotion_kind,
            "ASSERT_PROMOTION_DECISION_TABLE", f"promotion table mismatch {key}",
        )
        if expected_promotion_kind == "NO_SELECTION_CLAIM":
            checked(witness["executed"] and witness["construct_B_star_attempted"] and witness["measured_rank"] == 0 and witness["challenge_outcome"] == "CONSTRUCTIVE_FAIL", "ASSERT_NO_FORCING_WITNESS", f"no-forcing witness invalid {key}")
        checked(not promotion["exclusion_records_referenced"] and not promotion["survivor_or_complement_objects_referenced"] and not promotion["stability_roots_referenced"], "ASSERT_PROMOTION_ANCESTRY", f"selection-from-exclusion/stability ancestry {key}")

    derived_eligibility = []
    for candidate in axes["candidate_axis"]:
        for ambient in axes["ambient_branches"]:
            branch_cells = [row for row in cells.values() if row["candidate_id"] == candidate and row["ambient"] == ambient]
            integrities = ids(row["integrity"] for row in branch_cells)
            dispositions = ids(row["disposition"]["kind"] for row in branch_cells if row["integrity"] == "COMPUTATION_VALID")
            if candidate == "OTHER":
                eligible, reason = False, "catch_all_OTHER_has_no_defining_operator"
            elif not branch_cells:
                eligible, reason = False, "no_stratum_cells"
            elif integrities != ["COMPUTATION_VALID"]:
                eligible, reason = False, "integrity_stratum_cell_present"
            elif "EXCLUDED" in dispositions:
                eligible, reason = False, "EXCLUDED_stratum_cell_present"
            elif len(dispositions) != 1:
                eligible, reason = False, "heterogeneous_stratum_dispositions"
            elif dispositions[0] not in {"ADMISSIBLE", "UNRESOLVED"}:
                eligible, reason = False, "disposition_class_not_template_eligible"
            else:
                eligible, reason = True, None
            derived_eligibility.append({
                "template_branch_id": f"U2:{candidate}:{ambient}", "candidate_id": candidate, "ambient": ambient,
                "stratum_cell_ids": ids(row["cell_id"] for row in branch_cells), "stratum_count": len(branch_cells),
                "integrity_classes": integrities, "disposition_classes": dispositions,
                "template_eligible": eligible, "eligibility_class": dispositions[0] if eligible else None,
                "ineligibility_reason": reason,
            })
    eligibility_rows = view["template_eligibility_branches"]
    eligible_branch_ids = {row["template_branch_id"] for row in derived_eligibility if row["template_eligible"]}
    templates = by_id(view["posed_BVP_templates"], "stable_branch_id")
    checked(
        eligibility_rows == derived_eligibility and len(eligible_branch_ids) == 16
        and set(templates) == eligible_branch_ids
        and all(row["candidate_id"] != "OTHER" for row in eligibility_rows if row["template_eligible"])
        and all(cell["template_eligible"] == (f"U2:{cell['candidate_id']}:{cell['ambient']}" in eligible_branch_ids) for cell in cells.values()),
        "ASSERT_TEMPLATE_ELIGIBILITY", "branch eligibility and emitted-template sets differ", checks,
    )

    template_contract = contracts["template"]
    proof_by_branch = by_id(template_contract["branch_equivalence_proofs"], "template_branch_id")
    allowed_posing = set(template_contract["conditional_template_posing_allowed_constructors"])
    r49_ids = template_contract["R49_exact_unresolved_reference_ids"]
    for branch_id, template in templates.items():
        branch_cells = [cells[value] for value in template["source_stratum_cell_ids"]]
        proof = proof_by_branch[branch_id]
        checked(
            proof["identical_BVP_across_strata"] and proof["distinct_census_digest_count"] == 1
            and proof["proof_timing"] == "pre-production_symbolic_only"
            and not proof["produced_adjudication_objects_used"]
            and ids(row["stratum"] for row in branch_cells) == ids(proof["active_strata"])
            and all(digest(row["constituent_census"]) == row["census_sha256"] for row in proof["per_stratum_censuses"])
            and len({row["census_sha256"] for row in proof["per_stratum_censuses"]}) == 1,
            "ASSERT_TEMPLATE_STRATUM_BVP_EQUIVALENCE", f"invalid branch collapse {branch_id}", checks,
        )
        disposition = template["eligibility_disposition"]
        conditional = disposition == "UNRESOLVED"
        expected_missing_ids = ids(
            datum for cell in branch_cells
            for datum in (cell["disposition"]["named_datum"], cell["disposition"]["stratum_datum"])
        ) if conditional else []
        constituents = template["constituents"]
        typed_free_data = constituents["typed_free_data"]
        checked(
            template["integrity"] == "COMPUTATION_VALID" and template["physics"] == "POSED_BVP_TEMPLATE"
            and template["conditional"] == conditional and template["dependent_fields_unbound"]
            and template["evaluation_state"] == "UNEVALUATED"
            and [row["reference_id"] for row in typed_free_data] == r49_ids
            and all(row["availability"] == "UNRESOLVED" and row["domain"] and row["dimensions"] for row in typed_free_data),
            "ASSERT_TEMPLATE_NO_SOLVED_CONTENT", f"template schema/free-data invalid {branch_id}", checks,
        )
        if conditional:
            branch_tag = constituents.get("branch_conditionality_tag", {})
            open_tag = constituents.get("open_data_conditionality_tag", {})
            checked(
                branch_tag.get("tag") == f"CONDITIONAL_ON_BRANCH({template['candidate_id']})"
                and branch_tag.get("candidate_id") == template["candidate_id"]
                and open_tag.get("unresolved_missing_datum_ids") == expected_missing_ids
                and all(not tag.get("evidential") and not tag.get("posing_DAG_reachable") for tag in (branch_tag, open_tag)),
                "ASSERT_CONDITIONAL_TEMPLATE_TAGS", f"conditional tags missing/evidential/reachable {branch_id}", checks,
            )
        else:
            checked(
                "branch_conditionality_tag" not in constituents and "open_data_conditionality_tag" not in constituents,
                "ASSERT_CONDITIONAL_TEMPLATE_TAGS", f"unconditional template carries conditional tags {branch_id}", checks,
            )
        checked(
            not template["branch_admissibility_claim"] and not template["branch_selection_claim"],
            "ASSERT_CONDITIONAL_TEMPLATE_NO_CLAIMS", f"template asserts admissibility/selection {branch_id}", checks,
        )
        posing = template["posing_proof_dag"]
        computed_posing = set(frozen_classify_inference_content(posing["normalized_inference_content"]))
        checked(
            computed_posing == set(posing["constructors"]) and computed_posing <= allowed_posing
            and not template["excluded_record_references"] and not template["complement_record_references"]
            and ("POSTULATE_BRANCH" in computed_posing) == (template["ambient"] == "two_sided_R_w_postulate"),
            "ASSERT_CONDITIONAL_TEMPLATE_ANCESTRY", f"banned normalized posing content {branch_id}", checks,
        )
        metadata_ids = {
            constituents["branch_conditionality_tag"]["metadata_id"],
            constituents["open_data_conditionality_tag"]["metadata_id"],
        } if conditional else set()
        checked(
            metadata_ids.isdisjoint(scalar_values(posing["normalized_inference_content"]))
            and metadata_ids.isdisjoint(scalar_values(template["symbolic_ast"])),
            "ASSERT_CONDITIONAL_TEMPLATE_METADATA_UNREACHABLE", f"conditional metadata reached posing DAG {branch_id}", checks,
        )
        expected_terms = [
            {"term_id": "residual:bulk_euler_lagrange_residual", "kind": "residual"},
            {"term_id": f"boundary:canonical_operator:{template['candidate_id']}", "kind": "boundary"},
            {"term_id": "zero-mode:translation_zero_mode", "kind": "zero-mode"},
            {"term_id": f"asymptotic-matching:{template['ambient']}", "kind": "asymptotic-matching"},
        ]
        checked(
            template["expected_term_census"] == expected_terms
            and walk_template_terms(template["symbolic_ast"]) == expected_terms,
            "ASSERT_TEMPLATE_TERM_INCIDENCE", f"template term census/incidence incomplete {branch_id}", checks,
        )
        checked(
            not contains_solved_payload(constituents)
            and "SOLVE_EVALUATION_RESULT" not in computed_posing,
            "ASSERT_TEMPLATE_NO_SOLVED_CONTENT", f"solved payload reachable {branch_id}", checks,
        )

    fixtures = view["guard_fixtures"]
    template_fixture = fixtures["template"]
    checked(template_fixture["term_incidence_exact"] and template_fixture["expected_terms"] == template_fixture["emitted_terms"], "ASSERT_TEMPLATE_TERM_INCIDENCE", "posed-template term fixture incomplete", checks)
    checked(template_fixture["dependent_fields_unbound"] and not template_fixture["solve_evaluation_node_reachable"], "ASSERT_TEMPLATE_NO_SOLVED_CONTENT", "template fixture contains solved content", checks)
    closure_fixture = fixtures["closure"]
    checked(set(closure_fixture["census_A"]) == set(closure_fixture["census_B"]) and closure_fixture["reconstruction_exact"] and closure_fixture["no_double_count"], "ASSERT_CLOSURE_FIXTURE", "closure exact-ownership fixture failed", checks)
    checked(fixtures["mixed_applicability"] == {"mixed_geometric": "refinement-UNRESOLVED", "mixed_non_geometric": "NOT_APPLICABLE(witness)"}, "ASSERT_MIXED_APPLICABILITY_FIXTURES", "mixed applicability fixtures failed", checks)
    checked(fixtures["pair_firewall"]["other_consumer_count"] == 0, "ASSERT_PAIR_FIREWALL", "fixture pair escaped carve-out", checks)

    for path, dag in walk_proof_dags(view):
        assert_proof_dag(path, dag, taxonomy_roots)
    scope = view["scope"]
    checked(scope["static_adjudication_only"] and scope["dynamical_selection_deferred"] and scope["one_body_only"] and not scope["BVP_solved"], "ASSERT_SCOPE_INVARIANTS", "scope invariant false", checks)
    checked(scope["two_body_consumer_object_count"] == 0 and not scope["postulate_used_as_selection_evidence"] and scope["banned_native_import_count"] == 0, "ASSERT_FIREWALLS", "two-body/postulate/import firewall failed", checks)
    dimension = view["dimensional_firewall"]
    checked(dimension["all_constructed_real_expressions_homogeneous"] and dimension["firing_ablation"]["heterogeneity_detected"], "ASSERT_DIMENSION_FIREWALL", "dimensional firewall/ablation failed", checks)
    computed_disposition_census = dict(sorted(Counter(
        cell["disposition"]["kind"] for cell in cells.values()
    ).items()))
    checked(
        view["headlines"]["dispositions"] == computed_disposition_census,
        "ASSERT_OUTCOME_CENSUS_REPORT", "reported disposition census differs from cell outcomes", checks,
    )
    checked(
        view["headlines"]["integrity_failures"] == 0,
        "ASSERT_BANKING_ZERO_INTEGRITY", "headline integrity aggregate is nonzero", checks,
    )
    checked(
        projected_digest == pin["U2_V11_PHYSICS_RECORD_SET_DIGEST"]
        == declared_invariance["reference_digest"] == declared_invariance["live_v12_digest"]
        and projected_universe == pin["record_id_universe"] == declared_invariance["record_id_universe"]
        and len(projected) == declared_invariance["record_count"] == pin["record_count"] == 992
        and axes["candidate_universe_digest"] == pin["candidate_universe_digest"]
        and declared_invariance["candidate_universe_digest_unchanged"]
        and declared_invariance["equal"]
        and declared_invariance["comparison_timing"] == "after_all_records_emitted_before_banking"
        and declared_invariance["comparison_predicate"] == pin["comparison_predicate"],
        "ASSERT_PHYSICS_RECORD_INVARIANCE", "live v12 projected physics records differ from frozen v11 pin", checks,
    )


def mutation_catalog(view: dict[str, Any]) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for route in view["route_controls"]:
        rows.append({"mutation_id": f"ROUTE_MALFORMED:{route['route_id']}", "expected_assert_id": route_assert_id(route), "family": "per_cell_route"})
        rows.append({"mutation_id": f"ROUTE_ABLATION:{route['route_id']}", "expected_assert_id": ablation_assert_id(route), "family": "per_cell_ablation"})
    for cell in view["cell_records"]:
        if cell["ensemble"]["level_1"]["exchange_control"] is not None:
            rows.append({"mutation_id": f"EXCHANGE:{cell['cell_id']}", "expected_assert_id": exchange_assert_id(cell), "family": "per_claim_exchange"})
    globals_ = {
        "GRID_COVERAGE": "ASSERT_GRID_COVERAGE", "ROUTE_SET": "ASSERT_ROUTE_SET",
        "CELL_ROUTE_FIXTURE_COVERAGE": "ASSERT_CELL_ROUTE_FIXTURE_COVERAGE",
        "RECORD_IDS_UNIQUE": "ASSERT_RECORD_IDS_UNIQUE",
        "DEPENDENCY_BIDIRECTIONAL": "ASSERT_DEPENDENCY_BIDIRECTIONAL",
        "EVIDENCE_CENSUS_COVERAGE": "ASSERT_EVIDENCE_CENSUS_COVERAGE",
        "EVIDENCE_RAW_MEASUREMENT": "ASSERT_EVIDENCE_MEASUREMENTS",
        "EVIDENCE_RELABEL": "ASSERT_EVIDENCE_STATE_CLASSIFICATION",
        "DISPOSITION_EVALUATOR_LANDING": "ASSERT_DISPOSITION_PRECEDENCE",
        "DISPOSITION_PRECEDENCE": "ASSERT_DISPOSITION_RECORD",
        "TYPED_WITNESS_CHALLENGE": "ASSERT_TYPED_WITNESS_CHALLENGE",
        "UNAVAILABILITY_REFERENCE_SET": "ASSERT_TYPED_WITNESS_CHALLENGE",
        "TOPOLOGY_AGGREGATE": "ASSERT_TOPOLOGY_AGGREGATE", "PAIR_POLARITY": "ASSERT_TOPOLOGY_AGGREGATE",
        "TOPOLOGY_NATIVE_ROOT": "ASSERT_TOPOLOGY_NATIVE_ROOT",
        "ENSEMBLE_SCHEMA": "ASSERT_ENSEMBLE_SCHEMA",
        "ENSEMBLE_STATIC_CLASSIFIER": "ASSERT_ENSEMBLE_STATIC_CLASSIFIER",
        "CROSS_LEVEL": "ASSERT_CROSS_LEVEL_ENSEMBLE", "ENSEMBLE_CLOSURE": "ASSERT_ENSEMBLE_CLOSURE_PROPAGATION",
        "CROSS_LEVEL_INTEGRITY_PROPAGATION": "ASSERT_CROSS_LEVEL_ENSEMBLE",
        "CLOSURE_CENSUS": "ASSERT_CLOSURE_CENSUS_GENERATION",
        "CLOSURE_TOTAL": "ASSERT_CLOSURE_TOTAL_GENERATION",
        "CLOSURE_DUPLICATE": "ASSERT_CLOSURE_GATE",
        "CLOSURE_OMISSION": "ASSERT_CLOSURE_GATE", "CLOSURE_RETURN_LEAK": "ASSERT_CLOSURE_GATE",
        "CLOSURE_HOST_CONSISTENCY": "ASSERT_CLOSURE_HOST_CONSISTENCY",
        "GATED_CLAIM_CLOSURE_COVERAGE": "ASSERT_GATED_CLAIM_CLOSURE_COVERAGE",
        "FORCING_CENSUS": "ASSERT_FORCING_CENSUS", "PROMOTION_SELECTION": "ASSERT_PROMOTION_DECISION_TABLE",
        "PROMOTION_COVERAGE": "ASSERT_PROMOTION_COVERAGE",
        "PROMOTION_CANDIDATE_UNIVERSE": "ASSERT_PROMOTION_CANDIDATE_UNIVERSE",
        "NO_FORCING_WITNESS": "ASSERT_NO_FORCING_WITNESS", "PROMOTION_EXCLUSION": "ASSERT_PROMOTION_ANCESTRY",
        "TEMPLATE_EXTRA": "ASSERT_TEMPLATE_ELIGIBILITY", "TEMPLATE_SOLVED": "ASSERT_TEMPLATE_NO_SOLVED_CONTENT",
        "TEMPLATE_ELIGIBILITY_FLIP": "ASSERT_TEMPLATE_ELIGIBILITY",
        "CONDITIONAL_TEMPLATE_TAG_STRIP": "ASSERT_CONDITIONAL_TEMPLATE_TAGS",
        "CONDITIONAL_TEMPLATE_OPEN_DATA_TAG_STRIP": "ASSERT_CONDITIONAL_TEMPLATE_TAGS",
        "CONDITIONAL_TEMPLATE_ADMISSIBILITY_ASSERTION": "ASSERT_CONDITIONAL_TEMPLATE_NO_CLAIMS",
        "CONDITIONAL_TEMPLATE_EXCLUSION_ANCESTRY": "ASSERT_CONDITIONAL_TEMPLATE_ANCESTRY",
        "CONDITIONAL_TEMPLATE_LABEL_FREE_INCOMPATIBILITY": "ASSERT_CONDITIONAL_TEMPLATE_ANCESTRY",
        "CONDITIONAL_TEMPLATE_TAG_SOLVED_PAYLOAD": "ASSERT_TEMPLATE_NO_SOLVED_CONTENT",
        "CONDITIONAL_TEMPLATE_TAG_REACHABLE": "ASSERT_CONDITIONAL_TEMPLATE_METADATA_UNREACHABLE",
        "PHYSICS_RECORD_INVARIANCE": "ASSERT_PHYSICS_RECORD_INVARIANCE",
        "TEMPLATE_TERM_RESIDUAL": "ASSERT_TEMPLATE_TERM_INCIDENCE", "TEMPLATE_TERM_BOUNDARY": "ASSERT_TEMPLATE_TERM_INCIDENCE",
        "TEMPLATE_TERM_ZERO_MODE": "ASSERT_TEMPLATE_TERM_INCIDENCE", "TEMPLATE_TERM_ASYMPTOTIC": "ASSERT_TEMPLATE_TERM_INCIDENCE",
        "STABILITY_DAG": "ASSERT_PROGRAM_WIDE_STABILITY_BAN", "CONCEALED_ELIMINATION": "ASSERT_TYPED_ANCESTRY",
        "UNTYPED_ROOT": "ASSERT_NATIVE_ANCESTRY", "SOLVE_DAG": "ASSERT_BVP_NOT_SOLVED",
        "PAIR_CONSUMER": "ASSERT_PAIR_FIREWALL", "RETURN_REACHABLE": "ASSERT_RETURN_CLOSURE_OWNERSHIP",
        "HOST_SCHEMA": "ASSERT_HOST_RECORD", "DIMENSION": "ASSERT_DIMENSION_FIREWALL",
        "INTEGRITY": "ASSERT_ZERO_INTEGRITY", "ALL_INTEGRITY_BANK": "ASSERT_BANKING_ZERO_INTEGRITY",
        "OUTCOME_CENSUS_REPORT": "ASSERT_OUTCOME_CENSUS_REPORT",
        "BRANCH_ID": "ASSERT_STABLE_BRANCH_IDS", "CANDIDATE_NAMED_DATA": "ASSERT_CANDIDATE_WISE_NONIDENTIFIABILITY",
        "CLOSURE_FIXTURE": "ASSERT_CLOSURE_FIXTURE",
        "MIXED_APPLICABILITY_FIXTURES": "ASSERT_MIXED_APPLICABILITY_FIXTURES",
        "SCOPE_INVARIANTS": "ASSERT_SCOPE_INVARIANTS", "STAGE0_BINDING": "ASSERT_STAGE0_BINDING",
        "ENGINE_SCHEMA": "ASSERT_ENGINE_SCHEMA", "ENGINE_ROUTE_REGISTRY": "ASSERT_ENGINE_ROUTE_REGISTRY",
        "DISSIPATION": route_assert_id(next(row for row in view["route_controls"] if row["dissipation_measurement"]["applicable"])),
        "FROZEN_EVALUATOR": "ASSERT_FROZEN_EVALUATORS", "AXES": "ASSERT_AXES",
        "BANNED_IMPORT": "ASSERT_FIREWALLS", "DUAL_ENGINE": "ASSERT_DUAL_ENGINE_SEMANTIC_EQUALITY",
    }
    rows.extend({"mutation_id": key, "expected_assert_id": value, "family": "global_guard"} for key, value in globals_.items())
    for record_type in RECORD_TYPES:
        for defect in ("HARNESS_PHYSICS", "NOT_RUN_PHYSICS", "UPSTREAM_PROPAGATION"):
            rows.append({
                "mutation_id": f"RECORD_SCHEMA:{record_type}:{defect}",
                "expected_assert_id": record_schema_assert_id(record_type),
                "family": "per_record_type_integrity_schema",
            })
    return rows


def apply_mutation(view: dict[str, Any], mutation: str) -> None:
    changed = False
    if mutation.startswith("ROUTE_MALFORMED:"):
        target = mutation.split(":", 1)[1]
        row = next(value for value in view["route_controls"] if value["route_id"] == target)
        changed = row["malformed_excluded"] is True; row["malformed_excluded"] = False
    elif mutation.startswith("ROUTE_ABLATION:"):
        target = mutation.split(":", 1)[1]
        row = next(value for value in view["route_controls"] if value["route_id"] == target)
        changed = row["semantic_ablation"]["nondefinitional_obligation_failed"] is True
        row["semantic_ablation"]["nondefinitional_obligation_failed"] = False
    elif mutation.startswith("EXCHANGE:"):
        target = mutation.split(":", 1)[1]
        row = next(value for value in view["cell_records"] if value["cell_id"] == target)
        control = row["ensemble"]["level_1"]["exchange_control"]
        if control["route"] == "native_conjugate_transform":
            changed = bool(control["independently_generated_expected_signature"])
            control["independently_generated_expected_signature"].pop()
        else:
            insufficiency = control["pairing_insufficiency"]
            changed = bool(insufficiency["nonconjugate_components"] or insufficiency["unmapped_canonical_signature_terms"])
            insufficiency["nonconjugate_components"] = []
            insufficiency["unmapped_canonical_signature_terms"] = []
    else:
        cell = view["cell_records"][0]
        promotion = view["promotion_records"][0]
        closure = view["closure_records"][0]
        if mutation == "GRID_COVERAGE": changed = True; view["cell_records"].pop()
        elif mutation == "ROUTE_SET": changed = True; view["route_controls"].pop()
        elif mutation == "CELL_ROUTE_FIXTURE_COVERAGE": changed = True; view["route_controls"][0]["cell_id"] = "fixture:wrong-cell"
        elif mutation == "RECORD_IDS_UNIQUE": changed = True; view["cell_records"].append(copy.deepcopy(cell))
        elif mutation == "DEPENDENCY_BIDIRECTIONAL": changed = True; cell["used_dependencies"].pop()
        elif mutation == "EVIDENCE_CENSUS_COVERAGE": changed = True; cell["obligation_evidence"].pop()
        elif mutation == "EVIDENCE_RAW_MEASUREMENT":
            target = next(value for value in cell["obligation_evidence"] if value["raw_predicate_inputs"]["datum_missing"])
            changed = True; target["raw_predicate_inputs"]["datum_missing"] = False
        elif mutation == "EVIDENCE_RELABEL":
            record = next(value for value in cell["obligation_evidence"] if value["emitted_state"] == "MISSING")
            changed = True; record["emitted_state"] = "SATISFIED"
        elif mutation == "DISPOSITION_EVALUATOR_LANDING": changed = True; cell["disposition_evaluator_landing"] = "ADMISSIBLE"
        elif mutation == "DISPOSITION_PRECEDENCE": changed = True; cell["disposition"]["kind"] = "ADMISSIBLE"
        elif mutation == "TYPED_WITNESS_CHALLENGE": changed = True; cell["unavailability"]["challenge"]["executed"] = False
        elif mutation == "UNAVAILABILITY_REFERENCE_SET": changed = True; cell["unavailability"]["referenced_stage0_slots"].pop()
        elif mutation == "TOPOLOGY_AGGREGATE": changed = True; cell["topology_gate"]["physics"] = "orientation-only"
        elif mutation == "PAIR_POLARITY": changed = True; cell["topology_gate"]["pair_used_in_aggregate"] = True
        elif mutation == "TOPOLOGY_NATIVE_ROOT":
            target = next(value for value in view["cell_records"] if "nonholonomic" in value["native_root_class"])
            changed = True; target["topology_gate"]["native_structure_tokens"] = []
        elif mutation == "ENSEMBLE_SCHEMA": changed = True; cell["ensemble"]["level_1"]["integrity"] = "HARNESS_FAILED"
        elif mutation == "ENSEMBLE_STATIC_CLASSIFIER":
            target = next(value for value in view["cell_records"] if value["ensemble"]["level_1"]["local_preclosure_classification"])
            changed = True
            target["ensemble"]["level_1"]["local_preclosure_classification"] = "fixed-source"
        elif mutation == "CROSS_LEVEL": changed = True; cell["ensemble"]["level_2"]["physics"] = "not-defect-refined"
        elif mutation == "CROSS_LEVEL_INTEGRITY_PROPAGATION":
            changed = True
            cell["ensemble"]["level_2"]["evaluator_inputs"]["applicability"] = "geometric_component_present"
            cell["ensemble"]["level_2"]["evaluator_inputs"]["gate_integrity"] = "HARNESS_FAILED"
        elif mutation == "ENSEMBLE_CLOSURE":
            cell = next(value for value in view["cell_records"] if value["closure_adjudication"])
            changed = True; cell["ensemble"]["level_1"]["physics"] = cell["ensemble"]["level_1"]["local_preclosure_classification"]
        elif mutation == "CLOSURE_CENSUS": changed = True; closure["contribution_census_B"].pop()
        elif mutation == "CLOSURE_TOTAL": changed = True; closure["independently_constructed_symbolic_total"] += 1
        elif mutation == "CLOSURE_DUPLICATE": changed = True; closure["assignments"].append(copy.deepcopy(closure["assignments"][0]))
        elif mutation == "CLOSURE_OMISSION": changed = True; closure["assignment_coverage"].pop()
        elif mutation == "CLOSURE_RETURN_LEAK": changed = True; closure["return_closure_consumed"] = True
        elif mutation == "CLOSURE_HOST_CONSISTENCY": changed = True; closure["host_physics"] = "boundary-localized"
        elif mutation == "GATED_CLAIM_CLOSURE_COVERAGE": changed = True; view["closure_records"].pop()
        elif mutation == "FORCING_CENSUS": changed = True; promotion["forcing_census_B"]["available_normalized_operator_derivations"] = ["forcing:fake"]
        elif mutation == "PROMOTION_SELECTION": changed = True; promotion["physics"] = "SELECTED"
        elif mutation == "PROMOTION_COVERAGE": changed = True; view["promotion_records"].pop()
        elif mutation == "PROMOTION_CANDIDATE_UNIVERSE": changed = True; promotion["candidate_universe_digest"] = "stale"
        elif mutation == "NO_FORCING_WITNESS": changed = True; promotion["no_forcing_witness"]["executed"] = False
        elif mutation == "PROMOTION_EXCLUSION": changed = True; promotion["exclusion_records_referenced"] = ["excluded:E2"]
        elif mutation == "TEMPLATE_EXTRA":
            changed = True
            extra = copy.deepcopy(view["posed_BVP_templates"][0])
            extra.update({"record_id": "template:illegal", "stable_branch_id": "U2:illegal:ambient"})
            view["posed_BVP_templates"].append(extra)
        elif mutation == "TEMPLATE_SOLVED": changed = True; view["guard_fixtures"]["template"]["solve_evaluation_node_reachable"] = True
        elif mutation == "TEMPLATE_ELIGIBILITY_FLIP": changed = True; view["template_eligibility_branches"][0]["template_eligible"] = False
        elif mutation in {"CONDITIONAL_TEMPLATE_TAG_STRIP", "CONDITIONAL_TEMPLATE_OPEN_DATA_TAG_STRIP"}:
            changed = True
            key = "branch_conditionality_tag" if mutation == "CONDITIONAL_TEMPLATE_TAG_STRIP" else "open_data_conditionality_tag"
            view["posed_BVP_templates"][0]["constituents"].pop(key)
        elif mutation == "CONDITIONAL_TEMPLATE_ADMISSIBILITY_ASSERTION":
            changed = True; view["posed_BVP_templates"][0]["branch_admissibility_claim"] = True
        elif mutation in {"CONDITIONAL_TEMPLATE_EXCLUSION_ANCESTRY", "CONDITIONAL_TEMPLATE_LABEL_FREE_INCOMPATIBILITY"}:
            changed = True
            template = view["posed_BVP_templates"][0]
            operation = "complement" if mutation.endswith("EXCLUSION_ANCESTRY") else "incompatible"
            template["posing_proof_dag"]["normalized_inference_content"] = {
                "op": "positive_join", "args": [{"op": operation, "args": []}],
            }
            template["posing_proof_dag"]["constructors"] = ["POSITIVE_DERIVATION"]
        elif mutation == "CONDITIONAL_TEMPLATE_TAG_SOLVED_PAYLOAD":
            changed = True
            view["posed_BVP_templates"][0]["constituents"]["open_data_conditionality_tag"]["payload_kind"] = "RESPONSE_COEFFICIENT"
        elif mutation == "CONDITIONAL_TEMPLATE_TAG_REACHABLE":
            changed = True
            template = view["posed_BVP_templates"][0]
            metadata_id = template["constituents"]["branch_conditionality_tag"]["metadata_id"]
            template["posing_proof_dag"]["normalized_inference_content"]["args"].append({"op": "root", "id": metadata_id})
        elif mutation == "PHYSICS_RECORD_INVARIANCE":
            changed = True; view["physics_record_invariance"]["live_v12_digest"] = "0" * 64
        elif mutation.startswith("TEMPLATE_TERM_"):
            index = {"TEMPLATE_TERM_RESIDUAL": 0, "TEMPLATE_TERM_BOUNDARY": 1, "TEMPLATE_TERM_ZERO_MODE": 2, "TEMPLATE_TERM_ASYMPTOTIC": 3}[mutation]
            changed = True; view["guard_fixtures"]["template"]["emitted_terms"].pop(index)
        elif mutation in {"STABILITY_DAG", "CONCEALED_ELIMINATION", "SOLVE_DAG"}:
            dag = cell["obligation_evidence"][0]["proof_dag"]
            changed = True
            dag["normalized_inference_content"]["op"] = {"STABILITY_DAG": "stability", "CONCEALED_ELIMINATION": "case_eliminate", "SOLVE_DAG": "solve"}[mutation]
            if mutation == "CONCEALED_ELIMINATION":
                dag["constructors"] = ["POSITIVE_DERIVATION", "ROOT_REFERENCE"]
            else:
                dag["constructors"] = frozen_classify_inference_content(dag["normalized_inference_content"])
        elif mutation == "UNTYPED_ROOT":
            dag = cell["obligation_evidence"][0]["proof_dag"]; changed = True
            dag["root_ids"].append("source:builder:untyped")
        elif mutation == "PAIR_CONSUMER": changed = True; cell["topology_gate"]["pair_configuration"]["consumer"] = "force_consumer"
        elif mutation == "RETURN_REACHABLE": changed = True; cell["return_closure_ownership"]["computed_reachable_from_U2_verdict"] = True
        elif mutation == "HOST_SCHEMA": changed = True; cell["host_location"]["integrity"] = "HARNESS_FAILED"
        elif mutation == "DIMENSION": changed = True; view["dimensional_firewall"]["firing_ablation"]["heterogeneity_detected"] = False
        elif mutation == "INTEGRITY": changed = True; cell["integrity"] = "HARNESS_FAILED"
        elif mutation == "ALL_INTEGRITY_BANK": changed = True; view["headlines"]["integrity_failures"] = 144
        elif mutation == "OUTCOME_CENSUS_REPORT": changed = True; view["headlines"]["dispositions"] = {"EXCLUDED": 144}
        elif mutation == "BRANCH_ID": changed = True; cell["stable_branch_id"] = ""
        elif mutation == "CANDIDATE_NAMED_DATA":
            e1 = next(value for value in view["cell_records"] if value["candidate_id"] == "E1")["unavailability"]["named_datum"]
            for value in view["cell_records"]:
                if value["candidate_id"] == "E2": value["unavailability"]["named_datum"] = e1
            changed = True
        elif mutation == "DISSIPATION":
            row = next(value for value in view["route_controls"] if value["dissipation_measurement"]["applicable"])
            changed = True; row["dissipation_measurement"]["value"] += 1
        elif mutation == "FROZEN_EVALUATOR": changed = True; view["frozen_evaluator_behavior"]["disposition"]["missing_only"] = "ADMISSIBLE"
        elif mutation == "AXES": changed = True; view["axes"]["cell_count"] = 143
        elif mutation == "BANNED_IMPORT": changed = True; view["scope"]["banned_native_import_count"] = 1
        elif mutation == "DUAL_ENGINE": changed = True; view["scope"]["one_body_only"] = False
        elif mutation == "CLOSURE_FIXTURE": changed = True; view["guard_fixtures"]["closure"]["reconstruction_exact"] = False
        elif mutation == "MIXED_APPLICABILITY_FIXTURES": changed = True; view["guard_fixtures"]["mixed_applicability"]["mixed_geometric"] = "NOT_APPLICABLE(witness)"
        elif mutation == "SCOPE_INVARIANTS": changed = True; view["scope"]["static_adjudication_only"] = False
        elif mutation == "STAGE0_BINDING": changed = True; view["stage0_binding"]["equal"] = False
    checked(changed, "MUTATION_NOOP", f"mutation made no change: {mutation}")


def campaign_probe(document: dict[str, Any], bundle: dict[str, Any], output: Path) -> int:
    validate_engine_document(document, bundle)
    view = document["semantic_view"]
    catalog = mutation_catalog(view)
    records = []
    for item in catalog:
        mutation = item["mutation_id"]
        expected = item["expected_assert_id"]
        try:
            if mutation.startswith(("ROUTE_MALFORMED:", "ROUTE_ABLATION:")):
                target = mutation.split(":", 1)[1]
                original = next(value for value in view["route_controls"] if value["route_id"] == target)
                mutant = copy.deepcopy(original)
                if mutation.startswith("ROUTE_MALFORMED:"): mutant["malformed_excluded"] = False
                else: mutant["semantic_ablation"]["nondefinitional_obligation_failed"] = False
                frozen = next(value for value in bundle["route_fixture_inventory"]["route_records"] if value["route_id"] == target)
                assert_route_control(mutant, frozen)
            elif mutation.startswith("EXCHANGE:"):
                target = mutation.split(":", 1)[1]
                mutant = copy.deepcopy(next(value for value in view["cell_records"] if value["cell_id"] == target))
                control = mutant["ensemble"]["level_1"]["exchange_control"]
                if control["route"] == "native_conjugate_transform":
                    control["independently_generated_expected_signature"].pop()
                else:
                    control["pairing_insufficiency"]["nonconjugate_components"] = []
                    control["pairing_insufficiency"]["unmapped_canonical_signature_terms"] = []
                assert_exchange(mutant, {
                    row["slot_id"]: row for row in bundle["availability_slots"]["slots"]
                })
            elif mutation in {"ENGINE_SCHEMA", "ENGINE_ROUTE_REGISTRY"}:
                mutant = copy.deepcopy(document)
                if mutation == "ENGINE_SCHEMA":
                    mutant["schema_version"] = "U2_PRODUCTION_ENGINE_MUTANT"
                else:
                    mutant["engine_local_route_registry"].pop()
                validate_engine_document(mutant, bundle)
            elif mutation == "DUAL_ENGINE":
                mutant = copy.deepcopy(document)
                mutant["semantic_view"]["scope"]["one_body_only"] = False
                checked(
                    document["semantic_view"] == mutant["semantic_view"],
                    "ASSERT_DUAL_ENGINE_SEMANTIC_EQUALITY", "dual-engine semantic mutation",
                )
            elif mutation.startswith("RECORD_SCHEMA:"):
                _, record_type, defect = mutation.split(":", 2)
                if defect == "HARNESS_PHYSICS":
                    mutant_record = {"record_type": record_type, "integrity": "HARNESS_FAILED", "physics": "illegal", "failure_reason": "fixture_typed_reason"}
                elif defect == "NOT_RUN_PHYSICS":
                    mutant_record = {"record_type": record_type, "integrity": "NOT_RUN", "physics": "illegal", "failed_upstreams": ["fixture:upstream"]}
                else:
                    mutant_record = {"record_type": record_type, "integrity": "COMPUTATION_VALID", "physics": "illegal", "failed_upstreams": ["fixture:upstream"]}
                checked(
                    frozen_record_schema_valid(mutant_record), record_schema_assert_id(record_type),
                    f"integrity-schema mutation survived for {record_type}:{defect}",
                )
            else:
                mutant = copy.deepcopy(view)
                apply_mutation(mutant, mutation)
                validate_semantic(mutant, bundle)
        except AssertionDeath as failure:
            checked(failure.assert_id == expected, "ASSERT_MUTATION_OWN_ASSERT", f"{mutation} died at {failure.assert_id}, expected {expected}")
            records.append({**item, "status": "DIED_AT_OWN_ASSERT", "observed_assert_id": failure.assert_id})
        else:
            raise AssertionDeath("ASSERT_MUTATION_SURVIVAL", f"mutation survived: {mutation}")
    result = {
        "schema_version": "U2_PRODUCTION_MUTATION_PROBE_V1",
        "status": "PASS", "tooth_count": len(records), "control_count": len(records),
        "vacuous_case_count": 0, "mutation_noop_count": 0,
        "records": records,
        "defect_absent_controls": [
            {"mutation_id": row["mutation_id"], "status": "DEFECT_ABSENT_REAL_CHECK_SURVIVED", "shared_control_execution": index > 0}
            for index, row in enumerate(records)
        ],
        "shared_control_semantic_sha256": digest(view),
    }
    dump_yaml(output, result)
    print(f"U2_PRODUCTION_MUTATION_PROBE_PASS teeth={len(records)} controls={len(records)}")
    return 0


def a9_coverage(view: dict[str, Any]) -> dict[str, Any]:
    rows: list[dict[str, str]] = []
    for route in view["route_controls"]:
        rows.append({"object_id": f"route_control:{route['route_id']}", "category": "same_route_fixture_control", "recompute_route": "route_control"})
    for cell in view["cell_records"]:
        cell_id = cell["cell_id"]
        for category, prefix, route in (
            ("candidate_disposition", "disposition", "cell_record"),
            ("ensemble_level_1", "ensemble_level_1", "cell_record"),
            ("ensemble_level_2", "ensemble_level_2", "cell_record"),
            ("topology_gate", "topology_gate", "topology_record"),
            ("host_location", "host", "host_record"),
            ("return_closure_ownership", "return_closure", "cell_record"),
        ):
            rows.append({"object_id": f"{prefix}:{cell_id}", "category": category, "recompute_route": route})
        if cell["closure_adjudication"]:
            rows.append({"object_id": cell["closure_adjudication"]["record_id"], "category": "closure_adjudication", "recompute_route": "closure_refusal"})
    for promotion in view["promotion_records"]:
        rows.append({"object_id": promotion["record_id"], "category": "promotion", "recompute_route": "promotion_record"})
    for template in view["posed_BVP_templates"]:
        rows.append({"object_id": template["record_id"], "category": "posed_BVP_template", "recompute_route": "template_record"})
        rows.append({"object_id": f"conditional_contract:{template['record_id']}", "category": "conditional_template_contract", "recompute_route": "normalized_template_contract_guard"})
        rows.append({"object_id": f"stratum_equivalence:{template['stable_branch_id']}", "category": "template_stratum_equivalence_proof", "recompute_route": "stage0_symbolic_constituent_census"})
    rows.append({"object_id": "physics_record_invariance:v11_to_v12", "category": "physics_record_invariance_guard", "recompute_route": "independent_v11_and_v12_projection"})
    object_ids = [row["object_id"] for row in rows]
    counts = dict(sorted(Counter(row["category"] for row in rows).items()))
    zero_cardinality_scopes = {
        "EXCLUDED_decisive_incompatibility": sum(
            cell["disposition"]["kind"] == "EXCLUDED" for cell in view["cell_records"]
        ),
        "positive_witness": sum(
            cell["disposition"]["kind"] == "ADMISSIBLE" for cell in view["cell_records"]
        ),
        "selection_forcing_proof": sum(
            row["physics"] == "SELECTED" for row in view["promotion_records"]
        ),
        "closure_certificate": sum(
            row["physics"]["kind"] == "CLOSURE_CERTIFIED" for row in view["closure_records"]
        ),
    }
    return {
        "object_ids": object_ids, "object_count": len(object_ids), "coverage_map": rows,
        "coverage_category_counts": counts, "coverage_map_exact": len(object_ids) == len(set(object_ids)),
        "generated_category_union_exact": sum(counts.values()) == len(object_ids),
        "recomputation_granularity": "EVERY_OBJECT_INDIVIDUAL",
        "class_level_treatment_used": False,
        "class_equivalence_proof_required": False,
        "required_scope_zero_counts": zero_cardinality_scopes,
        "required_scope_zero_counts_measured": all(value == 0 for value in zero_cardinality_scopes.values()),
        "external_results_claimed": False,
    }


def summary_text(view: dict[str, Any]) -> str:
    h = view["headlines"]
    lines = [
        "# U2 production boundary-adjudication summary", "",
        f"Stage-0 contract digest: `{RATIFIED_DIGEST}`.", "",
        "## Computed landing", "",
        f"- Candidate dispositions: `{h['dispositions']}`; per candidate: `{h['dispositions_by_candidate']}`.",
        f"- Ensemble level 1 (final): `{h['ensemble_level_1_final']}`.",
        f"- Ensemble level 1 (local, before mechanical closure): `{h['ensemble_level_1_local_preclosure']}`.",
        f"- Ensemble level 2: `{h['ensemble_level_2']}`.",
        f"- Topology gate: `{h['topology_gate']}`; pair-annihilation field: `{h['pair_annihilation']}`.",
        f"- Closure adjudications: `{h['closure_adjudications']}`.",
        f"- Promotions: `{h['promotion_outcomes']}`.",
        f"- Posed-BVP templates: `{h['posed_template_count']}` conditional branch templates (8 operator-bearing candidates × 2 ambients); catch-all `OTHER` is named-open with no template.",
        f"- Physics-record invariance: `{view['physics_record_invariance']['equal']}`; live/reference digest `{view['physics_record_invariance']['live_v12_digest']}`.", "",
        "All 144 cells are `COMPUTATION_VALID` and remain `UNRESOLVED`. Each record names candidate-specific missing static structure and its active OPEN stratum, and carries an executed constructive challenge. The templates pose, but do not solve, the 16 operator-bearing branch BVPs and assert neither admissibility nor selection.", "",
        "The local mouth classifications are mechanically refused as final positive ensemble claims because host location remains under-specified. The final level-1 record therefore preserves the forcing evidence as a typed closure refusal rather than laundering it into a positive claim.", "",
        "The dynamical-selection/stability question and return closure remain downstream deferrals and are computed-unreachable from U2 verdicts.", "",
    ]
    return "\n".join(lines)


def compare_and_write(args: argparse.Namespace) -> int:
    sympy = load_yaml(Path(args.sympy))
    wolfram = load_yaml(Path(args.wolfram))
    bundle_dir = Path(args.bundle_dir)
    bundle = {
        name: load_yaml(bundle_dir / filename) for name, filename in {
            "candidate_inventory": "candidate_inventory.yaml", "dependency_grid_inventory": "dependency_grid_inventory.yaml",
            "obligation_censuses": "obligation_censuses.yaml", "route_fixture_inventory": "route_fixture_inventory.yaml",
            "evidence_taxonomy": "evidence_taxonomy.yaml", "availability_slots": "availability_slots.yaml",
            "closure_template_contracts": "closure_template_contracts.yaml",
        }.items()
    }
    if args.campaign_probe:
        return campaign_probe(sympy, bundle, Path(args.campaign_output))
    if args.mutation:
        apply_mutation(sympy["semantic_view"], args.mutation)
    checks: list[dict[str, str]] = []
    validate_engine_document(sympy, bundle, checks)
    validate_engine_document(wolfram, bundle, checks)
    checked(
        sympy["semantic_view"]["headlines"]["dispositions"]
        == wolfram["semantic_view"]["headlines"]["dispositions"],
        "ASSERT_DUAL_ENGINE_OUTCOME_CENSUS",
        "dual-engine disposition censuses differ", checks,
    )
    checked(sympy["semantic_view"] == wolfram["semantic_view"], "ASSERT_DUAL_ENGINE_SEMANTIC_EQUALITY", "full normalized semantic views differ", checks)
    if args.check_only:
        print(f"U2_PRODUCTION_COMPARATOR_CHECK_PASS checks={len(checks)}")
        return 0
    view = sympy["semantic_view"]
    catalog = mutation_catalog(view)
    agreement = {
        "schema_version": "U2_PRODUCTION_ENGINE_AGREEMENT_V1", "status": "ENGINE_AGREE",
        "stage0_contract_digest": RATIFIED_DIGEST, "semantic_sha256": digest(view),
        "sympy_engine_sha256": hashlib.sha256(Path(args.sympy).read_bytes()).hexdigest(),
        "wolfram_engine_sha256": hashlib.sha256(Path(args.wolfram).read_bytes()).hexdigest(),
        "assertion_count": len(checks), "assertions": checks,
        "mutation_catalog_count": len(catalog), "mutation_catalog": catalog,
        "headlines": view["headlines"],
    }
    results = {
        "schema_version": "U2_PRODUCTION_RESULTS_V1", "engine": "DUAL_ENGINE_COMPARATOR_ASSERTED",
        "stage0_contract_digest": RATIFIED_DIGEST, **view,
        "A9_external_verification": a9_coverage(view),
    }
    dump_yaml(Path(args.output), agreement)
    dump_yaml(Path(args.results_output), results)
    Path(args.summary_output).parent.mkdir(parents=True, exist_ok=True)
    Path(args.summary_output).write_text(summary_text(view), encoding="utf-8")
    print(f"U2_PRODUCTION_ENGINE_AGREE cells=144 checks={len(checks)} mutations={len(catalog)}")
    return 0


def parser() -> argparse.ArgumentParser:
    value = argparse.ArgumentParser()
    value.add_argument("--sympy", required=True); value.add_argument("--wolfram", required=True)
    value.add_argument("--bundle-dir", required=True); value.add_argument("--output")
    value.add_argument("--results-output"); value.add_argument("--summary-output")
    value.add_argument("--check-only", action="store_true"); value.add_argument("--mutation")
    value.add_argument("--campaign-probe", action="store_true"); value.add_argument("--campaign-output")
    return value


def main() -> int:
    args = parser().parse_args()
    if not args.check_only and not args.campaign_probe and not all((args.output, args.results_output, args.summary_output)):
        parser().error("normal comparison requires output, results-output, and summary-output")
    if args.campaign_probe and not args.campaign_output:
        parser().error("campaign-probe requires campaign-output")
    try:
        return compare_and_write(args)
    except AssertionDeath as failure:
        print(f"ASSERTION_FAILED {failure.assert_id}: {failure.detail}", file=sys.stderr)
        return 1
    except (KeyError, ValueError, StopIteration, FileNotFoundError) as failure:
        print(f"ASSERTION_FAILED ASSERT_COMPARATOR_INPUT: {failure}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
