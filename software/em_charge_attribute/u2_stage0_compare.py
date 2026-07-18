#!/usr/bin/env python3
"""Strict U2 stage-0 engine comparator and semantic mutation target."""

from __future__ import annotations

import argparse
import copy
import hashlib
import json
import sys
from pathlib import Path
from typing import Any

import yaml


SCHEMA = "U2_STAGE0_ENGINE_V1"
BASE_ENDPOINTS = {"E1", "E2", "E3", "E4", "E5"}
TILT_TYPES = {
    "indexed_density_tilt_profile", "indexed_flow_tilt_response",
    "indexed_h_tilt_profile", "indexed_phase_tilt_profile",
    "indexed_shear_tilt_profile", "indexed_sleeve_surface_normal_profile",
    "indexed_sleeve_tilt_profile", "indexed_uw_tilt_profile",
}
ALLOWED_WITNESS_KINDS = {
    "nonuniqueness/solvability failure",
    "absence of any typed producer in the complete authority census",
    "dimensional incompatibility of every candidate producer",
    "operator/domain well-posedness failure",
}


class AssertionDeath(RuntimeError):
    def __init__(self, assert_id: str, detail: str):
        super().__init__(f"{assert_id}: {detail}")
        self.assert_id = assert_id
        self.detail = detail


CHECK_ORDER: list[str] = []


def require(condition: bool, assert_id: str, detail: str) -> None:
    if assert_id not in CHECK_ORDER:
        CHECK_ORDER.append(assert_id)
    if not condition:
        raise AssertionDeath(assert_id, detail)


def load_yaml(path: Path) -> Any:
    with path.open("rb") as handle:
        return yaml.load(handle, Loader=yaml.CSafeLoader)


def canonical_bytes(value: Any) -> bytes:
    return json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(",", ":")).encode("utf-8")


def digest(value: Any) -> str:
    return hashlib.sha256(canonical_bytes(value)).hexdigest()


def by_id(rows: list[dict[str, Any]], key: str) -> dict[str, dict[str, Any]]:
    return {row[key]: row for row in rows}


def classify_content(content: dict[str, Any]) -> set[str]:
    mapping = {
        "root": "ROOT_REFERENCE", "positive_join": "POSITIVE_DERIVATION",
        "positive_equivalence": "POSITIVE_EQUIVALENCE", "static_force": "STATIC_COMMITTED_FORCING",
        "incompatible": "INCOMPATIBILITY", "not_member": "NEGATIVE_CANDIDATE_MEMBERSHIP",
        "case_eliminate": "CASE_ELIMINATION", "complement": "COMPLEMENT_SURVIVOR_COUNT",
        "exclude": "EXCLUSION_VERDICT", "postulate": "POSTULATE_BRANCH",
        "stability": "STABILITY_DYNAMICAL_CLASS", "solve": "SOLVE_EVALUATION_RESULT",
        "symbolic_collapse": "SYMBOLIC_EQUIVALENCE_COLLAPSE", "unavailability": "UNAVAILABILITY_WITNESS",
        "challenge": "DERIVABILITY_CHALLENGE", "evidence_state": "EVIDENCE_STATE_CLASSIFICATION",
    }
    operation = content.get("op")
    if operation not in mapping:
        return {"UNCLASSIFIED"}
    values = {mapping[operation]}
    for child in content.get("args", []):
        if isinstance(child, dict):
            values.update(classify_content(child))
    return values


def inference_root_ids(content: dict[str, Any]) -> set[str]:
    values = {content["id"]} if content.get("op") == "root" and isinstance(content.get("id"), str) else set()
    for child in content.get("args", []):
        if isinstance(child, dict):
            values.update(inference_root_ids(child))
    return values


def classify_evidence(raw: dict[str, Any]) -> str:
    if not raw["applicable"]:
        return "INAPPLICABLE"
    if raw["committed_incompatibility"] and raw["ancestry_complete_and_typed"]:
        return "DECISIVE_INCOMPATIBILITY"
    if raw["datum_missing"]:
        return "MISSING"
    return "SATISFIED"


def evaluate_disposition(case: dict[str, Any]) -> str:
    evidence = case["evidence_records"]
    states = [classify_evidence(row["raw_predicate_inputs"]) for row in evidence]
    if any(row["emitted_state"] != state for row, state in zip(evidence, states)):
        return "HARNESS_FAILED(contradictory_evidence)"
    if "DECISIVE_INCOMPATIBILITY" in states:
        return "EXCLUDED"
    applicable = [state for state in states if state != "INAPPLICABLE"]
    if applicable and all(state == "SATISFIED" for state in applicable):
        return "ADMISSIBLE"
    if "MISSING" in applicable:
        return "UNRESOLVED(named datum)"
    return "HARNESS_FAILED(contradictory_evidence)"


def evaluate_promotion(inputs: dict[str, Any]) -> str:
    if inputs["failed_required_upstreams"]:
        return "NOT_RUN(exact_set)"
    if inputs["uncanonicalized_overlap"]:
        return "HARNESS_FAILED(uncanonicalized_candidate_overlap)"
    forcing = inputs["forcing_records"]
    if not forcing:
        return "NO_SELECTION_CLAIM(witness,challenge)"
    if any(not row["candidate_in_current_axis"] for row in forcing):
        return "HARNESS_FAILED(forced_class_outside_axis)"
    if len({row["canonical_class"] for row in forcing}) > 1:
        return "HARNESS_FAILED(contradictory_forcing)"
    if forcing[0]["disposition"] == "EXCLUDED":
        return "HARNESS_FAILED(contradictory_committed_derivations)"
    if forcing[0]["disposition"] == "UNRESOLVED":
        return "PROMOTION_UNRESOLVED(admissibility_unresolved_refusal)"
    if inputs["closure_outcome"] == "CLOSURE_REFUSED":
        return "PROMOTION_UNRESOLVED(closure_refusal)"
    return "SELECTED"


def evaluate_cross_level(case: dict[str, Any]) -> str:
    applicability = case["applicability"]
    if applicability == "positively_non_geometric":
        return "NOT_APPLICABLE(witness)"
    if applicability == "missing":
        return "refinement-UNRESOLVED(named datum)"
    if case["gate_integrity"] != "COMPUTATION_VALID":
        return "NOT_RUN(exact_gate_set)"
    return {
        "topologically-distinct": "defect-refined",
        "orientation-only": "not-defect-refined",
        "UNRESOLVED": "refinement-UNRESOLVED",
    }[case["gate_outcome"]]


def evaluate_topology(case: dict[str, Any]) -> str:
    pair = (case["sector_disconnection"], case["interpolation"])
    if pair in {("DISCONNECTED", "INTERPOLABLE"), ("CONNECTED", "OBSTRUCTED")}:
        return "HARNESS_FAILED(inconsistent_subresults)"
    if pair == ("DISCONNECTED", "OBSTRUCTED"):
        return "topologically-distinct"
    if case["sector_disconnection"] == "CONNECTED" or case["interpolation"] == "INTERPOLABLE":
        return "orientation-only"
    return "UNRESOLVED(named data)"


def record_schema_valid(record: dict[str, Any]) -> bool:
    integrity = record.get("integrity")
    if integrity == "COMPUTATION_VALID":
        return record.get("physics") is not None and "failure_reason" not in record and "failed_upstreams" not in record
    if integrity == "HARNESS_FAILED":
        return record.get("physics") is None and bool(record.get("failure_reason")) and "failed_upstreams" not in record
    if integrity == "NOT_RUN":
        upstreams = record.get("failed_upstreams", [])
        return record.get("physics") is None and bool(upstreams) and upstreams == sorted(set(upstreams))
    return False


def matrix_residual(fixture: dict[str, Any]) -> list[int]:
    return [
        sum(value * candidate for value, candidate in zip(matrix_row, fixture["candidate"])) - rhs
        for matrix_row, rhs in zip(fixture["matrix"], fixture["rhs"])
    ]


def walk_template_terms(node: dict[str, Any]) -> list[dict[str, str]]:
    rows = []
    if node.get("op") == "term" and node.get("coefficient") != 0:
        rows.append({"term_id": node["term_id"], "kind": node["kind"]})
    for child in node.get("args", []):
        if isinstance(child, dict):
            rows.extend(walk_template_terms(child))
    return rows


def ast_contains_op(node: dict[str, Any], operation: str) -> bool:
    return node.get("op") == operation or any(
        isinstance(child, dict) and ast_contains_op(child, operation)
        for child in node.get("args", [])
    )


def expected_template_terms(inputs: dict[str, Any]) -> list[dict[str, str]]:
    return [
        *({"term_id": f"residual:{value}", "kind": "residual"} for value in inputs["static_field_equations"]),
        *({"term_id": f"boundary:{value}", "kind": "boundary"} for value in inputs["canonical_operator_witness"]),
        *({"term_id": f"zero-mode:{value}", "kind": "zero-mode"} for value in inputs["typed_free_data_ledger"]),
        *({"term_id": f"asymptotic-matching:{value}", "kind": "asymptotic-matching"} for value in inputs["ambient_branch"]),
    ]


def closure_fixture_checks(fixture: dict[str, Any]) -> dict[str, bool]:
    contributions = fixture["committed_root_contributions"]
    expected_ids = [row["contribution_id"] for row in contributions]
    certificate = fixture.get("certificate")
    assignments = certificate.get("assignments", []) if certificate else []
    assignment_ids = [row["contribution_id"] for row in assignments]
    contribution_by_id = {row["contribution_id"]: row for row in contributions}
    total = [sum(row["vector"][index] for row in contributions) for index in range(2)]
    return {
        "certificate_present": certificate is not None,
        "census_A_exact": fixture["census_A"] == expected_ids,
        "census_B_exact": set(fixture["census_B"]) == set(expected_ids),
        "assignment_exact_once": len(assignment_ids) == len(set(assignment_ids)) and set(assignment_ids) == set(expected_ids),
        "owner_consistent": all(
            row["contribution_id"] in contribution_by_id
            and row["owner"] == contribution_by_id[row["contribution_id"]]["expected_owner"]
            and row["computed_zero"] == (contribution_by_id[row["contribution_id"]]["vector"] == [0, 0])
            for row in assignments
        ),
        "total_exact": fixture["independently_constructed_total"] == total,
    }


def standard_predicate(standard_id: str, semantic: dict[str, Any]) -> bool:
    taxonomy = semantic["evidence_taxonomy"]
    routes = semantic["route_fixture_inventory"]["route_records"]
    slots = semantic["availability_slots"]
    unresolved = [row for row in slots if row["availability_outcome"] == "UNRESOLVED"]
    vocabulary = semantic["vocabulary_freeze"]
    predicates = {
        "S-1": lambda: all({"source_id", "category", "root_type"} <= set(row) for row in taxonomy["source_census"]["records"]),
        "S-2": lambda: vocabulary["evidence_classification"].startswith("transparent") and all(row["membership_predicate"].startswith("positive_") for row in semantic["candidate_inventory"]["candidate_records"] if row["candidate_id"] != "OTHER"),
        "S-3": lambda: len({tuple(tuple(value for value in matrix_row) for matrix_row in row["positive_fixture"]["matrix"]) for row in routes}) >= 6 and all(row["positive_fixture"]["nondegenerate_norm_squared"] > 0 for row in routes),
        "S-4": lambda: semantic["dimensional_firewall"]["base_pass"] and semantic["dimensional_firewall"]["ablation_fired"],
        "S-5": lambda: all(row["tooth_id"] == f"STANDARD_TOOTH:{row['standard_id']}" for row in semantic["standard_bindings"]),
        "S-6": lambda: len(semantic["template_contract"]["semantic_schema"]) == 6 and len(semantic["guard_fixtures"]["template"]["expected_term_census"]) == 4,
        "P-1": lambda: all(row.get("witness_id") == f"witness:{row['slot_id']}" for row in unresolved),
        "P-2": lambda: semantic["candidate_inventory"]["mixture_generator_A"] == semantic["candidate_inventory"]["mixture_generator_B"] and all(row["generator_A"] == row["generator_B"] for row in semantic["obligation_censuses"].values()),
        "P-3": lambda: "nonholonomic" in semantic["obligation_censuses"]["E4"]["native_root_class"] and "rayleigh" in semantic["obligation_censuses"]["E5"]["native_root_class"],
        "P-4": lambda: len(vocabulary["cross_level_decision_table"]) == 15 and len({(row["applicability"], row["gate_integrity"], row["gate_outcome"]) for row in vocabulary["cross_level_decision_table"]}) == 15,
        "P-5": lambda: set(semantic["guard_fixtures"]["closure"]["census_A"]) == set(semantic["guard_fixtures"]["closure"]["census_B"]),
        "P-6": lambda: "radiation/static-zero" in vocabulary["closure_channels"],
        "P-7": lambda: all("inherited_ratified_witness_reference" in row["witness"] for row in unresolved if row["slot_id"].startswith("template_free_data:")),
        "P-8": lambda: all(matrix_residual(row["positive_fixture"]) == row["positive_fixture"]["residual"] == [0] * len(row["positive_fixture"]["rhs"]) and any(value != 0 for value in matrix_residual(row["malformed_fixture"])) for row in routes),
        "P-9": lambda: set(vocabulary["topology_gate_enum"]) == {"topologically-distinct", "orientation-only", "UNRESOLVED(named datum)"},
        "P-10": lambda: taxonomy["guard_fixture_dags"]["pair_configuration"]["consumer"] == "topology_pair_annihilation_subquestion",
    }
    return bool(predicates[standard_id]())


def mutation_catalog(semantic: dict[str, Any]) -> list[dict[str, str]]:
    entries: list[tuple[str, str]] = [
        ("TOOTH_ENGINE_SCHEMA", "ASSERT_ENGINE_SCHEMA"),
        ("TOOTH_ENGINE_SEMANTIC", "ASSERT_ENGINE_SEMANTIC_AGREEMENT"),
        ("TOOTH_ENGINE_ROUTE_IDS", "ASSERT_ENGINE_ROUTE_IDS"),
        ("TOOTH_ENGINE_LOCAL_EXECUTION", "ASSERT_ENGINE_LOCAL_SYMBOLS"),
        ("TOOTH_COMPARATOR_OUTPUT", "ASSERT_COMPARATOR_OUTPUT"),
        ("TOOTH_SCOPE_INVARIANTS", "ASSERT_SCOPE_INVARIANTS"),
        ("TOOTH_MIXTURE_DUAL_GENERATORS", "ASSERT_MIXTURE_DUAL_GENERATORS"),
        ("TOOTH_CANDIDATE_AXIS_DISJOINT", "ASSERT_CANDIDATE_AXIS_DISJOINT"),
        ("TOOTH_UNCANONICALIZED_OVERLAP", "ASSERT_CANDIDATE_CANONICALIZATION"),
        ("TOOTH_CANDIDATE_UNIVERSE_DIGEST", "ASSERT_CANDIDATE_UNIVERSE_DIGEST"),
        ("TOOTH_CATCHALL_RULE", "ASSERT_CATCHALL_PRESENCE"),
        ("TOOTH_CONCRETE_OTHER_ALIAS", "ASSERT_CONCRETE_OTHER_CANONICAL"),
        ("TOOTH_MEMBERSHIP_PNF", "ASSERT_MEMBERSHIP_PREDICATES"),
        ("TOOTH_AMENDMENT_REKEY", "ASSERT_AMENDMENT_REKEY"),
        ("TOOTH_SOURCE_CENSUS", "ASSERT_SOURCE_CENSUS_COMPLETE"),
        ("TOOTH_TAXONOMY_CLOSED", "ASSERT_TAXONOMY_CLOSED"),
        ("CONCEALED_NEGATIVE_IN_POSITIVE_NODE", "ASSERT_CONTENT_CLASSIFICATION"),
        ("STABILITY_RELABELLED_STATIC_FORCING", "ASSERT_CONTENT_CLASSIFICATION"),
        ("TOOTH_STAGE0_ARTIFACT_COVERAGE", "ASSERT_STAGE0_ARTIFACT_COVERAGE"),
        ("STABILITY_IN_CANDIDATE_FORMATION", "ASSERT_STAGE0_ARTIFACT_ANCESTRY"),
        ("ELIMINATION_IN_CANONICALIZATION", "ASSERT_STAGE0_MEMBERSHIP_ANCESTRY"),
        ("SELECTION_FROM_NEGATIVE_LEMMA_NO_EXCLUDED_RECORD", "ASSERT_PROMOTION_ANCESTRY_ALLOWED"),
        ("STABILITY_IN_NONPROMOTION_VERDICT", "ASSERT_PROGRAM_WIDE_STABILITY_BAN"),
        ("TWO_BODY_PAIR_CONFIGURATION_MISCONSUMED", "ASSERT_ONE_BODY_FIREWALL"),
        ("HIDDEN_SOLVE_IN_TEMPLATE_ANCESTRY", "ASSERT_BVP_BOUNDARY"),
        ("POSTULATE_IN_PROMOTION_ANCESTRY", "ASSERT_POSTULATE_PROMOTION_BAN"),
        ("MAXWELL_NATIVE_ROOT_IMPORT", "ASSERT_NATIVE_ANCESTRY"),
        ("TOOTH_OBLIGATION_DUAL_GENERATORS", "ASSERT_OBLIGATION_DUAL_GENERATORS"),
        ("TOOTH_CENSUS_STRATUM_INPUT", "ASSERT_CENSUS_STRATUM_FREE"),
        ("TOOTH_DECORATIVE_ABLATION", "ASSERT_NONDEGENERACY_ANTI_ECHO"),
        ("TOOTH_OPEN_DEPENDENCY_JOIN", "ASSERT_OPEN_DEPENDENCY_JOIN"),
        ("TOOTH_OPEN_DEPENDENCY_INDEPENDENCE", "ASSERT_OPEN_DEPENDENCY_INDEPENDENCE"),
        ("TOOTH_ACTIVE_STRATA", "ASSERT_ACTIVE_STRATA"),
        ("TOOTH_RAW_GRID_CARDINALITY", "ASSERT_RAW_GRID_CARDINALITY"),
        ("TOOTH_GRID_COVERAGE", "ASSERT_GRID_COVERAGE"),
        ("TOOTH_COLLAPSE_TIMING", "ASSERT_COLLAPSE_TIMING"),
        ("TOOTH_COLLAPSED_CARDINALITY", "ASSERT_COLLAPSED_CARDINALITY"),
        ("TOOTH_PROMOTION_CONTEXT_COVERAGE", "ASSERT_PROMOTION_CONTEXT_COVERAGE"),
        ("TOOTH_ROUTE_PREPRODUCTION", "ASSERT_ROUTE_INVENTORY_PREPRODUCTION"),
        ("TOOTH_ROUTE_FIXTURE_COVERAGE", "ASSERT_ROUTE_FIXTURE_COVERAGE"),
        ("TOOTH_FIXTURE_POSITIVE", "ASSERT_FIXTURE_POSITIVE"),
        ("TOOTH_FIXTURE_MALFORMED", "ASSERT_FIXTURE_MALFORMED"),
        ("TOOTH_E4_E5_NATIVE_ROOT", "ASSERT_E4_E5_NATIVE_ROOT"),
        ("TOOTH_INTEGRITY_RULE", "ASSERT_VOCAB_INTEGRITY_RULE"),
        ("TOOTH_HARNESS_FAILED_PHYSICS_VALUE", "ASSERT_VOCAB_INTEGRITY_RULE"),
        ("TOOTH_NOT_RUN_PHYSICS_VALUE", "ASSERT_VOCAB_INTEGRITY_RULE"),
        ("TOOTH_UPSTREAM_PROPAGATION_SKIP", "ASSERT_VOCAB_INTEGRITY_RULE"),
        ("TOOTH_TYPED_FAILURE_REASON", "ASSERT_TYPED_FAILURE_REASONS"),
        ("TOOTH_EVIDENCE_ENUM", "ASSERT_EVIDENCE_STATE_FREEZE"),
        ("TOOTH_DISPOSITION_TABLE", "ASSERT_DISPOSITION_TABLE_TOTAL"),
        ("DISPOSITION_PRECEDENCE_DECISIVE_PLUS_MISSING_TO_UNRESOLVED", "ASSERT_DISPOSITION_PRECEDENCE"),
        ("EVIDENCE_DECISIVE_RELABELLED_OR_DROPPED", "ASSERT_EVIDENCE_STATE_CLASSIFICATION"),
        ("TOOTH_PROMOTION_TABLE_TOTAL", "ASSERT_PROMOTION_TABLE_TOTAL"),
        ("PROMOTION_REQUIRED_POSITIVE_TO_NO_SELECTION", "ASSERT_PROMOTION_DECISION_TABLE"),
        ("PROMOTION_CLOSURE_REFUSAL_TO_NO_SELECTION", "ASSERT_PROMOTION_DECISION_TABLE"),
        ("PROMOTION_FORCED_UNRESOLVED_WRONG_REFUSAL", "ASSERT_PROMOTION_DECISION_TABLE"),
        ("PROMOTION_FORCED_EXCLUDED_TO_PHYSICS", "ASSERT_PROMOTION_DECISION_TABLE"),
        ("PROMOTION_MULTI_FORCING_TO_SINGLE_CLASS", "ASSERT_PROMOTION_DECISION_TABLE"),
        ("PROMOTION_OUTSIDE_AXIS_TO_PHYSICS", "ASSERT_PROMOTION_DECISION_TABLE"),
        ("PROMOTION_ALIAS_TO_SELECTED", "ASSERT_PROMOTION_DECISION_TABLE"),
        ("PROMOTION_NOT_RUN_TO_PHYSICS", "ASSERT_PROMOTION_NOT_RUN_PROPAGATION"),
        ("TOOTH_ENSEMBLE_ENUMS", "ASSERT_ENSEMBLE_ENUMS"),
        ("TOOTH_CROSS_LEVEL_TABLE", "ASSERT_CROSS_LEVEL_TABLE"),
        ("MIXED_GEOMETRIC_TO_NOT_APPLICABLE", "ASSERT_MIXED_GEOMETRIC_FIXTURE"),
        ("MIXED_NONGEOMETRIC_TO_REFINEMENT", "ASSERT_MIXED_NONGEOMETRIC_FIXTURE"),
        ("ORIENTATION_ONLY_TO_REFINEMENT_UNRESOLVED", "ASSERT_ORIENTATION_NEGATIVE_PRESERVED"),
        ("GATE_INTEGRITY_TO_LEVEL2_PHYSICS", "ASSERT_GATE_INTEGRITY_PROPAGATION"),
        ("TOOTH_TOPOLOGY_TABLE_TOTAL", "ASSERT_TOPOLOGY_TABLE_TOTAL"),
        ("TOPOLOGY_INCONSISTENCY_TO_PHYSICS", "ASSERT_TOPOLOGY_INCONSISTENCY"),
        ("PAIR_ANNIHILATION_POLARITY_IN_AGGREGATE", "ASSERT_PAIR_ORTHOGONAL"),
        ("TOPOLOGY_E4_E5_BULK_ACTION_ONLY", "ASSERT_TOPOLOGY_NATIVE_ROOT"),
        ("TOOTH_ENSEMBLE_EXCHANGE", "ASSERT_ENSEMBLE_EXCHANGE_CONTRACT"),
        ("CLOSURE_MISSING_CERTIFICATE", "ASSERT_CLOSURE_GATE"),
        ("CLOSURE_DUPLICATE_OWNER", "ASSERT_CLOSURE_GATE"),
        ("CLOSURE_OMITTED_CONTRIBUTION", "ASSERT_CLOSURE_GATE"),
        ("CLOSURE_CENSUS_SIDE_OMISSION", "ASSERT_CLOSURE_GATE"),
        ("CLOSURE_INTEGRITY_LAUNDERING", "ASSERT_CLOSURE_GATE"),
        ("CLOSURE_INTEGRITY_TO_PHYSICS_REFUSAL", "ASSERT_CLOSURE_PROPAGATION"),
        ("TOOTH_TEMPLATE_SCHEMA", "ASSERT_TEMPLATE_SCHEMA"),
        ("TEMPLATE_REMOVE_RESIDUAL_TERM", "ASSERT_TEMPLATE_TERM_INCIDENCE"),
        ("TEMPLATE_ZERO_BOUNDARY_TERM", "ASSERT_TEMPLATE_TERM_INCIDENCE"),
        ("TEMPLATE_REMOVE_ZERO_MODE_TERM", "ASSERT_TEMPLATE_TERM_INCIDENCE"),
        ("TEMPLATE_ZERO_ASYMPTOTIC_TERM", "ASSERT_TEMPLATE_TERM_INCIDENCE"),
        ("TEMPLATE_HIDDEN_SOLVED_PAYLOAD", "ASSERT_TEMPLATE_NO_SOLVE"),
        ("RETURN_CLOSURE_CONSUMED_BY_U2_VERDICT", "ASSERT_RETURN_CLOSURE_NONOWNERSHIP"),
        ("TOOTH_AVAILABILITY_FLOOR", "ASSERT_AVAILABILITY_FLOOR"),
        ("TOOTH_WITNESS_SCHEMA", "ASSERT_WITNESS_SCHEMA"),
        ("TOOTH_CHALLENGE_SCHEMA", "ASSERT_CHALLENGE_SCHEMA"),
        ("TOOTH_WITNESS_CHALLENGE_LOCKSTEP", "ASSERT_WITNESS_CHALLENGE_LOCKSTEP"),
        ("TOOTH_WITNESS_RESTORE", "ASSERT_WITNESS_RESTORE"),
        ("OVERWRITE_DERIVED_WITH_UNRESOLVED", "ASSERT_DERIVED_OVERWRITE_CHALLENGE"),
        ("FIRST_TIME_CONSTRUCTION_HIDDEN_AS_UNRESOLVED", "ASSERT_FIRST_TIME_CONSTRUCTION"),
        ("TOOTH_DERIVED_OBJECT_CONTENT", "ASSERT_DERIVED_OBJECT_CONTENT"),
        ("TOOTH_STANDARD_COVERAGE", "ASSERT_STANDARD_COVERAGE"),
        ("TOOTH_DIMENSIONAL_FIREWALL", "ASSERT_DIMENSIONAL_FIREWALL"),
        ("TOOTH_BRANCH_ID_COVERAGE", "ASSERT_BRANCH_ID_COVERAGE"),
        ("TOOTH_STAGE0_NONZERO_INTEGRITY", "ASSERT_ZERO_STAGE0_INTEGRITY_FAILURES"),
        ("ALL_INTEGRITY_FAILURE_GRID_BANKED", "ASSERT_ALL_INTEGRITY_REJECTION"),
        ("TOOTH_MUTATION_CATALOG_COVERAGE", "ASSERT_MUTATION_CATALOG_COVERAGE"),
    ]
    unresolved = [row for row in semantic["availability_slots"] if row["availability_outcome"] == "UNRESOLVED"]
    entries.extend((f"ENTRY_WITNESS_DROP:{row['slot_id']}", f"ASSERT_ENTRY_WITNESS_RESOLVES:{row['slot_id']}") for row in unresolved)
    entries.extend((f"DERIVABILITY_CANARY:{row['slot_id']}", f"ASSERT_STAGE0_UNRESOLVED_REFUTED:{row['slot_id']}") for row in unresolved)
    entries.extend((f"STANDARD_TOOTH:{row['standard_id']}", row["reachable_check_id"]) for row in semantic["standard_bindings"])
    return [{"mutation_id": mutation_id, "expected_assert_id": assert_id} for mutation_id, assert_id in entries]


def find_artifact(semantic: dict[str, Any], artifact_id: str) -> dict[str, Any]:
    return next(row for row in semantic["evidence_taxonomy"]["physics_bearing_stage0_artifacts"] if row["artifact_id"] == artifact_id)


def apply_mutation(
    mutation: str | None, sympy: dict[str, Any], wolfram: dict[str, Any]
) -> None:
    if not mutation:
        return
    if mutation == "TOOTH_COMPARATOR_OUTPUT":
        return
    s = sympy["semantic_view"]
    w = wolfram["semantic_view"]
    both = (s, w)
    if mutation == "TOOTH_ENGINE_SCHEMA":
        sympy["schema_version"] = wolfram["schema_version"] = "MUTATED"
    elif mutation == "TOOTH_ENGINE_SEMANTIC":
        s["scope"]["one_body_only"] = False
    elif mutation == "TOOTH_ENGINE_ROUTE_IDS":
        sympy["engine_local_route_registry"][0]["semantic_route_id"] = "mutated_route"
    elif mutation == "TOOTH_ENGINE_LOCAL_EXECUTION":
        for engine in (sympy, wolfram):
            engine["engine_local_route_registry"][0]["executed"] = False
    elif mutation == "TOOTH_SCOPE_INVARIANTS":
        for value in both: value["scope"]["one_body_only"] = False
    elif mutation == "TOOTH_MIXTURE_DUAL_GENERATORS":
        for value in both: value["candidate_inventory"]["mixture_generator_B"].pop()
    elif mutation == "TOOTH_CANDIDATE_AXIS_DISJOINT":
        for value in both: value["candidate_inventory"]["candidate_axis"].append("E1")
    elif mutation == "TOOTH_UNCANONICALIZED_OVERLAP":
        for value in both: value["candidate_inventory"]["uncanonicalized_overlap_count"] = 1
    elif mutation == "TOOTH_CANDIDATE_UNIVERSE_DIGEST":
        for value in both: value["candidate_inventory"]["candidate_universe_digest"] = "0" * 64
    elif mutation == "TOOTH_CATCHALL_RULE":
        for value in both: value["candidate_inventory"]["basis_closure"]["catch_all_present"] = False
    elif mutation == "TOOTH_CONCRETE_OTHER_ALIAS":
        for value in both: value["candidate_inventory"]["concrete_other_candidates"] = [{"candidate_id": "OTHER(E1_alias)", "signature": ["normal_velocity_lock", "tangential_velocity_lock"]}]
    elif mutation == "TOOTH_MEMBERSHIP_PNF":
        for value in both:
            inventory = value["candidate_inventory"]
            inventory["candidate_records"][0]["membership_predicate"] = "negative_membership_by_elimination"
            inventory["candidate_universe_digest"] = digest(inventory["candidate_records"])
    elif mutation == "TOOTH_AMENDMENT_REKEY":
        for value in both: value["candidate_inventory"]["amendment_rekey_duty"] = ""
    elif mutation == "TOOTH_SOURCE_CENSUS":
        for value in both: value["evidence_taxonomy"]["source_census"]["records"].pop()
    elif mutation == "TOOTH_TAXONOMY_CLOSED":
        for value in both: value["evidence_taxonomy"]["constructor_grammar"].append("GENERIC")
    elif mutation in {"CONCEALED_NEGATIVE_IN_POSITIVE_NODE", "STABILITY_RELABELLED_STATIC_FORCING"}:
        operation = "not_member" if mutation.startswith("CONCEALED") else "stability"
        for value in both:
            artifact = find_artifact(value, "membership_equivalence_predicates")
            artifact["normalized_inference_content"] = {"op": "positive_equivalence", "args": [{"op": operation}]}
    elif mutation == "TOOTH_STAGE0_ARTIFACT_COVERAGE":
        for value in both:
            artifact = find_artifact(value, "obligation_censuses")
            artifact["reachable_root_ids"].pop()
            artifact["normalized_inference_content"]["args"].pop()
    elif mutation == "STABILITY_IN_CANDIDATE_FORMATION":
        for value in both:
            artifact = find_artifact(value, "mixture_inventory")
            artifact["normalized_inference_content"] = {"op": "stability", "args": [{"op": "root", "id": root_id} for root_id in artifact["reachable_root_ids"]]}
            artifact["computed_constructors"] = ["ROOT_REFERENCE", "STABILITY_DYNAMICAL_CLASS"]
    elif mutation == "ELIMINATION_IN_CANONICALIZATION":
        for value in both:
            artifact = find_artifact(value, "canonicalization_alias_proofs")
            artifact["normalized_inference_content"] = {"op": "case_eliminate", "args": [{"op": "root", "id": root_id} for root_id in artifact["reachable_root_ids"]]}
            artifact["computed_constructors"] = ["CASE_ELIMINATION", "ROOT_REFERENCE"]
    elif mutation == "SELECTION_FROM_NEGATIVE_LEMMA_NO_EXCLUDED_RECORD":
        for value in both: value["evidence_taxonomy"]["guard_fixture_dags"]["promotion_positive"]["constructors"].append("NEGATIVE_CANDIDATE_MEMBERSHIP")
    elif mutation == "STABILITY_IN_NONPROMOTION_VERDICT":
        for value in both: value["evidence_taxonomy"]["guard_fixture_dags"]["candidate_exclusion"]["constructors"].append("STABILITY_DYNAMICAL_CLASS")
    elif mutation == "TWO_BODY_PAIR_CONFIGURATION_MISCONSUMED":
        for value in both: value["evidence_taxonomy"]["guard_fixture_dags"]["pair_configuration"]["consumer"] = "two_body_force_consumer"
    elif mutation == "HIDDEN_SOLVE_IN_TEMPLATE_ANCESTRY":
        for value in both: value["evidence_taxonomy"]["guard_fixture_dags"]["posed_template"]["constructors"].append("SOLVE_EVALUATION_RESULT")
    elif mutation == "POSTULATE_IN_PROMOTION_ANCESTRY":
        for value in both: value["evidence_taxonomy"]["guard_fixture_dags"]["postulate_deferral_metadata"]["promotion_reachable"] = True
    elif mutation == "MAXWELL_NATIVE_ROOT_IMPORT":
        for value in both: value["evidence_taxonomy"]["guard_fixture_dags"]["candidate_exclusion"]["root_types"].append("MAXWELL_IMPORTED_ROOT")
    elif mutation == "TOOTH_OBLIGATION_DUAL_GENERATORS":
        for value in both: value["obligation_censuses"]["E1"]["generator_B"].pop()
    elif mutation == "TOOTH_CENSUS_STRATUM_INPUT":
        for value in both: value["obligation_censuses"]["E1"]["stratum_token_is_generation_input"] = True
    elif mutation == "TOOTH_DECORATIVE_ABLATION":
        for value in both: value["obligation_censuses"]["E1"]["semantic_ablation_criterion"] = "digest_changed"
    elif mutation == "TOOTH_OPEN_DEPENDENCY_JOIN":
        for value in both: value["open_dependency_relation"]["dependency_rows"][0]["generator_B"].pop()
    elif mutation == "TOOTH_OPEN_DEPENDENCY_INDEPENDENCE":
        for value in both:
            value["open_dependency_relation"]["generator_B_algorithm"] = value["open_dependency_relation"]["generator_A_algorithm"]
    elif mutation == "TOOTH_ACTIVE_STRATA":
        for value in both: value["open_dependency_relation"]["active_strata"].pop()
    elif mutation == "TOOTH_RAW_GRID_CARDINALITY":
        for value in both: value["grid_inventory"]["raw_ragged_cardinality"] -= 1
    elif mutation == "TOOTH_GRID_COVERAGE":
        for value in both: value["grid_inventory"]["grid_cells"].pop()
    elif mutation == "TOOTH_COLLAPSE_TIMING":
        for value in both: value["grid_inventory"]["collapse_proofs"][0]["stage0_objects_only"] = False
    elif mutation == "TOOTH_COLLAPSED_CARDINALITY":
        for value in both: value["grid_inventory"]["collapsed_cardinality"] -= 1
    elif mutation == "TOOTH_PROMOTION_CONTEXT_COVERAGE":
        for value in both: value["grid_inventory"]["promotion_contexts"][0]["candidate_cell_mappings"].pop()
    elif mutation == "TOOTH_ROUTE_PREPRODUCTION":
        for value in both: value["route_fixture_inventory"]["generated_preproduction"] = False
    elif mutation == "TOOTH_ROUTE_FIXTURE_COVERAGE":
        for value in both: value["route_fixture_inventory"]["route_records"].pop()
    elif mutation == "TOOTH_FIXTURE_POSITIVE":
        for value in both: value["route_fixture_inventory"]["route_records"][0]["positive_fixture"]["residual"] = [1, 0]
    elif mutation == "TOOTH_FIXTURE_MALFORMED":
        for value in both: value["route_fixture_inventory"]["route_records"][0]["malformed_fixture"]["residual"] = [0, 0]
    elif mutation == "TOOTH_E4_E5_NATIVE_ROOT":
        for value in both: value["obligation_censuses"]["E4"]["native_root_class"] = "variational_bulk"
    elif mutation == "TOOTH_INTEGRITY_RULE":
        for value in both: value["vocabulary_freeze"]["uniform_record_rule"]["record_types"].pop()
    elif mutation == "TOOTH_HARNESS_FAILED_PHYSICS_VALUE":
        for value in both:
            next(row for row in value["guard_fixtures"]["record_schema_cases"] if row["integrity"] == "HARNESS_FAILED")["physics"] = "ADMISSIBLE"
    elif mutation == "TOOTH_NOT_RUN_PHYSICS_VALUE":
        for value in both:
            next(row for row in value["guard_fixtures"]["record_schema_cases"] if row["integrity"] == "NOT_RUN")["physics"] = "UNRESOLVED"
    elif mutation == "TOOTH_UPSTREAM_PROPAGATION_SKIP":
        for value in both:
            value["guard_fixtures"]["upstream_propagation"]["emitted_dependent"].update({"integrity": "COMPUTATION_VALID", "physics": "ADMISSIBLE", "failed_upstreams": []})
    elif mutation == "TOOTH_TYPED_FAILURE_REASON":
        for value in both: value["vocabulary_freeze"]["typed_failure_reasons"]["free_text_reason"] = "builder prose"
    elif mutation == "TOOTH_EVIDENCE_ENUM":
        for value in both: value["vocabulary_freeze"]["obligation_evidence_state_enum"].pop()
    elif mutation == "TOOTH_DISPOSITION_TABLE":
        for value in both: value["vocabulary_freeze"]["disposition_precedence_table"].pop()
    elif mutation == "DISPOSITION_PRECEDENCE_DECISIVE_PLUS_MISSING_TO_UNRESOLVED":
        for value in both:
            next(row for row in value["guard_fixtures"]["disposition_cases"] if row["case_id"] == "decisive_plus_missing")["emitted_landing"] = "UNRESOLVED(named datum)"
    elif mutation == "EVIDENCE_DECISIVE_RELABELLED_OR_DROPPED":
        for value in both:
            case = next(row for row in value["guard_fixtures"]["disposition_cases"] if row["case_id"] == "decisive_plus_missing")
            case["evidence_records"][0]["emitted_state"] = "MISSING"
    elif mutation == "TOOTH_PROMOTION_TABLE_TOTAL":
        for value in both: value["vocabulary_freeze"]["promotion_decision_table"].pop()
    elif mutation.startswith("PROMOTION_") and mutation != "PROMOTION_NOT_RUN_TO_PHYSICS":
        changes = {
            "PROMOTION_REQUIRED_POSITIVE_TO_NO_SELECTION": ("required_positive", "NO_SELECTION_CLAIM"),
            "PROMOTION_CLOSURE_REFUSAL_TO_NO_SELECTION": ("closure_refusal", "NO_SELECTION_CLAIM"),
            "PROMOTION_FORCED_UNRESOLVED_WRONG_REFUSAL": ("forced_unresolved", "PROMOTION_UNRESOLVED(closure_refusal)"),
            "PROMOTION_FORCED_EXCLUDED_TO_PHYSICS": ("forced_excluded", "NO_SELECTION_CLAIM"),
            "PROMOTION_MULTI_FORCING_TO_SINGLE_CLASS": ("multi_forcing", "SELECTED"),
            "PROMOTION_OUTSIDE_AXIS_TO_PHYSICS": ("outside_axis", "SELECTED"),
            "PROMOTION_ALIAS_TO_SELECTED": ("alias", "SELECTED"),
        }
        key, landing = changes[mutation]
        for value in both:
            next(row for row in value["guard_fixtures"]["promotion_cases"] if row["case_id"] == key)["emitted_landing"] = landing
    elif mutation == "PROMOTION_NOT_RUN_TO_PHYSICS":
        for value in both:
            next(row for row in value["guard_fixtures"]["promotion_cases"] if row["case_id"] == "failed_upstream")["emitted_landing"] = "NO_SELECTION_CLAIM(witness,challenge)"
    elif mutation == "TOOTH_ENSEMBLE_ENUMS":
        for value in both: value["vocabulary_freeze"]["ensemble_level_1_enum"].pop()
    elif mutation == "TOOTH_CROSS_LEVEL_TABLE":
        for value in both: value["vocabulary_freeze"]["cross_level_decision_table"].pop()
    elif mutation in {"MIXED_GEOMETRIC_TO_NOT_APPLICABLE", "MIXED_NONGEOMETRIC_TO_REFINEMENT", "ORIENTATION_ONLY_TO_REFINEMENT_UNRESOLVED", "GATE_INTEGRITY_TO_LEVEL2_PHYSICS"}:
        targets = {
            "MIXED_GEOMETRIC_TO_NOT_APPLICABLE": ("geometric", "COMPUTATION_VALID", "topologically-distinct", "NOT_APPLICABLE(witness)"),
            "MIXED_NONGEOMETRIC_TO_REFINEMENT": ("positively_non_geometric", "COMPUTATION_VALID", "orientation-only", "defect-refined"),
            "ORIENTATION_ONLY_TO_REFINEMENT_UNRESOLVED": ("geometric", "COMPUTATION_VALID", "orientation-only", "refinement-UNRESOLVED"),
            "GATE_INTEGRITY_TO_LEVEL2_PHYSICS": ("geometric", "HARNESS_FAILED", None, "not-defect-refined"),
        }
        applicability, integrity, outcome, landing = targets[mutation]
        for value in both:
            case = next(row for row in value["guard_fixtures"]["cross_level_cases"] if (row["applicability"], row["gate_integrity"], row["gate_outcome"]) == (applicability, integrity, outcome))
            case["landing"] = landing
    elif mutation == "TOOTH_TOPOLOGY_TABLE_TOTAL":
        for value in both: value["vocabulary_freeze"]["topology_aggregate_table"].pop()
    elif mutation == "TOPOLOGY_INCONSISTENCY_TO_PHYSICS":
        for value in both:
            next(row for row in value["guard_fixtures"]["topology_cases"] if row["sector_disconnection"] == "DISCONNECTED" and row["interpolation"] == "INTERPOLABLE")["emitted_landing"] = "orientation-only"
    elif mutation == "PAIR_ANNIHILATION_POLARITY_IN_AGGREGATE":
        for value in both: value["guard_fixtures"]["topology_cases"][0]["pair_used_in_aggregate"] = True
    elif mutation == "TOPOLOGY_E4_E5_BULK_ACTION_ONLY":
        for value in both:
            route = next(
                row for row in value["route_fixture_inventory"]["route_records"]
                if "nonholonomic" in row["native_root_class"] or "rayleigh" in row["native_root_class"]
            )
            route["positive_fixture"]["native_structure_exercised"].update({"constraint_multiplier": False, "dissipation_bookkeeping": False})
    elif mutation == "TOOTH_ENSEMBLE_EXCHANGE":
        for value in both:
            value["obligation_censuses"]["E1"]["generator_A"].remove("ensemble_exchange_discharge")
            value["obligation_censuses"]["E1"]["generator_B"].remove("ensemble_exchange_discharge")
    elif mutation.startswith("CLOSURE_"):
        for value in both:
            fixture = value["guard_fixtures"]["closure"]
            if mutation == "CLOSURE_MISSING_CERTIFICATE": fixture.pop("certificate")
            elif mutation == "CLOSURE_DUPLICATE_OWNER": fixture["certificate"]["assignments"].append(copy.deepcopy(fixture["certificate"]["assignments"][0]))
            elif mutation == "CLOSURE_OMITTED_CONTRIBUTION": fixture["certificate"]["assignments"].pop()
            elif mutation == "CLOSURE_CENSUS_SIDE_OMISSION": fixture["census_A"].pop()
            elif mutation == "CLOSURE_INTEGRITY_LAUNDERING":
                fixture["dependent_record"] = {"integrity": "COMPUTATION_VALID", "physics": "PROMOTION_UNRESOLVED(closure_refusal)"}
            elif mutation == "CLOSURE_INTEGRITY_TO_PHYSICS_REFUSAL":
                fixture["closure_adjudication"] = {"integrity": "HARNESS_FAILED", "physics": None, "failure_reason": "schema_violation"}
                fixture["dependent_record"] = {"integrity": "COMPUTATION_VALID", "physics": "PROMOTION_UNRESOLVED(closure_refusal)"}
    elif mutation == "TOOTH_TEMPLATE_SCHEMA":
        for value in both: value["template_contract"]["semantic_schema"].pop()
    elif mutation.startswith("TEMPLATE_") and mutation != "TEMPLATE_HIDDEN_SOLVED_PAYLOAD":
        key = {
            "TEMPLATE_REMOVE_RESIDUAL_TERM": "residual", "TEMPLATE_ZERO_BOUNDARY_TERM": "boundary",
            "TEMPLATE_REMOVE_ZERO_MODE_TERM": "zero-mode", "TEMPLATE_ZERO_ASYMPTOTIC_TERM": "asymptotic-matching",
        }[mutation]
        for value in both:
            term = next(row for row in value["guard_fixtures"]["template"]["template_record"]["symbolic_ast"]["args"] if row["kind"] == key)
            if mutation.startswith("TEMPLATE_REMOVE_"): value["guard_fixtures"]["template"]["template_record"]["symbolic_ast"]["args"].remove(term)
            else: term["coefficient"] = 0
    elif mutation == "TEMPLATE_HIDDEN_SOLVED_PAYLOAD":
        for value in both:
            value["guard_fixtures"]["template"]["template_record"]["symbolic_ast"]["args"].append({"op": "solve", "args": []})
    elif mutation == "RETURN_CLOSURE_CONSUMED_BY_U2_VERDICT":
        for value in both: value["return_closure_ownership"]["U2_owned"] = True
    elif mutation == "TOOTH_AVAILABILITY_FLOOR":
        for value in both: value["availability_slots"].pop()
    elif mutation == "TOOTH_WITNESS_SCHEMA":
        for value in both:
            row = next(item for item in value["availability_slots"] if item["availability_outcome"] == "UNRESOLVED")
            row["witness"]["kind"] = "not found"
    elif mutation == "TOOTH_CHALLENGE_SCHEMA":
        for value in both:
            row = next(item for item in value["availability_slots"] if item["availability_outcome"] == "UNRESOLVED")
            row["challenge"]["empty_output"] = True
    elif mutation == "TOOTH_WITNESS_CHALLENGE_LOCKSTEP":
        for value in both:
            row = next(item for item in value["availability_slots"] if item["availability_outcome"] == "UNRESOLVED")
            row["challenge"]["kind"] = "operator/domain well-posedness failure"
    elif mutation == "TOOTH_WITNESS_RESTORE":
        for value in both:
            next(row for row in value["availability_slots"] if row["availability_outcome"] == "UNRESOLVED")["witness"]["counterfactual_restore_mutation"]["restored_status"] = "PASS_COMPUTED"
    elif mutation == "OVERWRITE_DERIVED_WITH_UNRESOLVED":
        for value in both:
            next(row for row in value["availability_slots"] if row["slot_id"].startswith("candidate_definition:") and row["availability_outcome"] == "DERIVED")["availability_outcome"] = "UNRESOLVED"
    elif mutation == "FIRST_TIME_CONSTRUCTION_HIDDEN_AS_UNRESOLVED":
        for value in both:
            next(row for row in value["availability_slots"] if row["slot_id"].startswith("mixture_law:"))["availability_outcome"] = "UNRESOLVED"
    elif mutation == "TOOTH_DERIVED_OBJECT_CONTENT":
        for value in both:
            next(row for row in value["availability_slots"] if row["availability_outcome"] == "DERIVED")["derived_object"].pop("canonical_signature", None)
    elif mutation == "TOOTH_STANDARD_COVERAGE":
        for value in both: value["standard_bindings"].pop()
    elif mutation.startswith("STANDARD_TOOTH:"):
        standard = mutation.split(":", 1)[1]
        for value in both:
            next(row for row in value["standard_bindings"] if row["standard_id"] == standard)["acceptance_predicate_id"] = "mutated:unbound"
    elif mutation == "TOOTH_DIMENSIONAL_FIREWALL":
        for value in both: value["dimensional_firewall"]["ablation_fired"] = False
    elif mutation == "TOOTH_BRANCH_ID_COVERAGE":
        for value in both: value["grid_inventory"]["grid_cells"][0]["stable_branch_id"] = ""
    elif mutation == "TOOTH_STAGE0_NONZERO_INTEGRITY":
        for value in both: value["availability_summary"]["integrity_failures"] = 1
    elif mutation == "ALL_INTEGRITY_FAILURE_GRID_BANKED":
        for value in both: value["guard_fixtures"]["banking_fixture"]["emitted_approval"] = True
    elif mutation == "TOOTH_MUTATION_CATALOG_COVERAGE":
        for value in both: value["production_contract"]["mutation_catalog_coverage_required"] = False
    elif mutation.startswith("ENTRY_WITNESS_DROP:"):
        slot_id = mutation.split(":", 1)[1]
        for value in both: next(row for row in value["availability_slots"] if row["slot_id"] == slot_id).pop("witness_id", None)
    elif mutation.startswith("DERIVABILITY_CANARY:"):
        slot_id = mutation.split(":", 1)[1]
        for value in both:
            row = next(item for item in value["availability_slots"] if item["slot_id"] == slot_id)
            row["challenge"].update({"outcome": "REFUTED", "defining_predicate_result": "PASS"})
    else:
        raise AssertionDeath("MUTATION_NOOP", f"unknown or unimplemented mutation {mutation}")


def run_checks(sympy: dict[str, Any], wolfram: dict[str, Any]) -> str:
    require(
        sympy.get("schema_version") == SCHEMA and wolfram.get("schema_version") == SCHEMA
        and "semantic_digest" not in sympy and "semantic_digest" not in wolfram,
        "ASSERT_ENGINE_SCHEMA", "engine schema mismatch or dead self-digest field present",
    )
    require(sympy["semantic_view"] == wolfram["semantic_view"], "ASSERT_ENGINE_SEMANTIC_AGREEMENT", "full normalized semantic views differ")
    semantic = sympy["semantic_view"]
    s_registry = sympy["engine_local_route_registry"]
    w_registry = wolfram["engine_local_route_registry"]
    s_registry_ids = {row["semantic_route_id"] for row in s_registry if row["exists"] and row["executed"]}
    w_registry_ids = {row["semantic_route_id"] for row in w_registry if row["exists"] and row["executed"]}
    require([row["semantic_route_id"] for row in s_registry] == [row["semantic_route_id"] for row in w_registry], "ASSERT_ENGINE_ROUTE_IDS", "semantic route IDs differ")
    require(
        all(row["exists"] and row["executed"] for row in s_registry + w_registry)
        and all(left["engine_local_function"] != right["engine_local_function"] for left, right in zip(s_registry, w_registry)),
        "ASSERT_ENGINE_LOCAL_SYMBOLS", "engine-local route names did not exist/execute independently",
    )
    scope = semantic["scope"]
    require(
        scope == {
            "static_adjudication_only": True, "dynamical_selection_deferred": True,
            "one_body_only": True, "BVP_solved": False,
            "pair_configuration_carveout": "pair_annihilation_subquestion_only",
        }, "ASSERT_SCOPE_INVARIANTS", "static/one-body/BVP scope changed",
    )

    candidates = semantic["candidate_inventory"]
    require(candidates["mixture_generator_A"] == candidates["mixture_generator_B"] and len(candidates["mixture_generator_A"]) == 3, "ASSERT_MIXTURE_DUAL_GENERATORS", "mixture generators differ")
    axis = candidates["candidate_axis"]
    require(len(axis) == len(set(axis)) == candidates["candidate_count"], "ASSERT_CANDIDATE_AXIS_DISJOINT", "candidate axis is not a disjoint union")
    require(candidates["uncanonicalized_overlap_count"] == 0 and candidates["canonical_signature_count"] == len(axis), "ASSERT_CANDIDATE_CANONICALIZATION", "uncanonicalized candidate overlap")
    require(candidates["candidate_universe_digest"] == digest(candidates["candidate_records"]), "ASSERT_CANDIDATE_UNIVERSE_DIGEST", "candidate universe digest stale")
    closure = candidates["basis_closure"]
    require(closure["full_current_canonical_union_used"] and closure["status"] == "UNRESOLVED" and closure["catch_all_present"] and axis[-1] == "OTHER", "ASSERT_CATCHALL_PRESENCE", "catch-all presence does not follow residual closure")
    concrete = candidates["concrete_other_candidates"]
    signatures = {tuple(row["canonical_signature"]) for row in candidates["candidate_records"]}
    require(all(tuple(row.get("signature", [])) not in signatures for row in concrete), "ASSERT_CONCRETE_OTHER_CANONICAL", "concrete OTHER aliases a canonical class")
    require(all(
        row["membership_predicate"] == "positive_componentwise_operator_signature_equivalence"
        for row in candidates["candidate_records"] if row["candidate_id"] != "OTHER"
    ) and not next(row for row in candidates["candidate_records"] if row["candidate_id"] == "OTHER")["promotion_membership_eligible"], "ASSERT_MEMBERSHIP_PREDICATES", "membership predicate not positive normal form")
    require("immutable_successors" in candidates["amendment_rekey_duty"], "ASSERT_AMENDMENT_REKEY", "candidate-universe amendment re-key duty absent")

    taxonomy = semantic["evidence_taxonomy"]
    source = taxonomy["source_census"]
    source_ids = [row["source_id"] for row in source["records"]]
    require(
        len(source["complete_category_list"]) == 10 and set(source["category_counts"]) == set(source["complete_category_list"])
        and len(source_ids) == len(set(source_ids)) == len(source["source_ids"])
        and set(source_ids) == set(source["source_ids"]) and source["source_to_root_exact"],
        "ASSERT_SOURCE_CENSUS_COMPLETE", "external source census and taxonomy coverage differ",
    )
    grammar = taxonomy["constructor_grammar"]
    required_constructors = {
        "INCOMPATIBILITY", "NEGATIVE_CANDIDATE_MEMBERSHIP", "CASE_ELIMINATION", "COMPLEMENT_SURVIVOR_COUNT",
        "EXCLUSION_VERDICT", "POSTULATE_BRANCH", "STABILITY_DYNAMICAL_CLASS", "SOLVE_EVALUATION_RESULT",
    }
    require(required_constructors <= set(grammar) and "GENERIC" not in grammar and not taxonomy["generic_or_unclassified_constructor_allowed"], "ASSERT_TAXONOMY_CLOSED", "constructor grammar open or incomplete")
    artifacts = taxonomy["physics_bearing_stage0_artifacts"]
    for artifact in artifacts:
        recomputed = classify_content(artifact["normalized_inference_content"])
        require(
            recomputed == set(artifact["computed_constructors"])
            and inference_root_ids(artifact["normalized_inference_content"]) == set(artifact["reachable_root_ids"]),
            "ASSERT_CONTENT_CLASSIFICATION", f"content classifier/real DAG-root walk mismatch {artifact['artifact_id']}",
        )
    require(all(
        set(row["expected_root_ids"]) == set(row["reachable_root_ids"])
        and row["reachable_generation"] == "structural_walk_of_emitted_artifact_computation"
        and row["expected_generation"] == "normative_projection_of_complete_source_census"
        for row in artifacts
    ), "ASSERT_STAGE0_ARTIFACT_COVERAGE", "artifact independently expected/reachable ancestry differs")
    require(all("STABILITY_DYNAMICAL_CLASS" not in row["computed_constructors"] for row in artifacts), "ASSERT_STAGE0_ARTIFACT_ANCESTRY", "stability shaped physics-bearing stage0 artifact")
    banned_membership = {"NEGATIVE_CANDIDATE_MEMBERSHIP", "CASE_ELIMINATION"}
    require(all(not (set(row["computed_constructors"]) & banned_membership) for row in artifacts if row["affects_promotion_membership_or_uniqueness"]), "ASSERT_STAGE0_MEMBERSHIP_ANCESTRY", "elimination/negative membership shaped promotion-affecting artifact")
    fixtures = taxonomy["guard_fixture_dags"]
    allowed_promotion = set(taxonomy["promotion_allowed_constructors"])
    require(set(fixtures["promotion_positive"]["constructors"]) <= allowed_promotion, "ASSERT_PROMOTION_ANCESTRY_ALLOWED", "promotion ancestry used elimination/incompatibility/complement")
    require(all("STABILITY_DYNAMICAL_CLASS" not in row.get("constructors", []) for row in fixtures.values()), "ASSERT_PROGRAM_WIDE_STABILITY_BAN", "stability reached a U2 physics output")
    pair = fixtures["pair_configuration"]
    require(pair["firewall_tag"] == "PAIR_ANNIHILATION_ONLY" and pair["consumer"] == "topology_pair_annihilation_subquestion", "ASSERT_ONE_BODY_FIREWALL", "static pair configuration escaped carve-out")
    require("SOLVE_EVALUATION_RESULT" not in fixtures["posed_template"]["constructors"] and fixtures["posed_template"]["fields_unbound"], "ASSERT_BVP_BOUNDARY", "solve/evaluation reached template")
    require("POSTULATE_BRANCH" not in fixtures["promotion_positive"]["constructors"] and not fixtures["postulate_deferral_metadata"]["promotion_reachable"], "ASSERT_POSTULATE_PROMOTION_BAN", "postulate promoted evidence")
    require(all("MAXWELL" not in value and "COULOMB" not in value for row in fixtures.values() for value in row.get("root_types", [])), "ASSERT_NATIVE_ANCESTRY", "non-native/Maxwell root imported")

    obligations = semantic["obligation_censuses"]
    require(all(row["generator_A"] == row["generator_B"] for row in obligations.values()), "ASSERT_OBLIGATION_DUAL_GENERATORS", "obligation generators differ")
    require(all(not row["stratum_token_is_generation_input"] and "stratum" not in row["generation_inputs"] for row in obligations.values()), "ASSERT_CENSUS_STRATUM_FREE", "census generation circularly consumed stratum")
    require(all(
        row["anti_echo_predicate"] == "witness_satisfies_independent_nondefinitional_obligation"
        and row["semantic_ablation_criterion"] == "fail_nondefinitional_obligation_or_change_canonical_operator_class"
        for row in obligations.values()
    ), "ASSERT_NONDEGENERACY_ANTI_ECHO", "witness ablation is decorative/self-echoing")

    dependency = semantic["open_dependency_relation"]
    require(all(row["generator_A"] == row["generator_B"] for row in dependency["dependency_rows"]), "ASSERT_OPEN_DEPENDENCY_JOIN", "expected OPEN-dependency joins differ")
    require(
        dependency["generator_A_algorithm"] == "obligation_to_authoritative_slot_relational_join"
        and dependency["generator_B_algorithm"] == "raw_field_and_endpoint_channel_schema_walk"
        and dependency["generator_A_algorithm"] != dependency["generator_B_algorithm"]
        and not dependency["shared_task_code_between_generators"],
        "ASSERT_OPEN_DEPENDENCY_INDEPENDENCE", "OPEN dependency join routes are not independently constructed",
    )
    require(set(dependency["active_strata"]) == TILT_TYPES, "ASSERT_ACTIVE_STRATA", "active OPEN strata differ from frozen ledger")
    grid = semantic["grid_inventory"]
    expected_raw = sum(grid["ambient_count"] * count for count in grid["raw_strata_per_candidate"].values())
    require(grid["raw_ragged_cardinality"] == expected_raw == 144, "ASSERT_RAW_GRID_CARDINALITY", "raw ragged cardinality wrong")
    cell_ids = [row["cell_id"] for row in grid["grid_cells"]]
    expected_cells = {
        f"candidate={candidate}|ambient={ambient}|stratum={stratum}"
        for candidate in axis for ambient in ("one_sided_pathA29", "two_sided_R_w_postulate") for stratum in TILT_TYPES
    }
    require(set(cell_ids) == expected_cells and len(cell_ids) == len(set(cell_ids)), "ASSERT_GRID_COVERAGE", "grid exact coverage failed")
    require(all(row["timing"] == "pre-production" and row["stage0_objects_only"] for row in grid["collapse_proofs"]), "ASSERT_COLLAPSE_TIMING", "collapse used produced objects")
    require(grid["collapsed_cardinality"] == len(grid["grid_cells"]) == 144 and grid["preproduction_collapse_count"] == 0, "ASSERT_COLLAPSED_CARDINALITY", "collapsed cardinality/proof mismatch")
    contexts = grid["promotion_contexts"]
    require(len(contexts) == 16 and all({row["candidate_id"] for row in context["candidate_cell_mappings"]} == set(axis) and len(context["candidate_cell_mappings"]) == len(axis) for context in contexts), "ASSERT_PROMOTION_CONTEXT_COVERAGE", "promotion common-refinement coverage failed")

    inventory = semantic["route_fixture_inventory"]
    require(inventory["generated_preproduction"] and inventory["executed_route_match_rule"].startswith("production_route_set_exactly"), "ASSERT_ROUTE_INVENTORY_PREPRODUCTION", "route inventory not frozen preproduction")
    route_rows = inventory["route_records"]
    require(len(route_rows) == inventory["route_count"] == inventory["fixture_count"] == len(expected_cells) and {row["route_id"].removeprefix("route:") for row in route_rows} == {value.removeprefix("candidate=").replace("|ambient=", "|ambient=") for value in expected_cells}, "ASSERT_ROUTE_FIXTURE_COVERAGE", "cell-route-fixture coverage differs")
    require(all(
        matrix_residual(row["positive_fixture"]) == row["positive_fixture"]["residual"]
        == [0] * len(row["positive_fixture"]["rhs"])
        and row["positive_fixture"]["nondegenerate_norm_squared"] > 0
        and row["positive_fixture"]["expected"] == "ADMISSIBLE"
        and all(row["positive_fixture"]["independent_row_equations_satisfied"])
        and row["positive_fixture"]["known_outcome_generator_A"] != row["positive_fixture"]["known_outcome_generator_B"]
        for row in route_rows
    ), "ASSERT_FIXTURE_POSITIVE", "positive native-route fixture/independent outcome invalid")
    require(all(
        matrix_residual(row["malformed_fixture"]) == row["malformed_fixture"]["residual"]
        and any(value != 0 for value in row["malformed_fixture"]["residual"])
        and row["malformed_fixture"]["expected"] == "EXCLUDED"
        and row["malformed_fixture"]["semantic_route_id"] == row["positive_fixture"]["semantic_route_id"]
        for row in route_rows
    ), "ASSERT_FIXTURE_MALFORMED", "malformed fixture did not exclude on identical semantic route")
    require("nonholonomic" in obligations["E4"]["native_root_class"] and "rayleigh" in obligations["E5"]["native_root_class"], "ASSERT_E4_E5_NATIVE_ROOT", "E4/E5 reduced to bulk action")

    vocabulary = semantic["vocabulary_freeze"]
    uniform = vocabulary["uniform_record_rule"]
    guard_fixtures = semantic["guard_fixtures"]
    propagation = guard_fixtures["upstream_propagation"]
    require(
        len(uniform["record_types"]) == 11 and uniform["approval_requires_zero_integrity_failures"]
        and all(record_schema_valid(row) for row in guard_fixtures["record_schema_cases"])
        and any(row["integrity"] != "COMPUTATION_VALID" for row in propagation["required_upstreams"])
        and record_schema_valid(propagation["emitted_dependent"])
        and propagation["emitted_dependent"]["integrity"] == "NOT_RUN",
        "ASSERT_VOCAB_INTEGRITY_RULE", "executable integrity-conditional record-schema/propagation guard failed",
    )
    allowed_failure_reasons = {
        "ancestry_incomplete", "challenge_error", "contradictory_committed_derivations", "contradictory_evidence",
        "contradictory_forcing", "environment_identity_mismatch", "evaluated_code_closure_failure",
        "forced_class_outside_axis", "inconsistent_subresults", "schema_violation", "stale_candidate_universe",
        "uncanonicalized_candidate_overlap",
    }
    require(set(vocabulary["typed_failure_reasons"]) == allowed_failure_reasons and all(value.startswith("predicate:") for value in vocabulary["typed_failure_reasons"].values()), "ASSERT_TYPED_FAILURE_REASONS", "failure reason is free text/outside frozen vocabulary")
    require(set(vocabulary["obligation_evidence_state_enum"]) == {"SATISFIED", "DECISIVE_INCOMPATIBILITY", "MISSING", "INAPPLICABLE"} and vocabulary["evidence_classification"].startswith("transparent"), "ASSERT_EVIDENCE_STATE_FREEZE", "evidence-state predicates not frozen/transparent")
    require(len(vocabulary["disposition_precedence_table"]) == 4 and {row["landing"].split("(")[0] for row in vocabulary["disposition_precedence_table"]} == {"HARNESS_FAILED", "EXCLUDED", "ADMISSIBLE", "UNRESOLVED"}, "ASSERT_DISPOSITION_TABLE_TOTAL", "disposition table non-total")
    disposition_cases = guard_fixtures["disposition_cases"]
    require(all(
        row["emitted_state"] == classify_evidence(row["raw_predicate_inputs"])
        for case in disposition_cases for row in case["evidence_records"]
    ), "ASSERT_EVIDENCE_STATE_CLASSIFICATION", "real evidence-state predicate result relabelled/dropped")
    precedence_case = next(row for row in disposition_cases if row["case_id"] == "decisive_plus_missing")
    require(
        precedence_case["emitted_landing"] == evaluate_disposition(precedence_case) == "EXCLUDED",
        "ASSERT_DISPOSITION_PRECEDENCE", "executable precedence evaluator laundered decisive exclusion",
    )
    required_promotion_conditions = {
        "failed_required_upstream_set_nonempty", "uncanonicalized_overlap", "one_forced_class_ADMISSIBLE_closure_certified",
        "one_forced_class_ADMISSIBLE_closure_refused_valid", "one_forced_class_UNRESOLVED", "one_forced_class_EXCLUDED",
        "multiple_non_equivalent_forced_classes", "forced_class_outside_axis_or_stale", "no_forcing_witness",
    }
    require({row["condition"] for row in vocabulary["promotion_decision_table"]} == required_promotion_conditions, "ASSERT_PROMOTION_TABLE_TOTAL", "promotion decision table non-total")
    promotion_cases = guard_fixtures["promotion_cases"]
    not_run_case = next(row for row in promotion_cases if row["case_id"] == "failed_upstream")
    require(
        not_run_case["emitted_landing"] == evaluate_promotion(not_run_case["inputs"]) == "NOT_RUN(exact_set)",
        "ASSERT_PROMOTION_NOT_RUN_PROPAGATION", "NOT_RUN promotion carried physics",
    )
    require(all(
        row["emitted_landing"] == evaluate_promotion(row["inputs"])
        for row in promotion_cases if row["case_id"] != "failed_upstream"
    ), "ASSERT_PROMOTION_DECISION_TABLE", "executable promotion decision-table landing bypassed")

    require(len(vocabulary["ensemble_level_1_enum"]) == 4 and len(vocabulary["ensemble_level_2_enum"]) == 4, "ASSERT_ENSEMBLE_ENUMS", "ensemble enums not total/disjoint")
    cross_table = vocabulary["cross_level_decision_table"]
    cross_keys = {(row["applicability"], row["gate_integrity"], row["gate_outcome"]) for row in cross_table}
    expected_cross_keys = {
        (applicability, integrity, outcome)
        for applicability in ("geometric", "positively_non_geometric", "missing")
        for integrity, outcome in (
            ("COMPUTATION_VALID", "topologically-distinct"),
            ("COMPUTATION_VALID", "orientation-only"),
            ("COMPUTATION_VALID", "UNRESOLVED"), ("HARNESS_FAILED", None), ("NOT_RUN", None),
        )
    }
    require(
        len(cross_table) == 15 and cross_keys == expected_cross_keys
        and all(row["landing"] == evaluate_cross_level(row) for row in cross_table),
        "ASSERT_CROSS_LEVEL_TABLE", "applicability-first cross-level table is not input-total/executable",
    )
    cross_cases = guard_fixtures["cross_level_cases"]
    def cross_case(applicability: str, integrity: str, outcome: str | None) -> dict[str, Any]:
        return next(row for row in cross_cases if (row["applicability"], row["gate_integrity"], row["gate_outcome"]) == (applicability, integrity, outcome))
    mixed_geo = cross_case("geometric", "COMPUTATION_VALID", "topologically-distinct")
    mixed_nongeo = cross_case("positively_non_geometric", "COMPUTATION_VALID", "orientation-only")
    orientation = cross_case("geometric", "COMPUTATION_VALID", "orientation-only")
    gate_failed = cross_case("geometric", "HARNESS_FAILED", None)
    require(mixed_geo["landing"] == evaluate_cross_level(mixed_geo), "ASSERT_MIXED_GEOMETRIC_FIXTURE", "mixed geometric record suppressed topology gate")
    require(mixed_nongeo["landing"] == evaluate_cross_level(mixed_nongeo), "ASSERT_MIXED_NONGEOMETRIC_FIXTURE", "mixed non-geometric record entered topology gate")
    require(orientation["landing"] == evaluate_cross_level(orientation), "ASSERT_ORIENTATION_NEGATIVE_PRESERVED", "orientation-only laundered into incompletion")
    require(gate_failed["landing"] == evaluate_cross_level(gate_failed), "ASSERT_GATE_INTEGRITY_PROPAGATION", "gate integrity record carried level-2 physics")
    topology = vocabulary["topology_aggregate_table"]
    require(len(topology) == 9 and len({(row["sector_disconnection"], row["interpolation"]) for row in topology}) == 9 and all(row["landing"] == evaluate_topology(row) for row in topology), "ASSERT_TOPOLOGY_TABLE_TOTAL", "topology table non-total/non-executable")
    topology_cases = guard_fixtures["topology_cases"]
    require(all(row["emitted_landing"] == evaluate_topology(row) for row in topology_cases), "ASSERT_TOPOLOGY_INCONSISTENCY", "topology record disagrees with aggregate evaluator")
    require(all(not row["pair_used_in_aggregate"] for row in topology_cases) and vocabulary["pair_annihilation_aggregate_role"].startswith("orthogonal"), "ASSERT_PAIR_ORTHOGONAL", "pair-annihilation improperly polarized aggregate")
    require(all(
        ("nonholonomic" not in row["native_root_class"] or row["positive_fixture"]["native_structure_exercised"]["constraint_multiplier"])
        and ("rayleigh" not in row["native_root_class"] or row["positive_fixture"]["native_structure_exercised"]["dissipation_bookkeeping"])
        and (row["candidate_id"] not in {"E1", "E2", "E3"} or (row["positive_fixture"]["native_structure_exercised"]["boundary_variation"] and row["positive_fixture"]["native_structure_exercised"]["jump_condition"]))
        for row in route_rows
    ), "ASSERT_TOPOLOGY_NATIVE_ROOT", "native fixture omitted constraint/dissipation/jump-boundary structure")
    require(all("ensemble_exchange_discharge" in row["generator_A"] for row in obligations.values()), "ASSERT_ENSEMBLE_EXCHANGE_CONTRACT", "positive ensemble route lacks exchange/no-pairing discharge")

    closure_contract = semantic["closure_contract"]
    closure_fixture = guard_fixtures["closure"]
    closure_checks = closure_fixture_checks(closure_fixture)
    require(
        all(closure_checks.values()) and closure_contract["no_double_count"]
        and closure_contract["radiation_static_zero_first_class"] and closure_contract["coverage"].startswith("exact")
        and record_schema_valid(closure_fixture["closure_adjudication"])
        and (
            closure_fixture["closure_adjudication"]["integrity"] != "COMPUTATION_VALID"
            or closure_fixture["dependent_record"] == {"integrity": "COMPUTATION_VALID", "physics": "SELECTED"}
        ),
        "ASSERT_CLOSURE_GATE", "real closure census/certificate/ownership guard failed",
    )
    closure_adjudication = closure_fixture["closure_adjudication"]
    dependent = closure_fixture["dependent_record"]
    require(
        (closure_adjudication["integrity"] == "COMPUTATION_VALID" and record_schema_valid(dependent))
        or (closure_adjudication["integrity"] != "COMPUTATION_VALID" and dependent["integrity"] == "NOT_RUN" and dependent.get("physics") is None),
        "ASSERT_CLOSURE_PROPAGATION", "closure integrity failure laundered into physics refusal",
    )
    template = semantic["template_contract"]
    require(len(template["semantic_schema"]) == 6 and len(template["R49_exact_unresolved_reference_ids"]) == 8 and template["dependent_fields_unbound"], "ASSERT_TEMPLATE_SCHEMA", "template semantic schema incomplete")
    template_fixture = guard_fixtures["template"]
    expected_terms = expected_template_terms(template_fixture["committed_inputs"])
    reachable_terms = walk_template_terms(template_fixture["template_record"]["symbolic_ast"])
    require(expected_terms == template_fixture["expected_term_census"] and set(map(digest, expected_terms)) == set(map(digest, reachable_terms)), "ASSERT_TEMPLATE_TERM_INCIDENCE", "real posed-template term census/incidence missing or zeroed")
    require(not ast_contains_op(template_fixture["template_record"]["symbolic_ast"], "solve") and template_fixture["template_record"]["evaluation_state"] == "UNEVALUATED" and template["posing_not_solving"] and template["forbidden_ancestry_constructors"] == ["SOLVE_EVALUATION_RESULT"], "ASSERT_TEMPLATE_NO_SOLVE", "hidden solved content in actual template AST")
    ownership = semantic["return_closure_ownership"]
    require(not ownership["U2_owned"] and ownership["preserved_terminal"] == "UNRESOLVED(return_closure)", "ASSERT_RETURN_CLOSURE_NONOWNERSHIP", "U2 consumed return closure")

    slots = semantic["availability_slots"]
    require(all(
        row["availability_outcome"] == "DERIVED"
        for row in slots if row["slot_id"].startswith("candidate_definition:") and row["candidate_id"] != "OTHER"
    ), "ASSERT_DERIVED_OVERWRITE_CHALLENGE", "derived committed endpoint object overwritten as unresolved")
    require(all(
        row["availability_outcome"] == "DERIVED"
        for row in slots if row["slot_id"].startswith("mixture_law:")
    ), "ASSERT_FIRST_TIME_CONSTRUCTION", "constructible mixture law hidden on first construction")
    unresolved = [row for row in slots if row["availability_outcome"] == "UNRESOLVED"]
    derived = [row for row in slots if row["availability_outcome"] == "DERIVED"]
    allowed_derived_kinds = {
        "canonical_endpoint_condition_record", "committed_structure_mixture_law_record",
        "boundary_action_variation_record", "posed_template_cell_data",
    }
    require(all(
        row.get("derived_object", {}).get("object_kind") in allowed_derived_kinds
        and row["value_digest"] == digest(row["derived_object"])
        and row["dual_engine_object_derivation"] == "independent_engine_walk_of_committed_endpoint_and_channel_structure"
        and len(canonical_bytes(row["derived_object"])) > len(canonical_bytes({"slot_id": row["slot_id"], "producer_set": row["producer_set"]}))
        for row in derived
    ), "ASSERT_DERIVED_OBJECT_CONTENT", "DERIVED value digest does not bind actual committed-structure object content")
    for row in unresolved:
        slot_id = row["slot_id"]
        require(
            row.get("witness_id") == f"witness:{slot_id}" and row.get("challenge_id") == f"challenge:{slot_id}"
            and row.get("witness", {}).get("datum_id") == slot_id,
            f"ASSERT_ENTRY_WITNESS_RESOLVES:{slot_id}", "entry witness does not resolve to its ratified UNRESOLVED slot",
        )
        require(
            not (row["challenge"].get("outcome") == "REFUTED" and row["challenge"].get("defining_predicate_result") == "PASS"),
            f"ASSERT_STAGE0_UNRESOLVED_REFUTED:{slot_id}", "constructive challenge derived an UNRESOLVED datum",
        )
    summary = semantic["availability_summary"]
    require(len(slots) == 84 and len(derived) == summary["DERIVED"] and len(unresolved) == summary["UNRESOLVED"] and len(slots) == summary["total"], "ASSERT_AVAILABILITY_FLOOR", "availability floor/summary differ")
    require(all(
        row["integrity_state"] == "COMPUTATION_VALID" and row["witness"]["kind"] in ALLOWED_WITNESS_KINDS
        and row["witness"]["insufficiency_certificate"]["status"] == "PASS_COMPUTED"
        and row["witness"]["insufficiency_certificate"]["executed"]
        and row["witness"]["insufficiency_certificate"]["executed_semantic_route_id"] in s_registry_ids & w_registry_ids
        and row["witness"]["insufficiency_certificate"]["dual_engine_execution_required"] == ["SymPy", "Wolfram"]
        and row["witness"]["complete_committed_input_closure"] == row["witness"]["producer_set"]
        for row in unresolved
    ), "ASSERT_WITNESS_SCHEMA", "unavailability witness empty/ill-typed/prose")
    require(all(
        row["challenge"]["outcome"] in {"CONSTRUCTIVE_FAIL", "REFUTED"}
        and not row["challenge"]["empty_output"] and not row["challenge"]["ill_typed_by_fiat"]
        and row["challenge"]["candidate_is_well_typed"] and row["challenge"]["dual_engine_certificate"]
        and row["challenge"]["executed_semantic_route_id"] in s_registry_ids & w_registry_ids
        and row["challenge"]["dual_engine_execution_required"] == ["SymPy", "Wolfram"]
        and row["challenge"]["dag_shared_nodes_are_committed_inputs_only"]
        for row in unresolved
    ), "ASSERT_CHALLENGE_SCHEMA", "challenge terminal/certificate invalid")
    require(all(row["witness"]["kind"] == row["challenge"]["kind"] for row in unresolved), "ASSERT_WITNESS_CHALLENGE_LOCKSTEP", "witness/challenge kind mismatch")
    require(all(row["witness"]["counterfactual_restore_mutation"]["baseline_status"] == "PASS_COMPUTED" and row["witness"]["counterfactual_restore_mutation"]["restored_status"] == "FAIL_COMPUTED" for row in unresolved), "ASSERT_WITNESS_RESTORE", "restore target did not kill insufficiency")
    dimensional = semantic["dimensional_firewall"]
    require(dimensional["base_pass"] and dimensional["ablation_fired"] and dimensional["surface_power_dimensions"] == dimensional["expected_power_dimensions"] == [2, -3, 1], "ASSERT_DIMENSIONAL_FIREWALL", "inline dimensional firewall/ablation failed")
    standards = semantic["standard_bindings"]
    expected_standards = {f"S-{index}" for index in range(1, 7)} | {f"P-{index}" for index in range(1, 11)}
    require(
        {row["standard_id"] for row in standards} == expected_standards and len(standards) == 16
        and len({row["reachable_check_id"] for row in standards}) == 16,
        "ASSERT_STANDARD_COVERAGE", "sixteen-standard binding incomplete/non-unique",
    )
    for row in standards:
        standard_id = row["standard_id"]
        expected_check = f"ASSERT_STANDARD_{standard_id.replace('-', '_')}"
        require(
            row == {
                "standard_id": standard_id, "authoritative_text": row["authoritative_text"],
                "acceptance_predicate_id": f"predicate:standard:{standard_id}",
                "reachable_check_id": expected_check,
                "tooth_id": f"STANDARD_TOOTH:{standard_id}",
                "evidence_id": f"evidence:standard:{standard_id}",
            } and standard_predicate(standard_id, semantic),
            expected_check, f"standard {standard_id} executable predicate/binding failed",
        )
    branch_ids = [row["stable_branch_id"] for row in grid["grid_cells"]]
    require(all(branch_ids) and len(branch_ids) == len(set(branch_ids)) == len(grid["grid_cells"]), "ASSERT_BRANCH_ID_COVERAGE", "stable branch IDs incomplete/duplicate")
    require(summary["integrity_failures"] == 0, "ASSERT_ZERO_STAGE0_INTEGRITY_FAILURES", "stage0 approval contains integrity failure")
    banking_fixture = guard_fixtures["banking_fixture"]
    require(
        banking_fixture["emitted_approval"] == all(row["integrity"] == "COMPUTATION_VALID" for row in banking_fixture["records"])
        and not banking_fixture["emitted_approval"],
        "ASSERT_ALL_INTEGRITY_REJECTION", "all-integrity failure records reached approval/banking",
    )

    catalog = mutation_catalog(semantic)
    catalog_asserts = {row["expected_assert_id"] for row in catalog}
    prospective = set(CHECK_ORDER) | {"ASSERT_MUTATION_CATALOG_COVERAGE"}
    require(semantic["production_contract"]["mutation_catalog_coverage_required"] and prospective <= catalog_asserts, "ASSERT_MUTATION_CATALOG_COVERAGE", f"checks without teeth: {sorted(prospective - catalog_asserts)}")
    return digest(semantic)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--sympy", required=True)
    parser.add_argument("--wolfram", required=True)
    parser.add_argument("--output")
    parser.add_argument("--mutation")
    args = parser.parse_args()
    try:
        sympy = copy.deepcopy(load_yaml(Path(args.sympy)))
        wolfram = copy.deepcopy(load_yaml(Path(args.wolfram)))
        apply_mutation(args.mutation, sympy, wolfram)
        semantic_digest = run_checks(sympy, wolfram)
        if args.mutation:
            raise AssertionDeath("MUTATION_NOOP", f"mutation {args.mutation} survived every assert")
        semantic = sympy["semantic_view"]
        require(bool(args.output), "ASSERT_COMPARATOR_OUTPUT", "missing output path")
        result = {
            "schema_version": "U2_STAGE0_ENGINE_AGREEMENT_V1", "status": "ENGINE_AGREE",
            "semantic_digest": semantic_digest, "comparison": "deep equality of full normalized semantic view",
            "checks": [{"assert_id": value, "status": "PASS", "reachable": True} for value in CHECK_ORDER],
            "mutation_catalog": mutation_catalog(semantic),
            "summary": {
                "candidate_count": semantic["candidate_inventory"]["candidate_count"],
                "raw_grid_cardinality": semantic["grid_inventory"]["raw_ragged_cardinality"],
                "collapsed_grid_cardinality": semantic["grid_inventory"]["collapsed_cardinality"],
                "route_count": semantic["route_fixture_inventory"]["route_count"],
                "availability": semantic["availability_summary"],
            },
        }
        output = Path(args.output)
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(yaml.safe_dump(result, sort_keys=False, allow_unicode=True, width=140), encoding="utf-8")
        print(f"U2_STAGE0_ENGINE_AGREE digest={semantic_digest} checks={len(CHECK_ORDER)} teeth={len(result['mutation_catalog'])}")
        return 0
    except AssertionDeath as failure:
        print(f"ASSERTION_FAILED {failure.assert_id}: {failure.detail}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
