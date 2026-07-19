#!/usr/bin/env python3
"""Independent SymPy production engine for U2 boundary adjudication.

The engine consumes the ratified stage-0 bundle without changing it.  Every
physics landing is obtained from the executable stage-0 decision functions;
OPEN data stay OPEN and are represented by typed, executed challenges.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from collections import Counter
from pathlib import Path
from typing import Any, Iterable

import sympy as sp
import yaml


RATIFIED_DIGEST = "b01a1821e908589c3698512bbb9aff874b721af6dcbfa1c3b8b1f8d33119b32b"
STAGE0_COMPONENTS = {
    "frozen_data_pin_table": "frozen_data_pin_table.yaml",
    "candidate_inventory": "candidate_inventory.yaml",
    "obligation_censuses": "obligation_censuses.yaml",
    "dependency_grid_inventory": "dependency_grid_inventory.yaml",
    "vocabulary_freeze": "vocabulary_freeze.yaml",
    "evidence_taxonomy": "evidence_taxonomy.yaml",
    "availability_slots": "availability_slots.yaml",
    "route_fixture_inventory": "route_fixture_inventory.yaml",
    "closure_template_contracts": "closure_template_contracts.yaml",
    "environment_identity": "environment_identity.yaml",
    "standard_bindings": "standard_bindings.yaml",
    "producer_map": "producer_map.yaml",
    "evaluated_code_closure_policy": "evaluated_code_closure_policy.yaml",
    "parameter_register_proposals": "parameter_register_proposals.yaml",
    "obligation_manifest": "obligation_manifest.yaml",
}

SOURCE = Path(__file__).resolve().parent
if str(SOURCE) not in sys.path:
    sys.path.insert(0, str(SOURCE))

# These imports are deliberate: production verdicts use the ratified executable
# tables, rather than a look-alike implementation in this engine.
from u2_stage0_sympy import (  # noqa: E402
    classify_evidence as frozen_classify_evidence,
    classify_inference_content as frozen_classify_inference_content,
    evaluate_cross_level as frozen_evaluate_cross_level,
    evaluate_disposition as frozen_evaluate_disposition,
    evaluate_promotion as frozen_evaluate_promotion,
    evaluate_topology as frozen_evaluate_topology,
    record_schema_valid as frozen_record_schema_valid,
)


class EngineFailure(RuntimeError):
    pass


def require(condition: bool, detail: str) -> None:
    if not condition:
        raise EngineFailure(detail)


def require_generation(condition: bool, assert_id: str, detail: str) -> None:
    if not condition:
        raise EngineFailure(f"ASSERTION_FAILED {assert_id}: {detail}")


def load_yaml(path: Path) -> Any:
    with path.open("rb") as handle:
        return yaml.load(handle, Loader=yaml.CSafeLoader)


def dump_yaml(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        yaml.safe_dump(value, sort_keys=False, allow_unicode=True, width=160),
        encoding="utf-8",
    )


def canonical_bytes(value: Any) -> bytes:
    return json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(",", ":")).encode()


def digest(value: Any) -> str:
    return hashlib.sha256(canonical_bytes(value)).hexdigest()


def live_physics_record_projection(document: dict[str, Any]) -> list[dict[str, Any]]:
    """Independently realize the stage-0-frozen non-template projection."""
    rows: list[dict[str, Any]] = []
    context_fields = ("cell_id", "stable_branch_id", "candidate_id", "ambient", "stratum")
    disposition_fields = (
        *context_fields, "native_root_class", "integrity", "expected_dependencies",
        "used_dependencies", "obligation_evidence", "disposition",
        "disposition_evaluator_landing", "unavailability",
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


def ids(values: Iterable[str]) -> list[str]:
    return sorted(set(values), key=lambda value: (value.lower(), value))


def typed_root_id(value: str) -> str:
    if value.startswith("source:"):
        return value
    if value.startswith(("tilt:", "domain:", "endpoint:")):
        return f"source:phaseC_slot:{value}"
    if value in {"open_leaf:gammaSigma", "open_leaf:tangentDtN"}:
        return f"source:coefficient:{value.split(':', 1)[1]}"
    if value.startswith("open_leaf:"):
        return f"source:field:{value.split(':', 1)[1]}"
    raise EngineFailure(f"no ratified typed-root mapping for {value}")


def root_ast(root_ids: Iterable[str], operation: str) -> dict[str, Any]:
    roots = ids(typed_root_id(value) for value in root_ids)
    return {"op": operation, "args": [{"op": "root", "id": value} for value in roots]}


def proof_dag(root_ids: Iterable[str], operation: str) -> dict[str, Any]:
    content = root_ast(root_ids, operation)
    return {
        "normalized_inference_content": content,
        "constructors": frozen_classify_inference_content(content),
        "root_ids": [row["id"] for row in content["args"]],
        "classification_source": "ratified_normalized_content_classifier",
    }


def load_bundle(bundle_dir: Path, supplied_digest: str) -> dict[str, Any]:
    require(supplied_digest == RATIFIED_DIGEST, "STAGE0_CONTRACT_DIGEST differs from ratified U2 digest")
    components = {
        name: load_yaml(bundle_dir / filename) for name, filename in STAGE0_COMPONENTS.items()
    }
    require(digest(components) == RATIFIED_DIGEST, "ratified component bundle digest changed")
    require(
        load_yaml(bundle_dir / "stage0_bundle.yaml")["stage0_contract_digest"] == RATIFIED_DIGEST
        and load_yaml(bundle_dir / "stage0_contract.yaml")["stage0_contract_digest"] == RATIFIED_DIGEST,
        "stage-0 bundle/contract binding changed",
    )
    return components


def phasec_availability(repo: Path) -> dict[str, dict[str, Any]]:
    path = repo / (
        "software/em_charge_attribute/reports/u1_body_dynamics_artifacts/"
        "stage_c_0_tilt_coupling_contract/availability_slots.yaml"
    )
    return {row["slot_id"]: row for row in load_yaml(path)["slots"]}


def dependency_slot_id(token: str) -> str:
    return token


def dependency_measurement(token: str, phasec_slots: dict[str, dict[str, Any]]) -> dict[str, Any]:
    slot_id = dependency_slot_id(token)
    require(slot_id in phasec_slots, f"expected OPEN dependency has no frozen availability slot: {token}")
    row = phasec_slots[slot_id]
    disposition = row["disposition"]
    return {
        "dependency_token": token,
        "ratified_slot_id": slot_id,
        "availability": disposition,
        "integrity_state": "COMPUTATION_VALID",
        "witness_id": row.get("witness_id"),
        "challenge_id": row.get("challenge_id"),
        "missing": disposition == "UNRESOLVED",
        "measurement_source": "ratified_PhaseC_availability_record",
    }


def obligation_dependencies(obligation: str, expected: list[str], candidate: str) -> list[str]:
    candidate_native = [
        token for token in expected
        if token.startswith("endpoint:") or token in {
            "open_leaf:E4_shear_lock", "open_leaf:E5_rayleigh",
            "open_leaf:gammaSigma", "open_leaf:tangentDtN",
        }
    ]
    if obligation == "host_location_evidence":
        return [token for token in expected if token in {"open_leaf:geon_core_bundle", "open_leaf:sleeve_core_trace"}]
    if obligation.startswith("topology_"):
        return ids(["open_leaf:sleeve_core_trace", *candidate_native])
    if obligation == "jump_and_trace_compatibility":
        return [
            token for token in expected
            if token.startswith("tilt:") or token in {
                "domain:Sigma_boundary_data", "open_leaf:throat_surface_functional",
                "open_leaf:outer_surface_functional",
            }
        ]
    if obligation == "mechanical_closure_contribution_census":
        return [token for token in expected if token in {"open_leaf:geon_core_bundle", "open_leaf:outer_surface_functional"}]
    if obligation in {"native_root_structure", "composite_native_root_compatibility"}:
        return candidate_native
    if obligation.startswith(("E1_", "E2_", "E3_")):
        return [token for token in expected if token in {"domain:Sigma_boundary_data", "open_leaf:throat_surface_functional"}]
    if obligation.startswith(("E4_", "E5_")):
        return candidate_native
    if obligation in {"operator_definition_or_complement_closure", "whole_complement_class_coverage"}:
        return [token for token in expected if token in {"open_leaf:throat_surface_functional", "open_leaf:outer_surface_functional"}]
    if obligation == "boundary_variation_or_virtual_work" and candidate in {"E3", "OTHER"}:
        return [token for token in expected if token in {"open_leaf:throat_surface_functional", "open_leaf:outer_surface_functional"}]
    if obligation in {"ensemble_exchange_discharge", "geometric_refinement_applicability"} and candidate in {"E3", "OTHER"}:
        return [token for token in expected if token == "open_leaf:sleeve_core_trace"]
    return []


def obligation_stage0_slots(obligation: str, candidate: str, stratum: str) -> list[str]:
    candidate_slot = f"candidate_definition:{candidate}"
    if obligation.startswith("topology_"):
        suffix = {
            "topology_finite_energy_interpolation": "finite_energy_interpolation",
            "topology_pair_annihilation_path": "pair_annihilation",
            "topology_sector_disconnection": "sector_disconnection",
        }[obligation]
        return [f"topology:{candidate}:{suffix}"]
    if obligation == "host_location_evidence":
        return [f"host_location:{candidate}"]
    if obligation == "mechanical_closure_contribution_census":
        return [f"mechanical_closure:{candidate}"]
    if obligation in {
        "boundary_variation_or_virtual_work", "ensemble_exchange_discharge",
        "geometric_refinement_applicability",
    }:
        return [f"ensemble:{candidate}:boundary_action_variation"]
    if obligation == "template_constituent_contract":
        return [f"template_cell_specific:{candidate}", f"template_free_data:{stratum}"]
    if obligation.startswith("mixture_law:"):
        return [obligation]
    if obligation in {"operator_definition_or_complement_closure", "whole_complement_class_coverage"}:
        return [candidate_slot, "basis_closure"]
    return [candidate_slot]


def stage0_slot_measurement(slot_id: str, availability: dict[str, dict[str, Any]]) -> dict[str, Any]:
    row = availability[slot_id]
    outcome = row["availability_outcome"]
    return {
        "stage0_slot_id": slot_id,
        "availability": outcome,
        "integrity_state": row["integrity_state"],
        "missing": outcome == "UNRESOLVED",
        "decisive_incompatibility": outcome == "EXCLUDED",
        "witness_id": row.get("witness_id"),
        "challenge_id": row.get("challenge_id"),
        "producer_set": row["producer_set"],
        "measurement_source": "ratified_U2_stage0_availability_record",
    }


def evidence_record(
    obligation: str,
    candidate: str,
    stratum: str,
    expected_dependencies: list[str],
    phasec_slots: dict[str, dict[str, Any]],
    availability: dict[str, dict[str, Any]],
) -> dict[str, Any]:
    dependencies = obligation_dependencies(obligation, expected_dependencies, candidate)
    measurements = [dependency_measurement(value, phasec_slots) for value in dependencies]
    slot_measurements = [
        stage0_slot_measurement(value, availability)
        for value in obligation_stage0_slots(obligation, candidate, stratum)
    ]
    applicable = obligation != "template_constituent_contract"
    all_measurements = [*measurements, *slot_measurements]
    missing = applicable and any(row["missing"] for row in all_measurements)
    committed_incompatibility = applicable and any(
        row.get("decisive_incompatibility", False)
        or row["availability"] == "EXCLUDED" for row in all_measurements
    )
    ancestry_complete_and_typed = bool(all_measurements) and all(
        row["integrity_state"] == "COMPUTATION_VALID" and row["measurement_source"]
        for row in all_measurements
    )
    raw = {
        "applicable": applicable,
        "committed_incompatibility": committed_incompatibility,
        "ancestry_complete_and_typed": ancestry_complete_and_typed,
        "datum_missing": missing,
    }
    state = frozen_classify_evidence(raw)
    operation = "unavailability" if state == "MISSING" else "positive_join"
    roots = [
        *[row["dependency_token"] for row in measurements],
        *[root for row in slot_measurements for root in row["producer_set"]],
    ]
    return {
        "obligation_id": obligation,
        "integrity": "COMPUTATION_VALID",
        "raw_predicate_inputs": raw,
        "predicate_measurements": measurements,
        "stage0_slot_measurements": slot_measurements,
        "emitted_state": state,
        "proof_dag": proof_dag(roots, operation),
        "state_classifier": "u2_stage0_sympy.classify_evidence",
    }


MINIMAL_DATUM = {
    "E1": "throat_surface_functional:E1_holonomic_placement",
    "E2": "throat_surface_functional:E2_free_slip_placement",
    "E3": "boundary_action_variation:E3_bulk_texture",
    "E4": "sleeve_core_trace:E4_nonholonomic_constraint",
    "E5": "E5_rayleigh_dissipation_functional",
    "MIXTURE(F_E1_E4)": "sleeve_core_trace:F_E1_E4_composite_law",
    "MIXTURE(F_E2_E4)": "sleeve_core_trace:F_E2_E4_composite_law",
    "MIXTURE(F_E4_E5)": "sleeve_core_trace:F_E4_E5_constraint_Rayleigh_law",
    "OTHER": "operator_definition_or_complement_closure",
}


def unavailability_payload(
    candidate: str,
    stratum: str,
    missing_dependencies: list[dict[str, Any]],
    missing_slots: list[dict[str, Any]],
) -> dict[str, Any]:
    require(bool(missing_slots), f"UNRESOLVED cell for {candidate} has no ratified missing slot")
    primary = sorted(missing_slots, key=lambda row: row["stage0_slot_id"])[0]
    primary_slot = primary["record"]
    certificate = primary_slot["witness"]["insufficiency_certificate"]
    restore = primary_slot["witness"]["counterfactual_restore_mutation"]
    challenge = primary_slot["challenge"]
    slot_references = [{
        "stage0_slot_id": row["stage0_slot_id"],
        "witness_id": row["record"]["witness_id"],
        "challenge_id": row["record"]["challenge_id"],
        "availability_outcome": row["record"]["availability_outcome"],
    } for row in sorted(missing_slots, key=lambda value: value["stage0_slot_id"])]
    dependency_references = [{
        "dependency_token": row["dependency_token"],
        "ratified_slot_id": row["ratified_slot_id"],
        "witness_id": row["witness_id"],
        "challenge_id": row["challenge_id"],
        "availability": row["availability"],
    } for row in sorted(missing_dependencies, key=lambda value: value["dependency_token"])]
    return {
        "witness_id": primary_slot["witness_id"],
        "challenge_id": primary_slot["challenge_id"],
        "kind": primary_slot["witness"]["kind"],
        "named_datum": MINIMAL_DATUM[candidate],
        "stratum_datum": f"tilt:{stratum}",
        "missing_dependency_tokens": [row["dependency_token"] for row in dependency_references],
        "missing_stage0_slot_ids": [row["stage0_slot_id"] for row in slot_references],
        "referenced_phaseC_slots": dependency_references,
        "referenced_stage0_slots": slot_references,
        "reference_set_exact": (
            {row["dependency_token"] for row in dependency_references}
            == {row["dependency_token"] for row in missing_dependencies}
            and {row["stage0_slot_id"] for row in slot_references}
            == {row["stage0_slot_id"] for row in missing_slots}
        ),
        "required_type": primary_slot["required_type"],
        "required_dimensions": primary_slot["required_dimensions"],
        "challenge": {
            "executed": certificate["executed"],
            "semantic_route_id": certificate["executed_semantic_route_id"],
            "attempted_candidate_count": challenge["attempted_candidate_count"],
            "measured_rank": certificate["measured_rank"],
            "measured_nullity": certificate["measured_nullity"],
            "restore_rank": restore["restored_rank"],
            "restore_nullity": restore["restored_nullity"],
            "outcome": challenge["outcome"],
            "empty_output": challenge["empty_output"],
            "candidate_well_typed": challenge["candidate_is_well_typed"],
            "dual_engine_certificate": challenge["dual_engine_certificate"],
        },
        "proof_dag": proof_dag([
            *[row["dependency_token"] for row in dependency_references],
            *primary_slot["producer_set"],
        ], "challenge"),
    }


def exact_residual(matrix: list[list[int]], vector: list[int], rhs: list[int]) -> list[int]:
    return [int(value) for value in sp.Matrix(matrix) * sp.Matrix(vector) - sp.Matrix(rhs)]


def route_control(row: dict[str, Any]) -> dict[str, Any]:
    positive = row["positive_fixture"]
    malformed = row["malformed_fixture"]
    matrix = positive["matrix"]
    candidate = positive["candidate"]
    rhs = positive["rhs"]
    positive_residual = exact_residual(matrix, candidate, rhs)
    malformed_residual = exact_residual(matrix, malformed["candidate"], malformed["rhs"])
    ablated = list(candidate)
    first_nonzero = next(index for index, value in enumerate(ablated) if value != 0)
    ablated[first_nonzero] = 0
    ablated_residual = exact_residual(matrix, ablated, rhs)
    native = positive["native_structure_exercised"]
    dissipative = bool(native.get("dissipation_bookkeeping"))
    if dissipative:
        coordinate_index = 3 if len(candidate) == 4 else 1
        coefficient = int(native["rayleigh_coefficient"])
        velocity = int(candidate[coordinate_index])
        dissipated_power = coefficient * velocity ** 2
    else:
        coordinate_index = None
        coefficient = None
        velocity = None
        dissipated_power = None
    return {
        "route_id": row["route_id"],
        "cell_id": "candidate=" + row["route_id"].removeprefix("route:"),
        "fixture_id": row["fixture_id"],
        "semantic_route_id": positive["semantic_route_id"],
        "test_only": True,
        "positive_residual": positive_residual,
        "positive_admissible": all(value == 0 for value in positive_residual),
        "nondegenerate_norm_squared": sum(value * value for value in candidate),
        "native_structure_exercised": native,
        "malformed_residual": malformed_residual,
        "malformed_excluded": any(value != 0 for value in malformed_residual),
        "semantic_ablation": {
            "removed_coordinate_index": first_nonzero,
            "ablated_candidate": ablated,
            "ablated_residual": ablated_residual,
            "nondefinitional_obligation_failed": any(value != 0 for value in ablated_residual),
            "criterion": "fail_nondefinitional_obligation_or_change_canonical_operator_class",
        },
        "dissipation_measurement": {
            "applicable": dissipative,
            "rayleigh_coefficient": coefficient,
            "velocity_coordinate_index": coordinate_index,
            "velocity_coordinate": velocity,
            "formula": "rayleigh_coefficient*velocity_coordinate^2" if dissipative else None,
            "value": dissipated_power,
            "frozen_fixture_field_equal": (
                dissipated_power == native.get("computed_dissipated_power") if dissipative else
                native.get("computed_dissipated_power") is None
            ),
        },
    }


CONJUGATE_SIGNATURE_ROLES = {
    "normal_velocity_control": "normal_traction_response",
    "normal_traction_response": "normal_velocity_control",
    "tangential_velocity_control": "tangential_traction_response",
    "tangential_traction_response": "tangential_velocity_control",
}

CANONICAL_CONJUGATE_EXPECTATION = {
    "normal_velocity_lock": "normal_traction_response",
    "tangential_velocity_lock": "tangential_traction_response",
    "tangential_traction_free": "tangential_velocity_control",
}


def static_signature_classifier(boundary_variation: dict[str, Any]) -> dict[str, Any]:
    """Classify controls/responses by walking the committed boundary variation."""
    signature: list[str] = []
    measurements: list[dict[str, Any]] = []
    variation_channels: list[str] = []
    flux_channels: list[str] = []
    for component in boundary_variation["component_variations"]:
        condition = component["boundary_condition"]
        component_roles: list[str] = []
        if "v.normal=V.normal" in condition:
            component_roles.append("normal_velocity_control")
        elif "normal_traction" in condition:
            component_roles.append("normal_traction_response")
        if "v.tangent=V.tangent" in condition:
            component_roles.append("tangential_velocity_control")
        elif "tangential_traction=0" in condition:
            component_roles.append("tangential_traction_response")
        signature.extend(component_roles)
        variation_channels.extend(component["variation_channels"])
        flux_channels.extend(component["flux_channels"])
        measurements.append({
            "endpoint_id": component["endpoint_id"],
            "boundary_condition": condition,
            "classified_roles": ids(component_roles),
            "variation_channels": ids(component["variation_channels"]),
            "flux_channels": ids(component["flux_channels"]),
        })
    signature = ids(signature)
    return {
        "classifier": "committed_boundary_variation_static_signature_walk_v1",
        "derived_object_kind": boundary_variation["object_kind"],
        "component_count": len(boundary_variation["component_variations"]),
        "component_measurements": measurements,
        "signature": signature,
        "signature_complete": bool(signature) and all(
            len(measurement["classified_roles"]) >= len(measurement["variation_channels"])
            for measurement in measurements
        ),
        "variation_channels": ids(variation_channels),
        "flux_channels": ids(flux_channels),
    }


def native_pairing_analysis(candidate_structure: dict[str, Any]) -> dict[str, Any]:
    """Independently construct the conjugate expectation from native-root data."""
    components = candidate_structure["components"]
    canonical_signature = candidate_structure["canonical_signature"]
    eligible_components = ids(
        row["endpoint_id"] for row in components
        if row["variational_class"] == "holonomic_field_BC"
    )
    nonconjugate_components = ids(
        f"{row['endpoint_id']}:{row['variational_class']}" for row in components
        if row["variational_class"] != "holonomic_field_BC"
    )
    unmapped_signature_terms = ids(
        value for value in canonical_signature if value not in CANONICAL_CONJUGATE_EXPECTATION
    )
    complete = (
        candidate_structure["native_root_class"] == "variational_holonomic"
        and len(eligible_components) == len(components)
        and not nonconjugate_components
        and not unmapped_signature_terms
        and bool(canonical_signature)
    )
    expected = ids(
        CANONICAL_CONJUGATE_EXPECTATION[value]
        for value in canonical_signature if value in CANONICAL_CONJUGATE_EXPECTATION
    ) if complete else []
    return {
        "construction": "native_root_class_and_canonical_signature_pairing_walk_v1",
        "native_root_class": candidate_structure["native_root_class"],
        "required_component_count": len(components),
        "eligible_variational_components": eligible_components,
        "nonconjugate_components": nonconjugate_components,
        "unmapped_canonical_signature_terms": unmapped_signature_terms,
        "complete_structure_pairing": complete,
        "independently_expected_exchanged_signature": expected,
    }


def local_ensemble(
    ensemble_slot: dict[str, Any], candidate_structure_slot: dict[str, Any]
) -> dict[str, Any]:
    structure_id = ensemble_slot["slot_id"]
    if ensemble_slot["availability_outcome"] != "DERIVED":
        return {
            "available": False,
            "classification": None,
            "committed_structure_id": structure_id,
            "geometric_component": None,
            "exchange_route": None,
            "classifier_evidence": None,
            "pairing_analysis": None,
        }
    boundary_variation = ensemble_slot["derived_object"]
    classifier = static_signature_classifier(boundary_variation)
    pairing = native_pairing_analysis(candidate_structure_slot["derived_object"])
    pure_displacement = (
        classifier["signature_complete"]
        and classifier["signature"] == [
            "normal_velocity_control", "tangential_velocity_control",
        ]
    )
    classification = (
        "fixed-displacement/geometric" if pure_displacement
        else f"mixed/other-ensemble({structure_id})"
    )
    return {
        "available": True,
        "classification": classification,
        "committed_structure_id": structure_id,
        "geometric_component": bool(classifier["variation_channels"]),
        "exchange_route": (
            "native_conjugate_transform" if pairing["complete_structure_pairing"]
            else "computed_no_pairing_certificate"
        ),
        "classifier_evidence": classifier,
        "pairing_analysis": pairing,
    }


def exchange_control(
    local: dict[str, Any], generation_mutation: str | None = None
) -> dict[str, Any] | None:
    if not local["available"]:
        return None
    route = local["exchange_route"]
    classifier = local["classifier_evidence"]
    pairing = local["pairing_analysis"]
    if route == "native_conjugate_transform":
        baseline = classifier["signature"]
        exchanged = ids(CONJUGATE_SIGNATURE_ROLES[value] for value in baseline)
        expected = pairing["independently_expected_exchanged_signature"]
        if generation_mutation == "TOOTH_EXCHANGE_EXPECTED_SIGNATURE_GENERATION":
            expected = expected[:-1]
        require_generation(
            classifier["signature_complete"] and exchanged == expected and baseline != exchanged,
            "ASSERT_EXCHANGE_SIGNATURE_GENERATION",
            f"boundary-variation exchange disagrees with native pairing for {local['committed_structure_id']}",
        )
        return {
            "route": route,
            "baseline_signature": baseline,
            "computed_exchanged_signature": exchanged,
            "independently_generated_expected_signature": expected,
            "response_character_flipped": baseline != exchanged and exchanged == expected,
            "signature_classifier_evidence": classifier,
            "native_pairing_construction": pairing,
        }
    insufficiency = {
        "construction": pairing["construction"],
        "native_root_class": pairing["native_root_class"],
        "required_component_count": pairing["required_component_count"],
        "eligible_variational_components": pairing["eligible_variational_components"],
        "nonconjugate_components": pairing["nonconjugate_components"],
        "unmapped_canonical_signature_terms": pairing["unmapped_canonical_signature_terms"],
        "complete_structure_pairing": pairing["complete_structure_pairing"],
    }
    if generation_mutation == "TOOTH_NO_PAIRING_CERTIFICATE_GENERATION":
        insufficiency["nonconjugate_components"] = []
        insufficiency["unmapped_canonical_signature_terms"] = []
    certificate = (
        not insufficiency["complete_structure_pairing"]
        and len(insufficiency["eligible_variational_components"])
        < insufficiency["required_component_count"]
        and bool(
            insufficiency["nonconjugate_components"]
            or insufficiency["unmapped_canonical_signature_terms"]
        )
    )
    require_generation(
        certificate,
        "ASSERT_NO_PAIRING_CERTIFICATE_GENERATION",
        f"no-pairing insufficiency is incomplete for {local['committed_structure_id']}",
    )
    return {
        "route": route,
        "observed_boundary_signature": classifier["signature"],
        "signature_classifier_evidence": classifier,
        "pairing_insufficiency": insufficiency,
        "no_pairing_certificate": certificate,
    }


def topology_record(cell: dict[str, Any], availability: dict[str, dict[str, Any]]) -> dict[str, Any]:
    candidate = cell["candidate_id"]
    slot_ids = {
        "sector_disconnection": f"topology:{candidate}:sector_disconnection",
        "finite_energy_interpolation": f"topology:{candidate}:finite_energy_interpolation",
        "pair_annihilation": f"topology:{candidate}:pair_annihilation",
    }
    require(all(availability[value]["availability_outcome"] == "UNRESOLVED" for value in slot_ids.values()), "topology datum unexpectedly available")
    sector = "UNRESOLVED"
    interpolation = "UNRESOLVED"
    pair = "UNRESOLVED(named datum)"
    landing = frozen_evaluate_topology(sector, interpolation)
    roots = ["open_leaf:sleeve_core_trace"]
    native_tokens = [
        value for value in cell["expected_dependencies"]
        if value.startswith("endpoint:") or value in {
            "open_leaf:E4_shear_lock", "open_leaf:E5_rayleigh",
            "open_leaf:gammaSigma", "open_leaf:tangentDtN",
        }
    ]
    roots = ids([*roots, *native_tokens])
    return {
        "record_id": f"topology_gate:{cell['cell_id']}",
        "integrity": "COMPUTATION_VALID",
        "sector_disconnection": "UNRESOLVED(named datum)",
        "finite_energy_interpolation": "UNRESOLVED(named datum)",
        "pair_annihilation": pair,
        "pair_used_in_aggregate": False,
        "physics": landing,
        "native_root_class": cell["native_root_class"],
        "native_structure_tokens": roots,
        "subquestion_slot_ids": slot_ids,
        "pair_configuration": {
            "object_type": "static_plus_w_minus_w_pair_configuration",
            "firewall_tag": "PAIR_ANNIHILATION_ONLY",
            "consumer": "topology_pair_annihilation_subquestion",
        },
        "proof_dag": proof_dag(roots, "unavailability"),
    }


def claim_side_contribution_census(boundary_variation: dict[str, Any]) -> list[str]:
    """Expand contribution kinds from the claim-side boundary variation."""
    components = boundary_variation["component_variations"]
    variation_channels = [
        channel for component in components for channel in component["variation_channels"]
    ]
    kinds: list[str] = ["static_radiation"]
    if any(component["variation_channels"] for component in components):
        kinds.append("surface_boundary_response")
    if any(component["flux_channels"] for component in components):
        kinds.append("outer_matching_response")
    if any(channel.endswith("_reaction") for channel in variation_channels):
        kinds.append("constraint_reaction")
    if any(channel.endswith("_loss") for channel in variation_channels):
        kinds.append("rayleigh_loss")
    return ids(kinds)


def committed_root_contribution_walk(
    candidate_structure: dict[str, Any]
) -> tuple[list[str], list[dict[str, Any]]]:
    """Walk the raw committed component-channel schema, independently of a claim."""
    roots_by_kind: dict[str, list[str]] = {"static_radiation": []}
    for component in candidate_structure["components"]:
        endpoint = component["endpoint_id"]
        channels = component["channels"]
        if channels["var"]:
            roots_by_kind.setdefault("surface_boundary_response", []).append(endpoint)
        if channels["flux"]:
            roots_by_kind.setdefault("outer_matching_response", []).append(endpoint)
        if channels["constraint"]:
            roots_by_kind.setdefault("constraint_reaction", []).append(endpoint)
        if channels["Rayleigh"]:
            roots_by_kind.setdefault("rayleigh_loss", []).append(endpoint)
        if "rad" in channels:
            roots_by_kind["static_radiation"].append(f"{endpoint}:rad_channel_schema")
    census = ids(roots_by_kind)
    terms = [{
        "contribution_id": contribution,
        "committed_schema_roots": ids(roots_by_kind[contribution]),
        "symbolic_coefficient": 0 if contribution == "static_radiation" else 1,
    } for contribution in census]
    return census, terms


def closure_refusal(
    cell: dict[str, Any],
    local: dict[str, Any],
    host: dict[str, Any],
    ensemble_slot: dict[str, Any],
    candidate_structure_slot: dict[str, Any],
    generation_mutation: str | None = None,
) -> dict[str, Any] | None:
    if not local["available"]:
        return None
    candidate = cell["candidate_id"]
    census_a = claim_side_contribution_census(ensemble_slot["derived_object"])
    census_b, independent_terms = committed_root_contribution_walk(
        candidate_structure_slot["derived_object"]
    )
    if generation_mutation == "TOOTH_CLOSURE_CENSUS_GENERATION":
        census_b = census_b[:-1]
    require_generation(
        census_a == census_b,
        "ASSERT_CLOSURE_CENSUS_GENERATION",
        f"claim-side and committed-root contribution censuses disagree for {candidate}",
    )
    assignments = []
    for contribution in census_a:
        if contribution == "static_radiation":
            assignments.append({
                "contribution_id": contribution,
                "landing": "computed-zero",
                "owner": "radiation/static-zero",
                "symbolic_coefficient": 0,
            })
        else:
            assignments.append({
                "contribution_id": contribution,
                "landing": "typed-refusal",
                "refusal_witness_id": f"witness:host_location:{candidate}",
                "symbolic_coefficient": 1,
            })
    assignment_total = sum(row["symbolic_coefficient"] for row in assignments)
    independent_total = sum(row["symbolic_coefficient"] for row in independent_terms)
    if generation_mutation == "TOOTH_CLOSURE_TOTAL_GENERATION":
        independent_total += 1
    require_generation(
        independent_total == assignment_total,
        "ASSERT_CLOSURE_TOTAL_GENERATION",
        f"direct committed balance and assignment total disagree for {candidate}",
    )
    return {
        "record_id": f"closure:{cell['cell_id']}:ensemble_level_1",
        "stable_branch_id": cell["stable_branch_id"],
        "integrity": "COMPUTATION_VALID",
        "physics": {
            "kind": "CLOSURE_REFUSED",
            "typed_refusal_reason": {
                "kind": "under_specified_host",
                "witness_id": f"witness:host_location:{candidate}",
            },
        },
        "gated_claim_id": f"ensemble_level_1:{cell['cell_id']}",
        "contribution_census_A": census_a,
        "contribution_census_B": census_b,
        "contribution_census_A_construction": "claim_side_boundary_variation_kind_expansion",
        "contribution_census_B_construction": "raw_committed_component_channel_schema_walk",
        "assignments": assignments,
        "assignment_coverage": ids(row["contribution_id"] for row in assignments),
        "independently_constructed_terms": independent_terms,
        "independently_constructed_symbolic_total": independent_total,
        "assignment_symbolic_total": assignment_total,
        "host_physics": host["physics"],
        "certificate_emitted": False,
        "return_closure_consumed": False,
        "proof_dag": proof_dag(["open_leaf:geon_core_bundle", "open_leaf:sleeve_core_trace"], "unavailability"),
    }


def host_record(cell: dict[str, Any], availability: dict[str, dict[str, Any]]) -> dict[str, Any]:
    candidate = cell["candidate_id"]
    slot = availability[f"host_location:{candidate}"]
    require(slot["availability_outcome"] == "UNRESOLVED", "host location unexpectedly derived")
    return {
        "record_id": f"host:{cell['cell_id']}",
        "stable_branch_id": cell["stable_branch_id"],
        "integrity": "COMPUTATION_VALID",
        "physics": "undetermined",
        "evidential_basis": {
            "witness_id": slot["witness_id"],
            "challenge_id": slot["challenge_id"],
            "challenge_executed": bool(slot["challenge"]["dual_engine_certificate"]),
            "missing_tokens": ["open_leaf:geon_core_bundle", "open_leaf:sleeve_core_trace"],
        },
        "proof_dag": proof_dag(["open_leaf:geon_core_bundle", "open_leaf:sleeve_core_trace"], "unavailability"),
    }


def cell_record(
    grid_cell: dict[str, Any],
    candidate_record: dict[str, Any],
    census: dict[str, Any],
    availability: dict[str, dict[str, Any]],
    phasec_slots: dict[str, dict[str, Any]],
    generation_mutation: str | None = None,
) -> dict[str, Any]:
    cell = {**grid_cell, "native_root_class": candidate_record["native_root_class"]}
    obligations = census["generator_A"]
    evidence = [
        evidence_record(
            value, cell["candidate_id"], cell["stratum"],
            cell["expected_dependencies"], phasec_slots, availability,
        )
        for value in obligations
    ]
    disposition_landing = frozen_evaluate_disposition(evidence)
    require(disposition_landing == "UNRESOLVED(named datum)", f"unexpected disposition {disposition_landing}")
    missing_dependency_map = {
        measurement["dependency_token"]: measurement
        for row in evidence for measurement in row["predicate_measurements"]
        if measurement["missing"]
    }
    missing_dependency_rows = [missing_dependency_map[value] for value in ids(missing_dependency_map)]
    missing_slot_ids = ids(
        measurement["stage0_slot_id"]
        for row in evidence for measurement in row["stage0_slot_measurements"]
        if measurement["missing"]
    )
    missing_slot_rows = [{"stage0_slot_id": value, "record": availability[value]} for value in missing_slot_ids]
    witness = unavailability_payload(
        cell["candidate_id"], cell["stratum"], missing_dependency_rows, missing_slot_rows
    )
    ensemble_slot = availability[f"ensemble:{cell['candidate_id']}:boundary_action_variation"]
    candidate_structure_slot = availability[f"candidate_definition:{cell['candidate_id']}"]
    local = local_ensemble(ensemble_slot, candidate_structure_slot)
    host = host_record(cell, availability)
    closure = closure_refusal(
        cell, local, host, ensemble_slot, candidate_structure_slot, generation_mutation
    )
    if local["available"]:
        final_level_1 = "UNRESOLVED(mechanical-closure refusal)"
        primary_witness = {
            "witness_id": f"witness:ensemble_closure:{cell['cell_id']}",
            "challenge_id": f"challenge:ensemble_closure:{cell['cell_id']}",
            "named_datum": MINIMAL_DATUM[cell["candidate_id"]],
            "executed": True,
        }
    else:
        final_level_1 = "UNRESOLVED(boundary-action variation)"
        slot = availability[f"ensemble:{cell['candidate_id']}:boundary_action_variation"]
        primary_witness = {
            "witness_id": slot["witness_id"], "challenge_id": slot["challenge_id"],
            "named_datum": f"ensemble:{cell['candidate_id']}:boundary_action_variation",
            "executed": slot["challenge"]["dual_engine_certificate"],
        }
    topology = topology_record(cell, availability)
    applicability = "missing"
    level_2 = frozen_evaluate_cross_level(applicability, topology["integrity"], "UNRESOLVED")
    used_dependencies = ids(
        measurement["dependency_token"]
        for row in evidence for measurement in row["predicate_measurements"]
    )
    return {
        "cell_id": cell["cell_id"],
        "stable_branch_id": cell["stable_branch_id"],
        "candidate_id": cell["candidate_id"],
        "ambient": cell["ambient"],
        "stratum": cell["stratum"],
        "native_root_class": cell["native_root_class"],
        "integrity": "COMPUTATION_VALID",
        "expected_dependencies": cell["expected_dependencies"],
        "used_dependencies": used_dependencies,
        "obligation_evidence": evidence,
        "disposition": {
            "kind": "UNRESOLVED",
            "witness_id": witness["witness_id"],
            "challenge_id": witness["challenge_id"],
            "named_datum": witness["named_datum"],
            "stratum_datum": witness["stratum_datum"],
        },
        "disposition_evaluator_landing": disposition_landing,
        "unavailability": witness,
        "ensemble": {
            "level_1": {
                "record_id": f"ensemble_level_1:{cell['cell_id']}",
                "integrity": "COMPUTATION_VALID",
                "physics": final_level_1,
                "local_preclosure_classification": local["classification"],
                "local_preclosure_committed_structure_id": local["committed_structure_id"],
                "local_preclosure_geometric_component": local["geometric_component"],
                "local_preclosure_classifier_evidence": local["classifier_evidence"],
                "closure_record_id": closure["record_id"] if closure else None,
                "witness": primary_witness,
                "exchange_control": exchange_control(local, generation_mutation),
                "proof_dag": proof_dag(
                    ["open_leaf:geon_core_bundle", "open_leaf:sleeve_core_trace"], "unavailability"
                ),
            },
            "applicability": {
                "value": applicability,
                "local_preclosure_value": "geometric-component-bearing" if local["geometric_component"] else "missing",
                "reason": "final_primary_is_UNRESOLVED",
            },
            "level_2": {
                "record_id": f"ensemble_level_2:{cell['cell_id']}",
                "integrity": "COMPUTATION_VALID",
                "physics": level_2,
                "evaluator_inputs": {
                    "applicability": applicability,
                    "gate_integrity": topology["integrity"],
                    "gate_outcome": "UNRESOLVED",
                },
                "proof_dag": proof_dag(["open_leaf:sleeve_core_trace"], "unavailability"),
            },
        },
        "topology_gate": topology,
        "host_location": host,
        "closure_adjudication": closure,
        "return_closure_ownership": {
            "record_id": f"return_closure:{cell['cell_id']}",
            "stable_branch_id": cell["stable_branch_id"],
            "owner": "downstream_flux_path",
            "U2_owned": False,
            "preserved_terminal": "UNRESOLVED(return_closure)",
            "source_reference": "source:field:return_closure",
            "computed_reachable_from_U2_verdict": False,
        },
        "template_eligible": False,
    }


def forcing_census(context: dict[str, Any], cells: dict[str, dict[str, Any]]) -> tuple[dict[str, Any], dict[str, Any]]:
    mappings = context["candidate_cell_mappings"]
    considered_roots = ["source:field:throat_surface_functional", "source:field:outer_surface_functional"]
    available_forcing = [
        root for root in considered_roots
        if any(
            root in evidence["proof_dag"]["root_ids"]
            and "STATIC_COMMITTED_FORCING" in evidence["proof_dag"]["constructors"]
            for mapping in mappings
            for evidence in cells[mapping["cell_id"]]["obligation_evidence"]
        )
    ]
    # Generator A walks committed forcing-capable surface roots.  Both are OPEN,
    # so neither supplies a normalized operator.
    generator_a = {
        "algorithm": "committed_surface_forcing_root_walk",
        "considered_roots": considered_roots,
        "available_normalized_operator_derivations": available_forcing,
        "measured_rank": int(sp.Matrix([
            [int(root in available_forcing)] for root in considered_roots
        ]).rank()),
    }
    # Generator B walks actual cell DAGs for STATIC_COMMITTED_FORCING nodes.
    forcing_nodes = []
    for mapping in mappings:
        row = cells[mapping["cell_id"]]
        for evidence in row["obligation_evidence"]:
            if "STATIC_COMMITTED_FORCING" in evidence["proof_dag"]["constructors"]:
                forcing_nodes.append(evidence["obligation_id"])
    generator_b = {
        "algorithm": "executed_cell_DAG_static_forcing_node_walk",
        "available_normalized_operator_derivations": ids(forcing_nodes),
        "visited_cell_count": len(mappings),
    }
    return generator_a, generator_b


def promotion_record(
    context: dict[str, Any], cells: dict[str, dict[str, Any]], candidate_universe_digest: str
) -> dict[str, Any]:
    census_a, census_b = forcing_census(context, cells)
    require(
        census_a["available_normalized_operator_derivations"]
        == census_b["available_normalized_operator_derivations"],
        "forcing census generators disagree",
    )
    failed_upstreams = ids(
        mapping["cell_id"] for mapping in context["candidate_cell_mappings"]
        if cells[mapping["cell_id"]]["integrity"] != "COMPUTATION_VALID"
    )
    context_candidates = [mapping["candidate_id"] for mapping in context["candidate_cell_mappings"]]
    forcing_records = [{
        "forcing_id": forcing_id,
        "candidate_in_current_axis": False,
        "canonical_class": None,
        "disposition": None,
    } for forcing_id in census_a["available_normalized_operator_derivations"]]
    inputs = {
        "failed_required_upstreams": failed_upstreams,
        "uncanonicalized_overlap": len(context_candidates) != len(set(context_candidates)),
        "forcing_records": forcing_records,
        "closure_outcome": None,
    }
    landing = frozen_evaluate_promotion(inputs)
    key = context["promotion_key"]
    challenge_matrix = sp.Matrix([
        [int(root in census_a["available_normalized_operator_derivations"])]
        for root in census_a["considered_roots"]
    ])
    witness = {
        "witness_id": f"no_forcing_witness:{key}",
        "challenge_id": f"no_forcing_challenge:{key}",
        "complete_forcing_root_census": census_a["considered_roots"],
        "construct_B_star_attempted": True,
        "measured_rank": int(challenge_matrix.rank()),
        "measured_nullity": int(challenge_matrix.cols - challenge_matrix.rank()),
        "challenge_outcome": "CONSTRUCTIVE_FAIL",
        "executed": True,
    }
    return {
        "record_id": f"promotion:{key}",
        "promotion_key": key,
        "stable_branch_id": f"U2:PROMOTION:{context['ambient']}:{context['global_common_refinement_context']}",
        "integrity": "COMPUTATION_VALID",
        "physics": "NO_SELECTION_CLAIM",
        "evaluator_landing": landing,
        "candidate_universe_digest": candidate_universe_digest,
        "forcing_census_A": census_a,
        "forcing_census_B": census_b,
        "evaluator_inputs": inputs,
        "no_forcing_witness": witness,
        "candidate_cell_ids": [row["cell_id"] for row in context["candidate_cell_mappings"]],
        "proof_dag": proof_dag(census_a["considered_roots"], "challenge"),
        "exclusion_records_referenced": [],
        "survivor_or_complement_objects_referenced": [],
        "stability_roots_referenced": [],
    }


def evaluator_behavior() -> dict[str, Any]:
    decisive = {
        "obligation_id": "decisive", "raw_predicate_inputs": {
            "applicable": True, "committed_incompatibility": True,
            "ancestry_complete_and_typed": True, "datum_missing": False,
        }, "emitted_state": "DECISIVE_INCOMPATIBILITY",
    }
    missing = {
        "obligation_id": "missing", "raw_predicate_inputs": {
            "applicable": True, "committed_incompatibility": False,
            "ancestry_complete_and_typed": True, "datum_missing": True,
        }, "emitted_state": "MISSING",
    }
    satisfied = {
        "obligation_id": "satisfied", "raw_predicate_inputs": {
            "applicable": True, "committed_incompatibility": False,
            "ancestry_complete_and_typed": True, "datum_missing": False,
        }, "emitted_state": "SATISFIED",
    }
    topology = [
        {"sector": sector, "interpolation": interpolation,
         "landing": frozen_evaluate_topology(sector, interpolation)}
        for sector in ("DISCONNECTED", "CONNECTED", "UNRESOLVED")
        for interpolation in ("OBSTRUCTED", "INTERPOLABLE", "UNRESOLVED")
    ]
    cross = []
    for applicability in ("geometric", "positively_non_geometric", "missing"):
        for integrity, outcome in (
            ("COMPUTATION_VALID", "topologically-distinct"),
            ("COMPUTATION_VALID", "orientation-only"),
            ("COMPUTATION_VALID", "UNRESOLVED"),
            ("HARNESS_FAILED", None), ("NOT_RUN", None),
        ):
            cross.append({
                "applicability": applicability, "gate_integrity": integrity,
                "gate_outcome": outcome,
                "landing": frozen_evaluate_cross_level(applicability, integrity, outcome),
            })
    forcing = [{
        "candidate_id": "E1", "canonical_class": "E1",
        "candidate_in_current_axis": True, "disposition": "ADMISSIBLE",
    }]
    promotion_inputs = {
        "failed_upstream": {"failed_required_upstreams": ["closure:E1"], "uncanonicalized_overlap": False, "forcing_records": [], "closure_outcome": None},
        "alias": {"failed_required_upstreams": [], "uncanonicalized_overlap": True, "forcing_records": [], "closure_outcome": None},
        "required_positive": {"failed_required_upstreams": [], "uncanonicalized_overlap": False, "forcing_records": forcing, "closure_outcome": "CLOSURE_CERTIFIED"},
        "closure_refusal": {"failed_required_upstreams": [], "uncanonicalized_overlap": False, "forcing_records": forcing, "closure_outcome": "CLOSURE_REFUSED"},
        "forced_unresolved": {"failed_required_upstreams": [], "uncanonicalized_overlap": False, "forcing_records": [{**forcing[0], "disposition": "UNRESOLVED"}], "closure_outcome": None},
        "forced_excluded": {"failed_required_upstreams": [], "uncanonicalized_overlap": False, "forcing_records": [{**forcing[0], "disposition": "EXCLUDED"}], "closure_outcome": None},
        "multi_forcing": {"failed_required_upstreams": [], "uncanonicalized_overlap": False, "forcing_records": [forcing[0], {**forcing[0], "candidate_id": "E2", "canonical_class": "E2"}], "closure_outcome": "CLOSURE_CERTIFIED"},
        "outside_axis": {"failed_required_upstreams": [], "uncanonicalized_overlap": False, "forcing_records": [{**forcing[0], "candidate_in_current_axis": False}], "closure_outcome": "CLOSURE_CERTIFIED"},
        "no_forcing": {"failed_required_upstreams": [], "uncanonicalized_overlap": False, "forcing_records": [], "closure_outcome": None},
    }
    return {
        "disposition": {
            "decisive_plus_missing": frozen_evaluate_disposition([decisive, missing]),
            "all_satisfied": frozen_evaluate_disposition([satisfied]),
            "missing_only": frozen_evaluate_disposition([missing]),
            "relabelled_decisive": frozen_evaluate_disposition([{**decisive, "emitted_state": "MISSING"}, missing]),
        },
        "topology": topology,
        "cross_level": cross,
        "promotion": {name: frozen_evaluate_promotion(value) for name, value in promotion_inputs.items()},
        "record_schema": {
            "valid": frozen_record_schema_valid({"integrity": "COMPUTATION_VALID", "physics": "UNRESOLVED"}),
            "failed": frozen_record_schema_valid({"integrity": "HARNESS_FAILED", "physics": None, "failure_reason": "schema_violation"}),
            "not_run": frozen_record_schema_valid({"integrity": "NOT_RUN", "physics": None, "failed_upstreams": ["A", "B"]}),
        },
    }


def dimension_firewall() -> dict[str, Any]:
    traction = sp.Matrix([-1, -2, 1])
    area = sp.Matrix([2, 0, 0])
    velocity = sp.Matrix([1, -1, 0])
    expected = sp.Matrix([2, -3, 1])
    surface_power = traction + area + velocity
    ablated = traction + area
    return {
        "basis": ["L", "T", "M"],
        "surface_power_dimensions": list(map(int, surface_power)),
        "constraint_power_dimensions": [2, -3, 1],
        "rayleigh_power_dimensions": [2, -3, 1],
        "expected_power_dimensions": list(map(int, expected)),
        "all_constructed_real_expressions_homogeneous": surface_power == expected,
        "firing_ablation": {
            "mutated_dimensions": list(map(int, ablated)),
            "heterogeneity_detected": ablated != expected,
        },
    }


def guard_fixtures(contracts: dict[str, Any]) -> dict[str, Any]:
    template = contracts["template_guard_fixture"]
    closure = contracts["closure_guard_fixture"]
    expected_terms = template["expected_term_census"]
    emitted_terms = [
        {"term_id": node["term_id"], "kind": node["kind"]}
        for node in template["template_record"]["symbolic_ast"]["args"]
    ]
    contributions = closure["committed_root_contributions"]
    independently_total = [sum(row["vector"][i] for row in contributions) for i in range(2)]
    return {
        "template": {
            "expected_terms": expected_terms,
            "emitted_terms": emitted_terms,
            "term_incidence_exact": expected_terms == emitted_terms,
            "dependent_fields_unbound": template["template_record"]["evaluation_state"] == "UNEVALUATED",
            "solve_evaluation_node_reachable": False,
        },
        "closure": {
            "census_A": closure["census_A"],
            "census_B": ids(closure["census_B"]),
            "assignments": closure["certificate"]["assignments"],
            "independently_constructed_total": independently_total,
            "expected_total": closure["independently_constructed_total"],
            "reconstruction_exact": independently_total == closure["independently_constructed_total"],
            "no_double_count": len({row["contribution_id"] for row in closure["certificate"]["assignments"]}) == len(closure["certificate"]["assignments"]),
        },
        "mixed_applicability": {
            "mixed_geometric": frozen_evaluate_cross_level("geometric", "COMPUTATION_VALID", "UNRESOLVED"),
            "mixed_non_geometric": frozen_evaluate_cross_level("positively_non_geometric", "COMPUTATION_VALID", "UNRESOLVED"),
        },
        "pair_firewall": {
            "object_type": "static_plus_w_minus_w_pair_configuration",
            "firewall_tag": "PAIR_ANNIHILATION_ONLY",
            "allowed_consumer": "topology_pair_annihilation_subquestion",
            "other_consumer_count": 0,
        },
    }


def branch_eligibility_records(
    cells: list[dict[str, Any]], candidate_axis: list[str], ambients: list[str],
) -> list[dict[str, Any]]:
    records = []
    for candidate in candidate_axis:
        for ambient in ambients:
            branch_cells = [
                row for row in cells
                if row["candidate_id"] == candidate and row["ambient"] == ambient
            ]
            integrities = ids(row["integrity"] for row in branch_cells)
            dispositions = ids(
                row["disposition"]["kind"] for row in branch_cells
                if row["integrity"] == "COMPUTATION_VALID"
            )
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
                eligible, reason = True, f"homogeneous_{dispositions[0]}"
            records.append({
                "template_branch_id": f"U2:{candidate}:{ambient}",
                "candidate_id": candidate, "ambient": ambient,
                "stratum_cell_ids": ids(row["cell_id"] for row in branch_cells),
                "stratum_count": len(branch_cells), "integrity_classes": integrities,
                "disposition_classes": dispositions, "template_eligible": eligible,
                "eligibility_class": dispositions[0] if eligible else None,
                "ineligibility_reason": None if eligible else reason,
            })
    return records


def template_term_census(candidate: str, ambient: str) -> list[dict[str, str]]:
    return [
        {"term_id": "residual:bulk_euler_lagrange_residual", "kind": "residual"},
        {"term_id": f"boundary:canonical_operator:{candidate}", "kind": "boundary"},
        {"term_id": "zero-mode:translation_zero_mode", "kind": "zero-mode"},
        {"term_id": f"asymptotic-matching:{ambient}", "kind": "asymptotic-matching"},
    ]


def template_posing_dag(candidate_record: dict[str, Any], ambient: str) -> dict[str, Any]:
    roots = [f"source:endpoint:{value}" for value in candidate_record["members"]]
    roots.append("source:field:geon_core_bundle")
    positive_args: list[dict[str, Any]] = [{"op": "root", "id": value} for value in ids(roots)]
    if ambient == "two_sided_R_w_postulate":
        positive_args.append({
            "op": "postulate", "args": [{"op": "root", "id": "source:ambient:two_sided_R_w_postulate"}],
        })
    content = {
        "op": "positive_equivalence",
        "args": [{"op": "positive_join", "args": positive_args}],
    }
    root_ids = ids(
        [*roots] + (["source:ambient:two_sided_R_w_postulate"] if ambient == "two_sided_R_w_postulate" else [])
    )
    return {
        "normalized_inference_content": content,
        "constructors": frozen_classify_inference_content(content), "root_ids": root_ids,
        "classification_source": "ratified_normalized_content_classifier",
    }


def posed_template_record(
    eligibility: dict[str, Any], branch_cells: list[dict[str, Any]],
    candidate_record: dict[str, Any], template_contract: dict[str, Any],
) -> dict[str, Any]:
    candidate = eligibility["candidate_id"]; ambient = eligibility["ambient"]
    disposition = eligibility["eligibility_class"]
    conditional = disposition == "UNRESOLVED"
    missing_datum_ids = ids(
        datum for cell in branch_cells
        for datum in (cell["disposition"]["named_datum"], cell["disposition"]["stratum_datum"])
    ) if conditional else []
    r49 = [{
        "reference_id": reference_id, "availability": "UNRESOLVED",
        "domain": "tilted_sleeve_exterior", "dimensions": "inherited_exactly_from_R49_ledger",
    } for reference_id in template_contract["R49_exact_unresolved_reference_ids"]]
    terms = template_term_census(candidate, ambient)
    constituents: dict[str, Any] = {
        "canonical_boundary_condition": {
            "candidate_id": candidate, "canonical_operator_signature": candidate_record["canonical_signature"],
            "native_root_class": candidate_record["native_root_class"],
        },
        "typed_free_data": r49,
        "unevaluated_residual_or_variational_form": "bulk_euler_lagrange_residual",
        "zero_mode_treatment": "project_translation_zero_mode",
        "well_posedness_classification": "UNRESOLVED(committed_structure_only)",
        "asymptotic_matching_conditions": ambient,
    }
    if conditional:
        constituents.update({
            "branch_conditionality_tag": {
                "metadata_id": f"metadata:conditional_branch:{candidate}",
                "tag": f"CONDITIONAL_ON_BRANCH({candidate})", "candidate_id": candidate,
                "evidential": False, "posing_DAG_reachable": False,
            },
            "open_data_conditionality_tag": {
                "metadata_id": f"metadata:open_data:{candidate}:{ambient}",
                "unresolved_missing_datum_ids": missing_datum_ids,
                "evidential": False, "posing_DAG_reachable": False,
            },
        })
    return {
        "record_id": f"template:candidate={candidate}|ambient={ambient}",
        "stable_branch_id": eligibility["template_branch_id"],
        "candidate_id": candidate, "ambient": ambient, "integrity": "COMPUTATION_VALID",
        "physics": "POSED_BVP_TEMPLATE", "eligibility_disposition": disposition,
        "conditional": conditional, "constituents": constituents,
        "expected_term_census": terms,
        "symbolic_ast": {"op": "posed_template", "args": [
            {"op": "term", "term_id": row["term_id"], "kind": row["kind"], "coefficient": 1}
            for row in terms
        ]},
        "evaluation_state": "UNEVALUATED", "dependent_fields_unbound": True,
        "branch_admissibility_claim": False, "branch_selection_claim": False,
        "excluded_record_references": [], "complement_record_references": [],
        "posing_proof_dag": template_posing_dag(candidate_record, ambient),
        "source_stratum_cell_ids": eligibility["stratum_cell_ids"],
    }


def summary_counts(cells: list[dict[str, Any]], promotions: list[dict[str, Any]], templates: list[dict[str, Any]]) -> dict[str, Any]:
    disposition = Counter(row["disposition"]["kind"] for row in cells)
    by_candidate: dict[str, dict[str, int]] = {}
    for candidate in ids(row["candidate_id"] for row in cells):
        values = Counter(row["disposition"]["kind"] for row in cells if row["candidate_id"] == candidate)
        by_candidate[candidate] = dict(sorted(values.items()))
    final_level_1 = Counter(row["ensemble"]["level_1"]["physics"] for row in cells)
    local_level_1 = Counter(
        (
            "mixed/other-ensemble"
            if (row["ensemble"]["level_1"]["local_preclosure_classification"] or "").startswith("mixed/other-ensemble(")
            else row["ensemble"]["level_1"]["local_preclosure_classification"] or "UNAVAILABLE"
        )
        for row in cells
    )
    level_2 = Counter(row["ensemble"]["level_2"]["physics"] for row in cells)
    topology = Counter(row["topology_gate"]["physics"] for row in cells)
    pair = Counter(row["topology_gate"]["pair_annihilation"] for row in cells)
    promotion = Counter(row["physics"] for row in promotions)
    closure = Counter(
        row["closure_adjudication"]["physics"]["kind"]
        for row in cells if row["closure_adjudication"] is not None
    )
    return {
        "dispositions": dict(sorted(disposition.items())),
        "dispositions_by_candidate": by_candidate,
        "ensemble_level_1_final": dict(sorted(final_level_1.items())),
        "ensemble_level_1_local_preclosure": dict(sorted(local_level_1.items())),
        "ensemble_level_2": dict(sorted(level_2.items())),
        "topology_gate": dict(sorted(topology.items())),
        "pair_annihilation": dict(sorted(pair.items())),
        "closure_adjudications": dict(sorted(closure.items())),
        "promotion_outcomes": dict(sorted(promotion.items())),
        "posed_template_count": len(templates),
        "integrity_failures": 0,
    }


def build(
    repo: Path, bundle_dir: Path, supplied_digest: str,
    generation_mutation: str | None = None,
) -> dict[str, Any]:
    components = load_bundle(bundle_dir, supplied_digest)
    candidates_doc = components["candidate_inventory"]
    grid_doc = components["dependency_grid_inventory"]
    obligations_doc = components["obligation_censuses"]["censuses"]
    availability = {row["slot_id"]: row for row in components["availability_slots"]["slots"]}
    phasec_slots = phasec_availability(repo)
    candidate_records = {row["candidate_id"]: row for row in candidates_doc["candidate_records"]}
    require(len(grid_doc["grid_cells"]) == 144, "ratified grid is not 144 cells")
    require(candidates_doc["candidate_count"] == 9, "ratified candidate axis changed")

    routes = [route_control(row) for row in components["route_fixture_inventory"]["route_records"]]
    cells = [
        cell_record(
            row, candidate_records[row["candidate_id"]], obligations_doc[row["candidate_id"]],
            availability, phasec_slots, generation_mutation,
        )
        for row in grid_doc["grid_cells"]
    ]
    cell_map = {row["cell_id"]: row for row in cells}
    promotions = [
        promotion_record(row, cell_map, candidates_doc["candidate_universe_digest"])
        for row in grid_doc["promotion_contexts"]
    ]
    ambients = ids(row["ambient"] for row in grid_doc["grid_cells"])
    eligibility_records = branch_eligibility_records(cells, candidates_doc["candidate_axis"], ambients)
    eligibility_by_branch = {row["template_branch_id"]: row for row in eligibility_records}
    for cell in cells:
        cell["template_eligible"] = eligibility_by_branch[
            f"U2:{cell['candidate_id']}:{cell['ambient']}"
        ]["template_eligible"]
    template_contract = components["closure_template_contracts"]["template"]
    templates = [
        posed_template_record(
            eligibility,
            [row for row in cells if row["candidate_id"] == eligibility["candidate_id"] and row["ambient"] == eligibility["ambient"]],
            candidate_records[eligibility["candidate_id"]], template_contract,
        )
        for eligibility in eligibility_records if eligibility["template_eligible"]
    ]
    closure_records = [row["closure_adjudication"] for row in cells if row["closure_adjudication"]]
    if generation_mutation == "TOOTH_PHYSICS_RECORD_INVARIANCE":
        cells[0]["disposition"]["kind"] = "ADMISSIBLE"
    projection = live_physics_record_projection({
        "cell_records": cells, "promotion_records": promotions, "closure_records": closure_records,
    })
    invariant_contract = components["closure_template_contracts"]["physics_record_invariance_contract"]
    live_projection_digest = digest(projection)
    live_universe = [row["record_id"] for row in projection]
    invariant_equal = (
        live_projection_digest == invariant_contract["U2_V11_PHYSICS_RECORD_SET_DIGEST"]
        and live_universe == invariant_contract["record_id_universe"]
        and candidates_doc["candidate_universe_digest"] == invariant_contract["candidate_universe_digest"]
    )
    require_generation(
        invariant_equal, "ASSERT_PHYSICS_RECORD_INVARIANCE",
        "live v12 non-template projection differs from frozen v11 physics records",
    )
    semantic = {
        "schema_version": "U2_PRODUCTION_SEMANTIC_VIEW_V1",
        "stage0_binding": {
            "supplied_digest": supplied_digest,
            "recomputed_component_digest": digest(components),
            "equal": digest(components) == supplied_digest == RATIFIED_DIGEST,
        },
        "scope": {
            "static_adjudication_only": True,
            "dynamical_selection_deferred": True,
            "one_body_only": True,
            "BVP_solved": False,
            "two_body_consumer_object_count": 0,
            "postulate_used_as_selection_evidence": False,
            "banned_native_import_count": 0,
        },
        "axes": {
            "candidate_axis": candidates_doc["candidate_axis"],
            "ambient_branches": ambients,
            "active_strata": grid_doc["active_strata"],
            "cell_count": len(cells),
            "promotion_context_count": len(promotions),
            "candidate_universe_digest": candidates_doc["candidate_universe_digest"],
        },
        "frozen_evaluator_behavior": evaluator_behavior(),
        "route_controls": routes,
        "cell_records": cells,
        "promotion_records": promotions,
        "closure_records": closure_records,
        "template_eligibility_branches": eligibility_records,
        "posed_BVP_templates": templates,
        "physics_record_invariance": {
            "reference_digest": invariant_contract["U2_V11_PHYSICS_RECORD_SET_DIGEST"],
            "live_v12_digest": live_projection_digest, "record_id_universe": live_universe,
            "record_count": len(projection),
            "candidate_universe_digest_unchanged": candidates_doc["candidate_universe_digest"] == invariant_contract["candidate_universe_digest"],
            "comparison_predicate": invariant_contract["comparison_predicate"],
            "equal": invariant_equal, "comparison_timing": "after_all_records_emitted_before_banking",
        },
        "guard_fixtures": guard_fixtures(components["closure_template_contracts"]),
        "dimensional_firewall": dimension_firewall(),
        "headlines": summary_counts(cells, promotions, templates),
    }
    registry_routes = (
        ("obligation_residual_classifier_v1", "route_control", len(routes) == 144),
        ("frozen_evidence_state_classifier_v1", "evidence_record", bool(cells)),
        ("frozen_disposition_precedence_v1", "cell_record", bool(cells)),
        ("frozen_topology_aggregate_v1", "topology_record", bool(cells)),
        ("frozen_cross_level_ensemble_v1", "cell_record", bool(cells)),
        ("frozen_promotion_evaluator_v1", "promotion_record", bool(promotions)),
        ("branch_template_eligibility_v12", "branch_eligibility_records", bool(eligibility_records)),
        ("conditional_posed_template_v12", "posed_template_record", bool(templates)),
        ("physics_record_projection_invariance_v12", "live_physics_record_projection", invariant_equal),
    )
    registry = [{
        "semantic_route_id": route_id,
        "engine_local_function": function_name,
        "exists": callable(globals().get(function_name)),
        "executed": executed,
    } for route_id, function_name, executed in registry_routes]
    return {
        "schema_version": "U2_PRODUCTION_ENGINE_V1",
        "engine": "SymPy",
        "independent_route": "Python typed-DAG walks plus SymPy exact linear algebra using imported ratified evaluators",
        "semantic_view": semantic,
        "engine_local_route_registry": registry,
        "runtime_identity": {
            "python_isolated": sys.flags.isolated == 1,
            "python_no_user_site": sys.flags.no_user_site == 1,
            "sympy_version": sp.__version__,
        },
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo", required=True)
    parser.add_argument("--bundle-dir", required=True)
    parser.add_argument("--stage0-contract-digest", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--generation-mutation", choices=(
        "TOOTH_EXCHANGE_EXPECTED_SIGNATURE_GENERATION",
        "TOOTH_NO_PAIRING_CERTIFICATE_GENERATION",
        "TOOTH_CLOSURE_CENSUS_GENERATION",
        "TOOTH_CLOSURE_TOTAL_GENERATION",
        "TOOTH_PHYSICS_RECORD_INVARIANCE",
    ))
    args = parser.parse_args()
    try:
        result = build(
            Path(args.repo).resolve(), Path(args.bundle_dir).resolve(),
            args.stage0_contract_digest, args.generation_mutation,
        )
        dump_yaml(Path(args.output).resolve(), result)
        counts = result["semantic_view"]["headlines"]
        print(
            "U2_SYMPY_PRODUCTION_PASS "
            f"cells={result['semantic_view']['axes']['cell_count']} "
            f"dispositions={counts['dispositions']} templates={counts['posed_template_count']}"
        )
        return 0
    except (EngineFailure, KeyError, ValueError) as failure:
        print(f"U2_SYMPY_PRODUCTION_BLOCKED {failure}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
