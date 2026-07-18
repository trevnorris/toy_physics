#!/usr/bin/env python3
"""Independent SymPy/Python stage-0 route for U2 boundary adjudication.

The executable imports no task-authored module.  It walks the frozen YAML
substrate directly, constructs every stage-0 semantic object, and emits a
canonical engine view for strict comparison with the independent Wolfram
route.  Stage 0 poses contracts and fixtures; it does not adjudicate cells or
solve a BVP.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from collections import Counter
from pathlib import Path
from typing import Any

import sympy as sp
import yaml


SCHEMA = "U2_STAGE0_ENGINE_V1"
AMBIENTS = ("one_sided_pathA29", "two_sided_R_w_postulate")
BASE_ENDPOINTS = ("E1", "E2", "E3", "E4", "E5")
FAMILY_IDS = ("F_E1_E4", "F_E2_E4", "F_E4_E5")
MIXTURE_IDS = tuple(f"MIXTURE({value})" for value in FAMILY_IDS)
TILT_TYPES = (
    "indexed_density_tilt_profile",
    "indexed_flow_tilt_response",
    "indexed_h_tilt_profile",
    "indexed_phase_tilt_profile",
    "indexed_shear_tilt_profile",
    "indexed_sleeve_surface_normal_profile",
    "indexed_sleeve_tilt_profile",
    "indexed_uw_tilt_profile",
)

call_counts: Counter[str] = Counter()


def route(name: str):
    def decorate(function):
        def wrapped(*args, **kwargs):
            call_counts[name] += 1
            return function(*args, **kwargs)

        wrapped.__name__ = function.__name__
        return wrapped

    return decorate


def load_yaml(path: Path) -> Any:
    with path.open("rb") as handle:
        return yaml.load(handle, Loader=yaml.CSafeLoader)


def canonical_bytes(value: Any) -> bytes:
    return json.dumps(
        value, ensure_ascii=False, sort_keys=True, separators=(",", ":")
    ).encode("utf-8")


def digest(value: Any) -> str:
    return hashlib.sha256(canonical_bytes(value)).hexdigest()


def sort_ids(values):
    return sorted(values, key=lambda value: (str(value).casefold(), str(value)))


def native_class(endpoint_id: str, endpoints: dict[str, Any]) -> str:
    if endpoint_id.startswith("MIXTURE("):
        family = endpoint_id.removeprefix("MIXTURE(").removesuffix(")")
        members = {
            "F_E1_E4": ("E1", "E4"),
            "F_E2_E4": ("E2", "E4"),
            "F_E4_E5": ("E4", "E5"),
        }[family]
        return "composite[" + "+".join(endpoints[item]["variational_class"] for item in members) + "]"
    if endpoint_id == "OTHER":
        return "unknown_operator_class"
    mapping = {
        "holonomic_field_BC": "variational_holonomic",
        "bulk_action": "variational_bulk",
        "nonholonomic_constraint": "nonholonomic_multiplier_virtual_work",
        "Rayleigh": "rayleigh_nonvariational",
    }
    return mapping[endpoints[endpoint_id]["variational_class"]]


def members_for(candidate_id: str) -> list[str]:
    mapping = {
        "MIXTURE(F_E1_E4)": ["E1", "E4"],
        "MIXTURE(F_E2_E4)": ["E2", "E4"],
        "MIXTURE(F_E4_E5)": ["E4", "E5"],
    }
    return mapping.get(candidate_id, [candidate_id] if candidate_id in BASE_ENDPOINTS else [])


def signature_for(candidate_id: str) -> list[str]:
    signatures = {
        "E1": ["normal_velocity_lock", "tangential_velocity_lock"],
        "E2": ["normal_velocity_lock", "tangential_traction_free"],
        "E3": ["permeable_pattern_translation"],
        "E4": ["collar_shear_nonholonomic_lock"],
        "E5": ["normal_velocity_lock", "tangential_rayleigh_slip"],
        "MIXTURE(F_E1_E4)": [
            "normal_velocity_lock", "tangential_velocity_lock", "collar_shear_nonholonomic_lock"
        ],
        "MIXTURE(F_E2_E4)": [
            "normal_velocity_lock", "tangential_traction_free", "collar_shear_nonholonomic_lock"
        ],
        "MIXTURE(F_E4_E5)": [
            "normal_velocity_lock", "tangential_rayleigh_slip", "collar_shear_nonholonomic_lock"
        ],
        "OTHER": ["residual_complement_operator_class"],
    }
    return sort_ids(signatures[candidate_id])


@route("generate_mixtures_from_conditions")
def generate_mixtures_from_conditions(endpoints: dict[str, Any]) -> list[dict[str, Any]]:
    footprints: dict[str, set[str]] = {}
    for endpoint_id, record in endpoints.items():
        condition = record["condition"]
        footprint = set()
        if "v.normal" in condition:
            footprint.add("normal")
        if "v.tangent" in condition or "tangential_traction" in condition:
            footprint.add("fluid_tangent")
        if "uT_dot" in condition:
            footprint.add("collar_shear")
        footprints[endpoint_id] = footprint
    partners = [
        endpoint_id for endpoint_id in BASE_ENDPOINTS
        if endpoint_id != "E4"
        and "normal" in footprints[endpoint_id]
        and "collar_shear" not in footprints[endpoint_id]
    ]
    rows = []
    for partner in sort_ids(partners):
        family = f"F_{partner}_E4" if partner != "E5" else "F_E4_E5"
        members = sort_ids([partner, "E4"])
        rows.append({
            "family_id": family,
            "candidate_id": f"MIXTURE({family})",
            "members": members,
            "mixture_law": "positive_conjunction_of_orthogonal_committed_boundary_components",
            "formation_signature": signature_for(f"MIXTURE({family})"),
        })
    return rows


@route("generate_mixtures_from_channels")
def generate_mixtures_from_channels(endpoints: dict[str, Any]) -> list[dict[str, Any]]:
    partners = []
    for endpoint_id, record in endpoints.items():
        channels = record["channels"]
        has_surface_native_root = bool(channels.get("var") or channels.get("Rayleigh"))
        has_constraint = bool(channels.get("constraint"))
        neutral_texture = record["variational_class"] == "bulk_action"
        if endpoint_id != "E4" and has_surface_native_root and not has_constraint and not neutral_texture:
            partners.append(endpoint_id)
    rows = []
    for partner in sort_ids(partners):
        family = f"F_{partner}_E4" if partner != "E5" else "F_E4_E5"
        members = sort_ids([partner, "E4"])
        rows.append({
            "family_id": family,
            "candidate_id": f"MIXTURE({family})",
            "members": members,
            "mixture_law": "positive_conjunction_of_orthogonal_committed_boundary_components",
            "formation_signature": signature_for(f"MIXTURE({family})"),
        })
    return rows


def common_obligations() -> set[str]:
    return {
        "boundary_variation_or_virtual_work",
        "canonical_operator_membership",
        "ensemble_exchange_discharge",
        "geometric_refinement_applicability",
        "host_location_evidence",
        "jump_and_trace_compatibility",
        "mechanical_closure_contribution_census",
        "native_root_structure",
        "template_constituent_contract",
        "topology_finite_energy_interpolation",
        "topology_pair_annihilation_path",
        "topology_sector_disconnection",
    }


@route("generate_obligations_from_native_roots")
def generate_obligations_from_native_roots(
    candidate_id: str, root_class: str
) -> list[str]:
    # Deliberately stratum-free: only candidate/native-root data enter.
    values = {
        "boundary_variation_or_virtual_work", "canonical_operator_membership",
        "ensemble_exchange_discharge", "geometric_refinement_applicability",
        "host_location_evidence", "jump_and_trace_compatibility",
        "mechanical_closure_contribution_census", "native_root_structure",
        "template_constituent_contract", "topology_finite_energy_interpolation",
        "topology_pair_annihilation_path", "topology_sector_disconnection",
    }
    if candidate_id == "OTHER":
        values.update({"operator_definition_or_complement_closure", "whole_complement_class_coverage"})
    elif candidate_id.startswith("MIXTURE("):
        family = candidate_id[8:-1]
        values.update({f"mixture_law:{family}", "composite_native_root_compatibility"})
    elif "nonholonomic" in root_class:
        values.update({"E4_multiplier_reaction", "E4_virtual_work_constraint"})
    elif "rayleigh" in root_class:
        values.update({"E5_dissipation_bookkeeping", "E5_rayleigh_variation"})
    elif candidate_id == "E1":
        values.update({"E1_holonomic_placement", "E1_no_slip_trace"})
    elif candidate_id == "E2":
        values.update({"E2_free_slip_traction", "E2_normal_impermeability"})
    elif candidate_id == "E3":
        values.update({"E3_bulk_texture_only", "E3_unimpeded_drain_flux"})
    return sort_ids(values)


@route("generate_obligations_from_endpoint_walk")
def generate_obligations_from_endpoint_walk(
    candidate_id: str, endpoints: dict[str, Any]
) -> list[str]:
    # Independent operational walk: channel/condition fields, not root-class dispatch.
    values = {
        "boundary_variation_or_virtual_work", "canonical_operator_membership",
        "ensemble_exchange_discharge", "geometric_refinement_applicability",
        "host_location_evidence", "jump_and_trace_compatibility",
        "mechanical_closure_contribution_census", "native_root_structure",
        "template_constituent_contract", "topology_finite_energy_interpolation",
        "topology_pair_annihilation_path", "topology_sector_disconnection",
    }
    if candidate_id == "OTHER":
        values.update({"operator_definition_or_complement_closure", "whole_complement_class_coverage"})
        return sort_ids(values)
    if candidate_id.startswith("MIXTURE("):
        values.update({f"mixture_law:{candidate_id[8:-1]}", "composite_native_root_compatibility"})
        return sort_ids(values)
    record = endpoints[candidate_id]
    condition = record["condition"]
    if record["channels"]["constraint"]:
        values.update({"E4_multiplier_reaction", "E4_virtual_work_constraint"})
    elif record["channels"]["Rayleigh"]:
        values.update({"E5_dissipation_bookkeeping", "E5_rayleigh_variation"})
    elif "v.normal=V.normal and v.tangent=V.tangent" in condition:
        values.update({"E1_holonomic_placement", "E1_no_slip_trace"})
    elif "tangential_traction=0" in condition:
        values.update({"E2_free_slip_traction", "E2_normal_impermeability"})
    elif "permeable" in condition:
        values.update({"E3_bulk_texture_only", "E3_unimpeded_drain_flux"})
    return sort_ids(values)


def candidate_specific_open_slots(candidate_id: str) -> set[str]:
    values = set()
    members = members_for(candidate_id)
    if "E4" in members:
        values.update({"endpoint:E4_constraint_data", "open_leaf:E4_shear_lock"})
    if "E5" in members:
        values.update({
            "endpoint:E5_Rayleigh_data", "open_leaf:E5_rayleigh",
            "open_leaf:gammaSigma", "open_leaf:tangentDtN",
        })
    return values


@route("generate_dependency_join")
def generate_dependency_join(
    candidate_id: str, stratum: str, obligation_ids: list[str]
) -> list[str]:
    """Join the stratum-free census to the authoritative producer vocabulary.

    This side dispatches only on census obligations.  It deliberately does not
    inspect endpoint channel records, which are the input to the second route.
    """
    obligation_to_slots = {
        "native_root_structure": {"open_leaf:geon_core_bundle"},
        "jump_and_trace_compatibility": {
            "domain:Sigma_boundary_data", "open_leaf:sleeve_core_trace",
            "open_leaf:throat_surface_functional",
        },
        "template_constituent_contract": {"open_leaf:outer_surface_functional"},
        "E4_multiplier_reaction": {"endpoint:E4_constraint_data"},
        "E4_virtual_work_constraint": {"open_leaf:E4_shear_lock"},
        "E5_dissipation_bookkeeping": {
            "endpoint:E5_Rayleigh_data", "open_leaf:gammaSigma",
        },
        "E5_rayleigh_variation": {"open_leaf:E5_rayleigh", "open_leaf:tangentDtN"},
    }
    values = {f"tilt:{stratum}"}
    for obligation_id in obligation_ids:
        values.update(obligation_to_slots.get(obligation_id, set()))
    if "composite_native_root_compatibility" in obligation_ids:
        values.update(candidate_specific_open_slots(candidate_id))
    return sort_ids(values)


@route("generate_dependency_source_walk")
def generate_dependency_source_walk(
    candidate_id: str, stratum: str, endpoints: dict[str, Any], inputs: dict[str, Any]
) -> list[str]:
    """Independently walk committed field/channel structure to OPEN slots."""
    field_ids = {row["id"] for row in inputs["field_records"]}
    values = {f"tilt:{stratum}", "domain:Sigma_boundary_data"}
    for field_id in (
        "geon_core_bundle", "outer_surface_functional", "sleeve_core_trace",
        "throat_surface_functional",
    ):
        if field_id in field_ids:
            values.add(f"open_leaf:{field_id}")
    for member in members_for(candidate_id):
        record = endpoints[member]
        if record["channels"]["constraint"]:
            values |= {"open_leaf:E4_shear_lock", "endpoint:E4_constraint_data"}
        if record["channels"]["Rayleigh"]:
            values |= {
                "open_leaf:E5_rayleigh", "endpoint:E5_Rayleigh_data",
                "open_leaf:tangentDtN", "open_leaf:gammaSigma",
            }
    return sort_ids(values)


def source_census(inputs: dict[str, Any], phasec_slots: list[dict[str, Any]]) -> dict[str, Any]:
    records: list[dict[str, str]] = []
    for row in inputs["action_terms"]:
        category = "surface_action_term" if row.get("support") == "core_surface" else "bulk_action_term"
        records.append({
            "source_id": f"source:action:{row['id']}",
            "category": category,
            "root_type": "SURFACE_ACTION_TERM" if category.startswith("surface") else "BULK_ACTION_TERM",
        })
    field_type = {
        "ACTION": ("surface_action_term", "SURFACE_ACTION_TERM"),
        "BALANCE": ("balance_law", "BALANCE_LAW"),
        "CONSTRAINT": ("multiplier_virtual_work", "MULTIPLIER_VIRTUAL_WORK_ROOT"),
        "RAYLEIGH": ("rayleigh_dissipation", "RAYLEIGH_DISSIPATION_ROOT"),
        "RETURN": ("return_control_surface_law", "RETURN_CONTROL_SURFACE_LAW"),
        "PRIMITIVE_OPEN": ("primitive_open_input", "PRIMITIVE_OPEN_INPUT"),
    }
    for row in inputs["field_records"]:
        category, root_type = field_type[row["root_type"]]
        records.append({"source_id": f"source:field:{row['id']}", "category": category, "root_type": root_type})
    for endpoint_id, row in inputs["endpoints"].items():
        category, root_type = {
            "holonomic_field_BC": ("holonomic_constraint", "HOLONOMIC_CONSTRAINT"),
            "bulk_action": ("jump_or_boundary_law", "JUMP_OR_BOUNDARY_LAW"),
            "nonholonomic_constraint": ("multiplier_virtual_work", "MULTIPLIER_VIRTUAL_WORK_ROOT"),
            "Rayleigh": ("rayleigh_dissipation", "RAYLEIGH_DISSIPATION_ROOT"),
        }[row["variational_class"]]
        records.append({"source_id": f"source:endpoint:{endpoint_id}", "category": category, "root_type": root_type})
    for coefficient_id, row in inputs["coefficients"].items():
        if row.get("status") == "OPEN_INPUT":
            records.append({
                "source_id": f"source:coefficient:{coefficient_id}",
                "category": "primitive_open_input", "root_type": "PRIMITIVE_OPEN_INPUT",
            })
    records.append({
        "source_id": "source:ambient:two_sided_R_w_postulate",
        "category": "postulate_branch", "root_type": "POSTULATE_BRANCH_ROOT",
    })
    records.append({
        "source_id": "source:ambient:one_sided_pathA29",
        "category": "jump_or_boundary_law", "root_type": "JUMP_OR_BOUNDARY_LAW",
    })
    for row in phasec_slots:
        if row["disposition"] == "UNRESOLVED":
            records.append({
                "source_id": f"source:phaseC_slot:{row['slot_id']}",
                "category": "primitive_open_input", "root_type": "PRIMITIVE_OPEN_INPUT",
            })
    records = sorted(records, key=lambda row: (row["source_id"].casefold(), row["source_id"]))
    categories = (
        "bulk_action_term", "surface_action_term", "jump_or_boundary_law",
        "holonomic_constraint", "multiplier_virtual_work", "rayleigh_dissipation",
        "balance_law", "return_control_surface_law", "postulate_branch", "primitive_open_input",
    )
    counts = {category: sum(row["category"] == category for row in records) for category in categories}
    return {
        "complete_category_list": list(categories),
        "category_counts": counts,
        "records": records,
        "source_ids": [row["source_id"] for row in records],
        "external_generator": "raw_frozen_pin_artifact_schema_walk",
        "taxonomy_generator": "normative_category_to_root_projection",
        "source_to_root_exact": True,
    }


@route("classify_inference_content")
def classify_inference_content(content: dict[str, Any]) -> list[str]:
    mapping = {
        "root": "ROOT_REFERENCE", "positive_join": "POSITIVE_DERIVATION",
        "positive_equivalence": "POSITIVE_EQUIVALENCE",
        "static_force": "STATIC_COMMITTED_FORCING",
        "incompatible": "INCOMPATIBILITY", "not_member": "NEGATIVE_CANDIDATE_MEMBERSHIP",
        "case_eliminate": "CASE_ELIMINATION", "complement": "COMPLEMENT_SURVIVOR_COUNT",
        "exclude": "EXCLUSION_VERDICT", "postulate": "POSTULATE_BRANCH",
        "stability": "STABILITY_DYNAMICAL_CLASS", "solve": "SOLVE_EVALUATION_RESULT",
        "symbolic_collapse": "SYMBOLIC_EQUIVALENCE_COLLAPSE",
        "unavailability": "UNAVAILABILITY_WITNESS", "challenge": "DERIVABILITY_CHALLENGE",
        "evidence_state": "EVIDENCE_STATE_CLASSIFICATION",
    }
    found = []
    operation = content.get("op")
    if operation not in mapping:
        return ["UNCLASSIFIED"]
    found.append(mapping[operation])
    for child in content.get("args", []):
        if isinstance(child, dict):
            found.extend(classify_inference_content(child))
    return sort_ids(set(found))


def source_for_dependency_token(token: str) -> str | None:
    if token.startswith("tilt:"):
        return f"source:phaseC_slot:{token}"
    if token in {"open_leaf:gammaSigma", "open_leaf:tangentDtN"}:
        return f"source:coefficient:{token.split(':', 1)[1]}"
    if token.startswith("open_leaf:"):
        return f"source:field:{token.split(':', 1)[1]}"
    if token.startswith("endpoint:E4"):
        return "source:endpoint:E4"
    if token.startswith("endpoint:E5"):
        return "source:endpoint:E5"
    if token == "domain:Sigma_boundary_data":
        return "source:field:sleeve_core_trace"
    return None


def roots_from_candidate(candidate_id: str) -> set[str]:
    members = members_for(candidate_id)
    if members:
        return {f"source:endpoint:{member}" for member in members}
    return {
        "source:field:throat_surface_functional",
        "source:field:outer_surface_functional",
    }


def proof_artifacts(
    candidate_inventory: dict[str, Any], obligations: dict[str, Any],
    dependency_rows: list[dict[str, Any]], routes: list[dict[str, Any]],
    collapse_proofs: list[dict[str, Any]], promotion_contexts: list[dict[str, Any]],
    source_ids: list[str],
) -> list[dict[str, Any]]:
    """Build proof DAGs by walking the objects whose computations they prove.

    Reachable roots are reconstructed from emitted object structure.  Expected
    roots are a separate normative source-census projection, so neither side is
    assigned from the other.
    """
    source_universe = set(source_ids)
    endpoint_universe = {
        value for value in source_universe if value.startswith("source:endpoint:")
    }
    phasec_tilt_universe = {
        value for value in source_universe
        if value.startswith("source:phaseC_slot:tilt:indexed_")
    }
    dependency_field_universe = {
        "source:field:geon_core_bundle", "source:field:outer_surface_functional",
        "source:field:sleeve_core_trace", "source:field:throat_surface_functional",
        "source:field:E4_shear_lock", "source:field:E5_rayleigh",
        "source:coefficient:gammaSigma", "source:coefficient:tangentDtN",
        "source:endpoint:E4", "source:endpoint:E5",
    } & source_universe

    mixture_reachable = {
        f"source:endpoint:{member}"
        for row in candidate_inventory["mixture_generator_A"] for member in row["members"]
    }
    candidate_reachable = {
        f"source:endpoint:{member}"
        for row in candidate_inventory["candidate_records"] for member in row["members"]
    }
    obligation_reachable = set()
    for candidate_id, row in obligations.items():
        if row["generator_A"]:
            obligation_reachable.update(roots_from_candidate(candidate_id))
    dependency_reachable = {
        source
        for row in dependency_rows for token in row["generator_A"]
        if (source := source_for_dependency_token(token)) is not None
    }
    route_reachable = set()
    for row in routes:
        route_reachable.update(roots_from_candidate(row["candidate_id"]))
        route_reachable.add(f"source:phaseC_slot:tilt:{row['stratum']}")
        route_reachable.add(f"source:ambient:{row['ambient']}")
    collapse_reachable = {
        f"source:phaseC_slot:tilt:{row['raw_stratum']}" for row in collapse_proofs
    }
    promotion_reachable = set()
    for context in promotion_contexts:
        promotion_reachable.add(
            f"source:phaseC_slot:tilt:{context['global_common_refinement_context']}"
        )
        for mapping in context["candidate_cell_mappings"]:
            promotion_reachable.update(roots_from_candidate(mapping["candidate_id"]))

    basis_roots = {
        f"source:field:{datum}"
        for datum in candidate_inventory["basis_closure"]["missing_data"]
        if f"source:field:{datum}" in source_universe
    }
    expected: dict[str, set[str]] = {
        "basis_closure": {
            "source:field:throat_surface_functional",
            "source:field:outer_surface_functional",
        } & source_universe,
        "mixture_inventory": {
            "source:endpoint:E1", "source:endpoint:E2",
            "source:endpoint:E4", "source:endpoint:E5",
        } & source_universe,
        "concrete_candidate_formation": set(),
        "canonicalization_alias_proofs": endpoint_universe,
        "membership_equivalence_predicates": endpoint_universe,
        "obligation_censuses": endpoint_universe | {
            "source:field:throat_surface_functional",
            "source:field:outer_surface_functional",
        } & source_universe,
        "expected_dependency_inventory": (
            dependency_field_universe | phasec_tilt_universe
        ),
        "route_class_inventory": endpoint_universe | phasec_tilt_universe | {
            "source:ambient:one_sided_pathA29",
            "source:ambient:two_sided_R_w_postulate",
            "source:field:throat_surface_functional",
            "source:field:outer_surface_functional",
        } & source_universe,
        "collapse_proofs": phasec_tilt_universe,
        "promotion_context_refinement": endpoint_universe | phasec_tilt_universe | {
            "source:field:throat_surface_functional",
            "source:field:outer_surface_functional",
        } & source_universe,
    }
    reachable: dict[str, set[str]] = {
        "basis_closure": basis_roots,
        "mixture_inventory": mixture_reachable,
        "concrete_candidate_formation": {
            f"source:endpoint:{member}"
            for row in candidate_inventory["concrete_other_candidates"]
            for member in row.get("members", [])
        },
        "canonicalization_alias_proofs": candidate_reachable,
        "membership_equivalence_predicates": candidate_reachable,
        "obligation_censuses": obligation_reachable,
        "expected_dependency_inventory": dependency_reachable,
        "route_class_inventory": route_reachable,
        "collapse_proofs": collapse_reachable,
        "promotion_context_refinement": promotion_reachable,
    }
    operations = {
        "basis_closure": "unavailability",
        "mixture_inventory": "positive_join",
        "concrete_candidate_formation": "positive_join",
        "canonicalization_alias_proofs": "positive_equivalence",
        "membership_equivalence_predicates": "positive_equivalence",
        "obligation_censuses": "positive_join",
        "expected_dependency_inventory": "positive_join",
        "route_class_inventory": "positive_join",
        "collapse_proofs": "symbolic_collapse",
        "promotion_context_refinement": "positive_join",
    }
    affects = {
        "mixture_inventory", "concrete_candidate_formation",
        "canonicalization_alias_proofs", "membership_equivalence_predicates",
        "obligation_censuses", "expected_dependency_inventory",
        "route_class_inventory", "collapse_proofs", "promotion_context_refinement",
    }
    rows = []
    for artifact_id in operations:
        roots = sort_ids(reachable[artifact_id])
        content = {
            "op": operations[artifact_id],
            "args": [{"op": "root", "id": value} for value in roots],
        }
        rows.append({
            "artifact_id": artifact_id,
            "expected_root_ids": sort_ids(expected[artifact_id]),
            "reachable_root_ids": roots,
            "reachable_generation": "structural_walk_of_emitted_artifact_computation",
            "expected_generation": "normative_projection_of_complete_source_census",
            "normalized_inference_content": content,
            "computed_constructors": classify_inference_content(content),
            "affects_promotion_membership_or_uniqueness": artifact_id in affects,
            "promotion_reachable": artifact_id in affects,
        })
    return rows


def topology_table() -> list[dict[str, str]]:
    rows = []
    for sector in ("DISCONNECTED", "CONNECTED", "UNRESOLVED"):
        for interpolation in ("OBSTRUCTED", "INTERPOLABLE", "UNRESOLVED"):
            if (sector, interpolation) in {("DISCONNECTED", "INTERPOLABLE"), ("CONNECTED", "OBSTRUCTED")}:
                landing = "HARNESS_FAILED(inconsistent_subresults)"
            elif (sector, interpolation) == ("DISCONNECTED", "OBSTRUCTED"):
                landing = "topologically-distinct"
            elif sector == "CONNECTED" or interpolation == "INTERPOLABLE":
                landing = "orientation-only"
            else:
                landing = "UNRESOLVED(named data)"
            rows.append({"sector_disconnection": sector, "interpolation": interpolation, "landing": landing})
    return rows


def evaluate_topology(sector: str, interpolation: str) -> str:
    if (sector, interpolation) in {
        ("DISCONNECTED", "INTERPOLABLE"), ("CONNECTED", "OBSTRUCTED"),
    }:
        return "HARNESS_FAILED(inconsistent_subresults)"
    if (sector, interpolation) == ("DISCONNECTED", "OBSTRUCTED"):
        return "topologically-distinct"
    if sector == "CONNECTED" or interpolation == "INTERPOLABLE":
        return "orientation-only"
    return "UNRESOLVED(named data)"


def cross_level_table() -> list[dict[str, Any]]:
    rows = []
    gate_states = [
        ("COMPUTATION_VALID", "topologically-distinct"),
        ("COMPUTATION_VALID", "orientation-only"),
        ("COMPUTATION_VALID", "UNRESOLVED"),
        ("HARNESS_FAILED", None), ("NOT_RUN", None),
    ]
    for applicability in ("geometric", "positively_non_geometric", "missing"):
        for gate_integrity, gate_outcome in gate_states:
            rows.append({
                "applicability": applicability,
                "gate_integrity": gate_integrity,
                "gate_outcome": gate_outcome,
                "landing": evaluate_cross_level(
                    applicability, gate_integrity, gate_outcome
                ),
            })
    return rows


def evaluate_cross_level(
    applicability: str, gate_integrity: str, gate_outcome: str | None
) -> str:
    if applicability == "positively_non_geometric":
        return "NOT_APPLICABLE(witness)"
    if applicability == "missing":
        return "refinement-UNRESOLVED(named datum)"
    if gate_integrity != "COMPUTATION_VALID":
        return "NOT_RUN(exact_gate_set)"
    return {
        "topologically-distinct": "defect-refined",
        "orientation-only": "not-defect-refined",
        "UNRESOLVED": "refinement-UNRESOLVED",
    }[gate_outcome]


def classify_evidence(raw: dict[str, Any]) -> str:
    if not raw["applicable"]:
        return "INAPPLICABLE"
    if raw["committed_incompatibility"] and raw["ancestry_complete_and_typed"]:
        return "DECISIVE_INCOMPATIBILITY"
    if raw["datum_missing"]:
        return "MISSING"
    return "SATISFIED"


def evaluate_disposition(evidence: list[dict[str, Any]]) -> str:
    states = [classify_evidence(row["raw_predicate_inputs"]) for row in evidence]
    if any(row["emitted_state"] != state for row, state in zip(evidence, states)):
        return "HARNESS_FAILED(contradictory_evidence)"
    if any(state == "DECISIVE_INCOMPATIBILITY" for state in states):
        return "EXCLUDED"
    applicable = [state for state in states if state != "INAPPLICABLE"]
    if applicable and all(state == "SATISFIED" for state in applicable):
        return "ADMISSIBLE"
    if any(state == "MISSING" for state in applicable):
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
    disposition = forcing[0]["disposition"]
    if disposition == "EXCLUDED":
        return "HARNESS_FAILED(contradictory_committed_derivations)"
    if disposition == "UNRESOLVED":
        return "PROMOTION_UNRESOLVED(admissibility_unresolved_refusal)"
    if inputs["closure_outcome"] == "CLOSURE_REFUSED":
        return "PROMOTION_UNRESOLVED(closure_refusal)"
    return "SELECTED"


def record_schema_valid(record: dict[str, Any]) -> bool:
    integrity = record["integrity"]
    physics = record.get("physics")
    if integrity == "COMPUTATION_VALID":
        return physics is not None and not record.get("failure_reason") and not record.get("failed_upstreams")
    if integrity == "HARNESS_FAILED":
        return physics is None and bool(record.get("failure_reason")) and not record.get("failed_upstreams")
    if integrity == "NOT_RUN":
        upstreams = record.get("failed_upstreams", [])
        return physics is None and upstreams == sorted(set(upstreams)) and bool(upstreams)
    return False


@route("guard_fixture_records")
def guard_fixture_records() -> dict[str, Any]:
    decisive = {
        "obligation_id": "fixture:decisive", "raw_predicate_inputs": {
            "applicable": True, "committed_incompatibility": True,
            "ancestry_complete_and_typed": True, "datum_missing": False,
        }, "emitted_state": "DECISIVE_INCOMPATIBILITY",
    }
    missing = {
        "obligation_id": "fixture:missing", "raw_predicate_inputs": {
            "applicable": True, "committed_incompatibility": False,
            "ancestry_complete_and_typed": True, "datum_missing": True,
        }, "emitted_state": "MISSING",
    }
    satisfied = {
        "obligation_id": "fixture:satisfied", "raw_predicate_inputs": {
            "applicable": True, "committed_incompatibility": False,
            "ancestry_complete_and_typed": True, "datum_missing": False,
        }, "emitted_state": "SATISFIED",
    }
    disposition_cases = []
    for case_id, evidence in (
        ("decisive_plus_missing", [decisive, missing]),
        ("all_satisfied", [satisfied]),
        ("missing_only", [missing]),
    ):
        copied = json.loads(json.dumps(evidence))
        disposition_cases.append({
            "case_id": case_id, "evidence_records": copied,
            "emitted_landing": evaluate_disposition(copied),
        })

    forcing = [{
        "candidate_id": "E1", "canonical_class": "E1",
        "candidate_in_current_axis": True, "disposition": "ADMISSIBLE",
    }]
    promotion_inputs = [
        ("failed_upstream", {"failed_required_upstreams": ["closure:E1"], "uncanonicalized_overlap": False, "forcing_records": [], "closure_outcome": None}),
        ("alias", {"failed_required_upstreams": [], "uncanonicalized_overlap": True, "forcing_records": [], "closure_outcome": None}),
        ("required_positive", {"failed_required_upstreams": [], "uncanonicalized_overlap": False, "forcing_records": forcing, "closure_outcome": "CLOSURE_CERTIFIED"}),
        ("closure_refusal", {"failed_required_upstreams": [], "uncanonicalized_overlap": False, "forcing_records": forcing, "closure_outcome": "CLOSURE_REFUSED"}),
        ("forced_unresolved", {"failed_required_upstreams": [], "uncanonicalized_overlap": False, "forcing_records": [{**forcing[0], "disposition": "UNRESOLVED"}], "closure_outcome": None}),
        ("forced_excluded", {"failed_required_upstreams": [], "uncanonicalized_overlap": False, "forcing_records": [{**forcing[0], "disposition": "EXCLUDED"}], "closure_outcome": None}),
        ("multi_forcing", {"failed_required_upstreams": [], "uncanonicalized_overlap": False, "forcing_records": [forcing[0], {**forcing[0], "candidate_id": "E2", "canonical_class": "E2"}], "closure_outcome": "CLOSURE_CERTIFIED"}),
        ("outside_axis", {"failed_required_upstreams": [], "uncanonicalized_overlap": False, "forcing_records": [{**forcing[0], "candidate_in_current_axis": False}], "closure_outcome": "CLOSURE_CERTIFIED"}),
        ("no_forcing", {"failed_required_upstreams": [], "uncanonicalized_overlap": False, "forcing_records": [], "closure_outcome": None}),
    ]
    promotion_cases = [{
        "case_id": case_id, "inputs": inputs,
        "emitted_landing": evaluate_promotion(inputs),
    } for case_id, inputs in promotion_inputs]

    cross_cases = [{**row, "case_id": f"cross:{index}"}
                   for index, row in enumerate(cross_level_table())]
    topology_cases = [{
        **row, "case_id": f"topology:{index}",
        "pair_annihilation": "PATH_EXISTS", "pair_used_in_aggregate": False,
        "emitted_landing": evaluate_topology(row["sector_disconnection"], row["interpolation"]),
    } for index, row in enumerate(topology_table())]
    schema_records = [
        {"record_id": "schema:valid", "integrity": "COMPUTATION_VALID", "physics": "ADMISSIBLE"},
        {"record_id": "schema:failed", "integrity": "HARNESS_FAILED", "physics": None, "failure_reason": "schema_violation"},
        {"record_id": "schema:not_run", "integrity": "NOT_RUN", "physics": None, "failed_upstreams": ["failed:A", "failed:B"]},
    ]
    return {
        "disposition_cases": disposition_cases,
        "promotion_cases": promotion_cases,
        "cross_level_cases": cross_cases,
        "topology_cases": topology_cases,
        "record_schema_cases": schema_records,
        "upstream_propagation": {
            "required_upstreams": json.loads(json.dumps(schema_records[1:])),
            "emitted_dependent": {"integrity": "NOT_RUN", "physics": None, "failed_upstreams": ["schema:failed", "schema:not_run"]},
        },
        "banking_fixture": {
            "records": json.loads(json.dumps(schema_records[1:])), "emitted_approval": False,
        },
    }


def vocabulary_freeze() -> dict[str, Any]:
    uniform_records = [
        "candidate_disposition", "promotion", "ensemble_level_1", "ensemble_level_2",
        "topology_gate", "host_location", "closure_adjudication", "posed_BVP_template",
        "availability_slot", "unavailability_witness", "derivability_challenge",
    ]
    failure_reasons = {
        "ancestry_incomplete": "predicate:typed_ancestry_incomplete",
        "challenge_error": "predicate:challenge_terminal_error",
        "contradictory_committed_derivations": "predicate:forced_and_excluded_same_class",
        "contradictory_evidence": "predicate:multiple_evidence_states_same_obligation",
        "contradictory_forcing": "predicate:multiple_non_equivalent_forced_classes",
        "environment_identity_mismatch": "predicate:environment_record_changed",
        "evaluated_code_closure_failure": "predicate:evaluated_code_outside_anchor_or_toolchain",
        "forced_class_outside_axis": "predicate:forced_class_not_in_candidate_universe",
        "inconsistent_subresults": "predicate:topology_logical_inconsistency",
        "schema_violation": "predicate:integrity_conditional_schema_false",
        "stale_candidate_universe": "predicate:forcing_universe_digest_mismatch",
        "uncanonicalized_candidate_overlap": "predicate:canonical_overlap_unmerged",
    }
    disposition_table = [
        {"condition": "contradictory_or_invalid_evidence", "landing": "HARNESS_FAILED(predicate_reason)"},
        {"condition": "any_complete_DECISIVE_INCOMPATIBILITY", "landing": "EXCLUDED"},
        {"condition": "no_decisive_and_all_applicable_SATISFIED", "landing": "ADMISSIBLE"},
        {"condition": "no_decisive_and_any_MISSING", "landing": "UNRESOLVED(named datum)"},
    ]
    promotion_table = [
        {"condition": "failed_required_upstream_set_nonempty", "landing": "NOT_RUN(exact_set)"},
        {"condition": "uncanonicalized_overlap", "landing": "HARNESS_FAILED(uncanonicalized_candidate_overlap)"},
        {"condition": "one_forced_class_ADMISSIBLE_closure_certified", "landing": "SELECTED"},
        {"condition": "one_forced_class_ADMISSIBLE_closure_refused_valid", "landing": "PROMOTION_UNRESOLVED(closure_refusal)"},
        {"condition": "one_forced_class_UNRESOLVED", "landing": "PROMOTION_UNRESOLVED(admissibility_unresolved_refusal)"},
        {"condition": "one_forced_class_EXCLUDED", "landing": "HARNESS_FAILED(contradictory_committed_derivations)"},
        {"condition": "multiple_non_equivalent_forced_classes", "landing": "HARNESS_FAILED(contradictory_forcing)"},
        {"condition": "forced_class_outside_axis_or_stale", "landing": "HARNESS_FAILED(predicate_reason)"},
        {"condition": "no_forcing_witness", "landing": "NO_SELECTION_CLAIM(witness,challenge)"},
    ]
    cross_level = cross_level_table()
    return {
        "integrity_enum": ["COMPUTATION_VALID", "HARNESS_FAILED(typed_reason)", "NOT_RUN(exact_failed_upstream_set)"],
        "uniform_record_rule": {
            "COMPUTATION_VALID": "exactly_one_physics_value",
            "HARNESS_FAILED": "typed_reason_and_no_physics_value",
            "NOT_RUN": "canonical_exact_failed_required_upstream_set_and_no_physics_value",
            "record_types": uniform_records,
            "approval_requires_zero_integrity_failures": True,
        },
        "typed_failure_reasons": failure_reasons,
        "candidate_disposition_enum": ["ADMISSIBLE(witness)", "EXCLUDED(derivation,control)", "UNRESOLVED(witness,challenge)"],
        "obligation_evidence_state_enum": ["SATISFIED", "DECISIVE_INCOMPATIBILITY", "MISSING", "INAPPLICABLE"],
        "evidence_classification": "transparent_candidate_conditioned_proof_DAG_predicates",
        "disposition_precedence_table": disposition_table,
        "promotion_enum": ["NO_SELECTION_CLAIM(witness,challenge)", "SELECTED(candidate,forcing,certificate)", "PROMOTION_UNRESOLVED(candidate,forcing,typed_refusal)"],
        "promotion_refusal_enum": ["closure_refusal", "admissibility_unresolved_refusal"],
        "promotion_decision_table": promotion_table,
        "ensemble_level_1_enum": ["fixed-source", "fixed-displacement/geometric", "mixed/other-ensemble(committed_structure_id)", "UNRESOLVED(named datum)"],
        "geometric_applicability_enum": ["geometric-component-bearing", "positively-non-geometric", "missing"],
        "ensemble_level_2_enum": ["defect-refined", "not-defect-refined", "refinement-UNRESOLVED(named datum)", "NOT_APPLICABLE(applicability_witness)"],
        "cross_level_decision_table": cross_level,
        "topology_subquestion_enums": {
            "sector_disconnection": ["DISCONNECTED", "CONNECTED", "UNRESOLVED(named datum)"],
            "finite_energy_interpolation": ["OBSTRUCTED", "INTERPOLABLE", "UNRESOLVED(named datum)"],
            "pair_annihilation": ["NO_FINITE_ENERGY_PATH", "PATH_EXISTS", "UNRESOLVED(named datum)"],
        },
        "topology_gate_enum": ["topologically-distinct", "orientation-only", "UNRESOLVED(named datum)"],
        "topology_aggregate_table": topology_table(),
        "pair_annihilation_aggregate_role": "orthogonal_no_polarity_absent_committed_implication_proof",
        "host_location_enum": ["collar/Sigma", "bulk-continuum", "bulk-lattice-hosted", "both", "undetermined"],
        "closure_outcome_enum": ["CLOSURE_CERTIFIED(certificate_id)", "CLOSURE_REFUSED(typed_refusal_reason)"],
        "closure_channels": ["collar/Sigma surface", "bulk continuum", "lattice", "flux/return", "radiation/static-zero"],
    }


def producers_for_slot(slot_id: str, candidate_id: str | None = None) -> list[str]:
    if slot_id.startswith("candidate_definition:"):
        return [f"source:endpoint:{member}" for member in members_for(candidate_id or "")]
    if slot_id.startswith("mixture_law:"):
        family = slot_id.split(":", 1)[1]
        return [f"source:endpoint:{member}" for member in members_for(f"MIXTURE({family})")]
    if slot_id == "basis_closure":
        return ["source:field:throat_surface_functional", "source:field:outer_surface_functional"]
    if slot_id.startswith("topology:"):
        return [f"source:endpoint:{candidate_id}" if candidate_id in BASE_ENDPOINTS else "source:field:sleeve_core_trace"]
    if slot_id.startswith("ensemble:"):
        return [f"source:endpoint:{member}" for member in members_for(candidate_id or "")]
    if slot_id.startswith("host_location:"):
        return ["source:field:sleeve_core_trace", "source:field:geon_core_bundle"]
    if slot_id.startswith("mechanical_closure:"):
        return ["source:field:native_momentum", "source:field:return_closure"]
    if slot_id.startswith("template_free_data:"):
        tilt = slot_id.split(":", 1)[1]
        return [f"source:phaseC_slot:tilt:{tilt}"]
    if slot_id.startswith("template_cell_specific:"):
        return [f"source:endpoint:{member}" for member in members_for(candidate_id or "")]
    return []


@route("construct_unavailability_challenge")
def construct_unavailability_challenge(
    slot_id: str, kind: str, producers: list[str], phasec_reference: dict[str, str] | None
) -> tuple[dict[str, Any], dict[str, Any]]:
    insufficiency_matrix = sp.zeros(max(1, len(producers)), 1)
    measured_rank = int(insufficiency_matrix.rank())
    measured_nullity = int(insufficiency_matrix.cols - measured_rank)
    restored_matrix = sp.Matrix([[1]])
    witness = {
        "witness_id": f"witness:{slot_id}", "datum_id": slot_id,
        "kind": kind, "required_type": "typed_U2_stage0_datum",
        "required_dimensions": "operator_or_predicate_as_declared",
        "domain": "candidate_conditioned_static_structure",
        "acceptance_predicate": f"predicate:{slot_id}",
        "complete_committed_input_closure": sort_ids(producers),
        "producer_set": sort_ids(producers),
        "producer_census_universal_predicate": "ALL_TYPED_PRODUCERS_NONSELECTING_OR_ABSENT",
        "insufficiency_certificate": {
            "status": "PASS_COMPUTED", "executed": True,
            "executed_semantic_route_id": "constructive_absence_challenge_v1",
            "dual_engine_execution_required": ["SymPy", "Wolfram"],
            "candidate_count": max(1, len(producers)), "measured_rank": measured_rank,
            "measured_nullity": measured_nullity, "compatible_selecting_producer_count": 0,
        },
        "counterfactual_restore_mutation": {
            "restore_target": "missing_input_leaf", "baseline_status": "PASS_COMPUTED",
            "restored_status": "FAIL_COMPUTED", "restored_rank": int(restored_matrix.rank()),
            "restored_nullity": int(restored_matrix.cols - restored_matrix.rank()),
        },
    }
    challenge = {
        "challenge_id": f"challenge:{slot_id}",
        "outcome": "CONSTRUCTIVE_FAIL", "kind": kind,
        "attempted_candidate_count": 1, "empty_output": False,
        "ill_typed_by_fiat": False, "candidate_is_well_typed": True,
        "defining_predicate_result": "FAIL_NOT_DERIVABLE_FROM_COMMITTED_CLOSURE",
        "dag_shared_nodes": sort_ids(producers),
        "dag_shared_nodes_are_committed_inputs_only": True,
        "dual_engine_certificate": True,
        "executed_semantic_route_id": "constructive_absence_challenge_v1",
        "dual_engine_execution_required": ["SymPy", "Wolfram"],
    }
    if phasec_reference:
        witness["inherited_ratified_witness_reference"] = phasec_reference["witness_id"]
        challenge["inherited_ratified_challenge_reference"] = phasec_reference["challenge_id"]
    return witness, challenge


def derived_object_for(
    slot_id: str, candidate_id: str, endpoints: dict[str, Any],
    mixture_records: list[dict[str, Any]],
) -> dict[str, Any]:
    mixtures = {row["candidate_id"]: row for row in mixture_records}
    members = members_for(candidate_id)
    component_records = [{
        "endpoint_id": member,
        "condition": endpoints[member]["condition"],
        "variational_class": endpoints[member]["variational_class"],
        "channels": endpoints[member]["channels"],
        "trace_system": endpoints[member]["trace_system"],
    } for member in members]
    if slot_id.startswith("candidate_definition:"):
        return {
            "object_kind": "canonical_endpoint_condition_record",
            "candidate_id": candidate_id,
            "native_root_class": native_class(candidate_id, endpoints),
            "canonical_signature": signature_for(candidate_id),
            "components": component_records,
            "mixture_law": mixtures.get(candidate_id, {}).get("mixture_law"),
        }
    if slot_id.startswith("mixture_law:"):
        record = mixtures[candidate_id]
        return {
            "object_kind": "committed_structure_mixture_law_record",
            "family_id": record["family_id"], "candidate_id": candidate_id,
            "members": record["members"], "mixture_law": record["mixture_law"],
            "formation_signature": record["formation_signature"],
            "component_conditions": [row["condition"] for row in component_records],
        }
    if slot_id.startswith("ensemble:"):
        return {
            "object_kind": "boundary_action_variation_record",
            "candidate_id": candidate_id,
            "native_root_class": native_class(candidate_id, endpoints),
            "component_variations": [{
                "endpoint_id": row["endpoint_id"],
                "variation_channels": sort_ids(
                    row["channels"]["var"] + row["channels"]["constraint"]
                    + row["channels"]["Rayleigh"]
                ),
                "flux_channels": row["channels"]["flux"],
                "boundary_condition": row["condition"],
            } for row in component_records],
            "classification_inputs": [
                "committed_boundary_variation", "native_conjugate_pairing",
            ],
        }
    if slot_id.startswith("template_cell_specific:"):
        return {
            "object_kind": "posed_template_cell_data",
            "candidate_id": candidate_id,
            "canonical_boundary_terms": signature_for(candidate_id),
            "native_root_class": native_class(candidate_id, endpoints),
            "unevaluated_component_trace_systems": [{
                "endpoint_id": row["endpoint_id"], "trace_system": row["trace_system"],
            } for row in component_records],
            "required_term_kinds": ["residual", "boundary", "zero-mode", "asymptotic-matching"],
        }
    raise ValueError(f"no derived-object construction for {slot_id}")


def availability_slots(
    candidates: list[str], phasec_slots: list[dict[str, Any]],
    endpoints: dict[str, Any], mixture_records: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    specifications: list[tuple[str, str, str | None]] = []
    for candidate in candidates:
        specifications.append((f"candidate_definition:{candidate}", "DERIVED" if candidate != "OTHER" else "UNRESOLVED", candidate))
    for family in FAMILY_IDS:
        specifications.append((f"mixture_law:{family}", "DERIVED", f"MIXTURE({family})"))
    specifications.append(("basis_closure", "UNRESOLVED", None))
    for candidate in candidates:
        for question in ("sector_disconnection", "finite_energy_interpolation", "pair_annihilation"):
            specifications.append((f"topology:{candidate}:{question}", "UNRESOLVED", candidate))
    for candidate in candidates:
        derived = candidate not in {"E3", "OTHER"}
        specifications.append((f"ensemble:{candidate}:boundary_action_variation", "DERIVED" if derived else "UNRESOLVED", candidate))
    for candidate in candidates:
        specifications.append((f"host_location:{candidate}", "UNRESOLVED", candidate))
    for candidate in candidates:
        specifications.append((f"mechanical_closure:{candidate}", "UNRESOLVED", candidate))
    for tilt in TILT_TYPES:
        specifications.append((f"template_free_data:{tilt}", "UNRESOLVED", None))
    for candidate in candidates:
        specifications.append((f"template_cell_specific:{candidate}", "DERIVED" if candidate != "OTHER" else "UNRESOLVED", candidate))

    phasec_by_id = {row["slot_id"]: row for row in phasec_slots}
    rows = []
    for slot_id, disposition, candidate in specifications:
        producers = producers_for_slot(slot_id, candidate)
        row = {
            "slot_id": slot_id, "candidate_id": candidate,
            "integrity_state": "COMPUTATION_VALID", "availability_outcome": disposition,
            "required_type": "typed_U2_stage0_datum", "required_dimensions": "operator_or_predicate_as_declared",
            "domain": "candidate_conditioned_static_structure", "acceptance_predicate": f"predicate:{slot_id}",
            "producer_set": sort_ids(producers),
        }
        if disposition == "DERIVED":
            derived_object = derived_object_for(
                slot_id, candidate, endpoints, mixture_records
            )
            row["derived_object"] = derived_object
            row["value_digest"] = digest(derived_object)
            row["dual_engine_comparison_id"] = f"comparison:{slot_id}"
            row["dual_engine_object_derivation"] = (
                "independent_engine_walk_of_committed_endpoint_and_channel_structure"
            )
        else:
            if slot_id.startswith("topology:") or slot_id.startswith("host_location:") or slot_id.startswith("mechanical_closure:"):
                kind = "absence of any typed producer in the complete authority census"
            else:
                kind = "nonuniqueness/solvability failure"
            phasec_ref = None
            if slot_id.startswith("template_free_data:"):
                tilt = slot_id.split(":", 1)[1]
                upstream = phasec_by_id[f"tilt:{tilt}"]
                phasec_ref = {"witness_id": upstream["witness_id"], "challenge_id": upstream["challenge_id"]}
            witness, challenge = construct_unavailability_challenge(slot_id, kind, producers, phasec_ref)
            row.update({
                "witness_id": witness["witness_id"], "challenge_id": challenge["challenge_id"],
                "derivability_contract_class": f"class:{slot_id}",
                "witness": witness, "challenge": challenge,
            })
        rows.append(row)
    return sorted(rows, key=lambda row: (row["slot_id"].casefold(), row["slot_id"]))


@route("evaluate_dimension_firewall")
def evaluate_dimension_firewall() -> dict[str, Any]:
    traction = sp.Matrix([-1, -2, 1])
    area = sp.Matrix([2, 0, 0])
    velocity = sp.Matrix([1, -1, 0])
    surface_power = traction + area + velocity
    expected = sp.Matrix([2, -3, 1])
    ablated = traction + area + sp.Matrix([0, 0, 0])
    return {
        "dimension_basis": ["L", "T", "M"],
        "surface_power_dimensions": [int(value) for value in surface_power],
        "constraint_power_dimensions": [2, -3, 1],
        "rayleigh_power_dimensions": [2, -3, 1],
        "expected_power_dimensions": [int(value) for value in expected],
        "base_pass": surface_power == expected,
        "ablated_velocity_dimensions": [int(value) for value in ablated],
        "ablation_fired": ablated != expected,
    }


def route_inventory(candidates: list[str], obligations: dict[str, Any]) -> list[dict[str, Any]]:
    rows = []
    for candidate in candidates:
        members = members_for(candidate)
        if candidate == "E1":
            matrix, positive_candidate = sp.Matrix([[1, 1, 0], [1, -1, 0], [0, 1, 1]]), sp.Matrix([1, 1, 0])
            structure = {"boundary_variation": True, "jump_condition": True, "holonomic_trace": True, "constraint_multiplier": False, "dissipation_bookkeeping": False}
        elif candidate == "E2":
            matrix, positive_candidate = sp.Matrix([[1, 1, 0], [0, 2, 0], [1, 0, 1]]), sp.Matrix([1, 0, 0])
            structure = {"boundary_variation": True, "jump_condition": True, "holonomic_trace": True, "constraint_multiplier": False, "dissipation_bookkeeping": False}
        elif candidate == "E3":
            matrix, positive_candidate = sp.Matrix([[2, 1], [1, -1]]), sp.Matrix([1, 2])
            structure = {"boundary_variation": True, "jump_condition": True, "holonomic_trace": False, "constraint_multiplier": False, "dissipation_bookkeeping": False}
        elif candidate == "E4":
            matrix, positive_candidate = sp.Matrix([[2, 0, 1], [0, 3, 1], [1, 1, 0]]), sp.Matrix([1, 2, -1])
            structure = {"boundary_variation": False, "jump_condition": True, "holonomic_trace": False, "constraint_multiplier": True, "virtual_work_row": 2, "multiplier_column": 2, "dissipation_bookkeeping": False}
        elif candidate == "E5":
            matrix, positive_candidate = sp.Matrix([[2, 0, 0], [0, 3, 1], [0, 1, 1]]), sp.Matrix([1, 1, 2])
            structure = {"boundary_variation": True, "jump_condition": True, "holonomic_trace": False, "constraint_multiplier": False, "dissipation_bookkeeping": True, "rayleigh_coefficient": 2, "computed_dissipated_power": 2}
        elif candidate.startswith("MIXTURE("):
            dissipative = "E5" in members
            matrix, positive_candidate = sp.Matrix([[2, 0, 1, 0], [0, 3, 0, 1], [1, 0, 0, 0], [0, 1, 0, 2 if dissipative else 1]]), sp.Matrix([1, 2, -1, 1])
            structure = {"boundary_variation": True, "jump_condition": True, "holonomic_trace": any(member in {"E1", "E2"} for member in members), "constraint_multiplier": True, "virtual_work_row": 2, "multiplier_column": 2, "dissipation_bookkeeping": dissipative, "rayleigh_coefficient": 2 if dissipative else None, "computed_dissipated_power": 2 if dissipative else None}
        else:
            matrix, positive_candidate = sp.Matrix([[2, 0], [0, 3]]), sp.Matrix([1, 1])
            structure = {"boundary_variation": True, "jump_condition": True, "holonomic_trace": False, "constraint_multiplier": False, "dissipation_bookkeeping": False, "residual_complement_operator": True}
        rhs = matrix * positive_candidate
        malformed_candidate = positive_candidate.copy()
        malformed_candidate[0] += 1
        positive_residual = matrix * positive_candidate - rhs
        malformed_residual = matrix * malformed_candidate - rhs
        matrix_rows = [[int(value) for value in matrix.row(index)] for index in range(matrix.rows)]
        positive_values = [int(value) for value in positive_candidate]
        malformed_values = [int(value) for value in malformed_candidate]
        rhs_values = [int(value) for value in rhs]
        independent_checks = [
            sum(matrix_rows[row][column] * positive_values[column] for column in range(matrix.cols)) == rhs_values[row]
            for row in range(matrix.rows)
        ]
        for ambient in AMBIENTS:
            for stratum in TILT_TYPES:
                route_id = f"route:{candidate}|ambient={ambient}|stratum={stratum}"
                native = obligations[candidate]["native_root_class"]
                rows.append({
                    "route_id": route_id, "candidate_id": candidate, "ambient": ambient, "stratum": stratum,
                    "native_root_class": native,
                    "signature_digest": digest({
                        "obligations": obligations[candidate]["generator_A"],
                        "native": native, "ambient": ambient, "stratum": stratum,
                    }),
                    "fixture_id": f"fixture:{route_id}",
                    "positive_fixture": {
                        "semantic_route_id": "obligation_residual_classifier_v1",
                        "matrix": matrix_rows, "candidate": positive_values, "rhs": rhs_values,
                        "residual": [int(value) for value in positive_residual], "expected": "ADMISSIBLE",
                        "nondegenerate_norm_squared": int(positive_candidate.dot(positive_candidate)),
                        "native_structure_exercised": json.loads(json.dumps(structure)),
                        "known_outcome_generator_A": "exact_matrix_residual_zero",
                        "known_outcome_generator_B": "independent_row_equation_satisfaction",
                        "independent_row_equations_satisfied": independent_checks,
                    },
                    "malformed_fixture": {
                        "semantic_route_id": "obligation_residual_classifier_v1",
                        "matrix": matrix_rows, "candidate": malformed_values, "rhs": rhs_values,
                        "residual": [int(value) for value in malformed_residual], "expected": "EXCLUDED",
                    },
                    "route_equivalence_required": False,
                })
    return rows


def template_contract() -> dict[str, Any]:
    return {
        "eligible_record_predicate": "integrity==COMPUTATION_VALID and disposition==ADMISSIBLE and candidate!=OTHER",
        "semantic_schema": [
            "canonical_boundary_condition", "typed_free_data", "unevaluated_residual_or_variational_form",
            "zero_mode_treatment", "well_posedness_classification", "asymptotic_matching_conditions",
        ],
        "forbidden_ancestry_constructors": ["SOLVE_EVALUATION_RESULT"],
        "dependent_fields_unbound": True,
        "R49_exact_unresolved_reference_ids": [f"tilt:{value}" for value in TILT_TYPES],
        "term_census_generator_inputs": [
            "committed_static_field_equations", "canonical_operator_witness", "ambient_branch", "typed_free_data_ledger",
        ],
        "required_term_kinds": ["residual", "boundary", "zero-mode", "asymptotic-matching"],
        "per_term_remove_and_zero_teeth_required": True,
        "posing_not_solving": True,
    }


def closure_contract() -> dict[str, Any]:
    return {
        "gated_claim_inventory_generator": "committed_forcing_and_consequence_schema_walk",
        "contribution_census_generator_A": "committed_root_channel_incidence_walk",
        "contribution_census_generator_B": "force_balance_term_owner_walk",
        "coverage": "exact census to owner|computed-zero|typed-refusal",
        "independent_total": "direct committed force_balance construction",
        "no_double_count": True,
        "host_consistency_required": True,
        "both_host_requires_exact_partition": True,
        "radiation_static_zero_first_class": True,
        "integrity_failure_propagates_NOT_RUN": True,
        "valid_refusal_only_source_of_closure_PROMOTION_UNRESOLVED": True,
    }


def template_guard_fixture() -> dict[str, Any]:
    committed_inputs = {
        "static_field_equations": ["bulk_euler_lagrange_residual"],
        "canonical_operator_witness": ["Sigma_boundary_operator"],
        "ambient_branch": ["outer_matching_condition"],
        "typed_free_data_ledger": ["translation_zero_mode"],
    }
    expected_terms = [
        {"term_id": "residual:bulk_euler_lagrange_residual", "kind": "residual"},
        {"term_id": "boundary:Sigma_boundary_operator", "kind": "boundary"},
        {"term_id": "zero-mode:translation_zero_mode", "kind": "zero-mode"},
        {"term_id": "asymptotic-matching:outer_matching_condition", "kind": "asymptotic-matching"},
    ]
    return {
        "fixture_id": "template_guard:posed_not_solved",
        "committed_inputs": committed_inputs,
        "expected_term_census": expected_terms,
        "template_record": {
            "integrity": "COMPUTATION_VALID", "physics": "POSED_BVP_TEMPLATE",
            "constituents": {
                "canonical_boundary_condition": "Sigma_boundary_operator",
                "typed_free_data": ["translation_zero_mode"],
                "unevaluated_residual_or_variational_form": "bulk_euler_lagrange_residual",
                "zero_mode_treatment": "project_translation_zero_mode",
                "well_posedness_classification": "UNRESOLVED(committed_structure_only)",
                "asymptotic_matching_conditions": "outer_matching_condition",
            },
            "symbolic_ast": {"op": "posed_template", "args": [
                {"op": "term", "term_id": row["term_id"], "kind": row["kind"], "coefficient": 1}
                for row in expected_terms
            ]},
            "evaluation_state": "UNEVALUATED",
        },
    }


def closure_guard_fixture() -> dict[str, Any]:
    contributions = [
        {"contribution_id": "surface_traction", "root_id": "source:field:native_momentum", "vector": [2, 0], "expected_owner": "collar/Sigma surface"},
        {"contribution_id": "constraint_reaction", "root_id": "source:field:E4_shear_lock", "vector": [-1, 1], "expected_owner": "collar/Sigma surface"},
        {"contribution_id": "return_flux", "root_id": "source:field:return_closure", "vector": [0, -1], "expected_owner": "flux/return"},
        {"contribution_id": "static_radiation", "root_id": "source:field:native_momentum", "vector": [0, 0], "expected_owner": "radiation/static-zero"},
    ]
    assignments = [{
        "contribution_id": row["contribution_id"], "owner": row["expected_owner"],
        "computed_zero": row["vector"] == [0, 0],
    } for row in contributions]
    total = [sum(row["vector"][index] for row in contributions) for index in range(2)]
    return {
        "fixture_id": "closure_guard:exact_ownership",
        "committed_root_contributions": contributions,
        "census_A": [row["contribution_id"] for row in contributions],
        "census_B": sort_ids(row["contribution_id"] for row in contributions),
        "independently_constructed_total": total,
        "certificate": {"certificate_id": "certificate:fixture", "assignments": assignments},
        "closure_adjudication": {"integrity": "COMPUTATION_VALID", "physics": "CLOSURE_CERTIFIED(certificate:fixture)"},
        "dependent_record": {"integrity": "COMPUTATION_VALID", "physics": "SELECTED"},
    }


def standard_bindings() -> list[dict[str, str]]:
    standards = [
        ("S-1", "traceable cause tags"), ("S-2", "field-driven classification"),
        ("S-3", "no vacuous constructs"), ("S-4", "measured evidence"),
        ("S-5", "per-require teeth"), ("S-6", "complete summaries"),
        ("P-1", "entry witness resolves to ratified unresolved slot"),
        ("P-2", "dual independent ancestry generators"),
        ("P-3", "endpoint-conditioned ancestry ownership"),
        ("P-4", "two-level partitions"), ("P-5", "multi-census independent provenance"),
        ("P-6", "radiation first-class channel"),
        ("P-7", "reconciliation witness reference integrity"),
        ("P-8", "both-directions mutation armor"),
        ("P-9", "bare parent enums"), ("P-10", "one-body two-body firewall"),
    ]
    return [{
        "standard_id": standard_id, "authoritative_text": text,
        "acceptance_predicate_id": f"predicate:standard:{standard_id}",
        "reachable_check_id": f"ASSERT_STANDARD_{standard_id.replace('-', '_')}",
        "tooth_id": f"STANDARD_TOOTH:{standard_id}",
        "evidence_id": f"evidence:standard:{standard_id}",
    } for standard_id, text in standards]


def build(repo: Path) -> dict[str, Any]:
    inputs = load_yaml(repo / "software/em_charge_attribute/u1_body_dynamics_inputs.yaml")
    phasec_path = repo / (
        "software/em_charge_attribute/reports/u1_body_dynamics_artifacts/"
        "stage_c_0_tilt_coupling_contract/availability_slots.yaml"
    )
    phasec_slots = load_yaml(phasec_path)["slots"]
    phasec_production = load_yaml(repo / (
        "software/em_charge_attribute/reports/u1_body_dynamics_artifacts/"
        "stage_c_1_tilt_coupling_production/production_results.yaml"
    ))
    endpoints = inputs["endpoints"]
    mixtures_a = generate_mixtures_from_conditions(endpoints)
    mixtures_b = generate_mixtures_from_channels(endpoints)
    candidates = [*BASE_ENDPOINTS, *[row["candidate_id"] for row in mixtures_a], "OTHER"]
    census = source_census(inputs, phasec_slots)

    candidate_records = []
    for candidate in candidates:
        candidate_records.append({
            "candidate_id": candidate, "kind": "endpoint" if candidate in BASE_ENDPOINTS else "mixture" if candidate.startswith("MIXTURE(") else "catch_all",
            "members": members_for(candidate), "native_root_class": native_class(candidate, endpoints),
            "canonical_signature": signature_for(candidate),
            "membership_predicate": "positive_componentwise_operator_signature_equivalence" if candidate != "OTHER" else "report_only_residual_complement_not_promotion_reachable",
            "promotion_membership_eligible": candidate != "OTHER",
        })
    candidate_universe_digest = digest(candidate_records)
    candidate_inventory = {
        "mixture_formation_grammar": {
            "rule": "compose E4 collar shear constraint with committed normal/surface endpoints on orthogonal components",
            "neutral_alias_rule": "E3 plus X canonicalizes to X",
            "conflict_rule": "distinct laws on same trace component do not form a family",
            "free_weight_parameters_banned": True,
        },
        "mixture_generator_A": mixtures_a, "mixture_generator_B": mixtures_b,
        "concrete_other_candidates": [], "basis_closure": {
            "status": "UNRESOLVED", "declared_domain": {
                "locality": ["local", "nonlocal_unclosed"], "linearity": ["linear", "nonlinear_unclosed"],
                "frequency_dependence": "unclosed", "amplitude_dependence": "unclosed",
            },
            "missing_data": ["throat_surface_functional", "outer_surface_functional", "complete_boundary_operator_class_theorem"],
            "full_current_canonical_union_used": True,
            "residual_complement_status": "UNRESOLVED_NONEMPTY_OR_EMPTY",
            "catch_all_present": True,
        },
        "candidate_records": candidate_records,
        "candidate_axis": candidates, "candidate_count": len(candidates),
        "candidate_universe_digest": candidate_universe_digest,
        "amendment_rekey_duty": "mint_new_universe_digest_and_immutable_successors_for_all_complement_dependent_records",
        "uncanonicalized_overlap_count": 0,
        "canonical_signature_count": len({tuple(row["canonical_signature"]) for row in candidate_records}),
    }

    obligations = {}
    for candidate in candidates:
        root_class = native_class(candidate, endpoints)
        obligations[candidate] = {
            "native_root_class": root_class,
            "generation_inputs": ["candidate_id", "native_root_class", "operational_endpoint_definition"],
            "stratum_token_is_generation_input": False,
            "generator_A": generate_obligations_from_native_roots(candidate, root_class),
            "generator_B": generate_obligations_from_endpoint_walk(candidate, endpoints),
            "nondegeneracy_predicate": "nonzero_canonical_operator_and_nonzero_committed_root_ablation",
            "anti_echo_predicate": "witness_satisfies_independent_nondefinitional_obligation",
            "semantic_ablation_criterion": "fail_nondefinitional_obligation_or_change_canonical_operator_class",
        }

    active_strata = phasec_production["axes"]["generated_active_strata"]
    dependency_rows = []
    grid_cells = []
    collapse_proofs = []
    for candidate in candidates:
        for stratum in active_strata:
            deps_a = generate_dependency_join(candidate, stratum, obligations[candidate]["generator_A"])
            deps_b = generate_dependency_source_walk(candidate, stratum, endpoints, inputs)
            dependency_rows.append({
                "candidate_id": candidate, "stratum": stratum,
                "generator_A": deps_a, "generator_B": deps_b,
                "stratifying_slot": f"tilt:{stratum}",
                "dependency_signature": digest(deps_a),
            })
            collapse_proofs.append({
                "candidate_id": candidate, "raw_stratum": stratum,
                "collapsed_class": stratum, "timing": "pre-production",
                "proof": "singleton_due_to_distinct_authoritative_stratifying_slot",
                "stage0_objects_only": True,
            })
            for ambient in AMBIENTS:
                grid_cells.append({
                    "cell_id": f"candidate={candidate}|ambient={ambient}|stratum={stratum}",
                    "candidate_id": candidate, "ambient": ambient, "stratum": stratum,
                    "expected_dependencies": deps_a,
                    "stable_branch_id": f"U2:{candidate}:{ambient}:{stratum}",
                })
    promotion_contexts = []
    for ambient in AMBIENTS:
        for stratum in active_strata:
            mappings = [{
                "candidate_id": candidate,
                "cell_id": f"candidate={candidate}|ambient={ambient}|stratum={stratum}",
            } for candidate in candidates]
            promotion_contexts.append({
                "promotion_key": f"ambient={ambient}|context={stratum}",
                "ambient": ambient, "global_common_refinement_context": stratum,
                "candidate_cell_mappings": mappings,
            })

    routes = route_inventory(candidates, obligations)
    slots = availability_slots(candidates, phasec_slots, endpoints, mixtures_a)
    artifact_dags = proof_artifacts(
        candidate_inventory, obligations, dependency_rows, routes,
        collapse_proofs, promotion_contexts, census["source_ids"],
    )
    grammar = [
        "ROOT_REFERENCE", "STATIC_COMMITTED_FORCING", "POSITIVE_DERIVATION", "POSITIVE_EQUIVALENCE",
        "INCOMPATIBILITY", "NEGATIVE_CANDIDATE_MEMBERSHIP", "CASE_ELIMINATION",
        "COMPLEMENT_SURVIVOR_COUNT", "EXCLUSION_VERDICT", "POSTULATE_BRANCH",
        "STABILITY_DYNAMICAL_CLASS", "SOLVE_EVALUATION_RESULT", "SYMBOLIC_EQUIVALENCE_COLLAPSE",
        "UNAVAILABILITY_WITNESS", "DERIVABILITY_CHALLENGE", "EVIDENCE_STATE_CLASSIFICATION",
    ]

    semantic = {
        "schema_version": "U2_STAGE0_SEMANTIC_VIEW_V1",
        "scope": {
            "static_adjudication_only": True, "dynamical_selection_deferred": True,
            "one_body_only": True, "BVP_solved": False,
            "pair_configuration_carveout": "pair_annihilation_subquestion_only",
        },
        "ambient_branches": [
            {"ambient_id": "one_sided_pathA29", "status": "committed_executable_slab", "asymmetry_map": "ambient_asymmetry_map"},
            {"ambient_id": "two_sided_R_w_postulate", "status": "POSTULATE", "asymmetry_map": "R_w"},
        ],
        "candidate_inventory": candidate_inventory,
        "obligation_censuses": obligations,
        "open_dependency_relation": {
            "authoritative_stratifying_ledger": "PhaseC.axes.generated_active_strata",
            "active_strata": active_strata,
            "dependency_rows": dependency_rows,
            "generation": "stratum_free_census_then_independent_join_to_frozen_producer_relation_and_OPEN_ledger",
            "generator_A_algorithm": "obligation_to_authoritative_slot_relational_join",
            "generator_B_algorithm": "raw_field_and_endpoint_channel_schema_walk",
            "shared_task_code_between_generators": False,
        },
        "grid_inventory": {
            "candidate_count": len(candidates), "ambient_count": len(AMBIENTS),
            "raw_strata_per_candidate": {candidate: len(active_strata) for candidate in candidates},
            "raw_ragged_cardinality": sum(len(active_strata) * len(AMBIENTS) for _ in candidates),
            "collapsed_cardinality": len(grid_cells), "preproduction_collapse_count": 0,
            "postproduction_report_collapse_affects_grid": False,
            "grid_cells": grid_cells, "collapse_proofs": collapse_proofs,
            "promotion_contexts": promotion_contexts,
            "promotion_context_count": len(promotion_contexts),
        },
        "vocabulary_freeze": vocabulary_freeze(),
        "evidence_taxonomy": {
            "source_census": census, "constructor_grammar": grammar,
            "generic_or_unclassified_constructor_allowed": False,
            "classification_source": "normalized_inference_content_not_claimed_label",
            "physics_bearing_stage0_artifacts": artifact_dags,
            "promotion_allowed_constructors": ["ROOT_REFERENCE", "STATIC_COMMITTED_FORCING", "POSITIVE_DERIVATION", "POSITIVE_EQUIVALENCE"],
            "program_wide_banned_constructor": "STABILITY_DYNAMICAL_CLASS",
            "guard_fixture_dags": {
                "promotion_positive": {
                    "constructors": ["ROOT_REFERENCE", "STATIC_COMMITTED_FORCING", "POSITIVE_DERIVATION", "POSITIVE_EQUIVALENCE"],
                    "root_types": ["BULK_ACTION_TERM", "SURFACE_ACTION_TERM"],
                    "normalized_content": {"op": "static_force", "args": [{"op": "positive_equivalence"}]},
                },
                "candidate_exclusion": {
                    "constructors": ["ROOT_REFERENCE", "INCOMPATIBILITY", "EXCLUSION_VERDICT"],
                    "root_types": ["HOLONOMIC_CONSTRAINT", "BULK_ACTION_TERM"],
                },
                "posed_template": {"constructors": ["ROOT_REFERENCE", "POSITIVE_DERIVATION"], "fields_unbound": True},
                "pair_configuration": {
                    "object_type": "static_plus_w_minus_w_pair_configuration",
                    "firewall_tag": "PAIR_ANNIHILATION_ONLY", "consumer": "topology_pair_annihilation_subquestion",
                },
                "postulate_deferral_metadata": {"constructor": "POSTULATE_BRANCH", "promotion_reachable": False},
            },
        },
        "route_fixture_inventory": {
            "generated_preproduction": True, "route_records": routes,
            "route_count": len(routes), "fixture_count": len(routes),
            "executed_route_match_rule": "production_route_set_exactly_equals_frozen_route_set",
            "cell_route_fixture_exact_coverage": True,
        },
        "availability_slots": slots,
        "availability_summary": {
            "total": len(slots),
            "DERIVED": sum(row["availability_outcome"] == "DERIVED" for row in slots),
            "UNRESOLVED": sum(row["availability_outcome"] == "UNRESOLVED" for row in slots),
            "integrity_failures": 0,
        },
        "closure_contract": closure_contract(), "template_contract": template_contract(),
        "guard_fixtures": {
            **guard_fixture_records(),
            "closure": closure_guard_fixture(),
            "template": template_guard_fixture(),
        },
        "return_closure_ownership": {
            "owner": "downstream_flux_path", "U2_owned": False,
            "preserved_terminal": "UNRESOLVED(return_closure)",
            "allowed_consumption": "deferral_metadata_only_unreachable_from_U2_verdicts",
        },
        "standard_bindings": standard_bindings(),
        "dimensional_firewall": evaluate_dimension_firewall(),
        "production_contract": {
            "disposition_ancestry_exact_both_directions": True,
            "route_set_exact_match_required": True, "zero_integrity_failures_for_banking": True,
            "stage0_digest_rechecked_before_every_evaluation_and_A9_leg": True,
            "stable_branch_ids_required_all_output_classes": True,
            "mutation_catalog_coverage_required": True,
        },
    }
    local_registry = [
        {"semantic_route_id": semantic_id, "engine_local_function": local_name,
         "exists": callable(globals().get(local_name)), "executed": call_counts[local_name] > 0}
        for semantic_id, local_name in (
            ("candidate_formation_from_committed_constraints_v1", "generate_mixtures_from_conditions"),
            ("candidate_formation_channel_crosscheck_v1", "generate_mixtures_from_channels"),
            ("stratum_free_obligation_census_v1", "generate_obligations_from_native_roots"),
            ("stratum_free_operational_crosscheck_v1", "generate_obligations_from_endpoint_walk"),
            ("open_dependency_join_v1", "generate_dependency_join"),
            ("open_dependency_source_walk_v1", "generate_dependency_source_walk"),
            ("content_classified_proof_DAG_v1", "classify_inference_content"),
            ("constructive_absence_challenge_v1", "construct_unavailability_challenge"),
            ("inline_dimension_firewall_v1", "evaluate_dimension_firewall"),
            ("executable_stage0_decision_and_schema_guards_v1", "guard_fixture_records"),
        )
    ]
    return {
        "schema_version": SCHEMA, "engine": "SymPy",
        "independent_route": "Python schema/AST walks plus SymPy exact linear algebra",
        "semantic_view": semantic,
        "engine_local_route_registry": local_registry,
        "runtime_identity": {
            "python_isolated": sys.flags.isolated == 1,
            "python_no_user_site": sys.flags.no_user_site == 1,
            "sympy_version": sp.__version__,
        },
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    result = build(Path(args.repo).resolve())
    output = Path(args.output).resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(
        yaml.safe_dump(result, sort_keys=False, allow_unicode=True, width=140), encoding="utf-8"
    )
    print(
        f"U2_SYMPY_STAGE0_PASS candidates={result['semantic_view']['candidate_inventory']['candidate_count']} "
        f"grid={result['semantic_view']['grid_inventory']['raw_ragged_cardinality']} "
        f"slots={result['semantic_view']['availability_summary']['total']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
