#!/usr/bin/env python3
"""Independent comparator for the U1 Phase-C stage-0 engines."""

from __future__ import annotations

import argparse
import ast
import copy
import hashlib
import json
import re
import sys
from collections import Counter
from pathlib import Path
from typing import Any

import yaml


SCHEMA = "U1_PHASE_C_STAGE0_ENGINE_V1"
EXPECTED_B1_TILT = {
    "indexed_density_tilt_profile",
    "indexed_flow_tilt_response",
    "indexed_h_tilt_profile",
    "indexed_phase_tilt_profile",
    "indexed_shear_tilt_profile",
    "indexed_sleeve_surface_normal_profile",
    "indexed_sleeve_tilt_profile",
    "indexed_uw_tilt_profile",
}
EXPECTED_FLOOR = {
    "common_drain",
    "orientation_sleeve",
    "post_mouth_axial_flow",
    "h",
    "u_T",
    "u_L",
    "E4_constraint_reaction",
    "outer_surface_flux_return",
    "wall_chi_modes",
}
EXPECTED_RADIATION_ROOTS = {
    "radiation:brane_normal_local",
    "radiation:brane_shear_transverse",
    "radiation:geon_open",
    "radiation:gnls_density_phase",
    "radiation:h_scalar",
    "radiation:throat_source_open",
    "radiation:wall_chi_u_coupled",
    "radiation:wall_mix_open",
}
ALLOWED_KINDS = {
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


def sha256_path(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def canonical_bytes(value: Any) -> bytes:
    return json.dumps(
        value, ensure_ascii=False, sort_keys=True, separators=(",", ":")
    ).encode("utf-8")


def digest(value: Any) -> str:
    return hashlib.sha256(canonical_bytes(value)).hexdigest()


def string_record_digest(records: list[str]) -> str:
    row_hashes = [hashlib.sha256(record.encode("utf-8")).hexdigest() for record in records]
    return hashlib.sha256("\n".join(sorted(row_hashes)).encode("ascii")).hexdigest()


def semantic_sink_records(engine: dict[str, Any]) -> list[str]:
    derivation = engine["native_derivation"]
    force = engine["force_term_census"]
    coupling = engine["coupling_source_census"]
    g8 = engine["g8_ablation_inventory"]
    projection = engine["projection_freeze"]
    return (
        [f"derivation_term:{row['id']}" for row in derivation["second_variation"]["term_records"]]
        + [
            f"hessian_pair:{row['term']}|{row['pair']}"
            for row in derivation["second_variation"]["nonzero_pair_records"]
        ]
        + [f"force:{value}" for value in force["expected_terms"]]
        + [f"coupling:{value}" for value in coupling["expected_entries"]]
        + [f"g8:{row['source_id']}" for row in g8["entries"]]
        + [
            f"slot:{row['slot_id']}|{row['disposition']}"
            for row in engine["availability_slots"]
        ]
        + [
            f"witness_measurement:{row['slot_id']}|"
            f"{row['witness']['insufficiency_certificate']['measurement_digest']}"
            for row in engine["availability_slots"]
            if row["disposition"] == "UNRESOLVED"
        ]
        + [
            f"challenge_measurement:{row['slot_id']}|{row['challenge']['measurement_digest']}"
            for row in engine["availability_slots"]
            if row["disposition"] == "UNRESOLVED"
        ]
        + [
            "hessian_constructive_challenge:"
            + engine["hessian_constructive_challenge"]["candidate_value_digest"]
        ]
        + [f"projection:{value}" for value in projection["predicates"]]
        + [
            f"reconciliation:{value}"
            for value in engine["frozen_assertions"]["reconciliation_expected_ids"]
        ]
    )


def by_key(rows: list[dict[str, Any]], key: str) -> dict[str, dict[str, Any]]:
    return {row[key]: row for row in rows}


def strip_engine_local(value: Any) -> Any:
    if isinstance(value, dict):
        out = {}
        for key, child in value.items():
            if key == "engine_local_function":
                continue
            out[key] = strip_engine_local(child)
        return out
    if isinstance(value, list):
        return [strip_engine_local(child) for child in value]
    return value


def normalized_derivation(derivation: dict[str, Any]) -> dict[str, Any]:
    s_body = strip_engine_local(derivation["S_body_Omega_c"])
    return {
        "S_body_Omega_c": s_body,
        "second_variation": {
            "status": derivation["second_variation"]["status"],
            "term_records": derivation["second_variation"]["term_records"],
            "nonzero_pair_records": derivation["second_variation"][
                "nonzero_pair_records"
            ],
            "chi_u_mixed_block_present": derivation["second_variation"][
                "chi_u_mixed_block_present"
            ],
            "value_digest": derivation["second_variation"]["value_digest"],
        },
    }


def normalized_census(census: dict[str, Any], kind: str) -> dict[str, Any]:
    value = strip_engine_local(census)
    if kind == "force":
        value["entries"] = sorted(value["entries"], key=lambda row: row["term_id"])
        value["root_incidence"] = sorted(
            value["root_incidence"], key=lambda row: row["root_id"]
        )
        value["expected_terms"] = sorted(value["expected_terms"])
        value["reachable_residual_terms"] = sorted(value["reachable_residual_terms"])
    elif kind == "coupling":
        value["sources"] = sorted(value["sources"], key=lambda row: row["source_id"])
        value["entries"] = sorted(value["entries"], key=lambda row: row["entry_id"])
        value["expected_entries"] = sorted(value["expected_entries"])
        value["reachable_entries"] = sorted(value["reachable_entries"])
        value["ordered_deltaO_entries"] = sorted(value["ordered_deltaO_entries"])
        value["coverage_checks"]["source_ids"] = sorted(
            value["coverage_checks"]["source_ids"]
        )
    else:
        value["entries"] = sorted(value["entries"], key=lambda row: row["source_id"])
        value["floor_coverage"] = {
            key: sorted(rows) for key, rows in sorted(value["floor_coverage"].items())
        }
    return value


def normalized_slots(rows: list[dict[str, Any]]) -> dict[str, dict[str, Any]]:
    return {
        value["slot_id"]: value
        for value in (strip_engine_local(row) for row in rows)
    }


def engine_route_certificates_valid(
    sympy: dict[str, Any], wolfram: dict[str, Any], repo: Path
) -> bool:
    """Verify shared semantics and each certificate's executed local symbol.

    Engine-local function names are intentionally excluded from cross-engine
    record equality, but they are not trusted as prose: each name is pinned to
    the expected engine, found in that engine's source, and carried only on an
    certificate emitted by the named function.  The route-specific witness,
    challenge, and Hessian assertions independently verify ``executed=True``
    so their mutation teeth still die at their own named assertions.
    """
    python_source = (
        repo / "software/em_charge_attribute/u1_phaseC_stage0_sympy.py"
    ).read_text(encoding="utf-8")
    wolfram_source = (
        repo / "software/em_charge_attribute/u1_phaseC_stage0_dual.wl"
    ).read_text(encoding="utf-8")
    python_functions = {
        node.name
        for node in ast.walk(ast.parse(python_source))
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef))
    }
    expected = {
        "witness": (
            "typed_witness_insufficiency_measurement_v1",
            "witness_insufficiency_measurement",
            "witnessInsufficiencyMeasurement",
        ),
        "challenge": (
            "constructive_derivation_challenge_v1",
            "constructive_derivation_challenge",
            "constructiveDerivationChallenge",
        ),
        "hessian": (
            "raw_action_hessian_challenge_v1",
            "hessian_challenge_from_raw_action",
            "hessianChallengeFromRawAction",
        ),
    }
    if not all(spec[1] in python_functions for spec in expected.values()):
        return False
    for _, _, local_name in expected.values():
        definition = rf"(?ms)^\s*{re.escape(local_name)}\[.*?\]\s*:=\s*Module\["
        if re.search(definition, wolfram_source) is None:
            return False

    sympy_rows = by_key(sympy["availability_slots"], "slot_id")
    wolfram_rows = by_key(wolfram["availability_slots"], "slot_id")
    if set(sympy_rows) != set(wolfram_rows):
        return False
    witness_spec = expected["witness"]
    challenge_spec = expected["challenge"]
    certified_slot_count = 0
    for slot_id, row_s in sympy_rows.items():
        row_w = wolfram_rows[slot_id]
        if row_s["disposition"] != "UNRESOLVED":
            continue
        if not all(key in row_s and key in row_w for key in ("witness", "challenge")):
            # The Hessian anti-dodge mutation temporarily turns a DERIVED slot
            # into a special unresolved candidate.  Its own guard validates
            # that record; it is not one of the 56 witness certificates.
            continue
        certified_slot_count += 1
        certificate_pairs = (
            (
                row_s["witness"]["insufficiency_certificate"]["engine_certificate"],
                row_w["witness"]["insufficiency_certificate"]["engine_certificate"],
                witness_spec,
            ),
            (
                row_s["witness"]["counterfactual_restore_mutation"][
                    "restored_certificate"
                ]["engine_certificate"],
                row_w["witness"]["counterfactual_restore_mutation"][
                    "restored_certificate"
                ]["engine_certificate"],
                witness_spec,
            ),
            (
                row_s["challenge"]["engine_certificate"],
                row_w["challenge"]["engine_certificate"],
                challenge_spec,
            ),
        )
        for certificate_s, certificate_w, (semantic, local_s, local_w) in certificate_pairs:
            if not (
                certificate_s.get("semantic_route_id")
                == certificate_w.get("semantic_route_id")
                == semantic
                and certificate_s.get("engine_local_function") == local_s
                and certificate_w.get("engine_local_function") == local_w
            ):
                return False

    hessian_s = sympy["hessian_constructive_challenge"]
    hessian_w = wolfram["hessian_constructive_challenge"]
    semantic, local_s, local_w = expected["hessian"]
    return (
        certified_slot_count
        == sympy["availability_summary"]["UNRESOLVED"]
        == wolfram["availability_summary"]["UNRESOLVED"]
        == 56
        and hessian_s.get("semantic_route_id")
        == hessian_w.get("semantic_route_id")
        == semantic
        and hessian_s.get("engine_local_function") == local_s
        and hessian_w.get("engine_local_function") == local_w
    )


def slot_schema(row: dict[str, Any]) -> dict[str, Any]:
    return {
        "required_type": row["required_type"],
        "required_dimensions": row["required_dimensions"],
        "domain": row["domain"],
    }


def one_column_rank(matrix: list[list[int]]) -> int:
    return int(any(cell for row in matrix for cell in row))


def expected_witness_status(
    certificate: dict[str, Any], kind: str
) -> tuple[str, bool]:
    measurements = certificate["producer_measurements"]
    matrix = [[int(row["acceptance_selecting_count"] > 0)] for row in measurements]
    rank = one_column_rank(matrix)
    nullity = 1 - rank
    selecting = sum(row["acceptance_selecting_count"] for row in measurements)
    domains = sum(row["domain_completion_count"] for row in measurements)
    if kind == "nonuniqueness/solvability failure":
        passes = nullity > 0
        predicate = "constraint_matrix_nullity_gt_zero"
    elif kind == "absence of any typed producer in the complete authority census":
        passes = selecting == 0
        predicate = "universal_typed_producer_elimination"
    elif kind == "operator/domain well-posedness failure":
        passes = domains == 0
        predicate = "no_closed_operator_domain_or_BC_completion"
    else:
        passes = selecting == 0
        predicate = "all_candidate_dimensions_incompatible"
    expected = "PASS_COMPUTED" if passes else "FAIL_COMPUTED"
    evidence_digests_valid = all(
        row["evidence_digest"] == digest(row["evidence_ids"]) for row in measurements
    )
    fields_valid = (
        certificate["executed"] is True
        and certificate["predicate"] == predicate
        and certificate["candidate_count"] == len(measurements)
        and certificate["constraint_matrix"] == matrix
        and certificate["measured_rank"] == rank
        and certificate["measured_nullity"] == nullity
        and certificate["compatible_selecting_producer_count"] == selecting
        and certificate["domain_completion_count"] == domains
        and certificate["measurement_digest"] == digest(measurements)
        and certificate["status"] == expected
        and certificate["engine_certificate"]["executed"] is True
        and certificate["engine_certificate"]["algorithm"]
        == "matrix_rank_plus_typed_authority_scan"
        and evidence_digests_valid
    )
    return expected, fields_valid


def expected_challenge_result(
    row: dict[str, Any]
) -> tuple[str, bool, bool]:
    challenge = row["challenge"]
    certificate = challenge["engine_certificate"]
    scans = certificate["raw_producer_scan"]
    matrix = [[int(scan["raw_selecting_count"] > 0)] for scan in scans]
    rank = one_column_rank(matrix)
    nullity = 1 - rank
    selecting = sum(scan["raw_selecting_count"] for scan in scans)
    domains = sum(scan["raw_domain_completion_count"] for scan in scans)
    kind = challenge["kind"]
    if kind == "nonuniqueness/solvability failure":
        result = "FAIL_NOT_UNIQUE" if nullity > 0 else "PASS"
    elif kind == "absence of any typed producer in the complete authority census":
        result = "FAIL_NO_APPROVED_INSTANTIATED_VALUE" if selecting == 0 else "PASS"
    elif kind == "operator/domain well-posedness failure":
        result = "FAIL_DOMAIN_NOT_CLOSED" if domains == 0 else "PASS"
    else:
        result = "FAIL_DIMENSIONAL_COMPATIBILITY" if selecting == 0 else "PASS"
    candidate = challenge["attempted_candidates"][0]
    schema_matches = (
        challenge["candidate_schema_pinned"] == slot_schema(row)
        and candidate["candidate_schema"] == challenge["candidate_schema_pinned"]
    )
    candidate_passes = schema_matches and result == "PASS"
    outcome = "REFUTED" if candidate_passes else "CONSTRUCTIVE_FAIL"
    fields_valid = (
        certificate["executed"] is True
        and certificate["algorithm"]
        == "independent_raw_authority_walk_plus_constraint_rank"
        and certificate["attempted_candidate_count"]
        == len(challenge["attempted_candidates"])
        == 1
        and certificate["empty_output"] is False
        and certificate["ill_typed_by_fiat"] is False
        and certificate["constraint_matrix"] == matrix
        and certificate["measured_rank"] == rank
        and certificate["measured_nullity"] == nullity
        and challenge["measurement_digest"] == digest(scans)
        and candidate["candidate_is_well_typed"] is schema_matches
        and candidate["defining_predicate_evaluated"] is True
        and candidate["defining_predicate_result"] == result
        and candidate["passes_defining_predicate"] is candidate_passes
        and challenge["outcome"] == outcome
    )
    return result, candidate_passes, fields_valid


def normalized_hessian_challenge(engine: dict[str, Any]) -> dict[str, Any]:
    return strip_engine_local(engine["hessian_constructive_challenge"])


def hessian_candidate_passes(engine: dict[str, Any], row: dict[str, Any]) -> bool:
    challenge = engine["hessian_constructive_challenge"]
    candidate = copy.deepcopy(challenge["candidate_value"])
    derived = normalized_derivation(engine["native_derivation"])["second_variation"]
    return (
        row["disposition"] == "UNRESOLVED"
        and row.get("challenge", {}).get("challenge_kind")
        == "hessian_constructive_candidate"
        and challenge["executed"] is True
        and challenge["candidate_schema"] == slot_schema(row)
        and challenge["candidate_is_well_typed"]
        == (challenge["candidate_schema"] == slot_schema(row))
        and challenge["defining_predicate_evaluated"] is True
        and challenge["defining_predicate_pass"]
        == (
            candidate["chi_u_mixed_block_present"]
            and len(candidate["term_records"])
            == len(derived["term_records"])
        )
        and challenge["candidate_value_digest"] == digest(challenge["candidate_value"])
        and candidate == derived
    )


def enforce_unresolved_guard(engine: dict[str, Any]) -> None:
    for row in engine["availability_slots"]:
        if row["disposition"] != "UNRESOLVED":
            continue
        challenge = row.get("challenge", {})
        if challenge.get("challenge_kind") == "hessian_constructive_candidate":
            passing = hessian_candidate_passes(engine, row)
        else:
            try:
                _, passing, fields_valid = expected_challenge_result(row)
                passing = passing and fields_valid
            except (KeyError, IndexError, TypeError):
                passing = False
        if passing:
            raise AssertionDeath("stage0_unresolved_refuted", row["slot_id"])


def plant_restored_witness(engine: dict[str, Any], slot_id: str) -> None:
    row = by_key(engine["availability_slots"], "slot_id")[slot_id]
    restored = copy.deepcopy(
        row["witness"]["counterfactual_restore_mutation"]["restored_certificate"]
    )
    row["witness"]["insufficiency_certificate"] = restored


def plant_refuting_candidate(engine: dict[str, Any], slot_id: str) -> None:
    row = by_key(engine["availability_slots"], "slot_id")[slot_id]
    challenge = row["challenge"]
    fixture = row["witness"]["counterfactual_restore_mutation"][
        "fixture_producer_measurement"
    ]
    scans = copy.deepcopy(challenge["engine_certificate"]["raw_producer_scan"])
    scans.append(
        {
            "producer_id": fixture["producer_id"],
            "raw_record_count": fixture["authority_record_count"],
            "raw_instantiated_count": fixture["instantiated_value_count"],
            "raw_domain_completion_count": fixture["domain_completion_count"],
            "raw_selecting_count": fixture["acceptance_selecting_count"],
        }
    )
    certificate = challenge["engine_certificate"]
    certificate["raw_producer_scan"] = scans
    matrix = [[int(scan["raw_selecting_count"] > 0)] for scan in scans]
    certificate["constraint_matrix"] = matrix
    certificate["measured_rank"] = one_column_rank(matrix)
    certificate["measured_nullity"] = 1 - certificate["measured_rank"]
    challenge["measurement_digest"] = digest(scans)
    candidate = challenge["attempted_candidates"][0]
    candidate["candidate_schema"] = copy.deepcopy(challenge["candidate_schema_pinned"])
    _, _, _ = expected_challenge_result(row)
    kind = challenge["kind"]
    selecting = sum(scan["raw_selecting_count"] for scan in scans)
    domains = sum(scan["raw_domain_completion_count"] for scan in scans)
    nullity = certificate["measured_nullity"]
    if kind == "nonuniqueness/solvability failure":
        result = "FAIL_NOT_UNIQUE" if nullity > 0 else "PASS"
    elif kind == "absence of any typed producer in the complete authority census":
        result = "FAIL_NO_APPROVED_INSTANTIATED_VALUE" if selecting == 0 else "PASS"
    elif kind == "operator/domain well-posedness failure":
        result = "FAIL_DOMAIN_NOT_CLOSED" if domains == 0 else "PASS"
    else:
        result = "FAIL_DIMENSIONAL_COMPATIBILITY" if selecting == 0 else "PASS"
    schema_matches = candidate["candidate_schema"] == slot_schema(row)
    candidate["candidate_is_well_typed"] = schema_matches
    candidate["defining_predicate_evaluated"] = True
    candidate["defining_predicate_result"] = result
    candidate["passes_defining_predicate"] = schema_matches and result == "PASS"
    challenge["outcome"] = (
        "REFUTED" if candidate["passes_defining_predicate"] else "CONSTRUCTIVE_FAIL"
    )
    challenge["mutation_evidence"] = {
        "source": "executed_counterfactual_restore_measurement",
        "fixture_evidence_digest": fixture["evidence_digest"],
    }


def plant_hessian_disposition_defect(engine: dict[str, Any], first_time: bool) -> None:
    slot_id = "support:complete_action_second_variation"
    rows = engine["availability_slots"]
    if first_time:
        engine["availability_slots"] = [row for row in rows if row["slot_id"] != slot_id]
        challenge = engine["hessian_constructive_challenge"]
        schema = challenge["candidate_schema"]
        row = {
            "slot_id": slot_id,
            "category": "derived_support_core",
            "required_type": schema["required_type"],
            "required_dimensions": schema["required_dimensions"],
            "domain": schema["domain"],
            "producer_set": ["action:*"],
            "acceptance_predicate": "ALL_ACTION_TERMS_INCIDENT_AND_CHI_U_MIXED_BLOCK_PRESENT",
        }
        engine["availability_slots"].append(row)
    else:
        row = by_key(rows, "slot_id")[slot_id]
        for key in ("value", "value_digest", "dual_engine_comparison_id"):
            row.pop(key, None)
    row["disposition"] = "UNRESOLVED"
    row["challenge"] = {
        "challenge_kind": "hessian_constructive_candidate",
        "candidate_ref": "hessian_constructive_challenge",
        "first_time_construction": first_time,
    }


def apply_mutation(sympy: dict[str, Any], wolfram: dict[str, Any], mutation: str | None) -> None:
    if not mutation:
        return
    s = sympy
    w = wolfram
    if mutation.startswith("ASSERT_WITNESS_INSUFFICIENCY:"):
        slot_id = mutation.split(":", 1)[1]
        plant_restored_witness(s, slot_id)
        plant_restored_witness(w, slot_id)
    elif mutation == "ASSERT_ENGINE_SCHEMA":
        s["schema_version"] = "MUTATED"
    elif mutation == "ASSERT_ENGINE_INDEPENDENCE":
        w["independent_route"] = s["independent_route"]
    elif mutation == "ASSERT_ENGINE_LOCAL_ROUTE_CERTIFICATES":
        by_key(s["availability_slots"], "slot_id")["J:h"]["challenge"][
            "engine_certificate"
        ]["engine_local_function"] = "nonexistent_route"
    elif mutation == "ASSERT_SOURCE_DIGESTS_LIVE":
        s["source_digests"]["phase_a_inputs"] = "0" * 64
        w["source_digests"]["phase_a_inputs"] = "0" * 64
    elif mutation == "ASSERT_B1_TILT_INDEX":
        s["frozen_assertions"]["b1_tilt_ingredient_types"].pop()
        w["frozen_assertions"]["b1_tilt_ingredient_types"].pop()
    elif mutation == "ASSERT_B1_PARTITION_TERMINAL":
        s["frozen_assertions"]["b1_partition_state"] = "closed"
        w["frozen_assertions"]["b1_partition_state"] = "closed"
    elif mutation == "ASSERT_B2_PARTITION_TERMINAL":
        s["frozen_assertions"]["b2_partition_state_sympy"] = "closed"
        w["frozen_assertions"]["b2_partition_state_sympy"] = "closed"
    elif mutation == "ASSERT_B2_DISPOSITION_PIN":
        s["frozen_assertions"]["b2_stage0_disposition_sha256"] = "0" * 64
        w["frozen_assertions"]["b2_stage0_disposition_sha256"] = "0" * 64
    elif mutation == "ASSERT_B2_TILT_PATH_EQUALITY":
        s["frozen_assertions"]["b2_tilt_paths_mathematica"].pop()
        w["frozen_assertions"]["b2_tilt_paths_mathematica"].pop()
    elif mutation == "ASSERT_RECONCILIATION_EXACT":
        s["frozen_assertions"]["reconciliation_expected_ids"].pop()
    elif mutation == "ASSERT_RECONCILIATION_CARDINALITY":
        fake = "B2_DEFERRED|invented_obligation"
        s["frozen_assertions"]["reconciliation_expected_ids"].append(fake)
        w["frozen_assertions"]["reconciliation_expected_ids"].append(fake)
    elif mutation == "ASSERT_ACTION_ROOT_SET":
        s["typed_root_inventory"].pop()
        w["typed_root_inventory"].pop()
    elif mutation == "ASSERT_SECOND_VARIATION_AGREE":
        s["native_derivation"]["second_variation"]["term_records"][0][
            "nonzero_hessian_pairs"
        ] = []
    elif mutation == "ASSERT_HESSIAN_CHALLENGE_EXECUTED":
        s["hessian_constructive_challenge"]["executed"] = False
        w["hessian_constructive_challenge"]["executed"] = False
    elif mutation == "ASSERT_HESSIAN_CHALLENGE_AGREE":
        s["hessian_constructive_challenge"]["candidate_value"]["term_records"].pop()
    elif mutation == "ASSERT_CHI_U_MIXED_BLOCK":
        s["native_derivation"]["second_variation"]["chi_u_mixed_block_present"] = False
        w["native_derivation"]["second_variation"]["chi_u_mixed_block_present"] = False
    elif mutation == "ASSERT_S_BODY_AGREE":
        s["native_derivation"]["S_body_Omega_c"]["canonical_terms"].pop()
    elif mutation == "ASSERT_FORCE_CENSUS_AGREE":
        s["force_term_census"]["entries"][0]["formal_expression"] = "MUTATED"
    elif mutation == "ASSERT_FORCE_INCIDENCE":
        s["force_term_census"]["root_incidence"][0]["incidence_complete"] = False
        w["force_term_census"]["root_incidence"][0]["incidence_complete"] = False
        s["force_term_census"]["coverage_checks"]["force_census_incidence_complete"] = False
        w["force_term_census"]["coverage_checks"]["force_census_incidence_complete"] = False
    elif mutation == "ASSERT_FORCE_EXPECTED_REACHABLE":
        s["force_term_census"]["reachable_residual_terms"].pop()
        w["force_term_census"]["reachable_residual_terms"].pop()
        s["force_term_census"]["coverage_checks"]["expected_reachable_exact_set_equal"] = False
        w["force_term_census"]["coverage_checks"]["expected_reachable_exact_set_equal"] = False
    elif mutation == "ASSERT_CHANNEL_EXACT_ONE":
        s["force_term_census"]["coverage_checks"]["exactly_one_channel"] = False
        w["force_term_census"]["coverage_checks"]["exactly_one_channel"] = False
    elif mutation == "ASSERT_COUPLING_CENSUS_AGREE":
        s["coupling_source_census"]["entries"][0]["components"].pop()
    elif mutation == "ASSERT_COUPLING_INCIDENCE":
        s["coupling_source_census"]["coverage_checks"][
            "coupling_census_incidence_complete"
        ] = False
        w["coupling_source_census"]["coverage_checks"][
            "coupling_census_incidence_complete"
        ] = False
    elif mutation == "ASSERT_COUPLING_EXPECTED_REACHABLE":
        s["coupling_source_census"]["reachable_entries"].pop()
        w["coupling_source_census"]["reachable_entries"].pop()
        s["coupling_source_census"]["coverage_checks"][
            "expected_reachable_exact_set_equal"
        ] = False
        w["coupling_source_census"]["coverage_checks"][
            "expected_reachable_exact_set_equal"
        ] = False
    elif mutation == "ASSERT_DELTAO_ORDERED_MATRIX":
        s["coupling_source_census"]["ordered_deltaO_entries"].pop()
        w["coupling_source_census"]["ordered_deltaO_entries"].pop()
        s["coupling_source_census"]["coverage_checks"]["all_ordered_deltaO_present"] = False
        w["coupling_source_census"]["coverage_checks"]["all_ordered_deltaO_present"] = False
    elif mutation == "ASSERT_G8_INDEPENDENT_PROVENANCE":
        s["g8_ablation_inventory"]["generator_provenance"]["route"] = s[
            "coupling_source_census"
        ]["generator_provenance"]["route"]
        w["g8_ablation_inventory"]["generator_provenance"]["route"] = w[
            "coupling_source_census"
        ]["generator_provenance"]["route"]
    elif mutation == "ASSERT_G8_INVENTORY_AGREE":
        s["g8_ablation_inventory"]["entries"][0]["floor_tags"] = ["MUTATED"]
    elif mutation == "ASSERT_RADIATION_CENSUS_COVERAGE":
        for engine in (s, w):
            coupling = engine["coupling_source_census"]
            coupling["sources"] = [
                row for row in coupling["sources"]
                if not row["source_id"].startswith("radiation:")
            ]
            coupling["entries"] = [
                row for row in coupling["entries"]
                if not row["source_id"].startswith("radiation:")
            ]
            coupling["expected_entries"] = [
                value for value in coupling["expected_entries"]
                if ":radiation:" not in value
            ]
            coupling["reachable_entries"] = [
                value for value in coupling["reachable_entries"]
                if ":radiation:" not in value
            ]
            coupling["coverage_checks"]["source_ids"] = [
                value for value in coupling["coverage_checks"]["source_ids"]
                if not value.startswith("radiation:")
            ]
            inventory = engine["g8_ablation_inventory"]
            inventory["entries"] = [
                row for row in inventory["entries"]
                if not row["source_id"].startswith("radiation:")
            ]
            inventory["floor_coverage"] = {
                floor: [value for value in values if not value.startswith("radiation:")]
                for floor, values in inventory["floor_coverage"].items()
            }
    elif mutation == "ASSERT_G8_FLOOR":
        s["g8_ablation_inventory"]["floor_coverage"]["common_drain"] = []
        w["g8_ablation_inventory"]["floor_coverage"]["common_drain"] = []
        s["g8_ablation_inventory"]["coverage_checks"]["floor_subset_or_certified"] = False
        w["g8_ablation_inventory"]["coverage_checks"]["floor_subset_or_certified"] = False
    elif mutation == "ASSERT_G8_TO_COUPLING":
        fake = {
            "source_id": "action:invented",
            "mediators": ["h"],
            "floor_tags": [],
            "level2_disposition": "entry_witness",
            "entry_witness_slot": "tilt:indexed_h_tilt_profile",
        }
        s["g8_ablation_inventory"]["entries"].append(copy.deepcopy(fake))
        w["g8_ablation_inventory"]["entries"].append(copy.deepcopy(fake))
        s["g8_ablation_inventory"]["coverage_checks"]["every_G8_maps_to_coupling_source"] = False
        w["g8_ablation_inventory"]["coverage_checks"]["every_G8_maps_to_coupling_source"] = False
    elif mutation == "ASSERT_G8_REVERSE_LEVEL1":
        # Named non-floor omission tooth: brane_shear_gradient has only a channel tag.
        for engine in (s, w):
            engine["g8_ablation_inventory"]["entries"] = [
                row
                for row in engine["g8_ablation_inventory"]["entries"]
                if row["source_id"] != "action:brane_shear_gradient"
            ]
            engine["g8_ablation_inventory"]["coverage_checks"][
                "level1_disjoint_union_exact"
            ] = False
    elif mutation == "ASSERT_G8_LEVEL2":
        s["g8_ablation_inventory"]["entries"][0]["level2_disposition"] = "none"
        w["g8_ablation_inventory"]["entries"][0]["level2_disposition"] = "none"
        s["g8_ablation_inventory"]["coverage_checks"]["level2_exactly_one_disposition"] = False
        w["g8_ablation_inventory"]["coverage_checks"]["level2_exactly_one_disposition"] = False
    elif mutation == "ASSERT_G8_ENTRY_WITNESS_SLOTS":
        for engine in (s, w):
            by_key(engine["g8_ablation_inventory"]["entries"], "source_id")[
                "input:native_momentum"
            ]["entry_witness_slot"] = "open_leaf:native_momentum"
    elif mutation == "ASSERT_SLOT_TAXONOMY":
        s["availability_slots"] = [
            row for row in s["availability_slots"] if row["slot_id"] != "deltaO:h:u_T"
        ]
        w["availability_slots"] = [
            row for row in w["availability_slots"] if row["slot_id"] != "deltaO:h:u_T"
        ]
    elif mutation == "ASSERT_SLOT_ENGINE_AGREE":
        by_key(s["availability_slots"], "slot_id")["J:h"]["required_type"] = "MUTATED"
    elif mutation == "ASSERT_DERIVED_VALUE_DIGESTS":
        for engine in (s, w):
            by_key(engine["availability_slots"], "slot_id")[
                "domain:S_body_Omega_c"
            ]["value_digest"] = "0" * 64
    elif mutation == "ASSERT_SLOT_EXACT_ONE":
        for engine in (s, w):
            row = by_key(engine["availability_slots"], "slot_id")["J:h"]
            row["value_digest"] = "f" * 64
    elif mutation == "ASSERT_PRODUCER_CENSUS":
        for engine in (s, w):
            by_key(engine["availability_slots"], "slot_id")["J:h"]["producer_set"] = []
    elif mutation == "ASSERT_S1_TRACEABLE_CAUSE":
        for engine in (s, w):
            by_key(engine["availability_slots"], "slot_id")["J:h"]["producer_set"] = [
                "unfrozen_builder_claim:J:h"
            ]
            by_key(engine["availability_slots"], "slot_id")["J:h"]["witness"][
                "producer_set"
            ] = ["unfrozen_builder_claim:J:h"]
    elif mutation == "ASSERT_COMMITTED_CLOSURE_EXACT":
        for engine in (s, w):
            by_key(engine["availability_slots"], "slot_id")["J:h"]["witness"][
                "complete_committed_input_closure_digest"
            ] = "0" * 64
    elif mutation == "ASSERT_CANDIDATE_SCHEMA_PINNED":
        for engine in (s, w):
            by_key(engine["availability_slots"], "slot_id")["J:h"]["challenge"][
                "candidate_schema_pinned"
            ]["required_type"] = "MUTATED"
    elif mutation == "ASSERT_WITNESS_EXECUTABLE":
        for engine in (s, w):
            by_key(engine["availability_slots"], "slot_id")["J:h"]["witness"][
                "insufficiency_certificate"
            ]["status"] = "PROSE_ONLY"
    elif mutation == "ASSERT_RESTORE_MUTATION":
        for engine in (s, w):
            by_key(engine["availability_slots"], "slot_id")["J:h"]["witness"][
                "counterfactual_restore_mutation"
            ]["restored_insufficiency_status"] = "PASS_COMPUTED"
    elif mutation == "ASSERT_CHALLENGE_TERMINAL":
        for engine in (s, w):
            by_key(engine["availability_slots"], "slot_id")["J:h"]["challenge"][
                "attempted_candidates"
            ] = []
    elif mutation == "ASSERT_CHALLENGE_LOCKSTEP":
        for engine in (s, w):
            by_key(engine["availability_slots"], "slot_id")["J:h"]["challenge"][
                "kind"
            ] = "operator/domain well-posedness failure"
    elif mutation == "ASSERT_CHALLENGE_NOT_REFUTED":
        for engine in (s, w):
            by_key(engine["availability_slots"], "slot_id")["J:h"]["challenge"][
                "outcome"
            ] = "REFUTED"
    elif mutation == "ASSERT_CHALLENGE_COMPUTED":
        for engine in (s, w):
            by_key(engine["availability_slots"], "slot_id")["J:h"]["challenge"][
                "measurement_digest"
            ] = "0" * 64
    elif mutation == "ASSERT_CONTRACT_CLASS_SHARING":
        for engine in (s, w):
            rows = by_key(engine["availability_slots"], "slot_id")
            rows["J:h"]["derivability_contract_class"] = rows["J:u_T"][
                "derivability_contract_class"
            ]
    elif mutation == "ASSERT_DERIVED_DUAL_ENGINE":
        for engine in (s, w):
            by_key(engine["availability_slots"], "slot_id")["domain:S_body_Omega_c"].pop(
                "dual_engine_comparison_id", None
            )
    elif mutation == "ASSERT_OVERWRITE_DERIVED_TOOTH":
        # The harness invokes this separately through --anti-dodge; retained here for listing.
        return
    elif mutation == "ASSERT_FIRST_CONSTRUCTION_TOOTH":
        return
    elif mutation == "ASSERT_PROJECTION_FREEZE":
        s["projection_freeze"]["projection"] = "POST_HOC_FIT"
    elif mutation == "ASSERT_DIMENSIONAL_SCHEMA":
        s["stage0_dimensional_firewall"]["all_declared_inputs_have_dimension_triples"] = False
        w["stage0_dimensional_firewall"]["all_declared_inputs_have_dimension_triples"] = False
    elif mutation == "ASSERT_BANNED_TILT_CONVENTION":
        s["guard_evidence"]["banned_collective_coordinate_symbol_present"] = True
        w["guard_evidence"]["banned_collective_coordinate_symbol_present"] = True
    elif mutation == "ASSERT_FORBIDDEN_ANCESTRY":
        s["guard_evidence"]["forbidden_ancestry_nodes"] = ["imported_point_current"]
        w["guard_evidence"]["forbidden_ancestry_nodes"] = ["imported_point_current"]
    elif mutation == "ASSERT_ONE_BODY_SCOPE":
        s["guard_evidence"]["two_body_objects_constructed"] = ["force_sign"]
        w["guard_evidence"]["two_body_objects_constructed"] = ["force_sign"]
    elif mutation == "ASSERT_FIELD_DRIVEN_CLASSIFICATION":
        s["guard_evidence"]["classification_from_runtime_fields"] = False
        w["guard_evidence"]["classification_from_runtime_fields"] = False
    elif mutation == "ASSERT_TYPED_SINK_LIVENESS_SYMPY":
        s["sink_digest"] = ""
    elif mutation == "ASSERT_TYPED_SINK_LIVENESS_WOLFRAM":
        w["sink_digest"] = ""
    else:
        raise ValueError(f"unknown comparator mutation: {mutation}")


def compare(
    sympy: dict[str, Any],
    wolfram: dict[str, Any],
    repo: Path,
    mutation: str | None = None,
) -> list[dict[str, str]]:
    apply_mutation(sympy, wolfram, mutation)
    checks: list[dict[str, str]] = []

    def checked(condition: bool, assert_id: str, detail: str) -> None:
        require(condition, assert_id, detail)
        checks.append({"assert_id": assert_id, "status": "PASS"})

    checked(
        sympy.get("schema_version") == SCHEMA and wolfram.get("schema_version") == SCHEMA,
        "ASSERT_ENGINE_SCHEMA",
        "engine schema mismatch",
    )
    checked(
        sympy.get("engine") == "SymPy"
        and wolfram.get("engine") == "Wolfram"
        and sympy.get("independent_route") != wolfram.get("independent_route"),
        "ASSERT_ENGINE_INDEPENDENCE",
        "routes are not independently identified",
    )
    source_paths = {
        "phase_a_inputs": repo / "software/em_charge_attribute/u1_body_dynamics_inputs.yaml",
        "b1": repo
        / "software/em_charge_attribute/reports/u1_body_dynamics_artifacts/b1_final_results_snapshot.yaml",
        "b2_contract": repo
        / "software/em_charge_attribute/reports/u1_body_dynamics_artifacts/stage_b2_0_intake_radiative_contract/stage0_contract.yaml",
        "b2_sympy": repo
        / "software/em_charge_attribute/reports/u1_body_dynamics_artifacts/stage_b2_1_intake_radiative_production/sympy_b2.yaml",
        "b2_mathematica": repo
        / "software/em_charge_attribute/reports/u1_body_dynamics_artifacts/stage_b2_1_intake_radiative_production/mathematica_b2.yaml",
        "stage3": repo / "software/stage1_solver/reports/pathA_39_stage3_operator_parity_results.yaml",
    }
    live_digests = {key: sha256_path(path) for key, path in source_paths.items()}
    checked(
        sympy["source_digests"] == wolfram["source_digests"] == live_digests,
        "ASSERT_SOURCE_DIGESTS_LIVE",
        "engine source digests do not match live frozen inputs",
    )
    sf = sympy["frozen_assertions"]
    wf = wolfram["frozen_assertions"]
    checked(
        set(sf["b1_tilt_ingredient_types"])
        == set(wf["b1_tilt_ingredient_types"])
        == EXPECTED_B1_TILT,
        "ASSERT_B1_TILT_INDEX",
        "eight canonical B1 ingredient types are not exact",
    )
    checked(
        sf["b1_partition_state"]
        == wf["b1_partition_state"]
        == "partition_open_pending_B2",
        "ASSERT_B1_PARTITION_TERMINAL",
        "B1 terminal was rewritten or misread",
    )
    checked(
        {
            sf["b2_partition_state_sympy"],
            sf["b2_partition_state_mathematica"],
            wf["b2_partition_state_sympy"],
            wf["b2_partition_state_mathematica"],
        }
        == {"UNRESOLVED(return_closure)"},
        "ASSERT_B2_PARTITION_TERMINAL",
        "B2 terminal was rewritten or misread",
    )
    checked(
        sf["b2_stage0_disposition_sha256"]
        == wf["b2_stage0_disposition_sha256"]
        == "d94341173b1f1ac643bb05cf52dbac2300668ce6a50b8b5042ee2c7fa35cc54f",
        "ASSERT_B2_DISPOSITION_PIN",
        "sealed B2 disposition digest mismatch",
    )
    sympy_tilt_paths = {row["path"] for row in sf["b2_tilt_paths_sympy"]}
    math_tilt_paths = {row["path"] for row in sf["b2_tilt_paths_mathematica"]}
    checked(
        sympy_tilt_paths
        == math_tilt_paths
        == {row["path"] for row in wf["b2_tilt_paths_sympy"]}
        == {row["path"] for row in wf["b2_tilt_paths_mathematica"]},
        "ASSERT_B2_TILT_PATH_EQUALITY",
        "B2 explicit tilt-row path sets differ",
    )
    reconciliation_s = set(sf["reconciliation_expected_ids"])
    reconciliation_w = set(wf["reconciliation_expected_ids"])
    checked(
        reconciliation_s == reconciliation_w,
        "ASSERT_RECONCILIATION_EXACT",
        "engine reconciliation inventories differ",
    )
    reconciliation_counts = Counter(row.split("|", 1)[0] for row in reconciliation_s)
    checked(
        reconciliation_counts
        == Counter({"B1_LEAF": 8, "B2_TILT_PATH": 77, "B2_DEFERRED": 14040})
        and len(reconciliation_s) == 14125,
        "ASSERT_RECONCILIATION_CARDINALITY",
        f"unexpected reconciliation categories {dict(reconciliation_counts)}",
    )
    typed_s = by_key(sympy["typed_root_inventory"], "id")
    typed_w = by_key(wolfram["typed_root_inventory"], "id")
    checked(
        strip_engine_local(typed_s) == strip_engine_local(typed_w)
        and len(typed_s) == 48
        and sum(row["root_type"] == "native_action_term" for row in typed_s.values()) == 15
        and sum(row["root_type"] == "radiative_channel" for row in typed_s.values()) == 8,
        "ASSERT_ACTION_ROOT_SET",
        "typed root inventories differ or are not complete against frozen native roots",
    )
    deriv_s = normalized_derivation(sympy["native_derivation"])
    deriv_w = normalized_derivation(wolfram["native_derivation"])
    checked(
        deriv_s["second_variation"] == deriv_w["second_variation"],
        "ASSERT_SECOND_VARIATION_AGREE",
        "independent second variations disagree",
    )
    checked(
        bool(deriv_s["second_variation"]["chi_u_mixed_block_present"])
        and bool(deriv_w["second_variation"]["chi_u_mixed_block_present"]),
        "ASSERT_CHI_U_MIXED_BLOCK",
        "wall chi-u mixed block missing",
    )
    checked(
        deriv_s["S_body_Omega_c"] == deriv_w["S_body_Omega_c"],
        "ASSERT_S_BODY_AGREE",
        "S_body termwise constructions disagree",
    )
    enforce_unresolved_guard(sympy)
    enforce_unresolved_guard(wolfram)
    hessian_challenge_s = normalized_hessian_challenge(sympy)
    hessian_challenge_w = normalized_hessian_challenge(wolfram)
    checked(
        hessian_challenge_s["executed"] is True
        and hessian_challenge_w["executed"] is True
        and hessian_challenge_s["defining_predicate_evaluated"] is True
        and hessian_challenge_w["defining_predicate_evaluated"] is True,
        "ASSERT_HESSIAN_CHALLENGE_EXECUTED",
        "fresh raw-action Hessian challenge did not execute",
    )
    checked(
        hessian_challenge_s == hessian_challenge_w
        and hessian_challenge_s["candidate_value"]
        == deriv_s["second_variation"]
        and sympy["hessian_constructive_challenge"]["candidate_value_digest"]
        == digest(sympy["hessian_constructive_challenge"]["candidate_value"])
        and wolfram["hessian_constructive_challenge"]["candidate_value_digest"]
        == digest(wolfram["hessian_constructive_challenge"]["candidate_value"]),
        "ASSERT_HESSIAN_CHALLENGE_AGREE",
        "fresh constructive Hessian routes disagree or do not reproduce the derived value",
    )
    force_s = normalized_census(sympy["force_term_census"], "force")
    force_w = normalized_census(wolfram["force_term_census"], "force")
    checked(force_s == force_w, "ASSERT_FORCE_CENSUS_AGREE", "force censuses disagree")
    checked(
        force_s["coverage_checks"]["force_census_incidence_complete"]
        and all(row["incidence_complete"] for row in force_s["root_incidence"]),
        "ASSERT_FORCE_INCIDENCE",
        "force root incidence incomplete",
    )
    checked(
        force_s["expected_terms"] == force_s["reachable_residual_terms"]
        and force_s["coverage_checks"]["expected_reachable_exact_set_equal"],
        "ASSERT_FORCE_EXPECTED_REACHABLE",
        "force expected/reachable sets differ",
    )
    checked(
        force_s["coverage_checks"]["exactly_one_channel"],
        "ASSERT_CHANNEL_EXACT_ONE",
        "force channel ownership is not exactly one",
    )
    coupling_s = normalized_census(sympy["coupling_source_census"], "coupling")
    coupling_w = normalized_census(wolfram["coupling_source_census"], "coupling")
    checked(
        coupling_s == coupling_w,
        "ASSERT_COUPLING_CENSUS_AGREE",
        "coupling censuses disagree",
    )
    checked(
        coupling_s["coverage_checks"]["coupling_census_incidence_complete"],
        "ASSERT_COUPLING_INCIDENCE",
        "coupling source incidence incomplete",
    )
    checked(
        coupling_s["expected_entries"] == coupling_s["reachable_entries"]
        and coupling_s["coverage_checks"]["expected_reachable_exact_set_equal"],
        "ASSERT_COUPLING_EXPECTED_REACHABLE",
        "coupling expected/reachable sets differ",
    )
    checked(
        len(coupling_s["ordered_deltaO_entries"]) == 16
        and coupling_s["coverage_checks"]["all_ordered_deltaO_present"],
        "ASSERT_DELTAO_ORDERED_MATRIX",
        "full ordered 4x4 deltaO matrix missing",
    )
    g8_s = normalized_census(sympy["g8_ablation_inventory"], "g8")
    g8_w = normalized_census(wolfram["g8_ablation_inventory"], "g8")
    checked(
        g8_s == g8_w,
        "ASSERT_G8_INVENTORY_AGREE",
        "independent engine G8 inventories disagree",
    )
    radiation_roots = {
        root_id for root_id in typed_s if root_id.startswith("radiation:")
    }
    checked(
        radiation_roots == EXPECTED_RADIATION_ROOTS
        and {
            row["source_id"] for row in coupling_s["sources"]
            if row["source_id"].startswith("radiation:")
        }
        == EXPECTED_RADIATION_ROOTS
        and {
            row["source_id"] for row in g8_s["entries"]
            if row["source_id"].startswith("radiation:")
        }
        == EXPECTED_RADIATION_ROOTS,
        "ASSERT_RADIATION_CENSUS_COVERAGE",
        "radiation typed roots are not covered by both independent coupling/G8 walks",
    )
    routes = {
        force_s["generator_provenance"]["route"],
        coupling_s["generator_provenance"]["route"],
        g8_s["generator_provenance"]["route"],
    }
    checked(
        len(routes) == 3
        and "coupling_source_census"
        in g8_s["generator_provenance"]["not_derived_from"],
        "ASSERT_G8_INDEPENDENT_PROVENANCE",
        "force/coupling/G8 generators are not provenance-distinct",
    )
    checked(
        set(g8_s["floor_coverage"]) == EXPECTED_FLOOR
        and all(g8_s["floor_coverage"][key] for key in EXPECTED_FLOOR)
        and g8_s["coverage_checks"]["floor_subset_or_certified"],
        "ASSERT_G8_FLOOR",
        "named G8 floor lacks coverage",
    )
    coupling_source_ids = set(coupling_s["coverage_checks"]["source_ids"])
    g8_source_ids = {row["source_id"] for row in g8_s["entries"]}
    checked(
        g8_source_ids.issubset(coupling_source_ids)
        and g8_s["coverage_checks"]["every_G8_maps_to_coupling_source"],
        "ASSERT_G8_TO_COUPLING",
        "G8 entry has no coupling-census source",
    )
    checked(
        coupling_source_ids
        == g8_source_ids
        | set(g8_s["certified_nonentries"])
        | set(g8_s["witnessed_nonentries"])
        and g8_s["coverage_checks"]["level1_disjoint_union_exact"],
        "ASSERT_G8_REVERSE_LEVEL1",
        "source-level reverse G8 coverage failed",
    )
    checked(
        all(
            row["level2_disposition"]
            in {"executed_ablation", "entry_certificate", "entry_witness"}
            for row in g8_s["entries"]
        )
        and g8_s["coverage_checks"]["level2_exactly_one_disposition"],
        "ASSERT_G8_LEVEL2",
        "G8 per-entry disposition invalid",
    )
    slots_s = normalized_slots(sympy["availability_slots"])
    slots_w = normalized_slots(wolfram["availability_slots"])
    raw_slots_s = by_key(sympy["availability_slots"], "slot_id")
    raw_slots_w = by_key(wolfram["availability_slots"], "slot_id")
    categories = Counter(row["category"] for row in slots_s.values())
    checked(
        set(slots_s) == set(slots_w)
        and len(slots_s) == 58
        and sum(slot_id.startswith("tilt:") for slot_id in slots_s) == 8
        and sum(slot_id.startswith("deltaO:") for slot_id in slots_s) == 16
        and sum(row["category"].startswith("7.5a") for row in slots_s.values()) == 3
        and sum(row["category"] == "declared_OPEN_leaf" for row in slots_s.values()) == 12,
        "ASSERT_SLOT_TAXONOMY",
        f"slot taxonomy malformed: {dict(categories)}",
    )
    checked(
        all(
            row["level2_disposition"] != "entry_witness"
            or (
                row.get("entry_witness_slot") in raw_slots_s
                and raw_slots_s[row["entry_witness_slot"]]["disposition"] == "UNRESOLVED"
            )
            for row in g8_s["entries"]
        )
        and all(
            row["level2_disposition"] != "entry_witness"
            or (
                row.get("entry_witness_slot") in raw_slots_w
                and raw_slots_w[row["entry_witness_slot"]]["disposition"] == "UNRESOLVED"
            )
            for row in g8_w["entries"]
        ),
        "ASSERT_G8_ENTRY_WITNESS_SLOTS",
        "G8 entry witness does not resolve to a ratified UNRESOLVED availability slot",
    )
    checked(
        engine_route_certificates_valid(sympy, wolfram, repo),
        "ASSERT_ENGINE_LOCAL_ROUTE_CERTIFICATES",
        "semantic route IDs or engine-local function certificates are invalid",
    )
    checked(slots_s == slots_w, "ASSERT_SLOT_ENGINE_AGREE", "slot records disagree")
    checked(
        all(
            row["value_digest"] == digest(row["value"])
            for row in [*raw_slots_s.values(), *raw_slots_w.values()]
            if row["disposition"] == "DERIVED"
        )
        and all(
            engine["native_derivation"]["S_body_Omega_c"]["value_digest"]
            == digest(
                engine["native_derivation"]["S_body_Omega_c"]["canonical_terms"]
            )
            and engine["native_derivation"]["second_variation"]["value_digest"]
            == digest(
                engine["native_derivation"]["second_variation"]["term_records"]
            )
            for engine in (sympy, wolfram)
        ),
        "ASSERT_DERIVED_VALUE_DIGESTS",
        "derived slot or native value digest is not bound to its canonical value",
    )
    checked(
        all(
            (
                row["disposition"] == "DERIVED"
                and "value_digest" in row
                and "witness_id" not in row
                and "challenge_id" not in row
            )
            or (
                row["disposition"] == "UNRESOLVED"
                and "value_digest" not in row
                and "value" not in row
                and "witness_id" in row
                and "challenge_id" in row
            )
            for row in [*raw_slots_s.values(), *raw_slots_w.values()]
        ),
        "ASSERT_SLOT_EXACT_ONE",
        "slot exact-one disposition failed",
    )
    unresolved = [row for row in slots_s.values() if row["disposition"] == "UNRESOLVED"]
    derived = [row for row in slots_s.values() if row["disposition"] == "DERIVED"]
    checked(
        all(
            row["producer_set"]
            and row["producer_set"] == row["witness"]["producer_set"]
            for row in unresolved
        ),
        "ASSERT_PRODUCER_CENSUS",
        "unresolved slot has empty or inconsistent producer census",
    )
    authoritative_prefixes = (
        "B1:",
        "B2:",
        "Stage3:",
        "action:",
        "action_incident:",
        "ambient:",
        "closure:",
        "input:",
        "parameter_register:",
        "second_variation_incidence:",
    )
    checked(
        all(
            row["witness"]["datum_id"] == row["slot_id"]
            and (
                (
                    row["category"] == "declared_OPEN_leaf"
                    and len(row["producer_set"]) == 2
                    and any(
                        producer.startswith("parameter_register:")
                        for producer in row["producer_set"]
                    )
                )
                or all(
                    producer == "action:*"
                    or producer.startswith(authoritative_prefixes)
                    for producer in row["producer_set"]
                )
            )
            for row in unresolved
        ),
        "ASSERT_S1_TRACEABLE_CAUSE",
        "unavailability cause does not trace to the frozen authority/producer census",
    )
    committed_closure_digest_s = string_record_digest(
        [f"{name}={value}" for name, value in sympy["source_digests"].items()]
    )
    committed_closure_digest_w = string_record_digest(
        [f"{name}={value}" for name, value in wolfram["source_digests"].items()]
    )
    checked(
        all(
            row["witness"]["complete_committed_input_closure_digest"]
            == committed_closure_digest_s
            for row in sympy["availability_slots"]
            if row["disposition"] == "UNRESOLVED"
        )
        and all(
            row["witness"]["complete_committed_input_closure_digest"]
            == committed_closure_digest_w
            for row in wolfram["availability_slots"]
            if row["disposition"] == "UNRESOLVED"
        ),
        "ASSERT_COMMITTED_CLOSURE_EXACT",
        "witness committed-input closure is not exact-set digest bound",
    )
    checked(
        all(
            row["challenge"]["candidate_schema_pinned"]
            == {
                "required_type": row["required_type"],
                "required_dimensions": row["required_dimensions"],
                "domain": row["domain"],
            }
            for row in unresolved
        ),
        "ASSERT_CANDIDATE_SCHEMA_PINNED",
        "challenge candidate type/domain schema is not pinned to the datum slot",
    )
    checked(
        all(
            row["witness"]["kind"] in ALLOWED_KINDS
            and row["witness"]["insufficiency_certificate"]["status"]
            in {"PASS_COMPUTED", "FAIL_COMPUTED"}
            and row["witness"]["insufficiency_certificate"]["executed"] is True
            and row["witness"]["datum_id"] == row["slot_id"]
            for row in unresolved
        ),
        "ASSERT_WITNESS_EXECUTABLE",
        "typed unavailability witness malformed",
    )
    for slot_id in sorted(slots_s, key=str.casefold):
        row_s = raw_slots_s[slot_id]
        if row_s["disposition"] != "UNRESOLVED":
            continue
        row_w = raw_slots_w[slot_id]
        assert_id = f"ASSERT_WITNESS_INSUFFICIENCY:{slot_id}"
        baseline_s, valid_s = expected_witness_status(
            row_s["witness"]["insufficiency_certificate"], row_s["witness"]["kind"]
        )
        baseline_w, valid_w = expected_witness_status(
            row_w["witness"]["insufficiency_certificate"], row_w["witness"]["kind"]
        )
        restored_s = row_s["witness"]["counterfactual_restore_mutation"][
            "restored_certificate"
        ]
        restored_w = row_w["witness"]["counterfactual_restore_mutation"][
            "restored_certificate"
        ]
        restored_status_s, restored_valid_s = expected_witness_status(
            restored_s, row_s["witness"]["kind"]
        )
        restored_status_w, restored_valid_w = expected_witness_status(
            restored_w, row_w["witness"]["kind"]
        )
        checked(
            valid_s
            and valid_w
            and restored_valid_s
            and restored_valid_w
            and baseline_s == baseline_w == "PASS_COMPUTED"
            and restored_status_s == restored_status_w == "FAIL_COMPUTED",
            assert_id,
            "executed insufficiency certificate did not fail under its typed restore fixture",
        )
    checked(
        all(
            row["witness"]["counterfactual_restore_mutation"]["measured_by_engine"]
            and row["witness"]["counterfactual_restore_mutation"][
                "candidate_schema_comparison"
            ]["equal"]
            and row["witness"]["counterfactual_restore_mutation"][
                "baseline_insufficiency_status"
            ]
            == "PASS_COMPUTED"
            and row["witness"]["counterfactual_restore_mutation"][
                "restored_insufficiency_status"
            ]
            == "FAIL_COMPUTED"
            and row["witness"]["counterfactual_restore_mutation"]["assert_id"]
            == f"ASSERT_WITNESS_INSUFFICIENCY:{row['slot_id']}"
            for row in unresolved
        ),
        "ASSERT_RESTORE_MUTATION",
        "witness restore mutation does not kill insufficiency assert",
    )
    checked(
        all(
            row["challenge"]["outcome"] in {"REFUTED", "CONSTRUCTIVE_FAIL", "ERROR"}
            and row["challenge"]["attempted_candidates"]
            and not row["challenge"]["engine_certificate"]["empty_output"]
            and not row["challenge"]["engine_certificate"]["ill_typed_by_fiat"]
            for row in unresolved
        ),
        "ASSERT_CHALLENGE_TERMINAL",
        "challenge empty, ill-typed-by-fiat, or non-exhaustive",
    )
    checked(
        all(row["challenge"]["kind"] == row["witness"]["kind"] for row in unresolved),
        "ASSERT_CHALLENGE_LOCKSTEP",
        "challenge kind differs from witness kind",
    )
    checked(
        all(row["challenge"]["outcome"] == "CONSTRUCTIVE_FAIL" for row in unresolved),
        "ASSERT_CHALLENGE_NOT_REFUTED",
        "REFUTED/ERROR challenge cannot be banked",
    )
    checked(
        all(
            expected_challenge_result(raw_slots_s[row["slot_id"]])[2]
            and expected_challenge_result(raw_slots_w[row["slot_id"]])[2]
            for row in unresolved
        ),
        "ASSERT_CHALLENGE_COMPUTED",
        "challenge certificate is not bound to its executed raw scan and predicate evaluation",
    )
    classes = [row["derivability_contract_class"] for row in unresolved]
    checked(
        len(classes) == len(set(classes)),
        "ASSERT_CONTRACT_CLASS_SHARING",
        "heterogeneous slots share a derivability class without equality proof",
    )
    checked(
        len(derived) == 2
        and all(row.get("dual_engine_comparison_id") for row in derived),
        "ASSERT_DERIVED_DUAL_ENGINE",
        "derived slot lacks comparator binding",
    )
    projection_s = strip_engine_local(sympy["projection_freeze"])
    projection_w = strip_engine_local(wolfram["projection_freeze"])
    checked(
        projection_s == projection_w
        and projection_s["range"] == "ordered_pair(c_sv,Delta_i)"
        and len(projection_s["predicates"]) == 7,
        "ASSERT_PROJECTION_FREEZE",
        "fixed projection artifact differs or is incomplete",
    )
    checked(
        sympy["stage0_dimensional_firewall"][
            "all_declared_inputs_have_dimension_triples"
        ]
        and wolfram["stage0_dimensional_firewall"][
            "all_declared_inputs_have_dimension_triples"
        ]
        and projection_s["dimension_firewall"]["subtraction_homogeneous"],
        "ASSERT_DIMENSIONAL_SCHEMA",
        "stage-0 dimensional schema/projection firewall failed",
    )
    checked(
        not sympy["guard_evidence"]["banned_collective_coordinate_symbol_present"]
        and not wolfram["guard_evidence"]["banned_collective_coordinate_symbol_present"],
        "ASSERT_BANNED_TILT_CONVENTION",
        "banned signed tilt coordinate entered emitted object graph",
    )
    checked(
        not sympy["guard_evidence"]["forbidden_ancestry_nodes"]
        and not wolfram["guard_evidence"]["forbidden_ancestry_nodes"]
        and not sympy["guard_evidence"]["maxwell_forms_used_as_ancestry"]
        and not wolfram["guard_evidence"]["maxwell_forms_used_as_ancestry"],
        "ASSERT_FORBIDDEN_ANCESTRY",
        "forbidden imported ancestry present",
    )
    checked(
        not sympy["guard_evidence"]["two_body_objects_constructed"]
        and not wolfram["guard_evidence"]["two_body_objects_constructed"],
        "ASSERT_ONE_BODY_SCOPE",
        "stage 0 constructed a two-body object",
    )
    checked(
        sympy["guard_evidence"]["classification_from_runtime_fields"]
        and wolfram["guard_evidence"]["classification_from_runtime_fields"],
        "ASSERT_FIELD_DRIVEN_CLASSIFICATION",
        "classification is not field-driven",
    )
    sink_records_s = semantic_sink_records(sympy)
    sink_records_w = semantic_sink_records(wolfram)
    checked(
        bool(sympy["sink_digest"])
        and sympy["sink_measurement_record_count"] == len(sink_records_s)
        and sympy["sink_digest"] == string_record_digest(sink_records_s),
        "ASSERT_TYPED_SINK_LIVENESS_SYMPY",
        "SymPy typed sink is dead or is not observation-bound",
    )
    checked(
        bool(wolfram["sink_digest"])
        and wolfram["sink_measurement_record_count"] == len(sink_records_w)
        and wolfram["sink_digest"] == string_record_digest(sink_records_w),
        "ASSERT_TYPED_SINK_LIVENESS_WOLFRAM",
        "Wolfram typed sink is dead or is not observation-bound",
    )
    return checks


def run_contract_class_canary(
    sympy: dict[str, Any], wolfram: dict[str, Any], slot_id: str, repo: Path
) -> None:
    """Plant executed restore measurements into real artifacts, then compare."""
    plant_refuting_candidate(sympy, slot_id)
    plant_refuting_candidate(wolfram, slot_id)
    compare(sympy, wolfram, repo)
    raise AssertionDeath("MUTATION_NOOP", f"canary survived:{slot_id}")


def anti_dodge(
    sympy: dict[str, Any], wolfram: dict[str, Any], mode: str, repo: Path
) -> None:
    """Plant the B2-pattern Hessian disposition defect, then really compare."""
    first_time = mode == "first_construction"
    plant_hessian_disposition_defect(sympy, first_time)
    plant_hessian_disposition_defect(wolfram, first_time)
    compare(sympy, wolfram, repo)
    raise AssertionDeath("MUTATION_NOOP", f"Hessian anti-dodge survived:{mode}")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--sympy", required=True)
    parser.add_argument("--wolfram", required=True)
    parser.add_argument("--repo", required=True)
    parser.add_argument("--output")
    parser.add_argument("--mutation")
    parser.add_argument("--anti-dodge", choices=("overwrite", "first_construction"))
    parser.add_argument(
        "--anti-dodge-control", choices=("overwrite", "first_construction")
    )
    parser.add_argument("--canary-slot")
    parser.add_argument("--canary-control")
    parser.add_argument("--control-assert")
    parser.add_argument("--list-asserts", action="store_true")
    args = parser.parse_args()
    sympy_path = Path(args.sympy).resolve()
    wolfram_path = Path(args.wolfram).resolve()
    repo = Path(args.repo).resolve()
    sympy = load_yaml(sympy_path)
    wolfram = load_yaml(wolfram_path)
    try:
        if args.anti_dodge:
            anti_dodge(sympy, wolfram, args.anti_dodge, repo)
        if args.canary_slot:
            run_contract_class_canary(sympy, wolfram, args.canary_slot, repo)
        if args.anti_dodge_control:
            slot_id = "support:complete_action_second_variation"
            require(
                by_key(sympy["availability_slots"], "slot_id")[slot_id]["disposition"]
                == by_key(wolfram["availability_slots"], "slot_id")[slot_id]["disposition"]
                == "DERIVED",
                "MUTATION_NOOP",
                f"anti-dodge control unexpectedly contains the guarded defect:{args.anti_dodge_control}",
            )
        if args.canary_control:
            control_s = by_key(sympy["availability_slots"], "slot_id")[args.canary_control]
            control_w = by_key(wolfram["availability_slots"], "slot_id")[args.canary_control]
            _, passing_s, valid_s = expected_challenge_result(control_s)
            _, passing_w, valid_w = expected_challenge_result(control_w)
            require(
                valid_s and valid_w and not passing_s and not passing_w,
                "MUTATION_NOOP",
                f"canary control unexpectedly contains a passing candidate:{args.canary_control}",
            )
        checks = compare(sympy, wolfram, repo, args.mutation)
        if args.control_assert:
            require(
                args.control_assert in {row["assert_id"] for row in checks},
                "MUTATION_NOOP",
                f"defect-absent control did not reach:{args.control_assert}",
            )
    except AssertionDeath as death:
        print(f"ASSERTION_FAILED {death.assert_id}: {death.detail}", file=sys.stderr)
        return 1
    if args.mutation:
        print("MUTATION_NOOP", file=sys.stderr)
        return 3
    if args.anti_dodge_control or args.canary_control or args.control_assert:
        control = args.anti_dodge_control or args.canary_control or args.control_assert
        print(f"PHASEC_REAL_COMPARATOR_GUARD_ABSENT control={control}")
        return 0
    if args.list_asserts:
        print(yaml.safe_dump({"assert_ids": [row["assert_id"] for row in checks]}, sort_keys=False))
        return 0
    if not args.output:
        parser.error("--output is required for normal comparison")
    slots = sympy["availability_slots"]
    unresolved = [row for row in slots if row["disposition"] == "UNRESOLVED"]
    agreement = {
        "schema_version": "U1_PHASE_C_STAGE0_ENGINE_AGREEMENT_V1",
        "status": "ENGINE_AGREE",
        "engine_artifacts": {
            "sympy": {"path": str(sympy_path), "sha256": sha256_path(sympy_path)},
            "wolfram": {"path": str(wolfram_path), "sha256": sha256_path(wolfram_path)},
        },
        "comparison_count": len(checks),
        "comparisons": checks,
        "semantic_digests": {
            "S_body_Omega_c": digest(
                normalized_derivation(sympy["native_derivation"])["S_body_Omega_c"]
            ),
            "complete_action_second_variation": digest(
                normalized_derivation(sympy["native_derivation"])["second_variation"]
            ),
            "force_term_census": digest(
                normalized_census(sympy["force_term_census"], "force")
            ),
            "coupling_source_census": digest(
                normalized_census(sympy["coupling_source_census"], "coupling")
            ),
            "G8_ablation_inventory": digest(
                normalized_census(sympy["g8_ablation_inventory"], "g8")
            ),
            "availability_slots": digest(normalized_slots(slots)),
            "committed_input_closure": string_record_digest(
                [f"{name}={value}" for name, value in sympy["source_digests"].items()]
            ),
            "witness_measurements": digest(
                {
                    row["slot_id"]: strip_engine_local(
                        row["witness"]["insufficiency_certificate"]
                    )
                    for row in unresolved
                }
            ),
            "challenge_measurements": digest(
                {
                    row["slot_id"]: strip_engine_local(row["challenge"])
                    for row in unresolved
                }
            ),
            "hessian_constructive_challenge": digest(
                normalized_hessian_challenge(sympy)
            ),
            "projection_freeze": digest(strip_engine_local(sympy["projection_freeze"])),
            "reconciliation_inventory": digest(
                sorted(sympy["frozen_assertions"]["reconciliation_expected_ids"])
            ),
        },
        "availability_summary": sympy["availability_summary"],
        "reconciliation_summary": {
            "B1_LEAF": 8,
            "B2_TILT_PATH": 77,
            "B2_DEFERRED": 14040,
            "total": 14125,
        },
        "challenge_dual_engine_certificates": [
            {
                "slot_id": row["slot_id"],
                "challenge_id": row["challenge_id"],
                "semantic_route_id": row["challenge"]["engine_certificate"][
                    "semantic_route_id"
                ],
                "sympy_local_function": row["challenge"]["engine_certificate"][
                    "engine_local_function"
                ],
                "wolfram_local_function": by_key(
                    wolfram["availability_slots"], "slot_id"
                )[row["slot_id"]]["challenge"]["engine_certificate"][
                    "engine_local_function"
                ],
                "semantic_outcome": "CONSTRUCTIVE_FAIL",
                "kind": row["witness"]["kind"],
                "sympy_measurement_sha256": row["challenge"]["measurement_digest"],
                "wolfram_measurement_sha256": by_key(
                    wolfram["availability_slots"], "slot_id"
                )[row["slot_id"]]["challenge"]["measurement_digest"],
            }
            for row in unresolved
        ],
    }
    output = Path(args.output).resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(
        yaml.safe_dump(agreement, sort_keys=False, allow_unicode=True, width=120),
        encoding="utf-8",
    )
    print(
        f"PHASEC_STAGE0_ENGINE_AGREE comparisons={len(checks)} "
        f"derived={agreement['availability_summary']['DERIVED']} "
        f"unresolved={agreement['availability_summary']['UNRESOLVED']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
