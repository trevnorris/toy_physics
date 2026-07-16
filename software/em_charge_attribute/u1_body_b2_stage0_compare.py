#!/usr/bin/env python3
"""Independent semantic comparator for directive-v48 B2 stage 0."""

from __future__ import annotations

import argparse
import copy
import itertools
import re
from pathlib import Path
from typing import Any

import sympy as sp

from u1_body_b2_common import ROOT, digest, dump_yaml, load_yaml, require, sha256_file


ENGINE_SCHEMA = "U1_PHASE_B2_STAGE0_ENGINE_V6"
AGREEMENT_SCHEMA = "U1_PHASE_B2_STAGE0_MATH_AGREEMENT_V6"
ANCHOR = "5f73be9f0738030bf73165fda5432644de5f4074"
ENDPOINTS = [f"E{i}" for i in range(1, 6)]
WITNESS_KINDS = {
    "nonuniqueness/solvability failure",
    "authority_census_producer_absence",
    "dimensional_incompatibility",
    "operator_domain_well_posedness_failure",
}
RESTORE_KINDS = {"missing_input_leaf", "missing_selector", "domain/BC completion", "dimensioned_branch_type", "fixture operator/domain swap"}
RESERVED = {
    "mass_sector_tolerance": "ir.g9_tolerance.mass_rate_scale",
    "energy_sector_tolerance": "ir.g9_tolerance.energy_power_scale",
    "momentum_residual_norm": "ir.g9_tolerance.momentum_residual_norm",
    "return_energy_closure": "field::dimensioned_energy_return_branch",
    "truncation_validity_band": "ir.truncation.validity_band",
    "truncation_tolerance": "ir.truncation.relative_tolerance",
}


def expr(raw: Any) -> sp.Expr:
    text = str(raw).replace("^", "**")
    return sp.sympify(text, locals={"I": sp.I, "Pi": sp.pi, "Abs": sp.Abs, "delta": sp.Function("delta")})


def equal_expr(left: Any, right: Any, tooth: str, detail: str) -> None:
    require(sp.simplify(expr(left) - expr(right)) == 0, tooth, detail)


def indexed(rows: list[dict[str, Any]], key: str = "id") -> dict[str, dict[str, Any]]:
    result = {str(row[key]): row for row in rows}
    require(len(result) == len(rows), "B2_S0_C_SCHEMA", f"duplicate {key}")
    return result


def compare_nested_expressions(left: Any, right: Any, tooth: str, detail: str) -> None:
    if isinstance(left, list) and isinstance(right, list):
        require(len(left) == len(right), tooth, f"{detail}:length")
        for index, (a, b) in enumerate(zip(left, right)):
            compare_nested_expressions(a, b, tooth, f"{detail}[{index}]")
    elif isinstance(left, dict) and isinstance(right, dict):
        require(set(left) == set(right), tooth, f"{detail}:keys")
        for key in sorted(left):
            compare_nested_expressions(left[key], right[key], tooth, f"{detail}.{key}")
    else:
        equal_expr(left, right, tooth, detail)


def normalize_wolfram_jet(raw: str) -> str:
    """Map Wolfram's deliberately distinct affine-jet names to the public jet basis."""
    text = raw
    replacements = [
        ("curlUd", "curl_ud_"), ("curlU", "curl_u_"),
        ("gradChi", "grad_chi_"), ("gradN", "grad_n_"), ("gradH", "grad_h_"),
        ("thetaT", "theta_t"), ("uwT", "uw_t"), ("uT", "u_t_"), ("hT", "h_t"),
    ]
    for old, new in replacements:
        text = text.replace(old, new)
    text = re.sub(r"d([LR])(curl_ud_|curl_u_|grad_chi_|grad_n_|grad_h_|u_t_)([A-Za-z]+)", r"d\1_\2\3", text)
    text = re.sub(r"d([LR])(theta_t|uw_t|h_t|chi|uw|n)", r"d\1_\2", text)
    text = re.sub(r"d([LR])v([wxyz])", r"d\1_v_\2", text)
    text = re.sub(r"\b(v)([wxyz])0\b", r"\g<1>_\g<2>0", text)
    text = re.sub(r"\b(grad_(?:n|h|chi)_)([wxyz])0\b", r"\g<1>\g<2>0", text)
    text = re.sub(r"\b(curl_(?:u|ud)_)(xy|xz|xw|yz|yw|zw)0\b", r"\g<1>\g<2>0", text)
    return text


def normalize_wolfram_primitives(raw: Any) -> str:
    text = str(raw)
    for old, new in {
        "chiGrad2": "chi_grad2", "hGrad2": "h_grad2", "hT2": "h_t2",
        "nGrad2": "n_grad2", "thetaT": "theta_t", "uCurl2": "u_curl2",
        "uT2": "u_t2", "uwT2": "uw_t2",
    }.items():
        text = text.replace(old, new)
    return text


def expected_operator_ancestry(phase: dict[str, Any]) -> dict[str, list[str]]:
    token_fields = {
        "n": "delta_n", "theta_t": "theta", "v2": "theta", "n_grad2": "delta_n",
        "chiB": "delta_chiB", "chi_grad2": "delta_chiB", "ud_curl2": "u_d_transverse",
        "f_throat": "native_throat_fields", "f_mix": "native_mixed_fields", "u_t2": "u_T",
        "u_curl2": "u_T", "uw_t2": "u_w", "uw": "u_w", "h_t2": "h", "h_grad2": "h",
    }
    coefficient_names = set(phase["coefficients"])
    roots: dict[str, set[str]] = {}
    for row in phase["action_terms"]:
        tokens = set(re.findall(r"\b[A-Za-z][A-Za-z0-9_]*\b", row["expression"])) - coefficient_names
        roots[row["id"]] = {token_fields[token] for token in tokens if token in token_fields}
    parent = {field: field for field in sorted(set().union(*roots.values()))}

    def find(value: str) -> str:
        while parent[value] != value:
            parent[value] = parent[parent[value]]
            value = parent[value]
        return value

    def union(a: str, b: str) -> None:
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[rb] = ra

    for fields in roots.values():
        ordered = sorted(fields)
        for other in ordered[1:]:
            union(ordered[0], other)
    groups: dict[str, set[str]] = {}
    for field in parent:
        groups.setdefault(find(field), set()).add(field)
    if any(row["id"] == "geon_core_bundle" for row in phase["field_records"]):
        groups["geon_core_bundle"] = {"geon_core_bundle"}
    names = {
        ("delta_n", "theta"): "gnls_density_phase",
        ("delta_chiB", "u_d_transverse"): "wall_chi_u_coupled",
        ("native_throat_fields",): "throat_source_open", ("native_mixed_fields",): "wall_mix_open",
        ("u_T",): "brane_shear_transverse", ("u_w",): "brane_normal_local", ("h",): "h_scalar",
        ("geon_core_bundle",): "geon_open",
    }
    expected: dict[str, list[str]] = {}
    for fields in groups.values():
        oid = names[tuple(sorted(fields))]
        expected[oid] = ["geon_core_bundle"] if oid == "geon_open" else sorted(root for root, targets in roots.items() if targets & fields)
    return expected


def flatten_leaf_ids(prefix: str, value: Any) -> list[str]:
    if isinstance(value, dict):
        return [leaf for key in sorted(value) for leaf in flatten_leaf_ids(f"{prefix}.{key}", value[key])]
    if isinstance(value, list):
        return [leaf for index, child in enumerate(value) for leaf in flatten_leaf_ids(f"{prefix}[{index}]", child)]
    return [prefix]


def slot_names(schema: dict[str, Any], branches: list[str]) -> list[str]:
    axes = copy.deepcopy(schema["axes"])
    axes["operator_branch"] = branches
    slots = set(schema["stage0_fixed_datum_slots"])
    for spec in schema["stage0_expanded_datum_slots"]:
        for values in itertools.product(*(axes[axis] for axis in spec["axes"])):
            slots.add(spec["name"] + "|" + "|".join(f"{axis}={value}" for axis, value in zip(spec["axes"], values)))
    return sorted(slots)


def split_slot(slot: str) -> tuple[str, dict[str, str]]:
    pieces = slot.split("|")
    return pieces[0], {key: value for key, value in (piece.split("=", 1) for piece in pieces[1:])}


def expected_producer(slot: str, action_ids: list[str], operators: dict[str, dict[str, Any]]) -> tuple[str, list[str]]:
    base, axes = split_slot(slot)
    operator = axes.get("operator_branch")
    if base == "full_action_second_variation":
        return "v48:hessian<-second_variation(EVERY_committed_action_term)+incidence_residual", [f"action::{name}" for name in action_ids]
    if base == "operator_hessian_block":
        producers = [f"action::{name}" for name in operators[operator]["action_ancestry"]]
        if operator == "geon_open":
            producers = ["field::geon_core_bundle_operator_functional"]
        return "v48:operator_block<-full_action_second_variation_component", producers
    if base in {"native_mass_current_law", "native_energy_current_law"}:
        return "v48:sector_current<-native_action_balance_roots", [*(f"action::{name}" for name in action_ids), "field::native_continuity", "field::native_momentum"]
    if base.startswith("integrated_"):
        return "v48:integrated_balance<-native_current+Reynolds+typed_channel_functionals", [*(f"field::{name}" for name in ["native_continuity", "native_momentum", "E4_shear_lock", "E5_rayleigh", "return_closure"]), "field::outer_control_flux_native_functional"]
    if base in RESERVED:
        return "v48:upstream_freeze<-ONLY_7.0/inputs_ledger_typed_leaf", [RESERVED[base]]
    if base == "causal_retarded_definition":
        return "v48:causal_condition<-directive_definition_of_radiative_residue", ["directive::section_3.1"]
    if base == "native_noether_flux":
        producers = ["field::geon_core_bundle_action"] if operator == "geon_open" else [f"action::{name}" for name in operators[operator]["action_ancestry"]]
        return "v48:noether_flux<-committed_action_block", producers
    if base == "retarded_green_operator":
        return "v48:green<-native_operator+endpoint_domain+retarded_definition", [*(f"action::{name}" for name in operators[operator]["action_ancestry"] if name != "geon_core_bundle"), f"endpoint::{axes['endpoint']}", "directive::section_3.1"]
    if base.startswith("functional_"):
        suffix = {"functional_test_space": "test_space", "functional_topology": "topology", "functional_contact_subtraction_set": "contact_subtraction_set", "functional_operator_norm": "operator_norm"}[base]
        return "v48:functional_datum<-IR_scheme_typed_leaf", [f"ir.functional_analytic.{operator}.{axes['endpoint']}.{suffix}"]
    if base in {"trajectory_frequency_representation", "fourier_convention", "retained_local_kernel_definition", "zero_denominator_policy"}:
        return "v48:definition<-directive_stage0_definition", ["directive::section_3.1" if base == "trajectory_frequency_representation" else "directive::section_3.4"]
    if base == "truncation_frequency_variable":
        return "v48:what<-geometry_a+committed_c_ref_coefficients", ["directive::section_3.4", "ir.g9_tolerance.force_scale"]
    if base == "complete_comoving_restriction_map":
        return "v48:restriction<-kinematics+moduli_fixed_base_fields+normalizations", ["kinematics::coordinate_map", "kinematics::control_volume", "kinematics::base_state_class", "ambient::medium_rest_velocity", "mechanics::indexed_coordinates", "B1::moduli_fixed_base_fields"]
    if base == "ownership_current_improvement_relation":
        return "v48:ownership_quotient<-typed_current_improvement_superpotential", ["field::current_improvement_superpotential"]
    raise AssertionError(f"unknown v48 producer class:{slot}")


def relevant_closure(slot: str, producers: list[str], inventory: list[str]) -> list[str]:
    base, axes = split_slot(slot)
    if base in RESERVED:
        return sorted(item for item in inventory if item.startswith("ir.") or (base == "return_energy_closure" and item.startswith("field::return_closure")))
    if base.startswith("functional_"):
        return sorted(item for item in inventory if item.startswith("ir."))
    if base in {"integrated_mass_balance_identity", "integrated_momentum_balance_identity", "integrated_energy_balance_identity"}:
        return sorted(item for item in inventory if item.startswith(("action::", "field::", "endpoint::", "kinematics::")))
    if base in {"native_mass_current_law", "native_energy_current_law", "full_action_second_variation"}:
        return sorted(item for item in inventory if item.startswith(("action::", "field::")))
    if base in {"operator_hessian_block", "native_noether_flux", "retarded_green_operator"}:
        actual = set(producers) & set(inventory)
        if "endpoint" in axes:
            actual.add(f"endpoint::{axes['endpoint']}")
        if base == "retarded_green_operator":
            actual.add("directive::section_3.1")
        return sorted(actual)
    return sorted(set(producers) & set(inventory))


def expanded_schema(schema: dict[str, Any], branches: list[str], dispositions: dict[str, str]) -> list[str]:
    axes = copy.deepcopy(schema["axes"])
    axes["operator_branch"] = branches
    axes["unordered_operator_pair"] = [f"{a}__PAIR__{b}" for a, b in itertools.combinations(branches, 2)]
    result = set(schema["fixed_products"])
    for spec in schema["expanded_products"]:
        for values in itertools.product(*(axes[axis] for axis in spec["axes"])):
            result.add(spec["name"] + "|" + "|".join(f"{axis}={value}" for axis, value in zip(spec["axes"], values)))
    for slot, disposition in dispositions.items():
        result.add(f"stage0_datum|slot={slot}|disposition={disposition}")
    return sorted(result)


def independent_kernel_application(serialized: dict[str, Any]) -> dict[str, sp.Expr]:
    operator, kernel = serialized["operator_ast"], serialized["kernel_ast"]
    if operator["kind"] == "multiplication":
        residual = sp.simplify(expr(operator["coefficient"]) * expr(kernel["coefficient"]) - 1)
        return {"contact_residual_coefficient": residual, "total_distributional_residual": residual}
    bop, z, q2 = expr(operator["B"]), expr(operator["z"]), expr(kernel["q_squared"])
    normalization = expr(kernel["normalization"])
    rp = sp.Symbol("rp", positive=True)
    dimension = int(kernel["dimension"])
    bulk = sp.simplify(z - bop * q2)
    jump = sp.simplify(-normalization / rp ** (dimension - 1))
    contact = sp.simplify(-bop * jump - 1 / rp ** (dimension - 1))
    fa, fpa, va, vpa = sp.symbols("f_a fp_a v_a vp_a")
    omega, epsilon, gamma = sp.symbols("omega epsilon gammaSigma")
    coefficient = sp.sympify(str(kernel["left_mixing_coefficient"]).replace("^", "**"), locals={"I": sp.I, "f_a": fa, "fp_a": fpa, "v_a": va, "vp_a": vpa, "omega": omega, "epsilon": epsilon, "gammaSigma": gamma})
    if operator["boundary_kind"] == "Dirichlet":
        endpoint = sp.simplify(fa + coefficient * va)
    elif operator["boundary_kind"] == "Neumann":
        endpoint = sp.simplify(fpa + coefficient * vpa)
    else:
        h = -sp.I * (omega + sp.I * epsilon) * gamma
        endpoint = sp.simplify(bop * (fpa + coefficient * vpa) + h * (fa + coefficient * va))
    return {
        "bulk_left_coefficient": bulk, "bulk_right_coefficient": bulk, "kernel_dispersion_residual": bulk,
        "continuity_jump": sp.Integer(0), "derivative_jump": jump,
        "delta_coefficient": sp.simplify(-bop * jump), "contact_residual_coefficient": contact,
        "endpoint_term_coefficient": endpoint, "total_distributional_residual": sp.simplify(2 * bulk + contact + endpoint),
    }


def assert_typecheck(record: dict[str, Any], expected: bool, tooth: str, detail: str) -> None:
    check = record["candidate_typecheck"]
    computed = all(check[key] is True for key in ["candidate_present", "type_matches_required_type", "dimensions_match_required_dimensions", "identity_domain_membership"])
    require(check["computed_result"] is computed and record["candidate_is_well_typed"] is computed and computed is expected, tooth, detail)


def unresolved_guard(record: dict[str, Any]) -> None:
    if record["disposition_candidate"] == "UNRESOLVED" and record["candidate_is_well_typed"] is True and record["defining_predicate_result"] == "PASS":
        require(False, "stage0_unresolved_refuted", record["slot"])


def validate_datum_banks(
    sym: dict[str, Any], math: dict[str, Any], phase: dict[str, Any], mechanics: dict[str, Any],
    schema: dict[str, Any], operators: dict[str, dict[str, Any]],
) -> tuple[list[str], dict[str, dict[str, Any]], list[dict[str, Any]]]:
    branches = sorted(operators)
    expected_slots = slot_names(schema, branches)
    action_ids = sorted(row["id"] for row in phase["action_terms"])
    field_ids = sorted(row["id"] for row in phase["field_records"])
    closure_inventory = sorted([
        *(f"action::{name}" for name in action_ids), *(f"field::{name}" for name in field_ids),
        *(f"endpoint::{name}" for name in sorted(mechanics["endpoint_functionals"])),
        *flatten_leaf_ids("ir", phase["ir_scheme"]), "kinematics::coordinate_map", "kinematics::control_volume",
        "kinematics::base_state_class", "ambient::medium_rest_velocity", "mechanics::indexed_coordinates",
        "B1::moduli_fixed_collective_tangents", "B1::moduli_fixed_base_fields",
        "directive::section_1.4", "directive::section_3.1", "directive::section_3.4",
    ])
    banks = [sym["stage0_datum_bank"], math["stage0_datum_bank"]]
    for engine, bank in zip(["SymPy", "Mathematica"], banks):
        require(bank["directive_version"] == 48, "B2_S0_C_DISPOSITION_FLOOR", f"{engine}:directive")
        require(set(bank["expected_slots"]) == set(bank["reachable_slots"]) == set(expected_slots) and len(bank["expected_slots"]) == len(bank["reachable_slots"]) == bank["record_count"] == len(expected_slots), "B2_S0_C_DISPOSITION_FLOOR", f"{engine}:expected/reachable exact set")
        require(bank["expected_reachable_exact_set_equal"] is True, "B2_S0_C_DISPOSITION_FLOOR", f"{engine}:set assert")
        require(set(bank["committed_input_closure_inventory"]) == set(closure_inventory) and len(bank["committed_input_closure_inventory"]) == len(closure_inventory) and len(bank["committed_input_closure_inventory_sha256"]) == 64, "B2_S0_C_WITNESS_CLOSURE", f"{engine}:mechanical closure")
    sym_rows, math_rows = indexed(banks[0]["records"], "slot"), indexed(banks[1]["records"], "slot")
    require(set(sym_rows) == set(math_rows) == set(expected_slots), "B2_S0_C_DISPOSITION_FLOOR", "dual record set")
    producer_map = []
    final_dispositions = []
    for slot in expected_slots:
        sr, mr = sym_rows[slot], math_rows[slot]
        rule, producers = expected_producer(slot, action_ids, operators)
        closure = relevant_closure(slot, producers, closure_inventory)
        for engine, row in [("SymPy", sr), ("Mathematica", mr)]:
            unresolved_guard(row)
            require(row["producer_rule"] == rule and set(row["producer_ids"]) == set(producers) and len(row["producer_ids"]) == len(producers), "B2_S0_C_PRODUCER_MAP", f"{engine}:{slot}")
            require(set(row["committed_input_closure"]) == set(closure) and len(row["committed_input_closure"]) == len(closure) and row["closure_exact_set_assert"] == "PASS", "B2_S0_C_WITNESS_CLOSURE", f"{engine}:{slot}")
            require(row["disposition_candidate"] in {"DERIVED", "UNRESOLVED"}, "B2_S0_C_DISPOSITION_FLOOR", f"{engine}:{slot}:one-of")
        require(sr["disposition_candidate"] == mr["disposition_candidate"], "B2_S0_C_DISPOSITION_FLOOR", f"dual disposition:{slot}")
        require(sr["required_type"] == mr["required_type"] and sr["required_dimensions"] == mr["required_dimensions"], "B2_S0_C_DISPOSITION_FLOOR", f"dual type:{slot}")
        producer_map.append({"slot": slot, "rule": rule, "producers": producers})
        if sr["disposition_candidate"] == "DERIVED":
            for engine, row in [("SymPy", sr), ("Mathematica", mr)]:
                assert_typecheck(row, True, "B2_S0_C_DERIVED_TYPECHECK", f"{engine}:{slot}")
                challenge = row["derivability_challenge"]
                assert_typecheck(challenge, True, "B2_S0_C_DERIVED_TYPECHECK", f"{engine}:{slot}:challenge")
                require(challenge["terminal"] == "REFUTED(well-typed PASS candidate)" and challenge["defining_predicate_result"] == "PASS", "B2_S0_C_DERIVED_TYPECHECK", f"{engine}:{slot}:terminal")
            comparison_id = "dual_engine::" + digest({"slot": slot, "candidate_ref": sr["candidate_ref"], "producer_rule": rule})
            final_dispositions.append({"slot": slot, "disposition": "DERIVED", "value_digest": digest({"slot": slot, "comparison_id": comparison_id}), "dual_engine_comparison_id": comparison_id})
        else:
            for engine, row in [("SymPy", sr), ("Mathematica", mr)]:
                assert_typecheck(row, False, "B2_S0_C_CHALLENGE", f"{engine}:{slot}:record typecheck")
                witness, challenge = row["unavailability_witness"], row["derivability_challenge"]
                require(witness["datum_id"] == slot and witness["required_type"] == row["required_type"] and witness["required_dimensions"] == row["required_dimensions"], "B2_S0_C_WITNESS_SCHEMA", f"{engine}:{slot}:identity")
                require(set(witness["enumerated_committed_inputs"]) == set(closure) and len(witness["enumerated_committed_inputs"]) == len(closure) and witness["complete_closure_exact_set_equal"] is True and set(witness["authoritative_roots"]) == set(producers) and len(witness["authoritative_roots"]) == len(producers) and witness["directive_generated_producer_rule"] == rule, "B2_S0_C_WITNESS_CLOSURE", f"{engine}:{slot}:witness")
                census = witness["producer_census"]
                require({entry["producer"] for entry in census} == set(producers) and len(census) == len(producers) and all((not entry["in_closure"]) or entry["type_compatible"] is False for entry in census), "B2_S0_C_WITNESS_KIND", f"{engine}:{slot}:executable census")
                require(witness["kind"] in WITNESS_KINDS and witness["executable_certificate_result"] == "PASS", "B2_S0_C_WITNESS_KIND", f"{engine}:{slot}:kind")
                restore = witness["counterfactual_restore_mutation"]
                require(restore["ingredient_kind"] in RESTORE_KINDS and restore["producer_to_restore"] in producers and restore["restored_type_compatible"] is True, "B2_S0_C_WITNESS_RESTORE", f"{engine}:{slot}:restore schema")
                restored = [dict(entry, in_closure=True, type_compatible=True) if entry["producer"] == restore["producer_to_restore"] else entry for entry in census]
                restored_absence = all((not entry["in_closure"]) or entry["type_compatible"] is False for entry in restored)
                require(not restored_absence and restore["certificate_after_restore"] == "FAIL" and restore["failed_at_own_assert"] == "B2_S0_WITNESS_RESTORE", "B2_S0_C_WITNESS_RESTORE", f"{engine}:{slot}:executable restore")
                require(set(challenge["same_committed_input_closure"]) == set(closure) and len(challenge["same_committed_input_closure"]) == len(closure) and challenge["dag_separated_from_witness"] is True and challenge["shared_only_committed_inputs"] is True, "B2_S0_C_CHALLENGE", f"{engine}:{slot}:DAG/closure")
                require(challenge["candidate_schema_pinned"] == {"type": row["required_type"], "dimensions": row["required_dimensions"]}, "B2_S0_C_CHALLENGE", f"{engine}:{slot}:pin")
                assert_typecheck(challenge, False, "B2_S0_C_CHALLENGE", f"{engine}:{slot}:typecheck")
                require(challenge["constructive_attempt_nonempty"] is True and challenge["attempt_records"] and {attempt["producer"] for attempt in challenge["attempt_records"]} == set(producers) and len(challenge["attempt_records"]) == len(producers), "B2_S0_C_CHALLENGE", f"{engine}:{slot}:constructive attempt")
                require(challenge["kind"] == witness["kind"] and challenge["terminal"] == f"CONSTRUCTIVE_FAIL({witness['kind']})", "B2_S0_C_CHALLENGE", f"{engine}:{slot}:kind lockstep")
                if split_slot(slot)[0] in RESERVED:
                    require(row["unresolved_tag"] == "UNRESOLVED_BY_UPSTREAM_FREEZE" and rule.startswith("v48:upstream_freeze"), "B2_S0_C_UPSTREAM_FREEZE", f"{engine}:{slot}")
                    forbidden = {"mdot", "mdot*c_ref**2", "candidate_momentum_norm", "directive::section_1.4"}
                    require(not (forbidden & set(producers)) and "approved_frozen_7.0_ledger_entry" in row["required_type"], "B2_S0_C_UPSTREAM_FREEZE", f"{engine}:{slot}:non-operative candidate exclusion")
            require(sr["unavailability_witness"]["kind"] == mr["unavailability_witness"]["kind"], "B2_S0_C_CHALLENGE", f"dual kind:{slot}")
            certificate = digest({"slot": slot, "sympy_attempts": sr["derivability_challenge"]["attempt_records"], "wolfram_attempts": mr["derivability_challenge"]["attempt_records"], "kind": sr["unavailability_witness"]["kind"]})
            final_dispositions.append({"slot": slot, "disposition": "UNRESOLVED", "witness_id": sr["unavailability_witness"]["witness_id"], "challenge_id": sr["derivability_challenge"]["challenge_id"], "dual_engine_challenge_certificate_sha256": certificate})
    for engine, bank in zip(["SymPy", "Mathematica"], banks):
        emitted = indexed(bank["directive_generated_producer_map"], "slot")
        require(set(emitted) == set(expected_slots), "B2_S0_C_PRODUCER_MAP", f"{engine}:complete directive map")
        for expected_row in producer_map:
            actual = emitted[expected_row["slot"]]
            require(actual["rule"] == expected_row["rule"] and set(actual["producers"]) == set(expected_row["producers"]) and len(actual["producers"]) == len(expected_row["producers"]), "B2_S0_C_PRODUCER_MAP", f"{engine}:{expected_row['slot']}:banked map")
        require(len(bank["directive_generated_producer_map_sha256"]) == 64, "B2_S0_C_PRODUCER_MAP", f"{engine}:producer map digest")
    return expected_slots, sym_rows, final_dispositions


def compare(sym: dict[str, Any], math: dict[str, Any], config: dict[str, Any], input_path: Path) -> tuple[list[str], list[dict[str, Any]]]:
    require(config["startup_contract_commit"] == ANCHOR and config["directive_version"] == 48, "B2_S0_C_INPUT", "v48 anchor")
    require(sym["schema_version"] == math["schema_version"] == ENGINE_SCHEMA, "B2_S0_C_SCHEMA", "engine schema")
    require(sym["engine"] == "SymPy" and math["engine"] == "Mathematica", "B2_S0_C_SCHEMA", "engine identities")
    require(sym["independent_representation"] != math["independent_representation"], "B2_S0_C_INDEPENDENCE", "distinct representations")
    for route in ["source_parser_ir", "operator_ir", "balance_route_A_ir", "balance_route_B_ir", "resolvent_ir", "ownership_ir"]:
        require(sym["route_separation"][route] != math["route_separation"][route], "B2_S0_C_INDEPENDENCE", route)
    require(not sym["route_separation"]["shared_derivation_helpers"] and not math["route_separation"]["shared_derivation_helpers"], "B2_S0_C_INDEPENDENCE", "no shared reduced helper")
    require(sym["input_sha256"] == math["input_sha256"] == sha256_file(input_path), "B2_S0_C_INPUT", "input digest")
    for source in ["phase_a_inputs", "phase_b1_inputs", "directive", "v48_obligation_schema"]:
        require(sym["source_digests"][source] == math["source_digests"][source], "B2_S0_C_INPUT", source)
    require(sym["source_digests"]["b1_artifact"] == sha256_file(ROOT / config["contracts"]["sympy_b1"]), "B2_S0_C_INPUT", "SymPy B1 artifact")
    require(math["source_digests"]["b1_artifact"] == sha256_file(ROOT / config["contracts"]["mathematica_b1"]), "B2_S0_C_INPUT", "Mathematica B1 artifact")
    require(sym["causal_definition"] == math["causal_definition"] == config["causal_definition"], "B2_S0_C_CAUSAL", "causal definition")
    require(sym["fourier_convention"] == math["fourier_convention"] == config["fourier_convention"], "B2_S0_C_FOURIER", "Fourier convention")

    phase = load_yaml(ROOT / config["contracts"]["phase_a_inputs"])
    mechanics = load_yaml(ROOT / config["contracts"]["phase_b1_inputs"])
    schema_path = ROOT / config["contracts"]["v48_obligation_schema"]
    schema = load_yaml(schema_path)
    expected = expected_operator_ancestry(phase)
    so, mo = indexed(sym["operator_inventory"]), indexed(math["operator_inventory"])
    require(set(so) == set(mo) == set(expected) and len(so) == 8, "B2_S0_C_OPERATORS", "connected action Hessian sectors")

    s2, m2 = sym["complete_action_second_variation"], math["complete_action_second_variation"]
    require(s2["status"] == m2["status"] and s2["wall_gate_block_status"] == m2["wall_gate_block_status"] == "DERIVED_SYMBOLIC_FROZEN_BASE_JET_NO_COUPLING_SUPPRESSED", "B2_S0_C_FULL_HESSIAN", "full second variation status")
    s2rows, m2rows = indexed(s2["termwise_records"]), indexed(m2["termwise_records"])
    action_ids = {row["id"] for row in phase["action_terms"]}
    require(set(s2rows) == set(m2rows) == action_ids, "B2_S0_C_FULL_HESSIAN", "every action term")
    for root in sorted(action_ids):
        require(s2rows[root]["status"] == m2rows[root]["status"], "B2_S0_C_FULL_HESSIAN", f"{root}:status")
        if s2rows[root]["status"].startswith("DERIVED_"):
            equal_expr(s2rows[root]["bilinear_second_variation"], normalize_wolfram_jet(m2rows[root]["bilinear_second_variation"]), "B2_S0_C_FULL_HESSIAN", f"{root}:bilinear")
        else:
            require(root in {"throat_source", "wall_mix"} and s2rows[root]["bilinear_second_variation"] == m2rows[root]["bilinear_second_variation"], "B2_S0_C_FULL_HESSIAN", f"{root}:honest open functional")
    gate = expr(s2rows["wall_shear_gate"]["bilinear_second_variation"])
    require(sp.diff(gate, sp.Symbol("dL_chi"), sp.Symbol("dR_curl_ud_xy")) != 0 and sp.diff(gate, sp.Symbol("dR_chi"), sp.Symbol("dL_curl_ud_xy")) != 0, "B2_S0_C_FULL_HESSIAN", "live chi-u block")
    require(sp.diff(gate, sp.Symbol("dL_curl_ud_xy"), sp.Symbol("dR_curl_ud_xy")).has(sp.Symbol("chi0")), "B2_S0_C_FULL_HESSIAN", "live u-u frozen chi0 block")

    for oid in sorted(so):
        sr, mr = so[oid], mo[oid]
        require(sr["field_block"] == mr["field_block"] and set(sr["action_ancestry"]) == set(mr["action_ancestry"]) == set(expected[oid]), "B2_S0_C_OPERATORS", f"{oid}:incidence")
        for key in ["family", "support_dimension", "spatial_order", "time_order", "status", "radial_kind"]:
            require(sr[key] == mr[key], "B2_S0_C_OPERATORS", f"{oid}:{key}")
        if str(sr["operator_fourier"]).startswith("UNRESOLVED"):
            require(sr["operator_fourier"] == mr["operator_fourier"], "B2_S0_C_OPERATORS", f"{oid}:open")
        else:
            equal_expr(sr["operator_fourier"], mr["operator_fourier"], "B2_S0_C_OPERATORS", f"{oid}:operator")
        if oid in {"gnls_density_phase", "wall_chi_u_coupled"}:
            if oid == "gnls_density_phase":
                compare_nested_expressions(sr["matrix_operator_fourier"], mr["matrix_operator_fourier"], "B2_S0_C_OPERATORS", f"{oid}:matrix")
            else:
                for key in ["H_chi_chi", "H_chi_u_reduced_polarization", "H_u_chi_reduced_polarization", "H_u_u_reduced_polarization"]:
                    equal_expr(sr["matrix_operator_fourier"][key], mr["matrix_operator_fourier"][key], "B2_S0_C_FULL_HESSIAN", f"wall:{key}")
                require(sr["frozen_symbolic_base_coefficients"] == mr["frozen_symbolic_base_coefficients"] == ["chi0", "curl_ud0"], "B2_S0_C_FULL_HESSIAN", "no Berry-only gate assumption")

    si, mi = indexed(sym["operator_action_incidence"], "action_root"), indexed(math["operator_action_incidence"], "action_root")
    require(set(si) == set(mi) == action_ids, "B2_S0_C_INCIDENCE", "every action root")
    for root in sorted(action_ids):
        require(si[root]["operator_ids"] == mi[root]["operator_ids"] and si[root]["incidence_count"] == mi[root]["incidence_count"] > 0, "B2_S0_C_INCIDENCE", root)
    require(sym["operator_action_incidence_residual"] == math["operator_action_incidence_residual"] == [], "B2_S0_C_INCIDENCE", "zero incidence residual")

    sn, mn = indexed(sym["native_noether_derivations"]), indexed(math["native_noether_derivations"])
    require(set(sn) == set(mn) == set(so), "B2_S0_C_NOETHER", "sector coverage")
    for oid in sorted(sn):
        sr, mr = sn[oid], mn[oid]
        require(sr["status"] == mr["status"] and sr["action_ancestry"] == mr["action_ancestry"] and sr["spatial_dimension"] == mr["spatial_dimension"], "B2_S0_C_NOETHER", f"{oid}:metadata")
        if sr["status"].startswith("DERIVED_FROM"):
            require(sr["component_labels"] == mr["component_labels"] and sr["fields"] == mr["fields"], "B2_S0_C_NOETHER", f"{oid}:basis")
            equal_expr(sr["native_lagrangian_density"], mr["native_lagrangian_density"], "B2_S0_C_NOETHER", f"{oid}:L")
            compare_nested_expressions(sr["euler_operators"], mr["euler_operators"], "B2_S0_C_NOETHER", f"{oid}:EL")
            compare_nested_expressions(sr["canonical_tensor_components"]["T_mu__nu"], mr["canonical_tensor_components"]["T_mu__nu"], "B2_S0_C_NOETHER", f"{oid}:tensor")
            require(all(expr(value) == 0 for value in sr["computed_component_residuals"] + mr["computed_component_residuals"]), "B2_S0_C_NOETHER", f"{oid}:identity residual")
            require(set(sr["action_contributions"]) == set(mr["action_contributions"]) == set(sr["action_ancestry"]), "B2_S0_C_NOETHER", f"{oid}:termwise roots")
        else:
            equal_expr(sr.get("native_source_expression"), mr.get("native_source_expression"), "B2_S0_C_NOETHER", f"{oid}:open source")

    scur, mcur = sym["current_derivations"], math["current_derivations"]
    require(scur["ambient_spatial_dimension"] == mcur["ambient_spatial_dimension"] == 4 and not scur["banned_particle_form_used"] and not mcur["banned_particle_form_used"], "B2_S0_C_CURRENTS", "native current metadata")
    for sector in ["mass", "momentum", "energy"]:
        sr, mr = scur["sectors"][sector], mcur["sectors"][sector]
        compare_nested_expressions(sr["comoving_current_components"], mr["comoving_current_components"], "B2_S0_C_CURRENTS", f"{sector}:current")
        require(expr(sr["computed_pullback_identity_residual"]) == expr(mr["computed_pullback_identity_residual"]) == 0, "B2_S0_C_CURRENTS", f"{sector}:pullback")
        require(sr["missing_native_current_laws"] == mr["missing_native_current_laws"] and sr["open_current_contributions"] == mr["open_current_contributions"], "B2_S0_C_CURRENTS", f"{sector}:open propagation")
    # Cross-engine agreement alone would accept the same jointly corrupted
    # current.  Reconstruct the fully derived U(1) mass current here from its
    # native action-variation records, independently of both emitted pullbacks.
    for engine, current in [("sympy", scur["sectors"]["mass"]), ("wolfram", mcur["sectors"]["mass"])]:
        ancestry = current["action_derivation"]
        equal_expr(current["native_density"], ancestry["mass_density"], "B2_S0_C_CURRENTS", f"{engine}:mass density ancestry")
        for axis, lab_flux, comoving in zip(["x", "y", "z", "w"], current["native_lab_flux_components"], current["comoving_current_components"]):
            equal_expr(lab_flux, ancestry["mass_current_components"][axis], "B2_S0_C_CURRENTS", f"{engine}:mass {axis} lab ancestry")
            reconstructed = expr(lab_flux) - sp.Symbol(f"V_{axis}") * expr(current["native_density"])
            require(sp.simplify(expr(comoving) - reconstructed) == 0, "B2_S0_C_CURRENTS", f"{engine}:mass {axis} co-moving reconstruction")

    sb, mb = sym["integrated_balance_identities"], math["integrated_balance_identities"]
    require(sb["status"] == mb["status"] == "UNRESOLVED(complete_integrated_balance_family)", "B2_S0_C_BALANCES", "honest family exit")
    for sector in ["mass", "momentum", "energy"]:
        sr, mr = sb["sectors"][sector], mb["sectors"][sector]
        require(sr["status"] == mr["status"] == f"UNRESOLVED(complete_{sector}_integrated_balance)" and sr["surrogate_symbol_list_used"] is mr["surrogate_symbol_list_used"] is False, "B2_S0_C_BALANCES", f"{sector}:no surrogate")
        ar, am = sr["route_A_native_reynolds"], mr["route_A_native_reynolds"]
        compare_nested_expressions(ar["native_comoving_current_components"], am["native_comoving_current_components"], "B2_S0_C_BALANCE_ROUTE_A", f"{sector}:actual current")
        current_row = scur["sectors"][sector]
        expected_route_a_sources = set(current_row.get("action_derivation", {}).get("action_ancestry", [])) | {cause.split(":", 1)[-1] for cause in current_row["missing_native_current_laws"]}
        require(ar["reynolds_transport_result"] == am["reynolds_transport_result"] and ar["missing_data"] == am["missing_data"] and set(ar["source_omission_premises"]) == set(am["source_omission_premises"]) == expected_route_a_sources and expected_route_a_sources, "B2_S0_C_BALANCE_ROUTE_A", f"{sector}:Reynolds/source premises")
        require(expr(ar["computed_pullback_residual"]) == expr(am["computed_pullback_residual"]) == 0, "B2_S0_C_BALANCE_ROUTE_A", f"{sector}:pullback")
        br, bm = sr["route_B_authenticated_typed_roots"], mr["route_B_authenticated_typed_roots"]
        broots, bmroots = indexed(br["authenticated_typed_roots"]), indexed(bm["authenticated_typed_roots"])
        require(broots == bmroots and br["complete_signed_expression"] is bm["complete_signed_expression"] is None and br["missing_data"] == bm["missing_data"], "B2_S0_C_BALANCE_ROUTE_B", f"{sector}:typed-root constructive exit")
        require(set(br["source_omission_premises"]) == set(bm["source_omission_premises"]) == {"native_continuity", "native_momentum", "return_closure", "E4_shear_lock", "E5_rayleigh", "outer_control_flux"}, "B2_S0_C_BALANCE_ROUTE_B", f"{sector}:source premises")
        require(sr["complete_signed_expression_residual"] is None and sr["reconstruction_residual_coefficients"] is None, "B2_S0_C_BALANCES", f"{sector}:no zero-by-construction")
        require(sr["g9_stage0_router_output"] == mr["g9_stage0_router_output"] and sr["g9_stage0_router_output"]["verdict"].startswith("UNRESOLVED("), "B2_S0_C_G9_ROUTER", sector)
    require(sb["router"]["fixture_controls"] == mb["router"]["fixture_controls"] and sb["router"]["fixture_controls"]["production_consulted"] is False and sb["router"]["fixture_controls"]["surrogate_removal_residual_used"] is False, "B2_S0_C_G9_ROUTER", "no fake balance evaluator")

    sg, mg = indexed(sym["endpoint_resolvent_cells"], "cell"), indexed(math["endpoint_resolvent_cells"], "cell")
    require(set(sg) == set(mg) and len(sg) == 5 * len(so), "B2_S0_C_RESOLVENT", "per-branch endpoint coverage")
    for cell in sorted(sg):
        sr, mr = sg[cell], mg[cell]
        operator_id, endpoint = sr["operator"], sr["endpoint"]
        require(operator_id == mr["operator"] and endpoint == mr["endpoint"], "B2_S0_C_RESOLVENT", f"{cell}:identity")
        endpoint_record = mechanics["endpoint_functionals"][endpoint]
        variational_class = endpoint_record["variational_class"]
        if operator_id == "brane_shear_transverse":
            if variational_class == "nonholonomic_constraint":
                expected_domain = "OPEN_projected_constraint"
            elif variational_class == "Rayleigh":
                expected_domain = "Rayleigh_Robin"
            elif "tangent" in endpoint_record.get("essential", {}):
                expected_domain = "Dirichlet"
            else:
                expected_domain = "Neumann"
        elif operator_id == "gnls_density_phase":
            expected_domain = "UNRESOLVED"
        elif operator_id in {"wall_chi_static", "h_scalar"}:
            expected_domain = "Neumann"
        elif operator_id == "brane_normal_local":
            expected_domain = "algebraic_contact"
        else:
            expected_domain = "UNRESOLVED"
        require(sr["endpoint_domain_kind"] == mr["endpoint_domain_kind"] == expected_domain, "B2_S0_C_RESOLVENT", f"{cell}:authenticated domain")
        a, b = sr["retarded_green_operator"], mr["retarded_green_operator"]
        require(a["status"] == b["status"] and a["representation"] == b["representation"], "B2_S0_C_RESOLVENT", f"{cell}:status")
        if a["representation"] is not None:
            require(a["executable_kernel"]["schema"] == b["executable_kernel"]["schema"], "B2_S0_C_RESOLVENT", f"{cell}:kernel schema")
            for engine, green in [("SymPy", a), ("Mathematica", b)]:
                applied = independent_kernel_application(green["executable_kernel"])
                for key, value in applied.items():
                    equal_expr(value, green["operator_applied_to_serialized_kernel_minus_identity"][key], "B2_S0_C_RESOLVENT", f"{cell}:{engine}:{key}")
                require(sp.simplify(applied["total_distributional_residual"]) == 0, "B2_S0_C_RESOLVENT", f"{cell}:{engine}:operator-on-kernel")
        control_a, control_b = sr["known_nonzero_control"], mr["known_nonzero_control"]
        require(control_a["status"] == control_b["status"], "B2_S0_C_NATIVE_CONTROL", f"{cell}:status")
        if control_a["status"] == "DERIVED_NONZERO_NATIVE_TENSOR_CONTROL":
            require(control_a["tensor_component_path"] == control_b["tensor_component_path"] and not control_a["route_separation"]["shared_reduced_helpers"] and not control_b["route_separation"]["shared_reduced_helpers"], "B2_S0_C_NATIVE_CONTROL", f"{cell}:route separation")
            equal_expr(control_a["sector_tensor_component"], control_b["sector_tensor_component"], "B2_S0_C_NATIVE_CONTROL_TENSOR_ROUTE", f"{cell}:actual T")
            math_production = str(control_b["production_route"]["flux"]).replace("amplitudeAbsSquared", "amplitude_abs_squared")
            math_oracle = str(control_b["oracle_route"]["flux"]).replace("amplitudeAbsSquared", "amplitude_abs_squared")
            amplitude, omega_symbol, k_symbol = sp.symbols("amplitude_abs_squared omega k")
            expected_production = expr(control_a["production_route"]["extracted_native_stiffness"]) * amplitude * omega_symbol * k_symbol / 2
            expected_oracle = expr(control_a["oracle_route"]["group_velocity"]) * expr(control_a["oracle_route"]["on_shell_energy_density"])
            equal_expr(control_a["production_route"]["flux"], expected_production, "B2_S0_C_NATIVE_CONTROL_TENSOR_ROUTE", f"{cell}:production premise")
            equal_expr(control_a["oracle_route"]["flux"], expected_oracle, "B2_S0_C_NATIVE_CONTROL_ORACLE_ROUTE", f"{cell}:oracle premise")
            equal_expr(control_a["production_route"]["flux"], math_production, "B2_S0_C_NATIVE_CONTROL_TENSOR_ROUTE", f"{cell}:production")
            equal_expr(control_a["oracle_route"]["flux"], math_oracle, "B2_S0_C_NATIVE_CONTROL_ORACLE_ROUTE", f"{cell}:oracle")
            require(expr(control_a["computed_route_residual"]) == expr(control_b["computed_route_residual"]) == 0 and expr(control_a["production_route"]["flux"]) != 0, "B2_S0_C_NATIVE_CONTROL", f"{cell}:nonzero agreement")
    require(sg["brane_shear_transverse|E4"]["retarded_green_operator"]["status"] == "UNRESOLVED(E4_shear_lock_constraint_domain)", "B2_S0_C_RESOLVENT", "honest E4 projected domain")

    sf, mf = indexed(sym["functional_analytic_test_data"], "cell"), indexed(math["functional_analytic_test_data"], "cell")
    require(set(sf) == set(mf) == set(sg), "B2_S0_C_FUNCTIONAL", "functional cell coverage")
    for cell in sf:
        for row in [sf[cell], mf[cell]]:
            require(row["status"].startswith("UNRESOLVED(") and row["test_space"] is row["topology"] is row["contact_subtraction_set"] is row["norm"] is None and "evades" in row["endpoint_honesty"], "B2_S0_C_FUNCTIONAL", cell)

    ssens, msens = indexed(sym["green_sensitivity_matrix"], "operator"), indexed(math["green_sensitivity_matrix"], "operator")
    require(set(ssens) == set(msens) == set(so), "B2_S0_C_SENSITIVITY", "every operator branch including brane-normal and wall")
    for oid in ssens:
        require(ssens[oid]["status"] == msens[oid]["status"] and len(ssens[oid]["mutations"]) == len(msens[oid]["mutations"]), "B2_S0_C_SENSITIVITY", f"{oid}:status")
        for srow, mrow in zip(ssens[oid]["mutations"], msens[oid]["mutations"]):
            require(srow["parameter"] == mrow["parameter"] and srow["must_change"] == mrow["must_change"] and srow["proved_invariant"] == mrow["proved_invariant"], "B2_S0_C_SENSITIVITY", f"{oid}:mutation")
            recomputed_derivative = sp.simplify(sp.diff(expr(so[oid]["operator_fourier"]), sp.Symbol(srow["parameter"])))
            require(recomputed_derivative != 0, "B2_S0_C_SENSITIVITY", f"{oid}:{srow['parameter']}:live parameter")
            equal_expr(srow["operator_derivative"], recomputed_derivative, "B2_S0_C_SENSITIVITY", f"{oid}:{srow['parameter']}:recomputed")
            equal_expr(srow["operator_derivative"], mrow["operator_derivative"], "B2_S0_C_SENSITIVITY", f"{oid}:{srow['parameter']}")

    sr, mr = sym["restriction_map"], math["restriction_map"]
    require(sr["id"] == mr["id"] == "B1_ACTUAL_EXPRESSION_P0_RESTRICTION_V4" and sr["status"] == mr["status"] and sr["tautological_reconstruction_used"] is mr["tautological_reconstruction_used"] is False, "B2_S0_C_RESTRICTION", "actual expression restriction")
    require(sr["substitution"] == mr["substitution"] == {"p_x": 0, "p_y": 0, "p_z": 0} and sr["frozen_comoving_data"] == mr["frozen_comoving_data"], "B2_S0_C_RESTRICTION", "frozen substitutions and normalization")
    equal_expr(sr["actual_expression_executions"]["action"]["source_expression"], normalize_wolfram_primitives(mr["actual_expression_executions"]["action"]["source_expression"]), "B2_S0_C_RESTRICTION", "actual action")
    require(all(row["computed_direct_substitution_residual"] == "0" for row in sr["actual_expression_executions"]["flux"] + mr["actual_expression_executions"]["flux"] + sr["actual_expression_executions"]["reconciliation"] + mr["actual_expression_executions"]["reconciliation"]), "B2_S0_C_RESTRICTION", "executed action/flux/reconciliation")

    sox, mox = sym["ownership_convention"], math["ownership_convention"]
    require(sox["status"] == mox["status"] == "UNRESOLVED(named_current_improvement_datum)" and sox["concrete_current_improvement_relations"] == mox["concrete_current_improvement_relations"] == [], "B2_S0_C_OWNERSHIP", "honest missing improvement relation")
    require(sox["surrogate_IBP_relation_used_as_current_improvement"] is mox["surrogate_IBP_relation_used_as_current_improvement"] is False and sox["overlap_subtraction_count"] == mox["overlap_subtraction_count"] == 0, "B2_S0_C_OWNERSHIP", "no surrogate quotient")
    scandidates, mcandidates = indexed(sox["candidate_terms"]), indexed(mox["candidate_terms"])
    require(set(scandidates) == set(mcandidates) and scandidates, "B2_S0_C_OWNERSHIP", "actual candidate currents")
    for candidate in scandidates:
        require(scandidates[candidate]["kind"] == mcandidates[candidate]["kind"] and scandidates[candidate]["action_root"] == mcandidates[candidate]["action_root"], "B2_S0_C_OWNERSHIP", candidate)
        equal_expr(scandidates[candidate]["concrete_expression"], mcandidates[candidate]["concrete_expression"], "B2_S0_C_OWNERSHIP", candidate)

    sw, mw = sym["shared_input_whitelist"], math["shared_input_whitelist"]
    require(sw == mw and "B1.indexed_cells.*.M_XX_p0_known" in sw["forbidden_shared_reduced_objects"], "B2_S0_C_WHITELIST", "raw-only sharing")
    allowed_routes = set(sw["route_partition_reconstruction"])
    forbidden_routes = set(sw["forbidden_shared_reduced_objects"])
    require(allowed_routes and allowed_routes.isdisjoint(forbidden_routes), "B2_S0_C_WHITELIST", "executable allowed/forbidden disjointness")
    slots, datum_rows, final_dispositions = validate_datum_banks(sym, math, phase, mechanics, schema, so)
    dispositions = {row["slot"]: row["disposition_candidate"] for row in datum_rows.values()}
    independent_floor = expanded_schema(schema, sorted(so), dispositions)
    sobl, mobl = sym["minimum_obligation_floor"], math["minimum_obligation_floor"]
    require(sobl["schema_sha256"] == mobl["schema_sha256"] == sha256_file(schema_path), "B2_S0_C_OBLIGATIONS", "v48 schema pin")
    require(set(sobl["expanded_records"]) == set(mobl["expanded_records"]) == set(independent_floor) and len(sobl["expanded_records"]) == len(mobl["expanded_records"]) == len(independent_floor) and set(sobl["stage0_datum_slots"]) == set(mobl["stage0_datum_slots"]) == set(slots), "B2_S0_C_OBLIGATIONS", "complete product and datum floor")
    require(sobl["expanded_record_count"] == mobl["expanded_record_count"] == len(independent_floor), "B2_S0_C_OBLIGATIONS", "floor count")
    require(sym["trajectory_representation"] == math["trajectory_representation"] and sym["validity_domain_construction"] == math["validity_domain_construction"], "B2_S0_C_VALIDITY", "trajectory/validity")

    required_checks = ["B2_S0_ACTION_PARSE", "B2_S0_ACTION_HESSIAN", "B2_S0_OPERATOR_INCIDENCE", "B2_S0_BRANE_NORMALIZATION", "B2_S0_NATIVE_NOETHER", "B2_S0_U1_CURRENT", "B2_S0_COMOVING", "B2_S0_BALANCE_HONEST_EXIT", "B2_S0_BALANCE_SLOTS", "B2_S0_G9_ROUTER", "B2_S0_RESOLVENT", "B2_S0_NATIVE_CONTROL", "B2_S0_RESTRICTION", "B2_S0_OWNERSHIP", "B2_S0_WHITELIST", "B2_S0_DISPOSITION_FLOOR", "B2_S0_UNAVAILABILITY_WITNESS", "B2_S0_OBLIGATIONS"]
    for check in required_checks:
        require(sym["checks"].get(check) == "PASS" and math["checks"].get(check) == "PASS", check, "engine-local executable premise")
    return [
        "complete dual-engine frozen-base action Hessian with live chi-u block",
        "full-dimensional native Noether/current derivations",
        "honest Reynolds/typed-root integrated-balance constructive exits",
        "dimensional endpoint kernels and actual-sector nonzero controls",
        "every-branch Green sensitivity matrix",
        "actual-expression co-moving restriction and honest ownership exit",
        "v48 exact datum floor with directive producer maps",
        "typed unavailability witnesses and constructive challenges",
        "v48 overwrite/first-construction unresolved guard",
        "expanded v48 full-key obligation floor",
    ], final_dispositions


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--sympy", type=Path, required=True)
    parser.add_argument("--mathematica", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    try:
        config = load_yaml(args.input)
        sym, math = load_yaml(args.sympy), load_yaml(args.mathematica)
        checks, dispositions = compare(sym, math, config, args.input)
        frozen_keys = [
            "causal_definition", "fourier_convention", "complete_action_second_variation", "operator_inventory",
            "operator_action_incidence", "operator_action_incidence_residual", "native_noether_derivations",
            "current_derivations", "integrated_balance_identities", "endpoint_resolvent_cells",
            "functional_analytic_test_data", "green_sensitivity_matrix", "restriction_map", "ownership_convention",
            "stage0_datum_bank", "shared_input_whitelist", "minimum_obligation_floor", "trajectory_representation",
            "validity_domain_construction", "route_separation",
        ]
        frozen = {key: sym[key] for key in frozen_keys}
        frozen["stage0_datum_dispositions"] = dispositions
        frozen["stage0_datum_disposition_sha256"] = digest(dispositions)
        artifact = {
            "schema_version": AGREEMENT_SCHEMA, "directive_version": 48, "startup_contract_commit": ANCHOR,
            "status": "ENGINE_AGREE", "checks": checks, "check_count": len(checks),
            "sympy_semantic_sha256": sym["semantic_payload_sha256"],
            "mathematica_semantic_sha256": math["semantic_payload_sha256"],
            "frozen_math": frozen, "frozen_math_sha256": digest(frozen),
        }
        dump_yaml(args.output, artifact)
        print(f"B2_STAGE0_COMPARE: PASS checks={len(checks)} operators={len(frozen['operator_inventory'])} slots={len(dispositions)}")
        return 0
    except Exception as exc:
        print(str(exc))
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
