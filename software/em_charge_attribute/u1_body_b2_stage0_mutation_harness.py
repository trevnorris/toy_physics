#!/usr/bin/env python3
"""Out-of-process one-tooth harness for the production stage-0 comparator."""

from __future__ import annotations

import argparse
import copy
import sys
from pathlib import Path
from typing import Any

from u1_body_b2_common import digest, load_yaml, require
from u1_body_b2_stage0_compare import compare
from u1_body_b2_stage0_mutations import comparator_cases, engine_assertion_cases

# Comparator failure is the expected outcome of this one-tooth subprocess.
# Do not route that deliberate exception through the host's Apport hook, which
# would add an unrelated package-manager/crash-reporter stack to first use.
sys.excepthook = sys.__excepthook__


def first_time_hessian_fixture(artifact: dict[str, Any]) -> None:
    """Plant UNRESOLVED before constructing the Hessian challenge candidate.

    This deliberately does not transform or read the prior datum record.  The
    record is removed first, then the challenge is constructed from the
    engine's independent complete-action second-variation route.  Thus the
    fixture exercises the no-prior-DERIVED case rather than being a second
    spelling of the overwrite tooth.
    """
    slot = "full_action_second_variation"
    bank = artifact["stage0_datum_bank"]
    bank["records"] = [row for row in bank["records"] if row["slot"] != slot]
    require(
        all(row["slot"] != slot for row in bank["records"]),
        "MUTATION_NOOP",
        f"{artifact['engine']}:first-time fixture has no prior datum record",
    )

    second_variation = artifact["complete_action_second_variation"]
    termwise = second_variation["termwise_records"]
    action_roots = sorted(row["action_root"] for row in artifact["operator_action_incidence"])
    require(
        second_variation["status"] == "DERIVED_COMPLETE_TERMWISE_SECOND_VARIATION_WITH_NAMED_BASE_JET"
        and second_variation["wall_gate_block_status"] == "DERIVED_SYMBOLIC_FROZEN_BASE_JET_NO_COUPLING_SUPPRESSED"
        and {row["id"] for row in termwise} == set(action_roots),
        "MUTATION_NOOP",
        f"{artifact['engine']}:first-time constructive Hessian route",
    )
    wall_gate = next(row for row in termwise if row["id"] == "wall_shear_gate")
    require(
        wall_gate["status"].startswith("DERIVED_")
        and "UNRESOLVED" not in str(wall_gate["bilinear_second_variation"]),
        "MUTATION_NOOP",
        f"{artifact['engine']}:first-time chi-u candidate",
    )

    producer_ids = [f"action::{root}" for root in action_roots]
    closure = sorted(
        item
        for item in bank["committed_input_closure_inventory"]
        if item.startswith(("action::", "field::"))
    )
    census = [
        {"producer": producer, "in_closure": producer in closure, "type_compatible": False}
        for producer in producer_ids
    ]
    required_type = "symmetric_bilinear_form_on_full_committed_field_jet"
    dimensions = "action_density/field^2 blockwise"
    typecheck = {
        "candidate_present": True,
        "type_matches_required_type": True,
        "dimensions_match_required_dimensions": True,
        "identity_domain_membership": True,
        "computed_result": True,
    }
    candidate_ref = "first_time_challenge::complete_action_second_variation"
    record = {
        "slot": slot,
        "disposition_candidate": "UNRESOLVED",
        "unresolved_tag": "UNRESOLVED",
        "required_type": required_type,
        "required_dimensions": dimensions,
        "defining_predicate": "directive_identity_and_schema_type_match",
        "candidate_ref": candidate_ref,
        "candidate_is_well_typed": True,
        "candidate_typecheck": copy.deepcopy(typecheck),
        "defining_predicate_result": "PASS",
        "producer_rule": "v48:hessian<-second_variation(EVERY_committed_action_term)+incidence_residual",
        "producer_ids": producer_ids,
        "committed_input_closure": closure,
        "closure_exact_set_assert": "PASS",
        "unavailability_witness": {
            "witness_id": f"witness::{slot}",
            "datum_id": slot,
            "required_type": required_type,
            "required_dimensions": dimensions,
            "acceptance_predicate": "directive_identity_and_schema_type_match",
            "authoritative_roots": producer_ids,
            "enumerated_committed_inputs": closure,
            "complete_closure_exact_set_equal": True,
            "directive_generated_producer_rule": "v48:hessian<-second_variation(EVERY_committed_action_term)+incidence_residual",
            "producer_census": census,
            "producer_census_predicate": "forall p in Producer(slot): p notin closure or type(p) incompatible",
            "kind": "authority_census_producer_absence",
            "executable_certificate_result": "PASS",
            "diagnostic": "planted apparently-passing absence contradicted by the independent Hessian challenge",
            "missing_typed_ingredient": "full_committed_action_second_variation",
            "counterfactual_restore_mutation": {
                "ingredient_kind": "missing_input_leaf",
                "ingredient": "full_committed_action_second_variation",
                "producer_to_restore": producer_ids[0],
                "restored_type_compatible": True,
                "fixture_type": required_type,
                "fixture_dimensions": dimensions,
                "certificate_after_restore": "FAIL",
                "failed_at_own_assert": "B2_S0_WITNESS_RESTORE",
            },
        },
        "derivability_challenge": {
            "challenge_id": f"challenge::{slot}",
            "status": "REFUTED",
            "same_committed_input_closure": closure,
            "dag_separated_from_witness": True,
            "shared_only_committed_inputs": True,
            "constructive_attempt": "fresh complete-action second variation after removal of the datum record",
            "constructive_attempt_nonempty": True,
            "candidate_ref": candidate_ref,
            "candidate_schema_pinned": {"type": required_type, "dimensions": dimensions},
            "candidate_typecheck": copy.deepcopy(typecheck),
            "candidate_is_well_typed": True,
            "defining_predicate_result": "PASS",
            "terminal": "REFUTED(well-typed PASS candidate)",
        },
        "first_time_construction_fixture": {
            "no_prior_datum_record_at_challenge_construction": True,
            "candidate_source": "complete_action_second_variation",
            "wall_chi_u_block_required": True,
        },
    }
    bank["records"].append(record)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--sympy", type=Path, required=True)
    parser.add_argument("--mathematica", type=Path, required=True)
    parser.add_argument("--case", required=True)
    args = parser.parse_args()
    sym, math, config = load_yaml(args.sympy), load_yaml(args.mathematica), load_yaml(args.input)
    cases = {case_id: (tooth, mutate) for case_id, tooth, mutate in comparator_cases() + engine_assertion_cases()}
    require(args.case in cases, "MUTATION_NOOP", f"unknown comparator case:{args.case}")
    _, mutate = cases[args.case]
    before = digest({"sym": sym, "math": math})
    if args.case == "first_time_Hessian_construction":
        first_time_hessian_fixture(sym)
        first_time_hessian_fixture(math)
    else:
        mutate(sym, math)
    require(before != digest({"sym": sym, "math": math}), "MUTATION_NOOP", args.case)
    compare(sym, math, config, args.input.resolve())
    print(f"MUTATION_SURVIVED:{args.case}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
