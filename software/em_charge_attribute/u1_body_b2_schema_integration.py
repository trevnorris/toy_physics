#!/usr/bin/env python3
"""Read-only v48 stage-0 schema and downstream-shape integration gate."""

from __future__ import annotations

import argparse
import copy
import itertools
from pathlib import Path
from typing import Any

from u1_body_b2_common import ROOT, digest, dump_yaml, load_yaml, require


def expand_slots(schema: dict[str, Any], branches: list[str]) -> list[str]:
    axes = copy.deepcopy(schema["axes"])
    axes["operator_branch"] = branches
    slots = set(schema["stage0_fixed_datum_slots"])
    for spec in schema["stage0_expanded_datum_slots"]:
        for values in itertools.product(*(axes[axis] for axis in spec["axes"])):
            slots.add(spec["name"] + "|" + "|".join(f"{axis}={value}" for axis, value in zip(spec["axes"], values)))
    return sorted(slots)


def expand_floor(schema: dict[str, Any], branches: list[str], dispositions: dict[str, str]) -> list[str]:
    axes = copy.deepcopy(schema["axes"])
    axes["operator_branch"] = branches
    axes["unordered_operator_pair"] = [f"{a}__PAIR__{b}" for a, b in itertools.combinations(branches, 2)]
    rows = set(schema["fixed_products"])
    for spec in schema["expanded_products"]:
        for values in itertools.product(*(axes[axis] for axis in spec["axes"])):
            rows.add(spec["name"] + "|" + "|".join(f"{axis}={value}" for axis, value in zip(spec["axes"], values)))
    rows.update(f"stage0_datum|slot={slot}|disposition={disposition}" for slot, disposition in dispositions.items())
    return sorted(rows)


def check(agreement: dict[str, Any], config: dict[str, Any]) -> dict[str, Any]:
    require(agreement["schema_version"] == "U1_PHASE_B2_STAGE0_MATH_AGREEMENT_V6" and agreement["directive_version"] == 48 and agreement["startup_contract_commit"] == config["startup_contract_commit"], "B2_S0_SCHEMA_INTEGRATION", "v48 agreement schema/anchor")
    frozen = agreement["frozen_math"]
    branches = sorted(row["id"] for row in frozen["operator_inventory"])
    require(len(branches) == 8 and "wall_chi_u_coupled" in branches, "B2_S0_SCHEMA_INTEGRATION", "connected operator inventory")
    require(frozen["complete_action_second_variation"]["wall_gate_block_status"] == "DERIVED_SYMBOLIC_FROZEN_BASE_JET_NO_COUPLING_SUPPRESSED", "B2_S0_SCHEMA_INTEGRATION", "complete chi-u Hessian")

    balances = frozen["integrated_balance_identities"]
    require(balances["status"] == "UNRESOLVED(complete_integrated_balance_family)" and balances["router"]["fixture_controls"]["surrogate_removal_residual_used"] is False, "B2_S0_SCHEMA_INTEGRATION", "honest balance family")
    for sector in ["mass", "momentum", "energy"]:
        row = balances["sectors"][sector]
        require(row["surrogate_symbol_list_used"] is False and "canonical_terms" not in row, "B2_S0_SCHEMA_INTEGRATION", f"{sector}:no hand-entered terms")
        require(row["route_A_native_reynolds"]["source_omission_premises"] and row["route_B_authenticated_typed_roots"]["source_omission_premises"], "B2_S0_SCHEMA_INTEGRATION", f"{sector}:dual source premises")
        require(row["route_B_authenticated_typed_roots"]["complete_signed_expression"] is None and row["complete_signed_expression_residual"] is None, "B2_S0_SCHEMA_INTEGRATION", f"{sector}:no constructed zero residual")

    schema = load_yaml(ROOT / config["contracts"]["v48_obligation_schema"])
    require(schema["schema_version"] == "U1_PHASE_B2_V48_OBLIGATION_SCHEMA_V3" and schema["directive_version"] == 48, "B2_S0_SCHEMA_INTEGRATION", "v48 obligation schema")
    expected_slots = expand_slots(schema, branches)
    bank = frozen["stage0_datum_bank"]
    dispositions = {row["slot"]: row for row in frozen["stage0_datum_dispositions"]}
    require(set(dispositions) == set(expected_slots) == set(bank["expected_slots"]) == set(bank["reachable_slots"]) and len(dispositions) == bank["record_count"], "B2_S0_SCHEMA_INTEGRATION", "expected/reachable/disposition exact set")
    for slot, row in dispositions.items():
        require(row["disposition"] in {"DERIVED", "UNRESOLVED"}, "B2_S0_SCHEMA_INTEGRATION", f"{slot}:one-of")
        if row["disposition"] == "DERIVED":
            require(set(row) == {"slot", "disposition", "value_digest", "dual_engine_comparison_id"}, "B2_S0_SCHEMA_INTEGRATION", f"{slot}:DERIVED shape")
        else:
            require({"slot", "disposition", "witness_id", "challenge_id", "dual_engine_challenge_certificate_sha256"} == set(row), "B2_S0_SCHEMA_INTEGRATION", f"{slot}:UNRESOLVED shape")
    candidate_dispositions = {row["slot"]: row["disposition_candidate"] for row in bank["records"]}
    require(candidate_dispositions == {slot: row["disposition"] for slot, row in dispositions.items()}, "B2_S0_SCHEMA_INTEGRATION", "bank/final disposition agreement")

    expected_floor = expand_floor(schema, branches, candidate_dispositions)
    floor = frozen["minimum_obligation_floor"]
    required_axes = {"endpoint", "parity", "closure", "open_stratum", "p_slice", "motion", "g9_sector", "deferred_gate", "operator_branch", "unordered_operator_pair"}
    require(required_axes <= set(floor["grid_axes"]), "B2_S0_SCHEMA_INTEGRATION", "full-key axes")
    require(set(floor["expanded_records"]) == set(expected_floor) and len(floor["expanded_records"]) == floor["expanded_record_count"] == len(expected_floor), "B2_S0_SCHEMA_INTEGRATION", "independent v48 floor expansion")
    require(set(floor["stage0_datum_slots"]) == set(expected_slots) and floor["category_counts"]["stage0_datum"] == len(expected_slots), "B2_S0_SCHEMA_INTEGRATION", "stage0 datum category")

    checks = {
        "v48_agreement_schema_and_anchor": "PASS",
        "complete_chi_u_Hessian": "PASS",
        "honest_dual_route_balance_exit": "PASS",
        "stage0_expected_reachable_exact_set": "PASS",
        "DERIVED_UNRESOLVED_disposition_shapes": "PASS",
        "v48_full_floor_independent_expansion": "PASS",
    }
    return {
        "schema_version": "U1_PHASE_B2_SCHEMA_INTEGRATION_V3", "directive_version": 48,
        "status": "PASS", "mode": "read_only_pre_contract_regeneration_gate", "checks": checks,
        "operator_count": len(branches), "stage0_datum_count": len(expected_slots),
        "stage0_datum_set_sha256": digest(expected_slots), "floor_count": len(expected_floor),
        "floor_set_sha256": digest(expected_floor),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--agreement", type=Path, required=True)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    try:
        dump_yaml(args.output, check(load_yaml(args.agreement), load_yaml(args.input)))
        print("B2_SCHEMA_INTEGRATION: PASS v48")
        return 0
    except Exception as exc:
        print(str(exc))
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
