#!/usr/bin/env python3
"""Out-of-process stage-0 mutations through production assertions."""

from __future__ import annotations

import argparse
import concurrent.futures
import copy
import os
import subprocess
import sys
from pathlib import Path
from typing import Any, Callable

import sympy as sp

from u1_body_b2_common import HERE, ROOT, digest, dump_yaml, load_yaml, require, sha256_file


Mutation = Callable[[dict[str, Any], dict[str, Any]], None]


def run(command: list[str], timeout: int = 180) -> subprocess.CompletedProcess[str]:
    return subprocess.run(command, cwd=HERE, check=False, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, timeout=timeout)


def expression(raw: Any) -> sp.Expr:
    return sp.sympify(str(raw).replace("^", "**"), locals={"I": sp.I, "delta": sp.Function("delta")})


def exact_equal(left: Any, right: Any) -> bool:
    return sp.simplify(expression(left) - expression(right)) == 0


def noether_row(artifact: dict[str, Any], ancestor: str) -> dict[str, Any] | None:
    return next((row for row in artifact["native_noether_derivations"] if ancestor in row["action_ancestry"]), None)


def current_structure(artifact: dict[str, Any]) -> dict[str, Any]:
    sectors = artifact["current_derivations"]["sectors"]
    return {
        "mass": sectors["mass"]["comoving_current_components"],
        "momentum": sectors["momentum"]["comoving_current_components"],
        "energy": sectors["energy"]["comoving_current_components"],
    }


def verify_exact_contribution_removal(baseline: dict[str, Any], altered: dict[str, Any], ancestor: str, engine: str) -> dict[str, Any]:
    before = noether_row(baseline, ancestor); after = noether_row(altered, ancestor)
    require(before is not None and after is not None, "MUTATION_NOOP", f"{ancestor}:{engine}:Noether path")
    if before["status"].startswith("DERIVED_FROM"):
        contribution = before["action_contributions"][ancestor]
        after_contribution = after["action_contributions"][ancestor]
        leaf_delta = sp.simplify(expression(after_contribution["L"]) - expression(contribution["L"]))
        require(leaf_delta != 0 and exact_equal(after_contribution["L"], "0"), "MUTATION_NOOP", f"{ancestor}:{engine}:exact leaf contribution removal")
        require(exact_equal(expression(after["native_lagrangian_density"]) - expression(before["native_lagrangian_density"]), leaf_delta), "MUTATION_NOOP", f"{ancestor}:{engine}:L exact delta")
        tensor_before = before["canonical_tensor_components"]["T_mu__nu"]
        tensor_after = after["canonical_tensor_components"]["T_mu__nu"]
        tensor_changed = []
        for i in range(len(tensor_before)):
            for j in range(len(tensor_before)):
                tensor_changed.append(sp.simplify(expression(tensor_after[i][j]) - expression(tensor_before[i][j])))
        require(any(value != 0 for value in tensor_changed), "MUTATION_NOOP", f"{ancestor}:{engine}:tensor exact contribution path")
        current_changed = current_structure(baseline) != current_structure(altered)
        require(current_changed, "MUTATION_NOOP", f"{ancestor}:{engine}:real current path")
        return {"engine": engine, "exit_code": 0, "certification": "isolated_native_term_removal_through_completed_production_path", "ancestor": ancestor, "removed_expression_is_zero": True, "lagrangian_contribution_delta": str(leaf_delta), "tensor_component_count": len(tensor_changed), "tensor_structure_changed": True, "current_structure_changed": True, "output_structure_or_exact_contribution_changed": True}
    after_source = after.get("native_source_expression", "0")
    require(before["status"].startswith("UNRESOLVED") and before["native_source_expression"] != after_source and after["status"] in {before["status"], "COMPUTED_ZERO(source_ablation)"}, "MUTATION_NOOP", f"{ancestor}:{engine}:OPEN native source expression path")
    delta = sp.simplify(expression(after_source) - expression(before["native_source_expression"]))
    require(delta != 0 and exact_equal(after_source, "0"), "MUTATION_NOOP", f"{ancestor}:{engine}:OPEN exact source removal")
    return {"engine": engine, "exit_code": 0, "certification": "isolated_OPEN_source_removal_through_completed_production_path", "ancestor": ancestor, "removed_expression_is_zero": True, "native_source_expression_delta": str(delta), "output_structure_or_exact_contribution_changed": True}


def source_ancestor_probes(args: argparse.Namespace, config: dict[str, Any], base_sym: dict[str, Any], base_math: dict[str, Any], batch: tuple[int, int] | None = None) -> list[dict[str, Any]]:
    phase = load_yaml(ROOT / config["contracts"]["phase_a_inputs"])
    ancestors = [row["id"] for row in phase["action_terms"]] + ["geon_core_bundle"]
    if batch is not None:
        index, total = batch
        ancestors = [ancestor for position, ancestor in enumerate(ancestors) if position % total == index]
    rows = []
    for ancestor in ancestors:
        case_dir = args.work / "source_ancestor_probes" / ancestor; case_dir.mkdir(parents=True, exist_ok=True)
        mutated = copy.deepcopy(phase)
        if ancestor == "geon_core_bundle":
            mutated["field_records"] = [row for row in mutated["field_records"] if row["id"] != ancestor]
        else:
            target = next(row for row in mutated["action_terms"] if row["id"] == ancestor)
            target["expression"] = "0"
            require(target["expression"] == "0", "MUTATION_NOOP", f"{ancestor}:source removal did not land")
        phase_mut = case_dir / "phase_a_mutated.yaml"; sym_out = case_dir / "sympy.yaml"; math_out = case_dir / "mathematica.yaml"; math_assert = case_dir / "mathematica_assert.log"
        math_assert.unlink(missing_ok=True)
        dump_yaml(phase_mut, mutated)
        commands = [
            [sys.executable, str(HERE / "u1_body_b2_stage0_sympy.py"), "--input", str(args.input), "--phase-input", str(phase_mut), "--mutation-projection", "--output", str(sym_out)],
            ["wolframscript", "-file", str(HERE / "u1_body_b2_stage0_dual.wl"), "--input", str(args.input), "--phase-input", str(phase_mut), "--mutation-projection", "--assert-log", str(math_assert), "--output", str(math_out)],
        ]
        with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
            sym_proc, math_proc = list(pool.map(lambda command: run(command, 240), commands))
        engine_rows = []
        for engine, proc, path, baseline in [("SymPy", sym_proc, sym_out, base_sym), ("Mathematica", math_proc, math_out, base_math)]:
            if proc.returncode != 0:
                own_markers = [
                    "ASSERT_FAIL:B2_S0_ACTION_HESSIAN:", "ASSERT_FAIL:B2_S0_OPERATOR_INCIDENCE:",
                    "ASSERT_FAIL:B2_S0_BRANE_NORMALIZATION:", "ASSERT_FAIL:B2_S0_NATIVE_NOETHER:",
                    "ASSERT_FAIL:B2_S0_U1_CURRENT:", "ASSERT_FAIL:B2_S0_RESOLVENT:",
                    "ASSERT_FAIL:B2_S0_NATIVE_CONTROL:",
                ]
                assertion_evidence = proc.stdout
                if engine == "Mathematica" and math_assert.is_file():
                    assertion_evidence += math_assert.read_text(encoding="utf-8")
                marker = next((item for item in own_markers if item in assertion_evidence), None)
                require(proc.returncode == 1 and marker is not None, "MUTATION_NOOP", f"{ancestor}:{engine}:source omission wrong failure:{assertion_evidence[-500:]}")
                engine_rows.append({"engine": engine, "exit_code": proc.returncode, "certification": "isolated_source_removal_fired_engine_executable_premise", "ancestor": ancestor, "removed_expression_is_zero": True, "expected_assert": marker.split(":", 2)[1], "killed_at_own_assert": True, "output_structure_or_exact_contribution_changed": True})
                continue
            altered = load_yaml(path)
            if ancestor == "geon_core_bundle":
                before_ids = {row["id"] for row in baseline["operator_inventory"]}; after_ids = {row["id"] for row in altered["operator_inventory"]}
                require("geon_open" in before_ids and "geon_open" not in after_ids, "MUTATION_NOOP", f"{ancestor}:{engine}:operator branch removal")
                require(baseline["minimum_obligation_floor"]["expanded_record_set_sha256"] != altered["minimum_obligation_floor"]["expanded_record_set_sha256"], "MUTATION_NOOP", f"{ancestor}:{engine}:full floor propagation")
                engine_rows.append({"engine": engine, "exit_code": 0, "certification": "completed_operator_branch_and_full_floor_change", "destroyed_structure": "geon_open", "output_structure_or_exact_contribution_changed": True})
            else:
                before_digest = digest({"operators": baseline["operator_inventory"], "noether": baseline["native_noether_derivations"], "currents": baseline["current_derivations"]})
                after_digest = digest({"operators": altered["operator_inventory"], "noether": altered["native_noether_derivations"], "currents": altered["current_derivations"]})
                require(before_digest != after_digest, "MUTATION_NOOP", f"{ancestor}:{engine}:isolated removal must change output")
                after_row = noether_row(altered, ancestor)
                if after_row is None:
                    engine_rows.append({"engine": engine, "exit_code": 0, "certification": "isolated_source_removal_changed_operator_Noether_structure", "ancestor": ancestor, "removed_expression_is_zero": True, "output_structure_or_exact_contribution_changed": True})
                else:
                    engine_rows.append(verify_exact_contribution_removal(baseline, altered, ancestor, engine))
        require(len(engine_rows) == 2, "MUTATION_NOOP", f"{ancestor}:dual route")
        rows.append({"id": f"physics_ancestor::{ancestor}", "kind": "dual_route_completed_source_to_output_probe", "ancestor": ancestor, "source_level": True, "mutation": "isolated source-term removal", "requires_completed_production_path_or_own_engine_assert": True, "engines": engine_rows})
    return rows


def comparator_cases() -> list[tuple[str, str, Mutation]]:
    def each(s: dict[str, Any], m: dict[str, Any], mutation: Callable[[dict[str, Any]], None]) -> None:
        mutation(s); mutation(m)

    def datum(artifact: dict[str, Any], slot: str) -> dict[str, Any]:
        return next(row for row in artifact["stage0_datum_bank"]["records"] if row["slot"] == slot)

    def schema(s: dict[str, Any], _: dict[str, Any]) -> None: s["schema_version"] = "CORRUPTED"
    def independence(s: dict[str, Any], m: dict[str, Any]) -> None: m["independent_representation"] = s["independent_representation"]
    def input_digest(s: dict[str, Any], m: dict[str, Any]) -> None: each(s, m, lambda row: row.__setitem__("input_sha256", "0" * 64))
    def causal(s: dict[str, Any], m: dict[str, Any]) -> None: each(s, m, lambda row: row["causal_definition"].__setitem__("prescription", "advanced"))
    def fourier(s: dict[str, Any], m: dict[str, Any]) -> None: each(s, m, lambda row: row["fourier_convention"]["derivative_map"].__setitem__("dt", "+I*omega"))
    def full_hessian(s: dict[str, Any], m: dict[str, Any]) -> None: each(s, m, lambda row: row["complete_action_second_variation"].__setitem__("wall_gate_block_status", "UNRESOLVED(coupled_chi_u_block)"))
    def operators(s: dict[str, Any], m: dict[str, Any]) -> None: each(s, m, lambda row: row["operator_inventory"].pop())
    def incidence(s: dict[str, Any], m: dict[str, Any]) -> None: each(s, m, lambda row: row.__setitem__("operator_action_incidence_residual", ["removed_root"]))
    def currents(s: dict[str, Any], m: dict[str, Any]) -> None: each(s, m, lambda row: row["current_derivations"]["sectors"]["mass"]["comoving_current_components"].__setitem__(0, "0"))
    def noether(s: dict[str, Any], m: dict[str, Any]) -> None: each(s, m, lambda row: next(x for x in row["native_noether_derivations"] if x["id"] == "h_scalar")["computed_component_residuals"].__setitem__(0, "1"))
    def balance_route_a(s: dict[str, Any], m: dict[str, Any]) -> None: each(s, m, lambda row: row["integrated_balance_identities"]["sectors"]["energy"]["route_A_native_reynolds"].__setitem__("source_omission_premises", []))
    def balance_route_b(s: dict[str, Any], m: dict[str, Any]) -> None: each(s, m, lambda row: row["integrated_balance_identities"]["sectors"]["energy"]["route_B_authenticated_typed_roots"].__setitem__("source_omission_premises", []))
    def balance_surrogate(s: dict[str, Any], m: dict[str, Any]) -> None: each(s, m, lambda row: row["integrated_balance_identities"]["sectors"]["mass"].__setitem__("surrogate_symbol_list_used", True))
    def g9(s: dict[str, Any], m: dict[str, Any]) -> None: each(s, m, lambda row: row["integrated_balance_identities"]["sectors"]["energy"]["g9_stage0_router_output"].__setitem__("verdict", "OK(exact)"))
    def kernel_ast(s: dict[str, Any], m: dict[str, Any]) -> None: each(s, m, lambda row: next(x for x in row["endpoint_resolvent_cells"] if x["cell"] == "brane_shear_transverse|E5")["retarded_green_operator"]["executable_kernel"]["kernel_ast"].__setitem__("normalization", "2/muR"))
    def e4_domain(s: dict[str, Any], m: dict[str, Any]) -> None: each(s, m, lambda row: next(x for x in row["endpoint_resolvent_cells"] if x["cell"] == "brane_shear_transverse|E4").__setitem__("endpoint_domain_kind", "Dirichlet"))
    def control_tensor(s: dict[str, Any], m: dict[str, Any]) -> None: each(s, m, lambda row: next(x for x in row["endpoint_resolvent_cells"] if x["cell"] == "brane_shear_transverse|E1")["known_nonzero_control"]["production_route"].__setitem__("flux", "0"))
    def control_oracle(s: dict[str, Any], m: dict[str, Any]) -> None: each(s, m, lambda row: next(x for x in row["endpoint_resolvent_cells"] if x["cell"] == "brane_shear_transverse|E1")["known_nonzero_control"]["oracle_route"].__setitem__("flux", "0"))
    def functional(s: dict[str, Any], m: dict[str, Any]) -> None: each(s, m, lambda row: row["functional_analytic_test_data"][0].__setitem__("test_space", "C_c_infinity"))
    def sensitivity(s: dict[str, Any], m: dict[str, Any]) -> None: each(s, m, lambda row: next(x for x in row["green_sensitivity_matrix"] if x["operator"] == "brane_normal_local")["mutations"][0].__setitem__("operator_derivative", "0"))
    def restriction(s: dict[str, Any], m: dict[str, Any]) -> None: each(s, m, lambda row: row["restriction_map"]["substitution"].__setitem__("p_x", 1))
    def ownership(s: dict[str, Any], m: dict[str, Any]) -> None: each(s, m, lambda row: row["ownership_convention"].__setitem__("surrogate_IBP_relation_used_as_current_improvement", True))
    def whitelist(s: dict[str, Any], m: dict[str, Any]) -> None: each(s, m, lambda row: row["shared_input_whitelist"]["route_partition_reconstruction"].append("B1.indexed_cells.*.M_XX_p0_known"))
    def obligations(s: dict[str, Any], m: dict[str, Any]) -> None: each(s, m, lambda row: row["minimum_obligation_floor"]["expanded_records"].pop())
    def disposition_omission(s: dict[str, Any], m: dict[str, Any]) -> None: each(s, m, lambda row: row["stage0_datum_bank"]["records"].pop())
    def producer_map(s: dict[str, Any], m: dict[str, Any]) -> None: each(s, m, lambda row: datum(row, "full_action_second_variation").__setitem__("producer_ids", ["action::bulk_berry"]))
    def closure(s: dict[str, Any], m: dict[str, Any]) -> None: each(s, m, lambda row: datum(row, "ownership_current_improvement_relation")["unavailability_witness"]["enumerated_committed_inputs"].append("builder::invented"))
    def restore(s: dict[str, Any], m: dict[str, Any]) -> None: each(s, m, lambda row: datum(row, "ownership_current_improvement_relation")["unavailability_witness"]["counterfactual_restore_mutation"].__setitem__("producer_to_restore", "not_a_producer"))
    def challenge_kind(s: dict[str, Any], m: dict[str, Any]) -> None: each(s, m, lambda row: datum(row, "ownership_current_improvement_relation")["derivability_challenge"].__setitem__("kind", "operator_domain_well_posedness_failure"))
    def derived_typecheck(s: dict[str, Any], m: dict[str, Any]) -> None:
        def mutate(row: dict[str, Any]) -> None:
            target = datum(row, "causal_retarded_definition"); target["candidate_typecheck"]["identity_domain_membership"] = False
        each(s, m, mutate)

    def refute_slot(slot: str, first_time: bool = False) -> Mutation:
        def mutate(s: dict[str, Any], m: dict[str, Any]) -> None:
            def one(artifact: dict[str, Any]) -> None:
                target = datum(artifact, slot)
                target["disposition_candidate"] = "UNRESOLVED"
                target["unresolved_tag"] = "UNRESOLVED"
                target["candidate_ref"] = target.get("candidate_ref") or "constructive_candidate"
                target["candidate_is_well_typed"] = True
                target["candidate_typecheck"] = {"candidate_present": True, "type_matches_required_type": True, "dimensions_match_required_dimensions": True, "identity_domain_membership": True, "computed_result": True}
                target["defining_predicate_result"] = "PASS"
                target["first_time_construction_fixture"] = first_time
                if "unavailability_witness" not in target:
                    template = copy.deepcopy(datum(artifact, "ownership_current_improvement_relation")["unavailability_witness"])
                    template.update({"witness_id": "witness::" + slot, "datum_id": slot, "required_type": target["required_type"], "required_dimensions": target["required_dimensions"], "authoritative_roots": target["producer_ids"], "enumerated_committed_inputs": target["committed_input_closure"], "directive_generated_producer_rule": target["producer_rule"], "producer_census": [{"producer": producer, "in_closure": producer in target["committed_input_closure"], "type_compatible": False if producer in target["committed_input_closure"] else None} for producer in target["producer_ids"]]})
                    target["unavailability_witness"] = template
                challenge = target["derivability_challenge"]
                challenge.update({"challenge_id": "challenge::" + slot, "candidate_ref": target["candidate_ref"], "candidate_is_well_typed": True, "candidate_typecheck": copy.deepcopy(target["candidate_typecheck"]), "defining_predicate_result": "PASS", "terminal": "REFUTED(well-typed PASS candidate)", "constructive_attempt_nonempty": True})
            each(s, m, one)
        return mutate

    cases: list[tuple[str, str, Mutation]] = [
        ("schema_tooth", "B2_S0_C_SCHEMA", schema), ("independence_tooth", "B2_S0_C_INDEPENDENCE", independence),
        ("input_tooth", "B2_S0_C_INPUT", input_digest), ("causal_tooth", "B2_S0_C_CAUSAL", causal),
        ("fourier_tooth", "B2_S0_C_FOURIER", fourier), ("full_Hessian_tooth", "B2_S0_C_FULL_HESSIAN", full_hessian),
        ("operators_tooth", "B2_S0_C_OPERATORS", operators), ("incidence_tooth", "B2_S0_C_INCIDENCE", incidence),
        ("Noether_tooth", "B2_S0_C_NOETHER", noether), ("currents_tooth", "B2_S0_C_CURRENTS", currents),
        ("balance_route_A_source_omission", "B2_S0_C_BALANCE_ROUTE_A", balance_route_a),
        ("balance_route_B_source_omission", "B2_S0_C_BALANCE_ROUTE_B", balance_route_b),
        ("balance_surrogate_tooth", "B2_S0_C_BALANCES", balance_surrogate), ("G9_tooth", "B2_S0_C_G9_ROUTER", g9),
        ("resolvent_kernel_tooth", "B2_S0_C_RESOLVENT", kernel_ast), ("E4_domain_tooth", "B2_S0_C_RESOLVENT", e4_domain),
        ("native_control_tensor_route", "B2_S0_C_NATIVE_CONTROL_TENSOR_ROUTE", control_tensor),
        ("native_control_oracle_route", "B2_S0_C_NATIVE_CONTROL_ORACLE_ROUTE", control_oracle),
        ("functional_tooth", "B2_S0_C_FUNCTIONAL", functional), ("sensitivity_tooth", "B2_S0_C_SENSITIVITY", sensitivity),
        ("restriction_tooth", "B2_S0_C_RESTRICTION", restriction), ("ownership_tooth", "B2_S0_C_OWNERSHIP", ownership),
        ("whitelist_tooth", "B2_S0_C_WHITELIST", whitelist), ("obligations_tooth", "B2_S0_C_OBLIGATIONS", obligations),
        ("disposition_floor_omission", "B2_S0_C_DISPOSITION_FLOOR", disposition_omission),
        ("producer_map_tooth", "B2_S0_C_PRODUCER_MAP", producer_map), ("witness_closure_tooth", "B2_S0_C_WITNESS_CLOSURE", closure),
        ("witness_restore_tooth", "B2_S0_C_WITNESS_RESTORE", restore), ("challenge_kind_tooth", "B2_S0_C_CHALLENGE", challenge_kind),
        ("derived_typecheck_tooth", "B2_S0_C_DERIVED_TYPECHECK", derived_typecheck),
        ("overwrite_Hessian_anti_dodge", "stage0_unresolved_refuted", refute_slot("full_action_second_variation")),
        ("first_time_Hessian_construction", "stage0_unresolved_refuted", refute_slot("full_action_second_variation", first_time=True)),
    ]
    for class_name, slot in [
        ("operator_Hessian", "operator_hessian_block|operator_branch=brane_shear_transverse"),
        ("native_current", "native_mass_current_law"),
        ("native_Noether", "native_noether_flux|operator_branch=h_scalar"),
        ("retarded_Green", "retarded_green_operator|operator_branch=brane_shear_transverse|endpoint=E1"),
        ("causal_definition", "causal_retarded_definition"),
        ("trajectory_definition", "trajectory_frequency_representation"),
        ("truncation_definition", "truncation_frequency_variable"),
        ("authority_absence_canary", "ownership_current_improvement_relation"),
        ("domain_failure_canary", "retarded_green_operator|operator_branch=gnls_density_phase|endpoint=E1"),
        ("upstream_freeze_canary", "mass_sector_tolerance"),
    ]:
        cases.append((f"derivability_class::{class_name}", "stage0_unresolved_refuted", refute_slot(slot)))
    return cases


def engine_assertion_cases() -> list[tuple[str, str, Mutation]]:
    keys = [
        "B2_S0_ACTION_PARSE", "B2_S0_ACTION_HESSIAN", "B2_S0_OPERATOR_INCIDENCE", "B2_S0_BRANE_NORMALIZATION",
        "B2_S0_NATIVE_NOETHER", "B2_S0_U1_CURRENT", "B2_S0_COMOVING", "B2_S0_BALANCE_HONEST_EXIT",
        "B2_S0_BALANCE_SLOTS", "B2_S0_G9_ROUTER", "B2_S0_RESOLVENT", "B2_S0_NATIVE_CONTROL",
        "B2_S0_RESTRICTION", "B2_S0_OWNERSHIP", "B2_S0_WHITELIST", "B2_S0_DISPOSITION_FLOOR",
        "B2_S0_UNAVAILABILITY_WITNESS", "B2_S0_OBLIGATIONS",
    ]
    result = []
    for key in keys:
        def mutate(s: dict[str, Any], m: dict[str, Any], check: str = key) -> None:
            for row in (s, m): row["checks"][check] = "FAIL"
        result.append((f"engine_assertion::{key}", key, mutate))
    return result


def run_comparator_cases(args: argparse.Namespace, base_sym: dict[str, Any], base_math: dict[str, Any], batch: tuple[int, int] | None = None) -> list[dict[str, Any]]:
    prepared = [(case_id, tooth) for case_id, tooth, _ in comparator_cases() + engine_assertion_cases()]
    if batch is not None:
        index, total = batch
        prepared = [item for position, item in enumerate(prepared) if position % total == index]

    def execute(item: tuple[str, str]) -> dict[str, Any]:
        case_id, tooth = item
        proc = run([sys.executable, str(HERE / "u1_body_b2_stage0_mutation_harness.py"), "--input", str(args.input), "--sympy", str(args.sympy), "--mathematica", str(args.mathematica), "--case", case_id], 240)
        marker = f"ASSERT_FAIL:{tooth}:"
        require(proc.returncode == 1 and marker in proc.stdout, "MUTATION_NOOP", f"{case_id}:expected {marker}:code={proc.returncode}:log={proc.stdout[-500:]}")
        return {"id": case_id, "kind": "out_of_process_production_comparator_function_tooth", "harness": "mutation applied after real artifact parse and before unchanged compare() call", "expected_assert": tooth, "exit_code": proc.returncode, "killed_at_own_assert": True, "log_tail": proc.stdout.strip().splitlines()[-1]}

    # These are independent production-comparator processes; Wolfram execution
    # remains strictly serial in source_ablations above.
    with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
        rows = list(pool.map(execute, prepared))
    return rows


def trace_semantics_teeth(work: Path) -> list[dict[str, Any]]:
    fixtures = [
        ("relative_dirfd_protected_rewrite", "B2_A1_PROTECTED_REWRITE"),
        ("preexisting_generated_laundering", "B2_TRACE_GENERATED_CREATE"),
    ]
    rows = []
    for fixture, tooth in fixtures:
        case = work / "trace_semantics" / fixture
        output = case / "unused.yaml"
        proc = run([sys.executable, str(HERE / "u1_body_b2_trace_audit.py"), "--mutation-fixture", fixture, "--fixture-root", str(case), "--output", str(output)], 60)
        marker = f"ASSERT_FAIL:{tooth}:"
        require(proc.returncode == 1 and marker in proc.stdout, "MUTATION_NOOP", f"{fixture}:expected {marker}:{proc.stdout[-500:]}")
        rows.append({"id": f"trace_semantics::{fixture}", "kind": "out_of_process_actual_trace_parser_tooth", "production_assertion_module": "u1_body_b2_trace_audit.py", "expected_assert": tooth, "exit_code": proc.returncode, "killed_at_own_assert": True, "log_tail": proc.stdout.strip().splitlines()[-1]})
    return rows


def parse_batch(raw: str | None) -> tuple[int, int] | None:
    if raw is None:
        return None
    index, total = (int(value) for value in raw.split("/", 1))
    require(total > 0 and 0 <= index < total, "MUTATION_NOOP", f"invalid batch:{raw}")
    return index, total


def batch_artifact(kind: str, batch: str, rows: list[dict[str, Any]]) -> dict[str, Any]:
    return {"schema_version": "U1_PHASE_B2_STAGE0_MUTATION_BATCH_V1", "status": "PASS", "kind": kind, "batch": batch, "case_count": len(rows), "cases": rows}


def aggregate_results(paths: list[Path]) -> dict[str, Any]:
    rows = []
    for path in paths:
        artifact = load_yaml(path)
        require(artifact["schema_version"] == "U1_PHASE_B2_STAGE0_MUTATION_BATCH_V1" and artifact["status"] == "PASS", "MUTATION_NOOP", f"batch:{path}")
        rows.extend(artifact["cases"])
    ids = [row["id"] for row in rows]
    require(len(ids) == len(set(ids)), "MUTATION_NOOP", "duplicate aggregated mutation id")
    physics = [row for row in rows if row["kind"] == "dual_route_completed_source_to_output_probe"]
    comparator = [row for row in rows if row["kind"] == "out_of_process_production_comparator_function_tooth"]
    trace = [row for row in rows if row["kind"] == "out_of_process_actual_trace_parser_tooth"]
    expected_comparator = len(comparator_cases()) + len(engine_assertion_cases())
    require(len(physics) == 16 and len(comparator) == expected_comparator and len(trace) == 2, "MUTATION_NOOP", f"aggregate coverage:{len(physics)},{len(comparator)}/{expected_comparator},{len(trace)}")
    return {"schema_version": "U1_PHASE_B2_STAGE0_MUTATIONS_V6", "directive_version": 48, "status": "PASS", "out_of_process_comparator_cases": True, "engines_mutation_unaware": True, "MUTATION_NOOP_sentinel": "armed", "case_count": len(rows), "physics_ancestor_probe_count": len(physics), "comparator_tooth_count": len(comparator), "trace_semantics_tooth_count": len(trace), "completed_production_path_or_own_engine_assert_required": True, "one_case_per_engine_assertion": True, "v48_overwrite_and_first_construction_teeth": True, "derivability_class_canaries": True, "cases": sorted(rows, key=lambda row: row["id"])}


def main() -> int:
    parser = argparse.ArgumentParser(); parser.add_argument("--input", type=Path); parser.add_argument("--sympy", type=Path); parser.add_argument("--mathematica", type=Path); parser.add_argument("--work", type=Path); parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--source-batch"); parser.add_argument("--comparator-batch"); parser.add_argument("--trace-only", action="store_true"); parser.add_argument("--aggregate-results", type=Path, nargs="+"); args = parser.parse_args()
    try:
        args.output = args.output.resolve()
        if args.aggregate_results:
            artifact = aggregate_results([path.resolve() for path in args.aggregate_results])
            dump_yaml(args.output, artifact); print(f"B2_STAGE0_MUTATIONS: PASS aggregate_cases={artifact['case_count']}"); return 0
        require(all(path is not None for path in [args.input, args.sympy, args.mathematica, args.work]), "MUTATION_NOOP", "non-aggregate paths")
        args.input, args.sympy, args.mathematica, args.work = (path.resolve() for path in (args.input, args.sympy, args.mathematica, args.work)); args.work.mkdir(parents=True, exist_ok=True)
        modes = sum([args.source_batch is not None, args.comparator_batch is not None, args.trace_only])
        require(modes == 1, "MUTATION_NOOP", "exactly one batch mode")
        config = load_yaml(args.input); base_sym, base_math = load_yaml(args.sympy), load_yaml(args.mathematica)
        if args.source_batch is not None:
            rows = source_ancestor_probes(args, config, base_sym, base_math, parse_batch(args.source_batch)); kind, batch = "physics_ancestor", args.source_batch
        elif args.comparator_batch is not None:
            rows = run_comparator_cases(args, base_sym, base_math, parse_batch(args.comparator_batch)); kind, batch = "comparator_assertion", args.comparator_batch
        else:
            rows = trace_semantics_teeth(args.work); kind, batch = "trace_semantics", "0/1"
        artifact = batch_artifact(kind, batch, rows)
        dump_yaml(args.output, artifact); print(f"B2_STAGE0_MUTATION_BATCH: PASS kind={kind} cases={len(rows)} batch={batch}"); return 0
    except Exception as exc:
        print(str(exc)); return 1


if __name__ == "__main__":
    raise SystemExit(main())
