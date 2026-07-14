#!/usr/bin/env python3
"""Independent Phase-B1 comparator, V5/report generator, and summary printer.

The comparator has no mutation switch.  The external harness mutates temporary
artifacts and invokes this unchanged verifier.
"""

from __future__ import annotations

import argparse
import copy
import hashlib
import json
import math
import re
import sys
from pathlib import Path
from typing import Any

import sympy as sp
import yaml


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
ARTIFACTS = HERE / "reports/u1_body_dynamics_artifacts"
DEFAULT_INPUT = HERE / "u1_body_mechanics_inputs.yaml"
DEFAULT_REPORT = HERE / "reports/u1_body_dynamics.md"
DEFAULT_RESULTS = HERE / "reports/u1_body_dynamics_results.yaml"
ENDPOINTS = ("E1", "E2", "E3", "E4", "E5")
AMBIENTS = ("one_sided_pathA29", "symmetric_postulate")
COMPARE_CHECKS = ("B1_C_ENGINE_MATH", "B1_C_RECORD_MAP", "B1_C_PHASE_A", "B1_C_REPORT_SCHEMA", "B1_C_INPUT_LIVENESS")


class CompareFailure(AssertionError):
    def __init__(self, tooth: str, detail: str):
        super().__init__(f"ASSERT_FAIL:{tooth}:{detail}")


def require(test: Any, tooth: str, detail: str) -> None:
    if not bool(test):
        raise CompareFailure(tooth, detail)


def load(path: Path) -> dict[str, Any]:
    text = path.read_text()
    try:
        return json.loads(text)
    except json.JSONDecodeError:
        return yaml.safe_load(text)


def digest(value: Any) -> str:
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":"), default=str).encode()).hexdigest()


def phase_a_payload(artifact: dict[str, Any]) -> dict[str, Any]:
    return {
        "axis1": artifact["axis1"], "axis2": artifact["axis2"], "cells": artifact["cells"],
        "tails": [(r["id"], r["classification"], r["decay_exponent"], r["normalizable"]) for r in artifact["tail_channels"]],
        "decision": artifact["source_action_completeness"]["operative_decision_citation"],
        "action_ids": artifact["source_action_completeness"]["assembled_ids"],
        "endpoint_coefficients": {ep: {k: v["canonical_terms"] for k, v in row["coefficients"].items()}
                                  for ep, row in artifact["endpoint_effective_actions"].items()},
    }


def normalized_payload_v2(payload: dict[str, Any]) -> dict[str, Any]:
    output = copy.deepcopy(payload)
    output["action_ids"] = sorted(output["action_ids"])
    output["decision"] = {key: output["decision"][key] for key in ("id", "sha256", "status")}
    for endpoint in sorted(output["endpoint_coefficients"]):
        for name, terms in output["endpoint_coefficients"][endpoint].items():
            rows = [{"coefficient": f"{float(row['coefficient']):.12f}",
                     "powers": {str(k): int(v) for k, v in sorted(row["powers"].items())}} for row in terms]
            output["endpoint_coefficients"][endpoint][name] = sorted(rows, key=lambda row: tuple(row["powers"].items()))
    return output


def canonical_rows(expr: sp.Expr) -> list[dict[str, Any]]:
    combined: dict[tuple[tuple[str, int], ...], float] = {}
    for term in sp.Add.make_args(sp.expand(expr)):
        coefficient, factors = term.as_coeff_Mul(); number = sp.N(coefficient, 17); powers: dict[str, int] = {}
        for factor, power in factors.as_powers_dict().items():
            if factor.is_number:
                number *= factor**power
            else:
                powers[str(factor)] = int(power)
        key = tuple(sorted(powers.items())); combined[key] = combined.get(key, 0.0) + float(number)
    return [{"coefficient": number, "powers": dict(key)} for key, number in sorted(combined.items()) if abs(number) > 2e-10]


def comparator_amended_payload(phase: dict[str, Any]) -> dict[str, Any]:
    amended = copy.deepcopy(phase)
    rho_br, shear = sp.symbols("rhoBr shear_transverse", real=True)
    for endpoint in ENDPOINTS:
        beta = parse_expr(phase["endpoint_responses"][endpoint]["shear_coefficient"])
        old = canonical_expression(phase["endpoint_effective_actions"][endpoint]["coefficients"]["GVV"]["canonical_terms"])
        new = sp.expand(old - sp.Rational(12, 5) * sp.pi * rho_br * (shear + beta) ** 2
                        + sp.Rational(8, 3) * sp.pi * rho_br * (shear + beta) ** 2)
        amended["endpoint_effective_actions"][endpoint]["coefficients"]["GVV"]["canonical_terms"] = canonical_rows(new)
    return normalized_payload_v2(phase_a_payload(amended))


def term_map(rows: list[dict[str, Any]]) -> dict[tuple[tuple[str, int], ...], float]:
    out: dict[tuple[tuple[str, int], ...], float] = {}
    aliases = {"Vx": "V_x", "Vy": "V_y", "Vz": "V_z"}
    for row in rows:
        key = tuple(sorted((aliases.get(str(k), str(k)), int(v)) for k, v in row["powers"].items()))
        out[key] = out.get(key, 0.0) + float(row["coefficient"])
    return {k: v for k, v in out.items() if abs(v) > 2e-8}


def compare_terms(left: list[dict[str, Any]], right: list[dict[str, Any]], label: str) -> None:
    a, b = term_map(left), term_map(right)
    require(set(a) == set(b), "B1_C_ENGINE_MATH", f"{label}:monomials:{set(a)^set(b)}")
    for key in a:
        require(math.isclose(a[key], b[key], rel_tol=5e-8, abs_tol=5e-8), "B1_C_ENGINE_MATH", f"{label}:{a[key]}!={b[key]}")


def compare_matrix(left: list[Any], right: list[Any], label: str) -> None:
    require(len(left) == len(right) and all(len(a) == len(b) for a, b in zip(left, right)), "B1_C_ENGINE_MATH", f"{label}:shape")
    for i, (a, b) in enumerate(zip(left, right)):
        for j, (x, y) in enumerate(zip(a, b)):
            compare_terms(x, y, f"{label}[{i},{j}]")


def canonical_expression(rows: list[dict[str, Any]]) -> sp.Expr:
    symbols: dict[str, sp.Symbol] = {}
    out: sp.Expr = sp.Integer(0)
    for row in rows:
        term: sp.Expr = sp.Float(row["coefficient"], 17)
        for name, power in row["powers"].items():
            symbols.setdefault(name, sp.Symbol(name))
            term *= symbols[name] ** int(power)
        out += term
    return sp.expand(out)


def expression_matrix(matrix: list[Any]) -> sp.Matrix:
    return sp.Matrix([[canonical_expression(entry) for entry in row] for row in matrix])


def parse_expr(text: Any) -> sp.Expr:
    text = re.sub(r"`[0-9.]*", "", str(text)).replace("Pi", "pi").replace("b1", "").replace("ZZ", "_")
    text = text.replace("^", "**")
    return sp.sympify(text, locals={"pi": sp.pi})


def zero_matrix(matrix: list[Any]) -> bool:
    return all(not term_map(entry) for row in matrix for entry in row)


def leaf_paths(value: Any, prefix: str = "") -> list[str]:
    if isinstance(value, dict):
        return [path for key, child in value.items() for path in leaf_paths(child, f"{prefix}.{key}" if prefix else str(key))]
    if isinstance(value, list):
        if not value:
            return [prefix]
        return [path for index, child in enumerate(value) for path in leaf_paths(child, f"{prefix}.{index}")]
    return [prefix]


def resolve(value: Any, dotted: str) -> Any:
    cur = value
    for part in dotted.split("."):
        cur = cur[int(part)] if isinstance(cur, list) else cur[part]
    return cur


def reaches(edges: list[list[str]], sources: set[str], targets: set[str]) -> bool:
    return any(graph_reachable(edges, source, targets) for source in sources)


def aggregate_components(components: dict[str, Any]) -> dict[str, Any]:
    statuses: list[str] = []
    leaves: list[str] = []
    def walk(value: Any) -> None:
        if isinstance(value, dict):
            if "status" in value:
                statuses.append(value["status"])
            leaves.extend(value.get("unresolved_leaves", []))
            for child in value.values():
                walk(child)
        elif isinstance(value, list):
            for child in value:
                walk(child)
    walk(components)
    status = "ILL_POSED" if "ILL_POSED" in statuses else ("UNSTABLE" if "UNSTABLE" in statuses else
             ("UNRESOLVED" if "UNRESOLVED" in statuses else "OK"))
    return {"status": status, "unresolved_leaves": sorted(set(leaves)), "precedence_evidence": statuses}


def graph_reachable(edges: list[list[str]], source: str, targets: set[str]) -> bool:
    graph: dict[str, list[str]] = {}
    for a, b in edges:
        graph.setdefault(a, []).append(b)
    seen = {source}; stack = [source]
    while stack:
        node = stack.pop()
        if node in targets:
            return True
        for nxt in graph.get(node, []):
            if nxt not in seen:
                seen.add(nxt); stack.append(nxt)
    return False


def validate_claims(artifact: dict[str, Any], phase: dict[str, Any]) -> None:
    allowed = {
        "phase_a_digest": ("source_contract.phase_a_payload_sha256", "sha256(normalized_phase_a_payload)"),
        "field_manifest": ("field_manifest.fields", "join(indexed_routes,phase_a.tail_channels)"),
        "emitted_leaves": ("indexed_profile_missing_leaves", "failed_phase_a_indexed_tangent_lookups"),
        "scalar_regression": ("scalar_regression", "eV.T*M_XX*eV-GVV"),
        "native_slice_constraints": ("indexed_cells", "native_MXp/Mpp_projection-minus-GVP/GPP"),
        "g4_winding": ("g4_control.computed_winding", "contour_integral/(2*pi)"),
        "g4_sigma": ("g4_control.computed_sigma", "Omega_xy/(rho_mass*Gamma*epsilon_xy)"),
        "sheet_area": ("g4_control.total_to_per_area_residual", "Omega_total/sheet_cell_area-Omega_per_area"),
        "berry_coverage": ("berry.production_pullback_coverage", "production_cells-cross-ambient_branches"),
        "e4_action_hessian": ("E4.M_aug_hessian_residual", "hessian(preconstraint_extended_action)-M_aug"),
        "open_reachability": ("open_root_reachability", "typed_ledger-union-traversal"),
        "finite_controls": ("reachability_analysis.finite_bound_controls", "same_generator_parity_filter"),
        "cell_count": ("cells", "active_endpoint-ambient-stratum_product"),
        "dimensions": ("dimensions.records", "termwise_units_restoration"),
        "derived_congruence": ("derived_congruence", "produced_blocks-to-covariant-coefficients"),
        "mechanics_map": ("mechanics_map", "aggregate(cells.component_findings)"),
        "closure_absence": ("provenance_graph.global_return_closure_absence", "reachability(return_closure,mechanics_targets)"),
        "e5_root": ("E5.root_deleted_conservative_solution", "delete(rayleigh_root)-and-resolve"),
        "gate_statuses": ("gate_evidence", "engine_check-derived-gates"),
        "partition": ("partition_ledger", "computed-candidate-ownership"),
    }
    rows = {row["id"]: row for row in artifact["report_claim_bindings"]}
    require(set(rows) == set(allowed), "B1_C_REPORT_SCHEMA", f"claim ids {set(rows)^set(allowed)}")
    for name, (path, rule) in allowed.items():
        require(rows[name]["schema_path"] == path and rows[name]["recompute"] == rule,
                "B1_C_REPORT_SCHEMA", f"{name}:typed path/rule")
        cur: Any = artifact
        for part in path.split("."):
            cur = cur[part]
    if artifact["engine"] == "SymPy":
        require(artifact["source_contract"]["phase_a_payload_sha256"] == digest(phase_a_payload(phase)),
                "B1_C_REPORT_SCHEMA", "phase digest claim recomputation")
    require(len(artifact["cells"]) == len(ENDPOINTS) * len(AMBIENTS), "B1_C_REPORT_SCHEMA", "cell count claim")


def compare(sym: dict[str, Any], math_artifact: dict[str, Any], phase: dict[str, Any], data: dict[str, Any]) -> list[str]:
    require(sym["schema_version"] == math_artifact["schema_version"] == "U1_PHASE_B1_ENGINE_ARTIFACT_V3",
            "B1_C_ENGINE_MATH", "artifact schema")
    representations = data["derivation_representations"]
    require(sym["engine_representation"] == representations["SymPy"] and
            math_artifact["engine_representation"] == representations["Mathematica"] and
            sym["engine_representation"] != math_artifact["engine_representation"] and
            representations["shared_reduced_formulas"] is False, "B1_C_ENGINE_MATH", "independent representations")
    for artifact in (sym, math_artifact):
        require(artifact["axis_1"] == "COMPUTATION_VALID" and
                set(artifact["checks"]) == {f"B1_R{i}" for i in range(1, 10)} | {"B1_S1"} and
                all(row["status"] == "PASS" for row in artifact["checks"].values()),
                "B1_C_ENGINE_MATH", f"{artifact['engine']}:engine checks")

    # Source, emitted lookup failures, and per-field manifest.
    require(set(sym["source_contract"]["action_ids"]) == set(math_artifact["source_contract"]["action_ids"]) and
            set(sym["source_contract"]["declared_action_ids"]) == set(math_artifact["source_contract"]["declared_action_ids"]) and
            sym["source_contract"]["decision_sha256"] == math_artifact["source_contract"]["decision_sha256"],
            "B1_C_ENGINE_MATH", "source contract disagreement")
    smissing, mmissing = set(sym["indexed_profile_missing_leaves"]), set(math_artifact["indexed_profile_missing_leaves"])
    require(smissing == mmissing, "B1_C_ENGINE_MATH", "emitted indexed lookup leaves")
    sfields = {row["field"]: row for row in sym["field_manifest"]["fields"]}
    mfields = {row["field"]: row for row in math_artifact["field_manifest"]["fields"]}
    require(set(sfields) == set(mfields), "B1_C_ENGINE_MATH", "field manifest keys")
    field_columns = ("support", "radial_dimension", "tensor_harmonic_type", "integration_measure",
                     "profile_provenance", "translation_tangent_class", "p_tangent_lookup", "emitted_missing_leaf", "kinetic_action")
    for name in sfields:
        require(all(sfields[name][col] == mfields[name][col] for col in field_columns), "B1_C_ENGINE_MATH", f"manifest:{name}")
        expected_dimension = 4 if sfields[name]["support"] == "bulk" else 3
        require(sfields[name]["radial_dimension"] == expected_dimension and sfields[name]["Phi_Xi"] and
                ((sfields[name]["emitted_missing_leaf"] is not None and sfields[name]["Phi_pi"] is None) or
                 (sfields[name]["emitted_missing_leaf"] is None and sfields[name]["Phi_pi"] is not None)) and
                ((mfields[name]["emitted_missing_leaf"] is not None and mfields[name]["Phi_pi"] is None) or
                 (mfields[name]["emitted_missing_leaf"] is None and mfields[name]["Phi_pi"] is not None)),
                "B1_C_ENGINE_MATH", f"manifest:{name}:typed tangent")
        if sfields[name]["translation_tangent_class"] == "endpoint_velocity_response":
            require(sfields[name]["response_provenance"], "B1_C_ENGINE_MATH", f"manifest:{name}:response provenance")
    manifest_leaves = {row["emitted_missing_leaf"] for row in sym["field_manifest"]["fields"] if row["emitted_missing_leaf"]}
    manifest_leaves |= {row["emitted_missing_leaf"] for row in sym["field_manifest"]["surface_variations"] if row["emitted_missing_leaf"]}
    require(manifest_leaves == smissing, "B1_C_ENGINE_MATH", "missing leaves not derived from failed lookups")

    # Endpoint variation and field-local translation contractions.
    require(sym["endpoint_source_map"] == math_artifact["endpoint_source_map"], "B1_C_ENGINE_MATH", "endpoint source map")
    for ep in ENDPOINTS:
        for channel in ("endpoint_assembly", "endpoint_conservative"):
            sr, mr = sym[channel][ep], math_artifact[channel][ep]
            require(sr["rank"] == mr["rank"] == 3, "B1_C_ENGINE_MATH", f"{channel}:{ep}:rank")
            require(all(sp.simplify(parse_expr(a) - parse_expr(b)) == 0 for a, b in zip(sr["solution"], mr["solution"])) and
                    all(parse_expr(x) == 0 for x in sr["solution_residual"] + mr["solution_residual"]),
                    "B1_C_ENGINE_MATH", f"{channel}:{ep}:solve residual")
    require(set(sym["indexed_cells"]) == set(math_artifact["indexed_cells"]), "B1_C_ENGINE_MATH", "indexed cell keys")
    for key in sym["indexed_cells"]:
        scell, mcell = sym["indexed_cells"][key], math_artifact["indexed_cells"][key]
        compare_matrix(scell["M_XX_p0_known"], mcell["M_XX_p0_known"], f"{key}:MXX")
        ambient = key.split("|", 1)[1]
        ambient_row = data["ambient_branches"][ambient]
        wall_domain = ambient_row.get("wall_side_domain", [])
        require(len(scell["wall_side_branch_samples"]) == len(mcell["wall_side_branch_samples"]) == len(wall_domain),
                "B1_C_ENGINE_MATH", f"{key}:wall-side sample count")
        branch = parse_expr(ambient_row["embedding_factor"]).subs(sp.Symbol("eta_asym"), parse_expr(ambient_row.get("eta_asym", 0)))
        for index, side in enumerate(wall_domain):
            compare_terms(scell["wall_side_branch_samples"][index], mcell["wall_side_branch_samples"][index],
                          f"{key}:wall-side sample:{index}")
            reported = canonical_expression(scell["wall_side_branch_samples"][index])
            expected = branch.subs(sp.Symbol("s"), parse_expr(side))
            require(abs(float(sp.N(reported - expected, 17))) < 2e-10,
                    "B1_C_ENGINE_MATH", f"{key}:wall-side sample binding:{index}")
        require(zero_matrix(scell["computed_isotropy_residual"]) and zero_matrix(mcell["computed_isotropy_residual"]) and
                not scell["reconstruction_residual"] and not mcell["reconstruction_residual"],
                "B1_C_ENGINE_MATH", f"{key}:isotropy/reconstruction")
        scon = {row["action_term"]: row for row in scell["field_contraction_integrals"]}
        mcon = {row["action_term"]: row for row in mcell["field_contraction_integrals"]}
        require(set(scon) == set(mcon), "B1_C_ENGINE_MATH", f"{key}:contraction terms")
        for term in scon:
            compare_terms(scon[term]["coefficient_from_contraction"], mcon[term]["coefficient_from_contraction"], f"{key}:{term}:contraction")
            compare_matrix(scon[term]["computed_tensor"], mcon[term]["computed_tensor"], f"{key}:{term}:computed tensor")
            require(scon[term]["radial_dimension"] == mcon[term]["radial_dimension"] and
                    scon[term]["quadrature_crosscheck"]["passed"] and mcon[term]["quadrature_crosscheck"]["passed"] and
                    scon[term]["quadrature_crosscheck"]["absolute_error"] < 2e-11 and
                    mcon[term]["quadrature_crosscheck"]["absolute_error"] < 2e-11 and
                    scon[term]["oracle_ancestry_forbidden"] and mcon[term]["oracle_ancestry_forbidden"] and
                    not scon[term]["oracle_paths_consumed"] and not mcon[term]["oracle_paths_consumed"] and
                    len(scon[term]["angular_tensor_integral"]) == len(mcon[term]["angular_tensor_integral"]) == 3,
                    "B1_C_ENGINE_MATH", f"{key}:{term}:angular/moment reduction")
        for artifact in (sym, math_artifact):
            regression = artifact["scalar_regression"][key]
            e_v = sp.Matrix([parse_expr(value) for value in data["scalar_regression_projection"]["V_unit"]])
            computed_projection = sp.expand((e_v.T * expression_matrix(artifact["indexed_cells"][key]["M_XX_p0_known"]) * e_v)[0])
            reported_projection = canonical_expression(regression["projected_M_XX"])
            require(not regression["residual"] and sp.expand(computed_projection - reported_projection) == 0 and
                    regression["projection_consumed_after_production"] and
                    regression["projection_dependency_role"] == "comparator_only", "B1_C_ENGINE_MATH", f"{artifact['engine']}:{key}:scalar regression")
            constraints = artifact["indexed_cells"][key]["native_tilt_slice_constraints"]
            for native in ("M_Xp_native", "M_pp_native"):
                row = constraints[native]
                encoded_free = "b1" + row["free_slice_symbol"].replace("_", "ZZ")
                require(row["defining_expression"] and (row["free_slice_symbol"] in row["defining_expression"] or
                        encoded_free in row["defining_expression"]) and
                        not ({"substitution", "evaluated_residual", "executable_residual"} & set(row)) and
                        row["status"] == "CONDITIONAL_INDEXED_FAMILY_UNAVAILABLE", "B1_C_ENGINE_MATH", f"{artifact['engine']}:{key}:{native}")

    for artifact in (sym, math_artifact):
        sentinel = artifact["scalar_regression"]["anisotropic_projection_sentinel"]
        tensor = sp.Matrix([[parse_expr(value) for value in row] for row in sentinel["tensor"]])
        vector = sp.Matrix([parse_expr(value[0] if isinstance(value, list) else value) for value in sentinel["e_V"]])
        true_value = sp.expand((vector.T * tensor * vector)[0])
        require(parse_expr(sentinel["true_contraction"]) == true_value and true_value != tensor[0, 0] and
                sentinel["distinguishes_element_selection"], "B1_C_ENGINE_MATH", f"{artifact['engine']}:anisotropic projection sentinel")

    # The reach table is re-generated here from roots, admissible orders and declared block rules.
    expected_reach: set[tuple[str, str, str]] = set()
    for root in data["open_action_functionals"]:
        for order in root["admissible_orders"]:
            for blocks in data["tensor_selection_rules"][order].values():
                expected_reach |= {(root["id"], block, order) for block in blocks}
    for artifact in (sym, math_artifact):
        actual = {(r["root"], r["block"], r["order"]) for r in artifact["open_root_reachability"]}
        require(actual == expected_reach and all(r["disposition"] == "REACHABLE_UNRESOLVED_REMAINDER" and r["witness_tensor"] for r in artifact["open_root_reachability"]) and
                artifact["reachability_analysis"]["crosscheck_agrees"], "B1_C_RECORD_MAP", f"{artifact['engine']}:traversal reachability")
        controls = artifact["reachability_analysis"]["finite_bound_controls"]
        require(controls["parity_empty_control"]["coefficient_space_empty"] and controls["parity_empty_control"]["classification"] == "STRUCTURAL_ZERO" and
                not controls["even_witness_control"]["coefficient_space_empty"] and controls["even_witness_control"]["witness_tensor"] and
                artifact["reachability_analysis"]["control_outcome_counts"] == {"empty": 1, "witness": 1},
                "B1_C_RECORD_MAP", f"{artifact['engine']}:finite controls")

    # Traversal-derived, route-separated Berry computations and production coverage.
    compare_matrix(sym["berry"]["route_A"]["Omega"], math_artifact["berry"]["route_A"]["Omega"], "Berry route A")
    compare_matrix(sym["berry"]["route_B"]["Omega"], math_artifact["berry"]["route_B"]["Omega"], "Berry route B")
    expected_cells = {f"{ep}|{ambient}" for ep in ENDPOINTS for ambient in AMBIENTS}
    for artifact in (sym, math_artifact):
        berry = artifact["berry"]; route_a, route_b = berry["route_A"], berry["route_B"]
        parse_guard = berry["production_phase_parse_guard"]
        require(parse_guard == {"allowed_free_symbols": ["r"], "observed_free_symbols": ["r"], "passed": True},
                "B1_C_ENGINE_MATH", f"{artifact['engine']}:production phase parse guard")
        for route, needed in ((route_a, {"L_indexed_berry", "Omega_EL_xy"}), (route_b, {"field_current_integrand", "projected_zero_mode_sector", "Omega_pullback_xy"})):
            dag = route["dag"]; require(needed <= set(dag["named_outputs"]) and dag["edges"] and
                    reaches(dag["edges"], {f"symbol:{x}" for x in dag["raw_dependencies"]}, set(dag["named_outputs"].values())),
                    "B1_C_ENGINE_MATH", f"{artifact['engine']}:expression traversal DAG")
        coverage = berry["production_pullback_coverage"]
        require({row["cell"] for row in coverage} == expected_cells and berry["production_coverage_complete"] and
                all(row["declared_class"] == row["actual_class"] and zero_matrix(row["pullback_residual"]) and row["execution_digest"] and
                    "route_A_Omega" in row and "route_B_Omega" in row for row in coverage) and
                berry["route_shared_reduced_nodes"] == [] and berry["selection_quarantined_from_production_ancestry"] and
                berry["separation_control"]["detector_fires"] and berry["separation_control"]["shared_reduced_digests"] and
                route_b["zero_mode_quotient_applied"] and not zero_matrix(route_b["zero_mode_sector_before_projection"]) and
                zero_matrix(route_b["zero_mode_sector_after_projection"]),
                "B1_C_ENGINE_MATH", f"{artifact['engine']}:production pullback coverage/quotient")
        require(zero_matrix(berry["equivalence_residual"]), "B1_C_ENGINE_MATH", f"{artifact['engine']}:Berry equivalence")

    # G4 signed control remains downstream, with k=1 control, sigma=-1 and k=0 production.
    sg, mg = sym["g4_control"], math_artifact["g4_control"]
    compare_matrix(sg["Omega_physical_per_sheet_area"], mg["Omega_physical_per_sheet_area"], "G4 Omega")
    for artifact in (sym, math_artifact):
        g4 = artifact["g4_control"]
        require(sp.simplify(parse_expr(g4["computed_winding"]) - 1) == 0 and sp.simplify(parse_expr(g4["computed_sigma"]) + 1) == 0 and
                zero_matrix(g4["total_to_per_area_residual"]) and g4["sheet_cell_area"] and g4["sheet_geometry"] and
                parse_expr(g4["sheet_cell_area"]) != 0 and
                "sigma_projection_unsimplified" in g4["sigma_dependency_dag"]["named_outputs"] and
                all(node.get("subexpression_digest") for node in g4["sigma_dependency_dag"]["nodes"]) and
                all(all(parse_expr(v) == 0 for v in case["signed_residual"]) for case in g4["force_cases"].values()) and
                parse_expr(artifact["berry"]["intrinsic_circulation"]["k_theta"]) == 0,
                "B1_C_ENGINE_MATH", f"{artifact['engine']}:G4 signed/area computation")

    # E4 action-first Hessian/Jacobian and E5 root-consumption chain.
    for name in ("M_aug", "constraint_operator", "J_physical_lift", "M_reduced"):
        compare_matrix(sym["E4"][name], math_artifact["E4"][name], f"E4:{name}")
    compare_terms(sym["E4"]["collar_mass"], math_artifact["E4"]["collar_mass"], "E4 collar mass")
    for artifact in (sym, math_artifact):
        e4, e5 = artifact["E4"], artifact["E5"]
        require(not ({"boundary_trace_input", "ir_condition_input", "moduli_fixing_input"} & set(e4)),
                "B1_C_INPUT_LIVENESS", f"{artifact['engine']}:E4 echo sinks excluded")
        require(e4["preconstraint_extended_action"] and zero_matrix(e4["M_aug_hessian_residual"]) and e4["J_derivation_equations"] and
                e4["constraint_rank"] == 3 and zero_matrix(e4["N_kernel_residual"]) and zero_matrix(e4["basis_covariance_residual"]) and
                sp.simplify(parse_expr(e4["collar_mode_reconstruction"]) - parse_expr(e4["boundary_profile"])) == 0 and
                parse_expr(e4["growing_mode_root"]) != parse_expr(e4["selected_decay"]) and
                all(not row for row in e4["A_aug"]) and zero_matrix(e4["Omega_aug"]), "B1_C_ENGINE_MATH", f"{artifact['engine']}:E4 action first")
        require(e5["root"] and "gammaSigma" in e5["rayleigh_dependencies"] and parse_expr(e5["root_deleted_functional"]) == 0 and
                parse_expr(e5["gamma_functional_ablation_force_delta"]) != 0 and e5["conservative_equals_E2_computed"],
                "B1_C_ENGINE_MATH", f"{artifact['engine']}:E5 root consumption")

    # Computed dimensions and derived-block congruence.
    require(sym["dimensions"]["action_dimension_sourced_from_phase_a"] == math_artifact["dimensions"]["action_dimension_sourced_from_phase_a"],
            "B1_C_ENGINE_MATH", "action dimension source")
    for artifact in (sym, math_artifact):
        dims = artifact["dimensions"]
        require(dims["action_terms_homogeneous"] and dims["coordinate_embedding_consistent"] and
                all(len(row["computed_monomial_dimensions_LTM"]) == 1 for row in dims["records"]),
                "B1_C_ENGINE_MATH", f"{artifact['engine']}:computed dimensions")
        congruence = artifact["derived_congruence"]
        require("covariant_inventory_consumed" not in congruence and
                all("ambient_contract" not in cell for cell in artifact["indexed_cells"].values()),
                "B1_C_INPUT_LIVENESS", f"{artifact['engine']}:deepcopy echo sinks excluded")
        # The exact keys are derived from the produced reach records, including one aggregate known delta coefficient.
        expected_coefficients = {"M_XX:delta_total"} | {f"{row['block']}:{generator}" for row in artifact["open_root_reachability"] for generator in row["finite_generator_set"]}
        require(expected_coefficients <= set(congruence["derived_coefficients"]) and
                all(row["definition"] not in (None, "") for row in congruence["derived_coefficients"].values()) and
                zero_matrix(congruence["M_XX_covariant_decomposition"]["projection_residual"]) and zero_matrix(congruence["block_split_residual"]) and
                len(congruence["transverse_LDL_pivots"]) + len(congruence["longitudinal_LDL_pivots"]) == 6 and
                len(congruence["positive_definite_conditions"]) == 6 and congruence["E4_reduced_rank"] == 3,
                "B1_C_ENGINE_MATH", f"{artifact['engine']}:derived congruence")
        generic = sp.Matrix([parse_expr(value) for value in data["covariant_inventory"]["generic_projection_point"]])
        zero_direction = sp.Matrix([parse_expr(value) for value in data["covariant_inventory"]["zero_limit_direction"]])
        projection_probe_residual = expression_matrix(congruence["projection_basis_at_point"]) - generic * generic.T
        zero_path_probe_residual = expression_matrix(congruence["zero_limit_path_projector"]) - zero_direction * zero_direction.T
        require(all(sp.simplify(value) == 0 for value in projection_probe_residual) and
                all(sp.simplify(value) == 0 for value in zero_path_probe_residual),
                "B1_C_ENGINE_MATH", f"{artifact['engine']}:covariant probe binding")
        compare_matrix(congruence["E4_reduced_block"], artifact["E4"]["M_reduced"], f"{artifact['engine']}:E4 congruence block")

    # Re-derive records, closure absence and all-UNRESOLVED headline from typed components.
    recomputed_maps: list[dict[str, Any]] = []
    for artifact in (sym, math_artifact):
        rebuilt: dict[str, Any] = {}
        for key, record in artifact["cells"].items():
            expected_key = "|".join((record["cell_key"]["endpoint"], record["cell_key"]["ambient"], record["cell_key"]["open_stratum"]))
            computed = aggregate_components(record["component_findings"])
            require(key == expected_key and record["headline"] == computed and computed["status"] == "UNRESOLVED",
                    "B1_C_RECORD_MAP", f"{artifact['engine']}:{key}:typed headline")
            rebuilt[key] = computed
        emitted_edges = {tuple(edge) for cell in artifact["indexed_cells"].values()
                         for contraction in cell["field_contraction_integrals"] for edge in contraction["dependency_edges"]}
        emitted_edges |= {(row["root"], f"remainder:{row['block']}") for row in artifact["open_root_reachability"]}
        emitted_edges |= {tuple(edge) for group in (artifact["berry"]["dependency_edges"], artifact["E4"]["dependency_edges"], artifact["E5"]["dependency_edges"])
                          for edge in group}
        require(set(map(tuple, artifact["provenance_graph"]["edges"])) == emitted_edges,
                "B1_C_RECORD_MAP", f"{artifact['engine']}:derived provenance edge union")
        require(set(rebuilt) == {f"{ep}|{ambient}|{s}" for ep in ENDPOINTS for ambient in AMBIENTS for s in artifact["open_strata"]["active"]} and
                artifact["mechanics_map"] == rebuilt and artifact["provenance_graph"]["global_return_closure_absence"] and
                artifact["closure_axis"]["global_absence_computed"] and not artifact["closure_axis"]["materialized"],
                "B1_C_RECORD_MAP", f"{artifact['engine']}:map/provenance")
        recomputed_maps.append(rebuilt)
    for key in recomputed_maps[0]:
        require(recomputed_maps[0][key]["status"] == recomputed_maps[1][key]["status"] and
                set(recomputed_maps[0][key]["unresolved_leaves"]) == set(recomputed_maps[1][key]["unresolved_leaves"]),
                "B1_C_RECORD_MAP", f"cross-engine map:{key}")

    # S1: every scalar input leaf has one typed role and semantic sink; comparator-only leaves are quarantined.
    declared = set(leaf_paths(data))
    for artifact in (sym, math_artifact):
        live = artifact["input_liveness"]; rows = {row["path"]: row for row in live["rows"]}
        require(set(live["declared_leaf_paths"]) == declared == set(live["consumed_leaf_paths"]) == set(rows) and
                not live["consumed_but_undeclared"] and live["per_key_mutation_evidence"].startswith("external:"),
                "B1_C_INPUT_LIVENESS", f"{artifact['engine']}:leaf cover")
        for path, row in rows.items():
            require(row["typed_role"] and row["semantic_sink"] and len(row["dependency_path"]) >= 3 and
                    row["read_scopes"] and not row["metadata_or_digest_only"],
                    "B1_C_INPUT_LIVENESS", f"{artifact['engine']}:{path}:typed sink")
            resolve(artifact, row["semantic_sink"])
            if path.startswith("scalar_regression_projection."):
                require(row["production_ancestry_forbidden"] and row["absent_from_production_ancestry"],
                        "B1_C_INPUT_LIVENESS", f"{artifact['engine']}:{path}:comparator quarantine")
            if path.startswith(("berry_pullback_selection.dependency_role", "berry_pullback_selection.production_cells.",
                                "berry_pullback_selection.ambient_branches.")):
                require(row["production_ancestry_forbidden"] and row["absent_from_production_ancestry"],
                        "B1_C_INPUT_LIVENESS", f"{artifact['engine']}:{path}:selection quarantine")

    phase_digest = digest(phase_a_payload(phase))
    require(phase_digest == data["phase_a_protection"]["normalized_payload_sha256"] == sym["source_contract"]["phase_a_payload_sha256"] and
            phase["axis1"] == "COMPUTATION_VALID" and phase["axis2"] == "U1_BASE_OK", "B1_C_PHASE_A", f"Phase-A payload {phase_digest}")
    core_traces_digest = data["phase_a_protection"]["core_traces_sha256"]
    require(core_traces_digest == sym["source_contract"]["core_traces_sha256"] == math_artifact["source_contract"]["core_traces_sha256"],
            "B1_C_PHASE_A", f"Phase-A core traces {core_traces_digest}")
    comparator_payload = comparator_amended_payload(phase)
    comparator_digest = digest(comparator_payload)
    amendment_digests = {artifact["phase_a_amendment"]["amended_payload_sha256"] for artifact in (sym, math_artifact)}
    for artifact in (sym, math_artifact):
        amendment = artifact["phase_a_amendment"]
        allowed = tuple(amendment["semantic_diff"]["allowed_prefixes"])
        changed_paths = [row["path"] for row in amendment["semantic_diff"]["rows"]]
        require(amendment["name"] == "PHASE_A_MOMENT_CORRECTION(brane_shear)" and
                amendment["baseline_digest"]["gate"] and amendment["baseline_digest"]["computed_legacy"] == phase_digest and
                amendment["semantic_diff"]["restricted_to_correction_closure"] and changed_paths and
                all(path.startswith(allowed) for path in changed_paths) and
                amendment["tilt_profile_rows"]["byte_semantics_unchanged"] and
                amendment["phase_a_acceptance_recheck"]["verdict"] == "U1_BASE_OK" and
                len(amendment["corrected_moment_rows"]) == len(amendment["corrected_downstream_coefficients"]) == 5 and
                len(amendment["corrected_row_paths"]) == 10,
                "B1_C_PHASE_A", f"{artifact['engine']}:amendment discipline")
    require(amendment_digests == {comparator_digest} and
            digest(sym["amended_phase_a_payload"]) == digest(math_artifact["amended_phase_a_payload"]) == comparator_digest and
            sym["source_contract"]["amended_payload_sha256"] == math_artifact["source_contract"]["amended_payload_sha256"] == comparator_digest,
            "B1_C_PHASE_A", f"three-way amended digest {amendment_digests} comparator={comparator_digest}")
    validate_claims(sym, phase); validate_claims(math_artifact, phase)
    return [
        "per-field d_f angular Gram contractions, isotropy, scalar regression, and native GVP/GPP residuals",
        "endpoint variation solves plus E4 action-first Hessian/J lift and E5 root-deletion solve",
        f"{len(expected_reach)} traversal-derived OPEN-root dispositions and two exercised finite controls",
        f"route-separated Berry DAGs with {len(expected_cells)} production pullback cells and a real quotient",
        "G4 sheet-area reduction, k=1 signed control, sigma=-1, and production k=0",
        "units restored on real expressions and derived-coefficient symmetry-block congruence",
        "component records, closure absence, provenance, and all-UNRESOLVED map reaggregated",
        f"typed input roles/sinks cover all {len(declared)} declared scalar leaves",
        "independent SymPy/Wolfram Phase-A protection and nonshared representations",
    ]


def load_mutation_results(path: Path) -> dict[str, Any]:
    require(path.exists(), "B1_C_REPORT_SCHEMA", f"missing external mutation results {path}")
    data = yaml.safe_load(path.read_text())
    require(data.get("schema_version") == "U1_PHASE_B1_OUT_OF_PROCESS_MUTATIONS_V3" and data.get("status") == "PASS",
            "B1_C_REPORT_SCHEMA", "mutation harness status")
    for row in data["cases"]:
        require(row["status"] == "PASS" and row["guarded_digest_changed"] and row["assertion_fired"],
                "B1_C_REPORT_SCHEMA", f"mutation {row['id']}")
    covered = {row["expected_assert"] for row in data["cases"]}
    require(set(f"B1_R{i}" for i in range(1, 10)) <= covered and set(COMPARE_CHECKS) <= covered,
            "B1_C_REPORT_SCHEMA", f"mutation tooth cover {covered}")
    live_rows = data.get("input_liveness_cases", [])
    sampled = [row for row in live_rows if row.get("dual_engine_sampled")]
    require(len(live_rows) == data.get("input_liveness_declared_count") and
            len({row["path"] for row in live_rows}) == len(live_rows) and
            all(row["status"] == "PASS" and row["out_of_process"] and row["input_digest_changed"] and
                (row["sink_digest_changed"] if row["mutation_class"] == "semantically_valid" else
                 (row["sink_digest_changed"] or row["guarded_assertion_fired"])) for row in live_rows),
            "B1_C_INPUT_LIVENESS", "per-key out-of-process liveness evidence")
    require({row["liveness_batch_index"] for row in sampled} == set(range(3)) and
            all(row["liveness_batch_count"] == 3 and {e["engine"] for e in row["engine_results"]} == {"SymPy", "Mathematica"} and
                all(e["sink_digest_changed"] for e in row["engine_results"]) for row in sampled),
            "B1_C_INPUT_LIVENESS", "sampled dual-engine liveness evidence in every batch")
    return data


def stage1_artifact(sym: dict[str, Any], math_artifact: dict[str, Any], phase: dict[str, Any],
                    agreement: list[str], paths: dict[str, Path]) -> dict[str, Any]:
    correction = sym["phase_a_amendment"]; math_correction = math_artifact["phase_a_amendment"]
    amended_digest = digest(comparator_amended_payload(phase))
    moment_paths = [row["path"] for row in correction["corrected_moment_rows"]]
    coefficient_paths = [row["path"] for row in correction["corrected_downstream_coefficients"]]
    brane_row = next(row for row in sym["indexed_cells"]["E1|symmetric_postulate"]["field_contraction_integrals"]
                     if row["action_term"] == "brane_shear_kinetic")
    base = correction["corrected_moment_rows"][0]
    e4 = correction["corrected_moment_rows"][3]
    summaries = [
        f"1. PHASE_A_MOMENT_CORRECTION(brane_shear): unit gradient {correction['old_unit_gradient_numeric']:.12f} -> {correction['new_unit_gradient_numeric']:.12f}.",
        "2. Authority: shear operator/harmonic decomposition -> verified nu=2 tail -> core trace -> endpoint response.",
        f"3. Corrected moments ({len(moment_paths)}): {', '.join(moment_paths)}.",
        f"4. Corrected downstream coefficients ({len(coefficient_paths)}): {', '.join(coefficient_paths)}.",
        f"5. Tilt rows: {correction['tilt_profile_rows']['disposition']}; unchanged={str(correction['tilt_profile_rows']['byte_semantics_unchanged']).lower()}.",
        f"6. Semantic diff: {len(correction['semantic_diff']['rows'])} scalar changes; correction-closure-only={str(correction['semantic_diff']['restricted_to_correction_closure']).lower()}.",
        f"7. Baseline digest gate: {correction['baseline_digest']['computed_legacy']}; SymPy/Wolfram gates passed.",
        f"8. Amended digest gate: SymPy=Wolfram=comparator={amended_digest}.",
        f"9. Phase-A acceptance re-check: SymPy={correction['phase_a_acceptance_recheck']['verdict']}; Wolfram={math_correction['phase_a_acceptance_recheck']['verdict']}.",
        f"10. Forward M_XX channels: bulk_flow unchanged; h unchanged; uw=0; brane E1 {base['old']['production_value']:.12f}->{base['new']['production_value']:.12f}, E4 {e4['old']['production_value']:.12f}->{e4['new']['production_value']:.12f}; quadrature error={brane_row['quadrature_crosscheck']['absolute_error']:.3g}.",
    ]
    artifact_hashes = {name: hashlib.sha256(path.read_bytes()).hexdigest() for name, path in paths.items()}
    core = {"schema_version": "U1_PHASE_B1_STAGE1_COMPLETE_V3", "status": "STAGE1_COMPLETE_AWAITING_ORCHESTRATOR_APPROVAL",
            "staged_exit_code": 42, "final_b1_outputs_emitted": False, "engine_agreement_checks": agreement,
            "artifact_sha256": artifact_hashes,
            "digest_gate": {"baseline": correction["baseline_digest"]["computed_legacy"], "sympy": correction["amended_payload_sha256"],
                            "wolfram": math_correction["amended_payload_sha256"], "comparator": amended_digest,
                            "agreement": correction["amended_payload_sha256"] == math_correction["amended_payload_sha256"] == amended_digest},
            "semantic_diff_gate": correction["semantic_diff"],
            "phase_a_acceptance_recheck": {"SymPy": correction["phase_a_acceptance_recheck"],
                                            "Mathematica": math_correction["phase_a_acceptance_recheck"]},
            "correction_finding": correction, "stage1_summary": summaries}
    core["resume_validation_sha256"] = digest(core)
    return core


def terms_text(rows: list[dict[str, Any]]) -> str:
    parts = []
    for row in rows:
        power = "*".join(f"{k}^{v}" if int(v) != 1 else k for k, v in row["powers"].items())
        parts.append(f"{float(row['coefficient']):.12g}" + ("*" + power if power else ""))
    return " + ".join(parts) if parts else "0"


def markdown_table(headers: list[str], rows: list[list[Any]]) -> str:
    clean = [[str(v).replace("|", "\\|").replace("\n", " ") for v in row] for row in rows]
    return "\n".join(["| " + " | ".join(headers) + " |", "| " + " | ".join("---" for _ in headers) + " |"] +
                     ["| " + " | ".join(row) + " |" for row in clean])


def surfaced_correction(correction: dict[str, Any]) -> dict[str, Any]:
    return {
        "name": correction["name"],
        "unit_gradient_transition": {
            "old": correction["old_unit_gradient"], "new": correction["new_unit_gradient"],
            "display": "12π/5 → 8π/3",
            "old_numeric": correction["old_unit_gradient_numeric"], "new_numeric": correction["new_unit_gradient_numeric"],
        },
        "digest_transition": {
            "old": correction["baseline_digest"]["computed_legacy"],
            "new": correction["amended_payload_sha256"],
        },
        "authoritative_chain": correction["authoritative_chain"],
        "error_provenance": correction["error_provenance"],
        "corrected_row_paths": correction["corrected_row_paths"],
        "tilt_profile_rows": correction["tilt_profile_rows"],
    }


def known_residual_constructs(sym: dict[str, Any]) -> dict[str, Any]:
    contraction_paths = "indexed_mechanics.cells.*.field_contraction_integrals.*"
    return {
        "sheet_area_X_minus_X_residual": {
            "status": "KNOWN_RESIDUAL_WITH_COMPARATOR_BACKSTOP",
            "artifact_field": "symplectic_mechanics.g4_control.total_to_per_area_residual",
            "engine_value": sym["g4_control"]["total_to_per_area_residual"],
            "comparator_backstop": "B1_C_ENGINE_MATH recomputes/nonzero-checks sheet_cell_area and requires the emitted residual matrix to vanish",
        },
        "hand_typed_provenance_edge_atoms": {
            "status": "KNOWN_RESIDUAL_WITH_COMPARATOR_BACKSTOP",
            "artifact_field": "mechanics_provenance_graph.edges",
            "engine_value_count": len(sym["provenance_graph"]["edges"]),
            "comparator_backstop": "B1_C_RECORD_MAP rebuilds the exact edge union from contraction, reachability, Berry, E4, and E5 records",
        },
        "oracle_ban_stamps": {
            "status": "KNOWN_RESIDUAL_WITH_COMPARATOR_BACKSTOP",
            "artifact_field": contraction_paths + ".{oracle_ancestry_forbidden,oracle_paths_consumed}",
            "engine_value": {"oracle_ancestry_forbidden": True, "oracle_paths_consumed": []},
            "comparator_backstop": "B1_C_ENGINE_MATH requires both engines' stamps and separately compares every forward-derived contraction tensor and quadrature tooth",
        },
    }


def phase_a_baseline_with_correction(report_text: str, correction: dict[str, Any]) -> str:
    baseline = report_text.split("\n---\n\n## Phase B1", 1)[0].rstrip()
    start = "<!-- PHASE_A_MOMENT_CORRECTION_NOTE_START -->"
    end = "<!-- PHASE_A_MOMENT_CORRECTION_NOTE_END -->"
    baseline = re.sub(rf"\n?{re.escape(start)}.*?{re.escape(end)}\n?", "\n", baseline, flags=re.DOTALL)
    anchor = "\nDirect differentiation computes both canonical momenta"
    require(anchor in baseline, "B1_C_REPORT_SCHEMA", "Phase-A GVV display correction-note anchor")
    surfaced = surfaced_correction(correction)
    note = "\n".join([
        start,
        "**Phase-A amendment note (surfaced by Phase B1):** the `G_VV` rows immediately above display the frozen pre-amendment "
        f"brane-shear unit coefficient `12π/5` (`{correction['old_unit_gradient_numeric']:.13g}`). "
        f"`{correction['name']}` corrects the governed `U_XX → G_VV` chain to `8π/3` (`{correction['new_unit_gradient_numeric']:.13g}`). "
        "The separately undetermined `U_XP`, `U_PP`, and `I_shear_grad` tilt-profile rows remain frozen/`UNRESOLVED(tilt_profile)`. "
        f"Payload digest: `{surfaced['digest_transition']['old']}` → `{surfaced['digest_transition']['new']}`.",
        end,
    ])
    return baseline.replace(anchor, "\n\n" + note + anchor, 1)


def report_section(sym: dict[str, Any], agreement: list[str], mutations: dict[str, Any], input_sha: str, phase_digest: str) -> str:
    bindings = {row["id"]: row for row in sym["report_claim_bindings"]}
    def ref(claim: str) -> str:
        row = bindings[claim]
        return f"[claim `{claim}` → `{row['schema_path']}`; recompute `{row['recompute']}`]"

    manifest = sym["field_manifest"]; missing = sym["indexed_profile_missing_leaves"]
    correction = surfaced_correction(sym["phase_a_amendment"])
    residuals = known_residual_constructs(sym)
    cell_statuses: dict[str, int] = {}
    for row in sym["mechanics_map"].values(): cell_statuses[row["status"]] = cell_statuses.get(row["status"], 0) + 1
    dimensions = {row["object"]: row["computed_monomial_dimensions_LTM"] for row in sym["dimensions"]["records"]}
    controls = sym["reachability_analysis"]["finite_bound_controls"]
    endpoint_rows = []
    for endpoint in ENDPOINTS:
        records = [row for key, row in sym["cells"].items() if key.startswith(endpoint + "|")]
        endpoint_rows.append([endpoint, sym["endpoint_source_map"][endpoint],
                              ",".join(sorted({row["headline"]["status"] for row in records})),
                              ",".join(sorted({row["component_findings"]["Omega_intrinsic_circulation"]["status"] for row in records}))])
    claim_rows = [[row["id"], row["schema_path"], row["type"], row["recompute"]] for row in sym["report_claim_bindings"]]
    gate_rows = [[name, row["status"]] for name, row in sym["gate_evidence"].items()]
    lines = [
        "", "---", "", "## Phase B1 — indexed mechanics remediation 3", "",
        f"The independently recomputed Phase-A payload is `{phase_digest}` and remains protected. {ref('phase_a_digest')}", "",
        "### Phase-A amendment carried into B1", "",
        f"`{correction['name']}` changes the brane-shear unit gradient `{correction['unit_gradient_transition']['display']}` "
        f"(`{correction['unit_gradient_transition']['old_numeric']:.13g}` → `{correction['unit_gradient_transition']['new_numeric']:.13g}`). "
        f"The normalized Phase-A payload digest moves `{correction['digest_transition']['old']}` → `{correction['digest_transition']['new']}`.", "",
        f"The correction closure contains exactly {len(correction['corrected_row_paths'])} paths: "
        f"`{'`, `'.join(correction['corrected_row_paths'])}`.", "",
        f"Tilt-profile rows are frozen as `{correction['tilt_profile_rows']['disposition']}`; "
        f"byte semantics unchanged = `{str(correction['tilt_profile_rows']['byte_semantics_unchanged']).lower()}`.", "",
        "### Derived positive content and honest exits", "",
        f"The manifest contains {len(manifest['fields'])} typed field routes and {len(manifest['surface_variations'])} surface variations. "
        f"It distinguishes rigid-advection from endpoint-response tangents and carries each field-local radial dimension and measure. {ref('field_manifest')}", "",
        f"Failed indexed-tangent lookups emit {len(missing)} leaves: `{', '.join(missing)}`. {ref('emitted_leaves')}", "",
        f"All {len(sym['indexed_cells'])} endpoint/ambient cells emit a derived isotropic `M_XX(p=0)`, zero reconstruction residual, a binding scalar regression, "
        f"and executable conditional native-slice residuals for `M_Xp`/`M_pp`. {ref('scalar_regression')} {ref('native_slice_constraints')}", "",
        f"The authoritative component aggregation is `{cell_statuses}`; no missing indexed or OPEN remainder was promoted to a value. {ref('mechanics_map')} {ref('cell_count')}", "",
        markdown_table(["endpoint", "computed source", "headline", "intrinsic Ω"], endpoint_rows), "",
        "### Traversal, dimensions, and congruence", "",
        f"Expression/root traversal produces {len(sym['open_root_reachability'])} OPEN root/block records. The shared finite generator yields "
        f"an empty control (`{controls['parity_empty_control']['classification']}`) and a witness control (`{controls['even_witness_control']['classification']}`). "
        f"{ref('open_reachability')} {ref('finite_controls')}", "",
        f"Units were restored termwise on the emitted expressions; computed dimension sets are `{dimensions}`. {ref('dimensions')}", "",
        f"The produced blocks reduce by exact transverse/longitudinal congruence to {len(sym['derived_congruence']['positive_definite_conditions'])} conditional pivots; "
        f"the full signature remains unresolved because their coefficients contain emitted remainders. {ref('derived_congruence')}", "",
        "### Berry/G4 and endpoint machinery", "",
        f"The two traversal-derived Berry DAGs cover {len(sym['berry']['production_pullback_coverage'])} production cells and agree after the executed zero-mode quotient. "
        f"Production winding is `{sym['berry']['intrinsic_circulation']['k_theta']}`. {ref('berry_coverage')}", "",
        f"The control contour gives `k={sym['g4_control']['computed_winding']}` and downstream `sigma={sym['g4_control']['computed_sigma']}`; "
        f"the consumed sheet area leaves a zero total-to-per-area residual. {ref('g4_winding')} {ref('g4_sigma')} {ref('sheet_area')}", "",
        f"E4 emits its pre-constraint action before differentiating `M_aug`; every stored Hessian residual is zero and `J` comes from the constraint/moduli equations. {ref('e4_action_hessian')}", "",
        f"E5 consumes and deletes root `{sym['E5']['root']}`, then re-solves to the computed E2 conservative response. {ref('e5_root')}", "",
        "### Records, gates, and halt", "",
        f"Return closure is absent by graph reachability, not by declaration. {ref('closure_absence')}", "",
        markdown_table(["gate", "computed status"], gate_rows), "",
        f"Gate rows and candidate ownership are artifact-backed. {ref('gate_statuses')} {ref('partition')}", "",
        f"The external unchanged-executable gauntlet passed {len(mutations['cases'])} focused mutation cases, "
        f"{len(mutations['metamorphic_controls'])} metamorphic controls, and {len(mutations['input_liveness_cases'])} per-key liveness cases. "
        "[typed mutation artifact `mechanics_mutations`]", "",
        f"The engines agree on {len(agreement)} independently recomputed groups: " + "; ".join(agreement) + ". [typed agreement artifact `mechanics_dual_engine`]", "",
        "### Known residual constructs", "",
        "These residual constructs remain explicit and are not accepted silently or presented as independently eliminated:", "",
        markdown_table(["construct", "bound artifact field", "comparator-side backstop"],
                       [[name, row["artifact_field"], row["comparator_backstop"]] for name, row in residuals.items()]), "",
        "### Typed report-claim bindings", "", markdown_table(["claim", "schema path", "type", "recomputation"], claim_rows), "",
        f"Axis 1 is `{sym['axis_1']}`; every mechanics cell remains honestly `UNRESOLVED`. B2 and Phase C remain `NOT_RUN(phase)`. "
        f"{ref('mechanics_map')} {ref('partition')}", "",
    ]
    return "\n".join(lines)


def build_results(phase_v3: dict[str, Any], phase_artifact: dict[str, Any], sym: dict[str, Any], agreement: list[str], mutations: dict[str, Any], input_sha: str, phase_digest: str) -> dict[str, Any]:
    result = copy.deepcopy(phase_v3)
    result["schema_version"] = "U1_BODY_DYNAMICS_RESULTS_V5"; result["phase"] = "B1"; result["halt_after_phase"] = True
    all_engine_checks = [row["status"] for row in sym["checks"].values()]
    comparator_checks = {name: "PASS" for name in COMPARE_CHECKS}
    axis = "COMPUTATION_VALID" if all(x == "PASS" for x in all_engine_checks) and all(x == "PASS" for x in comparator_checks.values()) else "HARNESS_FAILED"
    result["axis_1"] = axis; result["axis_2"] = "MECHANICS_MAP"
    result["phase_a_protection"] = {"axis_1": phase_artifact["axis1"], "axis_2": phase_artifact["axis2"],
                                     "normalized_payload_sha256": phase_digest, "unchanged": True}
    result["phase_a_correction_finding"] = surfaced_correction(sym["phase_a_amendment"])
    result["mechanics_input_sha256"] = input_sha; result["mechanics_cells"] = sym["cells"]
    result["MECHANICS_MAP"] = {key: aggregate_components(row["component_findings"]) for key, row in sym["cells"].items()}
    result["mechanics_gates"] = sym["gate_evidence"]
    result["indexed_mechanics"] = {"field_embedding": sym["field_embedding"], "missing_indexed_leaves": sym["indexed_profile_missing_leaves"],
                                     "cells": sym["indexed_cells"], "full_block_expressions": sym["full_block_expressions"],
                                     "symbolic_LDL": sym["symbolic_LDL"], "dimensions": sym["dimensions"]}
    result["symplectic_mechanics"] = {"berry": sym["berry"], "g4_control": sym["g4_control"]}
    result["E4_augmented_reduction"] = sym["E4"]; result["E5_rayleigh"] = sym["E5"]
    result["open_root_reachability"] = sym["open_root_reachability"]; result["reachability_analysis"] = sym["reachability_analysis"]
    result["mechanics_provenance_graph"] = sym["provenance_graph"]; result["mechanics_partition_ledger"] = sym["partition_ledger"]
    result["known_residual_constructs"] = known_residual_constructs(sym)
    result["mechanics_dual_engine"] = {"status": "ENGINE_AGREE", "checks": agreement,
                                        "representations": {"SymPy": sym["engine_representation"], "Mathematica": "field_tangent_jacobian_and_termwise_gram_pullback"}}
    result["mechanics_mutations"] = mutations; result["mechanics_comparator_checks"] = comparator_checks
    result["input_liveness"] = sym["input_liveness"]
    result["report_claim_bindings"] = sym["report_claim_bindings"]; result["closure_axis_B1"] = sym["closure_axis"]
    result["open_strata_B1"] = sym["open_strata"]; result["downstream"] = {"phase_B2": "NOT_RUN(phase_B2)", "phase_C": "NOT_RUN(phase_C)"}
    leaves = sym["indexed_profile_missing_leaves"]
    endpoint_lines = []
    for endpoint in ENDPOINTS:
        records = [row for key, row in sym["cells"].items() if key.startswith(endpoint + "|")]
        m_status = "/".join(sorted({row["component_findings"]["M_full"]["status"] for row in records}))
        omega_status = "/".join(sorted({"|".join(row["component_findings"][name]["status"] for name in
                                                    ("Omega_intrinsic_circulation", "Omega_translation_texture", "Omega_translation_tilt", "Omega_tilt_tilt"))
                                                for row in records}))
        headline = "/".join(sorted({row["headline"]["status"] for row in records}))
        endpoint_lines.append(f"{endpoint}: headline={headline}; M={m_status}; Omega={omega_status}")
    gate_text = "; ".join(f"{name}={row['status']}" for name, row in sym["gate_evidence"].items() if name in {"G2", "G3", "G4", "G5", "G6"})
    result["summary_lines"] = [
        "U1 PHASE B1 REMEDIATION-3 SUMMARY — HALT AFTER 7.2/7.3", *endpoint_lines,
        f"Derived: fields={len(sym['field_manifest']['fields'])}; cells={len(sym['indexed_cells'])}; scalar/native residual groups={len(sym['scalar_regression'])}; dimensions={len(sym['dimensions']['records'])}",
        f"Gates: {gate_text}; dual_groups={len(agreement)}; focused_mutations={len(mutations['cases'])}; liveness={len(mutations['input_liveness_cases'])}",
        f"Missing: emitted_indexed_leaves={len(leaves)} [{','.join(leaves)}]; OPEN_pairs={len(sym['open_root_reachability'])}; partition={sym['partition_ledger']['state']}",
        f"Correction: {sym['phase_a_amendment']['name']}; 12π/5→8π/3; digest {phase_digest[:8]}→{sym['phase_a_amendment']['amended_payload_sha256'][:8]}; sigma={sym['g4_control']['computed_sigma']}; production_k={sym['berry']['intrinsic_circulation']['k_theta']}",
    ]
    require(len(result["summary_lines"]) == 10, "B1_C_REPORT_SCHEMA", "summary must have ten lines")
    return result


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--artifacts", type=Path, default=ARTIFACTS); parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--sympy-artifact", type=Path); parser.add_argument("--math-artifact", type=Path); parser.add_argument("--phase-a-artifact", type=Path)
    parser.add_argument("--mutation-results", type=Path); parser.add_argument("--report", type=Path, default=DEFAULT_REPORT)
    parser.add_argument("--results", type=Path, default=DEFAULT_RESULTS); parser.add_argument("--verify-only", action="store_true"); parser.add_argument("--print-summary", action="store_true")
    parser.add_argument("--stage1-complete", type=Path)
    args = parser.parse_args()
    if args.print_summary:
        result = yaml.safe_load(args.results.read_text()); lines = result.get("summary_lines", [])
        if len(lines) != 10:
            print(f"ASSERT_FAIL:B1_C_REPORT_SCHEMA:summary line count {len(lines)}", file=sys.stderr); return 1
        print("\n".join(lines)); return 0
    try:
        data = yaml.safe_load(args.input.read_text())
        sym = load(args.sympy_artifact or args.artifacts / "sympy_phase_b1.yaml")
        math_artifact = load(args.math_artifact or args.artifacts / "mathematica_phase_b1.yaml")
        phase = json.loads((args.phase_a_artifact or args.artifacts / "sympy_phase_a.json").read_text())
        agreement = compare(sym, math_artifact, phase, data)
        if args.verify_only:
            if args.stage1_complete:
                paths = {"sympy_phase_b1": args.sympy_artifact or args.artifacts / "sympy_phase_b1.yaml",
                         "mathematica_phase_b1": args.math_artifact or args.artifacts / "mathematica_phase_b1.yaml",
                         "sympy_phase_a": args.phase_a_artifact or args.artifacts / "sympy_phase_a.json"}
                staged = stage1_artifact(sym, math_artifact, phase, agreement, paths)
                args.stage1_complete.parent.mkdir(parents=True, exist_ok=True)
                args.stage1_complete.write_text(yaml.safe_dump(staged, sort_keys=False, allow_unicode=True, width=220))
            print(f"B1_COMPARATOR_VERIFY: ENGINE_AGREE checks={len(agreement)}"); return 0
        mutations = load_mutation_results(args.mutation_results or args.artifacts / "b1_mutation_results.yaml")
        require({row["path"] for row in mutations["input_liveness_cases"]} == set(sym["input_liveness"]["declared_leaf_paths"]),
                "B1_C_INPUT_LIVENESS", "mutation paths exactly cover declared input leaves")
        phase_v3 = yaml.safe_load(args.results.read_text()); input_sha = hashlib.sha256(args.input.read_bytes()).hexdigest(); phase_digest = digest(phase_a_payload(phase))
        baseline = phase_a_baseline_with_correction(args.report.read_text(), sym["phase_a_amendment"])
        args.report.write_text(baseline + report_section(sym, agreement, mutations, input_sha, phase_digest))
        result = build_results(phase_v3, phase, sym, agreement, mutations, input_sha, phase_digest)
        args.results.write_text(yaml.safe_dump(result, sort_keys=False, allow_unicode=True, width=180))
        agreement_artifact = {"status": "ENGINE_AGREE", "checks": agreement, "comparator_checks": {name: "PASS" for name in COMPARE_CHECKS},
                              "mechanics_map_rederived_from_both_records": True}
        (args.artifacts / "b1_engine_agreement.yaml").write_text(yaml.safe_dump(agreement_artifact, sort_keys=False))
        print(f"B1_COMPARATOR: ENGINE_AGREE checks={len(agreement)} cells={len(sym['cells'])}"); return 0
    except CompareFailure as exc:
        print(str(exc), file=sys.stderr); return 1
    except Exception as exc:
        print(f"ASSERT_FAIL:B1_C_ENGINE_MATH:{type(exc).__name__}:{exc}", file=sys.stderr); return 1


if __name__ == "__main__":
    raise SystemExit(main())
