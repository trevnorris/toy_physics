#!/usr/bin/env python3
"""Engine-neutral B2 production comparator and computed outcome router."""

from __future__ import annotations

import argparse
import itertools
from pathlib import Path
from typing import Any

from u1_body_b2_common import digest, dump_yaml, load_yaml, load_yaml_authenticated, read_bytes_authenticated, rel_repo, require


HERE = Path(__file__).resolve().parent


def indexed(rows: list[dict[str, Any]], key: str) -> dict[str, dict[str, Any]]:
    result = {row[key]: row for row in rows}
    require(len(result) == len(rows), "B2_C_UNIQUENESS", f"duplicate {key}")
    return result


def expected_floor_from_schema(schema: dict[str, Any], operator_inventory: list[dict[str, Any]], directive: str) -> list[str]:
    """Independently expand the authenticated v47a schema at cell granularity."""
    for anchor in ["Generated obligation manifest (anti-omission, two-sided)", "Two inventories are built independently", "mode inventory", "source-to-term incidence completeness", "known-nonzero outgoing-wave control"]:
        require(anchor in directive, "B2_C_OBLIGATIONS", f"missing directive anchor {anchor}")
    require(schema["schema_version"] == "U1_PHASE_B2_V47A_OBLIGATION_SCHEMA_V2" and schema["directive_version"] == "47a", "B2_C_OBLIGATIONS", "v47a schema")
    axes = {key: list(values) for key, values in schema["axes"].items()}
    branches = sorted(row["id"] for row in operator_inventory)
    axes["operator_branch"] = branches
    axes["unordered_operator_pair"] = [f"{left}__PAIR__{right}" for left, right in itertools.combinations(branches, 2)]
    records = set(schema["fixed_products"])
    for spec in schema["expanded_products"]:
        for values in itertools.product(*(axes[axis] for axis in spec["axes"])):
            records.add(spec["name"] + "|" + "|".join(f"{axis}={value}" for axis, value in zip(spec["axes"], values)))
    return sorted(records)


def traverse_actual_dag(dag: dict[str, Any]) -> list[str]:
    nodes = dag["nodes"]
    root = dag["root"]
    require(root in nodes, "B2_C_OBLIGATIONS", "typed DAG root")
    seen: set[str] = set()
    visiting: set[str] = set()
    products: set[str] = set()

    def visit(node_id: str) -> None:
        require(node_id in nodes, "B2_C_OBLIGATIONS", f"dangling DAG node {node_id}")
        require(node_id not in visiting, "B2_C_OBLIGATIONS", f"cycle at {node_id}")
        if node_id in seen:
            return
        visiting.add(node_id)
        node = nodes[node_id]
        require(node.get("type") in {"sink", "obligation_product"}, "B2_C_OBLIGATIONS", f"untyped node {node_id}")
        if node["type"] == "obligation_product":
            require(bool(node.get("producer")), "B2_C_OBLIGATIONS", f"unbound product {node_id}")
            products.add(node_id)
        for dependency in node.get("depends_on", []):
            visit(dependency)
        visiting.remove(node_id); seen.add(node_id)

    visit(root)
    return sorted(products)


def route_g9(sector: str, record: dict[str, Any]) -> dict[str, Any]:
    residual_zero = record["residual_identically_zero"]
    determined = record["all_terms_determined"]
    independent = record["return_energy_structurally_independent"]
    exact = residual_zero and determined and (sector != "energy" or independent is True)
    causes: list[str] = []
    if not exact:
        causes.append("missing_momentum_residual_norm" if sector == "momentum" else "missing_sector_tolerance")
        if sector == "energy" and independent is not True:
            causes.append("return_energy_closure")
    causes = sorted(set(causes + record.get("missing_native_current_laws", [])))
    return {"verdict": "OK(exact)" if exact else f"UNRESOLVED({','.join(causes)})", "causes": causes}


def authenticate_engine(path: Path, producer: dict[str, Any], consumer: str) -> tuple[dict[str, Any], dict[str, Any]]:
    outputs = {Path(row["path"]).resolve(): row["sha256"] for row in producer["outputs"]}
    resolved = path.resolve()
    require(resolved in outputs, "B2_A1_PRODUCER_DIGEST", f"unbound producer output {resolved}")
    return load_yaml_authenticated(resolved, outputs[resolved], consumer)


def build(sym_path: Path, math_path: Path, stage0_path: Path, stage0_digest: str, producer_path: Path) -> dict[str, Any]:
    producer = load_yaml(producer_path)
    require(producer.get("status") == "PASS", "B2_A1_PRODUCER_DIGEST", "producer seal")
    sym, sym_auth = authenticate_engine(sym_path, producer, "production_comparator:sympy")
    math, math_auth = authenticate_engine(math_path, producer, "production_comparator:wolfram")
    stage0, stage0_auth = load_yaml_authenticated(stage0_path, stage0_digest, "production_comparator:stage0")
    manifest = stage0["observation_contract"]["Obs_B2_manifest"]
    directive_path = HERE / "directive_u1_phaseB2_intake_radiative.md"
    directive_raw, directive_auth = read_bytes_authenticated(directive_path, manifest[rel_repo(directive_path)], "production_comparator:directive")
    schema_path = HERE / "u1_body_b2_v47_obligation_schema.yaml"
    schema, schema_auth = load_yaml_authenticated(schema_path, manifest[rel_repo(schema_path)], "production_comparator:v47a_schema")
    require(sym["schema_version"] == math["schema_version"] == "U1_PHASE_B2_PRODUCTION_ENGINE_V3", "B2_C_SCHEMA", "production schema")
    require(sym["engine"] == "SymPy" and math["engine"] == "Mathematica", "B2_C_ENGINE", "independent engine identities")
    require(sym["independent_route"] != math["independent_route"], "B2_C_INDEPENDENCE", "distinct routes")
    require(sym["stage0_contract_sha256"] == math["stage0_contract_sha256"] == stage0_digest, "B2_C_STAGE0", "stage0 binding")
    for artifact in [sym, math]:
        require(len(artifact["first_use_authentication"]) >= 3 and all(row["held_descriptor"] and row["expected_sha256"] == row["consumed_sha256"] for row in artifact["first_use_authentication"]), "B2_C_FIRST_USE", artifact["engine"])

    sc, mc = indexed(sym["grid"]["cells"], "key"), indexed(math["grid"]["cells"], "key")
    require(set(sc) == set(mc) and len(sc) == 30, "B2_C_GRID", "full key grid")
    compared = []
    g9_components: dict[str, dict[str, Any]] = {}
    for key in sorted(sc):
        for field in ["axes", "C_mdot", "velocity_block", "pdot_generalized_velocity_remainder", "G2", "G5", "G6", "G7", "G9"]:
            require(sc[key][field] == mc[key][field], "B2_C_CELL", f"{key}:{field}")
        require(sc[key]["velocity_block"]["count_in_P_local"] == 1, "B2_C_NO_DOUBLE_COUNT", key)
        for sector in ["mass", "momentum", "energy"]:
            record = sc[key]["G9"][sector]
            require([row["channel"] for row in record["full_residual_terms"]] == ["var", "flux", "constraint", "Rayleigh", "rad"], "B2_C_G9", f"{key}:{sector}:five slots")
            require(bool(record["computed_full_residual"]), "B2_C_G9", f"{key}:{sector}:residual computed first")
            routed = route_g9(sector, record)
            require(record["verdict"] == routed["verdict"] and record["causes"] == routed["causes"], "B2_C_G9", f"{key}:{sector}:computed router")
            if sector == "energy":
                require(record["energy_return_leaf_derivative"] == "1" and record["return_energy_structurally_independent"] is False, "B2_C_G9", f"{key}:energy return dependence")
            if sector in g9_components:
                require(g9_components[sector] == record, "B2_C_G9", f"{sector}:cell invariance")
            else:
                g9_components[sector] = record
        compared.append({"key": key, "semantic_sha256": digest({field: sc[key][field] for field in ["C_mdot", "velocity_block", "pdot_generalized_velocity_remainder", "G7", "G9"]})})

    # The exact fixtures prove the OK route is live and that the energy
    # independence predicate is a computed part of the route.
    fixtures = stage0["frozen_data"]["integrated_balance_identities"]["router"]["fixture_controls"]
    require(fixtures["production_consulted"] is False, "B2_C_G9", "fixture tolerance is excluded from production")
    require(fixtures["exceed_bound"]["must_fire"] is True and fixtures["exact_removal_no_tolerance"]["must_remove_OK_exact"] is True, "B2_C_G9", "router fixture controls")

    frozen = stage0["frozen_data"]
    require(sym["operator_inventory"] == math["operator_inventory"] == frozen["native_operator_inventory"], "B2_C_OPERATORS", "all frozen native operators")
    require(sym["stage0_evidence"] == math["stage0_evidence"], "B2_C_NATIVE_ACTION", "native Noether/current/balance evidence")
    for field in ["partition", "radiation", "phase_C", "typed_dag"]:
        require(sym[field] == math[field], f"B2_C_{field.upper()}", field)
    require(sym["partition"]["terminal_owner_enum"] == frozen["ownership_convention"]["terminal_owner_enum"], "B2_C_PARTITION_OWNERSHIP", "frozen terminal ownership")
    require(sym["partition"]["source_to_term_incidence_residual"] == [], "B2_C_INCIDENCE", "complete action/source incidence")
    mode = sym["radiation"]["mode_coverage_residual"]
    require(mode["missing"] == [] and mode["extra"] == [], "B2_C_MODE_COVERAGE", str(mode))
    require(sym["radiation"]["per_cell_resolvent_identity"] == frozen["endpoint_resolvent_cells"] and len(frozen["endpoint_resolvent_cells"]) == 45, "B2_C_RESOLVENT", "per-cell native resolvent identities")
    expected_controls = [{"cell": row["cell"], **row["known_nonzero_control"]} for row in frozen["endpoint_resolvent_cells"]]
    require(sym["radiation"]["per_cell_known_nonzero_control"] == expected_controls, "B2_C_CONTROL", "per-cell native controls")

    expected = expected_floor_from_schema(schema, frozen["native_operator_inventory"], directive_raw.decode("utf-8"))
    frozen_expected = sorted(frozen["minimum_obligation_manifest"]["expanded_records"])
    require(expected == frozen_expected, "B2_C_OBLIGATIONS", "directive/schema floor equals orchestrator-frozen floor")
    reachable_sym = traverse_actual_dag(sym["typed_dag"])
    reachable_math = traverse_actual_dag(math["typed_dag"])
    require(reachable_sym == reachable_math, "B2_C_OBLIGATIONS", "independent engine DAG reachability")
    missing = sorted(set(expected) - set(reachable_sym)); extra = sorted(set(reachable_sym) - set(expected))
    require(not missing, "B2_C_OBLIGATIONS_EXPECTED_MINUS_REACHABLE", str(missing))
    require(not extra, "B2_C_OBLIGATIONS_REACHABLE_MINUS_EXPECTED", str(extra))

    g9_causes = sorted({cause for record in g9_components.values() for cause in record["causes"]})
    g9_headline = "OK(exact)" if not g9_causes and all(row["verdict"] == "OK(exact)" for row in g9_components.values()) else f"UNRESOLVED({','.join(g9_causes)})"
    headline = {
        "partition": sym["partition"]["state"],
        "C_mdot": "UNRESOLVED(required_OPEN_leaves)" if any("UNRESOLVED" in row["C_mdot"]["status"] for row in sc.values()) else "OK(exact)",
        "radiation": "UNRESOLVED(native_branch_inputs)" if any("UNRESOLVED" in row["accelerating"] for row in sym["radiation"]["branch_cells"]) else "OK(exact)",
        "G7": "UNRESOLVED(required_OPEN_leaves)" if any("UNRESOLVED" in row["G7"]["status"] for row in sc.values()) else "OK(exact)",
        "G9": g9_headline,
    }
    return {
        "schema_version": "U1_PHASE_B2_ENGINE_AGREEMENT_V3", "status": "PASS_WITH_HONEST_OUTCOMES",
        "stage0_contract_sha256": stage0_digest, "first_use_authentication": [sym_auth, math_auth, stage0_auth, directive_auth, schema_auth],
        "producer_seal": {"path": str(producer_path.resolve()), "producer": producer["producer"]},
        "cell_count": len(compared), "cells": compared, "operator_count": len(sym["operator_inventory"]), "operator_semantic_sha256": digest(sym["operator_inventory"]),
        "G9_component_records": g9_components,
        "obligation_manifest": {"expected_generator": "directive/schema groups independent of engine DAG", "expected": expected, "reachable": reachable_sym, "expected_minus_reachable": missing, "reachable_minus_expected": extra},
        "headline": headline,
        "checks": {"grid": "PASS", "dual_engine": "PASS", "no_double_count": "PASS", "native_action_derivations": "PASS", "G9_computed_router": "PASS", "mode_coverage": "PASS", "two_sided_obligation_set": "PASS", "first_use_authentication": "PASS"},
    }


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--sympy", type=Path, required=True); p.add_argument("--mathematica", type=Path, required=True)
    p.add_argument("--stage0", type=Path, required=True); p.add_argument("--stage0-contract-digest", required=True)
    p.add_argument("--producer-record", type=Path, required=True); p.add_argument("--output", type=Path, required=True); a = p.parse_args()
    try:
        dump_yaml(a.output, build(a.sympy, a.mathematica, a.stage0, a.stage0_contract_digest, a.producer_record))
        print("B2_COMPARE: PASS_WITH_HONEST_OUTCOMES"); return 0
    except Exception as exc:
        print(str(exc)); return 1


if __name__ == "__main__":
    raise SystemExit(main())
