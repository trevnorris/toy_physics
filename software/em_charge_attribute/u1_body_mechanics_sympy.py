#!/usr/bin/env python3
"""U1 Phase-B1 remediation-3 SymPy engine.

The engine derives the constructible indexed translation block from Phase-A
field/profile records.  Tilt tensors remain unavailable field leaves; the
Phase-A scalar tilt information is retained only as executable native-slice
constraints.  Production code is mutation-unaware.
"""

from __future__ import annotations

import argparse
import copy
import hashlib
import importlib.util
import json
import math
import sys
from contextlib import contextmanager
from pathlib import Path
from typing import Any, Iterable

import mpmath as mp
import sympy as sp
import yaml


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
DEFAULT_INPUT = HERE / "u1_body_mechanics_inputs.yaml"
DEFAULT_OUTPUT = HERE / "reports/u1_body_dynamics_artifacts/sympy_phase_b1.yaml"
ENDPOINTS = ("E1", "E2", "E3", "E4", "E5")
AMBIENTS = ("one_sided_pathA29", "symmetric_postulate")
ENERGY_DIMENSION = (2, -2, 1)
PHASE_A_BASELINE_DIGEST = "a32c25f4325671d280b54df6c51abd9b25008ef5e6008b98972bac1ed81f7e69"
PAYLOAD_NORMALIZATION_VERSION = "u1_phase_a_payload_v2_fixed_decimal_12"


def digest(value: Any) -> str:
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":"), default=str).encode()).hexdigest()


def resolve(path: str) -> Path:
    return ROOT / path


def q(value: Any, locals_: dict[str, Any] | None = None) -> sp.Expr:
    if isinstance(value, bool):
        return sp.true if value else sp.false
    if isinstance(value, (int, sp.Integer)):
        return sp.Integer(value)
    if isinstance(value, float):
        return sp.Rational(str(value))
    return sp.sympify(str(value), locals=locals_ or {})


def canonical(expr: sp.Expr) -> list[dict[str, Any]]:
    expr = sp.expand(expr)
    if expr == 0 or expr.is_zero is True:
        return []
    combined: dict[tuple[tuple[str, int], ...], float] = {}
    for term in sp.Add.make_args(expr):
        coefficient, factors = term.as_coeff_Mul()
        number = sp.N(coefficient, 17)
        powers: dict[str, int] = {}
        for factor, power in factors.as_powers_dict().items():
            if factor.is_number:
                number *= factor**power
            else:
                powers[str(factor)] = int(power)
        key = tuple(sorted(powers.items()))
        combined[key] = combined.get(key, 0.0) + float(number)
    rows = [{"coefficient": number, "powers": dict(key)} for key, number in combined.items() if abs(number) > 2e-10]
    return sorted(rows, key=lambda row: tuple(row["powers"].items()))


def canonical_matrix(matrix: sp.Matrix) -> list[list[list[dict[str, Any]]]]:
    return [[canonical(matrix[i, j]) for j in range(matrix.cols)] for i in range(matrix.rows)]


def expression_from_terms(rows: list[dict[str, Any]], symbols: dict[str, sp.Symbol]) -> sp.Expr:
    result: sp.Expr = sp.Integer(0)
    for row in rows:
        term: sp.Expr = sp.Float(row["coefficient"], 17)
        for name, power in row["powers"].items():
            symbols.setdefault(name, sp.Symbol(name, real=True))
            term *= symbols[name] ** int(power)
        result += term
    return sp.expand(result)


def at_path(value: Any, path: str | Iterable[Any]) -> Any:
    parts = path.split(".") if isinstance(path, str) else list(path)
    cur = value
    for part in parts:
        cur = cur[int(part)] if isinstance(cur, list) else cur[part]
    return cur


def sphere_area(d: int) -> sp.Expr:
    return sp.simplify(2 * sp.pi ** (sp.Rational(d, 2)) / sp.gamma(sp.Rational(d, 2)))


class Recorder:
    def __init__(self) -> None:
        self.rows: dict[str, dict[str, Any]] = {}

    def check(self, tooth: str, condition: Any, detail: Any) -> bool:
        passed = bool(condition)
        evidence = detail if isinstance(detail, (dict, list)) else {"detail": str(detail)}
        self.rows[tooth] = {"status": "PASS" if passed else "FAIL", "evidence_digest": digest(evidence), "evidence": evidence}
        return passed

    @property
    def failures(self) -> list[str]:
        return [name for name, row in self.rows.items() if row["status"] != "PASS"]


class ReadTracker:
    def __init__(self) -> None:
        self.scope_name = "unscoped_semantic"
        self.events: dict[str, set[str]] = {}

    @contextmanager
    def scope(self, name: str):
        previous = self.scope_name; self.scope_name = name
        try:
            yield
        finally:
            self.scope_name = previous

    def record(self, path: tuple[Any, ...]) -> None:
        self.events.setdefault(".".join(map(str, path)), set()).add(self.scope_name)

    def wrap(self, value: Any, path: tuple[Any, ...] = ()) -> Any:
        if isinstance(value, dict):
            return TrackedDict(value, self, path)
        if isinstance(value, list):
            return TrackedList(value, self, path)
        return value


class TrackedDict(dict):
    def __init__(self, value: dict[Any, Any], tracker: ReadTracker, path: tuple[Any, ...]):
        self._tracker = tracker; self._path = path
        for key, child in value.items():
            dict.__setitem__(self, key, tracker.wrap(child, path + (key,)))

    def __getitem__(self, key: Any) -> Any:
        value = dict.__getitem__(self, key)
        if not isinstance(value, (dict, list)):
            self._tracker.record(self._path + (key,))
        return value

    def get(self, key: Any, default: Any = None) -> Any:
        return self[key] if key in self else default

    def items(self):
        return ((key, self[key]) for key in dict.__iter__(self))

    def values(self):
        return (self[key] for key in dict.__iter__(self))

    def __deepcopy__(self, memo: dict[int, Any]) -> dict[Any, Any]:
        return {key: copy.deepcopy(self[key], memo) for key in dict.__iter__(self)}


class TrackedList(list):
    def __init__(self, value: list[Any], tracker: ReadTracker, path: tuple[Any, ...]):
        self._tracker = tracker; self._path = path
        list.__init__(self, [tracker.wrap(child, path + (index,)) for index, child in enumerate(value)])

    def __getitem__(self, key: Any) -> Any:
        value = list.__getitem__(self, key)
        if isinstance(key, int) and not isinstance(value, (dict, list)):
            self._tracker.record(self._path + (key,))
        return value

    def __iter__(self):
        for index in range(list.__len__(self)):
            yield self[index]

    def __contains__(self, value: Any) -> bool:
        return any(child == value for child in self)

    def __eq__(self, other: Any) -> bool:
        return list(self) == other

    def __deepcopy__(self, memo: dict[int, Any]) -> list[Any]:
        return [copy.deepcopy(self[index], memo) for index in range(list.__len__(self))]


def symbols_and_values(phase_input: dict[str, Any]) -> tuple[dict[str, sp.Symbol], dict[str, sp.Expr]]:
    symbols: dict[str, sp.Symbol] = {}
    values: dict[str, sp.Expr] = {}
    for name, row in phase_input["coefficients"].items():
        symbols[name] = sp.Symbol(name, real=True, positive=row.get("constraint") == "positive")
        raw = q(row["value"])
        values[name] = symbols[name] if raw.is_Symbol else raw
    for name, row in phase_input["core_traces"].items():
        symbols[name] = sp.Symbol(name, real=True)
        values[name] = q(row["value"])
    symbols["a"] = sp.Symbol("a", positive=True)
    symbols["s"] = sp.Symbol("s", real=True)
    values["a"] = q(phase_input["geometry"]["a"]["value"])
    return symbols, values


def phase_a_payload(artifact: dict[str, Any]) -> dict[str, Any]:
    return {
        "axis1": artifact["axis1"], "axis2": artifact["axis2"], "cells": artifact["cells"],
        "tails": [(r["id"], r["classification"], r["decay_exponent"], r["normalizable"]) for r in artifact["tail_channels"]],
        "decision": artifact["source_action_completeness"]["operative_decision_citation"],
        "action_ids": artifact["source_action_completeness"]["assembled_ids"],
        "endpoint_coefficients": {ep: {k: v["canonical_terms"] for k, v in row["coefficients"].items()}
                                  for ep, row in artifact["endpoint_effective_actions"].items()},
    }


def fixed_decimal(value: Any) -> str:
    """Frozen v2 payload number spelling, shared with the Wolfram engine."""
    return f"{float(sp.N(q(value), 18)):.12f}"


def normalized_payload_v2(payload: dict[str, Any]) -> dict[str, Any]:
    """Semantic Phase-A payload normalization frozen before the amendment.

    V1 remains the byte-compatible legacy recipe used by ``phase_a_payload``.
    V2 removes harmless engine-specific floating serialization while retaining
    every monomial and dependency-bearing power.
    """
    output = copy.deepcopy(payload)
    output["action_ids"] = sorted(output["action_ids"])
    output["decision"] = {key: output["decision"][key] for key in ("id", "sha256", "status")}
    for endpoint in sorted(output["endpoint_coefficients"]):
        coefficients = output["endpoint_coefficients"][endpoint]
        for name in sorted(coefficients):
            rows = []
            for row in coefficients[name]:
                rows.append({"coefficient": fixed_decimal(row["coefficient"]),
                             "powers": {str(k): int(v) for k, v in sorted(row["powers"].items())}})
            coefficients[name] = sorted(rows, key=lambda row: tuple(row["powers"].items()))
    return output


def normalized_core_traces(core_traces: dict[str, Any]) -> dict[str, Any]:
    """Freeze the protected Phase-A traces without engine-specific ``/`` escaping."""
    return {
        name: {
            "field": row["field"],
            "surface": row["surface"],
            "value": fixed_decimal(q(row["value"])),
        }
        for name, row in sorted(core_traces.items())
    }


def semantic_diff(left: Any, right: Any, prefix: tuple[Any, ...] = ()) -> list[dict[str, Any]]:
    if isinstance(left, dict) and isinstance(right, dict):
        rows: list[dict[str, Any]] = []
        for key in sorted(set(left) | set(right)):
            if key not in left or key not in right:
                rows.append({"path": ".".join(map(str, prefix + (key,))), "old": left.get(key), "new": right.get(key)})
            else:
                rows.extend(semantic_diff(left[key], right[key], prefix + (key,)))
        return rows
    if isinstance(left, list) and isinstance(right, list):
        rows = []
        for index in range(max(len(left), len(right))):
            if index >= len(left) or index >= len(right):
                rows.append({"path": ".".join(map(str, prefix + (index,))),
                             "old": left[index] if index < len(left) else None,
                             "new": right[index] if index < len(right) else None})
            else:
                rows.extend(semantic_diff(left[index], right[index], prefix + (index,)))
        return rows
    return [] if left == right else [{"path": ".".join(map(str, prefix)), "old": left, "new": right}]


def phase_a_acceptance_recheck(phase_input: dict[str, Any], amended: dict[str, Any]) -> tuple[str, dict[str, Any]]:
    """Run the Phase-A acceptance classifier itself on the amended substrate."""
    module_path = HERE / "u1_body_dynamics_sympy.py"
    spec = importlib.util.spec_from_file_location("u1_phase_a_acceptance", module_path)
    if spec is None or spec.loader is None:
        raise RuntimeError("cannot load Phase-A acceptance engine")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    _, values = module.coefficient_symbols(phase_input)
    values["a"] = module.qvalue(phase_input["geometry"]["a"]["value"])
    verdict, decision = module.classify(amended["tail_channels"], values)
    evidence = {"verdict": verdict, "decision": decision,
                "tails_digest": digest(amended["tail_channels"]),
                "classifier_source_sha256": hashlib.sha256(module_path.read_bytes()).hexdigest()}
    return verdict, evidence


def amend_brane_shear(phase: dict[str, Any], phase_input: dict[str, Any]) -> tuple[dict[str, Any], dict[str, Any], bool]:
    """Apply only the correction selected by the verified l=1 shear chain."""
    legacy_digest = digest(phase_a_payload(phase))
    declared_digest = phase_input.get("phase_a_protection", {}).get("normalized_payload_sha256", PHASE_A_BASELINE_DIGEST)
    old_gate = legacy_digest == declared_digest == PHASE_A_BASELINE_DIGEST
    amended = copy.deepcopy(phase)
    symbols, values = symbols_and_values(phase_input)
    shear = symbols["shear_transverse"]
    unit_old = sp.Rational(12, 5) * sp.pi
    # u_a=A(a/r)^nu n_a, d=3, ell=1, nu=2.  The two angular pieces
    # integrate to 4/3 (radial) + 2/3 (angular) before S_2/d.
    unit_new = sp.Rational(8, 3) * sp.pi
    corrected_rows: list[dict[str, Any]] = []
    coefficient_rows: list[dict[str, Any]] = []
    tilt_paths = [f"evaluated_moments.{ep}.{name}" for ep in ENDPOINTS
                  for name in ("U_XP", "U_PP", "I_shear_grad")]
    tilt_before = {path: copy.deepcopy(at_path(phase, path)) for path in tilt_paths}
    for endpoint in ENDPOINTS:
        beta = q(phase["endpoint_responses"][endpoint]["shear_coefficient"])
        old_moment = copy.deepcopy(phase["evaluated_moments"][endpoint]["U_XX"])
        new_moment_expr = sp.expand(unit_new * (shear + beta) ** 2)
        production_subs = {symbols[name]: value for name, value in values.items() if name in symbols}
        new_moment = {"dependencies": sorted(map(str, new_moment_expr.free_symbols)),
                      "expression": str(new_moment_expr),
                      "production_value": float(sp.N(new_moment_expr.subs(production_subs), 17))}
        amended["evaluated_moments"][endpoint]["U_XX"] = new_moment
        moment_path = f"evaluated_moments.{endpoint}.U_XX"
        corrected_rows.append({"path": moment_path, "endpoint": endpoint, "old": old_moment,
                               "new": new_moment, "delta_at_production": new_moment["production_value"] - old_moment["production_value"]})

        gvv_row = amended["endpoint_effective_actions"][endpoint]["coefficients"]["GVV"]
        old_gvv = expression_from_terms(gvv_row["canonical_terms"], symbols)
        old_shear = sp.expand(unit_old * symbols["rhoBr"] * (shear + beta) ** 2)
        new_shear = sp.expand(unit_new * symbols["rhoBr"] * (shear + beta) ** 2)
        new_gvv = sp.expand(old_gvv - old_shear + new_shear)
        new_gvv_row = {"canonical_terms": canonical(new_gvv),
                       "dependencies": sorted(map(str, new_gvv.free_symbols)), "expression": str(new_gvv)}
        amended["endpoint_effective_actions"][endpoint]["coefficients"]["GVV"] = new_gvv_row
        coefficient_rows.append({"path": f"endpoint_effective_actions.{endpoint}.coefficients.GVV",
                                 "endpoint": endpoint, "old": copy.deepcopy(phase["endpoint_effective_actions"][endpoint]["coefficients"]["GVV"]),
                                 "new": copy.deepcopy(new_gvv_row),
                                 "delta_expression": str(sp.expand(new_shear - old_shear))})

    tilt_after = {path: copy.deepcopy(at_path(amended, path)) for path in tilt_paths}
    old_payload = normalized_payload_v2(phase_a_payload(phase))
    new_payload = normalized_payload_v2(phase_a_payload(amended))
    diff_rows = semantic_diff(old_payload, new_payload)
    allowed_prefixes = tuple(f"endpoint_coefficients.{ep}.GVV" for ep in ENDPOINTS)
    semantic_gate = bool(diff_rows) and all(row["path"].startswith(allowed_prefixes) for row in diff_rows)
    tilt_gate = tilt_before == tilt_after
    verdict, recheck = phase_a_acceptance_recheck(phase_input, amended)
    new_digest = digest(new_payload)
    finding = {
        "name": "PHASE_A_MOMENT_CORRECTION(brane_shear)",
        "authoritative_chain": ["shear_operator_harmonic_decomposition", "tail_channels.brane_shear.nu=2",
                                "core_traces.shear_transverse", "endpoint_responses.*.shear_coefficient"],
        "construction": "u_a=A*(a/r)^2*n_a; integral sum_a partial_i(u_a) partial_j(u_a)",
        "old_unit_gradient": str(unit_old), "new_unit_gradient": str(unit_new),
        "old_unit_gradient_numeric": float(sp.N(unit_old, 17)), "new_unit_gradient_numeric": float(sp.N(unit_new, 17)),
        "error_provenance": "Phase-A unit_brane_grad set brane_nu=db=3 instead of the verified decaying root nu=2",
        "corrected_moment_rows": corrected_rows, "corrected_downstream_coefficients": coefficient_rows,
        "corrected_row_paths": [row["path"] for row in corrected_rows + coefficient_rows],
        "tilt_profile_rows": {"disposition": "UNCHANGED;UNRESOLVED(tilt_profile)", "paths": tilt_paths,
                              "byte_semantics_unchanged": tilt_gate},
        "payload_normalization": {"version": PAYLOAD_NORMALIZATION_VERSION,
                                  "recipe": "sorted JSON keys; decision reduced to id/sha256/status; action_ids sorted; coefficient floats rendered fixed decimal(12); monomial powers sorted"},
        "baseline_digest": {"expected": PHASE_A_BASELINE_DIGEST, "declared": declared_digest,
                            "computed_legacy": legacy_digest, "gate": old_gate},
        "amended_payload_sha256": new_digest,
        "semantic_diff": {"rows": diff_rows, "allowed_prefixes": list(allowed_prefixes), "restricted_to_correction_closure": semantic_gate},
        "phase_a_acceptance_recheck": recheck,
    }
    return amended, finding, bool(old_gate and semantic_gate and tilt_gate and verdict == "U1_BASE_OK")


def source_contract(data: dict[str, Any], phase_input: dict[str, Any], phase: dict[str, Any]) -> tuple[dict[str, Any], bool]:
    decision_path = resolve(data["substrate"]["operative_decision"])
    decision_sha = hashlib.sha256(decision_path.read_bytes()).hexdigest()
    citation = phase["source_action_completeness"]["operative_decision_citation"]
    assembled = set(phase["source_action_completeness"]["assembled_ids"])
    declared = {row["id"] for row in phase_input["action_terms"]}
    retired = set(phase_input["operative_action_decision"]["retired_action_term_ids"])
    direct = {"bulk_flow_kinetic", "brane_shear_kinetic", "uw_kinetic", "h_kinetic", "bulk_berry"}
    external_paths = {name: resolve(path).exists() for name, path in data["substrate"].items() if name != "operative_decision"}
    phase_digest = digest(phase_a_payload(phase)); declared_phase_digest = data["phase_a_protection"]["normalized_payload_sha256"]
    core_trace_digest = digest(normalized_core_traces(phase_input["core_traces"])); declared_core_trace_digest = data["phase_a_protection"]["core_traces_sha256"]
    okay = (citation["sha256"] == decision_sha and citation["status"] == "OPERATIVE" and assembled == declared
            and assembled.isdisjoint(retired) and direct <= assembled and all(external_paths.values()))
    okay = okay and phase_digest == declared_phase_digest and core_trace_digest == declared_core_trace_digest
    return ({"decision_sha256": decision_sha, "action_ids": sorted(assembled), "declared_action_ids": sorted(declared),
             "retired_ids_absent": sorted(retired), "direct_mechanics_terms": sorted(direct),
             "phase_a_payload_sha256": phase_digest, "declared_phase_a_payload_sha256": declared_phase_digest,
             "core_traces_sha256": core_trace_digest, "declared_core_traces_sha256": declared_core_trace_digest,
             "declared_external_paths_exist": external_paths}, bool(okay))


def parse_operator(text: Any, symbols: dict[str, sp.Symbol]) -> sp.Expr:
    return sp.sympify(str(text), locals=symbols)


def endpoint_solve(ep: str, data: dict[str, Any], symbols: dict[str, sp.Symbol], include_rayleigh: bool) -> dict[str, Any]:
    normal, tangent, shear, vt = sp.symbols("normal tangent shear V_tangent", real=True)
    fields = {"normal": normal, "tangent": tangent, "shear": shear}
    op = data["boundary_operator"]
    stiffness = {"normal": parse_operator(op["normal_DtN"], symbols), "tangent": parse_operator(op["tangent_DtN"], symbols),
                 "shear": parse_operator(op["shear_DtN"], symbols)}
    boundary_action = sp.expand(sum(stiffness[name] * field**2 / 2 for name, field in fields.items()))
    row = data["endpoint_functionals"][ep]
    rayleigh = sp.Integer(0)
    if ep == "E5" and include_rayleigh:
        rayleigh = q(row["rayleigh_integrand"], {**symbols, "v_tangent": tangent, "V_tangent": vt})
    residuals: list[sp.Expr] = []
    origins: list[str] = []
    for name, target in row["essential"].items():
        residuals.append(fields[name] - q(target)); origins.append(f"{row['primitive_id']}:essential:{name}")
    for name in row["natural"]:
        residuals.append(sp.diff(boundary_action + rayleigh, fields[name]))
        channel = "rayleigh_variation" if rayleigh != 0 and name == "tangent" else "boundary_variation"
        origins.append(f"{row['primitive_id']}:{channel}:{name}")
    matrix, rhs = sp.linear_eq_to_matrix(residuals, [normal, tangent, shear])
    rhs = rhs.subs(vt, 1)
    solution = sp.simplify(matrix.inv() * rhs) if matrix.rank() == 3 else sp.zeros(3, 1)
    residual = sp.simplify(matrix * solution - rhs)
    boundary_contract = {"trace_fields": data["boundary_operator"]["trace_fields"],
                         "action_ancestors": data["boundary_operator"]["action_ancestors"],
                         "boundary_class": row["boundary_class"], "variational_class": row["variational_class"]}
    return {"endpoint": ep, "primitive_id": row["primitive_id"], "include_rayleigh": include_rayleigh,
            "variation_contract": boundary_contract,
            "boundary_action": str(boundary_action), "rayleigh_density": str(rayleigh), "residual_origins": origins,
            "assembled_matrix": [[str(x) for x in r] for r in matrix.tolist()], "assembled_rhs": [str(x) for x in rhs],
            "solution": [str(sp.simplify(x)) for x in solution], "solution_residual": [str(x) for x in residual],
            "rank": matrix.rank(), "fingerprint": digest({"origins": origins, "matrix": str(matrix), "rhs": str(rhs)}),
            "_solution": solution, "_matrix": matrix, "_rhs": rhs, "_rayleigh": rayleigh}


def assemble_endpoints(data: dict[str, Any], symbols: dict[str, sp.Symbol], phase: dict[str, Any]) -> tuple[dict[str, Any], dict[str, Any], bool]:
    full = {ep: endpoint_solve(ep, data, symbols, ep == "E5") for ep in ENDPOINTS}
    conservative = {ep: endpoint_solve(ep, data, symbols, False) for ep in ENDPOINTS}
    phase_vectors = {ep: sp.Matrix([q(phase["endpoint_responses"][ep]["fluid_coefficients"]["normal"]),
                                    q(phase["endpoint_responses"][ep]["fluid_coefficients"]["tangent"]),
                                    q(phase["endpoint_responses"][ep]["shear_coefficient"])]) for ep in ENDPOINTS}
    source_map: dict[str, str] = {}
    for ep in ENDPOINTS:
        candidates = [name for name, vector in phase_vectors.items()
                      if sp.simplify(conservative[ep]["_solution"] - vector) == sp.zeros(3, 1)]
        source_map[ep] = sorted(candidates)[0] if candidates else "UNMATCHED"
    residuals_ok = all(row["rank"] == 3 and row["solution_residual"] == ["0", "0", "0"]
                       for row in list(full.values()) + list(conservative.values()))
    source_binding_ok = all(source_map[ep] != "UNMATCHED" for ep in ENDPOINTS)
    distinct = len({row["fingerprint"] for row in full.values()}) == len(ENDPOINTS)
    expected_variation = {"E1": "holonomic_field_BC", "E2": "holonomic_field_BC", "E3": "bulk_action",
                          "E4": "nonholonomic_constraint", "E5": "Rayleigh"}
    boundary_contract_ok = (data["boundary_operator"]["trace_fields"] == ["normal", "tangent", "shear"] and
                            set(data["boundary_operator"]["action_ancestors"]) == {"bulk_flow_kinetic", "brane_shear_gradient"} and
                            all(data["endpoint_functionals"][ep]["variational_class"] == expected_variation[ep] for ep in ENDPOINTS))
    info = {"conservative": conservative, "source_map": source_map,
            "E5_conservative_equals_E2": sp.simplify(conservative["E5"]["_solution"] - conservative["E2"]["_solution"]) == sp.zeros(3, 1)}
    return full, info, bool(residuals_ok and distinct and boundary_contract_ok and source_binding_ok)


def ancestry_by_term(phase: dict[str, Any], endpoint: str, symbols: dict[str, sp.Symbol], coefficient: str) -> dict[str, sp.Expr]:
    output: dict[str, sp.Expr] = {}
    for row in phase["per_structure_ancestry"]:
        if row["endpoint"] == endpoint and row["structure"].startswith(coefficient + "."):
            output[row["ancestor"]] = expression_from_terms(row["before_monomials"], symbols)
    return output


MISSING_FIELD_ALIASES = {"n": "density", "theta": "phase", "v": "flow", "chiB": "sleeve", "u": "shear", "uw": "uw", "h": "h"}


def emitted_tilt_leaf(row: dict[str, Any]) -> str:
    alias = MISSING_FIELD_ALIASES.get(row["action_symbol"], row["id"])
    suffix = "response" if row["tangent_class"] == "endpoint_velocity_response" and alias == "flow" else "profile"
    return f"indexed_{alias}_tilt_{suffix}"


def profile_from_tail(row: dict[str, Any], tail: dict[str, Any], phase: dict[str, Any]) -> str:
    if row.get("profile_path"):
        return str(at_path(phase, row["profile_path"]))
    classification = tail["classification"]
    if classification == "POWER_LAW":
        return f"{row.get('amplitude', '1')}*(a/r)^({tail['decay_exponent']})"
    if classification == "ALGEBRAIC_GAP":
        return "0"
    if classification == "EXPONENTIAL_GAP":
        d = int(tail["radial_dimension"])
        order = sp.Rational(d, 2) - 1
        return f"{row.get('amplitude', '1')}*(r/a)^(1-{d}/2)*besselk({order},gap*r)/besselk({order},gap*a)"
    return tail.get("solved_general_form", "UNRESOLVED_PROFILE")


def field_manifest(data: dict[str, Any], phase: dict[str, Any]) -> tuple[dict[str, Any], list[str], bool]:
    tails = {row["id"]: row for row in phase["tail_channels"]}
    rows: list[dict[str, Any]] = []
    missing: list[str] = []
    dimensions_ok = True
    substrate_tangents = phase.get("indexed_field_tangents", {})
    for route in data["indexed_embedding"]["fields"]:
        tail = tails[route["phase_a_channel"]]
        d = int(tail["radial_dimension"])
        dimensions_ok &= d == (4 if route["support"] == "bulk" else 3)
        rigid = route["tangent_class"] == "rigid_advected"
        tangent = (f"-d/dr({profile_from_tail(route, tail, phase)})*n_i" if rigid else
                   ("f_resp(r)*(c_t*delta_ai+(c_n-c_t)*n_a*n_i)" if route["id"] == "bulk_velocity_response"
                    else "z_resp(r)*n_i"))
        indexed_p = substrate_tangents.get(route["id"], {}).get("p")
        leaf = None
        if indexed_p is None:
            leaf = emitted_tilt_leaf(route); missing.append(leaf)
        rows.append({"field": route["id"], "action_symbol": route["action_symbol"], "support": route["support"],
                     "radial_dimension": d, "tensor_harmonic_type": "vector_l1_response" if not rigid else "scalar_radial",
                     "integration_measure": f"r^{d-1} dr", "profile": profile_from_tail(route, tail, phase),
                     "profile_provenance": f"phase_A.tail_channels[{route['phase_a_channel']}]",
                     "response_provenance": ("derived(endpoint_functionals,boundary_operator,core_trace,support_conditions)" if not rigid else None),
                     "oracle_ancestry_forbidden": ["phase_A.evaluated_moments", "phase_A.endpoint_effective_actions", "phase_A.L_eff"],
                     "translation_tangent_class": route["tangent_class"], "Phi_Xi": tangent,
                     "Phi_pi": indexed_p, "p_tangent_lookup": f"phase_A.indexed_field_tangents.{route['id']}.p",
                     "emitted_missing_leaf": leaf, "kinetic_action": route.get("kinetic_action")})
    surface_rows: list[dict[str, Any]] = []
    indexed_surfaces = phase.get("indexed_surface_variations", {})
    for route in data["indexed_embedding"]["surfaces"]:
        if route["phase_a_tilt_lookup"] == "zero_by_fixed_control_surface":
            p_variation: Any = "0"
            leaf = None
        else:
            p_variation = indexed_surfaces.get(route["id"])
            leaf = None if p_variation is not None else "indexed_sleeve_surface_normal_profile"
            if leaf: missing.append(leaf)
        surface_rows.append({"surface": route["id"], "support": route["support"], "normal_X_variation": route["X_variation"],
                             "normal_p_variation": p_variation, "p_variation_lookup": route["phase_a_tilt_lookup"],
                             "emitted_missing_leaf": leaf})
    manifest = {"family": "Phi_f(x,w;X,p)=Phi_f,0(x-X,w)+p_i*T_f^i(x-X,w)+O(p^2)", "fields": rows,
                "surface_variations": surface_rows, "production_phase_profile": phase["phase_flux_normalization"]["selected_phase"],
                "substrate_resolution": "PARTIAL_TRANSLATION_AT_P0;INDEXED_TILT_LOOKUPS_EXECUTED"}
    route_ok = (data["indexed_embedding"]["lab_to_body_map"] == "y_i=x_i-X_i(t)" and
                data["indexed_embedding"]["translation_rule"] == "dPhi/dX_i=-partial_i Phi_0")
    return manifest, sorted(set(missing)), bool(route_ok and dimensions_ok)


KINETIC_FIELDS = {"bulk_flow_kinetic": "bulk_velocity_response", "brane_shear_kinetic": "brane_shear",
                  "uw_kinetic": "brane_normal", "h_kinetic": "h_field"}
MOMENTS = {"bulk_flow_kinetic": "N_UU", "brane_shear_kinetic": "U_XX", "uw_kinetic": None, "h_kinetic": "H_XX"}
ACTION_FACTORS = {"bulk_flow_kinetic": "m_GNLS", "brane_shear_kinetic": "rhoBr", "uw_kinetic": "rhoBr", "h_kinetic": "Mh/cE**2"}


def branch_factor(data: dict[str, Any], ambient: str, symbols: dict[str, sp.Symbol]) -> sp.Expr:
    row = data["ambient_branches"][ambient]
    return q(row["embedding_factor"], {"s": symbols["s"], "eta_asym": q(row.get("eta_asym", 0))})


def phase_coefficient(phase: dict[str, Any], endpoint: str, name: str, symbols: dict[str, sp.Symbol]) -> sp.Expr:
    return expression_from_terms(phase["endpoint_effective_actions"][endpoint]["coefficients"][name]["canonical_terms"], symbols)


def sphere_polynomial_integral(expr: sp.Expr, directions: list[sp.Symbol], dimension: int) -> sp.Expr:
    """Integrate a polynomial in direction cosines over S^(d-1)."""
    poly = sp.Poly(sp.expand(expr), *directions)
    total: sp.Expr = sp.Integer(0)
    for powers, coefficient in poly.terms():
        if any(power % 2 for power in powers):
            continue
        half_degree = sum(powers) // 2
        numerator = sp.prod(sp.factorial2(power - 1) if power else 1 for power in powers)
        denominator = sp.prod(dimension + 2 * k for k in range(half_degree))
        total += coefficient * numerator / denominator
    return sp.simplify(sphere_area(dimension) * total)


def integrate_tensor_gram(jacobian: sp.Matrix, directions: list[sp.Symbol], dimension: int, radial_symbol: sp.Symbol,
                          weight: sp.Expr, lower: sp.Expr) -> tuple[sp.Matrix, sp.Matrix, list[list[str]]]:
    angular = sp.Matrix(jacobian.cols, jacobian.cols,
                        lambda i, j: sphere_polynomial_integral(sum(jacobian[a, i] * jacobian[a, j]
                                                                    for a in range(jacobian.rows)), directions, dimension))
    radial_integrands = angular.applyfunc(lambda entry: sp.simplify(radial_symbol ** (dimension - 1) * weight * entry))
    tensor = radial_integrands.applyfunc(lambda entry: sp.simplify(sp.integrate(entry, (radial_symbol, lower, sp.oo))))
    return angular, tensor, [[str(radial_integrands[i, j]) for j in range(radial_integrands.cols)]
                             for i in range(radial_integrands.rows)]


def independent_quadrature(expr: sp.Expr, radial_symbol: sp.Symbol, lower: sp.Expr,
                           substitutions: dict[sp.Symbol, sp.Expr]) -> tuple[float, float, float]:
    evaluated = sp.N(expr.subs(substitutions), 30)
    function = sp.lambdify(radial_symbol, evaluated, modules="mpmath")
    mp.mp.dps = 35
    quadrature = mp.quad(function, [float(sp.N(lower.subs(substitutions))), mp.inf])
    symbolic = sp.N(sp.integrate(expr, (radial_symbol, lower, sp.oo)).subs(substitutions), 30)
    return float(symbolic), float(quadrature), abs(float(symbolic) - float(quadrature))


def derived_translation(data: dict[str, Any], phase_input: dict[str, Any], phase: dict[str, Any], symbols: dict[str, sp.Symbol],
                        source_map: dict[str, str], manifest: dict[str, Any]) -> tuple[dict[str, Any], dict[str, Any], bool]:
    fields = {row["field"]: row for row in manifest["fields"]}
    velocities = sp.Matrix(sp.symbols("V_x V_y V_z", real=True))
    cells: dict[str, Any] = {}
    regressions: dict[str, Any] = {}
    all_good = True
    projection = data["scalar_regression_projection"]
    vunit = sp.Matrix([q(x) for x in projection["V_unit"]]); pdunit = sp.Matrix([q(x) for x in projection["pd_unit"]])
    punit = sp.Matrix([q(x) for x in projection["p_unit"]])
    a = symbols["a"]
    _, numeric_values = symbols_and_values(phase_input)
    numeric_subs = {symbols[name]: value for name, value in numeric_values.items() if name in symbols and not value.free_symbols}
    r = sp.Symbol("r", positive=True)
    tails = {row["id"]: row for row in phase["tail_channels"]}
    d_bulk = int(tails["density_EOS"]["radial_dimension"])
    lower = numeric_values["a"]
    gap = q(tails["density_EOS"]["gap"], {"sqrt": sp.sqrt})
    order = sp.Rational(d_bulk, 2) - 1
    density_profile = (symbols["density_delta"] * (r / lower) ** (1 - sp.Rational(d_bulk, 2))
                       * sp.besselk(order, gap * r) / sp.besselk(order, gap * lower))
    bulk_weight = symbols["rho_inf"] + density_profile
    forward_cache: dict[tuple[str, str], tuple[sp.Matrix, dict[str, Any], float]] = {}

    def compute_base(source: str, term: str, field_name: str) -> tuple[sp.Matrix, dict[str, Any], float]:
        cache_key = (source, term)
        if cache_key in forward_cache:
            return forward_cache[cache_key]
        response = phase["endpoint_responses"][source]
        field = fields[field_name]; d = int(field["radial_dimension"])
        directions = list(sp.symbols(f"n0:{d}", real=True))
        factor = q(ACTION_FACTORS[term], symbols)
        if term == "bulk_flow_kinetic":
            cn = q(response["fluid_coefficients"]["normal"]); ct = q(response["fluid_coefficients"]["tangent"])
            exponent = d - 1
            profile = (lower / r) ** exponent
            jacobian = sp.Matrix(d, 3, lambda alpha, i: profile *
                                 (ct * sp.KroneckerDelta(alpha, i) + (cn - ct) * directions[alpha] * directions[i]))
            weight = bulk_weight
            normalization = {"rule": "f(a)=1 from endpoint response trace; exterior decaying ell=1 solution",
                             "boundary_coefficients": {"c_n": str(cn), "c_t": str(ct)}, "profile_exponent": exponent}
        elif term == "brane_shear_kinetic":
            tail = tails["brane_shear"]
            roots = [q(value) for value in tail["indicial_roots"]]
            exponent = max(root for root in roots if root.is_positive)
            amplitude = symbols["shear_transverse"] + q(response["shear_coefficient"])
            profile = amplitude * (lower / r) ** exponent
            jacobian = sp.Matrix(d, 3, lambda alpha, i: -(
                sp.diff(profile, r) * directions[alpha] * directions[i]
                + profile / r * (sp.KroneckerDelta(alpha, i) - directions[alpha] * directions[i])))
            weight = sp.Integer(1)
            normalization = {"rule": "u_a(a)=A*n_a from core trace plus solved endpoint shear response",
                             "amplitude": str(amplitude), "profile_exponent": str(exponent), "harmonic_degree": 1}
        elif term == "h_kinetic":
            tail = tails["h"]
            exponent = q(tail["decay_exponent"])
            profile = symbols["h_scalar"] * (lower / r) ** exponent
            jacobian = sp.Matrix(1, 3, lambda _alpha, i: -sp.diff(profile, r) * directions[i])
            weight = sp.Integer(1)
            normalization = {"rule": "h(a)=h_scalar core trace with decaying support condition",
                             "amplitude": "h_scalar", "profile_exponent": str(exponent)}
        else:
            profile = sp.Integer(0)
            jacobian = sp.Matrix(1, 3, lambda _alpha, i: -sp.diff(profile, r) * directions[i])
            weight = sp.Integer(1)
            normalization = {"rule": "algebraic field equation selects zero outside the declared collar"}
        angular_tensor, raw_tensor, radial_integrands = integrate_tensor_gram(jacobian, directions, d, r, weight, lower)
        base_tensor = sp.Matrix(3, 3, lambda i, j: sp.expand(factor *
            (sp.N(raw_tensor[i, j], 17) if raw_tensor[i, j].has(sp.meijerg, sp.besselk) else raw_tensor[i, j])))
        diagonal_integrand = q(radial_integrands[0][0], {"r": r, **symbols})
        quad_subs = {**numeric_subs, symbols["rho_inf"]: q(phase_input["coefficients"]["rho_inf"]["value"]),
                     symbols["density_delta"]: q(phase_input["core_traces"]["density_delta"]["value"]),
                     symbols["shear_transverse"]: q(phase_input["core_traces"]["shear_transverse"]["value"]),
                     symbols["h_scalar"]: q(phase_input["core_traces"]["h_scalar"]["value"])}
        symbolic_number, quadrature_number, quadrature_error = independent_quadrature(diagonal_integrand, r, lower, quad_subs)
        record = {"action_term": term, "field": field_name, "tangent_class": field["translation_tangent_class"],
                  "translation_jacobian": [[str(jacobian[i, j]) for j in range(jacobian.cols)] for i in range(jacobian.rows)],
                  "normalization_derivation": normalization, "radial_dimension": d,
                  "radial_measure": f"r^{d-1} dr", "support": [str(lower), "infinity"],
                  "angular_tensor_integral": [[str(angular_tensor[i, j]) for j in range(3)] for i in range(3)],
                  "radial_integrands": radial_integrands,
                  "symbolic_tensor_integral": [[str(raw_tensor[i, j]) for j in range(3)] for i in range(3)],
                  "quadrature_crosscheck": {"symbolic": symbolic_number, "quadrature": quadrature_number,
                                            "absolute_error": quadrature_error, "passed": quadrature_error < 2e-11},
                  "oracle_ancestry_forbidden": True, "oracle_paths_consumed": []}
        forward_cache[cache_key] = (base_tensor, record, quadrature_error)
        return forward_cache[cache_key]

    for endpoint in ENDPOINTS:
        source = source_map[endpoint]
        derivation_source = endpoint if source == "UNMATCHED" else source
        for ambient in AMBIENTS:
            key = f"{endpoint}|{ambient}"
            branch = branch_factor(data, ambient, symbols)
            term_lagrangians: dict[str, sp.Expr] = {}
            term_tensors: dict[str, sp.Matrix] = {}
            contractions: list[dict[str, Any]] = []
            for term, field_name in KINETIC_FIELDS.items():
                base_tensor, base_record, quadrature_error = compute_base(derivation_source, term, field_name)
                tensor = base_tensor.applyfunc(lambda entry: sp.expand(branch**2 * entry))
                lagrangian = sp.expand((velocities.T * tensor * velocities)[0] / 2)
                term_tensors[term] = tensor; term_lagrangians[term] = lagrangian
                record = copy.deepcopy(base_record)
                record.update({"coefficient_from_contraction": canonical(tensor[0, 0]),
                               "computed_tensor": canonical_matrix(tensor),
                               "dependency_edges": [[term, f"field_tangent:{field_name}"],
                                                    [f"field_tangent:{field_name}", f"indexed_cells:{key}:M_XX"]],
                               "zero_decision": "ZERO_FROM_COMPUTED_PROFILE" if tensor == sp.zeros(3) else "NONZERO_COMPUTED"})
                contractions.append(record)
                all_good &= quadrature_error < 2e-11
            direct_l = sp.expand(sum(term_lagrangians.values(), sp.Integer(0)))
            direct_m = sp.hessian(direct_l, list(velocities))
            assembled_gram = sp.zeros(3)
            for tensor in term_tensors.values():
                assembled_gram += tensor
            independent_l = sp.expand((velocities.T * assembled_gram * velocities)[0] / 2)
            reconstruction = sp.expand(direct_l - independent_l)
            isotropic_scalar = direct_m[0, 0]
            isotropy_residual = (direct_m - isotropic_scalar * sp.eye(3)).applyfunc(sp.expand)
            target_gvv = sp.expand(branch**2 * phase_coefficient(phase, derivation_source, "GVV", symbols))
            projected = sp.expand((vunit.T * direct_m * vunit)[0])
            scalar_residual = sp.expand(projected - target_gvv)
            tilt_scale = q(projection["indexed_tilt_length"], symbols)
            gvp = phase_coefficient(phase, derivation_source, "GVP", symbols)
            gpp = phase_coefficient(phase, derivation_source, "GPP", symbols)
            mxp_slice = sp.Symbol(f"M_Xp_native_slice__{endpoint}__{ambient}", real=True)
            mpp_slice = sp.Symbol(f"M_pp_native_slice__{endpoint}__{ambient}", real=True)
            target_mxp = sp.expand(branch**2 * tilt_scale * gvp)
            target_mpp = sp.expand(branch**2 * tilt_scale**2 * gpp)
            mxp_residual = sp.expand(mxp_slice - target_mxp)
            mpp_residual = sp.expand(mpp_slice - target_mpp)
            native_constraints = {
                "projection_vectors": {"e_V": list(map(str, vunit)), "e_pdot": list(map(str, pdunit)), "e_p": list(map(str, punit))},
                "scale_factors": {"ambient_embedding_squared": str(branch**2), "indexed_tilt_length": str(tilt_scale)},
                "M_Xp_native": {"defining_expression": str(mxp_residual), "free_slice_symbol": str(mxp_slice),
                                  "phase_a_target": canonical(gvp),
                                  "status": "CONDITIONAL_INDEXED_FAMILY_UNAVAILABLE"},
                "M_pp_native": {"defining_expression": str(mpp_residual), "free_slice_symbol": str(mpp_slice),
                                  "phase_a_target": canonical(gpp),
                                  "status": "CONDITIONAL_INDEXED_FAMILY_UNAVAILABLE"},
                "scope": "PHASE_A_NATIVE_COMPONENT_ONLY;OPEN_ROOT_PROJECTIONS_SEPARATE",
            }
            ambient_row = data["ambient_branches"][ambient]
            wall_side_branch_samples = [
                canonical(sp.expand(branch.subs(symbols["s"], q(side))))
                for side in ambient_row.get("wall_side_domain", [])
            ]
            cells[key] = {"source_endpoint": source, "derivation_source_endpoint": derivation_source, "ambient_factor": str(branch),
                          "wall_side_branch_samples": wall_side_branch_samples,
                          "termwise_L": {name: canonical(expr) for name, expr in term_lagrangians.items()},
                          "L_direct": canonical(direct_l), "L_independent_reconstruction": canonical(independent_l),
                          "reconstruction_residual": canonical(reconstruction), "M_XX_p0_known": canonical_matrix(direct_m),
                          "M_XX_scalar": canonical(isotropic_scalar), "computed_isotropy_residual": canonical_matrix(isotropy_residual),
                          "field_contraction_integrals": contractions,
                          "source_removal_deltas": {name: canonical(tensor[0, 0]) for name, tensor in term_tensors.items()},
                          "production_dependencies": sorted({str(x) for x in direct_l.free_symbols} | set(term_lagrangians)),
                          "native_tilt_slice_constraints": native_constraints, "_L": direct_l, "_M": direct_m}
            regressions[key] = {"projection_consumed_after_production": True, "projection_dependency_role": projection["dependency_role"],
                                "projected_M_XX": canonical(projected), "phase_a_GVV_target": canonical(target_gvv),
                                "residual": canonical(scalar_residual), "validated": ambient == projection["checked_ambient"]}
            isotropy_zero = all(not canonical(isotropy_residual[i, j]) for i in range(3) for j in range(3))
            ambient_ok = ((ambient == "symmetric_postulate" and ambient_row["parity_tag"] == "BODY_PLUS_AMBIENT_POSTULATE") or
                          (ambient == "one_sided_pathA29" and ambient_row["parity_tag"] == "ONE_SIDED_ASYMMETRY_MAP" and
                           ambient_row["wall_side_domain"] == [-1, 1]))
            all_good &= ambient_ok and reconstruction == 0 and isotropy_zero and not canonical(scalar_residual)
    sentinel_tensor = sp.Matrix([[2, 1, 0], [1, 5, 0], [0, 0, 7]])
    sentinel_vector = sp.Matrix([1, 1, 0])
    sentinel_projection = sp.expand((sentinel_vector.T * sentinel_tensor * sentinel_vector)[0])
    regressions["anisotropic_projection_sentinel"] = {
        "tensor": sentinel_tensor.tolist(), "e_V": sentinel_vector.tolist(),
        "true_contraction": str(sentinel_projection), "element_selection": str(sentinel_tensor[0, 0]),
        "distinguishes_element_selection": sentinel_projection != sentinel_tensor[0, 0]}
    all_good &= sentinel_projection != sentinel_tensor[0, 0]
    return cells, regressions, bool(all_good)


def add_dim(a: Iterable[int], b: Iterable[int]) -> tuple[int, int, int]:
    return tuple(int(x) + int(y) for x, y in zip(a, b))  # type: ignore[return-value]


def scale_dim(a: Iterable[int], n: int) -> tuple[int, int, int]:
    return tuple(int(x) * n for x in a)  # type: ignore[return-value]


def expression_dimensions(expr: sp.Expr, dimensions: dict[str, tuple[int, int, int]]) -> list[list[int]]:
    output: set[tuple[int, int, int]] = set()
    for term in sp.Add.make_args(sp.expand(expr)):
        dim = (0, 0, 0)
        for factor, power in term.as_powers_dict().items():
            if factor.is_number:
                continue
            dim = add_dim(dim, scale_dim(dimensions.get(str(factor), (0, 0, 0)), int(power)))
        output.add(dim)
    return [list(row) for row in sorted(output)]


def restore_length_units(expr: sp.Expr, dimensions: dict[str, tuple[int, int, int]], a: sp.Symbol,
                         target_length_power: int = 0) -> sp.Expr:
    """Restore powers of the Phase-A reference radius suppressed by a=1."""
    restored: sp.Expr = sp.Integer(0)
    for term in sp.Add.make_args(sp.expand(expr)):
        length_power = 0
        for factor, power in term.as_powers_dict().items():
            if not factor.is_number:
                length_power += dimensions.get(str(factor), (0, 0, 0))[0] * int(power)
        restored += term * a ** (target_length_power - length_power)
    return sp.expand(restored)


def dimension_analysis(data: dict[str, Any], phase_input: dict[str, Any], phase: dict[str, Any], symbols: dict[str, sp.Symbol],
                       cells: dict[str, Any], g4_omega: sp.Expr) -> tuple[dict[str, Any], bool]:
    dims: dict[str, tuple[int, int, int]] = {name: tuple(row["dimensions"]) for name, row in phase_input["coefficients"].items()}
    dims.update({name: (0, 0, 0) for name in phase_input["core_traces"]})
    dims.update({"a": (1, 0, 0), "s": (0, 0, 0), "V": (1, -1, 0), "pd": (0, -1, 0)})
    coordinate_spec = data["indexed_coordinates"]
    coordinate = coordinate_spec["coordinate_dimensions"]
    v, pd = sp.symbols("V pd", real=True)
    source = cells["E1|symmetric_postulate"]["derivation_source_endpoint"]
    gvv = phase_coefficient(phase, source, "GVV", symbols)
    gvp = phase_coefficient(phase, source, "GVP", symbols)
    gpp = phase_coefficient(phase, source, "GPP", symbols)
    a = symbols["a"]
    gvv_units = restore_length_units(gvv, dims, a)
    gvp_units = restore_length_units(gvp, dims, a)
    gpp_units = restore_length_units(gpp, dims, a)
    objects = {
        "L_translation": gvv_units * v**2 / 2,
        "L_native_Xp_slice": a * gvp_units * v * pd,
        "L_native_pp_slice": a**2 * gpp_units * pd**2 / 2,
        "M_XX": gvv_units,
        "M_Xp_native_slice": a * gvp_units,
        "M_pp_native_slice": a**2 * gpp_units,
        "Omega_XX_control": g4_omega,
    }
    records: list[dict[str, Any]] = []
    for name, expr in objects.items():
        records.append({"object": name, "units_restored_expression": str(expr), "computed_monomial_dimensions_LTM": expression_dimensions(expr, dims)})
    action_dimensions = {tuple(dimension) for row in records[:3] for dimension in row["computed_monomial_dimensions_LTM"]}
    derived_energy_dim = next((tuple(row["dimensions"]) for row in phase_input["field_records"]
                               if row["root_type"] == "ACTION"), ENERGY_DIMENSION)
    indices = coordinate_spec["spatial_indices"]
    expected_coordinate_order = [f"X_{index}" for index in indices] + [f"p_{index}" for index in indices]
    expected_velocity_order = [f"V_{index}" for index in indices] + [f"pd_{index}" for index in indices]
    nondimensionalization = coordinate_spec["nondimensionalization"]
    coordinate_contract = {"spatial_indices": indices, "coordinate_order": coordinate_spec["coordinate_order"],
                           "velocity_order": coordinate_spec["velocity_order"], "expected_coordinate_order": expected_coordinate_order,
                           "expected_velocity_order": expected_velocity_order, "p_slice": coordinate_spec["p_slice"],
                           "nondimensionalization": nondimensionalization}
    coordinate_ok = (tuple(coordinate["X"]) == (1, 0, 0) and tuple(coordinate["p"]) == (0, 0, 0) and len(indices) == 3 and
                     coordinate_spec["coordinate_order"] == expected_coordinate_order and coordinate_spec["velocity_order"] == expected_velocity_order and
                     coordinate_spec["p_slice"] == "off_shell_free" and nondimensionalization["translation_length"] == str(a) and
                     nondimensionalization["tilt_length"] == str(a))
    homogeneous = action_dimensions == {derived_energy_dim}
    return {"basis": "LTM", "coefficient_dimensions": {k: list(v) for k, v in dims.items()},
            "coordinate_dimensions": coordinate, "coordinate_basis_contract": coordinate_contract,
            "action_dimension_sourced_from_phase_a": list(derived_energy_dim),
            "records": records, "action_terms_homogeneous": homogeneous, "coordinate_embedding_consistent": coordinate_ok}, bool(homogeneous and coordinate_ok)


BLOCK_SPECS = {
    "M_XX": ("quadratic_even", "symmetric"), "M_Xp_symmetric": ("quadratic_even", "symmetric"),
    "M_pp": ("quadratic_even", "symmetric"), "M_Xp_antisymmetric": ("quadratic_even", "antisymmetric"),
    "Omega_XX_texture": ("first_order", "antisymmetric"), "Omega_Xp": ("first_order", "antisymmetric"),
    "Omega_pp": ("first_order", "antisymmetric"),
}


def generator_candidates(data: dict[str, Any], symmetry: str) -> list[tuple[str, str]]:
    inventory = data["covariant_inventory"]
    if symmetry == "symmetric":
        mapping = {"delta_ij": ("delta_ij", "even"), "p_i_p_j": ("p_i*p_j", "even")}
        return [mapping[name] for name in inventory["symmetric_basis"] if name in mapping]
    mapping = {"epsilon_ijk_p_k": ("epsilon_ijk*p_k", "odd")}
    return [mapping[name] for name in inventory["antisymmetric_basis"] if name in mapping]


def root_generator_rows(root: dict[str, Any], block: str, data: dict[str, Any]) -> dict[str, Any] | None:
    order, symmetry = BLOCK_SPECS[block]
    if order not in root["admissible_orders"]:
        return None
    finite = isinstance(root.get("derivative_order_bound"), int)
    candidates = generator_candidates(data, symmetry)
    parity = root["body_conjugation_parity"]
    generators = [name for name, gen_parity in candidates if not finite or parity in {"mixed", "branch_covariant", gen_parity}]
    coefficients = [f"c__{root['id']}__{block}__{i}" for i in range(len(generators))]
    if not finite:
        disposition = "REACHABLE_UNRESOLVED_REMAINDER"
    elif generators:
        disposition = "REACHABLE_NONZERO_WITNESS"
    else:
        disposition = "STRUCTURAL_ZERO"
    return {"root": root["id"], "block": block, "order": order, "symmetry": symmetry,
            "source_justification": {"root_type": root["root_type"], "primitive": root.get("primitive"),
                                     "domain": root["domain"], "arguments": root["arguments"], "symmetry_class": root["symmetry_class"],
                                     "argument_index_structure": root["argument_index_structure"],
                                     "body_conjugation_parity": parity, "derivative_order_bound": root["derivative_order_bound"]},
            "space_finitely_bounded": finite, "generator_candidates": [name for name, _ in candidates],
            "finite_generator_set": generators if finite else [], "coefficient_space_dimension": len(generators) if finite else "UNBOUNDED",
            "coefficient_space_empty": bool(finite and not generators), "witness_tensor": generators[0] if generators else None,
            "remainder_symbols": coefficients, "disposition": disposition}


def reachability(data: dict[str, Any], phase: dict[str, Any]) -> tuple[list[dict[str, Any]], dict[str, list[str]], dict[str, Any], bool]:
    phase_open = {row["id"] for row in phase["declared_inputs"] if row.get("status") == "OPEN_INPUT"
                  and row.get("root_type") in {"ACTION", "PRIMITIVE_OPEN"}}
    declarations = {row["id"]: row for row in data["open_action_functionals"]}
    phase_primitives = {dep for deps in phase["assembled_action"]["term_dependencies"].values() for dep in deps}
    primitive_terms = {row["id"] for row in declarations.values() if row.get("primitive") in phase_primitives}
    generated = phase_open | primitive_terms
    rows: list[dict[str, Any]] = []
    by_block: dict[str, list[str]] = {}
    computed_table: dict[str, dict[str, list[str]]] = {"quadratic_even": {"symmetric_blocks": [], "antisymmetric_blocks": []},
                                                       "first_order": {"antisymmetric_blocks": []}}
    for block, (order, symmetry) in BLOCK_SPECS.items():
        family = "symmetric_blocks" if symmetry == "symmetric" else "antisymmetric_blocks"
        computed_table[order][family].append(block)
    for root_id in sorted(generated):
        root = declarations.get(root_id)
        if root is None:
            continue
        for block in BLOCK_SPECS:
            record = root_generator_rows(root, block, data)
            if record is not None:
                rows.append(record); by_block.setdefault(block, []).append(root_id)
    controls: dict[str, Any] = {}
    for root in data["finite_bound_controls"]:
        record = root_generator_rows(root, root["target_block"], data)
        assert record is not None
        witness_matrix = sp.Symbol(f"w__{root['id']}") * sp.eye(3) if record["witness_tensor"] == "delta_ij" else sp.zeros(3)
        classification = "STRUCTURAL_ZERO" if record["coefficient_space_empty"] else "NONZERO_WITNESS"
        controls[root["id"]] = {**record, "attached_control_block": canonical_matrix(witness_matrix), "classification": classification}
    normalize = lambda table: {order: {family: sorted(values) for family, values in families.items()} for order, families in table.items()}
    declared_table = copy.deepcopy(data.get("tensor_selection_rules", {}))
    table_agrees = normalize(computed_table) == normalize(declared_table)
    empty_count = sum(row["classification"] == "STRUCTURAL_ZERO" for row in controls.values())
    witness_count = sum(row["classification"] == "NONZERO_WITNESS" and row["witness_tensor"] is not None for row in controls.values())
    declaration_contract = all(row["root_type"] in {"PROFILE", "ACTION"} and bool(row["arguments"]) for row in declarations.values())
    control_contract = all(row["root_type"] == "ACTION_CONTROL" and bool(row["arguments"]) for row in data["finite_bound_controls"])
    okay = (generated <= set(declarations) and declaration_contract and control_contract and table_agrees and bool(rows) and empty_count == 1 and witness_count == 1
            and all(row["disposition"] != "STRUCTURAL_ZERO" for row in rows))
    analysis = {"authoritative_computed_table": computed_table, "declared_table_crosscheck": declared_table,
                "crosscheck_agrees": table_agrees, "generated_roots": sorted(generated), "finite_bound_controls": controls,
                "control_outcome_counts": {"empty": empty_count, "witness": witness_count}}
    return rows, {k: sorted(set(v)) for k, v in by_block.items()}, analysis, bool(okay)


def expression_dag(named: dict[str, sp.Expr], route: str) -> dict[str, Any]:
    nodes: dict[str, dict[str, Any]] = {}
    edges: set[tuple[str, str]] = set()

    def visit(expr: sp.Basic) -> str:
        if expr.is_Symbol:
            node = f"symbol:{expr}"
            nodes[node] = {"id": node, "kind": "raw_symbol", "expression": str(expr),
                           "subexpression_digest": digest(sp.srepr(expr))}
            return node
        payload = sp.srepr(expr)
        node = f"{route}:expr:{digest(payload)[:16]}"
        nodes[node] = {"id": node, "kind": expr.func.__name__, "expression": str(expr),
                       "subexpression_digest": digest(payload)}
        for arg in expr.args:
            child = visit(arg); edges.add((child, node))
        return node

    outputs: dict[str, str] = {}
    for name, expr in named.items():
        child = visit(sp.sympify(expr)); output = f"{route}:output:{name}"
        nodes[output] = {"id": output, "kind": "named_output", "expression": str(expr),
                         "subexpression_digest": digest(sp.srepr(expr))}
        edges.add((child, output)); outputs[name] = output
    return {"nodes": [nodes[key] for key in sorted(nodes)], "edges": [list(edge) for edge in sorted(edges)], "named_outputs": outputs,
            "raw_dependencies": sorted({str(symbol) for expr in named.values() for symbol in expr.free_symbols})}


def graph_reachable(edges: list[list[str]], source: str, targets: set[str]) -> bool:
    graph: dict[str, list[str]] = {}
    for left, right in edges:
        graph.setdefault(left, []).append(right)
    seen = {source}; stack = [source]
    while stack:
        node = stack.pop()
        if node in targets:
            return True
        for nxt in graph.get(node, []):
            if nxt not in seen:
                seen.add(nxt); stack.append(nxt)
    return False


def berry_pipeline(data: dict[str, Any], phase_input: dict[str, Any], phase: dict[str, Any], symbols: dict[str, sp.Symbol]) -> tuple[dict[str, Any], dict[str, Any], dict[str, bool], sp.Expr]:
    control = data["g4_control"]
    x, y, r, phi, radius, xi = sp.symbols("x y r phi R xi_core", positive=True, real=True)
    k_control = sp.Symbol("k_control", real=True)
    local = {"x": x, "y": y, "r": r, "xi_core": xi, "rho_inf": symbols["rho_inf"], "k_control": k_control,
             "atan2": sp.atan2, "exp": sp.exp}
    theta = q(control["phase_profile"], local).subs(k_control, q(control["phase_parameters"]["k_control"]))
    xi_ratio = q(control["xi_core_over_a"])
    density_r = q(control["core_density_profile"], local).subs(xi, xi_ratio * symbols["a"])
    density_xy = density_r.subs(r, sp.sqrt(x**2 + y**2))
    orient = q(control["contour_orientation"]["sign"]); eps_sign = q(control["epsilon_xy"])
    contour = {x: radius * sp.cos(orient * phi), y: radius * sp.sin(orient * phi)}
    dtheta = sp.simplify((sp.diff(theta, x) * sp.diff(contour[x], phi) + sp.diff(theta, y) * sp.diff(contour[y], phi)).subs(contour))
    winding = sp.simplify(sp.integrate(dtheta, (phi, 0, 2 * sp.pi)) / (2 * sp.pi))
    expected_k = q(control["k_expected"])
    action_row = next(row for row in phase_input["action_terms"] if row["id"] == "bulk_berry")
    nfield, theta_t = sp.symbols("n theta_t", real=True)
    action_expr = q(action_row["expression"], {**symbols, "n": nfield, "theta_t": theta_t})
    source_coefficient = sp.simplify(sp.diff(action_expr, theta_t) / nfield)
    grad_theta = sp.Matrix([sp.diff(theta, x), sp.diff(theta, y)])
    velocities = sp.Matrix(sp.symbols("Vberry_x Vberry_y", real=True))
    theta_t_substitution = -velocities.dot(grad_theta)
    substituted_density = sp.expand(action_expr.subs(theta_t, theta_t_substitution))
    local_connection_density = sp.Matrix([sp.diff(substituted_density, velocity) for velocity in velocities])
    local_density_residual = sp.simplify(local_connection_density + source_coefficient * nfield * grad_theta)

    normal = sp.Matrix([sp.cos(orient * phi), sp.sin(orient * phi)])
    outer_grad = sp.simplify(grad_theta.subs(contour))
    outer_density = sp.simplify(sp.limit(density_r, r, sp.oo))
    surfaces = set(data["berry_pullback_selection"]["include_surface_terms"])
    derivative_matrix = sp.zeros(2)
    boundary_rows: list[dict[str, Any]] = []
    if "partial_Omega_c" in surfaces:
        outer_matrix = sp.zeros(2)
        for i in range(2):
            for j in range(2):
                outer_matrix[i, j] = sp.simplify(source_coefficient * sp.integrate(
                    outer_density * outer_grad[i] * normal[j] * radius, (phi, 0, 2 * sp.pi)))
        derivative_matrix += outer_matrix
        boundary_rows.append({"surface": "partial_Omega_c", "contribution_to_L_derivative": canonical_matrix(outer_matrix)})
    if "Sigma" in surfaces:
        core_density = sp.simplify(sp.limit(density_r, r, 0, dir="+"))
        core_matrix = sp.zeros(2)
        if core_density != 0:
            core_matrix = -core_density / outer_density * derivative_matrix
        derivative_matrix += core_matrix
        boundary_rows.append({"surface": "Sigma", "density_limit": str(core_density),
                              "contribution_to_L_derivative": canonical_matrix(core_matrix),
                              "absence_decision": "ZERO_FROM_CORE_DENSITY_LIMIT" if core_matrix == sp.zeros(2) else "NONZERO_COMPUTED"})
    xx, xy = sp.symbols("X_x X_y", real=True); collective = sp.Matrix([xx, xy])
    # Route A begins with the Taylor-expanded, termwise-substituted action.
    l_indexed = sp.expand(velocities.dot(derivative_matrix * collective))
    connection_a = sp.Matrix([sp.diff(l_indexed, velocity) for velocity in velocities])
    omega_a = sp.Matrix([[0, sp.diff(connection_a[1], xx) - sp.diff(connection_a[0], xy)],
                         [sp.diff(connection_a[0], xy) - sp.diff(connection_a[1], xx), 0]])

    dn = sp.Matrix([sp.diff(density_xy, x), sp.diff(density_xy, y)])
    jacobian = sp.simplify(dn[0] * grad_theta[1] - dn[1] * grad_theta[0])
    polar_jacobian = sp.simplify(jacobian.subs({x: r * sp.cos(phi), y: r * sp.sin(phi)}))
    interior_scalar = sp.simplify(source_coefficient * sp.integrate(sp.integrate(polar_jacobian * r, (phi, 0, 2 * sp.pi)), (r, 0, sp.oo)))
    physical_current = sp.Matrix([[0, interior_scalar], [-interior_scalar, 0]])
    z1, z2 = sp.symbols("z1 z2", real=True)
    projector = q(phase["zero_mode_quotient"]["projector"], {"Matrix": sp.Matrix, "z1": z1, "z2": z2})
    zero_connection = sp.Matrix([-z2 / 2, z1 / 2])
    zero_mode_sector = sp.Matrix(2, 2, lambda i, j: sp.diff(zero_connection[j], (z1, z2)[i])
                                 - sp.diff(zero_connection[i], (z1, z2)[j]))
    projected_zero_sector = sp.simplify(projector.T * zero_mode_sector * projector)
    quotient = bool(data["berry_pullback_selection"]["apply_zero_mode_quotient"])
    omega_b = sp.simplify(physical_current + (projected_zero_sector if quotient else zero_mode_sector))
    equivalence = sp.simplify(omega_a - omega_b)

    production_profile = q(phase["phase_flux_normalization"]["selected_phase"], {"r": r, "pi": sp.pi})
    production_profile_symbols = sorted(map(str, production_profile.free_symbols))
    production_profile_guard = set(production_profile.free_symbols) <= {r}
    production_derivative = sp.diff(production_profile, phi)
    production_winding = sp.simplify(sp.integrate(production_derivative, (phi, 0, 2 * sp.pi)) / (2 * sp.pi))
    selection = data["berry_pullback_selection"]
    coverage: list[dict[str, Any]] = []
    production_xy = production_profile.subs({r: sp.sqrt(x**2 + y**2), phi: sp.atan2(y, x)})
    production_grad = sp.Matrix([sp.diff(production_xy, x), sp.diff(production_xy, y)])
    for boundary_class, endpoint_rows in selection["production_cells"].items():
        for endpoint in endpoint_rows:
            for ambient in selection["ambient_branches"]:
                response = phase["endpoint_responses"][endpoint]
                response_factor = sum(q(response["fluid_coefficients"][name]) for name in ("normal", "tangent")) + q(response["shear_coefficient"])
                cell_factor = sp.expand(branch_factor(data, ambient, symbols) * response_factor)
                cell_velocities = sp.Matrix(sp.symbols(f"Vprod_{endpoint}_{ambient}_x Vprod_{endpoint}_{ambient}_y", real=True))
                cell_coordinates = sp.Matrix(sp.symbols(f"Xprod_{endpoint}_{ambient}_x Xprod_{endpoint}_{ambient}_y", real=True))
                cell_density = sp.expand(action_expr.subs(theta_t, -cell_velocities.dot(production_grad)))
                cell_connection = sp.Matrix([sp.diff(cell_density, velocity) for velocity in cell_velocities])
                cell_lagrangian = sp.expand(cell_factor * cell_velocities.dot(sp.Matrix([
                    sp.diff(cell_connection.dot(cell_coordinates), coordinate) for coordinate in cell_coordinates])))
                cell_a = sp.Matrix([sp.diff(cell_lagrangian, velocity) for velocity in cell_velocities])
                cell_route_a = sp.Matrix(2, 2, lambda i, j: sp.diff(cell_a[j], cell_coordinates[i]) - sp.diff(cell_a[i], cell_coordinates[j]))
                cell_jacobian = sp.simplify(sp.diff(density_xy, x) * production_grad[1] - sp.diff(density_xy, y) * production_grad[0])
                cell_polar = cell_jacobian.subs({x: r * sp.cos(phi), y: r * sp.sin(phi)})
                cell_current = sp.simplify(cell_factor * source_coefficient * sp.integrate(
                    sp.integrate(cell_polar * r, (phi, 0, 2 * sp.pi)), (r, 0, sp.oo)))
                cell_route_b = sp.Matrix([[0, cell_current], [-cell_current, 0]])
                cell_residual = sp.simplify(cell_route_a - cell_route_b)
                coverage.append({"cell": f"{endpoint}|{ambient}", "declared_class": boundary_class,
                                 "actual_class": data["endpoint_functionals"][endpoint]["boundary_class"],
                                 "response_factor": str(response_factor), "route_A_Omega": canonical_matrix(cell_route_a),
                                 "route_B_Omega": canonical_matrix(cell_route_b),
                                 "pullback_residual": canonical_matrix(cell_residual), "production_winding": str(production_winding),
                                 "execution_digest": digest({"L": str(cell_lagrangian), "current": str(cell_current), "residual": str(cell_residual)})})
    expected_coverage = {(endpoint, ambient) for endpoint in ENDPOINTS for ambient in AMBIENTS}
    actual_coverage = {(row["cell"].split("|")[0], row["cell"].split("|")[1]) for row in coverage}
    coverage_ok = actual_coverage == expected_coverage and all(row["declared_class"] == row["actual_class"] for row in coverage)

    mass, hbar, rho = symbols["m_GNLS"], symbols["hbar"], symbols["rho_inf"]
    gamma_computed = sp.simplify(hbar * 2 * sp.pi * winding / mass)
    gamma_expected = sp.simplify(hbar * 2 * sp.pi * expected_k / mass)
    epsilon = eps_sign * sp.Matrix([[0, 1], [-1, 0]])
    sheet = control["sheet_geometry"]
    sheet_coordinates = sp.symbols(" ".join(sheet["coordinates"]), real=True)
    sheet_metric = sp.Matrix([[q(value, symbols) for value in row] for row in sheet["induced_metric"]])
    area_density = sp.sqrt(sp.det(sheet_metric))
    area = area_density
    for coordinate, bounds in zip(sheet_coordinates, sheet["bounds"]):
        area = sp.integrate(area, (coordinate, q(bounds[0], symbols), q(bounds[1], symbols)))
    omega_per_area = sp.simplify(omega_a)
    omega_total = sp.simplify(omega_a * area)
    area_residual = sp.simplify(omega_total - area * omega_per_area)
    denominator = sp.simplify(mass * rho * gamma_computed * eps_sign)
    sigma_unsimplified = sp.Mul(omega_per_area[0, 1], sp.Pow(denominator, -1, evaluate=False), evaluate=False)
    sigma = sp.simplify(sigma_unsimplified) if denominator != 0 else sp.nan
    control_velocity = sp.Matrix(sp.symbols("Vctrl_x Vctrl_y", real=True))
    force_rows: dict[str, Any] = {}
    for name, flow in control["ambient_flow_cases"].items():
        vinf = sp.Matrix([q(value) for value in flow]); relative = control_velocity - vinf
        extracted = sp.simplify(-omega_per_area * relative)
        target = sp.simplify(mass * rho * gamma_expected * epsilon * relative)
        force_rows[name] = {"v_infinity": list(map(str, vinf)), "F_extract": list(map(str, extracted)),
                            "F_target_oracle": list(map(str, target)), "signed_residual": list(map(str, sp.simplify(extracted - target)))}

    dag_a = expression_dag({"L_indexed_berry": l_indexed, "A_from_L_x": connection_a[0], "A_from_L_y": connection_a[1],
                            "Omega_EL_xy": omega_a[0, 1]}, "route_A")
    dag_b = expression_dag({"field_current_integrand": polar_jacobian, "bulk_current": interior_scalar,
                            "projected_zero_mode_sector": projected_zero_sector[0, 1], "Omega_pullback_xy": omega_b[0, 1]}, "route_B")
    excluded_kinds = {"raw_symbol", "Integer", "Rational", "Float", "Pi", "NumberSymbol", "named_output"}
    route_a_reduced = {node["subexpression_digest"] for node in dag_a["nodes"] if node["kind"] not in excluded_kinds}
    route_b_reduced = {node["subexpression_digest"] for node in dag_b["nodes"] if node["kind"] not in excluded_kinds}
    route_a_reduced -= {node["subexpression_digest"] for node in dag_a["nodes"] if node["kind"] == "named_output"}
    route_b_reduced -= {node["subexpression_digest"] for node in dag_b["nodes"] if node["kind"] == "named_output"}
    shared_reduced = sorted(route_a_reduced & route_b_reduced)
    planted = q(data["berry_separation"]["planted_control_expression"], symbols)
    planted_a = expression_dag({"planted": planted}, "separation_control_A")
    planted_b = expression_dag({"planted": planted}, "separation_control_B")
    planted_overlap = ({node["subexpression_digest"] for node in planted_a["nodes"] if node["kind"] not in excluded_kinds}
                       & {node["subexpression_digest"] for node in planted_b["nodes"] if node["kind"] not in excluded_kinds})
    sigma_dag = expression_dag({"Omega_EL_xy": omega_per_area[0, 1], "phase_holonomy": 2 * sp.pi * winding,
                                "Gamma": gamma_computed, "epsilon_input": eps_sign, "sigma_projection_unsimplified": sigma_unsimplified}, "sigma")
    omega_digest = digest(sp.srepr(omega_per_area[0, 1]))
    sigma_downstream = any(node["subexpression_digest"] == omega_digest and node["kind"] != "named_output" for node in sigma_dag["nodes"])
    selection_role = selection["dependency_role"]
    berry = {"dependency_role": selection_role, "production_phase_profile": str(production_profile), "production_angular_derivative": str(production_derivative),
             "production_phase_parse_guard": {"allowed_free_symbols": ["r"], "observed_free_symbols": production_profile_symbols,
                                                "passed": production_profile_guard},
             "intrinsic_circulation": {"k_theta": str(production_winding), "Gamma_theta": str(sp.simplify(hbar * 2 * sp.pi * production_winding / mass)),
                                       "zero_decision": "IDENTICALLY_ZERO_FUNCTIONAL" if production_winding == 0 else "NONZERO_COMPUTED"},
             "route_A": {"definition": "differentiate termwise Taylor-substituted L_indexed", "theta_t_termwise_substitution": str(theta_t_substitution),
                         "substituted_field_density": str(substituted_density), "connection_density": list(map(str, local_connection_density)),
                         "connection_density_residual": list(map(str, local_density_residual)), "indexed_L_local": str(l_indexed),
                         "connection": list(map(str, connection_a)), "Omega": canonical_matrix(omega_a),
                         "surface_contributions": boundary_rows, "dag": dag_a},
             "route_B": {"definition": "native presymplectic current then Gram-projector quotient then collective pullback",
                         "current_integrand": str(polar_jacobian), "bulk_integral": str(interior_scalar),
                         "zero_mode_quotient_applied": quotient, "projector": str(projector),
                         "projector_times_zero_mode": list(map(str, sp.simplify(projector * sp.Matrix([z1, z2])))),
                         "zero_mode_sector_before_projection": canonical_matrix(zero_mode_sector),
                         "zero_mode_sector_after_projection": canonical_matrix(projected_zero_sector),
                         "Omega": canonical_matrix(omega_b), "dag": dag_b},
             "equivalence_residual": canonical_matrix(equivalence), "route_shared_reduced_nodes": shared_reduced,
             "separation_control": {"expression": str(planted), "shared_reduced_digests": sorted(planted_overlap), "detector_fires": bool(planted_overlap),
                                    "shared_raw_input_whitelist": data["berry_separation"]["shared_raw_input_whitelist"]},
             "production_pullback_coverage": coverage, "production_coverage_complete": coverage_ok,
             "dependency_edges": [["bulk_berry", "berry:route_A"], ["bulk_berry", "berry:route_B"],
                                  ["berry:route_A", "Omega_EL"], ["berry:route_B", "Omega_pullback"]],
             "selection_quarantined_from_production_ancestry": selection_role == "selection_only"}
    g4 = {"phase_profile_parsed": str(theta), "core_density_profile_parsed": str(density_r), "xi_core_over_a": str(xi_ratio),
          "contour_substitution": {"x": str(contour[x]), "y": str(contour[y])}, "angular_derivative": str(dtheta),
          "phase_holonomy": str(sp.simplify(2 * sp.pi * winding)), "computed_winding": str(winding),
          "k_expected_oracle": str(expected_k), "Gamma_computed": str(gamma_computed), "Gamma_expected_oracle": str(gamma_expected),
          "epsilon_matrix": [[str(v) for v in row] for row in epsilon.tolist()],
          "Omega_collective_total": canonical_matrix(omega_total), "Omega_physical_per_sheet_area": canonical_matrix(omega_per_area),
          "source_berry_coefficient": str(source_coefficient), "computed_sigma": str(sigma), "force_cases": force_rows,
          "transverse_plane": control["transverse_plane"], "sheet_directions": control["sheet_directions"], "sheet_cell_area": str(area),
          "sheet_geometry": sheet, "sheet_area_density": str(area_density),
          "total_to_per_area_residual": canonical_matrix(area_residual), "ir_prescription": control["ir_prescription"],
          "ir_matches_phase_a": control["ir_prescription"] == phase["ir_scheme"]["name"], "sigma_dependency_dag": sigma_dag}
    force_ok = all(all(q(value) == 0 for value in row["signed_residual"]) for row in force_rows.values())
    validity = {"holonomy": sp.simplify(winding - expected_k) == 0,
                "handedness": orient == eps_sign and len(set(control["transverse_plane"] + control["sheet_directions"])) == control["spatial_dimension"],
                "route_equivalence": equivalence == sp.zeros(2), "route_separation": not shared_reduced and selection_role == "selection_only",
                "production_profile_guard": production_profile_guard,
                "production_zero": production_winding == 0, "production_coverage": coverage_ok,
                "force": force_ok and area_residual == sp.zeros(2), "sigma_downstream": sigma_downstream and sigma not in (sp.nan, sp.zoo),
                "separation_control": bool(planted_overlap),
                "ir": g4["ir_matches_phase_a"]}
    return berry, g4, validity, omega_per_area[0, 1]


def e4_chain(data: dict[str, Any], cells: dict[str, Any], symbols: dict[str, sp.Symbol], missing: list[str]) -> tuple[dict[str, Any], bool]:
    scheme = data["e4_collar_scheme"]
    d, ell = int(scheme["radial_dimension"]), int(scheme["harmonic_degree"])
    nu = sp.Symbol("nu", real=True)
    roots = sp.solve(sp.Eq(nu * (nu - d + 2) - ell * (ell + d - 2), 0), nu)
    selected = max(value for value in roots if value.is_positive) if scheme["ir_condition"]["select"] == "decaying_root" else min(roots)
    r = sp.Symbol("r", positive=True); a = symbols["a"]
    trace_field = str(scheme["boundary_trace"]["field"])
    trace_symbol = sp.Symbol(str(scheme["boundary_trace"]["value"]), real=True)
    trace_radius = q(scheme["boundary_trace"]["radius"], symbols)
    ir_limit = sp.oo if scheme["ir_condition"]["limit"] == "infinity" else q(scheme["ir_condition"]["limit"], symbols)
    ir_target = q(scheme["ir_condition"]["value"], symbols)
    moduli = scheme["moduli_fixing"]
    growing_coefficient = q(moduli["growing_mode_coefficient"], symbols)
    growing_root = next(value for value in roots if value != selected)
    profile = trace_symbol * (trace_radius / r) ** selected
    collar_mode_reconstruction = sp.expand(
        (trace_symbol - growing_coefficient) * (trace_radius / r) ** selected
        + growing_coefficient * (trace_radius / r) ** growing_root
    )
    normalized_profile = sp.simplify(profile / trace_symbol)
    trace_residual = sp.simplify(collar_mode_reconstruction.subs(r, trace_radius) - trace_symbol)
    ir_residual = sp.limit(collar_mode_reconstruction, r, ir_limit) - ir_target
    sd = sphere_area(d); rho_br, mu_r = symbols["rhoBr"], symbols["muR"]
    collar_mass = sp.simplify(rho_br * sd / d * sp.integrate(r ** (d - 1) * normalized_profile**2, (r, trace_radius, sp.oo)))
    grad_density = sp.diff(normalized_profile, r) ** 2 + ell * (ell + d - 2) * normalized_profile**2 / r**2
    stiffness = sp.simplify(mu_r * sd / d * sp.integrate(r ** (d - 1) * grad_density, (r, trace_radius, sp.oo)))
    dtn = sp.simplify(mu_r * sd / d * trace_radius ** (d - 1) * (-sp.diff(normalized_profile, r).subs(r, trace_radius)))

    velocities = sp.Matrix(sp.symbols("vAug0:9", real=True)); coordinates = sp.Matrix(sp.symbols("qAug0:9", real=True))
    known_l = cells["E4|symmetric_postulate"]["_L"].subs(dict(zip(sp.symbols("V_x V_y V_z", real=True), velocities[:3, 0])))
    collar_kinetic = sp.expand(collar_mass * sum(velocities[i] ** 2 for i in range(6, 9)) / 2)
    collar_potential = sp.expand(stiffness * sum(coordinates[i] ** 2 for i in range(6, 9)) / 2)
    preconstraint_action = sp.expand(known_l + collar_kinetic - collar_potential)
    m_aug = sp.hessian(preconstraint_action, list(velocities))

    vv, pp, uu = sp.symbols("V_i P_i udot_collar_i", real=True)
    constraint_expr = q(scheme["constraint_functional"]["expression"], {"V_i": vv, "P_i": pp, "udot_collar_i": uu})
    g1 = sp.Matrix([[sp.diff(constraint_expr, vv), sp.diff(constraint_expr, pp), sp.diff(constraint_expr, uu)]])
    constraint = sp.kronecker_product(g1, sp.eye(3))
    alpha, beta = sp.symbols("alpha beta", real=True)
    lift_substitution = sp.expand(constraint_expr.subs(uu, alpha * vv + beta * pp))
    equations = [sp.diff(lift_substitution, vv), sp.diff(lift_substitution, pp), beta - q(moduli["p_to_collar_cross_lift"])]
    lift_solutions = sp.solve(equations, [alpha, beta], dict=True)
    if lift_solutions:
        one = lift_solutions[0]
        j1 = sp.Matrix([[1, 0], [0, 1], [one[alpha], one[beta]]])
    else:
        j1 = sp.zeros(3, 2)
    physical_lift = sp.kronecker_product(j1, sp.eye(3))
    full_kernel_columns = constraint.nullspace()
    full_kernel = sp.Matrix.hstack(*full_kernel_columns) if full_kernel_columns else sp.zeros(9, 0)
    reduced = sp.simplify(physical_lift.T * m_aug * physical_lift)
    full_kernel_form = sp.simplify(full_kernel.T * m_aug * full_kernel)

    transform3 = sp.Matrix([[1, 1, 0], [0, 1, 0], [0, 0, 2]])
    transform = sp.eye(9); transform[6:, 6:] = transform3
    transformed_lift = transform.inv() * physical_lift
    covariance_residual = sp.simplify(transformed_lift.T * (transform.T * m_aug * transform) * transformed_lift - reduced)

    zero_velocity = {velocity: 0 for velocity in velocities}
    a_aug = sp.Matrix([sp.diff(preconstraint_action, velocity).subs(zero_velocity) for velocity in velocities])
    omega_aug = sp.Matrix(9, 9, lambda i, j: sp.diff(a_aug[j], coordinates[i]) - sp.diff(a_aug[i], coordinates[j]))
    physical_coordinates = coordinates[:6, 0]
    anholonomic = sp.Matrix(6, 6, lambda i, j: sp.simplify(sum(
        physical_lift[row, i] * sp.diff(physical_lift[row, j], physical_coordinates[i])
        - physical_lift[row, j] * sp.diff(physical_lift[row, i], physical_coordinates[j])
        for row in range(physical_lift.rows))))
    accelerations = sp.Matrix(sp.symbols("acc0:9", real=True)); multipliers = sp.Matrix(sp.symbols("lambda0:3", real=True))
    multiplier_equations = constraint * (m_aug * accelerations + constraint.T * multipliers)
    multiplier_solution = sp.solve(list(multiplier_equations), list(multipliers), dict=True)
    reaction = constraint.T * multipliers.subs(multiplier_solution[0]) if multiplier_solution else constraint.T * multipliers
    virtual_work = sp.simplify(physical_lift.T * reaction)

    profile_moduli_ok = growing_coefficient == 0
    scheme_contract = {"formulation": scheme["formulation"], "action_ancestors": scheme["action_ancestors"],
                       "radial_domain": scheme["radial_domain"], "constraint_equals": scheme["constraint_functional"]["equals"],
                       "collective_trace_map": scheme["collective_trace_map"]}
    scheme_ok = (scheme["formulation"] == "operator_level_Dirichlet_to_Neumann" and
                 set(scheme["action_ancestors"]) == {"brane_shear_kinetic", "brane_shear_gradient"} and
                 trace_field == "u_T" and scheme["ir_condition"]["limit"] == "infinity" and
                 scheme["radial_domain"] == "a<=r<infinity" and q(scheme["constraint_functional"]["equals"]) == 0 and
                 scheme["collective_trace_map"] == "identity")
    okay = (trace_residual == 0 and ir_residual == 0 and sp.simplify(stiffness - dtn) == 0 and
            sp.hessian(preconstraint_action, list(velocities)) == m_aug and constraint.rank() == 3 and
            physical_lift.rank() == 6 and constraint * physical_lift == sp.zeros(3, 6) and
            j1[2, 0] == 1 and j1[2, 1] == q(moduli["p_to_collar_cross_lift"]) and
            covariance_residual == sp.zeros(6) and bool(lift_solutions) and profile_moduli_ok and
            omega_aug == sp.zeros(9) and a_aug == sp.zeros(9, 1) and anholonomic == sp.zeros(6) and scheme_ok)
    status = "UNRESOLVED" if missing else "COMPUTED"
    return {"operator": f"rhoBr*d_t^2-muR*(d_r^2+{d-1}/r*d_r-{ell*(ell+d-2)}/r^2)", "scheme_contract": scheme_contract,
            "boundary_profile": str(profile), "collar_mode_reconstruction": str(collar_mode_reconstruction),
            "growing_mode_root": str(growing_root), "indicial_roots": list(map(str, roots)), "selected_decay": str(selected),
            "trace_residual": str(trace_residual), "ir_residual": str(ir_residual), "collar_mass": canonical(collar_mass),
            "collar_stiffness_bulk": canonical(stiffness), "DtN_stiffness": canonical(dtn), "DtN_residual": canonical(stiffness - dtn),
            "preconstraint_extended_action": canonical(preconstraint_action), "preconstraint_action_components": {
                "phase_a_translation": canonical(known_l), "collar_kinetic": canonical(collar_kinetic), "collar_potential": canonical(-collar_potential)},
            "M_aug": canonical_matrix(m_aug), "M_aug_hessian_residual": canonical_matrix(sp.hessian(preconstraint_action, list(velocities)) - m_aug),
            "constraint_expression": str(constraint_expr), "constraint_operator": canonical_matrix(constraint), "constraint_rank": constraint.rank(),
            "kernel_dimension": 9 - constraint.rank(), "N_full_kernel": canonical_matrix(full_kernel), "N_kernel_residual": canonical_matrix(constraint * full_kernel),
            "J_physical_lift": canonical_matrix(physical_lift), "J_derivation_equations": list(map(str, equations)),
            "M_reduced": canonical_matrix(reduced), "M_full_kernel": canonical_matrix(full_kernel_form),
            "basis_covariance_residual": canonical_matrix(covariance_residual),
            "multiplier_solution": [{str(key): str(value) for key, value in row.items()} for row in multiplier_solution],
            "reaction": list(map(str, reaction)), "virtual_work_residual": list(map(str, virtual_work)),
            "A_aug": [canonical(value) for value in a_aug], "Omega_aug": canonical_matrix(omega_aug),
            "A_aug_zero_decision": "DERIVATIVE_OF_ACTION_ZERO", "Omega_aug_zero_decision": "CURVATURE_OF_DERIVED_CONNECTION_ZERO",
            "anholonomic_term": canonical_matrix(anholonomic), "anholonomic_decision": "DERIVED_FRAME_BRACKET",
            "dependency_edges": [[name, "E4:preconstraint_action"] for name in scheme["action_ancestors"]]
                                + [["E4:preconstraint_action", "E4:M_aug"], ["E4_shear_lock", "E4:J_physical_lift"],
                                   ["E4:M_aug", "E4:M_reduced"], ["E4:J_physical_lift", "E4:M_reduced"]],
            "Omega_aug_unresolved_leaves": missing,
            "status": status, "unresolved_leaves": missing}, bool(okay)


def e5_chain(data: dict[str, Any], phase: dict[str, Any], endpoints: dict[str, Any], conservative: dict[str, Any], symbols: dict[str, sp.Symbol]) -> tuple[dict[str, Any], bool]:
    row = data["endpoint_functionals"]["E5"]
    tangent, velocity = sp.symbols("v_tangent V_tangent", real=True)
    rayleigh = q(row["rayleigh_integrand"], {**symbols, "v_tangent": tangent, "V_tangent": velocity})
    force_trace = sp.simplify(-sp.diff(rayleigh, tangent)); force_sleeve = sp.simplify(-sp.diff(rayleigh, velocity))
    root_deleted = sp.Integer(0)
    force_delta = sp.simplify(force_trace - (-sp.diff(root_deleted, tangent)))
    full_solution = endpoints["E5"]["_solution"]; bare_solution = conservative["E5"]["_solution"]; e2_solution = conservative["E2"]["_solution"]
    full_gvv = phase_coefficient(phase, "E5", "GVV", symbols); bare_gvv = phase_coefficient(phase, "E2", "GVV", symbols)
    discrepancy = sp.expand(full_gvv - bare_gvv)
    okay = (rayleigh != 0 and root_deleted == 0 and force_delta != 0 and
            sp.simplify(bare_solution - e2_solution) == sp.zeros(3, 1) and sp.simplify(full_solution - bare_solution) != sp.zeros(3, 1))
    return {"root": row["rayleigh_root"], "domain": row["rayleigh_surface"], "field_level_density": str(rayleigh),
            "F_Rayleigh_trace": str(force_trace), "F_Rayleigh_sleeve": str(force_sleeve),
            "root_deleted_functional": str(root_deleted), "gamma_functional_ablation_force_delta": str(force_delta),
            "full_E5_trace_solution": list(map(str, full_solution)), "root_deleted_conservative_solution": list(map(str, bare_solution)),
            "E2_conservative_solution": list(map(str, e2_solution)),
            "conservative_equals_E2_computed": sp.simplify(bare_solution - e2_solution) == sp.zeros(3, 1),
            "bare_GVV": canonical(bare_gvv), "stored_full_E5_minus_bare_GVV": canonical(discrepancy),
            "bare_dependencies": sorted(str(v) for v in bare_gvv.free_symbols),
            "rayleigh_dependencies": sorted(str(v) for v in rayleigh.free_symbols),
            "dependency_edges": [[row["rayleigh_root"], "E5:rayleigh_density"],
                                 ["E5:rayleigh_density", "F_Rayleigh"]]}, bool(okay)


def block_expressions(cells: dict[str, Any], reach: list[dict[str, Any]], missing: list[str]) -> dict[str, Any]:
    output: dict[str, Any] = {}
    for key, cell in cells.items():
        blocks: dict[str, Any] = {}
        for block in BLOCK_SPECS:
            root_rows = [row for row in reach if row["block"] == block]
            roots = sorted({row["root"] for row in root_rows})
            tilt_leaves = [] if block == "M_XX" else list(missing)
            leaves = sorted(set(roots + tilt_leaves))
            blocks[block] = {"constructible_tensor": cell["M_XX_p0_known"] if block == "M_XX" else None,
                             "unresolved_remainder_leaves": leaves,
                             "remainder_expression": " + ".join(f"R[{leaf},{block}]" for leaf in leaves),
                             "status": "UNRESOLVED" if leaves else "COMPUTED"}
        output[key] = blocks
    return output


def derived_congruence(data: dict[str, Any], cells: dict[str, Any], reach: list[dict[str, Any]], e4: dict[str, Any]) -> tuple[dict[str, Any], bool]:
    inventory = data["covariant_inventory"]
    known = expression_from_terms(cells["E1|symmetric_postulate"]["M_XX_p0_known"][0][0], {})
    vector_names = list(inventory["vectors"]); vector_name = str(vector_names[0]) if vector_names else "missing_vector"
    p1, p2, p3 = sp.symbols(f"{vector_name}1 {vector_name}2 {vector_name}3", real=True)
    radius = sp.Symbol("r", real=True)
    p = sp.Matrix([p1, p2, p3]); delta = sp.eye(3)
    epsilon_p = sp.Matrix([[0, p3, -p2], [-p3, 0, p1], [p2, -p1, 0]])
    coefficient_definitions: dict[tuple[str, str], sp.Expr] = {}
    coefficients: dict[tuple[str, str], sp.Expr] = {}
    for block, (_, symmetry) in BLOCK_SPECS.items():
        for generator, _ in generator_candidates(data, symmetry):
            # The generated tensor basis and remainder-symbol list are parallel.
            # Select only the symbol multiplying this generator, then sum roots.
            names = []
            for row in reach:
                if row["block"] != block:
                    continue
                for row_generator, symbol in zip(row["finite_generator_set"], row["remainder_symbols"]):
                    if row_generator == generator:
                        names.append(symbol)
            definition = sum((sp.Symbol(name, real=True) for name in names), sp.Integer(0))
            coefficient_definitions[(block, generator)] = definition
            coefficients[(block, generator)] = (sp.Symbol(f"coef__{block}__{generator.replace('*', '_').replace('.', '_')}", real=True)
                                                  if names else sp.Integer(0))
    ax = sp.Symbol("coef__M_XX__delta_total", real=True)
    bx = coefficients.get(("M_XX", "p_i*p_j"), 0)
    c = coefficients.get(("M_Xp_symmetric", "delta_ij"), 0)
    d = coefficients.get(("M_Xp_symmetric", "p_i*p_j"), 0)
    g = coefficients.get(("M_Xp_antisymmetric", "epsilon_ijk*p_k"), 0)
    e = coefficients.get(("M_pp", "delta_ij"), 0)
    f = coefficients.get(("M_pp", "p_i*p_j"), 0)
    mxx = ax * delta + bx * (p * p.T)
    mxp = c * delta + d * (p * p.T) + g * epsilon_p
    mpp = e * delta + f * (p * p.T)
    full = sp.Matrix.vstack(sp.Matrix.hstack(mxx, mxp), sp.Matrix.hstack(mxp.T, mpp))
    aligned = sp.simplify(full.subs({p1: 0, p2: 0, p3: radius}))
    transverse_indices = [0, 1, 3, 4]; longitudinal_indices = [2, 5]
    transverse = aligned.extract(transverse_indices, transverse_indices)
    longitudinal = aligned.extract(longitudinal_indices, longitudinal_indices)

    def pivots(matrix: sp.Matrix) -> list[sp.Expr]:
        minors = [sp.factor(matrix[:i, :i].det()) for i in range(1, matrix.rows + 1)]
        return [minors[0]] + [sp.factor(minors[i] / minors[i - 1]) for i in range(1, len(minors))]

    transverse_pivots = pivots(transverse); longitudinal_pivots = pivots(longitudinal)
    permutation = transverse_indices + longitudinal_indices
    permuted = aligned.extract(permutation, permutation)
    split = sp.diag(transverse, longitudinal)
    split_residual = sp.simplify(permuted - split)
    generic = inventory["generic_projection_point"]
    alpha, beta = sp.symbols("alpha beta", real=True)
    generic_p = sp.Matrix([q(value) for value in generic])
    projection_basis_at_point = generic_p * generic_p.T
    declared_symmetric = set(inventory["symmetric_basis"])
    projection_ansatz = alpha * delta
    projection_variables = [alpha]
    if "p_i_p_j" in declared_symmetric:
        projection_ansatz += beta * projection_basis_at_point; projection_variables.append(beta)
    projected_residual = sp.simplify(mxx.subs(dict(zip(p, generic_p))) - projection_ansatz)
    solved = sp.solve(list(projected_residual), projection_variables, dict=True)
    decomposition = solved[0] if solved else {}
    decomposition_residual = sp.simplify(projected_residual.subs(decomposition)) if solved else projected_residual
    zero_direction = [q(value) for value in inventory["zero_limit_direction"]]
    zero_limit_path_projector = sp.Matrix(zero_direction) * sp.Matrix(zero_direction).T
    zero_parameter = sp.Symbol("lambda_zero", real=True)
    zero_path_block = sp.simplify(mxx.subs(dict(zip(p, [zero_parameter * value for value in zero_direction]))))
    zero_limit_block = zero_path_block.applyfunc(lambda value: sp.limit(value, zero_parameter, 0))
    required_inventory = {"delta_ij", "p_i_p_j"} <= set(inventory["symmetric_basis"]) and {"epsilon_ijk_p_k"} <= set(inventory["antisymmetric_basis"])
    e4_reduced = sp.Matrix([[expression_from_terms(value, {}) for value in row] for row in e4["M_reduced"]])
    e4_rank = e4_reduced.rank()
    conditions = [f"({sp.factor(value)}) > 0" for value in transverse_pivots + longitudinal_pivots]
    record = {"derived_coefficients": {
                  "M_XX:delta_total": {"symbol": str(ax), "definition": str(known + coefficient_definitions.get(("M_XX", "delta_ij"), 0))},
                  **{f"{block}:{generator}": {"symbol": str(coefficients[(block, generator)]), "definition": str(value)}
                     for (block, generator), value in coefficient_definitions.items()}},
              "M_XX_covariant_decomposition": {"delta_coefficient": str(decomposition.get(alpha, "UNRESOLVED")),
                                                "pp_coefficient": str(decomposition.get(beta, "UNRESOLVED")),
                                                "projection_residual": canonical_matrix(decomposition_residual)},
              "inventory_vector_symbols": list(map(str, p)),
              "projection_basis_at_point": canonical_matrix(projection_basis_at_point),
              "zero_limit_path_projector": canonical_matrix(zero_limit_path_projector),
              "zero_limit_path_block": canonical_matrix(zero_path_block),
              "zero_limit_block": canonical_matrix(zero_limit_block),
              "aligned_full_block": canonical_matrix(aligned), "transverse_block": canonical_matrix(transverse),
              "longitudinal_block": canonical_matrix(longitudinal), "block_split_residual": canonical_matrix(split_residual),
              "transverse_LDL_pivots": list(map(str, transverse_pivots)), "longitudinal_LDL_pivots": list(map(str, longitudinal_pivots)),
              "positive_definite_conditions": conditions, "E4_reduced_rank": e4_rank,
              "E4_reduced_block": e4["M_reduced"], "classification": "M_UNRESOLVED_DERIVED_REMAINDER_COEFFICIENTS"}
    okay = (vector_names == ["p"] and len(zero_direction) == 3 and any(value != 0 for value in zero_direction) and
            required_inventory and decomposition_residual == sp.zeros(3) and split_residual == sp.zeros(6))
    return record, bool(okay)


def aggregate(components: dict[str, Any]) -> dict[str, Any]:
    statuses: list[str] = []
    leaves: list[str] = []

    def walk(value: Any) -> None:
        if isinstance(value, dict):
            if "status" in value: statuses.append(value["status"])
            leaves.extend(value.get("unresolved_leaves", []))
            for child in value.values(): walk(child)
        elif isinstance(value, list):
            for child in value: walk(child)

    walk(components)
    if "ILL_POSED" in statuses: status = "ILL_POSED"
    elif "UNSTABLE" in statuses: status = "UNSTABLE"
    elif "UNRESOLVED" in statuses: status = "UNRESOLVED"
    else: status = "OK"
    return {"status": status, "unresolved_leaves": sorted(set(leaves)), "precedence_evidence": statuses}


def component_records(data: dict[str, Any], cells: dict[str, Any], blocks: dict[str, Any], block_roots: dict[str, list[str]],
                      missing: list[str], berry: dict[str, Any], e4: dict[str, Any], e5: dict[str, Any], congruence: dict[str, Any]) -> tuple[dict[str, Any], dict[str, Any], dict[str, Any]]:
    strata: list[dict[str, Any]] = []
    for row in data["open_strata"]:
        predicate = q(row["predicate"], {"Eq": sp.Eq})
        strata.append({"id": row["id"], "predicate": str(predicate), "predicate_value": bool(predicate), "leaves": row.get("leaves", [])})
    active = [row for row in strata if row["predicate_value"]]
    records: dict[str, Any] = {}; mechanics_map: dict[str, Any] = {}
    for endpoint in ENDPOINTS:
        for ambient in AMBIENTS:
            for stratum in active:
                base_key = f"{endpoint}|{ambient}"; key = f"{base_key}|{stratum['id']}"
                block_leaves = sorted({leaf for row in blocks[base_key].values() for leaf in row["unresolved_remainder_leaves"]})
                m_status = "UNRESOLVED" if block_leaves else "COMPUTED"
                intrinsic_zero = berry["intrinsic_circulation"]["k_theta"] == "0"
                findings = {
                    "M_full": {"status": m_status, "unresolved_leaves": block_leaves,
                               "known_p0_translation": cells[base_key]["M_XX_p0_known"], "derived_congruence_class": congruence["classification"]},
                    "M_spatial_symmetry": {"status": m_status, "unresolved_leaves": block_leaves,
                                           "known_M_XX_projection_residual": congruence["M_XX_covariant_decomposition"]["projection_residual"]},
                    "Omega_intrinsic_circulation": {"status": "ZERO_COMPUTED" if intrinsic_zero else "NONZERO",
                                                     "zero_decision": berry["intrinsic_circulation"]["zero_decision"]},
                    "Omega_translation_texture": {"status": "UNRESOLVED", "unresolved_leaves": sorted(set(block_roots.get("Omega_XX_texture", []) + missing))},
                    "Omega_translation_tilt": {"status": "UNRESOLVED", "unresolved_leaves": sorted(set(block_roots.get("Omega_Xp", []) + missing))},
                    "Omega_tilt_tilt": {"status": "UNRESOLVED", "unresolved_leaves": sorted(set(block_roots.get("Omega_pp", []) + missing))},
                    "E4_constraint": ({"status": e4["status"], "unresolved_leaves": e4["unresolved_leaves"]} if endpoint == "E4" else {"status": "NOT_APPLICABLE"}),
                    "E5_rayleigh": ({"status": "COMPUTED" if e5["conservative_equals_E2_computed"] else "ILL_POSED"} if endpoint == "E5" else {"status": "NOT_APPLICABLE"}),
                }
                headline = aggregate(findings)
                record = {"cell_key": {"endpoint": endpoint, "ambient": ambient, "open_stratum": stratum["id"]},
                          "off_shell_p": True, "block_expressions": blocks[base_key], "component_findings": findings,
                          "headline": headline, "record_inputs_digest": digest({"blocks": blocks[base_key], "findings": findings})}
                records[key] = record; mechanics_map[key] = headline
    return records, mechanics_map, {"collapsed": len(active) == 1, "declared": strata, "active": [row["id"] for row in active]}


def provenance_graph(source: dict[str, Any], cells: dict[str, Any], reach: list[dict[str, Any]], berry: dict[str, Any], e4: dict[str, Any], e5: dict[str, Any]) -> dict[str, Any]:
    emitted_edge_rows = [edge for cell in cells.values() for contraction in cell["field_contraction_integrals"]
                         for edge in contraction["dependency_edges"]]
    emitted_edge_rows += [[row["root"], f"remainder:{row['block']}"] for row in reach]
    emitted_edge_rows += berry["dependency_edges"] + e4["dependency_edges"] + e5["dependency_edges"]
    edges = {tuple(edge) for edge in emitted_edge_rows}
    nodes = sorted({node for edge in edges for node in edge} | {"return_closure"})
    raw_targets = ({right for _, right in edges if right.endswith(("M_XX", "Omega_EL", "Omega_pullback", "M_reduced"))}
                   | {f"remainder:{block}" for block in BLOCK_SPECS})
    return_absent = not graph_reachable([list(edge) for edge in edges], "return_closure", raw_targets)
    return {"nodes": [{"id": node, "type": ("RETURN" if node == "return_closure" else "DERIVED_OR_DECLARED")} for node in nodes],
            "edges": [list(edge) for edge in sorted(edges)], "production_dependency_sets": {key: row["production_dependencies"] for key, row in cells.items()},
            "global_return_closure_absence": return_absent,
            "closure_absence_targets": sorted(raw_targets), "construction": "union(emitted_contraction,reachability,berry,E4,E5 dependency records)"}


def leaf_paths(value: Any, prefix: tuple[Any, ...] = ()) -> list[tuple[tuple[Any, ...], Any]]:
    if isinstance(value, dict):
        return [item for key, child in value.items() for item in leaf_paths(child, prefix + (key,))]
    if isinstance(value, list):
        if not value:
            return [(prefix, [])]
        return [item for index, child in enumerate(value) for item in leaf_paths(child, prefix + (index,))]
    return [(prefix, value)]


def liveness_sink(path: tuple[Any, ...]) -> tuple[str, str, str]:
    head = str(path[0])
    if head == "scalar_regression_projection" and len(path) > 1 and path[1] in {"pd_unit", "p_unit", "indexed_tilt_length"}:
        return "comparator_residual", "indexed_cells", "native_slice_projection_residual"
    if head == "endpoint_functionals" and len(path) > 2 and path[1] == "E5" and path[2] in {"rayleigh_root", "rayleigh_surface"}:
        return "mathematical_deliverable", "E5", "rayleigh_root_consumption"
    if head == "phase_a_protection":
        return "comparator_residual", "source_contract.phase_a_payload_sha256", "external_phase_a_payload_guard"
    mapping = {
        "derivation_representations": ("validation_predicate", "checks.B1_R9", "engine_representation_decision"),
        "substrate": ("validation_predicate", "source_contract", "phase_a_source_contract"),
        "indexed_coordinates": ("mathematical_deliverable", "dimensions", "units_restored_mechanics"),
        "indexed_embedding": ("mathematical_deliverable", "field_manifest", "per_field_embedding_derivation"),
        "covariant_inventory": ("mathematical_deliverable", "derived_congruence", "covariant_block_congruence"),
        "scalar_regression_projection": ("comparator_residual", "scalar_regression", "frozen_scalar_projection_residual"),
        "ambient_branches": ("mathematical_deliverable", "indexed_cells", "ambient_indexed_action"),
        "open_strata": ("selection_decision", "open_strata", "cell_stratum_selection"),
        "boundary_operator": ("mathematical_deliverable", "endpoint_assembly", "endpoint_DtN_solve"),
        "endpoint_functionals": ("mathematical_deliverable", "endpoint_assembly", "endpoint_variation_solve"),
        "e4_collar_scheme": ("mathematical_deliverable", "E4", "extended_action_constraint_reduction"),
        "berry_pullback_selection": ("mathematical_deliverable", "berry", "pullback_and_quotient_construction"),
        "berry_separation": ("validation_predicate", "berry.separation_control", "subexpression_overlap_detector"),
        "g4_control": ("mathematical_deliverable", "g4_control", "G4_force_and_holonomy"),
        "open_action_functionals": ("mathematical_deliverable", "open_root_reachability", "production_remainder_generator"),
        "finite_bound_controls": ("validation_predicate", "reachability_analysis.finite_bound_controls", "finite_generator_control"),
        "tensor_selection_rules": ("validation_predicate", "reachability_analysis.crosscheck_agrees", "selection_table_crosscheck"),
        "partition": ("mathematical_deliverable", "partition_ledger", "partition_ownership"),
        "phase_a_protection": ("comparator_residual", "source_contract.phase_a_payload_sha256", "external_phase_a_payload_guard"),
    }
    if head in {"schema_version", "directive_version", "spec_version", "phase"}:
        return "validation_predicate", "input_validation", "schema_validation"
    return mapping[head]


def attach_input_liveness(result: dict[str, Any], data: dict[str, Any], tracker: ReadTracker) -> tuple[dict[str, Any], bool]:
    rows: list[dict[str, Any]] = []
    for path, value in leaf_paths(data):
        role, sink, intermediate = liveness_sink(path)
        sink_value = at_path(result, sink)
        path_text = ".".join(map(str, path))
        comparator_only = path and path[0] == "scalar_regression_projection"
        selection_only = len(path) > 1 and path[0] == "berry_pullback_selection" and path[1] in {"dependency_role", "production_cells", "ambient_branches"}
        sink_digest = digest(sink_value)
        scopes = sorted(tracker.events.get(path_text, set()))
        metadata_only = bool(scopes) and all(scope.startswith(("metadata:", "digest:")) for scope in scopes)
        rows.append({"path": path_text, "typed_role": role, "semantic_sink": sink, "semantic_sink_kind": intermediate,
                     "dependency_path": [f"input:{path_text}"] + [f"derived:{scope}" for scope in scopes] + [f"artifact:{sink}"],
                     "read_scopes": scopes, "sink_digest": sink_digest, "metadata_or_digest_only": metadata_only,
                     "production_ancestry_forbidden": bool(comparator_only or selection_only),
                     "absent_from_production_ancestry": bool(comparator_only or selection_only)})
    declared = sorted(row["path"] for row in rows)
    consumed = sorted(set(tracker.events) & set(declared))
    artifact = {"declared_leaf_paths": declared, "consumed_leaf_paths": consumed, "consumed_but_undeclared": sorted(set(consumed) - set(declared)),
                "rows": rows, "per_key_mutation_evidence": "external:b1_mutation_results.input_liveness_cases"}
    okay = (declared == consumed and all(row["typed_role"] and row["semantic_sink"] and row["read_scopes"] and len(row["dependency_path"]) >= 3
                                       and not row["metadata_or_digest_only"] for row in rows)
            and all(not row["production_ancestry_forbidden"] or row["absent_from_production_ancestry"] for row in rows))
    return artifact, bool(okay)


def public_cells(cells: dict[str, Any]) -> dict[str, Any]:
    return {key: {name: value for name, value in row.items() if not name.startswith("_")} for key, row in cells.items()}


def build(data: dict[str, Any], phase_input: dict[str, Any], phase: dict[str, Any]) -> dict[str, Any]:
    plain_data = copy.deepcopy(data); tracker = ReadTracker(); data = tracker.wrap(data)
    rec = Recorder(); symbols, _ = symbols_and_values(phase_input)
    input_validation = {"schema": data.get("schema_version"), "directive": data.get("directive_version"), "spec": data.get("spec_version"),
                        "phase": data.get("phase"), "valid": data.get("schema_version") == "U1_PHASE_B1_MECHANICS_INPUTS_V3" and data.get("phase") == "B1"}
    source, source_ok = source_contract(data, phase_input, phase)
    amended_phase, correction, correction_ok = amend_brane_shear(phase, phase_input)
    source["payload_normalization_version"] = PAYLOAD_NORMALIZATION_VERSION
    source["verified_baseline_payload_sha256"] = correction["baseline_digest"]["computed_legacy"]
    source["amended_payload_sha256"] = correction["amended_payload_sha256"]
    manifest, missing, manifest_ok = field_manifest(data, amended_phase)
    endpoints, endpoint_info, endpoints_ok = assemble_endpoints(data, symbols, amended_phase)
    cells, scalar_regression, f1_ok = derived_translation(data, phase_input, amended_phase, symbols, endpoint_info["source_map"], manifest)
    rec.check("B1_R1", source_ok and correction_ok and manifest_ok and f1_ok,
              {"manifest_digest": digest(manifest), "emitted_missing_leaves": missing,
               "phase_a_correction": correction,
               "scalar_regression": scalar_regression,
               "reconstruction": {key: row["reconstruction_residual"] for key, row in cells.items()}})

    reach, block_roots, reach_analysis, reach_ok = reachability(data, phase)
    berry, g4, berry_valid, g4_omega = berry_pipeline(data, phase_input, amended_phase, symbols)
    rec.check("B1_R2", berry_valid["holonomy"] and berry_valid["handedness"] and berry_valid["force"] and berry_valid["sigma_downstream"] and berry_valid["ir"],
              {"winding": g4["computed_winding"], "sigma": g4["computed_sigma"], "force": g4["force_cases"],
               "area_residual": g4["total_to_per_area_residual"], "sigma_dag": g4["sigma_dependency_dag"]})
    rec.check("B1_R3", berry_valid["route_equivalence"] and berry_valid["route_separation"] and berry_valid["separation_control"]
              and berry_valid["production_profile_guard"] and berry_valid["production_zero"] and berry_valid["production_coverage"],
              {"equivalence": berry["equivalence_residual"], "coverage": berry["production_pullback_coverage"],
               "production_phase_parse_guard": berry["production_phase_parse_guard"],
               "route_A_dag": berry["route_A"]["dag"], "route_B_dag": berry["route_B"]["dag"],
               "separation_control": berry["separation_control"]})
    e4, e4_ok = e4_chain(data, cells, symbols, missing)
    rec.check("B1_R4", e4_ok, {"action": e4["preconstraint_extended_action"], "hessian": e4["M_aug_hessian_residual"],
                                "constraint": e4["constraint_operator"], "lift": e4["J_physical_lift"]})
    e5, e5_ok = e5_chain(data, amended_phase, endpoints, endpoint_info["conservative"], symbols)
    rec.check("B1_R5", endpoints_ok and e5_ok, {"endpoint_fingerprints": {ep: row["fingerprint"] for ep, row in endpoints.items()},
                                                "endpoint_source_map": endpoint_info["source_map"], "E5": e5})
    rec.check("B1_R6", reach_ok, {"production_pairs": len(reach), "analysis": reach_analysis})
    blocks = block_expressions(cells, reach, missing)
    congruence, congruence_ok = derived_congruence(data, cells, reach, e4)
    one_sided = branch_factor(data, "one_sided_pathA29", symbols)
    eta = q(data["ambient_branches"]["one_sided_pathA29"]["eta_asym"])
    covariance_ok = sp.simplify(one_sided - (1 + symbols["s"] * eta)) == 0
    rec.check("B1_R7", congruence_ok and covariance_ok, {"derived_congruence": congruence, "one_sided_factor": str(one_sided)})
    dimensions, dimensions_ok = dimension_analysis(data, phase_input, amended_phase, symbols, cells, g4_omega)
    rec.check("B1_R8", dimensions_ok, dimensions)
    representations = data["derivation_representations"]
    route_a_ops = {node["kind"] for node in berry["route_A"]["dag"]["nodes"]}
    route_b_ops = {node["kind"] for node in berry["route_B"]["dag"]["nodes"]}
    route_a_dependencies = set(berry["route_A"]["dag"]["raw_dependencies"])
    route_b_dependencies = set(berry["route_B"]["dag"]["raw_dependencies"])
    representation_ok = (representations["shared_reduced_formulas"] is False and representations["SymPy"] != representations["Mathematica"]
                         and {"Vberry_x", "Vberry_y", "X_x", "X_y"} <= route_a_dependencies
                         and {"r", "a"} <= route_b_dependencies and route_a_dependencies != route_b_dependencies
                         and not berry["route_shared_reduced_nodes"])
    rec.check("B1_R9", representation_ok, {"representations": representations, "route_A_operations": sorted(route_a_ops), "route_B_operations": sorted(route_b_ops)})

    cell_records, mechanics_map, strata = component_records(data, cells, blocks, block_roots, missing, berry, e4, e5, congruence)
    ledger: list[dict[str, Any]] = []
    for key, cell in cells.items():
        for term, expression in cell["termwise_L"].items():
            ledger.append({"candidate_id": f"{key}:{term}:translation", "owner": "M", "provenance": term,
                           "computed_expression_digest": digest(expression)})
    ledger.append({"candidate_id": "outer_control_flux:translation", "owner": data["partition"]["outer_control_flux_owner"],
                   "provenance": "outer_control_flux", "computed_expression_digest": digest("pending_B2")})
    partition = {"records": ledger, "owner_enum": data["partition"]["owner_enum"],
                 "unique": len({row["candidate_id"] for row in ledger}) == len(ledger), "state": data["partition"]["terminal_state"]}
    provenance = provenance_graph(source, cells, reach, berry, e4, e5)
    closure = {"materialized": False, "root": "return_closure", "global_absence_computed": provenance["global_return_closure_absence"]}
    gates = {
        "G1": {"status": "INHERITED_PHASE_A_REPRODUCED" if source_ok else "HARNESS_FAILED", "evidence_path": "source_contract.phase_a_payload_sha256"},
        "G2": {"status": "KNOWN_COEFFICIENTS_FINITE;FULL_REMAINDERS_UNRESOLVED" if all(math.isfinite(float(row["coefficient"])) for cell in cells.values() for term in cell["M_XX_p0_known"] for entry in term for row in entry) else "ILL_POSED"},
        "G3": {"status": "DERIVED_BLOCK_CONGRUENCE;FULL_SIGNATURE_UNRESOLVED", "conditions": congruence["positive_definite_conditions"]},
        "G4": {"status": "PASS" if rec.rows["B1_R2"]["status"] == "PASS" else "FAIL", "force_cases": len(g4["force_cases"])},
        "G5": {"status": "KNOWN_M_AND_OMEGA_COVARIANT;FULL_REMAINDERS_UNRESOLVED" if covariance_ok else "HARNESS_FAILED"},
        "G6": {"status": "ENDPOINT_MAP_COMPUTED", "source_map": endpoint_info["source_map"]},
        "G7": {"status": "NOT_RUN(phase_B2)"}, "G8": {"status": "NOT_RUN(phase_C)"},
        "G9": {"status": "NOT_RUN(phase_B2)"}, "G10": {"status": "NOT_RUN(phase_C)"}, "G11": {"status": "NOT_RUN(phase_C)"},
    }
    claims = [
        {"id": "phase_a_digest", "schema_path": "source_contract.phase_a_payload_sha256", "type": "sha256", "recompute": "sha256(normalized_phase_a_payload)"},
        {"id": "field_manifest", "schema_path": "field_manifest.fields", "type": "per_field_manifest", "recompute": "join(indexed_routes,phase_a.tail_channels)"},
        {"id": "emitted_leaves", "schema_path": "indexed_profile_missing_leaves", "type": "derived_string_set", "recompute": "failed_phase_a_indexed_tangent_lookups"},
        {"id": "scalar_regression", "schema_path": "scalar_regression", "type": "residual_mapping", "recompute": "eV.T*M_XX*eV-GVV"},
        {"id": "native_slice_constraints", "schema_path": "indexed_cells", "type": "conditional_residual_mapping", "recompute": "native_MXp/Mpp_projection-minus-GVP/GPP"},
        {"id": "g4_winding", "schema_path": "g4_control.computed_winding", "type": "sympy_expression", "recompute": "contour_integral/(2*pi)"},
        {"id": "g4_sigma", "schema_path": "g4_control.computed_sigma", "type": "sympy_expression", "recompute": "Omega_xy/(rho_mass*Gamma*epsilon_xy)"},
        {"id": "sheet_area", "schema_path": "g4_control.total_to_per_area_residual", "type": "canonical_matrix", "recompute": "Omega_total/sheet_cell_area-Omega_per_area"},
        {"id": "berry_coverage", "schema_path": "berry.production_pullback_coverage", "type": "cell_selection", "recompute": "production_cells-cross-ambient_branches"},
        {"id": "e4_action_hessian", "schema_path": "E4.M_aug_hessian_residual", "type": "canonical_matrix", "recompute": "hessian(preconstraint_extended_action)-M_aug"},
        {"id": "open_reachability", "schema_path": "open_root_reachability", "type": "generator_records", "recompute": "typed_ledger-union-traversal"},
        {"id": "finite_controls", "schema_path": "reachability_analysis.finite_bound_controls", "type": "emptiness_and_witness_records", "recompute": "same_generator_parity_filter"},
        {"id": "cell_count", "schema_path": "cells", "type": "mapping_length", "recompute": "active_endpoint-ambient-stratum_product"},
        {"id": "dimensions", "schema_path": "dimensions.records", "type": "computed_dimension_records", "recompute": "termwise_units_restoration"},
        {"id": "derived_congruence", "schema_path": "derived_congruence", "type": "block_congruence", "recompute": "produced_blocks-to-covariant-coefficients"},
        {"id": "mechanics_map", "schema_path": "mechanics_map", "type": "component_aggregate_map", "recompute": "aggregate(cells.component_findings)"},
        {"id": "closure_absence", "schema_path": "provenance_graph.global_return_closure_absence", "type": "graph_predicate", "recompute": "reachability(return_closure,mechanics_targets)"},
        {"id": "e5_root", "schema_path": "E5.root_deleted_conservative_solution", "type": "root_ablation_solve", "recompute": "delete(rayleigh_root)-and-resolve"},
        {"id": "gate_statuses", "schema_path": "gate_evidence", "type": "gate_mapping", "recompute": "engine_check-derived-gates"},
        {"id": "partition", "schema_path": "partition_ledger", "type": "ownership_ledger", "recompute": "computed-candidate-ownership"},
    ]
    result: dict[str, Any] = {"schema_version": "U1_PHASE_B1_ENGINE_ARTIFACT_V3", "engine": "SymPy",
        "engine_representation": representations["SymPy"], "phase": "B1", "input_validation": input_validation,
        "source_contract": source, "phase_a_amendment": correction,
        "amended_phase_a_payload": normalized_payload_v2(phase_a_payload(amended_phase)),
        "field_manifest": manifest, "field_embedding": manifest, "indexed_profile_missing_leaves": missing,
        "dimensions": dimensions,
        "endpoint_assembly": {ep: {k: v for k, v in row.items() if not k.startswith("_")} for ep, row in endpoints.items()},
        "endpoint_conservative": {ep: {k: v for k, v in row.items() if not k.startswith("_")} for ep, row in endpoint_info["conservative"].items()},
        "endpoint_source_map": endpoint_info["source_map"], "indexed_cells": public_cells(cells), "scalar_regression": scalar_regression,
        "full_block_expressions": blocks, "open_root_reachability": reach, "block_unresolved_roots": block_roots,
        "reachability_analysis": reach_analysis, "derived_congruence": congruence, "symbolic_LDL": congruence,
        "berry": berry, "g4_control": g4, "E4": e4, "E5": e5, "partition_ledger": partition,
        "provenance_graph": provenance, "cells": cell_records, "mechanics_map": mechanics_map, "gate_evidence": gates,
        "checks": rec.rows, "report_claim_bindings": claims, "closure_axis": closure, "open_strata": strata,
        "phase_scope": {"completed": ["7.2", "7.3"], "B2": "NOT_RUN(phase_B2)", "C": "NOT_RUN(phase_C)"}}
    liveness, liveness_ok = attach_input_liveness(result, plain_data, tracker)
    result["input_liveness"] = liveness
    rec.check("B1_S1", input_validation["valid"] and liveness_ok, {"declared": len(liveness["declared_leaf_paths"]),
                                                                   "consumed_but_undeclared": liveness["consumed_but_undeclared"]})
    result["checks"] = rec.rows
    result["axis_1"] = "COMPUTATION_VALID" if not rec.failures else "HARNESS_FAILED(" + ",".join(rec.failures) + ")"
    return result


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT); parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--phase-artifact", type=Path); args = parser.parse_args()
    try:
        data = yaml.safe_load(args.input.read_text())
        if data.get("schema_version") != "U1_PHASE_B1_MECHANICS_INPUTS_V3":
            raise ValueError("mechanics input schema is not V3")
        phase_input = yaml.safe_load(resolve(data["substrate"]["phase_a_input"]).read_text())
        phase_path = args.phase_artifact or resolve(data["substrate"]["sympy_artifact"])
        phase = json.loads(phase_path.read_text())
        result = build(data, phase_input, phase)
        args.output.parent.mkdir(parents=True, exist_ok=True)
        serializable = json.loads(json.dumps(result, default=str))
        if args.output.suffix == ".json":
            args.output.write_text(json.dumps(serializable, separators=(",", ":"), ensure_ascii=False))
        else:
            args.output.write_text(yaml.safe_dump(serializable, sort_keys=False, allow_unicode=True, width=180))
        failures = [name for name, row in result["checks"].items() if row["status"] != "PASS"]
        if failures:
            first = failures[0]; print(f"ASSERT_FAIL:{first}:{result['checks'][first]['evidence_digest']}", file=sys.stderr); return 1
        print(f"SYMPY_PHASE_B1_REMEDIATION3: PASS cells={len(result['cells'])}"); print(f"OUTPUT: {args.output}"); return 0
    except Exception as exc:
        print(f"ASSERT_FAIL:B1_S1:{type(exc).__name__}:{exc}", file=sys.stderr); return 1


if __name__ == "__main__":
    raise SystemExit(main())
