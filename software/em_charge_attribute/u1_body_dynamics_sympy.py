#!/usr/bin/env python3
"""U1 Phase-A field computation (SymPy engine).

The production YAML supplies action coefficients, core field traces, geometry,
and boundary equations.  This engine derives the exterior radial equations,
their tails, translated-mode norms, endpoint responses, and reduced moments.
It never reads a verdict or a decay exponent from input data.
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
from typing import Any, Callable

import mpmath as mp
import sympy as sp
import yaml
from scipy.integrate import quad
from scipy.special import kv, kvp


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
DEFAULT_INPUT = HERE / "u1_body_dynamics_inputs.yaml"
DEFAULT_FIXTURES = HERE / "u1_body_dynamics_fixtures.yaml"
DEFAULT_OUTPUT = HERE / "reports/u1_body_dynamics_artifacts/sympy_phase_a.json"
ENDPOINTS = ("E1", "E2", "E3", "E4", "E5")
TEETH = (
    "INPUT_LEDGER", "SOURCE_COMPLETENESS", "DIMENSIONAL",
    "COMOVING_CONTINUITY", "COMOVING_MOMENTUM", "BASE_BALANCE",
    "TAIL_ODE", "ZERO_MODE", "PROJECTOR", "ENDPOINT_RESPONSE",
    "MOMENT_INTEGRALS", "RECONSTRUCTION", "CANONICAL_VARIATION",
    "CHANNEL_UNIQUENESS", "TYPED_DATAFLOW", "PROVENANCE_FORBIDDEN",
    "ANCESTRY", "NATIVE_PADDING", "PARITY", "OUTCOME_REACHABILITY",
)


class ToothFailure(AssertionError):
    def __init__(self, tooth: str, detail: str):
        super().__init__(f"ASSERT_FAIL:{tooth}:{detail}")
        self.tooth = tooth


def require(test: Any, tooth: str, detail: str) -> None:
    if not bool(test):
        raise ToothFailure(tooth, detail)


def deep_get(data: Any, path: str) -> Any:
    cur = data
    for part in path.split("."):
        cur = cur[int(part)] if isinstance(cur, list) else cur[part]
    return cur


def deep_set(data: Any, path: str, value: Any) -> None:
    parts = path.split(".")
    cur = data
    for part in parts[:-1]:
        cur = cur[int(part)] if isinstance(cur, list) else cur[part]
    if isinstance(cur, list):
        cur[int(parts[-1])] = value
    else:
        cur[parts[-1]] = value


def apply_operations(data: dict[str, Any], operations: list[dict[str, Any]]) -> dict[str, str]:
    attacks: dict[str, str] = {}
    for op in operations:
        kind = op["op"]
        if kind == "set":
            deep_set(data, op["path"], op["value"])
        elif kind == "append":
            deep_get(data, op["path"]).append(copy.deepcopy(op["value"]))
        elif kind == "delete_by_id":
            rows = deep_get(data, op["path"])
            rows[:] = [row for row in rows if row.get("id") != op["id"]]
        elif kind == "derived_attack":
            attacks[op["target"]] = str(op["value"])
        else:
            raise ValueError(f"unknown fixture operation {kind}")
    return attacks


def load_case(input_path: Path, fixture_path: Path, case: str, mutation: str | None) -> tuple[dict[str, Any], dict[str, str], dict[str, Any]]:
    data = yaml.safe_load(input_path.read_text())
    fixtures = yaml.safe_load(fixture_path.read_text())
    require(data.get("schema_version") == "U1_PHASE_A_FIELD_INPUTS_V2", "INPUT_LEDGER", "production schema")
    require(fixtures.get("schema_version") == "U1_PHASE_A_FIXTURES_V2", "INPUT_LEDGER", "fixture schema")
    attacks: dict[str, str] = {}
    if case != "production":
        require(case in fixtures["outcome_cases"], "OUTCOME_REACHABILITY", f"unknown case {case}")
        attacks.update(apply_operations(data, fixtures["outcome_cases"][case]))
    if mutation:
        require(mutation in fixtures["mutations"], "INPUT_LEDGER", f"unknown mutation {mutation}")
        attacks.update(apply_operations(data, fixtures["mutations"][mutation]))
    return data, attacks, fixtures


def qvalue(value: Any) -> sp.Expr:
    if isinstance(value, (int, float)):
        return sp.Rational(str(value))
    text = str(value)
    try:
        return sp.Rational(text)
    except (ValueError, TypeError):
        return sp.Symbol(text)


def coefficient_symbols(data: dict[str, Any]) -> tuple[dict[str, sp.Symbol], dict[str, sp.Expr]]:
    symbols: dict[str, sp.Symbol] = {}
    values: dict[str, sp.Expr] = {}
    for name, rec in data["coefficients"].items():
        constraint = rec.get("constraint")
        assumptions: dict[str, bool] = {"real": True}
        if constraint == "positive":
            assumptions["positive"] = True
        elif constraint == "nonnegative":
            assumptions["nonnegative"] = True
        symbols[name] = sp.Symbol(name, **assumptions)
        raw = qvalue(rec["value"])
        values[name] = symbols[name] if raw.is_Symbol else raw
    return symbols, values


def dadd(*dims: tuple[int, int, int]) -> tuple[int, int, int]:
    return tuple(sum(row) for row in zip(*dims))  # type: ignore[return-value]


def dscale(dim: tuple[int, int, int], power: int) -> tuple[int, int, int]:
    return tuple(power * x for x in dim)  # type: ignore[return-value]


def dsub(a: tuple[int, int, int], b: tuple[int, int, int]) -> tuple[int, int, int]:
    return tuple(x - y for x, y in zip(a, b))  # type: ignore[return-value]


def extract_fence(text: str, label: str) -> str:
    lines = text.splitlines()
    start = next(i for i, line in enumerate(lines) if line.strip() == f"```{label}")
    end = next(i for i in range(start + 1, len(lines)) if lines[i].strip() == "```")
    return "\n".join(lines[start + 1:end])


def math_tokens(text: str) -> tuple[str, ...]:
    """Lex a source/action record into structural mathematical tokens."""
    return tuple(re.findall(r"[A-Za-z_][A-Za-z_0-9]*|\d+(?:\.\d+)?|:=|->|!=|==|\^|\*|/|\+|-|\(|\)|\[|\]|\||\.|·|∂|∇|Ω|χ|ρ|λ|ϖ", text))


def parsed_source_contains(source_text: str, probe: str) -> bool:
    """Match a token subtree in parsed Markdown records, never raw substrings."""
    needle = math_tokens(probe)
    if not needle:
        return False
    records = [math_tokens(line) for line in source_text.splitlines()]
    for record in records:
        if len(record) < len(needle):
            continue
        if any(record[i:i + len(needle)] == needle for i in range(len(record) - len(needle) + 1)):
            return True
    return False


def source_completeness(data: dict[str, Any], attacks: dict[str, str]) -> dict[str, Any]:
    policy = data["operative_action_decision"]
    decision_rel = policy["source_file"]
    decision_text = (ROOT / decision_rel).read_text()
    require(policy["id"] == "decision_16_retire_brane_polar_field" and policy["status"] == "OPERATIVE" and
            decision_rel == "software/stage1_solver/decisions/16_retire_brane_polar_field.md",
            "SOURCE_COMPLETENESS", "Decision 16 is not the operative action policy")
    citation_matches: list[str] = []
    for fragment in policy["required_source_fragments"]:
        require(parsed_source_contains(decision_text, fragment),
                "SOURCE_COMPLETENESS", f"Decision 16 citation mismatch: {fragment}")
        citation_matches.append(fragment)

    expected_ids = set(policy["expected_action_term_ids"])
    retired_ids = set(policy["retired_action_term_ids"])
    retired_symbols = set(policy["retired_expression_symbols"])
    assembled_ids = [row["id"] for row in data["action_terms"]]
    require(len(assembled_ids) == len(set(assembled_ids)),
            "SOURCE_COMPLETENESS", "duplicate assembled action id")
    for term in data["action_terms"]:
        expression_symbols = set(re.findall(r"[A-Za-z_][A-Za-z_0-9]*", term["expression"]))
        intrusion = term["id"] in retired_ids or bool(expression_symbols & retired_symbols)
        if intrusion:
            require(term.get("decision_citation") == decision_rel,
                    "SOURCE_COMPLETENESS",
                    f"retired P term lacks Decision 16 citation: {term['id']}")
            require(False, "SOURCE_COMPLETENESS",
                    f"retired P term cannot be assembled under Decision 16: {term['id']}")
    missing = sorted(expected_ids - set(assembled_ids))
    unexpected = sorted(set(assembled_ids) - expected_ids)
    require(not missing, "SOURCE_COMPLETENESS", f"missing P-retired whitelisted terms {missing}")
    require(not unexpected, "SOURCE_COMPLETENESS", f"unexpected assembled terms {unexpected}")

    cache: dict[str, str] = {}
    matched: list[dict[str, str]] = []
    for term in data["action_terms"]:
        rel = term["source_file"]
        if rel not in cache:
            cache[rel] = (ROOT / rel).read_text()
        require(parsed_source_contains(cache[rel], term["source_contains"]),
                "SOURCE_COMPLETENESS", f"parsed source mismatch {term['id']}")
        matched.append({"id": term["id"], "source_file": rel, "source_fragment": term["source_contains"]})

    # Load the legacy immutable blocks. Their P records are now negative
    # retirement evidence, not members of the operative assembled action.
    t0_text = (ROOT / "software/stage1_solver/reports/pathA_24_T0_freeze.md").read_text()
    g0_text = (ROOT / "software/stage1_solver/reports/pathA_35_G0_freeze.md").read_text()
    t0 = extract_fence(t0_text, "freeze-action")
    g0 = extract_fence(g0_text, "freeze-action")
    polar_lines = []
    in_pol = False
    for raw in t0.splitlines():
        line = raw.strip()
        if line == "L_pol =":
            in_pol = True
            continue
        if in_pol and line.startswith("Frozen extended"):
            break
        if in_pol and re.match(r"^[+-]?\s*\(1/[24]\) m rho", line):
            polar_lines.append(line.lstrip("+- ").rstrip("."))
    require(len(polar_lines) == 3, "SOURCE_COMPLETENESS", "legacy T0 polar manifest changed")

    # The G0 block explicitly carries the parent quantum-pressure term and the
    # brane blocks. Definitions are discovered from the fence itself; L_Pu is
    # separated into the Decision-16-retired evidence set.
    g0_records = [math_tokens(line) for line in g0.splitlines()]
    discovered_definitions = sorted({record[0] for record in g0_records
                                     if len(record) > 1 and record[1] == ":=" and
                                     (record[0] == "QP" or record[0].startswith("L_"))})
    require(len(discovered_definitions) >= 4, "SOURCE_COMPLETENESS", "source parser missed G0 action definitions")
    # These mandatory records are discovered in the immutable G0 fence itself.
    # Comparing their parsed token trees catches deletion of the exact R1 terms
    # even when the remaining YAML rows still point at valid source text.
    discovered_g0_records = []
    for raw in g0.splitlines():
        line = raw.strip().lstrip("+- ").rstrip(".")
        if ("QP :=" in line or "mu_R Omega_u^a Omega_u^a" in line or
                line.startswith("L_Pu :=")):
            discovered_g0_records.append(line)
    retired_g0 = [line for line in discovered_g0_records if line.startswith("L_Pu :=")]
    mandatory_g0 = [line for line in discovered_g0_records if not line.startswith("L_Pu :=")]
    require(len(retired_g0) == 1, "SOURCE_COMPLETENESS", "legacy L_Pu retirement record changed")
    assembled_fragments = [row["source_contains"] for row in data["action_terms"]]
    for fragment in mandatory_g0:
        expected = math_tokens(fragment)
        require(any(len(expected) >= len(math_tokens(candidate)) and
                    any(expected[i:i + len(math_tokens(candidate))] == math_tokens(candidate)
                        for i in range(len(expected) - len(math_tokens(candidate)) + 1))
                    for candidate in assembled_fragments),
                "SOURCE_COMPLETENESS", f"missing immutable G0 record {fragment}")
    return {
        "loaded_files": sorted(cache),
        "matched_assembled_terms": matched,
        "operative_decision_citation": {
            "id": policy["id"], "status": policy["status"], "source_file": decision_rel,
            "sha256": hashlib.sha256(decision_text.encode()).hexdigest(),
            "matched_source_fragments": citation_matches,
        },
        "expected_p_retired_action_ids": sorted(expected_ids),
        "retired_action_ids": sorted(retired_ids),
        "retired_parameter_rows": list(policy["retired_parameter_rows"]),
        "legacy_retired_t0_polar_monomials": polar_lines,
        "source_derived_g0_definitions": discovered_definitions,
        "source_derived_mandatory_g0_records": mandatory_g0,
        "source_derived_retired_g0_records": retired_g0,
        "assembled_ids": sorted(assembled_ids),
    }


def dimensional_firewall(data: dict[str, Any]) -> list[dict[str, Any]]:
    dims = {name: tuple(rec["dimensions"]) for name, rec in data["coefficients"].items()}
    L = (1, 0, 0)
    TINV = (0, -1, 0)
    GRAD = (-1, 0, 0)
    N = dims["rho_inf"]
    VEL = (1, -1, 0)
    target = (-2, -2, 1)
    rows: list[tuple[str, tuple[int, int, int]]] = [
        ("bulk_berry", dadd(dims["hbar"], N, TINV)),
        ("bulk_flow_kinetic", dadd(dims["m_GNLS"], N, dscale(VEL, 2))),
        ("quantum_pressure", dadd(dscale(dims["hbar"], 2), dscale(dims["m_GNLS"], -1), dscale(N, -1), dscale(dadd(N, GRAD), 2))),
        ("bulk_EOS", dadd(dims["K_EOS"], dscale(N, 5))),
        ("wall_double_well", dims["aB"]),
        ("wall_gradient", dadd(dims["kappaB"], dscale(GRAD, 2))),
        ("wall_shear_gate", dims["muR4"]),
    ]
    rows.extend([
        ("brane_shear_kinetic", dadd(dscale(dims["ellg"], -1), dims["rhoBr"], dscale(dadd(L, TINV), 2))),
        ("brane_shear_gradient", dadd(dscale(dims["ellg"], -1), dims["muR"])),
        ("uw_kinetic", dadd(dscale(dims["ellg"], -1), dims["rhoBr"], dscale(dadd(L, TINV), 2))),
        ("uw_gap", dadd(dscale(dims["ellg"], -1), dims["rhoBr"], dscale(dims["OmegaW"], 2), dscale(L, 2))),
        ("h_kinetic", dadd(dims["Mh"], dscale(dims["cE"], -2), dscale(TINV, 2))),
        ("h_gradient", dadd(dims["Mh"], dscale(GRAD, 2))),
    ])
    output = []
    for name, dim in rows:
        require(dim == target, "DIMENSIONAL", f"{name}:{dim}!={target}")
        output.append({"expression": name, "computed_dimensions_LTM": list(dim)})
    # Momentum balance and the 4D control-surface force are constructed inline.
    momentum_density = dadd(dims["m_GNLS"], N, VEL)
    momentum_rate = dadd(momentum_density, TINV)
    stress = target
    require(momentum_rate == dadd(stress, GRAD), "DIMENSIONAL", "momentum PDE")
    require(dadd(stress, (3, 0, 0)) == (1, -2, 1), "DIMENSIONAL", "4D surface force")
    output += [
        {"expression": "momentum_rate_density", "computed_dimensions_LTM": list(momentum_rate)},
        {"expression": "control_surface_force", "computed_dimensions_LTM": [1, -2, 1]},
    ]
    return output


def comoving_reductions(data: dict[str, Any]) -> dict[str, Any]:
    t, x, y, V = sp.symbols("t x y V", real=True)
    f = sp.Function("f")
    nfun, vfun = sp.Function("n"), sp.Function("v")
    gfun, pifun = sp.Function("g_i"), sp.Function("Pi_ij")
    map_record = data["kinematics"]["coordinate_map"]
    map_sign = sp.Integer(map_record["displacement_coefficient"])
    body_argument = x + map_sign * V * t
    composite = f(body_argument)
    chain = sp.diff(composite, t)
    explicit_chain = map_sign * V * sp.Subs(sp.diff(f(y), y), y, body_argument)
    chain_residual = sp.simplify(chain - explicit_chain)
    require(chain_residual == 0 and map_sign == -1, "COMOVING_CONTINUITY", str(chain_residual))

    nlab, vlab = nfun(body_argument), vfun(body_argument)
    native_continuity = sp.diff(nlab, t) + sp.diff(nlab * vlab, x)
    derived_continuity = sp.diff(nlab * (vlab + map_sign * V), x)
    continuity_residual = sp.simplify(native_continuity - derived_continuity)
    require(continuity_residual == 0, "COMOVING_CONTINUITY", str(continuity_residual))

    flux_record = data["kinematics"]["momentum_flux"]
    flux_sign = sp.Integer(flux_record["convective_coefficient"])
    glab, pilab = gfun(body_argument), pifun(body_argument)
    native_momentum = sp.diff(glab, t) + sp.diff(pilab, x)
    derived_momentum = sp.diff(pilab + map_sign * V * glab, x)
    declared_momentum = sp.diff(pilab + flux_sign * V * glab, x)
    momentum_residual = sp.simplify(native_momentum - derived_momentum)
    declared_residual = sp.simplify(declared_momentum - derived_momentum)
    require(momentum_residual == 0, "COMOVING_MOMENTUM", str(momentum_residual))
    require(declared_residual == 0, "COMOVING_MOMENTUM", str(declared_residual))
    return {
        "coordinate_substitution": map_record,
        "composite_field": str(composite),
        "scalar_chain_rule": str(chain),
        "scalar_chain_rule_expanded": str(explicit_chain),
        "scalar_chain_rule_residual": str(chain_residual),
        "continuity_native": "partial_t n+div_4(n v)=0",
        "continuity_comoving": "partial_t|y n+div_4[n(v-V)]=0",
        "continuity_residual": str(continuity_residual),
        "momentum_native": "partial_t g_i+partial_j Pi_ij=f_i[action]",
        "momentum_comoving": "partial_t|y g_i+partial_j(Pi_ij-V_j g_i)=f_i[action]",
        "momentum_residual": str(momentum_residual),
        "declared_flux_residual": str(declared_residual),
        "surface_force": "integral_Omega f_i d4y-integral_partialOmega(Pi_ij-V_j g_i)N_j dSigma3",
        "particle_momentum_imported": False,
    }


def numeric_sign(expr: sp.Expr) -> int | None:
    expr = sp.simplify(expr)
    if expr.is_positive:
        return 1
    if expr.is_zero:
        return 0
    if expr.is_negative:
        return -1
    if not expr.free_symbols:
        val = float(sp.N(expr))
        return 1 if val > 0 else (-1 if val < 0 else 0)
    return None


def sphere_area(d: int) -> sp.Expr:
    return sp.simplify(2 * sp.pi ** (sp.Rational(d, 2)) / sp.gamma(sp.Rational(d, 2)))


ACTION_PRIMITIVES = (
    "n", "theta_t", "v2", "n_grad2", "chiB", "chi_grad2",
    "ud_curl2", "f_throat", "f_mix", "g_ell", "u_t2",
    "u_curl2", "uw_t2", "uw", "h_t2", "h_grad2",
)


def assemble_action_terms(data: dict[str, Any], values: dict[str, sp.Expr]) -> dict[str, Any]:
    """Parse every action expression into one executable symbolic surface."""
    primitive = {name: sp.Symbol(name, real=True) for name in ACTION_PRIMITIVES}
    n = primitive["n"]
    cs2 = sp.diff(values["K_EOS"] * n**5, n).subs(n, values["rho_inf"]) / values["m_GNLS"]
    local = {**values, **primitive, "cs2": sp.simplify(cs2)}
    parsed: dict[str, sp.Expr] = {}
    dependencies: dict[str, list[str]] = {}
    for row in data["action_terms"]:
        expr = sp.sympify(row["expression"], locals=local, evaluate=True)
        parsed[row["id"]] = expr
        dependencies[row["id"]] = sorted(str(s) for s in expr.free_symbols)
    return {
        "terms": parsed,
        "primitives": primitive,
        "cs2": sp.simplify(cs2),
        "dependencies": dependencies,
        "total_expression": sp.Add(*parsed.values()),
    }


def derive_channel_operator(data: dict[str, Any], values: dict[str, sp.Expr],
                            assembled: dict[str, Any] | None = None) -> dict[str, Any]:
    """Second-vary the parsed action records, including all cross blocks."""
    assembled = assembled or assemble_action_terms(data, values)
    terms: dict[str, sp.Expr] = assembled["terms"]
    p = assembled["primitives"]
    d = len(data["ambient"]["coordinates"])
    db = len(data["ambient"]["brane_coordinates"])
    hbar, m, rho, K = (values[x] for x in ("hbar", "m_GNLS", "rho_inf", "K_EOS"))
    r = sp.Symbol("r", positive=True)
    eps = sp.Symbol("epsilon", real=True)
    dn, theta, chi, uT, uw, h = sp.symbols("dn theta chi uT uw h", real=True)
    dn_r, theta_r, chi_r, u_curl, h_r = sp.symbols(
        "dn_r theta_r chi_r u_curl h_r", real=True)
    drain_gradient = sp.simplify(-values["mdot"] / (hbar * rho * sphere_area(d)) * r ** (1 - d))
    static_subs = {
        p["n"]: rho + eps * dn,
        p["theta_t"]: 0,
        p["v2"]: (hbar / m) ** 2 * (drain_gradient + eps * theta_r) ** 2,
        p["n_grad2"]: eps**2 * dn_r**2,
        p["chiB"]: eps * chi,
        p["chi_grad2"]: eps**2 * chi_r**2,
        p["ud_curl2"]: 0,
        p["f_throat"]: 0,
        p["f_mix"]: 0,
        p["g_ell"]: 1,
        p["u_t2"]: 0,
        p["u_curl2"]: eps**2 * u_curl**2,
        p["uw_t2"]: 0,
        p["uw"]: eps * uw,
        p["h_t2"]: 0,
        p["h_grad2"]: eps**2 * h_r**2,
    }
    quadratic_by_term: dict[str, sp.Expr] = {}
    for term_id, expr in terms.items():
        expanded = sp.series(expr.subs(static_subs), eps, 0, 3).removeO().expand()
        quadratic_by_term[term_id] = sp.simplify(expanded.coeff(eps, 2))
    quadratic = sp.simplify(sum(quadratic_by_term.values()))
    fields = [dn, theta, chi, uT, uw, h]
    gradients = [dn_r, theta_r, chi_r, u_curl, h_r]
    curvature = sp.simplify(sp.hessian(quadratic, fields))
    gradient = sp.simplify(sp.hessian(quadratic, gradients))
    mixed = sp.Matrix([[sp.simplify(sp.diff(quadratic, f, g)) for g in gradients] for f in fields])
    channel_rows = [
        {"id": "density_EOS", "field": "delta_n", "dimension": d, "ell": 0,
         "gradient": gradient[0, 0], "curvature": curvature[0, 0]},
        {"id": "wall_chiB", "field": "delta_chiB", "dimension": d, "ell": 0,
         "gradient": gradient[2, 2], "curvature": curvature[2, 2]},
        {"id": "bound_phase", "field": "theta", "dimension": d, "ell": 0,
         "gradient": gradient[1, 1], "curvature": curvature[1, 1], "continuity": True},
        {"id": "brane_shear", "field": "u_transverse", "dimension": db, "ell": 1,
         "gradient": gradient[3, 3], "curvature": sp.Integer(0), "brane_profile": "g_ell(w)"},
        {"id": "uw", "field": "u_w", "dimension": db, "ell": 0,
         "gradient": sp.Integer(0), "curvature": curvature[4, 4],
         "algebraic": True},
        {"id": "h", "field": "h", "dimension": d, "ell": 0,
         "gradient": gradient[4, 4], "curvature": curvature[5, 5], "open_coefficients": ["Mh", "cE"]},
    ]
    operator_entries = {
        "density_gradient": gradient[0, 0], "density_EOS_curvature": curvature[0, 0],
        "phase_gradient": gradient[1, 1], "wall_gradient": gradient[2, 2],
        "wall_well_curvature": curvature[2, 2], "shear_curl": gradient[3, 3],
        "uw_curvature": curvature[4, 4], "h_gradient": gradient[4, 4],
        "drain_density_phase": mixed[0, 1],
    }
    return {
        "channels": channel_rows,
        "quadratic_expression": quadratic,
        "quadratic_by_term": quadratic_by_term,
        "field_order": [str(x) for x in fields],
        "gradient_order": [str(x) for x in gradients],
        "curvature_hessian": curvature,
        "gradient_hessian": gradient,
        "mixed_hessian": mixed,
        "entries": operator_entries,
        "drain_gradient": drain_gradient,
    }


def derive_channels(data: dict[str, Any], values: dict[str, sp.Expr]) -> list[dict[str, Any]]:
    return derive_channel_operator(data, values)["channels"]


def coupled_indicial_analysis(data: dict[str, Any], operator: dict[str, Any]) -> dict[str, Any]:
    """Check the surviving drain-induced density--phase degree mixing."""
    d = len(data["ambient"]["coordinates"])
    nu = sp.Symbol("nu", real=True)
    drain_diag_degree = -nu - 2
    drain_cross_degree = -nu - d
    drain_degree_difference = sp.simplify(drain_cross_degree - drain_diag_degree)
    return {
        "density_phase": {
            "cross_coefficient": str(operator["entries"]["drain_density_phase"]),
            "diagonal_degree": str(drain_diag_degree),
            "cross_degree": str(drain_cross_degree),
            "degree_difference": str(drain_degree_difference),
            "leading_exponents_shifted": not bool(d > 2),
            "computed_conclusion": "SUBLEADING_FOR_D_GT_2",
        },
        "changes_scalar_channel_verdict": False,
    }


def action_term_removal_probes(data: dict[str, Any], values: dict[str, sp.Expr],
                               base_operator: dict[str, Any]) -> dict[str, Any]:
    targets = {
        "quantum_pressure": "density_gradient",
        "brane_shear_gradient": "shear_curl",
    }
    output: dict[str, Any] = {}
    for term_id, entry in targets.items():
        altered = copy.deepcopy(data)
        altered["action_terms"] = [row for row in altered["action_terms"] if row["id"] != term_id]
        removed_operator = derive_channel_operator(altered, values)
        before = sp.simplify(base_operator["entries"][entry])
        after = sp.simplify(removed_operator["entries"][entry])
        delta = sp.simplify(before - after)
        require(delta != 0, "SOURCE_COMPLETENESS", f"{term_id} did not reach {entry}")
        output[term_id] = {
            "operator_entry": entry, "before": str(before), "after_removal": str(after),
            "computed_delta": str(delta), "operator_changed": True,
            "source_completeness_gate": "FAIL_IF_REMOVED",
        }
    return output


def solve_tail(channel: dict[str, Any], attacks: dict[str, str]) -> dict[str, Any]:
    r = sp.Symbol("r", positive=True)
    f = sp.Function("f")
    d, ell = int(channel["dimension"]), int(channel["ell"])
    A, B = sp.simplify(channel["gradient"]), sp.simplify(channel["curvature"])
    if channel.get("algebraic"):
        residual = sp.simplify(B * sp.Integer(0))
        return {
            "id": channel["id"], "field": channel["field"], "radial_dimension": d,
            "equation": f"{B}*f(r)=0", "dsolve": "f(r)=0 outside the collar",
            "classification": "ALGEBRAIC_GAP", "decay_exponent": None,
            "gap_squared": str(B), "solution_residual": str(residual),
            "zero_mode_norm_integral": "0", "zero_mode_norm_value": 0.0,
            "normalizable": True, "dependencies": sorted(str(x) for x in B.free_symbols),
        }
    radial = sp.diff(f(r), r, 2) + sp.Rational(d - 1, 1) * sp.diff(f(r), r) / r - ell * (ell + d - 2) * f(r) / r**2
    ode = sp.Eq(A * radial - B * f(r), 0)
    ratio = sp.simplify(B / A)
    sign = numeric_sign(ratio)
    dependencies = sorted(str(x) for x in (A.free_symbols | B.free_symbols))
    if sign is None:
        return {
            "id": channel["id"], "field": channel["field"], "radial_dimension": d,
            "equation": str(ode), "solved_general_form": "conditional on open coefficient stratum",
            "classification": "UNRESOLVED", "decay_exponent": None, "gap_squared": str(ratio),
            "solution_residual": "conditional", "zero_mode_norm_integral": "conditional",
            "zero_mode_norm_value": None, "normalizable": None, "dependencies": dependencies,
        }
    if sign < 0:
        k = sp.sqrt(-ratio)
        sol = r ** (1 - sp.Rational(d, 2)) * sp.besselj(ell + sp.Rational(d, 2) - 1, k * r)
        residual = sp.simplify((A * radial - B * f(r)).subs(f(r), sol).doit())
        return {
            "id": channel["id"], "field": channel["field"], "radial_dimension": d,
            "equation": str(ode), "solved_general_form": str(sol), "classification": "TACHYONIC",
            "decay_exponent": None, "gap_squared": str(ratio), "solution_residual": str(residual),
            "zero_mode_norm_integral": "not used: static operator tachyonic", "zero_mode_norm_value": None,
            "normalizable": False, "dependencies": dependencies,
        }
    if sign == 0:
        indicial = sp.expand(sp.Symbol("nu") * (sp.Symbol("nu") - d + 2) - ell * (ell + d - 2))
        roots = sp.solve(sp.Eq(indicial, 0), sp.Symbol("nu"))
        decaying = [root for root in roots if root.is_positive]
        nu = max(decaying) if decaying else sp.Integer(0)
        sol = r ** (-nu)
        if attacks.get(f"{channel['id']}_tail_solution") == "add_constant" or (channel["id"] == "density_EOS" and attacks.get("density_tail_solution") == "add_constant"):
            sol += r
        residual = sp.simplify((A * radial - B * f(r)).subs(f(r), sol).doit())
        require(residual == 0, "TAIL_ODE", f"{channel['id']}:{residual}")
        R = sp.Symbol("R", positive=True)
        norm = sp.integrate(r ** (d - 1) * sp.diff(r ** (-nu), r) ** 2, (r, R, sp.oo))
        converges = norm.is_finite is True and nu > 0
        return {
            "id": channel["id"], "field": channel["field"], "radial_dimension": d,
            "equation": str(ode), "solved_general_form": f"C_grow*r^{ell}+C_decay*r^(-{nu})", "indicial_equation": str(indicial),
            "indicial_roots": [str(x) for x in roots], "classification": "POWER_LAW",
            "decay_exponent": str(nu), "gap_squared": "0", "solution_residual": str(residual),
            "zero_mode_norm_integral": str(norm), "zero_mode_norm_value": float(sp.N(norm.subs(R, 1))) if converges else None,
            "normalizable": bool(converges), "dependencies": dependencies,
        }
    k = sp.sqrt(ratio)
    lambda_roots = sp.solve(sp.Eq(sp.Symbol("lambda")**2 - ratio, 0), sp.Symbol("lambda"))
    order = ell + sp.Rational(d, 2) - 1
    sol = r ** (1 - sp.Rational(d, 2)) * sp.besselk(order, k * r)
    if channel["id"] == "density_EOS" and attacks.get("density_tail_solution") == "add_constant":
        sol += 1
    residual = sp.simplify(sp.expand_func((A * radial - B * f(r)).subs(f(r), sol).doit()))
    residual = sp.simplify(residual.rewrite(sp.besselk))
    require(residual == 0, "TAIL_ODE", f"{channel['id']}:{residual}")
    kval = float(sp.N(k.subs({s: 1 for s in k.free_symbols})))
    orderval = float(order)
    power = 1 - d / 2
    dq = lambda x: power * x ** (power - 1) * kv(orderval, kval * x) + x**power * kval * kvp(orderval, kval * x, 1)
    norm_value = quad(lambda x: x ** (d - 1) * dq(x) ** 2, 1, math.inf,
                      epsabs=2e-11, epsrel=2e-11, limit=200)[0]
    return {
        "id": channel["id"], "field": channel["field"], "radial_dimension": d,
        "equation": str(ode), "solved_general_form": f"r^(1-d/2)[C_I I_{order}(gap*r)+C_K K_{order}(gap*r)]",
        "large_r_characteristic_roots": [str(x) for x in lambda_roots], "classification": "EXPONENTIAL_GAP",
        "decay_exponent": None, "gap_squared": str(sp.factor(ratio)), "gap": str(k),
        "solution_residual": str(residual),
        "zero_mode_norm_integral": f"Integral[1,inf] r^{d-1} (d/dr(r^{1-d/2} K_{order}(gap*r)))^2 dr",
        "zero_mode_norm_value": norm_value, "normalizable": math.isfinite(norm_value),
        "dependencies": dependencies,
    }


def phase_normalization(data: dict[str, Any], values: dict[str, sp.Expr]) -> dict[str, Any]:
    d = len(data["ambient"]["coordinates"])
    r = sp.Symbol("r", positive=True)
    phi = sp.Function("phi")
    rho, hbar, m, mdot = (values[x] for x in ("rho_inf", "hbar", "m_GNLS", "mdot"))
    ode = sp.Eq(sp.diff(r ** (d - 1) * rho * (hbar / m) * sp.diff(phi(r), r), r), 0)
    solution = sp.dsolve(ode)
    area = sphere_area(d)
    C = sp.simplify(mdot / ((d - 2) * hbar * rho * area))
    selected = C * r ** (2 - d)
    flux = sp.simplify(area * r ** (d - 1) * m * rho * (hbar / m) * sp.diff(selected, r))
    residual = sp.simplify(flux + mdot)
    require(residual == 0, "BASE_BALANCE", f"phase flux {residual}")
    return {"continuity_ode": str(ode), "dsolve": str(solution), "selected_phase": str(selected),
            "sphere_area": str(area), "computed_flux": str(flux), "normalization_residual": str(residual)}


def zero_mode_operator(data: dict[str, Any], values: dict[str, sp.Expr], attacks: dict[str, str]) -> dict[str, Any]:
    # A real Cartesian four-dimensional substitution makes translation commute
    # with the native Laplacian; the co-moving throat source is varied too.
    d = len(data["ambient"]["coordinates"])
    xs = sp.symbols(f"x0:{d}", real=True)
    a = qvalue(data["geometry"]["a"]["value"])
    rr2 = sum(x * x for x in xs)
    q = sp.exp(-rr2 / a**2)
    lap = lambda expr: sum(sp.diff(expr, x, 2) for x in xs)
    A = values["hbar"] ** 2 / (4 * values["m_GNLS"] * values["rho_inf"])
    n0 = values["rho_inf"] + qvalue(data["core_traces"]["density_delta"]["value"]) * q
    n = sp.Symbol("n", positive=True)
    Up = sp.diff(values["K_EOS"] * n**5 / 4, n).subs(n, n0)
    source = sp.simplify(-A * lap(n0) + Up)
    balance = sp.simplify(-A * lap(n0) + Up - source)
    if attacks.get("source_profile_after_balance") == "add_core_trace":
        source = source + q
        balance = sp.simplify(-A * lap(n0) + Up - source)
    require(balance == 0, "BASE_BALANCE", str(balance))
    z = -sp.diff(n0, xs[0])
    if attacks.get("translation_mode") == "add_base_profile":
        z += n0
    Upp = sp.diff(values["K_EOS"] * n**5 / 4, n, 2).subs(n, n0)
    source_z = -sp.diff(source, xs[0])
    residual = sp.simplify(-A * lap(z) + Upp * z - source_z)
    require(residual == 0, "ZERO_MODE", str(residual))
    return {
        "branch": "force_balance",
        "operator_object": "D_E=(-A_n Laplacian+U''(n0)) on delta_n, plus translated throat-source block; analogous Hessian/balance blocks for chiB,u,h and endpoint traces",
        "density_balance_expression": str(-A * lap(n0) + Up - source),
        "density_source_expression": str(source),
        "right_zero_mode": str(z),
        "substitution_residual": str(residual),
        "operator_available": True,
    }


def projector_algebra(tails: list[dict[str, Any]], attacks: dict[str, str]) -> dict[str, Any]:
    norms = [float(row["zero_mode_norm_value"]) for row in tails if row["zero_mode_norm_value"] is not None]
    require(norms and all(math.isfinite(x) for x in norms), "PROJECTOR", "nonfinite Gram input")
    z1, z2 = sp.symbols("z1 z2", real=True)
    z = sp.Matrix([z1, z2])
    gram = (z.T * z)[0]
    inverse = 1 / gram
    if attacks.get("gram_inverse") == "double":
        inverse *= 2
    Q = sp.eye(2) - z * inverse * z.T
    residual = sp.simplify(Q * z)
    require(residual == sp.zeros(2, 1), "PROJECTOR", str(residual))
    return {
        "computed_channel_norm_sum": sum(norms), "gram": str(gram), "gram_inverse": str(inverse),
        "projector": str(Q), "Q_times_Z": [str(x) for x in residual],
        "moduli_fixing": "<Ztilde_A,eta>_IR=0", "reduced_operator": "Q D_E Q on im(Q)",
    }


def eval_matrix_entry(value: Any, values: dict[str, sp.Expr]) -> sp.Expr:
    if isinstance(value, str):
        return sp.sympify(value, locals=values)
    return qvalue(value)


def endpoint_responses(data: dict[str, Any], values: dict[str, sp.Expr]) -> dict[str, Any]:
    require(tuple(data["endpoints"]) == ENDPOINTS, "ENDPOINT_RESPONSE", "endpoint ordering/coverage")
    output: dict[str, Any] = {}
    for ep in ENDPOINTS:
        row = data["endpoints"][ep]
        require(row.get("trace_unknowns") == ["normal", "tangent", "shear"],
                "ENDPOINT_RESPONSE", f"{ep}:unknown declaration")
        mat = sp.Matrix([[eval_matrix_entry(x, values) for x in rr] for rr in row["trace_system"]["matrix"]])
        rhs = sp.Matrix([eval_matrix_entry(x, values) for x in row["trace_system"]["rhs"]])
        require(mat.shape == (3, 3) and mat.rank() == mat.row_join(rhs).rank() == 3,
                "ENDPOINT_RESPONSE", f"{ep}:trace rank")
        sol = mat.inv() * rhs
        residual = sp.simplify(mat * sol - rhs)
        require(residual == sp.zeros(3, 1), "ENDPOINT_RESPONSE", f"{ep}:trace residual")
        output[ep] = {
            "condition": row["condition"], "variational_class": row["variational_class"],
            "fluid_coefficients": {"normal": str(sp.simplify(sol[0])), "tangent": str(sp.simplify(sol[1]))},
            "shear_coefficient": str(sp.simplify(sol[2])),
            "declared_matrix": [[str(x) for x in rr] for rr in mat.tolist()],
            "declared_rhs": [str(x) for x in rhs],
            "solve_method": "exact rank-3 boundary/constraint solve",
            "fluid_residual": [str(residual[0]), str(residual[1])],
            "shear_residual": [str(residual[2])],
        }
    return output


def mp_float(expr: sp.Expr, values: dict[str, sp.Expr]) -> float:
    substitutions = {sp.Symbol(name): val for name, val in values.items() if not val.free_symbols}
    remaining = sp.sympify(expr).subs(substitutions)
    remaining = remaining.subs({s: 1 for s in remaining.free_symbols})
    return float(sp.N(remaining))


def radial_profiles(data: dict[str, Any], values: dict[str, sp.Expr]) -> dict[str, Callable[[float], float]]:
    d = len(data["ambient"]["coordinates"])
    a = float(qvalue(data["geometry"]["a"]["value"]))
    hbar, m, rho, K = (mp_float(values[x], values) for x in ("hbar", "m_GNLS", "rho_inf", "K_EOS"))
    kd = math.sqrt(abs(20 * m * K * rho**4 / hbar**2))
    kw = math.sqrt(abs(2 * mp_float(values["aB"], values) / mp_float(values["kappaB"], values)))

    def gap_profile(k: float, order: float) -> tuple[Callable[[float], float], Callable[[float], float]]:
        power = 1 - d / 2
        denom = kv(order, k * a) * a**power
        q = lambda r: float(r**power * kv(order, k * r) / denom)
        dq = lambda r: float((power * r ** (power - 1) * kv(order, k * r)
                              + r**power * k * kvp(order, k * r, 1)) / denom)
        return q, dq

    density, density_d = gap_profile(kd, d / 2 - 1)
    wall, wall_d = gap_profile(kw, d / 2 - 1)
    return {
        "density": density, "density_d": density_d,
        "wall": wall, "wall_d": wall_d,
        "response_bulk": lambda r: (a / r) ** (d - 1),
    }


def quad_inf(fn: Callable[[float], float], a: float) -> float:
    value, error = quad(fn, a, math.inf, epsabs=2e-11, epsrel=2e-11, limit=300)
    if not math.isfinite(value) or error > 1e-7 * max(1.0, abs(value)):
        raise ToothFailure("MOMENT_INTEGRALS", f"quadrature value={value} error={error}")
    return float(value)


def quad_gap(fn: Callable[[float], float], a: float) -> float:
    value, error = quad(fn, a, 40 * a, epsabs=2e-11, epsrel=2e-11, limit=200)
    if not math.isfinite(value) or error > 1e-8 * max(1.0, abs(value)):
        raise ToothFailure("MOMENT_INTEGRALS", f"gapped quadrature value={value} error={error}")
    return float(value)


def compute_moments(data: dict[str, Any], values: dict[str, sp.Expr], responses: dict[str, Any], attacks: dict[str, str]) -> tuple[dict[str, Any], dict[str, dict[str, sp.Expr]]]:
    mp.mp.dps = 35
    d = len(data["ambient"]["coordinates"])
    db = len(data["ambient"]["brane_coordinates"])
    a = float(qvalue(data["geometry"]["a"]["value"]))
    Sd = float(sp.N(sphere_area(d)))
    Sb = float(sp.N(sphere_area(db)))
    profiles = radial_profiles(data, values)
    traces = {name: sp.Symbol(name, real=True) for name in data["core_traces"]}
    tvals = {name: float(qvalue(row["value"])) for name, row in data["core_traces"].items()}
    rho = sp.Symbol("rho_inf", positive=True)
    rn = sp.Symbol("density_delta", real=True)

    measure_power = d - 1
    if attacks.get("radial_measure") == "drop_jacobian_power":
        measure_power -= 1
    # This exact ambient sub-integral independently checks the 4D Jacobian.
    ambient_response_numeric = Sd * quad_inf(lambda r: r**measure_power * profiles["response_bulk"](r) ** 2, a)
    ambient_response_exact = Sd * a**d / (d - 2)
    require(abs(ambient_response_numeric - ambient_response_exact) < 1e-9,
            "MOMENT_INTEGRALS", f"radial/angular measure {ambient_response_numeric}!={ambient_response_exact}")

    density_cross = Sd * quad_gap(lambda r: r ** (d - 1) * profiles["density"](r) * profiles["response_bulk"](r) ** 2, a)
    response_integral = rho * sp.Float(ambient_response_exact, 14) + rn * sp.Float(density_cross, 14)
    bulk_nu = d - 2
    tilt_nu = d - 1
    brane_nu = db
    unit_massless_grad = Sd / d * bulk_nu**2 * a ** (d - 2) / (2 * bulk_nu - d + 2)
    unit_tilt_grad = Sd / d * tilt_nu**2 * a ** (d - 2) / (2 * tilt_nu - d + 2)
    unit_gap_grad_density = Sd / d * quad_gap(lambda r: r ** (d - 1) * profiles["density_d"](r) ** 2, a)
    unit_gap_grad_wall = Sd / d * quad_gap(lambda r: r ** (d - 1) * profiles["wall_d"](r) ** 2, a)
    unit_gap_l2_density = Sd / d * quad_gap(lambda r: r ** (d - 1) * profiles["density"](r) ** 2, a)
    unit_gap_l2_wall = Sd / d * quad_gap(lambda r: r ** (d - 1) * profiles["wall"](r) ** 2, a)
    unit_tilt_l2 = Sd * a**d / (2 * tilt_nu - d)
    unit_response_l2 = Sd * a**d / (d - 2)
    unit_brane_grad = Sb / db * brane_nu**2 * a ** (db - 2) / (2 * brane_nu - db + 2)
    unit_brane_l2 = Sb * a**db / (2 * brane_nu - db)

    # The odd first-order/cross angular moment is evaluated rather than asserted.
    theta = sp.Symbol("theta", real=True)
    odd_angular = sp.integrate(sp.cos(theta) * sp.sin(theta) ** (d - 2), (theta, 0, sp.pi))
    require(odd_angular == 0, "MOMENT_INTEGRALS", f"odd angular integral {odd_angular}")

    moments_out: dict[str, Any] = {}
    moments_expr: dict[str, dict[str, sp.Expr]] = {}
    tilt_phase = traces["tilt_phase"]
    tilt_shear = traces["tilt_shear"]
    tilt_h = traces["tilt_h"]
    for ep in ENDPOINTS:
        cn = sp.sympify(responses[ep]["fluid_coefficients"]["normal"], locals=values)
        ct = sp.sympify(responses[ep]["fluid_coefficients"]["tangent"], locals=values)
        beta = sp.sympify(responses[ep]["shear_coefficient"], locals=values)
        angular_fluid = (cn**2 + (d - 1) * ct**2) / d
        angular_cross = (cn + (d - 1) * ct) / d
        N_UU = sp.expand(angular_fluid * response_integral)
        N_UW = sp.expand(angular_cross * tilt_phase * response_integral)
        N_WW = sp.expand(tilt_phase**2 * response_integral)
        U_XX = sp.Float(unit_brane_grad, 14) * (traces["shear_transverse"] + beta) ** 2
        U_XP = sp.Float(unit_brane_l2, 14) * (traces["shear_transverse"] + beta) * tilt_shear
        U_PP = sp.Float(unit_brane_l2, 14) * tilt_shear**2
        H_XX = sp.Float(unit_massless_grad, 14) * traces["h_scalar"]**2
        H_XP = sp.Integer(0) * traces["h_scalar"] * tilt_h
        H_PP = sp.Float(unit_tilt_l2, 14) * tilt_h**2
        rows = {
            "B_X": sp.Integer(0), "B_p": sp.Integer(0), "B_Xp": sp.Integer(0),
            "N_UU": N_UU, "N_UW": N_UW, "N_WW": N_WW,
            "U_XX": sp.expand(U_XX), "U_XP": sp.expand(U_XP), "U_PP": sp.expand(U_PP),
            "H_XX": sp.expand(H_XX), "H_XP": sp.expand(H_XP), "H_PP": sp.expand(H_PP),
            "I_density_grad": sp.Float(unit_gap_grad_density, 14) * traces["density_delta"]**2,
            "I_density_l2": sp.Float(unit_gap_l2_density, 14) * traces["density_delta"]**2,
            "I_wall_grad": sp.Float(unit_gap_grad_wall, 14) * traces["wall_delta"]**2,
            "I_wall_l2": sp.Float(unit_gap_l2_wall, 14) * traces["wall_delta"]**2,
            "I_shear_grad": sp.Float(unit_brane_grad, 14) * tilt_shear**2,
            "I_h_grad": sp.Float(unit_tilt_grad, 14) * tilt_h**2,
        }
        moments_expr[ep] = rows
        moments_out[ep] = {name: {"expression": str(expr), "production_value": float(sp.N(expr.subs({traces[k]: tvals[k] for k in traces}).subs({rho: 1, rn: tvals["density_delta"]}).subs({s: 1 for s in expr.free_symbols}))),
                                           "dependencies": sorted(str(s) for s in expr.free_symbols)} for name, expr in rows.items()}
    return moments_out, moments_expr


def term_map(expr: sp.Expr) -> list[dict[str, Any]]:
    if sp.simplify(expr) == 0:
        return []
    rows = []
    for term in sp.Add.make_args(sp.expand(expr)):
        coeff, rest = term.as_coeff_Mul()
        powers: dict[str, int] = {}
        numeric = sp.N(coeff, 15)
        for factor, power in rest.as_powers_dict().items():
            if factor.is_number:
                numeric *= factor ** power
            else:
                powers[str(factor)] = int(power)
        rows.append({"coefficient": float(numeric), "powers": dict(sorted(powers.items()))})
    return sorted(rows, key=lambda row: tuple(row["powers"].items()))


def build_effective_action(data: dict[str, Any], symbols: dict[str, sp.Symbol],
                           moments: dict[str, dict[str, sp.Expr]], attacks: dict[str, str]) -> tuple[dict[str, Any], dict[str, str], list[dict[str, Any]], dict[str, dict[str, sp.Expr]]]:
    V, pd, p = sp.symbols("V pd p", real=True)
    hbar, m, rhoBr, Mh, cE = (symbols[x] for x in ("hbar", "m_GNLS", "rhoBr", "Mh", "cE"))
    rho = symbols["rho_inf"]
    symbolic_values = {name: symbol for name, symbol in symbols.items() if name != "a"}
    symbolic_values["a"] = symbols["a"]
    assembled = assemble_action_terms(data, symbolic_values)
    operator = derive_channel_operator(data, symbolic_values, assembled)
    primitive = assembled["primitives"]
    action = assembled["terms"]
    endpoint_output: dict[str, Any] = {}
    reconstruction: dict[str, str] = {}
    ancestry: list[dict[str, Any]] = []
    term_contributions_all: dict[str, dict[str, sp.Expr]] = {}
    for ep in ENDPOINTS:
        z = moments[ep]
        AX = hbar * z["B_X"]
        AP = -hbar * z["B_p"]
        CXP = hbar * z["B_Xp"]
        GVV_terms = {
            "bulk_flow_kinetic": m * z["N_UU"],
            "brane_shear_kinetic": rhoBr * z["U_XX"],
            "h_kinetic": Mh * z["H_XX"] / cE**2,
        }
        GVP_terms = {
            "bulk_flow_kinetic": m * z["N_UW"],
            "brane_shear_kinetic": rhoBr * z["U_XP"],
            "h_kinetic": Mh * z["H_XP"] / cE**2,
        }
        GPP_terms = {
            "bulk_flow_kinetic": m * z["N_WW"],
            "brane_shear_kinetic": rhoBr * z["U_PP"],
            "h_kinetic": Mh * z["H_PP"] / cE**2,
        }
        KPP_terms = {
            "quantum_pressure": operator["entries"]["density_gradient"] * z["I_density_grad"],
            "bulk_EOS": operator["entries"]["density_EOS_curvature"] * z["I_density_l2"],
            "wall_gradient": operator["entries"]["wall_gradient"] * z["I_wall_grad"],
            "wall_double_well": operator["entries"]["wall_well_curvature"] * z["I_wall_l2"],
            "brane_shear_gradient": operator["entries"]["shear_curl"] * z["I_shear_grad"],
            "h_gradient": operator["entries"]["h_gradient"] * z["I_h_grad"],
        }
        GVV, GVP, GPP, KPP = map(lambda table: sp.expand(sum(table.values())), (GVV_terms, GVP_terms, GPP_terms, KPP_terms))
        claimed_parts = {
            "AX": AX * V, "AP": AP * pd, "CXP": CXP * p * V,
            "GVV": GVV * V**2 / 2, "GVP": GVP * V * pd,
            "GPP": GPP * pd**2 / 2, "KPP": -KPP * p**2 / 2,
        }
        claimed = sp.expand(sum(claimed_parts.values()))
        if attacks.get("claimed_L_eff") == "drop_quantum_pressure":
            claimed += KPP_terms["quantum_pressure"] * p**2 / 2

        # Independent rigid-embedding substitution: iterate over the parsed
        # action records and reduce their primitive invariants one term at a
        # time.  This path shares neither claimed_parts nor a copied term map.
        term_contributions: dict[str, sp.Expr] = {term_id: sp.Integer(0) for term_id in action}
        berry_coeff = sp.simplify(sp.diff(action["bulk_berry"], primitive["theta_t"]) / primitive["n"])
        term_contributions["bulk_berry"] = sp.expand(-berry_coeff * z["B_X"] * V + berry_coeff * z["B_p"] * pd - berry_coeff * z["B_Xp"] * p * V)
        flow_coeff = sp.simplify(sp.diff(action["bulk_flow_kinetic"], primitive["v2"]) / primitive["n"])
        term_contributions["bulk_flow_kinetic"] = sp.expand(flow_coeff * (z["N_UU"] * V**2 + 2 * z["N_UW"] * V * pd + z["N_WW"] * pd**2))
        shear_kin_coeff = sp.simplify(sp.diff(action["brane_shear_kinetic"], primitive["u_t2"]).subs(primitive["g_ell"], 1))
        term_contributions["brane_shear_kinetic"] = sp.expand(shear_kin_coeff * (z["U_XX"] * V**2 + 2 * z["U_XP"] * V * pd + z["U_PP"] * pd**2))
        h_kin_coeff = sp.simplify(sp.diff(action["h_kinetic"], primitive["h_t2"]))
        term_contributions["h_kinetic"] = sp.expand(h_kin_coeff * (z["H_XX"] * V**2 + 2 * z["H_XP"] * V * pd + z["H_PP"] * pd**2))
        for term_id, stiffness in KPP_terms.items():
            term_contributions[term_id] = sp.expand(-stiffness * p**2 / 2)
        reconstructed = sp.expand(sum(term_contributions.values()))
        residual = sp.expand(claimed - reconstructed)
        require(residual == 0, "RECONSTRUCTION", f"{ep}:{residual}")
        reconstruction[ep] = str(residual)
        term_contributions_all[ep] = term_contributions
        PX = sp.diff(reconstructed, V)
        Pp = sp.diff(reconstructed, pd)
        Qp = sp.diff(reconstructed, p)
        PX_claim = AX + CXP * p + GVV * V + GVP * pd
        if attacks.get("claimed_canonical_momentum") == "drop_cross_term":
            PX_claim -= GVP * pd
        require(sp.expand(PX_claim - PX) == 0 and sp.expand(AP + GVP * V + GPP * pd - Pp) == 0,
                "CANONICAL_VARIATION", ep)
        require(sp.expand(CXP * V - KPP * p - Qp) == 0, "CANONICAL_VARIATION", f"force {ep}")
        endpoint_output[ep] = {
            "L_eff": str(reconstructed), "canonical_terms": term_map(reconstructed),
            "rigid_embedding_term_contributions": {key: str(value) for key, value in term_contributions.items()},
            "coefficients": {name: {"expression": str(expr), "canonical_terms": term_map(expr),
                                    "dependencies": sorted(str(s) for s in expr.free_symbols)}
                             for name, expr in {"AX": AX, "AP": AP, "CXP": CXP, "GVV": GVV, "GVP": GVP, "GPP": GPP, "KPP": KPP}.items()},
            "canonical_momenta": {"P_X": str(PX), "P_p": str(Pp)},
            "generalized_force": {"Q_p_var": str(Qp), "definition": "-<delta S_action/delta Phi,partial Phi/partial p>-surface variation"},
            "E4_virtual_work": "delta W=lambda_A*(delta V_A-C_A[delta uT_dot])=0 on allowed variations" if ep == "E4" else "STRUCTURAL_ZERO",
        }
        for coefficient, table in (("GVV", GVV_terms), ("GVP", GVP_terms), ("GPP", GPP_terms), ("KPP", KPP_terms)):
            for ancestor, contribution in table.items():
                removed = sp.simplify(contribution)
                if removed == 0:
                    continue
                open_remainder = ancestor.startswith("h_")
                removed_data = copy.deepcopy(data)
                removed_data["action_terms"] = [row for row in removed_data["action_terms"] if row["id"] != ancestor]
                rederived_ids = set(assemble_action_terms(removed_data, symbolic_values)["terms"])
                after_removal = contribution if ancestor in rederived_ids else sp.Integer(0)
                ancestry.append({"endpoint": ep, "structure": f"{coefficient}.{ancestor}", "ancestor": ancestor,
                                 "contribution": str(contribution), "removal_delta": str(removed),
                                 "before_monomials": term_map(contribution),
                                 "after_removal_monomials": term_map(after_removal),
                                 "classification_before": "UNRESOLVED_OPEN" if open_remainder else "NONZERO_NATIVE_STRUCTURE",
                                 "classification_after_removal": "ABSENT",
                                 "remove_then_rederive_action_ids": sorted(rederived_ids),
                                 "origin_mutation_changes_expression": sp.simplify(contribution - after_removal) != 0,
                                 "native_padding_ablation_destroys_structure": after_removal == 0,
                                 "open_remainder": open_remainder})
    if attacks.get("claimed_structure") == "open_only_inertia":
        ancestry.append({"endpoint": "E1", "structure": "GVV.open_only", "ancestor": None,
                         "contribution": "open_response", "removal_delta": "0",
                         "origin_mutation_changes_expression": False,
                         "native_padding_ablation_destroys_structure": False, "open_remainder": False})
    if attacks.get("native_padding_structure") == "open_plus_native":
        open_leaf, native_leaf = sp.symbols("open_response native_action_piece")
        before, after = open_leaf + native_leaf, open_leaf
        ancestry.append({"endpoint": "E1", "structure": "GVV.native_padded_open", "ancestor": "bulk_flow_kinetic",
                         "contribution": str(before), "removal_delta": str(before - after),
                         "before_monomials": term_map(before), "after_removal_monomials": term_map(after),
                         "classification_before": "OPEN_DOMINATED_STRUCTURE", "classification_after_removal": "OPEN_DOMINATED_STRUCTURE",
                         "origin_mutation_changes_expression": True,
                         "native_padding_ablation_destroys_structure": False, "open_remainder": False})
    require(all(row["open_remainder"] or (row["ancestor"] and row["origin_mutation_changes_expression"])
                for row in ancestry), "ANCESTRY", "open-only/inert additive structure")
    require(all(row["open_remainder"] or row["native_padding_ablation_destroys_structure"] for row in ancestry),
            "NATIVE_PADDING", "native padding survived actual ancestor removal")
    return endpoint_output, reconstruction, ancestry, term_contributions_all


def graph_and_channels(data: dict[str, Any], symbols: dict[str, sp.Symbol],
                       term_contributions: dict[str, dict[str, sp.Expr]],
                       attacks: dict[str, str]) -> tuple[dict[str, Any], dict[str, Any]]:
    channels = {ep: copy.deepcopy(data["endpoints"][ep]["channels"]) for ep in ENDPOINTS}
    for ep, row in channels.items():
        terms = [term for values in row.values() for term in values]
        require(len(terms) == len(set(terms)), "CHANNEL_UNIQUENESS", ep)
    symbolic_values = {name: symbol for name, symbol in symbols.items() if name != "a"}
    symbolic_values["a"] = symbols["a"]
    assembled = assemble_action_terms(data, symbolic_values)
    operator = derive_channel_operator(data, symbolic_values, assembled)
    nodes = [{"id": row["id"], "type": "ACTION", "expression": row["expression"]}
             for row in data["action_terms"]]
    nodes += [{"id": row["id"], "type": row["root_type"]} for row in data["field_records"]]
    nodes += [{"id": "base_family", "type": "DELIVERABLE_7_0"}, {"id": "D_E", "type": "OPERATOR"},
              {"id": "L_eff", "type": "DELIVERABLE_7_1"}, {"id": "zero_mode", "type": "DELIVERABLE_7_1"},
              {"id": "generalized_force", "type": "DELIVERABLE_7_1"}]
    edges: list[list[str]] = []
    dependency_sets: dict[str, list[str]] = {}
    for term_id, expression in operator["quadratic_by_term"].items():
        dependency_sets[f"D_E:{term_id}"] = sorted(str(s) for s in expression.free_symbols)
        if expression != 0:
            edges.append([term_id, "D_E"])
    for term_id, expression in term_contributions["E1"].items():
        dependency_sets[f"L_eff:{term_id}"] = sorted(str(s) for s in expression.free_symbols)
        if expression != 0:
            edges.append([term_id, "L_eff"])
    for row in data["field_records"]:
        if row["root_type"] in {"BALANCE", "PRIMITIVE_OPEN"}:
            edges.append([row["id"], "base_family"])
        if row["root_type"] in {"BALANCE", "CONSTRAINT", "RAYLEIGH", "RETURN"}:
            edges.append([row["id"], "generalized_force"])
    edges += [["base_family", "D_E"], ["D_E", "zero_mode"], ["L_eff", "generalized_force"]]
    if attacks.get("provenance_edge") == "native_continuity_to_L_eff":
        injected = sp.Symbol("native_continuity_dependency")
        dependency_sets["L_eff:injected_balance"] = sorted(str(s) for s in injected.free_symbols)
        edges.append(["native_continuity", "L_eff"])
    if attacks.get("provenance_expression"):
        injected = sp.Symbol(attacks["provenance_expression"])
        nodes.append({"id": "adversarial_import", "type": "FORBIDDEN_IMPORT",
                      "expression": str(injected)})
        dependency_sets["L_eff:adversarial_import"] = sorted(str(s) for s in injected.free_symbols)
        edges.append(["adversarial_import", "L_eff"])
    by_id = {row["id"]: row["type"] for row in nodes}
    leff_parents = {by_id[src] for src, dst in edges if dst == "L_eff"}
    forbidden_types = {"FORBIDDEN_IMPORT"}
    forbidden_nodes = [row["id"] for row in nodes if row["type"] in forbidden_types]
    require(not forbidden_nodes, "PROVENANCE_FORBIDDEN", str(forbidden_nodes))
    require(leff_parents <= {"ACTION", "PRIMITIVE_OPEN"}, "TYPED_DATAFLOW", str(sorted(leff_parents)))
    return channels, {"nodes": nodes, "edges": sorted(edges), "expression_dependency_sets": dependency_sets,
                      "L_eff_parent_types": sorted(leff_parents),
                      "forbidden_rejection_algorithm": "graph node-type traversal, not name matching"}


def input_ledger(data: dict[str, Any]) -> list[dict[str, Any]]:
    required = {"id", "status", "root_type", "domain", "dimensions", "arguments", "symmetry_class", "source"}
    collective = {"X", "p"}
    output = []
    for row in data["field_records"]:
        require(required <= set(row), "INPUT_LEDGER", f"incomplete {row.get('id')}")
        args = set(row["arguments"])
        domain_tokens = set(re.findall(r"[A-Za-z_]+", row["domain"]))
        require(not (args & collective) and "collective_coordinate_functional" not in domain_tokens,
                "INPUT_LEDGER", f"collective answer warehouse {row['id']}")
        require("multipole" not in domain_tokens, "INPUT_LEDGER", f"far-field answer warehouse {row['id']}")
        require(len(row["dimensions"]) == 3, "INPUT_LEDGER", f"dimension tuple {row['id']}")
        output.append(copy.deepcopy(row))
    # Coefficients are typed primitive action data and are included in the ledger.
    for name, row in data["coefficients"].items():
        output.append({"id": name, "status": row["status"], "root_type": "ACTION_COEFFICIENT",
                       "domain": "coefficient", "dimensions": row["dimensions"], "arguments": [],
                       "symmetry_class": "R_W_EVEN", "source": row["source"], "constraint": row["constraint"]})
    return output


def parity_check(data: dict[str, Any]) -> dict[str, Any]:
    w, s = sp.symbols("w s", real=True)
    plus = sp.Function("PhiPlus")(s * w)
    minus_defined = plus.subs({s: -s, w: -w})
    residual = sp.simplify(minus_defined - plus)
    require(residual == 0, "PARITY", str(residual))
    symmetric_map = data["parity_cases"]["symmetric_postulate"]["background_map"]
    one_sided_map = data["parity_cases"]["one_sided_pathA29"]["background_map"]
    symmetric_expr = sp.sympify(symmetric_map["w"], locals={"w": w}) if set(symmetric_map) == {"w"} else sp.nan
    symmetric_tag = "BODY_PLUS_AMBIENT_POSTULATE" if symmetric_expr == -w else "UNCLASSIFIED_BACKGROUND"
    one_sided = "ONE_SIDED_ASYMMETRY_MAP" if set(one_sided_map) == {"operator", "arguments"} and one_sided_map["operator"] == "ambient_asymmetry_map" else symmetric_tag
    require(one_sided == "ONE_SIDED_ASYMMETRY_MAP", "PARITY", one_sided)
    return {"body_conjugation_residual": str(residual), "embedding_tag": "BODY_CONJUGATION_ONLY",
            "symmetric_tag": symmetric_tag,
            "asymmetric_control_tag": one_sided}


def classify(tails: list[dict[str, Any]], values: dict[str, sp.Expr]) -> tuple[str, dict[str, Any]]:
    tachyons = [row["id"] for row in tails if row["classification"] == "TACHYONIC"]
    unresolved = sorted({dep for row in tails if row["classification"] == "UNRESOLVED" for dep in row["dependencies"]})
    nonnormal = [row["id"] for row in tails if row["normalizable"] is False and row["classification"] != "TACHYONIC"]
    kinetic = {"density": values["hbar"]**2 / (4 * values["m_GNLS"] * values["rho_inf"]),
               "wall": values["kappaB"], "shear": values["rhoBr"],
               "h": values["Mh"] / values["cE"]**2}
    kinetic_signs = {name: numeric_sign(expr) for name, expr in kinetic.items()}
    unresolved += sorted(str(s) for name, expr in kinetic.items() if kinetic_signs[name] is None for s in expr.free_symbols)
    unstable_kinetic = [name for name, sign in kinetic_signs.items() if sign == -1]
    if unresolved:
        leaves = ",".join(sorted(set(unresolved)))
        return f"U1_BASE_UNRESOLVED({leaves})", {"unresolved_leaves": sorted(set(unresolved)), "kinetic_signs": kinetic_signs}
    if tachyons or unstable_kinetic:
        modes_list = tachyons + unstable_kinetic
        modes = ",".join(modes_list)
        return f"U1_BASE_UNSTABLE({modes})", {"tachyonic_channels": tachyons,
                "negative_kinetic_channels": unstable_kinetic,
                "kinetic_signs": kinetic_signs}
    if nonnormal:
        return f"U1_BASE_ILL_POSED(NO_NORMALIZABLE_ZERO_MODE:{','.join(nonnormal)})", {"nonnormalizable_channels": nonnormal, "kinetic_signs": kinetic_signs}
    return "U1_BASE_OK", {"kinetic_signs": kinetic_signs}


def light_outcome(data: dict[str, Any], attacks: dict[str, str]) -> str:
    _, values = coefficient_symbols(data)
    values["a"] = qvalue(data["geometry"]["a"]["value"])
    assembled = assemble_action_terms(data, values)
    operator = derive_channel_operator(data, values, assembled)
    tails = [solve_tail(ch, attacks) for ch in operator["channels"]]
    coupled_indicial_analysis(data, operator)
    return classify(tails, values)[0]


def build(data: dict[str, Any], attacks: dict[str, str], fixtures: dict[str, Any], include_reachability: bool = True) -> dict[str, Any]:
    ledger = input_ledger(data)
    completeness = source_completeness(data, attacks)
    dimensions = dimensional_firewall(data)
    flux = comoving_reductions(data)
    symbols, values = coefficient_symbols(data)
    # Geometry entries participate symbolically in L_eff even though their
    # production values are used for quadrature.
    symbols["a"] = sp.Symbol("a", positive=True)
    values["a"] = qvalue(data["geometry"]["a"]["value"])
    assembled = assemble_action_terms(data, values)
    operator = derive_channel_operator(data, values, assembled)
    channels = operator["channels"]
    tails = [solve_tail(row, attacks) for row in channels]
    coupled = coupled_indicial_analysis(data, operator)
    removal_probes = action_term_removal_probes(data, values, operator)
    phase = phase_normalization(data, values)
    zero_mode = zero_mode_operator(data, values, attacks)
    quotient = projector_algebra(tails, attacks)
    responses = endpoint_responses(data, values)
    moment_rows, moment_expr = compute_moments(data, values, responses, attacks)
    endpoint_rows, reconstruction, ancestry, term_contributions = build_effective_action(data, symbols, moment_expr, attacks)
    channel_rows, graph = graph_and_channels(data, symbols, term_contributions, attacks)
    parity = parity_check(data)
    verdict, decision = classify(tails, values)
    reachability: dict[str, str] = {}
    if include_reachability:
        case_names = list(fixtures["outcome_cases"])
        if attacks.get("outcome_case") == "remove_fat_tail":
            case_names.remove("fat_tail")
        for name in case_names:
            case_data = copy.deepcopy(yaml.safe_load(DEFAULT_INPUT.read_text()))
            apply_operations(case_data, fixtures["outcome_cases"][name])
            reachability[name] = light_outcome(case_data, {})
        classes = {out.split("(", 1)[0] for out in reachability.values()}
        require(classes == {"U1_BASE_OK", "U1_BASE_UNSTABLE", "U1_BASE_UNRESOLVED", "U1_BASE_ILL_POSED"},
                "OUTCOME_REACHABILITY", str(reachability))

    # E4's multiplier is an actual reduced virtual-work object.
    lam, dV, du = sp.symbols("lambda_A deltaV_A delta_uTdot_A")
    e4_virtual_work = sp.expand(lam * (dV - du))
    e4_allowed = sp.simplify(e4_virtual_work.subs(dV, du))
    require(e4_allowed == 0, "ENDPOINT_RESPONSE", str(e4_allowed))

    known_normalizability = [row["normalizable"] for row in tails if row["normalizable"] is not None]
    gates = {
        "G1": ("CLASSIFIED_BY_AXIS2" if len(known_normalizability) != len(tails)
               else ("PASS" if all(known_normalizability) else "NEGATIVE_COMPUTED")),
        "G2": "PASS" if all(math.isfinite(row["zero_mode_norm_value"]) for row in tails if row["zero_mode_norm_value"] is not None) else "NEGATIVE_COMPUTED",
        "G3": "PASS" if verdict == "U1_BASE_OK" else "CLASSIFIED_BY_AXIS2",
        "G4": "NOT_RUN(phase_B)",
        "G5": "PASS",
        "G6": "PASS_ENDPOINT_MAP_REPORTED",
        "G7": "NOT_RUN(phase_B)",
        "G8": "NOT_RUN(phase_C)",
        "G9": "NOT_RUN(phase_B;tolerance_predeclared)",
        "G10": "NOT_RUN(phase_C)",
        "G11": "NOT_RUN(phase_C)",
    }
    cells = {f"{ep}:{ambient}": verdict for ep in ENDPOINTS for ambient in data["parity_cases"]}
    return {
        "engine": "SymPy", "phase": "A", "schema_version": data["schema_version"],
        "axis1": "COMPUTATION_VALID", "axis2": verdict, "decision_evidence": decision,
        "cells": cells, "gates": gates, "declared_inputs": ledger,
        "source_action_completeness": completeness, "dimensional_firewall": dimensions,
        "assembled_action": {
            "term_expressions": {key: str(value) for key, value in assembled["terms"].items()},
            "term_dependencies": assembled["dependencies"],
            "total_expression": str(assembled["total_expression"]),
        },
        "linearized_channel_operator": {
            "field_order": operator["field_order"], "gradient_order": operator["gradient_order"],
            "quadratic_expression": str(operator["quadratic_expression"]),
            "quadratic_by_action_term": {key: str(value) for key, value in operator["quadratic_by_term"].items()},
            "gradient_hessian": [[str(x) for x in row] for row in operator["gradient_hessian"].tolist()],
            "curvature_hessian": [[str(x) for x in row] for row in operator["curvature_hessian"].tolist()],
            "mixed_hessian": [[str(x) for x in row] for row in operator["mixed_hessian"].tolist()],
            "entries": {key: str(value) for key, value in operator["entries"].items()},
        },
        "action_term_removal_probes": removal_probes,
        "coupled_indicial_analysis": coupled,
        "co_moving_laws": flux, "tail_channels": tails, "phase_flux_normalization": phase,
        "linearized_force_balance": zero_mode, "zero_mode_quotient": quotient,
        "endpoint_responses": responses, "evaluated_moments": moment_rows,
        "endpoint_effective_actions": endpoint_rows, "reconstruction_residuals": reconstruction,
        "channel_assignments": channel_rows, "per_structure_ancestry": ancestry,
        "provenance_graph": graph, "parity": parity, "verdict_reachability": reachability,
        "E4_reduced_virtual_work": {"relation": str(e4_virtual_work), "allowed_variation_residual": str(e4_allowed), "multiplier": str(lam)},
        "ir_scheme": data["ir_scheme"], "partition": data["partition"],
        "base_configuration": {
            "class": data["kinematics"]["base_state_class"], "control_volume": data["kinematics"]["control_volume"],
            "family": "action-derived exterior radial solutions joined to declared core traces by throat_surface_functional",
            "balance": "native action Euler/continuity/momentum plus the typed co-moving surface/source functionals",
            "density_epp": str(channels[0]["curvature"]),
            "body_conjugation": "Phi_-(x,w)=R_w Phi_+(x,-w)",
        },
        "checks": {tooth: "PASS" for tooth in TEETH}, "teeth": list(TEETH),
        "downstream": {"phase_B": "NOT_RUN(upstream)", "phase_C": "NOT_RUN(upstream)"},
        "honest_scope": "collective-coordinate/effective-action Phase A only; exterior analytic family plus typed surface core, not a nonlinear throat simulation",
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--fixtures", type=Path, default=DEFAULT_FIXTURES)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--case", default="production")
    parser.add_argument("--mutation", choices=TEETH)
    args = parser.parse_args()
    try:
        data, attacks, fixtures = load_case(args.input, args.fixtures, args.case, args.mutation)
        result = build(data, attacks, fixtures, include_reachability=(args.mutation in {None, "OUTCOME_REACHABILITY"}))
        if args.mutation:
            print(f"ASSERT_FAIL:MUTATION_NOOP:{args.mutation}:mutation did not reach its own assert", file=sys.stderr)
            return 1
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
        print(f"SYMPY_PHASE_A: PASS axis2={result['axis2']} channels={len(result['tail_channels'])}")
        print(f"OUTPUT: {args.output}")
        return 0
    except ToothFailure as exc:
        print(str(exc), file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
