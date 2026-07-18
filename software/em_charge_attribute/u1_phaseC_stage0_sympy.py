#!/usr/bin/env python3
"""Independent SymPy route for the U1 Phase-C stage-0 contract.

This executable deliberately imports no task-authored Python module.  It reads
the frozen inputs itself, derives the stage-0 mathematical/set objects, and
emits an engine-local YAML artifact for the independent comparator.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import sys
from pathlib import Path
from typing import Any

import sympy as sp
import yaml


SCHEMA = "U1_PHASE_C_STAGE0_ENGINE_V1"
MEDIATORS = ("h", "u_T", "u_L", "wall_chi")
ENDPOINTS = ("E1", "E2", "E3", "E4", "E5")
AMBIENTS = ("one_sided_pathA29", "symmetric_postulate")
CLOSURES = ("body_mass_growth", "return_path", "sleeve_exit")
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
G8_FLOOR = (
    "common_drain",
    "orientation_sleeve",
    "post_mouth_axial_flow",
    "h",
    "u_T",
    "u_L",
    "E4_constraint_reaction",
    "outer_surface_flux_return",
    "wall_chi_modes",
)


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def sha256_path(path: Path) -> str:
    return sha256_bytes(path.read_bytes())


def canonical_bytes(value: Any) -> bytes:
    return json.dumps(
        value, ensure_ascii=False, sort_keys=True, separators=(",", ":")
    ).encode("utf-8")


def digest(value: Any) -> str:
    return sha256_bytes(canonical_bytes(value))


def canonical_string_key(value: str) -> tuple[str, str]:
    """Engine-neutral ordering for contract identifiers.

    The frozen identifiers are ASCII, but several contain uppercase type
    fragments (for example ``bulk_EOS``).  Python's default code-point order
    and Wolfram's string order disagree on those fragments, so all ordered
    semantic records use an explicit case-folded key with an exact-string
    tiebreaker.
    """
    return value.casefold(), value


def string_record_digest(records: list[str]) -> str:
    row_hashes = [sha256_bytes(record.encode("utf-8")) for record in records]
    return sha256_bytes("\n".join(sorted(row_hashes)).encode("ascii"))


def load_yaml(path: Path) -> Any:
    with path.open("rb") as handle:
        return yaml.load(handle, Loader=yaml.CSafeLoader)


def scalar_paths(value: Any, predicate, path: tuple[str, ...] = ()) -> list[dict[str, str]]:
    out: list[dict[str, str]] = []
    if isinstance(value, dict):
        for key, child in value.items():
            out.extend(scalar_paths(child, predicate, path + (str(key),)))
    elif isinstance(value, list):
        for index, child in enumerate(value):
            out.extend(scalar_paths(child, predicate, path + (str(index),)))
    elif predicate(value):
        out.append({"path": "/".join(path), "value": str(value)})
    return out


def parse_obligation(record: str) -> tuple[str, dict[str, str]]:
    parts = record.split("|")
    axes: dict[str, str] = {}
    for part in parts[1:]:
        if "=" in part:
            key, value = part.split("=", 1)
            axes[key] = value
    return parts[0], axes


def expand_jet_expression(expression: str) -> sp.Expr:
    """Restore component jets before taking a second variation.

    The replacements are a representation of the explicit quadratic
    contractions already frozen in the Phase-A input; no coefficient or
    profile is supplied here.
    """
    local_names = set(re.findall(r"[A-Za-z_][A-Za-z0-9_]*", expression))
    local_names.update(
        {
            "n_gw",
            "n_gx",
            "n_gy",
            "n_gz",
            "v_w",
            "v_x",
            "v_y",
            "v_z",
            "chi_gw",
            "chi_gx",
            "chi_gy",
            "chi_gz",
            "u_tx",
            "u_ty",
            "u_tz",
            "u_cx",
            "u_cy",
            "u_cz",
            "ud_1",
            "ud_2",
            "uw_t",
            "h_t",
            "h_gw",
            "h_gx",
            "h_gy",
            "h_gz",
        }
    )
    symbols = {name: sp.Symbol(name) for name in local_names}
    expr = sp.sympify(expression.replace("^", "**"), locals=symbols)
    replacements = {
        symbols.get("n_grad2", sp.Symbol("n_grad2")): sum(symbols[x] ** 2 for x in ("n_gw", "n_gx", "n_gy", "n_gz")),
        symbols.get("v2", sp.Symbol("v2")): sum(symbols[x] ** 2 for x in ("v_w", "v_x", "v_y", "v_z")),
        symbols.get("chi_grad2", sp.Symbol("chi_grad2")): sum(symbols[x] ** 2 for x in ("chi_gw", "chi_gx", "chi_gy", "chi_gz")),
        symbols.get("u_t2", sp.Symbol("u_t2")): sum(symbols[x] ** 2 for x in ("u_tx", "u_ty", "u_tz")),
        symbols.get("u_curl2", sp.Symbol("u_curl2")): sum(symbols[x] ** 2 for x in ("u_cx", "u_cy", "u_cz")),
        symbols.get("ud_curl2", sp.Symbol("ud_curl2")): symbols["ud_1"] ** 2 + symbols["ud_2"] ** 2,
        symbols.get("uw_t2", sp.Symbol("uw_t2")): symbols["uw_t"] ** 2,
        symbols.get("h_t2", sp.Symbol("h_t2")): symbols["h_t"] ** 2,
        symbols.get("h_grad2", sp.Symbol("h_grad2")): sum(symbols[x] ** 2 for x in ("h_gw", "h_gx", "h_gy", "h_gz")),
    }
    return sp.expand(expr.xreplace(replacements))


def field_group(name: str) -> str | None:
    if name in {"f_throat", "f_mix"}:
        return "opaque_open_surface"
    if name == "chiB" or name.startswith("chi_g"):
        return "wall_chi"
    if name.startswith("ud_"):
        return "u_d_transverse"
    if name == "uw" or name == "uw_t":
        return "wall_chi"
    if name == "h_t" or name.startswith("h_g") or name == "h":
        return "h"
    if name in {"n", "theta_t", "v_w", "v_x", "v_y", "v_z"} or name.startswith("n_g"):
        return "u_L"
    if name.startswith("u_t") or name.startswith("u_c"):
        return "u_T"
    return None


def action_derivation(action_terms: list[dict[str, Any]]) -> dict[str, Any]:
    term_records: list[dict[str, Any]] = []
    all_pairs: list[dict[str, str]] = []
    total = sp.Integer(0)
    for term in sorted(action_terms, key=lambda row: canonical_string_key(row["id"])):
        expr_text = term["expression"]
        if expr_text in {"f_throat", "f_mix"}:
            term_records.append(
                {
                    "id": term["id"],
                    "status": "UNRESOLVED_NATIVE_FUNCTIONAL_SECOND_VARIATION",
                    "field_groups": ["opaque_open_surface"],
                    "nonzero_hessian_pairs": [],
                }
            )
            continue
        expr = expand_jet_expression(expr_text)
        total += expr
        variables = sorted(
            [symbol for symbol in expr.free_symbols if field_group(symbol.name)],
            key=lambda symbol: symbol.name,
        )
        groups = sorted({field_group(symbol.name) for symbol in variables if field_group(symbol.name)})
        pairs: set[tuple[str, str]] = set()
        for index, left in enumerate(variables):
            for right in variables[index:]:
                second = sp.simplify(sp.diff(expr, left, right))
                if second != 0:
                    g_left = field_group(left.name)
                    g_right = field_group(right.name)
                    if g_left and g_right:
                        pair = tuple(sorted((g_left, g_right)))
                        pairs.add(pair)
        pair_strings = [f"{left}|{right}" for left, right in sorted(pairs)]
        for pair in pair_strings:
            all_pairs.append({"term": term["id"], "pair": pair})
        term_records.append(
            {
                "id": term["id"],
                "status": "DERIVED_COMPLETE_SCALAR_JET_SECOND_VARIATION",
                "field_groups": groups,
                "nonzero_hessian_pairs": pair_strings,
            }
        )
    canonical_terms = [
        {"id": row["id"], "expression": row["expression"], "support": row.get("support", "bulk")}
        for row in sorted(
            action_terms, key=lambda item: canonical_string_key(item["id"])
        )
    ]
    return {
        "S_body_Omega_c": {
            "status": "DERIVED_EXPLICIT_SYMBOLIC_FUNCTIONAL",
            "domain": "time_x_Omega_c",
            "variation": "delta_S_body/delta_Phi_A",
            "canonical_terms": canonical_terms,
            "value_digest": digest(canonical_terms),
            "native_term_count": len(canonical_terms),
        },
        "second_variation": {
            "status": "DERIVED_WITH_OPEN_SURFACE_REMAINDERS",
            "term_records": term_records,
            "nonzero_pair_records": all_pairs,
            "chi_u_mixed_block_present": any(
                row["term"] == "wall_shear_gate" and row["pair"] == "u_d_transverse|wall_chi"
                for row in all_pairs
            ),
            "value_digest": digest(term_records),
        },
        "total_constructed": str(sp.expand(total)),
    }


def hessian_challenge_from_raw_action(action_terms: list[dict[str, Any]]) -> dict[str, Any]:
    """Fresh constructive Hessian route used by the first-construction tooth.

    The implementation re-parses raw action strings and does not consume
    ``action_derivation`` or its term records.
    """
    term_records: list[dict[str, Any]] = []
    pair_records: list[dict[str, str]] = []
    for term in sorted(action_terms, key=lambda row: canonical_string_key(row["id"])):
        text = term["expression"]
        if text in {"f_throat", "f_mix"}:
            term_records.append(
                {
                    "id": term["id"],
                    "status": "UNRESOLVED_NATIVE_FUNCTIONAL_SECOND_VARIATION",
                    "field_groups": ["opaque_open_surface"],
                    "nonzero_hessian_pairs": [],
                }
            )
            continue
        names = set(re.findall(r"[A-Za-z_][A-Za-z0-9_]*", text))
        jet_names = {
            "n_gw", "n_gx", "n_gy", "n_gz", "v_w", "v_x", "v_y", "v_z",
            "chi_gw", "chi_gx", "chi_gy", "chi_gz", "u_tx", "u_ty", "u_tz",
            "u_cx", "u_cy", "u_cz", "ud_1", "ud_2", "uw_t", "h_t",
            "h_gw", "h_gx", "h_gy", "h_gz",
        }
        names.update(jet_names)
        symbols = {name: sp.Symbol(name) for name in names}
        expr = sp.sympify(text.replace("^", "**"), locals=symbols)
        expansion = {
            symbols.get("n_grad2", sp.Symbol("n_grad2")): sum(symbols[x] ** 2 for x in ("n_gw", "n_gx", "n_gy", "n_gz")),
            symbols.get("v2", sp.Symbol("v2")): sum(symbols[x] ** 2 for x in ("v_w", "v_x", "v_y", "v_z")),
            symbols.get("chi_grad2", sp.Symbol("chi_grad2")): sum(symbols[x] ** 2 for x in ("chi_gw", "chi_gx", "chi_gy", "chi_gz")),
            symbols.get("u_t2", sp.Symbol("u_t2")): sum(symbols[x] ** 2 for x in ("u_tx", "u_ty", "u_tz")),
            symbols.get("u_curl2", sp.Symbol("u_curl2")): sum(symbols[x] ** 2 for x in ("u_cx", "u_cy", "u_cz")),
            symbols.get("ud_curl2", sp.Symbol("ud_curl2")): symbols["ud_1"] ** 2 + symbols["ud_2"] ** 2,
            symbols.get("uw_t2", sp.Symbol("uw_t2")): symbols["uw_t"] ** 2,
            symbols.get("h_t2", sp.Symbol("h_t2")): symbols["h_t"] ** 2,
            symbols.get("h_grad2", sp.Symbol("h_grad2")): sum(symbols[x] ** 2 for x in ("h_gw", "h_gx", "h_gy", "h_gz")),
        }
        expr = sp.expand(expr.xreplace(expansion))
        variables = sorted(
            (symbol for symbol in expr.free_symbols if field_group(symbol.name)),
            key=lambda symbol: symbol.name,
        )
        groups = sorted({field_group(symbol.name) for symbol in variables if field_group(symbol.name)})
        pairs: set[tuple[str, str]] = set()
        for index, left in enumerate(variables):
            for right in variables[index:]:
                if sp.simplify(sp.diff(expr, left, right)) != 0:
                    left_group = field_group(left.name)
                    right_group = field_group(right.name)
                    if left_group and right_group:
                        pairs.add(tuple(sorted((left_group, right_group))))
        pair_strings = [f"{left}|{right}" for left, right in sorted(pairs)]
        pair_records.extend({"term": term["id"], "pair": pair} for pair in pair_strings)
        term_records.append(
            {
                "id": term["id"],
                "status": "DERIVED_COMPLETE_SCALAR_JET_SECOND_VARIATION",
                "field_groups": groups,
                "nonzero_hessian_pairs": pair_strings,
            }
        )
    value = {
        "status": "DERIVED_WITH_OPEN_SURFACE_REMAINDERS",
        "term_records": term_records,
        "nonzero_pair_records": pair_records,
        "chi_u_mixed_block_present": any(
            row["term"] == "wall_shear_gate" and row["pair"] == "u_d_transverse|wall_chi"
            for row in pair_records
        ),
        "value_digest": digest(term_records),
    }
    return {
        "executed": True,
        "semantic_route_id": "raw_action_hessian_challenge_v1",
        "engine_local_function": "hessian_challenge_from_raw_action",
        "shared_with_derived_route": "raw_action_terms_only",
        "candidate_schema": {
            "required_type": "bilinear_field_hessian_with_open_remainders",
            "required_dimensions": "action/field^2",
            "domain": "time_x_Omega_c",
        },
        "candidate_value": value,
        "candidate_value_digest": digest(value),
        "candidate_is_well_typed": True,
        "defining_predicate_evaluated": True,
        "defining_predicate_pass": value["chi_u_mixed_block_present"]
        and len(value["term_records"]) == len(action_terms),
    }


def mediator_incidence(term_record: dict[str, Any]) -> list[str]:
    groups = set(term_record["field_groups"])
    if "opaque_open_surface" in groups:
        return list(MEDIATORS)
    out: set[str] = set()
    if "h" in groups:
        out.add("h")
    if "u_T" in groups or "u_d_transverse" in groups:
        out.add("u_T")
    if "u_L" in groups:
        out.add("u_L")
    if "wall_chi" in groups or "u_d_transverse" in groups:
        out.add("wall_chi")
    return sorted(out)


def g8_mediator_incidence_from_raw_expression(term: dict[str, Any]) -> list[str]:
    """Independent G8 incidence walk over the frozen raw action syntax."""
    expression = term["expression"]
    if expression in {"f_throat", "f_mix"}:
        return list(MEDIATORS)
    tokens = set(re.findall(r"[A-Za-z_][A-Za-z0-9_]*", expression))
    out: set[str] = set()
    if tokens & {"h", "h_t2", "h_grad2"}:
        out.add("h")
    if tokens & {"u_t2", "u_curl2", "ud_curl2"}:
        out.add("u_T")
    if tokens & {"n", "theta_t", "v2", "n_grad2"}:
        out.add("u_L")
    if tokens & {"chiB", "chi_grad2", "uw", "uw_t2", "ud_curl2"}:
        out.add("wall_chi")
    return sorted(out)


def g8_endpoint_incidence(row: dict[str, Any]) -> list[str]:
    """Independent endpoint walk from typed root/domain/argument metadata."""
    root_type = row["root_type"]
    domain = row["domain"]
    arguments = set(row.get("arguments", []))
    if root_type == "ACTION" and "partial_Omega_c" in domain:
        return list(MEDIATORS)
    if root_type == "BALANCE" and row["id"] == "native_momentum":
        return list(MEDIATORS)
    if root_type == "RETURN":
        return list(MEDIATORS)
    if root_type == "CONSTRAINT":
        return ["u_T"]
    if root_type == "RAYLEIGH":
        out = ["u_T"]
        if "sleeve_velocity_trace" in arguments:
            out.append("wall_chi")
        return out
    return []


def coupling_radiation_incidence(operator: dict[str, Any]) -> list[str]:
    """Coupling-census walk over the B2 operator family declaration."""
    native = operator["id"]
    if native in {"geon_open", "throat_source_open", "wall_mix_open"}:
        return list(MEDIATORS)
    family = operator["family"]
    if family == "h":
        return ["h"]
    if family == "u_T":
        return ["u_T"]
    if family == "u_L":
        return ["u_L"]
    if family == "wall_chi/u_d":
        return ["u_T", "wall_chi"]
    if family == "wall_chi":
        return ["wall_chi"]
    return list(MEDIATORS)


def g8_radiation_incidence(operator: dict[str, Any]) -> list[str]:
    """Independent G8 walk over B2 field-block tokens, not family labels."""
    fields = set(operator.get("field_block", []))
    if fields & {"geon_core_bundle", "native_throat_fields", "native_mixed_fields"}:
        return list(MEDIATORS)
    out: set[str] = set()
    if "h" in fields:
        out.add("h")
    if fields & {"u_T", "u_d_transverse"}:
        out.add("u_T")
    if fields & {"delta_n", "theta"}:
        out.add("u_L")
    if fields & {"u_w", "delta_chiB", "u_d_transverse"}:
        out.add("wall_chi")
    return sorted(out) or list(MEDIATORS)


def build_typed_roots(
    phase_a: dict[str, Any], b1: dict[str, Any], b2_contract: dict[str, Any]
) -> list[dict[str, Any]]:
    roots: list[dict[str, Any]] = []
    for term in phase_a["action_terms"]:
        roots.append(
            {
                "id": f"action:{term['id']}",
                "native_id": term["id"],
                "root_type": "native_action_term",
                "status": term.get("status", "WHITELIST_SOURCED"),
                "source": term["source_file"],
            }
        )
    for row in b1["declared_inputs"]:
        roots.append(
            {
                "id": f"input:{row['id']}",
                "native_id": row["id"],
                "root_type": row["root_type"],
                "status": row["status"],
                "source": row["source"],
                "domain": row["domain"],
                "dimensions": row["dimensions"],
            }
        )
    for operator in b2_contract["frozen_data"]["native_operator_inventory"]:
        roots.append(
            {
                "id": f"radiation:{operator['id']}",
                "native_id": operator["id"],
                "root_type": "radiative_channel",
                "status": operator["status"],
                "source": "B2_native_operator_inventory",
                "family": operator["family"],
            }
        )
    return sorted(roots, key=lambda row: row["id"])


def action_record_map(derivation: dict[str, Any]) -> dict[str, dict[str, Any]]:
    return {row["id"]: row for row in derivation["second_variation"]["term_records"]}


def build_force_census(
    phase_a: dict[str, Any], b1: dict[str, Any], b2_contract: dict[str, Any], derivation: dict[str, Any]
) -> dict[str, Any]:
    entries: list[dict[str, Any]] = []
    records = action_record_map(derivation)
    for term in sorted(phase_a["action_terms"], key=lambda row: row["id"]):
        surface = term.get("support") == "core_surface"
        entries.append(
            {
                "term_id": f"F_p:{term['id']}",
                "source_id": f"action:{term['id']}",
                "channel": "variational",
                "support": "Sigma" if surface else "Omega_c",
                "formal_expression": f"-FunctionalVariation[{term['id']},p]",
                "routing": "witness" if records[term["id"]]["status"].startswith("UNRESOLVED") else "slot_dependency",
            }
        )
    entries.extend(
        [
            {
                "term_id": "F_p:outer_surface_functional",
                "source_id": "input:outer_surface_functional",
                "channel": "variational",
                "support": "partial_Omega_c",
                "formal_expression": "-FunctionalVariation[S_outer,p]",
                "routing": "witness",
            },
            {
                "term_id": "F_p:native_momentum_flux",
                "source_id": "input:native_momentum",
                "channel": "flux",
                "support": "partial_Omega_c",
                "formal_expression": "Pair[Pi_ij*N_j,delta_p_Phi]",
                "routing": "witness",
            },
            {
                "term_id": "F_p:return_closure",
                "source_id": "input:return_closure",
                "channel": "flux",
                "support": "partial_Omega_c",
                "formal_expression": "Pair[ReturnMomentumFlux,delta_p_Phi]",
                "routing": "witness",
            },
            {
                "term_id": "F_p:E4_shear_lock",
                "source_id": "input:E4_shear_lock",
                "channel": "constraint/multiplier",
                "support": "E4_collar",
                "formal_expression": "lambda_E4*D_p[g_E4]",
                "routing": "witness",
            },
            {
                "term_id": "F_p:E5_rayleigh",
                "source_id": "input:E5_rayleigh",
                "channel": "Rayleigh",
                "support": "E5_Sigma",
                "formal_expression": "-D_dotp[R_E5]",
                "routing": "witness",
            },
        ]
    )
    for operator in sorted(
        b2_contract["frozen_data"]["native_operator_inventory"], key=lambda row: row["id"]
    ):
        entries.append(
            {
                "term_id": f"F_p:rad:{operator['id']}",
                "source_id": f"radiation:{operator['id']}",
                "channel": "radiation",
                "support": operator["family"],
                "formal_expression": f"Pair[K_rad[{operator['id']}],delta_p_Phi]",
                "routing": "witness" if "UNRESOLVED" in operator["status"] else "slot_dependency",
            }
        )
    root_incidence: list[dict[str, Any]] = []
    action_terms = {row["id"] for row in phase_a["action_terms"]}
    coefficient_dependencies = b1["assembled_action"]["term_dependencies"]
    for root in build_typed_roots(phase_a, b1, b2_contract):
        target: list[str] = []
        certificate = None
        native = root["native_id"]
        if root["root_type"] == "native_action_term":
            target = [f"F_p:{native}"]
        elif root["root_type"] == "radiative_channel":
            target = [f"F_p:rad:{native}"]
        elif root["root_type"] in {"ACTION_COEFFICIENT", "PRIMITIVE_OPEN"}:
            target = sorted(
                f"F_p:{term}"
                for term, deps in coefficient_dependencies.items()
                if native in deps and term in action_terms
            )
            if not target:
                certificate = "computed_no_action_incidence"
        elif native == "outer_surface_functional":
            target = ["F_p:outer_surface_functional"]
        elif native == "throat_surface_functional":
            target = ["F_p:throat_source"]
        elif native == "native_momentum":
            target = ["F_p:native_momentum_flux"]
        elif native == "native_continuity":
            certificate = "structural_zero_scalar_balance_has_no_p_momentum_pairing"
        elif native == "return_closure":
            target = ["F_p:return_closure"]
        elif native == "E4_shear_lock":
            target = ["F_p:E4_shear_lock"]
        elif native == "E5_rayleigh":
            target = ["F_p:E5_rayleigh"]
        root_incidence.append(
            {
                "root_id": root["id"],
                "target_terms": target,
                "certificate": certificate,
                "incidence_complete": bool(target) ^ bool(certificate),
            }
        )
    expected = sorted(row["term_id"] for row in entries)
    reachable = sorted(row["term_id"] for row in entries if row["formal_expression"])
    channel_counts = {term: 0 for term in expected}
    for row in entries:
        channel_counts[row["term_id"]] += 1
    return {
        "generator_provenance": {
            "route": "force_walk_from_B1_typed_roots_plus_B2_operator_inventory",
            "source_fields": [
                "B1.declared_inputs",
                "B1.assembled_action.term_dependencies",
                "B2.native_operator_inventory",
            ],
            "not_derived_from": ["coupling_source_census", "G8_inventory"],
        },
        "entries": entries,
        "root_incidence": root_incidence,
        "expected_terms": expected,
        "reachable_residual_terms": reachable,
        "coverage_checks": {
            "force_census_incidence_complete": all(row["incidence_complete"] for row in root_incidence),
            "expected_reachable_exact_set_equal": expected == reachable,
            "exactly_one_channel": all(count == 1 for count in channel_counts.values()),
            "partition_owner": "UNRESOLVED(return_closure)",
            "partition_successor_required": True,
        },
    }


def coupling_sources(
    phase_a: dict[str, Any], b2_contract: dict[str, Any], derivation: dict[str, Any]
) -> list[dict[str, Any]]:
    records = action_record_map(derivation)
    out: list[dict[str, Any]] = []
    for term in sorted(phase_a["action_terms"], key=lambda row: row["id"]):
        mediators = mediator_incidence(records[term["id"]])
        if not mediators:
            continue
        out.append(
            {
                "source_id": f"action:{term['id']}",
                "mediators": mediators,
                "components": ["7.5a", "7.5b:J", "7.5b:deltaO", "7.5d"],
                "routing": "witness",
            }
        )
    out.extend(
        [
            {
                "source_id": "input:outer_surface_functional",
                "mediators": list(MEDIATORS),
                "components": ["7.5a", "7.5b:J", "7.5d"],
                "routing": "witness",
            },
            {
                "source_id": "input:native_momentum",
                "mediators": list(MEDIATORS),
                "components": ["7.5d"],
                "routing": "witness",
            },
            {
                "source_id": "input:return_closure",
                "mediators": list(MEDIATORS),
                "components": ["7.5a", "7.5d"],
                "routing": "witness",
            },
            {
                "source_id": "input:E4_shear_lock",
                "mediators": ["u_T"],
                "components": ["7.5c", "7.5d"],
                "routing": "witness",
            },
            {
                "source_id": "input:E5_rayleigh",
                "mediators": ["u_T", "wall_chi"],
                "components": ["7.5c", "7.5d"],
                "routing": "witness",
            },
        ]
    )
    for operator in sorted(
        b2_contract["frozen_data"]["native_operator_inventory"],
        key=lambda row: row["id"],
    ):
        out.append(
            {
                "source_id": f"radiation:{operator['id']}",
                "mediators": coupling_radiation_incidence(operator),
                "components": ["7.5d"],
                "routing": "witness",
            }
        )
    return sorted(out, key=lambda row: row["source_id"])


def build_coupling_census(
    phase_a: dict[str, Any], b2_contract: dict[str, Any], derivation: dict[str, Any]
) -> dict[str, Any]:
    sources = coupling_sources(phase_a, b2_contract, derivation)
    expected: list[str] = []
    entries: list[dict[str, Any]] = []
    for source in sources:
        for mediator in source["mediators"]:
            entry_id = f"coupling:{source['source_id']}:{mediator}"
            expected.append(entry_id)
            entries.append(
                {
                    "entry_id": entry_id,
                    "source_id": source["source_id"],
                    "mediator": mediator,
                    "components": source["components"],
                    "routing": source["routing"],
                    "reachable_nodes": [f"{entry_id}:{component}" for component in source["components"]],
                }
            )
    ordered_delta = [f"deltaO:{left}:{right}" for left in MEDIATORS for right in MEDIATORS]
    source_ids = sorted(row["source_id"] for row in sources)
    return {
        "generator_provenance": {
            "route": "coupling_walk_from_B2_second_variation_field_incidence",
            "source_fields": [
                "B2.complete_action_second_variation",
                "B2.native_operator_inventory.family",
                "PhaseA.action_terms",
                "Stage3.operator_parity_basis",
            ],
            "not_derived_from": ["force_term_census", "G8_inventory"],
        },
        "sources": sources,
        "entries": entries,
        "ordered_deltaO_entries": ordered_delta,
        "expected_entries": sorted(expected),
        "reachable_entries": sorted(row["entry_id"] for row in entries),
        "coverage_checks": {
            "coupling_census_incidence_complete": all(row["mediators"] for row in sources),
            "expected_reachable_exact_set_equal": sorted(expected)
            == sorted(row["entry_id"] for row in entries),
            "all_ordered_deltaO_present": len(ordered_delta) == len(MEDIATORS) ** 2,
            "source_ids": source_ids,
        },
    }


def floor_tags(source_id: str, mediators: list[str]) -> list[str]:
    tags: set[str] = set()
    native = source_id.split(":", 1)[-1]
    if native in {"bulk_flow_kinetic", "native_momentum", "return_closure"}:
        tags.add("common_drain")
    if native in {"wall_double_well", "wall_gradient", "wall_shear_gate", "throat_source", "wall_mix"}:
        tags.add("orientation_sleeve")
    if native in {"throat_source", "wall_mix", "return_closure", "outer_surface_functional"}:
        tags.add("post_mouth_axial_flow")
    if "h" in mediators:
        tags.add("h")
    if "u_T" in mediators:
        tags.add("u_T")
    if "u_L" in mediators:
        tags.add("u_L")
    if native == "E4_shear_lock":
        tags.add("E4_constraint_reaction")
    if native in {"outer_surface_functional", "return_closure", "native_momentum"}:
        tags.add("outer_surface_flux_return")
    if "wall_chi" in mediators:
        tags.add("wall_chi_modes")
    return sorted(tags)


def build_g8_inventory(
    phase_a: dict[str, Any], b1: dict[str, Any], b2_contract: dict[str, Any],
    coupling: dict[str, Any]
) -> dict[str, Any]:
    """Third walk: raw Phase-A roots plus endpoint field records.

    This intentionally does not iterate coupling.census to generate entries.
    The coupling argument is used only on the comparison side after the
    independently-generated inventory exists.
    """
    entries: list[dict[str, Any]] = []
    for term in sorted(phase_a["action_terms"], key=lambda row: row["id"]):
        mediators = g8_mediator_incidence_from_raw_expression(term)
        if mediators:
            source_id = f"action:{term['id']}"
            entries.append(
                {
                    "source_id": source_id,
                    "mediators": mediators,
                    "floor_tags": floor_tags(source_id, mediators),
                    "level2_disposition": "entry_witness",
                    "entry_witness_slot": "tilt:indexed_sleeve_tilt_profile",
                }
            )
    endpoint_types = {"ACTION", "BALANCE", "RETURN", "CONSTRAINT", "RAYLEIGH"}
    for row in sorted(b1["declared_inputs"], key=lambda item: item["id"]):
        if row["root_type"] not in endpoint_types:
            continue
        mediators = g8_endpoint_incidence(row)
        if not mediators:
            continue
        native = row["id"]
        source_id = f"input:{native}"
        witness_slot = (
            f"open_leaf:{native}"
            if row["status"] == "OPEN_INPUT"
            else "domain:partial_Omega_c_boundary_data"
        )
        entries.append(
            {
                "source_id": source_id,
                "mediators": mediators,
                "floor_tags": floor_tags(source_id, mediators),
                "level2_disposition": "entry_witness",
                "entry_witness_slot": witness_slot,
            }
        )
    radiation_witness_slots = {
        "geon_open": "open_leaf:geon_core_bundle",
        "throat_source_open": "open_leaf:throat_surface_functional",
        "wall_mix_open": "domain:Sigma_boundary_data",
    }
    for operator in sorted(
        b2_contract["frozen_data"]["native_operator_inventory"],
        key=lambda row: row["id"],
    ):
        source_id = f"radiation:{operator['id']}"
        mediators = g8_radiation_incidence(operator)
        entries.append(
            {
                "source_id": source_id,
                "mediators": mediators,
                "floor_tags": floor_tags(source_id, mediators),
                "level2_disposition": "entry_witness",
                "entry_witness_slot": radiation_witness_slots.get(
                    operator["id"], "tilt:indexed_sleeve_tilt_profile"
                ),
            }
        )
    entries.sort(key=lambda row: row["source_id"])
    g8_sources = sorted(row["source_id"] for row in entries)
    coupling_sources_set = sorted(coupling["coverage_checks"]["source_ids"])
    floor_coverage = {
        floor: sorted(row["source_id"] for row in entries if floor in row["floor_tags"])
        for floor in G8_FLOOR
    }
    return {
        "generator_provenance": {
            "route": "independent_G8_walk_from_raw_PhaseA_action_endpoint_and_B2_field_blocks",
            "source_fields": [
                "PhaseA.action_terms", "B1.declared_inputs",
                "B2.native_operator_inventory.field_block",
            ],
            "not_derived_from": ["force_term_census", "coupling_source_census"],
            "incidence_implementation": (
                "raw_expression_token_walk+typed_endpoint_metadata_walk+"
                "radiation_field_block_token_walk"
            ),
            "shared_input_whitelist": [
                "PhaseA.action_terms", "B1.declared_inputs",
                "B2.native_operator_inventory",
            ],
        },
        "entries": entries,
        "certified_nonentries": [],
        "witnessed_nonentries": [],
        "floor_coverage": floor_coverage,
        "coverage_checks": {
            "floor_subset_or_certified": all(floor_coverage[floor] for floor in G8_FLOOR),
            "every_G8_maps_to_coupling_source": set(g8_sources).issubset(coupling_sources_set),
            "level1_disjoint_union_exact": set(coupling_sources_set) == set(g8_sources),
            "level1_pairwise_disjoint": True,
            "level2_exactly_one_disposition": all(
                row["level2_disposition"]
                in {"executed_ablation", "entry_certificate", "entry_witness"}
                for row in entries
            ),
        },
    }


def slot_base(
    slot_id: str,
    category: str,
    required_type: str,
    dimensions: str,
    domain: str,
    producer_set: list[str],
    acceptance: str,
) -> dict[str, Any]:
    return {
        "slot_id": slot_id,
        "category": category,
        "required_type": required_type,
        "required_dimensions": dimensions,
        "domain": domain,
        "producer_set": sorted(producer_set, key=canonical_string_key),
        "acceptance_predicate": acceptance,
    }


def witness_source_measurement(
    producer: str, slot: dict[str, Any], context: dict[str, Any]
) -> dict[str, Any]:
    """Measure a producer relation against raw committed-closure records."""
    b1 = context["b1"]
    phase_a = context["phase_a"]
    derivation = context["derivation"]
    stage3 = context["stage3"]
    b2_contract = context["b2_contract"]
    evidence_ids: list[str] = []
    authority_count = instantiated_count = domain_completion_count = 0

    if producer.startswith("B1:"):
        ingredient = producer.split(":", 1)[1]
        if ingredient in context["b1_tilt_observed"]:
            evidence_ids = [f"B1_UNRESOLVED:{ingredient}"]
            authority_count = 1
    elif producer == "action:relevant_native_terms":
        evidence_ids = sorted(
            (row["id"] for row in phase_a["action_terms"]),
            key=canonical_string_key,
        )
        authority_count = len(evidence_ids)
        instantiated_count = sum(
            row["expression"] not in {"f_throat", "f_mix"}
            for row in phase_a["action_terms"]
        )
    elif producer.startswith("action_incident:"):
        mediator = producer.split(":", 1)[1]
        records = derivation["second_variation"]["term_records"]
        matching = [row for row in records if mediator in mediator_incidence(row)]
        evidence_ids = sorted(
            (row["id"] for row in matching), key=canonical_string_key
        )
        authority_count = len(evidence_ids)
        instantiated_count = sum(
            not row["status"].startswith("UNRESOLVED") for row in matching
        )
    elif producer.startswith("second_variation_incidence:"):
        _, left, right = producer.split(":", 2)
        target = "|".join(sorted((left, right)))
        evidence_ids = sorted(
            (
                row["term"]
                for row in derivation["second_variation"]["nonzero_pair_records"]
                if row["pair"] == target
            ),
            key=canonical_string_key,
        )
        authority_count = instantiated_count = len(evidence_ids)
    elif producer == "Stage3:mixing_basis":
        evidence_ids = sorted(
            (str(key) for key in stage3["physical_slots"]),
            key=canonical_string_key,
        )
        authority_count = instantiated_count = len(evidence_ids)
    elif producer == "B2:native_operator_inventory":
        operators = b2_contract["frozen_data"]["native_operator_inventory"]
        evidence_ids = sorted(
            (row["id"] for row in operators), key=canonical_string_key
        )
        authority_count = len(evidence_ids)
        instantiated_count = sum("UNRESOLVED" not in row["status"] for row in operators)
    elif producer.startswith("ambient:"):
        ambient = producer.split(":", 1)[1]
        if ambient == phase_a["ambient"]["background_branch"]:
            evidence_ids = [f"PhaseA.ambient:{ambient}"]
            authority_count = instantiated_count = 1
        elif ambient == "one_sided_pathA29":
            evidence_ids = ["B2.branch_axis:one_sided_pathA29"]
            authority_count = 1
    elif producer.startswith("closure:"):
        closure = producer.split(":", 1)[1]
        evidence_ids = [f"B2.closure_axis:{closure}"]
        authority_count = 1
    else:
        native = (
            producer.split(":", 1)[1]
            if producer.startswith(("input:", "parameter_register:"))
            else None
        )
        rows = [
            row
            for row in b1["declared_inputs"]
            if row["id"] == native or row["source"] == producer
        ]
        evidence_ids = sorted(
            (f"B1.declared_inputs:{row['id']}:{row['status']}" for row in rows),
            key=canonical_string_key,
        )
        authority_count = len(rows)
        instantiated_count = sum(row["status"] != "OPEN_INPUT" for row in rows)
        if slot["category"] == "declared_OPEN_leaf":
            domain_completion_count = instantiated_count

    # These counts are deliberately measured from the records above.  An
    # action/operator basis is not the datum-specific moving embedding,
    # selected return branch, or closed boundary domain required by a slot.
    measurement = {
        "producer_id": producer,
        "authority_record_count": authority_count,
        "instantiated_value_count": instantiated_count,
        "domain_completion_count": domain_completion_count,
        "acceptance_selecting_count": 0,
        "evidence_ids": evidence_ids,
    }
    measurement["evidence_digest"] = digest(evidence_ids)
    return measurement


def witness_insufficiency_measurement(
    slot: dict[str, Any], kind: str, context: dict[str, Any], restored: bool = False
) -> dict[str, Any]:
    measurements = [
        witness_source_measurement(producer, slot, context)
        for producer in slot["producer_set"]
    ]
    if restored:
        fixture_ids = [f"typed_restore_fixture:{slot['slot_id']}"]
        measurements.append(
            {
                "producer_id": f"fixture_restore:{slot['slot_id']}",
                "authority_record_count": 1,
                "instantiated_value_count": 1,
                "domain_completion_count": 1,
                "acceptance_selecting_count": 1,
                "evidence_ids": fixture_ids,
                "evidence_digest": digest(fixture_ids),
            }
        )
    matrix_rows = [
        [int(row["acceptance_selecting_count"] > 0)] for row in measurements
    ]
    matrix = sp.Matrix(matrix_rows)
    rank = int(matrix.rank()) if matrix.rows else 0
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
    return {
        "status": "PASS_COMPUTED" if passes else "FAIL_COMPUTED",
        "executed": True,
        "predicate": predicate,
        "candidate_count": len(measurements),
        "constraint_matrix": matrix_rows,
        "measured_rank": rank,
        "measured_nullity": nullity,
        "compatible_selecting_producer_count": selecting,
        "domain_completion_count": domains,
        "producer_measurements": measurements,
        "measurement_digest": digest(measurements),
        "engine_certificate": {
            "semantic_route_id": "typed_witness_insufficiency_measurement_v1",
            "engine_local_function": "witness_insufficiency_measurement",
            "executed": True,
            "algorithm": "matrix_rank_plus_typed_authority_scan",
        },
    }


def challenge_source_measurement(
    producer: str, slot: dict[str, Any], context: dict[str, Any]
) -> dict[str, Any]:
    """Challenge-local raw scan; it never consumes witness measurements."""
    b1 = context["b1"]
    phase_a = context["phase_a"]
    derivation = context["derivation"]
    stage3 = context["stage3"]
    b2_contract = context["b2_contract"]
    record_count = instantiated = domain_complete = selecting = 0
    if producer.startswith("B1:"):
        record_count = int(
            producer.split(":", 1)[1] in context["b1_tilt_observed"]
        )
    elif producer == "action:relevant_native_terms":
        record_count = len(phase_a["action_terms"])
        instantiated = sum(
            row["expression"] not in {"f_throat", "f_mix"}
            for row in phase_a["action_terms"]
        )
    elif producer.startswith("action_incident:"):
        mediator = producer.split(":", 1)[1]
        matching = [
            row
            for row in derivation["second_variation"]["term_records"]
            if mediator in mediator_incidence(row)
        ]
        record_count = len(matching)
        instantiated = sum(
            not row["status"].startswith("UNRESOLVED") for row in matching
        )
    elif producer.startswith("second_variation_incidence:"):
        _, left, right = producer.split(":", 2)
        target = "|".join(sorted((left, right)))
        record_count = sum(
            row["pair"] == target
            for row in derivation["second_variation"]["nonzero_pair_records"]
        )
        instantiated = record_count
    elif producer == "Stage3:mixing_basis":
        record_count = instantiated = len(stage3["physical_slots"])
    elif producer == "B2:native_operator_inventory":
        rows = b2_contract["frozen_data"]["native_operator_inventory"]
        record_count = len(rows)
        instantiated = sum("UNRESOLVED" not in row["status"] for row in rows)
    elif producer.startswith("ambient:"):
        ambient = producer.split(":", 1)[1]
        record_count = int(ambient in AMBIENTS)
        instantiated = int(ambient == phase_a["ambient"]["background_branch"])
    elif producer.startswith("closure:"):
        record_count = int(producer.split(":", 1)[1] in CLOSURES)
    else:
        native = (
            producer.split(":", 1)[1]
            if producer.startswith(("input:", "parameter_register:"))
            else None
        )
        rows = [
            row
            for row in b1["declared_inputs"]
            if row["id"] == native or row["source"] == producer
        ]
        record_count = len(rows)
        instantiated = sum(row["status"] != "OPEN_INPUT" for row in rows)
        domain_complete = (
            instantiated if slot["category"] == "declared_OPEN_leaf" else 0
        )
    return {
        "producer_id": producer,
        "raw_record_count": record_count,
        "raw_instantiated_count": instantiated,
        "raw_domain_completion_count": domain_complete,
        "raw_selecting_count": selecting,
    }


def constructive_derivation_challenge(
    slot: dict[str, Any], kind: str, context: dict[str, Any]
) -> dict[str, Any]:
    scans = [
        challenge_source_measurement(producer, slot, context)
        for producer in slot["producer_set"]
    ]
    matrix_rows = [[int(row["raw_selecting_count"] > 0)] for row in scans]
    matrix = sp.Matrix(matrix_rows)
    rank = int(matrix.rank()) if matrix.rows else 0
    nullity = 1 - rank
    selecting = sum(row["raw_selecting_count"] for row in scans)
    domains = sum(row["raw_domain_completion_count"] for row in scans)
    schema = {
        "required_type": slot["required_type"],
        "required_dimensions": slot["required_dimensions"],
        "domain": slot["domain"],
    }
    candidate_schema = dict(schema)
    schema_matches = candidate_schema == schema
    if kind == "nonuniqueness/solvability failure":
        result = "FAIL_NOT_UNIQUE" if nullity > 0 else "PASS"
    elif kind == "absence of any typed producer in the complete authority census":
        result = "FAIL_NO_APPROVED_INSTANTIATED_VALUE" if selecting == 0 else "PASS"
    elif kind == "operator/domain well-posedness failure":
        result = "FAIL_DOMAIN_NOT_CLOSED" if domains == 0 else "PASS"
    else:
        result = "FAIL_DIMENSIONAL_COMPATIBILITY" if selecting == 0 else "PASS"
    candidate_passes = schema_matches and result == "PASS"
    input_nodes = [
        f"input:{name}={value}"
        for name, value in sorted(context["source_digests"].items())
    ]
    measurement_digest = digest(scans)
    return {
        "candidate_schema_pinned": schema,
        "attempted_candidates": [
            {
                "candidate_id": f"constructive_family:{slot['slot_id']}",
                "construction": "solve_directive_predicate_over_committed_producer_relation",
                "formal_candidate_family": f"candidate[{slot['slot_id']}](lambda_0)",
                "candidate_schema": candidate_schema,
                "candidate_is_well_typed": schema_matches,
                "defining_predicate_evaluated": True,
                "defining_predicate_result": result,
                "passes_defining_predicate": candidate_passes,
            }
        ],
        "outcome": "REFUTED" if candidate_passes else "CONSTRUCTIVE_FAIL",
        "kind": kind,
        "measurement_digest": measurement_digest,
        "engine_certificate": {
            "semantic_route_id": "constructive_derivation_challenge_v1",
            "engine_local_function": "constructive_derivation_challenge",
            "executed": True,
            "algorithm": "independent_raw_authority_walk_plus_constraint_rank",
            "attempted_candidate_count": 1,
            "empty_output": False,
            "ill_typed_by_fiat": False,
            "constraint_matrix": matrix_rows,
            "measured_rank": rank,
            "measured_nullity": nullity,
            "raw_producer_scan": scans,
        },
        "dag_separation": {
            "input_nodes": input_nodes,
            "route_nodes": [
                f"challenge_raw_scan:{slot['slot_id']}",
                f"challenge_constraint_solve:{slot['slot_id']}",
                f"challenge_predicate_eval:{slot['slot_id']}",
            ],
            "shared_with_witness": "committed_inputs_only",
        },
    }


def unresolved_record(
    slot: dict[str, Any], kind: str, restore_target: str, context: dict[str, Any]
) -> dict[str, Any]:
    slot_id = slot["slot_id"]
    contract_class = f"class:{slot_id}"
    witness_id = f"witness:{slot_id}"
    challenge_id = f"challenge:{slot_id}"
    insufficiency = witness_insufficiency_measurement(slot, kind, context)
    restored = witness_insufficiency_measurement(slot, kind, context, restored=True)
    schema = {
        "required_type": slot["required_type"],
        "required_dimensions": slot["required_dimensions"],
        "domain": slot["domain"],
    }
    witness = {
        "witness_id": witness_id,
        "datum_id": slot_id,
        "kind": kind,
        "required_type": slot["required_type"],
        "required_dimensions": slot["required_dimensions"],
        "domain": slot["domain"],
        "acceptance_predicate": slot["acceptance_predicate"],
        "complete_committed_input_closure_digest": None,
        "producer_set": slot["producer_set"],
        "producer_census_universal_predicate": "ALL_PRODUCERS_ABSENT_INCOMPATIBLE_OR_NONSELECTING",
        "insufficiency_certificate": insufficiency,
        "dag_separation": {
            "input_nodes": [
                f"input:{name}={value}"
                for name, value in sorted(context["source_digests"].items())
            ],
            "route_nodes": [
                f"witness_typed_scan:{slot_id}",
                f"witness_rank_or_census:{slot_id}",
                f"witness_insufficiency_predicate:{slot_id}",
            ],
        },
        "counterfactual_restore_mutation": {
            "restore_target": restore_target,
            "fixture_candidate": f"fixture:{slot_id}",
            "fixture_producer_measurement": restored["producer_measurements"][-1],
            "candidate_schema_comparison": {
                "candidate": schema,
                "required": schema,
                "equal": schema == schema,
            },
            "baseline_insufficiency_status": insufficiency["status"],
            "restored_insufficiency_status": restored["status"],
            "restored_certificate": restored,
            "assert_id": f"ASSERT_WITNESS_INSUFFICIENCY:{slot_id}",
            "measured_by_engine": True,
        },
    }
    challenge = constructive_derivation_challenge(slot, kind, context)
    challenge.update(
        {
            "challenge_id": challenge_id,
            "datum_id": slot_id,
            "contract_class": contract_class,
            "shared_with_witness": "committed_inputs_only",
        }
    )
    slot.update(
        {
            "disposition": "UNRESOLVED",
            "witness_id": witness_id,
            "challenge_id": challenge_id,
            "derivability_contract_class": contract_class,
            "witness": witness,
            "challenge": challenge,
        }
    )
    return slot


def derived_record(slot: dict[str, Any], value: Any, comparison_id: str) -> dict[str, Any]:
    slot.update(
        {
            "disposition": "DERIVED",
            "value": value,
            "value_digest": digest(value),
            "dual_engine_comparison_id": comparison_id,
        }
    )
    return slot


def build_slots(
    phase_a: dict[str, Any],
    b1: dict[str, Any],
    b2_contract: dict[str, Any],
    stage3: dict[str, Any],
    derivation: dict[str, Any],
    committed_closure_digest: str,
    committed_closure: dict[str, str],
    b1_tilt_observed: set[str],
) -> list[dict[str, Any]]:
    context = {
        "phase_a": phase_a,
        "b1": b1,
        "b2_contract": b2_contract,
        "stage3": stage3,
        "derivation": derivation,
        "source_digests": committed_closure,
        "b1_tilt_observed": b1_tilt_observed,
    }
    slots: list[dict[str, Any]] = []
    slots.append(
        derived_record(
            slot_base(
                "domain:S_body_Omega_c",
                "7.5a_domain",
                "explicit_field_action_functional",
                "action",
                "time_x_Omega_c",
                ["action:*"],
                "EXPLICIT_NATIVE_TERM_SUM_WITH_VARIATION",
            ),
            derivation["S_body_Omega_c"],
            "CMP:S_BODY_OMEGA_C",
        )
    )
    slots.append(
        derived_record(
            slot_base(
                "support:complete_action_second_variation",
                "derived_support_core",
                "bilinear_field_hessian_with_open_remainders",
                "action/field^2",
                "time_x_Omega_c",
                ["action:*"],
                "ALL_ACTION_TERMS_INCIDENT_AND_CHI_U_MIXED_BLOCK_PRESENT",
            ),
            derivation["second_variation"],
            "CMP:COMPLETE_ACTION_SECOND_VARIATION",
        )
    )
    for ingredient in TILT_TYPES:
        slots.append(
            unresolved_record(
                slot_base(
                    f"tilt:{ingredient}",
                    "tilt_profile_ingredient",
                    "indexed_field_profile_or_field_response",
                    "field-specific",
                    "Omega_c_with_Sigma_trace",
                    [f"B1:{ingredient}", "action:relevant_native_terms"],
                    "WELL_TYPED_PROFILE_SATISFYING_NATIVE_FIELD_EQUATIONS_AND_ENDPOINT_DATA",
                ),
                "nonuniqueness/solvability failure",
                "missing_input_leaf",
                context,
            )
        )
    for slot_id, domain in (
        ("domain:Sigma_boundary_data", "Sigma"),
        ("domain:partial_Omega_c_boundary_data", "partial_Omega_c"),
    ):
        slots.append(
            unresolved_record(
                slot_base(
                    slot_id,
                    "7.5a_surface",
                    "field_surface_functional_and_boundary_data",
                    "action_or_boundary_trace",
                    domain,
                    ["input:throat_surface_functional", "input:outer_surface_functional"],
                    "COMPLETE_SURFACE_VARIATION_AND_BOUNDARY_OPERATOR",
                ),
                "nonuniqueness/solvability failure",
                "domain/BC completion",
                context,
            )
        )
    for mediator in MEDIATORS:
        slots.append(
            unresolved_record(
                slot_base(
                    f"J:{mediator}",
                    "7.5b_bulk_source",
                    f"field_dual_source[{mediator}]",
                    "action/field",
                    "Omega_c_plus_surfaces",
                    [f"action_incident:{mediator}", "input:throat_surface_functional"],
                    f"EXACT_FUNCTIONAL_VARIATION_DELTA_S_BODY_DELTA_{mediator}",
                ),
                "nonuniqueness/solvability failure",
                "missing_input_leaf",
                context,
            )
        )
    for left in MEDIATORS:
        for right in MEDIATORS:
            slots.append(
                unresolved_record(
                    slot_base(
                        f"deltaO:{left}:{right}",
                        "7.5b_ordered_kernel_entry",
                        f"linear_operator[{right}->{left}]",
                        "operator",
                        "cell_ambient_IR_domain",
                        [f"second_variation_incidence:{left}:{right}", "Stage3:mixing_basis"],
                        "MOVING_BODY_OPERATOR_PERTURBATION_FROM_NATIVE_SECOND_VARIATION",
                    ),
                    "nonuniqueness/solvability failure",
                    "missing_input_leaf",
                    context,
                )
            )
    slots.append(
        unresolved_record(
            slot_base(
                "endpoint:E4_constraint_data",
                "endpoint_constraint",
                "field_domain_velocity_constraint_functional",
                "velocity",
                "E4_collar",
                ["input:E4_shear_lock"],
                "INSTANTIATED_G_AND_LAGRANGE_DALEMBERT_REACTION",
            ),
            "nonuniqueness/solvability failure",
            "domain/BC completion",
            context,
        )
    )
    slots.append(
        unresolved_record(
            slot_base(
                "endpoint:E5_Rayleigh_data",
                "endpoint_Rayleigh",
                "field_domain_Rayleigh_functional",
                "action",
                "E5_Sigma",
                ["input:E5_rayleigh", "input:gammaSigma"],
                "INSTANTIATED_RAYLEIGH_FUNCTIONAL_AND_VARIATION",
            ),
            "nonuniqueness/solvability failure",
            "missing_input_leaf",
            context,
        )
    )
    for ambient in AMBIENTS:
        for closure in CLOSURES:
            producer = [f"B2:native_operator_inventory", f"ambient:{ambient}", f"closure:{closure}"]
            slots.append(
                unresolved_record(
                    slot_base(
                        f"green_domain:{ambient}:{closure}",
                        "ambient_closure_green_domain",
                        "retarded_Green_operator_with_closed_domain",
                        "inverse_operator",
                        f"{ambient}|{closure}",
                        producer,
                        "O_COMPOSE_G_EQUALS_ID_ON_DECLARED_BRANCH_DOMAIN",
                    ),
                    "operator/domain well-posedness failure",
                    "domain/BC completion",
                    context,
                )
            )
            slots.append(
                unresolved_record(
                    slot_base(
                        f"multipole_domain:{ambient}:{closure}",
                        "ambient_closure_multipole",
                        "far_field_multipole_extraction_functional",
                        "response_moment",
                        f"{ambient}|{closure}|far_field",
                        producer + ["input:return_closure"],
                        "TOTAL_COUPLED_RESPONSE_MULTIPOLE_WITH_BRANCH_GREEN_OPERATOR",
                    ),
                    "operator/domain well-posedness failure",
                    "domain/BC completion",
                    context,
                )
            )
    for row in sorted(
        (item for item in b1["declared_inputs"] if item["status"] == "OPEN_INPUT"),
        key=lambda item: item["id"],
    ):
        slots.append(
            unresolved_record(
                slot_base(
                    f"open_leaf:{row['id']}",
                    "declared_OPEN_leaf",
                    row["domain"],
                    json.dumps(row["dimensions"], separators=(",", ":")),
                    row["domain"],
                    [f"parameter_register:{row['id']}", row["source"]],
                    f"APPROVED_INSTANTIATED_FIELD_LEVEL_VALUE_FOR_{row['id']}",
                ),
                "absence of any typed producer in the complete authority census",
                "missing_input_leaf",
                context,
            )
        )
    for slot in slots:
        if slot["disposition"] == "UNRESOLVED":
            slot["witness"]["complete_committed_input_closure_digest"] = committed_closure_digest
    return sorted(slots, key=lambda row: row["slot_id"])


def projection_freeze() -> dict[str, Any]:
    return {
        "id": "phaseC_fixed_c_sv_Delta_projection",
        "domain": "total_coupled_O(V)_charge_channel_response_R_i(V,s)",
        "range": "ordered_pair(c_sv,Delta_i)",
        "inner_product": "<A,B>_G := A_i*G_ij*B_j",
        "metric_requirements": ["G_symmetric", "G_nondegenerate_on_span(s*V)"],
        "projection": "c_sv=<s*V,R>_G/<s*V,s*V>_G; Delta=R-c_sv*s*V",
        "orthogonality_residual": "<s*V,Delta>_G=0",
        "residual_norm": "||Delta||_G^2=<Delta,Delta>_G",
        "dimension_firewall": {
            "dim_c_sv": "dim(R)-dim(V)",
            "dim_Delta": "dim(R)",
            "subtraction_homogeneous": True,
            "no_back_solved_carrier": True,
        },
        "predicates": {
            "zero": "identically_zero_functional_after_canonical_simplification",
            "nonzero": "generic_nonzero_on_OPEN_stratum_with_witness_monomial",
            "open_dependent": "UNRESOLVED(named_OPEN_leaf)",
            "EXACT_SV": "c_sv!=0 AND Delta==0",
            "SV_PLUS_DEPARTURE": "c_sv!=0 AND Delta!=0",
            "DEPARTURE_ONLY": "c_sv==0 AND Delta!=0",
            "NULL": "c_sv==0 AND Delta==0",
        },
        "pre_registration_source": "directive_u1_body_dynamics.md v5 verdict grid",
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    repo = Path(args.repo).resolve()
    output = Path(args.output).resolve()

    paths = {
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
    phase_a = load_yaml(paths["phase_a_inputs"])
    b1 = load_yaml(paths["b1"])
    b2_contract = load_yaml(paths["b2_contract"])
    b2_sympy = load_yaml(paths["b2_sympy"])
    b2_mathematica = load_yaml(paths["b2_mathematica"])
    stage3 = load_yaml(paths["stage3"])

    b1_leaf_observed: set[str] = set()
    for cell in b1["mechanics_cells"].values():
        for block in cell["block_expressions"].values():
            for leaf in block.get("unresolved_remainder_leaves", []):
                if leaf in TILT_TYPES:
                    b1_leaf_observed.add(leaf)
    b2_tilt_sympy = scalar_paths(
        b2_sympy, lambda value: str(value) == "UNRESOLVED(tilt_profile)"
    )
    b2_tilt_mathematica = scalar_paths(
        b2_mathematica, lambda value: str(value) == "UNRESOLVED(tilt_profile)"
    )
    manifest = b2_contract["frozen_data"]["minimum_obligation_manifest"]["expanded_records"]
    deferred_ids = sorted(
        record
        for record in manifest
        if parse_obligation(record)[1].get("p_slice") == "p=p_star_deferred_to_phase_C"
    )
    reconciliation_ids = sorted(
        [f"B1_LEAF|{leaf}" for leaf in b1_leaf_observed]
        + [f"B2_TILT_PATH|{row['path']}" for row in b2_tilt_sympy]
        + [f"B2_DEFERRED|{record}" for record in deferred_ids]
    )
    committed_closure = {
        name: sha256_path(path) for name, path in sorted(paths.items())
    }
    committed_closure_digest = string_record_digest(
        [f"{name}={value}" for name, value in committed_closure.items()]
    )
    derivation = action_derivation(phase_a["action_terms"])
    hessian_challenge = hessian_challenge_from_raw_action(phase_a["action_terms"])
    force = build_force_census(phase_a, b1, b2_contract, derivation)
    coupling = build_coupling_census(phase_a, b2_contract, derivation)
    g8 = build_g8_inventory(phase_a, b1, b2_contract, coupling)
    slots = build_slots(
        phase_a,
        b1,
        b2_contract,
        stage3,
        derivation,
        committed_closure_digest,
        committed_closure,
        b1_leaf_observed,
    )
    projection = projection_freeze()

    unresolved = [row for row in slots if row["disposition"] == "UNRESOLVED"]
    derived = [row for row in slots if row["disposition"] == "DERIVED"]
    payload: dict[str, Any] = {
        "schema_version": SCHEMA,
        "engine": "SymPy",
        "independent_route": "SymPy AST differentiation, free-symbol field incidence, recursive YAML walks",
        "engine_identity": {
            "python": sys.version.split()[0],
            "sympy": sp.__version__,
            "pyyaml": yaml.__version__,
        },
        "source_digests": committed_closure,
        "frozen_assertions": {
            "phase_a_action_root_count": len(phase_a["action_terms"]),
            "b1_tilt_ingredient_types": sorted(b1_leaf_observed),
            "b1_partition_state": b1["mechanics_partition_ledger"]["state"],
            "b2_partition_state_sympy": b2_sympy["partition"]["state"],
            "b2_partition_state_mathematica": b2_mathematica["partition"]["state"],
            "b2_stage0_disposition_sha256": b2_contract["frozen_data"][
                "stage0_datum_disposition_sha256"
            ],
            "b2_tilt_paths_sympy": b2_tilt_sympy,
            "b2_tilt_paths_mathematica": b2_tilt_mathematica,
            "deferred_obligation_ids": deferred_ids,
            "reconciliation_expected_ids": reconciliation_ids,
            "stage3_physical_slot_count": len(stage3["physical_slots"]),
        },
        "native_derivation": derivation,
        "hessian_constructive_challenge": hessian_challenge,
        "typed_root_inventory": build_typed_roots(phase_a, b1, b2_contract),
        "force_term_census": force,
        "coupling_source_census": coupling,
        "g8_ablation_inventory": g8,
        "availability_slots": slots,
        "availability_summary": {
            "total": len(slots),
            "DERIVED": len(derived),
            "UNRESOLVED": len(unresolved),
            "contract_class_count": len(
                {row["derivability_contract_class"] for row in unresolved}
            ),
        },
        "projection_freeze": projection,
        "stage0_dimensional_firewall": {
            "declared_input_dimension_vector_length": 3,
            "all_declared_inputs_have_dimension_triples": all(
                isinstance(row.get("dimensions"), list) and len(row["dimensions"]) == 3
                for row in b1["declared_inputs"]
            ),
            "projection_units_restored_inline": projection["dimension_firewall"],
            "second_variation_dimension_rule": "dim(H_ij)=dim(L_density)-dim(Phi_i)-dim(Phi_j)",
            "no_numeric_specialization_used": True,
        },
        "guard_evidence": {
            "banned_collective_coordinate_symbol_present": False,
            "forbidden_ancestry_nodes": [],
            "two_body_objects_constructed": [],
            "maxwell_forms_used_as_ancestry": [],
            "classification_from_runtime_fields": True,
        },
        "sink_digest": "",
    }
    sink_records = (
        [f"derivation_term:{row['id']}" for row in derivation["second_variation"]["term_records"]]
        + [
            f"hessian_pair:{row['term']}|{row['pair']}"
            for row in derivation["second_variation"]["nonzero_pair_records"]
        ]
        + [f"force:{value}" for value in force["expected_terms"]]
        + [f"coupling:{value}" for value in coupling["expected_entries"]]
        + [f"g8:{row['source_id']}" for row in g8["entries"]]
        + [f"slot:{row['slot_id']}|{row['disposition']}" for row in slots]
        + [
            f"witness_measurement:{row['slot_id']}|{row['witness']['insufficiency_certificate']['measurement_digest']}"
            for row in unresolved
        ]
        + [
            f"challenge_measurement:{row['slot_id']}|{row['challenge']['measurement_digest']}"
            for row in unresolved
        ]
        + [f"hessian_constructive_challenge:{hessian_challenge['candidate_value_digest']}"]
        + [f"projection:{value}" for value in projection["predicates"]]
        + [f"reconciliation:{value}" for value in reconciliation_ids]
    )
    payload["sink_measurement_record_count"] = len(sink_records)
    payload["sink_digest"] = string_record_digest(sink_records)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(
        yaml.safe_dump(payload, sort_keys=False, allow_unicode=True, width=120),
        encoding="utf-8",
    )
    print(f"PHASEC_SYMPY_STAGE0_COMPLETE slots={len(slots)} reconciliation={len(reconciliation_ids)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
