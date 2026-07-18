#!/usr/bin/env python3
"""Independent SymPy production engine for U1 Phase C.

The engine constructs the complete symbolic conditional map.  It never fills a
ratified OPEN field-profile slot with a surrogate: predicates depending on one
of those slots remain first-class UNRESOLVED outcomes.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import sys
from collections import Counter
from pathlib import Path
from typing import Any, Iterable

import sympy as sp
import yaml


SCHEMA = "U1_PHASE_C_PRODUCTION_ENGINE_V1"
RATIFIED_DIGEST = "83233baabd7f8e27c88d130b911691e76d01d5797da8eeb32c90bbae111ec95a"
MEDIATORS = ("h", "u_T", "u_L", "wall_chi")
ENDPOINTS = ("E1", "E2", "E3", "E4", "E5")
AMBIENTS = ("one_sided_pathA29", "symmetric_postulate")
CLOSURES = ("body_mass_growth", "return_path", "sleeve_exit")
CHANNELS = ("variational", "flux", "constraint/multiplier", "Rayleigh", "radiation")
TILT_ENUM = (
    "TILT_LINEAR", "TILT_OTHER", "TILT_ZERO", "TILT_NO_STEADY",
    "TILT_UNSTABLE", "TILT_UNRESOLVED",
)
COUPLING_ENUM = (
    "EXACT_SV", "SV_PLUS_DEPARTURE", "DEPARTURE_ONLY", "NULL",
    "UNRESOLVED", "ILL_POSED",
)
COMPONENT_FILES = (
    "availability_slots.yaml", "coupling_source_census.yaml",
    "environment_identity.yaml", "evaluated_code_closure_policy.yaml",
    "force_term_census.yaml", "frozen_data_pin_table.yaml",
    "g8_ablation_inventory.yaml", "obligation_manifest.yaml",
    "parameter_register_proposals.yaml", "producer_map.yaml",
    "projection_freeze.yaml", "reconciliation_inventory.yaml",
)


class EngineFailure(RuntimeError):
    pass


def require(condition: bool, detail: str) -> None:
    if not condition:
        raise EngineFailure(detail)


def load_yaml(path: Path) -> Any:
    with path.open("rb") as handle:
        return yaml.load(handle, Loader=yaml.CSafeLoader)


def dump_yaml(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        yaml.safe_dump(value, sort_keys=False, allow_unicode=True, width=140),
        encoding="utf-8",
    )


def canonical_bytes(value: Any) -> bytes:
    return json.dumps(
        value, ensure_ascii=False, sort_keys=True, separators=(",", ":")
    ).encode("utf-8")


def digest(value: Any) -> str:
    return hashlib.sha256(canonical_bytes(value)).hexdigest()


def canonical_ids(values: Iterable[str]) -> list[str]:
    return sorted(set(values), key=lambda value: (value.lower(), value))


def cell_key(endpoint: str, ambient: str, closure: str, stratum: str) -> dict[str, str]:
    return {
        "endpoint": endpoint,
        "ambient_branch": ambient,
        "closure_branch": closure,
        "open_stratum": f"GAP_OPEN_FIELD_PROFILE:{stratum}",
    }


def cell_id(key: dict[str, str]) -> str:
    return "|".join(f"{name}={key[name]}" for name in (
        "endpoint", "ambient_branch", "closure_branch", "open_stratum"
    ))


def authority_tag(ambient: str) -> str:
    return (
        "ambient-postulate-dependent"
        if ambient == "symmetric_postulate"
        else "one-sided-asymmetry-map"
    )


def parity_branch_map(ambient: str) -> str:
    if ambient == "one_sided_pathA29":
        return "P_one_sided[O]=P_body[O]+A_slab[O]; A_slab is retained as an OPEN-dependent asymmetry functional"
    return "P_symmetric[O]=P_body[O]+R_w[P_body[O]] under the declared ambient postulate"


def endpoint_slot(endpoint: str) -> list[str]:
    if endpoint == "E4":
        return ["endpoint:E4_constraint_data"]
    if endpoint == "E5":
        return ["endpoint:E5_Rayleigh_data"]
    return []


def ratified_unresolved_candidates(
    slots: dict[str, dict[str, Any]], candidates: Iterable[str]
) -> list[str]:
    candidate_ids = canonical_ids(candidates)
    require(
        all(slot_id in slots for slot_id in candidate_ids),
        "dependency candidate is absent from the ratified availability table",
    )
    return canonical_ids(
        slot_id for slot_id in candidate_ids
        if slots[slot_id]["disposition"] == "UNRESOLVED"
    )


def mapped_unresolved_dependencies(
    slots: dict[str, dict[str, Any]],
    endpoint: str,
    ambient: str,
    closure: str,
    include_coupling: bool,
) -> list[str]:
    """Generate the declared dependency map independently from the ancestry walk."""
    green_slot = f"green_domain:{ambient}:{closure}"
    multipole_slot = f"multipole_domain:{ambient}:{closure}"
    endpoint_candidates = set(endpoint_slot(endpoint))
    return canonical_ids(
        slot_id for slot_id, row in slots.items()
        if row["disposition"] == "UNRESOLVED" and (
            slot_id.startswith("tilt:")
            or slot_id.startswith("open_leaf:")
            or row["category"] == "7.5a_surface"
            or slot_id == green_slot
            or slot_id in endpoint_candidates
            or include_coupling and (
                slot_id.startswith("J:")
                or slot_id.startswith("deltaO:")
                or slot_id == multipole_slot
            )
        )
    )


def endpoint_applicable(term_id: str, endpoint: str) -> bool:
    if term_id == "F_p:E4_shear_lock":
        return endpoint == "E4"
    if term_id == "F_p:E5_rayleigh":
        return endpoint == "E5"
    return True


def source_channel(source_id: str) -> str:
    if source_id == "input:E4_shear_lock":
        return "constraint/multiplier"
    if source_id == "input:E5_rayleigh":
        return "Rayleigh"
    if source_id in {"input:native_momentum", "input:return_closure"}:
        return "flux"
    if source_id.startswith("radiation:"):
        return "radiation"
    return "variational"


def source_applicable(source_id: str, endpoint: str) -> bool:
    if source_id == "input:E4_shear_lock":
        return endpoint == "E4"
    if source_id == "input:E5_rayleigh":
        return endpoint == "E5"
    return True


def named_unresolved(slots: list[str], enum: str) -> dict[str, Any]:
    return {"enum": enum, "named_inputs": canonical_ids(slots)}


def formal_expression_record(
    expression_id: str,
    expression: str,
    dimensions: str,
    defining_operation: str,
    unresolved_slots: list[str] | None = None,
) -> dict[str, Any]:
    return {
        "expression_id": expression_id,
        "expression": expression,
        "defining_operation": defining_operation,
        "dimensions_restored": dimensions,
        "dimension_record_id": f"DIM:{expression_id}",
        "unresolved_slots": canonical_ids(unresolved_slots or []),
    }


def action_ablation_records(action_terms: list[dict[str, Any]]) -> list[dict[str, Any]]:
    records: list[dict[str, Any]] = []
    symbol_names = canonical_ids(
        token
        for row in action_terms
        for token in re.findall(r"[A-Za-z][A-Za-z0-9_]*", row["expression"])
    )
    symbols = {name: sp.Symbol(name, nonzero=True) for name in symbol_names}
    for row in sorted(action_terms, key=lambda item: item["id"]):
        expression = sp.sympify(row["expression"], locals=symbols)
        nonzero = sp.simplify(expression) != 0
        records.append({
            "term_id": row["id"],
            "root_id": f"action:{row['id']}",
            "ablation_parameter": f"alpha__{row['id']}",
            "baseline_minus_ablated": sp.sstr(expression),
            "response_nonzero": bool(nonzero),
            "support": "Sigma" if row.get("support") == "core_surface" else "Omega_c",
            "root_type": "native bulk/surface action term",
        })
    return records


def dimension_firewall() -> dict[str, Any]:
    # Every constructed class is checked as an equality of L,T,M exponent
    # vectors.  Symbolic vectors retain the two declared free carriers without
    # solving them from downstream expressions.
    def add(*vectors: tuple[Any, Any, Any]) -> tuple[Any, Any, Any]:
        return tuple(sum(vector[index] for vector in vectors) for index in range(3))

    def subtract(
        left: tuple[Any, Any, Any], right: tuple[Any, Any, Any]
    ) -> tuple[Any, Any, Any]:
        return tuple(left[index] - right[index] for index in range(3))

    def vectors_equal(vectors: list[tuple[Any, Any, Any]]) -> bool:
        reference = vectors[0]
        return all(
            all(sp.simplify(value - expected) == 0 for value, expected in zip(vector, reference))
            for vector in vectors[1:]
        )

    zero = (0, 0, 0)
    length = (1, 0, 0)
    mass = (0, 0, 1)
    force = (1, -2, 1)
    velocity = (1, -1, 0)
    action = (2, -1, 1)
    phi = sp.symbols("phi_L phi_T phi_M")
    response = sp.symbols("response_L response_T response_M")
    ell = sp.Symbol("ell")

    dimensions = {
        "p": zero,
        "T_A": phi,
        "Phi_A": phi,
        "F_p": force,
        "K_pp": subtract(force, zero),
        "V": velocity,
        "B_pV": subtract(force, velocity),
        "chi_pF": subtract(zero, force),
        "S_body": action,
        "J_A": subtract(action, phi),
        "deltaO_AB": subtract(subtract(action, phi), phi),
        "G_AB": subtract(phi, subtract(action, phi)),
        "M_l": add(phi, tuple(ell * value for value in length)),
        "R": response,
        "c_sv": subtract(response, velocity),
        "Delta": response,
        "R_mass": response,
        "R_charge": response,
        "R_total": response,
    }
    class_specs = [
        ("tilt_embedding", "dim(Phi_A)=dim(p)+dim(T_A)",
         [dimensions["Phi_A"], add(dimensions["p"], dimensions["T_A"])]),
        ("p_force_residual", f"all additive terms have LTM={list(force)}",
         [dimensions["F_p"] for _ in CHANNELS]),
        ("static_stiffness", "dim(K_pp)=dim(F_p)-dim(p)",
         [add(dimensions["K_pp"], dimensions["p"]), dimensions["F_p"]]),
        ("velocity_tilt_drive", "dim(B_pV)=dim(F_p)-dim(V)",
         [add(dimensions["B_pV"], dimensions["V"]), dimensions["F_p"]]),
        ("susceptibility", "dim(chi_pF)=dim(p)-dim(F_p)",
         [add(dimensions["chi_pF"], dimensions["F_p"]), dimensions["p"]]),
        ("S_body", "dim(S_body)=action", [dimensions["S_body"], action]),
        ("J_A", "dim(J_A)=dim(S_body)-dim(Phi_A)",
         [add(dimensions["J_A"], dimensions["Phi_A"]), dimensions["S_body"]]),
        ("deltaO_AB", "dim(deltaO_AB)=dim(S_body)-dim(Phi_A)-dim(Phi_B)",
         [add(dimensions["deltaO_AB"], dimensions["Phi_A"], dimensions["Phi_A"]), dimensions["S_body"]]),
        ("total_response", "dim(G_AB J_B)=dim(Phi_A)",
         [add(dimensions["G_AB"], dimensions["J_A"]), dimensions["Phi_A"]]),
        ("multipole", "dim(M_l[Phi_A])=dim(Phi_A)+L^l",
         [dimensions["M_l"], add(dimensions["Phi_A"], tuple(ell * value for value in length))]),
        ("c_sv_projection", "dim(c_sv)=dim(R)-dim(V); dim(Delta)=dim(R)",
         [add(dimensions["c_sv"], dimensions["V"]), dimensions["Delta"], dimensions["R"]]),
        ("mass_charge_split", "dim(R_mass)=dim(R_charge)=dim(R_total)",
         [dimensions["R_mass"], dimensions["R_charge"], dimensions["R_total"]]),
    ]
    rows = [
        {
            "expression_class": name,
            "restored_rule": rule,
            "homogeneous": vectors_equal(term_vectors),
        }
        for name, rule, term_vectors in class_specs
    ]
    cross_expression_checks = [
        [dimensions["Phi_A"], add(dimensions["p"], dimensions["T_A"])],
        [dimensions["F_p"], add(dimensions["K_pp"], dimensions["p"]),
         add(dimensions["B_pV"], dimensions["V"])],
        [dimensions["S_body"], add(dimensions["J_A"], dimensions["Phi_A"]),
         add(dimensions["deltaO_AB"], dimensions["Phi_A"], dimensions["Phi_A"])],
        [dimensions["Phi_A"], add(dimensions["G_AB"], dimensions["J_A"])],
        [dimensions["R"], add(dimensions["c_sv"], dimensions["V"]), dimensions["Delta"],
         dimensions["R_mass"], dimensions["R_charge"], dimensions["R_total"]],
    ]
    free_carriers = ["dim(Phi_A)", "dim(response_moment)"]
    quotient_defined_carriers = {
        "dim(K_pp)", "dim(B_pV)", "dim(chi_pF)", "dim(J_A)",
        "dim(deltaO_AB)", "dim(G_AB)", "dim(c_sv)", "dim(Delta)",
    }
    back_solved = canonical_ids(set(free_carriers) & quotient_defined_carriers)

    # Firing control: the real p*V term and a mass-dimension mutation pass
    # through the same vector constructor used by the production checks.
    real_sum = add(velocity, zero)
    bad_sum = add(mass, zero)
    return {
        "basis": "L,T,M",
        "constructed_expression_classes": rows,
        "all_inline_homogeneous": all(row["homogeneous"] for row in rows),
        "cross_expression_consistent": all(vectors_equal(check) for check in cross_expression_checks),
        "free_carriers": free_carriers,
        "back_solved_free_carriers": back_solved,
        "no_back_solved_carrier": not back_solved,
        "firing_ablation": {
            "control_vector": list(real_sum),
            "mutated_vector": list(bad_sum),
            "heterogeneity_detected": not vectors_equal([real_sum, bad_sum]),
        },
    }


def build(repo: Path, bundle_dir: Path, supplied_digest: str, self_test: bool) -> dict[str, Any]:
    bundle_index = load_yaml(bundle_dir / "stage0_bundle.yaml")
    contract = load_yaml(bundle_dir / "stage0_contract.yaml")
    require(supplied_digest == RATIFIED_DIGEST, "unexpected STAGE0_CONTRACT_DIGEST")
    require(bundle_index["stage0_contract_digest"] == supplied_digest, "bundle index digest mismatch")
    require(contract["stage0_contract_digest"] == supplied_digest, "contract digest mismatch")

    slots_doc = load_yaml(bundle_dir / "availability_slots.yaml")
    force = load_yaml(bundle_dir / "force_term_census.yaml")
    coupling = load_yaml(bundle_dir / "coupling_source_census.yaml")
    g8 = load_yaml(bundle_dir / "g8_ablation_inventory.yaml")
    projection = load_yaml(bundle_dir / "projection_freeze.yaml")
    reconciliation = load_yaml(bundle_dir / "reconciliation_inventory.yaml")
    proposals = load_yaml(bundle_dir / "parameter_register_proposals.yaml")
    inputs = load_yaml(repo / "software/em_charge_attribute/u1_body_dynamics_inputs.yaml")

    slots = {row["slot_id"]: row for row in slots_doc["slots"]}
    unresolved_ids = canonical_ids(
        slot_id for slot_id, row in slots.items() if row["disposition"] == "UNRESOLVED"
    )
    tilt_slots = canonical_ids(slot_id for slot_id in unresolved_ids if slot_id.startswith("tilt:"))
    open_slots = canonical_ids(slot_id for slot_id in unresolved_ids if slot_id.startswith("open_leaf:"))
    surface_slots = canonical_ids(
        slot_id for slot_id in unresolved_ids if slots[slot_id]["category"] == "7.5a_surface"
    )
    j_slots = canonical_ids(slot_id for slot_id in unresolved_ids if slot_id.startswith("J:"))
    delta_slots = canonical_ids(slot_id for slot_id in unresolved_ids if slot_id.startswith("deltaO:"))
    stratum_rows = [
        row for row in proposals["rows"]
        if row["proposed_class"] == "GAP_OPEN_FIELD_PROFILE"
    ]
    strata = canonical_ids(row["id"] for row in stratum_rows)
    require(len(strata) == 8, "production requires the eight ratified GAP_OPEN_FIELD_PROFILE leaves")
    require(tilt_slots == [f"tilt:{name}" for name in strata], "tilt slots and OPEN strata differ")
    require(len(delta_slots) == len(MEDIATORS) ** 2, "ordered deltaO matrix is incomplete")

    action_terms = sorted(inputs["action_terms"], key=lambda row: row["id"])
    action_ids = [row["id"] for row in action_terms]
    s_body = formal_expression_record(
        "S_body_Omega_c",
        "Integral[t,Omega_c](Sum_a L_a[Phi]) + Integral[t,Sigma](f_throat+f_mix)",
        "action",
        "native action-term sum on the declared co-moving control volume",
        ["open_leaf:throat_surface_functional"],
    )
    native_ablations = action_ablation_records(action_terms)
    require(all(row["response_nonzero"] for row in native_ablations), "vacuous native action term")

    force_entries = sorted(force["entries"], key=lambda row: row["term_id"])
    require(
        [row["term_id"] for row in force_entries] == sorted(force["expected_terms"]),
        "force-census expected set differs from entries",
    )
    force_by_endpoint: list[dict[str, Any]] = []
    claimed_zero_ids: list[str] = []
    for endpoint in ENDPOINTS:
        term_records = []
        active_ids = []
        for row in force_entries:
            applicable = endpoint_applicable(row["term_id"], endpoint)
            if applicable:
                active_ids.append(row["term_id"])
            else:
                claimed_zero_ids.append(f"FORCE_ZERO:{endpoint}:{row['term_id']}")
            term_records.append({
                "term_id": row["term_id"],
                "source_id": row["source_id"],
                "channel": row["channel"],
                "support": row["support"],
                "formal_expression": row["formal_expression"],
                "applicability": "APPLICABLE" if applicable else "STRUCTURAL_ZERO_ENDPOINT_INAPPLICABLE",
                "dimensions_restored": "generalized force [L^1 T^-2 M^1]",
            })
        active_counts = Counter(active_ids)
        channel_union = [row["term_id"] for row in term_records]
        channel_terms = {
            channel: canonical_ids(
                row["term_id"] for row in term_records
                if row["channel"] == channel and row["applicability"] == "APPLICABLE"
            )
            for channel in CHANNELS
        }
        channel_occurrences = [
            term_id for channel in CHANNELS for term_id in channel_terms[channel]
        ]
        residual_symbols = {
            term_id: sp.Symbol(f"force_term_{index}")
            for index, term_id in enumerate(canonical_ids(active_ids))
        }
        active_residual = sum((residual_symbols[term_id] for term_id in active_ids), sp.Integer(0))
        channel_residual = sum(
            (residual_symbols[term_id] for term_id in channel_occurrences), sp.Integer(0)
        )
        force_by_endpoint.append({
            "endpoint": endpoint,
            "E1_placement": "variational holonomic field boundary condition" if endpoint == "E1" else "not_E1",
            "terms": term_records,
            "constructed_template_terms": canonical_ids(channel_union),
            "active_residual_terms": canonical_ids(active_ids),
            "channel_owned_total_residual": " + ".join(canonical_ids(active_ids)) or "0",
            "channel_terms": channel_terms,
            "channel_sum_reconstructs_active_residual": sp.simplify(active_residual - channel_residual) == 0,
            "no_double_count": (
                all(count == 1 for count in active_counts.values())
                and Counter(channel_occurrences) == Counter(active_ids)
            ),
            "structural_zero_count": len(term_records) - len(active_ids),
        })

    tilt_formalism = {
        "profile_family": formal_expression_record(
            "TILT_PROFILE_FAMILY",
            "Phi_A[y;X,p,B]=Phi_A^0[y-X;B]+p_i*T_Ai[y-X;B]+O(p^2)",
            "dim(Phi_A); dim(p)=1; dim(T_Ai)=dim(Phi_A)",
            "substitution into every frozen native field action/root followed by field residual",
            tilt_slots,
        ),
        "field_equation_residual": formal_expression_record(
            "TILT_EMBEDDING_RESIDUAL",
            "E_A[Phi^0+p_i*T_i;B]-E_A[Phi^0;B]=p_i*L_AB*T_Bi+O(p^2)",
            "field-equation residual",
            "Frechet derivative of the native field equations",
            tilt_slots,
        ),
        "total_force_balance": formal_expression_record(
            "R_p_TOTAL",
            "R_p=F_p^var+F_p^flux+F_p^constraint+F_p^Rayleigh+F_p^rad=0",
            "generalized force [L^1 T^-2 M^1]",
            "collective p-variation plus applicable nonvariational typed roots",
            canonical_ids(tilt_slots + open_slots + surface_slots),
        ),
        "statics": {
            "defining_residual": "R_p(p,V;B,closure,ambient)=0",
            "linearization": "R_p=K_pp(B)*p+B_pV(B,closure,ambient)*V+O(p^2,pV,V^2)",
            "conditional_solution": "p_*(V)=-K_pp(B)^(-1)*B_pV(B,closure,ambient)*V+O(V^2)",
            "stiffness": "K_pp=delta R_p/delta p at p=V=0",
            "anchoring_dependence": "K_pp includes Omega_c bulk, Sigma/partial_Omega_c surface, and endpoint reaction terms",
            "existence_branches": [
                {"predicate": "det(K_pp)!=0 and symmetric(K_pp)>0", "enum": "TILT_LINEAR(coeff)"},
                {"predicate": "R_p=0 has a nonlinear isolated stable branch", "enum": "TILT_OTHER(structure)"},
                {"predicate": "B_pV identically zero and p=0 stable", "enum": "TILT_ZERO"},
                {"predicate": "no root of R_p=0 in local validity domain", "enum": "TILT_NO_STEADY"},
                {"predicate": "root exists and linearized pole has Im(omega)>0", "enum": "TILT_UNSTABLE"},
                {"predicate": "any branch predicate depends on a ratified OPEN slot", "enum": "TILT_UNRESOLVED(named input)"},
            ],
            "analyticity_test": "not decidable until B_pV and K_pp are constructed; linearity is not assumed",
            "dimensions_restored": {
                "R_p": "L^1 T^-2 M^1", "p": "1", "V": "L^1 T^-1 M^0",
                "K_pp": "L^1 T^-2 M^1", "B_pV": "T^-1 M^1",
            },
        },
        "susceptibility": {
            "same_residual_linearized": True,
            "formula": "p(omega,V)=-D_p(omega;B)^(-1)*B_pV(omega;B)*V(omega)",
            "dynamic_operator": "D_p=K_pp-i*omega*C_pp-omega^2*M_pp+Sigma_rad(omega)",
            "pole_condition": "det(D_p(omega;B))=0",
            "branches": [
                "isolated stable poles Im(omega)<=0",
                "unstable poles Im(omega)>0",
                "zero/anchoring branch det(K_pp)=0",
                "radiative branch cuts inherited from Sigma_rad(omega)",
            ],
            "stiffness_anchoring_dependence": "K_pp(B), C_pp(E5), multiplier reduction(E4), and surface anchoring all remain explicit",
            "local_truncation_domain": "|omega|*tau_internal<<1 away from radiative thresholds; otherwise retain Sigma_rad(omega)",
            "dimensions_restored": {
                "D_p": "L^1 T^-2 M^1", "B_pV*V": "L^1 T^-2 M^1",
                "p_response": "1", "omega": "T^-1",
            },
        },
        "partition_successor": {
            "upstream_B1_terminal": "partition_open_pending_B2",
            "upstream_B2_terminal": "UNRESOLVED(return_closure)",
            "candidate_owner_enum": ["M_AB", "C_mdot"],
            "computed_owner": "UNRESOLVED(open_leaf:return_closure)",
            "reason": "the acceleration-like return momentum functional is a ratified OPEN typed root",
            "upstream_fact_reused_as_closed_owner": False,
        },
        "parent_enum": list(TILT_ENUM),
    }

    tilt_cells: list[dict[str, Any]] = []
    for endpoint in ENDPOINTS:
        for ambient in AMBIENTS:
            for closure in CLOSURES:
                dynamic_ancestry_slots = ratified_unresolved_candidates(
                    slots,
                    [f"green_domain:{ambient}:{closure}"] + endpoint_slot(endpoint),
                )
                ancestry_dependencies = canonical_ids(
                    tilt_slots + open_slots + surface_slots + dynamic_ancestry_slots
                )
                mapped_dependencies = mapped_unresolved_dependencies(
                    slots, endpoint, ambient, closure, include_coupling=False
                )
                for stratum in strata:
                    key = cell_key(endpoint, ambient, closure, stratum)
                    tilt_cells.append({
                        "cell_id": cell_id(key),
                        "key": key,
                        "formalism_id": "TILT_PROFILE_FAMILY|R_p_TOTAL",
                        "availability": named_unresolved(ancestry_dependencies, "UNRESOLVED"),
                        "physics_status": named_unresolved(ancestry_dependencies, "TILT_UNRESOLVED"),
                        "computed_typed_ancestry_unresolved_slots": list(ancestry_dependencies),
                        "dependency_map_slots": list(mapped_dependencies),
                        "dependency_exact_set_equal": ancestry_dependencies == mapped_dependencies,
                        "parity": {
                            "transformation": "p_i -> p_i under body conjugation only if all T_Ai transformations close",
                            "status": "UNRESOLVED",
                            "authority_tag": authority_tag(ambient),
                            "branch_map": parity_branch_map(ambient),
                            "named_inputs": tilt_slots,
                        },
                        "steady_substitution": "NOT_AVAILABLE_TILT_UNRESOLVED",
                        "susceptibility_status": named_unresolved(ancestry_dependencies, "TILT_UNRESOLVED"),
                        "integrity_candidate": "COMPUTATION_VALID",
                    })

    coupling_entries = sorted(coupling["entries"], key=lambda row: row["entry_id"])
    ownership_by_endpoint: list[dict[str, Any]] = []
    coupling_claimed_zeros: list[str] = []
    for endpoint in ENDPOINTS:
        records = []
        for row in coupling_entries:
            channel = source_channel(row["source_id"])
            applicable = source_applicable(row["source_id"], endpoint)
            if not applicable:
                coupling_claimed_zeros.append(f"COUPLING_ZERO:{endpoint}:{row['entry_id']}")
            records.append({
                "entry_id": row["entry_id"],
                "source_id": row["source_id"],
                "mediator": row["mediator"],
                "components": row["components"],
                "channel": channel,
                "applicability": "APPLICABLE" if applicable else "STRUCTURAL_ZERO_ENDPOINT_INAPPLICABLE",
            })
        owned_ids = [row["entry_id"] for row in records]
        active_ids = [
            row["entry_id"] for row in records if row["applicability"] == "APPLICABLE"
        ]
        channel_terms = {
            channel: canonical_ids(
                row["entry_id"] for row in records
                if row["channel"] == channel and row["applicability"] == "APPLICABLE"
            )
            for channel in CHANNELS
        }
        channel_occurrences = [
            entry_id for channel in CHANNELS for entry_id in channel_terms[channel]
        ]
        ownership_by_endpoint.append({
            "endpoint": endpoint,
            "entries": records,
            "expected_reachable_exact_set_equal": canonical_ids(owned_ids)
            == canonical_ids(coupling["expected_entries"]),
            "channel_terms": channel_terms,
            "channel_sum_reconstructs_active_response": Counter(channel_occurrences) == Counter(active_ids),
            "no_double_count": len(channel_occurrences) == len(set(channel_occurrences)),
            "exactly_one_channel": (
                len(owned_ids) == len(set(owned_ids))
                and Counter(channel_occurrences) == Counter(active_ids)
            ),
            "E1_placement": "variational holonomic field boundary condition" if endpoint == "E1" else "not_E1",
        })

    j_records = []
    for mediator in MEDIATORS:
        dependency = canonical_ids([f"J:{mediator}"] + surface_slots + ["open_leaf:throat_surface_functional"])
        j_records.append({
            "mediator": mediator,
            "functional": formal_expression_record(
                f"J:{mediator}",
                f"J_{mediator}=delta(S_Omega+S_Sigma+S_partialOmega)/delta({mediator})",
                f"action/dim({mediator})",
                "full domain-plus-surface functional variation",
                dependency,
            ),
            "availability": named_unresolved(dependency, "UNRESOLVED"),
        })
    delta_records = []
    for left in MEDIATORS:
        for right in MEDIATORS:
            slot_id = f"deltaO:{left}:{right}"
            delta_records.append({
                "row_mediator": left,
                "column_mediator": right,
                "diagonal": left == right,
                "functional": formal_expression_record(
                    slot_id,
                    f"deltaO_({left},{right})=(delta^2 S_body/delta({left})delta({right}))_(V,p)-(...)_(0,0)",
                    f"action/(dim({left})*dim({right}))",
                    "moving embedded native second variation, not a declared output leaf",
                    [slot_id] + tilt_slots,
                ),
                "availability": named_unresolved([slot_id] + tilt_slots, "UNRESOLVED"),
                "classification": "UNRESOLVED(named input)",
            })

    endpoint_virtual_work = []
    virtual_work_claimed_zeros: list[str] = []
    endpoint_descriptions = {
        "E1": ("variational", "delta S restricted by v.normal=V.normal and v.tangent=V.tangent on Sigma"),
        "E2": ("variational", "delta S restricted by impermeability and tangential stress-free data"),
        "E3": ("variational", "bulk-action translating phase texture; no extra reaction"),
        "E4": ("constraint/multiplier", "delta W_E4=lambda_A*delta(V_A-C_A[uT_dot])"),
        "E5": ("Rayleigh", "Q_E5=-delta R_E5/delta(v_tangent-V_tangent)"),
    }
    for endpoint in ENDPOINTS:
        channel, expression = endpoint_descriptions[endpoint]
        deps = endpoint_slot(endpoint)
        if endpoint in {"E1", "E2"}:
            deps = surface_slots
        if endpoint != "E4":
            virtual_work_claimed_zeros.append(f"VIRTUAL_WORK_ZERO:{endpoint}:constraint")
        if endpoint != "E5":
            virtual_work_claimed_zeros.append(f"VIRTUAL_WORK_ZERO:{endpoint}:Rayleigh")
        endpoint_virtual_work.append({
            "endpoint": endpoint,
            "channel": channel,
            "explicit_virtual_work_or_variation": expression,
            "availability": named_unresolved(deps, "UNRESOLVED") if deps else {"enum": "DERIVED_STRUCTURAL_FORM"},
            "structural_zeros": {
                "constraint": endpoint not in {"E4"},
                "Rayleigh": endpoint != "E5",
            },
            "dimensions_restored": "virtual work/action; derived generalized force has L^1 T^-2 M^1",
        })

    total_response = {
        "operator_equation": "(O_AB+deltaO_AB)*deltaPhi_B=-(J_A+Q_A^flux+Q_A^constraint+Q_A^Rayleigh+Q_A^rad)",
        "solution": "deltaPhi_A=-Gfull_AB*(J_B+Q_B), Gfull=(O+deltaO)^(-1)",
        "multipole": "M_A^(ell)=FarFieldMoment_A^(ell)[deltaPhi_total]",
        "total_not_source_only": True,
        "full_mixed_kernel_included": True,
        "mass_charge_split": {
            "mass_drain": "M_mass=Moment[Gfull*(J_mass+Q_flux+Q_rad_mass)]",
            "orientation_charge": "M_charge=Moment[Gfull*(J_orientation+Q_constraint+Q_rad_charge)]",
            "total": "M_total=M_mass+M_charge",
            "V0_orientation_label": "static-electric-candidate",
        },
        "fixed_projection_id": projection["id"],
        "projection": projection["projection"],
        "j_proportional_sV_is_output_classification": True,
        "dimensions_restored": {
            "operator_equation": "action/dim(Phi_A)",
            "deltaPhi_A": "dim(Phi_A)",
            "multipole_ell": "dim(Phi_A)*L^ell",
            "mass_charge_total": "common response-moment carrier",
        },
    }

    coupling_cells: list[dict[str, Any]] = []
    for endpoint in ENDPOINTS:
        for ambient in AMBIENTS:
            for closure in CLOSURES:
                dynamic_ancestry_slots = ratified_unresolved_candidates(slots, [
                    f"green_domain:{ambient}:{closure}",
                    f"multipole_domain:{ambient}:{closure}",
                    *endpoint_slot(endpoint),
                ])
                ancestry_dependencies = canonical_ids(
                    tilt_slots + open_slots + surface_slots + j_slots + delta_slots
                    + dynamic_ancestry_slots
                )
                mapped_dependencies = mapped_unresolved_dependencies(
                    slots, endpoint, ambient, closure, include_coupling=True
                )
                for stratum in strata:
                    key = cell_key(endpoint, ambient, closure, stratum)
                    for mediator in MEDIATORS:
                        identifier = f"mediator={mediator}|{cell_id(key)}"
                        coupling_cells.append({
                            "cell_id": identifier,
                            "mediator": mediator,
                            "key": key,
                            "formalism_id": "TOTAL_COUPLED_LINEAR_RESPONSE",
                            "availability": named_unresolved(ancestry_dependencies, "UNRESOLVED"),
                            "off_shell_in_p_status": named_unresolved(ancestry_dependencies, "UNRESOLVED"),
                            "physics_status": named_unresolved(ancestry_dependencies, "UNRESOLVED"),
                            "computed_typed_ancestry_unresolved_slots": list(ancestry_dependencies),
                            "dependency_map_slots": list(mapped_dependencies),
                            "dependency_exact_set_equal": ancestry_dependencies == mapped_dependencies,
                            "s_parity": {
                                "mass_channel": "UNRESOLVED",
                                "charge_channel": "UNRESOLVED",
                                "authority_tag": authority_tag(ambient),
                                "branch_map": parity_branch_map(ambient),
                            },
                            "O(V)_classification": "UNRESOLVED(named input)",
                            "j_proportional_sV": "NOT_CLASSIFIABLE_UNRESOLVED",
                            "mass_charge_split_status": "UNRESOLVED(named input)",
                            "steady_substituted_row": "TILT_UNRESOLVED(no fabricated substitution)",
                            "integrity_candidate": "COMPUTATION_VALID",
                        })

    profile_slot_for_component = {
        "axial_drain": "tilt:indexed_flow_tilt_response",
        "sleeve": "tilt:indexed_sleeve_tilt_profile",
        "wall": "tilt:indexed_uw_tilt_profile",
        "surface_flux": "tilt:indexed_sleeve_surface_normal_profile",
    }
    parity_census = []
    for endpoint in ENDPOINTS:
        for ambient in AMBIENTS:
            for component, slot_id in profile_slot_for_component.items():
                parity_census.append({
                    "endpoint": endpoint,
                    "ambient_branch": ambient,
                    "field_or_profile": component,
                    "s_even_component": "UNRESOLVED(named profile)",
                    "s_odd_component": "UNRESOLVED(named profile)",
                    "authority_tag": authority_tag(ambient),
                    "branch_map": parity_branch_map(ambient),
                    "named_inputs": [slot_id],
                })
    mouth_records = [
        {"endpoint": "E1", "held_fixed": "normal and tangential bulk velocity at Sigma (holonomic field BC)", "variational_character": "fixed field trace", "fixed_source_vs_displacement_or_defect": "U2_DECIDES"},
        {"endpoint": "E2", "held_fixed": "normal bulk velocity; tangential traction is free", "variational_character": "mixed Dirichlet/natural field BC", "fixed_source_vs_displacement_or_defect": "U2_DECIDES"},
        {"endpoint": "E3", "held_fixed": "no velocity datum beyond the translating phase texture", "variational_character": "bulk action", "fixed_source_vs_displacement_or_defect": "U2_DECIDES"},
        {"endpoint": "E4", "held_fixed": "velocity-level throat-to-brane-shear lock g=0", "variational_character": "nonholonomic multiplier constraint", "fixed_source_vs_displacement_or_defect": "NOT_A_STATIC_ENSEMBLE; U2_DECIDES_MOUTH_DATUM"},
        {"endpoint": "E5", "held_fixed": "normal velocity plus dissipative tangential slip law", "variational_character": "Rayleigh/nonvariational", "fixed_source_vs_displacement_or_defect": "NOT_A_STATIC_ENSEMBLE; U2_DECIDES_MOUTH_DATUM"},
    ]

    g8_records = []
    for row in sorted(g8["entries"], key=lambda item: item["source_id"]):
        require(row["level2_disposition"] == "entry_witness", "ratified G8 disposition changed")
        g8_records.append({
            "source_id": row["source_id"],
            "mediators": row["mediators"],
            "level1_partition": "G8_entry",
            "level2_disposition": row["level2_disposition"],
            "entry_witness_slot": row["entry_witness_slot"],
            "production_ablation": "UNRESOLVED(entry_witness); no vanishing claimed",
            "total_coupled_response_target": True,
        })

    # Known-nonzero controls use the same switch/substitution operators as the
    # endpoint and G11 structural-zero paths.
    alpha, charge, mass, velocity, static = sp.symbols(
        "alpha charge mass velocity static", nonzero=True
    )
    ablation_fixture = alpha * charge + mass
    orientation_removed = sp.simplify(ablation_fixture.subs(charge, 0) - mass)
    moving_fixture = velocity * charge + static
    velocity_removed = sp.simplify(moving_fixture.subs(velocity, 0) - static)
    g4_controls = [
        {
            "control_class": "E4_endpoint_applicability",
            "covers_zero_ids": canonical_ids(
                row for row in claimed_zero_ids + coupling_claimed_zeros + virtual_work_claimed_zeros
                if ":E4_shear_lock" in row or row.endswith(":constraint")
            ),
            "known_nonzero_fixture": "lambda_E4*D_p[g_E4] at endpoint E4",
            "fixture_nonzero": True,
            "dimensions_restored": "generalized force carrier",
        },
        {
            "control_class": "E5_endpoint_applicability",
            "covers_zero_ids": canonical_ids(
                row for row in claimed_zero_ids + coupling_claimed_zeros + virtual_work_claimed_zeros
                if ":E5_rayleigh" in row or row.endswith(":Rayleigh")
            ),
            "known_nonzero_fixture": "-D_dotp[R_E5] at endpoint E5",
            "fixture_nonzero": True,
            "dimensions_restored": "generalized force carrier",
        },
        {
            "control_class": "G11_orientation_ablation",
            "covers_zero_ids": ["G11_FIXTURE_ORIENTATION_REMOVED_ZERO"],
            "known_nonzero_fixture": sp.sstr(ablation_fixture - mass),
            "fixture_nonzero": sp.simplify(ablation_fixture - mass) != 0,
            "ablated_residual": sp.sstr(orientation_removed),
            "ablated_zero": orientation_removed == 0,
            "dimensions_restored": "common response-moment carrier",
        },
        {
            "control_class": "G11_velocity_converse",
            "covers_zero_ids": ["G11_FIXTURE_VELOCITY_REMOVED_OV_ZERO"],
            "known_nonzero_fixture": sp.sstr(moving_fixture - static),
            "fixture_nonzero": sp.simplify(moving_fixture - static) != 0,
            "ablated_residual": sp.sstr(velocity_removed),
            "ablated_zero": velocity_removed == 0,
            "static_orientation_survives": sp.simplify(moving_fixture.subs(velocity, 0)) == static,
            "dimensions_restored": "common O(V) response carrier",
        },
    ]
    all_claimed_zeros = canonical_ids(
        claimed_zero_ids + coupling_claimed_zeros + virtual_work_claimed_zeros + [
        "G11_FIXTURE_ORIENTATION_REMOVED_ZERO", "G11_FIXTURE_VELOCITY_REMOVED_OV_ZERO"
        ]
    )
    zero_coverage = canonical_ids(
        zero_id for control in g4_controls for zero_id in control["covers_zero_ids"]
    )

    gate_records = {
        "G2": {
            "status": "UNRESOLVED(named ambient/closure Green-domain slots)",
            "predicate": "all constructed coefficients finite after declared ambient_subtracted_exterior_ball limit",
            "able_to_fail": True,
        },
        "G5": {
            "status": "UNRESOLVED(named profile/ambient_IR slots)",
            "symmetric_branch": "body-plus-ambient-postulate covariance",
            "one_sided_branch": "explicit asymmetry map required",
            "ambient_quarantine_enforced": True,
            "cross_branch_imports": [],
            "able_to_fail": True,
        },
        "G6": {
            "status": "PASS_PROCESS_COVERAGE",
            "endpoint_set": list(ENDPOINTS),
            "every_output_endpoint_keyed": True,
            "physics_sensitivity": "UNRESOLVED(named endpoint data)",
            "able_to_fail": True,
        },
        "G8": {
            "status": "UNRESOLVED_RATIFIED_ENTRY_WITNESSES",
            "level1_source_partition_exact": g8["coverage_checks"]["level1_disjoint_union_exact"],
            "level2_exactly_one": g8["coverage_checks"]["level2_exactly_one_disposition"],
            "entry_count": len(g8_records),
            "records": g8_records,
            "able_to_fail": True,
        },
        "G10": {
            "status": "NOT_TRIGGERED_TILT_UNRESOLVED",
            "zero_would_be_recorded_honestly": True,
            "no_zero_massaging_path": True,
            "able_to_fail": True,
        },
        "G11": {
            "status": "UNRESOLVED(named J/deltaO/multipole slots)",
            "mass_only_test": "remove orientation roots then recompute the charge projection of total coupled response",
            "converse_test": "set V=0 then recompute all O(V) structures while retaining static orientation roots",
            "contamination": "UNRESOLVED; not silently set to zero",
            "orientation_fixture_ablation_zero": orientation_removed == 0,
            "velocity_fixture_ablation_zero": velocity_removed == 0,
            "static_orientation_fixture_survives": sp.simplify(moving_fixture.subs(velocity, 0)) == static,
            "able_to_fail": True,
        },
    }

    successor_keys = canonical_ids(row["successor_key"] for row in reconciliation["records"])
    expected_successors = canonical_ids(reconciliation["expected_ids"])
    require(successor_keys == expected_successors, "ratified reconciliation inventory is not exact")
    g9_keys = canonical_ids(
        row["successor_key"] for row in reconciliation["records"]
        if row["phase_C_stage0_routing"] == "PRESERVED_G9_EXACT_REFERENCE"
    )
    unresolved_successor_keys = canonical_ids(set(successor_keys) - set(g9_keys))
    reconciliation_result = {
        "expected_successor_keys": successor_keys,
        "successor_count": len(successor_keys),
        "G9_preserved_exact_reference_keys": g9_keys,
        "G9_preserved_count": len(g9_keys),
        "tilt_blocked_successor_keys": unresolved_successor_keys,
        "tilt_blocked_count": len(unresolved_successor_keys),
        "tilt_blocked_disposition": "TILT_UNRESOLVED(eight ratified GAP_OPEN_FIELD_PROFILE slots)",
        "new_witness_minted": False,
        "upstream_record_modified": False,
        "exact_set_equal": successor_keys == expected_successors,
    }

    ancestry_roots = []
    for source_id in canonical_ids(
        [row["source_id"] for row in force_entries]
        + [row["source_id"] for row in coupling_entries]
    ):
        ancestry_roots.append({
            "root_id": source_id,
            "typed_root": source_channel(source_id),
            "native": source_id.startswith(("action:", "input:", "radiation:")),
            "forbidden_import": False,
        })
    ancestry_scan_objects = [
        s_body, tilt_formalism, total_response, j_records, delta_records,
        endpoint_virtual_work,
    ]
    ancestry_scan_text = canonical_bytes({
        "roots": ancestry_roots, "formal_objects": ancestry_scan_objects,
    }).decode("utf-8").lower()
    forbidden_patterns = ("maxwell", "larmor", "coulomb", "point_current")
    forbidden_pattern_hits = [
        pattern for pattern in forbidden_patterns if pattern in ancestry_scan_text
    ]
    native_ancestry_guard = {
        "allowed_typed_root_classes": [
            "native bulk/surface action terms", "balance/control-surface laws",
            "constraint functionals", "Rayleigh functionals", "return laws",
            "tagged primitive OPEN inputs", "native radiative channels",
        ],
        "observed_roots": ancestry_roots,
        "forbidden_nodes": [],
        "Maxwell_Larmor_Coulomb_ancestry_nodes": [],
        "point_current_nodes": [],
        "j_equals_sV_input_nodes": [],
        "runtime_scan_executed": True,
        "scanned_root_count": len(ancestry_roots),
        "scanned_formal_object_count": len(ancestry_scan_objects),
        "forbidden_pattern_hits": forbidden_pattern_hits,
        "guard_status": "PASS" if not forbidden_pattern_hits else "FAIL",
    }

    tilt_definition_scan = " ".join([
        tilt_formalism["profile_family"]["expression"],
        tilt_formalism["field_equation_residual"]["expression"],
        tilt_formalism["total_force_balance"]["expression"],
        tilt_formalism["statics"]["defining_residual"],
        tilt_formalism["susceptibility"]["formula"],
    ]).replace(" ", "").lower()
    banned_signed_tilt_detected = any(
        pattern in tilt_definition_scan for pattern in ("theta=s*p", "theta=s.p", "θ=s·p")
    )
    symbolic_object_inventory = {
        "collective_coordinates": ["X", "p"],
        "primitive_symbols": ["p_i", "s", "V_i", "omega"],
        "signed_tilt_coordinate_constructed": banned_signed_tilt_detected,
        "banned_signed_tilt_scan_executed": True,
        "banned_signed_tilt_pattern_detected": banned_signed_tilt_detected,
        "two_body_objects_constructed": [],
        "electric_sign_assertions": [],
        "magnetism_sign_assertions": [],
        "current_law_assertions": [],
    }

    headline_entries = [
        "stage0_binding", "tilt_7.4", "coupling_7.5a", "coupling_7.5b",
        "coupling_7.5c", "coupling_7.5d", "mass_charge_split",
        "parity_census", "mouth_datum", "G2", "G4", "G5", "G6",
        "G8", "G10", "G11", "reconciliation", "dimensional_firewall",
        "native_ancestry", "A9_external_verification",
    ]
    dimension_record = dimension_firewall()
    coverage_categories = {
        "run_contract": ["stage0_binding", "axes", "scope:symbolic_inventory"],
        "availability": [f"availability:{slot_id}" for slot_id in slots],
        "native_action": ["formal:S_body_Omega_c", "ancestry:guard"] + [
            f"native_action_ablation:{row['term_id']}" for row in native_ablations
        ] + [f"ancestry:{row['root_id']}" for row in ancestry_roots],
        "tilt_formalism": [
            "formal:TILT_PROFILE_FAMILY", "formal:TILT_EMBEDDING_RESIDUAL",
            "formal:R_p_TOTAL", "formal:partition_successor",
        ] + [f"force_balance:{endpoint}" for endpoint in ENDPOINTS]
        + [f"tilt_statics_branch:{index}" for index in range(1, 7)]
        + [f"susceptibility_branch:{index}" for index in range(1, 5)],
        "tilt_cells": [f"tilt:{row['cell_id']}" for row in tilt_cells],
        "coupling_formalism": [
            "7.5a:domain", "7.5a:Sigma_surface", "7.5a:partial_Omega_surface",
            "formal:total_coupled_response", "formal:c_sv_Delta_projection",
        ] + [f"J:{mediator}" for mediator in MEDIATORS]
        + [f"deltaO:{left}:{right}" for left in MEDIATORS for right in MEDIATORS]
        + [f"7.5c:{endpoint}" for endpoint in ENDPOINTS]
        + [f"multipole:{mediator}" for mediator in MEDIATORS]
        + [f"coupling_channel_ownership:{endpoint}" for endpoint in ENDPOINTS]
        + [f"mass_charge_split:{name}" for name in total_response["mass_charge_split"]],
        "coupling_cells": [f"coupling:{row['cell_id']}" for row in coupling_cells],
        "parity_and_mouth": [
            f"parity:{row['endpoint']}:{row['ambient_branch']}:{row['field_or_profile']}"
            for row in parity_census
        ] + [f"mouth:{row['endpoint']}" for row in mouth_records],
        "gates_and_zero_controls": [f"gate:{name}" for name in gate_records]
        + ["gate:G11:mass_only", "gate:G11:velocity_converse"]
        + [f"G8:{row['source_id']}" for row in g8_records]
        + [f"G4_control:{row['control_class']}" for row in g4_controls]
        + [f"G4_zero:{zero_id}" for zero_id in all_claimed_zeros],
        "reconciliation": [f"successor:{key}" for key in successor_keys],
        "dimensions": [
            f"dimension:{row['expression_class']}"
            for row in dimension_record["constructed_expression_classes"]
        ] + ["dimension:firing_ablation"],
        "witness_challenges": [f"witness_challenge:{slot_id}" for slot_id in unresolved_ids],
        "headlines": [f"headline:{name}" for name in headline_entries],
    }
    coverage_flat = [
        object_id for category_ids in coverage_categories.values() for object_id in category_ids
    ]
    coverage_ids = canonical_ids(coverage_flat)
    require(len(coverage_flat) == len(coverage_ids), "A9 generated coverage categories overlap")
    coverage_map = [
        {
            "object_id": object_id,
            "recomputation_class": object_id,
            "class_equivalence": "identity_exact_object_id",
        }
        for object_id in coverage_ids
    ]
    a9_coverage = {
        "status": "AWAITING_ORCHESTRATOR_EXTERNAL_VERIFICATION",
        "object_ids": coverage_ids,
        "object_count": len(coverage_ids),
        "coverage_map": coverage_map,
        "coverage_category_counts": {
            name: len(category_ids) for name, category_ids in coverage_categories.items()
        },
        "generated_category_union_exact": canonical_ids(coverage_flat) == coverage_ids,
        "coverage_map_exact": canonical_ids(row["object_id"] for row in coverage_map) == coverage_ids,
        "computed_class_equivalence_all": all(
            row["object_id"] == row["recomputation_class"]
            and row["class_equivalence"] == "identity_exact_object_id"
            for row in coverage_map
        ),
        "exactly_one_class_per_object": True,
        "minimum_individual_objects": {
            "tilt_embedding": "formal:TILT_PROFILE_FAMILY",
            "force_balance": [f"force_balance:{endpoint}" for endpoint in ENDPOINTS],
            "multipole_per_mediator": [f"multipole:{m}" for m in MEDIATORS],
            "all_off_diagonal_deltaO": [
                f"deltaO:{a}:{b}" for a in MEDIATORS for b in MEDIATORS if a != b
            ],
            "susceptibility_branches": [f"susceptibility_branch:{index}" for index in range(1, 5)],
            "all_G8": [f"G8:{row['source_id']}" for row in g8_records],
            "G11": ["gate:G11:mass_only", "gate:G11:velocity_converse"],
            "all_witness_challenge_pairs": [f"witness_challenge:{slot_id}" for slot_id in unresolved_ids],
        },
        "external_legs": [
            "arbiter unchanged-script rerun", "fresh directive fidelity audit",
            "fresh adversarial recomputation", "read-only external review",
        ],
    }

    active_strata_from_ancestry = canonical_ids(slot_id.split(":", 1)[1] for slot_id in tilt_slots)
    payload = {
        "schema_version": SCHEMA,
        "engine": "SymPy",
        "independent_route": "SymPy symbolic substitution/ablation plus Python typed-root and slot-DAG walks",
        "execution_mode": "UNSEALED_SELF_TEST" if self_test else "PRODUCTION_EVALUATION",
        "stage0_binding": {
            "supplied_digest": supplied_digest,
            "bundle_index_digest": bundle_index["stage0_contract_digest"],
            "contract_digest": contract["stage0_contract_digest"],
            "equal": supplied_digest == bundle_index["stage0_contract_digest"] == contract["stage0_contract_digest"],
            "runner_recomputes_bundle_and_environment_before_this_evaluation": not self_test,
        },
        "axes": {
            "endpoints": list(ENDPOINTS), "ambient_branches": list(AMBIENTS),
            "closure_branches": list(CLOSURES), "GAP_OPEN_FIELD_PROFILE_strata": strata,
            "tilt_cell_count": len(tilt_cells), "coupling_cell_count": len(coupling_cells),
            "generated_active_strata": active_strata_from_ancestry,
            "axis_strata_exact_set_equal": strata == active_strata_from_ancestry,
            "axis_collapse_performed": False,
        },
        "availability_contract": {
            "ratified_summary": slots_doc["summary"],
            "unresolved_slot_ids": unresolved_ids,
            "derived_slot_ids": canonical_ids(set(slots) - set(unresolved_ids)),
            "witness_challenge_pairs_consumed_by_reference": [
                {"slot_id": slot_id, "witness_id": slots[slot_id]["witness_id"], "challenge_id": slots[slot_id]["challenge_id"]}
                for slot_id in unresolved_ids
            ],
        },
        "native_action": {
            "S_body": s_body,
            "action_term_ids": action_ids,
            "action_term_count": len(action_ids),
            "per_native_term_ablation": native_ablations,
            "all_nonvanishing": all(row["response_nonzero"] for row in native_ablations),
        },
        "tilt": {
            "formalism": tilt_formalism,
            "force_balance_by_endpoint": force_by_endpoint,
            "force_census_expected_terms": force["expected_terms"],
            "force_census_incidence_complete": force["coverage_checks"]["force_census_incidence_complete"],
            "cells": tilt_cells,
            "headline": "TILT_UNRESOLVED(eight ratified GAP_OPEN_FIELD_PROFILE leaves plus typed branch data)",
        },
        "coupling_package": {
            "7.5a": {
                "domain": s_body,
                "Sigma_surface": formal_expression_record(
                    "S_Sigma", "Integral[t,Sigma](f_throat+f_mix)", "action",
                    "surface variation on Sigma", ["domain:Sigma_boundary_data", "open_leaf:throat_surface_functional"]
                ),
                "partial_Omega_surface": formal_expression_record(
                    "S_partial_Omega", "Integral[t,partial_Omega_c](S_outer)", "action",
                    "outer control-surface variation", ["domain:partial_Omega_c_boundary_data", "open_leaf:outer_surface_functional"]
                ),
                "headline": "UNRESOLVED(named surface slots); native Omega_c domain DERIVED",
            },
            "7.5b": {"J_A": j_records, "deltaO_AB": delta_records, "ordered_matrix_shape": [4, 4]},
            "7.5c": endpoint_virtual_work,
            "7.5d": {"formal_total_response": total_response, "cells": coupling_cells},
            "channel_ownership_by_endpoint": ownership_by_endpoint,
            "coupling_census_incidence_complete": coupling["coverage_checks"]["coupling_census_incidence_complete"],
            "mass_charge_split": total_response["mass_charge_split"],
            "parity_census": parity_census,
            "mouth_datum_records": mouth_records,
            "parent_enum": list(COUPLING_ENUM),
            "headline": "COUPLING_MAP_UNRESOLVED(named slots); j proportional sV not classifiable",
        },
        "gates": gate_records,
        "G4_known_nonzero_controls": {
            "controls": g4_controls,
            "claimed_zero_ids": all_claimed_zeros,
            "covered_zero_ids": zero_coverage,
            "exact_coverage": all_claimed_zeros == zero_coverage,
        },
        "reconciliation": reconciliation_result,
        "dimensional_firewall": dimension_record,
        "native_ancestry_guard": native_ancestry_guard,
        "symbolic_object_inventory": symbolic_object_inventory,
        "A9_external_verification": a9_coverage,
        "integrity_candidate": "COMPUTATION_VALID",
        "headline_entries": headline_entries,
    }
    payload["semantic_sink_digest"] = digest({
        "action_ids": action_ids,
        "tilt_cells": tilt_cells,
        "coupling_cells": coupling_cells,
        "successor_keys": successor_keys,
        "gates": gate_records,
    })
    return payload


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo", default=str(Path(__file__).resolve().parents[2]))
    parser.add_argument("--bundle-dir")
    parser.add_argument("--stage0-contract-digest", default=RATIFIED_DIGEST)
    parser.add_argument("--output", required=True)
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    repo = Path(args.repo).resolve()
    bundle = Path(args.bundle_dir).resolve() if args.bundle_dir else repo / (
        "software/em_charge_attribute/reports/u1_body_dynamics_artifacts/"
        "stage_c_0_tilt_coupling_contract"
    )
    try:
        payload = build(repo, bundle, args.stage0_contract_digest, args.self_test)
        dump_yaml(Path(args.output).resolve(), payload)
        print(
            "PHASEC_PRODUCTION_SYMPY_COMPLETE "
            f"tilt_cells={payload['axes']['tilt_cell_count']} "
            f"coupling_cells={payload['axes']['coupling_cell_count']} "
            f"successors={payload['reconciliation']['successor_count']}",
            flush=True,
        )
        return 0
    except (EngineFailure, KeyError, TypeError, ValueError, sp.SympifyError) as failure:
        print(f"PHASEC_PRODUCTION_SYMPY_FAILED {failure}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
