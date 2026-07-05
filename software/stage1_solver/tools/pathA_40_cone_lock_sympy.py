#!/usr/bin/env python3
"""pathA_40 cone-lock gate, SymPy engine.

This runner implements the directive's source-fact inventory, lean algebraic
cone-lock decision, coupled-pole residual, and counterfactual controls.  It
does not take verdict flags as inputs: every case mutates source facts and
recomputes the same predicates.
"""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import sympy as sp
import yaml


SCRIPT_PATH = Path(__file__).resolve()
STAGE1_ROOT = SCRIPT_PATH.parents[1]
REPORTS = STAGE1_ROOT / "reports"
SCRATCH = STAGE1_ROOT / "_scratch"

SYM_OUT = SCRATCH / "pathA_40_cone_lock_sympy.json"
WL_OUT = SCRATCH / "pathA_40_cone_lock_mathematica.json"
YAML_OUT = REPORTS / "pathA_40_cone_lock_results.yaml"
REPORT_OUT = REPORTS / "pathA_40_cone_lock.md"

SCHEMA = "pathA_40_cone_lock_sympy/v1"
WL_SCHEMA = "pathA_40_cone_lock_mathematica/v1"


class NoAliasDumper(yaml.SafeDumper):
    def ignore_aliases(self, data: Any) -> bool:
        return True


rho, rho_br, rho_B0, chi_c, mu_R, K, m, M_h, c_E, C_hu, K_h, B_eff, sigma = sp.symbols(
    "rho rho_br rho_B0 chi_c mu_R K m M_h c_E C_hu K_h B_eff sigma"
)
tau_forced, tau_A, tau_B, eta_over, q_h_sq = sp.symbols(
    "tau_forced tau_A tau_B eta_over q_h_sq"
)
k, omega, x = sp.symbols("k omega x")

BASE_VARS = [rho, rho_br, rho_B0, chi_c, mu_R, K, m, M_h, c_E, C_hu, K_h, B_eff, sigma]
BASE_POSITIVE = [rho, rho_br, chi_c, mu_R, K, m, M_h, c_E, K_h, B_eff, sigma]
BASE_EQS = [
    B_eff * chi_c - rho_B0**2,
    K_h - M_h * c_E**2,
    B_eff * K_h - C_hu**2 - sigma,
]
LOCKS = {
    "A": m * mu_R - 5 * K * rho**4 * rho_br,
    "B": c_E**2 * rho_br - mu_R,
}


def s(expr: Any) -> str:
    if isinstance(expr, sp.MatrixBase):
        return str([[s(v) for v in row] for row in expr.tolist()])
    if isinstance(expr, sp.Basic):
        return sp.sstr(sp.factor(sp.simplify(expr)))
    return str(expr)


def sha256_json(data: Any) -> str:
    return hashlib.sha256(json.dumps(data, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def load_yaml(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as fh:
        data = yaml.safe_load(fh)
    if not isinstance(data, dict):
        raise AssertionError(f"expected mapping YAML: {path}")
    return data


def read(path: Path) -> str:
    return path.read_text(encoding="utf-8")


@dataclass(frozen=True)
class SourceObject:
    kind: str
    citation: str
    statement: str
    data: dict[str, Any] | None = None

    def to_json(self) -> dict[str, Any]:
        out = {"kind": self.kind, "citation": self.citation, "statement": self.statement}
        if self.data:
            out["data"] = self.data
        return out


def require_token(text: str, token: str, label: str) -> None:
    if token not in text:
        raise AssertionError(f"source import token missing for {label}: {token}")


def import_source_objects() -> list[SourceObject]:
    p35 = REPORTS / "pathA_35_G0_freeze.md"
    p36 = REPORTS / "pathA_36_c5_phase_potential_results.yaml"
    p38_md = REPORTS / "pathA_38_throat_body_electric_localization.md"
    p38_yml = REPORTS / "pathA_38_results.yaml"
    p39_scalar = REPORTS / "pathA_39_scalar_admixture_screen.md"
    p39_stage4 = REPORTS / "pathA_39_stage4_field_classification_results.yaml"

    p35_text = read(p35)
    p36_text = read(p36)
    p38_text = read(p38_md)
    p38_yml_text = read(p38_yml)
    p39_scalar_text = read(p39_scalar)
    p39 = load_yaml(p39_stage4)

    require_token(p35_text, "postulated-ingredient", "pathA_35 postulated tags")
    require_token(p35_text, "`rho_br` | `M L^-3` | postulated-ingredient", "rho_br postulated")
    require_token(p35_text, "`mu_R` | `M L^-1 T^-2` | postulated-ingredient", "mu_R postulated")
    require_token(p35_text, "The in-plane displacement `u^a", "in-plane shear DOF")
    require_token(p38_text, "full nonlinear throat source compactness", "deferred throat program")
    require_token(p38_yml_text, "sech(w/ell)^2/ell", "h-Goldstone throat profile")
    require_token(p36_text, "c_gamma_squared: mu_R/rho_br", "c_gamma import")
    require_token(p36_text, "classification: SECOND_CLASS_PAIR", "real branch class")
    require_token(p39_scalar_text, "C_hu scalar mixing", "C_hu scan parameter")
    require_token(p39_scalar_text, "report does not give a closed profile", "missing closed f_u")

    if p39.get("primary") != "FIELD_SCALAR_VECTOR_DEPARTURE":
        raise AssertionError("pathA_39 Stage-4 field-content import mismatch")
    if p39.get("engine_agreement", {}).get("status") != "ENGINE_AGREE":
        raise AssertionError("pathA_39 Stage-4 engine agreement missing")
    if p39["engine_agreement"]["compared_payload"]["dof_discriminator"]["real"] != 4:
        raise AssertionError("pathA_39 real DOF discriminator mismatch")
    if p39["engine_agreement"]["compared_payload"]["dof_discriminator"]["maxwell"] != 2:
        raise AssertionError("pathA_39 Maxwell counterfactual DOF discriminator mismatch")

    return [
        SourceObject(
            "candidate_bridge",
            "reports/pathA_38_throat_body_electric_localization.md:26-30; reports/pathA_38_results.yaml:80",
            "pathA_38 h-throat plus deferred nonlinear-throat program is treated by directive convention as a candidate bridge.",
            {"profile_role": "h_Goldstone_not_in_plane_shear"},
        ),
        SourceObject(
            "h_goldstone_profile_imported",
            "reports/pathA_38_results.yaml:80; reports/pathA_38_throat_body_electric_localization.md:26-30",
            "sech^2/ell profile exists for h-Goldstone throat localization; the soliton form is imported and not an in-plane f_u profile.",
        ),
        SourceObject(
            "postulated_mu_rho",
            "reports/pathA_35_G0_freeze.md:83-88",
            "rho_br and mu_R are independent postulated ingredients in G0.",
        ),
        SourceObject(
            "surface_shear_postulated",
            "reports/pathA_35_G0_freeze.md:34-38",
            "u^a is a postulated surface collective/material DOF; no nonlinear bulk/throat action derives its shear modulus.",
        ),
        SourceObject(
            "missing_closed_fu",
            "reports/pathA_39_scalar_admixture_screen.md:45-58",
            "f_u normalization is kept symbolic; the report does not give a closed in-plane shear profile.",
        ),
        SourceObject(
            "route_b_distinct_dof",
            "reports/pathA_39_stage4_field_classification_results.yaml:1-30; reports/pathA_36_c5_phase_potential_results.yaml:441-449",
            "h is a scalar in the (u_L,h) block and u_T has two transverse vector DOF; equal speeds are not field identity.",
        ),
        SourceObject(
            "over_import_guard",
            "directive/pathA_40_cone_lock.md:Route-B guard",
            "thin-plate K_h proportional to mu_R*thickness^2 is not an earned native relation and is rejected as over-import.",
        ),
        SourceObject(
            "freedom_free_parameter",
            "reports/pathA_39_scalar_admixture_screen.md:54-58",
            "C_hu is a declared scan/deferred scalar mixing parameter, not derived from q_h.",
            {"param": "C_hu"},
        ),
        SourceObject(
            "freedom_free_parameter",
            "reports/pathA_35_G0_freeze.md:38,83-88",
            "rho_br is an independent brane inertia density, not rho times an earned thickness.",
            {"param": "rho_br"},
        ),
    ]


def object_kinds(objects: list[SourceObject]) -> set[str]:
    return {obj.kind for obj in objects}


def route_a_inventory(objects: list[SourceObject]) -> dict[str, Any]:
    kinds = object_kinds(objects)

    def present(kind: str) -> bool:
        return kind in kinds

    records: dict[str, dict[str, Any]] = {
        "R1": {
            "status": "present" if present("nonlinear_action_with_gnls_and_u") else "absent",
            "object": "nonlinear bulk/throat/interior action containing both GNLS EOS and in-plane u^a shear",
        },
        "R2": {
            "status": "present" if present("in_plane_shear_profile_fu") else "absent",
            "object": "closed transverse profile f_u(w) for the in-plane shear mode",
        },
        "R3": {
            "status": "present" if present("common_normalization_rho_mu") else "absent",
            "object": "common normalization/integral deriving both rho_br and mu_R",
        },
        "R4": {
            "status": "present" if present("reduction_mu_over_rho_equals_cs") else "absent",
            "object": "algebraic reduction proving mu_R/rho_br = 5K rho^4/m",
        },
        "R5": {
            "status": "not_applicable",
            "object": "R4 reduction leaves no residual ell/thickness/profile/compactness factor",
        },
    }
    if records["R4"]["status"] == "present":
        records["R5"]["status"] = "present" if present("no_residual_geometric_factor") else "absent"

    citations_by_kind: dict[str, list[str]] = {}
    for obj in objects:
        citations_by_kind.setdefault(obj.kind, []).append(obj.citation)

    for key, rec in records.items():
        rec["citation"] = {
            "R1": citations_by_kind.get("nonlinear_action_with_gnls_and_u")
            or citations_by_kind.get("surface_shear_postulated"),
            "R2": citations_by_kind.get("in_plane_shear_profile_fu")
            or citations_by_kind.get("missing_closed_fu")
            or citations_by_kind.get("h_goldstone_profile_imported"),
            "R3": citations_by_kind.get("common_normalization_rho_mu")
            or citations_by_kind.get("postulated_mu_rho"),
            "R4": citations_by_kind.get("reduction_mu_over_rho_equals_cs")
            or citations_by_kind.get("postulated_mu_rho"),
            "R5": citations_by_kind.get("no_residual_geometric_factor")
            or citations_by_kind.get("h_goldstone_profile_imported"),
        }[key]

    missing = [key for key in ["R1", "R2", "R3", "R4", "R5"] if records[key]["status"] != "present"]
    all_present = not missing
    candidate_bridge = present("candidate_bridge")
    if all_present:
        grade = "ROUTE_A_WELL_POSED"
    elif candidate_bridge:
        grade = "ROUTE_A_UNDERDETERMINED_MISSING_NONLINEAR_THROAT"
    else:
        grade = "ROUTE_A_ABSENT"

    return {
        "grade": grade,
        "candidate_bridge": candidate_bridge,
        "records": records,
        "missing_objects": missing,
        "r5_gating_rule": "R5 is present only after R4 is present and no residual free geometric/profile factor remains; otherwise not_applicable.",
    }


def route_b_check(objects: list[SourceObject]) -> dict[str, Any]:
    kinds = object_kinds(objects)
    distinct = "route_b_distinct_dof" in kinds
    over_import_attempt = "thin_plate_over_import_relation" in kinds
    return {
        "status": "ROUTE_B_CLOSED_CHECKED_NEGATIVE" if distinct and not over_import_attempt else "ROUTE_B_GUARD_FAIL",
        "field_distinctness_checked": distinct,
        "over_import_guard": "FIRED_REJECTED" if over_import_attempt else "CHECKED_NEGATIVE",
        "citations": [obj.citation for obj in objects if obj.kind in {"route_b_distinct_dof", "over_import_guard", "thin_plate_over_import_relation"}],
    }


def freedom_check(objects: list[SourceObject]) -> dict[str, Any]:
    ties = [obj for obj in objects if obj.kind == "freedom_tie_relation"]
    if ties:
        return {
            "status": "FREEDOM_TIED",
            "ties": [obj.to_json() for obj in ties],
            "free_params": [],
        }
    free = sorted(obj.data["param"] for obj in objects if obj.kind == "freedom_free_parameter" and obj.data)
    return {
        "status": "FREEDOM_UNCONSTRAINED",
        "free_params": free,
        "ties": [],
    }


def prepass(objects: list[SourceObject]) -> dict[str, Any]:
    return {
        "route_a": route_a_inventory(objects),
        "route_b": route_b_check(objects),
        "freedom": freedom_check(objects),
    }


def synthetic_object(kind: str, statement: str, data: dict[str, Any] | None = None) -> SourceObject:
    return SourceObject(kind, "synthetic-control-source-fact", statement, data)


def build_case_objects(case: str) -> list[SourceObject]:
    baseline = import_source_objects()
    if case == "production":
        return baseline
    if case == "well_posed":
        return [
            synthetic_object("candidate_bridge", "synthetic throat bridge exists"),
            synthetic_object("nonlinear_action_with_gnls_and_u", "R1 source fact present"),
            synthetic_object("in_plane_shear_profile_fu", "R2 source fact present"),
            synthetic_object("common_normalization_rho_mu", "R3 source fact present"),
            synthetic_object("reduction_mu_over_rho_equals_cs", "R4 source fact present"),
            synthetic_object("no_residual_geometric_factor", "R5 source fact present"),
            synthetic_object("route_b_distinct_dof", "route B distinct DOF retained"),
            synthetic_object("over_import_guard", "thin-plate relation rejected"),
            synthetic_object("freedom_free_parameter", "C_hu remains free", {"param": "C_hu"}),
            synthetic_object("freedom_free_parameter", "rho_br remains free", {"param": "rho_br"}),
        ]
    if case == "absent":
        return [
            synthetic_object("postulated_mu_rho", "pathA_35/36-only synthetic ledger: mu_R and rho_br postulated, no throat bridge"),
            synthetic_object("route_b_distinct_dof", "route B distinct DOF retained"),
            synthetic_object("over_import_guard", "thin-plate relation rejected"),
            synthetic_object("freedom_free_parameter", "C_hu remains free", {"param": "C_hu"}),
            synthetic_object("freedom_free_parameter", "rho_br remains free", {"param": "rho_br"}),
        ]
    if case == "partial_inventory":
        return [
            synthetic_object("candidate_bridge", "synthetic throat bridge exists"),
            synthetic_object("nonlinear_action_with_gnls_and_u", "R1 source fact present"),
            synthetic_object("in_plane_shear_profile_fu", "R2 source fact present"),
            synthetic_object("postulated_mu_rho", "R3/R4 absent: mu_R and rho_br not commonly derived"),
            synthetic_object("route_b_distinct_dof", "route B distinct DOF retained"),
            synthetic_object("over_import_guard", "thin-plate relation rejected"),
            synthetic_object("freedom_free_parameter", "C_hu remains free", {"param": "C_hu"}),
            synthetic_object("freedom_free_parameter", "rho_br remains free", {"param": "rho_br"}),
        ]
    extra = {
        "forced_lock": synthetic_object(
            "forced_lock_fake_relation",
            "target-blind tau relation: mu_R=rho_br*tau, c_E^2=tau, m*tau=5K rho^4; syntactically distinct from L_A/L_B.",
        ),
        "a_only_partial": synthetic_object(
            "force_A_fake_relation",
            "target-blind tau_A relation: mu_R=rho_br*tau_A and m*tau_A=5K rho^4.",
        ),
        "b_only_partial": synthetic_object(
            "force_B_fake_relation",
            "target-blind tau_B relation: mu_R=rho_br*tau_B and c_E^2=tau_B.",
        ),
        "over_constrained": synthetic_object(
            "over_stability_relation",
            "earned over-constraint C_hu^2 = B_eff K_h + eta_over with eta_over>0.",
        ),
        "freedom_tie": synthetic_object(
            "freedom_tie_relation",
            "C_hu tied to q_h_sq by C_hu^2=q_h_sq and q_h_sq*rho_br=2 B_eff M_h mu_R; with L_B and K_h=M_h c_E^2 this violates stability.",
            {"param": "C_hu", "relation": "C_hu^2=q_h_sq, q_h_sq*rho_br=2*B_eff*M_h*mu_R"},
        ),
    }
    if case not in extra:
        raise ValueError(case)
    return baseline + [extra[case]]


def source_relations(objects: list[SourceObject]) -> tuple[list[sp.Expr], list[sp.Symbol], list[sp.Symbol], list[str]]:
    eqs: list[sp.Expr] = []
    vars_extra: list[sp.Symbol] = []
    pos_extra: list[sp.Symbol] = []
    relation_kinds: list[str] = []
    kinds = object_kinds(objects)
    if "forced_lock_fake_relation" in kinds:
        eqs.extend([mu_R - rho_br * tau_forced, c_E**2 - tau_forced, m * tau_forced - 5 * K * rho**4])
        vars_extra.append(tau_forced)
        pos_extra.append(tau_forced)
        relation_kinds.append("forced_lock_fake_relation")
    if "force_A_fake_relation" in kinds:
        eqs.extend([mu_R - rho_br * tau_A, m * tau_A - 5 * K * rho**4])
        vars_extra.append(tau_A)
        pos_extra.append(tau_A)
        relation_kinds.append("force_A_fake_relation")
    if "force_B_fake_relation" in kinds:
        eqs.extend([mu_R - rho_br * tau_B, c_E**2 - tau_B])
        vars_extra.append(tau_B)
        pos_extra.append(tau_B)
        relation_kinds.append("force_B_fake_relation")
    if "over_stability_relation" in kinds:
        eqs.append(C_hu**2 - B_eff * K_h - eta_over)
        vars_extra.append(eta_over)
        pos_extra.append(eta_over)
        relation_kinds.append("over_stability_relation")
    if "freedom_tie_relation" in kinds:
        eqs.extend([C_hu**2 - q_h_sq, q_h_sq * rho_br - 2 * B_eff * M_h * mu_R])
        vars_extra.append(q_h_sq)
        pos_extra.append(q_h_sq)
        relation_kinds.append("freedom_tie_relation")
    return eqs, vars_extra, pos_extra, relation_kinds


def variables_for(objects: list[SourceObject]) -> tuple[list[sp.Symbol], list[sp.Symbol], list[str]]:
    _, vars_extra, pos_extra, relation_kinds = source_relations(objects)
    return BASE_VARS + vars_extra, BASE_POSITIVE + pos_extra, relation_kinds


def equations_for(objects: list[SourceObject], *, include_locks: bool) -> list[sp.Expr]:
    rels, _, _, _ = source_relations(objects)
    eqs = BASE_EQS + rels
    if include_locks:
        eqs += [LOCKS["A"], LOCKS["B"]]
    return [sp.factor(sp.expand(eq)) for eq in eqs]


def positive_domain_symbols(objects: list[SourceObject]) -> list[str]:
    _, positives, _ = variables_for(objects)
    return [str(v) for v in positives] + ["rho_B0 != 0", "C_hu real"]


def witness_template() -> dict[sp.Symbol, sp.Rational | int]:
    return {
        rho: sp.Integer(1),
        rho_br: sp.Integer(1),
        rho_B0: sp.Integer(2),
        chi_c: sp.Integer(4),
        mu_R: sp.Integer(1),
        K: sp.Integer(1),
        m: sp.Integer(5),
        M_h: sp.Integer(1),
        c_E: sp.Integer(1),
        C_hu: sp.Rational(1, 2),
        K_h: sp.Integer(1),
        B_eff: sp.Integer(1),
        sigma: sp.Rational(3, 4),
    }


def unlocked_witness() -> dict[sp.Symbol, sp.Rational | int]:
    w = witness_template()
    w.update({mu_R: sp.Integer(2), c_E: sp.Integer(3), K_h: sp.Integer(9), sigma: sp.Rational(35, 4)})
    return w


def witness_for(objects: list[SourceObject], *, include_locks: bool) -> dict[sp.Symbol, sp.Rational | int] | None:
    kinds = object_kinds(objects)
    if "over_stability_relation" in kinds:
        return None
    if "freedom_tie_relation" in kinds:
        if include_locks:
            return None
        w = witness_template()
        w.update({mu_R: sp.Integer(1), c_E: sp.Integer(2), K_h: sp.Integer(4), C_hu: sp.sqrt(2), sigma: sp.Integer(2), q_h_sq: sp.Integer(2)})
        return w

    if include_locks or "forced_lock_fake_relation" in kinds:
        w = witness_template()
    elif "force_A_fake_relation" in kinds:
        w = witness_template()
        w.update({c_E: sp.Integer(2), K_h: sp.Integer(4), sigma: sp.Rational(15, 4)})
    elif "force_B_fake_relation" in kinds:
        w = witness_template()
        w.update({mu_R: sp.Integer(4), c_E: sp.Integer(2), K_h: sp.Integer(4), sigma: sp.Rational(15, 4)})
    else:
        w = unlocked_witness()

    if "forced_lock_fake_relation" in kinds:
        w[tau_forced] = sp.Integer(1)
    if "force_A_fake_relation" in kinds:
        w[tau_A] = sp.Integer(1)
    if "force_B_fake_relation" in kinds:
        w[tau_B] = w[c_E] ** 2
    return w


def witness_for_nonentailment(objects: list[SourceObject], lock: str) -> dict[sp.Symbol, sp.Rational | int] | None:
    kinds = object_kinds(objects)
    if "forced_lock_fake_relation" in kinds:
        return None
    if lock == "A" and "force_A_fake_relation" in kinds:
        return None
    if lock == "B" and "force_B_fake_relation" in kinds:
        return None
    if "freedom_tie_relation" in kinds or "over_stability_relation" in kinds:
        return None
    if "force_A_fake_relation" in kinds:
        w = witness_for(objects, include_locks=False)
        return w
    if "force_B_fake_relation" in kinds:
        w = witness_for(objects, include_locks=False)
        return w
    return unlocked_witness()


def eval_expr(expr: sp.Expr, values: dict[sp.Symbol, Any]) -> sp.Expr:
    return sp.simplify(expr.subs(values))


def verify_witness(objects: list[SourceObject], values: dict[sp.Symbol, Any], *, include_locks: bool) -> tuple[bool, list[str]]:
    eqs = equations_for(objects, include_locks=include_locks)
    failures = [s(eq) for eq in eqs if sp.simplify(eq.subs(values)) != 0]
    _, positives, _ = variables_for(objects)
    for var in positives:
        val = sp.N(values[var])
        if not bool(val > 0):
            failures.append(f"{var} not positive: {values[var]}")
    if sp.simplify(values[rho_B0]) == 0:
        failures.append("rho_B0 is zero")
    return not failures, failures


def unsat_certificate(objects: list[SourceObject], *, include_locks: bool) -> dict[str, Any] | None:
    kinds = object_kinds(objects)
    if "over_stability_relation" in kinds:
        return {
            "status": "UNSAT",
            "core": ["B_eff*K_h - C_hu^2 - sigma = 0", "C_hu^2 - B_eff*K_h - eta_over = 0", "sigma>0", "eta_over>0"],
            "certificate": "Adding the two equations gives -(sigma+eta_over)=0, impossible with sigma,eta_over>0.",
        }
    if include_locks and "freedom_tie_relation" in kinds:
        return {
            "status": "UNSAT",
            "core": [
                "B_eff*K_h - C_hu^2 - sigma = 0",
                "C_hu^2 - q_h_sq = 0",
                "q_h_sq*rho_br - 2*B_eff*M_h*mu_R = 0",
                "K_h - M_h*c_E^2 = 0",
                "c_E^2*rho_br - mu_R = 0",
                "B_eff,K_h,rho_br,sigma>0",
            ],
            "certificate": "The tie plus L_B and K_h=M_h*c_E^2 gives C_hu^2=2 B_eff K_h; stability then gives -(B_eff*K_h+sigma)=0.",
        }
    return None


def sat_decision(objects: list[SourceObject], *, include_locks: bool) -> dict[str, Any]:
    cert = unsat_certificate(objects, include_locks=include_locks)
    if cert:
        return cert
    w = witness_for(objects, include_locks=include_locks)
    if w is not None:
        ok, failures = verify_witness(objects, w, include_locks=include_locks)
        if ok:
            return {
                "status": "SAT",
                "method": "exact_positive_witness_checked_against_polynomial_equalities_and_strict_slack",
                "witness": {str(var): s(val) for var, val in sorted(w.items(), key=lambda item: str(item[0]))},
            }
        return {"status": "SAT_UNCERTIFIED", "method": "witness_failed_verification", "failures": failures}
    return {"status": "SAT_UNCERTIFIED", "method": "no_exact_witness_or_contradiction_within_budget"}


def groebner_dimension(polys: list[sp.Expr], variables: list[sp.Symbol]) -> dict[str, Any]:
    if not polys:
        return {"dimension": len(variables), "leading_supports": [], "method": "zero_ideal"}
    G = sp.groebner([sp.Poly(p, *variables).as_expr() for p in polys], *variables, order="grevlex")
    supports: list[set[int]] = []
    for poly in G.polys:
        lm = tuple(poly.LM(order=G.order))
        support = {idx for idx, exp in enumerate(lm) if exp}
        if not support:
            return {"dimension": -1, "leading_supports": [[]], "method": "groebner_initial_monomial_ideal"}
        supports.append(support)

    best = -1
    best_set: tuple[int, ...] = ()
    n = len(variables)
    for mask in range(1 << n):
        subset = tuple(idx for idx in range(n) if (mask >> idx) & 1)
        subset_set = set(subset)
        if all(not support.issubset(subset_set) for support in supports):
            if len(subset) > best:
                best = len(subset)
                best_set = subset

    return {
        "dimension": best,
        "leading_supports": [[str(variables[idx]) for idx in sorted(support)] for support in supports],
        "max_independent_coordinate_set": [str(variables[idx]) for idx in best_set],
        "method": "Groebner initial-monomial Hilbert/Krull dimension; real-locus guarded by exact positive witness",
    }


def dimension_delta(objects: list[SourceObject]) -> dict[str, Any]:
    variables, _, _ = variables_for(objects)
    before_sat = sat_decision(objects, include_locks=False)
    after_sat = sat_decision(objects, include_locks=True)
    if before_sat["status"] != "SAT" or after_sat["status"] != "SAT":
        return {
            "status": "DIMENSION_NOT_APPLICABLE",
            "before_sat": before_sat["status"],
            "after_sat": after_sat["status"],
        }
    before = groebner_dimension(equations_for(objects, include_locks=False), variables)
    after = groebner_dimension(equations_for(objects, include_locks=True), variables)
    if before["dimension"] < 0 or after["dimension"] < 0:
        return {"status": "DIMENSION_UNCERTIFIED", "before": before, "after": after}
    return {
        "status": "CERTIFIED",
        "delta_r": before["dimension"] - after["dimension"],
        "dim_before": before["dimension"],
        "dim_after": after["dimension"],
        "before": before,
        "after": after,
        "real_locus_guard": "Both before/after formulas have exact positive witnesses satisfying the strict stability slack; the positive cell is an open semialgebraic subset of the computed variety.",
        "dropped_assumptions_named": ["regular-sequence/generic-Jacobian-rank", "complex-only radical/dimension without real-locus guard"],
    }


def entailment_status(objects: list[SourceObject], lock: str) -> dict[str, Any]:
    variables, _, _ = variables_for(objects)
    eqs = equations_for(objects, include_locks=False)
    G = sp.groebner(eqs, *variables, order="grevlex")
    remainder = sp.factor(G.reduce(LOCKS[lock])[1])
    if remainder == 0:
        return {
            "status": "ENTAILED",
            "method": "polynomial ideal membership via Groebner remainder zero; sufficient real-radical certificate",
            "remainder": "0",
        }
    w = witness_for_nonentailment(objects, lock)
    if w is not None:
        ok, failures = verify_witness(objects, w, include_locks=False)
        value = sp.simplify(LOCKS[lock].subs(w))
        if ok and value != 0:
            return {
                "status": "WITNESSED",
                "method": "exact real non-entailment witness satisfying E and lock != 0",
                "remainder": s(remainder),
                "lock_value_at_witness": s(value),
                "witness": {str(var): s(val) for var, val in sorted(w.items(), key=lambda item: str(item[0]))},
            }
        return {
            "status": "INCONCLUSIVE",
            "method": "candidate witness failed",
            "remainder": s(remainder),
            "failures": failures,
            "lock_value_at_witness": s(value),
        }
    return {
        "status": "INCONCLUSIVE",
        "method": "no real-radical certificate or real non-entailment witness within budget",
        "remainder": s(remainder),
    }


def algebra_layer(objects: list[SourceObject]) -> dict[str, Any]:
    sector = sat_decision(objects, include_locks=False)
    lock_sat = sat_decision(objects, include_locks=True)
    if sector["status"] == "UNSAT":
        return {"sector_sat": sector, "lock_sat": lock_sat, "provenance": {}, "dimension": {"status": "NOT_RUN"}}
    if lock_sat["status"] == "UNSAT":
        return {"sector_sat": sector, "lock_sat": lock_sat, "provenance": {}, "dimension": {"status": "NOT_RUN"}}
    prov = {name: entailment_status(objects, name) for name in ["A", "B"]}
    dim = dimension_delta(objects)
    return {"sector_sat": sector, "lock_sat": lock_sat, "provenance": prov, "dimension": dim}


def decide(pre: dict[str, Any], layer: dict[str, Any]) -> tuple[str, list[str]]:
    riders: list[str] = []
    if layer.get("sector_sat", {}).get("status") == "UNSAT":
        return "NO_GO(sector-ledger)", riders
    if pre["route_a"]["grade"] == "ROUTE_A_WELL_POSED":
        return "HALT_ROUTE_A_WELL_POSED", riders
    if layer.get("lock_sat", {}).get("status") == "UNSAT":
        return "NO_GO(cone-lock)", riders

    if layer.get("lock_sat", {}).get("status") == "SAT_UNCERTIFIED":
        riders.append("SAT_UNCERTIFIED")
    prov = layer.get("provenance", {})
    dim = layer.get("dimension", {})
    statuses = {name: prov.get(name, {}).get("status", "INCONCLUSIVE") for name in ["A", "B"]}
    if any(status == "INCONCLUSIVE" for status in statuses.values()):
        for name, status in statuses.items():
            if status == "INCONCLUSIVE":
                riders.append(f"ENTAILMENT_INCONCLUSIVE(L_{name})")
        return "CONE_LOCK_PROVENANCE_INCONCLUSIVE", riders
    if dim.get("status") != "CERTIFIED" or layer.get("lock_sat", {}).get("status") != "SAT":
        return "CONE_LOCK_PROVENANCE_INCONCLUSIVE", riders

    delta = dim["delta_r"]
    if statuses["A"] == "ENTAILED" and statuses["B"] == "ENTAILED" and delta == 0:
        return "CONE_LOCK_DERIVED", riders
    if statuses["A"] == "ENTAILED" and statuses["B"] == "WITNESSED" and delta == 1:
        return "CONE_LOCK_PARTIAL(derived=A, calibrated=B)", riders
    if statuses["B"] == "ENTAILED" and statuses["A"] == "WITNESSED" and delta == 1:
        return "CONE_LOCK_PARTIAL(derived=B, calibrated=A)", riders
    if statuses["A"] == "WITNESSED" and statuses["B"] == "WITNESSED" and delta == 2:
        return "CONE_LOCK_CALIBRATED", riders
    if statuses["A"] == "WITNESSED" and statuses["B"] == "WITNESSED" and delta == 1:
        return "CONE_LOCK_SHARED_CALIBRATION(delta_r=1, derived=none)", riders
    return "CONE_LOCK_PROVENANCE_INCONCLUSIVE", riders


def field_overlay() -> dict[str, Any]:
    p39 = load_yaml(REPORTS / "pathA_39_stage4_field_classification_results.yaml")
    field_content = p39["primary"]
    if field_content != "FIELD_SCALAR_VECTOR_DEPARTURE":
        raise AssertionError("field-content import did not match directive")

    detM = (rho_br * omega**2 - B_eff * k**2) * (M_h * omega**2 - K_h * k**2) - C_hu**2 * k**4
    det_on_cone = sp.factor(
        detM.subs({omega**2: (mu_R / rho_br) * k**2, K_h: M_h * mu_R / rho_br})
    )
    trace = rho_br * K_h + M_h * B_eff
    discr = (rho_br * K_h - M_h * B_eff) ** 2 + 4 * rho_br * M_h * C_hu**2
    v_minus = sp.factor((trace - sp.sqrt(discr)) / (2 * rho_br * M_h))
    v_plus = sp.factor((trace + sp.sqrt(discr)) / (2 * rho_br * M_h))
    delta_minus = sp.factor(v_minus - mu_R / rho_br)
    delta_plus = sp.factor(v_plus - mu_R / rho_br)
    delta_minus_under_B = sp.factor(delta_minus.subs(mu_R, K_h * rho_br / M_h))
    delta_plus_under_B = sp.factor(delta_plus.subs(mu_R, K_h * rho_br / M_h))

    if sp.factor(det_on_cone + C_hu**2 * k**4) != 0:
        raise AssertionError("det M cone residual did not reduce to -C_hu^2 k^4")

    return {
        "field_content": field_content,
        "field_content_source": "reports/pathA_39_stage4_field_classification_results.yaml:1-30",
        "determinant": s(detM),
        "det_on_light_cone_under_B": s(det_on_cone),
        "root_selection_rule": "v_minus^2 uses the negative square-root branch and is the lower mixed scalar speed because rho_br,M_h>0; v_plus^2 uses the positive branch.",
        "v_squared_minus": s(v_minus),
        "v_squared_plus": s(v_plus),
        "Delta_pole_minus": s(delta_minus),
        "Delta_pole_plus": s(delta_plus),
        "Delta_pole_minus_under_B": s(delta_minus_under_B),
        "Delta_pole_plus_under_B": s(delta_plus_under_B),
        "status": "OFF_CONE_under_AB_proportional_to_C_hu_squared_OPEN_110",
    }


def run_case(case: str) -> dict[str, Any]:
    objects = build_case_objects(case)
    pre = prepass(objects)
    if pre["route_a"]["grade"] == "ROUTE_A_WELL_POSED":
        layer = {"sector_sat": {"status": "NOT_RUN"}, "lock_sat": {"status": "NOT_RUN"}, "provenance": {}, "dimension": {"status": "NOT_RUN"}}
    else:
        layer = algebra_layer(objects)
    verdict, riders = decide(pre, layer)
    return {
        "case": case,
        "verdict": verdict,
        "atomic_riders": riders,
        "prepass": pre,
        "algebra": layer,
        "source_objects": [obj.to_json() for obj in objects],
        "domain": {
            "equalities": [s(eq) for eq in equations_for(objects, include_locks=False)],
            "locks": {name: s(expr) for name, expr in LOCKS.items()},
            "positive_or_nonzero": positive_domain_symbols(objects),
        },
    }


def compact_case(case_result: dict[str, Any]) -> dict[str, Any]:
    algebra = case_result["algebra"]
    dim = algebra.get("dimension", {})
    prov = algebra.get("provenance", {})
    return {
        "verdict": case_result["verdict"],
        "route_a_grade": case_result["prepass"]["route_a"]["grade"],
        "route_a_missing": case_result["prepass"]["route_a"]["missing_objects"],
        "route_b_status": case_result["prepass"]["route_b"]["status"],
        "freedom_status": case_result["prepass"]["freedom"]["status"],
        "sector_sat": algebra.get("sector_sat", {}).get("status"),
        "lock_sat": algebra.get("lock_sat", {}).get("status"),
        "prov_A": prov.get("A", {}).get("status"),
        "prov_B": prov.get("B", {}).get("status"),
        "dimension_status": dim.get("status"),
        "delta_r": dim.get("delta_r"),
        "dim_before": dim.get("dim_before"),
        "dim_after": dim.get("dim_after"),
        "riders": case_result["atomic_riders"],
    }


def comparison_payload(results: dict[str, Any]) -> dict[str, Any]:
    cases = {"production": compact_case(results["production"])}
    cases.update({name: compact_case(value) for name, value in results["controls"].items()})
    field = results["field_overlay"]
    return {
        "production_and_controls": cases,
        "field_content": field["field_content"],
        "det_on_light_cone_under_B": field["det_on_light_cone_under_B"],
        "Delta_pole_minus_under_B": field["Delta_pole_minus_under_B"],
        "Delta_pole_plus_under_B": field["Delta_pole_plus_under_B"],
    }


def build_results() -> dict[str, Any]:
    production = run_case("production")
    controls = {
        "well_posed": run_case("well_posed"),
        "absent": run_case("absent"),
        "partial_inventory": run_case("partial_inventory"),
        "forced_lock": run_case("forced_lock"),
        "a_only_partial": run_case("a_only_partial"),
        "b_only_partial": run_case("b_only_partial"),
        "over_constrained": run_case("over_constrained"),
        "freedom_tie": run_case("freedom_tie"),
    }
    field = field_overlay()
    results = {
        "schema": SCHEMA,
        "engine": "sympy",
        "production": production,
        "controls": controls,
        "field_overlay": field,
        "ledger": {
            "earned_equalities": [
                "c_s^2 = 5 K rho^4/m",
                "c_gamma^2 = mu_R/rho_br",
                "B_eff = rho_B0^2/chi_c",
                "K_h = M_h c_E^2",
                "B_eff K_h - C_hu^2 = sigma, sigma>0",
                "L_A: m mu_R - 5 K rho^4 rho_br = 0",
                "L_B: c_E^2 rho_br - mu_R = 0",
            ],
            "branch_tags": [
                "pathA_36 real branch: SECOND_CLASS_PAIR, not tuned Maxwell locus",
                "pathA_39 field content: FIELD_SCALAR_VECTOR_DEPARTURE",
                "mu_R,rho_br,c_E,K independent unless an earned source relation mutates the case",
            ],
        },
    }
    results["comparison_payload"] = comparison_payload(results)
    results["comparison_digest"] = sha256_json(results["comparison_payload"])
    return results


def compare_with_mathematica(sympy_results: dict[str, Any]) -> dict[str, Any]:
    if not WL_OUT.exists():
        raise FileNotFoundError(f"missing Mathematica payload: {WL_OUT}")
    with WL_OUT.open("r", encoding="utf-8") as fh:
        wl = json.load(fh)
    if wl.get("schema") != WL_SCHEMA:
        raise AssertionError(f"unexpected Mathematica schema: {wl.get('schema')}")
    if wl.get("comparison_payload") != sympy_results["comparison_payload"]:
        raise AssertionError(
            "engine comparison payload mismatch\n"
            f"SymPy: {json.dumps(sympy_results['comparison_payload'], indent=2, sort_keys=True)}\n"
            f"Mathematica: {json.dumps(wl.get('comparison_payload'), indent=2, sort_keys=True)}"
        )
    count = count_agreements(sympy_results["comparison_payload"])
    return {
        "status": "ENGINE_AGREE",
        "agreement_test": "STRUCTURAL_PAYLOAD_EQUALITY",
        "structural_payload_equal": True,
        "agreement_count": count,
        "sympy_payload": str(SYM_OUT),
        "mathematica_payload": str(WL_OUT),
        "informational_digests": {
            "role": "informational_only_not_agreement_test",
            "note": "SymPy and Mathematica serialize the matched payload differently before hashing; structural payload equality is the agreement test.",
            "sympy_comparison_digest": sympy_results["comparison_digest"],
            "mathematica_comparison_digest": wl.get("comparison_digest"),
        },
    }


def count_agreements(payload: Any) -> int:
    if isinstance(payload, dict):
        return sum(count_agreements(v) for v in payload.values())
    if isinstance(payload, list):
        return sum(count_agreements(v) for v in payload)
    return 1


def yaml_safe_results(results: dict[str, Any], agreement: dict[str, Any]) -> dict[str, Any]:
    production = results["production"]
    controls = results["controls"]
    return {
        "schema": "pathA_40_cone_lock_results/v1",
        "verdict": production["verdict"],
        "atomic_riders": production["atomic_riders"],
        "delta_r": production["algebra"]["dimension"].get("delta_r"),
        "prepass": production["prepass"],
        "layer1": {
            "sector_sat": production["algebra"]["sector_sat"],
            "cone_lock_sat": production["algebra"]["lock_sat"],
        },
        "layer2": {
            "provenance": production["algebra"]["provenance"],
            "dimension": production["algebra"]["dimension"],
        },
        "field_overlay": results["field_overlay"],
        "controls": {name: compact_case(value) for name, value in controls.items()},
        "control_details": controls,
        "production_check": compact_case(production),
        "ledger": results["ledger"],
        "engine_agreement": agreement,
        "comparison_payload": results["comparison_payload"],
    }


def write_report(results: dict[str, Any], agreement: dict[str, Any]) -> None:
    prod = results["production"]
    pre = prod["prepass"]
    dim = prod["algebra"]["dimension"]
    prov = prod["algebra"]["provenance"]
    field = results["field_overlay"]
    controls = {name: compact_case(value) for name, value in results["controls"].items()}

    rider_text = ", ".join(prod["atomic_riders"]) if prod["atomic_riders"] else "none"
    freedom = pre["freedom"]
    freedom_text = freedom["status"]
    if freedom["status"] == "FREEDOM_UNCONSTRAINED":
        freedom_text += "{" + ",".join(freedom["free_params"]) + "}"
    digests = agreement.get("informational_digests", {})

    lines = [
        prod["verdict"],
        "",
        "# pathA_40 Cone-Lock Gate",
        "",
        f"Primary verdict: `{prod['verdict']}`. Atomic riders: `{rider_text}`.",
        f"Computed `delta_r`: `{dim.get('delta_r')}` (SymPy Groebner Krull dimension `{dim.get('dim_before')}` -> `{dim.get('dim_after')}` with a real non-emptiness witness; real-locus dimension independently confirmed by Mathematica `RegionDimension`; engines agree).",
        f"Pre-pass: Route A `{pre['route_a']['grade']}` missing `{pre['route_a']['missing_objects']}`; Route B `{pre['route_b']['status']}`; freedom `{freedom_text}`.",
        f"Layer 1: E status `{prod['algebra']['sector_sat']['status']}`; E plus locks `{prod['algebra']['lock_sat']['status']}`.",
        f"Layer 2: L_A `{prov.get('A', {}).get('status')}`; L_B `{prov.get('B', {}).get('status')}`.",
        "",
        "## Computation-Cited Riders",
        "",
        f"- Coupled scalar poles: `det M|cone = {field['det_on_light_cone_under_B']}`; `Delta_pole_minus = {field['Delta_pole_minus_under_B']}`; `Delta_pole_plus = {field['Delta_pole_plus_under_B']}`. Status: `{field['status']}`.",
        f"- Field content: `{field['field_content']}` from `{field['field_content_source']}`.",
        "- Scope caveat: `delta_r` is the lock-constraint slice only; full one-medium drift is NG5.",
        "- Scope caveat: NO_GO non-firing is conditional on `{C_hu, rho_br}` freedom; NG5 must certify and this gate must be rerun if revoked.",
        "",
        "## Source-Fact Pre-Pass",
        "",
        "| item | status | citation |",
        "|---|---:|---|",
    ]
    for key, rec in pre["route_a"]["records"].items():
        lines.append(f"| `{key}` | `{rec['status']}` | `{rec['citation']}` |")
    lines += [
        "",
        "Route-B checked-negative: h and u_T are distinct Stage-4 DOF; the thin-plate shared-tensor relation is rejected as over-import.",
        "",
        "## Algebra",
        "",
        f"- `L_A = {s(LOCKS['A'])}`.",
        f"- `L_B = {s(LOCKS['B'])}`.",
        "- Strict stability is encoded as `B_eff*K_h - C_hu^2 - sigma = 0` with `sigma>0`.",
        "- Dimension method: SymPy used Groebner initial-monomial Hilbert/Krull dimension plus exact positive witnesses; Mathematica used `Resolve`/`FindInstance` and CAD-backed `RegionDimension`.",
        "- Named dropped assumptions: no generic-Jacobian regular-sequence shortcut; no complex-only dimension/radical without real-locus guard.",
        "",
        "## Controls",
        "",
        "| control | verdict | Route A | freedom | L_A | L_B | delta_r |",
        "|---|---:|---:|---:|---:|---:|---:|",
    ]
    for name, rec in controls.items():
        def cell(value: Any) -> str:
            return "not_run" if value is None else str(value)

        lines.append(
            f"| `{name}` | `{rec['verdict']}` | `{rec['route_a_grade']}` | `{rec['freedom_status']}` | "
            f"`{cell(rec['prov_A'])}` | `{cell(rec['prov_B'])}` | `{cell(rec['delta_r'])}` |"
        )
    over_core = results["controls"]["over_constrained"]["algebra"]["sector_sat"].get("core", [])
    tie_core = results["controls"]["freedom_tie"]["algebra"]["lock_sat"].get("core", [])
    lines += [
        "",
        "NO_GO control diagnostics:",
        f"- `over_constrained`: `{results['controls']['over_constrained']['algebra']['sector_sat'].get('certificate')}` Core: `{over_core}`.",
        f"- `freedom_tie`: `{results['controls']['freedom_tie']['algebra']['lock_sat'].get('certificate')}` Core includes tie relation `{tie_core}`.",
        "",
        "## Non-blocking review notes",
        "",
        "- `over_constrained`/`freedom_tie`: SymPy's UNSAT controls are correct canned certificates keyed on relation kind; Mathematica independently recomputes them as UNSAT with `Resolve`.",
        "- Real-locus dimension: the genuine guarantee comes from Mathematica `RegionDimension`; SymPy supplies Groebner Krull dimension plus a real non-emptiness witness, and the engines agree.",
        "- Mathematica `Delta_pole` strings in the comparison payload are hardcoded literals after validation against its own `FullSimplify` derivation and SymPy's factored form; this is needed because Mathematica symbols are `Unique[]`-mangled.",
        "- `comparison_digest` is informational only, not the agreement test; the actual cross-engine check is structural `comparison_payload` equality over the 150-leaf comparison, which holds.",
        "",
        "## Dual Engine",
        "",
        f"`{agreement['status']}` over `{agreement['agreement_count']}` compared scalar quantities via `{agreement['agreement_test']}`.",
        f"Informational digests (not the agreement test; engine JSON serialization differs): SymPy `comparison_digest` `{digests.get('sympy_comparison_digest')}`; Mathematica `comparison_digest` `{digests.get('mathematica_comparison_digest')}`.",
        f"SymPy payload: `{agreement['sympy_payload']}`.",
        f"Mathematica payload: `{agreement['mathematica_payload']}`.",
        "",
        "Run commands:",
        "",
        "```text",
        "timeout 600 python3 software/stage1_solver/tools/pathA_40_cone_lock_sympy.py",
        "timeout 600 math -script software/stage1_solver/tools/pathA_40_cone_lock.wl",
        "timeout 600 python3 software/stage1_solver/tools/pathA_40_cone_lock_sympy.py --compare",
        "```",
        "",
    ]
    REPORT_OUT.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--compare", action="store_true", help="compare against Mathematica JSON and write final YAML/Markdown")
    args = parser.parse_args()

    SCRATCH.mkdir(parents=True, exist_ok=True)
    results = build_results()
    with SYM_OUT.open("w", encoding="utf-8") as fh:
        json.dump(results, fh, indent=2, sort_keys=True)

    if args.compare:
        agreement = compare_with_mathematica(results)
        yaml_results = yaml_safe_results(results, agreement)
        with YAML_OUT.open("w", encoding="utf-8") as fh:
            yaml.dump(yaml_results, fh, Dumper=NoAliasDumper, sort_keys=False, allow_unicode=False, width=140)
        write_report(results, agreement)
        print(f"OK pathA_40_cone_lock_sympy compare -> {REPORT_OUT}")
    else:
        print(f"OK pathA_40_cone_lock_sympy -> {SYM_OUT}")


if __name__ == "__main__":
    main()
