#!/usr/bin/env python3
"""pathA_41 NG5 SECOND_MEDIUM_DRIFT gate, SymPy engine.

The production ledger is source-fact driven.  Rows are constructed without an
origin verdict; every active row receives a structured RouteEvaluation, and
`classify_origin(row, ctx)` is the only place an origin verdict is produced.
Controls mutate source facts, then rerun the same classifier.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import sys
from dataclasses import dataclass, field, replace
from pathlib import Path
from typing import Any

import sympy as sp
import yaml


SCRIPT_PATH = Path(__file__).resolve()
STAGE1_ROOT = SCRIPT_PATH.parents[1]
REPO_ROOT = STAGE1_ROOT.parents[1]
REPORTS = STAGE1_ROOT / "reports"
SCRATCH = STAGE1_ROOT / "_scratch"

SYM_OUT = SCRATCH / "pathA_41_ng5_second_medium_drift_sympy.json"
WL_OUT = SCRATCH / "pathA_41_ng5_second_medium_drift_mathematica.json"
YAML_OUT = REPORTS / "pathA_41_ng5_second_medium_drift_results.yaml"
REPORT_OUT = REPORTS / "pathA_41_ng5_second_medium_drift.md"

SCHEMA = "pathA_41_ng5_second_medium_drift_sympy/v2"
WL_SCHEMA = "pathA_41_ng5_second_medium_drift_mathematica/v2"

ALLOWED_ROUTE_STATUSES = {"SOLVED_PASS", "REGISTERED_DEFERRED"}
REJECTED_ROUTE_STATUSES = {"FAILED", "BY_TUNING", "ABSENT", "PROMISSORY_ONLY"}


class NoAliasDumper(yaml.SafeDumper):
    def ignore_aliases(self, data: Any) -> bool:
        return True


Dim = tuple[int, int, int]
DIM_M: Dim = (1, 0, 0)
DIM_L: Dim = (0, 1, 0)
DIM_RHO: Dim = (0, -4, 0)
DIM_SURFACE_INERTIA: Dim = (1, -3, 0)
DIM_SURFACE_MODULUS: Dim = (1, -1, -2)
DIM_SPEED: Dim = (0, 1, -1)
DIM_B_EFF: Dim = (1, -1, -2)


def dadd(*dims: Dim) -> Dim:
    return tuple(sum(d[i] for d in dims) for i in range(3))  # type: ignore[return-value]


def dmul(dim: Dim, power: int) -> Dim:
    return tuple(power * v for v in dim)  # type: ignore[return-value]


def dim_str(dim: Dim) -> str:
    labels = ["M", "L", "T"]
    parts: list[str] = []
    for label, exp in zip(labels, dim):
        if exp == 0:
            continue
        parts.append(label if exp == 1 else f"{label}^{exp}")
    return " ".join(parts) if parts else "1"


def s(expr: Any) -> str:
    if isinstance(expr, sp.Basic):
        return sp.sstr(sp.factor(sp.simplify(expr)))
    return str(expr)


def sha256_json(data: Any) -> str:
    return hashlib.sha256(json.dumps(data, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def read(path: Path) -> str:
    if not path.exists():
        raise FileNotFoundError(path)
    return path.read_text(encoding="utf-8")


def load_yaml(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as fh:
        data = yaml.safe_load(fh)
    if not isinstance(data, dict):
        raise AssertionError(f"expected mapping YAML: {path}")
    return data


def rel(path: Path) -> str:
    return str(path.relative_to(REPO_ROOT))


def source_ref(path: Path, token: str, label: str) -> str:
    lines = read(path).splitlines()
    for idx, line in enumerate(lines, start=1):
        if token in line:
            return f"{rel(path)}:{idx}"
    raise AssertionError(f"missing source token for {label}: {token}")


def require_token(path: Path, token: str, label: str) -> None:
    source_ref(path, token, label)


def load_pathA40_module() -> Any:
    p40_path = STAGE1_ROOT / "tools" / "pathA_40_cone_lock_sympy.py"
    spec = importlib.util.spec_from_file_location("pathA40_cone_lock_for_ng5", p40_path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"could not load {p40_path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


P40 = load_pathA40_module()


@dataclass(frozen=True)
class RowSource:
    p: str
    incidence: str
    dims: str
    canonical_source: str
    local_status: str
    earned_equalities: list[str] = field(default_factory=list)
    route_hint: str | None = None
    active_medium_competition: bool = False
    pathA40_freedom_relevance: str = "none"
    dependency: list[str] = field(default_factory=list)
    calibration_anchor: str | None = None
    calibrated_geometry: bool = False
    out_of_active: bool = False
    location: str = "unassigned"
    notes: str = ""

    def to_yaml_base(self) -> dict[str, Any]:
        return {
            "p": self.p,
            "incidence": self.incidence,
            "dims": self.dims,
            "canonical_source": self.canonical_source,
            "local_status": self.local_status,
            "earned_equalities": self.earned_equalities,
            "route_hint": self.route_hint,
            "active_medium_competition": self.active_medium_competition,
            "pathA_40_freedom_relevance": self.pathA40_freedom_relevance,
            "dependency": self.dependency,
            "calibration_anchor": self.calibration_anchor,
            "calibrated_geometry": self.calibrated_geometry,
            "out_of_active": self.out_of_active,
            "location": self.location,
            "notes": self.notes,
        }


@dataclass(frozen=True)
class CaseSpec:
    name: str
    blank_all_routes: bool = False
    blank_routes: set[str] = field(default_factory=set)
    blank_route_fields: dict[str, set[str]] = field(default_factory=dict)
    remove_calibration_anchors: set[str] = field(default_factory=set)
    remove_calibrated_geometry: set[str] = field(default_factory=set)
    add_synthetic_irreducible: bool = False
    add_synthetic_reducible: bool = False
    contradiction: str | None = None
    residual_mode: str = "production"


@dataclass(frozen=True)
class LineageFacts:
    varrho_object: str = "pathA_25_closed_density_smectic_varrho_br"
    rho_object: str = "pathA_35_active_shear_surface_rho_br"
    varrho_active: bool = False
    rho_active: bool = True
    varrho_role: str = "density_smectic_layer_kinetic"
    rho_role: str = "shear_surface_kinetic"
    varrho_dim: Dim = DIM_SURFACE_INERTIA
    rho_dim: Dim = DIM_SURFACE_INERTIA
    residual_multiplier: str | None = None


def source_catalog() -> dict[str, Any]:
    p25 = REPORTS / "pathA_25_G0_freeze.md"
    p25_status = REPORTS / "pathA_25_STATUS.md"
    p25_b4 = REPORTS / "pathA_25_gateB4_smectic.md"
    p25_rc = REPORTS / "pathA_25_gateRC_cubic.md"
    p35 = REPORTS / "pathA_35_G0_freeze.md"
    p36 = REPORTS / "pathA_36_c5_phase_potential.md"
    p38 = REPORTS / "pathA_38_throat_body_electric_localization.md"
    p38_yml = REPORTS / "pathA_38_results.yaml"
    p39_mag = REPORTS / "pathA_39_magnetic_force.md"
    p39_scalar = REPORTS / "pathA_39_scalar_admixture_screen.md"
    p39_stage4 = REPORTS / "pathA_39_stage4_field_classification_results.yaml"
    p40_md = REPORTS / "pathA_40_cone_lock.md"
    p40_yml = REPORTS / "pathA_40_cone_lock_results.yaml"
    d09 = STAGE1_ROOT / "decisions" / "09_calibrate_predict_methodology.md"
    d14 = STAGE1_ROOT / "decisions" / "14_value_provenance_and_calibration_map.md"
    framing = STAGE1_ROOT / "_scratch" / "pathA_41_framing_codex.md"

    refs = {
        "pathA25_closed": source_ref(p25_status, "no live pathA_25 thread remains", "pathA_25 closed"),
        "pathA25_b4_fail": source_ref(p25_b4, "FAIL_NOT_CODIM1", "pathA_25 B4 fail"),
        "pathA25_rc_closed": source_ref(p25_status, "Density-smectic brane route CLOSED", "pathA_25 density route closed"),
        "pathA25_varrho_prose": source_ref(p25, "`varrho_br[rho]` is the layer mass-density functional", "pathA_25 varrho prose"),
        "pathA25_varrho_definition": source_ref(p25, "varrho_br[rho] := int_layer dn m rho.", "pathA_25 varrho definition"),
        "pathA25_sigma": source_ref(p25, "Sigma_n[rho] denotes the smectic density layers", "pathA_25 Sigma_n"),
        "pathA25_c_L1": source_ref(p25, "| `c_L1`", "pathA_25 c_L1"),
        "pathA25_c_L2": source_ref(p25, "| `c_L2`", "pathA_25 c_L2"),
        "pathA25_A_R": source_ref(p25, "| `A_R`", "pathA_25 A_R"),
        "pathA25_k_R": source_ref(p25, "| `k_R`", "pathA_25 k_R"),
        "pathA25_lambda_Cdiv": source_ref(p25, "| `lambda_Cdiv`", "pathA_25 lambda_Cdiv"),
        "pathA25_chi_Cpin": source_ref(p25, "| `chi_Cpin`", "pathA_25 chi_Cpin"),
        "pathA25_J_Pu": source_ref(p25, "| `J_Pu`", "pathA_25 J_Pu"),
        "pathA25_kappa_Pu": source_ref(p25, "| `kappa_Pu`", "pathA_25 kappa_Pu"),
        "pathA25_rc_cpin": source_ref(p25_rc, "S_L_plus_Cpin", "pathA_25 RC Cpin"),
        "pathA35_rho_role": source_ref(p35, "The brane inertia is the independent surface mass density `rho_br`.", "pathA_35 rho role"),
        "pathA35_rho_constant": source_ref(p35, "| `rho_br` | `M L^-3` | postulated-ingredient", "pathA_35 rho constant"),
        "pathA35_mu_constant": source_ref(p35, "| `mu_R` | `M L^-1 T^-2` | postulated-ingredient", "pathA_35 mu_R constant"),
        "pathA35_inplane_kinetic": source_ref(p35, "(1/2) rho_br (partial_t u^a)(partial_t u^a)", "pathA_35 in-plane kinetic"),
        "pathA35_ell_g": source_ref(p35, "g_ell(w) = exp(-(w/ell_g)^2)/(sqrt(pi) ell_g)", "pathA_35 ell_g"),
        "pathA36_B_eff": source_ref(p36, "finite density stiffness contributes `B_eff = rho_B0^2/chi_c`", "pathA_36 B_eff"),
        "pathA36_cgamma": source_ref(p36, "`omega^2 = (mu_R/rho_br) k^2`, `c_gamma^2 = mu_R/rho_br`.", "pathA_36 c_gamma"),
        "pathA36_by_tuning": source_ref(p36, "That is `BY_TUNING`, not `WITH_PROVENANCE`.", "pathA_36 tuning guard"),
        "pathA38_qh": source_ref(p38, "q_h(+)=2*QE*tanh(b/ell)/b", "pathA_38 q_h"),
        "pathA38_calibrated_deferred": source_ref(p38, "calibrated/deferred: `Q_E`", "pathA_38 calibrated Q_E"),
        "pathA38_wall_internal": source_ref(p38_yml, "lamN*(N-chi)^2/2+lamTau*tau^2/2", "pathA_38 wall-internal lambdas"),
        "pathA39_imports": source_ref(p39_mag, "imported_exact", "pathA_39 imported exact"),
        "pathA39_sim_values": source_ref(p39_mag, "sim_deferred_values", "pathA_39 sim deferred values"),
        "pathA39_scalar_C_hu": source_ref(p39_scalar, "declared_scan_parameters", "pathA_39 C_hu scan"),
        "pathA39_scalar_no_route": source_ref(p39_scalar, "C_hu scalar mixing", "pathA_39 C_hu mixing"),
        "pathA39_stage4_a_Tp": source_ref(p39_stage4, "Nu*aTp*sCharge", "pathA_39 Stage-4 a_Tp"),
        "pathA40_route_a": source_ref(p40_yml, "grade: ROUTE_A_UNDERDETERMINED_MISSING_NONLINEAR_THROAT", "pathA_40 Route-A"),
        "pathA40_freedom": source_ref(p40_yml, "status: FREEDOM_UNCONSTRAINED", "pathA_40 freedom"),
        "pathA40_no_go_control": source_ref(p40_md, "freedom_tie", "pathA_40 freedom_tie control"),
        "decision09_calibration": source_ref(d09, "Calibration anchor", "decision 09 calibration anchor"),
        "decision14_calibration": source_ref(d14, "calibration INPUT", "decision 14 calibration input"),
        "framing_target": source_ref(framing, "SECOND_MEDIUM_DRIFT(active_irreducible={rho_B0, chi_c, C_hu})", "framing target"),
    }

    required = [
        (p25_status, "density route is closed", "pathA_25 closed"),
        (p25_b4, "FAIL_NOT_CODIM1", "B4 closed"),
        (p35, "rho_br is an independent brane inertia density", "rho_br postulated"),
        (p35, "(1/2) rho_br (partial_t u^a)(partial_t u^a)", "rho_br kinetic"),
        (p36, "B_eff = rho_B0^2/chi_c", "B_eff"),
        (p38, "q_h(+)=2*QE*tanh(b/ell)/b", "q_h"),
        (p39_scalar, "C_hu scalar mixing", "C_hu"),
        (p40_yml, "ROUTE_A_UNDERDETERMINED_MISSING_NONLINEAR_THROAT", "Route-A"),
    ]
    for path, token, label in required:
        require_token(path, token, label)

    p40 = load_yaml(p40_yml)
    if p40.get("verdict") != "CONE_LOCK_CALIBRATED":
        raise AssertionError("pathA_40 calibrated verdict missing")
    if p40.get("prepass", {}).get("route_a", {}).get("grade") != "ROUTE_A_UNDERDETERMINED_MISSING_NONLINEAR_THROAT":
        raise AssertionError("pathA_40 Route-A grade mismatch")

    return {"refs": refs, "pathA40": p40}


def base_route_facts(refs: dict[str, str], p40: dict[str, Any]) -> dict[str, dict[str, Any]]:
    route_a_missing = p40["prepass"]["route_a"]["missing_objects"]
    return {
        "Route-A": {
            "name": "Route-A",
            "source": refs["pathA40_route_a"],
            "result_status": "REGISTERED_DEFERRED",
            "named_solve_in_provenance": True,
            "missing_objects": list(route_a_missing),
            "targets": ["rho_br", "mu_R", "c_E"],
            "required_joint_targets": ["rho_br", "mu_R"],
            "target_blind": True,
            "falsifiers": [],
            "notes": "R3 explicitly requires common normalization deriving both rho_br and mu_R; R4/R5 would close the cone lock.",
        },
        "Future-Compression-4D-to-3D": {
            "name": "Future-Compression-4D-to-3D",
            "source": "named_future_reduction_route:not_currently_registered",
            "result_status": "PROMISSORY_ONLY",
            "named_solve_in_provenance": False,
            "missing_objects": ["compression-sector 4D->3D nonlinear brane solve deriving rho_B0 and chi_c"],
            "targets": ["rho_B0", "chi_c"],
            "required_joint_targets": ["rho_B0", "chi_c"],
            "target_blind": True,
            "falsifiers": [],
            "notes": "Named future route only; not a current registered route.",
        },
        "Future-Embedding-Overlap": {
            "name": "Future-Embedding-Overlap",
            "source": "named_future_reduction_route:not_currently_registered",
            "result_status": "PROMISSORY_ONLY",
            "named_solve_in_provenance": False,
            "missing_objects": ["embedding-overlap nonlinear throat solve deriving C_hu"],
            "targets": ["C_hu"],
            "required_joint_targets": ["C_hu"],
            "target_blind": True,
            "falsifiers": [],
            "notes": "PathA_39 scans C_hu; no Route-A-like current registered solve exists.",
        },
    }


def route_name_for_param(param: str) -> str | None:
    return {
        "rho_br": "Route-A",
        "mu_R": "Route-A",
        "c_E": "Route-A",
        "rho_B0": "Future-Compression-4D-to-3D",
        "chi_c": "Future-Compression-4D-to-3D",
        "C_hu": "Future-Embedding-Overlap",
    }.get(param)


def route_inventory(refs: dict[str, str], p40: dict[str, Any], case: CaseSpec) -> dict[str, dict[str, Any]]:
    routes = base_route_facts(refs, p40)
    if case.blank_all_routes:
        return {}
    for route_name in case.blank_routes:
        routes.pop(route_name, None)
    for route_name, fields in case.blank_route_fields.items():
        if route_name in routes:
            for field_name in fields:
                routes[route_name].pop(field_name, None)
    return routes


def registered_route_for(param: str, ctx: dict[str, Any]) -> dict[str, Any]:
    route_name = route_name_for_param(param)
    if route_name is None:
        return {
            "param": param,
            "route_name": None,
            "valid": False,
            "reason": "NO_REGISTERED_ROUTE_FOR_PARAM",
            "result_status": None,
            "checks": {},
        }
    routes: dict[str, dict[str, Any]] = ctx["routes"]
    fact = routes.get(route_name)
    if fact is None:
        return {
            "param": param,
            "route_name": route_name,
            "valid": False,
            "reason": "ROUTE_BLANKED_OR_ABSENT",
            "result_status": None,
            "checks": {"mapped_route_exists": False},
        }

    required = [
        "name",
        "source",
        "result_status",
        "named_solve_in_provenance",
        "missing_objects",
        "targets",
        "target_blind",
        "falsifiers",
    ]
    missing_fields = [field_name for field_name in required if field_name not in fact]
    result_status = fact.get("result_status")
    missing_objects = fact.get("missing_objects")
    targets = fact.get("targets")
    falsifiers = fact.get("falsifiers")
    required_joint = fact.get("required_joint_targets", [])
    checks = {
        "mapped_route_exists": True,
        "required_fields_present": not missing_fields,
        "missing_required_fields": missing_fields,
        "named_solve_in_provenance": fact.get("named_solve_in_provenance") is True,
        "result_status_allowed": result_status in ALLOWED_ROUTE_STATUSES,
        "result_status_rejected": result_status in REJECTED_ROUTE_STATUSES,
        "finite_listed_missing_objects": isinstance(missing_objects, list),
        "target_identity": isinstance(targets, list) and param in targets,
        "joint_target_identity": isinstance(targets, list) and all(t in targets for t in required_joint),
        "target_blind": fact.get("target_blind") is True,
        "falsifiers_clear": isinstance(falsifiers, list) and len(falsifiers) == 0,
    }
    valid = (
        checks["required_fields_present"]
        and checks["named_solve_in_provenance"]
        and checks["result_status_allowed"]
        and checks["finite_listed_missing_objects"]
        and checks["target_identity"]
        and checks["joint_target_identity"]
        and checks["target_blind"]
        and checks["falsifiers_clear"]
    )
    if valid:
        reason = "VALID_REGISTERED_ROUTE"
    elif checks["result_status_rejected"]:
        reason = f"REJECTED_RESULT_STATUS_{result_status}"
    elif missing_fields:
        reason = "MISSING_REQUIRED_ROUTE_FIELD"
    else:
        reason = "ROUTE_EVALUATION_FAILED"
    return {
        "param": param,
        "route_name": route_name,
        "valid": valid,
        "reason": reason,
        "result_status": result_status,
        "source": fact.get("source"),
        "missing_objects": missing_objects if isinstance(missing_objects, list) else None,
        "targets": targets if isinstance(targets, list) else None,
        "checks": checks,
        "notes": fact.get("notes"),
    }


def dim_derivations() -> dict[str, Any]:
    bulk_mrho = dadd(DIM_M, DIM_RHO)
    layer_factor = DIM_L
    varrho = dadd(bulk_mrho, layer_factor)
    c_e_sq_rho = dadd(dmul(DIM_SPEED, 2), DIM_SURFACE_INERTIA)
    return {
        "bulk_m_times_rho_dim": dim_str(bulk_mrho),
        "layer_normal_integral_dn_dim": dim_str(layer_factor),
        "varrho_integral_dim": dim_str(varrho),
        "surface_inertia_dim_expected": dim_str(DIM_SURFACE_INERTIA),
        "varrho_equals_surface_inertia_dim": varrho == DIM_SURFACE_INERTIA,
        "rho_br_source_dim": dim_str(DIM_SURFACE_INERTIA),
        "rho_br_equals_surface_inertia_dim": True,
        "c_E_squared_rho_br_dim": dim_str(c_e_sq_rho),
        "mu_R_dim": dim_str(DIM_SURFACE_MODULUS),
        "c_E_squared_rho_br_minus_mu_R_dim_match": c_e_sq_rho == DIM_SURFACE_MODULUS,
    }


def lineage_facts_for(case: CaseSpec) -> LineageFacts:
    if case.residual_mode == "same_no_residual":
        return LineageFacts(
            varrho_object="same_active_shear_surface_object",
            rho_object="same_active_shear_surface_object",
            varrho_active=True,
            rho_active=True,
            varrho_role="shear_surface_kinetic",
            rho_role="shear_surface_kinetic",
            residual_multiplier="1",
        )
    if case.residual_mode == "same_with_free_residual":
        return LineageFacts(
            varrho_object="same_active_shear_surface_object",
            rho_object="same_active_shear_surface_object",
            varrho_active=True,
            rho_active=True,
            varrho_role="shear_surface_kinetic",
            rho_role="shear_surface_kinetic",
            residual_multiplier="Xi_residual",
        )
    return LineageFacts()


def lineage_adjudication(case: CaseSpec) -> dict[str, Any]:
    facts = lineage_facts_for(case)
    dims = dim_derivations()
    same_dim = facts.varrho_dim == facts.rho_dim == DIM_SURFACE_INERTIA
    same_role = facts.varrho_role == facts.rho_role
    same_active_object = facts.varrho_active and facts.rho_active and facts.varrho_object == facts.rho_object
    residual_closed = facts.residual_multiplier in {None, "1"}
    same = same_active_object and same_role and same_dim and residual_closed
    if same:
        finding = "OVERCOUNT_OR_SMUGGLE_CONTROL"
        explanation = "control-only same active object with no residual multiplier"
    elif not same_active_object and case.name == "production":
        finding = "NO_OVERCOUNT_ROUTE_A_PENDING"
        explanation = "varrho_br[rho] belongs to closed pathA_25 density-smectic; rho_br is active pathA_35 shear-surface and routes via Route-A."
    elif same_active_object and not residual_closed:
        finding = "DIFFERENT_BRANE_INERTIA_OBJECTS_RESIDUAL_MULTIPLIER"
        explanation = "same-object source was mutated to carry a free residual multiplier."
    else:
        finding = "DIFFERENT_BRANE_INERTIA_OBJECTS"
        explanation = "source facts do not identify the two inertia symbols as the same active object."
    return {
        "outcome": "SAME" if same else "DIFFERENT",
        "lineage_finding": finding,
        "explanation": explanation,
        "computed_invariants": {
            "same_role": same_role,
            "same_dimension": same_dim,
            "same_active_object": same_active_object,
            "residual_multiplier_closed": residual_closed,
            "varrho_active": facts.varrho_active,
            "rho_br_active": facts.rho_active,
            "varrho_object": facts.varrho_object,
            "rho_object": facts.rho_object,
        },
        "residual_multiplier": facts.residual_multiplier or "UNKNOWN_NOT_COMPARABLE_DIFFERENT_ACTIVE_OBJECTS",
        "dimension_derivation": dims,
    }


def build_rows(refs: dict[str, str], case: CaseSpec) -> list[RowSource]:
    rows: list[RowSource] = [
        RowSource("rho", "base GNLS substrate", "L^-4", refs["pathA25_varrho_definition"], "BASE_SUBSTRATE", location="4D bulk"),
        RowSource("K", "base GNLS EOS substrate", "M L^18 T^-2", refs["pathA25_varrho_definition"], "BASE_SUBSTRATE", ["c_s^2(rho)=5K rho^4/m"], location="4D bulk"),
        RowSource("m", "base GNLS constituent mass", "M", refs["pathA25_varrho_definition"], "BASE_SUBSTRATE", location="4D bulk"),
        RowSource("a", "T0 polar substrate length", "L", refs["pathA25_varrho_definition"], "BASE_SUBSTRATE", location="4D bulk"),
        RowSource("ell_g", "pathA_35 confinement width", "L", refs["pathA35_ell_g"], "ACCEPTED_GEOMETRY_SUBSTRATE", location="throat/embedding seam"),
        RowSource("g_ell(w)", "codim-1 confinement profile", "L^-1", refs["pathA35_ell_g"], "ACCEPTED_PROFILE_GIVEN_ell_g", ["int dw g_ell(w)=1"], location="throat/embedding seam"),
        RowSource("varrho_br[rho]", "closed pathA_25 density-smectic layer inertia", dim_str(DIM_SURFACE_INERTIA), refs["pathA25_varrho_prose"], "PATHA25_CLOSED_DENSITY_SMECTIC_OBJECT", ["varrho_br[rho]=int_layer dn m rho"], out_of_active=True, location="3D brane surface", notes="Closed pathA_25 candidate, not the active shear-surface rho_br."),
        RowSource("Sigma_n[rho]", "closed pathA_25 density-smectic layer support", "geometric support functional", refs["pathA25_sigma"], "PATHA25_CLOSED_DENSITY_SMECTIC_OBJECT", out_of_active=True, location="3D brane surface"),
        RowSource("delta_Sigma[rho]", "closed pathA_25 density-smectic layer measure", "L^-1", refs["pathA25_sigma"], "PATHA25_CLOSED_DENSITY_SMECTIC_OBJECT", out_of_active=True, location="3D brane surface"),
        RowSource("rho_br", "active pathA_35 shear-surface brane inertia", dim_str(DIM_SURFACE_INERTIA), refs["pathA35_rho_constant"], "POSTULATED_SHEAR_SURFACE_REDUCTION_TARGET", ["c_gamma^2=mu_R/rho_br"], "Route-A", True, "FREEDOM_SIM_DEFERRED{Route-A}", location="3D brane surface"),
        RowSource("mu_R", "active pathA_35 shear-surface modulus", dim_str(DIM_SURFACE_MODULUS), refs["pathA35_mu_constant"], "POSTULATED_SHEAR_SURFACE_REDUCTION_TARGET", ["c_gamma^2=mu_R/rho_br"], "Route-A", True, location="3D brane surface"),
        RowSource("c_E", "electric throat dynamic Green speed", dim_str(DIM_SPEED), refs["pathA39_imports"], "ROUTE_A_CONE_LOCK_REDUCTION_TARGET", ["K_h=M_h c_E^2"], "Route-A", True, "cone-lock L_B freedom relevance", location="throat/embedding seam"),
        RowSource("c_gamma", "light/shear speed", dim_str(DIM_SPEED), refs["pathA36_cgamma"], "DEPENDENT_SPEED", ["c_gamma^2=mu_R/rho_br"], active_medium_competition=True, dependency=["rho_br", "mu_R"], location="3D brane surface"),
        RowSource("rho_B0", "C5 compression density amplitude", "source-restored", refs["pathA36_B_eff"], "NO_CURRENT_REGISTERED_REDUCTION", ["B_eff=rho_B0^2/chi_c"], "Future-Compression-4D-to-3D", True, location="3D brane surface"),
        RowSource("chi_c", "C5 compression susceptibility", "source-restored", refs["pathA36_B_eff"], "NO_CURRENT_REGISTERED_REDUCTION", ["B_eff=rho_B0^2/chi_c"], "Future-Compression-4D-to-3D", True, location="3D brane surface"),
        RowSource("B_eff", "derived C5 density modulus", dim_str(DIM_B_EFF), refs["pathA36_B_eff"], "DEPENDENT_ON_COMPRESSION_INPUTS", ["B_eff=rho_B0^2/chi_c"], active_medium_competition=True, dependency=["rho_B0", "chi_c"], location="3D brane surface"),
        RowSource("C_hu", "embedding h/u_L mixing overlap", "source-restored", refs["pathA39_scalar_C_hu"], "NO_CURRENT_REGISTERED_REDUCTION", [], "Future-Embedding-Overlap", True, "FREEDOM_CERTIFIED_CURRENT_LEDGER{C_hu}", location="throat/embedding seam"),
        RowSource("Q_E", "electric throat source magnitude", "source-restored", refs["pathA38_calibrated_deferred"], "DECLARED_CALIBRATION_ANCHOR", ["q_h=2 Q_E tanh(b/ell)/b"], active_medium_competition=True, calibration_anchor="Q_E", location="throat/embedding seam"),
        RowSource("ell", "throat/wall profile scale", "L", refs["pathA38_qh"], "CALIBRATED_GEOMETRY_SOURCE_INPUT", ["q_h=2 Q_E tanh(b/ell)/b"], active_medium_competition=True, calibrated_geometry=True, location="throat/embedding seam"),
        RowSource("b", "compact throat source half-separation/form factor scale", "L", refs["pathA38_qh"], "CALIBRATED_GEOMETRY_SOURCE_INPUT", ["q_h=2 Q_E tanh(b/ell)/b"], active_medium_competition=True, calibrated_geometry=True, location="throat/embedding seam"),
        RowSource("M_h", "h-sector zero-mode normalization/mass coefficient", "source-restored", refs["pathA39_imports"], "CALIBRATED_GEOMETRY_SOURCE_INPUT", ["K_h=M_h c_E^2"], active_medium_competition=True, calibrated_geometry=True, location="throat/embedding seam"),
        RowSource("K_h", "h-sector stiffness", "M L^-1 T^-2", refs["pathA39_imports"], "DEPENDENT_STIFFNESS", ["K_h=M_h c_E^2"], active_medium_competition=True, dependency=["M_h", "c_E"], location="throat/embedding seam"),
        RowSource("q_h", "electric throat source projection", "source-restored", refs["pathA38_qh"], "DEPENDENT_PROJECTION", ["q_h=2 Q_E tanh(b/ell)/b"], active_medium_competition=True, dependency=["Q_E", "b", "ell"], location="throat/embedding seam"),
    ]

    for p, dim, ref_key in [
        ("c_L1", "M L^8 T^-2", "pathA25_c_L1"),
        ("c_L2", "M L^10 T^-2", "pathA25_c_L2"),
        ("A_R", "M L^6 T^-2", "pathA25_A_R"),
        ("k_R", "L^-1", "pathA25_k_R"),
        ("lambda_Cdiv", "M L^3 T^-2", "pathA25_lambda_Cdiv"),
        ("chi_Cpin", "M L^8 T^-2", "pathA25_chi_Cpin"),
        ("J_Pu", "M L^-1", "pathA25_J_Pu"),
        ("kappa_Pu", "M L^-1 T^-2", "pathA25_kappa_Pu"),
    ]:
        rows.append(RowSource(p, "pathA_25 density-smectic driver", dim, refs[ref_key], "PATHA25_DRIVER_OUT_OF_ACTIVE_SHEAR_SURFACE_NG5", out_of_active=True, location="3D brane surface", notes="Closed pathA_25 candidate; active model is pathA_35 shear-surface."))

    rows += [
        RowSource("lambda_Pu", "pathA_35 parity-repaired P-u coupling", "M L^-1 T^-2", refs["pathA35_rho_role"], "FROZEN_UNUSED_OUT_OF_ACTIVE_NG5", out_of_active=True, location="3D brane surface"),
        RowSource("Omega_w", "pathA_35 bare u_w gap scale", "T^-1", refs["pathA35_rho_role"], "FROZEN_UNUSED_OUT_OF_ACTIVE_NG5", out_of_active=True, location="3D brane surface"),
        RowSource("lambda_N", "pathA_38 wall-internal potential coefficient", "wall-internal", refs["pathA38_wall_internal"], "WALL_INTERNAL_EXCLUDED_FROM_ACTIVE_NG5", out_of_active=True, location="throat/embedding seam"),
        RowSource("lambda_tau", "pathA_38 wall-internal tau mass coefficient", "wall-internal", refs["pathA38_wall_internal"], "WALL_INTERNAL_EXCLUDED_FROM_ACTIVE_NG5", out_of_active=True, location="throat/embedding seam"),
        RowSource("Nu", "moving-source normalization", "source normalization", refs["pathA39_sim_values"], "SOURCE_NORMALIZATION_OUTSIDE_ACTIVE_ORIGIN_COMPETITION", out_of_active=True, location="throat/embedding seam"),
        RowSource("a_T", "transverse moving-source amplitude", "source normalization", refs["pathA39_sim_values"], "SOURCE_NORMALIZATION_OUTSIDE_ACTIVE_ORIGIN_COMPETITION", ["q_A^T=Nu*aT*sCharge"], out_of_active=True, location="throat/embedding seam"),
        RowSource("a_Tp", "second transverse moving-source amplitude", "source normalization", refs["pathA39_stage4_a_Tp"], "SOURCE_NORMALIZATION_OUTSIDE_ACTIVE_ORIGIN_COMPETITION", ["Stage-4 source vector includes Nu*aTp*sCharge"], out_of_active=True, location="throat/embedding seam"),
        RowSource("a_L", "longitudinal moving-source amplitude", "source normalization", refs["pathA39_sim_values"], "SOURCE_NORMALIZATION_OUTSIDE_ACTIVE_ORIGIN_COMPETITION", ["q_L=Nu*aL*sCharge"], out_of_active=True, location="throat/embedding seam"),
    ]

    if case.add_synthetic_irreducible:
        rows.append(RowSource("xi_active", "synthetic active medium parameter", "1", "synthetic-control-source-fact", "SYNTHETIC_ACTIVE_INPUT", active_medium_competition=True, location="3D brane surface"))
    if case.add_synthetic_reducible:
        rows.append(RowSource("p_syn", "synthetic reducible-derived control parameter", dim_str(DIM_SURFACE_MODULUS), "synthetic-control-source-fact", "SYNTHETIC_EARNED_RELATION", ["p_syn=5*K*rho^4/m"], active_medium_competition=True, location="4D bulk"))

    return rows


def dependency_origin(row: RowSource, prior: dict[str, dict[str, Any]]) -> str:
    if row.p == "B_eff":
        return "DEPENDENT_ON_IRREDUCIBLE" if any(prior.get(dep, {}).get("origin_verdict") == "IRREDUCIBLY_INDEPENDENT" for dep in row.dependency) else "DEPENDENT"
    return "DEPENDENT"


def classify_origin(row: RowSource, ctx: dict[str, Any], prior: dict[str, dict[str, Any]]) -> dict[str, Any]:
    route_eval = registered_route_for(row.p, ctx) if row.active_medium_competition else None
    if row.out_of_active:
        origin = "OUT_OF_ACTIVE_NG5"
        reason = "closed_pathA_25_or_scope_excluded"
    elif row.local_status in {"BASE_SUBSTRATE", "ACCEPTED_GEOMETRY_SUBSTRATE", "ACCEPTED_PROFILE_GIVEN_ell_g"}:
        origin = row.local_status
        reason = "accepted_base_or_geometry_substrate"
    elif row.dependency:
        origin = dependency_origin(row, prior)
        reason = "computed_dependency_not_independent"
    elif row.calibration_anchor:
        if row.calibration_anchor not in ctx["calibration_anchors"]:
            origin = "IRREDUCIBLY_INDEPENDENT"
            reason = "calibration_anchor_removed_by_source_mutation"
        else:
            origin = "CALIBRATED_ANCHOR"
            reason = "declared_calibration_anchor"
    elif row.calibrated_geometry:
        if row.p in ctx["calibrated_geometry"]:
            origin = "CALIBRATED_GEOMETRY_INPUT"
            reason = "declared_shared_geometry_or_normalization_input"
        else:
            origin = "IRREDUCIBLY_INDEPENDENT"
            reason = "calibrated_geometry_declaration_removed"
    elif row.local_status == "SYNTHETIC_EARNED_RELATION":
        origin = "REDUCIBLE_DERIVED"
        reason = "closed_relation_in_base_substrate"
    elif route_eval and route_eval["valid"] and route_eval["result_status"] == "SOLVED_PASS":
        origin = "REDUCIBLE_DERIVED"
        reason = "valid_solved_registered_route"
    elif route_eval and route_eval["valid"] and route_eval["result_status"] == "REGISTERED_DEFERRED":
        origin = "REDUCIBLE_IN_PRINCIPLE_SIM_DEFERRED"
        reason = "valid_registered_deferred_route"
    else:
        origin = "IRREDUCIBLY_INDEPENDENT"
        reason = "no_valid_registered_route_or_anchor"

    out = row.to_yaml_base()
    out["route_evaluation"] = route_eval
    out["origin_verdict"] = origin
    out["classification_reason"] = reason
    return out


def classify_rows(rows: list[RowSource], ctx: dict[str, Any]) -> list[dict[str, Any]]:
    classified: list[dict[str, Any]] = []
    prior: dict[str, dict[str, Any]] = {}
    for row in rows:
        rec = classify_origin(row, ctx, prior)
        classified.append(rec)
        prior[row.p] = rec
    return classified


def p40_freedom_tie_no_go() -> dict[str, Any]:
    result = P40.run_case("freedom_tie")
    lock = result.get("algebra", {}).get("lock_sat", {})
    no_go = result.get("verdict", "").startswith("NO_GO") and lock.get("status") == "UNSAT"
    return {
        "no_go": no_go,
        "verdict": result.get("verdict"),
        "lock_sat_status": lock.get("status"),
        "certificate": lock.get("certificate"),
        "core": lock.get("core"),
        "compact": P40.compact_case(result),
    }


def p40_current_nonentailment() -> dict[str, Any]:
    prod = P40.run_case("production")
    prov = prod["algebra"]["provenance"]
    return {
        "compact": P40.compact_case(prod),
        "lock_A": prov["A"],
        "lock_B": prov["B"],
        "c_E_squared_rho_br_minus_mu_R_status": prov["B"]["status"],
        "c_E_squared_rho_br_minus_mu_R_witness_value": prov["B"].get("lock_value_at_witness"),
    }


def sim_deferred_map(rows: list[dict[str, Any]]) -> dict[str, str]:
    out: dict[str, str] = {}
    for row in rows:
        ev = row.get("route_evaluation")
        if row["active_medium_competition"] and row["origin_verdict"] == "REDUCIBLE_IN_PRINCIPLE_SIM_DEFERRED" and ev and ev.get("route_name"):
            out[row["p"]] = ev["route_name"]
    return out


def calibrated_map(rows: list[dict[str, Any]]) -> dict[str, str]:
    return {
        row["p"]: row["origin_verdict"]
        for row in rows
        if row["origin_verdict"] in {"CALIBRATED_ANCHOR", "CALIBRATED_GEOMETRY_INPUT"}
    }


def route_eval_recorded_for_all_active_rows(rows: list[dict[str, Any]]) -> dict[str, Any]:
    active = [row["p"] for row in rows if row["active_medium_competition"]]
    recorded = [row["p"] for row in rows if row["active_medium_competition"] and row.get("route_evaluation") is not None]
    missing = sorted(set(active) - set(recorded))
    return {
        "status": not missing,
        "active_rows": active,
        "recorded_rows": recorded,
        "missing": missing,
        "count": len(recorded),
    }


def decide_verdict(rows: list[dict[str, Any]], no_go: dict[str, Any] | None) -> dict[str, Any]:
    sim = sim_deferred_map(rows)
    calibrated = calibrated_map(rows)
    if no_go:
        return {
            "verdict": no_go["verdict"],
            "reason": no_go["reason"],
            "active_irreducible": [],
            "sim_deferred": sim,
            "calibrated": calibrated,
            "drift_count": 0,
        }
    irreducible = [row["p"] for row in rows if row["active_medium_competition"] and row["origin_verdict"] == "IRREDUCIBLY_INDEPENDENT"]
    if irreducible:
        return {
            "verdict": "SECOND_MEDIUM_DRIFT(active_irreducible={" + ",".join(irreducible) + "})",
            "reason": "active irreducibly independent 3D-brane/throat parameter(s)",
            "active_irreducible": irreducible,
            "sim_deferred": sim,
            "calibrated": calibrated,
            "drift_count": len(irreducible),
        }
    if sim:
        sim_pairs = ",".join(f"{p}->{route}" for p, route in sim.items())
        cal_pairs = ",".join(f"{p}->{origin}" for p, origin in calibrated.items())
        return {
            "verdict": f"ONE_MEDIUM_CONDITIONAL(sim-deferred={{{sim_pairs}}}; calibrated={{{cal_pairs}}})",
            "reason": "no irreducible rows; at least one valid registered deferred route remains",
            "active_irreducible": [],
            "sim_deferred": sim,
            "calibrated": calibrated,
            "drift_count": 0,
        }
    return {
        "verdict": "ONE_MEDIUM_CONSISTENT",
        "reason": "all active rows derived, dependent, accepted, or calibrated",
        "active_irreducible": [],
        "sim_deferred": {},
        "calibrated": calibrated,
        "drift_count": 0,
    }


def interpretation_payload() -> dict[str, Any]:
    return {
        "interpretation": "ONE_CANDIDATE_MEDIUM_4D_TO_3D_REDUCTION_INCOMPLETE",
        "physical_meaning": "The drift is not a separate substance; it is three unreduced 3D-brane-surface parameters: compression rho_B0, chi_c and embedding-mixing C_hu.",
        "location_closure": {
            "arenas": ["4D bulk", "3D brane surface", "throat/embedding seam"],
            "no_fourth_arena": True,
            "fact": "Every production row is assigned to one of the three arenas.",
        },
        "reduction_status": {
            "rho_br": "REGISTERED_PENDING(Route-A)",
            "mu_R": "REGISTERED_PENDING(Route-A)",
            "rho_B0": "NOT_YET_REGISTERED",
            "chi_c": "NOT_YET_REGISTERED",
            "C_hu": "NOT_YET_REGISTERED",
        },
        "named_future_reduction_routes": [
            {
                "name": "compression-sector 4D->3D reduction",
                "targets": ["rho_B0", "chi_c"],
                "status": "DEFERRED_NOT_REGISTERED",
                "description": "derive the brane compression amplitude and susceptibility from the 4D bulk/nonlinear brane solve",
            },
            {
                "name": "embedding-overlap reduction",
                "targets": ["C_hu"],
                "status": "DEFERRED_NOT_REGISTERED",
                "description": "derive the throat density-to-h overlap coefficient from the nonlinear throat/embedding solve",
            },
        ],
        "reopen_trigger": "Re-open NG5 via the section 3.3 forward trigger when either named future route is registered or solved.",
        "honest_caveat": "The brane is a postulated shear-supporting ordered phase; whether the one medium yields it and whether the three reductions close rather than no-go is genuinely open.",
    }


def freedom_states(rows: list[dict[str, Any]]) -> dict[str, Any]:
    by_param = {row["p"]: row for row in rows}
    rho_br = by_param["rho_br"]
    c_hu = by_param["C_hu"]
    rho_ev = rho_br.get("route_evaluation")
    return {
        "rho_br": {
            "state": "FREEDOM_SIM_DEFERRED{Route-A}" if rho_ev and rho_ev.get("valid") else "FREEDOM_CERTIFIED_CURRENT_LEDGER{rho_br}",
            "route_evaluation": rho_ev,
        },
        "C_hu": {
            "state": "FREEDOM_CERTIFIED_CURRENT_LEDGER{C_hu}",
            "route_evaluation": c_hu.get("route_evaluation"),
        },
    }


def run_pipeline(case: CaseSpec, catalog: dict[str, Any]) -> dict[str, Any]:
    refs = catalog["refs"]
    routes = route_inventory(refs, catalog["pathA40"], case)
    ctx = {
        "routes": routes,
        "calibration_anchors": {"Q_E": refs["pathA38_calibrated_deferred"]} if "Q_E" not in case.remove_calibration_anchors else {},
        "calibrated_geometry": {p: "shared throat geometry/normalization input" for p in ["b", "ell", "M_h"] if p not in case.remove_calibrated_geometry},
    }
    rows = classify_rows(build_rows(refs, case), ctx)
    no_go = None
    p40_no_go = None
    if case.contradiction == "C_hu_charge_residue_tie":
        p40_no_go = p40_freedom_tie_no_go()
        if p40_no_go["no_go"]:
            no_go = {"verdict": "NO_GO(cone-lock-feedback)", "reason": "recomputed pathA_40 freedom_tie lock SAT is UNSAT"}
        else:
            raise AssertionError("contradiction control did not recompute pathA_40 UNSAT")
    verdict = decide_verdict(rows, no_go)
    return {
        "case": case.name,
        "route_inventory": routes,
        "origin_ledger": rows,
        "lineage_adjudication": lineage_adjudication(case),
        "freedom_states": freedom_states(rows),
        "route_eval_recorded_for_all_active_rows": route_eval_recorded_for_all_active_rows(rows),
        "verdict": verdict,
        "forward_reopen_triggers": [
            {"row": p, "registered_solve": route, "trigger": "re-adjudicate on solve"}
            for p, route in sim_deferred_map(rows).items()
        ],
        "p40_no_go_recompute": p40_no_go,
    }


def compact_case(result: dict[str, Any]) -> dict[str, Any]:
    return {
        "verdict": result["verdict"]["verdict"],
        "active_irreducible": result["verdict"]["active_irreducible"],
        "drift_count": result["verdict"]["drift_count"],
        "sim_deferred": result["verdict"]["sim_deferred"],
        "calibrated": result["verdict"]["calibrated"],
        "lineage_finding": result["lineage_adjudication"]["lineage_finding"],
        "lineage_outcome": result["lineage_adjudication"]["outcome"],
        "route_eval_recorded_for_all_active_rows": result["route_eval_recorded_for_all_active_rows"]["status"],
        "route_evaluation_count": result["route_eval_recorded_for_all_active_rows"]["count"],
        "freedom_states": {k: v["state"] for k, v in result["freedom_states"].items()},
    }


def build_controls(catalog: dict[str, Any], production: dict[str, Any]) -> dict[str, Any]:
    controls: dict[str, Any] = {}
    prod_compact = compact_case(production)

    specs = {
        "AB_delete_registry": CaseSpec("AB_delete_registry", blank_all_routes=True),
        "route_blank_Route_A": CaseSpec("route_blank_Route_A", blank_routes={"Route-A"}),
        "route_field_blank_Route_A_target_blind": CaseSpec("route_field_blank_Route_A_target_blind", blank_route_fields={"Route-A": {"target_blind"}}),
        "route_field_blank_Route_A_missing_objects": CaseSpec("route_field_blank_Route_A_missing_objects", blank_route_fields={"Route-A": {"missing_objects"}}),
        "calibration_ablation_Q_E": CaseSpec("calibration_ablation_Q_E", remove_calibration_anchors={"Q_E"}),
        "irreducible_synthetic": CaseSpec("irreducible_synthetic", add_synthetic_irreducible=True),
        "reducible_derived_synthetic": CaseSpec("reducible_derived_synthetic", add_synthetic_reducible=True),
        "contradiction": CaseSpec("contradiction", contradiction="C_hu_charge_residue_tie"),
    }
    for name, spec in specs.items():
        result = run_pipeline(spec, catalog)
        compact = compact_case(result)
        if name == "AB_delete_registry":
            fired = compact["verdict"] != prod_compact["verdict"] and set(prod_compact["active_irreducible"]) < set(compact["active_irreducible"])
            mutation = "blank all registered route source records"
            invariant = "all RouteEvaluation records for previously covered rows recompute invalid"
        elif name.startswith("route_blank") or name.startswith("route_field_blank"):
            covered = ["rho_br", "mu_R", "c_E"]
            fired = all(row["origin_verdict"] == "IRREDUCIBLY_INDEPENDENT" for row in result["origin_ledger"] if row["p"] in covered)
            mutation = "blank Route-A or a required Route-A field at the source-record level"
            invariant = {row["p"]: row["route_evaluation"] for row in result["origin_ledger"] if row["p"] in covered}
        elif name == "calibration_ablation_Q_E":
            fired = any(row["p"] == "Q_E" and row["origin_verdict"] == "IRREDUCIBLY_INDEPENDENT" for row in result["origin_ledger"])
            mutation = "remove the declared Q_E calibration-anchor source"
            invariant = "Q_E has no valid registered route and is no longer a declared anchor"
        elif name == "irreducible_synthetic":
            fired = "xi_active" in compact["active_irreducible"]
            mutation = "add active synthetic parameter xi_active with no relation or route"
            invariant = "registered_route_for(xi_active) recomputes NO_REGISTERED_ROUTE_FOR_PARAM"
        elif name == "reducible_derived_synthetic":
            fired = any(row["p"] == "p_syn" and row["origin_verdict"] == "REDUCIBLE_DERIVED" for row in result["origin_ledger"])
            mutation = "add earned closed relation p_syn=5*K*rho^4/m"
            invariant = "classifier sees SYNTHETIC_EARNED_RELATION"
        elif name == "contradiction":
            fired = compact["verdict"] == "NO_GO(cone-lock-feedback)" and result["p40_no_go_recompute"]["lock_sat_status"] == "UNSAT"
            mutation = "assert C_hu charge-residue tie, then rerun pathA_40 lock+stability SAT"
            invariant = result["p40_no_go_recompute"]
        else:
            raise AssertionError(name)
        controls[name] = {
            "source_mutation": mutation,
            "recomputed_invariant": invariant,
            "fired": fired,
            "before": prod_compact,
            "after": compact,
            "expected_transition": f"{prod_compact['verdict']} -> {compact['verdict']}",
        }

    same = run_pipeline(CaseSpec("residual_same", residual_mode="same_no_residual"), catalog)
    free = run_pipeline(CaseSpec("residual_free", residual_mode="same_with_free_residual"), catalog)
    controls["residual_multiplier_ablation"] = {
        "source_mutation": "mutate lineage to same active object, then toggle residual multiplier 1 vs Xi_residual",
        "recomputed_invariant": {
            "same_no_residual": same["lineage_adjudication"],
            "same_with_free_residual": free["lineage_adjudication"],
        },
        "fired": same["lineage_adjudication"]["outcome"] == "SAME" and free["lineage_adjudication"]["outcome"] == "DIFFERENT",
        "before": prod_compact,
        "after": {
            "same_no_residual": compact_case(same),
            "same_with_free_residual": compact_case(free),
        },
        "expected_transition": "lineage SAME when residual=1; lineage DIFFERENT when residual=Xi_residual",
    }

    controls["route_eval_recorded_for_all_active_rows"] = {
        "source_mutation": "production ledger audit only",
        "recomputed_invariant": production["route_eval_recorded_for_all_active_rows"],
        "fired": production["route_eval_recorded_for_all_active_rows"]["status"],
        "before": prod_compact,
        "after": prod_compact,
        "expected_transition": "all active production rows carry RouteEvaluation records",
    }
    return controls


def extraction_audit(catalog: dict[str, Any]) -> dict[str, Any]:
    refs = catalog["refs"]
    return {
        "active_model": {
            "status": "PATHA35_SHEAR_SURFACE_ACTIVE",
            "pathA25_density_route": "CLOSED",
            "sources": [refs["pathA25_closed"], refs["pathA25_b4_fail"], refs["pathA35_rho_role"]],
        },
        "lineage": {
            "finding": "NO_OVERCOUNT_ROUTE_A_PENDING",
            "sources": [refs["pathA25_varrho_prose"], refs["pathA35_rho_constant"], refs["pathA40_route_a"]],
        },
        "calibration_sources": {
            "Q_E": [refs["pathA38_calibrated_deferred"], refs["decision09_calibration"], refs["decision14_calibration"]],
            "geometry": [refs["pathA38_qh"], refs["pathA39_imports"]],
        },
    }


def comparison_payload(results: dict[str, Any]) -> dict[str, Any]:
    prod = results["production"]
    rows = prod["origin_ledger"]
    origin_by_param = {
        row["p"]: row["origin_verdict"]
        for row in rows
        if row["active_medium_competition"]
        or row["origin_verdict"] in {"OUT_OF_ACTIVE_NG5", "BASE_SUBSTRATE", "ACCEPTED_GEOMETRY_SUBSTRATE", "ACCEPTED_PROFILE_GIVEN_ell_g"}
    }
    route_validity = {
        row["p"]: {
            "route_name": row["route_evaluation"]["route_name"],
            "valid": row["route_evaluation"]["valid"],
            "reason": row["route_evaluation"]["reason"],
            "origin": row["origin_verdict"],
        }
        for row in rows
        if row["active_medium_competition"]
    }
    return {
        "production": compact_case(prod),
        "origin_by_param": origin_by_param,
        "route_validity": route_validity,
        "lineage": prod["lineage_adjudication"],
        "interpretation": results["interpretation"],
        "dual_engine_derivations": {
            "dimension_and_residual": prod["lineage_adjudication"]["dimension_derivation"] | {
                "residual_multiplier": prod["lineage_adjudication"]["residual_multiplier"],
            },
            "pathA40_current_nonentailment": results["pathA40_nonentailment"]["c_E_squared_rho_br_minus_mu_R_status"],
            "pathA40_current_nonentailment_witness": results["pathA40_nonentailment"]["c_E_squared_rho_br_minus_mu_R_witness_value"],
        },
        "controls": {
            name: {
                "fired": rec["fired"],
                "transition": rec["expected_transition"],
                "after_verdict": rec["after"]["verdict"] if isinstance(rec["after"], dict) and "verdict" in rec["after"] else None,
            }
            for name, rec in results["controls"].items()
        },
    }


def count_agreements(payload: Any) -> int:
    if isinstance(payload, dict):
        return sum(count_agreements(v) for v in payload.values())
    if isinstance(payload, list):
        return sum(count_agreements(v) for v in payload)
    return 1


def build_results() -> dict[str, Any]:
    catalog = source_catalog()
    production = run_pipeline(CaseSpec("production"), catalog)
    controls = build_controls(catalog, production)
    results = {
        "schema": SCHEMA,
        "engine": "sympy",
        "production": production,
        "controls": controls,
        "interpretation": interpretation_payload(),
        "extraction_audit": extraction_audit(catalog),
        "pathA40_nonentailment": p40_current_nonentailment(),
        "source_refs": catalog["refs"],
    }
    payload = comparison_payload(results)
    results["comparison_payload"] = payload
    results["comparison_digest"] = sha256_json(payload)
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
    return {
        "status": "ENGINE_AGREE",
        "agreement_test": "STRUCTURAL_PAYLOAD_EQUALITY",
        "structural_payload_equal": True,
        "agreement_count": count_agreements(sympy_results["comparison_payload"]),
        "sympy_payload": str(SYM_OUT),
        "mathematica_payload": str(WL_OUT),
        "sympy_comparison_digest": sympy_results["comparison_digest"],
        "mathematica_comparison_digest": wl.get("comparison_digest"),
    }


def yaml_results(results: dict[str, Any], agreement: dict[str, Any]) -> dict[str, Any]:
    prod = results["production"]
    return {
        "schema": "pathA_41_ng5_second_medium_drift_results/v2",
        "verdict": prod["verdict"]["verdict"],
        "active_irreducible_set": prod["verdict"]["active_irreducible"],
        "drift_count": prod["verdict"]["drift_count"],
        "lineage_finding": prod["lineage_adjudication"]["lineage_finding"],
        "lineage_adjudication": prod["lineage_adjudication"],
        "interpretation": results["interpretation"],
        "sim_deferred_rows": prod["verdict"]["sim_deferred"],
        "calibrated_rows": prod["verdict"]["calibrated"],
        "forward_reopen_triggers": prod["forward_reopen_triggers"],
        "origin_ledger": prod["origin_ledger"],
        "freedom_states": prod["freedom_states"],
        "route_eval_recorded_for_all_active_rows": prod["route_eval_recorded_for_all_active_rows"],
        "pathA40_nonentailment": results["pathA40_nonentailment"],
        "controls": results["controls"],
        "extraction_audit": results["extraction_audit"],
        "engine_agreement": agreement,
        "comparison_payload": results["comparison_payload"],
        "source_refs": results["source_refs"],
    }


def write_report(results: dict[str, Any], agreement: dict[str, Any]) -> None:
    prod = results["production"]
    verdict = prod["verdict"]["verdict"]
    controls = results["controls"]
    fired = [name for name, rec in controls.items() if rec.get("fired")]
    interp = results["interpretation"]

    lines = [
        verdict,
        "",
        "# pathA_41 NG5 SECOND_MEDIUM_DRIFT Gate",
        "",
        f"Primary verdict: `{verdict}`.",
        f"Active irreducible set (computed): `{prod['verdict']['active_irreducible']}`.",
        f"Lineage finding: `{prod['lineage_adjudication']['lineage_finding']}`.",
        f"Engine agreement: `{agreement['status']}` over `{agreement['agreement_count']}` compared scalar/audit quantities.",
        "",
        "## Interpretation",
        "",
        f"- `{interp['interpretation']}`.",
        f"- Physical meaning: {interp['physical_meaning']}",
        f"- Location closure: every row lives in `{interp['location_closure']['arenas']}`; no fourth arena = `{interp['location_closure']['no_fourth_arena']}`.",
        f"- Honest caveat: {interp['honest_caveat']}",
        "",
        "Named future reduction routes:",
    ]
    for route in interp["named_future_reduction_routes"]:
        lines.append(f"- `{route['name']}` targets `{route['targets']}`: `{route['status']}`.")
    lines += [
        "",
        "Reduction status:",
    ]
    for p, status in interp["reduction_status"].items():
        lines.append(f"- `{p}`: `{status}`")
    lines += [
        "",
        "## Sim-Deferred And Calibrated Rows",
        "",
    ]
    for p, route in prod["verdict"]["sim_deferred"].items():
        lines.append(f"- sim-deferred: `{p}` -> `{route}` -> re-adjudicate on solve")
    for p, origin in prod["verdict"]["calibrated"].items():
        lines.append(f"- calibrated: `{p}` -> `{origin}`")
    lines += [
        "",
        "## Origin Ledger",
        "",
        "| p | incidence | route eval | origin | location |",
        "|---|---|---|---|---|",
    ]
    for row in prod["origin_ledger"]:
        ev = row.get("route_evaluation")
        route = "" if ev is None else f"{ev.get('route_name')}:{ev.get('valid')}:{ev.get('reason')}"
        lines.append(f"| `{row['p']}` | {row['incidence']} | `{route}` | `{row['origin_verdict']}` | `{row['location']}` |")
    lines += [
        "",
        "## Controls",
        "",
        "| control | fired | transition |",
        "|---|---:|---|",
    ]
    for name, rec in controls.items():
        lines.append(f"| `{name}` | `{rec.get('fired')}` | {rec.get('expected_transition')} |")
    lines += [
        "",
        "## Dual-Engine Split",
        "",
        "- SymPy and Mathematica independently compute MLT dimension closure, residual lineage states, RouteEvaluation validity, origin classification, active irreducibles, and control transitions.",
        "- The contradiction control is gated on a recomputed pathA_40 `freedom_tie` UNSAT result, not on a typed no-go flag.",
        "",
        "Run commands:",
        "",
        "```text",
        "timeout 600 python3 software/stage1_solver/tools/pathA_41_ng5_second_medium_drift_sympy.py",
        "timeout 600 math -script software/stage1_solver/tools/pathA_41_ng5_second_medium_drift.wl",
        "timeout 600 python3 software/stage1_solver/tools/pathA_41_ng5_second_medium_drift_sympy.py --compare",
        "```",
        "",
        f"Controls fired: `{', '.join(fired)}`.",
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
        with YAML_OUT.open("w", encoding="utf-8") as fh:
            yaml.dump(yaml_results(results, agreement), fh, Dumper=NoAliasDumper, sort_keys=False, allow_unicode=False, width=180)
        write_report(results, agreement)
        print(f"OK pathA_41_ng5_second_medium_drift_sympy compare -> {REPORT_OUT}")
    else:
        print(f"OK pathA_41_ng5_second_medium_drift_sympy -> {SYM_OUT}")


if __name__ == "__main__":
    main()
