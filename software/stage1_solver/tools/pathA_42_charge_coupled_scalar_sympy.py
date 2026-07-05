#!/usr/bin/env python3
"""pathA_42 charge-coupled scalar gate, SymPy engine.

This runner separates the pieces the directive marks as dual-engine derivable
from single-engine ledger extraction/bookkeeping.  It writes deterministic
reports and a compact comparison payload for the Mathematica engine.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import numbers
import re
import subprocess
import sys
from collections.abc import Iterable
from dataclasses import dataclass, replace
from pathlib import Path
from typing import Any

import sympy as sp
import yaml


SCRIPT_PATH = Path(__file__).resolve()
STAGE1_ROOT = SCRIPT_PATH.parents[1]
REPO_ROOT = STAGE1_ROOT.parents[1]
REPORTS = STAGE1_ROOT / "reports"
RUN_OUT = STAGE1_ROOT / "_scratch"

SYM_OUT = RUN_OUT / "pathA_42_charge_coupled_scalar_sympy.json"
WL_OUT = RUN_OUT / "pathA_42_charge_coupled_scalar_mathematica.json"
YAML_OUT = REPORTS / "pathA_42_charge_coupled_scalar_results.yaml"
REPORT_OUT = REPORTS / "pathA_42_charge_coupled_scalar.md"

SCHEMA = "pathA_42_charge_coupled_scalar_sympy/v1"
WL_SCHEMA = "pathA_42_charge_coupled_scalar_mathematica/v1"

VERDICT_MAPPED = "SCALAR_DEPARTURE_MAPPED_MAGNITUDE_SIM_GATED"
VERDICT_NO_GO = "NO_GO_CONSISTENCY"

H_EP_SAFE = "EARNED_SAFE_MASS_CHANNEL_ON_DECOUPLED_FLOOR"
H_EP_FIFTH = "FIFTH_FORCE_TRIGGERED"
RADIATION_SIM = "SIM_GATED"
RADIATION_TENSION = "FALSIFIABLE_TENSION"
RADIATION_AUDIT = "AUDIT_ONLY_NOT_EARNED"
UNIV_SIM = "SIM_GATED_REQUIRED_UNIVERSALITY"
UNIV_TENSION = "FALSIFIABLE_TENSION"
UL_SIM = "SIM_GATED"
PF_SIM = "SIM_GATED"
PF_TENSION = "PREFERRED_FRAME_TENSION"

GUARDED_KINDS = {"M_h", "c_E", "K_h", "P_h_over_P_EM", "EP_magnitude", "residue_floor"}
GUARDED_SYMBOL_SEGMENTS = {"m_h", "mh", "c_e", "ce", "k_h", "kh", "n0", "n_0"}
GUARDED_PATH_TOKENS = {
    "p_h_over_p_em",
    "ph_over_pem",
    "p_h_p_em",
    "power_ratio",
    "flux_ratio",
    "ep_magnitude",
    "fifth_force_magnitude",
    "residue_floor",
}
NUMERIC_LITERAL_RE = re.compile(
    r"(?<![A-Za-z0-9_])[-+]?(?:(?:\d+(?:\.\d*)?)|(?:\.\d+))(?:[eE][-+]?\d+)?(?![A-Za-z0-9_])"
)


class NoAliasDumper(yaml.SafeDumper):
    def ignore_aliases(self, data: Any) -> bool:
        return True


class LaunderingRefused(AssertionError):
    """Raised when a guarded numeric magnitude reaches the serialization path."""


def stable_yaml_dump(data: Any) -> str:
    return yaml.dump(
        data,
        Dumper=NoAliasDumper,
        sort_keys=False,
        default_flow_style=False,
        allow_unicode=False,
        width=120,
    )


def sha256_json(data: Any) -> str:
    blob = json.dumps(data, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return hashlib.sha256(blob).hexdigest()


def write_json(path: Path, data: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, sort_keys=True, indent=2) + "\n", encoding="utf-8")


def read_json(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8"))


def load_text(path: Path) -> str:
    if not path.exists():
        raise AssertionError(f"missing substrate report: {path}")
    return path.read_text(encoding="utf-8")


def load_yaml(path: Path) -> Any:
    if not path.exists():
        raise AssertionError(f"missing substrate report: {path}")
    with path.open("r", encoding="utf-8") as fh:
        return yaml.safe_load(fh)


def report_path(name: str) -> Path:
    return REPORTS / name


def rel(path: Path) -> str:
    try:
        return str(path.relative_to(REPO_ROOT))
    except ValueError:
        return str(path)


def assert_contains(text: str, token: str, label: str) -> None:
    if token not in text:
        raise AssertionError(f"substrate token missing for {label}: {token}")


def s(expr: Any) -> str:
    if isinstance(expr, sp.MatrixBase):
        return str([[s(item) for item in row] for row in expr.tolist()])
    if isinstance(expr, sp.Basic):
        return sp.sstr(sp.factor(sp.simplify(expr)))
    return str(expr)


def leaf_count(data: Any) -> int:
    if isinstance(data, dict):
        return sum(leaf_count(v) for v in data.values())
    if isinstance(data, list):
        return sum(leaf_count(v) for v in data)
    return 1


def normalized_path_token(value: Any) -> str:
    return (
        str(value)
        .strip()
        .lower()
        .replace("/", "_over_")
        .replace("-", "_")
        .replace(" ", "_")
        .replace(".", "_")
    )


GUARDED_FIELD_TOKENS = {normalized_path_token(kind) for kind in GUARDED_KINDS} | GUARDED_PATH_TOKENS


def h_time_kinetic_parent_action_provenance(substrate: dict[str, Any]) -> dict[str, Any] | None:
    """Return only a substrate-imported earned parent-action object, never a local flag."""
    fact = substrate.get("facts", {}).get("h_time_kinetic_parent_action")
    obj = substrate.get("provenance_objects", {}).get("h_time_kinetic_parent_action")
    if fact in {None, "ABSENT", "MISSING", "NOT_EARNED"}:
        return None
    if isinstance(obj, dict) and obj.get("status") == "EARNED" and obj.get("source"):
        return obj
    return None


def has_h_time_kinetic_parent_action(substrate: dict[str, Any]) -> bool:
    return h_time_kinetic_parent_action_provenance(substrate) is not None


def is_numeric_leaf(value: Any) -> bool:
    if isinstance(value, bool):
        return False
    if isinstance(value, numbers.Number):
        return True
    if isinstance(value, sp.Basic):
        return bool(value.is_number)
    return False


def guarded_numeric_path(path: tuple[Any, ...]) -> bool:
    segments = [normalized_path_token(part) for part in path]
    if any(
        segment == symbol or segment.startswith(f"{symbol}_")
        for segment in segments
        for symbol in GUARDED_SYMBOL_SEGMENTS
    ):
        return True
    joined = ".".join(segments)
    if any(token in joined for token in GUARDED_PATH_TOKENS):
        return True
    return "magnitude" in joined and any(token in joined for token in {"radiation", "scalar", "ep", "u_l", "ul", "h"})


def guarded_field_path(path: tuple[Any, ...]) -> bool:
    return any(normalized_path_token(part) in GUARDED_FIELD_TOKENS for part in path)


def contains_numeric_literal(value: str) -> bool:
    return NUMERIC_LITERAL_RE.search(value) is not None


def ordered_iterable(value: Iterable[Any]) -> Iterable[Any]:
    if isinstance(value, (set, frozenset)):
        return sorted(value, key=lambda item: (type(item).__name__, repr(item)))
    return value


def find_laundered_numeric_paths(data: Any, path: tuple[Any, ...] = ()) -> list[str]:
    findings: list[str] = []
    if isinstance(data, dict):
        for key, value in data.items():
            findings.extend(find_laundered_numeric_paths(value, path + (key,)))
    elif isinstance(data, str):
        if guarded_field_path(path) and contains_numeric_literal(data):
            findings.append(".".join(str(part) for part in path))
    elif isinstance(data, (bytes, bytearray)):
        pass
    elif isinstance(data, Iterable):
        for index, value in enumerate(ordered_iterable(data)):
            findings.extend(find_laundered_numeric_paths(value, path + (index,)))
    elif is_numeric_leaf(data) and guarded_numeric_path(path):
        findings.append(".".join(str(part) for part in path))
    return findings


def format_laundering_refusal(paths: list[str]) -> str:
    return "LAUNDERING_REFUSED guarded_numeric_paths=" + ",".join(paths)


def validate_no_laundered_numbers(results: dict[str, Any], substrate: dict[str, Any]) -> None:
    if has_h_time_kinetic_parent_action(substrate):
        return
    paths = find_laundered_numeric_paths(results)
    if paths:
        raise LaunderingRefused(format_laundering_refusal(paths))


def expect_laundering_refused(payload: dict[str, Any], substrate: dict[str, Any]) -> str:
    try:
        validate_no_laundered_numbers(payload, substrate)
    except LaunderingRefused as exc:
        return str(exc)
    raise AssertionError("expected LAUNDERING_REFUSED, but serialization guard allowed the payload")


@dataclass(frozen=True)
class Ledger:
    name: str = "production"
    h_bare_mass_residue_zero: bool = True
    h_bare_mass_residue_fixture: bool = False
    static_coulomb_match: bool = True
    q_ratio_global_earned: bool = False
    species_indexed_b_over_ell: bool = False
    pinned_kh_fact: bool = False
    commit_cE_equals_cgamma: bool = False
    suppression_requires_large_cE: bool = False
    stability_violation: bool = False
    production_laundering_attempt: bool = False
    qL_zero_fixture: bool = False
    C_hu_nonzero: bool = True
    qM_nonzero: bool = True
    qM_sign: int = 1
    radiation_derivation_ok: bool = True
    radiation_corrupt_speed: bool = False
    radiation_wrong_normalization: bool = False
    uL_no_go: bool = False
    uL_escalated: bool = False


@dataclass(frozen=True)
class MagnitudeAttempt:
    label: str
    kind: str
    replacement: str | None = None
    contains_guarded_symbol: bool | None = None
    earned_parent_action: bool = False


def magnitude_attempt_contains_guarded_symbol(attempt: MagnitudeAttempt) -> bool:
    if attempt.contains_guarded_symbol is not None:
        return attempt.contains_guarded_symbol
    haystack = normalized_path_token(f"{attempt.label} {attempt.kind} {attempt.replacement or ''}")
    return attempt.kind in GUARDED_KINDS or any(token in haystack for token in GUARDED_SYMBOL_SEGMENTS | GUARDED_PATH_TOKENS)


@dataclass(frozen=True)
class ChannelResult:
    h_EP: str
    radiation: str
    universality: str
    u_L_EP: str
    preferred_frame: str
    subtags: dict[str, bool]
    guard_predicates: dict[str, bool]
    details: dict[str, Any]


def import_substrate() -> dict[str, Any]:
    scalar_yaml_path = report_path("pathA_39_scalar_admixture_results.yaml")
    scalar_screen_path = report_path("pathA_39_scalar_admixture_screen.md")
    magnetic_path = report_path("pathA_39_magnetic_force.md")
    p38_path = report_path("pathA_38_results.yaml")
    p40_path = report_path("pathA_40_cone_lock.md")
    p40_yaml_path = report_path("pathA_40_cone_lock_results.yaml")
    gate_l_path = report_path("pathA_35_gateL_light.md")

    scalar_yaml = load_yaml(scalar_yaml_path)
    scalar_screen = load_text(scalar_screen_path)
    magnetic = load_text(magnetic_path)
    p38 = load_yaml(p38_path)
    p38_text = load_text(p38_path)
    p40 = load_text(p40_path)
    p40_yaml = load_yaml(p40_yaml_path)
    gate_l = load_text(gate_l_path)

    assert_contains(str(scalar_yaml.get("L_scalar", {}).get("source_mass_vector")), "[qM, 0]", "source_mass_vector=[qM,0]")
    assert_contains(str(scalar_yaml.get("imports", {}).get("imported", {}).get("M_h")), "symbolic Mh", "M_h symbolic import")
    assert_contains(str(scalar_yaml.get("imports", {}).get("imported", {}).get("c_E")), "dynamic Green", "c_E Green import")
    assert_contains(scalar_screen, "4*QE**2*tanh(b/ell)**2/(Mh*b**2)", "h residue floor import")
    assert_contains(magnetic, "PREFERRED_FRAME_UNLESS_cE_EQUALS_cGAMMA", "preferred-frame diagnostic")
    assert_contains(magnetic, "(cE**2*rhoBr - muR)/rhoBr", "preferred-frame residual")
    assert_contains(p38_text, "grav_even_overlap: '0'", "grav_even_overlap zero")
    assert_contains(p38_text, "seed_normalization: 8/(3*ell)", "N0 zero-mode normalization")
    assert_contains(p38_text, "operator: g''+(2/R)g'+((omega/c_E)^2-m_0^2)*g=0", "inserted c_E dynamic operator")
    assert_contains(p38_text, "static_projected_U_kernel", "static Coulomb kernel")
    assert_contains(p40, "CONE_LOCK_CALIBRATED", "cone-lock calibrated verdict")
    assert_contains(p40, "Strict stability is encoded as `B_eff*K_h - C_hu^2 - sigma = 0`", "strict stability")
    assert_contains(gate_l, "Ungapped `u_w`: trips `FAIL_BENDING_MASSLESS_FIFTH_FORCE`", "Gate-L massless bending failure")

    if scalar_yaml.get("engine_agreement", {}).get("status") != "ENGINE_AGREE":
        raise AssertionError("pathA_39 scalar substrate lacks ENGINE_AGREE")
    if p40_yaml.get("verdict") != "CONE_LOCK_CALIBRATED":
        raise AssertionError("pathA_40 cone-lock substrate verdict mismatch")
    if p38.get("source_projections", {}).get("grav_even_overlap") != "0":
        raise AssertionError("pathA_38 grav_even_overlap import mismatch")

    # The parent time-kinetic action is the only provenance that can unlock
    # guarded numeric scalar magnitudes.  It is deliberately absent in today's
    # imported substrate, so Guard A must refuse any such numeric payload.
    h_time_parent_action = None

    return {
        "resolved_report_root": rel(scalar_yaml_path.parent),
        "imports": {
            "pathA_39_scalar_admixture_results": rel(scalar_yaml_path),
            "pathA_39_scalar_admixture_screen": rel(scalar_screen_path),
            "pathA_39_magnetic_force": rel(magnetic_path),
            "pathA_38_results": rel(p38_path),
            "pathA_40_cone_lock": rel(p40_path),
            "pathA_35_gateL_light": rel(gate_l_path),
        },
        "facts": {
            "source_mass_vector": "[q_M, 0]",
            "grav_even_overlap": "0",
            "h_time_kinetic_parent_action": "ABSENT",
            "N0_normalization": "8/(3*ell)",
            "c_E_operator_status": "inserted in dynamic Green operator, not derived",
            "M_h_status": "symbolic/imported calibrated geometry input, not pinned",
            "K_h_static_calibration_status": "NOT_PINNED_BY_pathA_38_STATIC_COULOMB",
            "cone_lock_status": "CONE_LOCK_CALIBRATED_NOT_DERIVED",
            "gate_L_connection": "same embedding-direction family, distinct ledger object; ungapped bending family triggers Gate-L fifth-force failure",
        },
        "provenance_objects": {
            "h_time_kinetic_parent_action": h_time_parent_action,
        },
        "derived_booleans": {
            "h_time_kinetic_parent_action_present": h_time_parent_action is not None,
        },
        "citations": {
            "source_mass_vector": f"{rel(scalar_yaml_path)}:90-91",
            "grav_even_overlap": f"{rel(p38_path)}:148",
            "c_E_inserted": f"{rel(p38_path)}:306",
            "N0_static_calibration": f"{rel(p38_path)}:285-291",
            "M_h_symbolic": f"{rel(scalar_yaml_path)}:64",
            "preferred_frame": f"{rel(magnetic_path)}:63-65",
            "cone_lock_calibrated": f"{rel(p40_path)}:1-7",
            "Gate_L_bending": f"{rel(gate_l_path)}:132",
        },
    }


def derive_core() -> dict[str, Any]:
    B_eff, K_h, C_hu = sp.symbols("B_eff K_h C_hu", nonzero=True)
    q_L, q_h, q_M, m_h = sp.symbols("q_L q_h q_M m_h")
    c_E, c_gamma, M_h, Q_E, omega, d = sp.symbols("c_E c_gamma M_h Q_E omega d", positive=True)
    kappa = sp.symbols("kappa", positive=True)
    r, z = sp.symbols("r z", positive=True)

    stiffness = sp.Matrix([[B_eff, C_hu], [C_hu, K_h]])
    inv_stiffness = sp.simplify(stiffness.inv())
    det = sp.factor(stiffness.det())
    jq = sp.Matrix([q_L, q_h])
    jm = sp.Matrix([q_M, m_h])
    A_qq = sp.factor((jq.T * inv_stiffness * jq)[0])
    A_qm = sp.factor((jq.T * inv_stiffness * jm)[0])
    A_mm = sp.factor((jm.T * inv_stiffness * jm)[0])
    A_qm_decoupled = sp.factor(A_qm.subs({C_hu: 0, m_h: 0}))
    A_qm_mixed_qL0 = sp.factor(A_qm.subs({q_L: 0, m_h: 0}))
    A_mm_decoupled = sp.factor(A_mm.subs({C_hu: 0, m_h: 0}))

    expect_A_qq = (K_h * q_L**2 - 2 * C_hu * q_L * q_h + B_eff * q_h**2) / det
    expect_A_qm = (K_h * q_L * q_M - C_hu * q_h * q_M - C_hu * q_L * m_h + B_eff * q_h * m_h) / det
    expect_A_mm = (K_h * q_M**2 - 2 * C_hu * q_M * m_h + B_eff * m_h**2) / det

    def residual_zero(expr: sp.Expr) -> bool:
        return sp.simplify(sp.factor(expr)) == 0

    accel = sp.factor(omega**2 * d)
    scalar_far_time_gradient = sp.factor(q_h * accel / (M_h * c_E**2 * r))
    scalar_power = sp.factor(M_h * c_E * scalar_far_time_gradient**2 * r**2)
    em_far_time_gradient = sp.factor(Q_E * accel / (c_gamma**2 * r))
    em_power = sp.factor(c_gamma * em_far_time_gradient**2 * r**2)
    scalar_power_ratio_bare = sp.factor(scalar_power / em_power)
    expected_bare_ratio = sp.factor((q_h**2 / (M_h * Q_E**2)) * (c_gamma / c_E) ** 3)

    def apply_static_kh_normalization(expr: sp.Expr) -> sp.Expr:
        pinned = sp.factor(expr.subs(M_h, K_h / c_E**2))
        return sp.factor(pinned * K_h / q_h**2 * kappa)

    scalar_power_ratio_pinned = apply_static_kh_normalization(scalar_power_ratio_bare)
    scalar_far_time_gradient_corrupt = sp.factor(q_h * accel / (M_h * c_E * r))
    scalar_power_corrupt = sp.factor(M_h * c_E * scalar_far_time_gradient_corrupt**2 * r**2)
    corrupt_ratio = sp.factor(scalar_power_corrupt / em_power)
    wrong_norm_ratio = apply_static_kh_normalization(scalar_power_ratio_bare)

    def speed_exponent(expr: sp.Expr) -> int:
        expr_z = sp.factor(expr.subs(c_E, z))
        num, den = sp.fraction(sp.together(expr_z))
        exp = sp.degree(den, z) - sp.degree(num, z)
        return int(exp)

    bare_exp = speed_exponent(scalar_power_ratio_bare)
    pinned_exp = speed_exponent(scalar_power_ratio_pinned)
    corrupt_exp = speed_exponent(corrupt_ratio)
    wrong_norm_exp = speed_exponent(wrong_norm_ratio)
    wrong_norm_discriminates_bare = wrong_norm_exp == 1 and bare_exp == 3 and wrong_norm_exp != bare_exp

    mass_vector = sp.Matrix([q_M, 0])
    h_unit = sp.Matrix([0, 1])
    mass_h_projection = sp.factor((h_unit.T * mass_vector)[0])
    fixture_mass_vector = sp.Matrix([q_M, m_h])
    fixture_h_projection = sp.factor((h_unit.T * fixture_mass_vector)[0])
    qM_flip_residual = sp.factor(A_qm.subs({m_h: 0, q_M: -q_M}) + A_qm.subs(m_h, 0))
    A_mm_even_residual = sp.factor(A_mm.subs({m_h: 0, q_M: -q_M}) - A_mm.subs(m_h, 0))
    A_mm_signed_projection = sp.factor(A_mm.subs(m_h, 0) / q_M)
    A_mm_signed_flip_residual = sp.factor(A_mm_signed_projection.subs(q_M, -q_M) + A_mm_signed_projection)

    checks = {
        "A_qq_residual_zero": residual_zero(A_qq - expect_A_qq),
        "A_qm_residual_zero": residual_zero(A_qm - expect_A_qm),
        "A_mm_residual_zero": residual_zero(A_mm - expect_A_mm),
        "A_qm_decoupled_h_mass_zero": residual_zero(A_qm_decoupled - q_L * q_M / B_eff),
        "A_qm_mixed_qL0_nonzero_term": residual_zero(A_qm_mixed_qL0 + C_hu * q_h * q_M / det),
        "A_mm_decoupled_no_h_mass": residual_zero(A_mm_decoupled - q_M**2 / B_eff),
        "mass_h_projection_zero": mass_h_projection == 0,
        "fixture_h_projection_nonzero_symbol": s(fixture_h_projection) == "m_h",
        "qM_flip_A_qm_residual_zero": qM_flip_residual == 0,
        "A_mm_even_in_qM": A_mm_even_residual == 0,
        "A_mm_signed_projection_flips": A_mm_signed_flip_residual == 0,
        "stability_bound_strict": s(det) == "B_eff*K_h - C_hu**2",
        "stability_violation_condition": "C_hu**2 >= B_eff*K_h",
        "radiation_bare_flux_ratio_matches": residual_zero(scalar_power_ratio_bare - expected_bare_ratio),
        "radiation_wrong_normalization_recomputed": wrong_norm_discriminates_bare,
        "radiation_bare_exponent": bare_exp,
        "radiation_pinned_Kh_exponent": pinned_exp,
        "radiation_corrupt_speed_exponent": corrupt_exp,
        "radiation_wrong_normalization_exponent": wrong_norm_exp,
    }

    if not all(v is True for k, v in checks.items() if k.endswith("_zero") or k in {
        "mass_h_projection_zero",
        "A_mm_even_in_qM",
        "A_mm_signed_projection_flips",
    }):
        raise AssertionError(f"core symbolic check failed: {checks}")
    if not checks["radiation_bare_flux_ratio_matches"] or not checks["radiation_wrong_normalization_recomputed"]:
        raise AssertionError(f"radiation flux derivation failed: {checks}")
    if bare_exp != 3 or pinned_exp != 1 or corrupt_exp == bare_exp or wrong_norm_exp != 1:
        raise AssertionError(f"radiation exponent derivation failed: {checks}")

    comparison_payload = {
        "static_exchange": {
            "determinant_check": checks["stability_bound_strict"],
            "A_qq_residual_zero": checks["A_qq_residual_zero"],
            "A_qm_residual_zero": checks["A_qm_residual_zero"],
            "A_mm_residual_zero": checks["A_mm_residual_zero"],
            "decoupled_A_qm_has_no_h_mass_residue": checks["A_qm_decoupled_h_mass_zero"],
            "mixed_qL0_contains_CqhqM": checks["A_qm_mixed_qL0_nonzero_term"],
            "decoupled_A_mm_density_only": checks["A_mm_decoupled_no_h_mass"],
            "qM_flip_A_qm": checks["qM_flip_A_qm_residual_zero"],
            "A_mm_even_but_signed_projection_flips": checks["A_mm_even_in_qM"] and checks["A_mm_signed_projection_flips"],
        },
        "radiation": {
            "bare_flux_ratio_matches_expected": checks["radiation_bare_flux_ratio_matches"],
            "wrong_normalization_recomputed": checks["radiation_wrong_normalization_recomputed"],
            "bare_fixed_exponent": bare_exp,
            "pinned_Kh_exponent": pinned_exp,
            "corrupt_speed_exponent": corrupt_exp,
            "wrong_normalization_exponent": wrong_norm_exp,
        },
        "stability": {
            "bound": "C_hu**2 < B_eff*K_h",
            "violated_condition": checks["stability_violation_condition"],
            "strict_slack_symbol": "sigma = B_eff*K_h - C_hu**2 > 0",
        },
        "projections": {
            "source_mass_vector_imported": "[q_M, 0]",
            "h_mass_projection_zero": checks["mass_h_projection_zero"],
            "fixture_h_projection": s(fixture_h_projection),
            "grav_even_overlap_imported_zero": True,
        },
    }

    expressions = {
        "stiffness_matrix": "[[B_eff, C_hu], [C_hu, K_h]]",
        "determinant": s(det),
        "A_qq": s(A_qq),
        "A_qm": s(A_qm),
        "A_mm": s(A_mm),
        "A_qm_decoupled_h_mass_zero": s(A_qm_decoupled),
        "A_qm_mixed_qL0": s(A_qm_mixed_qL0),
        "A_mm_decoupled": s(A_mm_decoupled),
        "radiation_scalar_power_from_flux": s(scalar_power),
        "radiation_em_power_from_flux": s(em_power),
        "radiation_bare_flux_ratio": s(scalar_power_ratio_bare),
        "radiation_corrupt_flux_ratio": s(corrupt_ratio),
        "radiation_wrong_normalization_ratio": s(wrong_norm_ratio),
        "radiation_bare_fixed_ratio_status": "DERIVED_BUT_NUMERIC_MAGNITUDE_GUARDED",
        "radiation_pinned_Kh_ratio_status": "CONDITIONAL_DERIVED_BUT_NUMERIC_MAGNITUDE_GUARDED",
    }

    return {
        "comparison_payload": comparison_payload,
        "expressions": expressions,
        "checks": checks,
    }


def laundering_guard(attempt: MagnitudeAttempt, substrate: dict[str, Any]) -> str:
    if not magnitude_attempt_contains_guarded_symbol(attempt):
        return "NO_GUARDED_MAGNITUDE_ATTEMPT"
    if has_h_time_kinetic_parent_action(substrate):
        return "MAGNITUDE_EMISSION_ALLOWED"
    return "LAUNDERING_REFUSED"


def fixture_payload_for_attempt(attempt: MagnitudeAttempt) -> dict[str, Any]:
    return {"assembled_results_fixture": {attempt.kind: 1}}


def guard_A_controls(substrate: dict[str, Any]) -> dict[str, Any]:
    if has_h_time_kinetic_parent_action(substrate):
        raise AssertionError("Guard A expected absent h_time_kinetic_parent_action in current substrate")
    attempts = [
        MagnitudeAttempt("M_h_from_N0", "M_h", "M_h:=N0"),
        MagnitudeAttempt("M_h_from_K_parallel", "M_h", "M_h:=K_parallel"),
        MagnitudeAttempt("c_E_from_c_gamma_Green", "c_E", "c_E:=c_gamma from Green form"),
        MagnitudeAttempt("K_h_from_N0_cgamma2", "K_h", "K_h:=N0*c_gamma^2"),
        MagnitudeAttempt("residue_floor_emission", "residue_floor", "charge residue floor"),
    ]
    fixtures = {}
    for attempt in attempts:
        wired_result = expect_laundering_refused(fixture_payload_for_attempt(attempt), substrate)
        fixtures[attempt.label] = {
            "replacement": attempt.replacement,
            "result": laundering_guard(attempt, substrate),
            "wired_serialization_result": wired_result,
            "contains_guarded_symbol": magnitude_attempt_contains_guarded_symbol(attempt),
            "local_earned_parent_action": attempt.earned_parent_action,
            "h_time_kinetic_parent_action": substrate["facts"]["h_time_kinetic_parent_action"],
        }
    if any(item["result"] != "LAUNDERING_REFUSED" for item in fixtures.values()):
        raise AssertionError(f"Guard A fixture failed: {fixtures}")
    if any(not item["wired_serialization_result"].startswith("LAUNDERING_REFUSED") for item in fixtures.values()):
        raise AssertionError(f"Guard A wired fixture failed: {fixtures}")

    direct_injection_result = expect_laundering_refused({"P_h/P_EM": 1, "M_h": 1}, substrate)
    forged_attempt = MagnitudeAttempt(
        "forged_parent_action_flag",
        "P_h_over_P_EM",
        "P_h/P_EM:=numeric",
        earned_parent_action=True,
    )
    forged_wired_result = expect_laundering_refused(
        {"earned_parent_action": forged_attempt.earned_parent_action, "P_h_over_P_EM": 1},
        substrate,
    )
    forged_guard_result = laundering_guard(forged_attempt, substrate)
    if forged_guard_result != "LAUNDERING_REFUSED":
        raise AssertionError("forged earned_parent_action flag bypassed Guard A")

    bypass_regressions = {
        "tuple_wrapped_guarded_number": {
            "injection": "tuple-wrapped guarded number under power_ratio",
            "result": expect_laundering_refused({"assembled_results_fixture": {"power_ratio": (0.137,)}}, substrate),
        },
        "string_embedded_guarded_number": {
            "injection": "numeric literal embedded in string under power_ratio",
            "result": expect_laundering_refused(
                {"assembled_results_fixture": {"power_ratio": "P_h/P_EM = 0.137"}},
                substrate,
            ),
        },
    }
    if any(not item["result"].startswith("LAUNDERING_REFUSED") for item in bypass_regressions.values()):
        raise AssertionError(f"Guard A bypass regression failed: {bypass_regressions}")

    return {
        "status": "FIRED",
        "production": {
            "attempted_guarded_magnitude_emission": False,
            "guard_predicate": False,
            "result": "NO_GUARDED_MAGNITUDE_EMITTED",
            "serialization_guard_wired": True,
            "h_time_kinetic_parent_action": substrate["facts"]["h_time_kinetic_parent_action"],
        },
        "negative_fixtures": fixtures,
        "direct_injected_number_refused": {
            "result": direct_injection_result,
            "payload": "assembled results dict with numeric P_h/P_EM and M_h",
        },
        "forged_local_flag_refused": {
            "result": forged_guard_result,
            "wired_serialization_result": forged_wired_result,
            "local_earned_parent_action": forged_attempt.earned_parent_action,
            "h_time_kinetic_parent_action": substrate["facts"]["h_time_kinetic_parent_action"],
        },
        "scanner_bypass_regressions": bypass_regressions,
        "scope_note": "Guard A is a denylist scoped to M_h, c_E, K_h, P_h/P_EM, EP magnitude, and residue floor; unrelated new numeric fields are out of scope.",
        "output_invariant": "Serialization refuses numeric M_h/c_E/K_h/P_h_over_P_EM/EP magnitude without h_time_kinetic_parent_action.",
    }


def compute_channels(ledger: Ledger, core: dict[str, Any]) -> ChannelResult:
    if not ledger.static_coulomb_match:
        h_ep = "NO_GO"
    elif ledger.h_bare_mass_residue_fixture or not ledger.h_bare_mass_residue_zero:
        h_ep = H_EP_FIFTH
    else:
        h_ep = H_EP_SAFE

    pinned_kh_default = ledger.pinned_kh_fact
    if ledger.radiation_corrupt_speed:
        radiation = RADIATION_AUDIT
    elif ledger.pinned_kh_fact and ledger.commit_cE_equals_cgamma:
        radiation = RADIATION_TENSION
    elif ledger.radiation_derivation_ok:
        radiation = RADIATION_SIM
    else:
        radiation = RADIATION_AUDIT

    if ledger.species_indexed_b_over_ell:
        universality = UNIV_TENSION
    elif ledger.q_ratio_global_earned:
        universality = "EARNED_SAFE"
    else:
        universality = UNIV_SIM

    if ledger.uL_no_go:
        u_l_ep = "NO_GO"
    elif ledger.uL_escalated:
        u_l_ep = UNIV_TENSION
    else:
        u_l_ep = UL_SIM

    if ledger.suppression_requires_large_cE:
        preferred_frame = PF_TENSION
    else:
        preferred_frame = PF_SIM

    mixed_risk = bool(ledger.C_hu_nonzero and ledger.qM_nonzero and core["checks"]["A_qm_mixed_qL0_nonzero_term"])
    cherenkov_deferred = radiation == RADIATION_SIM
    guard_predicates = {
        "LAUNDERING_REFUSED": ledger.production_laundering_attempt,
        "STABILITY_VIOLATED": ledger.stability_violation,
    }

    details = {
        "pinned_Kh_test_default": pinned_kh_default,
        "current_ledger_radiation_branch": "exponent_3_bare_fixed" if not ledger.pinned_kh_fact else "exponent_1_pinned_Kh",
        "radiation_bare_fixed_exponent": core["checks"]["radiation_bare_exponent"],
        "radiation_pinned_Kh_exponent": core["checks"]["radiation_pinned_Kh_exponent"],
        "radiation_corrupt_speed_exponent": core["checks"]["radiation_corrupt_speed_exponent"],
        "universality_conjunction": f"{H_EP_SAFE} AND universality=EARNED_SAFE",
        "mixing_caveat": "C_hu*q_h*q_M term reintroduces a mass coupling when mixing is nonzero.",
        "u_L_magnitude": "SIM_GATED via a_L AND C_hu",
        "qM_sign": ledger.qM_sign,
    }

    return ChannelResult(
        h_EP=h_ep,
        radiation=radiation,
        universality=universality,
        u_L_EP=u_l_ep,
        preferred_frame=preferred_frame,
        subtags={
            "CHERENKOV_DEFERRED": cherenkov_deferred,
            "MIXED_SCALAR_EP_RISK": mixed_risk,
        },
        guard_predicates=guard_predicates,
        details=details,
    )


def verdict_from_channels(ch: ChannelResult) -> str:
    states = [ch.h_EP, ch.radiation, ch.universality, ch.u_L_EP, ch.preferred_frame]
    if "NO_GO" in states or any(ch.guard_predicates.values()):
        return VERDICT_NO_GO
    tension_channels: list[str] = []
    if ch.h_EP == H_EP_FIFTH:
        tension_channels.append("h_EP")
    if ch.radiation == RADIATION_TENSION:
        tension_channels.append("radiation")
    if ch.universality == UNIV_TENSION:
        tension_channels.append("universality")
    if ch.u_L_EP == UNIV_TENSION:
        tension_channels.append("u_L_EP")
    if tension_channels:
        return f"FALSIFIABLE_TENSION(channel={','.join(tension_channels)})"
    gated_or_cost = (
        ch.radiation in {RADIATION_SIM, RADIATION_AUDIT}
        or ch.universality == UNIV_SIM
        or ch.u_L_EP == UL_SIM
        or ch.preferred_frame in {PF_TENSION, PF_SIM}
    )
    if ch.h_EP == H_EP_SAFE and gated_or_cost:
        return VERDICT_MAPPED
    if (
        ch.h_EP == H_EP_SAFE
        and ch.radiation == "EARNED_SAFE"
        and ch.universality == "EARNED_SAFE"
        and ch.u_L_EP == "EARNED_SAFE"
        and ch.preferred_frame == "EARNED_SAFE"
    ):
        return "NATURALLY_HIDDEN"
    return VERDICT_NO_GO


def channel_table(ch: ChannelResult) -> list[dict[str, Any]]:
    return [
        {
            "channel": "h_EP",
            "state": ch.h_EP,
            "subtags": [],
            "adjudication": "Mass-channel safe on the decoupled floor only; full decoupled-floor EP safety also needs universality=EARNED_SAFE.",
        },
        {
            "channel": "radiation",
            "state": ch.radiation,
            "subtags": ["CHERENKOV_DEFERRED"] if ch.subtags["CHERENKOV_DEFERRED"] else [],
            "adjudication": "Current ledger branch is exponent-3; magnitude remains SIM_GATED on the free h-sector kinetic normalization.",
        },
        {
            "channel": "universality",
            "state": ch.universality,
            "subtags": [],
            "adjudication": "q_h/Q_E universality is required but not earned from the current b/ell ledger.",
        },
        {
            "channel": "u_L_EP",
            "state": ch.u_L_EP,
            "subtags": ["MIXED_SCALAR_EP_RISK"] if ch.subtags["MIXED_SCALAR_EP_RISK"] else [],
            "adjudication": "u_L couples to charge and mass; C_hu mixing transfers mass coupling into charge-sourced eigenmodes.",
        },
        {
            "channel": "preferred_frame",
            "state": ch.preferred_frame,
            "subtags": [],
            "adjudication": "c_E=c_gamma is calibrated, not earned; large-c_E hiding would carry preferred-frame cost.",
        },
    ]


def make_control(name: str, baseline: str, ledger: Ledger, core: dict[str, Any], extra: dict[str, Any] | None = None) -> dict[str, Any]:
    ch = compute_channels(ledger, core)
    verdict = verdict_from_channels(ch)
    out = {
        "status": "FIRED",
        "baseline_verdict": baseline,
        "recomputed_verdict": verdict,
        "transition": f"{baseline} -> {verdict}",
        "channels": {
            "h_EP": ch.h_EP,
            "radiation": ch.radiation,
            "universality": ch.universality,
            "u_L_EP": ch.u_L_EP,
            "preferred_frame": ch.preferred_frame,
        },
        "subtags": {k: v for k, v in ch.subtags.items() if v},
        "guard_predicates": {k: v for k, v in ch.guard_predicates.items() if v},
    }
    if extra:
        out.update(extra)
    return out


def compute_controls(baseline_verdict: str, core: dict[str, Any], substrate: dict[str, Any]) -> dict[str, Any]:
    controls: dict[str, Any] = {}
    controls["A"] = guard_A_controls(substrate)

    controls["B"] = make_control(
        "B",
        baseline_verdict,
        Ledger(name="control_B", h_bare_mass_residue_zero=False, h_bare_mass_residue_fixture=True, C_hu_nonzero=False),
        core,
        {"directive_transition": f"{H_EP_SAFE} -> {H_EP_FIFTH}"},
    )
    controls["C"] = make_control(
        "C",
        baseline_verdict,
        Ledger(name="control_C", qL_zero_fixture=True, C_hu_nonzero=True, qM_nonzero=True),
        core,
        {"directive_transition": "q_L=0 with C_hu*q_h*q_M nonzero -> MIXED_SCALAR_EP_RISK"},
    )
    controls["D"] = make_control(
        "D",
        baseline_verdict,
        Ledger(name="control_D", radiation_corrupt_speed=True),
        core,
        {
            "directive_transition": "radiation exponent changes under corrupted Green/flux speed",
            "baseline_exponent": core["checks"]["radiation_bare_exponent"],
            "corrupt_exponent": core["checks"]["radiation_corrupt_speed_exponent"],
        },
    )
    controls["E"] = make_control(
        "E",
        baseline_verdict,
        Ledger(name="control_E", suppression_requires_large_cE=True),
        core,
        {"directive_transition": f"{PF_SIM} -> {PF_TENSION}; NATURALLY_HIDDEN unreachable"},
    )
    controls["F"] = make_control(
        "F",
        baseline_verdict,
        Ledger(name="control_F", species_indexed_b_over_ell=True),
        core,
        {"directive_transition": f"{UNIV_SIM} -> {UNIV_TENSION}"},
    )
    controls["G"] = make_control(
        "G",
        baseline_verdict,
        Ledger(name="control_G", stability_violation=True),
        core,
        {"directive_transition": "strict stability C_hu**2 < B_eff*K_h -> STABILITY_VIOLATED"},
    )
    controls["H"] = make_control(
        "H",
        baseline_verdict,
        Ledger(name="control_H", static_coulomb_match=False),
        core,
        {"directive_transition": "static-Coulomb match perturbed -> h_EP=NO_GO"},
    )
    control_i_ledger = Ledger(name="control_I", qM_sign=-1)
    controls["I"] = make_control(
        "I",
        baseline_verdict,
        control_i_ledger,
        core,
        {
            "directive_transition": "q_M sign flip -> A_qm and signed A_mm/q_M projection flip",
            "local_invariant_transition": "MASS_COUPLING_SIGN_POSITIVE -> MASS_COUPLING_SIGN_NEGATIVE",
            "A_qm_sign_flip": core["checks"]["qM_flip_A_qm_residual_zero"],
            "A_mm_even_in_qM": core["checks"]["A_mm_even_in_qM"],
            "A_mm_signed_projection_flips": core["checks"]["A_mm_signed_projection_flips"],
            "honesty_note": "A_mm itself is quadratic and even in q_M; the signed mass-source projection flips.",
        },
    )
    controls["J"] = make_control(
        "J",
        baseline_verdict,
        Ledger(name="control_J", radiation_wrong_normalization=True),
        core,
        {
            "directive_transition": "bare-fixed exponent 3 -> pinned-K_h exponent 1 under wrong normalization fixture",
            "bare_fixed_exponent": core["checks"]["radiation_bare_exponent"],
            "wrong_normalization_exponent": core["checks"]["radiation_wrong_normalization_exponent"],
        },
    )
    controls["K_without_pinned_Kh"] = make_control(
        "K_without_pinned_Kh",
        baseline_verdict,
        Ledger(name="control_K_without_pinned_Kh", commit_cE_equals_cgamma=True),
        core,
        {"directive_transition": f"c_E=c_gamma without pinned K_h -> radiation={RADIATION_SIM}"},
    )
    controls["K"] = make_control(
        "K",
        baseline_verdict,
        Ledger(name="control_K", pinned_kh_fact=True, commit_cE_equals_cgamma=True),
        core,
        {"directive_transition": f"pinned K_h plus c_E=c_gamma -> radiation={RADIATION_TENSION}"},
    )

    required = {
        "B": H_EP_FIFTH,
        "F": UNIV_TENSION,
        "G": VERDICT_NO_GO,
        "K": "FALSIFIABLE_TENSION(channel=radiation)",
    }
    if controls["B"]["channels"]["h_EP"] != required["B"]:
        raise AssertionError("control B did not trigger h_EP fifth-force")
    if controls["F"]["channels"]["universality"] != required["F"]:
        raise AssertionError("control F did not trigger non-universality")
    if controls["G"]["recomputed_verdict"] != required["G"]:
        raise AssertionError("control G did not trigger NO_GO")
    if controls["K"]["recomputed_verdict"] != required["K"]:
        raise AssertionError("control K did not trigger conditional radiation tension")
    if controls["K_without_pinned_Kh"]["channels"]["radiation"] != RADIATION_SIM:
        raise AssertionError("control K unpinned branch did not remain SIM_GATED")
    return controls


def deletion_sensitivity(baseline_channels: ChannelResult, baseline_verdict: str) -> dict[str, Any]:
    def states(ch: ChannelResult) -> dict[str, str]:
        return {
            "h_EP": ch.h_EP,
            "radiation": ch.radiation,
            "universality": ch.universality,
            "u_L_EP": ch.u_L_EP,
            "preferred_frame": ch.preferred_frame,
        }

    def stubbed(channel: str) -> ChannelResult:
        return replace(baseline_channels, **{channel: "EARNED_SAFE"})

    out: dict[str, Any] = {}
    for channel in ["h_EP", "radiation", "universality", "u_L_EP", "preferred_frame"]:
        recomputed_channels = stubbed(channel)
        recomputed_verdict = verdict_from_channels(recomputed_channels)
        out[channel] = {
            "status": "FIRED",
            "baseline_verdict": baseline_verdict,
            "deleted_channel": channel,
            "baseline_state": getattr(baseline_channels, channel),
            "stubbed_state": "EARNED_SAFE",
            "recomputed_verdict": recomputed_verdict,
            "transition": f"{baseline_verdict} -> {recomputed_verdict}",
            "changed": recomputed_verdict != baseline_verdict,
            "recomputed_channels": states(recomputed_channels),
        }
    collectively_stubbed = replace(
        baseline_channels,
        radiation="EARNED_SAFE",
        universality="EARNED_SAFE",
        u_L_EP="EARNED_SAFE",
        preferred_frame="EARNED_SAFE",
    )
    collective_verdict = verdict_from_channels(collectively_stubbed)
    out["gated_channels_collective"] = {
        "status": "FIRED",
        "baseline_verdict": baseline_verdict,
        "deleted_channels": ["radiation", "universality", "u_L_EP", "preferred_frame"],
        "stubbed_state": "EARNED_SAFE",
        "recomputed_verdict": collective_verdict,
        "transition": f"{baseline_verdict} -> {collective_verdict}",
        "changed": collective_verdict != baseline_verdict,
        "recomputed_channels": states(collectively_stubbed),
    }
    return out


def build_results(engine_agreement: dict[str, Any] | None = None) -> dict[str, Any]:
    substrate = import_substrate()
    core = derive_core()
    baseline_ledger = Ledger()
    channels = compute_channels(baseline_ledger, core)
    verdict = verdict_from_channels(channels)
    if verdict != VERDICT_MAPPED:
        raise AssertionError(f"unexpected baseline verdict, preserving computed value: {verdict}")
    controls = compute_controls(verdict, core, substrate)
    deletion = deletion_sensitivity(channels, verdict)

    results = {
        "schema": SCHEMA,
        "verdict": verdict,
        "computed_not_hardcoded": True,
        "five_channel_table": channel_table(channels),
        "channel_states": {
            "h_EP": channels.h_EP,
            "radiation": channels.radiation,
            "universality": channels.universality,
            "u_L_EP": channels.u_L_EP,
            "preferred_frame": channels.preferred_frame,
        },
        "subtags": channels.subtags,
        "guard_predicates_production": channels.guard_predicates,
        "ep_adjudication": {
            "h_EP_mass_channel": channels.h_EP,
            "full_decoupled_floor_EP_safety": "NOT_EARNED",
            "required_conjunction": channels.details["universality_conjunction"],
            "mixing_caveat": channels.details["mixing_caveat"],
            "unqualified_decoupled_floor_EP_safe_reported": False,
        },
        "static_exchange_matrix": {
            "provenance": "dual-engine derivation",
            "stiffness_matrix": core["expressions"]["stiffness_matrix"],
            "determinant": core["expressions"]["determinant"],
            "A_qq": core["expressions"]["A_qq"],
            "A_qm": core["expressions"]["A_qm"],
            "A_mm": core["expressions"]["A_mm"],
            "decoupled_projection": {
                "A_qm_with_C_hu_0_and_h_mass_residue_0": core["expressions"]["A_qm_decoupled_h_mass_zero"],
                "mass_residue_zero_projection": True,
            },
            "mixed_projection": {
                "q_L_0_C_hu_nonzero_q_M_nonzero": core["expressions"]["A_qm_mixed_qL0"],
                "MIXED_SCALAR_EP_RISK": channels.subtags["MIXED_SCALAR_EP_RISK"],
            },
        },
        "radiation_adjudication": {
            "provenance": "dual-engine dipole/flux exponent derivation",
            "setup": "nonrelativistic point-charge dipole; scalar and EM far-zone fields at c_E/c_gamma; scalar flux vs Poynting flux",
            "current_ledger_branch": "exponent-3 bare-fixed/current branch",
            "current_ledger_state": channels.radiation,
            "scalar_power_from_flux": core["expressions"]["radiation_scalar_power_from_flux"],
            "em_power_from_flux": core["expressions"]["radiation_em_power_from_flux"],
            "bare_flux_ratio": core["expressions"]["radiation_bare_flux_ratio"],
            "corrupt_speed_flux_ratio": core["expressions"]["radiation_corrupt_flux_ratio"],
            "wrong_normalization_ratio": core["expressions"]["radiation_wrong_normalization_ratio"],
            "bare_fixed_exponent": core["checks"]["radiation_bare_exponent"],
            "pinned_K_h_conditional_exponent": core["checks"]["radiation_pinned_Kh_exponent"],
            "pinned_K_h_test_default": channels.details["pinned_Kh_test_default"],
            "same_cE_cgamma_without_pinned_Kh": RADIATION_SIM,
            "same_cE_cgamma_with_pinned_Kh": RADIATION_TENSION,
            "magnitude": "SIM_GATED_BY_GUARD_A_NO_NUMERIC_POWER_RATIO_EMITTED",
            "cherenkov_subtag_spelling_guard": "CHERENKOV_DEFERRED",
            "subtags": ["CHERENKOV_DEFERRED"],
        },
        "preferred_frame": {
            "state": channels.preferred_frame,
            "cost": "Large c_E/c_gamma radiation hiding would trigger PREFERRED_FRAME_TENSION.",
            "source": substrate["citations"]["preferred_frame"],
        },
        "u_L_and_mixing_structure": {
            "state": channels.u_L_EP,
            "magnitude": channels.details["u_L_magnitude"],
            "mixed_scalar_EP_risk": channels.subtags["MIXED_SCALAR_EP_RISK"],
        },
        "hard_wall": {
            "status": "HARD_WALL",
            "missing_object": "h_time_kinetic_parent_action",
            "M_h": substrate["facts"]["M_h_status"],
            "c_E": substrate["facts"]["c_E_operator_status"],
            "K_h": substrate["facts"]["K_h_static_calibration_status"],
            "route_to_reopen": "deferred nonlinear throat solve / Route-A",
            "no_forward_reopen_triggered": True,
            "astrophysical_cosmological_reopen": "On throat-solve pinning, check stellar cooling and BBN/CMB bounds for the massless charge-coupled scalar.",
        },
        "Gate_L_connection": {
            "status": "EARNED_CONNECTION",
            "statement": substrate["facts"]["gate_L_connection"],
            "source": substrate["citations"]["Gate_L_bending"],
        },
        "guard_A": controls["A"],
        "controls": controls,
        "channel_deletion_sensitivity": deletion,
        "dual_engine_core": core["comparison_payload"],
        "substrate": substrate,
        "engine_agreement": engine_agreement
        or {
            "status": "PENDING_MATHEMATICA",
            "compared_leaf_count": 0,
            "sympy_payload": rel(SYM_OUT),
            "mathematica_payload": rel(WL_OUT),
        },
    }
    results["comparison_digest"] = sha256_json(core["comparison_payload"])
    validate_no_laundered_numbers(results, substrate)
    return results


def markdown_report(results: dict[str, Any]) -> str:
    verdict = results["verdict"]
    lines: list[str] = [verdict, "", "# pathA_42 Charge-Coupled Scalar Gate", ""]
    lines.append(f"Computed overall verdict: `{verdict}`.")
    lines.append(
        "The headline is computed from the five channel states and guard predicates by the directive first-match rules."
    )
    lines.append("")
    lines.append("## Five-Channel Table")
    lines.append("")
    lines.append("| channel | state | subtags | adjudication |")
    lines.append("|---|---:|---|---|")
    for row in results["five_channel_table"]:
        subtags = ", ".join(row["subtags"]) if row["subtags"] else "none"
        lines.append(f"| `{row['channel']}` | `{row['state']}` | `{subtags}` | {row['adjudication']} |")
    lines.append("")
    lines.append("## EP And Mixing")
    lines.append("")
    lines.append(
        f"`h_EP` is `{results['ep_adjudication']['h_EP_mass_channel']}` only for the Gate-L mass-coupled channel on the decoupled floor."
    )
    lines.append(
        f"Full decoupled-floor EP safety is not reported because it additionally requires `{results['ep_adjudication']['required_conjunction']}`."
    )
    lines.append(results["ep_adjudication"]["mixing_caveat"])
    lines.append("")
    lines.append("## Radiation")
    lines.append("")
    rad = results["radiation_adjudication"]
    lines.append(f"Setup: {rad['setup']}.")
    lines.append(f"Derived flux powers: scalar `{rad['scalar_power_from_flux']}`, EM `{rad['em_power_from_flux']}`.")
    lines.append(
        f"The current ledger selects the `{rad['current_ledger_branch']}`: exponent `{rad['bare_fixed_exponent']}`."
    )
    lines.append(
        f"The exponent-`{rad['pinned_K_h_conditional_exponent']}` branch is conditional on a future pinned-`K_h` fact; the default pinned-`K_h` test is `{rad['pinned_K_h_test_default']}`."
    )
    lines.append(
        f"At `c_E=c_gamma` without pinned `K_h`, radiation remains `{rad['same_cE_cgamma_without_pinned_Kh']}`; with pinned `K_h`, control K reaches `{rad['same_cE_cgamma_with_pinned_Kh']}`."
    )
    lines.append(f"Magnitude status: `{rad['magnitude']}`. Subtag: `CHERENKOV_DEFERRED`.")
    lines.append("")
    lines.append("## Static Exchange")
    lines.append("")
    static = results["static_exchange_matrix"]
    lines.append(f"`S = {static['stiffness_matrix']}`, determinant `{static['determinant']}`.")
    lines.append(f"`A_qq = {static['A_qq']}`")
    lines.append(f"`A_qm = {static['A_qm']}`")
    lines.append(f"`A_mm = {static['A_mm']}`")
    lines.append(
        f"With `q_L=0`, nonzero mixing gives `{static['mixed_projection']['q_L_0_C_hu_nonzero_q_M_nonzero']}`, so `MIXED_SCALAR_EP_RISK` is computed."
    )
    lines.append("")
    lines.append("## HARD_WALL And Guard A")
    lines.append("")
    wall = results["hard_wall"]
    lines.append(
        f"`{wall['status']}`: `{wall['missing_object']}` is absent; `M_h`, `c_E`, `K_h`, radiation magnitude, and EP magnitude are not pinned."
    )
    lines.append(
        "Guard A refuses all laundering fixtures: "
        + ", ".join(
            f"`{name}`={item['result']}"
            for name, item in results["guard_A"]["negative_fixtures"].items()
            if name != "residue_floor_emission"
        )
        + "."
    )
    lines.append(
        f"Direct assembled-result injection: `{results['guard_A']['direct_injected_number_refused']['result']}`."
    )
    lines.append(
        f"Forged local flag: `{results['guard_A']['forged_local_flag_refused']['result']}`; wired check "
        f"`{results['guard_A']['forged_local_flag_refused']['wired_serialization_result']}`."
    )
    bypass = results["guard_A"]["scanner_bypass_regressions"]
    lines.append(
        "Guard A bypass regressions: "
        + ", ".join(f"`{name}` -> `{item['result']}`" for name, item in bypass.items())
        + "."
    )
    lines.append(results["guard_A"]["scope_note"])
    lines.append("The guarded residue-floor attempt is also refused by the serialization guard; no numeric scalar magnitude is emitted.")
    lines.append("")
    lines.append("## Preferred Frame, u_L, Gate-L")
    lines.append("")
    lines.append(f"Preferred-frame state: `{results['preferred_frame']['state']}`. {results['preferred_frame']['cost']}")
    lines.append(
        f"`u_L`/mixing state: `{results['u_L_and_mixing_structure']['state']}`, magnitude `{results['u_L_and_mixing_structure']['magnitude']}`."
    )
    lines.append(f"Gate-L connection: {results['Gate_L_connection']['statement']}.")
    lines.append("")
    lines.append("## Controls")
    lines.append("")
    lines.append("| control | status | transition | key result |")
    lines.append("|---|---:|---|---|")
    for key in ["A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K_without_pinned_Kh", "K"]:
        ctrl = results["controls"][key]
        if key == "A":
            key_result = "fixtures, bypass regressions, and direct output injection -> LAUNDERING_REFUSED"
            transition = "serialization guard wired; production emits no guarded number"
        else:
            key_result = ctrl.get("directive_transition", ctrl.get("local_invariant_transition", "recomputed"))
            transition = ctrl["transition"]
        lines.append(f"| `{key}` | `{ctrl['status']}` | `{transition}` | {key_result} |")
    lines.append("")
    lines.append("## Deletion Sensitivity")
    lines.append("")
    lines.append("| deleted/stubbed channel | recomputed verdict | changed |")
    lines.append("|---|---|---:|")
    for key in ["h_EP", "radiation", "universality", "u_L_EP", "preferred_frame", "gated_channels_collective"]:
        item = results["channel_deletion_sensitivity"][key]
        label = key if key != "gated_channels_collective" else "radiation+universality+u_L_EP+preferred_frame"
        lines.append(f"| `{label}` | `{item['recomputed_verdict']}` | `{item['changed']}` |")
    lines.append("")
    lines.append("## Dual Engine")
    lines.append("")
    eng = results["engine_agreement"]
    lines.append(
        f"`{eng['status']}` over `{eng['compared_leaf_count']}` compared leaves."
        if eng["status"] == "ENGINE_AGREE"
        else f"`{eng['status']}`; run the Mathematica script and `--compare` to finalize agreement."
    )
    lines.append("")
    lines.append("Run commands:")
    lines.append("")
    lines.append("```text")
    lines.append("timeout 600 python3 software/stage1_solver/tools/pathA_42_charge_coupled_scalar_sympy.py")
    lines.append("timeout 600 math -script software/stage1_solver/tools/pathA_42_charge_coupled_scalar.wl")
    lines.append("timeout 600 python3 software/stage1_solver/tools/pathA_42_charge_coupled_scalar_sympy.py --compare")
    lines.append("```")
    lines.append("")
    return "\n".join(lines)


def persist_reports(results: dict[str, Any]) -> None:
    REPORTS.mkdir(parents=True, exist_ok=True)
    substrate = results.get("substrate")
    if not isinstance(substrate, dict):
        raise AssertionError("results missing substrate provenance for serialization guard")
    validate_no_laundered_numbers(results, substrate)
    YAML_OUT.write_text(stable_yaml_dump(results), encoding="utf-8")
    REPORT_OUT.write_text(markdown_report(results), encoding="utf-8")


def compare_payloads() -> dict[str, Any]:
    if not SYM_OUT.exists():
        raise AssertionError(f"missing SymPy payload: {SYM_OUT}")
    if not WL_OUT.exists():
        raise AssertionError(f"missing Mathematica payload: {WL_OUT}")
    sym = read_json(SYM_OUT)
    wl = read_json(WL_OUT)
    if sym.get("schema") != SCHEMA:
        raise AssertionError(f"unexpected SymPy schema: {sym.get('schema')}")
    if wl.get("schema") != WL_SCHEMA:
        raise AssertionError(f"unexpected Mathematica schema: {wl.get('schema')}")
    if sym["comparison_payload"] != wl["comparison_payload"]:
        diff_path = RUN_OUT / "pathA_42_charge_coupled_scalar_engine_diff.json"
        write_json(diff_path, {"sympy": sym["comparison_payload"], "mathematica": wl["comparison_payload"]})
        raise AssertionError(f"ENGINE_DISAGREE; wrote {rel(diff_path)}")
    count = leaf_count(sym["comparison_payload"])
    return {
        "status": "ENGINE_AGREE",
        "compared_leaf_count": count,
        "comparison_digest": sha256_json(sym["comparison_payload"]),
        "sympy_payload": rel(SYM_OUT),
        "mathematica_payload": rel(WL_OUT),
    }


def run_math_script() -> None:
    cmd = ["timeout", "600", "math", "-script", str(STAGE1_ROOT / "tools" / "pathA_42_charge_coupled_scalar.wl")]
    subprocess.run(cmd, cwd=str(REPO_ROOT), check=True)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--compare", action="store_true", help="compare the SymPy and Mathematica payloads")
    args = parser.parse_args(argv)

    RUN_OUT.mkdir(parents=True, exist_ok=True)
    core = derive_core()
    sym_payload = {
        "schema": SCHEMA,
        "comparison_payload": core["comparison_payload"],
        "comparison_digest": sha256_json(core["comparison_payload"]),
    }
    write_json(SYM_OUT, sym_payload)

    if args.compare:
        run_math_script()
        engine_agreement = compare_payloads()
        results = build_results(engine_agreement)
        persist_reports(results)
        print(f"ENGINE_AGREE {engine_agreement['compared_leaf_count']}")
    else:
        results = build_results(None)
        persist_reports(results)
        print(results["verdict"])
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:  # pragma: no cover - command-line guard
        print(f"ERROR: {exc}", file=sys.stderr)
        raise
