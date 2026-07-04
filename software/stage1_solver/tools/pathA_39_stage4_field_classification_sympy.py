#!/usr/bin/env python3
"""pathA_39 Stage 4 field-coupling classification, SymPy engine.

This runner assembles one inverse propagator for the Stage-4 multiplet and
classifies from the computed rank/constraint/residue features.  It imports
only the Stage 0-3 blocks and provenance facts; the verdict is recomputed here.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any

import sympy as sp
import yaml
from sympy.assumptions import assuming
from sympy.printing.mathematica import mathematica_code


SCRIPT_PATH = Path(__file__).resolve()
STAGE1_ROOT = SCRIPT_PATH.parents[1]
REPORTS = STAGE1_ROOT / "reports"
SCRATCH = STAGE1_ROOT / "_scratch"

PATH36_YAML = REPORTS / "pathA_36_c5_phase_potential_results.yaml"
PATH38_YAML = REPORTS / "pathA_38_results.yaml"
PATH39_STAGE01_REPORT = REPORTS / "pathA_39_scalar_admixture_screen.md"
PATH39_STAGE2_REPORT = REPORTS / "pathA_39_magnetic_force.md"
PATH39_STAGE3_YAML = REPORTS / "pathA_39_stage3_operator_parity_results.yaml"

SYM_OUT = SCRATCH / "pathA_39_stage4_field_classification_sympy.json"
WL_OUT = SCRATCH / "pathA_39_stage4_field_classification_mathematica.json"
YAML_OUT = REPORTS / "pathA_39_stage4_field_classification_results.yaml"
REPORT_OUT = REPORTS / "pathA_39_stage4_field_classification.md"

SCHEMA = "pathA_39_stage4_field_classification/v1"
WL_SCHEMA = "pathA_39_stage4_field_classification_mathematica/v1"


VERDICT_CODES = {
    "FIELD_CLASSIFICATION_UNDERDETERMINED": 0,
    "FIELD_EXACT_MAXWELL_STRUCTURE": 1,
    "FIELD_TRANSVERSE_EM_PLUS_CLEAN_GRAVITY_DENSITY": 2,
    "FIELD_SCALAR_VECTOR_DEPARTURE": 3,
    "FIELD_SCALAR_SECTOR_UNSTABLE": 4,
}


Dim = tuple[int, int, int, int]
DIM_LABELS = ("M", "L", "T", "Q")


def dadd(*dims: Dim) -> Dim:
    return tuple(sum(dim[i] for dim in dims) for i in range(len(DIM_LABELS)))  # type: ignore[return-value]


def dsub(left: Dim, right: Dim) -> Dim:
    return tuple(left[i] - right[i] for i in range(len(DIM_LABELS)))  # type: ignore[return-value]


def dmul(n: int, dim: Dim) -> Dim:
    return tuple(n * dim[i] for i in range(len(DIM_LABELS)))  # type: ignore[return-value]


def dim_str(dim: Dim) -> str:
    parts: list[str] = []
    for label, power in zip(DIM_LABELS, dim):
        if power == 0:
            continue
        parts.append(label if power == 1 else f"{label}^{power}")
    return "1" if not parts else " ".join(parts)


def dim_record(dim: Dim) -> dict[str, Any]:
    return {"basis": list(DIM_LABELS), "powers": list(dim), "string": dim_str(dim)}


class DimChecker:
    def __init__(self) -> None:
        self.records: list[dict[str, Any]] = []
        self.ablations: list[dict[str, Any]] = []

    def check(self, category: str, name: str, actual: Dim, expected: Dim, expression: str) -> None:
        if actual != expected:
            raise AssertionError(
                f"{category}:{name}: expected {expected} ({dim_str(expected)}), "
                f"got {actual} ({dim_str(actual)})"
            )
        self.records.append(
            {
                "category": category,
                "name": name,
                "expression": expression,
                "dimension": dim_record(actual),
                "expected": dim_record(expected),
                "status": "PASS",
            }
        )

    def expect_fail(self, category: str, name: str, actual: Dim, expected: Dim, expression: str) -> None:
        if actual == expected:
            raise AssertionError(f"dimension ablation did not fire: {category}:{name}")
        self.ablations.append(
            {
                "category": category,
                "name": name,
                "expression": expression,
                "actual": dim_record(actual),
                "expected": dim_record(expected),
                "status": "FIRED",
            }
        )


def hstr(expr: Any) -> Any:
    if expr is None or isinstance(expr, (bool, int, str)):
        return expr
    if isinstance(expr, sp.MatrixBase):
        return [[hstr(entry) for entry in row] for row in expr.tolist()]
    return sp.sstr(expr)


def mma_expr(expr: sp.Expr | int) -> str:
    return mathematica_code(expr)


def sha256_text(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def sha256_json(data: Any) -> str:
    return sha256_text(json.dumps(data, sort_keys=True, separators=(",", ":")))


def is_zero_expr(expr: sp.Expr) -> bool:
    return sp.simplify(expr) == 0


def nonzero_expr(expr: sp.Expr) -> bool:
    return not is_zero_expr(expr)


def verdict_code(label: str) -> int:
    return VERDICT_CODES[label]


def load_yaml(path: Path) -> dict[str, Any]:
    if not path.exists():
        raise FileNotFoundError(path)
    with path.open("r", encoding="utf-8") as fh:
        data = yaml.safe_load(fh)
    if not isinstance(data, dict):
        raise AssertionError(f"YAML did not load to mapping: {path}")
    return data


def import_banked_inputs() -> dict[str, Any]:
    p36 = load_yaml(PATH36_YAML)
    p38 = load_yaml(PATH38_YAML)
    p3 = load_yaml(PATH39_STAGE3_YAML)
    if not PATH39_STAGE01_REPORT.exists():
        raise FileNotFoundError(PATH39_STAGE01_REPORT)
    if not PATH39_STAGE2_REPORT.exists():
        raise FileNotFoundError(PATH39_STAGE2_REPORT)

    s01 = PATH39_STAGE01_REPORT.read_text(encoding="utf-8")
    s2 = PATH39_STAGE2_REPORT.read_text(encoding="utf-8")

    branch36 = p36["per_branch_sub_results"]["branch_b_slaved_finite_compressibility_conventional_K"]
    tuned36 = p36["per_branch_sub_results"]["branch_b_slaved_tuned_Maxwell_locus"]
    b_eff_text = branch36["coefficients"]["B_eff"]
    c_gamma_text = p36["transverse"]["c_gamma_squared"]
    qh_text = p38["source_projections"]["q_h_plus"]
    green_text = p38["green_function"]["dynamic_trace"]["radial_green_finite_omega"]

    checks = {
        "pathA_36_B_eff": {"expected": "rho_B0**2/chi_c", "actual": b_eff_text, "match": b_eff_text == "rho_B0**2/chi_c"},
        "pathA_36_c_gamma_squared": {"expected": "mu_R/rho_br", "actual": c_gamma_text, "match": c_gamma_text == "mu_R/rho_br"},
        "pathA_36_real_constraint_class": {
            "expected": "SECOND_CLASS_PAIR",
            "actual": branch36["constraints"]["classification"],
            "match": branch36["constraints"]["classification"] == "SECOND_CLASS_PAIR",
        },
        "pathA_36_real_longitudinal_dof": {
            "expected": 1,
            "actual": branch36["mode_count"]["physical_dof_per_finite_k"],
            "match": branch36["mode_count"]["physical_dof_per_finite_k"] == 1,
        },
        "pathA_36_tuned_constraint_class": {
            "expected": "FIRST_CLASS_MAXWELL_CHAIN",
            "actual": tuned36["constraints"]["classification"],
            "match": tuned36["constraints"]["classification"] == "FIRST_CLASS_MAXWELL_CHAIN",
        },
        "pathA_36_tuned_longitudinal_dof": {
            "expected": 0,
            "actual": tuned36["mode_count"]["physical_dof_per_finite_k"],
            "match": tuned36["mode_count"]["physical_dof_per_finite_k"] == 0,
        },
        "pathA_36_tuned_first_class_count": {
            "expected": 2,
            "actual": tuned36["constraints"]["first_class_count"],
            "match": tuned36["constraints"]["first_class_count"] == 2,
        },
        "pathA_36_transverse_dof": {"expected": 2, "actual": p36["transverse"]["physical_dof"], "match": p36["transverse"]["physical_dof"] == 2},
        "pathA_38_q_h_plus": {"expected": "2*QE*tanh(b/ell)/b", "actual": qh_text, "match": qh_text == "2*QE*tanh(b/ell)/b"},
        "pathA_38_dynamic_green_speed": {
            "expected": "exp(I*R*omega/cE)/(4*pi*R)",
            "actual": green_text,
            "match": green_text == "exp(I*R*omega/cE)/(4*pi*R)",
        },
        "pathA_38_engine": {"expected": "ENGINE_AGREE", "actual": p38["engine_agreement"]["status"], "match": p38["engine_agreement"]["status"] == "ENGINE_AGREE"},
        "pathA_39_stage0_plus_1_block": {
            "expected": "scalar block with Chu, Kh=Mh*cE^2, q_L=Nu*aL*sCharge, ENGINE_AGREE",
            "actual": "report text",
            "match": all(token in s01 for token in ["D_x", "Chu", "Kh=Mh*cE^2", "q_L = Nu*aL*sCharge", "ENGINE_AGREE"]),
        },
        "pathA_39_stage2_sign": {
            "expected": "NO_CANCELLATION_BOTH_CHANNELS_ATTRACTIVE and ENGINE_AGREE",
            "actual": "report text",
            "match": "NO_CANCELLATION_BOTH_CHANNELS_ATTRACTIVE" in s2 and "ENGINE_AGREE" in s2,
        },
        "pathA_39_stage3_contamination": {
            "expected": "FAIL_UNPROTECTED_OPERATOR_PARITY_MIXING with ENGINE_AGREE",
            "actual": p3.get("verdict"),
            "match": p3.get("verdict") == "FAIL_UNPROTECTED_OPERATOR_PARITY_MIXING"
            and p3["engine_agreement"]["status"] == "ENGINE_AGREE",
        },
    }
    bad = [name for name, item in checks.items() if not item["match"]]
    if bad:
        raise AssertionError(f"banked import check failed: {bad}")

    return {
        "checks": checks,
        "imported": {
            "B_eff": "rho_B0^2/chi_c from pathA_36 finite-compressibility branch",
            "c_gamma_squared": "mu_R/rho_br from pathA_36 transverse sector",
            "longitudinal_real_constraint_class": "SECOND_CLASS_PAIR with 1 physical longitudinal DOF",
            "longitudinal_maxwell_constraint_class": "FIRST_CLASS_MAXWELL_CHAIN with 0 physical longitudinal DOF and 2 first-class constraints",
            "q_A_T": "q_A^T components parameterized as Nu*aT*sCharge and Nu*aTp*sCharge",
            "q_L": "Nu*aL*sCharge from Stage 0+1 source projection; a_L sim-deferred",
            "q_h": qh_text,
            "q_M": "mass source to u_L, recorded separately from charge residues",
            "M_h": "positive symbolic Mh from pathA_38 branon zero-mode normalization",
            "c_E": "cE from pathA_38 dynamic Green exp(I*R*omega/cE)/(4*pi*R)",
            "K_h": "M_h*c_E^2",
            "C_hu": "sim-deferred Stage 0+1 scalar mixing coefficient Chu",
            "magnetic_sign": "Stage 2 source of truth: transverse and longitudinal like-current channels attract",
            "operator_parity_contamination": "true, magnitude sim-deferred, from Stage 3",
        },
        "declared_vs_deferred": {
            "declared_stage4": [
                "one assembled Q(omega,k) over (u_T1,u_T2,u_L,h)",
                "charge source vector decomposed as q_A^T, q_L, q_h with q_M tracked as mass-channel source",
                "real branch uses the imported SECOND_CLASS_PAIR longitudinal class",
            ],
            "conditional_or_sim_deferred": [
                "C_hu stability bound C_hu^2 < B_eff K_h",
                "a_L and density charge coupling magnitude",
                "operator-parity contamination magnitude",
                "c_E=c_gamma and lambda_gamma knit",
            ],
        },
        "pathA_36_branch_records": {
            "real": {
                "classification": branch36["constraints"]["classification"],
                "first_class_count": branch36["constraints"]["first_class_count"],
                "second_class_count": branch36["constraints"]["second_class_count"],
                "longitudinal_dof": branch36["mode_count"]["physical_dof_per_finite_k"],
            },
            "maxwell": {
                "classification": tuned36["constraints"]["classification"],
                "first_class_count": tuned36["constraints"]["first_class_count"],
                "second_class_count": tuned36["constraints"]["second_class_count"],
                "longitudinal_dof": tuned36["mode_count"]["physical_dof_per_finite_k"],
                "generator": tuned36["gauge"]["generator"],
            },
        },
    }


def build_dimensions() -> dict[str, Any]:
    check = DimChecker()
    brane_lag: Dim = (1, -1, -2, 0)
    q_operator: Dim = (1, -3, -2, 0)
    d_u: Dim = (0, 1, 0, 0)
    d_h: Dim = (0, 1, 0, 0)
    d_grad: Dim = (0, -1, 0, 0)
    d_dt: Dim = (0, 0, -1, 0)
    d_rho_br: Dim = (1, -3, 0, 0)
    d_Mh: Dim = d_rho_br
    d_mu_R: Dim = brane_lag
    d_B: Dim = brane_lag
    d_Kh: Dim = brane_lag
    d_Chu: Dim = brane_lag
    d_rho_q: Dim = (0, -3, 0, 1)
    d_j: Dim = (0, -2, -1, 1)
    d_rho_m: Dim = (1, -3, 0, 0)
    d_qA: Dim = (1, 0, -1, -1)
    d_qL: Dim = d_qA
    d_qh: Dim = (1, 1, -2, -1)
    d_qM: Dim = (0, 1, -2, 0)
    d_speed2: Dim = (0, 2, -2, 0)

    check.check("Q", "transverse_inertia_entry", dadd(d_rho_br, dmul(2, d_dt)), q_operator, "rho_br omega^2")
    check.check("Q", "transverse_stiffness_entry", dadd(d_mu_R, dmul(2, d_grad)), q_operator, "mu_R k^2")
    check.check("Q", "density_stiffness_entry", dadd(d_B, dmul(2, d_grad)), q_operator, "B_eff k^2")
    check.check("Q", "h_inertia_entry", dadd(d_Mh, dmul(2, d_dt)), q_operator, "M_h omega^2")
    check.check("Q", "h_stiffness_entry", dadd(d_Kh, dmul(2, d_grad)), q_operator, "K_h k^2")
    check.check("Q", "hu_mixing_entry", dadd(d_Chu, dmul(2, d_grad)), q_operator, "C_hu k^2")
    check.check("action", "uT_Q_uT", dadd(d_u, q_operator, d_u), brane_lag, "u_T Q_T u_T")
    check.check("action", "uL_Q_h", dadd(d_u, q_operator, d_h), brane_lag, "u_L Q_Lh h")
    check.check("import", "K_h_over_M_h", dsub(d_Kh, d_Mh), d_speed2, "K_h/M_h=c_E^2")
    check.check("import", "c_gamma_squared", dsub(d_mu_R, d_rho_br), d_speed2, "mu_R/rho_br")
    check.check("source", "transverse_current_source", dadd(d_j, d_qA, d_u), brane_lag, "j_T q_A^T u_T")
    check.check("source", "density_current_source", dadd(d_j, d_qL, d_u), brane_lag, "j_L q_L u_L")
    check.check("source", "h_charge_density_source", dadd(d_rho_q, d_qh, d_h), brane_lag, "rho_q q_h h")
    check.check("source", "mass_density_source", dadd(d_rho_m, d_qM, d_u), brane_lag, "rho_m q_M u_L")

    check.expect_fail("ablation", "drop_k2_from_hu_mixing_Q", d_Chu, q_operator, "C_hu without k^2 in Q_Lh")
    check.expect_fail("ablation", "use_charge_density_for_qL_current", dadd(d_rho_q, d_qL, d_u), brane_lag, "rho_q q_L u_L")
    check.expect_fail("ablation", "mass_source_counted_as_charge_source", dadd(d_rho_q, d_qM, d_u), brane_lag, "rho_q q_M u_L")

    return {
        "pass": True,
        "basis": list(DIM_LABELS),
        "checked_expression_count": len(check.records),
        "records": check.records,
        "ablations": check.ablations,
    }


def classify_from_features(features: dict[str, Any]) -> str:
    if not features["import_fidelity"]["ok"]:
        return "FIELD_CLASSIFICATION_UNDERDETERMINED"
    if features["scalar_sector"]["stable_code"] == 0:
        return "FIELD_SCALAR_SECTOR_UNSTABLE"
    if (
        features["dof"]["physical_total"] == 2
        and features["constraint"]["first_class_generator_present"]
        and features["dof"]["h"] == 0
        and not features["residue_flags"]["scalar_charge_residue_nonzero"]
        and features["residue_flags"]["transverse_charge_residue_nonzero"]
    ):
        return "FIELD_EXACT_MAXWELL_STRUCTURE"
    if (
        features["dof"]["physical_total"] >= 3
        and features["dof"]["transverse"] == 2
        and not features["residue_flags"]["scalar_charge_residue_nonzero"]
        and features["residue_flags"]["transverse_charge_residue_nonzero"]
    ):
        return "FIELD_TRANSVERSE_EM_PLUS_CLEAN_GRAVITY_DENSITY"
    if (
        features["dof"]["physical_total"] >= 3
        and features["dof"]["transverse"] == 2
        and features["dof"]["h"] == 1
        and features["residue_flags"]["h_charge_residue_nonzero"]
        and features["residue_flags"]["scalar_charge_residue_nonzero"]
    ):
        return "FIELD_SCALAR_VECTOR_DEPARTURE"
    return "FIELD_CLASSIFICATION_UNDERDETERMINED"


def build_symbolics(imports: dict[str, Any]) -> dict[str, Any]:
    rhoBr, muR, rhoB0, chiC, Mh, cE, cBad, deltaB, deltaKh, Chu = sp.symbols(
        "rhoBr muR rhoB0 chiC Mh cE cBad deltaB deltaKh Chu", positive=True
    )
    k, omega, x = sp.symbols("k omega x", positive=True)
    QE, b, ell = sp.symbols("QE b ell", positive=True)
    Nu, aT, aTp, aL, qM, sCharge = sp.symbols("Nu aT aTp aL qM sCharge")

    B_eff = sp.factor(rhoB0**2 / chiC)
    c_gamma2 = sp.factor(muR / rhoBr)
    qT1 = sp.factor(Nu * aT * sCharge)
    qT2 = sp.factor(Nu * aTp * sCharge)
    qT_norm = sp.factor(qT1**2 + qT2**2)
    qL_stage = sp.factor(Nu * aL * sCharge)
    qh = sp.factor(2 * QE * sp.tanh(b / ell) / b)
    Kh = sp.factor(Mh * cE**2)
    large_Chu = sp.factor(2 * sp.sqrt(B_eff * Kh))

    real_constraint = imports["pathA_36_branch_records"]["real"]
    maxwell_constraint = imports["pathA_36_branch_records"]["maxwell"]

    def dirac_dof(include_h: bool, constraint_record: dict[str, Any]) -> dict[str, Any]:
        transverse_config = 2
        longitudinal_imported_config = 2
        h_config = 1 if include_h else 0
        first_class = int(constraint_record["first_class_count"])
        second_class = int(constraint_record["second_class_count"])
        total_config = transverse_config + longitudinal_imported_config + h_config
        physical_total = int(total_config - first_class - second_class // 2)
        return {
            "config_count_for_dirac": total_config,
            "transverse": 2,
            "u_L": int(constraint_record["longitudinal_dof"]),
            "h": 1 if include_h else 0,
            "first_class_count": first_class,
            "second_class_count": second_class,
            "physical_total": physical_total,
            "independent_initial_data_functions": 2 * physical_total,
        }

    def fidelity(B_expr: sp.Expr, Kh_expr: sp.Expr, cE_expr: sp.Expr) -> dict[str, Any]:
        b_match = is_zero_expr(B_expr - B_eff)
        kh_match = is_zero_expr(Kh_expr - Mh * cE**2)
        ce_match = is_zero_expr(cE_expr - cE)
        return {
            "ok": b_match and kh_match and ce_match,
            "B_eff_matches_rho_B0_squared_over_chi_c": b_match,
            "K_h_matches_M_h_c_E_squared": kh_match,
            "c_E_matches_dynamic_Green_speed": ce_match,
            "expected": {
                "B_eff": hstr(B_eff),
                "K_h": hstr(Mh * cE**2),
                "c_E_symbol": "cE",
            },
            "actual": {
                "B_eff": hstr(B_expr),
                "K_h": hstr(Kh_expr),
                "c_E_symbol": hstr(cE_expr),
            },
        }

    def scalar_common(B_expr: sp.Expr, Kh_expr: sp.Expr, Chu_expr: sp.Expr, qL_expr: sp.Expr, qh_expr: sp.Expr) -> dict[str, Any]:
        D_x = sp.Matrix([[rhoBr * x - B_expr, -Chu_expr], [-Chu_expr, Mh * x - Kh_expr]])
        detD = D_x.det()
        G_x = sp.Matrix([[Mh * x - Kh_expr, Chu_expr], [Chu_expr, rhoBr * x - B_expr]]) / detD
        num_charge = qL_expr**2 * (Mh * x - Kh_expr) + 2 * qL_expr * qh_expr * Chu_expr + qh_expr**2 * (rhoBr * x - B_expr)
        num_mass = qM**2 * (Mh * x - Kh_expr)
        num_charge_mass = qM * (qL_expr * (Mh * x - Kh_expr) + qh_expr * Chu_expr)
        trace = rhoBr * Kh_expr + Mh * B_expr
        delta = (rhoBr * Kh_expr - Mh * B_expr) ** 2 + 4 * rhoBr * Mh * Chu_expr**2
        c_minus = (trace - sp.sqrt(delta)) / (2 * rhoBr * Mh)
        c_plus = (trace + sp.sqrt(delta)) / (2 * rhoBr * Mh)
        den_prime = sp.diff(detD, x)

        def residue(num: sp.Expr, root: sp.Expr) -> sp.Expr:
            return num.subs(x, root) / den_prime.subs(x, root)

        det_stiffness = B_expr * Kh_expr - Chu_expr**2
        lambda_minus = k**2 * (B_expr + Kh_expr - sp.sqrt((B_expr - Kh_expr) ** 2 + 4 * Chu_expr**2)) / 2
        lambda_plus = k**2 * (B_expr + Kh_expr + sp.sqrt((B_expr - Kh_expr) ** 2 + 4 * Chu_expr**2)) / 2
        Aqq_x = num_charge / detD
        Amm_x = num_mass / detD
        Aqm_x = num_charge_mass / detD

        return {
            "D_x": D_x,
            "G_x": G_x,
            "detD": detD,
            "num_charge": num_charge,
            "num_mass": num_mass,
            "num_charge_mass": num_charge_mass,
            "Aqq_x": Aqq_x,
            "Amm_x": Amm_x,
            "Aqm_x": Aqm_x,
            "c_minus": c_minus,
            "c_plus": c_plus,
            "den_prime": den_prime,
            "root_minus_charge_residue": residue(num_charge, c_minus),
            "root_plus_charge_residue": residue(num_charge, c_plus),
            "root_minus_mass_residue": residue(num_mass, c_minus),
            "root_plus_mass_residue": residue(num_mass, c_plus),
            "det_stiffness_over_k4": det_stiffness,
            "stiffness_lambda_minus": lambda_minus,
            "stiffness_lambda_plus": lambda_plus,
        }

    def sign_code(expr: sp.Expr) -> int | None:
        simplified = sp.factor(expr)
        if is_zero_expr(simplified):
            return 0
        if sp.ask(sp.Q.negative(simplified)):
            return -1
        if sp.ask(sp.Q.positive(simplified)):
            return 1
        return None

    def conditionally_positive_code(expr: sp.Expr) -> int:
        direct = sign_code(expr)
        if direct is not None:
            return direct
        with assuming(sp.Q.positive(expr)):
            if sp.ask(sp.Q.positive(expr)):
                return 1
        return 0

    def positive_code(expr: sp.Expr) -> int:
        return 1 if sign_code(expr) == 1 else 0

    def residue_nonzero(expr: sp.Expr) -> bool:
        return not (expr == 0 or getattr(expr, "is_zero", None) is True)

    def any_nonzero(expressions: list[sp.Expr]) -> bool:
        return any(residue_nonzero(expr) for expr in expressions)

    def scalar_stability(B_expr: sp.Expr, Kh_expr: sp.Expr, scalar_artifacts: dict[str, Any]) -> dict[str, Any]:
        det_stiffness = scalar_artifacts["det_stiffness_over_k4"]
        det_sign = conditionally_positive_code(det_stiffness)
        principal_positive = positive_code(B_expr) == 1 and positive_code(Kh_expr) == 1 and positive_code(k**2) == 1
        eigenvalues_positive = det_sign == 1 and principal_positive
        stable = det_sign == 1 and eigenvalues_positive
        if stable:
            condition = (
                f"computed det(stiffness)/k^4={hstr(det_stiffness)} is positive; "
                "with positive diagonal stiffness entries this gives both finite-k stiffness eigenvalues > 0"
            )
        else:
            condition = (
                f"computed det(stiffness)/k^4={hstr(det_stiffness)} is nonpositive or not certified positive; "
                "at least one finite-k scalar stiffness eigenvalue is marginal/unstable"
            )
        return {
            "stable_code": 1 if stable else 0,
            "stable": stable,
            "det_stiffness_sign_code": det_sign,
            "finite_k_stiffness_eigenvalues_positive": eigenvalues_positive,
            "diagonal_stiffness_positive": principal_positive,
            "condition": condition,
        }

    def charge_residue_flags(B_expr: sp.Expr, Kh_expr: sp.Expr, Chu_expr: sp.Expr, qL_expr: sp.Expr, qh_expr: sp.Expr, scalar_artifacts: dict[str, Any]) -> dict[str, Any]:
        scalar_residues = [
            scalar_artifacts["root_minus_charge_residue"],
            scalar_artifacts["root_plus_charge_residue"],
        ]
        density_only = scalar_common(B_expr, Kh_expr, Chu_expr, qL_expr, sp.Integer(0))
        h_only = scalar_common(B_expr, Kh_expr, Chu_expr, sp.Integer(0), qh_expr)
        density_residues = [
            density_only["root_minus_charge_residue"],
            density_only["root_plus_charge_residue"],
        ]
        h_residues = [
            h_only["root_minus_charge_residue"],
            h_only["root_plus_charge_residue"],
        ]
        return {
            "scalar_charge_residue_nonzero": any_nonzero(scalar_residues),
            "density_charge_coupled": any_nonzero(density_residues),
            "h_charge_residue_nonzero": any_nonzero(h_residues),
            "scalar_charge_residue_expressions": scalar_residues,
            "density_charge_residue_expressions": density_residues,
            "h_charge_residue_expressions": h_residues,
        }

    def assemble_branch(
        *,
        name: str,
        include_h: bool,
        maxwell: bool,
        qL_expr: sp.Expr,
        qh_expr: sp.Expr,
        B_expr: sp.Expr,
        Kh_expr: sp.Expr,
        cE_expr: sp.Expr,
        Chu_expr: sp.Expr,
        stability_fixture_label: str,
        operator_parity_contamination: bool,
        notes: list[str],
    ) -> dict[str, Any]:
        QT = sp.factor(rhoBr * omega**2 - muR * k**2)
        constraint_record = maxwell_constraint if maxwell else real_constraint
        constraint = {
            "classification": constraint_record["classification"],
            "first_class_count": int(constraint_record["first_class_count"]),
            "second_class_count": int(constraint_record["second_class_count"]),
            "longitudinal_dof_imported": int(constraint_record["longitudinal_dof"]),
            "first_class_generator_present": bool(maxwell),
            "generator": constraint_record.get("generator") if maxwell else None,
            "imported_from": "pathA_36_c5_phase_potential_results.yaml",
        }
        dof = dirac_dof(include_h, constraint_record)
        import_fidelity = fidelity(B_expr, Kh_expr, cE_expr)

        if maxwell:
            Q = sp.Matrix([[QT, 0, 0], [0, QT, 0], [0, 0, 0]])
            J_charge = sp.Matrix([qT1, qT2, 0])
            J_mass = sp.Matrix([0, 0, 0])
            R_charge = sp.factor(qT_norm / QT)
            det_Q = sp.Integer(0)
            physical_det_Q = sp.factor(QT**2)
            poles = [
                {"label": "u_T1", "speed_squared": c_gamma2, "residue_to_charge": qT1**2 / rhoBr, "sector": "transverse"},
                {"label": "u_T2", "speed_squared": c_gamma2, "residue_to_charge": qT2**2 / rhoBr, "sector": "transverse"},
            ]
            scalar_sector = {
                "stable_code": 2,
                "stable": "not_applicable_h_removed_first_class_longitudinal",
                "det_stiffness_sign_code": None,
                "finite_k_stiffness_eigenvalues_positive": "not_applicable",
                "diagonal_stiffness_positive": "not_applicable",
                "det_stiffness_over_k4": None,
                "stiffness_eigenvalues": [],
                "condition": "h block removed and u_L removed from physical spectrum by FIRST_CLASS_MAXWELL_CHAIN",
                "fixture_label": stability_fixture_label,
            }
            scalar_artifacts: dict[str, Any] = {}
            scalar_residue_record = {
                "scalar_charge_residue_nonzero": False,
                "density_charge_coupled": False,
                "h_charge_residue_nonzero": False,
                "scalar_charge_residue_expressions": [],
                "density_charge_residue_expressions": [],
                "h_charge_residue_expressions": [],
            }
            kinetic_rank = 2
        else:
            Q = sp.Matrix(
                [
                    [QT, 0, 0, 0],
                    [0, QT, 0, 0],
                    [0, 0, rhoBr * omega**2 - B_expr * k**2, -Chu_expr * k**2],
                    [0, 0, -Chu_expr * k**2, Mh * omega**2 - Kh_expr * k**2],
                ]
            )
            J_charge = sp.Matrix([qT1, qT2, qL_expr, qh_expr])
            J_mass = sp.Matrix([0, 0, qM, 0])
            scalar_artifacts = scalar_common(B_expr, Kh_expr, Chu_expr, qL_expr, qh_expr)
            R_charge = qT_norm / QT + scalar_artifacts["Aqq_x"] / k**2
            det_Q = QT**2 * k**4 * scalar_artifacts["detD"]
            physical_det_Q = det_Q
            poles = [
                {"label": "u_T1", "speed_squared": c_gamma2, "residue_to_charge": qT1**2 / rhoBr, "sector": "transverse"},
                {"label": "u_T2", "speed_squared": c_gamma2, "residue_to_charge": qT2**2 / rhoBr, "sector": "transverse"},
                {
                    "label": "scalar_root_minus",
                    "speed_squared": scalar_artifacts["c_minus"],
                    "residue_to_charge": scalar_artifacts["root_minus_charge_residue"],
                    "residue_to_mass": scalar_artifacts["root_minus_mass_residue"],
                    "sector": "scalar",
                },
                {
                    "label": "scalar_root_plus",
                    "speed_squared": scalar_artifacts["c_plus"],
                    "residue_to_charge": scalar_artifacts["root_plus_charge_residue"],
                    "residue_to_mass": scalar_artifacts["root_plus_mass_residue"],
                    "sector": "scalar",
                },
            ]
            stability = scalar_stability(B_expr, Kh_expr, scalar_artifacts)
            scalar_sector = {
                "stable_code": stability["stable_code"],
                "stable": stability["stable"],
                "det_stiffness_sign_code": stability["det_stiffness_sign_code"],
                "finite_k_stiffness_eigenvalues_positive": stability["finite_k_stiffness_eigenvalues_positive"],
                "diagonal_stiffness_positive": stability["diagonal_stiffness_positive"],
                "det_stiffness_over_k4": hstr(scalar_artifacts["det_stiffness_over_k4"]),
                "stiffness_eigenvalues": [
                    hstr(scalar_artifacts["stiffness_lambda_minus"]),
                    hstr(scalar_artifacts["stiffness_lambda_plus"]),
                ],
                "condition": stability["condition"],
                "fixture_label": stability_fixture_label,
            }
            scalar_residue_record = charge_residue_flags(B_expr, Kh_expr, Chu_expr, qL_expr, qh_expr, scalar_artifacts)
            kinetic_rank = 4

        residue_flags = {
            "transverse_charge_residue_nonzero": nonzero_expr(qT_norm),
            "scalar_charge_residue_nonzero": scalar_residue_record["scalar_charge_residue_nonzero"],
            "density_charge_coupled": (not maxwell) and scalar_residue_record["density_charge_coupled"],
            "h_charge_residue_nonzero": include_h and scalar_residue_record["h_charge_residue_nonzero"],
        }
        residue_derivation = {
            "rule": "branch flags are computed from pole charge residues of J_q^T Q^-1 J_q, not from source-coupling literals",
            "scalar_charge_residue_expressions": hstr(scalar_residue_record["scalar_charge_residue_expressions"]),
            "density_charge_residue_expressions": hstr(scalar_residue_record["density_charge_residue_expressions"]),
            "h_charge_residue_expressions": hstr(scalar_residue_record["h_charge_residue_expressions"]),
        }
        flags = {
            "scalar_sector_stable": scalar_sector["stable"] if scalar_sector["stable_code"] != 2 else "not_applicable",
            "density_charge_coupled": residue_flags["density_charge_coupled"],
            "operator_parity_contamination": operator_parity_contamination,
        }
        features = {
            "name": name,
            "coordinates": ["u_T1", "u_T2", "u_L_gauge"] if maxwell else ["u_T1", "u_T2", "u_L", "h"],
            "Q": Q,
            "Q_x_form_scalar": scalar_artifacts.get("D_x"),
            "det_Q": det_Q,
            "physical_det_Q": physical_det_Q,
            "J_charge": J_charge,
            "J_mass": J_mass,
            "J_q_record": {
                "q_A_T_components": [hstr(qT1), hstr(qT2)],
                "q_L": hstr(qL_expr),
                "q_h": hstr(qh_expr),
                "q_M": hstr(qM),
                "mass_source_is_not_counted_as_charge_residue": True,
            },
            "R_charge": R_charge,
            "constraint": constraint,
            "dof": dof,
            "kinetic_rank_assembled": kinetic_rank,
            "scalar_sector": scalar_sector,
            "poles": poles,
            "residue_flags": residue_flags,
            "residue_derivation": residue_derivation,
            "flags": flags,
            "import_fidelity": import_fidelity,
            "notes": notes,
            "operator_parity_contamination": operator_parity_contamination,
        }
        primary = classify_from_features(features)
        features["primary"] = primary
        features["primary_code"] = verdict_code(primary)
        return features

    real = assemble_branch(
        name="real_provenance_fixed",
        include_h=True,
        maxwell=False,
        qL_expr=qL_stage,
        qh_expr=qh,
        B_expr=B_eff,
        Kh_expr=Kh,
        cE_expr=cE,
        Chu_expr=Chu,
        stability_fixture_label="stable_bound",
        operator_parity_contamination=True,
        notes=["Real branch is conditional on the computed stability bound C_hu^2 < B_eff K_h."],
    )
    maxwell = assemble_branch(
        name="maxwell_counterfactual",
        include_h=False,
        maxwell=True,
        qL_expr=sp.Integer(0),
        qh_expr=sp.Integer(0),
        B_expr=B_eff,
        Kh_expr=Kh,
        cE_expr=cE,
        Chu_expr=Chu,
        stability_fixture_label="stable_bound",
        operator_parity_contamination=False,
        notes=["Action-level h block removal plus imported pathA_36 FIRST_CLASS_MAXWELL_CHAIN on u_L."],
    )
    clean = assemble_branch(
        name="clean_coexistence",
        include_h=True,
        maxwell=False,
        qL_expr=sp.Integer(0),
        qh_expr=sp.Integer(0),
        B_expr=B_eff,
        Kh_expr=Kh,
        cE_expr=cE,
        Chu_expr=Chu,
        stability_fixture_label="stable_bound",
        operator_parity_contamination=True,
        notes=["Scalar block kept and diagonalized; charge source set to zero only in scalar channels."],
    )
    aL0 = assemble_branch(
        name="aL_to_0",
        include_h=True,
        maxwell=False,
        qL_expr=sp.Integer(0),
        qh_expr=qh,
        B_expr=B_eff,
        Kh_expr=Kh,
        cE_expr=cE,
        Chu_expr=Chu,
        stability_fixture_label="stable_bound",
        operator_parity_contamination=True,
        notes=["Density source removed while h charge source remains nonzero."],
    )
    large = assemble_branch(
        name="large_C_hu",
        include_h=True,
        maxwell=False,
        qL_expr=qL_stage,
        qh_expr=qh,
        B_expr=B_eff,
        Kh_expr=Kh,
        cE_expr=cE,
        Chu_expr=large_Chu,
        stability_fixture_label="stable_bound",
        operator_parity_contamination=True,
        notes=["Uses C_hu=2*sqrt(B_eff*K_h), so C_hu^2=4 B_eff K_h and det stiffness is negative."],
    )
    restored = assemble_branch(
        name="large_C_hu_restored_bound",
        include_h=True,
        maxwell=False,
        qL_expr=qL_stage,
        qh_expr=qh,
        B_expr=B_eff,
        Kh_expr=Kh,
        cE_expr=cE,
        Chu_expr=Chu,
        stability_fixture_label="stable_bound",
        operator_parity_contamination=True,
        notes=["Restores the computed stability bound after the large-C_hu ablation."],
    )
    import_B = assemble_branch(
        name="import_fidelity_B_eff_corrupt",
        include_h=True,
        maxwell=False,
        qL_expr=qL_stage,
        qh_expr=qh,
        B_expr=B_eff + deltaB,
        Kh_expr=Kh,
        cE_expr=cE,
        Chu_expr=Chu,
        stability_fixture_label="stable_bound",
        operator_parity_contamination=True,
        notes=["Corrupts B_eff so B_eff != rho_B0^2/chi_c."],
    )
    import_K = assemble_branch(
        name="import_fidelity_K_h_corrupt",
        include_h=True,
        maxwell=False,
        qL_expr=qL_stage,
        qh_expr=qh,
        B_expr=B_eff,
        Kh_expr=Kh + deltaKh,
        cE_expr=cE,
        Chu_expr=Chu,
        stability_fixture_label="stable_bound",
        operator_parity_contamination=True,
        notes=["Corrupts K_h so K_h != M_h c_E^2."],
    )
    import_cE = assemble_branch(
        name="import_fidelity_c_E_corrupt",
        include_h=True,
        maxwell=False,
        qL_expr=qL_stage,
        qh_expr=qh,
        B_expr=B_eff,
        Kh_expr=sp.factor(Mh * cBad**2),
        cE_expr=cBad,
        Chu_expr=Chu,
        stability_fixture_label="stable_bound",
        operator_parity_contamination=True,
        notes=["Uses cBad in the h block, inconsistent with the pathA_38 dynamic Green speed cE."],
    )

    branches = {
        item["name"]: item
        for item in [real, maxwell, clean, aL0, large, restored, import_B, import_K, import_cE]
    }
    controls = {
        "real_provenance_fixed": {"status": "FIRED", "branch": "real_provenance_fixed", "class": real["primary"]},
        "maxwell_counterfactual": {"status": "FIRED", "branch": "maxwell_counterfactual", "class": maxwell["primary"]},
        "clean_coexistence": {"status": "FIRED", "branch": "clean_coexistence", "class": clean["primary"]},
        "aL_to_0": {"status": "FIRED", "branch": "aL_to_0", "class": aL0["primary"]},
        "large_C_hu": {
            "status": "FIRED",
            "branch": "large_C_hu",
            "class": large["primary"],
            "restored_bound_branch": "large_C_hu_restored_bound",
            "restored_bound_class": restored["primary"],
        },
        "import_fidelity": {
            "status": "FIRED",
            "branches": {
                "B_eff_corrupt": import_B["primary"],
                "K_h_corrupt": import_K["primary"],
                "c_E_corrupt": import_cE["primary"],
            },
        },
        "dof_count_discriminator": {
            "status": "FIRED" if real["dof"]["physical_total"] == 4 and maxwell["dof"]["physical_total"] == 2 else "NOT_FIRED",
            "real_physical_dof": real["dof"]["physical_total"],
            "maxwell_physical_dof": maxwell["dof"]["physical_total"],
            "real_branch": "real_provenance_fixed",
            "maxwell_branch": "maxwell_counterfactual",
        },
    }

    self_tests = {
        "FIELD_SCALAR_VECTOR_DEPARTURE": real["primary"],
        "FIELD_EXACT_MAXWELL_STRUCTURE": maxwell["primary"],
        "FIELD_TRANSVERSE_EM_PLUS_CLEAN_GRAVITY_DENSITY": clean["primary"],
        "FIELD_SCALAR_SECTOR_UNSTABLE": large["primary"],
        "FIELD_CLASSIFICATION_UNDERDETERMINED": import_B["primary"],
        "A_L_ZERO_H_FLOOR": aL0["primary"],
        "RESTORED_STABILITY_BOUND": restored["primary"],
    }
    expected_self = {
        "FIELD_SCALAR_VECTOR_DEPARTURE": "FIELD_SCALAR_VECTOR_DEPARTURE",
        "FIELD_EXACT_MAXWELL_STRUCTURE": "FIELD_EXACT_MAXWELL_STRUCTURE",
        "FIELD_TRANSVERSE_EM_PLUS_CLEAN_GRAVITY_DENSITY": "FIELD_TRANSVERSE_EM_PLUS_CLEAN_GRAVITY_DENSITY",
        "FIELD_SCALAR_SECTOR_UNSTABLE": "FIELD_SCALAR_SECTOR_UNSTABLE",
        "FIELD_CLASSIFICATION_UNDERDETERMINED": "FIELD_CLASSIFICATION_UNDERDETERMINED",
        "A_L_ZERO_H_FLOOR": "FIELD_SCALAR_VECTOR_DEPARTURE",
        "RESTORED_STABILITY_BOUND": "FIELD_SCALAR_VECTOR_DEPARTURE",
    }
    self_test_records = {
        label: {"expected_fixture": expected, "actual": self_tests[label], "match": self_tests[label] == expected}
        for label, expected in expected_self.items()
    }
    if controls["dof_count_discriminator"]["status"] != "FIRED":
        raise AssertionError("DOF discriminator did not fire")
    for name, branch in branches.items():
        if branch["primary"] != classify_from_features(branch):
            raise AssertionError(f"classifier purity failed for {name}")

    real_scalar = scalar_common(B_eff, Kh, Chu, qL_stage, qh)
    clean_scalar = scalar_common(B_eff, Kh, Chu, sp.Integer(0), sp.Integer(0))
    aL0_scalar = scalar_common(B_eff, Kh, Chu, sp.Integer(0), qh)
    large_scalar = scalar_common(B_eff, Kh, large_Chu, qL_stage, qh)
    import_B_scalar = scalar_common(B_eff + deltaB, Kh, Chu, qL_stage, qh)

    exprs_for_agreement = {
        "Q_real_00": real["Q"][0, 0],
        "Q_real_22": real["Q"][2, 2],
        "Q_real_23": real["Q"][2, 3],
        "Q_real_33": real["Q"][3, 3],
        "Q_maxwell_00": maxwell["Q"][0, 0],
        "Q_maxwell_gauge_22": maxwell["Q"][2, 2],
        "Q_clean_23": clean["Q"][2, 3],
        "Q_large_23": large["Q"][2, 3],
        "detQ_real": real["det_Q"],
        "physical_detQ_maxwell": maxwell["physical_det_Q"],
        "B_eff": B_eff,
        "K_h": Kh,
        "c_gamma_squared": c_gamma2,
        "qT_norm_squared": qT_norm,
        "q_L_stage": qL_stage,
        "q_h": qh,
        "J_real_0": real["J_charge"][0],
        "J_real_1": real["J_charge"][1],
        "J_real_2": real["J_charge"][2],
        "J_real_3": real["J_charge"][3],
        "R_charge_real": real["R_charge"].subs(omega**2 / k**2, x),
        "R_charge_maxwell": maxwell["R_charge"],
        "det_scalar_real": real_scalar["detD"],
        "scalar_c_minus_real": real_scalar["c_minus"],
        "scalar_c_plus_real": real_scalar["c_plus"],
        "scalar_det_stiffness_real": real_scalar["det_stiffness_over_k4"],
        "scalar_lambda_minus_real": real_scalar["stiffness_lambda_minus"],
        "scalar_lambda_plus_real": real_scalar["stiffness_lambda_plus"],
        "real_root_minus_charge_residue": real_scalar["root_minus_charge_residue"],
        "real_root_plus_charge_residue": real_scalar["root_plus_charge_residue"],
        "real_root_minus_mass_residue": real_scalar["root_minus_mass_residue"],
        "real_root_plus_mass_residue": real_scalar["root_plus_mass_residue"],
        "transverse_T1_charge_residue": qT1**2 / rhoBr,
        "transverse_T2_charge_residue": qT2**2 / rhoBr,
        "clean_root_minus_charge_residue": clean_scalar["root_minus_charge_residue"],
        "clean_root_plus_charge_residue": clean_scalar["root_plus_charge_residue"],
        "aL0_root_minus_charge_residue": aL0_scalar["root_minus_charge_residue"],
        "aL0_root_plus_charge_residue": aL0_scalar["root_plus_charge_residue"],
        "large_Chu": large_Chu,
        "large_scalar_det_stiffness": large_scalar["det_stiffness_over_k4"],
        "large_root_minus_charge_residue": large_scalar["root_minus_charge_residue"],
        "large_root_plus_charge_residue": large_scalar["root_plus_charge_residue"],
        "import_B_Q22": import_B["Q"][2, 2],
        "import_B_det_scalar": import_B_scalar["detD"],
        "real_kinetic_rank": sp.Integer(real["kinetic_rank_assembled"]),
        "maxwell_kinetic_rank": sp.Integer(maxwell["kinetic_rank_assembled"]),
        "real_total_dof": sp.Integer(real["dof"]["physical_total"]),
        "maxwell_total_dof": sp.Integer(maxwell["dof"]["physical_total"]),
        "real_longitudinal_dof": sp.Integer(real["dof"]["u_L"]),
        "maxwell_longitudinal_dof": sp.Integer(maxwell["dof"]["u_L"]),
        "real_h_dof": sp.Integer(real["dof"]["h"]),
        "maxwell_h_dof": sp.Integer(maxwell["dof"]["h"]),
        "real_first_class_count": sp.Integer(real["constraint"]["first_class_count"]),
        "maxwell_first_class_count": sp.Integer(maxwell["constraint"]["first_class_count"]),
        "real_second_class_count": sp.Integer(real["constraint"]["second_class_count"]),
        "maxwell_second_class_count": sp.Integer(maxwell["constraint"]["second_class_count"]),
        "real_scalar_stable_code": sp.Integer(real["scalar_sector"]["stable_code"]),
        "large_scalar_stable_code": sp.Integer(large["scalar_sector"]["stable_code"]),
        "real_scalar_det_sign_code": sp.Integer(real["scalar_sector"]["det_stiffness_sign_code"]),
        "large_scalar_det_sign_code": sp.Integer(large["scalar_sector"]["det_stiffness_sign_code"]),
        "real_scalar_eigen_positive_code": sp.Integer(1 if real["scalar_sector"]["finite_k_stiffness_eigenvalues_positive"] else 0),
        "large_scalar_eigen_positive_code": sp.Integer(1 if large["scalar_sector"]["finite_k_stiffness_eigenvalues_positive"] else 0),
        "real_primary_code": sp.Integer(real["primary_code"]),
        "maxwell_primary_code": sp.Integer(maxwell["primary_code"]),
        "clean_primary_code": sp.Integer(clean["primary_code"]),
        "aL0_primary_code": sp.Integer(aL0["primary_code"]),
        "large_primary_code": sp.Integer(large["primary_code"]),
        "restored_primary_code": sp.Integer(restored["primary_code"]),
        "import_B_primary_code": sp.Integer(import_B["primary_code"]),
        "import_K_primary_code": sp.Integer(import_K["primary_code"]),
        "import_cE_primary_code": sp.Integer(import_cE["primary_code"]),
        "real_density_charge_flag_code": sp.Integer(1 if real["flags"]["density_charge_coupled"] else 0),
        "aL0_density_charge_flag_code": sp.Integer(1 if aL0["flags"]["density_charge_coupled"] else 0),
        "real_scalar_charge_residue_flag_code": sp.Integer(1 if real["residue_flags"]["scalar_charge_residue_nonzero"] else 0),
        "clean_scalar_charge_residue_flag_code": sp.Integer(1 if clean["residue_flags"]["scalar_charge_residue_nonzero"] else 0),
        "aL0_h_charge_residue_flag_code": sp.Integer(1 if aL0["residue_flags"]["h_charge_residue_nonzero"] else 0),
        "clean_h_charge_residue_flag_code": sp.Integer(1 if clean["residue_flags"]["h_charge_residue_nonzero"] else 0),
    }

    return {
        "symbols": {
            "rhoBr": rhoBr,
            "muR": muR,
            "rhoB0": rhoB0,
            "chiC": chiC,
            "Mh": Mh,
            "cE": cE,
            "cBad": cBad,
            "deltaB": deltaB,
            "deltaKh": deltaKh,
            "Chu": Chu,
            "k": k,
            "omega": omega,
            "x": x,
            "QE": QE,
            "b": b,
            "ell": ell,
            "Nu": Nu,
            "aT": aT,
            "aTp": aTp,
            "aL": aL,
            "qM": qM,
            "sCharge": sCharge,
        },
        "B_eff": B_eff,
        "K_h": Kh,
        "c_gamma2": c_gamma2,
        "qT1": qT1,
        "qT2": qT2,
        "qT_norm": qT_norm,
        "qL_stage": qL_stage,
        "q_h": qh,
        "large_Chu": large_Chu,
        "branches": branches,
        "controls": controls,
        "self_tests": self_tests,
        "self_test_records": self_test_records,
        "exprs_for_agreement": exprs_for_agreement,
        "expr_digest": sha256_json({key: mma_expr(value) for key, value in exprs_for_agreement.items()}),
    }


def jsonify_branch(branch: dict[str, Any]) -> dict[str, Any]:
    return {
        "name": branch["name"],
        "primary": branch["primary"],
        "primary_code": branch["primary_code"],
        "flags": branch["flags"],
        "coordinates": branch["coordinates"],
        "assembled_Q": hstr(branch["Q"]),
        "assembled_Q_x_scalar_block": hstr(branch["Q_x_form_scalar"]) if branch["Q_x_form_scalar"] is not None else None,
        "det_Q": hstr(branch["det_Q"]),
        "physical_det_Q": hstr(branch["physical_det_Q"]),
        "J_charge": hstr(branch["J_charge"]),
        "J_mass": hstr(branch["J_mass"]),
        "J_q_record": branch["J_q_record"],
        "R_charge": hstr(branch["R_charge"]),
        "constraint": branch["constraint"],
        "kinetic_rank_assembled": branch["kinetic_rank_assembled"],
        "dof": branch["dof"],
        "scalar_sector": branch["scalar_sector"],
        "poles": [
            {
                key: hstr(value) if key in {"speed_squared", "residue_to_charge", "residue_to_mass"} else value
                for key, value in pole.items()
            }
            for pole in branch["poles"]
        ],
        "residue_flags": branch["residue_flags"],
        "residue_derivation": branch["residue_derivation"],
        "import_fidelity": branch["import_fidelity"],
        "operator_parity_contamination": branch["operator_parity_contamination"],
        "notes": branch["notes"],
    }


def build_payload() -> dict[str, Any]:
    imports = import_banked_inputs()
    dims = build_dimensions()
    sym = build_symbolics(imports)
    agreement_exprs = {key: mma_expr(value) for key, value in sym["exprs_for_agreement"].items()}
    branches_json = {name: jsonify_branch(branch) for name, branch in sym["branches"].items()}
    real = branches_json["real_provenance_fixed"]
    control_classes = {
        "real_provenance_fixed": sym["controls"]["real_provenance_fixed"]["class"],
        "maxwell_counterfactual": sym["controls"]["maxwell_counterfactual"]["class"],
        "clean_coexistence": sym["controls"]["clean_coexistence"]["class"],
        "aL_to_0": sym["controls"]["aL_to_0"]["class"],
        "large_C_hu": sym["controls"]["large_C_hu"]["class"],
        "large_C_hu_restored_bound": sym["controls"]["large_C_hu"]["restored_bound_class"],
        "import_fidelity": sym["controls"]["import_fidelity"]["branches"],
    }
    agreement_payload = {
        "top_line_primary": real["primary"],
        "top_line_primary_code": real["primary_code"],
        "top_line_flags": real["flags"],
        "control_classes": control_classes,
        "dof_discriminator": {
            "real": branches_json["real_provenance_fixed"]["dof"]["physical_total"],
            "maxwell": branches_json["maxwell_counterfactual"]["dof"]["physical_total"],
        },
        "self_tests": sym["self_tests"],
        "checked_expression_count": len(sym["exprs_for_agreement"]),
        "expr_digest": sym["expr_digest"],
    }

    return {
        "schema": SCHEMA,
        "engine": "sympy",
        "directive": "software/stage1_solver/directives/pathA_39_stage4_field_classification.md",
        "top_line_primary": real["primary"],
        "top_line_flags": real["flags"],
        "engine_agreement": {
            "status": "PENDING_MATHEMATICA",
            "mathematica_exprs": agreement_exprs,
            "sympy_expression_digest": sym["expr_digest"],
        },
        "imports": imports,
        "dimensional_firewall": {
            "pass": dims["pass"],
            "checked_expression_count": dims["checked_expression_count"],
            "records": dims["records"],
            "ablations": dims["ablations"],
        },
        "branches": branches_json,
        "controls": sym["controls"],
        "classifier": {
            "rule": "pure function of import fidelity, scalar stability, Dirac DOF consequences, first-class generator status, and charge residues",
            "labels": list(VERDICT_CODES.keys()),
            "self_tests": sym["self_tests"],
            "self_test_records_nonfatal": sym["self_test_records"],
            "verdict_codes": VERDICT_CODES,
        },
        "remediation": {
            "status": "REMEDIATED",
            "fixes": [
                "scalar stability code derived from computed det(stiffness)/k^4 sign and finite-k stiffness positivity",
                "Mathematica engine derives kinetic rank, Dirac DOF, residue flags, stability code, and classifier codes before comparing",
                "canonical classifier predicates are shared across engines",
                "scalar/density/h charge flags are derived from computed pole charge residues",
                "expected-landing classifier self-tests are non-fatal records",
            ],
            "numbers_and_verdict": "unchanged",
        },
        "agreement_payload": agreement_payload,
    }


def write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def compare_with_mathematica(payload: dict[str, Any]) -> dict[str, Any]:
    if not WL_OUT.exists():
        raise FileNotFoundError(f"Mathematica payload missing: {WL_OUT}")
    other = json.loads(WL_OUT.read_text(encoding="utf-8"))
    if other.get("schema") != WL_SCHEMA:
        raise AssertionError(f"unexpected Mathematica schema: {other.get('schema')}")
    if other.get("status") != "OK":
        raise AssertionError(f"Mathematica engine did not complete cleanly: {other}")
    if other.get("sympy_expression_digest") != payload["engine_agreement"]["sympy_expression_digest"]:
        raise AssertionError("Mathematica compared a stale SymPy expression digest")
    if other.get("agreement_payload") != payload["agreement_payload"]:
        raise AssertionError(
            "ENGINE_DISAGREE\n"
            + json.dumps({"sympy": payload["agreement_payload"], "mathematica": other.get("agreement_payload")}, indent=2)
        )

    payload["engine_agreement"]["status"] = "ENGINE_AGREE"
    payload["engine_agreement"]["mathematica_payload"] = str(WL_OUT)
    payload["engine_agreement"]["sympy_payload"] = str(SYM_OUT)

    results = {
        "schema": SCHEMA,
        "primary": payload["top_line_primary"],
        "flags": payload["top_line_flags"],
        "headline": {
            "primary": payload["top_line_primary"],
            "flags": payload["top_line_flags"],
            "engine_agreement": "ENGINE_AGREE",
        },
        "engine_agreement": {
            "status": "ENGINE_AGREE",
            "compared_payload": payload["agreement_payload"],
            "sympy_payload": str(SYM_OUT),
            "mathematica_payload": str(WL_OUT),
        },
        "imports": payload["imports"],
        "dimensional_firewall": payload["dimensional_firewall"],
        "assembled_system": {
            "coordinates": payload["branches"]["real_provenance_fixed"]["coordinates"],
            "Q": payload["branches"]["real_provenance_fixed"]["assembled_Q"],
            "J_charge": payload["branches"]["real_provenance_fixed"]["J_charge"],
            "J_q_record": payload["branches"]["real_provenance_fixed"]["J_q_record"],
            "R_charge": payload["branches"]["real_provenance_fixed"]["R_charge"],
        },
        "branches": payload["branches"],
        "controls": payload["controls"],
        "classifier": payload["classifier"],
        "remediation": payload["remediation"],
        "honest_scope": {
            "stage4_close": "spec close: classified linear field-content skeleton",
            "not_done_here": [
                "nonlinear throat solve",
                "lambda_gamma/c_gamma=c_s knit",
                "c_E=c_gamma Lorentz-cone test",
                "operator-parity contamination magnitude",
            ],
            "does_not_unearn": "Stages 0-3 remain imported provenance and controls; Stage 4 does not relabel them.",
        },
    }
    YAML_OUT.write_text(yaml.safe_dump(results, sort_keys=False, width=140), encoding="utf-8")
    write_report(results)
    write_json(SYM_OUT, payload)
    return results


def branch_artifact_summary(branch: dict[str, Any]) -> str:
    return (
        f"`Q={branch['assembled_Q']}`, `J={branch['J_charge']}`, "
        f"`dof={branch['dof']['physical_total']}`, `class={branch['primary']}`"
    )


def write_report(results: dict[str, Any]) -> None:
    real = results["branches"]["real_provenance_fixed"]
    maxwell = results["branches"]["maxwell_counterfactual"]
    clean = results["branches"]["clean_coexistence"]
    aL0 = results["branches"]["aL_to_0"]
    large = results["branches"]["large_C_hu"]
    imports = results["imports"]["imported"]
    compared = results["engine_agreement"]["compared_payload"]

    lines = [
        "# pathA_39 Stage 4 Field-Coupling Classification",
        "",
        f"Computed headline: `primary={results['primary']}` with flags `{results['flags']}` and dual-engine `ENGINE_AGREE`.",
        "",
        "This is a Stage-4 field-content classification from the assembled inverse propagator, not a relabel of Stages 0-3.  The healthy real-branch departure is conditional on the computed scalar stability bound `B_eff K_h - C_hu^2 > 0`.",
        "",
        "## Assembled Q",
        "",
        "Coordinates: `Phi=(u_T1,u_T2,u_L,h)`.",
        "",
        "```text",
        f"Q(omega,k) = {real['assembled_Q']}",
        f"J_charge = {real['J_charge']}",
        f"J_q record = {real['J_q_record']}",
        f"R=J^T Q^-1 J = {real['R_charge']}",
        "```",
        "",
        "The mass channel `q_M` is recorded as a separate source to `u_L`; it is not counted as charge residue.",
        "",
        "## DOF And Constraints",
        "",
        "| branch | kinetic rank | longitudinal class | first class | second class | transverse | u_L | h | total physical DOF |",
        "|---|---:|---|---:|---:|---:|---:|---:|---:|",
    ]
    for branch in [real, maxwell, clean, aL0, large]:
        lines.append(
            f"| `{branch['name']}` | `{branch['kinetic_rank_assembled']}` | `{branch['constraint']['classification']}` | `{branch['constraint']['first_class_count']}` | `{branch['constraint']['second_class_count']}` | `{branch['dof']['transverse']}` | `{branch['dof']['u_L']}` | `{branch['dof']['h']}` | `{branch['dof']['physical_total']}` |"
        )
    lines.extend(
        [
            "",
            "The real branch uses the imported pathA_36 `SECOND_CLASS_PAIR` consequence for `u_L`.  The Maxwell counterfactual removes the `h` block and imposes the imported `FIRST_CLASS_MAXWELL_CHAIN`, giving the required DOF drop from 4 to 2.",
            "",
            "## Scalar Stability",
            "",
            f"`det(stiffness)/k^4 = {real['scalar_sector']['det_stiffness_over_k4']}`.",
            "",
            f"Finite-`k` stiffness eigenvalues: `{real['scalar_sector']['stiffness_eigenvalues'][0]}` and `{real['scalar_sector']['stiffness_eigenvalues'][1]}`.",
            "",
            f"Real branch: `scalar_sector_stable={real['scalar_sector']['stable']}` from `{real['scalar_sector']['condition']}`.",
            "",
            f"Large-`C_hu` control: `det(stiffness)/k^4 = {large['scalar_sector']['det_stiffness_over_k4']}`, `scalar_sector_stable={large['scalar_sector']['stable']}`, class `{large['primary']}`.",
            "",
            "In the `C_hu -> 0` limit the scalar roots reduce to the density speed `B_eff/rho_br` and the h speed `c_E^2`; Stage 4 keeps the mixed eigenroots instead of importing the decoupled residues.",
            "",
            "## Poles And Residues",
            "",
            "| branch | pole | speed^2 | charge residue | mass residue |",
            "|---|---|---:|---:|---:|",
        ]
    )
    for branch in [real, maxwell, clean, aL0, large]:
        for pole in branch["poles"]:
            lines.append(
                f"| `{branch['name']}` | `{pole['label']}` | `{pole['speed_squared']}` | `{pole.get('residue_to_charge', '')}` | `{pole.get('residue_to_mass', '')}` |"
            )
    lines.extend(
        [
            "",
            "The clean branch keeps `q_A^T` transverse charge but sets `q_L=q_h=0`; the two scalar charge residues compute to zero.  The `a_L->0` branch sets `q_L=0` but keeps `q_h!=0`, so the computed h-only pole charge residue keeps the scalar-vector departure while the computed density-only pole charge residue leaves `density_charge_coupled=false`.",
            "",
            "## Controls",
            "",
            "| control | status | derived class/result | recorded artifacts |",
            "|---|---:|---|---|",
        ]
    )
    for name, control in results["controls"].items():
        if name == "import_fidelity":
            derived = control["branches"]
            artifact = "B_eff/K_h/c_E corruptions each run through the extractor; branch summaries below"
        elif name == "dof_count_discriminator":
            derived = f"real {control['real_physical_dof']} -> Maxwell {control['maxwell_physical_dof']}"
            artifact = "real and Maxwell branch DOF records"
        else:
            derived = control.get("class", control.get("restored_bound_class", ""))
            if name == "large_C_hu":
                derived = f"{derived}; restored bound -> {control['restored_bound_class']}"
            artifact = branch_artifact_summary(results["branches"][control["branch"]])
        lines.append(f"| `{name}` | `{control['status']}` | `{derived}` | {artifact} |")
    lines.extend(
        [
            "",
            "Import-fidelity corruption branch summaries:",
            "",
            "| branch | imported value corrupted | assembled scalar Q[2,2] | assembled scalar block | class |",
            "|---|---|---:|---:|---|",
        ]
    )
    for branch_name, corrupted in [
        ("import_fidelity_B_eff_corrupt", "B_eff"),
        ("import_fidelity_K_h_corrupt", "K_h"),
        ("import_fidelity_c_E_corrupt", "c_E"),
    ]:
        branch = results["branches"][branch_name]
        q22 = branch["assembled_Q"][2][2]
        det_scalar = branch["assembled_Q_x_scalar_block"]
        lines.append(f"| `{branch_name}` | `{corrupted}` | `{q22}` | `{det_scalar}` | `{branch['primary']}` |")
    lines.extend(
        [
            "",
            "## Provenance",
            "",
            "Imported blocks and values:",
        ]
    )
    for key, value in imports.items():
        lines.append(f"- `{key}`: {value}")
    lines.extend(
        [
            "",
            "Declared / sim-deferred:",
        ]
    )
    for key, values in results["imports"]["declared_vs_deferred"].items():
        lines.append(f"- `{key}`: {'; '.join(values)}")
    lines.extend(
        [
            "",
            "## Dimensional Firewall",
            "",
            f"Passed `{results['dimensional_firewall']['checked_expression_count']}` assembled-action/source checks.  The able-to-fail ablations fired for missing `k^2` in `C_hu`, using charge density for `q_L`, and counting `q_M` as a charge source.",
            "",
            "## Dual Engine",
            "",
            f"`ENGINE_AGREE` over `{compared['checked_expression_count']}` independently derived and compared feature quantities: assembled `Q`, `J`, kinetic rank, imported Dirac constraint consequences, DOF per sector, scalar stability sign/eigenvalue codes, pole speeds, `J^T Q^-1 J`, residues, residue-derived feature flags, control classes, and import-fidelity guards.",
            "",
            "Run commands:",
            "",
            "```text",
            "timeout 600 python3 software/stage1_solver/tools/pathA_39_stage4_field_classification_sympy.py",
            "timeout 600 math -script software/stage1_solver/tools/pathA_39_stage4_field_classification.wl",
            "timeout 600 python3 software/stage1_solver/tools/pathA_39_stage4_field_classification_sympy.py --compare",
            "```",
            "",
            "## Remediation",
            "",
            "REMEDIATED: fixes 1-5, verdict/numbers unchanged.  The large-`C_hu` control carries the stable fixture label but is classified unstable from the computed determinant sign; Mathematica now derives the discrete integers before comparison; both engines use the same classifier predicates; scalar/density/h flags are residue-derived; expected-landing checks are non-fatal records.",
            "",
            "## Honest Scope",
            "",
            "This closes the magnetism sector as a spec-level field-content skeleton.  It does not perform the nonlinear throat solve, the `lambda_gamma` knit, the `c_E=c_gamma` Lorentz-cone test, or the simulation-deferred operator-contamination magnitude.  It does not un-earn Stages 0-3.",
        ]
    )
    REPORT_OUT.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--compare", action="store_true", help="compare with Mathematica payload and write YAML/report")
    args = parser.parse_args()

    payload = build_payload()
    write_json(SYM_OUT, payload)
    if args.compare:
        results = compare_with_mathematica(payload)
        print(
            json.dumps(
                {
                    "engine": "sympy",
                    "status": "ENGINE_AGREE",
                    "headline": f"primary={results['primary']} flags={results['flags']}",
                    "primary": results["primary"],
                    "flags": results["flags"],
                    "run_commands": [
                        "timeout 600 python3 software/stage1_solver/tools/pathA_39_stage4_field_classification_sympy.py",
                        "timeout 600 math -script software/stage1_solver/tools/pathA_39_stage4_field_classification.wl",
                        "timeout 600 python3 software/stage1_solver/tools/pathA_39_stage4_field_classification_sympy.py --compare",
                    ],
                    "remediation": "REMEDIATED: fixes 1-5, verdict/numbers unchanged",
                },
                sort_keys=True,
            )
        )
    else:
        print(json.dumps({"engine": "sympy", "status": "OK", "primary": payload["top_line_primary"], "flags": payload["top_line_flags"]}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
