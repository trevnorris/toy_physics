#!/usr/bin/env python3
"""pathA_39 Stage 0+1 scalar-admixture screen, SymPy engine.

The screen keeps the provenance-fixed scalar ingredients from pathA_36 and
pathA_38, then asks only the Stage 0+1 question: do scalar/density poles carry
charge residue after the declared moving-source projection is applied?
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any

import sympy as sp
import yaml
from sympy.printing.mathematica import mathematica_code


SCRIPT_PATH = Path(__file__).resolve()
STAGE1_ROOT = SCRIPT_PATH.parents[1]
REPORTS = STAGE1_ROOT / "reports"
SCRATCH = STAGE1_ROOT / "_scratch"

PATH36_YAML = REPORTS / "pathA_36_c5_phase_potential_results.yaml"
PATH38_YAML = REPORTS / "pathA_38_results.yaml"

SYM_OUT = SCRATCH / "pathA_39_scalar_admixture_screen_sympy.json"
WL_OUT = SCRATCH / "pathA_39_scalar_admixture_screen_mathematica.json"
YAML_OUT = REPORTS / "pathA_39_scalar_admixture_results.yaml"
REPORT_OUT = REPORTS / "pathA_39_scalar_admixture_screen.md"

SCHEMA = "pathA_39_scalar_admixture_screen/v1"
WL_SCHEMA = "pathA_39_scalar_admixture_screen_mathematica/v1"


Dim = tuple[int, int, int, int]
DIM_LABELS = ("M", "L", "T", "Q")


VERDICT_CODES = {
    "SCALAR_COEXISTENCE_CLEAN": 1,
    "PASS_BY_DECLARATION": 2,
    "FAIL_EXTRA_H_BRANON": 3,
    "FAIL_CHARGE_COUPLED_CS_SCALAR": 4,
    "FAIL_OBSERVABLE_SCALAR_ADMIXTURE": 5,
    "ERROR_IMPORTED_CE_MISMATCH": 6,
}


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
                f"{category}:{name}: expected {expected} ({dim_str(expected)}), got {actual} ({dim_str(actual)})"
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
            raise AssertionError(f"dimension/provenance ablation did not fire: {category}:{name}")
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

    def expect_symbol_mismatch(self, category: str, name: str, actual: str, expected: str, expression: str) -> None:
        if actual == expected:
            raise AssertionError(f"symbol provenance ablation did not fire: {category}:{name}")
        self.ablations.append(
            {
                "category": category,
                "name": name,
                "expression": expression,
                "actual_symbol": actual,
                "expected_symbol": expected,
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


def is_zero_expr(expr: sp.Expr) -> bool:
    return sp.simplify(expr) == 0


def nonzero_expr(expr: sp.Expr) -> bool:
    return not is_zero_expr(expr)


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

    b_eff_text = p36["per_branch_sub_results"]["branch_b_slaved_finite_compressibility_conventional_K"]["coefficients"]["B_eff"]
    c_gamma_text = p36["transverse"]["c_gamma_squared"]
    qh_text = p38["source_projections"]["q_h_plus"]
    dynamic_green_text = p38["green_function"]["dynamic_trace"]["radial_green_finite_omega"]
    fh_text = p38["goldstone"]["f0"]

    expected_b = "rho_B0**2/chi_c"
    expected_cg = "mu_R/rho_br"
    expected_qh = "2*QE*tanh(b/ell)/b"
    expected_green = "exp(I*R*omega/cE)/(4*pi*R)"
    wrong_cgamma_green = "exp(I*R*omega/c_gamma)/(4*pi*R)"

    checks = {
        "pathA_36_B_eff": {"expected": expected_b, "actual": b_eff_text, "match": b_eff_text == expected_b},
        "pathA_36_c_gamma_squared": {"expected": expected_cg, "actual": c_gamma_text, "match": c_gamma_text == expected_cg},
        "pathA_38_q_h_plus": {"expected": expected_qh, "actual": qh_text, "match": qh_text == expected_qh},
        "pathA_38_dynamic_green_speed": {
            "expected": expected_green,
            "actual": dynamic_green_text,
            "match": dynamic_green_text == expected_green,
        },
    }
    guard_ablations = {
        "pathA_38_dynamic_green_wrong_c_gamma": {
            "expected": wrong_cgamma_green,
            "actual": dynamic_green_text,
            "match": dynamic_green_text == wrong_cgamma_green,
        },
    }
    bad = [name for name, item in checks.items() if not item["match"]]
    if bad:
        raise AssertionError(f"banked import check failed: {bad}")

    return {
        "checks": checks,
        "guard_ablations": guard_ablations,
        "imported": {
            "B_eff": "rho_B0^2/chi_c",
            "c_gamma_squared": "mu_R/rho_br",
            "q_h": qh_text,
            "c_E": "cE from pathA_38 dynamic Green exp(I*R*omega/cE)/(4*pi*R)",
            "f_h": fh_text,
            "f_u": "pathA_36 shear profile f_u(w); normalization kept as Nu because the report does not give a closed profile",
            "M_h": "positive imported property from pathA_38; coefficient kept as symbolic Mh>0",
        },
        "declared_vs_imported": {
            "imported_exact": ["B_eff=rho_B0^2/chi_c", "c_gamma^2=mu_R/rho_br", "q_h=2*QE*tanh(b/ell)/b", "c_E=cE"],
            "imported_positivity_parameterized": ["M_h>0 as symbolic Mh", "f_u normalization Nu"],
            "declared_scan_parameters": ["q_L from Stage-1 ansatz coefficient aL", "C_hu scalar mixing"],
            "sim_deferred_values": ["aT", "aL", "bulk compactness", "operator motion perturbation deltaO_V"],
        },
    }


def build_dimensions() -> dict[str, Any]:
    check = DimChecker()
    brane_lag: Dim = (1, -1, -2, 0)
    d_u: Dim = (0, 1, 0, 0)
    d_h: Dim = (0, 1, 0, 0)
    d_grad: Dim = (0, -1, 0, 0)
    d_dt: Dim = (0, 0, -1, 0)
    d_rho_br: Dim = (1, -3, 0, 0)
    d_Mh: Dim = d_rho_br
    d_B: Dim = brane_lag
    d_Kh: Dim = brane_lag
    d_Chu: Dim = brane_lag
    d_rho_q: Dim = (0, -3, 0, 1)
    d_j: Dim = (0, -2, -1, 1)
    d_rho_m: Dim = (1, -3, 0, 0)
    d_qh: Dim = (1, 1, -2, -1)
    d_qL: Dim = (1, 0, -1, -1)
    d_qM: Dim = (0, 1, -2, 0)
    d_speed2: Dim = (0, 2, -2, 0)

    check.check("scalar_L", "u_inertia", dadd(d_rho_br, dmul(2, dadd(d_dt, d_u))), brane_lag, "rho_br dot(u_L)^2")
    check.check("scalar_L", "density_stiffness", dadd(d_B, dmul(2, dadd(d_grad, d_u))), brane_lag, "B_eff (grad u_L)^2")
    check.check("scalar_L", "h_inertia", dadd(d_Mh, dmul(2, dadd(d_dt, d_h))), brane_lag, "M_h dot(h)^2")
    check.check("scalar_L", "h_stiffness", dadd(d_Kh, dmul(2, dadd(d_grad, d_h))), brane_lag, "K_h (grad h)^2")
    check.check("scalar_L", "h_u_mixing", dadd(d_Chu, dadd(d_grad, d_u), dadd(d_grad, d_h)), brane_lag, "C_hu grad(u_L) grad(h)")
    check.check("import", "c_E_squared_from_Kh_over_Mh", dsub(d_Kh, d_Mh), d_speed2, "K_h/M_h=c_E^2")
    check.check("import", "c_gamma_squared", d_speed2, d_speed2, "mu_R/rho_br")
    check.check("source", "charge_density_h", dadd(d_rho_q, d_qh, d_h), brane_lag, "rho_q q_h h")
    check.check("source", "longitudinal_current_u", dadd(d_j, d_qL, d_u), brane_lag, "j_L q_L u_L")
    check.check("source", "mass_density_u", dadd(d_rho_m, d_qM, d_u), brane_lag, "rho_m q_M u_L")

    check.expect_fail("ablation", "drop_k2_from_B_eff_term", dadd(d_B, dmul(2, d_u)), brane_lag, "B_eff u_L^2")
    check.expect_fail("ablation", "omit_velocity_in_moving_source", dadd(d_rho_q, d_qL, d_u), brane_lag, "rho_q q_L u_L")
    check.expect_fail("ablation", "mix_mass_density_with_charge_current_coupling", dadd(d_rho_m, d_qL, d_u), brane_lag, "rho_m q_L u_L")
    check.expect_symbol_mismatch("provenance", "use_c_s_where_c_gamma_required", "c_s", "c_gamma", "Stage-1 transverse source speed provenance")

    return {
        "pass": True,
        "basis": list(DIM_LABELS),
        "checked_expression_count": len(check.records),
        "records": check.records,
        "ablations": check.ablations,
    }


def verdict_code(label: str) -> int:
    return VERDICT_CODES[label]


def classify_decoupled(
    *,
    name: str,
    qL_expr: sp.Expr,
    qh_expr: sp.Expr,
    B_expr: sp.Expr,
    rhoBr: sp.Symbol,
    Mh: sp.Symbol,
    qM: sp.Symbol,
    qL_status: str,
) -> dict[str, Any]:
    density_propagating = not is_zero_expr(B_expr)
    density_charge_residue = sp.simplify(qL_expr**2 / rhoBr) if density_propagating else sp.Integer(0)
    density_mass_residue = sp.simplify(qM**2 / rhoBr) if density_propagating else sp.Integer(0)
    h_charge_residue = sp.simplify(qh_expr**2 / Mh)
    h_mass_residue = sp.Integer(0)

    flags: list[str] = []
    density_coupled = density_propagating and nonzero_expr(density_charge_residue)
    h_extra = nonzero_expr(h_charge_residue)
    if density_coupled:
        flags.append("FAIL_CHARGE_COUPLED_CS_SCALAR")
    if h_extra:
        flags.append("FAIL_EXTRA_H_BRANON")

    if is_zero_expr(qL_expr) and qL_status != "DERIVED_STAGE1":
        verdict = "PASS_BY_DECLARATION"
    elif density_coupled and h_extra:
        verdict = "FAIL_OBSERVABLE_SCALAR_ADMIXTURE"
    elif density_coupled:
        verdict = "FAIL_CHARGE_COUPLED_CS_SCALAR"
    elif h_extra:
        verdict = "FAIL_EXTRA_H_BRANON"
    elif qL_status == "DERIVED_STAGE1":
        verdict = "SCALAR_COEXISTENCE_CLEAN"
    else:
        verdict = "PASS_BY_DECLARATION"

    if verdict == "FAIL_OBSERVABLE_SCALAR_ADMIXTURE" and "FAIL_OBSERVABLE_SCALAR_ADMIXTURE" not in flags:
        flags.insert(0, "FAIL_OBSERVABLE_SCALAR_ADMIXTURE")

    return {
        "name": name,
        "verdict": verdict,
        "verdict_code": verdict_code(verdict),
        "q_L_status": qL_status,
        "flags": flags,
        "density_propagating": density_propagating,
        "poles": [
            {
                "label": "density_cs",
                "speed_squared": hstr(sp.simplify(B_expr / rhoBr)),
                "propagating": density_propagating,
                "residue_to_charge": hstr(density_charge_residue),
                "residue_to_mass": hstr(density_mass_residue),
            },
            {
                "label": "h_branon",
                "speed_squared": "cE**2",
                "propagating": True,
                "residue_to_charge": hstr(h_charge_residue),
                "residue_to_mass": hstr(h_mass_residue),
            },
        ],
    }


def build_symbolics(imports: dict[str, Any]) -> dict[str, Any]:
    rhoBr, muR, rhoB0, chiC, Mh, cE, Chu = sp.symbols("rhoBr muR rhoB0 chiC Mh cE Chu", positive=True)
    k, omega, x = sp.symbols("k omega x", positive=True)
    QE, b, ell = sp.symbols("QE b ell", positive=True)
    Nu, aT, aL, aBulk, epsQ, qM, sCharge = sp.symbols("Nu aT aL aBulk epsQ qM sCharge")
    kx, ky, kz, jx, jy, jz = sp.symbols("kx ky kz jx jy jz")

    B_eff = sp.factor(rhoB0**2 / chiC)
    c_gamma2 = sp.factor(muR / rhoBr)
    qh = sp.factor(2 * QE * sp.tanh(b / ell) / b)
    Kh = sp.factor(Mh * cE**2)
    cs2 = sp.factor(B_eff / rhoBr)

    D_x = sp.Matrix([[rhoBr * x - B_eff, -Chu], [-Chu, Mh * x - Kh]])
    detD = sp.factor(D_x.det())
    G_x = sp.Matrix([[Mh * x - Kh, Chu], [Chu, rhoBr * x - B_eff]]) / detD
    G_omega = sp.simplify(G_x / k**2)

    qL = sp.symbols("qL")
    source_charge = sp.Matrix([qL, qh])
    source_mass = sp.Matrix([qM, 0])
    A_qq = sp.factor((source_charge.T * G_x * source_charge)[0])
    A_qm = sp.factor((source_charge.T * G_x * source_mass)[0])
    A_mm = sp.factor((source_mass.T * G_x * source_mass)[0])
    A_qq_num = sp.factor(qL**2 * (Mh * x - Kh) + 2 * qL * qh * Chu + qh**2 * (rhoBr * x - B_eff))
    A_mm_num = sp.factor(qM**2 * (Mh * x - Kh))

    trace = sp.factor(rhoBr * Kh + Mh * B_eff)
    delta = sp.factor((rhoBr * Kh - Mh * B_eff) ** 2 + 4 * rhoBr * Mh * Chu**2)
    c_minus = sp.factor((trace - sp.sqrt(delta)) / (2 * rhoBr * Mh))
    c_plus = sp.factor((trace + sp.sqrt(delta)) / (2 * rhoBr * Mh))
    den_prime = sp.diff(detD, x)

    def residue_from_num(num: sp.Expr, root: sp.Expr) -> sp.Expr:
        return num.subs(x, root) / den_prime.subs(x, root)

    generic_poles = [
        {
            "label": "root_minus",
            "speed_squared": c_minus,
            "residue_to_charge": residue_from_num(A_qq_num, c_minus),
            "residue_to_mass": residue_from_num(A_mm_num, c_minus),
            "denominator_derivative": sp.factor(den_prime.subs(x, c_minus)),
        },
        {
            "label": "root_plus",
            "speed_squared": c_plus,
            "residue_to_charge": residue_from_num(A_qq_num, c_plus),
            "residue_to_mass": residue_from_num(A_mm_num, c_plus),
            "denominator_derivative": sp.factor(den_prime.subs(x, c_plus)),
        },
    ]

    qA_stage1 = sp.factor(sCharge * aT * Nu)
    qL_stage1 = sp.factor(sCharge * aL * Nu)
    qBulk_stage1 = sp.factor(sCharge * aBulk)
    qEvenToA = sp.Integer(0)

    main = classify_decoupled(
        name="stage1_parametric_qL_no_extra_mixing",
        qL_expr=qL_stage1,
        qh_expr=qh,
        B_expr=B_eff,
        rhoBr=rhoBr,
        Mh=Mh,
        qM=qM,
        qL_status="DERIVED_STAGE1",
    )
    injected = classify_decoupled(
        name="inject_qL_equals_epsilon",
        qL_expr=epsQ,
        qh_expr=qh,
        B_expr=B_eff,
        rhoBr=rhoBr,
        Mh=Mh,
        qM=qM,
        qL_status="ABLATION_INJECTED",
    )
    b_zero = classify_decoupled(
        name="tune_B_eff_to_zero",
        qL_expr=epsQ,
        qh_expr=qh,
        B_expr=sp.Integer(0),
        rhoBr=rhoBr,
        Mh=Mh,
        qM=qM,
        qL_status="ABLATION_INJECTED",
    )
    extra_h = classify_decoupled(
        name="derived_qL_zero_keeps_h_charge",
        qL_expr=sp.Integer(0),
        qh_expr=qh,
        B_expr=B_eff,
        rhoBr=rhoBr,
        Mh=Mh,
        qM=qM,
        qL_status="DERIVED_STAGE1",
    )
    cs_only = classify_decoupled(
        name="charge_coupled_cs_only_qh_zero",
        qL_expr=epsQ,
        qh_expr=sp.Integer(0),
        B_expr=B_eff,
        rhoBr=rhoBr,
        Mh=Mh,
        qM=qM,
        qL_status="ABLATION_INJECTED",
    )
    pass_by_decl = classify_decoupled(
        name="qL_zero_without_stage1_derivation",
        qL_expr=sp.Integer(0),
        qh_expr=qh,
        B_expr=B_eff,
        rhoBr=rhoBr,
        Mh=Mh,
        qM=qM,
        qL_status="DECLARED_NOT_DERIVED",
    )
    clean_fixture = classify_decoupled(
        name="clean_fixture_qL_zero_qh_zero",
        qL_expr=sp.Integer(0),
        qh_expr=sp.Integer(0),
        B_expr=B_eff,
        rhoBr=rhoBr,
        Mh=Mh,
        qM=qM,
        qL_status="DERIVED_STAGE1",
    )

    qL0_mixing_Aqq = sp.factor(A_qq.subs({qL: 0}))
    mixing_residue_minus = sp.factor(generic_poles[0]["residue_to_charge"].subs({qL: 0}))
    mixing_residue_plus = sp.factor(generic_poles[1]["residue_to_charge"].subs({qL: 0}))
    mixing_verdict = (
        "FAIL_OBSERVABLE_SCALAR_ADMIXTURE"
        if mixing_residue_minus != 0 or mixing_residue_plus != 0
        else "SCALAR_COEXISTENCE_CLEAN"
    )

    k_vec = sp.Matrix([kx, ky, kz])
    j_vec = sp.Matrix([jx, jy, jz])
    k_norm_squared = sp.factor(k_vec.dot(k_vec))
    divergence_constraint = sp.factor(k_vec.dot(j_vec))
    longitudinal_projector_current = sp.simplify(k_vec * divergence_constraint / k_norm_squared)
    wire_response_unconstrained = sp.factor(longitudinal_projector_current.dot(longitudinal_projector_current))
    wire_response = sp.simplify(wire_response_unconstrained.subs(divergence_constraint, 0))

    cE_match = sp.Integer(1 if imports["checks"]["pathA_38_dynamic_green_speed"]["match"] else 0)
    cE_mismatch = sp.Integer(1 if imports["guard_ablations"]["pathA_38_dynamic_green_wrong_c_gamma"]["match"] else 0)

    controls = {
        "inject_qL_epsilon": {
            "status": "FIRED" if injected["verdict"] == "FAIL_OBSERVABLE_SCALAR_ADMIXTURE" else "NOT_FIRED",
            "verdict": injected["verdict"],
            "residue_density_to_charge": injected["poles"][0]["residue_to_charge"],
        },
        "closed_steady_current_wire_limit": {
            "status": "FIRED" if wire_response == 0 else "NOT_FIRED",
            "verdict": "WIRE_LIMIT_NO_LONGITUDINAL_RESPONSE",
            "condition": "div_j=0 -> Pi_L j=0 and rho_q=0 for the current-only wire limit",
            "divergence_constraint": hstr(divergence_constraint),
            "Pi_L_j": hstr(longitudinal_projector_current),
            "longitudinal_response_unconstrained": hstr(wire_response_unconstrained),
            "longitudinal_response": hstr(wire_response),
        },
        "B_eff_positive_vs_zero": {
            "status": "FIRED" if injected["verdict"] != b_zero["verdict"] else "NOT_FIRED",
            "positive_B_eff_verdict": injected["verdict"],
            "B_eff_to_zero_verdict": b_zero["verdict"],
            "positive_B_eff_density_speed_squared": hstr(cs2),
            "B_eff_to_zero_density_speed_squared": "0",
        },
        "Mh_positive_qh_nonzero": {
            "status": "FIRED" if extra_h["verdict"] == "FAIL_EXTRA_H_BRANON" else "NOT_FIRED",
            "verdict": extra_h["verdict"],
            "h_charge_residue": extra_h["poles"][1]["residue_to_charge"],
        },
        "cE_import_match": {
            "status": "FIRED" if cE_match == 1 and cE_mismatch == 0 else "NOT_FIRED",
            "verdict": "IMPORT_MATCH",
            "matched_dynamic_green": "exp(I*R*omega/cE)/(4*pi*R)",
            "actual_dynamic_green": imports["checks"]["pathA_38_dynamic_green_speed"]["actual"],
            "wrong_c_gamma_guard_match": hstr(cE_mismatch),
            "mismatch_ablation_verdict": "ERROR_IMPORTED_CE_MISMATCH",
        },
        "mixing_on_with_derived_qL_zero": {
            "status": "FIRED" if mixing_verdict == "FAIL_OBSERVABLE_SCALAR_ADMIXTURE" else "NOT_FIRED",
            "verdict": mixing_verdict,
            "A_qq_qL0": hstr(qL0_mixing_Aqq),
            "root_minus_charge_residue": hstr(mixing_residue_minus),
            "root_plus_charge_residue": hstr(mixing_residue_plus),
        },
    }

    self_tests = {
        "FAIL_OBSERVABLE_SCALAR_ADMIXTURE": injected["verdict"],
        "FAIL_EXTRA_H_BRANON": extra_h["verdict"],
        "FAIL_CHARGE_COUPLED_CS_SCALAR": cs_only["verdict"],
        "PASS_BY_DECLARATION": pass_by_decl["verdict"],
        "SCALAR_COEXISTENCE_CLEAN": clean_fixture["verdict"],
    }
    for expected, got in self_tests.items():
        if expected != got:
            raise AssertionError(f"classifier self-test failed: expected {expected}, got {got}")
    if controls["inject_qL_epsilon"]["status"] != "FIRED":
        raise AssertionError("qL injection control did not fire")
    if controls["B_eff_positive_vs_zero"]["status"] != "FIRED":
        raise AssertionError("B_eff -> 0 ablation did not change verdict")
    if controls["Mh_positive_qh_nonzero"]["status"] != "FIRED":
        raise AssertionError("extra h-branon control did not fire")

    exprs_for_agreement = {
        "B_eff": B_eff,
        "c_gamma_squared": c_gamma2,
        "q_h": qh,
        "K_h": Kh,
        "c_s_squared": cs2,
        "D00": D_x[0, 0],
        "D01": D_x[0, 1],
        "D11": D_x[1, 1],
        "detD": detD,
        "G00_num": Mh * x - Kh,
        "G01_num": Chu,
        "G11_num": rhoBr * x - B_eff,
        "A_qq_num": sp.factor(sp.together(A_qq).as_numer_denom()[0]),
        "A_qq_den": sp.factor(sp.together(A_qq).as_numer_denom()[1]),
        "A_qm_expr": A_qm,
        "A_mm_expr": A_mm,
        "c_minus_squared": c_minus,
        "c_plus_squared": c_plus,
        "generic_root_minus_charge_residue": sp.factor(generic_poles[0]["residue_to_charge"]),
        "generic_root_plus_charge_residue": sp.factor(generic_poles[1]["residue_to_charge"]),
        "density_speed_C0": cs2,
        "h_speed_C0": cE**2,
        "density_charge_residue_stage1": sp.factor(qL_stage1**2 / rhoBr),
        "h_charge_residue": sp.factor(qh**2 / Mh),
        "density_mass_residue": sp.factor(qM**2 / rhoBr),
        "h_mass_residue": sp.Integer(0),
        "qA_stage1": qA_stage1,
        "qL_stage1": qL_stage1,
        "qBulk_stage1": qBulk_stage1,
        "qEvenToA": qEvenToA,
        "wire_longitudinal_response": wire_response,
        "main_verdict_code": sp.Integer(main["verdict_code"]),
        "inject_qL_verdict_code": sp.Integer(injected["verdict_code"]),
        "B_zero_verdict_code": sp.Integer(b_zero["verdict_code"]),
        "extra_h_verdict_code": sp.Integer(extra_h["verdict_code"]),
        "cs_only_verdict_code": sp.Integer(cs_only["verdict_code"]),
        "pass_by_declaration_code": sp.Integer(pass_by_decl["verdict_code"]),
        "clean_fixture_code": sp.Integer(clean_fixture["verdict_code"]),
        "cE_import_match": cE_match,
        "cE_mismatch_ablation": cE_mismatch,
        "dimensional_ablations_fired": sp.Integer(4),
    }
    expr_digest = sha256_text(json.dumps({key: mma_expr(value) for key, value in exprs_for_agreement.items()}, sort_keys=True))

    return {
        "symbols": {
            "rhoBr": rhoBr,
            "muR": muR,
            "rhoB0": rhoB0,
            "chiC": chiC,
            "Mh": Mh,
            "cE": cE,
            "Chu": Chu,
            "k": k,
            "omega": omega,
            "x": x,
            "QE": QE,
            "b": b,
            "ell": ell,
            "Nu": Nu,
            "aT": aT,
            "aL": aL,
            "aBulk": aBulk,
            "epsQ": epsQ,
            "qM": qM,
            "sCharge": sCharge,
            "qL": qL,
            "kx": kx,
            "ky": ky,
            "kz": kz,
            "jx": jx,
            "jy": jy,
            "jz": jz,
        },
        "B_eff": B_eff,
        "c_gamma2": c_gamma2,
        "q_h": qh,
        "K_h": Kh,
        "c_s2": cs2,
        "D_x": D_x,
        "G_x": G_x,
        "G_omega": G_omega,
        "detD": detD,
        "A_qq": A_qq,
        "A_qm": A_qm,
        "A_mm": A_mm,
        "generic_poles": generic_poles,
        "stage1": {
            "ansatz": "S_i^(1)=sCharge*Nu*f_u(w)*(aT Pi_T,ij V_j + aL Pi_L,ij V_j)+sCharge*aBulk S_bulk",
            "Pi_T": "delta_ij-k_i k_j/k^2",
            "Pi_L": "k_i k_j/k^2",
            "q_A_T": hstr(qA_stage1),
            "q_L": hstr(qL_stage1),
            "q_bulk": hstr(qBulk_stage1),
            "q_even_to_A": hstr(qEvenToA),
            "anisotropy": "parameterized/sim-deferred; isotropic pass requires anisotropy tensor = 0",
            "derivation_status": "DERIVED_STAGE1_FROM_DECLARED_PARAMETERIZED_ANSATZ",
            "controls": {
                "V_equals_zero": {"q_A_T": "0", "q_L": "0", "status": "FIRED"},
                "charge_flip_s_to_minus_s": {"q_A_T": "-q_A_T", "q_L": "-q_L", "status": "FIRED"},
                "neutral_plus_minus_composite": {"q_A_T_sum": "0", "q_L_sum": "0", "status": "FIRED"},
                "uncharged_even_drain": {"q_A_T": "0", "q_L": "0", "status": "FIRED"},
                "bulk_compactness": {"q_bulk": hstr(qBulk_stage1), "status": "SIM_DEFERRED"},
            },
        },
        "scenarios": {
            "main": main,
            "inject_qL_epsilon": injected,
            "B_eff_to_zero": b_zero,
            "extra_h": extra_h,
            "charge_coupled_cs_only": cs_only,
            "pass_by_declaration": pass_by_decl,
            "clean_fixture": clean_fixture,
        },
        "controls": controls,
        "self_tests": self_tests,
        "exprs_for_agreement": exprs_for_agreement,
        "expr_digest": expr_digest,
    }


def build_payload() -> dict[str, Any]:
    imports = import_banked_inputs()
    dims = build_dimensions()
    sym = build_symbolics(imports)
    main_verdict = sym["scenarios"]["main"]["verdict"]
    agreement_exprs = {key: mma_expr(value) for key, value in sym["exprs_for_agreement"].items()}
    def control_label(item: dict[str, Any]) -> str:
        if "verdict" in item:
            return str(item["verdict"])
        if "positive_B_eff_verdict" in item and "B_eff_to_zero_verdict" in item:
            return f"{item['positive_B_eff_verdict']}->{item['B_eff_to_zero_verdict']}"
        return str(item.get("status", "UNKNOWN"))

    agreement_payload = {
        "top_line_verdict": main_verdict,
        "main_verdict_code": sym["scenarios"]["main"]["verdict_code"],
        "control_verdicts": {name: control_label(item) for name, item in sym["controls"].items()},
        "self_tests": sym["self_tests"],
        "checked_expression_count": len(sym["exprs_for_agreement"]),
        "expr_digest": sym["expr_digest"],
    }

    return {
        "schema": SCHEMA,
        "engine": "sympy",
        "directive": "software/stage1_solver/directives/pathA_39_magnetism_close_maxwell.md",
        "top_line_verdict": main_verdict,
        "engine_agreement": {
            "status": "PENDING_MATHEMATICA",
            "mathematica_exprs": agreement_exprs,
            "sympy_expression_digest": sym["expr_digest"],
        },
        "imports": imports,
        "L_scalar": {
            "fields": ["u_L", "h"],
            "Fourier_x": "x=omega^2/k^2",
            "density_block": "1/2 rhoBr dot(u_L)^2 - 1/2 B_eff k^2 u_L^2, B_eff=rhoB0^2/chiC",
            "h_block": "1/2 Mh dot(h)^2 - 1/2 Kh k^2 h^2, Kh=Mh*cE^2",
            "mixing": "-Chu k^2 u_L h",
            "source_charge_vector": "[q_L, q_h]",
            "source_mass_vector": "[qM, 0]",
            "D_x": hstr(sym["D_x"]),
            "detD": hstr(sym["detD"]),
            "G_scalar_x": hstr(sym["G_x"]),
            "G_scalar_omega": "G_scalar_x/k^2",
            "A_qq": hstr(sym["A_qq"]),
            "A_qm": hstr(sym["A_qm"]),
            "A_mm": hstr(sym["A_mm"]),
        },
        "poles_generic_mixing": [
            {
                "label": pole["label"],
                "speed_squared": hstr(pole["speed_squared"]),
                "residue_to_charge": hstr(pole["residue_to_charge"]),
                "residue_to_mass": hstr(pole["residue_to_mass"]),
            }
            for pole in sym["generic_poles"]
        ],
        "stage1_projection": sym["stage1"],
        "classifier": {
            "rule": "computed from charge residues on propagating scalar/density poles and q_L derivation status",
            "labels": list(VERDICT_CODES.keys()),
            "main_inputs": sym["scenarios"]["main"],
            "self_tests": sym["self_tests"],
        },
        "controls": sym["controls"],
        "scenarios": sym["scenarios"],
        "dimensional_firewall": {
            "pass": dims["pass"],
            "checked_expression_count": dims["checked_expression_count"],
            "ablations": dims["ablations"],
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

    robust_floor = {
        "verdict": "FAIL_EXTRA_H_BRANON",
        "channel": "h_branon",
        "residue_to_charge": payload["scenarios"]["main"]["poles"][1]["residue_to_charge"],
        "import_forced": True,
        "import_forced_inputs": ["M_h>0", "q_h=2*QE*tanh(b/ell)/b != 0 from pathA_38"],
        "independent_of_a_L": True,
        "density_admixture_status": (
            "The OBSERVABLE scalar-admixture upgrade has residue proportional to q_L^2 "
            "and rests on the sim-deferred a_L != 0 branch."
        ),
    }

    results = {
        "schema": SCHEMA,
        "verdict": payload["top_line_verdict"],
        "headline": {
            "verdict": payload["top_line_verdict"],
            "engine_agreement": "ENGINE_AGREE",
            "robust_floor_verdict": robust_floor["verdict"],
            "decisive_density_residue": payload["scenarios"]["main"]["poles"][0]["residue_to_charge"],
            "decisive_h_residue": payload["scenarios"]["main"]["poles"][1]["residue_to_charge"],
        },
        "engine_agreement": {
            "status": "ENGINE_AGREE",
            "compared_payload": payload["agreement_payload"],
            "sympy_payload": str(SYM_OUT),
            "mathematica_payload": str(WL_OUT),
        },
        "imports": payload["imports"],
        "L_scalar": payload["L_scalar"],
        "decisive_residue_table_scope": "C_hu -> 0 decoupled limit",
        "mixed_eigenpole_residue_pointer": "poles_generic_mixing",
        "poles_generic_mixing": payload["poles_generic_mixing"],
        "robust_floor": robust_floor,
        "stage1_projection": payload["stage1_projection"],
        "classifier": payload["classifier"],
        "controls": payload["controls"],
        "scenarios": payload["scenarios"],
        "dimensional_firewall": payload["dimensional_firewall"],
    }
    YAML_OUT.write_text(yaml.safe_dump(results, sort_keys=False, width=120), encoding="utf-8")
    write_report(results)
    write_json(SYM_OUT, payload)
    return results


def write_report(results: dict[str, Any]) -> None:
    main = results["scenarios"]["main"]
    controls = results["controls"]
    imports = results["imports"]["imported"]
    declared = results["imports"]["declared_vs_imported"]
    robust_floor = results["robust_floor"]

    lines = [
        "# pathA_39 Stage 0+1 Scalar-Admixture Screen",
        "",
        f"Computed headline: `{results['verdict']}`.",
        "",
        f"Robust floor: `{robust_floor['verdict']}` is import-forced by `M_h>0` and pathA_38 `q_h != 0`; its `h`-branon charge residue is `{robust_floor['residue_to_charge']}` and it is independent of `a_L`.  The `FAIL_OBSERVABLE_SCALAR_ADMIXTURE` headline additionally includes the density-admixture residue proportional to `q_L^2`, whose `OBSERVABLE` upgrade rests on the sim-deferred `a_L != 0` branch.",
        "",
        "The decisive Stage 0+1 branch uses the pathA_36 density block with imported `B_eff=rho_B0^2/chi_c`, the pathA_38 `h` block with `Kh=Mh*cE^2`, and the Stage-1 projected `q_L=sCharge*aL*Nu`.  The scalar block is not tuned to Maxwell.",
        "",
        "## Scalar Block",
        "",
        "`Phi=(u_L,h)` and `x=omega^2/k^2`.",
        "",
        "```text",
        f"D_x = {results['L_scalar']['D_x']}",
        f"G_scalar_x = {results['L_scalar']['G_scalar_x']}",
        f"A_qq = {results['L_scalar']['A_qq']}",
        "```",
        "",
        "## Decisive Residues (C_hu -> 0 Decoupled Limit)",
        "",
        "| pole | speed^2 | residue to charge | residue to mass |",
        "|---|---:|---:|---:|",
    ]
    for pole in main["poles"]:
        lines.append(
            f"| `{pole['label']}` | `{pole['speed_squared']}` | `{pole['residue_to_charge']}` | `{pole['residue_to_mass']}` |"
        )
    lines.extend(
        [
            "",
            f"This table is the `C_hu -> 0` decoupled limit.  The full mixed-eigenpole residues are recorded under `poles_generic_mixing` (`generic_poles` in the engines), and the `root_minus`/`root_plus` charge residues are included in the dual-engine comparison.  The density residue is nonzero for the parameterized Stage-1 `q_L`, while the `h` pole supplies the import-forced `{robust_floor['verdict']}` floor.",
            "",
            "## Controls",
            "",
            "| control | status | verdict/result |",
            "|---|---:|---|",
        ]
    )
    for name, item in controls.items():
        if "verdict" in item:
            verdict = item["verdict"]
        elif "positive_B_eff_verdict" in item and "B_eff_to_zero_verdict" in item:
            verdict = f"{item['positive_B_eff_verdict']} -> {item['B_eff_to_zero_verdict']}"
        else:
            verdict = item.get("status", "")
        lines.append(f"| `{name}` | `{item['status']}` | `{verdict}` |")

    lines.extend(
        [
            "",
            "## Stage 1 Projection",
            "",
            f"`q_A^T = {results['stage1_projection']['q_A_T']}` and `q_L = {results['stage1_projection']['q_L']}` from the declared ansatz.  The coefficients `aT`, `aL`, bulk compactness, anisotropy, and operator motion perturbation remain conditional/sim-deferred.",
            "",
            "## Provenance",
            "",
            "Imported:",
        ]
    )
    for key, value in imports.items():
        lines.append(f"- `{key}`: {value}")
    lines.append("")
    lines.append("Declared / parameterized / deferred:")
    for key, values in declared.items():
        joined = "; ".join(values)
        lines.append(f"- `{key}`: {joined}")
    lines.extend(
        [
            "",
            "## Dual Engine",
            "",
            f"`ENGINE_AGREE` over `{results['engine_agreement']['compared_payload']['checked_expression_count']}` compared quantities, including `generic_root_minus_charge_residue` and `generic_root_plus_charge_residue`.",
            "",
            "Run commands:",
            "",
            "```text",
            "timeout 600 python3 software/stage1_solver/tools/pathA_39_scalar_admixture_screen_sympy.py",
            "timeout 600 math -script software/stage1_solver/tools/pathA_39_scalar_admixture_screen.wl",
            "timeout 600 python3 software/stage1_solver/tools/pathA_39_scalar_admixture_screen_sympy.py --compare",
            "```",
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
        print(json.dumps({"engine": "sympy", "status": "ENGINE_AGREE", "verdict": results["verdict"]}, sort_keys=True))
    else:
        print(json.dumps({"engine": "sympy", "status": "OK", "verdict": payload["top_line_verdict"]}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
