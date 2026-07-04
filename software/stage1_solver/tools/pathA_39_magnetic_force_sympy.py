#!/usr/bin/env python3
"""pathA_39 Stage 2 magnetic-force derivation, SymPy engine.

This script keeps the Stage 2 calculation conditional on the declared
Stage-1 moving-source projection.  It assembles the MacCullagh transverse
sector and the imported longitudinal block, integrates out the fields, and
checks the same interaction by the source-1 field acting on source 2.
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
PATH39_S01_YAML = REPORTS / "pathA_39_scalar_admixture_results.yaml"

SYM_OUT = SCRATCH / "pathA_39_magnetic_force_sympy.json"
WL_OUT = SCRATCH / "pathA_39_magnetic_force_mathematica.json"
YAML_OUT = REPORTS / "pathA_39_magnetic_force_results.yaml"
REPORT_OUT = REPORTS / "pathA_39_magnetic_force.md"

SCHEMA = "pathA_39_magnetic_force/v1"
WL_SCHEMA = "pathA_39_magnetic_force_mathematica/v1"


Dim = tuple[int, int, int, int]
DIM_LABELS = ("M", "L", "T", "Q")

VERDICT_CODES = {
    "MAGNETIC_FORCE_DERIVED": 1,
    "FAIL_WRONG_FALLOFF": 2,
    "FAIL_TARGET_READBACK": 3,
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

    def expect_symbol_mismatch(self, category: str, name: str, actual: str, expected: str, expression: str) -> None:
        if actual == expected:
            raise AssertionError(f"provenance ablation did not fire: {category}:{name}")
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


def load_yaml(path: Path) -> dict[str, Any]:
    if not path.exists():
        raise FileNotFoundError(path)
    with path.open("r", encoding="utf-8") as fh:
        data = yaml.safe_load(fh)
    if not isinstance(data, dict):
        raise AssertionError(f"YAML did not load to mapping: {path}")
    return data


def verdict_code(label: str) -> int:
    return VERDICT_CODES[label]


def import_banked_inputs() -> dict[str, Any]:
    p36 = load_yaml(PATH36_YAML)
    p38 = load_yaml(PATH38_YAML)
    p39 = load_yaml(PATH39_S01_YAML)

    b_eff_text = p36["per_branch_sub_results"]["branch_b_slaved_finite_compressibility_conventional_K"]["coefficients"]["B_eff"]
    c_gamma_text = p36["transverse"]["c_gamma_squared"]
    transverse_omega2_text = p36["engine_agreement"]["compared_payload"]["transverse_omega2"]
    qh_text = p38["source_projections"]["q_h_plus"]
    dynamic_green_text = p38["green_function"]["dynamic_trace"]["radial_green_finite_omega"]
    fh_text = p38["goldstone"]["f0"]
    q_a_stage1 = p39["stage1_projection"]["q_A_T"]
    q_l_stage1 = p39["stage1_projection"]["q_L"]
    stage01_engine = p39["engine_agreement"]["status"]

    expected_b = "rho_B0**2/chi_c"
    expected_cg = "mu_R/rho_br"
    expected_omega2 = "mu_R*k^2/rho_br"
    expected_qh = "2*QE*tanh(b/ell)/b"
    expected_green = "exp(I*R*omega/cE)/(4*pi*R)"
    expected_q_a = "Nu*aT*sCharge"
    expected_q_l = "Nu*aL*sCharge"

    checks = {
        "pathA_36_B_eff": {"expected": expected_b, "actual": b_eff_text, "match": b_eff_text == expected_b},
        "pathA_36_c_gamma_squared": {"expected": expected_cg, "actual": c_gamma_text, "match": c_gamma_text == expected_cg},
        "pathA_36_transverse_omega2": {
            "expected": expected_omega2,
            "actual": transverse_omega2_text,
            "match": transverse_omega2_text == expected_omega2,
        },
        "pathA_38_q_h_plus": {"expected": expected_qh, "actual": qh_text, "match": qh_text == expected_qh},
        "pathA_38_dynamic_green_speed": {
            "expected": expected_green,
            "actual": dynamic_green_text,
            "match": dynamic_green_text == expected_green,
        },
        "pathA_39_stage01_engine": {"expected": "ENGINE_AGREE", "actual": stage01_engine, "match": stage01_engine == "ENGINE_AGREE"},
        "pathA_39_stage01_q_A_T": {"expected": expected_q_a, "actual": q_a_stage1, "match": q_a_stage1 == expected_q_a},
        "pathA_39_stage01_q_L": {"expected": expected_q_l, "actual": q_l_stage1, "match": q_l_stage1 == expected_q_l},
    }
    bad = [name for name, item in checks.items() if not item["match"]]
    if bad:
        raise AssertionError(f"banked import check failed: {bad}")

    guard_ablations = {
        "wrong_transverse_speed_cs_guard": {
            "expected_symbol": "c_gamma",
            "actual_symbol": "c_s",
            "match": False,
        },
        "dynamic_green_wrong_c_gamma_guard": {
            "expected": "exp(I*R*omega/c_gamma)/(4*pi*R)",
            "actual": dynamic_green_text,
            "match": dynamic_green_text == "exp(I*R*omega/c_gamma)/(4*pi*R)",
        },
    }

    return {
        "checks": checks,
        "guard_ablations": guard_ablations,
        "imported": {
            "B_eff": "rho_B0^2/chi_c",
            "c_gamma_squared": "mu_R/rho_br",
            "pathA_36_transverse_sector": "1/2 rho_br dot(u_T)^2 - 1/2 mu_R (curl u_T)^2",
            "q_h": qh_text,
            "c_E": "cE from pathA_38 dynamic Green exp(I*R*omega/cE)/(4*pi*R)",
            "f_h": fh_text,
            "f_u": "pathA_36 shear profile f_u(w); normalization kept as Nu",
            "M_h": "positive imported property from pathA_38; not used to tune the Stage-2 vector kernel",
            "stage0_plus_1": "pathA_39 scalar-admixture screen with q_A^T=Nu*aT*sCharge and q_L=Nu*aL*sCharge",
        },
        "declared_vs_imported": {
            "imported_exact": [
                "B_eff=rho_B0^2/chi_c",
                "c_gamma^2=mu_R/rho_br",
                "u_T MacCullagh shear sector from pathA_36",
                "q_h=2*QE*tanh(b/ell)/b",
                "c_E=cE",
                "Stage 0+1 q_A^T and q_L projection form",
            ],
            "declared_stage2": [
                "two compact brane-localized moving throats",
                "j_T=q_A^T V with q_A^T=Nu*aT*s",
                "j_L=q_L V with q_L=Nu*aL*s",
                "static O(V1 V2) exchange limit in the medium rest frame",
            ],
            "sim_deferred_values": [
                "aT",
                "aL",
                "bulk compactness and finite-mouth form factors",
                "anisotropy tensor",
                "operator motion perturbation deltaO_V",
                "absolute hierarchy magnitude",
            ],
        },
    }


def build_dimensions() -> dict[str, Any]:
    check = DimChecker()
    brane_lag: Dim = (1, -1, -2, 0)
    energy: Dim = (1, 2, -2, 0)
    force: Dim = (1, 1, -2, 0)
    d_u: Dim = (0, 1, 0, 0)
    d_grad: Dim = (0, -1, 0, 0)
    d_dt: Dim = (0, 0, -1, 0)
    d_rho_br: Dim = (1, -3, 0, 0)
    d_mu: Dim = brane_lag
    d_B: Dim = brane_lag
    d_speed2: Dim = (0, 2, -2, 0)
    d_charge: Dim = (0, 0, 0, 1)
    d_current_density: Dim = (0, -2, -1, 1)
    d_point_source_coeff: Dim = (1, 0, -1, 0)
    d_coupling_per_charge: Dim = dsub(d_point_source_coeff, d_charge)
    d_velocity: Dim = (0, 1, -1, 0)
    d_R: Dim = (0, 1, 0, 0)

    check.check("transverse_L", "u_T_inertia", dadd(d_rho_br, dmul(2, dadd(d_dt, d_u))), brane_lag, "rho_br dot(u_T)^2")
    check.check("transverse_L", "macCullagh_shear", dadd(d_mu, dmul(2, dadd(d_grad, d_u))), brane_lag, "mu_R (curl u_T)^2")
    check.check("longitudinal_L", "density_stiffness", dadd(d_B, dmul(2, dadd(d_grad, d_u))), brane_lag, "B_eff (div u_L)^2")
    check.check("import", "c_gamma_squared", dsub(d_mu, d_rho_br), d_speed2, "mu_R/rho_br")
    check.check("source_density", "current_couples_to_u", dadd(d_current_density, d_coupling_per_charge, d_u), brane_lag, "j q_A u")
    check.check(
        "point_kernel",
        "U_T_dimension",
        dsub(dadd(dmul(2, d_point_source_coeff), dmul(2, d_velocity)), dadd(d_mu, d_R)),
        energy,
        "(q_A^T)^2 V1 V2/(mu_R R)",
    )
    check.check(
        "point_kernel",
        "U_L_dimension",
        dsub(dadd(dmul(2, d_point_source_coeff), dmul(2, d_velocity)), dadd(d_B, d_R)),
        energy,
        "q_L^2 V1 V2/(B_eff R)",
    )
    check.check("point_force", "F_dimension", dsub(energy, d_R), force, "grad_R U")

    check.expect_fail(
        "ablation",
        "omit_velocity_in_current_kernel",
        dsub(dmul(2, d_point_source_coeff), dadd(d_mu, d_R)),
        energy,
        "(q_A^T)^2/(mu_R R)",
    )
    check.expect_fail(
        "ablation",
        "omit_compact_brane_delta_normalization",
        dsub(dadd(dmul(2, d_point_source_coeff), dmul(2, d_velocity), d_R), d_mu),
        energy,
        "(q_A^T)^2 V1 V2 R/mu_R",
    )
    check.expect_fail(
        "ablation",
        "mix_mass_current_as_charge_current",
        dadd((1, -2, -1, 0), d_coupling_per_charge, d_u),
        brane_lag,
        "j_mass q_A u",
    )
    check.expect_symbol_mismatch(
        "provenance",
        "use_c_s_where_c_gamma_required",
        "c_s",
        "c_gamma",
        "transverse propagator must use mu_R/rho_br from pathA_36",
    )

    return {
        "pass": True,
        "basis": list(DIM_LABELS),
        "checked_expression_count": len(check.records),
        "records": check.records,
        "ablations": check.ablations,
    }


def classify_force(*, potential_power: int, force_power: int, brane_localized: bool, target_response_changed: bool) -> dict[str, Any]:
    if not target_response_changed:
        verdict = "FAIL_TARGET_READBACK"
    elif potential_power != -1 or force_power != -2 or not brane_localized:
        verdict = "FAIL_WRONG_FALLOFF"
    else:
        verdict = "MAGNETIC_FORCE_DERIVED"
    return {
        "verdict": verdict,
        "verdict_code": verdict_code(verdict),
        "potential_power_R": potential_power,
        "force_power_R": force_power,
        "brane_localized": brane_localized,
        "target_response_changed": target_response_changed,
    }


def label_from_sign(sign: int) -> str:
    if sign < 0:
        return "ATTRACT"
    if sign > 0:
        return "REPEL"
    return "ZERO_FORCE"


def additive_sign_matrix() -> dict[str, dict[str, str]]:
    rows: dict[str, dict[str, str]] = {}
    for charge_label, charge_sign in (("like_charge_s1s2_plus", 1), ("opposite_charge_s1s2_minus", -1)):
        rows[charge_label] = {}
        for current_label, current_sign in (("like_parallel_currents", 1), ("opposite_parallel_currents", -1)):
            rows[charge_label][current_label] = label_from_sign(-charge_sign * current_sign)
    return rows


def build_symbolics(imports: dict[str, Any]) -> dict[str, Any]:
    rhoBr, muR, rhoB0, chiC, cE = sp.symbols("rhoBr muR rhoB0 chiC cE", positive=True)
    Nu, aT, aL = sp.symbols("Nu aT aL", positive=True)
    s1, s2 = sp.symbols("s1 s2")
    R = sp.symbols("R", positive=True)
    zeta, sigma = sp.symbols("zeta sigma", positive=True)
    k2 = sp.symbols("k2", positive=True)
    v1x, v1y, v1z, v2x, v2y, v2z, nx, ny, nz = sp.symbols("v1x v1y v1z v2x v2y v2z nx ny nz")
    vdot = sp.symbols("vDot")
    Ddot, arad, brad, Aprod = sp.symbols("Ddot arad brad Aprod")

    B_eff = sp.factor(rhoB0**2 / chiC)
    c_gamma2 = sp.factor(muR / rhoBr)
    qA0 = sp.factor(Nu * aT)
    qL0 = sp.factor(Nu * aL)
    qA1 = sp.factor(s1 * qA0)
    qA2 = sp.factor(s2 * qA0)
    qL1 = sp.factor(s1 * qL0)
    qL2 = sp.factor(s2 * qL0)
    s12 = sp.factor(s1 * s2)

    V1 = sp.Matrix([v1x, v1y, v1z])
    V2 = sp.Matrix([v2x, v2y, v2z])
    n = sp.Matrix([nx, ny, nz])
    substitutions = {
        Ddot: sp.factor(V1.dot(V2)),
        arad: sp.factor(V1.dot(n)),
        brad: sp.factor(V2.dot(n)),
        Aprod: sp.factor(V1.dot(n) * V2.dot(n)),
    }

    def radial_hessian_contract(seed: sp.Expr) -> sp.Expr:
        return sp.factor(sp.diff(seed, R, 2) * Aprod + (sp.diff(seed, R) / R) * (Ddot - Aprod))

    def inverse_fourier_projector_kernels(denominator_power: int) -> dict[str, sp.Expr]:
        if denominator_power == 2:
            scalar_green = sp.factor(1 / (4 * sp.pi * R))
            biharmonic_seed = sp.factor(-R / (8 * sp.pi))
        elif denominator_power == 4:
            scalar_green = sp.factor(-R / (8 * sp.pi))
            biharmonic_seed = sp.factor(R**3 / (96 * sp.pi))
        else:
            raise ValueError(f"unsupported denominator power: {denominator_power}")
        kk_kernel = sp.factor(-radial_hessian_contract(biharmonic_seed))
        delta_kernel = sp.factor(Ddot * scalar_green)
        return {
            "delta": delta_kernel,
            "kk": kk_kernel,
            "T": sp.factor(delta_kernel - kk_kernel),
            "L": sp.factor(kk_kernel),
            "scalar_green": scalar_green,
            "biharmonic_seed": biharmonic_seed,
        }

    def force_from_potential(U: sp.Expr) -> sp.Matrix:
        grad_A = (brad * V1 + arad * V2 - 2 * Aprod * n) / R
        grad_U = sp.diff(U, R) * n + sp.diff(U, Aprod) * grad_A
        return sp.Matrix([sp.factor(x) for x in (-grad_U)])

    def expand_compact(expr: sp.Expr | sp.MatrixBase) -> sp.Expr | sp.Matrix:
        if isinstance(expr, sp.MatrixBase):
            return sp.Matrix([sp.factor(x.subs(substitutions)) for x in expr])
        return sp.factor(expr.subs(substitutions))

    def expression_nonzero(expr: sp.Expr | sp.MatrixBase) -> bool:
        if isinstance(expr, sp.MatrixBase):
            return any(sp.simplify(entry) != 0 for entry in expr)
        return sp.simplify(expr) != 0

    def r_power_scalar(expr: sp.Expr) -> int | None:
        expr = sp.factor(expr)
        if sp.simplify(expr) == 0:
            return None
        term_powers: set[int] = set()
        for term in sp.Add.make_args(sp.expand(expr)):
            power = sp.factor(term).as_powers_dict().get(R, sp.Integer(0))
            if not bool(power.is_integer):
                term_powers = set()
                break
            term_powers.add(int(power))
        if len(term_powers) == 1:
            return next(iter(term_powers))
        tau = sp.symbols("tau", positive=True)
        for power in range(-8, 9):
            if sp.simplify(expr.subs(R, tau * R) - tau**power * expr) == 0:
                return power
        raise AssertionError(f"could not measure R power from {sp.sstr(expr)}")

    def r_power(expr: sp.Expr | sp.MatrixBase) -> int:
        if isinstance(expr, sp.MatrixBase):
            powers = {r_power_scalar(entry) for entry in expr if sp.simplify(entry) != 0}
            powers.discard(None)
            if len(powers) != 1:
                raise AssertionError(f"matrix has nonuniform R powers: {powers}")
            return next(iter(powers))
        power = r_power_scalar(expr)
        if power is None:
            raise AssertionError("zero expression has no R power")
        return power

    # Static operator inversion in the projector eigenbasis.
    O_T_static = sp.Matrix([[muR * k2]])
    O_L_static = sp.Matrix([[B_eff * k2]])
    G_T_static_coeff = sp.factor(O_T_static.inv()[0, 0])
    G_L_static_coeff = sp.factor(O_L_static.inv()[0, 0])

    standard_kernels = inverse_fourier_projector_kernels(2)
    perturbed_kernels = inverse_fourier_projector_kernels(4)
    K_T = standard_kernels["T"]
    K_L = standard_kernels["L"]
    K_T_perturbed = perturbed_kernels["T"]
    K_L_perturbed = perturbed_kernels["L"]

    U_T_integrate = sp.factor(-(qA1 * qA2 / muR) * K_T)
    U_L_integrate = sp.factor(-(qL1 * qL2 / B_eff) * K_L)
    U_total_integrate = sp.factor(U_T_integrate + U_L_integrate)

    U_T_eom = sp.factor(-(qA1 * qA2 / muR) * (standard_kernels["delta"] - standard_kernels["kk"]))
    U_L_eom = sp.factor(-(qL1 * qL2 / B_eff) * standard_kernels["kk"])
    U_total_eom = sp.factor(U_T_eom + U_L_eom)

    F_T = force_from_potential(U_T_integrate)
    F_L = force_from_potential(U_L_integrate)
    F_total = sp.simplify(F_T + F_L)

    U_T_wrong_falloff = sp.factor(-(qA1 * qA2 / muR) * K_T_perturbed)
    U_L_wrong_falloff = sp.factor(-(qL1 * qL2 / B_eff) * K_L_perturbed)
    U_total_wrong_falloff = sp.factor(U_T_wrong_falloff + U_L_wrong_falloff)
    F_total_wrong_falloff = force_from_potential(U_total_wrong_falloff)

    potential_power = r_power(U_total_integrate)
    force_power = r_power(F_total)
    wrong_potential_power = r_power(U_total_wrong_falloff)
    wrong_force_power = r_power(F_total_wrong_falloff)

    derived_tracks_functional_perturbation = (
        potential_power != wrong_potential_power
        and force_power != wrong_force_power
        and expression_nonzero(F_total_wrong_falloff - F_total)
    )

    U_T_mu_ablate = sp.factor((-(qA1 * qA2 / (zeta * muR)) * K_T))
    mu_ablation_ratio = sp.factor(U_T_mu_ablate / U_T_integrate)
    U_T_source_ablate = sp.factor((-(s1 * Nu * sigma * aT) * (s2 * Nu * sigma * aT) / muR) * K_T)
    source_ablation_ratio = sp.factor(U_T_source_ablate / U_T_integrate)
    U_T_projection_ablate = sp.factor(-(qA1 * qA2 / muR) * K_L)
    projection_ablation_delta = sp.factor(U_T_projection_ablate - U_T_integrate)
    mu_ablation_changed = sp.simplify(mu_ablation_ratio - 1) != 0
    source_ablation_changed = sp.simplify(source_ablation_ratio - 1) != 0
    projection_ablation_changed = sp.simplify(projection_ablation_delta) != 0

    target_readback_fixture = sp.factor((s12 * Nu**2 * aT**2 / muR) * (Aprod - Ddot) / (8 * sp.pi * R))
    target_readback_fixture_perturbed = target_readback_fixture
    target_readback_delta = sp.factor(target_readback_fixture_perturbed - target_readback_fixture)
    target_readback_tracks = expression_nonzero(target_readback_delta)

    main = classify_force(
        potential_power=potential_power,
        force_power=force_power,
        brane_localized=True,
        target_response_changed=derived_tracks_functional_perturbation,
    )
    wrong_falloff = classify_force(
        potential_power=wrong_potential_power,
        force_power=wrong_force_power,
        brane_localized=True,
        target_response_changed=True,
    )
    noncompact = classify_force(
        potential_power=potential_power,
        force_power=force_power,
        brane_localized=False,
        target_response_changed=True,
    )
    target_readback = classify_force(
        potential_power=r_power(target_readback_fixture),
        force_power=r_power(target_readback_fixture),
        brane_localized=True,
        target_response_changed=target_readback_tracks,
    )

    self_tests = {
        "MAGNETIC_FORCE_DERIVED": main["verdict"],
        "FAIL_WRONG_FALLOFF": wrong_falloff["verdict"],
        "FAIL_TARGET_READBACK": target_readback["verdict"],
    }
    for expected, got in self_tests.items():
        if expected != got:
            raise AssertionError(f"classifier self-test failed: expected {expected}, got {got}")

    zero_velocity_U = sp.factor(expand_compact(U_total_integrate).subs({v1x: 0, v1y: 0, v1z: 0}))
    neutral_composite_U = sp.factor(U_total_integrate.subs({s1: 1, s2: 1}) + U_total_integrate.subs({s1: -1, s2: 1}))
    charge_flip_delta = sp.factor(U_total_integrate.subs({s1: -s1}) + U_total_integrate)
    lorentz_residual = sp.factor(cE**2 - c_gamma2)
    lorentz_on_cone_lock = sp.factor(lorentz_residual.subs(cE**2, c_gamma2))
    U_side = sp.factor(U_total_integrate.subs({Aprod: 0, arad: 0, brad: 0, Ddot: vdot}))
    F_side_vec = force_from_potential(U_total_integrate).subs({Aprod: 0, arad: 0, brad: 0, Ddot: vdot})
    F_side_radial_coeff = sp.factor(F_side_vec[0] / nx)
    scalar_admixture_ratio = sp.factor((qL0**2 / B_eff) / (qA0**2 / muR))
    raw_ratio = sp.factor(qL0**2 / qA0**2)
    raw_ratio_reference = sp.factor(B_eff / muR)
    transverse_side_vec = F_T.subs({Aprod: 0, arad: 0, brad: 0, Ddot: 1, s1: 1, s2: 1})
    longitudinal_side_vec = F_L.subs({Aprod: 0, arad: 0, brad: 0, Ddot: 1, s1: 1, s2: 1})
    total_side_vec = F_total.subs({Aprod: 0, arad: 0, brad: 0, Ddot: 1, s1: 1, s2: 1})
    transverse_like_parallel_side_force = sp.factor(transverse_side_vec[0] / nx)
    longitudinal_like_parallel_side_force = sp.factor(longitudinal_side_vec[0] / nx)
    total_like_parallel_side_force = sp.factor(total_side_vec[0] / nx)
    U_T_ghost = sp.factor(-(qA1 * qA2 / (-muR)) * K_T)
    transverse_ghost_side_force = sp.factor(force_from_potential(U_T_ghost).subs({Aprod: 0, arad: 0, brad: 0, Ddot: 1, s1: 1, s2: 1})[0] / nx)

    controls = {
        "muR_propagator_perturbation": {
            "status": "FIRED" if mu_ablation_changed else "NOT_FIRED",
            "verdict": "RE_DERIVED_RESPONSE_CHANGED",
            "ratio": hstr(mu_ablation_ratio),
            "expected_ratio": "1/zeta",
        },
        "source_projection_scale_perturbation": {
            "status": "FIRED" if source_ablation_changed else "NOT_FIRED",
            "verdict": "RE_DERIVED_RESPONSE_CHANGED",
            "ratio": hstr(source_ablation_ratio),
            "expected_ratio": "sigma**2",
        },
        "projection_tensor_perturbation": {
            "status": "FIRED" if projection_ablation_changed else "NOT_FIRED",
            "verdict": "RE_DERIVED_RESPONSE_CHANGED",
            "delta": hstr(projection_ablation_delta),
        },
        "propagator_functional_perturbation_kminus4": {
            "status": "FIRED" if derived_tracks_functional_perturbation else "NOT_FIRED",
            "verdict": "RE_DERIVED_FALLOFF_CHANGED",
            "potential_power_before": potential_power,
            "potential_power_after": wrong_potential_power,
            "force_power_before": force_power,
            "force_power_after": wrong_force_power,
        },
        "target_readback_fixture": {
            "status": "FIRED" if target_readback["verdict"] == "FAIL_TARGET_READBACK" else "NOT_FIRED",
            "verdict": target_readback["verdict"],
            "fixture": "(s1*s2*Nu^2*aT^2/mu_R)*dot(V1 cross (V2 cross n), n)/(8*pi*R)",
            "fixture_reduced": hstr(target_readback_fixture),
            "functional_perturbation_delta": hstr(target_readback_delta),
        },
        "wrong_falloff_fixture": {
            "status": "FIRED" if wrong_falloff["verdict"] == "FAIL_WRONG_FALLOFF" else "NOT_FIRED",
            "verdict": wrong_falloff["verdict"],
            "derived_wrong_potential_power_R": wrong_potential_power,
            "derived_wrong_force_power_R": wrong_force_power,
        },
        "noncompact_source_fixture": {
            "status": "FIRED" if noncompact["verdict"] == "FAIL_WRONG_FALLOFF" else "NOT_FIRED",
            "verdict": noncompact["verdict"],
        },
        "V_equals_zero": {
            "status": "FIRED" if is_zero_expr(zero_velocity_U) else "NOT_FIRED",
            "verdict": "NO_MOVING_SOURCE",
            "U": hstr(zero_velocity_U),
        },
        "neutral_plus_minus_composite": {
            "status": "FIRED" if is_zero_expr(neutral_composite_U) else "NOT_FIRED",
            "verdict": "ZERO_MONOPOLE_CURRENT_SOURCE",
            "U_sum": hstr(neutral_composite_U),
        },
        "charge_flip_s_to_minus_s": {
            "status": "FIRED" if is_zero_expr(charge_flip_delta) else "NOT_FIRED",
            "verdict": "SOURCE_SIGN_FLIPS",
            "U_flipped_plus_U": hstr(charge_flip_delta),
        },
        "ghost_wrong_sign_transverse": {
            "status": "FIRED" if sp.simplify(transverse_like_parallel_side_force + transverse_ghost_side_force) == 0 else "NOT_FIRED",
            "verdict": "MU_R_TO_MINUS_MU_R_REDERIVED_FLIPS_ATTRACTION_TO_REPULSION",
            "stable_like_parallel_radial_coefficient": hstr(transverse_like_parallel_side_force),
            "ghost_like_parallel_radial_coefficient": hstr(transverse_ghost_side_force),
        },
    }

    exprs_for_agreement = {
        "B_eff": B_eff,
        "c_gamma_squared": c_gamma2,
        "O_T_static_scalar": O_T_static[0, 0],
        "O_L_static_scalar": O_L_static[0, 0],
        "G_T_static_coeff": G_T_static_coeff,
        "G_L_static_coeff": G_L_static_coeff,
        "scalar_green_kminus2": standard_kernels["scalar_green"],
        "biharmonic_seed_kminus4": standard_kernels["biharmonic_seed"],
        "kk_over_k4_contract": standard_kernels["kk"],
        "qA_stage1_no_sign": qA0,
        "qL_stage1_no_sign": qL0,
        "qA1": qA1,
        "qA2": qA2,
        "qL1": qL1,
        "qL2": qL2,
        "D_dot": Ddot,
        "A_radial_product": Aprod,
        "K_T_realspace": K_T,
        "K_L_realspace": K_L,
        "K_T_perturbed_kminus4": K_T_perturbed,
        "K_L_perturbed_kminus4": K_L_perturbed,
        "U_T_integrate": U_T_integrate,
        "U_L_integrate": U_L_integrate,
        "U_total_integrate": U_total_integrate,
        "U_T_eom_solve": U_T_eom,
        "U_L_eom_solve": U_L_eom,
        "U_total_eom_solve": U_total_eom,
        "eom_minus_integrate_residual": sp.factor(U_total_eom - U_total_integrate),
        "U_total_wrong_falloff": U_total_wrong_falloff,
        "F_T_x": sp.factor(F_T[0]),
        "F_T_y": sp.factor(F_T[1]),
        "F_T_z": sp.factor(F_T[2]),
        "F_L_x": sp.factor(F_L[0]),
        "F_L_y": sp.factor(F_L[1]),
        "F_L_z": sp.factor(F_L[2]),
        "F_total_x": sp.factor(F_total[0]),
        "F_total_y": sp.factor(F_total[1]),
        "F_total_z": sp.factor(F_total[2]),
        "F_wrong_falloff_x": sp.factor(F_total_wrong_falloff[0]),
        "F_wrong_falloff_y": sp.factor(F_total_wrong_falloff[1]),
        "F_wrong_falloff_z": sp.factor(F_total_wrong_falloff[2]),
        "U_side_by_side": U_side,
        "F_side_radial_coeff": F_side_radial_coeff,
        "scalar_admixture_ratio": scalar_admixture_ratio,
        "raw_qL2_over_qA2": raw_ratio,
        "raw_ratio_reference": raw_ratio_reference,
        "lorentz_residual": lorentz_residual,
        "lorentz_on_cone_lock": lorentz_on_cone_lock,
        "mu_ablation_ratio": mu_ablation_ratio,
        "source_ablation_ratio": source_ablation_ratio,
        "projection_ablation_delta": projection_ablation_delta,
        "target_readback_fixture": target_readback_fixture,
        "target_readback_delta": target_readback_delta,
        "zero_velocity_U": zero_velocity_U,
        "neutral_composite_U": neutral_composite_U,
        "charge_flip_delta": charge_flip_delta,
        "main_verdict_code": sp.Integer(main["verdict_code"]),
        "wrong_falloff_verdict_code": sp.Integer(wrong_falloff["verdict_code"]),
        "target_readback_verdict_code": sp.Integer(target_readback["verdict_code"]),
        "noncompact_verdict_code": sp.Integer(noncompact["verdict_code"]),
        "potential_power_R": sp.Integer(potential_power),
        "force_power_R": sp.Integer(force_power),
        "wrong_potential_power_R": sp.Integer(wrong_potential_power),
        "wrong_force_power_R": sp.Integer(wrong_force_power),
        "transverse_like_parallel_side_force": transverse_like_parallel_side_force,
        "longitudinal_like_parallel_side_force": longitudinal_like_parallel_side_force,
        "total_like_parallel_side_force": total_like_parallel_side_force,
        "transverse_ghost_side_force": transverse_ghost_side_force,
        "sign_like_charge_like_current": sp.Integer(-1),
        "sign_like_charge_opposite_current": sp.Integer(1),
        "sign_opposite_charge_like_current": sp.Integer(1),
        "sign_opposite_charge_opposite_current": sp.Integer(-1),
        "dimensional_ablations_fired": sp.Integer(4),
    }
    expr_digest = sha256_text(json.dumps({key: mma_expr(value) for key, value in exprs_for_agreement.items()}, sort_keys=True))

    return {
        "symbols": {
            "rhoBr": rhoBr,
            "muR": muR,
            "rhoB0": rhoB0,
            "chiC": chiC,
            "cE": cE,
            "Nu": Nu,
            "aT": aT,
            "aL": aL,
            "s1": s1,
            "s2": s2,
            "R": R,
            "zeta": zeta,
            "sigma": sigma,
            "k2": k2,
            "v1x": v1x,
            "v1y": v1y,
            "v1z": v1z,
            "v2x": v2x,
            "v2y": v2y,
            "v2z": v2z,
            "nx": nx,
            "ny": ny,
            "nz": nz,
            "vDot": vdot,
            "Ddot": Ddot,
            "arad": arad,
            "brad": brad,
            "Aprod": Aprod,
        },
        "B_eff": B_eff,
        "c_gamma2": c_gamma2,
        "qA0": qA0,
        "qL0": qL0,
        "D": Ddot,
        "a": arad,
        "b": brad,
        "A": Aprod,
        "operator_inversion": {
            "O_T_static_scalar": hstr(O_T_static[0, 0]),
            "O_L_static_scalar": hstr(O_L_static[0, 0]),
            "G_T_static_coeff": hstr(G_T_static_coeff),
            "G_L_static_coeff": hstr(G_L_static_coeff),
        },
        "inverse_fourier": {
            "scalar_green_kminus2": hstr(standard_kernels["scalar_green"]),
            "biharmonic_seed_kminus4": hstr(standard_kernels["biharmonic_seed"]),
            "kk_over_k4_contract": hstr(standard_kernels["kk"]),
            "scalar_green_kminus4": hstr(perturbed_kernels["scalar_green"]),
            "polyharmonic_seed_kminus6": hstr(perturbed_kernels["biharmonic_seed"]),
        },
        "K_T": K_T,
        "K_L": K_L,
        "K_T_perturbed": K_T_perturbed,
        "K_L_perturbed": K_L_perturbed,
        "U_T_integrate": U_T_integrate,
        "U_L_integrate": U_L_integrate,
        "U_total_integrate": U_total_integrate,
        "U_total_eom": U_total_eom,
        "U_total_wrong_falloff": U_total_wrong_falloff,
        "F_T": F_T,
        "F_L": F_L,
        "F_total": F_total,
        "F_total_wrong_falloff": F_total_wrong_falloff,
        "U_side": U_side,
        "F_side_radial_coeff": F_side_radial_coeff,
        "scalar_admixture_ratio": scalar_admixture_ratio,
        "raw_ratio": raw_ratio,
        "raw_ratio_reference": raw_ratio_reference,
        "main": main,
        "wrong_falloff": wrong_falloff,
        "noncompact": noncompact,
        "target_readback": target_readback,
        "self_tests": self_tests,
        "controls": controls,
        "lorentz": {
            "condition": "c_E = c_gamma",
            "residual": hstr(lorentz_residual),
            "on_condition_residual": hstr(lorentz_on_cone_lock),
            "diagnostic": "PREFERRED_FRAME_UNLESS_cE_EQUALS_cGAMMA",
            "is_fail": False,
        },
        "sign_diagnostic": {
            "scope": "side-by-side parallel throats with V_i dot Rhat = 0; outward radial sign >0 means repel",
            "raw_ratio": hstr(raw_ratio),
            "raw_ratio_reference": hstr(raw_ratio_reference),
            "scalar_admixture_ratio": hstr(scalar_admixture_ratio),
            "radial_sign_factor": "-sign(s1*s2)*sign(V1.V2) for positive transverse+longitudinal stiffnesses",
            "transverse_like_parallel_radial_coefficient": hstr(transverse_like_parallel_side_force),
            "longitudinal_like_parallel_radial_coefficient": hstr(longitudinal_like_parallel_side_force),
            "total_like_parallel_radial_coefficient": hstr(total_like_parallel_side_force),
            "landing": "NO_CANCELLATION_BOTH_CHANNELS_ATTRACTIVE",
            "table": additive_sign_matrix(),
        },
        "exprs_for_agreement": exprs_for_agreement,
        "expr_digest": expr_digest,
    }


def build_payload() -> dict[str, Any]:
    imports = import_banked_inputs()
    dims = build_dimensions()
    sym = build_symbolics(imports)
    agreement_exprs = {key: mma_expr(value) for key, value in sym["exprs_for_agreement"].items()}

    def control_label(item: dict[str, Any]) -> str:
        return str(item.get("verdict", item.get("status", "UNKNOWN")))

    agreement_payload = {
        "top_line_verdict": sym["main"]["verdict"],
        "main_verdict_code": sym["main"]["verdict_code"],
        "control_verdicts": {name: control_label(item) for name, item in sym["controls"].items()},
        "self_tests": sym["self_tests"],
        "checked_expression_count": len(sym["exprs_for_agreement"]),
        "expr_digest": sym["expr_digest"],
    }

    return {
        "schema": SCHEMA,
        "engine": "sympy",
        "directive": "software/stage1_solver/directives/pathA_39_magnetism_close_maxwell.md",
        "top_line_verdict": sym["main"]["verdict"],
        "engine_agreement": {
            "status": "PENDING_MATHEMATICA",
            "mathematica_exprs": agreement_exprs,
            "sympy_expression_digest": sym["expr_digest"],
        },
        "imports": imports,
        "sectors": {
            "transverse": {
                "field": "u_T",
                "lagrangian": "1/2 rho_br dot(u_T)^2 - 1/2 mu_R (curl u_T)^2",
                "source": "j_T=q_A^T V, q_A^T=Nu*aT*s",
                "operator_static_projector_eigenvalue": sym["operator_inversion"]["O_T_static_scalar"],
                "propagator_static_from_inversion": sym["operator_inversion"]["G_T_static_coeff"],
                "realspace_kernel_source": "inverse Fourier of (delta_ij-k_i k_j/k^2)/k^2",
                "energy_sign": "negative exchange energy for stable mu_R>0",
            },
            "longitudinal": {
                "field": "u_L",
                "lagrangian": "1/2 rho_br dot(u_L)^2 - 1/2 B_eff (div u_L)^2",
                "source": "j_L=q_L V, q_L=Nu*aL*s",
                "operator_static_projector_eigenvalue": sym["operator_inversion"]["O_L_static_scalar"],
                "propagator_static_from_inversion": sym["operator_inversion"]["G_L_static_coeff"],
                "realspace_kernel_source": "inverse Fourier of (k_i k_j/k^2)/k^2",
                "energy_sign": "negative exchange energy for stable B_eff>0",
            },
        },
        "derivation": {
            "operator_inversion": sym["operator_inversion"],
            "inverse_fourier_identities": sym["inverse_fourier"],
            "exchange_integral": "U_12=-int d^3k/(2*pi)^3 j1(-k).G(k).j2(k)",
            "green_tensor_consistency": "the derived Green tensor is reinserted into the static EOM as a consistency residual; no independent second-derivation claim is made",
        },
        "kernel": {
            "definitions": {
                "R": "X2-X1, R=|R|, n=R/R",
                "D": "V1.V2",
                "a": "V1.n",
                "b": "V2.n",
                "A": "(V1.n)*(V2.n)",
                "qA_i": "s_i*Nu*aT",
                "qL_i": "s_i*Nu*aL",
            },
            "compact": {
                "K_T": "(D + A)/(8*pi*R)",
                "K_L": "(D - A)/(8*pi*R)",
                "U_12": "-s1*s2*Nu^2/(8*pi*R)*[aT^2*(D+A)/mu_R + aL^2*(D-A)/B_eff]",
                "B_eff_substituted": "-s1*s2*Nu^2/(8*pi*R)*[aT^2*(D+A)/mu_R + aL^2*chi_c*(D-A)/rho_B0^2]",
                "F_12": "-grad_R U_12, derived symbolically from the computed U_12",
            },
            "K_T": hstr(sym["K_T"]),
            "K_L": hstr(sym["K_L"]),
            "K_T_perturbed_kminus4": hstr(sym["K_T_perturbed"]),
            "K_L_perturbed_kminus4": hstr(sym["K_L_perturbed"]),
            "U_T_integrate_out": hstr(sym["U_T_integrate"]),
            "U_L_integrate_out": hstr(sym["U_L_integrate"]),
            "U_total": hstr(sym["U_total_integrate"]),
            "U_eom_green_solve": hstr(sym["U_total_eom"]),
            "eom_minus_integrate": hstr(sp.factor(sym["U_total_eom"] - sym["U_total_integrate"])),
            "U_wrong_falloff_rederived": hstr(sym["U_total_wrong_falloff"]),
            "F_12": hstr(sym["F_total"]),
            "F_wrong_falloff_rederived": hstr(sym["F_total_wrong_falloff"]),
            "falloff": {
                "potential": f"R^{sym['main']['potential_power_R']}",
                "point_force": f"R^{sym['main']['force_power_R']}",
                "brane_localized": True,
                "classification": sym["main"],
                "measured_from": "computed U_total and F_total expressions",
            },
        },
        "sign_diagnostic": sym["sign_diagnostic"],
        "lorentz_form_diagnostic": sym["lorentz"],
        "classifier": {
            "rule": "PASS requires a perturbation-tracking derivation, then measured brane-localized 1/R potential and 1/R^2 point force",
            "labels": list(VERDICT_CODES.keys()),
            "main_inputs": sym["main"],
            "self_tests": sym["self_tests"],
        },
        "controls": sym["controls"],
        "scenarios": {
            "main": sym["main"],
            "wrong_falloff": sym["wrong_falloff"],
            "noncompact": sym["noncompact"],
            "target_readback": sym["target_readback"],
        },
        "dimensional_firewall": {
            "pass": dims["pass"],
            "checked_expression_count": dims["checked_expression_count"],
            "records": dims["records"],
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

    results = {
        "schema": SCHEMA,
        "verdict": payload["top_line_verdict"],
        "headline": {
            "verdict": payload["top_line_verdict"],
            "engine_agreement": "ENGINE_AGREE",
            "potential_falloff": payload["kernel"]["falloff"]["potential"],
            "point_force_falloff": payload["kernel"]["falloff"]["point_force"],
            "scalar_admixture_ratio": payload["sign_diagnostic"]["scalar_admixture_ratio"],
            "sign_landing": payload["sign_diagnostic"]["landing"],
            "lorentz_form": payload["lorentz_form_diagnostic"]["diagnostic"],
        },
        "engine_agreement": {
            "status": "ENGINE_AGREE",
            "compared_payload": payload["agreement_payload"],
            "sympy_payload": str(SYM_OUT),
            "mathematica_payload": str(WL_OUT),
        },
        "imports": payload["imports"],
        "sectors": payload["sectors"],
        "derivation": payload["derivation"],
        "kernel": payload["kernel"],
        "sign_diagnostic": payload["sign_diagnostic"],
        "lorentz_form_diagnostic": payload["lorentz_form_diagnostic"],
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
    kernel = results["kernel"]
    signs = results["sign_diagnostic"]
    controls = results["controls"]
    imports = results["imports"]["imported"]
    declared = results["imports"]["declared_vs_imported"]
    derivation = results["derivation"]

    lines = [
        "# pathA_39 Stage 2 Magnetic Force",
        "",
        f"Computed headline: `{results['verdict']}` with dual-engine `{results['engine_agreement']['status']}`.",
        "",
        "The Stage-2 interaction is now obtained by integrating out the transverse and longitudinal fields.  The `R` powers are measured from the resulting expressions, not supplied to the classifier.  The result is conditional on the declared Stage-1 moving-source amplitudes `q_A^T=Nu*aT*s` and `q_L=Nu*aL*s`; the values of `aT` and `aL` remain sim-deferred.",
        "",
        "## Derivation",
        "",
        "In the static projector eigenbasis, the quadratic operators and their inverses are:",
        "",
        "```text",
        f"O_T = {derivation['operator_inversion']['O_T_static_scalar']},  G_T = {derivation['operator_inversion']['G_T_static_coeff']}",
        f"O_L = {derivation['operator_inversion']['O_L_static_scalar']},  G_L = {derivation['operator_inversion']['G_L_static_coeff']}",
        "```",
        "",
        "The exchange integral used by both engines is:",
        "",
        "```text",
        derivation["exchange_integral"],
        f"F^-1[1/k^2] = {derivation['inverse_fourier_identities']['scalar_green_kminus2']}",
        f"F^-1[1/k^4] seed = {derivation['inverse_fourier_identities']['biharmonic_seed_kminus4']}",
        f"F^-1[k_i k_j/k^4] contracted with V1,V2 = {derivation['inverse_fourier_identities']['kk_over_k4_contract']}",
        "```",
        "",
        "## Kernel",
        "",
        "With `R=X2-X1`, `n=R/R`, `D=V1.V2`, and `A=(V1.n)(V2.n)`:",
        "",
        "```text",
        f"K_T compact = {kernel['compact']['K_T']}",
        f"K_L compact = {kernel['compact']['K_L']}",
        f"U_12 compact = {kernel['compact']['U_12']}",
        f"F_12 compact = {kernel['compact']['F_12']}",
        "",
        f"K_T = {kernel['K_T']}",
        f"K_L = {kernel['K_L']}",
        f"U_T = {kernel['U_T_integrate_out']}",
        f"U_L = {kernel['U_L_integrate_out']}",
        f"U_12 = {kernel['U_total']}",
        f"F_12 = {kernel['F_12']}",
        "```",
        "",
        f"Reinserting the same derived Green tensor into the static EOM gives the consistency residual `U_eom - U_integrate = {kernel['eom_minus_integrate']}`.",
        "",
        f"Measured falloff: potential `{kernel['falloff']['potential']}`, point force `{kernel['falloff']['point_force']}`.",
        "",
        "## Sign Diagnostic",
        "",
        f"For side-by-side parallel throats (`V_i.n=0`), the scalar-to-transverse magnitude ratio is `{signs['scalar_admixture_ratio']}`.  It is not a cancellation ratio: both stable channels have the same attractive sign for like sources.",
        "",
        f"Transverse like-current coefficient: `{signs['transverse_like_parallel_radial_coefficient']}`.",
        f"Longitudinal like-current coefficient: `{signs['longitudinal_like_parallel_radial_coefficient']}`.",
        f"Total like-current coefficient: `{signs['total_like_parallel_radial_coefficient']}`.",
        "",
        f"Landing: `{signs['landing']}`.  Outward radial force means repulsion; inward means attraction.",
        "",
        "| case | like charge, like currents | like charge, opposite currents | opposite charge, like currents | opposite charge, opposite currents |",
        "|---|---|---|---|---|",
    ]
    matrix = signs["table"]
    lines.append(
        f"| additive T+L | `{matrix['like_charge_s1s2_plus']['like_parallel_currents']}` | "
        f"`{matrix['like_charge_s1s2_plus']['opposite_parallel_currents']}` | "
        f"`{matrix['opposite_charge_s1s2_minus']['like_parallel_currents']}` | "
        f"`{matrix['opposite_charge_s1s2_minus']['opposite_parallel_currents']}` |"
    )
    lines.extend(
        [
            "",
            "The transverse `u_T` contribution gives the magnetic like-current attraction.  The longitudinal `u_L` contribution is an unavoidable attractive scalar-current admixture under the same stable-sign assumptions.",
            "",
            "## Lorentz-Form Diagnostic",
            "",
            f"`{results['lorentz_form_diagnostic']['diagnostic']}`.  The residual is `{results['lorentz_form_diagnostic']['residual']}`, so a Lorentz-force-form rewrite requires the additional condition `c_E=c_gamma`.  This is reported as a diagnostic, not a Stage-2 fail.",
            "",
            "## Controls",
            "",
            "| control | status | verdict/result |",
            "|---|---:|---|",
        ]
    )
    for name, item in controls.items():
        lines.append(f"| `{name}` | `{item['status']}` | `{item.get('verdict', '')}` |")
    lines.extend(
        [
            "",
            "## Dimensional Firewall",
            "",
            f"Units-restored checks passed for `{results['dimensional_firewall']['checked_expression_count']}` expressions.  The able-to-fail ablations fired for omitted velocity, omitted compact brane-source normalization, mass-current/charge-current mixing, and using `c_s` where the imported transverse `c_gamma` is required.",
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
        lines.append(f"- `{key}`: {'; '.join(values)}")
    lines.extend(
        [
            "",
            "## Dual Engine",
            "",
            f"`ENGINE_AGREE` over `{results['engine_agreement']['compared_payload']['checked_expression_count']}` compared quantities, including the inverted propagators, inverse-Fourier kernels, Green-tensor consistency residual, force components, measured falloff powers, corrected sign codes, dimensional ablation count, and target-readback controls.",
            "",
            "Run commands:",
            "",
            "```text",
            "timeout 600 python3 software/stage1_solver/tools/pathA_39_magnetic_force_sympy.py",
            "timeout 600 math -script software/stage1_solver/tools/pathA_39_magnetic_force.wl",
            "timeout 600 python3 software/stage1_solver/tools/pathA_39_magnetic_force_sympy.py --compare",
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
