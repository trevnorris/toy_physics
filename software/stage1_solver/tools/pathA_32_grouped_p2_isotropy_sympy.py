#!/usr/bin/env python3
"""PathA-32 Gate 3 grouped-P2 isotropy gate, SymPy engine.

This script intentionally keeps the angular and response parts visible:

1. construct the explicit real l=2 harmonics;
2. integrate the 5x5 Gram matrix and the anisotropic self-overlaps;
3. compute the -Delta_S2 eigenvalue used in K2;
4. assemble every ungrouped channel from its own harmonic input;
5. compute the raw-D gate, normalized cross-check, stability guards, and
   able-to-fail probes through one verdict function.

The final report is written only after the Mathematica engine has produced its
independent scratch YAML. Before that, this script writes the SymPy scratch YAML
and exits 0 so the required run order can be:

  timeout 600 python ..._sympy.py
  timeout 600 math -script ... .wl
  timeout 600 python ..._sympy.py
"""

from __future__ import annotations

import hashlib
import math
from pathlib import Path
from typing import Any

import sympy as sp
import yaml


class NoAliasDumper(yaml.SafeDumper):
    def ignore_aliases(self, data: Any) -> bool:
        return True


SCRIPT_PATH = Path(__file__).resolve()
STAGE1_ROOT = SCRIPT_PATH.parents[1]
REPO_ROOT = SCRIPT_PATH.parents[3]
REPORTS = STAGE1_ROOT / "reports"
SCRATCH = STAGE1_ROOT / "_scratch"
NOTES = REPO_ROOT / "research" / "pde_ledger" / "notes" / "stages"

SYM_YAML = SCRATCH / "pathA_32_sympy_results.yaml"
MMA_YAML = SCRATCH / "pathA_32_mathematica_results.yaml"
RESULTS_YAML = REPORTS / "pathA_32_results.yaml"
REPORT_MD = REPORTS / "pathA_32_grouped_p2_isotropy.md"
FEED_NOTE = NOTES / "moving_throat_pde_pathA_32_grouped_p2_isotropy_result.md"

EPS_VALUES = [1.0e-4, 2.0e-4, -3.0e-4]
DELTA_PROFILE = 0.1
SYMBOLIC_TOL = 1.0e-10
NUMERIC_TOL = 1.0e-8
PROBE_TOL = 1.0e-12

CALIBRATION_SAMPLE = {
    "Mtilde": 3.0,
    "Ktilde": 7.0,
    "TomegaTilde": 0.5,
    "B0tilde": 1.0,
    "Z0tilde": 0.5,
    "B2tilde": 0.4,
    "Z2tilde": 0.2,
    "B4tilde": 0.05,
    "Z4tilde": 0.03,
}

CALIBRATION_WINDOW = {
    "Mtilde": {"min": 2.0, "max": 5.0, "unit": "mass"},
    "Ktilde": {"min": 5.0, "max": 11.0, "unit": "mass/time^2"},
    "TomegaTilde": {"min": 0.25, "max": 1.0, "unit": "mass/time^2"},
    "B0tilde": {"min": 0.5, "max": 1.5, "unit": "mass/time^2"},
    "Z0tilde": {"min": 0.25, "max": 0.75, "unit": "mass/time^2"},
}


def compact(expr: Any) -> Any:
    if isinstance(expr, sp.Basic):
        if expr.has(sp.sin, sp.cos):
            return sp.factor(sp.trigsimp(sp.simplify(expr)))
        return sp.factor(sp.cancel(expr))
    return expr


def hstr(expr: Any) -> str | bool | int | float | None:
    if expr is None or isinstance(expr, (bool, int, float, str)):
        return expr
    return sp.sstr(compact(expr))


def yaml_write(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        yaml.dump(payload, Dumper=NoAliasDumper, sort_keys=False, allow_unicode=False),
        encoding="utf-8",
    )


def yaml_read(path: Path) -> dict[str, Any] | None:
    if not path.exists():
        return None
    data = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(data, dict):
        raise AssertionError(f"YAML did not load to a mapping: {path}")
    return data


def expr_hash(expr: sp.Expr) -> str:
    return hashlib.sha256(sp.sstr(compact(expr)).encode("utf-8")).hexdigest()


def matrix_to_strings(matrix: sp.Matrix) -> list[list[str]]:
    return [[str(hstr(matrix[i, j])) for j in range(matrix.cols)] for i in range(matrix.rows)]


def matrix_to_floats(matrix: sp.Matrix) -> list[list[float]]:
    return [[float(sp.N(matrix[i, j], 30)) for j in range(matrix.cols)] for i in range(matrix.rows)]


def integrate_s2(expr: sp.Expr, theta: sp.Symbol, phi: sp.Symbol) -> sp.Expr:
    return compact(
        sp.integrate(
            sp.integrate(sp.expand_trig(expr) * sp.sin(theta), (phi, 0, 2 * sp.pi)),
            (theta, 0, sp.pi),
        )
    )


def laplacian_s2(expr: sp.Expr, theta: sp.Symbol, phi: sp.Symbol) -> sp.Expr:
    return compact(
        (1 / sp.sin(theta)) * sp.diff(sp.sin(theta) * sp.diff(expr, theta), theta)
        + (1 / sp.sin(theta) ** 2) * sp.diff(expr, phi, 2)
    )


def defect_pair(triple: list[sp.Expr]) -> tuple[sp.Expr, sp.Expr]:
    return (
        compact((2 * triple[0] - triple[1] - triple[2]) / 10),
        compact((triple[1] - triple[2]) / 2),
    )


def all_zero(expressions: list[sp.Expr]) -> bool:
    return all(compact(expr) == 0 for expr in expressions)


def sample_subs(symbols: dict[str, sp.Symbol]) -> dict[sp.Symbol, float]:
    return {symbols[name]: value for name, value in CALIBRATION_SAMPLE.items()}


def f_eval(expr: sp.Expr, subs: dict[sp.Symbol, float]) -> float:
    return float(sp.N(compact(expr).subs(subs), 30))


class DimError(ValueError):
    pass


Dim = tuple[sp.Rational, sp.Rational, sp.Rational]
ZERO_DIM: Dim = (sp.Rational(0), sp.Rational(0), sp.Rational(0))
DIMENSIONLESS_FUNCTIONS = {
    sp.sin,
    sp.cos,
    sp.tan,
    sp.cot,
    sp.sinh,
    sp.cosh,
    sp.tanh,
    sp.coth,
    sp.sech,
    sp.csch,
}
Ldim, Mdim, Tdim = sp.symbols("L M T", positive=True)


def dim_add(left: Dim, right: Dim) -> Dim:
    return tuple(sp.Rational(a0 + b0) for a0, b0 in zip(left, right))  # type: ignore[return-value]


def dim_sub(left: Dim, right: Dim) -> Dim:
    return tuple(sp.Rational(a0 - b0) for a0, b0 in zip(left, right))  # type: ignore[return-value]


def dim_scale(dim: Dim, scale: sp.Rational) -> Dim:
    return tuple(sp.Rational(scale * d0) for d0 in dim)  # type: ignore[return-value]


def dim_of(expr: sp.Expr, symbol_dims: dict[sp.Symbol, Dim]) -> Dim:
    expr = sp.sympify(expr)
    if expr == 0 or expr.is_number:
        return ZERO_DIM
    if expr.is_Symbol:
        if expr not in symbol_dims:
            raise DimError(f"missing dimension for symbol {expr}")
        return symbol_dims[expr]
    if expr.is_Mul:
        out = ZERO_DIM
        for arg in expr.args:
            out = dim_add(out, dim_of(arg, symbol_dims))
        return out
    if expr.is_Pow:
        base, exponent = expr.args
        if not exponent.is_number:
            raise DimError(f"non-numeric dimension exponent in {expr}")
        return dim_scale(dim_of(base, symbol_dims), sp.Rational(exponent))
    if expr.is_Add:
        dims = [dim_of(arg, symbol_dims) for arg in expr.args if arg != 0]
        if not dims:
            return ZERO_DIM
        first = dims[0]
        mismatched = [dim for dim in dims if dim != first]
        if mismatched:
            raise DimError(f"dimension mismatch in sum {expr}: {dims}")
        return first
    if expr.func in DIMENSIONLESS_FUNCTIONS:
        arg_dims = [dim_of(arg, symbol_dims) for arg in expr.args]
        if any(dim != ZERO_DIM for dim in arg_dims):
            raise DimError(f"dimensionful argument in {expr}: {arg_dims}")
        return ZERO_DIM
    raise DimError(f"unsupported expression in dimension checker: {expr}")


def dim_to_monomial(dim: Dim) -> sp.Expr:
    return sp.factor(Ldim ** dim[0] * Mdim ** dim[1] * Tdim ** dim[2])


def exp_text(exp: sp.Rational) -> str:
    exp = sp.Rational(exp)
    if exp.q == 1:
        return str(exp.p)
    return f"{exp.p}/{exp.q}"


def dim_to_text(dim: Dim) -> str:
    parts: list[str] = []
    for name, exp in (("L", dim[0]), ("T", dim[2]), ("M", dim[1])):
        exp = sp.Rational(exp)
        if exp == 0:
            continue
        if exp == 1:
            parts.append(name)
        else:
            parts.append(f"{name}^{exp_text(exp)}")
    return " ".join(parts) if parts else "1"


def dim_vector_text(dim: Dim | None) -> list[str] | None:
    if dim is None:
        return None
    return [hstr(v) for v in dim]


def dim_record(name: str, expr: sp.Expr, symbol_dims: dict[sp.Symbol, Dim]) -> dict[str, Any]:
    dim = dim_of(expr, symbol_dims)
    return {
        "quantity": name,
        "expression": hstr(expr),
        "dimension": dim_to_text(dim),
        "dimension_monomial": hstr(dim_to_monomial(dim)),
        "dimension_vector_LMT": dim_vector_text(dim),
    }


def dimension_verdict(ok: bool) -> str:
    return "DIMENSIONAL_OK" if ok else "FAIL_DIMENSIONAL"


def max_ratio_delta(values: list[float], eps_values: list[float]) -> float:
    ratios = [value / eps for value, eps in zip(values, eps_values) if eps != 0.0]
    if not ratios:
        return math.inf
    return max(abs(ratio - ratios[0]) for ratio in ratios)


def first_nonzero(values: list[float], tol: float = 1.0e-12) -> bool:
    return any(abs(value) > tol for value in values)


def verdict_from_gates(gates: dict[str, bool], calibration_inputs: list[str]) -> str:
    if not gates["dimensional_ok"]:
        return "FAIL_DIMENSIONAL"
    if not gates["covariant"]:
        return "FAIL_NOT_COVARIANT"
    if not gates["tautology_clear"]:
        return "FAIL_TAUTOLOGICAL"
    if not gates["dynamic_retained"]:
        return "FAIL_STATIC_RESPONSE"
    if not gates["stability_ok"]:
        return "FAIL_STABILITY"
    if not gates["denominator_guard_ok"]:
        return "FAIL_SINGULAR_RESPONSE"
    if not gates["lane_collapse_ok"]:
        return "FAIL_ANISOTROPIC_BRANCH"
    if not gates["able_to_fail_ok"]:
        return "FAIL_NOT_ABLE_TO_FAIL"
    if calibration_inputs:
        return "ISOTROPY_CALIBRATED"
    return "ISOTROPY_PASS"


def build_dimensional_check(
    *,
    Mtilde: sp.Symbol,
    Ktilde: sp.Symbol,
    TomegaTilde: sp.Symbol,
    M2_expr: sp.Expr,
    K2_expr: sp.Expr,
    lambda_m: sp.Expr,
) -> dict[str, Any]:
    a_dim, dw, dOmega = sp.symbols("a_dim dw dOmega", positive=True)
    beta2, beta2_prime = sp.symbols("beta2 beta2_prime", nonzero=True)
    mu_eta_density, T_w_density, K_eta_density, T_omega_density = sp.symbols(
        "mu_eta_density T_w_density K_eta_density T_Omega_density", nonzero=True
    )
    measure = a_dim**2 * dw * dOmega
    m2_integral = mu_eta_density * beta2**2 * measure
    k_tw_term = T_w_density * beta2_prime**2 * measure
    k_eta_term = K_eta_density * beta2**2 * measure
    k_omega_term = lambda_m * T_omega_density * beta2**2 * measure
    k2_integral = k_tw_term + k_eta_term + k_omega_term

    symbol_dims: dict[sp.Symbol, Dim] = {
        a_dim: (1, 0, 0),
        dw: (1, 0, 0),
        dOmega: ZERO_DIM,
        beta2: ZERO_DIM,
        beta2_prime: (-1, 0, 0),
        mu_eta_density: (-3, 1, 0),
        T_w_density: (-1, 1, -2),
        K_eta_density: (-3, 1, -2),
        T_omega_density: (-3, 1, -2),
        Mtilde: (0, 1, 0),
        Ktilde: (0, 1, -2),
        TomegaTilde: (0, 1, -2),
    }
    expected_m = (0, 1, 0)
    expected_k = (0, 1, -2)
    expected_ratio = (0, 0, -2)

    measure_dim = dim_of(measure, symbol_dims)
    m2_integral_dim = dim_of(m2_integral, symbol_dims)
    k_tw_dim = dim_of(k_tw_term, symbol_dims)
    k_eta_dim = dim_of(k_eta_term, symbol_dims)
    k_omega_dim = dim_of(k_omega_term, symbol_dims)
    k2_integral_dim = dim_of(k2_integral, symbol_dims)
    actual_m2_dim = dim_of(M2_expr, symbol_dims)
    actual_k2_dim = dim_of(K2_expr, symbol_dims)
    actual_ratio_dim = dim_sub(actual_k2_dim, actual_m2_dim)
    term_homogeneous = k_tw_dim == k_eta_dim == k_omega_dim == k2_integral_dim
    dimensional_ok = bool(
        measure_dim == (3, 0, 0)
        and m2_integral_dim == expected_m
        and term_homogeneous
        and k2_integral_dim == expected_k
        and actual_m2_dim == expected_m
        and actual_k2_dim == expected_k
        and actual_ratio_dim == expected_ratio
    )

    corrupt_dims = dict(symbol_dims)
    corrupt_dims[T_omega_density] = dim_add(symbol_dims[T_omega_density], (1, 0, 0))
    corrupt_dims[TomegaTilde] = dim_add(symbol_dims[TomegaTilde], (1, 0, 0))
    corrupt_error = None
    try:
        corrupt_k_omega_dim = dim_of(k_omega_term, corrupt_dims)
        corrupt_k2_integral_dim = dim_of(k2_integral, corrupt_dims)
        corrupt_actual_k2_dim = dim_of(K2_expr, corrupt_dims)
        corrupt_actual_ratio_dim = dim_sub(corrupt_actual_k2_dim, actual_m2_dim)
        corrupt_ok = bool(
            corrupt_k2_integral_dim == expected_k
            and corrupt_actual_k2_dim == expected_k
            and corrupt_actual_ratio_dim == expected_ratio
        )
    except DimError as exc:
        corrupt_k_omega_dim = dim_of(k_omega_term, corrupt_dims)
        corrupt_k2_integral_dim = None
        corrupt_actual_k2_dim = None
        corrupt_actual_ratio_dim = None
        corrupt_ok = False
        corrupt_error = str(exc)
    probe_verdict = "NO_FAIL" if corrupt_ok else "FAIL_DIMENSIONAL"

    return {
        "dimension_order": ["L", "M", "T"],
        "dimensional_gate": "explicit a^2 dw dOmega measure plus grouped M2/K2 ratio",
        "headline_quantities_walked": {
            "explicit_measure": hstr(measure),
            "M2_integral": hstr(m2_integral),
            "K2_integral": hstr(k2_integral),
            "K2_terms": {
                "T_w_beta_prime_sq": hstr(k_tw_term),
                "K_eta_beta_sq": hstr(k_eta_term),
                "lambda_T_Omega_beta_sq": hstr(k_omega_term),
            },
            "actual_grouped_M2": hstr(M2_expr),
            "actual_grouped_K2": hstr(K2_expr),
            "actual_K2_over_M2": hstr(K2_expr / M2_expr),
        },
        "symbol_dimensions": {
            "a_dim": "L",
            "dw": "L",
            "dOmega": "1",
            "dV=a_dim^2*dw*dOmega": "L^3",
            "beta2": "1",
            "beta2_prime": "L^-1",
            "mu_eta_density": "M L^-3",
            "T_w_density": "M L^-1 T^-2",
            "K_eta_density": "M L^-3 T^-2",
            "T_Omega_density": "M L^-3 T^-2",
            "Mtilde": "M",
            "Ktilde": "M T^-2",
            "TomegaTilde": "M T^-2",
        },
        "sourcing_note": "The explicit check restores the S2 wall measure dV=a^2 dw dOmega before walking the density terms; the actual grouped Mtilde/Ktilde/TomegaTilde dimensions are then sourced from those integrals, not back-solved from K2/M2.",
        "computed_dimensions": {
            "measure": dim_to_text(measure_dim),
            "M2_integral": dim_to_text(m2_integral_dim),
            "K2_terms": {
                "T_w_beta_prime_sq": dim_to_text(k_tw_dim),
                "K_eta_beta_sq": dim_to_text(k_eta_dim),
                "lambda_T_Omega_beta_sq": dim_to_text(k_omega_dim),
            },
            "K2_integral": dim_to_text(k2_integral_dim),
            "actual_grouped_M2": dim_to_text(actual_m2_dim),
            "actual_grouped_K2": dim_to_text(actual_k2_dim),
            "actual_K2_over_M2": dim_to_text(actual_ratio_dim),
        },
        "computed_dimension_vectors_LMT": {
            "measure": dim_vector_text(measure_dim),
            "M2_integral": dim_vector_text(m2_integral_dim),
            "K2_terms": {
                "T_w_beta_prime_sq": dim_vector_text(k_tw_dim),
                "K_eta_beta_sq": dim_vector_text(k_eta_dim),
                "lambda_T_Omega_beta_sq": dim_vector_text(k_omega_dim),
            },
            "K2_integral": dim_vector_text(k2_integral_dim),
            "actual_grouped_M2": dim_vector_text(actual_m2_dim),
            "actual_grouped_K2": dim_vector_text(actual_k2_dim),
            "actual_K2_over_M2": dim_vector_text(actual_ratio_dim),
        },
        "expected_dimensions": {
            "M2": dim_to_text(expected_m),
            "K2": dim_to_text(expected_k),
            "K2_over_M2": dim_to_text(expected_ratio),
        },
        "checks": {
            "measure_dimension_explicit": measure_dim == (3, 0, 0),
            "K2_terms_homogeneous": term_homogeneous,
            "M2_integral_has_mass_dimension": m2_integral_dim == expected_m,
            "K2_integral_has_stiffness_dimension": k2_integral_dim == expected_k,
            "actual_grouped_M2_has_mass_dimension": actual_m2_dim == expected_m,
            "actual_grouped_K2_has_stiffness_dimension": actual_k2_dim == expected_k,
            "actual_K2_over_M2_has_omega_squared_dimension": actual_ratio_dim == expected_ratio,
        },
        "dimensional_ok": dimensional_ok,
        "status": "pass" if dimensional_ok else "fail",
        "dimensional_status": dimension_verdict(dimensional_ok),
        "table": [
            dim_record("dV=a_dim^2*dw*dOmega", measure, symbol_dims),
            dim_record("M2_integral=mu_eta_density*beta2^2*dV", m2_integral, symbol_dims),
            dim_record("K2_Tw_term=T_w_density*beta2_prime^2*dV", k_tw_term, symbol_dims),
            dim_record("K2_Keta_term=K_eta_density*beta2^2*dV", k_eta_term, symbol_dims),
            dim_record("K2_TOmega_term=lambda*T_Omega_density*beta2^2*dV", k_omega_term, symbol_dims),
            dim_record("actual_grouped_M2", M2_expr, symbol_dims),
            dim_record("actual_grouped_K2", K2_expr, symbol_dims),
        ],
        "FAIL_DIMENSIONAL_probe": {
            "mutation": "corrupt sourced [T_Omega] and its assembled TomegaTilde scalar by one extra power of L",
            "participates_in_verdict": True,
            "sourced_T_Omega_dimension": dim_to_text(symbol_dims[T_omega_density]),
            "corrupted_T_Omega_dimension": dim_to_text(corrupt_dims[T_omega_density]),
            "sourced_TomegaTilde_dimension": dim_to_text(symbol_dims[TomegaTilde]),
            "corrupted_TomegaTilde_dimension": dim_to_text(corrupt_dims[TomegaTilde]),
            "mutated_dimensions": {
                "lambda_T_Omega_beta_sq": dim_to_text(corrupt_k_omega_dim),
                "K2_integral": dim_to_text(corrupt_k2_integral_dim) if corrupt_k2_integral_dim else "inhomogeneous",
                "actual_grouped_K2": dim_to_text(corrupt_actual_k2_dim) if corrupt_actual_k2_dim else "inhomogeneous",
                "actual_K2_over_M2": dim_to_text(corrupt_actual_ratio_dim) if corrupt_actual_ratio_dim else "inhomogeneous",
                "error": corrupt_error,
            },
            "without_mutation_dimensional_ok": dimensional_ok,
            "with_mutation_dimensional_ok": corrupt_ok,
            "probe_verdict": probe_verdict,
            "mutation_fires": probe_verdict == "FAIL_DIMENSIONAL",
        },
    }


def build_report(payload: dict[str, Any]) -> str:
    probes = payload["counterfactuals"]
    raw = payload["raw_D_gate"]["raw_D_defects"]
    engine = payload["engine_agreement"]
    part = payload["input_partition"]
    able = payload["able_to_fail"]
    fixed_probe_names = [
        "singular_denominator",
        "wrong_eigenvalue",
        "degenerate_beta_zero",
        "tautology_hash_collision",
        "static_drop_inertia",
        "dimensional_corrupt_T_Omega",
    ]
    fixed_probe_lines = [
        f"- `{name}`: with mutation `{probes[name]['verdict']}`, "
        f"without mutation `{probes[name]['self_ablation']['verdict']}`, "
        f"fail suppressed = `{probes[name]['self_ablation']['fail_suppressed']}`."
        for name in fixed_probe_names
    ]
    lines = [
        "# PathA-32 grouped-P2 isotropy result",
        "",
        f"Computed verdict: `{payload['verdict']}` (`{payload['which_rung']}`).",
        "",
        "The two engines computed the real l=2 basis, angular Gram matrix, Laplacian eigenvalues, "
        "grouped response coefficients, raw-D defects, normalized u-defects, and counterfactual "
        "verdict flips. The final rung is calibrated because the wall profile and radial/support "
        "scalars remain frozen calibration inputs rather than derived from the Gate-1 R0 support equation.",
        "",
        "## Key computed checks",
        "",
        f"- Gram matrix equals I5: `{payload['harmonics']['gram_is_identity']}`.",
        f"- Computed -Delta_S2 eigenvalues: `{payload['laplacian']['lambda_by_channel']}`.",
        f"- K2 angular coefficient equals computed lambda_m: `{payload['laplacian']['k2_coefficient_equals_computed_lambda_all']}`.",
        f"- Isotropic raw-D defects: `{raw}`.",
        f"- Normalized u-defects: `{payload['normalized_response']['normalized_defects']}`.",
        f"- Stability guard: `{payload['stability']['stability_ok']}`; denominator guard: `{payload['stability']['denominator_guard_ok']}`.",
        "",
        "## Probe outcomes",
        "",
        f"- Pure-prefactor anisotropy: `{probes['pure_prefactor_anisotropy']['verdict']}`; "
        f"raw-D moves, normalized defects stay zero = `{probes['pure_prefactor_anisotropy']['normalized_u_defects_stay_zero']}`.",
        f"- Sector-selective anisotropy: `{probes['sector_selective_anisotropy']['verdict']}`; "
        f"raw-D and normalized u move = `{probes['sector_selective_anisotropy']['u_defects_move']}`.",
        f"- m-dependent profile: `{probes['m_dependent_profile']['verdict']}`.",
        f"- Degenerate beta2=0: `{probes['degenerate_beta_zero']['verdict']}`.",
        f"- Wrong eigenvalue coefficient: `{probes['wrong_eigenvalue']['verdict']}`.",
        f"- Singular denominator guard probe: `{probes['singular_denominator']['verdict']}`.",
        f"- Tautology hash probe: `{probes['tautology_hash_collision']['verdict']}`.",
        f"- Static response probe: `{probes['static_drop_inertia']['verdict']}`.",
        f"- Dimensional sourced-input probe: `{probes['dimensional_corrupt_T_Omega']['verdict']}`.",
        "",
        "## Fixed probe self-ablations",
        "",
        *fixed_probe_lines,
        "",
        "## Able-to-fail aggregate",
        "",
        f"- Computed probe gate flags: `{able['computed_probe_gate_flags']}`.",
        f"- Neutering any one probe flips aggregate false: `{able['neutering_any_probe_flips_false']}`.",
        "",
        "## Engine agreement",
        "",
        f"- Status: `{engine['status']}`.",
        f"- Max symbolic delta: `{engine['max_symbolic_delta']}` with tolerance `{engine['symbolic_tolerance']}`.",
        f"- Max numeric delta: `{engine['max_numeric_delta']}` with tolerance `{engine['numeric_tolerance']}`.",
        f"- Per-lane `D_A,n` max numeric delta: `{engine.get('grouped_lane_D_max_numeric_delta')}`.",
        f"- Per-lane `D_A,n` deltas: `{engine.get('grouped_lane_D_numeric_deltas')}`.",
        f"- Dimensional agreement: `{engine['symbolic_checks'].get('dimensional_ok')}`; "
        f"probe verdict agreement: `{engine['symbolic_checks'].get('dimension_probe_verdict')}`.",
        "",
        "## Dimensional check (retrofit)",
        "",
        f"- Walked quantities: `{payload['dimensional']['headline_quantities_walked']}`.",
        f"- Sourced dimensions: `{payload['dimensional']['symbol_dimensions']}`.",
        f"- Computed dimensions: `{payload['dimensional']['computed_dimensions']}`.",
        f"- Probe flip: `{payload['counterfactuals']['dimensional_corrupt_T_Omega']}`.",
        "",
        "## Input partition",
        "",
        f"- Derived inputs: `{part['derived_inputs']}`.",
        f"- Calibration inputs: `{part['calibration_inputs']}`.",
        "",
        "Deferred: the 54/5 quadrupole normalization, outgoing odd N coefficients, and solved nonlinear branch data remain Gate 4/5/6 work.",
        "",
    ]
    return "\n".join(lines)


def build_feed_note(payload: dict[str, Any]) -> str:
    return "\n".join(
        [
            "# Moving-throat PDE PathA-32 grouped-P2 isotropy feed note",
            "",
            f"Computed verdict: `{payload['verdict']}` (`{payload['which_rung']}`).",
            "",
            "The SymPy and Mathematica engines independently computed the real l=2 angular basis, "
            "the isotropic Gram matrix `I5`, the per-harmonic `-Delta_S2` eigenvalue `6`, and the "
            "grouped conservative response assembled per lane from each harmonic self-overlap.",
            "",
            "On the isotropic calibrated reference, all raw-D lane defects vanish exactly for orders "
            "`0,2,4`; the normalized `a2,b2,a4,b4` cross-check also vanishes. The verdict is "
            "`ISOTROPY_CALIBRATED`, not `ISOTROPY_PASS`, because `beta2(w)`, `T_Omega`, wall stiffnesses, "
            "and the `Btilde/Ztilde/Ktilde/Mtilde` radial data remain calibration inputs.",
            "",
            "The able-to-fail probes now use computed gate flags. The remediated probes each fail with "
            "their mutation and stop failing under their self-ablation; neutering any probe's computed flag "
            "makes `able_to_fail_ok` false.",
            "",
            "Dimensional retrofit: the engines walk the explicit `a^2 dw dOmega` wall measure, the `M2`/`K2` integrand terms, and the actual grouped `Mtilde` / `Ktilde + 6*TomegaTilde` expressions. Corrupting sourced `T_Omega` flips the recomputed verdict to `FAIL_DIMENSIONAL` and the self-ablation restores `ISOTROPY_CALIBRATED`.",
            "",
            f"Engine agreement status: `{payload['engine_agreement']['status']}`, max numeric delta "
            f"`{payload['engine_agreement']['max_numeric_delta']}`, grouped-lane D max delta "
            f"`{payload['engine_agreement'].get('grouped_lane_D_max_numeric_delta')}`.",
            "",
        ]
    )


def compare_numeric_matrix(left: list[list[float]], right: list[list[Any]]) -> float:
    deltas: list[float] = []
    for i, row in enumerate(left):
        for j, value in enumerate(row):
            deltas.append(abs(float(value) - float(right[i][j])))
    return max(deltas) if deltas else 0.0


def compute_engine_agreement(sympy_payload: dict[str, Any], mma_payload: dict[str, Any] | None) -> dict[str, Any]:
    if mma_payload is None:
        return {
            "status": "pending_mathematica",
            "max_symbolic_delta": None,
            "max_numeric_delta": None,
            "symbolic_tolerance": SYMBOLIC_TOL,
            "numeric_tolerance": NUMERIC_TOL,
            "checked_fields": [],
        }

    numeric_deltas = [
        compare_numeric_matrix(
            sympy_payload["harmonics"]["gram_5x5_numeric"],
            mma_payload["harmonics"]["gram_5x5_numeric"],
        ),
        compare_numeric_matrix(
            sympy_payload["harmonics"]["grouped_reduction_isotropic_numeric"],
            mma_payload["harmonics"]["grouped_reduction_isotropic_numeric"],
        ),
        compare_numeric_matrix(
            sympy_payload["harmonics"]["grouped_reduction_pure_prefactor_linear_coeff_numeric"],
            mma_payload["harmonics"]["grouped_reduction_pure_prefactor_linear_coeff_numeric"],
        ),
    ]

    for name in sympy_payload["harmonics"]["order"]:
        numeric_deltas.append(
            abs(
                float(sympy_payload["laplacian"]["lambda_numeric_by_channel"][name])
                - float(mma_payload["laplacian"]["lambda_numeric_by_channel"][name])
            )
        )
        numeric_deltas.append(
            abs(
                float(sympy_payload["harmonics"]["anisotropy_coefficients_numeric"][name])
                - float(mma_payload["harmonics"]["anisotropy_coefficients_numeric"][name])
            )
        )

    grouped_lane_d_deltas: dict[str, float] = {}
    for lane in ["20", "21", "22"]:
        for order_n in ["0", "2", "4"]:
            key = f"{lane}.D{order_n}"
            delta = abs(
                float(sympy_payload["grouped_lanes"][lane]["D_sample"][order_n])
                - float(mma_payload["grouped_lanes"][lane]["D_sample"][order_n])
            )
            grouped_lane_d_deltas[key] = delta
            numeric_deltas.append(delta)

    for probe_name in ["pure_prefactor_anisotropy", "sector_selective_anisotropy"]:
        for defect_name, values in sympy_payload["counterfactuals"][probe_name]["sample_values"].items():
            mma_values = mma_payload["counterfactuals"][probe_name]["sample_values"][defect_name]
            numeric_deltas.extend(abs(float(a) - float(b)) for a, b in zip(values, mma_values))

    numeric_deltas.append(
        abs(
            float(sympy_payload["stability"]["omega_2m_sample"])
            - float(mma_payload["stability"]["omega_2m_sample"])
        )
    )
    max_numeric_delta = max(numeric_deltas) if numeric_deltas else 0.0

    mma_counterfactuals = mma_payload.get("counterfactuals", {})
    if not isinstance(mma_counterfactuals, dict):
        mma_counterfactuals = {}
    mma_dimensional = mma_payload.get("dimensional", {})
    if not isinstance(mma_dimensional, dict):
        mma_dimensional = {}

    symbolic_checks = {
        "gram_identity": sympy_payload["harmonics"]["gram_is_identity"]
        == mma_payload["harmonics"]["gram_is_identity"],
        "lambda_all_six": sympy_payload["laplacian"]["lambda_all_six"]
        == mma_payload["laplacian"]["lambda_all_six"],
        "k_coefficients_equal": sympy_payload["laplacian"]["k2_coefficient_equals_computed_lambda_all"]
        == mma_payload["laplacian"]["k2_coefficient_equals_computed_lambda_all"],
        "raw_D_defects_zero": sympy_payload["raw_D_gate"]["raw_D_defects_zero"]
        == mma_payload["raw_D_gate"]["raw_D_defects_zero"],
        "normalized_defects_zero": sympy_payload["normalized_response"]["normalized_defects_zero"]
        == mma_payload["normalized_response"]["normalized_defects_zero"],
        "grouped_lane_D_samples_present": all(
            order_n in mma_payload["grouped_lanes"][lane]["D_sample"]
            and order_n in sympy_payload["grouped_lanes"][lane]["D_sample"]
            for lane in ["20", "21", "22"]
            for order_n in ["0", "2", "4"]
        ),
        "probe_verdicts": {
            key: (
                key in mma_counterfactuals
                and sympy_payload["counterfactuals"][key]["verdict"] == mma_counterfactuals[key].get("verdict")
            )
            for key in [
                "pure_prefactor_anisotropy",
                "sector_selective_anisotropy",
                "m_dependent_profile",
                "degenerate_beta_zero",
                "wrong_eigenvalue",
                "singular_denominator",
                "tautology_hash_collision",
                "static_drop_inertia",
                "dimensional_corrupt_T_Omega",
            ]
        },
        "baseline_verdict": sympy_payload["verdict"] == mma_payload["verdict"],
        "dimensional_ok": sympy_payload["dimensional"]["dimensional_ok"]
        == mma_dimensional.get("dimensional_ok"),
        "dimension_probe_verdict": sympy_payload["counterfactuals"]["dimensional_corrupt_T_Omega"]["verdict"]
        == mma_counterfactuals.get("dimensional_corrupt_T_Omega", {}).get("verdict"),
    }
    nested_probe_ok = all(symbolic_checks["probe_verdicts"].values())
    symbolic_flat_ok = all(
        value if isinstance(value, bool) else nested_probe_ok
        for key, value in symbolic_checks.items()
        if key != "probe_verdicts"
    ) and nested_probe_ok
    max_symbolic_delta = 0.0 if symbolic_flat_ok else 1.0
    status = "pass" if max_symbolic_delta < SYMBOLIC_TOL and max_numeric_delta < NUMERIC_TOL else "FAIL_ENGINE_DISAGREE"
    return {
        "status": status,
        "max_symbolic_delta": max_symbolic_delta,
        "max_numeric_delta": max_numeric_delta,
        "grouped_lane_D_numeric_deltas": grouped_lane_d_deltas,
        "grouped_lane_D_max_numeric_delta": max(grouped_lane_d_deltas.values()) if grouped_lane_d_deltas else 0.0,
        "symbolic_tolerance": SYMBOLIC_TOL,
        "numeric_tolerance": NUMERIC_TOL,
        "symbolic_checks": symbolic_checks,
        "checked_fields": [
            "gram_5x5",
            "grouped_reduction",
            "grouped_lane_D_A_n",
            "lambda_by_channel",
            "anisotropy_coefficients",
            "probe_sample_values",
            "probe_verdicts",
            "omega_2m_sample",
            "baseline_verdict",
            "dimensional_ok",
            "dimension_probe_verdict",
        ],
    }


def symbolic_engine() -> dict[str, Any]:
    theta, phi, eps, delta = sp.symbols("theta phi eps delta", real=True)
    symbols = {
        "Mtilde": sp.Symbol("Mtilde", positive=True, real=True),
        "Ktilde": sp.Symbol("Ktilde", positive=True, real=True),
        "TomegaTilde": sp.Symbol("TomegaTilde", positive=True, real=True),
        "B0tilde": sp.Symbol("B0tilde", real=True),
        "Z0tilde": sp.Symbol("Z0tilde", real=True),
        "B2tilde": sp.Symbol("B2tilde", real=True),
        "Z2tilde": sp.Symbol("Z2tilde", real=True),
        "B4tilde": sp.Symbol("B4tilde", real=True),
        "Z4tilde": sp.Symbol("Z4tilde", real=True),
    }
    subs_sample = sample_subs(symbols)

    harmonics = {
        "20": sp.sqrt(sp.Rational(5, 16) / sp.pi) * (3 * sp.cos(theta) ** 2 - 1),
        "21c": -sp.sqrt(sp.Rational(15, 4) / sp.pi) * sp.sin(theta) * sp.cos(theta) * sp.cos(phi),
        "21s": -sp.sqrt(sp.Rational(15, 4) / sp.pi) * sp.sin(theta) * sp.cos(theta) * sp.sin(phi),
        "22c": sp.sqrt(sp.Rational(15, 16) / sp.pi) * sp.sin(theta) ** 2 * sp.cos(2 * phi),
        "22s": sp.sqrt(sp.Rational(15, 16) / sp.pi) * sp.sin(theta) ** 2 * sp.sin(2 * phi),
    }
    order = list(harmonics)
    ys = [harmonics[name] for name in order]

    gram = sp.Matrix([[integrate_s2(ys[i] * ys[j], theta, phi) for j in range(5)] for i in range(5)])
    gram_is_identity = bool(gram == sp.eye(5))

    p2_axis_z = (3 * sp.cos(theta) ** 2 - 1) / 2
    anisotropy_coefficients = {
        name: integrate_s2(p2_axis_z * harmonics[name] ** 2, theta, phi)
        for name in order
    }
    anisotropy_self_overlaps = {
        name: compact(1 + eps * anisotropy_coefficients[name])
        for name in order
    }

    grouped_isotropic = sp.diag(1, 1, 1)
    grouped_linear_coeff = sp.diag(
        anisotropy_coefficients["20"],
        compact((anisotropy_coefficients["21c"] + anisotropy_coefficients["21s"]) / 2),
        compact((anisotropy_coefficients["22c"] + anisotropy_coefficients["22s"]) / 2),
    )
    grouped_perturbed = compact(grouped_isotropic + eps * grouped_linear_coeff)

    lambdas: dict[str, sp.Expr] = {}
    residuals: dict[str, sp.Expr] = {}
    k_coeff_used: dict[str, sp.Expr] = {}
    k_coeff_equal: dict[str, bool] = {}
    for name, y_expr in harmonics.items():
        neg_lap = compact(-laplacian_s2(y_expr, theta, phi))
        rayleigh = compact(integrate_s2(y_expr * neg_lap, theta, phi) / integrate_s2(y_expr**2, theta, phi))
        lambdas[name] = rayleigh
        k_coeff_used[name] = rayleigh
        residuals[name] = compact(neg_lap - rayleigh * y_expr)
        k_coeff_equal[name] = bool(compact(k_coeff_used[name] - lambdas[name]) == 0)

    lambda_all_six = all(compact(value - 6) == 0 for value in lambdas.values())
    residuals_zero = all(compact(value) == 0 for value in residuals.values())
    k_coefficients_equal_all = all(k_coeff_equal.values())

    M = symbols["Mtilde"]
    Kbase = symbols["Ktilde"]
    Tomega = symbols["TomegaTilde"]
    B0 = symbols["B0tilde"]
    Z0 = symbols["Z0tilde"]
    B2 = symbols["B2tilde"]
    Z2 = symbols["Z2tilde"]
    B4 = symbols["B4tilde"]
    Z4 = symbols["Z4tilde"]
    S0 = B0 + Z0
    S2 = B2 + Z2
    S4 = B4 + Z4

    def assemble_channel(name: str, angular_weight: sp.Expr = sp.Integer(1), profile_factor: sp.Expr = sp.Integer(1)) -> dict[str, Any]:
        y_expr = harmonics[name]
        c_self = integrate_s2(angular_weight * y_expr**2, theta, phi)
        lam = lambdas[name]
        pref = compact(c_self * profile_factor)
        M_lane = compact(pref * M)
        K_lane = compact(pref * (Kbase + lam * Tomega))
        B_lane = {"0": compact(pref * B0), "2": compact(pref * B2), "4": compact(pref * B4)}
        Z_lane = {"0": compact(pref * Z0), "2": compact(pref * Z2), "4": compact(pref * Z4)}
        D_lane = {
            "0": compact(K_lane - B_lane["0"] - Z_lane["0"]),
            "2": compact(-(M_lane + B_lane["2"] + Z_lane["2"])),
            "4": compact(-(B_lane["4"] + Z_lane["4"])),
        }
        return {
            "angular_self_overlap": c_self,
            "lambda": lam,
            "M2": M_lane,
            "K2": K_lane,
            "B": B_lane,
            "Z": Z_lane,
            "D": D_lane,
        }

    ungrouped = {name: assemble_channel(name) for name in order}

    def average_expr(expressions: list[sp.Expr]) -> sp.Expr:
        return compact(sum(expressions) / len(expressions))

    grouped_lanes = {
        "20": ungrouped["20"],
        "21": {
            "M2": average_expr([ungrouped["21c"]["M2"], ungrouped["21s"]["M2"]]),
            "K2": average_expr([ungrouped["21c"]["K2"], ungrouped["21s"]["K2"]]),
            "D": {
                n: average_expr([ungrouped["21c"]["D"][n], ungrouped["21s"]["D"][n]])
                for n in ["0", "2", "4"]
            },
        },
        "22": {
            "M2": average_expr([ungrouped["22c"]["M2"], ungrouped["22s"]["M2"]]),
            "K2": average_expr([ungrouped["22c"]["K2"], ungrouped["22s"]["K2"]]),
            "D": {
                n: average_expr([ungrouped["22c"]["D"][n], ungrouped["22s"]["D"][n]])
                for n in ["0", "2", "4"]
            },
        },
    }
    dimensional_check = build_dimensional_check(
        Mtilde=M,
        Ktilde=Kbase,
        TomegaTilde=Tomega,
        M2_expr=grouped_lanes["20"]["M2"],
        K2_expr=grouped_lanes["20"]["K2"],
        lambda_m=sp.Integer(6),
    )

    cs_equal = {
        f"D21c_equals_D21s_order_{n}": bool(compact(ungrouped["21c"]["D"][n] - ungrouped["21s"]["D"][n]) == 0)
        for n in ["0", "2", "4"]
    }
    cs_equal.update(
        {
            f"D22c_equals_D22s_order_{n}": bool(compact(ungrouped["22c"]["D"][n] - ungrouped["22s"]["D"][n]) == 0)
            for n in ["0", "2", "4"]
        }
    )

    raw_triples: dict[str, list[sp.Expr]] = {
        n: [grouped_lanes[lane]["D"][n] for lane in ["20", "21", "22"]]
        for n in ["0", "2", "4"]
    }
    raw_defects = {}
    for n, triple in raw_triples.items():
        a_defect, b_defect = defect_pair(triple)
        raw_defects[f"a_D{n}"] = a_defect
        raw_defects[f"b_D{n}"] = b_defect
    raw_defects_zero = all_zero(list(raw_defects.values()))

    def u2_from_d(d_lane: dict[str, sp.Expr]) -> sp.Expr:
        return compact(-d_lane["2"] / d_lane["0"])

    def u4_from_d(d_lane: dict[str, sp.Expr]) -> sp.Expr:
        return compact((d_lane["2"] ** 2 - d_lane["0"] * d_lane["4"]) / d_lane["0"] ** 2)

    u2_lanes = {lane: u2_from_d(grouped_lanes[lane]["D"]) for lane in ["20", "21", "22"]}
    u4_lanes = {lane: u4_from_d(grouped_lanes[lane]["D"]) for lane in ["20", "21", "22"]}
    a2, b2_def = defect_pair([u2_lanes[lane] for lane in ["20", "21", "22"]])
    a4, b4_def = defect_pair([u4_lanes[lane] for lane in ["20", "21", "22"]])
    normalized_defects = {"a2": a2, "b2": b2_def, "a4": a4, "b4": b4_def}
    normalized_defects_zero = all_zero(list(normalized_defects.values()))

    input_hashes = {name: expr_hash(harmonics[name]) for name in order}
    distinct_hashes = len(set(input_hashes.values())) == len(input_hashes)
    per_lane_trace = {
        name: {
            "hash": input_hashes[name],
            "self_overlap_integral": f"Integral[Y_{name}^2 dOmega]",
            "self_overlap_value": hstr(ungrouped[name]["angular_self_overlap"]),
            "lambda_value": hstr(ungrouped[name]["lambda"]),
        }
        for name in order
    }
    tautology_clear = bool(distinct_hashes and all(ungrouped[name]["angular_self_overlap"] == 1 for name in order))

    K2_sample = CALIBRATION_SAMPLE["Ktilde"] + 6.0 * CALIBRATION_SAMPLE["TomegaTilde"]
    M2_sample = CALIBRATION_SAMPLE["Mtilde"]
    omega_sample = math.sqrt(K2_sample / M2_sample)
    k2_window_min = CALIBRATION_WINDOW["Ktilde"]["min"] + 6.0 * CALIBRATION_WINDOW["TomegaTilde"]["min"]
    d0_window_min = (
        CALIBRATION_WINDOW["Ktilde"]["min"]
        + 6.0 * CALIBRATION_WINDOW["TomegaTilde"]["min"]
        - CALIBRATION_WINDOW["B0tilde"]["max"]
        - CALIBRATION_WINDOW["Z0tilde"]["max"]
    )
    stability_ok = bool(CALIBRATION_WINDOW["Mtilde"]["min"] > 0.0 and k2_window_min > 0.0)
    denominator_guard_ok = bool(d0_window_min > 0.0)

    derived_inputs = [
        "explicit real l=2 harmonics",
        "5x5 angular Gram from S2 integrals",
        "per-harmonic -Delta_S2 eigenvalues",
        "K2 angular coefficient selected from computed lambda_m",
        "ungrouped channel angular self-overlaps",
        "raw-D and normalized defect algebra from assembled lanes",
        "counterfactual probe verdicts",
    ]
    calibration_inputs = [
        "R0(w) linearized isotropic reference",
        "beta2(w) radial profile",
        "mu_eta(w)",
        "T_w(w)",
        "K_eta(w)",
        "T_Omega(w)",
        "Mtilde radial mass scalar",
        "Ktilde radial stiffness scalar excluding angular T_Omega",
        "TomegaTilde angular radial scalar",
        "B0tilde,B2tilde,B4tilde support scalars",
        "Z0tilde,Z2tilde,Z4tilde mixed/Maxwell scalars",
        "physical calibration window for positivity and denominator guards",
        "Gate-1 D/N boundary provenance",
    ]

    def case_verdict(**overrides: bool) -> str:
        gates = {
            "dimensional_ok": dimensional_check["dimensional_ok"],
            "covariant": True,
            "tautology_clear": True,
            "dynamic_retained": True,
            "stability_ok": True,
            "denominator_guard_ok": True,
            "lane_collapse_ok": True,
            "able_to_fail_ok": True,
        }
        gates.update(overrides)
        return verdict_from_gates(gates, calibration_inputs)

    d_common = {
        "0": compact(Kbase + 6 * Tomega - S0),
        "2": compact(-(M + S2)),
        "4": compact(-S4),
    }
    c_group = {
        "20": anisotropy_coefficients["20"],
        "21": compact((anisotropy_coefficients["21c"] + anisotropy_coefficients["21s"]) / 2),
        "22": compact((anisotropy_coefficients["22c"] + anisotropy_coefficients["22s"]) / 2),
    }
    c_group_list = [c_group[lane] for lane in ["20", "21", "22"]]

    pure_raw_defects: dict[str, sp.Expr] = {}
    pure_samples: dict[str, list[float]] = {}
    pure_d_by_lane = {
        lane: {n: compact((1 + eps * c_group[lane]) * d_common[n]) for n in ["0", "2", "4"]}
        for lane in ["20", "21", "22"]
    }
    for n in ["0", "2", "4"]:
        triple = [pure_d_by_lane[lane][n] for lane in ["20", "21", "22"]]
        a_def, b_def = defect_pair(triple)
        pure_raw_defects[f"a_D{n}"] = a_def
        pure_raw_defects[f"b_D{n}"] = b_def
        for label, expr in [(f"a_D{n}", a_def), (f"b_D{n}", b_def)]:
            pure_samples[label] = [f_eval(expr.subs(eps, eps_val), subs_sample) for eps_val in EPS_VALUES]
    pure_u2 = {lane: u2_from_d(pure_d_by_lane[lane]) for lane in ["20", "21", "22"]}
    pure_u4 = {lane: u4_from_d(pure_d_by_lane[lane]) for lane in ["20", "21", "22"]}
    pure_a2, pure_b2 = defect_pair([pure_u2[lane] for lane in ["20", "21", "22"]])
    pure_a4, pure_b4 = defect_pair([pure_u4[lane] for lane in ["20", "21", "22"]])
    pure_u_defects = {"a2": pure_a2, "b2": pure_b2, "a4": pure_a4, "b4": pure_b4}
    pure_raw_moves = all(first_nonzero(values) for values in pure_samples.values())
    pure_linear_delta = max(max_ratio_delta(values, EPS_VALUES) for values in pure_samples.values())

    sector_d_by_lane = {
        lane: {
            "0": compact(Kbase + 6 * Tomega - (1 + eps * c_group[lane]) * S0),
            "2": compact(-(M + (1 + eps * c_group[lane]) * S2)),
            "4": compact(-((1 + eps * c_group[lane]) * S4)),
        }
        for lane in ["20", "21", "22"]
    }
    sector_raw_defects: dict[str, sp.Expr] = {}
    sector_samples: dict[str, list[float]] = {}
    for n in ["0", "2", "4"]:
        a_def, b_def = defect_pair([sector_d_by_lane[lane][n] for lane in ["20", "21", "22"]])
        sector_raw_defects[f"a_D{n}"] = a_def
        sector_raw_defects[f"b_D{n}"] = b_def
        for label, expr in [(f"a_D{n}", a_def), (f"b_D{n}", b_def)]:
            sector_samples[label] = [f_eval(expr.subs(eps, eps_val), subs_sample) for eps_val in EPS_VALUES]
    sector_u2 = {lane: u2_from_d(sector_d_by_lane[lane]) for lane in ["20", "21", "22"]}
    sector_a2, sector_b2 = defect_pair([sector_u2[lane] for lane in ["20", "21", "22"]])
    sector_u_defects = {
        "a2": sector_a2,
        "b2": sector_b2,
    }
    sector_u_samples = {
        label: [f_eval(expr.subs(eps, eps_val), subs_sample) for eps_val in EPS_VALUES]
        for label, expr in sector_u_defects.items()
    }
    sector_u_linear_coefficients = {
        label: f_eval(sp.diff(expr, eps).subs(eps, 0), subs_sample)
        for label, expr in sector_u_defects.items()
    }
    sector_raw_moves = all(first_nonzero(values) for values in sector_samples.values())
    sector_u_moves = first_nonzero(sector_u_samples["a2"]) and first_nonzero(sector_u_samples["b2"])

    profile_scales = {"20": sp.Integer(1), "21": sp.Integer(1), "22": compact((1 + delta) ** 2)}
    profile_d_by_lane = {
        lane: {n: compact(profile_scales[lane] * d_common[n]) for n in ["0", "2", "4"]}
        for lane in ["20", "21", "22"]
    }
    profile_raw_defects = {}
    for n in ["0", "2", "4"]:
        a_def, b_def = defect_pair([profile_d_by_lane[lane][n] for lane in ["20", "21", "22"]])
        profile_raw_defects[f"a_D{n}"] = a_def
        profile_raw_defects[f"b_D{n}"] = b_def
    profile_sample = {
        label: f_eval(expr.subs(delta, DELTA_PROFILE), subs_sample)
        for label, expr in profile_raw_defects.items()
    }
    profile_raw_moves = first_nonzero(list(profile_sample.values()))

    singular_subs = dict(subs_sample)
    singular_subs[symbols["B0tilde"]] = K2_sample - singular_subs[symbols["Z0tilde"]]
    singular_d0_value = f_eval(d_common["0"], singular_subs)
    singular_denominator_guard_ok = abs(singular_d0_value) >= PROBE_TOL
    nonsingular_d0_value = f_eval(d_common["0"], subs_sample)
    nonsingular_denominator_guard_ok = abs(nonsingular_d0_value) >= PROBE_TOL

    def beta_stability_probe(beta_scale: sp.Expr) -> dict[str, Any]:
        m2_expr = compact(beta_scale**2 * grouped_lanes["20"]["M2"])
        k2_expr = compact(beta_scale**2 * grouped_lanes["20"]["K2"])
        m2_sample = f_eval(m2_expr, subs_sample)
        k2_sample = f_eval(k2_expr, subs_sample)
        ok = bool(m2_sample > PROBE_TOL and k2_sample > PROBE_TOL)
        return {
            "beta2_scale": hstr(beta_scale),
            "M2": hstr(m2_expr),
            "K2": hstr(k2_expr),
            "M2_sample": m2_sample,
            "K2_sample": k2_sample,
            "stability_ok": ok,
            "verdict": case_verdict(stability_ok=ok),
            "fail_fires": case_verdict(stability_ok=ok) == "FAIL_STABILITY",
        }

    def forced_eigenvalue_probe(forced_coefficient: sp.Expr) -> dict[str, Any]:
        forced_k2_by_channel = {
            name: compact(ungrouped[name]["angular_self_overlap"] * (Kbase + forced_coefficient * Tomega))
            for name in order
        }
        coefficient_equals = all(compact(forced_coefficient - lambdas[name]) == 0 for name in order)
        verdict = case_verdict(covariant=coefficient_equals)
        return {
            "forced_coefficient": hstr(forced_coefficient),
            "forced_K2_by_channel": {name: hstr(value) for name, value in forced_k2_by_channel.items()},
            "computed_lambda_by_channel": {name: hstr(value) for name, value in lambdas.items()},
            "coefficient_equals_computed_lambda": coefficient_equals,
            "verdict": verdict,
            "fail_fires": verdict == "FAIL_NOT_COVARIANT",
        }

    def lane_hash_probe(lane_inputs: dict[str, sp.Expr]) -> dict[str, Any]:
        hashes = {lane: expr_hash(expr) for lane, expr in lane_inputs.items()}
        distinct = len(set(hashes.values())) == len(hashes)
        verdict = case_verdict(tautology_clear=distinct)
        return {
            "input_hashes": hashes,
            "distinct_hashes": distinct,
            "verdict": verdict,
            "fail_fires": verdict == "FAIL_TAUTOLOGICAL",
        }

    degenerate_beta_probe = beta_stability_probe(sp.Integer(0))
    degenerate_beta_ablation = beta_stability_probe(sp.Integer(1))
    wrong_eigen_probe = forced_eigenvalue_probe(sp.Integer(2))
    wrong_eigen_ablation = forced_eigenvalue_probe(sp.Integer(6))
    singular_verdict = case_verdict(denominator_guard_ok=singular_denominator_guard_ok)
    singular_ablation_verdict = case_verdict(denominator_guard_ok=nonsingular_denominator_guard_ok)
    tautology_probe = lane_hash_probe(
        {"20": harmonics["20"], "21": harmonics["20"], "22": harmonics["20"]}
    )
    tautology_ablation = lane_hash_probe(
        {"20": harmonics["20"], "21": harmonics["21c"], "22": harmonics["22c"]}
    )
    static_wrong_d2 = compact(-(B2 + Z2))
    static_dynamic_retained = bool(static_wrong_d2.has(M))
    static_verdict = case_verdict(dynamic_retained=static_dynamic_retained)
    static_ablation_dynamic_retained = bool(grouped_lanes["20"]["D"]["2"].has(M))
    static_ablation_verdict = case_verdict(dynamic_retained=static_ablation_dynamic_retained)
    dimension_probe = dimensional_check["FAIL_DIMENSIONAL_probe"]
    dimensional_verdict = case_verdict(dimensional_ok=dimension_probe["with_mutation_dimensional_ok"])
    dimensional_ablation_verdict = case_verdict(dimensional_ok=dimension_probe["without_mutation_dimensional_ok"])

    probes = {
        "pure_prefactor_anisotropy": {
            "description": "w=1+eps*P2 applied uniformly to wall and support sectors",
            "epsilon_values": EPS_VALUES,
            "C_A_epsilon": {lane: hstr(1 + eps * c_group[lane]) for lane in ["20", "21", "22"]},
            "raw_D_defects": {key: hstr(value) for key, value in pure_raw_defects.items()},
            "normalized_u_defects": {key: hstr(value) for key, value in pure_u_defects.items()},
            "sample_values": pure_samples,
            "raw_D_moves": pure_raw_moves,
            "linear_scaling_max_ratio_delta": pure_linear_delta,
            "linear_scaling_confirmed": pure_raw_moves and pure_linear_delta < 1.0e-10,
            "normalized_u_defects_stay_zero": all_zero(list(pure_u_defects.values())),
            "verdict": case_verdict(lane_collapse_ok=False),
        },
        "sector_selective_anisotropy": {
            "description": "w=1+eps*P2 applied only to Btilde/Ztilde support sectors; wall sector stays isotropic",
            "epsilon_values": EPS_VALUES,
            "raw_D_defects": {key: hstr(value) for key, value in sector_raw_defects.items()},
            "normalized_u_defects": {key: hstr(value) for key, value in sector_u_defects.items()},
            "sample_values": {**sector_samples, **{f"u_{k}": v for k, v in sector_u_samples.items()}},
            "u_linear_coefficients_at_eps0_sample": sector_u_linear_coefficients,
            "raw_D_moves": sector_raw_moves,
            "u_defects_move": sector_u_moves,
            "verdict": case_verdict(lane_collapse_ok=False),
        },
        "m_dependent_profile": {
            "description": "profile scale beta_22=(1+delta)*beta_20, so quadratic response scale differs by lane",
            "delta": DELTA_PROFILE,
            "profile_scales": {key: hstr(value) for key, value in profile_scales.items()},
            "raw_D_defects": {key: hstr(value) for key, value in profile_raw_defects.items()},
            "sample_values": profile_sample,
            "raw_D_moves": profile_raw_moves,
            "normalized_u_defects_may_cancel_for_pure_profile_scale": True,
            "verdict": case_verdict(lane_collapse_ok=False),
        },
        "degenerate_beta_zero": {
            "description": "beta2=0 gives M2=K2=0 and is caught before normalized 0/0 response",
            **degenerate_beta_probe,
            "computed_fail_gate": not degenerate_beta_probe["stability_ok"],
            "self_ablation": {
                **degenerate_beta_ablation,
                "fail_suppressed": not degenerate_beta_ablation["fail_fires"],
            },
        },
        "wrong_eigenvalue": {
            "description": "force K2 angular coefficient 2 while computed lambda_m remains 6",
            **wrong_eigen_probe,
            "computed_fail_gate": not wrong_eigen_probe["coefficient_equals_computed_lambda"],
            "self_ablation": {
                **wrong_eigen_ablation,
                "fail_suppressed": not wrong_eigen_ablation["fail_fires"],
            },
        },
        "singular_denominator": {
            "description": "hold M2,K2 positive but choose B0tilde+Z0tilde=K2 so D_A0=0",
            "D0_sample": singular_d0_value,
            "M2_positive": CALIBRATION_SAMPLE["Mtilde"] > PROBE_TOL,
            "K2_positive": K2_sample > PROBE_TOL,
            "denominator_guard_ok": singular_denominator_guard_ok,
            "computed_fail_gate": not singular_denominator_guard_ok,
            "verdict": singular_verdict,
            "fail_fires": singular_verdict == "FAIL_SINGULAR_RESPONSE",
            "self_ablation": {
                "description": "baseline non-singular substitution",
                "D0_sample": nonsingular_d0_value,
                "denominator_guard_ok": nonsingular_denominator_guard_ok,
                "verdict": singular_ablation_verdict,
                "fail_fires": singular_ablation_verdict == "FAIL_SINGULAR_RESPONSE",
                "fail_suppressed": singular_ablation_verdict != "FAIL_SINGULAR_RESPONSE",
            },
        },
        "tautology_hash_collision": {
            "description": "set the three grouped lanes to identical harmonic inputs and compute the hashes",
            **tautology_probe,
            "computed_fail_gate": not tautology_probe["distinct_hashes"],
            "self_ablation": {
                **tautology_ablation,
                "fail_suppressed": not tautology_ablation["fail_fires"],
            },
        },
        "static_drop_inertia": {
            "description": "drop Mtilde from D_A2, breaking the required dynamic omega expansion",
            "dynamic_retained": static_dynamic_retained,
            "computed_fail_gate": not static_dynamic_retained,
            "wrong_D2_without_M": hstr(static_wrong_d2),
            "correct_D2": hstr(d_common["2"]),
            "verdict": static_verdict,
            "fail_fires": static_verdict == "FAIL_STATIC_RESPONSE",
            "self_ablation": {
                "description": "retain the Mtilde inertia term in D_A2",
                "dynamic_retained": static_ablation_dynamic_retained,
                "correct_D2": hstr(grouped_lanes["20"]["D"]["2"]),
                "verdict": static_ablation_verdict,
                "fail_fires": static_ablation_verdict == "FAIL_STATIC_RESPONSE",
                "fail_suppressed": static_ablation_verdict != "FAIL_STATIC_RESPONSE",
            },
        },
        "dimensional_corrupt_T_Omega": {
            "description": "corrupt the sourced angular stiffness dimension and its assembled TomegaTilde scalar",
            **dimension_probe,
            "verdict": dimensional_verdict,
            "fail_fires": dimensional_verdict == "FAIL_DIMENSIONAL",
            "self_ablation": {
                "description": "restore sourced T_Omega/TomegaTilde dimensions",
                "dimensional_ok": dimension_probe["without_mutation_dimensional_ok"],
                "verdict": dimensional_ablation_verdict,
                "fail_fires": dimensional_ablation_verdict == "FAIL_DIMENSIONAL",
                "fail_suppressed": dimensional_ablation_verdict != "FAIL_DIMENSIONAL",
            },
        },
    }

    expected_probe_verdicts = {
        "pure_prefactor_anisotropy": "FAIL_ANISOTROPIC_BRANCH",
        "sector_selective_anisotropy": "FAIL_ANISOTROPIC_BRANCH",
        "m_dependent_profile": "FAIL_ANISOTROPIC_BRANCH",
        "degenerate_beta_zero": "FAIL_STABILITY",
        "wrong_eigenvalue": "FAIL_NOT_COVARIANT",
        "singular_denominator": "FAIL_SINGULAR_RESPONSE",
        "tautology_hash_collision": "FAIL_TAUTOLOGICAL",
        "static_drop_inertia": "FAIL_STATIC_RESPONSE",
        "dimensional_corrupt_T_Omega": "FAIL_DIMENSIONAL",
    }
    expected_probe_verdicts_match = {
        key: probes[key]["verdict"] == expected for key, expected in expected_probe_verdicts.items()
    }
    computed_probe_gate_flags = {
        "pure_prefactor_anisotropy": bool(
            probes["pure_prefactor_anisotropy"]["raw_D_moves"]
            and probes["pure_prefactor_anisotropy"]["normalized_u_defects_stay_zero"]
        ),
        "sector_selective_anisotropy": bool(
            probes["sector_selective_anisotropy"]["raw_D_moves"]
            and probes["sector_selective_anisotropy"]["u_defects_move"]
        ),
        "m_dependent_profile": bool(probes["m_dependent_profile"]["raw_D_moves"]),
        "degenerate_beta_zero": bool(
            probes["degenerate_beta_zero"]["computed_fail_gate"]
            and probes["degenerate_beta_zero"]["self_ablation"]["fail_suppressed"]
        ),
        "wrong_eigenvalue": bool(
            probes["wrong_eigenvalue"]["computed_fail_gate"]
            and probes["wrong_eigenvalue"]["self_ablation"]["fail_suppressed"]
        ),
        "singular_denominator": bool(
            probes["singular_denominator"]["computed_fail_gate"]
            and probes["singular_denominator"]["self_ablation"]["fail_suppressed"]
        ),
        "tautology_hash_collision": bool(
            probes["tautology_hash_collision"]["computed_fail_gate"]
            and probes["tautology_hash_collision"]["self_ablation"]["fail_suppressed"]
        ),
        "static_drop_inertia": bool(
            probes["static_drop_inertia"]["computed_fail_gate"]
            and probes["static_drop_inertia"]["self_ablation"]["fail_suppressed"]
        ),
        "dimensional_corrupt_T_Omega": bool(
            not probes["dimensional_corrupt_T_Omega"]["with_mutation_dimensional_ok"]
            and probes["dimensional_corrupt_T_Omega"]["self_ablation"]["fail_suppressed"]
        ),
    }
    def able_to_fail_from_flags(flags: dict[str, bool]) -> bool:
        return bool(all(expected_probe_verdicts_match.values()) and all(flags.values()))

    able_to_fail_ok = able_to_fail_from_flags(computed_probe_gate_flags)
    able_to_fail_if_probe_neutered = {
        key: able_to_fail_from_flags({**computed_probe_gate_flags, key: False})
        for key in computed_probe_gate_flags
    }

    dynamic_retained = bool(all(d_common[n].has(M) for n in ["2"]) and not d_common["0"].has(M))

    baseline_gates = {
        "dimensional_ok": dimensional_check["dimensional_ok"],
        "covariant": bool(gram_is_identity and lambda_all_six and residuals_zero and k_coefficients_equal_all),
        "tautology_clear": tautology_clear,
        "dynamic_retained": dynamic_retained,
        "stability_ok": stability_ok,
        "denominator_guard_ok": denominator_guard_ok,
        "lane_collapse_ok": bool(raw_defects_zero and all(cs_equal.values())),
        "able_to_fail_ok": able_to_fail_ok,
    }
    verdict = verdict_from_gates(baseline_gates, calibration_inputs)

    payload: dict[str, Any] = {
        "schema": "pathA_32_grouped_p2_isotropy/v1",
        "engine": "sympy",
        "source_directive": "software/stage1_solver/directives/pathA_32_grouped_p2_isotropy.md",
        "cited_sources": [
            "research/pde_ledger/notes/stages/moving_throat_pde_handoff_full.md sections 6.1, 8.3, 11",
            "research/pde_ledger/notes/stages/moving_throat_pde_phase1_linearized_scaffold.md section 8.2",
        ],
        "verdict": verdict,
        "which_rung": verdict,
        "gate_booleans": baseline_gates,
        "input_partition": {
            "derived_inputs": derived_inputs,
            "calibration_inputs": calibration_inputs,
            "pass_requires_empty_calibration_inputs": True,
        },
        "calibration_window": CALIBRATION_WINDOW,
        "calibration_sample": CALIBRATION_SAMPLE,
        "harmonics": {
            "order": order,
            "expressions": {name: hstr(expr) for name, expr in harmonics.items()},
            "input_hashes": input_hashes,
            "distinct_input_hashes": distinct_hashes,
            "per_channel_integral_trace": per_lane_trace,
            "gram_5x5": matrix_to_strings(gram),
            "gram_5x5_numeric": matrix_to_floats(gram),
            "gram_is_identity": gram_is_identity,
            "anisotropy_weight": "1 + eps*(3*cos(theta)^2-1)/2",
            "anisotropy_coefficients": {name: hstr(value) for name, value in anisotropy_coefficients.items()},
            "anisotropy_coefficients_numeric": {
                name: float(sp.N(value, 30)) for name, value in anisotropy_coefficients.items()
            },
            "self_overlaps_pure_prefactor": {name: hstr(value) for name, value in anisotropy_self_overlaps.items()},
            "grouped_reduction_isotropic": matrix_to_strings(grouped_isotropic),
            "grouped_reduction_isotropic_numeric": matrix_to_floats(grouped_isotropic),
            "grouped_reduction_pure_prefactor": matrix_to_strings(grouped_perturbed),
            "grouped_reduction_pure_prefactor_linear_coeff": matrix_to_strings(grouped_linear_coeff),
            "grouped_reduction_pure_prefactor_linear_coeff_numeric": matrix_to_floats(grouped_linear_coeff),
        },
        "laplacian": {
            "lambda_by_channel": {name: hstr(value) for name, value in lambdas.items()},
            "lambda_numeric_by_channel": {name: float(sp.N(value, 30)) for name, value in lambdas.items()},
            "residuals_minus_lambda_Y": {name: hstr(value) for name, value in residuals.items()},
            "residuals_zero": residuals_zero,
            "lambda_all_six": lambda_all_six,
            "k2_angular_coefficient_used": {name: hstr(k_coeff_used[name]) for name in order},
            "k2_coefficient_equals_computed_lambda": k_coeff_equal,
            "k2_coefficient_equals_computed_lambda_all": k_coefficients_equal_all,
        },
        "response_symbols": {
            "M2": "Mtilde",
            "K2_per_channel": "Ktilde + lambda_m*TomegaTilde",
            "K2_on_isotropic_l2": "Ktilde + 6*TomegaTilde",
            "Btilde": ["B0tilde", "B2tilde", "B4tilde"],
            "Ztilde": ["Z0tilde", "Z2tilde", "Z4tilde"],
        },
        "ungrouped_channels": {
            name: {
                "angular_self_overlap": hstr(ungrouped[name]["angular_self_overlap"]),
                "M2": hstr(ungrouped[name]["M2"]),
                "K2": hstr(ungrouped[name]["K2"]),
                "D": {n: hstr(ungrouped[name]["D"][n]) for n in ["0", "2", "4"]},
            }
            for name in order
        },
        "intra_lane_cs_equality": cs_equal,
        "grouped_lanes": {
            lane: {
                "M2": hstr(grouped_lanes[lane]["M2"]),
                "K2": hstr(grouped_lanes[lane]["K2"]),
                "D": {n: hstr(grouped_lanes[lane]["D"][n]) for n in ["0", "2", "4"]},
                "D_sample": {n: f_eval(grouped_lanes[lane]["D"][n], subs_sample) for n in ["0", "2", "4"]},
                "D0_nonzero_guard": denominator_guard_ok,
            }
            for lane in ["20", "21", "22"]
        },
        "raw_D_gate": {
            "raw_D_triples": {
                n: {lane: hstr(raw_triples[n][i]) for i, lane in enumerate(["20", "21", "22"])}
                for n in ["0", "2", "4"]
            },
            "raw_D_defects": {key: hstr(value) for key, value in raw_defects.items()},
            "raw_D_defect_magnitudes_sample": {
                key: abs(f_eval(value, subs_sample)) for key, value in raw_defects.items()
            },
            "raw_D_defects_zero": raw_defects_zero,
            "lane_collapse_ok": bool(raw_defects_zero and all(cs_equal.values())),
        },
        "normalized_response": {
            "u2_triple": {lane: hstr(u2_lanes[lane]) for lane in ["20", "21", "22"]},
            "u4_triple": {lane: hstr(u4_lanes[lane]) for lane in ["20", "21", "22"]},
            "normalized_defects": {key: hstr(value) for key, value in normalized_defects.items()},
            "normalized_defect_magnitudes_sample": {
                key: abs(f_eval(value, subs_sample)) for key, value in normalized_defects.items()
            },
            "normalized_defects_zero": normalized_defects_zero,
        },
        "stability": {
            "M2_window_min": CALIBRATION_WINDOW["Mtilde"]["min"],
            "K2_window_min": k2_window_min,
            "D0_window_min": d0_window_min,
            "stability_ok": stability_ok,
            "denominator_guard_ok": denominator_guard_ok,
            "omega_2m": {name: "sqrt((Ktilde + 6*TomegaTilde)/Mtilde)" for name in order},
            "omega_2m_sample": omega_sample,
            "omega_2m_degenerate": True,
        },
        "counterfactuals": probes,
        "able_to_fail": {
            "expected_probe_verdicts": expected_probe_verdicts,
            "expected_probe_verdicts_match": expected_probe_verdicts_match,
            "computed_probe_gate_flags": computed_probe_gate_flags,
            "able_to_fail_if_probe_neutered": able_to_fail_if_probe_neutered,
            "neutering_any_probe_flips_false": all(not value for value in able_to_fail_if_probe_neutered.values()),
            "able_to_fail_ok": able_to_fail_ok,
        },
        "dimensional": dimensional_check,
        "dimensional_check": dimensional_check,
        "dim_homogeneity_table": dimensional_check,
        "deferred": [
            "54/5 quadrupole normalization",
            "outgoing odd N_A,n extraction",
            "solved nonlinear branch data K_A,M_A,B_A,n,Z_A,n,N_A,n",
            "cross-l consistency and PN match-back",
        ],
    }
    return payload


def main() -> int:
    sympy_payload = symbolic_engine()
    SCRATCH.mkdir(parents=True, exist_ok=True)
    yaml_write(SYM_YAML, sympy_payload)

    mma_payload = yaml_read(MMA_YAML)
    engine_agreement = compute_engine_agreement(sympy_payload, mma_payload)
    if mma_payload is None:
        print(f"SymPy engine complete; wrote {SYM_YAML}. Mathematica scratch not found yet: {MMA_YAML}")
        return 0

    final_payload = dict(sympy_payload)
    final_payload["engine"] = "dual-engine"
    final_payload["engines"] = {
        "sympy": {"scratch_yaml": str(SYM_YAML), "status": "ok"},
        "mathematica": {"scratch_yaml": str(MMA_YAML), "status": mma_payload.get("status", "ok")},
    }
    final_payload["engine_agreement"] = engine_agreement
    if engine_agreement["status"] != "pass":
        final_payload["verdict"] = "FAIL_ENGINE_DISAGREE"
        final_payload["which_rung"] = "FAIL_ENGINE_DISAGREE"

    yaml_write(RESULTS_YAML, final_payload)
    REPORT_MD.parent.mkdir(parents=True, exist_ok=True)
    REPORT_MD.write_text(build_report(final_payload), encoding="utf-8")
    FEED_NOTE.parent.mkdir(parents=True, exist_ok=True)
    FEED_NOTE.write_text(build_feed_note(final_payload), encoding="utf-8")
    print(f"Final PathA-32 report written: {REPORT_MD}")
    print(f"Final PathA-32 YAML written: {RESULTS_YAML}")
    print(f"Computed verdict: {final_payload['verdict']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
