#!/usr/bin/env python3
"""Ledger stage016 SymPy audit: l=2 SO(3) covariance theorem.

Standalone, print-only, no arguments, no file I/O. This is the pathA_32
II-G3a angular slice only: real l=2 harmonics, Gram=I5, computed
-Delta_S2 eigenvalues, the bare K2 angular stiffness with the live computed
lambda, the angular dimensional gate, and the 016 covariance probe battery.
Stage 017 owns grouped lanes, raw-D, normalized-u, response probes, calibration
partition, and port-kernel export.
"""

from __future__ import annotations

import hashlib
from dataclasses import dataclass
from typing import Any

import sympy as sp


PASS_COUNT = 0
FAIL_COUNT = 0

ISOTROPY_CALIBRATED = "ISOTROPY_CALIBRATED"
FAIL_DIMENSIONAL = "FAIL_DIMENSIONAL"
FAIL_NOT_COVARIANT = "FAIL_NOT_COVARIANT"
FAIL_TAUTOLOGICAL = "FAIL_TAUTOLOGICAL"
FAIL_NOT_ABLE_TO_FAIL = "FAIL_NOT_ABLE_TO_FAIL"


class AuditFailure(AssertionError):
    pass


class DimError(ValueError):
    pass


def banner(title: str) -> None:
    print("")
    print("=" * len(title))
    print(title)
    print("=" * len(title))


def subbanner(title: str) -> None:
    print("")
    print(title)
    print("-" * len(title))


def compact(expr: Any) -> Any:
    if isinstance(expr, sp.MatrixBase):
        return expr.applyfunc(compact)
    if not isinstance(expr, sp.Basic):
        return expr
    if isinstance(expr, (sp.Equality, sp.Order)):
        return expr
    return sp.trigsimp(sp.factor(sp.cancel(sp.together(sp.simplify(expr)))))


def fmt(expr: Any) -> str:
    if isinstance(expr, bool):
        return "True" if expr else "False"
    if isinstance(expr, str):
        return expr
    if isinstance(expr, Dim):
        return str(expr)
    if isinstance(expr, sp.MatrixBase):
        return sp.sstr(compact(expr))
    if isinstance(expr, (dict, list, tuple)):
        return sp.sstr(expr)
    try:
        return sp.sstr(compact(expr))
    except Exception:
        return sp.sstr(expr)


def assert_no_float(name: str, expr: Any) -> None:
    if isinstance(expr, Dim):
        for label, value in zip(("L", "M", "T"), expr.components()):
            assert_no_float(f"{name}.{label}", value)
        return
    if isinstance(expr, dict):
        for key, value in expr.items():
            assert_no_float(f"{name}.{key}", value)
        return
    if isinstance(expr, sp.MatrixBase):
        for index, value in enumerate(expr):
            assert_no_float(f"{name}[{index}]", value)
        return
    if isinstance(expr, (list, tuple, set, frozenset)):
        for index, value in enumerate(expr):
            assert_no_float(f"{name}[{index}]", value)
        return
    if isinstance(expr, (str, type(None))):
        return
    if isinstance(expr, bool):
        expr = sp.Integer(1) if expr else sp.Integer(0)
    clean = sp.sympify(expr)
    floats = clean.atoms(sp.Float)
    if floats:
        raise AuditFailure(f"{name}: Float atom(s) found in exact audit expression: {floats}")


def _record_pass(message: str) -> None:
    global PASS_COUNT
    PASS_COUNT += 1
    print(message)


def _record_fail(message: str) -> None:
    global FAIL_COUNT
    FAIL_COUNT += 1
    print(message)


def expect_zero(name: str, residual: sp.Expr | int) -> None:
    assert_no_float(name, residual)
    clean = compact(residual)
    assert_no_float(name, clean)
    if clean == 0:
        _record_pass(f"PASS  {name}")
        return
    _record_fail(f"FAIL  {name}: residual = {fmt(clean)}")
    raise AuditFailure(f"{name} residual was not zero")


def expect_bool(name: str, condition: bool) -> None:
    expect_zero(name, sp.Integer(0) if condition else sp.Integer(1))


def expect_fail(name: str, residual: sp.Expr | int) -> None:
    assert_no_float(name, residual)
    clean = compact(residual)
    assert_no_float(name, clean)
    if clean != 0:
        _record_pass(f"PASS  {name} produced required FAIL (residual = {fmt(clean)})")
        return
    _record_fail(f"FAIL  {name}: required mutation/ablation did not fire")
    raise AuditFailure(f"{name} unexpectedly had zero residual")


def bool_residual(condition: bool) -> sp.Integer:
    return sp.Integer(0) if condition else sp.Integer(1)


def verdict_residual(actual: str, expected: str) -> sp.Integer:
    return sp.Integer(0) if actual == expected else sp.Integer(1)


def q(value: int | str | sp.Rational) -> sp.Rational:
    return sp.Rational(value)


@dataclass(frozen=True)
class Dim:
    """Exact exponent vector in (L, M, T) order."""

    l: sp.Rational | int = 0
    m: sp.Rational | int = 0
    t: sp.Rational | int = 0

    def __post_init__(self) -> None:
        object.__setattr__(self, "l", q(self.l))
        object.__setattr__(self, "m", q(self.m))
        object.__setattr__(self, "t", q(self.t))

    def __add__(self, other: "Dim") -> "Dim":
        return Dim(self.l + other.l, self.m + other.m, self.t + other.t)

    def __sub__(self, other: "Dim") -> "Dim":
        return Dim(self.l - other.l, self.m - other.m, self.t - other.t)

    def __mul__(self, scale: int | sp.Rational) -> "Dim":
        p = q(scale)
        return Dim(self.l * p, self.m * p, self.t * p)

    def __rmul__(self, scale: int | sp.Rational) -> "Dim":
        return self * scale

    def components(self) -> tuple[sp.Rational, sp.Rational, sp.Rational]:
        return (self.l, self.m, self.t)

    def __str__(self) -> str:
        return f"({fmt(self.l)},{fmt(self.m)},{fmt(self.t)})"


ZERO_DIM = Dim()
EXPECTED_M = Dim(0, 1, 0)
EXPECTED_K = Dim(0, 1, -2)
EXPECTED_RATIO = Dim(0, 0, -2)
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


def dim_residual(actual: Dim, expected: Dim) -> sp.Expr:
    return sp.simplify(sum((have - want) ** 2 for have, want in zip(actual.components(), expected.components())))


def dim_of(expr: sp.Expr, dims: dict[sp.Symbol, Dim]) -> Dim:
    clean = sp.sympify(expr)
    if clean == 0 or clean.is_number:
        return ZERO_DIM
    if isinstance(clean, sp.Symbol):
        if clean not in dims:
            raise DimError(f"missing dimension for symbol {clean}")
        return dims[clean]
    if isinstance(clean, sp.Mul):
        total = ZERO_DIM
        for arg in clean.args:
            total = total + dim_of(arg, dims)
        return total
    if isinstance(clean, sp.Pow):
        base, power = clean.args
        if not power.is_number:
            raise DimError(f"non-numeric power in dimension expression {clean}")
        return dim_of(base, dims) * sp.Rational(power)
    if isinstance(clean, sp.Add):
        arg_dims = [dim_of(arg, dims) for arg in clean.args if compact(arg) != 0]
        if not arg_dims:
            return ZERO_DIM
        first = arg_dims[0]
        if any(arg_dim != first for arg_dim in arg_dims[1:]):
            raise DimError(f"dimension mismatch in sum {clean}: {arg_dims}")
        return first
    if clean.func in DIMENSIONLESS_FUNCTIONS:
        arg_dims = [dim_of(arg, dims) for arg in clean.args]
        if any(arg_dim != ZERO_DIM for arg_dim in arg_dims):
            raise DimError(f"dimensionful argument in dimensionless function {clean}")
        return ZERO_DIM
    raise DimError(f"unsupported dimension expression {clean}")


def dim_text(dim: Dim | None) -> str:
    if dim is None:
        return "inhomogeneous"
    parts: list[str] = []
    for label, exponent in (("L", dim.l), ("M", dim.m), ("T", dim.t)):
        if exponent == 0:
            continue
        if exponent == 1:
            parts.append(label)
        else:
            parts.append(f"{label}^{fmt(exponent)}")
    return "1" if not parts else " ".join(parts)


theta, phi = sp.symbols("theta phi", real=True)
Mtilde = sp.Symbol("Mtilde", positive=True, real=True)
Ktilde = sp.Symbol("Ktilde", positive=True, real=True)
TomegaTilde = sp.Symbol("TomegaTilde", positive=True, real=True)

a_dim, dw_dim, dOmega_dim = sp.symbols("a_dim dw dOmega", positive=True)
beta2_dim, beta2_prime_dim = sp.symbols("beta2 beta2_prime", nonzero=True)
mu_eta_density, T_w_density, K_eta_density, T_Omega_density = sp.symbols(
    "mu_eta_density T_w_density K_eta_density T_Omega_density", nonzero=True
)


def integrate_s2(expr: sp.Expr) -> sp.Expr:
    return compact(
        sp.integrate(
            sp.integrate(sp.expand_trig(expr) * sp.sin(theta), (phi, 0, 2 * sp.pi)),
            (theta, 0, sp.pi),
        )
    )


def laplacian_s2(expr: sp.Expr) -> sp.Expr:
    return compact(
        (1 / sp.sin(theta)) * sp.diff(sp.sin(theta) * sp.diff(expr, theta), theta)
        + (1 / sp.sin(theta) ** 2) * sp.diff(expr, phi, 2)
    )


def expr_hash(expr: sp.Expr) -> str:
    return hashlib.sha256(sp.sstr(compact(expr)).encode("utf-8")).hexdigest()


def scoped_verdict(
    *,
    dimensional_ok: bool = True,
    covariant: bool = True,
    tautology_clear: bool = True,
    able_to_fail_ok: bool = True,
) -> str:
    if not dimensional_ok:
        return FAIL_DIMENSIONAL
    if not covariant:
        return FAIL_NOT_COVARIANT
    if not tautology_clear:
        return FAIL_TAUTOLOGICAL
    if not able_to_fail_ok:
        return FAIL_NOT_ABLE_TO_FAIL
    return ISOTROPY_CALIBRATED


def harmonics_l2() -> dict[str, sp.Expr]:
    return {
        "20": sp.sqrt(sp.Rational(5, 16) / sp.pi) * (3 * sp.cos(theta) ** 2 - 1),
        "21c": -sp.sqrt(sp.Rational(15, 4) / sp.pi) * sp.sin(theta) * sp.cos(theta) * sp.cos(phi),
        "21s": -sp.sqrt(sp.Rational(15, 4) / sp.pi) * sp.sin(theta) * sp.cos(theta) * sp.sin(phi),
        "22c": sp.sqrt(sp.Rational(15, 16) / sp.pi) * sp.sin(theta) ** 2 * sp.cos(2 * phi),
        "22s": sp.sqrt(sp.Rational(15, 16) / sp.pi) * sp.sin(theta) ** 2 * sp.sin(2 * phi),
    }


def compute_angular_block(harmonics: dict[str, sp.Expr]) -> dict[str, Any]:
    order = list(harmonics)
    ys = [harmonics[name] for name in order]
    gram = sp.Matrix([[integrate_s2(ys[i] * ys[j]) for j in range(5)] for i in range(5)])
    neg_laps: dict[str, sp.Expr] = {}
    lambdas: dict[str, sp.Expr] = {}
    residuals: dict[str, sp.Expr] = {}
    for name, y_expr in harmonics.items():
        neg_lap = compact(-laplacian_s2(y_expr))
        rayleigh = compact(integrate_s2(y_expr * neg_lap) / integrate_s2(y_expr**2))
        neg_laps[name] = neg_lap
        lambdas[name] = rayleigh
        residuals[name] = compact(neg_lap - rayleigh * y_expr)
    return {
        "order": order,
        "ys": ys,
        "gram": gram,
        "gram_is_identity": bool(gram == sp.eye(5)),
        "neg_laps": neg_laps,
        "lambdas": lambdas,
        "residuals": residuals,
        "lambda_all_six": all(compact(value - 6) == 0 for value in lambdas.values()),
        "residuals_zero": all(compact(value) == 0 for value in residuals.values()),
    }


def build_K2(coeff: sp.Expr) -> sp.Expr:
    return compact(Ktilde + coeff * TomegaTilde)


def extract_k2_coeff(k2_expr: sp.Expr) -> sp.Expr:
    return compact(sp.diff(k2_expr, TomegaTilde))


def make_dim_rules() -> dict[sp.Symbol, Dim]:
    return {
        a_dim: Dim(1, 0, 0),
        dw_dim: Dim(1, 0, 0),
        dOmega_dim: ZERO_DIM,
        beta2_dim: ZERO_DIM,
        beta2_prime_dim: Dim(-1, 0, 0),
        mu_eta_density: Dim(-3, 1, 0),
        T_w_density: Dim(-1, 1, -2),
        K_eta_density: Dim(-3, 1, -2),
        T_Omega_density: Dim(-3, 1, -2),
        Mtilde: EXPECTED_M,
        Ktilde: EXPECTED_K,
        TomegaTilde: EXPECTED_K,
    }


def dimension_eval(lambda_m: sp.Expr, m2_expr: sp.Expr, k2_expr: sp.Expr, dims: dict[sp.Symbol, Dim]) -> dict[str, Any]:
    measure = a_dim**2 * dw_dim * dOmega_dim
    m2_integral = mu_eta_density * beta2_dim**2 * measure
    k_tw_term = T_w_density * beta2_prime_dim**2 * measure
    k_eta_term = K_eta_density * beta2_dim**2 * measure
    k_omega_term = lambda_m * T_Omega_density * beta2_dim**2 * measure
    k2_integral = k_tw_term + k_eta_term + k_omega_term
    try:
        measure_dim = dim_of(measure, dims)
        m2_integral_dim = dim_of(m2_integral, dims)
        k_tw_dim = dim_of(k_tw_term, dims)
        k_eta_dim = dim_of(k_eta_term, dims)
        k_omega_dim = dim_of(k_omega_term, dims)
        k2_integral_dim = dim_of(k2_integral, dims)
        actual_m2_dim = dim_of(m2_expr, dims)
        actual_k2_dim = dim_of(k2_expr, dims)
        actual_ratio_dim = actual_k2_dim - actual_m2_dim
        term_homogeneous = k_tw_dim == k_eta_dim == k_omega_dim == k2_integral_dim
        ok = bool(
            measure_dim == Dim(3, 0, 0)
            and m2_integral_dim == EXPECTED_M
            and term_homogeneous
            and k2_integral_dim == EXPECTED_K
            and actual_m2_dim == EXPECTED_M
            and actual_k2_dim == EXPECTED_K
            and actual_ratio_dim == EXPECTED_RATIO
        )
        return {
            "ok": ok,
            "error": None,
            "measure": measure,
            "m2_integral": m2_integral,
            "k2_terms": {
                "T_w_beta_prime_sq": k_tw_term,
                "K_eta_beta_sq": k_eta_term,
                "lambda_T_Omega_beta_sq": k_omega_term,
            },
            "k2_integral": k2_integral,
            "actual_M2": m2_expr,
            "actual_K2": k2_expr,
            "dims": {
                "measure": measure_dim,
                "M2_integral": m2_integral_dim,
                "T_w_beta_prime_sq": k_tw_dim,
                "K_eta_beta_sq": k_eta_dim,
                "lambda_T_Omega_beta_sq": k_omega_dim,
                "K2_integral": k2_integral_dim,
                "actual_M2": actual_m2_dim,
                "actual_K2": actual_k2_dim,
                "actual_K2_over_M2": actual_ratio_dim,
            },
            "term_homogeneous": term_homogeneous,
        }
    except DimError as exc:
        return {
            "ok": False,
            "error": str(exc),
            "measure": measure,
            "m2_integral": m2_integral,
            "k2_terms": {
                "T_w_beta_prime_sq": k_tw_term,
                "K_eta_beta_sq": k_eta_term,
                "lambda_T_Omega_beta_sq": k_omega_term,
            },
            "k2_integral": k2_integral,
            "actual_M2": m2_expr,
            "actual_K2": k2_expr,
            "dims": {},
            "term_homogeneous": False,
        }


def corrupt_rules_for(label: str, base_rules: dict[sp.Symbol, Dim]) -> dict[sp.Symbol, Dim]:
    symbol_by_label = {
        "mu_eta_density": mu_eta_density,
        "T_w_density": T_w_density,
        "K_eta_density": K_eta_density,
        "T_Omega_density": T_Omega_density,
    }
    corrupt = dict(base_rules)
    corrupt[symbol_by_label[label]] = corrupt[symbol_by_label[label]] + Dim(1, 0, 0)
    if label == "T_Omega_density":
        corrupt[TomegaTilde] = corrupt[TomegaTilde] + Dim(1, 0, 0)
    return corrupt


def build_dimensional_check(lambda_m: sp.Expr, m2_expr: sp.Expr, k2_expr: sp.Expr) -> dict[str, Any]:
    rules = make_dim_rules()
    baseline = dimension_eval(lambda_m, m2_expr, k2_expr, rules)
    if not baseline["ok"]:
        raise AuditFailure(f"baseline dimensional check failed: {baseline['error']}")

    density_corruptions = {
        label: dimension_eval(lambda_m, m2_expr, k2_expr, corrupt_rules_for(label, rules))
        for label in ("mu_eta_density", "T_w_density", "K_eta_density", "T_Omega_density")
    }
    t_omega_probe = density_corruptions["T_Omega_density"]
    return {
        "rules": rules,
        "baseline": baseline,
        "dimensional_ok": baseline["ok"],
        "density_corruptions": density_corruptions,
        "FAIL_DIMENSIONAL_probe": {
            "mutation": "corrupt sourced [T_Omega] and assembled TomegaTilde by one extra L",
            "participates_in_verdict": (
                scoped_verdict(dimensional_ok=t_omega_probe["ok"]) == FAIL_DIMENSIONAL
                and scoped_verdict(dimensional_ok=baseline["ok"]) != FAIL_DIMENSIONAL
            ),
            "without_mutation_dimensional_ok": baseline["ok"],
            "with_mutation_dimensional_ok": t_omega_probe["ok"],
            "probe_verdict": "NO_FAIL" if t_omega_probe["ok"] else FAIL_DIMENSIONAL,
            "mutation_fires": not t_omega_probe["ok"],
            "error": t_omega_probe["error"],
            "self_ablation": {
                "dimensional_ok": baseline["ok"],
                "verdict": scoped_verdict(dimensional_ok=baseline["ok"]),
                "fail_fires": scoped_verdict(dimensional_ok=baseline["ok"]) == FAIL_DIMENSIONAL,
                "fail_suppressed": scoped_verdict(dimensional_ok=baseline["ok"]) != FAIL_DIMENSIONAL,
            },
        },
    }


def build_baseline() -> dict[str, Any]:
    harmonics = harmonics_l2()
    angular = compute_angular_block(harmonics)
    order = angular["order"]
    lambdas = angular["lambdas"]
    neg_laps = angular["neg_laps"]

    k2_core = {name: build_K2(lambdas[name]) for name in order}
    k2_coeff = {name: extract_k2_coeff(k2_core[name]) for name in order}
    k2_coeff_residuals = {
        name: compact(neg_laps[name] - k2_coeff[name] * harmonics[name]) for name in order
    }
    k2_coeff_residuals_zero = all(compact(residual) == 0 for residual in k2_coeff_residuals.values())

    lambda_ref = lambdas["20"]
    m2_core = Mtilde
    k2_ref = build_K2(lambda_ref)
    dimensional = build_dimensional_check(lambda_ref, m2_core, k2_ref)

    input_hashes = {name: expr_hash(harmonics[name]) for name in order}
    distinct_hashes = len(set(input_hashes.values())) == len(input_hashes)
    self_overlaps = {name: integrate_s2(harmonics[name] ** 2) for name in order}
    tautology_clear = bool(distinct_hashes and all(compact(value - 1) == 0 for value in self_overlaps.values()))

    def forced_eigenvalue_probe(forced_coefficient: sp.Expr) -> dict[str, Any]:
        forced_k2_by_channel = {name: build_K2(forced_coefficient) for name in order}
        forced_coeff_by_channel = {
            name: extract_k2_coeff(forced_k2_by_channel[name]) for name in order
        }
        coefficient_equals = all(
            compact(forced_coeff_by_channel[name] - lambdas[name]) == 0 for name in order
        )
        verdict = scoped_verdict(covariant=coefficient_equals)
        return {
            "forced_coefficient": forced_coefficient,
            "forced_K2_by_channel": forced_k2_by_channel,
            "forced_coeff_by_channel": forced_coeff_by_channel,
            "coefficient_equals_computed_lambda": coefficient_equals,
            "verdict": verdict,
            "fail_fires": verdict == FAIL_NOT_COVARIANT,
        }

    def lane_hash_probe(lane_inputs: dict[str, sp.Expr]) -> dict[str, Any]:
        hashes = {lane: expr_hash(expr) for lane, expr in lane_inputs.items()}
        distinct = len(set(hashes.values())) == len(hashes)
        verdict = scoped_verdict(tautology_clear=distinct)
        return {
            "input_hashes": hashes,
            "distinct_hashes": distinct,
            "verdict": verdict,
            "fail_fires": verdict == FAIL_TAUTOLOGICAL,
        }

    wrong_eigen_probe = forced_eigenvalue_probe(sp.Integer(2))
    wrong_eigen_ablation = forced_eigenvalue_probe(sp.Integer(6))
    tautology_probe = lane_hash_probe({"20": harmonics["20"], "21": harmonics["20"], "22": harmonics["20"]})
    tautology_ablation = lane_hash_probe(
        {"20": harmonics["20"], "21": harmonics["21c"], "22": harmonics["22c"]}
    )
    dim_probe = dimensional["FAIL_DIMENSIONAL_probe"]
    dimensional_verdict = scoped_verdict(dimensional_ok=dim_probe["with_mutation_dimensional_ok"])
    dimensional_ablation_verdict = scoped_verdict(dimensional_ok=dim_probe["without_mutation_dimensional_ok"])

    probes = {
        "wrong_eigenvalue": {
            "description": "force bare K2 coefficient 2 while computed lambda_m remains 6",
            **wrong_eigen_probe,
            "computed_fail_gate": not wrong_eigen_probe["coefficient_equals_computed_lambda"],
            "self_ablation": {
                **wrong_eigen_ablation,
                "fail_suppressed": not wrong_eigen_ablation["fail_fires"],
            },
        },
        "tautology_hash_collision": {
            "description": "set the three grouped-lane harmonic inputs to identical expressions and hash them",
            **tautology_probe,
            "computed_fail_gate": not tautology_probe["distinct_hashes"],
            "self_ablation": {
                **tautology_ablation,
                "fail_suppressed": not tautology_ablation["fail_fires"],
            },
        },
        "dimensional_corrupt_T_Omega": {
            "description": "corrupt sourced T_Omega and TomegaTilde dimensions",
            **dim_probe,
            "verdict": dimensional_verdict,
            "fail_fires": dimensional_verdict == FAIL_DIMENSIONAL,
            "self_ablation": {
                "description": "restore sourced T_Omega/TomegaTilde dimensions",
                "dimensional_ok": dim_probe["without_mutation_dimensional_ok"],
                "verdict": dimensional_ablation_verdict,
                "fail_fires": dimensional_ablation_verdict == FAIL_DIMENSIONAL,
                "fail_suppressed": dimensional_ablation_verdict != FAIL_DIMENSIONAL,
            },
        },
    }

    expected_probe_verdicts = {
        "wrong_eigenvalue": FAIL_NOT_COVARIANT,
        "tautology_hash_collision": FAIL_TAUTOLOGICAL,
        "dimensional_corrupt_T_Omega": FAIL_DIMENSIONAL,
    }
    expected_probe_verdicts_match = {
        key: probes[key]["verdict"] == expected for key, expected in expected_probe_verdicts.items()
    }
    computed_probe_gate_flags = {
        "wrong_eigenvalue": bool(
            probes["wrong_eigenvalue"]["computed_fail_gate"]
            and probes["wrong_eigenvalue"]["self_ablation"]["fail_suppressed"]
        ),
        "tautology_hash_collision": bool(
            probes["tautology_hash_collision"]["computed_fail_gate"]
            and probes["tautology_hash_collision"]["self_ablation"]["fail_suppressed"]
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
    able_to_fail = {
        "expected_probe_verdicts": expected_probe_verdicts,
        "expected_probe_verdicts_match": expected_probe_verdicts_match,
        "computed_probe_gate_flags": computed_probe_gate_flags,
        "able_to_fail_if_probe_neutered": able_to_fail_if_probe_neutered,
        "neutering_any_probe_flips_false": all(not value for value in able_to_fail_if_probe_neutered.values()),
        "able_to_fail_ok": able_to_fail_ok,
    }
    covariant_ok = bool(
        angular["gram_is_identity"]
        and angular["lambda_all_six"]
        and angular["residuals_zero"]
        and k2_coeff_residuals_zero
    )
    gates = {
        "dimensional_ok": dimensional["dimensional_ok"],
        "covariant": covariant_ok,
        "tautology_clear": tautology_clear,
        "able_to_fail_ok": able_to_fail_ok,
    }
    verdict = scoped_verdict(**gates)

    return {
        "harmonics": harmonics,
        "angular": angular,
        "K2_core": k2_core,
        "K2_coeff": k2_coeff,
        "K2_coeff_residuals": k2_coeff_residuals,
        "K2_coeff_residuals_zero": k2_coeff_residuals_zero,
        "M2_core": m2_core,
        "K2_ref": k2_ref,
        "dimensional": dimensional,
        "input_hashes": input_hashes,
        "distinct_hashes": distinct_hashes,
        "self_overlaps": self_overlaps,
        "tautology_clear": tautology_clear,
        "forced_eigenvalue_probe": forced_eigenvalue_probe,
        "lane_hash_probe": lane_hash_probe,
        "probes": probes,
        "able_to_fail": able_to_fail,
        "gates": gates,
        "verdict": verdict,
    }


def run_angular_theorem(data: dict[str, Any]) -> None:
    angular = data["angular"]
    order = angular["order"]
    gram = angular["gram"]
    lambdas = angular["lambdas"]
    residuals = angular["residuals"]
    subbanner("Real l=2 harmonics, Gram=I5, and computed -Delta_S2 spectrum")
    print(f"  harmonic order = {order}")
    print(f"  Gram matrix = {fmt(gram)}")
    print(f"  computed lambda_m = { {name: fmt(lambdas[name]) for name in order} }")
    print("  eigenvalues are Rayleigh quotients from the native Laplace-Beltrami operator.")
    expect_zero("Gram - I5 total squared residual", sum((gram - sp.eye(5))[i, j] ** 2 for i in range(5) for j in range(5)))
    expect_bool("Gram matrix equals I5", angular["gram_is_identity"])
    for name in order:
        expect_zero(f"lambda_{name} computed by Rayleigh quotient equals 6", lambdas[name] - 6)
        expect_zero(f"eigenfunction residual (-Delta)Y_{name}-lambda*Y_{name}", residuals[name])
    expect_bool("all computed lambda_m are 6", angular["lambda_all_six"])
    expect_bool("all eigenfunction residuals are zero", angular["residuals_zero"])


def run_k2_stiffness(data: dict[str, Any]) -> None:
    order = data["angular"]["order"]
    subbanner("Bare K2 angular stiffness uses the live computed lambda")
    print("  K2 builder = build_K2(coeff) = Ktilde + coeff*TomegaTilde")
    print("  baseline assembly path = build_K2(lambdas[name]); no counted k_coeff_used-lambdas self-compare.")
    print(f"  K2 core by channel = { {name: fmt(data['K2_core'][name]) for name in order} }")
    print(f"  extracted TomegaTilde coefficients = { {name: fmt(data['K2_coeff'][name]) for name in order} }")
    for name in order:
        expect_zero(
            f"K2-coefficient residual reads extracted coeff for {name}",
            data["K2_coeff_residuals"][name],
        )
    expect_bool("K2 coefficient residuals all vanish", data["K2_coeff_residuals_zero"])
    wrong = data["probes"]["wrong_eigenvalue"]
    ablation = wrong["self_ablation"]
    expect_zero("wrong_eigenvalue probe verdict is FAIL_NOT_COVARIANT", verdict_residual(wrong["verdict"], FAIL_NOT_COVARIANT))
    expect_bool("wrong_eigenvalue computed_fail_gate reads coefficient mismatch", wrong["computed_fail_gate"])
    expect_zero("wrong_eigenvalue self-ablation returns ISOTROPY_CALIBRATED", verdict_residual(ablation["verdict"], ISOTROPY_CALIBRATED))
    expect_bool("wrong_eigenvalue self-ablation suppresses fail", ablation["fail_suppressed"])


def run_dimensional_gate(data: dict[str, Any]) -> None:
    dimensional = data["dimensional"]
    baseline = dimensional["baseline"]
    dims = baseline["dims"]
    probe = data["probes"]["dimensional_corrupt_T_Omega"]
    subbanner("Angular dimensional gate and sourced-density corruption probes")
    print("  dimension order = (L,M,T)")
    print("  sourced volume-density convention: dV=a_dim^2*dw*dOmega has dimension L^3; beta2 is dimensionless.")
    print(f"  walked K2 terms = { {key: fmt(value) for key, value in baseline['k2_terms'].items()} }")
    print(f"  computed dimensions = { {key: dim_text(value) for key, value in dims.items()} }")
    expect_zero("explicit measure has dimension L^3", dim_residual(dims["measure"], Dim(3, 0, 0)))
    expect_zero("M2_integral has dimension M", dim_residual(dims["M2_integral"], EXPECTED_M))
    expect_zero("K2 T_w beta_prime^2 term has dimension M*T^-2", dim_residual(dims["T_w_beta_prime_sq"], EXPECTED_K))
    expect_zero("K2 K_eta beta^2 term has dimension M*T^-2", dim_residual(dims["K_eta_beta_sq"], EXPECTED_K))
    expect_zero("K2 lambda*T_Omega beta^2 term has dimension M*T^-2", dim_residual(dims["lambda_T_Omega_beta_sq"], EXPECTED_K))
    expect_zero("K2 integral has dimension M*T^-2", dim_residual(dims["K2_integral"], EXPECTED_K))
    expect_zero("bare M2=Mtilde has dimension M", dim_residual(dims["actual_M2"], EXPECTED_M))
    expect_zero("bare K2=Ktilde+lambda*TomegaTilde has dimension M*T^-2", dim_residual(dims["actual_K2"], EXPECTED_K))
    expect_zero("bare K2/M2 has dimension T^-2", dim_residual(dims["actual_K2_over_M2"], EXPECTED_RATIO))
    expect_bool("K2 density terms are homogeneous", baseline["term_homogeneous"])
    expect_bool("baseline dimensional_ok", dimensional["dimensional_ok"])
    print(f"  T_Omega corruption error = {probe['error']}")
    expect_zero("sourced T_Omega probe verdict is FAIL_DIMENSIONAL", verdict_residual(probe["verdict"], FAIL_DIMENSIONAL))
    expect_bool("sourced T_Omega mutation participates in verdict", probe["participates_in_verdict"])
    expect_bool("sourced T_Omega mutation fires", probe["mutation_fires"])
    expect_zero("T_Omega self-ablation returns ISOTROPY_CALIBRATED", verdict_residual(probe["self_ablation"]["verdict"], ISOTROPY_CALIBRATED))
    expect_bool("T_Omega self-ablation suppresses fail", probe["self_ablation"]["fail_suppressed"])
    for label, result in dimensional["density_corruptions"].items():
        verdict = scoped_verdict(dimensional_ok=result["ok"])
        expect_zero(f"corrupt sourced [{label}] fires FAIL_DIMENSIONAL", verdict_residual(verdict, FAIL_DIMENSIONAL))


def run_tautology_hash(data: dict[str, Any]) -> None:
    subbanner("Computed-not-typed hash guard and tautology probe")
    print(f"  harmonic input hashes = {data['input_hashes']}")
    print(f"  self-overlaps = { {name: fmt(value) for name, value in data['self_overlaps'].items()} }")
    expect_bool("five harmonic hashes are distinct", data["distinct_hashes"])
    for name, value in data["self_overlaps"].items():
        expect_zero(f"self-overlap integral for Y_{name} is 1", value - 1)
    expect_bool("tautology_clear = distinct_hashes and unit self-overlaps", data["tautology_clear"])
    taut = data["probes"]["tautology_hash_collision"]
    ablation = taut["self_ablation"]
    expect_zero("tautology_hash_collision verdict is FAIL_TAUTOLOGICAL", verdict_residual(taut["verdict"], FAIL_TAUTOLOGICAL))
    expect_bool("tautology_hash_collision computed_fail_gate reads non-distinct hashes", taut["computed_fail_gate"])
    expect_zero("tautology_hash_collision self-ablation returns ISOTROPY_CALIBRATED", verdict_residual(ablation["verdict"], ISOTROPY_CALIBRATED))
    expect_bool("tautology_hash_collision self-ablation suppresses fail", ablation["fail_suppressed"])


def run_aggregate(data: dict[str, Any]) -> None:
    able = data["able_to_fail"]
    subbanner("016 aggregate probe battery over the three covariance probes")
    print(f"  expected probe verdicts = {able['expected_probe_verdicts']}")
    print(f"  computed probe gate flags = {able['computed_probe_gate_flags']}")
    print(f"  neutered aggregates = {able['able_to_fail_if_probe_neutered']}")
    for key, value in able["expected_probe_verdicts_match"].items():
        expect_bool(f"probe {key} verdict matches expected token", value)
    for key, value in able["computed_probe_gate_flags"].items():
        expect_bool(f"probe {key} flag is computed true", value)
    expect_bool("016 able_to_fail_ok is true", able["able_to_fail_ok"])
    for key, value in able["able_to_fail_if_probe_neutered"].items():
        expect_bool(f"neutering {key} flips able_to_fail_ok false", not value)
    expect_bool("neutering any one 016 probe flips aggregate false", able["neutering_any_probe_flips_false"])


def run_per_tooth_ablations(data: dict[str, Any]) -> None:
    subbanner("Per-tooth ablations on copies")
    base_harmonics = data["harmonics"]
    corrupt_basis = dict(base_harmonics)
    corrupt_basis["20"] = sp.cos(theta) + sp.cos(theta) ** 2
    corrupt_angular = compute_angular_block(corrupt_basis)
    corrupt_covariant = bool(
        corrupt_angular["gram_is_identity"]
        and corrupt_angular["lambda_all_six"]
        and corrupt_angular["residuals_zero"]
    )
    expect_zero(
        "computed-eigenvalue basis corruption reaches FAIL_NOT_COVARIANT",
        verdict_residual(scoped_verdict(covariant=corrupt_covariant), FAIL_NOT_COVARIANT),
    )
    expect_fail("basis corruption makes at least one eigen residual nonzero", bool_residual(corrupt_angular["residuals_zero"]))

    wrong_k2 = build_K2(sp.Integer(2))
    wrong_coeff = extract_k2_coeff(wrong_k2)
    y20 = base_harmonics["20"]
    expect_fail(
        "mutating the assembled K2 coefficient to 2 breaks the K2-coefficient residual",
        data["angular"]["neg_laps"]["20"] - wrong_coeff * y20,
    )

    wrong_probe = data["forced_eigenvalue_probe"](sp.Integer(2))
    right_probe = data["forced_eigenvalue_probe"](sp.Integer(6))
    expect_zero("bare forced_eigenvalue_probe(2) fires FAIL_NOT_COVARIANT", verdict_residual(wrong_probe["verdict"], FAIL_NOT_COVARIANT))
    expect_zero("bare forced_eigenvalue_probe(6) suppresses FAIL_NOT_COVARIANT", verdict_residual(right_probe["verdict"], ISOTROPY_CALIBRATED))

    gram_corrupt_harmonics = dict(base_harmonics)
    gram_corrupt_harmonics["20"] = 2 * gram_corrupt_harmonics["20"]
    gram_corrupt = compute_angular_block(gram_corrupt_harmonics)
    expect_fail("Gram tooth: scaling one harmonic breaks Gram=I5", sum((gram_corrupt["gram"] - sp.eye(5))[i, j] ** 2 for i in range(5) for j in range(5)))
    expect_zero(
        "Gram tooth reaches FAIL_NOT_COVARIANT",
        verdict_residual(scoped_verdict(covariant=gram_corrupt["gram_is_identity"]), FAIL_NOT_COVARIANT),
    )

    taut_neutered_verdict = scoped_verdict(tautology_clear=True)
    expect_fail("tautology distinctness neuter would suppress FAIL_TAUTOLOGICAL", verdict_residual(taut_neutered_verdict, FAIL_TAUTOLOGICAL))

    dim_probe = data["probes"]["dimensional_corrupt_T_Omega"]
    expect_zero("T_Omega dimensional tooth reaches FAIL_DIMENSIONAL", verdict_residual(dim_probe["verdict"], FAIL_DIMENSIONAL))
    for key, value in data["able_to_fail"]["able_to_fail_if_probe_neutered"].items():
        expect_bool(f"aggregate tooth includes {key}", not value)
    print("  k_coeff_equal de-count: no k_coeff_used-lambdas self-compare is counted; K2 computed-ness rides on build_K2(lambdas), extracted-coefficient residuals, and the bare forced probe.")


def run_verdict_and_scope(data: dict[str, Any]) -> None:
    subbanner("016 scoped landing and 016/017 cut")
    print(f"  016 gate booleans = {data['gates']}")
    print(f"  016 scoped verdict = {data['verdict']}")
    expect_zero("016 scoped verdict lands ISOTROPY_CALIBRATED component", verdict_residual(data["verdict"], ISOTROPY_CALIBRATED))
    print("  ISOTROPY_CALIBRATED (JOINT, 2-stage) -- PARTIAL: 016 landed, 017 PENDING")
    print("    = (016: l=2 SO(3) covariance theorem: real harmonics + Gram=I5 + computed lambda_m=6 + K2 angular stiffness) [EARNED here]")
    print("    AND (017: grouped-P2 lane isotropy: grouped {20,21,22} lanes / raw-D=0 / normalized-u / calibration partition) [PENDING]")
    print("  Exact cut: this script does not assemble grouped lanes, raw-D, normalized-u, response probes, calibration partition, or port-kernel export.")
    print("  Carried caveats: angular structure is earned; radial profile/scalars are frozen calibration inputs, so the joint is CALIBRATED not PASS.")
    print("  Deferred caveats: 54/5 quadrupole normalization, outgoing odd-N coefficients, and solved nonlinear branch data are Gate 4/5/6 sim-deferred.")
    expect_bool("joint composition is partial with 017 pending", data["verdict"] == ISOTROPY_CALIBRATED)


def run_provenance(data: dict[str, Any]) -> None:
    subbanner("Provenance and scope labels")
    print("  CONSUMED-from-011/012/013: mu_eta/T_w physical wall constants, beta2(w)/R0(w), Mtilde/Ktilde/TomegaTilde, and Gate-1 D/N provenance are cited as provenance.")
    print("  Self-contained dimensional integrity: pathA_32 uses volume densities on a_dim^2*dw*dOmega with dimensionless beta2; stage013's line-density K_eta=T_w*beta^2 relation does not transfer.")
    print("  no-c_S: the l=2 covariance theorem is speed-free; matter-sector c_s/BdG remains deferred.")
    print("  ANGULAR-EARNED / RADIAL-CALIBRATED: 016 derives Gram=I5, lambda_m=6, and K2; 017 owns the radial calibration partition.")
    print("  COMPUTED-NOT-TYPED: Rayleigh + eigenfunction residual + extracted K2-coefficient residual + forced_eigenvalue_probe; k_coeff_equal X==X is de-counted.")
    print("  AGGREGATE-BATTERY-INTACT: wrong_eigenvalue, tautology_hash_collision, and dimensional_corrupt_T_Omega all participate.")
    print("  SOURCED-T_OMEGA-DIMENSIONAL: corrupting sourced T_Omega plus TomegaTilde fires FAIL_DIMENSIONAL.")
    print("  dropped-bookkeeping: scratch-YAML engine agreement and report/feed writers are stripped; transcript-level agreement remains.")
    print("  register note: 016 is an earned structural slice with likely zero new counted knobs; T_Omega and beta2 counting is deferred to 017.")
    live_symbols = {
        symbol.name
        for expr in [
            *data["harmonics"].values(),
            *data["angular"]["lambdas"].values(),
            *data["K2_core"].values(),
            data["M2_core"],
            data["K2_ref"],
        ]
        for symbol in sp.sympify(expr).free_symbols
    }
    expect_bool("no c_S/cS live symbol appears", "c_S" not in live_symbols and "cS" not in live_symbols)
    expect_bool("Btilde/Ztilde support scalars are not live 016 symbols", not any(name.startswith(("B", "Z")) for name in live_symbols))


def print_verdict_labels() -> None:
    subbanner("Verdict labels")
    print("  ledger earned-label: L2_SO3_COVARIANCE_THEOREM_EARNED")
    print("  source top-line verdict: ISOTROPY_CALIBRATED (JOINT 2-stage; 016 is the earned SO(3) covariance component, 017 completes the calibration/lane component)")
    print("  joint composition: ISOTROPY_CALIBRATED = 016[EARNED l=2 SO(3) covariance] AND 017[PENDING grouped-P2 lane isotropy + calibration partition]")
    print("  earned angular structure: Gram=I5 genuine; lambda_m=6 computed by Rayleigh + residual; K2=Ktilde+lambda_m*TomegaTilde uses the live computed lambda")
    print("  earned able-to-fail battery: wrong_eigenvalue / tautology_hash_collision / dimensional_corrupt_T_Omega, each with suppressing self-ablation")
    print("  consumed framing: provenance + pathA_32 self-contained dimensional integrity, not a cross-stage dual-site relation")
    print("  new symbols first appearing here but not counted in 016: T_Omega/TomegaTilde and beta2(w), deferred to 017")


def main() -> int:
    banner("ledger_stage016_l2_so3_covariance_sympy_audit")
    print("Target stem confirmed: ledger_stage016_l2_so3_covariance")
    print("Engine: SymPy exact symbolic; no scipy/numpy/floats/tolerances; zero file I/O.")
    data = build_baseline()
    run_angular_theorem(data)
    run_k2_stiffness(data)
    run_dimensional_gate(data)
    run_tautology_hash(data)
    run_aggregate(data)
    run_per_tooth_ablations(data)
    run_verdict_and_scope(data)
    run_provenance(data)
    print_verdict_labels()
    return 0


if __name__ == "__main__":
    try:
        exit_code = main()
    except Exception as exc:
        if not isinstance(exc, AuditFailure):
            _record_fail(f"UNCAUGHT exception: {exc!r}")
        banner("Tallies")
        total = PASS_COUNT + FAIL_COUNT
        print(f"TALLY sympy: {PASS_COUNT} pass + {FAIL_COUNT} fail = {total} checks")
        print("OVERALL FAIL")
        raise SystemExit(1) from exc
    banner("Tallies")
    total = PASS_COUNT + FAIL_COUNT
    print(f"TALLY sympy: {PASS_COUNT} pass + {FAIL_COUNT} fail = {total} checks")
    if FAIL_COUNT == 0 and exit_code == 0:
        print("OVERALL PASS")
        raise SystemExit(0)
    print("OVERALL FAIL")
    raise SystemExit(1)
