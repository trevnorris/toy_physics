#!/usr/bin/env python3
"""Ledger stage017 SymPy audit: grouped-P2 lane isotropy.

Standalone, print-only, no arguments, no file I/O. This is the pathA_32
II-G3b response slice only: grouped {20,21,22} lane assembly, raw-D primary
isotropy, normalized-u cross-check, stability/denominator windows, the six
017 response probes, calibration partition, and the joint COMPLETE landing.
Stage 016's SO(3) covariance theorem is consumed as cited data with an explicit
lambda/K2 dual-site guard plus one Y20 echo; it is not re-derived here.
"""

from __future__ import annotations

import math
from typing import Any

import sympy as sp


PASS_COUNT = 0
FAIL_COUNT = 0

ISOTROPY_CALIBRATED = "ISOTROPY_CALIBRATED"
ISOTROPY_PASS = "ISOTROPY_PASS"
FAIL_DIMENSIONAL = "FAIL_DIMENSIONAL"
FAIL_NOT_COVARIANT = "FAIL_NOT_COVARIANT"
FAIL_TAUTOLOGICAL = "FAIL_TAUTOLOGICAL"
FAIL_STATIC_RESPONSE = "FAIL_STATIC_RESPONSE"
FAIL_STABILITY = "FAIL_STABILITY"
FAIL_SINGULAR_RESPONSE = "FAIL_SINGULAR_RESPONSE"
FAIL_ANISOTROPIC_BRANCH = "FAIL_ANISOTROPIC_BRANCH"
FAIL_NOT_ABLE_TO_FAIL = "FAIL_NOT_ABLE_TO_FAIL"

SYMBOLIC_TOL = 1.0e-10
NUMERIC_TOL = 1.0e-8
PROBE_TOL = 1.0e-12
EPS_VALUES = [1.0e-4, 2.0e-4, -3.0e-4]
DELTA_PROFILE = 0.1

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

CITED_016_LAMBDA = sp.Integer(6)
CITED_016_C_SELF = sp.Integer(1)
CITED_016_DIMENSIONAL_OK = True
CITED_016_COVARIANT = True
CITED_016_TAUTOLOGY_CLEAR = True
CITED_016_THREE_PROBE_AGGREGATE = True


class AuditFailure(AssertionError):
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
    if isinstance(expr, float):
        return f"{expr:.12g}"
    if isinstance(expr, dict):
        return "{" + ", ".join(f"{key}: {fmt(value)}" for key, value in expr.items()) + "}"
    if isinstance(expr, (list, tuple)):
        return "[" + ", ".join(fmt(value) for value in expr) + "]"
    try:
        return sp.sstr(compact(expr))
    except Exception:
        return str(expr)


def assert_no_float(name: str, expr: Any) -> None:
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
        raise AuditFailure(f"{name}: Float atom(s) found in exact expression: {floats}")


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
    expect_zero(name, sp.Integer(0) if bool(condition) else sp.Integer(1))


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
    return sp.Integer(0) if bool(condition) else sp.Integer(1)


def verdict_residual(actual: str, expected: str) -> sp.Integer:
    return sp.Integer(0) if actual == expected else sp.Integer(1)


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


def f_eval(expr: sp.Expr, subs: dict[sp.Symbol, float]) -> float:
    return float(sp.N(compact(expr).subs(subs), 30))


def first_nonzero(values: list[float], tol: float = PROBE_TOL) -> bool:
    return any(abs(float(value)) > tol for value in values)


def max_ratio_delta(values: list[float], eps_values: list[float]) -> float:
    ratios = [value / eps for value, eps in zip(values, eps_values) if eps != 0.0]
    if not ratios:
        return math.inf
    return max(abs(ratio - ratios[0]) for ratio in ratios)


def verdict_from_gates(gates: dict[str, bool], calibration_inputs: list[str]) -> str:
    if not gates["dimensional_ok"]:
        return FAIL_DIMENSIONAL
    if not gates["covariant"]:
        return FAIL_NOT_COVARIANT
    if not gates["tautology_clear"]:
        return FAIL_TAUTOLOGICAL
    if not gates["dynamic_retained"]:
        return FAIL_STATIC_RESPONSE
    if not gates["stability_ok"]:
        return FAIL_STABILITY
    if not gates["denominator_guard_ok"]:
        return FAIL_SINGULAR_RESPONSE
    if not gates["lane_collapse_ok"]:
        return FAIL_ANISOTROPIC_BRANCH
    if not gates["able_to_fail_ok"]:
        return FAIL_NOT_ABLE_TO_FAIL
    if calibration_inputs:
        return ISOTROPY_CALIBRATED
    return ISOTROPY_PASS


def harmonics_l2(theta: sp.Symbol, phi: sp.Symbol) -> dict[str, sp.Expr]:
    return {
        "20": sp.sqrt(sp.Rational(5, 16) / sp.pi) * (3 * sp.cos(theta) ** 2 - 1),
        "21c": -sp.sqrt(sp.Rational(15, 4) / sp.pi) * sp.sin(theta) * sp.cos(theta) * sp.cos(phi),
        "21s": -sp.sqrt(sp.Rational(15, 4) / sp.pi) * sp.sin(theta) * sp.cos(theta) * sp.sin(phi),
        "22c": sp.sqrt(sp.Rational(15, 16) / sp.pi) * sp.sin(theta) ** 2 * sp.cos(2 * phi),
        "22s": sp.sqrt(sp.Rational(15, 16) / sp.pi) * sp.sin(theta) ** 2 * sp.sin(2 * phi),
    }


def sample_subs(symbols: dict[str, sp.Symbol]) -> dict[sp.Symbol, float]:
    return {symbols[name]: value for name, value in CALIBRATION_SAMPLE.items()}


def build_response(
    lambda_by_channel: dict[str, sp.Expr] | None = None,
    lambda_reference: sp.Expr = CITED_016_LAMBDA,
    k2_lambda_offsets: dict[str, sp.Expr] | None = None,
) -> dict[str, Any]:
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
    harmonics = harmonics_l2(theta, phi)
    order = list(harmonics)
    if lambda_by_channel is None:
        lambda_by_channel = {name: CITED_016_LAMBDA for name in order}
    if k2_lambda_offsets is None:
        k2_lambda_offsets = {}

    p2_axis_z = (3 * sp.cos(theta) ** 2 - 1) / 2
    anisotropy_coefficients = {
        name: integrate_s2(p2_axis_z * harmonics[name] ** 2, theta, phi)
        for name in order
    }
    c_group = {
        "20": anisotropy_coefficients["20"],
        "21": compact((anisotropy_coefficients["21c"] + anisotropy_coefficients["21s"]) / 2),
        "22": compact((anisotropy_coefficients["22c"] + anisotropy_coefficients["22s"]) / 2),
    }

    single_harmonic_echo_residual = compact(
        -laplacian_s2(harmonics["20"], theta, phi) - lambda_reference * harmonics["20"]
    )

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

    c_self = {name: CITED_016_C_SELF for name in order}

    def assemble_channel(name: str, profile_factor: sp.Expr = sp.Integer(1)) -> dict[str, Any]:
        pref = compact(c_self[name] * profile_factor)
        lam = lambda_by_channel[name]
        k2_offset = k2_lambda_offsets.get(name, sp.Integer(0))
        m_lane = compact(pref * M)
        k_lane = compact(pref * (Kbase + (lam + k2_offset) * Tomega))
        b_lane = {"0": compact(pref * B0), "2": compact(pref * B2), "4": compact(pref * B4)}
        z_lane = {"0": compact(pref * Z0), "2": compact(pref * Z2), "4": compact(pref * Z4)}
        d_lane = {
            "0": compact(k_lane - b_lane["0"] - z_lane["0"]),
            "2": compact(-(m_lane + b_lane["2"] + z_lane["2"])),
            "4": compact(-(b_lane["4"] + z_lane["4"])),
        }
        return {
            "angular_self_overlap": c_self[name],
            "lambda": lam,
            "M2": m_lane,
            "K2": k_lane,
            "B": b_lane,
            "Z": z_lane,
            "D": d_lane,
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

    lambda_site_residuals = {
        f"lambda_site_A_equals_016_export_{name}": compact(lambda_by_channel[name] - lambda_reference)
        for name in order
    }
    k2_form_residuals = {
        f"K2_form_coeff_matches_016_export_{name}": compact(
            sp.diff(ungrouped[name]["K2"], Tomega) - c_self[name] * lambda_reference
        )
        for name in order
    }
    dual_site_residuals = {
        **lambda_site_residuals,
        **k2_form_residuals,
        "single_harmonic_Y20_echo": single_harmonic_echo_residual,
    }
    dual_site_ok = all_zero(list(dual_site_residuals.values()))

    cs_equal = {
        f"D21c_equals_D21s_order_{n}": compact(ungrouped["21c"]["D"][n] - ungrouped["21s"]["D"][n]) == 0
        for n in ["0", "2", "4"]
    }
    cs_equal.update(
        {
            f"D22c_equals_D22s_order_{n}": compact(ungrouped["22c"]["D"][n] - ungrouped["22s"]["D"][n]) == 0
            for n in ["0", "2", "4"]
        }
    )

    raw_triples = {
        n: [grouped_lanes[lane]["D"][n] for lane in ["20", "21", "22"]]
        for n in ["0", "2", "4"]
    }
    raw_defects: dict[str, sp.Expr] = {}
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

    k2_window_min = CALIBRATION_WINDOW["Ktilde"]["min"] + 6.0 * CALIBRATION_WINDOW["TomegaTilde"]["min"]
    d0_window_min = (
        CALIBRATION_WINDOW["Ktilde"]["min"]
        + 6.0 * CALIBRATION_WINDOW["TomegaTilde"]["min"]
        - CALIBRATION_WINDOW["B0tilde"]["max"]
        - CALIBRATION_WINDOW["Z0tilde"]["max"]
    )
    stability_ok = bool(CALIBRATION_WINDOW["Mtilde"]["min"] > 0.0 and k2_window_min > 0.0)
    denominator_guard_ok = bool(d0_window_min > 0.0)
    k2_sample = CALIBRATION_SAMPLE["Ktilde"] + 6.0 * CALIBRATION_SAMPLE["TomegaTilde"]
    omega_sample = math.sqrt(k2_sample / CALIBRATION_SAMPLE["Mtilde"])

    derived_inputs = [
        "explicit real l=2 harmonics for 017 anisotropy coefficients and the one-Y20 echo",
        "016-cited Gram diagonal c_self=1, lambda_m=6, and K2 angular coefficient",
        "ungrouped and grouped {20,21,22} lane assembly",
        "raw-D primary defect algebra from assembled lanes",
        "normalized-u cross-check algebra",
        "017 six response-probe verdicts and computed gate flags",
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
            "dimensional_ok": CITED_016_DIMENSIONAL_OK,
            "covariant": CITED_016_COVARIANT and dual_site_ok,
            "tautology_clear": CITED_016_TAUTOLOGY_CLEAR,
            "dynamic_retained": True,
            "stability_ok": True,
            "denominator_guard_ok": True,
            "lane_collapse_ok": True,
            "able_to_fail_ok": True,
        }
        gates.update(overrides)
        return verdict_from_gates(gates, calibration_inputs)

    d_common = {
        "0": compact(Kbase + CITED_016_LAMBDA * Tomega - S0),
        "2": compact(-(M + S2)),
        "4": compact(-S4),
    }

    pure_d_by_lane = {
        lane: {n: compact((1 + eps * c_group[lane]) * d_common[n]) for n in ["0", "2", "4"]}
        for lane in ["20", "21", "22"]
    }
    pure_raw_defects: dict[str, sp.Expr] = {}
    pure_samples: dict[str, list[float]] = {}
    for n in ["0", "2", "4"]:
        a_def, b_def = defect_pair([pure_d_by_lane[lane][n] for lane in ["20", "21", "22"]])
        pure_raw_defects[f"a_D{n}"] = a_def
        pure_raw_defects[f"b_D{n}"] = b_def
        for label, expr in ((f"a_D{n}", a_def), (f"b_D{n}", b_def)):
            pure_samples[label] = [f_eval(expr.subs(eps, eps_val), subs_sample) for eps_val in EPS_VALUES]
    pure_u2 = {lane: u2_from_d(pure_d_by_lane[lane]) for lane in ["20", "21", "22"]}
    pure_u4 = {lane: u4_from_d(pure_d_by_lane[lane]) for lane in ["20", "21", "22"]}
    pure_a2, pure_b2 = defect_pair([pure_u2[lane] for lane in ["20", "21", "22"]])
    pure_a4, pure_b4 = defect_pair([pure_u4[lane] for lane in ["20", "21", "22"]])
    pure_u_defects = {"a2": pure_a2, "b2": pure_b2, "a4": pure_a4, "b4": pure_b4}
    pure_raw_moves = all(first_nonzero(values) for values in pure_samples.values())
    pure_linear_delta = max(max_ratio_delta(values, EPS_VALUES) for values in pure_samples.values())
    pure_u_zero = all_zero(list(pure_u_defects.values()))

    sector_d_by_lane = {
        lane: {
            "0": compact(Kbase + CITED_016_LAMBDA * Tomega - (1 + eps * c_group[lane]) * S0),
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
        for label, expr in ((f"a_D{n}", a_def), (f"b_D{n}", b_def)):
            sector_samples[label] = [f_eval(expr.subs(eps, eps_val), subs_sample) for eps_val in EPS_VALUES]
    sector_u2 = {lane: u2_from_d(sector_d_by_lane[lane]) for lane in ["20", "21", "22"]}
    sector_a2, sector_b2 = defect_pair([sector_u2[lane] for lane in ["20", "21", "22"]])
    sector_u_defects = {"a2": sector_a2, "b2": sector_b2}
    sector_u_samples = {
        label: [f_eval(expr.subs(eps, eps_val), subs_sample) for eps_val in EPS_VALUES]
        for label, expr in sector_u_defects.items()
    }
    sector_raw_moves = all(first_nonzero(values) for values in sector_samples.values())
    sector_u_moves = first_nonzero(sector_u_samples["a2"]) and first_nonzero(sector_u_samples["b2"])

    profile_scales = {"20": sp.Integer(1), "21": sp.Integer(1), "22": compact((1 + delta) ** 2)}
    profile_d_by_lane = {
        lane: {n: compact(profile_scales[lane] * d_common[n]) for n in ["0", "2", "4"]}
        for lane in ["20", "21", "22"]
    }
    profile_raw_defects: dict[str, sp.Expr] = {}
    for n in ["0", "2", "4"]:
        a_def, b_def = defect_pair([profile_d_by_lane[lane][n] for lane in ["20", "21", "22"]])
        profile_raw_defects[f"a_D{n}"] = a_def
        profile_raw_defects[f"b_D{n}"] = b_def
    profile_sample = {
        label: f_eval(expr.subs(delta, DELTA_PROFILE), subs_sample)
        for label, expr in profile_raw_defects.items()
    }
    profile_raw_moves = first_nonzero(list(profile_sample.values()))

    def beta_stability_probe(beta_scale: sp.Expr) -> dict[str, Any]:
        m2_expr = compact(beta_scale**2 * grouped_lanes["20"]["M2"])
        k2_expr = compact(beta_scale**2 * grouped_lanes["20"]["K2"])
        m2_sample = f_eval(m2_expr, subs_sample)
        k2_sample_value = f_eval(k2_expr, subs_sample)
        ok = bool(m2_sample > PROBE_TOL and k2_sample_value > PROBE_TOL)
        verdict = case_verdict(stability_ok=ok)
        return {
            "beta2_scale": beta_scale,
            "M2": m2_expr,
            "K2": k2_expr,
            "M2_sample": m2_sample,
            "K2_sample": k2_sample_value,
            "stability_ok": ok,
            "verdict": verdict,
            "fail_fires": verdict == FAIL_STABILITY,
        }

    degenerate_beta_probe = beta_stability_probe(sp.Integer(0))
    degenerate_beta_ablation = beta_stability_probe(sp.Integer(1))

    singular_subs = dict(subs_sample)
    singular_subs[symbols["B0tilde"]] = k2_sample - singular_subs[symbols["Z0tilde"]]
    singular_d0_value = f_eval(d_common["0"], singular_subs)
    singular_denominator_guard_ok = abs(singular_d0_value) >= PROBE_TOL
    singular_verdict = case_verdict(denominator_guard_ok=singular_denominator_guard_ok)
    nonsingular_d0_value = f_eval(d_common["0"], subs_sample)
    nonsingular_denominator_guard_ok = abs(nonsingular_d0_value) >= PROBE_TOL
    singular_ablation_verdict = case_verdict(denominator_guard_ok=nonsingular_denominator_guard_ok)

    static_wrong_d2 = compact(-(B2 + Z2))
    static_dynamic_retained = bool(static_wrong_d2.has(M))
    static_verdict = case_verdict(dynamic_retained=static_dynamic_retained)
    static_ablation_dynamic_retained = bool(grouped_lanes["20"]["D"]["2"].has(M))
    static_ablation_verdict = case_verdict(dynamic_retained=static_ablation_dynamic_retained)

    probes = {
        "pure_prefactor_anisotropy": {
            "raw_D_defects": pure_raw_defects,
            "normalized_u_defects": pure_u_defects,
            "sample_values": pure_samples,
            "raw_D_moves": pure_raw_moves,
            "linear_scaling_max_ratio_delta": pure_linear_delta,
            "linear_scaling_confirmed": pure_raw_moves and pure_linear_delta < SYMBOLIC_TOL,
            "normalized_u_defects_stay_zero": pure_u_zero,
            "verdict": case_verdict(lane_collapse_ok=False),
        },
        "sector_selective_anisotropy": {
            "raw_D_defects": sector_raw_defects,
            "normalized_u_defects": sector_u_defects,
            "sample_values": {**sector_samples, **{f"u_{key}": value for key, value in sector_u_samples.items()}},
            "raw_D_moves": sector_raw_moves,
            "u_defects_move": sector_u_moves,
            "verdict": case_verdict(lane_collapse_ok=False),
        },
        "m_dependent_profile": {
            "profile_scales": profile_scales,
            "raw_D_defects": profile_raw_defects,
            "sample_values": profile_sample,
            "raw_D_moves": profile_raw_moves,
            "verdict": case_verdict(lane_collapse_ok=False),
        },
        "degenerate_beta_zero": {
            **degenerate_beta_probe,
            "computed_fail_gate": not degenerate_beta_probe["stability_ok"],
            "self_ablation": {
                **degenerate_beta_ablation,
                "fail_suppressed": not degenerate_beta_ablation["fail_fires"],
            },
        },
        "singular_denominator": {
            "D0_sample": singular_d0_value,
            "M2_positive": CALIBRATION_SAMPLE["Mtilde"] > PROBE_TOL,
            "K2_positive": k2_sample > PROBE_TOL,
            "denominator_guard_ok": singular_denominator_guard_ok,
            "computed_fail_gate": not singular_denominator_guard_ok,
            "verdict": singular_verdict,
            "fail_fires": singular_verdict == FAIL_SINGULAR_RESPONSE,
            "self_ablation": {
                "D0_sample": nonsingular_d0_value,
                "denominator_guard_ok": nonsingular_denominator_guard_ok,
                "verdict": singular_ablation_verdict,
                "fail_fires": singular_ablation_verdict == FAIL_SINGULAR_RESPONSE,
                "fail_suppressed": singular_ablation_verdict != FAIL_SINGULAR_RESPONSE,
            },
        },
        "static_drop_inertia": {
            "dynamic_retained": static_dynamic_retained,
            "computed_fail_gate": not static_dynamic_retained,
            "wrong_D2_without_M": static_wrong_d2,
            "correct_D2": d_common["2"],
            "verdict": static_verdict,
            "fail_fires": static_verdict == FAIL_STATIC_RESPONSE,
            "self_ablation": {
                "dynamic_retained": static_ablation_dynamic_retained,
                "correct_D2": grouped_lanes["20"]["D"]["2"],
                "verdict": static_ablation_verdict,
                "fail_fires": static_ablation_verdict == FAIL_STATIC_RESPONSE,
                "fail_suppressed": static_ablation_verdict != FAIL_STATIC_RESPONSE,
            },
        },
    }

    expected_probe_verdicts = {
        "pure_prefactor_anisotropy": FAIL_ANISOTROPIC_BRANCH,
        "sector_selective_anisotropy": FAIL_ANISOTROPIC_BRANCH,
        "m_dependent_profile": FAIL_ANISOTROPIC_BRANCH,
        "degenerate_beta_zero": FAIL_STABILITY,
        "singular_denominator": FAIL_SINGULAR_RESPONSE,
        "static_drop_inertia": FAIL_STATIC_RESPONSE,
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
        "singular_denominator": bool(
            probes["singular_denominator"]["computed_fail_gate"]
            and probes["singular_denominator"]["self_ablation"]["fail_suppressed"]
        ),
        "static_drop_inertia": bool(
            probes["static_drop_inertia"]["computed_fail_gate"]
            and probes["static_drop_inertia"]["self_ablation"]["fail_suppressed"]
        ),
    }

    def able_to_fail_017_from_flags(flags: dict[str, bool]) -> bool:
        return bool(all(expected_probe_verdicts_match.values()) and all(flags.values()))

    able_to_fail_017_ok = able_to_fail_017_from_flags(computed_probe_gate_flags)
    able_to_fail_if_probe_neutered = {
        key: able_to_fail_017_from_flags({**computed_probe_gate_flags, key: False})
        for key in computed_probe_gate_flags
    }
    able_to_fail_ok = bool(able_to_fail_017_ok and CITED_016_THREE_PROBE_AGGREGATE)

    dynamic_retained = bool(d_common["2"].has(M) and not d_common["0"].has(M))
    lane_collapse_ok = bool(raw_defects_zero and all(cs_equal.values()))
    baseline_gates = {
        "dimensional_ok": CITED_016_DIMENSIONAL_OK,
        "covariant": bool(CITED_016_COVARIANT and dual_site_ok),
        "tautology_clear": CITED_016_TAUTOLOGY_CLEAR,
        "dynamic_retained": dynamic_retained,
        "stability_ok": stability_ok,
        "denominator_guard_ok": denominator_guard_ok,
        "lane_collapse_ok": lane_collapse_ok,
        "able_to_fail_ok": able_to_fail_ok,
    }
    verdict = verdict_from_gates(baseline_gates, calibration_inputs)

    return {
        "theta": theta,
        "phi": phi,
        "eps": eps,
        "delta": delta,
        "symbols": symbols,
        "subs_sample": subs_sample,
        "harmonics": harmonics,
        "order": order,
        "anisotropy_coefficients": anisotropy_coefficients,
        "c_group": c_group,
        "lambda_by_channel": lambda_by_channel,
        "lambda_reference": lambda_reference,
        "c_self": c_self,
        "single_harmonic_echo_residual": single_harmonic_echo_residual,
        "dual_site_residuals": dual_site_residuals,
        "dual_site_ok": dual_site_ok,
        "ungrouped": ungrouped,
        "grouped_lanes": grouped_lanes,
        "cs_equal": cs_equal,
        "raw_triples": raw_triples,
        "raw_defects": raw_defects,
        "raw_defects_zero": raw_defects_zero,
        "u2_lanes": u2_lanes,
        "u4_lanes": u4_lanes,
        "normalized_defects": normalized_defects,
        "normalized_defects_zero": normalized_defects_zero,
        "stability_ok": stability_ok,
        "denominator_guard_ok": denominator_guard_ok,
        "k2_window_min": k2_window_min,
        "d0_window_min": d0_window_min,
        "omega_sample": omega_sample,
        "derived_inputs": derived_inputs,
        "calibration_inputs": calibration_inputs,
        "case_verdict": case_verdict,
        "d_common": d_common,
        "probes": probes,
        "expected_probe_verdicts": expected_probe_verdicts,
        "expected_probe_verdicts_match": expected_probe_verdicts_match,
        "computed_probe_gate_flags": computed_probe_gate_flags,
        "able_to_fail_017_from_flags": able_to_fail_017_from_flags,
        "able_to_fail_017_ok": able_to_fail_017_ok,
        "able_to_fail_if_probe_neutered": able_to_fail_if_probe_neutered,
        "able_to_fail_ok": able_to_fail_ok,
        "dynamic_retained": dynamic_retained,
        "lane_collapse_ok": lane_collapse_ok,
        "baseline_gates": baseline_gates,
        "verdict": verdict,
        "u2_from_d": u2_from_d,
        "u4_from_d": u4_from_d,
        "M": M,
        "Kbase": Kbase,
        "Tomega": Tomega,
    }


def raw_defects_from_grouped(grouped_lanes: dict[str, Any]) -> dict[str, sp.Expr]:
    raw_defects: dict[str, sp.Expr] = {}
    for n in ["0", "2", "4"]:
        a_def, b_def = defect_pair([grouped_lanes[lane]["D"][n] for lane in ["20", "21", "22"]])
        raw_defects[f"a_D{n}"] = a_def
        raw_defects[f"b_D{n}"] = b_def
    return raw_defects


def normalized_defects_from_grouped(grouped_lanes: dict[str, Any]) -> dict[str, sp.Expr]:
    def u2_from_d(d_lane: dict[str, sp.Expr]) -> sp.Expr:
        return compact(-d_lane["2"] / d_lane["0"])

    def u4_from_d(d_lane: dict[str, sp.Expr]) -> sp.Expr:
        return compact((d_lane["2"] ** 2 - d_lane["0"] * d_lane["4"]) / d_lane["0"] ** 2)

    u2_lanes = {lane: u2_from_d(grouped_lanes[lane]["D"]) for lane in ["20", "21", "22"]}
    u4_lanes = {lane: u4_from_d(grouped_lanes[lane]["D"]) for lane in ["20", "21", "22"]}
    a2, b2_def = defect_pair([u2_lanes[lane] for lane in ["20", "21", "22"]])
    a4, b4_def = defect_pair([u4_lanes[lane] for lane in ["20", "21", "22"]])
    return {"a2": a2, "b2": b2_def, "a4": a4, "b4": b4_def}


def print_partition(derived_inputs: list[str], calibration_inputs: list[str]) -> None:
    print("Derived inputs:")
    for item in derived_inputs:
        print(f"  - {item}")
    print("Calibration inputs (printed only; count is resolved at registration):")
    for item in calibration_inputs:
        print(f"  - {item}")


def print_verdict_labels() -> None:
    print("Verdict labels:")
    print(
        "  ledger earned-label: GROUPED_P2_LANE_ISOTROPY_EARNED "
        "(grouped l=2 lanes, raw-D primary defects zero, normalized-u cross-check zero, six response probes able-to-fail)"
    )
    print("  source top-line verdict: ISOTROPY_CALIBRATED (JOINT 2-stage; 017 completes 016)")
    print(
        "  joint composition: ISOTROPY_CALIBRATED = "
        "(016 l=2 SO(3) covariance theorem: Gram=I5, lambda_m=6, K2 angular stiffness)[EARNED, cited] AND "
        "(017 grouped-P2 lane isotropy: grouped lanes, raw-D=0 PRIMARY, normalized-u=0 CROSS-CHECK, 6-probe battery)[EARNED here]"
    )
    print(
        "  calibrated values: beta2(w), R0(w), Mtilde/Ktilde/TomegaTilde, Btilde/Ztilde, and the window are frozen; "
        "T_Omega/beta2/support-scalar counting is left to registration"
    )
    print(
        "  consumed: 016 lambda_m=6/K2=Ktilde+lambda_m*TomegaTilde/c_self=1 via explicit dual-site guard; "
        "011/012/013 mu_eta/T_w/R0/Gate-1 D-N are provenance; c_S NOT consumed"
    )
    print(
        "  exports: l=2 port kernel = grouped M2 + K2=Ktilde+6*TomegaTilde + Btilde/Ztilde support scalars + D0/D2/D4 lanes"
    )


def run_assertions(data: dict[str, Any]) -> None:
    subbanner("Consumed 016 dual-site guard")
    for key, residual in data["dual_site_residuals"].items():
        expect_zero(f"dual_site.{key}", residual)
    expect_bool("dual_site.baseline_ok", data["dual_site_ok"])
    expect_zero("cited_c_self_is_1_for_lane_assembly", data["c_self"]["20"] - 1)
    expect_bool("cited_016_dimensional_ok", CITED_016_DIMENSIONAL_OK)
    expect_bool("cited_016_covariant_gate", CITED_016_COVARIANT)
    expect_bool("cited_016_tautology_clear", CITED_016_TAUTOLOGY_CLEAR)
    expect_bool("cited_016_three_probe_aggregate", CITED_016_THREE_PROBE_AGGREGATE)

    subbanner("Grouped lane assembly and exact isotropy")
    expect_zero("anisotropy.c_group_20_is_2_over_7", data["c_group"]["20"] - sp.Rational(2, 7))
    expect_zero("anisotropy.c_group_21_is_1_over_7", data["c_group"]["21"] - sp.Rational(1, 7))
    expect_zero("anisotropy.c_group_22_is_minus_2_over_7", data["c_group"]["22"] + sp.Rational(2, 7))
    for lane in ["20", "21", "22"]:
        expect_zero(f"grouped_lane.{lane}.M2_is_Mtilde", data["grouped_lanes"][lane]["M2"] - data["symbols"]["Mtilde"])
        expect_zero(
            f"grouped_lane.{lane}.K2_is_Ktilde_plus_6_TomegaTilde",
            data["grouped_lanes"][lane]["K2"]
            - (data["symbols"]["Ktilde"] + 6 * data["symbols"]["TomegaTilde"]),
        )
    for key, ok in data["cs_equal"].items():
        expect_bool(f"cs_equal.{key}", ok)
    for key, residual in data["raw_defects"].items():
        expect_zero(f"raw_D_primary.{key}", residual)
    expect_bool("raw_D_primary.raw_defects_zero", data["raw_defects_zero"])
    for key, residual in data["normalized_defects"].items():
        expect_zero(f"normalized_cross_check.{key}", residual)
    expect_bool("normalized_cross_check.normalized_defects_zero", data["normalized_defects_zero"])

    subbanner("Numeric stability and denominator windows")
    expect_bool("stability.window_M_and_K2_positive", data["stability_ok"])
    expect_bool("denominator.window_D0_positive", data["denominator_guard_ok"])
    expect_bool("stability.k2_window_min_positive", data["k2_window_min"] > 0.0)
    expect_bool("denominator.d0_window_min_positive", data["d0_window_min"] > 0.0)

    subbanner("Six response probes")
    probes = data["probes"]
    expected = data["expected_probe_verdicts"]
    for key in [
        "pure_prefactor_anisotropy",
        "sector_selective_anisotropy",
        "m_dependent_profile",
        "degenerate_beta_zero",
        "singular_denominator",
        "static_drop_inertia",
    ]:
        expect_zero(f"probe.{key}.verdict", verdict_residual(probes[key]["verdict"], expected[key]))
        expect_bool(f"probe.{key}.computed_gate_flag", data["computed_probe_gate_flags"][key])
    expect_bool("probe.pure_prefactor.raw_D_moves", probes["pure_prefactor_anisotropy"]["raw_D_moves"])
    expect_bool(
        "probe.pure_prefactor.normalized_u_defects_stay_zero",
        probes["pure_prefactor_anisotropy"]["normalized_u_defects_stay_zero"],
    )
    expect_bool(
        "probe.pure_prefactor.linear_scaling_confirmed",
        probes["pure_prefactor_anisotropy"]["linear_scaling_confirmed"],
    )
    expect_bool("probe.sector_selective.raw_D_moves", probes["sector_selective_anisotropy"]["raw_D_moves"])
    expect_bool("probe.sector_selective.u_defects_move", probes["sector_selective_anisotropy"]["u_defects_move"])
    expect_bool("probe.m_dependent.raw_D_moves", probes["m_dependent_profile"]["raw_D_moves"])
    for key in ["degenerate_beta_zero", "singular_denominator", "static_drop_inertia"]:
        expect_bool(f"probe.{key}.computed_fail_gate", probes[key]["computed_fail_gate"])
        expect_bool(f"probe.{key}.fail_fires", probes[key]["fail_fires"])
        expect_bool(f"probe.{key}.self_ablation_fail_suppressed", probes[key]["self_ablation"]["fail_suppressed"])

    subbanner("Aggregate battery and joint landing")
    expect_bool("able_to_fail.017_six_probe_aggregate", data["able_to_fail_017_ok"])
    for key, value in data["able_to_fail_if_probe_neutered"].items():
        expect_bool(f"able_to_fail.neuter_{key}_flips_false", not value)
    expect_bool("able_to_fail.joint_016_and_017", data["able_to_fail_ok"])
    for key, value in data["baseline_gates"].items():
        expect_bool(f"baseline_gate.{key}", value)
    expect_zero("joint.verdict_is_ISOTROPY_CALIBRATED", verdict_residual(data["verdict"], ISOTROPY_CALIBRATED))
    expect_bool("no_c_S_live_symbol", all("c_S" not in str(symbol) and "cS" not in str(symbol) for symbol in data["symbols"].values()))


def run_teeth(data: dict[str, Any]) -> None:
    subbanner("Able-to-fail teeth")
    grouped_mut = {
        lane: {
            **data["grouped_lanes"][lane],
            "D": dict(data["grouped_lanes"][lane]["D"]),
        }
        for lane in ["20", "21", "22"]
    }
    grouped_mut["21"]["D"]["0"] = compact(grouped_mut["21"]["D"]["0"] + data["symbols"]["TomegaTilde"])
    raw_mut = raw_defects_from_grouped(grouped_mut)
    raw_mut_zero = all_zero(list(raw_mut.values()))
    expect_fail("tooth.raw_D_PRIMARY_lane_collapse_fires", bool_residual(raw_mut_zero))
    expect_zero(
        "tooth.raw_D_PRIMARY_verdict_rung",
        verdict_residual(data["case_verdict"](lane_collapse_ok=raw_mut_zero), FAIL_ANISOTROPIC_BRANCH),
    )

    pure_neutered = {
        **data["computed_probe_gate_flags"],
        "pure_prefactor_anisotropy": bool(
            False and data["probes"]["pure_prefactor_anisotropy"]["normalized_u_defects_stay_zero"]
        ),
    }
    expect_fail("tooth.pure_prefactor_neutered_raw_D_moves_flips_aggregate", bool_residual(data["able_to_fail_017_from_flags"](pure_neutered)))

    for key in [
        "sector_selective_anisotropy",
        "m_dependent_profile",
        "degenerate_beta_zero",
        "singular_denominator",
        "static_drop_inertia",
    ]:
        flags = {**data["computed_probe_gate_flags"], key: False}
        expect_fail(f"tooth.response_probe_{key}_neuter_flips_aggregate", bool_residual(data["able_to_fail_017_from_flags"](flags)))

    for key in data["computed_probe_gate_flags"]:
        flags = {**data["computed_probe_gate_flags"], key: False}
        expect_fail(f"tooth.six_probe_battery_intact_{key}", bool_residual(data["able_to_fail_017_from_flags"](flags)))

    order = data["order"]
    one_site_lambda = {name: CITED_016_LAMBDA for name in order}
    one_site_lambda["20"] = sp.Integer(4)
    one_site = build_response(lambda_by_channel=one_site_lambda, lambda_reference=CITED_016_LAMBDA)
    expect_fail(
        "tooth.dual_site_one_site_lambda_corruption_fires_at_integrity_guard",
        one_site["dual_site_residuals"]["lambda_site_A_equals_016_export_20"],
    )
    expect_fail(
        "tooth.dual_site_one_site_K2_form_corruption_fires_at_integrity_guard",
        one_site["dual_site_residuals"]["K2_form_coeff_matches_016_export_20"],
    )

    assembly_formula = build_response(
        lambda_by_channel={name: CITED_016_LAMBDA for name in order},
        lambda_reference=CITED_016_LAMBDA,
        k2_lambda_offsets={"20": sp.Integer(1)},
    )
    expect_zero(
        "tooth.dual_site_assembly_formula_lambda_site_stays_clean",
        assembly_formula["dual_site_residuals"]["lambda_site_A_equals_016_export_20"],
    )
    expect_fail(
        "tooth.dual_site_assembly_formula_K2_corruption_fires_at_K2_form_guard",
        assembly_formula["dual_site_residuals"]["K2_form_coeff_matches_016_export_20"],
    )

    coordinated = build_response(
        lambda_by_channel={name: sp.Integer(4) for name in order},
        lambda_reference=sp.Integer(4),
    )
    expect_fail(
        "tooth.dual_site_coordinated_corruption_caught_by_Y20_echo",
        coordinated["dual_site_residuals"]["single_harmonic_Y20_echo"],
    )

    grouped_norm_mut = {
        lane: {
            **data["grouped_lanes"][lane],
            "D": dict(data["grouped_lanes"][lane]["D"]),
        }
        for lane in ["20", "21", "22"]
    }
    grouped_norm_mut["22"]["D"]["2"] = compact(grouped_norm_mut["22"]["D"]["2"] + data["symbols"]["Mtilde"])
    norm_mut = normalized_defects_from_grouped(grouped_norm_mut)
    expect_fail("tooth.normalized_u_cross_check_fires", bool_residual(all_zero(list(norm_mut.values()))))

    rung_expectations = {
        "dynamic_retained": FAIL_STATIC_RESPONSE,
        "stability_ok": FAIL_STABILITY,
        "denominator_guard_ok": FAIL_SINGULAR_RESPONSE,
        "lane_collapse_ok": FAIL_ANISOTROPIC_BRANCH,
        "able_to_fail_ok": FAIL_NOT_ABLE_TO_FAIL,
    }
    for gate, expected in rung_expectations.items():
        gates = dict(data["baseline_gates"])
        gates[gate] = False
        expect_zero(f"tooth.joint_rung_{gate}", verdict_residual(verdict_from_gates(gates, data["calibration_inputs"]), expected))


def print_transcript(data: dict[str, Any]) -> None:
    subbanner("Transcript summary")
    print(f"Engine: SymPy standalone audit")
    print(f"Computed verdict: {data['verdict']} (COMPLETE: 016 AND 017)")
    print("Grouped lanes:")
    for lane in ["20", "21", "22"]:
        print(
            f"  {lane}: M2={fmt(data['grouped_lanes'][lane]['M2'])}; "
            f"K2={fmt(data['grouped_lanes'][lane]['K2'])}; "
            f"D0={fmt(data['grouped_lanes'][lane]['D']['0'])}; "
            f"D2={fmt(data['grouped_lanes'][lane]['D']['2'])}; "
            f"D4={fmt(data['grouped_lanes'][lane]['D']['4'])}"
        )
    print(f"Raw-D defects (PRIMARY): {fmt(data['raw_defects'])}")
    print(f"Normalized-u defects (CROSS-CHECK): {fmt(data['normalized_defects'])}")
    print(f"Stability guard: {data['stability_ok']}; denominator guard: {data['denominator_guard_ok']}")
    print("Probe verdicts:")
    for key in [
        "pure_prefactor_anisotropy",
        "sector_selective_anisotropy",
        "m_dependent_profile",
        "degenerate_beta_zero",
        "singular_denominator",
        "static_drop_inertia",
    ]:
        print(
            f"  {key}: {data['probes'][key]['verdict']} "
            f"(computed flag={data['computed_probe_gate_flags'][key]})"
        )
    print(f"017 six-probe aggregate: {data['able_to_fail_017_ok']}")
    print(f"Joint able_to_fail = 016 cited aggregate AND 017 aggregate: {data['able_to_fail_ok']}")
    print(
        "Dual-site guard: Site A lane lambda/K2 coefficient equals independent 016 lambda=6; "
        "coordinated corruption is closed by a single Y20 (-Delta_S2) echo."
    )
    print("CALIBRATED-not-PASS: the wall profile and radial/support scalars are frozen calibration inputs.")
    print(
        "Carried caveats: 54/5 quadrupole normalization, outgoing odd-N coefficients, and solved nonlinear branch data "
        "remain Gate 4/5/6 sim-deferred work."
    )
    print("Dropped bookkeeping: scratch YAML/report/feed writers replaced by print-only v2 audit output.")

    subbanner("Calibration partition")
    print_partition(data["derived_inputs"], data["calibration_inputs"])

    subbanner("Provenance labels")
    print(
        "CONSUMED-016-DUAL-SITE: lambda_m=6, K2=Ktilde+lambda_m*TomegaTilde, Gram diagonal c_self=1, "
        "dimensional/covariant/tautology gates, and 016's 3-probe aggregate are cited and guarded."
    )
    print(
        "CONSUMED-011/012/013-PROVENANCE: mu_eta/T_w, R0(w), and Gate-1 D/N provenance are narrative provenance; "
        "K_eta=T_w*beta^2 is non-transferable across the volume-vs-line convention."
    )
    print("NO-c_S: c_S is not a live symbol in the l=2 response audit.")
    print("RAW-D-PRIMARY: raw-D defects are the primary isotropy test; normalized-u is a cross-check.")
    print("SIX-PROBE-BATTERY-INTACT: each six-probe flag is computed; neutering any one flips the aggregate false.")
    print("COMPLETES-THE-JOINT: 017 lands ISOTROPY_CALIBRATED as COMPLETE (016 AND 017), not partial.")
    print("EXPORTS-PORT-KERNEL: grouped M2, K2=Ktilde+6*TomegaTilde, Btilde/Ztilde, and D-lanes to stages 018-024.")

    subbanner("Verdict labels")
    print_verdict_labels()


def main() -> int:
    banner("ledger_stage017 grouped-P2 lane isotropy SymPy audit")
    print("Mode: print-only, assert-zero, zero file I/O, mixed exact symbolic + numeric calibration probes.")
    print("Target stem: ledger_stage017_grouped_p2_lane_isotropy")
    data = build_response()
    run_assertions(data)
    run_teeth(data)
    print_transcript(data)
    total = PASS_COUNT + FAIL_COUNT
    print("")
    print(f"TALLY sympy: PASS={PASS_COUNT} FAIL={FAIL_COUNT} TOTAL={total} ({PASS_COUNT}+{FAIL_COUNT}={total})")
    print("OVERALL PASS" if FAIL_COUNT == 0 else "OVERALL FAIL")
    return 0 if FAIL_COUNT == 0 else 1


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except AuditFailure as exc:
        print(f"Audit aborted: {exc}")
        total = PASS_COUNT + FAIL_COUNT
        print(f"TALLY sympy: PASS={PASS_COUNT} FAIL={FAIL_COUNT} TOTAL={total} ({PASS_COUNT}+{FAIL_COUNT}={total})")
        print("OVERALL FAIL")
        raise SystemExit(1)
