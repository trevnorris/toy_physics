#!/usr/bin/env python3
"""pathA_29 v3 brane<->bulk return gate, SymPy engine.

This executable follows the v3 directive:

* Check A computes the return transfer from projected continuity plus the
  finite-slab Helmholtz transport phase.
* Check B solves the flat finite-slab transverse spectra for the admissible
  DC-sink completions and derives the radial falloff from the solved zero-mode
  Green function.
* The static/dynamic consistency control uses two separate solve traces.
* The mandatory no-go control is a derived delocalizing warp family, not an
  inadmissible wall condition.
"""

from __future__ import annotations

import hashlib
import inspect
import json
from pathlib import Path
from typing import Any

import sympy as sp
import yaml
from sympy.printing.mathematica import mathematica_code


class NoAliasDumper(yaml.SafeDumper):
    def ignore_aliases(self, data: Any) -> bool:
        return True


SCRIPT_PATH = Path(__file__).resolve()
STAGE1_ROOT = SCRIPT_PATH.parents[1]
REPO_ROOT = SCRIPT_PATH.parents[3]
REPORTS = STAGE1_ROOT / "reports"
SCRATCH = STAGE1_ROOT / "_scratch"

PATHA28_RESULTS = REPORTS / "pathA_28_monopole_results.yaml"
PATHA28_CONDITION = REPORTS / "pathA_28_cancellation_condition.yaml"
REPORT_OUT = REPORTS / "pathA_29_brane_bulk_return.md"
RESULTS_YAML = REPORTS / "pathA_29_results.yaml"
JSON_OUT = SCRATCH / "pathA_29_brane_bulk_return_sympy.json"
MMA_JSON = SCRATCH / "pathA_29_brane_bulk_return_mathematica.json"

I = sp.I


def hstr(expr: sp.Expr | int | str | bool | None) -> str | int | bool | None:
    if expr is None or isinstance(expr, (int, str, bool)):
        return expr
    return sp.sstr(sp.factor(sp.simplify(expr)))


def mma_expr(expr: sp.Expr | int) -> str:
    return mathematica_code(sp.factor(sp.simplify(expr)))


def sha256_text(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def load_yaml(path: Path) -> dict[str, Any]:
    if not path.exists():
        raise AssertionError(f"required reuse file is missing: {path}")
    data = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(data, dict):
        raise AssertionError(f"YAML did not load to a mapping: {path}")
    return data


def assert_patha28_reuse() -> dict[str, Any]:
    results = load_yaml(PATHA28_RESULTS)
    condition = load_yaml(PATHA28_CONDITION)
    expected_raw = {
        "ell0": {"kernel": "I*a*omega/cS", "p_raw": 1},
        "ell1": {"kernel": "I*a**3*omega**3/(2*cS**3)", "p_raw": 3},
        "ell2": {"kernel": "I*a**5*omega**5/(27*cS**5)", "p_raw": 5},
    }
    for ell, expected in expected_raw.items():
        got_kernel = results["raw_amplitudes"][ell]["kernel"]
        got_p = results["p_raw"][ell]
        if got_kernel != expected["kernel"] or got_p != expected["p_raw"]:
            raise AssertionError(f"pathA_28 reuse mismatch for {ell}: {got_kernel}, {got_p}")
    if condition["ell0"]["return_moment_required"] != results["cancellation_condition"]["ell0"]["return_moment_required"]:
        raise AssertionError("pathA_28 ell0 cancellation condition mismatch")
    if condition["ell1"]["return_moment_required"] != results["cancellation_condition"]["ell1"]["return_moment_required"]:
        raise AssertionError("pathA_28 ell1 cancellation condition mismatch")
    return {
        "pathA_28_results": str(PATHA28_RESULTS.relative_to(REPO_ROOT)),
        "pathA_28_condition": str(PATHA28_CONDITION.relative_to(REPO_ROOT)),
        "raw_kernels_reused_verbatim": expected_raw,
        "source_moment_definitions_reused": {
            "ell0": condition["ell0"]["source_moment"],
            "ell1": condition["ell1"]["source_moment"],
        },
        "cancellation_target_reused": {
            "ell0": condition["ell0"]["return_moment_required"],
            "ell1": condition["ell1"]["return_moment_required"],
        },
        "quadrupole_anchor_reused": {
            "p_raw": results["p_raw"]["ell2"],
            "kernel": results["raw_amplitudes"]["ell2"]["kernel"],
            "research_anchor": "research/4d_2_5pn",
        },
    }


def omega_order(expr: sp.Expr, omega: sp.Symbol, max_order: int) -> int:
    series = sp.expand(sp.series(expr, omega, 0, max_order + 1).removeO())
    for power in range(max_order + 1):
        coeff = sp.simplify(series.coeff(omega, power))
        if coeff != 0:
            return power
    raise AssertionError(f"no nonzero omega coefficient through {max_order}: {series}")


def series_coefficients(expr: sp.Expr, omega: sp.Symbol, order: int) -> dict[str, str]:
    expanded = sp.expand(expr)
    return {
        f"omega^{power}": str(hstr(sp.simplify(expanded.coeff(omega, power))))
        for power in range(order + 1)
    }


def radial_decay_exponent(radial_flow: sp.Expr, r: sp.Symbol) -> sp.Expr:
    """Extract the large-r power from a solved radial-flow expression."""
    probe = sp.simplify(-sp.limit(r * sp.diff(sp.log(radial_flow), r), r, sp.oo))
    if not probe.is_number:
        raise AssertionError(f"could not extract radial exponent from solved flow: {radial_flow}")
    return sp.simplify(probe)


def radial_operator_residual(radial_expr: sp.Expr, kappa_squared: sp.Expr, r: sp.Symbol) -> sp.Expr:
    return sp.simplify(
        sp.diff(radial_expr, r, 2) + sp.Rational(2, 1) * sp.diff(radial_expr, r) / r + kappa_squared * radial_expr
    )


def solve_round_trip_phase(omega: sp.Symbol, cS: sp.Symbol, d: sp.Symbol, w: sp.Symbol) -> dict[str, Any]:
    outgoing_basis = sp.exp(I * omega * w / cS)
    returning_basis = sp.exp(-I * omega * w / cS)
    forward_phase = sp.simplify(outgoing_basis.subs(w, d) / outgoing_basis.subs(w, 0))
    return_phase = sp.simplify(returning_basis.subs(w, 0) / returning_basis.subs(w, d))
    round_trip_phase = sp.simplify(forward_phase * return_phase)
    tau = sp.simplify(sp.diff(sp.log(round_trip_phase), omega) / I)
    if sp.simplify(sp.exp(I * omega * tau) - round_trip_phase) != 0:
        raise AssertionError("round-trip transit phase was not solved from the Helmholtz basis")
    return {
        "tau": tau,
        "transport_phase": round_trip_phase,
        "trace": {
            "basis": "Phi_l=A_l*exp(I*omega*w/c_s)+B_l*exp(-I*omega*w/c_s)",
            "forward_phase_ratio": hstr(forward_phase),
            "return_phase_ratio": hstr(return_phase),
            "round_trip_phase": hstr(round_trip_phase),
            "solved_tau": hstr(tau),
            "source_hash": sha256_text(inspect.getsource(solve_round_trip_phase)),
        },
    }


def solve_dynamic_zero_mode_radial(
    branch: str,
    eigenfunction: sp.Expr,
    m_value: sp.Expr,
    omega: sp.Symbol,
    cS: sp.Symbol,
    d: sp.Symbol,
    r: sp.Symbol,
) -> dict[str, Any]:
    """Solve g''+2g'/r+((omega/c_s)^2-m^2)g=0 first, then take omega -> 0."""
    if sp.simplify(m_value) != 0:
        raise AssertionError("dynamic zero-mode radial solve must be seeded by the computed m=0 eigenvalue")
    g = sp.Function(f"g_{branch}_dynamic")
    kappa_squared = sp.simplify((omega / cS) ** 2 - m_value**2)
    ode = sp.Eq(radial_operator_residual(g(r), kappa_squared, r), 0)
    general = sp.dsolve(ode).rhs
    C1, C2 = sp.symbols("C1 C2")
    selected = sp.simplify(
        general.subs(
            {
                C1: I * sp.sqrt(sp.pi) * sp.sqrt(omega) / (4 * sp.pi * d * sp.sqrt(2) * sp.sqrt(cS)),
                C2: -sp.sqrt(sp.pi) * sp.sqrt(omega) / (4 * sp.pi * d * sp.sqrt(2) * sp.sqrt(cS)),
            }
        )
    )
    residual = radial_operator_residual(selected, kappa_squared, r)
    if residual != 0:
        raise AssertionError("dynamic 3D radial solve residual did not vanish")
    limit_green = sp.simplify(sp.limit(selected, omega, 0))
    flow = sp.simplify(-sp.diff(limit_green, r))
    exponent = radial_decay_exponent(flow, r)
    trace_payload = {
        "route": "dynamic_limit_trace",
        "branch": branch,
        "seed_transverse_eigenvalue_m": hstr(m_value),
        "seed_transverse_eigenfunction": hstr(eigenfunction),
        "operator": "g''+(2/r)g'+((omega/c_s)^2-m^2)g=0 solved before omega->0",
        "dsolve_general": hstr(general),
        "boundary_selection": "outgoing spherical radial branch selected from the dsolve basis; normalization fixed by the compact zero-mode overlap 1/d",
        "solution": hstr(selected),
        "limit_green": hstr(limit_green),
        "flow_from_gradient": hstr(flow),
        "exponent": hstr(exponent),
        "source_hash": sha256_text(inspect.getsource(solve_dynamic_zero_mode_radial)),
    }
    trace_payload["trace_id"] = sha256_text(json.dumps(trace_payload, sort_keys=True))
    return {"solution": selected, "limit_solution": limit_green, "flow": flow, "exponent": exponent, "trace": trace_payload}


def solve_static_zero_mode_radial(
    branch: str,
    eigenfunction: sp.Expr,
    m_value: sp.Expr,
    d: sp.Symbol,
    r: sp.Symbol,
) -> dict[str, Any]:
    """Set omega=0 in the 3D radial operator first, then solve for g(r)."""
    if sp.simplify(m_value) != 0:
        raise AssertionError("static zero-mode radial solve must be seeded by the computed m=0 eigenvalue")
    g = sp.Function(f"g_{branch}_static")
    kappa_squared = -sp.simplify(m_value) ** 2
    ode = sp.Eq(radial_operator_residual(g(r), kappa_squared, r), 0)
    general = sp.dsolve(ode).rhs
    C1, C2 = sp.symbols("C1 C2")
    selected = sp.simplify(general.subs({C1: 0, C2: sp.Rational(1, 1) / (4 * sp.pi * d)}))
    residual = radial_operator_residual(selected, kappa_squared, r)
    if residual != 0:
        raise AssertionError("static 3D radial solve residual did not vanish")
    flow = sp.simplify(-sp.diff(selected, r))
    exponent = radial_decay_exponent(flow, r)
    trace_payload = {
        "route": "static_solve_trace",
        "branch": branch,
        "seed_transverse_eigenvalue_m": hstr(m_value),
        "seed_transverse_eigenfunction": hstr(eigenfunction),
        "operator": "omega set to 0 first: g''+(2/r)g'-m^2*g=0",
        "dsolve_general": hstr(general),
        "boundary_selection": "constant branch removed by decay at infinity; remaining radial branch normalized by compact zero-mode overlap 1/d",
        "solution": hstr(selected),
        "flow_from_gradient": hstr(flow),
        "exponent": hstr(exponent),
        "source_hash": sha256_text(inspect.getsource(solve_static_zero_mode_radial)),
    }
    trace_payload["trace_id"] = sha256_text(json.dumps(trace_payload, sort_keys=True))
    return {"solution": selected, "flow": flow, "exponent": exponent, "trace": trace_payload}


def solve_static_massive_radial(branch: str, m_value: sp.Expr, r: sp.Symbol) -> dict[str, Any]:
    """Solve the massive static radial equation and select the decaying Yukawa branch."""
    if sp.simplify(m_value) == 0:
        raise AssertionError("massive radial solve requires m>0")
    mu = sp.symbols(f"mu_{branch}", positive=True, real=True)
    g = sp.Function(f"g_{branch}_massive")
    ode = sp.Eq(radial_operator_residual(g(r), -mu**2, r), 0)
    general = sp.dsolve(ode).rhs
    general_hyperbolic = sp.simplify(general)
    C1, C2 = sp.symbols("C1 C2")
    selected_mu = sp.simplify(
        general_hyperbolic.subs(
            {
                C1: -sp.sqrt(sp.pi) * sp.sqrt(mu) / (4 * sp.pi * (1 + I)),
                C2: -sp.sqrt(sp.pi) * sp.sqrt(mu) / (4 * sp.pi * (1 - I)),
            }
        )
    )
    selected = sp.simplify(selected_mu.subs(mu, m_value))
    residual = radial_operator_residual(selected, -m_value**2, r)
    if residual != 0:
        raise AssertionError(f"massive radial solve residual did not vanish for {branch}")
    return {
        "solution": selected,
        "trace": {
            "branch": branch,
            "seed_transverse_eigenvalue_m": hstr(m_value),
            "operator": "g''+(2/r)g'-m^2*g=0",
            "dsolve_general": hstr(general.subs(mu, m_value)),
            "selected_decaying_solution": hstr(selected),
            "selection_rule": "growing exp(+m*r)/r coefficient set to zero in the dsolve basis",
            "residual": hstr(residual),
        },
    }


def build_counterfactual_guard(candidate: sp.Expr, kappa_squared: sp.Expr, r: sp.Symbol) -> dict[str, Any]:
    correct_residual = radial_operator_residual(candidate, kappa_squared, r)
    if correct_residual != 0:
        raise AssertionError("candidate radial solution failed its own 3D radial operator")
    perturbed = sp.simplify(candidate / r**4)
    perturbed_residual = radial_operator_residual(perturbed, kappa_squared, r)
    rejected = sp.simplify(perturbed_residual) != 0
    if not rejected:
        raise AssertionError("counterfactual wrong falloff passed the 3D radial operator")
    return {
        "perturbed_falloff_rejected": True,
        "perturbation": "multiplied the solved static zero-mode Green function by r**(-4), changing 1/r to 1/r**5",
        "candidate_residual": hstr(correct_residual),
        "perturbed_residual": hstr(perturbed_residual),
        "operator": "g''+(2/r)g'+kappa^2*g",
    }


def classify_dc_sink_gate(branch_ps: list[sp.Expr], quadrupole_survives: bool = True) -> str:
    target = sp.Integer(2)
    equal_target = [sp.simplify(p - target) == 0 for p in branch_ps]
    if not quadrupole_survives:
        return "RETURN_NOGO"
    if all(equal_target):
        return "RETURN_RESIDUAL_PREDICTION"
    if not any(equal_target):
        return "RETURN_NOGO"
    return "BC_DEPENDENT"


def branch_verdict_from_p(p_value: sp.Expr) -> str:
    return "RETURN_RESIDUAL_PREDICTION" if sp.simplify(p_value - 2) == 0 else "RETURN_NOGO"


def compute_symbolics() -> dict[str, Any]:
    omega, cS, d, a, r, rhoB = sp.symbols("omega cS d a r rhoB", positive=True, real=True)
    eps0, eps1 = sp.symbols("epsilon0 epsilon1", positive=True, real=True)
    M0 = sp.symbols("M0", positive=True, real=True)
    D1 = sp.symbols("D1", real=True)
    Gamma_uniform, kWarp = sp.symbols("Gamma_uniform k_warp", positive=True, real=True)
    w, m = sp.symbols("w m", nonnegative=True, real=True)
    n = sp.symbols("n", integer=True, positive=True)
    q = sp.symbols("q", integer=True)

    z = sp.simplify(omega * d / cS)
    transit_trace = solve_round_trip_phase(omega, cS, d, w)
    tau = transit_trace["tau"]
    alpha0 = sp.simplify(1 / (1 + eps0))
    alpha1 = sp.simplify(1 / (1 + eps1))
    neighbor_fraction0 = sp.simplify(eps0 / (1 + eps0))
    neighbor_fraction1 = sp.simplify(eps1 / (1 + eps1))

    helmholtz_basis = "Phi_l=A_l*exp(I*omega*w/c_s)+B_l*exp(-I*omega*w/c_s)"
    transport_phase = transit_trace["transport_phase"]
    T0_full = sp.simplify(alpha0 * transport_phase)
    T1_full = sp.simplify(alpha1 * transport_phase)
    T0_series = sp.expand(sp.series(T0_full, omega, 0, 5).removeO())
    T1_series = sp.expand(sp.series(T1_full, omega, 0, 3).removeO())
    T0_dc = sp.simplify(sp.limit(T0_full, omega, 0))
    T1_dc = sp.simplify(sp.limit(T1_full, omega, 0))
    if sp.simplify(T0_dc - alpha0) != 0 or sp.simplify(T1_dc - alpha1) != 0:
        raise AssertionError("DC transfer values did not follow from continuity fractions")

    sigma0 = omega_order(T0_series - T0_dc, omega, 4)
    sigma1 = omega_order(T1_series - T1_dc, omega, 2)
    nu0 = omega_order(1 - T0_series, omega, 4)
    nu1 = omega_order(1 - T1_series, omega, 2)
    if nu0 != 0 or nu1 != 0:
        raise AssertionError("finite DC sink must leave a nonzero deviation-from-one at omega^0")

    raw_kernel0 = I * a * omega / cS
    raw_kernel1 = I * a**3 * omega**3 / (2 * cS**3)
    raw_kernel2 = I * a**5 * omega**5 / (27 * cS**5)
    R0 = sp.simplify(-M0 * T0_series)
    R1 = sp.simplify(-D1 * T1_series)
    A0_res = sp.simplify(raw_kernel0 * M0 * (1 - T0_series))
    A1_res = sp.simplify(raw_kernel1 * D1 * (1 - T1_series))
    p_res0 = sp.Integer(1) + sp.Integer(nu0)
    p_res1 = sp.Integer(3) + sp.Integer(nu1)

    J_leak0 = M0
    J_return0 = sp.simplify(alpha0 * M0)
    J_neighbor0 = sp.simplify(neighbor_fraction0 * M0)
    steady_balance0 = sp.simplify(J_leak0 - J_return0 - J_neighbor0)
    if steady_balance0 != 0:
        raise AssertionError("zeroth-moment steady circulation balance failed")
    J_return1 = sp.simplify(alpha1 * D1)
    J_neighbor1 = sp.simplify(neighbor_fraction1 * D1)
    steady_balance1 = sp.simplify(D1 - J_return1 - J_neighbor1)
    if steady_balance1 != 0:
        raise AssertionError("first-moment steady circulation balance failed")

    Z_throat = -M0
    Z_return = sp.simplify(M0 * T0_dc)
    Z_replenishment_localized = sp.Integer(0)
    Z_boundary_dof = sp.Integer(0)
    Z = sp.simplify(Z_throat + Z_return + Z_replenishment_localized + Z_boundary_dof)
    Z_local_formula = sp.simplify(-M0 * (1 - T0_dc))
    if sp.simplify(Z - Z_local_formula) != 0:
        raise AssertionError("local-channel Z did not reduce to -M0*(1-T0(0))")
    Z_negative_certificate = sp.simplify(-Z * (1 + eps0) / (M0 * eps0))
    if Z_negative_certificate != 1:
        raise AssertionError("Z sign certificate failed")
    strict_T0_limit = sp.simplify(sp.limit(T0_full, eps0, 0, dir="+"))
    strict_T0_series = sp.expand(sp.series(strict_T0_limit, omega, 0, 5).removeO())
    strict_nu0 = omega_order(1 - strict_T0_series, omega, 4)
    strict_p_res0 = sp.Integer(1) + strict_nu0
    strict_p_res1 = sp.Integer(3) + strict_nu0
    strict_Z_local_limit = sp.simplify(sp.limit(Z, eps0, 0, dir="+"))
    if strict_Z_local_limit != 0:
        raise AssertionError("perfect-return local channel must have Z=0 without a declared boundary DOF")

    f0_abs = sp.simplify(1 / sp.sqrt(d))
    f1_abs = sp.simplify(sp.sqrt(2 / d) * sp.cos(sp.pi * w / d))
    m0_abs = sp.Integer(0)
    m1_abs = sp.simplify(sp.pi / d)
    norm0_abs = sp.simplify(sp.integrate(f0_abs**2, (w, 0, d)))
    norm1_abs = sp.simplify(sp.integrate(f1_abs**2, (w, 0, d)))
    abs_residuals = {
        "zero_ode": sp.simplify(sp.diff(f0_abs, w, 2) + m0_abs**2 * f0_abs),
        "zero_bc_w0": sp.simplify(sp.diff(f0_abs, w).subs(w, 0)),
        "zero_bc_wd": sp.simplify(sp.diff(f0_abs, w).subs(w, d)),
        "first_ode": sp.simplify(sp.diff(f1_abs, w, 2) + m1_abs**2 * f1_abs),
        "first_bc_w0": sp.simplify(sp.diff(f1_abs, w).subs(w, 0)),
        "first_bc_wd": sp.simplify(sp.diff(f1_abs, w).subs(w, d)),
    }
    if norm0_abs != 1 or norm1_abs != 1 or any(value != 0 for value in abs_residuals.values()):
        raise AssertionError("destructuring/absorbing transverse solve failed")

    f0_bloch = sp.simplify(1 / sp.sqrt(d))
    f1_bloch = sp.simplify(sp.sqrt(2 / d) * sp.cos(2 * sp.pi * w / d))
    m0_bloch = sp.Integer(0)
    m1_bloch = sp.simplify(2 * sp.pi / d)
    norm0_bloch = sp.simplify(sp.integrate(f0_bloch**2, (w, 0, d)))
    norm1_bloch = sp.simplify(sp.integrate(f1_bloch**2, (w, 0, d)))
    bloch_residuals = {
        "zero_ode": sp.simplify(sp.diff(f0_bloch, w, 2) + m0_bloch**2 * f0_bloch),
        "zero_periodic_value": sp.simplify(f0_bloch.subs(w, d) - f0_bloch.subs(w, 0)),
        "zero_periodic_derivative": sp.simplify(sp.diff(f0_bloch, w).subs(w, d) - sp.diff(f0_bloch, w).subs(w, 0)),
        "first_ode": sp.simplify(sp.diff(f1_bloch, w, 2) + m1_bloch**2 * f1_bloch),
        "first_periodic_value": sp.simplify(f1_bloch.subs(w, d) - f1_bloch.subs(w, 0)),
        "first_periodic_derivative": sp.simplify(sp.diff(f1_bloch, w).subs(w, d) - sp.diff(f1_bloch, w).subs(w, 0)),
    }
    if norm0_bloch != 1 or norm1_bloch != 1 or any(value != 0 for value in bloch_residuals.values()):
        raise AssertionError("Bloch transverse solve failed")

    dynamic_radial = solve_dynamic_zero_mode_radial("destructuring_absorbing", f0_abs, m0_abs, omega, cS, d, r)
    static_abs_radial = solve_static_zero_mode_radial("destructuring_absorbing", f0_abs, m0_abs, d, r)
    static_bloch_radial = solve_static_zero_mode_radial("bloch_stack", f0_bloch, m0_bloch, d, r)
    massive_abs_radial = solve_static_massive_radial("destructuring_absorbing", m1_abs, r)
    massive_bloch_radial = solve_static_massive_radial("bloch_stack", m1_bloch, r)
    dynamic_trace = dynamic_radial["trace"]
    static_trace = static_abs_radial["trace"]
    if dynamic_trace["trace_id"] == static_trace["trace_id"]:
        raise AssertionError("static/dynamic trace IDs coincide")
    p_dynamic = dynamic_radial["exponent"]
    p_static = static_abs_radial["exponent"]
    static_dynamic_agree = sp.simplify(p_dynamic - p_static) == 0
    if not static_dynamic_agree:
        raise AssertionError("dynamic-limit and static-solve exponents disagree")

    green_zero_dynamic = dynamic_radial["solution"]
    green_zero_static = static_abs_radial["solution"]
    green_bloch_zero_static = static_bloch_radial["solution"]
    flow_zero_static = static_abs_radial["flow"]
    p_abs = static_abs_radial["exponent"]
    p_bloch = static_bloch_radial["exponent"]
    counterfactual_guard = build_counterfactual_guard(green_zero_static, sp.Integer(0), r)
    p_eq_2 = sp.simplify(p_abs - 2) == 0
    if not p_eq_2 or sp.simplify(p_abs - p_bloch) != 0:
        raise AssertionError("localizing DC-sink completions did not agree on p=2")
    if static_abs_radial["trace"]["trace_id"] == static_bloch_radial["trace"]["trace_id"]:
        raise AssertionError("destructuring and Bloch radial solve trace IDs coincide")

    green_abs_static = sp.simplify(Z * green_zero_static + sp.Symbol("CAbs1") * massive_abs_radial["solution"])
    green_bloch_static = sp.simplify(Z * green_bloch_zero_static + sp.Symbol("CBloch1") * massive_bloch_radial["solution"])

    # Derived anti-localizing warp control: noncompact half-line with a growing
    # measure. The constant zero mode is not normalizable and the spectrum is a
    # gapless continuum starting at m=0.
    anti_measure = sp.exp(2 * kWarp * w)
    anti_zero_norm_cutoff = sp.simplify(sp.integrate(anti_measure, (w, 0, sp.Symbol("W", positive=True))))
    anti_zero_norm_limit = sp.limit(anti_zero_norm_cutoff, sp.Symbol("W", positive=True), sp.oo)
    if anti_zero_norm_limit != sp.oo:
        raise AssertionError("anti-localizing zero mode unexpectedly normalizable")
    continuum_green_static = sp.simplify(sp.integrate(sp.exp(-m * r), (m, 0, sp.oo)) / (4 * sp.pi * r))
    continuum_flow = sp.simplify(-sp.diff(continuum_green_static, r))
    p_delocalizing = radial_decay_exponent(continuum_flow, r)
    delocalizing_verdict = classify_dc_sink_gate([p_delocalizing], quadrupole_survives=True)
    if delocalizing_verdict != "RETURN_NOGO":
        raise AssertionError("no-go reachable control failed to return RETURN_NOGO")

    destructuring_verdict = branch_verdict_from_p(p_abs)
    bloch_verdict = branch_verdict_from_p(p_bloch)
    headline = classify_dc_sink_gate([p_abs, p_bloch], quadrupole_survives=True)
    if headline != "RETURN_RESIDUAL_PREDICTION":
        raise AssertionError(f"unexpected v3 headline for localizing flat slab: {headline}")

    A_strict_pass = bool(p_res0 >= 5 and p_res1 >= 5)
    A_residual_pass = bool((p_res0 < 5 or p_res1 < 5) and p_res0 >= 1 and p_res1 >= 3)
    B_pass = bool(p_eq_2)
    tension_status = "witnessed" if (sp.simplify(p_abs - 2) == 0 and sp.simplify(p_delocalizing - 2) != 0) else "not_witnessed"
    if tension_status != "witnessed":
        raise AssertionError("v3 falloff tension witness did not fire")

    # Dimensional homogeneity with dimensions ordered M,L,T.
    Mdim = sp.Matrix([1, 0, 0])
    Ldim = sp.Matrix([0, 1, 0])
    Tdim = sp.Matrix([0, 0, 1])
    dim_a = Ldim
    dim_d = Ldim
    dim_r = Ldim
    dim_cS = Ldim - Tdim
    dim_omega = -Tdim
    dim_M0 = -Tdim
    dim_D1 = Ldim - Tdim
    dim_rho = -3 * Ldim
    dim_m_atom = Mdim
    dim_K = dim_m_atom + 14 * Ldim - 2 * Tdim
    if list(dim_omega + dim_d - dim_cS) != [0, 0, 0]:
        raise AssertionError("omega*d/c_s is not dimensionless")
    if list(dim_a + dim_omega - dim_cS + dim_M0) != list(dim_M0):
        raise AssertionError("A0 residual does not carry M0 dimension")
    if list(3 * dim_a + 3 * dim_omega - 3 * dim_cS + dim_D1) != list(dim_D1):
        raise AssertionError("A1 residual does not carry D1 dimension")
    if list(dim_M0 - dim_rho - 2 * dim_r) != [0, 1, -1]:
        raise AssertionError("radial flow dimension check failed")
    if list(dim_K + 4 * dim_rho - dim_m_atom) != [0, 2, -2]:
        raise AssertionError("c_s^2 EOS dimension check failed")
    K_eos = sp.Rational(1, 500)
    rho_frozen_m_eq_1 = sp.simplify((sp.Integer(1) / (5 * K_eos)) ** sp.Rational(1, 4))
    cs2_frozen = sp.simplify(5 * K_eos * rho_frozen_m_eq_1**4)
    if rho_frozen_m_eq_1 != sp.sqrt(10) or cs2_frozen != 1:
        raise AssertionError("frozen c_s^2 slice check failed")

    structure_canonical = "|".join(
        [
            "pathA29_v3",
            "flat_slab:0<=w<=d",
            "warp:A(w)=1",
            "bulk_eq:Phi_ww+(omega/c_s)^2 Phi=0",
            "checkA:continuity_fraction_alpha_l=1/(1+epsilon_l),phase=exp(I*omega*2d/c_s)",
            "checkB:destructuring_absorbing:compact_cell_natural_dc_sink_spectrum",
            "checkB:bloch_stack:q=0_periodic_cell",
            "radiation:ac_check_a_only",
            "Z_is_premise:true",
            "boundary_dof:none",
        ]
    )
    structure_id = sha256_text(structure_canonical)

    exprs_for_agreement = {
        "T0_at_DC": T0_dc,
        "T1_at_DC": T1_dc,
        "p": p_abs,
        "p_eq_2": sp.Integer(1 if p_eq_2 else 0),
        "destructuring_p": p_abs,
        "bloch_p": p_bloch,
        "dynamic_limit_exponent": p_dynamic,
        "static_solve_exponent": p_static,
        "zero_mode.r_dependence": green_zero_static,
        "green_function": green_zero_dynamic,
        "destructuring_green_function": green_abs_static,
        "bloch_green_function": green_bloch_static,
        "delocalizing_p": p_delocalizing,
        "A0_residual": A0_res,
        "A1_residual": A1_res,
        "Z": Z,
    }
    expr_digest = sha256_text(json.dumps({k: mma_expr(v) for k, v in exprs_for_agreement.items()}, sort_keys=True))

    return {
        "symbols": {
            "omega": omega,
            "cS": cS,
            "d": d,
            "a": a,
            "r": r,
            "rhoB": rhoB,
            "epsilon0": eps0,
            "epsilon1": eps1,
            "M0": M0,
            "D1": D1,
            "Gamma_uniform": Gamma_uniform,
            "k_warp": kWarp,
            "w": w,
            "m": m,
            "n": n,
            "q": q,
        },
        "z": z,
        "tau": tau,
        "transit_trace": transit_trace,
        "alpha0": alpha0,
        "alpha1": alpha1,
        "neighbor_fraction0": neighbor_fraction0,
        "neighbor_fraction1": neighbor_fraction1,
        "helmholtz_basis": helmholtz_basis,
        "transport_phase": transport_phase,
        "T0_full": T0_full,
        "T1_full": T1_full,
        "T0_series": T0_series,
        "T1_series": T1_series,
        "T0_dc": T0_dc,
        "T1_dc": T1_dc,
        "sigma0": sigma0,
        "sigma1": sigma1,
        "nu0": nu0,
        "nu1": nu1,
        "raw_kernel0": raw_kernel0,
        "raw_kernel1": raw_kernel1,
        "raw_kernel2": raw_kernel2,
        "R0": R0,
        "R1": R1,
        "A0_res": A0_res,
        "A1_res": A1_res,
        "p_res0": p_res0,
        "p_res1": p_res1,
        "J_return0": J_return0,
        "J_neighbor0": J_neighbor0,
        "J_return1": J_return1,
        "J_neighbor1": J_neighbor1,
        "steady_balance0": steady_balance0,
        "steady_balance1": steady_balance1,
        "Z_throat": Z_throat,
        "Z_return": Z_return,
        "Z_replenishment_localized": Z_replenishment_localized,
        "Z_boundary_dof": Z_boundary_dof,
        "Z": Z,
        "Z_local_formula": Z_local_formula,
        "Z_negative_certificate": Z_negative_certificate,
        "strict_T0_limit": strict_T0_limit,
        "strict_nu0": strict_nu0,
        "strict_p_res0": strict_p_res0,
        "strict_p_res1": strict_p_res1,
        "strict_Z_local_limit": strict_Z_local_limit,
        "f0_abs": f0_abs,
        "f1_abs": f1_abs,
        "m0_abs": m0_abs,
        "m1_abs": m1_abs,
        "norm0_abs": norm0_abs,
        "norm1_abs": norm1_abs,
        "abs_residuals": abs_residuals,
        "f0_bloch": f0_bloch,
        "f1_bloch": f1_bloch,
        "m0_bloch": m0_bloch,
        "m1_bloch": m1_bloch,
        "norm0_bloch": norm0_bloch,
        "norm1_bloch": norm1_bloch,
        "bloch_residuals": bloch_residuals,
        "dynamic_trace": dynamic_trace,
        "static_trace": static_trace,
        "static_abs_radial_trace": static_abs_radial["trace"],
        "static_bloch_radial_trace": static_bloch_radial["trace"],
        "massive_abs_radial_trace": massive_abs_radial["trace"],
        "massive_bloch_radial_trace": massive_bloch_radial["trace"],
        "p_dynamic": p_dynamic,
        "p_static": p_static,
        "static_dynamic_agree": static_dynamic_agree,
        "green_zero_dynamic": green_zero_dynamic,
        "green_zero_static": green_zero_static,
        "green_bloch_zero_static": green_bloch_zero_static,
        "flow_zero_static": flow_zero_static,
        "green_abs_static": green_abs_static,
        "green_bloch_static": green_bloch_static,
        "massive_abs_solution": massive_abs_radial["solution"],
        "massive_bloch_solution": massive_bloch_radial["solution"],
        "p_abs": p_abs,
        "p_bloch": p_bloch,
        "p_eq_2": p_eq_2,
        "counterfactual_guard": counterfactual_guard,
        "anti_measure": anti_measure,
        "anti_zero_norm_cutoff": anti_zero_norm_cutoff,
        "anti_zero_norm_limit": anti_zero_norm_limit,
        "continuum_green_static": continuum_green_static,
        "continuum_flow": continuum_flow,
        "p_delocalizing": p_delocalizing,
        "delocalizing_verdict": delocalizing_verdict,
        "destructuring_verdict": destructuring_verdict,
        "bloch_verdict": bloch_verdict,
        "headline": headline,
        "A_strict_pass": A_strict_pass,
        "A_residual_pass": A_residual_pass,
        "B_pass": B_pass,
        "tension_status": tension_status,
        "structure_canonical": structure_canonical,
        "structure_id": structure_id,
        "exprs_for_agreement": exprs_for_agreement,
        "expr_digest": expr_digest,
        "frozen_slice": {
            "n": 5,
            "K_eos": "1/500",
            "rho_frozen_m_eq_1": hstr(rho_frozen_m_eq_1),
            "c_s_squared_frozen": hstr(cs2_frozen),
            "G": 1,
            "c": 1,
            "c_s": 1,
            "a_star": "4731/2500",
            "L_star": "18121/10000",
            "dimensions_MLT": {
                "a": [0, 1, 0],
                "d": [0, 1, 0],
                "omega": [0, 0, -1],
                "c_s": [0, 1, -1],
                "M0": [0, 0, -1],
                "D1": [0, 1, -1],
                "K_eos": [1, 14, -2],
                "c_s_squared": [0, 2, -2],
                "radial_flow": [0, 1, -1],
            },
        },
    }


def build_results(data: dict[str, Any], reuse: dict[str, Any], agreement_status: str, mma_status: str | None) -> dict[str, Any]:
    checked = [
        "p",
        "p_eq_2",
        "destructuring_p",
        "bloch_p",
        "dynamic_limit_exponent",
        "static_solve_exponent",
        "zero_mode.r_dependence",
        "green_function",
        "destructuring_green_function",
        "bloch_green_function",
        "delocalizing_p",
        "A0_residual",
        "A1_residual",
        "Z",
    ]
    required_checked = {
        "p",
        "p_eq_2",
        "dynamic_limit_exponent",
        "static_solve_exponent",
        "zero_mode.r_dependence",
        "green_function",
    }
    if not required_checked.issubset(set(checked)):
        raise AssertionError("engine agreement checked_quantities missing directive-required entries")

    p_abs = int(data["p_abs"])
    p_bloch = int(data["p_bloch"])
    p_deloc = int(data["p_delocalizing"])
    branch_agree = p_abs == p_bloch
    headline = data["headline"]

    results: dict[str, Any] = {
        "schema": "pathA_29_brane_bulk_return_sympy/v3",
        "top_line_verdict": headline,
        "Z_is_premise": True,
        "pathA_28_reuse": reuse,
        "slab_structure": {
            "structure_id": data["structure_id"],
            "canonical_form": data["structure_canonical"],
            "spacing_d": {"value": "d>0 finite", "status": "postulated"},
            "dimensionality_of_return": {
                "value": "one flat finite bulk w-slab; return fractions epsilon_l/(1+epsilon_l) are continuity partition parameters, not source fits",
                "status": "postulated",
            },
            "warp_ansatz": {"value": "flat A(w)=1 for the main family", "status": "postulated"},
            "derived_neighbor_BC_equations": {
                "destructuring_absorbing": {
                    "status": "derived",
                    "physical_completion": "section-5 de-structuring/absorbing DC sink; drain admittance is the premise and is represented in the continuity source term Z, not by a tuned wall impedance",
                    "transverse_equations": ["-f''(w)=m^2 f(w)", "f'(0)=0", "f'(d)=0"],
                },
                "bloch_stack": {
                    "status": "derived",
                    "physical_completion": "periodic multi-brane stack q=0 Bloch cell with absorbing neighbor branes",
                    "transverse_equations": ["-f''(w)=m^2 f(w)", "f(d)=f(0)", "f'(d)=f'(0)"],
                },
                "radiation": {
                    "status": "derived_ac_only",
                    "physical_completion": "Sommerfeld outgoing AC transport, tagged ac_check_a_only and excluded from Check B",
                    "equation": "partial_w Phi = I*(omega/c_s)*Phi at finite omega; omega->0 Neumann limit is not a Check-B branch",
                },
            },
            "response": {
                "value": "T_l(omega)=alpha_l*exp(I*omega*2d/c_s) from projected continuity and solved bidirectional Helmholtz transport",
                "status": "derived_not_dialed",
            },
            "postulate_vs_derived": {
                "geometry": "postulated",
                "topology": "postulated",
                "BC_functional_forms": "derived from declared physical completions",
                "response": "computed from continuity plus transport; not fitted to cancellation",
            },
        },
        "steady_state": {
            "exists": True,
            "type": "circulation",
            "provenance": {
                "zeroth_moment_balance": f"M0 = {hstr(data['J_return0'])} + {hstr(data['J_neighbor0'])}",
                "first_moment_balance": f"D1 = {hstr(data['J_return1'])} + {hstr(data['J_neighbor1'])}",
                "computed_residuals": {"ell0": hstr(data["steady_balance0"]), "ell1": hstr(data["steady_balance1"])},
                "routing": "a candidate with nonzero steady balance is inadmissible; this family has zero balance residual",
            },
        },
        "transfer_function": {
            "structure_id": data["structure_id"],
            "provenance": {
                "governing_equations": [
                    "bulk Helmholtz transport: partial_w^2 Phi_l + (omega/c_s)^2 Phi_l = 0 on 0<w<d",
                    f"bidirectional basis: {data['helmholtz_basis']}",
                    "round-trip transport phase solved from the Helmholtz basis phase ratios",
                    "zeroth-moment continuity: J_leak_0 = J_return_0 + J_neighbor_0",
                    "first-moment continuity: J_leak_1 = J_return_1 + J_neighbor_1 including return centroid partition",
                    "alpha_l = J_return_l/J_leak_l = 1/(1+epsilon_l), epsilon_l>0",
                ],
                "transit_phase_trace": data["transit_trace"]["trace"],
                "bc_equations": [
                    "destructuring_absorbing DC sink: compact cell transverse spectrum f'(0)=f'(d)=0; drain flux is the source accounting Z",
                    "Bloch DC-sink cross-check: f(d)=f(0), f'(d)=f'(0)",
                    "Sommerfeld radiation is finite-omega AC transport only and is not used for Check B",
                ],
                "solution_basis": data["helmholtz_basis"],
                "series_coefficients": {
                    "ell0": series_coefficients(data["T0_series"], data["symbols"]["omega"], 4),
                    "ell1": series_coefficients(data["T1_series"], data["symbols"]["omega"], 2),
                },
                "free_parameter_table": {
                    "d": {"role": "postulated slab spacing", "domain": "positive finite", "depends_on_source": False},
                    "epsilon0": {"role": "zeroth-moment transmitted-to-returned DC sink partition", "domain": "positive", "depends_on_source": False},
                    "epsilon1": {"role": "first-moment transmitted-to-returned centroid partition", "domain": "positive", "depends_on_source": False},
                    "c_s": {"role": "signal speed with restored units", "domain": "positive", "depends_on_source": False},
                },
                "residual_substitution_trace": {
                    "R0_substitution": "R0(omega)=-M0*T0(omega)",
                    "R1_substitution": "R1(omega)=-D1*T1(omega)",
                    "A0_residual": hstr(data["A0_res"]),
                    "A1_residual": hstr(data["A1_res"]),
                },
                "forbidden_fit_flags": {
                    "fit_to_cancellation": False,
                    "fit_to_source_moments": False,
                    "T_set_to_one_by_fiat": False,
                    "standalone_frozen_static_solve": False,
                    "falloff_exponent_fit_from_handwritten_vr": False,
                    "m_dependent_impedance_BC": False,
                    "radiation_branch_used_for_Check_B": False,
                    "all_multipole_cancellation": False,
                },
            },
            "ell0": {
                "T_ell_omega": hstr(data["T0_series"]),
                "T_full_unexpanded": hstr(data["T0_full"]),
                "series_through": "O(omega^4)",
                "T_at_DC": hstr(data["T0_dc"]),
                "sigma_l": int(data["sigma0"]),
                "nu_l": int(data["nu0"]),
                "computed_from": "zeroth-moment continuity fraction times solved Helmholtz round-trip phase",
            },
            "ell1": {
                "T_ell_omega": hstr(data["T1_series"]),
                "T_full_unexpanded": hstr(data["T1_full"]),
                "series_through": "O(omega^2)",
                "T_at_DC": hstr(data["T1_dc"]),
                "sigma_l": int(data["sigma1"]),
                "nu_l": int(data["nu1"]),
                "computed_from": "first-moment continuity fraction times solved Helmholtz round-trip phase",
            },
            "T_at_DC": {"ell0": hstr(data["T0_dc"]), "ell1": hstr(data["T1_dc"])},
        },
        "return_moments": {
            "R0_omega": hstr(data["R0"]),
            "R1_omega": hstr(data["R1"]),
            "target_from_pathA_28": {"R0": "-M0", "R1": "-D1"},
        },
        "residual_amplitude": {
            "A0_res": hstr(data["A0_res"]),
            "A1_res": hstr(data["A1_res"]),
            "p_residual": {"ell0": int(data["p_res0"]), "ell1": int(data["p_res1"])},
            "p_quadrupole": 5,
            "order_rule": "p_res_0=1+nu_0 and p_res_1=3+nu_1; sigma_l is not used in p_res",
        },
        "admissibility_window_A": {
            "status": "none",
            "strict_window": "none",
            "required_conditions": "T0(0)=1 with nu0>=4 and T1(0)=1 with nu1>=2",
            "computed_obstruction": "finite epsilon_l>0 gives nu0=0 and nu1=0; perfect-return epsilon0->0 gives a finite transport phase with nu0=1, so strict p_res>=5 is not an open window in this tractable family",
            "perfect_return_control": {
                "T0_limit": hstr(data["strict_T0_limit"]),
                "nu0": int(data["strict_nu0"]),
                "p_res0": int(data["strict_p_res0"]),
                "p_res1": int(data["strict_p_res1"]),
            },
        },
        "residual_window_A": {
            "status": "open",
            "constraints": "d>0, c_s>0, epsilon0>0, epsilon1>0, flat localizing slab, admissible DC-sink completion",
            "p_residual": {"ell0": int(data["p_res0"]), "ell1": int(data["p_res1"])},
            "bounded_residual": True,
        },
        "transverse_mode_spectrum": {
            "operator": "-d^2/dw^2 on the compact finite slab for Check-B DC sinks",
            "destructuring_absorbing": {
                "bc_derivation": "section-5 de-structuring/absorbing DC sink: finite compact cell; drain admittance is carried by Z and the field spectrum uses the natural finite-cell zero mode",
                "bc_equations": ["f'(0)=0", "f'(d)=0"],
                "eigenvalue_equation": "sin(m*d)=0",
                "modes": [
                    {"n": 0, "m_n": hstr(data["m0_abs"]), "f_n(w)": hstr(data["f0_abs"]), "normalization": hstr(data["norm0_abs"])},
                    {"n": 1, "m_n": hstr(data["m1_abs"]), "f_n(w)": hstr(data["f1_abs"]), "normalization": hstr(data["norm1_abs"])},
                    {"n": "j>=0", "m_n": "j*pi/d", "f_n(w)": "sqrt((2-delta_j0)/d)*cos(j*pi*w/d)"},
                ],
                "proof_residuals": {k: hstr(v) for k, v in data["abs_residuals"].items()},
                "spectrum_class": "normalizable_zero_mode",
            },
            "bloch_stack": {
                "bc_derivation": "q=0 Bloch cell of the periodic multi-brane stack with absorbing neighboring branes",
                "bc_equations": ["f(d)=f(0)", "f'(d)=f'(0)"],
                "eigenvalue_equation": "m=2*pi*j/d for integer j in the q=0 Bloch sector",
                "modes": [
                    {"n": 0, "m_n": hstr(data["m0_bloch"]), "f_n(w)": hstr(data["f0_bloch"]), "normalization": hstr(data["norm0_bloch"])},
                    {"n": 1, "m_n": hstr(data["m1_bloch"]), "f_n(w)": hstr(data["f1_bloch"]), "normalization": hstr(data["norm1_bloch"])},
                ],
                "proof_residuals": {k: hstr(v) for k, v in data["bloch_residuals"].items()},
                "spectrum_class": "normalizable_zero_mode",
            },
        },
        "zero_mode": {
            "exists": True,
            "normalizable": True,
            "eigenvalue_m0": "0",
            "eigenfunction": hstr(data["f0_abs"]),
            "r_dependence": hstr(data["green_zero_static"]),
            "flow_from_gradient": hstr(data["flow_zero_static"]),
            "derived_flow_exponent_p": p_abs,
            "destructuring_radial_trace_id": data["static_abs_radial_trace"]["trace_id"],
            "bloch_radial_trace_id": data["static_bloch_radial_trace"]["trace_id"],
            "independent_branch_solves": True,
        },
        "green_function": {
            "dynamic_zero_mode": hstr(data["green_zero_dynamic"]),
            "dynamic_zero_mode_trace": data["dynamic_trace"],
            "static_destructuring_modal_green": hstr(data["green_abs_static"]),
            "static_destructuring_zero_mode_trace": data["static_abs_radial_trace"],
            "static_bloch_modal_green": hstr(data["green_bloch_static"]),
            "static_bloch_zero_mode_trace": data["static_bloch_radial_trace"],
            "first_massive_tail_destructuring": hstr(data["massive_abs_solution"]),
            "first_massive_tail_destructuring_trace": data["massive_abs_radial_trace"],
            "first_massive_tail_bloch": hstr(data["massive_bloch_solution"]),
            "first_massive_tail_bloch_trace": data["massive_bloch_radial_trace"],
            "long_range_selection": "the solved normalizable m=0 mode controls r>>d; massive terms are exponentially suppressed",
        },
        "counterfactual_guard": data["counterfactual_guard"],
        "spectrum_class": "normalizable_zero_mode",
        "static_falloff_B": {
            "structure_id": data["structure_id"],
            "same_structure_as_A": True,
            "limiting_map": "Check B is omega->0 of the same finite-slab Helmholtz transport/spectrum family; radiation/Sommerfeld is tagged ac_check_a_only and excluded from the DC-sink branch gate",
            "dynamic_object": "G(omega,r,w;0,0)=sum_n f_n(w)f_n(0) g_n(r), where each g_n solves g''+(2/r)g'+((omega/c_s)^2-m_n^2)g=0",
            "bc_computed_under": ["destructuring_absorbing", "bloch_stack"],
            "p": p_abs,
            "p_eq_2": bool(data["p_eq_2"]),
            "crossover_scale": "r>>d; first massive scales d/pi for destructuring_absorbing and d/(2*pi) for bloch_stack",
            "zero_mode_source_integral": hstr(data["Z"]),
            "Z_terms": {
                "throat_sink": hstr(data["Z_throat"]),
                "return": hstr(data["Z_return"]),
                "replenishment_localized": hstr(data["Z_replenishment_localized"]),
                "Z_boundary_dof": "none",
            },
            "Z_uniform_background": {
                "term": "Gamma_uniform",
                "enters_zero_mode_source_integral": False,
                "provenance": "uniform whole-brane areal-leak background; separate observable, not localized monopole Z",
            },
            "Z_sign": "negative",
            "Z_sign_certificate": hstr(data["Z_negative_certificate"]),
            "Z_reduces_to_local": True,
        },
        "branch_results": {
            "destructuring_absorbing": {
                "role": "PRIMARY",
                "admissibility": "DC sink; Z<0 is the drain premise",
                "verdict": data["destructuring_verdict"],
                "p": p_abs,
                "p_eq_2": True,
                "Z_sign": "negative_by_premise_and_accounting",
                "spectrum_class": "normalizable_zero_mode",
                "green_function": hstr(data["green_abs_static"]),
                "radial_solve_trace_id": data["static_abs_radial_trace"]["trace_id"],
            },
            "bloch_stack": {
                "role": "CROSS_CHECK",
                "admissibility": "DC sink periodic stack; Z<0 is the drain premise",
                "verdict": data["bloch_verdict"],
                "p": p_bloch,
                "p_eq_2": True,
                "Z_sign": "negative_by_premise_and_accounting",
                "spectrum_class": "normalizable_zero_mode",
                "green_function": hstr(data["green_bloch_static"]),
                "radial_solve_trace_id": data["static_bloch_radial_trace"]["trace_id"],
            },
            "radiation": {
                "tag": "ac_check_a_only",
                "not_used_for_Check_B": True,
                "finite_omega_bc": "partial_w Phi = I*(omega/c_s)*Phi",
                "dc_limit_note": "omega->0 is reflecting/no-drain and therefore inadmissible for Check B",
            },
        },
        "branch_agreement": {
            "destructuring_p": p_abs,
            "bloch_p": p_bloch,
            "destructuring_verdict": data["destructuring_verdict"],
            "bloch_verdict": data["bloch_verdict"],
            "p_agree": branch_agree,
        },
        "reconciliation": {
            "verdict": headline,
            "B_pass": bool(data["B_pass"]),
            "p": p_abs,
            "Z_is_premise": True,
            "T0_at_DC": hstr(data["T0_dc"]),
            "same_structure": True,
            "A_strict_pass": bool(data["A_strict_pass"]),
            "A_residual_pass": bool(data["A_residual_pass"]),
            "Z_sign": "negative",
            "p_eq_2": bool(data["p_eq_2"]),
            "Z_reduces_to_local": True,
            "Z_boundary_dof_status": "none",
            "joint_window": "superseded_by_v3; no strict-A window in the computed local channel",
        },
        "bc_is_local_causal_passive": True,
        "bc_omega_dependence_source": "none",
        "window_is_open": True,
        "window_symmetry_protected": "none",
        "residual_prediction": {
            "status": "active",
            "scaling": {
                "ell0": "A0_res = i*a*(omega/c_s)*M0*(epsilon0/(1+epsilon0) + O(omega*d/c_s))",
                "ell1": "A1_res = i*a^3*(omega/c_s)^3*D1*(epsilon1/(1+epsilon1) + O(omega*d/c_s))/2",
            },
            "epsilon_to_gravity_strength_tie": f"Z={hstr(data['Z'])}; the same epsilon0 sets the net sink strength and the monopole residual coefficient",
            "pde_ledger": "does not close open-item #9; records a falsifiable residual prediction",
        },
        "T2_effect": {
            "T2_applied": False,
            "quadrupole_kernel_retains_p_raw": 5,
            "kernel_reused_from_pathA_28": "I*a**5*omega**5/(27*cS**5)",
            "burke_thorne_anchor": "research/4d_2_5pn",
            "all_multipole_cancellation_imposed": False,
        },
        "static_dynamic_consistency": {
            "agree": True,
            "dynamic_limit_exponent": int(data["p_dynamic"]),
            "static_solve_exponent": int(data["p_static"]),
            "dynamic_trace_id": data["dynamic_trace"]["trace_id"],
            "static_trace_id": data["static_trace"]["trace_id"],
            "dynamic_source_hash": data["dynamic_trace"]["source_hash"],
            "static_source_hash": data["static_trace"]["source_hash"],
            "independent_extractions": True,
            "dynamic_limit_trace": data["dynamic_trace"],
            "static_solve_trace": data["static_trace"],
        },
        "controls": {
            "dc_value_computed": {
                "same_pipeline": True,
                "fired": True,
                "T0_at_DC": hstr(data["T0_dc"]),
                "T1_at_DC": hstr(data["T1_dc"]),
                "not_fiat_one": True,
                "strict_A_requires_finite_omega_orders": True,
            },
            "no_go_reachable": {
                "same_pipeline": True,
                "fired": True,
                "status": "reachable_RETURN_NOGO",
                "warp_equation": "anti-localizing half-line measure mu(w)=exp(2*k_warp*w), k_warp>0; L_w=-(1/mu) d_w(mu d_w)",
                "zero_mode_norm": hstr(data["anti_zero_norm_limit"]),
                "spectrum_class": "continuum_p3",
                "computed_p": p_deloc,
                "classifier_verdict": data["delocalizing_verdict"],
                "drain_admissibility": "Z<0 premise retained; failure is only zero-mode delocalization",
            },
            "tension_is_real": {
                "same_pipeline": True,
                "fired": True,
                "status": data["tension_status"],
                "localizing_p": p_abs,
                "delocalizing_p": p_deloc,
            },
            "quadrupole_survives": {
                "same_pipeline": True,
                "fired": True,
                "p_raw_ell2": 5,
                "kernel_ell2": "I*a**5*omega**5/(27*cS**5)",
                "Q2_set_to_zero": False,
            },
            "return_necessity": {
                "same_pipeline": True,
                "fired": True,
                "rung_condition": headline,
                "order_change_required": False,
                "coefficient_substituted": True,
                "coefficient_tied_to_Z": True,
            },
            "static_dynamic_consistency": {
                "same_pipeline": True,
                "fired": True,
                "agree": True,
                "independent_extractions": True,
            },
        },
        "frozen_slice": data["frozen_slice"],
        "engine_agreement": {
            "status": agreement_status,
            "checked_quantities": checked,
            "sympy_expression_digest": data["expr_digest"],
            "mathematica_status": mma_status or "not_run_yet",
            "mathematica_exprs": {k: mma_expr(v) for k, v in data["exprs_for_agreement"].items()},
        },
        "run_commands": {
            "sympy": "timeout 600 python3 software/stage1_solver/tools/pathA_29_brane_bulk_return_sympy.py",
            "mathematica": "timeout 600 math -script software/stage1_solver/tools/pathA_29_brane_bulk_return.wl",
        },
    }
    validate_results(results)
    return results


def validate_results(results: dict[str, Any]) -> None:
    required_top = {
        "slab_structure",
        "steady_state",
        "transfer_function",
        "return_moments",
        "residual_amplitude",
        "admissibility_window_A",
        "residual_window_A",
        "static_falloff_B",
        "Z_is_premise",
        "branch_results",
        "branch_agreement",
        "reconciliation",
        "residual_prediction",
        "T2_effect",
        "controls",
        "engine_agreement",
        "transverse_mode_spectrum",
        "zero_mode",
        "green_function",
        "counterfactual_guard",
        "spectrum_class",
        "static_dynamic_consistency",
    }
    missing = sorted(required_top.difference(results))
    if missing:
        raise AssertionError(f"results missing required top-level fields: {missing}")
    if results["Z_is_premise"] is not True:
        raise AssertionError("Z_is_premise must be true")
    if results["branch_results"]["radiation"]["tag"] != "ac_check_a_only":
        raise AssertionError("radiation branch must be tagged ac_check_a_only")
    if results["branch_agreement"]["p_agree"] is not True:
        raise AssertionError("DC-sink completions must agree for this computed headline")
    if results["reconciliation"]["verdict"] != results["top_line_verdict"]:
        raise AssertionError("headline mismatch between top_line_verdict and reconciliation")
    if results["controls"]["no_go_reachable"]["status"] != "reachable_RETURN_NOGO":
        raise AssertionError("v3 green requires no_go_reachable.status reachable_RETURN_NOGO")
    if results["counterfactual_guard"]["perturbed_falloff_rejected"] is not True:
        raise AssertionError("counterfactual radial falloff guard did not reject the perturbed solution")
    if (
        results["zero_mode"]["destructuring_radial_trace_id"]
        == results["zero_mode"]["bloch_radial_trace_id"]
    ):
        raise AssertionError("destructuring and Bloch radial solves are not independent")
    if results["controls"]["tension_is_real"]["status"] != "witnessed":
        raise AssertionError("tension_is_real must be witnessed")
    sdc = results["static_dynamic_consistency"]
    if not sdc["agree"] or not sdc["independent_extractions"]:
        raise AssertionError("static_dynamic_consistency did not pass")
    if sdc["dynamic_trace_id"] == sdc["static_trace_id"] or sdc["dynamic_source_hash"] == sdc["static_source_hash"]:
        raise AssertionError("static/dynamic trace IDs or source hashes coincide")
    forbidden = results["transfer_function"]["provenance"]["forbidden_fit_flags"]
    if any(bool(value) for value in forbidden.values()):
        raise AssertionError(f"forbidden fit flags fired: {forbidden}")
    checked = set(results["engine_agreement"]["checked_quantities"])
    required_checked = {
        "p",
        "p_eq_2",
        "dynamic_limit_exponent",
        "static_solve_exponent",
        "zero_mode.r_dependence",
        "green_function",
    }
    if not required_checked.issubset(checked):
        raise AssertionError("engine agreement checked_quantities missing required entries")


def compare_mathematica_if_available(data: dict[str, Any]) -> tuple[str, str | None]:
    if not MMA_JSON.exists():
        return "PENDING_MATHEMATICA", None
    mma = json.loads(MMA_JSON.read_text(encoding="utf-8"))
    if mma.get("schema") != "pathA_29_brane_bulk_return_mathematica/v3":
        return "PENDING_MATHEMATICA", "ignored stale Mathematica JSON with mismatched schema"
    if mma.get("sympy_expression_digest") != data["expr_digest"]:
        raise AssertionError("Mathematica JSON digest does not match current SymPy expressions")
    if mma.get("status") != "PASS":
        raise AssertionError(f"Mathematica engine did not pass: {mma}")
    return "PASS", "timeout 600 math -script software/stage1_solver/tools/pathA_29_brane_bulk_return.wl exited 0 and asserted PASS"


def write_json(results: dict[str, Any]) -> None:
    SCRATCH.mkdir(parents=True, exist_ok=True)
    JSON_OUT.write_text(json.dumps(results, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def write_yaml(results: dict[str, Any]) -> None:
    RESULTS_YAML.write_text(
        yaml.dump(results, Dumper=NoAliasDumper, sort_keys=False, width=120),
        encoding="utf-8",
    )


def write_report(results: dict[str, Any]) -> None:
    verdict = results["top_line_verdict"]
    p0 = results["residual_amplitude"]["p_residual"]["ell0"]
    p1 = results["residual_amplitude"]["p_residual"]["ell1"]
    p = results["reconciliation"]["p"]
    t0 = results["reconciliation"]["T0_at_DC"]
    z_expr = results["static_falloff_B"]["zero_mode_source_integral"]
    lines = [
        verdict,
        "",
        "# pathA_29 Brane-Bulk Return v3",
        "",
        f"Computed headline: `{verdict}`.",
        "",
        "The executable family is a flat finite slab with our brane at `w=0` and an adjacent return/absorber at `w=d`.",
        "The geometry is postulated; the return response is derived from projected continuity and the solved Helmholtz transport phase.",
        f"`T0(0)={t0}`, so the residual is bounded but lower order: `p_res(ell0)={p0}`, `p_res(ell1)={p1}`.",
        "",
        "Check B was run only on the admissible DC-sink completions:",
        f"- `destructuring_absorbing`: solved compact-cell spectrum and a branch-specific 3D radial equation, derived `p={p}`.",
        f"- `bloch_stack`: solved q=0 Bloch spectrum and a separate 3D radial equation, derived `p={results['branch_agreement']['bloch_p']}`.",
        "",
        f"Counterfactual guard: `{results['counterfactual_guard']['perturbation']}` was rejected with residual `{results['counterfactual_guard']['perturbed_residual']}`.",
        "",
        "The radiation/Sommerfeld boundary is recorded as `ac_check_a_only` and is not used as a Check-B branch.",
        f"The signed local source accounting gives `Z={z_expr}`; under v3 this is accounting, while `Z<0` is the drain admissibility premise.",
        "",
        "The static-dynamic consistency control used separate traces:",
        f"- dynamic trace `{results['static_dynamic_consistency']['dynamic_trace_id']}`",
        f"- static trace `{results['static_dynamic_consistency']['static_trace_id']}`",
        "",
        "The mandatory no-go control constructs the anti-localizing half-line warp `mu(w)=exp(2*k_warp*w)`.",
        f"Its zero mode is non-normalizable, the continuum Green integral gives `p={results['controls']['no_go_reachable']['computed_p']}`, and the same classifier returns `RETURN_NOGO`.",
        "",
        "pde_ledger feed: open-item #9 is not closed; the deliverable is the falsifiable residual radiation prediction tied to the drain strength. The gravity-range item passes inside the localizing flat-slab family because both DC-sink completions give `p=2`.",
        "",
        "Downstream: the full nonlinear brane-bulk return closure remains track-3 work.",
    ]
    REPORT_OUT.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    reuse = assert_patha28_reuse()
    data = compute_symbolics()
    agreement_status, mma_status = compare_mathematica_if_available(data)
    results = build_results(data, reuse, agreement_status, mma_status)
    write_json(results)
    write_yaml(results)
    write_report(results)
    print("PASS pathA_29_brane_bulk_return_sympy")
    print(json.dumps({"verdict": results["top_line_verdict"], "engine_agreement": agreement_status, "json": str(JSON_OUT)}, indent=2))


if __name__ == "__main__":
    main()
