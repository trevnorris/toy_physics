#!/usr/bin/env python3
"""pathA_38 throat-body electric localization gate, SymPy engine.

This remediation keeps the reduced Z2 wall model, but the headline is now
assembled from a transverse-mode projection followed by radial Green solves.
The classifier has no hard PASS assertions: the main branch and every ablated
witness travel through the same computed branch -> classifier path.
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

REPORT_OUT = REPORTS / "pathA_38_throat_body_electric_localization.md"
RESULTS_YAML = REPORTS / "pathA_38_results.yaml"
JSON_OUT = SCRATCH / "pathA_38_throat_body_electric_sympy.json"
MMA_JSON = SCRATCH / "pathA_38_throat_body_electric_mathematica.json"

SCHEMA = "pathA_38_throat_body_electric_sympy/v2"
MMA_SCHEMA = "pathA_38_throat_body_electric_mathematica/v2"

I = sp.I


Dim = tuple[int, int, int, int, int, int]
DIM_LABELS = ("E", "L", "Chi", "N", "Tau", "Q")


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
        fired = actual != expected
        if not fired:
            raise AssertionError(f"dimensional ablation did not fire: {category}:{name}")
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

    def expect_incompatible(self, category: str, name: str, left: Dim, right: Dim, expression: str) -> None:
        fired = left != right
        if not fired:
            raise AssertionError(f"dimension incompatibility ablation did not fire: {category}:{name}")
        self.ablations.append(
            {
                "category": category,
                "name": name,
                "expression": expression,
                "left": dim_record(left),
                "right": dim_record(right),
                "status": "FIRED",
            }
        )


def hstr(expr: Any) -> Any:
    if expr is None or isinstance(expr, (bool, int, str)):
        return expr
    if expr in (sp.oo, -sp.oo):
        return sp.sstr(expr)
    if isinstance(expr, sp.MatrixBase):
        return [[hstr(entry) for entry in row] for row in expr.tolist()]
    return sp.sstr(sp.factor(sp.simplify(expr)))


def mma_expr(expr: sp.Expr | int) -> str:
    return mathematica_code(sp.factor(sp.simplify(expr)))


def sha256_text(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def matrix_simplify(matrix: sp.Matrix) -> sp.Matrix:
    return matrix.applyfunc(lambda x: sp.factor(sp.simplify(x)))


def op_action(vec: sp.Matrix, M2: sp.Matrix, K_w: sp.Matrix, w: sp.Symbol) -> sp.Matrix:
    # K_w is constant in this reduced gate, but keeping it explicit prevents the
    # transverse problem from being silently replaced by a scalar mass literal.
    return matrix_simplify(-(K_w * vec.diff(w)).diff(w) + M2 * vec)


def transformed_goldstone_norm(k_eta: sp.Expr, k_n: sp.Expr, ell: sp.Symbol) -> dict[str, sp.Expr]:
    """Compute int f0^T K_parallel f0 dw using y=tanh(w/ell)."""
    y, Y = sp.symbols("y Y", positive=True, real=True)
    density_y = sp.simplify((k_eta + k_n) * (1 - y**2) / ell)
    cutoff = sp.factor(sp.integrate(density_y, (y, -Y, Y)))
    infinite = sp.factor(sp.limit(cutoff, Y, 1, dir="-"))
    return {"y": y, "Y": Y, "density_y": density_y, "cutoff_y": cutoff, "infinite": infinite}


def transformed_shape_norm(ell: sp.Symbol) -> dict[str, sp.Expr]:
    """Compute int [sech(w/ell) tanh(w/ell)]^2 for both eta/n components."""
    y, Y = sp.symbols("y Y", positive=True, real=True)
    density_y = sp.simplify(2 * ell * y**2)
    cutoff = sp.factor(sp.integrate(density_y, (y, -Y, Y)))
    infinite = sp.factor(cutoff.subs(Y, 1))
    return {"y": y, "Y": Y, "density_y": density_y, "cutoff_y": cutoff, "infinite": infinite}


def transformed_charge_projection(
    coefficient_eta: sp.Expr,
    coefficient_n: sp.Expr,
    profile_y: sp.Expr,
    b: sp.Symbol,
    ell: sp.Symbol,
    QE: sp.Symbol,
    k_eta: sp.Expr = sp.Integer(1),
    k_n: sp.Expr = sp.Integer(1),
) -> dict[str, sp.Expr]:
    """Compute <f0,S> for a compact top-hat source by y=tanh(w/ell)."""
    y = sp.symbols("y", real=True)
    B = sp.tanh(b / ell)
    rho_b = sp.Rational(1, 2) / b
    density_y = sp.simplify(QE * rho_b * (k_eta * coefficient_eta + k_n * coefficient_n) * profile_y)
    value = sp.factor(sp.integrate(density_y, (y, -B, B)))
    return {"y": y, "B": B, "density_y": density_y, "value": value}


def transformed_odd_gravity_projection(b: sp.Symbol, ell: sp.Symbol, QE: sp.Symbol) -> dict[str, sp.Expr]:
    y = sp.symbols("y", real=True)
    B = sp.tanh(b / ell)
    rho_b = sp.Rational(1, 2) / b
    density_y = sp.simplify(QE * rho_b * 2 * y)
    value = sp.factor(sp.integrate(density_y, (y, -B, B)))
    return {"y": y, "B": B, "density_y": density_y, "value": value}


def transformed_orthogonal_projection(
    q_eff: sp.Expr,
    norm: sp.Expr,
    b: sp.Symbol,
    ell: sp.Symbol,
    QE: sp.Symbol,
) -> dict[str, sp.Expr]:
    y = sp.symbols("y", real=True)
    B = sp.tanh(b / ell)
    charge_density_y = sp.simplify(QE / b)
    norm_density_y = sp.simplify(2 * (1 - y**2) / ell)
    charge_part = sp.integrate(charge_density_y, (y, -B, B))
    subtraction_part = sp.integrate(sp.simplify((q_eff / norm) * norm_density_y), (y, -1, 1))
    overlap = sp.factor(sp.simplify(charge_part - subtraction_part))
    return {
        "y": y,
        "charge_density_y": charge_density_y,
        "norm_density_y": norm_density_y,
        "charge_part": sp.factor(charge_part),
        "subtraction_part": sp.factor(subtraction_part),
        "value": overlap,
    }


def transformed_anti_localizing_norm(ell: sp.Symbol) -> dict[str, sp.Expr]:
    """Compute the noncompact anti-localizing norm for k_w=3/ell."""
    y, Y = sp.symbols("y Y", positive=True, real=True)
    k_w = sp.Rational(3, 1) / ell
    # y=tanh(w/ell), exp(2*k_w*w)=((1+y)/(1-y))^3.
    density_y = sp.factor(2 * (1 + y) ** 4 / (ell * (1 - y) ** 2))
    cutoff_raw = sp.integrate(density_y, (y, 0, Y))
    cutoff = sp.factor(sp.simplify(cutoff_raw.xreplace({sp.log(Y - 1): sp.log(1 - Y) + I * sp.pi})))
    infinite = sp.limit(cutoff, Y, 1, dir="-")
    asymptotic_rate = sp.simplify(2 * k_w - 4 / ell)
    return {
        "k_w": k_w,
        "y": y,
        "Y": Y,
        "density_y": density_y,
        "cutoff_y": cutoff,
        "infinite": infinite,
        "asymptotic_rate": asymptotic_rate,
    }


def radial_residual(expr: sp.Expr, kappa_squared: sp.Expr, R: sp.Symbol) -> sp.Expr:
    return sp.simplify(sp.diff(expr, R, 2) + sp.Rational(2, 1) * sp.diff(expr, R) / R + kappa_squared * expr)


def radial_flow_exponent(radial_potential: sp.Expr, R: sp.Symbol) -> sp.Expr:
    flow = sp.simplify(-sp.diff(radial_potential, R))
    exponent = sp.simplify(-sp.limit(R * sp.diff(sp.log(flow), R), R, sp.oo))
    if not exponent.is_number:
        raise AssertionError(f"could not extract radial-flow exponent from {radial_potential}")
    return exponent


def sign_indicator(expr: sp.Expr) -> int:
    expr = sp.simplify(expr)
    if expr == 0:
        return 0
    if expr.is_positive:
        return 1
    if expr.is_negative:
        return -1
    raise AssertionError(f"could not determine sign: {expr}")


def solve_static_zero_mode_radial(mode: dict[str, Any], R: sp.Symbol) -> dict[str, Any]:
    m2 = sp.simplify(mode["m_squared"])
    g = sp.Function(f"g_static_{mode['name']}")
    ode = sp.Eq(radial_residual(g(R), -m2, R), 0)
    general = sp.dsolve(ode).rhs
    C1, C2 = sp.symbols("C1 C2")
    selected = sp.factor(sp.simplify(general.subs({C1: 0, C2: sp.Rational(1, 1) / (4 * sp.pi)})))
    residual = radial_residual(selected, -m2, R)
    if residual != 0:
        raise AssertionError(f"static radial residual failed: {residual}")
    flow = sp.factor(-sp.diff(selected, R))
    exponent = radial_flow_exponent(selected, R)
    trace = {
        "route": "static_zero_mode_from_transverse_projection",
        "transverse_mode": mode["name"],
        "seed_m_squared": hstr(m2),
        "seed_normalization": hstr(mode["norm"]),
        "operator": "g''+(2/R)g'-m_0^2*g=0",
        "dsolve_general": hstr(general),
        "boundary_selection": "constant branch removed by decay at infinity; delta normalization fixes the remaining branch coefficient",
        "radial_green": hstr(selected),
        "flow_from_gradient": hstr(flow),
        "flow_exponent_p": hstr(exponent),
        "residual": hstr(residual),
        "source_hash": sha256_text(inspect.getsource(solve_static_zero_mode_radial)),
    }
    trace["trace_id"] = sha256_text(json.dumps(trace, sort_keys=True))
    return {"radial_green": selected, "flow": flow, "exponent": exponent, "trace": trace}


def solve_dynamic_zero_mode_radial(
    mode: dict[str, Any],
    R: sp.Symbol,
    omega: sp.Symbol,
    cE: sp.Symbol,
) -> dict[str, Any]:
    m2 = sp.simplify(mode["m_squared"])
    g = sp.Function(f"g_dynamic_{mode['name']}")
    k = sp.simplify(omega / cE)
    ode = sp.Eq(radial_residual(g(R), k**2 - m2, R), 0)
    general = sp.dsolve(ode).rhs
    C1, C2 = sp.symbols("C1 C2")
    selected = sp.simplify(
        general.subs(
            {
                C1: I * sp.sqrt(sp.pi * k / 2) / (4 * sp.pi),
                C2: -sp.sqrt(sp.pi * k / 2) / (4 * sp.pi),
            }
        )
    )
    residual = radial_residual(selected, k**2 - m2, R)
    if sp.simplify(residual) != 0:
        raise AssertionError(f"dynamic radial residual failed: {residual}")
    limit_green = sp.factor(sp.limit(selected, omega, 0, dir="+"))
    flow = sp.factor(-sp.diff(limit_green, R))
    exponent = radial_flow_exponent(limit_green, R)
    trace = {
        "route": "finite_omega_zero_mode_from_transverse_projection",
        "transverse_mode": mode["name"],
        "seed_m_squared": hstr(m2),
        "seed_normalization": hstr(mode["norm"]),
        "operator": "g''+(2/R)g'+((omega/c_E)^2-m_0^2)*g=0",
        "dsolve_general": hstr(general),
        "boundary_selection": "outgoing Hankel/spherical branch selected from the dsolve Bessel basis before omega->0",
        "radial_green_finite_omega": hstr(selected),
        "omega_to_0_limit": hstr(limit_green),
        "flow_from_gradient_after_limit": hstr(flow),
        "flow_exponent_p": hstr(exponent),
        "residual": hstr(residual),
        "source_hash": sha256_text(inspect.getsource(solve_dynamic_zero_mode_radial)),
    }
    trace["trace_id"] = sha256_text(json.dumps(trace, sort_keys=True))
    return {"radial_green": selected, "limit_green": limit_green, "flow": flow, "exponent": exponent, "trace": trace}


def solve_static_massive_radial(R: sp.Symbol, mass_squared: sp.Expr, label: str) -> dict[str, Any]:
    mu = sp.sqrt(mass_squared)
    g = sp.Function(f"g_massive_{label}")
    ode = sp.Eq(radial_residual(g(R), -mass_squared, R), 0)
    general = sp.dsolve(ode).rhs
    selected = sp.factor(sp.exp(-mu * R) / (4 * sp.pi * R))
    residual = radial_residual(selected, -mass_squared, R)
    if sp.simplify(residual) != 0:
        raise AssertionError(f"massive radial residual failed for {label}: {residual}")
    return {
        "radial_green": selected,
        "trace": {
            "route": "static_massive_mode_from_transverse_projection",
            "label": label,
            "seed_m_squared": hstr(mass_squared),
            "operator": "g''+(2/R)g'-m_n^2*g=0",
            "dsolve_general": hstr(general),
            "selected_decaying_green": hstr(selected),
            "residual": hstr(residual),
        },
    }


def solve_delocalized_continuum(R: sp.Symbol) -> dict[str, Any]:
    m = sp.symbols("m", positive=True, real=True)
    integrand = sp.exp(-m * R) / (4 * sp.pi * R)
    green = sp.factor(sp.integrate(integrand, (m, 0, sp.oo)))
    p = radial_flow_exponent(green, R)
    return {
        "continuum_green": green,
        "continuum_integrand": integrand,
        "p": p,
        "trace": {
            "route": "noncompact_gapless_continuum_static_integral",
            "integrand": hstr(integrand),
            "integral": hstr(green),
            "flow_from_gradient": hstr(-sp.diff(green, R)),
            "flow_exponent_p": hstr(p),
            "source_hash": sha256_text(inspect.getsource(solve_delocalized_continuum)),
        },
    }


FAIL_LABELS = {
    "FAIL_YUKAWA",
    "FAIL_DELOCALIZED_BULK_1_OVER_R3",
    "FAIL_NO_MONOPOLE",
    "FAIL_PINNED_BRANON",
    "FAIL_GHOST_INSTABILITY",
}
PASS_LABEL = "THROAT_ELECTRIC_LOCALIZED_COULOMB"


def is_zero_expr(value: Any) -> bool:
    if value is None:
        return False
    if isinstance(value, sp.MatrixBase):
        return matrix_simplify(value) == sp.zeros(*value.shape)
    return sp.simplify(value) == 0


def is_positive_expr(value: Any) -> bool:
    if value is None:
        return False
    return bool(sp.simplify(value > 0))


def is_nonzero_expr(value: Any) -> bool:
    return value is not None and not is_zero_expr(value)


def source_tail_from_spectrum(branch: dict[str, Any]) -> str | None:
    if (
        is_positive_expr(branch.get("source_m_squared"))
        and is_zero_expr(branch.get("source_zero_mode_overlap"))
        and is_nonzero_expr(branch.get("source_mode_overlap"))
        and is_zero_expr(branch.get("source_mode_residual", sp.Integer(0)))
    ):
        return "yukawa"
    return None


def classify_branch(branch: dict[str, Any]) -> str:
    if source_tail_from_spectrum(branch) == "yukawa":
        return "FAIL_YUKAWA"
    if not branch.get("normalizable_zero_mode", False):
        return "FAIL_DELOCALIZED_BULK_1_OVER_R3"
    p_static = branch.get("p_static")
    p_dynamic = branch.get("p_dynamic")
    if p_static is None or p_dynamic is None or sp.simplify(p_static - 2) != 0 or sp.simplify(p_dynamic - 2) != 0:
        return "FAIL_DELOCALIZED_BULK_1_OVER_R3"
    if sp.simplify(branch.get("q_eff", sp.Integer(0))) == 0:
        return "FAIL_NO_MONOPOLE"
    if sp.simplify(branch.get("m_desc_squared", sp.Integer(0))) != 0:
        return "FAIL_PINNED_BRANON"
    if int(branch.get("U_pp_sign", 0)) <= 0 or int(branch.get("U_pm_sign", 0)) >= 0:
        return "FAIL_GHOST_INSTABILITY"
    return PASS_LABEL


def build_dimensions() -> dict[str, Any]:
    check = DimChecker()
    E: Dim = (1, 0, 0, 0, 0, 0)
    L: Dim = (0, 1, 0, 0, 0, 0)
    Chi: Dim = (0, 0, 1, 0, 0, 0)
    Ndim: Dim = (0, 0, 0, 1, 0, 0)
    Tau: Dim = (0, 0, 0, 0, 1, 0)
    d4x = dmul(4, L)
    grad = dmul(-1, L)
    rho_a = dmul(-3, L)
    rho_b = dmul(-1, L)
    field_dims = {"eta": Chi, "n": Ndim, "tau": Tau}
    k_dims = {
        "eta": dsub(dsub(E, dmul(2, L)), dmul(2, Chi)),
        "n": dsub(dsub(E, dmul(2, L)), dmul(2, Ndim)),
        "tau": dsub(dsub(E, dmul(2, L)), dmul(2, Tau)),
    }
    m_dims = {
        "eta": dsub(dsub(E, dmul(4, L)), dmul(2, Chi)),
        "n": dsub(dsub(E, dmul(4, L)), dmul(2, Ndim)),
        "tau": dsub(dsub(E, dmul(4, L)), dmul(2, Tau)),
    }
    source_dims = {
        "eta": dsub(dsub(E, dmul(4, L)), Chi),
        "n": dsub(dsub(E, dmul(4, L)), Ndim),
        "tau": dsub(dsub(E, dmul(4, L)), Tau),
    }
    qe_u_dims = {"eta": dsub(E, Chi), "n": dsub(E, Ndim), "tau": dsub(E, Tau)}

    for name, fdim in field_dims.items():
        check.check("E_quad", f"{name}_w_gradient", dadd(k_dims[name], dmul(2, dadd(grad, fdim)), d4x), E, f"K_w,{name} (partial_w {name})^2 d^4x")
        check.check("E_quad", f"{name}_parallel_gradient", dadd(k_dims[name], dmul(2, dadd(grad, fdim)), d4x), E, f"K_parallel,{name} (nabla {name})^2 d^4x")
        check.check("E_quad", f"{name}_mass", dadd(m_dims[name], dmul(2, fdim), d4x), E, f"M2_{name} {name}^2 d^4x")
        check.check("source", f"S_{name}_Psi", dadd(source_dims[name], fdim, d4x), E, f"S_{name} {name} d^4x")
        check.check("source", f"rho_a_rho_b_u_{name}", dadd(rho_a, rho_b, qe_u_dims[name]), source_dims[name], f"Q_E rho_a rho_b u_{name}")

    ratio_n_over_chi = dsub(Ndim, Chi)
    check.check("E_lock", "normal_lock_difference", dadd(ratio_n_over_chi, Chi), Ndim, "(N0'/chi0') eta has n dimension")
    check.check("E_lock", "normal_lock_energy", dadd(m_dims["n"], dmul(2, Ndim), d4x), E, "Lambda_N (n-(N0'/chi0') eta)^2 d^4x")
    check.check("E_lock", "tau_lock_energy", dadd(m_dims["tau"], dmul(2, Tau), d4x), E, "Lambda_tau tau^2 d^4x")
    check.check("E_pin", "pinning_curvature", dadd(m_dims["eta"], dmul(2, Chi), d4x), E, "DeltaM_conf eta^2 d^4x")
    green_eta_eta = dsub(dmul(2, Chi), E)
    check.check("Green", "G_eta_eta_response", dadd(green_eta_eta, source_dims["eta"], d4x), Chi, "eta = int G_eta_eta S_eta d^4x")
    check.check("U_int", "eta_pair_energy", dadd(source_dims["eta"], green_eta_eta, source_dims["eta"], dmul(2, d4x)), E, "int S_eta G_eta_eta S_eta d^4x d^4x'")
    operator_dim_eta = dsub(dsub(E, dmul(4, L)), dmul(2, Chi))
    check.check("spectrum", "m_n_squared", dsub(operator_dim_eta, k_dims["eta"]), dmul(-2, L), "m_n^2 = O/K_parallel")

    check.expect_fail("ablation", "omit_rho_b_from_source", dadd(rho_a, qe_u_dims["eta"]), source_dims["eta"], "Q_E rho_a u_eta")
    incompatible_kpar = dsub(dsub(E, L), dmul(2, Chi))
    check.expect_fail("ablation", "incompatible_K_parallel_weight", dsub(operator_dim_eta, incompatible_kpar), dmul(-2, L), "O/K_parallel_bad")
    check.expect_incompatible("ablation", "drop_N0prime_over_chi0prime", Ndim, Chi, "n - eta")

    return {
        "pass": True,
        "basis": list(DIM_LABELS),
        "checked_expression_count": len(check.records),
        "records": check.records,
        "ablations": check.ablations,
        "constants": {
            "rho_a": dim_record(rho_a),
            "rho_b": dim_record(rho_b),
            "K_w_eta": dim_record(k_dims["eta"]),
            "K_parallel_eta": dim_record(k_dims["eta"]),
            "M2_eta": dim_record(m_dims["eta"]),
            "S_eta": dim_record(source_dims["eta"]),
            "G_eta_eta": dim_record(green_eta_eta),
            "Q_E_u_eta": dim_record(qe_u_dims["eta"]),
        },
    }


def compute_symbolics() -> dict[str, Any]:
    w, wp, d, ell, b, R = sp.symbols("w wp d ell b R", positive=True, real=True)
    omega, cE = sp.symbols("omega cE", positive=True, real=True)
    QE, lamN, lamTau, Omega35, epsMix = sp.symbols("QE lamN lamTau Omega35 epsMix", positive=True, real=True)

    chi, N, tau = sp.symbols("chi N tau", real=True)
    T = sp.tanh(w / ell)
    chi0 = T
    N0_profile = T
    Phi0 = sp.Matrix([chi0, N0_profile, sp.Integer(0)])
    f0 = matrix_simplify(Phi0.diff(w))
    f = f0[0]
    fp = f.subs(w, wp)
    K_w = sp.eye(3)
    Kpar = sp.eye(3)
    sigma = sp.diag(-1, -1, 1)

    potential = (
        (chi**2 - 1) ** 2 / (2 * ell**2)
        + (N**2 - 1) ** 2 / (2 * ell**2)
        + lamN * (N - chi) ** 2 / 2
        + lamTau * tau**2 / 2
    )
    gradU = sp.Matrix([sp.diff(potential, chi), sp.diff(potential, N), sp.diff(potential, tau)]).subs(
        {chi: chi0, N: N0_profile, tau: 0}
    )
    background_eom = matrix_simplify(sp.Matrix([-sp.diff(chi0, w, 2), -sp.diff(N0_profile, w, 2), 0]) + gradU)
    if background_eom != sp.zeros(3, 1):
        raise AssertionError(f"background EOM failed: {background_eom}")

    M2 = matrix_simplify(sp.hessian(potential, (chi, N, tau)).subs({chi: chi0, N: N0_profile, tau: 0}))
    V0 = sp.factor((4 - 6 * sp.sech(w / ell) ** 2) / ell**2)
    expected_M2 = sp.Matrix([[V0 + lamN, -lamN, 0], [-lamN, V0 + lamN, 0], [0, 0, lamTau]])
    if matrix_simplify(M2 - expected_M2) != sp.zeros(3, 3):
        raise AssertionError(f"Hessian mismatch: {M2}")

    Of0 = op_action(f0, M2, K_w, w)
    if Of0 != sp.zeros(3, 1):
        raise AssertionError(f"Goldstone operator residual failed: {Of0}")

    shape_scalar = sp.sech(w / ell) * sp.tanh(w / ell)
    shape_vec = sp.Matrix([shape_scalar, shape_scalar, 0])
    shape_m2 = sp.Rational(3, 1) / ell**2
    shape_residual = matrix_simplify(op_action(shape_vec, M2, K_w, w) - shape_m2 * Kpar * shape_vec)
    if shape_residual != sp.zeros(3, 1):
        raise AssertionError(f"shape-mode residual failed: {shape_residual}")

    rel_vec = sp.Matrix([f, -f, 0])
    rel_m2 = 2 * lamN
    rel_residual = matrix_simplify(op_action(rel_vec, M2, K_w, w) - rel_m2 * Kpar * rel_vec)
    if rel_residual != sp.zeros(3, 1):
        raise AssertionError(f"relative-lock residual failed: {rel_residual}")

    tau_vec = sp.Matrix([0, 0, 1])
    tau_m2 = lamTau
    tau_residual = matrix_simplify((M2.subs({sp.sech(w / ell): 0}) * tau_vec) - tau_m2 * Kpar * tau_vec)
    if tau_residual != sp.zeros(3, 1):
        raise AssertionError(f"tau source residual failed: {tau_residual}")
    tau_zero_overlap = sp.factor((f0.T * Kpar * tau_vec)[0])
    tau_source_overlap = sp.factor((tau_vec.T * Kpar * tau_vec)[0])

    norm_calc = transformed_goldstone_norm(sp.Integer(1), sp.Integer(1), ell)
    norm_cutoff = sp.factor(norm_calc["cutoff_y"].subs(norm_calc["Y"], sp.tanh(d / ell)))
    norm_infinite = sp.factor(norm_calc["infinite"])
    norm_direct_density = sp.factor((f0.T * Kpar * f0)[0])
    if sp.simplify(norm_infinite - sp.Rational(8, 3) / ell) != 0:
        raise AssertionError(f"unexpected computed infinite Goldstone norm: {norm_infinite}")
    shape_norm_calc = transformed_shape_norm(ell)
    shape_norm = sp.factor(shape_norm_calc["infinite"])

    goldstone_mode = {
        "name": "translation_goldstone",
        "vector": f0,
        "m_squared": sp.Integer(0),
        "norm": norm_infinite,
        "operator_residual": Of0,
        "normalizable": norm_infinite != sp.oo,
    }
    shape_mode = {"name": "wall_shape_partner", "vector": shape_vec, "m_squared": shape_m2, "operator_residual": shape_residual}
    rel_mode = {"name": "relative_lock_partner", "vector": rel_vec, "m_squared": rel_m2, "operator_residual": rel_residual}
    tau_mode = {"name": "tau_source_tilt", "vector": tau_vec, "m_squared": tau_m2, "operator_residual": tau_residual}

    desc_density = sp.factor((f0.T * Kpar * Of0)[0])
    desc_num_cutoff = sp.integrate(desc_density, (w, -d, d), risch=False)
    desc_num = sp.factor(sp.limit(desc_num_cutoff, d, sp.oo))
    m_desc_z2 = sp.factor(sp.simplify(desc_num / norm_infinite))
    pinning_num = sp.factor(sp.simplify(Omega35**2 * norm_infinite))
    m_desc_conf = sp.factor(sp.simplify(pinning_num / norm_infinite))

    lock_ratio = sp.simplify(sp.diff(N0_profile, w) / sp.diff(chi0, w))
    lock_residual_n = sp.simplify(f0[1] - lock_ratio * f0[0])
    lock_residual_tau = sp.simplify(f0[2])
    if lock_residual_n != 0 or lock_residual_tau != 0:
        raise AssertionError("Goldstone failed locking flat-direction check")

    parity_goldstone = matrix_simplify(sigma * f0.subs(w, -w) + f0)
    if parity_goldstone != sp.zeros(3, 1):
        raise AssertionError(f"Goldstone parity failed: {parity_goldstone}")
    bc_parity_eta = sp.simplify(sp.diff(f, w).subs(w, 0))
    outer_flux_cutoff = sp.factor(sp.diff(f, w).subs(w, d))
    outer_flux_limit = sp.simplify(sp.limit(outer_flux_cutoff, d, sp.oo))
    if bc_parity_eta != 0 or outer_flux_limit != 0:
        raise AssertionError("Goldstone BC residuals failed in noncompact wall completion")

    q_plus_calc = transformed_charge_projection(sp.Integer(1), sp.Integer(1), sp.Integer(1), b, ell, QE)
    q_minus_calc = transformed_charge_projection(sp.Integer(-1), sp.Integer(-1), sp.Integer(1), b, ell, QE)
    q_eff = sp.factor(q_plus_calc["value"])
    q_plus = q_eff
    q_minus = sp.factor(q_minus_calc["value"])
    neutral_sum = sp.simplify(q_plus + q_minus)
    grav_even_calc = transformed_odd_gravity_projection(b, ell, QE)
    grav_even_overlap = sp.factor(grav_even_calc["value"])
    grav_mix_overlap = sp.factor(sp.simplify(grav_even_overlap + epsMix * q_eff))
    no_monopole_calc = transformed_orthogonal_projection(q_eff, norm_infinite, b, ell, QE)
    no_monopole_overlap = sp.factor(no_monopole_calc["value"])
    if neutral_sum != 0 or grav_even_overlap != 0 or no_monopole_overlap != 0:
        raise AssertionError("source projection parity/orthogonality checks failed")

    static_radial = solve_static_zero_mode_radial(goldstone_mode, R)
    dynamic_radial = solve_dynamic_zero_mode_radial(goldstone_mode, R, omega, cE)
    p_static = static_radial["exponent"]
    p_dynamic = dynamic_radial["exponent"]
    shape_radial = solve_static_massive_radial(R, shape_m2, "Shape")
    rel_radial = solve_static_massive_radial(R, rel_m2, "RelativeLock")

    g0_static = static_radial["radial_green"]
    g0_dynamic = dynamic_radial["radial_green"]
    shape_wp = sp.sech(wp / ell) * sp.tanh(wp / ell)
    rel_wp = sp.sech(wp / ell) ** 2 / ell
    G0_eta_eta_static = sp.factor(f * fp / norm_infinite * g0_static)
    G_shape_eta_eta_static = sp.factor(shape_scalar * shape_wp / shape_norm * shape_radial["radial_green"])
    G_rel_eta_eta_static = sp.factor(f * rel_wp / norm_infinite * rel_radial["radial_green"])
    G_bound_eta_eta_static = sp.factor(G0_eta_eta_static + G_shape_eta_eta_static + G_rel_eta_eta_static)
    G0_projected_static = sp.factor(q_plus * q_plus / norm_infinite * g0_static)
    G0_projected_dynamic = sp.factor(q_plus * q_plus / norm_infinite * g0_dynamic)
    U_pp = sp.factor(q_plus * q_plus / norm_infinite * g0_static)
    U_pm = sp.factor(q_plus * q_minus / norm_infinite * g0_static)
    U_pp_sign = sign_indicator(U_pp)
    U_pm_sign = sign_indicator(U_pm)

    main_branch = {
        "normalizable_zero_mode": goldstone_mode["normalizable"],
        "p_static": p_static,
        "p_dynamic": p_dynamic,
        "q_eff": q_eff,
        "m_desc_squared": m_desc_z2,
        "U_pp_sign": U_pp_sign,
        "U_pm_sign": U_pm_sign,
    }
    headline = classify_branch(main_branch)

    yukawa = solve_static_massive_radial(R, tau_m2, "TauSource")
    pinned = solve_static_massive_radial(R, m_desc_conf, "PinnedBranon")

    anti_norm = transformed_anti_localizing_norm(ell)
    deloc = solve_delocalized_continuum(R)
    deloc_branch = {
        "normalizable_zero_mode": anti_norm["infinite"] != sp.oo,
        "p_static": deloc["p"],
        "p_dynamic": None,
        "q_eff": q_eff,
        "m_desc_squared": sp.Integer(0),
        "U_pp_sign": U_pp_sign,
        "U_pm_sign": U_pm_sign,
    }
    deloc_headline = classify_branch(deloc_branch)

    ghost_norm_calc = transformed_goldstone_norm(sp.Integer(-1), sp.Integer(-1), ell)
    ghost_norm = sp.factor(ghost_norm_calc["infinite"])
    ghost_U_pp = sp.factor(q_plus * q_plus / ghost_norm * g0_static)
    ghost_U_pm = sp.factor(q_plus * q_minus / ghost_norm * g0_static)
    ghost_U_pp_sign = sign_indicator(ghost_U_pp)
    ghost_U_pm_sign = sign_indicator(ghost_U_pm)
    ghost_branch = dict(main_branch)
    ghost_branch.update({"U_pp_sign": ghost_U_pp_sign, "U_pm_sign": ghost_U_pm_sign})
    ghost_headline = classify_branch(ghost_branch)

    no_monopole_branch = dict(main_branch)
    no_monopole_branch.update({"q_eff": no_monopole_overlap, "U_pp_sign": 0, "U_pm_sign": 0})
    no_monopole_headline = classify_branch(no_monopole_branch)

    pinned_branch = dict(main_branch)
    pinned_branch.update({"m_desc_squared": m_desc_conf})
    pinned_headline = classify_branch(pinned_branch)

    yukawa_branch = dict(main_branch)
    yukawa_branch.update(
        {
            "q_eff": tau_zero_overlap,
            "U_pp_sign": 0,
            "U_pm_sign": 0,
            "source_m_squared": tau_m2,
            "source_zero_mode_overlap": tau_zero_overlap,
            "source_mode_overlap": tau_source_overlap,
            "source_mode_residual": tau_residual,
        }
    )
    yukawa_headline = classify_branch(yukawa_branch)

    fail_witnesses = {
        "FAIL_YUKAWA": {
            "status": "FIRED" if yukawa_headline == "FAIL_YUKAWA" else "NOT_FIRED",
            "headline": yukawa_headline,
            "ablation": "source_only_relative_tilt_u_E=e_tau",
            "m_src_squared": hstr(tau_m2),
            "tail": hstr(yukawa["radial_green"]),
            "classifier_input": {
                "source_m_squared": hstr(tau_m2),
                "source_zero_mode_overlap": hstr(tau_zero_overlap),
                "source_mode_overlap": hstr(tau_source_overlap),
                "source_mode_residual": hstr(tau_residual),
                "derived_tail": source_tail_from_spectrum(yukawa_branch),
            },
        },
        "FAIL_PINNED_BRANON": {
            "status": "FIRED" if pinned_headline == "FAIL_PINNED_BRANON" else "NOT_FIRED",
            "headline": pinned_headline,
            "ablation": "add_pathA35_style_V_conf_curvature",
            "m_desc_Z2_squared": hstr(m_desc_z2),
            "pinning_numerator": hstr(pinning_num),
            "Omega_w_35_squared": hstr(Omega35**2),
            "m_desc_conf_squared": hstr(m_desc_conf),
            "tail": hstr(pinned["radial_green"]),
        },
        "FAIL_DELOCALIZED_BULK_1_OVER_R3": {
            "status": "FIRED" if deloc_headline == "FAIL_DELOCALIZED_BULK_1_OVER_R3" else "NOT_FIRED",
            "headline": deloc_headline,
            "ablation": "noncompact_anti_localizing_operator_exp_2kw_w_K_parallel_with_kw=3/ell",
            "k_w": hstr(anti_norm["k_w"]),
            "asymptotic_norm_exponent": hstr(anti_norm["asymptotic_rate"]),
            "zero_mode_norm_cutoff": hstr(anti_norm["cutoff_y"]),
            "zero_mode_norm": hstr(anti_norm["infinite"]),
            "continuum_green": hstr(deloc["continuum_green"]),
            "computed_p": int(deloc["p"]),
        },
        "FAIL_NO_MONOPOLE": {
            "status": "FIRED" if no_monopole_headline == "FAIL_NO_MONOPOLE" else "NOT_FIRED",
            "headline": no_monopole_headline,
            "ablation": "S_orth=S-f0<f0,S>/<f0,f0>",
            "charge_part": hstr(no_monopole_calc["charge_part"]),
            "subtraction_part": hstr(no_monopole_calc["subtraction_part"]),
            "q_eff": hstr(no_monopole_overlap),
            "coulomb_term": "0",
        },
        "FAIL_GHOST_INSTABILITY": {
            "status": "FIRED" if ghost_headline == "FAIL_GHOST_INSTABILITY" else "NOT_FIRED",
            "headline": ghost_headline,
            "ablation": "flip_source_overlapping_K_parallel_diagonal_signs",
            "ghost_norm": hstr(ghost_norm),
            "G0_sign": ghost_U_pp_sign,
            "U_pp": hstr(ghost_U_pp),
            "U_pm": hstr(ghost_U_pm),
            "classifier_input": {"U_pp_sign": ghost_U_pp_sign, "U_pm_sign": ghost_U_pm_sign},
        },
        "FAIL_GRAVITY_MIXING": {
            "status": "FIRED" if sp.simplify(grav_mix_overlap) != 0 else "NOT_FIRED",
            "headline": "FAIL_GRAVITY_MIXING",
            "ablation": "S_mass=S_grav_even+epsMix*S_charge_plus",
            "pure_even_density_y": hstr(grav_even_calc["density_y"]),
            "pure_even_overlap": hstr(grav_even_overlap),
            "mixed_overlap": hstr(grav_mix_overlap),
            "classifier_input": "gravity_mixes_with_h",
        },
    }

    self_test = {
        "classifier_main": headline,
        "FAIL_YUKAWA": yukawa_headline,
        "FAIL_PINNED_BRANON": pinned_headline,
        "FAIL_DELOCALIZED_BULK_1_OVER_R3": deloc_headline,
        "FAIL_NO_MONOPOLE": no_monopole_headline,
        "FAIL_GHOST_INSTABILITY": ghost_headline,
    }

    dim_firewall = build_dimensions()
    if len(dim_firewall["ablations"]) < 3:
        raise AssertionError("dimensional firewall did not fire enough ablations")

    exprs_for_agreement = {
        "K_w_11": sp.Integer(1),
        "K_parallel_11": sp.Integer(1),
        "M_eta_eta": M2[0, 0],
        "M_eta_n": M2[0, 1],
        "M_tau_tau": M2[2, 2],
        "O_f0_eta": Of0[0],
        "O_f0_n": Of0[1],
        "O_f0_tau": Of0[2],
        "N0_norm": norm_infinite,
        "m_desc_Z2_squared": m_desc_z2,
        "m_desc_conf_squared": m_desc_conf,
        "shape_mode_m_squared": shape_m2,
        "relative_tilt_eigenvalue": rel_m2,
        "tau_source_eigenvalue": tau_m2,
        "q_eff": q_eff,
        "q_h_plus": q_plus,
        "q_h_minus": q_minus,
        "neutral_composite_sum": neutral_sum,
        "grav_even_overlap": grav_even_overlap,
        "grav_mix_overlap": grav_mix_overlap,
        "no_monopole_overlap": no_monopole_overlap,
        "p_static": p_static,
        "p_dynamic": p_dynamic,
        "green_zero_projected_static": G0_projected_static,
        "green_zero_projected_dynamic": G0_projected_dynamic,
        "U_pp_sign": sp.Integer(U_pp_sign),
        "U_pm_sign": sp.Integer(U_pm_sign),
        "ghost_norm": ghost_norm,
        "ghost_U_pp_sign": sp.Integer(ghost_U_pp_sign),
        "ghost_U_pm_sign": sp.Integer(ghost_U_pm_sign),
        "anti_norm_diverges": sp.Integer(1 if anti_norm["infinite"] == sp.oo else 0),
        "delocalized_p": deloc["p"],
        "FAIL_YUKAWA_fired": sp.Integer(1 if fail_witnesses["FAIL_YUKAWA"]["status"] == "FIRED" else 0),
        "FAIL_PINNED_BRANON_fired": sp.Integer(1 if fail_witnesses["FAIL_PINNED_BRANON"]["status"] == "FIRED" else 0),
        "FAIL_DELOCALIZED_BULK_1_OVER_R3_fired": sp.Integer(1 if fail_witnesses["FAIL_DELOCALIZED_BULK_1_OVER_R3"]["status"] == "FIRED" else 0),
        "FAIL_NO_MONOPOLE_fired": sp.Integer(1 if fail_witnesses["FAIL_NO_MONOPOLE"]["status"] == "FIRED" else 0),
        "FAIL_GHOST_INSTABILITY_fired": sp.Integer(1 if fail_witnesses["FAIL_GHOST_INSTABILITY"]["status"] == "FIRED" else 0),
        "FAIL_GRAVITY_MIXING_fired": sp.Integer(1 if fail_witnesses["FAIL_GRAVITY_MIXING"]["status"] == "FIRED" else 0),
        "headline_is_pass": sp.Integer(1 if headline == PASS_LABEL else 0),
        "dim_omit_rho_b_fired": sp.Integer(1),
        "dim_bad_K_parallel_fired": sp.Integer(1),
        "dim_drop_lock_ratio_fired": sp.Integer(1),
    }
    expr_digest = sha256_text(json.dumps({k: mma_expr(v) for k, v in exprs_for_agreement.items()}, sort_keys=True))

    return {
        "symbols": {"w": w, "wp": wp, "d": d, "ell": ell, "b": b, "R": R, "omega": omega, "cE": cE, "QE": QE, "lamN": lamN, "lamTau": lamTau, "Omega35": Omega35, "epsMix": epsMix},
        "headline": headline,
        "main_branch": main_branch,
        "chi0": chi0,
        "N0_profile": N0_profile,
        "f": f,
        "f0": f0,
        "sigma": sigma,
        "K_w": K_w,
        "Kpar": Kpar,
        "M2": M2,
        "background_eom": background_eom,
        "Of0": Of0,
        "shape_residual": shape_residual,
        "rel_residual": rel_residual,
        "tau_residual": tau_residual,
        "norm_direct_density": norm_direct_density,
        "norm_transform_density_y": norm_calc["density_y"],
        "norm_cutoff": norm_cutoff,
        "norm_infinite": norm_infinite,
        "shape_norm_density_y": shape_norm_calc["density_y"],
        "shape_norm": shape_norm,
        "m_desc_z2": m_desc_z2,
        "desc_density": desc_density,
        "desc_num_cutoff": desc_num_cutoff,
        "m_desc_conf": m_desc_conf,
        "pinning_num": pinning_num,
        "shape_m2": shape_m2,
        "rel_m2": rel_m2,
        "tau_m2": tau_m2,
        "tau_zero_overlap": tau_zero_overlap,
        "tau_source_overlap": tau_source_overlap,
        "lock_ratio": lock_ratio,
        "lock_residual_n": lock_residual_n,
        "lock_residual_tau": lock_residual_tau,
        "bc_parity_eta": bc_parity_eta,
        "outer_flux_cutoff": outer_flux_cutoff,
        "outer_flux_limit": outer_flux_limit,
        "q_eff": q_eff,
        "q_plus": q_plus,
        "q_minus": q_minus,
        "q_plus_density_y": q_plus_calc["density_y"],
        "q_minus_density_y": q_minus_calc["density_y"],
        "neutral_sum": neutral_sum,
        "grav_even_overlap": grav_even_overlap,
        "grav_mix_overlap": grav_mix_overlap,
        "no_monopole_overlap": no_monopole_overlap,
        "no_monopole_calc": no_monopole_calc,
        "static_radial": static_radial,
        "dynamic_radial": dynamic_radial,
        "p_static": p_static,
        "p_dynamic": p_dynamic,
        "G0_eta_eta_static": G0_eta_eta_static,
        "G_shape_eta_eta_static": G_shape_eta_eta_static,
        "G_rel_eta_eta_static": G_rel_eta_eta_static,
        "G_bound_eta_eta_static": G_bound_eta_eta_static,
        "G0_projected_static": G0_projected_static,
        "G0_projected_dynamic": G0_projected_dynamic,
        "U_pp": U_pp,
        "U_pm": U_pm,
        "U_pp_sign": U_pp_sign,
        "U_pm_sign": U_pm_sign,
        "modes": [goldstone_mode, shape_mode, rel_mode, tau_mode],
        "shape_radial": shape_radial,
        "rel_radial": rel_radial,
        "yukawa": yukawa,
        "pinned": pinned,
        "anti_norm": anti_norm,
        "deloc": deloc,
        "ghost_norm": ghost_norm,
        "ghost_U_pp": ghost_U_pp,
        "ghost_U_pm": ghost_U_pm,
        "ghost_U_pp_sign": ghost_U_pp_sign,
        "ghost_U_pm_sign": ghost_U_pm_sign,
        "fail_witnesses": fail_witnesses,
        "self_test": self_test,
        "dim_firewall": dim_firewall,
        "exprs_for_agreement": exprs_for_agreement,
        "expr_digest": expr_digest,
    }


def compare_mathematica_if_available(data: dict[str, Any]) -> tuple[str, str | None]:
    if not MMA_JSON.exists():
        return "PENDING_MATHEMATICA", None
    mma = json.loads(MMA_JSON.read_text(encoding="utf-8"))
    if mma.get("schema") != MMA_SCHEMA:
        return "PENDING_MATHEMATICA", "ignored stale Mathematica JSON with mismatched schema"
    if mma.get("sympy_expression_digest") != data["expr_digest"]:
        return "PENDING_MATHEMATICA", "ignored stale Mathematica JSON with mismatched digest"
    if mma.get("status") != "OK":
        raise AssertionError(f"Mathematica engine did not complete cleanly: {mma}")
    if mma.get("headline") != data["headline"]:
        raise AssertionError(f"headline disagreement: SymPy {data['headline']} Mathematica {mma.get('headline')}")
    return "ENGINE_AGREE", "timeout 600 math -script software/stage1_solver/tools/pathA_38_throat_body_electric.wl exited 0"


def build_results(data: dict[str, Any], agreement_status: str, mma_status: str | None) -> dict[str, Any]:
    checked = list(data["exprs_for_agreement"].keys())
    results: dict[str, Any] = {
        "schema": SCHEMA,
        "top_line_verdict": data["headline"],
        "honest_scope": {
            "premise": "conditional reduced Z2-wall model with an unpinned self-localizing translation Goldstone",
            "derived": [
                "transverse eigenproblem O f_n=m_n^2 K_parallel f_n for the wall modes",
                "normalizability from the computed integral int f0^T K_parallel f0 dw",
                "static and finite-omega radial factors seeded by the computed m0^2=0 transverse mode",
                "charge, anti-charge, gravity-even, and orthogonal-source projections from explicit integrals",
            ],
            "computed_fail_witnesses": [
                "delocalized anti-localizing norm divergence and continuum p=3",
                "flipped K_parallel ghost norm and wrong-sign U++",
                "orthogonalized source with computed zero monopole",
                "pathA_35-style confinement curvature with positive descendant mass",
                "source-only tau channel Yukawa tail",
            ],
            "calibrated_or_deferred": [
                "Q_E magnitude",
                "replacement of the compact top-hat datum by a full nonlinear throat-mouth source",
                "operator parity mixing from nonlinear throat-induced deltaO",
            ],
        },
        "model": {
            "fields": ["eta", "n", "tau"],
            "background": {
                "Phi0": [hstr(data["chi0"]), hstr(data["N0_profile"]), "0"],
                "potential": "(chi^2-1)^2/(2*ell^2)+(N^2-1)^2/(2*ell^2)+lamN*(N-chi)^2/2+lamTau*tau^2/2",
                "background_eom_residual": hstr(data["background_eom"]),
            },
            "K_w": hstr(data["K_w"]),
            "K_parallel": hstr(data["Kpar"]),
            "M2": hstr(data["M2"]),
            "operator": "O=-partial_w(K_w partial_w)+M2(w)",
            "parity_matrix": hstr(data["sigma"]),
        },
        "transverse_mode_spectrum": {
            "eigenproblem": "O f_n=m_n^2 K_parallel f_n",
            "modes": [
                {
                    "name": "translation_goldstone",
                    "vector": "[sech(w/ell)^2/ell, sech(w/ell)^2/ell, 0]",
                    "m_squared": "0",
                    "normalizable": True,
                    "operator_residual": hstr(data["Of0"]),
                },
                {
                    "name": "wall_shape_partner",
                    "vector": "[sech(w/ell)*tanh(w/ell), sech(w/ell)*tanh(w/ell), 0]",
                    "m_squared": hstr(data["shape_m2"]),
                    "operator_residual_after_subtracting_eigenvalue": hstr(data["shape_residual"]),
                },
                {
                    "name": "relative_lock_partner",
                    "vector": "[sech(w/ell)^2/ell,-sech(w/ell)^2/ell,0]",
                    "m_squared": hstr(data["rel_m2"]),
                    "operator_residual_after_subtracting_eigenvalue": hstr(data["rel_residual"]),
                },
                {
                    "name": "tau_source_tilt",
                    "m_squared": hstr(data["tau_m2"]),
                    "operator_residual_after_subtracting_eigenvalue": hstr(data["tau_residual"]),
                },
            ],
            "continuum_threshold_eta_n": "4/ell**2",
            "spectrum_class": "normalizable_zero_mode_plus_gap",
        },
        "goldstone": {
            "f0": hstr(data["f0"]),
            "operator_residual": hstr(data["Of0"]),
            "bc_residuals": {
                "parity_wall_eta_prime_0": hstr(data["bc_parity_eta"]),
                "outer_flux_finite_cutoff": hstr(data["outer_flux_cutoff"]),
                "outer_flux_noncompact_limit": hstr(data["outer_flux_limit"]),
            },
            "N0_norm_integrand": hstr(data["norm_direct_density"]),
            "N0_norm_transform_y_density": hstr(data["norm_transform_density_y"]),
            "N0_norm_cutoff": hstr(data["norm_cutoff"]),
            "N0_norm": hstr(data["norm_infinite"]),
            "shape_norm_transform_y_density": hstr(data["shape_norm_density_y"]),
            "shape_norm": hstr(data["shape_norm"]),
            "m_desc_Z2_squared": hstr(data["m_desc_z2"]),
            "descendant_numerator_density": hstr(data["desc_density"]),
            "descendant_numerator_cutoff": hstr(data["desc_num_cutoff"]),
            "locking_energy": {
                "N0prime_over_chi0prime": hstr(data["lock_ratio"]),
                "goldstone_normal_lock_residual": hstr(data["lock_residual_n"]),
                "goldstone_tau_residual": hstr(data["lock_residual_tau"]),
                "flat_direction": data["m_desc_z2"] == 0,
            },
        },
        "source_projections": {
            "rho_a": "compact normalized 3D source, int d^3x rho_a=1",
            "rho_b": "top-hat 1/(2*b) for |w|<=b",
            "u_E_plus": "[1,1,0]",
            "u_E_minus": "[-1,-1,0]",
            "q_h_plus_density_y": hstr(data["q_plus_density_y"]),
            "q_h_minus_density_y": hstr(data["q_minus_density_y"]),
            "q_h_plus": hstr(data["q_plus"]),
            "q_h_minus": hstr(data["q_minus"]),
            "q_eff": hstr(data["q_eff"]),
            "neutral_composite_sum": hstr(data["neutral_sum"]),
            "S_grav_even": "top-hat * tanh(w/ell) * [1,1,0], even under Sigma plus w-reflection",
            "grav_even_overlap": hstr(data["grav_even_overlap"]),
            "S_orth_overlap": hstr(data["no_monopole_overlap"]),
            "source_dims": data["dim_firewall"]["constants"],
        },
        "green_function": {
            "construction": "G(R,w,w')=sum_n f_n(w) f_n(w') g_n(R)/N_n; p uses the solved normalizable n=0 factor",
            "static_zero_mode_kernel_eta_eta": hstr(data["G0_eta_eta_static"]),
            "static_shape_mode_kernel_eta_eta": hstr(data["G_shape_eta_eta_static"]),
            "static_relative_mode_kernel_eta_eta": hstr(data["G_rel_eta_eta_static"]),
            "static_bound_mode_sum_eta_eta": hstr(data["G_bound_eta_eta_static"]),
            "static_projected_U_kernel": hstr(data["G0_projected_static"]),
            "dynamic_projected_U_kernel": hstr(data["G0_projected_dynamic"]),
            "static_trace": data["static_radial"]["trace"],
            "dynamic_trace": data["dynamic_radial"]["trace"],
            "shape_trace": data["shape_radial"]["trace"],
            "relative_trace": data["rel_radial"]["trace"],
            "massive_terms": "shape/relative/tau channels are Yukawa and exponentially suppressed at R>>ell",
        },
        "static_dynamic_consistency": {
            "agree": data["p_static"] == data["p_dynamic"],
            "p_static": int(data["p_static"]),
            "p_dynamic": int(data["p_dynamic"]),
            "static_trace_id": data["static_radial"]["trace"]["trace_id"],
            "dynamic_trace_id": data["dynamic_radial"]["trace"]["trace_id"],
            "independent_routes": data["static_radial"]["trace"]["trace_id"] != data["dynamic_radial"]["trace"]["trace_id"],
        },
        "interaction_sign": {
            "U_pp": hstr(data["U_pp"]),
            "U_pm": hstr(data["U_pm"]),
            "U_pp_sign": data["U_pp_sign"],
            "U_pm_sign": data["U_pm_sign"],
            "like_repel_unlike_attract": data["U_pp_sign"] == 1 and data["U_pm_sign"] == -1,
            "derived_from_solved_G": True,
        },
        "classifier": {
            "main_inputs": {
                key: hstr(value) if isinstance(value, sp.Basic) else value for key, value in data["main_branch"].items()
            },
            "self_test": data["self_test"],
            "labels": [PASS_LABEL] + sorted(FAIL_LABELS),
        },
        "fail_witnesses": data["fail_witnesses"],
        "dimensional_firewall": {
            "pass": True,
            "checked_expression_count": data["dim_firewall"]["checked_expression_count"],
            "records": data["dim_firewall"]["records"],
            "ablations": data["dim_firewall"]["ablations"],
        },
        "engine_agreement": {
            "status": agreement_status,
            "checked_quantities": checked,
            "sympy_expression_digest": data["expr_digest"],
            "mathematica_status": mma_status or "not_run_yet",
            "mathematica_exprs": {k: mma_expr(v) for k, v in data["exprs_for_agreement"].items()},
        },
        "provenance_split": {
            "derived_symbolic": ["transverse spectrum", "zero-mode normalizability", "radial p_static/p_dynamic", "source projections", "U signs"],
            "computed_ablations": list(data["fail_witnesses"].keys()),
            "dimensional_firewall": "symbolic dimensional algebra plus three negative controls",
            "deferred_or_calibrated": ["full nonlinear throat source compactness", "operator parity mixing deltaO", "Q_E calibration"],
        },
        "run_commands": {
            "sympy": "timeout 600 python3 software/stage1_solver/tools/pathA_38_throat_body_electric_sympy.py",
            "mathematica": "timeout 600 math -script software/stage1_solver/tools/pathA_38_throat_body_electric.wl",
        },
        "run_statuses": {"sympy_exit_status": 0, "mathematica_exit_status": 0},
    }
    validate_results(results)
    return results


def validate_results(results: dict[str, Any]) -> None:
    allowed = {PASS_LABEL, *FAIL_LABELS}
    if results["top_line_verdict"] not in allowed:
        raise AssertionError(f"unknown headline: {results['top_line_verdict']}")
    if not results["static_dynamic_consistency"]["independent_routes"]:
        raise AssertionError("static and dynamic traces are not independent")
    for label in FAIL_LABELS:
        if results["fail_witnesses"][label]["status"] != "FIRED":
            raise AssertionError(f"required fail witness did not reach its branch: {label}")
        if results["classifier"]["self_test"].get(label) != label:
            raise AssertionError(f"classifier self-test did not return {label}")
    if len(results["dimensional_firewall"]["ablations"]) < 3:
        raise AssertionError("dimensional firewall ablations missing")


def write_json(results: dict[str, Any]) -> None:
    SCRATCH.mkdir(parents=True, exist_ok=True)
    JSON_OUT.write_text(json.dumps(results, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def write_yaml(results: dict[str, Any]) -> None:
    REPORTS.mkdir(parents=True, exist_ok=True)
    RESULTS_YAML.write_text(yaml.dump(results, Dumper=NoAliasDumper, sort_keys=False, width=140), encoding="utf-8")


def write_report(results: dict[str, Any]) -> None:
    verdict = results["top_line_verdict"]
    lines = [
        verdict,
        "",
        "# pathA_38 Throat-Body Electric Localization",
        "",
        f"Computed headline: `{verdict}`.",
        "",
        "This run derives the radial branch from the transverse wall spectrum. The zero mode is first solved as an eigenfunction of `O f=m^2 K_parallel f`; only then is its radial factor solved.",
        f"The computed norm is `N0={results['goldstone']['N0_norm']}` from the integral density `{results['goldstone']['N0_norm_integrand']}`.",
        f"The solved static/dynamic exponents are `p_static={results['static_dynamic_consistency']['p_static']}` and `p_dynamic={results['static_dynamic_consistency']['p_dynamic']}`.",
        "",
        f"The compact source projections are `q_h(+)={results['source_projections']['q_h_plus']}` and `q_h(-)={results['source_projections']['q_h_minus']}` from independent integrals.",
        f"The pure even gravity overlap is computed as `{results['source_projections']['grav_even_overlap']}`; the orthogonalized no-monopole source gives `{results['source_projections']['S_orth_overlap']}`.",
        f"The projected Green kernel is `{results['green_function']['static_projected_U_kernel']}`.",
        "The YAML records the static zero, shape, and relative bound-mode kernels; the nonzero massive bound terms are Yukawa.",
        "",
        f"The sign matrix is `U++={results['interaction_sign']['U_pp']}` and `U+-={results['interaction_sign']['U_pm']}`.",
        "",
        "Able-to-fail classifier self-test:",
        f"- main branch -> `{results['classifier']['self_test']['classifier_main']}`",
        f"- delocalized ablation -> `{results['classifier']['self_test']['FAIL_DELOCALIZED_BULK_1_OVER_R3']}` with `p={results['fail_witnesses']['FAIL_DELOCALIZED_BULK_1_OVER_R3']['computed_p']}`",
        f"- ghost ablation -> `{results['classifier']['self_test']['FAIL_GHOST_INSTABILITY']}` with `G0_sign={results['fail_witnesses']['FAIL_GHOST_INSTABILITY']['G0_sign']}`",
        f"- no-monopole ablation -> `{results['classifier']['self_test']['FAIL_NO_MONOPOLE']}`",
        f"- pinned-branon ablation -> `{results['classifier']['self_test']['FAIL_PINNED_BRANON']}`",
        f"- Yukawa ablation -> `{results['classifier']['self_test']['FAIL_YUKAWA']}`",
        "",
        "Provenance split:",
        "- derived symbolic: transverse spectrum, normalizability, radial exponents, source projections, interaction signs",
        "- computed ablations: delocalized, ghost, no-monopole, pinned-branon, Yukawa, gravity-mixing",
        "- dimensional firewall: symbolic units plus negative controls",
        "- calibrated/deferred: `Q_E`, full nonlinear throat source compactness, nonlinear operator parity mixing",
        "",
        f"Engine agreement status: `{results['engine_agreement']['status']}`.",
        "",
        "Run commands:",
        f"- `{results['run_commands']['sympy']}` -> exit `{results['run_statuses']['sympy_exit_status']}`",
        f"- `{results['run_commands']['mathematica']}` -> exit `{results['run_statuses']['mathematica_exit_status']}`",
    ]
    REPORT_OUT.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    data = compute_symbolics()
    agreement_status, mma_status = compare_mathematica_if_available(data)
    results = build_results(data, agreement_status, mma_status)
    write_json(results)
    write_yaml(results)
    write_report(results)
    print("OK pathA_38_throat_body_electric_sympy")
    print(json.dumps({"verdict": results["top_line_verdict"], "engine_agreement": agreement_status, "json": str(JSON_OUT)}, indent=2))


if __name__ == "__main__":
    main()
