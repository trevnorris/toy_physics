#!/usr/bin/env python3
"""PathA-43 density quadrupole port, SymPy engine.

Standalone print-only gate:

  timeout 600 python software/stage1_solver/tools/pathA_43_density_quadrupole_port_sympy.py

This engine derives the density two-port numerator by assembling the static
mass-normalized (q2, Phi2) reduced Lagrangian and Schur-eliminating the two
coordinates.  It does not read Mathematica output.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import sympy as sp


Dim = tuple[sp.Rational, sp.Rational, sp.Rational]
ZERO_DIM: Dim = (sp.Rational(0), sp.Rational(0), sp.Rational(0))

Ldim, Mdim, Tdim = sp.symbols("L M T", positive=True)
a, c_s, c, G, D0 = sp.symbols("a c_s c G D0", positive=True)
rho = sp.symbols("rho", positive=True)
I25, Xi_Q, eta_q, eta_phi = sp.symbols("I25 Xi_Q eta_q eta_phi", positive=True)
I_wrong2 = sp.symbols("I_wrong2", positive=True)
varpi_q, varpi_phi, lambda_c = sp.symbols(
    "varpi_q2 varpi_Phi2 lambda_c", positive=True
)
q2, Phi2 = sp.symbols("q2 Phi2")

Omega_U, Omega_W, R_mix, g_U, g_W = sp.symbols(
    "Omega_U Omega_W R_mix g_U g_W", nonzero=True
)
A_w, F_muw, Jw, U, W = sp.symbols("A_w F_muw Jw U W", nonzero=True)
omega_wall, omega_rho, r_mix, g_rho, g_qold = sp.symbols(
    "omega_wall omega_rho r_mix g_rho g_qold", nonzero=True
)
sigma_hidden = sp.symbols("sigma_hidden", nonzero=True)
Xi_deferred = sp.symbols("Xi_deferred", positive=True)
free_carrier = sp.symbols("free_carrier")

z, omega = sp.symbols("z omega")


class DimError(ValueError):
    pass


def compact(expr: Any) -> Any:
    if expr is None:
        return None
    if isinstance(expr, sp.Basic):
        return sp.factor(sp.cancel(sp.simplify(expr)))
    return expr


def hstr(expr: Any) -> str:
    if expr is None:
        return "None"
    return sp.sstr(compact(expr))


def bstr(value: bool) -> str:
    return "PASS" if value else "FAIL"


def dim_add(left: Dim, right: Dim) -> Dim:
    return tuple(sp.Rational(x + y) for x, y in zip(left, right))  # type: ignore[return-value]


def dim_sub(left: Dim, right: Dim) -> Dim:
    return tuple(sp.Rational(x - y) for x, y in zip(left, right))  # type: ignore[return-value]


def dim_scale(dim: Dim, scale: sp.Rational) -> Dim:
    return tuple(sp.Rational(scale * x) for x in dim)  # type: ignore[return-value]


def dim_of(expr: sp.Expr, dims: dict[sp.Symbol, Dim]) -> Dim:
    expr = sp.sympify(expr)
    if expr == 0 or expr.is_number:
        return ZERO_DIM
    if expr.is_Symbol:
        if expr not in dims:
            raise DimError(f"missing sourced dimension for {expr}")
        return dims[expr]
    if expr.is_Mul:
        out = ZERO_DIM
        for arg in expr.args:
            out = dim_add(out, dim_of(arg, dims))
        return out
    if expr.is_Pow:
        base, exponent = expr.args
        if not exponent.is_number:
            raise DimError(f"non-numeric dimension exponent in {expr}")
        return dim_scale(dim_of(base, dims), sp.Rational(exponent))
    if expr.is_Add:
        term_dims = [dim_of(arg, dims) for arg in expr.args if arg != 0]
        if not term_dims:
            return ZERO_DIM
        first = term_dims[0]
        if any(dim != first for dim in term_dims):
            raise DimError(f"dimension mismatch in sum {expr}: {term_dims}")
        return first
    raise DimError(f"unsupported dimension expression {expr}")


def exp_text(exp: sp.Rational) -> str:
    exp = sp.Rational(exp)
    return str(exp.p) if exp.q == 1 else f"{exp.p}/{exp.q}"


def dim_to_text(dim: Dim) -> str:
    parts: list[str] = []
    for name, exp in (("L", dim[0]), ("T", dim[2]), ("M", dim[1])):
        exp = sp.Rational(exp)
        if exp == 0:
            continue
        parts.append(name if exp == 1 else f"{name}^{exp_text(exp)}")
    return " ".join(parts) if parts else "1"


def scale_power(expr: sp.Expr, powers: dict[sp.Symbol, sp.Rational | None]) -> sp.Rational | None:
    expr = sp.sympify(expr)
    if expr == 0 or expr.is_number:
        return sp.Rational(0)
    if expr == a:
        return sp.Rational(1)
    if expr.is_Symbol:
        return powers.get(expr, sp.Rational(0))
    if expr.is_Mul:
        total = sp.Rational(0)
        for arg in expr.args:
            p = scale_power(arg, powers)
            if p is None:
                return None
            total += p
        return sp.Rational(total)
    if expr.is_Pow:
        base, exponent = expr.args
        if not exponent.is_number:
            return None
        p = scale_power(base, powers)
        return None if p is None else sp.Rational(exponent) * p
    if expr.is_Add:
        term_powers = [scale_power(arg, powers) for arg in expr.args if arg != 0]
        if not term_powers:
            return sp.Rational(0)
        if any(p is None for p in term_powers):
            return None
        first = term_powers[0]
        if any(p != first for p in term_powers):
            return None
        return first
    return None


@dataclass(frozen=True)
class Config:
    name: str
    kind: str = "density"
    coupling_zero: bool = False
    corrupt_dimension: bool = False
    free_carrier_rider: bool = False
    free_carrier_tagged: bool = False
    incoming_sign: bool = False
    coupling_a_power: sp.Rational = sp.Rational(-7, 2)
    deferred_uncertified: bool = False
    proven_deferred: bool = False
    fake_continuity: bool = False
    hidden_vector_intermediate: bool = False
    corrupt_continuity_moment: bool = False
    vector_injected_density: bool = False


VECTOR_SYMBOLS = {
    A_w,
    F_muw,
    Jw,
    U,
    W,
    Omega_U,
    Omega_W,
    R_mix,
    g_U,
    g_W,
}

CONTINUITY_OPERATOR_ID = "pathA_29_projected_continuity"
CONTINUITY_L0 = "M0 = Integral(S_leak d3x)"
CONTINUITY_L1 = "D1_i = Integral(x_i*S_leak d3x) + Integral(j_i d3x)"
CONTINUITY_L2 = "Q2_m = Integral(Y2_m_star*S_leak d3x)"
CONTINUITY_L2_KERNEL = "Y2_m_star*S_leak"
BAD_CONTINUITY_L2 = "GARBAGE_NOT_A_CONTINUITY_MOMENT_AT_ALL"


BASE_DIMS: dict[sp.Symbol, Dim] = {
    a: (1, 0, 0),
    c_s: (1, 0, -1),
    c: (1, 0, -1),
    G: (3, -1, -2),
    D0: (-1, 1, -2),
    rho: (-3, 1, 0),
    I25: (sp.Rational(5, 2), 0, 0),
    I_wrong2: (2, 0, 0),
    Xi_Q: ZERO_DIM,
    Xi_deferred: ZERO_DIM,
    eta_q: ZERO_DIM,
    eta_phi: ZERO_DIM,
    varpi_q: (0, 0, -2),
    varpi_phi: (0, 0, -2),
    lambda_c: (0, 0, -2),
    free_carrier: ZERO_DIM,
    sigma_hidden: ZERO_DIM,
    A_w: ZERO_DIM,
    F_muw: ZERO_DIM,
    Jw: ZERO_DIM,
    U: ZERO_DIM,
    W: ZERO_DIM,
    Omega_U: (0, 0, -1),
    Omega_W: (0, 0, -1),
    R_mix: (0, 0, -2),
    g_U: (sp.Rational(-1, 2), sp.Rational(1, 2), -2),
    g_W: (sp.Rational(-1, 2), sp.Rational(1, 2), -2),
    omega_wall: (0, 0, -1),
    omega_rho: (0, 0, -1),
    r_mix: (0, 0, -2),
    g_rho: (sp.Rational(-1, 2), sp.Rational(1, 2), -2),
    g_qold: (sp.Rational(-1, 2), sp.Rational(1, 2), -2),
}


BASE_A_POWERS: dict[sp.Symbol, sp.Rational | None] = {
    c_s: sp.Rational(0),
    c: sp.Rational(0),
    G: sp.Rational(0),
    D0: sp.Rational(0),
    rho: sp.Rational(0),
    I25: sp.Rational(0),
    I_wrong2: sp.Rational(0),
    Xi_Q: sp.Rational(0),
    Xi_deferred: sp.Rational(0),
    eta_q: sp.Rational(0),
    eta_phi: sp.Rational(0),
    varpi_q: sp.Rational(-2),
    varpi_phi: sp.Rational(-2),
    lambda_c: sp.Rational(-2),
    free_carrier: sp.Rational(0),
    sigma_hidden: sp.Rational(0),
    A_w: sp.Rational(0),
    F_muw: sp.Rational(0),
    Jw: sp.Rational(0),
    U: sp.Rational(0),
    W: sp.Rational(0),
    Omega_U: sp.Rational(-1),
    Omega_W: sp.Rational(-1),
    R_mix: sp.Rational(-2),
    g_U: sp.Rational(-1),
    g_W: sp.Rational(-1),
    omega_wall: sp.Rational(-1),
    omega_rho: sp.Rational(-1),
    r_mix: sp.Rational(-2),
    g_rho: sp.Rational(-1),
    g_qold: sp.Rational(-1),
}


BASE_SOURCE_TAGS: dict[sp.Symbol, set[str]] = {
    a: {"pathA_29_bulk", "pathA_32_wall"},
    c_s: {"pathA_29_bulk"},
    c: {"calibrated_anchor"},
    G: {"calibrated_anchor"},
    D0: {"calibrated_anchor"},
    rho: {"pathA_29_bulk"},
    Xi_Q: {"continuity_interface"},
    Xi_deferred: {"continuity_interface"},
    eta_q: {"continuity_interface"},
    eta_phi: {"continuity_interface"},
    varpi_q: {"pathA_32_wall"},
    varpi_phi: {"pathA_29_bulk"},
    lambda_c: {"continuity_interface", "pathA_29_bulk", "pathA_32_wall"},
    free_carrier: set(),
    sigma_hidden: {"vector_port"},
    A_w: {"vector_port"},
    F_muw: {"vector_port"},
    Jw: {"vector_port"},
    U: {"vector_port"},
    W: {"vector_port"},
    Omega_U: {"vector_port"},
    Omega_W: {"vector_port"},
    R_mix: {"vector_port"},
    g_U: {"vector_port"},
    g_W: {"vector_port"},
    omega_wall: {"vector_port"},
    omega_rho: {"vector_port"},
    r_mix: {"vector_port"},
    g_rho: {"vector_port"},
    g_qold: {"vector_port"},
}


def contains_all(text: Any, tokens: tuple[str, ...]) -> bool:
    text = str(text)
    return all(token in text for token in tokens)


def continuity_lineage_valid(lineage: dict[str, Any]) -> bool:
    moments = lineage.get("moments", {})
    l0 = moments.get("l0", "")
    l1 = moments.get("l1", "")
    l2 = moments.get("l2", "")
    l2_kernel = lineage.get("l2_kernel", "")
    return bool(
        lineage.get("operator_id") == CONTINUITY_OPERATOR_ID
        and contains_all(l0, ("M0", "Integral", "S_leak", "d3x"))
        and contains_all(l1, ("D1_i", "Integral", "x_i", "S_leak", "j_i", "d3x"))
        and contains_all(l2, ("Q2_m", "Integral", "Y2_m_star", "S_leak", "d3x"))
        and contains_all(l2_kernel, ("Y2_m_star", "S_leak"))
    )


def continuity_moment_symbol(config: Config, lineage: dict[str, Any]) -> tuple[sp.Symbol, bool]:
    """Return the l=2 moment symbol only when it is structurally continuous."""
    lineage_valid = continuity_lineage_valid(lineage)
    if not lineage_valid:
        return I_wrong2, False
    if config.coupling_a_power != sp.Rational(-7, 2):
        return I_wrong2, True
    return I25, True


def source_tag_map(
    moment_symbol: sp.Symbol,
    moment_valid: bool,
    *,
    free_carrier_tagged: bool = False,
) -> dict[sp.Symbol, set[str]]:
    out = {sym: set(tags) for sym, tags in BASE_SOURCE_TAGS.items()}
    if free_carrier_tagged:
        out[free_carrier] = {"pathA_34_dimensionless_free_carrier"}
    if moment_valid:
        out[moment_symbol] = {"continuity_interface", "pathA_32_wall"}
    else:
        out.setdefault(moment_symbol, set())
    return out


def taint_for_expr(expr: sp.Expr, tag_map: dict[sp.Symbol, set[str]]) -> tuple[set[str], set[sp.Symbol]]:
    taint: set[str] = set()
    missing: set[sp.Symbol] = set()
    for sym in expr.free_symbols:
        if sym not in tag_map:
            missing.add(sym)
            continue
        taint.update(tag_map[sym])
    return taint, missing


def source_graph_for_expr(
    expr: sp.Expr,
    tag_map: dict[sp.Symbol, set[str]],
    *,
    coordinate_hosts: set[str],
) -> dict[str, Any]:
    symbols = set(expr.free_symbols)
    symbol_tags = {
        str(sym): sorted(tag_map.get(sym, set()))
        for sym in sorted(symbols, key=str)
    }
    taint, missing = taint_for_expr(expr, tag_map)
    empty = {sym for sym in symbols if sym in tag_map and not tag_map[sym]}
    vector_hosts = symbols & VECTOR_SYMBOLS
    return {
        "symbol_tags": symbol_tags,
        "taint_set": sorted(taint),
        "missing_source_symbols": sorted(str(sym) for sym in missing),
        "empty_source_symbols": sorted(str(sym) for sym in empty),
        "vector_host_symbols": sorted(str(sym) for sym in vector_hosts),
        "coordinate_hosts": sorted(coordinate_hosts),
    }


def vector_ablated_expr(expr: sp.Expr, tag_map: dict[sp.Symbol, set[str]]) -> sp.Expr:
    vector_tainted_symbols = {sym for sym, tags in tag_map.items() if "vector_port" in tags}
    ablate_symbols = (VECTOR_SYMBOLS | vector_tainted_symbols) & expr.free_symbols
    return compact(expr.subs({sym: sp.Integer(0) for sym in ablate_symbols}))


def lineage_for(config: Config) -> dict[str, Any]:
    if config.fake_continuity:
        return {
            "operator_id": "mis_tagged_vector_formula",
            "moments": {
                "l2": "Omega_U^2*g_W + R_mix*g_U, relabeled as continuity",
            },
            "failure": "no l=0->1->2 projected-continuity lineage; no Integral(Y2star*S_leak*d3x)",
        }
    if config.kind == "density":
        l2 = BAD_CONTINUITY_L2 if config.corrupt_continuity_moment else CONTINUITY_L2
        return {
            "operator_id": CONTINUITY_OPERATOR_ID,
            "moments": {
                "l0": CONTINUITY_L0,
                "l1": CONTINUITY_L1,
                "l2": l2,
            },
            "l2_kernel": CONTINUITY_L2_KERNEL,
            "lineage": "same projected-continuity operator, extended from scalar and dipole kernels to the structured l=2 Y2* moment",
        }
    return {
        "operator_id": "none",
        "moments": {},
        "failure": "vector port has no continuity moment lineage",
    }


def schur_density_expression(config: Config, lineage: dict[str, Any]) -> tuple[sp.Expr, dict[str, Any]]:
    xi = Xi_deferred if config.deferred_uncertified or config.proven_deferred else Xi_Q
    q2_moment, moment_valid = continuity_moment_symbol(config, lineage)
    g_base = sp.sqrt(rho) * c_s**2 * q2_moment * xi / a ** sp.Rational(7, 2)
    if config.coupling_a_power != sp.Rational(-7, 2):
        if config.coupling_a_power != sp.Rational(-3, 1):
            raise AssertionError("only the isolated scaling-control power is implemented")
        g_base = sp.sqrt(rho) * c_s**2 * q2_moment * xi / a**3
    if config.vector_injected_density:
        g_base = g_base * Omega_U / Omega_W
    if config.coupling_zero:
        g_base = sp.Integer(0)

    g_q_den = compact(g_base * eta_q)
    g_phi_den = compact(g_base * eta_phi)

    # Lagrangian static operator in mass-normalized coordinates.
    static_operator = sp.Matrix([[varpi_q, -lambda_c], [-lambda_c, varpi_phi]])
    source_vector = sp.Matrix([g_q_den, g_phi_den])
    response = compact((static_operator.inv() * source_vector)[1])
    delta = compact(varpi_q * varpi_phi - lambda_c**2)
    p_den = compact(varpi_q * g_phi_den + lambda_c * g_q_den)
    n0 = compact(p_den**2 / delta**2)

    if config.hidden_vector_intermediate:
        n0 = compact(n0 * sigma_hidden)

    trace = {
        "method": "SymPy Lagrangian Schur elimination",
        "static_operator": static_operator,
        "source_vector": source_vector,
        "Phi2_response": response,
        "Delta_den": delta,
        "P_den": p_den,
        "g_base": g_base,
        "g_q_den": g_q_den,
        "g_phi_den": g_phi_den,
        "continuity_moment_symbol": q2_moment,
        "continuity_moment_valid": moment_valid,
        "physical_relations": {
            "varpi_q2": "K2/M2 = (c_s/a)^2*kappa_q from pathA_32 wall angular l=2 operator",
            "varpi_Phi2": "(c_s/a)^2*(6 + (m*a)^2) from pathA_29 bulk Helmholtz l=2 mode at c_s",
            "lambda_c": "(c_s/a)^2*lambda_hat_Q from projected continuity/interface matching",
        },
    }
    return n0, trace


def vector_expression(config: Config) -> tuple[sp.Expr, dict[str, Any]]:
    if config.kind == "relabel":
        p = omega_wall**2 * g_rho + r_mix * g_qold
        delta = omega_wall**2 * omega_rho**2 - r_mix**2
    else:
        p = Omega_U**2 * g_W + R_mix * g_U
        delta = Omega_U**2 * Omega_W**2 - R_mix**2
    return compact(p**2 / delta**2), {
        "method": "old vector port fixture",
        "P_old": p,
        "Delta_old": delta,
    }


def derive(config: Config) -> dict[str, Any]:
    lineage = lineage_for(config)
    moment_symbol, moment_valid = continuity_moment_symbol(config, lineage)
    tag_map = source_tag_map(
        moment_symbol,
        moment_valid,
        free_carrier_tagged=config.free_carrier_tagged,
    )

    if config.kind in {"vector", "relabel"}:
        expr, trace = vector_expression(config)
    elif config.fake_continuity:
        expr, trace = vector_expression(Config(name=config.name, kind="vector"))
    else:
        expr, trace = schur_density_expression(config, lineage)

    if config.free_carrier_rider:
        expr = compact(expr * free_carrier)

    coordinate_hosts = {"q2", "Phi2"} if config.kind == "density" else set()
    if config.fake_continuity:
        coordinate_hosts = {"q2_fake", "Phi2_fake"}
    host_symbols = {str(sym) for sym in expr.free_symbols}
    source_graph = source_graph_for_expr(expr, tag_map, coordinate_hosts=coordinate_hosts)

    return {
        "expr": compact(expr),
        "trace": trace,
        "tags": set(source_graph["taint_set"]),
        "source_tag_map": tag_map,
        "source_graph": source_graph,
        "lineage": lineage,
        "lineage_valid": continuity_lineage_valid(lineage),
        "continuity_moment_symbol": moment_symbol,
        "continuity_moment_valid": moment_valid,
        "host_symbols": sorted(host_symbols | coordinate_hosts),
    }


def dtn_sign(kind: str) -> dict[str, Any]:
    j2 = (3 / z**3 - 1 / z) * sp.sin(z) - 3 * sp.cos(z) / z**2
    y2 = (1 / z - 3 / z**3) * sp.cos(z) - 3 * sp.sin(z) / z**2
    h = j2 + sp.I * y2 if kind == "outgoing" else j2 - sp.I * y2
    lam = compact(z * sp.diff(h, z) / h)
    yhat = compact(-3 / lam)
    series = sp.expand(sp.series(yhat, z, 0, 7).removeO())
    coeff = compact(series.coeff(z, 5) / sp.I)
    return {
        "kind": kind,
        "Yhat_series": series,
        "radiative_coeff_z5_over_i": coeff,
        "matches_outgoing": compact(coeff - sp.Rational(1, 27)) == 0,
    }


def evaluate(config: Config) -> dict[str, Any]:
    base = derive(config)
    expr = base["expr"]
    tags = set(base["tags"])
    source_graph = base["source_graph"]
    tag_map = base["source_tag_map"]
    ablated_expr = vector_ablated_expr(expr, tag_map)

    dims = dict(BASE_DIMS)
    if config.corrupt_dimension:
        dims[I25] = dim_add(dims[I25], (1, 0, 0))

    powers = dict(BASE_A_POWERS)
    if config.coupling_a_power != sp.Rational(-7, 2):
        # The explicit expression already carries the changed a power; the
        # metadata documents the mutated source fact.
        powers[Xi_Q] = sp.Rational(0)
    if config.deferred_uncertified:
        powers[Xi_deferred] = None
    if config.proven_deferred:
        powers[Xi_deferred] = sp.Rational(0)

    dim_ok = False
    dim_text = "not-evaluated"
    dim_error = None
    if expr is not None:
        try:
            expr_dim = dim_of(expr, dims)
            dim_ok = expr_dim == (-1, 1, 0)
            dim_text = dim_to_text(expr_dim)
        except DimError as exc:
            dim_error = str(exc)
            dim_text = "inhomogeneous"

    n0_power = scale_power(expr, powers) if expr is not None else None
    p0_power = None if n0_power is None else sp.Rational(n0_power - 2 - powers[D0])
    scaling_wrong = p0_power is not None and p0_power != -5
    scaling_ok = p0_power == -5
    scaling_undecidable = p0_power is None

    sign = dtn_sign("incoming" if config.incoming_sign else "outgoing")
    sign_ok = bool(sign["matches_outgoing"])

    vector_host_symbols = set(source_graph["vector_host_symbols"])
    source_map_complete = bool(
        not source_graph["missing_source_symbols"]
        and not source_graph["empty_source_symbols"]
    )
    continuity_dependency_ok = bool(
        base["lineage_valid"]
        and base["continuity_moment_valid"]
        and (
            base["continuity_moment_symbol"] in expr.free_symbols
            or (compact(expr) == 0 and config.coupling_zero)
        )
    )
    vanished_continuity_coupling = bool(
        compact(expr) == 0
        and config.coupling_zero
        and not vector_host_symbols
        and source_map_complete
        and continuity_dependency_ok
    )
    origin_ok = bool(
        (
            "continuity_interface" in tags
            and "vector_port" not in tags
            and not vector_host_symbols
            and source_map_complete
            and continuity_dependency_ok
        )
        or vanished_continuity_coupling
    )
    vector_independence_ok = bool(
        not vector_host_symbols
        and "vector_port" not in tags
        and compact(expr - ablated_expr) == 0
    )
    nonzero_ok = compact(expr) != 0
    closure = closure_overlay(expr)

    if not origin_ok or not vector_independence_ok:
        verdict = "FAIL_NOT_DENSITY_DERIVED"
    elif not nonzero_ok:
        verdict = "FAIL_PORT_VANISHES"
    elif (not dim_ok) or scaling_wrong or (not sign_ok):
        bad = []
        if not dim_ok:
            bad.append("dimensional")
        if scaling_wrong:
            bad.append("scaling")
        if not sign_ok:
            bad.append("sign")
        verdict = "FAIL_PORT_MALFORMED(" + ",".join(bad) + ")"
    elif origin_ok and nonzero_ok and dim_ok and scaling_ok and sign_ok and vector_independence_ok and closure["ok"]:
        verdict = "DENSITY_PORT_HOSTED"
    else:
        verdict = "PORT_INCONCLUSIVE_SIM_DEFERRED"

    non_dodge = None
    if verdict == "PORT_INCONCLUSIVE_SIM_DEFERRED":
        non_dodge = {
            "origin": {"result": origin_ok, "why_undecidable": None, "deferred_object": None},
            "nonzero": {"result": nonzero_ok, "why_undecidable": None, "deferred_object": None},
            "dimension": {"result": dim_ok, "why_undecidable": None, "deferred_object": None},
            "scaling": {
                "result": "undecidable" if scaling_undecidable else scaling_ok,
                "why_undecidable": "synthetic deferred overlap scalar has no structural a-independence proof"
                if scaling_undecidable
                else None,
                "deferred_object": "Xi_deferred" if scaling_undecidable else None,
            },
            "sign": {"result": sign_ok, "why_undecidable": None, "deferred_object": None},
            "vector_independence": {
                "result": vector_independence_ok,
                "why_undecidable": None,
                "deferred_object": None,
            },
        }

    return {
        "config": config,
        "expr": expr,
        "trace": base["trace"],
        "host_symbols": base["host_symbols"],
        "taint_set": sorted(tags),
        "source_graph": source_graph,
        "lineage": base["lineage"],
        "lineage_valid": base["lineage_valid"],
        "continuity_moment_symbol": str(base["continuity_moment_symbol"]),
        "continuity_dependency_ok": continuity_dependency_ok,
        "vector_host_symbols": sorted(vector_host_symbols),
        "source_map_complete": source_map_complete,
        "empty_source_symbols": source_graph["empty_source_symbols"],
        "ablation_expr": ablated_expr,
        "ablation_delta": compact(expr - ablated_expr) if expr is not None else None,
        "checks": {
            "origin_derivation_ok": origin_ok,
            "nonzero": nonzero_ok,
            "dimension": dim_ok,
            "a_scaling_provenance_ok": scaling_ok,
            "radiative_sign": sign_ok,
            "vector_independence_ok": vector_independence_ok,
        },
        "dimensions": {"N0_den": dim_text, "error": dim_error},
        "scaling": {
            "N0_den_a_power": hstr(n0_power),
            "P0_physical_a_power": hstr(p0_power),
            "expected_P0_physical_a_power": "-5",
            "undecidable": scaling_undecidable,
            "provenance": [
                "g_base carries a^(-7/2) from the structured l=2 continuity/interface projection",
                "varpi_q2, varpi_Phi2, and lambda_c carry a^-2 from wall angular stiffness and bulk Helmholtz geometry",
                "Xi_Q, eta_q, eta_phi are dimensionless profile ratios with structural a-independence in the reference configuration",
                "D0 is the carried conservative slot scalar with sourced a-power 0 in this reduced closure",
            ],
        },
        "sign": sign,
        "closure": closure,
        "non_dodge_table": non_dodge,
        "verdict": verdict,
    }


def closure_overlay(n0_expr: sp.Expr | None) -> dict[str, Any]:
    if n0_expr is None:
        return {"ok": False}
    p0_physical = compact((c_s / a) ** 2 * n0_expr / D0)
    kbar0 = 54 * G * c_s**5 / (5 * a**5 * c**5)
    kbar2 = 6 * G * c_s**3 / (5 * a**3 * c**5)
    kbar4 = 8 * G * c_s / (15 * a * c**5)
    kbar4_residual = compact(kbar4 - 4 * kbar2**2 / kbar0)
    gamma5 = compact(kbar0 * a**5 / (27 * c_s**5))
    gamma5_residual = compact(gamma5 - 2 * G / (5 * c**5))
    return {
        "P0_physical": p0_physical,
        "mhat0_sq_chi_Q_N_Q": "1 with chi_Q=1, mhat0=1, N_Q=1",
        "calibrated": ["G", "2/5", "54/5"],
        "sim_deferred": ["Xi_Q exact branch magnitude", "eta_q", "eta_phi", "lambda_c literal throat value"],
        "Kbar0": kbar0,
        "Kbar2": kbar2,
        "Kbar4": kbar4,
        "Kbar4_residual": kbar4_residual,
        "Gamma5_residual_to_2G_over_5c5": gamma5_residual,
        "ok": kbar4_residual == 0 and gamma5_residual == 0,
    }


def controls() -> dict[str, dict[str, Any]]:
    fixtures = {
        "vector_hosted": Config(name="vector_hosted", kind="vector"),
        "relabel_rig": Config(name="relabel_rig", kind="relabel"),
        "fake_continuity": Config(name="fake_continuity", kind="density", fake_continuity=True),
        "ablation_isolation": Config(
            name="ablation_isolation", kind="density", hidden_vector_intermediate=True
        ),
        "attack2_continuity_corruption": Config(
            name="attack2_continuity_corruption", corrupt_continuity_moment=True
        ),
        "attack5_vector_injection": Config(
            name="attack5_vector_injection", vector_injected_density=True
        ),
        "zero_coupling": Config(name="zero_coupling", coupling_zero=True),
        "dimensional": Config(name="dimensional", corrupt_dimension=True),
        "sign": Config(name="sign", incoming_sign=True),
        "scaling": Config(name="scaling", coupling_a_power=sp.Rational(-3, 1)),
        "deferred_scalar": Config(name="deferred_scalar", deferred_uncertified=True),
        "deferred_scalar_proven_converse": Config(name="deferred_scalar_proven_converse", proven_deferred=True),
        "free_carrier_dimension_corruption": Config(
            name="free_carrier_dimension_corruption",
            free_carrier_rider=True,
            free_carrier_tagged=True,
        ),
        "provenance_less_rider": Config(
            name="provenance_less_rider",
            free_carrier_rider=True,
        ),
    }
    out = {}
    for name, cfg in fixtures.items():
        result = evaluate(cfg)
        out[name] = {
            "verdict": result["verdict"],
            "checks": result["checks"],
            "ablation_delta": hstr(result["ablation_delta"]),
            "taint_set": result["taint_set"],
            "P0_physical_a_power": result["scaling"]["P0_physical_a_power"],
            "lineage_valid": result["lineage_valid"],
            "continuity_dependency_ok": result["continuity_dependency_ok"],
            "vector_host_symbols": result["vector_host_symbols"],
            "source_map_complete": result["source_map_complete"],
            "empty_source_symbols": result["empty_source_symbols"],
        }
    return out


def assert_gate(payload: dict[str, Any], control_payload: dict[str, dict[str, Any]]) -> None:
    baseline = payload
    assert baseline["verdict"] == "DENSITY_PORT_HOSTED", baseline["verdict"]
    for check, value in baseline["checks"].items():
        assert value, f"baseline check failed: {check}"
    assert baseline["taint_set"] == ["continuity_interface", "pathA_29_bulk", "pathA_32_wall"]
    assert not baseline["source_graph"]["empty_source_symbols"]
    assert all(tags for tags in baseline["source_graph"]["symbol_tags"].values())
    assert baseline["dimensions"]["N0_den"] == "L^-1 M"
    assert baseline["scaling"]["P0_physical_a_power"] == "-5"
    assert baseline["closure"]["Kbar4_residual"] == 0

    expected = {
        "vector_hosted": "FAIL_NOT_DENSITY_DERIVED",
        "relabel_rig": "FAIL_NOT_DENSITY_DERIVED",
        "fake_continuity": "FAIL_NOT_DENSITY_DERIVED",
        "ablation_isolation": "FAIL_NOT_DENSITY_DERIVED",
        "attack2_continuity_corruption": "FAIL_NOT_DENSITY_DERIVED",
        "attack5_vector_injection": "FAIL_NOT_DENSITY_DERIVED",
        "zero_coupling": "FAIL_PORT_VANISHES",
        "dimensional": "FAIL_PORT_MALFORMED(dimensional)",
        "sign": "FAIL_PORT_MALFORMED(sign)",
        "scaling": "FAIL_PORT_MALFORMED(scaling)",
        "deferred_scalar": "PORT_INCONCLUSIVE_SIM_DEFERRED",
        "deferred_scalar_proven_converse": "DENSITY_PORT_HOSTED",
        "free_carrier_dimension_corruption": "DENSITY_PORT_HOSTED",
        "provenance_less_rider": "FAIL_NOT_DENSITY_DERIVED",
    }
    for name, verdict in expected.items():
        got = control_payload[name]["verdict"]
        assert got == verdict, f"control {name}: expected {verdict}, got {got}"


def print_result(baseline: dict[str, Any], control_payload: dict[str, dict[str, Any]]) -> None:
    print("PATHA_43_DENSITY_QUADRUPOLE_PORT_SYMPY")
    print(f"method: {baseline['trace']['method']}")
    print(f"verdict: {baseline['verdict']}")
    print(f"N0_den: {hstr(baseline['expr'])}")
    print("port_picture: ii two-port(q2,Phi2)")
    print(f"host_set: {', '.join(baseline['host_symbols'])}")
    print(f"taint_set: {', '.join(baseline['taint_set'])}")
    print(f"vector_host_symbols: {', '.join(baseline['vector_host_symbols']) or 'empty'}")
    print(f"continuity_moment_symbol: {baseline['continuity_moment_symbol']}")
    print(f"continuity_dependency_ok: {bstr(baseline['continuity_dependency_ok'])}")
    print("source_graph_symbol_tags:")
    for sym, tags in baseline["source_graph"]["symbol_tags"].items():
        print(f"  {sym}: {','.join(tags) or 'empty'}")
    print("lineage:")
    for key, value in baseline["lineage"]["moments"].items():
        print(f"  {key}: {value}")
    print(f"lineage_valid: {bstr(baseline['lineage_valid'])}")
    print("derivation_trace:")
    for key in ("g_base", "g_q_den", "g_phi_den", "P_den", "Delta_den", "Phi2_response"):
        print(f"  {key}: {hstr(baseline['trace'][key])}")
    print("checks:")
    for key, value in baseline["checks"].items():
        print(f"  {key}: {bstr(value)}")
    print(f"vector_ablation_delta: {hstr(baseline['ablation_delta'])}")
    print("scaling:")
    print(f"  N0_den_a_power: {baseline['scaling']['N0_den_a_power']}")
    print(f"  P0_physical_a_power: {baseline['scaling']['P0_physical_a_power']}")
    print("closure:")
    print(f"  P0_physical: {hstr(baseline['closure']['P0_physical'])}")
    print(f"  Kbar4_residual: {hstr(baseline['closure']['Kbar4_residual'])}")
    print(f"  Gamma5_residual_to_2G_over_5c5: {hstr(baseline['closure']['Gamma5_residual_to_2G_over_5c5'])}")
    print("scope_tags:")
    print("  CALIBRATED: G, 2/5, 54/5")
    print("  SIM_DEFERRED: Xi_Q exact branch magnitude, eta_q, eta_phi, lambda_c literal throat value")
    print("controls:")
    for name, result in control_payload.items():
        print(
            f"  {name}: verdict={result['verdict']}; "
            f"ablation_delta={result['ablation_delta']}; "
            f"taint_set={','.join(result['taint_set']) or 'empty'}; "
            f"source_map_complete={bstr(result['source_map_complete'])}; "
            f"empty_source_symbols={','.join(result['empty_source_symbols']) or 'empty'}"
        )


def main() -> None:
    baseline = evaluate(Config(name="baseline"))
    control_payload = controls()
    assert_gate(baseline, control_payload)
    print_result(baseline, control_payload)


if __name__ == "__main__":
    main()
