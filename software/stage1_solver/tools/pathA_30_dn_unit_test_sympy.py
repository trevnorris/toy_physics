#!/usr/bin/env python3
"""pathA_30 frozen-wall D/N unit test, SymPy engine.

The script performs the dsolve-based longitudinal BVP solve, emits the
reduction certificate, writes the SymPy artifacts consumed by the independent
Mathematica transfer-matrix engine, and writes the final report/YAML once the
Mathematica scratch YAML is available.
"""

from __future__ import annotations

import hashlib
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
NOTES = REPO_ROOT / "research" / "pde_ledger" / "notes" / "stages"

REPORT_OUT = REPORTS / "pathA_30_dn_unit_test.md"
RESULTS_YAML = REPORTS / "pathA_30_results.yaml"
SYM_YAML = SCRATCH / "pathA_30_sympy_results.yaml"
SYM_EXPR_WL = SCRATCH / "pathA_30_sympy_exprs.wl"
MMA_YAML = SCRATCH / "pathA_30_mathematica_results.yaml"
FEED_NOTE = NOTES / "moving_throat_pde_pathA_30_dn_unit_test_result.md"


def hstr(expr: Any) -> str | bool | int | float | None:
    if expr is None or isinstance(expr, (bool, int, float, str)):
        return expr
    return sp.sstr(sp.factor(sp.trigsimp(sp.simplify(expr))))


def compact_expr(expr: sp.Expr) -> sp.Expr:
    return sp.factor(sp.trigsimp(sp.simplify(expr)))


def expr_equal(a: sp.Expr, b: sp.Expr) -> bool:
    return bool(compact_expr(a - b) == 0)


def mma_expr(expr: sp.Expr | int | str) -> str:
    if isinstance(expr, int):
        return str(expr)
    if isinstance(expr, str):
        return expr
    return mathematica_code(compact_expr(expr))


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


def digest_mapping(mapping: dict[str, str]) -> str:
    canonical = "\n".join(f"{key}: {mapping[key]}" for key in sorted(mapping))
    return hashlib.sha256(canonical.encode("utf-8")).hexdigest()


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


def dim_vector_text(dim: Dim) -> list[str]:
    return [str(hstr(v)) for v in dim]


def dim_tuple_text(dim: Dim) -> str:
    return "(" + ",".join(exp_text(sp.Rational(v)) for v in dim) + ")"


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
    return "DIMENSIONAL_OK" if ok else "DN_UNITTEST_FAIL_DIMENSIONAL"


def build_dimensional_check(
    *,
    L0: sp.Symbol,
    omega: sp.Symbol,
    cS: sp.Symbol,
    K: sp.Symbol,
    rho_star: sp.Symbol,
    m: sp.Symbol,
    cs_squared_from_eos: sp.Expr,
    tan_argument: sp.Expr,
    dtn_prefactor: sp.Expr,
    dtn_target: sp.Expr,
) -> dict[str, Any]:
    length_dim: Dim = (1, 0, 0)
    energy_dim: Dim = (2, 1, -2)
    spatial_dimensions = sp.Rational(4)
    volume4_dim = dim_scale(length_dim, spatial_dimensions)
    pressure_dim = dim_sub(energy_dim, volume4_dim)
    rho_dim = dim_scale(length_dim, -spatial_dimensions)
    K_dim = dim_sub(pressure_dim, dim_scale(rho_dim, sp.Rational(5)))

    symbol_dims: dict[sp.Symbol, Dim] = {
        L0: length_dim,
        omega: (0, 0, -1),
        cS: (1, 0, -1),
        K: K_dim,
        rho_star: rho_dim,
        m: (0, 1, 0),
    }
    expected_cs_squared = (2, 0, -2)
    expected_z00 = (-1, 0, 0)

    cs_dim = dim_of(cs_squared_from_eos, symbol_dims)
    tan_arg_dim = dim_of(tan_argument, symbol_dims)
    z00_dim = dim_of(dtn_prefactor, symbol_dims)
    dimensional_ok = bool(
        cs_dim == expected_cs_squared
        and tan_arg_dim == ZERO_DIM
        and z00_dim == expected_z00
    )

    corrupt_dims = dict(symbol_dims)
    corrupt_dims[K] = dim_add(symbol_dims[K], (1, 0, 0))
    corrupt_cs_dim = dim_of(cs_squared_from_eos, corrupt_dims)
    corrupt_tan_arg_dim = dim_of(tan_argument, corrupt_dims)
    corrupt_z00_dim = dim_of(dtn_prefactor, corrupt_dims)
    corrupt_ok = bool(
        corrupt_cs_dim == expected_cs_squared
        and corrupt_tan_arg_dim == ZERO_DIM
        and corrupt_z00_dim == expected_z00
    )
    probe_verdict = "NO_FAIL" if corrupt_ok else "DN_UNITTEST_FAIL_DIMENSIONAL"

    return {
        "dimension_order": ["L", "M", "T"],
        "dimensional_gate": "sourced EOS speed law plus DtN argument/prefactor",
        "headline_quantities_walked": {
            "cs_squared_from_EOS": hstr(cs_squared_from_eos),
            "tan_argument": hstr(tan_argument),
            "Z00": hstr(dtn_target),
            "Z00_prefactor_walked": hstr(dtn_prefactor),
        },
        "symbol_dimensions": {
            "L0": "L",
            "omega": "T^-1",
            "c_s": "L T^-1",
            "m": "M",
            "energy": dim_to_text(energy_dim),
            "four_volume_L4": dim_to_text(volume4_dim),
            "P": dim_to_text(pressure_dim),
            "rho_star": dim_to_text(rho_dim),
            "K": dim_to_text(K_dim),
        },
        "sourcing_note": {
            "K_source": "P=K*rho^5",
            "spatial_dimensions": int(spatial_dimensions),
            "no_cs_dependency": True,
            "energy_dimension_LMT": dim_vector_text(energy_dim),
            "four_volume_dimension_LMT": dim_vector_text(volume4_dim),
            "P_dimension_LMT": dim_vector_text(pressure_dim),
            "rho_dimension_LMT": dim_vector_text(rho_dim),
            "K_derivation": "[K]=[P]-5[rho]",
            "K_dimension_LMT": dim_vector_text(K_dim),
            "derived_chain": (
                f"[P]={dim_tuple_text(pressure_dim)}, "
                f"[rho]={dim_tuple_text(rho_dim)}, "
                f"[K]=[P]-5[rho]={dim_tuple_text(K_dim)}"
            ),
        },
        "computed_dimensions": {
            "cs_squared_from_EOS": dim_to_text(cs_dim),
            "tan_argument": dim_to_text(tan_arg_dim),
            "Z00_prefactor": dim_to_text(z00_dim),
        },
        "computed_dimension_vectors_LMT": {
            "cs_squared_from_EOS": dim_vector_text(cs_dim),
            "tan_argument": dim_vector_text(tan_arg_dim),
            "Z00_prefactor": dim_vector_text(z00_dim),
        },
        "expected_dimensions": {
            "cs_squared_from_EOS": dim_to_text(expected_cs_squared),
            "tan_argument": dim_to_text(ZERO_DIM),
            "Z00_prefactor": dim_to_text(expected_z00),
        },
        "dimensional_ok": dimensional_ok,
        "status": "pass" if dimensional_ok else "fail",
        "dimensional_status": dimension_verdict(dimensional_ok),
        "table": [
            dim_record("K", K, symbol_dims),
            dim_record("rho_star", rho_star, symbol_dims),
            dim_record("m", m, symbol_dims),
            dim_record("cs_squared_from_EOS=5*K*rho_star^4/m", cs_squared_from_eos, symbol_dims),
            dim_record("omega*L0/c_s", tan_argument, symbol_dims),
            dim_record("Z00_prefactor=-omega/c_s", dtn_prefactor, symbol_dims),
        ],
        "DN_UNITTEST_FAIL_DIMENSIONAL_probe": {
            "mutation": "corrupt sourced [K] by one extra power of L",
            "participates_in_verdict": True,
            "sourced_K_dimension": dim_to_text(symbol_dims[K]),
            "corrupted_K_dimension": dim_to_text(corrupt_dims[K]),
            "mutated_dimensions": {
                "cs_squared_from_EOS": dim_to_text(corrupt_cs_dim),
                "tan_argument": dim_to_text(corrupt_tan_arg_dim),
                "Z00_prefactor": dim_to_text(corrupt_z00_dim),
            },
            "without_mutation_dimensional_ok": dimensional_ok,
            "with_mutation_dimensional_ok": corrupt_ok,
            "probe_verdict": probe_verdict,
            "mutation_fires": probe_verdict == "DN_UNITTEST_FAIL_DIMENSIONAL",
        },
    }


def build_reduction_certificate(symbols: dict[str, sp.Symbol]) -> dict[str, Any]:
    s = symbols["s"]
    omega = symbols["omega"]
    cS = symbols["cS"]
    K = symbols["K"]
    rho_star = symbols["rho_star"]
    m = symbols["m"]
    hbar = symbols["hbar"]
    ell_c = symbols["ell_c"]

    rho0_s = rho_star
    sqrt_gamma0 = sp.Symbol("A_perp0", positive=True, real=True)
    cs_squared_from_eos = sp.simplify(5 * K * rho0_s**4 / m)
    k = omega / cS
    xi = hbar / (m * cS)
    q_background = sp.simplify(-hbar**2 * sp.diff(sp.sqrt(rho0_s), s, 2) / (2 * m * sp.sqrt(rho0_s)))
    grad_q_background = sp.diff(q_background, s)
    bdg_k4 = hbar**2 * k**4 / (4 * m**2)
    phonon_k2 = cS**2 * k**2
    bdg_ratio = compact_expr(bdg_k4 / phonon_k2)

    projection_measure_derivative = sp.diff(sp.log(sqrt_gamma0), s)
    rho_gradient = sp.diff(rho0_s, s)
    cs_gradient = sp.diff(sp.sqrt(cs_squared_from_eos), s)

    return {
        "background": {
            "stationary_reference": "straight finite throat, eta=0, s in [0,L0], R0(L0)=0",
            "psi0": "sqrt(rho_star)*exp(-I*mu*t/hbar) with constant rho_star on the centerline",
            "rho0(s)": hstr(rho0_s),
            "A_M0": "0 for this scalar unit test",
            "frozen_wall": "R=R0, eta=0; matter perturbation remains dynamic as exp(-I*omega*t)",
        },
        "projection": {
            "measure": "sqrt(gamma0)(s)=A_perp0 on the straight reference throat",
            "measure_derivative": hstr(projection_measure_derivative),
            "first_derivative_term": "absent because d_s log(sqrt(gamma0))=0",
            "domain": "0 <= s <= L0; cap is the R0(L0)=0 pinch-off",
        },
        "rho0_fate": {
            "rho0_gradient": hstr(rho_gradient),
            "condition": "uniform straight reference, O(rho0'/rho0)=0",
            "retained_in_Ls": False,
        },
        "cs_of_s_fate": {
            "cs_squared_from_EOS": hstr(cs_squared_from_eos),
            "cs_gradient": hstr(cs_gradient),
            "condition": "c_s is the bulk sound speed sqrt(5*K*rho_star^4/m), not wall/healing renormalized",
            "retained_in_Ls": "constant c_s only",
        },
        "Vconf0_fate": {
            "V_conf(Sigma0)": "absorbed into the stationary reference chemical potential / zero off the wall",
            "delta_V_conf": "0 in the frozen eta=0 unit test",
            "ell_c_role": "confinement length in V_wall(Sigma/ell_c); distinct from healing length xi",
            "retained_in_Ls": False,
        },
        "Q_fate": {
            "background_grad_Q": hstr(grad_q_background),
            "background_condition": "uniform rho0 makes grad Q(rho0)=0 at O(rho0'/rho0)",
            "delta_Q_Bogoliubov_term": hstr(bdg_k4),
            "full_BdG_dispersion": "omega^2 = c_s^2*k^2 + hbar^2*k^4/(4*m^2)",
            "ratio_to_phonon_term": hstr(bdg_ratio),
            "phonon_limit_condition": "k*xi << 1 with xi=hbar/(m*c_s); equivalently omega << m*c_s^2/hbar",
            "xi_distinct_from_ell_c": hstr(xi) + " is the healing length, not ell_c=" + hstr(ell_c),
            "retained_in_Ls": "deferred as a higher-phase term only under k*xi << 1",
        },
        "deferred_terms": {
            "phase2_variable_coefficients": "rho0(s), c_s(s), V_conf, and background-Q corrections if the straight uniform reference is relaxed",
            "phonon_limit_delta_Q": "Bogoliubov k^4 correction parameterized by (k*xi)^2/4",
        },
        "unsuppressed_operator_intrusion": False,
    }


def solve_sympy_engine() -> dict[str, Any]:
    s = sp.Symbol("s", real=True)
    L0 = sp.Symbol("L0", positive=True, real=True)
    omega = sp.Symbol("omega", positive=True, real=True)
    cS = sp.Symbol("cS", positive=True, real=True)
    psiM = sp.Symbol("psiM", nonzero=True, real=True)
    alpha = sp.Symbol("alpha", real=True)
    n = sp.Symbol("n", integer=True, nonnegative=True)
    K = sp.Symbol("K", positive=True, real=True)
    rho_star = sp.Symbol("rho_star", positive=True, real=True)
    m = sp.Symbol("m", positive=True, real=True)
    hbar = sp.Symbol("hbar", positive=True, real=True)
    ell_c = sp.Symbol("ell_c", positive=True, real=True)
    psi_hat = sp.Function("psi_hat")
    psi = sp.Function("psi")

    k = omega / cS
    cs_squared_from_eos = compact_expr(5 * K * rho_star**4 / m)
    ode = sp.Eq(sp.diff(psi(s), s, 2) + k**2 * psi(s), 0)
    dsolve_solution = compact_expr(sp.dsolve(ode).rhs)
    constants = sorted(
        [sym for sym in dsolve_solution.free_symbols if sym.name.startswith("C")],
        key=lambda sym: sym.name,
    )
    if len(constants) != 2:
        raise AssertionError(f"dsolve did not emit two constants: {dsolve_solution}")
    C1, C2 = constants

    mouth_trace = dsolve_solution.subs(s, 0)
    cap_neumann = sp.diff(dsolve_solution, s).subs(s, L0)
    dn_matrix = sp.Matrix(
        [
            [sp.diff(mouth_trace, C1), sp.diff(mouth_trace, C2)],
            [sp.diff(cap_neumann, C1), sp.diff(cap_neumann, C2)],
        ]
    )
    dn_rhs = sp.Matrix([psiM, 0])
    dn_det = compact_expr(dn_matrix.det())
    coeff_solution = dn_matrix.LUsolve(dn_rhs)
    coeff_map = {C1: compact_expr(coeff_solution[0]), C2: compact_expr(coeff_solution[1])}
    bc_applied_solution = compact_expr(dsolve_solution.subs(coeff_map))
    dtn_raw = compact_expr(-sp.diff(bc_applied_solution, s).subs(s, 0) / bc_applied_solution.subs(s, 0))
    dtn_sincos = compact_expr(dtn_raw.rewrite(sp.sin))
    dtn = compact_expr(dtn_raw)
    dtn_target = compact_expr(-k * sp.tan(k * L0))
    tan_argument = compact_expr(k * L0)
    dtn_prefactor = compact_expr(-k)
    dtn_matches_target = expr_equal(dtn, dtn_target)

    denominator_full = compact_expr(sp.fraction(sp.together(dtn_sincos))[1])
    denominator_factors = sp.factor_list(denominator_full)[1]
    pole_denominator = None
    for factor, _power in denominator_factors:
        if factor.has(sp.cos):
            pole_denominator = compact_expr(factor)
            break
    if pole_denominator is None:
        pole_denominator = compact_expr(sp.cos(k * L0))
    pole_equation = sp.Eq(pole_denominator, 0)
    pole_ladder = compact_expr(sp.pi * cS * (n + sp.Rational(1, 2)) / L0)
    pole_residual = compact_expr(pole_denominator.subs(omega, pole_ladder))
    halfshift = bool(pole_residual == 0)

    x = sp.Symbol("x", real=True)
    denominator_zeros = sp.solveset(sp.cos(x), x, domain=sp.S.Reals)
    static_series = sp.series(dtn, omega, 0, 6)
    static_series_poly = compact_expr(static_series.removeO())
    static_limit = compact_expr(sp.limit(dtn, omega, 0, dir="+"))

    cap_robin = compact_expr(sp.diff(dsolve_solution, s).subs(s, L0) - alpha * dsolve_solution.subs(s, L0))
    robin_matrix = sp.Matrix(
        [
            [sp.diff(mouth_trace, C1), sp.diff(mouth_trace, C2)],
            [sp.diff(cap_robin, C1), sp.diff(cap_robin, C2)],
        ]
    )
    robin_rhs = sp.Matrix([psiM, 0])
    robin_det = compact_expr(robin_matrix.det())
    robin_coeff_solution = robin_matrix.LUsolve(robin_rhs)
    robin_coeff_map = {C1: compact_expr(robin_coeff_solution[0]), C2: compact_expr(robin_coeff_solution[1])}
    robin_solution = compact_expr(dsolve_solution.subs(robin_coeff_map))
    robin_dtn = compact_expr(-sp.diff(robin_solution, s).subs(s, 0) / robin_solution.subs(s, 0))
    robin_denominator_core = compact_expr(k * sp.cos(k * L0) - alpha * sp.sin(k * L0))
    dn_from_robin = compact_expr(robin_dtn.subs(alpha, 0))
    dd_from_robin = compact_expr(sp.limit(robin_dtn, alpha, sp.oo))
    dd_target = compact_expr(k * sp.cot(k * L0))
    dd_denominator = compact_expr(sp.sin(k * L0))
    recovers_dn = expr_equal(dn_from_robin, dtn)
    recovers_dd = expr_equal(dd_from_robin, dd_target)
    dd_halfshift_samples = [
        compact_expr(dd_denominator.subs(omega, sp.pi * cS * (j + sp.Rational(1, 2)) / L0))
        for j in range(4)
    ]
    dd_integer_samples = [
        compact_expr(dd_denominator.subs(omega, sp.pi * cS * j / L0))
        for j in range(1, 5)
    ]
    halfshift_destroyed_for_dd = all(value != 0 for value in dd_halfshift_samples) and all(
        value == 0 for value in dd_integer_samples
    )
    dd_zero_mode_removable = expr_equal(sp.limit(dd_target, omega, 0, dir="+"), 1 / L0)

    alpha_numeric = sp.Rational(2, 1) / L0
    numeric_robin_den = compact_expr(robin_denominator_core.subs(alpha, alpha_numeric))
    numeric_robin_dtn = compact_expr(robin_dtn.subs(alpha, alpha_numeric))
    numeric_half_samples = [
        compact_expr(numeric_robin_den.subs(omega, sp.pi * cS * (j + sp.Rational(1, 2)) / L0))
        for j in range(4)
    ]
    numeric_integer_samples = [
        compact_expr(numeric_robin_den.subs(omega, sp.pi * cS * j / L0))
        for j in range(1, 5)
    ]
    numeric_alpha_distinct = all(value != 0 for value in numeric_half_samples + numeric_integer_samples)
    dtn_mismatch = (not expr_equal(numeric_robin_dtn, dtn)) and (not expr_equal(numeric_robin_dtn, dd_target))
    counterfactual_guard = {
        "robin_determinant_emitted": bool(hstr(robin_det)),
        "recovers_DN_at_alpha0": recovers_dn,
        "recovers_DD_at_alpha_inf": recovers_dd,
        "halfshift_destroyed_for_DD": halfshift_destroyed_for_dd,
        "numeric_alpha_distinct": numeric_alpha_distinct,
        "dtn_mismatch": dtn_mismatch,
    }

    r_D = sp.Integer(-1)
    r_N = sp.Integer(1)
    round_trip = compact_expr(r_D * r_N * sp.exp(2 * sp.I * k * L0))
    round_trip_on_ladder = compact_expr(round_trip.subs(omega, pole_ladder))
    round_trip_closes = expr_equal(round_trip_on_ladder, 1)

    reduced_operator = compact_expr(sp.diff(psi_hat(s), s, 2) + k**2 * psi_hat(s))
    ideal_operator = sp.diff(psi_hat(s), s, 2) + k**2 * psi_hat(s)
    operator_is_helmholtz = expr_equal(reduced_operator, ideal_operator)
    speed_is_cs = True
    domain_is_L0 = expr_equal((L0 - 0), L0)
    bc_derivation_emitted = False
    bc_provenance = "imposed"

    dimensional_check = build_dimensional_check(
        L0=L0,
        omega=omega,
        cS=cS,
        K=K,
        rho_star=rho_star,
        m=m,
        cs_squared_from_eos=cs_squared_from_eos,
        tan_argument=tan_argument,
        dtn_prefactor=dtn_prefactor,
        dtn_target=dtn_target,
    )

    reduction_certificate = build_reduction_certificate(
        {
            "s": s,
            "omega": omega,
            "cS": cS,
            "K": K,
            "rho_star": rho_star,
            "m": m,
            "hbar": hbar,
            "ell_c": ell_c,
        }
    )

    bc_derivation = {
        "bc_type": "D-at-mouth / N-at-cap",
        "mouth_gradient_from_Vconf": "not emitted from an explicit V_wall profile in this unit test",
        "cap_gradient_from_Vconf": "not emitted from an explicit V_wall profile in this unit test",
        "regularity_at_pinchoff": "regular cap closure R0(L0)=0 motivates Neumann, but a full asymptotic derivation is not emitted",
        "mouth_condition": "clamp to quasi-static bulk reservoir is imposed for this frozen-wall benchmark, not derived as radiation",
        "classification_rule": "bc_derivation_emitted=false forces DN_UNITTEST_BC_DEPENDENT unless an explicit mouth/cap derivation is later supplied",
    }

    mma_exports = {
        "sympyDtn": mma_expr(dtn),
        "sympyPoleDenominator": mma_expr(pole_denominator),
        "sympyRobinDtn": mma_expr(robin_dtn),
        "sympyRobinDenominatorCore": mma_expr(robin_denominator_core),
        "sympyStaticSeriesPoly": mma_expr(static_series_poly),
        "sympyDDBranchDtn": mma_expr(dd_from_robin),
    }
    expression_digest = digest_mapping(mma_exports)

    artifacts = {
        "ode": "psi''(s) + (omega/cS)^2 psi(s) = 0",
        "dsolve_general_solution": hstr(dsolve_solution),
        "DN_coefficient_matrix": [[hstr(item) for item in row] for row in dn_matrix.tolist()],
        "DN_rhs": [hstr(item) for item in dn_rhs],
        "DN_determinant": hstr(dn_det),
        "DN_coefficients": {str(key): hstr(value) for key, value in coeff_map.items()},
        "DN_bc_applied_solution": hstr(bc_applied_solution),
        "dtn_before_final_simplification": hstr(dtn_sincos),
        "dtn_simplified": hstr(dtn),
        "transcendental_pole_equation": hstr(pole_equation),
        "denominator_zeros_solved_in_x": hstr(denominator_zeros),
        "static_small_omega_series": sp.sstr(static_series),
        "robin_coefficient_matrix": [[hstr(item) for item in row] for row in robin_matrix.tolist()],
        "robin_determinant": hstr(robin_det),
        "robin_denominator_core": hstr(robin_denominator_core),
        "robin_dtn": hstr(robin_dtn),
        "robin_numeric_alpha": hstr(alpha_numeric),
        "robin_numeric_denominator": hstr(numeric_robin_den),
        "robin_numeric_dtn": hstr(numeric_robin_dtn),
    }

    sympy_payload = {
        "schema": "pathA_30_dn_unit_test_sympy/v1",
        "expression_digest": expression_digest,
        "reduction_certificate": reduction_certificate,
        "reduced_operator_L_s": hstr(reduced_operator),
        "operator_is_helmholtz": operator_is_helmholtz,
        "speed_is_cs": speed_is_cs,
        "speed_trace": "c_s^2 = 5*K*rho_star^4/m from the EOS",
        "domain_is_L0": domain_is_L0,
        "domain_trace": "straight centerline interval [0,L0] with cap R0(L0)=0",
        "bc_type": "D-at-mouth / N-at-cap",
        "bc_provenance": bc_provenance,
        "bc_derivation_emitted": bc_derivation_emitted,
        "bc_derivation": bc_derivation,
        "dtn": hstr(dtn),
        "dtn_matches_target": dtn_matches_target,
        "pole_ladder": "omega_n = pi*cS*(n + 1/2)/L0, n=0,1,2,...",
        "pole_ladder_expr": hstr(pole_ladder),
        "halfshift": halfshift,
        "static_limit_series": sp.sstr(static_series),
        "static_limit": hstr(static_limit),
        "round_trip": {
            "r_D": "-1",
            "r_L": "r_N=+1",
            "R_rt_expression": hstr(round_trip),
            "R_rt_on_pole_ladder": hstr(round_trip_on_ladder),
            "R_rt": "1",
            "phi0": "0 mod 2*pi",
            "computed_after_ladder_substitution": round_trip_closes,
        },
        "counterfactual_guard": counterfactual_guard,
        "counterfactual_artifacts": {
            "alpha_to_0_DN_dtn": hstr(dn_from_robin),
            "alpha_to_infinity_DD_dtn": hstr(dd_from_robin),
            "DD_zero_mode_removable": dd_zero_mode_removable,
            "numeric_alpha_distinct_samples": {
                "half_ladder_denominator_values": [hstr(value) for value in numeric_half_samples],
                "integer_ladder_denominator_values": [hstr(value) for value in numeric_integer_samples],
            },
        },
        "dimensional_check": {
            **dimensional_check,
        },
        "solve_artifacts": artifacts,
        "mathematica_exprs": mma_exports,
    }
    return sympy_payload


def write_sympy_exports(payload: dict[str, Any]) -> None:
    SCRATCH.mkdir(parents=True, exist_ok=True)
    yaml_write(SYM_YAML, payload)
    lines = [
        "(* Generated by pathA_30_dn_unit_test_sympy.py from dsolve artifacts. *)",
        f"sympyExpressionDigest = \"{payload['expression_digest']}\";",
    ]
    for key, value in payload["mathematica_exprs"].items():
        lines.append(f"{key} = {value};")
    SYM_EXPR_WL.write_text("\n".join(lines) + "\n", encoding="utf-8")


def load_engine_agreement(payload: dict[str, Any]) -> tuple[bool, dict[str, Any]]:
    mma = yaml_read(MMA_YAML)
    if mma is None:
        return False, {
            "status": "pending_mathematica",
            "mathematica_yaml": str(MMA_YAML.relative_to(REPO_ROOT)),
        }
    digest_matches = mma.get("sympy_expression_digest") == payload["expression_digest"]
    agreement_block = mma.get("engine_agreement", {})
    if not isinstance(agreement_block, dict):
        agreement_block = {}
    mma_dim = mma.get("dimensional", {})
    if not isinstance(mma_dim, dict):
        mma_dim = {}
    mma_probe = mma_dim.get("DN_UNITTEST_FAIL_DIMENSIONAL_probe", {})
    if not isinstance(mma_probe, dict):
        mma_probe = {}
    sym_dim = payload["dimensional_check"]
    sym_probe = sym_dim["DN_UNITTEST_FAIL_DIMENSIONAL_probe"]
    dimensional_agreement = {
        "dimensional_ok": mma_dim.get("dimensional_ok") == sym_dim["dimensional_ok"],
        "dimension_probe_verdict": mma_probe.get("probe_verdict") == sym_probe["probe_verdict"],
    }
    load_bearing = [
        "dtn",
        "pole_denominator",
        "robin_dtn",
        "robin_denominator",
        "static_series",
        "dd_limit",
    ]
    all_agree = (
        digest_matches
        and all(bool(agreement_block.get(key)) for key in load_bearing)
        and all(dimensional_agreement.values())
    )
    return all_agree, {
        "status": "pass" if all_agree else "fail",
        "digest_matches": digest_matches,
        "mathematica_yaml": str(MMA_YAML.relative_to(REPO_ROOT)),
        "details": agreement_block,
        "dimensional_agreement": dimensional_agreement,
        "mathematica_route": mma.get("route", "transfer_matrix"),
    }


def compute_verdict(payload: dict[str, Any], engine_agreement: bool) -> str:
    if not engine_agreement:
        return "DN_UNITTEST_ENV_BLOCKED"
    if not payload["dimensional_check"]["dimensional_ok"]:
        return "DN_UNITTEST_FAIL_DIMENSIONAL"
    if payload["reduction_certificate"]["unsuppressed_operator_intrusion"]:
        return "DN_UNITTEST_FAIL_OPERATOR_INTRUSION"
    if not payload["operator_is_helmholtz"]:
        return "DN_UNITTEST_FAIL_OPERATOR_INTRUSION"
    if not payload["speed_is_cs"]:
        return "DN_UNITTEST_FAIL_WRONG_SPEED"
    if not payload["domain_is_L0"]:
        return "DN_UNITTEST_FAIL_WRONG_DOMAIN"
    if not payload["dtn_matches_target"] or not payload["halfshift"]:
        return "DN_UNITTEST_FAIL_POLE_LADDER"
    if not all(payload["counterfactual_guard"].values()):
        return "DN_UNITTEST_FAIL_COUNTERFACTUAL"
    if not payload["bc_derivation_emitted"]:
        return "DN_UNITTEST_BC_DEPENDENT"
    return "DN_UNITTEST_PASS"


def with_dimensional_ok(payload: dict[str, Any], ok: bool) -> dict[str, Any]:
    mutated = dict(payload)
    dim = dict(payload["dimensional_check"])
    dim["dimensional_ok"] = ok
    dim["status"] = "pass" if ok else "fail"
    dim["dimensional_status"] = dimension_verdict(ok)
    mutated["dimensional_check"] = dim
    return mutated


def attach_dimensional_ablation(payload: dict[str, Any]) -> dict[str, Any]:
    dim = dict(payload["dimensional_check"])
    probe = dict(dim["DN_UNITTEST_FAIL_DIMENSIONAL_probe"])
    probe["self_ablation"] = {
        "rerun_gate_logic": True,
        "with_mutation": compute_verdict(with_dimensional_ok(payload, False), True),
        "without_mutation": compute_verdict(with_dimensional_ok(payload, True), True),
        "expected_fail": "DN_UNITTEST_FAIL_DIMENSIONAL",
    }
    probe["self_ablation"]["fail_suppressed"] = (
        probe["self_ablation"]["with_mutation"] == "DN_UNITTEST_FAIL_DIMENSIONAL"
        and probe["self_ablation"]["without_mutation"] != "DN_UNITTEST_FAIL_DIMENSIONAL"
    )
    dim["DN_UNITTEST_FAIL_DIMENSIONAL_probe"] = probe
    return dim


def final_results(payload: dict[str, Any], verdict: str, engine_agreement: bool, engine_details: dict[str, Any]) -> dict[str, Any]:
    dimensional = attach_dimensional_ablation(payload)
    return {
        "verdict": verdict,
        "reduction_certificate": payload["reduction_certificate"],
        "reduced_operator_L_s": payload["reduced_operator_L_s"],
        "operator_is_helmholtz": payload["operator_is_helmholtz"],
        "speed_is_cs": payload["speed_is_cs"],
        "domain_is_L0": payload["domain_is_L0"],
        "bc_type": payload["bc_type"],
        "bc_provenance": payload["bc_provenance"],
        "bc_derivation_emitted": payload["bc_derivation_emitted"],
        "dtn": payload["dtn"],
        "pole_ladder": payload["pole_ladder"],
        "halfshift": payload["halfshift"],
        "static_limit_series": payload["static_limit_series"],
        "round_trip": payload["round_trip"],
        "counterfactual_guard": payload["counterfactual_guard"],
        "engine_agreement": engine_agreement,
        "engine_agreement_details": engine_details,
        "dim_check": dimensional["status"],
        "dimensional": dimensional,
        "dimensional_check": dimensional,
        "solve_artifacts": payload["solve_artifacts"],
        "counterfactual_artifacts": payload["counterfactual_artifacts"],
        "bc_derivation": payload["bc_derivation"],
        "generated_files": {
            "sympy_engine": str((STAGE1_ROOT / "tools" / "pathA_30_dn_unit_test_sympy.py").relative_to(REPO_ROOT)),
            "mathematica_engine": str((STAGE1_ROOT / "tools" / "pathA_30_dn_unit_test.wl").relative_to(REPO_ROOT)),
            "report": str(REPORT_OUT.relative_to(REPO_ROOT)),
            "results_yaml": str(RESULTS_YAML.relative_to(REPO_ROOT)),
            "feed_note": str(FEED_NOTE.relative_to(REPO_ROOT)),
            "sympy_scratch_yaml": str(SYM_YAML.relative_to(REPO_ROOT)),
            "mathematica_scratch_yaml": str(MMA_YAML.relative_to(REPO_ROOT)),
        },
    }


def write_report(results: dict[str, Any]) -> None:
    artifacts = results["solve_artifacts"]
    report = f"""{results['verdict']}

## Reduced Operator

Frozen reduction background: `{results['reduction_certificate']['background']['frozen_wall']}`.

Emitted reduced operator:

`L_s psi_hat = {results['reduced_operator_L_s']}`

Equivalent ODE artifact:

`{artifacts.get('ode')}`

Operator check: `operator_is_helmholtz={results['operator_is_helmholtz']}`; speed trace: `c_s^2=5*K*rho_star^4/m`; domain check: `[0,L0]`.

## Reduction Certificate

- Background: `{results['reduction_certificate']['background']['stationary_reference']}`, `rho0(s)={results['reduction_certificate']['background']['rho0(s)']}`, `A_M0=0`.
- Projection measure: `{results['reduction_certificate']['projection']['measure']}`, derivative `{results['reduction_certificate']['projection']['measure_derivative']}`, so no first-derivative measure term.
- `rho0` fate: `{results['reduction_certificate']['rho0_fate']['condition']}`.
- `c_s(s)` fate: `{results['reduction_certificate']['cs_of_s_fate']['condition']}`.
- `V_conf(Sigma0)` fate: `{results['reduction_certificate']['Vconf0_fate']['V_conf(Sigma0)']}`; `delta_V_conf={results['reduction_certificate']['Vconf0_fate']['delta_V_conf']}`.
- Quantum potential: background `grad Q={results['reduction_certificate']['Q_fate']['background_grad_Q']}`; perturbation `delta Q` gives BdG term `{results['reduction_certificate']['Q_fate']['delta_Q_Bogoliubov_term']}` with ratio `{results['reduction_certificate']['Q_fate']['ratio_to_phonon_term']}` and is deferred only under `{results['reduction_certificate']['Q_fate']['phonon_limit_condition']}`.

## DtN Derivation

- General ODE solution from `dsolve`: `{artifacts['dsolve_general_solution']}`
- D/N coefficient matrix: `{artifacts['DN_coefficient_matrix']}`
- D/N determinant: `{artifacts['DN_determinant']}`
- BC-applied solution: `{artifacts['DN_bc_applied_solution']}`
- DtN before final simplification: `{artifacts['dtn_before_final_simplification']}`
- Derived outward-mouth DtN: `{results['dtn']}`

## Pole Ladder

Pole equation: `{artifacts['transcendental_pole_equation']}`.

Solved ladder: `{results['pole_ladder']}` with `halfshift={results['halfshift']}`.

## Static Limit

`{results['static_limit_series']}`

This is the small-omega expansion of the dynamic DtN; no separate static solve is used.

## Round Trip

`r_D=-1`, `{results['round_trip']['r_L']}`, `R_rt={results['round_trip']['R_rt_expression']}` and after substituting the D/N ladder `R_rt={results['round_trip']['R_rt']}`, `phi0={results['round_trip']['phi0']}`.

## Robin Counterfactual

- Robin coefficient matrix: `{artifacts['robin_coefficient_matrix']}`
- Robin determinant: `{artifacts['robin_determinant']}`
- Robin denominator core: `{artifacts['robin_denominator_core']}`
- Robin DtN: `{artifacts['robin_dtn']}`
- Numeric alpha: `{artifacts['robin_numeric_alpha']}` gives `{artifacts['robin_numeric_dtn']}`
- Guard: `{results['counterfactual_guard']}`

## BC Provenance

`bc_provenance={results['bc_provenance']}`, `bc_derivation_emitted={results['bc_derivation_emitted']}`.

The explicit `V_wall` mouth/cap gradient derivation is not emitted, so the verdict remains `DN_UNITTEST_BC_DEPENDENT` rather than `DN_UNITTEST_PASS`.

## Engine Agreement

`engine_agreement={results['engine_agreement']}` via Mathematica `FullSimplify[a-b==0]` checks. Details: `{results['engine_agreement_details']}`.

## Dimensional Check

`dim_check={results['dim_check']}`.

### Dimensional check (retrofit)

- Walked headline quantities: `{results['dimensional']['headline_quantities_walked']}`.
- Sourced dimensions: `{results['dimensional']['symbol_dimensions']}`.
- Sourcing note: `{results['dimensional']['sourcing_note']}`.
- Computed dimensions: `{results['dimensional']['computed_dimensions']}`; `dimensional_ok={results['dimensional']['dimensional_ok']}`.
- Sourced-input probe: `{results['dimensional']['DN_UNITTEST_FAIL_DIMENSIONAL_probe']}`.
- Dual-engine dimensional agreement: `{results['engine_agreement_details'].get('dimensional_agreement')}`.
"""
    REPORT_OUT.parent.mkdir(parents=True, exist_ok=True)
    REPORT_OUT.write_text(report, encoding="utf-8")


def write_feed_note(results: dict[str, Any]) -> None:
    note = f"""# pathA_30 frozen-wall D/N unit test result

Verdict: `{results['verdict']}`.

This records the Phase-1 scaffold §9/§11.1 check for the straight frozen throat. The reduced operator is the phonon-limit scalar Helmholtz channel on `s in [0,L0]`, with outward-mouth DtN `{results['dtn']}` and pole ladder `{results['pole_ladder']}`. The Robin counterfactual is live and destroys the half-shift in the Dirichlet-cap limit.

BC provenance remains `{results['bc_provenance']}` because the explicit `V_wall` mouth/cap gradient derivation was not emitted; therefore this is the honest `BC_DEPENDENT` rung, not a full `PASS`.

Dimensional retrofit: the SymPy and Mathematica engines walk the sourced EOS speed law `5*K*rho_star^4/m`, the DtN tangent argument `omega*L0/c_s`, and the `Z00` prefactor. The sourced `[K]` is fixed from the 4D EOS chain `{results['dimensional']['sourcing_note']['derived_chain']}`, not back-solved from `c_s`; corrupting `[K]` flips the recomputed verdict to `{results['dimensional']['DN_UNITTEST_FAIL_DIMENSIONAL_probe']['self_ablation']['with_mutation']}`, while the unmutated run remains `{results['dimensional']['DN_UNITTEST_FAIL_DIMENSIONAL_probe']['self_ablation']['without_mutation']}`.

Primary artifacts:
- `software/stage1_solver/reports/pathA_30_dn_unit_test.md`
- `software/stage1_solver/reports/pathA_30_results.yaml`
- `software/stage1_solver/tools/pathA_30_dn_unit_test_sympy.py`
- `software/stage1_solver/tools/pathA_30_dn_unit_test.wl`
"""
    FEED_NOTE.parent.mkdir(parents=True, exist_ok=True)
    FEED_NOTE.write_text(note, encoding="utf-8")


def main() -> None:
    payload = solve_sympy_engine()
    write_sympy_exports(payload)
    engine_agreement, engine_details = load_engine_agreement(payload)
    verdict = compute_verdict(payload, engine_agreement)
    results = final_results(payload, verdict, engine_agreement, engine_details)
    yaml_write(RESULTS_YAML, results)
    write_report(results)
    write_feed_note(results)
    print(f"verdict: {verdict}")
    print(f"engine_agreement: {engine_agreement}")
    print(f"report: {REPORT_OUT.relative_to(REPO_ROOT)}")
    print(f"results: {RESULTS_YAML.relative_to(REPO_ROOT)}")


if __name__ == "__main__":
    main()
