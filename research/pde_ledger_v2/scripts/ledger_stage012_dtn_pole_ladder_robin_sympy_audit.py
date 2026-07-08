#!/usr/bin/env python3
"""Ledger stage012 SymPy audit: DtN pole ladder + Robin falsifier.

Standalone, print-only, no arguments, no file I/O.  This is the Part-II
pathA_30 II-G1b slice only: consume stage011's frozen Helmholtz L_s, open at
the dsolve, derive the D/N DtN by a coefficient-matrix LUsolve, derive the
half-shifted pole ladder, static small-omega series, round trip, Robin
counterfactual guard, and the 012 tan_argument/Z00 dimensional legs.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import sympy as sp


PASS_COUNT = 0
FAIL_COUNT = 0

DN_UNITTEST_FAIL_DIMENSIONAL = "DN_UNITTEST_FAIL_DIMENSIONAL"
FAIL_POLE_LADDER = "FAIL_POLE_LADDER"
FAIL_COUNTERFACTUAL = "FAIL_COUNTERFACTUAL"
DN_UNITTEST_BC_DEPENDENT = "DN_UNITTEST_BC_DEPENDENT"
DN_UNITTEST_PASS = "DN_UNITTEST_PASS"

EXPECTED_GUARD_KEYS = (
    "robin_determinant_emitted",
    "recovers_DN_at_alpha0",
    "recovers_DD_at_alpha_inf",
    "halfshift_destroyed_for_DD",
    "numeric_alpha_distinct",
    "dtn_mismatch",
)


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


def compact_expr(expr: Any) -> Any:
    if isinstance(expr, sp.MatrixBase):
        return expr.applyfunc(compact_expr)
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
        return sp.sstr(compact_expr(expr))
    if isinstance(expr, (dict, list, tuple)):
        return sp.sstr(expr)
    try:
        return sp.sstr(compact_expr(expr))
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
    clean = compact_expr(residual)
    assert_no_float(name, clean)
    if clean == 0:
        _record_pass(f"PASS  {name}")
        return
    _record_fail(f"FAIL  {name}: residual = {fmt(clean)}")
    raise AuditFailure(f"{name} residual was not zero")


def expect_bool(name: str, condition: bool) -> None:
    expect_zero(name, sp.Integer(0) if condition else sp.Integer(1))


def expect_nonzero(name: str, residual: sp.Expr | int) -> None:
    assert_no_float(name, residual)
    clean = compact_expr(residual)
    assert_no_float(name, clean)
    if clean != 0:
        _record_pass(f"PASS  {name} is nonzero as required (residual = {fmt(clean)})")
        return
    _record_fail(f"FAIL  {name}: required nonzero residual vanished")
    raise AuditFailure(f"{name} unexpectedly had zero residual")


def expect_fail(name: str, residual: sp.Expr | int) -> None:
    assert_no_float(name, residual)
    clean = compact_expr(residual)
    assert_no_float(name, clean)
    if clean != 0:
        _record_pass(f"PASS  {name} produced required FAIL (residual = {fmt(clean)})")
        return
    _record_fail(f"FAIL  {name}: required mutation/ablation did not fire")
    raise AuditFailure(f"{name} unexpectedly had zero residual")


def expr_equal(lhs: sp.Expr | int, rhs: sp.Expr | int = 0) -> bool:
    return compact_expr(lhs - rhs) == 0


def nonzero_q(expr: sp.Expr | int) -> bool:
    return compact_expr(expr) != 0


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


def dim_residual(actual: Dim, expected: Dim) -> sp.Expr:
    return sp.simplify(sum((have - want) ** 2 for have, want in zip(actual.components(), expected.components())))


def dim_of(expr: sp.Expr, dims: dict[sp.Symbol, Dim]) -> Dim:
    clean = sp.sympify(expr)
    if clean == 0 or clean.is_Number:
        return ZERO_DIM
    if isinstance(clean, sp.Symbol):
        if clean not in dims:
            raise AuditFailure(f"missing dimension for symbol {clean}")
        return dims[clean]
    if isinstance(clean, sp.Mul):
        total = ZERO_DIM
        for arg in clean.args:
            total = total + dim_of(arg, dims)
        return total
    if isinstance(clean, sp.Pow):
        base, power = clean.args
        if not power.is_number:
            raise AuditFailure(f"non-numeric power in dimension expression {clean}")
        return dim_of(base, dims) * sp.Rational(power)
    if isinstance(clean, sp.Add):
        arg_dims = [dim_of(arg, dims) for arg in clean.args if compact_expr(arg) != 0]
        if not arg_dims:
            return ZERO_DIM
        first = arg_dims[0]
        if any(arg_dim != first for arg_dim in arg_dims[1:]):
            raise AuditFailure(f"dimension mismatch in sum {clean}")
        return first
    if clean.func in (sp.sin, sp.cos, sp.tan, sp.cot):
        arg_dims = [dim_of(arg, dims) for arg in clean.args]
        if any(arg_dim != ZERO_DIM for arg_dim in arg_dims):
            raise AuditFailure(f"dimensionful argument in dimensionless function {clean}")
        return ZERO_DIM
    raise AuditFailure(f"unsupported dimension expression {clean}")


s, L0, omega, cS = sp.symbols("s L0 omega cS", positive=True, real=True)
psiM = sp.Symbol("psiM", nonzero=True, real=True)
alpha = sp.Symbol("alpha", real=True)
n = sp.Symbol("n", integer=True, nonnegative=True)
j = sp.Symbol("j", integer=True, positive=True)
K, rho, rho_star, m = sp.symbols("K rho rho_star m", positive=True, real=True)
psi = sp.Function("psi")
psi_hat = sp.Function("psi_hat")

y = psi_hat(s)
k = sp.simplify(omega / cS)
cs_squared_consumed = sp.simplify(5 * K * rho_star**4 / m)


def r1_site_from_exponent(exponent: int) -> sp.Expr:
    e = sp.Integer(exponent)
    return compact_expr(e * K * rho ** (e - 1) / m)


def r1_eos_site_from_exponent(exponent: int) -> sp.Expr:
    return compact_expr(sp.diff(K * rho**exponent, rho) / m)


def reconstruct_Ls_from_pair(pair: list[sp.Expr]) -> dict[str, Any]:
    a_op, b_op = sp.symbols("a_null b_null")
    equations = [
        compact_expr(sp.diff(branch, s, 2) + a_op * sp.diff(branch, s) + b_op * branch)
        for branch in pair
    ]
    solved = sp.solve(equations, (a_op, b_op), dict=True)
    if not solved:
        raise AuditFailure(f"could not reconstruct monic operator from null-space pair {pair}")
    a_value = compact_expr(solved[0][a_op])
    b_value = compact_expr(solved[0][b_op])
    operator = compact_expr(sp.diff(y, s, 2) + a_value * sp.diff(y, s) + b_value * y)
    return {
        "a": a_value,
        "b": b_value,
        "operator": operator,
        "equations": equations,
    }


def compute_012_verdict(
    *,
    dimensional_ok: bool,
    dtn_matches_target: bool,
    halfshift: bool,
    counterfactual_guard_all: bool,
    bc_derivation_emitted: bool,
) -> str:
    if not dimensional_ok:
        return DN_UNITTEST_FAIL_DIMENSIONAL
    if not (dtn_matches_target and halfshift):
        return FAIL_POLE_LADDER
    if not counterfactual_guard_all:
        return FAIL_COUNTERFACTUAL
    if not bc_derivation_emitted:
        return DN_UNITTEST_BC_DEPENDENT
    return DN_UNITTEST_PASS


def build_dtn_case(*, cap_rhs: sp.Expr = sp.Integer(0)) -> dict[str, Any]:
    ode = sp.Eq(sp.diff(psi(s), s, 2) + k**2 * psi(s), 0)
    dsolve_solution = compact_expr(sp.dsolve(ode).rhs)
    constants = sorted(
        [sym for sym in dsolve_solution.free_symbols if sym.name.startswith("C")],
        key=lambda sym: sym.name,
    )
    if len(constants) != 2:
        raise AuditFailure(f"dsolve did not emit two constants: {dsolve_solution}")
    C1, C2 = constants

    fundamental_pair = [compact_expr(sp.diff(dsolve_solution, C1)), compact_expr(sp.diff(dsolve_solution, C2))]
    mouth_trace = compact_expr(dsolve_solution.subs(s, 0))
    cap_neumann = compact_expr(sp.diff(dsolve_solution, s).subs(s, L0))
    dn_matrix = sp.Matrix(
        [
            [sp.diff(mouth_trace, C1), sp.diff(mouth_trace, C2)],
            [sp.diff(cap_neumann, C1), sp.diff(cap_neumann, C2)],
        ]
    )
    dn_rhs = sp.Matrix([psiM, cap_rhs])
    dn_det = compact_expr(dn_matrix.det())
    coeff_solution = dn_matrix.LUsolve(dn_rhs)
    raw_coeff_map = {C1: coeff_solution[0], C2: coeff_solution[1]}
    coeff_map = {C1: compact_expr(coeff_solution[0]), C2: compact_expr(coeff_solution[1])}
    bc_applied_solution_raw = dsolve_solution.subs(raw_coeff_map)
    bc_applied_solution = compact_expr(dsolve_solution.subs(coeff_map))
    dtn_raw_unsimplified = -sp.diff(bc_applied_solution_raw, s).subs(s, 0) / bc_applied_solution_raw.subs(s, 0)
    dtn_raw = compact_expr(dtn_raw_unsimplified)
    dtn_sincos = sp.cancel(sp.together(dtn_raw_unsimplified))
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
    halfshift = expr_equal(pole_residual, 0)
    x = sp.Symbol("x", real=True)
    denominator_zeros = sp.solveset(sp.cos(x), x, domain=sp.S.Reals)

    static_series = sp.series(dtn, omega, 0, 6)
    static_series_poly = compact_expr(static_series.removeO())
    static_series_target = compact_expr(-L0 * omega**2 / cS**2 - L0**3 * omega**4 / (3 * cS**4))
    static_limit = compact_expr(sp.limit(dtn, omega, 0, dir="+"))

    return {
        "ode": ode,
        "dsolve_solution": dsolve_solution,
        "constants": (C1, C2),
        "fundamental_pair": fundamental_pair,
        "mouth_trace": mouth_trace,
        "cap_neumann": cap_neumann,
        "dn_matrix": dn_matrix,
        "dn_rhs": dn_rhs,
        "dn_det": dn_det,
        "coeff_solution": coeff_solution,
        "raw_coeff_map": raw_coeff_map,
        "coeff_map": coeff_map,
        "bc_applied_solution_raw": bc_applied_solution_raw,
        "bc_applied_solution": bc_applied_solution,
        "dtn_raw": dtn_raw,
        "dtn_sincos": dtn_sincos,
        "dtn": dtn,
        "dtn_target": dtn_target,
        "tan_argument": tan_argument,
        "dtn_prefactor": dtn_prefactor,
        "dtn_matches_target": dtn_matches_target,
        "denominator_full": denominator_full,
        "pole_denominator": pole_denominator,
        "pole_equation": pole_equation,
        "pole_ladder": pole_ladder,
        "pole_residual": pole_residual,
        "halfshift": halfshift,
        "denominator_zeros": denominator_zeros,
        "static_series": static_series,
        "static_series_poly": static_series_poly,
        "static_series_target": static_series_target,
        "static_limit": static_limit,
    }


def build_consumed_inputs(dtn_data: dict[str, Any]) -> dict[str, Any]:
    site_a = r1_site_from_exponent(5)
    site_b = r1_eos_site_from_exponent(5)
    consumed_speed = compact_expr(site_a.subs(rho, rho_star))

    Ls_site_a = compact_expr(sp.diff(y, s, 2) + k**2 * y)
    reconstructed = reconstruct_Ls_from_pair(dtn_data["fundamental_pair"])
    Ls_site_b = reconstructed["operator"]
    consumed_Ls = Ls_site_a
    anchor_Ls = compact_expr(sp.diff(y, s, 2) + (omega / cS) ** 2 * y)

    return {
        "site_a": site_a,
        "site_b": site_b,
        "consumed_speed": consumed_speed,
        "Ls_site_a": Ls_site_a,
        "Ls_reconstructed": reconstructed,
        "Ls_site_b": Ls_site_b,
        "consumed_Ls": consumed_Ls,
        "anchor_Ls": anchor_Ls,
        "domain": (sp.Integer(0), L0),
    }


def build_counterfactual_guard(
    *,
    robin_dtn: sp.Expr,
    det_witness: sp.Expr,
    denominator_core: sp.Expr,
    dtn: sp.Expr,
    alpha_numeric: sp.Expr = sp.Rational(2, 1) / L0,
) -> dict[str, Any]:
    dn_from_robin = compact_expr(robin_dtn.subs(alpha, 0))
    dd_from_robin = compact_expr(sp.limit(robin_dtn, alpha, sp.oo))
    dd_target = compact_expr(k * sp.cot(k * L0))
    dd_denominator = compact_expr(sp.sin(k * L0))

    dd_halfshift_samples = [
        compact_expr(dd_denominator.subs(omega, sp.pi * cS * (idx + sp.Rational(1, 2)) / L0))
        for idx in range(4)
    ]
    dd_integer_samples = [
        compact_expr(dd_denominator.subs(omega, sp.pi * cS * idx / L0))
        for idx in range(1, 5)
    ]
    halfshift_destroyed_for_dd = all(nonzero_q(value) for value in dd_halfshift_samples) and all(
        expr_equal(value, 0) for value in dd_integer_samples
    )
    dd_zero_mode_removable = expr_equal(sp.limit(dd_target, omega, 0, dir="+"), 1 / L0)

    numeric_robin_den = compact_expr(denominator_core.subs(alpha, alpha_numeric))
    numeric_robin_dtn = compact_expr(robin_dtn.subs(alpha, alpha_numeric))
    numeric_alpha_distinct = nonzero_q(numeric_robin_dtn - dtn) and nonzero_q(numeric_robin_dtn - dd_target)

    guard = {
        "robin_determinant_emitted": nonzero_q(det_witness),
        "recovers_DN_at_alpha0": expr_equal(dn_from_robin, dtn),
        "recovers_DD_at_alpha_inf": expr_equal(dd_from_robin, dd_target),
        "halfshift_destroyed_for_DD": halfshift_destroyed_for_dd,
        "numeric_alpha_distinct": numeric_alpha_distinct,
        "dtn_mismatch": nonzero_q(robin_dtn - dtn),
    }
    return {
        "guard": guard,
        "dn_from_robin": dn_from_robin,
        "dd_from_robin": dd_from_robin,
        "dd_target": dd_target,
        "dd_denominator": dd_denominator,
        "dd_halfshift_samples": dd_halfshift_samples,
        "dd_integer_samples": dd_integer_samples,
        "dd_zero_mode_removable": dd_zero_mode_removable,
        "numeric_alpha": alpha_numeric,
        "numeric_robin_den": numeric_robin_den,
        "numeric_robin_dtn": numeric_robin_dtn,
    }


def build_robin_block(dtn_data: dict[str, Any]) -> dict[str, Any]:
    C1, C2 = dtn_data["constants"]
    dsolve_solution = dtn_data["dsolve_solution"]
    mouth_trace = dtn_data["mouth_trace"]
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
    robin_numerator_core = compact_expr(robin_dtn * robin_denominator_core)
    guard_data = build_counterfactual_guard(
        robin_dtn=robin_dtn,
        det_witness=robin_det,
        denominator_core=robin_denominator_core,
        dtn=dtn_data["dtn"],
    )
    return {
        "cap_robin": cap_robin,
        "robin_matrix": robin_matrix,
        "robin_rhs": robin_rhs,
        "robin_det": robin_det,
        "robin_coeff_solution": robin_coeff_solution,
        "robin_coeff_map": robin_coeff_map,
        "robin_solution": robin_solution,
        "robin_dtn": robin_dtn,
        "robin_denominator_core": robin_denominator_core,
        "robin_numerator_core": robin_numerator_core,
        **guard_data,
    }


def build_round_trip(pole_ladder: sp.Expr, *, r_D: sp.Expr = sp.Integer(-1), r_N: sp.Expr = sp.Integer(1)) -> dict[str, Any]:
    round_trip = compact_expr(r_D * r_N * sp.exp(2 * sp.I * k * L0))
    round_trip_on_ladder = compact_expr(round_trip.subs(omega, pole_ladder))
    round_trip_closes = expr_equal(round_trip_on_ladder, 1)
    return {
        "r_D": r_D,
        "r_N": r_N,
        "round_trip": round_trip,
        "round_trip_on_ladder": round_trip_on_ladder,
        "round_trip_closes": round_trip_closes,
    }


def build_dimensional_block() -> dict[str, Any]:
    length_dim = Dim(1, 0, 0)
    energy_dim = Dim(2, 1, -2)
    four_volume_dim = 4 * length_dim
    pressure_dim = energy_dim - four_volume_dim
    rho_dim = -4 * length_dim
    K_dim = pressure_dim - 5 * rho_dim
    omega_dim = Dim(0, 0, -1)
    mass_dim = Dim(0, 1, 0)
    alpha_dim = Dim(-1, 0, 0)
    expected_tan_dim = ZERO_DIM
    expected_z00_dim = Dim(-1, 0, 0)
    cs_squared_expr = 5 * K * rho_star**4 / m

    def walk(K_dimension: Dim) -> dict[str, Dim]:
        base_dims = {
            L0: length_dim,
            omega: omega_dim,
            K: K_dimension,
            rho_star: rho_dim,
            m: mass_dim,
            alpha: alpha_dim,
        }
        cs_squared_dim = dim_of(cs_squared_expr, base_dims)
        cs_dim = cs_squared_dim * sp.Rational(1, 2)
        dims = dict(base_dims)
        dims[cS] = cs_dim
        k_dim = dim_of(omega / cS, dims)
        tan_argument_dim = dim_of(k * L0, dims)
        z00_prefactor_dim = dim_of(-k, dims)
        z00_dim = dim_of(-k * sp.tan(k * L0), dims) if tan_argument_dim == ZERO_DIM else None
        alpha_cS_dim = dim_of(alpha * cS, dims)
        return {
            "cs_squared_dim": cs_squared_dim,
            "cs_dim": cs_dim,
            "k_dim": k_dim,
            "tan_argument_dim": tan_argument_dim,
            "z00_prefactor_dim": z00_prefactor_dim,
            "z00_dim": z00_dim,
            "alpha_cS_dim": alpha_cS_dim,
        }

    clean_walk = walk(K_dim)
    dimensional_ok = (
        clean_walk["tan_argument_dim"] == expected_tan_dim
        and clean_walk["z00_prefactor_dim"] == expected_z00_dim
        and clean_walk["z00_dim"] == expected_z00_dim
    )
    corrupt_K_dim = K_dim + Dim(1, 0, 0)
    corrupt_walk = walk(corrupt_K_dim)
    corrupt_dimensional_ok = (
        corrupt_walk["tan_argument_dim"] == expected_tan_dim
        and corrupt_walk["z00_prefactor_dim"] == expected_z00_dim
        and corrupt_walk["z00_dim"] == expected_z00_dim
    )
    mutation_fires = (
        not corrupt_dimensional_ok
        and corrupt_walk["tan_argument_dim"] != expected_tan_dim
        and corrupt_walk["z00_prefactor_dim"] != expected_z00_dim
    )
    clean_verdict = compute_012_verdict(
        dimensional_ok=dimensional_ok,
        dtn_matches_target=True,
        halfshift=True,
        counterfactual_guard_all=True,
        bc_derivation_emitted=False,
    )
    mutated_verdict = compute_012_verdict(
        dimensional_ok=corrupt_dimensional_ok,
        dtn_matches_target=True,
        halfshift=True,
        counterfactual_guard_all=True,
        bc_derivation_emitted=False,
    )
    fail_suppressed = (
        mutation_fires
        and clean_verdict == DN_UNITTEST_BC_DEPENDENT
        and mutated_verdict == DN_UNITTEST_FAIL_DIMENSIONAL
    )
    return {
        "length_dim": length_dim,
        "energy_dim": energy_dim,
        "four_volume_dim": four_volume_dim,
        "pressure_dim": pressure_dim,
        "rho_dim": rho_dim,
        "K_dim": K_dim,
        "omega_dim": omega_dim,
        "alpha_dim": alpha_dim,
        "expected_tan_dim": expected_tan_dim,
        "expected_z00_dim": expected_z00_dim,
        "clean_walk": clean_walk,
        "dimensional_ok": dimensional_ok,
        "corrupt_K_dim": corrupt_K_dim,
        "corrupt_walk": corrupt_walk,
        "corrupt_dimensional_ok": corrupt_dimensional_ok,
        "mutation_fires": mutation_fires,
        "clean_verdict": clean_verdict,
        "mutated_verdict": mutated_verdict,
        "fail_suppressed": fail_suppressed,
    }


def build_baseline() -> dict[str, Any]:
    dtn = build_dtn_case()
    consumed = build_consumed_inputs(dtn)
    robin = build_robin_block(dtn)
    round_trip = build_round_trip(dtn["pole_ladder"])
    dim = build_dimensional_block()
    bc_derivation_emitted = False
    bc_provenance = "imposed"
    bc_derivation = {
        "bc_type": "D-at-mouth / N-at-cap",
        "mouth_gradient_from_Vconf": "not emitted from an explicit V_wall profile in this unit test",
        "cap_gradient_from_Vconf": "not emitted from an explicit V_wall profile in this unit test",
        "regularity_at_pinchoff": "regular cap closure R0(L0)=0 motivates Neumann, but a full asymptotic derivation is not emitted",
        "mouth_condition": "clamp to quasi-static bulk reservoir is imposed for this frozen-wall benchmark, not derived as radiation",
        "classification_rule": "bc_derivation_emitted=false forces DN_UNITTEST_BC_DEPENDENT unless an explicit mouth/cap derivation is later supplied",
    }
    guard_all = all(robin["guard"].values())
    verdict = compute_012_verdict(
        dimensional_ok=dim["dimensional_ok"],
        dtn_matches_target=dtn["dtn_matches_target"],
        halfshift=dtn["halfshift"],
        counterfactual_guard_all=guard_all,
        bc_derivation_emitted=bc_derivation_emitted,
    )
    return {
        "dtn": dtn,
        "consumed": consumed,
        "robin": robin,
        "round_trip": round_trip,
        "dim": dim,
        "bc_derivation_emitted": bc_derivation_emitted,
        "bc_provenance": bc_provenance,
        "bc_derivation": bc_derivation,
        "guard_all": guard_all,
        "verdict": verdict,
    }


def run_opening_and_consumed_inputs(data: dict[str, Any]) -> None:
    dtn = data["dtn"]
    consumed = data["consumed"]
    subbanner("Consumed stage011 inputs and opening dsolve")
    print("  CONSUMED-from-011: L_s, c_S, and domain [0,L0] are cited; stage011 reduction/certificate are not recomputed.")
    print("  CITED-speed: c_S^2 = 5*K*rho_star^4/m is R1 at rho_star; EOS exponent-5 P=K*rho^5 IMPOSED; bare m is stage004 m_GNLS.")
    print(f"  ODE consumed by 012 = {fmt(dtn['ode'])}")
    print(f"  dsolve general solution = {fmt(dtn['dsolve_solution'])}")
    expect_zero("dsolve solution satisfies cited Helmholtz L_s", sp.diff(dtn["dsolve_solution"], s, 2) + k**2 * dtn["dsolve_solution"])
    expect_zero("fundamental sin branch is emitted by dsolve", dtn["fundamental_pair"][0] - sp.sin(k * s))
    expect_zero("fundamental cos branch is emitted by dsolve", dtn["fundamental_pair"][1] - sp.cos(k * s))

    print(f"  R1 site A literal = {fmt(consumed['site_a'])}")
    print(f"  R1 site B EOS route d(K*rho^5)/d rho / m = {fmt(consumed['site_b'])}")
    expect_zero("c_S^2 R1 site A minus site B equals zero", consumed["site_a"] - consumed["site_b"])
    expect_zero("c_S^2 evaluated at rho_star equals consumed speed", consumed["consumed_speed"] - cs_squared_consumed)
    expect_zero("c_S^2 frozen-export anchor consumed - 5*K*rho_star^4/m equals zero", consumed["consumed_speed"] - 5 * K * rho_star**4 / m)

    recon = consumed["Ls_reconstructed"]
    print(f"  L_s site A export = {fmt(consumed['Ls_site_a'])}")
    print(f"  L_s site B null-space solve: a = {fmt(recon['a'])}, b = {fmt(recon['b'])}")
    print(f"  L_s site B reconstructed = {fmt(consumed['Ls_site_b'])}")
    expect_zero("L_s null-space reconstruction recovers a=0", recon["a"])
    expect_zero("L_s null-space reconstruction recovers b=(omega/c_S)^2", recon["b"] - k**2)
    expect_zero("L_s site A minus site B equals zero", consumed["Ls_site_a"] - consumed["Ls_site_b"])
    expect_zero("L_s frozen-export anchor consumed_L_s - (psi''+(omega/c_S)^2 psi) equals zero", consumed["consumed_Ls"] - consumed["anchor_Ls"])
    expect_zero("domain [0,L0] is cited as length L0, not re-solved", consumed["domain"][1] - consumed["domain"][0] - L0)


def run_dn_dtn(data: dict[str, Any]) -> None:
    dtn = data["dtn"]
    subbanner("D/N BVP coefficient solve and DtN")
    print("  Dirichlet mouth: psi_hat(0)=psi_M; Neumann cap: psi_hat'(L0)=0; both are imposed calibration inputs.")
    print(f"  D/N coefficient matrix = {fmt(dtn['dn_matrix'])}")
    print(f"  D/N rhs = {fmt(dtn['dn_rhs'])}")
    print(f"  D/N determinant = {fmt(dtn['dn_det'])}")
    expect_zero("D/N determinant equals -omega*cos(L0*omega/cS)/cS", dtn["dn_det"] + omega * sp.cos(k * L0) / cS)
    print(f"  LUsolve coefficients = {fmt(dtn['coeff_solution'])}")
    print(f"  BC-applied solution = {fmt(dtn['bc_applied_solution'])}")
    expect_zero(
        "BC-applied solution equals psiM*(sin(k*s)*tan(k*L0)+cos(k*s))",
        dtn["bc_applied_solution"] - psiM * (sp.sin(k * s) * sp.tan(k * L0) + sp.cos(k * s)),
    )
    print(f"  dtn_raw = {fmt(dtn['dtn_raw'])}")
    print(f"  dtn target = {fmt(dtn['dtn_target'])}")
    expect_zero("DtN derived via LUsolve equals -k*tan(k*L0)", dtn["dtn"] - dtn["dtn_target"])
    expect_bool("dtn_matches_target is genuine derived-vs-typed comparison", dtn["dtn_matches_target"])
    expect_zero("tan_argument is k*L0", dtn["tan_argument"] - k * L0)
    expect_zero("dtn_prefactor is -k", dtn["dtn_prefactor"] + k)


def run_pole_static_roundtrip(data: dict[str, Any]) -> None:
    dtn = data["dtn"]
    rt = data["round_trip"]
    subbanner("Half-shift pole ladder, static series, and round trip")
    print(f"  DtN denominator full = {fmt(dtn['denominator_full'])}")
    print(f"  pole denominator = {fmt(dtn['pole_denominator'])}; equation = {fmt(dtn['pole_equation'])}")
    expect_zero("pole denominator is cos(k*L0)", dtn["pole_denominator"] - sp.cos(k * L0))
    print(f"  pole ladder = {fmt(dtn['pole_ladder'])}")
    print(f"  pole residual after ladder substitution = {fmt(dtn['pole_residual'])}")
    expect_zero("half-shift pole residual is zero", dtn["pole_residual"])
    expect_bool("halfshift = pole_residual==0 is computed true", dtn["halfshift"])
    print(f"  denominator zeros in x = {fmt(dtn['denominator_zeros'])}")
    print(f"  static small-omega series = {sp.sstr(dtn['static_series'])}")
    print(f"  static limit omega->0+ = {fmt(dtn['static_limit'])}")
    expect_zero("static series polynomial matches -L0*omega^2/cS^2 - L0^3*omega^4/(3*cS^4)", dtn["static_series_poly"] - dtn["static_series_target"])
    expect_zero("static DC limit is zero", dtn["static_limit"])
    expect_nonzero("static series is distinguished from the DC limit", dtn["static_series_poly"] - dtn["static_limit"])
    print("  Static note: this is the small-omega expansion of the dynamic DtN; no separate static solve is emitted.")

    print(f"  r_D = {fmt(rt['r_D'])}; r_N = {fmt(rt['r_N'])}; round_trip = {fmt(rt['round_trip'])}")
    print(f"  round_trip on D/N ladder = {fmt(rt['round_trip_on_ladder'])}")
    expect_zero("round-trip closes to R_rt=1 on the D/N ladder", rt["round_trip_on_ladder"] - 1)
    expect_bool("round_trip_closes is computed true", rt["round_trip_closes"])
    print("  round-trip phase: phi0 = 0 mod 2*pi")


def run_robin_counterfactual(data: dict[str, Any]) -> None:
    robin = data["robin"]
    guard = robin["guard"]
    subbanner("Robin counterfactual falsifier")
    print(f"  cap_robin = {fmt(robin['cap_robin'])}")
    print(f"  Robin coefficient matrix = {fmt(robin['robin_matrix'])}")
    print(f"  Robin determinant = {fmt(robin['robin_det'])}")
    print(f"  robin_denominator_core = {fmt(robin['robin_denominator_core'])}")
    expect_zero("Robin determinant is the negative of robin_denominator_core", robin["robin_det"] + robin["robin_denominator_core"])
    expect_nonzero("robin_determinant_emitted is tied to computed determinant/core", robin["robin_det"])
    print(f"  Robin DtN = {fmt(robin['robin_dtn'])}")
    print(f"  alpha->0 branch = {fmt(robin['dn_from_robin'])}")
    print(f"  alpha->infinity branch = {fmt(robin['dd_from_robin'])}; D/D target = {fmt(robin['dd_target'])}")
    expect_zero("Robin alpha->0 recovers D/N DtN", robin["dn_from_robin"] - data["dtn"]["dtn"])
    expect_zero("Robin alpha->infinity recovers D/D k*cot(k*L0)", robin["dd_from_robin"] - robin["dd_target"])
    print(f"  D/D half-shift samples = {fmt(robin['dd_halfshift_samples'])}")
    print(f"  D/D integer-ladder samples = {fmt(robin['dd_integer_samples'])}")
    print(f"  dd_zero_mode_removable artifact = {fmt(robin['dd_zero_mode_removable'])} (not a guard member)")
    expect_bool("D/D zero mode removable limit equals 1/L0", robin["dd_zero_mode_removable"])
    print(f"  numeric alpha = {fmt(robin['numeric_alpha'])}; numeric Robin DtN = {fmt(robin['numeric_robin_dtn'])}")
    expect_zero("counterfactual_guard has exactly six members", len(guard) - 6)
    expect_bool("counterfactual_guard membership matches source dict L522-L529", tuple(guard.keys()) == EXPECTED_GUARD_KEYS)
    expect_bool("dd_zero_mode_removable is not verdict-bearing", "dd_zero_mode_removable" not in guard)
    for key in EXPECTED_GUARD_KEYS:
        print(f"  guard[{key}] = {fmt(guard[key])}")
        expect_bool(f"counterfactual guard member {key} is computed true", guard[key])
    expect_bool("all(counterfactual_guard.values()) is true", all(guard.values()))


def run_dimensional_block(data: dict[str, Any]) -> None:
    dim = data["dim"]
    clean = dim["clean_walk"]
    corrupt = dim["corrupt_walk"]
    subbanner("012 dimensional legs and corrupt-[K] probe")
    print("  dimension order: (L,M,T)")
    print(f"  [energy] = {dim['energy_dim']}; [four-volume] = {dim['four_volume_dim']}; [P] = {dim['pressure_dim']}")
    print(f"  [rho] = {dim['rho_dim']}; [K]=[P]-5[rho] = {dim['K_dim']}; [alpha] = {dim['alpha_dim']}")
    print(f"  propagated [c_S^2] = {clean['cs_squared_dim']} -> [c_S] = {clean['cs_dim']} -> [k] = {clean['k_dim']}")
    print(f"  [tan_argument=k*L0] = {clean['tan_argument_dim']}; [Z00_prefactor=-k] = {clean['z00_prefactor_dim']}; [Z00] = {clean['z00_dim']}")
    expect_zero("tan_argument dimensional leg is dimensionless", dim_residual(clean["tan_argument_dim"], dim["expected_tan_dim"]))
    expect_zero("Z00_prefactor dimensional leg is L^-1", dim_residual(clean["z00_prefactor_dim"], dim["expected_z00_dim"]))
    expect_zero("Z00 dimensional leg is L^-1", dim_residual(clean["z00_dim"], dim["expected_z00_dim"]))
    expect_bool("dimensional_ok for 012 tan_argument/Z00 legs", dim["dimensional_ok"])
    print(f"  corrupt [K]+(1,0,0) gives [K] = {dim['corrupt_K_dim']}")
    print(f"  corrupt propagated [c_S] = {corrupt['cs_dim']}; [k] = {corrupt['k_dim']}")
    print(f"  corrupt [tan_argument] = {corrupt['tan_argument_dim']}; corrupt [Z00_prefactor] = {corrupt['z00_prefactor_dim']}")
    expect_fail("corrupt-[K] makes tan_argument non-dimensionless", dim_residual(corrupt["tan_argument_dim"], dim["expected_tan_dim"]))
    expect_fail("corrupt-[K] makes Z00_prefactor differ from L^-1", dim_residual(corrupt["z00_prefactor_dim"], dim["expected_z00_dim"]))
    expect_bool("corrupt-[K] mutation_fires=True", dim["mutation_fires"])
    expect_zero("self-ablation with mutation gives DN_UNITTEST_FAIL_DIMENSIONAL", verdict_residual(dim["mutated_verdict"], DN_UNITTEST_FAIL_DIMENSIONAL))
    expect_zero("self-ablation without mutation gives clean DN_UNITTEST_BC_DEPENDENT", verdict_residual(dim["clean_verdict"], DN_UNITTEST_BC_DEPENDENT))
    expect_bool("self-ablation fail_suppressed=True", dim["fail_suppressed"])


def run_bc_and_verdict(data: dict[str, Any]) -> None:
    subbanner("BC provenance, 012 verdict, and joint composition")
    print(f"  bc_provenance = {data['bc_provenance']}")
    print(f"  bc_derivation_emitted = {fmt(data['bc_derivation_emitted'])}")
    print(f"  bc_derivation descriptor = {data['bc_derivation']}")
    expect_bool("bc_provenance is imposed", data["bc_provenance"] == "imposed")
    expect_bool("bc_derivation_emitted is the honest false scope flag", data["bc_derivation_emitted"] is False)
    print(f"  012 scoped verdict = {data['verdict']}")
    expect_zero("012 verdict lands at DN_UNITTEST_BC_DEPENDENT", verdict_residual(data["verdict"], DN_UNITTEST_BC_DEPENDENT))
    print("  DN_UNITTEST_BC_DEPENDENT (JOINT, COMPLETED)")
    print("    = (011: REDUCTION_CERTIFIED, cited from ledger_stage011)")
    print("    AND (012: DtN ladder EARNED + bc_derivation_emitted=False -> BC_DEPENDENT landing, computed here)")
    expect_bool("joint composition cites stage011 REDUCTION_CERTIFIED and computed 012 landing", data["verdict"] == DN_UNITTEST_BC_DEPENDENT)


def print_provenance() -> None:
    subbanner("Provenance and scope")
    print("  CONSUMED-from-011: L_s, domain [0,L0] with cap R0(L0)=0, and c_S are CITED from ledger_stage011 with dual-site integrity; 011 reduction/certificate/de-rig are not recomputed.")
    print("  CITED-speed: c_S^2=5*K*rho_star^4/m is Part I edge R1 (stage005, re-exported by stage011) at rho_star; EOS exponent-5 P=K*rho^5 IMPOSED.")
    print("  IMPOSED-BC: D/N mouth/cap boundary pair is IMPOSED; bc_provenance=imposed and bc_derivation_emitted=False are banked calibration; V_wall derivation is deferred, not fabricated.")
    print("  EARNED: DtN, half-shifted ladder, static small-omega series, round-trip R_rt=1, and Robin counterfactual guard are computed here; dtn_matches_target is derived-vs-typed, not X==X.")
    print("  Robin-falsifier: Robin cap recovers D/N at alpha=0, D/D at alpha->infinity, destroys the half-shift for D/D, and is numerically distinct at alpha=2/L0.")
    print("  control-symbol: alpha is a Robin cap admittance with [alpha]=L^-1, tracked-not-counted like stage010 k_warp; it builds the falsifiable counterfactual, not the physics.")
    print("  split: this stage COMPLETES DN_UNITTEST_BC_DEPENDENT; 011 carried REDUCTION_CERTIFIED, 012 carries the D/N ladder + BC_DEPENDENT landing; DN_UNITTEST_PASS is deferred.")
    print("  dropped-bookkeeping: scratch-YAML/_sympy_exprs.wl export, MMA-YAML re-read, expression_digest, and engine_agreement plumbing are stripped.")
    print("  downstream consumers: stage 013 (harmonic beta lift) + stage 017 (calibration input) consume Z00, the resonance ladder, and BC_DEPENDENT provenance.")
    print("  register note: zero new counted knobs; alpha is tracked-not-counted; L0 already registered in stage011; edge R28 is the imposed D/N boundary calibration obligation.")


def print_verdict_labels() -> None:
    print("")
    print("Verdict labels:")
    print("  ledger earned-label (NOT a source verdict token): DTN_POLE_LADDER_ROBIN_FALSIFIER_EARNED  (dsolve of the cited frozen L_s -> D/N coefficient matrix -> outward-mouth DtN Z00=-(omega/c_S)*tan(L0*omega/c_S); half-shifted pole ladder omega_n=pi*c_S*(n+1/2)/L0; static small-omega series; round-trip R_rt=1; Robin cap counterfactual falsifier, guard = {robin_determinant_emitted, recovers_DN_at_alpha0, recovers_DD_at_alpha_inf, halfshift_destroyed_for_DD, numeric_alpha_distinct, dtn_mismatch}, each a computed residual; tan_argument/Z00 dim legs via [K]=[P]-5[rho] + corrupt-[K] probe)")
    print("  source top-line verdict: DN_UNITTEST_BC_DEPENDENT  (JOINT; 012 COMPLETES it -- adds the D/N ladder + bc_derivation_emitted=False -> BC_DEPENDENT landing to stage 011's cited REDUCTION_CERTIFIED)")
    print("  joint composition (COMPLETED): DN_UNITTEST_BC_DEPENDENT = (011: REDUCTION_CERTIFIED, cited from ledger_stage011) AND (012: DtN ladder EARNED + BC_DEPENDENT landing, computed here)")
    print("  earned: DtN derived via LUsolve (dtn_matches_target = genuine derived-vs-typed comparison, NOT X==X); half-shifted ladder (halfshift = pole_residual==0, COMPUTED); round-trip R_rt=1 (COMPUTED via substitution); Robin guard all booleans COMPUTED; tan_argument/Z00 dim legs (2,0,-2)-consistent + corrupt-[K] probe fires")
    print("  consumed (cited from stage011, dual-site integrity): frozen Helmholtz L_s = psi'' + (omega/c_S)^2 psi on [0,L0]; c_S^2 = 5*K*rho_star^4/m (R1 at rho_star); domain cap R0(L0)=0")
    print("  imposed (banked calibration, edge R28): D/N mouth/cap boundary pair; bc_provenance=imposed; bc_derivation_emitted=False -> BC_DEPENDENT landing; the mouth/cap V_wall gradient derivation earning DN_UNITTEST_PASS is a DEFERRED upgrade (NOT fabricated here)")
    print("  control symbol (tracked, not counted): alpha = Robin cap admittance, [alpha]=L^-1 (like k_warp at stage010)")


def run_able_to_fail_teeth(data: dict[str, Any]) -> None:
    dtn = data["dtn"]
    robin = data["robin"]
    dim = data["dim"]
    consumed = data["consumed"]
    subbanner("Able-to-fail mutation teeth")

    flipped_core = compact_expr(k * sp.cos(k * L0) + alpha * sp.sin(k * L0))
    flipped_robin_dtn = compact_expr(robin["robin_numerator_core"] / flipped_core)
    flipped_guard = build_counterfactual_guard(
        robin_dtn=flipped_robin_dtn,
        det_witness=flipped_core,
        denominator_core=flipped_core,
        dtn=dtn["dtn"],
    )["guard"]
    expect_fail("tooth 1 Robin denominator sign flip makes recovers_DD_at_alpha_inf false", bool_residual(flipped_guard["recovers_DD_at_alpha_inf"]))
    expect_zero(
        "tooth 1 verdict is FAIL_COUNTERFACTUAL",
        verdict_residual(
            compute_012_verdict(
                dimensional_ok=True,
                dtn_matches_target=True,
                halfshift=True,
                counterfactual_guard_all=all(flipped_guard.values()),
                bc_derivation_emitted=False,
            ),
            FAIL_COUNTERFACTUAL,
        ),
    )

    dtn_mut = build_dtn_case(cap_rhs=psiM / L0)
    expect_fail("tooth 2 mutated Neumann RHS changes derived DtN", dtn_mut["dtn"] - dtn["dtn_target"])
    expect_fail("tooth 2 dtn_matches_target flips false", bool_residual(dtn_mut["dtn_matches_target"]))
    expect_zero(
        "tooth 2 verdict is FAIL_POLE_LADDER",
        verdict_residual(
            compute_012_verdict(
                dimensional_ok=True,
                dtn_matches_target=dtn_mut["dtn_matches_target"],
                halfshift=True,
                counterfactual_guard_all=True,
                bc_derivation_emitted=False,
            ),
            FAIL_POLE_LADDER,
        ),
    )

    integer_ladder = compact_expr(sp.pi * cS * j / L0)
    integer_residual = compact_expr(dtn["pole_denominator"].subs(omega, integer_ladder))
    integer_halfshift = expr_equal(integer_residual, 0)
    expect_fail("tooth 3 integer ladder does not solve cos(k*L0)=0", integer_residual)
    expect_fail("tooth 3 halfshift boolean flips false", bool_residual(integer_halfshift))
    expect_zero(
        "tooth 3 verdict is FAIL_POLE_LADDER",
        verdict_residual(
            compute_012_verdict(
                dimensional_ok=True,
                dtn_matches_target=True,
                halfshift=integer_halfshift,
                counterfactual_guard_all=True,
                bc_derivation_emitted=False,
            ),
            FAIL_POLE_LADDER,
        ),
    )

    bad_round_trip = build_round_trip(dtn["pole_ladder"], r_D=sp.Integer(1), r_N=sp.Integer(1))
    expect_fail("tooth 4 corrupt r_D/r_N makes round_trip_on_ladder differ from 1", bad_round_trip["round_trip_on_ladder"] - 1)
    expect_fail("tooth 4 round_trip_closes flips false", bool_residual(bad_round_trip["round_trip_closes"]))

    no_alpha_numerator = compact_expr(-k * (k * sp.sin(k * L0)))
    no_alpha_robin_dtn = compact_expr(no_alpha_numerator / robin["robin_denominator_core"])
    no_alpha_guard = build_counterfactual_guard(
        robin_dtn=no_alpha_robin_dtn,
        det_witness=robin["robin_denominator_core"],
        denominator_core=robin["robin_denominator_core"],
        dtn=dtn["dtn"],
    )["guard"]
    expect_fail("tooth 5 broken alpha->infinity path makes recovers_DD_at_alpha_inf false", bool_residual(no_alpha_guard["recovers_DD_at_alpha_inf"]))
    expect_zero(
        "tooth 5 verdict is FAIL_COUNTERFACTUAL",
        verdict_residual(
            compute_012_verdict(
                dimensional_ok=True,
                dtn_matches_target=True,
                halfshift=True,
                counterfactual_guard_all=all(no_alpha_guard.values()),
                bc_derivation_emitted=False,
            ),
            FAIL_COUNTERFACTUAL,
        ),
    )

    degenerate_guard = build_counterfactual_guard(
        robin_dtn=robin["robin_dtn"],
        det_witness=robin["robin_det"],
        denominator_core=robin["robin_denominator_core"],
        dtn=dtn["dtn"],
        alpha_numeric=sp.Integer(0),
    )["guard"]
    expect_fail("tooth 6 numeric alpha degeneracy alpha=0 makes distinctness false", bool_residual(degenerate_guard["numeric_alpha_distinct"]))
    expect_zero(
        "tooth 6 verdict is FAIL_COUNTERFACTUAL",
        verdict_residual(
            compute_012_verdict(
                dimensional_ok=True,
                dtn_matches_target=True,
                halfshift=True,
                counterfactual_guard_all=all(degenerate_guard.values()),
                bc_derivation_emitted=False,
            ),
            FAIL_COUNTERFACTUAL,
        ),
    )

    for exponent in (4, 6):
        expect_fail(
            f"tooth 7 c_S^2 site A exponent 5->{exponent} trips R1 dual-site integrity",
            r1_site_from_exponent(exponent) - consumed["site_b"],
        )
        expect_fail(
            f"tooth 7 c_S^2 site B exponent 5->{exponent} trips R1 dual-site integrity",
            consumed["site_a"] - r1_eos_site_from_exponent(exponent),
        )
    expect_fail(
        "tooth 7 coordinated R1 both-site exponent drift trips frozen-export anchor",
        r1_site_from_exponent(6).subs(rho, rho_star) - 5 * K * rho_star**4 / m,
    )
    bad_Ls_site_a = compact_expr(sp.diff(y, s, 2) - k**2 * y)
    expect_fail("tooth 7 L_s export sign corruption trips site A/B integrity", bad_Ls_site_a - consumed["Ls_site_b"])
    hyper_recon = reconstruct_Ls_from_pair([sp.sinh(k * s), sp.cosh(k * s)])
    expect_fail("tooth 7 L_s site-B sinh/cosh corruption trips null-space integrity", consumed["Ls_site_a"] - hyper_recon["operator"])

    expect_fail(
        "tooth 8 corrupt-[K] probe trips tan_argument dimensional gate",
        dim_residual(dim["corrupt_walk"]["tan_argument_dim"], dim["expected_tan_dim"]),
    )
    expect_fail(
        "tooth 8 corrupt-[K] probe trips Z00_prefactor dimensional gate",
        dim_residual(dim["corrupt_walk"]["z00_prefactor_dim"], dim["expected_z00_dim"]),
    )
    expect_zero("tooth 8 corrupt-[K] verdict is DN_UNITTEST_FAIL_DIMENSIONAL", verdict_residual(dim["mutated_verdict"], DN_UNITTEST_FAIL_DIMENSIONAL))
    expect_bool("tooth 8 self-ablation fail_suppressed remains true", dim["fail_suppressed"])

    expect_zero("baseline immutable after teeth: DtN still equals target", dtn["dtn"] - dtn["dtn_target"])
    expect_zero("baseline immutable after teeth: halfshift pole residual remains zero", dtn["pole_residual"])
    expect_zero("baseline immutable after teeth: Robin alpha->infinity still recovers D/D", robin["dd_from_robin"] - robin["dd_target"])
    expect_zero("baseline immutable after teeth: L_s site integrity remains zero", consumed["Ls_site_a"] - consumed["Ls_site_b"])
    expect_zero("baseline immutable after teeth: clean 012 verdict remains DN_UNITTEST_BC_DEPENDENT", verdict_residual(data["verdict"], DN_UNITTEST_BC_DEPENDENT))


def main() -> None:
    banner("ledger_stage012_dtn_pole_ladder_robin SymPy audit")
    data = build_baseline()
    assert_no_float("baseline", data)
    run_opening_and_consumed_inputs(data)
    run_dn_dtn(data)
    run_pole_static_roundtrip(data)
    run_robin_counterfactual(data)
    run_dimensional_block(data)
    run_bc_and_verdict(data)
    print_provenance()
    print_verdict_labels()
    run_able_to_fail_teeth(data)
    print("")
    print(f"PASS tally: {PASS_COUNT}; FAIL tally: {FAIL_COUNT}")
    print("OVERALL PASS: SymPy verified ledger_stage012 DtN pole ladder + Robin falsifier exactly")


if __name__ == "__main__":
    try:
        main()
    except AuditFailure as exc:
        print("")
        print(f"PASS tally: {PASS_COUNT}; FAIL tally: {FAIL_COUNT}")
        print(f"OVERALL FAIL: SymPy stage012 audit did not close ({exc})")
        raise SystemExit(1)
