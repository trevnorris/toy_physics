#!/usr/bin/env python3
"""SymPy audit engine for S11 stray-longitudinal measurements.

The script streams one computed object per tag and writes the accumulated
Python ledger for the next step.  Its physics inputs are the equations in
``directives/S11_SHARED_PHYSICS.md``.
"""

from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations, combinations_with_replacement
from pathlib import Path
import hashlib
import math
import sys
import time

import sympy as sp
from sympy.core.symbol import Str
from sympy.polys.matrices import DomainMatrix

SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parents[2]
OUT_PATH = SCRIPT_DIR / "out" / "S11_stray_longitudinal_sympy_audit.out"
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from S10_exports import LEDGER as INCOMING_LEDGER  # noqa: E402


CLASS_TAGS = {"KNOB", "STRUCTURAL", "COORDINATE", "CONTROL", "PREMISE", "DERIVED"}
NOT_APPLICABLE = Str("NOT_APPLICABLE")
NOT_DEFINED_ON_COMPONENT = Str("NOT_DEFINED_ON_COMPONENT")
NOT_A_PURE_POWER = Str("NOT_A_PURE_POWER")
UNDEFINED_RATIO = Str("UNDEFINED_RATIO")
CODING_CONSISTENCY_ONLY = Str("CODING_CONSISTENCY_ONLY")
CROSS_STEP_CONSISTENCY_ONLY = Str("CROSS_STEP_CONSISTENCY_ONLY")
M_B_TOKEN = Str("M_B")


DECLARED_SYMBOLS: dict[sp.Symbol, dict[str, str]] = {}
ISSUES: list[str] = []
CURRENT_CELL: tuple[str, int] | None = None
CURRENT_ISSUE_CONTEXT: str | None = None
CUSTOM_QUANTITY_NAMES: set[str] = set()


def record_issue(message: str, *, package: str | None = None, n: int | None = None,
                 root: int | None = None, stratum: int | None = None,
                 context: str | None = None) -> None:
    if (package is None or n is None) and CURRENT_CELL is not None:
        package, n = CURRENT_CELL
    if package is not None and n is not None:
        prefix = f"{package} D{n}"
        if root is not None:
            prefix = f"{prefix} root {root}"
        if stratum is not None:
            prefix = f"{prefix} stratum {stratum}"
    else:
        prefix = context or CURRENT_ISSUE_CONTEXT or "GLOBAL"
    ISSUES.append(f"{prefix}: {message}")


def register_symbol(symbol: sp.Symbol, class_tag: str, description: str) -> sp.Symbol:
    if class_tag not in CLASS_TAGS:
        raise ValueError(f"unknown class tag {class_tag!r}")
    DECLARED_SYMBOLS[symbol] = {"class": class_tag, "description": description}
    return symbol


def declared_symbol(name: str, class_tag: str, description: str, **assumptions) -> sp.Symbol:
    return register_symbol(sp.Symbol(name, **assumptions), class_tag, description)


D = register_symbol(INCOMING_LEDGER["D"]["value"], "STRUCTURAL", "brane spatial dimension")
rho_br = register_symbol(INCOMING_LEDGER["rho_br"]["value"], "KNOB", "brane inertia density")
mu_R = register_symbol(INCOMING_LEDGER["mu_R"]["value"], "KNOB", "brane rotational stiffness")
B_comp = declared_symbol("B_comp", "KNOB", "brane compression stiffness", positive=True)
mu_br = declared_symbol("mu_br", "KNOB", "traceless stiffness", positive=True)
beta = declared_symbol("beta", "KNOB", "extra invariant coefficient", real=True)
s = declared_symbol("s", "CONTROL", "compression coefficient scale", positive=True)
s_rho = declared_symbol("s_rho", "CONTROL", "anisotropic inertia scale", positive=True)
c_s0 = declared_symbol("c_s0", "KNOB", "bulk scalar sound speed", positive=True)

omegaSquared = register_symbol(
    INCOMING_LEDGER.get("omegaSquared", {}).get("value", sp.Symbol("omegaSquared", real=True)),
    "COORDINATE",
    "squared frequency spectral coordinate",
)
lambdaScale = register_symbol(
    INCOMING_LEDGER.get("lambdaScale", {}).get("value", sp.Symbol("lambdaScale", positive=True)),
    "CONTROL",
    "wavevector scaling control",
)
kwSquared = declared_symbol("kwSquared", "COORDINATE", "normal wavevector square", real=True)
bulk_amplitude = declared_symbol("A", "COORDINATE", "bulk scalar amplitude", real=True)
phase = declared_symbol("phase", "COORDINATE", "plane-wave phase", real=True)
t = declared_symbol("t", "COORDINATE", "time coordinate", real=True)

X_ALL = tuple(declared_symbol(f"x{i}", "COORDINATE", "brane coordinate", real=True) for i in range(1, 6))
K_ALL = tuple(declared_symbol(f"k{i}", "COORDINATE", "wavevector component", real=True) for i in range(1, 6))
A_ALL = tuple(declared_symbol(f"a{i}", "COORDINATE", "brane displacement amplitude", real=True) for i in range(1, 6))
QG_ALL = tuple(declared_symbol(f"g_{i}", "COORDINATE", "row-major gradient coordinate", real=True) for i in range(1, 26))
G7_ALL = tuple(
    declared_symbol(f"g_{i}_{j}", "COORDINATE", "independent gradient coordinate", real=True)
    for i in range(1, 6)
    for j in range(1, 6)
)

L_DIM = sp.ImmutableMatrix([1, 0, 0])
T_DIM = sp.ImmutableMatrix([0, 1, 0])
M_DIM = sp.ImmutableMatrix([0, 0, 1])
ZERO_DIM = sp.ImmutableMatrix([0, 0, 0])
FIELD_DIM = L_DIM
ENERGY_DIM = sp.ImmutableMatrix([2, -2, 1])
ENERGY_DENSITY_DIM = sp.ImmutableMatrix([2 - D, -2, 1])
SPATIAL_DERIVATIVE_DIM = sp.ImmutableMatrix([-1, 0, 0])
TIME_DERIVATIVE_DIM = sp.ImmutableMatrix([0, -1, 0])
OMEGA_SQUARED_DIM = sp.ImmutableMatrix([0, -2, 0])
WAVEVECTOR_DIM = SPATIAL_DERIVATIVE_DIM


COEFFICIENT_NAME_MAP = {
    "rho_br": "rho_br",
    "mu_R": "mu_R",
    "B_comp": "B_comp",
    "mu_br": "mu_br",
    "beta": "beta",
    "s": "s",
    "s_rho": "s_rho",
    "c_s0": "c_s0",
}


PACKAGE_DIMS = {
    "MAIN": (2, 3, 4, 5),
    "XFORM_CURLONLY": (2, 3, 4, 5),
    "XFORM_EXTRA": (2, 3, 4, 5),
    "XFORM_DIVONLY": (3, 4),
    "XFORM_TRACELESS": (3, 4),
    "XCOEF_BSCALE": (3,),
    "XCOEF_BSIGN": (3,),
    "XKIN_ANISO": (2, 3, 4),
}
PACKAGE_ORDER = tuple(PACKAGE_DIMS)
PRIMARY_PACKAGE = "MAIN"


@dataclass(frozen=True)
class Term:
    factor: sp.Expr
    density_name: str
    density: sp.Expr

    @property
    def expression(self) -> sp.Expr:
        return sp.expand(self.factor * self.density)


@dataclass
class PackageBuild:
    kinetic_terms: tuple[Term, ...]
    stiffness_terms: tuple[Term, ...]
    coefficient_ordering: tuple[sp.Symbol, ...]
    lagrangian: sp.Expr


@dataclass
class RouteSelection:
    token: Str
    matrix: sp.ImmutableMatrix


@dataclass
class RootData:
    value: sp.Expr
    multiplicity: int


@dataclass
class StratumCandidate:
    source: str
    locus: dict[str, object]
    root: int | None


@dataclass
class CandidateRow:
    base_key: str
    value: object
    class_tag: str
    source_tag: str
    dimension_key: str | None = None
    description: str | None = None


def file_digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def casify(value: object) -> object:
    if isinstance(value, sp.MatrixBase):
        return sp.ImmutableMatrix(value)
    if isinstance(value, dict):
        return sp.Tuple(
            *[
                sp.Tuple(casify(key), casify(val))
                for key, val in value.items()
            ]
        )
    if isinstance(value, (list, tuple)):
        return sp.Tuple(*(casify(item) for item in value))
    if isinstance(value, str):
        return Str(value)
    if value is True:
        return sp.true
    if value is False:
        return sp.false
    if value is None:
        return Str("None")
    if isinstance(value, int):
        return sp.Integer(value)
    return value


def render(value: object) -> str:
    return sp.srepr(casify(value))


def display(value: object) -> str:
    return str(casify(value))


def tag_name(package: str | None, n: int | None, quantity: str, root: int | None = None,
             stratum: int | None = None, local: bool = False) -> str:
    prefix = "PY_S11_LOCAL" if local else "PY_S11"
    parts = [prefix]
    if package is not None:
        parts.append(package)
    if n is not None:
        parts.append(f"D{n}")
    if stratum is not None:
        parts.append(f"STRATUM{stratum}")
    if root is not None:
        parts.append(f"ROOT{root}")
    parts.append(quantity)
    return "_".join(parts)


def base_key(quantity: str, n: int | None = None, root: int | None = None,
             stratum: int | None = None) -> str:
    parts: list[str] = []
    if stratum is not None:
        parts.append(f"stratum{stratum}")
    if root is not None:
        parts.append(f"root{root}")
    parts.append(quantity.lower())
    key = "_".join(parts)
    if n is not None:
        key = f"{key}_d{n}"
    return key


class Emitter:
    def __init__(self) -> None:
        self.count = 0
        self.local_tags: list[str] = []
        self.export_candidates: list[CandidateRow] = []
        self.emitted: dict[str, object] = {}

    def emit(self, tag: str, payload: object, *, export_key: str | None = None,
             class_tag: str = "DERIVED", dimension_key: str | None = None,
             description: str | None = None) -> object:
        payload = casify(payload)
        print(f"{tag}: {render(payload)}", flush=True)
        self.count += 1
        self.emitted[tag] = payload
        if "_LOCAL_" in tag:
            self.local_tags.append(tag)
        if export_key is not None:
            self.export_candidates.append(
                CandidateRow(export_key, payload, class_tag, tag, dimension_key, description)
            )
        return payload


EMITTER = Emitter()


def emit_quantity(package: str, n: int, quantity: str, payload: object, *,
                  root: int | None = None, stratum: int | None = None,
                  local: bool = False, export: bool = True,
                  class_tag: str = "DERIVED", dimension_key: str | None = None) -> object:
    tag = tag_name(package, n, quantity, root=root, stratum=stratum, local=local)
    export_key = None
    if export and not local and package == PRIMARY_PACKAGE:
        export_key = base_key(quantity, n=n, root=root, stratum=stratum)
    return EMITTER.emit(tag, payload, export_key=export_key, class_tag=class_tag, dimension_key=dimension_key)


def emit_custom_name(name: str) -> None:
    CUSTOM_QUANTITY_NAMES.add(name)


def matrix_from_rows(rows: list[list[object]], cols: int) -> sp.ImmutableMatrix:
    if not rows:
        return sp.ImmutableMatrix.zeros(0, cols)
    return sp.ImmutableMatrix(rows)


def symbols_for_dimension(n: int) -> tuple[tuple[sp.Symbol, ...], tuple[sp.Symbol, ...], tuple[sp.Symbol, ...]]:
    return X_ALL[:n], K_ALL[:n], A_ALL[:n]


def derivative_placeholders(n: int) -> tuple[sp.ImmutableMatrix, tuple[sp.Symbol, ...]]:
    g = sp.ImmutableMatrix(
        n, n, lambda i, j: sp.Symbol(f"G_{i + 1}_{j + 1}", real=True)
    )
    v = tuple(sp.Symbol(f"V_{j + 1}", real=True) for j in range(n))
    return g, v


def u_functions(n: int) -> tuple[tuple[sp.Symbol, ...], tuple[sp.Expr, ...]]:
    xs = X_ALL[:n]
    args = (*xs, t)
    funcs = tuple(sp.Function(f"u{j + 1}")(*args) for j in range(n))
    return xs, funcs


def stiffness_densities(g: sp.ImmutableMatrix) -> dict[str, sp.Expr]:
    n = g.rows
    trace = sum(g[i, i] for i in range(n))
    sym = (g + g.T) / 2
    stl = sym - sp.eye(n) * trace / n
    return {
        "S_curl": sp.Rational(1, 2) * sum((g[i, j] - g[j, i]) ** 2 for i in range(n) for j in range(n)),
        "S_div": trace ** 2,
        "S_symtl": sum(stl[i, j] ** 2 for i in range(n) for j in range(n)),
    }


def coordinate_substitution(n: int, g: sp.ImmutableMatrix, v: tuple[sp.Symbol, ...]) -> dict[sp.Symbol, sp.Expr]:
    xs, funcs = u_functions(n)
    sub: dict[sp.Symbol, sp.Expr] = {}
    for i in range(n):
        for j in range(n):
            sub[g[i, j]] = sp.Derivative(funcs[j], xs[i], evaluate=False)
    for j in range(n):
        sub[v[j]] = sp.Derivative(funcs[j], t, evaluate=False)
    return sub


def to_coordinate(expr: sp.Expr, n: int, g: sp.ImmutableMatrix, v: tuple[sp.Symbol, ...]) -> sp.Expr:
    return sp.expand(expr.xreplace(coordinate_substitution(n, g, v)))


def euler_lagrange_from_placeholders(expr: sp.Expr, n: int, g: sp.ImmutableMatrix,
                                     v: tuple[sp.Symbol, ...]) -> tuple[sp.Expr, ...]:
    xs, funcs = u_functions(n)
    coord_sub = coordinate_substitution(n, g, v)
    system = []
    for j in range(n):
        time_part = sp.diff(expr, v[j])
        time_coord = sum(
            sp.diff(time_part, v[m]) * sp.Derivative(funcs[m], (t, 2), evaluate=False)
            for m in range(n)
        )
        space_coord = 0
        for i in range(n):
            grad_part = sp.diff(expr, g[i, j])
            for ell in range(n):
                for m in range(n):
                    coeff = sp.diff(grad_part, g[ell, m])
                    if coeff != 0:
                        space_coord += coeff * sp.Derivative(funcs[m], xs[i], xs[ell], evaluate=False)
        system.append(sp.expand((time_coord + space_coord).xreplace(coord_sub)))
    return tuple(system)


def route_a_matrix(expr: sp.Expr, n: int, g: sp.ImmutableMatrix, v: tuple[sp.Symbol, ...],
                   k: tuple[sp.Symbol, ...], a: tuple[sp.Symbol, ...]) -> tuple[sp.ImmutableMatrix, sp.Expr]:
    rows = []
    for j in range(n):
        time_part = sp.diff(expr, v[j])
        amp_eq = sum(sp.diff(time_part, v[m]) * (-omegaSquared * a[m]) for m in range(n))
        for i in range(n):
            grad_part = sp.diff(expr, g[i, j])
            for ell in range(n):
                for m in range(n):
                    coeff = sp.diff(grad_part, g[ell, m])
                    if coeff != 0:
                        amp_eq += coeff * (-k[i] * k[ell] * a[m])
        rows.append([sp.simplify(sp.diff(sp.expand(amp_eq), amp)) for amp in a])
    return sp.ImmutableMatrix(rows), sp.cos(phase)


def period_average(expr: sp.Expr) -> sp.Expr:
    sin_phase = sp.Symbol("sin_phase", real=True)
    tmp_omega = sp.Symbol("omega_linear", real=True)
    averaged = sp.expand(expr).xreplace({sin_phase ** 2: sp.Rational(1, 2)})
    averaged = averaged.xreplace({tmp_omega ** 2: omegaSquared})
    averaged = averaged.subs(sin_phase ** 2, sp.Rational(1, 2)).subs(tmp_omega ** 2, omegaSquared)
    return sp.simplify(averaged)


def route_b_matrix(expr: sp.Expr, n: int, g: sp.ImmutableMatrix, v: tuple[sp.Symbol, ...],
                   k: tuple[sp.Symbol, ...], a: tuple[sp.Symbol, ...]) -> tuple[sp.ImmutableMatrix, sp.Expr]:
    sin_phase = sp.Symbol("sin_phase", real=True)
    tmp_omega = sp.Symbol("omega_linear", real=True)
    sub: dict[sp.Symbol, sp.Expr] = {}
    for i in range(n):
        for j in range(n):
            sub[g[i, j]] = -a[j] * k[i] * sin_phase
    for j in range(n):
        sub[v[j]] = a[j] * tmp_omega * sin_phase
    averaged = period_average(sp.expand(expr.xreplace(sub)))
    return sp.ImmutableMatrix([[sp.simplify(sp.diff(averaged, ai, aj)) for aj in a] for ai in a]), averaged


def coefficient_ordering(terms: tuple[Term, ...]) -> tuple[sp.Symbol, ...]:
    coeffs: set[sp.Symbol] = set()
    known = {rho_br, mu_R, B_comp, mu_br, beta, s, s_rho, c_s0}
    for term in terms:
        coeffs.update(sym for sym in term.factor.free_symbols if sym in known)
    return tuple(sorted(coeffs, key=lambda sym: sym.name))


def package_build(package: str, n: int, q9_data: dict[str, object]) -> PackageBuild:
    g, v = derivative_placeholders(n)
    densities = stiffness_densities(g)
    pd_density = q9_data["PD_DENSITY_PLACEHOLDER"]
    kinetic_iso = (Term(rho_br / 2, "T_ISO_SUM", sum(vj ** 2 for vj in v)),)
    kinetic_aniso = (
        Term(rho_br / 2, "T_ANISO_TRANSVERSE_SUM", sum(v[j] ** 2 for j in range(1, n))),
        Term(s_rho * rho_br / 2, "T_ANISO_AXIS", v[0] ** 2),
    )
    stiffness_map = {
        "MAIN": (
            Term(mu_R / 2, "S_curl", densities["S_curl"]),
            Term(B_comp / 2, "S_div", densities["S_div"]),
        ),
        "XFORM_CURLONLY": (Term(mu_R / 2, "S_curl", densities["S_curl"]),),
        "XFORM_DIVONLY": (Term(B_comp / 2, "S_div", densities["S_div"]),),
        "XFORM_TRACELESS": (
            Term(mu_R / 2, "S_curl", densities["S_curl"]),
            Term(mu_br, "S_symtl", densities["S_symtl"]),
        ),
        "XFORM_EXTRA": (
            Term(mu_R / 2, "S_curl", densities["S_curl"]),
            Term(B_comp / 2, "S_div", densities["S_div"]),
            Term(beta / 2, "P_D", pd_density),
        ),
        "XCOEF_BSCALE": (
            Term(mu_R / 2, "S_curl", densities["S_curl"]),
            Term(s * B_comp / 2, "S_div", densities["S_div"]),
        ),
        "XCOEF_BSIGN": (
            Term(mu_R / 2, "S_curl", densities["S_curl"]),
            Term(-B_comp / 2, "S_div", densities["S_div"]),
        ),
        "XKIN_ANISO": (
            Term(mu_R / 2, "S_curl", densities["S_curl"]),
            Term(B_comp / 2, "S_div", densities["S_div"]),
        ),
    }
    kinetic_terms = kinetic_aniso if package == "XKIN_ANISO" else kinetic_iso
    stiffness_terms = stiffness_map[package]
    all_l_terms = kinetic_terms + tuple(Term(-term.factor, term.density_name, term.density) for term in stiffness_terms)
    return PackageBuild(
        kinetic_terms=kinetic_terms,
        stiffness_terms=stiffness_terms,
        coefficient_ordering=coefficient_ordering(all_l_terms),
        lagrangian=sp.expand(sum(term.expression for term in kinetic_terms) - sum(term.expression for term in stiffness_terms)),
    )


def q9_vector(poly: sp.Expr, variables: tuple[sp.Symbol, ...],
              pairs: tuple[tuple[int, int], ...]) -> list[sp.Expr]:
    pair_index = {pair: idx for idx, pair in enumerate(pairs)}
    row = [sp.Integer(0)] * len(pairs)
    p = sp.Poly(sp.expand(poly), *variables)
    for exponents, coeff in p.terms():
        indices: list[int] = []
        for idx, exponent in enumerate(exponents):
            indices.extend([idx] * int(exponent))
        if len(indices) == 2:
            row[pair_index[tuple(sorted(indices))]] = coeff
    return row


def q9_row_to_poly(row: sp.Matrix | list[object], monomials: tuple[sp.Expr, ...]) -> sp.Expr:
    entries = list(row)
    return sp.expand(sum(entries[i] * monomials[i] for i in range(len(monomials))))


def compute_q9(n: int) -> dict[str, object]:
    variables = QG_ALL[: n * n]
    pairs = tuple(combinations_with_replacement(range(n * n), 2))
    monomials = tuple(variables[p] * variables[q] for p, q in pairs)
    qg = sp.ImmutableMatrix(n, n, lambda i, j: variables[i * n + j])

    lie_equations: list[list[sp.Expr]] = []
    for a_idx in range(n):
        for b_idx in range(a_idx + 1, n):
            generator = sp.zeros(n)
            generator[a_idx, b_idx] = 1
            generator[b_idx, a_idx] = -1
            delta = generator * qg - qg * generator
            delta_vars = tuple(delta[i, j] for i in range(n) for j in range(n))
            for p, q in pairs:
                lie_equations.append(q9_vector(delta_vars[p] * variables[q] + variables[p] * delta_vars[q], variables, pairs))

    lie_matrix = matrix_from_rows(lie_equations, len(monomials))
    v1_null = lie_matrix.nullspace()
    v1_basis = matrix_from_rows([list(vec) for vec in v1_null], len(monomials)).rref()[0]
    v1_basis = sp.ImmutableMatrix(v1_basis)
    v1_pivots = v1_basis.rref()[1]

    r0 = sp.diag(-1, *([1] * (n - 1)))
    reflected_g = r0 * qg * r0.T
    reflected_vars = tuple(reflected_g[i, j] for i in range(n) for j in range(n))
    reflection_rows = [
        q9_vector(reflected_vars[p] * reflected_vars[q] - variables[p] * variables[q], variables, pairs)
        for p, q in pairs
    ]
    v2_matrix = sp.ImmutableMatrix(list(lie_equations) + reflection_rows)
    v2_null = v2_matrix.nullspace()
    v2_basis = matrix_from_rows([list(vec) for vec in v2_null], len(monomials)).rref()[0]
    v2_basis = sp.ImmutableMatrix(v2_basis)

    v1_polys = tuple(q9_row_to_poly(v1_basis.row(i), monomials) for i in range(v1_basis.rows))
    v4_reflected = tuple(
        sp.expand(poly.xreplace({variables[i]: reflected_vars[i] for i in range(n * n)}))
        for poly in v1_polys
    )
    v4_residual = tuple(sp.expand(v4_reflected[i] - v1_polys[i]) for i in range(len(v1_polys)))
    v4_sum = tuple(sp.expand(v4_reflected[i] + v1_polys[i]) for i in range(len(v1_polys)))

    operator_cols = []
    for reflected_poly in v4_reflected:
        reflected_row = q9_vector(reflected_poly, variables, pairs)
        operator_cols.append([reflected_row[pivot] for pivot in v1_pivots])
    if operator_cols:
        v6_operator = sp.ImmutableMatrix(operator_cols).T
        minus_null = (v6_operator + sp.eye(v6_operator.rows)).nullspace()
        v6_rows = []
        for coord in minus_null:
            row = (coord.T * v1_basis)
            v6_rows.append(list(row))
        v6_basis = matrix_from_rows(v6_rows, len(monomials)).rref()[0]
    else:
        v6_operator = sp.ImmutableMatrix.zeros(0, 0)
        v6_basis = sp.ImmutableMatrix.zeros(0, len(monomials))
    v6_basis = sp.ImmutableMatrix(v6_basis)

    g_placeholder, _ = derivative_placeholders(n)
    qg_to_placeholder = {variables[i * n + j]: g_placeholder[i, j] for i in range(n) for j in range(n)}
    pd_poly = sp.expand(sum(q9_row_to_poly(v6_basis.row(i), monomials) for i in range(v6_basis.rows)))
    pd_density = sp.expand(pd_poly.xreplace(qg_to_placeholder))

    return {
        "MONOMIAL_ORDERING": sp.Tuple(*monomials),
        "V1_BASIS": v1_basis,
        "V1_DIM": sp.Integer(v1_basis.rows),
        "V2_BASIS": v2_basis,
        "V2_DIM": sp.Integer(v2_basis.rows),
        "V3_DIFFERENCE": sp.Integer(v1_basis.rows - v2_basis.rows),
        "R0_MATRIX": sp.ImmutableMatrix(r0),
        "R0_ORTHOGONALITY_RESIDUAL": sp.ImmutableMatrix(r0.T * r0 - sp.eye(n)),
        "R0_DETERMINANT": sp.det(r0),
        "V4_REFLECTED": sp.Tuple(*v4_reflected),
        "V4_RESIDUAL": sp.Tuple(*v4_residual),
        "V4_SUM": sp.Tuple(*v4_sum),
        "V5_EULER_LAGRANGE": None,
        "V6_OPERATOR": v6_operator,
        "V6_BASIS": v6_basis,
        "V6_DIM": sp.Integer(v6_basis.rows),
        "V7_RESIDUAL": sp.Integer(v6_basis.rows - (v1_basis.rows - v2_basis.rows)),
        "PD_DENSITY_PLACEHOLDER": pd_density,
        "PD_POLY": pd_poly,
        "V1_POLYS": v1_polys,
        "QG_VARIABLES": variables,
    }


def q9_v5(n: int, q9_data: dict[str, object]) -> sp.Tuple:
    g, v = derivative_placeholders(n)
    qg_vars = q9_data["QG_VARIABLES"]
    sub = {qg_vars[i * n + j]: g[i, j] for i in range(n) for j in range(n)}
    out = []
    for poly in q9_data["V1_POLYS"]:
        density = sp.expand(poly.xreplace(sub))
        out.append(sp.Tuple(*euler_lagrange_from_placeholders(density, n, g, v)))
    return sp.Tuple(*out)


def emit_q9(package: str, n: int, q9_data: dict[str, object]) -> None:
    if q9_data["V5_EULER_LAGRANGE"] is None:
        q9_data["V5_EULER_LAGRANGE"] = q9_v5(n, q9_data)
    for quantity in (
        "MONOMIAL_ORDERING", "V1_BASIS", "V1_DIM", "V2_BASIS", "V2_DIM", "V3_DIFFERENCE",
        "R0_MATRIX", "R0_ORTHOGONALITY_RESIDUAL", "R0_DETERMINANT", "V4_REFLECTED",
        "V4_RESIDUAL", "V4_SUM", "V5_EULER_LAGRANGE", "V6_OPERATOR", "V6_BASIS",
        "V6_DIM", "V7_RESIDUAL",
    ):
        emit_quantity(package, n, quantity, q9_data[quantity])


def assumptions_for(package: str, n: int, include_bulk: bool = False) -> sp.Tuple:
    _, k, a = symbols_for_dimension(n)
    assumptions: list[object] = [
        sp.Gt(rho_br, 0, evaluate=False),
        sp.Gt(mu_R, 0, evaluate=False),
        sp.Gt(B_comp, 0, evaluate=False),
        sp.Gt(sum(ki ** 2 for ki in k), 0, evaluate=False),
        sp.Q.positive(D),
        sp.Q.integer(D),
    ]
    assumptions.extend(sp.Q.real(ki) for ki in k)
    assumptions.extend(sp.Q.real(ai) for ai in a)
    if package == "XFORM_TRACELESS":
        assumptions.append(sp.Gt(mu_br, 0, evaluate=False))
    if package == "XFORM_EXTRA":
        assumptions.append(sp.Q.real(beta))
    if package == "XCOEF_BSCALE":
        assumptions.extend([sp.Gt(s, 0, evaluate=False), sp.Ne(s, 1, evaluate=False)])
    if package == "XKIN_ANISO":
        assumptions.extend([sp.Gt(s_rho, 0, evaluate=False), sp.Ne(s_rho, 1, evaluate=False)])
    if include_bulk:
        assumptions.extend([sp.Q.real(bulk_amplitude), sp.Gt(c_s0, 0, evaluate=False), sp.Q.real(c_s0)])
    return sp.Tuple(*assumptions)


def assumption_text_inventory(package: str, n: int) -> sp.Tuple:
    qualitative = (
        Str("u is in-plane only"),
        Str("background brane at rest"),
        Str("no dissipation"),
        Str("linear quadratic response"),
        Str("frozen wall width"),
        Str("bulk scalar sound mode only"),
    )
    return sp.Tuple(assumptions_for(package, n), qualitative)


def symbolic_truth_status(test_object: object, operands: object) -> sp.Tuple:
    if test_object is True or test_object == sp.true:
        token = "PROVED_TRUE"
    elif test_object is False or test_object == sp.false:
        token = "PROVED_FALSE"
    else:
        token = "UNDECIDED"
    return sp.Tuple(sp.Tuple(Str("STATUS_TOKEN"), Str(token)), sp.Tuple(Str("TEST_OBJECT"), casify(test_object)), sp.Tuple(Str("OPERANDS"), casify(operands)))


def sign_payload(expr: sp.Expr, assumptions: sp.Tuple) -> sp.Tuple:
    expr = sp.factor(sp.simplify(expr))
    if expr == 0:
        token = "ZERO"
    elif expr.is_positive is True:
        token = "POSITIVE"
    elif expr.is_negative is True:
        token = "NEGATIVE"
    else:
        token = "UNDECIDED"
    return sp.Tuple(sp.Tuple(Str("SIGN_TOKEN"), Str(token)), sp.Tuple(Str("OPERAND"), casify(expr)))


def relation_residual(rel: object) -> sp.Expr:
    if isinstance(rel, sp.Equality):
        return sp.expand(rel.lhs - rel.rhs)
    return sp.sympify(rel)


def solve_result_to_cas(solution: object, variables: tuple[sp.Symbol, ...]) -> object:
    if isinstance(solution, list):
        branches = []
        for branch in solution:
            if isinstance(branch, dict):
                branches.append(sp.Tuple(*[sp.Tuple(var, casify(branch[var])) for var in variables if var in branch]))
            else:
                branches.append(casify(branch))
        return sp.Tuple(*branches)
    return casify(solution)


def exposed_branches(solution: object, variables: tuple[sp.Symbol, ...]) -> list[dict[sp.Symbol, object]] | None:
    if not isinstance(solution, list):
        return None
    branches = []
    for branch in solution:
        if not isinstance(branch, dict):
            return None
        branches.append(branch)
    return branches


def evaluate_premise(premise: object, branch: dict[sp.Symbol, object]) -> bool | None:
    substituted = premise.subs(branch) if hasattr(premise, "subs") else premise
    if substituted is True or substituted == sp.true:
        return True
    if substituted is False or substituted == sp.false:
        return False
    if isinstance(substituted, (sp.StrictGreaterThan, sp.StrictLessThan, sp.GreaterThan, sp.LessThan, sp.Equality, sp.Unequality)):
        simplified = sp.simplify(substituted)
        if simplified == sp.true:
            return True
        if simplified == sp.false:
            return False
        return None
    if getattr(substituted, "func", None).__name__ == "AppliedPredicate":
        argument = substituted.arguments[0] if getattr(substituted, "arguments", ()) else None
        predicate_name = getattr(substituted.function, "name", "")
        attribute_name = {
            "real": "is_real",
            "positive": "is_positive",
            "integer": "is_integer",
        }.get(predicate_name)
        if attribute_name is not None:
            decision = getattr(argument, attribute_name, None)
            if decision is True:
                return True
            if decision is False:
                return False
        return None
    return None


def admissibility_entry(branch: dict[sp.Symbol, object], assumptions: sp.Tuple) -> sp.Tuple:
    checks = [evaluate_premise(premise, branch) for premise in assumptions]
    test_terms = []
    for premise, check in zip(assumptions, checks):
        if check is True:
            test_terms.append(sp.true)
        elif check is False:
            test_terms.append(sp.false)
        else:
            test_terms.append(premise.subs(branch) if hasattr(premise, "subs") else premise)
    test_object = sp.And(*test_terms)
    if test_object == sp.false:
        token = "EXCLUDED"
    elif test_object == sp.true:
        token = "ADMISSIBLE"
    else:
        token = "UNDECIDED"
    return sp.Tuple(
        sp.Tuple(Str("BRANCH"), solve_result_to_cas([branch], tuple(branch.keys()))),
        sp.Tuple(Str("STATUS_TOKEN"), Str(token)),
        sp.Tuple(Str("TEST_OBJECT"), casify(test_object)),
        sp.Tuple(Str("OPERANDS"), sp.Tuple(casify(branch), assumptions)),
    )


def candidate_witness(variables: tuple[sp.Symbol, ...]) -> dict[sp.Symbol, sp.Expr]:
    values = {}
    for var in variables:
        if var in K_ALL:
            values[var] = sp.Integer(1) if var == variables[0] else sp.Integer(0)
        elif var in (rho_br, mu_R, B_comp, mu_br, c_s0):
            values[var] = sp.Integer(1)
        elif var == s or var == s_rho:
            values[var] = sp.Integer(2)
        elif var == beta:
            values[var] = sp.Integer(0)
        elif var == kwSquared:
            values[var] = sp.Integer(0)
        else:
            values[var] = sp.Integer(1)
    return values


def witness_valid(equations: tuple[object, ...], assumptions: sp.Tuple, witness: dict[sp.Symbol, sp.Expr]) -> bool:
    for eq in equations:
        residual = relation_residual(eq).subs(witness)
        if sp.simplify(residual) != 0:
            return False
    for premise in assumptions:
        result = evaluate_premise(premise, witness)
        if result is False:
            return False
    return True


def locus_conditionset(equations: tuple[object, ...], variables: tuple[sp.Symbol, ...]) -> sp.Set:
    condition = sp.And(*[sp.simplify(equation) for equation in equations]) if equations else sp.true
    if not variables:
        return sp.ConditionSet(sp.Tuple(), condition, sp.FiniteSet(sp.Tuple()))
    if len(variables) == 1:
        domain = sp.S.Complexes
        symbol: object = variables[0]
    else:
        domain = sp.ProductSet(*([sp.S.Complexes] * len(variables)))
        symbol = sp.Tuple(*variables)
    return sp.ConditionSet(symbol, condition, domain)


def canonical_locus(residuals: tuple[sp.Expr, ...], variables: tuple[sp.Symbol, ...],
                    *, package: str | None = None, n: int | None = None,
                    base_quantity: str = "CANONICAL_LOCUS", root: int | None = None,
                    stratum: int | None = None) -> object:
    if not residuals:
        return sp.Tuple()
    try:
        polys = [sp.Poly(sp.together(residual), *variables) for residual in residuals]
    except (sp.PolynomialError, ValueError):
        return NOT_APPLICABLE
    try:
        basis = sp.groebner([poly.as_expr() for poly in polys], *variables, order="lex")
        return sp.Tuple(*[sp.factor(expr) for expr in basis.polys])
    except Exception as exc:  # exact Groebner can be unavailable for some coefficient domains
        record_issue(
            f"{base_quantity}: canonical locus unavailable for variables {[v.name for v in variables]} "
            f"after {type(exc).__name__}; emitted NOT_APPLICABLE",
            package=package, n=n, root=root, stratum=stratum,
        )
        return NOT_APPLICABLE


def _single_quadratic_radical_locus(
    residuals: tuple[sp.Expr, ...], variables: tuple[sp.Symbol, ...],
) -> bool:
    if not residuals or not variables:
        return False
    numerators = []
    radical_nodes = []
    for residual in residuals:
        try:
            numerator, denominator = sp.fraction(sp.together(residual))
        except Exception:
            return False
        if quadratic_radical_nodes(denominator):
            return False
        numerators.append(numerator)
        radical_nodes.extend(quadratic_radical_nodes(numerator))

    radical_bases = {sp.srepr(node.base): node.base for node in radical_nodes}
    if len(radical_bases) != 1:
        return False
    radicand = next(iter(radical_bases.values()))
    radicand_variables = radicand.free_symbols & set(variables)
    if (
        quadratic_radical_nodes(radicand)
        or not radicand_variables
        or len(set(variables) - radicand_variables) != 2
    ):
        return False

    radical_symbol = sp.Dummy("quadratic_radical")
    replacements = {
        node: node.base ** ((int(node.exp.p) - 1) // 2) * radical_symbol
        for node in set(radical_nodes)
    }
    for numerator in numerators:
        lifted = numerator.xreplace(replacements)
        symbols = tuple(sorted(lifted.free_symbols - {radical_symbol}, key=sp.default_sort_key))
        if any(symbol not in DECLARED_SYMBOLS for symbol in symbols):
            return False
        try:
            sp.Poly(lifted, *(symbols + (radical_symbol,)), domain=sp.QQ)
        except (sp.PolynomialError, ValueError):
            return False
    return True


def _solve_single_quadratic_radical_locus(
    residuals: tuple[sp.Expr, ...], variables: tuple[sp.Symbol, ...],
) -> list[dict[sp.Symbol, sp.Expr]]:
    from sympy.solvers.solvers import _invert, _simple_dens, _vsolve, checksol

    dens = set()
    failed = []
    for residual in residuals:
        dens.update(_simple_dens(residual, variables))
        independent, dependent = _invert(residual, *variables)
        failed.append((dependent - independent).as_numer_denom()[0])

    result: list[dict[sp.Symbol, sp.Expr]] = [{}]
    legal = set(variables)

    def available_symbols(expr: sp.Expr, *, sort: bool = False) -> object:
        available = expr.free_symbols & legal
        if not sort:
            return available

        def key(symbol: sp.Symbol) -> object:
            poly = expr.as_poly(symbol)
            if poly is None:
                complexity = (sp.S.Infinity, sp.S.Infinity, sp.S.Infinity)
            else:
                coefficient_symbols = poly.LC().free_symbols
                complexity = (
                    poly.degree(),
                    len(coefficient_symbols & available),
                    len(coefficient_symbols),
                )
            return complexity + (sp.default_sort_key(symbol),)

        return sorted(available, key=key)

    check_flags = {"simplify": False, "manual": True}
    check_symbol = sp.Dummy()
    for equation in sp.ordered(failed, lambda expr: len(available_symbols(expr))):
        new_result: list[dict[sp.Symbol, sp.Expr]] = []
        bad_results: list[dict[sp.Symbol, sp.Expr]] = []
        hit = False
        for branch in result:
            solved_symbols = set()
            equation_substituted = equation.subs(branch)
            if branch:
                checked = checksol(
                    check_symbol, check_symbol, equation_substituted, minimal=True,
                )
                if checked is not None:
                    (new_result if checked else bad_results).append(branch)
                    continue
            candidates = available_symbols(equation_substituted, sort=True)
            if not candidates:
                if branch:
                    new_result.append(branch)
                break
            for symbol in candidates:
                try:
                    solutions = _vsolve(equation_substituted, symbol, **check_flags)
                except NotImplementedError:
                    continue
                for value in solutions:
                    if solved_symbols and any(
                        symbol in value.free_symbols for symbol in solved_symbols
                    ):
                        continue
                    expanded_branch = {
                        known_symbol: known_value.subs(symbol, value)
                        for known_symbol, known_value in branch.items()
                    }
                    expanded_branch[symbol] = value
                    expanded_items = set(expanded_branch.items())
                    for known_branch in new_result:
                        if len(known_branch) < len(expanded_items):
                            updated_items = {
                                (known_symbol, known_value.xreplace(expanded_branch))
                                for known_symbol, known_value in known_branch.items()
                            }
                            if not updated_items - expanded_items:
                                break
                    else:
                        new_result.append(expanded_branch)
                hit = True
                solved_symbols.add(symbol)
                break
            if not hit:
                raise NotImplementedError(f"could not solve {equation_substituted}")
        else:
            result = new_result
            for bad_result in bad_results:
                if bad_result in result:
                    result.remove(bad_result)

    result = [
        branch for branch in result
        if not any(checksol(denominator, branch, **check_flags) for denominator in dens)
    ]
    result = [
        branch for branch in result
        if not any(
            checksol(residual, branch, **check_flags) is False
            for residual in residuals
        )
    ]
    result = [branch for branch in result if branch]
    result = [
        {symbol: branch[symbol] for symbol in sp.ordered(branch)}
        for branch in result
    ]
    result.sort(key=sp.default_sort_key)
    return result


def emit_locus(package: str, n: int, base_quantity: str, residuals: tuple[sp.Expr, ...],
               variables: tuple[sp.Symbol, ...], assumptions: sp.Tuple,
               *, root: int | None = None, stratum: int | None = None,
               export: bool = True) -> dict[str, object]:
    equations = tuple(sp.Eq(sp.factor(residual), 0, evaluate=False) for residual in residuals)
    residual_exprs = tuple(relation_residual(eq) for eq in equations)
    canonical_payload = canonical_locus(
        residual_exprs, variables,
        package=package, n=n, base_quantity=base_quantity, root=root, stratum=stratum,
    )
    emit_quantity(package, n, f"{base_quantity}_EQUATIONS", sp.Tuple(*equations), root=root, stratum=stratum, export=export)
    try:
        solve_input = list(residual_exprs)
        solve_kwargs = {"dict": True, "simplify": False, "manual": True}
        if isinstance(canonical_payload, sp.Tuple) and canonical_payload:
            solve_input = [
                expr.as_expr() if isinstance(expr, sp.Poly) else expr
                for expr in canonical_payload
            ]
            solve_kwargs = {"dict": True, "simplify": False}
        if (
            canonical_payload == NOT_APPLICABLE
            and _single_quadratic_radical_locus(residual_exprs, variables)
        ):
            solution = _solve_single_quadratic_radical_locus(residual_exprs, variables)
        else:
            solution = sp.solve(solve_input, list(variables), **solve_kwargs)
    except Exception as exc:
        solution = Str(f"SOLVE_UNAVAILABLE_{type(exc).__name__}")
        record_issue(
            f"{base_quantity}: solve attempt raised {type(exc).__name__}: {exc!r}",
            package=package, n=n, root=root, stratum=stratum,
        )
    solution_payload = solve_result_to_cas(solution, variables)
    emit_quantity(package, n, f"{base_quantity}_SOLUTION", solution_payload, root=root, stratum=stratum, export=export)

    zero_statuses = [algebraic_zero_test(relation_residual(eq)) for eq in equations]
    if not equations:
        identical_object: object = sp.true
    elif all(status is True for status in zero_statuses):
        identical_object = sp.true
    elif any(status is False for status in zero_statuses):
        identical_object = sp.false
    else:
        identical_object = sp.And(*[
            sp.Eq(sp.simplify(relation_residual(eq)), 0, evaluate=False)
            for eq in equations
        ])
    ident_payload = symbolic_truth_status(identical_object, sp.Tuple(*equations))
    emit_quantity(package, n, f"{base_quantity}_IDENTICALLY_SATISFIED", ident_payload, root=root, stratum=stratum, export=export)

    equation_defined_set = locus_conditionset(equations, variables)
    inconsistent_object = sp.Eq(equation_defined_set, sp.EmptySet)
    inconsistent_operands = sp.Tuple(sp.Tuple(*equations), sp.Tuple(*variables), equation_defined_set)
    inconsistent_payload = symbolic_truth_status(inconsistent_object, inconsistent_operands)
    emit_quantity(package, n, f"{base_quantity}_INCONSISTENT", inconsistent_payload, root=root, stratum=stratum, export=export)

    branches = exposed_branches(solution, variables)
    admissible_entries: list[object] = []
    if branches is None:
        admissible_entries.append(
            sp.Tuple(
                sp.Tuple(Str("BRANCH"), solution_payload),
                sp.Tuple(Str("STATUS_TOKEN"), Str("UNDECIDED")),
                sp.Tuple(Str("TEST_OBJECT"), solution_payload),
                sp.Tuple(Str("OPERANDS"), sp.Tuple(solution_payload, assumptions)),
            )
        )
    else:
        for branch in branches:
            admissible_entries.append(admissibility_entry(branch, assumptions))
    real_admissible = sp.Tuple(*admissible_entries)
    emit_quantity(package, n, f"{base_quantity}_REAL_ADMISSIBLE", real_admissible, root=root, stratum=stratum, export=export)

    emit_quantity(
        package, n, f"{base_quantity}_CANONICAL_LOCUS", canonical_payload,
        root=root, stratum=stratum, export=export,
    )

    witness = candidate_witness(variables)
    if inconsistent_object == sp.true:
        real_status = Str("PROVED_EMPTY")
        real_witness = NOT_APPLICABLE
    elif witness_valid(equations, assumptions, witness):
        real_status = Str("PROVED_NONEMPTY")
        real_witness = sp.Tuple(*[sp.Tuple(var, witness[var]) for var in variables])
    else:
        real_status = Str("UNDECIDED")
        real_witness = NOT_APPLICABLE
    emit_quantity(package, n, f"{base_quantity}_REAL_STATUS", real_status, root=root, stratum=stratum, export=export)
    emit_quantity(package, n, f"{base_quantity}_REAL_WITNESS", real_witness, root=root, stratum=stratum, export=export)
    emit_quantity(package, n, f"{base_quantity}_REAL_STATUS_OPERANDS", sp.Tuple(sp.Tuple(*equations), sp.Tuple(*variables), assumptions), root=root, stratum=stratum, export=export)
    return {
        "equations": equations,
        "solution": solution,
        "real_admissible": real_admissible,
        "base_quantity": base_quantity,
        "variables": variables,
    }


def root_multiplicity_pairs(det_expr: sp.Expr) -> tuple[RootData, ...]:
    poly = sp.Poly(sp.factor(det_expr), omegaSquared)
    roots = sp.roots(poly.as_expr(), omegaSquared)
    if sum(roots.values()) != poly.degree():
        solved = sp.solve(sp.Eq(poly.as_expr(), 0, evaluate=False), omegaSquared)
        missing = [root for root in solved if root not in roots]
        for root in missing:
            roots[root] = 1
        if sum(roots.values()) != poly.degree():
            record_issue("root multiplicity solve did not account for every polynomial degree")
    return tuple(RootData(sp.factor(root), int(mult)) for root, mult in roots.items())


def distinct_roots(pairs: tuple[RootData, ...]) -> tuple[sp.Expr, ...]:
    roots: list[sp.Expr] = []
    for pair in pairs:
        if not any(sp.simplify(pair.value - existing) == 0 for existing in roots):
            roots.append(pair.value)
    return tuple(roots)


def reduced_expr(expr: sp.Expr) -> sp.Expr:
    return sp.factor(sp.cancel(expr))


def reduced_matrix(matrix: sp.MatrixBase) -> sp.Matrix:
    return sp.Matrix(matrix).applyfunc(reduced_expr)


def quadratic_radical_nodes(expr: sp.Expr) -> tuple[sp.Pow, ...]:
    return tuple(
        node
        for node in sp.preorder_traversal(expr)
        if isinstance(node, sp.Pow)
        and node.exp.is_Rational
        and int(node.exp.q) == 2
        and int(node.exp.p) % 2 != 0
    )


def _declared_polynomial_zero_test(expr: sp.Expr) -> bool | None:
    symbols = tuple(sorted(expr.free_symbols, key=lambda sym: sym.name))
    if any(symbol not in DECLARED_SYMBOLS for symbol in symbols):
        return None
    try:
        reduced = reduced_expr(expr)
        if not symbols:
            return reduced == 0 if isinstance(reduced, sp.Rational) else None
        poly = sp.Poly(reduced, *symbols, domain=sp.QQ)
    except Exception:
        return None
    return True if poly.is_zero else False


def _joint_zero_test_assumptions() -> sp.Expr:
    facts: list[sp.Expr] = [
        sp.Q.positive(rho_br),
        sp.Q.positive(mu_R),
        sp.Q.positive(B_comp),
        sp.Q.positive(D),
        sp.Q.integer(D),
    ]
    if CURRENT_CELL is not None:
        package, n = CURRENT_CELL
        facts.append(sp.Q.positive(sum(ki ** 2 for ki in K_ALL[:n])))
        if package == "XFORM_TRACELESS":
            facts.append(sp.Q.positive(mu_br))
        if package == "XCOEF_BSCALE":
            facts.extend((sp.Q.positive(s), sp.Q.nonzero(s - 1)))
        if package == "XKIN_ANISO":
            facts.extend((sp.Q.positive(s_rho), sp.Q.nonzero(s_rho - 1)))
    return sp.And(*facts)


def _linear_radical_zero_test(a: sp.Expr, b: sp.Expr, radicand: sp.Expr) -> bool | None:
    statuses = tuple(_declared_polynomial_zero_test(expr) for expr in (a, b, radicand))
    if any(status is None for status in statuses):
        return None
    a_zero, b_zero, radicand_zero = statuses
    if b_zero or radicand_zero:
        return a_zero
    norm = reduced_expr(a ** 2 - b ** 2 * radicand)
    norm_status = _declared_polynomial_zero_test(norm)
    if norm_status is False:
        return False
    if norm_status is None:
        return None
    quotient = reduced_expr(a / b)
    if quotient == 0 or quotient.is_zero is True:
        return True
    if quotient.is_positive is True or quotient.is_nonnegative is True:
        return False
    if quotient.is_nonpositive is True:
        return True
    assumptions = _joint_zero_test_assumptions()
    try:
        if sp.ask(sp.Q.positive(quotient), assumptions) is True:
            return False
        if sp.ask(sp.Q.nonnegative(quotient), assumptions) is True:
            return False
        if sp.ask(sp.Q.nonpositive(quotient), assumptions) is True:
            return True
    except Exception:
        pass
    return None


def algebraic_zero_test(entry: sp.Expr) -> bool | None:
    if entry == 0 or entry.is_zero is True:
        return True
    radical_nodes = quadratic_radical_nodes(entry)
    radical_bases = {sp.srepr(node.base): node.base for node in radical_nodes}
    if len(radical_bases) > 1:
        return None
    if not radical_bases:
        try:
            numerator, denominator = sp.fraction(sp.together(entry))
        except Exception:
            return None
        denominator_status = _declared_polynomial_zero_test(denominator)
        if denominator_status is not False:
            return None
        return _declared_polynomial_zero_test(numerator)

    radicand = next(iter(radical_bases.values()))
    if quadratic_radical_nodes(radicand):
        return None
    radical_symbol = sp.Dummy("quadratic_radical")
    replacements = {
        node: node.base ** ((int(node.exp.p) - 1) // 2) * radical_symbol
        for node in set(radical_nodes)
    }
    try:
        lifted = sp.together(entry.xreplace(replacements))
        numerator, denominator = sp.fraction(lifted)
        modulus = sp.Poly(radical_symbol ** 2 - radicand, radical_symbol, domain="EX")
        numerator_poly = sp.Poly(numerator, radical_symbol, domain="EX").rem(modulus)
        denominator_poly = sp.Poly(denominator, radical_symbol, domain="EX").rem(modulus)
    except Exception:
        return None
    numerator_status = _linear_radical_zero_test(
        numerator_poly.nth(0), numerator_poly.nth(1), radicand,
    )
    denominator_status = _linear_radical_zero_test(
        denominator_poly.nth(0), denominator_poly.nth(1), radicand,
    )
    return numerator_status if denominator_status is False else None


def exact_domain_matrix(matrix: sp.MatrixBase, *, already_reduced: bool = False) -> DomainMatrix:
    source = sp.Matrix(matrix) if already_reduced else reduced_matrix(matrix)
    return DomainMatrix.from_Matrix(source).to_field()


def matrix_rank(matrix: sp.MatrixBase) -> int:
    try:
        return int(exact_domain_matrix(matrix).rank())
    except MemoryError:
        raise
    except Exception as exc:
        record_issue(f"matrix rank exact-domain route unavailable after {type(exc).__name__}; used algebraic zero fallback")
        simplified = reduced_matrix(matrix)
        return int(simplified.rank(iszerofunc=algebraic_zero_test, simplify=True))


def zero_expr(expr: sp.Expr) -> bool:
    return algebraic_zero_test(expr) is True


def _single_quadratic_radical_determinant(matrix: sp.MatrixBase) -> sp.Expr | None:
    radical_nodes = tuple(
        node for entry in matrix for node in quadratic_radical_nodes(entry)
    )
    radical_bases = {sp.srepr(node.base): node.base for node in radical_nodes}
    if (
        len(radical_bases) != 1
        or any(node.exp < 0 for node in radical_nodes)
        or any(quadratic_radical_nodes(node.base) for node in radical_nodes)
    ):
        return None

    radicand = next(iter(radical_bases.values()))
    radical_symbol = sp.Dummy("quadratic_radical")
    replacements = {
        node: node.base ** ((int(node.exp.p) - 1) // 2) * radical_symbol
        for node in set(radical_nodes)
    }
    lifted_fractions: list[list[tuple[sp.Expr, sp.Expr]]] = []
    for i in range(matrix.rows):
        lifted_row = []
        for j in range(matrix.cols):
            lifted_entry = sp.together(matrix[i, j].xreplace(replacements))
            numerator, denominator = sp.fraction(lifted_entry)
            if radical_symbol in denominator.free_symbols:
                return None
            try:
                sp.Poly(numerator, radical_symbol, domain="EX")
            except sp.PolynomialError:
                return None
            lifted_row.append((numerator, denominator))
        lifted_fractions.append(lifted_row)

    numerator_matrix = sp.zeros(matrix.rows, matrix.cols)
    denominator_product = sp.Integer(1)
    for i, lifted_row in enumerate(lifted_fractions):
        row_denominator = sp.Integer(1)
        for _, denominator in lifted_row:
            row_denominator = sp.lcm(row_denominator, denominator)
        for j, (numerator, denominator) in enumerate(lifted_row):
            numerator_matrix[i, j] = sp.cancel(
                numerator * row_denominator / denominator
            )
        denominator_product *= row_denominator

    lifted_determinant = exact_domain_matrix(
        numerator_matrix, already_reduced=True,
    ).det().as_expr()
    modulus = sp.Poly(radical_symbol ** 2 - radicand, radical_symbol, domain="EX")
    reduced_determinant = sp.Poly(
        lifted_determinant, radical_symbol, domain="EX",
    ).rem(modulus).as_expr()
    return reduced_determinant.subs(
        radical_symbol, radicand ** sp.Rational(1, 2),
    ) / denominator_product


def exact_determinant(matrix: sp.MatrixBase) -> sp.Expr:
    simplified = reduced_matrix(matrix)
    if any(quadratic_radical_nodes(entry) for entry in simplified):
        det_expr = _single_quadratic_radical_determinant(simplified)
        if det_expr is None:
            det_expr = simplified.det(method="bareiss")
    else:
        try:
            det_expr = exact_domain_matrix(simplified, already_reduced=True).det().as_expr()
        except MemoryError:
            raise
        except Exception as exc:
            record_issue(f"determinant exact-domain route unavailable after {type(exc).__name__}; used SymPy determinant fallback")
            det_expr = simplified.det()
    return reduced_expr(det_expr)


def determinant_from_live_matrix(matrix: sp.MatrixBase, k: tuple[sp.Symbol, ...]) -> sp.Expr:
    n = matrix.rows
    coupling = None
    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            if not zero_expr(matrix[i, j]):
                coupling = sp.cancel(matrix[i, j] / (k[i] * k[j]))
                break
        if coupling is not None:
            break
    if coupling is None:
        coupling = sp.Integer(0)
    diagonal = [sp.cancel(matrix[i, i] - coupling * k[i] ** 2) for i in range(n)]
    reconstructed = sp.ImmutableMatrix(
        n, n,
        lambda i, j: diagonal[i] * (1 if i == j else 0) + coupling * k[i] * k[j],
    )
    if all(zero_expr(matrix[i, j] - reconstructed[i, j]) for i in range(n) for j in range(n)):
        product = sp.prod(diagonal)
        lemma_factor = 1 + coupling * sum(k[i] ** 2 / diagonal[i] for i in range(n))
        return reduced_expr(product * lemma_factor)
    if n <= 4:
        return exact_determinant(matrix)
    record_issue(f"D{n}: determinant fallback unavailable for unmatched matrix structure")
    return Str("DETERMINANT_UNAVAILABLE_UNMATCHED_STRUCTURE")


def all_minors(matrix: sp.MatrixBase, size: int) -> tuple[sp.Expr, ...]:
    if size == 0:
        return tuple()
    minors = []
    rows = range(matrix.rows)
    cols = range(matrix.cols)
    for row_set in combinations(rows, size):
        for col_set in combinations(cols, size):
            minors.append(exact_determinant(matrix.extract(row_set, col_set)))
    return tuple(minors)


def build_nullspace_from_pivot(mat: sp.Matrix, rank: int, selected_rows: tuple[int, ...],
                               selected_cols: tuple[int, ...]) -> list[sp.ImmutableMatrix]:
    selected_minor = mat.extract(selected_rows, selected_cols)
    pivot_determinant = exact_determinant(selected_minor)
    if algebraic_zero_test(pivot_determinant) is not False:
        raise ValueError("candidate pivot determinant was not proven nonzero")
    free_cols = [col for col in range(mat.cols) if col not in selected_cols]
    basis = []
    for free_col in free_cols:
        rhs = -sp.ImmutableMatrix([mat[row, free_col] for row in selected_rows])
        pivot_values = []
        for pivot_col in range(rank):
            replaced = sp.Matrix(selected_minor)
            replaced[:, pivot_col] = rhs
            pivot_values.append(reduced_expr(exact_determinant(replaced) / pivot_determinant))
        pivot_solution = sp.ImmutableMatrix(pivot_values)
        pivot_residual = selected_minor * pivot_solution - rhs
        if any(algebraic_zero_test(reduced_expr(entry)) is not True for entry in pivot_residual):
            raise ValueError("candidate pivot solve residual was not zero")
        vector_entries = [sp.Integer(0)] * mat.cols
        for idx, col in enumerate(selected_cols):
            vector_entries[col] = pivot_solution[idx]
        vector_entries[free_col] = sp.Integer(1)
        basis.append(sp.ImmutableMatrix(mat.cols, 1, vector_entries))
    if len(basis) != mat.cols - rank:
        raise ValueError("basis count did not match matrix nullity")
    residual_entries = []
    for vec in basis:
        residual_entries.extend(reduced_expr(entry) for entry in (mat * vec))
    if any(algebraic_zero_test(entry) is False for entry in residual_entries):
        raise ValueError("candidate pivot produced a nonzero nullspace residual")
    if any(algebraic_zero_test(entry) is None for entry in residual_entries):
        raise ValueError("candidate pivot residual was undecided")
    return basis


def generic_nullspace_vectors(matrix: sp.MatrixBase, rank: int, *, root: int | None = None) -> list[sp.ImmutableMatrix]:
    mat = reduced_matrix(matrix)
    cols = mat.cols
    if rank == cols:
        return []
    if rank == 0:
        return [sp.ImmutableMatrix(cols, 1, lambda i, _j, col=col: sp.Integer(1) if i == col else sp.Integer(0)) for col in range(cols)]

    undecided_candidates: list[tuple[tuple[int, ...], tuple[int, ...]]] = []
    for row_set in combinations(range(mat.rows), rank):
        for col_set in combinations(range(cols), rank):
            minor = exact_determinant(mat.extract(row_set, col_set))
            status = algebraic_zero_test(minor)
            if status is False:
                return build_nullspace_from_pivot(mat, rank, row_set, col_set)
            if status is None:
                undecided_candidates.append((row_set, col_set))

    last_error: Exception | None = None
    for row_set, col_set in undecided_candidates:
        try:
            basis = build_nullspace_from_pivot(mat, rank, row_set, col_set)
        except Exception as exc:
            last_error = exc
            continue
        record_issue(
            "generic nullspace source minor had undecided zero status; emitted basis validated by live matrix residual",
            root=root,
        )
        return basis

    detail = f"; last undecided candidate failed after {type(last_error).__name__}" if last_error is not None else ""
    record_issue(f"generic nullspace basis unavailable: no usable source minor found{detail}", root=root)
    return []


def scale_exponent_payload(ratio: object) -> object:
    if ratio == UNDEFINED_RATIO:
        return UNDEFINED_RATIO
    ratio = sp.simplify(ratio)
    if ratio == 1:
        return sp.Integer(0)
    powers = ratio.as_powers_dict()
    exponent = powers.get(lambdaScale)
    if exponent is None:
        return NOT_A_PURE_POWER
    residual = sp.simplify(ratio / (lambdaScale ** exponent))
    if residual == 1:
        return sp.simplify(exponent)
    return NOT_A_PURE_POWER


def expression_dimension(expr: object, dim_env: dict[sp.Symbol, sp.ImmutableMatrix],
                         g_dims: dict[sp.Symbol, sp.ImmutableMatrix] | None = None,
                         v_dims: dict[sp.Symbol, sp.ImmutableMatrix] | None = None) -> sp.ImmutableMatrix | None:
    g_dims = g_dims or {}
    v_dims = v_dims or {}
    expr = casify(expr)
    if isinstance(expr, sp.MatrixBase):
        entries = [expression_dimension(entry, dim_env, g_dims, v_dims) for entry in list(expr)]
        if not entries:
            return ZERO_DIM
        if any(entry is None for entry in entries):
            return None
        return entries[0] if all(sp.simplify(entries[0] - entry) == sp.zeros(3, 1) for entry in entries[1:]) else None
    if isinstance(expr, (sp.Integer, sp.Rational, sp.Float)):
        return ZERO_DIM
    if isinstance(expr, Str):
        return ZERO_DIM
    if isinstance(expr, sp.Symbol):
        if expr in dim_env:
            return sp.ImmutableMatrix(dim_env[expr])
        if expr in g_dims:
            return sp.ImmutableMatrix(g_dims[expr])
        if expr in v_dims:
            return sp.ImmutableMatrix(v_dims[expr])
        return ZERO_DIM
    if isinstance(expr, sp.Derivative):
        base = FIELD_DIM
        for var, count in expr.variable_count:
            if var == t:
                base += count * TIME_DERIVATIVE_DIM
            else:
                base += count * SPATIAL_DERIVATIVE_DIM
        return sp.ImmutableMatrix(base)
    if isinstance(expr, sp.Pow):
        base = expression_dimension(expr.base, dim_env, g_dims, v_dims)
        if base is None:
            return None
        exponent = expr.exp
        if exponent.free_symbols:
            return None
        return sp.ImmutableMatrix([sp.simplify(exponent * entry) for entry in base])
    if isinstance(expr, sp.Mul):
        total = ZERO_DIM
        for arg in expr.args:
            dim = expression_dimension(arg, dim_env, g_dims, v_dims)
            if dim is None:
                return None
            total += dim
        return sp.ImmutableMatrix([sp.simplify(entry) for entry in total])
    if isinstance(expr, sp.Add):
        dims = [expression_dimension(arg, dim_env, g_dims, v_dims) for arg in expr.args]
        if any(dim is None for dim in dims):
            return None
        if not dims:
            return ZERO_DIM
        if all(sp.simplify(dims[0] - dim) == sp.zeros(3, 1) for dim in dims[1:]):
            return dims[0]
        return None
    return ZERO_DIM


def dimension_solve(package: str, n: int, build: PackageBuild, g: sp.ImmutableMatrix,
                    v: tuple[sp.Symbol, ...]) -> dict[str, object]:
    coeffs = build.coefficient_ordering
    dim_symbols = {
        coeff: sp.ImmutableMatrix([
            sp.Symbol(f"dim_{coeff.name}_L"),
            sp.Symbol(f"dim_{coeff.name}_T"),
            sp.Symbol(f"dim_{coeff.name}_M"),
        ])
        for coeff in coeffs
        if coeff != s
    }
    dim_env: dict[sp.Symbol, sp.ImmutableMatrix] = {
        D: ZERO_DIM,
        omegaSquared: OMEGA_SQUARED_DIM,
        lambdaScale: ZERO_DIM,
        c_s0: sp.ImmutableMatrix([1, -1, 0]),
        kwSquared: sp.ImmutableMatrix([-2, 0, 0]),
        s: ZERO_DIM,
    }
    for ki in K_ALL[:n]:
        dim_env[ki] = WAVEVECTOR_DIM
    for ai in A_ALL[:n]:
        dim_env[ai] = FIELD_DIM
    for coeff, dim_vec in dim_symbols.items():
        dim_env[coeff] = dim_vec
    g_dims = {g[i, j]: FIELD_DIM + SPATIAL_DERIVATIVE_DIM for i in range(n) for j in range(n)}
    v_dims = {v[j]: FIELD_DIM + TIME_DERIVATIVE_DIM for j in range(n)}
    equations = []
    all_l_terms = list(build.kinetic_terms) + [Term(-term.factor, term.density_name, term.density) for term in build.stiffness_terms]
    for term in all_l_terms:
        dim = expression_dimension(term.expression, dim_env, g_dims, v_dims)
        if dim is None:
            record_issue(f"dimension walker found inhomogeneous action term {term.density_name}", package=package, n=n)
            continue
        equations.extend(sp.Eq(sp.simplify(dim[i]), ENERGY_DENSITY_DIM[i], evaluate=False) for i in range(3))
    unknowns = tuple(entry for coeff in coeffs if coeff != s for entry in dim_symbols[coeff])
    linear_exprs = [eq.lhs - eq.rhs for eq in equations]
    solution = sp.linsolve(linear_exprs, unknowns) if unknowns else sp.FiniteSet(())
    if solution:
        solution_tuple = tuple(next(iter(solution))) if unknowns else tuple()
        substitutions = dict(zip(unknowns, solution_tuple))
    else:
        substitutions = {}
    coeff_dimensions = []
    for coeff in coeffs:
        if coeff == s:
            coeff_dimensions.append(sp.Tuple(coeff, ZERO_DIM))
        else:
            coeff_dimensions.append(sp.Tuple(coeff, sp.ImmutableMatrix([sp.simplify(x.subs(substitutions)) for x in dim_symbols[coeff]])))
    if unknowns:
        coeff_matrix, rhs = sp.linear_eq_to_matrix(linear_exprs, unknowns)
        equation_count = int(coeff_matrix.rank())
    else:
        equation_count = 0
    unknown_count = len(coeffs) - (1 if s in coeffs else 0)
    if equation_count > 3 * unknown_count:
        determinacy = Str("OVER_DETERMINED")
    elif equation_count == 3 * unknown_count:
        determinacy = Str("EXACTLY_DETERMINED")
    else:
        determinacy = Str("UNDER_DETERMINED")
    return {
        "COEFFICIENT_ORDERING": sp.Tuple(*coeffs),
        "DIM_COEFFICIENTS": sp.Tuple(*coeff_dimensions),
        "DIM_EQUATIONS": sp.Tuple(*equations),
        "DIM_SOLUTION": casify(solution),
        "DIM_EQUATION_COUNT": sp.Integer(equation_count),
        "DIM_UNKNOWN_COUNT": sp.Integer(unknown_count),
        "DIM_COUNT_DIFFERENCE": sp.Integer(equation_count - 3 * unknown_count),
        "DIM_DETERMINACY": determinacy,
        "DIM_ENV": dim_env,
        "G_DIMS": g_dims,
        "V_DIMS": v_dims,
        "COEFF_DIMENSIONS_MAP": {entry[0]: entry[1] for entry in coeff_dimensions},
    }


def homogeneity_payload(expressions: tuple[object, ...], dim_env: dict[sp.Symbol, sp.ImmutableMatrix],
                        g_dims: dict[sp.Symbol, sp.ImmutableMatrix] | None = None,
                        v_dims: dict[sp.Symbol, sp.ImmutableMatrix] | None = None) -> sp.Tuple:
    rows = []
    for expr in expressions:
        expr = casify(expr)
        terms = expr.args if isinstance(expr, sp.Add) else (expr,)
        dims = [expression_dimension(term, dim_env, g_dims, v_dims) for term in terms]
        homogeneous = bool(dims and dims[0] is not None and all(dim is not None and sp.simplify(dim - dims[0]) == sp.zeros(3, 1) for dim in dims))
        rows.append(sp.Tuple(expr, sp.true if homogeneous else sp.false, sp.Tuple(*[dim if dim is not None else Str("UNDETERMINED") for dim in dims])))
    return sp.Tuple(*rows)


def q6r(package: str, n: int, coeff_order: tuple[sp.Symbol, ...], coeff_dims: dict[sp.Symbol, sp.ImmutableMatrix]) -> None:
    resolved = []
    unresolved = []
    for coeff in coeff_order:
        lookup_name = COEFFICIENT_NAME_MAP.get(coeff.name, coeff.name)
        coeff_row = INCOMING_LEDGER.get(lookup_name)
        dim_row = None
        if coeff_row is not None:
            dim_key = coeff_row.get("dimension_key")
            if dim_key is not None:
                dim_row = INCOMING_LEDGER.get(dim_key)
        if dim_row is None or "value" not in dim_row:
            unresolved.append(sp.Tuple(coeff, Str(lookup_name)))
        else:
            resolved.append(sp.Tuple(coeff, Str(lookup_name), Str(coeff_row["dimension_key"])))
            safe_name = coeff.name.upper().replace("_", "_")
            emit_custom_name(f"Q6R_{safe_name}_DERIVED_VECTOR")
            emit_custom_name(f"Q6R_{safe_name}_IMPORTED_VECTOR")
            emit_custom_name(f"Q6R_{safe_name}_DIFFERENCE")
            emit_custom_name(f"Q6R_{safe_name}_COEFFICIENT_PROVENANCE")
            emit_custom_name(f"Q6R_{safe_name}_DIMENSION_PROVENANCE")
            emit_quantity(package, n, f"Q6R_{safe_name}_DERIVED_VECTOR", coeff_dims[coeff], local=True, export=False)
            emit_quantity(package, n, f"Q6R_{safe_name}_IMPORTED_VECTOR", dim_row["value"], local=True, export=False)
            emit_quantity(package, n, f"Q6R_{safe_name}_DIFFERENCE", sp.simplify(coeff_dims[coeff] - dim_row["value"]), local=True, export=False)
            emit_quantity(
                package, n, f"Q6R_{safe_name}_COEFFICIENT_PROVENANCE",
                sp.Tuple(Str(str(coeff_row.get("class"))), Str(str(coeff_row.get("step")))),
                local=True, export=False,
            )
            emit_quantity(
                package, n, f"Q6R_{safe_name}_DIMENSION_PROVENANCE",
                sp.Tuple(
                    Str(str(dim_row.get("class"))),
                    Str(str(dim_row.get("step"))),
                    Str(str(dim_row.get("corroborated_steps"))),
                ),
                local=True, export=False,
            )
    emit_quantity(package, n, "Q6R_RESOLVED_COEFFICIENTS", sp.Tuple(*resolved), local=True, export=False)
    emit_quantity(package, n, "Q6R_UNRESOLVED_COEFFICIENTS", sp.Tuple(*unresolved), local=True, export=False)
    emit_quantity(package, n, "Q6R_RESIDUAL_SCOPE", CROSS_STEP_CONSISTENCY_ONLY, local=True, export=False)


def q7_objects(package: str, n: int, build: PackageBuild, q9_data: dict[str, object]) -> None:
    if n != 3:
        return
    g_symbols = sp.ImmutableMatrix(3, 3, lambda i, j: G7_ALL[i * 5 + j])
    densities = stiffness_densities(g_symbols)
    g_placeholder, _ = derivative_placeholders(n)
    sub = {g_placeholder[i, j]: g_symbols[i, j] for i in range(3) for j in range(3)}
    w_full = sp.expand(sum(term.expression.xreplace(sub) for term in build.stiffness_terms))
    curl_terms = [term for term in build.stiffness_terms if term.density_name == "S_curl"]
    curl_term = sp.expand(curl_terms[0].expression.xreplace(sub)) if curl_terms else sp.Integer(0)
    c_vec = []
    for i in range(3):
        component = 0
        for j in range(3):
            for k_idx in range(3):
                component += sp.LeviCivita(i, j, k_idx) * g_symbols[j, k_idx]
        c_vec.append(sp.expand(component))
    reference = sp.expand(sum(comp ** 2 for comp in c_vec))
    emit_quantity(package, n, "Q7_W_FULL", w_full)
    emit_quantity(package, n, "Q7_CURL_TERM", curl_term)
    emit_quantity(package, n, "Q7_CURL_DENSITY", densities["S_curl"])
    emit_quantity(package, n, "Q7_CURL_REFERENCE", reference)
    emit_quantity(package, n, "Q7_RESIDUAL", sp.expand(densities["S_curl"] - reference))


def q11_objects(package: str, n: int, roots: tuple[sp.Expr, ...], coeff_order: tuple[sp.Symbol, ...],
                assumptions: sp.Tuple, dim_env: dict[sp.Symbol, sp.ImmutableMatrix]) -> tuple[object, ...]:
    _, k, a = symbols_for_dimension(n)
    k_sq = sum(ki ** 2 for ki in k)
    bulk_assumptions = assumptions_for(package, n, include_bulk=True)
    kw_squared_values = []
    for idx, root in enumerate(roots, start=1):
        equation = sp.Eq(root, c_s0 ** 2 * (k_sq + kwSquared), evaluate=False)
        solved = sp.solve(equation, kwSquared, dict=True)
        kw_value = solved[0][kwSquared] if solved else Str("KW_SOLVE_UNAVAILABLE")
        kw_squared_values.append(kw_value)
        emit_quantity(package, n, "KW_EQUATION", equation, root=idx)
        emit_quantity(package, n, "KW_SQUARED", kw_value, root=idx)
        emit_quantity(package, n, "KW_SIGN", sign_payload(kw_value, bulk_assumptions), root=idx)
        emit_locus(
            package, n, "KW_ZERO_LOCUS", (sp.sympify(kw_value),),
            tuple(coeff_order) + (c_s0,), bulk_assumptions, root=idx,
        )
    c1_equations: tuple[object, ...] = tuple()
    c2_unknowns = sp.Tuple(*(a + (bulk_amplitude,)))
    c3_matrix = sp.zeros(0, len(c2_unknowns))
    emit_quantity(package, n, "C1_EQUATIONS", sp.Tuple(*c1_equations))
    emit_quantity(package, n, "C2_UNKNOWNS", c2_unknowns)
    emit_quantity(package, n, "C2_COUNT", sp.Integer(len(c2_unknowns)))
    emit_quantity(package, n, "C3_RANK", sp.Integer(c3_matrix.rank()))
    emit_quantity(package, n, "C4_DIFFERENCE", sp.Integer(len(c2_unknowns) - c3_matrix.rank()))
    return tuple(kw_squared_values)


def emit_q4(package: str, n: int, route: RouteSelection, roots: tuple[sp.Expr, ...],
            assumptions: sp.Tuple) -> dict[str, object]:
    _, k, _ = symbols_for_dimension(n)
    k_vector = sp.ImmutableMatrix(k)
    k_sq = sum(ki ** 2 for ki in k)
    root_matrices = []
    rank_data = []
    for idx, root in enumerate(roots, start=1):
        m_r = sp.ImmutableMatrix(route.matrix.subs(omegaSquared, root))
        rank = matrix_rank(m_r)
        nullity = n - rank
        stacked = sp.ImmutableMatrix.vstack(m_r, sp.ImmutableMatrix([list(k)]))
        stacked_rank = matrix_rank(stacked)
        transverse_nullity = n - stacked_rank
        m_dot_k = reduced_matrix(m_r * k_vector)
        basis = generic_nullspace_vectors(m_r, rank, root=idx)
        dot_k = sp.Tuple(*[reduced_expr((vec.T * k_vector)[0]) for vec in basis])
        residuals = sp.Tuple(*[
            reduced_matrix(k_sq * sp.ImmutableMatrix(vec) - ((vec.T * k_vector)[0]) * k_vector)
            for vec in basis
        ])
        basis_payload = sp.Tuple(*[sp.ImmutableMatrix(vec) for vec in basis])
        emit_quantity(package, n, "N1", m_r, root=idx)
        emit_quantity(package, n, "N2_RANK", sp.Integer(rank), root=idx)
        emit_quantity(package, n, "N2_NULLITY", sp.Integer(nullity), root=idx)
        emit_quantity(package, n, "N3_STACKED_RANK", sp.Integer(stacked_rank), root=idx)
        emit_quantity(package, n, "N3_TRANSVERSE_NULLITY", sp.Integer(transverse_nullity), root=idx)
        emit_quantity(package, n, "N4_NULLITY_DIFFERENCE", sp.Integer(nullity - transverse_nullity), root=idx)
        emit_quantity(package, n, "N5_M_DOT_K", m_dot_k, root=idx)
        emit_quantity(package, n, "N6_BASIS", basis_payload, root=idx)
        emit_quantity(package, n, "N6_DOT_K", dot_k, root=idx)
        emit_quantity(package, n, "N6_RESIDUAL", residuals, root=idx)
        emit_quantity(package, n, "N7_BASIS_COUNT", sp.Integer(len(basis)), root=idx)
        emit_quantity(package, n, "N7_RESIDUAL", sp.Integer(len(basis) - nullity), root=idx)
        root_matrices.append(m_r)
        rank_data.append((rank, stacked_rank, stacked))
    return {"root_matrices": tuple(root_matrices), "rank_data": tuple(rank_data)}


def emit_q5(package: str, n: int, roots: tuple[sp.Expr, ...]) -> None:
    _, k, _ = symbols_for_dimension(n)
    scale_sub = {ki: lambdaScale * ki for ki in k}
    for idx, root in enumerate(roots, start=1):
        scaled = sp.simplify(root.subs(scale_sub))
        emit_quantity(package, n, "SCALED", scaled, root=idx)
        emit_quantity(package, n, "UNSCALED", root, root=idx)
        if sp.simplify(root) == 0:
            ratio: object = UNDEFINED_RATIO
        else:
            ratio = sp.simplify(scaled / root)
        emit_quantity(package, n, "SCALE_RATIO", ratio, root=idx)
        emit_quantity(package, n, "SCALE_EXPONENT", scale_exponent_payload(ratio), root=idx)


LOCUS_PROTOCOL_SUFFIXES = (
    "EQUATIONS",
    "SOLUTION",
    "IDENTICALLY_SATISFIED",
    "INCONSISTENT",
    "REAL_ADMISSIBLE",
    "CANONICAL_LOCUS",
    "REAL_STATUS",
    "REAL_WITNESS",
    "REAL_STATUS_OPERANDS",
)


def record_field(record: object, field_name: str) -> object | None:
    if not isinstance(record, sp.Tuple):
        return None
    for field in record:
        if isinstance(field, sp.Tuple) and len(field) == 2 and field[0] == Str(field_name):
            return field[1]
    return None


def stratum_point_from_branch(variables: tuple[sp.Symbol, ...],
                               branch: dict[sp.Symbol, object]) -> dict[sp.Symbol, sp.Expr]:
    free_variables = tuple(variable for variable in variables if variable not in branch)
    free_wavevectors = tuple(variable for variable in free_variables if variable in K_ALL)
    values: dict[sp.Symbol, sp.Expr] = {}
    for variable in free_variables:
        if variable in K_ALL:
            values[variable] = sp.Integer(1) if variable == free_wavevectors[0] else sp.Integer(0)
        elif variable in (s, s_rho):
            values[variable] = sp.Integer(2)
        elif variable == beta:
            values[variable] = sp.Integer(0)
        else:
            values[variable] = sp.Integer(1)
    for _ in range(len(variables) + 1):
        updated = False
        for variable in variables:
            if variable not in branch:
                continue
            value = sp.simplify(sp.sympify(branch[variable]).subs(values))
            if variable not in values or values[variable] != value:
                values[variable] = value
                updated = True
        if not updated:
            break
    return {variable: values[variable] for variable in variables if variable in values}


def constancy_point_from_branch(variables: tuple[sp.Symbol, ...],
                                 branch: dict[sp.Symbol, object]) -> dict[sp.Symbol, sp.Expr]:
    free_variables = tuple(variable for variable in variables if variable not in branch)
    free_wavevectors = tuple(variable for variable in free_variables if variable in K_ALL)
    values: dict[sp.Symbol, sp.Expr] = {}
    for variable in free_variables:
        if variable in K_ALL:
            values[variable] = sp.Integer(1) if variable == free_wavevectors[0] else sp.Integer(0)
        elif variable in (s, s_rho):
            values[variable] = sp.Integer(2)
        elif variable == beta:
            values[variable] = sp.Integer(0)
        else:
            values[variable] = sp.Integer(1)
    for _ in range(len(variables) + 1):
        updated = False
        for variable in variables:
            if variable not in branch:
                continue
            value = sp.simplify(sp.sympify(branch[variable]).subs(values))
            if variable not in values or values[variable] != value:
                values[variable] = value
                updated = True
        if not updated:
            break
    return {variable: values[variable] for variable in variables if variable in values}


def exact_component_point(variables: tuple[sp.Symbol, ...], branch: dict[sp.Symbol, object],
                          equations: tuple[object, ...], assumptions: sp.Tuple) -> dict[sp.Symbol, sp.Expr] | None:
    point = stratum_point_from_branch(variables, branch)
    if tuple(point) != variables:
        return None
    if not witness_valid(equations, assumptions, point):
        return None
    return point


def component_free_parameters(n: int, coeff_order: tuple[sp.Symbol, ...],
                              branch: dict[sp.Symbol, object]) -> tuple[sp.Symbol, ...]:
    _, k, _ = symbols_for_dimension(n)
    ambient = tuple(dict.fromkeys(tuple(k) + tuple(coeff_order)))
    retained = [variable for variable in ambient if variable not in branch]
    extra_symbols = set()
    for value in branch.values():
        extra_symbols.update(sp.sympify(value).free_symbols)
    retained.extend(
        symbol for symbol in sorted(extra_symbols, key=lambda item: item.name)
        if symbol not in branch and symbol not in retained
    )
    return tuple(retained)


def admitted_strata(package: str, n: int, coeff_order: tuple[sp.Symbol, ...],
                     assumptions: sp.Tuple, candidates: list[StratumCandidate]) -> list[dict[str, object]]:
    components: list[dict[str, object]] = []
    component_indices: dict[str, int] = {}
    for candidate in candidates:
        variables = candidate.locus["variables"]
        solution = candidate.locus["solution"]
        branches = exposed_branches(solution, variables)
        entries = candidate.locus["real_admissible"]
        if branches is None:
            continue
        if len(branches) != len(entries):
            record_issue(
                f"{candidate.locus['base_quantity']}: solution branches and real-admissibility entries do not align",
                package=package, n=n, root=candidate.root,
            )
            continue
        source_tag = tag_name(
            package, n, candidate.locus["base_quantity"], root=candidate.root,
        )
        for branch, entry in zip(branches, entries):
            if record_field(entry, "STATUS_TOKEN") != Str("ADMISSIBLE"):
                continue
            component_payload = record_field(entry, "BRANCH")
            if component_payload is None:
                record_issue(
                    f"{candidate.locus['base_quantity']}: admitted branch has no BRANCH field",
                    package=package, n=n, root=candidate.root,
                )
                continue
            component_key = sp.srepr(component_payload)
            if component_key not in component_indices:
                defining_equations = tuple(
                    sp.Eq(variable, sp.sympify(branch[variable]), evaluate=False)
                    for variable in variables if variable in branch
                )
                component_indices[component_key] = len(components)
                components.append({
                    "component": component_payload,
                    "branch": branch,
                    "defining_equations": defining_equations,
                    "source_rows": [],
                })
            component = components[component_indices[component_key]]
            source_rows = component["source_rows"]
            if source_tag not in [row[2] for row in source_rows]:
                source_rows.append((candidate.source, candidate.locus, source_tag))

    promoted: list[dict[str, object]] = []
    _, k, _ = symbols_for_dimension(n)
    ambient = tuple(dict.fromkeys(tuple(k) + tuple(coeff_order)))
    for component in components:
        source_rows = component["source_rows"]
        point_variables = tuple(dict.fromkeys(
            variable
            for _source, locus, _tag in source_rows
            for variable in locus["variables"]
        ))
        source_equations = tuple(
            equation
            for _source, locus, _tag in source_rows
            for equation in locus["equations"]
        )
        point = exact_component_point(
            point_variables,
            component["branch"],
            tuple(component["defining_equations"]) + source_equations,
            assumptions,
        )
        if point is None:
            record_issue(
                "admitted candidate component had no exact point satisfying its defining equations and premises",
                package=package, n=n,
            )
            continue
        branch = component["branch"]
        component_substitution = {
            variable: sp.sympify(branch[variable])
            for variable in ambient if variable in branch
        }
        component["point_variables"] = point_variables
        component["point"] = point
        component["component_substitution"] = component_substitution
        component["free_parameters"] = component_free_parameters(n, coeff_order, branch)
        promoted.append(component)
    return promoted


def normalized_change_residuals(residuals: tuple[sp.Expr, ...]) -> tuple[sp.Expr, ...]:
    normalized = []
    for residual in residuals:
        expression = sp.factor(sp.together(sp.sympify(residual)))
        numerator = sp.factor(sp.fraction(expression)[0])
        if algebraic_zero_test(numerator) is not True:
            normalized.append(numerator)
    return tuple(normalized)


def substitution_is_regular(expressions: tuple[object, ...],
                            point: dict[sp.Symbol, sp.Expr]) -> bool:
    invalid_values = (sp.nan, sp.zoo, sp.oo, -sp.oo)
    if any(sp.sympify(value).has(*invalid_values) for value in point.values()):
        return False
    for expression in expressions:
        if not isinstance(expression, sp.Basic):
            continue
        for subexpression in sp.preorder_traversal(expression):
            if not isinstance(subexpression, sp.Pow) or subexpression.exp.is_negative is not True:
                continue
            denominator_value = sp.simplify(subexpression.base.subs(point))
            if denominator_value == 0 or denominator_value.has(*invalid_values):
                return False
    return True


def union_change_loci(left: tuple[sp.Expr, ...], right: tuple[sp.Expr, ...]) -> tuple[sp.Expr, ...]:
    left = normalized_change_residuals(left)
    right = normalized_change_residuals(right)
    if not left:
        return right
    if not right:
        return left
    return normalized_change_residuals(tuple(sp.factor(a * b) for a in left for b in right))


def component_count_record(status: str, value: object) -> sp.Tuple:
    return sp.Tuple(
        sp.Tuple(Str("STATUS_TOKEN"), Str(status)),
        sp.Tuple(Str("VALUE"), casify(value)),
    )


def emit_not_applicable_locus(package: str, n: int, base_quantity: str, *,
                              root: int | None, stratum: int) -> None:
    for suffix in LOCUS_PROTOCOL_SUFFIXES:
        emit_quantity(
            package, n, f"{base_quantity}_{suffix}", NOT_APPLICABLE,
            root=root, stratum=stratum,
        )


class ComponentConstancyAnalyzer:
    def __init__(self, free_parameters: tuple[sp.Symbol, ...], assumptions: sp.Tuple) -> None:
        self.free_parameters = free_parameters
        self.assumptions = assumptions
        self.witness_cache: dict[str, tuple[dict[sp.Symbol, sp.Expr], ...]] = {}

    def constant_certificate(self, residuals: tuple[sp.Expr, ...]) -> object | None:
        residuals = normalized_change_residuals(residuals)
        if not residuals:
            return sp.Tuple(sp.Tuple(*self.free_parameters), sp.EmptySet)
        free_set = set(self.free_parameters)
        for residual in residuals:
            if residual.free_symbols.isdisjoint(free_set) and algebraic_zero_test(residual) is False:
                derivatives = sp.Tuple(*[sp.diff(residual, parameter) for parameter in self.free_parameters])
                return sp.Tuple(sp.Tuple(*self.free_parameters), residual, derivatives)
            nonzero = sp.ask(sp.Q.nonzero(residual), sp.And(*self.assumptions))
            if nonzero is True:
                return sp.Tuple(
                    sp.Tuple(*self.free_parameters),
                    sp.Ne(residual, 0, evaluate=False),
                    self.assumptions,
                )
        return None

    def witnesses(self, residuals: tuple[sp.Expr, ...]) -> tuple[dict[sp.Symbol, sp.Expr], ...]:
        residuals = normalized_change_residuals(residuals)
        cache_key = sp.srepr(sp.Tuple(sp.Tuple(*self.free_parameters), sp.Tuple(*residuals)))
        if cache_key in self.witness_cache:
            return self.witness_cache[cache_key]
        if not residuals or not self.free_parameters:
            self.witness_cache[cache_key] = tuple()
            return tuple()
        equations = tuple(sp.Eq(residual, 0, evaluate=False) for residual in residuals)
        solve_branches = sp.solve(
            list(residuals), list(self.free_parameters),
            dict=True, simplify=False, manual=True,
        )
        branches = solve_branches if isinstance(solve_branches, list) else []
        points: list[dict[sp.Symbol, sp.Expr]] = []
        for branch in branches:
            if not isinstance(branch, dict):
                continue
            point = constancy_point_from_branch(self.free_parameters, branch)
            if tuple(point) != self.free_parameters:
                continue
            if not substitution_is_regular(equations + tuple(self.assumptions), point):
                continue
            if witness_valid(equations, self.assumptions, point):
                point_key = sp.srepr(casify(point))
                if point_key not in [sp.srepr(casify(existing)) for existing in points]:
                    points.append(point)
        self.witness_cache[cache_key] = tuple(points)
        return tuple(points)

    def analyze(self, residuals: tuple[sp.Expr, ...], generic_value: object,
                recompute: object) -> tuple[str, object, tuple[sp.Expr, ...]]:
        residuals = normalized_change_residuals(residuals)
        certificate = self.constant_certificate(residuals)
        if certificate is not None:
            return "CONSTANT", certificate, residuals
        for witness in self.witnesses(residuals):
            changed_value = recompute(witness)
            if sp.srepr(casify(changed_value)) != sp.srepr(casify(generic_value)):
                return "VARIES", NOT_APPLICABLE, residuals
        return "UNDECIDED", NOT_APPLICABLE, residuals


def emit_component_count(package: str, n: int, quantity: str, generic_value: object,
                         analyzer: ComponentConstancyAnalyzer,
                         change_residuals: tuple[sp.Expr, ...], recompute: object, *,
                         root: int | None, stratum: int) -> str:
    status, certificate, change_residuals = analyzer.analyze(
        change_residuals, generic_value, recompute,
    )
    payload_value = generic_value
    emit_quantity(
        package, n, quantity, component_count_record(status, payload_value),
        root=root, stratum=stratum,
    )
    emit_quantity(package, n, f"{quantity}_STATUS", Str(status), root=root, stratum=stratum)
    emit_quantity(
        package, n, f"{quantity}_CONSTANCY_CERTIFICATE", certificate,
        root=root, stratum=stratum,
    )
    change_base = f"{quantity}_CHANGE_LOCUS"
    if status == "VARIES":
        emit_locus(
            package, n, change_base, change_residuals,
            analyzer.free_parameters, analyzer.assumptions,
            root=root, stratum=stratum,
        )
    else:
        emit_not_applicable_locus(
            package, n, change_base, root=root, stratum=stratum,
        )
    return status


def scoped_quantity(quantity: str, point_evidence: bool) -> str:
    return f"POINT_EVIDENCE_{quantity}" if point_evidence else quantity


def scoped_polynomial_data(det_m: sp.Expr) -> tuple[sp.Poly, tuple[RootData, ...], tuple[sp.Expr, ...]]:
    poly = sp.Poly(sp.factor(det_m), omegaSquared)
    if not isinstance(poly.degree(), int):
        raise ValueError("component determinant does not have a finite omegaSquared degree")
    pairs = root_multiplicity_pairs(det_m)
    roots = distinct_roots(pairs)
    return poly, pairs, roots


def scoped_q3_count_values(det_m: sp.Expr) -> tuple[object, object, object, object]:
    if algebraic_zero_test(det_m) is True:
        return (
            NOT_DEFINED_ON_COMPONENT,
            NOT_DEFINED_ON_COMPONENT,
            NOT_DEFINED_ON_COMPONENT,
            NOT_DEFINED_ON_COMPONENT,
        )
    poly, pairs, roots = scoped_polynomial_data(det_m)
    root_count_all = sp.Integer(sum(pair.multiplicity for pair in pairs))
    degree = sp.Integer(poly.degree())
    return root_count_all, degree, degree - root_count_all, sp.Integer(len(roots))


def emit_scoped_q3(package: str, n: int, matrix: sp.ImmutableMatrix,
                   coeff_order: tuple[sp.Symbol, ...], assumptions: sp.Tuple, *,
                   stratum: int, free_parameters: tuple[sp.Symbol, ...],
                   point_evidence: bool,
                   analyzer: ComponentConstancyAnalyzer | None = None) -> tuple[sp.Expr, ...]:
    _, k, _ = symbols_for_dimension(n)
    det_m = exact_determinant(matrix)
    poly, pairs, roots = scoped_polynomial_data(det_m)
    prefix = lambda quantity: scoped_quantity(quantity, point_evidence)
    emit_quantity(package, n, prefix("DET_M"), sp.factor(det_m), stratum=stratum)
    emit_quantity(
        package, n, prefix("ROOT_MULTIPLICITY_PAIRS"),
        sp.Tuple(*[sp.Tuple(pair.value, sp.Integer(pair.multiplicity)) for pair in pairs]),
        stratum=stratum,
    )
    solution_set = sp.solveset(sp.Eq(det_m, 0, evaluate=False), omegaSquared)
    emit_quantity(package, n, prefix("ROOT_SOLUTION_SET"), solution_set, stratum=stratum)
    root_count_all = sp.Integer(sum(pair.multiplicity for pair in pairs))
    degree = sp.Integer(poly.degree())
    degree_residual = degree - root_count_all
    distinct_count = sp.Integer(len(roots))
    leading_change = normalized_change_residuals((poly.LC(),))
    distinct_factors = list(leading_change)
    distinct_factors.extend(sp.simplify(roots[p] - roots[q]) for p, q in combinations(range(len(roots)), 2))
    distinct_change = normalized_change_residuals(
        (sp.prod(distinct_factors),) if distinct_factors else tuple()
    )

    if point_evidence:
        emit_quantity(package, n, prefix("ROOT_COUNT_ALL"), root_count_all, stratum=stratum)
        emit_quantity(package, n, prefix("DET_M_DEGREE"), degree, stratum=stratum)
        emit_quantity(package, n, prefix("ROOT_DEGREE_RESIDUAL"), degree_residual, stratum=stratum)
        emit_quantity(package, n, prefix("ROOT_DISTINCT"), sp.Tuple(*roots), stratum=stratum)
        emit_quantity(package, n, prefix("ROOT_COUNT_DISTINCT"), distinct_count, stratum=stratum)
    else:
        if analyzer is None:
            raise ValueError("component Q3 count emission requires a constancy analyzer")
        emit_component_count(
            package, n, "ROOT_COUNT_ALL", root_count_all, analyzer, leading_change,
            lambda point: scoped_q3_count_values(det_m.subs(point))[0],
            root=None, stratum=stratum,
        )
        emit_component_count(
            package, n, "DET_M_DEGREE", degree, analyzer, leading_change,
            lambda point: scoped_q3_count_values(det_m.subs(point))[1],
            root=None, stratum=stratum,
        )
        emit_component_count(
            package, n, "ROOT_DEGREE_RESIDUAL", degree_residual, analyzer, leading_change,
            lambda point: scoped_q3_count_values(det_m.subs(point))[2],
            root=None, stratum=stratum,
        )
        emit_quantity(package, n, "ROOT_DISTINCT", sp.Tuple(*roots), stratum=stratum)
        emit_component_count(
            package, n, "ROOT_COUNT_DISTINCT", distinct_count, analyzer, distinct_change,
            lambda point: scoped_q3_count_values(det_m.subs(point))[3],
            root=None, stratum=stratum,
        )
    emit_quantity(package, n, prefix("ROOT_ORDERING"), sp.Tuple(*roots), stratum=stratum)
    for idx, root_value in enumerate(roots, start=1):
        emit_quantity(package, n, prefix("VALUE"), root_value, root=idx, stratum=stratum)
        emit_quantity(
            package, n, prefix("SIGN"), sign_payload(root_value, assumptions),
            root=idx, stratum=stratum,
        )

    remaining_k = tuple(variable for variable in k if variable in free_parameters)
    remaining_coeff = tuple(variable for variable in coeff_order if variable in free_parameters)
    if len(roots) < 2:
        emit_locus(
            package, n, prefix("ROOT_COINCIDENCE_K"), tuple(), remaining_k, assumptions,
            stratum=stratum,
        )
        emit_locus(
            package, n, prefix("ROOT_COINCIDENCE_COEFF"), tuple(), remaining_coeff, assumptions,
            stratum=stratum,
        )
    else:
        for p, q in combinations(range(len(roots)), 2):
            pair_name_k = prefix(f"ROOT_COINCIDENCE_R{p + 1}_R{q + 1}_K")
            pair_name_coeff = prefix(f"ROOT_COINCIDENCE_R{p + 1}_R{q + 1}_COEFF")
            emit_custom_name(pair_name_k)
            emit_custom_name(pair_name_coeff)
            residual = sp.simplify(roots[p] - roots[q])
            emit_locus(
                package, n, pair_name_k, (residual,), remaining_k, assumptions,
                stratum=stratum,
            )
            emit_locus(
                package, n, pair_name_coeff, (residual,), remaining_coeff, assumptions,
                stratum=stratum,
            )
    return roots


def emit_scoped_q4(package: str, n: int, matrix: sp.ImmutableMatrix,
                   roots: tuple[sp.Expr, ...], assumptions: sp.Tuple, *,
                   stratum: int, point_evidence: bool,
                   restriction: dict[sp.Symbol, sp.Expr],
                   analyzer: ComponentConstancyAnalyzer | None = None) -> list[str]:
    _, k, _ = symbols_for_dimension(n)
    prefix = lambda quantity: scoped_quantity(quantity, point_evidence)
    restricted_k = tuple(sp.sympify(ki).subs(restriction) for ki in k)
    k_vector = sp.ImmutableMatrix(restricted_k)
    k_sq = sum(ki ** 2 for ki in restricted_k)
    statuses: list[str] = []
    for idx, root_value in enumerate(roots, start=1):
        m_r = sp.ImmutableMatrix(matrix.subs(omegaSquared, root_value))
        rank = matrix_rank(m_r)
        nullity = n - rank
        stacked = sp.ImmutableMatrix.vstack(m_r, sp.ImmutableMatrix([list(restricted_k)]))
        stacked_rank = matrix_rank(stacked)
        transverse_nullity = n - stacked_rank
        m_dot_k = reduced_matrix(m_r * k_vector)
        basis = generic_nullspace_vectors(m_r, rank, root=idx)
        dot_k = sp.Tuple(*[reduced_expr((vec.T * k_vector)[0]) for vec in basis])
        residuals = sp.Tuple(*[
            reduced_matrix(k_sq * sp.ImmutableMatrix(vec) - ((vec.T * k_vector)[0]) * k_vector)
            for vec in basis
        ])
        basis_payload = sp.Tuple(*[sp.ImmutableMatrix(vec) for vec in basis])
        emit_quantity(package, n, prefix("N1"), m_r, root=idx, stratum=stratum)

        rank_change = normalized_change_residuals(all_minors(m_r, rank)) if rank else tuple()
        stacked_change = normalized_change_residuals(all_minors(stacked, stacked_rank)) if stacked_rank else tuple()
        combined_change = union_change_loci(rank_change, stacked_change)

        def rank_at(point: dict[sp.Symbol, sp.Expr]) -> sp.Integer:
            return sp.Integer(matrix_rank(sp.ImmutableMatrix(m_r.subs(point))))

        def stacked_rank_at(point: dict[sp.Symbol, sp.Expr]) -> sp.Integer:
            return sp.Integer(matrix_rank(sp.ImmutableMatrix(stacked.subs(point))))

        def basis_count_at(point: dict[sp.Symbol, sp.Expr]) -> sp.Integer:
            point_matrix = sp.ImmutableMatrix(m_r.subs(point))
            point_rank = matrix_rank(point_matrix)
            return sp.Integer(len(generic_nullspace_vectors(point_matrix, point_rank, root=idx)))

        def basis_residual_at(point: dict[sp.Symbol, sp.Expr]) -> sp.Integer:
            point_matrix = sp.ImmutableMatrix(m_r.subs(point))
            point_rank = matrix_rank(point_matrix)
            point_basis = generic_nullspace_vectors(point_matrix, point_rank, root=idx)
            return sp.Integer(len(point_basis) - (n - point_rank))

        count_rows = (
            ("N2_RANK", sp.Integer(rank), rank_change, rank_at),
            ("N2_NULLITY", sp.Integer(nullity), rank_change, lambda point: sp.Integer(n) - rank_at(point)),
            ("N3_STACKED_RANK", sp.Integer(stacked_rank), stacked_change, stacked_rank_at),
            ("N3_TRANSVERSE_NULLITY", sp.Integer(transverse_nullity), stacked_change, lambda point: sp.Integer(n) - stacked_rank_at(point)),
            ("N4_NULLITY_DIFFERENCE", sp.Integer(nullity - transverse_nullity), combined_change,
             lambda point: (sp.Integer(n) - rank_at(point)) - (sp.Integer(n) - stacked_rank_at(point))),
            ("N7_BASIS_COUNT", sp.Integer(len(basis)), rank_change, basis_count_at),
            ("N7_RESIDUAL", sp.Integer(len(basis) - nullity), rank_change, basis_residual_at),
        )
        if point_evidence:
            for quantity, value, _change, _recompute in count_rows[:5]:
                emit_quantity(package, n, prefix(quantity), value, root=idx, stratum=stratum)
        else:
            if analyzer is None:
                raise ValueError("component Q4 count emission requires a constancy analyzer")
            for quantity, value, change, recompute in count_rows[:5]:
                statuses.append(emit_component_count(
                    package, n, quantity, value, analyzer, change, recompute,
                    root=idx, stratum=stratum,
                ))
        emit_quantity(package, n, prefix("N5_M_DOT_K"), m_dot_k, root=idx, stratum=stratum)
        emit_quantity(package, n, prefix("N6_BASIS"), basis_payload, root=idx, stratum=stratum)
        emit_quantity(package, n, prefix("N6_DOT_K"), dot_k, root=idx, stratum=stratum)
        emit_quantity(package, n, prefix("N6_RESIDUAL"), residuals, root=idx, stratum=stratum)
        if point_evidence:
            for quantity, value, _change, _recompute in count_rows[5:]:
                emit_quantity(package, n, prefix(quantity), value, root=idx, stratum=stratum)
        else:
            for quantity, value, change, recompute in count_rows[5:]:
                statuses.append(emit_component_count(
                    package, n, quantity, value, analyzer, change, recompute,
                    root=idx, stratum=stratum,
                ))
    return statuses


def emit_stratum(package: str, n: int, stratum: int, component: dict[str, object],
                 route: RouteSelection, coeff_order: tuple[sp.Symbol, ...],
                 assumptions: sp.Tuple) -> None:
    source_rows = component["source_rows"]
    sources = tuple(
        Str(source)
        for source in ("RANK_DROP", "STACKED_DROP", "ROOT_COINCIDENCE")
        if source in [row[0] for row in source_rows]
    )
    source_tags = sp.Tuple(*[Str(row[2]) for row in source_rows])
    defining_equations = tuple(component["defining_equations"])
    free_parameters = tuple(component["free_parameters"])
    point = component["point"]
    point_payload = sp.Tuple(*[
        sp.Tuple(variable, point[variable]) for variable in component["point_variables"]
    ])
    emit_quantity(package, n, "SOURCES", sp.Tuple(*sources), stratum=stratum)
    emit_quantity(package, n, "SOURCE_LOCUS_TAGS", source_tags, stratum=stratum)
    emit_quantity(
        package, n, "DEFINING_EQUATIONS", sp.Tuple(*defining_equations), stratum=stratum,
    )
    emit_quantity(package, n, "FREE_PARAMETERS", sp.Tuple(*free_parameters), stratum=stratum)
    emit_quantity(package, n, "POINT", point_payload, stratum=stratum)

    component_substitution = component["component_substitution"]
    component_matrix = sp.ImmutableMatrix(route.matrix.subs(component_substitution))
    component_assumptions = sp.Tuple(*[
        premise.subs(component_substitution) if hasattr(premise, "subs") else premise
        for premise in assumptions
    ])
    analyzer = ComponentConstancyAnalyzer(free_parameters, component_assumptions)
    component_roots = emit_scoped_q3(
        package, n, component_matrix, coeff_order, component_assumptions,
        stratum=stratum, free_parameters=free_parameters,
        point_evidence=False, analyzer=analyzer,
    )
    count_statuses = emit_scoped_q4(
        package, n, component_matrix, component_roots, component_assumptions,
        stratum=stratum, point_evidence=False,
        restriction=component_substitution, analyzer=analyzer,
    )
    q3_count_statuses = [
        EMITTER.emitted[tag_name(package, n, f"{quantity}_STATUS", stratum=stratum)]
        for quantity in ("ROOT_COUNT_ALL", "DET_M_DEGREE", "ROOT_DEGREE_RESIDUAL", "ROOT_COUNT_DISTINCT")
    ]
    all_statuses = q3_count_statuses + [Str(status) for status in count_statuses]
    coverage = "COMPLETE_COVERAGE" if all(status != Str("UNDECIDED") for status in all_statuses) else "INCOMPLETE_COVERAGE"
    emit_quantity(package, n, "COMPONENT_Q3_Q4_COVERAGE", Str(coverage), stratum=stratum)

    point_matrix = sp.ImmutableMatrix(route.matrix.subs(point))
    point_assumptions = sp.Tuple(*[
        premise.subs(point) if hasattr(premise, "subs") else premise
        for premise in assumptions
    ])
    point_free_parameters = tuple(parameter for parameter in free_parameters if parameter not in point)
    point_roots = emit_scoped_q3(
        package, n, point_matrix, coeff_order, point_assumptions,
        stratum=stratum, free_parameters=point_free_parameters,
        point_evidence=True,
    )
    emit_scoped_q4(
        package, n, point_matrix, point_roots, point_assumptions,
        stratum=stratum, point_evidence=True, restriction=point,
    )


def emit_q8(package: str, n: int, roots: tuple[sp.Expr, ...], q4_data: dict[str, object],
            coeff_order: tuple[sp.Symbol, ...], assumptions: sp.Tuple,
            route: RouteSelection, root_coincidence_candidates: list[StratumCandidate]) -> None:
    _, k, _ = symbols_for_dimension(n)
    rank_candidates: list[StratumCandidate] = []
    stacked_candidates: list[StratumCandidate] = []
    for idx, root in enumerate(roots, start=1):
        m_r = q4_data["root_matrices"][idx - 1]
        rank, stacked_rank, stacked = q4_data["rank_data"][idx - 1]
        rank_minors = all_minors(m_r, rank)
        stacked_minors = all_minors(stacked, stacked_rank)
        emit_quantity(package, n, "RANK_DROP_MINORS", sp.Tuple(*rank_minors), root=idx)
        emit_quantity(package, n, "STACKED_DROP_MINORS", sp.Tuple(*stacked_minors), root=idx)
        for suffix, variables in (("K", k), ("COEFF", coeff_order), ("JOINT", tuple(k) + tuple(coeff_order))):
            locus = emit_locus(package, n, f"RANK_DROP_{suffix}", tuple(rank_minors), tuple(variables), assumptions, root=idx)
            rank_candidates.append(StratumCandidate("RANK_DROP", locus, idx))
        for suffix, variables in (("K", k), ("COEFF", coeff_order), ("JOINT", tuple(k) + tuple(coeff_order))):
            locus = emit_locus(package, n, f"STACKED_DROP_{suffix}", tuple(stacked_minors), tuple(variables), assumptions, root=idx)
            stacked_candidates.append(StratumCandidate("STACKED_DROP", locus, idx))
    stratum_candidates = rank_candidates + stacked_candidates + root_coincidence_candidates
    strata = admitted_strata(package, n, coeff_order, assumptions, stratum_candidates)
    emit_quantity(
        package, n, "STRATUM_ORDERING",
        sp.Tuple(*[component["component"] for component in strata]),
    )
    for stratum, component in enumerate(strata, start=1):
        emit_stratum(package, n, stratum, component, route, coeff_order, assumptions)


def emit_q10(package: str, n: int, roots: tuple[sp.Expr, ...], coeff_order: tuple[sp.Symbol, ...]) -> sp.ImmutableMatrix:
    jac = sp.ImmutableMatrix([[sp.simplify(sp.diff(root, coeff)) for coeff in coeff_order] for root in roots])
    emit_quantity(package, n, "ROOT_COEFFICIENT_JACOBIAN", jac)
    return jac


def emit_q3(package: str, n: int, det_m: sp.Expr, coeff_order: tuple[sp.Symbol, ...],
            assumptions: sp.Tuple,
            root_coincidence_candidates: list[StratumCandidate] | None = None) -> tuple[sp.Expr, ...]:
    _, k, _ = symbols_for_dimension(n)
    pairs = root_multiplicity_pairs(det_m)
    poly = sp.Poly(sp.factor(det_m), omegaSquared)
    roots = distinct_roots(pairs)
    emit_quantity(package, n, "DET_M", sp.factor(det_m))
    emit_quantity(package, n, "ROOT_MULTIPLICITY_PAIRS", sp.Tuple(*[sp.Tuple(pair.value, sp.Integer(pair.multiplicity)) for pair in pairs]))
    try:
        solution_set = sp.solveset(sp.Eq(det_m, 0, evaluate=False), omegaSquared)
    except Exception as exc:
        solution_set = Str(f"SOLVESET_UNAVAILABLE_{type(exc).__name__}")
        record_issue(f"root solution set unavailable after {type(exc).__name__}", package=package, n=n)
    emit_quantity(package, n, "ROOT_SOLUTION_SET", solution_set)
    root_count_all = sum(pair.multiplicity for pair in pairs)
    emit_quantity(package, n, "ROOT_COUNT_ALL", sp.Integer(root_count_all))
    emit_quantity(package, n, "DET_M_DEGREE", sp.Integer(poly.degree()))
    emit_quantity(package, n, "ROOT_DEGREE_RESIDUAL", sp.Integer(poly.degree() - root_count_all))
    emit_quantity(package, n, "ROOT_DISTINCT", sp.Tuple(*roots))
    emit_quantity(package, n, "ROOT_COUNT_DISTINCT", sp.Integer(len(roots)))
    emit_quantity(package, n, "ROOT_ORDERING", sp.Tuple(*roots))
    for idx, root in enumerate(roots, start=1):
        emit_quantity(package, n, "VALUE", root, root=idx)
        emit_quantity(package, n, "SIGN", sign_payload(root, assumptions), root=idx)
    if len(roots) < 2:
        locus_k = emit_locus(package, n, "ROOT_COINCIDENCE_K", tuple(), k, assumptions)
        locus_coeff = emit_locus(package, n, "ROOT_COINCIDENCE_COEFF", tuple(), coeff_order, assumptions)
        if root_coincidence_candidates is not None:
            root_coincidence_candidates.append(StratumCandidate("ROOT_COINCIDENCE", locus_k, None))
            root_coincidence_candidates.append(StratumCandidate("ROOT_COINCIDENCE", locus_coeff, None))
    else:
        for p, q in combinations(range(len(roots)), 2):
            pair_name_k = f"ROOT_COINCIDENCE_R{p + 1}_R{q + 1}_K"
            pair_name_coeff = f"ROOT_COINCIDENCE_R{p + 1}_R{q + 1}_COEFF"
            emit_custom_name(pair_name_k)
            emit_custom_name(pair_name_coeff)
            residual = sp.simplify(roots[p] - roots[q])
            locus_k = emit_locus(package, n, pair_name_k, (residual,), k, assumptions)
            locus_coeff = emit_locus(package, n, pair_name_coeff, (residual,), coeff_order, assumptions)
            if root_coincidence_candidates is not None:
                root_coincidence_candidates.append(StratumCandidate("ROOT_COINCIDENCE", locus_k, None))
                root_coincidence_candidates.append(StratumCandidate("ROOT_COINCIDENCE", locus_coeff, None))
    return roots


def run_cell(package: str, n: int, q9_cache: dict[int, dict[str, object]]) -> bool:
    q9_data = q9_cache[n]
    emit_q9(package, n, q9_data)
    emit_quantity(package, n, "PREMISE_INVENTORY", assumption_text_inventory(package, n), local=True, export=False, class_tag="PREMISE")
    g, v = derivative_placeholders(n)
    xs, k, a = symbols_for_dimension(n)
    build = package_build(package, n, q9_data)
    pd_term = to_coordinate(q9_data["PD_DENSITY_PLACEHOLDER"], n, g, v)
    emit_quantity(package, n, "PD_TERM", pd_term)

    kinetic_payload = sp.Tuple(*[to_coordinate(term.expression, n, g, v) for term in build.kinetic_terms])
    stiffness_payload = sp.Tuple(*[to_coordinate(term.expression, n, g, v) for term in build.stiffness_terms])
    lagrangian_coord = to_coordinate(build.lagrangian, n, g, v)
    el_system = sp.Tuple(*euler_lagrange_from_placeholders(build.lagrangian, n, g, v))
    emit_quantity(package, n, "LAGRANGIAN", lagrangian_coord, class_tag="PREMISE")
    emit_quantity(package, n, "KINETIC_TERMS", kinetic_payload, class_tag="PREMISE")
    emit_quantity(package, n, "STIFFNESS_TERMS", stiffness_payload, class_tag="PREMISE")
    emit_quantity(package, n, "EULER_LAGRANGE_SYSTEM", el_system)

    m_a, stripped = route_a_matrix(build.lagrangian, n, g, v, k, a)
    m_b, averaged_l = route_b_matrix(build.lagrangian, n, g, v, k, a)
    route = RouteSelection(M_B_TOKEN, m_b)
    emit_quantity(package, n, "M_ROUTE_A_STRIPPED_FACTOR", stripped)
    emit_quantity(package, n, "M_A", m_a)
    emit_quantity(package, n, "M_B", m_b)
    emit_quantity(package, n, "M_RESIDUAL", sp.simplify(m_a - m_b))
    emit_quantity(package, n, "M_RATIO", sp.simplify(m_a[0, 0] / m_b[0, 0]))
    emit_quantity(package, n, "M_ROUTE_RESIDUAL_SCOPE", CODING_CONSISTENCY_ONLY)
    emit_quantity(package, n, "M_ROUTE_USED", route.token)
    emit_quantity(package, n, "M_COEFFICIENT_JACOBIAN", sp.Tuple(*[sp.ImmutableMatrix([[sp.simplify(sp.diff(route.matrix[i, j], coeff)) for j in range(n)] for i in range(n)]) for coeff in build.coefficient_ordering]))

    assumptions = assumptions_for(package, n)
    det_m = determinant_from_live_matrix(route.matrix, k)
    root_coincidence_candidates: list[StratumCandidate] = []
    roots = emit_q3(
        package, n, det_m, build.coefficient_ordering, assumptions,
        root_coincidence_candidates,
    )
    q4_data = emit_q4(package, n, route, roots, assumptions)
    emit_q5(package, n, roots)

    dim_data = dimension_solve(package, n, build, g, v)
    for quantity in ("COEFFICIENT_ORDERING", "DIM_COEFFICIENTS", "DIM_EQUATIONS", "DIM_SOLUTION",
                     "DIM_EQUATION_COUNT", "DIM_UNKNOWN_COUNT", "DIM_COUNT_DIFFERENCE", "DIM_DETERMINACY"):
        emit_quantity(package, n, quantity, dim_data[quantity])
    k_sq = sum(ki ** 2 for ki in k)
    for idx, root in enumerate(roots, start=1):
        dim_expr = expression_dimension(sp.simplify(root / k_sq), dim_data["DIM_ENV"])
        emit_quantity(package, n, "DIM_OVER_KSQ", dim_expr if dim_expr is not None else Str("UNDETERMINED"), root=idx)

    action_terms = tuple(term.expression for term in build.kinetic_terms) + tuple(-term.expression for term in build.stiffness_terms)
    emit_quantity(package, n, "DIM_HOMOGENEITY_ACTION", homogeneity_payload(action_terms, dim_data["DIM_ENV"], dim_data["G_DIMS"], dim_data["V_DIMS"]))
    derived_for_homogeneity: list[object] = [det_m]
    derived_for_homogeneity.extend(roots)
    for m_dot in [EMITTER.emitted[tag_name(package, n, "N5_M_DOT_K", root=i)] for i in range(1, len(roots) + 1)]:
        derived_for_homogeneity.append(m_dot)
    for residual in [EMITTER.emitted[tag_name(package, n, "N6_RESIDUAL", root=i)] for i in range(1, len(roots) + 1)]:
        derived_for_homogeneity.append(residual)
    if n == 3:
        pass
    emit_quantity(package, n, "DIM_HOMOGENEITY_DERIVED", homogeneity_payload(tuple(derived_for_homogeneity), dim_data["DIM_ENV"]))
    q6r(package, n, build.coefficient_ordering, dim_data["COEFF_DIMENSIONS_MAP"])

    q7_objects(package, n, build, q9_data)
    emit_q8(
        package, n, roots, q4_data, build.coefficient_ordering, assumptions,
        route, root_coincidence_candidates,
    )
    emit_q10(package, n, roots, build.coefficient_ordering)
    kw_values = q11_objects(package, n, roots, build.coefficient_ordering, assumptions, dim_data["DIM_ENV"])
    if kw_values:
        emit_quantity(package, n, "DIM_HOMOGENEITY_DERIVED_Q11_KW", homogeneity_payload(tuple(kw_values), dim_data["DIM_ENV"]))
    return True


def compare_objects(left: object, right: object) -> tuple[str, object]:
    left = casify(left)
    right = casify(right)
    if isinstance(left, sp.Symbol) and isinstance(right, sp.Symbol):
        same = left.name == right.name and left.assumptions0 == right.assumptions0
        return ("PROVED_EQUAL" if same else "PROVED_DIFFERENT", sp.Tuple(left, right))
    if isinstance(left, Str) or isinstance(right, Str):
        same = sp.srepr(left) == sp.srepr(right)
        return ("PROVED_EQUAL" if same else "PROVED_DIFFERENT", sp.Tuple(left, right))
    if isinstance(left, sp.MatrixBase) and isinstance(right, sp.MatrixBase):
        if left.shape != right.shape:
            return "PROVED_DIFFERENT", sp.Tuple(left, right)
        residual = sp.simplify(left - right)
        if all(entry == 0 for entry in residual):
            return "PROVED_EQUAL", residual
        if any(entry != 0 and not entry.free_symbols for entry in residual):
            return "PROVED_DIFFERENT", residual
        return "UNDECIDED", residual
    if isinstance(left, sp.Tuple) and isinstance(right, sp.Tuple):
        if len(left) != len(right):
            return "PROVED_DIFFERENT", sp.Tuple(left, right)
        statuses = [compare_objects(l_item, r_item)[0] for l_item, r_item in zip(left, right)]
        if all(status == "PROVED_EQUAL" for status in statuses):
            return "PROVED_EQUAL", sp.Tuple(left, right)
        if any(status == "PROVED_DIFFERENT" for status in statuses):
            return "PROVED_DIFFERENT", sp.Tuple(left, right)
        return "UNDECIDED", sp.Tuple(left, right)
    if sp.srepr(left) == sp.srepr(right):
        return "PROVED_EQUAL", sp.Tuple(left, right)
    try:
        residual = sp.simplify(left - right)
        if residual == 0:
            return "PROVED_EQUAL", residual
        if residual != 0 and not getattr(residual, "free_symbols", set()):
            return "PROVED_DIFFERENT", residual
        return "UNDECIDED", residual
    except Exception as exc:
        record_issue(f"object comparison undecided after {type(exc).__name__}")
        return "UNDECIDED", sp.Tuple(left, right)


def record_source(record: dict[str, object]) -> str:
    lines = ["    {"]
    for key, value in record.items():
        if key in {"value", "f9_operands"}:
            lines.append(f"        {key!r}: _restore({sp.srepr(casify(value))!r}),")
        else:
            lines.append(f"        {key!r}: {value!r},")
    lines.append("    }")
    return "\n".join(lines)


def make_record(candidate: CandidateRow) -> dict[str, object]:
    record = {
        "display": display(candidate.value),
        "value": casify(candidate.value),
        "value_kind": "COMPUTED_OBJECT",
        "class": candidate.class_tag,
        "step": "S11",
    }
    if candidate.dimension_key is not None:
        record["dimension_key"] = candidate.dimension_key
    if candidate.description is not None:
        record["description"] = candidate.description
    return record


def add_symbol_candidates(candidates: list[CandidateRow], emitted_values: list[object],
                          dimension_keys: dict[sp.Symbol, str]) -> None:
    seen = {candidate.base_key for candidate in candidates}
    symbols = set()
    for value in emitted_values:
        value = casify(value)
        if isinstance(value, sp.MatrixBase):
            for entry in value:
                symbols.update(entry.free_symbols)
        elif hasattr(value, "free_symbols"):
            symbols.update(value.free_symbols)
        if isinstance(value, sp.Tuple):
            for item in value:
                if hasattr(item, "free_symbols"):
                    symbols.update(item.free_symbols)
    for symbol in sorted(symbols, key=lambda sym: sym.name):
        if symbol not in DECLARED_SYMBOLS:
            record_issue(f"unclassifiable free symbol {symbol}")
            continue
        if symbol.name in seen:
            continue
        metadata = DECLARED_SYMBOLS[symbol]
        candidates.append(
            CandidateRow(
                symbol.name,
                symbol,
                metadata["class"],
                "FREE_SYMBOL_POPULATION",
                dimension_keys.get(symbol),
                metadata["description"],
            )
        )
        seen.add(symbol.name)


def dimension_key_candidates(main_dim_data: dict[int, dict[str, object]]) -> tuple[list[CandidateRow], dict[sp.Symbol, str]]:
    rows: list[CandidateRow] = []
    symbol_dimension_keys: dict[sp.Symbol, str] = {}
    imported_dimension_rows = {
        key: record
        for key, record in INCOMING_LEDGER.items()
        if isinstance(record.get("value"), sp.MatrixBase) and str(key).endswith("dimension")
    }
    for n, dim_data in main_dim_data.items():
        for coeff, vector in dim_data["COEFF_DIMENSIONS_MAP"].items():
            match_key = None
            for key, record in imported_dimension_rows.items():
                status, _ = compare_objects(vector, record["value"])
                if status == "PROVED_EQUAL":
                    match_key = key
                    break
            if match_key is None:
                match_key = f"{coeff.name}_dimension"
            symbol_dimension_keys[coeff] = match_key
            rows.append(CandidateRow(match_key, vector, "DERIVED", f"DIMENSION_ROW_{coeff.name}_D{n}"))
    unique: dict[str, CandidateRow] = {}
    for row in rows:
        if row.base_key not in unique:
            unique[row.base_key] = row
    return list(unique.values()), symbol_dimension_keys


def merged_export(main_dim_data: dict[int, dict[str, object]],
                  *, emit_diagnostics: bool = True) -> dict[str, dict[str, object]]:
    candidates = list(EMITTER.export_candidates)
    dimension_rows, symbol_dimension_keys = dimension_key_candidates(main_dim_data)
    candidates.extend(dimension_rows)
    emitted_values = [candidate.value for candidate in candidates]
    add_symbol_candidates(candidates, emitted_values, symbol_dimension_keys)

    key_operands = []
    merged = {key: dict(record) for key, record in INCOMING_LEDGER.items()}
    for candidate in candidates:
        candidate_record = make_record(candidate)
        imported = INCOMING_LEDGER.get(candidate.base_key)
        if imported is None:
            write_key = candidate.base_key
            relation_token = Str("ABSENT")
            operands = sp.Tuple(Str(candidate.base_key), sp.Tuple())
        else:
            value_status, value_operand = compare_objects(candidate_record["value"], imported.get("value"))
            class_status = "PROVED_EQUAL" if candidate_record.get("class") == imported.get("class") else "PROVED_DIFFERENT"
            if value_status == "PROVED_EQUAL" and class_status == "PROVED_EQUAL":
                write_key = candidate.base_key
                relation_token = Str("PROVED_EQUAL")
                candidate_record["f9_operands"] = sp.Tuple(candidate_record["value"], imported.get("value"))
                imported_steps = []
                if imported.get("corroborated_steps"):
                    imported_steps.extend(imported["corroborated_steps"])
                elif imported.get("step"):
                    imported_steps.append(imported["step"])
                imported_steps.append("S11")
                candidate_record["corroborated_steps"] = tuple(dict.fromkeys(imported_steps))
            else:
                write_key = f"s11_{candidate.base_key}"
                relation_token = Str(value_status if value_status != "PROVED_EQUAL" else "CLASS_DIFFERENCE")
            operands = sp.Tuple(Str(candidate.base_key), Str(write_key), relation_token, casify(value_operand if imported is not None else sp.Tuple()))
        key_operands.append(sp.Tuple(Str(candidate.source_tag), Str(candidate.base_key), operands))
        merged[write_key] = candidate_record

    if emit_diagnostics:
        EMITTER.emit("PY_S11_LOCAL_EXPORT_CANDIDATE_KEY_OPERANDS", sp.Tuple(*key_operands), export_key=None)
    for n, dim_data in main_dim_data.items():
        export_coeff_rows = []
        for coeff in dim_data["COEFF_DIMENSIONS_MAP"]:
            lookup_name = COEFFICIENT_NAME_MAP.get(coeff.name, coeff.name)
            coeff_row = merged.get(lookup_name)
            status = Str("UNRESOLVED")
            dim_key_payload = NOT_APPLICABLE
            if coeff_row is not None and coeff_row.get("dimension_key") in merged:
                dim_record = merged[coeff_row["dimension_key"]]
                if "value" in dim_record:
                    status = Str("RESOLVED")
                    dim_key_payload = Str(coeff_row["dimension_key"])
            export_coeff_rows.append(sp.Tuple(coeff, Str(lookup_name), status, dim_key_payload))
            if emit_diagnostics:
                EMITTER.emit(
                    tag_name("MAIN", n, f"Q6R_EXPORT_LOOKUP_{coeff.name.upper()}", local=True),
                    sp.Tuple(coeff, Str(lookup_name), status, dim_key_payload),
                )
        if emit_diagnostics:
            EMITTER.emit(tag_name("MAIN", n, "Q6R_EXPORT_LOOKUP_ROWS", local=True), sp.Tuple(*export_coeff_rows))
    return merged


def write_exports(ledger: dict[str, dict[str, object]]) -> None:
    path = SCRIPT_DIR / "S11_exports.py"
    lines = [
        "# S11_exports.py -- GENERATED by S11_stray_longitudinal_sympy_audit.py. Do not edit.",
        "from types import MappingProxyType",
        "",
        "import sympy as sp",
        "from sympy.core.symbol import Str",
        "",
        "_RELATIONALS = {",
        "    'Equality': lambda left, right: sp.Eq(left, right, evaluate=False),",
        "    'Unequality': lambda left, right: sp.Ne(left, right, evaluate=False),",
        "    'StrictGreaterThan': lambda left, right: sp.Gt(left, right, evaluate=False),",
        "    'StrictLessThan': lambda left, right: sp.Lt(left, right, evaluate=False),",
        "    'GreaterThan': lambda left, right: sp.Ge(left, right, evaluate=False),",
        "    'LessThan': lambda left, right: sp.Le(left, right, evaluate=False),",
        "}",
        "",
        "def _restore(source):",
        "    return eval(source, {'__builtins__': {}, 'Str': Str, **vars(sp), **_RELATIONALS})",
        "",
        "BUILD_INPUT_DIGESTS = MappingProxyType({",
        f"    'S11_stray_longitudinal_sympy_audit.py': {file_digest(Path(__file__).resolve())!r},",
        f"    'S10_exports.py': {file_digest(SCRIPT_DIR / 'S10_exports.py')!r},",
        f"    'S11_SHARED_PHYSICS.md': {file_digest(SCRIPT_DIR.parent / 'directives' / 'S11_SHARED_PHYSICS.md')!r},",
        "})",
        "",
        "_LEDGER = {",
    ]
    for key in sorted(ledger):
        lines.append(f"    {key!r}: {record_source(ledger[key])},")
    lines.extend([
        "}",
        "LEDGER = MappingProxyType({",
        "    name: MappingProxyType(record) for name, record in _LEDGER.items()",
        "})",
        "del _LEDGER",
        "",
    ])
    path.write_text("\n".join(lines), encoding="utf-8")


def run_record_payloads(completed_pairs: list[tuple[str, int]]) -> tuple[sp.Tuple, sp.Tuple]:
    declared_pairs = [(package, n) for package in PACKAGE_ORDER for n in PACKAGE_DIMS[package]]
    skipped_pairs_list = [pair for pair in declared_pairs if pair not in completed_pairs]
    run_pairs_payload = sp.Tuple(*[sp.Tuple(Str(package), sp.Integer(n)) for package, n in completed_pairs])
    skipped_pairs_payload = sp.Tuple(*[sp.Tuple(Str(package), sp.Integer(n)) for package, n in skipped_pairs_list])
    return run_pairs_payload, skipped_pairs_payload


def main_is_complete(completed_pairs: list[tuple[str, int]]) -> bool:
    main_completed = {(package, n) for package, n in completed_pairs if package == PRIMARY_PACKAGE}
    main_declared = {(PRIMARY_PACKAGE, n) for n in PACKAGE_DIMS[PRIMARY_PACKAGE]}
    return main_completed == main_declared


def publish_time_cell_states(completed_pairs: list[tuple[str, int]],
                             attempted_pairs: list[tuple[str, int]]) -> sp.Tuple:
    completed = set(completed_pairs)
    attempted = set(attempted_pairs)
    rows = []
    for package in PACKAGE_ORDER:
        for n in PACKAGE_DIMS[package]:
            pair = (package, n)
            if pair in completed:
                state = "COMPLETED_AT_PUBLISH_TIME"
            elif pair in attempted:
                state = "ATTEMPT_FAILED_BEFORE_PUBLISH"
            else:
                state = "NOT_YET_ATTEMPTED_AT_PUBLISH_TIME"
            rows.append(sp.Tuple(Str(package), sp.Integer(n), Str(state)))
    return sp.Tuple(*rows)


def delete_stale_export() -> None:
    stale = SCRIPT_DIR / "S11_exports.py"
    if stale.exists():
        stale.unlink()


def run() -> None:
    global CURRENT_CELL, CURRENT_ISSUE_CONTEXT
    start = time.monotonic()
    q9_cache = {n: compute_q9(n) for n in (2, 3, 4, 5)}
    completed_pairs: list[tuple[str, int]] = []
    attempted_pairs: list[tuple[str, int]] = []
    main_dim_data: dict[int, dict[str, object]] = {}
    publish_error: Exception | None = None
    publish_attempted = False
    for package in PACKAGE_ORDER:
        for n in PACKAGE_DIMS[package]:
            attempted_pairs.append((package, n))
            CURRENT_CELL = (package, n)
            try:
                before_count = EMITTER.count
                run_cell(package, n, q9_cache)
                completed_pairs.append((package, n))
                if package == PRIMARY_PACKAGE:
                    dim_tag = tag_name(package, n, "DIM_COEFFICIENTS")
                    coeff_tag = tag_name(package, n, "COEFFICIENT_ORDERING")
                    coeff_dims = {}
                    for row in EMITTER.emitted[dim_tag]:
                        coeff_dims[row[0]] = row[1]
                    main_dim_data[n] = {
                        "COEFF_DIMENSIONS_MAP": coeff_dims,
                        "COEFFICIENT_ORDERING": tuple(EMITTER.emitted[coeff_tag]),
                    }
                if EMITTER.count == before_count:
                    record_issue("cell completed without emitted tags", package=package, n=n)
            except Exception as exc:
                EMITTER.emit(
                    tag_name(package, n, "CELL_EXCEPTION", local=True),
                    sp.Tuple(Str(package), sp.Integer(n), Str(type(exc).__name__), Str(repr(exc))),
                )
                record_issue(f"cell attempt failed after {type(exc).__name__}: {exc}", package=package, n=n)
            finally:
                CURRENT_CELL = None
            if not publish_attempted and main_is_complete(completed_pairs):
                publish_attempted = True
                publish_state = publish_time_cell_states(completed_pairs, attempted_pairs)
                EMITTER.emit(
                    "PY_S11_LOCAL_PUBLISH_TIME_CELL_STATES",
                    publish_state,
                    export_key="publish_time_cell_states",
                )
                CURRENT_ISSUE_CONTEXT = "PUBLISH"
                try:
                    ledger = merged_export(main_dim_data)
                    write_exports(ledger)
                except Exception as exc:
                    publish_error = exc
                    EMITTER.emit(
                        "PY_S11_LOCAL_PUBLISH_EXCEPTION",
                        sp.Tuple(Str(type(exc).__name__), Str(repr(exc))),
                    )
                    record_issue(f"S11_exports.py publish failed after {type(exc).__name__}: {exc}", context="PUBLISH")
                    try:
                        delete_stale_export()
                    except Exception as stale_exc:
                        record_issue(f"stale S11_exports.py deletion failed after {type(stale_exc).__name__}: {stale_exc}", context="PUBLISH")
                        if publish_error is None:
                            publish_error = stale_exc
                finally:
                    CURRENT_ISSUE_CONTEXT = None
    run_pairs_payload, skipped_pairs_payload = run_record_payloads(completed_pairs)
    EMITTER.emit("PY_S11_RUN_PAIRS", run_pairs_payload)
    EMITTER.emit("PY_S11_SKIPPED_PAIRS", skipped_pairs_payload)

    if not publish_attempted:
        try:
            delete_stale_export()
        except Exception as exc:
            publish_error = exc
            EMITTER.emit(
                "PY_S11_LOCAL_PUBLISH_EXCEPTION",
                sp.Tuple(Str(type(exc).__name__), Str(repr(exc))),
            )
            record_issue(f"stale S11_exports.py deletion failed after {type(exc).__name__}: {exc}", context="PUBLISH")
        record_issue("S11_exports.py not published because a declared MAIN cell did not complete", context="PUBLISH")

    EMITTER.emit("PY_S11_LOCAL_TAGS", sp.Tuple(*[Str(tag) for tag in EMITTER.local_tags]))

    runtime = sp.Float(time.monotonic() - start, 6)
    report = sp.Tuple(
        sp.Tuple(Str("SCRIPT"), Str(str(Path(__file__).resolve())), sp.Integer(len(Path(__file__).read_text(encoding="utf-8").splitlines()))),
        sp.Tuple(Str("EMITTED_TAGS"), sp.Integer(EMITTER.count)),
        sp.Tuple(Str("RUN_PAIRS"), run_pairs_payload),
        sp.Tuple(Str("SKIPPED_PAIRS"), skipped_pairs_payload),
        sp.Tuple(Str("RUNTIME_SECONDS"), runtime),
        sp.Tuple(Str("ISSUES"), sp.Tuple(*[Str(issue) for issue in ISSUES])),
        sp.Tuple(Str("CUSTOM_QUANTITY_NAMES"), sp.Tuple(*[Str(name) for name in sorted(CUSTOM_QUANTITY_NAMES)])),
    )
    EMITTER.emit("PY_S11_LOCAL_SECTION10_REPORT", report)
    if publish_error is not None:
        raise RuntimeError("S11 export publish failed; see PY_S11_LOCAL_SECTION10_REPORT") from publish_error


class StdoutTee:
    def __init__(self, *streams: object) -> None:
        self.streams = streams

    def write(self, text: str) -> int:
        for stream in self.streams:
            stream.write(text)
        return len(text)

    def flush(self) -> None:
        for stream in self.streams:
            stream.flush()


def run_with_stdout_tee() -> None:
    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    original_stdout = sys.stdout
    with OUT_PATH.open("w", encoding="utf-8") as out_stream:
        sys.stdout = StdoutTee(original_stdout, out_stream)
        try:
            run()
        finally:
            sys.stdout.flush()
            sys.stdout = original_stdout


if __name__ == "__main__":
    run_with_stdout_tee()
