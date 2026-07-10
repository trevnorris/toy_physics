#!/usr/bin/env python3
"""Ledger stage023 SymPy audit: native nullspace underdetermination.

Standalone, print-only, exact, and zero-file-I/O.  This is the pathA_34
II-G5b COMPLETING slice: it computes the genuine symbolic constraint rank,
the return-moving directions, the forward scalar/dipole residuals, the
dimensional and provenance firewalls, and the counterfactual selector
control.  The earned physics verdict is a characterized FAIL while the
audit process itself exits successfully when every tooth fires.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from functools import lru_cache
from typing import Any

import sympy as sp


PASS_COUNT = 0
FAIL_COUNT = 0

CROSS_L_RESIDUAL_PREDICTION = "CROSS_L_RESIDUAL_PREDICTION"
FAIL_UNDERDETERMINED = "FAIL_UNDERDETERMINED_NOT_PREDICTIVE"
FAIL_DECOUPLED = "FAIL_DECOUPLED"
FAIL_TAUTOLOGICAL = "FAIL_TAUTOLOGICAL"
FAIL_DIMENSIONAL = "FAIL_DIMENSIONAL"
FAIL_NO_CONSISTENT_RETURN = "FAIL_NO_CONSISTENT_RETURN"
FAIL_OVERCANCEL = "FAIL_OVERCANCEL"
FAIL_EPSILON_MISMATCH = "FAIL_EPSILON_MISMATCH"
NO_FAIL = "NO_FAIL"


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
    return sp.factor(sp.cancel(sp.together(sp.simplify(expr))))


def fmt(expr: Any) -> str:
    if isinstance(expr, bool):
        return "True" if expr else "False"
    if isinstance(expr, str):
        return expr
    try:
        return sp.sstr(compact(expr))
    except Exception:
        return sp.sstr(expr)


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
    floats = sp.sympify(expr).atoms(sp.Float)
    if floats:
        raise AuditFailure(f"{name}: Float atom(s) in exact audit expression: {floats}")


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


def report_self_consistency(name: str, condition: bool) -> None:
    """Report a transcription sanity check without counting it as an earned tooth."""
    status = "OK" if bool(condition) else "FAIL"
    print(f"SELF-CONSISTENCY {status} (not earned; not counted)  {name}")
    if not bool(condition):
        raise AuditFailure(f"self-consistency check failed: {name}")


def bool_zero(residual: sp.Expr | int) -> bool:
    assert_no_float("bool_zero", residual)
    return compact(residual) == 0


def verdict_residual(actual: str, expected: str) -> sp.Integer:
    return sp.Integer(0) if actual == expected else sp.Integer(1)


omega = sp.Symbol("omega", real=True)
a = sp.Symbol("a", positive=True, real=True)
c_s = sp.Symbol("c_s", positive=True, real=True)
M0 = sp.Symbol("M0", nonzero=True, real=True)
D1 = sp.Symbol("D1", nonzero=True, real=True)
R0, R1 = sp.symbols("R0 R1", real=True)
D0 = sp.Symbol("D0", nonzero=True, real=True)
K0c, Keta, TOmega = sp.symbols("K0c K_eta T_Omega", positive=True, real=True)
Z0ret, Z1ret = sp.symbols("Z0_ret Z1_ret", positive=True, real=True)
q_free = sp.Symbol("q_free", real=True)
gain0, gain1 = sp.symbols("gain0 gain1", positive=True, real=True)
eta_null = sp.Symbol("eta_null", positive=True, real=True)
OmegaU, OmegaW, Rmix, gU, gW = sp.symbols(
    "Omega_U Omega_W R_mix g_U g_W", nonzero=True, real=True
)

Ldim, Mdim, Tdim = sp.symbols("L M T", positive=True)
Dim = tuple[sp.Rational, sp.Rational, sp.Rational]
ZERO_DIM: Dim = (sp.Rational(0), sp.Rational(0), sp.Rational(0))


@dataclass(frozen=True)
class Mutation:
    # A transcript label only: it must not defeat caching of identical physics.
    name: str = field(default="baseline", compare=False)
    # A physics-inert, cache-bearing token used to prove a neutralized mutation
    # is independently rerun instead of aliasing the cached baseline object.
    inert_rerun_token: str = ""
    selector: str = "none"
    inject_null: bool = False
    wrong_sign_return: bool = False
    perfect_return: bool = False
    no_consistent_return: bool = False
    decouple_knobs: bool = False
    corrupt_dimension: bool = False
    corrupt_v1: bool = False
    corrupt_a1_order: bool = False
    assert_not_derive: bool = False
    emit_epsilon_magnitude_as_derived: bool = False
    corrupt_port_kernel: bool = False


GENERATOR_DOFS = [
    OmegaU,
    OmegaW,
    Rmix,
    gU,
    gW,
    D0,
    K0c,
    Keta,
    TOmega,
    Z0ret,
    Z1ret,
]


def build_port_kernel(mutation: Mutation) -> dict[str, sp.Expr]:
    omega_u = 2 * OmegaU if mutation.corrupt_port_kernel else OmegaU
    p_port = omega_u**2 * gW + Rmix * gU
    delta_port = omega_u**2 * OmegaW**2 - Rmix**2
    n0_port = compact(p_port**2 / delta_port**2)
    p0_raw = compact(n0_port / D0)
    return {
        "P_port": p_port,
        "Delta_port": delta_port,
        "N0_from_port": n0_port,
        "P0_raw": p0_raw,
        "P0_physical": compact((c_s / a) ** 2 * p0_raw),
    }


def selector_equations(mutation: Mutation) -> list[sp.Equality]:
    k1 = Keta + 2 * TOmega
    if mutation.selector in {"none", "neutered"}:
        return []
    if mutation.selector == "one_counterfactual_control":
        return [sp.Eq(Z0ret, K0c)]
    if mutation.selector in {"counterfactual_control", "asserted_unproven"}:
        return [sp.Eq(Z0ret, K0c), sp.Eq(Z1ret, k1)]
    raise ValueError(f"unknown selector mode: {mutation.selector}")


def selector_provenance(mutation: Mutation) -> dict[str, Any]:
    equations = selector_equations(mutation)
    asserted = mutation.selector == "asserted_unproven"
    control = mutation.selector in {
        "one_counterfactual_control",
        "counterfactual_control",
    }
    return {
        "present": bool(equations),
        "equations": [f"{fmt(eq.lhs)} = {fmt(eq.rhs)}" for eq in equations],
        "derived_from_named_pde": False,
        "control_only": control,
        "tautological": asserted,
        "source": (
            "counterfactual rank-collapse witness; exact Gate-6 equations are not earned"
            if control
            else "asserted selector without named PDE provenance"
            if asserted
            else "no Gate-5 selector equation"
        ),
    }


def selector_substitutions(mutation: Mutation) -> dict[sp.Symbol, sp.Expr]:
    return {eq.lhs: eq.rhs for eq in selector_equations(mutation)}


def positive_bounded_transfer(t_expr: sp.Expr, eps_expr: sp.Expr) -> bool:
    t_clean = compact(t_expr)
    one_minus = compact(1 - t_clean)
    return bool(
        eps_expr.is_positive is True
        and t_clean.is_positive is True
        and one_minus.is_positive is True
    )


def build_transfers(mutation: Mutation) -> dict[str, Any]:
    k0 = K0c
    k1 = Keta + 2 * TOmega
    substitutions = selector_substitutions(mutation)
    z0 = substitutions.get(Z0ret, Z0ret)
    z1 = substitutions.get(Z1ret, Z1ret)
    if mutation.perfect_return:
        z0 = sp.Integer(0)
        z1 = sp.Integer(0)
    elif mutation.wrong_sign_return:
        z0 = -z0
        z1 = -z1
    elif mutation.no_consistent_return:
        z0 = -2 * k0
        z1 = -2 * k1
    if mutation.inject_null:
        z0 = compact(z0 + eta_null * k0)
        z1 = compact(z1 + eta_null * k1)

    t0 = compact(k0 / (k0 + z0))
    t1 = compact(k1 / (k1 + z1))
    eps0 = compact(z0 / k0)
    eps1 = compact(z1 / k1)
    if mutation.decouple_knobs:
        t0 = compact(gain0 * t0)
        t1 = compact(gain1 * t1)

    relation_ok = (
        bool_zero(t0 - 1 / (1 + eps0))
        and bool_zero(t1 - 1 / (1 + eps1))
        and bool_zero(1 - t0 - eps0 / (1 + eps0))
        and bool_zero(1 - t1 - eps1 / (1 + eps1))
    ) if not mutation.decouple_knobs else False
    bounded = positive_bounded_transfer(t0, eps0) and positive_bounded_transfer(t1, eps1)
    overcancel = bool_zero(eps0) and bool_zero(eps1)
    no_consistent = bool_zero(eps0 + 2) and bool_zero(eps1 + 2) and not bounded
    return {
        "K0": k0,
        "K1": k1,
        "Z0": z0,
        "Z1": z1,
        "T0": t0,
        "T1": t1,
        "epsilon0": eps0,
        "epsilon1": eps1,
        "forward_relations_ok": relation_ok,
        "positive_bounded": bounded,
        "overcancel": overcancel,
        "no_consistent_return": no_consistent,
    }


def build_residuals(transfers: dict[str, Any], mutation: Mutation) -> dict[str, Any]:
    # CHECKABLE cited input from stage022; no fingerprint battery is rebuilt here.
    v0 = sp.Integer(1)
    v1 = sp.Rational(1, 3) if mutation.corrupt_v1 else sp.Rational(1, 2)
    power1 = 2 if mutation.corrupt_a1_order else 3
    t0, t1 = transfers["T0"], transfers["T1"]
    eps0, eps1 = transfers["epsilon0"], transfers["epsilon1"]
    a0 = compact(sp.I * v0 * (a * omega / c_s) * M0 * (1 - t0))
    a1 = compact(sp.I * v1 * (a * omega / c_s) ** power1 * D1 * (1 - t1))
    # Independently cited pathA_29 forward residual form.
    expected_a0 = compact(sp.I * a * omega * M0 * eps0 / (c_s * (1 + eps0)))
    expected_a1 = compact(
        sp.I * a**3 * omega**3 * D1 * eps1 / (2 * c_s**3 * (1 + eps1))
    )
    residual_a0 = compact(a0 - expected_a0)
    residual_a1 = compact(a1 - expected_a1)
    order0 = compact(a0).as_powers_dict().get(omega)
    order1 = compact(a1).as_powers_dict().get(omega)
    return {
        "v0": v0,
        "v1": v1,
        "A0": a0,
        "A1": a1,
        "expected_A0": expected_a0,
        "expected_A1": expected_a1,
        "A0_residual": residual_a0,
        "A1_residual": residual_a1,
        "A0_form": bool_zero(residual_a0),
        "A1_form": bool_zero(residual_a1),
        "A0_order": order0,
        "A1_order": order1,
    }


def rank_row(expr: sp.Expr, dofs: list[sp.Symbol]) -> list[sp.Expr]:
    return [compact(sp.diff(expr, dof)) for dof in dofs]


def matrix_rank(rows: list[list[sp.Expr]]) -> int:
    """The anti-v1 core: genuine native rank of the symbolic Jacobian rows."""
    return 0 if not rows else int(sp.Matrix(rows).rank())


def build_rank_audit(
    port: dict[str, sp.Expr], transfers: dict[str, Any], mutation: Mutation
) -> dict[str, Any]:
    dofs = list(GENERATOR_DOFS)
    if mutation.inject_null:
        dofs.append(eta_null)
    base_constraints = [port["P0_raw"], K0c, Keta + 2 * TOmega]
    selector_residuals = [eq.lhs - eq.rhs for eq in selector_equations(mutation)]
    constraints = base_constraints + selector_residuals
    rows = [rank_row(expr, dofs) for expr in constraints]
    rank0 = matrix_rank(rows)
    native_nullity = len(dofs) - rank0
    grad_t0 = rank_row(transfers["T0"], dofs)
    grad_t1 = rank_row(transfers["T1"], dofs)
    return_aug_rank = matrix_rank(rows + [grad_t0, grad_t1])
    return_moving_nullity = return_aug_rank - rank0

    witnesses: list[dict[str, Any]] = []
    for return_dof in (Z0ret, Z1ret):
        vector = [sp.Integer(1) if dof == return_dof else sp.Integer(0) for dof in dofs]
        constraint_deltas = [
            compact(sum(row_i * vec_i for row_i, vec_i in zip(row, vector)))
            for row in rows
        ]
        delta_t0 = compact(sum(x * y for x, y in zip(grad_t0, vector)))
        delta_t1 = compact(sum(x * y for x, y in zip(grad_t1, vector)))
        if all(bool_zero(value) for value in constraint_deltas):
            witnesses.append(
                {
                    "dof": fmt(return_dof),
                    "unit_vector": vector,
                    "delta_T0": delta_t0,
                    "delta_T1": delta_t1,
                    "moves_return": not (bool_zero(delta_t0) and bool_zero(delta_t1)),
                }
            )
    return {
        "dofs": dofs,
        "constraints": constraints,
        "constraint_jacobian": rows,
        "rank0": rank0,
        "native_nullspace_dimension": native_nullity,
        "return_augmented_rank": return_aug_rank,
        "return_moving_nullity": return_moving_nullity,
        "native_null_moves_return": return_moving_nullity > 0,
        "example_return_moving_directions": witnesses,
        "selector_provenance": selector_provenance(mutation),
        "pathA29_premise": {"Z_is_premise": True, "boundary_dof": "none"},
    }


def dim_add(left: Dim, right: Dim) -> Dim:
    return tuple(sp.Rational(x + y) for x, y in zip(left, right))  # type: ignore[return-value]


def dim_scale(dim: Dim, scale: sp.Rational) -> Dim:
    return tuple(sp.Rational(scale * value) for value in dim)  # type: ignore[return-value]


def dim_of(expr: sp.Expr, symbol_dims: dict[sp.Symbol, Dim]) -> Dim:
    expr = sp.sympify(expr)
    if expr == 0 or expr.is_number:
        return ZERO_DIM
    if expr.is_Symbol:
        if expr not in symbol_dims:
            raise DimError(f"missing dimension for symbol {expr}")
        return symbol_dims[expr]
    if expr.is_Mul:
        result = ZERO_DIM
        for argument in expr.args:
            result = dim_add(result, dim_of(argument, symbol_dims))
        return result
    if expr.is_Pow:
        base, exponent = expr.args
        if not exponent.is_number:
            raise DimError(f"non-numeric dimension exponent in {expr}")
        return dim_scale(dim_of(base, symbol_dims), sp.Rational(exponent))
    if expr.is_Add:
        dimensions = [dim_of(arg, symbol_dims) for arg in expr.args if arg != 0]
        if not dimensions:
            return ZERO_DIM
        if any(value != dimensions[0] for value in dimensions[1:]):
            raise DimError(f"dimension mismatch in sum {expr}: {dimensions}")
        return dimensions[0]
    raise DimError(f"unsupported expression in dimension checker: {expr}")


SOURCED_DIMS: dict[sp.Symbol, Dim] = {
    a: (1, 0, 0),
    c_s: (1, 0, -1),
    omega: (0, 0, -1),
    M0: (0, 1, -1),
    D1: (1, 1, -1),
    R0: (0, 1, -1),
    R1: (1, 1, -1),
    D0: (-1, 1, -2),
    K0c: (0, 1, -2),
    Keta: (0, 1, -2),
    TOmega: (0, 1, -2),
    Z0ret: (0, 1, -2),
    Z1ret: (0, 1, -2),
    OmegaU: (0, 0, -1),
    OmegaW: (0, 0, -1),
    Rmix: (0, 0, -2),
    gU: (sp.Rational(-1, 2), sp.Rational(1, 2), -2),
    gW: (sp.Rational(-1, 2), sp.Rational(1, 2), -2),
    eta_null: ZERO_DIM,
    gain0: ZERO_DIM,
    gain1: ZERO_DIM,
    q_free: ZERO_DIM,
}


EXPECTED_DIMS: dict[str, Dim] = {
    "A0": (0, 1, -1),
    "A1": (1, 1, -1),
    "T0": ZERO_DIM,
    "T1": ZERO_DIM,
    "epsilon0": ZERO_DIM,
    "epsilon1": ZERO_DIM,
    "P0_physical": ZERO_DIM,
}


def run_dimension_check(
    transfers: dict[str, Any],
    residuals: dict[str, Any],
    port: dict[str, sp.Expr],
    *,
    corrupt_sourced: bool = False,
    corrupt_free_carrier: bool = False,
) -> dict[str, Any]:
    dims = dict(SOURCED_DIMS)
    if corrupt_sourced:
        dims[M0] = dim_add(dims[M0], (1, 0, 0))
    if corrupt_free_carrier:
        dims[q_free] = (7, 0, 0)
    expressions = {
        "A0": residuals["A0"],
        "A1": residuals["A1"],
        "T0": transfers["T0"],
        "T1": transfers["T1"],
        "epsilon0": transfers["epsilon0"],
        "epsilon1": transfers["epsilon1"],
        "P0_physical": port["P0_physical"],
    }
    computed: dict[str, Dim] = {}
    errors: dict[str, str] = {}
    for name, expr in expressions.items():
        try:
            # An identically zero amplitude still occupies its independently
            # sourced ledger slot; this is the stage021/source dim_record rule.
            computed[name] = EXPECTED_DIMS[name] if sp.sympify(expr) == 0 else dim_of(expr, dims)
        except DimError as exc:
            errors[name] = str(exc)
    ok = not errors and all(computed[name] == expected for name, expected in EXPECTED_DIMS.items())
    return {
        "computed": computed,
        "expected": EXPECTED_DIMS,
        "errors": errors,
        "dimensional_ok": ok,
        "verdict": NO_FAIL if ok else FAIL_DIMENSIONAL,
        "checked_expressions_mention_q_free": any(q_free in expr.free_symbols for expr in expressions.values()),
    }


def build_provenance(mutation: Mutation, rank: dict[str, Any]) -> dict[str, Any]:
    items: dict[str, dict[str, Any]] = {
        "stage022_v0_v1_coefficients": {
            "tags": ["stage022_cross_l_fingerprints", "checkable_typed_input"],
            "computed_class": "cited_earned_input",
        },
        "forward_T0_T1_transfer_map": {
            "tags": ["continuity_partition", "forward_from_Z_over_K"],
            "computed_class": "derived_in_gate",
        },
        "forward_A0_A1_residual_map": {
            "tags": ["stage022_coefficients", "pathA29_independent_form"],
            "computed_class": "derived_in_gate",
        },
        "ell2_P0_map": {
            "tags": ["stage017_grouped_port_kernel"],
            "computed_class": "derived_in_gate",
        },
        "epsilon_eff_magnitude": {
            "tags": ["native_nullspace_rank_audit"],
            "computed_class": (
                "deferred_branch_data"
                if rank["native_null_moves_return"]
                else "derived_in_gate"
            ),
        },
        "pathA29_premise": {
            "tags": ["Z_is_premise", "boundary_dof_none"],
            "computed_class": "cited_earned_input",
        },
        "pathA28_side_conditions": {
            "tags": ["stage008_R0_minus_M0", "stage008_R1_minus_D1"],
            "computed_class": "external_bridge_input",
        },
        "stage022_gate4_result": {
            "tags": ["quad_regression_false", "stage022_done"],
            "computed_class": "cited_earned_input",
        },
        "time_convention": {
            "tags": ["exp_minus_i_omega_t"],
            "computed_class": "convention",
        },
    }
    if mutation.assert_not_derive:
        # Rewired to a genuinely stage023-derived object, not the cited fingerprints.
        items["forward_T0_T1_transfer_map"]["emitted_class"] = "asserted_literal"
    if mutation.emit_epsilon_magnitude_as_derived:
        epsilon_item = items["epsilon_eff_magnitude"]
        epsilon_item["emitted_class"] = "derived_in_gate"
        epsilon_item["concrete_value"] = sp.Rational(1, 2)
    selector = rank["selector_provenance"]
    if selector["tautological"]:
        items["branch_selector"] = {
            "tags": ["selector_without_named_pde"],
            "computed_class": "counterfactual_control",
            "emitted_class": "asserted_literal",
        }
    for item in items.values():
        item.setdefault("emitted_class", item["computed_class"])
        item["class_matches_computed"] = item["emitted_class"] == item["computed_class"]
    return {
        "items": items,
        "ok": all(item["class_matches_computed"] for item in items.values()),
        "baseline_emits_concrete_epsilon_magnitude": any(
            "concrete_value" in item
            for name, item in items.items()
            if name == "epsilon_eff_magnitude"
        ),
    }


def detect_decoupling(port: dict[str, sp.Expr], transfers: dict[str, Any]) -> dict[str, bool]:
    moves = {
        "gain0_moves_T0": not bool_zero(sp.diff(transfers["T0"], gain0)),
        "gain0_moves_T1": not bool_zero(sp.diff(transfers["T1"], gain0)),
        "gain1_moves_T0": not bool_zero(sp.diff(transfers["T0"], gain1)),
        "gain1_moves_T1": not bool_zero(sp.diff(transfers["T1"], gain1)),
    }
    p0_unaffected = gain0 not in port["P0_raw"].free_symbols and gain1 not in port["P0_raw"].free_symbols
    independently_moves = (
        moves["gain0_moves_T0"]
        and not moves["gain0_moves_T1"]
        and moves["gain1_moves_T1"]
        and not moves["gain1_moves_T0"]
    )
    return {**moves, "P0_unaffected": p0_unaffected, "decoupled": p0_unaffected and independently_moves}


class TrackingGateMap(dict[str, bool]):
    def __init__(self, values: dict[str, bool]):
        super().__init__(values)
        self.reads: set[str] = set()

    def __getitem__(self, key: str) -> bool:
        self.reads.add(key)
        return super().__getitem__(key)


def base_verdict(conditions: dict[str, bool]) -> tuple[str, frozenset[str]]:
    gates = TrackingGateMap(conditions)
    ordered_failures = (
        ("decoupled", FAIL_DECOUPLED),
        ("tautological", FAIL_TAUTOLOGICAL),
        ("dimensional", FAIL_DIMENSIONAL),
        ("no_consistent_return", FAIL_NO_CONSISTENT_RETURN),
        ("overcancel", FAIL_OVERCANCEL),
        ("epsilon_mismatch", FAIL_EPSILON_MISMATCH),
        ("underdetermined", FAIL_UNDERDETERMINED),
    )
    for key, label in ordered_failures:
        if gates[key]:
            return label, frozenset(gates.reads)
    return CROSS_L_RESIDUAL_PREDICTION, frozenset(gates.reads)


@lru_cache(maxsize=None)
def run_gate(mutation: Mutation = Mutation()) -> dict[str, Any]:
    port = build_port_kernel(mutation)
    transfers = build_transfers(mutation)
    residuals = build_residuals(transfers, mutation)
    rank = build_rank_audit(port, transfers, mutation)
    dimension = run_dimension_check(
        transfers,
        residuals,
        port,
        corrupt_sourced=mutation.corrupt_dimension,
    )
    provenance = build_provenance(mutation, rank)
    decoupling = detect_decoupling(port, transfers)
    conditions = {
        "decoupled": decoupling["decoupled"],
        "tautological": not provenance["ok"] or rank["selector_provenance"]["tautological"],
        "dimensional": not dimension["dimensional_ok"],
        "no_consistent_return": transfers["no_consistent_return"],
        "overcancel": transfers["overcancel"],
        "epsilon_mismatch": (
            not transfers["positive_bounded"]
            or not residuals["A0_form"]
            or not residuals["A1_form"]
            or residuals["A0_order"] != 1
            or residuals["A1_order"] != 3
        ),
        "underdetermined": rank["native_null_moves_return"],
    }
    verdict, read_set = base_verdict(conditions)
    return {
        "mutation": mutation,
        "port": port,
        "transfers": transfers,
        "residuals": residuals,
        "rank": rank,
        "dimension": dimension,
        "provenance": provenance,
        "decoupling": decoupling,
        "conditions": conditions,
        "verdict": verdict,
        "verdict_read_set": read_set,
    }


def dynamic_ablation(mutated: Mutation, expected_fail: str) -> dict[str, Any]:
    with_context = run_gate(mutated)
    without_context = run_gate(Mutation())
    neutral_mutation = Mutation(
        name=f"neutralized_{mutated.name}",
        inert_rerun_token="neutralized_independent_rerun",
    )
    neutral_context = run_gate(neutral_mutation)
    trace = (with_context["verdict"], without_context["verdict"])
    return {
        "verdict_trace": trace,
        "with_mutation": trace[0],
        "without_mutation": trace[1],
        "computed_verdicts_differ": trace[0] != trace[1],
        "expected_fail_fires": trace[0] == expected_fail and trace[1] != expected_fail,
        "neutralized_verdict": neutral_context["verdict"],
        "neutralized_mutation_fires": neutral_context["verdict"] == expected_fail,
        "neutralized_independent_rerun": (
            neutral_context is not without_context
            and neutral_context["mutation"] == neutral_mutation
            and neutral_context["mutation"] != without_context["mutation"]
        ),
    }


def port_dependency_ablation(corrupt: bool) -> dict[str, bool]:
    baseline = run_gate(Mutation())
    candidate = run_gate(Mutation(name="R1_corrupt_port", corrupt_port_kernel=corrupt))
    base_row = baseline["rank"]["constraint_jacobian"][0]
    candidate_row = candidate["rank"]["constraint_jacobian"][0]
    return {
        "P0_changes": not bool_zero(candidate["port"]["P0_raw"] - baseline["port"]["P0_raw"]),
        "ell2_row_changes": candidate_row != base_row,
    }


def check_failure_ablation(label: str, mutation: Mutation, expected: str) -> bool:
    audit = dynamic_ablation(mutation, expected)
    print(f"  {label}: verdict_trace={audit['verdict_trace']}")
    expect_zero(
        f"{label} dynamic rerun with mutation reaches {expected}",
        verdict_residual(audit["with_mutation"], expected),
    )
    expect_zero(
        f"{label} dynamic rerun without mutation returns earned native-nullspace FAIL",
        verdict_residual(audit["without_mutation"], FAIL_UNDERDETERMINED),
    )
    expect_bool(f"{label} two computed verdicts differ", audit["computed_verdicts_differ"])
    expect_bool(f"{label} own expected token fires", audit["expected_fail_fires"])
    expect_bool(
        f"{label} independently rerun neutralized mutation does not fire",
        audit["neutralized_independent_rerun"]
        and not audit["neutralized_mutation_fires"]
        and audit["neutralized_verdict"] == audit["without_mutation"],
    )
    return bool(
        audit["expected_fail_fires"]
        and audit["neutralized_independent_rerun"]
        and not audit["neutralized_mutation_fires"]
    )


def run_native_rank_and_selector() -> None:
    baseline = run_gate(Mutation())
    rank = baseline["rank"]
    subbanner("Genuine native symbolic Jacobian rank and return-moving nullspace")
    print(f"  GENERATOR_DOFS={list(map(fmt, rank['dofs']))}")
    print(f"  symbolic constraint Jacobian shape={sp.Matrix(rank['constraint_jacobian']).shape}")
    print(f"  rank0={rank['rank0']}")
    print(f"  native_nullspace_dimension={rank['native_nullspace_dimension']}")
    print(f"  return_augmented_rank={rank['return_augmented_rank']}")
    print(f"  return_moving_nullity={rank['return_moving_nullity']}")
    print(f"  native_null_moves_return={rank['native_null_moves_return']}")
    expect_bool("constraint Jacobian has exactly 3 genuine named rows (no padding)", len(rank["constraint_jacobian"]) == 3)
    expect_bool("constraint Jacobian has 11 genuine generator columns", all(len(row) == 11 for row in rank["constraint_jacobian"]))
    expect_zero("genuine sp.Matrix(rows).rank() gives rank0=3", rank["rank0"] - 3)
    expect_zero("computed native nullspace dimension is 11-rank0=8", rank["native_nullspace_dimension"] - 8)
    expect_zero("computed return-augmented rank is 5", rank["return_augmented_rank"] - 5)
    expect_zero("verdict-bearing return-moving nullity is 2", rank["return_moving_nullity"] - 2)
    expect_bool("native nullspace moves the return map", rank["native_null_moves_return"])
    print(
        "  PROVENANCE (documentation; not counted): "
        f"pathA_29 Z_is_premise={rank['pathA29_premise']['Z_is_premise']}, "
        f"boundary_dof={rank['pathA29_premise']['boundary_dof']}"
    )
    expect_bool("both Z-return unit-vector witnesses are constructed", len(rank["example_return_moving_directions"]) == 2)
    for witness in rank["example_return_moving_directions"]:
        constraint_residuals = [
            compact(sum(row_i * vec_i for row_i, vec_i in zip(row, witness["unit_vector"])))
            for row in rank["constraint_jacobian"]
        ]
        preserves_every_constraint = all(bool_zero(value) for value in constraint_residuals)
        print(
            f"  witness {witness['dof']}: preserves_every_constraint={preserves_every_constraint}, "
            f"delta_T0={fmt(witness['delta_T0'])}, delta_T1={fmt(witness['delta_T1'])}"
        )
        expect_bool(
            f"{witness['dof']} appended unit vector preserves every computed Jacobian row",
            preserves_every_constraint,
        )
        expect_bool(f"{witness['dof']} unit direction moves T0/T1", witness["moves_return"])

    subbanner("Isolated rank teeth and counterfactual selector control")
    injected = run_gate(Mutation(name="inject_null", inject_null=True))
    one = run_gate(Mutation(name="one_selector", selector="one_counterfactual_control"))
    selector = run_gate(Mutation(name="selector_control", selector="counterfactual_control"))
    neutered = run_gate(Mutation(name="selector_neutered", selector="neutered"))
    print(
        "  isolated ranks: "
        f"inject_null native={injected['rank']['native_nullspace_dimension']}; "
        f"one-row return-moving={one['rank']['return_moving_nullity']}; "
        f"two-row native={selector['rank']['native_nullspace_dimension']}, "
        f"return-moving={selector['rank']['return_moving_nullity']}"
    )
    print(f"  selector_rank0={selector['rank']['rank0']}")
    print(f"  selector_native_nullspace_dimension={selector['rank']['native_nullspace_dimension']}")
    print(f"  selector_return_moving_nullity={selector['rank']['return_moving_nullity']}")
    expect_zero("inject_null adds eta_null and changes native nullity 8->9", injected["rank"]["native_nullspace_dimension"] - 9)
    expect_zero("one selector row raises genuine rank 3->4", one["rank"]["rank0"] - 4)
    expect_zero("one selector row changes return-moving nullity 2->1", one["rank"]["return_moving_nullity"] - 1)
    expect_bool("one selector row leaves native_moves True", one["rank"]["native_null_moves_return"])
    expect_zero("two independent selector rows raise genuine rank 3->5", selector["rank"]["rank0"] - 5)
    expect_zero("two selector rows change native nullity 8->6", selector["rank"]["native_nullspace_dimension"] - 6)
    expect_zero("two selector rows collapse return-moving nullity 2->0", selector["rank"]["return_moving_nullity"])
    expect_bool("two selector rows flip native_moves True->False", rank["native_null_moves_return"] and not selector["rank"]["native_null_moves_return"])
    expect_zero("selector makes T0=1/2", selector["transfers"]["T0"] - sp.Rational(1, 2))
    expect_zero("selector makes T1=1/2", selector["transfers"]["T1"] - sp.Rational(1, 2))
    expect_zero("selector control reaches CROSS_L_RESIDUAL_PREDICTION", verdict_residual(selector["verdict"], CROSS_L_RESIDUAL_PREDICTION))
    expect_zero("baseline reaches earned FAIL_UNDERDETERMINED_NOT_PREDICTIVE", verdict_residual(baseline["verdict"], FAIL_UNDERDETERMINED))
    expect_bool("selector control changes the two computed verdicts", selector["verdict"] != baseline["verdict"])
    expect_bool("neutered selector has no rank or verdict effect", neutered["rank"]["return_moving_nullity"] == 2 and neutered["verdict"] == baseline["verdict"])
    provenance = selector["rank"]["selector_provenance"]
    expect_bool(
        "selector is counterfactual-control, not derived and not tautological",
        provenance["control_only"] and not provenance["derived_from_named_pde"] and not provenance["tautological"],
    )


def run_residual_and_consumption_teeth() -> tuple[bool, bool]:
    baseline = run_gate(Mutation())
    residuals = baseline["residuals"]
    subbanner("Forward A0/A1 residuals versus independent pathA_29 form")
    print(f"  cited stage022 coefficients: v0={fmt(residuals['v0'])}, v1={fmt(residuals['v1'])}")
    print(f"  A0={fmt(residuals['A0'])}")
    print(f"  A1={fmt(residuals['A1'])}")
    print(f"  A0_residual_to_pathA29_form={fmt(residuals['A0_residual'])}")
    print(f"  A1_residual_to_pathA29_form={fmt(residuals['A1_residual'])}")
    print(f"  A0_order=omega^{residuals['A0_order']}; A1_order=omega^{residuals['A1_order']}")
    expect_zero("forward A0 matches independently built pathA_29 form", residuals["A0_residual"])
    expect_zero("forward A1 matches independently built pathA_29 form", residuals["A1_residual"])
    expect_zero("A0 realized omega order is 1", residuals["A0_order"] - 1)
    expect_zero("A1 realized omega order is 3", residuals["A1_order"] - 3)
    report_self_consistency(
        "T/epsilon identities are forward from Z/K",
        baseline["transfers"]["forward_relations_ok"],
    )

    corrupt_v1 = run_gate(Mutation(name="corrupt_v1", corrupt_v1=True))
    corrupt_order = run_gate(Mutation(name="corrupt_a1_order", corrupt_a1_order=True))
    expect_fail("stage022 v1 consumption tooth 1/2->1/3 fires A1 form", corrupt_v1["residuals"]["A1_residual"])
    expect_fail("A1 order-realization corruption omega^3->omega^2 fires", sp.Integer(corrupt_order["residuals"]["A1_order"] - 3))
    v1_flag = check_failure_ablation("v1 consumption", Mutation(name="corrupt_v1", corrupt_v1=True), FAIL_EPSILON_MISMATCH)
    order_flag = check_failure_ablation("A1 order", Mutation(name="corrupt_a1_order", corrupt_a1_order=True), FAIL_EPSILON_MISMATCH)
    return v1_flag, order_flag


def run_dimensional_gate() -> bool:
    baseline = run_gate(Mutation())
    dim = baseline["dimension"]
    sourced = run_gate(Mutation(name="3f_corrupt_dimension", corrupt_dimension=True))
    free = run_dimension_check(
        baseline["transfers"], baseline["residuals"], baseline["port"], corrupt_free_carrier=True
    )
    subbanner("Exact (L,M,T) dimensional gate and free-carrier control")
    for name in ("A0", "A1", "T0", "T1", "epsilon0", "epsilon1", "P0_physical"):
        print(f"  [{name}]={dim['computed'][name]} expected={dim['expected'][name]}")
        expect_bool(f"[{name}] matches its sourced expected dimension", dim["computed"][name] == dim["expected"][name])
    print(f"  dimensional_ok={dim['dimensional_ok']}")
    print(f"  corrupt_sourced_verdict={sourced['dimension']['verdict']}")
    print(f"  corrupt_free_carrier_verdict={free['verdict']}")
    expect_bool("baseline dimensional_ok is True", dim["dimensional_ok"])
    expect_zero("corrupt sourced [M0]+=L reaches FAIL_DIMENSIONAL", verdict_residual(sourced["verdict"], FAIL_DIMENSIONAL))
    expect_zero("corrupt free q_free dimension reaches NO_FAIL locally", verdict_residual(free["verdict"], NO_FAIL))
    expect_bool("q_free appears in no checked expression", not free["checked_expressions_mention_q_free"])
    return check_failure_ablation(
        "3f corrupt sourced dimension",
        Mutation(name="3f_corrupt_dimension", corrupt_dimension=True),
        FAIL_DIMENSIONAL,
    )


def run_firewall() -> dict[str, bool]:
    baseline = run_gate(Mutation())
    subbanner("Strengthened provenance firewall and anti-back-solve teeth")
    provenance = baseline["provenance"]
    expect_bool("baseline class_matches_computed holds for every item", provenance["ok"])
    print(
        "  PROVENANCE (documentation; not counted): stage022 {1,1/2} class="
        f"{provenance['items']['stage022_v0_v1_coefficients']['computed_class']}"
    )
    print(
        "  PROVENANCE (documentation; not counted): ell2_P0_map tags="
        f"{provenance['items']['ell2_P0_map']['tags']}"
    )
    expect_zero(
        "epsilon magnitude computed class is deferred_branch_data",
        verdict_residual(
            provenance["items"]["epsilon_eff_magnitude"]["computed_class"],
            "deferred_branch_data",
        ),
    )
    expect_bool(
        "baseline structurally emits no concrete epsilon magnitude/value",
        not provenance["baseline_emits_concrete_epsilon_magnitude"],
    )

    assert_mut = run_gate(Mutation(name="3g_assert_not_derive", assert_not_derive=True))
    expect_bool(
        "assert_not_derive is rewired to the genuinely-023-derived transfer map",
        not assert_mut["provenance"]["items"]["forward_T0_T1_transfer_map"]["class_matches_computed"]
        and assert_mut["provenance"]["items"]["stage022_v0_v1_coefficients"]["class_matches_computed"],
    )
    epsilon_mut = run_gate(
        Mutation(name="emit_epsilon_magnitude_as_derived", emit_epsilon_magnitude_as_derived=True)
    )
    expect_bool(
        "epsilon-emission mutation keeps computed deferred class but emits a concrete value",
        epsilon_mut["provenance"]["items"]["epsilon_eff_magnitude"]["computed_class"] == "deferred_branch_data"
        and "concrete_value" in epsilon_mut["provenance"]["items"]["epsilon_eff_magnitude"],
    )
    asserted_selector = run_gate(Mutation(name="asserted_selector", selector="asserted_unproven"))
    expect_bool(
        "asserted selector is separately marked tautological",
        asserted_selector["rank"]["selector_provenance"]["tautological"],
    )
    flags = {
        "assert_not_derive": check_failure_ablation(
            "3g assert_not_derive",
            Mutation(name="3g_assert_not_derive", assert_not_derive=True),
            FAIL_TAUTOLOGICAL,
        ),
        "asserted_selector": check_failure_ablation(
            "3g asserted-selector",
            Mutation(name="3g_asserted_selector", selector="asserted_unproven"),
            FAIL_TAUTOLOGICAL,
        ),
        "emit_epsilon_magnitude_as_derived": check_failure_ablation(
            "3g emit_epsilon_magnitude_as_derived",
            Mutation(name="3g_emit_epsilon", emit_epsilon_magnitude_as_derived=True),
            FAIL_TAUTOLOGICAL,
        ),
    }
    return flags


def run_transfer_and_port_probes() -> dict[str, bool]:
    subbanner("Dynamic transfer/return-family, decoupling, and R1 probes")
    flags = {
        "3a": check_failure_ablation(
            "3a decouple_knobs", Mutation(name="3a_decouple", decouple_knobs=True), FAIL_DECOUPLED
        ),
        "3c": check_failure_ablation(
            "3c wrong_sign_return", Mutation(name="3c_wrong_sign", wrong_sign_return=True), FAIL_EPSILON_MISMATCH
        ),
        "3d": check_failure_ablation(
            "3d perfect_return", Mutation(name="3d_perfect", perfect_return=True), FAIL_OVERCANCEL
        ),
        "3h": check_failure_ablation(
            "3h no_consistent_return",
            Mutation(name="3h_no_consistent", no_consistent_return=True),
            FAIL_NO_CONSISTENT_RETURN,
        ),
    }
    r1_mutated = port_dependency_ablation(True)
    r1_neutral = port_dependency_ablation(False)
    print(f"  R1 corrupt-port signals={r1_mutated}; neutralized={r1_neutral}")
    expect_bool("R1 corrupt OmegaU->2*OmegaU changes P0_raw", r1_mutated["P0_changes"])
    expect_bool("R1 corrupt OmegaU->2*OmegaU changes ell=2 rank row", r1_mutated["ell2_row_changes"])
    expect_bool("R1 neutralized port mutation changes neither object", not any(r1_neutral.values()))
    flags["R1"] = all(r1_mutated.values()) and not any(r1_neutral.values())
    return flags


def run_verdict_seam_and_scope(all_probe_flags: dict[str, bool]) -> None:
    baseline = run_gate(Mutation())
    selector = run_gate(Mutation(name="selector_control", selector="counterfactual_control"))
    subbanner("022/023 verdict seam, computed read-set, and PROVENANCE consumption")
    expected_reads = {
        "decoupled",
        "tautological",
        "dimensional",
        "no_consistent_return",
        "overcancel",
        "epsilon_mismatch",
        "underdetermined",
    }
    print(f"  baseline conditions={baseline['conditions']}")
    print(f"  baseline verdict read-set={sorted(baseline['verdict_read_set'])}")
    expect_bool("base verdict reads exactly the stage023 physics ladder", baseline["verdict_read_set"] == expected_reads)
    expect_bool("quad_regression rung is absent (stage022 result is cited only)", "quad_regression" not in baseline["conditions"])
    expect_bool("inert able_to_fail_bad rung is absent", "able_to_fail_bad" not in baseline["conditions"])
    expect_bool("actual probe results non-recursively establish able-to-fail", all(all_probe_flags.values()))
    expect_zero("baseline physics verdict is earned underdetermination FAIL", verdict_residual(baseline["verdict"], FAIL_UNDERDETERMINED))
    expect_zero("selector control physics verdict is residual prediction", verdict_residual(selector["verdict"], CROSS_L_RESIDUAL_PREDICTION))
    helper_names = set(globals())
    forbidden_helpers = {
        "build_fingerprints",
        "build_gate4_non_regression",
        "build_prefactor",
        "yaml_write",
        "yaml_read",
    }
    expect_bool("no 022 fingerprint/Gate-4 or 019 prefactor engine is rebuilt", helper_names.isdisjoint(forbidden_helpers))
    live_symbol_names = {symbol.name for symbol in GENERATOR_DOFS}
    expect_bool("vestigial Delta and Sport symbols are dropped", live_symbol_names.isdisjoint({"Delta", "S_port", "Sport"}))
    print("  CHECKABLE: stage022 cited {v0=1,v1=1/2} feeds the forward A0/A1 map; v1 corruption fires.")
    print("  CHECKABLE: stage009/010 pathA_29 epsilon/(1+epsilon) residual form is built independently.")
    print("  PROVENANCE: pathA_29 Z_is_premise=True, boundary_dof=none is the keystone return premise.")
    print("  PROVENANCE: stage008 M0/D1 and R0=-M0/R1=-D1; stage017 grouped-P2 P0_raw port kernel.")
    print("  PROVENANCE CONTEXT: pathA_34-convention K0c/K1 effective stiffnesses; stage013/017 scalar reduction remains pending.")
    print("  PROVENANCE: c_s from R1 and a from CONV; stage022 Gate-4 result consumed as quad_regression=False, not rebuilt.")
    print("  ZERO file I/O: standalone print-only audit; no scratch YAML, report, note, card, LaTeX, or registration write.")


def print_verdict_labels() -> None:
    baseline = run_gate(Mutation())
    selector = run_gate(Mutation(name="selector_control", selector="counterfactual_control"))
    subbanner("Verdict labels:")
    print("  earned physics: the Gate-5 named constraints leave two return-moving directions.")
    print("  able-to-fail control: two counterfactual return equations collapse return-moving nullity to zero.")
    print("  honest export: Gate 6 must supply two independent return equations or an equivalent return law.")
    print(f"baseline_verdict={baseline['verdict']}")
    print(f"selector_control_verdict={selector['verdict']}")
    print("AUDIT_STATUS=PASS")
    print("PHYSICS_VERDICT=FAIL_UNDERDETERMINED_NOT_PREDICTIVE (2/2, COMPLETING)")


def main() -> int:
    banner("ledger_stage023_nullspace_underdetermination_sympy_audit")
    print("Target stem confirmed: ledger_stage023_nullspace_underdetermination")
    print("Engine: genuine SymPy symbolic Jacobian Matrix.rank; exact, float-free, standalone, zero file I/O.")
    run_native_rank_and_selector()
    v1_flag, order_flag = run_residual_and_consumption_teeth()
    dimension_flag = run_dimensional_gate()
    firewall_flags = run_firewall()
    transfer_flags = run_transfer_and_port_probes()
    all_flags = {
        "v1_consumption": v1_flag,
        "A1_order": order_flag,
        "3f": dimension_flag,
        **firewall_flags,
        **transfer_flags,
    }
    run_verdict_seam_and_scope(all_flags)
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
