#!/usr/bin/env python3
"""Ledger stage021 SymPy audit: mu_hat0-free dimensional closure.

Standalone, print-only, no arguments, and no file I/O.  This is the pathA_33
II-G4d dimensional-vector slice only: it computes the dimensionless physical
P0 gate from sourced N0/D0 dimensions, catches the dropped-frequency and
corrupt-port mutations, rejects the vacuous back-solved-mu_hat0 verdict, and
lands the completing 4/4 QUAD_CALIBRATED fold.  Adjacent stages are cited as
provenance only.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import sympy as sp


PASS_COUNT = 0
FAIL_COUNT = 0

DIMENSIONAL_OK = "DIMENSIONAL_OK"
FAIL_DIMENSIONAL = "FAIL_DIMENSIONAL"
NO_FAIL = "NO_FAIL"
QUAD_CALIBRATED = "QUAD_CALIBRATED"


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
    clean = sp.sympify(expr)
    floats = clean.atoms(sp.Float)
    if floats:
        message = f"{name}: Float atom(s) found in exact audit expression: {floats}"
        _record_fail(f"FAIL  {message}")
        raise AuditFailure(message)


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


def verdict_residual(actual: str, expected: str) -> sp.Integer:
    return sp.Integer(0) if actual == expected else sp.Integer(1)


a, c_s, c, G, omega = sp.symbols("a c_s c G omega", nonzero=True)
D0, N0 = sp.symbols("D0 N0", nonzero=True)
chi_Q, mtilde0, mu_hat0 = sp.symbols("chi_Q mtilde0 mu_hat0", nonzero=True)

Ldim, Mdim, Tdim = sp.symbols("L M T", positive=True)
Dim = tuple[sp.Rational, sp.Rational, sp.Rational]
ZERO_DIM: Dim = (sp.Rational(0), sp.Rational(0), sp.Rational(0))
SOURCED_N0_DIM: Dim = (-1, 1, 0)
SOURCED_D0_DIM: Dim = (-1, 1, -2)


def dim_add(left: Dim, right: Dim) -> Dim:
    return tuple(sp.Rational(a0 + b0) for a0, b0 in zip(left, right))  # type: ignore[return-value]


def dim_sub(left: Dim, right: Dim) -> Dim:
    return tuple(sp.Rational(a0 - b0) for a0, b0 in zip(left, right))  # type: ignore[return-value]


def dim_scale(dim: Dim, scale: sp.Rational) -> Dim:
    return tuple(sp.Rational(scale * d0) for d0 in dim)  # type: ignore[return-value]


def dim_of(expr: sp.Expr, symbol_dims: dict[sp.Symbol, Dim]) -> Dim:
    """Native recursive dimensional-vector evaluation in internal (L,M,T) order."""

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
        if any(dim != first for dim in dims):
            raise DimError(f"dimension mismatch in sum {expr}: {dims}")
        return first
    raise DimError(f"unsupported expression in dimension checker: {expr}")


def try_dim_of(expr: sp.Expr, symbol_dims: dict[sp.Symbol, Dim]) -> dict[str, Any]:
    """Turn DimError into data consumed by the named dimensional assertion."""

    try:
        return {"status": "OK", "dimension": dim_of(expr, symbol_dims), "error": None}
    except DimError as exc:
        return {"status": "DIM_ERROR", "dimension": None, "error": str(exc)}


class TrackingDimensions(dict[sp.Symbol, Dim]):
    """Record the dimension-map symbols genuinely read by dim_of."""

    def __init__(self, source: dict[sp.Symbol, Dim]) -> None:
        super().__init__(source)
        self.read_symbols: set[sp.Symbol] = set()

    def __getitem__(self, key: sp.Symbol) -> Dim:
        self.read_symbols.add(key)
        return super().__getitem__(key)


def dim_to_monomial(dim: Dim) -> sp.Expr:
    return compact(Ldim ** dim[0] * Mdim ** dim[1] * Tdim ** dim[2])


def exp_text(exp: sp.Rational) -> str:
    exp = sp.Rational(exp)
    if exp.q == 1:
        return str(exp.p)
    return f"{exp.p}/{exp.q}"


def dim_to_text(dim: Dim) -> str:
    # Internal order is (L,M,T); display order keeps the common L,T,M reading.
    parts: list[str] = []
    for name, exponent in (("L", dim[0]), ("T", dim[2]), ("M", dim[1])):
        exponent = sp.Rational(exponent)
        if exponent == 0:
            continue
        parts.append(name if exponent == 1 else f"{name}^{exp_text(exponent)}")
    return " ".join(parts) if parts else "1"


def dim_to_pretty_text(dim: Dim) -> str:
    superscript = str.maketrans("-0123456789/", "⁻⁰¹²³⁴⁵⁶⁷⁸⁹ᐟ")
    parts: list[str] = []
    for name, exponent in (("L", dim[0]), ("T", dim[2]), ("M", dim[1])):
        exponent = sp.Rational(exponent)
        if exponent == 0:
            continue
        parts.append(name if exponent == 1 else name + exp_text(exponent).translate(superscript))
    return " ".join(parts) if parts else "1"


def dim_record(name: str, expr: sp.Expr, symbol_dims: dict[sp.Symbol, Dim]) -> dict[str, Any]:
    dim = dim_of(expr, symbol_dims)
    return {
        "quantity": name,
        "expression": fmt(expr),
        "dimension": dim_to_text(dim),
        "dimension_monomial": fmt(dim_to_monomial(dim)),
        "dimension_vector_LMT": tuple(dim),
    }


def dim_vector_text(dim: Dim) -> tuple[str, str, str]:
    return tuple(fmt(value) for value in dim)  # type: ignore[return-value]


@dataclass(frozen=True)
class GateSignal:
    name: str
    value: bool
    upstream_reads: frozenset[str]


def compute_dimensional_gate(expr: sp.Expr, symbol_dims: dict[sp.Symbol, Dim]) -> dict[str, Any]:
    traced_dims = TrackingDimensions(symbol_dims)
    p0_dim = dim_of(expr, traced_dims)
    signal = GateSignal(
        name="dimensional_ok",
        value=p0_dim == ZERO_DIM,
        upstream_reads=frozenset({"p0_dim", *(symbol.name for symbol in traced_dims.read_symbols)}),
    )
    return {
        "dimension": p0_dim,
        "signal": signal,
        "symbol_reads": frozenset(symbol.name for symbol in traced_dims.read_symbols),
        "verdict_read_set": frozenset({signal.name, *signal.upstream_reads}),
    }


def scoped_dimensional_verdict(signal: GateSignal) -> str:
    return DIMENSIONAL_OK if signal.value else FAIL_DIMENSIONAL


def dynamic_dimensional_ablation(
    baseline_expr: sp.Expr,
    baseline_dims: dict[sp.Symbol, Dim],
    mutated_expr: sp.Expr,
    mutated_dims: dict[sp.Symbol, Dim],
) -> dict[str, Any]:
    """Dynamically rerun only the 021 dimensional gate with/without mutation."""

    with_gate = compute_dimensional_gate(mutated_expr, mutated_dims)
    without_gate = compute_dimensional_gate(baseline_expr, baseline_dims)
    verdict_trace: list[tuple[Dim, str]] = []

    def rerun_verdict(gate: dict[str, Any]) -> str:
        verdict = scoped_dimensional_verdict(gate["signal"])
        verdict_trace.append((gate["dimension"], verdict))
        return verdict

    with_mutation = rerun_verdict(with_gate)
    without_mutation = rerun_verdict(without_gate)
    rerun_gate_logic = (
        len(verdict_trace) == 2
        and len({dimension for dimension, _ in verdict_trace}) == 2
        and tuple(verdict for _, verdict in verdict_trace)
        == (with_mutation, without_mutation)
        and with_mutation != without_mutation
    )
    return {
        "rerun_gate_logic": rerun_gate_logic,
        "verdict_trace": tuple(verdict_trace),
        "with_mutation": with_mutation,
        "without_mutation": without_mutation,
        "expected_fail": FAIL_DIMENSIONAL,
        "fail_suppressed": (
            with_mutation == FAIL_DIMENSIONAL
            and without_mutation == DIMENSIONAL_OK
            and with_mutation != without_mutation
        ),
    }


def build_corruption_case(
    name: str,
    p0_raw: sp.Expr,
    p0_physical: sp.Expr,
    dims: dict[sp.Symbol, Dim],
) -> dict[str, Any]:
    raw_dim = dim_of(p0_raw, dims)
    gate = compute_dimensional_gate(p0_physical, dims)
    verdict = NO_FAIL if gate["signal"].value else FAIL_DIMENSIONAL
    return {
        "name": name,
        "P0_raw_dimension": raw_dim,
        "P0_physical_dimension": gate["dimension"],
        "dimensional_ok": gate["signal"].value,
        "verdict": verdict,
        "gate": gate,
    }


def back_solve_mutant(
    p0_physical: sp.Expr,
    target_rhs: sp.Expr,
    dims: dict[sp.Symbol, Dim],
) -> dict[str, Any]:
    """The rejected v1 verdict: solve mu anew, then test its defining identity."""

    p0_dim = dim_of(p0_physical, dims)
    rhs_dim = dim_of(target_rhs, dims)
    mu_dim = dim_scale(dim_sub(rhs_dim, p0_dim), sp.Rational(1, 2))
    homogeneity_pass = dim_add(dim_scale(mu_dim, sp.Rational(2)), p0_dim) == rhs_dim
    signal = GateSignal(
        name="homogeneity_pass",
        value=homogeneity_pass,
        upstream_reads=frozenset({"mu_dim", "p0_dim", "rhs_dim"}),
    )
    return {
        "p0_dim": p0_dim,
        "rhs_dim": rhs_dim,
        "mu_dim": mu_dim,
        "homogeneity_pass": homogeneity_pass,
        "wired_verdict": NO_FAIL if signal.value else FAIL_DIMENSIONAL,
        "verdict_read_set": frozenset({signal.name, *signal.upstream_reads}),
    }


def build_dimensions() -> dict[str, Any]:
    raw_symbol_dims: dict[sp.Symbol, Dim] = {
        a: (1, 0, 0),
        c_s: (1, 0, -1),
        c: (1, 0, -1),
        G: (3, -1, -2),
        omega: (0, 0, -1),
        D0: SOURCED_D0_DIM,
        N0: SOURCED_N0_DIM,
        chi_Q: ZERO_DIM,
        mtilde0: ZERO_DIM,
    }

    # PROVENANCE-FROZEN local fixture from stage018; no fingerprint machinery.
    u2_sourced = a**2 / (9 * c_s**2)
    u4_sourced = 4 * a**4 / (81 * c_s**4)
    v5_sourced = a**5 / (27 * c_s**5)

    # P0_raw is stage019's provenance-frozen definition; dimensions are 021's.
    p0_raw = N0 / D0
    frequency_normalization = (c_s / a) ** 2
    p0_physical = frequency_normalization * p0_raw
    yhat_physical = (
        1
        + u2_sourced * omega**2
        + u4_sourced * omega**4
        + sp.I * v5_sourced * omega**5
    )
    yhat_power_mutation = (
        1
        + u2_sourced * omega**3
        + u4_sourced * omega**4
        + sp.I * v5_sourced * omega**5
    )
    gamma5 = chi_Q * p0_physical * a**5 / (27 * c_s**5)
    # PROVENANCE-FROZEN assembled target from stage020; diagnostic only.
    target_rhs = 54 * G * c_s**5 / (5 * a**5 * c**5)

    p0_raw_dim = dim_of(p0_raw, raw_symbol_dims)
    frequency_norm_dim = dim_of(frequency_normalization, raw_symbol_dims)
    real_gate = compute_dimensional_gate(p0_physical, raw_symbol_dims)
    dimensional_ok = real_gate["signal"].value

    drop_norm_dim = dim_of(p0_raw, raw_symbol_dims)
    drop_norm_ok = drop_norm_dim == ZERO_DIM
    drop_norm_verdict = NO_FAIL if drop_norm_ok else FAIL_DIMENSIONAL

    corrupted_dims: dict[str, dict[sp.Symbol, Dim]] = {}
    for name, symbol in (("corrupt_N0", N0), ("corrupt_D0", D0), ("corrupt_G", G), ("corrupt_c_s", c_s)):
        mutation = dict(raw_symbol_dims)
        mutation[symbol] = ZERO_DIM
        corrupted_dims[name] = mutation

    truth_table = {
        name: build_corruption_case(name, p0_raw, p0_physical, dims)
        for name, dims in corrupted_dims.items()
    }
    truth_table["correct"] = build_corruption_case(
        "correct", p0_raw, p0_physical, raw_symbol_dims
    )

    yhat_try = try_dim_of(yhat_physical, raw_symbol_dims)
    yhat_mutant_try = try_dim_of(yhat_power_mutation, raw_symbol_dims)
    yhat_ok = yhat_try["status"] == "OK" and yhat_try["dimension"] == ZERO_DIM

    rhs_dim = dim_of(target_rhs, raw_symbol_dims)
    p0_dim = real_gate["dimension"]
    mu_dim = dim_scale(dim_sub(rhs_dim, p0_dim), sp.Rational(1, 2))
    symbol_dims = dict(raw_symbol_dims)
    symbol_dims[mu_hat0] = mu_dim
    lhs = (mu_hat0 * mtilde0) ** 2 * p0_physical
    lhs_raw_mutation = (mu_hat0 * mtilde0) ** 2 * p0_raw
    lhs_dim = dim_of(lhs, symbol_dims)
    lhs_raw_dim = dim_of(lhs_raw_mutation, symbol_dims)
    required_p0_dim = dim_sub(rhs_dim, dim_scale(mu_dim, sp.Rational(2)))
    homogeneity_pass = lhs_dim == rhs_dim and p0_dim == required_p0_dim
    diagnostic_participates_in_verdict = any(
        name in real_gate["verdict_read_set"]
        for name in {"mu_hat0", "mu_dim", "homogeneity_pass"}
    )

    backsolve_mutants = {
        name: back_solve_mutant(p0_physical, target_rhs, dims)
        for name, dims in corrupted_dims.items()
    }

    drop_ablation = dynamic_dimensional_ablation(
        p0_physical, raw_symbol_dims, p0_raw, raw_symbol_dims
    )
    corrupt_n0_ablation = dynamic_dimensional_ablation(
        p0_physical,
        raw_symbol_dims,
        p0_physical,
        corrupted_dims["corrupt_N0"],
    )

    local_pass_guard = (
        dimensional_ok
        and not drop_norm_ok
        and not truth_table["corrupt_N0"]["dimensional_ok"]
        and yhat_ok
        and truth_table["corrupt_G"]["dimensional_ok"]
    )

    return {
        "raw_symbol_dims": raw_symbol_dims,
        "symbol_dims": symbol_dims,
        "fixture": {"u2": u2_sourced, "u4": u4_sourced, "v5": v5_sourced},
        "P0_raw": p0_raw,
        "frequency_normalization": frequency_normalization,
        "P0_physical": p0_physical,
        "Yhat_physical": yhat_physical,
        "Yhat_power_mutation": yhat_power_mutation,
        "Gamma5": gamma5,
        "target_rhs": target_rhs,
        "P0_raw_dimension": p0_raw_dim,
        "frequency_normalization_dimension": frequency_norm_dim,
        "P0_physical_dimension": p0_dim,
        "dimensional_ok": dimensional_ok,
        "dimensional_status": scoped_dimensional_verdict(real_gate["signal"]),
        "dimensional_gate_depends_on_mu_hat0": "mu_hat0" in real_gate["verdict_read_set"],
        "real_gate": real_gate,
        "truth_table": truth_table,
        "corrupted_dims": corrupted_dims,
        "Yhat_try": yhat_try,
        "Yhat_mutant_try": yhat_mutant_try,
        "Yhat_ok": yhat_ok,
        "mu_hat0_homogeneity_diagnostic": {
            "label": "non-able-to-fail (mu_hat0 free carrier)",
            "participates_in_verdict": diagnostic_participates_in_verdict,
            "mu_dim_solved_from_rhs_minus_p0_over_2": True,
            "mu_hat0_dimension": mu_dim,
            "lhs_dimension": lhs_dim,
            "rhs_dimension": rhs_dim,
            "required_P0_dimension": required_p0_dim,
            "lhs_raw_without_frequency_normalization_dimension": lhs_raw_dim,
            "homogeneity_pass": homogeneity_pass,
        },
        "backsolve_mutants": backsolve_mutants,
        "mutation_drop_frequency_normalization": {
            "dimension": drop_norm_dim,
            "dimensional_ok": drop_norm_ok,
            "verdict": drop_norm_verdict,
            "mutation_fires": drop_norm_verdict == FAIL_DIMENSIONAL,
            "self_ablation": drop_ablation,
        },
        "mutation_corrupt_N0_dimension": {
            **truth_table["corrupt_N0"],
            "sourced_N0_dimension": SOURCED_N0_DIM,
            "corrupted_N0_dimension": ZERO_DIM,
            "mutation_fires": truth_table["corrupt_N0"]["verdict"] == FAIL_DIMENSIONAL,
            "self_ablation": corrupt_n0_ablation,
        },
        "local_pass_guard": local_pass_guard,
        "dimension_table": [
            dim_record("N0", N0, symbol_dims),
            dim_record("D0", D0, symbol_dims),
            dim_record("P0_raw=N0/D0", p0_raw, symbol_dims),
            dim_record("frequency_normalization=(c_s/a)^2", frequency_normalization, symbol_dims),
            dim_record("P0_physical=(c_s/a)^2*N0/D0", p0_physical, symbol_dims),
            dim_record("Yhat_out_physical_expansion", yhat_physical, symbol_dims),
            dim_record("Gamma5=chi_Q*P0*a^5/(27*c_s^5)", gamma5, symbol_dims),
            dim_record("target_lhs=(mu_hat0*mtilde0)^2*P0", lhs, symbol_dims),
            dim_record("target_rhs=54*G*c_s^5/(5*a^5*c^5)", target_rhs, symbol_dims),
        ],
    }


def run_mu_free_gate(data: dict[str, Any]) -> None:
    subbanner("Mu-hat0-free physical P0 dimensional gate")
    raw_dim = data["P0_raw_dimension"]
    norm_dim = data["frequency_normalization_dimension"]
    p0_dim = data["P0_physical_dimension"]
    read_set = data["real_gate"]["verdict_read_set"]
    print(f"  [P₀_raw]={dim_to_pretty_text(raw_dim)}")
    print(f"  [(c_s/a)²]={dim_to_pretty_text(norm_dim)}")
    print(f"  [P₀^phys]={dim_to_pretty_text(p0_dim)}")
    print(f"  dimensional_ok={data['dimensional_ok']} (COMPUTED by dim_of)")
    print(f"  dimension symbols actually read = {sorted(data['real_gate']['symbol_reads'])}")
    print(f"  verdict read-set = {sorted(read_set)}")
    print(
        "  dimensional_gate_depends_on_mu_hat0="
        + fmt(data["dimensional_gate_depends_on_mu_hat0"])
    )
    expect_bool("sourced [N0] is exactly (-1,1,0) in (L,M,T) order", data["raw_symbol_dims"][N0] == SOURCED_N0_DIM)
    expect_bool("sourced [D0] is exactly (-1,1,-2) in (L,M,T) order", data["raw_symbol_dims"][D0] == SOURCED_D0_DIM)
    expect_bool("computed [P0_raw] is T^2", raw_dim == (0, 0, 2))
    expect_bool("computed [(c_s/a)^2] is T^-2", norm_dim == (0, 0, -2))
    expect_bool("computed [P0_physical] is dimensionless", p0_dim == ZERO_DIM)
    expect_bool("dimensional_ok is computed True from [P0_physical]", data["dimensional_ok"])
    expect_bool("mu_hat0 is absent from the raw gate dimension map", mu_hat0 not in data["raw_symbol_dims"])
    expect_bool(
        "runtime gate trace reads exactly a,c_s,N0,D0",
        data["real_gate"]["symbol_reads"] == {"a", "c_s", "N0", "D0"},
    )
    expect_bool(
        "computed dimensional gate does not depend on mu_hat0",
        not data["dimensional_gate_depends_on_mu_hat0"],
    )


def run_yhat_structured_catch(data: dict[str, Any]) -> None:
    subbanner("Yhat dimensionless check with structured DimError catch")
    print(f"  local stage018 fixture = {data['fixture']}")
    print(f"  Yhat structured result = {data['Yhat_try']}")
    print(f"  Yhat omega-power mutant structured result = {data['Yhat_mutant_try']}")
    expect_bool(
        "Yhat dimensionless named assert consumes structured status and computed dimension",
        data["Yhat_try"]["status"] == "OK"
        and data["Yhat_try"]["dimension"] == ZERO_DIM
        and data["Yhat_ok"],
    )
    expect_fail(
        "Yhat dimensionless named assert catches u2*omega^3 sum mismatch",
        sp.Integer(1) if data["Yhat_mutant_try"]["status"] == "DIM_ERROR" else sp.Integer(0),
    )


def run_truth_table_and_probes(data: dict[str, Any]) -> None:
    subbanner("Corrupt-dimension scope truth-table and 3d/3d-prime probes")
    truth = data["truth_table"]
    expected_verdicts = {
        "corrupt_N0": FAIL_DIMENSIONAL,
        "corrupt_D0": FAIL_DIMENSIONAL,
        "corrupt_G": NO_FAIL,
        "corrupt_c_s": FAIL_DIMENSIONAL,
        "correct": NO_FAIL,
    }
    for name in ("corrupt_N0", "corrupt_D0", "corrupt_G", "corrupt_c_s", "correct"):
        case = truth[name]
        print(
            f"  {name}: [P₀_raw]={dim_to_pretty_text(case['P0_raw_dimension'])}; "
            f"[P₀^phys]={dim_to_pretty_text(case['P0_physical_dimension'])}; "
            f"vector(L,M,T)={dim_vector_text(case['P0_physical_dimension'])}; {case['verdict']}"
        )
        expect_zero(
            f"truth-table {name} reaches {expected_verdicts[name]} at its own assert",
            verdict_residual(case["verdict"], expected_verdicts[name]),
        )
    expect_bool(
        "corrupt-N0 keeps normalization: [P0_raw]=(1,-1,2)",
        truth["corrupt_N0"]["P0_raw_dimension"] == (1, -1, 2),
    )
    expect_bool(
        "corrupt-N0 keeps normalization: [P0_physical]=(1,-1,0)",
        truth["corrupt_N0"]["P0_physical_dimension"] == (1, -1, 0),
    )
    p0_free_symbol_names = {symbol.name for symbol in data["P0_physical"].free_symbols}
    print(f"  free_symbols(P0_physical) = {sorted(p0_free_symbol_names)}")
    expect_bool(
        "dependency-scope control computes G not in free_symbols(P0_physical)",
        "G" not in p0_free_symbol_names,
    )
    observed = tuple(truth[name]["verdict"] for name in expected_verdicts)
    expected = tuple(expected_verdicts.values())
    print(
        "  SUMMARY (not counted): full corrupt-dim truth-table "
        f"observed={observed}; expected={expected}; matches={observed == expected}"
    )

    drop = data["mutation_drop_frequency_normalization"]
    corrupt_n0 = data["mutation_corrupt_N0_dimension"]
    print(f"  3d_dimensional_break: drop (c_s/a)² -> {drop['verdict']}")
    print(f"  3d self_ablation = {drop['self_ablation']}")
    print(f"  3d_prime_corrupt_port_dimension: corrupt [N₀] -> {corrupt_n0['verdict']}")
    print(f"  3d_prime self_ablation = {corrupt_n0['self_ablation']}")
    expect_bool("drop-frequency-normalization mutation fires", drop["mutation_fires"])
    expect_zero("3d drop normalization reaches FAIL_DIMENSIONAL", verdict_residual(drop["verdict"], FAIL_DIMENSIONAL))
    expect_bool("corrupt-N0 dimension mutation fires", corrupt_n0["mutation_fires"])
    expect_zero("3d-prime corrupt port dimension reaches FAIL_DIMENSIONAL", verdict_residual(corrupt_n0["verdict"], FAIL_DIMENSIONAL))
    for probe_name, ablation in (("3d", drop["self_ablation"]), ("3d-prime", corrupt_n0["self_ablation"])):
        expect_bool(f"{probe_name} self-ablation dynamically reruns 021-local gate logic", ablation["rerun_gate_logic"])
        expect_zero(f"{probe_name} dynamic rerun with mutation is FAIL_DIMENSIONAL", verdict_residual(ablation["with_mutation"], FAIL_DIMENSIONAL))
        expect_zero(f"{probe_name} dynamic rerun without mutation is DIMENSIONAL_OK", verdict_residual(ablation["without_mutation"], DIMENSIONAL_OK))
        expect_bool(f"{probe_name} dynamic self-ablation changes the local verdict", ablation["fail_suppressed"])


def run_anti_v1_and_diagnostic(data: dict[str, Any]) -> None:
    subbanner("Decisive anti-v1 read-set tooth and non-verdict mu_hat0 diagnostic")
    diagnostic = data["mu_hat0_homogeneity_diagnostic"]
    real_read_set = data["real_gate"]["verdict_read_set"]
    forbidden = {"mu_hat0", "mu_dim", "homogeneity_pass"}
    print(f"  real verdict read-set = {sorted(real_read_set)}")
    print(f"  [mu_hat0]={dim_to_pretty_text(diagnostic['mu_hat0_dimension'])}")
    print(f"  [lhs]={dim_to_pretty_text(diagnostic['lhs_dimension'])}")
    print(f"  [rhs]={dim_to_pretty_text(diagnostic['rhs_dimension'])}")
    print(f"  homogeneity_pass={diagnostic['homogeneity_pass']} DIAGNOSTIC")
    print(f"  participates_in_verdict={diagnostic['participates_in_verdict']}")
    expect_bool(
        "computed real verdict read-set excludes mu_hat0/mu_dim/homogeneity_pass",
        real_read_set.isdisjoint(forbidden),
    )
    expect_bool(
        "real verdict read-set contains dimensional_ok and its computed p0_dim input",
        {"dimensional_ok", "p0_dim"} <= real_read_set,
    )
    expect_bool(
        "back-solved mu_hat0 has dimension (-1,-1/2,-1)",
        diagnostic["mu_hat0_dimension"] == (-1, sp.Rational(-1, 2), -1),
    )
    expect_bool("mu_hat0 homogeneity diagnostic computes True", diagnostic["homogeneity_pass"])
    expect_bool("mu_hat0 homogeneity diagnostic is explicitly non-verdict", diagnostic["participates_in_verdict"] is False)

    mutants = data["backsolve_mutants"]
    for name in ("corrupt_N0", "corrupt_D0", "corrupt_G", "corrupt_c_s"):
        mutant = mutants[name]
        print(
            f"  wired back-solve mutant {name}: homogeneity_pass="
            f"{mutant['homogeneity_pass']} -> {mutant['wired_verdict']}"
        )
        expect_zero(
            f"wired back-solve mutant stays NO_FAIL under {name}",
            verdict_residual(mutant["wired_verdict"], NO_FAIL),
        )
    mutant_read_set = frozenset().union(
        *(mutant["verdict_read_set"] for mutant in mutants.values())
    )
    wiring_policy_residual = sp.Integer(len(mutant_read_set & forbidden))
    expect_fail(
        "anti-v1 demotion policy rejects wiring homogeneity_pass/mu_dim into the verdict",
        wiring_policy_residual,
    )
    real_fail_count = sum(
        truth_case["verdict"] == FAIL_DIMENSIONAL
        for name, truth_case in data["truth_table"].items()
        if name != "correct"
    )
    mutant_fail_count = sum(
        mutant["wired_verdict"] == FAIL_DIMENSIONAL for mutant in mutants.values()
    )
    print(
        "  SUMMARY (not counted): real mu-free gate fail count="
        f"{real_fail_count}; wired back-solve mutant fail count={mutant_fail_count}"
    )


def run_scope_provenance_and_landing(data: dict[str, Any]) -> None:
    subbanner("021 scope, provenance-only consumption, and COMPLETING landing")
    expect_bool(
        "021-local pass guard excludes homogeneityPass and passes its dimensional teeth",
        data["local_pass_guard"],
    )
    print("  CONSUMES (PROVENANCE only): [N₀]=L⁻¹M from the density/continuity-port numerator; it genuinely enters the checked gate.")
    print("  CONSUMES (PROVENANCE only): [D₀]=L⁻¹T⁻²M from the carried reduced static conservative denominator D₀=K−B₀−Z₀.")
    print("  CONSUMES (PROVENANCE only): stage018 u₂/u₄/v₅ as a local frozen fixture for Yhat; no build_fingerprint machinery.")
    print("  CONSUMES (PROVENANCE only): stage019 P0=N0/D0 defines P₀_raw; no prefactor-series machinery.")
    print("  CONSUMES (PROVENANCE only): stage020 target_rhs enters only the non-verdict mu_hat0 diagnostic.")
    print("  018 fingerprint DONE; 019 prefactor DONE; 020 provenance partition + CALIBRATED classification DONE.")
    print("  QUAD_CALIBRATED (4/4) -- mu_hat0-free [P₀^phys]=1 dimensional closure EARNED (COMPLETING)")
    print("  joint QUAD_CALIBRATED COMPLETE = 018 DONE ∧ 019 DONE ∧ 020 DONE ∧ 021 dimensional closure EARNED.")
    print("  COMPLETE != PASS: token remains QUAD_CALIBRATED, not QUAD_PASS; 020 has G=GENUINE_BLOCKED.")
    print("  EXPORTS: completed Gate-4 dimensional closure to 022/027 and the Part-VII dimensional firewall.")
    print("  dropped-bookkeeping: scratch YAML, reports, feed writers, engine cross-reads, and all file I/O are absent.")


def print_verdict_labels() -> None:
    subbanner("Verdict labels")
    print("  ledger local gate label: DIMENSIONAL_OK (mu_hat0-free [P0_phys]=1 from sourced [N0]/[D0] and [(c_s/a)^2])")
    print("  source top-line verdict: QUAD_CALIBRATED (JOINT 4-stage; 021 is the COMPLETING 4/4 leg)")
    print("  joint composition: 018 fingerprint [DONE] AND 019 prefactor [DONE] AND 020 provenance partition + CALIBRATED classification [DONE] AND 021 dimensional closure [EARNED here] => COMPLETE")
    print("  COMPLETE != PASS: the joint token remains QUAD_CALIBRATED, not QUAD_PASS (020 provenance; G=GENUINE_BLOCKED)")
    print("  earned dimensional gate: [P0_raw]=T^2, [(c_s/a)^2]=T^-2, [P0_phys]=1; mu_hat0 is absent from raw_symbol_dims and the verdict read-set")
    print("  earned able-to-fail: 3d drop normalization and 3d-prime corrupt [N0] reach FAIL_DIMENSIONAL with DYNAMIC 021-local self-ablations; corrupt [G] is NO_FAIL scope control")
    print("  diagnostic only: back-solved [mu_hat0] and homogeneity_pass=True participate_in_verdict=False; the wired back-solve mutant is all-NO_FAIL and rejected as vacuous")


def main() -> int:
    banner("ledger_stage021_dimensional_closure_sympy_audit")
    print("Target stem confirmed: ledger_stage021_dimensional_closure")
    print("Engine: SymPy exact recursive (L,M,T) dimensional vectors; no floats/tolerances; zero file I/O.")
    data = build_dimensions()
    run_mu_free_gate(data)
    run_yhat_structured_catch(data)
    run_truth_table_and_probes(data)
    run_anti_v1_and_diagnostic(data)
    run_scope_provenance_and_landing(data)
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
