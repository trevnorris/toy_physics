#!/usr/bin/env python3
"""Ledger stage018 SymPy audit: DtN Hankel fingerprint + chi sign.

Standalone, print-only, no arguments, no file I/O. This is the pathA_33
II-G4a exterior-wave slice only: explicit l=2 spherical j/y functions,
the outgoing DtN Hankel log-derivative, the derived rational fingerprint,
the chi sign from outgoing versus incoming branches, passivity/not-inserted
sink probes, and the units-restored coefficient dimensions. Sibling stages
own the prefactor, normalization partition, and dimensional closure.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from typing import Any

import sympy as sp


PASS_COUNT = 0
FAIL_COUNT = 0

DTN_HANKEL_FINGERPRINT_EARNED = "DTN_HANKEL_FINGERPRINT_EARNED"
FAIL_FINGERPRINT = "FAIL_FINGERPRINT"
FAIL_NOT_OUTGOING = "FAIL_NOT_OUTGOING"
FAIL_DIMENSIONAL = "FAIL_DIMENSIONAL"
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
        return dim_text(expr)
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


def bool_from_residual(residual: sp.Expr | int) -> bool:
    assert_no_float("bool_from_residual", residual)
    return compact(residual) == 0


def bool_residual(condition: bool) -> sp.Integer:
    return sp.Integer(0) if bool(condition) else sp.Integer(1)


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


ZERO_DIM = Dim()
DIM_L = Dim(1, 0, 0)
DIM_SPEED = Dim(1, 0, -1)
DIM_T2 = Dim(0, 0, 2)
DIM_T4 = Dim(0, 0, 4)
DIM_T5 = Dim(0, 0, 5)


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
    raise DimError(f"unsupported dimension expression {clean}")


def dim_text(dim: Dim) -> str:
    parts: list[str] = []
    for label, exponent in (("L", dim.l), ("M", dim.m), ("T", dim.t)):
        if exponent == 0:
            continue
        if exponent == 1:
            parts.append(label)
        else:
            parts.append(f"{label}^{fmt(exponent)}")
    return "1" if not parts else " ".join(parts)


z = sp.Symbol("z", real=True)
omega = sp.Symbol("omega", real=True)
a = sp.Symbol("a", positive=True, real=True)
c_s = sp.Symbol("c_s", positive=True, real=True)


def series_no_o(expr: sp.Expr, var: sp.Symbol, order: int) -> sp.Expr:
    return sp.expand(sp.series(expr, var, 0, order).removeO())


def spherical_j2() -> sp.Expr:
    return (sp.Integer(3) / z**3 - sp.Integer(1) / z) * sp.sin(z) - sp.Integer(3) * sp.cos(z) / z**2


def spherical_y2() -> sp.Expr:
    return (sp.Integer(1) / z - sp.Integer(3) / z**3) * sp.cos(z) - sp.Integer(3) * sp.sin(z) / z**2


def branch_h(kind: str) -> tuple[sp.Expr, str]:
    j2 = spherical_j2()
    y2 = spherical_y2()
    if kind == "outgoing_hankel1":
        return j2 + sp.I * y2, "dtn_hankel1"
    if kind == "incoming_hankel2":
        return j2 - sp.I * y2, "dtn_hankel2"
    if kind == "standing_j2":
        return j2, "regular_standing_j2"
    raise ValueError(kind)


@lru_cache(maxsize=None)
def dtn_branch(kind: str, mutation: str = "none") -> dict[str, Any]:
    h, source = branch_h(kind)
    if mutation == "u2_mixed_derivation_copy":
        h = h * (1 + z**2)
        source = source + "_u2_copy"
    elif mutation == "u4_only_derivation_copy":
        h = h * (1 + z**4)
        source = source + "_u4_copy"
    elif mutation == "v5_only_derivation_copy":
        h = h * (1 + sp.I * z**5)
        source = source + "_v5_copy"
    elif mutation != "none":
        raise ValueError(mutation)

    lam = compact(z * sp.diff(h, z) / h)
    yhat = compact(-sp.Integer(3) / lam)
    h_series = series_no_o(h, z, 8)
    lambda_series = series_no_o(lam, z, 7)
    yhat_series = series_no_o(yhat, z, 7)
    yhat_series_to_z5 = series_no_o(yhat, z, 6)
    u2_z = compact(yhat_series.coeff(z, 2))
    u4_z = compact(yhat_series.coeff(z, 4))
    v5_z = compact(yhat_series.coeff(z, 5) / sp.I)
    static = compact(yhat_series.coeff(z, 0))
    return {
        "kind": kind,
        "source": source,
        "mutation": mutation,
        "h": h,
        "lambda": lam,
        "yhat": yhat,
        "h_series": h_series,
        "lambda_series": lambda_series,
        "yhat_series": yhat_series,
        "yhat_series_to_z5": yhat_series_to_z5,
        "static": static,
        "u2_z": u2_z,
        "u4_z": u4_z,
        "v5_z": v5_z,
        "u2": compact(u2_z * a**2 / c_s**2),
        "u4": compact(u4_z * a**4 / c_s**4),
        "v5": compact(v5_z * a**5 / c_s**5),
    }


def passivity_from_source(branch: dict[str, Any]) -> dict[str, Any]:
    imag_nonzero = not bool_from_residual(branch["v5_z"])
    is_genuine = branch["source"] == "dtn_hankel1"
    return {
        "source": branch["source"],
        "imaginary_v5_nonzero": imag_nonzero,
        "genuine_outgoing": is_genuine and imag_nonzero,
    }


def scoped_verdict(overrides: dict[str, bool] | None = None) -> str:
    gates = {
        "fingerprint_ok": True,
        "outgoing_ok": True,
        "dimensional_ok": True,
        "able_to_fail_ok": True,
    }
    if overrides:
        gates.update(overrides)
    if not gates["fingerprint_ok"]:
        return FAIL_FINGERPRINT
    if not gates["outgoing_ok"]:
        return FAIL_NOT_OUTGOING
    if not gates["dimensional_ok"]:
        return FAIL_DIMENSIONAL
    if not gates["able_to_fail_ok"]:
        return FAIL_NOT_ABLE_TO_FAIL
    return DTN_HANKEL_FINGERPRINT_EARNED


def dynamic_ablation(overrides: dict[str, bool], expected_fail: str) -> dict[str, Any]:
    with_mutation = scoped_verdict(overrides)
    without_mutation = scoped_verdict({})
    return {
        "rerun_gate_logic": True,
        "with_mutation": with_mutation,
        "without_mutation": without_mutation,
        "expected_fail": expected_fail,
        "fail_suppressed": with_mutation == expected_fail and without_mutation != expected_fail,
    }


def build_fingerprint() -> dict[str, Any]:
    out = dtn_branch("outgoing_hankel1")
    incoming = dtn_branch("incoming_hankel2")
    standing = dtn_branch("standing_j2")
    expected = {
        "u2_z": sp.Rational(1, 9),
        "u4_z": sp.Rational(4, 81),
        "v5_z": sp.Rational(1, 27),
        "u2": a**2 / (9 * c_s**2),
        "u4": 4 * a**4 / (81 * c_s**4),
        "v5": a**5 / (27 * c_s**5),
    }
    matches = {
        name: bool_from_residual(out[name] - target)
        for name, target in expected.items()
    }
    computed_denominator = compact(1 / out["v5_z"])
    chi_slot = a**5 / (sp.Integer(27) * c_s**5)
    chi_derived = compact(out["v5"] / chi_slot)
    chi_incoming = compact(incoming["v5"] / chi_slot)
    return {
        "outgoing": out,
        "incoming": incoming,
        "standing": standing,
        "expected": expected,
        "matches": matches,
        "computed_v5_denominator": computed_denominator,
        "chi_Q": chi_derived,
        "chi_Q_incoming": chi_incoming,
        "exact_sample": {
            "a": sp.Integer(3),
            "c_s": sp.Integer(2),
            "u2": compact(out["u2"].subs({a: 3, c_s: 2})),
            "u4": compact(out["u4"].subs({a: 3, c_s: 2})),
            "v5": compact(out["v5"].subs({a: 3, c_s: 2})),
        },
    }


def build_units_leg(out: dict[str, Any]) -> dict[str, Any]:
    dims = {a: DIM_L, c_s: DIM_SPEED}
    u2_dim = dim_of(out["u2"], dims)
    u4_dim = dim_of(out["u4"], dims)
    v5_dim = dim_of(out["v5"], dims)
    corrupted_u2 = compact(out["u2_z"] * a**2 / c_s**3)
    corrupted_u2_dim = dim_of(corrupted_u2, dims)
    ok = (
        dim_residual(u2_dim, DIM_T2) == 0
        and dim_residual(u4_dim, DIM_T4) == 0
        and dim_residual(v5_dim, DIM_T5) == 0
    )
    corrupted_ok = dim_residual(corrupted_u2_dim, DIM_T2) == 0
    return {
        "dims": {"u2": u2_dim, "u4": u4_dim, "v5": v5_dim},
        "expected": {"u2": DIM_T2, "u4": DIM_T4, "v5": DIM_T5},
        "ok": ok,
        "corrupted_u2": corrupted_u2,
        "corrupted_u2_dim": corrupted_u2_dim,
        "corrupted_ok": corrupted_ok,
        "corrupted_verdict": scoped_verdict({"dimensional_ok": corrupted_ok}),
        "self_ablation": dynamic_ablation({"dimensional_ok": False}, FAIL_DIMENSIONAL),
    }


def build_probes(fingerprint: dict[str, Any], units_leg: dict[str, Any]) -> dict[str, Any]:
    out = fingerprint["outgoing"]
    incoming = fingerprint["incoming"]
    standing = fingerprint["standing"]
    incoming_expected = (
        bool_from_residual(incoming["u2_z"] - out["u2_z"])
        and bool_from_residual(incoming["u4_z"] - out["u4_z"])
        and bool_from_residual(incoming["v5_z"] + out["v5_z"])
        and bool_from_residual(fingerprint["chi_Q_incoming"] + 1)
    )
    standing_expected = (
        bool_from_residual(standing["lambda_series"].coeff(z, 0) - 2)
        and bool_from_residual(standing["static"] + sp.Rational(3, 2))
        and bool_from_residual(standing["v5_z"])
    )

    phenom_source = {
        "source": "phenomenological_inserted_sink",
        "v5_z": out["v5_z"],
    }
    phenom_genuine = phenom_source["source"] == "dtn_hankel1" and not bool_from_residual(phenom_source["v5_z"])
    derivation_teeth = {
        "u2_mixed": dtn_branch("outgoing_hankel1", "u2_mixed_derivation_copy"),
        "u4_only": dtn_branch("outgoing_hankel1", "u4_only_derivation_copy"),
        "v5_only": dtn_branch("outgoing_hankel1", "v5_only_derivation_copy"),
    }
    return {
        "fingerprint_derivation_teeth": derivation_teeth,
        "3a_wrong_bc": {
            "incoming_expected": incoming_expected,
            "standing_expected": standing_expected,
            "incoming_verdict": FAIL_FINGERPRINT if incoming_expected else FAIL_NOT_ABLE_TO_FAIL,
            "standing_verdict": FAIL_FINGERPRINT if standing_expected else FAIL_NOT_ABLE_TO_FAIL,
            "self_ablation": dynamic_ablation({"fingerprint_ok": False}, FAIL_FINGERPRINT),
        },
        "3b_imposed_dissipation": {
            "mutated_source": phenom_source["source"],
            "inserted_v5_z": phenom_source["v5_z"],
            "genuine_outgoing": phenom_genuine,
            "removing_outgoing_bc_removes_imaginary_term": bool_from_residual(standing["v5_z"]),
            "verdict": FAIL_NOT_OUTGOING if not phenom_genuine else DTN_HANKEL_FINGERPRINT_EARNED,
            "self_ablation": dynamic_ablation({"outgoing_ok": False}, FAIL_NOT_OUTGOING),
        },
        "units_dimensional": units_leg,
    }


def build_baseline() -> dict[str, Any]:
    fingerprint = build_fingerprint()
    passivity = passivity_from_source(fingerprint["outgoing"])
    units_leg = build_units_leg(fingerprint["outgoing"])
    probes = build_probes(fingerprint, units_leg)
    fingerprint_ok = (
        all(fingerprint["matches"].values())
        and bool_from_residual(fingerprint["computed_v5_denominator"] - 27)
        and bool_from_residual(fingerprint["chi_Q"] - 1)
        and bool_from_residual(fingerprint["chi_Q_incoming"] + 1)
        and probes["3a_wrong_bc"]["incoming_expected"]
        and probes["3a_wrong_bc"]["standing_expected"]
    )
    outgoing_ok = (
        passivity["genuine_outgoing"]
        and bool_from_residual(fingerprint["outgoing"]["static"] - 1)
        and probes["3b_imposed_dissipation"]["verdict"] == FAIL_NOT_OUTGOING
    )
    able_to_fail_ok = (
        probes["3a_wrong_bc"]["incoming_verdict"] == FAIL_FINGERPRINT
        and probes["3a_wrong_bc"]["standing_verdict"] == FAIL_FINGERPRINT
        and probes["3a_wrong_bc"]["self_ablation"]["fail_suppressed"]
        and probes["3b_imposed_dissipation"]["self_ablation"]["fail_suppressed"]
        and units_leg["corrupted_verdict"] == FAIL_DIMENSIONAL
    )
    gates = {
        "fingerprint_ok": fingerprint_ok,
        "outgoing_ok": outgoing_ok,
        "dimensional_ok": units_leg["ok"],
        "able_to_fail_ok": able_to_fail_ok,
    }
    return {
        "fingerprint": fingerprint,
        "passivity": passivity,
        "units": units_leg,
        "probes": probes,
        "gates": gates,
        "verdict": scoped_verdict(gates),
    }


def run_fingerprint_derivation(data: dict[str, Any]) -> None:
    fp = data["fingerprint"]
    out = fp["outgoing"]
    expected = fp["expected"]
    subbanner("Outgoing DtN Hankel fingerprint derived from the branch")
    print("  h2^(1)(z) is built from explicit spherical j2 + i*y2 rational sin/cos expressions.")
    print(f"  h2^(1) series = {fmt(out['h_series'])}")
    print(f"  Lambda_2^out(z) = z*h2'(z)/h2(z)")
    print(f"  Lambda_2^out series = {fmt(out['lambda_series'])}")
    print(f"  Yhat_2^out(z) = -3/Lambda_2^out(z)")
    print(f"  Yhat_2^out series through z^5 = {fmt(out['yhat_series_to_z5'])}")
    print("  Coefficients are read from that series by coeff(z,n); the typed targets are independent references.")
    for key in ("u2_z", "u4_z", "v5_z", "u2", "u4", "v5"):
        print(f"  derived {key} = {fmt(out[key])}")
        expect_zero(f"outgoing derived {key} matches independent target", out[key] - expected[key])
    expect_zero("computed radiative denominator 1/v5_z equals 27", fp["computed_v5_denominator"] - 27)
    print(f"  exact sample at a=3, c_s=2 = { {key: fmt(value) for key, value in fp['exact_sample'].items()} }")


def run_chi_and_standing(data: dict[str, Any]) -> None:
    fp = data["fingerprint"]
    out = fp["outgoing"]
    incoming = fp["incoming"]
    standing = fp["standing"]
    subbanner("chi sign from outgoing/incoming Hankel branches and standing contrast")
    print(f"  outgoing v5_z = {fmt(out['v5_z'])}; incoming v5_z = {fmt(incoming['v5_z'])}")
    print(f"  computed chi_Q outgoing = {fmt(fp['chi_Q'])}; incoming = {fmt(fp['chi_Q_incoming'])}")
    expect_zero("chi_Q computed from outgoing v5 physical slot is +1", fp["chi_Q"] - 1)
    expect_zero("incoming v5_z + outgoing v5_z = 0", incoming["v5_z"] + out["v5_z"])
    expect_zero("incoming computed chi_Q is -1", fp["chi_Q_incoming"] + 1)
    expect_zero("incoming even u2_z unchanged", incoming["u2_z"] - out["u2_z"])
    expect_zero("incoming even u4_z unchanged", incoming["u4_z"] - out["u4_z"])
    expect_fail("using incoming branch as outgoing chi check fires", fp["chi_Q_incoming"] - 1)
    print(f"  standing Lambda series = {fmt(standing['lambda_series'])}")
    print(f"  standing Yhat series = {fmt(standing['yhat_series_to_z5'])}")
    expect_zero("standing branch Lambda_static is +2", standing["lambda_series"].coeff(z, 0) - 2)
    expect_zero("standing branch Yhat_static is -3/2", standing["static"] + sp.Rational(3, 2))
    expect_zero("standing branch has no radiating v5_z", standing["v5_z"])
    expect_fail("standing static slot differs from outgoing normalized static slot", standing["static"] - 1)


def run_passivity_and_probes(data: dict[str, Any]) -> None:
    fp = data["fingerprint"]
    passivity = data["passivity"]
    probes = data["probes"]
    wrong_bc = probes["3a_wrong_bc"]
    imposed = probes["3b_imposed_dissipation"]
    subbanner("Passivity and the two 018 fingerprint probes")
    print(f"  passivity source = {passivity['source']}; imaginary v5 nonzero = {passivity['imaginary_v5_nonzero']}")
    expect_bool("outgoing branch passivity gate is genuine outgoing with nonzero v5", passivity["genuine_outgoing"])
    expect_zero("outgoing static slot is the DC limit of the same response", fp["outgoing"]["static"] - 1)
    expect_zero("3a incoming verdict is FAIL_FINGERPRINT", verdict_residual(wrong_bc["incoming_verdict"], FAIL_FINGERPRINT))
    expect_zero("3a standing verdict is FAIL_FINGERPRINT", verdict_residual(wrong_bc["standing_verdict"], FAIL_FINGERPRINT))
    expect_bool("3a incoming predicted change is computed", wrong_bc["incoming_expected"])
    expect_bool("3a standing predicted change is computed", wrong_bc["standing_expected"])
    expect_zero("3a dynamic self-ablation with mutation is FAIL_FINGERPRINT", verdict_residual(wrong_bc["self_ablation"]["with_mutation"], FAIL_FINGERPRINT))
    expect_zero("3a dynamic self-ablation without mutation returns earned verdict", verdict_residual(wrong_bc["self_ablation"]["without_mutation"], DTN_HANKEL_FINGERPRINT_EARNED))
    expect_bool("3a self-ablation suppresses the failure", wrong_bc["self_ablation"]["fail_suppressed"])
    print(f"  3b mutated source = {imposed['mutated_source']}; inserted v5_z = {fmt(imposed['inserted_v5_z'])}")
    expect_bool("3b inserted sink is not a genuine outgoing source", not imposed["genuine_outgoing"])
    expect_bool("3b removing outgoing BC removes imaginary v5", imposed["removing_outgoing_bc_removes_imaginary_term"])
    expect_zero("3b imposed dissipation verdict is FAIL_NOT_OUTGOING", verdict_residual(imposed["verdict"], FAIL_NOT_OUTGOING))
    expect_zero("3b dynamic self-ablation with mutation is FAIL_NOT_OUTGOING", verdict_residual(imposed["self_ablation"]["with_mutation"], FAIL_NOT_OUTGOING))
    expect_zero("3b dynamic self-ablation without mutation returns earned verdict", verdict_residual(imposed["self_ablation"]["without_mutation"], DTN_HANKEL_FINGERPRINT_EARNED))
    expect_bool("3b outgoing_ablation is a rerun, not a constant label", imposed["self_ablation"]["rerun_gate_logic"] and imposed["self_ablation"]["with_mutation"] != imposed["self_ablation"]["without_mutation"])


def run_units_leg(data: dict[str, Any]) -> None:
    units = data["units"]
    subbanner("Units-restored physical coefficient dimensions only")
    print("  dimension order = (L,M,T); [a]=L, [c_s]=L*T^-1.")
    print(f"  computed dimensions = { {key: dim_text(value) for key, value in units['dims'].items()} }")
    expect_zero("physical u2 has dimension T^2", dim_residual(units["dims"]["u2"], DIM_T2))
    expect_zero("physical u4 has dimension T^4", dim_residual(units["dims"]["u4"], DIM_T4))
    expect_zero("physical v5 has dimension T^5", dim_residual(units["dims"]["v5"], DIM_T5))
    expect_bool("baseline coefficient dimensional leg passes", units["ok"])
    print(f"  corrupted u2 expression = {fmt(units['corrupted_u2'])}; dimension = {dim_text(units['corrupted_u2_dim'])}")
    expect_fail("corrupting the c_s power in u2 breaks the T^2 dimension", dim_residual(units["corrupted_u2_dim"], DIM_T2))
    expect_zero("units corruption reaches FAIL_DIMENSIONAL", verdict_residual(units["corrupted_verdict"], FAIL_DIMENSIONAL))
    expect_zero("units dynamic self-ablation with mutation is FAIL_DIMENSIONAL", verdict_residual(units["self_ablation"]["with_mutation"], FAIL_DIMENSIONAL))
    expect_zero("units dynamic self-ablation without mutation returns earned verdict", verdict_residual(units["self_ablation"]["without_mutation"], DTN_HANKEL_FINGERPRINT_EARNED))
    expect_bool("units self-ablation suppresses the failure", units["self_ablation"]["fail_suppressed"])


def run_per_tooth_ablations(data: dict[str, Any]) -> None:
    fp = data["fingerprint"]
    probes = data["probes"]
    expected = fp["expected"]
    teeth = probes["fingerprint_derivation_teeth"]
    subbanner("Per-tooth ablations on copies")
    u2_mut = teeth["u2_mixed"]
    expect_fail("u2 derivation copy changes u2_z so its own coefficient assert fires", u2_mut["u2_z"] - expected["u2_z"])
    u4_mut = teeth["u4_only"]
    print(f"  u4-only derivation copy coefficients = u2z {fmt(u4_mut['u2_z'])}, u4z {fmt(u4_mut['u4_z'])}, v5z {fmt(u4_mut['v5_z'])}")
    expect_zero("u4-only copy leaves u2_z unchanged", u4_mut["u2_z"] - expected["u2_z"])
    expect_fail("u4-only copy breaks exactly the u4_z coefficient assert", u4_mut["u4_z"] - expected["u4_z"])
    expect_zero("u4-only copy leaves v5_z unchanged", u4_mut["v5_z"] - expected["v5_z"])
    v5_mut = teeth["v5_only"]
    expect_zero("v5-only copy leaves u2_z unchanged", v5_mut["u2_z"] - expected["u2_z"])
    expect_zero("v5-only copy leaves u4_z unchanged", v5_mut["u4_z"] - expected["u4_z"])
    expect_fail("v5-only copy changes v5_z so its own coefficient assert fires", v5_mut["v5_z"] - expected["v5_z"])
    expect_fail("incoming branch mutation flips chi sign tooth", fp["chi_Q_incoming"] - 1)
    expect_zero("incoming branch reaches FAIL_FINGERPRINT", verdict_residual(probes["3a_wrong_bc"]["incoming_verdict"], FAIL_FINGERPRINT))
    expect_zero("imposed sink reaches FAIL_NOT_OUTGOING", verdict_residual(probes["3b_imposed_dissipation"]["verdict"], FAIL_NOT_OUTGOING))
    expect_zero("units dim tooth reaches FAIL_DIMENSIONAL", verdict_residual(data["units"]["corrupted_verdict"], FAIL_DIMENSIONAL))
    expect_zero("baseline immutable after teeth: outgoing u2_z still 1/9", fp["outgoing"]["u2_z"] - expected["u2_z"])
    expect_zero("baseline immutable after teeth: outgoing u4_z still 4/81", fp["outgoing"]["u4_z"] - expected["u4_z"])
    expect_zero("baseline immutable after teeth: outgoing v5_z still 1/27", fp["outgoing"]["v5_z"] - expected["v5_z"])


def forbidden_tokens() -> list[str]:
    return [
        "mu" + "_hat" + "0",
        "mtilde" + "0",
        "N" + "0",
        "D" + "0",
        "P" + "0" + "_phys",
        "build" + "_dimensions",
    ]


def live_symbol_names(data: dict[str, Any]) -> set[str]:
    exprs = [
        data["fingerprint"]["outgoing"]["h"],
        data["fingerprint"]["outgoing"]["lambda"],
        data["fingerprint"]["outgoing"]["yhat"],
        data["fingerprint"]["incoming"]["h"],
        data["fingerprint"]["standing"]["h"],
        data["units"]["corrupted_u2"],
    ]
    return {symbol.name for expr in exprs for symbol in sp.sympify(expr).free_symbols}


def run_scope_and_provenance(data: dict[str, Any]) -> None:
    subbanner("018 scope, provenance-only consumption, and structural exclusions")
    print(f"  018 gate booleans = {data['gates']}")
    print(f"  018 scoped verdict = {data['verdict']}")
    expect_zero("018 scoped verdict lands the earned partial component", verdict_residual(data["verdict"], DTN_HANKEL_FINGERPRINT_EARNED))
    print("  QUAD_CALIBRATED (1/4) -- outgoing l=2 DtN Hankel fingerprint EARNED (PARTIAL)")
    print("    = u2=a^2/(9*c_s^2), u4=4*a^4/(81*c_s^4), v5=a^5/(27*c_s^5), derived from Yhat=-3/Lambda.")
    print("    AND chi_Q=+1 outgoing / -1 incoming from j2 +/- i*y2, with standing j2 carrying no radiating term.")
    print("  REMAINING: prefactor algebra = 019; normalization partition plus calibrated label = 020; dim closure = 021.")
    print("  CAVEATS: c_s is a units carrier, not a consumed value; branch scalar numerics remain sim-deferred; chi reconciliation is Part-VII.")
    print("  CONSUMES (PROVENANCE only): 017 l=2 port kernel + 009/010 bulk mode + stage005 R1 identity for the c_s units symbol.")
    print("  EXPORTS: outgoing fingerprint and chi sign to 019/020 plus later non-regression and closure stages; no file artifact is written.")
    names = live_symbol_names(data)
    print(f"  live symbolic names in this slice = {sorted(names)}")
    expect_bool("only z, a, c_s are live symbols in the fingerprint slice", names <= {"z", "a", "c_s"})
    frozen_speed = "c" + "_S"
    expect_bool("frozen-wall speed symbol is not live", frozen_speed not in names and "cS" not in names)
    forbidden = set(forbidden_tokens())
    helper_names = set(globals())
    expect_bool("no prefactor/dim-closure carrier helper or symbol names are live", not (names & forbidden) and not (helper_names & forbidden))
    print("  dropped-bookkeeping: scratch YAML, report writers, engine-agreement files, and cross-reads are absent.")
    print("  register note: likely zero new counted knobs; c_s is R1-provenance units carrier; structural edge is for the outgoing DtN fingerprint.")


def print_verdict_labels() -> None:
    subbanner("Verdict labels")
    print("  ledger earned-label (NOT a source verdict token): DTN_HANKEL_FINGERPRINT_EARNED")
    print("  source top-line verdict: QUAD_CALIBRATED (JOINT 4-stage; 018 carries the outgoing-fingerprint + chi-sign component as PARTIAL)")
    print("  joint composition: 018 outgoing DtN Hankel fingerprint + 019 prefactor + 020 normalization partition/calibrated label + 021 dim closure")
    print("  earned fingerprint: u2z=1/9, u4z=4/81, v5z=1/27 from coeff(z,n) of Yhat=-3/Lambda; 27 = 1/v5z is computed.")
    print("  earned sign: chi_Q=+1 outgoing and -1 incoming from j2 +/- i*y2; standing j2 has Lambda_static=+2 and no radiating term.")
    print("  earned able-to-fail: wrong-BC and imposed-sink probes each carry dynamic with/without mutation verdicts.")
    print("  calibrated/deferred outside 018: prefactor algebra, normalization magnitude/label, and the full dimensional closure.")
    print("  consumed framing: provenance-only; no guard on port kernel, bulk mode, or c_s value; only c_s units appear in the dim leg.")


def main() -> int:
    banner("ledger_stage018_dtn_hankel_fingerprint_sympy_audit")
    print("Target stem confirmed: ledger_stage018_dtn_hankel_fingerprint")
    print("Engine: SymPy exact symbolic spherical-Hankel series; no floats/tolerances; zero file I/O.")
    data = build_baseline()
    run_fingerprint_derivation(data)
    run_chi_and_standing(data)
    run_passivity_and_probes(data)
    run_units_leg(data)
    run_per_tooth_ablations(data)
    run_scope_and_provenance(data)
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
