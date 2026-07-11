#!/usr/bin/env python3
"""Ledger stage027 SymPy audit: density-port checks, closure, and joint.

Standalone, print-only, no arguments, and zero file I/O.  Stage024's
canonical factored N0_den and the atomic stage025/026 certificates are cited,
not reconstructed.  This completing leg owns the DENSITY_PORT_HOSTED joint.
"""

from __future__ import annotations

import inspect
from dataclasses import dataclass, replace
from functools import lru_cache
from typing import Any, Callable

import sympy as sp


PASS_COUNT = 0
FAIL_COUNT = 0

JOINT_VERDICT = "DENSITY_PORT_HOSTED"
INCONCLUSIVE = "PORT_INCONCLUSIVE_SIM_DEFERRED"
FAIL_ORIGIN = "FAIL_NOT_DENSITY_DERIVED"
FAIL_VANISHES = "FAIL_PORT_VANISHES"


class AuditFailure(AssertionError):
    """A named stage027 audit assertion failed."""


class RigAssertion(AssertionError):
    """An expected control mutation reached its routed named assertion."""

    def __init__(self, name: str):
        super().__init__(name)
        self.name = name


class DimError(ValueError):
    """A sourced dimension is missing or an additive expression is mixed."""


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
    if not isinstance(expr, sp.Basic):
        return expr
    return sp.factor(sp.cancel(sp.together(sp.simplify(expr))))


def fmt(expr: Any) -> str:
    return expr if isinstance(expr, str) else sp.sstr(compact(expr))


def fmt_names(items: set[Any] | frozenset[Any]) -> str:
    return "{" + ", ".join(sorted(str(item) for item in items)) + "}"


def _record_pass(name: str) -> None:
    global PASS_COUNT
    PASS_COUNT += 1
    print(f"PASS  {name}")


def _record_fail(name: str, detail: str) -> None:
    global FAIL_COUNT
    FAIL_COUNT += 1
    print(f"FAIL  {name}: {detail}")


def expect_zero(name: str, residual: sp.Expr | int) -> None:
    clean = compact(sp.sympify(residual))
    if clean == 0:
        _record_pass(name)
        return
    _record_fail(name, f"residual = {fmt(clean)}")
    raise AuditFailure(name)


def expect_bool(name: str, condition: bool) -> None:
    expect_zero(name, 0 if bool(condition) else 1)


def routed_assert(name: str, condition: bool) -> None:
    if not bool(condition):
        raise RigAssertion(name)


def probe_assertion(
    assertion: Callable[[], None], expected_name: str
) -> tuple[bool, str]:
    try:
        assertion()
    except RigAssertion as exc:
        return exc.name == expected_name, exc.name
    return False, "NO_ASSERT_FIRED"


def assertion_passes(assertion: Callable[[], None]) -> bool:
    try:
        assertion()
    except RigAssertion:
        return False
    return True


def exercise_rig(
    label: str,
    assertion_name: str,
    rig_assertion: Callable[[], None],
    neutral_assertion: Callable[[], None],
    *,
    outcome: str,
) -> str:
    caught, fired_name = probe_assertion(rig_assertion, assertion_name)
    neutral_pass = assertion_passes(neutral_assertion)
    expect_bool(
        f"META {label} routed assertion fires and neutering stops it",
        caught and neutral_pass,
    )
    text = f"{outcome} at {fired_name}; neutralized=PASS"
    print(f"RIG {label}: {text}")
    return text


# Stage024's exact export symbols.  q2/Phi2 are coordinate metadata only.
a, c_s, rho_eff = sp.symbols("a c_s rho_eff", positive=True, real=True)
I25, I_wrong2, Xi_Q, Xi_deferred = sp.symbols(
    "I25 I_wrong2 Xi_Q Xi_deferred", nonzero=True, real=True
)
eta_q, eta_phi = sp.symbols("eta_q eta_phi", nonzero=True, real=True)
varpi_q2, varpi_Phi2 = sp.symbols(
    "varpi_q2 varpi_Phi2", positive=True, real=True
)
lambda_c = sp.Symbol("lambda_c", real=True)
q2, Phi2 = sp.symbols("q2 Phi2", real=True)
c, G, D0, z = sp.symbols("c G D0 z", positive=True, real=True)
foreign_subject = sp.Symbol("foreign_subject", nonzero=True, real=True)
delta4, delta_gamma = sp.symbols("delta4 deltaGamma", nonzero=True)

HOST_CONTRACT = {
    a,
    c_s,
    rho_eff,
    I25,
    Xi_Q,
    eta_q,
    eta_phi,
    varpi_q2,
    varpi_Phi2,
    lambda_c,
}

# Atomic stage025 export over the baseline cited N0_den.  It is not live-
# retainted on controls (especially not after zero_coupling specializes to 0).
BASELINE_TAGS = frozenset(
    {"continuity_interface", "pathA_29_bulk", "pathA_32_wall"}
)


@dataclass(frozen=True)
class Certificates:
    tags: frozenset[str] = BASELINE_TAGS
    vector_host_symbols: frozenset[str] = frozenset()
    source_map_complete: bool = True
    vector_free: bool = True
    moment_valid: bool = True
    lineage_valid: bool = True
    validated_I25: sp.Symbol = I25


@dataclass(frozen=True)
class Config:
    name: str
    coupling_zero: bool = False
    corrupt_dimension: bool = False
    corrupt_g_dimension: bool = False
    incoming_sign: bool = False
    coupling_a_power: sp.Rational = sp.Rational(-7, 2)
    deferred_uncertified: bool = False
    proven_deferred: bool = False
    second_subject_arm: bool = True
    kbar4_shift: sp.Expr = sp.Integer(0)
    gamma5_shift: sp.Expr = sp.Integer(0)


Dim = tuple[sp.Rational, sp.Rational, sp.Rational]
ZERO_DIM: Dim = (sp.Rational(0), sp.Rational(0), sp.Rational(0))
BASE_DIMS: dict[sp.Symbol, Dim] = {
    a: (1, 0, 0),
    c_s: (1, 0, -1),
    c: (1, 0, -1),
    G: (3, -1, -2),
    D0: (-1, 1, -2),
    rho_eff: (-3, 1, 0),
    I25: (sp.Rational(5, 2), 0, 0),
    I_wrong2: (2, 0, 0),
    Xi_Q: ZERO_DIM,
    Xi_deferred: ZERO_DIM,
    eta_q: ZERO_DIM,
    eta_phi: ZERO_DIM,
    varpi_q2: (0, 0, -2),
    varpi_Phi2: (0, 0, -2),
    lambda_c: (0, 0, -2),
    foreign_subject: (-3, 1, 0),
}
BASE_A_POWERS: dict[sp.Symbol, sp.Rational | None] = {
    c_s: 0,
    c: 0,
    G: 0,
    D0: 0,
    rho_eff: 0,
    I25: 0,
    I_wrong2: 0,
    Xi_Q: 0,
    Xi_deferred: 0,
    eta_q: 0,
    eta_phi: 0,
    varpi_q2: -2,
    varpi_Phi2: -2,
    lambda_c: -2,
    foreign_subject: 0,
}


def cited_N0_den() -> sp.Expr:
    """Cite stage024's canonical factored expression; never invert a matrix."""
    return (
        I25**2
        * Xi_Q**2
        * c_s**4
        * rho_eff
        * (eta_phi * varpi_q2 + eta_q * lambda_c) ** 2
        / (a**7 * (lambda_c**2 - varpi_Phi2 * varpi_q2) ** 2)
    )


def specialize_cited_port(config: Config) -> tuple[sp.Expr, sp.Symbol]:
    """Specialize the cited closed form for a control; do not re-solve 024."""
    xi = Xi_deferred if config.deferred_uncertified or config.proven_deferred else Xi_Q
    moment = I25 if config.coupling_a_power == sp.Rational(-7, 2) else I_wrong2
    expr = cited_N0_den().subs({I25: moment, Xi_Q: xi})
    # The source's g_base mutation -7/2 -> -3 changes its square by one
    # factor of a.  Apply that specialization to the cited export itself.
    if moment == I_wrong2:
        expr = expr * a
    return (sp.Integer(0) if config.coupling_zero else compact(expr), moment)


def dim_add(left: Dim, right: Dim) -> Dim:
    return tuple(sp.Rational(x + y) for x, y in zip(left, right))  # type: ignore[return-value]


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
        if any(value != term_dims[0] for value in term_dims):
            raise DimError(f"dimension mismatch in sum {expr}: {term_dims}")
        return term_dims[0]
    raise DimError(f"unsupported dimension expression {expr}")


def scale_power(
    expr: sp.Expr, powers: dict[sp.Symbol, sp.Rational | None]
) -> sp.Rational | None:
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
            value = scale_power(arg, powers)
            if value is None:
                return None
            total += value
        return sp.Rational(total)
    if expr.is_Pow:
        base, exponent = expr.args
        value = scale_power(base, powers)
        return None if value is None or not exponent.is_number else sp.Rational(exponent) * value
    if expr.is_Add:
        values = [scale_power(arg, powers) for arg in expr.args if arg != 0]
        if not values or any(value is None for value in values):
            return sp.Rational(0) if not values else None
        return values[0] if all(value == values[0] for value in values) else None
    return None


@lru_cache(maxsize=2)
def dtn_sign(kind: str) -> dict[str, Any]:
    """Compute chi's sign from j2 +/- i*y2, never from a typed chi value."""
    j2 = (3 / z**3 - 1 / z) * sp.sin(z) - 3 * sp.cos(z) / z**2
    y2 = (1 / z - 3 / z**3) * sp.cos(z) - 3 * sp.sin(z) / z**2
    hankel = j2 + sp.I * y2 if kind == "outgoing" else j2 - sp.I * y2
    response = compact(-3 / compact(z * sp.diff(hankel, z) / hankel))
    series = sp.expand(sp.series(response, z, 0, 7).removeO())
    coeff_over_i = compact(series.coeff(z, 5) / sp.I)
    return {
        "kind": kind,
        "series": series,
        "coeff_over_i": coeff_over_i,
        "chi_Q": 1 if coeff_over_i == sp.Rational(1, 27) else -1,
        "ok": coeff_over_i == sp.Rational(1, 27),
    }


def closure_overlay(
    n0_den: sp.Expr,
    *,
    kbar4_shift: sp.Expr = sp.Integer(0),
    gamma5_shift: sp.Expr = sp.Integer(0),
) -> dict[str, Any]:
    p0_physical = compact((c_s / a) ** 2 * n0_den / D0)
    kbar0 = 54 * G * c_s**5 / (5 * a**5 * c**5)
    kbar2 = 6 * G * c_s**3 / (5 * a**3 * c**5)
    kbar4 = 8 * G * c_s / (15 * a * c**5) + kbar4_shift
    gamma5 = compact(kbar0 * a**5 / (27 * c_s**5)) + gamma5_shift
    kbar4_residual = compact(kbar4 - 4 * kbar2**2 / kbar0)
    gamma5_residual = compact(gamma5 - 2 * G / (5 * c**5))
    return {
        "P0_physical": p0_physical,
        "Kbar0": kbar0,
        "Kbar2": kbar2,
        "Kbar4": kbar4,
        "Gamma5": gamma5,
        "Kbar4_residual": kbar4_residual,
        "Gamma5_residual": gamma5_residual,
        "ok": kbar4_residual == 0 and gamma5_residual == 0,
    }


def evaluate_port(config: Config, certificates: Certificates) -> dict[str, Any]:
    expr, fixture_moment = specialize_cited_port(config)
    port_moment = (
        certificates.validated_I25
        if config.coupling_a_power == sp.Rational(-7, 2)
        else I_wrong2
    )
    dims = dict(BASE_DIMS)
    if config.corrupt_dimension:
        dims[I25] = dim_add(dims[I25], (1, 0, 0))
    if config.corrupt_g_dimension:
        dims[G] = dim_add(dims[G], (1, 0, 0))
    try:
        expr_dim = dim_of(expr, dims)
        dim_ok = expr_dim == (-1, 1, 0)
        dim_error = None
    except DimError as exc:
        expr_dim = None
        dim_ok = False
        dim_error = str(exc)

    powers = dict(BASE_A_POWERS)
    if config.deferred_uncertified:
        powers[Xi_deferred] = None
    if config.proven_deferred:
        powers[Xi_deferred] = sp.Rational(0)
    n0_power = scale_power(expr, powers)
    p0_power = None if n0_power is None else sp.Rational(n0_power - 2 - powers[D0])
    scaling_wrong = p0_power is not None and p0_power != -5
    scaling_ok = p0_power == -5
    scaling_undecidable = p0_power is None

    sign = dtn_sign("incoming" if config.incoming_sign else "outgoing")
    nonzero_ok = compact(expr) != 0
    closure = closure_overlay(
        expr,
        kbar4_shift=config.kbar4_shift,
        gamma5_shift=config.gamma5_shift,
    )

    membership_arm = port_moment in expr.free_symbols
    vanished_arm = bool(
        config.second_subject_arm and compact(expr) == 0 and config.coupling_zero
    )
    subject_binding = membership_arm or vanished_arm
    origin_ok = bool(
        "continuity_interface" in certificates.tags
        and "vector_port" not in certificates.tags
        and not certificates.vector_host_symbols
        and certificates.source_map_complete
        and certificates.vector_free
        and certificates.lineage_valid
        and certificates.moment_valid
        and subject_binding
    )

    if not origin_ok:
        verdict = FAIL_ORIGIN
    elif not nonzero_ok:
        verdict = FAIL_VANISHES
    elif not dim_ok or scaling_wrong or not sign["ok"]:
        bad: list[str] = []
        if not dim_ok:
            bad.append("dimensional")
        if scaling_wrong:
            bad.append("scaling")
        if not sign["ok"]:
            bad.append("sign")
        verdict = "FAIL_PORT_MALFORMED(" + ",".join(bad) + ")"
    elif (
        origin_ok
        and nonzero_ok
        and dim_ok
        and scaling_ok
        and sign["ok"]
        and closure["ok"]
    ):
        verdict = JOINT_VERDICT
    else:
        verdict = INCONCLUSIVE

    return {
        "config": config,
        "expr": expr,
        "fixture_moment": fixture_moment,
        "port_moment": port_moment,
        "expr_dim": expr_dim,
        "dim_error": dim_error,
        "dim_ok": dim_ok,
        "n0_power": n0_power,
        "p0_power": p0_power,
        "scaling_wrong": scaling_wrong,
        "scaling_ok": scaling_ok,
        "scaling_undecidable": scaling_undecidable,
        "sign": sign,
        "nonzero_ok": nonzero_ok,
        "closure": closure,
        "membership_arm": membership_arm,
        "vanished_arm": vanished_arm,
        "subject_binding": subject_binding,
        "origin_ok": origin_ok,
        "verdict": verdict,
        "retirement_recorded": verdict == JOINT_VERDICT,
    }


def host_contract_ok(expr: sp.Expr) -> bool:
    return set(expr.free_symbols) == HOST_CONTRACT


def arity_scan(call_arity: int) -> bool:
    positional = [
        parameter
        for parameter in inspect.signature(evaluate_port).parameters.values()
        if parameter.kind
        in (parameter.POSITIONAL_ONLY, parameter.POSITIONAL_OR_KEYWORD)
    ]
    return call_arity == len(positional)


AUTHORED_LEAK = sp.Function("evaluate_port")


def leakage_scan(objects: list[Any]) -> bool:
    return not any(
        isinstance(obj, sp.Basic) and obj.has(AUTHORED_LEAK) for obj in objects
    )


def run_audit() -> dict[str, Any]:
    cert = Certificates()
    baseline_expr = compact(cited_N0_den())
    baseline = evaluate_port(Config("baseline"), cert)

    subbanner("I. consumed sibling exports and subject-integrity oracle")
    expect_bool(
        "SUBJECT-INTEGRITY cited N0_den free_symbols equals exact 10-symbol HOST_CONTRACT",
        host_contract_ok(baseline_expr),
    )
    expect_bool(
        "ASSEMBLY consumed 025 exact baseline taint set and atomic vector facts",
        cert.tags == BASELINE_TAGS
        and "continuity_interface" in cert.tags
        and "vector_port" not in cert.tags
        and not cert.vector_host_symbols
        and cert.source_map_complete
        and cert.vector_free,
    )
    expect_bool(
        "ASSEMBLY consumed 026 atomic moment_valid, validated_I25, and lineage_valid facts",
        cert.moment_valid and cert.validated_I25 == I25 and cert.lineage_valid,
    )

    subbanner("Six port checks and independent closure residuals")
    expect_bool("DIM [N0_den] == (-1,1,0) = L^-1 M", baseline["dim_ok"])
    expect_bool("SCALE P0_physical a-power == -5", baseline["p0_power"] == -5)
    expect_bool(
        "SIGN outgoing j2+i*y2 gives coeff/i == 1/27 and chi_Q=+1",
        baseline["sign"]["coeff_over_i"] == sp.Rational(1, 27)
        and baseline["sign"]["chi_Q"] == 1,
    )
    expect_bool("NONZERO cited density port is nonzero", baseline["nonzero_ok"])
    expect_zero(
        "CLOSURE Kbar4-4*Kbar2^2/Kbar0 standalone residual",
        baseline["closure"]["Kbar4_residual"],
    )
    expect_zero(
        "CLOSURE Gamma5-2*G/(5*c^5) standalone residual",
        baseline["closure"]["Gamma5_residual"],
    )
    expect_bool("POSITIVE assembled joint lands DENSITY_PORT_HOSTED", baseline["verdict"] == JOINT_VERDICT)

    subbanner("Per-control verdicts and routed per-tooth ablations")
    outcomes: dict[str, str] = {}
    controls = {
        "dimensional": evaluate_port(Config("dimensional", corrupt_dimension=True), cert),
        "corrupt_G_scope": evaluate_port(Config("corrupt_G_scope", corrupt_g_dimension=True), cert),
        "scaling": evaluate_port(Config("scaling", coupling_a_power=sp.Rational(-3)), cert),
        "sign": evaluate_port(Config("sign", incoming_sign=True), cert),
        "zero_coupling": evaluate_port(Config("zero_coupling", coupling_zero=True), cert),
        "zero_no_or": evaluate_port(
            Config("zero_no_or", coupling_zero=True, second_subject_arm=False), cert
        ),
        "deferred_scalar": evaluate_port(
            Config("deferred_scalar", deferred_uncertified=True), cert
        ),
        "deferred_scalar_proven_converse": evaluate_port(
            Config("deferred_scalar_proven_converse", proven_deferred=True), cert
        ),
        "closure_kbar4": evaluate_port(
            Config("closure_kbar4", kbar4_shift=delta4), cert
        ),
        "closure_gamma5": evaluate_port(
            Config("closure_gamma5", gamma5_shift=delta_gamma), cert
        ),
    }

    expected_controls = {
        "dimensional": "FAIL_PORT_MALFORMED(dimensional)",
        "corrupt_G_scope": JOINT_VERDICT,
        "scaling": "FAIL_PORT_MALFORMED(scaling)",
        "sign": "FAIL_PORT_MALFORMED(sign)",
        "zero_coupling": FAIL_VANISHES,
        "zero_no_or": FAIL_ORIGIN,
        "deferred_scalar": INCONCLUSIVE,
        "deferred_scalar_proven_converse": JOINT_VERDICT,
        "closure_kbar4": INCONCLUSIVE,
        "closure_gamma5": INCONCLUSIVE,
    }
    for label, expected in expected_controls.items():
        expect_bool(
            f"CONTROL {label} computed verdict == {expected}",
            controls[label]["verdict"] == expected,
        )
        print(f"CONTROL {label}: {controls[label]['verdict']}")

    expect_bool(
        "SCALE single-tag isolation keeps dim_ok True and only scaling_wrong fires",
        controls["scaling"]["dim_ok"]
        and controls["scaling"]["scaling_wrong"]
        and controls["scaling"]["verdict"] == "FAIL_PORT_MALFORMED(scaling)",
    )
    expect_bool(
        "NONZERO OR-arm meta flips zero_coupling FAIL_PORT_VANISHES to FAIL_NOT_DENSITY_DERIVED",
        controls["zero_coupling"]["verdict"] == FAIL_VANISHES
        and controls["zero_no_or"]["verdict"] == FAIL_ORIGIN,
    )
    expect_bool(
        "CLOSURE isolated Kbar4 mutation leaves Gamma5 residual zero",
        controls["closure_kbar4"]["closure"]["Kbar4_residual"] != 0
        and controls["closure_kbar4"]["closure"]["Gamma5_residual"] == 0,
    )
    expect_bool(
        "CLOSURE isolated Gamma5 mutation leaves Kbar4 residual zero",
        controls["closure_gamma5"]["closure"]["Gamma5_residual"] != 0
        and controls["closure_gamma5"]["closure"]["Kbar4_residual"] == 0,
    )

    rig_specs = [
        (
            "DIM dimensional",
            "DIM [N0_den]=L^-1 M named assert",
            lambda: routed_assert("DIM [N0_den]=L^-1 M named assert", controls["dimensional"]["dim_ok"]),
            lambda: routed_assert("DIM [N0_den]=L^-1 M named assert", baseline["dim_ok"]),
            controls["dimensional"]["verdict"],
        ),
        (
            "SCALE scaling-single-tag",
            "SCALE p0_power=-5 named assert",
            lambda: routed_assert("SCALE p0_power=-5 named assert", controls["scaling"]["scaling_ok"]),
            lambda: routed_assert("SCALE p0_power=-5 named assert", baseline["scaling_ok"]),
            controls["scaling"]["verdict"],
        ),
        (
            "SIGN incoming",
            "SIGN outgoing coeff/i=1/27 named assert",
            lambda: routed_assert("SIGN outgoing coeff/i=1/27 named assert", controls["sign"]["sign"]["ok"]),
            lambda: routed_assert("SIGN outgoing coeff/i=1/27 named assert", baseline["sign"]["ok"]),
            controls["sign"]["verdict"],
        ),
        (
            "NONZERO zero_coupling",
            "NONZERO N0_den!=0 named assert",
            lambda: routed_assert("NONZERO N0_den!=0 named assert", controls["zero_coupling"]["nonzero_ok"]),
            lambda: routed_assert("NONZERO N0_den!=0 named assert", baseline["nonzero_ok"]),
            controls["zero_coupling"]["verdict"],
        ),
        (
            "DEFERRED uncertified",
            "DEFERRED scaling-decidable named assert",
            lambda: routed_assert("DEFERRED scaling-decidable named assert", not controls["deferred_scalar"]["scaling_undecidable"]),
            lambda: routed_assert("DEFERRED scaling-decidable named assert", not controls["deferred_scalar_proven_converse"]["scaling_undecidable"]),
            controls["deferred_scalar"]["verdict"],
        ),
        (
            "CLOSURE Kbar4 isolated",
            "CLOSURE Kbar4 standalone named assert",
            lambda: routed_assert("CLOSURE Kbar4 standalone named assert", controls["closure_kbar4"]["closure"]["Kbar4_residual"] == 0),
            lambda: routed_assert("CLOSURE Kbar4 standalone named assert", baseline["closure"]["Kbar4_residual"] == 0),
            controls["closure_kbar4"]["verdict"],
        ),
        (
            "CLOSURE Gamma5 isolated",
            "CLOSURE Gamma5 standalone named assert",
            lambda: routed_assert("CLOSURE Gamma5 standalone named assert", controls["closure_gamma5"]["closure"]["Gamma5_residual"] == 0),
            lambda: routed_assert("CLOSURE Gamma5 standalone named assert", baseline["closure"]["Gamma5_residual"] == 0),
            controls["closure_gamma5"]["verdict"],
        ),
    ]
    for label, assertion, rig, neutral, outcome in rig_specs:
        outcomes[label] = exercise_rig(label, assertion, rig, neutral, outcome=outcome)

    # Each consumed atomic sibling fact is independently load-bearing.
    atomic_mutations = {
        "ASSEMBLY continuity_interface": replace(cert, tags=frozenset(set(cert.tags) - {"continuity_interface"})),
        "ASSEMBLY vector_port": replace(cert, tags=frozenset(set(cert.tags) | {"vector_port"})),
        "ASSEMBLY vector_host_symbols": replace(cert, vector_host_symbols=frozenset({"A_w"})),
        "ASSEMBLY source_map_complete": replace(cert, source_map_complete=False),
        "ASSEMBLY vector_free": replace(cert, vector_free=False),
        "ASSEMBLY moment_valid": replace(cert, moment_valid=False),
        "ASSEMBLY validated_I25": replace(cert, validated_I25=foreign_subject),
        "ASSEMBLY lineage_valid": replace(cert, lineage_valid=False),
    }
    for label, bad_cert in atomic_mutations.items():
        mutated = evaluate_port(Config(label), bad_cert)
        assert_name = label + " origin_ok named assert"
        expect_bool(label + " computed verdict == FAIL_NOT_DENSITY_DERIVED", mutated["verdict"] == FAIL_ORIGIN)
        outcomes[label] = exercise_rig(
            label,
            assert_name,
            lambda result=mutated, name=assert_name: routed_assert(name, result["origin_ok"]),
            lambda name=assert_name: routed_assert(name, baseline["origin_ok"]),
            outcome=mutated["verdict"],
        )

    outcomes["ASSEMBLY OR-arm"] = exercise_rig(
        "ASSEMBLY OR-arm",
        "ASSEMBLY zero-coupling subject-binding OR-arm named assert",
        lambda: routed_assert(
            "ASSEMBLY zero-coupling subject-binding OR-arm named assert",
            controls["zero_no_or"]["subject_binding"],
        ),
        lambda: routed_assert(
            "ASSEMBLY zero-coupling subject-binding OR-arm named assert",
            controls["zero_coupling"]["subject_binding"],
        ),
        outcome=f"{FAIL_VANISHES}->{FAIL_ORIGIN}",
    )

    corrupt_subject = compact(baseline_expr.subs(rho_eff, foreign_subject))
    expect_bool(
        "SUBJECT-INTEGRITY foreign replacement preserves dimension and a-power",
        dim_of(corrupt_subject, BASE_DIMS) == (-1, 1, 0)
        and scale_power(corrupt_subject, BASE_A_POWERS)
        == scale_power(baseline_expr, BASE_A_POWERS),
    )
    outcomes["SUBJECT-INTEGRITY"] = exercise_rig(
        "SUBJECT-INTEGRITY",
        "SUBJECT-INTEGRITY exact host-contract named assert",
        lambda: routed_assert(
            "SUBJECT-INTEGRITY exact host-contract named assert",
            host_contract_ok(corrupt_subject),
        ),
        lambda: routed_assert(
            "SUBJECT-INTEGRITY exact host-contract named assert",
            host_contract_ok(baseline_expr),
        ),
        outcome="CONTRACT_REJECTED_ONLY",
    )

    positive_flips = {
        name: result["verdict"] != JOINT_VERDICT
        for name, result in controls.items()
        if name in {"dimensional", "scaling", "sign", "zero_coupling", "deferred_scalar", "closure_kbar4"}
    }
    expect_bool(
        "POSITIVE baseline hosts and corruption of any of six checks flips the joint",
        baseline["verdict"] == JOINT_VERDICT and all(positive_flips.values()),
    )
    outcomes["POSITIVE"] = "DENSITY_PORT_HOSTED; any-check-corruption=FLIP"
    print(f"RIG POSITIVE: {outcomes['POSITIVE']}")

    outcomes["RETIREMENT-CONDITIONAL"] = exercise_rig(
        "RETIREMENT-CONDITIONAL",
        "RETIREMENT recorded only with hosted joint named assert",
        lambda: routed_assert(
            "RETIREMENT recorded only with hosted joint named assert",
            controls["sign"]["retirement_recorded"],
        ),
        lambda: routed_assert(
            "RETIREMENT recorded only with hosted joint named assert",
            baseline["retirement_recorded"],
        ),
        outcome="NOT_RECORDED_ON_FAIL",
    )

    outcomes["ARITY"] = exercise_rig(
        "ARITY scanner",
        "ARITY definition/call scanner named assert",
        lambda: routed_assert("ARITY definition/call scanner named assert", arity_scan(1)),
        lambda: routed_assert("ARITY definition/call scanner named assert", arity_scan(2)),
        outcome="SCANNER_CAUGHT",
    )
    actual_transcript = [baseline_expr, baseline["closure"]["Kbar4_residual"], baseline["verdict"]]
    leaked_transcript = actual_transcript + [AUTHORED_LEAK(baseline_expr)]
    outcomes["LEAKAGE"] = exercise_rig(
        "LEAKAGE scanner",
        "LEAKAGE unevaluated-helper scanner named assert",
        lambda: routed_assert("LEAKAGE unevaluated-helper scanner named assert", leakage_scan(leaked_transcript)),
        lambda: routed_assert("LEAKAGE unevaluated-helper scanner named assert", leakage_scan(actual_transcript)),
        outcome="SCANNER_CAUGHT",
    )
    expect_bool("LEAKAGE actual transcript is helper-free", leakage_scan(actual_transcript))

    return {
        "baseline_expr": baseline_expr,
        "cert": cert,
        "baseline": baseline,
        "controls": controls,
        "outcomes": outcomes,
    }


def emit(data: dict[str, Any]) -> None:
    baseline = data["baseline"]
    cert = data["cert"]
    closure = baseline["closure"]
    subbanner("Stage027 checks + closure transcript")
    print("consumes stage024 N0_den (canonical factored): " + fmt(data["baseline_expr"]))
    print("N0_den free_symbols: " + fmt_names(set(data["baseline_expr"].free_symbols)))
    print("HOST_CONTRACT (10): " + fmt_names(HOST_CONTRACT))
    print(f"free_symbols == contract: {host_contract_ok(data['baseline_expr'])}")
    print(
        "consumed stage025 atomics: tags=" + fmt_names(set(cert.tags))
        + f"; vector_port_not_in_tags={'vector_port' not in cert.tags}"
        + f"; vector_host_symbols={fmt_names(set(cert.vector_host_symbols))}"
        + f"; source_map_complete={cert.source_map_complete}; vector_free(P2)={cert.vector_free}"
    )
    print(
        "consumed stage026 atomics: moment_valid=True; validated_I25=I25; "
        "lineage_certificate=PASS; lineage_valid=True"
    )
    print("dimension tuple convention: (L,M,T)")
    print(f"CHECK DIM: [N0_den]={baseline['expr_dim']} = L^-1 M; dim_ok={baseline['dim_ok']}")
    print("CHECK SCALE: P0_physical=(c_s/a)^2*N0_den/D0")
    print(f"CHECK SCALE: p0_power={baseline['p0_power']}; scaling_ok={baseline['scaling_ok']}")
    print(
        "CHECK SIGN: outgoing j2+i*y2 coeff/i="
        f"{fmt(baseline['sign']['coeff_over_i'])}; chi_Q=+1; sign_ok={baseline['sign']['ok']}"
    )
    print(f"CHECK NONZERO: nonzero_ok={baseline['nonzero_ok']}")
    print(
        "CHECK DEFERRED: uncertified->PORT_INCONCLUSIVE_SIM_DEFERRED; "
        "proven converse->DENSITY_PORT_HOSTED"
    )
    print("CHECK CLOSURE: Kbar4-4*Kbar2^2/Kbar0=" + fmt(closure["Kbar4_residual"]))
    print("CHECK CLOSURE: Gamma5-2*G/(5*c^5)=" + fmt(closure["Gamma5_residual"]))
    print(
        "EXPORT Kbar moments: {Kbar0=" + fmt(closure["Kbar0"])
        + ", Kbar2=" + fmt(closure["Kbar2"])
        + ", Kbar4=" + fmt(closure["Kbar4"])
        + ", Gamma5=" + fmt(closure["Gamma5"]) + "}"
    )

    subbanner("Joint-owning verdict and completion")
    print(f"LOCAL_AUDIT_VERDICT: {baseline['verdict']} (CALIBRATED)")
    print(f"JOINT_VERDICT_OWNED_BY_027: {baseline['verdict']} (CALIBRATED, not PASS)")
    print("scope: STRUCTURE hosted; magnitude SIM_DEFERRED; G=GENUINE_BLOCKED")
    print("provenance: 27=derived_in_gate; 54/5 and Gamma5=external_bridge_input")
    print(
        "RETIREMENT_RECORD: EM A_w/U,W scaffold RETIRES; pathA_43 diagnostic sliver CLOSES "
        f"(conditional={baseline['retirement_recorded']})"
    )
    print(
        "COMPLETION: pathA_43 4-way split COMPLETE = 024 derivation AND 025 vector-freedom "
        "AND 026 lineage AND 027 checks+closure"
    )


def main() -> int:
    banner("ledger_stage027_port_checks_closure_sympy_audit")
    print("Target stem confirmed: ledger_stage027_port_checks_closure")
    print("Engine: exact SymPy runtime checks over cited exports; no floats; zero file I/O.")
    data = run_audit()
    emit(data)
    return 0


if __name__ == "__main__":
    try:
        exit_code = main()
    except Exception as exc:
        if not isinstance(exc, AuditFailure):
            _record_fail("UNCAUGHT exception", repr(exc))
        banner("Tallies")
        total = PASS_COUNT + FAIL_COUNT
        print(f"TALLY sympy: {PASS_COUNT} pass + {FAIL_COUNT} fail = {total} checks")
        print("OVERALL FAIL")
        raise SystemExit(1) from exc

    banner("Tallies")
    total = PASS_COUNT + FAIL_COUNT
    print(f"TALLY sympy: {PASS_COUNT} pass + {FAIL_COUNT} fail = {total} checks")
    if exit_code == 0 and FAIL_COUNT == 0:
        print("OVERALL PASS")
        raise SystemExit(0)
    print("OVERALL FAIL")
    raise SystemExit(1)
