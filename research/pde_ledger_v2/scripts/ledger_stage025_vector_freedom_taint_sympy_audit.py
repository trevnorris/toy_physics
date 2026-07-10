#!/usr/bin/env python3
"""Ledger stage025 SymPy audit: density-port vector-freedom taint proof.

Standalone, print-only, no arguments, and no file I/O.  Stage024's factored
N0_den is cited directly.  The decisive predicates are computed provenance
taint and source-map/host completeness; expression ablation is only a
redundant printed witness.  Continuity lineage and port closure belong to
stages026 and 027 respectively.
"""

from __future__ import annotations

import inspect
from dataclasses import dataclass
from typing import Any, Callable

import sympy as sp


PASS_COUNT = 0
FAIL_COUNT = 0

LOCAL_VERDICT = "DENSITY_PORT_VECTOR_FREE"
JOINT_PARTIAL = (
    "DENSITY_PORT_HOSTED (2/4, VECTOR-FREEDOM — N0_den is computationally "
    "vector-free; 024 = derivation, 026 = continuity lineage, "
    "027 = port checks + closure)"
)
BASELINE_TAGS = {
    "continuity_interface",
    "pathA_29_bulk",
    "pathA_32_wall",
}
TAGGED_CARRIER_TAG = "pathA_34_dimensionless_free_carrier"


class AuditFailure(AssertionError):
    """A named stage025 audit assertion failed."""


class RigAssertion(AssertionError):
    """An expected rig fired at its routed named assertion."""

    def __init__(self, name: str):
        super().__init__(name)
        self.name = name


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
    if isinstance(expr, str):
        return expr
    return sp.sstr(compact(expr))


def fmt_names(items: set[sp.Symbol] | set[str] | list[str]) -> str:
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


# Stage024 export symbols and the two coordinate-host metadata symbols.
a, c_s, rho_eff = sp.symbols("a c_s rho_eff", positive=True, real=True)
I25, Xi_Q = sp.symbols("I25 Xi_Q", nonzero=True, real=True)
eta_q, eta_phi = sp.symbols("eta_q eta_phi", nonzero=True, real=True)
varpi_q2, varpi_Phi2 = sp.symbols(
    "varpi_q2 varpi_Phi2", positive=True, real=True
)
lambda_c = sp.Symbol("lambda_c", real=True)
q2, Phi2 = sp.symbols("q2 Phi2", real=True)

# Retired vector coordinates and tag-only anti-fiat fixtures.
A_w, F_muw, Jw, U, W = sp.symbols("A_w F_muw Jw U W", nonzero=True)
Omega_U, Omega_W, R_mix, g_U, g_W = sp.symbols(
    "Omega_U Omega_W R_mix g_U g_W", nonzero=True
)
omega_wall, omega_rho, r_mix, g_rho, g_qold = sp.symbols(
    "omega_wall omega_rho r_mix g_rho g_qold", nonzero=True
)
sigma_hidden = sp.Symbol("sigma_hidden", nonzero=True)
free_carrier = sp.Symbol("free_carrier")
missing_rider = sp.Symbol("missing_rider")
foreign_subject = sp.Symbol("foreign_subject")

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
DENSITY_HOST_UNIVERSE = HOST_CONTRACT | {q2, Phi2}
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
RELABEL_SYMBOLS = {omega_wall, omega_rho, r_mix, g_rho, g_qold}

# Reshaped pathA_43 provenance ledger.  I25 is installed by source_tag_map
# from the cited moment_valid input, exactly as the source slice does.
BASE_SOURCE_TAGS: dict[sp.Symbol, set[str]] = {
    a: {"pathA_29_bulk", "pathA_32_wall"},
    c_s: {"pathA_29_bulk"},
    rho_eff: {"pathA_29_bulk"},
    Xi_Q: {"continuity_interface"},
    eta_q: {"continuity_interface"},
    eta_phi: {"continuity_interface"},
    varpi_q2: {"pathA_32_wall"},
    varpi_Phi2: {"pathA_29_bulk"},
    lambda_c: {
        "continuity_interface",
        "pathA_29_bulk",
        "pathA_32_wall",
    },
    free_carrier: set(),
    sigma_hidden: {"vector_port"},
    **{symbol: {"vector_port"} for symbol in VECTOR_SYMBOLS},
    **{symbol: {"vector_port"} for symbol in RELABEL_SYMBOLS},
}


def cited_N0_den() -> sp.Expr:
    """Cite, rather than re-derive, stage024's canonical factored export."""
    return (
        I25**2
        * Xi_Q**2
        * c_s**4
        * rho_eff
        * (eta_phi * varpi_q2 + eta_q * lambda_c) ** 2
        / (a**7 * (lambda_c**2 - varpi_Phi2 * varpi_q2) ** 2)
    )


def source_tag_map(
    moment_symbol: sp.Symbol,
    moment_valid: bool,
    *,
    free_carrier_tagged: bool = False,
) -> dict[sp.Symbol, set[str]]:
    tags = {symbol: set(values) for symbol, values in BASE_SOURCE_TAGS.items()}
    if free_carrier_tagged:
        tags[free_carrier] = {TAGGED_CARRIER_TAG}
    if moment_valid:
        tags[moment_symbol] = {"continuity_interface", "pathA_32_wall"}
    else:
        tags.setdefault(moment_symbol, set())
    return tags


def taint_for_expr(
    expr: sp.Expr, tag_map: dict[sp.Symbol, set[str]]
) -> tuple[set[str], set[sp.Symbol]]:
    """Union provenance tags over the expression's actual free symbols."""
    taint: set[str] = set()
    missing: set[sp.Symbol] = set()
    for symbol in expr.free_symbols:
        if symbol not in tag_map:
            missing.add(symbol)
        else:
            taint.update(tag_map[symbol])
    return taint, missing


@dataclass(frozen=True)
class SourceGraph:
    live_symbols: frozenset[sp.Symbol]
    taint_set: frozenset[str]
    missing_source_symbols: frozenset[sp.Symbol]
    empty_source_symbols: frozenset[sp.Symbol]
    vector_host_symbols: frozenset[sp.Symbol]

    @property
    def source_map_complete(self) -> bool:
        return not self.missing_source_symbols and not self.empty_source_symbols

    @property
    def vector_free(self) -> bool:
        return (
            self.source_map_complete
            and not self.vector_host_symbols
            and "vector_port" not in self.taint_set
        )


def source_graph_for_expr(
    expr: sp.Expr, tag_map: dict[sp.Symbol, set[str]]
) -> SourceGraph:
    symbols = frozenset(expr.free_symbols)
    taint, missing = taint_for_expr(expr, tag_map)
    empty = {symbol for symbol in symbols if symbol in tag_map and not tag_map[symbol]}
    return SourceGraph(
        live_symbols=symbols,
        taint_set=frozenset(taint),
        missing_source_symbols=frozenset(missing),
        empty_source_symbols=frozenset(empty),
        vector_host_symbols=frozenset(symbols & VECTOR_SYMBOLS),
    )


def vector_ablated_expr(
    expr: sp.Expr, tag_map: dict[sp.Symbol, set[str]]
) -> sp.Expr:
    vector_tainted = {
        symbol for symbol, tags in tag_map.items() if "vector_port" in tags
    }
    ablate = (VECTOR_SYMBOLS | vector_tainted) & expr.free_symbols
    return compact(expr.subs({symbol: 0 for symbol in ablate}))


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


def arity_scan(call_arity: int) -> bool:
    signature = inspect.signature(source_graph_for_expr)
    positional = [
        parameter
        for parameter in signature.parameters.values()
        if parameter.kind
        in (parameter.POSITIONAL_ONLY, parameter.POSITIONAL_OR_KEYWORD)
    ]
    return call_arity == len(positional)


AUTHORED_LEAK = sp.Function("source_graph_for_expr")


def leakage_scan(objects: list[Any]) -> bool:
    return not any(
        isinstance(obj, sp.Basic) and obj.has(AUTHORED_LEAK)
        for obj in objects
    )


def exercise_rig(
    label: str,
    assertion_name: str,
    rig_assertion: Callable[[], None],
    neutral_assertion: Callable[[], None],
) -> str:
    caught, fired_name = probe_assertion(rig_assertion, assertion_name)
    neutral_pass = assertion_passes(neutral_assertion)
    expect_bool(
        f"META {label} routed assertion fires and neutering stops it",
        caught and neutral_pass,
    )
    outcome = f"CAUGHT at {fired_name}; neutralized=PASS"
    print(f"RIG {label}: {outcome}")
    return outcome


def run_audit() -> dict[str, Any]:
    moment_valid = True  # Typed forward reference; stage026 earns this fact.
    n0_den = compact(cited_N0_den())
    tag_map = source_tag_map(I25, moment_valid)
    graph = source_graph_for_expr(n0_den, tag_map)

    subbanner("I. cited stage024 subject-integrity contract")
    expect_bool(
        "I cited N0_den live-symbol contract equals stage024's exact 10-symbol export",
        set(n0_den.free_symbols) == HOST_CONTRACT,
    )
    print(
        "DIAGNOSTIC (not counted): density-host metadata structural relationship: "
        "DENSITY_HOST_UNIVERSE minus HOST_CONTRACT = "
        f"{fmt_names(DENSITY_HOST_UNIVERSE - HOST_CONTRACT)} (10-vs-12)"
    )

    subbanner("P1/P2. computed provenance-taint proof")
    baseline_ancestry_ok = set(graph.taint_set) == BASELINE_TAGS
    expect_bool(
        "P1 baseline_ancestry_ok has exactly the three density tags",
        baseline_ancestry_ok,
    )
    expect_bool(
        "P2 source_map_complete has no missing or empty provenance",
        graph.source_map_complete,
    )
    expect_bool(
        "P2 vector_host_symbols is empty",
        not graph.vector_host_symbols,
    )
    expect_bool(
        "P2 computed taint excludes vector_port",
        "vector_port" not in graph.taint_set,
    )
    expect_bool(
        "P2 vector_free combines completeness, hosts, and computed taint",
        graph.vector_free,
    )
    print(
        "DIAGNOSTIC (not counted): moment_valid=True is a typed forward reference "
        "cited from stage026; the LOCAL verdict is conditional on it"
    )

    # Redundant witness only: deliberately not an expect_* gate.
    ablated = vector_ablated_expr(n0_den, tag_map)
    ablation_delta = compact(n0_den - ablated)

    subbanner("A-I. routed rig assertions and coupling meta-tests")
    outcomes: dict[str, str] = {}

    relabel_p = omega_wall**2 * g_rho + r_mix * g_qold
    relabel_delta = omega_wall**2 * omega_rho**2 - r_mix**2
    relabel_expr = compact(relabel_p**2 / relabel_delta**2)
    relabel_graph = source_graph_for_expr(relabel_expr, tag_map)
    relabel_neutered_tags = {symbol: set(tags) for symbol, tags in tag_map.items()}
    for symbol in RELABEL_SYMBOLS:
        relabel_neutered_tags[symbol] = {"continuity_interface"}
    relabel_neutral = source_graph_for_expr(relabel_expr, relabel_neutered_tags)
    outcomes["A"] = exercise_rig(
        "A relabel_rig",
        "A relabel_rig computed-taint assert: vector_port absent",
        lambda: routed_assert(
            "A relabel_rig computed-taint assert: vector_port absent",
            "vector_port" not in relabel_graph.taint_set,
        ),
        lambda: routed_assert(
            "A relabel_rig computed-taint assert: vector_port absent",
            "vector_port" not in relabel_neutral.taint_set,
        ),
    )

    hidden_expr = compact(n0_den * sigma_hidden)
    hidden_graph = source_graph_for_expr(hidden_expr, tag_map)
    hidden_neutered_tags = {symbol: set(tags) for symbol, tags in tag_map.items()}
    hidden_neutered_tags[sigma_hidden] = {"continuity_interface"}
    hidden_neutral = source_graph_for_expr(hidden_expr, hidden_neutered_tags)
    outcomes["B"] = exercise_rig(
        "B hidden_vector",
        "B hidden_vector computed-taint assert: vector_port absent",
        lambda: routed_assert(
            "B hidden_vector computed-taint assert: vector_port absent",
            "vector_port" not in hidden_graph.taint_set,
        ),
        lambda: routed_assert(
            "B hidden_vector computed-taint assert: vector_port absent",
            "vector_port" not in hidden_neutral.taint_set,
        ),
    )

    injected_graph = source_graph_for_expr(
        compact(n0_den * Omega_U / Omega_W), tag_map
    )
    outcomes["C"] = exercise_rig(
        "C vector_injection",
        "C vector_injection vector-host assert: vector_host_symbols empty",
        lambda: routed_assert(
            "C vector_injection vector-host assert: vector_host_symbols empty",
            not injected_graph.vector_host_symbols,
        ),
        lambda: routed_assert(
            "C vector_injection vector-host assert: vector_host_symbols empty",
            not graph.vector_host_symbols,
        ),
    )

    rider_expr = compact(n0_den * free_carrier)
    rider_graph = source_graph_for_expr(rider_expr, tag_map)
    tagged_map = source_tag_map(I25, moment_valid, free_carrier_tagged=True)
    tagged_graph = source_graph_for_expr(rider_expr, tagged_map)
    outcomes["D"] = exercise_rig(
        "D provenance_less_rider",
        "D provenance_less_rider source_map_complete assert",
        lambda: routed_assert(
            "D provenance_less_rider source_map_complete assert",
            rider_graph.source_map_complete,
        ),
        lambda: routed_assert(
            "D provenance_less_rider source_map_complete assert",
            tagged_graph.source_map_complete,
        ),
    )

    missing_expr = compact(n0_den * missing_rider)
    missing_graph = source_graph_for_expr(missing_expr, tag_map)
    repaired_map = {symbol: set(tags) for symbol, tags in tag_map.items()}
    repaired_map[missing_rider] = {"pathA_34_dimensionless_free_carrier"}
    repaired_graph = source_graph_for_expr(missing_expr, repaired_map)
    outcomes["E"] = exercise_rig(
        "E missing_symbol",
        "E missing_symbol source_map_complete assert",
        lambda: routed_assert(
            "E missing_symbol source_map_complete assert",
            missing_graph.source_map_complete,
        ),
        lambda: routed_assert(
            "E missing_symbol source_map_complete assert",
            repaired_graph.source_map_complete,
        ),
    )

    expected_four_tags = BASELINE_TAGS | {TAGGED_CARRIER_TAG}
    f_assert_name = "F tagged_carrier source_map_complete flip assert"
    f_flip_caught, f_flip_name = probe_assertion(
        lambda: routed_assert(f_assert_name, rider_graph.source_map_complete),
        f_assert_name,
    )
    expect_bool(
        "F tagged_carrier passes P2 with four tags and flips when tag is neutered",
        tagged_graph.vector_free
        and set(tagged_graph.taint_set) == expected_four_tags
        and not (set(tagged_graph.taint_set) == BASELINE_TAGS)
        and f_flip_caught,
    )
    outcomes["F"] = (
        "PASS P2; taint="
        + fmt_names(expected_four_tags)
        + f"; tag-neutered=FLIP at {f_flip_name}"
    )
    print(f"RIG F tagged_carrier: {outcomes['F']}")

    old_port = compact(Omega_U**2 * g_W + R_mix * g_U)
    old_port_graph = source_graph_for_expr(old_port, tag_map)
    outcomes["G"] = exercise_rig(
        "G raw_vector_port",
        "G raw_vector_port vector-host assert: vector_host_symbols empty",
        lambda: routed_assert(
            "G raw_vector_port vector-host assert: vector_host_symbols empty",
            not old_port_graph.vector_host_symbols,
        ),
        lambda: routed_assert(
            "G raw_vector_port vector-host assert: vector_host_symbols empty",
            not graph.vector_host_symbols,
        ),
    )

    outcomes["H_arity"] = exercise_rig(
        "H' arity_scanner",
        "H' definition/call arity scanner assert",
        lambda: routed_assert(
            "H' definition/call arity scanner assert", arity_scan(1)
        ),
        lambda: routed_assert(
            "H' definition/call arity scanner assert", arity_scan(2)
        ),
    )
    actual_transcript = [n0_den, ablated, sorted(graph.taint_set)]
    leaked_transcript = actual_transcript + [AUTHORED_LEAK(n0_den)]
    outcomes["H_leak"] = exercise_rig(
        "H' leakage_scanner",
        "H' unevaluated-leakage transcript scanner assert",
        lambda: routed_assert(
            "H' unevaluated-leakage transcript scanner assert",
            leakage_scan(leaked_transcript),
        ),
        lambda: routed_assert(
            "H' unevaluated-leakage transcript scanner assert",
            leakage_scan(actual_transcript),
        ),
    )

    corrupt_subject = compact(n0_den.subs(rho_eff, foreign_subject))
    outcomes["I"] = exercise_rig(
        "I subject_integrity",
        "I subject_integrity exact host-contract assert",
        lambda: routed_assert(
            "I subject_integrity exact host-contract assert",
            set(corrupt_subject.free_symbols) == HOST_CONTRACT,
        ),
        lambda: routed_assert(
            "I subject_integrity exact host-contract assert",
            set(n0_den.free_symbols) == HOST_CONTRACT,
        ),
    )

    expect_bool(
        "H' actual transcript has no unevaluated authored-helper leakage",
        leakage_scan(actual_transcript),
    )
    local_ok = bool(moment_valid and baseline_ancestry_ok and graph.vector_free)
    expect_bool("LOCAL DENSITY_PORT_VECTOR_FREE = P1 and P2", local_ok)

    return {
        "N0_den": n0_den,
        "live_symbols": set(n0_den.free_symbols),
        "graph": graph,
        "ablation_delta": ablation_delta,
        "ablation_witness_ok": ablation_delta == 0,
        "outcomes": outcomes,
        "moment_valid": moment_valid,
        "local_ok": local_ok,
    }


def emit(data: dict[str, Any]) -> None:
    graph: SourceGraph = data["graph"]
    subbanner("Computed density-port vector-freedom transcript")
    print("consumes: stage024 N0_den (cited canonical factored export)")
    print(f"N0_den (canonical factored): {fmt(data['N0_den'])}")
    print(f"N0_den live symbols: {fmt_names(data['live_symbols'])}")
    print(f"host contract (10): {fmt_names(HOST_CONTRACT)}")
    print(f"allowable density-host universe (12): {fmt_names(DENSITY_HOST_UNIVERSE)}")
    print(f"computed taint set: {fmt_names(set(graph.taint_set))}")
    print(f"missing_source_symbols: {fmt_names(set(graph.missing_source_symbols))}")
    print(f"empty_source_symbols: {fmt_names(set(graph.empty_source_symbols))}")
    print(f"vector_host_symbols: {fmt_names(set(graph.vector_host_symbols))}")
    print(f"source_map_complete: {graph.source_map_complete}")
    print(f"baseline_ancestry_ok (P1): {set(graph.taint_set) == BASELINE_TAGS}")
    print(f"vector_free (P2): {graph.vector_free}")
    print(f"ablation delta (redundant witness only): {fmt(data['ablation_delta'])}")
    print(f"ablation witness invariant: {data['ablation_witness_ok']} (de-counted)")
    print(
        "moment_valid: True (typed forward reference; LOCAL verdict is conditional "
        "until stage026 discharges it)"
    )

    subbanner("Verdict labels")
    print(f"LOCAL_AUDIT_VERDICT: {LOCAL_VERDICT}")
    print(f"JOINT_LANDING_LABEL (PARTIAL): {JOINT_PARTIAL}")


def main() -> int:
    banner("ledger_stage025_vector_freedom_taint_sympy_audit")
    print("Target stem confirmed: ledger_stage025_vector_freedom_taint")
    print(
        "Engine: computed SymPy provenance taint over actual free_symbols; "
        "zero file I/O."
    )
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
