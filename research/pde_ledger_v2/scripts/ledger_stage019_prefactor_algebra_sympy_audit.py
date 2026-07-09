#!/usr/bin/env python3
"""Ledger stage019 SymPy audit: squared-denominator prefactor algebra.

Standalone, print-only, no arguments, no file I/O. This is the pathA_33
II-G4b abstract-algebra slice only: series-extract P0/P2/P4 from
D0*N(omega)/Dcons(omega)^2, discriminate that object from plain N/D, and
exercise the local 3g failure with a dynamic rerun. Port moments and the
adjacent stage products are cited as provenance only.
"""

from __future__ import annotations

from typing import Any

import sympy as sp


PASS_COUNT = 0
FAIL_COUNT = 0

PREFACTOR_ALGEBRA_EARNED = "PREFACTOR_ALGEBRA_EARNED"
FAIL_PREFACTOR_ALGEBRA = "FAIL_PREFACTOR_ALGEBRA"
NO_FAIL = "NO_FAIL"


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


def bool_from_residual(residual: sp.Expr | int) -> bool:
    assert_no_float("bool_from_residual", residual)
    return compact(residual) == 0


def verdict_residual(actual: str, expected: str) -> sp.Integer:
    return sp.Integer(0) if actual == expected else sp.Integer(1)


omega = sp.Symbol("omega", real=True)
D0, D2, D4, N0, N2, N4 = sp.symbols(
    "D0 D2 D4 N0 N2 N4", real=True, nonzero=True
)


def series_no_o(expr: sp.Expr, var: sp.Symbol, order: int) -> sp.Expr:
    return sp.expand(sp.series(expr, var, 0, order).removeO())


def extract_coefficients(obj: sp.Expr) -> dict[str, Any]:
    """Derive the coefficient slots from the actual symbolic series."""

    derived_series = series_no_o(obj, omega, 6)
    return {
        "series": derived_series,
        "P0": compact(derived_series.coeff(omega, 0)),
        "P2": compact(derived_series.coeff(omega, 2)),
        "P4": compact(derived_series.coeff(omega, 4)),
    }


def scoped_verdict(prefactor_ok: bool, nd_self_check_ok: bool) -> str:
    """019-local gate only: coefficient match AND the N/D discriminator."""

    if not prefactor_ok or not nd_self_check_ok:
        return FAIL_PREFACTOR_ALGEBRA
    return PREFACTOR_ALGEBRA_EARNED


def dynamic_ablation(
    baseline_prefactor_ok: bool,
    mutated_prefactor_ok: bool,
    nd_self_check_ok: bool,
) -> dict[str, Any]:
    """Rerun the 019-local gate with and without the plain-N/D mutation."""

    with_mutation = scoped_verdict(mutated_prefactor_ok, nd_self_check_ok)
    without_mutation = scoped_verdict(baseline_prefactor_ok, nd_self_check_ok)
    return {
        "rerun_gate_logic": True,
        "with_mutation": with_mutation,
        "without_mutation": without_mutation,
        "expected_fail": FAIL_PREFACTOR_ALGEBRA,
        "fail_suppressed": (
            with_mutation == FAIL_PREFACTOR_ALGEBRA
            and without_mutation != FAIL_PREFACTOR_ALGEBRA
        ),
    }


def build_prefactor() -> dict[str, Any]:
    numerator = N0 + N2 * omega**2 + N4 * omega**4
    denominator = D0 + D2 * omega**2 + D4 * omega**4
    correct_obj = D0 * numerator / denominator**2
    plain_obj = numerator / denominator

    correct = extract_coefficients(correct_obj)
    plain = extract_coefficients(plain_obj)
    coeffs = {name: correct[name] for name in ("P0", "P2", "P4")}
    plain_coeffs = {name: plain[name] for name in ("P0", "P2", "P4")}
    expected = {
        "P0": N0 / D0,
        "P2": (D0 * N2 - 2 * D2 * N0) / D0**2,
        "P4": (
            D0**2 * N4
            - 2 * D0 * (D2 * N2 + D4 * N0)
            + 3 * D2**2 * N0
        )
        / D0**3,
    }
    residuals = {
        name: compact(coeffs[name] - expected[name])
        for name in ("P0", "P2", "P4")
    }
    matches = {
        name: bool_from_residual(residuals[name])
        for name in ("P0", "P2", "P4")
    }
    plain_p2_difference = compact(plain_coeffs["P2"] - expected["P2"])
    plain_equals_correct_p2 = bool_from_residual(plain_p2_difference)
    sample_subs = {
        D0: sp.Integer(19),
        D2: sp.Integer(23),
        D4: sp.Integer(29),
        N0: sp.Integer(11),
        N2: sp.Integer(13),
        N4: sp.Integer(17),
    }
    sample_values = {
        "P0": compact(coeffs["P0"].subs(sample_subs)),
        "P2": compact(coeffs["P2"].subs(sample_subs)),
        "P4": compact(coeffs["P4"].subs(sample_subs)),
        "plainP2": compact(plain_coeffs["P2"].subs(sample_subs)),
        "plainMinusCorrectP2": compact(plain_p2_difference.subs(sample_subs)),
    }
    return {
        "numerator": numerator,
        "denominator": denominator,
        "correct_object": correct_obj,
        "plain_object": plain_obj,
        "correct_series": correct["series"],
        "plain_series": plain["series"],
        "coefficients": coeffs,
        "plain_coefficients": plain_coeffs,
        "expected": expected,
        "residuals": residuals,
        "matches": matches,
        "ok": all(matches.values()),
        "plain_equals_correct_P2": plain_equals_correct_p2,
        "difference_plain_minus_correct_P2": plain_p2_difference,
        "sample_subs": sample_subs,
        "sample_values": sample_values,
    }


def build_probe(prefactor: dict[str, Any]) -> dict[str, Any]:
    expected = prefactor["expected"]
    plain = prefactor["plain_coefficients"]
    mutated_matches = {
        name: bool_from_residual(plain[name] - expected[name])
        for name in ("P0", "P2", "P4")
    }
    mutated_prefactor_ok = all(mutated_matches.values())
    nd_self_check_ok = not prefactor["plain_equals_correct_P2"]
    correct_equals_correct_p2 = bool_from_residual(
        prefactor["coefficients"]["P2"] - expected["P2"]
    )
    return {
        "mutated_object": prefactor["plain_object"],
        "mutated_prefactor_ok": mutated_prefactor_ok,
        "plain_equals_correct_P2": prefactor["plain_equals_correct_P2"],
        "verdict": (
            FAIL_PREFACTOR_ALGEBRA
            if not prefactor["plain_equals_correct_P2"]
            else NO_FAIL
        ),
        "correct_object_verdict": (
            NO_FAIL if correct_equals_correct_p2 else FAIL_PREFACTOR_ALGEBRA
        ),
        "self_ablation": dynamic_ablation(
            prefactor["ok"], mutated_prefactor_ok, nd_self_check_ok
        ),
    }


def build_baseline() -> dict[str, Any]:
    prefactor = build_prefactor()
    probe = build_probe(prefactor)
    nd_self_check_ok = not prefactor["plain_equals_correct_P2"]
    return {
        "prefactor": prefactor,
        "probe": probe,
        "gates": {
            "prefactor_match": prefactor["ok"],
            "N_over_D_self_check": nd_self_check_ok,
        },
        "verdict": scoped_verdict(prefactor["ok"], nd_self_check_ok),
    }


def run_prefactor_derivation(data: dict[str, Any]) -> None:
    prefactor = data["prefactor"]
    subbanner("Squared-denominator object and series-extracted coefficients")
    print(f"  N(omega) = {fmt(prefactor['numerator'])}")
    print(f"  Dcons(omega) = {fmt(prefactor['denominator'])}")
    print(f"  correctObj = {fmt(prefactor['correct_object'])}")
    print(f"  correctSeries = {fmt(prefactor['correct_series'])}")
    print("  P0/P2/P4 below are read by coeff(omega,n) from correctSeries; targets are independent typed references.")
    for name in ("P0", "P2", "P4"):
        print(f"  {name}(series-extracted) = {fmt(prefactor['coefficients'][name])}")
        expect_zero(
            f"series-extracted {name} matches independent expected {name}",
            prefactor["residuals"][name],
        )
    print(
        "  exact sample controls (D0,D2,D4,N0,N2,N4 only) = "
        + fmt(prefactor["sample_subs"])
    )
    print(f"  exact sample values (transcript only, non-tooth) = {fmt(prefactor['sample_values'])}")


def run_nd_self_check(data: dict[str, Any]) -> None:
    prefactor = data["prefactor"]
    plain = prefactor["plain_coefficients"]
    expected_plain_p2 = (D0 * N2 - D2 * N0) / D0**2
    expected_gap = D2 * N0 / D0**2
    subbanner("Plain N/D factor-of-2 discriminator")
    print(f"  plainObj = {fmt(prefactor['plain_object'])}")
    print(f"  plainSeries = {fmt(prefactor['plain_series'])}")
    print(f"  plainP2(series-extracted) = {fmt(plain['P2'])}")
    print(
        "  differencePlainMinusCorrectP2(computed) = "
        + fmt(prefactor["difference_plain_minus_correct_P2"])
    )
    print(
        "  plainEqualsCorrectP2 = "
        + fmt(prefactor["plain_equals_correct_P2"])
    )
    expect_zero(
        "plain N/D P2 is series-extracted with the single D2*N0 term",
        plain["P2"] - expected_plain_p2,
    )
    expect_zero(
        "computed plain-minus-correct P2 is D2*N0/D0^2",
        prefactor["difference_plain_minus_correct_P2"] - expected_gap,
    )
    expect_bool(
        "plainEqualsCorrectP2 is False (squared-denominator factor-of-2 detected)",
        not prefactor["plain_equals_correct_P2"],
    )


def run_probe(data: dict[str, Any]) -> None:
    probe = data["probe"]
    ablation = probe["self_ablation"]
    subbanner("Probe 3g wrong-prefactor object and dynamic 019-local self-ablation")
    print(f"  3g mutated object = {fmt(probe['mutated_object'])}")
    print(f"  3g plainEqualsCorrectP2 = {fmt(probe['plain_equals_correct_P2'])}")
    print(f"  3g wrong-object verdict = {probe['verdict']}")
    print(f"  3g correct-object verdict = {probe['correct_object_verdict']}")
    print(f"  3g dynamic self-ablation = {ablation}")
    expect_zero(
        "3g plain N/D reaches FAIL_PREFACTOR_ALGEBRA",
        verdict_residual(probe["verdict"], FAIL_PREFACTOR_ALGEBRA),
    )
    expect_zero(
        "3g correct squared-denominator object is NO_FAIL",
        verdict_residual(probe["correct_object_verdict"], NO_FAIL),
    )
    expect_zero(
        "3g dynamic rerun with plain-N/D mutation fails locally",
        verdict_residual(ablation["with_mutation"], FAIL_PREFACTOR_ALGEBRA),
    )
    expect_zero(
        "3g dynamic rerun without mutation returns the earned local verdict",
        verdict_residual(ablation["without_mutation"], PREFACTOR_ALGEBRA_EARNED),
    )
    expect_bool("3g self-ablation suppresses the failure", ablation["fail_suppressed"])


def run_per_tooth_ablations(data: dict[str, Any]) -> None:
    prefactor = data["prefactor"]
    expected = prefactor["expected"]
    subbanner("Per-tooth derivation ablations on independent series copies")
    for name, power in (("P0", 0), ("P2", 2), ("P4", 4)):
        mutated = extract_coefficients(prefactor["correct_object"] + omega**power)
        expect_fail(
            f"{name} derivation copy changes its own series-extracted coefficient",
            mutated[name] - expected[name],
        )
    correct_control = extract_coefficients(prefactor["correct_object"])
    correct_control_equals = bool_from_residual(
        correct_control["P2"] - expected["P2"]
    )
    expect_bool(
        "N/D discriminator comparison flips True when its plain control is replaced by the correct object",
        correct_control_equals,
    )
    expect_zero(
        "baseline remains immutable after derivation-copy ablations",
        sum(prefactor["residuals"].values(), sp.Integer(0)),
    )


def live_symbol_names(data: dict[str, Any]) -> set[str]:
    prefactor = data["prefactor"]
    expressions = [
        prefactor["numerator"],
        prefactor["denominator"],
        prefactor["correct_object"],
        prefactor["plain_object"],
        prefactor["correct_series"],
        prefactor["plain_series"],
        *prefactor["coefficients"].values(),
        *prefactor["plain_coefficients"].values(),
    ]
    return {
        symbol.name
        for expr in expressions
        for symbol in sp.sympify(expr).free_symbols
    }


def run_scope_and_provenance(data: dict[str, Any]) -> None:
    subbanner("019 scope, provenance-only consumption, and PARTIAL landing")
    print(f"  019 gate booleans = {data['gates']}")
    print(f"  019 scoped verdict = {data['verdict']}")
    expect_zero(
        "019 scoped verdict lands the earned partial component",
        verdict_residual(data["verdict"], PREFACTOR_ALGEBRA_EARNED),
    )
    names = live_symbol_names(data)
    print(f"  live symbolic names in the earned slice = {sorted(names)}")
    expect_bool(
        "earned algebra has exactly omega and D0..N4 as live symbols",
        names == {"omega", "D0", "D2", "D4", "N0", "N2", "N4"},
    )
    print("  QUAD_CALIBRATED (2/4) -- squared-denominator prefactor algebra EARNED (PARTIAL)")
    print("    = P0=N0/D0, P2=(D0*N2-2*D2*N0)/D0^2, P4=(D0^2*N4-2*D0*(D2*N2+D4*N0)+3*D2^2*N0)/D0^3.")
    print("    SERIES-EXTRACTED from P(omega)=D0*N/Dcons^2; the factor-of-2 is the squared-denominator signature.")
    print("    AND plain N/D reaches FAIL_PREFACTOR_ALGEBRA while the correct object is NO_FAIL.")
    print("  REMAINING: fingerprint = 018 (DONE); magnitude partition plus calibrated label = 020; dimension closure = 021.")
    print("  CAVEATS: D-lane and N-moment numerical branch scalars remain symbolic/Gate-6 deferred; the magnitude is calibrated in 020.")
    print("  CONSUMES (PROVENANCE only): 017 D-lanes + deferred port N-moments + 018 fingerprint context; no guard/dual-site.")
    print("  EXPORTS: the series-extracted prefactor algebra and N/D self-check to 020/021 and 027; no file artifact is written.")
    print("  PORT-MOMENTS PROVENANCE EXPORT -- LABELED NON-CHECK: concrete N-moments are deferred Gate-6 branch data, asserted against nothing.")
    print("  reduction certificate: FROZEN-INPUT symbolic D_n lanes and outgoing port-moment inputs/N_n; COMPUTED P0/P2/P4 plus N/D self-check; DEFERRED numerical branch scalars.")
    print("  dropped-bookkeeping: scratch-YAML agreement handoff and report/feed writers are absent; agreement is transcript-level.")
    print("  register note: likely zero new counted knobs; registration decides the structural edge.")


def print_verdict_labels() -> None:
    subbanner("Verdict labels")
    print("  ledger earned-label (NOT a source verdict token): PREFACTOR_ALGEBRA_EARNED  (the squared-denominator prefactor P(omega)=D0*N(omega)/Dcons(omega)^2, N=N0+N2*omega^2+N4*omega^4, Dcons=D0+D2*omega^2+D4*omega^4, series-expands to the coefficients P0=N0/D0, P2=(D0*N2-2*D2*N0)/D0^2, P4=(D0^2*N4-2*D0*(D2*N2+D4*N0)+3*D2^2*N0)/D0^3 -- SERIES-EXTRACTED, NOT typed; the -2*D2*N0 factor-of-2 is the squared-denominator signature (plain N/D gives only -D2*N0, provably wrong -> FAIL_PREFACTOR_ALGEBRA))")
    print("  source top-line verdict: QUAD_CALIBRATED  (JOINT 4-stage; 019 carries the squared-denominator prefactor-algebra component 2/4 and lands the token as a PARTIAL)")
    print("  joint composition (019 = the SECOND leg; 018 DONE, 020/021 remaining): QUAD_CALIBRATED = (018: outgoing DtN Hankel fingerprint + chi_Q sign)[EARNED, PARTIAL, DONE] AND (019: squared-denominator prefactor algebra P(omega)=D0*N/Dcons^2)[EARNED here, PARTIAL] AND (020: 54/5=2*27/5 provenance partition + the CALIBRATED label, G=GENUINE_BLOCKED) AND (021: mu_hat0-free [P0_phys]=1 dim closure)")
    print("  earned (the prefactor algebra): P0=N0/D0, P2=(D0*N2-2*D2*N0)/D0^2, P4=(D0^2*N4-2*D0*(D2*N2+D4*N0)+3*D2^2*N0)/D0^3 SERIES-EXTRACTED from P(omega)=D0*N/Dcons^2 (NOT typed); the N/D self-check (plain N/D gives -D2*N0 vs the correct -2*D2*N0 -> FAIL_PREFACTOR_ALGEBRA; the correct object does NOT fire)")
    print("  earned (able-to-fail): 3g_wrong_prefactor_object (plain N/D -> missing factor of 2 -> FAIL_PREFACTOR_ALGEBRA; correct object NO_FAIL), with a DYNAMIC self-ablation (re-run, NOT a constant); the prefactor coefficients are SERIES-EXTRACTED off the actual series (not a hardcode-and-compare, the section-4 firewall)")
    print("  calibrated / deferred (NOT 019): the fingerprint (018, DONE); the 54/5=2*27/5 magnitude + G (020, external_bridge_input, G=GENUINE_BLOCKED); the mu_hat0-free dim closure (021); the numerical (D_n,N_n) port scalars (Gate-6 sim-deferred, report :49)")
    print("  consumed (PROVENANCE only -- NO guard/dual-site; abstract port-agnostic algebra): 017's l=2 port kernel (D-lanes D0/D2/D4 + {Btilde,Ztilde}, carried as abstract symbols) + build_port_moments' concrete N-moments (deferred Gate-6 branch data, asserted-against-nothing) + 018's fingerprint context; NO c_s/a/G/mu_hat0")
    print("  exports (the prefactor algebra): P0=N0/D0 (-> 021's [P0_phys]=1 dim closure) + P2/P4 + the N/D self-check (-> 020's 54/5=2*27/5 partition context) => 020/021 + 027 (pathA_43 closure)")
    print("  new symbols first-appearing (019): none new-counted (the D0/D2/D4 are 017's exported D-lanes + the N0/N2/N4 are build_port_moments' deferred Gate-6 port N-moments, NOT new knobs); no units symbol (units-free abstract algebra); no counted knobs expected (an EARNED prefactor-algebra slice, like 018)")


def main() -> int:
    banner("ledger_stage019_prefactor_algebra_sympy_audit")
    print("Target stem confirmed: ledger_stage019_prefactor_algebra")
    print("Engine: SymPy exact symbolic Series/coeff derivation; no floats/tolerances; zero file I/O.")
    data = build_baseline()
    run_prefactor_derivation(data)
    run_nd_self_check(data)
    run_probe(data)
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
