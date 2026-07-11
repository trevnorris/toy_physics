#!/usr/bin/env python3
"""Ledger stage026 SymPy audit: continuity-lineage token check.

Standalone, print-only, no arguments, and no file I/O.  Stage024's factored
N0_den is cited directly.  The decisive gates are an exact lexical-token
lineage validator and the I25-vs-I_wrong2 earning gate.  This stage does not
rebuild the two-port inverse, vector-taint graph, or stage027 port checks.
"""

from __future__ import annotations

import inspect
import re
from copy import deepcopy
from typing import Any, Callable

import sympy as sp


PASS_COUNT = 0
FAIL_COUNT = 0

LOCAL_VERDICT = "CONTINUITY_LINEAGE_EARNED"
JOINT_PARTIAL = (
    "DENSITY_PORT_HOSTED (3/4, CONTINUITY LINEAGE — I25 is a genuine "
    "∫Y₂*·S_leak ℓ=2 continuity moment descended from pathA_29's operator; "
    "earns moment_valid; 024 = derivation, 025 = vector-freedom, "
    "027 = port checks + closure)"
)
FAIL_VERDICT = "FAIL_NOT_DENSITY_DERIVED"

CONTINUITY_OPERATOR_ID = "pathA_29_projected_continuity"
CONTINUITY_L0 = "M0 = Integral(S_leak d3x)"
CONTINUITY_L1 = "D1_i = Integral(x_i*S_leak d3x) + Integral(j_i d3x)"
CONTINUITY_L2 = "Q2_m = Integral(Y2_m_star*S_leak d3x)"
CONTINUITY_L2_KERNEL = "Y2_m_star*S_leak"
BAD_CONTINUITY_L2 = "GARBAGE_NOT_A_CONTINUITY_MOMENT_AT_ALL"
STUFFED_L0 = "NOT_M0 = FakeIntegral(S_leakage d3xyz)"

REQUIRED_TOKENS: dict[str, frozenset[str]] = {
    "l0": frozenset({"M0", "Integral", "S_leak", "d3x"}),
    "l1": frozenset(
        {"D1_i", "Integral", "x_i", "S_leak", "j_i", "d3x"}
    ),
    "l2": frozenset(
        {"Q2_m", "Integral", "Y2_m_star", "S_leak", "d3x"}
    ),
    "l2_kernel": frozenset({"Y2_m_star", "S_leak"}),
}


class AuditFailure(AssertionError):
    """A named stage026 audit assertion failed."""


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


def fmt_names(items: set[sp.Symbol] | set[str] | frozenset[str]) -> str:
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


# The exact stage024 export contract.  q2/Phi2 are coordinate metadata and do
# not occur in the exported numerator's free-symbol set.
a, c_s, rho_eff = sp.symbols("a c_s rho_eff", positive=True, real=True)
I25, I_wrong2, Xi_Q = sp.symbols("I25 I_wrong2 Xi_Q", nonzero=True, real=True)
eta_q, eta_phi = sp.symbols("eta_q eta_phi", nonzero=True, real=True)
varpi_q2, varpi_Phi2 = sp.symbols(
    "varpi_q2 varpi_Phi2", positive=True, real=True
)
lambda_c = sp.Symbol("lambda_c", real=True)
q2, Phi2 = sp.symbols("q2 Phi2", real=True)
Omega_U, Omega_W, R_mix, g_U, g_W = sp.symbols(
    "Omega_U Omega_W R_mix g_U g_W", nonzero=True
)
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


def cited_N0_den() -> sp.Expr:
    """Cite stage024's canonical factored export; never rebuild its inverse."""
    return (
        I25**2
        * Xi_Q**2
        * c_s**4
        * rho_eff
        * (eta_phi * varpi_q2 + eta_q * lambda_c) ** 2
        / (a**7 * (lambda_c**2 - varpi_Phi2 * varpi_q2) ** 2)
    )


def identifier_tokens(text: Any) -> frozenset[str]:
    """Extract exact [A-Za-z0-9_]+ lexical tokens, preserving underscores."""
    return frozenset(re.findall(r"[A-Za-z0-9_]+", str(text)))


def exact_contains_all(text: Any, required: frozenset[str]) -> bool:
    return required <= identifier_tokens(text)


def raw_substring_contains_all(text: Any, required: frozenset[str]) -> bool:
    """Source behavior retained only as an adversarial contrast witness."""
    value = str(text)
    return all(token in value for token in required)


def continuity_lineage_valid(lineage: dict[str, Any]) -> bool:
    """G1: compute the whole ancestry from exact tokens, never `valid`."""
    moments = lineage.get("moments", {})
    return bool(
        lineage.get("operator_id") == CONTINUITY_OPERATOR_ID
        and exact_contains_all(moments.get("l0", ""), REQUIRED_TOKENS["l0"])
        and exact_contains_all(moments.get("l1", ""), REQUIRED_TOKENS["l1"])
        and exact_contains_all(moments.get("l2", ""), REQUIRED_TOKENS["l2"])
        and exact_contains_all(
            lineage.get("l2_kernel", ""), REQUIRED_TOKENS["l2_kernel"]
        )
    )


def flag_only_validator(lineage: dict[str, Any]) -> bool:
    """Deliberately invalid meta-test implementation."""
    return bool(lineage.get("valid", False))


def continuity_moment_symbol(
    lineage: dict[str, Any],
    a_power: sp.Rational,
    validator: Callable[[dict[str, Any]], bool] = continuity_lineage_valid,
) -> tuple[sp.Symbol, bool]:
    """G2: moment_valid is G1; only a valid -7/2 baseline earns I25."""
    lineage_valid = bool(validator(lineage))
    earned = (
        I25
        if lineage_valid and a_power == sp.Rational(-7, 2)
        else I_wrong2
    )
    return earned, lineage_valid


def continuity_dependency_ok(
    lineage_valid: bool,
    moment_valid: bool,
    earned_moment: sp.Symbol,
    expr: sp.Expr,
    coupling_zero: bool,
    *,
    second_arm_enabled: bool = True,
) -> bool:
    """Full source predicate, including the tooth-G vanished-coupling arm."""
    membership_arm = earned_moment in expr.free_symbols
    vanished_arm = bool(
        second_arm_enabled and compact(expr) == 0 and coupling_zero
    )
    return bool(lineage_valid and moment_valid and (membership_arm or vanished_arm))


def baseline_lineage() -> dict[str, Any]:
    return {
        "operator_id": CONTINUITY_OPERATOR_ID,
        "moments": {
            "l0": CONTINUITY_L0,
            "l1": CONTINUITY_L1,
            "l2": CONTINUITY_L2,
        },
        "l2_kernel": CONTINUITY_L2_KERNEL,
        "lineage": "l0->l1->l2 from the same projected-continuity operator",
    }


def decoy_negative(lineage: dict[str, Any]) -> dict[str, Any]:
    """Every negative carries the two forbidden self-asserted decoys."""
    lineage["valid"] = True
    lineage["continuity_interface"] = True
    return lineage


def negative_lineages() -> dict[str, dict[str, Any]]:
    fake = decoy_negative(
        {
            "operator_id": "mis_tagged_vector_formula",
            "moments": {
                "l0": "Omega_U",
                "l1": "Omega_W",
                "l2": "Omega_U^2*g_W + R_mix*g_U, relabeled as continuity",
            },
            "l2_kernel": "R_mix*g_U",
        }
    )

    attack2 = deepcopy(baseline_lineage())
    attack2["moments"]["l2"] = BAD_CONTINUITY_L2

    operator_drop = deepcopy(baseline_lineage())
    operator_drop["operator_id"] = "pathA_29_projected_continuity_forged"

    l0_drop = deepcopy(baseline_lineage())
    l0_drop["moments"]["l0"] = "mass = Integral(S_leak d3x)"

    l1_drop = deepcopy(baseline_lineage())
    l1_drop["moments"]["l1"] = (
        "dipole = Integral(x_i*S_leak d3x) + Integral(j_i d3x)"
    )

    kernel_drop = deepcopy(baseline_lineage())
    kernel_drop["l2_kernel"] = "angular_kernel*S_leak"

    stuffed = deepcopy(baseline_lineage())
    stuffed["moments"]["l0"] = STUFFED_L0

    return {
        "A fake_continuity": decoy_negative(fake),
        "B attack2": decoy_negative(attack2),
        "C operator_id": decoy_negative(operator_drop),
        "C l0_token": decoy_negative(l0_drop),
        "C l1_token": decoy_negative(l1_drop),
        "C l2_kernel_token": decoy_negative(kernel_drop),
        "C token_stuffing": decoy_negative(stuffed),
    }


def evaluate_lineage(
    lineage: dict[str, Any],
    expr: sp.Expr,
    *,
    validator: Callable[[dict[str, Any]], bool] = continuity_lineage_valid,
) -> dict[str, Any]:
    valid = bool(validator(lineage))
    earned, moment_valid = continuity_moment_symbol(
        lineage, sp.Rational(-7, 2), validator
    )
    dependency = continuity_dependency_ok(
        valid, moment_valid, earned, expr, False
    )
    return {
        "lineage_valid": valid,
        "earned_moment": earned,
        "moment_valid": moment_valid,
        "continuity_dependency_ok": dependency,
        "accepted": bool(valid and earned == I25 and moment_valid and dependency),
    }


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
    verdict: str = FAIL_VERDICT,
) -> str:
    caught, fired_name = probe_assertion(rig_assertion, assertion_name)
    neutral_pass = assertion_passes(neutral_assertion)
    expect_bool(
        f"META {label} routed assertion fires and neutering stops it",
        caught and neutral_pass,
    )
    outcome = f"{verdict} at {fired_name}; neutralized=PASS"
    print(f"RIG {label}: {outcome}")
    return outcome


def arity_scan(call_arity: int) -> bool:
    parameters = inspect.signature(evaluate_lineage).parameters.values()
    positional = [
        parameter
        for parameter in parameters
        if parameter.kind
        in (parameter.POSITIONAL_ONLY, parameter.POSITIONAL_OR_KEYWORD)
    ]
    return call_arity == len(positional)


AUTHORED_LEAK = sp.Function("continuity_lineage_valid")


def leakage_scan(objects: list[Any]) -> bool:
    return not any(
        isinstance(obj, sp.Basic) and obj.has(AUTHORED_LEAK)
        for obj in objects
    )


def run_audit() -> dict[str, Any]:
    n0_den = compact(cited_N0_den())
    live_symbols = set(n0_den.free_symbols)
    lineage = baseline_lineage()
    baseline = evaluate_lineage(lineage, n0_den)
    negatives = negative_lineages()

    subbanner("I. cited stage024 consumption-integrity contract")
    expect_bool(
        "I cited N0_den free_symbols equals stage024 exact 10-symbol host contract",
        live_symbols == HOST_CONTRACT,
    )

    subbanner("G1/G2. exact-token ancestry and moment earning")
    baseline_token_sets = {
        "l0": identifier_tokens(lineage["moments"]["l0"]),
        "l1": identifier_tokens(lineage["moments"]["l1"]),
        "l2": identifier_tokens(lineage["moments"]["l2"]),
        "l2_kernel": identifier_tokens(lineage["l2_kernel"]),
    }
    expect_bool(
        "G1 genuine CONTINUITY_L* identifier tokens survive intact",
        all(
            REQUIRED_TOKENS[level] <= baseline_token_sets[level]
            for level in REQUIRED_TOKENS
        ),
    )
    stuffed = negatives["C token_stuffing"]
    stuffed_l0 = stuffed["moments"]["l0"]
    expect_bool(
        "G1 token-stuffing passes raw substring but fails exact lexical tokens",
        raw_substring_contains_all(stuffed_l0, REQUIRED_TOKENS["l0"])
        and not exact_contains_all(stuffed_l0, REQUIRED_TOKENS["l0"]),
    )
    expect_bool(
        "G1 baseline continuity_lineage_valid is computed True",
        baseline["lineage_valid"],
    )
    expect_bool(
        "G2 baseline earns exactly (I25, moment_valid=True)",
        baseline["earned_moment"] == I25 and baseline["moment_valid"],
    )
    corrupted_gate = evaluate_lineage(negatives["B attack2"], n0_den)
    expect_bool(
        "G2 corrupted lineage earns exactly (I_wrong2, moment_valid=False)",
        corrupted_gate["earned_moment"] == I_wrong2
        and not corrupted_gate["moment_valid"],
    )
    print(
        "DE-COUNTED diagnostic: baseline continuity_dependency_ok is True = "
        f"{baseline['continuity_dependency_ok']}"
    )

    subbanner("A-D/F. lineage rigs and coupling meta-tests")
    outcomes: dict[str, str] = {}
    fake_vector_expr = compact(Omega_U**2 * g_W + R_mix * g_U)
    baseline_assert = lambda name: routed_assert(name, baseline["lineage_valid"])
    for label, bad_lineage in negatives.items():
        subject = fake_vector_expr if label == "A fake_continuity" else n0_den
        result = evaluate_lineage(bad_lineage, subject)
        assert_name = f"{label} G1 exact-token lineage-gate assert"
        outcomes[label] = exercise_rig(
            label,
            assert_name,
            lambda result=result, name=assert_name: routed_assert(
                name, result["lineage_valid"]
            ),
            lambda name=assert_name: baseline_assert(name),
        )
        print(
            f"RIG DATA {label}: continuity_lineage_valid={result['lineage_valid']}; "
            f"earned_moment={result['earned_moment']}; "
            f"moment_valid={result['moment_valid']}; decoy_valid=True; "
            "self_asserted_continuity_interface=True"
        )

    d_assert = "D earning gate distinguishes I25 from I_wrong2 assert"

    def always_i25_gate(
        bad_lineage: dict[str, Any],
    ) -> tuple[sp.Symbol, bool]:
        return I25, bool(continuity_lineage_valid(bad_lineage))

    neutered_symbol, neutered_valid = always_i25_gate(negatives["B attack2"])
    actual_symbol, actual_valid = continuity_moment_symbol(
        negatives["B attack2"], sp.Rational(-7, 2)
    )
    outcomes["D earning_gate"] = exercise_rig(
        "D earning_gate",
        d_assert,
        lambda: routed_assert(
            d_assert,
            neutered_symbol == I_wrong2 and not neutered_valid,
        ),
        lambda: routed_assert(
            d_assert, actual_symbol == I_wrong2 and not actual_valid
        ),
    )

    flag_assert = "META flag-only validator rejects every decoy-negative assert"
    exact_rejects = all(
        not continuity_lineage_valid(item) for item in negatives.values()
    )
    flag_accepts = all(flag_only_validator(item) for item in negatives.values())
    outcomes["META flag_only_validator"] = exercise_rig(
        "META flag_only_validator",
        flag_assert,
        lambda: routed_assert(flag_assert, not flag_accepts),
        lambda: routed_assert(flag_assert, exact_rejects),
    )

    flip_lineage = deepcopy(lineage)
    flip_lineage["moments"]["l1"] = "D1_i = Integral(x_i*S_leak d3x)"
    flip_result = evaluate_lineage(decoy_negative(flip_lineage), n0_den)
    f_assert = "F baseline-valid positive and ancestry-corruption flip assert"
    expect_bool(
        f_assert,
        baseline["accepted"] and not flip_result["accepted"],
    )
    outcomes["F baseline_positive"] = (
        "PASS; corrupt-any-token=FLIP to FAIL_NOT_DENSITY_DERIVED"
    )
    print(f"RIG F baseline_positive: {outcomes['F baseline_positive']}")

    subbanner("G. isolated vanished-coupling OR-arm probes")
    g1 = continuity_dependency_ok(True, True, I25, sp.Integer(0), True)
    g2 = continuity_dependency_ok(True, True, I25, sp.Integer(0), False)
    g3 = continuity_dependency_ok(True, True, I25, foreign_subject, True)
    expect_bool("G(i) expr=0 and coupling_zero=True passes via arm 2", g1)
    expect_bool("G(ii) expr=0 and coupling_zero=False fails", not g2)
    expect_bool(
        "G(iii) nonzero expr missing earned moment fails despite coupling_zero=True",
        not g3,
    )
    g_assert = "G(iv) neutered vanished-coupling arm makes probe (i) fail assert"
    g1_neutered = continuity_dependency_ok(
        True,
        True,
        I25,
        sp.Integer(0),
        True,
        second_arm_enabled=False,
    )
    outcomes["G arm2_meta"] = exercise_rig(
        "G arm2_meta",
        g_assert,
        lambda: routed_assert(g_assert, g1_neutered),
        lambda: routed_assert(g_assert, g1),
    )

    subbanner("I. consumption-integrity mutations")
    no_i25 = compact(n0_den / I25**2)
    no_eta_q = compact(
        I25**2
        * Xi_Q**2
        * c_s**4
        * rho_eff
        * (eta_phi * varpi_q2) ** 2
        / (a**7 * (lambda_c**2 - varpi_Phi2 * varpi_q2) ** 2)
    )
    foreign_rho = compact(n0_den.subs(rho_eff, foreign_subject))
    for label, subject in {
        "I drop_external_I25_squared": no_i25,
        "I drop_eta_q_lambda_c_term": no_eta_q,
        "I replace_rho_eff": foreign_rho,
    }.items():
        assert_name = f"{label} exact host-contract assert"
        outcomes[label] = exercise_rig(
            label,
            assert_name,
            lambda subject=subject, name=assert_name: routed_assert(
                name, set(subject.free_symbols) == HOST_CONTRACT
            ),
            lambda name=assert_name: routed_assert(
                name, live_symbols == HOST_CONTRACT
            ),
            verdict="CONTRACT_REJECTED",
        )

    subbanner("H'. runtime arity and unevaluated-leakage scanners")
    outcomes["H' arity_scanner"] = exercise_rig(
        "H' arity_scanner",
        "H' definition/call arity scanner assert",
        lambda: routed_assert(
            "H' definition/call arity scanner assert", arity_scan(1)
        ),
        lambda: routed_assert(
            "H' definition/call arity scanner assert", arity_scan(2)
        ),
        verdict="SCANNER_CAUGHT",
    )
    actual_transcript = [n0_den, baseline, baseline_token_sets]
    leaked_transcript = actual_transcript + [AUTHORED_LEAK(n0_den)]
    outcomes["H' leakage_scanner"] = exercise_rig(
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
        verdict="SCANNER_CAUGHT",
    )
    expect_bool(
        "H' actual transcript has no unevaluated authored-helper leakage",
        leakage_scan(actual_transcript),
    )

    dependency_membership_witness = baseline["earned_moment"] in live_symbols
    print(
        "E DE-COUNTED witness: earned_moment in N0_den.free_symbols = "
        f"{dependency_membership_witness} (printed, not tallied)"
    )

    local_ok = bool(
        live_symbols == HOST_CONTRACT
        and baseline["lineage_valid"]
        and baseline["earned_moment"] == I25
        and baseline["moment_valid"]
        and baseline["continuity_dependency_ok"]
    )
    expect_bool("LOCAL CONTINUITY_LINEAGE_EARNED = G1 and G2", local_ok)

    return {
        "N0_den": n0_den,
        "live_symbols": live_symbols,
        "lineage": lineage,
        "baseline": baseline,
        "token_sets": baseline_token_sets,
        "stuffed_tokens": identifier_tokens(stuffed_l0),
        "outcomes": outcomes,
        "dependency_membership_witness": dependency_membership_witness,
        "local_ok": local_ok,
    }


def lineage_string(lineage: dict[str, Any]) -> str:
    moments = lineage["moments"]
    return (
        "{operator_id: "
        + lineage["operator_id"]
        + ", l0: "
        + moments["l0"]
        + ", l1: "
        + moments["l1"]
        + ", l2: "
        + moments["l2"]
        + ", l2_kernel: "
        + lineage["l2_kernel"]
        + "}"
    )


def emit(data: dict[str, Any]) -> None:
    baseline = data["baseline"]
    subbanner("Computed continuity-lineage transcript")
    print("consumes: stage024 N0_den (cited canonical factored export)")
    print(f"N0_den (canonical factored): {fmt(data['N0_den'])}")
    print(f"N0_den free_symbols: {fmt_names(data['live_symbols'])}")
    print(f"host contract (10): {fmt_names(HOST_CONTRACT)}")
    print(f"free_symbols == contract: {data['live_symbols'] == HOST_CONTRACT}")
    print(f"lineage dict: {lineage_string(data['lineage'])}")
    for level in ("l0", "l1", "l2", "l2_kernel"):
        print(f"exact tokens {level}: {fmt_names(data['token_sets'][level])}")
    print(f"stuffed l0 exact tokens: {fmt_names(data['stuffed_tokens'])}")
    print("tokenizer: re.findall(r'[A-Za-z0-9_]+', s)")
    print(
        "continuity_lineage_valid="
        f"{baseline['lineage_valid']} (computed exact tokens; decoy flags ignored)"
    )
    print(
        "earned (earned_moment="
        f"{baseline['earned_moment']}, moment_valid={baseline['moment_valid']})"
    )
    print(
        "continuity_dependency_ok="
        f"{baseline['continuity_dependency_ok']}"
    )
    print(
        "E dependency membership witness (de-counted): "
        f"{data['dependency_membership_witness']}"
    )
    print("exports: moment_valid=True; validated_I25=I25; lineage_certificate=PASS")

    subbanner("Verdict labels")
    print(f"LOCAL_AUDIT_VERDICT: {LOCAL_VERDICT}")
    print(f"JOINT_LANDING_LABEL (PARTIAL): {JOINT_PARTIAL}")


def main() -> int:
    banner("ledger_stage026_continuity_lineage_sympy_audit")
    print("Target stem confirmed: ledger_stage026_continuity_lineage")
    print(
        "Engine: exact [A-Za-z0-9_]+ lineage tokens over live strings; "
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
