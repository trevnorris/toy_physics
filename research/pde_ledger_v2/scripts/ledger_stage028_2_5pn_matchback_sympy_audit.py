#!/usr/bin/env python3
"""Ledger stage028 SymPy audit: calibrated 2.5PN match-back consistency.

Standalone, print-only, no arguments, and zero file I/O.  Stage027's
calibrated closure moments are restated locally as a provenance consumption;
the INV5 literal anchors are independent stage028-owned anti-rig constants.
"""

from __future__ import annotations

import inspect
from dataclasses import dataclass, replace
from typing import Any

import sympy as sp


PASS_COUNT = 0
FAIL_COUNT = 0
DIAGNOSTIC_PASS_COUNT = 0
DIAGNOSTIC_FAIL_COUNT = 0

LOCAL_VERDICT = "MATCHBACK_CONSISTENT"
R = sp.Rational
Z = sp.Integer


class AuditFailure(AssertionError):
    """A named stage028 audit assertion failed."""


class RigAssertion(AssertionError):
    """A computed mutation row missed its immutable caught-by oracle."""

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


def compact(expr: Any) -> sp.Expr:
    return sp.factor(sp.cancel(sp.together(sp.simplify(sp.sympify(expr)))))


def fmt(expr: Any) -> str:
    return sp.sstr(compact(expr))


def _record_pass(name: str) -> None:
    global PASS_COUNT
    PASS_COUNT += 1
    print(f"PASS  {name}")


def _record_fail(name: str, detail: str) -> None:
    global FAIL_COUNT
    FAIL_COUNT += 1
    print(f"FAIL  {name}: {detail}")


def _record_diagnostic(name: str, condition: bool, detail: str = "") -> None:
    global DIAGNOSTIC_PASS_COUNT, DIAGNOSTIC_FAIL_COUNT
    if condition:
        DIAGNOSTIC_PASS_COUNT += 1
        print(f"PASS  {name}")
        return
    DIAGNOSTIC_FAIL_COUNT += 1
    suffix = f": {detail}" if detail else ""
    print(f"FAIL  {name}{suffix}")
    raise AuditFailure(name)


def expect_zero(name: str, residual: sp.Expr | int) -> None:
    """Verdict-bearing zero diagnostic, deliberately not an EXIT-1 tooth tally."""
    clean = compact(residual)
    _record_diagnostic(name, clean == 0, f"residual = {sp.sstr(clean)}")


def expect_bool(name: str, condition: bool) -> None:
    expect_zero(name, Z(0) if bool(condition) else Z(1))


def routed_assert(name: str, condition: bool) -> None:
    if not bool(condition):
        raise RigAssertion(name)


def exercise_rig(
    mutation_name: str,
    actual: tuple[str, ...],
    expected: tuple[str, ...],
) -> None:
    """Route one runtime-computed caught-by row through its named tooth."""
    assertion_name = f"ROW {mutation_name} actual == EXPECTED_CAUGHT_BY"
    try:
        routed_assert(assertion_name, actual == expected)
    except RigAssertion as exc:
        _record_fail(exc.name, f"actual={actual!r}; expected={expected!r}")
        raise AuditFailure(exc.name) from exc
    _record_pass(assertion_name)


G, c_s, a, c = sp.symbols("G c_s a c", positive=True)

RESIDUAL_ORDER = (
    "INV1_moment_invariant",
    "INV2_pathA43_form_to_BT",
    "INV3_corpus_form_to_BT",
    "INV4_redundant_cross_form",
    "INV5_Kbar0_coeff_54_5",
    "INV5_Kbar2_coeff_6_5",
    "INV5_Kbar4_coeff_8_15",
    "INV5_BT_coeff_2_5",
    "INV5_pathA43_denominator_27",
    "INV5_corpus_factor_9",
    "INV5_exp_Kbar2_5_2",
    "INV5_exp_Kbar0_3_2",
)


@dataclass(frozen=True)
class MatchConfig:
    k0_coeff: sp.Expr
    k2_coeff: sp.Expr
    k4_coeff: sp.Expr
    bt_coeff: sp.Expr
    path_denominator: sp.Expr
    corpus_factor: sp.Expr
    exp_k2: sp.Expr
    exp_k0: sp.Expr
    k0_scale: sp.Expr = Z(1)
    k2_scale: sp.Expr = Z(1)
    k4_scale: sp.Expr = Z(1)
    bt_scale: sp.Expr = Z(1)


# Restated stage027 closure moments (PROVENANCE consumption, no runtime tie).
BASE = MatchConfig(
    k0_coeff=R(54, 5),
    k2_coeff=R(6, 5),
    k4_coeff=R(8, 15),
    bt_coeff=R(2, 5),
    path_denominator=Z(27),
    corpus_factor=Z(9),
    exp_k2=R(5, 2),
    exp_k0=R(3, 2),
)

# Stage028-owned INV5 RHS literals.  These are intentionally independent of
# BASE and are computed nowhere from MatchConfig or the moment expressions.
ANCHOR_K0 = R(54, 5)
ANCHOR_K2 = R(6, 5)
ANCHOR_K4 = R(8, 15)
ANCHOR_BT = R(2, 5)
ANCHOR_DENOMINATOR = Z(27)
ANCHOR_CORPUS_FACTOR = Z(9)
ANCHOR_EXP_K2 = R(5, 2)
ANCHOR_EXP_K0 = R(3, 2)

# Directive-v3 row order: coherent anti-rig, coupled control, nine singles.
MUTATIONS = (
    (
        "coherent_scale_all_moments_and_BT_x2",
        replace(BASE, k0_scale=Z(2), k2_scale=Z(2), k4_scale=Z(2), bt_scale=Z(2)),
    ),
    (
        "coupled_moments_x2_BT_fixed",
        replace(BASE, k0_scale=Z(2), k2_scale=Z(2), k4_scale=Z(2)),
    ),
    ("Kbar4_coeff_8_15_to_8_14", replace(BASE, k4_coeff=R(8, 14))),
    ("Kbar4_sign_flip", replace(BASE, k4_coeff=-R(8, 15))),
    ("Kbar2_coeff_6_5_to_7_5", replace(BASE, k2_coeff=R(7, 5))),
    ("Kbar0_coeff_54_5_to_55_5", replace(BASE, k0_coeff=R(55, 5))),
    ("pathA43_denominator_27_to_26", replace(BASE, path_denominator=Z(26))),
    ("corpus_factor_9_to_8", replace(BASE, corpus_factor=Z(8))),
    ("exp_Kbar2_5_2_to_3_2", replace(BASE, exp_k2=R(3, 2))),
    ("exp_Kbar0_3_2_to_1", replace(BASE, exp_k0=Z(1))),
    ("BT_coeff_2_5_to_3_5", replace(BASE, bt_coeff=R(3, 5))),
)

EXPECTED_CAUGHT_BY = {
    "coherent_scale_all_moments_and_BT_x2": (
        "INV5_Kbar0_coeff_54_5", "INV5_Kbar2_coeff_6_5",
        "INV5_Kbar4_coeff_8_15", "INV5_BT_coeff_2_5",
    ),
    "coupled_moments_x2_BT_fixed": (
        "INV2_pathA43_form_to_BT", "INV3_corpus_form_to_BT",
        "INV5_Kbar0_coeff_54_5", "INV5_Kbar2_coeff_6_5",
        "INV5_Kbar4_coeff_8_15",
    ),
    "Kbar4_coeff_8_15_to_8_14": (
        "INV1_moment_invariant", "INV5_Kbar4_coeff_8_15",
    ),
    "Kbar4_sign_flip": (
        "INV1_moment_invariant", "INV5_Kbar4_coeff_8_15",
    ),
    "Kbar2_coeff_6_5_to_7_5": (
        "INV1_moment_invariant", "INV3_corpus_form_to_BT",
        "INV4_redundant_cross_form", "INV5_Kbar2_coeff_6_5",
    ),
    "Kbar0_coeff_54_5_to_55_5": (
        "INV1_moment_invariant", "INV2_pathA43_form_to_BT",
        "INV3_corpus_form_to_BT", "INV4_redundant_cross_form",
        "INV5_Kbar0_coeff_54_5",
    ),
    "pathA43_denominator_27_to_26": (
        "INV2_pathA43_form_to_BT", "INV4_redundant_cross_form",
        "INV5_pathA43_denominator_27",
    ),
    "corpus_factor_9_to_8": (
        "INV3_corpus_form_to_BT", "INV4_redundant_cross_form",
        "INV5_corpus_factor_9",
    ),
    "exp_Kbar2_5_2_to_3_2": (
        "INV3_corpus_form_to_BT", "INV4_redundant_cross_form",
        "INV5_exp_Kbar2_5_2",
    ),
    "exp_Kbar0_3_2_to_1": (
        "INV3_corpus_form_to_BT", "INV4_redundant_cross_form",
        "INV5_exp_Kbar0_3_2",
    ),
    "BT_coeff_2_5_to_3_5": (
        "INV2_pathA43_form_to_BT", "INV3_corpus_form_to_BT",
        "INV5_BT_coeff_2_5",
    ),
}


def assert_no_float(label: str, expr: sp.Expr) -> None:
    floats = sp.sympify(expr).atoms(sp.Float)
    if floats:
        _record_fail("NO-FLOAT exact-rational guard", f"{label}: {sorted(map(str, floats))}")
        raise AuditFailure("NO-FLOAT exact-rational guard")


def build_residuals(cfg: MatchConfig) -> dict[str, sp.Expr]:
    k0 = cfg.k0_scale * cfg.k0_coeff * G * c_s**5 / (a**5 * c**5)
    k2 = cfg.k2_scale * cfg.k2_coeff * G * c_s**3 / (a**3 * c**5)
    k4 = cfg.k4_scale * cfg.k4_coeff * G * c_s / (a * c**5)
    bt = cfg.bt_scale * cfg.bt_coeff * G / c**5
    path_form = k0 * a**5 / (cfg.path_denominator * c_s**5)
    corpus_form = cfg.corpus_factor * k2**cfg.exp_k2 / k0**cfg.exp_k0

    raw = {
        "INV1_moment_invariant": k4 * k0 - Z(4) * k2**2,
        "INV2_pathA43_form_to_BT": path_form - bt,
        "INV3_corpus_form_to_BT": corpus_form - bt,
        "INV4_redundant_cross_form": path_form - corpus_form,
        "INV5_Kbar0_coeff_54_5": k0 * a**5 * c**5 / (G * c_s**5) - ANCHOR_K0,
        "INV5_Kbar2_coeff_6_5": k2 * a**3 * c**5 / (G * c_s**3) - ANCHOR_K2,
        "INV5_Kbar4_coeff_8_15": k4 * a * c**5 / (G * c_s) - ANCHOR_K4,
        "INV5_BT_coeff_2_5": bt * c**5 / G - ANCHOR_BT,
        "INV5_pathA43_denominator_27": cfg.path_denominator - ANCHOR_DENOMINATOR,
        "INV5_corpus_factor_9": cfg.corpus_factor - ANCHOR_CORPUS_FACTOR,
        "INV5_exp_Kbar2_5_2": cfg.exp_k2 - ANCHOR_EXP_K2,
        "INV5_exp_Kbar0_3_2": cfg.exp_k0 - ANCHOR_EXP_K0,
    }
    reduced: dict[str, sp.Expr] = {}
    for label in RESIDUAL_ORDER:
        assert_no_float(f"{label} raw", raw[label])
        reduced[label] = compact(raw[label])
        assert_no_float(f"{label} reduced", reduced[label])
    return reduced


def fired_labels(residuals: dict[str, sp.Expr]) -> tuple[str, ...]:
    return tuple(label for label in RESIDUAL_ORDER if residuals[label] != 0)


def vector_text(residuals: dict[str, sp.Expr]) -> str:
    return "[" + ", ".join(
        f"{label}={sp.sstr(residuals[label])}" for label in RESIDUAL_ORDER
    ) + "]"


def zero_test_text(residuals: dict[str, sp.Expr]) -> str:
    return "[" + ", ".join(
        f"{label}={'True' if residuals[label] == 0 else 'False'}"
        for label in RESIDUAL_ORDER
    ) + "]"


def arity_scan(function: Any, expected: int) -> bool:
    positional = [
        parameter
        for parameter in inspect.signature(function).parameters.values()
        if parameter.kind
        in (parameter.POSITIONAL_ONLY, parameter.POSITIONAL_OR_KEYWORD)
    ]
    return len(positional) == expected


AUTHORED_LEAK = sp.Function("build_residuals")


def leakage_scan(objects: list[Any]) -> bool:
    return not any(
        isinstance(obj, sp.Basic) and obj.has(AUTHORED_LEAK) for obj in objects
    )


def print_provenance() -> None:
    subbanner("PROVENANCE and A3 boundary")
    print("Kbar0,Kbar2,Kbar4 are CALIBRATED closure inputs restated from upstream stage027; moments are not derived here.")
    print("External calibrated literals: 2/5, 54/5, and G; G=GENUINE_BLOCKED; Gamma5/G is not derived here.")
    print("The pathA_43 denominator 27 is the upstream-EARNED density-Hankel fingerprint imported from stage018.")
    print("Corpus form 9*Kbar2^(5/2)/Kbar0^(3/2) is imported from 4d_2_5pn.tex:469, not re-derived.")
    print("Runtime isolation: no peer import/subprocess, no report/source/note/_scratch reads, no writes; zero file I/O.")
    print("A3 BOUNDARY: INV1/INV2 are SHARED re-expressions of stage027 closure (R46).")
    print("A3 BOUNDARY: INV3/INV4/INV5 and the coherent-rescale matrix are stage028 content.")
    print("A3 BOUNDARY: INV4 == INV2-INV3 identically; retained as a redundant diagnostic, not an independent tooth.")


def baseline_moments() -> tuple[sp.Expr, sp.Expr, sp.Expr, sp.Expr]:
    k0 = BASE.k0_coeff * G * c_s**5 / (a**5 * c**5)
    k2 = BASE.k2_coeff * G * c_s**3 / (a**3 * c**5)
    k4 = BASE.k4_coeff * G * c_s / (a * c**5)
    gamma5 = compact(k0 * a**5 / (BASE.path_denominator * c_s**5))
    return tuple(map(compact, (k0, k2, k4, gamma5)))  # type: ignore[return-value]


def run_audit() -> None:
    print_provenance()
    k0, k2, k4, gamma5 = baseline_moments()
    print("RESTATED stage027 closure_overlay moments:")
    print(f"  Kbar0={fmt(k0)}")
    print(f"  Kbar2={fmt(k2)}")
    print(f"  Kbar4={fmt(k4)}")
    print(f"  Gamma5=Kbar0*a^5/(27*c_s^5)={fmt(gamma5)}")
    print("  A3_FIDELITY: restated literals match stage027 closure_overlay (authoring-time comparison).")

    subbanner("Baseline exact-rational residuals and zero tests")
    baseline = build_residuals(BASE)
    print("BASELINE RESIDUAL VECTOR:")
    print("  " + vector_text(baseline))
    print("BASELINE ZERO TESTS:")
    print("  " + zero_test_text(baseline))
    for label in RESIDUAL_ORDER:
        expect_zero(f"ZERO-TEST {label} == 0", baseline[label])
    print("BASELINE ALL-ZERO: PASS (aggregate positive; DE-COUNTED)")

    subbanner("Runtime-computed mutation residuals")
    caught_matrix: dict[str, tuple[str, ...]] = {}
    transcript_objects: list[Any] = list(baseline.values())
    for name, cfg in MUTATIONS:
        residuals = build_residuals(cfg)
        transcript_objects.extend(residuals.values())
        caught = fired_labels(residuals)
        caught_matrix[name] = caught
        print(f"MUTATION {name}:")
        print("  residuals: " + vector_text(residuals))
        print("  zero_tests: " + zero_test_text(residuals))
        print("  caught_by: [" + ", ".join(caught) + "]")

    subbanner("Immutable EXPECTED caught-by matrix")
    for name, _cfg in MUTATIONS:
        expected = EXPECTED_CAUGHT_BY[name]
        actual = caught_matrix[name]
        print(
            f"ROW {name}: expected=[{', '.join(expected)}] "
            f"actual=[{', '.join(actual)}]"
        )
        exercise_rig(name, actual, expected)
    print("MUTATION PROBE: PASS (11 runtime-computed rows)")

    expect_bool(
        "ARITY build_residuals definition/call contract is 1",
        arity_scan(build_residuals, 1),
    )
    expect_bool(
        "LEAKAGE transcript has no unevaluated authored helper",
        leakage_scan(transcript_objects),
    )
    _record_pass("NO-FLOAT exact-rational guard over raw and reduced residuals")
    print("NO-FLOAT GUARD: PASS")

    subbanner("LOCAL calibrated-consistency verdict")
    print(f"LOCAL_AUDIT_VERDICT: {LOCAL_VERDICT}")
    print("scope: CALIBRATED consistency over closure moments; G=GENUINE_BLOCKED.")
    print("scope: NOT a first-principles Gamma5/G derivation; full 1PN->4PN from-throat is SIM-DEFERRED Gate 6.")


def main() -> int:
    banner("ledger_stage028_2_5pn_matchback_sympy_audit")
    print("Target stem confirmed: ledger_stage028_2_5pn_matchback")
    print("Engine: exact SymPy residual algebra; no floats; zero file I/O.")
    run_audit()
    return 0


if __name__ == "__main__":
    try:
        exit_code = main()
    except Exception as exc:
        if not isinstance(exc, AuditFailure):
            _record_fail("UNCAUGHT exception", repr(exc))
        banner("Tallies")
        print(
            f"TALLY sympy: {PASS_COUNT} pass + {FAIL_COUNT} fail = "
            f"{PASS_COUNT + FAIL_COUNT} independent EXIT-1 teeth; "
            f"diagnostics={DIAGNOSTIC_PASS_COUNT} pass + "
            f"{DIAGNOSTIC_FAIL_COUNT} fail"
        )
        print("OVERALL FAIL")
        raise SystemExit(1) from exc

    banner("Tallies")
    print(
        f"TALLY sympy: {PASS_COUNT} pass + {FAIL_COUNT} fail = "
        f"{PASS_COUNT + FAIL_COUNT} independent EXIT-1 teeth; "
        f"diagnostics={DIAGNOSTIC_PASS_COUNT} pass + "
        f"{DIAGNOSTIC_FAIL_COUNT} fail"
    )
    if exit_code == 0 and FAIL_COUNT == 0 and DIAGNOSTIC_FAIL_COUNT == 0:
        print("OVERALL PASS")
        raise SystemExit(0)
    print("OVERALL FAIL")
    raise SystemExit(1)
