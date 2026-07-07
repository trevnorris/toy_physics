#!/usr/bin/env python3
"""A3 2.5PN match-back verification, SymPy engine.

Self-contained by design: this script restates the calibrated literals locally,
does not read or write repo artifacts, and does not consume the Mathematica
peer in any way.
"""

from __future__ import annotations

from dataclasses import dataclass, replace

import sympy as sp


R = sp.Rational
Z = sp.Integer

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

MUTATIONS = (
    ("Kbar4_coeff_8_15_to_8_14", replace(BASE, k4_coeff=R(8, 14))),
    ("Kbar4_sign_flip", replace(BASE, k4_coeff=-R(8, 15))),
    ("Kbar2_coeff_6_5_to_7_5", replace(BASE, k2_coeff=R(7, 5))),
    ("Kbar0_coeff_54_5_to_55_5", replace(BASE, k0_coeff=R(55, 5))),
    ("pathA43_denominator_27_to_26", replace(BASE, path_denominator=Z(26))),
    ("corpus_factor_9_to_8", replace(BASE, corpus_factor=Z(8))),
    ("exp_Kbar2_5_2_to_3_2", replace(BASE, exp_k2=R(3, 2))),
    ("exp_Kbar0_3_2_to_1", replace(BASE, exp_k0=Z(1))),
    ("BT_coeff_2_5_to_3_5", replace(BASE, bt_coeff=R(3, 5))),
    (
        "coherent_scale_all_moments_and_BT_x2",
        replace(BASE, k0_scale=Z(2), k2_scale=Z(2), k4_scale=Z(2), bt_scale=Z(2)),
    ),
    (
        "coupled_moments_x2_BT_fixed",
        replace(BASE, k0_scale=Z(2), k2_scale=Z(2), k4_scale=Z(2)),
    ),
)

EXPECTED_CAUGHT_BY = {
    "Kbar4_coeff_8_15_to_8_14": (
        "INV1_moment_invariant",
        "INV5_Kbar4_coeff_8_15",
    ),
    "Kbar4_sign_flip": (
        "INV1_moment_invariant",
        "INV5_Kbar4_coeff_8_15",
    ),
    "Kbar2_coeff_6_5_to_7_5": (
        "INV1_moment_invariant",
        "INV3_corpus_form_to_BT",
        "INV4_redundant_cross_form",
        "INV5_Kbar2_coeff_6_5",
    ),
    "Kbar0_coeff_54_5_to_55_5": (
        "INV1_moment_invariant",
        "INV2_pathA43_form_to_BT",
        "INV3_corpus_form_to_BT",
        "INV4_redundant_cross_form",
        "INV5_Kbar0_coeff_54_5",
    ),
    "pathA43_denominator_27_to_26": (
        "INV2_pathA43_form_to_BT",
        "INV4_redundant_cross_form",
        "INV5_pathA43_denominator_27",
    ),
    "corpus_factor_9_to_8": (
        "INV3_corpus_form_to_BT",
        "INV4_redundant_cross_form",
        "INV5_corpus_factor_9",
    ),
    "exp_Kbar2_5_2_to_3_2": (
        "INV3_corpus_form_to_BT",
        "INV4_redundant_cross_form",
        "INV5_exp_Kbar2_5_2",
    ),
    "exp_Kbar0_3_2_to_1": (
        "INV3_corpus_form_to_BT",
        "INV4_redundant_cross_form",
        "INV5_exp_Kbar0_3_2",
    ),
    "BT_coeff_2_5_to_3_5": (
        "INV2_pathA43_form_to_BT",
        "INV3_corpus_form_to_BT",
        "INV5_BT_coeff_2_5",
    ),
    "coherent_scale_all_moments_and_BT_x2": (
        "INV5_Kbar0_coeff_54_5",
        "INV5_Kbar2_coeff_6_5",
        "INV5_Kbar4_coeff_8_15",
        "INV5_BT_coeff_2_5",
    ),
    "coupled_moments_x2_BT_fixed": (
        "INV2_pathA43_form_to_BT",
        "INV3_corpus_form_to_BT",
        "INV5_Kbar0_coeff_54_5",
        "INV5_Kbar2_coeff_6_5",
        "INV5_Kbar4_coeff_8_15",
    ),
}


def compact(expr: sp.Expr) -> sp.Expr:
    return sp.factor(sp.cancel(sp.simplify(sp.sympify(expr))))


def assert_no_float(label: str, expr: sp.Expr) -> None:
    floats = sp.sympify(expr).atoms(sp.Float)
    if floats:
        raise AssertionError(f"{label} contains Float atoms: {sorted(map(str, floats))}")


def build_residuals(cfg: MatchConfig) -> dict[str, sp.Expr]:
    k0 = sp.sympify(cfg.k0_scale) * cfg.k0_coeff * G * c_s**5 / (a**5 * c**5)
    k2 = sp.sympify(cfg.k2_scale) * cfg.k2_coeff * G * c_s**3 / (a**3 * c**5)
    k4 = sp.sympify(cfg.k4_scale) * cfg.k4_coeff * G * c_s / (a * c**5)
    bt = sp.sympify(cfg.bt_scale) * cfg.bt_coeff * G / c**5

    path_form = k0 * a**5 / (cfg.path_denominator * c_s**5)
    corpus_form = cfg.corpus_factor * k2**cfg.exp_k2 / k0**cfg.exp_k0

    raw = {
        "INV1_moment_invariant": k4 * k0 - Z(4) * k2**2,
        "INV2_pathA43_form_to_BT": path_form - bt,
        "INV3_corpus_form_to_BT": corpus_form - bt,
        "INV4_redundant_cross_form": path_form - corpus_form,
        "INV5_Kbar0_coeff_54_5": k0 * a**5 * c**5 / (G * c_s**5) - R(54, 5),
        "INV5_Kbar2_coeff_6_5": k2 * a**3 * c**5 / (G * c_s**3) - R(6, 5),
        "INV5_Kbar4_coeff_8_15": k4 * a * c**5 / (G * c_s) - R(8, 15),
        "INV5_BT_coeff_2_5": bt * c**5 / G - R(2, 5),
        "INV5_pathA43_denominator_27": cfg.path_denominator - Z(27),
        "INV5_corpus_factor_9": cfg.corpus_factor - Z(9),
        "INV5_exp_Kbar2_5_2": cfg.exp_k2 - R(5, 2),
        "INV5_exp_Kbar0_3_2": cfg.exp_k0 - R(3, 2),
    }

    reduced: dict[str, sp.Expr] = {}
    for label in RESIDUAL_ORDER:
        raw_expr = sp.sympify(raw[label])
        assert_no_float(f"{label} raw", raw_expr)
        reduced_expr = compact(raw_expr)
        assert_no_float(f"{label} reduced", reduced_expr)
        reduced[label] = reduced_expr
    return reduced


def fired_labels(residuals: dict[str, sp.Expr]) -> tuple[str, ...]:
    return tuple(label for label in RESIDUAL_ORDER if residuals[label] != 0)


def vector_text(residuals: dict[str, sp.Expr]) -> str:
    return "[" + ", ".join(f"{label}={sp.sstr(residuals[label])}" for label in RESIDUAL_ORDER) + "]"


def print_provenance() -> None:
    print("PROVENANCE:")
    print("  Kbar0,Kbar2,Kbar4 are CALIBRATED closure inputs; moments are not derived here.")
    print("  External calibrated literals: 2/5, 54/5, and G; G=GENUINE_BLOCKED.")
    print("  The pathA_43 denominator 27 is upstream-EARNED density-Hankel fingerprint, imported here.")
    print("  Corpus form 9*Kbar2^(5/2)/Kbar0^(3/2) is imported from 4d_2_5pn.tex, not re-derived.")
    print("  Runtime isolation: no peer-engine import/subprocess, no report/source/note/_scratch reads, no writes.")


def main() -> int:
    print_provenance()

    baseline = build_residuals(BASE)
    print("BASELINE RESIDUAL VECTOR (SymPy):")
    print("  " + vector_text(baseline))
    if any(expr != 0 for expr in baseline.values()):
        raise AssertionError("baseline residual vector is not all zero")

    print("MUTATED RESIDUAL VECTORS (SymPy):")
    caught_matrix: dict[str, tuple[str, ...]] = {}
    for name, cfg in MUTATIONS:
        residuals = build_residuals(cfg)
        caught = fired_labels(residuals)
        caught_matrix[name] = caught
        print(f"  {name}:")
        print("    residuals: " + vector_text(residuals))
        print("    caught_by: [" + ", ".join(caught) + "]")

    print("EXPECTED CAUGHT-BY MATRIX (SymPy):")
    for name, _cfg in MUTATIONS:
        expected = EXPECTED_CAUGHT_BY[name]
        actual = caught_matrix[name]
        print(f"  {name}: expected=[{', '.join(expected)}] actual=[{', '.join(actual)}]")
        if actual != expected:
            raise AssertionError(f"caught-by mismatch for {name}")
        if not actual:
            raise AssertionError(f"mutation was not caught: {name}")

    print("NO-FLOAT GUARD: PASS")
    print("BASELINE ALL-ZERO: PASS")
    print("MUTATION PROBE: PASS")
    print("SYMPY_MATCHBACK: PASS")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except AssertionError as exc:
        print(f"FAIL: {exc}")
        raise SystemExit(1)
