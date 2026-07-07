#!/usr/bin/env python3
"""Ledger stage001 SymPy audit: solid angles and second moments.

Print-only, standalone, no arguments.  The derived side is recomputed from the
parametrized S^2 and S^3 surface integrals used by pathA_21c.
"""

from __future__ import annotations

from dataclasses import dataclass

import sympy as sp


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


def compact(expr: sp.Expr) -> str:
    return sp.sstr(sp.factor(sp.cancel(sp.simplify(expr))))


def assert_no_float(name: str, expr: sp.Expr) -> None:
    floats = expr.atoms(sp.Float)
    if floats:
        raise AuditFailure(f"{name}: Float atom(s) found in exact audit expression: {floats}")


def expect_zero(name: str, residual: sp.Expr) -> None:
    assert_no_float(name, residual)
    clean = sp.simplify(residual)
    assert_no_float(name, clean)
    if clean == 0:
        print(f"PASS  {name}")
        return

    print(f"FAIL  {name}: residual = {compact(clean)}")
    raise AuditFailure(f"{name} residual was not zero")


@dataclass(frozen=True)
class SurfaceIntegrals:
    omega2: sp.Expr
    s2_n1_sq: sp.Expr
    s2_n1_n2: sp.Expr
    omega3: sp.Expr
    s3_n1_sq: sp.Expr
    s3_n1_n2: sp.Expr


def derive_surface_integrals() -> SurfaceIntegrals:
    theta, phi = sp.symbols("theta phi", real=True)
    chi = sp.symbols("chi", real=True)

    s2_n1 = sp.cos(theta)
    s2_n2 = sp.sin(theta) * sp.cos(phi)
    s2_measure = sp.sin(theta)
    s2_ranges = ((phi, 0, 2 * sp.pi), (theta, 0, sp.pi))
    omega2 = sp.integrate(s2_measure, *s2_ranges)

    s3_n1 = sp.cos(chi)
    s3_n2 = sp.sin(chi) * sp.cos(theta)
    s3_measure = sp.sin(chi) ** 2 * sp.sin(theta)
    s3_ranges = ((phi, 0, 2 * sp.pi), (theta, 0, sp.pi), (chi, 0, sp.pi))
    omega3 = sp.integrate(s3_measure, *s3_ranges)

    return SurfaceIntegrals(
        omega2=sp.simplify(omega2),
        s2_n1_sq=sp.simplify(sp.integrate(s2_n1**2 * s2_measure, *s2_ranges) / omega2),
        s2_n1_n2=sp.simplify(sp.integrate(s2_n1 * s2_n2 * s2_measure, *s2_ranges) / omega2),
        omega3=sp.simplify(omega3),
        s3_n1_sq=sp.simplify(sp.integrate(s3_n1**2 * s3_measure, *s3_ranges) / omega3),
        s3_n1_n2=sp.simplify(sp.integrate(s3_n1 * s3_n2 * s3_measure, *s3_ranges) / omega3),
    )


BASELINE_TARGETS: dict[str, sp.Expr] = {
    "Omega_2 from S^2 surface integral": 4 * sp.pi,
    "<n1^2>_S^2 normalized second moment": sp.Rational(1, 3),
    "<n1 n2>_S^2 normalized cross moment": sp.Integer(0),
    "Omega_3 from S^3 surface integral": 2 * sp.pi**2,
    "<n1^2>_S^3 normalized second moment": sp.Rational(1, 4),
    "<n1 n2>_S^3 normalized cross moment": sp.Integer(0),
}


def actual_claims(surface: SurfaceIntegrals) -> dict[str, sp.Expr]:
    return {
        "Omega_2 from S^2 surface integral": surface.omega2,
        "<n1^2>_S^2 normalized second moment": surface.s2_n1_sq,
        "<n1 n2>_S^2 normalized cross moment": surface.s2_n1_n2,
        "Omega_3 from S^3 surface integral": surface.omega3,
        "<n1^2>_S^3 normalized second moment": surface.s3_n1_sq,
        "<n1 n2>_S^3 normalized cross moment": surface.s3_n1_n2,
    }


def residuals(
    surface: SurfaceIntegrals,
    target_overrides: dict[str, sp.Expr] | None = None,
) -> dict[str, sp.Expr]:
    targets = dict(BASELINE_TARGETS)
    if target_overrides:
        targets.update(target_overrides)

    out: dict[str, sp.Expr] = {}
    for name, actual in actual_claims(surface).items():
        target = targets[name]
        assert_no_float(f"{name} actual", actual)
        assert_no_float(f"{name} target", target)
        out[name] = sp.simplify(actual - target)
    return out


@dataclass(frozen=True)
class MutationCase:
    label: str
    target_overrides: dict[str, sp.Expr]
    required_failure: str


def mutation_cases() -> list[MutationCase]:
    return [
        MutationCase(
            label="perturb Omega_2 target 4*pi -> 2*pi",
            target_overrides={"Omega_2 from S^2 surface integral": 2 * sp.pi},
            required_failure="Omega_2 from S^2 surface integral",
        ),
        MutationCase(
            label="perturb Omega_3 target 2*pi^2 -> 4*pi^2",
            target_overrides={"Omega_3 from S^3 surface integral": 4 * sp.pi**2},
            required_failure="Omega_3 from S^3 surface integral",
        ),
        MutationCase(
            label="perturb S^2 <n1^2> target 1/3 -> 1/2",
            target_overrides={"<n1^2>_S^2 normalized second moment": sp.Rational(1, 2)},
            required_failure="<n1^2>_S^2 normalized second moment",
        ),
        MutationCase(
            label="perturb S^3 <n1^2> target 1/4 -> 1/2",
            target_overrides={"<n1^2>_S^3 normalized second moment": sp.Rational(1, 2)},
            required_failure="<n1^2>_S^3 normalized second moment",
        ),
        MutationCase(
            label="cross-moment control: perturb S^2 <n1 n2> target 0 -> 1/5",
            target_overrides={"<n1 n2>_S^2 normalized cross moment": sp.Rational(1, 5)},
            required_failure="<n1 n2>_S^2 normalized cross moment",
        ),
    ]


def run_baseline(surface: SurfaceIntegrals) -> None:
    subbanner("Baseline exact residuals")
    for name, residual in residuals(surface).items():
        expect_zero(name, residual)


def assert_mutation_fails(surface: SurfaceIntegrals, case: MutationCase) -> None:
    mutated = residuals(surface, case.target_overrides)
    try:
        expect_zero(case.required_failure, mutated[case.required_failure])
    except AuditFailure:
        print(f"PASS  mutation probe: {case.label} produced the required FAIL")
        return

    print(f"FAIL  mutation probe: {case.label} survived")
    raise AuditFailure(f"mutation probe survived: {case.label}")


def run_mutation_probe(surface: SurfaceIntegrals) -> None:
    subbanner("Able-to-fail mutation probe")
    for case in mutation_cases():
        assert_mutation_fails(surface, case)


def main() -> None:
    banner("ledger_stage001_solid_angle_second_moment_primitives SymPy audit")
    surface = derive_surface_integrals()
    run_baseline(surface)
    run_mutation_probe(surface)
    print("")
    print("OVERALL PASS: SymPy derived all stage001 geometry primitives exactly")


if __name__ == "__main__":
    main()
