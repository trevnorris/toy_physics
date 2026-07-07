#!/usr/bin/env python3
"""Ledger stage002 SymPy audit: matter-stress force assembly.

Print-only, standalone, no arguments.  This stage consumes the stage001
geometry primitives as symbolic inputs:

  Omega_2 = 4*pi, Omega_3 = 2*pi^2, mu_d = <n1^2> = 1/d.

It does not re-integrate those surface primitives.  The force factors are
assembled from the Noether matter stress, Bernoulli pressure cross term, and
Gauss drain substitution.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import sympy as sp


Dim = tuple[sp.Rational, sp.Rational, sp.Rational]


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
    clean = sp.sympify(expr)
    floats = clean.atoms(sp.Float)
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


def q(value: int | str) -> sp.Rational:
    return sp.Rational(value)


def dim(l_power: int, t_power: int, m_power: int) -> Dim:
    return (q(l_power), q(t_power), q(m_power))


def dim_add(*items: Dim) -> Dim:
    out = [q(0), q(0), q(0)]
    for item in items:
        for i, power in enumerate(item):
            out[i] += power
    return (out[0], out[1], out[2])


def dim_neg(item: Dim) -> Dim:
    return (-item[0], -item[1], -item[2])


def dim_sub(left: Dim, right: Dim) -> Dim:
    return dim_add(left, dim_neg(right))


def dim_scale(scale: int | sp.Rational, item: Dim) -> Dim:
    factor = q(scale) if isinstance(scale, int) else scale
    return (factor * item[0], factor * item[1], factor * item[2])


def dimension_residual(terms: dict[str, Dim], expected: Dim) -> sp.Expr:
    residual = sp.Integer(0)
    for actual in terms.values():
        residual += sum((component - want) ** 2 for component, want in zip(actual, expected))
    return sp.simplify(residual)


@dataclass(frozen=True)
class Symbols:
    m_gnls: sp.Symbol
    n3: sp.Symbol
    rho4: sp.Symbol
    q1: sp.Symbol
    q2: sp.Symbol
    r12: sp.Symbol
    radius4: sp.Symbol
    v1: sp.Symbol
    k_eos: sp.Symbol
    rho: sp.Symbol
    vdot: sp.Symbol
    hbar: sp.Symbol
    mass_q: sp.Symbol
    x: sp.Symbol


def make_symbols() -> Symbols:
    return Symbols(
        m_gnls=sp.Symbol("m_GNLS", nonzero=True),
        n3=sp.Symbol("N", nonzero=True),
        rho4=sp.Symbol("rho4", nonzero=True),
        q1=sp.Symbol("Q1", nonzero=True),
        q2=sp.Symbol("Q2", nonzero=True),
        r12=sp.Symbol("r12", nonzero=True),
        radius4=sp.Symbol("R", nonzero=True),
        v1=sp.Symbol("v1"),
        k_eos=sp.Symbol("K", nonzero=True),
        rho=sp.Symbol("rho", nonzero=True),
        vdot=sp.Symbol("vdot"),
        hbar=sp.Symbol("hbar", nonzero=True),
        mass_q=sp.Symbol("m", nonzero=True),
        x=sp.Symbol("x", real=True),
    )


@dataclass(frozen=True)
class Stage001SurfaceInput:
    label: str
    ambient_dim: int
    omega_symbol: sp.Symbol
    omega_value: sp.Expr
    mu_symbol: sp.Symbol
    mu_value: sp.Expr
    shell_radius: sp.Symbol

    @property
    def stage001_subs(self) -> dict[sp.Symbol, sp.Expr]:
        return {
            self.omega_symbol: self.omega_value,
            self.mu_symbol: self.mu_value,
        }


def stage001_surfaces() -> tuple[Stage001SurfaceInput, Stage001SurfaceInput]:
    # ledger_stage001 owns the surface integrations.  This stage consumes the
    # resulting Omega_d and mu_d=<n1^2>=1/d values as live symbolic inputs.
    return (
        Stage001SurfaceInput(
            label="reduced-3D S2 boundary",
            ambient_dim=3,
            omega_symbol=sp.Symbol("Omega_2", positive=True),
            omega_value=4 * sp.pi,
            mu_symbol=sp.Symbol("mu_3"),
            mu_value=sp.Rational(1, 3),
            shell_radius=sp.Symbol("a", positive=True),
        ),
        Stage001SurfaceInput(
            label="bulk-4D S3 boundary",
            ambient_dim=4,
            omega_symbol=sp.Symbol("Omega_3", positive=True),
            omega_value=2 * sp.pi**2,
            mu_symbol=sp.Symbol("mu_4"),
            mu_value=sp.Rational(1, 4),
            shell_radius=sp.Symbol("A", positive=True),
        ),
    )


@dataclass(frozen=True)
class DerivationConfig:
    drain_sign: sp.Integer = sp.Integer(-1)
    force_flux_sign: sp.Integer = sp.Integer(-1)
    convective_moment3_override: sp.Expr | None = None


@dataclass(frozen=True)
class FluxSet:
    convective3: sp.Expr
    pressure3: sp.Expr
    total3: sp.Expr
    force3: sp.Expr
    total4: sp.Expr
    force4: sp.Expr


def derive_dimensions() -> dict[str, Dim]:
    l_unit = dim(1, 0, 0)
    t_unit = dim(0, 1, 0)
    m_unit = dim(0, 0, 1)

    action = dim_add(dim_scale(2, l_unit), dim_neg(t_unit), m_unit)
    energy = dim_sub(action, t_unit)
    force = dim_sub(energy, l_unit)
    velocity = dim_sub(l_unit, t_unit)
    rho3 = dim_scale(-3, l_unit)
    rho4 = dim_scale(-4, l_unit)
    number_rate = dim_neg(t_unit)
    q3_vol = dim_sub(number_rate, rho3)
    q4_vol = dim_sub(number_rate, rho4)
    mass_density3 = dim_add(m_unit, rho3)
    mass_density4 = dim_add(m_unit, rho4)
    stress3 = dim_sub(force, dim_scale(2, l_unit))
    stress4 = dim_sub(force, dim_scale(3, l_unit))
    force_density3 = dim_sub(force, dim_scale(3, l_unit))
    force_density4 = dim_sub(force, dim_scale(4, l_unit))
    energy_gradient = dim_sub(energy, l_unit)
    quantum_prefactor = dim_add(dim_scale(2, action), dim_neg(m_unit), dim_scale(-2, l_unit))

    return {
        "L": l_unit,
        "T": t_unit,
        "M": m_unit,
        "action": action,
        "energy": energy,
        "force": force,
        "velocity": velocity,
        "rho3": rho3,
        "rho4": rho4,
        "numberRate": number_rate,
        "q3Vol": q3_vol,
        "q4Vol": q4_vol,
        "massDensity3": mass_density3,
        "massDensity4": mass_density4,
        "stress3": stress3,
        "stress4": stress4,
        "forceDensity3": force_density3,
        "forceDensity4": force_density4,
        "energyGradient": energy_gradient,
        "quantumPrefactor": quantum_prefactor,
    }


def dim_claims(expected_overrides: dict[str, Dim] | None = None) -> dict[str, sp.Expr]:
    expected_overrides = expected_overrides or {}
    d = derive_dimensions()
    l_unit = d["L"]
    t_unit = d["T"]

    checks: list[tuple[str, dict[str, Dim], Dim]] = [
        (
            "reduced-3D momentum-balance terms",
            {
                "partial_t(m*N*v_i)": dim_sub(dim_add(d["massDensity3"], d["velocity"]), t_unit),
                "partial_j Pi_ij": dim_sub(d["stress3"], l_unit),
                "V_conf body force": dim_add(d["rho3"], d["energyGradient"]),
            },
            dim(-2, -2, 1),
        ),
        (
            "bulk-4D momentum-balance terms",
            {
                "partial_t(m*rho*v_i)": dim_sub(dim_add(d["massDensity4"], d["velocity"]), t_unit),
                "partial_J Pi_iJ": dim_sub(d["stress4"], l_unit),
                "V_conf body force": dim_add(d["rho4"], d["energyGradient"]),
            },
            dim(-3, -2, 1),
        ),
        (
            "reduced-3D Noether stress representatives",
            {
                "convective": dim_add(d["massDensity3"], dim_scale(2, d["velocity"])),
                "pressure": dim_add(dim_sub(d["energy"], dim_scale(4, l_unit)), l_unit),
                "quantum": dim_add(d["quantumPrefactor"], d["rho3"]),
            },
            dim(-1, -2, 1),
        ),
        (
            "Euler force-per-volume identity terms in reduced 3D",
            {
                "m*N*acceleration": dim_sub(dim_add(d["massDensity3"], d["velocity"]), t_unit),
                "N*grad h": dim_add(d["rho3"], d["energyGradient"]),
                "N*grad Q": dim_add(d["rho3"], d["energyGradient"]),
                "N*grad V_conf": dim_add(d["rho3"], d["energyGradient"]),
                "N*q(E+vB)": dim_add(d["rho3"], d["energyGradient"]),
            },
            dim(-2, -2, 1),
        ),
        (
            "reduced-3D surface traction integral gives force",
            {"stress3 + 2L": dim_add(d["stress3"], dim_scale(2, l_unit))},
            dim(1, -2, 1),
        ),
        (
            "bulk-4D surface traction integral gives force",
            {"stress4 + 3L": dim_add(d["stress4"], dim_scale(3, l_unit))},
            dim(1, -2, 1),
        ),
        (
            "reduced-3D Noether force coefficient m*N*Q1*Q2",
            {"massDensity3 + 2*q3Vol": dim_add(d["massDensity3"], dim_scale(2, d["q3Vol"]))},
            dim(3, -2, 1),
        ),
        (
            "bulk-4D Noether force coefficient m*rho4*Q1*Q2",
            {"massDensity4 + 2*q4Vol": dim_add(d["massDensity4"], dim_scale(2, d["q4Vol"]))},
            dim(4, -2, 1),
        ),
        (
            "V_conf body-force volume term dimension",
            {"forceDensity3 + 3L": dim_add(d["forceDensity3"], dim_scale(3, l_unit))},
            dim(1, -2, 1),
        ),
    ]

    return {
        name: dimension_residual(terms, expected_overrides.get(name, expected))
        for name, terms, expected in checks
    }


def derive_eos_identity(s: Symbols) -> sp.Expr:
    h_expr = sp.Rational(5, 4) * s.k_eos * s.rho**4
    p_expr = s.k_eos * s.rho**5
    return sp.diff(p_expr, s.rho) - s.rho * sp.diff(h_expr, s.rho)


def derive_pressure_cross_term(s: Symbols) -> sp.Expr:
    h_expr = sp.Rational(5, 4) * s.k_eos * s.rho**4
    p_expr = s.k_eos * s.rho**5
    dh_drho = sp.diff(h_expr, s.rho)
    dp_drho = sp.diff(p_expr, s.rho)
    delta_rho_cross = -s.m_gnls * s.vdot / dh_drho
    return sp.simplify(dp_drho * delta_rho_cross)


def derive_quantum_residual(s: Symbols) -> sp.Expr:
    u = sp.Function("u")
    ux = u(s.x)
    rho_expr = ux**2
    quantum_potential = -(s.hbar**2 / (2 * s.mass_q)) * sp.diff(ux, s.x, 2) / ux
    sigma_q = (s.hbar**2 / (4 * s.mass_q)) * (
        sp.diff(rho_expr, s.x) ** 2 / rho_expr - sp.diff(rho_expr, s.x, 2)
    )
    return sp.simplify(sp.diff(sigma_q, s.x) - rho_expr * sp.diff(quantum_potential, s.x))


def flux_for_surface(
    *,
    s: Symbols,
    surface: Stage001SurfaceInput,
    density: sp.Expr,
    separation: sp.Symbol,
    config: DerivationConfig,
    convective_moment_override: sp.Expr | None = None,
) -> tuple[sp.Expr, sp.Expr, sp.Expr, sp.Expr]:
    omega = surface.omega_symbol
    radius_power = surface.shell_radius ** (surface.ambient_dim - 1)
    area = omega * radius_power
    pressure_second_moment = area * surface.mu_symbol
    convective_mu = (
        surface.mu_symbol if convective_moment_override is None else convective_moment_override
    )
    convective_second_moment = area * convective_mu

    drain_prefactor = config.drain_sign * s.q2 / (omega * radius_power)
    convective_flux = sp.simplify(
        s.m_gnls * density * s.v1 * drain_prefactor * (area + convective_second_moment)
    )
    pressure_flux = sp.simplify(
        -s.m_gnls * density * s.v1 * drain_prefactor * pressure_second_moment
    )

    substituted = surface.stage001_subs
    convective_flux = sp.simplify(convective_flux.subs(substituted))
    pressure_flux = sp.simplify(pressure_flux.subs(substituted))
    total_flux = sp.simplify(convective_flux + pressure_flux)

    force_along_v1 = sp.simplify(config.force_flux_sign * total_flux)
    v1_from_gauss = -s.q1 / (
        surface.omega_symbol * separation ** (surface.ambient_dim - 1)
    )
    force_after_gauss = sp.simplify(
        force_along_v1.subs(s.v1, v1_from_gauss).subs(substituted)
    )
    return convective_flux, pressure_flux, total_flux, force_after_gauss


def derive_fluxes(
    s: Symbols,
    surface2: Stage001SurfaceInput,
    surface3: Stage001SurfaceInput,
    config: DerivationConfig,
) -> FluxSet:
    conv3, press3, total3, force3 = flux_for_surface(
        s=s,
        surface=surface2,
        density=s.n3,
        separation=s.r12,
        config=config,
        convective_moment_override=config.convective_moment3_override,
    )
    _conv4, _press4, total4, force4 = flux_for_surface(
        s=s,
        surface=surface3,
        density=s.rho4,
        separation=s.radius4,
        config=config,
    )
    return FluxSet(conv3, press3, total3, force3, total4, force4)


def expected_algebra(
    s: Symbols,
    surface2: Stage001SurfaceInput,
    surface3: Stage001SurfaceInput,
    expected_overrides: dict[str, sp.Expr] | None = None,
) -> dict[str, sp.Expr]:
    mu3 = surface2.mu_value
    omega2 = surface2.omega_value
    omega3 = surface3.omega_value
    expected: dict[str, sp.Expr] = {
        "EOS identity dP/drho = rho*dh/drho": sp.Integer(0),
        "Bernoulli/EOS pressure cross term": -s.m_gnls * s.rho * s.vdot,
        "Madelung quantum stress divergence": sp.Integer(0),
        "reduced-3D convective angular traction factor": -(1 + mu3)
        * s.m_gnls
        * s.n3
        * s.q2
        * s.v1,
        "reduced-3D pressure angular traction factor": mu3
        * s.m_gnls
        * s.n3
        * s.q2
        * s.v1,
        "reduced-3D total flux from Noether stress": -s.m_gnls * s.n3 * s.q2 * s.v1,
        "reduced-3D force structure after Gauss substitution": -s.m_gnls
        * s.n3
        * s.q1
        * s.q2
        / (omega2 * s.r12**2),
        "bulk-4D total flux from Noether stress": -s.m_gnls * s.rho4 * s.q2 * s.v1,
        "bulk-4D force structure after Gauss substitution": -s.m_gnls
        * s.rho4
        * s.q1
        * s.q2
        / (omega3 * s.radius4**3),
    }
    if expected_overrides:
        expected.update(expected_overrides)
    return expected


def actual_algebra(
    s: Symbols,
    surface2: Stage001SurfaceInput,
    surface3: Stage001SurfaceInput,
    config: DerivationConfig,
) -> dict[str, sp.Expr]:
    fluxes = derive_fluxes(s, surface2, surface3, config)
    return {
        "EOS identity dP/drho = rho*dh/drho": derive_eos_identity(s),
        "Bernoulli/EOS pressure cross term": derive_pressure_cross_term(s),
        "Madelung quantum stress divergence": derive_quantum_residual(s),
        "reduced-3D convective angular traction factor": fluxes.convective3,
        "reduced-3D pressure angular traction factor": fluxes.pressure3,
        "reduced-3D total flux from Noether stress": fluxes.total3,
        "reduced-3D force structure after Gauss substitution": fluxes.force3,
        "bulk-4D total flux from Noether stress": fluxes.total4,
        "bulk-4D force structure after Gauss substitution": fluxes.force4,
    }


def algebra_residuals(
    s: Symbols,
    surface2: Stage001SurfaceInput,
    surface3: Stage001SurfaceInput,
    config: DerivationConfig,
    expected_overrides: dict[str, sp.Expr] | None = None,
) -> dict[str, sp.Expr]:
    expected = expected_algebra(s, surface2, surface3, expected_overrides)
    actual = actual_algebra(s, surface2, surface3, config)
    return {name: sp.simplify(actual[name] - expected[name]) for name in expected}


def all_residuals(
    *,
    s: Symbols,
    surface2: Stage001SurfaceInput,
    surface3: Stage001SurfaceInput,
    config: DerivationConfig,
    dim_expected_overrides: dict[str, Dim] | None = None,
    algebra_expected_overrides: dict[str, sp.Expr] | None = None,
) -> dict[str, sp.Expr]:
    residuals = dim_claims(dim_expected_overrides)
    residuals.update(
        algebra_residuals(
            s=s,
            surface2=surface2,
            surface3=surface3,
            config=config,
            expected_overrides=algebra_expected_overrides,
        )
    )
    return residuals


@dataclass(frozen=True)
class MutationCase:
    label: str
    config: DerivationConfig
    required_failures: tuple[str, ...]
    dim_expected_overrides: dict[str, Dim] | None = None
    algebra_expected_overrides: dict[str, sp.Expr] | None = None


def mutation_cases(s: Symbols) -> list[MutationCase]:
    reduced_force = "reduced-3D force structure after Gauss substitution"
    bulk_force = "bulk-4D force structure after Gauss substitution"
    dim_force = "reduced-3D surface traction integral gives force"
    convective = "reduced-3D convective angular traction factor"
    total = "reduced-3D total flux from Noether stress"

    return [
        MutationCase(
            label="flip sign of reduced-3D force law",
            config=DerivationConfig(),
            required_failures=(reduced_force,),
            algebra_expected_overrides={
                reduced_force: s.m_gnls * s.n3 * s.q1 * s.q2 / (4 * sp.pi * s.r12**2)
            },
        ),
        MutationCase(
            label="change 4*pi -> 2*pi in reduced-3D force law",
            config=DerivationConfig(),
            required_failures=(reduced_force,),
            algebra_expected_overrides={
                reduced_force: -s.m_gnls * s.n3 * s.q1 * s.q2 / (2 * sp.pi * s.r12**2)
            },
        ),
        MutationCase(
            label="change 2*pi^2 -> 4*pi^2 in bulk-4D force law",
            config=DerivationConfig(),
            required_failures=(bulk_force,),
            algebra_expected_overrides={
                bulk_force: -s.m_gnls * s.rho4 * s.q1 * s.q2 / (4 * sp.pi**2 * s.radius4**3)
            },
        ),
        MutationCase(
            label="corrupt reduced-3D surface-traction dimensional exponent",
            config=DerivationConfig(),
            required_failures=(dim_force,),
            dim_expected_overrides={dim_force: dim(2, -2, 1)},
        ),
        MutationCase(
            label="corrupt consumed convective second moment 1/3 -> 1/2",
            config=DerivationConfig(convective_moment3_override=sp.Rational(1, 2)),
            required_failures=(convective, total, reduced_force),
        ),
        MutationCase(
            label="flip F = -int Pi.n force-definition sign",
            config=DerivationConfig(force_flux_sign=sp.Integer(1)),
            required_failures=(reduced_force,),
        ),
        MutationCase(
            label="flip v2 drain-inflow Gauss sign",
            config=DerivationConfig(drain_sign=sp.Integer(1)),
            required_failures=(reduced_force,),
        ),
    ]


def run_baseline(
    s: Symbols,
    surface2: Stage001SurfaceInput,
    surface3: Stage001SurfaceInput,
) -> None:
    residuals = all_residuals(
        s=s,
        surface2=surface2,
        surface3=surface3,
        config=DerivationConfig(),
    )

    subbanner("Baseline dimensional residuals")
    for name in list(dim_claims()):
        expect_zero(name, residuals[name])

    subbanner("Baseline algebraic residuals")
    for name in expected_algebra(s, surface2, surface3):
        expect_zero(name, residuals[name])

    print("")
    print("Verdict labels:")
    print("  leading matter-stress sign: FORCE_ATTRACTIVE_DERIVED")
    print("  full sign: SIGN_RESIDUAL_QUANTUM_VCONF_MAXWELL_PROFILE")
    print("  acceptance: PASS_WITH_NAMED_RESIDUALS")


def assert_mutation_fails(
    s: Symbols,
    surface2: Stage001SurfaceInput,
    surface3: Stage001SurfaceInput,
    case: MutationCase,
) -> None:
    residuals = all_residuals(
        s=s,
        surface2=surface2,
        surface3=surface3,
        config=case.config,
        dim_expected_overrides=case.dim_expected_overrides,
        algebra_expected_overrides=case.algebra_expected_overrides,
    )
    failed_required: list[str] = []
    for check_name in case.required_failures:
        try:
            expect_zero(check_name, residuals[check_name])
        except AuditFailure:
            failed_required.append(check_name)

    if tuple(failed_required) == case.required_failures:
        failed_list = ", ".join(case.required_failures)
        print(f"PASS  mutation probe: {case.label} produced required FAIL ({failed_list})")
        return

    missing = sorted(set(case.required_failures) - set(failed_required))
    print(f"FAIL  mutation probe: {case.label} survived for {', '.join(missing)}")
    raise AuditFailure(f"mutation probe survived: {case.label}")


def run_mutation_probe(
    s: Symbols,
    surface2: Stage001SurfaceInput,
    surface3: Stage001SurfaceInput,
) -> None:
    subbanner("Able-to-fail mutation probe")
    for case in mutation_cases(s):
        assert_mutation_fails(s, surface2, surface3, case)


def main() -> None:
    banner("ledger_stage002_matter_stress_force_assembly SymPy audit")
    s = make_symbols()
    surface2, surface3 = stage001_surfaces()
    run_baseline(s, surface2, surface3)
    run_mutation_probe(s, surface2, surface3)
    print("")
    print("OVERALL PASS: SymPy derived the stage002 matter-stress force assembly exactly")


if __name__ == "__main__":
    main()
