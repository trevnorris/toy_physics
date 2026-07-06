#!/usr/bin/env python3
"""SymPy companion for pathA_21c Noether-stress force checks.

Print-only engine:

  timeout 600 python3 software/stage1_solver/tools/pathA_21c_force_from_noether_stress_tensor_sympy.py

The script intentionally does not read the Mathematica harness, report, scratch
JSON, or any generated engine output.  Expected claims are local constants; the
derived side is recomputed from dimension bases, EOS/quantum formulas, and
explicit S2/S3 surface integrals.
"""

from __future__ import annotations

from dataclasses import dataclass
import random
import sys
from typing import Any

import sympy as sp


Dim = tuple[sp.Rational, sp.Rational, sp.Rational]


def q(value: int | str) -> sp.Rational:
    return sp.Rational(value)


def dim(l_power: int, t_power: int, m_power: int) -> Dim:
    return (q(l_power), q(t_power), q(m_power))


DIM0 = dim(0, 0, 0)


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
    scale = q(scale) if isinstance(scale, int) else scale
    return (scale * item[0], scale * item[1], scale * item[2])


def dim_text(item: Dim) -> str:
    pieces: list[str] = []
    for label, power in (("L", item[0]), ("T", item[1]), ("M", item[2])):
        if power == 0:
            continue
        if power == 1:
            pieces.append(label)
        else:
            pieces.append(f"{label}^{power}")
    return " ".join(pieces) if pieces else "1"


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
    m: sp.Symbol
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
        m=sp.Symbol("m", nonzero=True),
        x=sp.Symbol("x", real=True),
    )


@dataclass(frozen=True)
class SurfaceData:
    label: str
    ambient_dim: int
    area_integral: sp.Expr
    omega: sp.Expr
    moment11: sp.Expr
    cross12: sp.Expr

    @property
    def solid_angle_factor(self) -> sp.Expr:
        return sp.simplify(self.area_integral / self.omega)


def derive_s2_surface() -> SurfaceData:
    theta, phi = sp.symbols("theta phi", real=True)
    n1 = sp.cos(theta)
    n2 = sp.sin(theta) * sp.cos(phi)
    jacobian = sp.sin(theta)
    ranges = ((phi, 0, 2 * sp.pi), (theta, 0, sp.pi))

    omega = sp.integrate(jacobian, *ranges)
    moment11 = sp.simplify(sp.integrate(n1**2 * jacobian, *ranges) / omega)
    cross12 = sp.simplify(sp.integrate(n1 * n2 * jacobian, *ranges) / omega)
    return SurfaceData("S2", 3, sp.simplify(omega), sp.simplify(omega), moment11, cross12)


def derive_s3_surface() -> SurfaceData:
    chi, theta, phi = sp.symbols("chi theta phi", real=True)
    n1 = sp.cos(chi)
    n2 = sp.sin(chi) * sp.cos(theta)
    jacobian = sp.sin(chi) ** 2 * sp.sin(theta)
    ranges = ((phi, 0, 2 * sp.pi), (theta, 0, sp.pi), (chi, 0, sp.pi))

    omega = sp.integrate(jacobian, *ranges)
    moment11 = sp.simplify(sp.integrate(n1**2 * jacobian, *ranges) / omega)
    cross12 = sp.simplify(sp.integrate(n1 * n2 * jacobian, *ranges) / omega)
    return SurfaceData("S3", 4, sp.simplify(omega), sp.simplify(omega), moment11, cross12)


@dataclass(frozen=True)
class DerivationConfig:
    drain_sign: sp.Integer = sp.Integer(-1)
    force_flux_sign: sp.Integer = sp.Integer(-1)
    convective_moment3_override: sp.Expr | None = None


@dataclass(frozen=True)
class CheckResult:
    name: str
    kind: str
    passed: bool
    expected: Any
    actual: Any
    residual: sp.Expr | None = None
    note: str = ""


def compact(expr: Any) -> str:
    if isinstance(expr, sp.Basic):
        return sp.sstr(sp.factor(sp.cancel(sp.simplify(expr))))
    return str(expr)


def assert_no_float(expr: sp.Expr, name: str) -> None:
    if expr.atoms(sp.Float):
        raise AssertionError(f"{name}: Float atom found in checked expression {expr}")


def rational_samples(symbols: set[sp.Symbol], count: int = 3) -> list[dict[sp.Symbol, sp.Rational]]:
    ordered = sorted(symbols, key=lambda sym: sym.name)
    rng = random.Random(21021)
    samples: list[dict[sp.Symbol, sp.Rational]] = []
    for _ in range(count):
        point: dict[sp.Symbol, sp.Rational] = {}
        for sym in ordered:
            numerator = rng.randint(1, 17)
            denominator = rng.randint(1, 11)
            point[sym] = sp.Rational(numerator, denominator)
        samples.append(point)
    return samples


def expressions_equal(derived: sp.Expr, expected: sp.Expr) -> tuple[bool, sp.Expr, str]:
    assert_no_float(derived, "derived")
    assert_no_float(expected, "expected")
    residual = sp.simplify(derived - expected)
    if residual == 0:
        return True, residual, "symbolic"

    nonzero_points: list[str] = []
    for point in rational_samples(residual.free_symbols, 3):
        value = sp.simplify(residual.subs(point))
        if value == 0:
            continue
        numeric = sp.N(value, 60)
        if numeric != 0:
            nonzero_points.append(compact(value))
    if not nonzero_points:
        return True, residual, "symbolic residual survived simplify but vanished at 3 rational probes"
    return False, residual, f"nonzero rational probe residual(s): {', '.join(nonzero_points)}"


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


def dim_check(
    name: str,
    terms: dict[str, Dim],
    expected: Dim,
    expected_overrides: dict[str, Dim],
) -> CheckResult:
    expected = expected_overrides.get(name, expected)
    passed = all(actual == expected for actual in terms.values())
    return CheckResult(name, "dim", passed, expected, terms)


def derive_dim_checks(expected_overrides: dict[str, Dim] | None = None) -> list[CheckResult]:
    expected_overrides = expected_overrides or {}
    d = derive_dimensions()
    l_unit = d["L"]
    t_unit = d["T"]

    checks = [
        dim_check(
            "reduced-3D momentum-balance terms",
            {
                "partial_t(m*N*v_i)": dim_sub(dim_add(d["massDensity3"], d["velocity"]), t_unit),
                "partial_j Pi_ij": dim_sub(d["stress3"], l_unit),
                "V_conf body force": dim_add(d["rho3"], d["energyGradient"]),
            },
            dim(-2, -2, 1),
            expected_overrides,
        ),
        dim_check(
            "bulk-4D momentum-balance terms",
            {
                "partial_t(m*rho*v_i)": dim_sub(dim_add(d["massDensity4"], d["velocity"]), t_unit),
                "partial_J Pi_iJ": dim_sub(d["stress4"], l_unit),
                "V_conf body force": dim_add(d["rho4"], d["energyGradient"]),
            },
            dim(-3, -2, 1),
            expected_overrides,
        ),
        dim_check(
            "reduced-3D Noether stress representatives",
            {
                "convective": dim_add(d["massDensity3"], dim_scale(2, d["velocity"])),
                "pressure": dim_add(dim_sub(d["energy"], dim_scale(4, l_unit)), l_unit),
                "quantum": dim_add(d["quantumPrefactor"], d["rho3"]),
            },
            dim(-1, -2, 1),
            expected_overrides,
        ),
        dim_check(
            "Euler force-per-volume identity terms in reduced 3D",
            {
                "m*N*acceleration": dim_sub(dim_add(d["massDensity3"], d["velocity"]), t_unit),
                "N*grad h": dim_add(d["rho3"], d["energyGradient"]),
                "N*grad Q": dim_add(d["rho3"], d["energyGradient"]),
                "N*grad V_conf": dim_add(d["rho3"], d["energyGradient"]),
                "N*q(E+vB)": dim_add(d["rho3"], d["energyGradient"]),
            },
            dim(-2, -2, 1),
            expected_overrides,
        ),
        dim_check(
            "reduced-3D surface traction integral gives force",
            {"stress3 + 2L": dim_add(d["stress3"], dim_scale(2, l_unit))},
            dim(1, -2, 1),
            expected_overrides,
        ),
        dim_check(
            "bulk-4D surface traction integral gives force",
            {"stress4 + 3L": dim_add(d["stress4"], dim_scale(3, l_unit))},
            dim(1, -2, 1),
            expected_overrides,
        ),
        dim_check(
            "reduced-3D Noether force coefficient m*N*Q1*Q2",
            {"massDensity3 + 2*q3Vol": dim_add(d["massDensity3"], dim_scale(2, d["q3Vol"]))},
            dim(3, -2, 1),
            expected_overrides,
        ),
        dim_check(
            "bulk-4D Noether force coefficient m*rho4*Q1*Q2",
            {"massDensity4 + 2*q4Vol": dim_add(d["massDensity4"], dim_scale(2, d["q4Vol"]))},
            dim(4, -2, 1),
            expected_overrides,
        ),
        dim_check(
            "V_conf body-force volume term dimension",
            {"forceDensity3 + 3L": dim_add(d["forceDensity3"], dim_scale(3, l_unit))},
            dim(1, -2, 1),
            expected_overrides,
        ),
    ]
    return checks


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
    quantum_potential = -(s.hbar**2 / (2 * s.m)) * sp.diff(ux, s.x, 2) / ux
    sigma_q = (s.hbar**2 / (4 * s.m)) * (
        sp.diff(rho_expr, s.x) ** 2 / rho_expr - sp.diff(rho_expr, s.x, 2)
    )
    return sp.simplify(sp.diff(sigma_q, s.x) - rho_expr * sp.diff(quantum_potential, s.x))


@dataclass(frozen=True)
class FluxSet:
    convective3: sp.Expr
    pressure3: sp.Expr
    total3: sp.Expr
    force3: sp.Expr
    total4: sp.Expr
    force4: sp.Expr


def flux_for_surface(
    *,
    s: Symbols,
    surface: SurfaceData,
    density: sp.Expr,
    separation: sp.Symbol,
    config: DerivationConfig,
    convective_moment_override: sp.Expr | None = None,
) -> tuple[sp.Expr, sp.Expr, sp.Expr, sp.Expr]:
    convective_moment = (
        surface.moment11 if convective_moment_override is None else convective_moment_override
    )
    pressure_moment = surface.moment11
    solid_angle_part = surface.solid_angle_factor

    # v2 = drain_sign * Q2*n/(Omega*a^(d-1)); the surface factor cancels.
    convective_flux = sp.simplify(
        s.m_gnls
        * density
        * s.q2
        * s.v1
        * config.drain_sign
        * (solid_angle_part + convective_moment)
    )
    pressure_flux = sp.simplify(
        -config.drain_sign * s.m_gnls * density * s.q2 * s.v1 * pressure_moment
    )
    total_flux = sp.simplify(convective_flux + pressure_flux)

    force_along_v1 = sp.simplify(config.force_flux_sign * total_flux)
    v1_from_gauss = -s.q1 / (surface.omega * separation ** (surface.ambient_dim - 1))
    force_after_gauss = sp.simplify(force_along_v1.subs(s.v1, v1_from_gauss))
    return convective_flux, pressure_flux, total_flux, force_after_gauss


def derive_fluxes(
    s: Symbols,
    surface2: SurfaceData,
    surface3: SurfaceData,
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


def expected_algebra(s: Symbols, expected_overrides: dict[str, sp.Expr] | None = None) -> dict[str, sp.Expr]:
    expected: dict[str, sp.Expr] = {
        "EOS identity dP/drho = rho*dh/drho": sp.Integer(0),
        "Bernoulli/EOS pressure cross term": -s.m_gnls * s.rho * s.vdot,
        "Madelung quantum stress divergence": sp.Integer(0),
        "reduced-3D convective angular traction factor": -sp.Rational(4, 3)
        * s.m_gnls
        * s.n3
        * s.q2
        * s.v1,
        "reduced-3D pressure angular traction factor": sp.Rational(1, 3)
        * s.m_gnls
        * s.n3
        * s.q2
        * s.v1,
        "reduced-3D total flux from Noether stress": -s.m_gnls * s.n3 * s.q2 * s.v1,
        "reduced-3D force structure after Gauss substitution": -s.m_gnls
        * s.n3
        * s.q1
        * s.q2
        / (4 * sp.pi * s.r12**2),
        "bulk-4D total flux from Noether stress": -s.m_gnls * s.rho4 * s.q2 * s.v1,
        "bulk-4D force structure after Gauss substitution": -s.m_gnls
        * s.rho4
        * s.q1
        * s.q2
        / (2 * sp.pi**2 * s.radius4**3),
    }
    if expected_overrides:
        expected.update(expected_overrides)
    return expected


def actual_algebra(
    s: Symbols,
    surface2: SurfaceData,
    surface3: SurfaceData,
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


def derive_algebra_checks(
    s: Symbols,
    surface2: SurfaceData,
    surface3: SurfaceData,
    config: DerivationConfig,
    expected_overrides: dict[str, sp.Expr] | None = None,
) -> list[CheckResult]:
    expected = expected_algebra(s, expected_overrides)
    actual = actual_algebra(s, surface2, surface3, config)
    results: list[CheckResult] = []
    for name in expected:
        passed, residual, note = expressions_equal(actual[name], expected[name])
        results.append(CheckResult(name, "alg", passed, expected[name], actual[name], residual, note))
    return results


def surface_primitive_checks(surface2: SurfaceData, surface3: SurfaceData) -> list[CheckResult]:
    primitives = [
        ("Omega_2 from explicit S2 integral", surface2.omega, 4 * sp.pi),
        ("<n_x1^2>_S2 from explicit S2 integral", surface2.moment11, sp.Rational(1, 3)),
        ("<n_x1 n_x2>_S2 from explicit S2 integral", surface2.cross12, sp.Integer(0)),
        ("Omega_3 from explicit S3 integral", surface3.omega, 2 * sp.pi**2),
        ("<n_x1^2>_S3 from explicit S3 integral", surface3.moment11, sp.Rational(1, 4)),
        ("<n_x1 n_x2>_S3 from explicit S3 integral", surface3.cross12, sp.Integer(0)),
    ]
    results: list[CheckResult] = []
    for name, actual, expected in primitives:
        passed, residual, note = expressions_equal(actual, expected)
        results.append(CheckResult(name, "surface", passed, expected, actual, residual, note))
    return results


def run_claim_checks(
    *,
    s: Symbols,
    surface2: SurfaceData,
    surface3: SurfaceData,
    config: DerivationConfig,
    dim_expected_overrides: dict[str, Dim] | None = None,
    algebra_expected_overrides: dict[str, sp.Expr] | None = None,
) -> list[CheckResult]:
    return derive_dim_checks(dim_expected_overrides) + derive_algebra_checks(
        s, surface2, surface3, config, algebra_expected_overrides
    )


def print_result(result: CheckResult) -> None:
    if result.passed:
        suffix = "" if result.note in ("", "symbolic") else f" ({result.note})"
        print(f"{result.name} -> PASS{suffix}")
        return

    if result.kind == "dim":
        expected = dim_text(result.expected)
        actual_parts = ", ".join(
            f"{term}={dim_text(value)}" for term, value in result.actual.items()
        )
        print(f"{result.name} -> FAIL: expected {expected}; actual {actual_parts}")
        return

    print(
        f"{result.name} -> FAIL: expected {compact(result.expected)}; "
        f"derived {compact(result.actual)}; residual {compact(result.residual)}"
    )


@dataclass(frozen=True)
class MutationCase:
    name: str
    config: DerivationConfig
    required_failed_checks: tuple[str, ...]
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
            name="expected flip sign of reduced-3D force law",
            config=DerivationConfig(),
            required_failed_checks=(reduced_force,),
            algebra_expected_overrides={
                reduced_force: s.m_gnls * s.n3 * s.q1 * s.q2 / (4 * sp.pi * s.r12**2)
            },
        ),
        MutationCase(
            name="expected change 4*pi -> 2*pi in reduced-3D force law",
            config=DerivationConfig(),
            required_failed_checks=(reduced_force,),
            algebra_expected_overrides={
                reduced_force: -s.m_gnls * s.n3 * s.q1 * s.q2 / (2 * sp.pi * s.r12**2)
            },
        ),
        MutationCase(
            name="expected change 2*pi^2 -> 4*pi^2 in bulk-4D force law",
            config=DerivationConfig(),
            required_failed_checks=(bulk_force,),
            algebra_expected_overrides={
                bulk_force: -s.m_gnls * s.rho4 * s.q1 * s.q2 / (4 * sp.pi**2 * s.radius4**3)
            },
        ),
        MutationCase(
            name="expected corrupt one dimensional exponent",
            config=DerivationConfig(),
            required_failed_checks=(dim_force,),
            dim_expected_overrides={dim_force: dim(2, -2, 1)},
        ),
        MutationCase(
            name="derivation corrupt reduced-3D convective second moment 1/3 -> 1/2",
            config=DerivationConfig(convective_moment3_override=sp.Rational(1, 2)),
            required_failed_checks=(convective, total, reduced_force),
        ),
        MutationCase(
            name="derivation flip F = -Phi force-definition sign",
            config=DerivationConfig(force_flux_sign=sp.Integer(1)),
            required_failed_checks=(reduced_force,),
        ),
        MutationCase(
            name="derivation flip v2 drain inflow sign",
            config=DerivationConfig(drain_sign=sp.Integer(1)),
            required_failed_checks=(reduced_force,),
        ),
    ]


def run_mutation_probe(
    s: Symbols,
    surface2: SurfaceData,
    surface3: SurfaceData,
) -> bool:
    print("Mutation probe:")
    all_good = True
    for case in mutation_cases(s):
        results = run_claim_checks(
            s=s,
            surface2=surface2,
            surface3=surface3,
            config=case.config,
            dim_expected_overrides=case.dim_expected_overrides,
            algebra_expected_overrides=case.algebra_expected_overrides,
        )
        failed = {result.name for result in results if not result.passed}
        required = set(case.required_failed_checks)
        probe_passed = required.issubset(failed)
        all_good = all_good and probe_passed
        if probe_passed:
            failed_list = ", ".join(name for name in case.required_failed_checks)
            print(f"  {case.name} -> FAILED as required ({failed_list})")
        else:
            missing = ", ".join(sorted(required - failed)) or "none"
            observed = ", ".join(sorted(failed)) or "none"
            print(
                f"  {case.name} -> PROBE FAIL: required failure missing for {missing}; "
                f"observed failed checks: {observed}"
            )
    print(
        "Mutation probe summary: "
        + ("PASS (all corruptions failed)" if all_good else "FAIL (a corruption survived)")
    )
    return all_good


def main() -> int:
    s = make_symbols()
    surface2 = derive_s2_surface()
    surface3 = derive_s3_surface()

    surface_results = surface_primitive_checks(surface2, surface3)
    claim_results = run_claim_checks(
        s=s,
        surface2=surface2,
        surface3=surface3,
        config=DerivationConfig(),
    )

    print("Surface integral primitives:")
    for result in surface_results:
        print_result(result)

    print("Baseline claim checks:")
    for result in claim_results:
        print_result(result)

    surface_passed = all(result.passed for result in surface_results)
    baseline_passed = all(result.passed for result in claim_results)
    print(
        f"Baseline summary: {sum(result.passed for result in claim_results)}/"
        f"{len(claim_results)} claim checks passed; "
        f"{sum(result.passed for result in surface_results)}/{len(surface_results)} "
        "surface primitives passed"
    )

    mutation_passed = run_mutation_probe(s, surface2, surface3)

    overall = surface_passed and baseline_passed and mutation_passed
    print(f"pathA_21c SymPy companion: {'PASS' if overall else 'FAIL'}")
    return 0 if overall else 1


if __name__ == "__main__":
    sys.exit(main())
