"""Deterministic headless numeric verification report.

This module is the text-only verification surface required by
``notes/build_spec.md`` section 6.1.  All force and trajectory values come
from the existing render-agnostic :mod:`physics` package; this module only
samples those APIs, measures their outputs, and formats the results.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

if __package__:
    from .params import DEFAULT_PARAMS, PARAMETER_INFO, ModelParameters, value_of
    from .physics import charge, gravity, light, magnetism
else:  # Support ``python software/force_visualizer/report.py`` from the repo root.
    import sys

    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
    from software.force_visualizer.params import (  # type: ignore[no-redef]
        DEFAULT_PARAMS,
        PARAMETER_INFO,
        ModelParameters,
        value_of,
    )
    from software.force_visualizer.physics import (  # type: ignore[no-redef]
        charge,
        gravity,
        light,
        magnetism,
    )


OUTPUT_PATH = Path(__file__).resolve().parent / "output" / "verification_report.txt"

MASS1 = 1.0
MASS2 = 0.20
SEMI_MAJOR_AXIS = 3.0
ECCENTRICITY = 0.25
ORBIT_PERIODS = 5
STEPS_PER_PERIOD = 1800

PARAMETER_ORDER = (
    "c_s",
    "lambda_gamma",
    "c_gamma",
    "rho_br",
    "mu_R",
    "c_E",
    "rho_B0",
    "chi_c",
    "C_hu",
    "B_eff",
    "G",
    "Q_E",
    "aT",
    "aL",
    "N_u",
    "ell",
    "mouth_half_width_b",
    "M_h",
    "N0",
    "yukawa_mass",
    "yukawa_fraction",
    "epsilon0",
    "epsilon1",
    "radiation_length_a",
    "J_phase",
    "kappa_phase",
    "longitudinal_display_fraction",
    "radiation_reaction_benchmark_scale",
)

PARAMETER_UNITS = {
    "c_s": "length_unit/time_unit",
    "lambda_gamma": "1",
    "c_gamma": "length_unit/time_unit",
    "rho_br": "density_unit",
    "mu_R": "modulus_unit",
    "c_E": "length_unit/time_unit",
    "rho_B0": "density_unit",
    "chi_c": "susceptibility_unit",
    "C_hu": "coupling_unit",
    "B_eff": "modulus_unit",
    "G": "length_unit^3/(mass_unit*time_unit^2)",
    "Q_E": "charge_unit",
    "aT": "amplitude_unit",
    "aL": "amplitude_unit",
    "N_u": "normalization_unit",
    "ell": "length_unit",
    "mouth_half_width_b": "length_unit",
    "M_h": "inertia_unit",
    "N0": "1/length_unit",
    "yukawa_mass": "1/length_unit",
    "yukawa_fraction": "1",
    "epsilon0": "1",
    "epsilon1": "1",
    "radiation_length_a": "length_unit",
    "J_phase": "phase_coupling_unit",
    "kappa_phase": "phase_stiffness_unit",
    "longitudinal_display_fraction": "1",
    "radiation_reaction_benchmark_scale": "1",
}


@dataclass
class CheckTally:
    """Count report checks without hiding an individual failure."""

    passed: int = 0
    total: int = 0

    def mark(self, condition: bool) -> str:
        self.total += 1
        if condition:
            self.passed += 1
            return "✓"
        return "✗"


@dataclass(frozen=True)
class ReportResult:
    """Complete transcript and its machine-usable check totals."""

    transcript: str
    passed: int
    total: int


def _vector(values: np.ndarray) -> str:
    return "(" + ", ".join(f"{float(value):+.8e}" for value in values) + ")"


def _signed_radial_force(force: np.ndarray, outward: np.ndarray) -> float:
    return float(np.dot(force, outward / np.linalg.norm(outward)))


def _charge_force(params: ModelParameters, distance: float, orientation2: float) -> np.ndarray:
    return charge.force_on_first(
        [distance, 0.0],
        [0.0, 0.0],
        1.0,
        orientation2,
        params.Q_E,
        params.ell,
        params.mouth_half_width_b,
    )


def _magnetic_force(
    params: ModelParameters,
    distance: float,
    velocity2: tuple[float, float] = (0.0, 1.0),
    *,
    include_scalar: bool = True,
) -> np.ndarray:
    return magnetism.force_on_second(
        [-distance / 2.0, 0.0],
        [distance / 2.0, 0.0],
        [0.0, 1.0],
        velocity2,
        1.0,
        1.0,
        params.N_u,
        params.aT,
        params.aL,
        params.mu_R,
        params.B_eff,
        include_scalar,
    )


def _add_parameter_block(lines: list[str], params: ModelParameters) -> None:
    lines.extend(("PARAMETER BLOCK", "All quantities use deterministic visualization units."))
    for name in PARAMETER_ORDER:
        info = PARAMETER_INFO[name]
        lines.append(
            f"  [{info.status.value}] {name}={value_of(params, name):.10e} "
            f"{PARAMETER_UNITS[name]}"
        )


def _add_gravity(lines: list[str], tally: CheckTally, params: ModelParameters) -> None:
    total_mass = MASS1 + MASS2
    initial_position, initial_velocity = gravity.kepler_periapsis_state(
        SEMI_MAJOR_AXIS, ECCENTRICITY, total_mass, params.G
    )
    period = gravity.kepler_period(SEMI_MAJOR_AXIS, total_mass, params.G)
    dt = period / STEPS_PER_PERIOD
    orbit = gravity.integrate_relative_orbit(
        initial_position,
        initial_velocity,
        MASS1,
        MASS2,
        params.G,
        params.c_gamma,
        dt,
        ORBIT_PERIODS * STEPS_PER_PERIOD,
        pn_order=0.0,
    )
    radius = np.linalg.norm(orbit.positions, axis=1)
    speed_squared = np.sum(orbit.velocities**2, axis=1)
    reduced_mass = MASS1 * MASS2 / total_mass
    total_energy = reduced_mass * (
        0.5 * speed_squared - params.G * total_mass / radius
    )
    energy_drift = abs(
        float((total_energy[-1] - total_energy[0]) / total_energy[0])
    )
    energy_tolerance = 1.0e-9

    pn_orbit = gravity.integrate_relative_orbit(
        initial_position,
        initial_velocity,
        MASS1,
        MASS2,
        params.G,
        params.c_gamma,
        dt,
        int(8.2 * STEPS_PER_PERIOD),
        pn_order=1.0,
    )
    measured_precession = gravity.measured_precession(pn_orbit.positions)
    expected_precession = gravity.analytic_perihelion_precession(
        SEMI_MAJOR_AXIS,
        ECCENTRICITY,
        total_mass,
        params.G,
        params.c_gamma,
    )
    precession_ok = measured_precession > 0.0 and bool(
        np.isclose(measured_precession, expected_precession, rtol=0.035, atol=2.0e-4)
    )

    lines.extend(
        (
            "",
            "GRAVITY",
            f"  orbit setup: masses=({MASS1:.2f}, {MASS2:.2f}) mass_unit; "
            f"a={SEMI_MAJOR_AXIS:.2f} length_unit; e={ECCENTRICITY:.2f}; "
            f"N={ORBIT_PERIODS} periods; samples/period={STEPS_PER_PERIOD}",
            f"  sampled 0PN state at t={orbit.times[-1]:.8e} time_unit: "
            f"r={_vector(orbit.positions[-1])} length_unit; "
            f"v={_vector(orbit.velocities[-1])} length_unit/time_unit",
            f"  [{tally.mark(energy_drift <= energy_tolerance)}] 0PN relative total-energy drift: "
            f"measured={energy_drift:.8e} 1; expected=0.00000000e+00 1; "
            f"tolerance<={energy_tolerance:.1e} 1",
            f"  [{tally.mark(precession_ok)}] 1PN perihelion precession/orbit: "
            f"measured={measured_precession:.8e} rad/orbit; "
            f"expected={expected_precession:.8e} rad/orbit; tolerance=3.5%",
            "  optional 2.5PN inspiral-rate check: not included; the native normalization "
            "is GENUINE_BLOCKED and the configured scale is a CALIBRATED benchmark.",
        )
    )
    sample_point = np.array([2.0, 1.0])
    inflow = gravity.drain_inflow_field(
        sample_point, [[0.0, 0.0]], [MASS1], params.G
    )
    expected_direction = -sample_point / np.linalg.norm(sample_point)
    alignment = float(np.dot(inflow / np.linalg.norm(inflow), expected_direction))
    lines.append(
        f"  [{tally.mark(alignment > 1.0 - 1.0e-12)}] gravity drain-field direction: "
        f"sample={_vector(sample_point)} length_unit; measured={_vector(inflow)} field_unit; "
        f"toward-mass alignment={alignment:.8e}; expected=1.00000000e+00"
    )


def _add_light(lines: list[str], tally: CheckTally, params: ModelParameters) -> None:
    wavenumbers = np.array([1.0, 2.0, 3.0, 4.0])
    measured = np.array(
        [
            light.measured_transverse_angular_frequency(
                light.evolve_transverse_wave(
                    wavenumber,
                    params.c_gamma,
                    spatial_points=384,
                    periods=4.0,
                )
            )
            for wavenumber in wavenumbers
        ]
    )
    expected = params.c_gamma * wavenumbers
    lines.extend(("", "LIGHT", "  dispersion omega(k), measured from finite-difference time evolution:"))
    for k_value, measured_value, expected_value in zip(wavenumbers, measured, expected):
        ok = bool(np.isclose(measured_value, expected_value, rtol=4.0e-4, atol=1.0e-6))
        lines.append(
            f"    [{tally.mark(ok)}] k={k_value:.8e} 1/length_unit: "
            f"measured={measured_value:.8e} 1/time_unit; "
            f"expected(c_gamma*k)={expected_value:.8e} 1/time_unit; tolerance=0.04%"
        )

    wavevector = np.array([1.2, -0.7, 2.3])
    polarizations = light.transverse_polarizations(wavevector)
    transverse_count = polarizations.shape[0]
    transverse_residual = float(np.max(np.abs(polarizations @ wavevector)))
    polarization_ok = transverse_count == 2 and transverse_residual <= 1.0e-14
    lines.append(
        f"  [{tally.mark(polarization_ok)}] transverse polarization count: "
        f"measured={transverse_count:d} modes; expected=2 modes; "
        f"max(k dot epsilon)={transverse_residual:.8e} 1/length_unit"
    )

    impact = 2.5
    ray = light.trace_ray(impact, 20.0, 1.0, params.G, params.c_gamma, steps=3000)
    expected_deflection = light.analytic_signed_deflection(
        impact, 1.0, params.G, params.c_gamma
    )
    lens_ok = ray.signed_deflection < 0.0 and bool(
        np.isclose(ray.signed_deflection, expected_deflection, rtol=0.02, atol=1.0e-5)
    )
    lines.append(
        f"  [{tally.mark(lens_ok)}] lensing deflection (positive impact bends toward mass): "
        f"measured={ray.signed_deflection:.8e} rad; expected={expected_deflection:.8e} rad; "
        "expected_sign=negative; tolerance=2%"
    )


def _add_charge(lines: list[str], tally: CheckTally, params: ModelParameters) -> None:
    distances = np.array([0.75, 1.0, 2.0, 4.0])
    magnitudes = np.array(
        [np.linalg.norm(_charge_force(params, distance, 1.0)) for distance in distances]
    )
    measured_power = -float(np.polyfit(np.log(distances), np.log(magnitudes), 1)[0])
    lines.extend(("", "CHARGE", "  like-pair force samples:"))
    for distance, magnitude in zip(distances, magnitudes):
        lines.append(
            f"    R={distance:.8e} length_unit; measured |F|={magnitude:.8e} force_unit"
        )
    lines.append(
        f"  [{tally.mark(bool(np.isclose(measured_power, 2.0, rtol=0.0, atol=1.0e-12)))}] "
        f"fitted force exponent p: measured={measured_power:.8e} 1; "
        "expected=2.00000000e+00 1"
    )

    like_force = _signed_radial_force(_charge_force(params, 1.0, 1.0), np.array([1.0, 0.0]))
    unlike_force = _signed_radial_force(
        _charge_force(params, 1.0, -1.0), np.array([1.0, 0.0])
    )
    lines.extend(
        (
            f"  [{tally.mark(like_force > 0.0)}] like-pair radial sign: "
            f"measured_force={like_force:+.8e} force_unit, measured_sign={np.sign(like_force):+.0f}; "
            "expected_sign=+1 (repel)",
            f"  [{tally.mark(unlike_force < 0.0)}] unlike-pair radial sign: "
            f"measured_force={unlike_force:+.8e} force_unit, measured_sign={np.sign(unlike_force):+.0f}; "
            "expected_sign=-1 (attract)",
        )
    )
    sample_point = np.array([1.0, 0.0])
    positive_field = charge.electric_field(
        sample_point, [[0.0, 0.0]], [1.0], params.Q_E, params.ell, params.mouth_half_width_b
    )
    negative_field = charge.electric_field(
        sample_point, [[0.0, 0.0]], [-1.0], params.Q_E, params.ell, params.mouth_half_width_b
    )
    positive_radial = float(np.dot(positive_field, sample_point))
    negative_radial = float(np.dot(negative_field, sample_point))
    charge_field_ok = positive_radial > 0.0 and negative_radial < 0.0
    lines.append(
        f"  [{tally.mark(charge_field_ok)}] electric-field directions at +x: "
        f"positive-source E={_vector(positive_field)} field_unit (expected away, +x); "
        f"negative-source E={_vector(negative_field)} field_unit (expected toward, -x)"
    )


def _add_magnetism(lines: list[str], tally: CheckTally, params: ModelParameters) -> None:
    distances = np.array([0.8, 1.25, 2.0, 3.5, 6.0])
    magnitudes = np.array(
        [np.linalg.norm(_magnetic_force(params, distance)) for distance in distances]
    )
    measured_power = -float(np.polyfit(np.log(distances), np.log(magnitudes), 1)[0])
    lines.extend(("", "MAGNETISM", "  parallel-current force samples (scalar included):"))
    for distance, magnitude in zip(distances, magnitudes):
        lines.append(
            f"    R={distance:.8e} length_unit; measured |F|={magnitude:.8e} force_unit"
        )
    lines.append(
        f"  [{tally.mark(bool(np.isclose(measured_power, 2.0, rtol=0.0, atol=2.0e-12)))}] "
        f"fitted force exponent p: measured={measured_power:.8e} 1; "
        "expected=2.00000000e+00 1"
    )

    outward = np.array([1.0, 0.0])
    parallel_force = _signed_radial_force(_magnetic_force(params, 2.0), outward)
    antiparallel_force = _signed_radial_force(
        _magnetic_force(params, 2.0, (0.0, -1.0)), outward
    )
    transverse = np.linalg.norm(_magnetic_force(params, 2.0, include_scalar=False))
    full = np.linalg.norm(_magnetic_force(params, 2.0, include_scalar=True))
    scalar_magnitude = float(full - transverse)
    measured_ratio = scalar_magnitude / float(transverse)
    expected_ratio = magnetism.scalar_admixture_ratio(
        params.aT, params.aL, params.mu_R, params.B_eff
    )
    lines.extend(
        (
            f"  [{tally.mark(parallel_force < 0.0)}] parallel-current radial sign: "
            f"measured_force={parallel_force:+.8e} force_unit, "
            f"measured_sign={np.sign(parallel_force):+.0f}; expected_sign=-1 (attract)",
            f"  [{tally.mark(antiparallel_force > 0.0)}] antiparallel-current radial sign: "
            f"measured_force={antiparallel_force:+.8e} force_unit, "
            f"measured_sign={np.sign(antiparallel_force):+.0f}; expected_sign=+1 (repel)",
            f"  [{tally.mark(bool(np.isclose(measured_ratio, expected_ratio, rtol=1.0e-13, atol=1.0e-13)))}] "
            f"scalar-admixture contribution: measured_magnitude={scalar_magnitude:.8e} force_unit; "
            f"measured_ratio={measured_ratio:.8e} 1; expected_ratio={expected_ratio:.8e} 1",
        )
    )
    sample_point = np.array([1.0, 0.0])
    positive_field = magnetism.magnetic_field_on_brane(
        sample_point, [[0.0, 0.0]], [1.0], params.N_u, params.aT,
        params.aL, params.mu_R, params.B_eff, True
    )
    negative_field = magnetism.magnetic_field_on_brane(
        sample_point, [[0.0, 0.0]], [-1.0], params.N_u, params.aT,
        params.aL, params.mu_R, params.B_eff, True
    )
    circulation_ok = positive_field[1] > 0.0 and negative_field[1] < 0.0
    lines.append(
        f"  [{tally.mark(circulation_ok)}] magnetic-field circulation at +x: "
        f"forward-current B={_vector(positive_field)} field_unit (expected +y); "
        f"reverse-current B={_vector(negative_field)} field_unit (expected -y)"
    )


def _add_departures(lines: list[str], params: ModelParameters) -> None:
    total_mass = MASS1 + MASS2
    period = gravity.kepler_period(SEMI_MAJOR_AXIS, total_mass, params.G)
    orbital_frequency = 2.0 * np.pi / period
    gravity_departure = gravity.characterized_departure(params.epsilon0)
    monopole = abs(
        gravity.radiation_return_residual(
            0,
            orbital_frequency,
            1.0,
            params.epsilon0,
            params.radiation_length_a,
            params.c_s,
        )
    )
    dipole = abs(
        gravity.radiation_return_residual(
            1,
            orbital_frequency,
            1.0,
            params.epsilon1,
            params.radiation_length_a,
            params.c_s,
        )
    )

    light_departure = light.characterized_departure(
        params.rho_br,
        params.rho_B0,
        params.chi_c,
        params.J_phase,
        params.kappa_phase,
        params.longitudinal_display_fraction,
    )

    charge_departure = charge.characterized_departure(params.ell, params.yukawa_fraction)
    departure_distance = params.ell
    coulomb_force = np.linalg.norm(_charge_force(params, departure_distance, 1.0))
    corrected_force = np.linalg.norm(
        charge.force_on_first(
            [departure_distance, 0.0],
            [0.0, 0.0],
            1.0,
            1.0,
            params.Q_E,
            params.ell,
            params.mouth_half_width_b,
            include_yukawa=True,
            yukawa_fraction=params.yukawa_fraction,
        )
    )
    yukawa_force = float(corrected_force - coulomb_force)

    magnetic_departure = magnetism.characterized_departure(
        params.aT,
        params.aL,
        params.mu_R,
        params.B_eff,
        params.c_E,
        params.c_gamma,
    )
    transverse = np.linalg.norm(_magnetic_force(params, 2.0, include_scalar=False))
    full = np.linalg.norm(_magnetic_force(params, 2.0, include_scalar=True))

    lines.extend(
        (
            "",
            "CHARACTERIZED DEPARTURES",
            f"  gravity [{gravity_departure.code}] evaluation omega={orbital_frequency:.8e} "
            "1/time_unit, source_moment=1.00000000e+00:",
            "    [DERIVED-FORM] monopole residual law proportional to "
            "a*(omega/c_s)*M0*epsilon0/(1+epsilon0)",
            f"    [CALIBRATED-MAGNITUDE] monopole residual amplitude={monopole:.8e} amplitude_unit",
            "    [DERIVED-FORM] dipole residual law proportional to "
            "a^3*(omega/c_s)^3*D1*epsilon1/[2*(1+epsilon1)]",
            f"    [CALIBRATED-MAGNITUDE] dipole residual amplitude={dipole:.8e} amplitude_unit",
            f"    [CALIBRATED-MAGNITUDE] bounded drain factor="
            f"{float(gravity_departure.diagnostics['epsilon0/(1+epsilon0)']):.8e} 1",
            f"  light [{light_departure.code}]:",
            f"    [DERIVED-FORM] stray-longitudinal presence="
            f"{float(light_departure.diagnostics['longitudinal_dof']):.8e} mode",
            f"    [DERIVED-FORM] longitudinal speed="
            f"{float(light_departure.diagnostics['longitudinal_speed']):.8e} length_unit/time_unit",
            f"    [CALIBRATED-MAGNITUDE] display fraction="
            f"{float(light_departure.diagnostics['display_fraction']):.8e} 1",
            f"  charge [{charge_departure.code}] at R={departure_distance:.8e} length_unit:",
            f"    [DERIVED-FORM] Yukawa mass gap="
            f"{float(charge_departure.diagnostics['mass_gap']):.8e} 1/length_unit",
            f"    [CALIBRATED-MAGNITUDE] Yukawa force correction={yukawa_force:.8e} force_unit; "
            f"fraction_of_Coulomb={yukawa_force / float(coulomb_force):.8e} 1",
            f"  magnetism [{magnetic_departure.code}] at R=2.00000000e+00 length_unit:",
            f"    [DERIVED-FORM] scalar-current force contribution={float(full - transverse):.8e} "
            "force_unit",
            f"    [CALIBRATED-MAGNITUDE] scalar/transverse ratio="
            f"{float(magnetic_departure.diagnostics['scalar/transverse_ratio']):.8e} 1",
            f"    [CALIBRATED-MAGNITUDE] preferred-frame residual="
            f"{float(magnetic_departure.diagnostics['preferred_frame_residual']):.8e} "
            "length_unit^2/time_unit^2 (zero by calibrated c_E=c_gamma default)",
        )
    )


def build_report(params: ModelParameters = DEFAULT_PARAMS) -> ReportResult:
    """Run all headless measurements and return their deterministic transcript."""

    params.validate()
    tally = CheckTally()
    lines = [
        "FORCE VISUALIZER NUMERIC VERIFICATION REPORT",
        "Physics source: existing force_visualizer.physics core only; rendering: none.",
    ]
    _add_parameter_block(lines, params)
    _add_gravity(lines, tally, params)
    _add_light(lines, tally, params)
    _add_charge(lines, tally, params)
    _add_magnetism(lines, tally, params)
    _add_departures(lines, params)
    lines.extend(("", f"SUMMARY: {tally.passed} checks passed of {tally.total} total."))
    return ReportResult("\n".join(lines) + "\n", tally.passed, tally.total)


def write_report(result: ReportResult, output_path: Path = OUTPUT_PATH) -> Path:
    """Write a report result to the required deterministic output path."""

    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(result.transcript, encoding="utf-8")
    return output_path


def main() -> int:
    """Generate, print, and persist the verification transcript."""

    result = build_report()
    write_report(result)
    print(result.transcript, end="")
    return 0 if result.passed == result.total else 1


if __name__ == "__main__":
    raise SystemExit(main())
