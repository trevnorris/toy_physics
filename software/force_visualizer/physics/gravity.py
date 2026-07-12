"""Newtonian, 1PN, and optional benchmark 2.5PN gravity laws.

Sources
-------
* ``software/stage1_solver/reports/pathA_29_brane_bulk_return.md`` and
  ``pathA_29_results.yaml``: localized ``1/r`` Green function, ``1/r^2``
  attractive flow, and the bounded monopole/dipole return residual.
* ``research/4d_1pn_full/paper/4d_1pn_full.tex``, Secs. "Full two-body 1PN
  Lagrangian" and "Test-mass orbit reduction": EIH form and perihelion law.
* ``research/4d_2_5pn/paper/4d_2_5pn.tex``, Sec. "Burke--Thorne prototype":
  conditional local 2.5PN benchmark.  Its native toy normalization remains
  ``GENUINE_BLOCKED`` and is never represented here as derived.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import pi, sqrt
from typing import Tuple

import numpy as np
from numpy.typing import ArrayLike, NDArray

from .departures import CharacterizedDeparture
from .integrators import integrate_fixed

Vector = NDArray[np.float64]


def _separation(vector: ArrayLike) -> Tuple[Vector, float, Vector]:
    r_vec = np.asarray(vector, dtype=float)
    radius = float(np.linalg.norm(r_vec))
    if radius <= 0.0:
        raise ValueError("coincident bodies make the point-force law singular")
    return r_vec, radius, r_vec / radius


def newtonian_relative_acceleration(r: ArrayLike, total_mass: float, G: float) -> Vector:
    """Relative acceleration ``-G (m1+m2) r/r^3``.

    This is the point-law reduction of the ``1/r`` localized zero-mode Green
    function in pathA_29 ``zero_mode`` / ``static_falloff_B`` (derived form);
    ``G`` is an effective calibrated magnitude.
    Source: ``software/stage1_solver/reports/pathA_29_results.yaml``, result
    keys ``zero_mode.flow_from_gradient`` and ``static_falloff_B.p``.
    """

    r_vec, radius, _ = _separation(r)
    return -(G * total_mass / radius**3) * r_vec


def drain_inflow_field(
    sample_points: ArrayLike,
    source_positions: ArrayLike,
    source_masses: ArrayLike,
    G: float,
) -> Vector:
    """Sample the Newtonian core as the gravity-drain inflow field.

    Every returned vector is a sum of calls to
    :func:`newtonian_relative_acceleration`; no separate field or force law is
    introduced.  The renderer uses its direction for the in-brane medium
    inflow requested by ``software/force_visualizer/notes/build_spec.md`` §10
    and normalizes arrow/streamline speeds for display.  Thus the display does
    not claim a separately calibrated drain-velocity magnitude.

    The point-law source and inward sign are those of
    ``pathA_29_results.yaml:zero_mode.flow_from_gradient`` and
    ``static_falloff_B.p``.
    """

    points = np.asarray(sample_points, dtype=float)
    sources = np.asarray(source_positions, dtype=float)
    masses = np.asarray(source_masses, dtype=float)
    if points.ndim < 1 or points.shape[-1] != 2:
        raise ValueError("sample_points must end in two spatial coordinates")
    if sources.ndim != 2 or sources.shape[1] != 2 or masses.shape != (len(sources),):
        raise ValueError("source_positions must be (N,2) and source_masses must be (N,)")
    if np.any(masses <= 0.0):
        raise ValueError("source masses must be positive")

    flat_points = points.reshape(-1, 2)
    field = np.zeros_like(flat_points)
    for index, point in enumerate(flat_points):
        for position, mass in zip(sources, masses):
            if np.array_equal(point, position):
                continue  # singular source pixel is deliberately left blank
            field[index] += newtonian_relative_acceleration(point - position, mass, G)
    return field.reshape(points.shape)


def relative_acceleration_1pn(
    r: ArrayLike,
    v: ArrayLike,
    mass1: float,
    mass2: float,
    G: float,
    c: float,
) -> Vector:
    """Conservative EIH relative acceleration through first PN order.

    The coefficients are the equations-of-motion reduction of the exact EIH
    Lagrangian in ``4d_1pn_full.tex``, Eq. ``derived-eih``.  ``r=x1-x2`` and
    ``v=dr/dt``.  No calibrated coefficient modifies the derived 1PN form.
    Source path: ``research/4d_1pn_full/paper/4d_1pn_full.tex``, section
    ``Full two-body 1PN Lagrangian and exact EIH match``.
    """

    r_vec, radius, n = _separation(r)
    velocity = np.asarray(v, dtype=float)
    total_mass = mass1 + mass2
    if total_mass <= 0.0 or c <= 0.0:
        raise ValueError("masses and c must be positive")
    eta = mass1 * mass2 / total_mass**2
    mu = G * total_mass
    speed2 = float(np.dot(velocity, velocity))
    radial_speed = float(np.dot(n, velocity))
    correction_n = (
        (4.0 + 2.0 * eta) * mu / radius
        - (1.0 + 3.0 * eta) * speed2
        + 1.5 * eta * radial_speed**2
    )
    correction_v = (4.0 - 2.0 * eta) * radial_speed
    return (mu / (c**2 * radius**2)) * (correction_n * n + correction_v * velocity)


def relative_acceleration_25pn_benchmark(
    r: ArrayLike,
    v: ArrayLike,
    mass1: float,
    mass2: float,
    G: float,
    c: float,
    benchmark_scale: float = 1.0,
) -> Vector:
    """Burke--Thorne member of the local 2.5PN reaction family.

    Source: ``4d_2_5pn.tex``, "Burke--Thorne prototype", with
    ``alpha=4, beta=5``.  ``benchmark_scale`` is explicitly a calibrated
    display multiplier because the cited paper leaves the native response
    normalization ``mhat0^2 Gamma5`` genuinely blocked.  A value of one means
    the standard benchmark, not a toy-model derivation.
    Full path: ``research/4d_2_5pn/paper/4d_2_5pn.tex``, equations for
    ``A``, ``B`` and the ``alpha=4,beta=5`` member.
    """

    _, radius, n = _separation(r)
    velocity = np.asarray(v, dtype=float)
    total_mass = mass1 + mass2
    if total_mass <= 0.0 or c <= 0.0 or benchmark_scale < 0.0:
        raise ValueError("masses/c must be positive and benchmark_scale non-negative")
    eta = mass1 * mass2 / total_mass**2
    mu = G * total_mass
    speed2 = float(np.dot(velocity, velocity))
    radial_speed = float(np.dot(n, velocity))
    coefficient_a = 18.0 * speed2 + 2.0 * mu / (3.0 * radius) - 25.0 * radial_speed**2
    coefficient_b = -6.0 * speed2 + 2.0 * mu / radius + 15.0 * radial_speed**2
    prefactor = benchmark_scale * 8.0 * G**2 * total_mass**2 * eta / (5.0 * c**5 * radius**3)
    return prefactor * (coefficient_a * radial_speed * n + coefficient_b * velocity)


def relative_acceleration(
    r: ArrayLike,
    v: ArrayLike,
    mass1: float,
    mass2: float,
    G: float,
    c: float,
    pn_order: float = 0.0,
    radiation_reaction_scale: float = 1.0,
) -> Vector:
    """Assemble the Newtonian/1PN/optional benchmark-2.5PN ladder.

    Sources: ``pathA_29_results.yaml:static_falloff_B``;
    ``research/4d_1pn_full/paper/4d_1pn_full.tex:derived-eih``; and
    ``research/4d_2_5pn/paper/4d_2_5pn.tex:Burke--Thorne prototype``.
    """

    acceleration = newtonian_relative_acceleration(r, mass1 + mass2, G)
    if pn_order >= 1.0:
        acceleration = acceleration + relative_acceleration_1pn(r, v, mass1, mass2, G, c)
    if pn_order >= 2.5:
        acceleration = acceleration + relative_acceleration_25pn_benchmark(
            r, v, mass1, mass2, G, c, radiation_reaction_scale
        )
    return acceleration


def eih_lagrangian(
    x1: ArrayLike,
    x2: ArrayLike,
    v1: ArrayLike,
    v2: ArrayLike,
    mass1: float,
    mass2: float,
    G: float,
    c: float,
) -> float:
    """Exact conservative two-body 1PN EIH Lagrangian from the cited paper.

    Source: ``research/4d_1pn_full/paper/4d_1pn_full.tex``, equation
    ``derived-eih`` in section ``Full two-body 1PN Lagrangian``.
    """

    r_vec, radius, n = _separation(np.asarray(x1, dtype=float) - np.asarray(x2, dtype=float))
    del r_vec
    velocity1 = np.asarray(v1, dtype=float)
    velocity2 = np.asarray(v2, dtype=float)
    return float(
        eih_lagrangian_expression(
            radius, n, velocity1, velocity2, mass1, mass2, G, c
        )
    )


def eih_lagrangian_expression(
    radius,
    separation_direction,
    velocity1,
    velocity2,
    mass1,
    mass2,
    G,
    c,
):
    """Backend-neutral EIH expression used by numeric code and symbolic tests.

    Inputs may be ordinary numbers/arrays or symbolic scalar sequences.  This
    keeps the Euler--Lagrange consistency test tied to the exact production
    Lagrangian coefficients instead of maintaining a test-only transcription.
    """

    def dot(left, right):
        return sum(first * second for first, second in zip(left, right))

    v1_sq = dot(velocity1, velocity1)
    v2_sq = dot(velocity2, velocity2)
    cross = dot(velocity1, velocity2)
    radial_cross = dot(velocity1, separation_direction) * dot(
        velocity2, separation_direction
    )
    newtonian = mass1 * v1_sq / 2 + mass2 * v2_sq / 2 + G * mass1 * mass2 / radius
    kinetic = (mass1 * v1_sq**2 + mass2 * v2_sq**2) / (8 * c**2)
    pair = (G * mass1 * mass2 / (c**2 * radius)) * (
        3 * (v1_sq + v2_sq) / 2 - 7 * cross / 2 - radial_cross / 2
    )
    static = -(G**2 * mass1 * mass2 * (mass1 + mass2)) / (2 * c**2 * radius**2)
    return newtonian + kinetic + pair + static


def kepler_periapsis_state(
    semi_major_axis: float, eccentricity: float, total_mass: float, G: float
) -> Tuple[Vector, Vector]:
    """Newtonian periapsis initial data in the relative center-of-mass frame.

    Observable source: ``research/4d_1pn_full/paper/4d_1pn_full.tex``, section
    ``Test-mass orbit reduction``, relation ``ell^2=mu*a*(1-e^2)``.
    """

    if semi_major_axis <= 0.0 or total_mass <= 0.0 or not 0.0 <= eccentricity < 1.0:
        raise ValueError("require a>0, M>0 and 0<=e<1")
    periapsis = semi_major_axis * (1.0 - eccentricity)
    speed = sqrt(G * total_mass * (1.0 + eccentricity) / periapsis)
    return np.array([periapsis, 0.0]), np.array([0.0, speed])


def kepler_period(semi_major_axis: float, total_mass: float, G: float) -> float:
    """Newtonian period used only to choose integration/reporting windows.

    Law source: ``software/stage1_solver/reports/pathA_29_results.yaml``, key
    ``static_falloff_B.p=2``; this helper introduces no new force term.
    """

    return 2.0 * pi * sqrt(semi_major_axis**3 / (G * total_mass))


@dataclass(frozen=True)
class Orbit:
    """Deterministic relative trajectory and its simulation metadata.

    Architecture source: ``software/force_visualizer/notes/build_spec.md`` §2.
    """

    times: Vector
    positions: Vector
    velocities: Vector
    mass1: float
    mass2: float
    pn_order: float

    def body_positions(self) -> Tuple[Vector, Vector]:
        """Return center-of-mass-frame body tracks from the relative track.

        Convention source: ``research/4d_2_5pn/paper/4d_2_5pn.tex``, section
        ``Burke--Thorne prototype`` center-of-mass/source decomposition.
        """

        total = self.mass1 + self.mass2
        return (self.mass2 / total) * self.positions, -(self.mass1 / total) * self.positions


def integrate_relative_orbit(
    initial_position: ArrayLike,
    initial_velocity: ArrayLike,
    mass1: float,
    mass2: float,
    G: float,
    c: float,
    dt: float,
    steps: int,
    pn_order: float = 0.0,
    radiation_reaction_scale: float = 1.0,
) -> Orbit:
    """Integrate the render-agnostic relative two-body equations with RK4.

    Law sources: ``pathA_29_results.yaml:static_falloff_B`` and
    ``research/4d_1pn_full/paper/4d_1pn_full.tex:derived-eih``; deterministic
    integration contract: ``build_spec.md`` §§2,6.
    """

    position = np.asarray(initial_position, dtype=float)
    velocity = np.asarray(initial_velocity, dtype=float)
    if position.shape != velocity.shape or position.ndim != 1:
        raise ValueError("position and velocity must be same-shape vectors")
    dimension = position.size
    initial = np.concatenate((position, velocity))

    def rhs(_time: float, state: Vector) -> Vector:
        r = state[:dimension]
        v = state[dimension:]
        a = relative_acceleration(
            r,
            v,
            mass1,
            mass2,
            G,
            c,
            pn_order=pn_order,
            radiation_reaction_scale=radiation_reaction_scale,
        )
        return np.concatenate((v, a))

    times, states = integrate_fixed(rhs, initial, dt, steps)
    return Orbit(times, states[:, :dimension], states[:, dimension:], mass1, mass2, pn_order)


def periapsis_angles(positions: ArrayLike) -> Vector:
    """Find unwrapped angles at strict sampled radial minima.

    Acceptance source: ``software/force_visualizer/notes/build_spec.md`` §6,
    gravity perihelion measurement.
    """

    points = np.asarray(positions, dtype=float)
    if points.ndim != 2 or points.shape[1] != 2 or len(points) < 3:
        raise ValueError("positions must be an (N,2) array")
    radii = np.linalg.norm(points, axis=1)
    angles = np.unwrap(np.arctan2(points[:, 1], points[:, 0]))
    minima = np.flatnonzero((radii[1:-1] < radii[:-2]) & (radii[1:-1] <= radii[2:])) + 1
    refined = []
    for index in minima:
        y0, y1, y2 = radii[index - 1 : index + 2]
        denominator = y0 - 2.0 * y1 + y2
        fraction = 0.0 if denominator == 0.0 else 0.5 * (y0 - y2) / denominator
        local_angle = angles[index] + fraction * 0.5 * (angles[index + 1] - angles[index - 1])
        refined.append(local_angle)
    return np.asarray(refined, dtype=float)


def measured_precession(positions: ArrayLike) -> float:
    """Mean periapsis advance per radial cycle, subtracting ``2*pi``.

    Observable source: ``research/4d_1pn_full/paper/4d_1pn_full.tex``, section
    ``Perihelion coefficient``, ``Phi_r-2*pi``.
    """

    angles = periapsis_angles(positions)
    if len(angles) < 2:
        raise ValueError("trajectory contains fewer than two sampled periapses")
    return float(np.mean(np.diff(angles) - 2.0 * pi))


def analytic_perihelion_precession(
    semi_major_axis: float, eccentricity: float, total_mass: float, G: float, c: float
) -> float:
    """Derived EIH/GR 1PN shift ``6*pi*GM/(c^2*a*(1-e^2))``.

    Source: ``research/4d_1pn_full/paper/4d_1pn_full.tex``, equation
    ``deltaphi-model-beta3`` in section ``Perihelion coefficient``.
    """

    return 6.0 * pi * G * total_mass / (c**2 * semi_major_axis * (1.0 - eccentricity**2))


def radiation_return_residual(
    multipole: int,
    omega: float,
    source_moment: complex,
    epsilon: float,
    length_scale: float,
    c_s: float,
) -> complex:
    """Leading pathA_29 monopole/dipole return residual amplitude.

    ``multipole=0`` has power one and ``multipole=1`` power three.  The
    bounded magnitude factor ``epsilon/(1+epsilon)`` ties the falsifiable
    residual to the same drain strength as the static gravity source.
    Source: ``software/stage1_solver/reports/pathA_29_results.yaml``, keys
    ``residual_prediction.scaling`` and ``epsilon_to_gravity_strength_tie``.
    """

    if epsilon < 0.0 or c_s <= 0.0 or length_scale <= 0.0:
        raise ValueError("epsilon>=0 and positive length/speed are required")
    bounded = epsilon / (1.0 + epsilon)
    if multipole == 0:
        return 1j * length_scale * (omega / c_s) * source_moment * bounded
    if multipole == 1:
        return 0.5j * length_scale**3 * (omega / c_s) ** 3 * source_moment * bounded
    raise ValueError("pathA_29 characterizes only multipoles ell=0 and ell=1 here")


def characterized_departure(epsilon0: float) -> CharacterizedDeparture:
    """Return the mandatory honest-scope gravity departure record.

    Source: ``software/stage1_solver/reports/pathA_29_brane_bulk_return.md``
    headline ``RETURN_RESIDUAL_PREDICTION`` and results ``residual_prediction``.
    """

    bounded = epsilon0 / (1.0 + epsilon0)
    return CharacterizedDeparture(
        code="RETURN_RESIDUAL_PREDICTION",
        description="bounded monopole/dipole c_s-radiation, absent from GR gravitational waves",
        derived_form="A0 proportional to i*a*(omega/c_s)*M0*epsilon0/(1+epsilon0); p0=1, p1=3",
        calibrated_magnitude="epsilon0 and the effective drain/radiation scales are calibration inputs",
        diagnostics={"epsilon0/(1+epsilon0)": bounded, "monopole_power": 1.0, "dipole_power": 3.0},
    )
