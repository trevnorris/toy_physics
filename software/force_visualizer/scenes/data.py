"""Render-neutral scenario setup shared by the GIF and live front ends.

This module chooses initial conditions and samples the existing physics APIs.
It deliberately contains neither Matplotlib imports nor force-law equations.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import pi

import numpy as np

from ..params import DEFAULT_PARAMS, ModelParameters
from ..physics import charge, gravity, light, magnetism
from ..physics.integrators import integrate_fixed


@dataclass(frozen=True)
class GravitySceneData:
    newtonian: gravity.Orbit
    relativistic: gravity.Orbit
    body1: np.ndarray
    body2: np.ndarray
    expected_shift: float
    departure: object
    semi_major_axis: float
    eccentricity: float


def build_gravity_data(
    params: ModelParameters = DEFAULT_PARAMS,
    *,
    mass1: float = 0.03,
    mass2: float = 1.0,
    semi_major_axis: float = 3.0,
    eccentricity: float = 0.35,
    periods: float = 3.2,
    steps_per_period: int = 1000,
    include_25pn_benchmark: bool = False,
) -> GravitySceneData:
    """Sample the existing gravity core for the shared orbit scenario."""

    params.validate()
    if mass1 <= 0.0 or mass2 <= 0.0 or periods <= 0.0 or steps_per_period < 20:
        raise ValueError("masses/periods must be positive and steps_per_period >= 20")
    total_mass = mass1 + mass2
    initial_r, initial_v = gravity.kepler_periapsis_state(
        semi_major_axis, eccentricity, total_mass, params.G
    )
    period = gravity.kepler_period(semi_major_axis, total_mass, params.G)
    steps = max(1, int(periods * steps_per_period))
    dt = period / steps_per_period
    newtonian = gravity.integrate_relative_orbit(
        initial_r, initial_v, mass1, mass2, params.G, params.c_gamma, dt, steps, pn_order=0.0
    )
    pn_order = 2.5 if include_25pn_benchmark else 1.0
    relativistic = gravity.integrate_relative_orbit(
        initial_r,
        initial_v,
        mass1,
        mass2,
        params.G,
        params.c_gamma,
        dt,
        steps,
        pn_order=pn_order,
        radiation_reaction_scale=params.radiation_reaction_benchmark_scale,
    )
    body1, body2 = relativistic.body_positions()
    return GravitySceneData(
        newtonian=newtonian,
        relativistic=relativistic,
        body1=body1,
        body2=body2,
        expected_shift=gravity.analytic_perihelion_precession(
            semi_major_axis, eccentricity, total_mass, params.G, params.c_gamma
        ),
        departure=gravity.characterized_departure(params.epsilon0),
        semi_major_axis=semi_major_axis,
        eccentricity=eccentricity,
    )


@dataclass(frozen=True)
class LightSceneData:
    first_evolution: light.TransverseWaveEvolution
    second_evolution: light.TransverseWaveEvolution
    ray: light.RayTrace
    departure: object
    longitudinal_speed: float
    wavenumber: float
    c_gamma: float


def build_light_data(
    params: ModelParameters = DEFAULT_PARAMS,
    *,
    wavenumber: float = 2.0,
    c_gamma: float | None = None,
    spatial_points: int = 300,
    periods: float = 1.0,
    ray_steps: int = 2500,
) -> LightSceneData:
    """Sample the existing transverse-wave and ray-tracing core APIs."""

    params.validate()
    speed = params.c_gamma if c_gamma is None else float(c_gamma)
    first = light.evolve_transverse_wave(
        wavenumber, speed, spatial_points=spatial_points, periods=periods
    )
    second = light.evolve_transverse_wave(
        wavenumber, speed, spatial_points=spatial_points, periods=periods, phase=pi / 2.0
    )
    ray = light.trace_ray(
        impact_parameter=2.5,
        z_extent=20.0,
        mass=1.0,
        G=params.G,
        c_gamma=speed,
        steps=ray_steps,
    )
    departure = light.characterized_departure(
        params.rho_br,
        params.rho_B0,
        params.chi_c,
        params.J_phase,
        params.kappa_phase,
        params.longitudinal_display_fraction,
    )
    return LightSceneData(
        first_evolution=first,
        second_evolution=second,
        ray=ray,
        departure=departure,
        longitudinal_speed=float(departure.diagnostics["longitudinal_speed"]),
        wavenumber=float(wavenumber),
        c_gamma=speed,
    )


@dataclass(frozen=True)
class ChargeSceneData:
    trajectories: charge.ChargeTrajectories
    departure: object
    distances: np.ndarray
    coulomb_force: np.ndarray
    corrected_force: np.ndarray


def build_charge_data(
    params: ModelParameters = DEFAULT_PARAMS,
    *,
    orientations: tuple[float, float] = (1.0, -1.0),
    Q_E: float | None = None,
    show_departure: bool = True,
    dt: float = 0.003,
    steps: int = 420,
) -> ChargeSceneData:
    """Sample the existing localized Coulomb/Yukawa trajectory core."""

    params.validate()
    coupling_anchor = params.Q_E if Q_E is None else float(Q_E)
    trajectories = charge.integrate_charges(
        positions=[[-0.75, 0.0], [0.75, 0.0]],
        velocities=[[0.0, 0.0], [0.0, 0.0]],
        masses=[1.0, 1.0],
        orientations=orientations,
        Q_E=coupling_anchor,
        ell=params.ell,
        mouth_half_width=params.mouth_half_width_b,
        dt=dt,
        steps=steps,
        include_yukawa=show_departure,
        yukawa_fraction=params.yukawa_fraction,
    )
    departure = charge.characterized_departure(params.ell, params.yukawa_fraction)
    distances = np.geomspace(0.30, 5.0, 160)
    coupling = charge.electric_coupling(
        coupling_anchor, params.ell, params.mouth_half_width_b
    )
    coulomb_force = coupling / distances**2
    gap = float(departure.diagnostics["mass_gap"])
    corrected_force = coulomb_force * (
        1.0 + params.yukawa_fraction * np.exp(-gap * distances) * (1.0 + gap * distances)
    )
    return ChargeSceneData(
        trajectories=trajectories,
        departure=departure,
        distances=distances,
        coulomb_force=coulomb_force,
        corrected_force=corrected_force,
    )


@dataclass(frozen=True)
class MagnetismSceneData:
    trajectories: magnetism.CurrentTrajectories
    current_directions: tuple[float, float]
    departure: object
    scalar_ratio: float


def following_transverse_grid(
    source_positions: np.ndarray,
    *,
    columns: int = 17,
    rows: int = 11,
    padding: float = 0.55,
    minimum_half_height: float = 0.75,
) -> np.ndarray:
    """Return a transverse field grid whose extent follows all current markers."""

    sources = np.asarray(source_positions, dtype=float)
    if sources.ndim != 2 or sources.shape[1] != 2:
        raise ValueError("source_positions must have shape (N, 2)")
    if columns < 2 or rows < 2 or padding <= 0.0 or minimum_half_height <= 0.0:
        raise ValueError("grid dimensions and extents must be positive")
    center = sources.mean(axis=0)
    x_min = min(float(sources[:, 0].min() - padding), center[0] - padding)
    x_max = max(float(sources[:, 0].max() + padding), center[0] + padding)
    half_height = max(
        minimum_half_height,
        0.5 * float(np.ptp(sources[:, 1])) + padding,
    )
    grid_x, grid_y = np.meshgrid(
        np.linspace(x_min, x_max, columns),
        np.linspace(center[1] - half_height, center[1] + half_height, rows),
    )
    return np.stack((grid_x, grid_y), axis=-1)


def build_magnetism_data(
    params: ModelParameters = DEFAULT_PARAMS,
    *,
    current_directions: tuple[float, float] = (1.0, 1.0),
    aL: float | None = None,
    show_departure: bool = True,
    dt: float = 0.06,
    steps: int = 1500,
) -> MagnetismSceneData:
    """Integrate core-driven transverse drift for two end-on currents.

    The current vectors are fixed normal to the displayed transverse plane.
    Only the transverse response to :func:`magnetism.force_on_second` is
    integrated, so current direction is not mistaken for an in-plane body
    velocity.
    """

    params.validate()
    directions = np.asarray(current_directions, dtype=float)
    if directions.shape != (2,) or np.any(directions == 0.0):
        raise ValueError("current_directions must contain two nonzero values")
    if dt <= 0.0 or steps < 1:
        raise ValueError("dt and steps must be positive")
    longitudinal_amplitude = params.aL if aL is None else float(aL)
    initial_positions = np.array([[-1.30, 0.0], [1.30, 0.0]])
    initial_velocities = np.zeros((2, 2))
    initial_state = np.concatenate((initial_positions.ravel(), initial_velocities.ravel()))
    end_on_currents = np.column_stack((np.zeros((2, 2)), directions))

    def transverse_rhs(_time: float, state: np.ndarray) -> np.ndarray:
        positions = state[:4].reshape(2, 2)
        velocities = state[4:].reshape(2, 2)
        positions_3d = np.column_stack((positions, np.zeros(2)))
        force_on_second = magnetism.force_on_second(
            positions_3d[0],
            positions_3d[1],
            end_on_currents[0],
            end_on_currents[1],
            1.0,
            1.0,
            params.N_u,
            params.aT,
            longitudinal_amplitude,
            params.mu_R,
            params.B_eff,
            show_departure,
        )[:2]
        accelerations = np.stack((-force_on_second, force_on_second))
        return np.concatenate((velocities.ravel(), accelerations.ravel()))

    times, states = integrate_fixed(transverse_rhs, initial_state, dt, steps)
    trajectories = magnetism.CurrentTrajectories(
        times=times,
        positions=states[:, :4].reshape(steps + 1, 2, 2),
        velocities=states[:, 4:].reshape(steps + 1, 2, 2),
    )
    departure = magnetism.characterized_departure(
        params.aT,
        longitudinal_amplitude,
        params.mu_R,
        params.B_eff,
        params.c_E,
        params.c_gamma,
    )
    return MagnetismSceneData(
        trajectories=trajectories,
        current_directions=(float(directions[0]), float(directions[1])),
        departure=departure,
        scalar_ratio=float(departure.diagnostics["scalar/transverse_ratio"]),
    )
