"""Magnetism-sector sign and falloff analytic goldens."""

import numpy as np

from software.force_visualizer.params import DEFAULT_PARAMS as P
from software.force_visualizer.physics import magnetism
from software.force_visualizer.scenes.data import (
    build_magnetism_data,
    following_transverse_grid,
)


def _parallel_force(distance: float, include_scalar: bool = True) -> np.ndarray:
    return magnetism.force_on_second(
        [-distance / 2.0, 0.0],
        [distance / 2.0, 0.0],
        [0.0, 1.0],
        [0.0, 1.0],
        1.0,
        1.0,
        P.N_u,
        P.aT,
        P.aL,
        P.mu_R,
        P.B_eff,
        include_scalar,
    )


def test_parallel_like_currents_attract_in_both_channels() -> None:
    transverse_only = _parallel_force(2.0, include_scalar=False)
    full = _parallel_force(2.0, include_scalar=True)
    assert transverse_only[0] < 0.0  # body 2 at +x is pulled inward
    assert full[0] < transverse_only[0]  # scalar admixture is also attractive


def test_calibrated_scalar_admixture_is_visible_but_subdominant() -> None:
    ratio = magnetism.scalar_admixture_ratio(P.aT, P.aL, P.mu_R, P.B_eff)
    assert 0.1 <= ratio <= 0.3
    np.testing.assert_allclose(ratio, 0.2, rtol=0.0, atol=1e-15)


def test_measured_magnetic_force_exponent_is_two() -> None:
    distances = np.geomspace(0.8, 8.0, 40)
    forces = np.array([np.linalg.norm(_parallel_force(distance)) for distance in distances])
    measured_power = -np.polyfit(np.log(distances), np.log(forces), 1)[0]
    np.testing.assert_allclose(measured_power, 2.0, rtol=0.0, atol=2e-12)


def test_force_is_negative_gradient_of_report_potential() -> None:
    distance = 2.1
    delta = 1e-6
    plus = magnetism.current_potential(
        [0.0, 0.0], [distance + delta, 0.0], [0.0, 1.0], [0.0, 1.0],
        1.0, 1.0, P.N_u, P.aT, P.aL, P.mu_R, P.B_eff, True
    )
    minus = magnetism.current_potential(
        [0.0, 0.0], [distance - delta, 0.0], [0.0, 1.0], [0.0, 1.0],
        1.0, 1.0, P.N_u, P.aT, P.aL, P.mu_R, P.B_eff, True
    )
    numerical_force = -(plus - minus) / (2.0 * delta)
    np.testing.assert_allclose(_parallel_force(distance)[0], numerical_force, rtol=1e-9)


def test_brane_field_circulation_reverses_with_current_direction() -> None:
    points = np.array([[1.0, 0.0], [0.0, 1.0]])
    forward = magnetism.magnetic_field_on_brane(
        points, [[0.0, 0.0]], [1.0], P.N_u, P.aT, P.aL, P.mu_R, P.B_eff
    )
    reverse = magnetism.magnetic_field_on_brane(
        points, [[0.0, 0.0]], [-1.0], P.N_u, P.aT, P.aL, P.mu_R, P.B_eff
    )
    np.testing.assert_allclose(reverse, -forward, rtol=0.0, atol=0.0)
    assert forward[0, 1] > 0.0  # +x maps to +y for forward circulation
    assert forward[1, 0] < 0.0  # +y maps to -x


def test_end_on_scene_motion_is_transverse_core_attraction_or_repulsion() -> None:
    parallel = build_magnetism_data(P, current_directions=(1.0, 1.0), dt=0.25, steps=360)
    antiparallel = build_magnetism_data(
        P, current_directions=(1.0, -1.0), dt=0.25, steps=360
    )

    def separations(data) -> np.ndarray:
        positions = data.trajectories.positions
        return np.linalg.norm(positions[:, 1] - positions[:, 0], axis=1)

    parallel_separation = separations(parallel)
    antiparallel_separation = separations(antiparallel)
    assert parallel_separation[-1] < parallel_separation[0]
    assert antiparallel_separation[-1] > antiparallel_separation[0]

    # No current-streaming coordinate is integrated in the screen plane.
    np.testing.assert_allclose(parallel.trajectories.positions[:, :, 1], 0.0, atol=0.0)
    np.testing.assert_allclose(antiparallel.trajectories.positions[:, :, 1], 0.0, atol=0.0)


def test_following_field_grid_always_encloses_end_on_currents() -> None:
    data = build_magnetism_data(P, current_directions=(1.0, -1.0), dt=0.25, steps=360)
    for positions in data.trajectories.positions[::60]:
        grid = following_transverse_grid(positions)
        assert grid[..., 0].min() < positions[:, 0].min()
        assert grid[..., 0].max() > positions[:, 0].max()
        assert grid[..., 1].min() < positions[:, 1].min()
        assert grid[..., 1].max() > positions[:, 1].max()
