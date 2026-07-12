"""Charge-sector analytic goldens from build spec section 6."""

import numpy as np

from software.force_visualizer.params import DEFAULT_PARAMS as P
from software.force_visualizer.physics import charge


def _force_at(distance: float, second_orientation: float) -> np.ndarray:
    return charge.force_on_first(
        [distance, 0.0],
        [0.0, 0.0],
        1.0,
        second_orientation,
        P.Q_E,
        P.ell,
        P.mouth_half_width_b,
    )


def test_like_repels_and_unlike_attracts() -> None:
    assert _force_at(1.0, 1.0)[0] > 0.0
    assert _force_at(1.0, -1.0)[0] < 0.0


def test_measured_coulomb_force_exponent_is_two() -> None:
    distances = np.geomspace(0.7, 8.0, 40)
    forces = np.array([np.linalg.norm(_force_at(distance, 1.0)) for distance in distances])
    slope = np.polyfit(np.log(distances), np.log(forces), 1)[0]
    measured_power = -slope
    np.testing.assert_allclose(measured_power, 2.0, rtol=0.0, atol=1e-12)


def test_yukawa_is_short_range_and_toggleable() -> None:
    near = charge.yukawa_potential(0.5, 1.0, 1.0, P.Q_E, P.ell, P.mouth_half_width_b, P.yukawa_fraction)
    far = charge.yukawa_potential(5.0, 1.0, 1.0, P.Q_E, P.ell, P.mouth_half_width_b, P.yukawa_fraction)
    assert near > 0.0
    assert far / near < 1e-5


def test_electric_field_points_away_from_positive_and_toward_negative() -> None:
    point = np.array([1.2, -0.4])
    positive = charge.electric_field(
        point, [[0.0, 0.0]], [1.0], P.Q_E, P.ell, P.mouth_half_width_b
    )
    negative = charge.electric_field(
        point, [[0.0, 0.0]], [-1.0], P.Q_E, P.ell, P.mouth_half_width_b
    )
    assert np.dot(positive, point) > 0.0
    assert np.dot(negative, point) < 0.0
