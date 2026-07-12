"""Gravity-sector closed-orbit and 1PN perihelion analytic goldens."""

import numpy as np
import sympy as sp

from software.force_visualizer.params import DEFAULT_PARAMS as P
from software.force_visualizer.physics import gravity


MASS1, MASS2 = 1.0, 0.20
SEMI_MAJOR_AXIS, ECCENTRICITY = 3.0, 0.25


def _orbit(order: float, periods: float, steps_per_period: int = 1800):
    total = MASS1 + MASS2
    position, velocity = gravity.kepler_periapsis_state(
        SEMI_MAJOR_AXIS, ECCENTRICITY, total, P.G
    )
    period = gravity.kepler_period(SEMI_MAJOR_AXIS, total, P.G)
    return gravity.integrate_relative_orbit(
        position,
        velocity,
        MASS1,
        MASS2,
        P.G,
        P.c_gamma,
        period / steps_per_period,
        int(periods * steps_per_period),
        pn_order=order,
    )


def test_zero_pn_orbit_is_closed() -> None:
    orbit = _orbit(0.0, 5.2)
    measured = gravity.measured_precession(orbit.positions)
    assert abs(measured) < 2e-7


def test_1pn_precession_sign_and_magnitude() -> None:
    orbit = _orbit(1.0, 8.2)
    measured = gravity.measured_precession(orbit.positions)
    expected = gravity.analytic_perihelion_precession(
        SEMI_MAJOR_AXIS, ECCENTRICITY, MASS1 + MASS2, P.G, P.c_gamma
    )
    assert measured > 0.0
    np.testing.assert_allclose(measured, expected, rtol=0.035, atol=2e-4)


def test_return_residual_is_bounded_and_tied_to_drain() -> None:
    zero = gravity.radiation_return_residual(0, 0.2, 1.0, 0.0, 0.5, P.c_s)
    finite = gravity.radiation_return_residual(0, 0.2, 1.0, P.epsilon0, 0.5, P.c_s)
    assert zero == 0j
    assert 0.0 < abs(finite) < 0.5 * 0.2 / P.c_s


def test_drain_field_points_toward_each_isolated_mass() -> None:
    points = np.array([[2.0, 1.0], [-1.5, 0.7], [0.4, -2.2]])
    field = gravity.drain_inflow_field(points, [[0.0, 0.0]], [1.3], P.G)
    assert np.all(np.sum(field * points, axis=1) < 0.0)


def test_eih_euler_lagrange_equations_reproduce_relative_acceleration_1pn() -> None:
    """Symbolically reduce the production EIH Lagrangian through O(c^-2)."""

    x1, y1, x2, y2 = sp.symbols("x1 y1 x2 y2", real=True)
    u1, w1, u2, w2 = sp.symbols("u1 w1 u2 w2", real=True)
    mass1, mass2, coupling, pn = sp.symbols(
        "mass1 mass2 coupling pn", positive=True
    )
    coordinates = (x1, y1, x2, y2)
    velocities = (u1, w1, u2, w2)
    delta = (x1 - x2, y1 - y2)
    radius = sp.sqrt(sum(component**2 for component in delta))
    direction = tuple(component / radius for component in delta)

    lagrangian = gravity.eih_lagrangian_expression(
        radius,
        direction,
        (u1, w1),
        (u2, w2),
        mass1,
        mass2,
        coupling,
        1 / sp.sqrt(pn),
    )
    lagrangian_0pn = lagrangian.subs(pn, 0)
    lagrangian_1pn = sp.diff(lagrangian, pn).subs(pn, 0)

    def velocity_independent_el_terms(expression):
        return sp.Matrix(
            [
                sum(
                    sp.diff(expression, velocities[index], coordinate) * velocity
                    for coordinate, velocity in zip(coordinates, velocities)
                )
                - sp.diff(expression, coordinates[index])
                for index in range(4)
            ]
        )

    mass_matrix_0pn = sp.diag(mass1, mass1, mass2, mass2)
    acceleration_0pn = -mass_matrix_0pn.inv() * velocity_independent_el_terms(
        lagrangian_0pn
    )
    mass_matrix_1pn = sp.Matrix(
        4,
        4,
        lambda row, column: sp.diff(
            lagrangian_1pn, velocities[row], velocities[column]
        ),
    )
    acceleration_1pn = -mass_matrix_0pn.inv() * (
        mass_matrix_1pn * acceleration_0pn
        + velocity_independent_el_terms(lagrangian_1pn)
    )

    total_mass = mass1 + mass2
    separation, radial_velocity, transverse_velocity = sp.symbols(
        "separation radial_velocity transverse_velocity", positive=True
    )
    center_of_mass_frame = {
        x1: mass2 * separation / total_mass,
        x2: -mass1 * separation / total_mass,
        y1: 0,
        y2: 0,
        u1: mass2 * radial_velocity / total_mass,
        u2: -mass1 * radial_velocity / total_mass,
        w1: mass2 * transverse_velocity / total_mass,
        w2: -mass1 * transverse_velocity / total_mass,
    }
    derived_relative = sp.Matrix(
        [acceleration_1pn[0] - acceleration_1pn[2], acceleration_1pn[1] - acceleration_1pn[3]]
    ).subs(center_of_mass_frame)

    eta = mass1 * mass2 / total_mass**2
    mu = coupling * total_mass
    speed_squared = radial_velocity**2 + transverse_velocity**2
    expected_relative = mu / separation**2 * sp.Matrix(
        [
            (4 + 2 * eta) * mu / separation
            - (1 + 3 * eta) * speed_squared
            + sp.Rational(3, 2) * eta * radial_velocity**2
            + (4 - 2 * eta) * radial_velocity**2,
            (4 - 2 * eta) * radial_velocity * transverse_velocity,
        ]
    )
    for derived, expected in zip(derived_relative, expected_relative):
        assert sp.simplify(derived - expected) == 0

    sample = {
        mass1: 1.3,
        mass2: 0.4,
        coupling: 0.9,
        separation: 2.2,
        radial_velocity: 0.17,
        transverse_velocity: 0.63,
    }
    symbolic_sample = np.array(
        [float(component.subs(sample)) for component in derived_relative]
    )
    production_sample = gravity.relative_acceleration_1pn(
        [sample[separation], 0.0],
        [sample[radial_velocity], sample[transverse_velocity]],
        sample[mass1],
        sample[mass2],
        sample[coupling],
        1.0,
    )
    np.testing.assert_allclose(symbolic_sample, production_sample, rtol=2e-14, atol=1e-14)
