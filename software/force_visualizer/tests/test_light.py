"""Light-sector dispersion, polarization, and lensing analytic goldens."""

import numpy as np

from software.force_visualizer.params import DEFAULT_PARAMS as P
from software.force_visualizer.physics import light


def test_measured_dispersion_is_linear_with_c_gamma_slope() -> None:
    wavenumbers = np.array([1.0, 2.0, 3.0, 4.0])
    measured_omega = np.array(
        [
            light.measured_transverse_angular_frequency(
                light.evolve_transverse_wave(
                    wavenumber, P.c_gamma, spatial_points=384, periods=4.0
                )
            )
            for wavenumber in wavenumbers
        ]
    )
    slope, intercept = np.polyfit(wavenumbers, measured_omega, 1)
    np.testing.assert_allclose(slope, P.c_gamma, rtol=4e-4, atol=1e-6)
    assert abs(intercept) < 4e-3
    np.testing.assert_allclose(measured_omega, P.c_gamma * wavenumbers, rtol=4e-4)


def test_exactly_two_transverse_polarizations() -> None:
    wavevector = np.array([1.2, -0.7, 2.3])
    polarizations = light.transverse_polarizations(wavevector)
    assert polarizations.shape == (2, 3)
    np.testing.assert_allclose(polarizations @ wavevector, np.zeros(2), atol=1e-14)
    np.testing.assert_allclose(polarizations @ polarizations.T, np.eye(2), atol=1e-14)


def test_ray_bends_toward_mass_with_correct_weak_field_magnitude() -> None:
    impact = 2.5
    ray = light.trace_ray(impact, 20.0, 1.0, P.G, P.c_gamma, steps=3000)
    expected = light.analytic_signed_deflection(impact, 1.0, P.G, P.c_gamma)
    assert ray.signed_deflection < 0.0  # positive-x ray bends toward x=0
    np.testing.assert_allclose(ray.signed_deflection, expected, rtol=0.02, atol=1e-5)


def test_stray_longitudinal_mode_is_real_not_counted_as_transverse() -> None:
    analysis = light.analyze_longitudinal_mode(
        1.0, P.rho_br, P.rho_B0, P.chi_c, P.J_phase, P.kappa_phase
    )
    speed = light.stray_longitudinal_speed(
        P.rho_br, P.rho_B0, P.chi_c, P.J_phase, P.kappa_phase
    )
    positive_roots = [
        root.real
        for root in analysis.angular_frequency_roots
        if abs(root.imag) < 1e-12 and root.real > 0.0
    ]
    assert analysis.propagating_dof == 1
    assert len(positive_roots) == 1
    np.testing.assert_allclose(speed, positive_roots[0], rtol=0.0, atol=1e-14)
    departure = light.characterized_departure(
        P.rho_br,
        P.rho_B0,
        P.chi_c,
        P.J_phase,
        P.kappa_phase,
        P.longitudinal_display_fraction,
    )
    assert departure.diagnostics["longitudinal_dof"] == analysis.propagating_dof
