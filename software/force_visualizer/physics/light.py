"""Transverse brane waves, stray longitudinal mode, and refractive rays.

Sources
-------
* ``pathA_36_c5_phase_potential.md`` / results and read-only
  ``ledger_stage003_transverse_photons_stray_longitudinal.md``: exactly two
  transverse modes with ``omega^2=(mu_R/rho_br)k^2`` and one characterized
  provenance-fixed longitudinal second-class mode.
* read-only ``ledger_stage005_sound_speed_light_ratio.md`` and pathA_40:
  ``lambda_gamma=c_gamma/c_s`` is a free calibration, not a derived equality.
* ``research/1pn_optics/paper/1pn_optics.tex``, Sec. "Gravitational lensing":
  ``N(r)=1+2GM/(c_gamma^2 r)`` and Hamiltonian geometric-optics rays.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import ceil, pi, sqrt
from typing import Tuple

import numpy as np
from numpy.typing import ArrayLike, NDArray

from .departures import CharacterizedDeparture
from .integrators import integrate_fixed

Vector = NDArray[np.float64]


@dataclass(frozen=True)
class TransverseWaveEvolution:
    """Periodic finite-difference evolution of one transverse wave component."""

    coordinates: Vector
    times: Vector
    displacements: Vector
    wavenumber: float


@dataclass(frozen=True)
class LongitudinalModeAnalysis:
    """Roots and physical mode count obtained from the longitudinal pole."""

    wavenumber: float
    angular_frequency_roots: NDArray[np.complex128]
    propagating_dof: int


def photon_speed(mu_R: float, rho_br: float) -> float:
    """Derived transverse speed ``sqrt(mu_R/rho_br)``.

    Source: ``software/stage1_solver/reports/pathA_36_c5_phase_potential.md``,
    section ``Transverse Sector``.
    """

    if mu_R <= 0.0 or rho_br <= 0.0:
        raise ValueError("mu_R and rho_br must be positive")
    return sqrt(mu_R / rho_br)


def angular_frequency(wavenumber: ArrayLike, c_gamma: float) -> Vector | float:
    """Massless derived dispersion ``omega=c_gamma*|k|``.

    Source: ``research/pde_ledger_v2/notes/stages/
    ledger_stage003_transverse_photons_stray_longitudinal.md`` §2.
    """

    if c_gamma <= 0.0:
        raise ValueError("c_gamma must be positive")
    values = np.asarray(wavenumber, dtype=float)
    result = c_gamma * np.abs(values)
    return float(result) if result.ndim == 0 else result


def evolve_transverse_wave(
    wavenumber: float,
    c_gamma: float,
    *,
    spatial_points: int = 300,
    periods: float = 1.0,
    courant: float = 0.4,
    phase: float = 0.0,
) -> TransverseWaveEvolution:
    """Evolve ``u_tt=c_gamma^2*u_xx`` with the scene's leapfrog stepper.

    The periodic domain is ``[0, 2*pi)`` and therefore accepts integer Fourier
    wavenumbers.  The first step is a centered Taylor step using the requested
    right-moving initial velocity; all subsequent samples use the explicit
    second-order finite-difference update.
    """

    mode = int(round(abs(wavenumber)))
    if mode < 1 or not np.isclose(abs(wavenumber), mode, rtol=0.0, atol=1.0e-12):
        raise ValueError("wavenumber must be a nonzero integer on the 2*pi domain")
    if c_gamma <= 0.0 or spatial_points < 16 or periods <= 0.0:
        raise ValueError("require c_gamma>0, spatial_points>=16 and periods>0")
    if not 0.0 < courant <= 1.0:
        raise ValueError("the explicit wave stepper requires 0<courant<=1")

    coordinates = np.linspace(0.0, 2.0 * pi, spatial_points, endpoint=False)
    dx = 2.0 * pi / spatial_points
    duration = periods * 2.0 * pi / (c_gamma * mode)
    steps = max(2, int(ceil(duration / (courant * dx / c_gamma))))
    dt = duration / steps
    courant_actual = c_gamma * dt / dx

    initial = np.sin(wavenumber * coordinates + phase)
    initial_velocity = -c_gamma * wavenumber * np.cos(
        wavenumber * coordinates + phase
    )
    laplacian_stencil = np.roll(initial, -1) - 2.0 * initial + np.roll(initial, 1)
    first = initial + dt * initial_velocity + 0.5 * courant_actual**2 * laplacian_stencil

    displacements = np.empty((steps + 1, spatial_points), dtype=float)
    displacements[0] = initial
    displacements[1] = first
    for index in range(1, steps):
        current = displacements[index]
        displacements[index + 1] = (
            2.0 * current
            - displacements[index - 1]
            + courant_actual**2
            * (np.roll(current, -1) - 2.0 * current + np.roll(current, 1))
        )
    return TransverseWaveEvolution(
        coordinates=coordinates,
        times=np.linspace(0.0, duration, steps + 1),
        displacements=displacements,
        wavenumber=float(wavenumber),
    )


def measured_transverse_angular_frequency(evolution: TransverseWaveEvolution) -> float:
    """Measure temporal phase rotation from a finite-difference evolution."""

    sine = np.sin(evolution.wavenumber * evolution.coordinates)
    cosine = np.cos(evolution.wavenumber * evolution.coordinates)
    sine_projection = evolution.displacements @ sine
    cosine_projection = evolution.displacements @ cosine
    phase = np.unwrap(np.arctan2(-cosine_projection, sine_projection))
    slope, _ = np.polyfit(evolution.times, phase, 1)
    return abs(float(slope))


def transverse_polarizations(wavevector: ArrayLike) -> Vector:
    """Return exactly two orthonormal polarization vectors perpendicular to k.

    This realizes the two physical transverse DOF derived in pathA_36; it does
    not add a longitudinal Maxwell polarization.
    Source: ``software/stage1_solver/reports/
    pathA_36_c5_phase_potential_results.yaml``, result key
    ``transverse.physical_dof=2``.
    """

    k = np.asarray(wavevector, dtype=float)
    if k.shape != (3,):
        raise ValueError("the brane wavevector must be a three-vector")
    norm = float(np.linalg.norm(k))
    if norm <= 0.0:
        raise ValueError("a zero wavevector has no unique transverse basis")
    direction = k / norm
    reference = np.eye(3)[int(np.argmin(np.abs(direction)))]
    first = np.cross(direction, reference)
    first /= np.linalg.norm(first)
    second = np.cross(direction, first)
    second /= np.linalg.norm(second)
    return np.vstack((first, second))


def transverse_plane_wave(
    positions: ArrayLike,
    time: float,
    wavevector: ArrayLike,
    polarization: int,
    c_gamma: float,
    amplitude: float = 1.0,
    phase: float = 0.0,
) -> Vector:
    """Evaluate one of the two physical transverse displacement modes.

    Source: ``software/stage1_solver/reports/pathA_36_c5_phase_potential.md``,
    section ``Transverse Sector`` and its ``L_T`` wave equation.
    """

    x = np.asarray(positions, dtype=float)
    k = np.asarray(wavevector, dtype=float)
    if x.shape[-1] != 3:
        raise ValueError("positions must end in three brane coordinates")
    basis = transverse_polarizations(k)
    if polarization not in (0, 1):
        raise ValueError("polarization must be 0 or 1")
    omega = c_gamma * np.linalg.norm(k)
    scalar = amplitude * np.cos(x @ k - omega * time + phase)
    return scalar[..., None] * basis[polarization]


def stray_longitudinal_speed(
    rho_br: float,
    rho_B0: float,
    chi_c: float,
    J_phase: float,
    kappa_phase: float,
) -> float:
    """Derived provenance-fixed longitudinal pole speed from pathA_36.

    ``omega^2/k^2 = kappa_phase*rho_B0^2 /
    (chi_c*(J^2*rho_B0^2+kappa_phase*rho_br))``.
    Source: ``software/stage1_solver/reports/
    pathA_36_c5_phase_potential_results.yaml``, result key
    ``per_branch_sub_results.branch_b_slaved_finite_compressibility_conventional_K.dispersion``.
    """

    analysis = analyze_longitudinal_mode(
        1.0, rho_br, rho_B0, chi_c, J_phase, kappa_phase
    )
    if analysis.propagating_dof != 1:
        raise RuntimeError("the longitudinal dispersion has no real propagating mode")
    positive_root = max(
        root.real
        for root in analysis.angular_frequency_roots
        if abs(root.imag) <= 1.0e-10 and root.real > 0.0
    )
    return float(positive_root)


def analyze_longitudinal_mode(
    wavenumber: float,
    rho_br: float,
    rho_B0: float,
    chi_c: float,
    J_phase: float,
    kappa_phase: float,
) -> LongitudinalModeAnalysis:
    """Solve the longitudinal dispersion polynomial and count real modes.

    The pole equation is
    ``chi_c*(J^2*rho_B0^2+kappa*rho_br)*omega^2
    - kappa*rho_B0^2*k^2 = 0``.  A positive-frequency real root and its
    negative-frequency partner constitute one propagating physical DOF.
    """

    values = (rho_br, rho_B0, chi_c, kappa_phase)
    if any(value <= 0.0 for value in values) or wavenumber == 0.0:
        raise ValueError(
            "rho_br, rho_B0, chi_c and kappa_phase must be positive and k nonzero"
        )
    kinetic = chi_c * (J_phase**2 * rho_B0**2 + kappa_phase * rho_br)
    gradient = kappa_phase * rho_B0**2
    roots = np.asarray(
        np.roots([kinetic, 0.0, -gradient * wavenumber**2]), dtype=np.complex128
    )
    tolerance = 1.0e-10
    positive_real_roots = sum(
        abs(root.imag) <= tolerance and root.real > tolerance for root in roots
    )
    return LongitudinalModeAnalysis(
        wavenumber=float(wavenumber),
        angular_frequency_roots=roots,
        propagating_dof=int(positive_real_roots),
    )


def refractive_index(position: ArrayLike, mass: float, G: float, c_gamma: float) -> float:
    """Weak-field ``n(r)=1+2GM/(c_gamma^2 r)`` from 1PN optics.

    Source: ``research/1pn_optics/paper/1pn_optics.tex``, equation
    ``lensing-N`` in section ``Refractive index profile as input``.
    """

    point = np.asarray(position, dtype=float)
    radius = float(np.linalg.norm(point))
    if radius <= 0.0 or mass <= 0.0 or G <= 0.0 or c_gamma <= 0.0:
        raise ValueError("ray radius, mass, G and c_gamma must be positive")
    return 1.0 + 2.0 * G * mass / (c_gamma**2 * radius)


def refractive_index_gradient(position: ArrayLike, mass: float, G: float, c_gamma: float) -> Vector:
    """Analytic gradient of the cited weak-field refractive profile.

    Source: ``research/1pn_optics/paper/1pn_optics.tex``, section
    ``Analytic evaluation of the bending angle``, derivative of ``ln N``.
    """

    point = np.asarray(position, dtype=float)
    radius = float(np.linalg.norm(point))
    if radius <= 0.0:
        raise ValueError("the point-mass refractive gradient is singular at the origin")
    coefficient = 2.0 * G * mass / c_gamma**2
    return -coefficient * point / radius**3


def local_phase_speed_and_gradient(
    position: ArrayLike, mass: float, G: float, c_gamma: float
) -> Tuple[float, Vector]:
    """Return ``c_gamma/n`` and its gradient for Hamiltonian ray tracing.

    Source: ``research/1pn_optics/paper/1pn_optics.tex``, appendix section
    ``Ray-tracing algorithms``, Hamiltonian ``H=c/N*|k|``.
    """

    index = refractive_index(position, mass, G, c_gamma)
    grad_index = refractive_index_gradient(position, mass, G, c_gamma)
    speed = c_gamma / index
    grad_speed = -(c_gamma / index**2) * grad_index
    return speed, grad_speed


@dataclass(frozen=True)
class RayTrace:
    """A deterministic 2D geometric-optics ray.

    Source: ``research/1pn_optics/paper/1pn_optics.tex``, appendix
    ``Ray-tracing algorithms``.
    """

    positions: Vector
    wavevectors: Vector

    @property
    def signed_deflection(self) -> float:
        """Outgoing angle relative to the initial +z direction.

        Source: ``research/1pn_optics/paper/1pn_optics.tex``, appendix
        ray-tracing outgoing-direction deflection definition.
        """

        final = self.wavevectors[-1]
        initial = self.wavevectors[0]
        return float(
            np.arctan2(final[0], final[1]) - np.arctan2(initial[0], initial[1])
        )


def trace_ray(
    impact_parameter: float,
    z_extent: float,
    mass: float,
    G: float,
    c_gamma: float,
    steps: int = 4000,
) -> RayTrace:
    """Trace the cited Hamiltonian ray equations, reparameterized by z.

    The Hamiltonian is ``H=(c_gamma/n(x))*|k|``.  Reparameterization by the
    monotone coordinate ``z`` leaves the geometric ray unchanged and gives a
    fixed integration interval suitable for deterministic tests.
    Source: ``research/1pn_optics/paper/1pn_optics.tex``, sections
    ``Numerical ray-tracing`` and appendix ``Ray-tracing algorithms``.
    """

    if impact_parameter == 0.0 or z_extent <= 0.0 or steps < 2:
        raise ValueError("nonzero impact, positive extent, and steps>=2 are required")
    initial = np.array([impact_parameter, 0.0, 1.0], dtype=float)  # x, k_x, k_z
    dz = 2.0 * z_extent / steps

    def rhs(offset: float, state: Vector) -> Vector:
        z = -z_extent + offset
        x, kx, kz = state
        if kz <= 0.0:
            raise RuntimeError("ray ceased to be monotone in z")
        speed, grad_speed = local_phase_speed_and_gradient(
            np.array([x, z]), mass, G, c_gamma
        )
        wave_norm_squared = kx**2 + kz**2
        factor = -wave_norm_squared / (speed * kz)
        return np.array([kx / kz, factor * grad_speed[0], factor * grad_speed[1]])

    offsets, states = integrate_fixed(rhs, initial, dz, steps)
    z_values = -z_extent + offsets
    positions = np.column_stack((states[:, 0], z_values))
    wavevectors = states[:, 1:3]
    return RayTrace(positions, wavevectors)


def analytic_signed_deflection(
    impact_parameter: float, mass: float, G: float, c_gamma: float
) -> float:
    """Leading signed bend toward the mass, magnitude ``4GM/(|b|c^2)``.

    Source: ``research/1pn_optics/paper/1pn_optics.tex``, equation
    ``lensing-dtheta-result``.
    """

    if impact_parameter == 0.0:
        raise ValueError("impact parameter must be nonzero")
    return -np.sign(impact_parameter) * 4.0 * G * mass / (
        abs(impact_parameter) * c_gamma**2
    )


def characterized_departure(
    rho_br: float,
    rho_B0: float,
    chi_c: float,
    J_phase: float,
    kappa_phase: float,
    display_fraction: float,
) -> CharacterizedDeparture:
    """Return the mandatory stray-longitudinal honest-scope record.

    Source: ``research/pde_ledger_v2/notes/stages/
    ledger_stage003_transverse_photons_stray_longitudinal.md`` §§4-5.
    """

    analysis = analyze_longitudinal_mode(
        1.0, rho_br, rho_B0, chi_c, J_phase, kappa_phase
    )
    speed = stray_longitudinal_speed(rho_br, rho_B0, chi_c, J_phase, kappa_phase)
    return CharacterizedDeparture(
        code="FAIL_CAUCHY_STRAY_LONGITUDINAL",
        description="one physical second-class longitudinal propagating DOF; Maxwell has none",
        derived_form="omega^2=k^2*kappa*rho_B0^2/[chi_c*(J^2*rho_B0^2+kappa*rho_br)]",
        calibrated_magnitude="material/J/kappa magnitudes and visualization amplitude are not frozen",
        diagnostics={
            "longitudinal_speed": speed,
            "display_fraction": display_fraction,
            "longitudinal_dof": float(analysis.propagating_dof),
        },
    )
