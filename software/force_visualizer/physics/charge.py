"""Brane-localized electric point law and characterized Yukawa partner.

Authoritative source: ``software/stage1_solver/reports/
pathA_38_throat_body_electric_localization.md`` and
``pathA_38_results.yaml`` (``goldstone``, ``source_projections``,
``green_function``, ``interaction_sign`` and ``fail_witnesses.FAIL_YUKAWA``).
The source derives an odd ``+/-`` throat orientation, a positive zero-mode
kernel (like repel / unlike attract), and massive Yukawa partner kernels.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import pi, sqrt, tanh
from typing import Tuple

import numpy as np
from numpy.typing import ArrayLike, NDArray

from .departures import CharacterizedDeparture
from .integrators import integrate_fixed

Vector = NDArray[np.float64]


def zero_mode_norm(ell: float) -> float:
    """Derived Goldstone norm ``N0=8/(3 ell)``.

    Source: ``software/stage1_solver/reports/pathA_38_results.yaml``,
    result key ``goldstone.N0_norm``.
    """

    if ell <= 0.0:
        raise ValueError("ell must be positive")
    return 8.0 / (3.0 * ell)


def projected_charge(Q_E: float, ell: float, mouth_half_width: float, orientation: float) -> float:
    """Derived compact throat projection ``2 Q_E tanh(b/ell)/b`` times sign.

    Source: ``software/stage1_solver/reports/pathA_38_results.yaml``,
    result keys ``source_projections.q_h_plus/q_h_minus``.
    """

    if ell <= 0.0 or mouth_half_width <= 0.0:
        raise ValueError("ell and mouth_half_width must be positive")
    return orientation * 2.0 * Q_E * tanh(mouth_half_width / ell) / mouth_half_width


def electric_coupling(Q_E: float, ell: float, mouth_half_width: float) -> float:
    """Positive projected ``1/R`` coefficient from pathA_38.

    This is ``3 Q_E^2 ell tanh(b/ell)^2/(8 pi b^2)``.  ``Q_E`` is a
    calibrated anchor; every other factor is part of the derived projection.
    Source: ``software/stage1_solver/reports/pathA_38_results.yaml``, result
    key ``green_function.static_projected_U_kernel``.
    """

    if Q_E <= 0.0 or ell <= 0.0 or mouth_half_width <= 0.0:
        raise ValueError("Q_E, ell and mouth_half_width must be positive")
    return 3.0 * Q_E**2 * ell * tanh(mouth_half_width / ell) ** 2 / (
        8.0 * pi * mouth_half_width**2
    )


def coulomb_potential(
    distance: float,
    orientation1: float,
    orientation2: float,
    Q_E: float,
    ell: float,
    mouth_half_width: float,
) -> float:
    """Derived brane-localized Coulomb energy ``+k_E s1 s2/R``.

    Source: ``software/stage1_solver/reports/pathA_38_results.yaml``, result
    key ``interaction_sign``.
    """

    if distance <= 0.0:
        raise ValueError("distance must be positive")
    return electric_coupling(Q_E, ell, mouth_half_width) * orientation1 * orientation2 / distance


def yukawa_potential(
    distance: float,
    orientation1: float,
    orientation2: float,
    Q_E: float,
    ell: float,
    mouth_half_width: float,
    relative_amplitude: float,
) -> float:
    """Shape-partner correction ``alpha*k_E*s1*s2*exp(-sqrt(3)R/ell)/R``.

    The gap/form are derived by pathA_38.  The main compact-source residue is
    not universally fixed there, so ``relative_amplitude`` is calibrated and
    must not be described as derived.
    Source: ``software/stage1_solver/reports/pathA_38_results.yaml``, result
    keys ``transverse_mode_spectrum.wall_shape_partner`` and
    ``green_function.static_shape_mode_kernel_eta_eta``.
    """

    if relative_amplitude < 0.0:
        raise ValueError("relative_amplitude must be non-negative")
    mass_gap = sqrt(3.0) / ell
    return relative_amplitude * coulomb_potential(
        distance, orientation1, orientation2, Q_E, ell, mouth_half_width
    ) * np.exp(-mass_gap * distance)


def force_on_first(
    position1: ArrayLike,
    position2: ArrayLike,
    orientation1: float,
    orientation2: float,
    Q_E: float,
    ell: float,
    mouth_half_width: float,
    include_yukawa: bool = False,
    yukawa_fraction: float = 0.0,
) -> Vector:
    """Force on body 1; like orientations repel and unlike attract.

    Source: ``software/stage1_solver/reports/pathA_38_results.yaml``, result
    keys ``interaction_sign`` and ``static_dynamic_consistency.p_static``.
    """

    displacement = np.asarray(position1, dtype=float) - np.asarray(position2, dtype=float)
    distance = float(np.linalg.norm(displacement))
    if distance <= 0.0:
        raise ValueError("coincident point charges are singular")
    coupling = electric_coupling(Q_E, ell, mouth_half_width)
    multiplier = 1.0
    if include_yukawa:
        if yukawa_fraction < 0.0:
            raise ValueError("yukawa_fraction must be non-negative")
        gap = sqrt(3.0) / ell
        multiplier += yukawa_fraction * np.exp(-gap * distance) * (1.0 + gap * distance)
    return coupling * orientation1 * orientation2 * multiplier * displacement / distance**3


def electric_field(
    sample_points: ArrayLike,
    source_positions: ArrayLike,
    source_orientations: ArrayLike,
    Q_E: float,
    ell: float,
    mouth_half_width: float,
) -> Vector:
    """Evaluate the Coulomb field as force on a positive unit orientation.

    This is a sampling adapter over :func:`force_on_first`, not a new law.
    It deliberately excludes the optional Yukawa partner so the overlay is
    the long-range electric field requested by ``build_spec.md`` §10.  It is
    a throat-body interaction field, never a medium-flow velocity.

    Source: ``pathA_38_results.yaml``, ``interaction_sign`` and
    ``green_function.static_projected_U_kernel``.
    """

    points = np.asarray(sample_points, dtype=float)
    sources = np.asarray(source_positions, dtype=float)
    signs = np.asarray(source_orientations, dtype=float)
    if points.ndim < 1 or points.shape[-1] != 2:
        raise ValueError("sample_points must end in two spatial coordinates")
    if sources.ndim != 2 or sources.shape[1] != 2 or signs.shape != (len(sources),):
        raise ValueError("source_positions must be (N,2) and orientations must be (N,)")

    flat_points = points.reshape(-1, 2)
    field = np.zeros_like(flat_points)
    for index, point in enumerate(flat_points):
        for position, sign in zip(sources, signs):
            if np.array_equal(point, position):
                continue  # singular source pixel is deliberately left blank
            field[index] += force_on_first(
                point,
                position,
                1.0,
                sign,
                Q_E,
                ell,
                mouth_half_width,
            )
    return field.reshape(points.shape)


@dataclass(frozen=True)
class ChargeTrajectories:
    """Deterministic many-charge point trajectories.

    Architecture source: ``software/force_visualizer/notes/build_spec.md`` §2.
    """

    times: Vector
    positions: Vector
    velocities: Vector


def integrate_charges(
    positions: ArrayLike,
    velocities: ArrayLike,
    masses: ArrayLike,
    orientations: ArrayLike,
    Q_E: float,
    ell: float,
    mouth_half_width: float,
    dt: float,
    steps: int,
    include_yukawa: bool = False,
    yukawa_fraction: float = 0.0,
) -> ChargeTrajectories:
    """Integrate pairwise electric dynamics with deterministic fixed-step RK4.

    Law source: ``software/stage1_solver/reports/pathA_38_results.yaml``, keys
    ``green_function`` and ``interaction_sign``; integrator contract: spec §2.
    """

    x0 = np.asarray(positions, dtype=float)
    v0 = np.asarray(velocities, dtype=float)
    mass = np.asarray(masses, dtype=float)
    signs = np.asarray(orientations, dtype=float)
    if x0.ndim != 2 or v0.shape != x0.shape:
        raise ValueError("positions and velocities must be same-shape (N,D) arrays")
    bodies, dimension = x0.shape
    if mass.shape != (bodies,) or signs.shape != (bodies,) or np.any(mass <= 0.0):
        raise ValueError("masses/orientations must match bodies and masses must be positive")
    initial = np.concatenate((x0.ravel(), v0.ravel()))

    def rhs(_time: float, state: Vector) -> Vector:
        x = state[: bodies * dimension].reshape(bodies, dimension)
        v = state[bodies * dimension :].reshape(bodies, dimension)
        forces = np.zeros_like(x)
        for first in range(bodies):
            for second in range(first + 1, bodies):
                pair_force = force_on_first(
                    x[first],
                    x[second],
                    signs[first],
                    signs[second],
                    Q_E,
                    ell,
                    mouth_half_width,
                    include_yukawa,
                    yukawa_fraction,
                )
                forces[first] += pair_force
                forces[second] -= pair_force
        return np.concatenate((v.ravel(), (forces / mass[:, None]).ravel()))

    times, states = integrate_fixed(rhs, initial, dt, steps)
    positions_out = states[:, : bodies * dimension].reshape(steps + 1, bodies, dimension)
    velocities_out = states[:, bodies * dimension :].reshape(steps + 1, bodies, dimension)
    return ChargeTrajectories(times, positions_out, velocities_out)


def characterized_departure(ell: float, relative_amplitude: float) -> CharacterizedDeparture:
    """Return the mandatory charge-sector short-range departure record.

    Source: ``software/stage1_solver/reports/pathA_38_results.yaml``, keys
    ``green_function.massive_terms`` and ``fail_witnesses.FAIL_YUKAWA``.
    """

    gap = sqrt(3.0) / ell
    return CharacterizedDeparture(
        code="YUKAWA_PARTNER_CORRECTION",
        description="exponentially suppressed gapped wall-partner correction above Coulomb",
        derived_form="exp(-sqrt(3)*R/ell)/R with force factor exp(-mR)*(1+mR)/R^2",
        calibrated_magnitude="the compact main-source partner residue is not fixed; relative amplitude is calibrated",
        diagnostics={"mass_gap": gap, "relative_amplitude": relative_amplitude},
    )
