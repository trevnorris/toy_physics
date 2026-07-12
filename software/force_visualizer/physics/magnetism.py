"""Velocity-dependent throat current--current interaction.

Authoritative sources are ``pathA_39_magnetic_force.md`` and results
(``kernel.compact`` / ``sign_diagnostic``), plus
``pathA_39_stage4_field_classification.md`` and
``pathA_39_scalar_admixture_screen.md`` for the unavoidable scalar-vector
field-content departure.  The implemented potential is the report's exact
static ``O(V1.V2)`` compact-source kernel, not a textbook Biot--Savart import.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import pi

import numpy as np
from numpy.typing import ArrayLike, NDArray

from .departures import CharacterizedDeparture
from .integrators import integrate_fixed

Vector = NDArray[np.float64]


def scalar_admixture_ratio(aT: float, aL: float, mu_R: float, B_eff: float) -> float:
    """Derived side-by-side ratio ``(aL^2/B_eff)/(aT^2/mu_R)``.

    Source: ``software/stage1_solver/reports/
    pathA_39_magnetic_force_results.yaml``, result key
    ``sign_diagnostic.scalar_admixture_ratio``.
    """

    if aT == 0.0 or mu_R <= 0.0 or B_eff <= 0.0:
        raise ValueError("aT must be nonzero and stable moduli positive")
    return aL**2 * mu_R / (aT**2 * B_eff)


def _geometry(
    position1: ArrayLike, position2: ArrayLike, velocity1: ArrayLike, velocity2: ArrayLike
) -> tuple[Vector, float, Vector, Vector, Vector, float, float, float, float]:
    x1 = np.asarray(position1, dtype=float)
    x2 = np.asarray(position2, dtype=float)
    v1 = np.asarray(velocity1, dtype=float)
    v2 = np.asarray(velocity2, dtype=float)
    if x1.shape != x2.shape or v1.shape != x1.shape or v2.shape != x1.shape:
        raise ValueError("positions and velocities must be same-shape vectors")
    separation = x2 - x1  # report convention R=X2-X1
    distance = float(np.linalg.norm(separation))
    if distance <= 0.0:
        raise ValueError("coincident point currents are singular")
    n = separation / distance
    dot = float(np.dot(v1, v2))
    radial1 = float(np.dot(v1, n))
    radial2 = float(np.dot(v2, n))
    radial_product = radial1 * radial2
    return separation, distance, n, v1, v2, dot, radial1, radial2, radial_product


def current_potential(
    position1: ArrayLike,
    position2: ArrayLike,
    velocity1: ArrayLike,
    velocity2: ArrayLike,
    orientation1: float,
    orientation2: float,
    N_u: float,
    aT: float,
    aL: float,
    mu_R: float,
    B_eff: float,
    include_scalar: bool = True,
) -> float:
    """Exact derived ``1/R`` transverse plus longitudinal exchange energy.

    ``U=-s1*s2*N_u^2/(8*pi*R) * [aT^2(D+A)/mu_R +
    aL^2(D-A)/B_eff]``.  Setting ``include_scalar=False`` is an explicit
    counterfactual display toggle; the characterized model keeps it on.
    Source: ``software/stage1_solver/reports/
    pathA_39_magnetic_force_results.yaml``, result key ``kernel.compact.U_12``.
    """

    _, distance, _, _, _, dot, _, _, radial_product = _geometry(
        position1, position2, velocity1, velocity2
    )
    if mu_R <= 0.0 or B_eff <= 0.0:
        raise ValueError("stable transverse and longitudinal moduli must be positive")
    transverse = aT**2 * (dot + radial_product) / mu_R
    longitudinal = aL**2 * (dot - radial_product) / B_eff if include_scalar else 0.0
    return -orientation1 * orientation2 * N_u**2 * (transverse + longitudinal) / (
        8.0 * pi * distance
    )


def force_on_second(
    position1: ArrayLike,
    position2: ArrayLike,
    velocity1: ArrayLike,
    velocity2: ArrayLike,
    orientation1: float,
    orientation2: float,
    N_u: float,
    aT: float,
    aL: float,
    mu_R: float,
    B_eff: float,
    include_scalar: bool = True,
) -> Vector:
    """Derived point force ``-grad_R U`` on body 2 (report convention).

    For side-by-side like parallel currents this vector is inward, for both
    stable transverse and stable longitudinal channels, matching
    ``NO_CANCELLATION_BOTH_CHANNELS_ATTRACTIVE``.
    Source: ``software/stage1_solver/reports/
    pathA_39_magnetic_force_results.yaml``, keys ``kernel.F_12`` and
    ``sign_diagnostic``.
    """

    _, distance, n, v1, v2, dot, radial1, radial2, radial_product = _geometry(
        position1, position2, velocity1, velocity2
    )
    transverse_weight = aT**2 / mu_R
    longitudinal_weight = aL**2 / B_eff if include_scalar else 0.0
    sum_weight = transverse_weight + longitudinal_weight
    difference_weight = transverse_weight - longitudinal_weight
    angular_kernel = sum_weight * dot + difference_weight * radial_product
    angular_gradient = difference_weight * (
        radial2 * v1 + radial1 * v2 - 2.0 * radial_product * n
    )
    prefactor = orientation1 * orientation2 * N_u**2 / (8.0 * pi * distance**2)
    return prefactor * (angular_gradient - angular_kernel * n)


def magnetic_field_on_brane(
    sample_points: ArrayLike,
    source_positions: ArrayLike,
    current_directions: ArrayLike,
    N_u: float,
    aT: float,
    aL: float,
    mu_R: float,
    B_eff: float,
    include_scalar: bool = True,
) -> Vector:
    """Build the §10 circulating display field from the verified kernel.

    For each radius, the magnitude is obtained by evaluating
    :func:`force_on_second` in its verified side-by-side geometry.  Only that
    magnitude is used: the vector is turned into the tangential direction
    required for the brane manifestation of a moving throat's 4D-body swirl,
    with circulation reversed when the current direction reverses.  This is
    a render-neutral visualization adapter, not an added Biot--Savart or force
    law.  The literal swirl remains in the throat 4D body.

    Sources: ``pathA_39_magnetic_force_results.yaml:kernel`` for the sampled
    magnitude; ``pathA_39_stage4_field_classification.md`` and
    ``software/force_visualizer/notes/build_spec.md`` §10 for field content
    and circulation semantics.
    """

    points = np.asarray(sample_points, dtype=float)
    sources = np.asarray(source_positions, dtype=float)
    directions = np.asarray(current_directions, dtype=float)
    if points.ndim < 1 or points.shape[-1] != 2:
        raise ValueError("sample_points must end in two spatial coordinates")
    if sources.ndim != 2 or sources.shape[1] != 2 or directions.shape != (len(sources),):
        raise ValueError("source_positions must be (N,2) and current_directions must be (N,)")
    if np.any(directions == 0.0):
        raise ValueError("current directions must be nonzero")

    flat_points = points.reshape(-1, 2)
    field = np.zeros_like(flat_points)
    for index, point in enumerate(flat_points):
        for position, direction in zip(sources, directions):
            displacement = point - position
            radius = float(np.linalg.norm(displacement))
            if radius <= 0.0:
                continue  # singular source pixel is deliberately left blank
            reference_force = force_on_second(
                [-radius / 2.0, 0.0],
                [radius / 2.0, 0.0],
                [0.0, float(direction)],
                [0.0, 1.0],
                1.0,
                1.0,
                N_u,
                aT,
                aL,
                mu_R,
                B_eff,
                include_scalar,
            )
            magnitude = float(np.linalg.norm(reference_force))
            tangent = np.array([-displacement[1], displacement[0]]) / radius
            field[index] += np.sign(direction) * magnitude * tangent
    return field.reshape(points.shape)


@dataclass(frozen=True)
class CurrentTrajectories:
    """Deterministic trajectories for compact moving throat sources.

    Architecture source: ``software/force_visualizer/notes/build_spec.md`` §2.
    """

    times: Vector
    positions: Vector
    velocities: Vector


def integrate_currents(
    positions: ArrayLike,
    velocities: ArrayLike,
    masses: ArrayLike,
    orientations: ArrayLike,
    N_u: float,
    aT: float,
    aL: float,
    mu_R: float,
    B_eff: float,
    dt: float,
    steps: int,
    include_scalar: bool = True,
) -> CurrentTrajectories:
    """Integrate the pairwise report kernel with fixed deterministic RK4.

    Law source: ``software/stage1_solver/reports/
    pathA_39_magnetic_force_results.yaml``, result key ``kernel``; integrator
    contract: ``software/force_visualizer/notes/build_spec.md`` §2.
    """

    x0 = np.asarray(positions, dtype=float)
    v0 = np.asarray(velocities, dtype=float)
    mass = np.asarray(masses, dtype=float)
    signs = np.asarray(orientations, dtype=float)
    if x0.ndim != 2 or x0.shape != v0.shape:
        raise ValueError("positions and velocities must be same-shape (N,D) arrays")
    bodies, dimension = x0.shape
    if mass.shape != (bodies,) or signs.shape != (bodies,) or np.any(mass <= 0.0):
        raise ValueError("masses/orientations must match bodies and masses be positive")
    initial = np.concatenate((x0.ravel(), v0.ravel()))

    def rhs(_time: float, state: Vector) -> Vector:
        x = state[: bodies * dimension].reshape(bodies, dimension)
        v = state[bodies * dimension :].reshape(bodies, dimension)
        forces = np.zeros_like(x)
        for first in range(bodies):
            for second in range(first + 1, bodies):
                on_second = force_on_second(
                    x[first],
                    x[second],
                    v[first],
                    v[second],
                    signs[first],
                    signs[second],
                    N_u,
                    aT,
                    aL,
                    mu_R,
                    B_eff,
                    include_scalar,
                )
                forces[second] += on_second
                forces[first] -= on_second
        return np.concatenate((v.ravel(), (forces / mass[:, None]).ravel()))

    times, states = integrate_fixed(rhs, initial, dt, steps)
    positions_out = states[:, : bodies * dimension].reshape(steps + 1, bodies, dimension)
    velocities_out = states[:, bodies * dimension :].reshape(steps + 1, bodies, dimension)
    return CurrentTrajectories(times, positions_out, velocities_out)


def characterized_departure(
    aT: float,
    aL: float,
    mu_R: float,
    B_eff: float,
    c_E: float,
    c_gamma: float,
) -> CharacterizedDeparture:
    """Return the mandatory scalar-current/field-content departure record.

    Sources: ``software/stage1_solver/reports/
    pathA_39_stage4_field_classification.md``, sections ``DOF And Constraints``
    and ``Scalar Stability``; ``pathA_39_scalar_admixture_screen.md``, section
    ``Decisive Residues``.
    """

    ratio = scalar_admixture_ratio(aT, aL, mu_R, B_eff)
    return CharacterizedDeparture(
        code="FIELD_SCALAR_VECTOR_DEPARTURE",
        description="stable transverse vector plus uncancelable attractive scalar-current channel and h-branon",
        derived_form="U_L=-s1*s2*N_u^2*aL^2*(D-A)/(8*pi*B_eff*R); same attractive side-by-side sign",
        calibrated_magnitude="aT and aL are sim-deferred/calibrated; C_hu and scalar residues remain unreduced",
        diagnostics={
            "scalar/transverse_ratio": ratio,
            "physical_dof": 4.0,
            "preferred_frame_residual": c_E**2 - c_gamma**2,
            "c_E_equals_c_gamma": str(bool(np.isclose(c_E, c_gamma))),
        },
    )
