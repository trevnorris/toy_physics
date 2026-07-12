"""Small deterministic numerical integrators shared by the physics sectors.

Architecture source: ``software/force_visualizer/notes/build_spec.md`` §§2,6.
"""

from __future__ import annotations

from typing import Callable, Tuple

import numpy as np
from numpy.typing import ArrayLike, NDArray

State = NDArray[np.float64]
RHS = Callable[[float, State], State]


def rk4_step(rhs: RHS, t: float, state: ArrayLike, dt: float) -> State:
    """Advance one fixed classical fourth-order Runge--Kutta step.

    Architecture source: ``software/force_visualizer/notes/build_spec.md`` §2
    (builder-owned deterministic integrator; no rendering dependency).
    """

    y = np.asarray(state, dtype=float)
    k1 = np.asarray(rhs(t, y), dtype=float)
    k2 = np.asarray(rhs(t + 0.5 * dt, y + 0.5 * dt * k1), dtype=float)
    k3 = np.asarray(rhs(t + 0.5 * dt, y + 0.5 * dt * k2), dtype=float)
    k4 = np.asarray(rhs(t + dt, y + dt * k3), dtype=float)
    return y + (dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)


def integrate_fixed(rhs: RHS, initial_state: ArrayLike, dt: float, steps: int) -> Tuple[State, State]:
    """Return fixed-step times and states with no hidden randomness.

    Acceptance source: ``software/force_visualizer/notes/build_spec.md`` §6,
    determinism golden.
    """

    if dt <= 0.0 or steps < 1:
        raise ValueError("dt must be positive and steps must be at least one")
    initial = np.asarray(initial_state, dtype=float)
    states = np.empty((steps + 1,) + initial.shape, dtype=float)
    times = np.arange(steps + 1, dtype=float) * dt
    states[0] = initial
    for index in range(steps):
        states[index + 1] = rk4_step(rhs, times[index], states[index], dt)
    return times, states
