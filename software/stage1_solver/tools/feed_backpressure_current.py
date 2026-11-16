"""Steady, ungauged carrier on a flat segment with controlled reservoir phases.

Length 1/2, unit three-area, m = hbar = mu = 1, K = 4/5, V = 0.
Write psi(w,t) = (x + i y) exp(-i t). Solve for both fields with endpoint
amplitudes 1 and 9/10, then read the number current from the solution.
"""

import numpy as np
from scipy.integrate import solve_bvp


LENGTH = 0.5
OUTLET_AMPLITUDE = 0.9


def equations(w, state):
    x, y, dx, dy = state
    density = x**2 + y**2
    factor = 2 * (density**4 - 1)
    return np.vstack((dx, dy, factor * x, factor * y))


def solve(phase, grid_size, tolerance):
    outlet = OUTLET_AMPLITUDE * np.array([np.cos(phase), np.sin(phase)])
    feed = np.array([1.0, 0.0])

    def boundary(left, right):
        return np.concatenate((left[:2] - feed, right[:2] - outlet))

    grid = np.linspace(0, LENGTH, grid_size)
    slope = (outlet - feed) / LENGTH
    field = feed[:, None] + slope[:, None] * grid
    initial = np.vstack((field, np.broadcast_to(slope[:, None], field.shape)))
    result = solve_bvp(equations, boundary, grid, initial, tol=tolerance)
    if not result.success:
        raise RuntimeError(result.message)

    x, y, dx, dy = result.sol(np.linspace(0, LENGTH, 513))
    current = x * dy - y * dx
    error = np.max(np.abs(boundary(result.y[:, 0], result.y[:, -1])))
    spread = np.ptp(current)
    minimum_density = np.min(x**2 + y**2)
    print(f"phase={phase:.1f} tolerance={tolerance:.0e} "
          f"inlet={current[0]:.10f} outlet={current[-1]:.10f} "
          f"current_spread={spread:.2e} endpoint_error={error:.2e} "
          f"min_sampled_density={minimum_density:.8f}")
    if not (error <= 1e-12 and minimum_density > 0 and spread <= tolerance):
        raise ValueError("Endpoint, density or current-conservation check failed")
    return current[0]


if __name__ == "__main__":
    for phase in (0.0, 0.1, 0.2):
        coarse = solve(phase, 21, 1e-7)
        fine = solve(phase, 41, 1e-9)
        change = abs(coarse - fine)
        print(f"phase={phase:.1f} refinement_change={change:.2e}")
        if not change <= 1e-8:
            raise ValueError("Inlet current changed under refinement")
