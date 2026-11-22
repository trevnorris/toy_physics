"""Flat continuous-carrier path given the saved mouth's two reservoir phases.

The ordered and disordered sides use the same longitudinal Hamiltonian. This
is a current through a fixed plane in the additive neutral limit. The distance
between the baths is an explicit control rather than a lifetime or wall law.
"""

import json
from pathlib import Path

import numpy as np
from scipy.integrate import simpson, solve_bvp


BASE = Path(__file__).resolve().parents[3]
FLARE = BASE / "research/pde_ledger_v3/_measurements/mouth_free_radial_direction.json"


def saved_bath():
    record = json.loads(FLARE.read_text())
    rows = (json.loads(line) for run in record["runs"].values()
            for line in run["literal_stdout"].splitlines())
    target, = (row for row in rows if row.get("test") == "force"
               and row["mesh"] == 34 and row["q"] == 0.0)
    return target["bias"]


def solve(length, phase_drop):
    # psi=f*exp(i theta), j=f^2 theta'. Positive j points toward bulk.
    x = np.linspace(0.0, length, 301)
    weak_current = phase_drop / length

    def field(_, state, parameters):
        f, df, _ = state
        current = parameters[0]
        density = f**2
        phase_slope = current / density
        return np.array([df, f * ((density**4 - 1.0) / 2.0 + phase_slope**2),
                         phase_slope])

    def bath(left, right, _):
        return np.array([left[0] - 1.0, right[0] - 1.0, left[2], right[2] - phase_drop])

    guess = np.array([np.ones_like(x), np.zeros_like(x), weak_current * x])
    sol = solve_bvp(field, bath, x, guess, p=[weak_current], tol=2e-10)
    if not sol.success:
        raise ArithmeticError(sol.message)
    s = np.linspace(0.0, length, 4097)
    f, df, _ = sol.sol(s)
    n = f**2
    _, ddf, dtheta = sol.sol(s, 1)
    current = float(sol.p[0])
    pressure = n**5 / 5.0 + n * dtheta**2 + (df**2 - f * ddf) / 2.0
    row = dict(length=length, bulk_minus_brane_phase=phase_drop,
               lab_current=current, linear_current=weak_current,
               phase_integral=current * float(simpson(1.0 / n, x=s)),
               local_current_max_error=float(np.max(abs(n * dtheta - current))),
               normal_stress_range=float(np.ptp(pressure)), n_mid=float(n[len(s) // 2]),
               rms_residual_max=float(np.max(sol.rms_residuals)))
    print(json.dumps(row, allow_nan=False), flush=True)
    assert abs(row["phase_integral"] - phase_drop) < 1e-9
    assert row["local_current_max_error"] < 1e-9 and row["normal_stress_range"] < 5e-9
    return current


def main():
    phase_drop = saved_bath()
    for length in (36.0, 72.0):
        j_zero, j_plus, j_minus, j_half = [solve(length, scale * phase_drop)
                                           for scale in (0, 1, -1, 0.5)]
        assert abs(j_zero) < 1e-12 and abs(j_plus + j_minus) < 1e-10
        assert abs(2.0 * j_half / j_plus - 1.0) < 2e-4


if __name__ == "__main__":
    main()
