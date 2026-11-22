"""Stationary two-bath test of a held order front with passive reciprocal conversion.

The pressure difference is specified; the current and the force exerted by the
translation constraint are determined by the boundary-value problem. Units set
n0, constituent mass, hbar, a_B and delta_B to one (kappa_B = 2).
"""

import argparse
import json

import numpy as np
from scipy.integrate import simpson, solve_bvp
from scipy.special import expit


def measure(load, mobility, width, extent, spacing, tolerance, *, return_fields=False,
            endpoint_heads=None):
    # s increases toward the bulk. Projecting out translation holds the front;
    # the local order profile and the carrier density on both sides remain free.
    points = np.linspace(-extent, extent, int(round(2 * extent / spacing)) + 1)
    n_bulk = (1.0 + 5.0 * load) ** 0.2  # p(n) = n^5 / 5

    def fixture(s):
        q = expit(-s / width)
        return -q * (1.0 - q) / width

    def bath_head(n, current):
        return (n**4 - 1.0) / 4.0 + 0.5 * (current / n) ** 2

    def equations(s, y, params):
        f, df, q, dq, head, _ = y
        current, fixture_strength = params
        n = f**2
        slope_well = 2.0 * q * (1.0 - q) * (1.0 - 2.0 * q)
        affinity = -current * dq / (mobility * n)
        return np.array([
            df, 2.0 * f * (bath_head(n, current) - head),
            dq, (slope_well - fixture_strength * fixture(s) - affinity) / 2.0,
            affinity * dq / n, fixture(s) * (q - expit(-s)),
        ])

    def boundary(left, right, params):
        current = params[0]
        heads = endpoint_heads
        if heads is None:
            heads = (bath_head(1.0, current), bath_head(n_bulk, current))
        return np.array([
            left[0] - 1.0, right[0] - np.sqrt(n_bulk),
            left[2] - expit(extent), right[2] - expit(-extent),
            left[4] - heads[0], right[4] - heads[1],
            left[5], right[5],
        ])

    q = expit(-points)
    f = 1.0 + (np.sqrt(n_bulk) - 1.0) * (1.0 - q)
    guess = np.array([f, np.zeros_like(f), q, -q * expit(points),
                      (1.0 - q) * bath_head(n_bulk, 0), np.zeros_like(q)])
    result = solve_bvp(equations, boundary, points, guess, p=[0.0, 0.0],
                       tol=tolerance, max_nodes=30000)
    if not result.success:
        raise RuntimeError(result.message)

    s = np.linspace(-extent, extent, 4001)
    f, df, q, dq, head, _ = result.sol(s)
    ddf, ddq = result.sol(s, 1)[[1, 3]]
    current, strength = result.p
    n = f**2
    affinity = (2.0 * q * (1.0 - q) * (1.0 - 2.0 * q) - 2.0 * ddq
                - strength * fixture(s))
    force_external = simpson(-strength * fixture(s) * dq, x=s)
    # Normal stress from coherent flux plus square-gradient order stress.
    total_flux = (current**2 / n + n**5 / 5.0 + (df**2 - f * ddf) / 2.0
                  + dq**2 - (q * (1.0 - q))**2)
    from_baths = total_flux[-1] - total_flux[0]
    work = current * (head[0] - head[-1])
    dissipation = simpson(mobility * affinity**2, x=s)
    row = dict(bulk_pressure_excess=load, mobility=mobility, fixture_width=width,
                extent=extent, mesh_step=spacing, tolerance=tolerance,
                current_to_bulk=float(current), fixture_on_medium=float(force_external),
                normal_flux_difference=float(from_baths), bath_power=float(work),
                conversion_power=float(dissipation), nodes=len(result.x),
                rms_residual_max=float(max(result.rms_residuals)),
                boundary_residual_max=float(max(abs(boundary(result.y[:, 0], result.y[:, -1], result.p)))),
                order_crossing_position=float(np.interp(0.5, q[::-1], s[::-1])),
                center_density=float(n[len(s) // 2]))
    if return_fields:
        return row, result, fixture
    return row


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--extent", type=float, default=20.0)
    parser.add_argument("--spacing", type=float, default=0.125)
    parser.add_argument("--tolerance", type=float, default=1.0e-9)
    args = parser.parse_args()
    for mobility, width in [(0.5, 1.0), (1.0, 1.0), (0.5, 1.5)]:
        for load in (0.0, 0.0001, -0.0001, 0.00005, -0.00005):
            row = measure(load, mobility, width, args.extent, args.spacing, args.tolerance)
            print(json.dumps(row, sort_keys=True, allow_nan=False), flush=True)
            if abs(row["normal_flux_difference"] - row["fixture_on_medium"]) > 2e-9:
                raise ArithmeticError("Fixture and remote momentum flux disagree")
            if abs(row["bath_power"] - row["conversion_power"]) > 5e-13:
                raise ArithmeticError("Reservoir and conversion powers disagree")


if __name__ == "__main__":
    main()
