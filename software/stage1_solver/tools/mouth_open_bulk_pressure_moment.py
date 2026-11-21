"""Flat two-outlet control for the scalar pressure moment at a target outlet.

The upper 4D half-space has unit density, zero far velocity and no shear.
Each radius-one outlet prescribes the same smooth radial normal velocity.
The separately archived feeding fields supply only the comparison moments.
"""

import json
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss


RATE = 0.25
PEAK = RATE * 315 / (64 * np.pi)
FEED = Path("research/pde_ledger_v3/_measurements/mouth_front_pressure_dipole.json")


def unit_gauss(order):
    nodes, weights = leggauss(order)
    return (1 + nodes) / 2, weights / 2


def interior_rays(radius, order):
    """Integrate -grad(j)/(2 pi² |x-y|²) on rays starting at x."""
    t, wt = unit_gauss(order)
    cosine = np.r_[-t, t]
    r = radius[:, None]
    axial = r * cosine
    gap = (1 - r) * (1 + r)
    root = np.sqrt(axial**2 + gap)
    distance = np.where(cosine > 0, gap / (root + axial), root - axial)
    # The ray measure cancels the singular kernel. The remainder is a
    # polynomial of degree five in distance along each ray.
    s, ws = unit_gauss(3)
    result = np.zeros_like(distance)
    for node, weight in zip(s, ws, strict=True):
        step = distance * node
        shape = gap - step * (step + 2 * axial)
        result += weight * (r + step * cosine) * shape**2
    return (6 * PEAK / np.pi) * ((result * distance) @ np.tile(wt, 2))


def interior_shells(radius, order):
    """Project -grad(j) on source-centered shells and integrate their logarithm."""
    r = radius[:, None]
    nodes, weights = unit_gauss(order)
    velocity = np.zeros_like(radius)
    for span in (-r, 1 - r):
        offset = span * nodes**4  # resolve the logarithmic singularity at y=r
        y = r + offset
        a = r**2 + y**2
        z = 2 * r * y / a
        # Integral of t/(a-2ryt) dt over t in [-1,1]. Use the odd-power
        # expansion where the exact expression loses digits by cancellation.
        angular = 2 * (np.log((r + y) / abs(offset)) / z - 1) / (a * z)
        small = z < 0.025
        zs = z[small]
        angular[small] = (2 * zs / a[small]) * (
            1/3 + zs**2 * (1/5 + zs**2 * (1/7 + zs**2/9))
        )
        measure = abs(span) * 4 * nodes**3 * y**3 * (1 - y**2)**2
        velocity += (measure * angular) @ weights
    return (6 * PEAK / np.pi) * velocity


def exterior_series(radius, terms):
    """Differentiate the angularly integrated exterior Neumann potential."""
    field = np.zeros_like(radius)
    for k in range(terms):
        n = 2 * k
        moment = 48 / ((n+3) * (n+5) * (n+7) * (n+9))
        field += (n+2) / (n+1) * moment / radius**(n+3)
    return 2 * PEAK / np.pi * field


def exterior_shells(radius, order):
    """Integrate the observation gradient of the unsmoothed Neumann kernel."""
    y, wy = unit_gauss(order)
    t, wt = leggauss(order)
    r = radius[..., None]
    field = np.zeros_like(radius)
    for shell, weight in zip(y, wy, strict=True):
        squared = r**2 + shell**2 - 2 * r * shell * t
        angular = ((r - shell * t) / squared**2) @ wt
        field += weight * shell**2 * (1 - shell**2)**3 * angular
    return 2 * PEAK / np.pi * field


def archived_feed():
    run = json.loads(FEED.read_text())["runs"]["mouth_front_pressure_dipole02"]
    if run["guard_exit"] or run["child_exit"] or not run["limits_verified"]:
        raise ValueError("Archived feed run failed")
    rows = [json.loads(line) for line in run["stdout_lines"] if line.startswith("{")]
    return {row["d"]: row["moment"]["pressure"] for row in rows if row["grid"] == [96, 6, 24]}


def measure(order, d):
    # +xi points from target to companion, which is at (d,0,0).
    r, wr = unit_gauss(order)
    mu, wm = leggauss(order)
    axial = r[:, None] * mu
    transverse = r[:, None] * np.sqrt(1 - mu**2)
    dist = np.hypot(axial - d, transverse)
    weights = 2 * np.pi * (wr * r**2)[:, None] * wm
    normal = PEAK * (1 - r**2)**3
    source_normal = PEAK * np.maximum(1 - dist**2, 0)**3
    # Ray/series pressure uses the expanded Bernoulli difference.
    self_ray = interior_rays(r, order)[:, None]
    other_series = exterior_series(dist, 12)
    cross = -self_ray * other_series * (r[:, None] - d * mu) / dist
    other_pressure = -other_series**2 / 2
    cross_moment = np.sum(weights * axial * cross)
    source_square_moment = np.sum(weights * axial * other_pressure)
    # Independent radial velocities and Cartesian sum of squares, retaining
    # all four components (the two transverse axes combine by symmetry).
    self_shell = interior_shells(r, order)[:, None]
    other_shell = exterior_shells(dist, 16)
    vx = self_shell * mu + other_shell * (axial - d) / dist
    vt = self_shell * np.sqrt(1 - mu**2) + other_shell * transverse / dist
    pressure = -(vx**2 + vt**2 + (normal[:, None] + source_normal)**2) / 2
    pressure_alone = -(self_shell**2 + normal[:, None]**2) / 2
    shell_moment = np.sum(weights * axial * (pressure - pressure_alone))
    return dict(
        order=order, d=d, rate_per_outlet=float(np.sum(weights * normal[:, None])),
        companion_normal_speed_on_target_max=float(np.max(source_normal)),
        radial_velocity_difference_inside=float(np.max(abs(self_ray - self_shell))),
        radial_velocity_difference_outside=float(np.max(abs(other_series - other_shell))),
        cross_moment=float(cross_moment), companion_square_moment=float(source_square_moment),
        ray_series_moment=float(cross_moment + source_square_moment),
        shell_moment=float(shell_moment),
    )


def main():
    feed = archived_feed()
    bulk = {}
    slopes = {}
    for d in (8, 16):
        rows = [measure(order, d) for order in (24, 48, 96)]
        for row in rows:
            print(json.dumps(row))
        last = rows[-1]
        bulk[d] = last["ray_series_moment"]
        mirror = measure(48, -d)["shell_moment"]
        jump = feed[d] - bulk[d]
        slopes[d] = jump / (4 * np.pi / 3)
        print(json.dumps(dict(d=d, feed_moment=feed[d], bulk_moment=bulk[d],
                              reflected_global_x_moment=mirror, jump_moment=jump,
                              bulk_over_feed=bulk[d]/feed[d], sigma_s1=slopes[d])))
        for row in rows:
            np.testing.assert_allclose(row["rate_per_outlet"], RATE, atol=1e-14, rtol=0)
            assert row["companion_normal_speed_on_target_max"] == 0
        np.testing.assert_allclose(last["shell_moment"], bulk[d], rtol=2e-9, atol=0)
        np.testing.assert_allclose(rows[-2]["ray_series_moment"], bulk[d], rtol=2e-9, atol=0)
        np.testing.assert_allclose(rows[-2]["shell_moment"], -mirror, rtol=2e-9, atol=0)
        assert last["radial_velocity_difference_inside"] < 1e-9
        assert last["radial_velocity_difference_outside"] < 1e-15
    print(json.dumps(dict(bulk_doubling_ratio=bulk[16]/bulk[8],
                          feed_doubling_ratio=feed[16]/feed[8],
                          sigma_s1_doubling_ratio=slopes[16]/slopes[8])))


if __name__ == "__main__":
    main()
