"""Integrate the finite source and target directly, with no point-source substitution."""

import numpy as np
from numpy.polynomial.legendre import leggauss


def rings(radius, profile, order):
    nodes, weights = leggauss(order)
    r = radius * (1 + nodes[:, None]) / 2
    axial = r * nodes[None, :]
    transverse = r * np.sqrt(1 - nodes[None, :]**2)
    # Integrate full spherical rings. Each profile has unit total volume rate.
    volume_weights = 2 * np.pi * r**2 * (radius * weights[:, None] / 2) * weights[None, :]
    density = 3 / (4 * np.pi * radius**3)
    if profile == "taper":
        density *= 2.5 * (1 - (r/radius)**2)
    elif profile == "dipole":
        density *= 1 + axial/(2 * radius)
    elif profile != "uniform":
        raise ValueError(profile)
    return axial.ravel(), transverse.ravel(), (volume_weights * density).ravel()


def pair_coefficient(target, source, separation, order):
    # Integrate (source_x - target_x) / |source - target|^3 over two
    # normalized mouths. Axial symmetry leaves one relative azimuth.
    x_t, h_t, w_t = target
    x_s, h_s, w_s = source
    dx = separation + x_s[None, :] - x_t[:, None]
    squared = dx**2 + h_t[:, None]**2 + h_s[None, :]**2
    cross = 2 * h_t[:, None] * h_s[None, :]
    weight = w_t[:, None] * w_s[None, :]
    cosines = np.cos(2 * np.pi * np.arange(2 * order) / (2 * order))
    axial = [np.sum(weight * dx / (squared - cross * cosine)**1.5) for cosine in cosines]
    return np.mean(axial)


def main():
    # Output 4*pi*Delta F_x/(rho*Q_t*Q_s) in the length units of this test.
    for order in (10, 18):
        for target_name, source_name in (("uniform", "taper"), ("taper", "uniform"),
                                         ("dipole", "taper")):
            target = rings(1, target_name, order)
            source = rings(1.5, source_name, order)
            readings = [pair_coefficient(target, source, d, order) for d in (4, 8)]
            rates = [target[2].sum(), source[2].sum()]
            print(f"order={order} profiles={target_name},{source_name} "
                  f"normalized rates={rates[0]:.12f},{rates[1]:.12f} "
                  f"scaled loads(d=4,8)={readings[0]:.12f},{readings[1]:.12f}", flush=True)
            if max(abs(rate - 1) for rate in rates) > 1e-12:
                raise ArithmeticError("Rate changed with the profile or quadrature")


if __name__ == "__main__":
    main()
