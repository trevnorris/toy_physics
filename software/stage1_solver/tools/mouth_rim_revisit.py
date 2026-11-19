"""Chance of finding a locally deposited quantum when the next parcel arrives.

The parcel arrivals sample the unprobed survival probability. The added bulk
branch is an illustrative energy-receiving channel with no transverse shear.
Only one excitation is present, including on the branched channel geometry.
"""
import argparse
from pathlib import Path
import re

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq


BOUND = 37 / 6
Z = 5 / 9
RIM = 5.5


def surface(z):
    # For Im z>0 the two principal roots give the retarded branch on both bands.
    return 2 / (z + np.sqrt(z - 2 + 0j) * np.sqrt(z + 2 + 0j))


def resolvent(energy, coupling):
    return 1 / (energy - RIM - surface(energy - 4) - coupling**2 * surface(energy - 6))


def next_arrival(rate, coupling):
    def integrand(energy):
        return rate * abs(resolvent(energy + .5j * rate, coupling))**2 / (2 * np.pi)
    # Split at both bands and at the original surface pole. The Laplace frequency
    # keeps the integrand finite even for zero coupling to the second channel.
    breaks = [-np.inf, 2., 4., 6., BOUND, 8., np.inf]
    pieces = [quad(integrand, a, b, epsabs=3e-12, epsrel=3e-10)
              for a, b in zip(breaks, breaks[1:])]
    return sum(value for value, _ in pieces), sum(error for _, error in pieces)


def one_excitation_spectrum(coupling):
    # G has a point mass Z at BOUND only when the second port is shut.
    edges = [2., 4., 6., BOUND, 8.]
    def density(energy):
        return -resolvent(energy + 0j, coupling).imag / np.pi
    if coupling == 0:
        continuum = quad(density, 2., 6., epsabs=3.e-12)[0]
        return continuum, Z
    continuum = sum(quad(density, a, b, epsabs=5.e-12, epsrel=1e-10)[0]
                    for a, b in zip(edges, edges[1:]))
    return continuum, 0.


def test(rates):
    for coupling in (0., .025, .05, .1, .2):
        continuum, pole = one_excitation_spectrum(coupling)
        # With both ports open, the outer-edge denominators exclude exterior poles.
        low = (1 / resolvent(2. + 0j, coupling)).real
        high = (1 / resolvent(8. + 0j, coupling)).real
        print(f"spectrum second_hop={coupling:g} continuous_weight={continuum:.12g} "
              f"pole_weight={pole:.12g} lower_edge_den={low:.12g} upper_edge_den={high:.12g}")
        assert low < 0 < high and abs(continuum + pole - 1) < 1e-9
    for rate in rates:
        base, _ = next_arrival(rate, 0.)
        for coupling in (0., .025, .05, .1, .2):
            exact, error = next_arrival(rate, coupling)
            width = -2 * Z * coupling**2 * surface(BOUND-6).imag
            estimate = Z**2 * rate / (rate + width)
            print(f"revisit rate={rate:.12g} mean_wait={1/rate:.12g} second_hop={coupling:g} "
                  f"mean_local={exact:.12g} fraction_of_one_band={exact/base:.12g} "
                  f"weak_width={width:.12g} pole_only_estimate={estimate:.12g} quad_error={error:.2g}")
            assert 0 < exact < 1 and error < 1e-8
        h = brentq(lambda b: next_arrival(rate, b)[0] - base / 2, 0, .25)
        print(f"half_reference rate={rate:.12g} second_hop={h:.12g} "
              f"mean_local={next_arrival(rate,h)[0]:.12g}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("earlier_rate_output", type=Path)
    source = parser.parse_args().earlier_rate_output.read_text()
    rates = [float(x) for x in re.search(r"mesh level=2 edge=132 .*?rates=\[([^\]]+)\]", source).group(1).split(",")]
    test(rates)
