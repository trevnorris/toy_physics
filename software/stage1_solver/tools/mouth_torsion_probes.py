"""Free mouth-site orientations and compare two mechanical torque impulses."""

from itertools import product
from math import asin, cos, pi, sin

import numpy as np

from mouth_torsion_power import MASS, spectrum


def finite_site_test():
    """Energy bounds in units of the individual face contact stiffness."""
    site_strength, link_strength, tilt = 400, 4, .1
    site_cost = site_strength * (1 - cos(tilt))
    overlap_error = 4 * sin(tilt / 4)
    for radius in (5, 10):
        trial = (72 + 48 / radius + 12 / radius**2) / (16 * (radius - 1))
        h = cos(1 / (2 * radius))
        slip = (pi - 1 / (2 * radius**2)) / 4 - asin(overlap_error / h)
        bad_bonds = 2 * link_strength * (h - overlap_error) * (1 - cos(slip))
        print(f"R={radius} free site strength={site_strength} "
              f"trial_upper={trial:.12f} single_tilt_floor={site_cost:.12f} "
              f"wrong_rim_two_link_floor={bad_bonds:.12f} slip_radians={slip:.9f}")
        assert slip > 0 and trial < min(site_cost, bad_bonds)


def impulses(k):
    # Fourier transform of the physical impulses on the four shared edge
    # rotors bounding the positively oriented xy square at the origin.
    kick = np.zeros((15, 2), complex)
    kick[0, 0] = 1
    kick[:2, 1] = [1 - np.exp(-1j * k[1]), np.exp(-1j * k[0]) - 1]
    return kick


def classify(k):
    values, coordinates = spectrum(k)
    kick = impulses(k)
    modal = abs(coordinates.conj().T @ kick) ** 2
    norm = np.sum(abs(kick) ** 2 / MASS[:, None], axis=0)
    assert np.allclose(modal.sum(axis=0), norm, atol=1e-13)
    return values, modal, norm


def main():
    finite_site_test()
    k = .04 * np.array([1., 2., 3.]) / np.sqrt(14)
    values, modal, norm = classify(k)
    for i, name in enumerate(("single edge", "oriented square")):
        print(f"one harmonic k=0.04 ray=123 probe={name} "
              f"zero_fraction={sum(modal[values < 1e-9, i])/norm[i]:.12g} "
              f"two_acoustic_fraction={sum(modal[4:6,i])/norm[i]:.12g}")
    for side in (5, 9, 15):
        locked = np.zeros(2)
        total = np.zeros(2)
        for k in product(2 * pi * np.fft.fftfreq(side), repeat=3):
            values, modal, norm = classify(k)
            locked += modal[values < 1e-9].sum(axis=0)
            total += norm
        fraction = locked / total
        print(f"finite periodic L={side} localized impulse zero-frequency kinetic fraction "
              f"single_edge={fraction[0]:.12f} oriented_square={fraction[1]:.3g}")
        assert .05 < fraction[0] < .5 and fraction[1] < 1e-20


if __name__ == "__main__":
    main()
