"""Survival of an initially local rim excitation coupled to a bulk band.

One excitation, initially on the rim, evolves under a time-independent
Hermitian hopping law. This is an added receiver model. Energy and time
are measured in units of the bulk hopping and its inverse (hbar=1).
The continuum integral is for a true half-line. A finite-chain calculation
provides an independent short-time check, not a late-time extrapolation.
"""
from fractions import Fraction

import numpy as np
from scipy.integrate import quad_vec
from scipy.linalg import eigh_tridiagonal
from scipy.optimize import brentq
from scipy.special import j0, j1


HOPPING = 1.0
BULK_CENTER = 4.0
TIMES = np.array([0., .33, .66, 1., 2., 4., 8., 16., 32., 64.])


def equal_hopping_pole(offset):
    """For g=J, a tail proportional to lambda**n fixes lambda=J/offset."""
    if abs(offset) <= HOPPING:
        return None, 0.
    decay = HOPPING / offset
    return offset + HOPPING * decay, 1 - decay**2


def bulk_measure(angle, offset):
    # The retarded terminal self-energy is J exp(-ik), E-E0=2 J cos k.
    # At |offset|=J use the cancelled expression at the threshold.
    if abs(offset) == HOPPING:
        return (1 + np.sign(offset) * np.cos(angle)) / np.pi
    denominator = abs(HOPPING * np.exp(1j * angle) - offset)**2
    return 2 * HOPPING**2 * np.sin(angle)**2 / (np.pi * denominator)


def spectral_response(offset):
    def integrand(angle):
        energy = 2 * HOPPING * np.cos(angle)
        moment_and_propagator = np.r_[BULK_CENTER + energy, np.exp(-1j * energy * TIMES)]
        return bulk_measure(angle, offset) * moment_and_propagator

    integrated, _ = quad_vec(integrand, 0., np.pi, epsabs=2e-13, epsrel=2e-13)
    continuum_energy, amplitude = integrated[0].real, integrated[1:]
    pole, weight = equal_hopping_pole(offset)
    bound_energy = 0.
    if pole is not None:
        amplitude += weight * np.exp(-1j * pole * TIMES)
        bound_energy = weight * (BULK_CENTER + pole)
    return amplitude, continuum_energy, bound_energy


def finite_chain(offset, contact):
    diagonal = np.full(360, BULK_CENTER)
    diagonal[0] += offset
    off_diagonal = np.full(len(diagonal)-1, HOPPING)
    off_diagonal[0] = contact
    energy, modes = eigh_tridiagonal(diagonal, off_diagonal)
    weights = abs(modes[0])**2
    return energy, weights, weights @ np.exp(-1j * np.outer(energy - BULK_CENTER, TIMES))


def main():
    print("one incoming local excitation; J=1, hbar=1, bulk band=[2,6]; "
          "amplitudes below omit the common phase exp(-4 i t).")
    for offset in (0., 1., 1.5, -1.5):
        amp, continuum_energy, bound_energy = spectral_response(offset)
        where, weight = equal_hopping_pole(offset)
        spectrum, local_weights, observed = finite_chain(offset, HOPPING)
        isolated = (spectrum < BULK_CENTER-2*HOPPING-1e-9) | (spectrum > BULK_CENTER+2*HOPPING+1e-9)
        pole_energy = None if where is None else where + BULK_CENTER
        energy_error = abs(bound_energy+continuum_energy-(BULK_CENTER+offset))
        print(f"g=1 offset={offset:g} pole_absolute={pole_energy} total_bound={weight:.12g} "
              f"late_local={weight**2:.12g} bound_energy={bound_energy:.12g} "
              f"continuum_energy={continuum_energy:.12g} energy_residual={energy_error:.3g} "
              f"normalization_residual={abs(1-amp[0]):.3g} finite_chain_max_amplitude_error="
              f"{max(abs(observed-amp)):.3g} finite_chain_bound_weight={sum(local_weights[isolated]):.12g}")
        for time, local in zip(TIMES, abs(amp)**2):
            print(f"delay={time:g} offset={offset:g} local_population={local:.12g}")
        assert spectrum[0] > 0 and energy_error < 1e-12 and abs(amp[0]-1)<1e-12
        assert max(abs(observed-amp)) < 2e-12 and abs(sum(local_weights[isolated])-weight)<1e-12
        if where is not None:
            assert isolated.sum() == 1 and abs(spectrum[isolated][0]-pole_energy) < 1e-12
        if offset == 0.:
            exact = np.r_[1., j1(2 * TIMES[1:]) / TIMES[1:]]
            first_half = brentq(lambda t: (j1(2*t)/t)**2 - .5, .01, 1.)
            # Plancherel: pi integral rho(E)^2 dE for the semi-circle measure.
            residence = 8 / (3 * np.pi * HOPPING)
            print(f"unshifted Bessel_max_error={max(abs(amp-exact)):.3g} "
                  f"first_half_population_time={first_half:.12g} integrated_local_population={residence:.12g}")
            assert max(abs(amp-exact))<1e-12
        if offset == HOPPING:
            threshold = j0(2*TIMES)-1j*j1(2*TIMES)
            print(f"at_threshold Bessel_max_error={max(abs(amp-threshold)):.3g} "
                  "no normalizable bound state; late probability ~1/(pi J t).")
            assert max(abs(amp-threshold)) < 1e-12

    bare_offset = Fraction(3, 2)
    tail = 1 / bare_offset
    retained = 1 - tail**2
    print(f"g=1 offset=3/2 exact_tail_ratio={tail} pole_absolute={4+bare_offset+tail} "
          f"total_bound={retained} late_local={retained**2} bound_cloud={retained*(1-retained)} "
          f"total_escape={1-retained}")
    _, _, isolated_amp = finite_chain(float(bare_offset), 0.)
    print(f"g=0 site_isolated offset=1.5 max_survival_error={max(abs(abs(isolated_amp)**2-1)):.3g}")
    assert max(abs(abs(isolated_amp)**2-1))<1e-12


if __name__ == "__main__":
    main()
