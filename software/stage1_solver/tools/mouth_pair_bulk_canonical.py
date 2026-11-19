"""Equilibrium twists in 3D boxes at fixed constituent densities on the two faces.

The ideal mixture has a mobile bosonic atom on each face and one bosonic
cross-face dimer species. Its band minimum is -binding relative to two atoms.
Equilibrium sums over every composition consistent with the two atom counts.
Twists are imposed probes, with no independent dynamical gauge fields here.
"""

from dataclasses import dataclass

import numpy as np
from scipy.optimize import brentq
from scipy.special import ive, logsumexp, zeta


ATOM_HOP = 1.0
PAIR_HOP = 1/8
FILLING = 1/8


@dataclass
class Canonical:
    log_partition: np.ndarray
    twist_curvature: np.ndarray  # Total energy curvature per uniform bond angle


def single_species(length, number, hop, temperature, twist=0.0):
    """Canonical Bose trace by permutation cycles. twist is the total x boundary angle."""
    cycles = np.arange(1, number+1)
    momentum = 2*np.pi*np.arange(length)/length
    x = 2*hop*cycles/temperature
    exponent = x[:, None]*(np.cos(momentum)-1)
    log_axis = logsumexp(exponent, axis=1)
    shifted = x[:, None]*(np.cos(momentum+twist/length)-1)
    log_trace = 2*log_axis + logsumexp(shifted, axis=1)

    # At zero twist, winding numbers avoid subtracting large almost equal
    # diamagnetic and thermal terms. Check the winding cutoff against z itself.
    response = np.zeros(number)
    if twist == 0:
        winds = np.arange(1, int(np.ceil(15*np.sqrt(x[-1])/length))+12)
        weights = ive(length*winds[None, :], x[:, None])
        weight_sum = ive(0, x)+2*weights.sum(axis=1)
        assert np.allclose(length*weight_sum, np.exp(log_axis), atol=1e-13, rtol=3e-12)
        response = temperature*length**2*2*(weights @ (winds**2))/weight_sum

    log_z = np.zeros(number+1)
    curvature = np.zeros(number+1)
    for n in range(1, number+1):
        terms = log_trace[:n]+log_z[n-1::-1]
        norm = logsumexp(terms)
        log_z[n] = norm-np.log(n)
        if twist == 0:
            curvature[n] = np.exp(terms-norm) @ (response[:n]+curvature[n-1::-1])
    return Canonical(log_z, curvature)


def mixture(length, temperature, binding, face_angles=(0.0, 0.0)):
    volume = length**3
    number = int(volume*FILLING)
    assert number == volume*FILLING
    top, bottom = face_angles
    atom_t = single_species(length, number, ATOM_HOP, temperature, top)
    atom_b = atom_t if bottom == top else single_species(length, number, ATOM_HOP, temperature, bottom)
    pair = single_species(length, number, PAIR_HOP, temperature, top+bottom)
    populations = np.arange(number+1)
    logs = (pair.log_partition + binding*populations/temperature
            + atom_t.log_partition[::-1] + atom_b.log_partition[::-1])
    total = logsumexp(logs)
    probabilities = np.exp(logs-total)
    unpaired = probabilities @ (number-populations)/volume
    # All species/composition first derivatives vanish at zero twist.
    independent = (probabilities @ (atom_t.twist_curvature[::-1]+atom_b.twist_curvature[::-1]))/volume
    common = independent + 4*(probabilities @ pair.twist_curvature)/volume
    return -temperature*total, unpaired, common, independent


def excited_density(hop, temperature, chemical_offset, terms):
    """Infinite-lattice noncondensed density with a controlled band-edge tail."""
    n = np.arange(1, terms+1)
    slope = 2*hop/temperature
    density = np.sum(np.exp(n*chemical_offset/temperature)*ive(0, n*slope)**3)
    if chemical_offset == 0:
        leading = (2*np.pi*slope)**(-1.5)
        density += leading*(zeta(1.5, terms+1)+3/(8*slope)*zeta(2.5, terms+1)
                            + 33/(128*slope**2)*zeta(3.5, terms+1))
    return density


def thermodynamic(temperature, binding, terms=12000):
    # Equal face chemical potentials. At saturation the lowest of the atomic
    # threshold mu=0 and the dimer threshold 2 mu=-binding is reached first.
    assert binding != 0
    mu = min(0., -binding/2)
    atoms = excited_density(ATOM_HOP, temperature, mu, terms)
    pairs = excited_density(PAIR_HOP, temperature, 2*mu+binding, terms)
    excess = max(0., FILLING-atoms-pairs)
    atom_condensed = excess if binding < 0 else 0.
    pair_condensed = excess if binding > 0 else 0.
    saturation = atoms+pairs
    if saturation > FILLING:
        def density_residual(value):
            return (excited_density(ATOM_HOP,temperature,value,terms)
                    + excited_density(PAIR_HOP,temperature,2*value+binding,terms) - FILLING)
        mu = brentq(density_residual, mu-40*temperature, mu, xtol=5e-15)
        atoms = excited_density(ATOM_HOP,temperature,mu,terms)
    # The periodic Brillouin-zone derivative integrates to zero for the
    # thermal species. The k=0 condensate contributes its band curvature.
    opposite = 4*ATOM_HOP*atom_condensed
    common = opposite + 8*PAIR_HOP*pair_condensed
    return saturation, atom_condensed, pair_condensed, common, opposite, atoms+atom_condensed


def critical_temperature(binding):
    mu = min(0.,-binding/2)
    def critical_residual(temperature):
        return (excited_density(ATOM_HOP,temperature,mu,12000)
                + excited_density(PAIR_HOP,temperature,2*mu+binding,12000) - FILLING)
    return brentq(critical_residual,0.01,4.)


def verify_one_plus_one():
    # One constituent per face on eight sites has just two composition sectors.
    # Trace their explicit Bloch bands independently of the cycle recursion.
    modes = np.indices((2,2,2)).reshape(3,-1)*np.pi
    def boltzmann(hop, twist):
        shifted = modes.copy()
        shifted[0] += twist/2
        return np.exp(-2*hop*np.sum(1-np.cos(shifted),axis=0)/0.2).sum()
    top, bottom = 0.13,-0.07
    direct = boltzmann(ATOM_HOP,top)*boltzmann(ATOM_HOP,bottom)
    direct += np.exp(0.25/0.2)*boltzmann(PAIR_HOP,top+bottom)
    free = mixture(2,0.2,0.25,(top,bottom))[0]
    print("one-plus-one direct canonical free-energy residual:", free+0.2*np.log(direct))
    assert abs(free+0.2*np.log(direct)) < 1e-13


def verify_curvatures(length, temperature, binding):
    free, _, common, opposite = mixture(length,temperature,binding)
    for sign, expected in ((1,common),(-1,opposite)):
        estimates = []
        for theta in (0.04,0.02):
            forward = mixture(length,temperature,binding,(theta,sign*theta))[0]
            backward = mixture(length,temperature,binding,(-theta,-sign*theta))[0]
            estimates.append((forward+backward-2*free)/(theta**2*length))
        extrapolated = (4*estimates[1]-estimates[0])/3
        print(f"direct fixed-number free-energy curvature L={length} sign={sign}"
              f" trace={expected:.12g} finite-difference={extrapolated:.12g}")
        assert abs(expected-extrapolated) < 3e-8


def run(label, temperature, binding):
    print(f"{label}: t_atom={ATOM_HOP} t_pair={PAIR_HOP} T={temperature} binding={binding}"
          f" constituents per face per site={FILLING}")
    for length in (4,6,8,12,16,24,32):
        _, unpaired, common, opposite = mixture(length,temperature,binding)
        print(f"L={length:2d} N/face={int(length**3*FILLING):4d}"
              f" unpaired/face/site={unpaired:.12g} Y_common={common:.12g} Y_opposite={opposite:.12g}")
    result = thermodynamic(temperature,binding)
    assert np.allclose(result,thermodynamic(temperature,binding,3000),rtol=1e-10,atol=1e-13)
    print("bulk saturation density, atom condensate/face, pair condensate, Y_common, Y_opposite, unpaired/face/site:",
          " ".join(f"{x:.12g}" for x in result))
    verify_curvatures(6,temperature,binding)
    return result


if __name__ == "__main__":
    verify_one_plus_one()
    transition = critical_temperature(0.25)
    print(f"fixed density critical T at binding 0.25: {transition:.12g}")
    paired = run("below pair transition",0.12,0.25)
    normal = run("same Hamiltonian above pair transition",0.20,0.25)
    atomic = run("atomic-favored detuning control",0.12,-1.0)
    assert 0.12 < transition < 0.20
    assert paired[3] > 0 and paired[4] == 0
    assert normal[3] == normal[4] == 0
    assert atomic[3] == atomic[4] > 0
