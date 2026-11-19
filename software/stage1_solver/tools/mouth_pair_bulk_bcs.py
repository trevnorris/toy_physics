"""Independent face twists in a balanced 3D attractive Fermi-Hubbard model.

Uniform onsite pairing and Hartree mean field, t = lattice spacing = k_B = 1.
The two fermion species are labelled by face. Each has one orbital per site;
their numbers are separately conserved by the microscopic Hamiltonian. The
uniform saddle is taken on the zero-winding branch near zero total twist.
The two temperatures use identical hopping, attraction and fixed mean filling.
Finite periodic cubes discretize this saddle, not the exact many-body Gibbs state.
"""

import numpy as np
from scipy.optimize import brentq
from scipy.special import expit, logsumexp


U = 4.0                 # positive magnitude of the onsite attraction
FILLING = 0.8           # both face species combined, half filling would be 1
TEMPERATURES = (0.15, 1.0)
LENGTHS = (6, 8, 12, 16, 24, 32, 48, 64)


class Cube:
    def __init__(self, length):
        # Keep all kx. Collapse only exact ky,kz reflection/permutation orbits.
        half = np.arange(length // 2 + 1)
        cosines = np.cos(2 * np.pi * half / length)
        degeneracy = np.where((half == 0) | (2 * half == length), 1, 2)
        y, z = np.triu_indices(len(half))
        self.transverse = -2 * (cosines[y] + cosines[z])
        self.weights = degeneracy[y] * degeneracy[z] * np.where(y == z, 1, 2)
        self.volume = length ** 3
        self.kx = (2 * np.pi * np.arange(length) / length)[:, None]
        assert length * self.weights.sum() == self.volume

    def average(self, values):
        return np.sum(values * self.weights) / self.volume

    def bands(self, top_angle, bottom_angle):
        top = self.transverse - 2 * np.cos(self.kx + top_angle)
        bottom = self.transverse - 2 * np.cos(-self.kx + bottom_angle)
        return (top + bottom) / 2, (top - bottom) / 2


def thermal_kernel(xi, mismatch, gap, temperature):
    """Return [1-f(E+d)-f(E-d)]/(2E), using its continuous value at E=0."""
    energy = np.hypot(xi, gap)
    top = expit(-(energy + mismatch) / temperature)
    bottom = expit(-(energy - mismatch) / temperature)
    f = expit(-mismatch / temperature)
    kernel = f * (1 - f) / temperature
    return np.divide(1 - top - bottom, 2 * energy, out=kernel,
                     where=energy > 1e-8 * temperature)


def equilibrium(cube, top_angle, bottom_angle, temperature):
    average_band, mismatch = cube.bands(top_angle, bottom_angle)

    def chemical_potential(gap):
        def excess(mu):
            xi = average_band - mu
            kernel = thermal_kernel(xi, mismatch, gap, temperature)
            return cube.average(1 - 2 * xi * kernel) - FILLING

        return brentq(excess, -12, 12, xtol=2e-13)

    def pairing_residual(gap, mu):
        kernel = thermal_kernel(average_band - mu, mismatch, gap, temperature)
        return 1 / U - cube.average(kernel)

    normal_mu = chemical_potential(0)
    normal_residual = pairing_residual(0, normal_mu)
    if normal_residual >= 0:
        gap, mu = 0., normal_mu
    else:
        gap = brentq(lambda d: pairing_residual(d, chemical_potential(d)),
                     0, U, xtol=2e-13)
        mu = chemical_potential(gap)
    xi = average_band - mu
    filling = cube.average(1 - 2 * xi * thermal_kernel(xi, mismatch, gap, temperature))
    return gap, mu, filling, gap * pairing_residual(gap, mu), normal_mu, normal_residual


def fock_free_energy(cube, top, bottom, temperature, gap, mu):
    """Fixed-filling Helmholtz energy/site from an explicit four-state Fock trace."""
    # For each k, |empty> and |top,bottom> mix, whereas either single is diagonal.
    # Assemble from individual hopping dispersions without the quasiparticle kernel.
    top_xi = cube.transverse - 2 * np.cos(cube.kx + top) - mu
    bottom_xi = cube.transverse - 2 * np.cos(-cube.kx + bottom) - mu
    block = np.zeros((*top_xi.shape, 2, 2))
    block[..., 0, 1] = block[..., 1, 0] = -gap
    block[..., 1, 1] = top_xi + bottom_xi
    even = np.linalg.eigvalsh(block)
    levels = np.stack((even[..., 0], even[..., 1], top_xi, bottom_xi), axis=-1)
    grand = gap ** 2 / U - temperature * cube.average(logsumexp(-levels / temperature, axis=-1))
    # mu is the Hartree-shifted chemical potential. At equal filling per species
    # the uniform Hartree energy is constant throughout either twist preparation.
    return grand + mu * FILLING - U * FILLING ** 2 / 4


def curvature(cube, temperature, gap, mu):
    # Inversion makes the twist derivatives of the density and gap equations
    # vanish at the origin, so the explicit second derivative is canonical here.
    band, _ = cube.bands(0., 0.)
    xi = band - mu
    energy = np.hypot(xi, gap)
    speed_squared = (2 * np.sin(cube.kx)) ** 2
    band_curvature = 2 * np.cos(cube.kx)
    f = expit(-energy / temperature)
    fprime = -f * (1 - f) / temperature
    if gap == 0:
        diam = cube.average(band_curvature * 2 * expit(-xi / temperature))
        param = cube.average(2 * speed_squared * fprime)
        return diam, param, param, 0.
    tanh = np.tanh(energy / (2 * temperature))
    occupation = 1 - xi * tanh / energy
    dn_dxi = (2 * fprime * (xi / energy) ** 2
              - tanh * gap ** 2 / energy ** 3)
    diam = cube.average(band_curvature * occupation)
    common_param = cube.average(2 * speed_squared * fprime)
    opposite_param = cube.average(speed_squared * dn_dxi)
    difference = cube.average(speed_squared * gap ** 2 / energy ** 3
                              * (tanh + 2 * energy * fprime))
    return diam, common_param, opposite_param, difference


def curvature_from_fock(cube, temperature, sign, f0):
    # Vary a constituent phase per x bond. Resolve density and gap at all four
    # perturbations and differentiate the spectrum-based free energy.
    def central(h):
        f = []
        for angle in (-h, h):
            gap, mu, *_ = equilibrium(cube, angle, sign * angle, temperature)
            f.append(fock_free_energy(cube, angle, sign * angle, temperature, gap, mu))
        return ((f[0] - f0) + (f[1] - f0)) / h ** 2

    h = 0.006
    return (4 * central(h / 2) - central(h)) / 3


def immobile_control():
    # Exact grand-canonical onsite trace, with fixed mean filling. With t=0 the
    # microscopic Hamiltonian has no twist-bearing term in any dimension.
    temperature = TEMPERATURES[0]
    occupation = np.array([0., 1., 1., 2.])
    interaction = np.array([0., 0., 0., -U])

    def probabilities(mu):
        logits = (mu * occupation - interaction) / temperature
        return np.exp(logits - logsumexp(logits))

    mu = brentq(lambda m: probabilities(m) @ occupation - FILLING, -12., 12.)
    p = probabilities(mu)
    print(f"atomic t=0 T={temperature:g} mu={mu:.12g} n={(p @ occupation):.12g}"
          f" pair_probability={p[3]:.12g} single_probability={p[1]+p[2]:.12g}")
    assert abs(p @ occupation - FILLING) < 1e-11


def main():
    print("3D periodic cubic uniform Hartree-BCS saddle; t=a=kB=1,"
          f" U={U:g}, total n={FILLING:g}, mu_eff=mu_physical+Un/2")
    print("C,R=d2(F/site)/dA2 for (A,+A),(A,-A) per x bond;"
          " for total angle theta, D=L^-1*d2F_total/dtheta2")
    print("T L Delta mu_eff filling Delta*gap_residual normal_inverse_susceptibility"
          " Fpaired_minus_Fnormal diam paramC paramR C R CminusR_identity")
    for temperature in TEMPERATURES:
        for length in LENGTHS:
            cube = Cube(length)
            gap, mu, filling, residual, normal_mu, linear = equilibrium(cube, 0., 0., temperature)
            f0 = fock_free_energy(cube, 0., 0., temperature, gap, mu)
            normal = fock_free_energy(cube, 0., 0., temperature, 0., normal_mu)
            diam, pc, pr, difference = curvature(cube, temperature, gap, mu)
            common, relative = diam + pc, diam + pr
            values = (gap, mu, filling, residual, linear, f0 - normal, diam, pc, pr,
                      common, relative, difference)
            print(f"{temperature:g} {length}", *(f"{v:.12g}" for v in values), flush=True)
            assert abs(filling - FILLING) < 1e-12 and abs(residual) < 1e-12
            assert abs(common - relative - difference) < 2e-13
            assert f0 <= normal + 3e-13
            if length in (6, 16, 32):
                for sign, analytic in ((1, common), (-1, relative)):
                    numeric = curvature_from_fock(cube, temperature, sign, f0)
                    print(f"fock_FD T={temperature:g} L={length} bottom_sign={sign:+d}"
                          f" analytic={analytic:.12g} numeric={numeric:.12g}"
                          f" residual={numeric-analytic:.12g}", flush=True)
                    assert abs(numeric - analytic) < 1e-7
    immobile_control()


if __name__ == "__main__":
    main()
