"""Three-dimensional Brillouin-zone integrals for an ideal atom/dimer mixture.

Each face has fixed filling NU counting free atoms and one atom per dimer.
The two free atom bands coincide. Each has chemical potential mu, and each
dimer consumes one particle from each face. Its band bottom is -binding.
Temperature, hoppings and energy are expressed in the same lattice units.
"""

from dataclasses import dataclass

import numpy as np
from scipy.optimize import brentq
from scipy.special import roots_legendre


NU = 1 / 8
ATOM_HOP = 1.0
DIMER_HOP = 1 / 8


class CubicZone:
    def __init__(self, order):
        nodes, weights = roots_legendre(order)
        nodes = (nodes + 1) / 2
        weights = weights / 2
        # Even and cubic integrands: in one octant, take x as the largest
        # coordinate and multiply by three. k=(pi*r,pi*r*u,pi*r*v).
        # The normalized measure 3 r^2 dr du dv cancels the Bose singularity.
        r = nodes[:, None, None]
        kx = np.pi * r
        ky = kx * nodes[None, :, None]
        kz = kx * nodes[None, None, :]
        self.weights = 3 * r**2 * (
            weights[:, None, None] * weights[None, :, None] * weights[None, None, :]
        )
        self.energy = 4 * sum(np.sin(k / 2)**2 for k in (kx, ky, kz))
        # Average the three diagonal entries for a single Cartesian stiffness
        # before restricting to the x-largest subdomain.
        self.curvature = 2 * sum(np.cos(k) for k in (kx, ky, kz)) / 3
        self.velocity_squared = 4 * sum(np.sin(k)**2 for k in (kx, ky, kz)) / 3

    def occupation(self, temperature, hop, gap):
        x = (hop * self.energy + gap) / temperature
        return np.exp(-x) / -np.expm1(-x)

    def integrate(self, values):
        return float(np.sum(self.weights * values))

    def excited_density(self, temperature, hop, gap):
        return self.integrate(self.occupation(temperature, hop, gap))

    def thermal_curvature(self, occupation, temperature, hop):
        diamagnetic = self.integrate(hop * self.curvature * occupation)
        current = -self.integrate(
            hop**2 * self.velocity_squared * occupation * (1 + occupation) / temperature
        )
        return diamagnetic, current


@dataclass
class Mixture:
    zone: CubicZone
    binding: float

    @property
    def saturation_mu(self):
        return min(0.0, -self.binding / 2)

    def total_excited(self, temperature, mu):
        return self.zone.excited_density(temperature, ATOM_HOP, -mu) + (
            self.zone.excited_density(temperature, DIMER_HOP, -self.binding - 2 * mu)
        )

    def critical_temperature(self):
        def excess(temperature):
            return self.total_excited(temperature, self.saturation_mu) - NU

        return brentq(excess, 0.02, 4.0, xtol=5e-15)

    def equilibrium(self, temperature):
        mu = self.saturation_mu
        if self.total_excited(temperature, mu) > NU:
            mu = brentq(lambda value: self.total_excited(temperature, value) - NU,
                        mu - 8.0, mu, xtol=5e-15)
        atom = self.zone.occupation(temperature, ATOM_HOP, -mu)
        dimer = self.zone.occupation(temperature, DIMER_HOP, -self.binding - 2 * mu)
        na = self.zone.integrate(atom)
        nd = self.zone.integrate(dimer)
        condensate = max(NU - na - nd, 0.0) if mu == self.saturation_mu else 0.0
        # Each tested binding has one strictly lower threshold in mu.
        n0a = condensate if self.binding < 0 else 0.0
        n0d = condensate if self.binding > 0 else 0.0
        da, ca = self.zone.thermal_curvature(atom, temperature, ATOM_HOP)
        dd, cd = self.zone.thermal_curvature(dimer, temperature, DIMER_HOP)
        # Each atomic species responds to its own twist, while the dimer
        # responds to their sum. These include both thermal-current terms.
        ya = 2 * ATOM_HOP * n0a + da + ca
        yd = 2 * DIMER_HOP * n0d + dd + cd
        return (mu, na, nd, n0a, n0d, 2 * ya + 4 * yd, 2 * ya,
                da, ca, dd, cd, na + nd + n0a + n0d - NU)


def main():
    lower_order = {}
    for order in (24, 48, 72):
        zone = CubicZone(order)
        print(f"Gauss order={order} normalized zone volume={zone.integrate(1.0):.15g}")
        for binding in (0.25, -1.0):
            mixture = Mixture(zone, binding)
            tc = mixture.critical_temperature()
            print(f"binding={binding:+.2f} saturation mu={mixture.saturation_mu:.8f} Tc={tc:.14f}")
            for temperature in (0.12, 0.20):
                result = mixture.equilibrium(temperature)
                mu, na, nd, n0a, n0d, same, opposite, da, ca, dd, cd, balance = result
                print(f"  T={temperature:.2f} mu={mu:.14g} excited(atom,dimer)=({na:.14g},{nd:.14g}) "
                      f"condensed(atom,dimer)=({n0a:.14g},{n0d:.14g}) "
                      f"same,opposite=({same:.14g},{opposite:.14g})")
                print(f"    thermal diamagnetic,current: atom=({da:.14g},{ca:.14g}) "
                      f"dimer=({dd:.14g},{cd:.14g}) balance={balance:.3g}")
                vector = np.array((tc, *result))
                key = (binding, temperature)
                if order == 72:
                    diff = np.max(np.abs(vector - lower_order[key]))
                    print(f"    max absolute change from order 48={diff:.6g}")
                    if diff > 1e-9 or max(abs(balance), abs(da + ca), abs(dd + cd)) > 1e-9:
                        raise ValueError("Unconverged Brillouin-zone integral or number balance")
                lower_order[key] = vector


if __name__ == "__main__":
    main()
