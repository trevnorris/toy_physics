"""Force on a movable shear-order profile and an ideal reflecting waist, per ray.

The smooth archived profile is sampled in thin uniform slices. The exterior
exponential continuation extends to infinity. Inertia and modulus both follow
order. A uniform upstream section changes length when the front moves.
"""

import argparse

import numpy as np
from scipy.optimize import brentq

from mouth_shear_tail_impedance import digest, emit, profile


class Ray:
    def __init__(self, mesh, cut, cells):
        sample, cross, self.source = profile(1.0, mesh)
        self.start, self.half, self.end = (cross(level) for level in (.99, .5, cut))
        self.ell = (self.end - cross(10 * cut)) / np.log(10.)
        self.cut = cut
        self.step = (self.end - self.start) / cells
        self.layers = sample(self.start + self.step * (np.arange(cells) + .5))
        if np.min(self.layers) <= 0 or self.start <= .52:
            raise ValueError("Profile or available ordered lead is invalid")

    def mode(self, omega, gap_change, boundary, with_energy=False):
        k = np.asarray(omega)
        decay = np.sqrt(1.0 - (2.0 * self.ell * k)**2) / (2.0 * self.ell)
        # Finite-energy branch of (chi y')' + omega^2 chi y = 0.
        beta = k**2 / (1.0 / (2.0 * self.ell) + decay)
        y = np.ones_like(k)
        traction = self.cut * beta
        energy = self.cut * (k**2 + beta**2) / (8.0 * decay) if with_energy else None
        cos, sin = np.cos(k * self.step), np.sin(k * self.step)
        for chi in self.layers[::-1]:
            if with_energy:
                energy += .25 * self.step * (chi * (k*y)**2 + traction**2/chi)
            y, traction = (cos*y - sin*traction/(k*chi),
                           k*chi*sin*y + cos*traction)
        lead = .99
        length = self.start - .5 + gap_change
        pressure = .25 * (lead*(k*y)**2 + traction**2/lead)
        cos, sin = np.cos(k*length), np.sin(k*length)
        y, traction = cos*y - sin*traction/(k*lead), k*lead*sin*y + cos*traction
        if with_energy:
            action = (energy + length*pressure)/omega
            return pressure/action  # momentum flux at the ideal mirror, per unit action
        return y if boundary == "fixed" else traction/k

    def roots(self, gap_change, boundary):
        omega = np.linspace(1e-6, (1.0 - 1e-10) / (2.0*self.ell), 801)
        values = self.mode(omega, gap_change, boundary)
        flips = np.flatnonzero(np.signbit(values[:-1]) != np.signbit(values[1:]))
        residual = lambda w: self.mode(w, gap_change, boundary)
        return [brentq(residual, omega[i], omega[i+1], xtol=2e-14) for i in flips]


def main(args):
    ray = Ray(args.mesh, args.cut, args.cells)
    emit(test="preparation", mesh=args.mesh, cells=args.cells, cut=args.cut,
         source=str(ray.source), sha256=digest(ray.source), lead_order=.99, waist=.5,
         lead_start=ray.start, half_order=ray.half, tail_join=ray.end, decay_length=ray.ell,
         cutoff=1/(2*ray.ell), mass_law="rho_T=chi; shear modulus=chi")
    for boundary in ("fixed", "free"):
        # Track each positive-frequency bound root through all signed changes.
        # The zero-frequency uniform displacement for the free end is excluded.
        shifts = [-.02, -.01, 0., .01, .02]
        roots = [ray.roots(q, boundary) for q in shifts]
        counts = [len(row) for row in roots]
        emit(test="bound_counts", boundary=boundary, shifts=shifts, counts=counts)
        for branch, frequencies in enumerate(zip(*roots), start=1):
            forces = [ray.mode(w, q, boundary, with_energy=True)
                      for q, w in zip(shifts, frequencies)]
            full_work = -(frequencies[4] - frequencies[0])/.04
            half_work = -(frequencies[3] - frequencies[1])/.02
            work = (4.0*half_work-full_work)/3.0
            emit(test="fixed_action", boundary=boundary, branch=branch, omega=frequencies,
                 toward_bulk_per_action=forces, waist_force_per_action=-forces[2],
                 center_work_full=full_work, center_work_half=half_work, center_work_extrap=work,
                 tangent_force_slope=(forces[3]-forces[1])/.02,
                 empty_wave_force=0.)
            if abs(work-forces[2]) > .005 * max(1.,abs(work)):
                raise ValueError("Energy derivative and local momentum flux disagree")
        if min(counts) < 1:
            raise ValueError("No bound branch could be followed")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mesh", type=int, required=True)
    parser.add_argument("--cut", type=float, required=True)
    parser.add_argument("--cells", type=int, required=True)
    main(parser.parse_args())
