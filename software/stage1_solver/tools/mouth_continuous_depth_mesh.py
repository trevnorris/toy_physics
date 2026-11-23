"""Move the joined terminal section without translating the flat bath interface.

The interior knots are a coordinate choice. A log-depth perturbation has a
rapidly vanishing mesh velocity towards the stationary plane w=0. No order
values or material points are pinned to interior knots.
"""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.special import expit

from mouth_breathing_mobile_patch import KAPPA, WELL, WIDTH
from mouth_exchange_carrier import mapped_nodes
from mouth_last_rings_fluid import GRIDS


def axes_for_depth(length, reference_depth, power=24):
    nr, _, nz = GRIDS["double"]
    width = .12
    end, origin = length-width, reference_depth-width
    if min(end, origin) <= 0:
        raise ValueError("Terminal section outside the open mesh chart")
    rr, default_z = mapped_nodes(nr, length, 28., 22.)
    s = np.linspace(0., 1., nr+1)
    # z^(-p) - z_ref^(-p) is constant for every positive interior knot.
    # At the reference geometry these are exactly the old uniform knots.
    lower = (origin*s)/(1. + np.expm1(power*np.log(origin/end))*s**power)**(1./power)
    lower[-1] = end
    near = np.linspace(end, length, nz+1)
    upper = np.linspace(length, length+width, nz+1)
    far = length+width+(22.-length-width)*(np.sinh(2*s)/np.sinh(2.))**2
    return rr, (default_z[0], np.r_[lower[:-1], near], np.r_[upper[:-1], far])


def neutral_force(length, reference, step, power):
    """False depth force on a neutral planar profile for two choices of lift."""
    def knots(l):
        _, (a, b, c) = axes_for_depth(l, reference, power)
        return np.concatenate([a, b[1:], c[1:]])

    left, original, right = (knots(length*np.exp(sign*step)) for sign in (-1, 0, 1))
    # Accumulate signed differences cell by cell, including the radial integral.
    # The force calculation in the physical problem transports the Q1 coefficients.
    # Resampling chi(w) gives a separate check independent of that lift.
    q = (1.+np.array([-1., 1.])/np.sqrt(3.))/2.
    def energies(w, chi):
        h = np.diff(w)
        sampled = chi[:-1,None]*(1.-q)+chi[1:,None]*q
        return np.stack((.5*KAPPA*np.diff(chi)**2/h,
                         h*np.mean(WELL*(sampled*(1.-sampled))**2, axis=1)))
    result = {}
    for lift in ("transported", "resampled"):
        cl, cr = (expit(-original/WIDTH),)*2 if lift == "transported" else (
            expit(-left/WIDTH), expit(-right/WIDTH))
        diff = (energies(left, cl)-energies(right, cr)).sum(axis=1)*((4*np.pi/3)*28**3/(2*step))
        result[lift] = dict(gradient=diff[0], well=diff[1], total=float(np.sum(diff)))
    return result


if __name__ == "__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("endpoint", type=Path)
    d=json.loads(p.parse_args().endpoint.read_text())
    if d["pole_radius"][0] != 4.:
        raise ValueError("The test requires the radius-four input")
    reference=d["depth"]
    print(json.dumps([dict(length=l, reference=reference, exponent=power, step=h,
                           **neutral_force(l,reference,h,power))
                      for l in (reference, .95*reference, 1.05*reference)
                      for power in (16, 24, 32) for h in (1.e-4, 5.e-5)], allow_nan=False))
