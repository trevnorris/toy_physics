"""Test direct phase matching of the saved carrier to an ideal open 4D receiver.

Use a flat aperture into w>0, zero normal flow elsewhere, constant unit
receiver number density, and m=hbar=1. The coordinate contrast used below
annihilates both a constant potential and any uniform far velocity.
"""

import json
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.special import erf

from mouth_carrier_momentum import read_state
from mouth_release_input_check import manifest_fields


def cell_kernel_surface(offset, dx, order):
    # Integral of 1/r^2 over each source cell by div(r/r^2) = 1/r^2.
    # Probes sit on grid vertices. Faces through a probe have zero normal
    # projection, and an infinitesimal sphere around it contributes zero.
    nodes, weights = leggauss(order)
    total = np.zeros(len(offset))
    for axis in range(3):
        a, b, c = offset[:, axis], offset[:, (axis+1) % 3], offset[:, (axis+2) % 3]
        for sign in (-1, 1):
            face = a + sign * dx/2
            active = np.abs(face) > dx*1e-8
            fixed = face[active]
            for u, wu in zip(nodes, weights):
                first = (b[active] + dx*u/2)**2 + fixed**2
                for v, wv in zip(nodes, weights):
                    total[active] += sign*fixed*wu*wv/(first + (c[active] + dx*v/2)**2)
    return total*(dx/2)**2


def cell_kernel_gaussian(offset, dx, order):
    # Independent volume integral using 1/r^2 = integral exp(-s*r^2) ds.
    # Set sqrt(s) = t/((1-t)*dx). Integrate each Gaussian over its
    # rectangular cell analytically before the remaining scalar quadrature.
    nodes, weights = leggauss(order)
    result = np.zeros(len(offset))
    low, high = (offset-dx/2)/dx, (offset+dx/2)/dx
    low = np.rint(low)
    high = np.rint(high)
    for t, wt in zip((nodes+1)/2, weights/2):
        ratio = t/(1-t)
        lengths = np.sqrt(np.pi)/(2*ratio) * (erf(high*ratio)-erf(low*ratio))
        result += wt * dx * 2*t/(1-t)**3 * np.prod(lengths, axis=1)
    return result


def contrasts(values):
    # At each of x=-1/2 and x=+1/2, compare y=0 with the mean at y=±1/2.
    # z stays zero. Return the x-even and x-odd parts of the two contrasts.
    left, right = values[:, 1] - (values[:, 0] + values[:, 2])/2
    return np.array([(left+right)/2, (right-left)/2])


def measure_level(n, base_dir, driven_dir):
    scratch = Path("_scratch/s11c")
    records = Path("docs/_measurements")
    coords = (np.arange(n)+0.5-n/2)*4/n
    positions = np.stack(np.meshgrid(coords, coords, coords, indexing="ij"), axis=-1)
    mask = np.sum(positions**2, axis=-1) < 1
    points = np.array([[x, y, 0.] for x in (-0.5, 0.5) for y in (-0.5, 0., 0.5)])
    dx = 4/n
    kernels = []
    stencils = []
    for point in points:
        offset = positions[mask] - point
        np.testing.assert_allclose(offset/dx-0.5, np.rint(offset/dx-0.5), atol=2e-14, rtol=0)
        surface = cell_kernel_surface(offset, dx, 8)
        volume = cell_kernel_gaussian(offset, dx, 64)
        kernels.append((surface, volume))
        stencil = np.all(np.abs(positions-point) < 0.51*dx, axis=-1)
        if stencil.sum() != 8 or not mask[stencil].all():
            raise ValueError("Interpolation must stay inside the opening")
        stencils.append(stencil)
    kernels = np.array(kernels)
    print(f"grid={n} max cell integral difference={abs(kernels[:,0]-kernels[:,1]).max():.3g}")
    for magnitude in (0., 0.02, 0.01):
        directory = scratch / (driven_dir if magnitude else base_dir)
        outcome = json.loads((directory/"outcome.json").read_text())
        if not outcome["limitsVerified"] or outcome["exitCode"] or outcome["childOutcome"]["exitCode"]:
            raise ValueError(f"Invalid original run: {directory}")
        manifest = "mouth_intake_redistribution.json" if magnitude else "mouth_feed_weak_run.json"
        field, current, field_dx, _ = read_state(directory/"fields", n, "restoring_slope",
                                                  magnitude, manifest_fields(records/manifest))
        if field_dx != dx or np.any(current[~mask]) or not np.all(current[mask] > 0):
            raise ValueError("Incorrect geometry or inward current")
        face = field[-1]
        for sign in ((1, -1) if magnitude else (1,)):
            speed = sign*magnitude
            envelope = face if sign == 1 else face[::-1]
            j = current if sign == 1 else current[::-1]
            rate = dx**3 * j.sum()
            phase_samples = np.array([envelope[s].mean() for s in stencils])
            if np.min(phase_samples.real) <= 0:
                raise ValueError("Expected the archived unwrapped phase branch")
            feed = contrasts((np.angle(phase_samples) + speed*points[:, 0]).reshape(2, 3))
            # Neumann Green function on the plane for flow into the half-space.
            potentials = -(kernels @ j[mask])/(2*np.pi**2)
            recv = contrasts(potentials[:,0].reshape(2, 3))
            check = contrasts(potentials[:,1].reshape(2, 3))
            print(f"n={n} U={speed:+.2f} J={rate:.12g} "
                  f"feed_even={feed[0]:+.12g} bulk_even={recv[0]:+.12g} "
                  f"feed_odd={feed[1]:+.12g} bulk_odd={recv[1]:+.12g} "
                  f"quadrature_difference={np.max(abs(check-recv)):.3g}", flush=True)
            np.testing.assert_allclose(rate, 0.25, rtol=0, atol=1e-12)
            np.testing.assert_allclose(potentials[:,0], potentials[:,1], rtol=1e-7, atol=1e-12)


if __name__ == "__main__":
    measure_level(24, "mouth_feed_weak01", "mouth_redistribution_coarse01")
    measure_level(32, "mouth_feed_depth01", "mouth_redistribution_depth01")
