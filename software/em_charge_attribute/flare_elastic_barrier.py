"""Elastic work and standing-wave energy on the previously tested tilt path."""
import json
from pathlib import Path
import re

import numpy as np
from numpy.polynomial.legendre import leggauss

from wave_at_level_patch import Collar, transition


def wall(tilt, nodes):
    # All lengths in units R. Return area and integral of mean curvature
    # squared separately on the fixed-neck join and the conversion region.
    q, weights = leggauss(nodes)
    mu = q[None, :]
    sine = np.sqrt(1-mu**2)
    pieces = []
    for start in (1, 2):
        r = start + (q[:, None]+1)/2
        x = np.clip(r-1, 0, 1)
        f = transition(x)
        df = 30*x**2*(1-x)**2
        ddf = 60*x*(1-x)*(1-2*x)
        # Gradient and Hessian of h=epsilon*(r^2/2-tilt*z*f)
        # in an orthonormal radial, polar, azimuthal frame.
        pr = 0.2*(r-tilt*mu*(f+r*df))
        pt = 0.2*tilt*f*sine
        rr = 0.2*(1-tilt*mu*(2*df+r*ddf))
        rt = 0.2*tilt*df*sine
        tt = 0.2*(1-tilt*mu*df)
        w = np.sqrt(1+pr**2+pt**2)
        projected = rr*pr**2 + 2*rt*pr*pt + tt*pt**2
        mean = (rr+2*tt-projected/w**2)/(3*w)
        measure = np.pi*r**2*w*weights[:, None]*weights[None, :]
        pieces.append([measure.sum(), (measure*mean**2).sum()])
    return np.array(pieces)


def archived_frequencies(cells):
    manifest = json.loads(Path("docs/_measurements/wave_at_level_patch.json").read_text())
    values = {}
    for line in manifest["runs"]["sympy"]["stdout_claim_lines"]:
        match = re.match(r"cells=(\d+) angular=\d+ tilt=([\d.]+) omega=([\d.]+)", line)
        if match and int(match[1]) == cells:
            values[float(match[2])] = float(match[3])
    if not {0, 2, 2.49, 2.5, 2.51, 3, 4} <= values.keys():
        raise ValueError("Incomplete reference spectrum")
    return values


def main():
    # Reuse completed spectral readings. Fill the gaps up to the first complete loss.
    cache = archived_frequencies(24)
    solver = Collar(24, 10)
    grid = sorted(set(np.arange(61)/20) | {0.025, 0.0125, 2.49, 2.51, 4.0})
    baseline = wall(0, 64)
    print("baseline join/convert [A,B] =", baseline.tolist())
    rows = []
    for beta in grid:
        parts = wall(beta, 64)
        coarse = wall(beta, 40)
        if np.max(abs(parts-coarse)) > 5e-8:
            raise ArithmeticError("Elastic quadrature unresolved")
        if beta not in cache:
            cache[beta] = solver.solve(beta)[0]
        da, db = (parts-baseline).sum(axis=0)
        dw = cache[beta]-cache[0]
        print(f"beta={beta:g} delta_area={da:.10f} delta_bend={db:.10f} "
              f"delta_frequency={dw:+.10f}", flush=True)
        if beta in (2, 2.5, 3, 4):
            print("  join/convert differences =", (parts-baseline).tolist())
        if beta > 0 and (dw >= 0 or da <= 0):
            raise ArithmeticError("Area or frequency change requires investigation")
        rows.append([beta, da, db, dw])
    rows = np.array(rows)
    for beta in (2, 2.5, 3):
        _, da, db, dw = rows[rows[:, 0] == beta][0]
        print(f"beta={beta} equal-energy action: {da/-dw:.8f} + {db/-dw:.8f} k")
    # A check free of an arbitrarily selected J: allow only a common scale R
    # of this isolated energy, hold shape, action and dimensional moduli fixed.
    a0, b0 = baseline.sum(axis=0)
    j_sigma, j_bend = 3*a0/cache[0], b0/cache[0]
    print(f"isolated scale balance j={j_sigma:.9f}+{j_bend:.9f} k")
    active = rows[rows[:, 0] <= 3]
    b, a, h, dw = active.T
    offset, slope = a+j_sigma*dw, h+j_bend*dw
    for beta in (0.0125,0.025,0.05,2,2.5,3):
        idx = np.where(b == beta)[0][0]
        print(f"scale balanced beta={beta} offset={offset[idx]:+.10f} "
              f"slope={slope[idx]:+.10f} k_equal={-offset[idx]/slope[idx]:.8f}")
    for k in (0.0, 1.0, 4.0, 10.0):
        energies = offset + k*slope
        peak = energies.argmax()
        print(f"scale balanced k={k:g} sampled maximum={energies[peak]:+.9f} "
              f"at beta={b[peak]:g} energy_at_loss={energies[-1]:+.9f}")
    if np.all(slope[1:] > 0):
        ratios = -offset[1:]/slope[1:]
        print(f"first positive sampled energy above base at k={min(ratios):.8f} "
              f"beta={b[1:][ratios.argmin()]:g}")


if __name__ == "__main__":
    main()
