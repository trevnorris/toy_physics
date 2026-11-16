"""Integrate companion-directed advection and mixed stress on each 3D mouth.

Restore the full periodic x axis from the reflected solution. Compare each
paired mouth to the isolated field translated to the same aperture cells.
The two other reflections give a factor of four in the integration measure.
"""

import argparse
from pathlib import Path

import numpy as np

from mouth_feed_weak_check import dx4


def read_faces(path):
    with np.load(path, allow_pickle=False) as p:
        layers = p["layers"]
        face = np.concatenate((layers[:, ::-1], layers), axis=1)
        j = np.concatenate((p["current"][::-1], p["current"]))
        weight = np.concatenate((p["weight"][::-1], p["weight"]))
        dx, dz = float(p["dx"]), float(p["dz"])
        d, width, depth = (int(p[key]) for key in ("separation", "width", "depth"))
    size = j.shape[0]//2
    x = (np.arange(2*size)-size+0.5)*dx
    y = (np.arange(size)+0.5)*dx
    r2 = ((abs(x)[:, None, None]-d/2)**2 + y[None, :, None]**2 + y[None, None, :]**2)
    if not np.allclose(weight, np.maximum(1-r2, 0)**3, atol=2e-14, rtol=0):
        raise ValueError("Archived rim or profile differs from the intended geometry")
    n = abs(face[-1])**2
    vx = (face[-1].conj()*dx4(face[-1], dx)).imag/n
    slope = -2*weight*(1-1/np.sqrt(n))
    def normal(data):
        return (11*data[-1]-18*data[-2]+9*data[-3]-2*data[-4])/(6*dz)
    j_fd = (face[-1].conj()*normal(face)).imag
    slope_fd = normal(2*np.log(abs(face)))
    # The slope imposed at each mouth vanishes at the rim. Integration by
    # parts on that ball gives the slope * grad(n) expression. The direct
    # readout differentiates the reconstructed slope and cuts at the rim.
    currents = (j*vx, slope*dx4(n, dx)/4, j_fd*vx, -n*dx4(slope_fd, dx)/4)
    return (width, depth, dx, dz), x, r2 < 1, j, n, currents


def run(directory):
    grid, x, mask0, j0, n0, alone = read_faces(directory/"d0.npz")
    width, depth, dx, dz = grid
    element = 4*dx**3
    baseline = np.array([element*f[mask0].sum() for f in alone])
    base_rate = element*j0[mask0].sum()
    if np.max(abs(baseline)) > 1e-12 or abs(base_rate-0.25) > 2e-13:
        raise ArithmeticError("Isolated mouth has a lateral load or incorrect intake")
    print(f"width={width} depth={depth} dx={dx:.12g} dz={dz:.12g} "
          f"alone_load_x={(baseline[0]+baseline[1]):.12g} alone_rate={base_rate:.12g}")
    totals = {}
    for d in (3, 4, 6):
        grid_d, _, mask, j, n, pair = read_faces(directory/f"d{d}.npz")
        if grid_d != grid:
            raise ValueError("Pair and isolated grids differ")
        sides = []
        for side in (-1, 1):
            cut = mask & (side*x[:, None, None] > 0)
            shift = side*round(d/(2*dx))
            if not np.array_equal(cut, np.roll(mask0, shift, axis=0)):
                raise ValueError("Pair and isolated openings do not coincide after translation")
            # Recover the normal speed profile up to its one allowed scalar.
            speed0 = np.roll(j0/n0, shift, axis=0)[cut]
            ratio = (j[cut]/n[cut])/speed0
            if (abs(element*j[cut].sum()-0.25) > 2e-13 or np.ptp(ratio) > 2e-12
                    or not np.isfinite(ratio).all()):
                raise ArithmeticError("Per-mouth rate or normal-speed shape changed")
            sides.append(-side*element*np.array([
                np.sum(f[cut] - np.roll(f0, shift, axis=0)[cut]) for f, f0 in zip(pair, alone)]))
        if not np.isfinite(sides).all() or np.max(abs(sides[1]-sides[0])) > 1e-12:
            raise ArithmeticError("Mismatch under reflection across the pair midpoint")
        carried, stress, carried_fd, stress_fd = sides[1]
        total = carried + stress
        totals[d] = total
        print(f"  d={d} carried={carried:.12g} stress={stress:.12g} total={total:.12g} "
              f"carried_fd={carried_fd:.12g} stress_fd={stress_fd:.12g} "
              f"total_fd={carried_fd+stress_fd:.12g} d2_total={d*d*total:.12g}")
    if totals[3] != 0:
        print(f"  F6_over_F3={totals[6]/totals[3]:.12g}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directories", nargs="+", type=Path)
    for directory in parser.parse_args().directories:
        run(directory)
