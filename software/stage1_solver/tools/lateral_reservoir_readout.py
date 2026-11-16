"""Read signed one-drain currents and paired-minus-alone mouth loads from saved fields."""

import argparse
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.interpolate import RegularGridInterpolator

from mouth_feed_weak_check import dx4


def depth_sum(field, dz):
    return dz*(field[1:-1].sum(axis=0)+(field[0]+field[-1])/2)


def sphere_currents(f, x, r, dz):
    dx = x[1]-x[0]
    # Currents on cell links conserve number in the steady difference equation.
    jx = depth_sum((f[:, :-1].conj()*f[:, 1:]).imag/dx, dz)
    jr = depth_sum((f[:, :, :-1].conj()*f[:, :, 1:]).imag/dx, dz)
    # Also sample a geometrically spherical surface, interpolating staggered
    # components separately. Reflection at the radial axis enforces regularity.
    axial = RegularGridInterpolator(((x[1:]+x[:-1])/2, np.r_[-r[0],r]),
                                    np.column_stack((jx[:, 0], jx)))
    radial = RegularGridInterpolator((x, np.r_[0,(r[1:]+r[:-1])/2]),
                                     np.column_stack((np.zeros(len(x)),jr)))
    for radius in (2, 3, 4, 6, 8):
        inside = (x[:, None]**2+r**2 < radius**2).astype(int)
        # Face area on each axial annulus or radial cylindrical strip.
        out_x = np.sum((inside[:-1]-inside[1:])*jx*(2*np.pi*r*dx))
        out_r = np.sum((inside[:, :-1]-inside[:, 1:])*jr*(2*np.pi*(r[:-1]+dx/2)*dx))
        cut_rate = -(out_x+out_r)
        readings = []
        for order in (96, 192):
            mu, w = leggauss(order)
            sine = np.sqrt(1-mu**2)
            points = np.column_stack((radius*mu, radius*sine))
            readings.append(-2*np.pi*radius**2*np.dot(w, axial(points)*mu+radial(points)*sine))
        if (not np.isfinite([cut_rate,*readings]).all() or abs(cut_rate-0.25) > 1e-8):
            raise ArithmeticError(f"Isolated drain failed the shell balance at {radius=}")
        print(f"  R={radius} link_cut={cut_rate:.12g} sphere96={readings[0]:.12g} "
              f"sphere192={readings[1]:.12g} mean_inward_j={readings[1]/(4*np.pi*radius**2):.12g}")


def mouth_load(layers, currents, x, r, dz, expected_currents):
    dx = x[1]-x[0]
    dv = 2*np.pi*r*dx**2
    mask = x[:, None]**2+r**2 < 1
    j = currents[-1]
    if (len(currents) != expected_currents or np.any(j[~mask] != 0)
            or abs(np.sum(dv*j)-0.25) > 2e-13):
        raise ValueError("Target rate, position or mouth count changed")
    n = abs(layers[-1])**2
    vx = (layers[-1].conj()*dx4(layers[-1],dx)).imag/n
    weight = np.maximum(1-x[:, None]**2-r**2, 0)**3
    slope = -2*weight*(1-1/np.sqrt(n))
    def normal(values):
        return (11*values[-1]-18*values[-2]+9*values[-3]-2*values[-4])/(6*dz)
    j_fd = (layers[-1].conj()*normal(layers)).imag
    slope_fd = normal(2*np.log(abs(layers)))
    # The source is on the negative x axis. Use the same signed reference
    # on the isolated ball, at the same location in the same reservoir.
    def inward(values):
        return -np.sum((dv*values)[mask])
    return np.array([inward(j*vx), inward(slope*dx4(n,dx)/4),
                     inward(j_fd*vx), inward(-n*dx4(slope_fd,dx)/4)])


def run(directory):
    with np.load(directory/"d0.npz", allow_pickle=False) as packet:
        x, r, dz = packet["x"], packet["r"], float(packet["dz"])
        f = packet["field"]
        print(f"extent={r[-1]+r[0]:.12g} dx={2*r[0]:.12g} dz={dz:.12g}")
        sphere_currents(f, x, r, dz)
        baseline = mouth_load(f[-4:], packet["currents"], x, r, dz, 1)
        base_speed = packet["currents"][0]/abs(f[-1])**2
    if not np.isfinite(baseline).all() or np.max(abs(baseline)) > 1e-11:
        raise ArithmeticError(f"Unexpected isolated sideways load {baseline}")
    print(f"  alone_load={baseline[:2].sum():.12g}")
    forces = {}
    for d in (3,4,6,8):
        with np.load(directory/f"d{d}.npz", allow_pickle=False) as p:
            if (int(p["d"]) != d or float(p["dz"]) != dz or not np.array_equal(p["x"],x)
                    or not np.array_equal(p["r"],r)):
                raise ValueError("Pair is not on the isolated reference grid")
            layers, currents = p["layers"], p["currents"]
        for current, center in zip(currents, (-d,0)):
            ref = np.roll(base_speed, round(center/(x[1]-x[0])), axis=0)
            if (np.any((ref>0) != (current>0)) or abs(np.sum(2*np.pi*r*(2*r[0])**2*current)-0.25)>2e-13
                    or np.ptp(current[ref>0]/(abs(layers[-1][ref>0])**2*ref[ref>0]))>2e-11):
                raise ValueError("Changed normal-speed profile or rate")
        parts = mouth_load(layers, currents, x, r, dz, 2)-baseline
        if not np.isfinite(parts).all():
            raise ArithmeticError("Nonfinite force")
        total = parts[:2].sum()
        forces[d] = total
        print(f"  d={d} carried={parts[0]:.12g} stress={parts[1]:.12g} total={total:.12g} "
              f"carried_fd={parts[2]:.12g} stress_fd={parts[3]:.12g} "
              f"total_fd={parts[2:].sum():.12g} d2_total={d*d*total:.12g}")
    print(f"  F6_over_F3={forces[6]/forces[3]:.12g} F8_over_F4={forces[8]/forces[4]:.12g}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directories", nargs="+", type=Path)
    for directory in parser.parse_args().directories:
        run(directory)
