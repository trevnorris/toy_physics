"""Linear carrier reaction in an accelerating sleeve with comoving baths.

The feed is fixed in the accelerating frame; the far bath can impose a known
lateral phase gradient. The frame potential is +m*A*x. The outputs below
are derivatives at A=0; acceleration is never inferred from a force.
"""

import argparse
from pathlib import Path

import numpy as np
from scipy.sparse import bmat, diags
from scipy.sparse.linalg import spsolve

from stage1_solver.config import TensorGridSpec
from stage1_solver.open_sleeve import OpenSleeve


def carrier_dipole(sleeve, psi, far_phase_slope):
    r = sleeve.r_index
    x, y = psi.real, psi.imag
    n = abs(psi)**2
    volume = sleeve.volume
    h = 5 * sleeve.eos_k / 4
    # The radial half-cell face conductances are unchanged for an l=1
    # perturbation. The angular kinetic term integrates 2/r^2 over a shell.
    angular = sleeve.kinetic * 8 * np.pi * sleeve.dr * sleeve.dw
    base = sleeve.operator + diags(angular + volume*h*(n**4 - sleeve.bath_density**4))
    slope = 8 * volume * h * n**3
    mixed = diags(slope*x*y)
    jac = bmat([[base + diags(slope*x*x), mixed],
                [mixed, base + diags(slope*y*y)]], format="csc")
    radial_center = (r + 0.5)*sleeve.dr
    rhs = -sleeve.constituent_mass * radial_center * volume * psi
    # Optional measured/controlled relative outlet velocity per acceleration
    # (time units), imposed as a phase gradient at the far receiver only.
    rhs += (1j * sleeve.kinetic * sleeve.outlet * np.sqrt(sleeve.bath_density)
            * radial_center * far_phase_slope)
    vector = np.r_[rhs.real, rhs.imag]
    solved = spsolve(jac, vector)
    residual = np.max(abs(jac @ solved - vector))
    real, imag = np.split(solved, 2)
    return real + 1j*imag, residual


def cut_currents(logn, vel_potential, b, g, *, z, spacing, dr, m, hbar):
    """Flux along +w across an interior flat cut at z (a face index).

    The cells on both sides of this cut lie in the straight sleeve. H and Phi
    are log-density and velocity potential; b and g are their cos(theta)
    dipole perturbations per unit acceleration in the comoving frame.
    """
    def face(arr):
        return (arr[:, z-1] + arr[:, z]) / 2

    def dw(arr):
        return (arr[:, z] - arr[:, z-1]) / spacing

    def dr_neumann(arr):
        padded = np.r_[arr[0], arr, arr[-1]]
        return (padded[2:]-padded[:-2]) / (2*dr)

    def dr_dipole(arr):
        # Odd regularity through the origin, even across the reflecting wall.
        padded = np.r_[-arr[0], arr, arr[-1]]
        return (padded[2:]-padded[:-2]) / (2*dr)

    nr = logn.shape[0]
    r = (np.arange(nr) + 0.5)*dr
    shells = 4*np.pi/3 * np.diff((dr*np.arange(nr+1))**3)
    n = np.exp(face(logn))
    vw = dw(vel_potential)
    vr = dr_neumann(face(vel_potential))
    vel1x = (dr_dipole(face(g)) + 2*face(g)/r) / 3
    advective = m*n*(vw*vel1x + vr*(dw(g) + face(b)*vw)/3)
    quantum = -hbar**2*n/(12*m)*(face(b)*dr_neumann(dw(logn))
                                  + dr_dipole(dw(b)) + 2*dw(b)/r)
    return np.dot(shells, advective), np.dot(shells, quantum)


def side_stress(logn, b, first, last, *, dr, dw, m, hbar, eos_k):
    """Outward flux to the straight reflecting wall on [first,last) cells.

    At the wall H_r=B_r=0. The angular average of the x-normal quantum
    term retains both the radial Hessian and the spherical cross derivative.
    """
    def boundary_value(arr):
        return (9*arr[-1] - arr[-2]) / 8

    def second_radial(arr):
        return (arr[-2] - arr[-1]) / dr**2

    radius = logn.shape[0] * dr
    h_wall, b_wall = boundary_value(logn), boundary_value(b)
    density = np.exp(h_wall)
    pressure = 5 * eos_k * density**5 * b_wall / 3
    quantum = -hbar**2*density/(12*m) * (
        b_wall * second_radial(logn) + second_radial(b) - 2*b_wall/radius**2
    )
    area = 4*np.pi * radius**2 * dw
    return area*pressure[first:last].sum(), area*quantum[first:last].sum()


def bath_quantum(psi, eta, bath, *, dr, dw, hbar, m, at_exit):
    """Differentiate mixed density stress at a flat, constant-density bath.

    For H=ln n and its dipole B, H and B at the bath are spatial constants
    (B=0). Isotropic pressure has zero xw component. The angular/radial
    integral of mixed density-gradient stress is proportional to [r^2 B_w] at
    the rim. Estimate the interior corner trace from inside the carrier.
    """
    b = 2*(eta/psi).real
    # Homogeneous density dipole trace. Sign is +w at either end.
    if at_exit:
        bw = -(9*b[:, -1] - b[:, -2])/(3*dw)
    else:
        bw = (9*b[:, 0] - b[:, 1])/(3*dw)
    rim_bw = (9*bw[-1] - bw[-2]) / 8
    radius = psi.shape[0]*dr
    return -hbar**2*bath/(12*m) * 4*np.pi*radius**2 * rim_bw


def examine(path, outlet_velocity_per_acceleration):
    with np.load(path, allow_pickle=False) as d:
        nr, nw = int(max(d["r_index"])+1), int(max(d["w_index"])+1)
        m, hbar, k, bath = (float(d[key]) for key in ("constituent_mass", "hbar", "eos_k", "bath_density"))
        dr, dw = float(d["dr"]), float(d["dw"])
        a, length = int(d["mouth_cells"]), int(d["sleeve_cells"])
        if length % 4 or length < 8:
            raise ValueError("The quarter-length readout requires a resolved, face-aligned window")
        domain = OpenSleeve(TensorGridSpec(nr=nr, nw=nw, r_max=nr*dr, w_min=0, w_max=nw*dw),
                            a, length, hbar, m, k, bath)
        if not (np.array_equal(domain.r_index, d["r_index"]) and
                np.array_equal(domain.w_index, d["w_index"])):
            raise ValueError("Incompatible saved geometry or field ordering")
        psi = d["psi"]
        rate = float(d["rate"])
        feed_phase = float(d["feed_phase"])
    entry, exit_rate, far, _ = domain.budget(psi, feed_phase)
    if np.max(abs(np.array([entry, exit_rate, far]) - rate)) > 1e-10 or min(abs(psi)) <= 0:
        raise ArithmeticError("Background carrier current or density inconsistent")
    eta, residual = carrier_dipole(domain, psi, m*outlet_velocity_per_acceleration/hbar)
    # All collar derivatives use only occupied straight-sleeve cells.
    mask = domain.w_index < length
    shape = (a, length)
    n = (abs(psi[mask])**2).reshape(shape)
    h = np.log(n)
    phi = hbar/m * np.unwrap(np.unwrap(np.angle(psi[mask].reshape(shape)), axis=0), axis=1)
    relative = (eta[mask]/psi[mask]).reshape(shape)
    b = 2*relative.real
    g = hbar/m * relative.imag
    first, last = length//4, 3*length//4
    ka = dict(spacing=dw, dr=dr, m=m, hbar=hbar)
    incoming = cut_currents(h, phi, b, g, z=first, **ka)
    outgoing = cut_currents(h, phi, b, g, z=last, **ka)
    wall = side_stress(h, b, first, last, dr=dr, dw=dw, m=m, hbar=hbar, eos_k=k)
    receiver = ~mask
    receiver_shape = (nr, nw-length)
    bath_args = dict(dr=dr, dw=dw, hbar=hbar, m=m)
    feed_q = bath_quantum(psi[mask].reshape(shape), eta[mask].reshape(shape), bath,
                          at_exit=False, **bath_args)
    exit_q = bath_quantum(psi[receiver].reshape(receiver_shape), eta[receiver].reshape(receiver_shape),
                          bath, at_exit=True, **bath_args)
    far_carried = m * far * outlet_velocity_per_acceleration
    mass = m * np.sum(domain.volume[mask].reshape(shape)[:, first:last] * n[:, first:last])
    # d/dA of integral div Pi = -m integral n, since the fictitious body
    # force is -m*n*A e_x. Hydrostatic isotropic pressure has zero xw flux.
    balance = sum(outgoing) - sum(incoming) + sum(wall) + mass
    row = np.r_[incoming, outgoing, wall, mass, feed_q, exit_q, far_carried, balance]
    print(f"rate={rate:g} level={a} rate_in/out/far={entry:.11g},{exit_rate:.11g},{far:.11g} "
          f"max_linear_equation_residual={residual:.3e}")
    print(f"  sleeve_quarter: carried={incoming[0]:+.12g} quantum={incoming[1]:+.12g} "
          f"sleeve_three_quarters: carried={outgoing[0]:+.12g} quantum={outgoing[1]:+.12g}")
    print(f"  feed_quantum={feed_q:+.12g} far_receiver_quantum={exit_q:+.12g} "
          f"far_receiver_carried={far_carried:+.12g}")
    print(f"  wall_pressure={wall[0]:+.12g} wall_quantum={wall[1]:+.12g} "
          f"wall_total={sum(wall):+.12g} carrier_inventory_mass={mass:.12g} "
          f"balance_residual={balance:+.12g} relative_residual={abs(balance)/mass:.6g}", flush=True)
    if not np.all(np.isfinite(np.r_[row, residual])) or residual > 1e-11 or mass <= 0:
        raise ArithmeticError("Invalid field solution or lateral readings")
    return rate, a, row


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("fields", nargs="+", type=Path)
    parser.add_argument("--outlet-velocity-per-acceleration", type=float, default=0.0,
                        help="controlled sideways bath speed per imposed acceleration, in time units")
    args = parser.parse_args()
    if not np.isfinite(args.outlet_velocity_per_acceleration):
        parser.error("Outlet velocity per acceleration must be finite")
    print(f"outlet_velocity_per_acceleration={args.outlet_velocity_per_acceleration:g}")
    cases = {}
    for path in args.fields:
        rate, level, values = examine(path, args.outlet_velocity_per_acceleration)
        if (rate, level) in cases:
            raise ValueError("Duplicate field resolution and rate")
        cases[rate, level] = values
    for rate, level in sorted(cases):
        if rate and (0, level) in cases:
            diff = cases[rate, level] - cases[0, level]
            print(f"drain_minus_zero level={level} sleeve_three_quarters_carried={diff[2]:+.9g} "
                  f"sleeve_three_quarters_quantum={diff[3]:+.9g} wall_total={sum(diff[4:6]):+.9g} "
                  f"far_receiver_quantum={diff[8]:+.9g} far_receiver_carried={diff[9]:+.9g} "
                  f"balance_residual={diff[-1]:+.9g}", flush=True)
