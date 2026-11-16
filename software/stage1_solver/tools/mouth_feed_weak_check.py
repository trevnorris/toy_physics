"""Check the saved brane-side solutions without importing the solver.

Use the mixed Madelung stress, a one-sided polynomial at the end of the
feeding slab, and virtual integration by parts on the full aperture. Under
the restoring law the normal log-density derivative is compactly supported
inside the opening and vanishes at the rim. There is no rim term.
"""

import argparse
from pathlib import Path

import numpy as np
from scipy.interpolate import RegularGridInterpolator


def dx4(field, dx):
    return (np.roll(field, 2, axis=0) - 8*np.roll(field, 1, axis=0)
            + 8*np.roll(field, -1, axis=0) - np.roll(field, -2, axis=0))/(12*dx)


def check_file(path, law):
    with np.load(path) as packet:
        f, current = packet["field"], packet["current"]
        dx, dz, speed = (float(packet[name]) for name in ("dx", "dz", "speed"))
        mouth = packet["mouth"]
    count = f.shape[1]
    coord = (np.arange(count) + 0.5 - count/2)*dx
    x, y, z = np.meshgrid(coord, coord, coord, sparse=True, indexing="ij")
    radius_squared = x**2 + y**2 + z**2
    aperture = radius_squared < 1
    weight = np.maximum(1-radius_squared, 0)**3
    if not np.array_equal(aperture, mouth) or abs(current.sum()*dx**3 - 0.25) > 1e-12:
        raise ArithmeticError("Incorrect saved aperture or intake")
    if abs(f[0]-1).max() > 1e-12:
        raise ArithmeticError("Feed envelope is not the calibrated plane wave")
    # Verify the stationary equation with real-space neighbours, including
    # the last plane where a ghost enforces the normal density and flux laws.
    amplitude = abs(f[-1])
    restoring = weight if law == "restoring_slope" else np.zeros_like(weight)
    slope = -2*restoring*(1-1/amplitude)
    wall_derivative = (slope/2 + 1j*current/amplitude**2)*f[-1]
    ghost = f[-2] + 2*dz*wall_derivative
    w_with_ghost = np.concatenate((f, ghost[None, ...]))
    lap = np.diff(w_with_ghost, n=2, axis=0)/dz**2
    bulk = f[1:]
    for axis in (1, 2, 3):
        lap += (np.roll(bulk, 1, axis) - 2*bulk + np.roll(bulk, -1, axis))/dx**2
    adv = (np.roll(bulk, -1, axis=1) - np.roll(bulk, 1, axis=1))/(2*dx)
    residual = -lap/2 - 1j*speed*adv + (abs(bulk)**8-1)*bulk/4
    max_residual = np.max(abs(residual))
    if max_residual > 2e-8:
        raise ArithmeticError(f"Large steady-equation residual {max_residual}")
    # Independently read the flux on every feed-depth link and the imposed
    # lateral plane wave from its phase change between two feed samples.
    link_rates = np.imag(f[:-1].conj()*f[1:]).sum(axis=(1, 2, 3))*dx**3/dz
    feed_phase = f[0,:,0,0]*np.exp(1j*speed*coord)
    incident = np.angle(np.vdot(feed_phase[:-1], feed_phase[1:]))/dx
    if np.max(abs(link_rates-0.25)) > 2e-7 or abs(incident-speed) > 1e-12:
        raise ArithmeticError("Current or feed drive disagrees with input")
    n = amplitude**2
    face = f[-1]
    vx = speed + np.imag(face.conj()*dx4(face, dx))/n
    carried = dx**3 * np.sum(current*vx)
    stress_parts = dx**3 * np.sum(slope*dx4(n, dx))/4
    # Read the same physical integral from the full three-area ball with a
    # one-sided four-point field derivative. This is a separate diagnostic
    # with its own spatial and rim-integration errors.
    logn = 2*np.log(abs(f[-4:]))
    slope_fd = (11*logn[-1]-18*logn[-2]+9*logn[-3]-2*logn[-4])/(6*dz)
    stress_fd = -n*dx4(slope_fd, dx)/4
    j_fd = np.imag(face.conj()*(11*face-18*f[-2]+9*f[-3]-2*f[-4]))/(6*dz)
    # Project the solved face density onto a dipole. Give the axial
    # half-difference at selected radii as a simple centre-to-rim readout.
    interp = RegularGridInterpolator((coord,)*3, n, bounds_error=True)
    radii = np.array([0, 0.25, 0.5, 0.75, 1])
    points = np.column_stack((radii, np.zeros((len(radii), 2))))
    dipole = (interp(points)-interp(-points))/2
    return (max_residual, np.max(abs(link_rates-0.25)), incident, carried, stress_parts,
            dx**3 * np.sum(j_fd[mouth]*vx[mouth]), dx**3*stress_fd[mouth].sum(), dipole)


def main(directory, levels):
    for n in levels:
        for law in ("flat_density_slope", "restoring_slope"):
            rows = [check_file(directory/f"n{n}_{law}_u{u:g}.npz", law)
                    for u in (0., 0.01, 0.02)]
            for u, row in zip((0.01, 0.02), rows[1:]):
                vals = [(row[i]-rows[0][i])/u for i in range(3, 7)]
                profile = (row[-1]-rows[0][-1])/u
                print(f"n={n} law={law} U={row[2]:g} max_residual={row[0]:.3g} "
                      f"max_flux_error={row[1]:.3g} "
                      "slopes(flow,stress_parts,flow_fd,stress_fd)="
                      + ",".join(f"{v:.10g}" for v in vals)
                      + " density_half_difference/U(r=0,.25,.5,.75,1)="
                      + ",".join(f"{v:.8g}" for v in profile), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    parser.add_argument("--levels", type=int, nargs="+", default=[16, 24, 32])
    args = parser.parse_args()
    main(args.directory, args.levels)
