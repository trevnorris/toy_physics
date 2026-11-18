"""Carrier normal momentum on the two reflecting faces outside one active drain."""

import json
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.interpolate import RegularGridInterpolator
import sympy as sp

from lateral_source_rate import verified_packet


def symbolic_checks():
    b, w0 = sp.symbols("B W0", positive=True)
    c, kw = sp.symbols("C kW", real=True)
    z, w, density, t_plus, t_minus = sp.symbols("Z W theta t_plus t_minus", real=True)
    # t_plus and t_minus are both signed along laboratory +w.
    work = t_plus*(z+w/2) + t_minus*(z-w/2)
    print("center and width loads =", sp.diff(work, z), sp.diff(work, w))
    energy = b*density**2/2 + c*density*w + kw*w**2/2
    relaxed = sp.solve(sp.diff(energy, density), density)[0]
    kr = sp.diff(energy.subs(density, relaxed), w, 2)
    kc = sp.diff(energy.subs(density, -w/w0), w, 2)
    print("closed-column minus relaxed stiffness =", sp.factor(kc - kr))
    assert sp.simplify(kc - kr - b*(1/w0-c/b)**2) == 0
    r, amp, stiffness, grad = sp.symbols("r A K D", positive=True)
    def lap(f):
        return sp.diff(r**2*sp.diff(f, r), r)/r**2
    # A negative r^-4 tail in the actual net applied load is an assumption.
    gapped = -amp/(stiffness*r**4)
    soft = amp/(2*grad*r**2)
    print("width operator on driven tails =", sp.expand(stiffness*gapped-grad*lap(gapped)),
          sp.simplify(-grad*lap(soft)))
    assert sp.simplify(-grad*lap(soft) + amp/r**4) == 0


def normal_stress(face, next_layer, second_layer, radial, dx, dz):
    """Subtract the uniform bath stress. Use only where the face reflects."""
    n = abs(face)**2
    pressure = np.expm1(5*np.log(n))/5
    d1 = abs(next_layer)**2 - n
    d2 = abs(second_layer)**2 - n
    nww = (16*d1 - d2)/(6*dz**2)
    direct = pressure - nww/4

    # Independently recover f_ww from the steady carrier equation using
    # only in-plane neighbors of this face. Outside the aperture f_w = 0.
    lap_x = (face[2:] - 2*face[1:-1] + face[:-2])/dx**2
    radial_flux = (radial[:-1]+dx/2)*np.diff(face[1:-1], axis=1)/dx
    lap_r = np.diff(radial_flux, axis=1)/(dx*radial[1:-1])
    axis_value = radial_flux[:, 0]/(dx*radial[0])
    lap = lap_x[:, :-1] + np.column_stack((axis_value, lap_r))
    core = face[1:-1, :-1]
    potential = np.expm1(4*np.log(n[1:-1, :-1]))/4
    fww = 2*potential*core - lap
    pde = pressure[1:-1, :-1] - (core.conj()*fww).real/2
    return direct[1:-1, :-1], pde


def sphere_means(values, x, radial, radii, order):
    """Average a scalar over the in-brane directions of a sphere at fixed w."""
    interp = RegularGridInterpolator((x, np.r_[-radial[0], radial]),
                                     np.column_stack((values[:, 0], values)))
    mu, weights = leggauss(order)
    directions = np.column_stack((mu, np.sqrt(1-mu**2)))
    return np.array([np.dot(weights, interp(radius*directions))/2 for radius in radii])


def read_grid(tag, hashes):
    path = Path(f"_scratch/s11c/lateral_bath_{tag}_01/fields/d0.npz")
    with verified_packet(path, hashes) as packet:
        f, currents = packet["field"], packet["currents"]
        x, radial, dz = packet["x"], packet["r"], float(packet["dz"])
        isolated = int(packet["d"]) == 0 and len(currents) == 1
    dx = x[1] - x[0]
    rate = np.sum(currents[0]*2*np.pi*radial*dx**2)
    width = (len(f)-1)*dz
    radii = np.array([2, 4, 8, 16])
    print(f"grid={tag} intake={rate:.12g} prescribed_width={width:.12g}", flush=True)
    if not isolated or abs(rate-0.25) > 1e-12 or abs(width-1) > 1e-12:
        raise ValueError("Unexpected baseline drain")
    if np.any(currents[0][x[:, None]**2+radial**2 >= 1]):
        raise ValueError("Drain current extends onto a sampled reflecting face")
    face_means = []
    for face, first, second in ((f[-1], f[-2], f[-3]), (f[0], f[1], f[2])):
        both = normal_stress(face, first, second, radial, dx, dz)
        face_means.append([sphere_means(v, x[1:-1], radial[:-1], radii, 192) for v in both])
        coarse = sphere_means(both[0], x[1:-1], radial[:-1], radii, 96)
        if not np.allclose(coarse, face_means[-1][0], rtol=.005, atol=5e-13):
            raise ValueError("Shell sampling changed with angular resolution")
    upper, lower = np.array(face_means)
    width_load = (upper+lower)/2
    center_load = upper-lower
    for i, radius in enumerate(radii):
        print(f"  r={radius} upper={upper[0,i]:+.12g} lower={lower[0,i]:+.12g} "
              f"width={width_load[0,i]:+.12g} width_pde={width_load[1,i]:+.12g} "
              f"center={center_load[0,i]:+.12g}")
    consistent = all(np.allclose(direct, pde, rtol=.05, atol=5e-13) for direct, pde in face_means)
    if not np.isfinite(face_means).all() or not consistent:
        raise ValueError("The two normal-stress reconstructions disagree")


if __name__ == "__main__":
    symbolic_checks()
    runs = json.loads(Path("docs/_measurements/lateral_reservoir_pair.json").read_text())["runs"]
    hashes = {Path(path): digest for run in runs for path, digest in run.get("fields", {}).items()}
    for grid in ("e24_c4", "e24_c6", "e48_c4"):
        read_grid(grid, hashes)
