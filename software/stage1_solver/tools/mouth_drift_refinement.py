"""Fixed-extent radius-one grid refinement of the slow sideways mouth load."""

import argparse
from pathlib import Path

import numpy as np

from lateral_mouth_radius import SpeedTranslation, Target
from lateral_radius_momentum import slow_tangent
from lateral_reservoir_field import LateralSlab
from mouth_drift_flux_check import read_flux
from mouth_slow_drift import read_drift


def full_face_flux(packet):
    """Read the mixed flux across the entire top face from the four brane layers."""
    f, u = packet["field"], packet["u"]
    dx, dz = 2*packet["r"][0], packet["dz"]
    dv = 2*np.pi*packet["r"]*dx**2

    def dx_at_bath(q, bath):
        # Reflect about the two Dirichlet end faces to fill two ghost cells.
        extended = np.concatenate((2*bath-q[1::-1], q, 2*bath-q[:-3:-1]))
        return (extended[:-4]-8*extended[1:-3]+8*extended[3:-1]-extended[4:])/(12*dx)

    def normal(q):
        return (11*q[-1]-18*q[-2]+9*q[-3]-2*q[-4])/(6*dz)

    product = (dx_at_bath(u[-1], 0).conj()*normal(f)
               + dx_at_bath(f[-1], 1).conj()*normal(u)).real
    dn = 2*(f.conj()*u).real
    # The normal derivative of the density perturbation vanishes at each
    # endcap because its density is fixed at every depth.
    divergence = dx_at_bath(normal(dn), 0)/4
    return np.sum(dv*(product-divergence)), np.sum(dv*product)


def run(cells, path, normal_cells=None):
    if normal_cells is None:
        normal_cells = 4*cells
    s = LateralSlab(24, cells, normal_cells)
    f, currents, steps = s.solve(0)
    n0 = abs(f[-1])**2
    mean_speed = np.sum(s.dv*currents[0]/n0)/(4*np.pi/3)
    motion = SpeedTranslation(Target(s, 1., n0, mean_speed), f, currents[0])
    p, rp = motion.static_tangent()
    u, ru = slow_tangent(motion, 1j*p)
    packet = dict(x=s.x, r=s.r, dz=s.dz, radius=1., field=f[-4:], u=u[-4:],
                  normal=motion.normal_factor*motion.face)
    print(f"extent=24 cells={cells} normal_cells={normal_cells} steps={steps} "
          f"static_residual={rp:.12g} velocity_residual={ru:.12g}", flush=True)
    with path.open("xb") as stream:
        np.savez_compressed(stream, **packet)

    drift = read_drift(packet, 1.)
    flux = read_flux(packet)
    whole, product = full_face_flux(packet)
    lab = drift["lab"]
    comparison = dict(whole_layer=whole, whole_product=product,
                      gap_layer_halo=abs(flux["layer_halo"]-lab),
                      gap_layer_whole=abs(whole-lab),
                      gap_rim=abs(flux["layer_halo"]-flux["layer_cut"]))
    for group in (drift, flux, comparison):
        print(" ".join(f"{k}={v:.12g}" for k,v in group.items()), flush=True)

    if (not np.isfinite([*drift.values(), *flux.values(), *comparison.values()]).all()
            or abs(drift["rate"]-0.25) > 3e-14 or abs(drift["rate_slope"]) > 1e-10
            or abs(drift["rest"]) > 1e-11):
        raise ArithmeticError("Preparation, rest symmetry or velocity response failed")
    if not np.allclose([drift["fd_1em4"], drift["fd_5em5"]], lab, rtol=1e-6, atol=2e-13):
        raise ArithmeticError("Centered flux derivative differs from analytic derivative")
    if (abs(drift["interior"]-flux["layer_cut"]) > 5e-15
            or abs(flux["law_link"]-flux["law_halo"]) > 5e-15
            or abs(whole-product) > 5e-15):
        raise ArithmeticError("Mixed flux or integrated derivative identity failed")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("cells", type=int)
    parser.add_argument("output", type=Path)
    parser.add_argument("--normal-cells", type=int, help="Cells through the slab (default: 4 per lateral cell per radius)")
    args = parser.parse_args()
    if args.cells < 4:
        parser.error("At least four lateral cells per mouth radius required")
    if args.normal_cells is not None and args.normal_cells < 3:
        parser.error("The normal derivative needs four grid planes")
    run(args.cells, args.output, args.normal_cells)
