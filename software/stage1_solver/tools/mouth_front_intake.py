"""Pressure-selected throughput and observable shape controls for a permeable order front."""
from itertools import product
from pathlib import Path

import numpy as np
import sympy as sp
from numpy.polynomial.legendre import leggauss

from mouth_carrier_momentum import read_state
from mouth_pinned_front_check import Material, branch
from mouth_release_input_check import manifest_fields


def mixed_drive():
    """Test a rest calibration with new normal heads, without prescribing a current."""
    material = Material()
    head, step = .002, .000025
    flow_low = branch(material, head-step, 22)
    flow_high = branch(material, head+step, 22)
    rest_slope = (flow_high-flow_low)/(2*step)
    print(f"rest_head={head:g} static_dD_dH={rest_slope:.12g}")
    # These are the two distinct heads among four signed trials for W=0 under
    # full transmission of the ideal exterior head. Solve each new head once.
    for u, v in ((.01, .02), (.005, .01)):
        aligned = branch(material, head + u*u/2 - u*v, 22)
        opposed = branch(material, head + u*u/2 + u*v, 22)
        measured = (aligned-opposed)/(2*u*v)
        print(f"U={u:g} V={v:g} same_sign_D={aligned:.12g} opposite_sign_D={opposed:.12g} "
              f"mixed_D_over_UV={measured:.12g} rest_slope_prediction={-rest_slope:.12g} "
              "fixed_local_head_mixed=0")
        assert abs(measured + rest_slope) < .001


def aperture_shape():
    # Dimensionless radius, R=1. For prescribed positive mass flux j(r), let
    # b=rho*nu*h. Then -(r² b')'=r² j, regular at the center and b(1)=0.
    r, s = sp.symbols("r s", positive=True)
    for description, density in (("uniform", sp.Integer(1)), ("archived_weight", (1-s*s)**3)):
        integral = sp.integrate(s**2*density, (s, 0, r))
        slope = sp.cancel(-integral/r**2)
        center = -sp.integrate(slope, (r, 0, 1))
        exact_rate = 4*sp.pi*integral.subs(r,1)
        center_rate = 8*sp.pi*center
        edge_rate = -4*sp.pi*slope.subs(r,1)
        print(f"radial_profile={description} center={center} true_rate={exact_rate} "
              f"center_estimate_over_true={sp.factor(center_rate/exact_rate)} "
              f"edge_estimate_over_true={sp.factor(edge_rate/exact_rate)}")
        assert sp.simplify(edge_rate-exact_rate)==0

    # Normal pressure traction on a graph (feed minus bulk, normal downstream).
    # For a fixed edge h=0 integration by parts gives Fx=integral h*d_x(dp).
    p0, gradient, radius, sigma = sp.symbols("p0 gradient R sigma", positive=True)
    height = p0*(radius**2-r**2)/(6*sigma)
    mean_height = sp.integrate(3*r**2*height/radius**3,(r,0,radius))
    print("normal_load_only Fx/(area*gradient*center_sag) =", sp.factor(mean_height/height.subs(r,0)))


def archived_shape():
    # Evaluate the exact continuum ball Green function at its center for the
    # archived *input profile* held constant in each grid cell. Clip boundary
    # cells at R=1. Also evaluate the enclosed rate on this identical support.
    hashes = manifest_fields(Path("docs/_measurements/mouth_feed_weak_run.json"))
    for n, dirname in ((24,"mouth_feed_weak01"), (32,"mouth_feed_depth01")):
        folder = Path("_scratch/s11c") / dirname / "fields"
        _, flux, dx, _ = read_state(folder, n, "restoring_slope", 0., hashes)
        axis = (np.arange(n) + .5 - n/2)*dx
        positions = np.stack(np.meshgrid(axis, axis, axis, indexing="ij"), axis=-1)
        mask = flux > 0
        cells = positions[mask]
        weights = flux[mask]
        total = dx**3 * weights.sum()
        radius2 = np.sum(positions**2, axis=-1)
        profile = np.maximum(1-radius2, 0)**3
        profile /= dx**3 * profile.sum()
        assert abs(total-.25) < 1e-12 and np.allclose(flux, total*profile, rtol=1e-12, atol=1e-14)
        for order in (4, 8):
            points, gauss = leggauss(order)
            mass, center_moment = 0., 0.
            for i,j,k in product(range(order), repeat=3):
                offset = dx/2 * points[[i,j,k]]
                r = np.linalg.norm(cells + offset, axis=1)
                dv = (dx/2)**3 * gauss[i]*gauss[j]*gauss[k]
                inside = r < 1
                weighted_volume = dv*weights[inside]
                mass += weighted_volume.sum()
                center_moment += weighted_volume @ (1/r[inside]-1)
            # rho*nu*h(0)=center_moment/(4*pi). Fixed-edge uniform estimate is
            # 8*pi*rho*nu*h(0) for R=1. True rate on the clipped ball is mass.
            ratio = 2*center_moment/mass
            print(f"archived_grid={n} cubature={order} saved_rate={total:.12g} "
                  f"ball_rate={mass:.12g} uniform_sag_estimate_ratio={ratio:.12g}")
            assert abs(mass/total-1) < .002 and 2.7 < ratio < 3.05


if __name__ == "__main__":
    mixed_drive()
    aperture_shape()
    archived_shape()
