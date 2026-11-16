"""Measure the steady paired-minus-alone mouth load across an exterior cylinder."""

import argparse
import json
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.interpolate import RectBivariateSpline

from lateral_radius_momentum import link_budget
from lateral_reservoir_field import LateralSlab
from lateral_reservoir_readout import mouth_load
from lateral_source_rate import check, verified_packet
from mouth_drift_volume_balance import FaceReconstruction
from mouth_feed_weak_check import dx4


REACHES = (1.5, 2., 2.5)


def panels(lower, upper, grid):
    """Gauss rule on each piece of the cubic surface interpolant."""
    edges = np.r_[lower, grid[(grid > lower) & (grid < upper)], upper]
    nodes, weights = leggauss(16)
    half = np.diff(edges)[:, None]/2
    return ((edges[1:]+edges[:-1])[:, None]/2+half*nodes).ravel(), (half*weights).ravel()


def surface_flux(s, f):
    """Outward +x momentum through the endcaps and side, from exterior fields."""
    cuts = {}
    for R in REACHES:
        r, wr = panels(0., R, s.r)
        x, wx = panels(-R, R, s.x)
        cuts[R] = (r, 2*np.pi*r*wr, x, 2*np.pi*R*wx)
    result = {R: np.zeros(5) for R in REACHES}
    radial_axis = np.r_[-s.r[::-1], s.r]
    for layer, wz in zip(f, s.dz*s.sqrt_z**2):
        even = np.column_stack((layer[:, ::-1], layer))
        re = RectBivariateSpline(s.x, radial_axis, even.real)
        im = RectBivariateSpline(s.x, radial_axis, even.imag)

        def field(x, r, dx=0, dr=0):
            return re.ev(x,r,dx=dx,dy=dr)+1j*im.ev(x,r,dx=dx,dy=dr)

        for R, (r, wr, x, wx) in cuts.items():
            parts = result[R]
            for sign in (-1, 1):
                xc = np.full_like(r, sign*R)
                value, fx, fxx = field(xc,r), field(xc,r,1), field(xc,r,2)
                n = abs(value)**2
                carried = (value.conj()*fx).imag**2/n
                gradient = (abs(fx)**2-(value.conj()*fxx).real)/2
                parts[0] += sign*wz*np.dot(wr,(n**5-1)/5)
                parts[1] += sign*wz*np.dot(wr,carried)
                parts[2] += sign*wz*np.dot(wr,gradient-carried)
            rc = np.full_like(x,R)
            value, fx, fr = field(x,rc), field(x,rc,1), field(x,rc,0,1)
            carried = (value.conj()*fx).imag*(value.conj()*fr).imag/abs(value)**2
            shear = ((fx.conj()*fr).real-(value.conj()*field(x,rc,1,1)).real)/2
            parts[3] += wz*np.dot(wx,carried)
            parts[4] += wz*np.dot(wx,shear-carried)
    return result


def mouth_normal(s, face, currents, centers):
    g = np.zeros_like(face.real)
    for j, center in zip(currents, centers, strict=True):
        weight = np.maximum(1-((s.x[:,None]-center)**2+s.r**2),0)**3
        if np.any((j > 0) != (weight > 0)) or np.any(j < 0):
            raise ValueError("Archived intake has the wrong support")
        g += weight
    return (-g*(1-1/abs(face))+1j*currents.sum(axis=0)/abs(face)**2)*face


def link_difference(alone, pair, normal0, normal1, s):
    """Integrate the existing link-budget differential between the two fields.

    Bilinear fluxes are quadratic in the field. Pressure and potential terms
    have degree at most ten. Five Gauss points give the exact endpoint difference.
    Only the two endpoints represent physical steady solutions.
    """
    ix = np.flatnonzero(abs(s.x)<4)
    xs = slice(ix[0],ix[-1]+1)
    rs = slice(0,np.count_nonzero(s.r<4))
    f = alone[:,xs,rs]
    df = pair[:,xs,rs]-f
    g = normal0[xs,rs]
    dg = normal1[xs,rs]-g
    zero = np.zeros_like(f)
    budgets = {(R,m): np.zeros(6) for R in REACHES for m in (1,2)}
    nodes, weights = leggauss(5)
    for node, weight in zip((1+nodes)/2, weights/2):
        f_node, g_node = f+node*df, g+node*dg
        for (R,m), values in budgets.items():
            values += weight*link_budget(f_node,zero,df,g_node,dg,s.x[xs],s.r[rs],
                                         s.dx,s.dz,R,m)
    return budgets


def target_link(s, face, current):
    g = np.maximum(1-(s.x[:,None]**2+s.r**2),0)**3
    top = (-g*(1-1/abs(face))+1j*current/abs(face)**2)*face
    return -np.sum(s.dv*(dx4(face,s.dx).conj()*top).real)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("reference", type=Path)
    args = parser.parse_args()
    archive = json.loads(Path("docs/_measurements/lateral_reservoir_pair.json").read_text())
    hashes = {Path(p):h for row in archive["runs"] for p,h in row.get("fields",{}).items()}
    with verified_packet(args.reference/"d0.npz",hashes) as packet:
        x, r, dz = packet["x"], packet["r"], float(packet["dz"])
        s = LateralSlab(round(r[-1]+r[0]),round(1/(2*r[0])),round(1/dz))
        if (int(packet["d"]) or dz!=s.dz or not np.array_equal(x,s.x)
                or not np.array_equal(r,s.r)):
            raise ValueError("Mismatched isolated grid")
        alone, currents0 = packet["field"], packet["currents"]
    normal0 = mouth_normal(s,alone[-1],currents0,[0])
    res0, bath0, rates0 = check(s,alone,normal0,currents0,np.array([.25]))
    flux0 = surface_flux(s,alone)
    read0 = mouth_load(alone[-4:],currents0,x,r,dz,1)[:2].sum()
    link0 = target_link(s,alone[-1],currents0[0])
    symmetry = max(np.max(abs(v)) for v in flux0.values())
    print(f"extent={r[-1]+r[0]:g} cells={1/s.dx:g} normal_cells={1/dz:g} "
          f"alone_intake={rates0[0]:.12g} alone_bath={bath0:.12g} alone_residual={res0:.12g} "
          f"alone_mouth={read0:.12g} alone_surface_max_component={symmetry:.12g}",flush=True)
    if max(abs(read0),abs(link0),symmetry)>1e-11:
        raise ArithmeticError("Isolated target failed reflection control")
    reconstruction = FaceReconstruction(s)
    forces = {}
    for d in (4,8):
        with verified_packet(args.reference/f"d{d}.npz",hashes) as packet:
            if (int(packet["d"])!=d or float(packet["dz"])!=dz or
                    not np.array_equal(packet["x"],x) or not np.array_equal(packet["r"],r)):
                raise ValueError("Pair differs from isolated reference grid")
            saved, currents = packet["layers"], packet["currents"]
        pair, steps = reconstruction.field(saved[-1],1.,lambda f: -(abs(f)**8-1)*f/4)
        withheld = np.max(abs(pair[-4:-1]-saved[:-1]))
        normal = mouth_normal(s,pair[-1],currents,[-d,0])
        residual, bath, rates = check(s,pair,normal,currents,np.array([.25,.25]))
        at_mouth = mouth_load(saved,currents,x,r,dz,2)[:2].sum()-read0
        at_link = target_link(s,pair[-1],currents[-1])-link0
        print(f"d={d} steps={steps} withheld={withheld:.12g} residual={residual:.12g} "
              f"source_rate={rates[0]:.12g} target_rate={rates[1]:.12g} bath={bath:.12g} "
              f"mouth={at_mouth:.12g} mouth_link4={at_link:.12g}",flush=True)
        if not np.isfinite(withheld) or withheld>3e-11:
            raise ArithmeticError("Reconstruction failed withheld-layer check")
        surface = surface_flux(s,pair)
        links = link_difference(alone,pair,normal0,normal,s)
        for R in REACHES:
            if not (1 < R-s.dx and R+s.dx < d-1):
                raise ValueError("Cylinder links must isolate the target")
            parts = surface[R]-flux0[R]
            tensor = parts.sum()
            print(f"d={d} R={R:g} cap_pressure={parts[0]:.12g} cap_carried={parts[1]:.12g} "
                  f"cap_gradient_stress={parts[2]:.12g} side_carried={parts[3]:.12g} "
                  f"side_gradient_stress={parts[4]:.12g} tensor={tensor:.12g}",flush=True)
            for order, b in ((2,links[R,1]),(4,(4*links[R,1]-links[R,2])/3)):
                _, pressure, gradient, side, mouth_x, body = b
                outer = pressure+gradient+side
                corrected = outer-body
                error = corrected+mouth_x
                print(f"d={d} R={R:g} order={order} outer={outer:.12g} grid_body={body:.12g} "
                      f"corrected={corrected:.12g} direct_link={-mouth_x:.12g} closure={error:.12g}",
                      flush=True)
                if (not np.isfinite([*b,*parts]).all() or abs(error)>2e-11 or
                        (order==4 and abs(mouth_x+at_link)>2e-12)):
                    raise ArithmeticError("Finite-grid stationary momentum balance failed")
                if order==4:
                    forces[d,R] = np.array([outer,tensor,corrected,at_mouth])
    for R in REACHES:
        ratios = forces[8,R]/forces[4,R]
        print(f"R={R:g} ratio_outer={ratios[0]:.12g} ratio_tensor={ratios[1]:.12g} "
              f"ratio_corrected={ratios[2]:.12g} ratio_mouth={ratios[3]:.12g} "
              f"outer_percent_from_quarter={100*(4*ratios[0]-1):.12g} "
              f"tensor_percent_from_quarter={100*(4*ratios[1]-1):.12g}",flush=True)


if __name__ == "__main__":
    main()
