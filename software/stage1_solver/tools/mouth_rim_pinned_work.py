"""Independent R4 fixed-source envelope derivative and bond compatibility audit.

Separately measure swept material at the sloping wall, prescribed order work,
and work from remapping each registered crystal sample at its fixed location.
"""
import argparse
import json
from pathlib import Path

import numpy as np

from mouth_bonded_front_adjoint import angular_weights
from mouth_closed_rim_surface_review import quadrature
from mouth_disk_pillbox import dot
from mouth_last_ring_solid import interpolation
from mouth_rim_bank_work import Bank
from mouth_rim_disk_work_review import Disk
from mouth_rim_pinned_sources import Pinned, amplitude
from mouth_supercore_joined_front import KAPPA, ORDER_INERTIA, WELL, file_sha, state_from


def trace_sweep(shape, knots, eps):
    # r is a radius in three directions, so cross-sectional volume is 4πr³/3.
    lo, mid, hi = knots
    breaks = np.unique(np.r_[np.linspace(lo, mid, 65), np.linspace(mid, hi, 65)])
    z, w = quadrature(breaks, 8)
    hat = np.minimum((z-lo)/(mid-lo), (hi-z)/(hi-mid))
    dz = np.where(z < mid, 1/(mid-lo), -1/(hi-mid))
    b, slope = shape.radius(z)
    other = shape.radius(z + eps*hat)[0]
    volume = 4*np.pi/3 * (b-other) * (b*b+b*other+other*other) * (1+eps*dz)
    normal = -eps * hat * slope / np.sqrt(1+slope**2)
    return dict(step=eps, exposed_volume=w@np.maximum(volume, 0.),
                excluded_volume=w@np.maximum(-volume, 0.),
                sampled_normal_motion_min=normal.min(), sampled_normal_motion_max=normal.max(),
                common_bond_axial_mismatch=abs(eps))


def pin_order_work(disk, pins, side, edges, knots, steps):
    m = disk.model
    xaxis, zaxis = m.rr[side], m.z[1]
    if np.ndim(edges):
        bounds = np.unique(np.r_[edges, xaxis[(xaxis >= min(edges)) if not side else (xaxis <= max(edges))]])
    else:
        bounds = np.unique(np.r_[edges, xaxis[(xaxis >= edges) if not side else (xaxis <= edges)]])
    # Integrate the independent nodal pin derivative against the driven and
    # resting physical stresses, on the same native two-point Gauss rule.
    x, wx = quadrature(bounds, 2)
    a, wa = quadrature(m.a, 2)
    zn, wz = quadrature(np.unique(np.r_[knots, zaxis[(zaxis > knots[0]) & (zaxis < knots[-1])]]), 2)
    x, a = np.meshgrid(x, a, indexing='ij')
    wx = wx[:, None] * wa
    base = Bank(disk, side, knots, edges[0] if np.ndim(edges) else edges)
    changes = [pins.ref_change(side, edges, knots, h).T/h for h in steps]
    effects = np.zeros((len(steps), 2))
    for z, weight in zip(zn, wz):
        r0, slope = disk.sh.radius(z)
        width = 28.-r0 if side else r0
        r = r0+x*width if side else x*width
        measure = (4*np.pi*weight*width*wx*r**2*np.sin(a)).ravel()
        lift = []
        for delta in changes:
            v, dx, dz = interpolation(zaxis, xaxis, delta, z, x.ravel())
            vr = dx/width
            vz = dz - slope*(1-x.ravel() if side else x.ravel())*vr
            lift.append((v, np.column_stack([vr, np.zeros_like(vr), vz])))
        for sign, (state, *_) in zip((1., -1.), base.physical(x, a, z)):
            _, gp, gc, n, _, c = state
            p = ORDER_INERTIA*dot(gp, gc)/(n+ORDER_INERTIA*dot(gc, gc))
            current = gp-p[:, None]*gc
            gorder = 2*WELL*c*(1-c)*(1-2*c)
            for k, (v, grad) in enumerate(lift):
                effects[k,0] -= sign*(measure@(gorder*v + KAPPA*dot(gc, grad)))
                effects[k,1] += sign*(measure@(n*p*dot(current, grad)))
    return effects


def measure(a):
    disk = Disk(a.state, a.wall, a.fields)
    m = disk.model
    with np.load(a.reference) as ref:
        knots = ref['duct_z'][[-9,-5,-1]]
        xin, xout = (ref[k] for k in ('inside_x','outside_x'))
    native_order = state_from(a.state)[1]
    pins = Pinned(m, a.state, a.wall, native_order, np.zeros(m.size))
    archive = json.loads(a.work.read_text())['data'][disk.level]
    if (any(archive['source'].get(k) != v for k,v in disk.source_hash.items())
            or archive['source']['reference'] != file_sha(a.reference)):
        raise ValueError('Mismatched archived work')
    steps = (-2e-5,2e-5,-1e-5,1e-5)
    with np.load(a.fields) as f:
        diff = f['wound'][2] - f['rest'][2]
        mean = np.einsum('raj,a->jr', diff[m.ids[1,0]], angular_weights(m.a))
    if np.any(mean[:, -1]):
        raise ValueError('Duct drive and resting order differ on the wall')
    rows=[]
    for side, widths in ((0,(4,8)),(1,(4,)),(1,(8,))):
        edges = tuple((xout[w] if side else xin[-1-w]) for w in widths)
        edge = edges if len(edges)>1 else edges[0]
        prior = [next(r for r in archive['results'] if r['bank']==('outside' if side else 'duct')
                      and r['strip']=='near' and r['quadrature']==2 and r['radial_cells']==w) for w in widths]
        differences=np.array([1.,-1.] if len(edges)>1 else [1.])
        carrier = differences@np.array([r['total_native']-r['solid_reaction'] for r in prior])
        advected = differences@np.array([r['solid_reaction'] for r in prior])
        order_work = pin_order_work(disk, pins, side, edge, knots, steps)
        analytic_crystal = 0.; fixed_crystal = np.zeros(len(steps))
        if side == 0:
            r, z, weight = pins.material_sites()
            R, slope = disk.sh.radius(z)
            x = r/R
            eta = amplitude(side, edge, knots, x, z)
            _, dcx, dcz = interpolation(m.z[1], m.rr[0], mean, z, x)
            analytic_crystal = weight@(eta*(dcz - x*slope*dcx/R))
            for k,h in enumerate(steps):
                delta, _ = pins.crystal_change(edge,knots,h)
                fixed_crystal[k] = -(delta@diff)/h
            if np.max(abs(fixed_crystal-analytic_crystal))>4e-5:
                raise ArithmeticError('Fixed-site analytic derivative and signed nodal load differ')
        rows.append(dict(bank=('duct','outside')[side], widths=widths, signed_steps=steps,
                         frozen_carrier_reaction=carrier, convected_crystal=advected,
                         fixed_site_crystal=analytic_crystal, signed_crystal_check=fixed_crystal.tolist(),
                         prescribed_order_static=order_work[:,0].tolist(), prescribed_order_carrier=order_work[:,1].tolist(),
                         predicted_pin_corrected_reaction=(carrier+analytic_crystal+order_work.sum(axis=1)).tolist()))
    return dict(source=disk.source_hash, signed_steps=steps, depth_knots=knots.tolist(),
                bond_sweep=[trace_sweep(disk.sh,knots,h) for h in steps], banks=rows,
                work_sha=file_sha(a.work), state_sha=file_sha(a.state), reference_sha=file_sha(a.reference))


if __name__ == '__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    for key in ('state','wall','fields','reference','work','output'):
        parser.add_argument(key,type=Path)
    args=parser.parse_args()
    args.output.write_text(json.dumps(measure(args),indent=2,allow_nan=False)+'\n')
