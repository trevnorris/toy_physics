"""Compare two laboratory-frame representations of the R4 receiver order.

The measured Q1 coefficient at zero displacement is identical in both.
Only its representation after a change to the carrier grid differs.
"""

import argparse
import json
from pathlib import Path

import numpy as np

from mouth_closed_rim_surface_review import quadrature
from mouth_last_ring_solid import interpolation
from mouth_rim_bank_work import Bank
from mouth_rim_disk_work_review import Disk
from mouth_rim_exterior_pin import projection_energy


def coefficient_gradient(model, values, r, z):
    r, z = np.broadcast_arrays(r, z)
    radius, slope = model.shape.radius(z)
    x = (r - radius) / (28. - radius)
    _, cx, cz = interpolation(model.z[1], model.rr[1], values, z, x)
    radial = cx / (28. - radius)
    return np.stack((radial, np.zeros_like(radial), cz - (1-x)*slope*radial), axis=-1)


def pull(gradient, radial, axial, h):
    result = gradient.copy()
    result[:, 2] /= 1 + h*axial
    result[:, 0] -= h*radial*result[:, 2]
    return result


def measurement(disk, edge, knots):
    m = disk.model
    bank = Bank(disk, 1, knots, edge)
    x, z = np.meshgrid(m.rr[1], m.z[1])
    radius = m.shape.radius(z)[0]
    r = radius + x*(28.-radius)
    eta = bank.motion(r,z)[0]
    outside = m.terminal_profile.outside
    baseline = outside(r,z,False)
    for field in disk.reader.fields[1,1]:
        assert np.max(abs(field.values[...,2] - baseline.T[:,None,:])) < 1e-12
    steps = (-2e-5,2e-5,-1e-5,1e-5)
    nodal = [outside(r,z+h*eta,False) for h in steps]
    xs, wx = quadrature(np.unique(np.r_[m.rr[1][m.rr[1]<edge],edge]),2)
    aa, wa = quadrature(m.a,2)
    x, angle = np.meshgrid(xs,aa,indexing='ij')
    wx = wx[:,None]*wa
    lo, _, hi = knots
    zs, wz = quadrature(np.unique(np.r_[lo,m.z[1][(m.z[1]>lo)&(m.z[1]<hi)],hi]),2)
    results = np.zeros((len(steps),2))
    max_base = 0.
    for z, w in zip(zs,wz):
        b = m.shape.radius(z)[0]
        r = b + x*(28.-b)
        weight = (4*np.pi*w*(28.-b)*r*r*np.sin(angle)*wx).ravel()
        g,gr,gz = (arr.ravel() for arr in bank.motion(r,z))
        rflat = r.ravel()
        zero_grad = coefficient_gradient(m,baseline,rflat,z)
        for sign,(state,*_) in zip((1.,-1.),bank.physical(x,angle,z)):
            gp,gc,n=state[1:4]
            max_base=max(max_base,float(np.max(abs(zero_grad-gc))))
            for i,h in enumerate(steps):
                p = pull(gp,gr,gz,h)
                gradients = (pull(coefficient_gradient(m,nodal[i],rflat,z),gr,gz,h),
                             coefficient_gradient(m,baseline,rflat,z+h*g))
                carried = projection_energy(p,pull(zero_grad,gr,gz,h),n)
                for j,c in enumerate(gradients):
                    dk = projection_energy(p,c,n)-carried
                    results[i,j] -= sign*((weight*(1+h*gz))@dk)/h
    if max_base>2e-10:
        raise ArithmeticError(f'Receiver gradient at zero differs: {max_base}')
    return dict(steps=steps,baseline_gradient_error=max_base,
                reaction_per_step=results.tolist(),
                columns=['resampled_to_moved_nodes','archived_Q1_fixed_at_physical_quadrature'])


if __name__ == '__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('state','wall','fields','axes'):
        p.add_argument(name,type=Path)
    a=p.parse_args()
    d=Disk(a.state,a.wall,a.fields)
    with np.load(a.axes) as f:
        edges=f['outside_x'][[4,8]]
        knots=f['duct_z'][[-9,-5,-1]]
    print(json.dumps(dict(mesh=d.level,source=d.source_hash,
                          by_width={k:measurement(d,e,knots) for k,e in zip(('4','8'),edges)}),
                     allow_nan=False))
