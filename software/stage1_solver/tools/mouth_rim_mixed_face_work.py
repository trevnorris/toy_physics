"""First-order work in the thin strips swept between carrier and order cells.

The order gradient jumps at element faces. A carrier Gauss point can stay in
its original order element under both signs even as a layer of material between
the two meshes changes element. Integrate the limiting change at each internal
face from its two native Q1 traces. Use the lab-frame wall and disk geometry.
"""
import argparse
import json
from itertools import product
from pathlib import Path

import numpy as np

from mouth_closed_rim_surface_review import Reader
from mouth_exchange_carrier import mapped_nodes
from mouth_rim_disk_work_review import Disk
from mouth_supercore_joined_front import ORDER_INERTIA, file_sha, state_from


def traces(model, values, ijk, coords, side, r, slope):
    x, angle, z = coords
    kx, ka, kz = ijk
    grids = (model.rr[side], model.a, model.z[1 if z.max() <= model.shape.length else 2])
    spacing = [np.diff(g)[i] for g, i in zip(grids, ijk)]
    fraction = [(v - g[i]) / step for g, i, v, step in zip(grids, ijk, coords, spacing)]
    total = np.zeros((x.size, 3))
    gradient = np.zeros((x.size, 2, 3))
    for i, j, k in product((0, 1), repeat=3):
        bits = (i, j, k)
        widths = [f if bit else 1-f for bit, f in zip(bits, fraction)]
        triple = values[kx+i, ka+j, kz+k]
        total += (widths[0]*widths[1]*widths[2])[:, None]*triple
        for d in range(3):
            coefficient = (2*bits[d]-1)*widths[(d+1)%3]*widths[(d+2)%3] / spacing[d]
            gradient[:, :, d] += coefficient[:, None]*triple[:, 1:]
    radius = model.shape.radius(z)[0]
    gradient[:, :, 0] /= (28-radius if side else radius)[:, None]
    gradient[:, :, 1] /= r[:, None]
    gradient[:, :, 2] -= (slope*(1-x if side else x))[:, None]*gradient[:, :, 0]
    return total, gradient


def projection(n, phase_gradient, order_gradient):
    pc = np.einsum("ij,ij->i", phase_gradient, order_gradient)
    cc = np.einsum("ij,ij->i", order_gradient, order_gradient)
    return -.5*ORDER_INERTIA*n*pc**2 / (n+ORDER_INERTIA*cc)


def relative_motions(x, z, slab, side, bounds, heights, exterior):
    """Carrier displacement minus material displacement in solver units."""
    tent = np.interp(z, heights, [0., 1., 0., 1., 0.], left=0., right=0.)
    edges = bounds[side]
    g = [np.clip(1-x/b if side else (x-b)/(1-b), 0., 1.) for b in edges]
    plateau = np.clip((exterior[1]-x)/(exterior[1]-exterior[0]), 0., 1.) if side else 1.
    stationary = tent*(g[0]-g[1]) if (slab, side)==(1, 0) else np.zeros_like(x)
    return [stationary, tent*(g[0]-plateau), tent*(g[1]-plateau)]


def integrate(disk, bounds, heights, exterior, order):
    nodes, weights = np.polynomial.legendre.leggauss(order)
    nodes, weights = (1+nodes)/2, weights/2
    m = disk.model
    # case, side of step (-,+), face axis (radial,depth), slab, bank
    result = np.zeros((3, 2, 2, 2, 2))
    for slab, side, axis in product((1,2), (0,1), (0,2)):
        grids = (m.rr[side], m.a, m.z[slab])
        indices = [np.arange(len(g)-1) for g in grids]
        indices[axis] = np.arange(1, len(grids[axis])-1)
        zz = grids[2]
        selected = indices[2]
        indices[2] = selected[(zz[selected+int(axis!=2)] > heights[0]) & (zz[selected] < heights[-1])]
        ii = tuple(v.ravel() for v in np.meshgrid(*indices, indexing="ij"))
        left = tuple(v-int(d==axis) for d, v in enumerate(ii))
        for a, b in product(range(order), repeat=2):
            free = iter((a, b))
            coords = []
            measure = np.ones(len(ii[0]))
            for d, (g, i) in enumerate(zip(grids, ii)):
                if d == axis:
                    coords.append(g[i])
                else:
                    q = next(free)
                    width = np.diff(g)[i]
                    measure *= width*weights[q]
                    coords.append(g[i]+nodes[q]*width)
            x, angle, z = coords
            radius, slope = disk.sh.radius(z)
            r = radius + x*(28-radius) if side else x*radius
            # Flux area projected on the axial direction, oriented along
            # increasing x or z, including the two brane angular dimensions.
            normal = -slope*(1-x if side else x) if axis==0 else np.ones_like(x)
            cross = 1. if axis==0 else (28-radius if side else radius)
            measure *= 4*np.pi*r*r*np.sin(angle)*cross
            motions = relative_motions(x,z,slab,side,bounds,heights,exterior)
            jumps = np.zeros((2,len(r)))
            F = disk.reader.vortex(r*np.sin(angle))[0]
            for sign, interpolator in zip((1.,-1.),disk.reader.fields[slab,side]):
                v0, g0 = traces(m, interpolator.values, left, coords, side, r, slope)
                v1, g1 = traces(m, interpolator.values, ii, coords, side, r, slope)
                if np.max(abs(v0-v1)) > 2e-10:
                    raise ArithmeticError("Two traces disagree in field value on an internal face")
                density = (F*v0[:,0])**2
                for index, gp in enumerate((g0[:,0],g1[:,0])):
                    jumps[index] += sign*(projection(density,gp,g1[:,1])-projection(density,gp,g0[:,1]))
            for k, velocity in enumerate(motions):
                relative = normal*velocity
                for j, step in enumerate((-1.,1.)):
                    # Positive relative sweep brings material from the plus
                    # side into the carrier cell on the minus side.
                    jump = np.where(step*relative>0, jumps[0], jumps[1])
                    result[k,j,int(axis==2),slab-1,side] -= measure @ (relative*jump)
    return dict(gauss=order, cases={name:dict(negative=v[0].sum(),positive=v[1].sum(),
                centered=v.sum()/2, breakdown=v.tolist()) for name,v in
                zip(("fixed_wall_difference","rim_four","rim_eight"),result)})


def run(a):
    previous = json.loads(a.record.read_text())
    case = next(c for c in previous["cases"] if c["inputs"]["grid"] == a.zero.parent.name)
    meta = case["inputs"]
    if (file_sha(a.zero) != case["zero_state_sha256"] or file_sha(a.fields) != meta["field_sha256"]
            or file_sha(a.axes) != meta["reference_axes_sha256"]):
        raise ValueError("Unmatched state or axes")
    disk = Disk(a.state,a.wall,a.fields)
    with np.load(a.zero) as data:
        disk.reader = Reader(disk.model,state_from(a.state)[1],data["wound"],data["rest"])
    bounds = (meta["duct_bounds"],meta["receiver_bounds"])
    with np.load(a.axes) as f:
        heights = np.r_[f["duct_z"][[-9,-5,-1]], f["receiver_z"][[4,8]]]
    exterior = mapped_nodes(13,disk.sh.length,28.,22.)[0][1][4:6]
    return dict(field=file_sha(a.zero), wall=file_sha(a.wall), heights=heights.tolist(), bounds=bounds,
                exterior=exterior.tolist(), readings=[integrate(disk,bounds,heights,exterior,q) for q in (2,4)])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for key in ("state","wall","fields","axes","zero","record"):
        parser.add_argument(key,type=Path)
    print(json.dumps(run(parser.parse_args()),allow_nan=False))
