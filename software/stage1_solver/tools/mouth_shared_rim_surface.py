"""Independent closed caps and sides for a bonded displacement of the whole R4 disk.

The slit traction and the historical boundary-node work are additional
measurements. Keeping the old wall-node hat fixed in physical space makes the
change in that work separable from changing the width of its support.
"""
import argparse
import json
from itertools import product
from pathlib import Path

import numpy as np

from mouth_closed_rim_endpoint import wall_reaction
from mouth_closed_rim_surface_review import Reader, quadrature
from mouth_exchange_carrier import mapped_nodes
from mouth_rim_disk_work_review import Disk
from mouth_supercore_joined_front import file_sha, state_from
from mouth_terminal_element_recovery import SECTORS, add_solid, add_wall, volume


def ramp(z, window):
    lo, left, right, hi = window
    return np.clip((z-lo)/(left-lo), 0, 1) * np.clip((hi-z)/(hi-right), 0, 1)


def segment(axis, lo, hi, degree):
    return quadrature(np.r_[lo, axis[(axis > lo) & (axis < hi)], hi], degree)


def outer_faces(reader, window, radial_window, degree):
    """Average actual caps and constant-chart-radius sides through each FE taper.

    At fixed exterior radial coordinate x the actual radius has derivative
    (1-x)R'(z); both components of the oriented side normal enter.
    """
    m = reader.model
    q1, q2 = radial_window
    a, wa = quadrature(m.a, degree)
    faces = {}
    for slab, height, sign in [(1, window[:2], -1.), (2, window[2:], 1.)]:
        zc, wz = segment(m.z[slab], *height, degree)
        for outer in (0, 1):
            xs, wx = segment(m.rr[outer], 0., q2 if outer else 1., degree)
            weight_x = np.clip((q2-xs)/(q2-q1), 0, 1) if outer else np.ones_like(xs)
            x, theta = np.meshgrid(xs, a, indexing='ij')
            cap = np.zeros(2)
            for z, w in zip(zc, wz):
                b = m.shape.radius(z)[0]
                width = 28.-b if outer else b
                r = b+x*width if outer else x*b
                area = 4*np.pi*r**2*np.sin(theta)*width*(weight_x*wx)[:,None]*wa
                flux = reader.sample(slab, outer, r, theta, z, np.array([0., 0., sign]))
                cap += w/(height[1]-height[0]) * (area.ravel()@flux)
            faces[f'cap_{slab}_{outer}'] = cap
        lo, hi = (window[0], m.shape.length) if slab == 1 else (m.shape.length, window[-1])
        zside, wz = segment(m.z[slab], lo, hi, degree)
        xs, wx = segment(m.rr[1], q1, q2, degree)
        x, theta = np.meshgrid(xs, a, indexing='ij')
        side = np.zeros(2)
        for z, w in zip(zside, wz):
            b, slope = m.shape.radius(z)
            r = b+x*(28.-b)
            area = 4*np.pi*r*r*np.sin(theta)*wx[:,None]*wa/(q2-q1)
            radial = reader.sample(slab, 1, r, theta, z, np.array([1., 0., 0.]))
            axial = reader.sample(slab, 1, r, theta, z, np.array([0., 0., 1.]))
            flux = radial - ((1-x)*slope).reshape(-1, 1)*axial
            side += w*ramp(z, window) * (area.ravel()@flux)
        faces[f'side_{slab}'] = side
    return faces


def material_nodes(disk, state, wall, zero, q2, degree):
    m = disk.model
    lo, hi = disk.windows[1][0], disk.windows[1][-1]
    regions = {}
    for slab, bank in product((1, 2), (0, 1)):
        z = m.z[slab]
        depth = np.flatnonzero((z[1:] > lo+1e-13) & (z[:-1] < hi-1e-13))
        radial = np.arange(len(m.rr[0])-1) if bank == 0 else np.flatnonzero(m.rr[1][:-1] < q2-1e-13)
        nodes = volume(disk, slab, bank, radial, depth, degree)
        if slab == 1:
            add_wall(disk, bank, depth, degree, nodes)
        regions[slab, bank] = nodes
    add_solid(disk, state, wall, zero, regions[1, 0])
    return regions


def nodal_works(m, regions, window, radial_window, reference_axes):
    """All-side reaction and two physical collar-hat profiles on this mesh."""
    q1, q2 = radial_window
    result = {key: np.zeros(len(SECTORS)) for key in ('all', 'native_hat', 'doubled_hat')}
    for (slab, bank), forces in regions.items():
        x, z = m.rr[bank], m.z[slab]
        p = np.clip((q2-x)/(q2-q1), 0., 1.) if bank else np.ones_like(x)
        a = ramp(z, window)
        result['all'] += np.einsum('i,j,ijk->k', p, a, forces)
        for key, axes in [('native_hat', (m.rr, m.z[2])), ('doubled_hat', reference_axes)]:
            rr, downstream = axes
            dz_tip = downstream[1] - downstream[0]
            inside = a if slab == 1 else np.clip(1-(z-m.shape.length)/dz_tip, 0., 1.)
            g = np.clip(1-x/rr[1][1], 0., 1.) if bank else np.clip((x-rr[0][-2])/(1-rr[0][-2]), 0., 1.)
            result[key] += np.einsum('i,j,ijk->k', g, inside, forces)
    return result


def measure(args):
    disk = Disk(args.state, args.wall, args.base)
    model = disk.model
    checks = json.loads(args.verified.read_text())
    source = next(c for c in checks['cases'] if c['inputs']['grid'] == disk.level)
    if (file_sha(args.zero) != source['zero_state_sha256'] or
            disk.source_hash['field'] != source['inputs']['field_sha256'] or
            file_sha(args.reference) != source['inputs']['reference_axes_sha256']):
        raise ValueError('Mismatched zero-state fields or comparison axes')
    with np.load(args.zero, allow_pickle=False) as saved:
        pair = [saved[k] for k in ('wound', 'rest')]
    disk.reader = Reader(model, state_from(args.state)[1], *pair)
    with np.load(args.reference, allow_pickle=False) as coarse:
        reference_axes = ((coarse['inside_x'], coarse['outside_x']), coarse['receiver_z'])
    radial_window = mapped_nodes(13, model.shape.length, 28., 22.)[0][1][4:6]
    # Direct solid stress differences are zero on the exterior caps and side.
    # The prescribed external order is identical in the driven/resting cases.
    dc = pair[0][2] - pair[1][2]
    inner = np.broadcast_to(model.z[1], model.ids[1,0].shape) < model.shape.length-model.h-1e-13
    upper = np.broadcast_to(model.z[2], model.ids[2,0].shape) > model.shape.length+model.h+1e-13
    support_check = max(np.max(abs(dc[model.ids[1,0][inner]])), np.max(abs(dc[model.ids[2,0][upper]])),
                        *(np.max(abs(dc[model.ids[slab,1]])) for slab in (1, 2)))
    if support_check != 0.:
        raise ValueError('Paired material stress extends into an exterior taper')
    results = []
    for degree in (2, 4):
        regions = material_nodes(disk, args.state, args.wall, args.zero, radial_window[-1], degree)
        for name, window in zip(('near', 'far'), disk.windows):
            faces = outer_faces(disk.reader, window, radial_window, degree)
            outgoing = sum(faces.values())
            slit = wall_reaction(disk.reader, lambda z: ramp(z, window), window[0], 8)
            work = nodal_works(model, regions, window, radial_window, reference_axes)
            if max(abs(work['all'][4]), abs(outgoing[0])) > 4e-5:
                raise ValueError('Unexpected net crystalline work or net outward current')
            results.append(dict(degree=degree, window=name, faces={k:v.tolist() for k,v in faces.items()},
                                outward=outgoing.tolist(), slit=slit.tolist(), sectors=SECTORS,
                                works={k:v.tolist() for k,v in work.items()},
                                common_reaction=float(work['all'].sum()),
                                common_plus_exterior=float(work['all'].sum()+outgoing[1]),
                                native_hat_plus_exterior=float(work['native_hat'].sum()+outgoing[1]),
                                doubled_hat_plus_exterior=float(work['doubled_hat'].sum()+outgoing[1]),
                                slit_plus_exterior=float(slit[:,1].sum()+outgoing[1])))
    return dict(grid=disk.level, source=disk.source_hash, zero_sha256=file_sha(args.zero),
                reference_sha256=file_sha(args.reference), max_external_order_difference=support_check,
                copies_of_common_rim=1, results=results)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    for key in ('state', 'wall', 'base', 'zero', 'reference', 'verified'):
        parser.add_argument(key, type=Path)
    print(json.dumps(measure(parser.parse_args()), allow_nan=False))
