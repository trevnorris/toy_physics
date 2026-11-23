"""Rematch the fine R4 field and wound rest with the order pins moved to existing planes.

Only the free order mask changes. The original functional, coordinates, fixed
material and source assembly are used on the same instance before and after.
"""
import argparse
import json
from math import fsum
from pathlib import Path
import sys

import numpy as np
from scipy.sparse.linalg import factorized

from mouth_whole_r4_order_edges import charts, digest, profile, read


def planes(model, target):
    L = model.shape.length
    return {s: int(np.argmin(abs(model.z[s]-L-sign*target))) for s, sign in ((1, -1), (2, 1))}


def expand_mask(model, new_edges):
    from mouth_breathing_mobile_patch import KAPPA, WELL
    old = model.masks[2].copy()
    mask = np.zeros_like(old)
    for slab, edge in new_edges.items():
        window = slice(edge+1, None) if slab == 1 else slice(None, edge)
        mask[model.ids[slab, 0][:-1, :, window]] = True
    if np.any(old & ~mask) or not np.any(mask & ~old):
        raise ArithmeticError('Expanded mask must strictly contain the original')
    if any(np.any(mask[model.ids[s, 1]]) for s in (0, 1, 2)):
        raise ArithmeticError('Wider order crossed a prescribed lateral bank')
    if model.z[1][new_edges[1]] <= .62*model.shape.length:
        raise ArithmeticError('Wider support reaches the polar deformation')

    # The order metric in Joined.matrix is formed before the polar pullback.
    # Its restriction to these nodes lies entirely beyond the polar support.
    V = model.V[:, mask]
    metric = V.T@V.multiply((2*WELL*model.hq)[:, None])
    for G in model.G:
        grad = G[:, mask]
        metric += grad.T@grad.multiply((KAPPA*model.hq)[:, None])
    model.masks[2] = mask
    model.slices = np.r_[0, np.cumsum(model.masks.sum(axis=1))]
    model.solvers[2] = factorized(metric.tocsc())
    return dict(narrow=int(old.sum()), wider=int(mask.sum()), released=int((mask & ~old).sum()))


def edge_grad(model, gr, edges):
    output = {}
    for slab, j in edges.items():
        ids = np.unique(model.ids[slab, 0][:-1, :, j])
        t = gr[2, ids]
        output[str(slab)] = dict(peak_abs_node=float(abs(t).max()), signed_sum=float(t.sum()))
    return output


def limits(model, chart, edges):
    values = []
    for slab, j in edges.items():
        window = slice(j, None) if slab == 1 else slice(None, j+1)
        values.append(chart[slab][1]['wound'][:, :, window, 2])
    return [float(min(c.min() for c in values)), float(max(c.max() for c in values))]


def row_profile(model, chart, annulus, old_edges, new_edges):
    output = {}
    for place, frac in (('axis', 0.), ('annulus', annulus), ('near_rim', .98), ('bonded_rim', 1.)):
        row = profile(model, chart, frac)
        subset = row[abs(row[:, 0]) <= .25]
        # charts() receives identical arrays in both slots for a single state.
        order = subset[:, 3]
        spans = (subset[:-1, 0] >= model.z[1][new_edges[1]]-model.shape.length-1.e-13) & (
            subset[1:, 0] <= model.z[2][new_edges[2]]-model.shape.length+1.e-13)
        half = []
        for i in np.flatnonzero(spans):
            z0, z1 = subset[i:i+2, 0]
            y0, y1 = order[i:i+2]-.5
            if abs(y0)<5.e-14 and abs(y1)<5.e-14:
                if half and abs(half[-1][1]-z0)<1.e-13:
                    half[-1][1] = z1
                else:
                    half.append([z0, z1])
            elif y0*y1 < 0:
                cross = z0-(z1-z0)*y0/(y1-y0)
                half.append([cross, cross])
            elif abs(y0)<5.e-14 or abs(y1)<5.e-14:
                cross = z0 if abs(y0)<5.e-14 else z1
                if not half or abs(half[-1][1]-cross)>1.e-13:
                    half.append([cross, cross])
        edge = {}
        for key, plane in (('original', old_edges), ('wider', new_edges)):
            edge[key] = [float(order[np.argmin(abs(subset[:, 0]-(model.z[s][j]-model.shape.length)))])
                         for s, j in plane.items()]
        output[place] = dict(radial_fraction=float(frac), old_and_new_edge_chi=edge,
                             half_order_offset_ranges=half, columns=['w_minus_L', 'native_reference','total_chi'],
                             profile=subset[:, [0, 1, 3]].tolist())
    return output


def run(args):
    from mouth_exchange_carrier import BIAS
    from mouth_free_shape_solid import Solid
    from mouth_joined_front_audit import Front
    from mouth_last_rings_fluid import GRIDS
    from mouth_mobile_joined_fields import JointFields
    from mouth_rim_pinned_comparison import validate_guard
    folder = args.run/'fine'
    validate_guard(folder)
    report = read(folder/'endpoint.json')
    audit = read(args.edge_readback)['results']['fine']
    checkpoint = folder/report['file']
    controls = folder/'endpoint_controls.npz'
    if (digest(folder/'endpoint.json') != audit['endpoint_sha256'] or
        digest(checkpoint) != audit['state_sha256'] or digest(controls) != audit['controls_sha256']):
        raise ValueError('The endpoint and order-edge readback have different inputs')
    command = read(folder/'invocation.json')['command']
    state = Path(command[command.index('-u')+2])
    if digest(state) != audit['preparation_sha256']:
        raise ValueError('Order preparation changed')
    with np.load(checkpoint) as f, np.load(controls) as c:
        x, q, narrow = f['x'], f['q'], dict(wound=f['field'], rest=c['rest'])
    base = Front(state, GRIDS['double'], 1)
    f = JointFields(base, state, x, precondition=True)
    m = f.model
    o = audit['solid_discretization']
    solid = Solid(state, m.shape, m.polar_shift, o['n'], o['k'], o['angles'], m,
                  degree=o['degree'], subdivisions=o['subdivisions'])
    es, _, load, _ = solid.evaluate(q, gradient=False)
    static = fsum((f.static, es))
    old_edges, new_edges = planes(m, .12), planes(m, .20)
    narrow_energy = {k: f.energy(v, load) for k, v in narrow.items()}
    narrow_stationarity = {k: f.residual(v, load) for k, v in narrow.items()}
    if (abs(fsum((static,narrow_energy['wound']))-report['energy']) > 2.e-9 or
        max(narrow_stationarity.values()) >= 2.e-7):
        raise ArithmeticError('Narrow wound energy or stationary controls do not match')
    counts = expand_mask(m, new_edges)
    start_energy = {k: f.energy(v, load) for k, v in narrow.items()}
    if start_energy != narrow_energy:
        raise ArithmeticError('The expanded mask changed the starting energy')
    log = dict(command=sys.argv, reader_sha256=digest(Path(__file__)),
               profile_reader_sha256=digest(Path(__file__).with_name('mouth_whole_r4_order_edges.py')),
               source_manifest_sha256=digest(args.run/'sources.json'), edge_readback_sha256=digest(args.edge_readback),
               state_sha256=audit['state_sha256'], controls_sha256=audit['controls_sha256'],
               preparation_sha256=audit['preparation_sha256'], static_energy=static,
               order_nodes=counts, narrow_energy=narrow_energy, wider_starting_energy=start_energy,
               narrow_stationarity=narrow_stationarity, offsets={k: [float(m.z[s][j]-m.shape.length) for s,j in p.items()]
                  for k,p in (('narrow',old_edges),('wider',new_edges))}, comparisons={})
    for name, v in narrow.items():
        gradient = f.gradient(v, load)
        old = edge_grad(m, gradient, old_edges)
        rest = name=='rest'
        w, sol = f.solve(v, load, rest=rest)
        np.savez_compressed(args.scratch/(name+'.npz'), field=w)
        residual = f.residual(w, load)
        gradient = f.gradient(w, load)
        released = edge_grad(m, gradient, old_edges)
        current = m.measure(w)
        chart = charts(m, base.c, w, w)
        bounds = limits(m, chart, new_edges)
        zero_energy, narrow_g, _, _ = solid.evaluate(q, v[2], load=False)
        new_energy, wide_g, new_load, _ = solid.evaluate(q, w[2])
        if not np.array_equal(new_load, load):
            raise ArithmeticError('The fixed material deposition changed')
        narrow_solid = solid.residual(q, narrow_g)
        wide_solid = solid.residual(q, wide_g)
        dg_solid = solid.residual(q, wide_g-narrow_g)
        if not rest and abs(narrow_solid['material_force_max']-report['material']['material_force_max']) > 1.e-12:
            raise ArithmeticError('Narrow material load does not reproduce the saved endpoint')
        mu = float(np.mean(audit['annulus_fraction_cell']))
        comparison = dict(solver=sol, stationarity=residual, field_energy_change=sol['energy']-start_energy[name],
            max_order_correction_change=float(abs(w[2]-v[2]).max()), total_order_nodal_range=bounds,
            old_edge_narrow=old, old_edge_released=released, new_edge_held=edge_grad(m, gradient, new_edges),
            narrow_feed=report['controls'][name]['feed_current'], wider_feed=current['feed_current'],
            wider_bulk=current['bulk_current'], wider_weak_disk=current['weak_disk_current'],
            feed_change=current['feed_current']-report['controls'][name]['feed_current'],
            solid_energy_change=new_energy-zero_energy, narrow_solid=narrow_solid, wider_solid=wide_solid,
            solid_force_change_at_fixed_q=dg_solid, wider_profiles=row_profile(m, chart, mu, old_edges, new_edges))
        checks = dict(stationary=residual<2.e-7, bounded=bounds[0]>=-2.e-12 and bounds[1]<=1+2.e-12,
            descending=sol['energy']<=start_energy[name], ports=np.all(w[1,m.feed]==0.) and np.all(w[1,m.bulk]==(0. if rest else BIAS)),
            balance=current['port_error']<3.e-5 and abs(current['weak_disk_current']-current['feed_current'])<3.e-5,
            held_banks=np.count_nonzero(w[2,~m.masks[2]])==0,
            old_edge_unloaded=all(released[k]['peak_abs_node']<old[k]['peak_abs_node'] for k in old))
        comparison['checks'] = {k: bool(t) for k,t in checks.items()}
        log['comparisons'][name] = comparison
        (args.scratch/'readback.json').write_text(json.dumps(log, allow_nan=False)+'\n')
    log['accepted'] = all(all(c['checks'].values()) for c in log['comparisons'].values())
    args.output.write_text(json.dumps(log, indent=2, allow_nan=False)+'\n')
    if not log['accepted']:
        raise ArithmeticError('Held-width controls failed reported endpoint checks; see the recorded measurements')


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    for key in ('run', 'edge_readback', 'scratch', 'output'):
        p.add_argument(key, type=Path)
    args = p.parse_args()
    pin = args.run/'pin'
    for key, sha in read(args.run/'sources.json').items():
        if digest(pin/key) != sha:
            raise ValueError(f'Changed pinned input {key}')
    sys.path.insert(0, str(pin/'software/stage1_solver/tools'))
    args.scratch.mkdir(parents=True, exist_ok=False)
    run(args)
