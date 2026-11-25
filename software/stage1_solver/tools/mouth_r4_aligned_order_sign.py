"""Independent local +1e-9 energy sweep on both fixed-label R4 field rematches."""
import argparse
import json
from itertools import product
from pathlib import Path
import sys

import numpy as np
from scipy import sparse

from mouth_r4_aligned_rematch import FORMER, STRENGTH, digest, material_nodes, pin_sources, save_json


def gauge_polynomial(t, labels, chi):
    """Independently integrate each joined basis's exact cubic weight difference."""
    from mouth_cp10_label_condition import physical_frame
    from mouth_cp10_label_gauge import VERTICES, gradient_basis
    from mouth_r4_aligned_stencil import barrier_means
    from mouth_r4_eulerian_release import block_axes
    axes = block_axes(t.b)
    sites = t.number[t.b.ids.transpose(0, 2, 1)]
    cells = np.indices(tuple(n-1 for n in sites.shape)).reshape(3, -1).T
    nodes = np.column_stack([sites[tuple((cells+v).T)] for v in VERTICES])
    corners = np.stack([labels[tuple((cells+v).T)] for v in VERTICES], axis=1)
    lo = np.column_stack([ax[cells[:,i]] for i,ax in enumerate(axes)])
    widths = np.column_stack([np.diff(ax)[cells[:,i]] for i,ax in enumerate(axes)])
    offsets = corners - np.stack([np.column_stack([a[cells[:,j]+v[j]] for j,a in enumerate(axes)])
                                  for v in VERTICES], axis=1)
    _, det_energy = barrier_means(corners, widths)
    rows = np.repeat(np.arange(len(cells)), 8)
    coeff = np.zeros((3, len(chi)), dtype=np.longdouble)
    unit = np.longdouble(0)
    q1d = (.5-.5/np.sqrt(3), .5+.5/np.sqrt(3))
    for q in product(q1d, repeat=3):
        frame, hoop = physical_frame(t.b, lo + np.asarray(q)*widths)
        volume = np.pi/2 * widths.prod(axis=1) * hoop * np.linalg.det(frame)
        chart = np.einsum('nvc,vd->ncd', offsets, gradient_basis(q)) / widths[:,None,:]
        gradient = chart @ np.linalg.inv(frame)
        density = volume * (det_energy + .5*np.sum(gradient**2, axis=(1,2)))
        phi = np.where(VERTICES, q, 1-np.asarray(q)).prod(axis=1)
        N = sparse.coo_matrix((np.tile(phi,len(cells)), (rows,nodes.ravel())),
                              shape=(len(cells),len(chi)), dtype=np.longdouble).tocsr()
        v = np.clip(N @ chi, 0, 1)
        unit += np.sum(density*(1-v)**2*(1+2*v))
        coeff[0] += N.T @ (density*(-6*v*(1-v)))
        coeff[1] += N.power(2).T @ (density*(6*v-3))
        coeff[2] += N.power(3).T @ (2*density)
    return unit, coeff


def order_work(t, labels, u, eps):
    """Sweep every joined free zero at the supplied fields and physical model.material_load.

    Return the strict sign readout and full arrays, with local physical energy and
    a separate exact cubic gauge difference. The same entry applies after crystal motion.
    """
    from mouth_r4_smooth_order_boundary import PhysicalOrderEnergy
    physical_energy, physical_gradient = t.model.evaluate(u)
    physical_gradient = t.pull.T @ physical_gradient[2]
    chi = u[2, t.first]
    if not np.array_equal(u[2], chi[t.number]):
        raise ValueError('The joined order registry has split')
    ids = np.flatnonzero((chi == 0) & t.free)
    increment = np.longdouble(1e-9)
    physical = PhysicalOrderEnergy(t, u)
    slope = np.empty(len(ids), dtype=np.longdouble)
    change = np.empty_like(slope)
    for k,i in enumerate(ids):
        local = physical.at(i)
        slope[k] = local(0)
        change[k] = increment*local(increment)
    unit, coeff = gauge_polynomial(t, labels, chi)
    a,b,c = coeff[:,ids]
    gauge = eps*increment*(a+increment*(b+increment*c))
    total = change + gauge
    good = change > 0
    ratios = np.divide(abs(gauge), change, where=good, out=np.zeros_like(change))
    max_error = float(np.max(abs(slope-physical_gradient[ids])))
    fields = dict(joined=ids, physical=change, eps_gauge=gauge, total=total, physical_slope=slope)
    def row(k):
        node = ids[k]
        return dict(node=int(node), material_grid=material_nodes(t,node), depth=float(t.depth[node]),
                    physical=float(change[k]), gauge=float(gauge[k]), total=float(total[k]),
                    physical_slope=float(slope[k]), gauge_slope=float(eps*a[k]))
    def counts(v):
        return dict(positive=int((v>0).sum()), zero=int((v==0).sum()), negative=int((v<0).sum()),
                    nonfinite=int((~np.isfinite(v)).sum()))
    former = []
    for i in FORMER:
        local = physical.at(i)
        z = float(chi[i])
        item = dict(node=i, after=z,
                    material_grid=material_nodes(t,i), physical_slope=float(local(0)),
                    gauge_slope=float(eps*coeff[0,i]), total_slope=float(local(0)+eps*coeff[0,i]))
        if z+increment <= 1:
            item['physical_plus_step'] = float(increment*local(increment))
            x,y,z3 = coeff[:,i]
            item['gauge_plus_step'] = float(eps*increment*(x+increment*(y+increment*z3)))
            item['total_plus_step'] = item['physical_plus_step']+item['gauge_plus_step']
        former.append(item)
    result = dict(former=former, increment=float(increment), free_zero=len(ids), physical=counts(change), total=counts(total),
                  gauge=counts(gauge), long_double_mantissa_bits=int(np.finfo(np.longdouble).nmant),
                  min_physical=row(np.argmin(change)), min_total=row(np.argmin(total)),
                  worst_fraction=row(np.argmax(ratios)), max_abs_gauge_over_positive_physical=float(np.max(ratios)),
                  nonpositive=[row(k) for k in np.flatnonzero((total<=0) | ~good)],
                  native_physical_slope_max_error=max_error, physical_energy=float(physical_energy), unit_gauge=float(unit))
    result['passed'] = bool(len(ids) and np.isfinite(total).all() and np.all(good) and np.all(total>0)
                            and np.all(gauge<=0) and max_error<=5e-10)
    return result, fields


def sweep(study, case):
    prior = pin_sources(study)
    from mouth_continuous_depth_mesh import axes_for_depth
    from mouth_continuous_shape_relax import ShapeTrial
    from mouth_exact_checkpoint_control import CHECKPOINT_SHA
    from mouth_r4_eulerian_solid import SpatialCrystal
    from mouth_r4_newton_stop import read_state
    from mouth_r4_smooth_fields import material_load
    from mouth_r4_smooth_held import SmoothHeld

    folder = study/case
    record = json.loads((folder/'rematched.json').read_bytes())
    path = folder/'rematched.npz'
    eps = STRENGTH[case]
    if not record['field_matched'] or record['sha256'] != digest(path) or record['epsilon'] != eps:
        raise ValueError('Unverified rematched state')
    driver = SmoothHeld(study.parent, folder, 1e-8)
    source, old, labels = read_state(prior/'primary', 'last_active', driver, CHECKPOINT_SHA)
    with np.load(path, allow_pickle=False) as f:
        if (str(f['source_sha256']) != source['sha256'] or str(f['parent_sha256']) != CHECKPOINT_SHA
                or float(f['epsilon']) != eps or int(f['stencil_points']) != 39
                or not np.array_equal(labels, f['labels']) or not np.array_equal(driver.x, f['x'])):
            raise ValueError('Rematch changed the fixed step-108 preparation')
        u = f['field']
    t = ShapeTrial(driver.base, driver.prep, driver.x, 8, 24, 6, 3, 2, precondition=False,
                   axes=axes_for_depth(driver.shape(driver.x).length, driver.reference_depth))
    t.model.extended_energy_sum = True
    t.model.material_load, _ = material_load(t, SpatialCrystal(t.b), labels, u)
    result, fields = order_work(t, labels, u, eps)
    archive = folder/'order_work.npz'
    with archive.open('xb') as stream:
        np.savez_compressed(stream, **fields, increment=result['increment'], epsilon=eps, source_sha256=record['sha256'])
    for row in result['former']:
        row['before'] = float(old[2, t.first[row['node']]])
    result.update(command=sys.orig_argv, source_sha256=source['sha256'], rematched_sha256=record['sha256'],
                  script_sha256=digest(Path(__file__)), archive_sha256=digest(archive), epsilon=eps,
                  unit_gauge_energy_error=result['unit_gauge']-record['final']['unit_gauge'],
                  native_physical_energy_error=result['physical_energy']-record['final']['physical_energy'])
    result['passed'] &= abs(result['unit_gauge_energy_error'])<=1e-11 and abs(result['native_physical_energy_error'])<=1e-10
    save_json(folder/'gate1.json', result)
    print(json.dumps(result,allow_nan=False), flush=True)
    if not result['passed']:
        raise ArithmeticError('The fixed-label one-sided order sweep failed')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('study', type=Path)
    parser.add_argument('case', choices=STRENGTH)
    args = parser.parse_args()
    sweep(args.study.resolve(), args.case)
