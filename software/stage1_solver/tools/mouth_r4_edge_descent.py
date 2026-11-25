"""Probe frozen step 108 along a registered entrance coordinate and a melted edge normal."""
import argparse
from fractions import Fraction as Q
from hashlib import file_digest
from itertools import product
import json
from pathlib import Path
import sys

import mpmath as mp
import numpy as np


STUDY_PIN = '22a06d53f4b64bd6d6be19fc2c8c2f8d038909ebd5f7357c9a327941d386820b'
STATE = 'a2210361aa5b38712de7db618b5220985e12f39f3a089351dfd3767bef0a7178'
CELL = (25, 11, 10)
ENTRANCE = (24, 2, 10)
BITS = list(product((0, 1), repeat=3))
NODES = [tuple(i+b for i, b in zip(CELL, bits)) for bits in BITS]


def digest(path):
    with path.open('rb') as stream:
        return file_digest(stream, 'sha256').hexdigest()


def cross(x, y):
    return tuple(x[(i+1)%3]*y[(i+2)%3] - x[(i+2)%3]*y[(i+1)%3] for i in range(3))


def dot(x, y):
    return sum(a*b for a, b in zip(x, y))


def edge_value(vertices, widths, t):
    """At (0,0,t), form Q1 chart columns and cofactors in exact binary-input arithmetic."""
    basis = []
    for bits in BITS:
        factors = [v if side else 1-v for v, side in zip((0, 0, t), bits)]
        basis.append(tuple((2*bits[j]-1)*factors[(j+1)%3]*factors[(j+2)%3]/widths[j]
                           for j in range(3)))
    columns = [[sum(v[i]*w[j] for v, w in zip(vertices, basis)) for i in range(3)]
               for j in range(3)]
    cof = [cross(columns[(j+1)%3], columns[(j+2)%3]) for j in range(3)]
    normal = [[sum(w[j]*cof[j][i] for j in range(3)) for i in range(3)] for w in basis]
    return dot(columns[0], cof[0]), normal


def edge_data(labels, widths, *, check_normal=False):
    vertices = [[Q(float(x)) for x in labels[node]] for node in NODES]
    p0, pm, p1 = [edge_value(vertices, widths, t)[0] for t in (Q(0), Q(1, 2), Q(1))]
    a = 2*(p1+p0-2*pm)
    b = p1-p0-a
    t = -b/(2*a) if a > 0 and 0 < -b < 2*a else Q(int(p1 < p0))
    minimum, normal = edge_value(vertices, widths, t)
    if minimum != a*t*t+b*t+p0:
        raise ArithmeticError('Quadratic edge reconstruction failed')
    if check_normal:
        if not 0 < t < 1:
            raise ArithmeticError('Saved limiting edge no longer has a strict interior minimum')
        # Changing one scalar node coefficient is a rank-one update at fixed t.
        for k in range(len(vertices)):
            for j in range(3):
                vertices[k][j] += 1
                difference = edge_value(vertices, widths, t)[0]-minimum
                vertices[k][j] -= 1
                if difference != normal[k][j]:
                    raise ArithmeticError('Cofactor derivative disagrees with exact rank-one update')
    return t, minimum, normal, (p0, b, a)


def incidence(joint, normal):
    rows = []
    for node, ds in zip(NODES, normal):
        incident = [tuple(i-b for i, b in zip(node, bits)) for bits in BITS]
        incident = [c for c in incident if all(0 <= k < n for k, n in zip(c, joint.opt.active.shape))]
        active = [c for c in incident if joint.opt.active[c]]
        groups = []
        flat = np.ravel_multi_index(node, joint.order.shape)
        for axis, cm in enumerate(joint.components):
            g = cm.groups[flat]
            state = 'held' if cm.held[g] else 'free_active' if cm.active[g] else 'free_melt_only'
            groups.append(dict(component=axis, group=int(g), status=state,
                               registered_nodes=int(cm.C.getcol(g).nnz)))
        rows.append(dict(node=node, chi=float(joint.order[node]), incident_cell_count=len(incident),
                         incident_positive_order_cells=active, groups=groups,
                         minimum_derivative=[float(v) for v in ds]))
    return rows


def run(study, area):
    pinned = area/'pin/software/stage1_solver/tools'
    manifest = study/'pinned.json'
    if digest(manifest) != STUDY_PIN or any(digest(pinned/k) != v for k, v in json.loads(manifest.read_bytes()).items()):
        raise ValueError('Required step-108 evaluator source pin differs')
    sys.path.insert(0, str(pinned))
    from mouth_exact_checkpoint_control import CHECKPOINT_SHA
    from mouth_r4_eulerian_corner_audit import point, profile
    from mouth_r4_eulerian_relax_labels import InfeasibleLabels
    from mouth_r4_eulerian_solid import InvalidSpatialCrystal, SpatialCrystal
    from mouth_r4_newton_held import build_joint
    from mouth_r4_refinement_review import _collar_reference
    from mouth_r4_smooth_held import SmoothHeld

    source = study/'primary/last_active.npz'
    saved = json.loads(source.with_suffix('.json').read_bytes())
    if digest(source) != STATE or saved['sha256'] != STATE or saved['outer'] != 108:
        raise ValueError('Expected the preserved 108th accepted step')
    driver = SmoothHeld(study.parent, area/'evaluation', 1.e-8)
    with np.load(source, allow_pickle=False) as f:
        if (str(f['parent_sha256']) != CHECKPOINT_SHA or float(f['epsilon']) != driver.epsilon
                or int(f['radius']) != 4 or int(f['winding']) != 1 or int(f['iteration']) != driver.iteration
                or not np.array_equal(f['x'], driver.x)):
            raise ValueError('Incorrect shape, gauge or checkpoint ancestry')
        driver.labels, driver.fields = f['labels'], f['field']
    labels = driver.labels
    trial = driver.trial()
    driver.current_trial = trial
    solid = SpatialCrystal(trial.b)
    joint, (value, grad, before) = build_joint(driver, trial, solid, registry_update=False)
    site = before['worst_free_active']
    if (tuple(site['node']) != ENTRANCE or site['component'] != 1
            or abs(site['signed_force']-saved['material']['worst_free_active']['signed_force']) > 1.e-11):
        raise ArithmeticError('Frozen entrance residual differs from step 108')
    i = int(np.argmax(np.where(joint.active, abs(grad), -np.inf)))
    widths = [Q(float(x[k+1]))-Q(float(x[k])) for x, k in zip(joint.opt.axes, CELL)]
    t, initial_minimum, nodal_normal, polynomial = edge_data(labels, widths, check_normal=True)
    nodal = np.zeros_like(labels)
    for node, vals in zip(NODES, nodal_normal):
        nodal[node] = list(map(float, vals))
    normal = joint.reduce(nodal)
    native_normal, native_det = joint.corner_normal(labels, CELL, (0., 0., float(t)))
    error = float(np.max(abs(native_normal-normal)))
    if error > 1.e-11*max(abs(normal)) or abs(native_det-float(initial_minimum)) > 1.e-13:
        raise ArithmeticError('Exact edge differential disagrees with the native Q1 cofactor')
    with mp.workdps(80):
        physical = point(trial.b, joint.opt.axes, labels, joint.order, CELL,
                         (0., 0., float(t)), _collar_reference)
    rows = incidence(joint, nodal_normal)
    # At the fixed minimizing chart location, Jinv = det(Dxi) * ref(xi)/phys(p).
    # Differentiate the reference volume, including angular hoop, on its native side.
    # The spatial volume is held fixed under a label variation at this site.
    with mp.workdps(75):
        mt = mp.mpf(t.numerator)/t.denominator
        weights = [1-mt, mt, 0, 0, 0, 0, 0, 0]
        xi = [sum(w*mp.mpf(float(labels[n][k])) for n, w in zip(NODES, weights)) for k in range(3)]
        def log_reference(x, depth, theta):
            if x > mp.mpf(float(trial.b.cut)):
                c = _collar_reference(trial.b, x, depth)
                volume = abs(mp.det(c))*c[2,2]
            else:
                r = profile(trial.b, depth, True)[0]*mp.mpf(float(trial.b.data.core_ratio))
                volume = r*mp.mpf(float(trial.b.data.ref_length))*(x*r)**2
            return mp.log(volume*mp.sin(theta))
        dlog = [mp.diff(log_reference, tuple(xi), tuple(int(j == k) for j in range(3))) for k in range(3)]
        ratio = physical['Jinv']/physical['chart_det']
        for row, ds, w in zip(rows, nodal_normal, weights):
            row['physical_Jinv_fixed_site_derivative'] = [
                ratio*float(n) + physical['Jinv']*float(w*q) for n, q in zip(ds, dlog)]

    unit = np.zeros_like(grad)
    unit[i] = 1.
    column = joint.hvp(labels, unit)
    curvature = float(column[i])
    if curvature == 0:
        raise ArithmeticError('Entrance coordinate has zero measured diagonal')
    coordinate = -grad[i]/abs(curvature)*unit

    # Work in a registered Q1 neighborhood of the entrance. A negative Hii makes
    # gradient descent on one coordinate increase its own residual infinitesimally.
    # The angular bisector of the two negative gradients (energy and |g_i|) gives
    # a common descent direction if their restrictions are not antiparallel.
    near = np.zeros(len(grad), bool)
    cursor = 0
    for cm in joint.components:
        mobile = np.flatnonzero(~cm.held)
        active = np.flatnonzero(cm.free[mobile])
        nodes = np.array(np.unravel_index(cm.source[mobile[active]], joint.order.shape)).T
        near[cursor+active] = (abs(nodes-ENTRANCE).max(axis=1) <= 1)
        cursor += len(mobile)
    scaled_energy = np.where(near, joint.width*grad, 0.)
    scaled_force = np.where(near, joint.width*np.sign(grad[i])*column, 0.)
    ehat = scaled_energy/np.linalg.norm(scaled_energy)
    fhat = scaled_force/np.linalg.norm(scaled_force)
    cosine = float(ehat@fhat)
    entrance = -(ehat+fhat)*joint.width
    if not np.any(entrance) or grad@entrance >= 0:
        raise ArithmeticError('No strict local bisector of the measured energy and force gradients')
    curvature_response = joint.hvp(labels, entrance)
    initial_rate = np.sign(grad[i])*curvature_response[i]
    if initial_rate < 0:
        scale = min(abs(.5*grad[i]/initial_rate), .2/np.max(abs(entrance/joint.width)))
    else:
        # A nonsmooth moving interface can give different column and row secants.
        # Keep the independently bounded finite probe and let the exact functional decide.
        scale = .2/np.max(abs(entrance/joint.width))
    entrance *= scale
    if normal@entrance != 0 or any(np.any(joint.lift(entrance)[node]) for node in NODES):
        raise ArithmeticError('Entrance-neighborhood tangent should leave the limiting cell unchanged')

    # The minimum-norm correction in cell-width-scaled packed coordinates satisfies
    # n.p = m at unit fraction, so the first-order edge margin becomes twice its held value.
    melt_metric_normal = np.where(joint.active, 0., joint.width**2*normal)
    capacity = float(normal@melt_metric_normal)
    correction = melt_metric_normal * (float(initial_minimum)/capacity) if capacity > 0 else np.zeros_like(grad)
    if np.any(joint.lift(correction)[joint.opt.stencil]):
        raise ArithmeticError('Melt-only correction moves the active stencil')
    directions = [('single_coordinate_control', coordinate), ('entrance_neighborhood', entrance)]
    if capacity > 0:
        directions.extend([('melt_only_control', correction), ('entrance_plus_melt', entrance+correction)])

    def summary(m):
        return {k: m[k] for k in ('elastic_energy', 'gauge_energy', 'max_free_label_force',
                                   'worst_free_active', 'gradient_l2', 'min_Jinv_global', 'min_Jinv_active',
                                   'max_unit_melt_residual', 'chart')}

    tests = {}
    def test(name, direction):
        move = joint.lift(direction)
        slope = float(grad@direction)
        dn = sum(dot(ds, map(lambda x: Q(float(x)), move[node])) for node, ds in zip(NODES, nodal_normal))
        record = dict(gradient_dot_direction=slope, edge_minimum_derivative=float(dn),
                      changed_packed_active=int(np.count_nonzero(direction[joint.active])),
                      changed_packed_melt_only=int(np.count_nonzero(direction[~joint.active])), trials=[])
        tests[name] = record
        for k in range(14):
            alpha = 2.**-k
            candidate = labels+alpha*move
            u, m, _, _ = edge_data(candidate, widths)
            row = dict(alpha=alpha, edge_minimum=float(m), edge_minimum_delta=float(m-initial_minimum),
                       minimizing_edge_fraction=float(u))
            record['trials'].append(row)
            try:
                new_value, new_gradient, after = joint.evaluate(candidate)
            except (InfeasibleLabels, InvalidSpatialCrystal) as exc:
                row.update(refusal=str(exc), witness=getattr(exc, 'witness', getattr(exc, 'diagnostics', {})))
                continue
            de = after['elastic_energy']-before['elastic_energy']
            dg = float(driver.epsilon*joint.weighted.work(labels, candidate))
            delta = de+dg
            with mp.workdps(80):
                pp = point(trial.b, joint.opt.axes, candidate, joint.order, CELL,
                           (0., 0., float(u)), _collar_reference)
            unchanged_active = np.array_equal(candidate[joint.opt.stencil], labels[joint.opt.stencil])
            row.update(elastic_work=de, gauge_work=dg, joint_work=delta,
                       direct_joint_difference=new_value-value, physical_Jinv_at_edge_minimum=pp['Jinv'],
                       entrance_signed_force=-float(new_gradient[i]), after=summary(after),
                       active_stencil_unchanged=unchanged_active, armijo=bool(slope < 0 and delta <= 1.e-4*alpha*slope))
            if abs((new_value-value) - delta) > 32*np.spacing(value):
                raise ArithmeticError('Native joint work disagrees with signed component work')
            if name == 'single_coordinate_control':
                return
            if name == 'melt_only_control':
                if not unchanged_active or de != 0:
                    raise ArithmeticError('Strictly melted control changes crystal energy or stencil')
                return
            if (row['armijo'] and delta < 0 and m >= initial_minimum
                    and abs(new_gradient[i]) < abs(grad[i]) and after['max_free_label_force'] < before['max_free_label_force']
                    and after['chart']['chart_det_global_lower'] > 0 and after['min_Jinv_global'] > 0):
                record['successful_fraction'] = alpha
                return

    for name, direction in directions:
        test(name, direction)
    # If the single-node probe transfers the peak to another entrance group, try
    # the restricted steepest descent over every active group on these depth planes.
    if not any('successful_fraction' in r for r in tests.values()):
        mask = np.zeros_like(grad, dtype=bool)
        cursor = 0
        for cm in joint.components:
            mobile = np.flatnonzero(~cm.held)
            active = cm.free[mobile]
            idx = np.unravel_index(cm.source[mobile[active]], joint.order.shape)
            mask[cursor+np.flatnonzero(active)] = np.abs(idx[1]-ENTRANCE[1]) <= 1
            cursor += len(mobile)
        d = -grad*mask
        hd = joint.hvp(labels, d)
        if d@hd > 0:
            d *= -(grad@d)/(d@hd)
        if any(np.any(joint.lift(d)[node]) for node in NODES):
            raise ArithmeticError('Restricted entrance support touches the limiting cell')
        test('entrance_layer_plus_melt', d+correction)

    base_inputs = [driver.prep, study.parent/'r4_whole_mouth01/coarse/endpoint.json',
                   study.parent/'r4_melted_survival01/coarse/state_00010.npz',
                   study.parent/'r4_cp10_labels01/labels.npz', source, source.with_suffix('.json')]
    return dict(step=108, state_sha256=STATE, source_manifest_sha256=STUDY_PIN,
                governing_script_sha256=digest(Path(__file__)), source_commit='1f46c4ad',
                inputs_sha256={str(f.relative_to(study.parent)):digest(f) for f in base_inputs},
                epsilon=driver.epsilon, entrance_node=ENTRANCE, packed_entrance_index=i,
                entrance_diagonal_curvature=curvature, single_coordinate_displacement=float(coordinate[i]),
                local_energy_abs_force_gradient_cosine=cosine,
                local_group_count=int(near.sum()), entrance_unit_fraction_displacement=float(entrance[i]),
                entrance_abs_force_linear_rate_from_Hd=float(initial_rate*scale),
                entrance_abs_force_linear_rate_from_He_i=float(np.sign(grad[i])*(column@entrance)),
                cell=CELL, edge_polynomial_ascending=[float(z) for z in polynomial],
                exact_edge_polynomial_ascending=list(map(str, polynomial)), exact_edge_minimum=str(initial_minimum),
                exact_minimizing_edge_fraction=str(t), minimizing_edge_fraction=float(t), edge_minimum=float(initial_minimum),
                physical_Jinv_at_edge_minimum=physical['Jinv'], reference_side=physical['material_side'],
                vertices=rows, packed_normal_check_max_error=error, native_chart_det_minus_exact=native_det-float(initial_minimum),
                width_scaled_melt_normal_squared=capacity, before=summary(before), tests=tests)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('study', type=Path)
    parser.add_argument('area', type=Path)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    result = run(args.study.resolve(), args.area.resolve())
    with args.output.open('x') as stream:
        json.dump(result, stream, indent=2, allow_nan=False)
        stream.write('\n')
    print(json.dumps(result, allow_nan=False))
