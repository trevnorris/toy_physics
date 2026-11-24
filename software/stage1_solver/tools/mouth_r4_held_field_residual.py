"""Evaluate the saved R4 fields with the initial and last accepted crystal labels."""
import argparse
import hashlib
import json
from pathlib import Path
import sys
from tempfile import TemporaryDirectory

import numpy as np


def read(path):
    return json.loads(path.read_bytes())


def digest(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def peak(array):
    return float(np.max(np.abs(array), initial=0.))


def readout(area, output):
    pin = area / 'pin/software/stage1_solver/tools'
    hashes = read(area / 'sources.json')
    if any(digest(pin/name) != checksum for name, checksum in hashes.items()):
        raise ValueError('Changed held-solver source')
    sys.path.insert(0, str(pin))
    from mouth_continuous_depth_mesh import axes_for_depth
    from mouth_continuous_order import KAPPA, WELL
    from mouth_continuous_shape_relax import ShapeTrial
    from mouth_exact_checkpoint_control import CHECKPOINT_SHA
    from mouth_r4_eulerian_relax_labels import LabelRelaxation
    from mouth_r4_eulerian_release import Driver, order_grid
    from mouth_r4_eulerian_solid import SpatialCrystal
    from mouth_r4_melt_barrier import MeltBarrier
    from mouth_r4_slaved_melt import FLOW

    scratch = area.parent
    source = area/'primary/held_stop.npz'
    earlier = scratch/'r4_slaved_melt_barrier01/primary/fixed_active.npz'
    expected = read(area/'primary/held_stop.json')
    limits = read(area/'primary_guard/limit-validation.json')
    outcome = read(area/'primary_guard/outcome.json')
    control = read(scratch/'r4_barrier_work01/PRECISION_GATE.json')['cases']['primary']
    if (digest(source) != expected['sha256'] or digest(source) != digest(area/'primary/last_active.npz')
            or digest(earlier) != control['fixed_trace_sha256'] or not limits['verified']
            or int(limits['actual']['memory.swap.max']) != 0 or not outcome['limitsVerified']
            or outcome['exitCode'] != 1 or outcome['childOutcome']['guardReason'] is not None):
        raise ValueError('Unverified held study')
    with np.load(source) as st, np.load(earlier) as prior:
        if int(st['radius']) != 4 or int(st['winding']) != 1 or str(st['parent_sha256']) != CHECKPOINT_SHA:
            raise ValueError('Wrong scientific input')
        fields, labels, q = (st[k].copy() for k in ('field','labels','x'))
        old_labels, epsilon = prior['labels'].copy(), float(st['epsilon'])
    driver = Driver(scratch, area/'primary', epsilon)
    if not np.array_equal(q, driver.x) or epsilon != 1e-8:
        raise ValueError('Changed held shape or gauge')
    body_shape = driver.shape(q)
    trial = ShapeTrial(driver.base, driver.prep, q, 8, 24, 6, 3, 2, precondition=False,
                       axes=axes_for_depth(body_shape.length, driver.reference_depth))
    model = trial.model
    model.extended_energy_sum = True
    if not np.all(fields[1, model.feed] == 0.) or not np.all(fields[1, model.bulk] == .4510084282857874):
        raise ValueError('Changed phase reservoirs')
    order = order_grid(trial, fields)
    crystal = SpatialCrystal(trial.b)
    opt = LabelRelaxation(trial.b, labels, order, epsilon)
    chi = fields[2, trial.first]
    # The production obstacle projection uses the diagonal of
    # kappa G.T H G + 2 A V.T H V on the joined order space.
    trial.diag = sum(g.power(2).T @ (KAPPA*model.hq) for g in trial.G) + trial.V.power(2).T @ (2*WELL*model.hq)
    def evaluate(label):
        elastic, grad, order_source, _ = crystal.evaluate(label, order)
        inner = np.bincount(trial.b.ids.ravel(), weights=order_source.transpose(0,2,1).ravel(), minlength=model.size)
        model.material_load = trial.ext + inner
        energy, raw = model.evaluate(fields)
        joined = trial.pull.T @ raw[2]
        feasible = np.where(chi == 0., np.minimum(joined, 0.), np.where(chi == 1., np.maximum(joined, 0.), joined))
        site = int(np.argmax(np.where(trial.free, np.abs(feasible), 0.)))
        record = dict(physical_energy=float(energy), crystal_energy=elastic, source_exchange=float(elastic-inner @ fields[2]),
                      phase_raw=peak(raw[1,model.masks[1]]), amplitude_raw=peak(raw[0,model.masks[0]]),
                      order_raw=peak(joined[trial.free]), order_kkt=peak(feasible[trial.free]),
                      order_projected=peak(trial.projected(chi,joined)), largest_kkt_group=site,
                      largest_kkt_gradient=float(joined[site]), order_at_largest_kkt=float(chi[site]))
        return record, grad, raw
    before, _, pre_grad = evaluate(old_labels)
    after, material_grad, end_grad = evaluate(labels)
    earlier_report = expected['field']
    for key in ('phase_raw','amplitude_raw','order_raw','order_kkt','order_projected'):
        if abs(before[key] - earlier_report[key]) > 1e-9:
            raise ArithmeticError(f'Prior fixed-field residual failed reproduction at {key}')
    carrier_change = peak(pre_grad[:2]-end_grad[:2])
    if carrier_change != 0.:
        raise ArithmeticError('Carrier work changed with crystal label at fixed order and fixed fields')

    biggest = []
    for component, cm in enumerate(opt._components()):
        gradient = cm.sum(material_grad[...,component])
        group = int(np.flatnonzero(cm.free)[np.argmax(abs(gradient[cm.free]))])
        nodes = np.argwhere((cm.groups==group).reshape(order.shape))
        physical_nodes = []
        for i,j,k in nodes:
            x, t, angle = (ax[n] for ax,n in zip(opt.axes,(i,j,k)))
            w = float(t*body_shape.length)
            alpha = (t-.06)/.56
            bend = 16*alpha**2*(1.-alpha)**2 if 0 < alpha < 1 else 0.
            r = float(x*body_shape.radius(w)[0]+x**3*bend*trial.polar_shift*(3*np.cos(angle)**2-1)/2)
            physical_nodes.append(dict(index=[int(i),int(j),int(k)], radius=r, depth=w,
                                       theta=float(angle), chi=float(order[i,j,k])))
        biggest.append(dict(component=('radial','axial','angular')[component], group=group,
                            signed_label_force=float(-gradient[group]), nodes=physical_nodes))
    with TemporaryDirectory(dir=area.parent) as tmp:
        barrier = MeltBarrier(trial.b, opt, Path(tmp))
        unit, gauge_grad, _, _ = barrier.sample(labels)
        gauge = epsilon*float(unit)
        inner_res = peak(barrier.extract(gauge_grad,'dead'))
        active = peak(barrier.extract(material_grad+epsilon*gauge_grad,'free'))
    currents = trial.measure(fields,model.material_load)
    max_flow_error = max(abs(currents[key]-expected['current'][key]) for key in FLOW)
    m = expected['stop']['material']
    if (abs(after['physical_energy']-currents['excess_energy'])>1e-8 or max_flow_error>1e-11
            or abs(active-m['max_free_label_force'])>1e-11 or abs(inner_res-m['max_unit_melt_residual'])>1e-10
            or abs(gauge-m['gauge_energy'])>1e-15 or abs(after['crystal_energy']-m['elastic_energy'])>1e-12):
        raise ArithmeticError('Saved energy, force or current changed in stationary readback')
    result = dict(input_sha256=digest(source), original_fixed_trace_sha256=digest(earlier),
                  source_manifest_sha256=digest(area/'sources.json'), epsilon=epsilon, radius=4,
                  before_label_steps=before, after_label_steps=after, fixed_field_carrier_gradient_change=carrier_change,
                  carrier_metric_from_identical_fields=earlier_report['carrier_metric'], largest_free_forces_by_component=biggest,
                  total_fixed_field_energy=after['physical_energy']+gauge, gauge_energy=gauge,
                  gauge_to_crystal=gauge/after['crystal_energy'], unscaled_melt_residual=inner_res,
                  free_label_residual=active, currents={key:currents[key] for key in FLOW},
                  max_abs_change_from_saved_current=max_flow_error)
    with output.open('x') as fp:
        json.dump(result, fp, indent=2, allow_nan=False)
        fp.write('\n')
    print('Final physical order KKT, projected, gauge energy, total energy:',after['order_kkt'],after['order_projected'],gauge,result['total_fixed_field_energy'])


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('area', type=Path)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    readout(args.area.resolve(), args.output.resolve())
