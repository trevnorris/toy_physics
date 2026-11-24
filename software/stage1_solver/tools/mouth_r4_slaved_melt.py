"""Re-solve the frozen R4 harmonic gauge and measure work, orientation and currents."""
import argparse
import json
from pathlib import Path

import mpmath as mp
import numpy as np

from mouth_continuous_depth_mesh import axes_for_depth
from mouth_continuous_shape_relax import ShapeTrial
from mouth_exact_checkpoint_control import CHECKPOINT_SHA, digest, read
from mouth_r4_eulerian_corner_audit import point
from mouth_r4_eulerian_relax_labels import LabelRelaxation, absmax, chart_orientation
from mouth_r4_eulerian_release import Driver, order_grid
from mouth_r4_eulerian_solid import SpatialCrystal
from mouth_r4_refinement_review import _collar_reference


FLOW = ('feed_current', 'bulk_current', 'weak_disk_current', 'incident_slice_current',
        'receiving_slice_current', 'positive_slice_current', 'port_error', 'disk_minus_feed_current',
        'dense_total_subsonic_share', 'total_axial_stress')


def observe(trial, fields, opt, spatial, labels, order, cell, fractions):
    energy, label_gradient, order_load, diagnostics = spatial.evaluate(labels, order, True)
    inner_load = np.bincount(trial.b.ids.ravel(), weights=order_load.transpose(0,2,1).ravel(),
                             minlength=trial.model.size)
    fluid = trial.measure(fields, trial.ext + inner_load)
    offsets = (labels - opt.base).reshape(-1,3)
    gauge_gradient = opt.K @ offsets
    physical, joint, dead = [], [], []
    for k,component in enumerate(opt._components()):
        elastic = component.sum(label_gradient[...,k])
        harmonic = component.sum(gauge_gradient[:,k])
        physical.append(elastic[component.free])
        joint.append((elastic + opt.epsilon*harmonic)[component.free])
        dead.append(absmax(harmonic[component.dead]))
    gauge = opt.gauge_energy(labels)
    with mp.workdps(65):
        corner = point(trial.b, opt.axes, labels, order, cell, fractions, _collar_reference)
    return (dict(crystal_energy=energy, gauge_energy=gauge,
                 total_energy=float(fluid['excess_energy'])+gauge,
                 max_free_force=absmax(np.concatenate(joint)),
                 unit_harmonic_residual=dead, corner=corner,
                 whole_cell_chart=chart_orientation(opt.axes,labels,opt.active),
                 min_physical_at_quadrature=diagnostics['min_Jinv_global'],
                 min_physical_at_ordered_quadrature=diagnostics['min_Jinv_active'],
                 currents={name:float(fluid[name]) for name in FLOW}),
            np.concatenate(physical), np.concatenate(joint), order_load, label_gradient)


def audit(scratch, case):
    run = scratch/'r4_eulerian_stage2_01'
    base = read(Path(__file__).resolve().parents[3]/'research/pde_ledger_v3/_measurements/mouth_r4_eulerian_stage2.json')
    saved = base['cases'][case]
    guard = read(run/(case+'_guard')/'outcome.json')
    snapshot = run/case/'label_stall.npz'
    if (not guard['limitsVerified'] or guard['childOutcome']['guardReason'] is not None
            or guard['exitCode'] != 1 or digest(snapshot) != saved['saved_state_sha256']):
        raise ValueError('Unverified held state')
    with np.load(snapshot) as st:
        if int(st['radius']) != 4 or int(st['winding']) != 1 or str(st['parent_sha256']) != CHECKPOINT_SHA:
            raise ValueError('Wrong preparation')
        labels, fields, x = (st[k].copy() for k in ('labels','field','x'))
    driver = Driver(scratch,run/case,saved['epsilon'])
    if not np.array_equal(x,driver.x):
        raise ValueError('Shape changed from checkpoint 10')
    shape = driver.shape(x)
    trial = ShapeTrial(driver.base,driver.prep,x,8,24,6,3,2,precondition=False,
                       axes=axes_for_depth(shape.length, driver.reference_depth))
    trial.model.extended_energy_sum = True
    order = order_grid(trial,fields)
    opt = LabelRelaxation(trial.b,labels,order,saved['epsilon'])
    # Preserve the absolute binary trace as well as the production axis/fixture groups.
    lifted = opt._lift(opt._pack(labels))
    constrained = np.stack([~c.dead[c.groups].reshape(order.shape)
                            for c in opt._components()], axis=-1)
    updated = np.where(constrained, labels, lifted)
    spatial = SpatialCrystal(trial.b)
    corner = saved['accepted_corner']
    def measure(label):
        return observe(trial,fields,opt,spatial,label,order,corner['cell'],corner['local'])
    before, old_grad, old_joint, old_load, old_nodes = measure(labels)
    after, new_grad, new_joint, new_load, new_nodes = measure(updated)
    delta = dict(crystal_energy=after['crystal_energy']-before['crystal_energy'],
                 gauge_energy=after['gauge_energy']-before['gauge_energy'],
                 total_energy=after['total_energy']-before['total_energy'],
                 active_labels=absmax((updated-labels)[opt.stencil]),
                 melt_only_labels=absmax((updated-labels)[~constrained]),
                 lifted_trace_roundoff=absmax((lifted-labels)[constrained]),
                 all_active_elastic_nodal_forces=absmax((new_nodes-old_nodes)[opt.stencil]),
                 all_free_elastic_forces=absmax(new_grad-old_grad),
                 all_free_combined_forces=absmax(new_joint-old_joint),
                 order_conjugate_load=absmax(new_load-old_load),
                 currents={k:after['currents'][k]-before['currents'][k] for k in FLOW})
    if (delta['active_labels'] != 0 or delta['all_active_elastic_nodal_forces'] != 0
            or delta['crystal_energy'] != 0 or any(delta['currents'].values())):
        raise AssertionError(f'Frozen physical state changed: {delta}')
    if (abs(before['total_energy']-saved['fixed_field_total_energy'])>1.e-8
            or abs(before['max_free_force']-saved['final_free_label_force'])>1.e-12):
        raise AssertionError('Starting energy or force differs from reviewed endpoint')
    return dict(case=case,epsilon=opt.epsilon,state_sha256=digest(snapshot), active_nodes=int(opt.stencil.sum()),
                melt_only_nodes=int((~opt.stencil).sum()), free_active_dofs=len(old_grad),
                eliminated_dofs=sum(int(c.dead.sum()) for c in opt._components()), before=before,after=after,
                delta=delta, last_accepted_alpha_from_production=saved['final_alpha'])


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('scratch',type=Path)
    p.add_argument('output',type=Path)
    args=p.parse_args()
    rows={}
    for case in ('primary','secondary'):
        rows[case]=audit(args.scratch.resolve(),case)
        print(case, json.dumps(rows[case]['delta'],allow_nan=False),flush=True)
    with args.output.open('x') as f:
        json.dump(dict(origin='12e02d90',phase_drop=.4510084282857874,
                       phase='held R4, unchanged active trace, repeated harmonic solution',
                       cases=rows), f, indent=2, allow_nan=False)
        f.write('\n')
