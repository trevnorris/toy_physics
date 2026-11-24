"""Test the opt-in precise inner solver with the mouth and every active label fixed."""
import argparse
from hashlib import sha256
import importlib.util
import json
from pathlib import Path
import sys

import numpy as np


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def run(area, direction, audit, work, output):
    pin = area/'pin/software/stage1_solver/tools'
    sys.path.insert(0,str(pin))
    from mouth_cp10_label_grid import get_cp10
    from mouth_r4_eulerian_relax_labels import LabelRelaxation, chart_orientation
    from mouth_r4_eulerian_solid import SpatialCrystal
    from mouth_r4_stage2_outcome import verify_guard

    for directory in ('precision_guard','direction_guard'):
        verify_guard(direction.parent.parent/directory,0)
    verify_guard(area/'primary_guard',1)
    fixed = json.loads(audit.read_bytes())
    fork = Path(__file__).with_name('mouth_r4_melt_barrier.py')
    spec = importlib.util.spec_from_file_location('audited_inner',fork)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    if digest(direction)!=fixed['historical']['direction_sha256']:
        raise ValueError('Wrong historical direction')
    state = area/'primary/failed_trial.npz'
    if digest(state)!=fixed['inner_state_sha256']:
        raise ValueError('Wrong archived iterate')
    b,_,axes = get_cp10(area.parent)
    with np.load(state,allow_pickle=False) as f:
        labels = f['labels']
        order = np.ascontiguousarray(f['field'][2,b.ids].transpose(0,2,1))
    opt = LabelRelaxation(b,labels,order,1e-8)
    barrier = module.MeltBarrier(b,opt,work)
    with np.load(direction,allow_pickle=False) as f:
        step = f['step']
    alpha = fixed['historical']['first_extended_descent']['alpha']
    delta = float(barrier.local_work(labels,labels+alpha*step))
    derivative = fixed['historical']['slope']
    signed_reference = fixed['historical']['first_extended_descent']['extended_rounded']
    if abs(delta-signed_reference)>1e-23 or not delta<1e-4*alpha*derivative:
        raise ArithmeticError('Native precise work disagrees with the independent measurement')
    loose, early = barrier.solve(labels,'far_outer',work,outer_force=0.07797064397780097)
    if not np.array_equal(loose,labels) or early['accepted_steps'] or not early['stationary']:
        raise ArithmeticError('Adaptive far-outer target did not admit the starting melt field')
    polished, result = barrier.solve(labels,'terminal_outer',work,outer_force=.000030)
    journal = [json.loads(v) for v in (work/'terminal_outer.jsonl').read_bytes().splitlines()]
    increments = [v['last'] for v in journal[1:]]
    solid = SpatialCrystal(b)
    e0,g0,_,_ = solid.evaluate(labels,order)
    e1,g1,_,_ = solid.evaluate(polished,order)
    _,gr,parts,_ = barrier.sample(polished)
    metrics = dict(max_trace_motion=float(np.abs((polished-labels)[opt.stencil]).max()),
                   crystal_energy_delta=float(e1-e0),
                   active_physical_force_max_difference=float(np.abs((g1-g0)[opt.stencil]).max()),
                   direct_residual=float(np.abs(barrier.extract(gr,'dead')).max()),
                   chart=chart_orientation(axes,polished,opt.active),
                   initial_to_polished_unit_work=float(barrier.local_work(labels,polished)),
                   final_corner_min=parts['corner_min'])
    if metrics['max_trace_motion'] or metrics['crystal_energy_delta'] or metrics['active_physical_force_max_difference']:
        raise ArithmeticError('Fixed-trace physical state changed')
    if result['residual']!=metrics['direct_residual'] or any(v['delta']>=0 or v['chart_lower']<=0 for v in increments):
        raise ArithmeticError('Inner-step readout disagrees')
    np.savez_compressed(work/'polished_labels.npz',labels=polished)
    record=dict(command=sys.orig_argv,input_sha256={str(p):digest(p) for p in (state,direction,audit)},
        source_sha256={p.name:digest(p) for p in (Path(__file__),fork,fork.with_suffix('.cpp'))},
        authority='Diagnostic only, fixed order and active labels; output is not a resumed held state',
        reference_alpha=alpha,work_native=delta,work_independent=signed_reference,
        work_absolute_discrepancy=abs(delta-signed_reference),analytic_slope=derivative,
        loose_scale=.07797064397780097,loose=early,terminal_scale=.000030,terminal=result,checks=metrics,
        journal=journal,polished_sha256=digest(work/'polished_labels.npz'))
    with output.open('x') as f:
        json.dump(record,f,indent=2,allow_nan=False)
        f.write('\n')
    print(json.dumps({k:v for k,v in record.items() if k not in ('journal','command','input_sha256','source_sha256')}),flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    for key in ('area','direction','audit','work','output'):
        parser.add_argument(key,type=Path)
    args=parser.parse_args()
    args.work.mkdir(exist_ok=False)
    run(args.area.resolve(),args.direction.resolve(),args.audit.resolve(),args.work.resolve(),args.output)
