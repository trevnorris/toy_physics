"""Check the two signed depth derivatives using promoted bath-order quadrature.

Each moving shape still uses its original stationary inner-label solution. Only
terms independent of those labels are replaced in its signed energy difference.
"""
import argparse
from copy import deepcopy
from pathlib import Path

import numpy as np

from mouth_r4_barrier_gradient import evidence
from mouth_exact_checkpoint_control import digest, read
from mouth_r4_stage2_outcome import verify_guard
from mouth_whole_r4_relax import atomic_json


def collect(area, durable=None):
    rawfile, hpfile = area/'GATE.json', area/'depth_precision.json'
    original, hp = read(rawfile), read(hpfile)
    verify_guard(area/'collection_guard', int(not original['passed']))
    verify_guard(area/'depth_precision_readout_guard', 0)
    if hp['schemas'] != ['order_gradient','double_well']:
        raise ValueError('Wrong depth quadrature columns')
    root = Path(__file__).resolve().parent
    for name, sha in original['sources'].items():
        if digest(root/name) != sha:
            raise RuntimeError(f'Source changed before re-read: {name}')
    gate = deepcopy(original)
    failures = [reason for reason in original['failed'] if not (
        len(reason) == 4 and reason[1] == 'log_depth' and reason[3] == 'central check')]
    gate.update(raw_gate_sha256=digest(rawfile),precision_measurement_sha256=digest(hpfile),
                raw_passed=original['passed'],raw_failed=original['failed'],
                corrections='Only the two extensive bath order energies are reevaluated in promoted precision, at both signs and both widths; no field or label was re-optimized.',
                guards=['depth_precision_readout_guard','precision_collection_guard','collection_guard'])
    gate['sources'].update({name:digest(root/name) for name in ('mouth_r4_depth_precision.py',Path(__file__).name)})
    tol = gate['absolute_tolerance']
    if not (tol==1.e-7 and hp['longdouble_epsilon']<1.e-18 and
            np.max(np.abs(np.subtract(hp['analytic'],hp['independent_80bit_analytic'])))<1.e-8 and
            max(map(abs,hp['energy_anchor_delta']))<2.e-9):
        raise RuntimeError('Cannot transfer the high-precision raw depth work to the anchored analytic evaluator')
    for case, group in gate['cases'].items():
        state, repair, analytic = evidence(area.parent,case)
        if digest(state)!=group['source_sha256'] or digest(repair)!=group['fixed_trace_sha256']:
            raise RuntimeError('Changed original or repaired source: '+case)
        for name, mode in group['modes'].items():
            work=mode['work']; r=area/(case+'_'+name)
            verify_guard(r.with_name(r.name+'_guard'),0)
            if mode != read(r/'readout.json'):
                raise RuntimeError('Original guard readout differs from the raw table')
            if name=='log_depth' and hp['analytic']!=analytic['modes'][name]['analytic'][3:5]:
                raise RuntimeError('The two bath contributions do not belong to this depth table')
            for row in work['steps']:
                banks=row['banks']
                if name=='log_depth':
                    hpwidth=hp['banks'][str(row['h'])]['signed']
                    row['raw_central']=row['central'];row['raw_central_minus_analytic']=row['central_minus_analytic']
                    for sign, bank in banks.items():
                        d=hpwidth[sign]
                        correction=float(np.sum(np.subtract(d['extended_work'],d['native_work']),dtype=np.longdouble))
                        bank['one_sided_raw']=bank['one_sided']
                        bank['order_work_correction']=correction
                        bank['one_sided']+=correction
                    row['central']=float(np.mean([b['one_sided'] for b in banks.values()],dtype=np.longdouble))
                    row['central_minus_analytic']=row['central']-work['analytic']
                if (set(banks)!={'-1','1'} or not all(b['valid'] and b['min_chart']>0 and
                        b['inner_residual']<=3e-8 for b in banks.values()) or
                        not abs(row['central_minus_analytic'])<=tol):
                    failures.append([case,name,row['h'],'signed reduced gradient'])
                if name!='active':
                    row['unit_gauge_central_minus_analytic']=float(np.mean([b['one_sided_unit_gauge_force'] for b in banks.values()])-mode['analytic_unit_gauge'])
        group['max_abs_centered_error']=max(abs(row['central_minus_analytic']) for mode in group['modes'].values() for row in mode['work']['steps'])
        group['max_abs_unit_gauge_centered_error']=max(abs(row['unit_gauge_central_minus_analytic']) for k,mode in group['modes'].items() if k!='active' for row in mode['work']['steps'])
    gate['passed']=not failures;gate['failed']=failures
    atomic_json(area/'PRECISION_GATE.json',gate)
    if durable is not None:
        atomic_json(durable,{'gate':gate,'depth_quadrature':hp})
    print('passed',gate['passed'],'failures',failures,flush=True)
    for key,item in gate['cases'].items():
        print(key,'errors',[(s['raw_central_minus_analytic'],s['central_minus_analytic']) for s in item['modes']['log_depth']['work']['steps']],flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('area',type=Path)
    p.add_argument('--durable',type=Path)
    args=p.parse_args()
    collect(args.area.resolve(),args.durable)
