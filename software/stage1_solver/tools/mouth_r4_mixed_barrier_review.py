"""Verify the signed reduced work of the mixed three-component R4 label probe."""
import argparse
import json
from math import fsum, isfinite, sqrt
from pathlib import Path

import numpy as np

from mouth_r4_barrier_work_review import digest, guard, read


WORKER_SHA='f5306f7197a4a97b755188083f49c50e1fcf7b316920cc9fcb01e131c629ccf3'


def check(area):
    source=Path(__file__).resolve().with_name('mouth_r4_mixed_barrier_work.py')
    bankmap=read(area/'BANKS.json')
    if digest(source)!=WORKER_SHA or bankmap['worker_source_sha256']!=WORKER_SHA:
        raise ValueError('The executed direction differs from the independently reviewed source')
    rows={}; guards={}
    for name,eps in [('primary',1.e-8),('secondary',1.e-10)]:
        run=area/name
        log=area/(name+'_guard')
        guard(log)
        guards[name]=digest(log/'outcome.json')
        raw=read(run/'readout.json')
        source_npz=area.parent/'r4_eulerian_stage2_01'/name/'label_stall.npz'
        repair=area.parent/'r4_slaved_melt_barrier01'/name/'fixed_active.npz'
        deposit=np.load(run/'packed_direction.npy',allow_pickle=False)
        refs=bankmap['cases'][name]
        if not (refs['guard_sha256']==guards[name] and refs['readout_sha256']==digest(run/'readout.json')
                and refs['source_sha256']==raw['source_sha256']==digest(source_npz)
                and refs['fixed_trace_sha256']==raw['fixed_trace_sha256']==digest(repair)
                and raw['packed_direction_sha256']==digest(run/'packed_direction.npy')
                and raw['source_test_sha256']==WORKER_SHA and raw['epsilon']==eps and raw['case']==name
                and raw['seed']==9256324 and raw['baseline_inner']<=3e-8 and raw['passed']):
            raise ValueError(f'Unlinked held state or failed input in {name}')
        components=raw['component_work']
        if [c['component'] for c in components]!=['radial','axial','angular']:
            raise ValueError('Incomplete virtual displacement')
        i=0
        for c in components:
            v=deposit[i:i+c['free']]; i+=c['free']
            n=np.count_nonzero(v)
            if not(n==c['moved'] and n>.5*len(v) and c['cross_coupled_moved']>0
                   and np.isclose(sqrt(fsum(float(x)**2 for x in v)),c['l2'],rtol=1e-14)
                   and np.max(abs(v))==c['max_abs']):
                raise ValueError(f'Packed group direction disagrees with {name} component record')
        if i!=len(deposit):
            raise ValueError('Packed direction length mismatch')
        analytic=fsum(c['analytic_physical']+eps*c['analytic_unit_gauge'] for c in components)
        if max(abs(analytic-raw['analytic_from_nodes']),abs(analytic-raw['work']['analytic']))>1e-12:
            raise ValueError('Nonconjugate packed and nodal derivatives')
        errors=[]; bounds=[]; station=[]; gains=[]; active=[]
        for j,step in enumerate(raw['work']['steps']):
            h=step['h']
            if j>1 or h!=(2.e-6,1.e-6)[j]:
                raise ValueError('Unexpected mixed probe width')
            banks=step['banks']
            if set(banks)!={'-1','1'}:
                raise ValueError('Missing physical sign')
            for sign,p in banks.items():
                state=run/(('plus' if sign=='1' else 'minus')+('_half' if j else '_full')+'.npz')
                if not (digest(state)==p['sha256'] and p['valid'] and p['iterations']>0
                        and p['min_Jinv']>0 and p['min_Jinv_active']>0 and p['min_chart']>0
                        and p['inner_residual']<=3e-8 and p['melted_label_change']>0):
                    raise ValueError('Unverified signed minimizer')
                bounds.append(p['min_chart']);station.append(p['inner_residual'])
                gains.append(p['melted_label_change']);active.append(p['min_Jinv_active'])
            eplus,eminus=banks['1'],banks['-1']
            fd=fsum([eplus['elastic']-eminus['elastic'],eps*(eplus['gauge_unit']-eminus['gauge_unit'])])/(2*h)
            error=fd-analytic
            if not(isfinite(error) and abs(error)<=1e-7 and abs(error-step['central_minus_analytic'])<=1e-12):
                raise ArithmeticError(f'Mixed endpoint energies do not reproduce the reduced derivative: {name} {h} {error}')
            errors.append(error)
        if len(errors)!=2:
            raise ValueError('Missing full or half step')
        rows[name]=dict(analytic=analytic,full_half_errors=errors,min_chart=min(bounds),
                        min_active_physical_Jinv=min(active),max_inner_residual=max(station),
                        min_melt_adjustment=min(gains),moved_groups=[c['moved'] for c in components],
                        coupled_groups=[c['cross_coupled_moved'] for c in components])
    return dict(passed=True,banks_sha256=digest(area/'BANKS.json'),worker_source_sha256=WORKER_SHA,
                reviewer_source_sha256=digest(Path(__file__)),cases=rows,guards=guards,
                independent_review='All three coordinate components participate in the same admissible active-label variation. At both widths and both gauge weights, the four signed stationary endpoint energies reproduce the analytic directional derivative within 1e-7.',
                limitations='One reproducible dense linear combination checks registry, seam and aggregate graded work. It does not establish componentwise agreement for every basis vector or stationarity at the held shape.')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('area',type=Path)
    parser.add_argument('destination',type=Path)
    args=parser.parse_args()
    result=check(args.area)
    with args.destination.open('x') as stream:
        json.dump(result,stream,allow_nan=False,indent=2)
        stream.write('\n')
    print(json.dumps(result))
