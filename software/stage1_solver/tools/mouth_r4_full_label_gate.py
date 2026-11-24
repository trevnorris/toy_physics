"""Combine untouched broad derivatives and the preselected narrow held-out coordinates."""
import argparse
from pathlib import Path

import numpy as np

from mouth_exact_checkpoint_control import digest, read
from mouth_r4_full_label_work import STATE_SHA
from mouth_r4_stage2_outcome import verify_guard
from mouth_whole_r4_relax import atomic_json


def run(area,measurements):
    for name in ('admitted_guard','precision_guard','quad_guard'):
        verify_guard(area/name,0)
    names=['mouth_r4_full_label_'+s+'.json' for s in ('work','precision','quad')]
    raw,hp,quad=[read(measurements/n) for n in names]
    if raw['source_sha']!=STATE_SHA or hp['state_sha']!=quad['state_sha'] or quad['state_sha']!=STATE_SHA:
        raise ValueError('State pins differ')
    if hp['raw_measurement_sha']!=digest(measurements/names[0]) or quad['raw_measurement_sha']!=hp['raw_measurement_sha']:
        raise ValueError('Different full trial')
    widths=2.**-np.arange(25,39,2)
    output={};counts=[0,0]
    for epsilon in (1e-8,1e-10):
        peaks=[];melted=[]
        for c,row in enumerate(raw['components']):
            with np.load(area/'result'/f'coordinate_{c}.npz') as f:
                free=f['free'];rich=f['richardson'];p=f['analytic_crystal'];gg=f['analytic_gauge']
                error=(rich[0]-p)+epsilon*(rich[1]-gg)
                active_count=row['active_count'];dead_count=row['melted_count']
                if len(free)!=active_count+dead_count:raise ValueError('Incomplete broad trial')
                if epsilon==1e-8:
                    counts[0]+=active_count;counts[1]+=dead_count
                flagged=free[abs((rich[0,free]-p[free])+1e-8*(rich[1,free]-gg[free]))>5e-11]
                if len(flagged)!=quad['components'][c]['examined']:raise ValueError('Wrong narrow holdout')
                if len(flagged):
                    if np.any(gg[flagged]!=0):raise ValueError('Nonzero gauge held out')
                    with np.load(area/'quad'/f'coordinate_{c}.npz') as q:
                        if not np.array_equal(q['group'],flagged) or not np.array_equal(q['widths'],widths):raise ValueError('Unexpected quad bank')
                        # Require both final raw centered widths and their mutual agreement;
                        # additionally keep each of the seven independent raw widths in evidence.
                        d=q['signed'].mean(axis=1)-p[flagged]
                        if max(abs(d[-2:]).ravel())>1e-10 or max(abs(d[-1]-d[-2]))>1e-10:
                            raise ArithmeticError('Promoted signed derivatives unresolved')
                        error[flagged]=d[-1]
                peaks.append(float(max(abs(error[free]))))
                melted.append(row[str(epsilon)]['fine_richardson']['melted'])
        output[str(epsilon)]=dict(max_abs_combined=max(peaks),components=peaks,max_abs_melt=max(melted),passed=max(peaks)<=1e-10)
    result=dict(source_sha256=STATE_SHA,checksums={n:digest(measurements/n) for n in names},
        accepted_dof_counts=dict(active=counts[0],melted=counts[1],all=sum(counts)),
        binary128_held_out_count=sum(row['examined'] for row in quad['components']),
        complete_minimum_certified_chart_volume=raw['signed_banks_min_chart_orientation_bound'],
        binary128_promoted_raw_max_errors=quad['components'][0]['raw_max_abs_errors'],
        width_rule='All registry-free coordinates used coarse, half and quarter centered signed pairs. Marginal or failing broad extrapolations (>5e-11) used disjoint new 80-bit widths and seven new 128-bit widths. Both final raw binary128 widths and their mutual agreement must pass independently of any intermediate 80-bit width or post-hoc choice.',
        provenance='Binary128 widens all original inputs and subtracts cell anchor energy before narrowing. Its 120 bisections resolve reference plane roots below binary128 machine precision. The three narrowest raw 80-bit widths lose precision with 52-bisection roots. All componentwise wide and narrow signed results are retained.',
        results=output,passed=all(v['passed'] for v in output.values()))
    atomic_json(measurements/'mouth_r4_full_label_gate.json',result)
    return result['passed']


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('area',type=Path);parser.add_argument('measurements',type=Path)
    args=parser.parse_args()
    raise SystemExit(0 if run(args.area,args.measurements) else 1)
