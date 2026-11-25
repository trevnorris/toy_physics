"""Read accepted held R4 states and count recrossings by plane and observation site."""
import argparse
from collections import Counter
import json
from pathlib import Path
import sys

import numpy as np

from mouth_exact_checkpoint_control import digest, read
from mouth_r4_aligned_derivatives import load
from mouth_r4_eulerian_release import order_grid
from mouth_r4_segmented_held import CollarFaces
from mouth_r4_stage2_outcome import verify_guard


def measure(area):
    verify_guard(area/'primary_guard',1)
    old=area.parent/'r4_rematch39_01'
    primary=area/'primary'
    saved=sorted(primary.glob('step_*.json'))
    if not saved or saved[-1].stem!='step_00166' or (primary/'HELD.json').exists():
        raise ValueError('Expected the recorded fixed-shape stop at step 167')
    driver,trial,solid,order=load(old,'primary',area)
    monitor=CollarFaces(solid,driver.labels,order)
    term=read(saved[-1])['material']['worst_free_active']['node']
    target=int(np.argmin(np.abs(monitor.planes-1/12)))
    if monitor.planes[target]!=1/12: raise ValueError('The registered depth grid lacks the requested reference face')
    sources={}
    history={}
    rows=[]
    events=[]
    for step,recfile in enumerate(saved,109):
        paths=[recfile]
        matched=primary/f'matched_{step:05d}.json'
        if matched.exists(): paths.append(matched)
        for path in paths:
            rec=read(path)
            state=path.with_suffix('.npz')
            if step!=rec['outer'] or digest(state)!=rec['sha256']: raise ValueError(str(state))
            sources[path.name]=rec['sha256']
            with np.load(state,allow_pickle=False) as s:
                if int(s['radius'])!=4 or int(s['winding'])!=1 or not np.array_equal(s['x'],driver.x):
                    raise ValueError('Unexpected held state')
                labels=s['labels']
                order=order_grid(trial,s['field'])
            before=monitor.last
            after=monitor.banks(labels,order)
            peak=rec['material']['worst_free_active']['node']
            for key in list(history):
                bank,i,j,k,_=key
                if not after[bank][1][i,j,k]: del history[key]
            stage='fixed_order' if path==recfile else 'field_rematch'
            for bank,((ra,ma),(rb,mb)) in enumerate(zip(before,after)):
                for site in np.argwhere(ma & mb & (ra!=rb)):
                    node=tuple(site)
                    a,b=int(ra[node]),int(rb[node])
                    direction=1 if b>a else -1
                    for plane in range(min(a,b),max(a,b)):
                        key=(bank,*map(int,node),plane)
                        prev=history.get(key)
                        same=bool(prev and prev[0]==-direction)
                        history[key]=(direction,step,stage)
                        events.append(dict(step=step,stage=stage,kind='ray' if bank else 'node',
                            rule_bank=bank,site=site.tolist(),plane=plane,xi_t=float(monitor.planes[plane]),
                            direction=direction,recross=same,previous_step=prev[1] if same else None,
                            near_peak=bool(np.all(np.abs(site-peak)<=1)),near_last_peak=bool(np.all(np.abs(site-term)<=1))))
            native=monitor.observe(labels,order,peak if stage=='fixed_order' else None)
            if stage=='fixed_order' and native!=rec['reference_crossings']:
                raise ArithmeticError(f'Archived crossing detector differs at {step}')
    def counts(group):
        return dict(Counter((e['kind']+'_recross' if e['recross'] else e['kind']+'_first_cross') for e in group))
    def number(group):
        result=counts(group)
        result['total']=len(group)
        result['distinct_site_plane']=len({(e['rule_bank'],*e['site'],e['plane']) for e in group})
        return result
    for n in range(109,167):
        local=[e for e in events if e['step']==n]
        rows.append(dict(step=n,all=counts(local),near_peak_recross=sum(e['recross'] and e['near_peak'] for e in local),
                         one_twelfth_recross=sum(e['recross'] and e['plane']==target for e in local)))
    one=[e for e in events if e['plane']==target]
    reversals=[e for e in events if e['recross']]
    near=[e for e in reversals if e['near_peak']]
    last=[e for e in reversals if e['near_last_peak']]
    planes=Counter(e['plane'] for e in reversals)
    evidence=dict(sha256=digest(Path(__file__)),runner_sha256=digest(Path(__file__).with_name('mouth_r4_segmented_held.py')),
        source_step_108=digest(old/'primary/rematched.npz'),accepted_states=sources)
    return dict(command=sys.orig_argv,evidence=evidence,definition='A recrossing is opposite signed traversal of the same reference depth by the same continuously tracked ordered collar node or quadrature-ray endpoint. Near means at most one material-grid index away in each direction. Between saved endpoints only net changes are identifiable.',
                all=number(events),per_plane=[dict(xi_t=float(monitor.planes[i]),recross=n) for i,n in sorted(planes.items())],
                during_field_rematch=number([e for e in events if e['stage']=='field_rematch']),
                one_twelfth=number(one),recross_near_instantaneous_peak=number(near),
                recross_near_last_peak=number(last),one_twelfth_recross_near_peak=[e for e in near if e['plane']==target],
                last_reversal=reversals[-1] if reversals else None,
                last_cross=events[-1] if events else None,
                one_twelfth_last_reversal=next((e for e in reversed(one) if e['recross']),None),
                distinct_recrossed_sites=len({(e['rule_bank'],*e['site']) for e in reversals}),
                final_peak_node=term,counts_by_step=rows)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('study',type=Path)
    parser.add_argument('output',type=Path)
    args=parser.parse_args()
    report=measure(args.study.resolve())
    with args.output.open('x') as out: json.dump(report,out,indent=2,allow_nan=False);out.write('\n')
    print(json.dumps({k:report[k] for k in ('all','one_twelfth','recross_near_instantaneous_peak',
            'recross_near_last_peak','last_reversal','last_cross','one_twelfth_last_reversal')},allow_nan=False),flush=True)
