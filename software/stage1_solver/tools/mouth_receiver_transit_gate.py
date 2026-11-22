"""Test persistence and paired shear/compression during one conditional receiver transit."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from mouth_receiver_material_screen import signature
from mouth_receiver_pressure_readout import fraction, interval, weighted


def read(path):
    return json.loads(path.read_bytes())


def phase(state,ordered):
    blocks=state['blocks']; late=blocks[-8:]; assert len(blocks)==16 and all('state' in b for b in blocks)
    level=weighted(late[-4:],'density'); previous=weighted(late[:4],'density')
    dense=all(b['mean']['density']>=1.2 and b['structure']['max_cluster']>=116 for b in late)
    prepared=signature(state['prepared'])
    # A fully dispersed configuration has already failed density. Avoid a divergent reciprocal-grid count.
    observed=[signature(b['state']) if b['end_reading']['density']>=.3 else None for b in late]
    def expected(s):
        return s is not None and (min(s['four'])>=.2 if ordered else s['peak']<=.2)
    end=blocks[-1]['structure'];lost=end['absent']/end['original'] if end['original'] else None
    mobile=lost is not None and lost>=.5 and end['nonaffine_rms']>=1.
    accepted=(dense and abs(level/previous-1)<=.1 and expected(prepared)
              and all(expected(s) for s in observed) and (mobile or ordered))
    return dict(pass_state=bool(accepted),ordered_control=ordered,
        mean_densities=[float(b['mean']['density']) for b in late],min_late_cluster=min(b['structure']['max_cluster'] for b in late),
        drift_q4_over_q3=float(level/previous-1),pressure_q4=weighted(blocks[-4:],'pressure'),T_q4=weighted(blocks[-4:],'temperature'),
        absent_fraction=lost,rms=end['nonaffine_rms'],prepared_four=prepared['four'],prepared_peak=prepared['peak'],
        late_peaks=[s['peak'] if s is not None else None for s in observed],late_four=[s['four'] if s is not None else None for s in observed])


def gate(root,stage):
    plan=read(root/'design.json'); run_index=read(root/('ZERO_RUNS.json' if stage=='zero' else 'SIGNED_RUNS.json'))['results']
    states={};conditions={}
    for k,c in plan['cases'].items():
        if stage=='zero' and c['label']!='zero':continue
        raw=(root/'runs'/(k+'.json')).read_bytes()
        if hashlib.sha256(raw).hexdigest()!=run_index[k]: raise ValueError(f'Changed data: {k}')
        r=json.loads(raw); h=c['id']; label=c['label']; ordered=c['ordered']
        assert (r['history'],r['seed'],r['strain'],r['pressure_external'],r['temperature'],r['axis_index'],r['warm'],r['measure'])==(h,c['seed'],c['strain'],c['pressure'],1.088,1,0.,plan['horizon_native'])
        assert r['source_temperature']==1.088 and r['blocks'][-1]['state']==r['final']
        conditions[k]=phase(r,ordered);states.setdefault(h,{})[label]=r
    assert len(run_index)==len(conditions)
    states_ok=all(x['pass_state'] for x in conditions.values())
    record=dict(stage=stage,state_gate=states_ok,conditions=conditions)
    if stage=='signed':
        response={}
        for h,arms in states.items():
            z=arms['zero']; assert all(r['prepared']==z['prepared'] for r in arms.values())
            def signed(key):
                p,m=arms[key+'p'],arms[key+'m']; den=p['strain']-m['strain'];
                early=(p['initial_reading']['shear']-m['initial_reading']['shear'])/den
                gp=weighted(p['blocks'][-4:],'shear');gm=weighted(m['blocks'][-4:],'shear');gz=weighted(z['blocks'][-4:],'shear')
                return dict(prompt=early,late=(gp-gm)/den,odd_late_over_prompt=(gp-gm)/(den*early),even_scaled=(gp+gm-2*gz)/den)
            v={'full':signed('g'),'half':signed('h')};v['full_minus_half_late']=v['full']['late']-v['half']['late']
            if 'pp' in arms:
                p,m=arms['pp'],arms['pm'];den=p['pressure_external']-m['pressure_external']
                v['compliance']=(weighted(m['blocks'][-4:],'log_volume')-weighted(p['blocks'][-4:],'log_volume'))/den
                v['realized_pressure_difference']=weighted(p['blocks'][-4:],'pressure')-weighted(m['blocks'][-4:],'pressure')
            response[h]=v
        train=[response[h] for h in (0,1,2)]
        ratios={k:fraction([x[k]['prompt'] for x in train],[x[k]['late'] for x in train]) for k in ('full','half')}
        small=all(b is not None and b[0]>=-.25 and b[1]<=.25 for b in [r['bounded_ratio_t95'] for r in ratios.values()])
        compression=interval([x['compliance'] for x in train]);pressure=interval([x['realized_pressure_difference'] for x in train]);diff=interval([x['full_minus_half_late'] for x in train]);zero_P=interval([conditions[f'id{h}_zero']['pressure_q4'] for h in (0,1,2)])
        matched_pressure=abs(zero_P['mean'])<=.25 and pressure['t95'][0]>0 and abs(pressure['mean']-.5)<=.25
        held=response[4]; held_ok=all(abs(held[k]['odd_late_over_prompt'])<=.25 and held[k]['prompt']>0 for k in ('full','half')) and held['compliance']>0 and held['realized_pressure_difference']>0
        ordered_ok=all(response[h][k]['prompt']>0 and response[h][k]['late']>0 for h in (40,42) for k in ('full','half'))
        record.update(response=response,seed87_ratio=ratios,seed87_compliance=compression,seed87_pressure_diff=pressure,seed87_zero_pressure=zero_P,seed87_full_minus_half_late=diff,
            pressure_resolved=matched_pressure,seed41_single_path_meets=held_ok,ordered_control_resists=ordered_ok,
            eligible_for_interface=bool(states_ok and small and compression['t95'][0]>0 and matched_pressure and held_ok and ordered_ok and diff['t95'][0]<=0<=diff['t95'][1]))
    return record


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('root',type=Path);p.add_argument('stage',choices=('zero','signed'));args=p.parse_args();r=gate(args.root,args.stage)
    out=args.root/('ZERO_GATE.json' if args.stage=='zero' else 'SIGNED_GATE.json');out.write_text(json.dumps(r,indent=2,allow_nan=False)+'\n')
    for n,c in r['conditions'].items(): print(json.dumps(dict(case=n,pass_state=c['pass_state'],min_late_cluster=c['min_late_cluster'],absent_fraction=c['absent_fraction'],density_end=c['mean_densities'][-1],final_peak=c['late_peaks'][-1])))
    if args.stage=='signed': print(json.dumps({k:v for k,v in r.items() if k!='conditions'},allow_nan=False))
