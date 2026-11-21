"""Compare sign controls and the measured modulus-curvature requirement for the first bend."""
import json
from pathlib import Path
import numpy as np

from mouth_order_wave_pulse import digest, emit

ROOT=Path('_scratch/s11c')


def records(name):
    folder=ROOT/name
    guard=json.loads((folder/'outcome.json').read_text())
    limits=json.loads((folder/'limit-validation.json').read_text())
    if guard['exitCode'] or guard['childOutcome']['exitCode'] or not guard['limitsVerified'] or not limits['verified']:
        raise ValueError(f'{name}: incomplete guarded measurement')
    rows=[json.loads(line) for line in (folder/'stdout').read_text().splitlines()]
    return folder,rows


def read(name):
    folder,rows=records(name)
    trace=folder/'bend_trace.npz'
    if digest(trace)!=rows[-1]['trace_sha256']:
        raise ValueError(f'{name}: wave trace mismatch')
    with np.load(trace) as data:
        fields=dict(zip(data['columns'],data['readings'].T,strict=True))
    return fields,rows[-1]


def avg(fields,row,lo,hi):
    phase=(fields['t']-.5*row['dt']-lo)/(hi-lo)
    w=np.where((phase>=0)&(phase<=1),np.sin(np.pi*phase)**2,0.)
    return {key:float(w@fields[key]/w.sum()) for key in ('flare_force','metric_force','energy')}


def main():
    conditions=('pos','neg','half','quiet','finedt','sub2','ramp4','deep')
    before=('positive','negative','unbent','half','quiet','medium','finedt')
    rows={}
    for prefix,cases in (('ringing_after_',conditions),('ringing_bend_',before)):
        for case in cases:
            data,row=read(prefix+case+'_01')
            if data['t'][-1]<48:
                continue
            win=avg(data,row,24.,48.)
            rows[prefix+case]=data,row,win
            emit(test='finite_bend',case=prefix+case,energy=win['energy'],work=row['work'],
                 flare=win['flare_force'],entrance=-win['flare_force'],metric=win['metric_force'],
                 response_per_offset=win['flare_force']/row['bend'] if row['bend'] else None)
    positive,pr,pw=rows['ringing_after_pos']
    for case in ('neg','half','quiet'):
        data,r,win=rows['ringing_after_'+case]
        prior=positive['t']<=pr['start']
        common=np.flatnonzero(prior)
        # The opposite and half-wave launches compare squared amplitudes. Check the
        # two signed stress meters as well as energy on the common initial history.
        ratio=r['amplitude']/pr['amplitude']
        errs={k:float(np.max(abs(data[k][common]-ratio*positive[k][common]))) for k in ('collar_stress','deep_stress')}
        errs['energy']=float(np.max(abs(data['energy'][common]-ratio**2*positive['energy'][common])))
        expect=(r['bend']/pr['bend'])*ratio**2*pw['flare_force']
        emit(test='same_preparation',case=case,pre_bend_errors=errs,linear_force_prediction=expect,actual=win['flare_force'])
        if max(errs.values())>1e-13:
            raise ValueError('initial wave histories differ')
    for case,sub in (('pos',1),('sub2',2)):
        _,cr=records('ringing_bend_modulus_sub'+str(sub)+'_01')
        c=cr[-1]
        d,r,win=rows['ringing_after_'+case]
        if abs(r['work']-c['work'])>1e-9*c['work']:
            raise ValueError('incremental modulus test used a different straight excitation')
        for label,weight in c['stiffness_per_lambda'].items():
            lo,hi=map(float,label.split('-'))
            f=avg(d,r,lo,hi)['flare_force']/r['bend']
            emit(test='condition',sub=sub,window=label,flare_per_offset_at_constant_modulus=f,
                 restoring_per_lambda=weight,lambda_at_balance=f/weight,
                 wave_flare_per_offset_trial_lambda_point_one=f-.1*weight)
    # Compare the unbent edge field to the earlier axisymmetric dual-field
    # calculation on the same type of stroke. Its time step and remote bath differ.
    _,edge=read('ringing_bend_unbent_01')
    old=json.loads((Path('research/pde_ledger_v3/_measurements/mouth_formation_wall_work.json')).read_text())
    raw=old['guarded_runs']['formation_wall_ax2_01']['literal_stdout']
    reference=json.loads(raw.splitlines()[-1])['signed_wall_work']
    emit(test='launch_comparison',edge=edge['work'],dual=reference,
         relative_difference=edge['work']/reference-1)
    if abs(edge['work']/reference-1)>.06:
        raise ValueError('straight pulse differs from the archived wave')


if __name__=='__main__':
    main()
