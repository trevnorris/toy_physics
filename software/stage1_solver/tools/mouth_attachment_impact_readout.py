"""Reduce guarded impact trajectories using the axial or face arrival denominator."""
import argparse
import json
from pathlib import Path

import numpy as np

from mouth_attachment_impact_trial import sha
from mouth_curved_site_pair import validated
from neutral_patch_sites import emit


def read(folder):
    prep,archive=validated(folder)
    path=Path(archive['path'])
    if sha(path)!=archive['sha256'] or sha(Path(prep['section']))!=prep['section_sha256']:
        raise ValueError('Impact input or output changed')
    rows=[json.loads(line) for line in path.read_text().splitlines()]
    if len(rows)!=prep['trajectories'] or len(rows)!=archive['traces']:
        raise ValueError('Missing launches')
    return prep,rows


def keys(row):
    r=row['result']
    before_exit=r['first_joint_time']>=0 and not r['gate_before_binding']
    censored=r['first_joint_time']<0 and r['gate_exit_time']<0
    return (before_exit,censored,r['first_joint_time']>=0,r['joint_end'],r['joint_both'],
            r['residence15'],r['joint_residence'],r['gated_joint_residence'],
            r['energy_tag_final']-r['energy_tag_initial'],r['host_kinetic'],r['host_strain'],
            r['work_to_host'],r['vz0']-r['vz_end'])


COLS=('captured','censored','ever','at_end','both_times','near_time','joint_time','gate_time',
      'tracer_energy_change','host_kinetic','host_strain','host_work','normal_impulse')


def mean(rows):
    weights=np.array([r['weight'] for r in rows])
    data=np.array([keys(r) for r in rows])
    return dict(zip(COLS,np.average(data,axis=0,weights=weights),strict=True))


def report(name,rows,axial):
    avg=mean(rows)
    diag=[r['result'] for r in rows]
    capture=np.array([keys(r)[0] for r in rows])
    times=[r['first_joint_time'] for r in diag if r['first_joint_time']>=0]
    emit(test='zero_flow',group=name,arrivals=len(rows),
         axial_gate_applies=axial,early_count=int(sum(capture)) if axial else None,
         ever_count=sum(r['first_joint_time']>=0 for r in diag),
         end_count=sum(r['joint_end'] for r in diag),
         half_and_end_count=sum(r['joint_both'] for r in diag),
         unresolved_gate_count=int(sum(keys(r)[1] for r in rows)) if axial else None,
         initial_energy_min=min(r['energy_tag_initial'] for r in diag),first_joint_time_range=
         [min(times),max(times)] if times else None,energy_drift_max=max(r['max_energy_error'] for r in diag),
         work_error_max=max(r['max_work_error'] for r in diag),
         tangent_momentum_error_max=max(r['max_momentum_error'] for r in diag),
         plane_crossing_count=sum(r['plane_crossing_time']>=0 for r in diag),
         weighted_means={k:v for k,v in avg.items() if axial or k not in ('captured','censored','gate_time')})
    if max(r['max_energy_error'] for r in diag)>1.e-3 or max(r['max_momentum_error'] for r in diag)>1.e-8:
        raise ArithmeticError(f'Unresolved trajectory budget in {name}')
    if any(r['plane_crossing_time']>=0 for r in diag):
        raise ArithmeticError('Trajectory crossed the supported plane')


def signed(name,rows):
    for speed in (.002,.001):
        for idx,quantity in [(0,'early_capture'),(7,'truncated_gate_residence')]:
            moments=[]
            groups=[]
            for U in (speed,-speed):
                subset=[r for r in rows if r['U']==U]
                if not subset: return
                groups.append(subset)
                values=np.array([r['weight']*keys(r)[idx] for r in subset])
                angles=np.array([r['cosine'] for r in subset])
                moments.append(np.dot(angles,values)/values.sum() if values.sum()>0 else None)
            if [v['sample'] for v in groups[0]]!=[v['sample'] for v in groups[1]]:
                raise ValueError('Unpaired signed inputs')
            changed=sum(keys(a)[0]!=keys(b)[0] for a,b in zip(*groups,strict=True))
            emit(test='signed',group=name,quantity=quantity,U=speed,
                 angular_means=moments,normalized_signed_slope=
                 None if None in moments else 3*(moments[0]-moments[1])/(2*speed),
                 paired_outcome_switches=changed,paired_comparisons=len(groups[0]))


def run(root,sets):
    all_zero={}
    inputs={}
    for name in sets:
        prep,rows=read(root/name)
        inputs[name]=sha(root/name/'stdout')
        zero=[r for r in rows if r['U']==0]
        all_zero[name]=zero
        report(name,zero,prep['preparation']!='beam')
        signed(name,rows)
    for name in ('ax15_l0','ax15_l1','ax20_l0','ax20_l1','beam_l1'):
        labels=[f'{name}_b{i}' for i in range(4)]
        if all(l in all_zero for l in labels):
            group=[r for l in labels for r in all_zero[l]]
            report(name+'_pooled',group,not name.startswith('beam'))
            emit(test='scrambles',group=name,gate_capture_by_scramble=
                 [mean(all_zero[l])['captured'] for l in labels] if not name.startswith('beam') else None,
                 end_by_scramble=[mean(all_zero[l])['at_end'] for l in labels])
    for primary,control in [('ax15_l1_b0','ax15_fine'),('ax15_l1_b0','ax15_large'),
                            ('ax15_l1_b0','ax15_rigid_control'),('ax15_l1_b0','ax15_long'),
                            ('beam_l1_b0','beam_fine'),('beam_l1_b0','beam_large'),
                            ('beam_l1_b0','beam_rigid_control'),('beam_l1_b0','beam_long')]:
        if primary not in all_zero or control not in all_zero:continue
        x,y=all_zero[primary],all_zero[control]
        if [r['initial'] for r in x]!=[r['initial'] for r in y]:
            raise ValueError('Control changed launch states')
        emit(test='control',primary=primary,control=control,rows=len(x),
             early_state_disagreements=sum(keys(a)[0]!=keys(b)[0] for a,b in zip(x,y,strict=True)),
             end_state_disagreements=sum(keys(a)[3]!=keys(b)[3] for a,b in zip(x,y,strict=True)),
             max_abs_gate_time_difference=max(abs(keys(a)[7]-keys(b)[7]) for a,b in zip(x,y,strict=True)))
    emit(test='source_outputs',sha256=inputs)


def durations(root,sets):
    """Report observed bound time, separately for early and late acquisitions.

    Quantiles are weighted inverse empirical CDFs. Bound time can consist of
    several visits; a bound last sample is a censored observation of departure.
    """
    for label in sets:
        names=[f'{label}_b{i}' for i in range(4)] if not (root/label).is_dir() else [label]
        rows=[]
        for name in names:
            prep,batch=read(root/name)
            rows.extend(r for r in batch if r['U']==0)
        clock=prep['model_time']/prep['native_time']
        for kind in ('ever','in_slice'):
            if kind=='in_slice' and prep['preparation']=='beam':
                continue
            selected=[r for r in rows if (r['result']['first_joint_time']>=0 if kind=='ever' else keys(r)[0])]
            if not selected: continue
            weight=np.array([r['weight'] for r in selected])
            results=[r['result'] for r in selected]
            measured={}
            for key in ('first_joint_time','joint_residence','gated_joint_residence','energy_tag_final',
                        'work_to_host','host_kinetic','host_strain'):
                if key=='gated_joint_residence' and prep['preparation']=='beam':continue
                vals=np.array([r[key] for r in results],dtype=float)
                if key in ('first_joint_time','joint_residence','gated_joint_residence'):
                    vals*=clock
                ix=np.argsort(vals)
                quantiles=vals[ix][np.searchsorted(np.cumsum(weight[ix])/sum(weight),[.1,.5,.9])]
                measured[key]=dict(mean=np.average(vals,weights=weight),q10_q50_q90=quantiles.tolist())
            emit(test='observed_duration',group=label,selection=kind,acquisitions=len(results),
                 still_bound_at_end=sum(r['joint_end'] for r in results),
                 time_unit='mouth model time',energy_unit='native well energy per tracer mass',
                 horizon=prep['model_time'],measurements=measured)
    for early,later in [('ax15_l1_b0','ax15_long'),('beam_l1_b0','beam_long')]:
        if later not in sets: continue
        _,short=read(root/early)
        _,long=read(root/later)
        short=[r for r in short if r['U']==0]
        long=[r for r in long if r['U']==0]
        assert [r['initial'] for r in short]==[r['initial'] for r in long]
        pairs=np.array([[r['result']['joint_end'],s['result']['joint_end']] for r,s in zip(short,long,strict=True)])
        emit(test='end_states',early=early,late=later,neither=int(sum(~pairs.any(axis=1))),
             at_both=int(sum(pairs.all(axis=1))),early_only=int(sum(pairs[:,0]&~pairs[:,1])),
             late_only=int(sum(~pairs[:,0]&pairs[:,1])))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('root',type=Path)
    parser.add_argument('runs',nargs='+')
    parser.add_argument('--durations',action='store_true')
    args=parser.parse_args()
    (durations if args.durations else run)(args.root,args.runs)
