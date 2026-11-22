"""Validate the guarded coupled order and projected time readouts without replaying trajectories."""
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
BASE = ROOT / '_scratch/s11c/live_waist_dynamics01'


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def record(name):
    p=BASE/name
    read=lambda s: json.loads((p/s).read_text())
    inv,lim,duration,child=map(read,('invocation.json','limit-validation.json',
                                     'duration-validation.json','child-outcome.json'))
    adopted = not (p/'outcome.json').exists()
    if adopted:
        if name != 'r1_tight_48_support':
            raise ValueError(f'Missing native guard outcome for {name}')
        receipt=read('terminal_worker_verification.json')
        if not (receipt['worker_outcome']==child and receipt['worker_exit']==0
                and receipt['guard_exit'] is None and not receipt['guard_ended_normally']
                and receipt['adopted_authorized_worker'] and receipt['limits_verified']
                and receipt['unit']==inv['unit']+'.service'
                and all(digest(p/f)==h for f,h in receipt['evidence_sha256'].items())):
            raise ValueError('Incomplete original worker terminal verification')
        endfile=p/'child-outcome.json'
        guard_exit=None
        wall_seconds=child['wallSeconds']
    else:
        endfile=p/'outcome.json'
        out=read('outcome.json')
        guard_exit=out['exitCode']
        wall_seconds=out['wallSeconds']
        if guard_exit!=child['exitCode'] or out['childOutcome']!=child or not out['limitsVerified']:
            raise ValueError(f'Native guard outcome disagrees with worker in {name}')
    if not (lim['verified'] and duration['verified'] and inv['wallDeadlineSeconds'] is None
            and inv['swapMax']==0 and lim['actual']['memory.swap.max']=='0'
            and duration['actual']['RuntimeMaxUSec']=='infinity'):
        raise ValueError(f'Guard containment does not close for {name}')
    samples=list(map(json.loads,(p/'resource-samples.jsonl').read_text().splitlines()))
    if not samples or any(q['memory.swap.current']!='0' for q in samples):
        raise ValueError(f'Invalid memory samples for {name}')
    f=p/'stdout'
    raw=f.read_text()
    return dict(command=inv['command'],cwd=inv['cwd'], started_utc=inv['startedUtc'],
                finished_record_utc=datetime.fromtimestamp(endfile.stat().st_mtime,timezone.utc).isoformat(),
                duration=wall_seconds,guard_exit=guard_exit,worker_exit=child['exitCode'],
                guard_ended_normally=not adopted,adopted_authorized_worker=adopted,
                source_of_terminal_metadata=str(endfile.relative_to(ROOT)),
                adoption_evidence=receipt if adopted else None,
                limits_verified=True,memory_cap=inv['memoryMax'],swap=0,wall_deadline=None,
                peak_memory_sample_bytes=max(int(a['memory.peak']) for a in samples),
                stdout=str(f.relative_to(ROOT)),stdout_bytes=f.stat().st_size,stdout_sha256=digest(f),
                literal_stdout=raw, stderr_tail=(p/'stderr').read_text().splitlines()[-2:])


def main():
    names=['pilot']
    resolved = {}
    for r in (1,2):
        for point in ('broad','mid','tight'):
            key=f'r{r}_{point}'
            special = r==1 and point!='broad'
            for n in (48,72):
                if special and n==48:
                    names.extend((key+'_48',key+'_48_bound'))
                case=f'{key}_{n}'+('_support' if special else '')
                names.append(case)
                resolved[case]=(r,point,n)
            names.append(f'time{r}_{point}_48'+('_final' if special else ''))
            if (r,point) in ((1,'broad'),(2,'tight')):
                names.append(f'time{r}_{point}_72')
    runs={name:record(name) for name in names}
    fields,summary={},{}
    for case,(r,point,n) in resolved.items():
        out=runs[case]
        rows=list(map(json.loads,out['literal_stdout'].splitlines()))
        accepted=[q for q in rows if q['kind']=='coupled']
        entry=dict(last_returned_stationary_action=max((q['action'] for q in accepted),default=None))
        if out['worker_exit']==0 and (out['guard_ended_normally'] or out['adopted_authorized_worker']):
            center,=(q for q in accepted if q['action']==1.)
            signed=[q for q in rows if q['kind']=='signed']
            if len(signed)!=8 or center['axis_order']<-.0001 or center['axis_order']>1.0001:
                raise ArithmeticError(f'Incomplete signed result for {case}')
            entry.update(axis_order=center['axis_order'], min_order=center['min_order'],max_order=center['max_order'],
                         wave_frequency=center['wave_frequency'],current=center['current'])
            for coord in ('waist','depth'):
                energies={q['offset']:q['total_inside_and_drives'] for q in signed if q['coordinate']==coord}
                energies[0.]=center['total_inside_and_drives']
                entry[coord]={name:(energies[b]-energies[a])/(b-a) for name,a,b in (
                    ('minus',-.002,0.),('zero',-.001,.001),('plus',0.,.002),('double',-.002,.002))}
            cmd=out['command']
            file=ROOT/cmd[cmd.index('--save')+1]
            fields[str(file.relative_to(ROOT))]=digest(file)
        else:
            entry['diagnostic']=out['stderr_tail']
        summary[case]=entry
    src=ROOT/'software/stage1_solver/tools'
    scripts=['mouth_live_waist_'+s+'.py' for s in (
        'order','fields','bounded','support','active_wave','motion','wave','transient','order_readout')]
    scripts += ['mouth_exchange_carrier.py','mouth_exchange_carrier_density.py','mouth_exchange_shear.py',
                'mouth_exchange_geometry.py','mouth_conversion_shape_mode.py','mouth_neck_order_check.py',
                'mouth_sealed_shape_relax.py','mouth_exchange_interior.py']
    sources=[src/s for s in scripts]
    sources.extend(ROOT/'research/pde_ledger_v3/_measurements'/s for s in
                   ('mouth_conversion_shape.json','mouth_sealed_shape_relax.json'))
    sources.append(BASE/'transient_positive.py')
    times={}
    for name in (name for name in names if name.startswith('time')):
        result=runs[name]
        rows=list(map(json.loads,result['literal_stdout'].splitlines()))
        traces=[q for q in rows if q['kind']=='time' and 'relative_trace_error' in q]
        if (result['guard_exit']!=0 or len(rows)!=5 or len(traces)!=2
                or {q['phase'] for q in traces}!={'pure_strain','pure_kinetic'}
                or any(q['reached_bound'] or q['relative_trace_error']>.004 or
                       q['relative_energy_drift']>1e-7 or q['time']!=q['goal_time'] for q in traces)):
            raise ValueError(f'Incomplete phase-resolved control for {name}')
        times[name]=dict(axis_min=min(q['axis_order_min'] for q in traces),
                         axis_max=max(q['axis_order_max'] for q in traces),
                         timestep_error=max(q['relative_trace_error'] for q in traces),
                         energy_error=max(q['relative_energy_drift'] for q in traces),
                         execution_source='mouth_live_waist_transient.py' if name.endswith('_final') else 'transient_positive.py')
    print(json.dumps(dict(time_controls=times,base='bca2e63c',
        scope='Interior order and two active drives. Exterior order energy is a separate measured contribution. Time control is a fixed-density, fixed-phase, one-order-coordinate and phase-resolved single-wave Galerkin Hamiltonian.',
        input_sha256={str(p.relative_to(ROOT)):digest(p) for p in sources},
        saved_field_sha256=fields, summary=summary, runs=runs),indent=2))


if __name__=='__main__':
    main()
