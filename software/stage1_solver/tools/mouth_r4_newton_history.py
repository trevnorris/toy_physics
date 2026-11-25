"""Recover held R4 timings, accepted work and packed force norms from sealed histories.

Original journal birth and final-write times enclose only the iterative part
of a melt solve. The original worker had no per-call instrumentation.
"""
import argparse
import hashlib
import json
import os
from datetime import datetime, timezone
from pathlib import Path
import subprocess

import numpy as np


def read(path):
    return json.loads(path.read_bytes())


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def time_ns(text):
    base, fraction = text.split('.')
    digits, zone = fraction.split()
    return int(datetime.strptime(base + zone, '%Y-%m-%d %H:%M:%S%z').timestamp()) * 10**9 + int(digits)


def extents(paths):
    output = subprocess.check_output(['stat', '--printf=%w\t%y\n', '--', *paths],
        text=True, env={**os.environ,'TZ':'UTC','LC_ALL':'C'})
    return {p:tuple(map(time_ns,line.split('\t'))) for p,line in zip(paths,output.splitlines(),strict=True)}


def central_norms(path, expected):
    with np.load(path, allow_pickle=False) as data:
        g = data['gradient']
        maximum = float(np.max(abs(g)))
        if maximum != expected:
            raise ValueError(f'Saved gradient has unexpected state: {path}')
        return dict(linf=maximum,l2=float(np.linalg.norm(g)),packed_free_active=int(g.size))


def audit(folder, baseline):
    run = folder/'primary'
    guard = folder/'primary_guard'
    source = folder/'pin/software/stage1_solver/tools'
    inputs = [source/(name+'.py') for name in ('mouth_r4_projected_held','mouth_r4_smooth_held',
              'mouth_r4_barrier_held','mouth_r4_melt_barrier')]
    pins = read(folder/'sources.json')
    if any(sha(p) != pins[p.name] for p in inputs):
        raise ValueError('Pinned loop differs from the saved production manifest')
    saved = read(run/'last_active.json')
    shape_sha = sha(run/'last_active.npz')
    if (saved['outer'] != 29 or saved['sha256'] != shape_sha or
            shape_sha != 'd2e764212c08bdcacfbc18121812ceb484f2eaffcd0b44f46991a9cf03120202'):
        raise ValueError('Expected positively oriented held step 29')
    invocation = read(guard/'invocation.json')
    start = datetime.fromisoformat(invocation['startedUtc']).astimezone(timezone.utc)
    start_ns = int(start.timestamp())*10**9 + start.microsecond*1000
    events = [json.loads(x) for x in (guard/'stdout').read_bytes().splitlines() if x.strip()]
    steps = {e['outer']:e for e in events if e['kind']=='smooth_active_step'}
    rematches = {e['outer']:e for e in events if e['kind']=='smooth_rematched'}
    count = saved['outer']
    if sorted(steps) != list(range(1,count+1)) or sorted(rematches)!=list(range(0,count,4)):
        raise ValueError('Unexpected acceptance or field-rematch sequence')
    directions = sorted(run.glob('direction_*.npz'))
    inner = sorted(run.glob('inner_*.jsonl'))
    op_count = count+len(rematches)
    if len(directions)!=count or not all(p.name==f'inner_{i:07d}.jsonl' for i,p in enumerate(inner)):
        raise ValueError('Missing directions or a skipped melt operation')
    if not op_count<=len(inner)<=op_count+1:
        raise ValueError('Melt calls must be classified explicitly before using this readout')
    times=extents([*directions,*inner])
    op=0
    previous_end=start_ns
    previous_energy=np.longdouble(rematches[0]['energy'])
    rows=[]
    artifacts=[*inputs,*directions,run/'last_active.json',guard/'invocation.json',guard/'stdout']
    for k in range(1,count+1):
        rec=steps[k]
        direction=directions[k-1]
        if direction.name!=f'direction_{k-1:06d}.npz' or rec['family']!='local_diagonal_harmonic':
            raise ValueError('Unexpected direction marker or search family')
        end=times[direction][1]
        interval=(end-previous_end)/1e9
        if not previous_end<end:
            raise ValueError('Reversed file-time interval')
        phase=rematches.get(k-1)
        journals=[]
        for role in (['fixed_trace','accepted_trial'] if phase else ['accepted_trial']):
            p=inner[op]
            begin,finish=times[p]
            if not previous_end<=begin<=finish<=times[direction][0] or (journals and begin<journals[-1]['end_ns']):
                raise ValueError(f'Overlapping/incorrectly assigned melt journal at {k}: {p}')
            records=[json.loads(x) for x in p.read_bytes().splitlines()]
            iterations=len(records)-1
            if ([r['i'] for r in records]!=list(range(iterations+1)) or
                    records[-1]['melt_residual']>1e-5):
                raise ValueError(f'Unfinished or malformed melt journal {p}')
            journals.append(dict(role=role,op=op,begin_ns=begin,end_ns=finish,
                journal_seconds=(finish-begin)/1e9,iterations=iterations))
            artifacts.append(p)
            op+=1
        log2=-np.log2(rec['alpha'])
        if not log2.is_integer() or log2>=14:
            raise ValueError('Unexpected outer fraction')
        # For the published trial family this final harmonic pass tried
        # 1, 1/2, ... alpha. Any prior normal-set refinements are unlogged.
        minimum=int(log2)
        exact=len(saved['search']['rejections']) if k==count else None
        if exact is not None and exact<minimum:
            raise ValueError('Recorded rejection list shorter than final search')
        before=(phase if phase else steps[k-1])['active_force']
        pre=central_norms(direction,before)
        active_E=previous_energy+np.longdouble(rec['energy_change'])
        matched=rematches.get(k)
        final_E=np.longdouble(matched['energy']) if matched else active_E
        post=(matched if matched else rec)['active_force']
        after=central_norms(directions[k],post) if k<count else dict(linf=post,l2=None)
        step=dict(step=k,wall_seconds=interval,accepted_marker_ns=end,journals=journals,
            other_seconds=interval-sum(j['journal_seconds'] for j in journals),
            rematched_before=bool(phase),rematched_after=bool(matched),
            rematch_field_time_upper_s=(journals[0]['begin_ns']-previous_end)/1e9 if phase else 0.,
            projection_and_early_backtrack_upper_s=(journals[-1]['begin_ns']-
                (journals[0]['end_ns'] if phase else previous_end))/1e9,
            outer_step_alpha=rec['alpha'],outer_backtracks_at_least=minimum,outer_backtracks_exact=exact,
            pre_label_force=pre,post_label_force=after,immediate_post_label_linf=rec['active_force'],
            fixed_field_total_energy_change=rec['energy_change'],
            subsequent_rematch_energy_change=float(final_E-active_E),
            net_total_energy_change=float(final_E-previous_energy),total_energy_after=float(final_E),
            minimum_sampled_Jinv=rec['min_Jinv'])
        rows.append(step)
        previous_end=end
        previous_energy=final_E
    if op!=op_count:
        raise ValueError('Unaccounted completed melt calls')
    top_sites=None
    if baseline:
        control=read(baseline)
        # The baseline evaluator reports all physical forces on the accepted
        # snapshot. Never substitute the gradient before step 29.
        done=read(baseline.parent.parent/'baseline_guard/outcome.json')
        if not done['limitsVerified'] or done['exitCode'] or done['childOutcome']['guardReason'] is not None:
            raise ValueError('Unverified endpoint force readback')
        if control['state_sha256']!=shape_sha:
            raise ValueError('Endpoint L2 diagnostic has a different state')
        row=rows[-1]['post_label_force']
        if abs(row['linf']-control['material']['max_free_label_force'])>2e-13:
            raise ValueError('Endpoint active-force anchor differs')
        row['l2']=control['active_gradient_l2']
        row['packed_free_active']=control['active_gradient_count']
        top_sites=control['top5_nodes_by_vector_norm']
        artifacts.append(baseline)
    def window(chunk):
        duration=sum(s['wall_seconds'] for s in chunk)
        out=dict(steps=[chunk[0]['step'],chunk[-1]['step']],step_count=len(chunk),
            elapsed_s=duration,mean_minutes_per_step=duration/60/len(chunk),
            total_energy_change=sum(s['net_total_energy_change'] for s in chunk),
            accepted_fixed_field_delta_E=sum(s['fixed_field_total_energy_change'] for s in chunk),
            outer_backtracks_at_least=sum(s['outer_backtracks_at_least'] for s in chunk),
            rematches=sum(s['rematched_before'] for s in chunk))
        for role in ('fixed_trace','accepted_trial'):
            jobs=[j for s in chunk for j in s['journals'] if j['role']==role]
            sec=sum(j['journal_seconds'] for j in jobs)
            out[role]=dict(completed=len(jobs),inner_updates=sum(j['iterations'] for j in jobs),
                           journal_s=sec,journal_share=sec/duration)
        other=sum(s['other_seconds'] for s in chunk)
        out.update(rejected_completed_inner_solves=0,rejected_inner_journal_s=0.,
                   other_s=other,other_share=other/duration,
                   field_rematch_time_upper_s=sum(s['rematch_field_time_upper_s'] for s in chunk),
                   projection_and_early_backtrack_upper_s=sum(s['projection_and_early_backtrack_upper_s'] for s in chunk))
        return out
    interrupted=None
    if len(inner)>op:
        p=inner[-1]
        journal=[json.loads(x) for x in p.read_bytes().splitlines() if x.endswith(b'}')]
        if journal:
            interrupted=dict(journal=p.name,recorded_updates=journal[-1]['i'],
                last_recorded_residual=journal[-1]['melt_residual'],journal_seconds=(times[p][1]-times[p][0])/1e9)
    report=dict(state_sha256=shape_sha,source_sha256=sha(Path(__file__)),
        schema=dict(clock='Invocation to acceptance markers. A direction-file final write follows successful fixed-order line search and precedes outer state serialization. Startup is included in step 1.',
        inner='Original file birth to final write encloses iteration records, excluding initial sampling and return overhead. Sequential calls: one fixed-trace journal per rematch and one accepted-trial journal per accepted update. Contiguous operation numbers show no additional rejected inner solves.',
        other='All non-journalled work. Outer quadrature, determinant checks/projections and order/phase iterations had no individual clocks. Upper bounds overlap this category and must not be added to it.',
        backtracks='Minimum counts are exact halvings in the successful harmonic pass. Earlier failed direction-refinement passes were overwritten, except in the final checkpoint. No extra full inner solve was launched for a rejected trial.',
        residuals='Euclidean norm and largest component of the packed free active gradient. Pre is the accepted direction archive. Post includes any same-step field rematch. Immediate post-step infinity norms precede rematches. L2 at step 29 requires a separate saved-state gradient readout.',
        energy='Fixed-field crystal-plus-gauge change equals the total change during that substep. The net change also includes any following same-step rematch. Re-anchored to logged float64 total energies at rematches.',
        tail='Time after accepted step 29, including any unfinished inner operation, is excluded from both windows.'),
        top5_paused_nodes_by_vector_norm=top_sites,
        input_hashes={str(p.relative_to(folder) if p.is_relative_to(folder) else p.name):sha(p) for p in artifacts},
        windows=[window(rows[:21]),window(rows[21:])],steps=rows,unaccepted_next_trial=interrupted)
    return report


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('run',type=Path)
    parser.add_argument('output',type=Path)
    parser.add_argument('--baseline',type=Path)
    args=parser.parse_args()
    result=audit(args.run.resolve(),args.baseline)
    with args.output.open('x') as f:
        json.dump(result,f,indent=2,allow_nan=False)
        f.write('\n')
    print(json.dumps(result['windows'],indent=2),flush=True)
