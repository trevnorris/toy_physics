"""Read the archived held R4 work, field rematches and unsuccessful step 167."""
import argparse
from collections import Counter
import hashlib
import json
import math
from pathlib import Path
import sys


def read(path):
    return json.loads(path.read_bytes())


def sha(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def require(test, message):
    if not test:
        raise ValueError(message)


def measure(study):
    primary = study/'primary'
    old = study.parent/'r4_rematch39_01/primary'
    guard = study/'primary_guard'
    receipt, child = read(guard/'outcome.json'), read(guard/'child-outcome.json')
    limits = read(guard/'limit-validation.json')
    require(receipt['exitCode']==child['exitCode']==1 and receipt['childOutcome']==child
            and child['guardReason'] is None and receipt['limitsVerified'] and limits['verified']
            and int(limits['actual']['memory.swap.max'])==0, 'Unexpected guarded exit')
    require('primary: guard returned 1' in read(study/'PAIRED_ISSUE.json')['reason'], 'Coordinator outcome')
    absent = ('PAIRED_DONE.json','paired_held.json','secondary','secondary_guard','pair_guard',
              'primary/HELD.json','primary/endpoint.json','primary/stationary_state.npz')
    require(all(not (study/name).exists() for name in absent), 'Unexpected stationary or paired output')
    review = read(study/'gate_review.json')
    manifest = read(study/'held_pin/manifest.json')
    pin = study/'held_pin/software/stage1_solver/tools'
    require(review['passed'] and all(sha(pin/name)==key for name,key in manifest.items()), 'Solver source pin')
    initial_sha = sha(old/'rematched.npz')
    require(initial_sha==read(old/'rematched.json')['sha256']==review['evidence'][
        str((old/'rematched.npz').relative_to(study.parents[2]))], 'Unverified restart state')
    events = [json.loads(line) for line in (guard/'stdout').read_text().splitlines()]
    rows, rematches = [], []
    start_energy = previous = events[0]['total_energy']
    step = 108
    for event in events:
        if event['kind']=='held_fresh_fields':
            require(event['step']==step and event['total_energy']<=previous, f'Rematch energy at {step}')
            folder = primary if rows else old
            state = read(folder/(f'matched_{step:05}.json' if rows else 'rematched.json'))
            sign = read(folder/(f'order_work_{step:05}.json' if rows else 'gate1.json'))
            workfile = folder/(f'order_work_{step:05}.npz' if rows else 'order_work.npz')
            energies = state['total_energy'] if rows else state['total']
            require(state['sha256']==sha(folder/state['file']) and energies==event['total_energy']
                    and sign['passed'] and sign['increment']==1e-9 and sign['total']['positive']==sign['free_zero']
                    and sign.get('work_sha256',sign.get('archive_sha256'))==sha(workfile)
                    and max(event['field'].values())<=1e-8, f'Incomplete fresh field state at {step}')
            rematches.append(dict(step=step,energy=energies,field_work=energies-previous,
                       crystal_residual=event['native_label'],native_field_residuals=event['field'],
                       front=event['front'], flow=event['current'],state_sha256=state['sha256'],
                       free_zero=sign['free_zero'],sign_counts=sign['total'],
                       max_abs_gauge_over_physical=sign['max_abs_gauge_over_positive_physical']))
            previous = energies
            continue
        step += 1
        require(event['kind']=='held_crystal_step' and event['outer']==step, 'Discontinuous accepted trajectory')
        state = read(primary/f'step_{step:05}.json')
        s, m = state['search'], state['material']
        energy = state['total_energy']
        require(state['source_sha256']==initial_sha and sha(primary/state['file'])==state['sha256']
                and state['labels_fixed_order'] and event['shape_locked'], f'Step {step} inputs')
        require(event['work']==s['energy_change']<0 and energy<previous and
                abs((energy-previous)-s['energy_change'])<=math.ulp(energy)
                and event['alpha']==s['alpha']>1e-4 and m['chart']['chart_det_global_lower']>0
                and m['min_Jinv_global']>0, f'Step {step} acceptance')
        rows.append(dict(step=step,state_sha256=state['sha256'],total_energy=energy,work=s['energy_change'],
                    alpha=s['alpha'],crystal_force=m['max_free_label_force'],force_l2=m['gradient_l2'],
                    chart_lower=m['chart']['chart_det_global_lower'],sampled_Jinv=m['min_Jinv_global'],
                    cg_reason=s['cg_reason'],cg_iterations=s['cg_iterations']))
        previous = energy
    final = read(primary/f'step_{step:05}.json')
    issue = read(primary/'ISSUE.json')
    require(step==166 and len(rows)==58 and [v['step'] for v in rematches]==[108,128,148]
            and final['sha256']=='4b84a0d2ad53889a1d619e68960df70233ceb589fc4fe6c336c6e59e832a2694'
            and issue['witness']['material']==final['material'] and issue['type']=='InfeasibleLabels',
            'Endpoint differs from the held primary stop')
    witness = issue['witness']
    tests = witness['tests']
    searches = {}
    for name, attempts in [('Newton',tests[:14]), ('projected Jacobi',tests[14:])]:
        require([t['alpha'] for t in attempts]==[2.**-j for j in range(14)], name+' step fractions')
        feasible = [t for t in attempts if t['reason']=='Insufficient resolved signed work']
        require(len(feasible)==9 and all(t['site']['delta']>0 for t in feasible), name+' feasible work')
        searches[name] = dict(reason_counts=dict(Counter(t['reason'] for t in attempts)),
                              smallest_test=attempts[-1],tests=attempts)
    record = dict(solver_commit='05a1b666',starting_state=initial_sha,command=sys.orig_argv,
                tool_sha256=sha(Path(__file__)),source_manifest_sha256=sha(study/'held_pin/manifest.json'),
                guard=receipt, held_shape=dict(entrance=4.,depth=final['depth'],exit=final['exit_radius']),
                accepted_steps=rows,freshly_equilibrated=rematches,absent_outputs=list(absent),
                total_energy_change=previous-start_energy,fixed_field_work_sum=math.fsum(v['work'] for v in rows),
                best_active_force=min(rows,key=lambda v:v['crystal_force']),
                accepted_cg_reasons=dict(Counter(v['cg_reason'] for v in rows)),
                minimum_accepted_fraction=min(v['alpha'] for v in rows),
                minimum_accepted_chart_lower=min(v['chart_lower'] for v in rows),
                lowest_sampled_physical_Jinv=min(v['sampled_Jinv'] for v in rows),
                final_material=final['material'],final_sha256=final['sha256'],failure=issue['reason'],
                trial_cg_reason=witness['cg_reason'],trial_cg_iterations=witness['cg_iterations'],searches=searches)
    return record


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('study',type=Path)
    parser.add_argument('output',type=Path)
    args=parser.parse_args()
    data=measure(args.study.resolve())
    with args.output.open('x') as stream:
        json.dump(data,stream,indent=2,allow_nan=False)
        stream.write('\n')
    print(json.dumps({k:data[k] for k in ('final_sha256','accepted_cg_reasons','total_energy_change',
                     'best_active_force','minimum_accepted_fraction','minimum_accepted_chart_lower')},allow_nan=False))
