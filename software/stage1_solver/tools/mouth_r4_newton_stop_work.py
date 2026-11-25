"""Read the signed work history and the refusal at the end of the held R4 solve."""
import argparse
from collections import Counter
import hashlib
import json
import math
from pathlib import Path


def read(path):
    return json.loads(path.read_bytes())


def sha(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def require(condition, description):
    if not condition:
        raise ValueError(description)


def reduce_run(study):
    root = study.parents[2]
    primary = study / 'primary'
    pin = study / 'pin/software/stage1_solver/tools'
    entry = read(study / 'entry.json')
    sources = read(study / 'pinned.json')
    require(entry['passed'] and all(sha(root / k) == v for k, v in entry['files'].items()), 'entry inputs')
    require(all(sha(pin / k) == v for k, v in sources.items()), 'pinned worker sources')

    for name in ('baseline', 'oracle_guard_serial', 'timing', 'hessian', 'benchmark', 'primary'):
        folder = study / (name if name == 'oracle_guard_serial' else name + '_guard')
        limits = read(folder / 'limit-validation.json')
        outcome = read(folder / 'outcome.json')
        child = read(folder / 'child-outcome.json')
        require(outcome['limitsVerified'] and limits['verified'] and int(limits['actual']['memory.swap.max']) == 0,
                f'{name} resource limits')
        require(child == outcome['childOutcome'] and child['guardReason'] is None and
                child['exitCode'] == outcome['exitCode'] == int(name == 'primary'), f'{name} outcome')

    review = read(study / 'benchmark_review.json')
    bench = read(study / 'benchmark/result.json')
    require(review['passed'] and bench['passed'] and
            review['benchmark_sha256'] == read(study / 'BENCH_DONE.json')['sha256'] == sha(study / 'benchmark/result.json'),
            'benchmark eligibility')
    require(bench['saved'] == review['saved_state_sha256'] == sha(study / 'benchmark/last_active.npz'),
            'primary starting state')
    failure = read(study / 'PAIRED_ISSUE.json')
    issue = read(primary / 'ISSUE.json')
    require('primary stopped (1)' in failure['error'] and issue['error'] == 'InfeasibleLabels', 'primary issue')
    absent = ('primary/HELD.json', 'primary/stationary_state.npz', 'secondary', 'secondary_guard',
              'pair', 'pair_guard', 'paired_held.json', 'PAIRED_DONE.json')
    require(all(not (study / path).exists() for path in absent), 'unexpected pair or stationary certificate')
    accepted = read(primary / 'last_active.json')
    rematched = read(primary / 'rematched.json')
    for rec in (accepted, rematched):
        require(sha(primary / rec['file']) == rec['sha256'], rec['file'])

    journal = study / 'primary_guard/stdout'
    initial, *events = [json.loads(line) for line in journal.read_text().splitlines()]
    require(initial['kind'] == 'newton_initial' and initial['source_sha256'] == accepted['source_sha256'] == bench['saved'],
            'step-34 ancestry')
    columns = ['outer', 'energy', 'work', 'alpha', 'max_force', 'gradient_l2', 'min_Jinv', 'cg_iterations',
               'cg_reason', 'backtracks']
    rows, rematches = [], []
    last_energy = initial['energy']
    step = bench['outer']
    for item in events:
        energy = item['energy']
        if item['kind'] == 'newton_step':
            step += 1
            require(item['outer'] == step and energy < last_energy and item['work'] < 0 and
                    abs(energy - last_energy - item['work']) <= math.ulp(energy), f'crystal work at {step}')
            require(item['min_Jinv'] > 0 and item['alpha'] > 1e-4, f'accepted orientation and fraction at {step}')
            rows.append([item[k] for k in columns])
        else:
            require(item['kind'] == 'newton_rematched' and item['outer'] == step and
                    energy <= last_energy and item['field_residual'] <= 1e-8, f'field rematch at {step}')
            rematches.append(dict(step=step, energy=energy, work=energy-last_energy,
                                  max_field_residual=item['field_residual']))
        last_energy = energy
    require(step == accepted['outer'] and last_energy == accepted['energy'] and
            rematches[-1]['step'] == rematched['outer'] and rematches[-1]['energy'] == rematched['total_energy'],
            'saved endpoints versus journal')
    witness = issue['witness']
    require(witness['material'] == accepted['material'], 'rejected trial anchored to accepted state')
    trials = witness['tests']
    fractions = [2.**-i for i in range(14)]
    require([q['alpha'] for q in trials] == fractions * 2, 'two complete dyadic line searches')
    refusals = []
    for name, group in zip(('projected Newton', 'projected Jacobi'), (trials[:14], trials[14:])):
        reasons = Counter(q['reason'] for q in group)
        require(set(reasons) == {'Nonpositive or roundoff-unresolved chart determinant',
                                 'Label leaves the saved material domain'}, 'feasibility-only refusals')
        last = group[-1]
        require(last['site']['chart_det'] < -last['site']['roundoff_margin'], 'resolved negative determinant')
        refusals.append(dict(direction=name, reasons=reasons, final=last))
    work = [r[2] for r in rows]
    total_change = last_energy - initial['energy']
    closure = total_change - math.fsum(work + [v['work'] for v in rematches])
    require(abs(closure) <= len(events) * math.ulp(last_energy), 'signed work closure')
    return dict(scope='Fixed geometry, primary epsilon=1e-8; trial 109 refused. No stationary pair.',
                command='python3 software/stage1_solver/tools/mouth_r4_newton_stop_work.py '
                        '_scratch/s11c/r4_newton_held01',
                source_sha256={str(f.relative_to(study)): sha(f) for f in
                               (journal, primary / 'ISSUE.json', primary / 'last_active.npz', primary / 'rematched.npz')},
                start=dict(step=bench['outer'], energy=initial['energy'], active_Linf=initial['force'],
                           joint_L2=initial['material']['gradient_l2']),
                columns=columns, accepted_steps=rows, field_rematches=rematches,
                fixed_field_work=math.fsum(work), total_energy_change=total_change, work_closure_error=closure,
                attempted_step=step+1, rejected_directions=refusals)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('study', type=Path)
    print(json.dumps(reduce_run(parser.parse_args().study.resolve()), indent=2, allow_nan=False))
