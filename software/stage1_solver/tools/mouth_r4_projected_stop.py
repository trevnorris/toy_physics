"""Reduce the archived R4 projected active steps and the stopped melt-only solve."""
import argparse
import json
from math import fsum, isclose
from pathlib import Path

import numpy as np

from mouth_r4_stage2_outcome import digest, read, verify_guard


def read_lines(path):
    return [json.loads(line) for line in path.read_bytes().splitlines() if line]


def read_endpoint(folder, name):
    record = read(folder/(name + '.json'))
    path = folder/record['file']
    if digest(path) != record['sha256']:
        raise ValueError(f'Changed {name} checkpoint')
    with np.load(path, allow_pickle=False) as file:
        if (int(file['radius']), int(file['winding']), float(file['epsilon'])) != (4, 1, 1.e-8):
            raise ValueError('Unexpected physical state')
        return record, {key: file[key] for key in ('x', 'labels', 'field')}


def measurements(area):
    folder = area/'primary'
    guard = area/'primary_guard'
    verify_guard(guard, 1)
    verify_guard(area/'status15_guard', 0)
    eligibility = read(area/'eligibility.json')
    pins = read(area/'sources.json')
    root = area.parents[2]
    if (not eligibility['passed'] or any(digest(root/path) != sha
            for path,sha in eligibility['artifacts'].items())):
        raise ValueError('Eligibility or frozen inputs changed')
    if any(digest(area/'pin/software/stage1_solver/tools'/path) != sha for path,sha in pins.items()):
        raise ValueError('Pinned evaluator changed')
    invocation = read(guard/'invocation.json')
    start = root/eligibility['starting_state']
    if (Path(invocation['command'][-1]) != start or digest(start) != eligibility['starting_sha256']
            or read(folder/'ISSUE.json')['error'] != 'HeldStop' or not (area/'ISSUE.json').is_file()):
        raise ValueError('Incorrect starting state or termination')
    absent = [area/'DONE.json',area/'secondary_guard',area/'secondary',area/'paired_held.json',
              folder/'HELD.json',folder/'signed_forces.json',*folder.glob('continuation_*')]
    if any(p.exists() for p in absent):
        raise ValueError('Unexpected stationary or continued branch')
    records, states = {}, {}
    for key in ('matched','last_active','held_stop','failed_trial'):
        records[key], states[key] = read_endpoint(folder,key)
    with np.load(start) as saved:
        for s in states.values():
            if not np.array_equal(saved['x'],s['x']):
                raise ValueError('Held shape moved')
    for a,b,same in [('matched','last_active','field'),('last_active','held_stop','labels'),
                     ('held_stop','failed_trial','field')]:
        change = 'field' if same == 'labels' else 'labels'
        if (not np.array_equal(states[a][same],states[b][same])
                or np.array_equal(states[a][change],states[b][change])):
            raise ValueError('Invalid provenance of the rematched or trial state')

    steps = read_lines(guard/'stdout')
    last = records['last_active']
    directions = sorted(folder.glob('direction_*.npz'))
    if (len(directions)!=len(steps) or last['outer']+1!=len(steps)
            or records['matched']['outer']!=last['outer']):
        raise ValueError('Incomplete active trajectory')
    for i,(path,s) in enumerate(zip(directions,steps)):
        if (s['outer']!=i or path.stem!=f'direction_{i:06d}' or s['kind']!='projected_active_step'
                or s['alpha']<=1.e-4 or s['energy_change']>=0 or s['chart_lower']<=0):
            raise ValueError(f'Invalid active step at {i}')
        with np.load(path) as f:
            slope = float(f['gradient']@f['direction'])
        if not (slope<0 and s['energy_change']<=1.e-4*s['alpha']*slope):
            raise ValueError('Reduced-energy Armijo inequality failed')
        s.update(slope=slope)
    if (steps[-1]['material_force']!=last['search']['active_force']
            or steps[-1]['energy_change']!=last['search']['energy_change']
            or not isclose(steps[-1]['slope'],last['search']['slope'],rel_tol=1.e-13)):
        raise ValueError('Final saved step differs from log')
    stop = records['held_stop']
    decision = stop['stop']['inner']
    if (decision!=records['failed_trial']['inner'] or decision['stationary']
            or decision['reason']!='No resolved Armijo descent'):
        raise ValueError('Inner solve disposition differs')
    logs = sorted(folder.glob('inner_*.jsonl'))
    inner_target = 3.e-8
    converged = []
    for f in logs:
        rows = read_lines(f)
        final = rows[-1]
        if any(r['i']!=i or r['last']['delta']>=0 or r['last']['chart_lower']<=0
               for i,r in enumerate(rows[1:],1)):
            raise ValueError('Inconsistent accepted inner step')
        if f!=logs[-1]:
            if final['melt_residual']>inner_target:
                raise ValueError('Earlier inner solution is unresolved')
            converged.append(final['melt_residual'])
    if (final['melt_residual']!=decision['residual'] or final['i']!=decision['accepted_steps']
            or final['melt_residual']<=inner_target):
        raise ValueError('Inner stop disagrees with archived history')
    physical = stop['current']['excess_energy']
    gauge = stop['field']['order_step_energy']['new_gauge']
    print(f"active_steps={len(steps)} last_active_force={steps[-1]['material_force']:.12g} "
          f"inner_steps={final['i']} inner_residual={final['melt_residual']:.12g} "
          f"inner_target={inner_target:g} front={stop['front']:.12g}")
    return dict(solver_commit='d72f1d27',parent_commit='4c3c0a59',epsilon=1.e-8,
                held_shape=dict(depth=stop['depth'],exit_radius=stop['exit_radius']),
                state_sha256={key:r['sha256'] for key,r in records.items()},
                starting_sha256=eligibility['starting_sha256'],pre_descent_gradient_gate=eligibility['artifacts'][
                    'research/pde_ledger_v3/_measurements/mouth_r4_full_label_gate.json'],
                shape_held_in_all_saves=True,accepted_active_steps=steps,
                fixed_field_active_energy_changes_sum=fsum(s['energy_change'] for s in steps),
                min_active_fraction=min(s['alpha'] for s in steps),
                physical_force_target=3.e-5,rematched_physical_fields=stop['field'],
                physical_energy=physical,crystal_energy=stop['elastic_energy'],gauge_energy=gauge,
                total_energy=fsum((physical,gauge)),gauge_over_crystal=gauge/stop['elastic_energy'],
                current=stop['current'],fully_melted_depth=stop['front'],
                earlier_converged_inner_solves=len(converged),earlier_max_terminal_unit_residual=max(converged),
                failed_inner=dict(target=inner_target,reported=decision,first=rows[0],last=final,
                    last_rejected_fraction=2*decision['alpha'],unaccepted_exit_fraction=decision['alpha'],
                    same_reported_residual_since=next(r['i'] for r in reversed(rows)
                        if r['melt_residual']!=final['melt_residual'])+1),
                stationary_shape_forces=None,second_gauge='unadmitted',stronger_gauge_initializers=0,
                paired_comparison=None,shape_release=False)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('study',type=Path)
    parser.add_argument('output',type=Path)
    args=parser.parse_args()
    result=measurements(args.study.resolve())
    with args.output.open('x') as stream:
        json.dump(result,stream,indent=2,allow_nan=False)
        stream.write('\n')
