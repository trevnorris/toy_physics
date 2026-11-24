"""Assemble the two R4 fixed-shape label stops from archived states and read-only checks."""
import argparse
import hashlib
import json
from math import fsum, isclose
from pathlib import Path

import numpy as np


def read(path):
    return json.loads(path.read_bytes())


def digest(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def verify_guard(folder, expected_exit):
    result = read(folder / 'outcome.json')
    child = result['childOutcome']
    limits = read(folder / 'limit-validation.json')
    if (not result['limitsVerified'] or not limits['verified']
            or result['exitCode'] != expected_exit or child['exitCode'] != expected_exit
            or child['guardReason'] is not None or int(limits['actual']['memory.swap.max']) != 0):
        raise ValueError(f'Unverified outcome at {folder}')


def current_budget(current):
    keys = ('feed_current', 'bulk_current', 'port_error', 'weak_disk_current',
            'incident_slice_current', 'receiving_slice_current', 'dense_total_subsonic_share',
            'disk_minus_feed_current', 'total_axial_stress')
    return {k: current[k] for k in keys}


def collect(scratch, case, epsilon, baseline):
    run = scratch / 'r4_eulerian_stage2_01'
    folder = run / case
    guard = run / (case + '_guard')
    verify_guard(guard, 1)
    if read(folder / 'issue.json')['error'] != 'LabelStall':
        raise ValueError(f'Unexpected endpoint at {folder}')
    for name in ('checkpoint_00011.json', 'physical_stop.json', 'equilibrium.json'):
        if (folder / name).exists():
            raise ValueError(f'Already advanced to {name}')
    state = read(folder / 'label_stall.json')
    sha = digest(folder / state['file'])
    force_dir = scratch / 'r4_spatial_stall_force01' / (
        'secondary_fixed_field' if case == 'secondary' else case)
    corner_dir = scratch / 'r4_eulerian_corner01'
    verify_guard(force_dir.with_name(force_dir.name + '_guard'), 0)
    verify_guard(corner_dir / (case + '_guard'), 0)
    review_guard = 'review_primary_margin_guard' if case == 'primary' else 'review_secondary_guard'
    verify_guard(run / review_guard, 0)
    forces = read(force_dir / 'readout.json')
    corner = read(corner_dir / (case + '.json'))
    review = read(run / ('review_' + case + '.json'))
    if not (sha == state['sha256'] == forces['state_sha256'] == corner['state_sha256']
            == review['saved_sha256'] and forces.get('complete') and state['iteration'] == 10):
        raise ValueError(f'Inconsistent fixed-state diagnostics for {case}')
    rows = [json.loads(line) for line in (guard / 'stdout').read_bytes().splitlines()]
    start = next(row for row in rows if row.get('event') == 'label_start')
    accepted = [row for row in rows if row.get('event') == 'label_accepted']
    label = state['label']
    if not (rows[0]['checkpoint'] == baseline['sha256'] and rows[0]['phase_drop'] == .4510084282857874
            and rows[0]['epsilon'] == epsilon == label['epsilon']
            and [r['iteration'] for r in accepted] == list(range(1,review['accepted_steps']+1))
            and all(row['energy_delta'] < 0 for row in accepted)):
        raise ValueError(f'History validation failed for {case}')
    limiting = corner['accepted_limiting_cell']
    if not isclose(limiting['Jinv'], review['accepted_corner']['Jinv'], rel_tol=1.e-8):
        raise ValueError(f'Independent one-sided physical determinants disagree for {case}')
    modes = rows[0]['coordinates']
    if set(modes) != set(forces['modes']) or len(modes) != 6:
        raise ValueError('Missing shape modes')
    derivatives = {}
    for name in modes:
        m = forces['modes'][name]
        signed = []
        for p in m['signed']:
            banks = p['banks']
            # Differences were taken before serializing the large common bath energy.
            one_sided = {s:b.get('one_sided') for s,b in banks.items()}
            central = fsum(one_sided.values())/2 if p['valid_central'] else None
            signed.append(dict(h=p['h'], feasible_banks={s:b['feasible'] for s,b in banks.items()},
                               one_sided=one_sided, central=central,
                               analytic_minus_central=None if central is None else m['force']-central))
        derivatives[name] = dict(analytic=m['force'],signed=signed)
    if label['stationary'] or len(accepted) == 0 or min(row['alpha'] for row in accepted) <= 0:
        raise ValueError(f'Unexpected label-search record {case}')
    current = state['held_field_current']
    total = fsum(forces['energy'])
    fixed_delta = fsum(row['energy_delta'] for row in accepted)
    before = fsum((current['excess_energy'], start['gauge_energy']))
    # The pre-label energy uses the rematched order and the original crystal load.
    if abs((total-before)-fixed_delta) > 2.e-8:
        raise ValueError(f'Accumulated energy changes disagree with direct signed-work base at {case}')
    return dict(epsilon=epsilon, input_sha256=baseline['sha256'], saved_state_sha256=sha,
                geometry={k:state[k] for k in ('depth','exit_radius','wall_polar_amplitude')},
                fixed_field_total_energy=total, initial_rematched_total_energy=before,
                initial_crystal_energy=start['elastic_energy'], crystal_energy=label['elastic_energy'],
                gauge_energy=label['gauge_energy'], gauge_over_crystal=label['gauge_to_elastic'],
                strict_label_energy_change=fixed_delta, accepted_steps=len(accepted),
                initial_free_label_force=start['max_free_label_force'],
                final_free_label_force=label['max_free_label_force'],
                final_free_force_by_component=label['free_force_by_component'],
                harmonic_residual=label['max_unit_harmonic_residual'],
                first_alpha=accepted[0]['alpha'], final_alpha=accepted[-1]['alpha'],
                final_energy_change=accepted[-1]['energy_delta'],
                last_rejected_alpha=review['last_rejected_alpha'], last_rejected=label['last_rejection'],
                chart_lower_all=label['min_chart_Jinv_bound'],
                chart_lower_active_cells=label['min_chart_Jinv_active_bound'],
                lowest_sampled_physical=label['worst_global'],
                lowest_sampled_ordered_physical=label['worst_active'],
                accepted_corner=limiting, active_incident_cells=[
                    c for c in corner['accepted_incident_cells'] if c['positive_vertices'] > 0],
                flow=current_budget(current), sampled_fully_melted_depth=state['frozen_field_front'],
                fields_before_label_descent=state['field_before_labels'],
                frozen_shape_work=derivatives, force_order=modes,
                diagnostic_hashes={p.name+'.'+p.parent.name:digest(p) for p in (
                    folder/'label_stall.json', force_dir/'readout.json', corner_dir/(case+'.json'),
                    run/('review_'+case+'.json'))})


def main(scratch, output):
    run = scratch / 'r4_eulerian_stage2_01'
    serial = read(run / 'SERIAL_ISSUE.json')
    if 'endpoint codes primary 1, secondary 1' not in serial['error']:
        raise ValueError(serial['error'])
    if any((run/name).exists() for name in ('gauge_11.json','gauge_final.json','SERIAL_DONE.json')):
        raise ValueError('A stationary-state gauge comparison or terminal success was recorded')
    history = scratch / 'r4_melted_survival01/coarse'
    base = read(history / 'latest.json')
    if base['iteration'] != 10 or digest(history / base['file']) != base['sha256']:
        raise ValueError('Wrong input checkpoint')
    cases = {case:collect(scratch,case,eps,base) for case,eps in [('primary',1.e-8),('secondary',1.e-10)]}
    with np.load(history / base['file']) as cp, np.load(run/'primary/label_stall.npz') as a, \
            np.load(run/'secondary/label_stall.npz') as b:
        if not all(np.array_equal(s['x'],cp['x']) and str(s['parent_sha256'])==base['sha256']
                   and int(s['radius'])==4 and int(s['winding'])==1 for s in (a,b)):
            raise ValueError('Shape or source changed')
        state_delta = {name:float(np.max(np.abs(a[name]-b[name]))) for name in ('x','field','labels')}
    primary, secondary = cases.values()
    comparison = dict(state_max_abs_difference=state_delta,
                      physical_energy_difference=(primary['crystal_energy']-secondary['crystal_energy']
                          if state_delta['field']==0 else None),
                      total_energy_difference=primary['fixed_field_total_energy']-secondary['fixed_field_total_energy'],
                      shape_force_differences={k:primary['frozen_shape_work'][k]['analytic']
                          - secondary['frozen_shape_work'][k]['analytic'] for k in primary['force_order']})
    cp_front = base['depth']*base['current']['shaft_order']['terminal_from_fraction']['disk']['0.0']
    corner_dir = scratch / 'r4_eulerian_corner01'
    for name,case in cases.items():
        verify_guard(corner_dir/(name + '_support_guard'), 0)
        support = read(corner_dir/(name + '_support.json'))
        if support['state_sha256'] != case['saved_state_sha256']:
            raise ValueError('Melt-only support readback belongs to another state')
        case['vertices_of_limiting_cell'] = support['vertices']
    verify_guard(run/'review_flux_guard', 0)
    currents = read(run/'review_common_currents.json')
    flux_names = dict(feed='feed_current',bulk='bulk_current',disk='weak_disk_current',
                      subsonic_share='dense_total_subsonic_share')
    for case in cases.values():
        for key, recorded in flux_names.items():
            if abs(currents[key]-case['flow'][recorded]) > 1.e-10:
                raise ValueError(f'Independent current disagreement for {key}')
    result = dict(command='mouth_r4_stage2_outcome.py _scratch/s11c OUTPUT', phase='fixed-shape Stage 2 stop',
                  release_solver_commit='01622952', stage1_clearance_commit='fb0347a8',
                  frozen_source_manifest_sha256=digest(run/'sources.json'),
                  force_convention='Minus the partial shape derivative at fixed saved fields, with production nodal transport. No stationary released derivative.',
                  checkpoint_11=False, stationary_gauge_comparison=False, reference=dict(
                  sha256=base['sha256'], energy=base['energy'], six_shape_forces=base['shape_force'],
                  flow=current_budget(base['current']), sampled_fully_melted_depth=cp_front),
                  cases=cases, primary_minus_secondary=comparison, independent_common_flow=currents)
    with output.open('x') as dest:
        json.dump(result, dest, indent=2, allow_nan=False)
        dest.write('\n')
    for case,r in cases.items():
        print(case, r['fixed_field_total_energy'], r['accepted_steps'], r['final_free_label_force'],
              r['final_alpha'], r['accepted_corner']['Jinv'], r['flow']['feed_current'], flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('scratch', type=Path)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    main(args.scratch, args.output)
