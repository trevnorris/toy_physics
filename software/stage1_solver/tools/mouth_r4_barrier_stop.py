"""Reduce the first held R4 determinant-gauge run, without advancing its fields."""
import argparse
from collections import Counter
from math import fsum
from pathlib import Path
import json

from mouth_r4_stage2_outcome import digest, read, verify_guard


def measure(scratch):
    area = scratch / 'r4_barrier_held01'
    work = scratch / 'r4_barrier_work01'
    folder = area / 'primary'
    run = area / 'primary_guard'
    issue = read(area / 'ISSUE.json')
    reduced = read(area / 'gauge_held.json')
    verify_guard(run, 1)
    verify_guard(area / 'held_gate_guard', 1)
    verify_guard(area / 'precision_eligibility_guard', 0)
    verify_guard(scratch / 'r4_mixed_barrier_work01' / 'independent_guard', 0)
    approved = read(area / 'precision_eligibility.json')
    corrected = read(work / 'PRECISION_GATE.json')
    raw_sha = digest(work / 'GATE.json')
    gate_sha = digest(work / 'PRECISION_GATE.json')
    if (read(work / 'GATE.json')['passed'] or not corrected['passed'] or corrected['failed']
            or not read(scratch / 'r4_mixed_barrier_work01' / 'VERIFIED.json')['passed']
            or raw_sha != approved['raw_gate_sha256'] or raw_sha != corrected['raw_gate_sha256']
            or gate_sha != approved['gate_sha256'] or gate_sha != issue['gate_sha256']
            or any(digest(work / p) != h for p, h in approved['guards_sha256'].items())):
        raise ValueError('Unverified derivative measurement or admission')
    if (issue['outcomes'] != {'primary': 1} or (area / 'secondary_guard').exists()
            or reduced['passed'] or reduced['continuation_allowed'] or (folder / 'HELD.json').exists()):
        raise ValueError('The case disposition differs from a single held stop')

    row = reduced['rows']['primary']
    record = read(folder / 'held_stop.json')
    last = read(folder / 'last_active.json')
    state_sha = digest(folder / record['file'])
    if not (state_sha == record['sha256'] == last['sha256'] == digest(folder / last['file'])
            == row['accepted_state']['sha256']):
        raise ValueError('Last accepted state and held stop differ')

    steps = [json.loads(line) for line in (run / 'stdout').read_bytes().splitlines()]
    steps = [s for s in steps if s['kind'] == 'barrier_active_step']
    stop = record['stop']
    if ([s['step'] for s in steps] != list(range(1, stop['accepted_steps'] + 1))
            or steps != row['accepted_active_steps'] or steps[-1]['alpha'] != last['alpha']
            or any(s['alpha'] < 1e-4 or s['energy_change'] >= 0 for s in steps)):
        raise ValueError('Accepted steps do not match the energy and feasibility record')
    probes = stop['probes']
    floor = [p for p in probes if p['alpha'] == 1e-4]
    location = floor[-1]['witness']
    profile = row['accepted_state']
    corner = next(p['one_sided_on_this_saved_state'] for p in profile['read_at_rejection_coordinates']
                  if p['witness']['cell'] == location['cell']
                  and p['one_sided_on_this_saved_state']['local'] == location['chart_fractions'])
    field = record['field']
    current = record['current']
    parent = read(scratch / 'r4_eulerian_stage2_01' / 'primary' / 'label_stall.json')
    material = stop['material']
    return dict(
        radius=4, winding=1, epsilon=1e-8, solver_base='091e8c4d', saved_state_sha256=state_sha,
        held_solver_sha256=approved['held_source_sha256'], verified_derivative_gate_sha256=gate_sha,
        held_readback_sha256=digest(area / 'gauge_held.json'), order_and_flow_rematch=field['order_step_energy'],
        stop=stop['reason'], held_shape=dict(depth=record['depth'], exit_radius=record['exit_radius']),
        prior_held_fields=dict(front=parent['frozen_field_front'], current=parent['held_field_current']['feed_current']),
        sampled_melted_front_depth=record['front'],
        current={k:current[k] for k in ('feed_current','bulk_current','weak_disk_current','port_error',
                                        'disk_minus_feed_current','dense_total_subsonic_share','minimum_total_mach')},
        energies=dict(physical=current['excess_energy'], crystal=record['elastic_energy'],
                      gauge=material['gauge_energy'], total=fsum((current['excess_energy'], material['gauge_energy'])),
                      gauge_to_crystal=material['gauge_to_elastic'],
                      fixed_fields_active_steps_delta=fsum(s['energy_change'] for s in steps)),
        free_crystal_force=material['max_free_label_force'], material_target=3e-5,
        physical_force_components=profile['physical_force_components'],
        melt_inner_unit_residual=material['max_unit_melt_residual'], accepted_steps=steps,
        last_step_fraction=stop['last_accepted_alpha'], min_step_fraction=stop['min_alpha'],
        rejected_predictors=dict(Counter(p['reason'] for p in probes)), trials_at_floor=floor,
        accepted_corner=corner, active_quadrature_minimum=profile['worst_active'],
        global_quadrature_minimum=profile['worst_global'], chart_bounds=profile['chart'],
        pre_label_field_residuals={k:field[k] for k in ('phase_raw','amplitude_raw','carrier_metric',
                                                        'order_raw','order_kkt','order_projected')},
        stationary_shape_forces=row['six_stationary_forces'], second_gauge='not admitted',
        paired_held_comparison=None, shape_release=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('scratch', type=Path)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    result = measure(args.scratch)
    with args.output.open('x') as stream:
        json.dump(result, stream, indent=2, allow_nan=False)
        stream.write('\n')
    print(result['stop'], result['free_crystal_force'], result['accepted_corner']['Jinv'])
