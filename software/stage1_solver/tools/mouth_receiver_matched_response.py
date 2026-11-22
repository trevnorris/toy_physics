"""Measure shear and compression after an independent geometry meets the precursor test."""
import argparse
import hashlib
import json
from pathlib import Path

from mouth_receiver_material_screen import signature
from mouth_receiver_pressure_readout import weighted


def load(path):
    return json.loads(path.read_bytes())


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def measure(root):
    done = load(root / 'DONE.json')
    classification = root / 'classification.json'
    assert digest(classification) == done['readout_sha256']
    choice = load(classification)['selected']
    assert choice == done['selected'] and choice is not None
    selected_hash = digest(root / (choice + '_back.json'))
    groups = load(root / 'design.json')['jobs']
    assert set(done['results']) == {c['name'] for group in groups.values() for c in group}
    result = {}
    for name, cases in groups.items():
        trajectories = {}
        for c in cases:
            path = root / 'runs' / (c['name'] + '.json')
            assert digest(path) == done['results'][c['name']]
            row = load(path)
            assert row['source_sha256'] == selected_hash
            assert (row['history'], row['seed'], row['strain'], row['pressure_external'],
                    row['requested_dt'], row['axis_index'], row['temperature']) == (
                        c['history'], c['seed'], c['strain'], c['p_ext'], c['dt'], 1, 1.088)
            assert row['warm'] == c['warm'] and row['measure'] == c['duration']
            assert len(row['blocks']) == 16
            trajectories[c['label']] = row
        zero = trajectories['zero']
        prepared = signature(zero['prepared'])['peak']
        arms = {}
        for key, row in trajectories.items():
            assert row['prepared'] == zero['prepared']
            blocks = row['blocks']
            last = blocks[-1]['structure']
            density = weighted(blocks[-4:], 'density')
            drift = density / weighted(blocks[-8:-4], 'density') - 1
            current = {
                'prompt_stress': row['initial_reading']['shear'],
                'late_stress': weighted(blocks[-4:], 'shear'),
                'density': density, 'relative_density_change': drift,
                'minimum_cluster': min(b['structure']['max_cluster'] for b in blocks[-4:]),
                'pressure': weighted(blocks[-4:], 'pressure'),
                'log_volume': weighted(blocks[-4:], 'log_volume'),
                'old_neighbors_absent': last['absent'] / last['original'],
                'nonaffine_rms': last['nonaffine_rms'],
                'end_peak': signature(row['final'])['peak'],
            }
            current['meets_structure_test'] = (
                density >= 1.2 and abs(drift) <= .10 and current['minimum_cluster'] >= 116
                and max(prepared, current['end_peak']) <= .20 and current['old_neighbors_absent'] >= .5
                and current['nonaffine_rms'] >= 1.)
            arms[key] = current
        response = {}
        for label, key, step in [('full', 'g', .02), ('half', 'h', .01)]:
            plus, minus = arms[key+'p'], arms[key+'m']
            initial = (plus['prompt_stress'] - minus['prompt_stress']) / (2*step)
            late = (plus['late_stress'] - minus['late_stress']) / (2*step)
            response[label] = dict(prompt=initial, late=late, fraction=late/initial if initial else None)
        compression = (arms['pm']['log_volume'] - arms['pp']['log_volume']) / .5
        result[name] = dict(prepared_peak=prepared, response=response, compliance=compression, arms=arms)
        print(json.dumps(dict(group=name, prepared_peak=prepared, force=response, compliance=compression,
                              number_passing_structure=sum(v['meets_structure_test'] for v in arms.values()),
                              density_range=[min(v['density'] for v in arms.values()),max(v['density'] for v in arms.values())],
                              peak_range=[min(v['end_peak'] for v in arms.values()),max(v['end_peak'] for v in arms.values())]),
                         allow_nan=False))
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('root', type=Path)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    args.output.write_text(json.dumps(measure(args.root), indent=2, allow_nan=False) + '\n')
