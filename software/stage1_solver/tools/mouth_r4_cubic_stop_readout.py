"""Read the last accepted cubic R4 state without a field solve or trial step."""
import argparse
from pathlib import Path

import numpy as np

from mouth_exact_checkpoint_control import digest, read
from mouth_r4_collar166_state import load_held, readout
from mouth_r4_cubic_held import approved_reference, reviewed_step288
from mouth_r4_cubic_physical_orientation import physical_certificate
from mouth_r4_stage2_outcome import verify_guard
from mouth_whole_r4_relax import atomic_json


def measure(study, output):
    case = study / 'primary'
    verify_guard(study / 'primary_guard', 1)
    state_record = sorted(case.glob('step_*.json'))[-1]
    record = read(state_record)
    saved = case / record['file']
    gate = read(study / 'gate_review.json')
    if 'step288_sha256' in gate:
        _, source_sha = reviewed_step288(study)
    else:
        source_sha = gate['reviewed']['primary']['fresh_fields']['state_sha256']
    if digest(saved) != record['sha256'] or record['source_sha256'] != source_sha:
        raise ValueError('Accepted state or rematched source hash differs')
    output.mkdir(exist_ok=False)
    driver, trial, solid, _ = load_held(study, output, 1e-8, solve_fields=True)
    approved_reference(trial.b)
    with np.load(saved, allow_pickle=False) as state:
        if int(state['radius']) != 4 or int(state['winding']) != 1 or not np.array_equal(state['x'], driver.x):
            raise ValueError('Unexpected held geometry')
        driver.labels, driver.fields = state['labels'], state['field']
    fresh = readout(driver, trial, solid)
    if (fresh['chart'] != record['material']['chart']
            or fresh['crystal_energy'] != record['material']['elastic_energy']
            or fresh['peak_active'] != record['material']['max_free_label_force']
            or abs(fresh['total_energy'] - record['total_energy']) > 2 * np.spacing(fresh['total_energy'])):
        raise ArithmeticError('Read-only evaluation does not reproduce the saved state')
    orientation = physical_certificate(trial.b, driver.labels, fresh['chart']['chart_det_global_lower'])
    if not 0 < orientation['physical_global_lower'] <= fresh['min_Jinv']:
        raise ArithmeticError('Physical determinant bound failed at the accepted state')
    atomic_json(output / 'result.json', dict(step=record['outer'], source_sha256=record['sha256'],
                record_sha256=digest(state_record), readout=fresh, physical_certificate=orientation))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('study', type=Path)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    measure(args.study.resolve(), args.output.resolve())
