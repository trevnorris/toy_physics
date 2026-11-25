"""Separate physical bath and numerical gauge work from the R4 log-depth crystal probe."""
import argparse
import json
from pathlib import Path
import sys

import numpy as np

from mouth_cp10_field_shape import Bath
from mouth_r4_aligned_derivatives import SignedGauge, load
from mouth_r4_smooth_endpoint import LD, promoted_bath
from mouth_r4_stage2_outcome import digest, read, verify_guard


def measure(study, case, work):
    archive = study/'derivatives'/(case+'_shape')/'record.json'
    verify_guard(study/'derivatives'/(case+'_shape_guard'), 1)
    native = read(archive)
    reference = native['rows'][0]
    state_hash = digest(study/case/'rematched.npz')
    if (reference['mode'] != 'log_depth' or native['state_sha256'] != state_hash
            or native['passed'] or reference['central'][0]['minus_analytic'] < 1e-7):
        raise ValueError('Expected the held step-108 log-depth stop')
    driver, trial, _, order = load(study,case,work)
    gauge = SignedGauge(trial,driver.reference_depth,order)
    bath = Bath(trial,driver.base,driver.fields,driver.reference_depth)
    unit = np.eye(6,dtype=LD)[0]
    h = LD(2)**-19
    ratio = LD(driver.epsilon)
    def at(delta):
        return (promoted_bath(bath,trial,driver.reference_depth,delta*unit),
                gauge.evaluate(driver.labels,delta*unit))
    neg, pos = at(-h), at(h)
    pieces = (neg[0]-pos[0])/(2*h)
    mid = dict(h=float(h), bath=float(pieces.sum(dtype=LD)), physical_bath_force_terms=pieces.tolist(),
               unit_gauge=float(np.sum(neg[1]-pos[1],dtype=LD)/(2*h)))
    mid['physical_bath_force_terms'] = [float(v) for v in pieces]
    records = []
    for item in (reference['central'][0],mid,reference['central'][1]):
        bath_error = item['bath']-reference['bath']
        gauge_error = float(ratio*(LD(item['unit_gauge'])-LD(reference['unit_gauge'])))
        result = dict(h=item['h'], bath_force=item['bath'], bath_force_error=bath_error,
                      unit_gauge_force=item['unit_gauge'], weighted_gauge_force_error=gauge_error,
                      noncell_force_error=bath_error+gauge_error)
        if 'crystal' in item:
            result.update(archived_crystal_force=item['crystal'], archived_total_force_error=item['minus_analytic'],
                          archived_crystal_force_error=item['crystal']-reference['elastic'])
            result['unaccounted_force_error'] = (item['minus_analytic']-
                (result['archived_crystal_force_error']+result['noncell_force_error']))
            if abs(result['unaccounted_force_error']) > 2e-15:
                raise ArithmeticError('Component-sum disagreement with the archived force')
        records.append(result)
    return dict(state_sha256=state_hash, wider_and_narrower_record_sha256=digest(archive),
                gauge_strength=float(ratio), bath_analytic_force=reference['bath'],
                gauge_unit_analytic_force=reference['unit_gauge'], ladder=records,
                middle_physical_bath_force_terms=mid['physical_bath_force_terms'])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('study',type=Path)
    parser.add_argument('output',type=Path)
    args = parser.parse_args()
    out = args.output.resolve()
    out.mkdir()
    report = dict(command=sys.orig_argv, script_sha256=digest(Path(__file__)), cases={},
                  convention='Forces equal minus the energy derivative in log depth. Errors are centered minus analytic forces.',
                  physical_bath_terms=['amplitude gradient','excess fluid potential','coherent phase',
                  'order gradient','order potential','wall','exterior crystal'])
    for key in ('primary','secondary'):
        report['cases'][key] = measure(args.study.resolve(),key,out)
    if any(len(c['middle_physical_bath_force_terms'])!=len(report['physical_bath_terms']) for c in report['cases'].values()):
        raise ArithmeticError('Unexpected bath energy terms')
    with (out/'noncell.json').open('x') as f:
        json.dump(report,f,indent=2,allow_nan=False)
        f.write('\n')
    print(json.dumps(report,allow_nan=False),flush=True)


if __name__ == '__main__':
    main()
