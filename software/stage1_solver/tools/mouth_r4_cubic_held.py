"""Paired held-mouth Newton solves on the approved two-metric cubic collar."""
import argparse
import hashlib
from pathlib import Path

import numpy as np

from mouth_exact_checkpoint_control import digest, read
from mouth_r4_aligned_rematch import FIELDS, FIELD_TOL, STRENGTH
from mouth_r4_collar166_state import SOURCE as STEP166_SHA
from mouth_r4_segmented_held import run
from mouth_r4_smooth_collar_reference import enable_smooth_collar
from mouth_r4_smooth_endpoint import certified_state, compare
from mouth_r4_stage2_outcome import verify_guard
from mouth_whole_r4_relax import atomic_json


R4_COEFFICIENT_SHA = '4c09ddff651e5b501e4b58014bd75b109be36878ed5e8d27b954709482ef7b12'


def approved_reference(body):
    proof = enable_smooth_collar(body)
    coefficients = np.ascontiguousarray(body.buffers['collar_coeff'], dtype='<f8')
    if hashlib.sha256(coefficients).hexdigest() != R4_COEFFICIENT_SHA:
        raise ValueError('The material reference differs from the approved radius-4 cubic')
    return proof


def reviewed_inputs(study):
    repo = study.parents[2]
    audit = study.parent / 'r4_crease167_01'
    gate = read(study / 'gate_review.json')
    if gate['passed'] is not True or not gate['verdict']:
        raise ValueError('Independent review did not clear cubic held descent')
    evidence = gate['evidence']

    def check(path):
        if digest(path) != evidence[str(path.relative_to(repo))]:
            raise ValueError(f'Reviewed evidence differs at {path}')

    for name in ('derivatives_lean', 'frozen_controls'):
        check(audit / name / 'DONE.json')
        if read(audit / name / 'DONE.json')['passed'] is not True:
            raise ValueError(f'Uncleared pre-descent controls in {name}')
    for case, strength in STRENGTH.items():
        case_dir = audit / 'frozen_controls' / ('rematch_' + case)
        verify_guard(case_dir.with_name(case_dir.name + '_guard'), 0)
        check(case_dir / 'rematched.npz')
        check(case_dir / 'result.json')
        record = read(case_dir / 'result.json')
        sign = record['sign']
        if (record['state_sha256'] != digest(case_dir / 'rematched.npz') or record['source_sha256'] != STEP166_SHA
                or record['epsilon'] != strength or record['outer'] != 166 or not record['passed']
                or not all(0 <= record['rematched']['field'][k] <= FIELD_TOL for k in FIELDS)
                or not sign['passed'] or sign['total']['positive'] != sign['free_zero'] or sign['free_zero'] == 0):
            raise ValueError(f'Invalid fixed-label rematch or one-sided order work for {case}')
    return audit / 'frozen_controls/rematch_primary'


def solve(study, case):
    origin = reviewed_inputs(study)
    record = read(origin / 'result.json')
    if case == 'pair':
        compare(study, primary_initial=record['state_sha256'], aligned=True, configure_reference=approved_reference)
        return
    if case == 'primary':
        start = origin / 'rematched.npz', 166
        initial_match = record
    else:
        primary, _, _ = certified_state(study, 'primary', STRENGTH['primary'], record['state_sha256'])
        start = study / 'primary/stationary_state.npz', primary['outer']
        initial_match = None
    run(study, case, start=start, configure_reference=approved_reference, reviewed_match=initial_match)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('study', type=Path)
    parser.add_argument('case', choices=[*STRENGTH, 'pair'])
    args = parser.parse_args()
    study = args.study.resolve()
    try:
        solve(study, args.case)
    except Exception as error:
        folder = study / args.case
        if folder.is_dir():
            witness = {k:getattr(error,k) for k in ('diagnostics','details','witness') if hasattr(error,k)}
            atomic_json(folder / 'ISSUE.json', dict(type=type(error).__name__, reason=str(error), **witness))
        raise
