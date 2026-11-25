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


def reviewed_step288(study):
    """Read the accepted production replay from the independently reviewed gates."""
    from mouth_r4_step288_cells import STEP287
    review = read(study / 'gate_review.json')
    if review['passed'] is not True or not review['verdict']:
        raise ValueError('Repaired projection has not passed independent review')
    repo = study.parents[2]
    evidence = review['evidence']
    replay = study / 'replay'
    state = replay / 'step_00288.npz'
    for path in (study/'permutations.json', replay/'result.json', replay/'step_00288.json', state):
        if digest(path) != evidence[str(path.relative_to(repo))]:
            raise ValueError(f'Reviewed gate artifact changed: {path}')
    for name in ('permutations_guard', 'replay_guard'):
        verify_guard(study/name, 0)
    for name in ('mouth_r4_corner_projection.py', 'mouth_r4_newton_crystal.py'):
        if digest(Path(__file__).with_name(name)) != evidence[f'software/stage1_solver/tools/{name}']:
            raise ValueError(f'The held solver differs from the reviewed step: {name}')
    permutations = read(study/'permutations.json')
    trial = read(replay/'result.json')
    step = read(replay/'step_00288.json')
    if not (permutations['passed'] and len(permutations['permutations'])==24 and trial['passed']
            and all(trial['checks'].values()) and trial['outer']==step['outer']==288
            and trial['epsilon']==1e-8 and trial['initial_sha256']==step['source_sha256']==STEP287
            and trial['state_sha256']==step['sha256']==digest(state)
            and trial['record_sha256']==digest(replay/'step_00288.json')
            and trial['orientation']['physical_global_lower']>=8.54e-8):
        raise ValueError('Unverified repaired step-288 production replay')
    return state, trial['state_sha256']


def solve(study, case, *, projector_replay=False):
    if projector_replay:
        path, primary_initial = reviewed_step288(study)
        first_step = 288
        initial_match = None
    else:
        origin = reviewed_inputs(study)
        initial_match = read(origin/'result.json')
        primary_initial = initial_match['state_sha256']
        first_step = 166
        path = origin/'rematched.npz'
    if case == 'pair':
        comparison = compare(study, primary_initial=primary_initial, aligned=True,
                             configure_reference=approved_reference)
        if projector_replay:
            endpoints = {}
            for name in STRENGTH:
                endpoint = read(study/name/'endpoint.json')
                sign, flow = endpoint['order_boundary'], endpoint['current']
                endpoints[name] = dict(epsilon=endpoint['epsilon'], step=endpoint['outer'],
                        total_energy=endpoint['total_energy'], crystal_energy=endpoint['crystal_energy'],
                        gauge_energy=endpoint['gauge_energy'], gauge_to_crystal=endpoint['gauge_to_crystal'],
                        shape_forces=dict(zip(comparison['shape_force_delta'], endpoint['shape_forces'], strict=True)),
                        sampled_melted_front=endpoint['sampled_melted_front'], native_residuals=endpoint['native_residuals'],
                        sign_positive=sign['total']['positive'], free_zero=sign['free_zero'],
                        order_increment=sign['increment'], max_abs_gauge_over_physical=sign['max_abs_gauge_over_positive_physical'],
                        feed_current=flow['feed_current'], dense_subsonic_share=flow['dense_total_subsonic_share'],
                        feed_minus_bulk=flow['port_error'],
                        physical_Jinv_lower=endpoint['orientation']['physical_certificate']['physical_global_lower'])
            comparison['stationary_endpoints'] = endpoints
            atomic_json(study/'paired_held.json', comparison)
        return
    if case == 'primary':
        start = path, first_step
    else:
        primary, _, _ = certified_state(study, 'primary', STRENGTH['primary'], primary_initial)
        start = study / 'primary/stationary_state.npz', primary['outer']
        initial_match = None
    run(study, case, start=start, configure_reference=approved_reference, reviewed_match=initial_match,
        rematch_start=projector_replay)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('study', type=Path)
    parser.add_argument('case', choices=[*STRENGTH, 'pair'])
    parser.add_argument('--projector-replay', action='store_true',
                        help='Start the primary from the independently verified repaired step 288')
    args = parser.parse_args()
    study = args.study.resolve()
    try:
        solve(study, args.case, projector_replay=args.projector_replay)
    except Exception as error:
        folder = study / args.case
        if folder.is_dir():
            witness = {k:getattr(error,k) for k in ('diagnostics','details','witness') if hasattr(error,k)}
            atomic_json(folder / 'ISSUE.json', dict(type=type(error).__name__, reason=str(error), **witness))
        raise
