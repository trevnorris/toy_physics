"""Independent reduction of the signed bonded-rim and stationary-wall tests.

Use only completed fields and per-step records. The physical surface target
was calculated from the saved unshifted fields before any signed fields were
produced. Keep the two one-sided derivatives in case their centered average
masks a cusp where the two finite-element grids separate.
"""
import argparse
import json
from math import fsum, ulp
from pathlib import Path

import numpy as np

from mouth_rim_pinned_comparison import digest, read, validate_guard


def derivatives(rows, key, baseline):
    pos_w, pos_r = [s[key] for s in rows[0]['pair']]
    neg_w, neg_r = [s[key] for s in rows[1]['pair']]
    base_w, base_r = baseline
    h = rows[0]['step']
    return dict(positive=fsum((base_w, -base_r, -pos_w, pos_r))/h,
                negative=fsum((neg_w, -neg_r, -base_w, base_r))/h,
                centered=fsum((neg_w, -neg_r, -pos_w, pos_r))/(2*h))


def grid(folder, previous_zero, physical, prior):
    validate_guard(folder)
    result = read(folder/'result.json')
    meta = read(folder/'input.json')
    level = folder.name
    if (meta != result['input'] or level != meta['grid'] or
            meta['zero_sha'] != physical['zero_sha256'] or
            digest(previous_zero/'zero.npz') != meta['zero_sha'] or
            prior['zero_state_sha256'] != meta['zero_sha']):
        raise ValueError('Inconsistent input or zero-state archives')
    for key, other in [('state_sha','state_sha256'), ('wall_sha','wall_sha256'),
                       ('axes_sha','reference_axes_sha256'), ('fields','field_sha256'),
                       ('load_sha','load_sha256')]:
        if meta[key] != prior['inputs'][other]:
            raise ValueError(f'{level} input changed: {key}')
    with np.load(previous_zero/'zero.npz', allow_pickle=False) as saved:
        zero = {k: saved[k] for k in ('wound', 'rest')}
    e0 = [s['energy'] for s in meta['zero']]
    for row, old in zip(meta['zero'],prior['zero'],strict=True):
        if abs(row['energy']-old['rechecked_energy']) > 2.e-12 or row['stationarity'] >= meta['tolerance']:
            raise ValueError(f'{level} different or unconverged zero')
    by_name = {key: {} for key in ('four', 'eight', 'stationary')}
    rows_checked, gradients, port_errors, roundoffs = [], [], [], []
    for row in result['signed']:
        filename = Path(row['file'])
        if filename.name != row['file'] or row != read(folder/filename.with_suffix('.json')):
            raise ValueError(f'Summary disagrees with signed record: {filename}')
        if digest(folder/filename) != row['sha256'] or row['jacobian'] < .95:
            raise ValueError(f'Unusable signed field or map: {filename}')
        with np.load(folder/filename, allow_pickle=False) as saved:
            for key, state in zip(('wound', 'rest'), row['pair'], strict=True):
                delta = np.max(abs(saved[key]-zero[key]),axis=1)
                if not np.array_equal(delta,state['change_max']) or not np.isfinite(saved[key]).all():
                    raise ValueError(f'Saved state disagrees with the recorded field change: {filename}')
                gradients.append(state['stationarity'])
                port_errors.append(abs(state['feed']-state['bulk']))
                if key=='wound' and min(state['feed'],state['bulk']) <= 0:
                    raise ValueError(f'Inward drainage lost at {filename}')
            roundoffs.append(sum(ulp(s['rechecked_energy']) for s in row['pair'])/abs(row['step']))
        table = by_name[row['case']]
        if row['step'] in table:
            raise ValueError(f'Duplicate signed pair: {filename}')
        table[row['step']] = row
        rows_checked.append(row['sha256'])
    if max(gradients) >= meta['tolerance'] or max(port_errors) >= 3.e-5:
        raise ValueError(f'Invalid stationarity or current continuity on {level}')
    near, far = [next(row for row in physical['results'] if row['window']==label and row['degree']==2)
                 for label in ('near','far')]
    complement = near['common_reaction'] - far['common_reaction']
    refs = dict(near=-near['outward'][1], far=-far['outward'][1],
                complementary_material_work=complement, far_on_near=-far['outward'][1]+complement)
    slopes = {}
    for name, table in by_name.items():
        step = meta['step']
        if set(table) != {step, -step, step/2, -step/2}:
            raise ValueError(f'Missing signed histories: {level} {name}')
        fields = {}
        for label, h in [('full',step), ('half',step/2)]:
            pair = [table[h], table[-h]]
            work = {k: derivatives(pair, k, z) for k,z in (
                ('rechecked_energy',e0), ('frozen_energy',e0),
                ('swept_layer_energy',(0.,0.)), ('crystal_strain_increment',(0.,0.)))}
            if name != 'stationary':
                f = work['rechecked_energy']['centered']
                work['minus_near_surface'] = f-refs['near']
                work['minus_far_surface_on_near_motion'] = f-refs['far_on_near']
            fields[label] = work
        slopes[name] = fields
    gap = {h:{k: slopes['four'][h]['rechecked_energy'][k]-slopes['eight'][h]['rechecked_energy'][k]
              for k in ('positive','negative','centered')} for h in ('full','half')}
    return dict(grid=level, references=refs, steps=slopes, four_minus_eight=gap,
                versus_old_fixed_gap=dict(old=prior['four_minus_eight']['duct_relative']['half']['centered'],
                    new=slopes['stationary']['half']['rechecked_energy']['centered']),
                max_stationarity=max(gradients), max_current_error=max(port_errors),
                energy_representation_scale=max(roundoffs), signed_state_sha256=rows_checked,
                result_sha256=digest(folder/'result.json'), zero_sha256=meta['zero_sha'])


def run(a):
    if read(a.batch/'DONE.json') != {'double':0,'local':0} or (a.batch/'ISSUE.json').exists():
        raise ValueError('Batch incomplete')
    frozen = read(a.batch/'src_hashes.json')
    if any(digest(a.batch/'src'/name) != sha for name,sha in frozen.items()):
        raise ValueError('Frozen solver source changed')
    old = read(a.prior)
    ref = read(a.surfaces)
    if digest(a.prior) != ref['record_before_shift_sha256']:
        raise ValueError('The prior zero-state manifest changed')
    rows = []
    for level in ('double', 'local'):
        p = next(k for k in old['cases'] if k['inputs']['grid']==level)
        s = next(k for k in ref['cases'] if k['grid']==level)
        rows.append(grid(a.batch/level, a.previous/level, s, p))
    drift = {name: {h:rows[1]['steps'][name][h]['rechecked_energy']['centered']
                      -rows[0]['steps'][name][h]['rechecked_energy']['centered'] for h in ('full','half')}
             for name in ('stationary','four','eight')}
    return dict(cases=rows, fine_minus_doubled=drift, surface_sha256=digest(a.surfaces), source_sha256=frozen)


if __name__ == '__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ('batch','previous','prior','surfaces'):
        parser.add_argument(name,type=Path)
    print(json.dumps(run(parser.parse_args()),allow_nan=False,indent=2))
