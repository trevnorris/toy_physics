"""Recover the missing inner search vector from the fixed accepted fractions.

This is an arithmetic readback of the last failed fixed-field inner call.
There is no line search, equilibrium solve, or adoption of any physical state.
A single differing journal scalar or endpoint label invalidates the vector.
"""
import argparse
import hashlib
import json
import sys
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np


def sha(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def registered_start(opt, labels):
    result = labels.copy()
    for c, cm in enumerate(opt._components()):
        old = (labels - opt.base)[..., c].ravel()
        value = old[cm.representative].copy()
        value[cm.free] = old[cm.source[cm.free]]
        value[cm.held] = 0
        m = ~opt.stencil
        registered = opt.base[..., c] + value[cm.groups].reshape(m.shape)
        result[..., c][m] = registered[m]
    return result


def run(area, output):
    pin = area/'pin/software/stage1_solver/tools'
    sys.path.insert(0, str(pin))
    from mouth_cp10_label_grid import get_cp10
    from mouth_r4_eulerian_relax_labels import LabelRelaxation, absmax, chart_orientation
    from mouth_r4_melt_barrier import MeltBarrier
    from mouth_r4_stage2_outcome import verify_guard

    manifest = area/'sources.json'
    if sha(manifest) != '625c07dc0ed83f93df23959c4310c336c49e3ef61cd751fdbe6df0f9dbf37c49':
        raise ValueError('Changed source manifest')
    for name, digest in json.loads(manifest.read_bytes()).items():
        if sha(pin/name) != digest:
            raise ValueError(name)
    verify_guard(area/'primary_guard', 1)
    primary = area/'primary'
    path = primary/'inner_0000034.jsonl'
    rows = [json.loads(line) for line in path.read_bytes().splitlines()]
    states = []
    for name, digest in [('held_stop','ae9dd7452d5291c6af95cdcfa23ae1cb6d0705fe29fc7d32cfc596608854f008'),
                         ('failed_trial','23172a606edb2b0707f3ebd132a42bae06e9bb336e4712ad34b154fe0d7839d1')]:
        state_path = primary/(name+'.npz')
        if sha(state_path) != digest:
            raise ValueError('Unverified inner input')
        with np.load(state_path, allow_pickle=False) as data:
            states.append((data['labels'], data['field']))
    (start, fields), (endpoint, failed_fields) = states
    if not np.array_equal(fields, failed_fields):
        raise ArithmeticError('Field changed during inner operation')
    b, _, axes = get_cp10(area.parent)
    order = np.ascontiguousarray(fields[2,b.ids].transpose(0,2,1))
    opt = LabelRelaxation(b, start, order, 1e-8)
    labels = registered_start(opt, start)
    with TemporaryDirectory() as work:
        barrier = MeltBarrier(b, opt, Path(work))
        value, grad, parts, det = barrier.sample(labels)
        gradient = barrier.extract(grad, 'dead')
        last = None
        pairs = []
        checks = 0
        for index, row in enumerate(rows):
            actual = dict(i=index, unit_energy=float(value), melt_residual=absmax(gradient), last=last, **parts)
            if actual != row:
                raise ArithmeticError(f'Archived inner row {index} disagrees: {actual!r} versus {row!r}')
            checks += 1
            y = gradient.copy()
            dots = []
            for step, change, inv in reversed(pairs):
                dot = inv*np.dot(step,y)
                y -= dot*change
                dots.append(dot)
            y = barrier.pre(y)
            for (step,change,inv), dot in zip(pairs,reversed(dots)):
                y += step*(dot-inv*np.dot(change,y))
            direction = -y
            if np.dot(direction,gradient) >= 0:
                pairs.clear()
                direction = -barrier.pre(gradient)
            step = barrier.melt_direction(direction)
            scale = max(1., *(absmax(step[...,c])/(.25*np.median(np.diff(ax))) for c,ax in enumerate(axes)))
            direction /= scale
            step /= scale
            if index == len(rows)-1:
                break
            alpha = rows[index+1]['last']['alpha']
            trial = labels+alpha*step
            value, next_grad, parts, new_det = barrier.sample(trial)
            next_gradient = barrier.extract(next_grad, 'dead')
            dstep = trial-labels
            midpoint = ((labels-opt.base)+.5*dstep).reshape(-1,3)
            dv = np.sum(barrier.weight[:,None]*((new_det-det)-np.log1p((new_det-det)/det)), dtype=np.longdouble)
            dv += np.sum(dstep.reshape(-1,3)*(opt.K@midpoint), dtype=np.longdouble)
            orientation = chart_orientation(axes,trial,opt.active)
            last = dict(alpha=alpha, search_normalization=scale, delta=float(dv),
                        chart_lower=orientation['chart_det_global_lower'])
            displacement, change = alpha*direction, next_gradient-gradient
            sy = np.dot(displacement,change)
            if sy > 1e-12*np.linalg.norm(displacement)*np.linalg.norm(change):
                pairs.append((displacement,change,1./sy))
                pairs = pairs[-10:]
            labels, det, gradient = trial, new_det, next_gradient
    if not np.array_equal(labels.view(np.uint64),endpoint.view(np.uint64)) or checks != 106:
        raise ArithmeticError('Final fixed path differs from archived failed inner state')
    if np.any(step[opt.stencil]) or gradient@direction>=0:
        raise ArithmeticError('Recovered vector changes the active trace or fails the descent test')
    npz_path = output/'direction.npz'
    np.savez_compressed(npz_path, packed=direction, step=step, labels_sha256=hashlib.sha256(endpoint.tobytes(order='C')).hexdigest())
    record = dict(command=sys.orig_argv, source_sha256=sha(Path(__file__)), input_manifest_sha256=sha(manifest),
                  inner_journal_sha256=sha(path), accepted_fractions_used=len(rows)-1,
                  journal_rows_all_scalars_equal=checks, final_labels_bit_exact=True, history_length=len(pairs),
                  max_active_motion=absmax(step[opt.stencil]), slope=float(gradient@direction),
                  largest_unit_gradient=absmax(gradient), largest_label_step=absmax(step),
                  direction_sha256=sha(npz_path))
    (output/'reconstruction.json').write_text(json.dumps(record, indent=2)+'\n')
    print(json.dumps(record),flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('area',type=Path)
    parser.add_argument('output',type=Path)
    args = parser.parse_args()
    args.output.mkdir(exist_ok=False)
    run(args.area.resolve(),args.output.resolve())
