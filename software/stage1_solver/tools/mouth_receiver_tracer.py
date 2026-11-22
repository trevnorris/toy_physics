"""Reconstruct bulk tagged-site motion from the undriven pressure-hold endpoints."""
import argparse
from hashlib import sha256
import json
from pathlib import Path
from statistics import mean

import numpy as np


def reconstruct(batch):
    design = json.loads((batch / 'design.json').read_text())
    hashes = json.loads((batch / 'DONE.json').read_text())['results']
    histories = {}
    for group, cases in design['jobs'].items():
        case, = [c for c in cases if c['label'] == 'zero']
        raw = (batch / 'runs' / (case['name'] + '.json')).read_bytes()
        assert sha256(raw).hexdigest() == hashes[case['name']]
        r = json.loads(raw)
        assert r['strain'] == r['pressure_external'] == 0 and r['n'] == 128
        assert (r['history'], r['seed']) == (case['history'], case['seed'])
        start, end = r['prepared'], r['final']
        # The solver writes unwrapped Cartesian positions in a box with row vectors.
        # Map the final fractional positions back to the initial physical box.
        frac = np.linalg.solve(np.asarray(end['box']).T, np.asarray(end['pos']).T).T
        dx = frac @ np.asarray(start['box']) - np.asarray(start['pos'])
        dx -= dx.mean(axis=0)
        squares = dx*dx
        distances = np.sqrt(squares.sum(axis=1))
        rms = float(np.sqrt(squares.sum(axis=1).mean()))
        last = r['blocks'][-1]['structure']
        np.testing.assert_allclose(rms, last['nonaffine_rms'], atol=1e-10, rtol=1e-10)
        row = dict(native_duration=r['measure'], rms=rms,
                   w_rms=float(np.sqrt(squares[:,3].mean())),
                   parallel_rms=float(np.sqrt(squares[:,:3].sum(axis=1).mean())),
                   over_one_bond=int(np.count_nonzero(distances >= 1.0)),
                   lost_fraction=last['absent']/last['original'])
        histories[group] = row
        print(json.dumps(dict(history=group, **row), allow_nan=False))
    assert len(histories) == 6
    summary = dict(mean_rms=mean(h['rms'] for h in histories.values()),
                   mean_fraction_original_absent=mean(h['lost_fraction'] for h in histories.values()),
                   min_over_one_bond=min(h['over_one_bond'] for h in histories.values()),
                   max_over_one_bond=max(h['over_one_bond'] for h in histories.values()))
    print(json.dumps(summary, allow_nan=False))
    return dict(histories=histories, summary=summary)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('batch', type=Path)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    args.output.write_text(json.dumps(reconstruct(args.batch), indent=2, allow_nan=False) + '\n')
