"""Follow independent density-wave directions selected on the pressure-prepared 4D receiver."""
import argparse
from hashlib import sha256
from itertools import permutations, product
import json
from pathlib import Path

import numpy as np


LEDGER = Path('research/pde_ledger_v3/_measurements/mouth_bulk_temperature_steps.json')
BATCH = Path('_scratch/s11c/thermal87_batch02')
COLD = Path('_scratch/s11c/iso4_pressure_glass87_02/native/endpoint.npz')
Q_MAX = 12.


def determinant(matrix):
    """Exact for integer Gram matrices of size at most four."""
    n = len(matrix)
    total = 0
    for p in permutations(range(n)):
        parity = sum(p[i] > p[j] for i in range(n) for j in range(i+1, n))
        term = -1 if parity % 2 else 1
        for i, j in enumerate(p):
            term *= int(matrix[i, j])
        total += term
    return total


def grid(box):
    extent = np.floor(Q_MAX*np.linalg.norm(box, axis=1)/(2*np.pi)).astype(int)
    indices = np.array([h for h in product(*(range(-k,k+1) for k in extent)) if h > (0,0,0,0)])
    q = np.linalg.norm(2*np.pi*indices @ np.linalg.inv(box).T, axis=1)
    return indices[q<=Q_MAX]


def intensity(indices, frac):
    return abs(np.exp(2j*np.pi*(indices @ frac.T)).mean(axis=1))**2


def independent_peaks(indices, power):
    selected = []
    for i in np.argsort(-power, kind='stable'):
        candidates = indices[selected+[i]]
        if determinant(candidates @ candidates.T):
            selected.append(i)
            if len(selected) == 4:
                return selected
    raise ValueError('Reciprocal grid did not span all four directions')


def saved_reading(name, manifest):
    folder = BATCH / name
    receipt = json.loads((folder/'outcome.json').read_text())
    child = json.loads((folder/'child-outcome.json').read_text())
    limits = json.loads((folder/'limit-validation.json').read_text())
    raw = (folder/'stdout').read_bytes()
    if not (receipt['exitCode']==child['exitCode']==0 and receipt['childOutcome']==child
            and child['guardReason'] is None and limits['verified']
            and int(limits['actual']['memory.swap.max'])==0
            and manifest['runs'][name]['stdout_sha256']==sha256(raw).hexdigest()):
        raise ValueError(f'Invalid source for {name}')
    setup,*frames = map(json.loads,raw.splitlines())
    if [f['frame'] for f in frames] != list(range(49)):
        raise ValueError('Incomplete source')
    return setup,frames


def state(path):
    with np.load(path,allow_pickle=False) as f:
        if 'pos' in f:
            box=f['box']; frac=f['pos'] @ np.linalg.inv(box)
        else:
            box=f['cell']; frac=f['y']
    assert frac.shape==(128,4) and box.shape==(4,4)
    return frac,box


def compare(start, finish, cold, source, endpoint):
    pos,box = start; final,final_box = finish
    h = grid(box)
    a = intensity(h,pos)
    ix = independent_peaks(h,a)
    chosen = h[ix]
    b = intensity(chosen,final)
    # Same torus indices register homogeneous changes in the volume. Absolute intensities are translation invariant.
    assert abs(max(a) - source['diffraction']) < 1e-10
    final_grid = grid(final_box)
    final_power = intensity(final_grid,final)
    assert abs(max(final_power) - endpoint['diffraction']) < 1e-10
    cold_power = intensity(chosen,cold[0])
    qend = np.linalg.norm(2*np.pi*chosen @ np.linalg.inv(final_box).T,axis=1)
    return dict(pre_ranked_indices=chosen.tolist(),pre_intensity=a[ix].tolist(),post_same_intensity=b.tolist(),
                cold_same_intensity=cold_power.tolist(),post_same_q=qend.tolist(),post_strongest=float(max(final_power)),
                original_neighbors=endpoint['initial'],original_missing=endpoint['original_missing'],
                lost_fraction=endpoint['original_missing']/endpoint['initial'],
                largest_cluster=endpoint['largest_cluster'],starting_reciprocal_pairs=len(h))


def main(out):
    record = json.loads(LEDGER.read_bytes())
    setup,_ = saved_reading('r0h0_prepare',record)
    assert sha256(COLD.read_bytes()).hexdigest()==setup['source_sha256']
    cold = state(COLD)
    h=grid(cold[1]); power=intensity(h,cold[0]); selected=independent_peaks(h,power)
    result = {'q_max':Q_MAX, 'normalization':'|sum exp(2 pi i h dot y)/N|^2',
              'uniform_independent_fixed_index_mean':1/128,
              'cold_reference':{'indices':h[selected].tolist(),'intensities':power[selected].tolist(),'sha256':setup['source_sha256']},
              'histories':[]}
    for hist in range(3):
        name=f'r4h{hist}_'; psetup,prep=saved_reading(name+'prepare',record); usetup,un=saved_reading(name+'undriven',record)
        assert psetup['temperature']==usetup['temperature']==1.024 and prep[-1]['diffraction']==un[0]['diffraction']
        source=BATCH/(name+'prepare')/'native/end.npz'; end=BATCH/(name+'undriven')/'native/end.npz'
        assert sha256(source.read_bytes()).hexdigest()==usetup['source_sha256']
        a=compare(state(source),state(end),cold,un[0],un[-1]);a.update(history=hist,files={str(f):sha256(f.read_bytes()).hexdigest() for f in (source,end)})
        result['histories'].append(a)
    out.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps(result,allow_nan=False))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('output',type=Path)
    main(parser.parse_args().output)
