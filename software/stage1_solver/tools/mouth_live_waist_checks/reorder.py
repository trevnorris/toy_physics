"""One-sided order variations within an empty patch of the held radius-one throat."""
import argparse
import hashlib
import json
from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from mouth_live_waist_active_wave import ActiveWave
from mouth_live_waist_order import BIAS, Meridian, Order, load_for, state


def directions(chi, mask):
    """Keep a whole-patch control and separate the released depth and radial halves."""
    yield 'whole', mask
    zz, rr = np.nonzero(mask)
    z, r = np.indices(chi.shape)
    for label, values, grid in (('depth', zz, z), ('radius', rr, r)):
        mid = np.median(values)
        yield label + '_low', mask & (grid <= mid)
        yield label + '_high', mask & (grid > mid)



def run(source, name, mesh, field_mesh):
    with np.load(source) as saved:
        chi, pole, vector = saved['order'], complex(saved['pole']), saved['vector']
    sample = state(1, name)
    shape = Meridian(1, sample['coordinates'])
    core = Order(shape, mesh)
    fields = ActiveWave.install(load_for(shape, mesh, field_mesh, pole.real, (pole, vector)))

    def energy(c):
        # Restore the exact branch at each test. Local order variations can reactivate null amplitudes.
        fields.wave.target, fields.wave.seed = pole.real, (pole, vector.copy())
        eo, _ = core.energy_gradient(c)
        ef, _, info = fields.evaluate(c, 1., BIAS)
        return eo + ef, info

    ref, initial = energy(chi)
    mask = (chi <= 1e-12) & core.free
    assert mask.any()
    print(json.dumps({'source': str(source), 'sha256': hashlib.sha256(source.read_bytes()).hexdigest(),
                      'reference_energy': ref, 'empty_nodes': int(mask.sum()),
                      'zero_cells': initial['zero_cells'], 'inactive': initial['inactive_nodes']}), flush=True)
    for zone, where in directions(chi, mask):
        if not where.any():
            continue
        results = []
        for h in (5e-5, 1e-4):
            c = chi.copy()
            c[where] += h
            power, details = energy(c)
            results.append({'step': h, 'energy_slope': (power - ref) / h,
                            'overlap': details['mode_overlap'], 'pole_residual': details['pole_residual']})
        print(json.dumps({'zone': zone, 'sites': int(where.sum()), 'slopes': results}), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('source', type=Path)
    parser.add_argument('name', choices=['mid', 'tight'])
    parser.add_argument('mesh', type=int)
    parser.add_argument('field_mesh', type=int)
    args = parser.parse_args()
    run(args.source, args.name, args.mesh, args.field_mesh)
