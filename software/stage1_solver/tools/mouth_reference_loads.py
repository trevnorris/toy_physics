"""Resolve the two sides of each archived force probe before adding crystal prestress."""
import hashlib
import json
from math import fsum, isclose
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
SOURCE = ROOT / 'research/pde_ledger_v3/_measurements/mouth_unified_strip_hold.json'
ENERGIES = ('order_inside', 'order_outside', 'carrier_two_sided', 'stored_wave')


def signed_work(first, last):
    x0, left = first
    x1, right = last
    if x1 <= x0:
        raise ValueError('Reversed or degenerate coordinate interval')
    return {component: (left[component] - right[component]) / (x1 - x0)
            for component in ENERGIES}


def main():
    source = json.loads(SOURCE.read_text())['runs']['shape']
    raw = '\n'.join(source['literal_stdout']) + '\n'
    digest = hashlib.sha256(raw.encode()).hexdigest()
    if (digest != source['sha256_stdout'] or source['guard_and_child_exit'] != 0
            or source['swap_cap_bytes'] != 0):
        raise ValueError('Unverified matched-shape input')
    rows = [json.loads(line) for line in raw.splitlines()]
    results = []
    for radius in (1, 2):
        for mesh in (48, 72):
            case = f'r{radius}_broad_{mesh}'
            for axis in ('radial', 'axial'):
                pairs = [r for r in rows if (r['case'], r['axis']) == (case, axis)]
                if sorted(r['offset'] for r in pairs) != [.001, .002]:
                    raise ValueError('Missing signed offsets')
                if len({(r['input_wave_sha256'], r['input_order_sha256']) for r in pairs}) != 1:
                    raise ValueError('Inputs differ between steps')
                endpoints = {}
                for row in pairs:
                    h = row['offset']
                    for sign, leg in zip((-1, 1), row['legs'], strict=True):
                        offset = sign * h
                        x = offset if axis == 'radial' else leg['depth']
                        endpoints[offset] = (x, leg)
                    center = signed_work(endpoints[-h], endpoints[h])
                    if not isclose(fsum(center.values()), row['held_zero_crystal_prestress_total'],
                                   rel_tol=1e-11):
                        raise ValueError('Native center-work control failed')
                for a, b in ((-.002, -.001), (.001, .002)):
                    force = signed_work(endpoints[a], endpoints[b])
                    results.append(dict(radius=radius, mesh=mesh, axis=axis, offsets=[a,b],
                                        coordinate_interval=[endpoints[a][0], endpoints[b][0]],
                                        components=force, force_without_crystal=fsum(force.values())))
    print(json.dumps(dict(source_stdout_sha256=digest, secants=results), indent=2))


if __name__ == '__main__':
    main()
