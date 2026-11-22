"""Read normal and tangential pressure on the prospective brane plane from saved bulk stresses."""
import argparse
import hashlib
import json
from math import isclose
from pathlib import Path


def read_pressures(root):
    manifest = json.loads((root / 'ZERO_RUNS.json').read_bytes())
    result = {}
    for name, checksum in manifest['results'].items():
        raw = (root / 'runs' / (name + '.json')).read_bytes()
        assert hashlib.sha256(raw).hexdigest() == checksum
        run = json.loads(raw)
        assert run['strain'] == run['pressure_external'] == 0 and run['axis_index'] == 1
        blocks = run['blocks'][-4:]
        count = sum(b['steps'] for b in blocks)
        diag = [-sum(b['steps'] * b['mean']['stress_tensile'][i][i] for b in blocks) / count
                for i in range(4)]
        pressure = sum(b['steps'] * b['mean']['pressure'] for b in blocks) / count
        normal, in_plane = diag[3], sum(diag[:3])/3
        assert isclose((normal+3*in_plane)/4, pressure, rel_tol=1e-10, abs_tol=1e-11)
        result[name] = dict(scalar=pressure, normal=normal, mean_in_plane=in_plane)
        print(json.dumps(dict(sample=name, **result[name]), allow_nan=False))
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('root', type=Path)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    args.output.write_text(json.dumps(read_pressures(args.root), indent=2, allow_nan=False) + '\n')
