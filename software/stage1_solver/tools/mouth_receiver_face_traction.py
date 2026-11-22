"""Measure the mean traction on the periodic face that would be used to join two parent cells."""
import argparse
from hashlib import sha256
import json
from pathlib import Path

import numpy as np


def face_frame(box):
    # Preserve the directions of the three ordered tangents and point n toward row four.
    basis = []
    for vector in box:
        component = vector.copy()
        for axis in basis:
            component -= np.dot(axis, component) * axis
        basis.append(component / np.linalg.norm(component))
    frame = np.array(basis)
    np.testing.assert_allclose(frame @ frame.T, np.eye(4), atol=1e-12)
    return frame


def read_faces(batch):
    hashes = json.loads((batch / 'ZERO_RUNS.json').read_text())['results']
    result = {}
    for key, checksum in sorted(hashes.items()):
        data = (batch / 'runs' / (key+'.json')).read_bytes()
        assert sha256(data).hexdigest() == checksum
        r = json.loads(data)
        assert r['strain'] == r['pressure_external'] == 0.
        frame = face_frame(np.asarray(r['prepared']['box']))
        n = frame[-1]
        blocks = r['blocks'][-4:]
        for block in blocks:
            np.testing.assert_allclose(face_frame(np.asarray(block['state']['box'])), frame, atol=1.e-12)
        stress = np.average([b['mean']['stress_tensile'] for b in blocks],
                            weights=[b['steps'] for b in blocks], axis=0)
        components = frame @ stress @ n
        np.testing.assert_allclose(np.linalg.norm(components), np.linalg.norm(stress @ n))
        result[key] = dict(normal_lab=n.tolist(),
                           angle_deg_from_w=float(np.degrees(np.arccos(np.clip(n[-1],-1.,1.)))),
                           face_pressure=float(-components[-1]), tangential_stress=components[:3].tolist(),
                           tangential_magnitude=float(np.linalg.norm(components[:3])))
        print(json.dumps(dict(sample=key, **result[key]), allow_nan=False))
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('batch', type=Path)
    parser.add_argument('out', type=Path)
    args = parser.parse_args()
    args.out.write_text(json.dumps(read_faces(args.batch), indent=2, allow_nan=False)+'\n')
