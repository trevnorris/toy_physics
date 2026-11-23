"""Rematch the phase, amplitude and order on a nested carrier mesh at the held R4 endpoint."""
import argparse
import gc
import hashlib
import json
from pathlib import Path
import sys
from types import SimpleNamespace

import numpy as np


def digest(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def read(path):
    return json.loads(path.read_text())


def pin(run, edges):
    manifest = run/'sources.json'
    references = read(edges)
    if digest(manifest) != references['pin_manifest_sha256']:
        raise ValueError('Endpoint pin manifest changed')
    for name, sha in read(manifest).items():
        if digest(run/'pin'/name) != sha:
            raise ValueError(f'Changed pinned file {name}')
    sys.path.insert(0, str((run/'pin/software/stage1_solver/tools').resolve()))
    return references['results']


def endpoint(run, mesh, reference):
    from mouth_rim_pinned_comparison import validate_guard
    job = run/mesh
    validate_guard(job)
    report = read(job/'endpoint.json')
    checkpoint = job/report['file']
    controls = job/'endpoint_controls.npz'
    command = read(job/'invocation.json')['command']
    state = Path(command[command.index('-u')+2])
    for path, key in ((state, 'preparation_sha256'), (checkpoint, 'state_sha256'),
                      (job/'endpoint.json', 'endpoint_sha256'), (controls, 'controls_sha256')):
        if digest(path) != reference[key]:
            raise ValueError(f'Unexpected saved input {path}')
    with np.load(checkpoint) as f, np.load(controls) as c:
        return state, f['x'], f['q'], f['field'], c['rest'], report


def refine(args):
    ref = pin(args.run, args.edges)['fine']
    from mouth_disk_mesh_rematch import remap
    from mouth_free_shape_solid import Solid
    from mouth_joined_front_audit import Front
    from mouth_last_rings_fluid import GRIDS
    from mouth_mobile_joined_fields import JointFields

    state, x, q, wound, rest, previous = endpoint(args.run, 'fine', ref)
    base = Front(state, GRIDS['double'], 1)
    original = JointFields(base, state, x, precondition=False)
    # Only the interpolation topology and receiver are needed from this mesh.
    old = SimpleNamespace(rr=original.model.rr, a=original.model.a, z=original.model.z,
                          ids=original.model.ids, size=original.model.size)
    del original
    gc.collect()
    new = JointFields(base, state, x, grid='local')
    m = new.model
    if any(np.max(np.min(abs(a[:, None]-b), axis=0)) > 2e-13
           for a, b in zip((*m.rr, *m.z, m.a), (*old.rr, *old.z, old.a))):
        raise ArithmeticError('Original nodes are missing from the refined chart')
    s = ref['solid_discretization']
    solid = Solid(state, m.shape, m.polar_shift, s['n'], s['k'], s['angles'], m,
                  degree=s['degree'], subdivisions=s['subdivisions'])
    es, _, load, info = solid.evaluate(q, gradient=False)
    args.saved.mkdir(parents=True, exist_ok=False)
    report = dict(command=sys.argv, script_sha256=digest(Path(__file__)),
                  source_manifest_sha256=digest(args.run/'sources.json'),
                  edge_readback_sha256=digest(args.edges), original_sha256=ref['state_sha256'],
                  controls_sha256=ref['controls_sha256'], state_sha256=ref['preparation_sha256'],
                  solid_discretization=s, coordinates=x.tolist(), solid_reference_energy=es,
                  deposition_sum=float(load.sum()), deposition_l1=float(abs(load).sum()), solid=info,
                  grid=GRIDS['local'], comparisons={})
    for name, seed in [('wound', wound), ('rest', rest)]:
        remapped = remap(old, m, seed)
        u, sol = new.solve(remapped, load, rest=name=='rest')
        path = args.saved/(name+'.npz')
        np.savez_compressed(path, field=u)
        measured = new.measure(u, load)
        sol['independent_residual'] = new.residual(u, load)
        report['comparisons'][name] = dict(sha256=digest(path), original=previous['controls'][name],
                                            rematched=measured, solver=sol)
        (args.saved/'readback.json').write_text(json.dumps(report, allow_nan=False)+'\n')
        if (sol['independent_residual'] >= 2e-7 or measured['port_error'] >= 3e-5
                or abs(measured['disk_minus_feed_current']) >= 3e-5
                or min(measured['nodal_order_range_in_band']) < -2e-12
                or max(measured['nodal_order_range_in_band']) > 1+2e-12):
            raise ArithmeticError(f'Unresolved held {name} field rematch')
    args.output.write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
    print(json.dumps(dict(output_sha256=digest(args.output), field_sha256={
        name: row['sha256'] for name, row in report['comparisons'].items()})))


if __name__=='__main__':
    p = argparse.ArgumentParser(description=__doc__)
    for key in ('run', 'edges', 'saved', 'output'):
        p.add_argument(key, type=Path)
    refine(p.parse_args())
