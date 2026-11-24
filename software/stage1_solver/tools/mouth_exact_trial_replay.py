"""Reconstruct only the original full shape proposal following radius-four update 10.

The archived trajectory contains the accepted shape, shape derivatives and every
outer rematch rejection. Recover the selected inverse shape Hessian by matching
each intermediate accepted shape bit for bit, then run the producer's Armijo
selection and field match for the one missing trial. The archived gate is
observed before any change to its verdict. The candidate is never promoted here.
"""
import argparse
import hashlib
import json
from pathlib import Path
import sys
from types import SimpleNamespace

import numpy as np


CP10 = '51727c2d48eea199c29bacbfd4bf07b7b11a3e997a737ebcc8c2f69925cd112b'


def read(path):
    return json.loads(path.read_bytes())


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def shape_history(run):
    folder = run/'coarse'
    accepted, residual, rejects = {}, {}, {}
    tail = None
    for line in (folder/'stdout').read_text().splitlines():
        row = json.loads(line)
        if row['kind']=='accepted':
            accepted[row['iteration']] = row
        if row['kind']=='residuals':
            tail = row['iteration']
            residual[tail] = row
            rejects[tail] = []
        if row['kind']=='rejected_shape_rematch':
            rejects[tail].append(row)
    if accepted[10]['sha256'] != CP10 or rejects[10][0]['scale']!=1.:
        raise ValueError('Incorrect archived trial')
    g_old = x_old = metric = None
    history = []
    for i in range(7,11):
        record = accepted[i]
        if digest(folder/record['file'])!=record['sha256']:
            raise ValueError(f'Checkpoint hash mismatch at {i}')
        x = np.array(record['x'])
        g = -np.array(residual[i]['shape_force'])
        diag = np.array(residual[i]['fixed_material_shape_curvatures'])
        diagonal = np.diag(1./np.maximum(abs(diag),1.))
        if metric is None:
            metric = diagonal
        if g_old is not None:
            dx, dg = x-x_old, g-g_old
            prod = dx@dg
            if prod>1.e-10*np.linalg.norm(dx)*np.linalg.norm(dg):
                b = np.eye(len(dx))-np.outer(dx,dg)/prod
                metric = b@metric@b.T+np.outer(dx,dx)/prod
        if i==10:
            return metric,history,accepted[i],residual[i],rejects[i][0]
        scales = [r['scale'] for r in rejects[i]]
        if scales != [2.**-j for j in range(len(scales))]:
            raise ValueError('Noncontiguous rematch backtracking history')
        rematch_scale = 2.**-len(scales)
        next_x = np.array(accepted[i+1]['x'])
        matches = []
        for label,candidate in [('carried',metric),('diagonal',diagonal)]:
            direction = -(candidate@g)
            direction *= min(1., .035/max(1.e-20,abs(direction).max()))
            for k in range(26):
                if np.array_equal(x+rematch_scale*(direction*2.**-k),next_x):
                    matches.append((label,k,candidate))
        if not matches or any(not np.array_equal(m[2],matches[0][2]) for m in matches):
            raise ValueError(f'Archived accepted trial {i+1} does not fix the shape metric')
        label,k,metric = matches[0]
        history.append(dict(iteration=i,metric=label,armijo_scale=2.**-k,
                            rematch_scale=rematch_scale,bitwise_next_x=True))
        g_old,x_old = g,x


def main(run, out):
    from mouth_continuous_shape_relax import ContinuousRelease
    from mouth_joined_front_audit import Front
    from mouth_last_rings_fluid import GRIDS
    from mouth_melted_collar_survival import MeltTrial, ReorderingUnresolved
    from mouth_whole_r4_relax import Release, atomic_json
    from mouth_exchange_carrier import BIAS

    out.mkdir(parents=True,exist_ok=True)
    if any((out/name).exists() for name in ('proposal.json','trial.npz','replay.json')):
        raise FileExistsError('Use a new output directory for a distinct diagnostic')
    for relative,sha in read(run/'sources.json').items():
        if digest(run/'pin'/relative)!=sha:
            raise ValueError(f'Producer pin changed at {relative}')
    metric,history,accepted,residual,expected = shape_history(run)
    source = run/'coarse'/accepted['file']
    cpcmd = read(run/'coarse/invocation.json')['command']
    origin = Path(cpcmd[-3])/'coarse'
    sourcecmd = read(origin/'invocation.json')['command']
    previous = Path(sourcecmd[-5])/'coarse'
    proof = read(Path(sourcecmd[-3]))['results']['coarse']
    prep = Path(read(previous/'invocation.json')['command'][5])
    if digest(prep)!=proof['preparation_sha256'] or digest(previous/'endpoint.json')!=proof['endpoint_sha256'] or BIAS!=.4510084282857874:
        raise ValueError('Original radius-four model inputs differ')
    with np.load(source,allow_pickle=False) as saved:
        q,x,u = [saved[key].copy() for key in ('q','x','field')]
        if int(saved['radius']) != 4 or int(saved['winding']) != 1 or str(saved['state_sha'])!=digest(prep):
            raise ValueError('Saved preparation mismatch')
    class Capture(MeltTrial):
        def check_settled_fields(self, q, u, gradient):
            try:
                super().check_settled_fields(q,u,gradient)
            except ReorderingUnresolved as exc:
                self.rejection = str(exc)
            else:
                raise ValueError('Original full-step rejection failed to reproduce')
            self.final_gradient = gradient.copy()
            return {}
    class Shape:
        build = ContinuousRelease.build
        shape_step = Release.shape_step
        trial_type = Capture
    owner=Shape()
    owner.base = Front(prep,GRIDS['double'],1)
    owner.args = SimpleNamespace(state=prep)
    owner.reference_depth = read(previous/'endpoint.json')['depth']
    owner.discretization=proof['solid_discretization']
    owner.x,owner.shape_metric = x,metric
    g = -np.array(residual['shape_force'])
    dx,selected = owner.shape_step(g,np.array(residual['fixed_material_shape_curvatures']),q,u,accepted['energy'])
    if dx is None:
        raise ValueError('Pinned Armijo selector did not find the full prior proposal')
    proposal=dict(input_sha=CP10,phase_drop=BIAS,shape_history=history,previous_x=x.tolist(),
                  full_dx=dx.tolist(),x=(x+dx).tolist(),metric=selected.tolist(),expected_rejection=expected)
    atomic_json(out/'proposal.json',proposal)
    print(json.dumps(proposal,allow_nan=False),flush=True)
    trial=owner.build(x+dx,True)
    new_u,e,gq,load,info=trial.match(q,u)
    node=51603
    joined=trial.pull.T@trial.final_gradient[2]
    state_file=out/'trial.npz'
    tmp=out/'trial.partial.npz'
    np.savez_compressed(tmp,field=new_u,gradient=trial.final_gradient,
                        joined_order_gradient=joined,q=q,x=trial.x,solid=trial.b.unpack(q),
                        gq=gq,load=load,energy=e,shape_params=trial.shape.params,
                        shape_extra=trial.shape.extra,state_sha=digest(prep),radius=4,winding=1)
    tmp.replace(state_file)
    result=dict(proposal_sha256=digest(out/'proposal.json'),snapshot_sha256=digest(state_file),
                archived=expected['detail'],replayed=trial.rejection,
                identical_rejection=(expected['detail']==trial.rejection),
                order=float(new_u[2,trial.first[node]]),lower_bound_order_derivative=float(joined[node]),
                field=info['field'],energy=e,solid_order_work_error=info['exchange_error'],
                full_dx=dx.tolist(),entrance_radius=trial.shape.neck,exit_radius=trial.shape.exit,depth=trial.shape.length,
                current=trial.measure(new_u,load))
    atomic_json(out/'replay.json',result)
    print(json.dumps({k:result[k] for k in ('identical_rejection','replayed','order','lower_bound_order_derivative','energy')},allow_nan=False),flush=True)
    if not result['identical_rejection']:
        raise ValueError('Original rejection differs, do not admit this archived full-step replay')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('run',type=Path)
    p.add_argument('output',type=Path)
    args=p.parse_args()
    run=args.run.resolve()
    sys.path.insert(0,str(run/'pin/software/stage1_solver/tools'))
    main(run,args.output.resolve())
