"""Check compatibility of a smooth collar metric with the stored boundary traces."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.interpolate import CubicSpline


def audit(root):
    prep = root/'supercore_joined_front01/adopt309/state.npz'
    step = root/'r4_segmented_held01/primary/step_00166.npz'
    def sha(path):
        return hashlib.sha256(path.read_bytes()).hexdigest()
    if sha(step) != json.loads(step.with_suffix('.json').read_text())['sha256']:
        raise ValueError('Step 166 checksum mismatch')
    with np.load(root/'r4_melted_survival01/coarse/state_00010.npz',allow_pickle=False) as old:
        if sha(prep) != str(old['state_sha']):
            raise ValueError('The saved spatial solver uses a different registration')
    with np.load(prep,allow_pickle=False) as reg:
        if int(reg['radius'])!=4:
            raise ValueError('Only the held radius-four registry is admitted')
        points=reg['labels'].copy()
    points[...,0] += 4
    t=np.linspace(0,1,len(points))
    eta=np.linspace(-.5,.5,points.shape[1])
    i=points.shape[1]//2
    boundaries=dict(weld=(t,points[:,-1]),lateral_wall=(t,points[:,i]),
                    entrance=(eta[i:],points[0,i:]),exit=(eta[i:],points[-1,i:]))
    out={}
    for name,(q,z) in boundaries.items():
        d=np.diff(z,axis=0)/np.diff(q)[:,None]
        u=d/np.linalg.norm(d,axis=1)[:,None]
        jump=np.linalg.norm(np.diff(d,axis=0),axis=1)
        k=int(np.argmax(jump))+1
        p=q[:-1,None]+np.diff(q)[:,None]*np.linspace(0,1,33)
        spline=CubicSpline(q,z,axis=0,bc_type=((1,d[0]),(1,d[-1])))
        baseline=z[:-1,None,:]+(p-q[:-1,None])[...,None]*d[:,None,:]
        placed=spline(p)
        out[name]=dict(number_segments=len(d),anchor=z[k].tolist(),label=float(q[k]),
                       derivative_left=d[k-1].tolist(),derivative_right=d[k].tolist(),
                       max_jump=float(jump[k-1]),minimax_tangent_discrepancy=float(jump[k-1]/2),
                       max_unit_tangent_jump=float(np.max(np.linalg.norm(np.diff(u,axis=0),axis=1))),
                       stored_curve_is_straight=bool(np.max(np.linalg.norm(np.diff(u,axis=0),axis=1))<1.e-14),
                       fixed_label_position_delta_clamped_cubic=float(np.max(np.linalg.norm(placed-baseline,axis=-1))),
                       max_anchor_position_delta_clamped_cubic=float(np.max(np.abs(spline(q)-z))),
                       parameter_direction='eta' if name in ('entrance','exit') else 'normalized axial label')
    return dict(step166_sha256=sha(step),registry_sha256=sha(prep),boundaries=out,
                chart='S=4+stored radial label; T=stored axial label, from the production natural reference',
                exact_parameterized_C1_compatible=all(o['max_jump']<1.e-14 for o in out.values()))


if __name__=='__main__':
    a=argparse.ArgumentParser(description=__doc__)
    a.add_argument('root',type=Path)
    a.add_argument('destination',type=Path)
    opt=a.parse_args()
    report=audit(opt.root)
    opt.destination.parent.mkdir(parents=True,exist_ok=True)
    opt.destination.write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    print(json.dumps(report,allow_nan=False))
