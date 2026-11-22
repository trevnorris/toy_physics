"""Finite simultaneous shape/order work at the last resolved free-site state."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from mouth_poly_shallow_release import Trial,emit
from mouth_sealed_shape_relax import Meridian


def run(radius,mesh,path):
    parent = json.loads((path.parent/'outcome.json').read_bytes())
    if parent['exitCode'] or parent['childOutcome']['exitCode'] or not parent['limitsVerified']:
        raise ValueError('Parent readout not validated')
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    rows = [json.loads(line) for line in (path.parent/'stdout').read_text().splitlines()]
    end, = [row for row in rows if row.get('kind')=='material_endpoint' and row['snapshot_sha256']==digest]
    with np.load(path,allow_pickle=False) as f:
        if f['radius'].item()!=radius or f['mesh'].item()!=mesh:
            raise ValueError('Endpoint identity changed')
        shifts,chi,seed = f['shifts'],f['chi'],(complex(f['pole'].item()),f['vector'])
    trial = Trial(radius,mesh,bonded_core=True)
    def read(sh,order):
        val,g,own,loads,info,_ = trial.evaluate(sh,order,seed)
        if (not np.isfinite(val) or info['mode_overlap']<.94 or info['pole_residual']>1e-7
                or abs(info['current']-info['outflow'])>1e-6*info['current'] or info['current']<=0):
            raise ArithmeticError('Coupled variation lost resolved wave or mass balance')
        info['exit_radius']=Meridian(radius,trial.x+sh[:4]).exit
        return val,g,own,info
    en,g,owner,info = read(shifts,chi)
    if max(abs(info['energy'][k]-v) for k,v in end['energy'].items())>3e-7:
        raise ArithmeticError('Source state was not recovered')
    # A common physical interpolation in height and radius on both grids keeps
    # previously silent cells exactly silent. All six geometric works have the
    # same central step used by the independent signed witnesses.
    t = np.linspace(0,1,chi.shape[0])[:,None]
    z = np.linspace(0,1,chi.shape[1])[None,:]
    patch = chi*(1-chi)*np.exp(-((t-.7)**2+(z-.55)**2)/(2*.12**2))
    patch[~owner.free] = 0.
    slope = np.sum(g*patch)
    patch *= -np.sign(slope)
    grad = []
    for i in range(6):
        delta = np.zeros(6)
        delta[i] = .0005
        plus = read(shifts+delta,chi)[0]
        minus = read(shifts-delta,chi)[0]
        grad.append((plus-minus)/.001)
    grad = np.array(grad)
    direction = -grad/max(abs(grad))
    emit(kind='center',source_sha256=digest,source_status=end['status'],structural_gradient=grad.tolist(),
         direction=direction.tolist(),patch_slope_before_orientation=slope,joint_slope=grad@direction-abs(slope),**info)
    for h in (-.005,-.0025,.0025,.005):
        try:
            value,_,_,result = read(shifts+h*direction,chi+h*patch)
        except (ArithmeticError,ValueError,RuntimeError) as err:
            emit(kind='unresolved_probe',step=h,reason=str(err))
            continue
        emit(kind='joint',step=h,delta_energy=value-en,**result)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--radius',type=int,choices=(1,2),required=True)
    p.add_argument('--mesh',type=int,choices=(48,72),required=True)
    p.add_argument('--snapshot',type=Path,required=True)
    a=p.parse_args()
    run(a.radius,a.mesh,a.snapshot)
