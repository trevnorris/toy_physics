"""Locate the material residual of a saved, fixed R4 wall and test its angular dependence."""
import argparse
import json
from pathlib import Path

import numpy as np

from mouth_bonded_front_adjoint import set_front
from mouth_fixed_wall_angular import AngularBonded
from mouth_joined_front_audit import Front
from mouth_joint_polar_shape import Trial
from mouth_supercore_joined_front import check_guard,emit,file_sha


def scan(state, snapshot):
    check_guard(state)
    folder=snapshot.parent
    outcome=json.loads((folder/'outcome.json').read_text())
    limits=json.loads((folder/'limit-validation.json').read_text())
    rows=[json.loads(s) for s in (folder/'stdout').read_text().splitlines() if s.startswith('{')]
    if not outcome['limitsVerified'] or not limits['verified'] or int(limits['actual']['memory.swap.max']):
        raise ValueError('Unverified guarded path')
    with np.load(snapshot,allow_pickle=False) as z:
        x,d,field,s=z['x'],z['displacement'],z['field'],float(z['internal_amplitude'])
        data=json.loads(z['readout'].item())
        if z['state_sha'].item()!=file_sha(state) or z['input_sha'].item()!=rows[0]['snapshot']:
            raise ValueError('Source does not match the recorded preparation')
    last=next(row for row in reversed(rows) if row['kind']=='fixed_wall_sweep')
    if last['sweep']!=data['sweep'] or last['energy']!=data['energy']:
        raise ValueError('Expected the last fully saved sweep')
    base=Front(state,(13,10,14),1)
    n,k=(d.shape[1]-1)//3,d.shape[0]-1
    trial=Trial(base,state,x,n,k,False)
    b=AngularBonded(state,n,k,shape=trial.shape,epsilon=x[5])
    set_front(b,trial.model,field)
    b.internal_amplitude=s
    q=d[b.mask]*b.scale
    b.energy(q)
    g=np.where(b.mask,b.native,0.).copy()
    index=np.unravel_index(np.argmax(abs(g)), g.shape)
    iz,ix,axis=index
    peak=float(abs(g[index]))
    if abs(peak-data['native_meridional_residual'])>1e-9:
        raise ArithmeticError('Cannot recover the saved native material residual')
    b.internal_amplitude=0.
    b.energy(q)
    z0=b.t[iz]*trial.shape.length
    r0=b.x[ix]*trial.shape.radius(z0)[0]
    site=[r0+d[iz,ix,0],z0+d[iz,ix,1]]
    emit(kind='saved_peak',snapshot=file_sha(snapshot),source=file_sha(state),path_exit=outcome['exitCode'],
         last_complete_sweep=data['sweep'], grid=[3*n,k], angular_amplitude=s,
         native_force_signed=float(-g[index]),peak=peak,angular_frozen_zero_force=float(-b.native[index]),
         component=('radial','axial')[axis], depth_fraction=float(b.t[iz]),
         radius_fraction=float(b.x[ix]), deformed_meridional_site=site,
         actual_depth_fraction=site[1]/trial.shape.length,angular_depth_band=[.06,.62],
         radial_max=float(abs(g[:,:,0]).max()),axial_max=float(abs(g[:,:,1]).max()))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('state',type=Path)
    parser.add_argument('snapshot',type=Path)
    args=parser.parse_args()
    scan(args.state,args.snapshot)
