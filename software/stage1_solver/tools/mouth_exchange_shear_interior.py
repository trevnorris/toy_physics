"""Same outgoing degree-four mode under an independent held-opening redistribution."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.sparse.linalg import eigsh

from mouth_exchange_shear import assemble, solve
from mouth_exchange_interior import InteriorShape


def run(mesh, base):
    guard=json.loads((base/'outcome.json').read_text())
    if guard['exitCode'] or guard['childOutcome']['exitCode'] or not guard['limitsVerified']:
        raise ValueError('Reference did not pass the guard')
    data=(base/'stdout').read_bytes()
    zero, = (row for row in map(json.loads,data.splitlines())
             if row['test']=='mode' and row['step']==0 and row['fade']==0 and row['mesh']==mesh and row['port']==3.)
    prior=complex(zero['frequency'],-zero['amplitude_loss'])
    # Inspect separation of the next *localized clamped seed*. The open brane
    # itself has a continuum, so this is not a total spectral gap.
    points,k,m,_=assemble(InteriorShape(0.),mesh,3.,0.)
    pocket=points[:,1]>0.
    values=eigsh(k[pocket][:,pocket],k=3,M=m[pocket][:,pocket],sigma=0.,
                 return_eigenvectors=False,v0=np.ones(sum(pocket)),tol=1e-10)
    print(json.dumps(dict(test='reference',mesh=mesh,frequency=prior.real,
                          reference_sha256=hashlib.sha256(data).hexdigest(),
                          pocket_seed_frequencies=np.sqrt(np.sort(values)).tolist())),flush=True)
    rows=[]
    for step in (-.005,.005):
        out,_,_=solve(InteriorShape(step),mesh,3.,0.,prior)
        print(json.dumps(dict(test='mode',step=step,mesh=mesh,**out)),flush=True)
        rows.append(out)
    slope=(rows[1]['frequency']-rows[0]['frequency'])/.01
    work=(rows[1]['stored_energy_slope']+rows[0]['stored_energy_slope'])/2.
    print(json.dumps(dict(test='second_direction',mesh=mesh,frequency_first=slope,
                          mean_virtual_energy_work=work)),flush=True)
    if abs(slope-work)>.002*max(.01,abs(slope)):
        raise ArithmeticError('Stored energy and spectral first-order slope disagree')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--mesh',type=int,required=True)
    p.add_argument('--base',type=Path,required=True)
    a=p.parse_args()
    run(a.mesh,a.base)
