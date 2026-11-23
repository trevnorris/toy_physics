"""Compare quadrature summation orders on the two saved, stationary R4 fields."""
import argparse
import json
import math
from pathlib import Path
import sys

import numpy as np


def totals(t, u, load):
    from mouth_breathing_mobile_patch import KAPPA, ORDER_INERTIA, WELL
    m=t.model
    v, phase, c = u
    wt, f = m.hq, m.D[:,4]
    if np.any(m.ref) or np.any(m.D[:,10]):
        raise ValueError("The saved state requires absolute order and one material deposition")
    test=m.V@v
    n=(f*test)**2
    y=test**2
    chi=m.V@c
    gv,gp,gc=(np.column_stack([g@a for g in m.G]) for a in u)
    square=lambda x: np.einsum('ij,ij->i',x,x)
    dot=np.einsum('ij,ij->i',gp,gc)
    den=n+ORDER_INERTIA*square(gc)
    field_phase=.5*np.einsum('ij,ij->i',gp,n[:,None]*(gp-(ORDER_INERTIA*dot/den)[:,None]*gc))
    phase_kin=.5*n*square(gp)
    coupling=-.5*ORDER_INERTIA*n*(dot**2/den)
    grad=.5*KAPPA*square(gc)
    well=WELL*(chi*(1.-chi))**2
    kinetic=.5*f**2*square(gv)
    density=f**10*(y-1.)**2*(y**3+2*y**2+3*y+4)/20.
    combined=kinetic+density+field_phase+(grad+well)
    def precise(a):
        return (wt.astype(np.longdouble)*a).sum(dtype=np.longdouble)
    boundary=.5*((v-1)@(m.B@(v+1)))
    solid=load@c
    grouped=float(precise(combined)+boundary+solid)
    parts=float(sum(precise(a) for a in (kinetic,density,phase_kin,coupling,grad,well))+boundary+solid)
    native=(wt@combined+boundary)+solid
    return dict(native_grouped=native,precise_grouped=grouped,precise_components=parts,
                max_local_phase_identity=float(abs(phase_kin+coupling-field_phase).max()))


def main(held, original, ring, destination):
    # Validate the same material, geometry and fields used in both archived sums.
    sys.path[:0]=[str(held/'pin'),str(original/'pin/software/stage1_solver/tools')]
    from mouth_continuous_evidence import held_receipt
    from mouth_continuous_order import ContinuousTrial
    from mouth_joined_front_audit import Front
    from mouth_last_rings_fluid import GRIDS
    from mouth_rim_pinned_comparison import digest,read
    receipt=held_receipt(held/'held')
    rec=read(held/'saved/record.json')
    evidence=read(ring)
    state=Path(read(original/'fine/invocation.json')['command'][5])
    if digest(state)!=rec['state'] or np.finfo(np.longdouble).eps>=np.finfo(float).eps:
        raise ValueError("Incompatible input or accumulator precision")
    report=dict(source=digest(Path(__file__)),ring_readout=digest(ring),held_receipt=digest(receipt),
                longdouble_mantissa=np.finfo(np.longdouble).nmant, cases={})
    for name in ('rest','wound'):
        path=held/'saved'/(name+'.npz')
        expected=rec['solutions'][name]
        if digest(path)!=expected['sha256']: raise ValueError(f'Changed field: {name}')
        with np.load(path) as data:
            if name=='rest':
                trial=ContinuousTrial(Front(state,GRIDS['double'],1),state,data['x'],**rec['solid'],precondition=False)
            row=totals(trial,data['field'],data['load'])
        archived=evidence['cases'][name]
        row.update(saved=expected['energy'],separate=float(math.fsum(archived['energy_terms'].values())),
                   field_sha256=expected['sha256'])
        row['precise_difference']=row['precise_components']-row['precise_grouped']
        row['saved_minus_native']=row['saved']-row['native_grouped']
        report['cases'][name]=row
    destination.write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    for row in report['cases'].values():
        if abs(row['saved_minus_native'])>2e-8 or abs(row['precise_difference'])>2e-8:
            raise ArithmeticError("Unresolved native or independently accumulated energy")


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for arg in ('held','original','ring','destination'):
        p.add_argument(arg,type=Path)
    a=p.parse_args()
    main(a.held,a.original,a.ring,a.destination)
