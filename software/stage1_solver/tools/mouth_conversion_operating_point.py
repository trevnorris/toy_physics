"""Use measured order position to trace the neutral branch through its drive turning point.

Position parameterizes solutions. The original residual vanishes at all free
order nodes; no extra anchor is installed on the measured level set.
"""
import argparse
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
from scipy.sparse.linalg import splu

from mouth_conversion_supply import SuppliedOrder
from mouth_supplied_contact_work import verified
from mouth_carrier_lateral_response import prepare
from mouth_empty_bend_stiffness import quadrature
from mouth_curved_phase import quadratic
from mouth_order_wave_pulse import digest, emit


def main(a):
    records=verified(a.source)
    setup=next(row for row in records if row['test']=='calibration')
    source=Path(records[0]['source']).parent
    dom,psi=prepare(SimpleNamespace(source=source))
    start=next(row for row in records if row.get('multiple')==.8)
    file=Path(start['path'])
    if digest(file)!=start['sha256']:
        raise ValueError('Starting order changed')
    with np.load(file) as data:
        x=data['chi']
    model=SuppliedOrder(dom,psi,setup['width'],setup['mobility'])
    pg,vg,cells=dom.blocks[2]
    h=start['head']
    _,_,dr,dw=quadrature(dom)
    f=dom.N@psi
    ur,uw=(np.imag(g@(psi-1)/f) for g in (dr,dw))
    for z in (1.29,1.28,1.275,1.27,1.265,1.26,1.255,1.25,1.24,1.23,1.22):
        i=np.searchsorted(pg,z)-1
        nodes=cells[i*(len(vg)-1),[0,3,6]]
        trace=quadratic(2*(z-pg[i])/(pg[i+1]-pg[i])-1)[0]
        for it in range(24):
            c=(dom.N@x)/model.quad_sum
            d=(dom.N@(1-x))/model.quad_sum
            rhs=model.Nf.T@(model.wt*setup['mobility']*model.density*6*c*d)
            residual,jac=model.measure(x)
            residual+=h*rhs
            jac+=model.mass_weight(model.wt*setup['mobility']*h*model.density*6*(d-c))
            factor=splu(jac.tocsc())
            updates=np.zeros((len(x),2))
            updates[dom.free]=factor.solve(np.column_stack([-residual,rhs]))
            u,v=updates.T
            dh=(trace@(x[nodes]+u[nodes])-.5)/(trace@v[nodes])
            step=u-dh*v
            x+=step
            h+=dh
            if max(np.max(abs(step)),abs(dh))<2e-11:
                break
        else:
            raise ArithmeticError('Could not resolve position-parametrized neutral equilibrium')
        # Negative curvature on this explicit free variation proves instability.
        direction=v[dom.free]
        rq=(direction@(jac@direction))/(direction@(model.mass@direction))
        rate=ur*(dr@x)+uw*(dw@x)
        loss=dom.weight@(rate**2)/setup['mobility']
        emit(test='branch',mesh=setup['mesh'],position=z,head=h,head_power=h*setup['intake'],
             loss=loss,unpaid_loss=loss-h*setup['intake'],free_variation_curvature=rq,
             free_residual=np.max(abs(residual)),level_error=abs(trace@x[nodes]-.5),iterations=it+1)
        if max(np.max(abs(residual)),abs(trace@x[nodes]-.5))>1e-8:
            raise ArithmeticError('Constrained residual larger than numerical tolerance')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source',type=Path,required=True)
    main(p.parse_args())
