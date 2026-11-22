"""Equilibrium of the fraction field on its declared 0 to 1 interval."""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.optimize import minimize

from mouth_live_waist_order import BIAS, Order, Meridian, emit, load_for, state


def bounded(order, chi, loads, action):
    chi = chi.copy()
    chi[~order.free] = order.seed[~order.free]
    scale = np.sqrt((2*order.stiff.diagonal()+1.5*order.lump)[order.flatfree])
    def sources(x):
        chi[order.free] = x/scale
        eo, go = order.energy_gradient(chi)
        ef, gf, _ = loads.evaluate(chi,action=action,bias=BIAS)
        return eo+ef, (go+gf)[order.free]/scale
    x = np.maximum(0., np.minimum(1.,chi[order.free]))*scale
    answer = minimize(sources,x,jac=True,method='L-BFGS-B',bounds=list(zip(np.zeros_like(scale),scale)),
                      options={'gtol':6e-7,'ftol':1e-13,'maxiter':1500,'maxls':30,'maxcor':18})
    chi[order.free] = answer.x/scale
    eo, go = order.energy_gradient(chi)
    ef, gf, info = loads.evaluate(chi, action=action, bias=BIAS)
    g = (go+gf)[order.free]
    values = chi[order.free]
    bed = (values<1e-10)&(g>0); top = (values>1-1e-10)&(g<0)
    proj = np.where(bed|top, 0., g)/scale
    norm = max(abs(proj))
    m = {k:float(v) for k,v in info.items() if np.isscalar(v) and isinstance(v,(int,float,np.number))}
    m.update(total_inside_and_drives=eo+ef,order_energy=eo,axis_order=order.axis(chi),
             min_order=float(chi.min()),max_order=float(chi.max()),iterations=answer.nit,
             order_residual=norm, waist_radius=float(order.shape.radius(order.shape.length*order.t_waist)[0]),
             axial_depth=order.shape.length, fraction_at_zero=int(sum(bed)), fraction_at_one=int(sum(top)),
             outward_reaction_at_zero=sum(g[bed]),outward_reaction_at_one=-sum(g[top]),
             convergence_message=str(answer.message))
    if answer.success and norm<2e-5:
        return chi, m
    emit(kind='unresolved',action=action,**m)
    raise ArithmeticError(f'Bounded order balance unresolved at scaled KKT residual {norm}')


def run(a):
    s = state(a.radius,a.point)
    params=np.array(s['coordinates'])
    shape=Meridian(a.radius,params); order=Order(shape,a.mesh)
    loads=load_for(shape,a.mesh,a.field_mesh,s['wave_frequency'])
    chi, flow=order.condition(order.seed,loads,0.,BIAS)
    emit(kind='flow',**flow)
    for action in (.125,.25,.5,.75,1.):
        chi,m=bounded(order,chi,loads,action)
        emit(kind='coupled',action=action,**m)
    np.savez_compressed(a.save,order=chi,frequency=m['wave_frequency'],
                        vector=loads.wave.seed[1],pole=loads.wave.seed[0])
    for coord in ('waist','depth'):
        for h in (-.002,-.001,.001,.002):
            p=params.copy(); index=2 if coord=='waist' else 0
            p[index] += h/(a.radius if coord=='waist' else 1.)
            q=Meridian(a.radius,p)
            fields=load_for(q,a.mesh,a.field_mesh,s['wave_frequency'],loads.wave.seed)
            _,shifted=bounded(Order(q,a.mesh),chi,fields,1.)
            emit(kind='signed',coordinate=coord,offset=h,**shifted)


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--mesh',type=int,required=True)
    ap.add_argument('--field-mesh',type=int,required=True)
    ap.add_argument('--radius',type=int,choices=(1,2),required=True)
    ap.add_argument('--point',choices=('broad','mid','tight'),required=True)
    ap.add_argument('--save',type=Path,required=True)
    run(ap.parse_args())
