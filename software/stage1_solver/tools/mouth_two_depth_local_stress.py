"""Independent in-place renewal check using stress as the memory variable.

For spatially fixed renewal, S mdot = L bdot - B G m/tau. Rest offsets
which produce no stress drop out. The positive-rate search excludes zero.
"""
import argparse

import numpy as np
from scipy import linalg, sparse
from scipy.optimize import brentq
from scipy.sparse.linalg import splu

from mouth_bend_renewal import ROOT, volume_cost, wall_data
from mouth_bend_static_shear import ANGULAR, archive_inputs
from mouth_order_wave_pulse import emit
from mouth_two_depth_renewal import CoupledMouth, DEPTHS, GAUGES, PORTS, pairs


def run(args):
    args.reference='local'
    model=CoupledMouth(args)
    order=ROOT/f'mouth_free_order_cubic_n{args.mesh}_d04/order_seed1.55.npz'
    geo,chi,_=archive_inputs(order,None)
    with np.load(ROOT/f'mouth_flared_cubic_n{args.mesh}/field.npz') as data:
        field=data['psi']
    edges,keep,bulk=volume_cost(geo,chi,args.sub,5,1e-8)
    contact,load,B,G,_,tau,_=wall_data(geo,chi,field,edges,keep,5)
    S=bulk+contact
    decay=B@G/tau
    stress_dofs=len(S)
    n=model.nc
    scale=1/np.sqrt(S.diagonal())
    def by_bend(rate):
        metric=S+decay/rate
        return scale[:,None]*linalg.solve(scale[:,None]*metric*scale,scale[:,None]*load)
    def force_matrix(rate):
        z=model.carrier.force_matrix(rate).copy()
        z[1:,1:]+=ANGULAR*load.T@by_bend(rate)
        return z
    # Pick up the material force on the actual wall from waist to middle.
    middle=DEPTHS[1]
    lo=np.searchsorted(edges,middle)-1
    frac=(middle-edges[lo])/(edges[lo+1]-edges[lo])
    band=np.zeros(len(edges))
    band[lo:lo+2]=[1-frac,frac]
    band=-ANGULAR*geo.radius(np.array([middle]))[0][0]*band[keep]
    s=.30+.40j
    q=linalg.solve(force_matrix(s),PORTS)
    emit(test='local_stress_resolvent',**pairs('gauge',GAUGES@q),
         **pairs('bend_holding',ANGULAR*load.T@by_bend(s)),
         upstream_band_gain=abs(band@by_bend(s)@q[1:]).tolist())
    empty=sparse.csc_matrix((n-2,stress_dofs))
    drive=sparse.vstack((empty,load.T),format='csc')
    energy=sparse.bmat([[model.car_energy,ANGULAR*drive],[None,decay]],format='csc')
    time=sparse.bmat([[model.car_time,None],[drive.T,-S]],format='csc')
    ports=np.zeros((energy.shape[0],2))
    ports[model.q]=PORTS
    dt=args.dt
    advance=splu((energy/2-time/dt).tocsc())
    state=np.zeros_like(ports)
    work=np.zeros(2)
    renewal=np.zeros(2)
    scale_work=np.zeros(2)
    worst=0.
    for i in range(round(tau/dt)):
        mid=(i+.5)*dt
        f=np.sin(np.pi*mid/8)**2 if mid<8 else 0.
        delta=advance.solve(f*ports-energy@state)
        m=(state+delta/2)[n:]
        step_renewal=dt*ANGULAR*np.sum(m*(decay@m),axis=0)
        step_work=f*np.sum(PORTS*delta[model.q],axis=0)
        renewal+=step_renewal
        work+=step_work
        scale_work+=abs(step_work)+abs(step_renewal)
        state+=delta
        t=(i+1)*dt
        if t in (4.,8.,16.,round(tau/dt)*dt):
            stored=(np.sum(state[:n]*(model.car_energy@state[:n]),axis=0)
                    +ANGULAR*np.sum(state[n:]*(S@state[n:]),axis=0))/2
            error=abs(stored+renewal-work)/scale_work
            worst=max(worst,max(error))
            emit(test='local_stress_pulse',t=t,gauge_per_peak=(GAUGES@state[model.q]).tolist(),
                 upstream_band_force_per_peak=(band@state[n:]).tolist(),
                 stored_energy=stored.tolist(),relaxation_work=renewal.tolist(),
                 input_work=work.tolist(),relative_work_error=error.tolist())
    if not np.isfinite(worst) or worst>2e-5:
        raise ArithmeticError('Stress-history pulse work check failed')
    left,right=.035,.08
    def det(rate):
        return np.linalg.det(force_matrix(rate))
    bounds=[det(left),det(right)]
    root=brentq(det,left,right,xtol=1e-13)
    _,_,vh=np.linalg.svd(force_matrix(root))
    q=vh[-1]
    carrier=model.carrier
    rhs=(carrier.bend+root*carrier.motion)@q
    h=splu((carrier.operator-root*carrier.skew).tocsc()).solve(rhs)
    m=by_bend(root)@q[1:]
    ref=G@m/(tau*root)
    v=np.r_[h,q,ref]
    av,tv=model.energy@v,model.time@v
    err=np.linalg.norm(av-root*tv)/(np.linalg.norm(av)+root*np.linalg.norm(tv))
    gauges=GAUGES@q
    emit(test='local_positive_mode',bracket=[left,right],determinants=bounds,rate=root,
         original_relative_residual=err,gauges=(gauges/gauges[np.argmax(abs(gauges))]).tolist())
    if not np.isfinite(err) or err>1e-6:
        raise ArithmeticError('Positive mode failed original coupled equations')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--mesh',type=int,choices=[26,34],required=True)
    p.add_argument('--extension',choices=['compact','wide'],required=True)
    p.add_argument('--dt',type=float,choices=[.25,.125],required=True)
    p.add_argument('--sub',type=int,choices=[1,2,3,4],required=True)
    run(p.parse_args())
