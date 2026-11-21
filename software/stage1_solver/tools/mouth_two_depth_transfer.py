"""Spatially integrated attachment traction, pulse work and a free-mode check.

Reassemble the contact law before eliminating its dual stress. This also gives
the local stress trace and complementary energy needed for the readouts.
"""
import argparse
from pathlib import Path

import numpy as np
from scipy import linalg
from scipy.sparse.linalg import LinearOperator, eigs, splu

from mouth_bend_renewal import ROOT, volume_cost, wall_data
from mouth_bend_static_shear import ANGULAR, archive_inputs
from mouth_order_wave_pulse import emit
from mouth_two_depth_renewal import CoupledMouth, DEPTHS, GAUGES, PORTS, pairs


def free_modes(model):
    shift=.1
    fac=splu((model.energy-shift*model.time).tocsc())
    inv=LinearOperator(model.energy.shape,matvec=lambda v: fac.solve(model.time@v),dtype=float)
    seed=np.random.default_rng(31927).normal(size=model.energy.shape[0])
    vals,vectors=eigs(inv,k=4,which='LM',v0=seed,tol=1e-10)
    for value,v in zip(vals,vectors.T,strict=True):
        rate=shift+1/value
        av,tv=model.energy@v,model.time@v
        error=np.linalg.norm(av-rate*tv)/(np.linalg.norm(av)+abs(rate)*np.linalg.norm(tv))
        motion=GAUGES@v[model.q]
        motion=motion/motion[np.argmax(abs(motion))]
        emit(test='free_mode',rate=[float(rate.real),float(rate.imag)],
             relative_residual=float(error),**pairs('motion',motion))
        if not np.isfinite(error) or error>3e-5:
            raise ArithmeticError('Unforced field and attachment mode failed')


def run(args):
    model=CoupledMouth(args)
    order=ROOT/f'mouth_free_order_cubic_n{args.mesh}_d04/order_seed1.55.npz'
    geo,chi,_=archive_inputs(order,None)
    with np.load(ROOT/f'mouth_flared_cubic_n{args.mesh}/field.npz') as data:
        field=data['psi']
    edges,keep,bulk=volume_cost(geo,chi,args.sub,5,1e-8)
    contact,load,B,G,drift,tau,_=wall_data(geo,chi,field,edges,keep,5)
    if args.reference=='local':
        drift=np.zeros_like(drift)
    metric=bulk+contact
    # Material-on-wall resultant between two depth cuts is -ANGULAR [a*m].
    # The endpoint tractions vanish. The middle gauge divides the whole guide.
    middle=DEPTHS[1]
    i=np.searchsorted(edges,middle)-1
    fraction=(middle-edges[i])/(edges[i+1]-edges[i])
    cut=np.zeros(len(edges))
    cut[i:i+2]=[1-fraction,fraction]
    band=-ANGULAR*geo.radius(np.array([middle]))[0][0]*cut[keep]
    s=.30+.40j
    memory=linalg.solve(s*np.eye(len(drift))+drift,G/tau)
    cost=metric+B@memory
    scale=1/np.sqrt(np.diag(metric))
    mper=scale[:,None]*linalg.solve(scale[:,None]*cost*scale,scale[:,None]*load)
    holding=ANGULAR*load.T@mper
    full=model.carrier.force_matrix(s).copy()
    full[1:,1:]+=holding
    q=linalg.solve(full,PORTS)
    mode=mper@q[1:]
    emit(test='material_resolvent',**pairs('bend_holding',holding),
         **pairs('gauge',GAUGES@q),**pairs('upstream_band_force',band@mode),
         upstream_band_gain=abs(band@mode).tolist())

    dt=args.dt
    # Trapezoidal midpoint update in next-state form. Carrier energy excludes
    # renewing reference coordinates. Their transport and relaxation have a
    # separate signed work account.
    advance=splu((model.energy/2-model.time/dt).tocsc())
    past=-model.energy/2-model.time/dt
    state=np.zeros_like(model.forces)
    drive_work=np.zeros(2)
    abs_input=np.zeros(2)
    relief=np.zeros(2)
    transported=np.zeros(2)
    scale_work=np.zeros(2)
    worst=0.
    for j in range(round(tau/dt)):
        t=(j+.5)*dt
        strength=np.sin(np.pi*t/8)**2 if t<8 else 0.
        future=advance.solve(past@state+strength*model.forces)
        mid=(future+state)/2
        qp=mid[model.q]
        ref=mid[model.nc:]
        dual=model.dual_bend@qp[1:]-model.dual_ref@ref
        conjugate=ANGULAR*(B.T@dual)
        renewal_step=dt*np.sum(conjugate*(G@dual/tau),axis=0)
        transport_step=-dt*np.sum(conjugate*(drift@ref),axis=0)
        relief+=renewal_step
        transported+=transport_step
        applied=strength*np.sum(PORTS*(future[model.q]-state[model.q]),axis=0)
        drive_work+=applied
        abs_input+=abs(applied)
        scale_work+=abs(applied)+abs(renewal_step)+abs(transport_step)
        state=future
        end=(j+1)*dt
        if end not in (4.,8.,16.,round(tau/dt)*dt):
            continue
        bare=state[:model.nc]
        m=model.dual_bend@state[model.q][1:]-model.dual_ref@state[model.nc:]
        kinetic=.5*np.sum(bare*(model.car_energy@bare),axis=0)
        elastic=ANGULAR*.5*np.sum(m*(metric@m),axis=0)
        imbalance=kinetic+elastic+relief+transported-drive_work
        error=abs(imbalance)/scale_work
        worst=max(worst,max(error))
        emit(test='pulse_budget',t=end,upstream_band_force_per_peak=(band@m).tolist(),
             gauge_per_peak=(GAUGES@state[model.q]).tolist(),
             carrier_energy=kinetic.tolist(),material_energy=elastic.tolist(),
             relaxation_work=relief.tolist(),reference_transport_work=transported.tolist(),
             input_work=drive_work.tolist(),relative_work_error=error.tolist(),
             abs_input_work=abs_input.tolist())
    if not np.isfinite(worst) or worst>2e-5:
        raise ArithmeticError('Carrier, material and reference energy accounts do not close')
    free_modes(model)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--mesh',type=int,choices=[26,34],required=True)
    p.add_argument('--extension',choices=['compact','wide'],required=True)
    p.add_argument('--dt',type=float,choices=[.25,.125],required=True)
    p.add_argument('--sub',type=int,choices=[1,2,3,4],required=True)
    p.add_argument('--reference',choices=['travelling','local'],required=True)
    run(p.parse_args())
