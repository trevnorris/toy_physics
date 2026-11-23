"""Wave-induced order work and the response of the calibrated order coordinate on a held exit."""
import argparse
import hashlib
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
from scipy import sparse
from scipy.optimize import brentq
from scipy.special import expit,logit

from mouth_exchange_carrier import BIAS
from mouth_exit_release import Release
from mouth_gallery_action import archived
from mouth_live_waist_order import WELL
from mouth_live_waist_wave import transfer
from mouth_poly_shallow_release import emit
from mouth_trapping_modes import Sector,pole,form
from mouth_strained_wave import cell_tangents
from mouth_trapping_angular_tail import recorded_seed


def guarded_action(radius,n,digest):
    path=Path('_scratch/s11c/mouth_gallery_action01')/f'r{radius}n{n}'
    meta=json.loads((path/'outcome.json').read_bytes())
    raw=(path/'stdout').read_bytes()
    rows=list(map(json.loads,raw.splitlines()))
    if (meta['exitCode'] or meta['childOutcome']['exitCode'] or not meta['limitsVerified']
            or rows[0]['digest']!=digest or rows[-1]['kind']!='exit_balance'):
        raise ValueError('Invalid matched action measurement')
    return rows[-1]['action'], hashlib.sha256(raw).hexdigest()


def run(radius,n):
    trial,sh,chi,old,digest,target=archived(radius,n)
    action,source=guarded_action(radius,n,digest)
    problem=Release(radius,n)
    _,_,owner,loads,initial,elastic=trial.evaluate(sh,chi,old)
    active=problem.active
    logits=logit(chi[active])
    frozen_reservoir=initial['energy']['static_outside']+initial['energy']['outer_order']
    op=Sector(trial,sh,chi,12)
    seed,_=recorded_seed(op,target)
    cells_per_depth=len(op.ids)//op.mesh
    projection=sparse.kron(transfer(n,np.linspace(0,1,op.mesh+1)),
                           transfer(2*n,np.linspace(0,1,cells_per_depth+1)),format='csr')
    holder=SimpleNamespace(mesh=op.mesh,ids=op.ids,kc=op.mc,mc=op.mc)
    kc,mc=cell_tangents(holder,op.shape,op.coeff,degree=12,core_reference=op.reference)
    def qfield(q):
        c=chi.copy()
        c[active]=expit(logits+q*problem.rate)
        return c
    def state(c,wseed,J):
        e,g=owner.energy_gradient(c)
        ec,gc,_=loads.carrier.response(c,BIAS)
        value=e+ec+frozen_reservoir+(elastic*c).sum()
        g+=gc+elastic
        if J==0:
            return value,g,wseed,None
        obj=Sector(trial,sh,c,12)
        if not np.array_equal(op.live,obj.live) or not np.array_equal(op.labels,obj.labels):
            raise ArithmeticError('Shear support changed along order variation')
        choice,info=pole(obj,wseed)
        if info['seed_overlap']<.98:
            raise ArithmeticError('Lost the prepared occupied pole')
        omega,v=choice
        full=np.zeros(len(obj.live),complex); full[obj.live]=v
        vcell=full[obj.ids]
        dcell=np.einsum('ci,cij,cj->c',vcell,kc-omega**2*mc,vcell)
        dport=(obj.port.block(omega+2e-6)-obj.port.block(omega-2e-6))/4e-6
        edge=v[obj.port.ids]
        denom=2*omega*(v@(obj.M@v))+edge@(dport@edge)
        dwdc=(projection.T@(dcell/denom).real).reshape(n+1,2*n+1).T
        return value+J*omega.real,g+J*dwdc,choice,dwdc
    def state_q(q,wseed,J):
        c=qfield(q)
        value,g,choice,gw=state(c,wseed,J)
        dg=g[active]@(c[active]*(1-c[active])*problem.rate)
        return value,dg,choice,gw,c

    v0,g0,_,_=state(chi,seed,0.)
    baseline=state_q(0.,seed,action)
    value,slope,_,gw,_=baseline
    if abs(value-v0-action*target['frequency'])>1e-7:
        raise ArithmeticError('Unmatched wave energy')
    for h in (.0002,.0001):
        low=state_q(-h,seed,action)[0]
        high=state_q(h,seed,action)[0]
        if abs((high-low)/(2*h)-slope)>3e-5*max(1,abs(slope)):
            raise ArithmeticError('Unresolved order work on the represented coordinate')
        emit(kind='coordinate_work',h=h,analytic=slope,direct=(high-low)/(2*h))
    # Signed work at an independently movable site with room for both perturbations.
    # The signed step is a probe, not an inferred equilibrium change or failure strain.
    smooth=owner.free & (chi>.02) & (chi<.98)
    volume=owner.lump.reshape(chi.shape)/WELL
    i=np.unravel_index(np.argmax(np.where(smooth,abs(gw)/volume,-1)),chi.shape)
    for h in (.001,.0005):
        d=np.zeros_like(chi); d[i]=h
        lo=state(chi-d,seed,action); hi=state(chi+d,seed,action)
        dw=(hi[2][0]-lo[2][0]).real/(2*h)
        if abs(dw-gw[i])>2e-4*max(abs(gw[i]),1e-6):
            raise ArithmeticError('Unresolved local incremental order work')
        emit(kind='local_order_work',radius=radius,n=n,h=h,index=list(map(int,i)),local_chi=chi[i],
             wave_slope=action*gw[i],direct_wave_slope=action*dw,material_slope=g0[i],
             changes=[lo[0]-value,hi[0]-value],
             wave_work_per_volume=action*gw[i]/volume[i],max_homogeneous_well_slope=WELL/(3*np.sqrt(3)))
    emit(kind='origin',radius=radius,n=n,order_field_sha=digest,threshold_sha=source,action=action,
         old_wave=old[0].real,degree12=target['frequency'],peak_first_order_chi=0.)
    # Quasistatic response along exactly the archived order coordinate, at fixed geometry.
    # Roots bracketed by descent, positive curvature checked at both step sizes.
    final={}
    for label,J in [('no_wave',0.),('threshold',action)]:
        q=0.
        value,slope,branch,_,c=state_q(q,seed,J)
        start=value
        sign=-np.sign(slope)
        if sign==0:
            raise ArithmeticError('Exactly zero initial slope requires curvature test')
        dq=sign*.005
        while True:
            proposed=q+dq
            trial_c=qfield(proposed)
            if np.max(abs(trial_c-c))>.025:
                dq/=2
                continue
            if np.max(abs(trial_c-c))<1e-12:
                emit(kind='order_unresolved',label=label,reason='No resolvable change in order',q=q)
                break
            try:
                val,dg,next_seed,_,next_c=state_q(proposed,branch,J)
            except (ValueError,ArithmeticError,RuntimeError) as err:
                dq/=2
                if np.max(abs(qfield(q+dq)-c))<1e-8:
                    emit(kind='order_unresolved',label=label,reason=str(err),q=q);break
                continue
            if slope*dg<=0:
                bracket=sorted((proposed,q))
                def projected_grad(x):return state_q(x,branch,J)[1]
                root=brentq(projected_grad,*bracket,xtol=2e-10)
                ev,dlast,selected,_,res=state_q(root,branch,J)
                curv=[(projected_grad(root+h)-projected_grad(root-h))/(2*h) for h in (.0002,.0001)]
                if min(curv)<=0 or abs(dlast)>1e-6 or ev>start+1e-8:
                    raise ArithmeticError('Projected minimum failed tests')
                final[label]=res
                di=np.unravel_index(np.argmax(abs(res-chi)),chi.shape)
                emit(kind='order_minimum',label=label,q=root,peak_change=float(abs(res-chi).max()),
                     signed_change_here=res[di]-chi[di],index=list(map(int,di)),chi_before=chi[di],
                     chi_after=res[di],energy_change=ev-start,curvature=curv,
                     stored_action=J,frequency=selected[0].real if J else None)
                break
            q,value,slope,branch,c=proposed,val,dg,next_seed,next_c
            dq*=1.4
    if len(final)==2:
        diff=final['threshold']-final['no_wave']
        i=np.unravel_index(np.argmax(abs(diff)),chi.shape)
        emit(kind='paired_response',peak=float(abs(diff).max()),signed=diff[i],index=list(map(int,i)),
             chi_off=final['no_wave'][i],chi_on=final['threshold'][i],reference=chi[i])


if __name__=='__main__':
    cli=argparse.ArgumentParser(description=__doc__)
    cli.add_argument('--radius',required=True,type=int,choices=(1,2))
    cli.add_argument('--mesh',required=True,type=int,choices=(48,72))
    a=cli.parse_args();run(a.radius,a.mesh)
