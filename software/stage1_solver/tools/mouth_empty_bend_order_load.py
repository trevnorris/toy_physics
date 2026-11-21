"""Measure the passive conversion field's separate bare-energy traction on the tilted sidewall.

The carrier drives order, with no reciprocal update to carrier momentum. The
integrating factor is used only to solve the stationary equation. Physical
wall traction is evaluated from the unweighted order energy and checked against
volume stress plus the force on the advected material.
"""
import argparse
import json
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.sparse.linalg import splu

from mouth_curved_phase_stability import weighted_gram
from mouth_empty_bend_stiffness import emit,prepare,quadrature,shapes,trace
from mouth_free_order_diffuse import NU,SIGMA
from mouth_order_wave_pulse import digest


def main(a):
    dom,psi=prepare(a)
    archive=json.loads(Path('research/pde_ledger_v3/_measurements/mouth_free_order_depth.json').read_text())
    if digest(a.order)!=archive['diffuse_orders_sha256'][str(a.order)]:
        raise ValueError('Order does not match the archived preparation')
    for p in (a.tangent,a.order.parent):
        limits=json.loads((p/'limit-validation.json').read_text())
        o=json.loads((p/'outcome.json').read_text())
        if not limits['verified'] or o['exitCode'] or o['childOutcome']['exitCode']:
            raise ValueError('Invalid carrier tangent or stationary order')
    order_meta=json.loads((a.order.parent/'stdout').read_text().splitlines()[0])
    tangent_meta=[json.loads(line) for line in (a.tangent/'stdout').read_text().splitlines()]
    source_sha=digest(a.source/'field.npz')
    if (order_meta['source_sha256']!=source_sha or tangent_meta[0]['sha256']!=source_sha
            or tangent_meta[-1]['extension']!='compact' or order_meta['nu']!=NU or order_meta['sigma']!=SIGMA):
        raise ValueError('Order and carrier tangent must share the material, carrier field and compact chart')
    with np.load(a.tangent/'tangent.npz') as data:
        np.testing.assert_array_equal(dom.points,data['points'])
        np.testing.assert_array_equal(psi,data['field'])
        v=data['tangents']
    with np.load(a.order) as z: c=z['chi']
    r,w,gr,gw=quadrature(dom)
    n=dom.N; q=dom.weight
    h,sr,sw=shapes(r,w,'compact',dom.far,dom.roof)
    f=n@psi; fr=gr@(psi-1); fw=gw@(psi-1)
    vq=n@v
    phase=np.angle(f); phase1=(vq/f[:,None]).imag
    pr=(fr/f).imag; pw=(fw/f).imag
    p1r=((gr@v-vq*(fr/f)[:,None])/f[:,None]).imag
    p1w=((gw@v-vq*(fw/f)[:,None])/f[:,None]).imag
    d=order_meta['delta']
    kappa=6*SIGMA*d
    mobility=NU/kappa
    c0=n@c; rem=n@(1-c); cr=gr@c; cw=gw@c
    well=c0*rem*(rem-c0)/d**2
    weight=q*np.exp(-phase/NU)
    coord=sr*cr[:,None]+sw*cw[:,None]
    def source(weight,relative_phase):
        return (gr.T@(weight[:,None]*(coord+relative_phase*cr[:,None]))
                +gw.T@(weight[:,None]*(sw*cr[:,None]-sr*cw[:,None]+relative_phase*cw[:,None]))
                +n.T@(weight[:,None]*(2*coord/r[:,None]+(relative_phase-sr)*well[:,None])))
    rhs=source(weight,phase1/NU)
    regular=dom.free&(dom.points[:,0]>0)
    jac=(weighted_gram(gr,weight)+weighted_gram(gw,weight)
         +weighted_gram(n, weight*((1-6*c0*rem)/d**2+2/r**2)))
    y=np.zeros_like(v.real)
    y[regular]=splu(jac[regular][:,regular].tocsc()).solve(rhs[regular])
    error=np.max(abs((jac@y-rhs)[regular]))
    # Direct physical sidewall integral with natural n·grad chi=0.
    pg=dom.blocks[2][0]; x,weights=leggauss(10)
    hp=np.diff(pg)[:,None]/2; p=(pg[:-1,None]+hp*(1+x)).ravel()
    rr,ww,_,_,_,_,val,rrgrad,wwgrad=trace(dom,2,p,np.zeros_like(p),np.column_stack([c,y]))
    hh,srr,sww=shapes(rr,ww,'compact',dom.far,dom.roof)
    surface=rrgrad[:,0,None]*(rrgrad[:,1:]-srr*rrgrad[:,0,None])+wwgrad[:,0,None]*(wwgrad[:,1:]-sww*rrgrad[:,0,None])
    surface+=(val[:,0]*(1-val[:,0])*(1-2*val[:,0])/d**2)[:,None]*val[:,1:]
    wall=kappa/3*hh.T@((4*np.pi*rr**2*(hp*weights).ravel())[:,None]*surface)
    # Independent shape derivative over volume plus integral of s_i * mu_chi * d_x chi.
    dd=q*cr**2
    vol=kappa/3*(sr.T@(dd[:,None]*sr)+sw.T@(dd[:,None]*sw)-source(q,0.).T@y)
    dq=n@y; dr=gr@y; dw=gw@y
    mu=-(pr*cr+pw*cw)/mobility
    m1=-((p1r-sr*pr[:,None])*cr[:,None]+(p1w-sw*pr[:,None])*cw[:,None]+
         (dr-sr*cr[:,None])*pr[:,None]+(dw-sw*cr[:,None])*pw[:,None])/mobility
    body=h.T@((q/3)[:,None]*(m1*cr[:,None]+mu[:,None]*(dr+2*dq/r[:,None])))
    emit(test='order_traction',mesh=len(dom.blocks[2][1])-1,order_sha256=digest(a.order),
         tangent_sha256=digest(a.tangent/'tangent.npz'),trial_interface_tension=SIGMA,
         kappa=kappa,mobility=mobility,width=d,order_file=str(a.order),
         direct_wall=wall.tolist(),bare_volume_virtual_work=vol.tolist(),conversion_body_term=body.tolist(),
         reconstructed_wall=(vol+body).tolist(),max_wall_disagreement=np.max(abs(wall-vol-body)),
         dipole_residual=error,maximum_order_tangent=np.max(abs(y),axis=0).tolist())
    if error>1e-8 or not np.isfinite(wall).all(): raise ArithmeticError('Order tangent failed')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source',type=Path,required=True)
    p.add_argument('--tangent',type=Path,required=True)
    p.add_argument('--order',type=Path,required=True)
    main(p.parse_args())
