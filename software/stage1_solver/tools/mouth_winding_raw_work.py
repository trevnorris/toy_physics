"""Original-field energy at signed exit radii on independently held vortex solutions.

Reuse the central stationary fields. Differentiate at fixed residual-amplitude
and meridional-phase coordinates. This is an envelope derivative for a stationary
discretization and a directional work diagnostic when changing the quadrature.
"""
import argparse
import hashlib
import itertools
import json
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss

from mouth_circulation_material import profile
from mouth_exchange_carrier import REGIONS, assemble, mapped_nodes
from mouth_exit_witness import inputs
from mouth_live_waist_fields import ORDER_INERTIA, WIDTH
from mouth_poly_shallow_release import ROOT, Trial, emit
from mouth_sealed_shape_relax import Meridian


def basis(frac):
    bits=np.array(list(itertools.product((0,1),repeat=len(frac))))
    terms=np.where(bits,np.array(frac),1-np.array(frac))
    shape=terms.prod(axis=1)
    grad=[(2*bits[:,k]-1)*np.delete(terms,k,axis=1).prod(axis=1) for k in range(len(frac))]
    return shape,grad


def quadrature(count):
    x,w=leggauss(count)
    return list(zip((x+1)/2,w/2))


def topology(shape,mesh,polar,outer,roof):
    _,_,grids,*_=assemble(0.,mesh,WIDTH/12,outer,roof,lambda _:shape,True,True)
    nodes=np.arange((1+max(g.max() for g in grids.values()))*(polar+1)).reshape(-1,polar+1)
    axis=np.unique(np.concatenate([grids[k,0][:,0] for k in range(3)]))
    nodes[axis]=nodes[axis,:1]
    _,flat=np.unique(nodes,return_inverse=True)
    phase={k:flat.reshape(nodes.shape)[g] for k,g in grids.items()}
    density=np.arange(flat.max()+1)
    density[phase[2,0][0].ravel()]=phase[1,0][-1].ravel()
    _,dids=np.unique(density,return_inverse=True)
    return phase,dids


def amplitudes_at(p,shape,vertices,args,labels,did,order,fields):
    radius,mesh,polar,outer,roof,degree,nq=args
    xs,zs=mapped_nodes(mesh,shape.length,outer,roof)
    angles=np.pi/2*np.linspace(0,1,polar+1)**1.35
    total=np.zeros(4)
    by_region={}
    points=quadrature(nq)
    for slab,ws in enumerate(zs):
        for side,xs_side in enumerate(xs):
            axes=(ws,xs_side,angles)
            k,i,l=np.indices(tuple(len(u)-1 for u in axes))
            x0,dx=xs_side[i].ravel(),np.diff(xs_side)[i].ravel()
            w0,dw=ws[k].ravel(),np.diff(ws)[k].ravel()
            a0,da=angles[l].ravel(),np.diff(angles)[l].ravel()
            grid=labels[slab,side]
            ids=np.stack([grid[k+bb,i+aa,l+cc] for aa,bb,cc in itertools.product((0,1),repeat=3)],axis=-1).reshape(-1,8)
            vv=vertices[did[ids]]; pp=p[ids]; center=vv.mean(axis=1)
            # Replicate the cell-central order gradient only for the force control
            # with matched quadrature. Raw fields and their derivatives stay pointwise.
            if slab==1 and side==0:
                x=x0+dx/2; w=w0+dw/2
                b,sl=shape.radius(w); t=w/shape.length
                xi=x*(order.shape[1]-1);tj=t*(order.shape[0]-1)
                ii=xi.astype(int);jj=tj.astype(int); a=xi-ii; h=tj-jj
                gr=(order.shape[1]-1)*((1-h)*(order[jj,ii+1]-order[jj,ii])+h*(order[jj+1,ii+1]-order[jj+1,ii]))/b
                gw=(order.shape[0]-1)*((1-a)*(order[jj+1,ii]-order[jj,ii])+a*(order[jj+1,ii+1]-order[jj,ii+1]))/shape.length-x*sl*gr
            else:
                gr=gw=0.
            gr2=gr**2+gw**2
            def coords(fr,fw,fa):
                x=x0+fr*dx; w=w0+fw*dw; angle=a0+fa*da
                b,slope=shape.radius(w); width=outer-b if side else b
                r=b+x*width if side else x*b
                tilt=slope*((1-x) if side else x)
                sin=np.sin(angle);cos=np.cos(angle)
                dj=4*np.pi*r*r*sin*width*dx*dw*da
                return width,r,tilt,sin,cos,dj
            if slab==1 and side==0:
                vol=np.zeros(len(ids));mean0=np.zeros(len(ids))
                for (a,wa),(h,wh),(c,wc) in itertools.product(quadrature(3),repeat=3):
                    _,r,_,sin,_,dj=coords(a,h,c)
                    vol+=dj*wa*wh*wc
                    mean0+=dj*wa*wh*wc*fields(r*sin)[0]**2
                mean0/=vol
            else:
                mean0=1.
            sums=np.zeros(4)
            for (a,wa),(h,wh),(c,wc) in itertools.product(points,repeat=3):
                width,r,tilt,sin,cos,dj=coords(a,h,c)
                shp,g=basis([a,h,c])
                v=vv@shp
                f0,fp=fields(r*sin)
                n0=f0*f0;n=n0*v*v
                vr=(vv@g[0])/(dx*width)
                vw=(vv@g[1])/dw-tilt*vr
                va=(vv@g[2])/(da*r)
                tr=(pp@g[0])/(dx*width)
                tw=(pp@g[1])/dw-tilt*tr
                ta=(pp@g[2])/(da*r)
                kinetic=tr*tr+tw*tw+ta*ta
                dot=gr*tr+gw*tw
                original=.5*n*(kinetic-ORDER_INERTIA*dot**2/(n+ORDER_INERTIA*gr2))
                lumped_phase=.5*n0*center**2*(kinetic-ORDER_INERTIA*dot**2/(mean0*center**2+ORDER_INERTIA*gr2))
                # The neutral vortex gradient points in cylindrical e_s, whose
                # spherical components are (sin(angle), cos(angle), 0 along w).
                grad=.5*((v*v-1)*fp*fp+2*f0*fp*v*(sin*vr+cos*va)+n0*(vr*vr+vw*vw+va*va))
                angular=.5*(n-n0)*degree**2/(r*sin)**2
                interaction=((n**5-n0**5)-5*(n-n0))/20
                # Exact correction for the native cell-averaged v in the
                # integrated-by-parts local potential, separated for comparison.
                xl=center**2; xf=v*v
                quadrature_change=n0**5*((xl**5-xf**5)-5*(xl-xf))/20
                weight=dj*wa*wh*wc
                sums+=np.array([weight@(grad+angular+interaction),weight@original,
                                weight@quadrature_change,weight@(lumped_phase-original)])
            total+=sums
            by_region[REGIONS[2*slab+side]]=sums.tolist()
    # Terminal layer joins common N*phi at fixed w, so only meridional phase
    # appears in the jump. The face has a physical three-volume.
    x0,ang=np.meshgrid(xs[0][:-1],angles[:-1],indexing='ij')
    dx,da=np.meshgrid(np.diff(xs[0]),np.diff(angles),indexing='ij')
    ids=labels[1,0][-1]; ids2=labels[2,0][0]
    i,k=np.indices((mesh,polar))
    lower=np.stack([ids[i+aa,k+cc] for aa,cc in itertools.product((0,1),repeat=2)],axis=-1).reshape(-1,4)
    upper=np.stack([ids2[i+aa,k+cc] for aa,cc in itertools.product((0,1),repeat=2)],axis=-1).reshape(-1,4)
    v=vertices[did[lower]]
    phase=p[upper]-p[lower]
    cap=np.zeros(2)
    for (a,wa),(c,wc) in itertools.product(points,repeat=2):
        r=shape.exit*(x0+a*dx)
        sin=np.sin(ang+c*da)
        S=basis([a,c])[0]
        measure=(4*np.pi*r*r*sin*shape.exit*dx*da*wa*wc).ravel()
        fv=fields((r*sin).ravel())[0]
        energy=.5*fv**4*(phase@S)**2/(WIDTH/12)
        cap[0]+=measure@((v@S)**4*energy)
        cap[1]+=measure@(v.mean(axis=1)**4*energy)
    return dict(continuous_raw=total[0]+total[1]+cap[0],
                original_with_native_lumping=total.sum()+cap[1],regions=by_region,cap=cap.tolist())


def check(name, args):
    folder=ROOT/name
    guard=json.loads((folder/'outcome.json').read_text())
    lines=list(map(json.loads,(folder/'stdout').read_text().splitlines()))
    assert guard['limitsVerified'] and guard['exitCode']==guard['childOutcome']['exitCode']==0
    conf=lines[0]
    with np.load(folder/'fields.npz',allow_pickle=False) as archive:
        v=archive['amplitude_ratio']; phase=archive['phase']
        assert archive['material_sha256'].item()==conf['material_sha256']
    radius,degree=conf['radius'],conf['winding']
    assert radius in (1,2)
    trial=Trial(radius,72,bonded_core=True)
    path=ROOT/f'_scratch/s11c/mouth_exit_release01/r{radius}n72/endpoint.npz'
    shifts,chi,_,digest,status,*_=inputs(path,radius,72,trial)
    assert digest==conf['material_sha256'] and status=='balanced'
    cfg=trial.x+shifts[:4]; center=Meridian(radius,cfg)
    b=center.exit
    n,nt,outer,roof=(conf[k] for k in ['mesh','polar_cells','outer','roof'])
    phase_ids,dids=topology(center,n,nt,outer,roof)
    fields=profile(degree,max(2*outer,80.),2e-9)[0] if degree else lambda s:(np.ones_like(s),np.zeros_like(s))
    emit(kind='input',radius=radius,winding=degree,material_sha256=digest,mesh=n,polar=nt,
         center_field_sha256=hashlib.sha256((folder/'fields.npz').read_bytes()).hexdigest())
    for nq in args.quad:
        for step in args.steps:
            work=[]
            for sign in (-1,1):
                shape_cfg=cfg.copy();shape_cfg[1]+=np.log1p(sign*step)
                shape=Meridian(radius,shape_cfg)
                work.append(amplitudes_at(phase,shape,v,(radius,n,nt,outer,roof,degree,nq),
                                          phase_ids,dids,chi,fields))
            div=2*step*b
            emit(kind='fixed_central_field_work',shape_step=step,quadrature=nq,
                 continuous_raw_force=(work[0]['continuous_raw']-work[1]['continuous_raw'])/div,
                 matched_lumping_raw_force=(work[0]['original_with_native_lumping']-work[1]['original_with_native_lumping'])/div,
                 signed_raw_energies=[w['continuous_raw'] for w in work],region_terms=[w['regions'] for w in work])


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('directories',type=Path,nargs='+')
    parser.add_argument('--quad',type=int,nargs='+',default=[3,5])
    parser.add_argument('--steps',type=float,nargs='+',default=[.003])
    args=parser.parse_args()
    for directory in args.directories:
        check(directory,args)
