"""Reordering work at previously identified radius-one sites, on both held mouths."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss

from mouth_circulation_material import material_parts, profile
from mouth_live_waist_fields import ORDER_INERTIA
from mouth_poly_shallow_release import emit
from mouth_wound_carrier import order_slopes


def contained_edges(material,coarse,low,high):
    return np.unique(np.r_[low,high,material[(material>low)&(material<high)],
                           coarse[(coarse>low)&(coarse<high)]])


def q1_and_slopes(values,x,t,shape,b,slope):
    nt,nr=np.subtract(values.shape,1)
    i=np.minimum(np.floor(x*nr).astype(int),nr-1)
    j=np.minimum(np.floor(t*nt).astype(int),nt-1)
    a,beta=x*nr-i,t*nt-j
    grid=values
    lo=(1-a)*grid[j,i]+a*grid[j,i+1]
    hi=(1-a)*grid[j+1,i]+a*grid[j+1,i+1]
    value=(1-beta)*lo+beta*hi
    gr=nr*((1-beta)*(grid[j,i+1]-grid[j,i])+beta*(grid[j+1,i+1]-grid[j+1,i]))/b
    gw=nt*(hi-lo)/shape.length-x*slope*gr
    return value,gr,gw


def reconstructed(z,shape,chi,vectors,level):
    # Resolve the footprint of a material test, even when its two opposite
    # slopes are narrower than a cell of the meridional carrier. Coarse
    # phase and amplitude are interpolated, with no claim of new fluid resolution.
    xs,ts,angs=(z[key] for key in ('core_x_nodes','core_height_nodes','core_angle_nodes'))
    amp=z['core_amp']
    phase=z['core_phase']
    m=int(z['winding'])
    vortex=(profile(m,64.,1.e-9)[0] if m else lambda s:(np.ones_like(s),np.zeros_like(s)))
    q,wt=leggauss(level)
    probes=[]
    rq=np.linspace(0,1,chi.shape[1]);tq=np.linspace(0,1,chi.shape[0])
    for key,direction in vectors.items():
        js,ix=np.nonzero(direction)
        ex=contained_edges(rq,xs,rq[max(ix.min()-1,0)],rq[min(ix.max()+1,len(rq)-1)])
        ez=contained_edges(tq,ts,tq[max(js.min()-1,0)],tq[min(js.max()+1,len(tq)-1)])
        def quadrature(edges):
            return (edges[:-1,None]+np.diff(edges)[:,None]*(1+q)/2).ravel(),(np.diff(edges)[:,None]*wt/2).ravel()
        x,xw=quadrature(ex);t,tw=quadrature(ez)
        t,x=np.meshgrid(t,x,indexing='ij')
        drdt=tw[:,None]*xw[None,:]
        i=np.searchsorted(xs,x)-1
        j=np.searchsorted(ts,t)-1
        dx,dt=xs[i+1]-xs[i],ts[j+1]-ts[j]
        xx,tt=(x-xs[i])/dx,(t-ts[j])/dt
        b,slope=shape.radius(t*shape.length)
        r=x*b
        gr,gw=order_slopes(chi,x,t,b,slope,shape.length)
        _,ar,aw=q1_and_slopes(direction,x,t,shape,b,slope)
        density_weight=4*np.pi*r*r*b*shape.length*drdt
        result,finite,abswork=0.,0.,0.
        # All prescribed order tests are axisymmetric. The polar component of
        # the order gradient is zero, including after the variation.
        for k,(start,stop) in enumerate(zip(angs[:-1],angs[1:])):
            for c,wc in zip((1+q)/2,wt/2):
                angle=start+(stop-start)*c
                vol=density_weight*(np.sin(angle)*(stop-start)*wc)
                def slice_values(grid):
                    return (1-c)*grid[...,k]+c*grid[...,k+1]
                am=slice_values(amp)
                a0=(1-xx)*am[j,i]+xx*am[j,i+1]
                a1=(1-xx)*am[j+1,i]+xx*am[j+1,i+1]
                n=(vortex(r*np.sin(angle))[0]*((1-tt)*a0+tt*a1))**2
                p=slice_values(phase)
                p0=(1-xx)*p[j,i]+xx*p[j,i+1]
                p1=(1-xx)*p[j+1,i]+xx*p[j+1,i+1]
                pr=((1-tt)*(p[j,i+1]-p[j,i])+tt*(p[j+1,i+1]-p[j+1,i]))/(dx*b)
                pw=(p1-p0)/(dt*shape.length)-x*slope*pr
                pdot=gr*pr+gw*pw
                adot=ar*pr+aw*pw
                gdot=gr*ar+gw*aw
                den=n+ORDER_INERTIA*(gr*gr+gw*gw)
                work=-vol*n*ORDER_INERTIA*pdot/den*(adot-ORDER_INERTIA*pdot*gdot/den)
                result+=work.sum()
                abswork+=abs(work).sum()
                def e(delta):
                    dot=pdot+delta*adot
                    denom=n+ORDER_INERTIA*((gr+delta*ar)**2+(gw+delta*aw)**2)
                    return -.5*vol*n*ORDER_INERTIA*dot**2/denom
                h=1e-6
                finite+=np.sum((e(h)-e(-h))/(2*h))
        probes.append(dict(probe=key,carrier_work=float(result),finite_difference_check=float(finite),
                           integrated_absolute_work=float(abswork),quadrature_order=level))
    return probes


def run(path):
    guard=json.loads((path.parent/'outcome.json').read_bytes())
    if not guard['limitsVerified'] or guard['exitCode'] or guard['childOutcome']['exitCode']:
        raise ValueError('Unverified source solve')
    raw=(path.parent/'stdout').read_bytes()
    rows=[json.loads(l) for l in raw.splitlines()]
    inp,=[r for r in rows if r['kind']=='input']
    if not any(r['kind']=='geometry' and r['relative_exit']==0 for r in rows):
        raise ValueError('No stationary central result')
    with np.load(path) as z:
        radius=int(z['radius'])
        if str(z['material_sha256'])!=inp['material_sha256'] or int(z['winding'])!=inp['winding'] or radius!=inp['radius']:
            raise ValueError('Mismatch with source solve')
        shape,chi,owner,probes,elastic,order=material_parts(radius)
        mesh=inp['mesh']
        emit(kind='input',radius=radius,winding=inp['winding'],mesh=mesh,polar=inp['polar_cells'],
             material_sha256=inp['material_sha256'],field_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
             source_stdout_sha256=hashlib.sha256(raw).hexdigest())
        x,t=z['core_x'],z['core_t'];b,slope=shape.radius(t*shape.length)
        for key,direction in probes.items():
            gr,gw=order_slopes(direction,x,t,b,slope,shape.length)
            carrier=gr@z['core_order_gr_work']+gw@z['core_order_gw_work']
            static=float(np.sum((elastic+order)*direction))
            sample=np.count_nonzero((gr!=0)|(gw!=0))
            emit(kind='discrete_derivative',probe=key,material_slope=static,carrier_slope=float(carrier),
                 sampled_core_cells=int(sample),sum=static+carrier)
        if 'core_amp' in z:
            for level in (2,3):
                for result in reconstructed(z,shape,chi,probes,level):
                    signed=result['carrier_work']
                    check=result['finite_difference_check']
                    static=float(np.sum((elastic+order)*probes[result['probe']]))
                    emit(kind='field_reconstruction',material_slope=static,sum=static+signed,**result)
                    if abs(signed-check)>2e-8*(1+abs(signed)):
                        raise ArithmeticError('Derivative disagrees with native Hamiltonian variation')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('central_fields',type=Path,nargs='+')
    for path in p.parse_args().central_fields:
        run(path)
