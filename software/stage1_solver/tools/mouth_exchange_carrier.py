"""Weak reversible carrier flow through the same perturbed meridian, with a surrounding bulk.

The thin Route A layer is a conservative sheet storing phase energy; sides
and cap transmit normal mass current. The only sealed internal boundary is
the intact feeding face outside the aperture. All integrals include 4*pi*r².
"""
import argparse
import json

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy import sparse
from scipy.sparse.linalg import spsolve

from mouth_exchange_geometry import ExchangeShape


BIAS = .4510084282857874
REGIONS = ('feed_inner','feed_outer','core','bulk_beside','bulk_above_inner','bulk_above_outer')


def mapped_nodes(n, length, outer, roof):
    u = np.linspace(0., 1., n+1)
    radial = (1-(1-u)**2, (np.sinh(2*u)/np.sinh(2))**2)
    axial = (-(1-u)**2, length*(1-np.cos(np.pi*u))/2,
             length+(roof-length)*(np.sinh(2*u)/np.sinh(2))**2)
    return radial, axial


def assemble(step, n, sheet, outer, roof, geometry=ExchangeShape, return_system=False):
    shape = geometry(step)
    delta = 2.e-5
    up, down = geometry(step+delta), geometry(step-delta)
    radial, axial = mapped_nodes(n,shape.length,outer,roof)
    q, weights = leggauss(3)
    roots, weights = (q+1)/2, weights/2
    pieces, nodes, indices = [], {}, {}
    volumes = {}
    for slab in range(3):
        for side in range(2):
            grid = np.empty((n+1,n+1),int)
            for z,r in np.ndindex(grid.shape):
                iz,ir = slab*n+z,side*n+r
                ordered = slab==0 or (slab==1 and side==0)
                bank = int(ordered) if sheet else int(iz==n and ir>n and slab==0)
                key = (iz,ir,bank)
                if key not in nodes:
                    nodes[key]=len(nodes)
                grid[z,r]=nodes[key]
            indices[slab,side]=grid

            dy,dw = np.diff(radial[side])[None,:],np.diff(axial[slab])[:,None]
            corners = np.stack((grid[:-1,:-1],grid[:-1,1:],grid[1:,:-1],grid[1:,1:]),axis=-1).reshape(-1,4)
            K=np.zeros((n,n,4,4))
            dK=np.zeros_like(K)
            cell_volume=np.zeros((n,n))
            cell_volume_derivative=np.zeros_like(cell_volume)
            for b,wb in zip(roots,weights):
                w=axial[slab][:-1,None]+b*dw
                radius,slope=shape.radius(w)
                move=(up.radius(w)[0]-down.radius(w)[0])/(2*delta)
                move_slope=(up.radius(w)[1]-down.radius(w)[1])/(2*delta)
                width=outer-radius if side else radius
                for a,wa in zip(roots,weights):
                    y=radial[side][None,:-1]+dy*a
                    alpha=1-y if side else y
                    r=radius+y*width if side else y*radius
                    jac=4*np.pi*r*r*width*dy*dw*(wa*wb)
                    basis_y=np.array([b-1,1-b,-b,b])/(dy[...,None]*width[...,None])
                    basis_w=np.array([a-1,-a,1-a,a])/dw[...,None]-(alpha*slope)[...,None]*basis_y
                    def gram(left,right):
                        return left[...,None]*right[...,None,:]
                    Grr=gram(basis_y,basis_y)
                    Gww=gram(basis_w,basis_w)
                    div_r=(1-2*side)*move/width
                    dflow_dw=alpha*(move_slope-div_r*slope)
                    div=div_r+2*alpha*move/r
                    cell_volume+=jac
                    cell_volume_derivative+=jac*div
                    K+=jac[...,None,None]*(Grr+Gww)
                    dK+=jac[...,None,None]*((div-2*div_r)[...,None,None]*Grr
                          +div[...,None,None]*Gww
                          -dflow_dw[...,None,None]*(gram(basis_y,basis_w)+gram(basis_w,basis_y)))
            volumes[REGIONS[2*slab+side]]=(cell_volume.ravel(), cell_volume_derivative.ravel())
            pieces.append((REGIONS[2*slab+side],corners,K.reshape(-1,4,4),dK.reshape(-1,4,4)))

    for cap in (False,True):
        if not sheet:
            continue
        if cap:
            inner,bulk=indices[1,0][-1,:], indices[2,0][0,:]
        else:
            inner,bulk=indices[1,0][:,-1], indices[1,1][:,0]
        corners=np.stack((inner[:-1],inner[1:],bulk[:-1],bulk[1:]),axis=-1)
        K=np.zeros((n,4,4))
        dK=np.zeros_like(K)
        for t,wt in zip(roots,weights):
            jump=np.array([t-1,-t,1-t,t])
            if cap:
                y=radial[0][:-1]+t*np.diff(radial[0])
                r,dr=shape.exit*y,shape.exit*np.diff(radial[0])
                darea_t=3*(up.exit-down.exit)/(2*delta*shape.exit)
                factor=4*np.pi*r*r*dr
            else:
                w=axial[1][:-1]+t*np.diff(axial[1])
                r,rw=shape.radius(w)
                dr=(up.radius(w)[0]-down.radius(w)[0])/(2*delta)
                drw=(up.radius(w)[1]-down.radius(w)[1])/(2*delta)
                factor=4*np.pi*r*r*np.hypot(1.,rw)*np.diff(axial[1])
                darea_t=2*dr/r+rw*drw/(1+rw*rw)
            local= (wt*factor/sheet)[:,None,None]*np.outer(jump,jump)
            K+=local
            dK+=np.asarray(darea_t)[...,None,None]*local
        pieces.append(('cap_layer' if cap else 'side_layer',corners,K,dK))
    count=len(nodes)
    rows=np.concatenate([np.broadcast_to(ids[:,:,None],mat.shape).ravel() for _,ids,mat,_ in pieces])
    cols=np.concatenate([np.broadcast_to(ids[:,None,:],mat.shape).ravel() for _,ids,mat,_ in pieces])
    data=np.concatenate([mat.ravel() for _,_,mat,_ in pieces])
    stiffness=sparse.coo_matrix((data,(rows,cols)),shape=(count,count)).tocsr()
    feed=indices[0,1][:,-1]
    bulk=np.unique(np.r_[indices[1,1][:,-1], indices[2,1][:,-1], indices[2,1][-1],indices[2,0][-1]])
    free=np.ones(count,bool)
    free[feed]=free[bulk]=False
    if return_system:
        return pieces, volumes, indices, stiffness, feed, bulk, free
    theta=np.zeros(count)
    theta[bulk]=1.
    theta[free]=spsolve(stiffness[free][:,free],-(stiffness@theta)[free])
    port=stiffness@theta
    # Use a differential within each element before quadratic forms to avoid
    # subtracting near-equal energies from a large uniform reservoir phase.
    energy,work={},{}
    crossing={}
    for region,ids,K,dK in pieces:
        local=theta[ids]-theta[ids][:,:1]
        energy[region]=float(.5*np.einsum('ni,nij,nj->',local,K,local))
        work[region]=float(.5*np.einsum('ni,nij,nj->',local,dK,local))
        if region.endswith('layer'):
            crossing[region]=float(np.einsum('nij,nj->',K[:,:2,:],local))
    conductance=2*sum(energy.values())
    if max(abs(-port[feed].sum()/conductance-1),abs(port[bulk].sum()/conductance-1))>5.e-7:
        raise ArithmeticError('Energy and two port fluxes disagree')
    if sheet and abs(-sum(crossing.values())/conductance-1)>5.e-7:
        raise ArithmeticError('Side and cap fluxes did not close the receiving balance')
    return dict(step=step,mesh=n,linear_phase_bias=BIAS,sheet_length=sheet,outer=outer,roof=roof,
                theta_energy=sum(energy.values())*BIAS**2,port_current=conductance*BIAS,
                relative_port_mismatch=abs(port[feed].sum()+port[bulk].sum())/conductance,
                component_energy={k:v*BIAS**2 for k,v in energy.items()},
                component_work={k:v*BIAS**2 for k,v in work.items()},
                work=sum(work.values())*BIAS**2,
                outward_current={k:-v*BIAS for k,v in crossing.items()})


def run(a):
    rows=[assemble(s,a.mesh,a.sheet,a.outer,a.roof) for s in (-.01,-.005,0.,.005,.01)]
    print(json.dumps(dict(inputs=vars(a),states=rows),allow_nan=False),flush=True)
    checks=[]
    for i,h in ((0,.01),(1,.005)):
        m,z,p=rows[i],rows[2],rows[4-i]
        slope=(p['theta_energy']-m['theta_energy'])/(2*h)
        curvature=(p['theta_energy']+m['theta_energy']-2*z['theta_energy'])/h**2
        parts={key:(p['component_work'][key]-m['component_work'][key])/(2*h)
               for key in z['component_work']}
        print(json.dumps(dict(half_width=h,slope=slope,direct_work=z['work'],stiffness=curvature,
                              stiffness_from_work=sum(parts.values()),components=parts)),flush=True)
        checks.append((sum(parts.values())-curvature,z['work']-slope))
    richardson=(4*np.asarray(checks[1])-np.asarray(checks[0]))/3
    if abs(richardson[0])>2e-4 or abs(richardson[1])>1e-5:
        raise ArithmeticError('Port energy and integrated variational stress disagree after step extrapolation')


if __name__=='__main__':
    ap=argparse.ArgumentParser()
    ap.add_argument('--mesh',type=int,required=True)
    ap.add_argument('--sheet',type=float,default=.04/6)
    ap.add_argument('--outer',type=float,default=28.)
    ap.add_argument('--roof',type=float,default=22.)
    run(ap.parse_args())
