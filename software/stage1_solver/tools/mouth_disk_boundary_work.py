"""Independent boundary diagnostics on the same R4 fields and bonded lattice.

The direct cap Cauchy stress uses the deformed material chart on its held top
face. Separate FE tests decompose shared cap-node reactions by neighboring bank.
"""
import argparse
import copy
from pathlib import Path

import numpy as np

from mouth_bonded_front_adjoint import angular_weights
from mouth_bonded_front_feedback import Relocated
from mouth_disk_mesh_rematch import inputs, digest
from mouth_disk_pillbox import Pillbox, gauss_cells, lift, dot
from mouth_supercore_joined_front import emit, ORDER_INERTIA, KAPPA, WELL


def crystal_cap(base,pill,compact,n):
    m,body,d = pill.m,pill.b,pill.d
    chart=copy.copy(base)
    chart.sh=m.shape
    geom=Relocated(chart,dict(radial_fraction=body.x,depth_fraction=body.t,displacement=d.copy()))
    if not np.array_equal(base.co,base.co0):
        raise ValueError('The reference spacing used in the recovered core strain needs a separate evaluation')
    pressure=np.exp(base.co[0])
    # Split the radial integral at the faces of all three finite element grids.
    label_edges=np.exp(-geom.pair.a*pressure)
    edges=np.unique(np.r_[body.x,m.rr[0],label_edges[(label_edges<1)&(label_edges>0)],.75,.95])
    x,w=gauss_cells(edges,n)
    t=np.ones_like(x)
    r,z,motion,B,b,bw=geom.mapping(x,t)
    if abs(motion).max()>1e-14:raise ValueError('Direct crystal trace requires a held material cap')
    F,h,J=geom.initial_strain(x,t,r,b,bw)
    J*=np.linalg.det(B)
    F=B@F
    Kzz=dot(F[:,1],F[:,1])
    tr=(Kzz-1+.5*J*(J-1))/J
    raw=(pill.pair[0][2]-pill.pair[1][2])[m.ids[1,0][:,:,-1]]
    dc=np.interp(x,m.rr[0],raw@angular_weights(m.a))
    eta=lift(m.shape,r,z,.12,compact)[0]
    return float(-np.sum(w*4*np.pi*b*r*r*eta*dc*tr))


def cap_node_work(pill):
    m=pill.m
    wt,r,a,z,F,*_=m.D.T
    cap=np.unique(m.ids[1,0][:,:,-1])
    L=m.shape.length
    core=(z<L)&(r<m.shape.radius(z)[0])&(z>0)
    receiver=(z>L)&(r<m.shape.exit)
    masks=[core,receiver,~(core|receiver)]
    # Element test functions for a shared cap-node displacement of each
    # primary field. The rim belongs to the outer bank as well as the two caps.
    columns=[]
    for U in pill.pair:
        n=(F*(m.V@U[0]))**2
        chi=m.ref[:,0]+m.V@U[2]
        gp=np.stack([G@U[1] for G in m.G],axis=1)
        gc=m.ref[:,1:]+np.stack([G@U[2] for G in m.G],axis=1)
        den=n+ORDER_INERTIA*dot(gc,gc)
        coeff=ORDER_INERTIA*dot(gp,gc)/den
        vel=gp-coeff[:,None]*gc
        j=n[:,None]*vel
        pi=n*coeff
        gc_flux=KAPPA*gc-pi[:,None]*vel
        dw=2*WELL*chi*(1-chi)*(1-2*chi)
        banks=[]
        for take in masks:
            w=wt*take
            phase=sum(G.T@(w*j[:,k]) for k,G in enumerate(m.G))
            order=sum(G.T@(w*gc_flux[:,k]) for k,G in enumerate(m.G))+m.V.T@(w*dw)
            banks.append(np.stack([phase[cap],order[cap]]))
        columns.append(np.stack(banks))
    # The same fixed-lattice W load in both controls cancels from the
    # drive-minus-rest order reaction, even on prescribed order nodes.
    delta=columns[0]-columns[1]
    total=delta.sum(axis=0)
    oc=m.masks[2,cap]
    return dict(phase_bank_sums=delta[:,0].sum(axis=1).tolist(),phase_net=total[0].sum(),
                order_bank_sums=delta[:,1].sum(axis=1).tolist(),order_net=total[1].sum(),
                order_pinned_resultant=total[1,~oc].sum(),order_free_resultant=total[1,oc].sum(),
                phase_max_node=abs(total[0]).max(),order_free_max_node=abs(total[1,oc]).max())


def local_slopes(axis,values,location):
    i=np.clip(np.searchsorted(axis,location,side='right')-1,0,len(axis)-2)
    h=np.diff(axis)
    slopes=np.diff(values,axis=2)/h[None,None,:]
    if np.isclose(location,axis[i],atol=1e-13,rtol=0.) and i>0:
        return .5*(slopes[:,:,i-1]+slopes[:,:,i])
    return slopes[:,:,i]


def rim_order_work(pill,h,from_side):
    m=pill.m
    # A fixed-r axial lift through the rim needs a slope extension at the
    # junction. Use the inner disk incident/receiving limits separately.
    # At remaining cylinder nodes use the adjacent two cells' slope average.
    work=np.zeros(2)
    L=m.shape.length
    for slab in (1,2):
        for iz,z in enumerate(m.z[slab]):
            if abs(z-L)>=h-1e-12 or (slab==2 and iz==0):continue
            face=from_side if abs(z-L)<1e-12 else slab
            axis=m.z[face]
            ids=m.ids[slab,0][-1,:,iz]
            if np.any(m.masks[2,ids]):raise ValueError('Rim is not fully constrained')
            b,slope=m.shape.radius(z)
            if face==2:slope=0.
            eta=lift(m.shape,b,z,h,False)[0]
            j=np.argmin(abs(axis-z))
            for k,vals in enumerate(pill.faces.nodal[face]):
                c=vals[2]
                dr=(c[-1,:,j]-c[-2,:,j])/(np.diff(m.rr[0])[-1]*b)
                der=local_slopes(axis,c[-1:,:,:],z)[0]-slope*dr
                work[k]+=eta*(pill.order_grad[k][ids]@der)
    return work[0]-work[1]


def run(args):
    base,body,d,pair=inputs(args.state,args.joint,args.fields)
    model=base.joined(base.sh,False)
    pill=Pillbox(model,base.c,body,d,pair)
    emit(kind='preparation',source_sha256=digest(args.state),joint_sha256=digest(args.joint),
         fields_sha256=digest(args.fields),grid=base.grid,winding=base.N)
    emit(kind='shared_cap_variation',banks=['duct','over_disk','outside_rim'],**cap_node_work(pill))
    for compact in (False,True):
        for q in (6,10):
            emit(kind='sharp_solid',radial='interior_to_.95' if compact else 'full',quad=q,
                 crystal_axial=crystal_cap(base,pill,compact,q))
    for h in (.12*2/7,.12*4/7,.12):
        emit(kind='constrained_rim',h=h,gedot_chi_z_in= float(rim_order_work(pill,h,1)),
             gedot_chi_z_out=float(rim_order_work(pill,h,2)))
        for compact in (False,True):
            emit(**pill.audit(h,compact,8))
    emit(kind='surface_quadrature_replica',quad=12,readout=pill.sides(1,.12,False,12)[0].tolist(),
         receiving=pill.sides(2,.12,False,12)[0].tolist())


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__)
    for key in ('state','joint','fields'):ap.add_argument(key,type=Path)
    run(ap.parse_args())
