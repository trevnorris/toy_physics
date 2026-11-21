"""Leading change of passive conversion power from the strain dependence of order.

Physical gradients of bulk and wall energies are taken at fixed material strain
and contact reference. The complementary curl field gives those physical strains.
The monopole suffices for the integrated power derivative about a round mouth.
"""
import argparse
import json
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy import sparse
from scipy.linalg import solve
from scipy.sparse.linalg import splu

from mouth_bend_renewal import ROOT, volume_cost, wall_data, gram
from mouth_bend_static_shear import ANGULAR, archive_inputs, stiffness_matrix
from mouth_carrier_lateral_response import prepare
from mouth_curved_phase import quadratic,tensor
from mouth_curved_phase_stability import weighted_gram
from mouth_empty_bend_stiffness import quadrature
from mouth_order_wave_pulse import mesh_geometry, product, linear, emit


def patch(geo, block, grid, order, sub, quad, floor):
    pg,vg,nodes=grid
    _,vo,old=geo.order_blocks[block]
    i,j=[a.ravel() for a in np.indices((len(pg)-1,len(vg)-1))]
    parent=old[(i//sub)*(len(vo)-1)+j//sub]
    root,weights=leggauss(quad)
    xi,eta=[v.ravel() for v in np.meshgrid(root,root,indexing='ij')]
    def original_basis(x,y):
        return tensor(quadratic(2*((i%sub)[:,None]+(x+1)/2)/sub-1)[0],
                      quadratic(2*((j%sub)[:,None]+(y+1)/2)/sub-1)[0])
    probe=np.array([-1.,0.,1.])
    xprobe,yprobe=[v.ravel() for v in np.meshgrid(probe,probe,indexing='ij')]
    to_bernstein=np.array([[1.,0,0],[-.5,2.,-.5],[0,0,1.]])
    sample=(original_basis(xprobe,yprobe)@order[parent][...,None])[...,0].reshape(-1,3,3)
    live=(to_bernstein@sample@to_bernstein.T).min(axis=(1,2))>floor
    orig=original_basis(xi,eta)
    chi=np.clip((orig@order[parent][...,None])[...,0],0,1)
    hp,hv=np.diff(pg)[i,None]/2,np.diff(vg)[j,None]/2
    r,_,rp,rv,zp,zv=geo.chart(block,pg[i,None]+hp*(1+xi),vg[j,None]+hv*(1+eta))
    det=rp*zv-rv*zp
    derivative=np.broadcast_to([-.5,.5],(len(xi),2))
    dp=product(derivative,linear(eta))/hp[:,:,None]
    dv=product(linear(xi),derivative)/hv[:,:,None]
    dr=(zv[:,:,None]*dp-zp[:,:,None]*dv)/det[:,:,None]
    dw=(rp[:,:,None]*dv-rv[:,:,None]*dp)/det[:,:,None]
    val=np.broadcast_to(product(linear(xi),linear(eta)),dr.shape)
    measure=abs(det)*hp*hv*np.outer(weights,weights).ravel()
    inverse=np.divide(1.,chi,out=np.zeros_like(chi),where=live[:,None])
    return parent,orig,live,r,dr,dw,val,measure,inverse


def energy_sources(geo,order,psi,edges,boundary_keep,boundary_modes,args):
    n=args.sub
    points,grids=mesh_geometry(geo,geo.grids,n,False)
    matrix=sparse.csc_matrix((len(points),len(points)))
    active=np.ones(len(points),bool)
    sides=[]
    for block,grid in enumerate(grids):
        pg,vg,nodes=grid
        _,_,live,r,dr,dw,val,measure,inv=patch(geo,block,grid,order,n,args.quad,args.floor)
        active[nodes[~live].ravel()]=False
        A=measure*inv
        stiffness=(gram(A,dr,dr)+gram(A,dw,dw))/2+gram(A/r**2,val,val)
        matrix+=stiffness_matrix(stiffness,nodes,len(points))
        boxes=nodes.reshape(len(pg)-1,len(vg)-1,4)
        sides.extend((boxes[0,:,[0,1]].T,boxes[-1,:,[2,3]].T,boxes[:,0][:,[0,2]],boxes[:,-1][:,[1,3]]))
        if block==2:
            wall_pairs=boxes[:,0][:,[0,2]]
            wall=np.r_[wall_pairs[:,0],wall_pairs[-1,1]]
    pairs,count=np.unique(np.sort(np.vstack(sides),axis=1),axis=0,return_counts=True)
    allowed={tuple(pair) for pair in np.sort(wall_pairs,axis=1)}
    for pair in pairs[count==1]:
        xy=points[pair]
        if tuple(pair) not in allowed and not (np.all(xy[:,0]==geo.far) or np.all(xy[:,1]==geo.roof)):
            active[pair]=False
    active[points[:,0]==0]=False
    np.testing.assert_array_equal(active[wall],boundary_keep)
    interior=active.copy()
    interior[wall]=False
    m=np.zeros((len(points),2))
    m[wall[boundary_keep]]=boundary_modes
    inner=matrix[interior][:,interior]
    scale=sparse.diags(1/np.sqrt(inner.diagonal()))
    m[interior]=-scale@splu((scale@inner@scale).tocsc()).solve(scale@(matrix@m)[interior])
    bulk_stored=ANGULAR/2*np.sum(m*(matrix@m),axis=0)
    bulk_source=np.zeros((len(order),2))
    for block,grid in enumerate(grids):
        nodes=grid[2]
        parent,orig,_,r,dr,dw,val,measure,inv=patch(geo,block,grid,order,n,args.quad,args.floor)
        phase=np.angle((orig@psi[parent][...,None])[...,0])
        for col in range(2):
            local=m[nodes,col,None]
            strain_squared=(np.square(dr@local)+np.square(dw@local))[...,0]/2+np.square((val@local)[...,0]/r)
            derivative=ANGULAR/2*strain_squared*measure*inv**2*np.exp(-phase/.4)
            np.add.at(bulk_source[:,col],parent,np.einsum('eq,eqi->ei',derivative,orig))
    # Same two-sphere angular integral for tangential meridional and polar wall strains.
    xi,wt=leggauss(args.quad)
    length=np.diff(edges)[:,None]
    z=edges[:-1,None]+length*(1+xi)/2
    po,vo,old=geo.order_blocks[2]
    parent=old[(np.arange(len(length))//n)*(len(vo)-1)][:,[0,3,6]]
    q=quadratic(2*((np.arange(len(length))%n)[:,None]+(xi+1)/2)/n-1)[0]
    chi=(q@order[parent][...,None])[...,0]
    weights=length*wt/2
    phase=np.angle((q@psi[parent][...,None])[...,0])
    r,rp=geo.radius(z)
    arc=np.hypot(1.,rp)
    wall_source=np.zeros_like(bulk_source)
    wall_stored=[]
    for col in range(2):
        mode=m[wall,col]
        at=mode[:-1,None]+np.diff(mode)[:,None]*(1+xi)/2
        slope=np.diff(mode)[:,None]/length
        valid=(at!=0)|(slope!=0)
        inv=np.divide(1.,chi,out=np.zeros_like(chi),where=valid)
        elastic=ANGULAR/2*weights*(arc*(at/r)**2+slope**2/(2*arc))
        wall_stored.append(np.sum(elastic*inv))
        derivative=elastic*inv**2*np.exp(-phase/.4)
        np.add.at(wall_source[:,col],parent,np.einsum('eq,eqi->ei',derivative,q))
    return bulk_source,wall_source,bulk_stored,np.array(wall_stored)


def run(a):
    order_file=ROOT/f'mouth_free_order_delayed_n{a.mesh}_d04/order_seed1.55.npz'
    geo,order,_=archive_inputs(order_file,ROOT/f'bend_static_shear_delcarrier{a.mesh}_01')
    meta=json.loads((order_file.parent/'stdout').read_text().splitlines()[0])
    a.source=Path(meta['source']).parent
    dom,psi=prepare(a)
    np.testing.assert_array_equal(dom.points,geo.points)
    nu=meta['nu']
    if nu!=.4:
        raise ValueError('Wrong curvature mobility')
    kappa=nu/meta['mobility']
    edges,keep,bulk=volume_cost(geo,order,a.sub,a.quad,a.floor)
    contact,load,B,G,drift,clock,_=wall_data(geo,order,psi,edges,keep,a.quad)
    operator=bulk+contact+B@solve(clock*drift,G)
    scale=1/np.sqrt((bulk+contact).diagonal())
    inlet=np.zeros(len(drift))
    inlet[len(inlet)//2:]=1.
    straight,influence=(scale[:,None]*solve(scale[:,None]*operator*scale,
                            scale[:,None]*np.column_stack([load[:,1],B@inlet]))).T
    modes=np.column_stack([straight,straight-(straight[0]/influence[0])*influence])
    body,wall,body_energy,wall_energy=energy_sources(geo,order,psi,edges,keep,modes,a)
    for matrix,estimate in ((bulk,body_energy),(contact,wall_energy)):
        np.testing.assert_allclose(estimate,ANGULAR/2*np.sum(modes*(matrix@modes),axis=0),rtol=2e-8)
    _,_,dr,dw=quadrature(dom)
    chi=dom.N@order
    f=dom.N@psi
    ur,uw=(np.imag(d@(psi-1)/f) for d in (dr,dw))
    adv=ur*(dr@order)+uw*(dw@order)
    w=dom.weight*np.exp(-np.angle(f)/nu)
    H=weighted_gram(dr,w)+weighted_gram(dw,w)+weighted_gram(dom.N,w*(1-6*chi*(1-chi))/meta['delta']**2)
    factor=splu(H[dom.free][:,dom.free].tocsc())
    loss_slope=(dr.T@(dom.weight*adv*ur)+dw.T@(dom.weight*adv*uw))*2/meta['mobility']
    adjoint=factor.solve(loss_slope[dom.free],trans='T')
    for j,name in enumerate(('straight','traction_free')):
        rhs=-np.column_stack([body[:,j],wall[:,j]])/kappa
        shifts=np.zeros((len(order),2))
        shifts[dom.free]=factor.solve(rhs[dom.free])
        derivative=loss_slope@shifts
        rate_shift=ur[:,None]*(dr@shifts)+uw[:,None]*(dw@shifts)
        direct=2/meta['mobility']*((dom.weight*adv)@rate_shift)
        np.testing.assert_allclose(direct,derivative,rtol=1e-9,atol=3e-12)
        np.testing.assert_allclose(adjoint@rhs[dom.free],derivative,rtol=1e-9,atol=3e-12)
        central=(chi>.05)&(chi<.95)
        order_shift=(dom.N@shifts).sum(axis=1)*.01**2
        projected=chi+order_shift
        emit(test='order_feedback',mesh=a.mesh,sub=a.sub,floor=a.floor,inlet=name,
             derivative_loss_per_bend_squared=derivative.tolist(),
             change_in_loss_at_b01=float(derivative.sum()*.01**2),
             term_order=['volume_shear','wall_attachment'],renewal_at_b01=float(2*wall_energy[j]/clock*.01**2),
             stored_bulk_wall_per_bend_squared=[body_energy[j],wall_energy[j]],
             peak_monopole_shift_at_b01=float(np.max(abs(order_shift))),
             peak_shift_inside_transition=float(np.max(abs(order_shift[central]))),
             baseline_order_range=[float(chi.min()),float(chi.max())],
             extrapolated_order_range=[float(projected.min()),float(projected.max())],
             min_shifted_inside_transition=float(projected[central].min()),
             highest_order_where_shift_exceeds_base=float(chi[(abs(order_shift)>abs(chi))].max(initial=0)),
             norm_volume_wall_sources=[float(np.linalg.norm(rhs[:,i])) for i in (0,1)],
             adjoint_check=(adjoint@rhs[dom.free]).tolist())


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--mesh',type=int,choices=(34,48),required=True)
    p.add_argument('--sub',type=int,required=True)
    p.add_argument('--quad',type=int,default=4)
    p.add_argument('--floor',type=float,default=1e-9)
    run(p.parse_args())
