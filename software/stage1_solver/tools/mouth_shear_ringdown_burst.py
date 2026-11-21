"""Resonant burst with a separate readout of momentum during a virtual depth change."""
import argparse
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy import sparse

from mouth_shear_ringdown import (EXIT, assemble, axial_stiffness, digest, emit, energies,
                                 linear, prepared, product, quadratic, spot, tensor)


def momentum_and_mass(dom, blocks, chi, count, floor):
    """B_ij = int rho xi N_i d_w N_j dV with torsional angular and radial weights."""
    roots,weights=leggauss(3)
    xi,eta=(x.ravel() for x in np.meshgrid(roots,roots,indexing='ij'))
    shape=product(linear(xi),linear(eta))
    der=np.tile([-.5,.5],(len(xi),1))
    dp,dv=product(der,linear(eta)),product(linear(xi),der)
    mass=np.zeros((6,count))
    pieces=[]
    for ib,(pg,vg,ids) in enumerate(blocks):
        po,vo,old=dom.order_blocks[ib]
        ii,jj=(x.ravel() for x in np.indices((len(pg)-1,len(vg)-1)))
        for start in range(0,len(ii),1536):
            i,j=ii[start:start+1536],jj[start:start+1536]
            dofs=ids[start:start+1536]
            hp,hv=np.diff(pg)[i,None]/2,np.diff(vg)[j,None]/2
            p,v=pg[i,None]+hp*(1+xi),vg[j,None]+hv*(1+eta)
            oi=np.searchsorted(po,(pg[i]+pg[i+1])/2)-1
            oj=np.searchsorted(vo,(vg[j]+vg[j+1])/2)-1
            ax=quadratic(2*(p-po[oi,None])/np.diff(po)[oi,None]-1)[0]
            ay=quadratic(2*(v-vo[oj,None])/np.diff(vo)[oj,None]-1)[0]
            order=np.clip(np.einsum('eqi,ei->eq',tensor(ax,ay),chi[old[oi*(len(vo)-1)+oj]]),0,1)
            r,w,rp,rv,wp,wv=dom.chart(ib,p,v)
            det=rp*wv-rv*wp
            rho=floor+(1-floor)*order
            wt=(8*np.pi/3)*np.outer(weights,weights).ravel()*r**4*abs(det)*hp*hv*rho
            groups=np.select([w<0,w>=EXIT,order>=.9,order>=.5,order>=.1],[0,5,1,2,3],default=4)
            for region in np.unique(groups):
                local=(wt*(groups==region))@shape
                np.add.at(mass[region],dofs.ravel(),local.ravel())
            if ib in (1,5):
                continue
            dz=(rp[:,:,None]*dv/hv[:,:,None]-rv[:,:,None]*dp/hp[:,:,None])/det[:,:,None]
            temp=np.einsum('eq,qi,eqj->eij',wt*np.clip(w,0.,EXIT),shape,dz)
            pieces.append((dofs,temp))
    ids=np.concatenate([q for q,_ in pieces])
    values=np.concatenate([v for _,v in pieces])
    mat=sparse.coo_matrix((values.ravel(),(np.repeat(ids,4,axis=1).ravel(),
                          np.tile(ids,(1,4)).ravel())),shape=(count,count)).tocsr()
    mat.eliminate_zeros()
    return mass,mat


def run(a):
    dom,chi,intake,path=prepared(a)
    cut=None
    if a.cut is not None:
        from mouth_shear_cut import truncate,impedance
        dom,chi,cut=truncate(dom,chi,a.cut)
    pts,blocks,old_m,parts,damp,_=assemble(dom,chi,a.sub,a.tail,a.feed_spacing/a.radius)
    regional_m,bmat=momentum_and_mass(dom,blocks,chi,len(pts),a.floor)
    if a.floor==1:
        np.testing.assert_allclose(regional_m,old_m,rtol=4e-12,atol=4e-14)
    mass=regional_m.sum(axis=0)
    k=sum(parts).tocsr()
    kw=axial_stiffness(dom,blocks,chi,len(pts))
    tail_load=np.zeros(len(pts))
    cut_properties=None
    if cut is not None:
        cut_z,zmin,zmax=impedance(dom,blocks,pts,chi,cut,a.floor)
        if a.cut_load=='matched':
            tail_load=cut_z
        cut_properties=dict(w=cut,minimum_order=zmin,maximum_order=zmax,
                            termination=a.cut_load,matched_impedance=cut_z.sum())
        active=mass>0
        if np.any(k[~active].data) or np.any(tail_load[~active]):
            raise ValueError('Removed node carries shear stiffness or applied traction')
        pts,old_m,regional_m=pts[active],old_m[:,active],regional_m[:,active]
        mass,damp,tail_load=mass[active],damp[active],tail_load[active]
        parts=[op[active][:,active] for op in parts]
        k,kw,bmat=(op[active][:,active] for op in (k,kw,bmat))
    # Same spatial distribution of applied torque in both inertia laws.
    ref_mass=old_m.sum(axis=0)
    f=ref_mass*spot(pts,a.radius,a.source_depth,.16)
    f/=f.sum()
    near=ref_mass*spot(pts,a.radius,.30,.16)
    near/=near.sum()
    probe_indices=np.flatnonzero(near)
    g=ref_mass*spot(pts,a.radius,.88,.16)
    g/=g.sum()
    d=np.flatnonzero(f)
    o=np.flatnonzero(g)
    boundary=np.flatnonzero(damp)
    edge=np.flatnonzero(tail_load)
    sample=.04
    limit=np.max(np.asarray(abs(k).sum(axis=1)).ravel()/mass)
    stride=int(np.ceil(sample*np.sqrt(limit)/1.4))
    dt=sample/stride
    duration=2*np.pi*6/a.frequency
    emit(test='burst_preparation', radius=a.radius,mesh=a.mesh,sub=a.sub,tail=a.tail,
         feed_spacing=a.feed_spacing,outer=a.outer,order_file=str(path),order_hash=digest(path),
         current=intake,inertia_floor=a.floor,angular_frequency=a.frequency,
         drive_duration=duration,end=a.end,dt=dt,wave_nodes=len(pts),torque_amplitude=.0001,
         full_mass_law='floor+(1-floor)*chi', earliest_outer_echo_at_drive=2*(a.outer-.47*a.radius),cut=cut_properties,source_depth=a.source_depth)
    u=np.zeros_like(mass)
    vel=np.zeros_like(mass)
    damping=damp+tail_load
    inv=dt/(mass+.5*dt*damping)
    ratio=(mass-.5*dt*damping)/(mass+.5*dt*damping)
    work=loss=cut_export=0.
    time_series=[]
    for step in range(round(a.end/sample)*stride):
        t=step*dt
        applied=(1e-4*np.sin(np.pi*t/duration)**2*np.cos(a.frequency*(t-duration/2))
                 if t<duration else 0.)
        residual=k@u
        residual[d]-=applied*f[d]
        vn=ratio*vel-inv*residual
        if applied:
            work+=.5*dt*applied*(f[d]@(vel[d]+vn[d]))
        loss+=.25*dt*(damp[boundary]@(vel[boundary]+vn[boundary])**2)
        if len(edge):
            cut_export+=.25*dt*(tail_load[edge]@(vel[edge]+vn[edge])**2)
        un=u+dt*vn
        if (step+1)%stride==0:
            kinetic,strain=energies(regional_m,parts,vn,u,un)
            total=kinetic+strain
            frozen=kinetic[1:5].sum()-strain[1:5].sum()+u@(kw@un)
            transport=vn@(bmat@(.5*(u+un)))
            time_series.append([(step+.5)*dt,work,loss,*total,frozen,near[probe_indices]@vn[probe_indices],g[o]@vn[o],
                                total.sum()+loss+cut_export-work,transport,cut_export])
        u,vel=un,vn
    trace=np.array(time_series)
    columns=['time','work','loss','feed','ordered','upper','lower','tail','receiver',
             'axial_force','near_velocity','far_velocity','balance','transport_momentum','cut_export']
    out=a.outdir/'trace.npz'
    np.savez_compressed(out,columns=columns,trace=trace)
    err=np.max(abs(trace[:,columns.index('balance')]))/work
    emit(test='completed',path=str(out),sha256=digest(out),input_work=work,max_balance_over_final_work=err)
    if work<=0 or not np.isfinite(trace).all() or err>1e-8:
        raise ArithmeticError('Burst energy ledger failed')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--radius',type=float,required=True)
    p.add_argument('--mesh',type=int,default=34)
    p.add_argument('--state',type=Path)
    p.add_argument('--sub',type=int,required=True)
    p.add_argument('--tail',type=int,required=True)
    p.add_argument('--feed-spacing',type=float,default=.05)
    p.add_argument('--outer',type=float,default=21.)
    p.add_argument('--source-depth',type=float,default=.30)
    p.add_argument('--frequency',type=float,default=12.2)
    p.add_argument('--floor',type=float,default=1.)
    p.add_argument('--cut',type=float)
    p.add_argument('--cut-load',choices=('free','matched'),default='matched')
    p.add_argument('--end',type=float,default=22.)
    p.add_argument('--outdir',type=Path,required=True)
    args=p.parse_args()
    if not 0<=args.floor<=1 or (args.floor==0 and args.cut is None) or args.frequency<=0:
        p.error('Zero inertia floor requires an explicit tail cut; drive frequency must be positive')
    run(args)
