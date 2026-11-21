"""Measure virtual work and wave momentum from a local burst after switch-off.
The first variation at zero offset is evolved independently of the pole solve.
"""
import argparse

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy import sparse

from mouth_native_wave_bend import (Bend,bend_forms,emit,geometry,linear,product,
                                    quadratic,taper,tensor,x_derivative)
from mouth_order_wave_pulse import mesh_geometry


def momentum_forms(g,blocks,interval):
    """Order-weighted integrals of H(r N_i) times d_x(r N_j T_j), volume and cut."""
    q,weights=leggauss(3)
    xi,eta=(v.ravel() for v in np.meshgrid(q,q,indexing='ij'))
    n=product(linear(xi),linear(eta))
    der=np.tile([-.5,.5],(len(xi),1))
    dp_ref,dv_ref=product(der,linear(eta)),product(linear(xi),der)
    pieces=[[] for _ in range(6)]
    def add(group,ids,wt,left,right):
        val=np.einsum('eq,eqi,eqj->eij',wt,left,right)
        pieces[group].append((ids,val))
    for ib,(pg,vg,cells) in enumerate(blocks):
        if ib in (3,4): continue
        po,vo,old=g.dom.order_blocks[ib]
        ii,jj=(v.ravel() for v in np.indices((len(pg)-1,len(vg)-1)))
        for start in range(0,len(ii),1000):
            sl=slice(start,start+1000)
            i,j,ids=ii[sl],jj[sl],cells[sl]
            hp,hv=np.diff(pg)[i,None]/2,np.diff(vg)[j,None]/2
            p,v=pg[i,None]+hp*(1+xi),vg[j,None]+hv*(1+eta)
            oi=np.searchsorted(po,(pg[i]+pg[i+1])/2)-1
            oj=np.searchsorted(vo,(vg[j]+vg[j+1])/2)-1
            ax=quadratic(2*(p-po[oi,None])/np.diff(po)[oi,None]-1)[0]
            ay=quadratic(2*(v-vo[oj,None])/np.diff(vo)[oj,None]-1)[0]
            order=np.clip(np.einsum('eqi,ei->eq',tensor(ax,ay),g.chi[old[oi*(len(vo)-1)+oj]]),0,1)
            r,w,rp,rv,wp,wv=g.dom.chart(ib,p,v)
            det=rp*wv-rv*wp
            rgrad=r[:,:,None]*(wv[:,:,None]*dp_ref/hp[:,:,None]-wp[:,:,None]*dv_ref/hv[:,:,None])/det[:,:,None]
            rn=r[:,:,None]*n
            ns=np.broadcast_to(n,rn.shape)
            wt=(8*np.pi/3)*r**2*abs(det)*hp*hv*order*np.outer(weights,weights).ravel()
            entrance=taper(w,*interval)[0]
            for s,h in enumerate((entrance,1-entrance)):
                if np.any(wt*h):
                    add(2*s,ids,wt*h,rn,rgrad)
                    add(2*s+1,ids,wt*h,rn,ns)
    # Plane w=cut. The radial derivative lies within that plane. No n_x term.
    pg,vg,cells=blocks[2]
    po,vo,old=g.dom.order_blocks[2]
    row=np.searchsorted(pg,g.cut)-1
    j=np.arange(len(vg)-1)
    hv=np.diff(vg)[:,None]/2
    v=vg[j,None]+hv*(1+q)
    p=np.full_like(v,g.cut)
    r,_,_,rv,_,_=g.dom.chart(2,p,v)
    oj=np.searchsorted(vo,(vg[:-1]+vg[1:])/2)-1
    oi=np.searchsorted(po,g.cut)-1
    oy=quadratic(2*(v-vo[oj,None])/np.diff(vo)[oj,None]-1)[0]
    order=np.clip(np.einsum('eqi,ei->eq',oy, g.chi[old.reshape(len(po)-1,len(vo)-1,3,3)[oi,oj,-1,:]]),0,1)
    ns=np.broadcast_to(linear(q),(*r.shape,2))
    dr=r[:,:,None]*np.array([-.5,.5])/(hv*rv)[:,:,None]
    wt=(8*np.pi/3)*r**2*abs(rv)*hv*order*weights
    ids=cells.reshape(len(pg)-1,len(vg)-1,2,2)[row,:,-1,:]
    add(4,ids,wt,r[:,:,None]*ns,dr)
    add(5,ids,wt,r[:,:,None]*ns,ns)
    mats=[]
    for group in pieces:
        ids=np.concatenate([entry[0] for entry in group])
        vals=np.concatenate([entry[1] for entry in group])
        nloc=ids.shape[1]
        mat=sparse.coo_matrix((vals.ravel(),(np.repeat(ids,nloc,axis=1).ravel(),np.tile(ids,(1,nloc)).ravel())),shape=(g.size,g.size)).tocsr()
        mats.append(g.restrict(mat))
    return mats


def run(a):
    g=geometry(a)
    _,blocks=mesh_geometry(g.dom,g.dom.grids,a.sub,False,feed_spacing=a.feed_spacing)
    interval=(.55,1.15) if a.profile=='collar' else (1.2,2.)
    b=Bend([g.restrict(m) for m in bend_forms(g.dom,g.block,g.chi,g.size,interval)])
    channels=(3,4,5)
    k=[g.k+(l*(l+1)-2)*g.ang for l in channels]
    adjacent=[b.linear(l,4) for l in (3,5)]
    curvature=b.square(4,4)
    forms=momentum_forms(g,blocks,interval)
    pairs=((0,1),(2,1),(1,0),(1,2))
    angular_ports=[]
    for i,j in pairs:
        c,coeff=x_derivative(channels[i],channels[j])
        angular_ports.append([c*(forms[s]+coeff*forms[s+1]) for s in (0,2,4)])
    r,w=g.points.T
    force=np.sqrt(g.mass)*r**3*np.exp(-(r/.7)**2)*np.maximum(0.,1-((w-.95)/.18)**2)**3
    force/=np.linalg.norm(force)
    forced=np.flatnonzero(force)
    outer=g.points[:,0]==a.outer
    sample=.12
    limit=max(np.max(np.asarray(abs(m).sum(axis=1))) for m in k)
    stride=int(np.ceil(sample*np.sqrt(limit)/1.45))
    dt=sample/stride
    inv=dt/(1+.5*dt*g.damp)
    ratio=(1-.5*dt*g.damp)/(1+.5*dt*g.damp)
    u=np.zeros((3,len(r))); vel=np.zeros_like(u)
    work=absorbed=far_absorbed=0.
    duration=12*np.pi/6.446
    trace=[]
    for step in range(round(a.end/sample)*stride):
        t=step*dt
        residual=np.array([m@v for m,v in zip(k,u,strict=True)])
        residual[0]+=adjacent[0]@u[1]
        residual[2]+=adjacent[1]@u[1]
        amp=1e-4*np.sin(np.pi*t/duration)**2*np.cos(6.446*(t-duration/2)) if t<duration else 0.
        if amp: residual[1,forced]-=amp*force[forced]
        vn=ratio*vel-inv*residual
        un=u+dt*vn
        midpoint=.5*(vel[1]+vn[1])
        if amp: work+=dt*amp*(force[forced]@midpoint[forced])
        absorbed+=dt*(g.damp@midpoint**2)
        far_absorbed+=dt*(g.damp[outer]@midpoint[outer]**2)
        if (step+1)%stride==0:
            mid=.5*(u+un)
            en=.5*(vn[1]@vn[1]+u[1]@(k[1]@un[1]))
            configurational=-mid[1]@(curvature@mid[1])
            for z,mat in zip((0,2),adjacent,strict=True):
                configurational-=mid[z]@(mat@mid[1])
            ports=np.zeros(3)
            for (i,j),port in zip(pairs,angular_ports,strict=True):
                for n,mat in enumerate(port):
                    ports[n]+=vn[i]@(mat@mid[j])
            trace.append([(step+.5)*dt-duration,work,en,absorbed,far_absorbed,configurational,*ports])
        u,vel=un,vn
    d=np.array(trace)
    err=np.max(abs(d[:,2]+d[:,3]-d[:,1]))/work
    emit(test='release',profile=a.profile,sub=a.sub,cut=g.cut,outer=a.outer,dt=dt,work=work,
         ledger_error=err,far_work_per_input=far_absorbed/work,omega_drive=6.446,
         time_end_after_switch=d[-1,0],intake=g.intake)
    if err>1e-8 or not np.isfinite(d).all(): raise ArithmeticError('Straight burst energy failed')
    hcut=1-taper(g.cut,*interval)[0]
    for lo,hi in ((8,24),(20,44),(32,54),(44,70),(56,80),(70,100)):
        if hi>d[-1,0]: continue
        sel=(d[:,0]>=lo)&(d[:,0]<=hi)
        phase=2*np.pi*(d[sel,0]-lo)/(hi-lo)
        weight=(1-np.cos(phase))/2
        deriv=np.pi*np.sin(phase)/(hi-lo)
        energy=weight@d[sel,2]
        geom=weight@d[sel,5]/energy
        # +d/dt int rho H q_t q_x, averaged by parts using a window
        # vanishing at both ends. The bath term is +int Z H q_t q_x.
        storage=-deriv@d[sel,6:8]/energy
        bath=(weight@d[sel,8])/energy*np.array([1-hcut,hcut])
        totals=np.array([-geom,geom])+storage+bath
        emit(test='window',window=[lo,hi],mean_energy_per_input=energy/(weight.sum()*work),
             deeper_config_force_per_energy_per_offset=geom,storage_correction=storage.tolist(),
             cut_bath_correction=bath.tolist(),both_reactions=totals.tolist(),
             sign_order='entrance, deeper flare; positive along the imposed deeper offset')
    if a.trace:
        np.savez_compressed(a.trace, values=d,columns=['time_after_off','work','energy','absorbed',
            'remote_work','config_flare','minus_momentum_entrance','minus_momentum_flare','cut_bath'])


if __name__=='__main__':
    from pathlib import Path
    p=argparse.ArgumentParser(description=__doc__)
    p.set_defaults(radius=1.,state=None,boundary='open',mesh=34)
    p.add_argument('--sub',type=int,default=1)
    p.add_argument('--cut',type=float,default=2.15)
    p.add_argument('--outer',type=float,default=91.)
    p.add_argument('--feed-spacing',type=float,default=.12)
    p.add_argument('--end',type=float,default=87.)
    p.add_argument('--profile',choices=('collar','deep'),required=True)
    p.add_argument('--trace',type=Path)
    a=p.parse_args()
    if a.outer<a.end+2: p.error('Remote termination must lie beyond the source light cone.')
    run(a)
