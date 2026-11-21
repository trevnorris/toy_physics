"""Drive the material torsional field, switch off the torque, and measure its ringdown.

The material is at rest for transverse propagation on a prescribed draining
carrier/order state. Three brane coordinates and w are retained through exact
angular integration. Constant transverse inertia, stiffness proportional to chi.
"""
import argparse
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy import sparse

from mouth_order_wave_pulse import (ARCHIVE, ROOT, EXIT, assemble, digest, emit, energies,
                                    exit_transit_bound, linear, product, quadratic, tensor)
from mouth_shear_ringdown_state import domain


def prepared(a):
    dom=domain(a.radius,a.mesh)
    if a.state:
        header=[json.loads(x) for x in (a.state/'stdout').read_text().splitlines()]
        saved=header[-1]
        if (not json.loads((a.state/'limit-validation.json').read_text())['verified']
                or json.loads((a.state/'outcome.json').read_text())['exitCode']
                or header[0]['radius']!=a.radius or header[0]['mesh']!=a.mesh):
            raise ValueError('Unverified or mismatched carrier and layer')
        path=a.state/'state.npz'
        if digest(path)!=saved['sha256']:
            raise ValueError('State digest mismatch')
        with np.load(path) as field:
            np.testing.assert_array_equal(dom.points,field['points'])
            chi=field['order']
        intake=header[0]['feed']
    else:
        if a.radius!=1:
            raise ValueError('Native archive has unit radius')
        archive=json.loads(ARCHIVE.read_text())
        c=ROOT/f'mouth_flared_cubic_n{a.mesh}'
        carrier=c/'field.npz'
        if digest(carrier)!=archive['input_carriers_sha256'][str(carrier)]:
            raise ValueError('Carrier digest mismatch')
        with np.load(carrier) as field:
            np.testing.assert_allclose(dom.points,field['points'],rtol=0,atol=1.1e-11)
        path=ROOT/f'mouth_free_order_cubic_n{a.mesh}_d04'/'order_seed1.55.npz'
        if digest(path)!=archive['diffuse_orders_sha256'][str(path)]:
            raise ValueError('Order digest mismatch')
        with np.load(path) as field:
            chi=field['chi']
        intake=json.loads((c/'stdout').read_text().splitlines()[0])['feed']
    order_blocks=list(dom.blocks.values())
    grids=[(p,v) for p,v,_ in order_blocks]
    # Continue only the ordered feeding slab beyond its original r=9 edge;
    # carrier and order in the original calculation are exactly as prepared.
    pg=np.linspace(9.,a.outer,int(np.ceil((a.outer-9.)/a.feed_spacing))+1)
    vg=grids[1][1]
    grids.append((pg,vg))
    order_blocks.append((pg,vg,np.full(((len(pg)-1)*(len(vg)-1),9),len(chi),dtype=int)))
    def chart(block,p,v):
        if block==5:
            one=np.ones_like(p)
            zero=np.zeros_like(p)
            return p,-v,one,zero,zero,-one
        return dom.chart(block,p,v)
    config=SimpleNamespace(grids=grids,order_blocks=order_blocks,chart=chart,far=a.outer)
    return config,np.r_[chi,1.],intake,path


def axial_stiffness(dom,blocks,chi,count):
    """Physical-depth shear-strain energy of the confined region alone."""
    roots,weights=leggauss(3)
    xi,eta=(x.ravel() for x in np.meshgrid(roots,roots,indexing='ij'))
    deriv=np.tile([-.5,.5],(len(xi),1))
    sp,sv=product(deriv,linear(eta)),product(linear(xi),deriv)
    pieces=[]
    for block in (0,2):
        pg,vg,ids=blocks[block]
        po,vo,old=dom.order_blocks[block]
        rows,cols=(x.ravel() for x in np.indices((len(pg)-1,len(vg)-1)))
        for start in range(0,len(rows),1536):
            i,j=rows[start:start+1536],cols[start:start+1536]
            hp,hv=np.diff(pg)[i,None]/2,np.diff(vg)[j,None]/2
            p,v=pg[i,None]+hp*(1+xi),vg[j,None]+hv*(1+eta)
            oi=np.searchsorted(po,(pg[i]+pg[i+1])/2)-1
            oj=np.searchsorted(vo,(vg[j]+vg[j+1])/2)-1
            ax=quadratic(2*(p-po[oi,None])/np.diff(po)[oi,None]-1)[0]
            ay=quadratic(2*(v-vo[oj,None])/np.diff(vo)[oj,None]-1)[0]
            order=np.clip(np.einsum('eqi,ei->eq',tensor(ax,ay),
                          chi[old[oi*(len(vo)-1)+oj]]),0,1)
            r,w,rp,rv,wp,wv=dom.chart(block,p,v)
            det=rp*wv-rv*wp
            dz=(rp[:,:,None]*sv/hv[:,:,None]-rv[:,:,None]*sp/hp[:,:,None])/det[:,:,None]
            wt=(8*np.pi/3)*r**4*abs(det)*hp*hv*np.outer(weights,weights).ravel()*order*(w>0)*(w<EXIT)
            pieces.append((ids[start:start+1536],np.einsum('eq,eqi,eqj->eij',wt,dz,dz)))
    ids=np.concatenate([q for q,_ in pieces])
    values=np.concatenate([v for _,v in pieces])
    mat=sparse.coo_matrix((values.ravel(),(np.repeat(ids,4,axis=1).ravel(),
                          np.tile(ids,(1,4)).ravel())),shape=(count,count)).tocsr()
    mat.eliminate_zeros()
    return mat


def spot(points,radius,center,halfwidth):
    r,w=points.T
    return (np.maximum(0,1-(r/(.47*radius))**2)*
            np.maximum(0,1-((w-center)/halfwidth)**2))**3


def run(a):
    dom,chi,intake,state=prepared(a)
    emit(test='preparation',radius=a.radius,mesh=a.mesh,intake=intake,order_file=str(state),
         order_hash=digest(state),outer=a.outer,**exit_transit_bound(dom,chi))
    pts,blocks,regional_m,parts,damp,_=assemble(dom,chi,a.sub,a.tail,a.feed_spacing/a.radius)
    mass=regional_m.sum(axis=0)
    k=sum(parts).tocsr()
    kw=axial_stiffness(dom,blocks,chi,len(mass))
    spatial=spot(pts,a.radius,.30,.16)
    force=mass*spatial
    force/=force.sum()
    receiver=spot(pts,a.radius,.88,.16)*mass
    receiver/=receiver.sum()
    drives=np.flatnonzero(force)
    receivers=np.flatnonzero(receiver)
    remote=np.flatnonzero(damp)
    if len(remote)==0 or np.min(pts[remote,0])!=a.outer:
        raise ValueError('Missing outer radiation load')
    sample=.04
    bound=np.max(np.asarray(abs(k).sum(axis=1)).ravel()/mass)
    stride=int(np.ceil(sample*np.sqrt(bound)/1.4))
    dt=sample/stride
    end=a.end
    nt=round(end/dt)
    if nt%stride or not np.isclose(nt*dt,end):
        raise ValueError('Sampling grid must cover the measurement window')
    duration=a.duration
    lo=1e-4
    def torque(t):
        if not 0<t<duration:
            return 0.
        envelope=np.sin(np.pi*t/duration)**2
        return lo*envelope*(np.cos(a.frequency*(t-duration/2)) if a.frequency else 1.)
    emit(test='drive',dofs=len(pts),dt=dt,radius=a.radius,frequency=a.frequency,
         duration=duration,end=end,torque_scale=lo,mesh=a.mesh,sub=a.sub,tail=a.tail,
         cfl_bound=bound,feed_spacing=a.feed_spacing,
         causal_return_lower_bound=2*(a.outer-.47*a.radius))
    u=np.zeros_like(mass)
    v=np.zeros_like(mass)
    inv=dt/(mass+.5*dt*damp)
    stay=(mass-.5*dt*damp)/(mass+.5*dt*damp)
    f=force[drives]
    b=damp[remote]
    work=loss=0.
    rows=[]
    for step in range(nt):
        field_force=k@u
        applied=torque(step*dt)
        field_force[drives]-=applied*f
        next_v=stay*v-inv*field_force
        if applied:
            work+=dt*applied*.5*(f@(v[drives]+next_v[drives]))
        loss+=dt*(b@((v[remote]+next_v[remote])*.5)**2)
        next_u=u+dt*next_v
        if (step+1)%stride==0:
            kin,strain=energies(regional_m,parts,next_v,u,next_u)
            energy=kin+strain
            axial=.5*u@(kw@next_u)
            wall=kin[1:5].sum()-strain[1:5].sum()+2*axial
            rows.append([(step+.5)*dt,work,loss,*energy,wall,force[drives]@next_v[drives],
                         receiver[receivers]@next_v[receivers],(energy.sum()+loss)-work])
        u,v=next_u,next_v
    cols=['time','work','loss','feed','ordered','upper','lower','tail','receiver',
          'axial_force','near_velocity','far_velocity','balance']
    trace=np.asarray(rows)
    out=a.outdir/'trace.npz'
    np.savez_compressed(out,columns=cols,trace=trace)
    maximum=float(np.max(abs(trace[:,-1]))/work)
    emit(test='completed',path=str(out),sha256=digest(out),input_work=work,
         max_balance_over_final_work=maximum)
    if work<=0 or not np.isfinite(trace).all() or maximum>1e-8:
        raise ValueError('Forced energy balance failed')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--radius',type=float,required=True)
    p.add_argument('--mesh',type=int,required=True)
    p.add_argument('--state',type=Path)
    p.add_argument('--sub',type=int,default=1)
    p.add_argument('--tail',type=int,required=True)
    p.add_argument('--feed-spacing',type=float,default=.055)
    p.add_argument('--outer',type=float,default=21.)
    p.add_argument('--frequency',type=float,default=0.)
    p.add_argument('--duration',type=float,default=.48)
    p.add_argument('--end',type=float,default=20.)
    p.add_argument('--outdir',type=Path,required=True)
    run(p.parse_args())
