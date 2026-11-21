"""Direct meridional shear work during a weak prescribed stroke of the radius-one collar.

For X=(U e_r,V e_w), C=V_r-U_w. Write b=r^2 int mu*C dt=r^3 H.
From rest, (v_r,v_w)=(-b_w,b_r)/(rho*r^2). The natural boundary drive
for b is the material's tangential velocity. Work is 4 pi int b_t v_t ds.
Use mu=rho=chi, the archived shear-order map and the frozen drain geometry.
"""
import argparse
import json
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy import sparse

from mouth_order_wave_pulse import linear, mesh_geometry, product, quadratic, tensor
from mouth_shear_ringdown import digest, emit, prepared

CENTER = .88
WIDTH = .22
SAMPLE = .1


def profile(w):
    return np.maximum(0., 1-((w-CENTER)/WIDTH)**2)**3


def boundary_rows(a):
    dom, chi, intake, path = prepared(a)
    # Preserve the four-volume quadrature on the rounded feed, slab and conduit.
    original = (0, 1, 2, 5)
    grids = [dom.grids[i] for i in original]
    p, v = grids[2]
    cut_index = int(np.argmin(abs(p-a.cut)))
    cut = p[cut_index]
    grids[2] = p[:cut_index+1], v

    class Charts:
        @staticmethod
        def chart(index, p, v):
            return dom.chart(original[index], p, v)

    pts, blocks = mesh_geometry(Charts(), grids, a.sub, False)
    nq = 4
    x, weights = leggauss(nq)
    xi, eta = [a.ravel() for a in np.meshgrid(x, x, indexing='ij')]
    shape = product(linear(xi), linear(eta))
    deriv = np.tile([-.5,.5], (nq*nq,1))
    dp, dv = product(deriv, linear(eta)), product(linear(xi), deriv)
    count = len(pts)
    mass = np.zeros((3,count))
    stiff = [sparse.csr_matrix((count,count)) for _ in mass]
    damping = np.zeros((2,count))
    drives = np.zeros((2,count))
    pressure_probes = np.zeros((2,count))
    fixed = pts[:,1]==-1.
    order_min = 1.
    for block,(pg,vg,ids) in enumerate(blocks):
        which = original[block]
        po, vo, archived = dom.order_blocks[which]
        ii, jj = [t.ravel() for t in np.indices((len(pg)-1,len(vg)-1))]
        hp, hv = np.diff(pg)[ii,None]/2, np.diff(vg)[jj,None]/2
        p, v = pg[ii,None]+hp*(1+xi), vg[jj,None]+hv*(1+eta)
        # Subdivision respects every original chart-cell boundary.
        oi, oj = ii//a.sub, jj//a.sub
        old_ids = archived[oi*(len(vo)-1)+oj]
        def get_order(p,v):
            ax = quadratic(2*(p-po[oi,None])/np.diff(po)[oi,None]-1)[0]
            ay = quadratic(2*(v-vo[oj,None])/np.diff(vo)[oj,None]-1)[0]
            return np.einsum('eqi,ei->eq', tensor(ax,ay), chi[old_ids])

        order = get_order(p,v)
        order_min = min(order_min,float(order.min()))
        if np.min(order)<=0:
            raise ValueError('dual field encountered zero shear stiffness before the receiving cut')
        r,w,rp,rv,wp,wv = dom.chart(which,p,v)
        jac = rp*wv-rv*wp
        gp,gv = dp/hp[:,:,None],dv/hv[:,:,None]
        gr = (wv[:,:,None]*gp-wp[:,:,None]*gv)/jac[:,:,None]
        gz = (rp[:,:,None]*gv-rv[:,:,None]*gp)/jac[:,:,None]
        wt = 4*np.pi*abs(jac)*hp*hv*np.outer(weights,weights).ravel()/order
        dr = 3*shape+r[:,:,None]*gr
        dz = r[:,:,None]*gz
        group = np.where(w<0,0,np.where(order>=.1,1,2))
        for g in np.unique(group):
            use = wt*(group==g)
            np.add.at(mass[g], ids.ravel(), ((use*r**4)@shape).ravel())
            grad = use*r*r
            loc = (np.einsum('eq,eqi,eqj->eij',grad,dr,dr)
                   +np.einsum('eq,eqi,eqj->eij',grad,dz,dz))
            stiff[g] += sparse.coo_matrix((loc.ravel(),(np.repeat(ids,4,axis=1).ravel(),
                               np.tile(ids,(1,4)).ravel())),shape=(count,count)).tocsr()
        edge_ids = ids.reshape(len(pg)-1,len(vg)-1,2,2)
        if which==0:
            fixed[edge_ids[:,0,:,0].ravel()] = True   # freely sliding rounded lip
        if which in (1,5):
            fixed[edge_ids[:,0,:,0].ravel()] = True   # unloaded planar faces
            fixed[edge_ids[:,-1,:,-1].ravel()] = True
        if which in (2,5):
            pe = np.full((len(ii),nq),pg[-1])
            ve = vg[jj,None]+hv*(1+x)
            rr,_,_,rrv,_,wwv = dom.chart(which,pe,ve)
            last = ii==len(pg)-2
            # Evaluate only the element row owning the cut; other rows have
            # different polynomial charts and are irrelevant to the integral.
            olast = oi[last]
            u = quadratic(2*(pe[last]-po[olast,None])/np.diff(po)[olast,None]-1)[0]
            z = quadratic(2*(ve[last]-vo[oj[last],None])/np.diff(vo)[oj[last],None]-1)[0]
            ch = np.einsum('eqi,ei->eq',tensor(u,z),chi[old_ids[last]])
            if ch.min()<=0:
                raise ValueError('zero material impedance on numerical port')
            area = 4*np.pi*rr[last]**4*np.hypot(rrv[last],wwv[last])*hv[last]*weights/ch
            port_nodes = edge_ids[-1,:,-1,:]
            np.add.at(damping[int(which==2)],port_nodes.ravel(), (area@linear(x)).ravel())
        if which==2:
            wall = edge_ids[:,0,:,0]
            h = np.diff(pg)[:,None]/2
            z = pg[:-1,None]+h*(1+x)
            r,_,slope,_,_,_ = dom.chart(which,z,np.zeros_like(z))
            # ds v_t for unit w-directed stroke is h(w) dw. A radial stroke
            # with identical normal speed gives -a'(w)^2 h(w) dw.
            axial = 4*np.pi*r**3*h*weights*profile(z)
            for drive,weight in zip(drives,(axial,-slope**2*axial),strict=True):
                np.add.at(drive,wall.ravel(),(weight@linear(x)).ravel())
            # Signed average of physical shear stress b_t/r² = r H_t.
            for sensor,depth in zip(pressure_probes,(.88,1.12),strict=True):
                patch = np.maximum(0.,1-((z-depth)/.045)**2)**2*h*weights
                np.add.at(sensor,wall.ravel(),((r*patch)@linear(x)).ravel()/patch.sum())
            if a.free_wall:
                nodes=np.unique(wall)
                fixed[nodes[abs(pts[nodes,1]-CENTER)>=WIDTH]]=True
    live = ~fixed
    damp, drives, sensors = damping[:,live],drives[:,live],pressure_probes[:,live]
    mass, stiff = mass[:,live], [part[live][:,live] for part in stiff]
    emit(test='geometry',radius=1.,intake=intake,order_file=str(path),sha256=digest(path),
         cut=cut,outer=a.outer,mesh=a.mesh,subdivision=a.sub,free_outside_patch=a.free_wall,
         minimum_sampled_order=order_min,nodes=int(live.sum()),boundary_drive_sums=drives.sum(axis=1).tolist(),
         region_mass_sums=mass.sum(axis=1).tolist(),port_sums=damp.sum(axis=1).tolist())
    if mass.sum(axis=0).min()<=0 or np.any(damp.sum(axis=1)<=0) or drives[0].sum()<=0 or drives[1].sum()>=0:
        raise ValueError('invalid meridional energy or boundary load')
    return mass,stiff,damp,drives,sensors


def run(a):
    regions,parts,ports,drives,sensors = boundary_rows(a)
    force = {'normal':0.*drives[0],'axial':drives[0],'radial':drives[1]}[a.motion]
    mass,k,damp = regions.sum(axis=0),sum(parts).tocsr(),ports.sum(axis=0)
    # Centered, staggered velocities and exact discrete work balance.
    bound = np.max(np.asarray(abs(k).sum(axis=1)).ravel()/mass)
    stride = int(np.ceil(SAMPLE*np.sqrt(bound)/1.35))*a.time_refine
    dt = SAMPLE/stride
    update = dt/(mass+.5*dt*damp)
    decay = (mass-.5*dt*damp)/(mass+.5*dt*damp)
    pos,vel = np.zeros_like(mass),np.zeros_like(mass)
    ends = [np.flatnonzero(p) for p in ports]
    driven = np.flatnonzero(force)
    work, pos_work, neg_work = 0.,0.,0.
    lost = np.zeros(2)
    snapshots = []
    for i in range(round(a.end/SAMPLE)*stride):
        clock = i*dt
        x = clock/a.duration
        speed = a.amplitude*30*x*x*(1-x)**2/a.duration if 0<x<1 else 0.
        residual = -(k@pos)
        if speed:
            residual[driven] += speed*force[driven]
        new = decay*vel+update*residual
        next_pos = pos+dt*new
        for j,edge in enumerate(ends):
            lost[j] += .25*dt*(ports[j,edge]@(new[edge]+vel[edge])**2)
        if speed:
            added = .5*dt*speed*(force[driven]@(vel[driven]+new[driven]))
            work += added
            pos_work += max(added,0.)
            neg_work += min(added,0.)
        if (i+1)%stride==0:
            energy = .5*(regions@(new**2)+np.array([pos@(mat@next_pos) for mat in parts]))
            snapshots.append([(i+.5)*dt,work,*energy,*lost,energy.sum()+lost.sum()-work,
                              *(sensors@new)])
        vel,pos = new,next_pos
    values=np.array(snapshots)
    names=['time','work','feeding_field','ordered_throat','weak_throat','remote_feed','receiving_cut',
           'balance_error','collar_stress','deep_stress']
    np.savez_compressed(a.outdir/'trace.npz',columns=names,readings=values)
    windows={}
    for lo,hi in ((4,8),(20,24),(44,48)):
        window=values[(values[:,0]>=lo)&(values[:,0]<hi)]
        if len(window):
            means=window.mean(axis=0)
            windows[f'{lo}-{hi}']=dict(zip(names[1:7],means[1:7].tolist(),strict=True))
    emit(test='stroke',motion=a.motion,amplitude=a.amplitude,duration=a.duration,dt=dt,
         signed_wall_work=work,positive_work=pos_work,negative_work=neg_work,peak_stress=np.max(abs(values[:,8:]),axis=0).tolist(),
         windows=windows,last=dict(zip(names,values[-1].tolist(),strict=True)),
         max_balance_error=np.max(abs(values[:,7])),trace_sha256=digest(a.outdir/'trace.npz'))
    if np.max(abs(values[:,7]))>1e-8*max(abs(work),1.e-30) or not np.isfinite(values).all():
        raise ValueError('work, port losses and stored energy fail to balance')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.set_defaults(radius=1.,state=None)
    p.add_argument('--mesh',type=int,default=34)
    p.add_argument('--sub',type=int,default=1)
    p.add_argument('--cut',type=float,default=1.85)
    p.add_argument('--outer',type=float,default=14.)
    p.add_argument('--feed-spacing',type=float,default=.1)
    p.add_argument('--end',type=float,default=48.)
    p.add_argument('--duration',type=float,default=3.)
    p.add_argument('--time-refine',type=int,default=1)
    p.add_argument('--amplitude',type=float,default=.01)
    p.add_argument('--free-wall',action='store_true')
    p.add_argument('--motion',choices=('axial','normal','radial'),required=True)
    p.add_argument('--outdir',type=Path,required=True)
    a=p.parse_args()
    if not (a.end>=a.duration>0 and a.sub>=1 and a.time_refine>=1 and a.outer>9. and 1.2<a.cut<3.05):
        p.error('use a resolved throat, two ports and a complete stroke')
    run(a)
