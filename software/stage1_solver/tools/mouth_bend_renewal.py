"""Hold the two bends while wall contacts renew in place or travel with a tagged flow.

The shear field equilibrates quasi-statically. Its dipole dual uses the whole
feeding and receiving domain. The two rest-slip components are resolved in the
orthonormal wall frame. The angular area factor is restored only for readouts.
"""
import argparse
import json
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy import linalg, sparse
from scipy.sparse.linalg import splu

from mouth_bend_static_shear import ANGULAR, archive_inputs, stiffness_matrix
from mouth_carrier_lateral_inertia import taper
from mouth_curved_phase import quadratic, tensor
from mouth_empty_bend_stiffness import INTERVALS
from mouth_order_wave_pulse import digest, emit, linear, mesh_geometry, product

ROOT = Path('_scratch/s11c')


def gram(weights, left, right):
    return np.einsum('eq,eqi,eqj->eij', weights, left, right, optimize=True)


def volume_cost(geo, order, sub, quad, floor):
    """Complementary curl energy. Free tangential boundaries have m=0."""
    points, blocks = mesh_geometry(geo, geo.grids, sub, False)
    roots, weights = leggauss(quad)
    xi, eta = [x.ravel() for x in np.meshgrid(roots,roots,indexing='ij')]
    shape = product(linear(xi),linear(eta))
    derivative = np.tile([-.5,.5], (len(xi),1))
    sp, sv = product(derivative,linear(eta)), product(linear(xi),derivative)
    samples = np.array([0.,.5,1.])
    bx, by = [x.ravel() for x in np.meshgrid(samples,samples,indexing='ij')]
    transform = np.array([[1,0,0],[-.5,2,-.5],[0,0,1]])
    count=len(points)
    active=np.ones(count,bool)
    matrix=sparse.csc_matrix((count,count))
    all_edges=[]
    for block,((pg,vg,cells),(po,vo,archived)) in enumerate(zip(blocks,geo.order_blocks,strict=True)):
        i,j=[x.ravel() for x in np.indices((len(pg)-1,len(vg)-1))]
        hp,hv=np.diff(pg)[i,None]/2,np.diff(vg)[j,None]/2
        pp, vv=pg[i,None]+hp*(1+xi),vg[j,None]+hv*(1+eta)
        original=order[archived[(i//sub)*(len(vo)-1)+j//sub]]
        def interpolate(x,y):
            px=quadratic(2*((i%sub)[:,None]+x)/sub-1)[0]
            vy=quadratic(2*((j%sub)[:,None]+y)/sub-1)[0]
            return np.einsum('eqk,ek->eq',tensor(px,vy),original)
        coeff=transform@interpolate(bx,by).reshape(-1,3,3)@transform.T
        live=coeff.min(axis=(1,2))>floor
        active[cells[~live].ravel()]=False
        chi=np.clip(interpolate((xi+1)/2,(eta+1)/2),0,1)
        r,_,rp,rv,zp,zv=geo.chart(block,pp,vv)
        det=rp*zv-rv*zp
        dp,dv=sp/hp[:,:,None],sv/hv[:,:,None]
        gr=(zv[:,:,None]*dp-zp[:,:,None]*dv)/det[:,:,None]
        gz=(rp[:,:,None]*dv-rv[:,:,None]*dp)/det[:,:,None]
        weight=abs(det)*hp*hv*np.outer(weights,weights).ravel()
        inv=np.divide(weight,chi,out=np.zeros_like(chi),where=live[:,None])
        values=np.broadcast_to(shape,gr.shape)
        local=(gram(inv,gr,gr)+gram(inv,gz,gz))/2+gram(inv/r**2,values,values)
        matrix+=stiffness_matrix(local,cells,count)
        grid=cells.reshape(len(pg)-1,len(vg)-1,4)
        all_edges.extend((grid[0,:,[0,1]].T,grid[-1,:,[2,3]].T,
                          grid[:,0][:,[0,2]],grid[:,-1][:,[1,3]]))
        if block==2:
            guide_edges=grid[:,0][:,[0,2]]
            zedges=pg
    pairs,counts=np.unique(np.sort(np.vstack(all_edges),axis=1),axis=0,return_counts=True)
    guide_set={tuple(xx) for xx in np.sort(guide_edges,axis=1)}
    for pair in pairs[counts==1]:
        exterior=points[pair]
        far=np.all(exterior[:,0]==geo.far) or np.all(exterior[:,1]==geo.roof)
        if not far and tuple(pair) not in guide_set:
            active[pair]=False
    active[points[:,0]==0]=False
    wall=np.r_[guide_edges[:,0],guide_edges[-1,1]]
    boundary=wall[active[wall]]
    interior=active.copy()
    interior[wall]=False
    coupling=matrix[interior][:,boundary]
    inner=matrix[interior][:,interior]
    scale=sparse.diags(1/np.sqrt(inner.diagonal()))
    solver=splu((scale@inner@scale).tocsc())
    cross=scale@coupling
    schur=matrix[boundary][:,boundary].toarray()-cross.T@solver.solve(cross.toarray())
    return zedges, active[wall], (schur+schur.T)/2


def wall_data(geo, order, field, edges, keep, quad):
    root,weights=leggauss(quad)
    t=(1+root)/2
    width=np.diff(edges)
    z=edges[:-1,None]+width[:,None]*t
    parent,_,nodes=geo.order_blocks[2]
    idx=np.searchsorted(parent,z,side='right')-1
    cells=nodes[idx*(len(geo.order_blocks[2][1])-1)]
    def wall_field(data,depth):
        q,dq=quadratic(2*(depth-parent[idx])/np.diff(parent)[idx]-1)
        wall=data[cells[:,:,[0,3,6]]]
        return (q*wall).sum(axis=-1),(dq*wall).sum(axis=-1)*2/np.diff(parent)[idx]
    def speed(depth):
        psi,dpsi=wall_field(field,depth)
        return (dpsi/psi).imag/(1+geo.radius(depth)[1]**2)
    velocity=speed(z)
    if velocity.min()<=0:
        raise ValueError('Non-inward carrier speed on the guide')
    chi=np.clip(wall_field(order,z)[0],0,1)
    radius,slope=geo.radius(z)
    arc=np.sqrt(1+slope**2)
    def age_at(w):
        span=np.clip(w-edges[:-1],0,width)
        # Each edge is inside an original element; collapsed terminal cells
        # use that same element's one-sided velocity with zero integration weight.
        zp=edges[:-1,None]+span[:,None]*t
        return float(np.sum(span[:,None]*weights/2/speed(zp)))
    clock=age_at(2.)
    last=np.flatnonzero(keep)[-1]+1
    length=width[:last]
    chi,radius,arc,velocity,z=[x[:last] for x in (chi,radius,arc,velocity,z)]
    n=len(length)
    ids=np.arange(n)[:,None]+np.array([0,1])
    wall_keep=keep[:n+1]
    line=length[:,None]*weights/2
    val=np.broadcast_to(np.stack((1-t,t),axis=-1),(n,len(t),2))
    deriv=np.broadcast_to(np.array([-1.,1.]),val.shape)/length[:,None,None]
    used=wall_keep[ids].any(axis=1)
    if chi[used].min()<=0:
        raise ValueError('Active dual trace encountered zero order')
    inv=np.divide(line,chi,out=np.zeros_like(line),where=used[:,None])
    local=gram(inv*arc/radius**2,val,val)+gram(inv/(2*arc),deriv,deriv)
    contact=stiffness_matrix(local,ids,n+1)[wall_keep][:,wall_keep].toarray()
    load=np.stack([radius*taper(z,*band)[1] for band in INTERVALS],axis=-1)
    rhs=np.zeros((n+1,2))
    np.add.at(rhs,ids,gram(line,val,load))
    # DG1 rest-slip coefficients, discontinuous between cells. B integrates
    # physical slips against m and m_z; G returns local extension = traction/chi.
    time_mass=gram(line/velocity,val,val)
    rest_nodes=np.arange(2*n).reshape(n,2)
    def scatter_rect(elem):
        return sparse.coo_matrix((elem.ravel(),
                    (np.repeat(ids,2,axis=1).ravel(),np.tile(rest_nodes,(1,2)).ravel())),
                    shape=(n+1,2*n)).toarray()[wall_keep]
    bz=scatter_rect(gram(line*arc,val,val))
    ba=scatter_rect(gram(line*radius,deriv,val))
    gz=gram(inv/(velocity*radius**2),val,val)
    ga=gram(inv/(2*velocity*radius*arc),val,deriv)
    def scatter_rate(elem):
        solved=np.linalg.solve(time_mass,elem)
        return scatter_rect(solved.transpose(0,2,1)).T
    G=np.vstack((scatter_rate(gz),scatter_rate(ga)))
    B=np.hstack((bz,ba))
    a_cell=-gram(line,deriv,val)
    a_cell[:,1,1]+=1
    advection=stiffness_matrix(a_cell,rest_nodes,2*n).tolil()
    for j in range(1,n):
        advection[2*j,2*j-1]=-1
    mass_inverse=sparse.block_diag(list(np.linalg.inv(time_mass)),format='csr')
    upwind=(mass_inverse@advection.tocsr()).toarray()
    return (contact,rhs[wall_keep],B,G,linalg.block_diag(upwind,upwind),clock,
            {'projected_clock_to_deep_end':clock,'to_exit':age_at(edges[-1]),
             'min_sampled_wall_speed':float(velocity.min()),'active_wall_end':float(edges[n])})


def response(force):
    stiffness=force+force.T
    return np.linalg.eigvalsh(stiffness/2).tolist()


def run(args):
    label=f'delayed_n{args.mesh}'
    order_path=ROOT/f'mouth_free_order_{label}_d04/order_seed1.55.npz'
    flow=ROOT/f'bend_static_shear_delcarrier{args.mesh}_01'
    geo,order,carrier=archive_inputs(order_path,flow)
    psi_path=ROOT/f'mouth_flared_{label}/field.npz'
    source=json.loads((order_path.parent/'stdout').read_text().splitlines()[0])
    if digest(psi_path)!=source['source_sha256']:
        raise ValueError('Changed carrier preparation')
    with np.load(psi_path) as data:
        np.testing.assert_array_equal(geo.points,data['points'])
        field=data['psi']
    edges,keep,bulk=volume_cost(geo,order,args.sub,args.quad,args.floor)
    contact,load,B,G,drift,clock,settings=wall_data(geo,order,field,edges,keep,args.quad)
    matrix=bulk+contact
    # Normalize the complementary elastic solve before applying its response.
    scale=1/np.sqrt(matrix.diagonal())
    factor=linalg.cho_factor(scale[:,None]*matrix*scale)
    def solve(rhs):
        return scale[:,None]*linalg.cho_solve(factor,scale[:,None]*rhs)
    initial=solve(load)
    inverseB=solve(B)
    read=ANGULAR*load.T
    frozen=read@initial
    eig,V=linalg.eigh(scale[:,None]*contact*scale,scale[:,None]*matrix*scale)
    if eig.min() < -1e-9 or eig.max()>1+1e-8:
        raise ArithmeticError('Invalid local renewal spectrum')
    mode=V.T@(scale[:,None]*load)
    def site(t):
        return ANGULAR*mode.T@(np.exp(-eig*t)[:,None]*mode)
    emit(test='initial',**settings,sub=args.sub,quad=args.quad,positive_order_floor=args.floor,
         carrier=carrier.tolist(),frozen=frozen.tolist(),net_eigs=response(carrier+frozen),
         site_spectrum=[eig.min(),eig.max()],transport_speed_over_carrier=args.transport_speed)
    samples=(0.,.25,.5,1.,2.,4.,8.,16.)
    for t in samples:
        s=site(t)
        emit(test='fixed_sites',time_over_lifetime=t,shear=s.tolist(),net_eigs=response(carrier+s))
    from scipy.optimize import brentq
    root=brentq(lambda t:response(carrier+site(t))[0],0,1e4)
    filtered=solve(load)
    direct=linalg.solve(bulk+2*contact,load,assume_a='pos')
    laplace_error=np.max(abs(read@direct-ANGULAR*mode.T@(mode/(1+eig[:,None]))))
    emit(test='site_check',first_zero_over_lifetime=root,filter_error=laplace_error)
    if laplace_error>1e-8 or np.max(abs(site(0)-read@filtered))>1e-8:
        raise ArithmeticError('Local renewal elastic correspondence failed')
    # Steady streaming admits direct integration of the same DG transport
    # from the inflow. The wall complementary solve is smaller than the
    # history system, allowing a separately imposed inlet traction condition.
    reference_drift=args.transport_speed*clock*drift
    history=linalg.solve(reference_drift,G)
    advected_cost=B@history
    constant=np.zeros((len(drift),2))
    constant[:len(drift)//2,0]=1
    constant[len(drift)//2:,1]=1
    if not keep[1] or keep[0]:
        raise ValueError('Expected a free joining lip and first active guide node')
    def saturated(tau, stress_free):
        norm=scale[:,None]*(matrix+advected_cost/tau)*scale
        m=scale[:,None]*linalg.solve(norm,scale[:,None]*load)
        entry=np.zeros((2,2))
        if stress_free:
            adjust=scale*linalg.solve(norm,scale*(B@constant[:,1]))
            # Anchor the meridional inlet reference to the straight feed.
            # m=0 at the lip sets its traction to zero. Zero angular
            # traction sets the one-sided derivative of m to zero as well.
            entry[1]=m[0]/adjust[0]
            m-=adjust[:,None]*entry[1]
        return m,entry
    def minimum(tau, stress_free):
        m,_=saturated(tau,stress_free)
        return response(carrier+read@m)[0]
    for free in (False,True):
        left,right=.03,.25
        if minimum(left,free)*minimum(right,free)<0:
            tc=brentq(lambda tau:minimum(tau,free),left,right)
            emit(test='dc_crossing',stress_free_entry=free,bracket=[left,right],
                 lifetime_over_clock=tc,model_time=tc*clock)
    for tau in args.lifetimes:
        # t and tau in units of the independent wall transit clock. Both
        # meridional and angular amplitudes are parallel transported along a meridian.
        feedback=G@inverseB/tau
        drift_op=reference_drift
        operator=drift_op+feedback
        source_term=G@initial/tau
        dc=linalg.solve(operator,source_term)
        mdc=initial-inverseB@dc
        residual=drift_op@dc-G@mdc/tau
        diss=(B.T@mdc).T@(G@mdc)/tau
        balance=(B.T@mdc).T@residual
        static=read@mdc
        emit(test='travelling_dc',lifetime_over_clock=tau,shear=static.tolist(),net_eigs=response(carrier+static),
             renewal_power_per_clock=(ANGULAR*diss).tolist(),power_balance_error=ANGULAR*np.max(abs(balance)),
             relative_dc_residual=np.linalg.norm(residual)/max(1.,np.linalg.norm(source_term)))
        for free in (False,True):
            m,entry=saturated(tau,free)
            if not free and np.max(abs(read@(m-mdc)))>1e-8:
                raise ArithmeticError('Direct transit integration disagrees with DG steady solve')
            material=read@m
            total=material+carrier
            exact_power=ANGULAR*m.T@contact@m/tau
            transported=history@m/tau+constant@entry
            renewal_term=ANGULAR*(B.T@m).T@(G@m)/tau
            stream_term=ANGULAR*(B.T@m).T@(drift_op@(transported-constant@entry))
            emit(test='inlet_comparison',stress_free_entry=free,lifetime_over_clock=tau,
                 incoming_reference=entry.tolist(),net_eigs=response(total),net=total.tolist(),
                 entrance_test_deep_holder=float(.01*total[1,0]),
                 deep_test_entrance_holder=float(.01*total[0,1]),
                 ccw_holder_cycle_work=float((total[1,0]-total[0,1])*.01**2),
                 positive_reference_loss_per_clock=exact_power.tolist(),
                 relative_loss_projection_error=float(np.linalg.norm(renewal_term-exact_power)/np.linalg.norm(exact_power)),
                 relative_stream_balance_error=float(np.linalg.norm(stream_term-renewal_term)/np.linalg.norm(renewal_term)))
        dt=1/args.steps
        ident=np.eye(len(dc))
        first=linalg.lu_factor(ident/dt+operator)
        bdf=linalg.lu_factor(1.5*ident/dt+operator)
        ref=np.zeros_like(dc); prev=ref.copy()
        sample_steps={round(t*args.steps) for t in samples}
        for step in range(1,16*args.steps+1):
            if step==1:
                updated=linalg.lu_solve(first,source_term)
            else:
                updated=linalg.lu_solve(bdf,source_term+(2*ref-.5*prev)/dt)
            prev,ref=ref,updated
            if step in sample_steps:
                stiffness=read@(initial-inverseB@ref)
                emit(test='travelling_hold',time_over_clock=step*dt,lifetime_over_clock=tau,
                     shear=stiffness.tolist(),net_eigs=response(carrier+stiffness),
                     ideal_feed_marker_reaches_deep_end=bool(step*dt>=1))
        if np.linalg.norm(residual)>1e-7*max(1.,np.linalg.norm(source_term)):
            raise ArithmeticError('Advective stationary residual failed')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--mesh',type=int,choices=(34,48),default=48)
    p.add_argument('--sub',type=int,default=1)
    p.add_argument('--quad',type=int,default=4)
    p.add_argument('--floor',type=float,default=1e-8)
    p.add_argument('--steps',type=int,default=32)
    p.add_argument('--lifetimes',type=float,nargs='+',default=[.25,1,4])
    p.add_argument('--transport-speed',type=float,default=1.,help='Reference transport speed relative to the projected carrier')
    args=p.parse_args()
    if args.steps<=0 or args.transport_speed<=0 or min(args.lifetimes)<=0:
        p.error('Steps, traveling reference speed and renewal times must be positive')
    run(args)
