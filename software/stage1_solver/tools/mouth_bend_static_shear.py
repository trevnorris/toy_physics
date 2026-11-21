"""Two variational estimates of static curl elasticity for an attached draining mouth.

On the conduit r=a(z)s, subtract the gradient of the angular displacement
potential from a dipole displacement one-form. P ds+T dz is its remaining
meridional covector, multiplied by the direction cosine Y. The conjugate scalar
m is the angular dipole of the radial/depth curl multiplied by chi*r^2.
Both trial fields extend into the reservoirs with zero curl or zero m.
"""
import argparse
import json
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy import sparse
from scipy.linalg import eigvalsh
from scipy.sparse.linalg import splu

from mouth_carrier_lateral_inertia import taper
from mouth_curved_phase import quadratic, tensor
from mouth_empty_bend_stiffness import INTERVALS
from mouth_order_wave_pulse import ARCHIVE, Geometry, digest, emit

FLOW = Path('research/pde_ledger_v3/_measurements/mouth_empty_bend_stiffness.json')
ANGULAR = 4*np.pi/3


def archive_inputs(order_path, flow_run):
    record = json.loads(ARCHIVE.read_text())
    state = json.loads((order_path.parent/'stdout').read_text().splitlines()[0])
    if digest(order_path) != record['diffuse_orders_sha256'][str(order_path)]:
        raise ValueError('Order hash disagrees with the reference record')
    if json.loads((order_path.parent/'outcome.json').read_text())['exitCode']:
        raise ValueError('Order source did not complete')
    if flow_run is None:
        prior = json.loads(FLOW.read_text())
        src = prior['runs'][f"_scratch/s11c/empty_bend_c{state['n']}_01"]
        if not src['verified_limits'] or src['worker_exit'] or src['guarded_exit']:
            raise ValueError('Unverified reference carrier')
        lines = src['literal_stdout']
        source_record = FLOW
    else:
        outcome=json.loads((flow_run/'outcome.json').read_text())
        limits=json.loads((flow_run/'limit-validation.json').read_text())
        if (outcome['exitCode'] or outcome['childOutcome']['exitCode'] or not outcome['limitsVerified']
                or not limits['verified'] or int(limits['actual']['memory.swap.max'])):
            raise ValueError('Unverified flow calculation')
        source_record=flow_run/'stdout'
        lines=source_record.read_text().splitlines()
    meta,carrier=[json.loads(line) for line in lines]
    if (state['source_sha256']!=meta['sha256'] or state['n']!=meta['mesh'] or
            state['profile']!=meta['profile'] or carrier['extension']!='compact'):
        raise ValueError('Unmatched carrier or bend geometry')
    geo = Geometry(meta['mesh'],meta['profile'],meta['roof'],meta['far'])
    with np.load(order_path) as data:
        order = data['chi']
    emit(test='preparation',mesh=state['n'],rate=meta['rate'],radius=1.,order=str(order_path),
         sha256=digest(order_path),flow_record_sha256=digest(source_record))
    return geo, order, np.array(carrier['holding_stiffness'])


def stiffness_matrix(local, nodes, count):
    n = nodes.shape[1]
    return sparse.coo_matrix((local.ravel(), (np.repeat(nodes,n,axis=1).ravel(),
                                              np.tile(nodes,(1,n)).ravel())),shape=(count,count)).tocsc()


def stationary_solve(matrix, rhs, active):
    dofs=np.flatnonzero(active)
    sub=matrix[dofs][:,dofs]
    scale=sparse.diags(1/np.sqrt(sub.diagonal()))
    answer=np.zeros_like(rhs)
    answer[dofs]=scale@splu((scale@sub@scale).tocsc()).solve(scale@rhs[dofs])
    return answer


def subdivide(edges,n):
    return np.r_[(edges[:-1,None]+np.diff(edges)[:,None]*np.arange(n)/n).ravel(),edges[-1]]


def calculate(geo,order,sub,quad,deep_cutoff,contacts):
    zz, vv, parent_nodes = geo.order_blocks[2]
    zedge, sedge = subdivide(zz,sub),subdivide(1-vv[::-1],sub)
    nz, ns = len(zedge)-1,len(sedge)-1
    iz, js = [v.ravel() for v in np.indices((nz,ns))]
    shape=(nz+1,ns+1)
    count=(nz+1)*(ns+1)
    nodes=(iz*(ns+1)+js)[:,None]+np.array([0,1,ns+1,ns+2])
    vector_nodes=np.column_stack((nodes,nodes+count))
    dz,ds=np.diff(zedge)[iz,None],np.diff(sedge)[js,None]
    base_order=order[parent_nodes[(iz//sub)*(len(vv)-1)+len(vv)-2-js//sub]]
    root,weights=leggauss(quad)
    unit=(1+root)/2
    zzq,ssq=[x.ravel() for x in np.meshgrid(unit,unit,indexing='ij')]
    wz,ws=[x.ravel() for x in np.meshgrid(weights/2,weights/2,indexing='ij')]
    z=zedge[iz,None]+dz*zzq
    s=sedge[js,None]+ds*ssq
    a,ap=geo.radius(z)
    r=a*s
    measure=a*dz*ds*wz*ws
    def chi_at(t,u):
        # The archive's radial chart runs from wall to axis.
        qz=quadratic(2*((iz%sub)[:,None]+t)/sub-1)[0]
        qv=quadratic(1-2*((js%sub)[:,None]+u)/sub)[0]
        return np.einsum('eqk,ek->eq',tensor(qz,qv),base_order)
    chi=np.clip(chi_at(zzq,ssq),0,1)
    # A tensor quadratic is bounded by its Bernstein coefficients. Retain
    # dual basis supports only on cells with a strictly positive lower bound.
    grid=np.array([0.,.5,1.])
    tz,ts=[x.ravel() for x in np.meshgrid(grid,grid,indexing='ij')]
    transform=np.array([[1,0,0],[-.5,2,-.5],[0,0,1]])
    bern=transform@chi_at(tz,ts).reshape(-1,3,3)@transform.T
    dual_cells=bern.min(axis=(1,2))>deep_cutoff

    lin=lambda x:np.stack((1-x,x),axis=-1)
    deriv=np.tile([-1.,1.],(len(ssq),1))
    def bilinear(x,y):
        return (x[..., :,None]*y[...,None,:]).reshape(-1,4)
    val=bilinear(lin(zzq),lin(ssq))
    vz=bilinear(deriv,lin(ssq))/dz[:,:,None]
    vs=bilinear(lin(zzq),deriv)/ds[:,:,None]
    blank=np.zeros_like(vs)
    def rows(first,second):
        return np.concatenate((np.broadcast_to(first,vs.shape),np.broadcast_to(second,vs.shape)),axis=2)
    curl=rows(-vz,vs)/a[:,:,None]
    radial=rows(val,blank)/a[:,:,None]
    axial=rows(-ap[:,:,None]*s[:,:,None]/a[:,:,None]*val,val)
    slopes=np.stack([-taper(z,*interval)[1] for interval in INTERVALS],axis=-1)
    # Exact tangential boundary data from the physical section translation:
    # T(1,z)=-a*h'. The lift P=0, T=-r*h' has curl=-h'.
    curl_lift=-slopes
    axial_lift=-r[:,:,None]*slopes
    wt=measure*chi
    gram=lambda m,x,y:np.einsum('eq,eqi,eqj->eij',m,x,y,optimize=True)
    # A tiny positive floor supplies a legal primal trial for the softer
    # exact material. Evaluate and report its energy with the actual chi.
    solve_wt=measure*np.maximum(chi,1e-11)
    def blocks(w):
        return (gram(w*r*r,curl,curl)+2*(gram(w,radial,radial)+gram(w,axial,axial)),
                gram(w*r*r,curl,curl_lift)+2*gram(w,axial,axial_lift))
    kl, force=blocks(solve_wt)
    rhs=np.zeros((2*count,2))
    np.add.at(rhs,vector_nodes,-force)
    # On both ends of the conduit the tangential radial covector P vanishes.
    # T can jump normally there. On the wall T's correction vanishes.
    free_p=np.ones(shape,bool); free_p[[0,-1],:]=False; free_p[:,0]=False
    free_t=np.ones(shape,bool); free_t[:,[0,-1]]=False
    v=stationary_solve(stiffness_matrix(kl,vector_nodes,2*count),rhs,np.r_[free_p.ravel(),free_t.ravel()])
    loc=v[vector_nodes]
    cs=curl@loc+curl_lift
    ar=radial@loc
    az=axial@loc+axial_lift
    summ=lambda w,x:gram(w,x,x).sum(axis=0)*ANGULAR
    upper=summ(wt*r*r,cs)+2*(summ(wt,ar)+summ(wt,az))
    locked=summ(wt*r*r,curl_lift)+2*summ(wt,axial_lift)
    # Conforming complementary field. m=0 at both axial ends and on the axis.
    active=np.ones(shape,bool)
    active[[0,-1],:]=False; active[:,0]=False
    active.ravel()[np.unique(nodes[~dual_cells])]=False
    dr=vs/a[:,:,None]
    dw=vz-ap[:,:,None]*s[:,:,None]*dr
    inv=np.zeros_like(wt)
    inv[dual_cells]=measure[dual_cells]/chi[dual_cells]
    dual_local=(gram(inv,dr,dr)+gram(inv,dw,dw))/2+gram(inv/r**2,np.broadcast_to(val,vs.shape),np.broadcast_to(val,vs.shape))
    side=sedge[-1]
    assert side==1
    wzq=zedge[:-1,None]+np.diff(zedge)[:,None]*unit
    radius,_=geo.radius(wzq)
    load=np.stack([radius*taper(wzq,*it)[1] for it in INTERVALS],axis=-1)
    line=np.einsum('eq,qi,eqj->eij',np.diff(zedge)[:,None]*weights/2,lin(unit),load)
    boundary=np.arange(nz)[:,None]*(ns+1)+ns+np.array([0,ns+1])
    drhs=np.zeros((count,2)); np.add.at(drhs,boundary,line)
    dual_matrix=stiffness_matrix(dual_local,nodes,count)
    m=stationary_solve(dual_matrix,drhs,active.ravel())
    lower=ANGULAR*(drhs.T@m+m.T@drhs-m.T@(dual_matrix@m))
    contact_results=[]
    if contacts:
        # Finite tangential spring on ordered material, bulk modulus set to one.
        # Dual m vanishes at both ends. Integration by parts pairs its wall
        # derivative with b-a*h, the angular slip of the material.
        wall=(js==ns-1)
        chi_wall=np.clip(chi_at(unit,np.ones_like(unit))[wall],0,1)
        radius,slope=geo.radius(wzq)
        arc=np.sqrt(1+slope**2)
        line_weight=np.diff(zedge)[:,None]*weights/2
        pos=dual_cells[wall]
        inverse=np.zeros_like(chi_wall)
        inverse[pos]=line_weight[pos]/chi_wall[pos]
        shape1=np.broadcast_to(lin(unit),(*wzq.shape,2))
        derivative=np.broadcast_to(np.array([-1.,1.]),shape1.shape)/np.diff(zedge)[:,None,None]
        local_contact=gram(inverse*arc/radius**2,shape1,shape1)+gram(inverse/(2*arc),derivative,derivative)
        contact=stiffness_matrix(local_contact,boundary,count)
        # An admissible gradient follows the normal wall motion and angular
        # translation. Only its tangential depth component stretches the bond.
        trial=ANGULAR*gram(line_weight*chi_wall*radius**2/arc,load,load).sum(axis=0)
        for strength in contacts:
            cost=dual_matrix+contact/strength
            response=stationary_solve(cost,drhs,active.ravel())
            bound=ANGULAR*(drhs.T@response+response.T@drhs-response.T@(cost@response))
            contact_results.append((strength,bound,strength*trial))
    return lower,upper,locked,contact_results


def run(args):
    dom,order,k=archive_inputs(args.order,args.flow_run)
    low,high,rigid,contacts=calculate(dom,order,args.sub,args.quad,args.cutoff,args.contacts)
    def marginal(shear):
        # For a singular shear estimate, a finite sufficient threshold may be unavailable.
        if np.linalg.eigvalsh(shear)[0]<=0:
            return None
        return eigvalsh(-k,shear)[-1]
    emit(test='static_shear',sub=args.sub,quadrature=args.quad,dual_order_floor=args.cutoff,
         lower=low.tolist(),upper=high.tolist(),rigid_sections=rigid.tolist(),carrier=k.tolist(),
         necessary_mu=marginal(high),sufficient_mu=marginal(low),
         rigid_necessary_mu=marginal(rigid),
         unit_modulus_total_lower=np.linalg.eigvalsh(k+low).tolist(),
         unit_modulus_total_upper=np.linalg.eigvalsh(k+high).tolist(),
         slip_total=np.linalg.eigvalsh(k).tolist())
    for strength,bound,trial in contacts:
        low_eigs=np.linalg.eigvalsh(k+bound)
        slip_eigs=np.linalg.eigvalsh(k+trial)
        emit(test='finite_contact',ordered_modulus=1.,contact_stiffness=strength,
             lower=bound.tolist(),gradient_upper=trial.tolist(),total_lower=low_eigs.tolist(),
             total_gradient_upper=slip_eigs.tolist())
        if np.linalg.eigvalsh(low-bound)[0]<-1e-8 or np.linalg.eigvalsh(trial-bound)[0]<-1e-8:
            raise ArithmeticError('Finite spring variational order failed')
    if np.linalg.eigvalsh(high-low)[0]<-1e-8 or np.linalg.eigvalsh(rigid-high)[0]<-1e-8:
        raise ArithmeticError('Variational ordering failed')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--order',type=Path,required=True)
    p.add_argument('--flow-run',type=Path,help='Guarded carrier bend result for this geometry')
    p.add_argument('--sub',type=int,default=1)
    p.add_argument('--quad',type=int,default=5)
    p.add_argument('--contacts',type=float,nargs='+',default=[],help='Ordered contact stiffnesses at unit bulk modulus')
    p.add_argument('--cutoff',type=float,default=1e-10)
    run(p.parse_args())
