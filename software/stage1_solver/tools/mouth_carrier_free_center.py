"""Solve for the unassigned rigid-center velocity and carrier field together.

Use laboratory carrier phase on a translated chart (no phase boost). The
moving-wall Neumann flux and advection of the background are a weak source.
The last equation sets the total carrier force plus applied push to zero.
"""
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy import sparse
from scipy.sparse.linalg import splu

from mouth_carrier_lateral_response import Dipole, emit, parser, prepare, weighted_gram
from mouth_curved_phase import quadratic, tensor


def basis(dom,block,p,v,i,j):
    pg,vg,dofs=dom.blocks[block]
    hp,hv=np.diff(pg)[i]/2,np.diff(vg)[j]/2
    bp,dp=quadratic((p-pg[i])/hp-1)
    bv,dv=quadratic((v-vg[j])/hv-1)
    n=tensor(bp,bv)
    gp=tensor(dp,bv)/hp[...,None]
    gv=tensor(bp,dv)/hv[...,None]
    r,w,rp,rv,wp,wv=dom.chart(block,p,v)
    det=rp*wv-rv*wp
    gr=(wv[...,None]*gp-wp[...,None]*gv)/det[...,None]
    gw=(rp[...,None]*gv-rv[...,None]*gp)/det[...,None]
    return dofs[i*(len(vg)-1)+j],r,wp,abs(det),n,gr,gw


def motion_source(dom,field):
    # Integrating the translated Neumann work by parts leaves the
    # antisymmetrized volume derivative, including both angular directions.
    vector=np.zeros(len(field),complex)
    root,weights=leggauss(4)
    u,v=[x.ravel() for x in np.meshgrid(root,root,indexing='ij')]
    for block,(pg,vg,dofs) in dom.blocks.items():
        for lo in range(0,len(dofs),256):
            i,j=np.unravel_index(np.arange(lo,min(lo+256,len(dofs))),(len(pg)-1,len(vg)-1))
            i,j=i[:,None],j[:,None]
            hp,hv=np.diff(pg)[i]/2,np.diff(vg)[j]/2
            p,s=pg[i]+hp*(1+u),vg[j]+hv*(1+v)
            idx,r,_,det,n,gr,_=basis(dom,block,p,s,i,j)
            nodes=field[idx]
            f=(n*nodes).sum(axis=-1)
            fr=(gr*nodes).sum(axis=-1)
            wt=4*np.pi*r**2*hp*hv*det*np.outer(weights,weights).ravel()
            contribution=.5j*wt[...,None]*((gr+2*n/r[...,None])*f[...,None]-n*fr[...,None])
            np.add.at(vector,idx[:,0,:].ravel(),contribution.sum(axis=1).ravel())
    return np.r_[vector.real,vector.imag]


def wall_row(dom,field):
    """Historical direct pressure trace at the curved reflecting wall; diagnostic only."""
    static=np.zeros(len(field),complex)
    time=np.zeros_like(static)
    boost=0.
    root,weights=leggauss(8)
    for block in (0,2):
        pg,_,_=dom.blocks[block]
        i=np.arange(len(pg)-1)[:,None]
        hp=np.diff(pg)[:,None]/2
        p=pg[i]+hp*(1+root)
        idx,r,wp,_,n,gr,gw=basis(dom,block,p,np.zeros_like(p),i,np.zeros_like(i))
        nodes=field[idx]
        f=(n*nodes).sum(axis=-1)
        fr=(gr*nodes).sum(axis=-1)
        fw=(gw*nodes).sum(axis=-1)
        wt=-(4*np.pi/3)*r**2*wp*hp*weights
        coefficients=fr[...,None]*gr+fw[...,None]*gw+((abs(f)**8-1)*f/2)[...,None]*n
        np.add.at(static,idx[:,0,:].ravel(),(wt[...,None]*coefficients).sum(axis=1).ravel())
        np.add.at(time,idx[:,0,:].ravel(),(wt[...,None]*1j*f[...,None]*n).sum(axis=1).ravel())
        boost-=np.sum(wt*(f.conj()*fr).imag)
    def real_pair(x):
        return np.r_[x.real,x.imag]
    return real_pair(static),real_pair(time),boost



def momentum_force_rows(dom, field, source):
    """External force on a massless rigid mouth from lab momentum conservation.

    For dipole field h on the translating chart the force is (A+iωB)h.
    A is the stationary momentum flux at the two remote reservoirs. B includes
    lab momentum stored in the volume and the temporal part of reservoir stress.
    The explicit base-current force from the lab time derivative cancels the
    swept base momentum on each moving reservoir face.
    """
    # The same quadrature used for the translation source differentiates
    # P_x = integral Im(psi* d_x psi). The angular dipole average is 1/3.
    static = np.zeros(len(field), complex)
    temporal_port = np.zeros_like(static)
    root, weights = leggauss(8)
    faces = ((1, 'v', dom.far, 0), (4, 'p', 1., 0),
             (3, 'v', dom.roof, 1), (4, 'v', dom.roof, 1))
    for block, along, fixed, outward_axis in faces:
        pg, vg, _ = dom.blocks[block]
        edges = pg if along == 'p' else vg
        half = np.diff(edges)[:, None] / 2
        q = edges[:-1, None] + half * (root + 1)
        p, v = (q, np.full_like(q, fixed)) if along == 'p' else (np.full_like(q, fixed), q)
        i = np.clip(np.searchsorted(pg, p, side='right') - 1, 0, len(pg) - 2)
        j = np.clip(np.searchsorted(vg, v, side='right') - 1, 0, len(vg) - 2)
        idx, r, _, _, n, gr, gw = basis(dom, block, p, v, i, j)
        _, _, rp, rv, wp, wv = dom.chart(block, p, v)
        tr, tw = (rp, wp) if along == 'p' else (rv, wv)
        nr, nw = tw, -tr
        orientation = nr if outward_axis == 0 else nw
        if np.any(np.abs(orientation) < 1e-12):
            raise ValueError('Ambiguous outward reservoir normal')
        nr, nw = np.sign(orientation)[..., None] * nr[..., None], np.sign(orientation)[..., None] * nw[..., None]
        nodal = field[idx]
        f = np.sum(n * nodal, axis=-1, keepdims=True)
        fr = np.sum(gr * nodal, axis=-1, keepdims=True)
        fw = np.sum(gw * nodal, axis=-1, keepdims=True)
        normal_gradient = nr * fr + nw * fw
        grad_stress = fr * (nr * gr + nw * gw) + normal_gradient * (gr + 2*n/r[..., None])
        energy = fr*gr + fw*gw + .5*(abs(f)**8-1)*f*n
        wt = (4*np.pi/3 * r**2 * half * weights)[..., None]
        np.add.at(static, idx.ravel(), (wt * (grad_stress - nr*energy)).ravel())
        np.add.at(temporal_port, idx.ravel(), (-wt * nr * 1j * f * n).ravel())
    return (np.r_[static.real, static.imag],
            (2./3)*source + np.r_[temporal_port.real, temporal_port.imag])


def response(dom,field,model,omega):
    m=weighted_gram(dom.N,dom.weight)
    skew=sparse.bmat([[None,-m],[m,None]],format='csr')
    op=model.operator.astype(complex)-1j*omega*skew
    rigid=motion_source(dom,field)
    static,timerow=momentum_force_rows(dom,field,rigid)
    # Last row uses the force ON the mouth, opposite to the actuator.
    force=-(static+1j*omega*timerow)
    trace_static,trace_time,trace_boost=wall_row(dom,field)
    mask=model.free
    # Independent far-brane phase drive in lab coordinates. The far bulk rests.
    drive=1j*dom.points[:,0]*field*dom.feed
    boundary=np.r_[drive.real,drive.imag]
    rhs=np.column_stack([np.r_[-(op@boundary)[mask],-force@boundary],
                         np.r_[np.zeros(mask.sum()),-1.]])
    bordered=sparse.bmat([[op[mask][:,mask],sparse.csr_matrix(-rigid[mask,None])],
                          [sparse.csr_matrix(force[None,mask]),None]],format='csc')
    fac=splu(bordered)
    result=fac.solve(rhs)
    residual=np.max(abs(bordered@result-rhs),axis=0)
    # Mixed counterdrive checks superposition using the same two source columns.
    ratio=.002
    sign=-.5
    mixed=fac.solve(sign*rhs[:,0]+ratio*rhs[:,1])
    if max(abs(mixed-sign*result[:,0]-ratio*result[:,1]))>2e-8:
        raise ArithmeticError('Mixed-drive superposition failed')
    for col,(label,incident,push) in enumerate([('incident',1.,0.),('push',0.,1.)]):
        z=(incident*boundary).astype(complex)
        z[mask]=result[:-1,col]
        velocity=result[-1,col]
        total=force@z+push
        trace_phase=(trace_static+1j*omega*trace_time)@z
        direct_load=trace_phase+trace_boost*velocity
        acceleration=1j*omega*velocity
        emit(test='free_response',drive=label,omega=omega,velocity=[velocity.real,velocity.imag],
             acceleration=[acceleration.real,acceleration.imag],force_balance=abs(total),
             matrix_residual=residual[col],rim_force_from_phase=abs(trace_phase),
             uniform_moving_wall_term=trace_boost,force_method='volume_momentum_and_remote_stress',
             direct_trace_force=[direct_load.real,direct_load.imag],
             direct_trace_balance=abs(direct_load+push))
        if max(abs(total),residual[col])>2e-8:
            raise ArithmeticError('Coupled force-free center or field failed')


if __name__=='__main__':
    p=parser()
    p.add_argument('--omega',type=float,nargs='+',required=True)
    args=p.parse_args()
    dom,field=prepare(args)
    model=Dipole(dom,field)
    for omega in args.omega:
        response(dom,field,model,omega)
