"""Sideways relative translation of two parts of the radius-one native shear guide.

The wave is azimuthal about the direction of translation. The section at w is
translated by X_entrance + (X_flare-X_entrance)*h(w), together with its order.
The straight wave has degree 4. Compute the complete second variation within
this axisymmetric polarization, allowing the wave to change its spatial shape.
"""
import argparse
from pathlib import Path
from types import SimpleNamespace

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy import sparse
from scipy.sparse.linalg import eigsh, splu

from mouth_carrier_lateral_inertia import taper
from mouth_shear_pattern_release import (angular_stiffness, assemble, digest, emit, impedance,
    linear, momentum_and_mass, prepared, product, quadratic, tensor, truncate)


def bend_forms(dom, block, chi, size, interval):
    """Five radial/depth integrals common to all angular harmonics."""
    roots, weights = leggauss(3)
    xi, eta = (v.ravel() for v in np.meshgrid(roots, roots, indexing='ij'))
    n = product(linear(xi), linear(eta))
    deriv = np.tile([-.5,.5], (len(xi),1))
    np_ref, nv_ref = product(deriv,linear(eta)), product(linear(xi),deriv)
    pg, vg, cells = block
    po, vo, old = dom.order_blocks[2]
    ii, jj = (v.ravel() for v in np.indices((len(pg)-1,len(vg)-1)))
    take = (pg[ii+1] > interval[0]) & (pg[ii] < interval[1])
    i, j, ids = ii[take], jj[take], cells[take]
    hp, hv = np.diff(pg)[i,None]/2, np.diff(vg)[j,None]/2
    p, v = pg[i,None]+hp*(1+xi), vg[j,None]+hv*(1+eta)
    oi = np.searchsorted(po,(pg[i]+pg[i+1])/2)-1
    oj = np.searchsorted(vo,(vg[j]+vg[j+1])/2)-1
    ax = quadratic(2*(p-po[oi,None])/np.diff(po)[oi,None]-1)[0]
    ay = quadratic(2*(v-vo[oj,None])/np.diff(vo)[oj,None]-1)[0]
    order = np.clip(np.einsum('eqi,ei->eq',tensor(ax,ay),chi[old[oi*(len(vo)-1)+oj]]),0,1)
    r,w,rp,rv,wp,wv = dom.chart(2,p,v)
    det = rp*wv-rv*wp
    dp, dv = np_ref/hp[:,:,None], nv_ref/hv[:,:,None]
    wr = r[:,:,None]*(wv[:,:,None]*dp-wp[:,:,None]*dv)/det[:,:,None]
    wz = r[:,:,None]*(rp[:,:,None]*dv-rv[:,:,None]*dp)/det[:,:,None]
    wt = (8*np.pi/3)*r**2*abs(det)*hp*hv*order*np.outer(weights,weights).ravel()
    hprime = -taper(w,*interval)[1]
    wn = np.broadcast_to(n,wr.shape)
    rows = np.repeat(ids,4,axis=1).ravel()
    cols = np.tile(ids,(1,4)).ravel()
    def form(left,right,weight):
        local = np.einsum('eq,eqi,eqj->eij',weight,left,right)
        mat = sparse.coo_matrix((local.ravel(),(rows,cols)),shape=(size,size)).tocsr()
        mat.eliminate_zeros()
        return mat
    one, two = wt*hprime, wt*hprime**2
    return [form(wz,wr,one),form(wz,wn,one),form(wr,wr,two),form(wr,wn,two),form(wn,wn,two)]


def x_derivative(out_degree, degree):
    """Coefficients of r f_r and f in the indicated angular channel of d_x(r f T_l)."""
    top = max(out_degree,degree)
    c = np.sqrt((top**2-1)/(4*top**2-1))
    return c, (degree+2 if out_degree < degree else 1-degree)


class Bend:
    def __init__(self, forms):
        self.p,self.q,self.t,self.u,self.v = forms

    def linear(self, i, j):
        if abs(i-j) != 1:
            return None
        c,a = x_derivative(i,j)
        _,b = x_derivative(j,i)
        return -c*(self.p+a*self.q+(self.p+b*self.q).T)

    def square(self, i, j):
        terms = []
        for n in sorted({i-1,i+1}&{j-1,j+1}):
            c,a = x_derivative(n,i)
            d,b = x_derivative(n,j)
            terms.append(c*d*(self.t+b*self.u+a*self.u.T+a*b*self.v))
        return sum(terms) if terms else None


def geometry(a):
    dom,chi,intake,src = prepared(a)
    dom,chi,cut = truncate(dom,chi,a.cut)
    pts,blocks,_,parts,outer,_ = assemble(dom,chi,a.sub,1,a.feed_spacing)
    masses,_ = momentum_and_mass(dom,blocks,chi,len(pts),0.)
    ang = sum(angular_stiffness(dom,blocks,chi,len(pts)))
    cut_z,_,_ = impedance(dom,blocks,pts,chi,cut,0.)
    mass = masses.sum(axis=0)
    live = (mass>0)&(pts[:,0]>0)
    if a.boundary=='waist':
        live &= pts[:,1]>.5000000001
    m = mass[live]
    scale = sparse.diags(1/np.sqrt(m))
    def restrict(mat):
        return (scale@mat[live][:,live]@scale).tocsr()
    drag = (outer+cut_z)[live]/m if a.boundary=='open' else np.zeros_like(m)
    return SimpleNamespace(dom=dom,chi=chi,block=blocks[2],size=len(pts),points=pts[live],
        k=restrict(sum(parts)),ang=restrict(ang),mass=m,cut=cut,damp=drag,
        restrict=restrict,order_sha=digest(src),intake=intake)


def pencil(k, damp, omega):
    return (k-sparse.diags(omega**2+1j*omega*damp)).tocsc()


def find_pole(k, damp, start, probe):
    """Nonlinear inverse iteration; compact localized probe keeps the same branch."""
    omega = complex(start)
    u = probe.astype(complex)
    for it in range(18):
        y = splu(pencil(k,damp,omega)).solve(-(2*omega+1j*damp)*u)
        shift = 1/(probe@y)
        u = y*shift
        omega -= shift
        if abs(shift)<3e-12:
            break
    residual = np.linalg.norm(pencil(k,damp,omega)@u)/(abs(omega)**2*np.linalg.norm(u))
    if not np.isfinite(residual) or residual > 2e-8:
        raise ArithmeticError(f'Open branch residual {residual}')
    return omega,u/np.linalg.norm(u),residual


def coupled_system(k, b, degrees, delta):
    blocks=[]
    for i in degrees:
        row=[]
        for j in degrees:
            terms=[]
            if i==j: terms.append(k[i])
            l,h = b.linear(i,j), b.square(i,j)
            if l is not None: terms.append(delta*l)
            if h is not None: terms.append(delta**2*h)
            row.append(sum(terms) if terms else None)
        blocks.append(row)
    return sparse.bmat(blocks,format='csr')


def run(a):
    g=geometry(a)
    k={ell:g.k+(ell*(ell+1)-2)*g.ang for ell in range(2,7)}
    behind=g.points[:,1]>.5000000001
    vals,vec=eigsh(k[4][behind][:,behind], k=1, sigma=0,
                   v0=np.ones(behind.sum()),tol=1e-12)
    seed=np.zeros_like(g.mass)
    seed[behind]=vec[:,0]
    ref=np.sqrt(vals[0])
    omega,u,res=find_pole(k[4],g.damp,ref+1e-4j,seed)
    energy=(abs(omega)**2*np.vdot(u,u).real+np.vdot(u,k[4]@u).real)/4
    emit(test='straight',mesh=a.mesh,sub=a.sub,outer=a.outer,cut=g.cut,boundary=a.boundary,
         seed_frequency=ref,open_frequency=omega.real,amplitude_loss=-omega.imag,
         seed_overlap_squared=abs(seed@u)**2,relative_residual=res,nodes=len(seed),
         intake=g.intake,order_sha256=g.order_sha)
    for interval in ((.55,1.15),(1.2,2.0)):
        b=Bend([g.restrict(m) for m in bend_forms(g.dom,g.block,g.chi,g.size,interval)])
        direct=b.square(4,4)@u
        mixed=[]; responses={}
        for l in (3,5):
            l4=b.linear(l,4)
            forcing=l4@u
            responses[l]=-splu(pencil(k[l],g.damp,omega)).solve(forcing)
            mixed.append(l4.T@responses[l])
        coeff=direct+sum(mixed)
        # Symmetric analytic pencil uses transpose; physical mean work uses conjugation.
        curvature=2*(u@coeff)/(u@((2*omega+1j*g.damp)*u))
        stress=np.vdot(u,coeff).real/(2*energy)
        emit(test='bend',interval=interval,k_over_energy_pole=curvature.real/omega.real,
             amplitude_loss_second_derivative=-curvature.imag,virtual_work_over_energy=stress,
             frozen_over_energy=np.vdot(u,direct).real/(2*energy),
             readjustment_over_energy=[np.vdot(u,v).real/(2*energy) for v in mixed],
             shape_norm_derivatives=[np.linalg.norm(responses[l]) for l in (3,5)])
        if a.fd:
            deg=list(range(2,7)) if a.five else [3,4,5]
            probe=np.zeros((len(deg),len(seed)))
            probe[deg.index(4)]=seed
            branch_start=omega
            for d in (.02,-.02,.01,.005):
                kk=coupled_system(k,b,deg,d)
                freq,shape,err=find_pole(kk,np.tile(g.damp,len(deg)),branch_start+1e-6j,probe.ravel())
                branch_start=freq
                ds=shape.reshape(len(deg),-1)
                derivative=0j
                for i,ui in zip(deg,ds,strict=True):
                    for j,uj in zip(deg,ds,strict=True):
                        ll,hh=b.linear(i,j),b.square(i,j)
                        if ll is not None: derivative+=np.vdot(ui,ll@uj)
                        if hh is not None: derivative+=2*d*np.vdot(ui,hh@uj)
                ed=(abs(freq)**2+np.vdot(shape,kk@shape).real)/4
                emit(test='displaced',interval=interval,offset=d,degrees=deg,frequency=freq.real,
                     loss=-freq.imag,frequency_secant=2*(freq.real-omega.real)/(omega.real*d*d),
                     flare_force_over_energy=-derivative.real/(4*ed),
                     curvature_from_direct_stress=derivative.real/(4*ed*d),
                     straight_overlap_squared=abs(np.vdot(u,ds[deg.index(4)]))**2,residual=err)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.set_defaults(radius=1., state=None)
    parser.add_argument('--mesh',type=int,default=34)
    parser.add_argument('--sub',type=int,default=1)
    parser.add_argument('--outer',type=float,default=12.)
    parser.add_argument('--feed-spacing',type=float,default=.09)
    parser.add_argument('--cut',type=float,default=2.15)
    parser.add_argument('--boundary',choices=('open','waist','reflect'),default='open')
    parser.add_argument('--fd',action='store_true')
    parser.add_argument('--five',action='store_true')
    run(parser.parse_args())
