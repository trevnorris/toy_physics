"""Forced DC monopole of the moving 0PN carrier at second order in rim speed.

Fixed complex Dirichlet traces. Unit sinusoidal velocity is along x, and mu=x/r.
Angular/time mean of each squared dipole cosine or sine is 1/6. The order-two
monopole determines the whole area-integrated mean particle current on each bath.
"""
import argparse
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy import sparse
from scipy.sparse.linalg import spsolve

from mouth_carrier_free_center import basis
from mouth_carrier_lateral_response import prepare, trace, weighted_gram
from mouth_flared_receiver import emit
from mouth_occupied_push_baths import digest


def second_order_source(dom, field, hc, hs):
    f, c, s = (dom.N @ a for a in (field, hc, hs))
    n = abs(f)**2
    q = np.zeros_like(f)
    for h in (c, s):
        v = 2*(f.conj()*h).real
        q += (n**3*(v*h+abs(h)**2*f) + 1.5*n**2*v*v*f)/6
    interaction = dom.N.T @ (dom.weight*q)
    advection = np.zeros_like(field)
    roots, weights = leggauss(4)
    u, v = [a.ravel() for a in np.meshgrid(roots, roots, indexing='ij')]
    for block, (pg, vg, dofs) in dom.blocks.items():
        for start in range(0,len(dofs),128):
            elem = np.arange(start, min(start+128, len(dofs)))
            i,j = np.unravel_index(elem,(len(pg)-1,len(vg)-1))
            i,j = i[:,None], j[:,None]
            hp,hv = np.diff(pg)[i]/2,np.diff(vg)[j]/2
            idx,r,_,det,b,gr,_ = basis(dom,block,pg[i]+hp*(1+u),vg[j]+hv*(1+v),i,j)
            ch, cr = (b*hc[idx]).sum(axis=-1), (gr*hc[idx]).sum(axis=-1)
            wt = 4*np.pi*r*r*hp*hv*det*np.outer(weights,weights).ravel()
            term = gr*ch[...,None] - b*(cr+2*ch/r)[...,None]
            np.add.at(advection,idx[:,0,:].ravel(),
                      (1j/12*wt[...,None]*term).sum(axis=1).ravel())
    return interaction, advection


def dc_state(dom, f, interaction, advection):
    b = dom.N
    fq = b @ f
    n=abs(fq)**2
    scalar = .5*dom.stiff + weighted_gram(b,dom.weight*(n**4-1)/4)
    w = 2*dom.weight*n**3
    cross = weighted_gram(b,w*fq.real*fq.imag)
    linear = sparse.bmat([[scalar+weighted_gram(b,w*fq.real**2),cross],
                          [cross,scalar+weighted_gram(b,w*fq.imag**2)]],format='csr')
    force = advection-interaction
    rhs=np.r_[force.real,force.imag]
    free=np.tile(dom.free,2)
    out=np.zeros_like(rhs)
    out[free]=spsolve(linear[free][:,free].tocsc(),rhs[free])
    residual=linear @ out-rhs
    k=np.asarray(out[:len(f)]+1j*out[len(f):])
    error=residual[:len(f)]+1j*residual[len(f):]
    outward=2*(f.conj()*error).imag
    emit(test='dc', free_residual=max(abs(error[dom.free])),
         max_mean_field_per_speed2=max(abs(k)), feed_out=outward[dom.feed].sum(),
         bulk_out=outward[dom.bulk].sum(), weak_sum=outward[~dom.free].sum())
    if max(abs(error[dom.free]))>1e-8:
        raise ArithmeticError('Mean-field solve failed')
    return k


def check_sections(dom, fields):
    # Direct totals through physical horizontal disks independent of nodal phase reactions.
    roots, weights=leggauss(8)
    vgrid=dom.blocks[2][1]
    half=np.diff(vgrid)/2
    v=(vgrid[:-1,None]+half[:,None]*(roots+1)).ravel()
    dv=(half[:,None]*weights).ravel()
    values=[]
    for w in (.7,1.4,2.4,3.0):
        r,_,_,rv,_,_,g,_,gw=trace(dom,2,np.full_like(v,w),v,fields)
        f,c,s,k=g.T
        fw,cw,sw,kw=gw.T
        mean_field=(k.conj()*fw+f.conj()*kw).imag
        carried=((c.conj()*cw+s.conj()*sw)/6).imag
        measure=4*np.pi*r*r*abs(rv)*dv
        values.append(dict(w=w, mean_field=measure @ mean_field, dipole=measure @ carried,
                           sum=measure @ (mean_field+carried)))
    emit(test='cross_section_outward_toward_bulk', values=values)


def finite_ports(dom, fields, speed0, omega, force0):
    """Direct current without its quadratic expansion on the computed trial state.

    g=f+A mu(c cos+s sin)+A² k. Higher angular and temporal harmonics of order
    two vanish from area-integrated DC particle flux to this order; they and
    cubic dipoles can affect quartic power, so these are observable checks on
    a truncated field, not finite-amplitude nonlinear trajectory solutions.
    """
    angular, angular_w = leggauss(8)
    roots, weights=leggauss(8)
    phases=np.arange(24)*(2*np.pi/24)
    for side, block, along, val, direction in (
        ('feed',1,'v',dom.far,(1,0)),('bulk_side',4,'p',1.,(1,0)),
        ('roof_inner',3,'v',dom.roof,(0,1)),('roof_outer',4,'v',dom.roof,(0,1))):
        pg,vg,_=dom.blocks[block]
        edges=pg if along=='p' else vg
        half=np.diff(edges)/2
        q=(edges[:-1,None]+half[:,None]*(roots+1)).ravel()
        fixed=np.full_like(q,val)
        p,v=(q,fixed) if along=='p' else (fixed,q)
        r,_,rp,rv,wp,wv,f,fr,fw=trace(dom,block,p,v,fields)
        tr,tw=(rp,wp) if along=='p' else (rv,wv)
        nr,nw=tw,-tr
        sg=np.sign(nr*direction[0]+nw*direction[1])
        norm=np.hypot(nr,nw)
        nr,nw=sg*nr/norm,sg*nw/norm
        measure=4*np.pi*r*r*norm*(half[:,None]*weights).ravel()
        dn=fr*nr[:,None]+fw*nw[:,None]
        b,c,s,k=dn.T
        base=f[:,0]
        if max(abs(f[:,1:]).ravel())>1e-10:raise ValueError('Far trace varied')
        branch=[]
        for fraction in (0,1,-1,.5,-.5):
            a=speed0*fraction
            j0=pwr=0.
            for ph in phases:
                co,si=np.cos(ph),np.sin(ph)
                for mu,wmu in zip(angular, angular_w/2):
                    delta=a*mu*(c*co+s*si)+a*a*k
                    grad=b+delta
                    # On a constant Dirichlet surface only the normal gradient
                    # survives: in the moving frame e-flux = V n_x(K-W).
                    nx=mu*nr
                    energy=(.5*abs(grad)**2-(abs(base)**10-5*abs(base)**2+4)/20)
                    pwr+=wmu/len(phases)*a*co*(measure @ (nx*energy))
                    j0+=wmu/len(phases)*(measure @ ((base.conj()*delta).imag-a*co*nx*abs(base)**2))
            branch.append(dict(fraction=fraction,force=force0*fraction,excess_number=j0,
                               shifted_energy=pwr))
        emit(test='truncated_field_finite_ports', side=side, omega=omega, values=branch)


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--source',type=Path,required=True)
    parser.add_argument('--response',type=Path,required=True)
    parser.add_argument('--force',type=float,default=.00001)
    args=parser.parse_args()
    dom,f=prepare(args)
    with np.load(args.response) as a:
        np.testing.assert_array_equal(dom.points,a['points'])
        assert str(a['background_sha256'])==digest(args.source/'field.npz')
        ar,ai=a['harmonic_real'],a['harmonic_imag']
        hc,hs=ar.real+1j*ai.real,-ar.imag-1j*ai.imag
        speed=args.force*abs(a['velocity_per_force'].item())
        omega=a['omega'].item()
    emit(test='input', response_sha256=digest(args.response), speed=speed, omega=omega)
    interaction, advection=second_order_source(dom,f,hc,hs)
    k=dc_state(dom,f,interaction,advection)
    fields=np.column_stack((f,hc,hs,k))
    check_sections(dom,fields)
    finite_ports(dom,fields,speed,omega,args.force)


if __name__=='__main__':
    main()
