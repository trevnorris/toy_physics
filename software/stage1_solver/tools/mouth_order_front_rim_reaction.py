"""Wave reaction on a locally actuated rounded lip while the deeper order front shifts.

The rim actuator carries its local ordered material in w. Its compactly supported
map deforms both geometry and the material kinetic term. Remote ports and the
entire displaced conversion front remain fixed during each rim virtual step.
"""
import argparse
from types import SimpleNamespace

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy import sparse
from scipy.sparse.linalg import eigsh

from mouth_order_front_wave_force import Front, smooth
from mouth_native_wave_bend import find_pole
from mouth_shear_pattern_release import (angular_stiffness, assemble, digest, emit,
    impedance, linear, momentum_and_mass, prepared, product, quadratic, tensor, truncate)


def window(x, start, flat_start, flat_end, end):
    def slope(z):
        return 30*np.clip(z,0,1)**2*(1-np.clip(z,0,1))**2
    a = (x-start)/(flat_start-start)
    b = (end-x)/(end-flat_end)
    return (smooth(a)*smooth(b),
            slope(a)*smooth(b)/(flat_start-start)-smooth(a)*slope(b)/(end-flat_end))


def actuator(r, w):
    radial, dr = window(r,.25,.48,1.10,1.38)
    axial, dw = window(w,-.32,-.04,.48,.58)
    return radial*axial, dr*axial, radial*dw


def deformed(dom, y):
    def chart(block,p,v):
        r,w,rp,rv,wp,wv = dom.chart(block,p,v)
        h,hr,hw = actuator(r,w)
        if np.min(1+y*hw) <= 0:
            raise ValueError('Non-invertible virtual displacement')
        return r,w+y*h,rp,rv,(1+y*hw)*wp+y*hr*rp,(1+y*hw)*wv+y*hr*rv
    return SimpleNamespace(**dict(vars(dom),chart=chart))


def rim_work(dom, blocks, chi, wave, omega, energy, degree):
    roots, weights = leggauss(3)
    xi, eta = (x.ravel() for x in np.meshgrid(roots,roots,indexing='ij'))
    shape = product(linear(xi),linear(eta))
    dlin = np.tile([-.5,.5],(len(xi),1))
    sp, sv = product(dlin,linear(eta)),product(linear(xi),dlin)
    strain_deriv = kinetic_deriv = 0.
    for b in (0,1,2):
        p,v,cells = blocks[b]
        p0,v0,old = dom.order_blocks[b]
        i,j = (x.ravel() for x in np.indices((len(p)-1,len(v)-1)))
        if b==1:
            take = p[i]<1.38
        elif b==2:
            take = p[i]<.58
        else:
            take = np.ones_like(i,dtype=bool)
        i,j,ids = i[take],j[take],cells[take]
        hp,hv = np.diff(p)[i,None]/2,np.diff(v)[j,None]/2
        pp,vv = p[i,None]+hp*(1+xi),v[j,None]+hv*(1+eta)
        oi = np.searchsorted(p0,(p[i]+p[i+1])/2)-1
        oj = np.searchsorted(v0,(v[j]+v[j+1])/2)-1
        aa = quadratic(2*(pp-p0[oi,None])/np.diff(p0)[oi,None]-1)[0]
        bb = quadratic(2*(vv-v0[oj,None])/np.diff(v0)[oj,None]-1)[0]
        order = np.clip(np.einsum('eqi,ei->eq',tensor(aa,bb),
                                  chi[old[oi*(len(v0)-1)+oj]]),0,1)
        r,w,rp,rv,wp,wv = dom.chart(b,pp,vv)
        _, hr, hw = actuator(r,w)
        det = rp*wv-rv*wp
        local = wave[ids]
        f = local@shape.T
        fp,fv = (local@sp.T)/hp,(local@sv.T)/hv
        fr,fw = (wv*fp-wp*fv)/det,(rp*fv-rv*fp)/det
        radial_curl = 2*f+r*fr
        strain = degree*(degree+1)*abs(f)**2+abs(radial_curl)**2+r**2*abs(fw)**2
        # At fixed pulled-back wave, d f_r = -h_r f_w, d f_w = -h_w f_w.
        ds = hw*strain-2*r*hr*(radial_curl.conj()*fw).real-2*r*r*hw*abs(fw)**2
        weight = (8*np.pi/3)*r**2*np.abs(det)*hp*hv*np.outer(weights,weights).ravel()*order
        strain_deriv += np.sum(weight*ds)
        kinetic_deriv += np.sum(weight*hw*r**2*(np.abs(local)**2@shape.T)*abs(omega)**2)
    scale = omega.real/(4*energy)
    return -scale*strain_deriv, scale*kinetic_deriv


def mode(dom, chi, cut, a, compute_work):
    pts, blocks, _, parts, feed, _ = assemble(dom,chi,a.sub,1,a.feed_spacing)
    masses,_ = momentum_and_mass(dom,blocks,chi,len(pts),0.)
    k = sum(parts)+(a.degree*(a.degree+1)-2)*sum(angular_stiffness(dom,blocks,chi,len(pts)))
    outlet,_,_ = impedance(dom,blocks,pts,chi,cut,0.)
    full_mass = masses.sum(axis=0)
    live = (full_mass>0)&(pts[:,0]>0)
    m = full_mass[live]
    scale = sparse.diags(1/np.sqrt(m))
    op = scale@k[live][:,live]@scale
    behind = pts[live,1]>.5000000001
    val,vec = eigsh(op[behind][:,behind],k=1,sigma=0,tol=1e-12,v0=np.ones(behind.sum()))
    seed = np.zeros_like(m)
    seed[behind] = vec[:,0]
    omega, u, residual = find_pole(op,(feed+outlet)[live]/m,np.sqrt(val[0])+1.e-4j,seed)
    if residual>2e-8 or abs(seed@u)**2<.90 or not (0<=-omega.imag<.03*omega.real):
        raise ArithmeticError('Confined branch identification failed')
    data = dict(frequency=omega.real,loss=-omega.imag,residual=residual,
                waist_seed_overlap=abs(seed@u)**2)
    if compute_work:
        wave = np.zeros(len(pts),dtype=complex)
        wave[live] = u/np.sqrt(m)
        energy = (np.vdot(u,op@u).real+abs(omega)**2)/4
        shear, inertia = rim_work(dom,blocks,chi,wave,omega,energy,a.degree)
        data.update(rim_force_per_action=shear+inertia,shear_part=shear,kinetic_part=inertia)
    return data


def run(a):
    raw, base, intake, src = prepared(a)
    front = Front(raw,base)
    emit(test='rounded_lip', radius=a.radius, degree=a.degree, mesh=a.mesh, sub=a.sub,
         intake=intake, order_source=str(src), order_sha256=digest(src), actual_lip_r=1.,
         actuator_radial_window=[.25,.48,1.10,1.38], actuator_axial_window=[-.32,-.04,.48,.58],
         outer=a.outer, feed_spacing=a.feed_spacing, cut_request=a.cut, rim_step=a.rim_step)
    rows = {}
    for q in (0.,-.02,-.01,.01,.02):
        field,_ = front.fields(q)
        dom, chi, cut = truncate(raw,field,a.cut)
        signed = []
        for y in (-a.rim_step, 0., a.rim_step):
            measurement = mode(deformed(dom,y),chi,cut,a,y==0.)
            emit(test='wave',q=q,rim_virtual_shift=y,**measurement)
            signed.append(measurement)
        fn = -(signed[2]['frequency']-signed[0]['frequency'])/(2*a.rim_step)
        direct = signed[1]['rim_force_per_action']
        emit(test='reaction',q=q,cut=cut,stress_force=direct,frequency_force=fn,
             empty_wave_force=0.,fixture_force=-direct)
        if not np.isclose(direct,fn,rtol=.035,atol=5e-7):
            raise ArithmeticError('Rim stress and frequency slope disagree')
        rows[q] = (direct,fn)
    for h in (.02,.01):
        change = (np.array(rows[h])-np.array(rows[-h]))/(2*h)
        emit(test='cross_reaction',front_excursion=h,
             rim_force_change_per_front_shift=change[0],frequency_mixed_force_derivative=change[1],
             base_rim_force=rows[0.][0])


if __name__=='__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.set_defaults(radius=1.,degree=4,state=None)
    p.add_argument('--mesh',type=int,default=34)
    p.add_argument('--sub',type=int,default=1)
    p.add_argument('--outer',type=float,default=12.)
    p.add_argument('--cut',type=float,default=2.15)
    p.add_argument('--feed-spacing',type=float,default=.09)
    p.add_argument('--rim-step',type=float,default=.002)
    run(p.parse_args())
