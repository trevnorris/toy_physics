"""Compare wall traction, interior momentum balance, and low-frequency acceleration.

Both remote material reservoirs rest in the lab. The prescribed rigid mouth and
finite bath geometry translate at common velocity V. In the accelerating frame
the extra scalar potential is a*x and remote phases acquire -V*x.
"""
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy import sparse
from scipy.sparse.linalg import spsolve

from mouth_carrier_lateral_response import (Dipole, emit, parser, prepare, surfaces,
                                            trace, weighted_gram)
from mouth_flared_receiver import EXIT


def taper(x,start,end):
    t=np.clip((x-start)/(end-start),0,1)
    return 1-t**3*(10-15*t+6*t*t), -30*t*t*(1-t)**2/(end-start)


def volume(dom,fields,window,accelerating,*,rest_energy=False):
    """Virtual translation supported around the reflecting throat. No corner trace."""
    root,weights=leggauss(4)
    u,v=[x.ravel() for x in np.meshgrid(root,root,indexing='ij')]
    applied=storage=kinetic=0.
    for block,(pg,vg,dofs) in dom.blocks.items():
        for start in range(0,len(dofs),160):
            i,j=np.unravel_index(np.arange(start,min(start+160,len(dofs))),
                                 (len(pg)-1,len(vg)-1))
            hp,hv=np.diff(pg)[i,None]/2,np.diff(vg)[j,None]/2
            p=pg[i,None]+hp*(1+u)
            s=vg[j,None]+hv*(1+v)
            r,w,rp,rv,wp,wv,f,fr,fw=trace(dom,block,p,s,fields)
            wt=4*np.pi/3*r*r*hp*hv*abs(rp*wv-rv*wp)*np.outer(weights,weights).ravel()
            chi_r,dr=taper(r,*window[:2])
            chi_w,dw=taper(w,*window[2:])
            n=abs(f[...,0])**2
            lag=(fr[...,0].conj()*fr[...,1]+fw[...,0].conj()*fw[...,1]).real
            lag+=(n**4-1)/2*(f[...,0].conj()*f[...,1]).real
            if accelerating:
                lag+=n*r+(f[...,0].conj()*f[...,2]).imag
                current=(f[...,0].conj()*(fr[...,2]+2*f[...,2]/r)+f[...,2].conj()*fr[...,0]).imag
                storage+=np.sum(wt*chi_r*chi_w*(3*n+current))
            grad=fr[...,1]+2*f[...,1]/r
            radial=(fr[...,0].conj()*fr[...,1]+grad.conj()*fr[...,0]).real-lag
            axial=(fr[...,0].conj()*fw[...,1]+grad.conj()*fw[...,0]).real
            applied+=np.sum(wt*(radial*dr*chi_w+axial*chi_r*dw))
            # Independent positive kinetic norm in a homogeneous, nondraining bath.
            if rest_energy:
                potential=f[...,1].imag
                kinetic+=np.sum(wt*((1+fr[...,1].imag)**2+fw[...,1].imag**2+2*(1+potential/r)**2))
    reading=dict(virtual_wall_force=applied-storage,stress_term=applied,lab_momentum_rate=storage)
    if rest_energy:
        reading['rest_potential_norm']=kinetic
    return reading


WINDOWS=[(1.7,4.3,EXIT+.4,EXIT+2.5),(2.3,6.,EXIT+.7,EXIT+3.4)]


def independent_frequency(model,velocity,omega):
    dom=model.dom
    mass=weighted_gram(dom.N,dom.weight)
    skew=sparse.bmat([[None,-mass],[mass,None]],format='csr')
    mat=model.operator.astype(complex)-1j*omega*skew
    bath= -1j*model.field*dom.points[:,0]*(dom.feed|dom.bulk)
    out=np.r_[bath.real,bath.imag].astype(complex)
    load=dom.N.T@(dom.weight*model.rq*(dom.N@model.field))
    rhs= -1j*omega*np.r_[load.real,load.imag]
    rhs-=mat@out
    out[model.free]=spsolve(mat[model.free][:,model.free].tocsc(),rhs[model.free])
    vr,vi=np.split(out,2)
    in_phase=vr.real+1j*vi.real
    quadrature=vr.imag+1j*vi.imag
    fields=np.column_stack([model.field,quadrature/omega,in_phase])
    # At this frequency Im(dot h)/omega is the in-phase velocity-drive field.
    wall,_,_=surfaces(dom,fields,accelerating=True)
    projected=volume(dom,fields,WINDOWS[0],True)
    return dict(omega=omega,actuator_mass_wall=-sum(x[1] for x in wall.values()),
                actuator_mass_volume=-projected['virtual_wall_force'],
                field_deviation_from_dc=float(max(abs(in_phase-velocity))))


def experiment(args):
    dom,field=prepare(args)
    model=Dipole(dom,field)
    feed,_=model.steady(1.,0.)
    bulk,_=model.steady(0.,1.)
    velocity,_=model.steady(-1.,-1.)
    for name,h in [('feed',feed),('bulk',bulk),('tow',velocity)]:
        fields=np.column_stack([field,h])
        windows=[volume(dom,fields,w,False) for w in WINDOWS]
        wall,outer,_=surfaces(dom,fields)
        emit(test='steady_balance',drive=name,direct_wall=sum(x[1] for x in wall.values()),
             outer_cut=-sum(outer.values()),interior=[x['virtual_wall_force'] for x in windows])
    source=dom.N.T@(dom.weight*(1j*(dom.N@velocity)-model.rq*(dom.N@field)))
    acceleration,err=model.solve(np.zeros_like(field),source)
    fields=np.column_stack([field,acceleration,velocity])
    wall,outer,cuts=surfaces(dom,fields,accelerating=True)
    windows=[volume(dom,fields,w,True) for w in WINDOWS]
    mass=-windows[0]['virtual_wall_force']
    direct=-sum(x[1] for x in wall.values())
    emit(test='slow_acceleration',actuator_mass_volume=mass,actuator_mass_wall=direct,
         windows=windows,wall_parts={k:-v[1] for k,v in wall.items()},
         dipole_residual=err,section_flux=cuts,outer_cut=outer)
    if max(abs(field-1))<1e-14:
        norm=volume(dom,np.column_stack([field,velocity]),WINDOWS[0],False,rest_energy=True)['rest_potential_norm']
        emit(test='no_drain_kinetic_check',positive_kinetic_norm=norm,
             wall_mass=direct,volume_mass=mass)
        if max(abs(norm-direct),abs(norm-mass))>.02*norm:
            raise ArithmeticError('Resting inertia disagrees with independent kinetic norm')
    for omega in args.omega:
        emit(test='harmonic',**independent_frequency(model,velocity,omega))


if __name__=='__main__':
    p=parser()
    p.add_argument('--omega',type=float,nargs='*',default=[])
    experiment(p.parse_args())
