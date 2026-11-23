"""Axisymmetric shape stiffness, registered lattice mass and open sound channels on both reference mouths."""
import json
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss, legval
from scipy.interpolate import RegularGridInterpolator
from scipy.linalg import eigh
from scipy.special import spherical_jn

from mouth_exit_witness import inputs
from mouth_poly_shallow_release import CUT, ROOT, Trial, emit, source
from mouth_sealed_shape_relax import Meridian
from mouth_waist_outer_order import exterior


def coordinates(trial, shifts, b):
    x = shifts.copy()
    base = Meridian(trial.radius, trial.x + shifts[:4])
    x[1] += np.log(b / base.exit)
    return x, Meridian(trial.radius, trial.x + x[:4])


def load_map(trial, shifts, chi):
    shape = Meridian(trial.radius, trial.x+shifts[:4])
    c = trial.co.copy()
    c[0] += shifts[4]
    cr = trial.cr.copy()
    cr[0] += shifts[5]

    def receiver_load(w, x):
        return source(shape, cr, w/shape.length, x, True)

    energy, (w, x, order) = exterior(shape, trial.n, CUT, return_side_field=True,
                                     side_order_work=receiver_load)
    inside = RegularGridInterpolator((np.linspace(0,1,2*trial.n+1), np.linspace(0,1,trial.n+1)), chi)
    outside = RegularGridInterpolator((w/shape.length, x),order)
    return shape,c,cr,inside,outside,energy


def measure_mass(trial, fields, exit_scale, n):
    shape,c,cr,inside,outside,_ = fields
    original,co = trial.core_reference
    z,wt = leggauss(n)
    y,wy = leggauss(n)
    t=(z[:,None]+1)/2
    y=(y[None,:]+1)/2
    weight=wt[:,None]*wy[None,:]/4
    current_b = shape.profile(t)[0]
    cut0 = original.profile(t)[0]*np.exp(-.5*np.exp(legval(2*t-1,co)))
    matrices={}
    for name in ('core','inner_collar','outer_collar'):
        if name=='core':
            x = np.exp(-.5*np.exp(legval(2*t-1,c)))*y
            rest_r = cut0*y
            measure = weight*4*np.pi*rest_r**2*cut0*original.length
        else:
            # Signed a is distance from the midplane of the initially flat slab.
            side = 1 if name=='inner_collar' else -1
            a = side*y/2
            x = np.exp(-np.exp(legval(2*t-1,c if side==1 else cr))*a)
            measure = weight*2*np.pi*(trial.radius+a)**2
        th = np.broadcast_to(t, x.shape)
        if name=='outer_collar':
            xouter=(x-1)*current_b/(CUT-current_b)
            order=outside(np.stack([th,xouter],axis=-1))
        else:
            order=inside(np.stack([th,x],axis=-1))
        coeff=measure*order
        # Basis {log(depth), exit/initial_exit, middle}. b(t) depends linearly
        # on the physical exit and middle; material depth is affine in depth.
        sides=np.stack([np.zeros_like(x),x*t*exit_scale,x*4*trial.radius*t*(1-t)],axis=-1)
        mass=np.einsum('mn,mni,mnj->ij',coeff,sides,sides)
        mass[0,0]+=np.sum(coeff*(shape.length*t)**2)
        matrices[name]=mass
    return matrices


def volume_derivatives(shape, exit_scale):
    nodes,w=leggauss(128)
    t=(nodes+1)/2
    b=shape.profile(t)[0]
    V=shape.length*(w@(4*np.pi*b**3/3))/2
    dV=shape.length* (w*b*b*2*np.pi) @ np.column_stack([exit_scale*t,4*shape.neck*t*(1-t)])
    return np.r_[V,dV]


def exterior_probe(omega, radius):
    # Positive-frequency Bogoliubov mode in the homogeneous resting bath.
    # k^2 solves Omega^2=k^2+k^4/4. The normalized source is a uniform
    # chemical-potential perturbation on a 3-ball in a 4D infinite bath.
    # The angular average of the on-shell squared Fourier source is the
    # emitted power relative to a point source of the same integrated strength.
    k=np.sqrt(2*omega**2/(1+np.sqrt(1+omega**2)))
    mat=np.array([[1+k*k/2,1],[-1,-1-k*k/2]])
    if abs(max(np.linalg.eigvals(mat))-omega)>1e-10:
        raise ArithmeticError('Acoustic spectral check failed')
    def form(n):
        x,w=leggauss(n)
        ang=(x+1)*np.pi/2
        q=k*radius*np.sin(ang)
        disk=3*spherical_jn(1,q)/q
        return float(np.sum(w*(np.sin(ang)*disk)**2))
    near,far=form(80),form(160)
    if abs(near-far)>1.e-11 or not 0 < far <= 1+1e-12:
        raise ArithmeticError('Open sound angular integral failed')
    return dict(k_bulk_retarded=k,bulk_unit_disk_monopole_relative_power=far)


def run(radius,n):
    trial=Trial(radius,n,bonded_core=True)
    path=ROOT/f'_scratch/s11c/mouth_exit_release01/r{radius}n{n}/endpoint.npz'
    shifts,chi,seed,sha,status,*_=inputs(path,radius,n,trial)
    if status!='balanced':
        raise ValueError('Unbalanced source in its four released coordinates')
    rows=list(map(json.loads,(path.parent/'stdout').read_text().splitlines()))
    init,endpoint=rows[0],rows[-1]
    stiffness=np.array([row['hessian'] for row in rows if row['kind']=='four_coordinate_curvature'])
    sch=init['exit_radius']
    base=Meridian(radius,trial.x+shifts[:4])
    # Distinct physical exit probes; the original four-coordinate released point
    # supplies the center energy and all symmetric degree-four comparisons.
    center=sum(endpoint['energy'].values())
    cache={}
    def state(offset):
        if offset in cache: return cache[offset]
        shifted,shape=coordinates(trial,shifts,base.exit*(1+offset))
        fields=load_map(trial,shifted,chi)
        masses=measure_mass(trial,fields,sch,112)
        check=measure_mass(trial,fields,sch,160)
        m=sum(check.values())
        mass_error=np.linalg.norm(m-sum(masses.values()))/np.linalg.norm(m)
        e=center if offset==0 else trial.evaluate(shifted,chi,seed)[0]
        if mass_error>.003 or min(np.linalg.eigvalsh(m))<=0:
            raise ArithmeticError('Lattice mass quadrature unresolved')
        cache[offset]=e,m
        emit(kind='basis_state',radius=radius,n=n,relative_exit=offset,energy=e,
             mass=m.tolist(),relative_quadrature_error=mass_error,shape_volume_slopes=volume_derivatives(shape,sch).tolist(),
             receiving_order_residual=fields[-1]['residual'])
        return e,m

    E,M=state(0.)
    hessian=stiffness[-1,:3,:3]
    vals,vec=eigh(hessian,M)
    if min(vals)<=0: raise ArithmeticError('Nonpositive projected shape stiffness')
    grad_vol=volume_derivatives(base,sch)
    emit(kind='setup',radius=radius,n=n,endpoint_sha256=sha,physical_exit=base.exit,
         depth=base.length,cap_three_area=4*np.pi*base.exit**3/3,
         order4_frequency=float(seed[0].real),current=endpoint['current'],outflow=endpoint['outflow'])
    for index,(value,mode) in enumerate(zip(vals,vec.T)):
        freq=np.sqrt(value)
        emit(kind='three_shape_trial',radius=radius,n=n,index=index,freq=freq,
             basis=mode.tolist(),abs_dVdq_per_sqrt_mass=abs(grad_vol@mode),
             physical_exit_amplitude_per_sqrt_mass=abs(mode[1]*sch),
             **exterior_probe(freq,base.exit))

    # Report both step sizes; measurements are finite differences of the stored
    # potential and participating lattice mass. No pressure impedance is fitted.
    for h in (.012,.006):
        freq=[]
        for sign in (-1,0,1):
            u=sign*h
            eo,m=state(u)
            ep=state((sign+1)*h)[0]
            em=state((sign-1)*h)[0]
            curv=(ep+em-2*eo)/(base.exit*h)**2
            omega=np.sqrt(curv/(m[1,1]/sch**2))
            if not np.isfinite(omega): raise ArithmeticError('Nonpositive breathing scale')
            freq.append(omega)
        emit(kind='one_exit_trial',radius=radius,n=n,step=h,frequencies=freq,
             domega_dexit=(freq[2]-freq[0])/(2*h*base.exit),
             archived_center=np.sqrt(hessian[1,1]/M[1,1]),
             **exterior_probe(freq[1],base.exit))


if __name__=='__main__':
    for radius in (1,2):
        for n in (48,72):
            run(radius,n)
