"""Signed incoming material in the upstream radial well, followed to the true outlet plane.

Full spherical angular slip, prescribed inward w clock, unmodified quadratic
order on its original elements. The mathematical well is two-sided and has
no unilateral wall collision. Aperture admission is explicitly diagnosed.
"""
import argparse

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.integrate import solve_ivp

from mouth_calibrated_flare_coating import read_inputs
from mouth_flared_receiver import EXIT
from mouth_middle_tag_grip import cubic_trace
from mouth_segment_contact_hold import PITCH
from neutral_patch_sites import BETWEEN, TEMPERATURE, emit


def field_on_element(tr, element, w):
    start, end = tr.edges[element:element+2]
    f=(w-start)/(end-start)
    def quadratic(nodes):
        l,m,r=nodes
        b=4*m-3*l-r
        a=2*(l+r-2*m)
        return l+f*(b+f*a),(b+2*a*f)/(end-start)
    R,slope=tr.geometry.radius(w)
    psi,gradient=quadratic(tr.psi[element])
    speed=(gradient/psi).imag/(1+slope*slope)
    order=quadratic(tr.chi[element])[0]
    if order<=0 or speed<=0:
        raise ValueError(f'Nonpositive native order or clock at {w}: {order}, {speed}')
    return float(R),float(order),float(speed)


def integrate(args):
    tr,intake=cubic_trace(args.mesh)
    mat,_,_=read_inputs(args.mesh)
    scale=2*mat['shared_scale']
    D,k,rest=BETWEEN
    sigma=np.sqrt(scale*TEMPERATURE)
    lo=1.175
    R0,R0p=map(float,tr.geometry.radius(lo))
    base_speed=R0p*tr.speed(lo)
    mu,angular_weight=leggauss(args.angles)
    angular_weight=angular_weight[:,None]/2
    sin0=np.sqrt(1-mu**2)[:,None]
    mu=mu[:,None]
    U=args.speed
    n=args.velocities
    if args.launch=='cold':
        normal_vel=np.zeros((1,1))
        weight=angular_weight
    else:
        grid,normal_weight=leggauss(n)
        normal_vel=9*sigma*grid[None,:]
        # On a lab-velocity grid the shifted radial launch distribution implements
        # an incident U*e, including its radial part. Controls also integrate
        # the identical physical preparation on an unshifted, source-velocity grid.
        exponent=normal_vel if args.launch=='source' else normal_vel-U*mu
        weight=angular_weight*normal_weight[None,:]*np.exp(-.5*(exponent/sigma)**2)
        weight/=weight.sum()
    vr0=base_speed+normal_vel+(U*mu if args.launch!='lab' else 0.)
    vr0=np.broadcast_to(vr0,weight.shape)
    L=R0*U*sin0
    state=np.zeros((5,*weight.shape))
    state[0]=R0
    state[1]=vr0
    initial_vx=vr0*mu+U*sin0**2
    mean_in=np.sum(weight*initial_vx)
    if abs(mean_in-U)>2e-14:
        raise ArithmeticError(f'Incorrect incident weak drift ({mean_in}, {U})')
    gates=np.unique([tr.crossing(.01),1.55,1.75,2.,2.25,2.5,2.75,3.,tr.edges[-1]])
    if abs(gates[-1]-EXIT)>1e-12: raise ArithmeticError('Wrong physical exit')
    edges=np.unique(np.r_[lo,tr.edges[(tr.edges>lo)&(tr.edges<gates[-1])],gates])
    left_support=np.zeros(weight.shape,dtype=bool)
    negative_length=np.zeros(weight.shape,dtype=bool)
    min_length=rest
    emit(test='signed_setup',mesh=args.mesh,intake=intake,U=U,launch=args.launch,
         radial_nodes=weight.shape[1],angle_nodes=args.angles,R0=R0,incoming_vx=mean_in,
         thermal_sigma=sigma,actual_exit=EXIT,far_radius=tr.geometry.far)
    for a,b in zip(edges[:-1],edges[1:],strict=True):
        element=np.searchsorted(tr.edges,(a+b)/2)-1
        def rhs(w,y):
            r,v,theta,_,_=y.reshape(state.shape)
            R,order,speed=field_on_element(tr,element,w)
            # Compute chi*exp(-2*k*gap/ell) directly in log form in the deep
            # tail, retaining both the attractive and the compressive terms.
            log_x=-k*(r-R)/PITCH
            log_chi=np.log(order)
            acc=2*D*k*scale/PITCH*(np.exp(log_chi+2*log_x)-np.exp(log_chi+log_x))
            ang=L/r**2
            return (np.array((v,acc+L*ang/r,ang,acc*np.cos(theta),acc*np.sin(theta)))/speed).ravel()
        sol=solve_ivp(rhs,(a,b),state.ravel(),method='DOP853',rtol=args.rtol,atol=2e-13)
        if not sol.success or not np.isfinite(sol.y).all():
            raise ArithmeticError(sol.message)
        state=sol.y[:,-1].reshape(state.shape)
        # Sampled excursions only; endpoints in particular are directly observed.
        radii=sol.y[:weight.size].reshape(*weight.shape,-1)
        gaps=radii-tr.geometry.radius(sol.t)[0]
        left_support |= np.max(gaps,axis=-1)>1e-8
        d=rest+np.min(gaps,axis=-1)/PITCH
        negative_length |= d<0
        min_length=min(min_length,d.min())
        if b not in gates: continue
        r,v,theta,Jc,Js=state
        R,order,_=field_on_element(tr,element,b)
        cos,sin=np.cos(theta),np.sin(theta)
        nx=mu*cos+sin0*sin
        tx=sin0*cos-mu*sin
        vx=v*nx+(L/r)*tx
        on_well=-np.sum(weight*(mu*Jc+sin0*Js))
        inside=r<=R
        out_in=np.sum(weight*vx*inside)
        out_beyond=np.sum(weight*vx*~inside)
        balance=mean_in-out_in-out_beyond-on_well
        emit(test='signed_discharge_cut',depth=b,radius=R,chi=order,incoming_vx=mean_in,
             on_well_per_tag_flux=on_well,external_holder_per_tag_flux=-on_well,
             flux_gain_estimate=mean_in-np.sum(weight*vx),momentum_residual=balance,
             normalized_on_well=on_well/U if U else None,inside_site_flux=np.sum(weight*inside),
             vx_flux_inside=out_in,vx_flux_outside=out_beyond,
             minimum_gap=np.min(r-R),maximum_gap=np.max(r-R),
             ever_outside_at_sampled_points=np.sum(weight*left_support),
             ever_negative_pair_length=np.sum(weight*negative_length),minimum_pair_length=min_length)
        if abs(balance)>1e-9:
            raise ArithmeticError(f'Tag momentum identity fails at {b}: {balance}')
    if not np.any(state[0]<=R):
        emit(test='receiver_admission',entered_nodes=0,formal_continuation='outside the modeled opening',
             actual_downstream_tag_force=None)
    else:
        emit(test='receiver_admission',entered_nodes=int((state[0]<=R).sum()),
             entering_weight=np.sum(weight*(state[0]<=R)),actual_downstream_tag_force=None)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--mesh',type=int,required=True)
    p.add_argument('--speed',type=float,required=True)
    p.add_argument('--velocities',type=int,required=True)
    p.add_argument('--angles',type=int,default=4)
    p.add_argument('--launch',choices=['cold','lab','source'],default='lab')
    p.add_argument('--rtol',type=float,default=1e-10)
    integrate(p.parse_args())
