"""Native coherent vortex on a neutral order interface and order/strain probes on the held mouths."""
import json

import numpy as np
from numpy.polynomial.legendre import leggauss, legder, legval
from scipy.integrate import solve_bvp
from scipy.optimize import brentq

from mouth_curved_steady import SteadyMouth
from mouth_exit_witness import inputs
from mouth_finite_bore_support import graph_values
from mouth_live_waist_fields import BIAS, Carrier, ORDER_INERTIA, WIDTH
from mouth_poly_bonded_core import core_source
from mouth_poly_shallow_release import ROOT, Order, Trial, emit, source
from mouth_sealed_shape_relax import Meridian


def profile(m, outer, tolerance):
    # f=(s/sqrt(1+s²))^m g. The first regular core correction and leading
    # analytic bath correction replace arbitrary core radii or imposed speed caps.
    m=abs(m)
    if m<1:
        raise ValueError('Use the homogeneous state as the zero-winding control')
    inner=1e-5
    s=np.unique(np.r_[np.geomspace(inner,1.,90),np.linspace(1,outer,180)])
    def envelope(x):
        return (x/np.hypot(1,x))**m
    def rhs(x,y):
        g,dg=y
        h=envelope(x)
        correction=m*(2*(1+m)+m*x*x)/(1+x*x)**2
        return np.array([dg,(.5*((h*g)**8-1)+correction)*g
                         -(1+2*m/(1+x*x))*dg/x])
    def ends(lo,hi):
        return np.array([lo[1]-inner*(m-1/(4*(m+1)))*lo[0],
                         hi[0]*envelope(outer)-(1-m*m/(4*outer**2))])
    guess=(np.hypot(1,s)/np.hypot(1.65*m,s))**m
    sol=solve_bvp(rhs,ends,s,np.array([guess,np.gradient(guess,s)]),tol=tolerance,max_nodes=25000)
    if not sol.success:
        raise ArithmeticError(sol.message)
    def fields(x):
        x=np.asarray(x)
        g,dg=sol.sol(np.maximum(x,inner))
        f=envelope(x)*g
        df=envelope(x)*(dg+m*g/(np.maximum(x,inner)*(1+x*x)))
        return f,df
    return fields,float(max(sol.rms_residuals))


def stress_check(fields,m):
    # On annuli away from the core s T_ss|_a^b = integral T_phiphi ds.
    q,v=leggauss(64)
    def stress(x):
        f,df=fields(x)
        n=f*f
        ddf=(m*m/x**2+.5*(n**4-1))*f-df/x
        p=(n**5-1)/5
        return p+.5*(df*df-f*ddf),p+n*m*m/x**2-.5*f*df/x
    residuals=[]
    for lo,hi in ((.1,.5),(.5,1.),(1.,2.),(2.,5.),(5.,12.)):
        x=(hi+lo+(hi-lo)*q)/2
        residuals.append(abs(hi*stress(hi)[0]-lo*stress(lo)[0]-(hi-lo)*(v@stress(x)[1])/2))
    return max(residuals)


def held(radius):
    mesh=72
    trial=Trial(radius,mesh,bonded_core=True)
    path=ROOT/f'_scratch/s11c/mouth_exit_release01/r{radius}n72/endpoint.npz'
    shifts,chi,_,digest,status,*_=inputs(path,radius,mesh,trial)
    manifest=ROOT/'research/pde_ledger_v3/_measurements/mouth_breathing_exit_support.json'
    if status!='balanced' or json.loads(manifest.read_text())['physical_preparations'][f'r{radius}n72']!=digest:
        raise ValueError('Input differs from the archived breathing-force preparation')
    sh=Meridian(radius,trial.x+shifts[:4])
    co=trial.co.copy()
    cr=trial.cr.copy()
    co[0]+=shifts[4]
    cr[0]+=shifts[5]
    return trial,sh,chi,co,cr,digest


def supported_order_probes(radius,chi,owner):
    """Named axisymmetric positive reordering tests on the 145 × 73 native mesh."""
    t=np.linspace(0,1,chi.shape[0])[:,None]
    x=np.linspace(0,1,chi.shape[1])[None,:]
    center=(106/144,43/72)
    # A compact C1 test of an extended patch centered on the archived melted
    # R1 node.  The width is fixed in the native material coordinates at both
    # radii, and the (1-chi) factor preserves bounds for a positive probe.
    distance=((t-center[0])/.12)**2+((x-center[1])/.12)**2
    patch=np.maximum(0.,1-distance)**2*(1-chi)
    patch[~owner.free]=0
    point=np.zeros_like(chi)
    point[106,43]=1.
    partial=np.zeros_like(chi)
    partial[61,48]=1.
    return {'archived_site':point,'partial_order_site':partial,'extended_patch':patch}


def material_parts(radius):
    trial,shape,chi,co,cr,digest=held(radius)
    owner=Order(shape,trial.n)
    t=np.linspace(0.,1.,chi.shape[0])
    x=np.linspace(0.,1.,chi.shape[1])
    elastic=source(shape,co,t,x,False)+core_source(shape,co,*trial.core_reference,t,x)
    order=owner.energy_gradient(chi)[1]
    return shape,chi,owner,supported_order_probes(radius,chi,owner),elastic,order


def collar(shape,coeffs,inside):
    # Both principal meridional stretches and both degenerate hoop stretches.
    # Determinant and stretches refer to the unloaded flat crystalline bore.
    a,t=np.meshgrid(np.linspace(0,.5 if inside else -.5,81),np.linspace(0,1,513))
    _,rs,rz,h,j,_=graph_values(a,t,shape,legval(2*t-1,coeffs),2*legval(2*t-1,legder(coeffs)))
    trace=rs*rs+rz*rz+shape.length**2
    deter=(rs*shape.length)**2
    large=np.sqrt((trace+np.sqrt(np.maximum(0,trace**2-4*deter)))/2)
    small=np.sqrt(deter)/large
    return dict(min_J=j.min(),max_J=j.max(),min_stretch=np.minimum(small,h).min(),
                max_stretch=np.maximum(large,h).max())


def direct_angular_work(shape,chi,mesh):
    steady=SteadyMouth(shape,chi,mesh)
    # Brane-polar equator, with density and meridional phase held fixed.
    t=np.linspace(.02,.98,16)[:,None]
    x=np.linspace(.03,.98,15)
    b,slope=shape.radius(t*shape.length)
    r=x*b
    n,ur,uw=steady.base('core',r,t*shape.length)
    # Independent test vectors in the meridian. The cancellation holds for
    # every axisymmetric order profile and every differentiable test direction.
    gr=np.sin(3*t+2*x)/WIDTH
    gw=np.cos(2*t-x)/WIDTH
    dr=np.cos(t-3*x)
    dw=np.sin(t+x)
    h=.0001
    def work(m,shift):
        g1,g2=gr+shift*dr,gw+shift*dw
        a,b=ur,uw
        inner=a*g1+b*g2
        return .5*n*(a*a+b*b+(m/r)**2-ORDER_INERTIA*inner**2/(n+ORDER_INERTIA*(g1*g1+g2*g2)))
    def slope(m):
        return (work(m,h)-work(m,-h))/(2*h)
    # Use the same arbitrary meridional test phase gradient in all three
    # probes. This is a functional identity, not a wound background solution.
    plus,minus,zero=slope(1),slope(-1),slope(0)
    angular=n/(2*r*r)
    return dict(max_signed_angle_energy_parity_error=float(np.max(abs(work(1,0)-work(-1,0)))),
                max_added_density_energy_error=float(np.max(abs((work(1,0)-work(0,0))-angular))),
                max_abs_order_slope_difference=float(np.max(abs(plus-zero))),
                max_abs_reversed_slope_difference=float(np.max(abs(plus-minus))))


def main():
    fields=None
    for m in (1,2):
        for outer,tol in ((32.,1e-7),(64.,1e-9)):
            f,residual=profile(m,outer,tol)
            rhalf=brentq(lambda s:f(s)[0]**2-.5,.001,outer/2)
            sonic=brentq(lambda s:s*f(s)[0]**4-m,.001,outer/2)
            err=stress_check(f,m)
            emit(kind='neutral_control',winding=m,outer=outer,collocation_tol=tol,
                 max_rms_residual=residual,stress_integral_error=err,half_density_radius=rhalf,
                 low_k_sonic_radius=sonic,core_density=float(f(0.)[0]**2))
            if err>2e-7:
                raise ArithmeticError('Native vortex violates independent cylindrical stress balance')
            if m==1 and outer==64: fields=f
    q,v=leggauss(80)
    for radius in (1,2):
        trial,sh,chi,co,cr,digest=held(radius)
        b=sh.exit
        # Ball B³ sliced into cylindrical shells; s=b sin u integrates the
        # regular volume 4π s sqrt(b²-s²) ds without a polar-axis cutoff.
        u=np.pi*(1+q)/4
        x=b*np.sin(u)
        weight=np.pi*v/4*np.sin(u)*np.cos(u)**2*3
        feq=float(fields(b)[0]); neq=feq**2
        volume=np.sum(weight*fields(x)[0]**2)
        shape,chi,owner,directions,elastic,order=material_parts(radius)
        coefficients=direct_angular_work(shape,chi,trial.mesh)
        carrier=Carrier(sh,trial.n,2*trial.n,trial.mesh)
        baseline_grad=carrier.response(chi,BIAS)[1]
        direction_readings={}
        for name,direction in directions.items():
            g0=np.sum(order*direction)
            ge=np.sum(elastic*direction)
            gc=np.sum(baseline_grad*direction)
            direction_readings[name]=dict(order_work=g0,elastic_work=ge,unwound_carrier_work=gc,
                unwound_total=g0+ge+gc,max_direction=float(direction.max()),
                centered_site_order=float(chi[106,43] if name!='partial_order_site' else chi[61,48]))
        emit(kind='held_material',entrance_radius=radius,material_sha256=digest,exit_radius=b,
             neutral_equator_density=neq,neutral_equator_low_k_Mach=1/(b*neq**2),
             neutral_B3_mean_density=volume,neutral_B3_fraction_above_low_k_speed=np.sum(weight*(1/x>fields(x)[0]**4)),
             neutral_shear_order_half_level=.5,registered_inner=collar(sh,co,True),
             registered_outer=collar(sh,cr,False),fixed_state_angular=coefficients,
             baseline_order_directions=direction_readings)
        if abs(weight.sum()-1)>1e-12 or coefficients['max_abs_order_slope_difference']>2e-6:
            raise ArithmeticError('Quadrature or direct angular order identity failed')


if __name__=='__main__':
    main()
