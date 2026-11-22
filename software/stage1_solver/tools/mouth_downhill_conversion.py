"""Test a per-constituent downhill order energy on a coherent two-sided front.

A planar front has an analytic selected relative current. In four dimensions
solve the freely positioned spherical order front between two reservoirs,
using the planar fields only for an initial Newton iterate. The bath densities
and homogeneous specific energies are fixed; the phase drop and flow follow.
"""
import argparse
import math

import numpy as np
from scipy.integrate import cumulative_simpson, cumulative_trapezoid, simpson, solve_bvp
from scipy.optimize import brentq
from scipy.special import expit

from mouth_free_front_wave_head import bath_head, emit, fraction
from mouth_curved_front_release import LocalHolder


KAPPA = 2.
MOBILITY = .5
EPSILON = .015
OUTER = 400.


def front(inner, outer, eps, guess, interval, tol, holder=None, planar=False):
    # q=1 in the feed, q=0 in the bulk. In a planar cell fix translation by
    # reporting the representative with the half-order level at the origin.
    width = outer-inner
    grid = np.linspace(0, 1, math.ceil(width/(2*interval))+1)
    def scale(r):
        return 1.+0.*r if planar else (r/inner)**3

    def gradients(r, fields, flux):
        f, df, q, dq, Bc, _ = fields  # Bc excludes the new specific energy
        n = f*f
        j = flux/scale(r)
        curve = 0. if planar else 3/r
        affinity = -j*dq/(MOBILITY*n)
        g = 6*q*(1-q)
        potential = 2*q*(1-q)*(1-2*q)
        applied = 0. if holder is None else holder.affinity(r)
        return np.array([df,2*f*(bath_head(n,j)-Bc)-curve*df, dq,
                         (potential+eps*n*g+applied-affinity)/KAPPA-curve*dq,
                         (affinity/n-eps*g)*dq,j/n])

    def ode(t, y, params):
        flux, mid = params
        a, b = mid-inner, outer-mid
        return np.vstack([a*gradients(inner+a*t, y[:6], flux),
                          b*gradients(mid+b*t, y[6:], flux)])

    def ends(a, b, params):
        flux, mid = params
        left, right = a[:6], b[6:]
        # Pure remote ends have bulk head defined thermodynamically. Planar
        # translation is gauged, so drop the redundant right head condition
        # and check it after solving instead.
        right_head = (mid-guess if planar else
                      right[4]-bath_head(1.,flux/scale(outer)))
        return np.r_[left[0]-1, left[2]-1, left[4]-bath_head(1.,flux),
                     left[5], right[0]-1, right[2], right_head,
                     b[:6]-a[6:], b[2]-.5]

    j0 = 3*eps if planar else scale(guess)*(3*eps+3/guess)
    r = np.r_[inner+(guess-inner)*grid,guess+(outer-guess)*grid[1:]]
    q = expit(guess-r)
    v = j0/scale(r)
    B0 = bath_head(1.,j0)*q + bath_head(1.,j0/scale(outer))*(1-q)
    f = np.ones_like(q)
    start=np.array([f,0*f,q,-q*(1-q),B0,cumulative_trapezoid(v,r,initial=0)])
    result=solve_bvp(ode,ends,grid,np.vstack([start[:,:len(grid)],start[:,len(grid)-1:]]),
                     p=[j0,guess],tol=tol,max_nodes=22000)
    if not result.success:
        raise ArithmeticError(result.message)
    flux,mid=result.p
    if not inner<mid<outer:
        raise ArithmeticError('Half-order surface outside two pure phases')
    t=np.linspace(0,1,5001)
    z=result.sol(t); z1=result.sol(t,1)
    f,df,q,dq,Bc,phase=np.c_[z[:6],z[6:,1:]]
    dr=np.c_[z1[:6]/(mid-inner),z1[6:,1:]/(outer-mid)]
    r=np.r_[inner+(mid-inner)*t,mid+(outer-mid)*t[1:]]
    area=scale(r); n=f*f; j=flux/area; H=fraction(q); g=6*q*(1-q)
    curve = 0. if planar else 3/r
    applied=0.*r if holder is None else holder.affinity(r)
    mu=2*q*(1-q)*(1-2*q)+eps*n*g-KAPPA*(dr[3]+curve*dq)+applied
    loss=simpson(area*MOBILITY*mu**2,x=r)
    bk=bath_head(1.,flux)-bath_head(1.,flux/scale(outer))
    fuel_power=eps*flux
    work=flux*(Bc[0]-Bc[-1])+fuel_power
    local_order_error=max(abs(MOBILITY*mu+j*dq/n))
    bare=(q*(1-q))**2
    kin=j*j/n; pres=n**5/5
    normal=kin+pres+(df**2-f*dr[1])/2+KAPPA*dq*dq/2-bare
    tangent=pres-KAPPA*dq*dq/2-bare
    if not planar:
        tangent-=f*df/(2*r)
    # Integrate the angular stress independently on each smooth chart.
    # Trapezoid integration of the rest pressure can dominate the tiny
    # response; Simpson quadrature is exact on the quadratic radial weight.
    force_density=area*(curve*tangent+applied*dq)
    ns=len(t)-1
    left_sum=cumulative_simpson(force_density[:ns+1],x=r[:ns+1],initial=0.)
    right_sum=cumulative_simpson(force_density[ns:],x=r[ns:],initial=0.)
    ledger=area*normal-np.r_[left_sum,left_sum[-1]+right_sum[1:]]
    end_traction=KAPPA*max(dq[[0,-1]]**2)/2
    out=dict(test='planar' if planar else 'annular',epsilon=eps,inner=inner,outer=outer,
             guess=guess,half_order=mid,ordered_current=flux,front_current=flux/scale(mid),
             front_density=z[0,-1]**2,chemical_drop=eps,end_kinetic_drop=bk,
             total_drop=eps+Bc[0]-Bc[-1],phase_difference=phase[-1]-phase[0],
             fuel_power=fuel_power,net_bath_work=work,conversion_loss=loss,
             max_total_normal_balance=np.ptp(ledger),end_order_stress=end_traction,
             right_head_error=Bc[-1]-bath_head(1.,flux/scale(outer)),
             local_order_error=local_order_error,external_holder=holder is not None,
             mesh_nodes=len(result.x),residual=max(result.rms_residuals))
    if holder is not None:
        # Only the holder affinity is switched off at t=0. No carrier impulse.
        unheld=mu-applied
        k=len(t)-1
        speed=j[k]/n[k]+MOBILITY*unheld[k]/dq[k]
        out.update(probe_offset=mid-holder.center,velocity_on_release=speed,
                   removed_affinity_speed=-MOBILITY*applied[k]/dq[k])
    emit(**out)
    if (max(abs(z[:6,-1]-z[6:,0]))>2e-8 or min(n)<=0 or max(np.diff(q))>1e-7 or
            end_traction>2e-8 or np.ptp(ledger)>2e-7 or abs(work-loss)>4e-9 or
            abs(out['right_head_error'])>1e-8 or local_order_error>2e-7):
        raise ArithmeticError('Phase, energy, local order or total momentum validation failed')
    return mid, out


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('case',choices=('planar','small','large'))
    parser.add_argument('--refined',action='store_true')
    args=parser.parse_args()
    if args.case=='planar':
        for eps in (0., EPSILON, -EPSILON):
            _,row=front(-24.,24.,eps,0.,.32,5e-9,planar=True)
            if abs(row['front_current']-6*MOBILITY*eps)>2e-7:
                raise ArithmeticError('Resolved relative current fails analytic value')
        return
    inner=100. if args.case=='small' else 200.
    # Latitude-free thin-front estimate is used strictly to initialize the
    # root solve. The diffuse curved carrier and BVP set the final value.
    def estimate(r):
        flux=(r/inner)**3*(3*EPSILON+3/r)
        return .5*flux**2*(1-(inner/OUTER)**6)-1/r
    guess=brentq(estimate,inner,OUTER)
    spacing,tol=(.3,1e-8) if args.refined else (.55,4e-8)
    mid,_=front(inner,OUTER,EPSILON,guess,spacing,tol)
    if not args.refined:
        for strength in (-1e-6,1e-6,-5e-7,5e-7):
            front(inner,OUTER,EPSILON,mid,spacing,tol,LocalHolder(mid,strength))


if __name__=='__main__':
    main()
