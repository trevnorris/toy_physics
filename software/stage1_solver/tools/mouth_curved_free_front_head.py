"""Stationary radial two-phase front, 4D capillary tension and a supplied head.

Baths and the shear-only inner mirror are fixed; the internal order level is
free. Two collocation intervals share an unknown radius and match both values
and derivatives there. The radiation closure is a radial-ray hypothesis with
E_w per unit inner area = J*pi/(2Y), Y=int h(chi) dr. All three equal angular
stresses are included in the momentum check.
"""
import argparse

import numpy as np
from scipy.integrate import cumulative_trapezoid, simpson, solve_bvp
from scipy.special import expit

from mouth_free_front_wave_head import bath_head, emit, fraction
from mouth_flat_phase_control import saved_bath


def solve(inner, span, load, action, initial_depth, step, tol, phase_target=None, holder=None,
          radial_mobility=None):
    kappa = 2.
    def mobility(r):
        return .5 if radial_mobility is None else radial_mobility(r)
    outer = inner+span
    wave_constant = action*np.pi/2
    n_bulk = (1-5*load)**.2
    tau = np.linspace(0.,1.,max(30,round(span/step/2))+1)

    def area(r):
        return (r/inner)**3

    def holder_affinity(r):
        return 0. if holder is None else holder.affinity(r)

    def local(r,y,j,p):
        f,df,q,dq,head,_,_ = y
        n = f**2
        j_local = j/area(r)
        affinity = -j_local*dq/(mobility(r)*n)
        wave_affinity = -p*6*q*(1-q)/area(r)
        return np.array([df, 2*f*(bath_head(n,j_local)-head) - 3*df/r,
                         dq, (2*q*(1-q)*(1-2*q)+wave_affinity+holder_affinity(r)-affinity)/kappa-3*dq/r,
                         affinity*dq/n, fraction(q), j_local/n])

    def ode(t,y,parameters):
        j,p,mid = parameters[:3]
        left = mid-inner
        right = outer-mid
        return np.vstack([left*local(inner+left*t,y[:7],j,p),
                          right*local(mid+right*t,y[7:],j,p)])

    def bc(a,b,parameters):
        j,p,mid = parameters[:3]
        left,right = a[:7],b[7:]
        far_density = n_bulk if phase_target is None else parameters[3]
        extra = [] if phase_target is None else [right[6]-phase_target]
        return np.r_[
            left[0]-1,left[2]-1,left[4]-bath_head(1,j),left[5],left[6],
            right[0]-np.sqrt(far_density),right[2],
            right[4]-bath_head(far_density,j/area(outer)),p*right[5]**2-wave_constant,
            b[:7]-a[7:],b[2]-.5,extra]

    # This estimate initializes Newton only. Neither depth nor current occurs
    # as prescribed data at the common internal section.
    mid=inner+initial_depth
    j=6*mobility(mid)*load*area(mid)
    r=np.r_[inner+initial_depth*tau,mid+(outer-mid)*tau[1:]]
    q=expit(mid-r)
    f=1+(np.sqrt(n_bulk)-1)*(1-q)
    head=bath_head(1,j)*q+bath_head(n_bulk,j/area(outer))*(1-q)
    guess=np.array([f,0*q,q,-q*(1-q),head,
                    cumulative_trapezoid(fraction(q),r,initial=0),
                    cumulative_trapezoid(j/area(r)/f**2,r,initial=0)])
    guess=np.vstack([guess[:,:len(tau)],guess[:,len(tau)-1:]])
    parameters = [j,wave_constant/initial_depth**2,mid]
    if phase_target is not None:
        parameters.append(n_bulk)
    sol=solve_bvp(ode,bc,tau,guess,p=parameters,
                  tol=tol,max_nodes=28000)
    if not sol.success:
        raise ArithmeticError(sol.message)
    j,p,mid=sol.p[:3]
    far_density = n_bulk if phase_target is None else sol.p[3]
    pressure_drop = (1-far_density**5)/5
    if not inner < mid < outer:
        raise ArithmeticError("Order level outside physical domain")
    t=np.linspace(0,1,4001)
    z=sol.sol(t)
    dz=sol.sol(t,1)
    states=np.hstack([z[:7],z[7:,1:]])
    gradients=np.hstack([dz[:7]/(mid-inner),dz[7:,1:]/(outer-mid)])
    r=np.r_[inner+(mid-inner)*t, mid+(outer-mid)*t[1:]]
    f,df,q,dq,B,width,phase=states
    ddf,ddq=gradients[[1,3]]
    n=f*f
    a=area(r)
    W=(q*(1-q))**2
    net_mu=2*q*(1-q)*(1-2*q)-kappa*(ddq+3*dq/r)-p*6*q*(1-q)/a
    net_mu = net_mu + holder_affinity(r)
    heat=simpson(a*mobility(r)*net_mu**2,x=r)
    bath_power=j*(B[0]-B[-1])
    prr=((j/a)**2/n+n**5/5+(df**2-f*ddf)/2
         +kappa*dq**2/2-W+p*fraction(q)/a)
    ptt=n**5/5-f*df/(2*r)-kappa*dq**2/2-W
    # Radial divergence includes the applied order-conjugate external force.
    holder_force = holder_affinity(r)*dq
    ledger=a*prr-cumulative_trapezoid(3*a*ptt/r+a*holder_force,r,initial=0)
    weighted_end=np.max(abs(kappa*dq[[0,-1]]**2/2))
    bare_end=(W-kappa*dq**2/2)[-1]-(W-kappa*dq**2/2)[0]
    curvature=-simpson(3*kappa*dq**2/r,x=r)
    wave_push=-simpson(p/a*6*q*(1-q)*dq,x=r)
    conversion=-j*simpson(dq**2/(mobility(r)*a*n),x=r)
    holder_projection=simpson(holder_force,x=r)
    emit(test='free_spherical_front' if holder is None else 'held_spherical_front', inner=inner, outer=outer, brane_pressure_excess=pressure_drop, phase_target=phase_target,
         bulk_density=far_density,
         trial_action_per_inner_area=action, seed_depth=initial_depth, interface_radius=mid,
         effective_length=width[-1], current_per_inner_area=j, current_at_front=j/area(mid),
         density_at_front=sol.sol(1)[0]**2, mirror_pressure=p,
         pressure_at_front_from_sharp_ray=p/area(mid), curvature_projection=curvature,
         wave_projection=wave_push, end_order_projection=bare_end, conversion_projection=conversion,
         holder_projection=holder_projection, holder_force_per_inner_area=simpson(a*holder_force,x=r),
         phase_bulk_minus_brane=phase[-1], head_brane_minus_bulk=B[0]-B[-1],
         reservoir_power=bath_power, conversion_loss=heat, radial_momentum_residual=np.ptp(ledger),
         boundary_order_traction_max=weighted_end, seam_error=max(abs(z[:7,-1]-z[7:,0])),
         rms_residual=max(sol.rms_residuals), nodes=len(sol.x))
    if max(abs(z[:7,-1]-z[7:,0]))>5e-9 or min(n)<=0 or max(np.diff(q))>1e-8:
        raise ArithmeticError("Internal matching or phase ordering failed")
    if weighted_end > 2e-5*abs(pressure_drop) or np.ptp(ledger)>2e-8:
        raise ArithmeticError("Boundary stress or radial momentum contamination")
    if abs(bath_power-heat)>2e-10 or abs(bare_end+curvature+wave_push+holder_projection-conversion)>2e-8:
        raise ArithmeticError("Energy or unconstrained translation balance failed")
    if holder is not None:
        # Evaluate both sides at the shared level. The switch removes only the
        # order-conjugate affinity. Density, local velocity and wave load are
        # continuous at that instant under a finite, impulse-free release.
        speeds = []
        for state, derivative in ((z[:7,-1],dz[:7,-1]/(mid-inner)),
                                  (z[7:,0],dz[7:,0]/(outer-mid))):
            amp, _, c, slope = state[:4]
            unheld_mu = (2*c*(1-c)*(1-2*c) - kappa*(derivative[3]+3*slope/mid)
                         - p*6*c*(1-c)/area(mid))
            speeds.append(j/(area(mid)*amp**2)+mobility(mid)*unheld_mu/slope)
        speed_from_removed_force = -mobility(mid)*holder_affinity(mid)/z[3,-1]
        offset = mid-holder.center
        emit(test='initial_release', action=action, strength=holder.strength,
             holder_center=holder.center, width=holder.width, radius=mid, offset=offset,
             carrier_speed=j/(area(mid)*z[0,-1]**2), left_speed=speeds[0],
             right_speed=speeds[1], speed_from_removed_force=speed_from_removed_force,
             order_slope=z[3,-1], offset_times_speed=offset*speeds[0])
        if max(abs(np.array(speeds)-speed_from_removed_force))>1e-10:
            raise ArithmeticError("Two-sided release velocity and removed force disagree")
        return offset, speeds[0]
    if radial_mobility is not None:
        integrated_head = j*simpson(dq**2/(a*n*n*mobility(r)), x=r)
        emit(test='radial_mobility_balance', minimum_mobility=np.min(mobility(r)),
             front_mobility=mobility(mid), maximum_mobility=np.max(mobility(r)),
             independently_integrated_head=integrated_head, bath_head=B[0]-B[-1])
        if abs(integrated_head-(B[0]-B[-1]))>2e-9:
            raise ArithmeticError("Integrated diffuse-layer resistance disagrees with the bath head")


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--spacing',type=float,default=.3)
    p.add_argument('--tolerance',type=float,default=1e-8)
    p.add_argument('--seed-depth',type=float,default=17.)
    p.add_argument('--lower-head-seed',type=float,default=24.)
    p.add_argument('--cases',choices=('all','lower-head'),default='all')
    p.add_argument('--phase-match',action='store_true')
    args=p.parse_args()
    span = 32. if args.phase_match else 50.
    emit(test='preparation', space_dimensions=4, radius_unit='delta_B', kappa=2.,a_B=1.,
         mobility=.5, inner=200., outer=200.+span, tension=1/3,
         control='phase bias with head supplied as solved' if args.phase_match else 'pressure control',
         wave_energy_per_inner_area='J*pi/(2*integral h(chi)dr)', h='chi^2(3-2chi)',
         order_anchor='none', wave_mirror='fixed inner reservoir only')
    if args.phase_match:
        for action in (0.,.001):
            solve(200.,span,.0045,action,args.seed_depth,args.spacing,args.tolerance,saved_bath())
        return
    if args.cases=='all':
        for action in (0.,.001,.002):
            solve(200.,50.,.0045,action,args.seed_depth,args.spacing,args.tolerance)
    solve(200.,50.,.00435,.001,args.lower_head_seed,args.spacing,args.tolerance)


if __name__=='__main__':
    main()
