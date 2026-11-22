"""Free two-sided planar order front under a conservative shear-cavity load.

Use a thin, perfectly reflecting shear front for the wave energy only. The
carrier and diffuse order field share the whole interval with no internal
boundary or imposed order level. A fixed shear mirror at x=0 takes the opposite
wave load. The coherent carrier can cross that mirror plane. A steady solution
is driven by a specified mechanical head from two prepared baths.
"""

import argparse
import json

import numpy as np
from scipy.integrate import cumulative_trapezoid, simpson, solve_bvp
from scipy.optimize import brentq
from scipy.special import expit


def emit(**values):
    print(json.dumps(values, allow_nan=False), flush=True)


def fraction(q):
    return q*q*(3.0 - 2.0*q)


def bath_head(n, j):
    return (n**4 - 1.0)/4.0 + .5*(j/n)**2


def solve(length, load, action, seed_factor, tol, spacing):
    # a_B=delta_B=n0=m=hbar=c_T=1, kappa_B=2, M=1/2.
    # The clamped/free fundamental shear wave has E=J pi/(2Y);
    # Y is the effective length of ordered material.
    mobility = .5
    n_bulk = (1.0+5.0*load)**.2
    wave_constant = action*np.pi/2
    x = np.linspace(0, length, round(length/spacing)+1)

    def equations(_, states, parameters):
        f, df, q, dq, head, _, _ = states
        j, pressure = parameters
        n = f*f
        bare = 2*q*(1-q)*(1-2*q)
        wave_affinity = -pressure*6*q*(1-q)
        affinity = -j*dq/(mobility*n)
        return np.asarray([
            df, 2*f*(bath_head(n,j)-head),
            dq, .5*(bare+wave_affinity-affinity),
            affinity*dq/n, fraction(q), j/n+0*q,
        ])

    def baths(left, right, parameters):
        j, pressure = parameters
        return np.asarray([
            left[0]-1, right[0]-np.sqrt(n_bulk), left[2]-1, right[2],
            left[4]-bath_head(1,j), right[4]-bath_head(n_bulk,j),
            left[5], right[5]**2*pressure-wave_constant, left[6],
        ])

    guess_depth = seed_factor*np.sqrt(wave_constant/load)
    q = expit(guess_depth-x)
    f = 1+(np.sqrt(n_bulk)-1)*(1-q)
    initial = np.asarray([f,0*q,q,-q*(1-q),(1-q)*bath_head(n_bulk,0),
                          cumulative_trapezoid(fraction(q),x,initial=0),0*q])
    sol = solve_bvp(equations,baths,x,initial,p=[0.,load],
                    tol=tol,max_nodes=20000)
    if not sol.success:
        raise ArithmeticError(sol.message)
    grid = np.linspace(0,length,6001)
    f, df, q, dq, head, width, phase = sol.sol(grid)
    ddf, ddq = sol.sol(grid,1)[[1,3]]
    j, pressure = sol.p
    n=f**2
    half=brentq(lambda y:sol.sol(y)[2]-.5,0,length)
    # Compute affinity from the energy, and compare its power with the baths.
    mu = 2*q*(1-q)*(1-2*q)-2*ddq-pressure*6*q*(1-q)
    friction = simpson(mobility*mu**2,x=grid)
    bath_work = j*(head[0]-head[-1])
    order_boundary = (dq**2-(q*(1-q))**2)[[0,-1]]
    # For this wave approximation the radiation stress changes from P to zero
    # across the front. The fixed mirror feels -P and does no mechanical work.
    flux = j*j/n + n**5/5 + (df**2-f*ddf)/2 + dq**2-(q*(1-q))**2 + pressure*fraction(q)
    phase_only_head_needed = (j/mobility)*simpson((dq/n)**2,x=grid)
    row = dict(test='head_selected', length=length, bulk_pressure_excess=load,
               action=action, seed_factor=seed_factor, nodes=len(sol.x),
               effective_length=width[-1], half_order=half, relative_current=j,
               incident_side='bulk' if j<0 else 'brane', fixed_plane_current=j,
               wave_force_on_front=pressure, mirror_force=-pressure,
               carrier_phase_bulk_minus_brane=phase[-1], head_brane_minus_bulk=head[0]-head[-1],
               integrated_head_requirement=phase_only_head_needed,
               bath_work=bath_work, conversion_loss=friction,
               total_normal_flux_range=float(np.ptp(flux)), density_min=n.min(),
               boundary_order_traction_max=float(max(abs(order_boundary))),
               maximum_bvp_residual=max(sol.rms_residuals),
               phase_integral=j*simpson(1/n,x=grid),
               max_order_equation_error=max(abs(j*dq/n+mobility*mu)))
    emit(**row)
    if row['boundary_order_traction_max']>1e-5*load or np.ptp(flux)>3e-8:
        raise ArithmeticError("End order tractions or momentum imbalance contaminate front")
    if abs(bath_work-friction)>1.5e-12 or abs(phase[-1]-row['phase_integral'])>2e-8:
        raise ArithmeticError("Reservoir work, conversion or carrier phase failed")
    if abs(row['head_brane_minus_bulk']-phase_only_head_needed)>3e-10:
        raise ArithmeticError("Head drop derived from conversion disagrees with reservoirs")


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--length',type=float,default=46.)
    p.add_argument('--seed-factor',type=float,default=1.)
    p.add_argument('--tolerance',type=float,default=1e-9)
    p.add_argument('--spacing',type=float,default=.25)
    a=p.parse_args()
    emit(test='preparation', mobility=.5,a_B=1.,delta_B=1.,shear_speed=1.,
         wave_energy='J*pi/(2*Y)', Y='integral chi^2(3-2chi) dx',
         pin='none', baths='remote pure phases, homogeneous specific heads, bulk pressure excess')
    for action,load in ((.02,.00010),(.01,.00010),(.02,.00005)):
        solve(a.length,load,action,a.seed_factor,a.tolerance,a.spacing)


if __name__=='__main__':
    main()
