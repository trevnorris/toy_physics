"""Conversion mobility of supported planar and free radial order fronts at fixed bath pressures.

b is the Euclidean norm of the three principal curvatures times the neutral
interface width. Each regular order sphere has b=sqrt(3)/r in these units.
All other material coefficients and the two pressure preparations are fixed.
"""
import argparse

import numpy as np
from scipy.integrate import simpson, solve_bvp
from scipy.special import expit

from mouth_curved_free_front_head import solve
from mouth_free_front_wave_head import bath_head, emit


def mobile_fraction(bend, onset):
    t = np.clip(bend/onset - 1., 0., 1.)
    return t*t*(3.-2.*t)


def locked_flat(pressure_drop, spacing):
    # A material face at rest has j=0. A compactly concentrated conservative
    # order load supplies its pressure support. No division by M occurs.
    n_bulk = (1.-5.*pressure_drop)**.2
    exterior = 25.
    grid = np.linspace(-exterior, exterior, round(2*exterior/spacing)+1)

    def shape(s):
        return expit(-s)*expit(s)

    def ode(s, y, p):
        f, df, head = y
        return np.array([df, 2*f*(bath_head(f*f, 0.)-head), p[0]*shape(s)**2/(f*f)])

    def bc(left, right, _):
        return np.array([left[0]-1., right[0]-np.sqrt(n_bulk),
                         left[2], right[2]-bath_head(n_bulk, 0.)])

    guess = np.array([1.+(np.sqrt(n_bulk)-1.)*expit(grid), grid*0,
                      bath_head(n_bulk,0.)*expit(grid)])
    sol = solve_bvp(ode, bc, grid, guess, p=[-6*pressure_drop], tol=2e-9)
    if not sol.success:
        raise ArithmeticError(sol.message)
    s = np.linspace(-exterior, exterior, 8001)
    f, df, head = sol.sol(s)
    ddf = sol.sol(s, 1)[1]
    stress = f**10/5 + (df**2-f*ddf)/2  # kink has zero net normal order stress
    applied_force = sol.p[0]*simpson(shape(s)**2, x=s)
    emit(test='locked_planar_face', pressure_brane_minus_bulk=pressure_drop,
         head_brane_minus_bulk=head[0]-head[-1], number_current_relative=0.,
         per_constituent_work=None, reservoir_power=0., fixture_force=applied_force,
         stress_difference=stress[-1]-stress[0], no_fixture_momentum_residual=stress[-1]-stress[0],
         chemical_affinity_coefficient=sol.p[0], equal_well_energy_difference=0.,
         rms_residual=max(sol.rms_residuals))
    if abs(applied_force-(stress[-1]-stress[0])) > 2e-9 or np.min(f)<=0:
        raise ArithmeticError('Locked face traction and positive-density checks failed')


def radial(onset, action, step, tol, seed_depth):
    inner, span = 200., 50.
    def mobility(radius):
        return .5*mobile_fraction(np.sqrt(3.)/radius, onset)
    floor, ceiling = mobility(inner+span), mobility(inner)
    emit(test='radial_law', onset=onset, action=action, inner=inner, outer=inner+span,
         brane_pressure_excess=.0045, face_mobility=.5*mobile_fraction(0.,onset),
         rim_side_mobility=ceiling, bulk_side_mobility=floor)
    if ceiling == 0:
        # j chi'/n = -M mu, and any profile connecting distinct pure phases
        # has at least one nonzero derivative. Constant radial flux must vanish.
        emit(test='fully_locked_radial_law', current_allowed_by_order=0.,
             pure_phase_order_change=1., maximum_mobility=ceiling,
             static_stress_balance_solved=False)
    elif floor == 0:
        emit(test='partly_locked_tail', active_up_to_radius=np.sqrt(3.)/onset,
             stationary_rate=None, reason='A degenerate free-boundary solve is needed to decide whether the order is constant throughout the locked region.')
    elif floor == .5:
        emit(test='uniform_operator', same_as_archived=True, minimum_mobility=floor)
    else:
        solve(inner,span,.0045,action,seed_depth,step,tol,radial_mobility=mobility)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--onset', type=float, default=.006)
    parser.add_argument('--action', type=float, default=.001)
    parser.add_argument('--spacing', type=float, default=.22)
    parser.add_argument('--tol', type=float, default=2e-9)
    parser.add_argument('--seed-depth', type=float, default=20.)
    parser.add_argument('--flat', action='store_true')
    args=parser.parse_args()
    if args.flat:
        for load in (.0045, -.0045, 0.):
            locked_flat(load,args.spacing)
    else:
        if args.onset<=0:
            raise ValueError('A positive dimensionless curvature onset is required')
        radial(args.onset,args.action,args.spacing,args.tol,args.seed_depth)


if __name__ == '__main__':
    main()
