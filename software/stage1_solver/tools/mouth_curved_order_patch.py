"""Held neutral half-order surfaces on a coherent material patch in four dimensions.

Units: length δ_B=sqrt(kappa_B/(2*a_B)), energy density a_B. The homogeneous
carrier is at rest at unit density on both sides, with no shear displacement
or wave and the minimal additive carrier/order energy. Moving the pin in the
energy check changes the interface location while both bath surfaces stay put.
"""
import argparse
import json

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.integrate import solve_bvp
from scipy.special import expit


KAPPA = 2.
TRAIN = 1/50
HELD_OUT = 1/100


def well(chi):
    return (chi*(1-chi))**2


def emit(**row):
    print(json.dumps(row, allow_nan=False), flush=True)


def half_profile(left, right, ends, pin, curvature, spacing, tol):
    """Equilibrium in the physical annulus, with fixed endpoint order values."""
    n = int(np.ceil((right-left)/spacing))
    nodes = np.linspace(left, right, n+1)
    guess = expit(pin-nodes)

    def equation(s, y):
        chi, slope = y
        local_curvature = 3*curvature/(1+curvature*s)
        return np.vstack([slope, chi*(1-chi)*(1-2*chi)-local_curvature*slope])

    def boundary(lo, hi):
        return np.array([lo[0]-ends[0], hi[0]-ends[1]])

    solution = solve_bvp(equation, boundary, nodes, np.vstack([guess,-guess*(1-guess)]),
                         tol=tol, max_nodes=10000)
    if not solution.success:
        raise ArithmeticError(solution.message)
    return solution


def integrate_half(solution, curvature, rule):
    # Gauss integration over each adaptive collocation interval, normalized
    # to unit physical area at the original central surface.
    x, weights = rule
    edges = solution.x
    half_width = np.diff(edges)[:, None]/2
    s = edges[:-1,None]+half_width*(1+x)
    chi, grad = solution.sol(s.ravel()).reshape(2,*s.shape)
    if np.any(grad>2e-10) or np.min(chi)<-2e-10 or np.max(chi)>1+2e-10:
        raise ArithmeticError('Order profile left the monotone connecting branch')
    measure = weights*half_width
    jac = (1+curvature*s)**3
    stress_anisotropy = KAPPA*grad**2
    energy = np.sum(measure*jac*(well(chi)+stress_anisotropy/2))
    radial_force = np.sum(measure*3*curvature/(1+curvature*s)*stress_anisotropy)
    return energy, np.sum(measure*stress_anisotropy), radial_force


def relax(curvature, extent, pin, spacing, tol):
    if abs(curvature)*extent >= 1:
        raise ValueError('A radial chart reached its symmetry center')
    lo = half_profile(-extent,pin,(1.,.5),pin,curvature,spacing,tol)
    hi = half_profile(pin,extent,(.5,0.),pin,curvature,spacing,tol)
    left_slope, right_slope = lo.sol(pin)[1], hi.sol(pin)[1]
    lo_stress = well(.5)-.5*KAPPA*left_slope**2
    hi_stress = well(.5)-.5*KAPPA*right_slope**2
    far_stresses = [-.5*KAPPA*lo.sol(-extent)[1]**2, -.5*KAPPA*hi.sol(extent)[1]**2]
    rule = leggauss(8)
    energy, anisotropy, integrated = np.add(integrate_half(lo,curvature,rule),
                                            integrate_half(hi,curvature,rule))
    return dict(energy=energy, anisotropy=anisotropy, integrated_force=integrated,
                far_stresses=far_stresses, left_stress=lo_stress, right_stress=hi_stress,
                traction=lo_stress-hi_stress, left_order_slope=left_slope,
                right_order_slope=right_slope, nodes=[len(lo.x),len(hi.x)],
                max_rms_residual=max(max(lo.rms_residuals),max(hi.rms_residuals)))


def run(args):
    emit(test='preparation', carrier_density=1, carrier_phase=0, carrier_current=0,
         shear_displacement=0, wave_action=0, carrier_pressure_jump=0,
         bare_well_depth_jump=0, f_mix=0, f_throat=0, static_carrier_equation_residual=(abs(1.)**8-1.)/4,
         end_conditions='fixed orders one and zero; endpoint slopes measured', delta_B=1, a_B=1, kappa_B=KAPPA, orientation='ordered to disordered',
         extent=args.extent, spacing=args.spacing,tolerance=args.tolerance,
         training_K=3*TRAIN, held_out_K=3*HELD_OUT)
    pressure = {}
    energy_step = .004
    for c in (0., TRAIN, -TRAIN, HELD_OUT, -HELD_OUT):
        center = relax(c,args.extent,0.,args.spacing,args.tolerance)
        derivatives = []
        for step in (energy_step, energy_step/2):
            minus = relax(c,args.extent,-step,args.spacing,args.tolerance)['energy']
            plus = relax(c,args.extent,+step,args.spacing,args.tolerance)['energy']
            derivatives.append(dict(step=step, minus_energy=minus, plus_energy=plus,
                                    energy_derivative=(plus-minus)/(2*step)))
        boundary_jump = center['far_stresses'][1]-center['far_stresses'][0]
        stress_integral = center['integrated_force']-boundary_jump
        emit(test='patch', curvature=3*c, **center, stress_integral=stress_integral,
             displaced_pins=derivatives)
        np.testing.assert_allclose(center['traction'],stress_integral,rtol=7e-6,atol=5e-10)
        np.testing.assert_allclose(center['traction'],derivatives[-1]['energy_derivative'],
                                   rtol=7e-6,atol=5e-10)
        pressure[c] = center['traction']
        if c==0:
            planar_energy = center['energy']
            planar_anisotropy = center['anisotropy']
            np.testing.assert_allclose(planar_anisotropy,planar_energy,rtol=1e-8)
    coefficient = (pressure[TRAIN]-pressure[-TRAIN])/(6*TRAIN)
    for c in (HELD_OUT, -HELD_OUT):
        predicted = pressure[0.]+3*c*coefficient
        emit(test='prediction', curvature=3*c, planar_energy=planar_energy,
             planar_anisotropy=planar_anisotropy, tension_signed_training=coefficient,
             planar_prediction=3*c*planar_energy, signed_pair_prediction=predicted,
             measured=pressure[c], fractional_error=(predicted-pressure[c])/pressure[c])


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--extent', type=float, required=True)
    parser.add_argument('--spacing', type=float, required=True)
    parser.add_argument('--tolerance', type=float, required=True)
    run(parser.parse_args())
