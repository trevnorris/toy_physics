"""Free turning neck and flare of a 4D membrane with axial wave work.

Both end radii and the end separation vary. Pressure acts on the occupied
volume and a wave of fixed action has energy alpha/L. This sharp-interface
test omits side loading by the wave and nonuniform carrier stresses.
"""
import json
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.integrate import solve_ivp

SIGMA = 1.
EXIT_RATIO = 1.5


def profile(exit_ratio):
    # b is measured in waist radii, w in the same units. p is p_phys R/sigma.
    p = 3*(exit_ratio+1)/(exit_ratio**2+exit_ratio+1)
    tension_flux = 1-p/3

    def derivative(_, z):
        b, slope = z
        return [slope, 2*(1+slope**2)/b - p*(1+slope**2)**1.5]

    def end(_, z):
        return z[1]
    end.direction = -1.
    end.terminal = True
    sol = solve_ivp(derivative, (0., 12*exit_ratio), [1., 0.], events=end,
                    method='DOP853', dense_output=True, max_step=.06,
                    rtol=2e-12, atol=2e-14)
    if not sol.success or len(sol.t_events[0]) != 1:
        raise ArithmeticError('Missing flared turning section')
    length = sol.t_events[0][0]
    np.testing.assert_allclose(sol.sol(length), [exit_ratio,0.], atol=1e-10, rtol=0)

    # Independent first-integral quadrature; sine substitution resolves the two turning points.
    root, weight = leggauss(80)
    angle = (root+1)*np.pi/4
    b = 1+(exit_ratio-1)*np.sin(angle)**2
    cosine = tension_flux/b**2 + p*b/3
    length_check = np.pi/4 * weight @ ((exit_ratio-1)*np.sin(2*angle)*cosine/np.sqrt(1-cosine**2))
    np.testing.assert_allclose(length_check, length, rtol=1e-8)
    return p, tension_flux, length, length_check, sol.sol


def measure(radius, dimensionless):
    p, tension_flux, unit_length, integral_length, solution = dimensionless
    load = p*SIGMA/radius
    length = radius*unit_length
    axial_force = 4*np.pi*SIGMA*radius**2*tension_flux
    alpha = axial_force*length**2
    x, weight = leggauss(100)
    dw = length*weight/2
    b, slope = solution(unit_length*(1+x)/2)
    b = b*radius
    arc = np.hypot(1.,slope)
    k_angle = 1/(b*arc)
    k_meridian = load/SIGMA-2*k_angle

    def potential(normal_step=0., stretch=1.):
        # Parallel surface obtained by displacing along the physical unit normal.
        new_b = b+normal_step/arc
        axial_jacobian = 1+normal_step*k_meridian
        if min(axial_jacobian)<=0 or min(new_b)<=0:
            raise ArithmeticError('Normal displacement crossed a focal point')
        area = dw @ (4*np.pi*new_b**2*axial_jacobian*np.hypot(stretch,slope))
        volume = dw @ ((4*np.pi/3)*new_b**3*stretch*axial_jacobian)
        return SIGMA*area-load*volume+alpha/(length*stretch)

    base = potential()
    shape_curvature = -SIGMA*dw @ (4*np.pi*b**2*arc*(2*k_angle**2+k_meridian**2))
    rows = []
    for fraction in (.01, .005):
        eps = radius*fraction
        minus, plus = potential(-eps), potential(eps)
        curvature = (plus+minus-2*base)/eps**2
        # The parallel-surface volume is quartic, so the central second difference has O(eps²) error.
        np.testing.assert_allclose(curvature,shape_curvature,rtol=2e-4)
        if not (plus < base and minus < base):
            raise ArithmeticError('Both signed free radial perturbations must be checked')
        rows.append(dict(fraction=fraction, inward_energy_change=minus-base,
                         outward_energy_change=plus-base, second_difference=curvature))
    eps = .0001
    axial_derivative = (potential(stretch=1+eps)-potential(stretch=1-eps))/(2*eps)
    if abs(axial_derivative)>2e-6*max(1.,abs(base)):
        raise ArithmeticError('Axial wave thrust failed independent stretch-work check')
    return dict(waist=radius,flare=EXIT_RATIO*radius,free_length=length,
                length_from_first_integral=radius*integral_length,
                sigma=SIGMA,required_uniform_pressure=load,required_axial_wave_constant=alpha,
                wave_thrust=axial_force,axial_stretch_energy_derivative=axial_derivative,
                normal_stiffness_from_curvature=shape_curvature,
                normal_displacement_controls=rows)


def main():
    dimensionless = profile(EXIT_RATIO)
    training = measure(1., dimensionless)
    diagnostic = measure(2., dimensionless)
    pressure = training['required_uniform_pressure']
    # At a free orthogonal waist b'=0, b''=2/R-p/sigma. A strict
    # interior flaring minimum requires positive curvature in this convention.
    frozen = [dict(waist=r, waist_second_derivative=2/r-pressure/SIGMA)
              for r in (1.,2.)]
    print(json.dumps(dict(equal_pressure_turning_end_condition='b_exit²=b_waist²',
                          fitted_reference=training, frozen_pressure_and_action=frozen,
                          retuned_for_diagnosis_only=diagnostic),indent=2,allow_nan=False))


if __name__=='__main__':
    main()
