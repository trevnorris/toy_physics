"""Reciprocal held-mouth test of material registration in the two-state contact.

Impose a calibrated dipolar angular displacement y=Y sin(theta) on one
isolated, held cut. Compare the same exterior preparation to carried sites.
Only the difference in the external force needed to keep the center fixed is
computed. The angular waveform is in phase with peak displacement at t=0.
"""
import numpy as np
from numpy.polynomial.legendre import leggauss

from mouth_pair_constituent_response import BOND, EXCHANGE, PIN, analytic_site, site


def virtual_work(radius, eta, carry, ratio, amplitude):
    """Sum the change of microscopic pin and bond energies on displaced tags.

    Spin and exterior state are held instantaneous during each virtual move.
    The same exterior terms cancel against a carried-registration control.
    Taking the part odd in drive removes the undriven energy of the dimer.
    """
    u, weights = leggauss(64)
    sin = np.sqrt(1-u**2)
    _, gap = analytic_site(EXCHANGE)
    pol, chi, _ = site(EXCHANGE, eta, ratio*gap)
    angular_compliance = chi.real/(2*radius*pol)**2
    torque = BOND*amplitude/(1+BOND*angular_compliance)
    spin_angle = angular_compliance*torque
    staggered = 2*radius*pol*spin_angle*sin
    def energy(displacement):
        # Original site R*n travels carry*displacement, the marked rim by displacement.
        shift = (carry-1)*displacement/radius
        length = np.hypot(sin, u+shift)
        normal_z = (u+shift)/length
        sin_turn = -shift*sin/length
        # Difference of actual signed solid angles of the two latitude loops.
        preferred = eta*(u-normal_z)/(2*radius*sin)
        pin_odd = -eta*PIN*sin_turn*staggered
        bond_odd = BOND*(amplitude-spin_angle)*sin*preferred
        return 2*np.pi*radius**2*(weights@(pin_odd+bond_odd))
    step = 1e-4*radius
    force = (8*(energy(step)-energy(-step))-energy(2*step)+energy(-2*step))/(12*step)
    return force, torque, spin_angle


def main():
    for radius in (5., 10.):
        pol, _ = analytic_site(EXCHANGE)
        h = 4*radius**2*pol*PIN
        for carry in (0., .4, 1.):
            force, torque, slip = virtual_work(radius, 1, carry, .65, 1e-4)
            # Invert measurements only from this positive, held preparation.
            unit_force = (4*np.pi/3)*(torque-h*slip)
            lag = force/unit_force
            print(f"held R={radius:g} carry={carry:g} plus_force_per_Y={force/1e-4:.10g} "
                  f"inferred_lag={lag:.10g}")
            assert abs(lag-(1-carry))<1e-8
            for ratio, eta, amplitude in ((.30,-1,5e-5), (1.15,-1,-5e-5), (1.15,1,1e-4)):
                observed, load, tilt = virtual_work(radius, eta, carry, ratio, amplitude)
                prediction = eta*lag*(4*np.pi/3)*(load-h*tilt)
                error = abs(observed-prediction)/abs(amplitude)
                print(f"held control R={radius:g} carry={carry:g} gap_ratio={ratio:g} eta={eta:+d} "
                      f"Y={amplitude:g} observed_F_over_Y={observed/amplitude:.10g} "
                      f"predicted_F_over_Y={prediction/amplitude:.10g} error_per_Y={error:.3g}")
                assert error<1e-8


if __name__ == "__main__":
    main()
