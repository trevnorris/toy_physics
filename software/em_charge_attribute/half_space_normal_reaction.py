"""Excess pressure force on the intact wall of a stationary ideal half-space."""

import numpy as np
import sympy as sp


def outside_force(peak_speed, terms):
    # Angular integration of the 4D Neumann kernel gives the logarithm
    # log((r+y)/(r-y)) on the intact face r>1. Expand it there. The kth
    # term of d(phi)/dr is coefficient[k] / r**(2*k+3), with positive terms.
    k = np.arange(terms, dtype=float)
    power = 2*k
    moment = 48/((power+3)*(power+5)*(power+7)*(power+9))
    coefficient = 2*peak_speed/np.pi * (power+2)/(power+1)*moment
    # In the receiver, Bernoulli gives p-p_inf = -rho * v_t**2 / 2
    # on the intact plane. Its force on the plane points toward +w.
    integrated_square = sum(c*np.sum(coefficient/(p+power+3))
                            for c, p in zip(coefficient, power))
    return 2*np.pi*integrated_square  # rho=1 in the benchmark


if __name__ == "__main__":
    y = sp.symbols("y", nonnegative=True)
    shape = (1-y**2)**3
    measure = 4*sp.pi*y**2
    norm = sp.integrate(measure*shape, (y, 0, 1))
    peak = sp.Rational(1, 4)/norm  # n_b = m = radius = 1
    carriage = sp.integrate(measure*(peak*shape)**2, (y, 0, 1))
    print("normal volume rate =", sp.simplify(peak*norm), "; aperture carriage =", carriage)
    last = None
    for terms in (512, 1024, 2048):
        force = outside_force(float(peak), terms)
        print(f"exterior terms={terms} receiver_on_intact={force:+.12g}")
        if last is not None and not (0 <= force-last < 1e-12):
            raise ArithmeticError("Exterior positive series did not stabilize")
        last = force
    if not (0 < force < float(carriage)/2):
        raise ArithmeticError("Bernoulli pressure or orientation mismatch")
