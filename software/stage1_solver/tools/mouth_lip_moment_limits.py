"""Exact thin-front moments of the symmetric wave-free order functional in kink units."""
import json

import sympy as sp


def main():
    t, s, a, y = sp.symbols("t s a y", real=True)
    # With t=chi_0, the magnitude of chi_0' is t*(1-t).
    speed = t*(1-t)
    # Integral of |chi_0'|^2 to the outer bath on the right half.
    tail = sp.integrate(y*(1-y), (y, 0, t))
    tension = 2*sp.integrate(speed, (t, 0, 1))
    bending = -4*sp.integrate(sp.cancel(tail**2/speed**3), (t, 0, sp.Rational(1, 2)))
    # The second moment is twice d^2/da^2 of integral t^(1-a)*(1-t)^(1+a) dt.
    tilted = sp.gamma(2-a)*sp.gamma(2+a)/sp.gamma(4)
    second = 2*sp.diff(tilted, a, 2).subs(a, 0)
    # A three-tangent symmetric strain film provides a separate wave-off bend example.
    k, H, young = sp.symbols('k H Y', positive=True)
    strip_energy = sp.integrate(young*(k*s)**2/2, (s, -H/2, H/2))
    print(json.dumps(dict(tension=str(tension), bending_over_tension=str(bending/tension),
                          jacobian_second_over_tension=str(sp.simplify(second/tension)),
                          independent_strip_modulus=str(sp.diff(strip_energy, k, 2)))))


if __name__ == '__main__':
    main()
