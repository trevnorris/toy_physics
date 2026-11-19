"""Exact initial torque on a standing wave before the split, with zero velocity."""

import sympy as sp

x, w, eta = sp.symbols("x w eta", real=True)
V = x**2 + w**4/4
B = x + eta*(w**2-sp.Rational(1, 2))  # the omitted amplitude factor is exp(-r^2/2)


def gradient(poly, axis):
    return sp.diff(poly, axis) - axis*poly


def xmoment(power):
    return 0 if power % 2 else sp.gamma(sp.Rational(power+1, 2))


def half(poly, side):
    return sp.factor(sum(coef*xmoment(i)*side**j*sp.gamma(sp.Rational(j+1, 2))/2
                         for (i,j), coef in sp.Poly(poly, x, w).terms()))


ux, uw = gradient(B, x), gradient(B, w)
acc = gradient(ux, x) + gradient(uw, w) - V*B
energy = (ux**2 + uw**2 + V*B**2)/2
dpx, dpw = -acc*ux, -acc*uw
for side in (1, -1):
    total = half(energy, side)
    raw = half(x*dpw-w*dpx, side)
    origin_motion = (half(x*energy, side)*half(dpw, side)
                     - half(w*energy, side)*half(dpx, side))/total
    internal = sp.factor(raw-origin_motion)
    well_torque = half(-B**2*(x*V.diff(w)-w*V.diff(x))/2, side)
    cut_stress = sp.Poly((x*(uw**2-energy)).subs(w, 0), x)
    cut_torque = side*sum(c*xmoment(i) for (i,), c in cut_stress.terms())
    print("side", side, "raw, trap and cut torque =", raw, well_torque, sp.factor(cut_torque))
    print("side", side, "intrinsic torque linear coefficient =", sp.factor(internal.diff(eta).subs(eta, 0)))
    for tilt in (sp.Rational(1,100), -sp.Rational(1,100), 0):
        print("side", side, "tilt", tilt, "initial intrinsic torque =", float(internal.subs(eta, tilt)))
    assert sp.simplify(raw-well_torque-cut_torque) == 0
    assert sp.simplify(internal + internal.subs(eta, -eta)) == 0
