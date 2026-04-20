#!/usr/bin/env python3
"""
Stage 77 SymPy audit

Checks the exact l=0 / l=2 decoupling on the isotropic wall branch by explicit
real-harmonic integrals on S^2.
"""

from __future__ import annotations
import sympy as sp

pi = sp.pi
th, ph = sp.symbols('th ph', real=True)


def domega(expr):
    return sp.simplify(sp.integrate(sp.integrate(sp.simplify(expr*sp.sin(th)), (ph, 0, 2*sp.pi)), (th, 0, sp.pi)))


def lap_s2(Y):
    return sp.simplify((1/sp.sin(th))*sp.diff(sp.sin(th)*sp.diff(Y, th), th) + (1/sp.sin(th)**2)*sp.diff(Y, ph, 2))


Y00 = sp.Rational(1, 2)/sp.sqrt(pi)
Y20 = sp.sqrt(sp.Rational(5, 16)/pi) * (3*sp.cos(th)**2 - 1)
Y21c = sp.sqrt(sp.Rational(15, 4)/pi) * sp.sin(th)*sp.cos(th)*sp.cos(ph)
Y21s = sp.sqrt(sp.Rational(15, 4)/pi) * sp.sin(th)*sp.cos(th)*sp.sin(ph)
Y22c = sp.sqrt(sp.Rational(15, 16)/pi) * sp.sin(th)**2 * sp.cos(2*ph)
Y22s = sp.sqrt(sp.Rational(15, 16)/pi) * sp.sin(th)**2 * sp.sin(2*ph)
Y2 = {'20': Y20, '21c': Y21c, '21s': Y21s, '22c': Y22c, '22s': Y22s}

print('Y00 normalization =', domega(Y00*Y00))
for name, Y in Y2.items():
    print(f'Y{name} normalization =', domega(sp.simplify(Y*Y)))
    print(f'<Y00|Y{name}> =', domega(sp.simplify(Y00*Y)))
    print(f'(-Delta)Y{name} - 6Y{name} =', sp.simplify(-lap_s2(Y) - 6*Y))
    grad_cross = sp.simplify(sp.diff(Y00, th)*sp.diff(Y, th) + (1/sp.sin(th)**2)*sp.diff(Y00, ph)*sp.diff(Y, ph))
    print(f'<grad Y00 . grad Y{name}> =', domega(grad_cross))
    print(f'<Y00|(-Delta)Y{name}> =', domega(sp.simplify(Y00*(-lap_s2(Y)))))
    assert domega(sp.simplify(Y00*Y)) == 0
    assert domega(grad_cross) == 0
    assert domega(sp.simplify(Y00*(-lap_s2(Y)))) == 0

mu, Tw, TOm, K = sp.symbols('mu Tw TOm K', real=True)
# Generic isotropic quadratic cross coefficient between l=0 and one l=2 mode.
# Time and w derivatives factor out and multiply the same angular orthogonality integrals.
for name, Y in Y2.items():
    I_mass = domega(Y00*Y)
    I_grad = domega(sp.diff(Y00, th)*sp.diff(Y, th) + (1/sp.sin(th)**2)*sp.diff(Y00, ph)*sp.diff(Y, ph))
    I_lap = domega(sp.simplify(Y00*(-lap_s2(Y))))
    I_pot = domega(Y00*Y)
    Ccross = sp.simplify(mu*I_mass - Tw*I_mass - TOm*I_lap - K*I_pot)
    print(f'Generic isotropic cross coefficient C_0,{name} =', Ccross)
    assert sp.simplify(Ccross) == 0

print('\nStage 77 theorem verified: isotropic l=0 <-> l=2 cross terms vanish exactly.')
