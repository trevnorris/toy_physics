"""Check normal-current matching and what stationary measurements fix about a breathing pole."""
import json

import sympy as sp


def check():
    n, inertia, gradient, width = sp.symbols('n I g delta', positive=True)
    slope, speed, velocity, pi = sp.symbols('p V u pi', real=True)
    # The normal points along grad(chi); g=|grad(chi)|. V is the speed of
    # the chi level in that same direction. Both equations are instantaneous.
    sol = sp.solve([n*velocity - n*slope + pi*gradient,
                    pi - inertia*gradient*(velocity-speed)], (velocity, pi), dict=True)[0]
    relative = sp.factor(n*(sol[velocity]-speed))
    expected = n**2*(slope-speed)/(n+inertia*gradient**2)
    assert sp.simplify(relative-expected) == 0
    # Radius-independent translation of half of the planar logistic front.
    # |d chi/d ell|=chi(1-chi)/width, in monotonically increasing orientation.
    chi = sp.symbols('chi', real=True)
    excess = sp.integrate(inertia*chi*(1-chi)/width, (chi, 0, sp.Rational(1, 2)))
    assert sp.simplify(excess.subs(inertia,width**2)-width/12) == 0

    # Counterexample at the level of static data, never a fitted mouth law.
    # The two test coordinates have independent unit inertias, an additional
    # mixed inertia s, and the same passive outgoing admittance on coordinate y.
    # Eliminate y. The unit impedance is the load of a normalized outgoing line.
    w, b, s = sp.symbols('omega b s', real=True)
    mech = 3-w**2
    bath = b+2-w**2-sp.I*w
    dynamic = sp.Matrix([[mech,-s*w**2],[-s*w**2,bath]])
    reduced = sp.factor(dynamic.det()/bath)
    assert sp.simplify(reduced-(mech-s**2*w**4/bath)) == 0
    assert sp.simplify(reduced.subs(w,0)-3) == 0
    loss, reactive, derivative = [], [], []
    for mixed_mass in (0,sp.Rational(1,2)):
        point = {b:1,w:1,s:mixed_mass}
        assert 1-mixed_mass**2 > 0
        z = sp.expand_complex(reduced.subs(point))
        loss.append(str(-sp.im(z)))
        reactive.append(str(sp.re(z)))
        derivative.append(str(sp.re(sp.diff(reduced,b).subs(point).expand(complex=True))))
    assert reactive==['2','19/10'] and loss==['0','1/20']
    assert derivative==['0','3/100']
    return {'relative_current':str(relative),'order_momentum':str(sp.factor(sol[pi])),
            'integrated_half_front_weight':str(excess),
            'at_I_equal_width_squared':str(excess.subs(inertia,width**2)),
            'example':'dimensionless passive two-port, e^(-i omega t), unit frequency and b',
            'mixed_inertias':['0','1/2'],'reduced_real_drive_stiffness':reactive,
            'minus_imaginary_drive_stiffness':loss,
            'exit_derivative_of_real_drive_stiffness':derivative}


if __name__ == '__main__':
    print(json.dumps(check(), sort_keys=True))
