"""Eliminate the far displacement of a translating rim and compare field loads."""

import sympy as sp


def check(label, actual, expected):
    actual = sp.factor(actual)
    print(f"{label} = {actual}", flush=True)
    if sp.simplify(actual - expected) != 0:
        raise ArithmeticError(label)


def main():
    a, b, c, d, load, load2 = sp.symbols("A B C D K1 K2")
    x, y = sp.symbols("X Y", nonzero=True)
    y_free = sp.solve(c*x + (d+load)*y, y)[0]
    rim = sp.factor((a*x+b*y_free)/x)
    check("force per rim displacement", rim, a-b*c/(d+load))
    check("change for a different end load", rim.subs(load, load2)-rim,
          b*c*(load2-load)/((d+load)*(d+load2)))

    # For the passive wall, integrate its Lagrange identity between endpoints.
    # p0, ps are unit-rim-displacement solutions at s=0 and s>0 respectively.
    w = sp.symbols("w")
    s, tau, z1, z2 = sp.symbols("s tau Z1 Z2", positive=True)
    p0, ps, tension, density, pin = [sp.Function(n)(w) for n in
                                       ("p0","ps","T","mu","H")]
    stress0, stress_s = tension*p0.diff(w), tension*ps.diff(w)
    divergence = p0*stress_s.diff(w) - ps*stress0.diff(w)
    check("interior concomitant derivative", divergence.xreplace({
        stress_s.diff(w): (pin+density*s**2)*ps,
        stress0.diff(w): pin*p0}), density*s**2*p0*ps)

    # In each unpinned layer, [displacement, tension * spatial derivative]
    # propagates from its left end to its right end. The last force is zero.
    h, ch = sp.symbols("h ch", positive=True)  # tanh(s tau), cosh(s tau)
    def transfer(z):
        return ch*sp.Matrix([[1,h/(s*z)], [s*z*h,1]])

    def free_end(left, right):
        system = transfer(right)*transfer(left)
        return sp.factor(system[1,0]/system[1,1])

    near1, near2 = free_end(z1,z2), free_end(z2,z1)
    check("end stiffness with layer 1 nearest", near1,
          s*h*(z1+z2)/(1+h**2*z2/z1))
    check("end stiffness with layer 2 nearest", near2,
          s*h*(z1+z2)/(1+h**2*z1/z2))
    print("end stiffness difference =", sp.factor(near1-near2))
    for response in (near1,near2):
        exact = response.subs(h,sp.tanh(s*tau))
        check("static stiffness", sp.limit(exact,s,0), 0)
        check("slow force per acceleration", sp.limit(exact/s**2,s,0), tau*(z1+z2))
    check("first differing coefficient (s^4)",
          sp.limit(((near1-near2)/s**4).subs(h,sp.tanh(s*tau)),s,0),
          tau**3*(z1+z2)*(z1/z2-z2/z1))


if __name__ == "__main__":
    main()
