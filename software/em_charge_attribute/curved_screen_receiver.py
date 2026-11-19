"""Net tangential momentum across a curved-screen receiver's injection cut.

Integrate the excess Bernoulli pressure on the *intact* graph. Compare a
point probe with a finite uniform spherical source of the prescribed height.
"""

from math import log2, pi, sqrt

import sympy as sp
from scipy.integrate import quad


def check(label, value, expected):
    print(f"{label} = {sp.factor(value)}", flush=True)
    if sp.simplify(value - expected) != 0:
        raise AssertionError(label)


def symbolic():
    eps, px, slope, local, pressure, dp = sp.symbols("eps x gx g0 p0 p1", real=True)
    tang = sp.Matrix([1,eps*local])
    wall_area_normal = sp.Matrix([eps*slope,-1])
    load = (pressure+eps*dp)*wall_area_normal.dot(tang)
    check("screen flux per projected area at first order", load.diff(eps).subs(eps,0),
          pressure*(slope-local))
    check("common tilt", load.subs(slope,local), 0)
    r,d,a,rho,rate = sp.symbols("r d a rho D", positive=True)
    e = sp.symbols("e",real=True)
    # Reference upper-hemisphere solution restricted to its intact plane.
    q = rate/(sp.pi**2*rho)
    potential = -q/(2*r**2) + q*e*a**2/(4*r**4)
    p0 = -rho*potential.diff(r)**2/2
    tail = -4*sp.pi*sp.integrate(sp.expand(r**2*p0),(r,d,sp.oo))
    check("positive pressure-tail weight", tail,
          2*rate**2/(sp.pi**3*rho*d**3)*(sp.Rational(1,3)-2*e*a**2/(5*d**2)+e**2*a**4/(7*d**4)))
    height_charge = sp.symbols("A", positive=True)
    injection = tail * height_charge/d**2
    check("coefficient of d^-2 on the inner cut", sp.limit(d**2*injection,d,sp.oo),0)
    check("coefficient of d^-5 on the inner cut", sp.limit(d**5*injection,d,sp.oo),
          2*height_charge*rate**2/(3*sp.pi**3*rho))
    # For a uniform ball of radius b, the second moment per axis is b^2/5.
    y,z,b = sp.symbols("y z b", real=True)
    rr = sp.sqrt(px**2+y**2+z**2)
    component = px/rr**6
    laplacian = sum(component.diff(v,2) for v in (px,y,z))
    check("finite uniform source second-order coefficient relative to leading term",
          (b**2*laplacian/(10*component)).subs({px:d,y:0,z:0}) - 6*e*a**2/(5*d**2),
          (9*b*b-6*e*a*a)/(5*d*d))


def mean_gradient(r, d, b):
    """Direct angular integral of the gradient of a uniform source ball.

    Separation t runs over the shell around the *target*. Inside the
    source the field is linear, outside it is inverse square.
    """
    low, high = abs(d-r), d+r
    def integrand(t):
        return (d*d-r*r+t*t)*t/(4*r*d*d*max(b,t)**3)
    edges = [low] + ([b] if low < b < high else []) + [high]
    return sum(quad(integrand,left,right,epsabs=2e-17,epsrel=2e-13)[0]
               for left,right in zip(edges,edges[1:]))


def finite_probe(d, b, a, e):
    # A=D=rho=1; the small-slope parameter is differentiated at zero.
    center_gradient = 1/d**2
    def weighted_pressure(r):
        return -2*(1-e*(a/r)**2)**2/(pi**3*r**4)
    def tail(radius):
        return -quad(weighted_pressure,radius,float("inf"),epsabs=1e-18)[0]
    assert abs(mean_gradient(0.5*(a+d-b),d,b)/center_gradient-1) < 1e-12
    assert abs(mean_gradient(d+2*b,d,b)/center_gradient) < 1e-12
    def screen(r):
        return weighted_pressure(r) * (mean_gradient(r,d,b)-center_gradient)
    partial = sum(quad(screen,left,right,epsabs=1e-21,epsrel=2e-12)[0]
                  for left,right in ((d-b,d),(d,d+b)))
    current = partial + center_gradient*tail(d+b)
    if not tail(d+b) <= current/center_gradient <= tail(d-b):
        raise AssertionError("positive compact probe violates the pressure-tail bounds")
    return current


if __name__ == "__main__":
    symbolic()
    for e in (0, 0.25):
        previous = None
        for d in (8,16,32,64):
            current = finite_probe(d,b=1,a=1,e=e)
            scaled = current*d**5*3*pi**3/2
            print(f"shape {e:g}, d {d}, normalized d^5 J = {scaled:.12f}")
            if previous is not None:
                print(f"shape {e:g}, d {d//2}->{d}, force exponent = {log2(previous/current):.10f}")
            previous = current
        if abs(scaled-1) > 0.001:
            raise AssertionError("normalized asymptote")
