"""Mixed stream/motion energy budget at fixed mass intake and far heads.

Use signed surface fluxes for ideal outer domains, then separately test a
speed-independent, passive head-loss law. All power signs are laboratory signs.
The normal-jet example is a constitutive control, not a measured throat.
"""
import math
from itertools import product

import sympy as sp


def cut(dimension, q, stream, center, susceptibility, radius):
    mu = sp.Symbol("mu", real=True)
    rho, area = sp.symbols("rho area", positive=True)
    far_head = sp.Symbol("H", real=True)  # Far p/rho. A constant internal energy adds q*e.
    radial = q / (rho * area)
    dipole = susceptibility * (stream - center) / radius**dimension
    vn = radial + (stream - (dimension - 1)*dipole)*mu
    tangent = stream + dipole
    speed2 = vn**2 + tangent**2 * (1-mu**2)
    vx = vn*mu + tangent*(1-mu**2)
    crossing = vn - center*mu
    pressure = rho*(far_head + (stream**2 - speed2)/2 + center*(vx-stream))

    def integrate(expr):
        # Spherical first three angular moments are sufficient for these cuts.
        poly = sp.Poly(expr, mu)
        if poly.degree() > 3:
            raise ValueError("Cut needs a higher angular moment")
        return sp.factor(area*(poly.nth(0) + poly.nth(2)/dimension))

    energy = integrate(rho*speed2*crossing/2)
    work = integrate(pressure*vn)
    impulse = integrate(rho*vx*crossing + pressure*mu)
    return energy, work, impulse, far_head


def odd_both(fn, u, v):
    return sum(a*b*fn(a*u, b*v) for a, b in product((1, -1), repeat=2))/4


def exterior_check():
    u, w, v, d = sp.symbols("U W V D", real=True)
    for radius, response in ((1, 0), (1, sp.Rational(1,4)), (2, sp.Rational(1,4))):
        brane = cut(3, -d, u, v, response, radius)
        bulk = cut(4, d, w, v, -response, radius)
        energy_terms = [sp.factor(brane[i] + bulk[i]) for i in (0,1)]
        force = sp.factor(brane[2] + bulk[2])
        net_each = [sp.factor(sum(exterior[:2])) for exterior in (brane,bulk)]
        print(f"cut={radius} susceptibility={response} "
              f"outward_flux_per_signed_rate={[sp.factor(net_each[0]/-d), sp.factor(net_each[1]/d)]} "
              f"holder={force}")
        ratios = [sp.factor(odd_both(lambda a,b: term.subs({u:a, v:b, w:0}, simultaneous=True), u,v)
                             / (d*u*v)) for term in (*energy_terms, force*v)]
        print(f"mixed_per_DUV transport={ratios[0]} pressure={ratios[1]} actuator={ratios[2]}")
        assert sum(ratios[:2]) == 0 and ratios[2] == -1
        assert net_each[0].expand() == (-d*(brane[3]+u**2/2)).expand()
        assert net_each[1].expand() == (d*(bulk[3]+w**2/2)).expand()
        assert force == d*(w-u) or sp.expand(force-d*(w-u)) == 0


def head_control():
    # An ordinary jet lost on mixing gives a rest loss D²/[2(rho*A)²].
    # Fix rho*A=1 and D0=1/4. The stationary slope is independently available
    # from the rest head curve. Compare its tangent prediction to the exact
    # root, without fitting the four moving states.
    area_density, d0 = 1., .25
    def loss(rate):
        return (rate/area_density)**2/2
    baseline = loss(d0)
    step = d0 / 16
    slope = (loss(d0+step)-loss(d0-step))/(2*step)
    print(f"jet rho_area={area_density:g} rest_D={d0:g} rest_head={baseline:g} "
          f"rest_head_slope={slope:g}")
    for u, v in ((.01,.02), (.005,.01)):
        def released_rate(stream, center):
            return area_density * math.sqrt(2*baseline+stream**2-2*stream*center)
        for a,b in product((1,-1), repeat=2):
            stream, center = a*u, b*v
            need = stream*center - stream**2/2
            flow = released_rate(stream, center)
            incoming = d0*stream**2/2
            delivered = d0*stream*center
            print(f"U={stream:+g} V={center:+g} incremental_arrival={incoming:.12g} "
                  f"delivered={delivered:.12g} fixed_D_extra_power={d0*need:.12g} "
                  f"fixed_D_extra_head={need:.12g} fixed_head_jet_D={flow:.12g}")
            assert math.isclose(delivered-incoming,d0*need,abs_tol=1e-17)
            assert math.isclose(loss(flow),baseline-need,abs_tol=1e-16)
        observed = odd_both(released_rate, u, v)/(u*v)
        predicted = -1/slope
        print(f"U_abs={u:g} V_abs={v:g} fixed_head_jet_mixed_D_over_UV={observed:.12g} "
              f"rest_slope_prediction={predicted:g}")
        assert abs(observed-predicted) < .005


if __name__ == "__main__":
    exterior_check()
    head_control()
