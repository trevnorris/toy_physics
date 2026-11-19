"""Leading one-mouth response in an incident bend, with explicit momentum ports."""

import json
import re
from pathlib import Path

import sympy as sp


def verify(label, actual, expected):
    value = sp.factor(actual)
    print(f"{label} = {value}", flush=True)
    if sp.simplify(value - expected) != 0:
        raise ArithmeticError(label)


def scalar_cut():
    # The unit normal is in the x,w-independent in-brane section; integrate
    # in three spatial dimensions. cos(theta) is measured from the gradient.
    mu, a, source, k, g = sp.symbols("mu r L k g", real=True, nonzero=True)
    q = -source/(4*sp.pi*k*a**2)
    # Linear term of grad(h)_x grad(h)_n - n_x |grad h|^2/2.
    cross = k*((q*mu)*(g*mu) + g*q - mu*q*g*mu)
    incoming = -2*sp.pi*a**2*sp.integrate(cross, (mu, -1, 1))
    verify("inward scalar momentum in the gradient direction", incoming, source*g)


def moving_stream():
    tilt, ell, cin, cout, slope, s = sp.symbols("eta ell c_in c_out g s", real=True)
    speed, mass_rate, velocity = sp.symbols("v D U", real=True)
    hprime = ell*tilt*slope
    tangent = sp.Matrix([1, hprime])
    guide_velocity = velocity*tangent

    def material_velocity(c):
        return s*speed*sp.Matrix([-c*hprime, 1])/sp.sqrt(1+(c*hprime)**2)

    vin, vout = material_velocity(cin), material_velocity(cout)
    impulse = mass_rate*(vin-vout)
    load = impulse.dot(tangent)
    coeff = sp.diff(load, tilt).subs(tilt, 0)
    verify("tangent stream force, linear in incident slope", coeff,
           s*mass_rate*speed*ell*(cout-cin)*slope)
    for a, b in ((1, 0), (0, 1), (1, 1)):
        print(f"in/upstream follows {a}, out/downstream follows {b}:",
              sp.factor(coeff.subs({cin:a, cout:b})))
    power = mass_rate*((vin+guide_velocity).dot(vin+guide_velocity)
                        - (vout+guide_velocity).dot(vout+guide_velocity))/2
    verify("kinetic power delivered minus stream force times guide velocity",
           sp.factor(sp.expand(power - velocity*load)), 0)
    verify("resting guide kinetic power", power.subs(velocity, 0), 0)
    # At any nonzero slope a purely normal load does zero tangent work.
    verify("unit-following normal motion dotted with tangent", material_velocity(1).dot(tangent), 0)


def footprint():
    # Potential of a monopole at distance d, averaged over a spherical
    # sampling surface of radius b < d. Use the endpoints of the primitive
    # so the absolute-value branch is explicit.
    d, b = sp.symbols("d b", positive=True)
    mu = sp.symbols("mu", real=True)
    primitive = -sp.sqrt(d*d+b*b-2*d*b*mu)/(d*b)
    verify("shell integral primitive", primitive.diff(mu), 1/sp.sqrt(d*d+b*b-2*d*b*mu))
    gap = sp.symbols("gap", positive=True)
    ends = [sp.sqrt(sp.factor((d*d+b*b-2*d*b*end).subs(d,b+gap))) for end in (-1,1)]
    average = ((ends[0]-ends[1]).subs(gap,d-b))/(2*d*b)
    verify("shell-averaged field at either port", average, 1/d)
    verify("shell-averaged incident radial derivative", average.diff(d), -1/d**2)
    # A fixed displacement along the line of centers is deliberately
    # anisotropic. Its first differential signal has one extra power of 1/d.
    offset = sp.symbols("delta", real=True)
    verify("shifted gradient derivative at zero offset",
           sp.diff(-1/(d+offset)**2,offset).subs(offset,0), 2/d**3)


def bulk_cut():
    # 4D irrotational, constant-density steady receiver. On a sphere use
    # phi1=p*(n dot ex)/r^3. The monopole v0 is radial. Keeping the pressure
    # cross term is necessary before taking the hemisphere moment.
    r, p, rho, q, mu = sp.symbols("R p rho Q mu", positive=True)
    v0 = q/r**3
    v1n = -3*p*mu/r**4
    v1x = p*(1-4*mu**2)/r**4
    advective = rho*v0*(v1x + mu*v1n)
    pressure = -rho*v0*v1n*mu
    def hemisphere(expr):
        # mu is a tangential direction cosine; half of S^3 selected by w>0.
        return sp.expand(expr*r**3*sp.pi**2).subs(mu**2,sp.Rational(1,4))
    verify("hemisphere advective x flux for a decaying bulk dipole", hemisphere(advective),
           -3*sp.pi**2*p*q*rho/(4*r**4))
    verify("hemisphere pressure x flux for a decaying bulk dipole", hemisphere(pressure),
           3*sp.pi**2*p*q*rho/(4*r**4))
    verify("hemisphere total x flux", hemisphere(advective+pressure), 0)
    verify("uniform incident bulk velocity, surviving coefficient per unit velocity",
           hemisphere(rho*v0*(1+mu**2)-rho*v0*mu**2), sp.pi**2*q*rho)
    verify("independent dipolar pressure, coefficient for p*n_x/R^3",
           hemisphere(p*mu**2/r**3), sp.pi**2*p/4)


def radii():
    manifest = json.loads(Path("docs/_measurements/unequal_mouth_electric.json").read_text())
    readings = manifest["runs"]["annular_mouth_radius_sympy"]["stdout_claim_lines"]
    loads = {}
    for line in readings:
        result = re.fullmatch(r"radius (\d+), escape = ([\d.]+), self response = [\d.]+", line)
        if result:
            loads[int(result[1])] = sp.Rational(result[2])
    ratio = loads[2]/loads[1]  # positive sources of equal bare strength
    print("same target, emitted strengths 2a/a =", ratio.evalf(12))
    print("equalized diagonal products, ratio of the two ordered mixed pairs =", (ratio**2).evalf(12))
    rate, rho, radius, shape = sp.symbols("D rho a beta", positive=True)
    area = 4*sp.pi*radius**3/3
    thrust = shape*rate**2/(rho*area)
    verify("three-ball thrust radius ratio with same density and shape",
           thrust.subs({radius:2*radius,rate:sp.Symbol('D2')}) / thrust,
           sp.Symbol('D2')**2/(8*rate**2))


if __name__ == "__main__":
    scalar_cut()
    moving_stream()
    footprint()
    bulk_cut()
    radii()
