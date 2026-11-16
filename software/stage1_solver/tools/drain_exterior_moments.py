"""Axisymmetric exterior potential driven by inward flux on a finite sphere.

The pole points towards the target. The calculation concerns the steady
source-only excess in a constant-density, unbounded 3D potential exterior.
"""
import sympy as sp


r, a, q, d = sp.symbols("r a Q d", positive=True)
mu = sp.Symbol("mu", real=True)
p, s = sp.symbols("p s", real=True)
normal_speed = -q * (1 + p*mu + s*(3*mu**2 - 1)/2) / (4*sp.pi*a**2)
inward_flux = -2*sp.pi*a**2 * sp.integrate(normal_speed, (mu, -1, 1))
print("integrated inward source flux =", sp.factor(inward_flux))

# Determine the exterior radial power from Laplace's equation for each
# Legendre mode, omitting the constant/growing solution at infinity.
power = sp.Symbol("power", real=True)
potential = 0
for degree in range(3):
    angular = sp.legendre(degree, mu)
    angular_laplacian = sp.diff((1-mu**2)*sp.diff(angular, mu), mu)
    equation = sp.simplify((power*(power+1)*angular + angular_laplacian) / angular)
    exponent, = [root for root in sp.solve(equation, power) if root < 0]
    projection = (2*degree+1)*sp.integrate(normal_speed*angular, (mu, -1, 1))/2
    coefficient = sp.simplify(projection / (exponent * a**(exponent-1)))
    print("degree, radial exponent, potential coefficient =", degree, exponent, coefficient)
    potential += coefficient * r**exponent * angular

speed = sp.diff(potential, r)
boundary_error = sp.simplify(speed.subs(r, a) - normal_speed)
flux = sp.factor(-2*sp.pi*r**2 * sp.integrate(speed, (mu, -1, 1)))
# Inward speed at the target is the negative outward radial component.
scale = 4*sp.pi*d**2/q
axial = sp.factor(-scale * speed.subs({r: d, mu: 1}))
equatorial = sp.factor(-scale * speed.subs({r: d, mu: 0}))
print("source BC residual, exterior sphere intake =", boundary_error, flux)
print("scaled radial inflow [pole, equator] =", axial, equatorial)
cases = {"round": (0, 0), "dipole": (sp.Rational(1,2), 0),
         "quadrupole": (0, sp.Rational(1,2))}
for label, (dipole, quadrupole) in cases.items():
    print(label, "scaled pole readings d/a=2,5,10 =",
          [axial.subs({p: dipole, s: quadrupole, a: 1, d: distance})
           for distance in (2, 5, 10)])
round_speed = -speed.subs({mu: 1, p: 0, s: 0, r: d})
print("round, doubled intake, speed ratio =", round_speed, round_speed.subs(q, 2*q),
      sp.simplify(round_speed.subs(q, 2*q) / round_speed))
if boundary_error != 0 or sp.simplify(flux - inward_flux) != 0:
    raise AssertionError("Exterior solution does not reproduce the prescribed intake")
