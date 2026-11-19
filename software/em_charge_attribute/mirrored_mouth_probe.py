"""Weak temporal probe, reflection, and removal of carrier residence charge."""

import sympy as sp


def record(name, value, expected):
    result = sp.factor(value)
    print(f"{name} = {result}", flush=True)
    if sp.simplify(result - expected) != 0:
        raise ArithmeticError(name)


q, probe, n, hbar, phase_rate, f = sp.symbols("q probe n hbar phase_rate f", real=True)
covariant_rate = sp.I * (phase_rate + q * probe * f / hbar)
temporal_action = sp.I*hbar*n*(covariant_rate-sp.conjugate(covariant_rate))/2
record("weak static source", -temporal_action.diff(probe), q*n*f)

w, lam = sp.symbols("w lambda", positive=True)
z = sp.exp(-(w/lam)**2)
for profile in (sp.Integer(1), 2*w/lam):
    eigenvalue = -sp.diff(z*profile.diff(w), w)/(z*profile)
    record(f"depth eigenvalue for {profile}", eigenvalue, 0 if profile == 1 else 2/lam**2)

# Independent variables for response of any reflected transit-subtracted carrier
# profile at the two points +w and -w. No even-density assumption per mouth.
p, m, fplus, fminus = sp.symbols("p m fplus fminus")
upper = q*(p*fplus + m*fminus)
lower = q*(m*fplus + p*fminus)
odd = (upper - lower)/2
record("odd source per mirror pair", odd, q*(fplus-fminus)*(p-m)/2)
record("flat orientation difference", odd.subs({fplus: 1, fminus: 1}), 0)
record("first depth moment", odd.subs({fplus: 2*w/lam, fminus: -2*w/lam}), 2*q*w*(p-m)/lam)

# Positive control: spread a unit over the angular sphere and a smooth
# conversion interval. Keep its normalization fixed while the radius varies.
r, u, g = sp.symbols("R u g", positive=True)
chi = 1-3*u**2+2*u**3
layer_areal_density = -g*chi.diff(u)/(4*sp.pi*r**2)
layer_response = sp.integrate(4*sp.pi*r**2 * layer_areal_density, (u, 0, 1))
for radius in (1, 2):
    record(f"added layer oriented residual R={radius}", (odd.subs({fplus: 1, fminus: 1})
           + layer_response.subs(r, radius)), g)
