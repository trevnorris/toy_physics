"""Angular averages on an upright, stationary cut across an open sleeve.

Take a radial isolated state and its first in-plane dipole perturbation in
an ungauged condensate. The profiles are free functions, evaluated on the
positive-density part of the cut. No core boundary condition is selected.
"""
import sympy as sp


r = sp.Symbol("r", positive=True)
w, mu, amplitude = sp.symbols("w mu epsilon", real=True)
m, hbar = sp.symbols("m hbar", positive=True)

# Log density and velocity potential. mu is the cosine to the companion.
log_base, log_dipole, phi_base, phi_dipole = (
    sp.Function(name)(r, w) for name in ("H", "B", "F", "G"))
log_density = log_base + amplitude*mu*log_dipole
density = sp.exp(log_density)
potential = phi_base + amplitude*mu*phi_dipole


def along_companion(field):
    return mu*sp.diff(field, r) + (1-mu**2)*sp.diff(field, mu)/r


mass_current = m*density*sp.diff(potential, w)
flow_current = mass_current*along_companion(potential)
quantum_current = -hbar**2*density*along_companion(sp.diff(log_density, w))/(4*m)

for name, current in (("mass", mass_current), ("carried", flow_current),
                      ("density-gradient stress", quantum_current)):
    isolated = current.subs(amplitude, 0)
    linear = sp.diff(current, amplitude).subs(amplitude, 0)
    for label, term in (("isolated", isolated), ("dipole derivative", linear)):
        average = sp.factor(sp.integrate(term, (mu, -1, 1))/2)
        print(name, label, "angular mean =", average)
