"""Excess energy of a charged patch on an existing planar order wall."""
import sympy as sp


# The normal coordinate is traded for the monotone order of a logistic wall.
chi = sp.symbols("chi", nonnegative=True)
a, kappa, ell, area, spectator, charge, density_scale, coulomb, radius, tension = (
    sp.symbols("a kappa ell A L Q nu0 C r T", positive=True)
)
well = chi**2 * (1 - chi)**2
slope = chi * (1 - chi) / ell


def energy(well_coefficient):
    density = kappa * slope**2 / 2 + well_coefficient * well
    return sp.integrate(sp.cancel(density / slope), (chi, 0, 1))


neutral = energy(a)
delta = sp.solve(neutral.diff(ell), ell)[0]
surface = sp.simplify(neutral.subs(ell, delta))
print("neutral column energy, minimizing width, tension =", neutral, delta, surface)
fixed_width = sp.symbols("ell_fixed", positive=True)
fixed_excess = sp.factor((neutral.subs(ell, fixed_width) - surface).subs(kappa, 2*a*ell**2))
print("imposed-width excess per wall three-area (ell = neutral width) =", fixed_excess)

# Counterexample: the well softens in response to the local areal charge
# magnitude. No change of sign or disappearance of the double well is involved.
area_scale = charge / density_scale
ratio = area_scale / area
softened = energy(a / (1 + ratio)**2)
optimum = sp.solve(softened.diff(ell), ell)[0]
relaxed = area * spectator * (softened.subs(ell, optimum) - surface)
frozen = area * spectator * (softened.subs(ell, delta) - surface)
prefactor = spectator * surface
scaled = sp.factor(relaxed / prefactor)
scaled_frozen = sp.factor(frozen / prefactor)
print("relaxed/neutral width, relaxed excess/(L*tension) =",
      sp.factor(optimum / delta), scaled)
print("area derivative and infinite-area limit of scaled relaxed excess =",
      sp.factor(scaled.diff(area)), sp.limit(scaled, area, sp.oo))
print("scaled frozen-trial excess, infinite-area limit =",
      scaled_frozen, sp.limit(scaled_frozen, area, sp.oo))
print("scaled excess above the dilute limit =", sp.factor(scaled + area_scale))

disk_energy = coulomb * charge**2 / radius + sp.pi * radius**2 * tension
disk_radius = sp.solve(disk_energy.diff(radius), radius)[0]
curvature = sp.simplify(disk_energy.diff(radius, 2).subs(radius, disk_radius))
print("positive excess-area cost: minimizing radius and curvature =", disk_radius, curvature)
assert sp.factor(scaled.diff(area)).is_negative
assert sp.simplify(sp.limit(scaled_frozen-scaled, area, sp.oo)) == 0
assert curvature.is_positive
