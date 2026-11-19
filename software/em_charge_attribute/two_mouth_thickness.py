"""Work of moving two mirror-related wave-supported mouths on an elastic slab."""

import sympy as sp


def check(label, actual, expected):
    print(f"{label} = {sp.factor(actual)}", flush=True)
    if sp.simplify((actual - expected).rewrite(sp.exp)) != 0:
        raise AssertionError(label)


# Uniform, unit-integral spherical surface weights of radius a in the brane.
# The angular average at distance d from one shell uses rho d(rho)/(2*a*d).
# Kernels below use units 4*pi*kappa=1; each physical mode has its own kappa.
a, d, mu, rho = sp.symbols("a d mu rho", positive=True)
form = sp.sinh(mu*a)/(mu*a)
one = sp.integrate(sp.exp(-mu*rho), (rho, d-a, d+a))/(2*a*d)
check("one shell outside", one, form*sp.exp(-mu*d)/d)
self_kernel = sp.integrate(sp.exp(-mu*rho), (rho, 0, 2*a))/(2*a**2)
check("shell self response", self_kernel, form*sp.exp(-mu*a)/a)
pair = form**2*sp.exp(-mu*d)/d
check("nonoverlapping double shell average", form*one, pair)
check("massless shell pair limit", sp.limit(pair,mu,0), 1/d)
check("unrestored thickness, equal stiffness, opposite-face response",
      sp.limit(pair-1/d,mu,0), 0)
check("far oriented energy coefficient divided by isolated squared load",
      sp.limit(-d*((1/d+pair)-(-1/d+pair))/2, d, sp.oo), -1)

# Outward depth on the loaded face is s*m+t. In the symmetric pair the
# stationary field energy is z^2/(S+g), before relaxing the two mouths.
S, g, kb, A, T = sp.symbols("S g k_b A T", positive=True)
z, Y, rear, P = sp.symbols("z Y Y_ref P", real=True)
length = Y-z
energy = z**2/(S+g) + 2*(T*length+A/length) + kb*(Y-rear)**2
residuals = [energy.diff(q)/2 for q in (z,Y)]
on_branch = {z:(S+g)*P}
check("stationary energy derivative in mutual response", energy.diff(g).subs(on_branch), -P**2)
for s1, s2 in [(1,1),(1,-1),(-1,1),(-1,-1)]:
    force = P**2*sp.diff(s1*s2/d+pair, d)
    check(f"outward force for {s1:+d} {s2:+d}", force,
          -P**2*(s1*s2+(1+mu*d)*form**2*sp.exp(-mu*d))/d**2)

# The projected gap and the gap measured perpendicular to two locally
# parallel faces agree to linear order in small displacements and slopes.
e, slope, thick, H = sp.symbols("eps slope t H", real=True)
gap_change = (H+2*e*thick)/sp.sqrt(1+e**2*slope**2)-H
check("leading squared normal gap change", sp.diff(gap_change**2,e,2).subs(e,0)/2, 4*thick**2)

# Same rear reference, action, wall tension and radius in all eight solves.
# Energies and forces are divided by P0^2 in these dimensionless units.
p0 = sp.Rational(1,20)
self_value = 1+self_kernel.subs({a:1,mu:1})
material = {S:self_value, kb:1, A:1, T:1+p0, rear:1+(1+self_value)*p0}
single_pair = {z:self_value*p0, Y:1+self_value*p0, g:0}
base = energy.subs(material).subs(single_pair)
print("s1 s2 distance excess_energy/P0^2 outward_force/P0^2 load/P0", flush=True)
for sa, sb in [(1,1),(1,-1),(-1,1),(-1,-1)]:
    kernel = sa*sb/d + pair.subs({a:1,mu:1})
    for dist in (8,4):
        mutual = kernel.subs(d, dist)
        params = {**material, g:mutual}
        equations = [eq.subs(params) for eq in residuals]
        rim, tip = sp.nsolve(equations,(z,Y),(0.07,1.07),prec=45)
        state = {z:rim, Y:tip}
        load = rim/(self_value+mutual)
        excess = (energy.subs(params).subs(state)-base)/p0**2
        force = load**2*kernel.diff(d).subs(d,dist)/p0**2
        print(sa, sb, dist, *(f"{float(x):+.10f}" for x in (excess,force,load/p0)), flush=True)
        if tip <= rim or any(abs(eq.subs(state)) > 1e-35 for eq in equations):
            raise ArithmeticError("invalid equilibrium")
