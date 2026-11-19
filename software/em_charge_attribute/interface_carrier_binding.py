"""Weak stationary number response to the existing single-trace contact law."""
import sympy as sp


def check(label, value, expected):
    result = sp.factor(value)
    print(f"{label} = {result}", flush=True)
    if sp.simplify(result - expected) != 0:
        raise AssertionError(label)


# The stationary GNLS amplitude at fixed number current per true area.
m, hb, eos, n0, j = sp.symbols("m hbar K n0 j", positive=True)
z, eps = sp.symbols("z eps", real=True)
eta = sp.Function("dn")(z)
n = n0 + eps*eta
amp = sp.sqrt(n)
mu = lambda density: 5*eos*density**4/4
head = -hb**2*amp.diff(z, 2)/(2*m*amp) + mu(n) + m*j**2/(2*n**2)
linear = sp.expand(head.diff(eps).subs(eps, 0))
gradient = hb**2/(4*m*n0)
stiffness = 5*eos*n0**3-m*j**2/n0**3
check("linear density coefficient at fixed current", linear.coeff(eta), stiffness)
check("density second derivative coefficient", linear.coeff(eta.diff(z, 2)), -gradient)

# In the continuous-trace surface law c [amp'] = (partial_n g) amp,
# weak coupling gives gradient*[dn'] = h, where h=(partial_n g)(n0).
jump = sp.Symbol("jump_dn_prime", real=True)
normal_amplitude_jump = jump/(2*sp.sqrt(n0))
check("surface slope from normal density jump",
      hb**2*normal_amplitude_jump/(2*m*sp.sqrt(n0)), gradient*jump)

# Use positive placeholders for the subsonic stiffness and gradient energy.
k, a, area = sp.symbols("k a area", positive=True)
h, q = sp.symbols("h q", real=True)
length = sp.sqrt(a/k)
amplitude = sp.Symbol("amplitude", real=True)
right = amplitude*sp.exp(-z/length)
left = amplitude*sp.exp(z/length)
h_at_contact = a*(right.diff(z)-left.diff(z)).subs(z, 0)
scale, = sp.solve(h_at_contact-h, amplitude)
check("solution away from contact", -a*right.diff(z,2)+k*right, 0)
per_area = (sp.integrate(left,(z,-sp.oo,0)) +
            sp.integrate(right,(z,0,sp.oo))).subs(amplitude,scale)
check("stationary excess constituents per area", per_area, -h/k)
radius, energy = sp.symbols("R E", positive=True)
common_slope = sp.Function("H")(radius,energy)
odd_slope = sp.Function("Hodd")(radius,energy)
number = area*per_area
upper = number.subs(h,common_slope+odd_slope)
lower = number.subs(h,common_slope-odd_slope)
check("charge half-difference for a supplied odd contact",
      q*(upper-lower)/2, -q*area*odd_slope/k)
check("charge half-difference for the reflected scalar contact",
      (q*(upper-lower)/2).subs(odd_slope,0), 0)

# All constituents in this single transmitting stream cross the patch.
# Carry out the time integrals before comparing to the density inventory.
delay = (sp.integrate(left/j,(z,-sp.oo,0)) +
         sp.integrate(right/j,(z,0,sp.oo))).subs(amplitude,scale)
check("excess crossing time", delay, -h/(j*k))
check("excess number beyond changed transit inventory", number-area*j*delay, 0)
