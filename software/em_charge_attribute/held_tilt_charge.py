"""Moving versus held aperture: relative carrier flux, storage and signed count."""
import sympy as sp


def check(label, actual, expected):
    value = sp.factor(actual)
    print(f"{label} = {value}", flush=True)
    if sp.simplify(value - expected) != 0:
        raise AssertionError(label)


# Orthonormal (radial, polar, w) components. The azimuthal direction is
# perpendicular to all three. The axisymmetric graph can have any tilt rate.
radial, polar, ht = sp.symbols("h_r h_theta_over_r h_t", real=True)
r, R, n0, intake = sp.symbols("r R n0 D", positive=True)
q, sign = sp.symbols("q s", real=True)
tr = sp.Matrix([1, 0, radial])
tang = sp.Matrix([0, 1, polar])
orthogonal = tr - tang * tr.dot(tang)/tang.dot(tang)
w = sp.sqrt(1 + radial**2 + polar**2)
# |orthogonal| = W / |tang|; area per unit solid angle is r^2 |tang|.
check("squared norm of cross-section conormal",
      orthogonal.dot(orthogonal), w**2/tang.dot(tang))
directed_area = r**2 * tang.dot(tang) * orthogonal / w
density = n0/w
vr = intake/(4*sp.pi*n0*r**2)
cut_velocity = sp.Matrix([0, 0, ht])
moving_velocity = cut_velocity + vr*tr
held_velocity = vr*tr
moving_flux = q*density*(moving_velocity-cut_velocity).dot(directed_area)
held_flux = q*density*held_velocity.dot(directed_area)
check("moving charge crossing per solid angle", moving_flux, q*intake/(4*sp.pi))
check("held charge crossing per solid angle", held_flux, q*intake/(4*sp.pi))
check("unwinding minus held integrated crossing rate", 4*sp.pi*(moving_flux-held_flux), 0)
check("uncorrected sweep term per solid angle",
      q*density*moving_velocity.dot(directed_area)-held_flux, q*n0*r**2*ht*radial/w**2)
stored = sp.integrate(4*sp.pi*r**2*q*density*w, (r, 2*R, 3*R))
check("carrier charge in conversion layer, either branch", stored, 76*sp.pi*q*n0*R**3/3)

# The geometric completion moves a weighted point. In units R=1 its
# oriented flux at a cut is A(beta)*beta_dot*delta(cut-beta). Integrate
# against monotone d beta and keep the bulk edge as a one-sided limit.
beta, u = sp.symbols("beta u", real=True)
lost = sp.integrate(30*u**2*(1-u)**2, (u, 0, beta-2))
for c in (sp.Rational(5, 2), sp.Rational(8, 3)):
    val = sp.integrate(sign*lost*sp.DiracDelta(beta-c), (beta, 2, 3))
    check(f"completed geometric dose at cut r/R={c}", val, sign*lost.subs(beta, c))
check("completed geometric dose approaching edge", sp.limit(sign*lost, beta, 3, dir="-"), sign)

# Cumulative accounting after a putative charge Q arrives at the edge.
# Released fraction p and propagation time tau are unspecified. This tests
# what two separate probes would read for different allowed continuations.
t, tau, received = sp.symbols("t tau Q", real=True)
p = sp.Function("p")
at_edge = received*(1-p(t))
in_flight = received*(p(t)-p(t-tau))
past_cut = received*p(t-tau)
check("edge, in-flight and downstream total", at_edge+in_flight+past_cut, received)
for final in (0, 1):
    print(f"conditional late (edge, flight, passed cut), released fraction {final}: "
          f"{tuple(a.subs({p(t): final, p(t-tau): final}) for a in (at_edge,in_flight,past_cut))}")

# Continuity permits a transient number inventory even with fixed inlet.
# This control uses the normalized section coordinate x in [0,1].
x = sp.symbols("x", real=True)
stock = sp.Function("B")(t)
number_per_x = sp.diff(x*stock,x)
number_flux = -sp.diff(x*stock,t)
check("transit control continuity", number_per_x.diff(t)+number_flux.diff(x), 0)
check("transit control extra inlet", number_flux.subs(x,0), 0)
check("transit control extra outlet charge rate", q*number_flux.subs(x,1), -q*stock.diff(t))
check("transit control extra stored charge", q*sp.integrate(number_per_x,(x,0,1)), q*stock)
