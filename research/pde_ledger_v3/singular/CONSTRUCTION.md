# S15 Singular construction

**Status: unfinished and unreviewed.** The trials were stopped on the user's
2026-11-10 instruction. This describes the attempted construction, not an
accepted build; see the [unfinished-build report](../directives/_measurements/S15_separation_departure_singular_run_record.md).

The inputs are the supplied equations and labelled definitions in the S15 and
S15-reading decision lists. No other engine's source, exports or output is used.
The Python program constructs labelled density/velocity products, integrates
them, and forms the geometric separation derivative. Singular alone computes
the conditions that its output is identically zero. This document specifies the
method, not a result or a comparison.

For each surface the program constructs all seven nonempty label sets for the
carried current, centre transport, radius transport, viscosity and first-order
pressure expansion. For the region it constructs the seven momentum groups and
their time rates. It applies the scope rule to the unintegrated fields, so
empty groups are retained. Whole pieces and the named rest/transport/crossing
sums contain general fields and are outside this item-4 directive. Pressure is
never replaced by its first-order expansion. Viscosity's active single-label
groups contain the general viscosities; its other groups are empty by the
linearity of the supplied stress in velocity derivatives.

Write `h=D/d`, `x=X2+r*e`. Azimuthal tensor moments of ranks zero through three
reduce every in-scope integrand to scalar integrals of
`mu^j*(d^2-2*d*r*mu+r^2)^(s/2)`. Substituting the distance for `mu` gives finite
sums of `H(p,r)`, the primitive of `u^p` between `d-r` and `d+r`. Radial
integration gives `K(t,p)`, the integral of `r^t*H(p,r)`, and the elementary
moment `T(t)`. These are printed as finite closed expressions using real powers,
logs, Appell F1 and the dilogarithm. Euler's Appell integral is used on its real
positive-base continuation; a Taylor truncation is never used. The primitive
cases `p+1=0` and `t+1=0` are printed explicitly. They are cases needed to form
the exact objects, not assumed generic exponent restrictions.

The rate of the region momentum differentiates amplitudes and separation, with
`Ddot=dot_X1-dot_X2`, and includes both radius endpoint terms. Exponents are
constant in time. A sufficient C1 condition for the time-dependent data and the
compact-domain analytic differentiation premise accompany the objects.

Introduce finite state variables for each required H, K and endpoint power U,
including their derivative dependencies. Their d derivatives satisfy a linear
system with denominators only in `d` and `d +/- radius`. Radius powers and T are
held constants with explicit definitions. The program forms `E=d*Q_d+2*Q` and
clears only nonzero geometric denominators and the common nonzero factor pi.

Singular reduces the polynomial in h modulo `h1^2+h2^2+h3^2-1`, then extracts
all its direction coefficients. Each coefficient is linear in the analytic
states. It is represented by a row in a free module with basis `(1,y_0,...)`,
over the polynomial ring in the other quantities. The source checks that the
rows contain no state variables as coefficients and reconstruct their outputs.
The generated ring also names the states for input/output, but no module
coefficient contains them. Its coefficient field is Q.

Write the printed state system as `delta(y)=M*y`, with a constant coordinate
adjoined if needed, where `delta=P(d)*partial_d`. Python forms P as the least
common multiple of the actual state-system denominators for this object; its
factors are nonzero geometric quantities on the domain. Singular computes the
row connection `nabla(v)=P*partial_d(v)+v*M`. It repeatedly computes a module
standard basis and adjoins connection remainders. A whole row may be divided
by powers of `d,R,rin,d-R,d+R,d-rin,d+rin` only when every entry is divisible.
These factors are nonzero throughout the domain. Divisions and remainders are
printed. Normalized basis rows are also adjoined until no further row is added.
There is no round or time limit.

Exponents, amplitudes, rates, velocities, radii and lifted constants are
**polynomial variables**, not rational-function field parameters. No generic
parameter division or exceptional-exponent exclusion occurs. A finite free
module over this Noetherian ring is Noetherian, so the ascending module chain
stabilizes. Each adjoined row is a polynomial combination of derivatives of the
original outputs, allowing only inverses of the printed nonzero geometric
factors; hence an identically vanishing original output satisfies it. At
closure, connection derivatives of basis rows are polynomial combinations of
basis rows. Their values along the analytic state trajectory satisfy a finite
homogeneous linear ODE, regular on the domain because P is nonzero. Zero initial
values at `d=2*R` therefore imply vanishing throughout the connected domain.
This gives necessity and sufficiency, including parameter specializations.
Real-sphere polynomial density supplies the direction quantifier. These
uncomputed mathematical facts are printed as premises at their use sites.

The previous polynomial-ideal construction is retained as
`differential_ideal.sing` to identify attempt01. The current executable uses
`differential_module.sing`. The module method preserves linearity in the analytic
states throughout the closure computation; nonlinear products of those states
are not introduced by its standard-basis operations.

The conditions are a finite conjunction of equations in the named quantities,
using the printed special functions and all their primitive cases. They are
not a list of numerically fitted roots or a generic-exponent answer. After
Singular computes them, Python only renders its equations, substituting the
explicit state values at `d=2*R`. Under those conditions Q is printed as
`(2*R/d)^2*Q(2*R,h)`, with its reference value expanded in named quantities.
No other engine is compared, and no force interpretation is made.
