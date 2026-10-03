# Packet-action preflight: adapters and analytic truncation plan

2026-10-03. Concrete prerequisite of the wavelength-resolved two-packet method.
No scientific payload has been restored during preparation. There is no READY
gate, launch or computed packet action. This build is not the future numerical
evaluator. The new analytic bounds below require independent assessment and
actual guarded operand joins before they can be used as truncation certificates.

The method remains the two real-frequency-3 e_W-input/THETA-test complex-bilinear
Gaussian pairings, with s=8, centres -5/2 and +5/2, carriers sqrt(595)/10 and 0,
cs=sqrt(6)/2, rest bulk, LAB_HELD/RHO4_CONSTANT and tangents 1/5,1/10. The matching
speed comes from the saved LEFT point, not a primitive calibration. The optional
finite-origin display is not evaluated. No response integral, Fourier transform,
local action, finite matrix, inverse, current or loss is computed by this worker.

## Restored inputs and new adapters

The manifest pins 222 complete JSON files, including the full original THETA
address row, all 400 local cells, the reviewed 544/16 projections, all 34 field
polynomial inputs and zero returns, 17 distinct inherited complete-factor proof
sets, the whole definitions, global bounds, physical context and saved match.
The old producer's `operandSha256` is a tuple digest; it is not confused with a
file receipt. Byte receipts are separate. Copies are byte-identical and indexed
incrementally. All original files are checked again after the guarded run.

Selection uses the actual row/field/jet metadata and joins the full originals
against the reviewed projections. It retains 102 formal, 106 zero-source and
336 zero-consumer addresses. The 16 ordered pressure grade triples are different
from the 16 local derivative-order/grade cells. Epsilon is recorded once on
formal entries and zero on exact-zero entries. All 48 selected normal-height
addresses have zero consumer. No numerical smallness removes an address.

For pressure fields, the old polynomial numerator coefficients are divided by
the saved constant denominator and compared to the saved quotient coefficients.
Their original field and physical-tanh reconstruction operands join the literal
old cancelled-zero return; that proof is restored, not repeated. Local records
already store quotient coefficient vectors. Their source sums, physical-tanh
proofs and derivatives remain inherited. No old summation or derivative is
recomputed. The new coefficient vector conversion checks exact complex rational
entries and calculates an explicit real/imaginary coefficient norm.

The actual source context fixes real omega3 and the native tangents; the original
input still says omega1/cs10. Only the previously declared frequency/speed use is
carried forward. Source coefficients are not rebound to an auxiliary regulator.
Native profile L=10 and raw absence of source/consumer speed symbols are restored.
The preflight inherits the source and four THETA consumer unit joins and the
1/(2pi) Fourier normalization. It explicitly does not claim a new post-binding
unit derivation for every local and pressure summand. That remains a requirement
of the future evaluator before any action assembly.

Every scalar address joins its actual complete mapped factor and normal sign.
Twenty new numerical templates cover both faces, both slots and all five
components. Each function call is checked against the full original signature:
q(k), q(l), h(l-k), j(l-k), Hwhole(l-k,1,10), and the face-specific Jwhole/Dwhole
(l,k,3,cs,1/5,1/10,1,10). Whole values remain independent formal tags. Flat support
uses output q(l) in the entire pole. There is no added resolvent and no second
middle integral of a whole direct term. Reuse of a template requires actual
identical complete factors, not an ID alone.

The new J and D numerical density formulae join the actual saved definitions at
omega3. Their unrestricted depth symbols remain independent during the identity
check: qm=q(k+t), qh=q(k+t), qs=q(l-t), qi=q(k), qo=q(l). In particular qs is not
replaced by qh. H's contact and subtracted integrand join separately. Physical
h=(1+tanh(x/10))/4 and j=(1-tanh(x/10)^2)/4 have h'=j/10. The Fourier derivative,
sech-squared transform and product/convolution identity are assessed analytic
lemmas, not symbolic numerical integrations. The stored removable A(0)=1/(2pi)
joins the same A(t)=5t/(2sinh(5pi t)). A future evaluator must implement that
removable value without sampling an invented 0/0 value.

## Rational contour envelopes

For any coefficient polynomial P(T), let N(P) be the sum over its coefficients
of |Re c|+|Im c|. All actual coefficients must be complex rational. On |Im z|<=5,
|tanh(z/10)|<=1, so |P(tanh(z/10))|<=N(P). This follows from 5/10<pi/4; the closest
poles are at |Im z|=5pi.

Write y=|Re z-x_u|, c=5 and use |p0|<=3 for both declared carriers. Nonnegative
rational polynomials M_n majorize the carrier-inclusive input derivative:

    M_0(y)=1,
    M_(n+1)(y)=M_n'(y)+(3+(y+5)/64) M_n(y).

This is the coefficient triangle bound on the actual Gaussian derivative
recurrence, not a replacement of the actual numerical derivative. The native
jet magnitude also contains 3^nt/5^n2/10^n3, with spatial derivatives at the
original source wave momentum before coefficient multiplication.

For two-sided absolute Gaussian moments at s=8 use rational upper bounds
I_0<3s, I_1=2s^2 and I_n=(n-1)s^2 I_(n-2). Let I(P) be the resulting positive
linear functional on a nonnegative polynomial P(y). Since exp(25/128)<2,
1/(2pi)<1/6 and sqrt(2pi)<3, the following constants safely dominate:

    CX = nativeJetMagnitude * N(source) * I(M_n)/3,
    CY = 2 * N(consumer) * I(1),
    CYprime = 2 * N(consumer) * I(y+5/2+5).

Thus |X(k)|<=CX exp(-5|k-p0|) and |Y(l)|<=CY exp(-5|l-p0|).
CYprime bounds the Fourier derivative Y' with the same decay; |z| is bounded
by y+|x_v|+5 on the shifted test contour. Y includes its required factor 2pi.
The signed contour offsets are k-p0 for X and p0-l for Y in the method. The
bound uses their absolute values and does not change those signs in evaluation.
These envelopes are not assertions that a nonconstant coefficient times a
Gaussian has a pure Gaussian Fourier transform.

## Source-bound pressure truncation budget

Use the already certified global bounds with b=3000/11101 and
kappa_min=sqrt(879)/20>1. Cq_saved=36/sqrt(kappa_min)<36. Exact new scalar joins
compare the saved beta law, profile, shift-root bound, J/D constants and both
height constants to these same operands. Put

    KJ=(4/5)*121*36/b^3,
    KD=18*121*72/b^2,
    C=1+4*KJ+4*KD+400/b^2,
    A0=8/(5b), A1=16/(5b^2).

Here C is a coarse common ordinary-kernel bound C*(1+|k|+|l|)^3, including normal
slots, flat, slope, H and the whole J/D terms. H<=100 is inherited; its C(l,k)
factor and normal factor are included. A0 bounds the pressure height magnitude
by A0*(1+|k|). A1 is an upper bound on the saved Holder constant after using
sqrt(3)<2. The selected live height slots must actually all be pressure slots.
A future expansion to live normal height is not silently covered by this test.

Define F_n=integral_R (1+|x|)^n exp(-5|x|) dx, and E_n(K) the same integral
outside |x|<=K. They have finite positive polynomial formulae. The implementation
uses exact rational integration-by-parts coefficients and replaces exp(-rK) by
2^(-rK), for positive integer r and nonnegative integer K. This only enlarges
the tail. Bounds exp(15)<3^15 and exp(30)<3^30 absorb the two |p0|<=3 offsets.
There is no floating estimate of a claimed analytic inequality.

For each live ordinary off-diagonal address, the outer square tail is at most

    2*C*CX*CY*3^30*F_3*E_3(K).

P^3<=(1+|k|)^3(1+|l|)^3 and the union bound over the two tails proves this; overlap
is harmless overcount. On a flat address only the diagonal integral remains:
C*CX*CY*3^30*E_3(K) is safe after dropping one decaying exponential.

For height use the method's positive-half-line paired PV, with the contact
included. Set

    Ch=(A0/4+11*A0)*CY + (2*A1*CY+A0*CYprime)/6.

For 0<Q<=1, |f(Q)-f(-Q)| is bounded by
2*A1*Pk^(3/2)*CY*sqrt(Q)+2*A0*Pk*CYprime*Q. With A<=1/(2pi)<1/6 and W/2=1/2,
the near contribution is the last term in Ch. For Q>1, A<=11 exp(-Q) and the
two-sided difference give the 11*A0 term; W/4 contact gives A0/4. Therefore the
height inner action is bounded by Ch*Pk^2. Its outer tail is

    Ch*CX*3^15*E_2(K).

The separate discarded Q>U tail, U>=1, is bounded over all k by
11*A0*CX*CY*3^15*F_1*2^(-U). No square-domain tail is substituted for this
(k,Q) domain, and no cancellation is used to pay the budget.

For middle truncation require T>=K+4 with kappa<3. Every endpoint of q(k+t) and
q(l-t) lies strictly inside |t|<T, and the depth moduli outside exceed one.
The saved numerator estimates then imply density tails before normal factors:

    J: (4/5)*121/b^3 * P^2 * (1+|t|)^2 exp(-|t|),
    D: 36*121/b^2 * P^2 * (1+|t|) exp(-|t|).

The D bound treats the SUM of reflected/height terms; it introduces no product
of roots. Multiply each by 4P (safe for pressure as well as normal slots) and
by the packet envelopes. Integrating k/l over all real momenta gives the factor
4*CX*CY*3^30*F_3^2. The remaining middle tails use the rational formula above
with rate1 and degree2 or1. H's subtracted middle tail is bounded by
(55/3)*exp(-T): A<=11 exp(-|t|), |j(Q-t)-j(Q)|<=5/pi, W/2=1/2 and T>=1.
Its outer C/normal multiplier is bounded by (4/b)P^2. This gives the additional
middle tail (4/b)*CX*CY*3^30*F_2^2*(55/3)*2^(-T) for native iteration addresses.

The worker sums POSITIVE per-address outer, height-Q and middle tails. Starting
at K=4,U=1,T=8, it increases only radii whose allotted total exceeds one third
of 1e-11, always maintaining T>=K+4. Every iteration is saved before a capacity
check. Predeclared capacity K<=256,U<=512,T<=512 records impracticality if missed.
These are spatial/momentum truncation capacities, not time limits, and no
quadrature or post-failure retry occurs. If the sum over all 102 live addresses
is below 1e-11, each unweighted grade/component is below that budget by
positivity. The bound covers both carriers independently. Local x-integral
truncation is a separate future obligation and is not claimed here.

## Execution and what remains

Twenty-two standard-library tests check metadata, original/projection identity,
all byte receipts, argument schemas, fail-closed gate/invocation, exact-true
predicate behavior, evidence persistence, containment order and launch route.
They do not restore a scientific expression or check the mathematics numerically.
All mathematical imports and payload restoration happen only after the pinned
gate and the unchanged shared guard/supervisor containment check. Resources are
4GiB native/cgroup in the16GiB pool, zero swap, one CPU/thread,32 tasks,4GiB host
reserve, desktop priority and no deadline. The hook is armed before launch.

The four numerical controls are only preselected from live address metadata.
This build cannot declare them responsive. A future evaluator must still join
all units; construct actual Fourier derivative/radius/panel arrays; pass both
independent transform routes; implement every collision cell, PV prescription
and branch limit; allocate local tails; evaluate both complete action routes;
and run the four numerical controls. Preflight success is not quadrature
readiness. New exact adapter failures remain evidence and stop the run; no
unknown is turned into true, no tolerance is relaxed, and no automatic retry
or source/producer replay is enabled.
