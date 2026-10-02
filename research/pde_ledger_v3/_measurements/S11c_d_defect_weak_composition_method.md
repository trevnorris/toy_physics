# Global momentum bounds and weak retained composition: proposed method

2026-10-02. Author derivation for independent assessment; no new science has run.
The completed source/consumer inventory at 7455ee78 supplies actual ordered
addresses, both-face source joins, coefficient fields and response arguments.
This proposal addresses its explicit global test-space gap. It does not evaluate
any convolution, restore old producers, compute a finite inverse or start a
near-unity defect sweep. The conclusions below are proposed, not established by
writing this document or by the earlier compact-domain certificates.

## 1. Target and exclusions

Keep real omega=3, saved conserved edges (1/5,1/10), strict rest bulk,
LAB_HELD/RHO4_CONSTANT, W=1, L=10, rho=1/10, Lambda=1/100, tau=1/10 and the
same native coefficients. Keep eta and sigma independent in G={00,10,01,11}.
The effective cs family remains [1,2]; this is not primitive calibration.

The new domain extension is only in internal response momenta: k,l range over
all real numbers, without a projection. The proposed result is a continuous
bilinear weak action on Schwartz test fields in the one varying coordinate,
with a continuous real outgoing limit at modal speed matching. This is weaker
than a bounded inverse, an outgoing scattering solution, a current estimate,
or smooth/differentiable leakage. Incoming plane waves are not Schwartz. There
is no automatic extension to them, to the old finite basis or to cutoffs 4/6.
The untouched local slab part is outside this pressure-response proposition.

The saved physical source/consumer coefficients remain bound to real omega=3.
For proofs only, regularize the RESPONSE law by Omega=3+i delta, 0<delta<=1/10,
and take delta down to zero. The test fields and source/consumer coefficients
are held at real frequency. This auxiliary regularization is not a claim that
the physically composed operator has been analytically continued to complex
frequency, and does not reuse a conjugated power law. It selects the same
native response contact/PV prescription. No finite-regulator solve is proposed.

Restore the complete saved response definitions and actual new inventory.
Do not call old boundary, triangular solve, source, trace, grade or row functions.
New work consists only of global weight/test-space certificates and their joins.
All physical response formulas stay unchanged. A failed join or a distribution
product that this prescription cannot justify is a substantive finding.

## 2. Source fields, order and weak Fourier normalization

The inventory contains 34 distinct coefficient fields. Their displayed forms
are finite polynomials in T=tanh(x/10), with exact constant coefficients; this
must be certified for EVERY actual field, including zeros and consumers, not
inferred from a few controls. Save the original constructor, polynomial,
constant denominators and literal reconstruction residual. Require no leftover
x, unbound symbol or other function in a coefficient. Reuse the existing field
IDs and every address-to-field map. Do not recompute the source grades.

For a polynomial P(T), differentiation uses the exact recurrence

    P_0=P,   P_(n+1)=(1-T^2) P_n'(T)/10.

Since |T|<=1, every derivative is bounded by the finite sum of absolute
coefficients of P_n. This is an analytic induction, not an attempt to compute
infinitely many derivatives. Multiplication by each saved field therefore maps
Schwartz space S(R) continuously into itself. The existing wave-jet map gives
D_j=(-i*3)^n_t (i/5)^n_2 (i/10)^n_3 partial_x^n_1.
Persist each actual jet and order; no guessed maximum derivative order.
In particular the source is b(x) D_j u(x), NOT D_j[b(x)u(x)].

Use the already joined convention

    hat f(p)=(1/(2pi)) integral exp(-ipx) f(x) dx,
    f(x)=integral exp(ipx) hat f(p) dp.

For each inventory address with source coefficient b and consumer coefficient c,
put X(k)=hat[b D_j u](k), Y(l)=2pi hat[c v](-l), where u,v are Schwartz.
The pairing is complex bilinear (no conjugation implicit); an antilinear test
pairing would replace v by its conjugate explicitly. Neither pairing is power.
Define the addressed pressure contribution as

    B(v,u)=integral Y(l) F(l,k) X(k) dl dk,                 (1)

interpreted by the distribution prescription below. The native epsilon factor
is restored once as recorded by the address. Flat F includes delta(l-k), so
both its entire pole and its jet reduce to q(l)=q(k) on support.

This is the weak form of the recorded r-l, k-p (flat: l-p) coefficient convolutions;
source differentiation remains at the ORIGINAL wave momentum p before source
multiplication. Absorbing b and c into X,Y is a test-space construction, not a
change of their Fourier arguments and not evaluation of their distributional
Fourier transforms. General coefficient transforms must not be multiplied as
arbitrary distributions. The inverse-Fourier duality above supplies the precise
meaning. Constant coefficients give the saved delta rules.

For normal slots use the separately joined kernel N_f(l,k)=i*f*q(l)*F(l,k),
f=+1/-1. Do not assume q(l)Y(l) is a smooth test function at a branch point.
Prove bounds for the whole normal kernel, including its height-PV coefficient.

## 3. Global outgoing-depth and profile estimates

Use the saved unrestricted native laws

    q_delta(p)=sqrt(Omega^2/cs^2 - 1/20 - p^2),
    mu=rho*Omega, a=Lambda/[rho^2(1-i*Omega*tau)], beta=a*mu.

All depths and beta lie in the first quadrant; beta is strictly interior.
Reuse the source-joined constants b=3000/11101, |mu|<=2/5, |a|<=1,
|beta|<=2/5 and Re beta>=b. These do not depend on p. Derive the missing global
upper bound |q_delta(p)|<=|p|+4<=4(1+|p|), replacing the OLD compact bound 5.
For same-quadrant roots,

    |q+beta|>=b, |q/(q+beta)|<=1, |qi/(qm+qi)|<=1,
    |q(p)-q(k)|<=sqrt(|p-k|*|p+k|).

Let a_* = sqrt(879)/20 and kappa_delta=sqrt((9-delta^2)/cs^2-1/20).
The saved lower-root certificate gives a_*<=kappa_delta<=3. For all real p,

    1/|q_delta(p)| <= a_*^(-1/2) sum_(s=+-1) |p-s*kappa_delta|^(-1/2). (2)

At zero denominators this is an almost-everywhere envelope, not a finite
pointwise value. At delta>0 the real comparison endpoints are not called
complex branch points. Endpoint motion and coincidences retain a SUM of simple
inverse-square-root envelopes, not a product.

Reuse the actual tanh factor A(s)=L*s/[4*sinh(pi*L*s/2)] with removable
A(0)=1/(2pi), h(s)=(W/2)[delta(s)/2+PV(A(s)/(i*s))], j(s)=L*A(s)/2.
The existing global estimates A<=1/(2pi), |A'|<=L/4 and |j'|<=L^2/8 imply
bounded j and its derivative. One deliberately loose NEW useful envelope is

    |A(s)|<=11 exp(-|s|),
    |A(t)A(Q-t)|<=121 exp(-|t|).                          (3)

For |s|<=1 use A<=1/(2pi)<=11/e; for |s|>=1 use the inherited
A<=10|s| exp(-5pi|s|) and |s|exp(-(5pi-1)|s|)<=1. Check these constants and
actual L join. No large-Q exponential prefactor is permitted in the final bound.

A shift-uniform singular integral estimate closes the tail issue. With
w(t)=(1+|t|)^2 exp(-|t|), sup w<=2 and integral w=10. Splitting at |t-z|=1,

    integral w(t)|t-z|^(-1/2) dt <= 4 sup(w)+integral(w) <=18,
    integral w(t)/|q_delta(k+t)| dt <= Cq=36/sqrt(a_*).     (4)

The reflected route q_delta(l-t) has the same bound, uniformly in every k,l.
These are analytic bounds on integrals, not numerical evaluations of the
response convolutions. The later instrument checks algebra/constants/routes;
it does not claim a symbolic program proves measure theory.

## 4. Whole mixed responses and polynomial growth

Restore the saved reference kernels (real specialization and unrestricted law):

    F0=mu/(qi+beta) delta(l-k),
    Fh=-i*mu*qi/(qo+beta) h(l-k),
    Fj=mu*k/[(qo+beta)(qi+beta)] j(l-k),
    F11_iter=C(l,k) H(l-k) + integral J(t;l,k) dt,
    C=-i*mu*k*qo/[(qo+beta)(qi+beta)].

H is the inherited contact plus subtracted PV whole convolution. Its existing
bounds on j,j' and A give |H(Q)|<=100 for ALL Q: on |s|<=1 use the derivative
bound; on |s|>1 use |A(s)/s|<=10 exp(-5pi|s|) and |j(Q-s)-j(Q)|<=5/pi.
Keep the nonzero contact (W/4)j(Q). The old restriction |Q|<=6 is unnecessary
for these particular estimates; this extension must be assessed and certified.

The inherited remaining iteration density, m=k+t, is

    J=(a*mu^2*W*L/4) A(t)A(l-k-t)
       * m(m+k)*qi/[qm*(qo+beta)*(qm+beta)*(qi+beta)*(qm+qi)].

With P=1+|k|+|l|, use |m(m+k)|<=2 P^2(1+|t|)^2 and the quadrant ratio
qi/(qm+qi). From (3)-(4),

    integral |J| dt <= KJ P^2,
    KJ=(4/5)*121*Cq/b^3.                                 (5)

The separately saved complete direct density is

    GD=(W*L/(4i)) A(t)A(l-k-t) Bc,
    Bc=-i*mu/[(qi+beta)(qo+beta)] *
       [ k(2l-t)/(qs+qo) + k(t+2k)*qi/[qh*(qh+qi)] + qi^2/qh ],
    qh=q(k+t), qs=q(l-t).

The four-depth common sheet and both native faces must still join exactly.
Do not replace qs by qh. Use |qi|^2<=16 P^2, the quadrant ratio and
1/|qs+qo|<=1/|qs| to obtain

    integral |GD| dt <= KD P^2,
    KD=18*121*(2*Cq)/b^2.                                (6)

In (6), |mu|*W*L/4<=1 and each polynomial numerator is bounded by
18 P^2(1+|t|); (4) covers both routes. These constants are deliberately loose
existence bounds, not numerical error bars. Check all factors independently
against the saved density, not against a newly typed substitute alone.
The native normal multiplier has modulus <=4P, so the corresponding J/direct
normal bounds have growth P^3. C*H grows at most linearly, its normal jet at
most quadratically; the slope kernel has the same or smaller polynomial growth.
The diagonal flat pressure and normal multipliers are bounded.

For external parameters in ANY fixed compact k,l set, the same moving-endpoint
uniform-integrability argument as before, with that compact set's constants,
gives continuity of J and direct whole integrals in (k,l,cs,delta), including
simultaneous external grazing and coincident internal endpoints. To prove this,
retain (i) uniform absolute continuity on bounded t intervals, (ii) a uniform
tail for that compact external set, and (iii) almost-everywhere convergence.
A global pointwise dominating function is not assumed. Equations (5)-(6) then
give the polynomial bound needed to pair with Schwartz X,Y over ALL momenta.
No assertion of differentiability at matching follows.

The inherited Dwhole is a complete closed convolution. It enters grade11 ONCE;
no new external resolvent, reference inverse or second middle integration.
Bare Rprod, factored E and Dwhole remain distinct. H and J retain their own
actual signatures, transfer variables and native first-shape iteration once.

## 5. Global height-PV action, including the normal jet

For fixed k define B(l,k)=-i*mu*qi/(qo+beta),
f_k(Q)=B(k+Q,k)Y(k+Q), and chi(Q)=1 on |Q|<=1, zero otherwise. The inherited
prescription for the inner height action is

    (W/4)f_k(0) + W/(2i) integral A(Q)[f_k(Q)-f_k(0)chi(Q)]/Q dQ. (7)

Use the symmetric subtraction, with PV integral A(Q)chi(Q)/Q=0. Its cutoff is
a local subtraction convention, NOT a momentum cutoff on the operator.
On |Q|<=1, the global root-difference bound implies, with Pk=1+|k|,

    |B| <= (8/(5b))*Pk,
    |B(k+Q,k)-B(k,k)| <= (8*sqrt(3)/(5b^2))*Pk^(3/2)*sqrt(|Q|).

The product with Y adds a term (8/(5b))*Pk*||Y'||_infinity*|Q|. After division
by Q, these bounds are integrable. Outside |Q|<=1 the exponential A/Q tail
and bounded Y suffice. Thus the absolute subtracted action is bounded by
C*Pk^(3/2)*(||Y||_infinity+||Y'||_infinity), uniformly in cs and delta.
Multiplication by Schwartz X(k) makes the outer integral absolutely convergent.

For the normal slot, Bn=i*f*qo*B=f*mu*qi*qo/(qo+beta). Its magnitude is
<= (8/5)Pk and its variation is bounded by
(16*sqrt(3)/(25b^2))*Pk^(3/2)*sqrt(|Q|), using |beta|<=2/5.
Use this direct coefficient bound; never argue that q(l)Y is C1 at a root.
The same subtracted-action proof applies with the actual lower sign.

Pointwise convergence of the subtracted numerator for Q!=0, these local bounds,
the exponential Q tail, and Schwartz decay in k give the outgoing distributional
limit by dominated convergence. The contact is included before taking limits.
The regularized response is only an auxiliary proof family; no physical
positive-regulator composed source law or current is asserted.

## 6. Proposed composition conclusion and bounded instrument

Equations (1)-(7), if supported, define each retained pressure/normal block as
a continuous bilinear form on S(R)xS(R). Polynomial bounds of order at most
three for ordinary whole kernels and the height estimate supply a finite
Schwartz-seminorm bound. The actual native source derivatives are finite order,
and multiplication by every certified coefficient preserves S. The sum over
the finite saved inventory is consequently a continuous weak retained pressure
contribution, with weak continuity in cs through both selected matches.
This is a distribution/test-space conclusion, not a physical field value,
coercivity, positivity, invertibility, an operator-norm scattering bound or a
claim that loss varies smoothly. The complete physical model's drain remains off.

After substantive independent assessment, a task-local bounded certificate
instrument should restore only the required saved operands and identities. It
must not rerun the source inventory or response derivations. Required NEW outputs:

- A certificate for all 34 actual coefficient fields, polynomial-in-tanh
  reconstruction, derivative recurrence and actual source-jet/address coverage.
- Source-joined global q/beta/profile/shift-uniform envelope constants, actual
  J/GD numerator estimates and separate pressure/normal polynomial degrees.
- Full contact/PV and whole-tag argument joins, both faces, native signs and
  independent eta/sigma order. No hidden source/consumer grade collapse.
- Source-bound weak Fourier-duality/ordering and normal-height certificates;
  the universal functional analysis is an assessed analytic argument, not a
  finite set of symbolic samples promoted to a theorem.
- Responsive addressed controls: interchange source differentiation and profile
  multiplication on an ACTUAL nonconstant source field and jet; use q(r) in
  place of q(l) through a normal-slot address with a nonconstant consumer;
  remove the surviving H contact; corrupt the reflected direct root route in
  the actual physical dispersion identity at applicable nonzero momenta.
  Reuse prior response controls where identical; derive only genuinely missing
  global/composition controls. Persist operands and applicability before guards.

Do not evaluate the response integral to make the controls. Formal coefficient
or algebraic-identity movement must be labelled as such, not a loss signal.
Unknown or nonzero exact joins stop with saved evidence; no numerical fallback,
relaxed tolerance or automatic retry. At completion inspect actual resources,
all operands/returns, copies and hashes. Review clearance is not result acceptance.

The scope of this method ends at global weak applicability of the retained
pressure contribution. A finite scattering method would still need a justified
trial/test space and asymptotic excitation, numerical action/endpoint controls,
current and finite-inverse sensitivity, regulator policy and applicable reviews.
No finite solve, full Green/FORM/A11/A12, calibrated light claim or defect sweep
is authorized by this proposal alone.
