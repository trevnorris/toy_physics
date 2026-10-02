# Remaining retained face response at grazing: proposed bounded method

2026-10-02. This is a source-based author derivation for independent assessment,
not an executed identity, a cleared instrument, an integral value or a defect
result. The completed direct-kernel certificate is preserved at 30b0fbaa. This
proposal addresses the remaining zero/first-order face response and the existing
mixed first-shape iteration, including its reference-trace conversion. It does
not repeat the direct boundary construction or silently extend its conclusion
to the full slab operator.

## 1. Fixed scope and the objects that must join

Keep strict-rest-bulk LAB_HELD/RHO4_CONSTANT, omega=3, saved edge momenta
(1/5,1/10), W=1, L=10, rho=1/10, Lambda=1/100 and tau=1/10. Keep eta and
sigma_W independent through the retained rectangle (0,0), (1,0), (0,1), (1,1).
There is no eta^2 or sigma_W^2 result. The effective speed cs is in [1,2]; input
and output profile momenta k,l are in [-3,3]. This is a mathematical compact
domain containing the selected modal matches, not a primitive calibration or
a frequency/defect sweep. Physical permeability and memory remain present.

Use Omega=3+i delta, 0<delta<=1/10, and then the outgoing limit delta->0+:

    q(p) = sqrt(Omega^2/cs^2 - 1/20 - p^2),
    mu = rho Omega,
    a = Lambda/[rho^2(1-i Omega tau)],
    beta = a mu,                 R(p) = q(p)/(q(p)+beta).

The square root lies in the first quadrant; at real frequency it is the
positive real or positive imaginary outgoing root. A zero depth is a limiting
value, never a pointwise division by zero. This is the same analytic
limiting-absorption prescription as the direct-kernel certificate; it does not
remove physical damping or analytically continue a conjugated power form.

The saved law gives beta=(30+9i)/109 at delta=0. Reuse, with actual symbol and
source joins, the certificates

    |q(p)+beta| >= b0 = 3000/11101,
    |a| <= 1, |mu| <= 2/5, |q(k)|,|q(l)| <= 5,
    kappa_delta = sqrt((9-delta^2)/cs^2-1/20) >= sqrt(879)/20.

Only the latter two elementary upper bounds are new if not supplied by the
saved certificates. No kappa=0 or beta=0/impermeable limit is covered.

The native sources are c1 `shape_source`, `dtn_first_kernel`, c2 `kernel_bridge`,
`reference_pressure_kernels`, `build_face`, `kernel_apply`, and the saved
both-face closure and trace operands. Inspect the actual sources in the packet.
At the lower face the lab height and lab normal multiplier both change sign:
their product in the trace operator has the same sign as at the upper face.
Do not also replace the outgoing depth by its negative. The final reference
normal jet remains i f q(l) times reference pressure, f=+1 or -1.

The scope is the response kernel per unit native source, before source-field
and slab-row composition. It supplies no new arbitrary/lower-face Fourier
source map, incoming transverse mode or complete retained slab operator. A
finite formal jet multiplier is not a physical source state.

## 2. Native first shape, closure and trace algebra

With Q=l-k, use the previously joined normalized one-dimensional transforms

    A(s)=L s/[4 sinh(pi L s/2)],        A(0)=1/(2 pi),
    h(s)=(W/2)[delta(s)/2 + PV(A(s)/(i s))],
    j(s)=(L/2) A(s).

Here h and j have their eta and sigma_W factors stripped separately. Do not
bind sigma_W=eta W/L before grade selection. The saved edge deltas remove the
two conserved edge integrations. All remaining compositions use plain dm; do
not add another 2 pi factor or a second full-weight transfer assignment.

The proposed one-dimensional forms of native c1 are

    Z0(k)=mu/q(k),
    Zh(l,k)=i mu [q(l)-q(k)]/q(l) h(Q),
    Zj(l,k)=mu k/[q(l)q(k)] j(Q).                         (1)

They follow from the actual on-shell first-shape numerator, including the
unchanged edge momentum contribution. They must join both native face operands
and the independently saved lower first-shape joins; equality on independent
formal depths is not a substitute for the physical dispersion join.

Native closure is (I+a Z)^(-1) Z. The physical first-order response is Rl Z1 Rk.
The physical mixed first-shape iteration is, with m the middle momentum,

    -a Rl [Zh(l,m) Rm Zj(m,k) + Zj(l,m) Rm Zh(m,k)] Rk.    (2)

There is no additional symmetrization. The separate saved direct mixed kernel
is added once, outside (2). Its already-whole convolution must never be passed
through a second middle integration.

The native reference trace is T=I+T1 with

    T1(l,k)=i q(k) h(Q).                                  (3)

This uses lab height f h times native lab normal i f q(k), for each face. In
the retained rectangle, converting a physical response to reference pressure
subtracts T1 F0 at grade (1,0) and T1 Fj at grade (1,1). The T1^2 term is eta^2,
outside this rectangle. Keep that exclusion explicit rather than asserting a
full second-order trace inverse.

The proposed reference kernels at the first three grades are

    F0(k)=mu/[q(k)+beta],                                 (4)
    Fh(l,k)=-i mu q(k)/[q(l)+beta] h(Q),                   (5)
    Fj(l,k)=mu k/[(q(l)+beta)(q(k)+beta)] j(Q).             (6)

F0 is a diagonal multiplier, with its delta understood. `build_face` first
constructs PHYSICAL pressure from `pmat`, then uses the trace-derived reference
normal jet in the affine reference-value solve. Its final slot must join (4–6)
and the mixed formula below after projection. Do not apply the trace inverse
twice, or merely check an intermediate matrix with the final native slot
unexamined. No producer or completed triangular solve is to be replayed.

## 3. Mixed iteration and the surviving trace contact

Let qi=q(k), qm=q(m), qo=q(l), and

    D3=(qo+beta)(qm+beta)(qi+beta).

Before combining terms, preserve the three native mixed contributions:

    physical height-left/slope-right:
      -i a mu^2 k (qo-qm)/D3 h(l-m) j(m-k),
    physical slope-left/height-right:
      -i a mu^2 m qi (qm-qi)/(qm D3) j(l-m) h(m-k),
    reference-trace subtraction:
      -i mu k qm/[(qm+beta)(qi+beta)] h(l-m) j(m-k).        (7)

The proposed exact combination of the first and third terms is

    C(l,k) h(l-m) j(m-k),
    C(l,k)=-i mu k qo/[(qo+beta)(qi+beta)].                 (8)

The middle-depth factor cancels, but the height contact DOES NOT generically
cancel. Its value is C(l,k) (W/4) j(Q). Zero k or qo can remove this particular
contact; neither is imposed for the general kernel. The direct kernel's zero
contact result must not be transferred to this expression.

The whole convolution H(Q)=(h*j)(Q) can be represented, without evaluating it,
by a contact plus an ordinary absolutely convergent subtracted integral:

    H(Q) = (W/4) j(Q)
         + W/(2i) integral_R A(s)[j(Q-s)-j(Q)]/s ds.       (9)

At finite delta this is exactly the original symmetric PV distribution,
because A is even and PV integral A(s)/s ds=0. Preserve its original contact/PV
form, the actual variable change s=l-m and its orientation, and the subtraction
identity. Equation (9) is a definition/representation, not a computed value.

For the remaining second term in (7), set t=m-k. On the common sheet

    qm-qi = -t(m+k)/(qm+qi).

At finite delta its height contact at t=0 is exactly zero. Cancelling t against
the PV profile gives the proposed ordinary density

    J(t;l,k) = a mu^2 W L/4 * A(t) A(Q-t)
        * m(m+k) qi/[qm (qo+beta)(qm+beta)(qi+beta)(qm+qi)],
    m=k+t.                                               (10)

Thus the reference mixed ITERATION is C(l,k) H(Q)+integral J dt. Add the saved
closed direct convolution once for the full mixed face-response block. No
integral is evaluated in this proposed instrument. Normal jets are obtained
separately with the native i f qo multiplier after the same grade selection.

The factors, signs, independent grades, two faces and contact are all proposed
source joins, not established by this hand derivation. A disagreement is a
substantive finding; do not modify the native physical inputs or soften a
residual to force agreement.

## 4. Proposed grazing bounds and strength of the limit

Depths and beta share the first quadrant. Consequently
|qi/(qm+qi)|<=1, |qo/(qo+beta)|<=1, and the three shifted denominators in (10)
are bounded below by b0. With |k|<=3, |a mu^2 W L/4|<=2/5,

    |J| <= (2/(5 b0^3)) |A(t) A(Q-t)|
                  (|t|+3)(|t|+6)/|q(k+t)|.              (11)

Reuse the exact simple-root and moving-endpoint argument from the completed
direct-kernel work, with the actual new route joined: 1/|q(k+t)| is bounded by
kappa_min^(-1/2) times the SUM of inverse square-root distances to
-k+/-kappa_delta. It is not a product. Its integral over a measurable set of
measure u is bounded by 4 sqrt(2u)/sqrt(kappa_min). This supplies local uniform
integrability even when a profile contact coincides with an endpoint, or when
both external depths vanish. The factor qi/(qm+qi) must remain together for
this argument; separately bounding it by 1/|qm+qi| would create a false pole.

For |t|>=12, |Q|<=6, use |q(k+t)|>=|t|/sqrt(2),
(|t|+3)(|t|+6)<=15 |t|^2/8 and
|A(t)A(Q-t)|<=150 exp(30 pi) |t|^2 exp(-10 pi |t|).
A deliberately loose explicit tail is therefore

    |J| <= (200/b0^3) exp(30 pi) |t|^3 exp(-10 pi |t|).    (12)

Record the inequalities and the inherited exponential-tail primitive, not a
numerical integral. Together with pointwise convergence off the finite endpoint
set, local uniform integrability and this tail propose an L1 limit for J along
any path in the compact parameter domain. In particular, qi->0 makes J tend to
zero in L1; no hidden contact may then be inferred from individual summands.
The same/opposite external matches l=k and l=-k are included. This does not
prove differentiability in the speed or a bound on a finite inverse.

For (9), A is even, smooth with A(0)=1/(2pi), and exponentially decreasing.
One explicit bound is |A'|<=L/4, hence |j'|<=L^2/8. To assess it, write
A=(1/(2pi)) x/sinh(x), x=pi L s/2. For x>=0,
0<=x cosh(x)-sinh(x)<=sinh(x)^2: the upper gap has derivative
sinh(x)[2 cosh(x)-x]>=0 and zero value at zero. The bound extends by parity.
Then |A(s)[j(Q-s)-j(Q)]/s|<=L^2/(16pi) locally, uniformly in Q. For |s|>=1,
|A(s)/s|<=L exp(-pi L |s|/2) and |j(Q-s)-j(Q)|<=L/(2pi), giving an explicit
integrable tail independent of Q. The contact and (9) are finite and continuous
on Q in [-6,6]. With |C|<=6/(5 b0), (8–9) are continuous at grazing; their
normal-jet multiplier is bounded. At qo=0 this combined contribution vanishes.

## 5. First-order height is a distribution, not an ordinary L1 density

Equations (4) and (6) are bounded diagonal/ordinary kernels on the compact
domain. Equation (5) retains its delta and PV. Do not call the entire retained
face response an L1 function. Its proposed limiting statement is distributional
action on fixed phi in C_c^1((-3,3)) in OUTPUT momentum l, with k in [-3,3].
This states a kernel-distribution result, not a uniform bound on the untruncated
solver, arbitrary input fields, source Fourier transforms or parameter derivatives.

Put B(l,k)=-i mu qi/[q(l)+beta], f_delta(Q)=B(k+Q,k) phi(k+Q), and let chi be
the even indicator of |Q|<=1. The height action is

    (W/4) f_delta(0)
      + W/(2i) integral_R A(Q)[f_delta(Q)-f_delta(0)chi(Q)]/Q dQ.  (13)

The jump at |Q|=1 is harmless. The removed PV integral A(Q)chi(Q)/Q is zero.
For roots in the same first quadrant,

    |q(p)-q(k)|^2 <= |q(p)^2-q(k)^2|=|p-k| |p+k|.

On |Q|<=1 and |k|<=3 the latter gives sqrt(7)|Q|^(1/2). With phi extended
smoothly by zero outside its support,

    |f_delta(Q)-f_delta(0)|
      <= (2/b0)||phi'||_infinity |Q|
       + (2 sqrt(7)/b0^2)||phi||_infinity |Q|^(1/2).       (14)

This makes the local subtracted integrand bounded by a constant plus an
integrable |Q|^(-1/2), uniformly in delta, cs and k. Outside that neighborhood
the denominator is bounded away from zero and the test function has compact
support; |Q|>6 contributes zero. The resulting dominated-convergence argument
proposes a unique limiting distribution, even when transfer zero and a grazing
point coincide. qi->0 makes the height action vanish by the same bound with
its original |qi| factor retained.

The normal-jet height coefficient is i f qo B=f mu qi qo/(qo+beta). Its
supremum is at most 2. Its variation with output depth is bounded by
(4 sqrt(7)/(5 b0^2)) |Q|^(1/2), since |beta|<=2/5. Combining this with the
test-function Lipschitz term gives the same integrable subtraction argument.
Do not assume q(l)phi(l) is C^1 at a root; bound the actual jet coefficient
instead. This is a separate native sign/normal-jet join.

These are proposed analytic arguments to assess, followed by exact algebraic
certificates. They are not a claim that a symbolic program proves measure
theory. Actual finite-matrix actions on their own trial functions will require
their own source/basis/measure and convergence checks.

## 6. Bounded instrument and falsifiable acceptance

After substantive method assessment, prepare one task-local instrument. Restore
the saved two-face closures, original first-shape operands and trace/sign maps,
and completed direct-kernel domain certificates without calling their old
functions. Do not repeat `kernel_bridge`, `reference_pressure_kernels`, a
triangular solve, boundary derivation, profile transform or prior control.
Extract/derive only missing coefficients, source joins and bounds. Unrestricted
frequency/depth maps must be explicit; a symbol with nonzero assumptions cannot
silently represent a depth that reaches zero.

The saved numerical closure matrix is an omega=3 anchor, not evidence of the
complex-frequency continuation. For the new continuation, bind the unrestricted
native definitions Z0=mu/q and I+a Z, verify the proposed generic coefficients
by multiplication in those defining equations (including T Fref=Fphysical),
and then join their omega=3 specialization coefficient by coefficient to the
saved matrices. Do not replace every occurrence of a numeric 3 in an old
expression, assume a frequency-independent ratio, or rerun the old inverse
constructor. These new unrestricted-law identities and the actual saved anchor
are both necessary for the limiting prescription.

Required new evidence:

- Both native first-height/slope coefficients, (1), with the common dispersion,
  both face signs, actual eta/sigma factors, Fourier/delta/measure and units.
- Coefficients of the SAVED closure matrix joined to (2); actual native trace
  and final pressure/normal-jet slot routing joined to (3–8). Preserve each
  uncombined term, not only a simplified result reconstructed from the proposal.
- Full contact and PV operands, both transfer assignments, their variable
  changes and exact residuals for (9–10). Preserve the surviving contact.
- Exact root/denominator/Holder/local/tail certificates and limiting formulas,
  including qi=0, qo=0, both matches and contact-endpoint coincidence. Classify
  ordinary L1 densities, diagonal deltas and height PV distributions separately.
- Responsive controls through the actual expressions: remove the native trace
  subtraction (must change the combined contact and first height), corrupt the
  saved lower lab height while retaining its normal multiplier (must change the
  native trace response), remove the middle resolvent in (2), and use a wrong
  native middle outgoing root in the closed iterated expression. Use declared
  nonzero physical momenta with nongrazing depths for finite response controls;
  also save the failed algebraic identities. Controls demonstrate sensitivity,
  not independent physical correctness. Persist operands/results before guards.
- A grade census separating direct(1,1), iterated(1,1), and trace(1,1), with
  first-shape iteration once and no double integration of the whole direct term.

Unknown or nonzero source/identity joins stop with evidence. No numerical
tolerance, sampling of a limit or new physical input may replace these joins.
No integral evaluation, full source/consumer contraction, finite matrix, producer
regeneration, production overwrite, old benchmark rerun or defect sweep is
authorized by this method document. Guarded execution/readiness and later
result inspection remain separate from method or build review.

## 7. What this would and would not remove from the critical path

If supported and executed, this would establish the retained RESPONSE-kernel
grazing prescriptions alongside the already saved direct block. It would still
leave full physical source/consumer grade and Fourier joins, the corrected
operator's source-address inventory, quadrature implementation and finite-current
controls before a useful defect calculation. The current finite engine uses
momentum cutoffs 4/6, outside the present [-3,3] certificate; no transfer to those
settings is claimed. A later finite plan must explicitly resolve this domain
issue and the existing matrix-route/Gaussian/momentum-zero review debts. Its
untruncated inverse is not bounded by graded-kernel regularity.

No old benchmark, review verdict, accepted uniform result or failed computation
is overwritten. This is still rest bulk, effective speed and LAB_HELD, not the
calibrated/draining model, a leakage bound or a claim that loss varies smoothly
through unity. If source joins or limiting distributions fail, report the actual
obstruction rather than starting a larger exterior-method campaign silently.
