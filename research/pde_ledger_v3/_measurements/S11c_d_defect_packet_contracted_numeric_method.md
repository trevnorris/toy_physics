# Numerical J and direct packet contributions on the certified finite windows

2026-10-04. New method proposal for independent assessment. No new numerical
worker, integral, constant-transform certificate or execution gate exists.
The completed contraction certificate is at 31146c6a. Its original operands,
57 new exact-zero returns and 92 inherited returns are supplied unchanged.
This proposal asks for a concrete numerical subset of the two original packet
actions, with a deliberate change to the constant-field Fourier prescription.
It does not claim an evaluated action or a leakage factor.

## 1. Target, subset and unchanged dependencies

Keep real omega=3, strict rest bulk, LAB_HELD/RHO4_CONSTANT, saved edge momenta
1/5 and 1/10, cs=sqrt(6)/2, kappa=sqrt(595)/10, W=1, L=10, width s=8,
x_u=-5/2, x_v=5/2 and carriers p0=kappa and p0=0. The packets are

    u(x)=exp(-(x-x_u)^2/(2s^2))*exp(i*p0*x),
    v(x)=exp(-(x-x_v)^2/(2s^2))*exp(-i*p0*x).

The pairing is complex bilinear, with no extra conjugation. Amplitudes and
units stay those of the original e_W input and dual THETA_BALANCE test. The
coefficient of epsilon is extracted once. Independent eta/sigma grades remain
independent, with no finite-origin weighting needed by this numerical subset.

Evaluate ONLY the ordinary J part of NATIVE_MIXED_ITERATION and the three
ADDED terms of INHERITED_DIRECT_WHOLE_OFF_DIAGONAL. The former is not the full
native mixed response: its H term remains pending. Flat, height/contact/PV,
slope, H and their numerical controls remain pending under the original
method. The already computed local values are preserved, not recomputed or
added to this incomplete pressure subtotal. Do not call the new J subtotal
the complete mixed contribution or the J+D subtotal the packet action.

The full selected inventory still has 544 addresses. Its saved metadata has
64 J/direct-address entries, of which 20 have formal nonzero-eligible source
and consumer fields. All 20 are pressure slots with source/consumer grades
00 and target/response grade 11. Ten are J entries and ten are direct entries,
with five native jets on each face: e_W, e_W_d1d1, e_W_d2d2, e_W_d3d3 and
e_W_t. The other 44 keep their inherited exact-zero reasons. These are source
metadata observations, not a new mathematical reduction or numerical result.
The future worker must restore and join all 544 entries, their original full
row, all grade maps, 20 complete templates, accepted units and complete old
factor proofs. A different runtime census is fatal, not an opportunity to
change the selected subset. In particular no off-diagonal normal contribution
is inferred from this pressure-only subset.

All prior source, grade, field, response, geometry, unit, rule and numerical
bank operations are inherited without calling their original functions.
New constant-transform identities, request maps and quadrature are new work.
The old finite Fourier/inner banks are not a uniform numerical accuracy proof.

## 2. Explicit proposed exception for exactly constant coefficient fields

Every one of the 20 eligible entries has a saved degree-zero source polynomial
and a saved degree-zero consumer polynomial. Runtime eligibility must join the
full saved polynomial and denominator, reconstruction input and zero return,
source/consumer original expression, constantValue, field ID, actual native
jet, accepted unit input/return, normal factor and complete response map. A
constant flag, field hash or value at x=0 alone is insufficient. A zero source
or consumer is retained as a proved zero, not a manufactured transform request.
No unknown/nonconstant field is routed through this exception.

Let b and c be the two actual constant coefficients, n=n1 the profile-axis
derivative order, and P_j=(-3i)^nt*(i/5)^n2*(i/10)^n3. The proposed exact formulas
in the original Fourier convention are

    G_u(k)=s*sqrt(2*pi)*exp(-s^2*(k-p0)^2/2-i*(k-p0)*x_u),
    G_v(l)=s*sqrt(2*pi)*exp(-s^2*(l-p0)^2/2+i*(l-p0)*x_v),
    X(k)=b*P_j*(i*k)^n*G_u(k)/(2*pi),
    Y(l)=c*G_v(l).

Derivatives act on u before multiplication by b. Integration by parts is
allowed here because b is exactly constant and the Gaussian and its derivatives
decay at both ends. The carrier derivative is i*k after Fourier transformation,
not i*p0. The Y sign follows its actual argument -l and test carrier -p0;
no conjugation of c, P_j or the response is introduced. The physical derivative
units belong to k and the native jet; numerical magnitudes do not erase them.

This is a substantive change from requiring new GL24/48/B50 Fourier quadrature
at EVERY request. For this exactly joined constant subset, replace those
Fourier integrals by the analytic formulas, including their all-real Gaussian
tails exactly. There is no finite x-radius or Fourier truncation error for the
analytic identity. Its floating evaluation still has roundoff. No such change
is proposed for any nonconstant field or for H/profile integrals. Do not invoke
the original Fourier request/rule constructor merely to recreate old tests.

Two numerical implementations remain separate. Route A uses the displayed
Fourier derivative identity at 30 digits. Route B at 50 digits independently
forms the physical Gaussian derivative polynomial Q_n(w), starting with 1 and
using Q_next=Q_n'+(i*p0-w/s^2)Q_n, and evaluates its Gaussian moments:

    M0(nu)=s*sqrt(2*pi)*exp(-s^2*nu^2/2),
    M1(nu)=-i*s^2*nu*M0(nu),
    M(r+1)=-i*s^2*nu*Mr+r*s^2*M(r-1).

Use nu=k-p0 and multiply exp(-i*nu*x_u) and the source normalization for X.
For Y, use the test's actual Fourier argument -l, carrier -p0 and center x_v;
its nu is p0-l, and there is no 1/(2*pi) factor after the declared 2*pi in Y.
The future build derives and saves the moment/derivative algebra for every
actual n and checks exact equality to the displayed form before quadrature.
The Gaussian transform theorem is assessed mathematics, not machine analysis.

At each new numerical request compare the two formula evaluations at that
request's EXACT same represented momentum (conversion preserves its MP tuple,
not a newly rounded intended node), using a private checking context. Record
the actual value/difference and require absolute difference <1e-12 in the
fixed transform unit. The independent-route values are not replaced by the
check's values or shared between routes. These are algebra/roundoff checks,
not independent numerical Fourier integrals. The two complete action routes
therefore share an assessed analytic Gaussian identity; state that limitation.
Persist both formula operands before refusing any mismatch. A future worker
must also compare original saved bank constants only when full arguments and
settings match, without rerunning them; an unmatched bank point is not a miss.

## 3. Addressed contraction and permitted sharing

Restore the certified four finite-window formulas and their full old/new
operands. Put E(p)=q(p)+beta, q(p)=sqrt(kappa^2-p^2) on the positive-real /
positive-imaginary outgoing branch, a=100/109+30i/109, mu=3/10,
beta=30/109+9i/109 in their established physical units. Keep
Cj=a*mu^2*W*L/4 and Cd=(-i*mu)*W*L/(4i), with original phase provenance.
A(z)=L*z/[4*sinh(pi*L*z/2)], A(0)=1/(2*pi), is unchanged.

For each address retain x=X/E and y=N(l)*Y/E, with the actual native N(l).
Runtime N=1 is required for the selected nonzero subset; do not silently remove
normal factors from the rest of the inventory. The certified definitions are

    Y0=integral_[-K,K] A(l-m)*y(l) dl,
    X1=integral_[-K,K] A(m-k)*k*x(k) dk,
    XrT=integral_I(m) A(m-k)*k^r*x(k) dk, r=1,2,
    X02T=integral_I(m) A(m-k)*q(k)^2*x(k) dk,
    CrT=integral_I(m) A(m-k)*k^r*q(k)*x(k)/(q(m)+q(k)) dk, r=0,1,2,
    YrCT=integral_I(m) A(l-m)*l^r*y(l)/(q(m)+q(l)) dl, r=0,1,

where I(m)=[max(-K,m-T),min(K,m+T)], M=K+T, m in [-M,M]. Then

    BJ=Cj*integral Y0*m/[q(m)*E(m)]*(C1T+m*C0T) dm,
    BDh=Cd*integral Y0/q(m)*(C2T+m*C1T) dm,
    BDq=Cd*integral Y0/q(m)*X02T dm,
    BDr=Cd*integral X1*(Y1CT+m*Y0CT) dm.

Both original windows (K,T)=(27,122) and (29,124) are evaluated. The full outer
wings remain: M=149 or 153, with clipping transitions at |m|=95. J/Dh/Dq clip
INPUT k; Dr clips OUTPUT l. The actual reflected root is q(l-t); it is not
converted to q(k+t). All four primitives are kept separately before addition.
Neither finite quadrature sums nor individual numerical errors are inferred
from the exact continuous identity.

Constant coefficients can be factored out ONLY after new exact full-expression
joins. A candidate normalized family uses x_n(k)=k^n*G_u(k)/(2*pi*E(k)) and
y0(l)=G_v(l)/E(l), with multiplier alpha=b*c*P_j*i^n. Actual n is 0 or 2.
This is exact linearity, not low-rank approximation, field pruning or averaging.
Sharing across faces additionally requires the actual complete J or D template
identity and units, not equal face labels. Keep every address and alpha; never
cancel them to pay a numerical error budget. If a family equality is not exact,
the build must keep separate families. Sharing evaluated values is allowed only
within the same route, precision, carrier, window, complete mathematical and
numerical request identity. A24, A48 and B50 do not share evaluated transforms,
contractions, adaptive choices or action values. Their immutable exact source
and previously constructed rule objects may be shared.

## 4. Complete panel geometry, including inner grazing layers

Build a new exact plan from the saved dependencies; the old square mesh is
ancestry, not a mesh to replay or an excuse to integrate absent root products.
Transport every original collision label per primitive as certified. Labels
absent from a primitive remain recorded with their actual dependency reason.
All new interval endpoints and order changes must be explicit before nodes.

For an inner variable z=k or l, start with the appropriate entire [-K,K] or
I(m). Include fixed cuts z=+-kappa and z=p0+d for d=0,+-1/s,+-2/s,+-4/s,+-8/s.
Include z=m+d for d=0,+-1/L,+-2/L,+-4/L,+-8/L. To resolve q(m)+q(z) near
simultaneous grazing, also include z=+-kappa+-d_g(m), where
d_g(m)=min(|m-kappa|,|m+kappa|). These are resolution cuts, not extra physical
singularities. Their purpose is to resolve the scale |q(z)| approximately
|q(m)| on each side of either simple root, including m approaching a root.
Retain both branches and all coalesced labels. There is no pointwise 0/0 value.

The m plan includes endpoints +-M, clipping transitions +-(T-K), +-kappa, 0,
and the same carrier and profile offset cuts around p0. On each of the four
sectors divided by -kappa,0,kappa, all the above inner-boundary functions are
affine in m: fixed cuts, m+d, the appropriate clipped window endpoints, and
the two signs of d_g. Include EVERY intersection of these boundary functions
inside the outer sector, including crossings with +-K and the active clipped
endpoints. Coalesce only exactly equal points, retaining their labels. Then
verify order constancy, disjoint open intervals and complete inner coverage on
each m slab. There is no omitted narrow strip or clipped wing. This is a new
plan certificate in the future numerical worker, not a completed geometry job.

On each original interval, Route A splits at the midpoint and uses open
z^2 endpoint maps on each half with positive Jacobians, followed by the saved
GL24 or GL48 rule. Apply this to m and the actual m-dependent inner intervals.
Both A levels use the same exact interval plan, but their actual nodes remain
different. Do not substitute a uniform tiny spacing across the entire wings;
the stated exact resolution/branch/window plan defines all initial intervals.
No posterior interval removal or extra A refinement campaign is authorized.

Route B uses the same mathematical boundaries and coverage obligation with
independent physical-variable G7/K15 adaptive panels at 50 digits. No squared
Route A nodes, values or subdivisions are used. Embedded errors are empirical.
Select the active leaf with largest normalized error, split at its physical
midpoint, replace that leaf in the global error sum, and preserve both parent
and child records. This global criterion avoids the old pathological local
tolerance-halving policy. Rounding alias, precision stagnation, nonfinite or
unsupported endpoint behavior is a failure; do not change precision or rules
after a miss. Root endpoints are never nodes. The profile's removable value
and existing finite sinhc expansion keep their original bound and source join.

## 5. Nested error accounting and comparison rules

No inner error is treated as zero because the outer rules agree. Store every
scalar or vector evaluation as value plus separately named empirical indicators
and analytic tails. For a product use the nonnegative propagation expression

    e_product=|a|*e_b+|b|*e_a+e_a*e_b,

and for a linear sum use the sum of absolute coefficient times each indicator.
Apply it before cancellations, including Cj/Cd, external resolvents, m powers,
alpha, rule weights and Jacobians. Source Gaussian formula discrepancies and
profile-approximation bounds enter the same path. This propagation of empirical
indicators is not a rigorous forward-error enclosure. Do not label it one.

For each carrier/window use a numerical budget epsilon=1e-11/80 per original
eligible address and primitive (20*4 is conservative). Zero/nonapplicable
primitives do not donate their budgets. A shared normalized family must satisfy
the tightest scaled budget of ALL its original addresses; multiplying a family
back by a small alpha may not hide a failed other address. Controls have their
own identical allocation. This is internal resolution bookkeeping; the final
physical tolerance is unchanged.

For B, solve the inner vector integrals jointly at a given m until the absolute
propagated density indicator for every affected addressed primitive is at most

    epsilon*w(m)/(16*W_M),  w(m)=1+1/|q(m)|, W_M=3*M+4.

Here m and q denote magnitudes in the fixed original momentum unit, so the
addition in w is dimensionless. This must be explicitly nondimensionalized in
the build. Since kappa>2 and M>kappa,
integral_[-M,M] 1/|q(m)| dm=pi+2*acosh(M/kappa) <4+M,
giving integral w < W_M. This is a new analytic allocation lemma to assess and
join, not an evaluated response integral. It avoids demanding vanishing inner
error merely because an integrable outer 1/q factor grows near a root.

Start every inner active leaf from its prescribed G7/K15 pair. Propagate all
active-leaf indicators through the complete addressed products. Refine the
leaf with the largest contribution to the worst normalized density indicator
until all density criteria hold. Decisions use no outer or other-route result.
This adaptive loop is the predeclared integrator, not a retry of a failed run.
At each outer B leaf, also include the inner indicators multiplied by the
absolute Kronrod AND embedded Gauss weights when assessing K-G. Stop the
outer global refinement only when all addressed sums satisfy their epsilon
budget, and persist separate outer K-G and propagated-inner totals. Actual
weighted inner totals must be <=epsilon/4; the analytic allocation lemma alone
does not certify a discrete sum of w.

For each A level, at its OWN m nodes evaluate the inner GL24/48 pair with the
same m and mathematical intervals. Its base/refined difference is the local
inner indicator; the A24 outer value uses its own inner24 values, and A48 its
own inner48 values. The extra inner pair is for error accounting, not sharing
an evaluated A24 request with A48. Sum the propagated indicators using the
actual positive outer weights and require <=epsilon/4 for each addressed
primitive. If this fails, preserve it and stop without adding orders or panels.
Independently sum the finished A24/A48/B50 values and keep the full differences.

For EACH primitive, face and independent grade, and for the subsequent J/direct
subtotals, retain the original comparisons:

    tolerance = 1e-9 + 1e-7*abs(A48 reference value).

Require A24 versus A48, A48 versus B50, and corresponding base/enlarged-window
changes to pass. An unresolved Dr, Dh or Dq is not hidden by cancellation in D.
Keep actual per-address indicators even where final required comparisons are
per-face/component/grade. Both carriers must be assessed separately.

Restore and join the original positive per-address all-real tail budgets for
these same J/D integrands. Reordering the IDENTICAL finite windows creates no
new truncation. Do not use the compact Fubini C/|q(m)| bound as an all-real tail
estimate. No cancellation pays a tail budget. Save base/enlarged results and
tail allocations separately. Conditional analytic tails plus propagated
empirical indicators and route/refinement/window differences form a declared
empirical envelope, not a certified total error bar, flux uncertainty or loss.

## 6. Numerical controls, identities and preservation

Use the certified wrong-reflected-root formula on the actual address with the
smallest eligible direct ID, separately for each carrier. Replacing q(l-t) by
q(k+t) requires m=k+t, a clipped INPUT, an unclipped OUTPUT, and

    Dr_wrong=Cd*integral {2*X1T*Y1C+(X2T-m*X1T)*Y0C} dm,

where YrC is unclipped. Do not reuse correct Dr's clipping. Persist the complete
baseline/mutant values and their difference with both routes/window checks.
Require finite measured movement >10 times their summed empirical envelopes;
silence is unestablished coverage, not permission to choose another address.
This control supplies no normal-slot/H-contact/Leibniz coverage: those remain
required when the other pressure contributions are evaluated.

As a Fourier-adapter control use the smallest eligible n1=2 address, replacing
(ik)^2 by (ip0)^2 before the SAME complete J assembly. This is a new derivative
routing mutant, not the nonconstant-coefficient Leibniz control. The carrier-zero
mutant can be exactly zero while the baseline need not be. Use the same finite
movement/envelope criterion per carrier; no asserted outcome in advance.

Each numerical request identity includes the full exact field/jet/unit/template
and source-proof identities, native coefficients, carrier/centers/width,
component and mutant status, window/clipping variable, exact parent interval,
node descriptor AND represented MP tuple, route, precision, saved rule receipt,
panel/adaptive settings, profile branch algorithm and applicable error target.
Hash equality alone is not equality of full operands. If a completed request
is reused, persist its full input join and original return receipt. A looser
completed error target is not silently accepted as a tighter one. Old fixed
bank values are not requested again; unmatched settings simply prohibit reuse.

Use the proposed lossless SQLite FULL-transaction journal only after concrete
build assessment. Its existing stdlib tests are tooling evidence, not scientific
integration clearance. Keep route namespaces separate. Store all numerical
nodes, weights, function operands, returned vectors, active/refined leaf ancestry,
inner-to-outer error propagation, sums and controls as exact encoded MP tuples
with byte/hash receipts. Commit input before evaluation and the returned panel
before any scalar guard. Persist a complete prefix on failure. No lossy decimal
serialization, output pruning, old-bank rewriting or destructive migration.

Process at most one inner panel record (<=48 nodes) at a time; release its live
arrays after durable storage. Keep adaptive leaves and immutable request indices
in SQLite rather than retaining every evaluated node in RAM. Bound an optional
route-local LRU to 64 MiB of explicitly accounted stored payload; evicting a
payload must allow exact immutable retrieval, not scientific recomputation.
Keep SQLite's small fixed cache and file-backed temporary tables. A build must
test this orchestration with synthetic records before source assessment.
Check available disk before launch and before bounded record batches, reserving
20 GiB; stop preserving complete records on resource refusal. This is a resource
safeguard, not an elapsed-time or operation-count deadline. Total adaptive work,
storage and runtime remain unknown; no speedup or completion-time forecast.

## 7. Build obligations and result boundary

After method assessment, prepare ONE concrete guarded numerical build containing
the new constant-field/moment/source joins, exact panel plan, nested accounting,
two independent implementations and lossless record orchestration. These joins
can run as the preamble to the bounded numerical instrument after build review;
another standalone certificate campaign is not proposed. Pin all actual source,
saved inputs, old rule receipts, helpers, method/build records and authority.
No method verdict alone launches an unspecified worker. Preserve every old
review/failure and do not relabel inherited proofs as fresh computations.

The run remains one 4 GiB native/cgroup worker under the unchanged shared guard
and existing supervisor, within the 16 GiB pool, zero swap, one assigned CPU
and thread, 32 tasks, 4 GiB host reserve and desktop priority. No wall/native/
CPU/inactivity deadline, scheduler change or unguarded fallback. Arm the existing
local completion hook first. No automatic scientific retry or model polling.

A passing result would be two carriers' addressed J and three-direct primitive
tables on both finite windows with conditional all-real tail/resolution evidence.
It would leave the other pressure components and complete pairings pending.
No full pressure action, incoming mode, inverse, current normalization,
scattering, leakage factor, loss, drain/calibration or defect/speed sweep follows.
