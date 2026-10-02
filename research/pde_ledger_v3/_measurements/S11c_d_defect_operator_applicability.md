# Defect operator applicability and the next bounded correction

Revised 2026-10-01. This is a proposed method and source-based decision, not a
new calculation or independent clearance. The first reviewed draft is preserved.
The explicit rules below are for assessment before bounded implementation.
The objective remains a useful rest-bulk near-unity defect test, with the drain
and primitive calibration limitations stated. Repeating completed uniform or
selected mixed-term calculations is unnecessary.

## What is established, and what remains different

The selected uniform near-unity schedule completed on both sides of each end's
actual modal/acoustic match. Its selected transverse states and currents were
finite, with zero selected face drives and joined grazing limits. That result
concerns constant ends, not propagation through the profile.

The upper-face tanh diagnostic supplied a nonzero direct mixed integrand, with
a sign/endpoint argument, at profile momenta 0 to 1/10, omega=3, c_s0=10.
No numerical integral was evaluated. The saved boundary coefficient
before that specialization keeps the input momentum and both transfers live.
It passed the original boundary and native linear-coefficient identities, but
has not been completed as a physical operator at arbitrary external momenta.
The selected addition survives the closed-face and reference-pressure maps.
It is distinct from the existing iterated first-shape contribution.

The source/consumer diagnostic now supplies complete selected controls. Its
zero-grade source annihilates arbitrary transverse curl on both faces; its
scalar and longitudinal sources and upper-face scalar consumers are nonzero.
The retained contribution is kernel(1,1) times source(0,0) times consumer(0,0).
It does not directly supply transverse forcing at that tested grade. The native
U pressure-slot absence is a complete selected-row census, not a claim about
every other coupling. The full source expressions and their higher grades remain
saved. None of these results constructs the missing lower-face operator or a
complete nonuniform response.

The packet supplies actual operands, original methods and addressed controls.
They are author-run evidence, with inherited scope and provenance limitations.
No earlier independent report or outside commentary is supplied to this review.

## Why the numerical benchmark is not automatically protected

Native c2 `retained_shape` keeps grades (0,0), (1,0), (0,1), (1,1) after closing
the slab rows. The numerical finite library then binds contrast, assembles a
finite matrix and applies a full LU solve, cross-checked with SVD. It does not
truncate the inverse solution to that rectangle. Thus a vanishing direct
transverse-source column does not by itself make a missing scalar-sector
matrix entry irrelevant to its finite-contrast response.

There is a useful *conditional* formal statement, proposed for assessment only.
For a regular graded boundary-value problem, changing only L11 gives
L00 delta_psi11 = -delta_L11 psi00, with any boundary/incoming-data changes also
included. If psi00 is transverse, delta_L11 annihilates it, the boundary data
are unchanged, and the reference problem has a uniquely specified inverse,
then delta_psi11 is zero. Missing joins include the actual reference solution,
both-face/source/consumer applicability, regularity and the exterior/incoming
boundary prescription. This is not yet a source-joined theorem for S11c-d.

Even if established, that statement does not bound the untruncated finite solve.
At finite contrast, an induced scalar field can feed the missing scalar block;
an exterior pole or threshold can defeat an informal small-order argument.
No physical-current or power equivalence follows merely from a field-grade
identity. Do not spend a separate scientific job proving only this formal
statement unless reviewers identify a concrete way it removes required work
for the requested finite-defect result. Do not relabel the old 0.12 benchmark
as a corrected-operator or calibrated-model answer.

## Recommended next construction: one-dimensional direct increment

Prepare a task-local, source-bound correction for the actual one-dimensional
tanh profile, LAB_HELD/RHO4_CONSTANT, fixed conserved edge momenta (1/5,1/10),
with both faces. Keep the profile-direction input/output momenta live. Carry the
effective bulk speed as a declared parameter for later use near the two actual
modal matches; this is not a primitive calibration and the bulk remains at rest.
Do not expand to arbitrary profiles, anchoring cases, density rules or a full
three-dimensional shape campaign. Do not overwrite or regenerate the production
exports merely to construct this increment.

### 1. Fixed arguments, physical sheet and both faces

Let k be the profile input momentum, Q=k_out-k, H=t, S=Q-t, and
b^2=(1/5)^2+(1/10)^2. Keep omega=3 and c_s>0 explicit; c_s is the effective
parameter, not a derived primitive calibration. Put Kappa^2=omega^2/c_s^2-b^2.
For real p and omega>0 use

    q(p) = sqrt(Kappa^2-p^2)       if Kappa^2-p^2 > 0,
           i*sqrt(p^2-Kappa^2)    if Kappa^2-p^2 < 0.

This is the positive outgoing/decaying root. At a zero radicand the limiting
sheet is omega -> omega+i*0^+ with c_s fixed, not a pointwise division by zero.
Do not use solver ordering to select a root or use the opposite root below.
The four distinct depths are

    qi=q(k), qh=q(k+t), qs=q(k+Q-t), qo=q(k+Q).

In native middle coordinates m=k+t, the reflected route has momentum
k+k_out-m. Give its depth a separate symbol with the same rule; it is not
MIDDLE_Q. No independent depth assignments may stand in for this shared sheet.

Reuse the saved generic upper coefficient, its boundary equations, modes and
residuals without replaying that derivation. Reconstruct only missing lower-face
or general-argument work. For face f in {+1,-1}, use the actual graph reference
f*W_0/2, signed lab displacement f*h, lab slope f*s, and native unit normal
(-f*grad(h_lab), f)/sqrt(1+|grad(h_lab)|^2). Here h and s are outward height
and slope amplitudes whose native factors are joined below. Keep these lab and
outward conventions distinct: merely changing f at a fixed lab slope is not the
same as the actual mirrored thickness profile. The outgoing extension is
exp(i*f*q*(w-f*W_0/2)), pressure is the native i*omega*rho_m times potential,
and prescribed normal velocity uses this native normal. Do not also send q to
-q. Use the actual minus pressure/normal-jet slots for f=-1.

Extract only the independent height-slope grade and join zero/first orders to
native c1 for each face. Save the original boundary residuals, normal/location
and pressure convention, and both sign maps before checking. The normal's
quadratic slope normalization cannot silently supply or remove mixed terms.
An outward-coordinate mirror comparison is an extra check, not the derivation.
All new claims retain independent eta and sigma_W until grade selection.

### 2. One ordered convolution and its contact/PV pieces

The saved upper coefficient, per unit height and slope and prescribed normal
velocity, is

    C+(k;H,S) = i*omega*rho_m * N/(qh*qi*qo),
    N = -H*qi^2 + k*qh*qi + k*qh*qo - k*qh*qs - k*qi^2,

where qh=q(k+H), qs=q(k+S), qo=q(k+H+S). It is an ordered coefficient.
The native shape factors per grade are

    h_hat(t) = (W_0/2) * w_hat(t),
    j_hat(u) = (1/2) * jet_hat(u),
    w_hat(t) = delta(t)/2 + PV[A(t)/(i*t)],
    A(t) = L*t/[4*sinh(pi*L*t/2)], A(0)=1/(2*pi),
    jet_hat(u) = L*A(u), L=L_W.

Join these factors to actual native shape_source/profile_bindings and preserve
both edge deltas until the reduced measure is established. For each derived
face coefficient Cf the ordered object is exactly

    Df(k_out,k) = integral_R Cf(k;t,Q-t) h_hat(t) j_hat(Q-t) dt.

It is a whole-convolution *definition*, not a request to evaluate the integral.
Do not add a second full-weight swapped coefficient. The reference computation
uses this ordered form. If a comparison is symmetrized, only
[F(t)+F(Q-t)]/2, with F the entire ordered integrand including its hats and depth
arguments, is a change of integration variable; averaging Cf alone is not.
Persist an addressed times-two control. No symmetry of Cf is assumed.

Before any cancellation, retain the explicit candidate contact and PV parts:

    contact_f = (W_0/4) Cf(k;0,Q) j_hat(Q),
    pv_f = (W_0/2) PV integral_R [Cf(k;t,Q-t) A(t) j_hat(Q-t)/(i*t)] dt.

The contact uses the complete coefficient, not separately H-independent-looking
summands. For the upper coefficient and regular external depths, the common
sheet gives qh=qi, qs=qo at H=0, hence N=0. More explicitly,

    qh^2-qi^2 = -H*(H+2*k),
    qo^2-qs^2 = -H*(2*k+2*Q-H),

and, where the denominators are nonzero, the proposed exact identity is

    C+ = H*B+,
    B+ = -i*omega*rho_m * [
        qi/(qh*qo)
        + k*(H+2*k)/(qh*qo*(qh+qi))
        + k*(2*k+2*Q-H)/(qi*qo*(qs+qo)) ].

This is an author hand derivation to be checked against the saved expression,
not a new executed symbolic result. Persist its original operands, both
physical dispersion identities and exact residual. For regular external qi,qo,
B+ is regular at t=0, so the upper contact is zero and t*PV(1/t)=1 gives the
ordinary raw integrand (W_0/2)*B+*A(t)*j_hat(Q-t)/i. Preserve the uncancelled
contact/PV definition alongside this result and its domain. Derive and check the
corresponding rule for the actual lower coefficient; do not assume it.

This regular-domain cancellation is not a distribution product at a coincident
grazing singularity. If qi or qo vanishes, or a factor sum vanishes, retain the
limiting expression and domain refusal until section 5's closed limit is
established. No numerical smallness or unknown-to-zero promotion is allowed.

### 3. Native measure and direct insertion

Native profile_bindings defines each three-dimensional profile hat with
exp(-i delta_k*y)/(2*pi)^3; its two constant-edge transforms therefore give
one delta in each edge direction. Source fields use the unnormalized forward
transform and the inverse in kernel_apply carries (2*pi)^-3. These are distinct
objects. Join the executable definitions, not an ambiguous comment describing
only one of them. Their composition is

    f(x) = integral f_hat(k) exp(i*k*x) dk

for the normalized profile hat, so profile multiplication contributes a plain
middle integral with no extra 2*pi. In the native second slot the two profile
hats give two sets of edge deltas; integrating over the two middle edge momenta
leaves exactly the output-input edge deltas. The reduced transfer Jacobian
m_profile=k+t is +1. Keep source transform/inverse factors as in kernel_apply.
Persist this reduction and compare its normalization and dimensions to the
native iterated first-shape second slot. A disagreement is a method finding,
not permission to insert a compensating constant.

Represent the new direct slot by the raw product in section 2 using native
profile hats before reduction, or its explicitly joined reduced counterpart.
Keep any contact separate until its reduction is proven. Add it once at
z_three[0,2] for each face with eta*sigma_W; the native first-shape iteration
remains once. Never insert the inherited whole D into kernel_apply(second=...),
which already integrates the middle leg. Compare to the saved k_in=0 integrand
and whole-D closure at their actual saved point without recalculating them.
No new integral, finite matrix or baseline regeneration is part of this step.

### 4. Closure, consumers and reuse

Use the actual source/density/memory/epsilon, resolvent and reference maps on
both faces. At zero grade write z0(p)=rho_m*omega/q(p) and let a_f be the actual
coefficient of Z in the inherited inverse I+a_f*Z. The direct closed difference
must join the two external factors

    R_f(k_out) Cf R_f(k),  R_f(p)=q(p)/(q(p)+a_f*rho_m*omega).

Keep leg-dependent factors explicit if the actual source gives them. Verify
this from the native ordered closure, not by assuming that a numerical factor
saved at one point is a general function. In particular, before pointwise
external grazing use, require the rational identity that removes explicit
1/qi and 1/qo by these R factors, with nonzero resolvent denominators stated.
This algebra alone is not an integrable uniform grazing limit. The reference
trace inverse and i*f*qo normal derivative must be carried too.

Contract the increment into the full native slab rows before extracting grades
(0,0),(1,0),(0,1),(1,1). Preserve pressure and normal jet separately, even where
a selected zero-grade jet consumer vanishes. Carry actual units, epsilon and
source/row addresses. The target is the **(1,1)-complete rectangle model**;
pure eta^2 and sigma_W^2 are not included. It is neither a full second-order
shape answer nor a physical loss estimate. No relative magnitude of omitted
pure grades has been established by the selected diagnostic.

List exactly which coefficient records change and establish actual both-face
constant-end/zero-jet specialization. Only these joins allow unchanged end
modes, currents and trace maps to be reused; a changed file hash alone does not
require replay. Keep local/first-shape entries unchanged after actual operand
joins. New scalar entries mean old finite matrices remain uncorrected. Later
source-address and term inventories must explicitly include new entries; their
375/160 guards cannot be bypassed. The present output is a separate increment
and its provenance, not an overwrite of accepted production exports.

### 5. Branch/contact coincidence and bounded acceptance

When Kappa^2>0 the potential middle branch locations are

    t = -k +/- Kappa, and t = k+Q +/- Kappa.

The latter route may be a numerator cusp in the unfactored expression; determine
orders from the complete physical coefficient rather than from an artificially
split summand. Retain contact t=0, removable profile factors and coincident
locations. At Kappa^2<=0 record the corresponding domain and any zero endpoint
instead of using a real square root without qualification.

The saved c_s=10, k=0 sign/tail calculation does not transfer to near-match
speeds. As a nonzero k approaches input grazing, the approaching root has
|t_b|=|Kappa-k|=|qi^2/(Kappa+k)|, of order |qi|^2/(2*|k|) for the positive
branch. The exact root location, not that asymptotic scale, sets any future
endpoint partition. Negative momentum uses the corresponding root.

Before accepting exact-match applicability, save a symbolic limit of the
*closed* combination including the reference maps on the same outgoing sheet,
from both sides. Resolve the shrinking contact/branch layer, for example using
t=t_b*u when t_b!=0, with orientation/Jacobian and the other branch locations
carried. A pointwise finite limit alone does not justify exchange with the
middle integral: supply an integrable local bound or a separate limiting
integration prescription, plus tail behavior and pole exclusions. If this
requires a larger method, record the exact obstruction and stop; do not silently
invent a regulator or a new exterior method. Exact-match finite matrices stay
blocked if this step is unresolved. A nongrazing raw increment may still be
reported with its restricted domain, never as uniform near-match applicability.

Necessary addressed checks are the native boundary/linear joins, physical
four-depth dispersion, contact/factorization residuals, normalization/units,
closed external factors, lower-face orientation and both-face source/consumer
joins. Direct-slot omission, native slope omission, wrong sheet, times-two
insertion and wrong face sign must respond where applicable. Save inputs and
results before checking; exact rational normalization may decide a residual but
must preserve the raw form. A failed control or unjoined source is unresolved,
not a zero increment. No producer, completed calculation or original integral
is replayed to make these checks.

The construction stops at this increment and its applicability evidence. It
contains no evaluated middle integral, finite solve, defect sweep or loss
number. Before finite matrices change, every new contact/PV/branch piece needs
an explicit supported integration recipe. A nonintegrable or unjoined piece
is a substantive method finding, not a reason to loosen a numerical threshold.

## Review decision requested

Assess the foregoing route and the already saved selected evidence, with the
following concrete outcomes distinguished:

- Which selected source/consumer conclusions are supported by the actual
  operands, and which gaps still block that narrow result?
- Does a justified retained-solution exemption eliminate any required work for
  the *untruncated* numerical defect result? If not, state that directly rather
  than requesting an otherwise unhelpful extra formal diagnostic.
- Is the proposed one-dimensional, both-face raw-increment route sufficient in
  scope and representation? Identify the smallest missing physical operand or
  distributional rule before implementation; do not demand an unrelated full
  model certification.
- Which saved uniform objects and numerical machinery can be reused, under
  which actual joins, and which matrices/solves must change?

This is a source/evidence method assessment, not a claim that a future worker,
lower-face calculation, new integral or physical response has passed. A clear
method verdict would permit faithful implementation and its necessary guarded
checks under the user's continuing direction. Actual readiness still requires
pinned code/inputs, checkpoint preservation and enforced resources. A substantive
method change would require its applicable assessment, not author clearance.

No new science has run during this revision. No future worker has been cleared. Scientific work will use the
existing no-deadline pooled guard/supervisor, zero swap, native memory limits,
disjoint CPU assignments and aggregate reservation, with the local completion
hook armed first. No replay, time cutoff, scheduler change, automatic retry or
new external export is authorized by a process exit. Costs of general contact/
endpoint handling remain uncertain; the five-second selected continuation is
not a runtime estimate for this construction.
