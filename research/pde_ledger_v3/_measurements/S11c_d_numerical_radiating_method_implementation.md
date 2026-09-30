# Numerical radiating pilot — user-authorized implementation method

2026-09-30. The reviewed round-2 draft is preserved at `ba15c87b` and in
`S11c_d_numerical_radiating_method_r2.md`. Claude literally cleared that method
with the five small corrections below. Grok delivered no final report. The
user then expressly directed: “No other grok pass make the minor corrections
and run it.” This is user authority to proceed without the missing report,
not independent paired clearance. No further review is queued.

The user also withdrew the automatic two-day/16-hour authoring stop and will
decide the time budget. This does not alter the scientific method, the finite
case list, physical inputs, failed-method stop, or ordinary per-job containment.
Exact RIGHT stays parked. Actual results remain uncomputed.

## Deliverable and finite scope

One case, LAB_HELD/RHO4_CONSTANT, tangential momenta `(1/5,1/10)`, frequencies
**3, then 2.3, then 4**. The near-threshold point is above the saved acoustic
threshold `sqrt(5)`. The analog-light frequency calibration remains open.
Material constants, physical face permeability and memory, and the saved
profiles remain those in `S11c_d_variable_profile_development_input.json`.

Let `a=0,1/4,1/2,1` multiply the complete saved contrast: `eta_bg=a/100` and
`sigma_W=eta_bg*W_0/L_W`. Both original thickness and density-profile
perturbations follow this binding. Independent eta/sigma source grades remain
distinct until binding; this is not thickness-only forcing. At `a=0` the same
assembly is bound to its reference background, without deleting rows by hand.

For each of four transverse incident columns report incident, reflected and
transmitted transverse current and

    D_T = 1 - (J_T,reflected + J_T,transmitted)/J_T,incident.

Reflection counts as surviving transverse output. Contract the actual
polarized quadratic current matrices, including cross terms within each
output block, instead of summing squared field amplitudes. This is a finite
numerical solution of the saved retained operator, not a Born expansion.

Two statuses travel separately with every number:

- Numerical status: the signed finite-model deficit and whether the declared
  numerical controls resolve it.
- Physical interpretation: **RETAINED_ORDER_LOSS_INTERPRETATION_UNRESOLVED**.
  This pilot does not establish that the deficit equals physical leakage;
  convergence and contrast-squared scaling alone cannot do that. It does not
  separate absorption, radiation and conversion or certify a loss bound.

This distinction is part of the requested energy-balance report, not a reason
to suppress a successfully computed finite-model number. A boundary, pairing
or integration failure instead makes that numerical number unresolved too.
If the central case fails a required method/control check, preserve its data
and stop before the other frequencies. No automatic retry or larger exterior
method follows a failed point.

## Sources and numerical binding

Restore the accepted frequency-source packet's original `liveFrequencyAndGrades`
records, addresses, ordered integrals and unit maps. Its 375 records cover
the local matrices, cells, 80 nonlocal rows and 35 source amplitudes. Bind
frequency, material and tangential values before new numerical compilation.
Do not reconstruct a producer, general frequency derivatives, epsilon Poly
carrier extraction, exact root census or the earlier failed premise work.

Reuse finite Chebyshev/source/profile integration and linear-algebra structure.
Reassemble each changed numerical operator. Reuse completed arrays only when
the actual inputs match. Changed-contrast end tables are bound from the saved
original algebraic pencils and mappings; later frequency-only tables already
fixed at nonzero contrast are not substituted for them.

Keep physical depth momentum distinct from the scaled algebraic radical.
The old `abs(omega-1)<=1/4` chart, its negative-radicand proof and its
`i*sqrt(-radicand)` replacement are not used at these real radiating targets.

## Real bulk branch and endpoint integration

c1 supplies flat impedance `Z0=rho_m*omega/q`. The closure uses the **memory
kernel** `lambda_A(omega)=Lambda_A_0/(1-i*omega*tau_A)` and

    R(q) = 1/(1 + lambda_A(omega)*Z0/rho_m**2).

It does not replace that kernel by the static parameter `Lambda_A_0`. Physical
permeability remains nonzero. c2's ordered second scattering retains its output,
middle and input legs. No sponge, changed material loss or imaginary-frequency
regulator is added to the real-frequency field operator.

After numerical binding, inventory every radical-bearing momentum leg, its
radicand, physical scale and real zeros. The source at the saved tangents
suggests `b**2=omega**2/100-1/20`; use a common pair `(-b,b)` only after the
actual leg-wise coefficients/zeros join to it. An unexpected radicand is a
stop with its operand, not an invented general radical solver. Profile xi
and non-radical variables retain their existing rules.

On each matching real Fourier leg the source branch is
`q=+sqrt(b**2-k**2)` inside the band and
`q=+i*sqrt(k**2-b**2)` outside it. Join this to the native real-axis Piecewise
and its physical scale. Split the finite momentum interval at both endpoints:

    interior: k=b*sin(theta), -pi/2<theta<pi/2;
    exterior: k=+/-b*cosh(u), 0<u<acosh(K/b).

Keep Jacobians, the negative-exterior orientation and the original nested
Fourier measures. Endpoints are panel bounds, not evaluation nodes. Map the
existing transfer/Abel panel centres from physical k into the transformed
coordinates, including centres dependent on outer legs. Source and profile
cutoffs remain finite; no interchange of infinite integrals is claimed.

The closed resolvents can cancel inverse-q factors; the middle leg of ordered
second scattering can retain an integrable 1/q singularity. That is a source-
reasoned candidate, not a completed row census. Inventory the actual bound
denominators and powers, inspect complete sums rather than isolated summands,
and save the transformed weighted values on two-sided sequences approaching
each endpoint. Keep physical lambda_A and all coupled denominators live in
this check. Stop on an uncancelled nonintegrable term, a coupled zero on the
physical integration path or unsupported source correspondence. No finite part,
dropped term, fabricated removable value or extra damping repairs it silently.

Independent checks are concrete:

1. First inventory the actual ordered variable/limit position of the middle
   radical in every three-momentum row. c2 appends it last, which makes it
   outermost under SymPy's limit convention; do not assume the reduced row
   keeps a different order. On at least one nonzero such row, compare the
   complete middle-leg integral with adaptive Gauss–Kronrod directly in k,
   split at radical zeros and Abel centres. The reference uses neither the
   sin/cosh map nor its weights. In the nested row comparison, replace the
   rule at the actual middle-leg position. If it is outermost, rebuild all
   dependent inner panels for every adaptive middle node. Shared remaining
   rules are disclosed and do not establish their independent accuracy.
   Tolerances are relative 1e-6 and absolute 1e-8 in the recorded coefficient
   units. Save operands and reference error estimates; refine the main rule
   separately. Unsupported native ordering stops instead of testing the
   wrong leg. This implements C4.
2. At zero contrast, apply the assembled uniform operator to Gaussian packets
   in a face-driving field direction (including eW), with momentum centres0
   and3b and width8 in the saved length frame. These are not claimed band-
   limited. Compare interior values with independent adaptive 1-D integration
   of the original REFERENCE pencil on the physical branch times the packet's
   Fourier transform. Keep normalization and component units. Account for
   finite-window derivative boundary terms, or verify their contribution is
   below the absolute comparison tolerance for every derivative used. Compare
   the actual finite-cutoff problem first; a cutoff change tests its tail.
   A zero transverse drive does not substitute for this responsive check.
   Use relative 1e-6 and absolute 1e-8 comparison tolerances. Inspect the actual
   zero-contrast operator for regulator dependence: either establish absence
   or use exactly the same regulator in the independent reference integral.
   A mismatch is a failed check, not silently removed as a regulator effect (C2).
3. Retain a one-sided wrong-sheet sign and omitted-Jacobian control, saving
   movements of actual operator values. Compare native term accumulation with
   assembled row matrices. Responsive mutations are not correctness proofs.

These are calculations in the guarded instrument, not extra solve campaigns.
An independent middle check or the wavepacket comparison that cannot be
implemented within this bounded numerical route is a method stop.

## Per-contrast end states and finite boundary

Use the original full rational pencil and wave relation, with invariant pairs
R,K,Q and whole transverse doublets. Saved candidates are 18 records per end;
the old selected set has five clusters/seven directions, of which five are
outgoing and two incoming. These counts do not prove the selected set stays
appropriate at changed contrast or frequency.

The continuation order is explicit:

1. Start from saved `omega=1,a=1` LEFT/RIGHT candidate records. Reuse their
   existing right bases, k, q and full-pencil checks; do not recalculate saved
   roots. For `a=1/2,1/4,0`, continue the same algebraic end pencil in a at real
   omega1 using invariant-pair Newton with maximum contrast step0.05 and
   comparison step0.025. At a0 compare against the independently saved
   REFERENCE pencil and candidate subspaces. Compare projectors, not printed
   expressions or vector phases. A collision, ambiguous correspondence or
   changed nullity is unresolved; no new root search follows.
2. For each completed contrast, transport **all 18 saved candidates**, not
   just the old selected five clusters, along positive-imaginary frequency
   paths of heights0.05 and0.025 to the target, then descend to real frequency.
   Horizontal steps are at most0.05 and0.025 respectively; descent steps at
   most0.01 and0.005. This path is for boundary selection only. Do not put its
   complex frequency into the real field operator or reported currents.
3. Preserve actual R,K,Q, full original-pencil/wave/commutator residuals,
   denominator margins, rank and conditioning at every completed point.
   Use saved derivatives where available, otherwise a secant predictor;
   no new general symbolic frequency derivative. Numerical corrections solve
   the existing invariant-pair equations.
4. Independently transport each paired radical along its actual joint
   (omega,k) path from the recorded seed, retaining sheet history rather than
   resetting a principal square root. Compare the corrected Q with that
   transported lift. Record physical Im(q), |q| and closure-denominator
   margins. A branch-locus encounter, ambiguous joint eigenpair/doublet lift,
   path dependence or selected-root sheet change makes BOUNDARY_UNRESOLVED.
   Im(q) alone is not a sheet classifier on a complex joint path.
   Any failed Newton correction, denominator margin or joint-lift check for
   any of the 18 candidates makes BOUNDARY_UNRESOLVED, including non-selected
   candidates. Sheet change means disagreement of corrected Q with the
   independently transported lift, not a change in Im(q) along a complex
   path. At the real target require Im(q)>0, or real q>0, on the physical
   branch; finite comparison tolerances and raw values are recorded (C3).

At each **real** target enumerate the continued candidates with their physical
q and spatial k. Report proper-depth/outgoing-radiation membership, spatial
decay/growth and current sign separately. For real k use the source's outgoing
q and actual current sign. For complex k require the selected boundary state
to decay towards that end and have a supported physical-depth sheet; a depth-
growing leaky-pole state is not silently used as a standalone physical mode.
If the five selected outgoing trace directions cannot be supported, or another
continued candidate changes their counted outgoing span, stop. Never discard
a direction until a 5x5 inverse works. This is a conservative eligibility check
on the saved candidate set, not proof that it exhausts all possible channels.

Require two incoming transverse directions per end with real k to numerical
tolerance, resolved current sign and vanishing face/bulk drives. Retain the
other outgoing directions for interior matching. Use the existing trace-map
structure, compare both frequency paths/steps and record full-pencil residuals.
Starting tolerances are scaled1e-9 for pencil/wave/trace residuals,1e-8 for
projector agreement and real-k/zero-drive checks; save raw values and scales.
Ill-conditioning that defeats these comparisons is unresolved.

This is an approximate finite modal end condition. The radiation continuum
outside the finite interval is not completed by 18-candidate continuation.
Domain variation probes that omission but does not bound it. Failure to obtain
an eligible central end map ends the pilot; do not construct a new transparent
boundary or extend the thickness-classification programme.

## Current, finite-end interference and retained-order interpretation

Evaluate the source-polarized matrices after numerical binding. Source
linearity in epsilon and bilinear contractions justify epsilon-squared
homogeneity; retain the saved coefficient residuals as provenance. Evaluate
at epsilon1 and check rescalings1/2 and2. Numeric samples alone are not an
identity proof. Unsupported dependence stops; no full Poly extraction.

Define the reported transverse current directly by projecting the **slab
current matrix** onto the actual selected transverse bases, conditional on
their zero loss-side drives. Apply orientation and incident normalization.
Check actual face lifts and projected bulk-density matrices; if the transverse
reduction is not supported, the reported balance is unresolved. This definition
never computes an undefined infinite-depth integral and multiplies by zero.
Do not add interface/port power to slab current: they enter different places
and have different units in the source balance. Radiating non-transverse
output is not included in J_T and needs no invented infinite-depth pairing.

At each finite matching end additionally report anchored slab-current cross
terms between the transverse and remaining modes, and between incoming and
outgoing transverse states. Use the actual solved amplitudes and phases.
Because a transverse leg has zero acoustic lift, the relevant bulk mixed
contraction must also vanish directly; otherwise that partition is unresolved.
Normalize cross contributions to incident current. If they exceed the larger
of the numerical envelope and10% of |D_T|, mark PARTITION_UNRESOLVED and do not
interpret the separate incident/reflected/transmitted currents as a resolved
loss partition. Preserve Hermitian congruence and a sign/off-diagonal control.

The saved rectangular retained model keeps1,eta,sigma,eta*sigma, not all pure
second grades. A finite solve of that model is not the separate retained
response series; neither should be misdescribed. Nonetheless, the source only
checks some balance identities after grade projection. This proposal supplies
no bound on an unprojected O(a**2) conservation defect. It does not restore
missing parent-order terms or infer their absence from quadratic scaling.

Therefore a numerically resolved D_T remains a **finite-model transverse
deficit, not an established physical leakage rate**. Report this limitation
next to the number, not only in an appendix. A sub-threshold nonzero-contrast
control could include real face absorption and would not bound this error at
omega3; no extra frequency or two extra solves is silently added. If a physical
loss value separated from retained-order error is required, this pilot does
not supply it and must stop rather than build a parent-order/power campaign.

## Finite list, controls and output

Before assembly record every selected real transverse |k_T| and bulk endpoint
b. Require |k_T|+0.5 inside the base momentum cutoff. Require at least four
points per transverse wavelength at the largest Chebyshev-node gap. Determine
one common base basis size for all contrasts at the frequency before solving:
the smallest odd N >=257 meeting that spacing gate at the maximum selected
|k_T| and source bound. Refinement uses odd N >=max(385,ceil(1.5*N_base)).
Record actual sizes and scales; failure to cover the spectrum stops instead of
accepting agreement between two under-resolved settings (C1). This is a
predeclared finite-discretization adjustment, not a change in physical inputs.

Starting base settings:257 coefficients/field, source/matching interval[-64,64],
profile interval[-14,14], momentum cutoff4, Abel0.1, source/profile orders
512/512, transformed outer order24 and inner orders8/8. These are starting
values, not asserted adequacy. Compile/bind one numerical case at a time and
bound quadrature caches; do not accumulate all cases in memory.

| Setting at each frequency | Contrast a | Purpose |
|---|---|---|
| Base | 0,1/4,1/2,1 | Uniform control and contrast scaling |
| Basis385, source/profile768/768, outer32, inner12/12, cutoff6 | 0,1 | Joint numerical refinement |
| Matching/source interval[-80,80], otherwise refined | 0,1 | Domain sensitivity |
| Abel0.05, otherwise larger-domain refined | 0,1 | Regulator sensitivity |

Maximum10 solves/frequency,30 total. Central omega3 comes first; measure its
base assembly/solve cost before committing the refined cases. All controls
retain physical permeability/memory. Abel is a numerical Fourier weight.
All completed case arrays and failed operands persist. No automatic retry.

For each incident column, retain raw signed D_T, raw matched D_T(0), their
difference (explicitly a contrast diagnostic), and every sensitivity. Let
E_num be computed separately for each incident column as the maximum of
1e-6 and all observed refinement/domain/regulator changes over a in {0,1}.
Apply that same per-column envelope to a=1/4 and 1/2; those contrasts do not
have separate refinement solves (C5). The uniform floor is |D_T(0)|, shown separately; the comparison
envelope is max(E_num, uniform floor). Neither is a physical error bound or
a known subtractable absorption floor.

- A uniform deficit outside1e-6 is explicitly flagged. With the selected
  lossless transverse premise this is a failed background control, not proof
  that merely opening the acoustic band makes uniform transverse light leak.
  Report the central data and stop before 2.3/4; do not hide it by raising
  a resolution label and proceeding (C5).
- Any deficit below minus its numerical comparison envelope is
  NEGATIVE_DEFICIT_UNRESOLVED. Preserve the negative value; do not clip or
  call it gain. Report the comparison against E_num as well.
- Report contrast-halving ratios only for resolved values. A ratio near4
  supports quadratic scaling but also fits certain artifacts. Use
  SCALING_UNRESOLVED otherwise; do not fit an exponent to numerical noise.
- FINITE_MODEL_DEFICIT_NUMERICALLY_RESOLVED requires a positive deficit over
  three times the comparison envelope, at most20% observed sensitivity,
  supported current partition and all source/boundary/integration controls.
  It does not remove RETAINED_ORDER_LOSS_INTERPRETATION_UNRESOLVED.
- Otherwise return actual finite values with the appropriate unresolved
  statuses. No result inherits the old omega1 reporting resolution by fiat.

Report central costs, peak resources and visited coverage. Frequencies2.3/4
are attempted only if the central numerical method/control gates pass; the
retained-order interpretation label alone does not conceal or erase numerical
coverage. Three frequencies do not give a threshold map or no-leakage region.

## Implementation authority and method stop

The user authorized proceeding without a further Grok pass. Preserve the
literal Claude verdict and missing Grok report; no paired method or build
clearance is asserted. C1–C5 above implement Claude's small pre-implementation
corrections. The prior reviews and exact drafts remain unchanged.

Execute the central candidate-continuation gate first, using only accepted
original end pencils and whole saved mode bases. Save all new contrast/path
operands and returns; reuse a completed result when the actual source is
unchanged. This first worker does not claim a face/current check, a full
boundary map, field assembly or deficit. If continuation passes, continue
implementation of the remaining current/face and integration/finite-solve
parts within this same user-authorized goal, using saved completed work.
If a required central method gate fails, report the concrete evidence and stop;
no new transparent boundary, automatic retry or exact-track restart follows.

Ordinary shared guard around the existing supervisor, native limits, one
scientific worker and silent completion/error hook remain mandatory. Per job:
900s outer, 840s native, 2GiB, zero swap, one CPU, 32 tasks and one native
thread. No scheduler changes or inherited exact-track resource exceptions.
The user decides the overall authoring duration; the former 16-hour automatic
stop no longer applies. Pure tooling corrections get local tests and rationale.
Scratch is runtime storage only; canonical code/status is committed outside it.
Stop with actual finite-model results and limitations or a concrete method
blocker. No physical leakage, A11/A12, full Green/FORM or optical calibration
claim is supplied by a successful continuation or finite solve alone.
