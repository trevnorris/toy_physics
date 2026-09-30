# Bounded numerical radiating pilot — method for independent review

2026-09-29 local. **Proposed method, not clearance or computed results.**
The user has directed a numerical pilot with two independent method reviews
and a hard stop if boundary/integration authoring becomes a larger project.
The exact-symbolic omega=1 premise track stays parked. This document makes
the numerical choices concrete for that one review round before construction.
No scientific payload has been restored in preparing it.

## Object, scope and stopping point

Compute the finite-model transverse-current deficit for LAB_HELD/RHO4_CONSTANT
at the saved tangents `(1/5,1/10)` and frequencies **3, then 2.3, then 4**.
All material constants, permeability/memory, profiles and reference units stay
as in `S11c_d_variable_profile_development_input.json`. The lower point is
just above the recorded bulk threshold `sqrt(5)`; it is not at grazing.
No claim maps these frequencies to a calibrated optical band.

Let `a` scale the complete saved profile contrast:
`eta_bg=a/100`, `sigma_W=eta_bg*W_0/L_W`, with `a=0,1/4,1/2,1`.
Keep independent eta/sigma dependence in the saved operands until this
numerical binding. Scale both the saved thickness and density-profile
perturbations as their original source requires; do not replace the problem
by thickness-only forcing. At zero contrast the same operator assembly binds
those amplitudes to zero, yielding the reference uniform background. It does
not delete difficult rows by hand or change the physical face permeability.

For each of the four transverse incident columns, evaluate

    D_T = 1 - (J_T,reflected + J_T,transmitted) / J_T,incident.

Both output directions count as surviving transverse light. Use physical
quadratic current matrices and their cross terms, not squared field amplitudes
or a first-order transmission deficit. This is a direct finite solve of the
saved retained operator at finite contrast, not a Born-rate construction.
It inherits the operator's retained-order limitations.

The deliverable is a signed value, numerical sensitivity and status at each
attempted point. It does not partition radiation, face absorption and
non-transverse conversion; certify total physical loss; close A11/A12/FORM;
or claim that availability of bulk radiation guarantees significant leakage.
If an essential method issue remains after two authoring days, stop. If the
central numerical case is unresolved, return it rather than start the other
frequencies, more tooling cycles or an exterior simulation project.

## Reuse without replay

Use the saved frequency-source packet's `liveFrequencyAndGrades` records,
addresses, original ordered integrals and unit maps. Its 375 records cover
the local matrices, cells, factors and 35 source amplitudes of all 80 nonlocal
rows. Bind numerical frequency and material/tangent values before new
numerical compilation. No frequency derivative, coefficient-domain Poly,
epsilon-carrier extraction, exact source producer or old accepted calculation
is required for the new field solve.

Use the existing Chebyshev collocation, finite source/profile integration,
linear solve, residual and physical-current machinery. Reassemble matrices
at each changed frequency/contrast; old omega=1 matrices are not those
matrices. Completed numerical arrays may be reused only when their actual
inputs agree. The saved frequency-live end pencils retain their original
algebraic sources and mappings; rebuild a changed-contrast numerical binding
from those, not the already specialized right-end matrix.

The source chart with `abs(omega-1)<=1/4` is not used outside that domain.
In particular its `i*sqrt(-radicand)` replacement and negative-radicand proofs
are not transported to these real radiating frequencies. Native real-frequency
branch expressions, physical face closure and units are the governing inputs.

## Bulk radiation and branch-endpoint integration

The native bulk elimination already contains an outgoing acoustic boundary
relation. Its flat impedance is `rho_m*omega/q`; c2 supplies the real-frequency
outgoing/decaying root on each of the output, input and middle momentum legs.
The pilot adds **no sponge, physical damping change or artificial complex
frequency to the real-frequency field operator**.

At the saved tangents define `b^2=omega^2/100-1/20>0`. For a real Fourier
normal momentum `k`, use the original root `q=sqrt(b^2-k^2)` with positive
real value inside the acoustic band and positive imaginary value outside.
The source uses both physical depth momentum and a scaled radical: retain
the original scale map, rather than identifying their units or coefficients.

For every momentum leg, split the actual finite integration interval at
`-b,+b`. Proposed changes of variables are

    interior: k=b*sin(theta), -pi/2<theta<pi/2;
    exterior: k=+/- b*cosh(u), 0<u<acosh(K/b).

Weights include the actual Jacobian and reversed orientation of the negative
exterior interval. This removes the elementary 1/q endpoint singularity.
Use interior Gauss nodes, retaining exact endpoint locations as panel bounds.
Map existing transfer-centered Abel panels into these coordinates, including
centers depending on already fixed outer legs. Preserve the native nesting,
Fourier factors and integration measures. The source/profile quadratures
remain finite; this is not an interchange of infinite integrals or an Abel
limit. Increase panel order for a separate numerical comparison.

**The 1/q premise needs review and runtime evidence.** The c1 port closure
has `R(q)=1/(1+Lambda_A*rho_m*omega/(rho_m^2*q))`. For nonzero physical
permeability it can cancel exterior-leg inverse-q factors in the closed
kernel; c2 also has an ordered three-leg second scattering term. This is
motivation, not proof that all assembled native coefficients have only the
allowed singularity. Review the supplied `port_matrix` and `kernel_bridge`
source. In construction, inventory every native denominator/power on each
leg after physical binding, preserve full sums, and check source-specific
cancellation and the transformed weighted integrand. Isolated summands may
have worse behavior than the complete coefficient.

Use finite numerical endpoint approaches on both sides at multiple distances
and compare quadrature refinements, saving the actual weighted arrays. These
are bounded numerical diagnostics, not a general integrability theorem. An
uncancelled nonintegrable endpoint, a denominator zero away from the branch
endpoints, or an unsupported source identity stops the pilot. No finite-part
prescription, contour deformation of the field integral, dropped singular
row, made-up removable value or extra regulator is silently introduced.

Controls: compare physical q against the original source's real-axis branch;
retain a wrong-sheet sign control in the bulk operand and an omitted-Jacobian
control in quadrature. Save movements in actual matrices/current outputs
where applicable. An identically zero transverse drive alone is not a
responsive radiation-boundary control. Independent scalar quadrature with
the same source integrand checks selected endpoint-bearing row actions.

## Finite end boundaries and transverse current

Use numerical continuation of the saved full end invariant subspaces, keeping
transverse doublets together. The existing R,K,Q equations, full original
pencil residual and five-dimensional trace map are starting structure. Keep
all five outgoing field directions, not just transverse waves; the other
directions affect the interior field even though the headline observable
counts transverse output only. No complete thickness-channel map is required.

Proposed outgoing selection: continue the already selected omega=1 subspaces
on a positive-imaginary-frequency path to each target and then to its real
boundary value. Compare two path heights (0.05 and 0.025 in reference units)
and step sizes (0.05 and 0.025); coefficients use the full rational end pencil
with physical relaxation retained. The auxiliary path is for numerical end
selection only; the field solve and reported currents have real frequency.
The positive-half-plane prescription itself is a method-review question,
not an established global analyticity assertion for the whole pencil.

Use saved frequency derivatives where already available; otherwise use a
numerical secant predictor, with corrected residuals and a smaller-step
comparison. Do not generate new general symbolic frequency derivatives for
the predictor. All changed-contrast end bindings precede numerical work.

At every retained point save actual K,Q,R, wave/full-pencil residuals,
denominator margins, subspace dimension, conditioning and branch history.
Compare subspace projectors and final trace maps rather than individual
basis-vector phases. The numerical chart must not reset q to a principal
root along a complex joint path. Step failure, an unresolved collision,
dimension change or material path dependence stops; no new root census or
exact discriminant campaign is opened.

At each real target/contrast require an independent full-pencil check and
the expected transverse subspace, real normal momenta within declared
numerical tolerance, resolved physical current sign and vanishing loss-side
end drives at that numerical tolerance. Complex non-transverse roots are
not labelled closed. Require their selected spatial behavior to be outgoing
under the reviewed continuation; a growing/ambiguous direction is unresolved,
not removed until a 5x5 inverse happens to work.

Construct the finite trace condition from those complete numerical bases
using the existing `frequency_end.maps` structure. This is an approximate
nonlocal end condition, not a transparent-boundary theorem. Check domain
sensitivity and a source-level boundary-selection mutation. Do not reuse
the old `BULK_DECAY_DISK_CERTIFIED` filter as the radiating selector.

For the selected transverse inputs/outputs, evaluate the original physical
current form after numerical binding and project onto the actual bases.
Use its existing polarized bilinear construction to justify removing the
common epsilon_shape-squared amplitude normalization. Evaluate this known
quadratic form at epsilon=1 and check its 1/2 and 2 rescalings numerically;
do not construct a general Poly coefficient domain. The independent review
must check that homogeneity against the supplied current source, rather than
treating three numeric samples as a polynomial-identity proof. Face rows can
similarly be evaluated on amplitude basis vectors with saved linearity checks.
If the actual payload has unsupported amplitude dependence, stop rather
than resume the old exact carrier-extraction track.
If loss-side face/bulk amplitudes vanish to tolerance, justify the remaining
slab-current contraction directly; never form an undefined infinite-depth
integral and multiply it by a nominal zero. If this reduction is unsupported,
stop. Apply Hermitian current congruence, preserve cross terms and incident
normalization, and test an addressed sign/off-diagonal mutation. The old
omega=1 current normalizations and labels are not transferred unchanged.

## Finite experiment and acceptance labels

Proposed starting settings: 257 coefficients per field, source interval
[-64,64], profile interval [-14,14], momentum cutoff 4, Abel regulator 0.1,
source/profile orders 512/512, outer momentum order 24 and inner panel orders
8/8 in the transformed coordinates. These are planning choices, not measured
adequacy. The larger basis allows for the higher frequency than the old
129-coefficient omega=1 runs. No assertion guarantees that it resolves omega=4.

For each frequency, the finite list is:

| Setting | Contrast multipliers a | Purpose |
|---|---|---|
| Base | 0, 1/4, 1/2, 1 | Uniform control and small-contrast scaling |
| Joint numerical refinement: basis 385, source/profile 768/768, outer 32, inner 12/12, cutoff 6 | 0, 1 | Check discretization, endpoint integration and finite momentum extent |
| Domain/source interval [-80,80], otherwise refined | 0, 1 | Expose finite end/source-window sensitivity |
| Abel regulator 0.05 on the larger-domain refined setting | 0, 1 | Expose regulator sensitivity |

That is at most **10 solves per frequency, 30 total**, sequentially, rather
than the earlier 6–10-solve sketch. The added count implements the user's
requested controls at every frequency; it is one reviewed instrument/experiment,
not thirty preparation/review stages. Central omega=3 comes first. A central
method, sign, convergence or uniform-control failure stops later frequencies.
This is not a Cartesian parameter sweep or an automatic retry policy.

Report every raw signed deficit, the full matrices/vectors behind its current
contraction, and the matched uniform controls. Define a displayed empirical
resolution floor per incident column from the declared `1e-6` floor, largest
absolute matched-uniform deficit, and observed refinement/domain/regulator
changes. Keep those components separate as well as the envelope. The floor
is not a certified bound and does not automatically subtract from D_T.

Use the following reporting rules, subject to the independent method review:

- A uniform deficit outside `1e-6` is flagged; do not accept a profile-loss
  claim merely by redefining the floor to hide it. Explain a larger measured
  floor and stop with unresolved if it defeats the signal.
- Flag any deficit below minus its resolution as `NEGATIVE_DEFICIT_UNRESOLVED`.
  Do not clip it, call it loss, or attribute it to gain without a separate
  physical explanation. Numerical controls that fail remain visible.
- Show raw D_T at all three nonzero contrasts, and also differences from
  their common zero-contrast control, with separate labels. Report ratios
  against the expected quadratic scaling only when both values are resolved.
  A ratio near four on halving is supportive of leading quadratic behavior;
  artifacts can also scale quadratically. Do not force that law or use it
  as a substitute for boundary/regulator checks. Unresolved small-contrast
  points give `SCALING_UNRESOLVED`, not a measured exponent.
- A positive deficit is called numerically resolved only if it clearly
  exceeds the empirical envelope and moves little under the finite checks.
  Propose a conservative reporting margin of three times that envelope and
  no more than 20% sensitivity of the resolved deficit. These are practical
  pilot criteria, not physics laws or guaranteed coverage. Otherwise report
  the signed finite results with an unresolved loss interpretation.

The uniform comparison is a numerical/background diagnostic, not proof of
the regulator's absolute error; profile-dependent artifacts need not cancel.
The three selected frequencies do not establish a threshold curve or a
no-leakage region. No resolution achieved at omega=1 is inherited at 3.

## Implementation, review and budget

Fresh Claude and Grok receive this method, the same source excerpts and
identical questions, with no peer reports or external commentary. This is
one method-review round. It covers the integration prescription, auxiliary
end-selection path, physical current and experiment/claim criteria before
implementation. Clearance is not an executed result. A substantive blocker
returns a stop/decision; no automatic review-until-clear campaign is planned.

Implementation must follow the cleared mathematics, with the numerical
functions, inputs and checks tied to that review. Any substantive departure
is a method change, not a tooling exemption. Ordinary loader/JSON/format
fixes get local tests and a short rationale; they neither require separate
review rounds nor turn failed physics/integrity evidence into acceptance.
Do not compare symbolic display strings as mathematical identity. Keep raw
input hashes and actual numerical/source provenance instead. Save arrays and
case-level checkpoints without thousands of per-scalar rendering receipts.

One scientific worker at a time under the existing guard/supervisor, with
native limits, complete failure preservation and the silent local completion
hook. Split the finite list into bounded sequential batches as necessary;
completed cases are not rerun. A timeout/failed scientific case does not
automatically authorize a retry. No scheduler changes. No new scientific
run is started during method-packet preparation or review.

The two-authoring-day hard stop is a maximum of 16 active authoring hours;
review/user wait does not count as authoring. Record active work in the one
preparation record, do not create an unattended task to use up the budget.
The first authoring session began 2026-09-29 local; at the first durable UTC
mark, 2026-09-30T02:25:17Z, charge 0.5 h already used. Stop sooner if the source
shows that the simple endpoint/end-map route is unsupported. The earlier
1–4-hour compute estimate is uncertain and predated the full 30-solve control
list; do not promise that runtime. Measure the first central cases before
spending the remaining finite budget. Stop with data and limits, not a new
exact-method or simulation project.
