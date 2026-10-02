# Defect operator applicability and the next bounded correction

Prepared 2026-10-01 after the user renewed continued work. This is a proposed
method and source-based decision, not a new calculation or independent clearance.
The objective remains a useful rest-bulk near-unity defect test, with the drain
and primitive calibration limitations stated. Repeating completed uniform or
selected mixed-term calculations is unnecessary.

## What is established, and what remains different

The selected uniform near-unity schedule completed on both sides of each end's
actual modal/acoustic match. Its selected transverse states and currents were
finite, with zero selected face drives and joined grazing limits. That result
concerns constant ends, not propagation through the profile.

The upper-face tanh diagnostic established a nonzero direct mixed bare action
at profile momenta 0 to 1/10, omega=3, c_s0=10. The saved boundary coefficient
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

1. Reuse the saved generic upper boundary coefficient and its actual modes,
   input equations and native graph/linear joins. Derive the lower-face change
   from the actual native face normal, location, pressure convention and outgoing
   depth sign. A guessed face-sign copy is insufficient. Join zero/first-order
   coefficients to native c1 on each face; require the original boundary
   residual before accepting the new coefficient. Preserve eta and sigma_W as
   independent grades and retain the specified rectangle, not full pure second
   order. Reconstruct only a genuinely missing face/argument calculation.
2. Impose one physical outgoing dispersion on input, both intermediate routes
   and output, with the conserved edge variables joined. Insert both assignments
   of height and slope to the two transfers using the saved native transform
   convention. At general input momentum, rederive the delta/contact and
   principal-value contributions: the special k_in=0 cancellation is not a
   universal license to discard them. The two assignments must not double a
   coefficient already generated by the ordered expansion. Retain explicit
   distributional and domain conditions instead of numerically sampling them.
3. Keep the result as a declared sum of raw middle-momentum integrands and any
   separately supported contact terms. Native `kernel_apply(second=...)`
   integrates its second slot over the middle leg. An already-whole convolution
   cannot be fed unchanged into that slot. The selected inherited whole D can
   be used only as a comparison at its original point, without recomputing that
   completed integral or treating it as the arbitrary-momentum kernel.
4. Add the direct contribution once to each face's native ordered response,
   retaining the existing first-shape iteration exactly once. Reuse the complete
   original source/density/memory/epsilon and physical/reference-pressure maps.
   Contract into the actual full slab rows before grade restriction. Verify
   units, face orientation and source/row addresses. Keep pressure and normal
   jet changes separate; do not infer their generic behavior from the selected
   zero-grade jet coefficient alone.
5. At source/operand level, identify exactly which local/nonlocal coefficient
   records change. Establish actual constant-end/zero-jet specialization so
   unchanged end modes/currents can be reused. A source-file hash change alone
   does not require replay. Compare the selected new increment to saved evidence
   and persist the difference, not a regenerated baseline pipeline. Addressed
   direct-slot, slope, sheet, face-orientation and single-integration controls
   must be responsive where applicable. Controls establish sensitivity rather
   than physical acceptance.

This first construction ends at the incremental operator and its applicability
evidence. It does not include a defect sweep or claim a new loss number. Before
changing finite matrices, require an explicit integration recipe for every new
contact/PV/branch piece. In particular, endpoint integrability at c_s0=10 and the
selected k_in=0 witness does not automatically cover grazing/branch coincidence
near modal matching. A new nonintegrable or unjoined term is a real method
finding, not an occasion to loosen a numerical threshold.

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

No new science has run during this preparation. Scientific work will use the
existing no-deadline pooled guard/supervisor, zero swap, native memory limits,
disjoint CPU assignments and aggregate reservation, with the local completion
hook armed first. No replay, time cutoff, scheduler change, automatic retry or
new external export is authorized by a process exit. Costs of general contact/
endpoint handling remain uncertain; the five-second selected continuation is
not a runtime estimate for this construction.
