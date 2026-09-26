# A9 construction boundary: decision before a new solver

Status: **source-grounded dependency finding and bounded proposal**. No new
physics, independent clearance or change to the accepted Option B scope.
The document review is finished. This is the concrete stop anticipated by
the amendment, not a request for a tenth amendment review.

## What the implementation inspection established

The existing finite response is reusable on its recorded domain. It is not
yet the general radiating response needed to publish all retained FORM roots.
Two specific dependencies now replace the earlier vague “general FORM” task.
The [inspection record](S11c_d_FORM_boundary_source_inspection.json) pins the
whole source functions and exact local excerpts used below. No scientific
module was imported and no saved calculation was repeated.

**1. The slab-to-face map needs a centre-motion disposition.** S11c-a supplies
separate `DELTA_W` and `ZETA_C` source directions. Its
`build_material_face_source` sets the centre contribution to zero when
constructing the thickness direction; this is an individual source direction,
not a computed elimination of the independent centre degree of freedom.
S11c-b computes and retains `CENTER_FACE_GENERALIZED_ROW` beside the slab
rows. The actual c2 caller sets `REPRESENTATION = 'DELTA_W'`, consumes that
face-velocity direction, and closes only `U`, `THETA` and `E_W` in
`expanded_rows`/`build_case`. d's `ReducedPencil` has the resulting five fields.
The eight previously indexed saved c2 velocity identifications confirm this
particular route; their values are not being recomputed here.

The alternative d `ClosedAcousticEnergy.construct` does not supply the missing
general centre map: its per-face displacement is explicitly the signed
thickness amplitude. It is an end construction, not a profile-wide exterior
solution. The inherited b/c1 specifications retain an independent centre DOF;
the accepted amendment requires an actual elimination/face-drive map and
forbids substituting a parity assumption. **No centre elimination or complete
invariant-sector argument has been established by this inspected chain.**
This is not proof that none exists elsewhere, that the centre must be excited,
or that the computed five-field results are numerically wrong.

**2. The available end closure has a decaying-bulk domain.** In the actual
`continuum_boundary.construct_end`, modal selection requires
`SHEET_MEMBERSHIP` and `BULK_DECAY_DISK_CERTIFIED`; the current-pair loop requires
positive imaginary bulk momenta on both legs. It constructs five outgoing
trace columns and two incident columns at each end. `continuum_response`
then assembles a finite coefficient system and four incident right-hand sides.
Those are actual domain and representation choices, not general channel
counts. `EndResolventAudit` constructs fixed-frequency inverse-pencil/local
Laurent data; that caller does not build a profile-wide outgoing Green kernel.

Consequently these inspected constructors do not establish A12 on nonempty
radiating support, nor continuous-spectrum coverage for a general outgoing
FORM. The located c1 far-field construction is useful existing evidence, but
its restricted drive and first-shape power expression do not close either
dependency. This finding does not invalidate a saved result on its actual
domain or prescribe a nonzero-leakage example.

## Concrete construction proposal and the limit of that proposal

If the five-field sector is justified, a candidate generic construction is
the full reduced reference-symbol Green kernel on its actual outgoing sheet,
followed by the retained independent-grade response recursion using the saved
local/nonlocal profile kernels. The baseline symbol must keep the computed
off-diagonal couplings. End-mode, incident-source, observation and current
normalization variations must enter the same recursion. The exported inverse
must be constructed from the actual symbol, with its contour/domain and
continuous contribution; an opaque inverse name is not an implementation.
Nondecaying profile endpoints also require the actual asymptotic source/phase
terms rather than treating the whole perturbation as compactly supported.

This is a **method proposal, not a computed Green kernel or established domain**.
It is new work beyond wrapping the saved finite solver. A finite-basis formula
would be an explicitly approximate alternative and would still need the
missing radiation coverage; changing the representation does not waive it.
Neither method is launched or pre-approved by this report.

## Recommended next authorization: one bounded dependency resolution

Prepare a focused b/c1/c2-to-d compatibility task before commissioning the
general radiating solver. Its sole deliverable is a disposition of the centre
drive for the actual retained scattering sector:

1. Trace the saved centre generalized-force operand and both independent face
   source directions through the supplied kinetic, constraint and closure
   data. Determine whether there is a saved complete reduction, a separately
   prescribed input, or an unresolved upstream equation. A zero face-force
   term alone is not a complete centre dynamics or decoupling proof.
2. If a genuinely missing algebraic compatibility check suffices, specify it
   from those exact operands for the four cases and independent retained
   grades. Reuse every completed source/closure operation. Persist both
   operands, any residual and the actual assumptions. Do not add a mode,
   impose centre zero or rebuild any response to make the check pass.
3. Stop with either the usable reconstruction/sector restriction and its
   evidence, or a precise upstream input/repair request and affected claims.
   Do not turn that result into automatic authority for a new radiation method.

**Bound and rough cost:** medium authoring/review effort for this focused
dependency task; at most one new guarded algebraic job (900 s, 2 GiB, one CPU)
after the applicable build gate, with no automatic retry or numerical solve.
The time cap is a budget, not a predicted completion time. If a missing
equation or new degree of freedom is required, stop within this task instead
of beginning a full upstream rebuild. No new external review has been sent.

The subsequent generic radiating response is **large/uncertain** work: an
outgoing representation and supported domain, retained response/observable
construction, then its independent/Wolfram/T7 checks. A reliable total runtime
cannot be inferred from the inexpensive existing finite solves. Approving the
bounded dependency task need not authorize that larger method.

Publishing only the existing numerical annex would be smaller, but would not
fulfil A9/A11/A12 under the accepted amendment. It would require an explicit
scope change; none is proposed as a silent substitute here. Pole deferrals,
all saved results and failure histories remain unchanged. A6/A7 and the
located first-jet routes stay resolved; no repeat mapping campaign is queued.

The reason for the stop is the accepted amendment's explicit rule: “If no
upstream map determines a required face drive, record an upstream c1/c2
dependency finding and its implications for any affected end-channel result;
do not close the gap by adding or fixing a mode in d.” Its separate
radiation-method decision rule also applies. This report records an unresolved
dependency, not a finding that every earlier result must be regenerated.
