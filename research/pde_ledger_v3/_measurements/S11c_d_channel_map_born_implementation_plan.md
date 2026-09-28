# Channel map and conditional Born: proposed implementation plan

2026-09-28. **Prepared for user go/no-go; nothing has been implemented or run.**
Both fresh round-2 reviewers returned CLEAR FOR THIS PREPARATION AMENDMENT on
the unchanged [amendment](../directives/S11c_d_CHANNEL_MAP_BORN_AMENDMENT_DRAFT.md).
Their [literal reports and receipts](S11c_d_channel_map_born_review_r2_record.json)
clear that scope only, not future workers, currents, results or historical builds.
This plan specifies work bounds within that amendment; it is not a new governing
physics prescription or launch gate.

The pilot stays LAB_HELD/RHO4_CONSTANT, with saved material/profile/damping inputs.
It asks where channels are available, then attempts weak end conversion where
supported. Bulk availability and bulk power remain separate. No full FORM,
four-case response, centre project or A11/A12 completion is promised.

## (a) Map at the saved azimuth

Use the original accepted `uniform-source.pickle`/strong symbols with frequency
and both tangential components live. The compact frequency-live pencils froze
tangential momentum and are only the fixed-momentum fallback. Parameterize
`(k1,k2) = (2*kappa/sqrt(5), kappa/sqrt(5))`, where `kappa = |k_parallel|`, with
actual source units. This is the saved 2:1 direction, not all azimuths. No new
symmetry campaign is required.

Proposed window in the saved `L_ref,T_ref,M_ref` coefficient frame:
**0.1 ≤ omega ≤ 4; 0 ≤ kappa ≤ 0.4**. It includes the saved frequency-one point
and the frequency-three feasibility question, without claiming all relevant
frequencies are covered. Initial 16 frequency rows:
`[.1,.2,.25,.275,.3,.4,.5,.75,1,1.5,2,2.25,2.5,3,3.5,4]`;
12 momentum columns:
`[0,.025,.05,.075,.1,.15,.2,sqrt(1/20),.25,.3,.35,.4]`.
Allow at most 64 additional transition points: **256 parameter pairs maximum**,
each with REFERENCE/LEFT/RIGHT. Zero-momentum direction/degeneracy gets its own
label. Sampled claims stay sampled; no interpolated no-leakage theorem.

Before evaluating points, emit source/unit/endpoint/origin joins and the actual
`COMPUTED_BRANCH_BINDINGS`; check for a hidden saved tangential substitution or
sheet restriction. Preserve each original pencil and radical relation. Reuse
accepted seed results at identical settings, not completed mode jobs. Share
LEFT/REFERENCE work only after a source/origin join. No producer replay.

At each visited point retain candidate roots, denominator/opposite-sheet
exclusions, original full-pencil/nullspace residuals, complete degenerate
blocks, physical current and computed sector classification. A root count or
field-coordinate norm is insufficient. Keep physical damping. Unsupported
sheet, current, rank, threshold or leaky classification is **unresolved**, not
zero or lossless. A complete per-point root census and coverage between points
are distinct. No complex-pole search.

Within the region with an admissible incident transverse channel, report two
separate layers: thickness-like channels at each end and incidence; acoustic
depth availability with its allowed output normal-momentum interval. Keep
excluded points and reasons. Frequency and edge-parallel momentum are conserved;
edge-normal transfer is allowed. Bulk availability is reported even without
open end-thickness channels.

**Cost stop:** measure at most the first 12 new point evaluations, saving each
return. Continue only within the same 900-second job. If 2D work is too costly,
prioritize the saved-momentum column within the remaining budget and label that
restriction. No automatic second grid, wider window or longer runtime. If a
physical-current classifier needs a new unsupported method, preserve the
candidate/availability evidence and return that blocker instead of calling it
a completed physical-channel map.

## (b) Test actual reference transverse face drive

In the map consumer, join reference modes to the saved c2 `REFERENCE_TRACE_MAP`,
`REFERENCE_PRESSURE`, `IDENTIFICATIONS`, `NORMAL_JET` and both face orientations.
Use [the exact dependency routes](S11c_d_A9_saved_dependency_routes.json) and
[eight-face index](S11c_d_A9_saved_face_drives.json). Emit operands, units and
contractions, with an addressed physical-term omission or channel change that
makes the check responsive. A transverse/shear label does not establish zero.

Zero reference drive does **not** remove first-order induced thickness fields,
face-map changes or tilt terms. Closed/evanescent end-thickness channels do not
imply zero bulk radiation. Missing centre/exterior dependencies limit the
observable; they do not start a centre-mechanics project.

Feed this check back into (a). A uniform mode's exterior root must be evaluated
at its **own normal momentum** as well as its frequency and conserved momentum.
Some radiating output momentum being available to the edge does not alone make
that uniform mode leaky. Nonzero drive into its supported propagating exterior
channel must agree with the full-pencil/current classification, or be labeled
leaky/unresolved. This finalizes a provisional map; it is not a new radiation
method or a complete bulk-power calculation.

## (c) Born end conversion only on supported open-thickness points

If (a)/(b) supply regular flux-carrying channels, use **at most two** open-thickness
parameter pairs, selected by reported threshold/gap/current margins rather than
a desired nonzero answer. Include all incoming transverse polarizations from
**both incident ends** and reflected/transmitted thickness-like channels at
**both outgoing ends**. Closed modes remain matching data.

Consume the complete own reduced coupling, including supported zero/first-jet
and mixed-grade content, phases, transforms and physical current normalization.
Keep eta/sigma_W independent before the saved homotopy. Keep profile-width and
momentum-transfer dependence explicit; no width sweep, uncontrolled sharp-step
limit or presumed exponential law. Old evanescent-domain currents are not an
open-bulk normalization. Save operands before checks, including an actual
term-omission control and the reduced/full residual.

Report four separate statuses: symbolic baseline zeros, regular-domain checks,
reduced/full comparison, open-domain flux validation. Baseline end zeros do not
mean zero interior coupling. Leading power order follows the actual expansion;
a nonzero eta-squared coefficient is not supplied. If no useful anchor exists,
use **unanchored leading-order estimate** with every missing status. The label
and qualifications travel with all handoffs; they do not discharge §2 or A11/A12.

Bulk power stays pending without the actual outgoing exterior solution on both
half-spaces from the same solved face state, signed far-field flux, measure and
incident denominator. The restricted c1 identity and slab-current deficit do
not supply it. No full radiating solve follows automatically if Born fails.

## (d) Inspect saved contrast vectors before proposing new anchor runs

For this preparation, **16 selected saved artifacts and their checkpoint checks
were hash-verified; no scientific payload was restored**. Exact routes/hashes
are in the round-2 record. No new anchor solve is proposed.

A later guarded saved-output consumer uses:

- Response `coefficient-solutions.pickle` / `solved['coefficients']` and matching
  field/trace maps; `formal-remainders.pickle` contains `direct`, `retained`,
  `difference`, `directEquationResidual` at the three existing contrast pairs.
- `channel-selectors.pickle` and `continuum-currents.pickle → closedAmplitudes`
  identify matching thickness components and the possible baseline boundary floor.

Inspect all saved incident columns; compare matched first-order coefficients
before tiny physical amplitudes. Restrict direct/remainder vectors to the same
observable: their full-vector maxima are not thickness scaling or a Born/full
residual. Preserve profile, units, regulator, finite boundaries and modal
conventions on both sides. Aim for roughly 1% only on a resolved nonzero quantity,
with a declared absolute floor near zero. The old 1e-4 physical-amplitude floor
is not a coefficient arithmetic-error bound. Attainable comparison accuracy is
unmeasured; reviewers' coefficient-size/precision suggestions are not new results.
A closed-domain match tests assembly, not open flux or bulk power. If no useful
comparison is possible, keep the unanchored label and identify the mismatch.
No saved solve replay or new contrast values; a new open-point anchor needs a
separately priced proposal after this inspection.

## Proposed budget and stopping point

**Recommend go for the map/face implementation, with Born conditional.** Each
new instrument retains applicable independent build review before execution;
amendment review does not pay that gate. No implementation review is submitted.

Planning estimate, not measured: **1–2 working days** to build/review the map and
face test; **another 1–3 days** for supported Born/saved-anchor work. Confidence
is low until the first pilot measures cost and current-domain support. A concrete
method blocker stops the pilot rather than expanding that estimate indefinitely.

Later proposed compute ceiling: one map/face job, one conditional Born/saved-anchor
job, and at most one saved-output validator per completed job if needed —
**four ordinary 15-minute job envelopes maximum**, not a completion guarantee.
Review elapsed time is separate. No retry or deadline extension is included.

Use the unchanged shared guard around the normalization supervisor: 900 s outer,
840 s native, 2 GiB, zero swap, one CPU, nice15, 32 tasks, one native thread;
one job at a time, hook first. Preserve operands/returns, strict stderr and
stdout/checks identity, enforced resources and posthashes. No producer replay,
general certification or optional-review campaign. Scratch stays uncommitted;
canonical sources/reports live outside it. Lean/S11_lean/shared guard, protected
suffix and incident/failure history remain untouched.

**Stop here for user go/no-go. No worker or scientific job has been launched.**
