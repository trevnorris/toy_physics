# Next response step: the saved localized thickness forcing

Status: **concrete implementation target from inspected saved operands**.
The user directed continued work after the end-lift/forcing stage. The guarded
[saved-payload inspection](../_measurements/S11c_d_two_asymptote_saved_inspection_report.md)
is complete. This plan selects the smallest useful next response calculation;
it is not an external submission, scientific launch receipt or new acceptance
claim. The implementation and exact method-review packet must exist before
any new submission. No new external submission has been authorized.

## Why this is the next calculation

The two saved thickness-grade `(0,1)` cases, blocks 16 and 17, have **zero end
lifts and no profile Abel regulator**. Each has 36 retained profile-forcing
records, with complete source terms and ordered integrals saved. The density
grade `(1,0)` instead has the two nonzero RIGHT end fields and 18 Abel-bearing
source records per case. Starting with thickness therefore produces a useful
part of the retained response without making the unresolved density-grade
half-line/plane-jet issue a prerequisite for every coefficient.

This changes execution order, not retained scope. Density and mixed grades,
other cases, physical incoming/current normalization, general bindable FORM
and supported A11/A12 duties remain open. The unchanged slice is evanescent;
this calculation is not a search for a radiating witness.

## One bounded deliverable

At the unchanged LAB_HELD/RHO4_CONSTANT input, construct the **two thickness
forcing transforms and their action under the existing reference outgoing
prescription**, in the saved redundant residue-column frame. Reuse the full
5-by-5 reference prescription, branch, units and signed-current pole blocks.
Do not replace it with a scalar acoustic inverse or diagonal approximation.

Since the saved thickness lift is zero, the target is `L0 V01 = F01`, with
`F01` equal to the already saved full native forcing. The deliverable should
retain its continuous Fourier contribution and explicit inherited pole terms.
An unevaluated but properly defined source-specific Fourier/PV representation
is useful; it must not be labelled as a numerically integrated response or a
flux-normalized scattering coefficient.

## Implementation work, kept together

1. Restore only the two saved thickness forcing bundles, relevant native
   operands/units, residue-frame data and accepted outgoing-prescription
   artifacts. Transfer their existing checks through exact saved-object
   joins. Do not reassemble producers, end fields, modes, roots or inverses.
2. Form a small transform dictionary for the **actual localized tanh-derived
   factors present in these operands**. Retain the source coordinate scale,
   phase and Fourier mass. Translate the native plane-character/source
   integrals into momentum forcing with all local and nonlocal terms kept.
   The implementation must state the distributional identities used and
   justify them on this fixed localized source; it cannot silently reverse
   the original ordered integrals because an ordinary Fubini argument would
   be convenient. Save the original term, transform rule and transformed
   return before checking their joins.
3. Combine each complete forcing transform before applying the inverse.
   Check regularity at both inherited real poles and the actual large-real-
   momentum behavior needed for this pairing. Apparent transform singularities
   at zero momentum transfer need their removable values; no numerical zero
   may replace a missing value. Use the fixed explicit profile and saved
   branch/denominator evidence, not a global profile or operator theorem.
4. Apply the saved symmetric PV and coupled signed-current delta prescription
   to these source transforms, with the original measure and order retained.
   Save a full source-specific response representation and its equation/pole
   operands. A smooth-source pairing does not require inventing a pointwise
   kernel value on the diagonal. Any unsupported multiplication or source
   domain stops at its explicit operand, without a new method on the fly.
5. Include only checks that can change this claim: source/sign/unit/measure
   joins, selected independent profile-transform evaluations, the equation
   on this source and one actual responsive thickness-term control. Use a
   small declared set of spectral probes; no adaptive sweep or root queue.
   Do not rerun the earlier density-grade mutation as a substitute for the
   new transform check.

The new transform and outgoing-action rules need applicable independent
method/build review. The review question is concrete: does this native,
fixed-profile thickness forcing lie in the class required by the saved
coupled prescription, and does the proposed representation implement its
action with the correct phases and normalization? It is not a request to
certify all profiles, all coupled operators or complex-frequency causality.
Literal prior review verdicts remain history, not automatic clearance of
these new rules. Optional wording is not grounds for a review loop.

## Bounds, precision and stopping point

Budget one ordinary guarded worker: **900s outer / 840s native**, 2GiB,
zero swap, one CPU, nice15, 32 tasks and one native thread. Persist complete
inputs/returns and reuse saved operations if a future separately approved
continuation becomes necessary. No automatic retry, indefinite exception,
overlap or unguarded fallback. Use a new directory and the existing silent
completion hook for session `01a0e01b-ef84-7192-817f-584cda5d339b`.

The cost of new symbolic contractions is not yet measured. Avoid expanding
large residue constants unnecessarily; keep complete factored operands and
use selected checks. Once numerical response observables are reported, use
the governing practical roughly-1% stability and declared absolute reporting
goals, not arbitrary universal intermediate tolerances. Exact branch/unit/
normalization identities required by the representation still need to hold.

Success would close the localized fixed-input thickness action only. An
honest partial outcome is the complete forcing transform with a precisely
identified pairing issue. Neither outcome closes density/mixed response,
full FORM/export, A11/A12, current normalization, diagonal extension or
retarded equivalence. No general tail/Abel certification campaign, new centre
mechanics, new physical input, complex rows or poles is added to the queue.
