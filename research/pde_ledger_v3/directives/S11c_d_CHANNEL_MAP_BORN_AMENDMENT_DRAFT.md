# Draft: bounded channel map and conditional leading-order leakage

**DRAFT — pending independent review, not execution authority.** 2026-09-28. This proposes a
cheaper scientific target following the user's quoted discussion. It changes
neither accepted source bytes nor historical results. The immediate task is
preparation and the two amendment reviews requested in the user's clarification;
stop afterward for the user's go/no-go. No physics run, implementation
repair loop or response construction is authorized by this draft.

## 1. Proposed scope

For LAB_HELD/RHO4_CONSTANT, answer where propagating thickness-like end
conversion or bulk-depth radiation is available, then estimate the leading
weak conversion only where the existing §2 reduction is justified. Retain
the original material/closure values, density convention and full coupling;
do not remove physical damping to obtain open channels.

This is a scoped numerical/leading-order result. It does not complete the
four-case general FORM/export, strong-edge conversion, pole branch or all
A11/A12 obligations. Report end-channel conversion and bulk-depth escape
separately. Supported leading-order coefficients may be handed off with their
domain and normalization; missing bulk terms remain unavailable.

## 2. Channel map

Prefer a bounded (omega,k_parallel) map, restricted to points with an actual
open incoming transverse branch of the full end pencil and nonzero physical
incident current. The saved original uniform-end strong symbols retain
frequency and both tangential components; the smaller frequency-live pencils
have tangential momentum frozen. First assess whether using the former needs
only a small binding step. If 2D binding/coverage is not cheap, use the user's
explicit fallback: vary omega at saved tangential momenta (1/5,1/10) and limit
every conclusion to that slice. No full source-producer replay is needed just
to recover symbols already saved.

Here k_parallel has two components. A scalar-magnitude map must establish
azimuthal equivalence from the actual uniform symbols or fix the saved 2:1
tangential direction and label that restriction. Do not assume rotational
symmetry. Fixed-angle curves may be overlaid later if the map supports them.
The later execution proposal must pin a finite frequency/momentum domain and
coverage/work budget before launch. A live parameter does not make a 2D grid
free; no unbounded search or claim that it costs the same as 1D is authorized.

Use the saved source operators/pencils, both ends, physical field units and
actual branch/current classifiers. Roots of a cleared polynomial are only
candidates: retain denominator artifacts, opposite sheets, full nullspaces,
thresholds, ambiguous rank and leaky complex modes. Do not classify thickness
by the norm of selected field coordinates or turn complex roots into open
real-normal-momentum channels. Preserve physical damping and distinguish it
from numerical regulators. Use the full coupled pencil even when a baseline
zero permits subsequent sector simplification.

Bulk availability uses the separate outgoing acoustic depth root and allowed
output normal-momentum continuum at the conserved frequency and edge-parallel
momentum. Edge-normal momentum transfer is allowed. A channel map is not a
leakage magnitude or proof that coupling to an available channel is nonzero.
A sampled map supports sampled claims; interval-wide absence requires actual
coverage evidence. Regions without an admissible incident branch or with
unresolved classification are labelled unavailable/unresolved, not lossless.

## 3. Conditional Born calculation

The [baseline assessment](../_measurements/S11c_d_born_shortcut_baseline_assessment.md)
records computed zero K0/Kminus/Kplus in both directions for the saved case
and endpoints. Consume those actual operands with source/unit joins. Also
establish §2's regular-domain and reduced/full comparison premises at the
claimed settings. Exclude thresholds, gap closures, unresolved resonance
enhancement and accumulated repeated conversion. Zero baselines alone do not
clear the shortcut.

Review whether the zeros are symbolic in omega/k_parallel or only pointwise.
The later saved all-case inventory also reports zero coupling blocks in all
twelve case/end backgrounds; assess its scope without assuming identical
end currents or adding three more response calculations to this proposal.

Construct the first-order amplitude from the complete own reduced kernel,
including tilt, modulus-gradient/advection and supported zero-jet/first-jet
terms. Preserve independent eta and sigma_W until the declared physical
homotopy is applied. Derive the actual transfer/profile dependence and keep
phase-matched/removable values and step-tail distributions where they matter.
Do not insert a preselected scalar “coupling squared” formula or assume
universal exponential suppression/nonzero sharp-step leakage. Width changes
must remain within the retained smallness/regularity domain.

Normalize with the actual incident and outgoing current forms. Bulk
availability uses the radiation-selected acoustic depth root. Bulk power
additionally requires the actual outgoing exterior solution driven by the
same slab face state on both disconnected half-spaces, its signed far-field
flux and measure, and the incident-current denominator. The restricted c1
far-field check alone does not establish this map; absent the supported map,
report availability only and leave bulk power pending. Neither a slab-current
deficit nor a decaying-mode norm is bulk escape. An unavailable radiating
normalization stops that observable; it does not trigger the full construction
automatically. Report a flux ratio or power coefficient, not a temporal decay
rate without its additional map.

First compare with matching saved omega=1 coefficient/field operands under
the same units, grades, profiles, regulator and boundary convention. Reuse
the solve; do not rerun it. Emit the actual reduced/full residual and a
responsive physical-term omission control. Use practical approximately 1%
checks where the reference resolves the observable and a declared absolute
tolerance near zero. The saved approximately 1e-4 thickness amplitudes and
1e-4 reported resolution may make that comparison inconclusive. An inconclusive
anchor stops validation rather than forcing a precision campaign. A match
there validates only the compared closed-domain construction; open currents
still require their own domain-appropriate evidence.

A contrast-scaling check is an alternative to a percent-level amplitude match:
test the leading amplitude and full-minus-leading remainder along the declared
eta/sigma homotopy where the relevant nonzero signal is resolved. Reuse saved
evaluations first. Existing coefficient-polynomial remainder checks are not
automatically a Born/full comparison. Any new contrast points need a later
bounded plan; do not repeat completed points. If no useful anchor exists,
label the result “unanchored leading-order estimate” and report baseline-zero
evidence, regular-domain checks, reduced/full comparison and open-domain flux
validation separately, identifying what is established and what remains
unestablished. This label does not discharge §2's comparison requirement or
A11/A12. The amendment reviewers must assess this explicit scope exception;
it is not an author-issued waiver or a reason to force a precision campaign.

## 4. Budget, review and stop

After a scope go, channel-map and Born implementations each need the applicable
independent two-leg review before their results are trusted. Their execution
plans must give concrete budgets and stop conditions, reuse completed inputs
and preserve failed operands. Ordinary shared guard/supervisor limits apply:
900 s outer/840 s native, 2 GiB, zero swap, one CPU, nice15, 32 tasks, one native
thread, one job at a time, hook first, no automatic retry or replay. Runtime
scratch stays uncommitted. No old indefinite-time exception carries over.

Review this draft specifically for: the sufficiency and scope of the saved
baseline zeros; the path/conservation and channel definitions; whether a
useful reduced/full comparison is possible at saved precision; and the minimum
supported open-current/bulk map. Identify a concrete blocker if the cheap
route cannot answer its restricted question. Optional wording does not create
a review loop. Report both literal verdicts and the resulting go/no-go to the
user, then stop before implementation. A failed shortcut does not authorize
the full radiating method or a PDE simulation as an automatic fallback.
