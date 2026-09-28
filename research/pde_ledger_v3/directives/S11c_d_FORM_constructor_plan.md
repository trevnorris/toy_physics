# S11c-d FORM constructor: response representation and first method decision

Status: **proposed construction, not implementation or physics clearance**.
Prepared after resuming session `01a08c53-295c-70f0-baaf-353db1745867` in a
replacement conversation on 2026-09-26. No scientific job or external review
has been launched. The accepted
[amendment](S11c_d_SCATTERING_FORM_AMENDMENT.md),
[current build instructions](S11c_d_FORM_build_directive.md), and completed
[centre claim disposition](../_measurements/S11c_d_centre_claim_disposition.md)
govern this proposal. This document does not reopen centre mechanics.

## Decision in concrete terms

Develop the missing response from the **full saved reference symbol's outgoing
Fourier Green kernel**, retaining its continuous contribution, then construct
the two-asymptote independent-grade response from the saved reduced operators.
Begin with one bounded reference-kernel construction in the existing
`LAB_HELD/RHO4_CONSTANT` case. Do not commission the whole radiation solver or
a new physical-frequency response as part of this first stage.

This is the new boundary/continuum-method decision anticipated in amendment
§2. The previous source report proposed this family of methods but did not
authorize or construct it. Here the first stage, reusable inputs, downstream
dependencies and stopping rule are specified for a decision. The scoped
implementation/instrument gate follows route approval; no new whole-amendment
review is proposed.

| Alternative | What it would do | Disposition |
|---|---|---|
| Full-reference outgoing spectral kernel | Construct the actual inverse symbol and its outgoing integral, including branch/continuous content; use it in retained-grade response construction. | Recommended first route, with the bounded stage below. It directly addresses the missing general response ingredient without recomputing accepted finite responses. |
| New exterior-domain numerical method | Replace the current end closure with a representation that resolves exterior radiation and its coupling to the slab. | Larger new method and numerical acceptance burden. No mesh, absorbing-layer or continuum-discretization campaign is proposed now. |
| Existing finite numerical annex only | Publish the already evaluated finite solutions and their current bookkeeping on their recorded domains. | Useful existing evidence, but not completion of A9/A11/A12. Selecting this as the final deliverable requires an explicit scope reduction. |

## What is available, and why wrapping the old solver is insufficient

The [source receipt](../_measurements/S11c_d_FORM_constructor_source_receipt.json)
pins the selected definitions and accepted packet addresses without importing
science modules or restoring pickles.

| Consumer | Reuse and actual boundary |
|---|---|
| `uniform_source.construct` | Four saved `uniform-source.pickle` packets contain `records[REFERENCE/LEFT/RIGHT]`, full `strong`/`weak` symbols, `coupling`, units, and Fourier mass. Consume them; do not call this constructor again. The reference is not a manually diagonalized bare-sector model. |
| `ReducedPencil` and saved action/grade packets | Five fields `(u1,u2,u3,theta,eW)`, full local derivatives and nonlocal kernels. Reuse the native reduced action, Fourier measures, profile tails, grade coefficients and units. No `EdgeReduction` replay. |
| `continuum_boundary.construct_end` | The existing selection requires sheet membership **and bulk-decay certification**. Its current-pair contraction requires both bulk momenta to have positive imaginary parts. Five outgoing and two incoming trace columns per end are observed finite-case choices, not general completeness statements. |
| `continuum_response.systems/solve` | The finite problem replaces endpoint rows with derivative-minus-trace equations and retains incident insertions; its recursion includes the mixed operator forcing. Reuse these dependency identities and saved evaluations, not their fixed arrays as a generic inverse. |
| `continuum_response.channels` | Extraction subtracts incident traces, applies outgoing-coordinate inverses and origin phases, and varies incoming/outgoing current normalization. These dependencies must enter the FORM. |
| `continuum_currents` | Current products retain interference entries, slab/bulk parts, and a varying incident denominator. End-normal current containing a depth-integrated bulk contribution is not bulk-depth escape. |
| c1/c2 face and radiation maps | Already indexed saved fields, closure maps, reference traces and normal jets feed the separately typed bulk observable. The five-field `DELTA_W` claim boundary and inherited c1/c2 debts remain attached. |

One useful restriction follows directly from the inspected rest-acoustic wave
equation and the saved development-input numbers. For real normal momentum
`k_n`, that equation gives

```text
q² = (omega/c_s0)² - k_parallel,1² - k_parallel,2² - k_n².
```

The saved point has `omega=1`, `c_s0=10`, and tangential components `1/5,1/10`.
Thus this rest-acoustic expression is `q² = -1/25 - k_n²`: its real-momentum
bulk spectrum is evanescent. This is a hand substitution into the inspected
source equation, not a new calculated response or a general dispersion result.
Tangential momentum is conserved by the specified profile ansatz. Consequently
testing a new kernel only at this point cannot provide A12's nonempty
radiating-support evidence. It also does not decide whether a different
admissible input supplies A11's thickness-like end channel. No input parameter
is changed by this proposal.

## Constructed representation, not a named inverse

Write `P0(k; omega,k_parallel)` for the actual saved five-field reference
symbol, with its physical row/field units and its source-derived acoustic
branch. The candidate reference response is an outgoing inverse Fourier
integral of `adj(P0)/det(P0)` with the exact inherited Fourier measure and
phase convention. This is a design equation, not an already emitted result.

The constructor must instantiate the actual scalar entries and calculate the
cofactor numerators and determinant from those entries. Preserve shared
subexpressions and unevaluated, explicit finite sums where expansion would
inflate them. A downstream decoder must recover scalar algebra and the actual
integral integrand; a string named `P0Inverse`, a fingerprint, or an unbound
matrix-inverse node does not finish this construction. Do not invert the
gauge-redundant weak-sector matrix instead of the physical strong symbol.

The integration prescription is part of the output, not metadata that can be
left unspecified. It must identify the inherited time/Fourier sign, outgoing
bulk branch, boundary value or contour, real-axis singularity treatment,
measure and regularity domain. Keep a profile Abel regulator distinct from any
limiting-absorption prescription; substituting one for the other is not allowed.
Do not choose a square-root sign independently at each integration point.
Branch-cut/continuous contributions cannot be replaced by the existing finite
mode list. Any needed analytic-continuation or stability premise must be
identified rather than inferred from a successful matrix inversion.

No frequency-pole search, contour census of profile poles, global domain
certificate or inverse-norm theorem is added. If the actual outgoing
prescription cannot be constructed under the supplied assumptions, preserve
the inverse-symbol work and stop with that precise dependency. An algebraic
inverse alone is then a completed ingredient, not a completed outgoing kernel.

## How this ingredient enters the retained FORM

The response constructor will preserve the background-grade rectangle
`(eta,sigma) = (0,0),(1,0),(0,1),(1,1)` and the native epsilon bookkeeping.
For an actual constructed outgoing boundary-value problem `A U = F`,
coefficient convolution gives the following dependency pattern:

```text
A00 U00 = F00
A00 U10 = F10 - A10 U00
A00 U01 = F01 - A01 U00
A00 U11 = F11 - A10 U01 - A01 U10 - A11 U00.
```

Here `A` includes the full operator and the actual outgoing/incident conditions;
it is not just the local slab matrix. These equations specify required
dependencies, not permission to export an unconstructed action. A Green-kernel
implementation must construct their forcing and extraction maps, including
the variations of end modes, phases, current forms and incident normalization.
It may not discard reference off-diagonal terms or replace the response with
a bare `K1` Born matrix element without the governing computed premises.

Because the thickness profile has unequal end limits, `L-L0` is not simply a
compactly supported perturbation. The subsequent two-end construction must
retain the actual left/right asymptotic fields and their grade variations,
and the source's half-line/distributional tail terms. Naively integrating a
nondecaying forcing against the reference kernel is not the completed method.
An auxiliary lifting or tail subtraction, if used, must carry the full
nonlocal operator and its reconstruction/cancellation checks; it cannot erase
cross-interface terms. This is a downstream method obligation, not work hidden
inside the first-stage budget.

For outgoing amplitude `S_g` and current `J_g`, form all coefficient products
`S_a^* J_b S_c` before applying the requested output-grade selection. Then
derive each incident-state fraction with the actual incident denominator,
retaining its variations. Only afterward impose the declared physical
eta/sigma homotopy. The full total, baseline/interference contributions and
induced-field quadratic remain separate outputs. A retained-model quadratic
does not establish missing parent-theory second-shape power.

| Prospective export root | Dependency on the new construction |
|---|---|
| `s11cdRetainedScatteringForm` | Completed two-end outgoing response, forcing/extraction maps, all four cases and declared domains. The first-stage reference kernel alone does not complete it. |
| `s11cdEndChannelConversionForm` | Same response, actual thickness-like outgoing subspace, full current pairing and transverse incident denominator; A11 remains open until its nonempty-domain check. |
| `s11cdTransverseSurvivalForm` | Same response with the transverse observable; preserve its restricted interpretation. |
| `s11cdSupportedBulkEscapeForm` | Solved-state-to-face maps, outgoing exterior fields and depth flux/measure, including supported scattered–scattered terms; A12 remains open. |
| `s11cdRetainedWeakCoefficients` | Actual retained response and observable convolutions, with units, independent grades and capability distinctions. |
| `s11cdFormCapabilities` | Completed versus unavailable case/domain/grade claims and inherited debts; never a replacement for missing expressions. |

## Bounded first stage proposed for approval

**Deliverable:** a source-bound reference-kernel artifact for the existing
`LAB_HELD/RHO4_CONSTANT` case: full saved reference symbol, explicit inverse
entries, source-derived outgoing integral prescription where established,
units and bind dependencies, actual residual/control operands, and an honest
domain/obstruction record. It is a reusable construction ingredient, not a
new scattering result, radiating witness, published FORM root or clearance.

**Reuse:** select the baseline packet through the accepted
`S11c_d_remaining_case_end_sources_checkpoint.json`; it is byte-identical to
the original baseline uniform packet. Consume the saved reduction/branch/unit
context and existing physical-input JSON. Keep frequency and tangential
momentum symbolic in the constructed expressions. Only the already approved
physical point may be used for a bounded algebraic binding check in this stage.
Normal momentum is the kernel's integration variable, not a newly selected
incident channel. No old closure, mode, current, Taylor or response is rebuilt.

**Checks to specify in the implementation gate:** left and right inverse
residuals with physical units; an independent direct scalar/linear solve at a
few regular probe momenta in the existing input frame, away from saved modal
singularities; source-derived outgoing branch/phase evidence; a one-sided
mutation of an actual nonzero symbol entry with the constructed inverse held
fixed. Compute residuals and control responses rather than typing expected
values. Preserve failures and non-applicable controls. Algebraic inversion
checks do not establish the integral's radiation/domain coverage.

**Cost and stop:** at most one new guarded symbolic/algebraic worker, 900 s,
2 GiB, zero swap, one CPU, nice 15, 32 tasks, one native thread, after the
scoped implementation/instrument gate. This is an execution cap, not a runtime
prediction. Persist each completed numerator/denominator/prescription/check
before later guards. No automatic retry or move to a larger memory profile.
Authoring/review time and the cost of the eventual full four-case FORM remain
unmeasured. Do not infer a total-project estimate from the small saved packet.

If the stage finishes, report exactly which ingredient is available and the
measured cost before expanding to the two-asymptote construction. If it reaches
a missing branch premise, expression-size obstruction or the resource limit,
report that obstruction and saved partial artifacts. A11/A12 witnesses and
new physical bindings remain separately planned duties, not automatic follow-on
jobs. This proposal does not request a new centre proof, old pole queue,
external review submission, or unrestricted radiation campaign.
