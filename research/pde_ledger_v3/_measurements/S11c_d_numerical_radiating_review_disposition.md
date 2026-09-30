# Numerical radiating method review — stop before implementation

2026-09-29 local. Both fresh reviewers returned **NEEDS REVISION**. The
submitted method is not independently cleared. No new implementation, scientific
run, payload restoration or second review was performed. This is a source-only
adjudication, not a numerical result or an author clearance.

The numerical route is not shown impossible. Both reviewers consider its bulk
integration approach plausible and propose bounded corrections. Nevertheless,
the unresolved items change boundary selection, independent integration checks
or what a deficit can mean; they are not loader/formatting fixes. The submitted
method's one-round stop rule (lines 246–259) applies. Preserve the exact draft
and return this decision, rather than start another repair/review campaign.

## What finished

The 19-file, 311,706-byte packet
`625116f563970d763870b5d1dfdcd987140ef738db5dae5be6053686ce18bde8`
and its archive, all 18 original source records, prompt and transport hashes
match. Claude finished in 644.796 s, Grok in 569.707 s, in separate fresh sessions
with the same prompt and no peer reports. These are review times, not compute
cost estimates. Both produced usable literal reports. Claude/coordinator stderr
is empty; Grok's 3,598-byte stderr contains CLI configuration/context warnings.
The completion hook finished. The canonical review record preserves both entire
literal reports, their receipts, stderr, authorization and packet inventory.

Neither reviewer executed a derivation, ablation or scientific calculation.
Their endpoint reasoning is prospective evidence, not a completed 80-row census.

## Substantive findings that survive source inspection

1. **End selection needs an explicit recipe at every contrast and frequency.**
   `S11c_d_frequency_source.py:174` specializes LEFT/RIGHT to the nonzero
   contrast and REFERENCE to zero. `S11c_d_frequency_end.py:56` requires the
   resulting coefficients to depend only on frequency. Its `seeds` method
   (lines 119–135) retains five clusters containing seven directions from 18
   candidate records; five directions are outgoing and two incoming. Changing
   the contrast binding alone does not supply new seeds. The proposed method
   needs contrast continuation/selection and an explicit physical-sheet and
   outgoing audit on the frequency path and at its endpoint. Small pencil/wave
   residuals do not by themselves distinguish a wrong sheet. Continuing the
   saved 18 candidates is Claude's suggested bounded check, not an established
   completeness theorem or an already approved correction. This accepts the
   substance of Claude B1/B2 and Grok C1, with the cluster/direction distinction.

2. **The independent integration check must cover the middle momentum leg.**
   The c2 `kernel_bridge` at lines 404–417 contains ordered second scattering
   through a middle radical. Both reviewers support investigating the closed
   endpoint cancellation. Existing `CachedMomentum.adaptive_outer` only varies
   the outer leg and `inner_at_outer` retains fixed inner panels; it does not
   supply an independent check of the singular middle leg. A selected nested
   check using a different rule there is a real missing validation choice.
   The source also supports an assembled uniform non-transverse test against
   the constant-end symbol (Claude B3): the zero-drive transverse control alone
   is blind to this operator. Neither test has been implemented or executed.

3. **Physical leakage and retained-model deficit remain distinct.**
   `UniformSlabCurrent.retained` at engine lines 2605–2609 keeps the rectangular
   grades 1, eta, sigma and eta*sigma. The slab balance residual is projected
   again at line 2853, and c2 line 417 retains the mixed second scattering only.
   Along the proposed joint scaling, omitted pure second grades and the
   leading deficit can have the same contrast order. Thus numerical convergence
   plus approximately quadratic scaling cannot alone rule out retained-order
   non-conservation. This accepts Claude B4 as an unresolved interpretation
   risk, **not proof that an artifact exists or its size is known**. No source-
   only assessment here proves or disproves exact conservation of the finite
   retained model. This finding is grounded in its source; it does not transfer
   the separate retained-response-series omission wholesale onto the old full
   finite solve.

   Claude suggests two sub-threshold, nonzero-contrast solves (32 rather than
   30 total). Those could be diagnostics, but physical face absorption remains
   present below threshold and an error there need not bound an error at omega3.
   They are not a demonstrated cure and are not added to an execution queue.

4. **Finite-end transverse partition needs a cross-current check.**
   `S11c_d_finite_scattering.py:384–391` forms the reported flux from the selected
   open outputs. A finite-end transverse/non-transverse cross term is not
   automatically included in that block or negligible. Its size or vanishing
   must be checked using the actual physical pairing before interpreting the
   selected deficit. This extends the method's cross-term validation; it is
   not a request for a complete outgoing Green operator.

## Findings qualified or rejected

- **Grok B1: reject the blanket undefined-headline conclusion and proposed
  addition of port power to current.** The draft counts selected transverse
  input/output only and explicitly requires a vanishing-loss-drive reduction
  or a stop (lines 169–185). Radiating Fourier legs elsewhere do not establish
  that these selected modal pairs require a divergent depth integral.
  `ClosedCurrentPairing.construct`, engine lines 3374–3396, has
  `total_current = slab_current + depth_integral*bulk_current`; interface power
  enters the balance separately. `ModalCurrentSubspaces` lines 4136–4138 also
  gives current a length factor relative to power. Adding the raw interface/
  port-power matrix to slab current is not the source's flux definition.
  Accept the obligation to verify actual pair domains and reject undefined
  infinite-depth products; do not implement Grok's replacement formula.
- **Grok B2: accept a leg-by-leg source/scale join, not a demonstrated wrong
  threshold.** c2 lines 441–446 use the same wave relation on each ordered leg,
  while the compact chart verifies three separate radical arguments. Actual
  bound radicands and physical scale maps must agree with the proposed split.
  No source inspection here has found unequal endpoint locations. Profile
  coordinates are not radical-bearing Fourier legs.
- **Grok B3: accept explicit memory-kernel notation and endpoint evidence.**
  c1 `lambda_kernels` uses `Lambda_A_0/(1-i*omega*tau_A)`, and `port_matrix`
  consumes that kernel. The draft's `Lambda_A` is underspecified; it must not
  become the static JSON amplitude. This is a specification clarification,
  not evidence that an executed worker froze physical damping. Both reviews'
  cancellation arguments still require actual source-bound checks.
- **Grok C2: reject bulk availability as proof of uniform transverse loss.**
  The draft requires the incoming transverse state to have vanishing loss-side
  drives at each target. If that premise fails, it stops. An open acoustic band
  alone does not invalidate a zero-contrast control under that premise. Keep
  raw deficit, matched-uniform difference and separate numerical sensitivities
  visible; the uniform result is neither an absolute error bound nor a known
  subtractable absorption floor. The draft already carries those limitations.
- Optional loss-matrix eigenvalues, extra path sampling and regulator/bandwidth
  reporting are not converted into a new prerequisite campaign. No optional
  review cycle follows this disposition.

## Decision and retained result

**NO GO for implementing/running the submitted version.** This is not a claim
that radiating numerics require the full exact machinery, or that two authoring
days have elapsed. Substantial budget remains; the stop is the agreed method-
blocker rule. A revised method would need a new explicit decision and applicable
independent review; none is queued. Exact RIGHT stays parked. No new science
gate exists and no calculation is running from this proposal.

The honest physics status is unchanged: the saved omega1 finite calculation
resolves no deficit at its declared 1e-6 current resolution, without supplying
a physical upper bound. Selected REFERENCE/LEFT premise evidence survives;
RIGHT is unresolved from the summary guard. The proposed above-threshold
frequencies remain untested. Calibration of analog-light frequency, total loss,
full Green/FORM and A11/A12 remain open. Do not spend the remaining budget merely
to obtain a finite number whose physical interpretation these gaps leave open.
