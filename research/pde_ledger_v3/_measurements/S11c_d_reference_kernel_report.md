# Accepted bounded reference inverse ingredient

The approved first stage is complete for `LAB_HELD/RHO4_CONSTANT`.
The [acceptance checkpoint](S11c_d_reference_kernel_checkpoint.json) records
the exact saved artifact routes and hashes, input joins, resource records,
review disposition and validation results. Acceptance covers the algebraic
ingredient below; it does not close the full FORM task.

Available in
`_scratch/s11c/s11c-d-reference-kernel-20260926/production/complete`:

- The full saved five-by-five strong reference symbol's explicit scalar
  determinant, transposed cofactors and inverse entries, conditional on a
  nonzero determinant. Frequency and tangential momentum remain symbolic.
- The regular Fourier spectral density with the saved normalization and
  phase, physical row/field/inverse/kernel units, and branch/profile context.
- Fifty zero generic adjugate-identity certificates with exact source-entry
  bindings, physical left/right residual expressions, three fixed numerical
  probe comparisons and one-sided source mutations.
- All 115 operation input/result receipts and 497 artifact files
  (3,551,424 bytes). Every constructor artifact and all 16 consumed source
  routes remain unchanged.

The checkpoint's `primaryArtifacts` identifies the determinant, inverse,
density, certificate, context and units without needing to scan the payloads.
The `artifactSha256` inventory covers the complete result directory.

## Validation and cost

The user explicitly approved the corrected saved-output validation after the
first validator stopped on expression layout. SymPy reconstructs the saved
Add/Mul expressions canonically; the original comparisons incorrectly required
the unevaluated layout to survive serialization. The separate corrected
validator passed all **17 structural checks** and all **three numerical probe
checks**, without reconstructing the inverse or repeating LU solves. The
original failed run, validator and diagnosis remain preserved.

At normal momenta `0, 1/10, 1/5` in the unchanged physical input frame, both
inverse residual norms are below `9e-80`; differences from the independently
constructed LU inverse are below `8e-81`. Condition estimates range from
10.48 to 11.90. The source-entry sign mutations produce response norms from
3.22 to 3.38. The corrected validator reproduces these metrics from the saved
probe matrices. Stored-residual roundtrip differences are below `5e-80`.
These are reference-frame algebraic checks, not physical error estimates or
radiation-domain coverage.

| Run | Measured time | Peak whole-job memory | Outcome |
|---|---:|---:|---|
| Constructor | 19.314 s worker; 19.583 s guard | 109,998,080 bytes (104.9 MiB) | Completed |
| Original saved-output validator | 5.036 s supervisor | 82,608,128 bytes (78.8 MiB) | Preserved layout-check failure |
| Explicitly approved corrected validator | 4.394 s supervisor; 4.522 s guard | 95,674,368 bytes (91.2 MiB) | Passed |

All runs used the required 2 GiB, zero-swap, single-CPU, nice-15, 32-task,
single-native-thread containment. No swap or memory-limit event occurred.
Successful runs have empty scientific stderr, matching stdout/checks and
verified source posthashes. Completion hooks returned to the existing session;
there is no active job or automatic follow-on. These measurements do not
estimate the cost of the remaining FORM construction.

Claude's build verdict remains CLEAR FOR THIS BOUNDED INGREDIENT. Grok's
remains NEEDS REVISION; its sole scalar-index request was adopted through the
documented local disposition. No new independent CLEAR or independent
physical validation is claimed.

## Remaining boundary

The actual result status remains
`REGULAR_SPECTRAL_DENSITY_BUILT_OUTGOING_PRESCRIPTION_PENDING`.
The pointwise inverse does not prescribe the treatment of real normal-momentum
singularities of the coupled operator. No completed outgoing Green operator,
two-asymptote response, final FORM root, A11 witness or A12 radiating-domain
coverage is available from this stage. The existing physical point has no
real-momentum radiating bulk support. Centre mechanics and the deferred pole
queue were not reopened.

The next method decision is how to construct the coupled outgoing prescription
from the retained branch and reference-mode/current evidence, then handle the
unequal profile end limits and actual forcing/extraction maps. The constructor
proposal requires a concrete scope and cost decision before expanding into
that work. This report and checkpoint close the approved ingredient stage;
they authorize no new calculation, external review or retry. Lean files,
shared guard changes, accepted older science and the protected builder-report
suffix remain untouched.
