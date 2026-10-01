# Useful work alongside the uniform continuation

2026-10-01. User asked what can be done in parallel to speed completion.
This is source-only author preparation, not a scientific result or independent
review. Two local analysis agents worked independently on the smallest direct
mixed-term test and on its downstream dependencies. Neither imported scientific
code, restored payloads, launched a probe, changed equations or exported files.
The running uniform job and its source pins are untouched.

## Resources and priority

The one-time resource snapshot at 19:28:42 UTC records about 401 MiB used by the
uniform job, an 8 GiB cgroup/native reservation, CPU 15, no swap or memory events,
and about 20.5 GiB host memory available. The authorized aggregate pool is 16 GiB,
so 8 GiB remains available for independent guarded jobs. This could accommodate
two 4 GiB workers on separate CPUs, subject to the guard's live reservation and
host-reserve checks. No computation deadline applies.

The current worker uses a serial schedule of exact symbolic point checks. Extra
cores do not automatically accelerate those calls. Preserve it and its completed
work; do not stop it to retrofit parallel scheduling. The useful independent
work is the mixed-term applicability question, which does not depend on the
new-speed uniform results. Another defect sweep would depend on its answer.

## Small input, specific question

Determine whether the candidate direct `eta*sigma_W` impedance contribution
survives the actual profile/jet relation and common physical dispersion/sheet
constraints. The independent-q expression in the delivered Grok report does not
answer that question. A test must include both assignments of height and slope
to the momentum transfers and distinguish the direct term from the already
retained iterated first-shape product.

The necessary source is small: c1's literal `dtn_kernel` at
`scripts/S11c_c1_exports.py:87` is 17,868 bytes; its flat symbol at line 79 is
2,260 bytes. The face normal/shift/velocity excerpts are already in the upstream
packet. Join them to c2 `kernel_bridge:367`, `fourier_profiles:429`, and
`profile_bindings:850`. The large finite-balance database is unnecessary.

The development profile is `(1+tanh(xi))/2`, with
`sigma_W=eta_bg*W0/L_W`. A Gaussian or periodic witness must be labelled as a
separate generic-operator diagnostic, not substituted for that benchmark.
Enforce one physical bulk dispersion on input, intermediate and output legs;
exclude unproved grazing/pole values; carry the native Fourier convention and
both face orientations. Report any unresolved source join instead of silently
choosing it. A nonzero direct coefficient is not a leakage magnitude.

### Candidate smallest diagnostic, not yet computed

First assess one upper-face scalar bare-impedance action on the actual tanh
profile, at omega=3, c_s0=10, L_W=10, saved edge tangents, profile-direction input momentum 0
and output momentum 1/10. The common outgoing dispersion then fixes input
depth 1/5, output depth sqrt(3)/10 and middle depth
`sqrt(1/25-t^2)`, with the outgoing continuation on the exterior intervals.
These are a source-level diagnostic of the existing benchmark parameters,
not a new end mode or a changed physical pilot.

The saved review expression suggests that setting the profile-direction input
momentum to zero leaves a factor of the height transfer. That might annihilate
the tanh constant-end delta and remove the transfer singularity. A sign/phase
argument on the propagating middle interval might then distinguish a nonzero
direct action from the zero native slot without evaluating a large integral.
**This is an unverified algebraic lead.** Independently rederive the actual
boundary coefficient, including both height/slope assignments, before using it.
Derive the Fourier normalization and distribution products from the native
profile definition; do not hard-code a recalled transform. If those joins or
endpoint integrability fail, report that missing dependency and stop.
Require the original boundary residual, native linear-kernel correspondence,
rigid-translation/zero-jet limits, addressed slope/sheet controls and units.
An upper-face witness would not determine cancellation in the closed two-face
slab response or its physical power.

A weighted tanh kernel value at two nonzero transfers would close only a
pointwise off-shell/profile gap: it cannot by itself exclude convolution or
contact-term cancellations. An exact finite-Fourier periodic witness could
instead test the generic operator with a finite sum, but would use a different
profile and establish no tanh-benchmark correction. Do not silently exchange
these claims. Neither route computes a closed slab response or loss fraction.

The scalar tanh route is the preparation priority. It needs a small immutable
source snapshot and one guarded worker; provisionally reserve 4 GiB on a
different CPU, with the existing pool/no-deadline controls. Runtime is unmeasured.
The exit deliverable is a source-joined direct-term witness, a proved
cancellation on that selected action, or the exact unresolved join. Full c2
repair, producer regeneration, two-dimensional sweeps and power calculations
are outside this diagnostic. It is not ready or launched by this note.

## Reuse map if the term is confirmed

| Object | Required action |
| --- | --- |
| Running selected uniform check | Continue with its existing conditional applicability. Verify the candidate's actual constant-end specialization before promoting that condition. |
| Saved 0.12 benchmark | Preserve as the result of its implemented operator. Its numerical residuals do not bound the effect of a missing term. |
| Old end pencils, modes, currents and face maps | Reuse only where zero-jet specialization proves the corrected contribution vanishes and actual bound operands match. A changed source-file hash alone does not require recomputing them. |
| Nonuniform closed rows/kernel and numerical inputs | Rebuild affected coefficient/row/cell operands, domain/unit/grade joins, quadrature arrays and matrices. Do not inherit old address counts if the correction adds terms. |
| Numerical machinery | Retain basis, quadrature, solver and reporting code where unchanged. Recompute solves/current balances whose matrices change. |

The source chain is c2's zero direct entry at lines 409–417, `build_face` at
555–604, closed slab rows at 640–655, then the d consumer at
`scripts/S11c_d_mixing_scattering_sympy_audit.py:616`. The numerical integration
joins 375 saved addresses and 160 occurrences at
`_measurements/S11c_d_numerical_radiating_integration.py:63–104`; these establish
identity with the saved operator, not completeness of its physics.

This preparation narrows a future guarded test and avoids unnecessary rebuilds.
No new external packet, physics correction, producer run or defect sweep was
started. Applicable method review and exact-packet consent remain separate from
this source analysis. Stop at the selected uniform evidence and upstream
findings under the current scope.
