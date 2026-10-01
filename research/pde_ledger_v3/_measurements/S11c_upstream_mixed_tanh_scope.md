# Focused direct height–slope diagnostic

2026-10-01. Prepared after the user's “great. let's keep going,” followed by
“k. approved. keep going.” This authorizes the focused next test described in
`S11c_upstream_parallel_next_work.md`. It does not authorize a defect sweep or
change the production equations. The selected uniform result is preserved.

The new instrument asks whether the **bare upper-face impedance coefficient at
wave order one and background grade eta*sigma_W** survives the actual tanh
profile, its derivative identity, and the common outgoing bulk dispersion.
The prescribed normal-velocity input has unit amplitude; this is not an
incident transverse slab mode. A nonzero result would identify a native
operator-completeness issue, not establish excitation by light or a loss value.

## Concrete inputs and calculation

Use LAB_HELD/RHO4_CONSTANT, zero drain, omega=3, rho_m=1/10, c_s0=10,
W0=1, L_W=10 and saved conserved edge momenta (1/5,1/10). Keep the original
development-input file unchanged; omega=3 is the explicit diagnostic override.
Choose profile-direction input momentum zero and output momentum 1/10.
This is a diagnostic source action, not a new end mode. All three depth legs
obey q(p)^2=1/25-p^2 with positive real outgoing q in the propagating interval
and positive imaginary q outside. Input and output are nongrazing; the two
middle-leg branch points must be included rather than excluded from the action.

1. Restore only c1's literal flat/kernel records through their serialized
   constructor expressions. Extract and execute the small native shape-source
   and first-kernel functions; compare actual operands. Import no producer and
   replay no completed uniform, path, matrix or integral calculation.
2. Construct the rectangular h/s boundary Taylor coefficient from the original
   graph normal and shifted pressure. Keep both orderings, arbitrary input
   momentum in the derivation, and the input/output/intermediate depth legs.
   Verify the original boundary residuals, native first-height/first-slope
   coefficients, amplitude homogeneity and rigid-height/zero-jet limits.
3. Bind the actual `(1+tanh(y/L_W))/2` profile using c2's own local-jet and
   Fourier-definition functions. Derive the derivative transform by the
   logistic substitution and Euler beta/reflection identity, including its
   removable zero value. Keep the constant delta and principal-value height
   transform explicit. A nonzero total transfer alone does not remove these.
4. Form the mixed convolution with both native half-height/half-slope factors.
   Multiplication by the height transfer must remove the constant delta and
   principal-value singularity exactly. Show the fully symmetrized transfer
   expression with its factor 1/2; do not double the full convolution.
5. Save physical-sheet expressions, endpoint numerators, tail behavior,
   sign/phase ingredients, units and responsive native slope/sheet controls.
   Assess whether these imply noncancellation without computing a large
   integral. If any source, distribution or endpoint join is missing, report
   it as unresolved. No numerical smallness or unknown-to-zero conversion.

Pure eta^2 and sigma_W^2 terms are outside the existing rectangular retained
grade. The output must not be called the full second-order shape correction
after sigma_W=eta*W0/L_W is inserted. The bare coefficient is also distinct
from the closed permeable response's already-retained iterated first-shape
product, and from the inverse map to reference pressure.

## Controls, persistence and scope

Each identity saves both operands and the computed residual before guarding.
The worker writes immutable JSON with exact expression representations, stage
inputs/returns, source posthashes, strict stderr and final checks. It restores
no pickle. The first-harmonic convention, native Fourier normalization, two
conserved edge deltas and four-spatial-dimensional density units remain explicit.
The source-native direct three-leg entry is inspected as an AST literal; the
c2 producer is not executed or changed. No closed slab response is constructed.

One 4 GiB worker in the existing 16 GiB pool, one assigned CPU/thread, zero
swap, 32 tasks, 4 GiB host reserve, shared guard and existing normalization
supervisor. No wall/native/CPU/inactivity deadline; desktop-managed priority,
no scheduler changes. Use a fresh pinned gate and the existing local completion
hook before execution. Stop after this diagnostic; no automatic retry, loss
pilot, new producer or defect sweep.

## Review and authorization

This is a physics-bearing instrument, not a formatting fix. CLAUDE.md E1/G1
requires a Codex-written diagnostic to receive fresh Claude and Grok review
before its output is trusted. This packet is one prospective build/method
review of the complete small instrument and its question. Reviewers receive
native sources and physical inputs, no earlier peer reports, agent derivations
or external commentary. Review clearance is not runtime acceptance. Exact
external-packet consent is separate from the user's already-given next-test
authorization; no repetitive science permission is needed after substantive
clearance and actual readiness.

The possible exit records are a source-joined direct-term witness, a proved
cancellation on the selected action, or the exact unresolved dependency.
Any production correction and its downstream rebuild remain separate work.
