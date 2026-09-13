# S11c-d SymPy builder checkpoint

The approved acoustic/face-lift split is implemented. The focused checkpoint
has **802 literal zero scalar residuals in 34 families**, with complete
metadata and no unresolved dimension constraints. No new upstream repair is
indicated by these results. The completed mechanical-load repair remains
committed as `c643112a`; this continuation is uncommitted.

The native `ClosedCurrentPairing` builder verifies the local acoustic identity
before closing the face amplitudes, reconstructs the full face rows and
bulk/port matrices, and composes their finite-depth balance with the slab
balance. The equal-depth branch has a separate wave reduction. Clearing the
original inverse carriers resolves the unequal-depth composition without
expanding the full response algebra; its denominator domain remains explicit.
The original import wiring, Fourier reduction and current classes are intact.

Regular-sheet normal/frequency derivative operands retain the radical chain
rule, row-power maps, interface exchange and upper-boundary current. Their
local wave and source product-rule checks are complete. Full-subspace mode
contraction, the modal current/frequency-pairing relation and flux normalization
remain next in the [continuation plan](S11c_d_modal_current_plan.md).

The [current report](S11c_d_modal_current_report.md) and
[checkpoint inventory](S11c_d_modal_current_checkpoint.json) specify the
symbolic and reference-material-slice coverage. The new focused transcript
under `scripts/out/` is 1,386,029 bytes, with 117 object/metadata pairs and 232
indexed source assignments. Heavy objects and report hashes use fingerprints
and DAG digests. The full calculation took 221.17 seconds; source-checked
publication after correcting zero-object unit tags took 65.26 seconds.

The four-case production transcript remains the committed repair checkpoint;
no new full run or incomplete export was produced. All ten broad TODOs remain,
including mode/flux normalization, scattering, profile-frequency bound poles,
survival, bookkeeping, weak coefficients, controls and export. Threshold,
defective, sheet, denominator and convergence limitations and supplied
upstream debts remain explicit. Constant-end momentum poles remain distinct
from §3b's profile-dependent frequency poles. No authority change, S10/Lean
edit, review leg, comparator, Wolfram, downstream step, commit or push occurred.

## Retained user-approved solver/export contract


1. Preserve `EdgeReduction`, the positional three-parent fold, and the exact
   direct-lookup manifest. All numerical assembly must consume the computed
   reduced rows, including the full nonlocal terms and full coupling vertex.
2. Accept explicit, independent dimensionless profile functions w(xi), m(xi),
   their derivatives, asymptotic limits and tail information. A selected smooth
   step with an independently adjustable localized modulus bump is a numerical
   instance; it does not replace the interface class. Store the profile formula
   and digest in every case record. The current preflight input is recorded in `S11c_d_channel_preflight_input.json`.
3. Require a complete parameter map in a declared L/T/M unit frame, real
   continuum frequency and tangential momentum, small contrast, and
   sigma_W = eta_bg W_0/L_W for evaluations on the physical homotopy. Retain
   independent eta/sigma grades in the symbolic calculation. Test actual
   reference/end channel availability before attempting flux normalization.
   Do not manufacture an incident channel by assigning a sector label.
4. Compute both-end modes, left/right normalization, the S11b-derived current,
   and the variable-profile matching problem. Re-expand the continuum response
   to the retained rectangle; do not present a finite-contrast numerical
   solution as a higher-order continuum prediction. Retain evanescent matching
   modes, channel degeneracies and domain failures explicitly.
5. Numerical pole searches have an explicit profile, parameter map, sheet,
   bounded search region and isolating contours. Evaluate the retained operator
   without the continuum re-expansion. Record boundary/quadrature resolution,
   domain size, precision, root residuals, contour-count evidence, and changes
   under refinement. A bounded search does not establish a global pole set.
   An unsuccessful or inconclusive search is unresolved, not an empty pole set.
   Compute residues/projectors and sheet/decay/width/closure tests only for
   actually resolved candidates; emit spectral overlap, not capture probability.
6. Separate transparent symbolic expressions from evaluated numerical records.
   Symbolic operator/continuum/weak-coefficient exports remain differentiable
   SymPy expressions, compacted with algebraic equivalence checks. Numerical
   mode and pole datasets retain their input bindings, domain, convergence
   evidence, dimensions and truncated-model status. They are not stand-ins for
   a generic symbolic profile-dependent root function. This is the export
   distinction motivating the user-approved contract; the downstream consumer
   will need to bind the appropriate representation explicitly.
7. Fingerprints summarize already constructed/evaluated objects. They do not
   evaluate nonlocal integrals or replace a spectral solve. Algebraic PIT and
   physical numerical evaluation are separate records. All completed roots use
   fresh lowerCamel write-keys and the existing bind-closure/minimal-delta guards.
8. Finish the one-case path, then implement controls and bookkeeping, then run
   all four cases once and write the complete export. No review legs,
   comparator, Wolfram engine, downstream stage, or commit belongs to this lane.
