# S11c-d SymPy builder checkpoint

The LAB_HELD/RHO4_CONSTANT right-end source calculation now runs on this box
with all material parameters symbolic. Native `ClosedAcousticEnergy` (coefficient collection at
2082; mass/mechanical comparisons at 2228–2229) uses exact structural
collection before material cancellation. The resumable instrument preserves
raw operands and denominator exclusions; proof arithmetic combines fractions
before sparse cross-multiplication. **All 127 exact reconstruction residuals
are zero.** No physical input, source formula or extra truncation was changed.

The computed RIGHT source joins have **10 nonzero scalar entries**: five mass
and five mechanical, each at `(epsilon,eta,sigma)=(0,1,0)`, lambda order one.
The other 23 top-level source residual scalars are zero. This is a retained
source-consistency discrepancy whose origin remains undetermined. It is no
longer an uncomputed dependency deferred for machine capacity. See the
[source report](S11c_d_end_current_source_report.md) and
[runtime plan](S11c_d_right_current_runtime_plan.md).

The final source/proof replay and emission took 345.70 seconds at 201436 KiB
peak RSS and produced 9520298 bytes. This excludes prior construction/proof
attempts. Exit 1 followed complete emission through the computed-residual
guard. Both transcripts are validated and atomically published. The fresh LEFT run
took 123.03 seconds at 201636 KiB: all 33 source residual scalars and 82 exact
arithmetic identities are zero. All 85 source objects agree with the preserved
baseline (83 structurally; the complete 52-entry parameter map by key; two
face-load expression comparisons with zero residual). See the
[RIGHT inventory](S11c_d_end_current_source_right_checkpoint.json),
[LEFT regression inventory](S11c_d_end_current_source_left_runtime_checkpoint.json)
and [runtime inventory](S11c_d_right_current_runtime_checkpoint.json).

The independent **both-end frequency construction is complete** for the
LAB_HELD/RHO4_CONSTANT development binding: 36 isolated-root/lift candidates,
44 full basis directions, 90 exact-zero reconstruction scalars and 3920
numerical algebraic scalars with maximum norm 2.535e-13. The 1800 finite-frequency
checks and complete subspace/rank/domain records remain unchanged. See the
[frequency report](S11c_d_end_frequency_report.md),
[LEFT inventory](S11c_d_end_frequency_left_checkpoint.json) and
[RIGHT inventory](S11c_d_end_frequency_right_checkpoint.json). These algebraic
frequency projectors do not constitute section 3b bound-pole/Riesz data.

The runtime work is complete. Stop and assess the
first-contrast source discrepancy before changing any physical formula. Trace
it through the chemical driver, mass-rate normalization and face reconstruction
to determine whether repair belongs here or upstream. Verified RIGHT physical
current and flux normalization retain this dependency; the preserved end-pencil
and frequency results keep their existing domains.

All ten broad TODOs remain. Section 1 is supplied/unfalsifiable, and c2
cross-engine operand debt, global/sheet/exceptional coverage, complete
scattering, profile-frequency bound poles, controls/bookkeeping and export
remain open. No full four-case regeneration, export, S10/Lean/authority edit,
review/comparator/Wolfram/downstream run or commit was made. New outputs remain
uncommitted for DataLad/git-annex storage at a requested checkpoint.

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
