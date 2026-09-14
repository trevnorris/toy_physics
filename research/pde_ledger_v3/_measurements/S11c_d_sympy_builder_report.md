# S11c-d SymPy builder checkpoint

The user controls continuation. Scheduled follow-ups have been removed; do
not create or recreate automations without an explicit scheduling request.
The already-authorized calculation was left running.

The thickness-coordinate repair and native b/c1/c2/d regeneration are committed
through 8b2e3cf2. All eight fresh endpoint source/frequency/pairing stages are
now validated, published and committed through **432db7e7**. RIGHT, LEFT and
REFERENCE two-frequency pairings each contain 802 retained residual scalars,
all zero. RIGHT is committed at 7c2001b6; LEFT at 86e61755. The earlier three
RIGHT retained discrepancies are preserved in their historical checkpoint.

The next stage is fresh full-subspace current and adjoint normalization at both
ends, with a new reference regression. The prepared runner, validator and
serial-job guard are committed at **5bb2c811**; see the
[normalization plan](S11c_d_end_normalization_plan.md). The native emitters now
accept the actual end/case context; an AST comparison verifies that their
calculation bodies are unchanged. RIGHT normalization is running. Its results
are not yet validated or published.

All important working runs are stored under the repository's `_scratch/s11c/`.
Published `.out` files use DataLad/git-annex with post-save full-hash checks;
ordinary sources, reports and inventories use Git. Historical `/tmp` paths
are compatibility links; see the [storage report](S11c_storage_report.md).
The [execution checkpoint](S11c_thickness_coordinate_execution_checkpoint.json)
points to the live durable normalization record.

The complete two-ended variable-profile S-matrix, profile-frequency bound
poles/residues/overlap, survival, flux bookkeeping, weak coefficients, section 5
controls and final own-row export remain program work. Existing per-root and
full-subspace checks retain their stated sheet and exceptional-domain limits.
The supplied physical premises and c2 operand debt remain inherited premises,
not independently verified by this builder. No final S11c_d_exports.py exists.

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
