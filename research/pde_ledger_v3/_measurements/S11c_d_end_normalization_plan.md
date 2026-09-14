# S11c-d current and adjoint normalization after the thickness repair

The user authorized continuing and committing each substantive step on
2026-09-14, with a stop before another upstream physical repair. This extends
step 6 of the thickness-coordinate repair plan and the existing native
`ModalCurrentSubspaces` and `AdjointCurrentMap` constructions.

Queue disposition, 2026-09-14: LEFT is validated and committed at b271c71b.
The [Mathematica audit](S11c_wolfram_repair_audit_report.md) found a duplicate
physical/reference pressure shift in c2 N6. Its [repair plan](S11c_wolfram_pressure_trace_repair_plan.md)
is user approved. The repair passed its focused tests and native four-case
regeneration/validation: 218 exact focused residuals and 9,968,256 native
numerator evaluations are zero under the adopted covariance criterion.
The native annex publication is committed at `5acfdf30` and hash-verified.
REFERENCE construction and validation completed successfully. RIGHT, LEFT and
REFERENCE normalization are now accepted on the supplied case; publish/commit
the REFERENCE checkpoint before variable-profile matching. The license allows at most two
Mathematica scripts across sessions; heavy CAS jobs remain serial.

Prerequisites completed: all eight endpoint stages are validated, published
and committed through 432db7e7. RIGHT pairing is committed at 7c2001b6, LEFT
at 86e61755, and REFERENCE at 432db7e7. Each has 802 retained residual scalars,
all zero. The native source context extension follows this completed queue.
Read `_scratch/s11c/s11c-end-normalization-20260914/active.json` at the
repository root for the current normalization command and durable live logs.

1. Add an explicit context parameter to the two existing modal/adjoint
   emitters, preserving their reference defaults and calculation bodies.
   No energy, balance, mode or derivative formula changes in this step.
2. Source-check the new pairing packet and native spectrum. Join the complete
   physical pencil, material/profile bindings, grade origin, isolated disks
   and both normal-momentum lifts. Reuse the computed raw objects together
   with their independent-grade retained/remainder checkpoint; do not replace
   a physical endpoint with the reference calculation.
3. Run RIGHT, then LEFT, then a fresh REFERENCE regression, sequentially.
   Use the supplied finite-contrast homotopy at physical ends and the zero
   background origin only at REFERENCE. Compute every full modal subspace,
   frequency/current form and adjoint field map using the existing engine.
   Compare full coordinate projectors and the frequency derivative to the
   independently validated end-frequency packet. Basis rotations are allowed;
   missing directions or rank/property mismatches are not.
4. Apply the specified outward orientations to computed signed currents and
   emit incoming/outgoing basis-column records. Preserve every closed,
   evanescent, wrong-sheet and unresolved candidate as a domain record.
   The reference regression has no physical end orientation.
5. Validate literal residuals, heavy-object fingerprints, every emitted
   metadata path, full-basis ranks, source joins and transcript census.
   A finite-contrast raw balance discrepancy must be investigated against its
   recorded truncation remainder; it cannot be silently accepted as closure.
   RIGHT's completed investigation finds only eta-squared coefficients in the
   saved discarded balance and its regular-sheet derivatives. Contract these
   independent operands through every full right and adjoint basis, emit the
   raw residual, remainder and difference separately, and require the
   difference to meet the numerical diagnostic threshold. Preserve the raw
   nonzero values in the accepted inventory. Apply the same calculation to
   LEFT and REFERENCE; do not infer their results from RIGHT.
   Publish each accepted transcript atomically and save with DataLad/git-annex;
   commit ordinary source, plan, report and inventory files with Git. Verify
   full payload hashes after annexing before the next stage.

All new runs live under the physical repository path
`/var/projects/toy_physics/_scratch/s11c/s11c-end-normalization-20260914/`.
Keep failed attempts and their source snapshots. The native b/c1/c2/d producer
and existing endpoint packets remain reusable subject to constructor/input
joins; a metadata or emission adapter does not require regenerating them.

After accepted normalization, continue one-case variable-profile matching
under the unchanged solver/export contract. This checkpoint does not claim
global exceptional coverage, a continuum re-expansion, a complete S-matrix,
or section 3b profile-dependent frequency poles. Those remain program work.

REFERENCE execution preparation: atomic modal and adjoint construction packets
are now saved before emission, alongside the established final packets. The
validator checks their exact structural equality. Progress identifies emission
stages and resource use; periodic stack diagnostics stay in the run directory.
The accepted REFERENCE input/current/pairing hashes passed a read-only preflight.
The fresh REFERENCE scientific run and independent remainder validation have
completed; pre/post-emission construction packets agree exactly.
Use a local completion/error watcher for each owned long-running script;
no recurring checks or model polling while a healthy job runs.
