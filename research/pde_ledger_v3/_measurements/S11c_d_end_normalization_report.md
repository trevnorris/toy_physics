# S11c-d endpoint normalization execution

RIGHT normalization is validated and published for the supplied LAB_HELD /
RHO4_CONSTANT case. The 31.8-minute constructor computed 18 isolated-root/lift
candidates, all 22 basis directions, and 18 invertible adjoint field maps. Two
two-dimensional subspaces have defined physical current normalization: two
incoming and two outgoing basis directions at RIGHT. Independent frequency
and full-coordinate-projector joins agree within 4.0e-14.

Validation covers 5,990 original tags, 28,146 metadata paths and 5,850 numerical
residual scalars. The 250 coefficient and 50 symbolic map residual scalars are
zero. Nine raw balance/reconstruction families remain nonzero, up to 0.01342
in the declared numerical unit frame. The independent pairing remainder and
its regular-sheet derivatives contain only eta-squared coefficients. All 275
remainder decomposition/projection/limit residual scalars are zero. Contracting
those operands through every full right and adjoint subspace accounts for the
raw residuals, with a largest difference of 2.15e-13. This establishes the stated
retained-order accounting; finite-contrast balance is not claimed to be exact.

The original normalization transcript and a separate dimensioned, graded
remainder transcript are published together. Their fingerprints, hashes,
source snapshots, complete per-mode comparisons and raw diagnostics are in
[S11c_d_end_normalization_right_thickness_repair_checkpoint.json](S11c_d_end_normalization_right_thickness_repair_checkpoint.json).
Frozen physical producer files remain unchanged. Checker repairs handle JSON
representation and input-map order; the validator checks every binding before
replaying the original ordering. The launcher now attaches logs only for the
actual constructor and persists child completion before log bookkeeping.

All eight fresh source/frequency/pairing prerequisites remain committed through
432db7e7, with 802 zero retained pairing scalars at each of RIGHT, LEFT and
REFERENCE. RIGHT is committed at 0da0f746 and preparation at 8ca2f3b0.
LEFT construction finished successfully in 178.18 seconds and validation in
56.45 seconds, both with empty stderr. Its 18 candidates, 22 basis directions,
18 adjoint field maps, and two normalized two-dimensional subspaces passed the
full checks. The original transcript has 5,990 tags, 28,146 metadata paths and
5,850 numerical residual scalars. No raw norm exceeds the 1e-8 diagnostic
threshold. Its independently computed discarded balance matrices vanish;
275 exact remainder checks are zero and the largest contraction difference
is 8.21e-13. The LEFT original and remainder transcripts are published with
source pins in [the checkpoint](S11c_d_end_normalization_left_thickness_repair_checkpoint.json).

LEFT is committed at b271c71b. The [Mathematica audit](S11c_wolfram_repair_audit_report.md)
found a duplicate physical/reference pressure shift in c2 N6; the three mechanical
issues are absent on the audited b domain. The [repair plan](S11c_wolfram_pressure_trace_repair_plan.md)
is user approved. The complete focused repair check has 218 zero residuals;
native four-case validation has 9,968,256 zero required numerator evaluations.
The repaired main transcript is committed at `5acfdf30` and its complete annex
payload hash is verified. REFERENCE construction started at 23:03:14 UTC using
the prepared input/current/pairing packets. A refreshed preflight verified all
16 joins, including unchanged pairing constructors and prepared sources. Its
completion watcher delivered the successful constructor exit after 166.28 seconds
(202,188 KiB peak RSS). Construction records 18 candidates, all 22 basis directions,
18 invertible adjoint field maps and two current-normalized subspaces. Its 250
coefficient and 50 symbolic residual scalars are zero; the pre/post-emission
packet hashes match. Full-subspace validation and independent remainder accounting
are now launched with a separate completion/error watcher. These constructor
results are provisional until validation passes.
Accepted SymPy results are unchanged.
Existing sheet/exceptional-domain limits remain.
Owned scripts use user-authorized local completion/error wake-ups; no recurring
checks or model polling are active.
