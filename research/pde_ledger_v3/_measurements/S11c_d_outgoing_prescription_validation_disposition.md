# Saved-output validator Boolean comparison correction

The outgoing-prescription constructor completed successfully. Its 1069 saved
files are unchanged. The first validator then failed because its own
`tail['decays'] is True` test rejected a saved SymPy `BooleanTrue` object.
This is a validator identity-comparison error, not evidence that the inverse
tail fails to decay. The user explicitly approved the corrected attempt;
it has now passed. The original failure remains unchanged.

The [failure checkpoint](S11c_d_outgoing_prescription_validation_failure_checkpoint.json)
preserves all 217 completed validation operations and 901 files / 1,197,482
bytes, including every input and return. The exact failed validator, manifest,
gate and launcher are frozen in
`_scratch/s11c/s11c-d-outgoing-prescription-20260927/validation-failed-source`.
Original `saved-validation` output is untouched. The failure was inspected
during the constructor-completion turn; a later duplicate watcher notice
must not cause a retry or new launch.

## Actual evidence

The run started at 00:28:19 UTC on 2026-09-28 and stopped after 2.878 supervisor
seconds (3.004 guard seconds). All resource controls were verified. Peak
memory was 93,073,408 bytes, with zero swap and memory events. Infrastructure
stderr is empty; scientific stderr contains the 838-byte failed-check trace.
Scientific stdout is empty and final checks/posthash/index receipts were not
reached. The checkpoint labels externally collected hashes and operations.
No job PID or service cgroup remains.

All 30 structural/source/unit/domain checks and all 26 saved denominator
certificate checks passed. On `inverse-0-0-tail-minus`, these six checks passed:
the rational operand join, tail side, positive reciprocal coordinate,
positive integer power (two), saved leading term, and saved series leading
term. Only `savedDecay` was false in the validator's result.

Read-only `pickletools` opcode inspection of the saved tail value identified
the `decays` key followed by `sympy.logic.boolalg`, `BooleanTrue`, `STACK_GLOBAL`
and `NEWOBJ`. No scientific pickle was restored for this diagnosis. This
matches the constructor's expression: its last `power > 0` comparison returns
a SymPy Boolean, while the constructor correctly used its truth value. The
validator incorrectly demanded identity with Python's distinct singleton.

## Explicitly approved corrected attempt, completed

The approved [v2 worker](S11c_d_outgoing_prescription_validate_v2.py) changes
the saved flag test to membership in `(True, sp.true)`. It does not accept
arbitrary truthy or unresolved expressions, and it retains all six mathematical
tail checks. For the already completed failed check, it restores the original
return and performs one new operation that changes only this flag result,
records the actual class name, and confirms the other check values are
unchanged.

All 217 prior complete validation operations are restored in order with their
functions bypassed. The failed tail-check return is preserved, not overwritten.
After the separate literal-Boolean correction, the validator continues the
remaining 49 tails, both block residual/mutation checks and explicit PV/delta
assembly. It reuses the exact constructor outputs; it does not reconstruct
derivatives/tails, import producers, solve modes/roots/LU, or evaluate integrals.

The [preparation receipt](S11c_d_outgoing_prescription_validation_v2_preparation.json)
and [11 metadata/static checks](S11c_d_outgoing_prescription_validation_v2_static_checks.json)
pin 1069 unchanged candidate files and 901 unchanged first-validation files.
Remaining denominator/block/assembly validation functions are AST-identical;
native containment is unchanged. The launcher was checked to refuse its
pending gate before creating a run directory or arming a job.

- Worker SHA256: `45fd003457e42ec067e23f54d19b263fe60a41b8420c57827d908c9146c821d8`.
- Manifest SHA256: `b298db3dd446374f8b7a9603e2195f3dd410e448c006d6942bc5e074610e40c8`.
- Completed root: `_scratch/s11c/s11c-d-outgoing-prescription-20260927/saved-validation-v2`.
- Limits: one existing shared guard/supervisor, 900 seconds outer / 840 native,
  2 GiB, zero swap, one CPU, nice 15, 32 tasks, one native thread.
- Hook: armed before launch, targeting session
  `01a0e01b-ef84-7192-817f-584cda5d339b`; no model polling or automatic retry.

The user replied **Approved** to this prepared correction. The
[authorization](S11c_d_outgoing_prescription_validation_v2_authorization.json)
and [gate](S11c_d_outgoing_prescription_validation_v2_gate.json) pin the single
corrected attempt. The exact pending gate is preserved in scratch history.
This was an explicitly approved attempt, not an automatic retry.

The [launch receipt](S11c_d_outgoing_prescription_validation_v2_launch.json)
and [checkpoint](S11c_d_outgoing_prescription_validation_v2_checkpoint.json)
record completion at 01:38:06 UTC on 2026-09-28. All 217 prior return files
were restored byte-for-byte; 224 new operations completed. The separate first
tail correction records `sympy.logic.boolalg.BooleanTrue`, the original
false flag, all other checks unchanged, and all seven corrected checks true.
All 30 structural checks, 26 denominator certificates, 50 tails, both full
block residual/mutation checks and seven assembly checks passed.

The validator took 3.930 worker seconds, with 89,706,496-byte peak whole-job
memory, zero swap/events and all required containment verified. Scientific
and infrastructure stderr are empty; stdout and checks are byte-identical.
All 1069 candidate files, 901 prior validation files and 498 constructor
input routes remain unchanged. All job PIDs and the cgroup are gone. The
final checkpoint preserves 441 complete operations / 1854 files / 3,490,903
bytes, exact source snapshots, logs and hashes. No new physics campaign or
external review ran. The [current report](S11c_d_outgoing_prescription_candidate_report.md)
and [current checkpoint](S11c_d_outgoing_prescription_checkpoint.json) close
only the bounded saved-output prescription stage.

No fresh independent CLEAR, full Green/FORM/A11/A12 clearance, diagonal
extension or radiation witness is claimed. Shared guard, Lean and the
protected builder suffix are unchanged.
