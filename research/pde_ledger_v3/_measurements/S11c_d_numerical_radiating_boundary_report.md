# Central boundary gate: saved progress, native timeout

The first numerical radiating boundary gate stopped at its native 840-second
limit. There is no recorded mathematical check failure, but the boundary gate
is incomplete. Nothing is running from this launch, and no continuation or
scientific retry has been launched.

At omega = 3, 17 of 18 LEFT candidates completed both prescribed upper-half
frequency paths. Their largest recorded original-pencil residual is
4.2891e-13; the largest projector disagreement between paths is 2.8922e-13.
All 18 LEFT zero-contrast/reference joins passed. The literal contrast-independent
LEFT pencil reused 54 candidate records without Newton replay.

Candidate 17 completed its coarse path. Its fine path has a fully joined saved
point at 1.85 + 0.025i and a subsequent completed Newton check at 1.875 + 0.025i.
The alarm interrupted persistence of that next joint-sheet return. The return
is absent and remains unaccepted; its complete input is preserved. RIGHT was
not reached. There are zero finite solves, current/face results or transverse
loss estimates, and no result at omega = 2.3 or 4.

The saved LEFT transverse candidate 16 has real normal momentum approximately
2.43926218 (imaginary roundoff 3.8e-17) at omega = 3. Its current sign has not
been evaluated. This is useful partial continuation evidence, not a lossless-end
premise or radiation/loss measurement.

## Integrity and cost

Worker time was 840.270613 seconds. Peak cgroup memory was 264,159,232 bytes
(251.9 MiB), with zero swap and zero memory events. The ordinary 2 GiB, one-CPU,
32-task, one-thread limits were enforced; the native alarm, not a memory or
host-reserve stop, ended this job. All scientific and guard stderr was empty;
stdout and checks were byte-identical.

Lightweight byte/hash inspection verified all 71 input posthashes, 67 source
snapshots, four copied accepted packets, and all 5,001 complete journal
input/return receipts. The archive contains 10,004 unique, accounted-for members;
no scientific payload was deserialized during inspection. The complete result
has 36 files / 30,729,309 bytes. Exact artifact routes and hashes, partial census,
checks, resources and pending operation are recorded in
`S11c_d_numerical_radiating_boundary_checkpoint.json`.

The preserved archive is
`_scratch/s11c/s11c-d-numerical-radiating-20260929/boundary-01/complete/operations.zip`,
SHA256 `0de8a6e1246a7c125f75d15cd3394cfc24f166f51b55945962602f7091a1ff86`.
It and every earlier result remain in place; scratch is not committed.

## Storage correction and next dependency

The writer reopened the ZIP archive for every operand and return, repeatedly
scanning and rewriting its growing central directory. This was avoidable
instrumentation overhead in my implementation. The timeout stack is inside
that index scan; total time by cause was not profiled, so it would be wrong to
attribute the entire runtime to storage.

A separate immutable-byte SQLite helper is prepared, with durable transactions,
unique names and hash checks. Synthetic tests passed for 512 byte roundtrips,
duplicate refusal, interrupted-transaction recovery, read-only reopening,
existing-file refusal and corruption detection. See `blob_store_checks.json`.
This is a storage-only correction, requiring no physics review. It is not yet
integrated into a continuation and has restored no scientific payload.

The concrete next step would be one guarded continuation that restores the
5,001 complete returns and their actual arguments, resumes the incomplete
joint-sheet operation, then finishes the remaining LEFT and RIGHT work. It
must preserve this run and avoid replaying completed calculations. Runtime for
RIGHT remains uncertain; the storage fix does not promise completion inside
another 840 seconds. The completion event forbids an automatic scientific
retry, so no new run is authorized by this inspection.

The user's removal of the authoring-day cap remains effective; it did not
change this job's native runtime limit. No further Grok pass was made. Literal
round-2 Claude CLEAR and the missing Grok report remain recorded without paired
independent clearance. Physical-loss interpretation, analog-light calibration,
Green/FORM and A11/A12 remain open. Lean, shared guard, incident/review history
and the protected builder suffix are unchanged.
