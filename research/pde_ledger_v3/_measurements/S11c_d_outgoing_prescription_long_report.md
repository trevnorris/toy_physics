# Longer outgoing-prescription continuation: partial completion

The explicitly approved longer run ended at its 3540-second native deadline
on 2026-09-27, with **217 complete journal operations: 207 restored and 10
new**. It finished the previously interrupted block-17 adjugate derivative,
saved that block's residue and full-block checks, then timed out in the
corresponding block-16 derivative evaluation. No outgoing kernel is accepted.
No retry, external review or new scientific validation was launched.

The [checkpoint](S11c_d_outgoing_prescription_long_checkpoint.json) is the
complete byte-hash inventory and inspection receipt. Its SHA256 is
`382d7e474aca66487da7dc49e323efbf5de0d6a68b0d2f9c5f9295108ed3ab40`.
The run remains at
`_scratch/s11c/s11c-d-outgoing-prescription-20260927/extended`.
All earlier production, continuation and reference artifacts retain their
original paths. Exact source/manifest/gate/guard/authorization snapshots are
preserved separately at sibling `extended-timeout-source`.

## Actual outcome and resource cost

The worker started at 19:46:34 UTC and the coordinator finished at 20:45:34 UTC.
Worker time was 3540.028 seconds; supervisor time 3540.288 seconds; guard time
3540.409 seconds. The native alarm raised `TimeoutError` and saved the failure
receipt. Guard, supervisor and coordinator report exit 1; the guard's stop
reason is null because the worker's native deadline fired before the
3600-second outer ceiling. This is a timeout, not a failed residual test.

The 1771 resource samples record a peak of 137,863,168 bytes (131.48 MiB),
zero swap and zero memory events. At most three processes/tasks were observed.
Actual containment was 2 GiB memory, zero swap, 32 tasks, CPU affinity `[15]`,
nice 15 and all six native thread settings equal to one. The task-local
duration exception is pinned and recorded in both guard invocation and native
containment. The shared guard hash is unchanged.

Final supervisor invocation and `active.json` agree on failure; guard stdout
matches its outcome JSON, and resource-guard stdout matches the supervisor
receipt. Scientific stdout is empty. Scientific stderr is a 5636-byte timeout
traceback; infrastructure stderr is empty. There is no `checks.json`, so a
successful stdout/checks identity cannot be asserted. Worker final posthash
and operation-index receipts were not reached. The checkpoint explicitly
labels its externally collected inventory and posthashes instead. The hook
finished and queued this completion event. All four recorded PIDs are gone
and the service cgroup is absent.

## Preserved progress and its limits

All 207 prior journal entries are labelled `RESTORED_PRIOR_COMPLETE_RETURN`
and link exactly to the previous checkpoint; their operation functions were
bypassed. The new block-17 derivative evaluation took 1721.252 seconds
(28 minutes 41 seconds). Its residue construction, source-at-pole evaluation,
physical/native bulk-coordinate joins and subsequent full-block checks
completed. Block 17 has nullity two and current direction minus one.

The worker's saved block-17 summary reports all eleven residual norms below
`1e-8`, with the largest `5.744387398727955e-16`. In particular, the relative
cofactor-versus-projected-block residue norm is that largest value; the actual
pencil join is `4.488947401037362e-16`. The direction mutation residual is
`2.8284271247461903`, and changing the delta sign moves its coefficient by norm
`3.7717975547957288`. Source inspection and reaching block 16 establish that
the block-17 determinant-order, Laurent-order, residual and mutation guards
passed. These are preserved worker results, not a fresh independent scientific
validation. Inspection read JSON and hashes only; no scientific pickle was
restored or residual recomputed.

All 26 denominator-proof summaries and 50 tail summaries match the previous
continuation. The stored denominator coefficients are nonnegative with at
least one positive term in the selected real or imaginary component. The
determinant certificate uses imaginary coefficients
`[5760000000000, 5760000000000, 0]`. The tail reciprocal-power counts remain
1:2, 2:18, 3:14, 4:16; only entry (2,3) has order one on both tails. The
18-candidate census still selects blocks 16 and 17 and excludes the opposite
real-sheet lifts 14 and 15. Restored inverse/adjugate and excluded-sheet joins
remain preserved, without being recomputed during inspection.

New operations 213–216 saved block 16's determinant values through order two
and its zeroth adjugate value. Its later order guards were **not reached**.
Operation `217-block-16-adj-value-1` began at 20:15:22 UTC and has a complete
70,328-byte input but no return. The traceback again locates the cost in
general `simplify` zero/sign assumptions, minimal-polynomial construction and
integer polynomial factorization after exact pole substitution. There is no
entry-level completion record inside that interrupted matrix operation.

The saved-density join, explicit `Ne(z,zp)` domain artifact, PV integral and
delta-term assembly were not reached. The first block's success does not
complete the two-block prescription, a coincident-point distributional
extension, a full Green operator, FORM, two-asymptote response or A11/A12.
The unchanged physical slice remains evanescent, with no new radiation witness.
The earlier literal reviews remain Claude NEEDS REVISION and Grok CLEAR FOR
THIS BOUNDED STAGE; no new independent CLEAR is claimed.

## Exact saved-result routes

Paths below are relative to `extended/complete`; every operand and return is
indexed with bytes and SHA256 in the checkpoint.

| Artifact | SHA256 |
|---|---|
| `operations/207-block-17-adj-value-1/value.pickle` | `6f8577a1d8fed69833c54a04754e44079528a8ee409e67c2f60e80a2e682e1b9` |
| `operations/208-block-17-exact-residue/value.pickle` | `4c263e7ad9e66dbe247567f97f434b92a1bc1143e032341389355216d3857ca2` |
| `block-17-projected-operands.pickle` | `e19bb4c4f2f27cb7bf5f4e50f5a061170b83f7d6269d2bf98ef96bc78d6452d5` |
| `block-17-checks-and-mutation.pickle` | `fc5fe893dd829ff4c6568d53a1fe763d98a5adcc80614730c0f1e07bfcad4e5d` |
| `block-17-checks.json` | `01dcdc3d13085ccfb7bf7f38ca6a0e6c1ccf8670fcc0e57734a27614dd631240` |
| `operations/217-block-16-adj-value-1/input.pickle` | `57f648984b31ea08667e46b20c8b8dbe502741797f306c1f7b723ad1fa72c4fc` |

The checkpoint inventories 1024 result files / 4,883,752 bytes. All 470 source
routes, canonical paths, byte counts and hashes match the consumed manifest.
Inspection also verified all 181 original failed-production artifacts, 980
previous-continuation artifacts and 497 accepted reference artifacts against
their existing checkpoints. The protected builder suffix remains
`f01024512c9a028e7e524e5a95b5c455be42e125d10c36f94e20493109b013c2`.
The new 1,414,132-byte checkpoint has an explicit `.gitattributes` diff
exclusion. Large binary scratch artifacts remain at their original paths.

## Next execution boundary

The first unfinished scientific operation is block 16's adjugate derivative
at its inherited exact pole. More time did complete block 17, so the previous
timeout is not evidence that this calculation cannot finish. The remaining
cost is still unknown. No further duration extension or run is inferred from
this completion notice.

Any separately approved continuation must reuse all 217 journal returns and
also restore the complete saved block-17 contraction/check/mutation bundle.
Those block checks occur outside `Journal.op`; merely increasing the journal
resume prefix would repeat them. Add an explicit saved-block restore boundary
before a future launch. Continue only block 16 and the unfinished density,
domain and integral work, with unchanged inputs and full checks. A narrowly
scoped rational-radical evaluation/persistence change remains an implementation
option, not an authorized new method or a change already made. Any scientific
validation must use its own guarded saved-output job. No script remains active.
