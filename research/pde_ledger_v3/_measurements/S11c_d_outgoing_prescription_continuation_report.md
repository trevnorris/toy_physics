# Outgoing prescription continuation: partial results preserved at time limit

Status: **unfinished; no outgoing kernel accepted**. The explicitly approved
continuation reached its 900-second wall-time limit on 2026-09-27. The original
denominator-check failure and both independent review reports remain preserved.
There has been no further attempt, review submission or scientific replay.

The [checkpoint](S11c_d_outgoing_prescription_continuation_checkpoint.json)
pins all **980 saved files / 3,883,186 bytes**, the 207 complete operation
receipts and the one incomplete operation. All **137 consumed file hashes**
remain unchanged. Inspection used JSON, hashes and source text; it did not
restore scientific pickles or independently re-evaluate saved expressions.

## Completed and retained

All first **43 operations** are labelled `RESTORED_PRIOR_COMPLETE_RETURN`,
with their prior records matching the pinned first-run operations. Their
functions were bypassed. The continuation then completed **164 new operations**.

The saved worker summaries report:

| Work | Saved outcome |
|---|---|
| Determinant and 25 adjugate denominator proofs | All 26 passed the positive-radical nonzero certificate. Actual real/imaginary polynomial coefficients and signs are saved. |
| Both tails of all 25 inverse entries | All 50 passed the decay requirement: 2 tails at reciprocal-momentum order 1, 18 at order 2, 14 at order 3, and 16 at order 4. |
| Full inherited candidate classification | All 18 recorded; real on-sheet blocks 16 and 17 selected, each nullity two. |
| Excluded real-normal candidates 14 and 15 | Exact opposite-sheet joins completed before the subsequent work. |
| Determinant derivatives through order two and adjugate first derivative | Complete operands and returns saved. |
| First pole block's determinant values through order two and adjugate value | Saved; the later full residue guard was not reached. |

These are retained partial computation results. Their raw operands and returns
are available for the eventual saved-output validation; this metadata inspection
does not turn them into a completed outgoing-kernel result. The two order-one
tails are entry `(2,3)` on the negative and positive real axes. The separated-point
domain remains the intended scope; no diagonal distributional extension is claimed.

## Exact stopping point

Operation **207 — `block-17-adj-value-1`** started at 19:14:12.444163 UTC.
Its complete 68,063-byte input is retained at:

```text
_scratch/s11c/s11c-d-outgoing-prescription-20260927/continuation/complete/
operations/207-block-17-adj-value-1/input.pickle
SHA256 5e3acd42972a5ddd61766b3e4a329f91f2c961c66a211981ee3c47cf450a3b93
```

There is a start receipt and no completed return. The worker was evaluating
the saved adjugate derivative at an exact pole using
`a.subs(kn, b).applyfunc(sp.simplify)`. The traceback shows `simplify` querying
`expr.is_zero`, entering SymPy's algebraic-number minimal-polynomial fallback,
and spending the remaining time in polynomial factorization. The native
900-second alarm wrote `failure.json`; the outer guard also stopped the job
for its wall-time limit. This operation occupied the last approximately nine
minutes. The largest earlier completed operation took about 24 seconds.

No block residue, block comparison/mutation report, saved-density join,
prescription-domain artifact or final candidate integral was produced. The
scheduled determinant/adjugate order guards after this incomplete matrix
evaluation were not reached and must not be reported as passing.

## Resources and outcome integrity

The actual cgroup enforced 2 GiB, zero swap and 32 tasks, with CPU affinity
`[15]`, nice 15 and all six native thread limits set to one. Across 452 resource
samples, peak whole-job memory was **176,340,992 bytes (about 168 MiB)**; swap
and all memory events remained zero. The guard exited 124 after 900.349 seconds,
and its child outcome explicitly records `wall-time limit`. Native failure
runtime was 900.014 seconds.

Scientific stdout is empty; stderr is the 5,104-byte timeout traceback. There
is no successful `checks.json`, so stdout/checks identity is unavailable.
The guard's stdout matches its saved outcome. Final supervisor invocation and
worker operation-index/posthash receipts are missing: the guard ended the
supervisor before final bookkeeping. Its `active.json` retains a stale
`running` value. Direct process/cgroup inspection found all recorded run PIDs
gone and the service cgroup absent; **no job remains active**.

The checkpoint reconstructs the operation inventory from the durable per-call
receipts and records externally collected posthashes. It does not overwrite
the missing worker receipts or imply a clean scientific finish. Worker, input
manifest and gate snapshots are preserved in `continuation-timeout-source/`.
The protected builder suffix remains
`f01024512c9a028e7e524e5a95b5c455be42e125d10c36f94e20493109b013c2`.

## Next implementation boundary

The bottleneck is now specific: exact evaluation of the already saved adjugate
derivative at the first pole. The next bounded repair should normalize/reduce
those saved expressions in the existing rational radical coordinates before
algebraic point substitution, and preserve entry-level work so an expensive
entry cannot discard a complete matrix operation. It should consume the 207
complete returns and the saved interrupted-operation input, without replaying
the denominator, tail, census or derivative computations.

That repair is not yet implemented or authorized to execute. It must preserve
the exact residue and independent full-block comparison requirements, use the
same physical point and resource limits, and receive the applicable concrete
gate and explicit continuation approval before a new worker. There is no
automatic retry, limit increase, new mode/root/current/producer/profile job,
or external review cycle.

The accepted reference inverse and regular density remain available. Full
outgoing Green, two-asymptote response, FORM, A11 and A12 remain unfinished.
The literal reviews remain Claude **NEEDS REVISION** and Grok **CLEAR FOR THIS
BOUNDED STAGE**; no fresh independent CLEAR is claimed.
