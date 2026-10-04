CLEAR FOR THIS EXPLICIT NUMERICAL-PREFIX RECOVERY METHOD

This clears the method as a basis for a separately assessed concrete build. It does not clear any worker, and it gives no science result, action, gate or acceptance. The conditions C1–C9 below must be met by that build. No blocker applies to the method itself.

## What I read
- **Method and prototype:** `guide.md`, `method.md`, `recovery.py`, `tooling-tests.py` and `tooling-tests.log`. The log shows 54 passing tests, which I did not re-run.
- **Source files:** `original-numeric.py` in full, and `request-index.py` and `evidence-store.py` in full.
- **Records and exports:**
  - Record 2 in full.
  - Records 3 and 486 by targeted field extraction, not a full read. I checked route, precision, kind, cell, selectedRule, ownRoute, side and ruleIndex.
  - `contract.json` header and pins (lines 1–12 and 10928–10939).
  - `export-index.json` head, and the PENDING count in `requests.json`.
  - Not read in full: records 0, 1, 485 and 1449, the full `requests.json`/`receipts.json`/`namespaces.json`, the views, the `failure/` tree, `original-worker.py`, `numeric-evidence-fix.py` and the tail/geometry files.
- **Not done:** I did not compare the large context blobs across records 0, 1, 3 and 486 byte for byte. That equality rests on what `audit_prefix` and `join_frame` would enforce when run, not on my own check. Nothing was executed.

## Real-record joins
- **Record 2** is named `1/outer-node-input`, in the A24 baseline namespace `9d3e6…`, with side=0, ruleIndex=0, K27/T122, n=0, slabId 0, selectedRule and ownRoute both A24.
- **Records 3 and 486** are both `inner-panel` requests with route A24 and precision 30. Both have cell `{cellId:0, clipped:true}`, an outer parent with side 0 and ruleIndex 0, and ownRoute A24. Record 3 selects A24 and record 486 selects A48. This matches `join_frame` and the A48-inside-A24 claim.
- **Serial arithmetic** checks out independently:
  - 1 (node) + 144 (inner-node-input) + 144 (root-operands) + 288 (formula-checks) + 3 (batches) = 580.
  - Total records: 1 + 2×434 + 1 pending input + 1 node + 579 other puts = 1450.
  - So node serial 1 is skipped, the new evaluator must start at 580, and the pair becomes `581/inner-positive-panel-pair`. `after_put` checks exactly this name.
- **PENDING count:** `requests.json` has one PENDING row.

## The scalar-accumulation boundary
This boundary is sound and needs no different mechanism.
- **What actually happens:** The design re-enters the unfinished outer callback in a fresh process. The completed panels, nodes and transforms are replaced by byte-exact lookups, and everything above them is recomputed.
- **Why the frame is pinned:** Each lookup requires a byte-equal descriptor. These descriptors embed the encoded MP tuples for m, carrier, a, b, outerParent and cell. The node is compared byte for byte against record 2. A diverged reconstruction therefore refuses before any callback.
- **What is not restored:** The `local` and `result` accumulators were never published, so recomputing them is reconstruction. The protocol labels it that way: `restoredPublishedReturn:False` and `completedCallbacksInvoked:False`.
- **Scope of the repeated arithmetic:** Only `paired_indicator` and the `Estimate` addition on already decoded saved values are repeated. That is a deterministic add, subtract and abs at 30 digits.
- **Limit:** No original pair record exists, so there is no oracle for the reconstruction. Its correctness rests on determinism and on an unchanged source expression (C2, C6).

## Concerns found
**Prototype defects:**
- **P1. No refusal latch.** After a raised refusal, `state` stays PARENT, CHILD24 and so on. A caller that catches `ValueError` can retry with the correct descriptor in the same process. The method says refusal stops the run. The protocol should latch to a terminal FAILED state at `recovery.py:330-428`.
- **P2. Pair contributions are checked for key set only** (`recovery.py:408`). Some checks need no arithmetic. For `ownOrder` 24, `contribution[k].value` must equal the A24 `components[k].value` encoding exactly. Each contribution must also have the encoded value/error shape.
- **P3. Loose sequence check.** `after_put` requires `sequence >= records` (`recovery.py:420`). It should require the exact next sequence, given the provenance record's position (C3).
- **P4. Hard-coded flag.** `originalLeafTableWasEmpty: True` in `prefix_preserved` (`recovery.py:251`) is a constant, not a measurement.
- **P5. Incomplete preservation check.** `prefix_preserved` does not compare `sqlite_master`/schema. It does not audit the new suffix chain, which should link from head `97034b…` with no orphan chunks. It also does not report new PENDING rows.
- **P6. Test gaps.**
  - No test mutates the pair's `pairedPanels`, `ownOrder` or component set, or the name or serial.
  - There is no test of a kill between commit and `after_put`.
  - There is no test of NEW-state reuse after recovery.
  - Fixtures use serial 1 rather than 580.

**Future concrete-build duties:**
- **C1. Decode/encode round trip.** `before_put` requires `encode(decode(saved)) == saved` for the pairedPanels. If that fails, recovery would fail late, with the clone stuck in PAIR and the parent PENDING. Pre-verify that round trip on records 485 and 1449 before the clone is created.
- **C2. Wrapper integration.** The numeric wrapper needs real call sites, and a fresh AST certificate, since the delta is no longer only the evidence boundary.
  - It must not emit the original `exact-request-reuse` put in recovery states. `before_put` would refuse it. In NEW the intended behavior must be specified, because the protocol records reuse in memory only, not as the original's serial-consuming record.
  - The evaluator serial must start at 580.
  - `contracted_leaves` must not be re-created.
- **C3. Execution provenance and durable ledger.** The provenance record needs a fixed location.
  - It cannot go through the protocol, which refuses every put before NEW.
  - Commit it directly to the clone, in its own namespace, before `claim_parent`.
  - Commit a durable recovery-events record immediately after the pair. The `events` and `attempts` are in memory until `final_record`. An OOM kill under the 4 GiB guard would lose them.
- **C4. Path identity.** The descriptors embed absolute rule paths. For example, record 2's `ruleReceipt.path` is `/var/projects/.../inner-preparation/A-GL24.json`. The contexts must be restored from record 0, never regenerated from new paths, or every lookup will refuse.
- **C5. Namespace configuration.** The index's `maximum`, `reserve` and `cache` values (8388608, 21474836480, 0) feed into the namespace descriptor. They must reproduce the original namespaces exactly. The prototype does enforce this via `exact()` in `namespace`.
- **C6. Numerical environment.** The concrete build should pin the mpmath version and backend, per the disclosed backend duty. It should add typed-geometry checks and a full-AST delta limited to the evidence key and the reviewed recovery hooks.
- **C7. Final audit.** Do a full `EvidenceReader.audit` of the clone at finalization, plus the hashes in P5.
- **C8. No resume after a crash.** A crash in NEW leaves nested PENDING rows that ordinary rules will refuse. That matches "no generic automatic resume". It also means any later failure needs its own reviewed method.
- **C9. Pair-contribution verification.** Add the P2 check, and have the wrapper verify the pair contributions against the original expression.

**Mathematical concerns:** None new. The design leaves the equations, windows, precision, rules, tolerances and tails untouched.
- The A24 and A48 paired estimates are not independent, and the Gaussian identities are shared between routes.
- Error figures are empirical indicators, not rigorous bounds.
- H/flat/heightPV/slope are pending, and J/direct is not a full packet or a leakage/loss/current.

## Other checks
- **Latent integer-key sites:** In `original-numeric.py` I found `pairedPanels` to be the only dict with integer keys that reaches `put`. The `adaptive`, `outer_decision`, `addressed` and `action` dicts all use string keys. This came from reading the source only. The concrete build should still confirm it.
- **Clone path:** The checks run in the right order (admission audit, then disk reserve, then exclusive `xb` copy, then file fsync, directory fsync and a three-way SHA check). A partial copy is preserved and never reopened.
- **Request-index behavior:** `claim_parent` raises before any `begin` on a mismatch, and only one claim is allowed. Ordinary `lookup` and `begin` still refuse PENDING rows. Unreachable callbacks hold for both child mismatches and exact matches, because the child branch never calls `callback`.
- **Prefix immutability:** Completed request rows are compared exactly. The only allowed change is the one-way PENDING-to-COMPLETE transition, with input fields unchanged and an output sequence of at least 1450. New adaptive leaves are allowed at finalization.