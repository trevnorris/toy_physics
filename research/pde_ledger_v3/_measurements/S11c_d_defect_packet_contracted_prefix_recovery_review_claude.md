NEEDS REVISION

I'm not clearing this. I found no equation or tolerance defect. But the method and prototype do not yet stop two silent-recompute paths, and the build spec leaves one duplicate-record decision open. All three need fixing before a separately assessed build can rely on them.

## Verdict

The recovery design is sound in its core. The two completed child panels are reused from the store, and the unfinished parent's pre-request arithmetic is re-executed, not restored. The blockers below are method-text and prototype gaps that can be fixed without redesigning the recovery.

## What I verified against the records

- **Failure site.** `failure/defect_packet_contracted_continue_fix.stderr` matches `original-numeric.py:394`. The `put` of `inner-positive-panel-pair` fails at `canonical` ("string JSON keys required") because `pairedPanels=panels` has integer keys. `self.serial+=1` runs before `encode`, so no record 581 exists. `contract.json:10933` gives `lastCommittedNumericSerial: 580`, and `receipts.json:17381` shows the last serial name is `580/inner-node-batch`.
- **Parent and child joins.**
  - Record 0001 is the pending parent: A24, baseline, precision 30, outer slab `a=-149`.
  - Record 0002 holds the first outer-node input: side 0, ruleIndex 0.
  - Records 0003 (A24 rule) and 0486 (A48 rule) are both `route:"A24"`, `kind:"inner-panel"`, precision 30, baseline.
  - Both carry an identical `outerParent`: the same `originalA/B`, jacobian, `representedPoint`, `ruleNode`, K27/T122, carrier, epsilon and geometry receipt (seq 12424, `new-geometry-27-1-full`).
  - Both inner panels use cell `{cellId:0, clipped:true}`.
  - The inner `m` in both equals the `representedPoint` in 0002, and both inner arguments start at `a=-27` (`-K`).
  - Only `settings.selectedRule` differs (A24 vs A48). That confirms the A48 rule sits inside the A24 own route, not the independent A48 action route.
- **Census.**
  - `requests.json` has 435 rows and exactly one PENDING.
  - `contract.json` pins 1450 records, 434 complete requests, and the reserve, record-size and cache settings to the original values.
  - The mathematical-context namespace `872bea…` is separate from the request namespace `9d3e6b…`.
- **First action.** `original-worker.py:164-171` runs baseline, then plan key, then route, then n, so this was the first action. No earlier outer leaf exists, and `contracted_leaves` is empty. That is why re-entering `panel.work()` from its start reproduces the failed state with no loop-index recovery. The accumulators were empty at failure: `sums={}`, `batch` empty, so the `except` batch flush did nothing.

## The unpublished pair accumulation

The method's boundary is sound and does not need a different mechanism.

- In the original run, `request()` returned `decode(c, encode(value))` for fresh computations. Decoding the stored panels is therefore bit-identical to what the failed run held in memory. The `mpf` tuples round-trip exactly, and `decode` asserts `_mpf_==original`.
- `paired_indicator` plus the first `Estimate(0,0)+local` is a pure function of those two decoded returns. Re-applying it is reconstruction. The method labels it that way and does not claim a restored return.
- No saved frame exists to restore, so reconstruction is the only option. No panel, node or transform callback needs to run again.

## Blockers

1. **A miss on a completed child must refuse, not fall through to `begin`.**
   - `method.md:40` says the wrapper uses ordinary lookup for both inner panels "and otherwise" original new-request logic. `request()` calls `begin` on any miss.
   - A miss would silently recompute a completed panel and its 144 nodes and 288 transforms, mixing old and new work. This could happen if a reconstructed namespace or descriptor differed by even one byte.
   - The ledger requirement ("before continuing beyond the first pair") fires too late.
   - The method must say that the two inner-panel descriptors (0003, 0486) and the parent are pre-registered. A lookup miss on any of them must raise before `begin`.
2. **`ExplicitParentIndex.claim_parent` returns `None` on any non-matching key** (`recovery.py:205-206`).
   - A reconstructed outer-panel descriptor that differed in any byte would fall through to `begin`. That creates a new nearby parent and leaves the original PENDING forever.
   - This conflicts with the method's own "first parent must be the exact allowlisted parent" rule.
   - The prototype or method should require that the first request through the wrapper is claimed, and refuse otherwise.
3. **Record 0002 would be duplicated on re-entry, and the method does not say what to do.** Re-running `work()` calls `self.put(..., 'outer-node-input', parent)` again at original-numeric.py:273 for the same node 0002 already records. Choose one of two options:
   - **Skip.** Verify that the canonical encoding of the reconstructed `parent` equals record 0002, and skip the put.
   - **Re-put.** Re-put it under a serial of 581 or higher, labelled as a recovery re-execution.

   Either is acceptable. Leaving it open is not.

## Prototype defects (`recovery.py`)

- `audit_prefix` checks only that the two witnesses share the global `context` and route. It does not check that each witness's `outerParent`, `arguments`, cell and `representedPoint` join to the pending parent. A manufactured contract with panels from a different node or cell would pass. I did that join by hand above, but the prototype should do it.
- `ExplicitCloneStore.__init__` does not check free space before the copy, and it does not fsync the directory. A failed copy leaves a partial destination that the exclusive `xb` then refuses to reuse. That is arguably correct preservation, but it should be stated.
- `prefix_preserved` is a method the caller must invoke. It does not check that `contracted_leaves` is still empty or that no original namespace changed beyond the prefix check.
- `import os` sits inside the `with` block in `__init__`. This is a style point only.

## Future concrete-build duties

- **Namespace allowlist.** `store.namespace()` inserts a new namespace on any digest mismatch. The build must assert that every A24 baseline and `formula-check:baseline` namespace equals the existing digest, with no new namespaces for those routes.
- **Serial counter.** Derive the serial from record names (max leading integer is 580). Do not hard-code it.
- **Fixed-init path.** `Evaluator.__init__` runs `CREATE TABLE contracted_leaves` and will fail on the clone. The build needs a fixed recovery-only initialization, with an AST delta certificate.
- **Typed geometry.** `self.line` and `self.quad` read `.slope`, `.intercept.a/.b/.d`. The typed geometry plans must be rebuilt from the saved `new-geometry-*` JSON and checked against it. Only cell 0 is joined against record 0003. The exact `rule receipt` path string (an absolute path in `settings.ruleReceipt`) must come from record 0 and not be recomputed.
- **Backend pinning.** Pin the mpmath version and backend for the re-applied pair arithmetic.
- **Execution provenance.** The provenance record should list which sequence ranges are old and which are new, including the reconstructed pair record and the parent's return. Otherwise a later lookup of the parent return would not show that it mixes old-code children with new-code continuation.
- **Failure handling.** `prefix_preserved` and the posthashes must run in a `finally` block. A failed clone is preserved and never resumed. A later attempt makes a new clone from the original.
- **Other int-keyed dicts.** I saw no other int-keyed dict reaching `put`, but later paths (fixed outer budget, completed-action) have not been exercised.

## Mathematical concerns (all inherited, none new)

- **A24/A48 pairing.** The A24 and A48 inner estimates are paired-order and not independent.
- **Error indicators.** The error values are empirical, not rigorous bounds.
- **Shared Gaussian identity.** The two full routes share one analytic Gaussian identity.
- **Scope.** J/direct is not a full packet action, and no leakage or current follows.

## Coverage

- **Read in full:**
  - `guide.md`
  - `method.md`
  - `recovery.py`
  - `request-index.py`
  - `evidence-store.py`
  - `original-numeric.py`
  - the stderr
  - `original-worker.py` lines 140-240
  - `tooling-tests.py` lines 17-90
- **Read in part:** `numeric-evidence-fix.py` (grep only: the import and the `panel_evidence(panels)` call at line 395), `contract.json` (grep and one excerpt), `export-index.json`, `receipts.json`, and `requests.json` (counts). Records 0001 and 0003 were truncated by the viewer, and I did not see the rest of record 0001 (including the end of its `exactFamilyContexts`) or record 0003 past its first 40,000 characters. Records 0002, 0003, 0486 and 1449 were checked by field-level grep for the joins above.
- **Not read:**
  - Returns 0485 and 1449 beyond the route and rule fields.
  - Individual node and transform rows.
  - `tooling-tests.log` and `initial-tooling-tests.log`, so I did not check the 33-test claim.
  - Geometry, tail, wing and census control files.
  - `original-worker.py` lines 1-139.
  - The `original-*` prepare, geometry and resume files.
  - Hashes: the database and geometry hashes are contract assertions I could not recompute.
- **Not done:** nothing was executed, and the A24/A48 pair values were not compared numerically.