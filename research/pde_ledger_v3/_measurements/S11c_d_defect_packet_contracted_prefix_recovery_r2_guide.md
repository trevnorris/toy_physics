# Review map: exact fixed-prefix state machine

Read `method.md`, `recovery.py`, `tooling-tests.py` and `tooling-tests.log` first. This is an explicit recovery-method/storage-prototype review, not a worker/build or numerical-result review. No scientific worker or READY gate for recovery exists.

`contract.json` binds the actual failed database, complete census, pending parent and both panel witnesses. `sqlite-export/` exposes every record as its exact original canonical JSON bytes. No scientific type was restored to make these exports. Original `failure/complete/numerical.sqlite` is also present.

Key raw records:

| Record path | Meaning |
|---|---|
| `sqlite-export/sqlite-records/0000.json` | Complete original contexts/rules/source/tail operands |
| `sqlite-export/sqlite-records/0001.json` | The single PENDING outer-panel full descriptor |
| `sqlite-export/sqlite-records/0002.json` | First outer-node caller input, actual represented point |
| `sqlite-export/sqlite-records/0003.json` | Complete A24-rule inner-panel input |
| `sqlite-export/sqlite-records/0485.json` | Its complete return |
| `sqlite-export/sqlite-records/0486.json` | Complete A48-rule inner-panel input, still own route A24 |
| `sqlite-export/sqlite-records/1449.json` | Its complete return and final chain head |
| `sqlite-export/requests.json` | All 435 full request descriptors/states/links |
| `sqlite-export/receipts.json`, `namespaces.json`, `export-index.json` | Exact chain, namespace and byte-origin receipts |

`original-numeric.py` is the unchanged evaluator whose `inner` evidence write failed. Inspect `Evaluator.__init__`, `ns`, `put`, `request`, `panel`, `inner` and `action`, plus `paired_indicator`, `encode` and `decode`. `original-worker.py` has the unchanged numerical tail; its prepare/geometry prefix must not execute again. The proposed worker path is specified in the method but is intentionally not implemented or cleared by this packet.

`numeric-evidence-fix.py` is a new exact copy with only the helper import and the pairedPanels evidence-boundary expression changed. The source tests compare its entire AST to the original after reversing those two changes. It has no database-attachment or pending-parent recovery path yet. It is not the future complete worker. The supplied final tooling log contains 54 distinct tests; the initial revision log is also preserved and explicitly includes duplicate inherited test executions. Neither test run restored native science.

`request-index.py` and `evidence-store.py` are unchanged shared libraries. Their normal refusal of pending/existing work remains. `recovery.py` now includes strict `join_frame`, clone admission, `ExplicitParentIndex` and `ExplicitRecoveryProtocol`. Manufactured tests exercise complete encoded callback routing with callback spies, required-child misses, exact first-node skip and committed-pair release, plus storage/prefix controls. Actual numerical wrapper integration and guarded runtime frame arithmetic are still future build duties.

`failure/` supplies the failed tree's complete scientific data and execution records. Historical reports/adjudications/authority and hook prose are omitted from this independent packet with exact opaque receipts in `excluded-historical-receipts.json`; original files are preserved and remain mandatory for the future runtime full-tree copy. `packet-index.json` maps every supplied path/hash/size. Some exact scientific files recur under older byte copies; no scientific array is shortened.

The four new geometry plans and their full input/return records are `failure/complete/new-geometry-{27,29}-{1,0}-{input,full,return}.json`. The four tail groups are `failure/complete/new-separate-analytic-tail-sums-{27,29}-{matching,zero}.json`; final guard records and all census/wing controls are in the same directory. The actual error is in `failure/defect_packet_contracted_continue_fix.stderr`, with context in `failure/complete/failure.json` and `numerical-journal-final-receipt.json`.

The main judgment is whether this explicit recovery preserves the reviewed numerical meaning while distinguishing complete child reuse from reconstruction of unpublished parent scalars. If it does not, explain the concrete stopping condition or required design change. A full action, numerical error pass or leakage cannot be inferred from the stored prefix. State what you actually read and any coverage gaps.

For readable full contexts, `views/record-0001/root.json` points to complete constituent values such as `context.exactFamilyContexts.entries.0.json`; its index carries JSON pointers and raw-source hashes. The same lossless view exists for records 0,2,3,485,486,1449. Each split object/array lists every child and its original key/index. Follow all entries needed to assess context rather than assuming a long viewer truncation is the complete record. Raw originals and receipts are unchanged.

The required sequence is PARENT -> NODE -> CHILD24 -> CHILD48 -> PAIR -> NEW. Lookups cannot use nullable fallback in the first four request states. No new namespace or request is available until the pair has actually committed. Record 2 is verified and skipped; serial increment is skipped too. The original mathematical callback remains unintegrated in the boundary-only numeric copy. This is still method plus prototype assessment, not full build clearance.
