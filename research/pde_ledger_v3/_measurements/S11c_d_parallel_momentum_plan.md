# S11c-d isolated numerical workers

The user authorizes four independent numerical workers in one supervised job.
This replaces serial scheduling for the current source/profile quadrature;
retain one BLAS/OpenMP thread per worker and keep unrelated heavy CAS jobs out.
The physical engine and existing source/profile checker stay byte-identical.

1. Each worker owns one (test field, quadrature setting), its logs and checkpoint
   directory. Preserve the native node sequence and every within-integral floating
   addition. The coordinator waits on process sentinels and merges completed
   operands in the original field/setting order, irrespective of completion order.
2. Resume an interrupted serial grid from its last verified cumulative values,
   mutated values, mass, node/batch counts and exact source/field/domain joins.
   Insert only the saved accumulator seed before the unchanged native group loop;
   prove the remainder of that function is AST-identical. Skip the corresponding
   original batches without repeating their integrals. Replay the literal native
   source-frequency prefix on those nodes to reconstruct complete frequency counts.
   Retain the legacy combined workspace estimate as a conservative envelope for
   each split estimate and state this scope. It is distinct from measured RSS.
3. Before handoff, verify resume against an uninterrupted bounded serial prefix,
   including exact numerical arrays, mass, node and source-frequency counts.
   Test four-process scheduling on actual operands. For end-to-end aggregation,
   compute small, explicitly underresolved full quadratures, changing only numeric
   rule orders in the test constructor. Compare independent serial construction
   with parallel construction, all native terms/actions, emitted payloads and
   metadata; no accuracy/convergence conclusion follows from these focused inputs.
4. Preserve all original run files. Disarm only its owned completion watcher,
   capture and hash its saved operands, interrupt only its verified child process,
   and record the intentional handoff separately from physics/runtime failures.
   Re-read stable inventories after it exits. Completed grids can be resumed at
   their final accumulator; no accepted numerical work is discarded.
5. Launch at most four workers under the existing single-job supervisor. Each
   has a 2 GiB address-space ceiling and one native numerical thread. A worker
   error stops peers after preserving their checkpoints and wakes the assistant.
   Each worker saves result operands before guards; no failed result is accepted.
6. Replay complete worker results through the unchanged serial constructor's
   comparisons and guard logic, then its full metadata/emission validator. Keep
   structural manifests lossless. Verify source, worker, partial and pre/post
   emission packet hashes. Publish only validated output through DataLad/git-annex
   and commit each implementation/launch/publication checkpoint.

All work stays under repository scratch. Use the same silent local completion/
error watcher, never model polling or recurring checks. Keep the approved inputs,
full symbolic parameter dependence and retained solver/export contract. This
changes execution scheduling only; no new physical result or limit is assumed.
