# stage023 step (h) — orchestrator ablation, re-run 2026-07-28

⚠ **Supplementary evidence.** An orchestrator ablation does not replace the fresh-agent adversarial leg,
which owns the verdict (`docs/development_pipeline.md`, Roles/Ablations row). What it establishes is the
per-target behaviour recorded below — each row backed by its own capture at run time, and by its
`results.tsv` row in this commit (see Files).

## Method

Include-list: `include_list.tsv`, 51 rows, orchestrator-owned (`make_lists.py` derives the rows
mechanically from the two source blocks; the *choice* of targets — every declaration, every emitted
record, one at a time — is the orchestrator's). Driver: `ablate.py`. Per mutant, in this order:

1. restore the stage source from `PRISTINE.py` (`cp`, never `git checkout`);
2. **delete the emitted sidecar** `…_sympy_audit.dimensions.txt`;
3. `rm -rf scripts/__pycache__`;
4. apply the mutation **on disk at the real path** (the sidecar's `source_sha256` freshness binding and
   the module pin only stay honest if the real file at the real path is what changes);
5. run the stage (`timeout 600`), capture stdout to `captures/<axis>.<key>.stdout`;
6. run `compare_dimension_artifacts.py 023`, capture to `captures/<axis>.<key>.cmp`;
7. append one flushed TSV row to `results.tsv`.

Restore at the end was proved by digest rather than asserted, in `run.log`'s last two lines:
`stage_sha256=10f125307d0a3a0bd09b19d32e0f4f66ce6512a4fc417f58ac2fca11c670edc4`,
`sidecar_sha256=cba3a9de579073d60fec8d11268cb5a626017349591290f0d5d274b7f637f096`, both equal to
pristine; `git hash-object` = `6ecc9fd00a1879c62fc1710d2f2c5c41ce741324`. ⚠ `run.log` does not survive
the commit (see Files), so **for a reader of the commit these three digests are a transcription, not a
check they can re-run**; what is independently checkable is the committed script's own pair, next.

⚠ **That digest is the ablation-time source, not the committed one — establish the difference before
reading the tables.** The committed stage script is
`sha256=46f90b7c937d7d229727821db0ba5c080a93240ff506afad41ba156de7dde5c1`
(`git hash-object` `26995ae877f092ba447cac7edd8300bfbe776439`); step-(h) remediation moved it on after
this run. Measured by `diff` of the run-time pristine copy (`PRISTINE.py`, likewise not retained)
against the committed file, the **whole** difference is one statement-order change in `main()`: `emit_dimension_sidecar(...)` now precedes
`print_verdict_labels()` instead of following it (`_scratch/stage023_h/REPORT_REMEDIATE_H.md` R2). No
other byte differs — round-1 R1 also edited `fmt` and `dimension_records`, and round-2 R4 reverted both,
so they cancel.

**That difference cannot affect any row of either table.** The two reordered statements are the last two
in `main()`, after every check and every print these tables record, and sidecar emission itself prints
nothing — so neither the stdout the captures hold, nor the sidecar bytes the comparator reads, nor the
tally, verdict or exit code changes with the order. The order becomes observable only when the emission
**raises**, which is the injected fault R2 was written for; neither axis injects one. The 16 caught rows
abort before `main()`'s tail under either order, and every completing row runs both statements. The
tables therefore stand as run.

**Two defects in the previous run are closed, and both were closed by construction rather than by
watchfulness.** (i) The sidecar is deleted before every iteration, so no verdict is decided by a residual
artifact — in the previous run 16 of 22 comparator verdicts were `COMPARISON_SKIPPED` freshness failures
rather than value comparisons. (ii) Every mutant kept its own `.cmp` and `.stdout` capture — 51 of each,
none overwritten — so no row's verdict was read off another's. ⚠ **What survives into the commit is the
tables, not the captures:** all 51 rows carry their exact mutation in `include_list.tsv`
(`old_text`/`new_text`) and their recorded `stage_exit`, tally, `first_fail`, `sidecar_written`,
`record_moved`, `cmp_exit`, `cmp_status` and `mismatch_names` in `results.tsv` — enough to re-derive
every claim below, and enough to re-run any row from scratch. What a reader **cannot** do from the commit
is check those recorded values against the raw capture each was transcribed from (see Files).
**A third defect is closed by choosing a different decoy:** axis 1 previously used `Dim(7, 0, 0)`, which
is the value the stage's own free-carrier control substitutes at `.py:509`, so during the `q_free` row
that control was silently a no-op. The decoy here is `Dim(13, -7, 5)`, which occurs nowhere in the source.

## Axis 1 — DECLARATION: replace one `SOURCED_DIMS` value with the decoy (22 rows)

**Every one of the 22 is detected — but by two different instruments, and the split is the finding.**

- **16 declarations are caught by the stage itself** (exit 1). No sidecar is emitted, so the comparator
  reports the artifact absent (exit 2) rather than a value mismatch: a failing producer leaves the
  committed sidecar to go stale, which is the freshness binding's job, not the comparator's.
- **6 declarations leave the stage green at 111 PASS** — `R0`, `R1`, `eta_null`, `gain0`, `gain1`,
  `q_free` — and the **comparator is the sole detector**, naming exactly the mutated record. `q_free` is
  a first measurement rather than a repeat, for the collision reason above.

⛔ **In none of the 16 caught rows does the dimensional gate's own assertion fire first.** `base_verdict`
ranks `dimensional` third in the verdict ladder (`.py:673-680`), so a corrupted declaration flips the
gate verdict before `run_dimensional_gate()` is ever reached. The assertion that actually fires is
`selector control reaches CROSS_L_RESIDUAL_PREDICTION` (11 rows: `a`, `c_s`, `omega`, `M0`, `D1`, `D0`,
`Omega_U`, `Omega_W`, `R_mix`, `g_U`, `g_W`) or `baseline reaches earned
FAIL_UNDERDETERMINED_NOT_PREDICTIVE` (5 rows: `K0c`, `K_eta`, `T_Omega`, `Z0_ret`, `Z1_ret`). The
detection is real; what it is *not* is the stage's dimensional assertions doing the detecting.

## Axis 2 — BINDING: repoint one emitted record at a different-valued source (29 rows)

`dimension_records()` is an explicit per-name dict literal (`.py:539+`) whose return value is passed to
the sidecar emitter and to nothing else — **no assertion and no verdict reads it** — so every row runs to
completion.

**29 of 29: stage exit 0, 111 PASS, `OVERALL PASS` — and the comparator names exactly the mis-bound
record, one mismatch, no collateral.** The decoy is `SOURCED_DIMS[a]` (`SOURCED_DIMS[c_s]` for the `a`
record itself); the `record_moved` column is computed per row against the pristine sidecar, so
"the decoy really is distinct from the true value" is measured for all 29 rather than assumed.

## What neither axis reaches

25 of the 29 records share their exact triple with a sibling (only `a`, `c_s`, `D0`, `R_mix` are unique),
so a mis-binding **within** a triple-class emits bytes identical to a correct one. Both engines and the
comparator are blind to it by construction; only a source read reaches it. That is the fidelity leg's
territory, not this instrument's, and it is why 29/29 here is a statement about *cross-class* bindings.

## Files

**Committed in this directory — the whole of what a reader can audit:** `ABLATION_SUMMARY.md` (this
file) · `include_list.tsv` (51 rows) · `results.tsv` (header + 51 rows).

⚠ **NOT retained — told about, not auditable.** The run's working set lived in
`research/pde_ledger_v2/_scratch/stage023_h/`, which `.gitignore:96` (`research/**/_scratch/`) excludes,
so none of it survives the commit: `captures/` (51 `.cmp` + 51 `.stdout` — the raw evidence each
`results.tsv` row was transcribed from; Method (ii) describes the run, not this commit) · `run.log`
(source of the restore digests quoted above) · `ablate.py` · `make_lists.py` · `PRISTINE.py` ·
`PRISTINE.sidecar` · `base_stdout.txt` · `base_cmp.txt`. From the commit alone, every per-row claim here
is checkable **against `results.tsv`** and against no capture.
