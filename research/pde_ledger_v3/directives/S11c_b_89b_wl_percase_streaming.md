# S11c-b #89b (WL) — per-case STREAMING so the full run fits in memory (finish the build)

## 0 · Role and single deliverable

The engine `research/pde_ledger_v3/mathematica/S11c_b_brane_operator_mathematica_audit.wl` is **already
complete and correct** — the physics (operator un-freeze, tractable activation, per-summand Series) and every
tag are done and were verified out-of-band. The ONE defect is **memory**: the full run reproducibly OOMs the
30 GB machine — the WolframKernel peaks ~15–16 GB and the run is killed at the late equivalence controls. Make
the run **complete on this box** by holding **only one case's large operators in memory at a time**, changing
**nothing** about what any tag computes or its value. Re-run to a fully-regenerated `.out` with every tag.
This is a memory/architecture fix, **not physics**. ⛔ Do not commit; the orchestrator commits after review.

⛔ Every emitted payload must be value-identical to the pre-fix computation on the coverage it reports. ⛔ No
operator/kernel/residual/control MEANING changes.

## 1 · The measured cause (SUPPLIED — verified)

`mainModels` (built at ~L2146) stores, per case, the **un-reduced live operator**
`mainModels[key]["KERNEL_SOURCE_OPERATOR"] = operatorLive` (~L1390/1403) — a huge object (`LeafCount` ~50k+/row).
All **four** cases' `KERNEL_SOURCE_OPERATOR` (and `frozenModels`, ~L2429) stay resident from the operator emit
(~L2154) through the late controls (~L2531–2560), where each control additionally builds `liveOperator`,
`activateSpatialDivergences[liveOperator]`, `activateOperatorPerRowTopDownReference[liveOperator]`, and their
expanded difference. Four cases' un-reduced operators + those transients = the ~16 GB peak. Codex's earlier fix
freed the control *transients* but NOT the resident four-case `KERNEL_SOURCE_OPERATOR` baseline, so it still
OOM'd. `$HistoryLength=0` and per-case `ClearSystemCache[]` are already present.

## 2 · What to do — hold ONE case at a time (you choose the mechanism)

The invariant: **at no point are more than one case's un-reduced/large operators resident.** Reach it however
is cleanest — e.g.:

- **Case-outer streaming.** Compute one case's model, contribute its share to every tag's (small, reduced)
  emit-ready payload AND run every control for that case accumulating only the small results, then
  `Clear`/drop that case's `KERNEL_SOURCE_OPERATOR` and all large transients (`KeyDropFrom`/set the heavy
  fields to a placeholder) before the next case. Emit the assembled tag-keyed payloads at the end. The emitted
  `.out` format and per-tag content stay identical (still `TAG: <|case1->…, case4->…|>`); only the *build
  order* and *memory residency* change.
- **Or** keep the current tag-outer emit but, for the memory-heavy controls, iterate cases one at a time and
  free each case's un-reduced source operator immediately after that case's control result is captured — so the
  four-case un-reduced baseline is never simultaneously live during the controls. The reduced emitted operators
  (small) may stay.
- Keep the small reduced emitted operators; it is only the **un-reduced** `KERNEL_SOURCE_OPERATOR` and the
  control transients that must be streamed/freed.

Add at the end: `WriteString[First[$Output], "S11CB_MAX_MEMORY_USED_BYTES: " <> ToString[MaxMemoryUsed[]] <>
"\n"];` (a computed diagnostic — fine to emit).

## 3 · Allowed fallback if one case STILL will not fit (with a REQUIRED note)

If, after per-case streaming, a **single** case's heaviest equivalence control still peaks too high to finish:
you MAY compute that one control's full-operator equivalence residual (`FULL_OPERATOR_RESIDUAL_*`) on **one
representative case only** — it is a **case-independent code property** (the same activation/Series functions
run for every case) — and for the other cases emit an explicit marker
`"VALIDATED_ON_REPRESENTATIVE_CASE" -> <repr key>` INSTEAD of silently dropping the residual. ⛔ Do this ONLY
for the genuinely case-independent full-operator equivalence residual, and ONLY if per-case streaming alone
does not fit. State in an emitted note which controls (if any) you scoped this way and why. Every other
residual stays per-case. If you use this fallback, the full-fidelity all-cases in-band run is DEFERRED to a
bigger box — the orchestrator records that separately.

## 4 · Acceptance + report

- The run **completes**: `.out` regenerates with **all** ~41 `WL_S11CB_*` tags (operator/kernel/§5.E + the
  three equivalence controls), ends cleanly, and does NOT OOM. Emit `S11CB_MAX_MEMORY_USED_BYTES`.
- ⛔ Every payload value-identical to pre-fix on its reported coverage; ⛔ no physics/meaning change.
- Run it yourself: `wolframscript -file <this engine> > <its out.out> 2> <stderr>` (⛔ no `timeout` on the full
  run; the `.out` is a plain file, safe to overwrite). Watch memory; if a run is going to OOM, it is better to
  let it be killed and report than to leave a 15 GB kernel — but the point of the fix is that it should not.
- Report: how you streamed (case-outer / per-control free), the final `MaxMemoryUsed[]`, the total tag count,
  whether you used the §3 fallback (and on which controls), and confirm no residual's meaning changed.
