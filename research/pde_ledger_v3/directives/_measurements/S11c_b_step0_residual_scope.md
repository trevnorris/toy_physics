# STEP 0 — cross-engine residual memory/time scope on the 30 GB box (2026-09-03)

Question (resume-prompt STEP 0): does regenerating the objects the cross-engine `row_residual` needs — a
single-case un-frozen WL SLAB operator **and** the primary COUPLING `FINAL_KERNEL` — fit on this 30 GB box, or
is it ≥64 GB-gated like the full in-band run? The 0.9 GB figure from 2026-09-02 was the U-row operator ONLY; the
kernel footprint was unmeasured (Codex compact-prep correction).

## What was measured
A guarded single-case WL run for `EULERIAN / LAB_HELD / RHO4_CONSTANT` that builds ONLY `evaluatedModel[…]`
(the operator) and `extractCouplingData[…]["FINAL_KERNEL"]` (the primary kernel) — i.e. exactly the objects the
residual compares, and NONE of the per-case loop's tower-depth control variants (`operatorTruncated/Extended`,
`kernelTruncated/Extended`) that dominate the full run's ~16 GB/case peak.

- Probe: `scratchpad/step0_wl_op_kernel.wl` (the proven StringTake-between-markers loader, FULL 40-invariant basis;
  emits `MaxMemoryUsed[]` per phase; exports the operator U-row + `FINAL_KERNEL`).
- Runner (memory watchdog): `scratchpad/step0_run.sh` — samples WolframKernel RSS + system MemAvailable every 4 s,
  records peak RSS and min MemAvailable, kills by pid if MemAvailable < 2500 MB.
- Command: `bash scratchpad/step0_run.sh` (backgrounded; no timeout).
- Evidence: `~/.s11_build/S11c_b_step0_scope/{run.log, mem_samples.log, wl_op_urow.txt, wl_final_kernel.txt}`.

## Literal result (`~/.s11_build/S11c_b_step0_scope/run.log`)
```
PHASE definitions-loaded basisReps=40
PHASE model-evaluated memInUse_GB=0.192535 maxUsed_GB=0.217158
OPERATOR_ROWS uRow_leaf=24601 thickness_leaf=13597 mass_leaf=5111 muTheta_leaf=2745
PHASE operator-rows-extracted memInUse_GB=0.192538 maxUsed_GB=0.217158
FINAL_KERNEL kernel_leaf=250223
PHASE final-kernel-built memInUse_GB=0.926238 maxUsed_GB=0.987058
PHASE exported memInUse_GB=0.926343 maxUsed_GB=0.987058
DONE
RUN_END rc=0 peak_kernel_RSS_GB=7.95 min_MemAvailable_GB=14.94 2026-09-03 01:17:06
```
Wall time (samples): 00:50:43 → 01:17:02 ≈ **26.3 min** for one case (operator + FINAL_KERNEL).

## Finding
- **In-kernel high-water (`MaxMemoryUsed[]`): 0.99 GB.** **External kernel RSS peak: 7.95 GB.** **min MemAvailable:
  14.94 GB** — never near the floor. The single-case PRIMARY operator + kernel fits the 30 GB box with ~2× headroom
  on RSS, ~30× on the in-kernel measure.
- The ~16 GB/case in `DEFERRED_HEAVY_RUNS.md` is the per-case loop's TOWER-DEPTH CONTROL VARIANTS
  (`operatorTruncated/Extended`, `kernelTruncated/Extended`), built unconditionally at `…mathematica_audit.wl:2204-2231`
  and NOT needed by the residual. ⇒ the cross-engine residual is NOT ≥64 GB-gated. `DEFERRED_HEAVY_RUNS.md` is
  correct for the full 4-case in-band `.out` regen (with all tower/heavy controls); it does not bound the residual.
- ⚠ OPERATIONAL DISCOVERY: the engine's NORMAL emit path always builds the tower variants (no env gate skips them,
  only `S11CB_SKIP_HEAVY_CONTROLS` gates the two equivalence controls). So a comparator-PARSEABLE single-case `.out`
  at this ~8 GB footprint requires a **residual-mode single-case emit** — restrict to one case and call the engine's
  OWN `emitShared`/`modelRecord`/`kernelRecord` on the primary objects, skipping the tower-variant build. Same
  objects, same emit functions; it drops only the controls the residual does not compare. This is an added piece of
  the integration pass (both engines).
- PY side: the committed PY `.out` was generated PRIMARIES_ONLY (completed once), and the #90 build emitted one
  folded+#90 case (564 s smoke). A single-case folded PY emit is known-tractable; a fresh single-case PY `.out`
  regen is the analogous residual-mode emit. (To be reconfirmed when regenerated.)

## Verdict
The cross-engine single-case residual (validating the constraint-fold's U/E_W rows and #90's face+response kernel)
is DOABLE on this 30 GB box. The ≥64 GB box stays required only for the full 4-case in-band `.out` regen with tower
+ heavy controls (belt-and-suspenders, already proven out-of-band per `DEFERRED_HEAVY_RUNS.md`).
