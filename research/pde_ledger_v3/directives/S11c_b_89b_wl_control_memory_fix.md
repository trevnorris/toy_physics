# S11c-b #89b (WL) — TARGETED memory fix: the equivalence controls OOM the machine (peak ~15 GB)

## 0 · Role and single deliverable

The engine `research/pde_ledger_v3/mathematica/S11c_b_brane_operator_mathematica_audit.wl` is **already
complete and correct** (it is your own prior build — do NOT re-derive the operator, the kernel, the
tractability activation, the Series-linearity fix, or any physics; those are done and verified to emit). The
ONLY defect is **memory**: the full run reproducibly OOMs the 30 GB machine — the WolframKernel balloons to
~15 GB and the run is killed — always at the **late equivalence controls**, right after
`WL_S11CB_CONTROL_BACKGROUND_JET_TOWER_OPERANDS` emits. Make the run **fit in memory and complete**, changing
**nothing** about what any tag computes or its numeric/structural value. Re-run the engine to completion so the
`.out` is fully regenerated with every tag. This is a build-completion / memory-hygiene fix, **not physics**.

⛔ Do not change any operator, kernel, residual, or control's MEANING. Every emitted payload must be
value-identical to what it would be without this fix (only computed with a lower memory peak).

## 1 · The measured cause (SUPPLIED — verified by direct timing/RSS)

The two controls below build several copies of the large live operator (`LeafCount` ~50k/row) for **all four
`mainModels` cases at once**, while `mainModels` + `frozenModels` (8 full operators) are already held globally:

- **`CONTROL_SURFACE_BACKGROUND_JET_DIFFERENCE_ATOMS`** (`surfaceHessianWitnessPayload`, ~L2500–2532):
  `AssociationMap` over all `Keys[mainModels]`, referencing both `mainModels[key]` and `frozenModels[key]`
  surface/face/divergence-source expressions.
- **`CONTROL_TRACTABLE_ACTIVATION_EQUIVALENCE`** (`tractabilityPayload`, ~L2538–2560): for each key it builds
  `liveOperator`, `tractableOperator = activateSpatialDivergences[liveOperator]`, and `referenceOperator =
  activateOperatorPerRowTopDownReference[liveOperator]`, then `mapDifference[…]` — three large operators plus
  their expanded difference, held simultaneously across the `AssociationMap`.

The kernel peaks at ~15 GB there; with other machine consumers that exhausts RAM and the run dies. Everything
BEFORE these controls (operator, kernel, term origins, freeze diagnostic, rep-invariance, independence,
Hessian-zero, jet-tower operands) emits fine at low memory.

## 2 · What to do (memory only — you choose the mechanism)

Reduce the peak so the run completes on ~20 GB of headroom. Keep every residual's meaning. Options you may
combine (your call — you are the builder):

- **Stream, don't accumulate.** Compute these controls one key (case) at a time and one surface/row at a time,
  emitting or accumulating only the small final payload, and `Clear[…]`/`ClearSystemCache[]` the large
  transient operators (`liveOperator`, `tractableOperator`, `referenceOperator`, per-surface expressions)
  **before** moving to the next key/surface. Do not hold all four cases' big operators live at once.
- **Do not hold both `mainModels` and `frozenModels` in full if you can recompute or stream** the per-surface
  difference; free each surface's big expression once its small atom-set difference is extracted.
- **`MemoryConstrained` is a guard, not a solution** — if you wrap a sub-computation in it, a hit must be
  surfaced as an explicit emitted marker, ⛔ never silently drop a residual.
- ⚠ If — and only if — the full-operator equivalence residual is genuinely a **case-independent code
  property** (the same activation/Series functions run for every case), you MAY compute the heavy
  `FULL_OPERATOR_RESIDUAL_*` on **one representative case** and, for the other cases, emit a cheap marker
  stating it was validated on the representative case. If a residual could differ per case, keep it per case.
  State in a one-line emitted note which choice you made and why. ⛔ Do not silently reduce coverage.

Also set/confirm `$HistoryLength = 0` is in force (it is, lines 1/3) and that the per-case emit loop for the
operator/kernel already `ClearSystemCache[]`s between cases (it does) — extend the same discipline into these
controls.

## 3 · Acceptance

- The run **completes**: the `.out` regenerates with **all** ~41 `WL_S11CB_*` tags present (including
  `CONTROL_SURFACE_BACKGROUND_JET_DIFFERENCE_ATOMS`, `CONTROL_TRACTABLE_ACTIVATION_EQUIVALENCE`,
  `CONTROL_SERIES_LINEARITY_EQUIVALENCE`, and the §5.E `ENERGY_BASIS_NEW_INVARIANT_DIMENSION_DERIVATION`), and
  the file ends cleanly.
- Peak `MaxMemoryUsed[]` stays well under the machine limit (target: the run does not OOM). You may print
  `MaxMemoryUsed[]` at the end as a diagnostic (a computed number, fine to emit).
- ⛔ Every equivalence/residual payload is value-identical to the pre-fix computation on whatever coverage it
  reports. ⛔ No physics, operator, kernel, or control MEANING changed.

## 4 · Run + report

Run the engine to completion via `wolframscript -file <this engine> > <its out.out> 2> <stderr>` (⛔ no
`timeout` on the full run; it is silence/RSS-bounded). The `.out` is a plain file now (safe to overwrite).
Report: which controls you restructured and how; the final `MaxMemoryUsed[]`; the total tag count in the
regenerated `.out`; and confirm no residual's meaning changed. ⛔ Do not commit; the orchestrator commits after
review.
