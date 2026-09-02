# S11c-b #89b (WL) — complete the run by deferring the two heavy in-band controls (PY-precedent fallback)

## 0 · Role and single deliverable

The engine `research/pde_ledger_v3/mathematica/S11c_b_brane_operator_mathematica_audit.wl` is complete and
correct (physics done; per-case streaming already added). It STILL OOMs the 30 GB box: even one case's two
heaviest equivalence controls peak ~14–15 GB. Per the PY-side precedent (`S11CB_PRIMARIES_ONLY` skips
build-heavy in-band controls, proven out-of-band), make the run **complete** by putting those two controls
behind a skip switch and running WITH the switch set. Everything else — the operator, kernel, μ_θ, term
origins, §5.E, and the cheap controls — must emit. Change **nothing** about what any emitted tag computes.
⛔ Do not commit; the orchestrator commits after review.

## 1 · Exactly what to guard

Wrap ONLY these two controls (and their heavy upstream builds) in an env switch
`S11CB_SKIP_HEAVY_CONTROLS` (skip when the variable is set):
- `CONTROL_SURFACE_BACKGROUND_JET_DIFFERENCE_ATOMS` (~L2594)
- `CONTROL_TRACTABLE_ACTIVATION_EQUIVALENCE` (~L2598)

When skipped, **still emit each tag** with a deferred marker instead of the heavy payload, e.g.
`emitShared["CONTROL_TRACTABLE_ACTIVATION_EQUIVALENCE", <|"DEFERRED_HEAVY_CONTROL" -> True, "REASON" ->
"OOM on 30 GB box; full in-band run deferred to a bigger box (DEFERRED_HEAVY_RUNS.md); verified OUT-OF-BAND by
the #89b build legs + tractability decision legs", "RUN_FULL_WITH" -> "unset S11CB_SKIP_HEAVY_CONTROLS"|>];` —
so the tag is present and self-documenting, ⛔ never silently absent. ⛔ Do NOT guard any other control
(Hessian-zero, re-freeze regression, independence, rep-invariance, jet-tower operands, Series-linearity all
stay live) and ⛔ do NOT touch the operator/kernel/μ_θ/§5.E emits.

Keep the `S11CB_MAX_MEMORY_USED_BYTES` diagnostic at the end.

⚠ **PLACEMENT — the heavy builds live INSIDE dual-purpose `Do` loops, so you must gate the tractability /
Hessian-witness statements ONLY, ⛔ never a whole loop iteration** (skipping the iteration would drop core
emits):
- `tractabilityPayload` (→ `CONTROL_TRACTABLE_ACTIVATION_EQUIVALENCE`) is assembled inside the SAME `Do`
  at ~L2179 that builds the core `mainModels` / `mainKernels` / `mainKernelOrigins` /
  `mainKernelActivationPostconditions`. Gate ONLY the `tractableOperator` / `referenceOperator` /
  `tractabilityCase` build and its `AssociateTo[tractabilityPayload, …]`; ⛔ KEEP every
  `mainModels` / `mainKernels` / origins / activation-postcondition assignment LIVE and unguarded.
- `surfaceHessianWitnessPayload` (→ `CONTROL_SURFACE_BACKGROUND_JET_DIFFERENCE_ATOMS`) is assembled inside
  the `Do` at ~L2497 that ALSO streams the must-keep `CONTROL_OPERATOR_REFREEZE_REGRESSION_OPERAND`
  (its `beginAssociationEmission` at ~L2491, `endAssociationEmission` at ~L2590). Gate ONLY the
  surface-Hessian-witness statements and their `AssociateTo[surfaceHessianWitnessPayload, …]`; ⛔ KEEP the
  refreeze-regression streaming (`beginAssociationEmission` / per-case emit / `endAssociationEmission`) LIVE.
- When gated, leave the payload association as (or set it to) the `DEFERRED_HEAVY_CONTROL` marker so the
  downstream `emitShared` at ~L2594 / ~L2598 still emits a self-documenting tag.

## 2 · Do NOT run the engine — parse-check only, then STOP

⛔ Do NOT run the full engine and ⛔ do NOT generate the `.out` — the orchestrator runs it via
`wolframscript` (an unrelated background session in this repo kills `codex` processes on a timer, which
would abort a long in-process run; a `wolframscript` launched by the orchestrator is immune). Your job is
the edit ONLY.

After making the guard edit, do only a fast SYNTAX check that does NOT execute any cell:
```
cd /var/projects/toy_physics
wolframscript -code 'r = ToExpression[Import["research/pde_ledger_v3/mathematica/S11c_b_brane_operator_mathematica_audit.wl","Text"], InputForm, Hold]; Print[If[Head[r]===Hold && FreeQ[r, $Failed], "PARSE_OK", "PARSE_FAIL"]]'
```
It must print `PARSE_OK` (this parses the whole file but evaluates nothing). Then STOP.

## 3 · Acceptance + report (edit only)

- Both controls' heavy BUILDS **and** their `emitShared` are gated by `S11CB_SKIP_HEAVY_CONTROLS`; when the
  variable is set each still emits its `DEFERRED_HEAVY_CONTROL` marker (⛔ never silently absent). When unset,
  the engine is byte-behaviour-identical to before (full heavy payloads).
- ⛔ Every OTHER build/emit is untouched — operator, kernel, μ_θ, §5.E, refreeze regression, Hessian-zero,
  independence, rep-invariance, jet-tower operands, Series-linearity, and `S11CB_MAX_MEMORY_USED_BYTES`.
  ⛔ No operator/kernel/physics change; ⛔ no whole-loop-iteration skip (per §1 PLACEMENT).
- The syntax check prints `PARSE_OK`.
- Report: the unified `git diff` of the engine, and confirm exactly which two controls are gated and that
  no other emit's build was moved or removed. ⛔ Do not commit.
