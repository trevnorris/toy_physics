---
title: Stage1 Verdict
type: source
status: current
sources:
- software/stage1_solver/STAGE1_VERDICT.md
- software/stage1_solver/frozen/m1c/834835c999efd572e7aba46e04831a201ef51a3a7c88820459964b44b303d5e8/forbidden_R_norm_validation_report.json
- software/stage1_solver/frozen/m1c/834835c999efd572e7aba46e04831a201ef51a3a7c88820459964b44b303d5e8/freeze_hash.txt
- software/stage1_solver/frozen/m1c/834835c999efd572e7aba46e04831a201ef51a3a7c88820459964b44b303d5e8/freeze_sheet.json
- software/stage1_solver/frozen/m1c/834835c999efd572e7aba46e04831a201ef51a3a7c88820459964b44b303d5e8/m1c_physical_derived_diagnostics.json
- software/stage1_solver/frozen/m1c/834835c999efd572e7aba46e04831a201ef51a3a7c88820459964b44b303d5e8/m1c_physical_derived_report.md
- software/stage1_solver/frozen/m1c/834835c999efd572e7aba46e04831a201ef51a3a7c88820459964b44b303d5e8/m1c_physical_summary.json
- software/stage1_solver/frozen/m1c/834835c999efd572e7aba46e04831a201ef51a3a7c88820459964b44b303d5e8/m1c_v2_22b_physical_frozen_packet.json
- software/stage1_solver/frozen/m1c/834835c999efd572e7aba46e04831a201ef51a3a7c88820459964b44b303d5e8/repro_snapshot_run1.json
- software/stage1_solver/frozen/m1c/834835c999efd572e7aba46e04831a201ef51a3a7c88820459964b44b303d5e8/repro_snapshot_run2.json
- software/stage1_solver/frozen/m1c/834835c999efd572e7aba46e04831a201ef51a3a7c88820459964b44b303d5e8/reproducibility_report.json
- software/stage1_solver/frozen/m1c/834835c999efd572e7aba46e04831a201ef51a3a7c88820459964b44b303d5e8/v2_21_manifest.json
- software/stage1_solver/frozen/m1c/834835c999efd572e7aba46e04831a201ef51a3a7c88820459964b44b303d5e8/v2_22a_profile_manifest.json
- software/stage1_solver/frozen/m1c/834835c999efd572e7aba46e04831a201ef51a3a7c88820459964b44b303d5e8/v2_22c_observable_packet.json
- software/stage1_solver/frozen/m1c/834835c999efd572e7aba46e04831a201ef51a3a7c88820459964b44b303d5e8/v2_22c_pipeline_report.json
- software/stage1_solver/frozen/m1c/834835c999efd572e7aba46e04831a201ef51a3a7c88820459964b44b303d5e8/v2_22c_tolerance_budget.json
- software/stage1_solver/frozen/m1c/834835c999efd572e7aba46e04831a201ef51a3a7c88820459964b44b303d5e8/wp1_background_10x8.json
- software/stage1_solver/frozen/m1c/834835c999efd572e7aba46e04831a201ef51a3a7c88820459964b44b303d5e8/wp1_background_6x4.json
- software/stage1_solver/frozen/m1c/834835c999efd572e7aba46e04831a201ef51a3a7c88820459964b44b303d5e8/wp1_background_8x6.json
- software/stage1_solver/frozen/m1c/834835c999efd572e7aba46e04831a201ef51a3a7c88820459964b44b303d5e8/wp1_background_resolution.json
- software/stage1_solver/frozen/m1c/834835c999efd572e7aba46e04831a201ef51a3a7c88820459964b44b303d5e8/wp1_solve_manifest.json
last_updated: '2026-08-25'
---

## Purpose and scope

### source-stage1-verdict--target-blind-effective-closure-test — Target-blind effective-closure test

Status: `lifecycle=current` · `evidence=provisional` · `memory_review=ai_draft`

The result source describes a target-blind, dual-engine Stage-1 test of one frozen moving-throat effective-closure branch. It says the workflow solves a self-consistent WP1 background, derives a BdG-plus-wall and vector-harmonic Maxwell outgoing transfer, passes that transfer through the existing Python V2 chain, and judges a GR-normalized residual. This scope is conditional on the frozen decision-07 posits.

Sources:

- `software/stage1_solver/STAGE1_VERDICT.md` — heading `## The test`

## Source-unit map

- Result and interpretation: `software/stage1_solver/STAGE1_VERDICT.md`
- Readable diagnostic measurement report: `software/stage1_solver/frozen/m1c/834835c999efd572e7aba46e04831a201ef51a3a7c88820459964b44b303d5e8/m1c_physical_derived_report.md`
- Authoritative V2 handoff packet and supporting validation, reproducibility, summary, manifest, background, and tolerance artifacts: identity-only members listed in `source_unit.members`

## Key statements

### source-stage1-verdict--robust-miss — Provisional effective-closure miss

Status: `lifecycle=current` · `evidence=provisional` · `memory_review=ai_draft`

The verdict reports a robust magnitude miss for this branch: `R_norm = -10.7999993`, with outgoing transfer about seven orders of magnitude below the stated GR target. This remains provisional because the prepared unit supplies neither the exact recorded invocation nor readable contents of the authoritative V2 packet.

Sources:

- `software/stage1_solver/STAGE1_VERDICT.md` — heading `## The result — a robust MISS`

### source-stage1-verdict--diagnostic-preview — Direct-formula values are diagnostic only

Status: `lifecycle=current` · `evidence=provisional` · `memory_review=ai_draft`

The readable derived report gives post-freeze diagnostic values including `D0 = 4.055722430990979`, `R_norm = -10.799999340551249`, `R_pole = -0.8310527171393044`, `P2 = -2.0665687783479044e-7`, and `P4 = 6.419188997684018e-8`. The report explicitly labels these as physical-run diagnostics and says V2 remains authoritative; the authoritative packet is not semantically readable here.

Sources:

- `software/stage1_solver/frozen/m1c/834835c999efd572e7aba46e04831a201ef51a3a7c88820459964b44b303d5e8/m1c_physical_derived_report.md` — heading `## Direct Formula Preview`
- `software/stage1_solver/frozen/m1c/834835c999efd572e7aba46e04831a201ef51a3a7c88820459964b44b303d5e8/m1c_v2_22b_physical_frozen_packet.json` — `anchor-unavailable`

### source-stage1-verdict--magnitude-and-resolution-guards — Reported magnitude and resolution guards

Status: `lifecycle=current` · `evidence=provisional` · `memory_review=ai_draft`

The verdict characterizes the discrepancy as a magnitude gap rather than a sign flip: it says `N0` is nonnegative, the spatial-gauge sign flip changes `N0` by zero, and only transfer `+10.8` would make `R_norm = 0`, whereas the reported transfer is about `6.6e-7`. It also claims the miss is not an `R_norm`-path resolution artifact: the large relative floor concerns near-null `min_density`, while the stated transfer-mesh, interpolation, and absolute-residual budgets remain far smaller than the `10.8` gap. These are source-reported guards, not independently verified measurements in this capsule.

Sources:

- `software/stage1_solver/STAGE1_VERDICT.md` — heading `## The result — a robust MISS`

### source-stage1-verdict--interpretation-scope — Branch result does not decide the full framework

Status: `lifecycle=current` · `evidence=provisional` · `memory_review=ai_draft`

The verdict limits the reported miss to one pre-registered branch of a deliberately incomplete effective-closure model that omits the matter/gauge-to-wall return coupling. It does not present the result as falsifying the overall framework, and it rejects fitting free-choice values after seeing the target.

Sources:

- `software/stage1_solver/STAGE1_VERDICT.md` — heading `## Interpretation — what it means (and does NOT mean)`

## Computed evidence represented by the source

### source-stage1-verdict--available-evidence-chain — Prepared evidence boundary

Status: `lifecycle=current` · `evidence=provisional` · `memory_review=ai_draft`

The readable verdict and derived report supply a prose workflow description, diagnostic operands and residuals, and an interpretation, but the exact recorded invocation is not supplied. The authoritative V2 packet, forbidden-residual validation report, and reproducibility report are identity-only members. Their literal contents, values, tags, and success states were not readable and are not asserted here.

Sources:

- `software/stage1_solver/STAGE1_VERDICT.md` — heading `## The test`
- `software/stage1_solver/frozen/m1c/834835c999efd572e7aba46e04831a201ef51a3a7c88820459964b44b303d5e8/m1c_physical_derived_report.md` — heading `## Direct Formula Preview`
- `software/stage1_solver/frozen/m1c/834835c999efd572e7aba46e04831a201ef51a3a7c88820459964b44b303d5e8/m1c_v2_22b_physical_frozen_packet.json` — `anchor-unavailable`
- `software/stage1_solver/frozen/m1c/834835c999efd572e7aba46e04831a201ef51a3a7c88820459964b44b303d5e8/forbidden_R_norm_validation_report.json` — `anchor-unavailable`
- `software/stage1_solver/frozen/m1c/834835c999efd572e7aba46e04831a201ef51a3a7c88820459964b44b303d5e8/reproducibility_report.json` — `anchor-unavailable`

## Assumptions, exclusions, and open questions

### source-stage1-verdict--full-path-a-open — Full Path-A model remains open

Status: `lifecycle=current` · `evidence=open` · `memory_review=ai_draft`

Whether the full self-consistent Path-A model reproduces GR remains explicitly open. The verdict defers derivation of the matter/gauge-to-wall return coupling `S_η^(ψ,A)` to a new, separately pre-registered subprogram and does not promise that it closes the reported gap.

Sources:

- `software/stage1_solver/STAGE1_VERDICT.md` — heading `## Open after Stage-1 (deferred)`

### source-stage1-verdict--diagnostic-resolution-limits — Resolution and interpolation limits

Status: `lifecycle=current` · `evidence=provisional` · `memory_review=ai_draft`

The derived report labels the physical frozen decision-07 result as conditional on the frozen target-blind posits and modest CPU resolution. Separately, its order-2 interpolation comparison is a diagnostic re-extraction on the same frozen background rather than a background-resolution study; the reported maximum relative absolute interpolation-method delta is `0.009519365305710407`.

Sources:

- `software/stage1_solver/frozen/m1c/834835c999efd572e7aba46e04831a201ef51a3a7c88820459964b44b303d5e8/m1c_physical_derived_report.md` — exact text `LOUD LABEL: PHYSICAL frozen decision-07 effective-closure branch. The freeze hash was written before the WP1 solve; the result is conditional on the frozen target-blind posits and modest CPU resolution.`
- `software/stage1_solver/frozen/m1c/834835c999efd572e7aba46e04831a201ef51a3a7c88820459964b44b303d5e8/m1c_physical_derived_report.md` — heading `## Interpolation-Method Sensitivity`

## Revision and supersession relationships

## Related topics and scripts

- `memory/sources/software/stage1-solver.md`