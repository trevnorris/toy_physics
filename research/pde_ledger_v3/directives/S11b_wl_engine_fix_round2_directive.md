# S11b WL engine — fix round 2 (F-WL-3b/3c dead controls + chi5-tuned grazing)

## Authority and boundary
Repair `research/pde_ledger_v3/mathematica/S11b_interface_coupling_law_mathematica_audit.wl` IN PLACE
(committed baseline `a5186dce`). `CLAUDE.md` binds; `directives/S11b_SHARED_PHYSICS.md` is the sole physics
authority (⛔ do not re-open). ⛔ Add no expected value (rule 5). This is a NARROW fix: F-WL-1, F-WL-2, and
F-WL-3a are already genuine and correct and must stay behaviourally unchanged; so must every Scope object
(impedance/regimes, added mass, `S11B_ZPERM_SLICE`, transverse, breathing, longitudinal dispersion, roots).
Two defects remain, both "a residual/control that cannot fail" (the corollary-3 tautology).

## F-WL-3c — `GRAZING_MODE_CLASSIFICATION` classifies a chi5-TUNED block, and its "unrelated" control is a dead alias
Sites: `thresholdLeadingClassification` (~L1604-1615), `thresholdUnrelatedClassification` (~L1649).
- The classification is `classifyLeadingBlock[stratumBlock]` where `stratumBlock = block /. First[Solve[
  dispersion == 0, chi5]]` — i.e. `chi5` is SOLVED to force `Det = 0`, so the classified block is
  degenerate BY CONSTRUCTION and invariant to the block's actual form (a form change to the block is
  cancelled by `chi5` re-solving to keep `Det = 0`).
  Fix: classify the **ACTUAL, un-tuned sound-cone block** — `classifyLeadingBlock[block]` on the un-tuned
  `thresholdSoundConeBlocks` entries — so `RANK`/`NULLSPACE` report the assembled system's real leading-order
  structure (grazing shows up as a genuine rank drop only where the physics produces one, ⛔ not forced by
  tuning `chi5`). ⛔ Do not substitute a `Solve[Det==0, chi5]` stratum into the classified block. If the
  dispersion/stratum is still wanted, emit it under a DISTINCT field/tag as a computed object, ⛔ never as the
  classification's input.
- `thresholdUnrelatedClassification = thresholdLeadingClassification` is a DEAD alias ⇒ its `RANK_RESIDUALS`/
  `NULLITY_RESIDUALS` are `{0,0}` by construction and cannot fail. Fix: **DROP the unrelated-classification
  control entirely** (the built-in one-sided block-ablation control at ~L1633, which genuinely re-classifies
  an ablated block, is the real specificity control — keep it). ⛔ Do not leave an `A−A` residual emitted as
  a control.

## F-WL-3b — `CAUSALITY_CHECK` unrelated control is a dead A−A aggregation
Site: `unrelatedCausalityAggregation` (~L1237-1238, 1259-1262). The removed-record leg is genuine (presence
count and its residual bite). The "unrelated" leg re-aggregates the UNMODIFIED `kernelOrientationIdentities`
(identical to the baseline aggregation ⇒ `A−A ≡ 0`, cannot fail). Fix: **DROP the dead unrelated-aggregation
control** (the removed-record presence leg is the genuine control — keep it). ⛔ Do not emit an `A−A` residual.

## The three script clauses / corollary 3 (verbatim)
1. PRINT computed objects; ⛔ no prose conclusion. 2. PRINT the residual; ⛔ do not `assert` it zero. 3.
Interpretation → step record. ⛔ No tautological residual: an emitted residual/control must be able to be
nonzero for SOME input; where no genuine second route exists, emit the object and say so, ⛔ never an `A−A`.

## Run discipline (Mathematica)
One kernel at a time (2-seat licence); orchestrator writes the committed `.out` after review; demonstration
runs to scratch. Kill a demo kernel at 600 s no-new-output or RSS > 6 GB (~0.5 GB is normal). ⚠ Keep
demonstration runs FOCUSED — this repair touches only the grazing/causality tails; a truncated run that
reaches those tags suffices. (Long full runs have been spuriously killed here.)

## Acceptance — executable, no expected values (rule 5)
- **F-WL-3c**: the emitted `GRAZING_MODE_CLASSIFICATION.RANK`/`NULLSPACE` now reflect the ACTUAL sound-cone
  block — demonstrate that a FORM change to that block (a sign+off-diagonal flip) MOVES the emitted rank/
  nullspace (it no longer stays fixed via chi5-tuning), and that no `A−A` unrelated-classification residual
  remains. ⛔ Do not state the resulting rank value.
- **F-WL-3b**: no `A−A` unrelated aggregation residual remains; the removed-record presence control still
  fails on a dropped record.
- **No regression**: F-WL-1/2/3a and every Scope object stay behaviourally unchanged (byte-identical tags).

## Report (§13) — under 15 lines
The edits (tag:line), the FORM-ablation demonstration that the grazing classification now moves, confirmation
both `A−A` controls are gone (and which genuine control remains for each), and that no Scope object changed.
