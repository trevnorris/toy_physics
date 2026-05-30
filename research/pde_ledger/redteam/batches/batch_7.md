---
batch: 7
pass: IV.x/V.1 orchestrator-direct integrity remediation
range: [148, 150, 157, 166]
total_stages: 4
verified: 4
findings_count: 7
findings_resolved: 7
findings_blocked_legitimate: 0
material_change_count: 1
material_change_stages: [148]
clean_stages: []
status_only: []
dirty_stages: [148, 150, 157, 166]
checkpoints: []
status_only_candidates: []
consult: redteam/codex_reviews/_consult_batch7.md
audit_date: 2026-05-29
verify_date: 2026-05-29
status: closed
---

# Red-team remediation batch 7 — stages 148, 150, 157, 166

## Summary

Seventh batch of the IV.x/V.1 orchestrator-direct integrity remediation. Four stages
whose first-pass fixes had been applied orchestrator-direct (Codex bypassed) were
Codex-reconciled and re-verified:

- **148** (IV.5, `representative_positive_families`)
- **150** (IV.5, `full_profile_residual`)
- **157** (IV.6, `core_mouth_coevolution_status`)
- **166** (V.1, `bundle_inversion_four_drifts`)

All four are dual-engine computational stages (SymPy + Mathematica). **NO checkpoints
and NO status-only units** in the batch-7 range — `is_checkpoint: false` and
`is_status_only_candidate: false` confirmed in `redteam/MANIFEST.yaml` for all four
(IV.5/IV.6/V.1 carry no checkpoints; nearest checkpoints are 096 and 105).

7 findings across the four stages — **all 7 resolved, 0 blocked**. All four
REMEDIATION-`verified` end-to-end; both engines exit 0/0; all committed outputs fresh.

**Material change:** 1 (stage 148), with **ZERO downstream propagation**. 148 carried a
genuine live bug — the Mathematica `dT` was *wrong* (a defective mirror dropped a chain
term) while SymPy and the paper were correct. The fix CORRECTED the Mathematica `dT` to
match the already-correct paper/SymPy `dT`; the correct `A_T`/`B_T`/`dT` did NOT move, so
no downstream stage is stale. 150/157/166 are `material_change: false`.

## Consult

One Claude+Codex read-only consult (`redteam/codex_reviews/_consult_batch7.md`,
Q1–Q6): **5 of 6 unconditional CONCUR + 1 conditional CONCEPTUAL-ESCALATE (Q3/157)
RESOLVED against escalation**. No DISPUTES. No item changes any stage's published
(paper-card) conceptual claim, so nothing was escalated to the user.

| Q | Stage / finding | Verdict |
|---|---|---|
| Q1 | 148-F2 live cross-engine divergence | CONCUR |
| Q2 | 148-F1 ξ_* bridge + tolerance | CONCUR (with scope) |
| Q3 | 157 canonical-even tautology | CONCUR — conditional escalate **resolved: NOT conceptual** |
| Q4 | 157 assumptions domain `0<σ<1` | CONCUR |
| Q5 | 166 matrix round-trip tautology | CONCUR — REPLACE (not keep) |
| Q6 | 150 compact display | CONCUR |

The Q3 escalation condition (escalate *iff* the stage TEXT claims an in-stage proof of
δC=0 from family motion) was resolved AGAINST escalation: `paper/stages/stage_157.tex`
is tagged `\StatusNumerical/\StatusOpen`, already calls the deviation-to-normalization
map "the next task" (deferred to Stage 158), and "imposes" the even-preservation
constraint rather than "proving from motion." Only the SCRIPT docstring item 6 overstated
relative to its own card — a scripts-only labeling fix. Paper card NOT edited.

## Per-stage findings tally

| Stage | Status | Findings | material_change | Notes |
|-------|--------|----------|-----------------|-------|
| 148 | dirty | 2 | **true** (zero downstream) | F2 LIVE bug: Mathematica `dSigmaOfDeltas`/`dTOfDeltas` dropped the S-follows-Π chain term → wrong `dTU=0.4976…`; de-mirrored to independent `D[]` autodiff + both engines anchor `A_T`/`B_T` to paper literals. F1 ξ_* bridge raised to exact symbolic zero on `100π²` (stale `168π²` directive typo purged). |
| 150 | dirty | 1 | false | F1 DISPLAY-only: compact slope `S_q(Π)=Aq·k−Cq·Π` printed from free placeholders then `.subs`; source slope already correct/committed. |
| 157 | dirty | 3 | false | F1/F2 canonical-even tautology de-mirrored in both engines → parallel det non-degeneracy `det([[1,−9σ],[5,−72σ]])=−27σ≠0`; F3 added branch domain `0<σ<1`; docstring item 6 + a carry-forward banner corrected to the card's existing Stage-158 deferral. |
| 166 | dirty | 1 | false | F1 vacuous matrix round-trip (`Mmat·Inverse[Mmat]·v−v≡0`, X−X) replaced by hand-typed forward-transcription check from notes §1 boxed laws. |

**Totals:** 7 findings, 4 dirty stages, 0 clean, 0 status-only, 0 blocked.

## Findings closed (detail)

### 148 — the live cross-engine bug

- **F2 (insufficient_verification, the LIVE bug):** the Mathematica first-order
  traction-shift route (`dSigmaOfDeltas`/`dTOfDeltas`, wl:43-47) routed `dG` only through
  `dPi=−dG/gPrimeStar` and **silently dropped the S-follows-Π chain term** that SymPy's
  `AT` carries, so Mathematica computed a WRONG `dTU=0.4976…` while SymPy was correct at
  `dTU=0.5087…` — and **nothing asserted cross-engine agreement**, so both engines passed
  while disagreeing. Fixed by DELETING that block and having Mathematica derive `aT`/`bT`
  by its OWN `D[]` autodiff of `Tm[p]=Sqrt[(9/20)(p/(1−sFormula/4))]` along the S(Π) curve
  (regenerating the dropped `sFormula'` chain term) — NOT a hand-port of SymPy's algebra.
  BOTH engines now anchor `A_T`/`B_T` to the published paper literals
  `A_T=−4.27263956256927`, `B_T=0.134875005736706`
  (`paper/appendices/stage_appendix_part04.tex:846,848`) as the EXTERNAL cross-engine
  anchor; no baked SymPy literal in the `.wl`; `dTU`/`dTD` now agree across engines to ~16
  digits (old divergent `0.4976…`/`−0.1144…` gone).
- **F1 (paper_misalignment / insufficient_assertion):** the same-source ξ_* bridge check
  `(1−λ_{Π,0})−ξ_*` was raised to an EXACT symbolic zero (`exact_resid == 0`, built from
  the exact `gminus`) on the rF1-forced `100π²` radical (`12·(37/20)²=4107/100` ⇒
  `rF1²=(4107−100π²)/(100π²)`). The stale `168π²` was a directive-doc typo (a Codex
  false-positive) — the scripts always used `100π²`; the directive was corrected.

**material_change rationale:** the buggy Mathematica audit was corrected to match the
already-correct paper/SymPy values; the correct `A_T`/`B_T`/`dT` (paper literals) did NOT
change, so **no downstream stage is stale**.

### 150 — display-only

- **F1 (insufficient_verification, display corroboration):** both transcripts now print
  the compact slope `S_q(Π)=Aq·k−Cq·Π` built from FREE coefficient placeholders, then
  `.subs`/`/.` to the concrete definitions that feed the load-bearing `T_q'(0)−S_q==0`
  assertion — so the printed form is provably the real slope (can-fail self-test: deleting
  the `.subs` step leaves a free-symbol form and the assert FAILs). The source slope was
  already correct and committed; this is a display-only change with no derived result moved.

### 157 — de-tautologized canonical-even check

- **F1/F2 (insufficient_verification / transliteration):** the duplicate homogeneous
  re-solve (SymPy) and the mirrored literal 9/72/5 system (Mathematica) were replaced by a
  parallel determinant non-degeneracy assertion `det([[1,−9σ],[5,−72σ]])=−27σ≠0` (a genuine
  fail mode — a mistyped −72→−71 gives `−36σ`, residual `−9σ≠0`, FAIL).
- **F3 (symbol_assumption_error):** added the physical branch assumption `0<σ<1` scoped to
  the Section-3 even-preservation block.
- Plus the SymPy docstring item 6 corrected (was overclaiming "tangent motion kills delta
  C") to match the published card's EXISTING deferral of the deviation-to-normalization map
  to Stage 158, and a carry-forward banner fixed ("138-139"→"155-156"). Resolved against
  escalation (see Consult, Q3); paper card NOT edited.

### 166 — forward-transcription replaces vacuous round-trip

- **F1 (tautological_check):** the vacuous Mathematica matrix "round-trip"
  `Total[(Mmat·Inverse[Mmat]·v − v)^2]` (an X−X self-cancellation true for any invertible
  `Mmat`) was **replaced** (not supplemented) by a genuine forward-transcription check
  `Total[(Mmat.{drho,da,dcs,dZ} − fwdLaws)^2]==0`, with
  `fwdLaws={2 drho, drho+2 da, dZ+2 dcs−2 da, 5(dcs−da)}` HAND-TYPED from the notes §1
  boxed laws (not built from `Mmat`/`Solve`/`Inverse`). A wrong `Mmat` coefficient now
  propagates into the residual and FAILs.

## Orchestrator catches (directive-review, stage 148)

Two orchestrator catches on the 148 directive, both during directive review before
`codex-invoke`, both Codex-agreed:

1. **Fragile baked-literal → symmetric paper anchor.** An early directive route would have
   anchored the Mathematica `dTU`/`dTD` against a BAKED SymPy literal in the `.wl` (a
   fragile, non-independent cross-engine pin). Replaced with a symmetric EXTERNAL anchor:
   BOTH engines check their `A_T`/`B_T` against the published paper literals
   (`stage_appendix_part04.tex:846,848`), so neither engine imports the other's value.
2. **Flawed exact-route code → fully-symbolic construction.** The exact ξ_* bridge route
   was reworked to build the bias-neutral point from the exact `gminus`
   (`rF1_exact`/`gminus_exact`/`one_minus_lam_exact`) and assert `exact_resid == 0`
   fully symbolically, rather than chasing a loose numeric tolerance on a residual that
   does not actually need `Pi_star`.

## Verification

All four verification files written under `redteam/verifications/stage_{148,150,157,166}.md`.
Final verdicts:

- `verified` (4): 148, 150, 157, 166.
- `needs_rework` (0); `blocked_unfixable` (0).

Material change: 1 (stage 148 — defective Mathematica `dT` corrected to the already-correct
paper/SymPy value; ZERO downstream propagation, no stage stale). 150/157/166 are
`material_change: false`.

## Tracker sync

Six prose trackers synced for batch 7:

- `MATHEMATICA_MIRROR_POLICY` — 148 (`dT` DE-MIRRORED to independent `D[]` autodiff +
  paper-literal anchor; prior mirror was DEFECTIVE, dropped a term), 157 (canonical-even
  DE-MIRRORED to parallel det non-degeneracy), 166 (round-trip → forward-transcription,
  X−X self-cancellation removed), 150 (no mirror issue — display only).
- `CHECKPOINT_TRUST_AUDIT` + `CHECKPOINT_CONSTANT_PROVENANCE` — no checkpoints in range;
  148's `A_T`/`B_T` provenance = paper literals (`stage_appendix_part04.tex:846,848`) now
  anchored in BOTH engines.
- `PAPER_CLEANUP_TRACKER` — **P4-48**: NO new paper-card items (paper `A_T`/`B_T` correct,
  scripts fixed to match); flagged the THREE misfiled notes-review artifacts
  (`notes/stages/review/stage_{148,150,157}_review.md` bodies point at pre-renumber
  `stage029`/`stage031`/`stage038` source files) for separate orchestrator/notes repair.
- `EM_PROJECTED_INTEGRATION_TRACKER` — no impact (out of EM-projected range).
- `STAGE_VERIFICATION_COVERAGE` — 148/150/157/166 → verified (remediation); 148
  `material_change: true` (zero downstream propagation), 150/157/166 `material_change: false`.

## Cumulative

Range 001-175 paper-aligned at v2 depth and remediation-hardened through batch 7. Next
remediation batch (sequential-audit-chunks rule, awaits explicit user authorization):
the remaining FINDINGS stage **175**.
