# Claude+Codex consult — BATCH 7 (stages 148, 150, 157, 166)

Date: 2026-05-29. Mode: read-only ephemeral (`codex-chat -s read-only`). Raw deleted after
distillation (bloat policy). Prompt: `redteam/tmp_prompts/consult_batch7.md`.

## SUMMARY
**5 of 6 unconditional CONCUR; Q3 (157) flagged CONCEPTUAL-ESCALATE *conditionally* — condition
resolved by the orchestrator AGAINST escalation (see Q3).** No DISPUTES. No item changes any
stage's published (paper-card) conceptual claim.

| Q | Stage / finding | Verdict |
|---|---|---|
| Q1 | 148-R2 live cross-engine divergence | CONCUR |
| Q2 | 148-R1 ξ_* bridge + tolerance | CONCUR (with scope) |
| Q3 | 157-R1/R2 canonical-even tautology | CONCUR — option (ii); conditional escalate **resolved: NOT conceptual** |
| Q4 | 157-R3 assumptions domain | CONCUR |
| Q5 | 166-R1 matrix round-trip tautology | CONCUR |
| Q6 | 150-R1 compact display | CONCUR |

## Resolutions to encode in the rewritten directives

### Q1 — 148-R2 (the live bug)
- **SymPy is paper-consistent; Mathematica `dSigmaOfDeltas`/`dTOfDeltas` (wl:43-47) is the bug** —
  it routes `dG` only through `dPi=-dG/gPrimeStar` and never uses `sPrimeStar`, dropping the
  S-follows-Π chain term that SymPy's `AT` (py:32-35) carries. Paper anchors `A_T=-4.27263956256927`,
  `B_T=0.134875005736706` (part04:846/848); model `A_T(ḡ-g_*)+B_T(S̄-S_*)` (part04:839-840).
- **Fix:** Mathematica computes `aT`/`bT` INDEPENDENTLY via `D[]` autodiff of
  `Tm[p_]:=Sqrt[(9/20)*(p/(1-sFormula[p]/4))]` (S substituted as the Π-function, à la Stage 147) —
  NOT a hand-port of SymPy's algebra — then `dT=aT*dg+bT*dS`. ADD hard assertions: `aT≈-4.27263956256927`
  and `bT≈0.134875005736706` (paper literals, external anchor) AND `dTU`/`dTD` agree with SymPy to ≥25 digits
  (currently NOTHING asserts cross-engine agreement → both pass while disagreeing).
- **Q1b — no double-counting** (CONCUR): `A_T` is the paper's total-derivative g-projection coefficient;
  `B_T` multiplies the full projected S-deviation in the same validated coordinate model (the projection
  identity + source-centering already verified at Stage 147 / batch 6). NOT a conceptual issue.

### Q2 — 148-R1 (ξ_* bridge)
- Keep the `(1-λ_{Π,0}) - ξ_*` check as a **same-source closed-form consistency check** (acceptable; both
  sides descend from `gMinus`). Raise the SymPy side to a real `assert`.
- **`100π²` is CORRECT** (`12·(37/20)²=4107/100` ⇒ `rF1²=(4107-100π²)/(100π²)`); purge the directive's
  stale `168π²`. Do NOT touch wl:36 `MaxIterations->100`.
- **Tolerance:** Codex's preferred route — an EXACT symbolic comparison (`gLam==gMinus` reduction) is cleaner
  than chasing `<1e-25`, since the bridge does not actually need `Pi_star`. Acceptable alternative: raise
  `nsolve` working precision (`prec=50`) so the residual honestly drops below `1e-25`. Either closes the
  review's tolerance finding; prefer the exact route if SymPy reduces it cleanly, else the prec-50 numeric.

### Q3 — 157-R1/R2 (canonical-even tautology) — **escalation condition RESOLVED, not conceptual**
- py:107-110 already solves+asserts `{dE2=0,dE4=0}=={δC:0,δκ:0}`; py:112-114 RE-SOLVES the identical
  homogeneous pair (det `[[1,-9σ],[5,-72σ]]=-27σ≠0`) ⇒ `δC≡0` by construction → tautological (R1).
  wl:102-105 re-states the same literal numerators → mirror (R2).
- **Codex: option (ii) is the honest fix** — collapse to a SINGLE non-degeneracy/trivial-kernel assertion +
  RELABEL honestly as "canonical-even constraint pins δC=δκ=0 (carried-coefficient consistency)", deferring
  the genuine tangent→defect (deviation-to-normalization) map to Stage 158 (py:116-123 / wl:107-115 is
  already the Stage-158 expansion point). Option (i) [reconstruct 9/72/5 from an upstream Galerkin source
  in-script] only if that source is actually present here — it is NOT (no in-stage provenance).
- **Codex flagged CONCEPTUAL-ESCALATE *iff* the stage text claims an in-stage proof of δC=0 from family
  motion.** ORCHESTRATOR DETERMINATION (resolved against escalation), evidence from `paper/stages/stage_157.tex`:
  - `:5` stage is `\StatusNumerical/\StatusOpen` (not an unconditional theorem);
  - `:16` "leaves the deviation-to-normalization map as **the next task**" (map already declared deferred);
  - `:23` Checks: even-preservation constraints are **"imposed"** (= the kernel check), NOT "proven from motion";
  - `:27` "the card is a derivation ledger entry, **not an unconditional actual-branch theorem**."
  The published claim is the weaker "constraint imposed ⇒ δC=0 + map deferred". Only the SCRIPT DOCSTRING
  item 6 ("Tangent motion kills delta C") overstates relative to its own card. Option (ii) brings the script
  INTO alignment with the already-deferring card → it does NOT change Stage 157's conceptual claim.
  **Therefore: how-it's-checked/labeling fix, resolved Claude+Codex; NO user escalation.** (Docstring item 6
  is a .py comment → scripts-only fix-loop territory; the paper card is NOT edited and already says the right thing.)

### Q4 — 157-R3 (assumptions)
- Add `0 < sigmaStar < 1` to wl:93 `$Assumptions` (claim invalid at σ=0 for δκ; `(1-σ)` singular at σ=1).
- Minimal safe form: simplify the residual under that assumption before `expectZero`; if a
  `ConditionalExpression[0, cond]` appears, unwrap only after verifying `cond` is implied by `0<σ<1`
  (per the Mathematica-idiom memo).

### Q5 — 166-R1 (matrix round-trip)
- Replace wl:73-76 (`Total[(Mmat.solVec - v)^2]`, vacuous since `solVec=Inverse[Mmat].v`) with a
  HAND-TYPED forward-transcription check `Total[(Mmat.{drho,da,dcs,dZ} - fwdLaws)^2]==0`,
  `fwdLaws={2 drho, drho+2 da, dZ+2 dcs-2 da, 5(dcs-da)}` from notes §1 boxed laws (NOT built from
  Mmat/Solve/Inverse). Confirmed `Mmat` rows 3,4 = `{0,-2,2,1}`,`{0,-5,5,0}` in column order
  `(drho,da,dcs,dZ)` (wl:61-66) are correct. **Replace** the round-trip (don't keep the misleading line).

### Q6 — 150-R1 (compact display)
- Source slope `Sq=Aq*k-Cq*Pi` already committed; only the DISPLAY expands (Aq/Cq concretized). Print a
  compact placeholder form (build from free `Aq,Cq,k`), THEN `.subs` the concrete definitions for the
  load-bearing `T_q'(0)-S_q==0` assert, so the printed form is provably the real slope. Approach (b).
