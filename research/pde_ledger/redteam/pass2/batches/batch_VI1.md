---
batch: VI.1
pass: 2
range: 201-218
total_stages: 18
verified: 18
findings_count: 1   # one load-bearing finding (211 transliteration/port); the rest are stale_output (refresh-only, P4-52) + a 16-card card-text-lag cluster (user-deferred to P4-51) + one numbering cross-ref (206, deferred)
material_change_count: 0
clean_stages: [201, 204, 207, 208, 209, 215, 216, 217]
status_only: []
dirty_stages: [211]   # the only stage with a source-code change (the .wl re-author); all others had refresh-only stale_output and/or the deferred card-lag
checkpoints: [203, 218]
value_recon_total: 222
value_recon_misaligned: 0
audit_date: 2026-06-09
verify_date: 2026-06-09
status: closed
---

# Red-team PASS-2 batch VI.1 — Explicit realization, scalar slice, ray ranking

## Summary

Pass-2 re-audit of VI.1 (`Part VI.1 — Explicit realization, scalar slice, ray
ranking`), stages 201–218 (18), isolated under `redteam/pass2/`. **All 18 reached
`verified`; `material_change: false` on all 18** (the only source change is the 211
`.wl` re-author — a method/route change; no derived value, constant, identity target,
or paper number moved). 0 stop-cold, 0 blocked, 0 Codex deviations; 211's codex-invoke
exit 0 on iter-1 (no iter-2). 36 orchestrator exec runs exit 0 (reliability gate);
slowest 217-mma 148s, 214-mma 72s, 218-mma 36s, 202-mma 14s (the first-pass iter-2
timeout-rework stage — confirmed deterministic, well under the 600s cap).

Value reconciliation (pass-2 augmentation): **222 deliverable values checked
batch-wide, 0 misaligned** (201=11, 202=11, 203=13, 204=16, 205=15, 206=10, 207=7,
208=16, 209=13, 210=18, 211=14, 212=9, 213=12, 214=6, 215=11, 216=14, 217=11, 218=15).
No genuine `paper_misalignment` → **ZERO substantive paper/notes edits this batch.**

## Headline — 211 `.wl` re-author (the lone port; first-pass retrofit was insufficient)

VI.1 was the first pass's dual-engine RETROFIT-heavy batch (16 new `.wl` + 218
re-authored from a transliteration + 203 strengthened). The pass-2 heads-up flagged
the standing risk — realized on V.3-200 — that a first-pass retrofit can be
insufficient. Here it was realized on **one non-checkpoint stage, 211**.

- **The per-stage audit agent FLAGGED 211** as `mathematica_transliteration` (medium).
- **The orchestrator ground-truth `.wl`-vs-`.py` read CONFIRMED the port:** the `.wl`'s
  M1 computed `Numerator[Together[D[Phi,r]]]` and checked it against the *same posited*
  `Nr` (identical to the `.py`'s `sp.diff(Phi,r)` check); M2/M3 posited the *same*
  `Ccross`/`Sr`/`Ss` and asserted the *same* algebraic identities; M4 the same Bézout
  product; M5/M6 the same iso/sym substitutions. **No independent extraction operation
  anywhere** — exactly the V.3-200 / V.2-port signature.
- **CALIBRATION-AGENT technique (V.3 method) across the ray-ranking family** (207, 208,
  209, 210, 213, 214, 216) with 211 as the explicit known-port standard → **PORTS
  BESIDES 211: NONE.** Each sibling clears the bar via a genuine independent op the
  `.py` lacks: `Resultant` lift-elimination (209, 214), Lagrange-`Solve` (210, 213,
  216), `Solve`+`Reduce[…,Reals]` branch-select (207), `Maximize`/`Eigenvalues`
  (213/216), `Solve`+uniqueness (208). (Caveat logged: 214's M1 numerator block IS a
  211-style autodiff-vs-posit port, but 214's *theorem* — the eliminant — is
  independently derived via `Resultant`, so 214 is independent overall.)
- **Re-author-vs-accept = USER-LEVEL → ⭐ USER AUTHORIZED re-author** (matching the
  V.2/V.3 precedent and `feedback-dual-engine-required`).
- **Codex (designs+writes) re-authored ONLY the `.wl`** (`.py` untouched): the
  load-bearing eliminants are now DERIVED by `Resultant` of the differentiated
  stationary numerators — `crossDerived = Resultant[numRq, numSq, q]` (quartic),
  `srDerived/ssDerived = Resultant[numRq/numSq, q^2 - Delta, q]` (sextics), with
  `numR/numS = Numerator[Together[D[Phi,var]]]` and `Sqrt[Delta]->q`. The SymPy closed
  forms survive ONLY as labeled "SymPy comparison targets" (`derived cross / SymPy
  target = -1`, scaled residual `= 0`). Verifier-confirmed GENUINELY INDEPENDENT;
  derived cross degree 4, derived S_r/S_s degree 6, Bézout 24 from the *derived*
  degrees, diag-iso (M5) and symmetric equal-mix `Nr(1,1)=Ns(1,1)=0` (M6) preserved;
  committed `.wl` output refreshed to the derived forms. **211 ADDED to the
  Independent-Mirror Set.** `material_change: false`.

## Checkpoints — both CLEARED (203 & 218)

- **203 (free quintuple scalar closure slice + crossing theorem) — CLEARED.** The §VI
  crossing theorem is re-derived in-script (carried Stage-197 closure scalar composed
  with the concrete Stage-202 path → `32^(2τ-1)-1`, sign change `-31/32 → 31`, unique
  root `τ=1/2`); both engines substantive; the `.wl` is independent (log-additive/posit
  + `Reduce` + target-monomial invariance vs the `.py`'s power-multiplicative/derive +
  `solveset`). 13-identity reconciliation, 0 misaligned.
- **218 (full support-≤5 completion + local mixed-ray search closure) — CLEARED by
  orchestrator ground-truth read** (the hardest stage; pass-1 re-author from a
  transliteration was **SUFFICIENT**, unlike V.3-200). M1 boundary incidence via
  `Subsets`/`ContainsAll`/`Boole`/`Tally`/`2^5-2` vs `.py` itertools; **M2-M3 splice via
  `Resolve[ForAll,…],Reals]` real-QE vs the `.py`'s `simplify_logic` finite
  propositional contradiction**; M4 regimes via independently-GENERATED witness windows
  (`makeBoundaryWindows`) vs hand-listed; M5 benign shared budget arithmetic. Budget
  ledger (162→324, 750→1500, 1140+324=1464, 1140+1500=2640) computed in both engines,
  internally consistent, and appendix-exact.

## Reconfirmations (the 5 first-pass user-resolved paper_misalignment — all HOLD)

- **217 PUBLISHED `179/230 → 162` ✅ HOLDS.** Both engines DERIVE `162` (product of the
  computed degree pattern `(3,3,3,3,2)` = `3^4·2`, never hardcoded); a grep for `179`/
  `230` across all 7 artifacts (card, notes, part-06 appendix, both scripts, both
  outputs) returned ZERO matches. Arithmetic chain intact: `162 → 324 (=2·162) → 1464
  (=1140+324)`; fallback `750 → 1500`, `1140+1500 = 2640`.
- **212 notes typos HOLD** — 188→120 fixed; the 246/243/245/247→212/209/211/213
  renumber present; budget 120/480/600 in both engines; no stray tokens.
- **214 notes typos HOLD** — projected bound `150 = 5·5·6`, lifted `54`; no stray
  218/162/230.
- **206** — first-pass scope finding was resolved by ADDITIVE both-engine checks; no
  value mismatch (confirmed; the `.py` strict-sign/`Resolve` checks the `.wl` adds).

## Non-blocking classes (established deferral policy)

- **stale_output (P4-52 class)** — 13 committed SymPy outputs (201–204, 206–210, 213,
  214, 217, 218) + 203-mma carried pre-renumber `STAGE 184–201` banners (SOURCE banners
  already canonical). The orchestrator reliability-gate re-run + sed-refresh regenerated
  them; all 18 self-banners now read `STAGE NNN`; 205/212/215/216 SymPy + the fresh mma
  outputs byte-identical. Arbiter grep CLEAN (no stale self-epoch banner, no `168`
  class, no FAIL).
- **Card-text lag (16 cards: 201–217 except 203; 203 & 218 already cite their `.wl`)** —
  `\stagefield{Verification}` says "Mathematica audit: none yet" though a passing `.wl`
  exists. A stale STATUS annotation, NOT a value/identity mismatch → **USER-DECIDED
  2026-06-09 to DEFER to PAPER_CLEANUP P4-51** (appended to the existing card-text-lag
  entry; fixed in the later dedicated paper pass, Codex-applied + Claude-reviewed).
- **Numbering (content-keyed, never offset-sweep → NUMBERING_SCRIPT_OUTPUT_BAND_PLAN)** —
  206 `.py` collapse-target label "Stage 239" should be the Stage-205 quad-log predictor
  (math correct); 218 dead `source_stage:198/200` dict fields + "Stage 249" comment
  (never asserted/read); 203 `.wl` `chiFromStage180`/`closureNumStage180` var names vs
  the Stage-197 attribution. All cosmetic, deferred.
- **+17 `\stagefield{Purpose}` card self-label class (IV.3 discovery)** — VI.1 cards
  CLEAN (218–235 absent).

## Infra

- 36 orchestrator exec runs exit 0 (reliability gate). 15 committed outputs refreshed
  (13 stale SymPy banners + 203-mma banner + 211-mma re-authored content); the rest
  byte-identical. Arbiter grep CLEAN.
- Seat policy: 211 = 1 `.wl`-touching Codex session solo (held one math seat); the
  orchestrator exec re-run ran sequentially AFTER Codex done (no overlap). 202's exec
  re-confirmed deterministic <600s.
- A wrong-project-root path-typo (`toy_*`, not `toy_physics`) an audit agent leaked into `directives/stage_215.md` was
  purged; the wrong-root path-typo grep is empty (the only hit is a gitignored prior-batch
  `codex_logs/053_iter1.txt`).
- 6 trackers synced (PAPER_CLEANUP **P5-19** = ZERO substantive paper/notes edits;
  211 added to the Independent-Mirror Set). Pass-1 `MANIFEST.yaml` untouched.
