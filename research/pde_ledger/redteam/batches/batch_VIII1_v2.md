---
batch: VIII.1
range: 243-253
total_stages: 11
verified: 11
findings_count: 16
findings_resolved: 16
findings_blocked_legitimate: 0
material_change_count: 0
clean_stages: []
status_only: []
dirty_stages: [243, 244, 245, 246, 247, 248, 249, 250, 251, 252, 253]
checkpoints: [243, 248, 253]
consult: none (248 F2 reframed via Claude math-coverage resolution after Codex correctly blocked a Wolfram-version dead-end Reduce/ToRules route; non-conceptual, no paper edit)
audit_date: 2026-06-03
verify_date: 2026-06-03
status: closed
---

# Red-team batch VIII.1 — Relaxed branch, dynamic event chain, cold survival

> `findings_count: 16` is a CURATED THEMATIC count (consistent with the
> convention — VII.2 used 18). The total per-directive findings across the 11
> stages sum to ~38; this field tracks the distinct thematic defect classes
> closed, not the raw per-directive total.

## Summary

11-stage audit unit for VIII.1 (`Part VIII.1 — Relaxed branch, dynamic event
chain, cold survival`), **the FINAL forward first-pass batch**, forward
first-pass under the v2 paper-grounded auditor **WITH the dual-engine rule** in
force. **Three checkpoints this batch — 243, 248, 253** (all at the higher
verify bar); no status-only units. All 11 reached `verified`; **`material_change:
false` on all 11**; **0 stop-cold, 0 ultimately-blocked, 0 needs_rework left
open.** Every change is a strengthening / route change / notes typo correction;
no derived value, constant, identity target, or PUBLISHED paper number on the
SCRIPT side moved, so no `upstream_stale` propagation.

**10 of 11 codex-invoke runs exited 0 on iteration 1; one iteration-2 (stage 248
— see the reframe narrative below).**

### Milestone — first end-to-end red-team pass COMPLETE

**Reaching stage 253 COMPLETES the first end-to-end red-team pass.** Cumulative:
**242/253 → 253/253 stages red-team verified (100%)**; the entire range 001–253
is now paper-aligned at v2 depth. There are no `pending` stages remaining. The
planned full end-to-end **second pass** remains a later cross-check.

### Headline — full dual-engine coverage with zero sanctioned mirrors

VIII.1 ran with the dual-engine rule (a Mathematica `.wl` is **REQUIRED wherever
Mathematica CAN independently verify** — the test is "is it possible," not "is
it necessary") in force. **Every stage now has an INDEPENDENT second engine; 0
sanctioned mirrors were accepted in VIII.1.**

- **8 stages got a NEW independent-route `.wl`**: 244, 245, 246, 247, 249, 250,
  251, 252.
- **All 3 checkpoint `.wl` were caught as line-by-line transliterations and
  RE-AUTHORED to genuinely independent routes**: 243, 248, 253.

Every new/re-authored `.wl` was confirmed independent by a clean verify agent —
native primitives via a DIFFERENT decomposition than the SymPy `.py`.

The labor split was strictly enforced: **Claude reviews** (audit + verify);
**Codex writes ALL script code** (designs and writes the new/re-authored `.wl`);
the directives stated only the requirement + acceptance criteria, never script
code.

## Checkpoint findings (the higher bar, no rubber-stamp)

All THREE VIII.1 checkpoints had EXISTING `.wl` that were line-by-line
transliterations of their `.py`; all THREE were RE-AUTHORED to independent
routes and all THREE cleared the higher bar. **No checkpoint constant CHANGED —
the route-only checkpoint rewrites land on the SAME values.**

### 243 — relaxed-branch lift / leakage-work lane / non-rigid solve / recovery / short-range firewall

Cleared the higher checkpoint bar:

- **`.wl` RE-AUTHORED** from a transliteration to a genuinely independent route:
  IBP closure / native `LinearSolve` / `TrigExpand` / `Series`-at-∞ for the
  exact Gaussian leakage/work lane and the linear `(U,V)` non-rigid solve, with
  the hardcoded `expected*` residuals DELETED (so the checks are now genuine
  re-derivations, not equality-to-a-baked-literal).
- In-file stale-label cosmetic rode the fix loop: `.wl` banner `STAGE 226→243`
  (single-file, NOT a batch renumber).

**Checkpoint bar MET.**

### 248 — dynamic event chain / energy conservation / threshold speeds / Coulomb reference / near-top action

Cleared the higher checkpoint bar (the one iter-2 of the batch):

- **`.wl` RE-AUTHORED** from a transliteration to a genuinely independent route:
  the §II `Solve`-mirror was replaced by a native **SATISFACTION route** — the
  compiler closed forms are verified to satisfy their defining energy equalities
  via substitution + `FullSimplify`, plus a **non-vacuity guard** and a
  **positive-branch guard**.
- **ITER-2 reframe (see the dedicated narrative below):** iter-1 prescribed a
  `Reduce`/`ToRules` route which Codex correctly **BLOCKED** as a Wolfram-version
  dead end; the orchestrator reframed to the satisfaction route above.
- In-file stale-label cosmetic rode the fix loop: `.wl` banner `STAGE 231→248`
  (single-file).
- A NOTES-ONLY benchmark typo (notes:506 `×168% → ×100%`) was corrected to the
  already-correct script (notes' own 23.3128% + script ×100); recurs the stale
  "168" seen at 148/232.

**Checkpoint bar MET.**

### 253 — lattice-turnover / calibration recovery / stiffness map / temperature ceiling / screening ratios

Cleared the higher checkpoint bar:

- **`.wl` RE-AUTHORED** from a transliteration to a genuinely independent route:
  native `D[Log[V[r]]]` + 5 `Solve` energy/force balances + regrouped threshold
  / Korringa / screening blocks, versus the `.py`'s back-substituted
  `r_turn_phys` expr−expr round-trip.
- In-file stale-label cosmetic rode the fix loop: `.wl` banner → 253
  (single-file).
- Two NOTES-ONLY benchmark typos were corrected to the already-correct
  script/card: notes:274 `187.23361317 → 119.23361317` (notes' own
  65.45193926/0.5489386551 = 119.2336 + cross-engine `.wl`) and notes:419 a_int
  `10.95423247 → 10.95423248` (= 4·K_turn, cross-engine).

**Checkpoint bar MET.**

## The dominant defect theme — the "variable-independence self-test trap" (+ tautological round-trips)

Continuing the VII.2 theme, the recurring VIII.1 defect was support-blindness /
independence "verified" by differentiating an expression w.r.t. variables it
never contains (vacuous), plus tautological / round-trip checks that pass for any
input. Hit across both the new-`.wl` stages and the re-authored checkpoints:

- **244-F1, 245-F1 (variable-independence self-test trap)** — support-blindness
  "verified" by differentiating w.r.t. an absent variable → vacuous; fixed via
  free-symbol containment + positive control (244), or a live-channel positive
  control (245).
- **246** — a σ_min-vs-its-own-piecewise tautology → tied to an independent
  derivation.
- **247** — an inverse round-trip + a back-substituted `lambda_L` (X−X) →
  de-tautologized.
- **249** — subtraction-by-linearity over blank placeholders (a vacuous
  difference) → de-vacuized.
- **251** — a quintic X−X self-cancellation → de-tautologized.
- **252** — a `gamma_safe_eq` X−X → de-tautologized.
- **253** — the `r_turn_phys` expr−expr round-trip (the checkpoint defect; see
  above).
- **250** — a global claim tested only at a single sample point → fixed to a
  GLOBAL strict-monotonicity proof via `Resolve[ForAll]`.

### Stage 248 iter-2 (notable — the one iteration-2 this batch)

On iter-1, the F2 directive prescribed a `Reduce`/`ToRules` route to extract the
event-chain compiler closed forms. Codex correctly **BLOCKED** it: on this
Wolfram version `ToRules` rejects the domain-predicate conjunction the route
requires — a genuine version-level dead end, not a refusal to do the work. The
**ORCHESTRATOR REFRAMED F2** (Claude math-coverage resolution, NON-conceptual, no
paper edit) to a native **SATISFACTION route**: rather than re-deriving the
closed forms by `Reduce`, verify that the declared compiler closed forms
**satisfy their defining energy equalities** under substitution + `FullSimplify`,
and gate it with (a) a **non-vacuity guard** (the substitution does not collapse
the equality trivially) and (b) a **positive-branch guard** (the physical
turning-point / threshold branch is selected). Codex applied on iter-2; the
independent verify confirmed the route is genuine and non-vacuous.

## Notes-only paper_misalignment resolutions (5 typos across 4 stages, user-resolved)

**5 notes-only `paper_misalignment` numerical/coefficient typos, user-resolved
(direction: correct the notes to the script/canonical value; PUBLISHED paper
cards UNAFFECTED; each cross-engine or internally corroborated):**

- **244** — notes:366 `196√2 → 128√2` (consistent with the script's `E0 = 16`
  structure).
- **247** — notes:406 Δ `210.17750000 → 142.17750000` (the notes' own formula
  9·16 − 1.35² = 142.1775, plus the adjacent D0 = 3.76481862, both corroborate).
- **248** — notes:506 `×168% → ×100%` (the notes' own 23.3128% + the script's
  ×100 give the consistent value; recurs the stale "168" previously corrected at
  148 IV.5 and 232 VII.2).
- **253** — notes:274 benchmark `187.23361317 → 119.23361317` (the notes' own
  65.45193926/0.5489386551 = 119.2336, plus the cross-engine `.wl`) AND notes:419
  a_int `10.95423247 → 10.95423248` (= 4·K_turn, cross-engine).

**Card-misattribution catch:** the audit agents initially flagged 247 and 253 as
PUBLISHED-CARD misalignments (`stage_247.tex:407`, `stage_253.tex:274`). The
orchestrator verified the cards directly — the 247 card is 93 lines and the 253
card is 140 lines, and **neither contains the flagged value** — so the typos were
NOTES-ONLY in both cases. All 5 corrections were Codex-applied + Claude-reviewed.

## Mathematica mirror policy — full dual-engine, zero sanctioned mirrors

VIII.1 ran with the dual-engine rule in force. **0 sanctioned mirrors were
accepted.** All 8 new `.wl` (244, 245, 246, 247, 249, 250, 251, 252) are GENUINE
INDEPENDENT routes; all 3 checkpoint `.wl` (243, 248, 253), each caught as a
transliteration, were **RE-AUTHORED** to independent routes. All 11 are recorded
in the Independent-Mirror Set in `MATHEMATICA_MIRROR_POLICY.md`.

## Per-stage findings tally

| Stage | Status | Findings | Notes |
|-------|--------|----------|-------|
| 243 | dirty (ckpt) | `.wl` re-authored | **Checkpoint.** Transliteration `.wl` RE-AUTHORED to an independent route: IBP closure + native `LinearSolve` + `TrigExpand` + `Series`-at-∞ for the leakage-work lane and the non-rigid `(U,V)` solve, with the hardcoded `expected*` residuals DELETED. In-file stale-label cosmetic: `.wl` banner `STAGE 226→243` (single-file). Checkpoint bar MET. Mathematica 43-44 PASS |
| 244 | dirty | 1 + `.wl` | New independent-route `.wl`. **F1 (variable-independence self-test trap):** support-blindness "verified" by differentiating w.r.t. an absent variable → fixed via free-symbol containment + positive control. **paper_misalignment (notes-only):** notes:366 `196√2→128√2` (script E0=16 structure). Mathematica 24 PASS |
| 245 | dirty | 1 + `.wl` | New independent-route `.wl`. **F1 (variable-independence self-test trap):** support-blindness → fixed via a live-channel positive control. Mathematica 38 PASS |
| 246 | dirty | 1 + `.wl` | New independent-route `.wl`. Tautology fix: σ_min-vs-its-own-piecewise → tied to an independent derivation. Mathematica 20 PASS |
| 247 | dirty | 1 + `.wl` | New independent-route `.wl`. Tautology fix: inverse round-trip + back-substituted `lambda_L` (X−X) → de-tautologized. **paper_misalignment (notes-only):** notes:406 Δ `210.17750000→142.17750000` (notes' own 9·16−1.35²=142.1775 + adjacent D0=3.76481862; card confirmed clean, 93 lines). Mathematica 19 PASS |
| 248 | dirty (ckpt) | `.wl` re-authored | **Checkpoint. ITER-2.** Transliteration `.wl` RE-AUTHORED: §II `Solve`-mirror → native SATISFACTION route (compiler closed forms verified to satisfy their defining energy equalities via substitution + `FullSimplify` + non-vacuity guard + positive-branch guard). Iter-1's `Reduce`/`ToRules` route correctly BLOCKED by Codex (Wolfram-version dead end — `ToRules` rejects the domain-predicate conjunction); ORCHESTRATOR REFRAMED to the satisfaction route (Claude math-coverage resolution, non-conceptual, no paper edit), Codex applied iter-2, verified. In-file stale-label cosmetic: `.wl` banner `STAGE 231→248`. **paper_misalignment (notes-only):** notes:506 `×168%→×100%` (recurs the stale "168"). Checkpoint bar MET. Mathematica 30-31 PASS |
| 249 | dirty | 1 + `.wl` | New independent-route `.wl`. Vacuity fix: subtraction-by-linearity over blank placeholders → de-vacuized. Mathematica 13 PASS |
| 250 | dirty | 1 + `.wl` | New independent-route `.wl`. **Single-sample-point trap:** a global claim was tested at one sample point only → fixed to a GLOBAL strict-monotonicity proof via `Resolve[ForAll]`. Mathematica 18-19 PASS |
| 251 | dirty | 1 + `.wl` | New independent-route `.wl`. Tautology fix: quintic X−X self-cancellation → de-tautologized. Mathematica 31-32 PASS |
| 252 | dirty | 1 + `.wl` | New independent-route `.wl`. Tautology fix: `gamma_safe_eq` X−X → de-tautologized. Mathematica 28 PASS |
| 253 | dirty (ckpt) | `.wl` re-authored | **Checkpoint.** Transliteration `.wl` RE-AUTHORED to an independent route: native `D[Log[V[r]]]` + 5 `Solve` energy/force balances + regrouped threshold / Korringa / screening blocks (vs the `.py`'s back-substituted `r_turn_phys` expr−expr round-trip). In-file stale-label cosmetic: `.wl` banner → 253. **paper_misalignment (notes-only):** notes:274 `187.23361317→119.23361317` (notes' own 65.45193926/0.5489386551=119.2336 + cross-engine `.wl`) AND notes:419 a_int `10.95423247→10.95423248` (=4·K_turn); card confirmed clean (140 lines). Checkpoint bar MET. Mathematica 26 PASS |

**Totals (thematic):** 16 distinct thematic defect classes closed (the
missing-`.wl` dual-engine gap on the 8 non-checkpoint stages; the three
checkpoint transliterations re-authored; the variable-independence self-test-trap
fixes on 244/245; the tautological/round-trip de-tautologizations on
246/247/249/251/252/253; the single-sample-point → global-monotonicity fix on
250; and the 5 notes-only paper_misalignment typos on 244/247/248/253), 0
blocked, 0 status-only. **8 new independent-route `.wl` (244, 245, 246, 247, 249,
250, 251, 252)** + **3 re-authored checkpoint `.wl` (243, 248, 253)**. (Per-stage
directive findings sum to ~38; the 16 above is the curated thematic count.)

## Paper / notes edits (Codex-applied, Claude-reviewed)

This batch APPLIED notes edits (per the file-ownership contract: Codex owns
`paper/*.tex` + `notes/stages/*.md` edits, Claude reviews). **ALL
paper_misalignment items were NOTES-ONLY; the PUBLISHED paper cards were
UNAFFECTED — they carry abstract forms.** Each correction was cross-engine or
internally corroborated. All were orchestrator-reviewed: correct, isolated, no
collateral.

**Numerical/coefficient typos:**

- **244:** notes:366 `196√2 → 128√2`.
- **247:** notes:406 Δ `210.17750000 → 142.17750000`.
- **248:** notes:506 `×168% → ×100%`.
- **253:** notes:274 benchmark `187.23361317 → 119.23361317` and notes:419 a_int
  `10.95423247 → 10.95423248`.

**In-file stale-label `.wl` fixes (single-file, NOT a batch renumber, rode the
fix loop):**

- **243** `.wl` banner `STAGE 226 → 243`.
- **248** `.wl` banner `STAGE 231 → 248`.
- **253** `.wl` banner → 253.

## Out-of-scope residuals — NO new numbering residual this batch

**VIII.1 introduced NO new numbering residual.** All 11 notes/stages H1 titles
are already canonical (Stage 243…253), and the only stale `.wl` banners were
fixed in-loop (243/248/253, single-file). The deferred PROJECT-WIDE stem-keyed
renumber reconciliation already logged for earlier batches — VII.1 notes drift
(219, 221, 222, 223-title, 227, 228, 229) + `.wl` banners (221, 227) at P4-53;
VII.2 notes-title drift (232, 234, 235, 236) at P4-54 — **STILL STANDS for the
post-253 dedicated paper-cleanup pass** per the project numbering-drift policy
(build a `stem→canonical` lookup ONCE from `scripts/`, map every citation
deterministically; Codex applies, Claude reviews). **NEVER offset-sweep** —
offsets are inconsistent across the EM-extension realignment.

## Capture-diff limitation

The committed batch patches OMIT the untracked new `.wl` files (244, 245, 246,
247, 249, 250, 251, 252) and the notes edits captured outside the diff snapshot;
the **verifiers checked the working tree directly** (reading each new/re-authored
`.wl` and each corrected notes line in place) rather than relying on the patch
capture, so the dual-engine independence and the notes corrections were
confirmed against the on-disk state, not the (incomplete) diff.

## Iteration-2 reworks

One: **stage 248** (the F2 `Reduce`/`ToRules` dead-end → native satisfaction
route reframe — see the iter-2 narrative above). The other 10 codex-invoke runs
exited 0 on iteration 1.

## Consult

None escalated to the user. The 248 F2 reframe was a Claude math-coverage
resolution (NON-conceptual, no paper edit — a faithful native satisfaction route
substituted after Codex correctly blocked a Wolfram-version dead-end
`Reduce`/`ToRules` route). The 5 notes-only paper_misalignment numerical typos
(244/247/248/253) were USER-RESOLVED 2026-06-03 — corrected to the already-correct
SymPy scripts, each cross-engine or internally corroborated.

## Verification

All 11 verification files under `redteam/verifications/stage_243.md` …
`stage_253.md`. Final verdicts:
- `verified` (11): 243–253.
- `needs_rework` → reworked → re-`verified`: 248 (iter-2).
- `blocked_unfixable` (0).

Per-stage Mathematica PASS (all exit 0, 0 FAIL): 243=43-44, 244=24, 245=38,
246=20, 247=19, 248=30-31, 249=13, 250=18-19, 251=31-32, 252=28, 253=26. SymPy
all exit 0.

Material change: **0** (`material_change: false` on all 11 — new/re-authored
second engine, de-tautologized / de-vacuized checks, and notes-only
numerical-typo corrections that align prose to the already-correct script; no
SCRIPT-side derived value, constant, identity target, or PUBLISHED paper number
moved; the three checkpoint re-authorings land on the SAME values).

## Cumulative

Range 001-253 paper-aligned at v2 depth. **253/253 stages red-team verified
(100%)** (was 242 after VII.2; VIII.1 adds 11 across 243–253). **THE FIRST
END-TO-END RED-TEAM PASS IS COMPLETE** — there are no `pending` stages remaining.
The FINAL forward batch: zero stop-cold, zero material change, full dual-engine
coverage with zero sanctioned mirrors, and ALL THREE checkpoints (243, 248, 253)
cleared the higher bar (each via a transliteration `.wl` re-authored to an
independent route, no checkpoint constant changed).

Next: the planned full end-to-end **second pass** (a later cross-check, now that
the first pass reaches stage 253), and the post-253 dedicated paper-cleanup pass
(the deferred project-wide stem-keyed renumber reconciliation, P4-53/P4-54).
