---
batch: VIII.1
pass: 2
range: 243-253
total_stages: 11
verified: 11
findings_count: 5
findings_resolved: 0
findings_deferred: 5
findings_blocked_legitimate: 0
material_change_count: 0
script_corrections: 0
codex_invocations: 0
clean_stages: [243, 245, 247, 248, 251, 253]
findings_stages: [244, 246, 249, 250, 252]
status_only: []
checkpoints: [243, 248, 253]
ports_found: 0
audit_date: 2026-06-10
verify_date: 2026-06-10
status: closed
pass_complete: true
---

# Red-team PASS-2 batch VIII.1 — Relaxed branch / dynamic event chain / cold survival (243–253) — THE FINAL BATCH

## Headline — ⭐ THE SECOND PASS IS COMPLETE (253/253)

**Closing VIII.1 completes the full end-to-end SECOND PASS — all 253 stages
re-audited at v2 depth with an independent dual-engine check.** VIII.1 is the
**THIRD consecutive ZERO-SCRIPT-CORRECTION pass-2 batch** (after VII.1, VII.2):
no `.py`/`.wl` edit anywhere, **0 Codex invocations**, all 11 captured diffs 0
bytes. All 11 stages `verified`, **`material_change: false` on all 11.**

VIII.1 was a first-pass **RETROFIT** batch (8 NEW independent-route `.wl` +
3 checkpoint `.wl` re-authored FROM transliterations). Pass-2 re-confirmed all 11
`.wl` genuinely independent — **0 ports, 0 borderline.**

## Independence — EARNED, not relayed (the standing user-distrust standard)

- **Orchestrator ground-truth `.wl`-vs-`.py` read IN FULL on ALL THREE
  checkpoints (243, 248, 253)** — the highest-risk re-authors (where V.3-200 was
  caught). All three confirmed genuinely independent (see below).
- The 11 per-stage audit agents each supplied operation-level, line-referenced
  independence evidence; no port was flagged → no CALIBRATION agent needed
  (as VII.1/VII.2).

## Checkpoints (the higher bar) — all THREE re-authors confirmed SUFFICIENT

The VI.1-218 / VII.1-221 / VII.2-239&242 outcome (sufficient), NOT V.3-200.

### 243 — relaxed constraint branch / short-range open-system compiler
- `.py`: `S_leak` via `sp.integrate(diff(W,w)·j_w)` (py:40), U/V via `sp.solve`
  (py:76), short-range limits via `sp.limit` (py:39).
- `.wl`: ADDS an **IBP-closure cross-check** `expectZero["IBP closure",
  Sleak + ibpInterior - boundary]` (wl:101, computes the integral a second way
  via `Integrate[W D[jw,w]]`) — a route absent from `.py`; derives U/V via BOTH
  `Solve` (wl:112) AND `LinearSolve[{{kU,-chiLam},{-chiLam,kV}},{fU,0}]` (wl:116)
  then asserts agreement (wl:132); short-range via `Limit` + an independent
  `Series[…,{x,Infinity,0}]` asymptotic route (wl:194-251). **Banner canonical
  (STAGE 243; the 226→243 fix landed). Checkpoint bar MET.**

### 248 — dynamic event chain / WKB threshold compiler
- `.py`: threshold via `sp.solve(Eq(E_launch_new, Vpeak), v0)` (py:73) — SOLVES.
- `.wl`: POSITS `vcritNew = Sqrt[2(Vpeak-V0)/ms]` (wl:97) then verifies the
  **satisfaction route** `EAtVcrit = ElaunchNew /. v0->vcritNew` with an explicit
  **non-vacuity guard** `TrueQ[FullSimplify[deltaNew==0]] || FreeQ[deltaNew,v0]
  → fail` (wl:110-111) — would FAIL if `deltaNew` were vacuously 0 or v0-free.
  One SOLVES, the other POSITS-and-verifies-satisfaction non-vacuously — opposite
  information flow (the pass-1 iter-2 reframe holds). **ZERO surviving 168** (the
  `×168%→×100%` notes fix holds; arbiter grep clean). Banner canonical (STAGE
  248; 231→248). **Checkpoint bar MET.**

### 253 — physical calibration / material threshold (THE FINAL STAGE)
- `.py`: `gamma_lat_turn_phys = gamma_lat_safe_eq/(Upsilon_lat·t_star)` (py:51) —
  DIVISION; chi_lambda hardcoded `2/r_phys`.
- `.wl`: `Solve`s the physical balance `gammaLatTurnPhys == zetaEp·lambdaEpOmegaD`
  for the threshold (wl:96-97); derives chi_lambda via `D[Log[V[r]], r]` with
  `V[rad_]:=(1/2) kEff rad^2` → `2/r` (wl:157-159); force-match `Solve` for
  `kEffReq` (wl:166-167). DERIVES what the `.py` DIVIDES/hardcodes. Benchmark
  `119.23361317476524` + `10.95423248` hold; ZERO surviving `187.2336`/
  `136.2336`/`10.95423247`. Published card clean/abstract (pass-1 misattribution
  re-overruled). Banner canonical. **Checkpoint bar MET.**

## Reconfirmations — 5 notes-only paper_misalignment fixes HOLD

Cross-engine-corroborated; **published cards/appendices UNAFFECTED** (notes-only):
- **244** — `196√2 → 128√2` (notes:366; both engines emit `128√2`).
- **247** — Δ `210.17750000 → 142.17750000` (notes:406; `210.17750000` purged
  everywhere; published card clean/abstract — pass-1's card misattribution
  re-overruled).
- **248** — figure-of-merit `×168% → ×100%` (zero surviving 168; the recurring
  stale-168 family 148/232/248 is purged here).
- **253** — benchmark `187.23361317 → 119.23361317` (notes:274) + a_int
  `10.95423247 → 10.95423248` (notes:419); the stale `136.23361317476524`
  caught in the first-pass tracker is gone (engines compute `119.…`).

## Reconfirmations — script-side fixes HOLD

- **Variable-independence self-test trap** (244-F1, 245-F1): positive controls
  (differentiate an expression where the variable genuinely appears, assert
  nonzero) present and load-bearing in both engines.
- **250 single-sample-point → GLOBAL**: the goldilocks-window / monotonicity
  claim is established GLOBALLY via `Resolve[ForAll]` / `Reduce` (M1/M3/M4),
  NOT a point evaluation — the pass-1 fix is intact.
- **Round-trip / tautology de-tauts** (246 independent `Min`/`MinValue`;
  247 falsifiable closure on the paper literal; 249 Möbius inverses genuinely
  `Solve`-derived; 251 difference-quotient certificate + independent 50-digit
  `NSolve`; 252 X−X → can-fail `safe_combo/sc == Γ₃sc²+Γ₅sc⁴`) — all hold.

## Findings (5) — all the standing-deferred CARD-TEXT-LAG class

5 `paper_missing_script_claim`: the card `\stagefield{Verification}` says
"Mathematica audit: none yet" despite a passing `.wl` (244, 246, 249, 250, 252).
Stale STATUS annotation, NOT a value/identity mismatch → **DEFERRED to
PAPER_CLEANUP P4-51** per the standing user decision; stages still `verified`.

## Numbering

- **Cards CLEAN** of the `+17 \stagefield{Purpose}` class (would be 260–270, all
  beyond 253 → absent).
- The **250/253 filenames legitimately embed** `…from_the_stage248_event_chain…`
  / `…from_the_stage252_export_and_cold_survival…` — BUILT-FROM cross-references
  in the descriptive title, NOT numbering errors.

## Infrastructure

- **22 DIRECT exec runs** (11 `python3` + 11 `math -script`, ≤2 mma concurrent in
  pairs, NOT parallel `$RT exec-*`) — **all exit 0, FAIL=0, deterministic, ALL
  byte-identical to committed** → ZERO output refreshes needed.
- 11 captured diffs all 0 bytes → `material_change: false` confirmed mechanically.
- 6 prose trackers synced (PAPER_CLEANUP **P5-22** = ZERO substantive paper/notes
  edits; card-text-lag deferral appended to P4-51).

## Per-stage verdicts

| Stage | Verdict | Note |
|-------|---------|------|
| 243 | verified (ckpt) | bar MET; `.wl` IBP-closure + Solve∧LinearSolve agreement + Series route |
| 244 | verified (finding) | card-text-lag → P4-51; 128√2 holds; F1 self-test-trap fix load-bearing |
| 245 | verified (clean) | `.wl` independent (hand-Hessian + Solve); F1 positive-control fix load-bearing |
| 246 | verified (finding) | card-text-lag → P4-51; round-trip backed by independent `MinValue` |
| 247 | verified (clean) | 142.17750000 holds, 210 purged; card clean/abstract; tautology fixes hold |
| 248 | verified (ckpt) | bar MET; satisfaction route non-vacuous (FreeQ guard); zero surviving 168 |
| 249 | verified (finding) | card-text-lag → P4-51; Möbius inverses genuinely Solve-derived |
| 250 | verified (finding) | card-text-lag → P4-51; global goldilocks via `Resolve[ForAll]` (not point) |
| 251 | verified (clean) | `.wl` Resolve[ForAll] + diff-quotient certificate + 50-digit NSolve |
| 252 | verified (finding) | card-text-lag → P4-51; X−X resolved to can-fail identity |
| 253 | verified (ckpt) | bar MET (FINAL); `.wl` Solves balance + `D[Log[V]]`; 119.23361317476524 holds |

## Cumulative — SECOND PASS COMPLETE

Pass-2 = **253/253 verified (100%)**. VIII.1 is the last batch; closing it
completes the full end-to-end second pass. Three consecutive zero-script-
correction batches (VII.1, VII.2, VIII.1) close it out — the first-pass retrofit
work held everywhere it was re-checked. Pass-1 `MANIFEST.yaml` untouched.

**Remaining gated follow-ups (all await the user):** the post-253 stem-keyed
numbering reconciliation (Phase 1 deterministic done), the dual-engine-required
backfill (121/122/123), and the script/output-band numbering plan.
