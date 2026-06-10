---
batch: VII.2
pass: 2
range: 231-242
total_stages: 12
verified: 12
findings_count: 5
findings_resolved: 0
findings_deferred: 5
findings_blocked_legitimate: 0
material_change_count: 0
script_corrections: 0
codex_invocations: 0
clean_stages: [232, 233, 234, 236, 238, 239, 242]
findings_stages: [231, 235, 237, 240, 241]
status_only: []
checkpoints: [239, 242]
ports_found: 0
audit_date: 2026-06-09
verify_date: 2026-06-09
status: closed
---

# Red-team PASS-2 batch VII.2 — Rigid-mouth orbit-lock / branch-dressing / twin-support (231–242)

## Headline

**The SECOND pass-2 batch needing ZERO script corrections** (after VII.1). No
`.py`/`.wl` edit anywhere; **0 Codex invocations.** All 12 stages reached
`verified`, **`material_change: false` on all 12** (working tree carried only new
report/directive artifacts — every captured diff was 0 bytes). Cumulative pass-2
= **242/253**; only VIII.1 (243–253) remains.

VII.2 was a first-pass **RETROFIT** batch (10 NEW independent-route `.wl` +
2 checkpoint `.wl` re-authored FROM transliterations), the same class as
V.3/VI.1/VII.1. The pass-2 mandate was therefore: **confirm both engines present
AND re-check `.wl` independence with the orchestrator ground-truth read** —
because a first-pass retrofit can be insufficient (V.3-200 proved it on a
checkpoint). This pass it was **sufficient everywhere**: 12/12 `.wl` confirmed
genuinely independent, **0 ports, 0 borderline.**

## Independence — EARNED, not relayed

Per the standing distrust of an all-clean result (the VII.1 precedent), the
orchestrator did NOT relay the audit agents' clean verdicts:

- **Orchestrator ground-truth `.wl`-vs-`.py` read IN FULL on BOTH checkpoints
  (239, 242)** — the highest-risk re-authors (where V.3-200 was caught).
- **Orchestrator spot-check of the two highest port-risk non-checkpoints
  (235, 236)** — the pure-linear-algebra stages where structural parallelism is
  hardest to distinguish from a port.
- The 12 per-stage audit agents each supplied operation-level, line-referenced
  independence evidence (different primitives, not relabels).

No port was flagged by any audit agent, so (as in VII.1) no CALIBRATION agent was
needed.

## Checkpoints (the higher bar) — both re-authors confirmed SUFFICIENT

The VI.1-218 / VII.1-221 outcome (sufficient), NOT V.3-200 (insufficient).

### 239 — rigid-mouth physical normal form / Cartesian orbit-lock
- `.py` **POSITS** the compiler matrix `S_rm_dep = [[0,0],[0,-1],[1,-1]]`
  (py:174) and the left-inverse `L_phys_dep = [[0,-1,1],[0,-1,0]]` (py:206) as
  literals, then checks consistency; orbit-lock via `sp.solve` (py:310-316).
- `.wl` **DERIVES** them: `compilerJacobian = Table[D[dependentDelta[[row]],
  physicalCoordinates[[col]]], …]` (forward Jacobian of the boxed dependent
  vector, wl:135-146); left-inverse via native `PseudoInverse[compilerJacobian]`
  (wl:150) cross-checked to reproduce `{DeltaMu-DeltaKeta, -DeltaKeta}`;
  orbit-lock via `Reduce[…]` + `Equivalent[orbitLockReduced, U==0 && V==0]`
  (wl:229-235).
- **Opposite information flow** (one DERIVES what the other POSITS) — not a
  transliteration. Banner canonical (`STAGE 239`; carried labels Stage238/236).
  **Checkpoint bar MET.**

### 242 — twin-support strict inclusion / coherent orbit-lock
- Load-bearing claim = strict twin-window inclusion `C_mix < Pi_tr < 2·C_mix`
  (the ρ_α-style 4/3 region).
- `.py` (py:94-97): `ratio = nsimplify(Pi_tr/C_mix)`; `assert ratio == 4/3`;
  `assert ratio > 1 and ratio < 2` — scalar rational compare, STRICT.
- `.wl` (wl:152-165): 4/3 independently DERIVED via
  `FullSimplify[traceLoad/mixedCapacity]`, then a `Resolve[ForAll,{lambdaWin,
  epsilonWin}, …, Reals]` strict-inequality QE certificate over the admissible
  family — STRICT (`1 < … < 2`). De-transliterated devices: abstract-`Function`
  → `D[closedObservablePacket, zeta]` on real closed forms; `Exp[t·d]` →
  `logDrift` total-log-differential.
- Genuinely different operations (scalar compare vs `Resolve` QE), both strict,
  each fails at the boundary. Banner canonical (`STAGE 242`, no `225` residue).
  **Checkpoint bar MET.**

## Reconfirmations — all 3 notes-only paper_misalignment fixes HOLD

Cross-engine-corroborated; **published cards/appendices UNAFFECTED** (notes-only):

- **231** — notes `dF/dξ` numerator coeffs `240·δ²ξ→189·δ²ξ`, `189·ξ³→121·ξ³`.
  Both engines emit 189/121; the `.wl` M1 output prints the numerator factor
  `… + 189*d^2*x + 297*d*x^2 + 121*x^3`. Notes match. HOLDS.
- **232** — figure-of-merit prefactor `168→100`. **ZERO surviving `168`** in any
  VII.2 script/note/card (arbiter grep clean; only incidental digit-substrings
  inside residual mantissas). Scripts emit `c=100`, reproducing the notes'
  decimals Ξ_χ≈5.5548e5 / Ξ_J≈1.2664e5. HOLDS. (The recurring stale-"168"
  family — cf. 148/232/248 — is purged here.)
- **241** — `ϱ_WΛ` bound `193/369→125/369`. Notes line 577 carries `125/369≈
  0.338753`; `.wl` M7 computes `ϱ_WΛ|_{β=2/11}=125/369`; the ϱ windows
  `1/3, 125/369, 2/3, 250/441` all reconcile. HOLDS.

## Reconfirmations — script-side fixes HOLD (the variable-independence self-test trap)

The dominant first-pass defect theme (`D[]`/`sp.diff` w.r.t. an ABSENT variable
→ vacuous ≡0 pass). All fixes present and load-bearing in BOTH engines:

- **237-F2** — live-channel negative control + `Not[FreeQ[#,Derivative]]` leak
  detector; every `D[]`/`diff` differentiates w.r.t. a variable that appears.
- **238-F1/F4** (the pass-1 iter-2 reframe) — negative control `∂_ζ M_tr ≠ 0`
  (genuine nonzero witnesses) + leak detector `∂_ζ(R_tr·M_tr/M_mix) ≠ 0` +
  structural exclusion of `{ζ, M_mix, M_tr}` from the reduced observables; would
  genuinely FAIL if support-blindness were faked.
- **240-F1** — weights extracted from `Y_support` (which carries `Omega_Q` via
  its pole, asserted by presence guards), exercising the independence claim
  against a variable-bearing object, not a constant; `ρ_α=4/3`, `ζ_req=1/3`,
  `Π_tr=(4/3)C_mix` reconcile.

## Findings (5) — all the standing-deferred CARD-TEXT-LAG class

5 `paper_misalignment` / `paper_missing_script_claim`: the card
`\stagefield{Verification}` says "Mathematica audit: none yet" despite a passing
`.wl`. Formalized on **231, 235, 237, 240, 241** (clean stages 232/234/236
noted it informally; 238/239/242 cards similar). This is a stale STATUS
annotation, NOT a value/identity mismatch — an inconsistent-threshold
formalization, the documented agent under-call. **DEFERRED to PAPER_CLEANUP
P4-51** per the standing user decision; stages still go `verified` (non-blocking,
later Codex-applied + Claude-reviewed paper pass).

## Numbering

- **Cards CLEAN** of the `+17 \stagefield{Purpose}` self-label class
  (+17 = 248–259 absent on 231–242).
- **233 `.py` stale forward-renumber comment/print labels** (Stage 239/240/241
  where the source is canonical Stage 188 / budgets 224 / compat-point 223) —
  the `.wl` is canonical; comment-only, touches no computed value →
  `NUMBERING_SCRIPT_OUTPUT_BAND_PLAN.md` (content-keyed, never offset-sweep).
- Notes-TITLE drift (232 "Stage 249", 234 "251", 235 "251/252/253", 236 "253")
  → deferred post-253 stem-keyed cleanup.

## Infrastructure

- **24 DIRECT exec runs** (12 `python3` + 12 `math -script`, ≤2 mma concurrent in
  pairs, NOT parallel `$RT exec-*` — MANIFEST race avoidance) — **all exit 0,
  FAIL=0, deterministic.** **All 24 byte-identical to committed** → ZERO output
  refreshes needed (cleaner than VII.1's 2 label-only refreshes; banners already
  canonical).
- 12 captured diffs all 0 bytes → `material_change: false` confirmed mechanically.
- Seat policy held (no `.wl`-touching Codex this batch; mma exec ≤2 concurrent).
- 6 prose trackers synced (PAPER_CLEANUP **P5-21** = ZERO substantive
  paper/notes edits; card-text-lag deferral appended to P4-51).

## Per-stage verdicts

| Stage | Verdict | Note |
|-------|---------|------|
| 231 | verified (finding) | card-text-lag → P4-51; 189/121 fix holds; `.wl` independent (Resolve[ForAll] + NSolve) |
| 232 | verified (clean) | 168→100 holds, zero surviving 168; `.wl` independent (native Integrate vs mp.quad) |
| 233 | verified (clean) | `.wl` independent (Coefficient vs subs); `.py` stale comment-labels → numbering band |
| 234 | verified (clean) | `.wl` independent (Series/Coefficient + Minimize vs diff+Lagrange); M4 cancellation real |
| 235 | verified (finding) | card-text-lag → P4-51; `.wl` independent (MatrixPower involution + Det + Solve-cardinality) |
| 236 | verified (clean) | `.wl` independent (left-inverse two ways: LinearSolve AND Solve, agree; +eigen/rank/Equivalent) |
| 237 | verified (finding) | card-text-lag → P4-51; self-test-trap guard present+load-bearing both engines |
| 238 | verified (clean) | iter-2 reframe holds (neg-control + leak-detector + exclusion, both engines) |
| 239 | verified (ckpt) | checkpoint bar MET; re-author SUFFICIENT (D[]-Jacobian + PseudoInverse + Reduce vs posited literals) |
| 240 | verified (finding) | card-text-lag → P4-51; F1 fix holds (weights from Omega_Q-bearing Y_support pole) |
| 241 | verified (finding) | card-text-lag → P4-51; 125/369 fix holds; `.wl` independent (Solve vs hardcoded targets) |
| 242 | verified (ckpt) | checkpoint bar MET; re-author SUFFICIENT (Resolve QE + FullSimplify-derived 4/3 vs scalar compare) |

## Cumulative

Pass-2 = **242/253** verified (was 230 after VII.1; VII.2 adds 12). Remaining:
**VIII.1 (243–253, 11 stages)** — awaits the explicit user gate
(sequential-audit-chunks). Pass-1 `MANIFEST.yaml` untouched.
