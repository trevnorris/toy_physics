---
unit_id: 213
batch: VI.1
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-06-09T00:00:00Z
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 2
findings_total: 2
material_change: false
---

# Verification — unit 213

Both findings in the original report are non-script (F1 = paper-side card-text lag — "Mathematica audit: none yet" is now false; F2 = stale committed SymPy `.txt` carrying the pre-renumber "STAGE 196" banner). No Codex source edit was directed and none was needed: the directive's front-matter is `applied: false / needs_user_resolution: true`, and there is no `stage_213_diff.patch` (confirmed absent), so nothing was touched in `scripts/`, `paper/`, or `notes/`. Verification confirms (A) the output refresh landed clean (SymPy banner now "STAGE 213"), (B) the audit disposition still holds on the refreshed artifacts, (C) F1 is correctly classed as paper-prose lag routed to user (P4-51), and (D) `material_change: false`.

## Per-finding outcomes

### F1 — paper_misalignment (card-text lag: "Mathematica audit: none yet")

**Classification:** resolved (correctly routed — non-script, USER-DEFERRED to P4-51)

**What changed:**
Nothing in scripts. The directive (`directives/stage_213.md`, `## F1`, "Resolve before fix_loop") routes `paper/stages/stage_213.tex:11` to the user — Codex must not edit `paper/`. The card's `\stagefield{Verification}` line still reads "Mathematica audit: none yet", but a substantive, independently-derived `.wl` is present and passes. Correctly classified as a paper-prose coverage lag with no math impact, deferred to P4-51.

**Assessment:**
Legitimately deferred. The `.wl` exists and passes (Mathematica output lines 19–178 all PASS, no FAIL), so the card text understates coverage but contradicts no result. Non-blocking and outside the scripts-only scope. The disposition is correct: this is card-text lag from the first-pass dual-engine retrofit, not a value mismatch (reconciliation is 12/12 MATCH; the lone discrepancy is this prose line).

### F2 — stale_output (committed SymPy `.txt` carried pre-renumber "STAGE 196" banner)

**Classification:** resolved

**What changed:**
No source change. The orchestrator's independent re-run regenerated the committed outputs. `scripts/output/...stage213..._sympy_audit.txt` now reads `STAGE 213 — FOUR-COORDINATE MIXED SIMPLEX AND THE SUPPORT-CARDINALITY-4 GATE` (line 3) and `STAGE 213 SYMPY AUDIT COMPLETED SUCCESSFULLY` (line 228); the stale "STAGE 196" banner (the known +17 offset, 196+17=213) is gone. mtime is Jun 9 16:51, newer than the `.py` (Jun 3 15:59). The Mathematica output was already fresh and was likewise refreshed (mtime Jun 9 16:51).

**Assessment:**
Correct and complete; the refresh is the prescribed remedy. Every SymPy result line reduces to `0`/`True` (combinatorial ledger lines 19–34; face norm/slope reductions lines 119–126; gradient-optimal norm/slope/ratios + synergy gaps lines 159–167; `w_Σ` identity + Cauchy slack + 3/2/1 lines 174–181; diagonal-neutral + discriminant + ratio-tau + face collapses lines 192–212; the five brute-force gate theorems lines 226–230). No regression introduced.

## Disposition re-confirmation (post-refresh)

- **Genuinely independent `.wl`:** confirmed on the refreshed artifacts. The three load-bearing objects are DERIVED, not posited: (1) gradient ray via Lagrange `Solve` of the KKT system — Mathematica output line 48 reports `M3 positive Lagrange point = {ki/Sqrt[...], ...}` computed from stationarity, vs the `.py` posited `avec_grad`; (2) leverage bound via `Maximize` over the constrained nonneg sphere — line 106 `M5 Maximize result = {3, {x1->1/2,...,x4->1/2}}` returns both value and maximizer, vs the `.py` posited barycenter; (3) ten discriminant coefficients by `Coefficient`-extraction + matrix-kernel comparison — line 129 extracted set matches the boxed A..J (lines 130–151 all PASS), vs the `.py` posited scalar literals. Different operations on the load-bearing objects → INDEPENDENT, not a port. This is the V.3-style independence re-check the batch requires.
- **Gate theorems non-vacuous:** confirmed. The five interval-gate theorems are exhaustive integer sampling (1500625 ordered samples for the boundary splice; 924 each for screen dominance, non-improvement, support-4 improvement, support-4 non-improvement — SymPy output lines 226–230) — real can-fail brute-force checks, not trivially-true identities.
- **Both engines agree:** the combinatorial ledger (5 quadruples, incidence 2/4), face reductions, gradient ray and `||k||`, ratios, synergy gaps `k_missing²`, `w_Σ` max 3 at barycenter (2/1 at triple/pair), diagonal-neutral reduction, ten coeffs A..J, ratio tau bracket, and all face collapses land identically in both transcripts. `Maximize` returns exactly `{3,{1/2,1/2,1/2,1/2}}` matching the SymPy posited barycenter value 3.
- **0 reconciliation misalignments:** confirmed — 12/12 deliverable values MATCH per the report's reconciliation table; the only discrepancy is the F1 verification-prose line (a coverage statement, not a result value).

## Exec log assessment

**SymPy:** exit=0. Notable lines: `#quadruples - binomial(5,4) = 0` (L19); `gradient-optimal slope value = 0` (L160); `Kgrad^2 - Kijk^2 - k_l^2 = 0` (L164); `w_Sigma(four-way equal mix) - 3 = 0` (L179); `discriminant numerator reduction = 0` (L206); `verified four-face boundary splice theorem on 1500625 ordered integer samples` (L226). All equalities hold to zero; `STAGE 213 SYMPY AUDIT COMPLETED SUCCESSFULLY` (L233), `# exit_code: 0` (L235).

**Mathematica:** exit=0. Notable lines: `PASS: M3 positive Lagrange normalization` (L51); `PASS: M3 maximum value from positive branch minus ||k||` (L53); `M5 Maximize result = {3, {x1 -> 1/2, x2 -> 1/2, x3 -> 1/2, x4 -> 1/2}}` (L106) → `PASS: M5 constrained maximum value - 3` (L112); `PASS: M7 coefficient A matches boxed set` (L131) through `PASS: M7 reconstructed extracted Delta sharp equals numerator` (L151); `All Stage 213 Mathematica audit checks passed.` (L178), `# exit_code: 0` (L179). Every check prints PASS; no FAIL.

**Output freshness:** confirmed. Both `.txt` outputs carry mtime Jun 9 16:51, newer than the `.py` (Jun 3 15:59) and `.wl` (Jun 2 10:19). The SymPy banner is now canonical "STAGE 213" (the F2 stale-196 banner is cleared). No `stage_213_diff.patch` exists, consistent with no source edit having been directed.

## Material-change assessment

`material_change`: false. No source code changed; the only edits were regenerated committed `.txt` outputs (banner relabel 196→213 + transcript refresh). No derived result changed, so no downstream unit (> 213) is affected.

## Side observations (non-blocking)

None beyond the two reported findings. The (6) `S_quad` packet "partial" entry in the report is correctly a labeling artifact (a bookkeeping tuple; both interior screen points and the four-face interval Min are exercised) and is not a finding. I concur and add nothing.

## Verdict justification

`verified`. Both findings are non-script and resolved: F2's stale SymPy banner is cleared by the orchestrator re-run — the refreshed `.txt` carries the canonical "STAGE 213" banner with every SymPy check reducing to 0/True, and the Mathematica output prints PASS on every check with no FAIL (both exit 0); F1 is a paper-side card-text lag ("Mathematica audit: none yet") correctly USER-DEFERRED to P4-51. The audit disposition holds on the refreshed artifacts: the `.wl` is genuinely independent (Lagrange `Solve` ray, `Maximize` leverage bound, `Coefficient`-extracted discriminant coeffs — all absent from the posit-and-verify `.py`), the gate theorems are non-vacuous brute-force samples, and reconciliation is 12/12 MATCH. No diff patch (no source edit), no regressions; `material_change: false`.
