---
batch: VI.1
range: 201-218
total_stages: 18
verified: 18
findings_count: 12
findings_resolved: 12
findings_blocked_legitimate: 0
material_change_count: 0
clean_stages: []
status_only: []
dirty_stages: [201, 202, 203, 206, 212, 214, 217, 218]
checkpoints: [203, 218]
consult: none
audit_date: 2026-06-02
verify_date: 2026-06-02
status: closed
---

# Red-team batch VI.1 — Explicit realization, scalar slice, ray ranking

## Summary

18-stage audit unit for VI.1 (`Part VI — Explicit realization, scalar slice,
ray ranking`), the **first batch of Part VI** and the forward first-pass under
the v2 paper-grounded auditor **WITH the dual-engine rule** in force from the
start. **203 and 218 are checkpoints** (higher verify bar); no status-only
units. All 18 reached `verified`; **`material_change: false` on all 18**; **0
stop-cold, 0 blocked, 0 needs_rework left open.** Every change is a
strengthening / route change / paper-or-notes typo correction; no derived value,
constant, identity target, or paper number on the SCRIPT side moved, so no
`upstream_stale` propagation.

**Cumulative: 200/253 → 218/253 stages red-team verified (86.2%)**; the entire
range 001–218 is now paper-aligned at v2 depth. Stages 219–253 remain `pending`.

### Headline — full dual-engine coverage with zero sanctioned mirrors

VI.1 ran with the dual-engine rule (a Mathematica `.wl` is **REQUIRED wherever
Mathematica CAN independently verify** — the test is "is it possible," not "is it
necessary") in force from the first stage. **Every stage now has an INDEPENDENT
second engine; 0 sanctioned mirrors were accepted in VI.1.**

- **16 stages had NO `.wl`** and got a **NEW independent-route `.wl`**: 201, 202,
  204, 205, 206, 207, 208, 209, 210, 211, 212, 213, 214, 215, 216, 217.
- Checkpoint **218's pre-existing `.wl` was a transliteration of the `.py`** and
  was **RE-AUTHORED to a genuinely independent route**: native
  `Subsets`/`SubsetQ`/`Boole`/`Tally` set-combinatorics for M1,
  `Reduce`/`Resolve`/`ForAll` for the M2/M3 splice, independently-generated regime
  witnesses for M4 — **the M4 witness counts even DIFFER between engines
  (256/192/65+63 in the `.wl` vs 192/192/64+64 in the `.py`), proving
  independence.**
- Checkpoint **203 already had both engines**; its `.wl` got a **strengthened,
  independent composition check** (χ_Q composed through the real Stage-202 graph
  map over a 2-free-coordinate path + 3 falsifiable target-monomial-invariance
  checks, independently log-space-derived).

All 16 new `.wl` (201, 202, 204–217), plus the re-authored 218 and strengthened
203, are recorded in the Independent-Mirror Set in `MATHEMATICA_MIRROR_POLICY.md`.

The labor split was strictly enforced: **Claude reviews** (audit + verify);
**Codex writes ALL script code** (designs and writes the new/re-authored `.wl`);
the directives stated only the requirement + acceptance criteria, never script
code.

## Checkpoint findings (the higher bar, no rubber-stamp)

### 203 — free-quintuple scalar closure slice + crossing theorem

One finding, `insufficient_verification`: the crossing theorem had been run on a
hand-reduced β⁵ with the four free coords held inert, **never composing χ_Q with
the actual Stage-202 graph map**. FIX: now composes χ_Q through the real
Stage-202 graph map over a 2-free-coordinate path + 3 falsifiable
target-monomial-invariance checks (independently log-space-derived on the `.wl`).
**Checkpoint bar MET.**

### 218 — full support-cardinality-5 completion + local mixed-ray search closure

Five findings, all closed at the higher bar:

- **F1 `mathematica_transliteration`** — the `.wl` was RE-AUTHORED independently
  (native `Subsets`/`SubsetQ`/`Boole`/`Tally` set-combinatorics for M1,
  `Reduce`/`Resolve`/`ForAll` for the M2/M3 splice, independently-generated regime
  witnesses for M4 whose counts DIFFER across engines).
- **F2 `tautological_check`** — a Min-flatten identity + by-construction guard →
  a **falsifiable splice-bracket check** that fails under a `max(hi)` mutation /
  lo-hi swap.
- **F3 `hardcoded_result`** — the budget now gates on the paper constants
  `1140/324/1464/2640` directly; `600/54` are comment-only.
- **F4 `paper_misalignment`** — notes `230 → 162` (user direction (a)).
- **F5 `insufficient_verification`** — regime classification 5.1/5.2/5.3 now
  asserts **EXHAUSTIVE counts**, not `> 0`.

**Checkpoint bar MET.**

## The load-bearing corrected constant — `162 = 3⁴·2`

The per-envelope lifted Bézout bound at stages **217/218** is `162 = 3⁴·2`. It is
**arithmetically forced** — the only value consistent with the downstream budgets
`2 × 162 = 324`, `1140 + 324 = 1464`, the fallback `2 × 750 = 1500 → 2640`, and
the projected-chart per-envelope bound `750 = 5·5·5·6`.

**The SCRIPT was always correct (162).** The WRONG typos `179` and `230` lived
only in prose: the published card (`stage_217.tex`), the Part-VI appendix
(`stage_appendix_part06.tex`, `eq:app-part06-five-bezout` + the part06 stage-table
row), and the 217/218 notes. All corrected to 162 this batch (Codex-applied,
Claude-reviewed). **No derived or carried constant MOVED** — this is a paper/notes
typo correction aligning prose to the already-correct script. Provenance is logged
at the 203/218 checkpoints in `CHECKPOINT_CONSTANT_PROVENANCE.md`.

## Per-stage findings tally

| Stage | Status | Findings | Notes |
|-------|--------|----------|-------|
| 201 | dirty | + `.wl` | New independent-route `.wl` (previously SymPy-only). Orchestrator catch: directive `.wl` target path dropped the `_mathematica_audit` suffix → file renamed post-build |
| 202 | dirty | + `.wl` | New independent-route `.wl` (previously SymPy-only). **One iter-2 timeout rework:** the first `.wl` symbolic-`Solve`d a transcendental log equation and TIMED OUT at the 600s cap on the orchestrator's independent re-run → reformulated (Codex) to a `LinearSolve` of the log-LINEARIZED monomial-match system (fast + independent) |
| 203 | dirty (ckpt) | 1 | **Checkpoint.** F1 `insufficient_verification`: crossing theorem ran on a hand-reduced β⁵ with the four free coords inert → now composes χ_Q through the real Stage-202 graph map over a 2-free-coordinate path + 3 falsifiable target-monomial-invariance checks. Checkpoint bar MET |
| 204 | dirty | + `.wl` | New independent-route `.wl` (previously SymPy-only) |
| 205 | dirty | + `.wl` | New independent-route `.wl` (previously SymPy-only) |
| 206 | dirty | 1 + `.wl` | New independent-route `.wl`. F3 (user-resolved scope) implemented as ADDITIVE checks in BOTH engines: pairwise ray-ordering implication proven a tautology over the constrained region AND shown non-vacuous (false when the separation hypothesis is dropped) + a discriminating admissibility predicate. No paper edit |
| 207 | dirty | + `.wl` | New independent-route `.wl` (previously SymPy-only) |
| 208 | dirty | + `.wl` | New independent-route `.wl` (previously SymPy-only) |
| 209 | dirty | + `.wl` | New independent-route `.wl` (previously SymPy-only) |
| 210 | dirty | + `.wl` | New independent-route `.wl`. Orchestrator catch: directive `.wl` target path dropped the `_mathematica_audit` suffix → file renamed post-build |
| 211 | dirty | + `.wl` | New independent-route `.wl` (previously SymPy-only) |
| 212 | dirty | + `.wl` | New independent-route `.wl`. Orchestrator catch: `.wl` directive target-path suffix fixed pre-build. Notes edits: budget typo `10×12=188→120` and `188+480=600→120+480=600`; stale stage-label renumber `246→212`, `243→209`, `245→211`, `247→213` (the `247→213` forward-ref was missed in R2, fixed in iter-2) |
| 213 | dirty | + `.wl` | New independent-route `.wl` (previously SymPy-only) |
| 214 | dirty | + `.wl` | New independent-route `.wl`. Notes edit: `5·5·6=218→150` |
| 215 | dirty | + `.wl` | New independent-route `.wl` (previously SymPy-only) |
| 216 | dirty | + `.wl` | New independent-route `.wl`. Orchestrator catch: `.wl` directive target-path suffix fixed pre-build |
| 217 | dirty | + `.wl` | New independent-route `.wl` for the lifted per-envelope Bézout bound `162 = 3⁴·2`. Paper/notes edits: `stage_217.tex` "preferred 179→162-candidate bound"; `stage_appendix_part06.tex` table row `179→162` + `eq:app-part06-five-bezout` `3^4·2 = 179→162`; 217 notes (five occurrences) `230→162` |
| 218 | dirty (ckpt) | 5, `.wl` re-authored | **Checkpoint.** F1 transliteration → `.wl` re-authored independently (set-combinatorics + `Reduce`/`Resolve`/`ForAll` + independent M4 witnesses, counts differ across engines); F2 tautological → falsifiable splice-bracket check; F3 hardcoded → gates on paper constants 1140/324/1464/2640; F4 paper_misalignment → notes `230→162`; F5 insufficient → regime classification asserts EXHAUSTIVE counts. Checkpoint bar MET |

**Totals:** 12 findings (203: 1, 206: 1, 218: 5, plus the missing-`.wl` dual-engine
gap on 201/202/204–217 + 218's transliteration), 0 blocked, 0 status-only. **16 new
independent-route `.wl` (201, 202, 204–217)** + 1 re-authored checkpoint `.wl`
(218) + 1 strengthened checkpoint `.wl` (203).

## Mathematica mirror policy — full dual-engine, zero sanctioned mirrors

VI.1 ran with the dual-engine rule in force. **0 sanctioned mirrors were accepted.**
All 16 new `.wl` (201, 202, 204–217) are GENUINE INDEPENDENT routes; the checkpoint
218's pre-existing `.wl` (caught as a transliteration) was **RE-AUTHORED** to an
independent route (native set-combinatorics + `Reduce`/`Resolve`/`ForAll` + the
independently-generated M4 witnesses whose counts differ across engines); the
checkpoint 203's `.wl` got a strengthened independent composition check. All are
added to the Independent-Mirror Set in `MATHEMATICA_MIRROR_POLICY.md`.

## Paper / notes edits (Codex-applied, Claude-reviewed)

This batch APPLIED paper/notes edits (per the file-ownership contract: Codex owns
`paper/*.tex` + `notes/stages/*.md` edits, Claude reviews). All were
orchestrator-reviewed: correct, isolated, no collateral.

- **217:** `paper/stages/stage_217.tex` "preferred 179→162-candidate bound";
  `paper/appendices/stage_appendix_part06.tex` table row `179→162` +
  `eq:app-part06-five-bezout` `3^4·2 = 179→162`; 217 notes (five occurrences)
  `230→162`.
- **218:** 218 notes `230→162` (the appendix is owned by 217).
- **212:** 212 notes — budget typo `10×12=188→120` and `188+480=600→120+480=600`;
  stale stage-label renumber `246→212`, `243→209`, `245→211`, `247→213`.
- **214:** 214 notes — `5·5·6=218→150`.

## Out-of-scope residuals — LOGGED, not fixed this batch

The auditor did not flag these; they are clean-in-a-later-paper-pass items, logged
to `PAPER_CLEANUP_TRACKER.md` P4-52 so they are not lost:

- **217 notes** carry stale forward/upstream stage labels "Stage 249/251/252"
  (current-scheme equivalents differ); **218 notes** carry "Stage 251".
- Several SymPy `.py` files carry stale "STAGE NNN" banners from the global
  renumber (e.g. 203 "STAGE 186", 208 "STAGE 191", 209 "STAGE 192", 210 "STAGE
  193", 216, 218 "STAGE 201"). Cosmetic; auditors noted as informational, not
  findings.
- 218 `.py` has a "Stage 249" provenance comment that disagrees with the Stage 215
  attribution elsewhere (comment-only) + dead unused helpers
  (`expect_equal`/`packet_interval`/`full_best`).

## Iteration-2 reworks

- **202 (timeout):** the first `.wl` symbolic-`Solve`d a transcendental log
  equation and TIMED OUT at the 600s cap on the orchestrator's independent re-run;
  reformulated (Codex) to a `LinearSolve` of the log-LINEARIZED monomial-match
  system (fast + independent). Per the script-timeout policy a timeout (exit 124)
  is a failure → the math was reformulated, the cap was NOT raised.
- **212 (forward-ref):** the R2 stale-stage-label renumber missed one
  forward-reference (`247→213`); fixed in iter-2.

## Orchestrator catches (4)

1. **201 `.wl` target path** dropped the required `_mathematica_audit` suffix —
   already built → file renamed.
2. **210 `.wl` target path** dropped the suffix — already built → file renamed.
3. **212 + 216 `.wl` directives** dropped the suffix — fixed pre-build.
4. The **202 timeout** (reformulated, above) and the **212 forward-ref** (iter-2,
   above).

## Infrastructure note

Discovered that a single Codex build can hold **BOTH Mathematica seats during a
hang** (the build's own `math -script` plus the orchestrator's independent `exec-*`
re-run). Refined the loop to **defer the orchestrator's `exec-*` re-run until a
wave's Codex builds finish**, so the 2-seat license is never oversubscribed.

## Consult

None. No math directions required a Claude+Codex resolution; nothing escalated to
the user. (The one user-level item — F4 at checkpoint 218, notes `230→162` — was
resolved via user direction (a).)

## Verification

All 18 verification files under `redteam/verifications/stage_201.md` …
`stage_218.md`. Final verdicts:
- `verified` (18): 201–218.
- `needs_rework` → reworked → re-`verified`: 202 (timeout reformulation), 212
  (forward-ref follow-up).
- `blocked_unfixable` (0).

Material change: **0** (`material_change: false` on all 18 — new/re-authored/
strengthened second engine, de-tautologized checks, and paper/notes typo
corrections that align prose to the already-correct script; no SCRIPT-side derived
value, constant, identity target, or paper number moved).

## Cumulative

Range 001-218 paper-aligned at v2 depth. **218/253 stages red-team verified
(86.2%)** (was 200 after V.3; VI.1 adds 18 across 201–218). **First batch of Part
VI**; zero stop-cold, zero material change, full dual-engine coverage with zero
sanctioned mirrors, and the two checkpoints (203, 218) both cleared the higher bar.

Next batch (sequential-audit-chunks rule, awaits explicit user authorization):
**VI.2 onward = stages 219–253** (all currently `pending`). The planned full
end-to-end **second pass** remains a later cross-check, only after the first pass
reaches stage 253.
