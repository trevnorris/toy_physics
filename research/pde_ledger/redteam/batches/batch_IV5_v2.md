---
batch: IV.5
range: 139-150
total_stages: 12
verified: 12
findings_count: 31
findings_resolved: 30
findings_blocked_legitimate: 1
material_change_count: 0
clean_stages: [141, 145, 149]
status_only: [141, 145, 149]
dirty_stages: [139, 140, 142, 143, 144, 146, 147, 148, 150]
checkpoints: []
audit_date: 2026-05-27
verify_date: 2026-05-27
status: closed
---

# Red-team batch IV.5 — Susceptibility, branch, defect transport

## Summary

12-stage audit unit for IV.5 (`Part IV.5 — Susceptibility, branch, defect transport`). 9 dual-engine stages (139, 140, 142, 143, 144, 146, 147, 148, 150) and 3 status-only stages (141, 145, 149). No checkpoints.

31 findings distributed across 9 dirty stages — 30 resolved + 1 blocked_legitimate (stage 144 F4 transliteration policy, accepted as written per user gate). 12 stages verified end-to-end; both engines exit 0; outputs fresh.

Zero material change to any derived constant or downstream-cited quantity. Eleventh consecutive zero-redirection batch.

## Per-stage findings tally

| Stage | Status | Findings | Engines | Notes |
|-------|--------|----------|---------|-------|
| 139 | dirty | 4 | SymPy + Mathematica | F1 insufficient (all 6 paper deliverables), F2 transliteration, F3 hardcoded provenance, F4 banner |
| 140 | dirty | 2 | SymPy + Mathematica | F1 transliteration+banner, F2 insufficient (boxed numerics) |
| 141 | clean | 0 | — | Status-only; verdict clean |
| 142 | dirty | 5 | SymPy + Mathematica | F1 tautological R_q(g_-), F2 insufficient (5 canonical-point anchors), F3 transliteration, F4 hardcoded, F5 banner (paper_misalignment subtype, resolved by Cluster A) |
| 143 | dirty | 3 | SymPy + Mathematica | F1 insufficient (8 paper deliverables), F2 hardcoded sInf=1, F3 transliteration |
| 144 | dirty | 4 | SymPy + Mathematica | F1 banner (Cluster A), F2 insufficient (7 numerical targets), F3 paper-card Checks downgrade (Cluster C), F4 transliteration policy (blocked_legitimate, user-accepted) |
| 145 | clean | 0 | — | Status-only; verdict clean |
| 146 | dirty | 3 | SymPy + Mathematica | F1 tautological affine-law, F2 missing symbolic moments, F3 banner |
| 147 | dirty | 3 | SymPy + Mathematica | F1 SymPy zero-asserts, F2 Mathematica tautological + banner, F3 moment-stability missing |
| 148 | dirty | 3 | SymPy + Mathematica | F1 tautological D4 / hardcoded xi_star + directive-typo catch, F2 subsumed, F3 transliteration + banner |
| 149 | clean | 0 | — | Status-only; verdict clean |
| 150 | dirty | 1 | SymPy + Mathematica | F1 tautological T_q'(0) - S_q |

**Totals:** 31 findings, 9 dirty stages.

## Three user-gate clusters (all `(Recommended)`)

### Cluster A — Mass renumbering (mechanical, 21 edits)

- 11 `.wl` banner edits off by −17 (139→122, 140→123, 142→125 + LEDGER, 143→126 + LEDGER, 144→127 + LEDGER, 146→129, 147→130, 148→131).
- 6 `.py` banner edits off by −17 (142, 143, 144 with LEDGER variants).
- 4 notes H1 edits off by +102 (146→248, 147→249, 148→250, 149→251).

### Cluster B — Body-text forward-stage citation re-attribution (22 edits)

22 pre-renumber citations across 11 of 12 notes. Two offsets: −51 for stages 188-199, −102 for stages 220-251. Re-attributed to current numbering. Notes are now self-consistent.

### Cluster C — Stage 144 paper-card Checks downgrade

Items (i) outlet-consistency and (ii) self-matched-susceptibility rewritten as carry-forward citations of `\ref{stage:135}` and `\ref{stage:140}`, mirroring IV.4's stage 134 downgrade pattern.

## Orchestrator catches (rework loop)

1. **Stage 148 directive-prescribed closed form was wrong.** Auditor copied `4107 - 168*pi^2` from stage 148 notes (which itself was a typo); correct form is `4107 - 100*pi^2` (confirmed against stage 126 upstream). Fixed in both engines + stage 148 notes typo also corrected.
2. **Stage 139 pitfall #13 recurrence**: `Pi_*)` substring in Mathematica comment parses as comment-terminator. Rewrote with ASCII-safe names.
3. **Stage 142 SymPy tolerance**: 1e-20 too tight for nsolve's 30-digit precision. Loosened to 1e-15.
4. **Stage 147 SymPy precision**: `sp.N(AT)` default-truncates to 15 digits; assigned `AT_30 = sp.N(AT, 30)` explicitly.
5. **Stage 146 SymPy `Integrate` fallbacks**: both F1 (eps-sample) and F2 (Pi-sample) needed numeric fallbacks mirroring the Mathematica path because `sp.simplify` / `sp.integrate` couldn't reduce the residuals symbolically.

## Verification

All 12 verification files written under `redteam/verifications/stage_*.md`. Final verdicts:

- `verified` (12): all stages
- `needs_rework` (0)
- `blocked_unfixable` (0)

## Cumulative

Range 001-150 paper-aligned at v2 depth. 162/253 stages red-team verified (64.0%). Eleventh consecutive zero-redirection batch.

Next batch (sequential-audit-chunks rule, awaits explicit user authorization): **IV.6 = stages 151-163** ("Correction, coevolution, traction, off-family"), 13 stages.
