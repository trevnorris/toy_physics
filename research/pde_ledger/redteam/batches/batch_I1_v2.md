---
batch_id: I.1-v2
label: Part I.1 — Geometry lift, BdG coupling, projected Maxwell setup (paper-grounded re-audit)
started: 2026-05-25
completed: 2026-05-25
stages_total: 12
stages_verified: 12
stages_blocked: 0
material_change_any: true
---

# Batch I.1 v2 — Final summary (paper-grounded re-audit)

This is the **first batch processed under the v2 paper-grounded auditor**. The auditor now reads `paper/stages/stage_NNN.tex`, the per-stage notes under `notes/stages/`, and the part-level appendix BEFORE opening the scripts. A 10th finding category `paper_misalignment` was added with subtypes (`target_mismatch`, `value_mismatch`, `script_missing_paper_claim`, `paper_missing_script_claim`, `notes_contradicts_script`). Findings of that category route to user resolution rather than codex — codex never edits paper/notes.

All 12 stages reached `verified` status. 2 stages flagged `material_change: true` (001 sign flips, 004 Faraday/Bianchi restructure).

## Per-stage outcomes

| Stage | Audit findings | Resolution path | Verifier verdict | material_change |
|---|---|---|---|---|
| 001 | 2 paper_misalignment (Q1 source sign, Q2 gauge-fix sign) | Codex apply per resolutions (a)+(b) | verified | **true** |
| 002 | 0 (clean) | n/a | verified | false |
| 003 | 0 (clean) | n/a | verified | false |
| 004 | 1 tautological_check (Faraday/Bianchi block) | fix_loop + orchestrator inline-F fix | verified | **true** |
| 005 | 0 (clean) | n/a | verified | false |
| 006 | 1 paper_misalignment (Q3 gauge_μ placeholder) | Codex apply per resolution (c) | verified | false |
| 007 | 1 paper_misalignment (Q4 xi_eff^proj missing) + 2 script-side (obviated post-Q4) | Codex apply per resolution (a) + v2 re-audit | verified | false |
| 008 | 1 tautological_check (Integrate[W*Z] self-cancel) | fix_loop | verified | false |
| 009 | 0 (clean) | n/a | verified | false |
| 010 | 2 paper_misalignment (Q5 δu_n vs δP_n, Q6 7 unanchored clusters) | Codex apply per resolutions (c)+(a) | verified | false |
| 011 | 1 paper_misalignment (Q7 5 unanchored clusters) | Codex apply per resolution (c) — trim only | verified | false |
| 012 | 1 tautological_check (M1 block) | fix_loop | verified | false |

## Findings breakdown

- `paper_misalignment`: 7 (across stages 001, 006, 007, 010, 011)
- `tautological_check`: 4 (stages 004, 008, 012; one additional that didn't recur on 007 post-Q4)
- `insufficient_verification`: 0 in final v2 audit (one from initial 007 audit, obviated by H(w) extension)
- `stale_output`: 1 (stage 007 informational, resolved by output refresh)

V2-introduced findings missed by v1: **5 substantive** (007 F1 paper_misalignment, 010 F1+F2, 011 F1, plus the new script-side tautologies on 004, 008, 012). The 011 finding is particularly notable — v1 would have shipped 5 script clusters with no paper anchor; v2 caught that they should live in stages 022-024 instead.

## Resolution methodology

All 7 paper_misalignment items were routed to **Codex** as the math authority via two markdown-file sessions:
1. **Questions session**: `redteam/resolutions/batch_I1_paper_alignment.md` + `redteam/resolutions/codex_prompt_batch_I1.md`. Codex read paper, notes, scripts, and downstream stages, then filled `## Recommendation` blocks with `direction: a|b|c|skip` and rationales (cited file:line). Codex answered all 7, skipped 0.
2. **Apply session**: `redteam/resolutions/codex_apply_batch_I1.md`. Codex re-read its own recommendations, double-checked, applied edits across scripts and paper (authorized scope was paper/stages/, paper/appendices/, scripts/, mathematica/; NOT notes/ or redteam/), ran scripts to verify exit 0. All 7 applied (`0 revised, 0 blocked`).

User reviewed Codex's recommendations and concurred on all 7 before the apply session was launched.

## New toolchain pitfall documented

Added pitfall #5 to `prompts/codex.md`: **`Part[]`-indexing on pattern parameters inside `Do[Module[{locals = list}, ...]]` silently drops half the body.** Discovered on stage 004 when Codex's fix for the Faraday/Bianchi tautology used `fieldStrength[mu_, nu_] := D[potentialList[[nu+1]], ...]` inside a `Do[Module[...]]` loop. Mathematica's evaluator runs the Module body once during analysis (before `{alpha, beta, gamma} = triple` fires), passing the gensym'd locals to `fieldStrength`. Even `_Integer` pattern guards are insufficient — the partial evaluation leaks a malformed expression. Fix: precompute all needed `F_ij` as immediate-valued expressions BEFORE any Do/Module scope opens, then iterate over precomputed `(label, expr)` pair lists. Documented with full before/after code blocks.

## Process observations

1. **First batch under v2 auditor.** The new "paper first" reading order surfaced 7 paper_misalignment findings that v1 never could have caught — v1 was explicitly forbidden from reading paper/notes ("doc alignment is out of scope"). These are real defects: 011 alone had 5 script-verified clusters with no paper anchor, all of which actually live in stages 022-024.

2. **Codex as math authority pattern worked.** Routing 7 paper-vs-script disagreements to Codex via a structured questions-markdown file and a 3-option (a/b/c/skip) format produced 7/7 answers with citations. Codex's recommendations included cross-file evidence the auditor didn't have (notation firewall for metric signature, stage 008 already handling H(w) for stage 007, stages 022-024 publishing stage 011's clusters). User-as-reviewer + Codex-as-authority + Claude-as-orchestrator is a clean separation of concerns.

3. **Two acceptable controlled deviations** during apply (carried from earlier in batch I.1 v1; not new this v2 sweep).

4. **One mid-batch orchestrator hot-fix** (stage 004 Mathematica Part-indexing). Codex tried 3 iterations without solving it; the orchestrator applied the inline-F precomputation fix and documented the pattern in codex.md so future codex sessions know. Pattern is the same kind of "Mathematica gotcha" as the `ConditionalExpression` strip and `Limit` infinity test from III.2.

5. **Exec_logs and post-edit `.txt` outputs** were stale on the paper-alignment apply stages (001, 006, 007, 010, 011) because codex ran scripts but didn't persist outputs to the canonical paths. Manually refreshed all 5 sympy + mathematica outputs serially (Mathematica single-seat rule respected). Fix_loop stages (004, 008, 012) had fresh outputs (the fix_loop's "refresh sympy/mathematica" step handles this).

6. **`paper_misalignment` resolution loop is slower than fix_loop.** Each paper_misalignment requires: auditor flags it → orchestrator surfaces question → Codex recommends → user reviews → Codex applies → verifier confirms. Five separate sessions (audit, recommend, review, apply, verify) vs the fix_loop's two (codex applies + verifier). Worth tracking the per-stage time delta as more batches go through.

7. **Material change cascade**. Stage 001's sign flips affect downstream stages 003 + 005-009 in principle, but those stages were all just verified in this same v2 sweep against the corrected stage 001, so no separate upstream-stale check needed. Stage 004's Bianchi restructure has no downstream-quoted constant — no cascade.

## Tracker updates landed in this commit

- `notes/MATHEMATICA_MIRROR_POLICY.md`: snapshot bumped to v2 sweep close; no Independent-Mirror Set additions (the 5 paper-alignment apply stages didn't change Mathematica derivation routes, only added content).
- `notes/CHECKPOINT_TRUST_AUDIT.md`: snapshot bumped; no checkpoint stages in I.1.
- `notes/CHECKPOINT_CONSTANT_PROVENANCE.md`: snapshot bumped; no new constants surfaced.
- `notes/PAPER_CLEANUP_TRACKER.md`: new P4-38 row for batch I.1 v2 sweep; change log entry.
- `notes/EM_PROJECTED_INTEGRATION_TRACKER.md`: v2 paper-grounded re-audit row added (stages 004-021 are in scope for this tracker; v2 sweep confirmed alignment).
- `notes/STAGE_VERIFICATION_COVERAGE.md`: snapshot bumped to "batch I.1 v2 close"; cumulative count unchanged (84/253 — same stages, deeper verification).
