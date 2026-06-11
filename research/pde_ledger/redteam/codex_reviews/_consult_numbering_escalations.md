# Consult — SCRIPT/OUTPUT numbering-pass open escalations (Claude ↔ Codex)

**Date:** 2026-06-10
**Mode:** Claude (orchestrator) + Codex read-only consult (`codex_session_id: 019eb3de-ad92-7870-a23b-bfc8ac568e30`, model gpt-5.5, sandbox read-only, reasoning xhigh). Raw transcript: `redteam/codex_reviews/_consult_numbering_escalations.raw`. Prompt: `redteam/tmp_prompts/consult_numbering_escalations.md`.
**Trigger:** user asked Claude+Codex to jointly decide the OPEN, strictly-math/content-resolvable questions left by the (now committed+closed) script/output numbering pass, so they need not adjudicate math themselves.
**Method:** Claude posed 4 questions with a preliminary read; Codex opened the files itself and returned independent verdicts; Claude then re-verified every disputed claim against the `notes/` + card sources (the standing backstop — agents under-call AND can fabricate corroboration). **All 4 agreed after verification.**

## Iron rule (held throughout)
Canonical = filename `stageNNN_<stem>` / paper card. +17 (old PDE K→K+17) is corroboration only, never the decision rule. Content (which stage OWNS the cited deliverable) decides. A numbering-pass FIX must be **digit-only** (strip-digits proof); anything that changes a stem, removes a token, or renames a code identifier is OUT of the label-only pass and routes to its proper applier.

---

## Q1 — `stage182.py:45` `# Stage-30 coherent branch definitions.` → **FIX 30→047** (CONCUR; Codex disputed Claude's lean-LEAVE, verification upheld Codex)
- **Owner = stage047** (`coherent_kernel_map`). Decisive content (Claude-verified):
  - `notes/stages/...stage182_...md:83` — "**Stage 047** already expressed the coherent local D/N branch in terms of the microscopic couplings" (ε_η, χ₀, ε_W).
  - `notes/stages/...stage182_...md:138` — "Differentiating the **Stage-047 coherent ratios** gives the physical Stage 181 drift variables directly." (The .py block below the comment defines exactly those drift vars: zetaZ, chi1, eta1, varepsW, …)
  - `notes/stages/...stage047_coherent_kernel_map.md:92` — "Exact coherent dimensionless ratios" (eps_eta, eps_W).
- So "Stage-30" is a stale OLD-epoch cross-ref (30+17=47), NOT a loose pointer. Earlier left-as-ambiguous only because the orchestrator hadn't opened the notes.
- **Disposition:** digit-only comment cross-ref → in label-only scope → **orchestrator applies** `# Stage-30` → `# Stage-047` (3-digit cross-ref convention) in `stage182.py:45`; re-run sympy to refresh (comment → output byte-identical). Codex found NO `Stage-30` in the stage182 `.wl`.

## Q2 — `stage203.wl` `chiFromStage180`/`closureNumStage180` → **RENAME →197** (CONCUR, airtight)
- The `.wl` (lines 289/290/294/296) and the `.py` twin define **byte-identical** expressions; the `.py` names them `chi_from_stage197` / `closure_num_stage197` (`stage203.py:312/313`).
- **Owner = stage197** (`conditional_packetA_closure_theorem`): `notes/stages/...stage197_...md:229` owns `chi_Q = 3(S β^5 + 9 Σ5)/(3S − Σ0)`; `:239` owns the closure numerator. stage180 = `effective_transfer_shape_collapse` owns none of it. 197−17=180 → stale old number baked into a CODE IDENTIFIER.
- **Disposition:** a code-identifier rename (4 occ: `chiFromStage180`→`chiFromStage197` @289/294; `closureNumStage180`→`closureNumStage197` @290/296) is OUT of the label-only pass (zero-variable-bytes rule) → **Codex applies** + re-run (value-inert; output unchanged). Per [[feedback-claude-reviews-codex-codes]]/[[feedback-codex-is-fix-applier]], orchestrator does NOT hand-edit code identifiers.

## Q3 — `stage100` chi_Q family cites `stage 097 (single_normalization_defect)` → **LEAVE / flag for content review** (CONCUR-escalate)
- `grep chi_Q` in stage097 `.py`+notes is **EMPTY** — 097 has no χ_Q; it closes the single normalization defect (N_Q=1). The exact χ_Q=1 identification lives at **105** = `chiQ_fix_from_outgoing_dtn` (`stage105.py` solves the outgoing-DtN match; `notes:67` yields χ_Q=1).
- BUT the cited number AND its parenthetical stem are internally consistent (both name 097), and 097 IS the genuine *precursor defect* that χ_Q later fixes — so whether the citation is a mis-attribution (→105) or an intentional "the defect from 097" pointer is an **author-intent / semantic** call. Changing it would alter the stem, not just a digit.
- **Disposition:** NOT a numbering fix. Both engines concur: **LEAVE** for this pass; record as a candidate for a future content/red-team review of stage100/101's chi_Q provenance.
- **⮑ REOPENED + RESOLVED 2026-06-10** (user clarified there is no external author — Claude+Codex wrote the ledger): see `_consult_q3_chiQ_provenance.md`. Decided **FIX** — the attribution is a genuine error (097=branch context; **104**=exact fingerprint; **105**=χ_Q=1 fix, downstream; the "upstream 097" wording is self-contradictory + inconsistent with stage106). Also surfaced + fixed a **paper-card overclaim** in `stage_100.tex:13`/`stage_101.tex:13`. Applied across 6 files (scripts+notes+cards) + committed; math untouched, 3/3 re-run exit 0.

## Q4 — `stage105.py:8/31` compound `Stage 088/074` → **FIX (drop /074 → `Stage 088`)** (CONCUR; Codex disputed Claude's lean-LEAVE, verification upheld Codex)
- Quantity = pole scale `Ω_Q = 3c_s/(2a)` + the minimal isotropic conservative quadrupole module.
- **Owner = stage088** (`loading_ratio_from_minimal_module`): card `stage_088.tex:7` Inputs = the minimal precursor `Y_Q^cons = 3/4 + (1/4)/(1−ω²/Ω_Q²)`; `:24` **boxes** `Ω_Q = 3c_s/(2a)`. stage074 = `family1_healing_lock` (card Purpose :6 "healing-length lock"; outputs χ_s=37/2, κ=9/5 Λ², α — **no Ω_Q**). So "/074" is a wrong/stale co-citation.
- **Disposition:** removing the `/074` token is a content correction (NOT digit-only, so out of the label-only pass) → **Codex applies** `Stage 088/074` → `Stage 088` at `stage105.py:8,31` + re-run (comment → output unchanged).

---

## Joint outcome
| Q | Verdict | Owner | Edit | Applier | In numbering-pass scope? |
|---|---------|-------|------|---------|--------------------------|
| Q1 | FIX | 047 | `Stage-30`→`Stage-047` comment, stage182.py:45 | orchestrator | yes (digit-only) |
| Q2 | RENAME | 197 | `…Stage180`→`…Stage197` ×4, stage203.wl | Codex | no (code identifier) |
| Q3 | LEAVE/flag | (105 vs 097) | none | — | no (semantic attribution) |
| Q4 | FIX | 088 | drop `/074`→`Stage 088` ×2, stage105.py | Codex | no (content/token removal) |

All four are cosmetic (comments + one value-inert rename): **zero equation/value/logic change**; every script re-runs exit 0 with byte-identical numerical output. Q1/Q2/Q4 are safe to apply; Q3 is a deliberate leave.

## ✅ APPLIED 2026-06-10 (user-approved "apply all three"; committed in this commit)
- **Q1** — orchestrator-applied directly: `stage182.py:45` `# Stage-30` → `# Stage-047` (digit-only comment; strip-digits proof PASS).
- **Q2** — Codex-applied (`redteam/codex_logs/apply_numbering_consult_q2q4.txt`): `stage203.wl` 4 identifier occurrences `…Stage180`→`…Stage197` (RHS expressions untouched).
- **Q4** — Codex-applied (same log): `stage105.py:8,31` `Stage 088/074` → `Stage 088` (comment-text token removal; intentionally NOT digit-only).
- **Q3** — no edit (deliberate leave; flagged for a future content/red-team review of stage100/101 chi_Q provenance).
- **Orchestrator verification (the backstop):** `git diff` line-exact on all 3 files (no equation/value/literal/logic byte changed; residual scan for `Stage180`/`088/074`/`Stage-30` empty); re-ran sympy 182+105 + mathematica 203 → **3/3 exit 0**, outputs byte-identical (stage182's `free_symbols` SET-ordering re-run noise reverted as spurious, same as band 181-253); `git grep` for the wrong-root path typo empty. 3 source files changed, 0 output files changed.
