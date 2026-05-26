---
batch_id: II.1-v2
created: 2026-05-25
auditor_version: v2 (paper-grounded)
total_paper_misalignment_items: 3
codex_role: math-authority + cross-stage research; provide direction (a/b/c/skip) + rationale per Q
user_role: review Codex recommendations, approve or redirect before any apply session
---

# Batch II.1 v2 — paper_misalignment consolidation

The v2 audit of stages 024–036 (Overlap isotropy through continuum kernel) produced **3 paper_misalignment items**: stage 029 F1, stage 029 F4, stage 035 F1. Script-side findings (covered separately via fix_loop) are NOT in this file.

For each item Codex must fill in the `## Recommendation` block with `direction:` set to (a), (b), (c), or skip, plus rationale citing the relevant file:line evidence. Codex MUST NOT apply any change in this session — recommendations only.

---

## Q1 — Stage 029 F1 (paper_misalignment / notes_contradicts_script, low)

**Item:** SymPy docstring lines 3 and 5 reference "Stage 12" / "moving_throat_pde_stage12_..."; Mathematica banner at line 33 prints `STAGE 012 — DYNAMIC LOADING`. Paper card title is "Stage 029: Coupled Response Operator" (`paper/stages/stage_029.tex:1`). The notes file `notes/stages/moving_throat_pde_stage029_dynamic_loading.md` mixes legacy "Stage 11/12" references with current "Stage 029".

**Audit finding location:** `redteam/reports/stage_029.md` § F1

**Context:** Pure labeling drift from an earlier numbering scheme; math content unaffected. The closing line of the Mathematica script (line 232) already correctly says "Stage 029 Mathematica audit passed." So the file knows both numbers.

This is borderline as a paper_misalignment finding — it's a script-side metadata error, not a paper↔script substantive disagreement. (Compare stage 022/028 where similar legacy-label drift was explicitly NOT raised as paper_misalignment, treated as cosmetic.) The auditor flagged it under paper_misalignment but noted "Codex-applicable script-side label fix" in its summary. So Codex needs to advise whether to:

**Options:**
- **(a) Re-categorize as script-side cosmetic** — treat as not requiring user gate; relabel docstring/banner under normal fix_loop without paper change. (This matches the I.1/I.2 v2 precedent for stage 022/028 legacy labels.)
- **(b) Relabel scripts AND propagate "Stage 029" into the notes file** to fully eliminate legacy numbering. (More work; notes are read-only for Codex but user could approve.)
- **(c) Acknowledge as known limitation** — leave the legacy labels intact (descriptive of the original numbering); add a one-line script comment noting the rename.
- **skip** — Codex defers entirely.

## Recommendation

direction: a
rationale: Stage 029 is the current paper stage (`paper/stages/stage_029.tex:1`), while the SymPy docstring still names `moving_throat_pde_stage12...` and "Stage 12" (`scripts/moving_throat_pde_stage029_dynamic_loading_sympy_audit.py:3`, `scripts/moving_throat_pde_stage029_dynamic_loading_sympy_audit.py:5`) and the Mathematica banner still prints `STAGE 012` (`mathematica/moving_throat_pde_stage029_dynamic_loading_mathematica_audit.wl:33`). The same Mathematica file already ends with the current-stage success message (`mathematica/moving_throat_pde_stage029_dynamic_loading_mathematica_audit.wl:232`), and the notes mix current title text with legacy Stage 11/12 prose (`notes/stages/moving_throat_pde_stage029_dynamic_loading.md:1`, `notes/stages/moving_throat_pde_stage029_dynamic_loading.md:5`, `notes/stages/moving_throat_pde_stage029_dynamic_loading.md:297`, `notes/stages/moving_throat_pde_stage029_dynamic_loading.md:303`). The precedent points the same way: Stage 022 called stale labels documentation drift rather than `paper_misalignment` (`redteam/reports/stage_022.md:129`), and Stage 028 explicitly treated Stage 11/028 drift as a renumbering artifact rather than a math/paper disagreement (`redteam/reports/stage_028.md:164`). Re-categorize this as script-side cosmetic and relabel the script identifiers under the normal fix loop, without paper or notes changes.

---

## Q2 — Stage 029 F4 (paper_misalignment / paper_missing_script_claim, low)

**Item:** Both scripts verify a closed form for the refined softening threshold

```
alpha_crit = K0t * K1t / (K1t * kappa_0^2 + K0t * kappa_1^2)
```

at `scripts/moving_throat_pde_stage029_dynamic_loading_sympy_audit.py:189-194` and `mathematica/moving_throat_pde_stage029_dynamic_loading_mathematica_audit.wl:159-166`. The notes file (`notes/stages/moving_throat_pde_stage029_dynamic_loading.md:173-186`) defines and discusses `alpha_crit` in §3.1. But the paper card (`paper/stages/stage_029.tex`) does NOT mention `alpha_crit` in either `\stagefield{Output}` (line 106) or `\stagefield{Checks}` (lines 108-113). The card's Output enumerates only Schur split, static branch data, β₀, and selected odd coefficient.

**Audit finding location:** `redteam/reports/stage_029.md` § F4

**Cross-stage research needed:** check downstream stages (030, 031, 032, 033, 034 — all consume α_crit per the v1 audits already reviewed) to determine which stage actually OWNS the α_crit derivation in the paper. Per "each stage builds on prior" principle, look up the chain in `paper/stages/stage_030.tex` through `stage_033.tex` for `\stagefield{Inputs}` lines that reference α_crit. Whichever paper card first introduces α_crit is its owner; the others import it.

**Options:**
- **(a) Add α_crit to stage 029 paper card** — expand `\stagefield{Output}` and `\stagefield{Checks}` to include the closed-form derivation. Notes-source authority preserved.
- **(b) Trim α_crit from stage 029 scripts** — remove the assertion if a downstream stage (030/031) actually owns the derivation in its paper card. Verify destination ownership via grep before recommending. Apply [[feedback-redteam-priorities]] "destination-verification guardrail".
- **(c) Acknowledge as known scope mismatch** — leave as-is; the script-side check is a consistency anchor that the paper card chooses not to highlight. Add `\stagefield{Notes}` clarification in the card.
- **skip**

## Recommendation

direction: b
rationale: Stage 029's paper card lists only the Schur split, static branch data, `beta_0`, and selected odd coefficient in Output and Checks (`paper/stages/stage_029.tex:106`, `paper/stages/stage_029.tex:108`), while the Stage 029 audits assert `alpha_crit` anyway (`scripts/moving_throat_pde_stage029_dynamic_loading_sympy_audit.py:191`, `scripts/moving_throat_pde_stage029_dynamic_loading_sympy_audit.py:194`, `mathematica/moving_throat_pde_stage029_dynamic_loading_mathematica_audit.wl:164`, `mathematica/moving_throat_pde_stage029_dynamic_loading_mathematica_audit.wl:166`). The downstream owner is Stage 031: its Inputs introduce the stable interval `0 <= alpha_0 < alpha_crit` (`paper/stages/stage_031.tex:9`), it boxes the refined threshold `alpha_{\rm crit}=AB/(B\kappa_0^2+A\kappa_1^2)` (`paper/stages/stage_031.tex:39`, `paper/stages/stage_031.tex:43`), and its Output explicitly includes the refined threshold (`paper/stages/stage_031.tex:65`). Stage 031's own audits verify the threshold and determinant/softening identities (`scripts/moving_throat_pde_stage031_selected_branch_reachability_sympy_audit.py:86`, `scripts/moving_throat_pde_stage031_selected_branch_reachability_sympy_audit.py:94`, `mathematica/moving_throat_pde_stage031_selected_branch_reachability_mathematica_audit.wl:57`, `mathematica/moving_throat_pde_stage031_selected_branch_reachability_mathematica_audit.wl:61`). Trim the Stage 029 `alpha_crit` assertion and scope bullet instead of expanding the Stage 029 card.

---

## Q3 — Stage 035 F1 (paper_misalignment / target_mismatch, medium)

**Item:** The boxed equation `eq:app-stage035-F-derivative` in `paper/stages/stage_035.tex:65-73` and the notes file `notes/stages/moving_throat_pde_stage035_dimensionless_normalization_locus.md:85-87` state the bracket polynomial of ∂F/∂ξ as:

```
81 delta^3 + 72 delta^2 + 206 delta^2 xi + 297 delta xi^2 + 138 xi^3
```

Both engines, however, verify the polynomial:

```
81 delta^3 + 72 delta^2 + 189 delta^2 xi + 297 delta xi^2 + 121 xi^3
```

(scripts at `scripts/moving_throat_pde_stage035_..._sympy_audit.py:65-70` and `mathematica/...wl:64-69`). Independent hand derivation (quotient rule on F = (9δ+11ξ)⁴ / [81(1-ξ)M²] with M = 9δ² + 18δξ + 11ξ²) confirms the scripts: `189` and `121` are correct. A numerical check at δ=0 gives `dF/dξ = 121 / [81(1-ξ)²]`, matching the script form (the paper form gives `138 / [81(1-ξ)²]` which is wrong).

**Audit finding location:** `redteam/reports/stage_035.md` § F1

**Context:** Positivity (the monotonicity / existence-uniqueness theorem) survives either way because both polynomials have strictly positive coefficients on δ,ξ > 0. But `eq:app-stage035-F-derivative` is a load-bearing boxed identity referenced from `\stagefield{Output}`, so the wrong coefficients undermine the exact-closure claim and could propagate to any downstream consumer that quotes the polynomial.

**Options:**
- **(a) Fix the paper and notes** — update `paper/stages/stage_035.tex:71` and `notes/stages/moving_throat_pde_stage035_...md:86` so the bracket polynomial reads `81 delta^3 + 72 delta^2 + 189 delta^2 xi + 297 delta xi^2 + 121 xi^3`. Scripts unchanged. This is the auditor's recommended direction.
- **(b) Reverse the scripts** — if there's any chance the paper coefficients are correct and the auditor's hand derivation is wrong. (Auditor cross-checked numerically; this is very unlikely.)
- **(c) Add a paper-side note** acknowledging the discrepancy and committing to recompute later. (Not advisable for a load-bearing identity.)
- **skip**

## Recommendation

direction: a
rationale: The paper and notes currently use the `206 delta^2 xi` and `138 xi^3` coefficients (`paper/stages/stage_035.tex:71`, `notes/stages/moving_throat_pde_stage035_dimensionless_normalization_locus.md:86`), while both audit engines use `189 delta^2 xi` and `121 xi^3` (`scripts/moving_throat_pde_stage035_dimensionless_normalization_locus_sympy_audit.py:65`, `scripts/moving_throat_pde_stage035_dimensionless_normalization_locus_sympy_audit.py:68`, `mathematica/moving_throat_pde_stage035_dimensionless_normalization_locus_mathematica_audit.wl:64`, `mathematica/moving_throat_pde_stage035_dimensionless_normalization_locus_mathematica_audit.wl:66`). Using the source formula `F=(9\delta+11\xi)^4/[81(1-\xi)(9\delta^2+18\delta\xi+11\xi^2)^2]` (`paper/stages/stage_035.tex:47`, `paper/stages/stage_035.tex:48`), set `S=9 delta+11 xi`, `M=9 delta^2+18 delta xi+11 xi^2`, and `M'=2S`; the quotient-rule/log-derivative form is `S^3[44(1-xi)M + SM - 4(1-xi)S^2]/[81(1-xi)^2 M^3]`. Expanding the bracket gives `44(1-xi)M = 396 delta^2 + 792 delta xi + 484 xi^2 - 396 delta^2 xi - 792 delta xi^2 - 484 xi^3`, `SM = 81 delta^3 + 261 delta^2 xi + 297 delta xi^2 + 121 xi^3`, and `-4(1-xi)S^2 = -324 delta^2 - 792 delta xi - 484 xi^2 + 324 delta^2 xi + 792 delta xi^2 + 484 xi^3`, summing to `81 delta^3 + 72 delta^2 + 189 delta^2 xi + 297 delta xi^2 + 121 xi^3`. The sanity check at `delta=0` gives `F=121/[81(1-xi)]` and `dF/dxi=121/[81(1-xi)^2]`, so the paper's `138` coefficient cannot be correct; this derivative is a paper Output item (`paper/stages/stage_035.tex:131`), so fix the paper and notes and leave scripts unchanged.

---

## End of questions

After Codex fills in `## Recommendation` blocks above, halt and present to user for review. Apply session is a separate Codex invocation with explicit per-question scope.

## Apply log

### Q1 applied
- direction: a
- files modified: `scripts/moving_throat_pde_stage029_dynamic_loading_sympy_audit.py:3,5`; `mathematica/moving_throat_pde_stage029_dynamic_loading_mathematica_audit.wl:33`
- destination_verified: n/a
- post-edit checks: `python3 scripts/moving_throat_pde_stage029_dynamic_loading_sympy_audit.py` exit 0; `math -script mathematica/moving_throat_pde_stage029_dynamic_loading_mathematica_audit.wl` exit 0
- notes: checks were run after the combined Q1/Q2 Stage 029 script edits.

### Q2 applied
- direction: b
- files modified: `scripts/moving_throat_pde_stage029_dynamic_loading_sympy_audit.py:17,189-194`; `mathematica/moving_throat_pde_stage029_dynamic_loading_mathematica_audit.wl:159-166`
- destination_verified: yes — `paper/stages/stage_031.tex:43`, `paper/stages/stage_031.tex:65`; `scripts/moving_throat_pde_stage031_selected_branch_reachability_sympy_audit.py:87`, `scripts/moving_throat_pde_stage031_selected_branch_reachability_sympy_audit.py:94`, `scripts/moving_throat_pde_stage031_selected_branch_reachability_sympy_audit.py:116`; `mathematica/moving_throat_pde_stage031_selected_branch_reachability_mathematica_audit.wl:59`, `mathematica/moving_throat_pde_stage031_selected_branch_reachability_mathematica_audit.wl:61`, `mathematica/moving_throat_pde_stage031_selected_branch_reachability_mathematica_audit.wl:71`
- post-edit checks: `python3 scripts/moving_throat_pde_stage029_dynamic_loading_sympy_audit.py` exit 0; `math -script mathematica/moving_throat_pde_stage029_dynamic_loading_mathematica_audit.wl` exit 0
- notes: the mandatory grep confirmed Stage 031 script ownership; the Stage 031 paper card uses TeX spelling `\alpha_{\rm crit}`, so its owner lines were verified by direct line readback.

### Q3 applied
- direction: a
- files modified: `paper/stages/stage_035.tex:71`; `notes/stages/moving_throat_pde_stage035_dimensionless_normalization_locus.md:86`
- destination_verified: n/a
- post-edit checks: n/a — paper/notes only; targeted readback confirmed the corrected bracket polynomial and intact surrounding displayed-equation delimiters.
- notes: Stage 035 scripts were not modified or rerun, per scope.
