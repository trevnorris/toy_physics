# Codex recommendation pass — batch III.3 paper-alignment items

You are filling in `## Recommendation` blocks in `redteam/resolutions/batch_III3_paper_alignment.md`. Read that file in full first; it contains two questions (Q1 = stage 062 F1, Q2 = stage 062 F2) with full context, options, and per-question destination-verification research guidance.

For each question:

1. Read the relevant paper card (`paper/stages/stage_062.tex`) and notes (`notes/stages/moving_throat_pde_stage062_parent_action_gain.md`).
2. Read the audit report (`redteam/reports/stage_062.md`) for the full finding text.
3. For destination-verification: grep the codebase (paper/stages/stage_063.tex through stage_072.tex, all scripts/*.py and mathematica/*.wl in batch III.3, plus stage 064 specifically for Q2) for usage of the contested quantity:
   - Q1: search for `N_φφ · C_{σφ}²`, `C_sp`, `C_sigma_phi`, `O_sp²/(N_ss N_pp)`, Cauchy-Schwarz language
   - Q2: search for explicit usage of `σ_*(φ)` expression (not just `G_micro`); look at stage 064 ("equilibrium source profile") in detail.
4. Cite specific file:line evidence in the rationale.
5. Append a `## Recommendation` block under each question with the schema:

   ```markdown
   ## Recommendation

   - direction: <a | b | c>
   - rationale: <one paragraph citing file:line evidence from destination verification>
   - downstream_impact: <which downstream stages consume the contested quantity, if any>
   - notes: <optional caveats>
   ```

Do NOT edit any paper, notes, or scripts. Recommendations only — application is a separate pass.

Do NOT auto-apply. The user reviews recommendations and approves before any edits land.
