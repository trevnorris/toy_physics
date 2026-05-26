You are acting as math authority and cross-stage researcher for the red-team v2 paper-grounded audit of batch III.1 (stages 037–048) of the moving-throat PDE derivation.

# Your task

Fill in the `## Recommendation` block beneath each `## Q<N>` heading in
`/var/projects/toy_physics/research/pde_ledger/redteam/resolutions/batch_III1_paper_alignment.md`.

For each question:
1. Read the question fully (item description, audit finding location, context, options).
2. Read the cited files (paper TeX, notes, scripts) to verify the auditor's claims.
3. Write your recommendation:
   ```
   ## Recommendation

   direction: <a|b|c|skip>
   rationale: <2–6 sentences citing file:line evidence>
   ```

# Per-question research guidance

**Q1 (043 F2 D_phi sign)**: Check whether `D_phi` is referenced in any downstream stage (044, 045, 046, 047, 048). Grep for `D_phi`, `dPhi`, `D_{\phi}`, `dphi` across `paper/stages/stage_04[3-8].tex`, `scripts/moving_throat_pde_stage04[4-8]_*.py`, `mathematica/moving_throat_pde_stage04[4-8]_*.wl`, and `notes/stages/moving_throat_pde_stage04[3-8]_*.md`. Find which convention downstream consumers use. Whichever side has more downstream callers should determine the canonical direction.

**Q2 (045 F3 Stage-27 residual)**: "Stage 27" in 045's notes is current paper numbering Stage 044 ("Continuum-Selected Rank-2 Closure"). Read `paper/stages/stage_044.tex` (especially the boxed normalization-residual equation), `notes/stages/moving_throat_pde_stage044_continuum_selected_rank2_closure.md`, and `scripts/moving_throat_pde_stage044_continuum_selected_rank2_closure_sympy_audit.py`. Identify the canonical `R_target` residual expression. Then check whether substituting `R_phi = R_U`, `lam_0 = 2/9` into it yields the `F_tr` form in 045's notes. If yes, direction (a) is mechanically applicable. If 044 doesn't have such a residual or it's in a different form, consider (c) move to downstream.

**Q3 (045 F4 script label drift)**: Note the I.2/II.1 precedent (stage 022, 028, 029 F1 were treated as cosmetic-only without user gate). Recommend on that basis. Also check `notes/stages/moving_throat_pde_stage045_coherent_local_tracking.md` for "Stage 28" mentions — if present, that's an extra fix target for option (b).

**Q4 (046 F1 notes coefficient typos)**: Independently verify by hand expansion that at least one of the disputed coefficients (e.g., `P_R[R delta^3] = 162` vs `230`) is correct. Use the script's `F_tr` definition at `scripts/moving_throat_pde_stage046_tracking_branch_bounds_sympy_audit.py:44-48`. Expansion sketch: `P_R = (1-xi) * (denominators)^2 * d/dR [stuff]` — work out the term-by-term contribution to `R delta^3`. Cite the Mathematica output line 14 and 20 evidence as cross-check.

# Critical rules

- **Recommendations ONLY.** Do not edit any script, paper TeX, or notes file in this session. Apply is a separate session.
- **No fake commentary scripts.** Do not write `python3 -c` commentary blocks. Read and reason directly.
- **Cite file:line for every claim** in your rationale.
- **For Q4 coefficient check** — write out the expansion step-by-step in your rationale so the user can audit your work.
- **For Q1 / Q2 downstream check** — explicitly cite the destination stage and its file:line if you find one.

# Working directory

`/var/projects/toy_physics/research/pde_ledger`

# When you are done

Print a one-line summary per question: `Q<N>: direction=<x>, key evidence: <file:line>`. Do not modify any other file.
