You are acting as math authority and cross-stage researcher for the red-team v2 paper-grounded audit of batch III.2 (stages 049–060, "Tracking, zeta thresholds, asymmetry, boost") of the moving-throat PDE derivation.

# Your task

Fill in the `## Recommendation` block beneath each `## Q<N>` heading in
`/var/projects/toy_physics/research/pde_ledger/redteam/resolutions/batch_III2_paper_alignment.md`.

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

**Q1 (050 F2 enhancement-ceiling theorem `S_n^(max)`)**: The scripts assert `S_n^(twin) < S_n^(max) := 1 + (1-eps)/((2n+1)^2 - eps)` (`scripts/moving_throat_pde_stage050_zeta_threshold_comparison_sympy_audit.py:90,94,104-107`; `mathematica/moving_throat_pde_stage050_zeta_threshold_comparison_mathematica_audit.wl:83`). Notes describe `S_n^(max)` in Section 5 (`notes/stages/moving_throat_pde_stage050_zeta_threshold_comparison.md:172-174`). Paper card body has no boxed equation for `S_n^(max)` (`paper/stages/stage_050.tex:11-42`) and Output line doesn't mention a ceiling (`paper/stages/stage_050.tex:44`).

Research questions to answer:
- Does any downstream stage (051, 052, 053, 054, 055) consume `S_n^(max)` as an upstream input? Grep `S_n.*max`, `enhancement.*ceiling`, `S_n_max`, `sNMax`, `S_max` across `paper/stages/stage_05[1-9].tex`, `notes/stages/moving_throat_pde_stage05[1-9]_*.md`, and `scripts/moving_throat_pde_stage05[1-9]_*_sympy_audit.py`.
- If downstream stages use it, the ceiling needs to live somewhere — option (a) advertising in stage 050's paper card is the natural home (since notes Section 5 already covers it).
- If no downstream stage uses `S_n^(max)`, the ceiling is locally derivable and (b) "ceiling belongs to a different stage" may not have a natural target — meaning (a) is still the answer or it should be demoted to an unboxed remark in the paper.
- Stage 050's paper card already mentions "higher-harmonic exclusion/softness thresholds" in the Output — is the ceiling a softness threshold? If so, advertising it explicitly is option (a).

**Q2 (057 F1 `partial_Pe zeta > 0` carry-forward)**: Paper card claims monotonicity in Pe (`paper/stages/stage_057.tex:23`), notes anchor it to Stage 056 (`notes/stages/moving_throat_pde_stage057_physical_parameter_map.md:140-148`). The scripts (sympy:62-73, wl:56-67) only compute `partial_kappa` and `partial_y`; no Pe derivative.

Research questions to answer:
- Does Stage 056 actually prove `dOmega_Pe/dPe > 0`? Read `paper/stages/stage_056.tex`, `notes/stages/moving_throat_pde_stage056_*.md`, and `scripts/moving_throat_pde_stage056_*_sympy_audit.py`. Look for an explicit sign check on `D[Omega_Pe, Pe]` (or `dOmega_Pe/dPe`, or `partial_Pe Omega_Pe`).
- If Stage 056 has a clean sign check, the carry-forward is sound and option (b) is conservative (paper/notes edit, no script change).
- If Stage 056 lacks it, the chain is broken — option (a) (re-verify locally in Stage 057) is mandatory.
- Stage 057's `zeta_phys = Omega_Pe^2 * (kappa+pi^2/4)/(kappa+y^2)` — note that `Omega_Pe^2` factors out, so `partial_Pe zeta_phys = 2 Omega_Pe (dOmega_Pe/dPe) * (kappa+pi^2/4)/(kappa+y^2)`. Positivity reduces to `Omega_Pe > 0` AND `dOmega_Pe/dPe > 0` on the constructive branch (`Pe >= 0`). The first is from a separate stage; the second is the actual Pe-monotonicity claim.
- Comment on the destination-verification guardrail: if you recommend (b), explicitly cite the Stage 056 line that proves `dOmega_Pe/dPe > 0`.

# Critical rules

- **Recommendations ONLY.** Do not edit any script, paper TeX, or notes file in this session. Apply is a separate session.
- **No fake commentary scripts.** Do not write `python3 -c` commentary blocks. Read and reason directly.
- **Cite file:line for every claim** in your rationale.
- **For Q1 downstream grep** — explicitly state which downstream stages (if any) consume `S_n^(max)`. Cite the file:line if found.
- **For Q2 destination check** — read the Stage 056 source files and cite the specific assertion (or its absence). If absent, lean toward (a). If present, cite the file:line evidence as the foundation for (b).

# Working directory

`/var/projects/toy_physics/research/pde_ledger`

# When you are done

Print a one-line summary per question: `Q<N>: direction=<x>, key evidence: <file:line>`. Do not modify any other file.
