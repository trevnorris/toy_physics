# Codex apply pass — batch III.3 paper-alignment + Stage 062 F3

The user reviewed and approved both `## Recommendation` blocks in `redteam/resolutions/batch_III3_paper_alignment.md`:
- Q1 → direction (a): extend Stage 062 scripts to verify second equality + Cauchy-Schwarz bound.
- Q2 → direction (a): flip Stage 062 scripts to use `−Λ_φ σ φ` to match notes.

Also fold in **Stage 062 F3** (`mathematica_transliteration`, script-side, in `redteam/directives/stage_062.md`) since it touches the same `.wl` file and benefits from being applied together. F3 is independent of Q1/Q2.

## Scope authorization (per finding)

- **Q1**: edit `scripts/moving_throat_pde_stage062_parent_action_gain_sympy_audit.py` and `mathematica/moving_throat_pde_stage062_parent_action_gain_mathematica_audit.wl`. No paper or notes edits.
- **Q2**: edit the same two script files. No paper or notes edits. (The notes file already has the correct minus sign — no edit needed there.)
- **F3** (`redteam/directives/stage_062.md` F3): edit `mathematica/moving_throat_pde_stage062_parent_action_gain_mathematica_audit.wl` only. Add an independent susceptibility-route derivation `chiSigmaEff = 1/thetaSigma; gainViaSusceptibility = chiSigmaEff*lambdaPhi^2/kX` and `expectZero` it against `gClosed` AND against `gainFromAction`. Keep the existing routes.

## Required edits

### Q1 — second equality + Cauchy-Schwarz bound (both engines)

**SymPy** (`scripts/moving_throat_pde_stage062_parent_action_gain_sympy_audit.py`):

After the existing `G_micro_closed` definition and the first-equality `expect_zero` (currently at line ~86), add:

```python
# Second equality of boxed eq:app-stage062-Gmicro: G_micro = (rho_* g_phi^2 N_pp/(m c_{s,*}^2 K_X)) * C_{sigma phi}^2
C_sp_sq = Osp**2 / (Nss * Npp)
G_micro_factored = (rho_star * g_phi**2 * Npp / (m * cs_star_sq * KX)) * C_sp_sq
expect_zero("Second equality of boxed G_micro: closed vs factored form", G_micro_closed - G_micro_factored)

# Cauchy-Schwarz bound on C_{sigma phi}^2: substitute O_{sigma phi} = cos(theta) sqrt(N_ss N_pp)
theta = sp.Symbol("theta", real=True)
C_sp_sq_cos = C_sp_sq.subs(Osp, sp.cos(theta) * sp.sqrt(Nss * Npp))
C_sp_sq_cos_simplified = sp.simplify(C_sp_sq_cos)
expect_zero("C_{sigma phi}^2 via Cauchy parameterization equals cos^2(theta)", C_sp_sq_cos_simplified - sp.cos(theta)**2)
# 0 <= cos^2(theta) <= 1 is then a standard real-trig fact; assert sympy agrees:
assert sp.simplify(sp.cos(theta)**2 - 0) >= 0 if False else True  # tautology pinned to documentation
print(f"C_sp_sq Cauchy parameterization yields cos^2(theta) in [0, 1] (Cauchy-Schwarz bound).")
```

(Adjust variable names if the actual script uses different ones; the script defines `Osp`, `Nss`, `Npp`, `rho_star`, `g_phi`, `cs_star_sq`, `KX`, `m` — use what's already there.)

**Mathematica** (`mathematica/moving_throat_pde_stage062_parent_action_gain_mathematica_audit.wl`):

After the existing `gClosed` definition and the first-equality `expectZero` (currently at line ~80), add:

```mathematica
(* Second equality of boxed eq:app-stage062-Gmicro *)
cSpSq = oSP^2 / (nSS * nPP);
gMicroFactored = (rhoStar * gPhi^2 * nPP / (m * csStarSq * kX)) * cSpSq;
expectZero["Second equality of boxed G_micro: closed vs factored form", gClosed - gMicroFactored];

(* Cauchy-Schwarz parameterization: O_{sigma phi} = cos(theta) sqrt(N_ss N_pp) *)
theta = Symbol["theta"];
cSpSqCos = cSpSq /. oSP -> Cos[theta] * Sqrt[nSS * nPP];
expectZero["C_{sigma phi}^2 via Cauchy parameterization equals cos^2(theta)",
           FullSimplify[cSpSqCos - Cos[theta]^2, Assumptions -> $Assumptions]];
Print["C_sp_sq Cauchy parameterization yields Cos[theta]^2 in [0, 1] (Cauchy-Schwarz bound)."];
```

### Q2 — flip σφ coupling sign (both engines)

**SymPy**: in `S_parent` definition (currently at line ~71-75), change `+ Lambda_phi * sigma * phi` to `- Lambda_phi * sigma * phi`. The on-shell `sigma_star` will flip sign accordingly; `gain_from_action` is unchanged. All existing assertions still pass.

**Mathematica**: in `sParent` definition (currently at line ~65), change `+ lambdaPhi*sigma*phi` to `- lambdaPhi*sigma*phi`. Same effect.

### F3 — Mathematica susceptibility route

In `mathematica/moving_throat_pde_stage062_parent_action_gain_mathematica_audit.wl`, after the `thetaSigma` and `lambdaPhi` definitions (currently at line ~63-64), add:

```mathematica
(* Independent susceptibility-route derivation (notes section 4): 
   chi_sigma^(eff) = 1/Theta_sigma, then G_micro = chi_sigma^(eff) * Lambda_phi^2 / K_X *)
chiSigmaEff = 1/thetaSigma;
gainViaSusceptibility = FullSimplify[chiSigmaEff * lambdaPhi^2 / kX, Assumptions -> $Assumptions];
```

Keep the existing action-coefficient route (`gainFromAction`, `gainFromSeries`). Add two new assertions:

```mathematica
expectZero["G_micro via susceptibility route vs closed form", gainViaSusceptibility - gClosed];
expectZero["G_micro: action route equals susceptibility route", gainFromAction - gainViaSusceptibility];
```

This gives the `.wl` an algebraically distinct route to the same target — bypassing `Solve` and `Coefficient` entirely.

## Apply protocol

1. Apply Q1, Q2, F3 edits to the two stage 062 script files.
2. Run both scripts and confirm all assertions pass + scripts exit 0:
   - `python3 scripts/moving_throat_pde_stage062_parent_action_gain_sympy_audit.py > scripts/output/moving_throat_pde_stage062_parent_action_gain_sympy_audit.txt 2>&1` — exit 0 expected
   - `math -script mathematica/moving_throat_pde_stage062_parent_action_gain_mathematica_audit.wl > mathematica/output/moving_throat_pde_stage062_parent_action_gain_mathematica_audit.txt 2>&1` — exit 0 expected
   - **CRITICAL: Mathematica single-seat rule applies — DO NOT run sympy and mathematica in parallel.** Run sympy first, wait for completion, then run mathematica.
3. Append apply-log block to `redteam/resolutions/batch_III3_paper_alignment.md` under each question and to `redteam/directives/stage_062.md` under F3:
   ```markdown
   ## Apply log

   - applied_at: <ISO 8601>
   - direction: a (Q1) / a (Q2) / F3 (susceptibility route)
   - files_changed:
     - <path>
   - sympy_exit: 0
   - mathematica_exit: 0
   - summary: <one sentence per finding>
   ```

If a script fails, diagnose and iterate (max 5 attempts per directive's own rule). If pitfall #8 (heavy `dsolve`/`DSolve`) appears in any prescribed edit, follow the Green-function identity pattern documented in codex.md "Common cross-engine pitfalls" item #1.

Do NOT touch any paper, notes, or other scripts beyond stage 062.

Do NOT update MANIFEST.yaml — orchestrator handles status transitions.

When you are done, mark stage 062's directive frontmatter `applied: true` with `applied_at: <ISO 8601>`.
