---
unit_id: 062
batch: III.3
created_at: 2026-05-26T17:40:40Z
findings_count: 3
stop_cold: null
applied: true
applied_at: 2026-05-26T18:40:40Z
verification_status: pending
needs_user_resolution: true
---

# Codex directive — unit 062

Apply each non-paper_misalignment finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

For `paper_misalignment` findings (F1, F2 below), do nothing — the orchestrator is holding for user resolution. Do not edit paper.tex, notes/, or scripts to "fix" a paper_misalignment unless the user has explicitly chosen a direction in a follow-up directive.

If a non-paper_misalignment finding's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named.

After editing, RUN the affected scripts (`python3 <path>` for SymPy, `math -script <path>` for Mathematica) and iterate until they exit 0 with all in-file checks passing. Getting the scripts to run cleanly is your job; the orchestrator independently re-runs afterward.

Do NOT touch paper.tex, notes/, or any prose documents. The red-team only modifies scripts (except when a follow-up directive explicitly authorizes a paper-side edit after user resolution).

## F1 — paper_misalignment

**Subtype:** `script_missing_paper_claim`

**Paper side:**
- `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_062.tex:33-46` quote (boxed eq. `app-stage062-Gmicro` plus eq. `app-stage062-Csigma`):
  > `G_micro = rho_* g_phi^2 O_{sigma phi}^2 / (m c_{s,*}^2 K_X N_{sigma sigma}) = (rho_* g_phi^2 N_{phi phi} / (m c_{s,*}^2 K_X)) * C_{sigma phi}^2,` with `C_{sigma phi}^2 := O_{sigma phi}^2 / (N_{sigma sigma} N_{phi phi}), 0 <= C_{sigma phi}^2 <= 1.`

**Script side:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage062_parent_action_gain_sympy_audit.py:88` quote:
  > `print("Coherence factor (definition):  C_(sigma phi)^2 := Osp^2 / (Nss Npp)")`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage062_parent_action_gain_mathematica_audit.wl:82` quote:
  > `Print["Coherence factor (definition):  C_(sigma phi)^2 := O_sp^2 / (N_ss N_pp)"];`

## Resolve before fix_loop

The paper's boxed `G_micro` has two equalities and a bound. The scripts only assert the first equality and print the coherence-factor *definition*; the second equality `G_micro = (rho_* g_phi^2 N_{phi phi}/(m c_{s,*}^2 K_X)) * C_{sigma phi}^2` and the Cauchy-Schwarz bound `0 <= C_{sigma phi}^2 <= 1` have no script-side check. Should the script be extended to cover both, or should the paper card be trimmed to only the first equality?

Possible directions (the user picks one):
- (a) Script should cover both — Codex will add to both engines: (1) a `C_sp_sq = Osp**2/(Nss*Npp)` definition followed by `expect_zero("Second equality of boxed G_micro", G_micro_closed - (rho_star*g_phi**2*Npp/(m*cs_star_sq*KX))*C_sp_sq)`, and (2) substitute `Osp -> cos(theta) * sqrt(Nss*Npp)` for a symbolic `theta` and verify `C_sp_sq` reduces to `cos(theta)**2` (which is manifestly in `[0,1]`). Re-run sympy+mathematica.
- (b) Paper card is overspecified — trim eq. `app-stage062-Gmicro` to its first equality and remove eq. `app-stage062-Csigma`. Notes §5 can keep the factorization as commentary but not as a verified deliverable. No script change.
- (c) Keep paper as-is and document in the notes that the second equality and bound are reader-level algebra, not auditor-verified. No script or paper change beyond a notes annotation.

The orchestrator will not invoke Codex on this unit (for F1) until the user has chosen a direction.

## F2 — paper_misalignment

**Subtype:** `notes_contradicts_script`

**Paper side:**
- `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage062_parent_action_gain.md:44` quote:
  > `F_red[sigma,phi] = int_0^L ds [ (Theta_sigma/2) sigma^2 - Lambda_phi sigma phi + (T_X/2) phi_s^2 + (K_X/2) phi^2 ] + (K_m/2) phi(0)^2,`
- (Origin of the sign: same file, line 114-122: `delta V_conf(s,y) = - g_phi chi_phi(y) phi(s)` and line 135: `- int_0^L ds sigma(s) phi(s) g_phi int d^3y chi_sigma chi_phi`.)

**Script side:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage062_parent_action_gain_sympy_audit.py:71-75` quote:
  > `S_parent = (sp.Rational(1, 2) * Theta_sigma * sigma**2 + Lambda_phi * sigma * phi + sp.Rational(1, 2) * KX * phi**2)`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage062_parent_action_gain_mathematica_audit.wl:65` quote:
  > `sParent = (1/2)*thetaSigma*sigma^2 + lambdaPhi*sigma*phi + (1/2)*kX*phi^2;`

## Resolve before fix_loop

The notes derive `F_red` with `-Lambda_phi sigma phi` (minus on the cross-coupling), tracing back to the physical perturbation `delta V_conf = -g_phi chi_phi phi`. Both scripts use `+Lambda_phi sigma phi`. The gain magnitude `Lambda^2/(Theta K)` is invariant under this sign flip and all assertions still pass, but the on-shell `sigma_star(phi)` prints with the opposite sign from what the notes' `F_red` would give. Which sign is canonical?

Possible directions (the user picks one):
- (a) Notes' sign is canonical (matches the physical derivation `delta V_conf = -g_phi chi_phi phi`). Codex will change `S_parent` / `sParent` to `(1/2)*Theta_sigma*sigma**2 - Lambda_phi*sigma*phi + (1/2)*KX*phi**2` in both engines; the gain assertion will still pass, the printed `sigma_star` will flip sign, and downstream consumers of `sigma_star(phi)` will now match notes.
- (b) Script's sign is the working convention used downstream (perhaps because of how `Lambda_phi` is defined to be positive). Codex will not edit scripts; user must edit notes line 44 (and any related text) to drop the minus and add a remark explaining the convention.
- (c) The two sides should be kept consistent by changing both the notes derivation chain (replace `delta V_conf = +g_phi chi_phi phi` so the minus sign disappears at its source) and the scripts to match — a deeper notes/paper review the user must authorize.

The orchestrator will not invoke Codex on this unit (for F2) until the user has chosen a direction.

## F3 — mathematica_transliteration

**Target:** `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage062_parent_action_gain_mathematica_audit.wl:63-80`

**Issue:**

The Mathematica script is a line-by-line port of the SymPy script's algebraic chain: same `S_parent` recipe, same `Solve` for `sigma_star`, same `Coefficient[..., phi, 2]` extraction, same closed-form target. The only extra check (`gainFromSeries` at line 71, asserted at line 79) compares two coefficient-extraction methods on the *same* `sEff` object, which is internal-consistency, not an independent derivation. Both engines walk identical algebra.

**Required change:**

Add a second, structurally distinct derivation of `G_micro` in the Mathematica file using the susceptibility chain `chi_sigma^(eff) = 1/Theta_sigma`, `G_micro = chi_sigma^(eff) * Lambda_phi^2 / K_X` (notes §4). This route bypasses `Solve` and `Coefficient` entirely.

Edit the Mathematica file between lines 63 and 80:

**Before** (current lines 63-80):

```
thetaSigma = (m*csStarSq/rhoStar)*nSS;
lambdaPhi = gPhi*oSP;
sParent = (1/2)*thetaSigma*sigma^2 + lambdaPhi*sigma*phi + (1/2)*kX*phi^2;
sigmaStar = sigma /. First[Solve[D[sParent, sigma] == 0, sigma]];
sigmaStar = sigmaStar /. ConditionalExpression[e_, _] :> e;
sEff = Expand[sParent /. sigma -> sigmaStar];
phiQuadraticCoeff = FullSimplify[Coefficient[sEff, phi, 2], Assumptions -> $Assumptions];
gainFromAction = FullSimplify[(kX - 2*phiQuadraticCoeff)/kX, Assumptions -> $Assumptions];
gainFromSeries = FullSimplify[(kX - 2*SeriesCoefficient[Series[sEff, {phi, 0, 2}], 2])/kX, Assumptions -> $Assumptions];
gClosed = FullSimplify[rhoStar*gPhi^2*oSP^2/(m*csStarSq*kX*nSS), Assumptions -> $Assumptions];

Print["Theta_sigma = ", fmt[thetaSigma]];
Print["Lambda_phi = ", fmt[lambdaPhi]];
Print["sigmaStar = ", fmt[sigmaStar]];
Print["S_eff = ", fmt[sEff]];
Print["G_micro from action = ", fmt[gainFromAction]];
expectZero["Mathematica two-route consistency", gainFromAction - gainFromSeries];
expectZero["gMicro from parent action vs closed form", gainFromAction - gClosed];
```

**After** (insert susceptibility-route block; keep the action-coefficient route as cross-check):

```
thetaSigma = (m*csStarSq/rhoStar)*nSS;
lambdaPhi = gPhi*oSP;

(* Susceptibility route (notes section 4):
   chi_sigma^(eff) = 1/Theta_sigma; G_micro = chi_sigma^(eff) * Lambda_phi^2 / K_X.
   Independent algebraic path from the action-coefficient route below. *)
chiSigmaEff = FullSimplify[1/thetaSigma, Assumptions -> $Assumptions];
gainViaSusceptibility = FullSimplify[chiSigmaEff*lambdaPhi^2/kX, Assumptions -> $Assumptions];

(* Action-coefficient route (kept as redundant second route). *)
sParent = (1/2)*thetaSigma*sigma^2 + lambdaPhi*sigma*phi + (1/2)*kX*phi^2;
sigmaStar = sigma /. First[Solve[D[sParent, sigma] == 0, sigma]];
sigmaStar = sigmaStar /. ConditionalExpression[e_, _] :> e;
sEff = Expand[sParent /. sigma -> sigmaStar];
phiQuadraticCoeff = FullSimplify[Coefficient[sEff, phi, 2], Assumptions -> $Assumptions];
gainFromAction = FullSimplify[(kX - 2*phiQuadraticCoeff)/kX, Assumptions -> $Assumptions];
gainFromSeries = FullSimplify[(kX - 2*SeriesCoefficient[Series[sEff, {phi, 0, 2}], 2])/kX, Assumptions -> $Assumptions];
gClosed = FullSimplify[rhoStar*gPhi^2*oSP^2/(m*csStarSq*kX*nSS), Assumptions -> $Assumptions];

Print["Theta_sigma = ", fmt[thetaSigma]];
Print["Lambda_phi = ", fmt[lambdaPhi]];
Print["chi_sigma^(eff) = ", fmt[chiSigmaEff]];
Print["G_micro via susceptibility route = ", fmt[gainViaSusceptibility]];
Print["sigmaStar = ", fmt[sigmaStar]];
Print["S_eff = ", fmt[sEff]];
Print["G_micro from action = ", fmt[gainFromAction]];
expectZero["G_micro via susceptibility route vs closed form", gainViaSusceptibility - gClosed];
expectZero["G_micro: action route equals susceptibility route", gainFromAction - gainViaSusceptibility];
expectZero["Mathematica two-route consistency", gainFromAction - gainFromSeries];
expectZero["gMicro from parent action vs closed form", gainFromAction - gClosed];
```

The susceptibility-route block (`chiSigmaEff`, `gainViaSusceptibility`) uses only the polytrope-derived coefficients `thetaSigma`, `lambdaPhi` and arithmetic — no `Solve`, no `Coefficient`, no `Series`. It is a genuinely independent algebraic path from the action-coefficient route. The two new `expectZero` calls (`... vs closed form` and `... action route equals susceptibility route`) verify both that the susceptibility route arrives at `gClosed` independently and that the two engines' internal routes agree.

Do not touch the SymPy file. Do not touch the EOS block (lines 28-55) or the `Xi_micro` block (lines 84-93).

**Verification command:**

After Codex applies, the verifier will run `redteam exec-mathematica 062` and confirm:
1. The script exits 0.
2. The new transcript contains the lines `PASS: G_micro via susceptibility route vs closed form` and `PASS: G_micro: action route equals susceptibility route` in addition to the existing PASS lines.
3. The script source contains the new symbol `chiSigmaEff` and the variable `gainViaSusceptibility`.

If the Mathematica engine produces a non-zero residual for `gainViaSusceptibility - gClosed`, that means `chi_sigma^(eff)` or `lambda_phi` is being constructed differently — escalate to user; do not auto-paper over.

## Apply log

- applied_at: 2026-05-26T18:40:40Z
- direction: F3 (susceptibility route)
- files_changed:
  - mathematica/moving_throat_pde_stage062_parent_action_gain_mathematica_audit.wl
- sympy_exit: 0
- mathematica_exit: 0
- summary: Added the independent Mathematica susceptibility-route derivation and assertions against both the closed form and the action route.
