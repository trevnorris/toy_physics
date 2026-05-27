---
unit_id: 062
batch: III.3
auditor_model: claude-opus-4-7-1m
audit_date: 2026-05-26T17:40:40Z
verdict: findings
stop_cold: null
findings_count: 3
paper_alignment: partial
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files:
    - /var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage062_parent_action_gain.md
  paper_appendix: present
---

# Audit unit 062 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_062.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage062_parent_action_gain.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part03.tex` (this stage appears as `\input{stages/stage_062}` on line 242; no inline summary row beyond that)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage062_parent_action_gain_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage062_parent_action_gain_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage062_parent_action_gain_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage062_parent_action_gain_mathematica_audit.txt`

## What the paper claims

Stage 062 projects the microscopic support/source gain into parent-action variables for the frozen `n=5` GNLS EOS. The card's `\stagefield{Output}` is the boxed equation eq. `app-stage062-Gmicro`:

> `G_micro = rho_* g_phi^2 O_{sigma phi}^2 / (m c_{s,*}^2 K_X N_{sigma sigma}) = [rho_* g_phi^2 N_{phi phi} / (m c_{s,*}^2 K_X)] * C_{sigma phi}^2`

with `C_{sigma phi}^2 := O_{sigma phi}^2 / (N_{sigma sigma} N_{phi phi})` and `0 <= C_{sigma phi}^2 <= 1` (eq. `app-stage062-Csigma`). The card also fixes the EOS identities `h'(rho_*) = 5 K rho_*^3 = m c_{s,*}^2/rho_*` (eq. `app-stage062-hprime`) and the overlap-invariant definitions (eq. `app-stage062-overlaps`). The notes file enumerates further deliverables: the projected coefficients `Theta_sigma = h'(rho_*) N_{sigma sigma}` and `Lambda_phi = g_phi O_{sigma phi}`, the reduced action `F_red[sigma,phi] = int_0^L ds [(Theta_sigma/2)sigma^2 - Lambda_phi sigma phi + (T_X/2)phi_s^2 + (K_X/2)phi^2] + (K_m/2)phi(0)^2` (note the **minus** on the sigma-phi coupling), the effective susceptibility `chi_sigma^(eff) = rho_*/(m c_{s,*}^2 N_{sigma sigma})`, and `Xi_micro = kappa G_micro` with `kappa = K_X L^2/T_X`.

## What the script claims to verify

The SymPy and Mathematica scripts both run the same audit list: (1) the general polytrope identity `h'(rho) = m c_s^2/rho` with an explicit `n=5` specialization and a non-tautology probe using `n+1` in `c_s^2`; (2) the projected coefficients `Theta_sigma = (m c_{s,*}^2/rho_*) N_{sigma sigma}` and `Lambda_phi = g_phi O_{sigma phi}`; (3) construction of a local quadratic action `S_parent = (1/2) Theta_sigma sigma^2 + Lambda_phi sigma phi + (1/2) K_X phi^2`, on-shell elimination of `sigma`, and extraction of the quadratic-in-phi coefficient yielding `(K_X - 2*coeff)/K_X`, asserted equal to the first-form closed expression `rho_* g_phi^2 O_{sigma phi}^2/(m c_{s,*}^2 K_X N_{sigma sigma})`; (4) solving `Xi_micro - Xi_target = 0` for `kappa` and asserting the solution equals `K_X L^2/T_X`. The coherence factor `C_{sigma phi}^2` is *printed as a definition* but never used in any assertion.

## Paper ↔ script cross-check

| Paper deliverable | Script-side check | Status |
|---|---|---|
| `h'(rho_*) = m c_{s,*}^2/rho_*` (eq. app-stage062-hprime) | `expect_zero("h'(rho) = m c_s^2 / rho (general polytrope)", ...)` and n=5 specialization (sympy lines 42-56) / wl lines 39-47 | match |
| Overlap-invariant definitions `N_{ss}, N_{pp}, O_{sp}` (eq. app-stage062-overlaps) | Declared as positive real symbols only (sympy line 65; wl lines 57-61); not exercised analytically since they are abstract integrals | match (definitions only; nothing to verify beyond carrying them through) |
| First equality of boxed eq.: `G_micro = rho_* g_phi^2 O_sp^2/(m c_{s,*}^2 K_X N_ss)` | `expect_zero("G_micro from parent action vs closed form", gain_from_action - G_micro_closed)` (sympy line 86; wl lines 70-80) | match |
| Second equality of boxed eq.: `G_micro = [rho_* g_phi^2 N_pp/(m c_{s,*}^2 K_X)] * C_sp^2` | None - `print("Coherence factor (definition): ...")` only (sympy line 88; wl line 82) | missing |
| Cauchy-Schwarz bound `0 <= C_sp^2 <= 1` (paper line 45, notes §5) | None - `Osp, Nss, Npp` declared positive but with no Cauchy-Schwarz tie | missing |
| `Xi_micro = kappa G_micro` with `kappa = K_X L^2/T_X` (notes §6) | `kappa_solved == [KX*L**2/TX]` (sympy line 95; wl lines 87-93) | match |
| Notes §2 reduced action `F_red` with **-Lambda_phi sigma phi** coupling | Script builds `S_parent` with **+Lambda_phi sigma phi** (sympy line 73; wl line 65) | mismatch (sign convention) |
| Notes §3 effective susceptibility `chi_sigma^(eff) = rho_*/(m c_{s,*}^2 N_ss)` | Not exposed by a dedicated assertion; only implicit in the gain identity | partial (gain identity contains it, but no standalone check) |

Dominant pattern is partial alignment: the headline gain identity (first equality) is verified, but the second equality of the boxed paper result and the Cauchy-Schwarz bound are not exercised, and the sigma-phi coupling sign in the script's reduced action disagrees with the notes.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 42-45 | `expect_zero("h'(rho) = m c_s^2 / rho (general polytrope)", hprime_general - m*cs_sq_general/rho)` | eq. app-stage062-hprime (general polytrope route) | yes (general n, then specialized) |
| A2 | sympy | 53-56 | `expect_zero("n=5 specialization ...", ...)` | eq. app-stage062-hprime at n=5 | yes |
| A3 | sympy | 60 | `assert sp.simplify(residual_wrong) != 0` | non-tautology probe: differentiating the wrong exponent yields a non-zero residual at n=5 | yes (the residual `K rho^3 (5 - 6 rho)` is a non-constant polynomial in `rho`, so `simplify` cannot collapse it to 0; this is a real check that A1/A2 weren't trivially passing) |
| A4 | sympy | 86 | `expect_zero("G_micro from parent action vs closed form", gain_from_action - G_micro_closed)` | first equality of boxed eq. app-stage062-Gmicro | yes |
| A5 | sympy | 95 | `assert kappa_solved == [KX*L**2/TX]` | `kappa = K_X L^2/T_X` from notes §6 | yes |
| A6 | mathematica | 39 | `expectZero["h'(rho) = m c_s^2 / rho (general polytrope)", ...]` | mirrors A1 | yes |
| A7 | mathematica | 47 | `expectZero["n=5 specialization ...", ...]` | mirrors A2 | yes |
| A8 | mathematica | 52-55 | `If[FullSimplify[...] === 0, fail[...], Print[...]]` | mirrors A3 | yes |
| A9 | mathematica | 79 | `expectZero["Mathematica two-route consistency", gainFromAction - gainFromSeries]` | cross-check of `Coefficient[..., phi, 2]` vs `SeriesCoefficient[Series[..., {phi,0,2}], 2]` for same `sEff` | partial (this is internal-consistency: both routes operate on the same `sEff`, so they must agree by construction of polynomial coefficient extraction; the check verifies the two extraction methods agree, not that the gain is independently right) |
| A10 | mathematica | 80 | `expectZero["gMicro from parent action vs closed form", gainFromAction - gClosed]` | mirrors A4 | yes |
| A11 | mathematica | 90-93 | `If[FullSimplify[kappaSolved == kX*ell^2/tX] === True, pass[...], fail[...]]` | mirrors A5 | yes |

No script-side check anchors to the second equality of eq. app-stage062-Gmicro, to the Cauchy-Schwarz bound, or to the sign convention of the sigma-phi coupling in F_red.

## Findings

### F1 — paper_misalignment

**Subtype:** `script_missing_paper_claim`

**Severity:** medium

**Files:**
- `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_062.tex:33-46` (boxed eq. `app-stage062-Gmicro`, two equalities; `app-stage062-Csigma` with the bound `0 <= C_{sigma phi}^2 <= 1`)
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage062_parent_action_gain_sympy_audit.py:88` (`print("Coherence factor (definition): ...")` only - no assertion)
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage062_parent_action_gain_mathematica_audit.wl:82` (same pattern - print only)

**What's wrong:**

The paper card's `\stagefield{Output}` is the boxed equation eq. `app-stage062-Gmicro`, which contains **two** equalities:

(i) `G_micro = rho_* g_phi^2 O_{sp}^2/(m c_{s,*}^2 K_X N_{ss})` - checked by A4 / A10.

(ii) `G_micro = [rho_* g_phi^2 N_{pp}/(m c_{s,*}^2 K_X)] * C_{sp}^2`, with `C_{sp}^2 := O_{sp}^2/(N_{ss} N_{pp})` - **not** asserted anywhere. Both scripts only `Print["Coherence factor (definition): C_(sigma phi)^2 := O_sp^2 / (N_ss N_pp)"]` and stop.

The card also states `0 <= C_{sigma phi}^2 <= 1` (paper line 45-46; notes §5) - also not exercised by any assertion.

Both omissions are paper deliverables that have no corresponding script-side check.

**Why this matters:**

The paper's boxed equation has two faces. The script verifies one face and prints the definition of the second; it never actually shows the algebraic equivalence of the two forms or carries the Cauchy-Schwarz bound through. A reader auditing the paper from the script transcript cannot confirm the boxed equation's second equality nor the bound on `C_{sigma phi}^2` is mathematically realized. This is the kind of carried-forward identity that downstream stages (063 amplitude/coherence thresholds, 064 equilibrium source profile) will lean on - if a downstream stage cites the `N_{pp} * C_{sp}^2` form, it currently rests on a printed definition, not a verified identity.

**Resolution direction (user decides):**

The fix is small if the user wants script-side coverage: introduce `C_sp_sq = Osp**2/(Nss*Npp)` and a second closed form `G_micro_factored = (rho_star*g_phi**2*Npp/(m*cs_star_sq*KX))*C_sp_sq`, then `expect_zero("Second equality of boxed G_micro", G_micro_closed - G_micro_factored)`. For the bound, the symbolic check is to substitute `O_{sigma phi} = cos(theta) * sqrt(N_ss * N_pp)` and confirm `C_sp_sq = cos(theta)^2 in [0,1]`. Alternatively, the paper card could be trimmed to only the first equality. Codex must not auto-pick a direction.

### F2 — paper_misalignment

**Subtype:** `notes_contradicts_script`

**Severity:** low

**Files:**
- Notes `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage062_parent_action_gain.md:44` quote:
  > `F_red[sigma,phi] = int_0^L ds [ (Theta_sigma/2) sigma^2 - Lambda_phi sigma phi + (T_X/2) phi_s^2 + (K_X/2) phi^2 ] + (K_m/2) phi(0)^2,`
- Script `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage062_parent_action_gain_sympy_audit.py:71-75` quote:
  > `S_parent = (sp.Rational(1, 2) * Theta_sigma * sigma**2 + Lambda_phi * sigma * phi + sp.Rational(1, 2) * KX * phi**2)`
- Script `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage062_parent_action_gain_mathematica_audit.wl:65` quote:
  > `sParent = (1/2)*thetaSigma*sigma^2 + lambdaPhi*sigma*phi + (1/2)*kX*phi^2;`

**What's wrong:**

The notes derive `F_red` with `-Lambda_phi sigma phi` (minus sign), tracing back to the physical perturbation `delta V_conf = -g_phi chi_phi(y) phi(s)` so that the source-support coupling enters the energy with a minus sign. The scripts build `S_parent` with **+Lambda_phi sigma phi**. The gain `G_micro = Lambda^2/(Theta K)` is invariant under this sign flip (quadratic in Lambda), which is why all assertions still pass, but the on-shell solution prints with the opposite sign:

```
sigma_star = -O_sp*g_phi*phi*rho_star/(N_ss*cs_star_sq*m)
```

versus what the notes' `F_red` would give, namely

```
sigma_star = +O_sp*g_phi*phi*rho_star/(N_ss*cs_star_sq*m).
```

**Why this matters:**

For the gain identity alone, this is cosmetic. But the script docstring claims to verify the projected reduced action against the notes, and the action it builds is sign-inverted on the cross-coupling. If a downstream unit (e.g., Stage 064 equilibrium source profile) imports `sigma_star(phi)` rather than just the gain, it will pick up the wrong sign.

**Resolution direction (user decides):**

(a) The notes/paper sign is the physically correct one (matches `delta V_conf = -g_phi chi_phi phi` derivation) - change scripts to `... - Lambda_phi * sigma * phi + ...` on `S_parent`/`sParent`. The gain assertion will still pass; the printed `sigma_star` will flip sign to match the notes. **OR**
(b) The script's sign is the working convention used downstream - update the notes to drop the minus on `-Lambda_phi sigma phi` (and explain the convention). This requires a paper-side edit which the user must authorize.

### F3 — mathematica_transliteration

**Severity:** medium

**Files:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage062_parent_action_gain_mathematica_audit.wl:34-47` vs `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage062_parent_action_gain_sympy_audit.py:36-56`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage062_parent_action_gain_mathematica_audit.wl:63-80` vs `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage062_parent_action_gain_sympy_audit.py:64-86`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage062_parent_action_gain_mathematica_audit.wl:84-93` vs `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage062_parent_action_gain_sympy_audit.py:91-95`

**What's wrong:**

The Mathematica script is a line-by-line port of the SymPy script:

| SymPy | Mathematica |
|---|---|
| `U_general = K * rho**n_poly / (n_poly - 1)` | `uGeneral = capitalK*rho^nPoly/(nPoly - 1);` |
| `h_general = sp.diff(U_general, rho)` | `hGeneral = D[uGeneral, rho]` |
| `hprime_general = sp.diff(h_general, rho)` | `hPrimeGeneral = D[hGeneral, rho]` |
| `cs_sq_general = sp.diff(K * rho**n_poly, rho) / m` | `csSqGeneral = (1/m)*D[capitalK*rho^nPoly, rho]` |
| `cs_sq_wrong = sp.diff(K * rho**(n_poly + 1), rho) / m` | `csSqWrong = (1/m)*D[capitalK*rho^(nPoly + 1), rho]` |
| `Theta_sigma = (m * cs_star_sq / rho_star) * Nss; Lambda_phi = g_phi * Osp` | `thetaSigma = (m*csStarSq/rhoStar)*nSS; lambdaPhi = gPhi*oSP` |
| `S_parent = (1/2)*Theta_sigma*sigma**2 + Lambda_phi*sigma*phi + (1/2)*KX*phi**2` | `sParent = (1/2)*thetaSigma*sigma^2 + lambdaPhi*sigma*phi + (1/2)*kX*phi^2` |
| `sigma_star = sp.solve(sp.diff(S_parent, sigma), sigma)[0]` | `sigmaStar = sigma /. First[Solve[D[sParent, sigma] == 0, sigma]]` |
| `gain_from_action = (KX - 2*S_eff_phi.coeff(phi, 2))/KX` | `gainFromAction = (kX - 2*Coefficient[sEff, phi, 2])/kX` |
| solve `Xi_micro - Xi_target = 0` for `kappa` and compare to `KX*L**2/TX` | `Solve[xiMicro == xiTarget, kappa]; compare to kX*ell^2/tX` |

The Mathematica script adds *one* auxiliary cross-check (`gainFromSeries` via `SeriesCoefficient`) on line 71 / asserted on line 79 (A9), but that compares `Coefficient[sEff, phi, 2]` against `SeriesCoefficient[Series[sEff, ...], 2]` for **the same `sEff` object** - an internal-consistency check between two algebraically equivalent extraction methods, not an independent re-derivation of the gain. Both engines walk the same path: same polytrope normalization, same on-shell elimination of `sigma` with the same `+ Lambda sigma phi` sign convention, same closed-form target.

**Why this matters:**

The framework's second-engine policy requires independent re-derivations so that a sign/branch/normalization error in one engine doesn't propagate to the other. Here the same sign error in F2 appears identically in both engines because the Mathematica script copied the SymPy script's `S_parent` expression rather than rederiving from `delta V_conf = -g_phi chi_phi phi` and `delta rho = sigma chi_sigma`. A real second-engine derivation would, for example: derive `G_micro` directly from `chi_sigma^(eff) * Lambda_phi^2/K_X` (the notes §4 chain) rather than via the quadratic-coefficient route; the SymPy script already does the coefficient route, so this gives Mathematica an independent algebraic path.

**Required change:**

Rewrite the Mathematica gain derivation so it follows the susceptibility chain, not the action-coefficient chain. Concrete edit (preferred):

1. After defining `thetaSigma` and `lambdaPhi` (line 63-64), introduce a `chiSigmaEff` directly:
   ```
   chiSigmaEff = 1/thetaSigma;
   ```
2. Compute the gain via the susceptibility chain (notes §4):
   ```
   gainViaSusceptibility = FullSimplify[chiSigmaEff*lambdaPhi^2/kX, Assumptions -> $Assumptions];
   ```
3. Keep the existing action-coefficient route (`gainFromAction`) and its internal-consistency `gainFromSeries` cross-check (the current A9), then add an *independent-route* check:
   ```
   expectZero["G_micro via susceptibility route vs closed form", gainViaSusceptibility - gClosed];
   expectZero["G_micro: action route equals susceptibility route", gainFromAction - gainViaSusceptibility];
   ```
4. The susceptibility route bypasses `Solve` and `Coefficient` entirely; it uses only `1/thetaSigma` and a product. This is a structurally distinct algebraic path from the SymPy script's action-coefficient route.

Keep the `Xi_micro` block unchanged (its `Solve` for `kappa` is already a real check).

**Verification:**

After Codex applies, the verifier will run the Mathematica script and confirm a new `PASS: G_micro via susceptibility route vs closed form` line appears in the transcript (alongside the existing `gMicro from parent action vs closed form`), and the script exits 0. The diff should show no introduction of new physical claims, only a second algebraic path to the same `gClosed` target.

## Independent-derivation check (Mathematica)

The `.wl` script does not currently independently rederive `G_micro` from the physical premises. As laid out in F3, its derivation is a syntactic port of the SymPy script's chain (`U -> h -> h' -> c_s^2`; `S_parent -> sigma_star -> S_eff -> coefficient -> gain`; `kappa` from `xiMicro = xiTarget`). The single `gainFromSeries` cross-check is an internal coefficient-extraction sanity, not an independent route.

## Engine cross-check

Both engines produce the same simplified outputs:

| Quantity | SymPy | Mathematica |
|---|---|---|
| `h'(rho)` (general) - `m c_s^2/rho` | 0 | 0 |
| `h'(rho)` at n=5 | `5*K*rho**3` | `5*capitalK*rho^3` |
| `c_s^2(rho)` at n=5 | `5*K*rho**4/m` | `(5*capitalK*rho^4)/m` |
| `Theta_sigma` | `N_ss*cs_star_sq*m/rho_star` | `(csStarSq*m*nSS)/rhoStar` |
| `Lambda_phi` | `O_sp*g_phi` | `gPhi*oSP` |
| `sigma_star` | `-O_sp*g_phi*phi*rho_star/(N_ss*cs_star_sq*m)` | `-((gPhi*oSP*phi*rhoStar)/(csStarSq*m*nSS))` |
| `G_micro` from action | `O_sp**2*g_phi**2*rho_star/(K_X*N_ss*cs_star_sq*m)` | `(gPhi^2*oSP^2*rhoStar)/(csStarSq*kX*m*nSS)` |
| `kappa` solved | `K_X*L**2/T_X` | `(ell^2*kX)/tX` |
| inconsistency probe residual | `K*rho**3*(5 - 6*rho)` | `capitalK*(5 - 6*rho)*rho^3` |

All outputs agree symbol-for-symbol after renaming. The agreement is not surprising given F3 (the engines share the same derivation), but at the level of the final identities they do match. `engines_agree: true` in the front matter reflects pure output agreement; the underlying issue is *how* they agree (F3), not *whether* they agree.

Output transcripts are both newer than their respective scripts (sympy script mtime 1779500305 < output 1779500439; wl script mtime 1779500366 < output 1779500444), so `outputs_fresh: true`.

## Verdict justification

The script's headline assertion - the first equality of the boxed `G_micro` paper formula - holds up: it is non-tautological (the gain is solved from a quadratic on-shell elimination, not pre-baked) and the n=5 specialization is anchored by a real inconsistency probe. The polytrope EOS identity is similarly real. The `kappa = K_X L^2/T_X` solve is non-trivial. So at the level of the load-bearing identity, the unit is verified.

What does not hold up: (i) the **second** equality of the boxed paper equation and the Cauchy-Schwarz bound on `C_{sigma phi}^2` are paper-side deliverables with no script-side check (F1, paper_misalignment / script_missing_paper_claim); (ii) the sigma-phi coupling sign in the reduced action `S_parent` disagrees with the notes' `F_red` (F2, paper_misalignment / notes_contradicts_script) - invariant for the gain magnitude but a real convention disagreement that could bite downstream consumers of `sigma_star(phi)`; (iii) the Mathematica audit is a syntactic transliteration of the SymPy audit (F3, mathematica_transliteration), so the two-engine policy is currently a no-op for this unit.

Attacks attempted and what they found:
- **Perturb the gain identity.** Replacing `Lambda_phi = g_phi * Osp` with `Lambda_phi = g_phi * Osp**2` would cause `gain_from_action - G_micro_closed` to evaluate to a non-zero residual, since `gain_from_action` is computed from `sp.solve` of `S_parent` and would now have an `Osp**4` factor against the hand-typed `G_micro_closed`. A4/A10 catches this. **Real check.**
- **Perturb the polytrope index.** Hardcoding `n_poly = 6` in `U_general` only (not in `cs_sq_general`) would cause the general-polytrope `expect_zero` to fail. The inconsistency probe A3/A8 demonstrates this catch path concretely. **Real check.**
- **Sign flip on Lambda_phi.** The gain assertion would still pass (quadratic in Lambda), so this is invisible to the script - this is the F2 issue.
- **Drop a factor on N_pp in the second equality.** The script doesn't carry the second equality, so this is invisible - this is the F1 issue.

Verdict: `findings` (3). Not `stop_cold` - none of these would invalidate the headline `G_micro` identity, and F1/F2 are paper-side resolutions for the user while F3 is a script-side rewrite Codex can apply.

## Self-test notes

For F1's proposed second-equality check, I mentally substituted `Osp**2/(Nss*Npp)` into the factored form and confirmed it reduces exactly to `G_micro_closed` - non-trivial in symbolic display but algebraically a single substitution, so the new `expect_zero` will fire correctly. For F3's susceptibility-route check `gainViaSusceptibility = (1/thetaSigma)*lambdaPhi^2/kX = (rhoStar/(m*csStarSq*nSS))*(gPhi*oSP)^2/kX`, this equals `(gPhi^2 * oSP^2 * rhoStar)/(csStarSq * kX * m * nSS) = gClosed`, so the new `expectZero` will pass. The directive for F3 names `mathematica/` for the `.wl` path. F2 is a `## Resolve before fix_loop` block with no Codex edit prescribed.
