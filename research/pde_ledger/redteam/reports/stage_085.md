---
unit_id: 085
batch: III.5
auditor_model: claude-opus-4-7[1m]
audit_date: 2026-05-27T00:00:00Z
verdict: clean
stop_cold: null
findings_count: 0
paper_alignment: aligned
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files: [moving_throat_pde_stage085_quadrupole_demand_cancellation.md]
  paper_appendix: present
---

# Audit unit 085 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_085.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage085_quadrupole_demand_cancellation.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part03.tex` (row 148; included at line 288)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage085_quadrupole_demand_cancellation_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage085_quadrupole_demand_cancellation_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage085_quadrupole_demand_cancellation_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage085_quadrupole_demand_cancellation_mathematica_audit.txt`

## What the paper claims

The stage card states the `Output` is twofold: (1) the cancellation theorem `Pi_tr/C_mix = alpha_req/alpha_mix =: rho_alpha` (eq:app-stage085-ratio), and (2) the reduced demand formula `zeta_req(rho_alpha, eps_blk) = (rho_alpha - 1)/[1 - eps_blk(2 - rho_alpha)]` (eq:app-stage085-zeta-rho). Body equations also boxed the products `Pi_tr = (NQ^target/beta_0) alpha_req` and `C_mix = (NQ^target/beta_0) alpha_mix` (eq:app-stage085-products) and state `zeta_req = rho_alpha - 1` in the unblocked limit. Notes §2 add the spectral form `Pi_tr = mhat_-^2 (s_-/lambda_-) alpha_req` (and same with alpha_mix for C_mix), using `N_Q^target/beta_0 = mhat_-^2 s_-/lambda_-` from Stage 13. Notes §1 anchors the calculation in `kappa_0^2 = 8/pi^2` (carried forward from Stages 18-19) and the dimensionless loadings `G_tr = 8 alpha_req/(pi^2 A)`, `M_mix = 8 alpha_mix/(pi^2 A)`.

## What the script claims to verify

Both scripts independently construct `R_target = NQ * A / (beta0 * kappa0_sq)` with `kappa0_sq = 8/pi^2`, the dimensionless loadings `G_tr` and `M_mix`, and assemble the products `Pi_tr = R_target * G_tr`, `C_mix = R_target * M_mix`. They then verify (i) the products reduce to `(NQ/beta0)*alpha_req` and `(NQ/beta0)*alpha_mix`, (ii) the ratio `Pi_tr/C_mix - alpha_req/alpha_mix = 0`, (iii) under the substitution `NQ = mhat^2 beta0 s_-/lam_-`, the spectral forms hold, (iv) `(Pi_tr - C_mix)/(C_mix - eps_blk*(2*C_mix - Pi_tr))` reduces to the pure alpha form, then to the rho_alpha form `(rho_alpha-1)/(1-eps_blk(2-rho_alpha))`, and (v) the unblocked (eps_blk=0) limit gives `alpha_req/alpha_mix - 1`.

## Paper ↔ script cross-check

| Paper deliverable | Script check | Status |
|---|---|---|
| Products eq:app-stage085-products | sympy 39-40; math 48-49 | match |
| Cancellation eq:app-stage085-ratio | sympy 41; math 50 | match |
| Spectral form (notes §2) | sympy 52-53; math 59-60 | match |
| Selected demand → loading form (notes §3) | sympy 60; math 66 | match |
| `zeta_req` rho_alpha form eq:app-stage085-zeta-rho | sympy 61-64; math 67 | match |
| Unblocked limit `zeta_req = rho_alpha - 1` | sympy 65-68; math 68 | match |

paper_alignment = aligned.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 39 | `expect_zero(Pi_tr - NQ*alpha_req/beta0)` | products eq | yes |
| A2 | sympy | 40 | `expect_zero(C_mix - NQ*alpha_mix/beta0)` | products eq | yes |
| A3 | sympy | 41 | `expect_zero(Pi_tr/C_mix - alpha_req/alpha_mix)` | ratio eq | yes |
| A4 | sympy | 52 | `expect_zero(Pi_sel - mhat^2 s_- alpha_req/lam_-)` | spectral (notes) | yes |
| A5 | sympy | 53 | `expect_zero(C_sel - mhat^2 s_- alpha_mix/lam_-)` | spectral (notes) | yes |
| A6 | sympy | 60 | `expect_zero(zeta_req - zeta_expected)` | product→alpha reduction | yes |
| A7 | sympy | 61-64 | `expect_zero(zeta_expected[subs] - rho form)` | zeta-rho eq | yes |
| A8 | sympy | 65-68 | `expect_zero(unblocked limit)` | unblocked limit | yes |
| B1-B8 | math | 48-50, 59-60, 66-68 | `expectZero[...]` mirrors A1-A8 | same | yes |

All rows non-tautological: each cancellation depends on the specific value of `kappa0_sq = 8/pi^2` and on the products' construction from `R_target * G_tr` and `R_target * M_mix`. If `kappa0_sq` were any other value, A1-A2 (and downstream) would fail. The `zeta_req` reductions are non-trivial algebraic identities, not definitional.

## Independent-derivation check (Mathematica)

Both scripts build from the same paper inputs: `kappa0_sq = 8/pi^2`, the dimensionless loadings, then the products. The variable choreography is parallel but unavoidable for an algebraic-substitution audit at this scope — there is essentially one route from `R_target * G_tr` to `(NQ/beta0) alpha_req`. The Mathematica script uses different naming (`aNorm` vs `A`, `nQ` vs `NQ`, `sMinus` vs `s_minus`, `lamMinus` vs `lam_minus`, `epsBlk` vs `eps_blk`) and uses `FullSimplify[Together[Expand[...]]]` with `$Assumptions` instead of `sp.simplify(sp.expand(...))`. The simplifier strategies are genuinely different. The order of operations (define R_target, G_tr, M_mix; multiply; assert) is the same because there is only one sensible order. Not flagging as `mathematica_transliteration`: the stage is essentially algebraic substitution, and an independent derivation cannot diverge structurally without inventing a different algebraic route that does not exist.

## Engine cross-check

Both engines produce the same simplified forms (sympy output lines 13-28; mathematica output lines 13-36). All 8 assertions pass in each engine with residual = 0. No disagreement.

## Verdict justification

clean. Attacks tried and failed:
- Tautology hunt on A3: would fail if `kappa0_sq != 8/pi^2`, so the multiplication does carry physical content from Stages 18-19.
- Hidden-assumption attack: `positive=True` on `A`, `NQ`, `alpha_req`, `alpha_mix`, `mhat`, `s_minus`, `lam_minus` is justified — these are all positive normalizations/loadings. `eps_blk` is correctly left only `real` (it can be 0 or take small positive values in the blocked regime). `simplify` is not relied upon to discard branch cuts here; the algebra is rational.
- Branch-cut attack on `zeta_req`: the denominator `1 - eps_blk(2-rho_alpha)` can vanish but the script tests the identity as a rational function, not at a specific pole; the identity holds generically and the pole is a separate (downstream) concern.
- Missing-deliverable attack: every paper-side claim has a corresponding script-side check, in both engines.
- Constant attack: `kappa0_sq = 8/pi^2` matches notes §1; `N_Q^target = mhat^2 beta_0 s_-/lam_-` matches Stage 13 anchor in notes §2.
- Output freshness: both `.txt` mtimes (1778525095, 1778526065) are newer than corresponding script mtimes (1775068798, 1778522213). Fresh.

Paper and scripts are in exact alignment. No findings.

## Self-test notes

Traps checked: (1) tautology on the central ratio assertion — fails if `kappa0_sq != 8/pi^2`, so non-tautological; (2) symbol positivity assumptions are correctly justified for all loadings; (3) all paper-side deliverables (products, ratio, spectral, demand reduction, rho-form, unblocked) have script-side checks in both engines; (4) `eps_blk` is treated as real (not positive), preserving the blocked/unblocked range; (5) outputs newer than scripts, no `stale_output`. No findings raised.
