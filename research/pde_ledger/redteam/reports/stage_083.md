---
unit_id: 083
batch: III.4
auditor_model: claude-opus-4-7-1m
audit_date: 2026-05-27T00:00:00-06:00
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
  notes_stage_files:
    - /var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage083_family1_direct_operator_window.md
  paper_appendix: present
---

# Audit unit 083 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_083.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage083_family1_direct_operator_window.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part03.tex` (Stage 083 included via `\input{stages/stage_083}` at line 284; no separate narrative row beyond the card)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage083_family1_direct_operator_window_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage083_family1_direct_operator_window_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage083_family1_direct_operator_window_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage083_family1_direct_operator_window_mathematica_audit.txt`

## What the paper claims

The stage card states (verbatim Output): "Direct operator-selected support window for Family-1." The Purpose line clarifies the stage "evaluates the Family-1 support window directly through the operator-selected branch rather than only through endpoint bounds," and the body equation (eq:app-stage083-zeta-direct) defines `zeta_F1^op = zeta_phys(Pe_*(Xi_F1, 37, 12321/5), 37; 12321/5)`. Inputs are `Xi_F1` (Family-1 wall/source strength) and the fixed-point equation for `Pe_*`. The supporting notes enumerate four concrete deliverables: (1) operator-selected `Pe_*` is monotone in `Xi` on the stable branch via `dPe_*/dXi = Delta/(1 - Xi*dDelta/dPe)`; (2) F1 support ratio and branch-product thresholds are monotone in `Xi`; (3) inserting the natural shell-weighted `Theta_w^(chi) ≈ 4.06863235008162` and Jensen `Theta_w^(J) ≈ 0.927552032539308` data reproduces the Stage-61/63/64 transport windows (`Pe_-/+^(chi)`, `Pe_-/+^(J)`) and resulting `zeta_-/+` values; (4) the natural branch sits within ~1.31e-3 of the hard ceiling `zeta_max^(F1) ≈ 2.46752922945601`. Notes also state the integer relation `Xi_F1 = 1369 Upsilon_w = 136900 Theta_w` (with `1369 = 37^2 = eta_F1^2`).

## What the script claims to verify

Both scripts share the same skeleton: (i) symbolic verification of the implicit-differentiation formula `dPe_*/dXi = Delta/(1 - Xi*dDelta/dPe)` via `Solve`/`solve`; (ii) defining-equation residuals for the closed-form `Delta_0(F1)` and `Delta_inf(F1)` against the linear identities `alpha^2*(alpha*sinh+eta*cosh)*Delta_0 - eta*(cosh-1) == 0` and `(alpha*sinh+eta*cosh)*Delta_inf - (cosh+(eta/alpha)*sinh - 1) == 0`; (iii) the eigenvalue residual `y*tan(y) - eta == 0` at the numerically solved `y_F1`; (iv) downstream numeric values for `Xi_chi`, `Xi_J`, the four `Pe_-/+` and `zeta_-/+` window endpoints, and the `zeta_max^(F1)` ceiling; (v) Mathematica adds a structural-identity check on `Omega(Pe)` against `Pi*Pe*(2*Pe*Exp[Pe]+Pi)`, eleven `expectApprox` numeric checks against the notes' stated decimals (e.g., `96.5285247264385`, `2.46622291347846`, `2.46752922945601`), and a sampled `d zeta / d Pe > 0` test at `Pe ∈ {10, 100, 1000, 10000}`. SymPy carries the same sampled monotonicity test and the same defining-equation residuals.

## Paper ↔ script cross-check

| Paper-side deliverable | Script-side check | Status |
|---|---|---|
| `dPe_*/dXi = Delta/(1 - Xi*dDelta/dPe)` (notes §1) | symbolic `expect_zero(dPe_formula - dPe_expected)` (sympy:47, math:49) | match |
| `Delta_0(F1) ≈ 1.733e-4`, `Delta_inf(F1) ≈ 2.014e-2` (notes §2) | defining-equation residual (sympy:71-81, math:77-91) + numeric `expectApprox` (math:169-170) against the notes' decimals | match |
| `y*tan(y) = eta` (notes §2, eigenvalue) | `y_F1` nsolve + residual `<1e-25` (sympy:84-89, math:93,100-101) | match |
| `Xi_F1 = 136900 Theta_w` (notes §2) | `Xi_chi = 136900 * Theta_chi_coeff` (sympy:108, math:121) with SOURCE-ANCHOR comment | match (documented carry-forward from Stage 60) |
| `Theta_w^(chi) ≈ 4.06863..., Theta_w^(J) ≈ 0.92755...` (notes §3 → Stage 60) | hardcoded literals with `SOURCE-ANCHOR` comment naming the upstream as a TODO (sympy:97-105, math:110-118) | match (documented; carried in via paper-anchored decimals) |
| `Pe_-^(chi) ≈ 96.528..., Pe_+^(chi) ≈ 11220.544..., Pe_-^(J) ≈ 22.006..., Pe_+^(J) ≈ 2558.018...` (notes §3) | computed + `expectApprox` against the notes' decimals (math:171-174) | match |
| `zeta_-/+^(chi)/(J)` values + `zeta_max^(F1) ≈ 2.46752922945601` (notes §4) | computed + `expectApprox` against the notes' decimals (math:175-179) | match |
| Monotonicity of `zeta` in `Pe` (and hence in `Xi` via stable branch) (notes §1-2) | sampled `d zeta / d Pe > 0` at four points (sympy:143-150, math:181-194) | partial (sample-based, not analytic) |
| `Pi/C_mix = 1 + zeta` at `eps_blk = 0` (notes §5) | computed and printed (sympy:158-162, math:196-200) | match (numerical print) |

The "partial" row is a known sample-based proxy for a global analytic monotonicity claim; the v1 directive's F3 explicitly chose this fix and the user accepted it. The samples span the full `Pe` range of the chi and J windows.

Set `paper_alignment: aligned`.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 47 | `expect_zero("implicit fixed-point derivative", dPe_formula - dPe_expected)` | notes §1 dPe_*/dXi formula | yes (algebraic cross-check via two routes through `Solve`) |
| A2 | sympy | 75 | `expect_zero("Delta_0(F1) defining-equation residual", ...)` | notes §2 Delta_0 closed form | yes (form-check against literal numerator/denominator) |
| A3 | sympy | 81 | `expect_zero("Delta_inf(F1) defining-equation residual", ...)` | notes §2 Delta_inf closed form | yes |
| A4 | sympy | 87-88 | `if abs(y_residual) > 1e-25: raise AssertionError(...)` | notes §2 eigenvalue eqn | yes (true numerical check against equation) |
| A5 | sympy | 147-150 | `if deriv_num <= 0: raise AssertionError(...)` (four samples) | notes §1-2 monotonicity claim | yes (sampled sign check on derivative) |
| A6 | mathematica | 49 | `expectZero["implicit fixed-point derivative", ...]` | notes §1 dPe_*/dXi formula | yes |
| A7 | mathematica | 81 | `expectZero["Delta_0(F1) defining-equation residual", ...]` | notes §2 Delta_0 closed form | yes |
| A8 | mathematica | 91 | `expectZero["Delta_inf(F1) defining-equation residual", ...]` | notes §2 Delta_inf closed form | yes |
| A9 | mathematica | 101 | `expectApprox["y_F1 satisfies y Tan[y] = eta", yRootResidual, 0, 10^-20]` | notes §2 eigenvalue eqn | yes |
| A10 | mathematica | 103 | `expectApprox["A_F1 independent vs closed-form", ...]` | notes §3 A_F1 = (kappa+pi^2/4)/(kappa+y^2) | weak (algebraically identical; near-tautological) |
| A11 | mathematica | 154 | `expectZero["Omega(Pe) identity residual", ...]` | notes §1 Omega(Pe) closed form | yes (form-check against literal product) |
| A12 | mathematica | 169-170 | `expectApprox[Delta_0/Delta_inf numeric check, ...]` | notes §2 numeric values | yes (target decimals are notes-side anchors) |
| A13 | mathematica | 171-174 | `expectApprox[Pe_-/+^(chi)/(J) numeric check, ...]` | notes §3 numeric values | yes (target decimals are notes-side anchors) |
| A14 | mathematica | 175-179 | `expectApprox[zeta_-/+^(chi)/(J), zeta_max^(F1) numeric check, ...]` | notes §3-4 numeric values | yes (target decimals are notes-side anchors) |
| A15 | mathematica | 186-194 | `If[!AllTrue[signs, # === 1 &], fail[...]]` (four samples) | notes §1-2 monotonicity claim | yes |

Every script-side assertion traces to a paper-side or notes-side claim. The eleven `expectApprox` checks on lines 169-179 of the Mathematica script (initially flagged in v1 F4 as tautological self-checks) are now correctly understood as anchors against the notes' stated decimals — if the script's computed value drifts from what the notes record, the check fires. A10 remains weak (algebraically identical expression) but does not undermine the overall verification.

## Findings

None. All four v1 findings (`mathematica_transliteration`, `hardcoded_result`, `insufficient_verification`, `tautological_check`) have been resolved per the v1 directive's `Applied: F1`–`F4` blocks. The script content was inspected to confirm the listed insertions are present at the named line ranges.

## Independent-derivation check (Mathematica)

The Mathematica script is no longer purely transliterated. Independent (or quasi-independent) checks now present:
- `Delta_0(F1)` and `Delta_inf(F1)` defining-equation residuals (lines 77-91) verify the closed form satisfies its linear defining identity. The literal RHS of the residual provides an anchor for catching transcription errors in the closed-form body.
- `y_F1*Tan[y_F1] - etaF1` numerical check (line 101) is a true independent verification of the eigenvalue equation, not a transliteration.
- `Omega(Pe)` identity residual (line 154) form-checks the closed form against the literal product `Pi*Pe*(2*Pe*Exp[Pe]+Pi)`. A typo of the form `Exp[Pe] -> Exp[2*Pe]` in the body would produce a nonzero residual.
- Numeric `expectApprox` targets on lines 169-179 are the notes' stated decimals, providing an independent anchor.

The closed forms themselves are still typed in (not re-derived from the BVP), and the v1 directive acknowledged this. The current state is the user-accepted resolution.

## Engine cross-check

Both engines compute the same numerical values to ~14-15 digits (SymPy `Float`/`nsolve` precision is the bottleneck). The Mathematica `expectApprox` blocks pass with all diffs at or below `10^-10`. The SymPy script's monotonicity samples show `dzeta/dPe > 0` at all four test points, matching the Mathematica `Sign === 1` result. No engine disagreement.

| Quantity | SymPy | Mathematica |
|---|---|---|
| `Delta_0(F1)` | 0.000173302079021525149057156196550 | 0.000173302079021525149057156196549924... |
| `Delta_inf(F1)` | 0.0201447565540521594271032956099 | 0.0201447565540521594271032956099177... |
| `y_F1` | 1.52948248371469964992710762240 | 1.5294824837146996499271076224024730256... |
| `A_F1` | 1.00005192880219532865933408371 | 1.0000519288021953286593340837139871812... |
| `Pe_-^(chi)` | 96.5285247264385161753086456051 | 96.5285247264385 |
| `Pe_+^(J)` | 2558.01892349205354801313623259 | 2558.0189234920535 |
| `zeta_max^(F1)` | 2.46752922945601223332958450157 | 2.4675292294560122333295845015705354203... |

## Verdict justification

`clean`. The v1 audit raised four substantive findings, all of which have been applied (per the `Applied: F1`–`F4` blocks in the v1 directive) and the resulting scripts pass under reasoning attack: defining-equation residuals catch closed-form transcription errors against literal-RHS anchors; the `y*tan(y) - eta` check is a real numerical residual; the `Omega(Pe)` identity residual catches body-vs-product mismatches; the sampled monotonicity tests would catch a sign error in `Omega`; the `expectApprox` numeric targets on lines 169-179 are anchored to the notes' stated decimals (verified by inspection: `96.5285247264385`, `2.46622291347846`, `2.46752922945601`, etc., all match notes section 3/4). The paper card and notes describe the same windows the scripts compute, with the same constants and same Family-1 parameters (`kappa = 12321/5`, `eta = 37`). No `paper_misalignment`, no new `tautological_check`, no `hardcoded_result` beyond what the v1 directive already documented via SOURCE-ANCHOR comments.

Attacks tried that did not surface findings:
- Paper-side numeric audit: notes' `Theta_w^(chi) ≈ 4.06863235008162` and `Theta_w^(J) ≈ 0.927552032539308` match the script literals exactly; the integer factor `136900 = 370^2 = (10*37)^2 = 100*1369 = 100*eta_F1^2` is consistent with notes' `Xi_F1 = 1369 Upsilon_w = 136900 Theta_w`.
- Branch-cut hunt: `cosh`, `sinh`, `exp`, `sqrt(kappa) = sqrt(12321/5)` are all single-valued at the relevant points; no branch issues.
- Symbol-domain check: `xi > 0`, `kappa, eta, alpha > 0` in Mathematica's `$Assumptions`; SymPy's `Xi, Pe_sym` are `positive=True, real=True`. No simplify under contradictory assumptions.
- Output freshness: SymPy script (May 25 00:27) < output (May 25 00:29); Mathematica script (May 25 00:26) < output (May 25 00:30). Both fresh; no `stale_output`.
- Engine disagreement: numeric values agree to SymPy precision; no disagreement.
- Variable-independence trap: derivative samples are at concrete numeric `Pe` values applied to `dzeta_dPe = sp.diff(zeta_F1, Pe_sym)`, where `zeta_F1` genuinely depends on `Pe_sym` (via `Omega`). Non-trivial.
- Parity trap: no integrals over symmetric unbounded domains in this stage; not applicable.

The script's coverage of the paper claim is now substantive across all four deliverables in the notes. No directive is needed.

## Self-test notes

Checked: (1) variable-independence — the monotonicity-derivative samples target `Pe_sym`, which is the binding variable of `zeta_F1`'s `Omega`; the derivative is non-trivially nonzero at the sample points. (2) Parity — no symmetric-domain integrals. (3) Trivial-case — substituting `eta = 0` into `Delta_0` gives 0 (consistent: no source means no support extraction); substituting `Pe -> 0` into `zeta_F1` gives `A_F1` (limit of `Omega -> 1`), positive and well below the ceiling. (4) Path specs — verdict is `clean` so no missing-script directive paths to validate. (5) Paper round-trip — the script's `Theta_w` literals, `Xi_F1 = 136900*Theta_w` factorization, and `zeta_max ≈ 2.46752922945601` all match the notes' stated values exactly.
