---
unit_id: 151
batch: IV.6
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-08T00:00:00Z
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
  notes_stage_files: [/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage151_first_order_selected_correction.md]
  paper_appendix: present
---

# Audit unit 151 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_151.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage151_first_order_selected_correction.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part04.tex` (only an `\input{stages/stage_151}` include at line 1336; no dedicated status row)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage151_first_order_selected_correction_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage151_first_order_selected_correction_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage151_first_order_selected_correction_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage151_first_order_selected_correction_mathematica_audit.txt`

## What the paper claims

Stage 151 projects the exact full-profile mouth residual \(R_*(x)\) onto the Family-1 rigidity formulas to obtain the **actual** first-order non-exponential source correction. The card's `\stagefield{Purpose}` is terse ("a finite mouth-profile corrections ledger step. Its audit target is the verification output quoted below") and the in-card claim is the quote "Actual selected deformation is the centered residual \(R_*\); correction is fixed by two covariances." The notes are authoritative on intent and enumerate the deliverables: (1) the normalized first-order correction \(\Sigma_{\rm act}=\Sigma_*[1-\widetilde R_*]\) with the centered residual \(\widetilde R_*=R_*-\langle R_*\rangle_*\) (§1); (2) the two moment shifts \(\delta\mathfrak g_{\rm act}=-\operatorname{Cov}_*(c,R_*)\) and \(\delta\mathcal S_{\rm act}=-\operatorname{Cov}_*(K_q,R_*)\) (§2); (3) the bias retuning \(\delta\Pi_{\rm act}=\operatorname{Cov}_*(c,R_*)/\mathfrak g_*'\) and the traction retuning \(\delta\widehat T_{m,{\rm act}}=-A_T\operatorname{Cov}_*(c,R_*)-B_T\operatorname{Cov}_*(K_q,R_*)\) (§3). The kernels are fixed: \(c(x)=\cos(\pi x/2)\), \(K_q(x)=\cosh(\tfrac{\pi}{2}(1-x))/\cosh(\pi/2)\). No numeric values are pinned at this stage — the notes explicitly defer "actual numerical covariances on the explicit Family-1 branch" to a downstream stage.

## What the script claims to verify

Both scripts build the normalized full mouth source \(\Sigma_{\rm full}(x)=e^{-\Pi_* x-\epsilon(r_1 x+r_2 x^2)}/Z\), expand to first order in \(\epsilon\), and extract \(\Sigma_*\) (order 0) and \(\delta\Sigma\) (order 1). They then assert: the centered-residual identity \(\delta\Sigma + \Sigma_*(R-\langle R\rangle)=0\) (M1), normalization \(\int\Sigma_*=1\) (M2, sympy only), centering \(\int\delta\Sigma=0\) (M3 / Mathematica's "centering" check), the two moment-shift identities \(\delta g_{\rm int}+\operatorname{Cov}(c,R)=0\) (M4) and \(\delta S_{\rm int}+\operatorname{Cov}(K,R)=0\) (M5), and the two retuning identities \(\delta\Pi-\operatorname{Cov}(c,R)/g'=0\) (M6) and \(\delta T+A_T\operatorname{Cov}(c,R)+B_T\operatorname{Cov}(K,R)=0\) (M7). \(R(x)=r_1 x+r_2 x^2\) is a generic positive-family residual ansatz; \(r_1,r_2,A_T,B_T,g'\) stay symbolic in both engines. SymPy verifies M1–M7 exactly at five rational \(\Pi_*\in\{1/2,1,3/2,2,5/3\}\); Mathematica verifies the same identities once, fully symbolic in \(\Pi_*>0\).

## Paper ↔ script cross-check

| Paper deliverable | Script check | Status |
|---|---|---|
| §1 \(\Sigma_{\rm act}=\Sigma_*[1-\widetilde R_*]\), \(\widetilde R_*=R_*-\langle R_*\rangle_*\) | M1 (both engines): `delta_Sigma + Sigma_star*(R - <R>) == 0` | match |
| §1 normalization \(\langle 1\rangle_*=1\) | M2 (sympy): `int Sigma_star - 1 == 0` | match |
| §1 centering \(\int\delta\Sigma=0\) | M3 (sympy) / `<deltaSigma>_* (centering)` (wl) | match |
| §2 \(\delta\mathfrak g_{\rm act}=-\operatorname{Cov}_*(c,R_*)\) | M4: `dg + CovcR == 0` | match |
| §2 \(\delta\mathcal S_{\rm act}=-\operatorname{Cov}_*(K_q,R_*)\) | M5: `dS + CovKR == 0` | match |
| §3 \(\delta\Pi_{\rm act}=\operatorname{Cov}_*(c,R_*)/\mathfrak g_*'\) | M6: `deltaPi - CovcR/gprime == 0` | match |
| §3 \(\delta\widehat T_{m,{\rm act}}=-A_T\operatorname{Cov}_*(c,R_*)-B_T\operatorname{Cov}_*(K_q,R_*)\) | M7: `deltaT + A_T*CovcR + B_T*CovKR == 0` | match |

Every paper-side boxed deliverable maps to a substantive script-side assertion in both engines (M2 sympy-only, but the equivalent normalization is implicit in Mathematica's `SigmaFull = unnorm/Z`). `paper_alignment: aligned`.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored? |
|---|---|---|---|---|---|
| A1 | sympy | 153–156 | `expect_zero(delta_Sigma + Sigma_star*(pert - Rbar))` | §1 centered-residual law | yes |
| A2 | sympy | 157–160 | `expect_zero(int Sigma_star - 1)` | §1 normalization | yes |
| A3 | sympy | 161–164 | `expect_zero(int delta_Sigma)` | §1 centering | yes |
| A4 | sympy | 168–171 | `expect_zero(dg + CovcR)` | §2 \(\delta\mathfrak g\) | yes |
| A5 | sympy | 172–175 | `expect_zero(dS + CovKR)` | §2 \(\delta\mathcal S\) | yes |
| A6 | sympy | 179–182 | `expect_zero(deltaPi - CovcR/gprime)` | §3 \(\delta\Pi\) | yes (algebraically downstream of A4) |
| A7 | sympy | 183–186 | `expect_zero(deltaT + AT*CovcR + BT*CovKR)` | §3 \(\delta\widehat T_m\) | yes (lin. comb. of A4,A5) |
| B1 | wl | 53 | `expectZero(int deltaSigma)` | §1 centering | yes |
| B2 | wl | 54 | `expectZero(deltaSigma + SigmaStar*(R-rBar))` | §1 centered-residual law | yes |
| B3 | wl | 58 | `expectZero(deltaGInt + covCR)` | §2 \(\delta\mathfrak g\) | yes |
| B4 | wl | 59 | `expectZero(deltaSInt + covKR)` | §2 \(\delta\mathcal S\) | yes |
| B5 | wl | 63 | `expectZero(deltaPi - covCR/gPrime)` | §3 \(\delta\Pi\) | yes (downstream of B3) |
| B6 | wl | 64 | `expectZero(deltaT + aT*covCR + bT*covKR)` | §3 \(\delta\widehat T_m\) | yes (lin. comb. of B3,B4) |

Note on A6/A7 (and B5/B6): these are algebraically downstream of A4/A5 (M6 = −(dg+CovcR)/g'; M7 = A_T(dg+CovcR)+B_T(dS+CovKR)). They are NOT tautological X−X checks — `deltaPi`/`deltaT` are built from the independently-integrated `dg`/`dS`, while `CovcR`/`CovKR` are built from the separate covariance route, so each can fail if the two routes disagree. They are redundant with M4/M5 but they each correspond to a distinct boxed paper deliverable (\(\delta\Pi\), \(\delta\widehat T_m\)), so retaining them is correct and not a finding.

## Findings

None.

The two strongest attack vectors were checked and held:

1. **Tautology on M4/M5 (the load-bearing covariance identities).** `dg = ∫ c·delta_Sigma dx` is computed from the eps-series-derived, renormalized `delta_Sigma` (sympy lines 138–141 via `num/(Z0+eps*Z1)`; wl lines 35–37 via `Series[SigmaFull]`). `CovcR = mean(c·pert) − cbar·Rbar` is computed independently from `Sigma_star`-weighted moments (sympy lines 146–151). The two are built from structurally different expressions, so `dg + CovcR == 0` is a real identity, not a definitional rearrangement. It can fail (e.g., if the normalization derivative term `−Z1/Z0²` in `delta_Sigma` were dropped, M4 would break while M3 might still pass). Not tautological.

2. **mathematica_transliteration.** The engines are genuinely independent, not a line-by-line port. SymPy installs a bespoke recursive integration-by-parts integrator (`_poly_exp_moment`, `_expand_linear_exponentials`, `_integrate_exp_poly_term`, lines 55–119) that rewrites cos/cosh to exponentials and integrates exp×polynomial over [0,1] in closed form, then runs it at five rational `Pi_star` samples. Mathematica keeps `piStar` fully symbolic and calls native `Integrate[..., Assumptions -> piStar > 0]` once. Different integration machinery, different Pi_star treatment (5 rationals vs. all-Pi_star symbolic), different series plumbing. This is the deliberate dual-engine split, not transliteration.

## Independent-derivation check (Mathematica)

Genuinely independent. Mathematica derives \(\Sigma_*\)/\(\delta\Sigma\) from `Series[Exp[-Phi[x]]/Z, {epsilon,0,1}]` with `Z = Integrate[Exp[-Phi[x]], {x,0,1}]` evaluated symbolically in `piStar` (wl lines 31–37). SymPy cannot do this symbolic integral (it hangs — confirmed by the anti-footgun comment, lines 13–23) and instead uses its own IBP-based moment routine at concrete rational `Pi_star` (lines 55–119, 130–141). The covariance assembly is conceptually the same (it must be — it's the same physics), but the underlying integration is done by completely different code paths. No shared intermediate expressions are imported across engines; each computes the kernels' moments from scratch. This is the correct dual-engine design, and the anti-footgun comment forbidding re-attempts at symbolic-Pi_star SymPy integration is honestly documented and is not itself a defect.

## Engine cross-check

Both engines assert the same seven identities (Mathematica omits the explicit `∫Sigma_star=1` because normalization is enforced by construction in `SigmaFull = unnorm/Z`) and both report all checks `= 0` / `PASS` with `exit_code: 0`. SymPy: M1–M7 each `= 0` at all five `Pi_star` samples (sympy output lines 6–49). Mathematica: all six `expectZero` checks `= 0` / `PASS` symbolically (wl output lines 5–16). The five SymPy rational samples are exact special cases of the Mathematica symbolic result, so the two engines agree at the level claimed (exact-zero residuals). No disagreement.

## Verdict justification

Clean. I read the paper card, the notes (authoritative on intent), and the appendix include, then attacked the scripts. Every boxed deliverable in the notes (§1 centered-residual law + normalization + centering; §2 the two covariance moment shifts; §3 the bias and traction retunings) maps to a substantive, non-tautological assertion in both engines, verified exactly. The custom SymPy integrator is a correct recursive-IBP implementation, not a shortcut that bakes in the answer; the cos/cosh kernels are genuinely integrated against the exponential weight. The two engines are independent (custom rational-sample integrator vs. native symbolic Integrate), so this is not a transliteration. No hardcoded numeric answers, no `168π²`/`100π²` or `√(4107−100π²)/(10π)` constants appear (correctly — this stage is symbolic and defers numerics downstream). No symbol-assumption errors (`piStar>0`, reals, `gprime≠0` all physically justified). Outputs are newer than their scripts (fresh). Attacks tried and failed: (a) M4/M5 tautology — refuted, the two sides are built from distinct expressions; (b) transliteration — refuted, independent integration code paths; (c) hidden hardcoded/stale constants — none found; (d) symbol-domain abuse — none.

## Self-test notes

Checked: (1) Variable independence — no `sp.diff(EXPR, VAR)`/`D[expr,var]` with VAR absent from EXPR; the only derivatives are inside the integrator's `sp.diff(exponent, x)` (line 95) which is by construction linear-in-x, so `slope` is x-free as the code's own guard requires. (2) Parity/symmetry — integrals are over the finite [0,1], not unbounded, so the odd/even-vanishing trap does not apply; vanishing of `∫delta_Sigma` is enforced by the renormalization, not by parity, and is correctly the centering check. (3) Trivial-case pre-check — substituting `r1=r2=0` (zero residual) makes `delta_Sigma=0`, `CovcR=CovKR=0`, so all of M1,M3–M7 reduce to `0=0`; the non-triviality comes from generic symbolic `r1,r2`, which the script retains. (4) Path specs — N/A (no missing-script findings). (5) Paper round-trip — N/A (no fixes prescribed).

## Value Reconciliation (pass-2 augmentation)

This stage is intentionally **symbolic**: the scripts pin no closed-form numeric figure-of-merit. The deliverables are symbolic identities (the boxed equations of notes §1–§3) over the generic residual ansatz \(R(x)=r_1 x + r_2 x^2\), with `r1, r2, A_T, B_T, gprime` and (for SymPy) `Pi_star` left as free parameters. The notes explicitly defer all "actual numerical covariances on the explicit Family-1 branch" to a downstream stage, so no Family-1 numeric (e.g., a `Pi_star ≈ 1.508…` or `√(4107−100π²)/(10π)` radius) is expected here, and none appears — correctly.

The only emitted RESULT values are the asserted symbolic identities (each printed as `... = 0`), which are the verification scaffolding for the boxed deliverables rather than standalone reportable numbers. Each maps to a boxed equation in the notes:

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| \(\delta\Sigma=-\Sigma_*\widetilde R_*\) identity (M1/B2 residual `=0`) | py L153–156 / py.txt L7; wl L54 / wl.txt L7 | notes §1 boxed `\Sigma_{\rm act}=\Sigma_*[1-\widetilde R_*]` (md L28–35) | MATCH |
| normalization \(\langle1\rangle_*=1\) (M2 `=0`) | py L157–160 / py.txt L8 | notes §1 `\langle f\rangle_*:=\int_0^1\Sigma_* f` (md L38–40) | MATCH (def. carrier) |
| centering \(\int\delta\Sigma=0\) (M3/B1 `=0`) | py L161–164 / py.txt L9; wl L53 / wl.txt L5 | notes §1 \(\widetilde R_*=R_*-\langle R_*\rangle_*\) (md L34) | MATCH |
| \(\delta\mathfrak g=-\operatorname{Cov}_*(c,R)\) (M4/B3 `=0`) | py L168–171 / py.txt L10; wl L58 / wl.txt L9 | notes §2 boxed (md L57–65) | MATCH |
| \(\delta\mathcal S=-\operatorname{Cov}_*(K_q,R)\) (M5/B4 `=0`) | py L172–175 / py.txt L11; wl L59 / wl.txt L11 | notes §2 boxed (md L66–74) | MATCH |
| \(\delta\Pi=\operatorname{Cov}_*(c,R)/\mathfrak g_*'\) (M6/B5 `=0`) | py L179–182 / py.txt L12; wl L63 / wl.txt L13 | notes §3 boxed (md L95–102) | MATCH |
| \(\delta\widehat T_m=-A_T\operatorname{Cov}_*(c,R)-B_T\operatorname{Cov}_*(K_q,R)\) (M7/B6 `=0`) | py L183–186 / py.txt L13; wl L64 / wl.txt L15 | notes §3 boxed (md L106–116) | MATCH |

INTERNAL scaffolding (accounted for, no finding): `Pi_star` sample set `{1/2,1,3/2,2,5/3}`, `Z0`/`Z1` normalization partials, `Sigma_star`/`delta_Sigma` series coefficients, `cbar`/`Kbar`/`Rbar` moments, `CovcR`/`CovKR` covariance intermediates, the `_poly_exp_moment` recursion values, PASS/FAIL flags, `exit_code`.

reconciliation: complete; 7 deliverable identities checked, 0 misaligned.
