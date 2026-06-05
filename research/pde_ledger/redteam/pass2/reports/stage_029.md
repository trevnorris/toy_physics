---
unit_id: 029
batch: II.1
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-04T00:00:00Z
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
  notes_stage_files: [/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage029_dynamic_loading.md]
  paper_appendix: present
---

# Audit unit 029 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_029.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage029_dynamic_loading.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part02.tex` (row 48 references this stage)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage029_dynamic_loading_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage029_dynamic_loading_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage029_dynamic_loading_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage029_dynamic_loading_mathematica_audit.txt`

## What the paper claims

Stage 029 ("Coupled Response Operator", anchor MTDC-T5) gives the Stage-028 loading coefficient `alpha` a microscopic origin. Its `\stagefield{Output}` reads: *"Stage~029 outputs the Schur-complement split \eqref{eq:app-stage029-schur-split}, the static branch data \eqref{eq:app-stage029-static-data}, the outgoing transfer coefficient \eqref{eq:app-stage029-beta0}, and the selected odd coefficient \eqref{eq:app-stage029-selected-odd}."* The distinct deliverables are: (1) the exact Schur-complement split `Sigma_wall(omega) = Xi(omega) I_2 + alpha(omega) vv^T` with `Xi = lambda_U^2/A_U` and the boxed `alpha(omega)` rational form (eq. Xi-alpha); (2) the conservative static data `Xi_0 = lambda_U^2/Omega_U^2`, `Delta_0 = Omega_U^2 Omega_W^2 - lambda_R^2 sigma`, and `alpha_0` (eq. static-data); (3) the outgoing transfer coefficient `beta_0 = (Omega_U^2 lambda_W + lambda_R lambda_U)^2 / Delta_0^2 >= 0` (eq. beta0), with the first-order dressing `alpha = alpha_cons + beta Pi_out + O(Pi_out^2)`; (4) the selected-mode odd coefficient `delta D_-^odd = -i (a^5/27 c_s^5) beta_0 (v.e_-)^2 omega^5 + O(omega^7)` (eq. selected-odd). The notes add: the constant `sigma = 88/(9 pi^2)` (= Stage-027 max-coupling), the refined angle law `tan(2 theta_-) = 2 alpha_0 kappa_0 kappa_1 / (Delta K_ax + alpha_0 (kappa_0^2 - kappa_1^2))`, the limiting identities `(v.e_-)^2 = kappa(theta_-)^2`, and `beta_5 = Gamma_2^port beta_0` with `Gamma_2^port = a^5/(27 c_s^5)`.

## What the script claims to verify

Both engines verify the same six chained results. The SymPy docstring states it checks: the full Schur-complement elimination of the coupled wall/U/W/phi block, the exact `Xi I_2 + alpha vv^T` decomposition, the conservative static limits `Xi_0` and `alpha_0`, the refined angle law with the isotropic shift, the first-order outgoing expansion of `alpha(omega)`, and the selected-mode odd coefficient projected onto the conservative lower eigenvector. Concretely the load-bearing assertions are: (A) the matrix `Sigma = C^T M_int^{-1} C` equals `Xi I_2 + alpha vv^T` (sympy line 126 / mathematica line 116 via independent sequential elimination); (B) the overlap constants `sigma=88/(9 pi^2)`, `xi=56/(9 pi^2)`, `eta=-8 sqrt(2)/(3 pi^2)` (mathematica lines 117-119); (C) the stationarity / angle-law identity (sympy 183 / math 148,149); (D) `beta = beta_clean` (sympy 208 / math 166) and the `O(Pi)` expansion (sympy 212 / math 170); (E) `beta_5` extraction from the `omega^5` series coefficient (sympy 224 / math 179,180); (F) Hellmann-Feynman `kappa_sel^2` agreeing with a direct eigenvector projection, plus the al→0 (→kappa0^2) and al→∞ (→sigma) limits, and the paper-formula selected-odd identity (sympy 282-303 / math 214-232).

## Paper ↔ script cross-check

| Paper deliverable | Script-side check | Status |
|---|---|---|
| (1) Schur split `Sigma = Xi I_2 + alpha vv^T`, `Xi=lambda_U^2/A_U`, boxed `alpha(omega)` | sympy 119-126 (`C^T M_int^{-1} C` vs `Xi I + alpha vv^T`); math 90-116 (independent sequential elimination) | match |
| (2) static data `Xi_0`, `Delta_0`, `alpha_0` | sympy 140-148 prints all three; math 121-129 | match |
| (3) `beta_0 = P_0^2/Delta_0^2 >= 0` and `alpha = alpha_cons + beta Pi + O(Pi^2)` | sympy 199-212 (`beta - beta_clean=0`, `O(Pi)` residual=0); math 159-170; `beta_5 = Gamma_port beta_0` math 174-179 | match |
| (4) selected odd `delta D_-^odd = -i Gamma beta_0 (v.e_-)^2 omega^5` | sympy 291-303 (`delta_D_script - delta_D_paper=0`); math 222-232 | match |
| sigma=88/(9 pi^2) (notes/tex) | math 117 `sigma - 88/(9 Pi^2)=0`; sympy prints 88/(9*pi**2) | match |
| tan(2 theta_-) angle law (notes) | sympy 180-186 stationarity + tan2theta; math 142-157 | match |
| `(v.e_-)^2 = kappa(theta_-)^2`, limits kappa0^2 / sigma (notes) | sympy 281-289; math 203,216,217 | match |

`paper_alignment: aligned` — every paper-side deliverable maps to a non-tautological script-side check; no script-side check is orphaned from the paper.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored? |
|---|---|---|---|---|---|
| A1 | sympy | 126 | `expect_zero(Sigma - (Xi I + alpha vv^T))` | claim 1 (Schur split) | yes |
| A2 | sympy | 154 | `expect_zero(DeltaK_tilde - DeltaK_ax)` | none (self-labeled typo guard) | no (acknowledged tautology) |
| A3 | sympy | 183 | `expect_zero(dE/dtheta - stationarity/2)` | angle law (claim 2) | yes |
| A4 | sympy | 208 | `expect_zero(beta - beta_clean)` | claim 3 (transfer factor) | yes |
| A5 | sympy | 212 | `expect_zero(alpha_series - (alpha_cons + beta_clean Pi))` | claim 3 (O(Pi) dressing) | yes |
| A6 | sympy | 224 | `expect_zero(extracted beta_5 - beta_5)` | claim 3/4 (omega^5 coeff) | yes |
| A7 | sympy | 282 | `expect_zero(kappa_sel_template - kappa_sel_direct)` | claim 4 (selected overlap) | yes |
| A8 | sympy | 288 | `expect_zero(kappa_sel(0) - kappa0^2)` | claim 4 limit | yes |
| A9 | sympy | 289 | `expect_zero(limit(kappa_sel, oo) - sigma)` | claim 4 limit | yes |
| A10 | sympy | 300 | `expect_zero(delta_D_script - delta_D_paper)` | claim 4 (selected odd) | yes |
| B1 | math | 116 | `expectMatrixZero(sigmaWallSeq - sigmaExpected)` | claim 1 (independent elimination) | yes |
| B2 | math | 117-119 | `expectZero(sigma/xi/eta - literals)` | overlap constants | yes |
| B3 | math | 132 | `expectZero((k1t-k0t)-DeltaKax)` | none (typo guard) | no (acknowledged) |
| B4 | math | 148-149 | stationarity + stationarity@theta_- | angle law | yes |
| B5 | math | 166,170 | `beta - beta_clean`, `O(portPi)` | claim 3 | yes |
| B6 | math | 179,180 | `beta_5 - target`, extracted coeff | claim 3/4 | yes |
| B7 | math | 214,216,217 | kappa_sel cross-check + limits | claim 4 | yes |
| B8 | math | 229-232 | `deltaDScript - deltaDPaper` | claim 4 | yes |

A2/B3 are the only "no" rows; both are explicitly labeled in-comment as algebraically-forced typo guards (sympy 151-153, math 130-131), not as physics checks, and they coexist with substantive checks for the same quantities. They are not findings.

## Findings

None. See verdict justification.

## Independent-derivation check (Mathematica)

The Mathematica script is NOT a transliteration of the SymPy script; it derives the Schur complement by a genuinely different route. SymPy forms the full 4x4 internal mass matrix and inverts it in one step:
- sympy 100-119: `Mint = [[AU,0,-lambda_R k0,0],[0,AU,-lambda_R k1,0],[-lambda_R k0,-lambda_R k1,AW,0],[0,0,0,Aphi]]`; `Sigma = C.T * Mint.inv() * C`.

Mathematica instead eliminates the internal fields sequentially and sums three separately-derived rank pieces:
- math 90: `sigmaPhi = (lambdaB^2/aphi) Outer[Times,v,v]` (integrate out phi);
- math 97-101: `uMassInv = Inverse[diag(aU,aU)] + (lambdaR^2/(aU (aU aW - lambdaR^2 sigma))) Outer[Times,v,v]` (Woodbury-style W elimination);
- math 106-112: `sigmaU = cU^T uMassInv cU`, `sigmaW = ((aU lambdaW^2 + 2 lambdaR lambdaU lambdaW)/(aU aW - lambdaR^2 sigma)) Outer[Times,v,v]`, then `sigmaWallSeq = sigmaU + sigmaW + sigmaPhi`.

This is a structurally independent derivation (sequential block elimination + Woodbury identity) reaching the same boxed `Xi I_2 + alpha vv^T`, exactly what the second-engine policy requires. Likewise the `kappa_sel^2` check is computed two independent ways within each engine (Hellmann-Feynman `-d lambda_-/d alpha` vs. direct eigenvector nullspace/Eigensystem projection), and the engines use different selection mechanics (SymPy `nullspace`, Mathematica `Eigensystem` + position-match). No `mathematica_transliteration` finding.

## Engine cross-check

Both engines agree on every emitted form. Spot comparisons of the saved outputs:
- `Xi_0`: sympy `lambda_U**2/Omega_U**2` (txt 356) vs math `lambdaU^2/OmegaU^2` (txt 17) — agree.
- `Delta_0`: sympy `Omega_U**2*Omega_W**2 - 88*lambda_R**2/(9*pi**2)` (txt 357) vs math `OmegaU^2*OmegaW^2 - (88*lambdaR^2)/(9*Pi^2)` (txt 18) — agree.
- `beta_5`: sympy `81*pi**4*Gamma_port*(...)^2/(9*pi**2*Omega_U**2*Omega_W**2-88*lambda_R**2)**2` (txt 819) vs math `GammaPort*(...)^2/(OmegaU^2*OmegaW^2-(88*lambdaR^2)/(9*Pi^2))^2` (txt 34) — these are the same value (multiply num+den of the math form by `81 pi^4` to reach the sympy form). Both also pass `extracted beta_5 - expected beta_5 = 0` and the sympy `delta_D_script - delta_D_paper = 0` / math `deltaDScript - deltaDPaper = 0`.
- `sigma/xi/eta`: sympy `88/(9 pi^2)`, `56/(9 pi^2)`, `-8 sqrt(2)/(3 pi^2)` (txt 345-347) vs math literal checks all PASS (txt 11-16).

`engines_agree: true`.

## Verdict justification

`clean`. I read the paper card, the notes, and the appendix row before opening the scripts, then attacked the assertions. Attacks tried and failed: (1) I checked whether the Schur check is tautological — it is not, because `Sigma` is built by matrix inversion of `M_int` (sympy) / sequential Woodbury elimination (math) while `Sigma_expected` is the independently-written boxed formula, so the residual genuinely tests the elimination. (2) I checked whether the two acknowledged typo-guards (A2/B3) substitute for real physics checks — they do not; they sit alongside substantive checks and are commented as guards. (3) I verified the Pi- and al-derivatives are over real dependencies (`alpha` depends on `Pi` through `A_W`; `lambda_minus_template` depends on `al`), so neither `diff` is identically zero and the `O(Pi)` / Hellmann-Feynman checks are non-trivial. (4) I checked the two limit assertions (al→0 → kappa0^2, al→∞ → sigma) are real boundary tests of the eigenvector overlap, and that the strong-loading limit (sympy `sp.limit(..., al, oo)`, math `Limit[..., al->Infinity, Assumptions->DeltaKax>0]`) is properly conditioned. (5) I confirmed both saved outputs ran to completion with all `expect_zero`/`expectZero` residuals = 0 and no AssertionError/FAIL, and that outputs are newer than scripts. The script's claim matches the paper's claim on all four deliverables plus the notes' angle-law and limiting identities. The only blemish is cosmetic: the SymPy docstring (lines 9, 16, 21, 134) calls this "Stage-11" loading, the pre-renumber name for Stage 028 (11+17=28, the known +17 renumber drift); this is a comment-only label, touches no assertion, and is not a `paper_misalignment` of any verified value — noted, not filed.

## Value Reconciliation (pass-2 augmentation)

reconciliation: complete; 12 deliverable values checked, 0 misaligned.

Every RESULT/deliverable value the scripts emit is accounted for below. Values were read from the script source and the committed saved outputs (sympy `.txt`, math `.txt`); nothing was executed. Outputs are fresh (sympy/math `.txt` mtime May 26 > script mtime May 25), so the saved transcripts reflect the current scripts.

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| `sigma = 88/(9 pi^2)` | py txt 345; wl txt 11 (`PASS sigma - 88/(9 Pi^2)`) | tex 35-36 `\sigma=...=\frac{88}{9\pi^2}`; md 50 | MATCH |
| `Xi(omega) = lambda_U^2/A_U` | py line 88; wl txt 6 | tex 50 (boxed); md 111 | MATCH |
| `alpha(omega)` rational form | py line 89; wl txt 7 | tex 52-53 (boxed); md 115-117 | MATCH |
| `Xi_0 = lambda_U^2/Omega_U^2` | py txt 356; wl txt 17 | tex 59 (boxed); md 145 | MATCH |
| `Delta_0 = Omega_U^2 Omega_W^2 - lambda_R^2 sigma` | py txt 357; wl txt 18 | tex 61 (boxed); md 147 | MATCH |
| `alpha_0` static load | py txt 358; wl txt 19 | tex 63-64 (boxed); md 149-150 | MATCH |
| `beta(omega) = (A_U lambda_W + lambda_R lambda_U)^2 / Delta_cons^2` | py line 202; wl txt 28 | tex 84-86 (boxed `beta`, `Delta_cons`); md 208-212 | MATCH |
| `beta_0 = (Omega_U^2 lambda_W + lambda_R lambda_U)^2 / Delta_0^2` | py line 216 (beta5 at Gamma=1); wl txt 34 | tex 92 (boxed `beta_0`); md 218 | MATCH |
| `beta_5 = Gamma_port * beta_0` | py txt 819; wl txt 34 | md 244-250 (`beta_5 = Gamma_2^port beta_0`); tex 96-103 via selected-odd | MATCH |
| `tan(2 theta_-) = 2 alpha_0 eta / (Delta K_ax + alpha_0 xi)` | py txt 515; wl txt 26 | md 167-168 | MATCH |
| `kappa_sel^2 = (v.e_-)^2` (HF = direct, limits kappa0^2 / sigma) | py 264-289; wl 203-217 | md 267 `(v.e_-)^2 = kappa(theta_-)^2`; tex 103,112 `(v\cdot e_-)^2` | MATCH |
| `delta D_-^odd = -i Gamma beta_0 (v.e_-)^2 omega^5` | py 299-303; wl txt 48-50 | tex 102-103 (boxed selected-odd); md 271-272 | MATCH |

Internal scaffolding (accounted for, no finding expected): `xi = 56/(9 pi^2)` and `eta = -8 sqrt(2)/(3 pi^2)` are intermediate combinations of `kappa_0=2 sqrt(2)/pi`, `kappa_1=-4/(3 pi)` (md 46-48; tex 35 `\kappa_0,\kappa_1`); they feed `tan(2 theta_-)` (where `kappa_0^2-kappa_1^2` and `kappa_0 kappa_1` appear symbolically in md 168) rather than being separately-boxed deliverables. Also internal-only: `D_bare` (bare operator, given as input eq. D-bare tex 18-19), `A_phi/A_U/A_W` kernels (definitions, tex 22-28 / md 73-78), `Delta_cons` (named in tex 86), `alpha_cons(omega)`, `Delta_UW(omega)` (tex 33), `lambda_-/lambda_+` eigenvalues, `kappa(theta)^2` symbolic form, `K_eff^(0)` matrix (tex 69-70), and all verification scaffolding (residual-zero flags, `DeltaK_tilde - DeltaK_ax` typo guard, series-truncation intermediates). All such items either appear correctly in the docs as definitions/inputs or are genuine internal scaffolding; none is a missing deliverable.

All 12 deliverable values reconcile against the `.tex` card and/or the `.md` notes. No MISMATCH, no MISSING-DELIVERABLE. Standard verdict `clean` is retained.

## Self-test notes

I checked the variable-independence trap: `sp.diff(alpha, Pi)` (sympy 200) is nonzero because `alpha` depends on `Pi` via `A_W`→`Delta_UW`; `sp.diff(lam_minus_template, al)` (sympy 264) is nonzero because the template depends on `al` — so neither the `O(Pi)` expansion nor the Hellmann-Feynman `kappa_sel^2` collapses to a trivial 0=0. I checked the trivial-case/limit pre-checks: al→0 gives `kappa0^2` and al→∞ gives `sigma`, both confirmed as residual=0 in both engines, which are real boundary exercises of the eigenvector overlap. I confirmed no integrals over symmetric domains are involved (this stage is algebraic/eigenvalue, no parity trap). Conclusion: assertions are substantive and the two acknowledged typo-guards do not substitute for physics checks.
