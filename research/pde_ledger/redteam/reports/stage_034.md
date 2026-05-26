---
unit_id: 034
batch: II.1
auditor_model: claude-opus-4-7[1m]
audit_date: 2026-05-25T00:00:00Z
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
    - /var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage034_softening_depth_normal_form.md
  paper_appendix: present
---

# Audit unit 034 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_034.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage034_softening_depth_normal_form.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part02.tex` (only the `\input{stages/stage_034}` line at L106 references this unit — there is no separate appendix row beyond inclusion)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage034_softening_depth_normal_form_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage034_softening_depth_normal_form_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage034_softening_depth_normal_form_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage034_softening_depth_normal_form_mathematica_audit.txt`

## What the paper claims

Quoting `\stagefield{Output}` (stage_034.tex:93): *"Stage~034 outputs the softening-depth variable [x:=A-lambda_-], the loading law [alpha_0(x) = x(x+DeltaK_ax) / (kappa_0^2(x+DeltaK_ax) + kappa_1^2 x)], the selected overlap [s_-(x)], the normalization product [N_-(x)], and the required support formula [g_{B,req}^2/varpi^2]."*

Five distinct deliverables, all derived from the rank-one two-level secular equation `1 = alpha_0[kappa_0^2/x + kappa_1^2/(x + DeltaK_ax)]` (eq. app-stage034-secular):

1. **Softening-depth variable** (eq:app-stage034-x): `x := A - lambda_-`, stable branch `0 <= x < A`.
2. **Loading law** (eq:app-stage034-alpha-x): `alpha_0(x) = x(x+DeltaK_ax) / [kappa_0^2(x+DeltaK_ax) + kappa_1^2 x]`.
3. **Selected overlap** (eq:app-stage034-s-x): `s_-(x) = [kappa_0^2(x+DeltaK_ax) + kappa_1^2 x]^2 / [kappa_0^2(x+DeltaK_ax)^2 + kappa_1^2 x^2]`.
4. **Normalization product** (eq:app-stage034-N-x): `N_-(x) = beta_0 [kappa_0^2(x+DeltaK_ax) + kappa_1^2 x]^4 / { kappa_0^2(A-x) [kappa_0^2(x+DeltaK_ax)^2 + kappa_1^2 x^2]^2 }`.
5. **Required support** (eq:app-stage034-support-req): `g_{B,req}^2/varpi^2 = x(x+DeltaK_ax) / [kappa_0^2(x+DeltaK_ax) + kappa_1^2 x] - Chi^2 / (Omega_U^2 Delta_0)`.

The paper also asserts (eq:app-stage034-dalpha-dx) the manifestly positive form `d alpha_0/dx = [kappa_0^2(x+DeltaK_ax)^2 + kappa_1^2 x^2] / [kappa_0^2(x+DeltaK_ax) + kappa_1^2 x]^2 > 0`, and the checklist item that `s_- = dx/d alpha_0` follows from `x = A - lambda_-` and `s_- = -d lambda_-/d alpha_0`. The notes file at lines 27-31 mirrors the same five-item list and additionally states monotonicity as deliverable (4) in its enumeration. Inputs are explicitly Stage 032's rank-one secular equation and selected product `N_-`.

## What the script claims to verify

Both engines build the lambda-form expressions `alpha(lambda) = 1/(kappa0_sq/(A-lambda) + kappa1_sq/(A+DeltaK-lambda))`, `s_-(lambda) = S1^2/(dS1/dlambda)` with `S1 = kappa0_sq/(A-lambda) + kappa1_sq/(A+DeltaK-lambda)`, and `N_-(lambda) = beta_0 s_-(lambda)^2 / (kappa0_sq * lambda)`. They then independently write the paper-stated x-form closed forms `alpha_x`, `s_x`, `N_x` (sympy lines 48-50; Mathematica 46-58) and assert that substituting `lambda -> A - x` into each lambda-form reproduces the corresponding x-form symbolically (residual = 0 after simplify). They verify `d alpha_x/dx` matches the manifestly positive form (sympy:70; WL:75), check the reciprocity `s_x * d alpha_x/dx == 1` (sympy:71; WL:76) — which encodes `s_- = dx/d alpha_0`, and check that solving the support-partition equation `gB^2/varpi^2 + alpha_mix == alpha_x` in the x-variable agrees, under `lambda -> A - x`, with the same solve done in the lambda-variable (sympy:84-92; WL:85-97).

## Paper <-> script cross-check

| Paper deliverable | Script-side check(s) | Status |
|---|---|---|
| (1) softening-depth definition `x := A - lambda_-` | Used as substitution rule `lambda -> A - x` in all three substitution checks (sympy:61-63; WL:65-67) | match |
| (2) loading law `alpha_0(x)` | `alpha_lam.subs(lam, A-x) - alpha_x == 0` (sympy:61; WL:65) | match |
| (3) selected overlap `s_-(x)` | `s_lam.subs(lam, A-x) - s_x == 0` (sympy:62; WL:66) | match |
| (4) normalization product `N_-(x)` | `N_lam.subs(lam, A-x) - N_x == 0` (sympy:63; WL:67) | match |
| (5) required support `g_{B,req}^2/varpi^2` | x-form constructed by solve; lambda-form solved independently and verified to agree under `lambda -> A - x` (sympy:78,86-91; WL:85-87,91-96) | match |
| (body) `d alpha_0/dx` manifestly positive form | `dalpha_dx - dalpha_target == 0` (sympy:70; WL:75) | match |
| (checklist) `s_- = dx/d alpha_0` | `s_x * dalpha_dx - 1 == 0` (sympy:71; WL:76) | match |

`paper_alignment` set to `aligned` — every paper-listed deliverable has at least one substantive script-side check; no script-side check is orphan to a missing paper claim.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 61 | `simplify(alpha_lam.subs(lam, A-x) - alpha_x) == 0` | claim 2 (loading law) | yes |
| A2 | sympy | 62 | `simplify(s_lam.subs(lam, A-x) - s_x) == 0` | claim 3 (selected overlap) | yes |
| A3 | sympy | 63 | `simplify(N_lam.subs(lam, A-x) - N_x) == 0` | claim 4 (normalization product) | yes |
| A4 | sympy | 70 | `simplify(dalpha_dx - dalpha_target) == 0` | body eq dalpha/dx (manifestly positive) | yes |
| A5 | sympy | 71 | `simplify(s_x * dalpha_dx - 1) == 0` | checklist `s_- = dx/d alpha_0` | yes |
| A6 | sympy | 89-92 | `simplify(gBreq_lambda.subs(lam, A-x) - gBreq_sq_over_varpi2) == 0` | claim 5 (required support) | yes |
| A7 | mathematica | 65 | `expectZero[(alphaLambda /. lambda -> A-x) - alphaX]` | claim 2 | yes |
| A8 | mathematica | 66 | `expectZero[(sLambda /. lambda -> A-x) - sX]` | claim 3 | yes |
| A9 | mathematica | 67 | `expectZero[(nLambda /. lambda -> A-x) - nX]` | claim 4 | yes |
| A10 | mathematica | 75 | `expectZero[dAlphaDx - dAlphaTarget]` | body eq dalpha/dx | yes |
| A11 | mathematica | 76 | `expectZero[sX*dAlphaDx - 1]` | checklist `s_- = dx/d alpha_0` | yes |
| A12 | mathematica | 94-97 | `expectZero[(gBReqLambda /. lambda -> A-x) - gBReqSqOverVarpi2]` | claim 5 | yes |

Every assertion is non-tautological: each constructs the lambda-side expression from the secular form, the x-side expression from the paper-stated closed form, and verifies their equivalence symbolically. The CAS simplifier must do real algebraic work; if either the closed form or the differentiation/solve step were wrong, the residual would be a nonzero rational expression.

## Findings

None. Verdict is `clean`.

(The prior v1 audit at this same path flagged a `tautological_check` on the support-loading rearrangement; that finding was already resolved — the current script at sympy:86-92 and WL:91-97 implements the v1 directive's recommended lambda-form-vs-x-form cross-check rather than the original by-construction residual.)

## Independent-derivation check (Mathematica)

The Mathematica file is structurally parallel to the SymPy file. Variable names mirror across engines (`kappa0Sq`/`ks0`, `alphaLambda`/`alpha_lam`, `s1`/`S1`, `s1p`/`S1p`, `sLambda`/`s_lam`, `nLambda`/`N_lam`, `alphaX`/`alpha_x`, `sX`/`s_x`, `nX`/`N_x`, `dAlphaDx`/`dalpha_dx`, `dAlphaTarget`/`dalpha_target`, `gBReqLambda`/`gBreq_lambda`). The same algebraic choreography is followed in the same order: build `alphaLambda` and `s1`/`s1p`, derive `sLambda = s1^2/s1p`, build `nLambda`, hardcode `alphaX`/`sX`/`nX` from the paper, run the three substitution checks, run the derivative check, run the reciprocity check, solve for `gB` in both forms.

Three corresponding sections to justify the call:

- SymPy lines 42-45 vs WL lines 38-44:
  ```
  S1 = sp.simplify(ks0 / (A - lam) + ks1 / (A + DeltaK - lam))
  S1p = sp.simplify(sp.diff(S1, lam))
  s_lam = sp.simplify(S1**2 / S1p)
  N_lam = sp.simplify(beta0 * s_lam**2 / (ks0 * lam))
  ```
  vs
  ```
  s1 = FullSimplify[kappa0Sq/(A - lambda) + kappa1Sq/(A + DeltaK - lambda), ...];
  s1p = FullSimplify[D[s1, lambda], ...];
  sLambda = FullSimplify[s1^2/s1p, ...];
  nLambda = FullSimplify[beta0*sLambda^2/(kappa0Sq*lambda), ...];
  ```
- SymPy lines 48-50 vs WL lines 46-58: same x-form closed forms, hardcoded identically with engine-syntactic differences only.
- SymPy lines 61-63 vs WL lines 65-67: identical three-line block of substitution residual checks, same order.

That said, the load-bearing computation in each assertion is the engine's own simplifier (`sp.simplify(sp.expand(...))` in the SymPy `expect_zero` helper at line 25; `FullSimplify[Together[Expand[...]]]` in the WL `expectZero` at line 21) operating on a residual rational expression. Each engine reduces the residual to 0 using its native rational-function machinery; a bug in either simplifier would surface here independently. Additionally, the verification target IS an algebraic-equivalence claim — the paper deliverable is exactly "these two parametrizations agree under the substitution `lambda -> A - x`", so writing the substitution explicitly in both engines is not "doing the same algebra"; it is doing the same TEST. There is no obviously richer independent reformulation available without restating the test itself (e.g., substituting `alpha_x` into the secular `1 = alpha_x*(ks0/x + ks1/(x+DeltaK))` would be an equivalent but no-stronger check). On balance I do not flag `mathematica_transliteration`; the borderline parallelism is appropriate given the nature of the claim.

## Engine cross-check

Both engines compute identical printed closed forms:

| Quantity | SymPy printout | Mathematica printout |
|---|---|---|
| `alpha(x)` | `x*(DeltaK + x) / (kappa_0_sq*(DeltaK+x) + kappa_1_sq*x)` (txt L9-12) | `(x*(DeltaK + x))/(DeltaK*kappa0Sq + (kappa0Sq + kappa1Sq)*x)` (txt L6) — same after distribution |
| `s_-(x)` | `(kappa_0_sq*(DeltaK+x)+kappa_1_sq*x)^2 / (kappa_0_sq*(DeltaK+x)^2 + kappa_1_sq*x^2)` (txt L13-18) | `(DeltaK*kappa0Sq + (kappa0Sq+kappa1Sq)*x)^2 / (kappa1Sq*x^2 + kappa0Sq*(DeltaK+x)^2)` (txt L7) — same |
| `N_-(x)` | `beta_0 * (...)^4 / [kappa_0_sq*(A-x)*(...)^2]` (txt L19-25) | `beta0*(...)^4 / (kappa0Sq*(A-x)*(...)^2)` (txt L8) — same |
| `g_{B,req}^2/varpi^2` | `(-Chi^2*(DeltaK*kappa0Sq + (kappa0Sq+kappa1Sq)*x) + Delta0*OmegaU^2*x*(DeltaK+x)) / (Delta0*OmegaU^2*(...))` (txt L50-55) | `-(Chi^2/(Delta0*OmegaU^2)) + x*(DeltaK+x)/(...)` (txt L21) — algebraically equal |

All six residuals in each script reduce to 0; both scripts print `All Stage 17 checks passed.` / `Stage 034 Mathematica audit passed.` No disagreement.

Output freshness: sympy script mtime `May 21 17:32`, output mtime `May 21 17:33` (fresh). Mathematica script mtime `May 21 17:32`, output mtime `May 21 17:33` (fresh). Both outputs reflect the current scripts.

## Verdict justification

The unit is `clean`. I attacked the audit in five ways and all attacks failed:

1. **Tried to break the closed forms by hand.** Substituting `lambda = A - x` into `alpha_lam = 1/(ks0/(A-lam) + ks1/(A+DeltaK-lam))` gives `1/(ks0/x + ks1/(x+DeltaK)) = x(x+DeltaK)/(ks0(x+DeltaK) + ks1*x)`, which is exactly `alpha_x` (sympy:48). Same for `s_lam = S1^2/S1'` -> `s_x` and `N_lam = beta_0 s_lam^2 / (ks0*lam)` -> `N_x`. All three substitutions reduce to the paper-stated closed forms exactly. No factor errors, no sign flips, no missing terms.

2. **Tried to find a missing paper deliverable.** All five `\stagefield{Output}` items plus the body `d alpha/dx` form plus the checklist `s_- = dx/d alpha_0` identity have script-side counterparts. No `script_missing_paper_claim` finding.

3. **Tried to find an extra unverified script claim.** Every assertion traces to a paper deliverable or paper-listed checklist item. No `paper_missing_script_claim`.

4. **Tried to find domain-assumption errors.** SymPy declares `lam` `positive=True` but does not impose `lam < A`; `x` is `nonnegative=True` without imposing `x < A`. The substitution checks don't require the upper bound (they're purely symbolic identity checks). Mathematica is stricter, imposing `0 < lambda < A` and `0 <= x < A` in `$Assumptions`. No `simplify` step requires the upper bound to be valid. No `symbol_assumption_error`.

5. **Tried to find a load-bearing tautology.** The lambda-form / x-form solve cross-check (sympy:86-91, WL:91-96) is genuinely non-tautological: `gBreq_lambda` is built by SymPy/Mathematica solving `gB + alpha_mix == alpha_lam` for `gB`, yielding `alpha_lam - alpha_mix`. The residual `(alpha_lam - alpha_mix).subs(lam, A-x) - (alpha_x - alpha_mix)` simplifies to `alpha_lam.subs(lam, A-x) - alpha_x`, which is the load-bearing identity from check A1. So A6/A12 effectively re-exercise A1/A7 with the additional solve-step in the path; not stronger than A1 but not tautological either.

The previously-flagged tautological_check (v1 audit) on the support-loading rearrangement is resolved: lines 84-91 of sympy and 91-97 of WL now implement the lambda-form vs x-form cross-check rather than the original `gBreq - (alpha_x - alpha_mix) == 0` residual.

I confirm I read the paper card (stage_034.tex), the notes file (moving_throat_pde_stage034_softening_depth_normal_form.md), and the appendix inclusion line in stage_appendix_part02.tex before opening the scripts. The script's verified claims match the paper's stated claims item-for-item.

## Self-test notes

- **Variable independence**: The two derivative calls (`sp.diff(S1, lam)` at sympy:43 and `sp.diff(alpha_x, x)` at sympy:66; corresponding `D[s1, lambda]` and `D[alphaX, x]` in WL) target variables that genuinely appear in the differentiand — `S1` is a function of `lam` and `alpha_x` is a function of `x`. No identically-zero derivatives.
- **Symmetry/parity**: N/A — no integrals here, only rational-function identities.
- **Trivial-case pre-check**: For each substitution residual, at `ks0=1, ks1=1, A=2, DeltaK=1, x=1` (so `lam=1`): `alpha_lam = 1/(1/1 + 1/2) = 2/3`; `alpha_x = 1*2/(1*2 + 1*1) = 2/3`. Match. `s_lam = S1^2/S1'` with `S1=3/2`, `dS1/dlam = 1/(2-lam)^2 + 1/(3-lam)^2 = 1 + 1/4 = 5/4`, so `s_lam = (9/4)/(5/4) = 9/5`; `s_x = (1*2+1*1)^2/(1*4 + 1*1) = 9/5`. Match. The residuals genuinely reduce to 0 at concrete values, confirming the symbolic check is not vacuously passing.
- **Paper round-trip**: Each script-side check name and form maps onto a paper equation label (eq:app-stage034-alpha-x for A1/A7, eq:app-stage034-s-x for A2/A8, eq:app-stage034-N-x for A3/A9, eq:app-stage034-dalpha-dx for A4/A10, the body sentence "s_- = dx/d alpha_0" in the checklist for A5/A11, eq:app-stage034-support-req for A6/A12). No new paper-misalignment introduced.
