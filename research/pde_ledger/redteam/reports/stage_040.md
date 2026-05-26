---
unit_id: 040
batch: III.1
auditor_model: claude-opus-4-7[1m]
audit_date: 2026-05-26T00:00:00Z
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
  notes_stage_files: ["/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage040_generalized_selected_branch.md"]
  paper_appendix: present
---

# Audit unit 040 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_040.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage040_generalized_selected_branch.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part03.tex` (rows for stage 040 only)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage040_generalized_selected_branch_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage040_generalized_selected_branch_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage040_generalized_selected_branch_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage040_generalized_selected_branch_mathematica_audit.txt`

## What the paper claims

Stage 040 replaces the flat one-vector selected branch by a source/loading mismatch law. Verbatim `\stagefield{Output}`: "The generalized selected functions \eqref{eq:app-stage040-Fqeta}--\eqref{eq:app-stage040-Gq} and their split--$U$ specialization \eqref{eq:app-stage040-FU}--\eqref{eq:app-stage040-GU}." Concretely, the stage proves:
1. The required rank-1 loading is `alpha_req = A0*xi*(xi+delta)/[z0^2*(delta+(1+q^2)*xi)]` for `0 <= xi < 1` (eq. app-stage040-alpha-req).
2. The selected eigenvector has ratio `e1/e0 = q*xi/(delta+xi)` (eq. app-stage040-evec).
3. The generalized normalization function `F_{q,eta_s}(xi,delta) = [delta+(1+q^2)xi]^2 [delta+(1+eta_s)xi]^2 / ((1-xi)[(delta+xi)^2+q^2*xi^2]^2)` with `eta_s = s1*z1/(s0*z0)` (eq. app-stage040-Fqeta).
4. The deformed loading function `G_q(xi,delta) = xi*(xi+delta)/[delta+(1+q^2)xi]` (eq. app-stage040-Gq).
5. The split-U specialization `q = -sqrt(2)/3 * R_U`, `eta_s = (2/9) R_U` (eq. app-stage040-qeta-RU), yielding `F_U(xi,delta;R_U)` and `G_U(xi,delta;R_U)` (eqs. app-stage040-FU and app-stage040-GU).
6. The flat branch is recovered at `R_U = 1` (matches Stage-18/19 functions). The notes additionally enumerate (Section 5) the first-order deformation about the flat-U limit with explicit `H_F`, `H_G`.

## What the script claims to verify

The SymPy script (and its Mathematica counterpart) verify exactly the items above:
1. The exact `alpha_req` closed form (sympy L59) and the eigenvector residual `(M - alpha_req z z^T - lam I) e_- = 0` for `e_- = (1, q*xi/(delta+xi))^T` (sympy L73-79).
2. Construction of `F_general` from normalized overlaps and `A0/lam_minus`, with equality to paper's literal `F_expected` (sympy L89, L95-98).
3. Construction of `G_general = (z0^2/A0)*alpha_req` with equality to paper's `G_expected` (sympy L90, L96, L99).
4. Substitution `q -> -sqrt(2/9)*R_U`, `eta -> (2/9)*R_U` to form `F_U`, `G_U`, with verification that `F_U(R_U=1) - F_stage18 = 0` and `G_U(R_U=1) - G_stage19 = 0` (sympy L105-122).
5. First-order deformation `H_F`, `H_G` derived by two independent routes (differentiating `F_U` after R_U substitution vs. differentiating `F_general` after eps-parametrized (q,eta) substitution) and cross-checked equal (sympy L124-149).
The Mathematica script independently derives `alpha_req` via `Solve[Det[mPert[alpha] - lam I] == 0, alpha]` and the eigenvector via `NullSpace[...]`, not by copying SymPy's algebraic substitution path.

## Paper - script cross-check

| Paper deliverable | Script-side check | Status |
|---|---|---|
| `alpha_req` (eq. alpha-req) | sympy L59 closed form + L73-79 residual; math L52-55 Solve | match |
| `e1/e0 = q*xi/(delta+xi)` (eq. evec) | sympy L64-66, L77-79 residual; math L58-65 NullSpace, L74-77 residual | match |
| `F_{q,eta_s}` (eq. Fqeta) | sympy L89, L95-98; math L90, L99-105 | match |
| `G_q` (eq. Gq) | sympy L90, L96, L99; math L91, L104, L106 | match |
| Split-U `q, eta_s` (eq. qeta-RU) | sympy L105-106; math L110-111 | match |
| `F_U`, `G_U` (eqs. FU, GU) via R_U=1 recovery | sympy L107-122; math L112-129 | match |
| Notes Sec. 5 `H_F`, `H_G` first-order deformation | sympy L128-149 two-path; math L133-150 two-path | match (extra detail beyond .tex Output but anchored in notes) |

`paper_alignment`: `aligned`.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 66 | `expect_zero("e1/e0 closed form", r - xi*q/(delta+xi))` | claim 2 (evec ratio) | yes |
| A2 | sympy | 78 | `expect_zero("eigenvector residual row 0", eig_residual[0])` | claims 1+2 (alpha_req and evec consistent with lam_minus) | yes |
| A3 | sympy | 79 | `expect_zero("eigenvector residual row 1", eig_residual[1])` | claims 1+2 | yes |
| A4 | sympy | 98 | `expect_zero("F_general - expected", F_general - F_expected)` | claim 3 (F_{q,eta_s} closed form) | yes |
| A5 | sympy | 99 | `expect_zero("G_general - expected", G_general - G_expected)` | claim 4 (G_q closed form) | yes |
| A6 | sympy | 121 | `expect_zero("F_U(R_U=1) - Stage18 F", ...)` | claim 6 (flat-U recovery of F) | yes |
| A7 | sympy | 122 | `expect_zero("G_U(R_U=1) - Stage19 G", ...)` | claim 6 (flat-U recovery of G) | yes |
| A8 | sympy | 148 | `expect_zero("H_F cross-check (F_U vs F_general)", HF - HF_direct)` | notes Sec. 5 H_F | yes |
| A9 | sympy | 149 | `expect_zero("H_G cross-check (G_U vs G_general)", HG - HG_direct)` | notes Sec. 5 H_G | yes |
| B1 | math | 65 | `expectZero["e1/e0 closed form", r - xi q/(delta + xi)]` | claim 2 | yes |
| B2 | math | 76 | `expectZero["eigenvector residual row 0", eigResidual[[1]]]` | claims 1+2 | yes |
| B3 | math | 77 | `expectZero["eigenvector residual row 1", eigResidual[[2]]]` | claims 1+2 | yes |
| B4 | math | 105 | `expectZero["F_general - expected", fGeneral - fExpected]` | claim 3 | yes |
| B5 | math | 106 | `expectZero["G_general - expected", gGeneral - gExpected]` | claim 4 | yes |
| B6 | math | 128 | `expectZero["F_U(R_U=1) - Stage18 F", (fU /. rU -> 1) - fStage18]` | claim 6 | yes |
| B7 | math | 129 | `expectZero["G_U(R_U=1) - Stage19 G", (gU /. rU -> 1) - gStage19]` | claim 6 | yes |
| B8 | math | 149 | `expectZero["H_F cross-check (F_U vs F_general)", hF - hFDirect]` | notes Sec. 5 H_F | yes |
| B9 | math | 150 | `expectZero["H_G cross-check (G_U vs G_general)", hG - hGDirect]` | notes Sec. 5 H_G | yes |

All assertions are non-tautological and trace to a specific paper-side or notes-side deliverable.

## Findings

None.

Attacks attempted and dismissed:

- **alpha_req mismatch attack**: Paper's denominator is `z0^2*[delta + (1+q^2)*xi]`. SymPy L59 writes `(z0**2*(delta + xi) + (q*z0)**2*xi)` which expands to `z0^2*(delta+xi+q^2*xi) = z0^2*(delta+(1+q^2)xi)`. Match. Mathematica independently solves `Det[mPert[alpha] - lamMinus I] == 0` for alpha and gets the same expression. No finding.
- **eta vs eta_s naming**: Paper uses `eta_s`, script uses `eta`. Both defined as `s1*z1/(s0*z0)`, and in script comment L47 `eta := (s1/s0)*q`, which is identical since `q=z1/z0`. Cosmetic rename, not misalignment.
- **xi domain attack**: Paper says `0 <= xi < 1`. SymPy declares `xi` as `positive=True` (excludes `xi=0`); does not encode `xi < 1`. The factor `(1-xi)` appears symbolically; no division is performed during the algebraic check. The boundary `xi=0` is a trivial endpoint not used by any of the stage's theorem statements. Not a finding.
- **R_U positivity / eps positivity**: SymPy declares both as `positive=True`. Paper notes (Sec. 5) allow `R_U = 1 + eps` with eps potentially negative on the natural branch (`rho_0 > 0` -> `R_U < 1`). However, the script's checks are symbolic identities. The derivative-at-eps=0 check does not depend on the sign of eps, and the substitution `R_U -> 1` for Stage-18/19 recovery is also sign-independent. Not a finding.
- **Hardcoded Stage-18/19 forms (sympy L118-119; math L120-124)**: These are literals, but their comments cite upstream verification scripts (`stage035_*_sympy_audit.py` lines 46-58 and `stage036_*_sympy_audit.py` lines 53-70). They are used as targets in `F_U(R_U=1) - F_stage18 = 0`, a non-tautological reduction-to-limit check. This is the standard carry-forward pattern with explicit provenance; not `hardcoded_result`.
- **Mathematica transliteration check**: SymPy directly substitutes the closed-form `alpha_req` and computes `r` from a sign-fixed algebraic expression. Mathematica solves `Solve[Det[...] == 0, alpha]` and extracts the eigenvector via `NullSpace[...]`. These are genuinely different derivation paths converging to the same physics. Not transliteration.
- **Tautology attack on H_F cross-check**: `HF` is `diff(F_U.subs(R_U, 1+eps), eps).subs(eps, 0) / F_stage18`. `HF_direct` is `diff(F_general.subs({q: q_U_eps, eta: eta_U_eps}), eps).subs(eps, 0) / F_stage18`. The two paths use different intermediate expressions (`F_U` already has q,eta substituted; `F_general` retains the symbolic q,eta and substitutes the eps-parametrized versions at differentiation time). The assertion `HF - HF_direct == 0` is a real cross-check of two algebra paths, not a CAS tautology.
- **Tautology attack on F_general assertion**: `F_general` is `(A0/lam_minus)*z_overlap_sq*s_overlap_sq` constructed from analytically derived overlaps, with `r = (A0*xi - alpha_req*z0^2)/(alpha_req*z0*(q*z0))` (sympy L64) — a route through the eigenvalue equation. `F_expected` is the paper's literal closed form. They are connected only through nontrivial algebraic identities. Equality is a real check.
- **Engine disagreement attack**: Both engines produce identical `alpha_req`, identical `e1/e0`, identical overlaps and F/G, identical F_U/G_U, identical H_F/H_G. All `expectZero`/`expect_zero` checks PASS in both transcripts. No disagreement.
- **Stage labeling drift**: Script docstring and banner say "Stage 23"; file/paper say "Stage 040". Notes file uses "Stage 040" in heading but refers internally to "Stage 22/23". This is a known internal-vs-paper numbering convention, not a math finding.

## Independent-derivation check (Mathematica)

The `.wl` is **not** a transliteration. Key independent moves:

- **Section 1, alpha derivation (lines 52-55)**:
  ```
  charEq = Det[mPert[alpha] - lamMinus IdentityMatrix[2]] == 0;
  alphaSol = Solve[charEq, alpha];
  alphaReq = FullSimplify[alpha /. alphaSol[[1]], ...];
  ```
  This actually solves the characteristic polynomial for `alpha` given `lam = a0(1-xi)`. The SymPy script's L59 writes the closed form directly. Different derivation paths to the same target.

- **Section 1, eigenvector via NullSpace (lines 57-65)**:
  ```
  nsVec = NullSpace[mPert[alphaReq] - lamMinus IdentityMatrix[2]];
  eMinusRaw = FullSimplify[nsVec[[1]], ...];
  eMinus = FullSimplify[eMinusRaw/eMinusRaw[[1]], ...];
  ```
  Computes the null space of the loaded matrix at the selected eigenvalue and normalizes. SymPy constructs `r` algebraically from `(A0*xi - alpha_req*z0^2)/(alpha_req*z0*(q*z0))`. Different derivations.

- **Section 2, overlaps (lines 83-91)**: Mathematica uses `sVec = {1, eta/q}` and `zVec = {z0, q z0}`, computes `(z.e)^2 / (z0^2 ||e||^2)` and `(s.e)^2 / (1 * ||e||^2)`. SymPy assembles them via the analytic identities `(1+qr)^2/(1+r^2)` and `(1+eta*xi/(delta+xi))^2/(1+r^2)`. Mathematica's approach starts from vectors directly; SymPy's from rationalized algebra. Different and convergent.

No `mathematica_transliteration` finding.

## Engine cross-check

| Quantity | SymPy output | Mathematica output | Agree? |
|---|---|---|---|
| `alpha_req` | `A0*xi*(delta + xi)/(z0**2*(delta + q**2*xi + xi))` | `(a0*xi*(delta + xi))/((delta + xi + q^2*xi)*z0^2)` | yes |
| `e1/e0` | `q*xi/(delta + xi)` | `(q*xi)/(delta + xi)` | yes |
| `(z.e_-)^2/z0^2` (normalized) | `(delta+q**2*xi+xi)**2/(delta**2+2*delta*xi+q**2*xi**2+xi**2)` | `(delta+xi+q^2*xi)^2/(delta^2+2*delta*xi+(1+q^2)*xi^2)` | yes |
| `(s.e_-)^2/s0^2` (normalized) | `(delta+eta*xi+xi)**2/(delta**2+2*delta*xi+q**2*xi**2+xi**2)` | `(delta+xi+eta*xi)^2/(delta^2+2*delta*xi+(1+q^2)*xi^2)` | yes |
| `F_(q,eta)` | matches `F_expected` (residual 0) | matches `fExpected` (residual 0) | yes |
| `G_q` | `xi*(delta+xi)/(delta+q**2*xi+xi)` | `(xi*(delta+xi))/(delta+xi+q^2*xi)` | yes |
| `F_U(R_U=1) - F_stage18` | 0 | 0 | yes |
| `G_U(R_U=1) - G_stage19` | 0 | 0 | yes |
| `H_F` | `4*xi*(27*delta**2+36*delta*xi+11*xi**2)/((9*delta+11*xi)*(9*delta**2+18*delta*xi+11*xi**2))` | `(4*xi*(27*delta^2+36*delta*xi+11*xi^2))/((9*delta+11*xi)*(9*delta^2+18*delta*xi+11*xi^2))` | yes |
| `H_G` | `-4*xi/(9*delta+11*xi)` | `(-4*xi)/(9*delta+11*xi)` | yes |

All `expectZero` / `expect_zero` checks return 0 in both transcripts. Engines fully agree.

Output mtimes both newer than corresponding script mtimes (SymPy: script 1779474522 < output 1779474740; Mathematica: script 1779474673 < output 1779474749). Outputs fresh.

## Verdict justification

`clean`. The paper card states six exact deliverables for Stage 040; both scripts verify all six, plus the notes-only first-order deformation in Section 5, by non-tautological algebraic checks that derive the targets via genuinely different routes (closed-form substitution and explicit matrix residual in SymPy vs. `Solve[Det == 0]` and `NullSpace` in Mathematica). The engines produce identical final expressions for every quantity and identical PASS outcomes for every `expectZero`. I attempted attacks on: `alpha_req` algebraic form (matches paper exactly after expansion), `eta` vs `eta_s` naming (identical definition), `xi`/`R_U`/`eps` positivity declarations (do not affect symbolic identities tested), hardcoded Stage-18/19 forms (cited upstream sources, used as reduction targets not as answers), Mathematica transliteration (independent derivations confirmed), tautological F_general construction (real algebraic identity), tautological H_F cross-check (two genuinely different paths), and engine drift (none). The stage-23 internal label vs. stage-040 paper label is cosmetic. The Stage-040 audit holds up.

## Self-test notes

Walked through the math by hand: (i) `z0^2*(delta+xi) + (q*z0)^2*xi = z0^2*(delta+(1+q^2)*xi)`, so script L59 alpha_req equals paper's eq. alpha-req; (ii) the split-U substitution `q -> -sqrt(2)/3 R_U`, `eta -> (2/9) R_U` into `F_expected` yields paper's `F_U` after pulling factors of 9 through (numerator gains `1/6561`, denominator gains `1/81`, leaving the explicit `1/(81*(1-xi))` factor in paper's eq. FU); (iii) at R_U=1, `9+2R_U^2 = 11 = 9+2R_U`, giving `(9d+11xi)^4` in F_U numerator and the Stage-18 denominator, matching `F_stage18`. Verified F_general parity in q (depends only on `q^2` via `(1+q^2)xi` and `q^2*xi^2`), so the negative sign of `q_U` is harmless. The H_F cross-check uses two different intermediate expressions (one with q,eta already substituted, one parametrized by eps via q_U(1+eps), eta_U(1+eps)), so `expect_zero(HF - HF_direct)` is substantive. No trap from variable-independence in `sp.diff(EXPR, eps)`: in both paths `eps` enters the expression nontrivially before differentiation.
