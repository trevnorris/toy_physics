---
unit_id: 030
batch: II.1
auditor_model: claude-opus-4-7[1m]
audit_date: 2026-05-25T00:00:00Z
verdict: findings
stop_cold: null
findings_count: 1
paper_alignment: partial
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files: ["/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage030_selected_mode_normalization.md"]
  paper_appendix: present
---

# Audit unit 030 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_030.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage030_selected_mode_normalization.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part02.tex` (row 50 references stage 030; input at line 98; Part-II `N_Q^target` definition at lines 72-81)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage030_selected_mode_normalization_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage030_selected_mode_normalization_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage030_selected_mode_normalization_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage030_selected_mode_normalization_mathematica_audit.txt`

(Note: a previous v1 report at this same path raised three findings — two tautological_check items and one mathematica_transliteration. Inspection of the current scripts shows all three have been addressed: the tautological assertions at the old sympy:113 and sympy:144 have been removed and replaced with explanatory comments (sympy lines 113-117 and 148-155; mathematica lines 102-106 and 127-130), and the Mathematica script now derives lamMinus / lamPlus from `Eigenvalues[mMat]` with a basis-reordered explicit `2x2` matrix (mathematica lines 55-73), with an explicit inline comment at lines 48-54 calling out the de-transliteration. This new audit is the v2 paper-grounded pass on the current state.)

## What the paper claims

The stage card `\stagefield{Output}` says: "Stage~030 outputs the lower eigenvalue \eqref{eq:app-stage030-lambda-minus}, the selected overlap \eqref{eq:app-stage030-s-minus}, the Hellmann--Feynman identity \eqref{eq:app-stage030-HF}, the selected prefactor \eqref{eq:app-stage030-P0minus}, and the target \eqref{eq:app-stage030-selected-target}." The deliverables are:

1. Closed form `lambda_- = (A + B - alpha_0 sigma - R_lambda)/2` with `B = A + Delta K_ax` and `R_lambda = sqrt((Delta K_ax + alpha_0 delta_kappa)^2 + 4 alpha_0^2 KappaProd)`.
2. Closed form `s_- = (v.e_-)^2 = (1/2)[sigma + ((Delta K_ax + alpha_0 delta_kappa) delta_kappa + 4 alpha_0 KappaProd)/R_lambda]`.
3. Hellmann-Feynman identity `s_- = (v.e_-)^2 = -d lambda_-/d alpha_0`. The card's `\stagefield{Checks}` item 1 enumerates this verbatim: "Differentiating the eigenvalue with respect to alpha_0 gives `d lambda_-/d alpha_0 = -e_-^T v v^T e_- = -s_-`."
4. Normalized-response coefficients `u_{2,-} = -D_{-2}/D_{-0}`, `u_{4,-} = (D_{-2}^2 - D_{-0} D_{-4})/D_{-0}^2`, `Gamma_{5,-} = C_{5,-}/D_{-0}`.
5. Selected odd-coefficient identity `Gamma_{5,-} = (a^5/(27 c_s^5)) beta_0 s_-/lambda_-`.
6. Selected static prefactor `P_{0,-} = beta_0 s_-/lambda_- = -beta_0 d ln lambda_-/d alpha_0`.
7. Selected target `mhat_-^2 P_{0,-} = N_Q^target` with `N_Q^target = 54 G c_s^5/(5 a^5 c^5)` (Part-II definition at appendix line 78).

The notes file expands all of these and adds in section 6 the exact determinant identity `lambda_- lambda_+ = A B - alpha_0(B kappa_0^2 + A kappa_1^2)`, which the paper card does not box but the script verifies.

## What the script claims to verify

Both scripts share a four-part structure. **Part I** posits the abstract operator `D_-(omega) = D_0 + D_2 omega^2 + D_4 omega^4 - i C_5 omega^5`, series-expands `Y_- = D_0/D_-` to order omega^5, and checks that the omega^2, omega^4, and imaginary omega^5 coefficients reproduce the textbook closed forms (paper deliverable 4). **Part II** produces `lambda_-` (SymPy: typed directly at line 73; Mathematica: derived from `Eigenvalues[mMat]` lines 55-68 with branch selection against the typed expression), computes `s_minus_hf = -d lambda_-/d alpha`, defines `s_minus_closed` from the paper's explicit overlap formula, asserts `s_minus_hf - s_minus_closed = 0`, and verifies the alpha to 0 limit `s_- -> kappa_0^2`. **Part III** defines `C5_sel = beta5 s_-`, `Gamma5_sel = C5_sel/lambda_-`, `P0_sel = beta0 s_-/lambda_-`, asserts `P0_sel + beta0 d ln lambda_-/d alpha = 0`, and checks the determinant identity `lambda_- lambda_+ = A(A+DK) - alpha((A+DK)x0 + A x1)`. **Part IV** introduces the physical constants `G5_phys = a^5/(27 c_s^5)`, `N_Q^target = 54 G c_s^5/(5 a^5 c^5)`, and asserts the equivalence `cond1 - G5_phys cond2 = 0` where `cond1 = mhat^2 G5_phys P0_sel - 2G/(5 c^5)` and `cond2 = mhat^2 P0_sel - N_Q^target`.

## Paper to script cross-check

| Paper deliverable | Script coverage |
|---|---|
| (1) `lambda_-` closed form | SymPy: assumed (line 73 types the closed form). Mathematica: derived from `Eigenvalues[mMat]` (lines 55-68) and matched against the typed form via branch selection. Cross-engine `match`. |
| (2) `s_-` closed form | partial — neither engine derives the closed form from the eigenvector. Both define `s_minus_closed` (sympy line 82; math line 79) and verify it against the derivative (HF). The closed-form-vs-derivative loop is verified; the closed-form-vs-`(v.e_-)^2` anchor is not. |
| (3) HF identity `s_- = (v.e_-)^2 = -d lambda_-/d alpha` | partial — both engines verify `closed_form - (-d lambda_-/d alpha) = 0` (sympy line 87; math line 83). Neither engine computes `(v.e_-)^2` directly from the eigenvector of the loaded wall block, so the leftmost equality `(v.e_-)^2 = -d lambda_-/d alpha` (the actual HF theorem statement listed verbatim in the paper's `\stagefield{Checks}`) is never exercised. See F1. |
| (4) `u_{2,-}, u_{4,-}, Gamma_{5,-}` | match — Part I series expansion + coefficient check (sympy lines 53-55, math lines 39-41). |
| (5) `Gamma_{5,-} = (a^5/(27 c_s^5)) beta_0 s_-/lambda_-` | match (by composition) — Part III uses abstract `G5` (`beta5 = G5*beta0`, `Gamma5_sel = beta5 s/lambda_-`), Part IV substitutes `G5_phys = a^5/(27 c_s^5)`. No single-line assertion of deliverable 5 stands alone, but the definitional chain produces it; the script acknowledges this with the inline note at sympy lines 113-117 / math lines 102-106. |
| (6) `P_{0,-} = beta_0 s/lambda_- = -beta_0 d ln lambda_-/d alpha` | match — definition + HF identity check at sympy line 118, math line 107. |
| (7) `mhat_-^2 P_{0,-} = N_Q^target` equivalent to `mhat_-^2 Gamma_{5,-} = 2G/(5 c^5)` | match (algebraic) — sympy line 142, math line 124 establish `cond1 - G5_phys cond2 = 0`. The check reduces to the prefactor identity `G5_phys * N_Q^target = 2G/(5 c^5)` (i.e., `54/(27*5) = 2/5`); this is a substantive consistency check on the numeric prefactors `27` and `54` themselves but does not exercise `P0_sel` or `mhat`. |
| (extra) determinant identity `lambda_- lambda_+ = AB - alpha(B kappa_0^2 + A kappa_1^2)` | match — sympy line 125, math line 110. This is a notes-side claim (notes section 6); the paper appendix references the structure implicitly. Useful sanity check, not an `extra` mismatch. |

Paper alignment is **partial** because deliverable 3 (the actual Hellmann-Feynman theorem statement `(v.e_-)^2 = -d lambda_-/d alpha`) is only half-verified.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored? |
|---|---|---|---|---|---|
| A1 | sympy | 53 | `Ysel.coeff(w,2) - (-D2/D0) == 0` | (4) `u_{2,-}` | yes |
| A2 | sympy | 54 | `Ysel.coeff(w,4) - (D2^2-D0*D4)/D0^2 == 0` | (4) `u_{4,-}` | yes |
| A3 | sympy | 55 | `im(Ysel).coeff(w,5) - C5/D0 == 0` | (4) `Gamma_{5,-}` | yes |
| A4 | sympy | 87 | `s_minus_hf - s_minus_closed == 0` | (3) HF: closed-form vs derivative leg | partial (eigenvector leg not exercised) |
| A5 | sympy | 91 | `s_minus_closed.subs(alpha,0) - x0 == 0` | (2) weak-loading limit `s_- -> kappa_0^2` | yes (non-trivial limit) |
| A6 | sympy | 118-121 | `P0_sel + beta0 * d log(lam_-)/d alpha == 0` | (6) `P_{0,-} = -beta_0 d ln lambda_-/d alpha` | yes (algebraically equivalent to A4) |
| A7 | sympy | 125 | `lam_- lam_+ - (A(A+DK) - alpha T0) == 0` | (notes 6) determinant identity | yes |
| A8 | sympy | 142-145 | `cond1 - G5_phys cond2 == 0` | (7) target equivalence | partial (collapses to constant identity `G5_phys*N_Q^target = 2G/(5c^5)`; independent of `P0_sel` and `mhat`) |
| B1 | math | 39 | `Coefficient[ySel,w,2] - (-d2/d0) == 0` | (4) `u_{2,-}` | yes |
| B2 | math | 40 | `Coefficient[ySel,w,4] - (d2^2-d0 d4)/d0^2 == 0` | (4) `u_{4,-}` | yes |
| B3 | math | 41 | `Coefficient[Im[ySel],w,5] - c5/d0 == 0` | (4) `Gamma_{5,-}` | yes |
| B4 | math | 83 | `sMinusHF - sMinusClosed == 0` | (3) HF closed-form vs derivative | partial (same as A4) |
| B5 | math | 85 | `(sMinusClosed/.alpha->0) - x0 == 0` | (2) weak-loading limit | yes |
| B6 | math | 107 | `p0Sel + beta0 D[Log[lamMinus],alpha] == 0` | (6) `P0_-` HF identity | yes |
| B7 | math | 110 | `lamMinus*lamPlus - (a(a+dK) - alpha t0) == 0` | (notes 6) det identity | yes |
| B8 | math | 124 | `cond1 - g5Phys*cond2 == 0` | (7) target equivalence | partial (same as A8) |

The most load-bearing identities (A4/B4 and A6/B6) are algebraically equivalent — they each say `s_minus_closed = -d lambda_-/d alpha`. A8/B8 reduce to a pure prefactor identity once the `P0_sel`-dependent terms cancel. The strongest substantive assertions are the Part I series-expansion checks (A1-A3, B1-B3), the alpha to 0 limit (A5/B5), and the determinant identity (A7/B7). The Mathematica-side independent derivation of `lamMinus` from `Eigenvalues[mMat]` (lines 55-73) is the strongest evidence that the typed closed form is correct.

## Findings

### F1 — insufficient_verification

**Severity:** medium

**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage030_selected_mode_normalization_sympy_audit.py:73-91`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage030_selected_mode_normalization_mathematica_audit.wl:55-85`

**What's wrong:**
The paper card lists as its first explicit check (`\stagefield{Checks}` item 1) the Hellmann-Feynman identity *in full*: "Differentiating the eigenvalue with respect to alpha_0 gives `d lambda_-/d alpha_0 = -e_-^T v v^T e_- = -s_-`." The chain has three sides: the derivative, the eigenvector projection `(v.e_-)^2`, and the closed-form `s_-`. The stage card `\stagefield{Output}` likewise boxes `s_- = (v.e_-)^2 = (1/2)[sigma + ...]` (eq:app-stage030-s-minus); the leftmost equality is asserted as a paper-side output, not a definition.

Both scripts verify only the rightmost two sides:
- sympy:87 `expect_zero("selected overlap: HF - closed form", s_minus_hf - s_minus_closed)` where `s_minus_hf = -sp.diff(lam_minus, alpha)` (line 81) and `s_minus_closed` is the typed closed form (lines 82-85).
- math:83 `expectZero["selected overlap: HF - closed form", sMinusHF - sMinusClosed]` with the same structure (lines 78-82).

Neither engine ever constructs `e_-` (the eigenvector of the loaded wall block) or computes the inner product `(v.e_-)^2` directly. The Mathematica script already builds the explicit matrix `mMat` (lines 55-56) and calls `Eigenvalues[mMat]` (line 63); adding one `Eigenvectors[mMat]` call there to verify `(v.e_-)^2 = sMinusClosed` would close the gap with minimal change. The SymPy script does not construct the matrix at all, but can do so following the Mathematica basis ordering.

The leftmost equality `(v.e_-)^2 = -d lambda_-/d alpha` is the *physical* content of paper item 3 (the actual Hellmann-Feynman theorem applied to this specific operator); the rightmost equality (closed form vs derivative) is a smaller consistency check between two formulas that are both downstream of the same algebra.

**Why this matters:**
Without the eigenvector-side check, the audit cannot detect, e.g., a sign or normalization error in the paper's identification of `s_-` with `(v.e_-)^2`. The paper invokes HF as a theorem rather than verifying it on the explicit operator; the audit should plug that gap, especially since `lambda_-` itself is independently derived from `mMat` on the Mathematica side. Stage 031 then uses `s_-` to prove monotonicity of `P_{0,-}` (per the stage card's `\stagefield{Downstream use}`), so a hidden mis-identification of `s_-` with `(v.e_-)^2` would propagate. The risk is low (the closed-form-vs-derivative check would catch most realistic algebra mistakes) but the gap is real and inexpensive to close.

**Required change:**
Add an explicit eigenvector projection check on each side. The basis ordering of `mMat` at math:55-56 places the kappa_1-mode in row 1 (diagonal `a + dK - alpha*x1 = B - alpha*kappa_1^2`) and the kappa_0-mode in row 2 (diagonal `a - alpha*x0 = A - alpha*kappa_0^2`); the off-diagonal `-alpha*Sqrt[x0*x1] = -alpha*kappa_0*kappa_1`. In this basis, the loading direction is `v = {kappa_1, kappa_0} = {Sqrt[x1], Sqrt[x0]}`. See directive F1 for the precise insertion blocks.

**Verification:**
After re-running, the new line `HF eigenvector check = 0` (sympy) and `PASS: HF eigenvector check` (mathematica) should appear in the captured transcripts. The simplified `(v.e_-)^2` expression should equal `sMinusClosed`/`s_minus_closed`. All existing assertions remain unchanged.

## Independent-derivation check (Mathematica)

The `.wl` script is **not** a line-by-line transliteration of the `.py` script. Evidence:

1. Mathematica Part II (lines 55-73) constructs the explicit loaded `2x2` wall matrix `mMat = {{a + dK - alpha*x1, -alpha*Sqrt[x0*x1]}, {-alpha*Sqrt[x0*x1], a - alpha*x0}}` and routes `lamMinus`/`lamPlus` through `Eigenvalues[mMat]`, with branch selection by matching against the typed closed form. SymPy (line 73) directly types `lam_minus = (2A + DK - alpha sigma - R)/2`. The Mathematica version is a genuinely independent re-derivation from the operator. The inline comment at math lines 48-54 explicitly flags this as the de-transliteration step.

2. Verification of the matrix: trace `2a + dK - alpha*(x0+x1)` matches paper `A + B - alpha sigma`; det `a*(a+dK) - alpha*((a+dK)*x0 + a*x1) = AB - alpha(B kappa_0^2 + A kappa_1^2)` matches the notes section-6 determinant identity (and is independently verified by the script's own A7/B7 check). The discriminant `T^2 - 4D` reduces to `dK^2 + 2 dK alpha (x0-x1) + alpha^2 sigma^2 = R^2`, matching the paper's `R_lambda^2`. So `mMat` is the paper's `K_eff^(0) = diag(A,B) - alpha v v^T` up to a 1-to-2 basis reordering (row 1 is the kappa_1-mode, row 2 is the kappa_0-mode), and `Eigenvalues[mMat]` independently reproduces the closed form.

Elsewhere (Part I, Part III, Part IV), the two scripts share the same algebraic structure because the claims themselves are algebraic identities on abstract symbols (D0, D2, D4, C5, beta0, G5, etc.) — there is no first-principles re-derivation either engine could mount that the other could avoid. The Mathematica script's basis-reordered matrix in Part II is sufficient to clear the `mathematica_transliteration` bar for this stage.

## Engine cross-check

The two engines' final outputs agree on every check:

| sympy check | sympy output | math check | math output | residual |
|---|---|---|---|---|
| u2 coefficient check | line 11 | u2 coefficient check | line 6 (PASS line 7) | 0 |
| u4 coefficient check | line 12 | u4 coefficient check | line 8 (PASS line 9) | 0 |
| Gamma5 coefficient check | line 13 | Gamma5 coefficient check | line 10 (PASS line 11) | 0 |
| selected overlap: HF - closed form | line 30 | selected overlap: HF - closed form | line 18 (PASS line 19) | 0 |
| weak-loading overlap limit | line 47 | weak-loading overlap limit | line 21 (PASS line 22) | 0 |
| P0_sel + beta0 d log lam/d alpha | line 100 | (same) | line 30 (PASS line 31) | 0 |
| det identity | line 101 | det identity | line 32 (PASS line 33) | 0 |
| normalization equivalence | line 106 | normalization equivalence | line 38 (PASS line 39) | 0 |

Closed forms also agree: SymPy `lambda_-` (output lines 18-23) and Mathematica `lambda_-` (output line 16) are the same expression in different surface forms; `s_-` (sympy output lines 31-46) and `s_-` (math output line 20) likewise; `lambda_req` (sympy output lines 107-122) and `lambda_req` (math output line 40) likewise.

Outputs are fresh: sympy script and output have the same mtime (1779405337 / 1779405424; output newer); Mathematica script and output (1779405337 / 1779405433; output newer). No stale_output finding.

No engine disagreement.

## Verdict justification

Verdict is `findings` with one medium-severity `insufficient_verification` issue. The scripts faithfully execute the algebra surrounding the closed forms for `lambda_-`, `s_-`, `P_{0,-}`, and the target equivalence; the Mathematica script independently re-derives `lambda_-` from the explicit wall matrix. All numeric prefactors in Part IV (`G5_phys = a^5/(27 c_s^5)`, `N_Q^target = 54 G c_s^5/(5 a^5 c^5)`) match the paper card and the Part-II appendix definition (eq:app-part02-NQ-target). The single gap is that the paper-side equality `(v.e_-)^2 = -d lambda_-/d alpha` (Hellmann-Feynman theorem on the explicit `K_eff`) is never exercised — both scripts only verify the closed-form-vs-derivative leg. This is fixable in a handful of lines on each side and does not invalidate downstream stages, so `stop_cold` is null.

Attacks attempted that failed: (i) sign flip on `R` inside `lambda_-`; the explicit Mathematica derivation pins the sign via branch selection. (ii) factor-of-2 error in `s_-`; the alpha to 0 limit `s_- -> kappa_0^2` would not satisfy a wrong overall factor (a factor of 2 error would give `2 kappa_0^2` or `kappa_0^2/2`). (iii) wrong constant in `N_Q^target`; Part IV would report a nonzero residue `2G/(5 c^5) - G5_phys N_Q^target` if either prefactor were off. (iv) wrong matrix entries in Mathematica; trace, determinant, and discriminant all match the paper's formulas exactly. (v) `simplify` masking errors under positivity assumptions; the symbols `A, DK, alpha, x0, x1` are restricted but the assumptions match the paper's stated setup (positive wall stiffnesses, nonnegative loading, positive kappa^2 values). (vi) basis-reordering trap in Mathematica; trace/det invariant under basis swap, so eigenvalues are independent of the row ordering.

## Self-test notes

- **Variable independence:** `sp.diff(lam_minus, alpha)` and `D[lamMinus, alpha]` are well-formed; `lam_minus(A, DK, alpha, x0, x1)` genuinely depends on `alpha` via `-alpha*sigma` and the `R(alpha)` term. The HF check is non-trivial; no zero-by-construction trap. The proposed F1 eigenvector check projects `v = {Sqrt[x1], Sqrt[x0]}` (in the Mathematica basis) onto an eigenvector that does depend on `alpha, x0, x1, A, DK`, so the resulting expression is a genuine function of those symbols.
- **Symmetry/parity:** no oscillatory integrals; not applicable.
- **Trivial-case pre-check:** at `alpha = 0`, `R = DK` (since `DK > 0`), `s_- = (sigma + DK*(x0-x1)/DK)/2 = (sigma + delta_kappa)/2 = x0`; matches the existing A5/B5 check. Also at `alpha = 0`, the matrix `mMat` becomes `diag(a+dK, a) = diag(B, A)`, eigenvectors are `{1,0}` and `{0,1}`, so `(v.e_-)^2` with `v = {Sqrt[x1], Sqrt[x0]}` and `e_-` being the eigenvector at the lower eigenvalue `a` (i.e., `{0,1}`) gives `(Sqrt[x0])^2 = x0`. The proposed F1 check reduces correctly at `alpha = 0`.
- **Path specifications:** directive points to existing files; no new script files prescribed. `.py` lives in `scripts/`, `.wl` lives in `mathematica/`.
- **Paper round-trip:** the proposed eigenvector check uses `v = {Sqrt[x1], Sqrt[x0]}` in the Mathematica basis (and `v = {sqrt(x1), sqrt(x0)}` in SymPy if it adopts the same matrix ordering) consistent with the paper's `v = (kappa_0, kappa_1)^T` after the explicit 1-to-2 basis swap; no new constants introduced; no new paper_misalignment risk.
- **Codex pre-check:** the directive only specifies *adding* lines, not modifying existing ones; no risk of clobbering working checks.
