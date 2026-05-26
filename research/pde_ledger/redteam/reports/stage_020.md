---
unit_id: 020
batch: I.2
auditor_model: claude-opus-4-7[1m]
audit_date: 2026-05-25T00:00:00Z
verdict: findings
stop_cold: null
findings_count: 2
paper_alignment: partial
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files: []
  paper_appendix: present
---

# Audit unit 020 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_020.tex`
- notes: (none)
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part01.tex` (row for stage 020 at line 62)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage020_parent_throat_action_weak_axisym_packet_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage020_parent_throat_action_weak_axisym_packet_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage020_parent_throat_action_weak_axisym_packet_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage020_parent_throat_action_weak_axisym_packet_mathematica_audit.txt`

## What the paper claims

The paper card states this stage proves the weak-axisymmetric parent throat-action packet. Specifically, it defines (eq. stage020-wall-slopes) the wall-slope integrals `delta M_Sigma = int dw delta mu_eta beta_2^2` and `delta K_Sigma = int dw [delta T_w (beta_2')^2 + (delta K_eta + 6 delta T_Omega) beta_2^2]`; (eq. stage020-denominator-slopes) the denominator slopes `D_{01} = delta K_Sigma - B_{01} - Z_{01}`, `D_{21} = -(delta M_Sigma + B_{21} + Z_{21})`, `D_{41} = -(B_{41} + Z_{41})`; (eq. stage020-gates) the live gates `K_1 = D_{21} + D_{01}/9` and `H_even = D_{41} - (2/3) D_{21} - D_{01}/27`; (eq. stage020-xi) `Xi_1 = N_{01}/N_0 - D_{01}/D_0`. The exact wall-slope solve (eq. stage020-wall-slope-solve) is `delta K_Sigma = B_{01} + Z_{01} + 27(B_{41}+Z_{41})` and `delta M_Sigma = -(B_{21}+Z_{21}) + 3(B_{41}+Z_{41})`; the residual prefactor (eq. stage020-residual-xi) is `Xi_1 = N_{01}/N_0 - 27(B_{41}+Z_{41})/(K_Sigma - B_0 - Z_0)`. Output paragraph: "Stage~020 exports the weak-axisymmetric parent packet \eqref{eq:stage020-wall-slopes}--\eqref{eq:stage020-residual-xi}." The part-appendix row at line 62 of `stage_appendix_part01.tex` repeats this as "Weak-axisymmetric packet exported from the parent throat-action bundle." No spherical-harmonic / Gaunt lane ratios are mentioned in the paper card or appendix row.

## What the script claims to verify

The script docstring states it "verifies the exact wall-slope solve of the even gates once the parent throat action is combined with the already-derived support and mixed slope packets." Concretely, the script (a) defines the same `D_{01}, D_{21}, D_{41}, K_1, H_even, Xi_1` as the paper; (b) computes Y20-squared overlap ratios `lambda_{20}, lambda_{21}, lambda_{22}` (and same-sign vanishing checks) from `sympy.physics.wigner.gaunt`; (c) asserts the even-gate Jacobian determinant equals `1/27`; (d) solves `K_1 = H_even = 0` for `dKSigma, dMSigma` and asserts the closed forms match `B_{01}+Z_{01}+27(B_{41}+Z_{41})` and `-(B_{21}+Z_{21})+3(B_{41}+Z_{41})`; (e) checks that on the compensated branch, `D_{01} = 27(B_{41}+Z_{41})`, `D_{21} = -3(B_{41}+Z_{41})`, `D_{41} = -(B_{41}+Z_{41})`, and `Xi_1 = N_{01}/N_0 - 27(B_{41}+Z_{41})/(K_Sigma - B_0 - Z_0)`. The Mathematica script does the same with independently-built `GauntIntegral` from `ThreeJSymbol`.

## Paper <-> script cross-check

| Paper deliverable | Script-side check | Status |
|---|---|---|
| Wall-slope integral definitions (eq14) `delta M_Sigma, delta K_Sigma` | Treated as free symbols `dKSigma, dMSigma`; integral form not exercised | partial (the integral form is a definition, not an identity to verify; treating as free symbols is acceptable) |
| Denominator slopes (eq22) `D_{01}, D_{21}, D_{41}` | sympy lines 36-39; wl lines 77-80 | match |
| Live gates (eq28) `K_1, H_even` | sympy lines 41-42; wl lines 81-82 | match |
| Live gate `Xi_1` (eq32) | sympy line 43; wl line 83 | match |
| Wall-slope solve (eq42) `delta K_Sigma, delta M_Sigma` | sympy lines 58, 65-66; wl lines 96, 110, 114 | match |
| Residual `Xi_1` (eq49) | sympy line 70; wl lines 140-149 | match |
| Even-gate Jacobian det = 1/27 | sympy line 56; wl line 89 | extra (not stated in paper card, but a clean implication; auxiliary) |
| Y20 lane ratios `lambda_{20}=1, lambda_{21}=1/2, lambda_{22}=-1` and same-sign vanishing | sympy lines 45-50; wl lines 34-67 | extra (paper card and appendix do not mention Y20/Gaunt) |

The dominant pattern is that the algebraic-core claims of the paper are all faithfully exercised. The script additionally verifies Y20 overlap ratios that the paper card does not anchor — this is the source of the `paper_alignment: partial` flag.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 48 | `assert_zero('Y20 overlap lane 20', lam20 - 1)` | none (lam20 is identically 1 by construction: gaunt(2,2,2,0,0,0)/gaunt(2,2,2,0,0,0)) | no (tautological) |
| A2 | sympy | 49 | `assert_zero('Y20 overlap lane 21', lam21 - 1/2)` | none (paper card omits Y20) | yes within script's own claim, but no paper anchor |
| A3 | sympy | 50 | `assert_zero('Y20 overlap lane 22', lam22 + 1)` | none (paper card omits Y20) | yes within script's own claim, but no paper anchor |
| A4 | sympy | 25 (inside `real_y20_square_ratio`) | `if same_sign != 0: raise` | none (paper card omits Y20) | yes (non-trivial selection-rule check) |
| A5 | sympy | 56 | `assert_zero('even-gate solve determinant', coeff_matrix.det() - 1/27)` | implicit; supports eq28 invertibility | yes |
| A6 | sympy | 65 | `assert_zero('dKSigma closed form', sol[dKSigma] - (B01+Z01+27(B41+Z41)))` | eq42 left | yes |
| A7 | sympy | 66 | `assert_zero('dMSigma closed form', sol[dMSigma] - (-(B21+Z21)+3(B41+Z41)))` | eq42 right | yes |
| A8 | sympy | 67 | `assert_zero('compensated D01', D01_comp - 27(B41+Z41))` | eq22+eq42 corollary, consistent with eq49 numerator | yes |
| A9 | sympy | 68 | `assert_zero('compensated D21', D21_comp + 3(B41+Z41))` | eq22+eq42 corollary | yes |
| A10 | sympy | 69 | `assert_zero('compensated D41', D41_comp + B41 + Z41)` | eq22 | yes |
| A11 | sympy | 70 | `assert_zero('compensated Xi1', Xi1_comp - (N01/N0 - 27(B41+Z41)/(KSigma-B0-Z0)))` | eq49 | yes |
| M1 | wl | 49-51 | `If[!TrueQ[FullSimplify[lambda0 - 1] === 0], ... Exit[1]]` | none (tautological: lambda0 = GauntIntegral[2,2,2,0,0,0]/overlapBase where overlapBase = GauntIntegral[2,2,2,0,0,0]) | no (tautological) |
| M2 | wl | 53-55 | `... lambda1 - 1/2 ...` | none (paper card omits Y20) | yes within script claim |
| M3 | wl | 57-59 | `... lambda2 - (-1) ...` | none (paper card omits Y20) | yes within script claim |
| M4 | wl | 61-63 | `... crossOne - 0 ...` | none (paper card omits Y20); selection-rule | yes |
| M5 | wl | 65-67 | `... crossTwo - 0 ...` | none (paper card omits Y20); selection-rule | yes |
| M6 | wl | 91-93 | even-gate Jacobian determinant - 1/27 | implicit | yes |
| M7 | wl | 98-100 | unique-solution count | implicit | yes |
| M8 | wl | 110-112 | dKSigma closed form | eq42 | yes |
| M9 | wl | 114-116 | dMSigma closed form | eq42 | yes |
| M10 | wl | 127-137 | D01, D21, D41 deficits | eq22+eq42 corollary | yes |
| M11 | wl | 146-149 | Xi1 residual | eq49 | yes |

## Findings

### F1 — paper_misalignment

**Subtype:** paper_missing_script_claim

**Severity:** medium

**Files:**
- `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_020.tex:1-56` (no Y20 / Gaunt / `lambda` lane mention anywhere in card)
- `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part01.tex:62` (appendix row: "Weak-axisymmetric packet exported from the parent throat-action bundle" — no Y20 mention)
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage020_parent_throat_action_weak_axisym_packet_sympy_audit.py:45-50`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage020_parent_throat_action_weak_axisym_packet_mathematica_audit.wl:34-67`

**What's wrong:**

The script (both engines) computes and asserts the Y20-squared lane overlap ratios
```
lambda_{20} = 1,   lambda_{21} = 1/2,   lambda_{22} = -1
```
and the vanishing of same-sign m-cross terms. The paper card for stage 020 makes no reference to spherical harmonics, Y20 overlaps, Gaunt integrals, or these lane ratios. The card's content is purely the one-dimensional wall-slope solve in the `beta_2(w)` profile and the algebraic identities listed in eq. stage020-wall-slopes through eq. stage020-residual-xi.

Paper card quote (entire Output paragraph, `stage_020.tex:53-55`):
> "Stage~020 exports the weak-axisymmetric parent packet \eqref{eq:stage020-wall-slopes}--\eqref{eq:stage020-residual-xi}."

Script quote (`sympy_audit.py:45-50`):
```
lam20 = real_y20_square_ratio(0)
lam21 = real_y20_square_ratio(1)
lam22 = real_y20_square_ratio(2)
assert_zero('Y20 overlap lane 20', lam20 - 1)
assert_zero('Y20 overlap lane 21', lam21 - sp.Rational(1, 2))
assert_zero('Y20 overlap lane 22', lam22 + 1)
```

**Why this matters:**

Either (a) the paper card should anchor these Y20 ratios (the "weak-axisymmetric" content of the stage title suggests they are conceptually load-bearing somewhere), or (b) the script should drop them (they belong to a different stage). Without an explicit paper anchor, the script's load-bearing assertion of `lambda_{21} = 1/2` and `lambda_{22} = -1` cannot be cross-checked against the derivation. The direction of resolution is the user's call.

**Required change:**
(See directive `## Resolve before fix_loop` — codex does not act on this finding.)

**Verification:**
After user resolves: either (a) the paper card is amended to add a paragraph or equation citing the Y20 lane ratios, in which case the script checks remain and the row in this report flips to `match`, or (b) the script section at lines 45-50 (sympy) and lines 34-67 (wl) is removed, in which case the report row drops out entirely.

### F2 — tautological_check

**Severity:** low

**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage020_parent_throat_action_weak_axisym_packet_sympy_audit.py:20-26, 45, 48`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage020_parent_throat_action_weak_axisym_packet_mathematica_audit.wl:34-35, 49-51`

**What's wrong:**

In `real_y20_square_ratio(m)` (sympy, lines 20-26):
```
base = sp.simplify(gaunt(2, 2, 2, 0, 0, 0))
...
return sp.simplify((sp.Integer(-1) ** m) * gaunt(2, 2, 2, 0, m, -m) / base)
```
For `m = 0`, the return value is `gaunt(2,2,2,0,0,0) / gaunt(2,2,2,0,0,0) = 1` by construction. Then line 48 asserts `lam20 - 1 == 0`, which is identically true regardless of what `gaunt` returns. The assertion cannot detect any error in the Gaunt code.

The same defect is present in Mathematica (`mathematica_audit.wl:34-35, 49-51`):
```
overlapBase = GauntIntegral[2, 2, 2, 0, 0, 0];
lambda0 = FullSimplify[(-1)^0 GauntIntegral[2, 2, 2, 0, 0, 0]/overlapBase];
...
If[!TrueQ[FullSimplify[lambda0 - 1] === 0], Print["FAIL: M1 lambda_0"]; Exit[1]];
```
`lambda0` is `overlapBase / overlapBase = 1` by construction.

**Why this matters:**

The "M1 lambda_0 residual = 0" line in the saved Mathematica output and the silently-passing sympy assertion mask the fact that nothing is being checked. If the contributor were to ever swap or alter the `gaunt(2,2,2,0,0,0)` definition, the `lam20`/`lambda0` assertion would still pass.

**Required change:**

Either (a) replace the tautological self-ratio with a concrete check of `gaunt(2,2,2,0,0,0)` against its closed-form value `sqrt(5/(4*pi)) * 2/7` (this anchors the base), or (b) restructure `lam20` so it is not the ratio of an expression to itself — for instance, compute it from a direct angular integral and check the resulting numeric ratio. If the easier fix is preferred, anchor `base` to its closed value before forming the ratios. NB: this finding is contingent on F1 — if the user decides the Y20 block should be removed from the script entirely (F1 direction b), F2 disappears automatically.

**Verification:**
After Codex applies, `lam20 - 1` (sympy) and `lambda0 - 1` (Mathematica) should each be the result of an assertion that *could fail* if `gaunt(2,2,2,0,0,0)` were perturbed, e.g. by anchoring `base` to its closed-form numerical value.

## Independent-derivation check (Mathematica)

The `.wl` script does NOT transliterate the SymPy `.py`. SymPy imports `gaunt` directly from `sympy.physics.wigner` (one library call); Mathematica builds `GauntIntegral` from first principles using `ThreeJSymbol` and a `gauntWeight = Sqrt[(2 la + 1)(2 lb + 1)(2 lc + 1)/(4 Pi)]` factor (`mathematica_audit.wl:21-30`). The angular machinery is genuinely independent. The algebraic core (D-slopes, K_1, H_even, Xi_1, Jacobian determinant, Solve) is structurally parallel because the underlying algebra is small and there is essentially one natural way to verify a 2x2 linear solve — this is constrained-by-the-problem similarity, not transliteration. Not flagged.

## Engine cross-check

Both engines emit `STATUS: PASS`. All algebraic-core residuals (D01, D21, D41 deficits, dKSigma, dMSigma, Xi1) reduce to zero in both engines. The Jacobian determinant equals `1/27` in both. The Y20 lane ratios `(1, 1/2, -1)` agree (modulo the tautology in `lambda_0` flagged in F2). No disagreement.

## Verdict justification

The script's algebraic-core content (eq. stage020-wall-slopes through eq. stage020-residual-xi) is fully verified by both engines and matches the paper card exactly. I attempted to break the wall-slope solve by hand: substituting `K_1 = 0` for `dKSigma` and back into `H_even = 0` reproduces `dMSigma = 3(B_{41}+Z_{41}) - (B_{21}+Z_{21})`, matching paper eq42. The Jacobian determinant is `(1/9)(2/3) - (-1)(-1/27) = 2/27 - 1/27 = 1/27`, matching the script. `Xi_1` on the compensated branch evaluates by direct substitution to `N_{01}/N_0 - 27(B_{41}+Z_{41})/(K_Sigma - B_0 - Z_0)`, matching paper eq49. Outputs are fresh (sympy script May 21 13:39 / output May 21 15:00; wl script May 21 13:41 / output May 21 15:02). The two findings are: (F1) the Y20 lane block at sympy:45-50 / wl:34-67 verifies content that has no anchor in the paper card or appendix row — `paper_misalignment` routes to user; (F2) within that block, `lambda_0 = 1` is a self-ratio and the assertion is tautological. Verdict: `findings`; not stop-cold because the algebraic-core math is sound and the Y20 question is a paper/script scope decision rather than a downstream-breaking math error.

## Self-test notes

- **Variable independence (Self-test 1):** Checked `D[K1, dKSigma]` and similar in both engines. `K1 = D21 + D01/9` depends on `dKSigma` through `D01` and on `dMSigma` through `D21`; `H_even` depends on both through `D21` (chain) and `D01`; `D41` is independent. All four derivatives are well-defined and non-trivial.
- **Symmetry/parity:** No unbounded-domain integrals in this script — all integrals are abstracted as free symbols `dKSigma, dMSigma`. Not applicable.
- **Trivial-case pre-check:** Set `B_{n1} = Z_{n1} = 0` for all `n`: closed forms collapse to `dKSigma = 0, dMSigma = 0`; `Xi_1 = N_{01}/N_0`. Consistent.
- **Paper round-trip:** Fix for F2 (if applied) does not alter any paper-side claim — it only strengthens the script-side anchor for the Y20 base value. F1 itself is paper-side and routes to user, not to Codex.
