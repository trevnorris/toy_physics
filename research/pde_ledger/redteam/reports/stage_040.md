---
unit_id: 040
batch: III.1
auditor_model: claude-opus-4-7[1m]
audit_date: 2026-05-22T00:00:00-06:00
verdict: findings
stop_cold: null
findings_count: 4
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: true
---

# Audit unit 040 red-team report

## Files reviewed

- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage040_generalized_selected_branch_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage040_generalized_selected_branch_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage040_generalized_selected_branch_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage040_generalized_selected_branch_mathematica_audit.txt`

## What the script claims to verify

The scripts claim five things about a generalized selected-branch normalization for a 2x2 baseline `diag(A0, A0(1+delta))` perturbed by a rank-1 loading `alpha z z^T` with loading ratio `q = z1/z0` and a separate source-direction ratio captured by `eta`. The docstring asserts: (1) exact closed-form eigenvalue `lam_minus = A0(1-xi)` and eigenvector `(1, r)` with `r = q xi/(delta+xi)`; (2) a two-vector normalization function `F_{q,eta}` that reduces to the Stage-18 form when source and loading are aligned; (3) a one-parameter family `F_U(xi,delta;R_U), G_U(xi,delta;R_U)` for the split-U continuum; (4) `R_U = 1` recovers Stage-18/19 closed forms; (5) the first-order Taylor expansion about `R_U = 1` is exact. The assertions verify algebraic identities between posited expressions and hardcoded closed forms; the eigenvalue/eigenvector statement is never tested against the actual matrix equation, and the Taylor-expansion check is a structural identity that cannot fail.

## Assertion inventory

| # | Script | Line | Form | Anchored to claim? |
|---|---|---|---|---|
| A1 | sympy | 66 | `simplify(r - xi*q/(delta+xi)) == 0` | partial — verifies algebra of posited `r`, not eigenvalue equation |
| A2 | sympy | 85 | `simplify(F_general - F_expected) == 0` | partial — both sides hardcoded/derived in script, no eigenvalue check feeds in |
| A3 | sympy | 86 | `simplify(G_general - G_expected) == 0` | partial — same as A2 |
| A4 | sympy | 103 | `simplify(F_U(R_U=1) - F_stage18) == 0` | partial — `F_stage18` is hardcoded, not imported |
| A5 | sympy | 104 | `simplify(G_U(R_U=1) - G_stage19) == 0` | partial — `G_stage19` is hardcoded |
| A6 | sympy | 118-121 | `expand(F_ratio - (1 + eps*HF)) == 0` | no — tautological: `F_ratio` is the Taylor series, `HF` is the derivative, equality holds by definition |
| A7 | sympy | 122-125 | `expand(G_ratio - (1 + eps*HG)) == 0` | no — same as A6 |
| B1 | mathematica | 50 | `expectZero[r - xi*q/(delta+xi)]` | partial — mirror of A1 |
| B2 | mathematica | 70 | `expectZero[fGeneral - fExpected]` | partial — mirror of A2 |
| B3 | mathematica | 71 | `expectZero[gGeneral - gExpected]` | partial — mirror of A3 |
| B4 | mathematica | 88 | `expectZero[(fU /. rU->1) - fStage18]` | partial — mirror of A4 |
| B5 | mathematica | 89 | `expectZero[(gU /. rU->1) - gStage19]` | partial — mirror of A5 |
| B6 | mathematica | 102 | `expectZero[fRatio - (1 + eps*hF)]` | no — same tautology as A6 |
| B7 | mathematica | 103 | `expectZero[gRatio - (1 + eps*hG)]` | no — same tautology as A7 |

## Findings

### F1 — insufficient_verification

**Severity:** high
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage040_generalized_selected_branch_sympy_audit.py:58-66`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage040_generalized_selected_branch_mathematica_audit.wl:44-50`

**What's wrong:**
The docstring (sympy lines 6-8) claims "the selected lower wall branch for a diagonal 2x2 baseline plus a rank-1 loading vector z has exact closed-form eigenvalue and eigenvector formulas." The script posits `lam_minus = A0*(1 - xi)` (line 58) and `alpha_req = A0*xi*(delta + xi)/(z0**2*(delta + xi) + (q*z0)**2*xi)` (line 59), then defines `r` as a particular algebraic combination (line 64) and asserts `r == q*xi/(delta+xi)` (line 66). At no point does the script:

- construct the 2x2 perturbed matrix `M = diag(A0, A0*(1+delta))` (with or without the rank-1 term `alpha z z^T`), or
- write the determinantal equation `det(M_perturbed - lam_minus * I) == 0`, or
- write the eigenvector residual `(M_perturbed - lam_minus * I) @ e_minus == 0`.

The assertion at line 66 verifies only an algebraic identity between two posited symbolic expressions. With the convention `M_perturbed = M - alpha z z^T`, the eigenvalue and eigenvector claims are in fact correct (I verified numerically at `q=1, delta=1, xi=1/2, z0=1, A0=1`: the matrix `[[1/8, -3/8], [-3/8, 9/8]]` has determinant 0 and null vector `(1, 1/3)`). But this is left implicit; the script never checks it. The Mathematica script mirrors the same omission.

**Why this matters:**
The "selected-branch" claim is the entire physical premise of stages 22-23. The script can pass while the underlying matrix relationship is silently wrong; any sign or convention error in `alpha_req` or `lam_minus` would slip through because the assertion only tests internal algebra. The downstream `F_general`, `G_general`, and split-U specializations all sit on this unverified foundation.

**Required change:**
Insert before line 68 (sympy) and before subbanner "2." (mathematica) an explicit eigenvalue/eigenvector residual check. Build the perturbed matrix and the eigenvector explicitly using SymPy `Matrix`/`MatMul` (or Mathematica equivalents), then assert each component of `(M_perturbed - lam_minus * I) @ e_minus` simplifies to zero. The convention is `M_perturbed = M - alpha_req * z * z^T` (subtracted rank-1), based on the script's sign choice in `r`.

**Verification:**
After the fix, the sympy script will have a new `expect_zero` block printing two residual components for the eigenvector equation, both reducing to 0; the mathematica script gains the analogous `expectZero` calls. If the residuals do not simplify to 0, the underlying premise is wrong and a separate finding is warranted.

### F2 — tautological_check

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage040_generalized_selected_branch_sympy_audit.py:106-125`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage040_generalized_selected_branch_mathematica_audit.wl:91-103`

**What's wrong:**
Section 23.4 ("Exact small-deformation expansion about the flat-U limit") computes:

```
F_ratio = sp.series(F_U.subs(R_U, 1 + eps) / F_stage18, eps, 0, 2).removeO()
HF      = sp.diff(F_U.subs(R_U, 1 + eps), eps).subs(eps, 0) / F_stage18
expect_zero("F_ratio - (1 + eps * HF)", sp.expand(F_ratio - (1 + eps * HF)))
```

By definition, `sp.series(f, eps, 0, 2).removeO() = f(0) + f'(0)*eps`. Here `f(0) = F_U(R_U=1)/F_stage18 = 1` (already checked in 23.3) and `f'(0) = HF`. Therefore `F_ratio = 1 + HF*eps` is true by construction of `series` and `diff`; the assertion `F_ratio - (1 + eps*HF) == 0` is an algebraic identity that cannot fail regardless of the physics. The same applies to the `G_ratio` check at lines 122-125 (and Mathematica lines 102-103). The substantive Taylor coefficient `H_F` itself is just *defined* as the derivative and is never compared against an independently derived value.

**Why this matters:**
The output advertises that the "first-order deformation around the flat-U limit is exact," but the assertion contributes zero physical information — it only verifies that SymPy's `series` and `diff` agree on the same expression, which is a CAS sanity check. A genuinely substantive verification would derive `H_F` from an independent route (e.g., a direct first-order perturbation of the eigenvalue problem when `R_U = 1 + eps`) and confirm the two computations match.

**Required change:**
Replace the tautological assertion at sympy lines 118-125 (and mathematica lines 102-103) with a substantive cross-check. Compute `HF_direct` and `HG_direct` by an *independent* path: take the `q_U`, `eta_U` substitutions at `R_U = 1 + eps`, expand `F_general` and `G_general` (the eigenvector-overlap construction from section 23.2) to first order in `eps` with `q -> q_U(1+eps)`, `eta -> eta_U(1+eps)`, and assert the resulting first-order coefficients equal the script's `HF`, `HG`. If an independent path is not feasible within the unit's scope, delete the tautological assertions (and the corresponding `H_F`, `H_G` prints) entirely and demote section 23.4 to an informational print of the leading-order ratio, since it does not exercise any claim.

**Verification:**
After the fix, the sympy output should show either (a) two new `HF - HF_direct == 0` and `HG - HG_direct == 0` checks both passing, where `HF_direct` is built from `F_general`/`G_general` rather than from `F_U` itself; or (b) section 23.4 contains only informational prints with no `expect_zero` calls on `F_ratio - (1 + eps*HF)`.

### F3 — mathematica_transliteration

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage040_generalized_selected_branch_mathematica_audit.wl:33-103`
- (compared against `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage040_generalized_selected_branch_sympy_audit.py:44-125`)

**What's wrong:**
The Mathematica script is a line-by-line transliteration of the SymPy script. Quoting three corresponding pairs:

Pair 1 — same posited `alpha_req`:
- SymPy line 59: `alpha_req = sp.simplify(A0 * xi * (delta + xi) / (z0**2 * (delta + xi) + (q * z0)**2 * xi))`
- Mathematica line 45: `alphaReq = FullSimplify[a0 xi (delta + xi)/(z0^2 (delta + xi) + (q z0)^2 xi), ...]`

Pair 2 — same hardcoded "expected" closed form for F:
- SymPy line 82: `F_expected = sp.simplify((delta + (1 + q**2) * xi)**2 * (delta + (1 + eta) * xi)**2 / ((1 - xi) * ((delta + xi)**2 + q**2 * xi**2)**2))`
- Mathematica lines 59-63: `fExpected = FullSimplify[(delta + (1 + q^2) xi)^2 (delta + (1 + eta) xi)^2/((1 - xi) ((delta + xi)^2 + q^2 xi^2)^2), ...]`

Pair 3 — same hardcoded Stage-18 target and same comparison structure:
- SymPy line 100: `F_stage18 = sp.simplify((9*delta + 11*xi)**4 / (81 * (1 - xi) * (9*delta**2 + 18*delta*xi + 11*xi**2)**2))`
- Mathematica lines 80-83: `fStage18 = FullSimplify[(9 delta + 11 xi)^4/(81 (1 - xi) (9 delta^2 + 18 delta xi + 11 xi^2)^2), ...]`

Every intermediate variable name is the same modulo casing (`alpha_req`↔`alphaReq`, `F_general`↔`fGeneral`, `z_overlap_sq`↔`zOverlapSq`, `F_stage18`↔`fStage18`, `H_F`↔`hF`, etc.). The algebraic ordering is identical, the same closed forms are hardcoded as targets in both engines, and the section banners (`23.1`/`1.`, `23.2`/`2.`) map one-to-one. Neither engine independently derives `alpha_req` from the eigenvalue equation (see F1), so both are echoing the same algebra. This violates the second-engine policy: the engines must derive results from the physical premise (the eigenvalue problem) independently, not import each other's algebra.

**Why this matters:**
A second engine that echoes the first engine's algebra cannot catch algebraic errors in the first engine. If `alpha_req` had a sign error or a misplaced factor, both scripts would carry the same error and both would "pass." The current pair only verifies that two CAS engines agree on simplification, not that the physics is right.

**Required change:**
Restructure the Mathematica script to derive the closed forms from the eigenvalue equation directly, using Mathematica's `Eigenvalues`/`Eigenvectors`/`Solve` on the explicit perturbed matrix. Concretely: build `Mmat = DiagonalMatrix[{a0, a0 (1 + delta)}] - alpha {z0, q z0}.Transpose[{{z0, q z0}}]`, then `Solve[Det[Mmat - lambda IdentityMatrix[2]] == 0 /. lambda -> a0 (1 - xi), alpha]` to derive `alphaReq` independently, and use `NullSpace[Mmat /. alpha -> alphaReq - a0 (1 - xi) IdentityMatrix[2]]` (or equivalent) to derive the eigenvector. Then `fExpected` and `gExpected` need not be hardcoded — they fall out of substituting the derived eigenvector into the overlap construction. Hardcoded targets `fStage18`, `gStage19` may remain (they reference upstream units' results), but the intermediate algebra path must differ from the SymPy script's.

**Verification:**
After the fix, the Mathematica script's structure differs visibly: it contains `Solve[Det[...] == 0, alpha]` and `NullSpace[...]` or analogous matrix-eigenvalue calls, derives `alphaReq` rather than positing it, and arrives at `fGeneral`/`gGeneral` without relying on `fExpected`/`gExpected` as the comparison target (the comparison becomes between the two engines' derived `F_general` forms rather than against a shared hardcoded closed form).

### F4 — hardcoded_result

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage040_generalized_selected_branch_sympy_audit.py:100-101`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage040_generalized_selected_branch_mathematica_audit.wl:80-84`

**What's wrong:**
The "Stage-18 F" and "Stage-19 G" closed forms are hardcoded in both scripts:

- SymPy lines 100-101:
  ```
  F_stage18 = sp.simplify((9 * delta + 11 * xi)**4 / (81 * (1 - xi) * (9 * delta**2 + 18 * delta * xi + 11 * xi**2)**2))
  G_stage19 = sp.simplify(9 * xi * (delta + xi) / (9 * delta + 11 * xi))
  ```
- Mathematica lines 80-84:
  ```
  fStage18 = FullSimplify[(9 delta + 11 xi)^4/(81 (1 - xi) (9 delta^2 + 18 delta xi + 11 xi^2)^2), ...]
  gStage19 = FullSimplify[9 xi (delta + xi)/(9 delta + 11 xi), ...]
  ```

These literal expressions are stated without any in-script derivation and without a comment citing the specific upstream script and line where these forms are verified. The checks `F_U(R_U=1) - F_stage18 == 0` and `G_U(R_U=1) - G_stage19 == 0` then verify the current script's specialization against a hardcoded snippet of unstated provenance.

**Why this matters:**
If the upstream Stage-18/Stage-19 scripts later change their canonical form (sign, factor, or simplification representative), this stage's hardcoded copies will silently drift out of sync. The check would still pass on its own terms but would no longer reflect "Stage-18 F" as that file currently defines it. Provenance comments protect against this drift.

**Required change:**
Add a comment immediately above the `F_stage18` and `G_stage19` definitions in both scripts naming the specific upstream script file and the line range where these closed forms are verified. The form is informational only (no code change to the expressions themselves). Example for sympy line 99:

```
# F_stage18 and G_stage19 reproduce the closed forms verified in
# scripts/<filename>.py (lines NN-MM). If those upstream forms change,
# update here.
```

The auditor cannot supply the exact upstream filename and lines because notes/ and other units' scripts are out of scope for this audit; Codex must locate the upstream Stage-18 and Stage-19 sympy scripts (by filename containing "stage018" or "stage19" / "stage_18" / similar) and cite them.

**Verification:**
After the fix, both scripts have a non-blank comment immediately preceding the `F_stage18` / `fStage18` and `G_stage19` / `gStage19` definitions, citing an upstream script path and line range. No assertion is added or removed; this is a provenance fix.

## Independent-derivation check (Mathematica)

The Mathematica script is **not** an independent derivation. As detailed in F3, it transliterates the SymPy script: same hardcoded `alpha_req` (Mathematica line 45 ↔ SymPy line 59), same hardcoded comparison targets `fExpected` (line 59-63 ↔ line 82), same `fStage18` literal (line 80-83 ↔ line 100), and same algebraic recipe in the same order. Neither script derives `alpha_req` from the eigenvalue equation; both posit it and verify algebraic simplifications about `r`. Mathematica's `Eigenvalues`, `Eigenvectors`, `Solve`, and `NullSpace` primitives are not used anywhere.

## Engine cross-check

The two engines produce numerically and symbolically equivalent outputs:

| Quantity | SymPy output | Mathematica output |
|---|---|---|
| `alpha_req` | `A0*xi*(delta+xi)/(z0^2*(delta+q^2*xi+xi))` | `(a0*xi*(delta+xi))/((delta+xi+q^2*xi)*z0^2)` |
| `e1/e0` | `q*xi/(delta+xi)` | `(q*xi)/(delta+xi)` |
| `F_(q,eta)` | `-(delta+eta*xi+xi)^2*(delta+q^2*xi+xi)^2/((xi-1)*(...)^2)` | `((delta+xi+eta*xi)^2*(delta+xi+q^2*xi)^2)/((1-xi)*(...)^2)` (same up to factoring (xi-1)= -(1-xi)) |
| `G_q` | `xi*(delta+xi)/(delta+q^2*xi+xi)` | `(xi*(delta+xi))/(delta+xi+q^2*xi)` |
| `H_F` | `4*xi*(27*delta^2 + 36*delta*xi + 11*xi^2)/((9*delta+11*xi)*(9*delta^2+18*delta*xi+11*xi^2))` | `(4*xi*(27*delta^2 + 36*delta*xi + 11*xi^2))/((9*delta+11*xi)*(9*delta^2+18*delta*xi+11*xi^2))` |
| `H_G` | `-4*xi/(9*delta+11*xi)` | `(-4*xi)/(9*delta+11*xi)` |

All seven `expect_zero`/`expectZero` checks pass (`= 0`) in both engines. Engines agree at the symbolic level. This agreement, however, is consistent with the F3 transliteration finding: identical inputs to identical algebra produce identical outputs, which is not the same as two independent derivations corroborating each other.

## Verdict justification

The scripts pass every assertion they make, but the assertions are weaker than the docstring's claims. The eigenvalue/eigenvector existence statement of section 23.1 is asserted in the docstring and used as the foundation for everything downstream, yet no script writes the perturbed matrix or checks the eigenvalue residual — the check at line 66 is a pure algebraic identity about posited expressions (F1). Section 23.4's "exact first-order expansion" check is tautological by construction: it verifies that `sp.series` and `sp.diff` agree on the same expression (F2). The Mathematica script is a line-by-line transliteration of the SymPy script rather than an independent re-derivation (F3); both engines posit `alpha_req` and `r` rather than deriving them from the matrix equation, so the cross-engine agreement provides limited corroboration. Stage-18/19 reference forms are hardcoded without provenance (F4, low severity). Attacks I tried that failed: (a) testing whether `lam_minus = A0(1-xi)` could actually be an eigenvalue under the rank-1-subtracted convention — it can, I verified at `q=1, delta=1, xi=1/2, A0=z0=1` and got a singular matrix `[[1/8, -3/8], [-3/8, 9/8]]` with null vector `(1, 1/3)` matching the script's `r`; (b) checking parity of `F_general, G_general` under `q -> -q` (both depend only on `q^2`, so the negative-q substitution for split-U is safe); (c) checking that the `removeO()` truncation in `sp.series(..., 2)` correctly produces the linear-order expansion (it does). Verdict: `findings` (not stop-cold) — the math holds where it's exercised, but the verification surface is too narrow and the second-engine independence is absent.

## Self-test notes

**Variable independence (F1 directive):** the proposed eigenvalue residual `(M_perturbed - lam_minus * I) @ e_minus` depends on `A0, delta, xi, q, z0` through every component. No `sp.diff(EXPR, VAR)` with `VAR` outside `EXPR`'s symbol set. Trivial-case check: at `q=1, delta=1, xi=1/2, A0=z0=1` I worked the arithmetic out by hand (`alpha_req = 3/8`, `r = 1/3`, residual = `(1/8, -3/8)*(1,1/3) - 1/2*(1,1/3) = (1/8 - 1/8, -3/8 + 3/8) - (1/2, 1/6) = (0, 0) - ...` — wait, recomputed: `M_perturbed - lam I = [[1/8, -3/8], [-3/8, 9/8]]`, applied to `(1, 1/3)` gives `(1/8 - 1/8, -3/8 + 3/8) = (0, 0)`. ✓ Residual is genuinely zero, so the proposed `expect_zero` will pass after the directive is applied. **Parity (F4 / split-U):** confirmed `F_general`, `G_general` are even in `q`, so the negative-q substitution for `q_U` is safe — not raised as a finding. **Tautology probe (F2):** confirmed by reading the script that `F_ratio` and `1 + eps*HF` are constructed from the same expression via series/diff, with `f(0) = 1` already established in section 23.3; so the assertion is structurally true and the finding is real. **Path specifications:** no `missing_verification_script` findings here, so no path-disambiguation traps. The F3 directive targets a file in `mathematica/` (correct directory for `.wl`). The F1 directive edits both `scripts/` (`.py`) and `mathematica/` (`.wl`) at the cited line ranges, both confirmed to exist in this audit.
