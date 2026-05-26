---
unit_id: 037
batch: III.1
created_at: 2026-05-26T00:00:00Z
findings_count: 1
stop_cold: null
applied: true
applied_at: 2026-05-26T07:35:25Z
findings_applied: 1
findings_blocked: 0
verification_status: pending
needs_user_resolution: false
---

# Codex directive - unit 037

Apply each non-paper_misalignment finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

If a finding's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead - skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named.

Do NOT run python or mathematica. Only edit files.

Do NOT touch paper.tex, notes/, or any prose documents. The red-team only modifies scripts.

## F1 - insufficient_verification

**Target:** `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage037_continuum_kernel_mathematica_audit.wl:114-127`

**Issue:**

The Mathematica audit's Schur-complement step (lines 114-125 in the current `.wl`) verifies only the *rank-one+identity shape* of `Sigma_wall`: it computes `sigmaWall = cMat . LinearSolve[bMat, Transpose[cMat]]`, solves backwards for `xi` and `alpha` from two entries of the computed Sigma, and then checks one consistency equation on the third (the `(2,2)` entry). It does NOT assert that the recovered `xiSolved` equals `g_U^2/A_U` or that `alphaSolved` equals the closed form in the notes:

```
Xi(omega) = g_U^2 / A_U(omega),
alpha(omega) = g_B^2 / A_phi(omega) + ( A_U(omega) g_W + g_R g_U )^2 / [ A_U(omega) Delta_UW(omega) ],
Delta_UW(omega) = A_U(omega) A_W(omega) - g_R^2 sigma.
```

(See `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage037_continuum_kernel_extraction.md` section 4 for the verbatim statement.)

The SymPy audit verifies these closed forms directly via a full 2x2 matrix comparison (`scripts/moving_throat_pde_stage037_continuum_kernel_sympy_audit.py` line 144). The Mathematica side must also verify the closed forms - the shape-only check is strictly weaker and would not catch a regression in either `Xi` or `alpha` that preserved the rank-one+identity structure.

**Required change:**

In `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage037_continuum_kernel_mathematica_audit.wl`, edit Section 3 only (Schur-complement factorization). After the existing `sigmaWall = FullSimplify[cMat . LinearSolve[bMat, Transpose[cMat]], Assumptions -> $Assumptions];` line (current line 114), and before the existing `alphaSolved = ...` line (current line 116), insert independent closed-form definitions and a full 2x2 matrix assertion. The new block must:

1. Define
```
deltaUW = aU*aW - gR^2 * sigma;
xiExpected = gU^2 / aU;
alphaExpected = gB^2/aPhi + (aU*gW + gR*gU)^2 / (aU * deltaUW);
sigmaWallExpected = FullSimplify[
  xiExpected * IdentityMatrix[2]
    + alphaExpected * {{kappa0^2, kappa0*kappa1}, {kappa0*kappa1, kappa1^2}},
  Assumptions -> $Assumptions
];
```

2. Assert the full 2x2 matrix residual vanishes:
```
expectMatrixZero["Sigma_wall - [Xi I + alpha v v^T]", sigmaWall - sigmaWallExpected];
```

3. Leave the existing solve-then-consistency choreography (current lines 116-127) in place. It is still useful: it provides a *second* route through the Schur form via backward solve, and the printed `xiSolved` / `alphaSolved` lines are diagnostically valuable. The new matrix-equality check is added in addition, not in replacement.

Concretely, the inserted lines come immediately after `sigmaWall = ...;` (current line 114) and before `v = {{kappa0}, {kappa1}};` (current line 115). The function `expectMatrixZero` is already defined at the top of the file (current lines 33-38) and accepts the same call pattern.

The edit affects only lines 114-115 (insertion); the existing assertion at the current line 125 (`expectZero["Sigma_wall (2,2) consistency with ansatz", ...]`) and surrounding code remain untouched.

**Claim manifest:** (not applicable - script exists; this is an added assertion, not a missing-script reconstruction)

**Verification command:**

After Codex applies, the verifier runs the Mathematica script and confirms:
1. The new assertion `Sigma_wall - [Xi I + alpha v v^T]` appears in `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage037_continuum_kernel_mathematica_audit.txt` with a `PASS:` line (residual `0`).
2. The script exits 0.
3. The existing `Sigma_wall (2,2) consistency with ansatz` PASS line remains present (the patch added a check; it did not replace one).
4. No other assertions or printed values change between the prior transcript and the new transcript (apart from the inserted PASS line).

## Applied: F1

- files_changed:
  - `mathematica/moving_throat_pde_stage037_continuum_kernel_mathematica_audit.wl`
- summary: Added the closed-form Schur complement matrix assertion before the existing recovered `Xi`/`alpha` consistency check.
- deviation: none
