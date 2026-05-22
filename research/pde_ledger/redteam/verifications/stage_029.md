---
unit_id: 029
batch: II.1
verifier_model: claude-opus-4-7[1m]
verify_date: 2026-05-22T00:00:00Z
verdict: verified
sympy_exit: n/a
mathematica_exit: n/a
findings_resolved: 3
findings_total: 3
material_change: false
---

# Verification — unit 029

## Per-finding outcomes

### F1 — tautological_check

**Classification:** resolved

**What changed:**
`mathematica/moving_throat_pde_stage029_dynamic_loading_mathematica_audit.wl:149-156` — the prior tautological residual `-tan2Theta - 2*alpha0*(-eta)/(DeltaKax + alpha0*xiConst)` (label "manifestly positive form") was replaced with a substantive stationarity-root check:

```
expectZero[
  "stationarity at theta_-",
  FullSimplify[
    ((DeltaKax + alpha0*xiConst)*Sin[2*theta] - 2*alpha0*eta*Cos[2*theta])
      /. theta -> ArcTan[2*alpha0*eta/(DeltaKax + alpha0*xiConst)]/2,
    Assumptions -> $Assumptions
  ]
];
```

The `tan2Theta` definition on line 142 and the `Print["tan(2 theta_-) = ", ...]` on line 157 are preserved per directive.

**Assessment:**
Correct and non-tautological. Substituting `theta -> ArcTan[2 a0 eta/(DKax + a0 xi)]/2` gives `2*theta = ArcTan[2 a0 eta/(DKax + a0 xi)]`. The residual then equals `(DKax + a0 xi)*sin(arctan(x)) - 2 a0 eta*cos(arctan(x))` where `x = 2 a0 eta/(DKax + a0 xi)`. Using `sin(arctan(x)) = x/sqrt(1+x^2)`, `cos(arctan(x)) = 1/sqrt(1+x^2)`, this becomes `(DKax + a0 xi)*x/sqrt(1+x^2) - 2 a0 eta/sqrt(1+x^2) = [2 a0 eta - 2 a0 eta]/sqrt(1+x^2) = 0`. Vanishing requires Mathematica to actually evaluate the arctan-sin/cos identity; it is not a `-x - (-x)` algebraic identity. The saved output (line 24-25) records `stationarity at theta_- = 0` followed by `PASS: stationarity at theta_-`, and the phrase "manifestly positive form" no longer appears anywhere in the output. The label change is honest: it now describes a true stationarity-root check rather than an identity rewriting.

### F2 — tautological_check

**Classification:** resolved

**What changed:**
- `scripts/moving_throat_pde_stage029_dynamic_loading_sympy_audit.py:152-154` — replaced the single comment `# Verify the isotropic shift cancels from the stiffness splitting.` with the prescribed three-line typo-guard comment. The `expect_zero("DeltaK_tilde - DeltaK_bare (Xi_0 cancellation)", DeltaK - DeltaK_ax)` call on line 155 is unchanged.
- `mathematica/moving_throat_pde_stage029_dynamic_loading_mathematica_audit.wl:130-131` — inserted the prescribed two-line typo-guard comment immediately before the existing `expectZero["DeltaK_tilde - DeltaK_ax", ...]` on line 132. The assertion line itself is unchanged.

**Assessment:**
Edits match the directive verbatim. Both assertions retained as typo guards (per directive instruction "Do NOT delete the assertion itself"). Output transcripts confirm both still PASS: SymPy output line 359 (`DeltaK_tilde - DeltaK_bare (Xi_0 cancellation) = 0`), Mathematica output lines 20-21 (`DeltaK_tilde - DeltaK_ax = 0` / `PASS`). Honest relabeling, no over-claim remains.

### F3 — mathematica_transliteration

**Classification:** resolved

**What changed:**
- `mathematica/moving_throat_pde_stage029_dynamic_loading_mathematica_audit.wl:87-116` — Schur block restructured. The one-shot `sigmaWall = Transpose[cMat].LinearSolve[mint, cMat]` was deleted; replaced by a sequential-elimination construction: `sigmaPhi = (lambdaB^2/aphi)*vv^T`, `uMassInv = Inverse[diag(aU,aU)] + (lambdaR^2/(aU*deltaUW))*vv^T` (the Schur-inverted U-block including the W back-reaction), `sigmaU = cU^T uMassInv cU`, and a separate `sigmaW` term capturing the direct lambdaW-vertex contribution through aW. The check is now `Sigma_seq - (Xi I + alpha vv^T) = 0`.
- `mathematica/moving_throat_pde_stage029_dynamic_loading_mathematica_audit.wl:194-226` — Hellmann–Feynman block restructured. Added a direct `Eigensystem[keffAl]` derivation that picks the lower eigenvector by `Position` against `(tr - Sqrt[disc])/2` and computes `kappaSelSqDirect = (vecMinusUnit.v)^2`. Both the closed-form `kappaSelSq = -D[lambdaMinusTemplate, al]` and `kappaSelSqDirect` are printed (output lines 43-44), and an explicit cross-engine check `expectZero["kappa_sel^2 closed-form vs eigenvector projection", kappaSelSq - kappaSelSqDirect]` is added. The weak/strong limit checks are repointed at `kappaSelSqDirect` so the limit tests exercise the independent derivation path.

**Codex deviation:** Codex's `## Applied: F3` block states "Adjusted the prescribed sigmaW split to avoid double-counting the u-block correction already present in sigmaU." The directive proposed `sigmaW = (aU*lambdaW + lambdaR*lambdaU)^2/(aU*deltaUW) * vv^T`; Codex used `sigmaW = (aU*lambdaW^2 + 2*lambdaR*lambdaU*lambdaW)/deltaUW * vv^T`. The difference is exactly `lambdaR^2*lambdaU^2/(aU*deltaUW)`, which is the `lambdaR^2/(aU*deltaUW)` term carried into `sigmaU` via the second piece of `uMassInv` (scaled by `lambdaU^2`). Without the deviation the directive's literal recipe would double-count the U-W mixing piece and the assertion would FAIL. Codex's correction is mathematically necessary, not a back-fit.

**Assessment:**
The sequential-elimination path is genuinely independent from the SymPy one-shot `C^T M^{-1} C`. It writes `Sigma_wall` as a sum of three rank-1 vv^T contributions plus an isotropic Xi*I piece, derived by serial elimination of phi -> W -> U with explicit accounting of the U-W Schur back-reaction. The algebraic recombination that lands on `Xi I + alpha vv^T` is non-trivial: coefficient of vv^T sums to `lambdaR^2 lambdaU^2/(aU deltaUW) + (aU lambdaW^2 + 2 lambdaR lambdaU lambdaW)/deltaUW + lambdaB^2/aphi`, which factors to `(aU lambdaW + lambdaR lambdaU)^2/(aU deltaUW) + lambdaB^2/aphi = alpha` (matches expected); coefficient of I sums to `lambdaU^2/aU = Xi`. A sign or transposition error in the 4x4 `mint` matrix would NOT propagate identically into the sequential path — the bug would have to be reproduced separately in `aphi`, `uMassInv`, or `sigmaW`. The second-engine policy is now genuinely satisfied.

The eigensystem path is also independent: `Eigensystem[keffAl]` constructs eigenvectors via row-reduction of `keffAl - lambda*I`, returning a different parametrisation than the SymPy `-d lambda_- / d alpha` template. The cross-check `kappaSelSq - kappaSelSqDirect = 0` (output line 45-46) is a real consistency test of two derivations rather than a transliteration; a bug in either would surface. The weak (`al -> 0 => kappa0^2`) and strong (`al -> Infinity => sigma`) limits are now exercised on the eigenvector-derived expression and both PASS (output lines 47-50).

No regressions in the diff: deleted code is the previous one-shot Schur formula and the previous limit-checks-on-closed-form; both are subsumed by the new independent paths plus cross-check. The constant-overlap PASS lines and downstream `kappa(theta)^2`, `delta D_-^(odd)` are unaffected.

## Exec log assessment

**SymPy:** exit=n/a. No `stage_029_sympy.log` exists in `redteam/exec_logs/`; only `stage_029_diff.patch` is present there. The saved transcript `scripts/output/moving_throat_pde_stage029_dynamic_loading_sympy_audit.txt` (mtime 17:13, vs the script's mtime 17:08) is post-fix. Notable lines from that saved output:
- L359: `DeltaK_tilde - DeltaK_bare (Xi_0 cancellation) = 0`
- L910: `weak-loading kappa_sel^2 -> kappa0^2 = 0`
- L911: `strong-loading kappa_sel^2 -> sigma = 0`
No FAIL lines anywhere. Final tail (not quoted in the grep) terminates with the standard pass banner.

**Mathematica:** exit=n/a. No `stage_029_mathematica.log` exists in `redteam/exec_logs/`. The saved transcript `mathematica/output/moving_throat_pde_stage029_dynamic_loading_mathematica_audit.txt` (mtime 17:14, vs the script's mtime 17:10) is post-fix. Notable lines:
- L9-10: `Sigma_seq - (Xi I + alpha vv^T) = {{0, 0}, {0, 0}}` / `PASS: Sigma_seq - (Xi I + alpha vv^T)`
- L24-25: `stationarity at theta_- = 0` / `PASS: stationarity at theta_-`
- L43-44: prints both `kappa_sel^2 (closed-form, Hellmann-Feynman) = ...` and `kappa_sel^2 (direct eigenvector projection) = ...` (their printed forms differ literally, confirming Codex did not simply alias them)
- L45-46: `kappa_sel^2 closed-form vs eigenvector projection = 0` / `PASS: kappa_sel^2 closed-form vs eigenvector projection`
- L53: `Stage 029 Mathematica audit passed.`
No FAIL anywhere, no "manifestly positive form" string anywhere.

**Output freshness:** Confirmed. SymPy output mtime 17:13 > SymPy script mtime 17:08. Mathematica output mtime 17:14 > Mathematica script mtime 17:10. Outputs are post-edit, so the PASS lines reflect the new substantive assertions.

## Material-change assessment

`material_change`: false.

The edits change the *form* of the Mathematica derivation (sequential vs one-shot Schur; Eigensystem vs `-D[lam_minus,al]`) but the final symbolic results (`Xi`, `alpha`, `kappa_sel^2`, `beta_5`, `alpha_crit`) printed in the transcripts are identical to the pre-fix values. The SymPy edit was a comment-only change. No downstream unit consumes the Mathematica intermediate `sigmaWall`; downstream stages consume the verified closed forms (`Xi`, `alpha`, `alpha_0`, `kappa_sel^2`, `beta_5`), which are unchanged. No regression risk.

## Side observations (non-blocking)

- Line 103-105 comment in the `.wl` file references a `cWphi` variable that was never instantiated in Codex's deviation (the directive's outline mentioned splitting the coupling into `cU` and `cWphi`, but Codex collapsed the W and phi contributions into the scalar coefficient on `sigmaW` rather than via a `cWphi` matrix). The comment is therefore mildly stale prose. Purely cosmetic, no functional impact.
- The `sigmaW` coefficient `(aU*lambdaW^2 + 2*lambdaR*lambdaU*lambdaW)/deltaUW` looks unusual at first glance (missing a `lambdaR^2*lambdaU^2/aU` term relative to the closed form), but as shown in the F3 assessment this is intentional and required for consistency with `uMassInv`. A future reviewer who reads the script without the directive context might want a one-line comment near line 108 saying "(missing lambdaR^2 lambdaU^2/aU piece is already inside uMassInv)". Not blocking.

## Verdict justification

All three findings are resolved with substance, not paper-overs. F1's replacement is a real arctan-sin/cos trig identity, not an algebraic rewrite. F2's edits are comment-only relabeling as the directive requested; the typo-guard assertion remains active in both engines. F3's restructuring delivers a genuinely independent Schur path (sequential elimination phi -> W -> U with explicit U-W Schur back-reaction) and a genuinely independent kappa_sel path (Eigensystem-based eigenvector projection), and crucially adds a cross-engine PASS that compares the two Mathematica derivations directly. Codex's `sigmaW` deviation from the directive's literal text was mathematically required to avoid double-counting and is correctly disclosed. Output transcripts are post-fix (mtime newer than script mtime in both engines), all PASS lines hold, no "manifestly positive form" residue, and no FAIL lines appear. No new findings noted that would justify needs_rework.
