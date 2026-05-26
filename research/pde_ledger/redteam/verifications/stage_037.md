---
unit_id: 037
batch: III.1
verifier_model: claude-opus-4-7[1m]
verify_date: 2026-05-26T00:00:00Z
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 1
findings_total: 1
material_change: false
---

# Verification — unit 037 (batch III.1 v2)

## Per-finding outcomes

### F1 — insufficient_verification

**Classification:** resolved

**What changed:**

In `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage037_continuum_kernel_mathematica_audit.wl:115-123`, Codex inserted nine new lines immediately after the existing `sigmaWall = FullSimplify[cMat . LinearSolve[bMat, Transpose[cMat]], Assumptions -> $Assumptions];` (line 114) and before the existing `v = {{kappa0}, {kappa1}};` (now line 124):

```
deltaUW = aU*aW - gR^2*sigma;
xiExpected = gU^2/aU;
alphaExpected = gB^2/aPhi + (aU*gW + gR*gU)^2/(aU*deltaUW);
sigmaWallExpected = FullSimplify[
  xiExpected*IdentityMatrix[2]
    + alphaExpected*{{kappa0^2, kappa0*kappa1}, {kappa0*kappa1, kappa1^2}},
  Assumptions -> $Assumptions
];
expectMatrixZero["Sigma_wall - [Xi I + alpha v v^T]", sigmaWall - sigmaWallExpected];
```

The pre-existing solve-then-consistency choreography is preserved unchanged: `v = {{kappa0}, {kappa1}};` (now line 124), the two `Solve` calls for `alphaSolved` (lines 125-128) and `xiSolved` (lines 129-132), the original `expectZero["Sigma_wall (2,2) consistency with ansatz", sigmaWall[[2, 2]] - xiSolved - alphaSolved*kappa1^2];` (now line 134), and the diagnostic `Print["Xi (recovered) = ...]` / `Print["alpha (recovered) = ...]` (lines 135-136). The diff at `redteam/exec_logs/stage_037_diff.patch` confirms the patch is a clean nine-line insertion with no other edits.

**Assessment:**

The edit matches the directive line-for-line. The four definitions are the literal closed forms from the notes (`Xi(omega) = g_U^2 / A_U`, `alpha(omega) = g_B^2 / A_phi + (A_U g_W + g_R g_U)^2 / (A_U Delta_UW)`, `Delta_UW = A_U A_W - g_R^2 sigma`), defined independently of the computed `sigmaWall = cMat . LinearSolve[bMat, Transpose[cMat]]`. The assertion uses the pre-existing `expectMatrixZero` helper (lines 33-38) which applies `FullSimplify[Map[Together[Expand[#]] &, expr, {2}], Assumptions -> $Assumptions]` and tests element-wise strict equality to a zero 2x2 matrix — so the residual is computed on all four entries of the 2x2 symmetric matrix, strictly stronger than the prior `(2,2)`-only consistency check.

Non-tautology check: `sigmaWall` is computed from `bMat`, `cMat` via `LinearSolve`, while `sigmaWallExpected` is built directly from `xiExpected`, `alphaExpected`, `kappa0`, `kappa1`. The two paths share only the input symbols (`aU`, `aW`, `aPhi`, `gU`, `gW`, `gB`, `gR`, `kappa0`, `kappa1`, `sigma`). A sign-flip or missing-coupling regression in either `Xi` or `alpha` would leave a nonzero rational residual that `FullSimplify` would not zero out — the new check has real teeth.

The Mathematica exec log captures the post-fix run with exit code 0; the new assertion prints `Sigma_wall - [Xi I + alpha v v^T] = {{0, 0}, {0, 0}}` followed by `PASS: Sigma_wall - [Xi I + alpha v v^T]` (log lines 54-55). The pre-existing `Sigma_wall (2,2) consistency with ansatz` PASS line is preserved (log lines 56-57). The recovered-`Xi`/`alpha` diagnostic prints (log lines 58-59) are unchanged. All Section 4 / Section 5 PASS lines are unchanged in content. No collateral edits, no deviation.

The new check substantively closes the gap F1 identified: the Mathematica engine now verifies the explicit closed forms for `Xi` and `alpha` from the notes (section 4), bringing it to parity with the SymPy assertion at `scripts/moving_throat_pde_stage037_continuum_kernel_sympy_audit.py:144`. The directive's success criteria are met (new PASS line present, residual zero, script exits 0, prior PASS line retained, no other assertions changed).

## Exec log assessment

**SymPy:** exit=0. The SymPy script was not edited under this directive. Notable lines (from `redteam/exec_logs/stage_037_sympy.log`):

```
Sigma_wall - [Xi I + alpha v v^T] =
⎡0  0⎤
⎣0  0⎦
A continuum formula = 0
delta continuum formula = 0
# exit_code: 0
```

All twelve assertions across Sections 1-5 pass; final exit 0.

**Mathematica:** exit=0. Notable lines (from `redteam/exec_logs/stage_037_mathematica.log`):

```
Sigma_wall - [Xi I + alpha v v^T] = {{0, 0}, {0, 0}}
PASS: Sigma_wall - [Xi I + alpha v v^T]
Sigma_wall (2,2) consistency with ansatz = 0
PASS: Sigma_wall (2,2) consistency with ansatz
Xi (recovered) = gU^2/aU
alpha (recovered) = (-88*aU*gB^2*gR^2 + 9*(aPhi*gR^2*gU^2 + aU*(aU*aW*gB^2 + 2*aPhi*gR*gU*gW + aPhi*aU*gW^2))*Pi^2)/(aPhi*aU*(-88*gR^2 + 9*aU*aW*Pi^2))
...
Stage 037 Mathematica audit passed.
# exit_code: 0
```

Both the new matrix-equality PASS and the pre-existing `(2,2)`-consistency PASS appear; the script terminates with the standard success banner and exit 0.

**Engine cross-check:** `Xi (recovered) = gU^2/aU` matches SymPy's `Xi = gU^2/A_U` symbolically. The `alpha (recovered)` rational form is the closed-form `alpha` from the notes after substituting `sigma = 88/(9*Pi^2)` and clearing common denominators. The new closed-form matrix-equality check in Mathematica directly mirrors the SymPy `expect_zero("Sigma_wall - [Xi I + alpha v v^T]", Sigma - Sigma_expected)` at `scripts/moving_throat_pde_stage037_continuum_kernel_sympy_audit.py:144`. Both engines now verify the same content.

**Output freshness:** The Mathematica `.wl` script mtime is 2026-05-26 01:35 (post-edit); the exec log timestamp header is 2026-05-26T01:49:53, after the edit, so the post-fix code was actually run. The saved transcript at `mathematica/output/moving_throat_pde_stage037_continuum_kernel_mathematica_audit.txt` still has mtime 2026-05-22 12:18 — older than the `.wl` script. The transcript on disk was NOT refreshed by the orchestrator; the canonical post-fix record lives in `redteam/exec_logs/stage_037_mathematica.log`, which is the artifact this prompt directs the verifier to consult. The SymPy script (`scripts/...sympy_audit.py`, mtime Apr 1 12:39) is unchanged; its saved output (May 22 12:17) is newer than the script and consistent with the unchanged exec log. Flag for the orchestrator: if the saved `.txt` is the downstream-published artifact, it should be regenerated so it matches the new exec log.

## Material-change assessment

`material_change`: false.

The patch adds a strictly stronger assertion in the Mathematica engine; it does not change any computed quantity, closed form, or numerical value. `A`, `DeltaK_ax`, `alpha_mix`, `beta_0`, `M_mix`, `delta`, `Xi`, `alpha`, `Delta_0`, `Chi` are unchanged in both engines. Downstream units (stage 038+) that consume these closed forms see identical values. The fix tightens the two-engine cross-check on the Schur factorization without altering any carry-forward data. No substantive reason to mark units > 037 stale.

## Side observations (non-blocking)

- The saved Mathematica transcript at `mathematica/output/moving_throat_pde_stage037_continuum_kernel_mathematica_audit.txt` was not regenerated post-fix (its mtime predates the `.wl` edit). The exec log captures the correct post-fix output, so verification is sound; this is a tooling concern for the orchestrator if downstream consumers read the saved `.txt` rather than the exec log.
- Stale internal banner: the `.wl` `banner["STAGE 020 — CONTINUUM-KERNEL EXTRACTION"]` (line 40) and SymPy banner `STAGE 20 — CONTINUUM-KERNEL EXTRACTION AUDIT` are holdovers from pre-renumbering. The original auditor explicitly classified this as cosmetic, not a finding. No action.
- Auditor's B4 observation (Mathematica `aDerived = FullSimplify[Together[k0 - gUCont^2/omegaU2]]` standalone identity, algebraically tautological as a standalone check but covered by the substantive B5 closed-form check) is unchanged by this patch. The auditor deliberately did not file it as a finding; the verifier accepts that classification and does not introduce a new one.
- The variable `v = {{kappa0}, {kappa1}}` at line 124 is defined but never used in the file after that point (the existing checks index `sigmaWall` directly). Harmless dead scaffolding; outside F1 scope.

## Verdict justification

Codex applied F1 exactly as the directive prescribed: a nine-line insertion in Section 3 of the `.wl` script defining independent closed-form values for `xiExpected`, `alphaExpected`, `deltaUW`, and `sigmaWallExpected`, then asserting the full 2x2 matrix residual `sigmaWall - sigmaWallExpected` is zero via the pre-existing `expectMatrixZero` helper. The pre-existing solve-then-consistency choreography remains in place as the directive required. The Mathematica exec log shows the new PASS line with residual `{{0, 0}, {0, 0}}`, the prior `(2,2)` consistency PASS line is preserved, all other assertions are unchanged, and the script exits 0. No collateral edits in the diff, no regressions, and the Mathematica engine now verifies the same closed-form content for `Xi` and `alpha` as the SymPy engine — closing the second-engine gap that F1 identified. Verdict: verified.
