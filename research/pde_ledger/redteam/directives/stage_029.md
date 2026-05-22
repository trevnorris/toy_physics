---
unit_id: 029
batch: II.1
created_at: 2026-05-21T00:00:00Z
findings_count: 3
stop_cold: null
applied: true
applied_at: 2026-05-21T23:11:15Z
findings_applied: 3
findings_blocked: 0
verification_status: pending
---

# Codex directive — unit 029

Apply each finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

If a finding's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named.

Do NOT run python or mathematica. Only edit files.

Do NOT touch paper.tex, notes/, or any prose documents. The red-team only modifies scripts.

## F1 — tautological_check

**Target:** `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage029_dynamic_loading_mathematica_audit.wl:120-123`

**Issue:** The current Mathematica check
```
expectZero[
  "-tan(2 theta_-) - manifestly positive form",
  -tan2Theta - 2*alpha0*(-eta)/(DeltaKax + alpha0*xiConst)
];
```
expands to `-x - (-x) = 0` after substituting `tan2Theta = 2*alpha0*eta/(DeltaKax + alpha0*xiConst)`. The residual is algebraically forced to vanish regardless of physics. The label claims a "manifestly positive form" but the check verifies an identity that holds for any expression. SymPy has no corresponding check, so this is a Mathematica-side defect.

**Required change:**

Replace the existing `expectZero[...]` block on lines 120-123 with a substantive stationarity-root check. The replacement should verify that the closed-form `tan(2 theta_-) = 2 alpha0 eta / (DeltaK_ax + alpha0 xi)` is indeed a root of the stationarity equation `(DeltaK_ax + alpha0 xi) Sin[2 theta] - 2 alpha0 eta Cos[2 theta] = 0`.

Before (lines 120-123):
```
expectZero[
  "-tan(2 theta_-) - manifestly positive form",
  -tan2Theta - 2*alpha0*(-eta)/(DeltaKax + alpha0*xiConst)
];
```

After:
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

Do NOT remove the `tan2Theta = FullSimplify[...]` definition on line 113 — it is still used by the `Print["tan(2 theta_-) = ", fmt[tan2Theta]]` on line 124.

Do NOT touch lines 113 (`tan2Theta` definition) or 124 (`Print["tan(2 theta_-) = ..."]`).

**Verification command:**
After Codex applies, the verifier will run `redteam exec-mathematica 029` and confirm the new check labeled `stationarity at theta_-` appears in the captured output with `PASS`, that the line containing the phrase `manifestly positive form` is gone, and that the script exits 0.

## Applied: F1

- files_changed:
  - `mathematica/moving_throat_pde_stage029_dynamic_loading_mathematica_audit.wl`
- summary: Replaced the tautological tan-form check with a stationarity-root check at theta_-.
- deviation: none

## F2 — tautological_check

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage029_dynamic_loading_sympy_audit.py:152`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage029_dynamic_loading_mathematica_audit.wl:102` (insert a new comment line immediately before the existing line 103 `expectZero[...]`)

**Issue:**
The check `(K1t - K0t) - DeltaK_ax = 0` is algebraically forced by the definitions on the two preceding lines (`K1 = K0 + DeltaK_ax`, `K0t = K0 - Xi0`, `K1t = K1 - Xi0`). The current comment in SymPy (`# Verify the isotropic shift cancels from the stiffness splitting.`) overclaims what the check accomplishes. The assertion has nonzero value as a typo guard (catches e.g. `K1t = K1 + Xi0`) but should be honestly labeled.

**Required change:**

(1) In SymPy at line 152, replace the single comment line:
```python
    # Verify the isotropic shift cancels from the stiffness splitting.
```
with the following three-line comment block:
```python
    # Sanity check: K0t = K0 - Xi0, K1t = K1 - Xi0, K1 = K0 + DeltaK_ax, so
    # (K1t - K0t) - DeltaK_ax = 0 is algebraically forced by the construction.
    # Kept here as a typo guard, not as a physics check.
```
Do NOT touch line 153 (the `expect_zero(...)` call itself); keep the assertion.

(2) In Mathematica, insert the following two-line comment block immediately before the existing `expectZero["DeltaK_tilde - DeltaK_ax", ...]` on line 103 (i.e. between current line 102 (which prints `alpha_0`) and current line 103):
```
(* Sanity check: K0t = K0 - Xi0, K1t = K1 - Xi0, K1 = K0 + DeltaKax, so
   (K1t - K0t) - DeltaKax = 0 is algebraically forced. Kept as typo guard. *)
```
Do NOT touch the existing `expectZero[...]` line itself; keep the assertion.

**Verification command:**
After Codex applies, the verifier will run both `redteam exec-sympy 029` and `redteam exec-mathematica 029`. Both should still exit 0, and both should still report PASS on `DeltaK_tilde - DeltaK_(ax|bare)`. Codex's `## Applied: F2` block should confirm the comments were edited but the assertions were not changed.

## Applied: F2

- files_changed:
  - `scripts/moving_throat_pde_stage029_dynamic_loading_sympy_audit.py`
  - `mathematica/moving_throat_pde_stage029_dynamic_loading_mathematica_audit.wl`
- summary: Relabeled the DeltaK cancellation checks as algebraic typo guards while leaving both assertions unchanged.
- deviation: none

## F3 — mathematica_transliteration

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage029_dynamic_loading_mathematica_audit.wl:74-89` (Schur block)
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage029_dynamic_loading_mathematica_audit.wl:161-168` (Hellmann–Feynman block)

**Issue:**
The Mathematica script reuses SymPy's exact derivation path: same `Mint` matrix layout with internal-field ordering `(u0, u1, W, phi)`, same coupling matrix `C` with the same row order, same one-shot Schur formula `C^T M^{-1} C`, same Hellmann–Feynman template `-D[lam_minus, al]`, same two limit checks. Surface syntax differs but the algebraic content is line-by-line identical. This violates the second-engine policy: a sign or transposition error in `Mint` or `C` would be reproduced verbatim in both engines and would PASS in both.

The fix restructures the Mathematica derivation to use a genuinely independent path. The same SymPy script remains unchanged — only Mathematica is restructured. The final assertions (which compare the independently-derived results to the same `Xi I + alpha vv^T` target and the same SymPy-computed `-d lambda_- / d al` template) become a real cross-engine check rather than a transliteration confirmation.

**Required change:**

(1) Schur block, lines 74-89. Keep the `mint`, `cMat` definitions on lines 74-85 for reference (they are also used by the printed output on line 87 to keep the diagnostic message). REPLACE the body of the Schur computation on lines 87-89 with a sequential-elimination derivation.

Insert immediately after line 85 (which ends `cMat = {...};`):

```
(* Independent derivation: eliminate internal fields one at a time. *)
(* Step 1: integrate out phi. phi couples to q via lambdaB v; Aphi phi = lambdaB v.q
   gives sigmaPhi = lambdaB^2/Aphi * Outer[Times, v, v]. *)
sigmaPhi = FullSimplify[(lambdaB^2/aphi)*Outer[Times, v, v], Assumptions -> $Assumptions];

(* Step 2: integrate out W. After eliminating phi, the W equation is
   aW*W = lambdaR (kappa0 u0 + kappa1 u1) + lambdaW v.q.
   Solving for W and substituting back contributes (per the U-W coupling) to both
   the diagonal Xi piece and the vv^T alpha piece. We compute the contribution
   sigmaW from the wall-q -> W -> (u-block, wall-q) closed loop. *)
uMassInv = FullSimplify[
  Inverse[{{aU, 0}, {0, aU}}] +
    (lambdaR^2/(aU*(aU*aW - lambdaR^2*sigma)))*Outer[Times, v, v],
  Assumptions -> $Assumptions
];

(* Step 3: contract C against the inverted internal block reconstructed from
   the sequential elimination. cU is the (u0,u1) part of the coupling; cWphi the
   (W,phi) part. *)
cU = {{lambdaU, 0}, {0, lambdaU}};
sigmaU = FullSimplify[Transpose[cU].uMassInv.cU, Assumptions -> $Assumptions];
sigmaW = FullSimplify[
  ((aU*lambdaW + lambdaR*lambdaU)^2 / (aU*(aU*aW - lambdaR^2*sigma)))*
    Outer[Times, v, v],
  Assumptions -> $Assumptions
];
sigmaWallSeq = FullSimplify[sigmaU + sigmaW + sigmaPhi, Assumptions -> $Assumptions];

sigmaExpected = FullSimplify[xiShift*i2 + alpha*Outer[Times, v, v], Assumptions -> $Assumptions];
expectMatrixZero["Sigma_seq - (Xi I + alpha vv^T)", sigmaWallSeq - sigmaExpected];
```

DELETE the original lines 87-89:
```
sigmaWall = FullSimplify[Transpose[cMat].LinearSolve[mint, cMat], Assumptions -> $Assumptions];
sigmaExpected = FullSimplify[xiShift*i2 + alpha*Outer[Times, v, v], Assumptions -> $Assumptions];
expectMatrixZero["Sigma - (Xi I + alpha vv^T)", sigmaWall - sigmaExpected];
```

The constant-overlap checks on lines 90-92 (`sigma - 88/(9 Pi^2)`, etc.) are independent of the Schur block and are kept unchanged.

(2) Hellmann–Feynman block, lines 161-168. REPLACE with a derivation that uses `Eigensystem` on the symbolic 2×2 `keff0` directly, rather than the closed-form `(tr - Sqrt[disc])/2` template.

Insert immediately before line 161:

```
(* Independent derivation: diagonalise the al-dependent effective stiffness
   matrix directly, take the eigenvalue with the smaller real part, and read
   off kappa_sel^2 = |projection of v onto the lower eigenvector|^2. *)
keffAl = {
  {k0t - al*kappa0Sq, -al*kappa0*kappa1},
  {-al*kappa0*kappa1, k1t - al*kappa1Sq}
};
{eigvals, eigvecs} = Eigensystem[keffAl];
(* The lower eigenvalue (with the - sign on Sqrt[disc]) corresponds to
   tr/2 - Sqrt[(tr/2)^2 - det]. Eigensystem orders eigenvalues by Mathematica's
   internal heuristic; we select by demanding the (tr - lam) / 2 piece matches. *)
lowerIdx = First[
  Position[Simplify[eigvals - ((k0t + k1t - al*sigma)/2 - Sqrt[(DeltaKax + al*xiConst)^2 + 4*al^2*eta^2]/2)],
    0, Infinity, Heads -> False]
];
lamMinusDirect = eigvals[[First[lowerIdx]]];
vecMinusDirect = eigvecs[[First[lowerIdx]]];
vecMinusUnit = vecMinusDirect/Sqrt[vecMinusDirect.vecMinusDirect];
kappaSelSqDirect = FullSimplify[(vecMinusUnit.v)^2, Assumptions -> $Assumptions];
```

Then REPLACE the existing lines 161-168:

```
discTemplate = FullSimplify[(DeltaKax + al*xiConst)^2 + 4*al^2*eta^2, Assumptions -> $Assumptions];
trTemplate = FullSimplify[k0t + k1t - al*sigma, Assumptions -> $Assumptions];
lambdaMinusTemplate = FullSimplify[(trTemplate - Sqrt[discTemplate])/2, Assumptions -> $Assumptions];
kappaSelSq = FullSimplify[-D[lambdaMinusTemplate, al], Assumptions -> $Assumptions];

Print["kappa_sel^2 = ", fmt[kappaSelSq]];
expectZero["weak-loading kappa_sel^2 - kappa0^2", (kappaSelSq /. al -> 0) - kappa0Sq];
expectZero["strong-loading kappa_sel^2 - sigma", FullSimplify[Limit[kappaSelSq, al -> Infinity], Assumptions -> DeltaKax > 0] - sigma];
```

with:

```
discTemplate = FullSimplify[(DeltaKax + al*xiConst)^2 + 4*al^2*eta^2, Assumptions -> $Assumptions];
trTemplate = FullSimplify[k0t + k1t - al*sigma, Assumptions -> $Assumptions];
lambdaMinusTemplate = FullSimplify[(trTemplate - Sqrt[discTemplate])/2, Assumptions -> $Assumptions];
kappaSelSq = FullSimplify[-D[lambdaMinusTemplate, al], Assumptions -> $Assumptions];

Print["kappa_sel^2 (closed-form, Hellmann-Feynman) = ", fmt[kappaSelSq]];
Print["kappa_sel^2 (direct eigenvector projection)  = ", fmt[kappaSelSqDirect]];

(* Cross-engine check: the two independent derivations must agree. *)
expectZero["kappa_sel^2 closed-form vs eigenvector projection", kappaSelSq - kappaSelSqDirect];

expectZero["weak-loading kappa_sel^2 - kappa0^2", (kappaSelSqDirect /. al -> 0) - kappa0Sq];
expectZero["strong-loading kappa_sel^2 - sigma", FullSimplify[Limit[kappaSelSqDirect, al -> Infinity], Assumptions -> DeltaKax > 0] - sigma];
```

Note: the weak and strong limit checks are repointed at the new `kappaSelSqDirect` so that the limit test exercises the independent derivation path. The closed-form `kappaSelSq` is retained and cross-checked against `kappaSelSqDirect`.

If `Eigensystem[keffAl]`'s symbolic output ordering is non-deterministic enough that `lowerIdx` cannot be unambiguously selected (e.g. the `Position[...]` match fails to find exactly one index), STOP and write a `## Blocked: F3` block asking the auditor for an alternative selection rule (e.g. `Sort[eigvals, Less /. _Symbol -> 0][[1]]`). Do not attempt to guess.

**Verification command:**
After Codex applies, the verifier will run `redteam exec-mathematica 029` and confirm:
- the captured output contains the new line `Sigma_seq - (Xi I + alpha vv^T) = {{0, 0}, {0, 0}}` followed by a `PASS: Sigma_seq - ...` line;
- the captured output contains both `kappa_sel^2 (closed-form, Hellmann-Feynman) = ...` AND `kappa_sel^2 (direct eigenvector projection) = ...`;
- the captured output contains a PASS line for `kappa_sel^2 closed-form vs eigenvector projection`;
- the script exits 0.

## Applied: F3

- files_changed:
  - `mathematica/moving_throat_pde_stage029_dynamic_loading_mathematica_audit.wl`
- summary: Replaced the transliterated Mathematica Schur and selected-mode checks with sequential-elimination and direct-eigenvector derivations.
- deviation: Adjusted the prescribed sigmaW split to avoid double-counting the u-block correction already present in sigmaU.
