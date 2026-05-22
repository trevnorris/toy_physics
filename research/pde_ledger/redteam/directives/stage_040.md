---
unit_id: 040
batch: III.1
created_at: 2026-05-22T00:00:00-06:00
findings_count: 4
stop_cold: null
applied: true
applied_at: 2026-05-22T12:29:39-06:00
findings_applied: 4
findings_blocked: 0
verification_status: pending
---

# Codex directive — unit 040

Apply each finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

If a finding's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named.

Do NOT run python or mathematica. Only edit files.

Do NOT touch paper.tex, notes/, or any prose documents. The red-team only modifies scripts.

## F1 — insufficient_verification

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage040_generalized_selected_branch_sympy_audit.py:66` (insert new check immediately after line 66, before subbanner on line 68)
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage040_generalized_selected_branch_mathematica_audit.wl:50` (insert new check immediately after line 50, before subbanner on line 52)

**Issue:** The docstring claims a selected-branch eigenvalue and eigenvector formula, but neither script writes the perturbed 2x2 matrix or verifies `(M_perturbed - lam_minus * I) @ e_minus == 0`. The current assertion only validates an algebraic identity between two posited expressions. Add an explicit matrix-residual check so the eigenvalue claim is actually exercised.

**Required change:**

For the SymPy script, insert after line 66 (and before the blank line / subbanner on lines 67-68):

```python

# Explicit eigenvalue/eigenvector residual check.
# Convention (consistent with the sign chosen for r above): the perturbed
# matrix is M - alpha z z^T, where M = diag(A0, A0*(1+delta)) and
# z = (z0, q*z0)^T. The selected-branch claim is that lam_minus = A0*(1-xi)
# is an eigenvalue and (1, r)^T with r = q*xi/(delta+xi) is the eigenvector.
M_base = sp.Matrix([[A0, 0], [0, A0 * (1 + delta)]])
z_vec = sp.Matrix([z0, q * z0])
M_perturbed = M_base - alpha_req_q * (z_vec * z_vec.T)
e_minus = sp.Matrix([1, q * xi / (delta + xi)])
eig_residual = sp.simplify(M_perturbed * e_minus - lam_minus * e_minus)
expect_zero("eigenvector residual row 0", eig_residual[0])
expect_zero("eigenvector residual row 1", eig_residual[1])
```

For the Mathematica script, insert after line 50 (and before the blank line / subbanner on lines 51-52):

```mathematica

(* Explicit eigenvalue/eigenvector residual check. *)
(* Convention: perturbed matrix is M - alpha z z^T, with             *)
(* M = DiagonalMatrix[{a0, a0 (1+delta)}] and z = {z0, q z0}.         *)
mBase = DiagonalMatrix[{a0, a0 (1 + delta)}];
zVec = {z0, q z0};
mPerturbed = mBase - alphaReq Outer[Times, zVec, zVec];
eMinus = {1, q xi/(delta + xi)};
eigResidual = FullSimplify[mPerturbed.eMinus - lamMinus eMinus,
  Assumptions -> $Assumptions];
expectZero["eigenvector residual row 0", eigResidual[[1]]];
expectZero["eigenvector residual row 1", eigResidual[[2]]];
```

**Verification command:**
After Codex applies, the verifier will run `redteam exec-sympy 040` and `redteam exec-mathematica 040` and confirm both new `eigenvector residual` checks appear in the output and both scripts exit 0.

## Applied: F1

- files_changed:
  - `scripts/moving_throat_pde_stage040_generalized_selected_branch_sympy_audit.py`
  - `mathematica/moving_throat_pde_stage040_generalized_selected_branch_mathematica_audit.wl`
- summary: Added explicit perturbed-matrix eigenvector residual checks for the selected lower branch.
- deviation: none

## F2 — tautological_check

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage040_generalized_selected_branch_sympy_audit.py:106-125`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage040_generalized_selected_branch_mathematica_audit.wl:91-103`

**Issue:** The assertions `F_ratio - (1 + eps * HF) == 0` and `G_ratio - (1 + eps * HG) == 0` are tautological by construction — `F_ratio` is the linear Taylor expansion (via `sp.series`), and `HF` is the first derivative at `eps=0` (via `sp.diff`), so equality holds for any underlying expression. The assertions verify a CAS identity, not a physical claim. Replace with a substantive cross-check: derive `HF` and `HG` from an independent path (the section-23.2 overlap construction with `q -> q_U(1+eps)`, `eta -> eta_U(1+eps)`) and assert the two derivations agree.

**Required change:**

For the SymPy script, replace lines 106-125 (the entire 23.4 block from the subbanner through both `expect_zero` calls) with the following block. Keep the subbanner.

```python
subbanner("23.4 — Independent cross-check of first-order deformation about flat-U limit")

# Define R_U-dependent q and eta substitutions parametrized by eps,
# so that R_U = 1 + eps.
q_U_eps = sp.simplify(-sp.sqrt(lam0) * (1 + eps))
eta_U_eps = sp.simplify(lam0 * (1 + eps))

# Path 1 (script's original route): differentiate F_U with respect to R_U.
HF = sp.simplify(sp.diff(F_U.subs(R_U, 1 + eps), eps).subs(eps, 0) / F_stage18)
HG = sp.simplify(sp.diff(G_U.subs(R_U, 1 + eps), eps).subs(eps, 0) / G_stage19)

# Path 2 (independent): substitute the eps-parametrized (q, eta) into
# the general two-vector functions F_general, G_general (section 23.2),
# then expand to first order in eps and read off the coefficient.
F_general_eps = F_general.subs({q: q_U_eps, eta: eta_U_eps})
G_general_eps = G_general.subs({q: q_U_eps})
HF_direct = sp.simplify(sp.diff(F_general_eps, eps).subs(eps, 0) / F_stage18)
HG_direct = sp.simplify(sp.diff(G_general_eps, eps).subs(eps, 0) / G_stage19)

print("H_F (via F_U)        =", sp.factor(HF))
print("H_F (via F_general)  =", sp.factor(HF_direct))
print("H_G (via G_U)        =", sp.factor(HG))
print("H_G (via G_general)  =", sp.factor(HG_direct))

expect_zero("H_F cross-check (F_U vs F_general)", sp.simplify(HF - HF_direct))
expect_zero("H_G cross-check (G_U vs G_general)", sp.simplify(HG - HG_direct))
```

For the Mathematica script, replace lines 91-103 (the entire section 4 block from the subbanner through both `expectZero` calls) with the following. Keep the subbanner.

```mathematica
subbanner["4. Independent cross-check of first-order deformation about flat-U limit"];

qUEps = FullSimplify[-Sqrt[lambda0] (1 + eps), Assumptions -> $Assumptions];
etaUEps = FullSimplify[lambda0 (1 + eps), Assumptions -> $Assumptions];

hF = FullSimplify[(D[fU /. rU -> 1 + eps, eps] /. eps -> 0)/fStage18, Assumptions -> $Assumptions];
hG = FullSimplify[(D[gU /. rU -> 1 + eps, eps] /. eps -> 0)/gStage19, Assumptions -> $Assumptions];

fGeneralEps = fGeneral /. {q -> qUEps, eta -> etaUEps};
gGeneralEps = gGeneral /. q -> qUEps;
hFDirect = FullSimplify[(D[fGeneralEps, eps] /. eps -> 0)/fStage18, Assumptions -> $Assumptions];
hGDirect = FullSimplify[(D[gGeneralEps, eps] /. eps -> 0)/gStage19, Assumptions -> $Assumptions];

Print["H_F (via F_U)        = ", fmt[hF]];
Print["H_F (via F_general)  = ", fmt[hFDirect]];
Print["H_G (via G_U)        = ", fmt[hG]];
Print["H_G (via G_general)  = ", fmt[hGDirect]];

expectZero["H_F cross-check (F_U vs F_general)", hF - hFDirect];
expectZero["H_G cross-check (G_U vs G_general)", hG - hGDirect];
```

**Verification command:**
After Codex applies, the verifier will run `redteam exec-sympy 040` and `redteam exec-mathematica 040` and confirm the outputs contain `H_F cross-check` and `H_G cross-check` lines that both report `= 0` and PASS, and that the old `F_ratio - (1 + eps H_F)` / `G_ratio - (1 + eps H_G)` lines no longer appear.

## Applied: F2

- files_changed:
  - `scripts/moving_throat_pde_stage040_generalized_selected_branch_sympy_audit.py`
  - `mathematica/moving_throat_pde_stage040_generalized_selected_branch_mathematica_audit.wl`
- summary: Replaced the tautological first-order ratio checks with independent F/G deformation cross-checks through the general two-vector formulas.
- deviation: none

## F3 — mathematica_transliteration

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage040_generalized_selected_branch_mathematica_audit.wl:42-71`

**Issue:** The Mathematica script transliterates the SymPy algebra: same posited `alphaReq`, same hardcoded `fExpected`/`gExpected` comparison targets, same intermediate variable names. To satisfy the second-engine independence policy, the Mathematica script must derive `alphaReq` and the eigenvector from Mathematica primitives (`Solve`, `Det`, `NullSpace`, or `Eigenvectors`) on the explicit perturbed matrix, then construct `fGeneral`/`gGeneral` from that derivation rather than positing it.

**Required change:**

Replace lines 42-71 (the entire sections 1 and 2 blocks from the subbanner through `expectZero["G_general - expected", ...]`) with the following block. Keep the subbanners.

```mathematica
subbanner["1. Exact 2x2 selected-branch solve with generic loading ratio"];

(* Independent derivation: build the perturbed matrix, solve the         *)
(* characteristic equation det(M - alpha z z^T - lam I) = 0 for alpha    *)
(* at lam = a0 (1 - xi), then derive the eigenvector via NullSpace.     *)
mBase = DiagonalMatrix[{a0, a0 (1 + delta)}];
zVec = {z0, q z0};
mPert[alphaVal_] := mBase - alphaVal Outer[Times, zVec, zVec];
lamMinus = a0 (1 - xi);

charEq = Det[mPert[alpha] - lamMinus IdentityMatrix[2]] == 0;
alphaSol = Solve[charEq, alpha];
alphaReq = FullSimplify[alpha /. alphaSol[[1]], Assumptions -> $Assumptions];
Print["alpha_req (from Det = 0) = ", fmt[alphaReq]];

(* Eigenvector via NullSpace of (M - alpha_req z z^T - lam I).          *)
nsVec = NullSpace[mPert[alphaReq] - lamMinus IdentityMatrix[2]];
eMinusRaw = FullSimplify[nsVec[[1]], Assumptions -> $Assumptions];
(* Normalize so the first component is 1. *)
eMinus = FullSimplify[eMinusRaw/eMinusRaw[[1]], Assumptions -> $Assumptions];
r = eMinus[[2]];
Print["e1/e0 (from NullSpace) = ", fmt[r]];
expectZero["e1/e0 closed form", r - xi q/(delta + xi)];

subbanner["2. Exact overlap formulas and generalized F,G functions"];

(* Build z and s overlaps from the derived eigenvector, NOT from a      *)
(* posited r ratio. Use s = (1, eta/q) so that s1/s0 * q = eta.         *)
sVec = {1, eta/q};
zDotE = zVec.eMinus;
sDotE = sVec.eMinus;
normESq = eMinus.eMinus;
zOverlapSq = FullSimplify[zDotE^2/(z0^2 normESq), Assumptions -> $Assumptions];
sOverlapSq = FullSimplify[sDotE^2/(1^2 normESq), Assumptions -> $Assumptions];
fGeneral = FullSimplify[(a0/lamMinus) zOverlapSq sOverlapSq, Assumptions -> $Assumptions];
gGeneral = FullSimplify[(z0^2/a0) alphaReq, Assumptions -> $Assumptions];

Print["(z.e_-)^2 / z0^2 = ", fmt[zOverlapSq]];
Print["(s.e_-)^2 / s0^2 = ", fmt[sOverlapSq]];
Print["F_(q,eta) = ", fmt[fGeneral]];
Print["G_q = ", fmt[gGeneral]];

(* Cross-check against the SymPy script's claimed closed forms. *)
fExpected = FullSimplify[
  (delta + (1 + q^2) xi)^2 (delta + (1 + eta) xi)^2/
    ((1 - xi) ((delta + xi)^2 + q^2 xi^2)^2),
  Assumptions -> $Assumptions
];
gExpected = FullSimplify[xi (delta + xi)/(delta + (1 + q^2) xi), Assumptions -> $Assumptions];
expectZero["F_general - expected", fGeneral - fExpected];
expectZero["G_general - expected", gGeneral - gExpected];
```

Note: the existing `lamMinus = a0 (1 - xi)` definition on line 44 is folded into the new block. The `alpha` symbol must be cleared at the top of the file if it isn't already; verify line 35's `Clear[a0, delta, z0, q, eta, xi, rU, eps]` covers it — if `alpha` is not listed there, add it: change line 35 to `Clear[a0, delta, z0, q, eta, xi, rU, eps, alpha]`.

**Verification command:**
After Codex applies, the verifier will run `redteam exec-mathematica 040` and confirm: (a) the output contains `alpha_req (from Det = 0)`, (b) the output contains `e1/e0 (from NullSpace)`, (c) the existing `e1/e0 closed form`, `F_general - expected`, `G_general - expected` checks still pass, and (d) the script exits 0.

## Applied: F3

- files_changed:
  - `mathematica/moving_throat_pde_stage040_generalized_selected_branch_mathematica_audit.wl`
- summary: Reworked the Mathematica selected-branch block to derive alpha via Det/Solve and the eigenvector via NullSpace before constructing the overlap formulas.
- deviation: Stored/restored the NullSpace-derived eigenvector around the F1 residual check so the overlap construction remains derivation-based.

## F4 — hardcoded_result

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage040_generalized_selected_branch_sympy_audit.py:99-101`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage040_generalized_selected_branch_mathematica_audit.wl:79-84`

**Issue:** The Stage-18 `F` and Stage-19 `G` closed forms are hardcoded with no in-script provenance comment. If the upstream scripts later canonicalize these expressions differently, this stage's hardcoded copies will silently drift out of sync.

**Required change:**

Codex must first locate the upstream files that verify the Stage-18 and Stage-19 closed forms (search the `scripts/` directory for filenames containing `stage018`, `stage_018`, `stage19`, or `stage_019` and matching the moving-throat-pde topic). Once located, add a comment immediately above the `F_stage18` / `G_stage19` definitions in both scripts citing those files and the relevant line ranges where the closed forms are verified.

For the SymPy script, insert a comment block immediately before line 100 (i.e., between lines 99 and 100). Suggested form (Codex fills in `<upstream_file>` and `<line_range>`):

```python
# F_stage18 reproduces the Stage-18 closed-form normalization F(xi, delta)
# verified in scripts/<upstream_file>.py lines <line_range>.
# G_stage19 reproduces the Stage-19 closed-form loading G(xi, delta)
# verified in scripts/<upstream_file>.py lines <line_range>.
# Keep these literals in sync with the upstream source of truth.
```

For the Mathematica script, insert an analogous comment immediately before line 80 (i.e., between lines 79 and 80):

```mathematica
(* fStage18 reproduces the Stage-18 closed-form F(xi, delta) verified  *)
(* in mathematica/<upstream_file>.wl lines <line_range>.                *)
(* gStage19 reproduces the Stage-19 closed-form G(xi, delta) verified  *)
(* in mathematica/<upstream_file>.wl lines <line_range>.                *)
(* Keep these literals in sync with the upstream source of truth.      *)
```

If Codex cannot locate the upstream files unambiguously (e.g., zero matches or multiple equally plausible matches), append `## Blocked: F4` with a question describing what was found and requesting clarification — do not invent file names.

**Verification command:**
After Codex applies, the verifier opens both scripts and confirms a non-empty comment block citing concrete upstream file paths and line numbers appears immediately above the `F_stage18`/`fStage18` and `G_stage19`/`gStage19` definitions. No assertions are added or removed.

## Applied: F4

- files_changed:
  - `scripts/moving_throat_pde_stage040_generalized_selected_branch_sympy_audit.py`
  - `mathematica/moving_throat_pde_stage040_generalized_selected_branch_mathematica_audit.wl`
- summary: Added provenance comments above the Stage-18 F and Stage-19 G literals citing the upstream audit files and line ranges.
- deviation: none
