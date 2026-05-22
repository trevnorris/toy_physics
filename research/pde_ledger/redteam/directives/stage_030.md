---
unit_id: 030
batch: II.1
created_at: 2026-05-21T00:00:00Z
findings_count: 3
stop_cold: null
applied: true
applied_at: 2026-05-21T17:15:49-06:00
findings_applied: 3
findings_blocked: 0
verification_status: pending
---

# Codex directive — unit 030

Apply each finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

If a finding's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named.

Do NOT run python or mathematica. Only edit files.

Do NOT touch paper.tex, notes/, or any prose documents. The red-team only modifies scripts.

## F1 — tautological_check

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage030_selected_mode_normalization_sympy_audit.py:113`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage030_selected_mode_normalization_mathematica_audit.wl:83`

**Issue:**
The assertion `Gamma5_sel - G5*P0_sel == 0` is forced to zero by the three preceding definition lines (`beta5 = G5*beta0`, `C5_sel = beta5*s`, `Gamma5_sel = C5_sel/lam_-`, `P0_sel = beta0*s/lam_-`). The check cannot fail no matter what symbolic inputs are chosen — it does not verify any physical claim, only the algebraic consequence of the four definitions written immediately above it. Replace it with a check that derives `Gamma5_sel` by an independent route — namely, the generic Part-I series expansion result `Im[Y(omega)]/omega^5 |_{omega -> 0} = C5/D0` instantiated at `D0 = lam_-` and `C5 = C5_sel`.

**Required change:**

The relation `Gamma5_sel = G5 * P0_sel` is a definitional consequence of how `Gamma5_sel`, `C5_sel`, and `P0_sel` are written in this script (specifically: `C5_sel := G5*beta0*s`, `Gamma5_sel := C5_sel/lam_-`, `P0_sel := beta0*s/lam_-`, so `Gamma5_sel = G5*beta0*s/lam_- = G5 * P0_sel` follows by substitution). There is no separate physical derivation visible in this script that the assertion could falsify. The honest fix is to drop the assertion and add an explanatory comment so the script does not claim to verify a relation it cannot test.

*SymPy* (`scripts/moving_throat_pde_stage030_selected_mode_normalization_sympy_audit.py`)

Replace the single line at `113`:

```python
expect_zero("Gamma5_sel - G5*P0_sel", Gamma5_sel - G5 * P0_sel)
```

with the explanatory comment (insert in place of line 113, leaving lines 100-112 and 114-117 unchanged):

```python
# Note: Gamma5_sel - G5*P0_sel == 0 follows by construction from the
# definitions on lines 101-104 (C5_sel := G5*beta0*s, Gamma5_sel := C5_sel/lam_-,
# P0_sel := beta0*s/lam_-), so it is not verified separately here. The physical
# content of this relation is the identification Gamma5 = C5/D0 at the selected
# mode, which is verified in generic form by the Part-I expansion (line 55).
```

*Mathematica* (`mathematica/moving_throat_pde_stage030_selected_mode_normalization_mathematica_audit.wl`)

Replace the single line at `83`:

```mathematica
expectZero["Gamma5_sel - G5*P0_sel", gamma5Sel - g5*p0Sel];
```

with:

```mathematica
(* Note: gamma5Sel - g5*p0Sel == 0 follows by construction from lines 74-77
   (c5Sel := g5*beta0*sMinusClosed, gamma5Sel := c5Sel/lamMinus,
   p0Sel := beta0*sMinusClosed/lamMinus). Not verified separately. The
   physical content (Gamma5 = C5/D0 at the selected mode) is verified in
   generic form by Part I (line 41). *)
```

**Verification command:**
After Codex applies, the verifier runs `redteam exec-sympy 030` then `redteam exec-mathematica 030`. The output transcripts must no longer contain the lines `"Gamma5_sel - G5*P0_sel = 0"` / `PASS: Gamma5_sel - G5*P0_sel`. Both scripts must still exit 0. The verifier confirms the removal by checking that the assertion text no longer appears in the captured output.

## Applied: F1

- files_changed:
  - `scripts/moving_throat_pde_stage030_selected_mode_normalization_sympy_audit.py`
  - `mathematica/moving_throat_pde_stage030_selected_mode_normalization_mathematica_audit.wl`
- summary: Removed the tautological Gamma5/P0 assertions and replaced them with the directive's explanatory comments.
- deviation: none

## F2 — tautological_check

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage030_selected_mode_normalization_sympy_audit.py:143-146`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage030_selected_mode_normalization_mathematica_audit.wl:103-109`

**Issue:**
The "spectral condition rewrite" assertion expands to `0 = 0` by direct substitution of `P0_sel = beta0*s/lam_-` and `lambda_req = mhat^2*beta0*s/NQ_target`. No spectral structure is exercised. The check cannot fail for any choice of `lam_-`, `s`, `mhat`, `beta0`, `NQ_target`. Replace it with a non-tautological check that exercises the actual spectral condition: at the physical fixed point `mhat^2*P0_sel = NQ_target`, the equality `lam_- = lambda_req` should hold as an identity. Implement this by solving `mhat^2*P0_sel = NQ_target` for `mhat^2` and substituting into `lambda_req`, then asserting `lambda_req` reduces to `lam_-`.

**Required change:**

Like F1, the "spectral condition rewrite" identity follows mechanically from the literal definitions `lambda_req := mhat^2*beta0*s/NQ_target` and `P0_sel := beta0*s/lam_-`. Once these are substituted, the asserted residual collapses to `0` regardless of the physical content of `lam_-`. The honest fix is to drop the assertion and add an explanatory comment.

*SymPy* (`scripts/moving_throat_pde_stage030_selected_mode_normalization_sympy_audit.py`)

Replace lines `143-146`:

```python
lambda_req = sp.simplify(mhat**2 * beta0 * s_minus_closed / NQ_target)
expect_zero("spectral condition rewrite", sp.simplify(lam_minus - lambda_req) + sp.simplify((mhat**2 * P0_sel - NQ_target) * lam_minus / NQ_target))
print("lambda_req =")
sp.pprint(lambda_req)
```

with:

```python
lambda_req = sp.simplify(mhat**2 * beta0 * s_minus_closed / NQ_target)
# Note: the algebraic rearrangement
#   (lam_- - lambda_req) + (mhat^2*P0_sel - NQ_target)*lam_-/NQ_target = 0
# follows by substituting the definitions P0_sel := beta0*s/lam_- and
# lambda_req := mhat^2*beta0*s/NQ_target; both sides collapse to 0 without
# exercising the physical content of lam_- or s. It is therefore not
# verified separately. The substantive content (that lam_- = lambda_req
# when the constraint mhat^2 P0_sel = NQ_target holds) is recorded here
# as a definitional identity, not as a check.
print("lambda_req =")
sp.pprint(lambda_req)
```

*Mathematica* (`mathematica/moving_throat_pde_stage030_selected_mode_normalization_mathematica_audit.wl`)

Replace lines `103-109`:

```mathematica
lambdaReq = FullSimplify[mhat^2*beta0*sMinusClosed/nQTarget, Assumptions -> $Assumptions];
expectZero[
  "spectral condition rewrite",
  FullSimplify[lamMinus - lambdaReq, Assumptions -> $Assumptions] +
  FullSimplify[(mhat^2*p0Sel - nQTarget)*lamMinus/nQTarget, Assumptions -> $Assumptions]
];
Print["lambda_req = ", fmt[lambdaReq]];
```

with:

```mathematica
lambdaReq = FullSimplify[mhat^2*beta0*sMinusClosed/nQTarget, Assumptions -> $Assumptions];
(* Note: (lamMinus - lambdaReq) + (mhat^2 p0Sel - nQTarget) lamMinus / nQTarget
   == 0 follows from the definitions p0Sel := beta0*sMinusClosed/lamMinus and
   lambdaReq := mhat^2*beta0*sMinusClosed/nQTarget by substitution, with no
   physical content of lamMinus required. Not verified separately. *)
Print["lambda_req = ", fmt[lambdaReq]];
```

**Verification command:**
After Codex applies, `redteam exec-sympy 030` and `redteam exec-mathematica 030`. The output transcripts must no longer contain `"spectral condition rewrite = 0"` / `PASS: spectral condition rewrite`. Both scripts must still exit 0. The verifier confirms the removal by checking that the assertion text no longer appears in the captured output, and that the explanatory comment is present in the source.

## Applied: F2

- files_changed:
  - `scripts/moving_throat_pde_stage030_selected_mode_normalization_sympy_audit.py`
  - `mathematica/moving_throat_pde_stage030_selected_mode_normalization_mathematica_audit.wl`
- summary: Removed the tautological spectral condition rewrite checks and kept lambda_req printing with explanatory comments.
- deviation: none

## F3 — mathematica_transliteration

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage030_selected_mode_normalization_mathematica_audit.wl:43-66` (Part II eigenvalue/overlap block)

**Issue:**
The Mathematica script is a line-by-line port of the SymPy script with cosmetic renamings (uppercase Greek to lowercase, `_` to camelCase). Part I, Part II, Part III, and Part IV all retrace the same variable choreography and intermediate expressions as the SymPy file. To satisfy the second-engine policy without rewriting every part, restructure Part II so that `lamMinus` and `lamPlus` are produced by Mathematica's `Eigenvalues[]` operating on an explicit 2x2 matrix `M`, rather than by typing the closed-form expressions verbatim. The matrix must reproduce the trace `2 a + dK - alpha*(x0+x1)` and the determinant `a*(a+dK) - alpha*((a+dK)*x0 + a*x1)` that the existing det identity (wl:87) checks against. The matrix below satisfies both:

```mathematica
mMat = {{a + dK - alpha*x1, -alpha*Sqrt[x0*x1]},
        {-alpha*Sqrt[x0*x1], a - alpha*x0}};
```

Verification of the matrix (already done by the auditor; included here so Codex does not need to re-derive):

- `Tr[mMat] = (a + dK - alpha*x1) + (a - alpha*x0) = 2 a + dK - alpha*(x0+x1)`. ✓
- `Det[mMat] = (a + dK - alpha*x1)*(a - alpha*x0) - alpha^2*x0*x1 = a^2 + a*dK - alpha*((a+dK)*x0 + a*x1)*alpha + alpha^2*x0*x1 - alpha^2*x0*x1 = a*(a+dK) - alpha*((a+dK)*x0 + a*x1)`. ✓ (matches `t0 = (a+dK)*x0 + a*x1` on existing wl:86.)
- Therefore `Eigenvalues[mMat] = {(Tr ± Sqrt[Tr^2 - 4 Det])/2}`, and the discriminant `Tr^2 - 4 Det = dK^2 + 2 alpha dK (x0-x1) + alpha^2 (x0+x1)^2 = (dK + alpha (x0-x1))^2 + 4 alpha^2 x0 x1`. ✓ matches the existing `r` in wl:51.

The loading direction `v` is NOT defined in the existing script, so the auditor cannot mechanically derive `sMinusClosed` from an eigenvector overlap. Keep the existing closed-form expression for `sMinusClosed` and the HF cross-check; only `lamMinus`/`lamPlus` are re-sourced from `Eigenvalues[mMat]`.

**Required change:**

In `mathematica/moving_throat_pde_stage030_selected_mode_normalization_mathematica_audit.wl`, replace lines `43-66` exactly with the block below. Do not change Parts I, III, IV (lines 26-41 and 68 onward).

```mathematica
banner["PART II — EXACT SELECTED LOWER EIGENVALUE AND OVERLAP"];

Clear[a, dK, alpha, x0, x1];
$Assumptions = Element[{a, dK, alpha, x0, x1}, Reals] && a > 0 && dK > 0 && alpha >= 0 && x0 > 0 && x1 > 0;

(* Independent derivation: define the 2x2 wall block explicitly and let
   Mathematica's Eigenvalues[] produce the spectrum. The matrix is chosen so
   that Tr[mMat] = 2 a + dK - alpha (x0+x1) and
   Det[mMat] = a (a+dK) - alpha ((a+dK) x0 + a x1) -- the same trace/det as the
   determinant identity verified later in this script (see wl:87). This breaks
   the line-by-line transliteration of the SymPy file by routing lamMinus /
   lamPlus through Eigenvalues[mMat] rather than through a typed closed form. *)
mMat = {{a + dK - alpha*x1, -alpha*Sqrt[x0*x1]},
        {-alpha*Sqrt[x0*x1], a - alpha*x0}};

sigma = x0 + x1;
deltaKappa = x0 - x1;
kappaProd = x0*x1;
r = Sqrt[(dK + alpha*deltaKappa)^2 + 4*alpha^2*kappaProd];

eigVals = Eigenvalues[mMat];
lamMinus = FullSimplify[
  First[Select[eigVals,
    FullSimplify[# - ((2*a + dK - alpha*sigma) - r)/2, Assumptions -> $Assumptions] === 0 &]],
  Assumptions -> $Assumptions
];
lamPlus = FullSimplify[
  First[Select[eigVals,
    FullSimplify[# - ((2*a + dK - alpha*sigma) + r)/2, Assumptions -> $Assumptions] === 0 &]],
  Assumptions -> $Assumptions
];

Print["lambda_- = ", fmt[lamMinus]];
Print["lambda_+ = ", fmt[lamPlus]];

sMinusHF = FullSimplify[-D[lamMinus, alpha], Assumptions -> $Assumptions];
sMinusClosed = FullSimplify[
  (sigma + ((dK + alpha*deltaKappa)*deltaKappa + 4*alpha*kappaProd)/r)/2,
  Assumptions -> $Assumptions
];
expectZero["selected overlap: HF - closed form", sMinusHF - sMinusClosed];
Print["s_- = (v.e_-)^2 = ", fmt[sMinusClosed]];
expectZero["weak-loading overlap limit", (sMinusClosed /. alpha -> 0) - x0];
```

Notes for Codex:

- The `Select` + `First` idiom is used because Mathematica's `Eigenvalues[]` may return the eigenvalues in either order; we select the one matching the minus-root closed-form to bind to `lamMinus`. If `Select` returns an empty list (which would mean the matrix is wrong), the script will error out at `First[]`, which is the correct failure behavior.
- The `r`, `sigma`, etc. assignments are kept because they are used by `sMinusClosed`; do not delete them.
- Do not change `t0` (line 86) or the `det identity` assertion (line 87) — those work unchanged with the new `lamMinus`/`lamPlus`.

**Verification command:**
After Codex applies, the verifier runs `redteam exec-mathematica 030`. Expected outcome:
1. Output transcript still prints `lambda_- = ...` and `lambda_+ = ...` with the same canonical forms as before (Mathematica's `FullSimplify` of `Eigenvalues[mMat]` reduces to the same expression).
2. All four Part II assertions (`selected overlap: HF - closed form`, `weak-loading overlap limit`) plus the downstream Part III det identity continue to pass with residual 0.
3. The script exits 0.
4. The source file no longer contains `(2*a + dK - alpha*sigma - r)/2` typed as an *assignment* to `lamMinus` / `lamPlus`; instead, `Eigenvalues[mMat]` is the source.

## Applied: F3

- files_changed:
  - `mathematica/moving_throat_pde_stage030_selected_mode_normalization_mathematica_audit.wl`
- summary: Re-sourced lamMinus and lamPlus from Eigenvalues[mMat] in the Mathematica Part II block while preserving the existing overlap checks.
- deviation: none
