---
unit_id: 039
batch: III.1
created_at: 2026-05-26T00:00:00Z
findings_count: 2
stop_cold: null
applied: true
applied_at: 2026-05-26T07:36:40Z
findings_applied: 2
findings_blocked: 0
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 039 (v2)

Apply each finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

If a finding's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named.

Do NOT run python or mathematica. Only edit files.

Do NOT touch paper.tex, notes/, or any prose documents. The red-team only modifies scripts.

Note: an earlier v1 directive for stage 039 has already been applied (it introduced the `*_derived` / `*_postulated` ordering in the Mathematica script and replaced the old `z1/z0 - (kappa1/kappa0) R_U` check with the explicit `(1+rho0)` / `(1+rho0/(1+deltaU))` form). Both findings below address residual issues not covered by that earlier round.

## F1 — tautological_check

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage039_split_u_sector_sympy_audit.py:138-139`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage039_split_u_sector_mathematica_audit.wl:118-119`

**Issue:** Section 22.4 of the scripts defines `M_mix_split = 8 Z_W (1+rho0)^2 / (pi^2 (1-eps_eta)(1-eps_W_split))` by *literally substituting* `eps_W_split` into the flat formula, and then asserts that `M_mix_split == M_mix_flat.subs(eps_W, eps_W_split)`. By construction of `M_mix_split`, the residual is identically zero regardless of the physics. The same holds for `R_target_split` vs. `R_target_flat.subs(eps_W, eps_W_split)`. The paper's `\stagefield{Output}` does not require these placement-map identities (only the notes mention them), and the assertion as written carries no verification content.

**Required change:**

In the SymPy script `moving_throat_pde_stage039_split_u_sector_sympy_audit.py`:

1. Delete the two `expect_zero` calls at lines 138 and 139 entirely.

   Specifically, remove:

   ```python
   expect_zero("M_mix split is M_mix_flat under eps_W -> eps_W_split", M_mix_split - M_mix_flat.subs(eps_W, eps_W_split))
   expect_zero("R_target split is R_target_flat under eps_W -> eps_W_split", R_target_split - R_target_flat.subs(eps_W, eps_W_split))
   ```

   Keep the `print("M_mix^(split U) =", M_mix_split)`, `print("R_target^(split U) =", R_target_split)`, and `print("product =", product)` lines at lines 135-137 unchanged — they are documentation prints, not assertions, and so are not tautological.

2. After deleting the two assertions, the definitions of `M_mix_flat` and `R_target_flat` at lines 129-130 are still used to print/document but no longer referenced by any assertion. That is fine; leave them in place for documentation continuity.

In the Mathematica script `moving_throat_pde_stage039_split_u_sector_mathematica_audit.wl`:

1. Delete the two `expectZero` calls at lines 118 and 119 entirely.

   Specifically, remove:

   ```mathematica
   expectZero["M_mix split is M_mix_flat under epsW -> epsWSplit", mMixSplit - (mMixFlat /. epsW -> epsWSplit)];
   expectZero["R_target split is R_target_flat under epsW -> epsWSplit", rTargetSplit - (rTargetFlat /. epsW -> epsWSplit)];
   ```

   Keep `Print["M_mix^(split U) = ", fmt[mMixSplit]]`, `Print["R_target^(split U) = ", fmt[rTargetSplit]]`, and `Print["product = ", fmt[product]]` at lines 115-117 unchanged.

2. The definitions of `mMixFlat` and `rTargetFlat` at lines 109-110 may stay as documentation.

**Verification command:**

After Codex applies, the verifier will run `redteam exec-sympy 039` and `redteam exec-mathematica 039` and confirm:
- The SymPy output no longer contains the lines `M_mix split is M_mix_flat under eps_W -> eps_W_split = 0` and `R_target split is R_target_flat under eps_W -> eps_W_split = 0` (visible in the current saved output at lines 43-44).
- The Mathematica output no longer contains `PASS: M_mix split is M_mix_flat under epsW -> epsWSplit` or `PASS: R_target split is R_target_flat under epsW -> epsWSplit` (visible in the current saved output at lines 50, 52).
- Both scripts still exit `0`.
- All other assertions (A1-A5 SymPy, A8-A13 Mathematica per the report) still pass with residual `0`.

## Applied: F1

- files_changed:
  - `scripts/moving_throat_pde_stage039_split_u_sector_sympy_audit.py`
  - `mathematica/moving_throat_pde_stage039_split_u_sector_mathematica_audit.wl`
- summary: Removed the tautological split-vs-flat placement-map assertions while preserving the documentation prints and flat definitions.
- deviation: none

## F2 — insufficient_verification

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage039_split_u_sector_sympy_audit.py:124` (immediately after the existing `expect_zero("direction-splitting invariant", D_dir - D_dir_expected)` at line 122 and before the print at line 124)
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage039_split_u_sector_mathematica_audit.wl:105-106` (immediately after the existing `expectZero["direction-splitting invariant derived matches postulated", ...]` at line 105 and before the `subbanner` at line 107)

**Issue:** The paper's `\stagefield{Output}` (line 54 of `paper/stages/stage_039.tex`) lists the collinearity iff theorem (eq:app-stage039-collinearity) as a top-level deliverable: `D_dir = 0 <=> delta_U = 0 or rho_0 = 0`. The current scripts derive the closed form of `D_dir` and print a comment claiming the iff theorem, but issue no formal assertion of either direction. The "if" direction (substitute `deltaU = 0` or `rho0 = 0` and confirm `D_dir` vanishes) is trivial to add; the "only if" direction requires confirming that the numerator of `D_dir` has `rho0 * deltaU` as a factor and that the residual after pulling out that factor is nonvanishing (i.e., not zero as a function of `rho0`, `deltaU`).

**Required change:**

In the SymPy script `moving_throat_pde_stage039_split_u_sector_sympy_audit.py`:

1. Locate the line `expect_zero("direction-splitting invariant", D_dir - D_dir_expected)` (line 122) and the immediately-following `print("Collinearity theorem: D_dir = 0 iff deltaU = 0 or rho0 = 0 (equivalently g_R g_U = 0).")` (line 124).

2. Between line 122 and line 124, insert exactly the following block (preserving the existing 4-space indentation level, which in this script is column 0):

   ```python
   # Collinearity iff theorem (paper eq:app-stage039-collinearity).
   expect_zero("collinearity if-leg: D_dir(deltaU=0) = 0", D_dir.subs(deltaU, 0))
   expect_zero("collinearity if-leg: D_dir(rho0=0) = 0", D_dir.subs(rho0, 0))
   # Only-if leg: the numerator of D_dir must factor exactly as rho0 * deltaU * <nonzero>.
   num, den = sp.fraction(sp.together(D_dir))
   num_reduced = sp.simplify(num / (rho0 * deltaU))
   print("Numerator(D_dir) / (rho0*deltaU) =", num_reduced)
   if num_reduced.has(rho0) or num_reduced.has(deltaU):
       raise AssertionError(
           f"D_dir numerator still depends on rho0 or deltaU after dividing by rho0*deltaU: {num_reduced}"
       )
   if sp.simplify(num_reduced) == 0:
       raise AssertionError("D_dir numerator is identically zero after factoring out rho0*deltaU")
   ```

   The `if-leg` assertions are non-tautological because `D_dir` was derived at line 119 from `kappa0 * z1 - kappa1 * z0`, where `z0` and `z1` were independently constructed; the substitution `rho0 -> 0` collapses both `z0` and `z1` to their bare overlap terms, and `kappa0 z1 - kappa1 z0` then reduces to `kappa0 kappa1 g_W - kappa1 kappa0 g_W = 0`. Similarly `deltaU -> 0` collapses `(1 + rho0/(1+deltaU)) -> (1+rho0)`, after which `z1 = kappa1 g_W (1+rho0)` and `kappa0 z1 - kappa1 z0 = kappa0 kappa1 g_W (1+rho0) - kappa1 kappa0 g_W (1+rho0) = 0`. Both `if`-legs vanish for non-trivial algebraic reasons. The only-if leg uses `sp.fraction(sp.together(D_dir))` to get the numerator, divides by the conjectured factor `rho0 * deltaU`, and checks the quotient is free of `rho0` and `deltaU` (so they were the only `rho0`/`deltaU`-dependent factors in the numerator). Cross-check: the printed `D_dir` is `8*sqrt(2)*c_etaW*deltaU*rho0/(3*pi**2*sqrt(mu_W)*sqrt(mu_eta)*(deltaU + 1))` (from the saved sympy output at line 33), so `sp.fraction(sp.together(D_dir))` returns numerator `8*sqrt(2)*c_etaW*deltaU*rho0` and denominator `3*pi**2*sqrt(mu_W)*sqrt(mu_eta)*(deltaU + 1)`. Then `num / (rho0*deltaU) = 8*sqrt(2)*c_etaW / (3*pi**2*sqrt(mu_W)*sqrt(mu_eta))`, which is `has(rho0) == False`, `has(deltaU) == False`, and is nonzero (`c_etaW > 0`). All three new checks pass.

In the Mathematica script `moving_throat_pde_stage039_split_u_sector_mathematica_audit.wl`:

1. Locate `expectZero["direction-splitting invariant derived matches postulated", dDirDerived - dDirPostulated];` at line 105.

2. Immediately after line 105 (before the blank line at line 106 and the `subbanner` at line 107), insert exactly the following block (at column 0):

   ```mathematica
   (* Collinearity iff theorem (paper eq:app-stage039-collinearity). *)
   expectZero["collinearity if-leg: D_dir(deltaU=0) = 0", dDir /. deltaU -> 0];
   expectZero["collinearity if-leg: D_dir(rho0=0) = 0", dDir /. rho0 -> 0];
   (* Only-if leg: numerator of D_dir must factor as rho0 * deltaU * <nonzero, no rho0, no deltaU>. *)
   dDirTogether = Together[dDir];
   dDirNum = Numerator[dDirTogether];
   dDirNumReduced = FullSimplify[dDirNum/(rho0 deltaU), Assumptions -> $Assumptions];
   Print["Numerator(D_dir) / (rho0*deltaU) = ", fmt[dDirNumReduced]];
   If[!FreeQ[dDirNumReduced, rho0] || !FreeQ[dDirNumReduced, deltaU],
     fail["collinearity only-if: numerator still depends on rho0 or deltaU after factoring", dDirNumReduced]];
   If[TrueQ[FullSimplify[dDirNumReduced == 0, Assumptions -> $Assumptions]],
     fail["collinearity only-if: numerator is identically zero after factoring rho0*deltaU", dDirNumReduced]];
   pass["collinearity only-if: residual factor is nonzero and independent of rho0, deltaU"];
   ```

   Cross-check against the saved Mathematica output line 39 (`D_dir = (8*Sqrt[2]*cEtaW*deltaU*rho0)/(3*(1 + deltaU)*Sqrt[muEta*muW]*Pi^2)`): `Numerator[Together[dDir]] = 8*Sqrt[2]*cEtaW*deltaU*rho0`, so `dDirNum/(rho0 deltaU) = 8*Sqrt[2]*cEtaW/(3*Sqrt[muEta*muW]*Pi^2)`. `FreeQ[..., rho0] == True`, `FreeQ[..., deltaU] == True`, and the expression is not zero under the `$Assumptions` (`cEtaW > 0`, `muEta > 0`, `muW > 0`). All three new checks pass.

**Verification command:**

After Codex applies, the verifier will run `redteam exec-sympy 039` and `redteam exec-mathematica 039` and confirm:
- SymPy output contains three new lines at section 22.3:
  - `collinearity if-leg: D_dir(deltaU=0) = 0 = 0`
  - `collinearity if-leg: D_dir(rho0=0) = 0 = 0`
  - `Numerator(D_dir) / (rho0*deltaU) = 8*sqrt(2)*c_etaW/(3*pi**2*sqrt(mu_W)*sqrt(mu_eta))` (or an equivalent form)
- Mathematica output contains three new lines at section 3:
  - `PASS: collinearity if-leg: D_dir(deltaU=0) = 0`
  - `PASS: collinearity if-leg: D_dir(rho0=0) = 0`
  - `Numerator(D_dir) / (rho0*deltaU) = (8*Sqrt[2]*cEtaW)/(3*Sqrt[muEta*muW]*Pi^2)` followed by `PASS: collinearity only-if: residual factor is nonzero and independent of rho0, deltaU`
- Both scripts still exit `0`.
- All previously-passing assertions still pass.

## Applied: F2

- files_changed:
  - `scripts/moving_throat_pde_stage039_split_u_sector_sympy_audit.py`
  - `mathematica/moving_throat_pde_stage039_split_u_sector_mathematica_audit.wl`
- summary: Added formal collinearity if-leg checks and only-if numerator factor checks to both audit scripts.
- deviation: none
