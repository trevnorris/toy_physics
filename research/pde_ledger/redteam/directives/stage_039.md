---
unit_id: 039
batch: III.1
created_at: 2026-05-22T07:59:49Z
findings_count: 3
stop_cold: null
applied: true
applied_at: 2026-05-22T18:24:37Z
findings_applied: 3
findings_blocked: 0
verification_status: pending
---

# Codex directive — unit 039

Apply each finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

If a finding's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named.

Do NOT run python or mathematica. Only edit files.

Do NOT touch paper.tex, notes/, or any prose documents. The red-team only modifies scripts.

## F1 — tautological_check

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage039_split_u_sector_sympy_audit.py:107-114`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage039_split_u_sector_mathematica_audit.wl:82-89`

**Issue:** The assertion `expect_zero("z1/z0 - (kappa1/kappa0) R_U", z1/z0 - (kappa1/kappa0)*R_U)` (sympy:114) and its Mathematica counterpart (mma:89) are tautological by construction. `R_U` is defined as `(1 + rho0/(1+deltaU))/(1+rho0)`, and `z0`, `z1` are constructed with `(1+rho0)` and `(1+rho0/(1+deltaU))` respectively, so `z1/z0` literally equals `(kappa1/kappa0)*R_U` by substitution of the literal definitions. The check can't fail and contributes no verification content.

**Required change:**

In the SymPy script `moving_throat_pde_stage039_split_u_sector_sympy_audit.py`:

1. At line 114, replace the existing assertion

   ```
   expect_zero("z1/z0 - (kappa1/kappa0) R_U", sp.simplify(z1 / z0 - (kappa1 / kappa0) * R_U))
   ```

   with

   ```
   expect_zero(
       "z1*(1+rho0) - (kappa1/kappa0)*z0*(1+rho0/(1+deltaU))",
       sp.simplify(z1 * (1 + rho0) - (kappa1 / kappa0) * z0 * (1 + rho0 / (1 + deltaU))),
   )
   ```

   The new check substitutes `z0`, `z1` against their *explicit* `(1+rho0)` / `(1+rho0/(1+deltaU))` rho-structure rather than against the named symbol `R_U`. After simplification the residual must reduce to `0` only if `z0`, `z1` have the right kappa-rho construction. Keep the `print("R_U =", R_U)` (line 113) as informational.

In the Mathematica script `moving_throat_pde_stage039_split_u_sector_mathematica_audit.wl`:

1. At line 89, replace

   ```
   expectZero["z1/z0 - (kappa1/kappa0) R_U", z1/z0 - (kappa1/kappa0) rU];
   ```

   with

   ```
   expectZero[
     "z1*(1+rho0) - (kappa1/kappa0)*z0*(1+rho0/(1+deltaU))",
     z1*(1 + rho0) - (kappa1/kappa0)*z0*(1 + rho0/(1 + deltaU))
   ];
   ```

   Keep `Print["R_U = ", fmt[rU]]` (line 88) as informational.

**Verification command:**

After Codex applies, the verifier will run `redteam exec-sympy 039` and `redteam exec-mathematica 039` and confirm the new assertion lines appear, the scripts exit `0`, and the output prints the new label `"z1*(1+rho0) - (kappa1/kappa0)*z0*(1+rho0/(1+deltaU)) = 0"` instead of the old label.

## Applied: F1

- files_changed:
  - `scripts/moving_throat_pde_stage039_split_u_sector_sympy_audit.py`
  - `mathematica/moving_throat_pde_stage039_split_u_sector_mathematica_audit.wl`
- summary: Replaced the ratio-vs-R_U check with an explicit z0/z1 kappa-rho structure residual in both audit scripts.
- deviation: none

## F2 — tautological_check

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage039_split_u_sector_sympy_audit.py:126-133`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage039_split_u_sector_mathematica_audit.wl:98-105`

**Issue:** The product-law assertion `expect_zero("product law", M_mix_split * R_target_split - 8*Lambda*(1-eps_W_split)/pi**2)` (sympy:133) and its Mathematica counterpart (mma:105) are algebraic identities by direct factor cancellation. `Z_W`, `(1+rho0)^2`, `(1-eps_eta)`, and one factor of `(1-eps_W_split)` cancel between `M_mix_split` and `R_target_split` literally, leaving `8*Lambda*(1-eps_W_split)/pi^2` by elementary algebra regardless of any physics. The script docstring item 5 claims this verifies that "the Stage-21 continuum placement map survives at the scalar level with `eps_W -> eps_W_split` and `delta -> delta_split`", which is a substitution claim, not a product-law cancellation.

**Required change:**

In the SymPy script `moving_throat_pde_stage039_split_u_sector_sympy_audit.py`:

1. Between line 125 and line 126, insert the flat-U baseline definitions:

   ```
   M_mix_flat = sp.simplify(8 * Z_W * (1 + rho0)**2 / (pi**2 * (1 - eps_eta) * (1 - eps_W)))
   R_target_flat = sp.simplify(Lambda * (1 - eps_eta) * (1 - eps_W)**2 / (Z_W * (1 + rho0)**2))
   ```

2. At line 133, replace

   ```
   expect_zero("product law", product - 8 * Lambda * (1 - eps_W_split) / pi**2)
   ```

   with the two substitution checks

   ```
   expect_zero("M_mix split is M_mix_flat under eps_W -> eps_W_split", M_mix_split - M_mix_flat.subs(eps_W, eps_W_split))
   expect_zero("R_target split is R_target_flat under eps_W -> eps_W_split", R_target_split - R_target_flat.subs(eps_W, eps_W_split))
   ```

   Keep `print("product =", product)` (line 132) as informational so the closed form remains visible.

In the Mathematica script `moving_throat_pde_stage039_split_u_sector_mathematica_audit.wl`:

1. Between line 97 (the `subbanner` for section 4) and line 98 (the `mMixSplit` definition), insert the flat-U baseline definitions:

   ```
   mMixFlat = FullSimplify[8 zW (1 + rho0)^2/(Pi^2 (1 - epsEta) (1 - epsW)), Assumptions -> $Assumptions];
   rTargetFlat = FullSimplify[lambda (1 - epsEta) (1 - epsW)^2/(zW (1 + rho0)^2), Assumptions -> $Assumptions];
   ```

2. At line 105, replace

   ```
   expectZero["product law", product - 8 lambda (1 - epsWSplit)/Pi^2];
   ```

   with

   ```
   expectZero["M_mix split is M_mix_flat under epsW -> epsWSplit", mMixSplit - (mMixFlat /. epsW -> epsWSplit)];
   expectZero["R_target split is R_target_flat under epsW -> epsWSplit", rTargetSplit - (rTargetFlat /. epsW -> epsWSplit)];
   ```

   Keep `Print["product = ", fmt[product]]` (line 104) as informational.

**Verification command:**

After Codex applies, the verifier will run `redteam exec-sympy 039` and `redteam exec-mathematica 039` and confirm both new substitution assertions appear and pass (residuals print as `0`), and that the old `"product law"` assertion is gone. The script must still exit `0`.

## Applied: F2

- files_changed:
  - `scripts/moving_throat_pde_stage039_split_u_sector_sympy_audit.py`
  - `mathematica/moving_throat_pde_stage039_split_u_sector_mathematica_audit.wl`
- summary: Added flat-U baseline definitions and replaced the product-law cancellation check with split-vs-flat substitution checks in both audit scripts.
- deviation: none

## F3 — mathematica_transliteration

**Target:** `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage039_split_u_sector_mathematica_audit.wl:55-94`

**Issue:** The Mathematica script postulates the same closed-form answers (`deltaSplit`, `epsWSplit`, `dDirExpected`) as the SymPy script and then verifies them against direct constructions whose algebra mirrors the SymPy choreography step for step. As a result, the second-engine check cannot surface any disagreement other than a transcription typo. The independent-derivation policy requires the Mathematica engine to obtain at least the key derived quantities from its own algebra rather than echoing SymPy's postulates.

**Required change:**

In `moving_throat_pde_stage039_split_u_sector_mathematica_audit.wl`:

1. **deltaSplit (lines 59, 61, 67)**: Replace the postulated form with a Mathematica-derived form.

   At line 59, replace

   ```
   deltaSplit = FullSimplify[(delta0 + epsEta deltaU/(1 + deltaU))/(1 - epsEta), Assumptions -> $Assumptions];
   ```

   with

   ```
   a1Direct = FullSimplify[(a1 /. cEtaU^2 -> epsEta kU kEtaEff), Assumptions -> $Assumptions];
   deltaSplitDerived = FullSimplify[a1Direct/a0Expected - 1, Assumptions -> $Assumptions];
   deltaSplitPostulated = (delta0 + epsEta deltaU/(1 + deltaU))/(1 - epsEta);
   deltaSplit = deltaSplitDerived;
   ```

   At line 67 (the existing `expectZero["A1 direct - expected", ...]`), leave as is. Immediately after line 67, append a new check:

   ```
   expectZero["delta_split derived matches postulated", deltaSplitDerived - deltaSplitPostulated];
   ```

2. **epsWSplit (lines 73, 77)**: Replace the postulated form with a Mathematica-derived form.

   At line 73, replace

   ```
   epsWSplit = FullSimplify[epsW (1 - (2/11) deltaU/(1 + deltaU)), Assumptions -> $Assumptions];
   ```

   with

   ```
   epsWSplitDerived = FullSimplify[(epsWDirect /. cUW^2 -> epsW kU kWEff/sigma), Assumptions -> $Assumptions];
   epsWSplitPostulated = epsW (1 - (2/11) deltaU/(1 + deltaU));
   epsWSplit = epsWSplitDerived;
   ```

   At line 77, replace the existing assertion

   ```
   expectZero["eps_W direct - split formula", (epsWDirect /. cUW^2 -> epsW kU kWEff/sigma) - epsWSplit];
   ```

   with

   ```
   expectZero["eps_W_split derived matches postulated", epsWSplitDerived - epsWSplitPostulated];
   ```

3. **dDirExpected (lines 91-94)**: Replace the postulated form with a Mathematica-derived form.

   At line 91, replace

   ```
   dDir = FullSimplify[kappa0 z1 - kappa1 z0, Assumptions -> $Assumptions];
   dDirExpected = FullSimplify[-kappa0 kappa1 gW rho0 deltaU/(1 + deltaU), Assumptions -> $Assumptions];
   ```

   with

   ```
   dDir = FullSimplify[kappa0 z1 - kappa1 z0, Assumptions -> $Assumptions];
   dDirDerived = dDir;
   dDirPostulated = FullSimplify[-kappa0 kappa1 gW rho0 deltaU/(1 + deltaU), Assumptions -> $Assumptions];
   dDirExpected = dDirDerived;
   ```

   At line 94, replace

   ```
   expectZero["direction-splitting invariant", dDir - dDirExpected];
   ```

   with

   ```
   expectZero["direction-splitting invariant derived matches postulated", dDirDerived - dDirPostulated];
   ```

**Claim manifest:** Not applicable — `.wl` file already exists; this is a restructuring, not a missing-script finding.

**Verification command:**

After Codex applies, the verifier will run `redteam exec-mathematica 039` and confirm: (a) the three new `derived matches postulated` assertions appear in the output and print residual `0`; (b) the original `A0 direct - expected`, `A1 direct - expected` assertions remain and still pass; (c) the script exits `0`. The verifier should also confirm via grep that `deltaSplitDerived`, `epsWSplitDerived`, `dDirDerived` symbols appear as left-hand sides of definitions before being asserted equal to their postulated counterparts.

## Applied: F3

- files_changed:
  - `mathematica/moving_throat_pde_stage039_split_u_sector_mathematica_audit.wl`
- summary: Introduced derived Mathematica forms for deltaSplit, epsWSplit, and dDir before checking each against its postulated counterpart.
- deviation: none
