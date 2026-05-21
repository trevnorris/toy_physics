---
unit_id: 002
batch: I.1
created_at: 2026-05-20T00:00:00Z
findings_count: 2
stop_cold: null
applied: true
applied_at: 2026-05-21T00:44:31-06:00
findings_applied: 2
findings_blocked: 0
verification_status: pending
---

# Codex directive — unit 002

Apply each finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

If a finding's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named.

Do NOT run python or mathematica. Only edit files.

Do NOT touch paper.tex, notes/, or any prose documents. The red-team only modifies scripts.

## F1 — mathematica_transliteration

**Target:** `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage002_breathing_reduction_mathematica_audit.wl:82-187`

**Issue:** The Mathematica script reproduces the SymPy script's algebra line-by-line: identical intermediate variable names (`y00`, `avgY00`, `mouthAvg`, `angularPrefactor`, `Mintegrand`, `Kintegrand`), identical 2x2 hand-built matrix structure with `4 Pi` pulled outside, identical hand-implemented Euler-Lagrange formula, and byte-identical `expectZero[...]` name strings. This is a structural transliteration rather than an independent second-engine derivation.

**Required change:**
Make three local edits to remove the line-by-line correspondence with the SymPy script. Do not change which physical statements are checked — only how Mathematica obtains them.

1. **Section I — use the built-in spherical harmonic.**
   - At WL line 82, replace `y00 = 1/(2 Sqrt[Pi]);` with:
     ```
     y00 = FullSimplify[SphericalHarmonicY[0, 0, theta, phi]];
     ```
   - The existing `expectZero["norm(Y00) - 1", normY00 - 1]` on line 94 then becomes a non-trivial check that Mathematica's built-in `Y_00` normalizes correctly, rather than a self-check on a hardcoded literal.

2. **Section II — extract M and K via `Coefficient`, do not hand-build them.**
   - At WL lines 133-145, replace the hand-built `Mintegrand`, `Kintegrand`, and `lwTarget` block with an extraction:
     ```
     MaaExt = (1/(dadt^2)) Coefficient[lw, dadt, 2];
     MLLExt = (1/(dLdt^2)) Coefficient[lw, dLdt, 2];
     MaLExt = (1/2) Coefficient[Coefficient[lw, dadt], dLdt];
     KaaExt = -(1/(da^2)) Coefficient[lw, da, 2];
     KLLExt = -(1/(dL^2)) Coefficient[lw, dL, 2];
     KaLExt = -(1/2) Coefficient[Coefficient[lw, da], dL];
     MintegrandExtracted = {{MaaExt, MaLExt}, {MaLExt, MLLExt}};
     KintegrandExtracted = {{KaaExt, KaLExt}, {KaLExt, KLLExt}};
     MintegrandBoxed = 4 Pi {
       {muEta alphaA^2, muEta alphaA alphaL},
       {muEta alphaA alphaL, muEta alphaL^2}
     };
     KintegrandBoxed = 4 Pi {
       {Tw D[alphaA, w]^2 + K0 alphaA^2, Tw D[alphaA, w] D[alphaL, w] + K0 alphaA alphaL},
       {Tw D[alphaA, w] D[alphaL, w] + K0 alphaA alphaL, Tw D[alphaL, w]^2 + K0 alphaL^2}
     };
     ```
   - Then on line 146 (the existing `expectZero` for "reduced Lagrangian density from the action - target overlap form"), replace it with two checks:
     ```
     expectZero["extracted M matches boxed M (4 Pi overlap form)", MintegrandExtracted - MintegrandBoxed];
     expectZero["extracted K matches boxed K (4 Pi overlap form)", KintegrandExtracted - KintegrandBoxed];
     ```
   - This way Mathematica derives `M` and `K` from the action by coefficient extraction and then compares to the boxed form, rather than postulating both and checking they agree (which is what the SymPy script does and what the current WL ports).

3. **Section II — use `VariationalMethods` for the Euler-Lagrange equations.**
   - Add `Needs["VariationalMethods\`"];` near the top of the script (after line 2 is fine).
   - At WL lines 177-178 (the calls to the local helper `eulerLagrange1D`), replace with:
     ```
     {elAEq, elLEq} = EulerEquations[lredTime, {qa, qLfun}, t];
     elA = elAEq[[1]] - elAEq[[2]];
     elL = elLEq[[1]] - elLEq[[2]];
     ```
     (`EulerEquations` returns `{lhs == rhs, ...}`; rearrange to a residual.)
   - The local helper `eulerLagrange1D` defined on lines 59-62 can stay for use by Section III's single-component EL check on line 277.

**Verification command:**
After Codex applies, the verifier will run `redteam exec-mathematica 002` and confirm:
- The new `extracted M matches boxed M (4 Pi overlap form)` and `extracted K matches boxed K (4 Pi overlap form)` PASS lines appear in the output.
- `Y_00` is no longer hardcoded as `1/(2 Sqrt[Pi])`; line 82 now references `SphericalHarmonicY`.
- The Euler-Lagrange equation checks still PASS, but now via `EulerEquations` rather than the hand-built formula.
- The script exits 0.

## Applied: F1

- files_changed:
  - `mathematica/moving_throat_pde_stage002_breathing_reduction_mathematica_audit.wl`
- summary: Replaced the hardcoded monopole, hand-built two-mode matrix extraction, and local two-mode EL derivation with Mathematica built-in `SphericalHarmonicY`, coefficient extraction, and `EulerEquations`.
- deviation: Used Mathematica's coefficient factors and residual sign convention instead of the literal snippet so the extracted matrices and EL residuals reduce to zero.

## F2 — insufficient_verification

**Target:** `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage002_breathing_reduction_mathematica_audit.wl:200-267`

**Issue:** Section III claims to verify the full five-component P2 multiplet (Y20, Y21c, Y21s, Y22c, Y22s) is orthonormal with common angular-stiffness eigenvalue 6 and reduces to a degenerate five-component isotropic system. The Mathematica script only exercises three of the five components (Y20, Y21c, Y22c). The `s` companions Y21s and Y22s are defined but only checked as `phi`-rotations of their `c` siblings, not directly integrated.

**Required change:**

1. **Extend the test basis to all five components.**
   - At WL line 200, change
     ```
     basis = {y20, y21c, y22c};
     basisNames = {"Y20", "Y21c", "Y22c"};
     ```
     to
     ```
     basis = {y20, y21c, y21s, y22c, y22s};
     basisNames = {"Y20", "Y21c", "Y21s", "Y22c", "Y22s"};
     ```
   - The existing `Do[..., {i, 1, Length[basis]}]` loops on lines 206-219 and 221-227 will then automatically pick up Y21s and Y22s.

2. **Add explicit 5x5 matrix-form orthonormality and stiffness assertions.**
   - After the existing `Do[ ... , {i, 1, Length[basis]}]` block that ends at WL line 219 (the per-component norm and angular-energy checks), insert:
     ```
     normMatrix5 = Table[
       FullSimplify[
         Integrate[basis[[i]] basis[[j]] dOmega, {phi, 0, 2 Pi}, {theta, 0, Pi}],
         Assumptions -> $Assumptions
       ],
       {i, 1, Length[basis]}, {j, 1, Length[basis]}
     ];
     gradMatrix5 = Table[
       FullSimplify[
         Integrate[gradS2Inner[basis[[i]], basis[[j]], theta, phi] dOmega, {phi, 0, 2 Pi}, {theta, 0, Pi}],
         Assumptions -> $Assumptions
       ],
       {i, 1, Length[basis]}, {j, 1, Length[basis]}
     ];
     expectZero["real P2 norm matrix - IdentityMatrix[5]", normMatrix5 - IdentityMatrix[5]];
     expectZero["real P2 angular stiffness matrix - 6 IdentityMatrix[5]", gradMatrix5 - 6 IdentityMatrix[5]];
     ```

3. **Extend the per-component reduced-density loop to all five.**
   - At WL line 235 change
     ```
     qvec = {q20, q21c, q22c};
     qdotvec = {q20d, q21cd, q22cd};
     ```
     to
     ```
     qvec = {q20, q21c, q21s, q22c, q22s};
     qdotvec = {q20d, q21cd, q21sd, q22cd, q22sd};
     ```
   - The `Do[...]` loop on lines 238-267 will then iterate over all five components.

4. **Leave the phase-shift identities alone.**
   - WL lines 203-204 (`expectZero["phase shift: Y21s(theta,phi) - Y21c(theta,phi-Pi/2)", ...]` and the Y22 analogue) stay as supplementary checks. They are not the problem; they should not be removed.

**Verification command:**
After Codex applies, the verifier will run `redteam exec-mathematica 002` and confirm:
- Two new top-level assertions `real P2 norm matrix - IdentityMatrix[5]` and `real P2 angular stiffness matrix - 6 IdentityMatrix[5]` appear in the output and PASS.
- Per-component `norm(Y21s) - 1`, `angular energy(Y21s) - 6`, `-Delta_S2 Y21s - 6 Y21s`, `norm(Y22s) - 1`, `angular energy(Y22s) - 6`, `-Delta_S2 Y22s - 6 Y22s` PASS lines appear.
- Per-component `single-component reduced density for Y21s` and `... for Y22s` PASS lines appear.
- The script exits 0.

## Applied: F2

- files_changed:
  - `mathematica/moving_throat_pde_stage002_breathing_reduction_mathematica_audit.wl`
- summary: Extended the real P2 checks and reduced-density loop to all five components, adding explicit 5x5 norm and angular-stiffness matrix assertions.
- deviation: none
