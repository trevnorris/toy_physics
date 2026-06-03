---
unit_id: 239
batch: VII.2
created_at: 2026-06-02T00:00:00-06:00
findings_count: 2
stop_cold: null
applied: true
applied_at: 2026-06-03T08:38:12-06:00
findings_applied: 2
findings_blocked: 0
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 239

Apply each finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

If a finding's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes beyond what is named. Do NOT touch paper.tex, notes/, or any prose documents.

After editing, RUN the affected scripts (`python3 <path>` for SymPy, `math -script <path>` for Mathematica) and iterate until they exit 0 with all in-file checks passing.

## F1 — mathematica_transliteration

**Target:** `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage239_rigid_mouth_physical_normal_form_exact_physical_to_microscopic_correction_compiler_and_cartesian_orbit_lock_theorem_mathematica_audit.wl` (re-author in place; KEEP the `_mathematica_audit.wl` filename suffix)

**Issue:** The current `.wl` is a mechanical line-by-line port of the `.py`: identical camelCased variable names, identical assertion ordering, identical verbatim check-label strings, the identical hardcoded `SrmDep . Mphys` compiler decomposition, the identical `Solve`-then-compare-to-a-literal-`expectedInverse` choreography, and the identical first-order `D[...,h]/.h->0` trick. It therefore provides no independent second-engine confirmation — it would reproduce any algebra error the `.py` made because it makes the identical structural choices. A checkpoint requires a genuinely independent route.

**Required change:**
Re-derive the SAME physical claims (manifest M1-M9 below) from the paper's premises using a DIFFERENT decomposition and native Mathematica primitives. The acceptance bar is *independence of route*, not a different answer. Concretely, you MUST NOT reuse:
- the `SrmDep = {{0,0},{0,-1},{1,-1}}` then `CphysDep = SrmDep . Mphys` two-step construction. Instead build the physical-to-microscopic map a different way — e.g. define the dependent correction directly as the symbolic vector `{0, -V, U - V}` from the paper's boxed result and then *recover* the compiler matrix by `D[#, {{U, V}}] &` (Jacobian / `CoefficientArrays` / `Grad`), and separately confirm it factors through `M_phys = I_2`;
- the `Solve[...] === expectedInverse` literal-list-comparison pattern for the left inverse. Instead obtain the left inverse natively (e.g. `PseudoInverse` of the integer compiler matrix, or `LinearSolve`/`Inverse` on the appropriate 2x2 reduction) and verify `Lphys . Cphys == IdentityMatrix[2]` and that it reproduces `U = Delta_mu - Delta_Keta`, `V = -Delta_Keta` symbolically;
- the verbatim English check-label strings copied from the `.py`. Use your own labels.

For the orbit-lock claim, prefer a native `Reduce[yDep[[2]] == 0 && yDep[[3]] == 0, {U, V}]` (give `Reduce`, not the `.py`'s `Solve`-then-`=!=`-literal compare) and confirm it returns `U == 0 && V == 0`. For support-blindness, the independent route may use abstract functions with `D[..., zeta]` and explicit `Derivative` zero rules (this primitive is shared by necessity), but the surrounding choreography and labels must be your own.

Keep the existing `expectZero` / `pass` / `fail` harness if convenient, but the derivation BODY must be independent.

**Claim manifest** (each must be independently verified):
- M1: Chart definitions `Exp[U] == T2/T2ref`, `Exp[V] == epsEta/epsEtaRef`; diagonal packet `q_nt = U`, `q_eta = V` with `M_phys = I_2`.
- M2: Target ratio `R_target/R_target_ref == ((1 - epsEtaRef Exp[V])/(1 - epsEtaRef)) Exp[-U]`, consistent with the selected-branch identity `R_target T^2 == Lambda0 (1 - epsEta)`.
- M3: Complementary projectors `P_T = diag(1,0)`, `P_eta = diag(0,1)`: idempotent, orthogonal, sum to `I_2`; and the transfer/dressing legs factor the target ratio (`ratioFromUV == ratioTransfer * ratioDressing`).
- M4: Physical-to-microscopic dependent compiler `(Delta_T, Delta_Keta, Delta_mu) == (0, -V, U - V)`.
- M5: Left inverse `U == Delta_mu - Delta_Keta`, `V == -Delta_Keta`, with `Lphys . Cphys == I_2` (obtained by an independent native method, not literal-list compare).
- M6: Axis images `(U,0) -> (0,0,U)` and `(0,V) -> -V (0,1,1)`; dependent correction splits as their sum.
- M7: Correction compilers `Delta_static == (0,0,-U)`, `Delta_eta_rest == (0,V,V)`, `Delta_orbit == (0,V,V-U)`, with `Delta_orbit == Delta_static + Delta_eta_rest` and `yDep + Delta_orbit == 0`.
- M8: Support-blindness: the Stage-238 branch formula for `T^2` has zero `zeta`/`M_mix` derivative, and this propagates (chain rule) to `U`, `V`, the dependent correction, and the static/orbit corrections.
- M9: Cartesian orbit lock `yDep[[2]]==0 && yDep[[3]]==0` <=> `U==0 && V==0`; `T2==T2ref => U==0`; `epsEta==epsEtaRef => V==0`; target ratio at `U=V=0` equals `1`; first-order form `(0, -dlnEps, dlnT2 - dlnEps)` via `D[..., h] /. h -> 0`.

**Verification command:** the verifier runs `redteam exec-mathematica 239` and confirms the re-authored `.wl` (a) no longer shares the `SrmDep . Mphys` construction, the `expectedInverse` literal-compare, or the `.py`'s verbatim labels, (b) covers M1-M9, and (c) exits 0.

## Applied: F1

- files_changed:
  - `mathematica/moving_throat_pde_stage239_rigid_mouth_physical_normal_form_exact_physical_to_microscopic_correction_compiler_and_cartesian_orbit_lock_theorem_mathematica_audit.wl`
- summary: Re-authored the Mathematica audit to recover the dependent compiler by Jacobian from the boxed correction vector, use `PseudoInverse` for the left inverse, and verify M1-M9 with native Mathematica checks.
- deviation: none

## F2 — stale stage label (script-side)

**Target:** `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage239_rigid_mouth_physical_normal_form_exact_physical_to_microscopic_correction_compiler_and_cartesian_orbit_lock_theorem_mathematica_audit.wl:59` and the `Stage221` variable suffixes (wl:141-203)

**Issue:** The banner reads `"STAGE 222 — RIGID-MOUTH PHYSICAL NORMAL FORM"` and the carried Stage-238 branch variables are suffixed `Stage221` (`T2Stage221`, `RtargetStage221`, `qNtRigidStage221`, `UStage221`, `VStage221`). Canonical numbering is 239 (paper card, script filename, appendix row 90) and the carried branch is Stage 238 (card `\stagefield{Inputs}`, notes §0). This is the known project-wide incomplete-renumber drift surfacing on the script side. The "STAGE 222" banner also propagates into the saved output header.

**Required change:**
- Change the banner string to `"STAGE 239 — RIGID-MOUTH PHYSICAL NORMAL FORM"`.
- Rename the `Stage221` variable suffixes to `Stage238` (consistently, including the `Clear[...]` list if affected and all use sites): `T2Stage221 -> T2Stage238`, `RtargetStage221 -> RtargetStage238`, `qNtRigidStage221 -> qNtRigidStage238`, `UStage221 -> UStage238`, `VStage221 -> VStage238`, `xPhysStage221 -> xPhysStage238`, `yDepStage221 -> yDepStage238`.
- This is a pure rename; the symbolic content is unchanged. If you re-author the `.wl` for F1, simply use the correct "Stage 239" banner and "Stage238" suffixes in the new file and mark F2 satisfied by F1.
- Do NOT batch-renumber any other file. This is a single-file label fix only. Do NOT touch the `.py` comments (their "Stage 236/221" references are the same drift but out of scope for this directive — leave them).

**Verification command:** `grep -i "STAGE 222\|Stage221" <wl>` returns nothing after the fix; `redteam exec-mathematica 239` exits 0 and the output header reads "STAGE 239".

## Applied: F2

- files_changed:
  - `mathematica/moving_throat_pde_stage239_rigid_mouth_physical_normal_form_exact_physical_to_microscopic_correction_compiler_and_cartesian_orbit_lock_theorem_mathematica_audit.wl`
  - `mathematica/output/moving_throat_pde_stage239_rigid_mouth_physical_normal_form_exact_physical_to_microscopic_correction_compiler_and_cartesian_orbit_lock_theorem_mathematica_audit.txt`
- summary: The re-authored audit and regenerated saved output use the `STAGE 239` banner and the Stage238 carried-branch variable suffixes throughout.
- deviation: none
