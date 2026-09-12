# Independent statement-fidelity review — S10 Lean contract

## 1. Reviewed packet

| Field | Value |
|---|---|
| Revision | `S10 compact fidelity contract v1` (MANIFEST.json:2) |
| Packet SHA-256 (as pinned) | `96d2c9d90a32b26d1c031e43ccc8580a8f8a102e18425491b231d4ac0cbc3b35` |
| Reviewed artifact | `lean/s10/COVERAGE.md` §"S10 work contract" (lines 29–177), sha256 `013cf4df…c1c5db` |
| Policy | `lean/FORMALIZATION_POLICY.md`, sha256 `065a3f04…9e58c8`, policy commit `c2bdb663` |
| Proof checkpoint | `f21e6459` |
| Author | Codex |
| Reviewer | independent non-author; this is one of the two L4 review legs |

I verified file contents directly against the manifest hashes for the files I opened; I could not recompute the aggregate `packet_sha256` (no shell in this review environment).

## 2. Verdict

**CLEAR** for the stated Lean scope (the conditional real plane-wave classification of the six supplied quadratic action families, its compact fidelity link, coverage contract, and mutation controls).

No substantive blocker found. All findings below are editorial. This clearance does not cover the second required review leg, nor the separate CAS production / comparator-export / ledger-paper tasks that the contract explicitly excludes.

## 3. What I checked and what it resolved to

### 3.1 Six families, action identity, parameter and index maps

The contract's common form `rho/2 · vᵀWv − mu·c/2 · S(J)` is exactly `S10Audit.packageAction_uses_stiffness` (`S10Audit/Packages.lean:38-45`), with `packageStiffness`/`stiffnessCoefficient`/`packageKinetic` at lines 16–36. Each branch matches its underlying definition:

- MAIN: `S10Pilot.stiffness J = (1/2)·Σ_ij (J i.succ j − J j.succ i)^2` (`S10Pilot/Action.lean:22-23`).
- FULLGRAD: `Σ_ij (J i.succ j)^2`; DIVONLY: `(Σ_i J i.succ i)^2` (`S10Controls/Action.lean:16-18`).
- SIGNFLIP / XCOEF_SCALE: `S10ScalarControls.lagrangian c rho mu J = rho/2·normSq(J 0) − c·mu/2·stiffness J` (`S10Controls/Scalar.lean:14-15`), instantiated at `c = -1` and `c = scale` (`Packages.lean:26,28`).
- ANISO: `kinetic e sigma v = Σ_i (if i = e then sigma else 1)·v_i^2` (`S10Anisotropic/Action.lean:13-17`) — i.e. only `W_ee` changes, all other inertia entries stay 1, as claimed.

Against the CAS builders:

- `stiffness_density` (`scripts/S10_brane_mode_spectrum_sympy_audit.py:338-348`) and `buildPackage` (`mathematica/S10_brane_mode_spectrum_mathematica_audit.wl:154-165`) give the same three densities with the same index convention `gradient[i][j] = ∂_{x_i} u_j`, matching Lean's `J i.succ j`, and `J 0` = velocity (confirmed via `fieldJet`/`waveCovector`, `S10Pilot/PlaneWave.lean:20-24`).
- `build_action` (py:374-395) / `buildPackage` (wl:167-194): `inertial_coefficients[0] = s_rho·rho_br`, others `rho_br`; `stiffness_coefficient = s·mu_R` for XCOEF_SCALE; `stiffness_sign = +1` for SIGNFLIP. The SIGNFLIP rewriting `L = T + (μ/2)S = T − ((−1)μ/2)S` is exactly `c = −1`; no sign slip.
- `PACKAGES` (py:85-92) is the six-entry selector the contract names, with the same stiffness/sign/aniso/scale flags as the Wolfram `Switch` blocks.
- Parameter map `rho_br|rhoBr → rho`, `mu_R|muR → mu`, `s_rho|sRho → sigma`, `s|coefficientScale → c` is correct.
- The CAS joint assumptions (py:351-371, wl:102-104) — `rho>0`, `mu>0`, `|k|²>0`, real `k`, real `a`, `sigma>0 ∧ sigma≠1`, `c>0` — map onto the contract's declared domain without slack. Lean's domain (arbitrary `D`, arbitrary `e : Fin D`) is a strict superset of the emitted dimensions `{2,3,4,5}` and of the CAS's fixed first component.

The contract's plane-wave operator `rho·z·W − mu·c·B(k)` with `B = K I − k kᵀ / K I / k kᵀ` is confirmed by `packageOperator` and `actionMatrix`/`actionMatrix_mulVec` (`S10Audit/MatrixTrees.lean:28-72`) and by each family's `modalOperator`.

### 3.2 Stationarity notion, derivative and normalization

`ActionStationary` is stationarity of a *relative* action against every smooth compactly supported vector field (`S10Pilot/Variation.lean:84-86,162-163`), upgraded to the pointwise Euler–Lagrange identity by `actionStationary_iff_eulerLagrange` (`Variation.lean:188`) and specialized by `actionStationary_planeWave_iff` (`Variation.lean:217`). The contract's sentences "against all admissible compact variations, not only variations within that ansatz" and "a relative action avoids assuming finite total action for a nondecaying plane wave" are accurate; `S10Pilot/FiniteAction.lean:30-41` supplies the finite-action compatibility separately. `phaseAverage_eq` (`S10Pilot/PhaseAverage.lean:29-36`, and the `S10Controls`/`S10Anisotropic`/`S10ScalarControls` copies) gives the exact `0..2π` integral with the explicit factor `1/2`.

### 3.3 Physical vs normalized squared frequency

The contract's distinction holds in the sources. `S10Anisotropic/Spectrum.lean:11-29` defines `normalizedOperator … z …` with `modalOperator_normalized` fixing `z = rho·omega^2/mu`; `normalized_frequency_iff` and `extraConeValue = (mu/rho)·extraValue` (`S10Anisotropic/Certificate.lean:13-24`) convert back. So `R = (mu/rho)K` and `E = (mu/rho)(p² + q/sigma)` are the physical branches, and "ANISO's internal normalized variable is `rho·z/mu`" is exactly right. `S10ScalarControls.nonzero_root_iff` (`Scalar.lean:93-116`) treats `z` as an arbitrary real, so the signed-`z` reading for SIGNFLIP is genuinely proved rather than inferred.

### 3.4 Complete amplitude classification and kernel dimensions

Every table cell resolves to a theorem with matching content:

| Contract cell | Source |
|---|---|
| MAIN `R`/`T`, `D-1/D-1`; static `L`, `1/0` | `s10_variational_certificate` (`VariationalCertificate.lean:12`), `propagating_variational_mode_iff`, `transverseSpace_finrank`, `longitudinalSpace_finrank`, `S10Controls.longitudinal_inf_transverse` |
| FULLGRAD `R`/⊤, `D/D-1`; static `⊥`, `0/0` | `full_on_cone_space`, `full_zero_space`, `full_mode_counts` (`S10Controls/Certificate.lean:28-36`) |
| DIVONLY `R`/`L`, `1/0`; static `T`, `D-1/D-1` | `div_on_cone_space`, `div_zero_space`, `div_mode_counts` (`Certificate.lean:38-46`) |
| SIGNFLIP `−R`/`T`, `D-1/D-1` | `signflip_counts`, `coneValue_neg`, `negative_control_no_real_wave` (`Scalar.lean:123-199`) |
| XCOEF_SCALE `cR`/`T`, `D-1/D-1` | `cone_modeSpace`, `cone_counts`, `zero_counts` (`Scalar.lean:136-207`) |
| ANISO parallel `R=E`, `O=T`, `D-1/D-1` | `parallel_kernel_iff`, `parallel_modeSpace`, `ordinarySpace_parallel(_finrank)`, `parallel_counts` |
| ANISO perpendicular `D-2/D-2` and `1/1` | `perpendicular_counts`, `extra_perpendicular_inf_transverse` (`Census.lean:74-110`) |
| ANISO oblique `D-2/D-2` and `1/0` | `oblique_counts`, `extra_oblique_inf_transverse` (`Census.lean:67-96`) |
| ANISO static `L`, `1/0` | `zero_counts`, `zero_variational_iff` |

`N2`/`N3` are `Module.finrank` of the kernel submodule and of its meet with `transverseSpace k` — genuinely dimensions, not displayed-basis counts, as the contract insists. "Elsewhere, frequencies outside the listed branches have no nonzero amplitude" is discharged family by family (`propagating_mode_iff`, `full_stationary_iff`, `div_propagating_iff`, `nonzero_root_iff`, `split_propagating_iff`, `parallel_kernel_iff`), each stated with `a ≠ 0` and nonzero frequency, so candidate roots and actual modes stay distinguished.

### 3.5 Anisotropic strata: exhaustiveness, disjointness, degeneracies

`q = perpSq e k ≥ 0` (`Geometry.lean:40`) and `perpSq_zero_iff` (`Geometry.lean:46`) give the parallel case; the residual `q > 0` splits on `k_e = 0`. The three cases are mutually exclusive and exhaustive on `k ≠ 0` (parallel ∧ perpendicular would force `K = 0`). This trichotomy is actually *exercised* inside Lean, e.g. `nonzero_root_positive` (`Certificate.lean:27-38`) branches on it. `frequency_coincidence_iff` (`Certificate.lean:113`) is an iff giving exactly `q=0` for the merger; `extra_exactly_transverse_iff` (`Certificate.lean:127`) is an iff giving exactly `p=0` for the extra mode's transversality — these are the two load-bearing exceptional-stratum claims and both are biconditional, not one-directional. No `GenericChart` hypothesis appears anywhere in the census (it appears only in the CAS basis bindings), matching the contract's explicit disclaimer.

Low-dimensional edges: `two_le_dimension` (`Geometry.lean:84`) derives `D ≥ 2` from `q ≠ 0`; at `D = 2` the ordinary split branch has `finrank = D-2 = 0` and `split_total_count` (`Census.lean:129-136`) still sums to `D-1` — the contract's caveat is correct. `dimension_one_no_propagating_mode` (`Checks.lean:77`) is correctly placed outside the contract domain.

### 3.6 Compact object identification and the limited D3 link

`Bindings.lean` proves `matrixB = (1/2)•referenceMatrix` and `matrixB_action : matrixB rho mu sigma (omega^2) k = (1/2)•actionMatrix .anisotropic 0 rho mu sigma 1 omega k` (`Bindings.lean:35-36,72-76,147-148,184-187`), with `routes : matrixA = (−2)•matrixB` (lines 67-70, 179-182). `referenceMatrix` (`CAS/Support.lean:65-67`) is literally `rho·z·W − mu(K I − k kᵀ)` with `W₀₀ = sigma`, and `referenceMatrix_action` ties it to the action-derived `actionMatrix`. So the identification is an action identity, not spectral coincidence, exactly as the contract claims — and the contract's caveat that a nonzero scalar preserves the kernel but not resolvent normalization or residues is the correct L2 caveat.

I independently corroborated the transcription at the source level: `PY.lean` declares input SHA-256 `10ebaecb…` = the in-packet `scripts/out/S10_anisotropic_strata_sympy_audit.out`. Line 32 of that transcript emits
`Q2_MATRIX_A = [[k2²μ_R + k3²μ_R − ω²ρ_br·s_rho, −k1k2μ_R, …], …]`,
and `PY.n13/n18/n20/n24` (`PY.lean:108-197`) evaluate to exactly those entries under `k1,k2,k3 → k 0, k 1, k 2`. Line 42's `Q2_MATRIX_B` is `−A/2`, consistent with `matrixB = (1/2)·referenceMatrix` and `A = (−1)·referenceMatrix`. The `M_A = −2 M_B` claim is therefore true of the actual emitted artifacts, not only of the Lean transcription.

The focused comparator's characterization is also accurate and appropriately hedged: `S10_anisotropic_strata_comparator.py:100-135` recomputes ranks from the emitted matrices, checks the stacked matrix, checks that the emitted roots exhaust the emitted determinant, and — importantly for L2 — checks `matrix·basis = 0 ∧ rank(basis) = n2 ∧ N7 = len(vectors) = n2`, i.e. a *complete* nullspace basis, not one residual-zero vector. Its outputs (`scripts/out/S10_anisotropic_strata_comparator.json`) agree numerically with the Lean census at D=3,4 for both exceptional directions (parallel `1/0` and `D-1/D-1`; perpendicular `1/0`, `D-2/D-2`, `1/1`; extra root `muR/(rhoBr·sRho)` per unit `K` = `E` at `p=0`). Coverage is honestly labelled "sampled".

### 3.7 Mutation controls

`_measurements/S10_lean_contract_check.py` never edits canonical sources (asserts `before == after`, line 119), requires the unmutated file to compile first, and records command, canonical hash, replacement text and full diagnostic. The recorded canonical hashes in `S10_lean_contract_checks.json` match this packet's MANIFEST hashes, so the run is on the reviewed revision.

I inspected every diagnostic. All seven rejections fail for the intended mathematical reason, not environmentally:

- `curl_normalization`, `fullgrad_replaced_by_divergence`, `divergence_sum_of_squares`: `ring` cannot close the altered quadratic-expansion identities (`lagrangian_increment`, `stiffnessPair_self`).
- `anisotropy_all_inertias`: unsolved `⊢ ∑ i, sigma·v i^2 = ∑ i, v i·v i + (sigma−1)·v e^2` inside `kinetic_eq` — the one-axis structure is genuinely load-bearing.
- `ignore_scalar_control`: the coefficient-scale definitional equalities collapse (`lagrangian c rho mu J` no longer defeq `S10Pilot.lagrangian rho (c·mu) J`), cascading through `variational_planeWave_iff`.
- `phase_average_missing_half`: `⊢ modalAction = 2·modalAction`.
- `negative_branch_reported_positive`: `⊢ False` from asserting the SIGNFLIP cone value is positive.

Six passing counterparts are present (five canonical originals plus `scalar_admissible_control`), matching the contract's "seven rejections and six passing controls".

The retained CAS suite (`S10_lean_cas_bridge_checks.json`) contains every mutation the contract names, each rejected with a substantive goal (`⊢ False`, `tauto` failures exhibiting the wrong stratum predicate, a `c 0 + c 1 = 0` vs `c 0 = 0` mismatch for the dependent-basis control), each with a paired positive control (`parallel_basis_control`, `root_multiplicity_control`, `generic_perpendicular_boundary_control`, `denominator_control`, …). Its recorded canonical hashes also match this packet.

**Vacuity.** I looked specifically for hypothesis sets that could be empty or self-defeating. They are not: `concrete_oblique_wave`, `concrete_perpendicular_wave`, `concrete_perpendicular_transverse_dimension`, `concrete_one_axis_kinetic` (`S10Anisotropic/Checks.lean:32-73`) are closed instances; `extra_wave_exists` and `ordinary_wave_exists` (with the necessary `3 ≤ D`) supply nonemptiness in general `D`; `positive_variational_certificate` and `negative_root_exists` (`Scalar.lean:165-221`) do the same for the coefficient family; and `positive_root_without_nonzero_wavevector` is precisely a hypothesis-necessity control showing `k ≠ 0` cannot be dropped.

**Build/axiom.** `CONTRACT_BUILD_VERIFICATION.txt`: exit 0, 3834 jobs, no warnings/admissions/custom axioms, observed axioms only `propext, Classical.choice, Quot.sound`, 30 axiom-free declarations, and the line "Proof build and mutation results do not constitute statement-fidelity review" — the correct L4 posture.

### 3.8 C1–C4

- **C1 (contract mapped to existing theorems): met.** Every anchor name in the evidence table exists at the cited namespace and states what the table says. The combinations the contract declines to re-prove (e.g. MAIN's static `N3 = 0`) are single rewrites off cited lemmas.
- **C2 (compact fidelity connection at the stated level): met.** The in-Lean action identity, the reviewed source correspondence for all six constructors, and the limited D3 imported-matrix link are all present, and the contract's boundary language ("not a claim that Lean executes or certifies either builder") matches what is actually proved.
- **C3 (controls with passing counterparts and substantive failures): met**, per §3.7.
- **C4 (two reviews resolved): partially.** This report discharges one leg; the second remains open, as the contract itself states. Build/axiom and provenance evidence is retained.

## 4. Findings

All editorial. None changes the action/operator identification, so none blocks the fidelity link.

**E1 — "first spatial component" is imprecise (`COVERAGE.md:113-114`).** The anisotropy in both builders sits on the first *displacement/velocity component* (`inertial_coefficients[0]` multiplying `(∂_t u₁)²`; wl:170,190-192), i.e. an amplitude index, not a coordinate direction — consistent with Lean's `kinetic e sigma (J 0)` weighting `v_e`. Minimal correction: "Both CAS engines distinguish their first displacement component (the amplitude index); this is Lean index `0`."

**E2 — `sigma`'s dimensionlessness is omitted from the parameter map (`COVERAGE.md:47-48`).** The contract declares `c` dimensionless but not `sigma`, although both builders declare it so (`py:466` `s_rho: ZERO_DIM`; `wl:204` `dimensionOf[sRho] == {0,0,0}`), and the Q5/Q6 row leans on it. Add "`sigma` is likewise declared dimensionless" to the ANISO domain sentence.

**E3 — the XCOEF_SCALE domain is narrower than what is proved (`COVERAGE.md:47`).** `c ≠ 1` is not needed for the table row: `cone_modeSpace`/`cone_counts` require only `c ≠ 0`, and the CAS assumes only `s > 0`. This understates coverage rather than overstating it. Either drop `c ≠ 1` from the domain or note that it is required only by `coefficient_changes_frequency`.

**E4 — two table cells lack a named anchor (`COVERAGE.md:61,64`).** MAIN's static `1 / 0` and DIVONLY's static `N3 = D-1` follow trivially but are not cited. Suggest adding `S10ScalarControls.zero_counts` (at `c = 1`) and `S10Controls.longitudinal_inf_transverse` to the MAIN row.

**E5 — the control instrument's rejection predicate scans the whole build output (`_measurements/S10_lean_contract_check.py:74-76`).** Two of the seven mutants (`anisotropy_all_inertias`, `ignore_scalar_control`) additionally emit unused-variable/unused-simp-arg linter errors under `-DwarningAsError=true`. I confirmed both still fail mathematically *inside the required declaration*, and the line→declaration mapping already drops the pre-theorem linter errors, so the controls are valid. Still, adding `set_option linter.unusedVariables false` / `linter.unusedSimpArgs false` to the scratch mutants, or scoping the diagnostic regex to the required declaration's error, would make the record unambiguous on its face.

**E6 — mutant text is not retained in the older CAS suite.** `S10_lean_cas_bridge_checks.json` records only `source_sha256` for each mutant, while the new instrument records the `replacement` text and the scratch `source`. The diagnostics do exhibit the false goal, so provenance is adequate, but the two suites are inconsistent in self-containedness.

## 5. Limitations and review coverage

**Reviewed in full:** `FORMALIZATION_POLICY.md`; the S10 work contract (COVERAGE.md:29-177); `S10Audit/Packages.lean`, `MatrixTrees.lean`, `Curl.lean`; `S10Pilot/{Action, PlaneWave, Spectrum, Variation, FiniteAction, PhaseAverage, Certificate, VariationalCertificate}.lean`; `S10Controls/{Action, Spectrum, Certificate, Scalar, PhaseAverage}.lean`; `S10Anisotropic/{Action, Geometry, Spectrum, Census, Certificate, Checks, PhaseAverage}.lean`; `S10Audit/CAS/Bindings.lean` and `RootBindings.lean`; `_measurements/S10_lean_contract_check.py` and its JSON result (every diagnostic); `scripts/S10_brane_mode_spectrum_sympy_audit.py` (`PACKAGES`, `stiffness_density`, `build_joint_assumptions`, `build_action`) and `mathematica/…audit.wl` (`buildPackage`); `scripts/S10_anisotropic_strata_comparator.py` and its result; both verification transcripts.

**Spot-checked, not exhaustively verified:** `S10Audit/CAS/PY.lean` node definitions (I traced `n0`–`n24` and the `matrixA/B/Residual` cells against the in-packet transcript); `S10_lean_cas_bridge_checks.json` (I read ~15 of 50 records in full and the full name/outcome list); `CAS_BRIDGE_VERIFICATION.txt` axiom listing (header plus a sample).

**Not reviewed (out of the stated scope, per the contract's exclusions):** the 916-expression bridge, minor/locus/metadata/rerun bindings beyond the roots and matrices above, the S9 packet, `scripts/S10_cross_engine_comparator.py`, the broad comparator/export pipeline, and ledger/paper reconciliation. Also not reviewed: whether the supplied physics is the intended physics — that remains a premise, as the contract states.

**Limitations of this review:**

1. **Generator not in packet.** `scripts/S10_lean_cas_bridge.py` (the `--check` provenance tool named in `PY.lean:1` and `CAS_BRIDGE_VERIFICATION.txt:21`) is not included, so the transcript → `PY.lean`/`WL.lean` transcription is not machine-verifiable from this packet. My spot-check of the D3 `Q2_MATRIX_A/B` cells against `scripts/out/S10_anisotropic_strata_sympy_audit.out:32,42` matched exactly, which is direct corroboration for the cells carrying the "stronger D3 link", but it is not a proof of the transcription as a whole. The contract and verification file both already disclose this ("Parsers/transcript translation remain tested software; generated proofs are kernel checked"), so I treat it as a disclosed limitation, not a finding.
2. **No re-execution.** I did not run `lake build`, the mutation instruments, the CAS engines, or the comparator; I read their recorded commands, hashes and diagnostics. The recorded canonical hashes match this packet's MANIFEST, so the evidence corresponds to the reviewed revision.
3. **No recomputation of `packet_sha256`** (no shell available here).
4. **Multiplicity claims.** `RootBindings.multiplicities` is D=3 only and carries explicit `GenericChart` / fixed-locus hypotheses (generic: all simple; parallel: `0` simple, merged root multiplicity 2). The contract's hedge — "certify multiplicities on their declared scope, not on every possible emitted package" — accurately describes what is there.

I found no statement in the contract that overstates what the sources prove, no parameter-map or normalization mismatch against either builder, no vacuous hypothesis set, and no mutation control that rejects for a non-mathematical reason.
