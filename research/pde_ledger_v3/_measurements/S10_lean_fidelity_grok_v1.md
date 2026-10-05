# S10 Lean contract — independent statement-fidelity review

**Reviewer:** Grok 4.6 (non-author). Author of the packet: Codex.

**No other reviewer findings were used.**

## 1. Packet identity

Pinned in `MANIFEST.json` (not independently recomputed; a prior hash command was cancelled):

| Field | Value |
|---|---|
| Revision | `S10 compact fidelity contract v1` |
| `packet_sha256` | `96d2c9d90a32b26d1c031e43ccc8580a8f8a102e18425491b231d4ac0cbc3b35` |
| Policy commit | `c2bdb663` |
| Proof checkpoint | `f21e6459` |
| Author | Codex |
| Contract file (pinned) | `lean/s10/COVERAGE.md` → `013cf4df7b9cfcb55881550dca3bc9d9effef417343d488c0cda917123c1c5db` |

Artifact under review: the **S10 work contract** in `lean/s10/COVERAGE.md` (lines 29–177), with `_measurements/S10_lean_contract_check.py` and `_measurements/S10_lean_contract_checks.json`. Kernel-build success in `lean/s10/CONTRACT_BUILD_VERIFICATION.txt` was treated only as a compilation record, not as fidelity evidence.

## 2. Verdict

**CLEAR** for the stated Lean scope (the compact six-family contract, C1–C3).

No statement-fidelity blocker was found. C4 (two independent non-author reviews of this fixed revision) is a process obligation and is not claimed complete by the contract; this report is one of those two legs. This is not clearance of the full S10 ledger, CAS production, comparator/export, or paper tasks.

## 3. Findings

No substantive findings. The classification table, action/operator identification, parameter/index maps, physical vs normalized frequency, N2/N3 counts, anisotropic case split, D=2 edge, compact CAS construction link, and limited D3 matrix identity match the cited definitions and hypotheses. Mutation diagnostics fail on the intended identities and have passing counterparts. Details of that check follow; they are confirmation, not a defect list.

### Object identification (not counts alone)

`S10Audit.packageAction_uses_stiffness` (`lean/s10/S10Audit/Packages.lean:38–45`) is a definitional identity

\[
\texttt{packageAction}=\frac{\rho}{2}\,\texttt{packageKinetic}-\frac{\texttt{stiffnessCoefficient}}{2}\,\texttt{packageStiffness},
\]

not a spectral inference. That matches the contract formula \(\rho/2\,v^TWv-\mu c/2\,S(J)\) (`COVERAGE.md:92–99`).

| Package | Lean constructor | CAS constructor |
|---|---|---|
| MAIN | `S10Pilot.lagrangian`: \(S=\frac12\sum_{ij}(\partial_i u_j-\partial_j u_i)^2\) (`S10Pilot/Action.lean:22–25`) | `stiffness_density("curl")` / `curlStiffness` with factor \(1/2\) (`sympy_audit.py:338–343`, `mathematica_audit.wl:154–159`) |
| FULLGRAD | `fullGradientStiffness=\sum_{ij}J_{i.succ\,j}^2` (`S10Controls/Action.lean:16`) | `kind=="fullgrad"` / `fullGradientStiffness` |
| DIVONLY | `divergenceOnlyStiffness=(\sum_i J_{i.succ\,i})^2` (`S10Controls/Action.lean:17–18`) | `divonly` / `Total[Diagonal[gradient]]^2` |
| SIGNFLIP | `S10ScalarControls.lagrangian (-1)` ⇒ \(c=-1\) (`Scalar.lean:14–15`, `Packages.lean:26`) | `stiffness_sign=1` / `stiffnessMultiplier=1` ⇒ \(T+\mu/2\,S\), same density |
| XCOEF_SCALE | `lagrangian scale` with \(c=\texttt{scale}\) (`Packages.lean:28`) | `s * mu_R` / `coefficientScale muR` |
| ANISO | `kinetic`: only \(W_{ee}=\sigma\), others \(1\) (`S10Anisotropic/Action.lean:13–17`) | first inertial slot \(s_\rho\rho\) (`sympy_audit.py:384–385`, `mathematica_audit.wl:167–172`) |

Jet convention: `J i.succ j=\partial_i u_j`, `J 0` is velocity (`COVERAGE.md:96–97`; `fieldJet` in `S10Pilot/PlaneWave.lean:20–21`; CAS `gradient[i][j]=\partial u_j/\partial x_i`). Parameter map \(\rho_{br}/\rho Br\mapsto\rho\), \(\mu_R/\mu R\mapsto\mu\), \(s_\rho/sRho\mapsto\sigma\), \(s/\texttt{coefficientScale}\mapsto c\) is used consistently. CAS distinguished axis = Lean `Fin` index `0`; Lean ANISO theorems allow arbitrary `e : Fin D`.

Plane-wave operator \(\rho z W-\mu c B(k)\) with \(B=KI-kk^T\) (curl), \(KI\) (FULLGRAD), \(kk^T\) (DIVONLY) matches:

- MAIN: `(rho*ω²-μK)•a+(μ k·a)•k` (`S10Pilot/Action.lean:111–113`)
- FULLGRAD: `(rho*ω²-μK)•a` (`S10Controls/Action.lean:98`)
- DIVONLY: `(rho*ω²)•a+(-μ k·a)•k` (`S10Controls/Action.lean:99`)
- ANISO: MAIN operator plus \((\rho\omega^2(\sigma-1)a_e)•e\) (`S10Anisotropic/Action.lean:69–70`)
- `stiffnessMatrixTree` (`S10Audit/MatrixTrees.lean:13–19`) is the same \(B\).

The extra \(1/2\) in curl \(S\) is load-bearing: `antisym_mode_pair` supplies a factor \(2\), so the modal operator has \(B=KI-kk^T\) without a leftover \(1/2\) (`S10Pilot/Action.lean:123–145`). Phase average is a further \(\langle\sin^2\rangle=1/2\) (`S10Pilot/PhaseAverage.lean:29–35` and the Controls/ANISO/Scalar copies). Stationarity is against all compact `TestField`s, not the ansatz (`S10Pilot/Variation.lean:188–220`; relative action, not a finite total action).

### Physical vs normalized squared frequency

Contract \(z\) is physical \(\omega^2\). ANISO Lean `z` in `normalizedOperator` is \(\rho\omega^2/\mu\) (`S10Anisotropic/Spectrum.lean:3,11–28`). The table uses physical \(R=(\mu/\rho)K\) and \(E=(\mu/\rho)(p^2+q/\sigma)\). That matches `extraConeValue=(μ/ρ)*extraValue` with `extraValue=(q+σp²)/σ` (`Certificate.lean:13–14`, `Geometry.lean:114–117`). For \(\mu\neq0\), kernels of the physical operator at \(R\) and \(E\) equal `modeSpace` at `normSq k` and `extraValue`. `split_variational_iff` / `split_variational_certificate` (`Certificate.lean:47–111`) state the census in physical \(\omega^2\).

### Classification, N2/N3, ANISO cases, edges

N2 = `Module.finrank` of the kernel; N3 = finrank of kernel \(\cap T\), \(T=\{a:k\cdot a=0\}\) (`COVERAGE.md:56–57`). Not a count of displayed basis vectors.

| Claim | Anchor | Hypotheses |
|---|---|---|
| MAIN: \(R\), \(T\), \(D-1/D-1\); static \(L\), \(1/0\) | `s10_variational_certificate`, `propagating_variational_mode_iff`; static N3 from `S10Controls.longitudinal_inf_transverse` / scalar `zero_counts` | \(\rho,\mu>0\), \(k\neq0\) |
| FULLGRAD: \(R\), whole space, \(D/D-1\); static \(0/0\) | `full_mode_counts`, `full_zero_space` | \(\rho\neq0\), \(\mu\neq0\), \(k\neq0\) |
| DIVONLY: \(R\), \(L\), \(1/0\); static \(T\), \(D-1/D-1\) | `div_mode_counts`, `div_zero_space` | same, plus \(\omega\neq0\) on the cone |
| SIGNFLIP: \(-R\), \(T\), \(D-1/D-1\); no real \(\omega\neq0\) | `signflip_counts`, `negative_control_no_real_wave`, `negative_root_exists` | \(c=-1<0\), \(D\ge2\) for a nonzero transverse vector |
| XCOEF_SCALE: \(cR\), \(T\), \(D-1/D-1\) | `coneValue`, `cone_counts`, `coefficient_changes_frequency` | \(c>0\), \(c\neq1\) for distinctness from MAIN |
| ANISO parallel \(q=0\): one root \(R=E\), \(O=T\), \(D-1/D-1\) | `parallel_counts`, `frequency_coincidence_iff`, `ordinarySpace_parallel` | \(\sigma\neq0\), \(k\neq0\), \(q=0\) (then \(p\neq0\)) |
| Perp \(p=0,q>0\): \(R\) on \(O\) \(D-2/D-2\); \(E\) on \(\mathrm{span}\{w\}\) \(1/1\) | `perpendicular_counts`, `extra_exactly_transverse_iff` | \(\sigma>0\), \(\sigma\neq1\) |
| Oblique \(p\neq0,q>0\): same ordinary; extra \(1/0\) | `oblique_counts`, `extra_oblique_inf_transverse` | same |
| Static ANISO: \(L\), \(1/0\) | `zero_counts`, `zero_frequency_iff` (kinetic correction vanishes at \(\omega=0\)) | \(k\neq0\) |

\(w=(q+\sigma p^2)e-\sigma p\,k\) is `extraVector` (`Geometry.lean:119–120`). Cases are exhaustive and disjoint on \(k\neq0\): `perpSq_nonneg`, `perpSq_zero_iff` (`Geometry.lean:40–57`); remainder splits on \(k_e=0\). No `GenericChart` hypothesis on the census. At \(D=2\) the ordinary split kernel has rank \(D-2=0\) (`Census.lean:127–128`; `ordinary_wave_exists` requires \(D\ge3\), `Checks.lean:19–21`). \(D=1\) theorems are outside the contract domain (`S10Pilot/EdgeCases.lean:36–40`; `S10Anisotropic/Checks.lean:77–93`). FULLGRAD at \(\omega=0\) is \(\bot\) (`full_zero_space`), so \(z=0\) is not a determinant root on this domain. Completeness for nonzero frequency is by iff theorems (`propagating_mode_iff`, `full_stationary_iff`, `div_propagating_iff`, `split_propagating_iff`, `parallel_kernel_iff`, `nonzero_root_iff`), not only the chosen positive root.

### Limited D3 matrix link

`S10Audit.CAS.PY/WL.matrixB_action` (`Bindings.lean:72–75,184–187`): imported route-B at \(\omega^2\) is \(\frac12\) of `actionMatrix .anisotropic 0`. `routes`: \(M_A=-2M_B\). `matrixB_kernel`: nonzero scalar, same kernel. `referenceMatrix` (`Support.lean:65–77`) is \(\rho z W-\mu(KI-kk^T)\) with \(e=0\). The contract correctly does **not** claim resolvent/residue equality. Focused comparator (`scripts/out/S10_anisotropic_strata_comparator.json`) is sampled D=3,4 parallel/perpendicular N2/N3; general coverage is the theorem census.

### Mutation controls and vacuity

Closure instrument (`S10_lean_contract_check.py`, result status `PASS`): 7 mathematical rejections, 6 passing controls, canonical SHA-256 unchanged. Failures are unsolved goals or type mismatch on the mutated identity, not missing imports/timeouts:

- `curl_normalization`: drop \(1/2\) in \(S\); `lagrangian_increment` / `modalAction_eq` unsolved (`S10_lean_contract_checks.json:38–47`).
- `fullgrad_replaced_by_divergence` / `divergence_sum_of_squares`: form swap; `stiffnessPair_self` unsolved.
- `anisotropy_all_inertias`: all inertias \(\sigma\); `kinetic_eq` unsolved (unused-`e` warning-as-error is a side effect; the identity still fails).
- `ignore_scalar_control`: drop \(c\); `lagrangian_eq` is not `rfl`.
- `phase_average_missing_half`: claims `modalAction = 2*modalAction`.
- `negative_branch_reported_positive`: `coneValue (-1) … > 0` reduces to `False`; passing twin `scalar_admissible_control` proves `< 0` and scale \(2\).

Named older controls in `S10_lean_cas_bridge_checks.json` match the table (`missing_parallel_branch`, `missing_perpendicular_branch`, `incomplete_parallel_kernel`, `generic_count_beyond_chart`, `collapsed_parallel_multiplicity`, `wrong_extra_root_sign`, `missing_chart_denominator`, …) with unsolved goals / failed `tauto` on false equivalences. Passing counterparts exist (`denominator_control`, `parallel_basis_control`, and concrete waves in `S10Anisotropic/Checks.lean:32–69`). Census hypotheses are not contradictory: `q\neq0` forces \(D\ge2\); `ordinary_wave_exists` uses \(D\ge3\); `negative_root_exists` uses \(D\ge2\).

### C1–C4 (non-review deliverables)

| Item | Status |
|---|---|
| **C1** contract mapped to existing theorems | Met. Cited names exist in the stated namespaces (table in §3). No new physics. |
| **C2** compact action/operator link at the stated level | Met. Construction-boundary identity to both CAS builders, plus the explicitly limited D3 `matrixB_action`/`routes` link. |
| **C3** essential controls with passing counterparts and substantive failure | Met. Closure instrument plus mapped older records. |
| **C4** two independent reviews resolved | **Open.** The contract does not claim it (`COVERAGE.md:160–163`). This report is one leg of that obligation for this scoped Lean contract, not ledger-wide clearance. |

## 4. Limitations and coverage

**Not independently verified:** packet `sha256` (quoted from `MANIFEST.json`); re-execution of `lake build` or the mutation script (diagnostics and hashes in the supplied JSON were inspected instead); any numeric CAS residual beyond reading constructors and the focused comparator’s reported N2/N3.

**Editorial, not blockers** (no contract change required for fidelity):

- Scalar SIGNFLIP/XCOEF_SCALE inherit EL/compact stationarity from MAIN via `lagrangian_eq` / `actionStationary_eq` (`Scalar.lean:17–41`), not under the identical names `actionStationary_iff_eulerLagrange` / `eulerLagrange_planeWave`. The evidence table already lists scalar `lagrangian_eq` and `variational_planeWave_iff` separately (`COVERAGE.md:140`).
- MAIN static N3 \(=0\) is not restated in `s10_variational_certificate`; it follows from kernel \(=L\) and \(L\cap T=\bot\). DIVONLY static N3 is implicit because the static space is \(T\) itself (`div_zero_space`).
- CAS `XCOEF_SCALE` joint assumptions are \(s>0\) without \(s\neq1\) (`sympy_audit.py:369–370`); the Lean distinctness theorem and SHARED PHYSICS use \(c\neq1\). That is a production-assumption difference, not a Lean identification error.

**Out of scope (as required):** exhaustive per-output CAS bridge; production Q7 / stratum handling / broad comparator-export / ledger-paper reconciliation; deriving the action or selecting physical \(D=3\); explicit exponentially growing solutions.

**Coverage of this review:** `FORMALIZATION_POLICY.md`; S10 work contract; closure instrument and result; `Packages.lean`, `ActionTrees.lean`, `MatrixTrees.lean`, `CAS/Bindings.lean`, `CAS/Support.lean`; all six family `Action`/`Spectrum`/`Certificate`/`Variation`/`PhaseAverage`/`PlaneWave` cores; ANISO `Geometry.lean`/`Census.lean`/`Checks.lean`/`Scaling.lean`; Scalar; EdgeCases; Specialization; Curl Q7; both CAS `build_action`/`stiffness_density`/`buildPackage` and `PACKAGES`; focused comparator JSON (sampled counts); mapped mutation records; SHARED PHYSICS action/ansatz; `CONTRACT_BUILD_VERIFICATION.txt`. Historical `*_RESULT.md` files were not treated as extra acceptance criteria.
