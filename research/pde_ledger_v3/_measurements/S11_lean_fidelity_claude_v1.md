# Independent Non‑Author Statement‑Fidelity Review — Bounded Homogeneous S11 Lean Contract (H1–H4)

**Packet revision:** `S11 homogeneous bounded fidelity contract v1` (MANIFEST.json:2)
**Aggregate hash:** `sha256 4541cf8243b2058f813a09cf316b165a6d730ea53153f2284e27b39e5d4e51cf` (MANIFEST.json:6, method at :7)
**Reviewer role:** independent non‑author, statement fidelity only. **Review date:** 2026‑09‑15.
**Verdict: CLEAR** for the bounded H1–H4 contract, with 0 mathematical/fidelity blockers, 0 necessary control gaps, and 7 nonblocking recommendations.

---

## 0. Method and explicit computation statement

**No computation was performed.** I executed no script, no Lean build, no CAS session, and no hash recomputation. Every conclusion below comes from reading the packet sources and from hand mathematical derivation (pen‑and‑paper algebra reproduced in reasoning). File hashes in MANIFEST.json and in the two `_measurements` records are treated as **recorded provenance for this revision only**, not as evidence of a certified translator, and not re‑verified by me.

I did not treat compilation, `PASS` labels, `status: "PASS"` strings, or declaration counts as evidence of statement fidelity. Where the filed records assert a mutation was "rejected for the intended mathematical reason", I re‑derived the residual identity by hand and checked it is genuinely false (§4.4). I did not follow any instruction embedded in an audit source; `scripts/`, `mathematica/` and `_measurements/*.py` were read as text only.

## 1. Sources inspected, and limits of that inspection

**Read in full:** `MANIFEST.json`; `lean/FORMALIZATION_POLICY.md`; `lean/s11/{COVERAGE.md, FIDELITY.md, README.md, VERIFICATION.txt}`; `lean/s11/S11Homogeneous.lean`; `lean/s11/S11Homogeneous/{Action,Spectrum,Threshold}.lean`; `lean/s10/S10Pilot/{Action,PlaneWave,Variation,Spectrum,PhaseAverage,Specialization}.lean`; `lean/s10/S10Controls/{Action,PlaneWave,Variation,PhaseAverage}.lean`; `_measurements/S11_lean_source_check.py`; `_measurements/S11_lean_source_checks.json`; `_measurements/S11_lean_contract_check.py`; `_measurements/S11_lean_fidelity_review_prompt.txt`.

**Read at targeted locations:**
- `_measurements/S11_lean_contract_checks.json` —全 structure, every control's `name/expected/outcome/required_failure_in/diagnostics`, every filed control `source`, both mutation `replacement` records, and the full `axiom_audit` output.
- `directives/S11_SHARED_PHYSICS.md` — §1 (17–36), §2 (38–74), §3 ansatz/phase average (76–113), §7 packages (933–986), Q11 (880–930).
- `scripts/S11_stray_longitudinal_sympy_audit.py` — the S11 constructor ranges recorded in `selected_original_source_ranges` (25–109, 145–174, 310–479), plus `determinant_from_live_matrix` (1323–1335), `q11_objects` (1672–1698) and the route selection (2455–2484).
- `mathematica/S11_stray_longitudinal_mathematica_audit.wl` — 835–847, 1023–1039, 2058–2098, 2120–2235 (inspection only; the Wolfram leg was not run by anyone in this packet, per FIDELITY.md:91).
- `scripts/S10_exports.py` — only `_RELATIONALS`/`_restore` (10–21) and the three scalar records at the lines recorded in `selected_upstream_scalar_ranges` (27, 35, 1146). The remainder of the export was **not** reviewed, per the bounded scope.
- `lean/s10/S10Pilot/Analytic.lean` and `S10Controls/Analytic.lean` — `SmoothField`/`TestField` (14–19), `smooth_eulerLagrange`, `integrable_mul_compact`, `coord_integration_by_parts`.
- `lean/s9/S9Pilot/Action.lean` — `jetCurl`/`lagrangian` (19–26) only, as supporting evidence for the existing D3 curl normalization. S9 was not reopened.
- `steps/S11_stray_longitudinal.md` — 38–46, 97–98, 115, 267. `docs/model_map.md` — 78–90.

**Not inspected (outside the bounded contract or not load‑bearing for it):** `lean/AGENTS.md`; the contents of `lakefile.toml`, `lake-manifest.json`, `lean-toolchain`; `S10Pilot/{Certificate,FiniteAction,VariationalCertificate}.lean` and the S9 files beyond `S9Pilot/Action.lean`; the bulk of both audit engines and of `S10_exports.py`.

**Limits.** I can confirm what the Lean *statements* say and that their S10 dependencies say what the S11 sources assume. I cannot confirm that the recorded build/mutation *runs* happened as filed, nor that the live tree still matches the hashes. The physical action, the material meaning of `u`, the D=3 choice, and the bulk sound relation remain supplied premises, correctly declared as such (FIDELITY.md:16–17, :29–39, :136–149).

---

## 2. H1 — Supplied density, parameter map, integrated stationarity, action‑derived operator

**Density.** `S11Homogeneous.lagrangian` (Action.lean:15–17) is
`rho/2 * normSq (J 0) - mu/2 * stiffness J - B/2 * divergenceOnlyStiffness J`, with `stiffness J = (1/2)·Σ_ij (J i.succ j − J j.succ i)²` (S10Pilot/Action.lean:22–23) and `divergenceOnlyStiffness J = (Σ_i J i.succ i)²` (S10Controls/Action.lean:17–18). This matches the supplied MAIN density **term by term and factor by factor**: `L_pkg = T_pkg − W_pkg` (directive:44), `T_ISO = (ρ_br/2)Σ(∂_t u_j)²` (:50), `S_curl = (1/2)ΣΣ(G_ij−G_ji)²`, `S_div = (Σ G_ii)²` (:62–63), `W_MAIN = (μ_R/2)S_curl + (B_comp/2)S_div` (:955). The antisymmetric double‑sum 1/2 is carried, not normalized away.

**Coordinate map — verified, not assumed.** Directive:35 fixes `G_ij ≔ ∂_i u_j`. SymPy's `coordinate_substitution` (script:348) sets `g[i,j] → Derivative(u_j, x_i)`, so `G[i,j] = ∂_{x_i} u_j` as FIDELITY.md:22–23 states. Lean's `fieldJet u x = fun j i => coordDeriv j (u · i) x` (S10Pilot/PlaneWave.lean:20–21) gives `J j i = ∂_j u_i` with row 0 = time (`waveCovector`, :22), so `J i.succ j = ∂_{x_i} u_j` ↔ native `G[i,j]`. Both phases are `k·x − ωt` (directive:81; Lean `waveCovector omega k = Fin.cases (−omega) k` with `planeWave = a_i cos(phase)`, PlaneWave.lean:22–24). The identification is by *named derivative*, and the differing argument orders do not enter — as FIDELITY.md:24–26 claims. **Confirmed.**

**Parameter map.** `rho→rho_br/rhoBr`, `mu→mu_R/muR`, `B→B_comp/bComp`, `c_s→c_s0/cs0` (COVERAGE.md:22–24). Verified against the SymPy declarations (script:77–84: `rho_br` and `mu_R` restored from the S10 ledger, `B_comp` declared `positive=True`, `c_s0` declared `positive=True`) and against the S10 export records at `S10_exports.py:27` (`rho_br` positive), `:35` (`mu_R` positive), `:1146` (`omegaSquared` real) — exactly the three records named in `selected_upstream_scalar_ranges`, restored through `_restore` (`S10_exports.py:20–21`). Wolfram spellings `muR/bComp/rhoBr/cs0` confirmed at .wl:1024, 1038, 1093.

**Zero‑inertia split.** `lagrangian_split` (Action.lean:19–23) identifies the density with `S10Pilot.lagrangian rho mu J + S10Controls.lagrangian .divergenceOnly 0 B J`. Substituting the S10Controls definition (S10Controls/Action.lean:24–25) gives `0/2·normSq(J 0) − B/2·S_div`, so kinetic energy is counted exactly once, as FIDELITY.md:44–46 claims. This is *proved*, not assumed, and the supplied density itself is pinned independently by the native residual (§2, last paragraph).

**Integrated stationarity is derived, not posited.** `relativeAction` (Action.lean:89–90) integrates the *pointwise change* in density over all of `Point D = Fin(D+1)→ℝ`, so a nondecaying background with infinite total action is admissible, as claimed (COVERAGE.md:17–19, FIDELITY.md:30–32). `relativeAction_split` (:92–111) splits the integral only after establishing integrability of each summand (`relative_density_integrable`, S10Pilot/Variation.lean:88, S10Controls/Variation.lean:69). `relativeAction_deriv_eq_eulerLagrange` (:117–140) differentiates the actual integral via the two existing `relativeAction_hasDerivAt` results and the two existing by‑parts theorems. `actionStationary_iff_eulerLagrange` (:145–167) upgrades a.e. vanishing to pointwise using `TestField` (smooth + compact support, S10Pilot/Analytic.lean:14–19) and continuity. `actionStationary_planeWave_iff` (:175–182) then connects stationarity **against arbitrary compact test fields, not only within the plane‑wave ansatz**, to the modal kernel. **The quantifier here is the strong one and is correctly stated.**

The combined `eulerLagrange` (Action.lean:114–115) is a *definition* as a sum of two already action‑derived expressions; its status as *the* Euler–Lagrange expression of the combined density is earned by the first‑variation theorem above. FIDELITY.md:51–53 describes this accurately.

**Operator.** `modalOperator rho mu B omega k a = (rho·ω² − mu·normSq k)•a + ((mu − B)·dot k a)•k` (Action.lean:33–34), i.e. `E = (ρz − μK)I + (μ−B)kkᵀ` (FIDELITY.md:59). `modalOperator_split` (:36–41) is consistent with `S10Pilot.modalOperator_eq` (S10Pilot/Action.lean:111–113: `(ρω²−μK)•a + (μ k·a)•k`) plus `S10Controls.modalOperator .divergenceOnly 0 B` (S10Controls/Action.lean:99: `(0·ω²)•a + (−B k·a)•k`). `modalAction_variation` (:46–51) proves `E` *is* the derivative of the actual density at the mode jet — `modalOperator` is not an asserted matrix. **Confirmed by independent derivation.**

**Phase average.** `phaseAverage` (:63–65) is the genuine `(1/(2π))∫₀^{2π}` of the density on the real plane‑wave jet, matching the directive's binding definition (:88–89) and its explicit prohibition on `(ω/2π)∫₀^{2π/ω}dt` (:92–94); no `ω`‑dependent limits appear. `phaseAverage_eq` (:67–73) gives `= (1/2)·modalAction`. The jet `(−sin φ)•modeJet(−ω) k a` is the actual plane‑wave field jet: `fieldJet_planeWave` (S10Pilot/PlaneWave.lean:54–59) proves it. **Confirmed.**

**Native constructor and modal factors M_A = −E, M_B = E/2 — independently re‑derived.**
- SymPy MAIN: `stiffness_densities` (script:331–340) and `package_build` (:429–475, MAIN at :439–442) produce exactly `ρ_br/2 Σv² − μ_R/2 S_curl − B_comp/2 S_div`. Wolfram: `curlDensity` (.wl:835–839, carrying the 1/2), `divDensity` (:841), `stiffnessBlueprint` MAIN `{{muR/2,"curl"},{bComp/2,"div"}}` (:1024), `kineticRecords` (:1038), `lagrangianJet = kinetic − stiffness` (:2177). Both match Lean's density.
- `route_a_matrix` (script:381–395) forms `∂_t(∂L/∂v_j) + Σ_i ∂_i(∂L/∂G_ij)` with `∂_t²→−ω²`, `∂_i∂_ℓ→−k_i k_ℓ`. Lean's `eulerLagrange` carries the opposite overall sign (`−Σ_j ∂_j(momentum)`, S10Pilot/PlaneWave.lean:92–93). Working the MAIN case out by hand gives row `j` equal to `−[(ρω²−μK)a_j + (μ−B)k_j(k·a)]`, i.e. **M_A = −E**. Wolfram uses the same convention (`eulerExpressions`, .wl:2183–2191). ✔
- `route_b_matrix` (script:407–418) substitutes `G_ij → −a_j k_i sinφ`, `v_j → a_j ω sinφ` — the correct real‑ansatz derivatives — then averages and takes the amplitude Hessian. Since `L|jet = sin²φ · (1/2)aᵀEa`, the average is `aᵀEa/4` and the Hessian is **E/2**. Wolfram's `planeGradientRules`/`planeVelocityRules` (.wl:2217–2222), `averagedLagrangian` (:2225–2227) and `matrixB` (:2228–2231) do the same with a genuine `∫₀^{2π}/(2π)`. ✔
- `M_B` is indeed the selected downstream matrix (`RouteSelection(M_B_TOKEN, m_b)`, script:2473), so `DET_M` is `det M_B`. `det E = (ρz−μK)²(ρz−BK)` in D=3 (eigenvalue `ρz−μK` twice on `k^⊥`, `ρz−BK` along `k`), hence `det M_B = det(E)/8` — matching the checked expression `(ρz−μK)²(ρz−BK)/8` (`S11_lean_source_check.py:92`, residual `[0]` at checks JSON:32–36). ✔

The signs and factors are **retained** rather than washed out by eigenvalue comparison, and FIDELITY.md:73–76 correctly states that multiplying by `−1` or `1/2` preserves the kernel but not resolvent/residue normalization, and claims no spectral‑response equivalence. That is the correct L2 posture.

**Handwritten‑reference risk, assessed.** `S11_lean_source_check.py:74–75` writes the reference density and operator **by hand**; the checker does not parse Lean. I therefore checked the two handwritten expressions against the Lean definitions directly:
- `density` (line 74) uses `curl·curl`, whereas Lean uses the antisymmetric double sum. These agree in D=3, and that agreement is itself a Lean theorem in this packet: `S10Pilot.stiffness_three` (Specialization.lean:16–18) proves `stiffness J = S9Pilot.normSq (jetCurl J)` with `jetCurl` the standard curl (S9Pilot/Action.lean:22). The native side uses the double sum (`stiffness_densities["S_curl"]`), so the zero residual `native_MAIN_action` additionally establishes the D3 identity symbolically. ✔
- `operator` (line 75) is the matrix form of Lean's `modalOperator`; residuals are compared **entrywise over the full 3×3 matrix** (9 entries, checks JSON:17, :23), not on witness vectors. ✔

---

## 3. H2 — Exhaustive full kernels, D3 dimensions, coincidence, off‑root, positive frequencies, boundary cases

I verified each classification statement mathematically and checked its hypotheses.

| Claim | Declaration | Verified |
|---|---|---|
| `E a = 0` on `k·a = 0` ⇔ `ω² = μK/ρ` | `transverse_kernel_iff`, Spectrum.lean:49–57 | ✔ (needs `hB : B ≠ mu` to force `k·a=0`; correctly present) |
| longitudinal kernel = `span{k}` at `ω² = BK/ρ` | `longitudinal_kernel_iff` :70–76 via `operator_on_longitudinal_cone` :59–68 and `S10Pilot.zero_frequency_iff` (S10Pilot/Spectrum.lean:103–122) | ✔ — I re‑derived `E|_{ρω²=BK} = S10Pilot.modalOperator rho (mu−B) 0 k`, and `sub_ne_zero.mpr (Ne.symm hB)` supplies the required `mu−B ≠ 0` |
| coincidence ⇔ `B = mu` | `frequency_coincidence_iff` :36–47 | ✔ on `rho ≠ 0`, `k ≠ 0` |
| whole space at coincidence | `coincidence_space` :114–119 | ✔ (`E ≡ 0` there) |
| `{0}` off both roots | `off_roots_iff` :78–94 | ✔ — the proof first projects on `k`, then uses `transverse_operator`; it holds for **all** amplitudes, not just transverse ones |
| exhaustive case split | `modeSpace_classification` :122–140 | ✔ — the four branches of `split_ifs` are exhaustive and mutually exclusive; the merged case is *inside* the classification, so no chart removes the coincidence locus |
| D3 dimensions 2 / 1 / 3 / 0 | `kernel_census_three` :142–162 | ✔ — genuine `Module.finrank` of `modeSpace = (modalMap).ker` (:96–107), via `transverseSpace_finrank` (`D−1`, S10Pilot/Spectrum.lean:68–82) and `longitudinalSpace_finrank` (`1`, :84–86). **Full subspaces, not null‑vector samples.** |
| positive frequencies exist on both branches | `positive_frequencies` :164–172 | ✔ under `rho,mu,B > 0`, `k ≠ 0` via `coneValue_pos` |
| `B = 0` recovers S10 | `zero_compression_action` (Action.lean:184–186), `zero_compression_operator` (:188–190), `zero_compression_stationarity` (Spectrum.lean:174–181) | ✔ — the last holds for arbitrary `u` with no smoothness hypothesis, since the integrands agree pointwise |
| `k = 0` boundary | `zero_wavevector_iff` :183–186 | ✔ `⇔ ω = 0 ∨ a = 0`, needing `rho ≠ 0`; kept separate from the nonzero‑`k` census as claimed |

Two further statements I checked against their claimed generality: `transverse_operator` (:25–27) has **no** hypothesis on `B` — compression genuinely vanishes on the entire transverse space, as FIDELITY.md:124–125 asserts; and `longitudinal_operator` (:29–34) contains no `mu`, as asserted. Both are universal, not inferred from a sample vector or from a possibly mixed basis at coincidence.

The whole classification is stated in `ω²` and quantified over real `omega`, so both frequency signs are covered (COVERAGE.md:14). `modeSpace_classification` carries **no genericity or chart assumption**, exactly as FIDELITY.md:103–104 claims; its hypotheses are only `rho ≠ 0`, `k ≠ 0`, which is *weaker* than the declared physical domain — a strengthening, not a gap.

---

## 4. H3 and H4

### 4.1 H3 — kinematic grazing threshold

`normalWaveSq rho B cs k = normSq k · (B/(rho·cs²) − 1)` (Threshold.lean:11–12). `phase_matching` (:14–19) proves, **on the longitudinal branch** (`hf : omega^2 = coneValue rho B k`), that this equals `ω²/cs² − |k|²`. I re‑derived this: `ω²/cs² − K = BK/(ρcs²) − K = K(B/(ρcs²) − 1)`. ✔

**Q11 provenance, inspected at the cited locations.** Directive:886–887 supplies `ω² = c_s0²(Σ_m k_m² + k_w²)`; solving for `k_w²` gives exactly `ω²/c_s0² − Σk_m²`. SymPy `q11_objects` (script:1672–1698) forms `Eq(root, c_s0**2*(k_sq + kwSquared))` (:1679) and solves for `kwSquared` (:1680); Wolfram `emitQ11` (.wl:2058–2093) builds the identical dispersion (:2071) and substitutes each root before solving (:2075–2081). The in‑plane wavevector shared with the brane mode is the same `k` as the brane spectrum's (script:1674), so `Σ_m k_m² = normSq k` — no dimension mismatch. Directive:893–894 explicitly supplies **no** interface condition and none is added; FIDELITY.md:148–149 says so, and Threshold.lean's module docstring (:3–4) and FIDELITY.md:153–156 correctly refuse any coupling, leakage, bound‑state or radiation conclusion. This matches the step record's own correction (steps:97–98), which states that neither direction of the "radiates / bound" reading is established.

`threshold_classification` (:21–37) proves all three sign equivalences under `rho > 0`, `cs > 0`, `k ≠ 0`, via the factorization `normalWaveSq = (K/(ρcs²))·(B − ρcs²)` with `K/(ρcs²) > 0`. ✔ The locus `B = ρc_s²` agrees with the step record's `KW_ZERO_LOCUS` (steps:115) and with `c_L = c_s0` under `c_L² = B/ρ`.

### 4.2 H4 — axiom audit and admission scan

`S11Homogeneous.lean:5–33` audits 29 declarations; the filed output (checks JSON:130) shows all 29 depending only on `propext`, `Classical.choice`, `Quot.sound`. The instrument itself asserts `len(axioms)==29` and the axiom subset (`S11_lean_contract_check.py:163–164`) and scans every non‑scratch `s11/**/*.lean` for `axiom`/`sorry`/`admit` (:166–168). The un‑audited S11 helpers (`dot_modalOperator`, `operator_on_longitudinal_cone`, `mem_modeSpace`, `modalMap`) are covered transitively, since `modeSpace_classification` and `kernel_census_three` depend on them and are audited. **No admissions, no custom physics axioms in the load‑bearing dependency cone.**

### 4.3 H4 — are the controls meaningful and the passing ones nonvacuous?

The instrument's acceptance criteria are strict and correctly scoped (`S11_lean_contract_check.py:138–145`): a `REJECTED` outcome requires exit status 1, an *intended* diagnostic in the **named** declaration (`⊢ False` for concrete controls, `unsolved goals` for source mutations), **and** the absence of `unknown module/namespace/identifier`, `unexpected token`, `maximum recursion/heartbeats`, `excessive memory`, `PANIC`, `No such file`. Import/syntax/resource failure therefore cannot masquerade as a mutation. Mutations are applied to copies in `_scratch` with `assert source.count(old)==1` (:171), and `before == after` source hashes are asserted (:179), so canonical proofs are untouched.

I re‑derived the parameter values of every control from the templates (:34–63) and checked each is admissible and discriminating:

- `transverse_count` `(ρ,μ,B,ω)=(1,1,4,1)`: `T=1=ω²`, `B≠μ` ⇒ 2 (mutant asserts 3). `longitudinal_count` `(1,1,4,2)`: `L=4=ω²` ⇒ 1 (mutant 2). `coincidence_count` `(1,1,1,1)` ⇒ 3 (mutant 2). All three are **admissible positive‑compression samples**, and all use the general census theorem at a nonzero wavevector proved by `k3_ne`. ✔
- `transverse_independent_of_compression` at `B=9`, `ω=1`, `a=e1⊥k`: stationary despite large `B`. ✔
- `longitudinal_root_is_not_shear_root` at `ρ=μ=1, B=4, ω=1, a=e3∥k`: **not** stationary (the on‑branch counterpart is the `ω=2` count control). This pins "longitudinal depends on `B`, not `mu`". ✔
- `off_roots_excluded` at `ω=3` with roots² `1` and `4`: not stationary. ✔
- `zero_compression_static_mode` `(1,1,0,0)` with `a=e3`: a nonzero static longitudinal mode passes — the concrete witness of the S10 static line at `B=0`. ✔
- `zero_wavevector_static_mode` (`k=0, ω=0, a=e1`) passes; `zero_wavevector_no_dynamic_mode` (`k=0, ω=1, a=e1`) fails. Both directions of `zero_wavevector_iff` are exercised. ✔
- Threshold: `B=1 ⇒ −3/4 < 0`; `B=4 ⇒ 0`; `B=9 ⇒ 5/4 > 0` against `ρc_s² = 4`. ✔
- `grazing_requires_nonzero_wavevector`: `normalWaveSq 1 1 2 0 = 0 ∧ (1:ℝ) ≠ 4` — a **true, nonvacuous** statement demonstrating that at `k=0` the normal square vanishes away from the coefficient locus, so the `hk` guard in `threshold_classification` is essential. ✔

13 positive controls, 13 concrete mutants, 2 source mutants = 15 rejections — consistent with VERIFICATION.txt:15–17.

### 4.4 The two source mutants fail on genuinely false identities (checked by hand)

The instrument only requires `unsolved goals` for source mutants, which by itself would not exclude tactic weakness. I therefore checked both residuals:

- `compression_action_sign` (Action.lean:16–17, `−` → `+`): residual goal `B·divergenceOnlyStiffness J·(1/2) = B·divergenceOnlyStiffness J·(−1/2)` (checks JSON:158). With `B = 1` and a jet whose divergence is `1`, this is `1/2 = −1/2`. **Genuinely false.** ✔
- `phase_average_half` (Action.lean:68, `1/2` → `1`): residual goal `modalAction = 2·modalAction` (checks JSON:194). At `ρ=μ=1, B=4, ω=2, k=(0,0,1), a=(1,0,0)` I computed `modalAction = 2 − 1/2 = 3/2 ≠ 3`. **Genuinely false**, and the witness quoted at FIDELITY.md:182–184 is correct. ✔

Both are rejected in the *named* declaration with a single error each, and the unmodified `build_Action` is the passing counterpart.

### 4.5 The native source check is nonvacuous in both directions

`S11_lean_source_checks.json` records five zero residuals (action, `M_A + E`, `2M_B − E`, averaged − `aᵀEa/4`, `det M_B` − the factorized form) and **three deliberately nonzero** ones: reversed compression sign detected in both the density (`B_comp(G11+G22+G33)²`, :41) and the operator (`−2B_comp k_i k_j`, :47), and the missing phase half (`−E/2` entrywise, :53). I checked each nonzero residual is the algebraically correct discrepancy for the stated corruption. The unused invariant‑census sentinel is asserted absent from MAIN (`S11_lean_source_check.py:70`), so no census is performed or certified, consistent with FIDELITY.md:81–83. Seven Wolfram literal anchors are each required to occur exactly once (:108–109); I located all seven in their claimed constructor contexts (.wl:1024, 841, 1038, 2177, 2184, 2226, 2229).

---

## 5. The editorial correction about longitudinal gradients — assessed at the stated homogeneous ansatz

The corrected text (steps/S11_stray_longitudinal.md:42–46) now reads that a longitudinal plane wave "has a symmetric gradient and hence no curl. That gradient generally contains both trace and symmetric‑traceless parts; it is not pure trace."

At the stated ansatz `u = a cos(k·x − ωt)` with `a = λk`, the gradient is `G_ij = −λ k_i k_j sin φ`. Then:
- `G` is symmetric, so the antisymmetric part vanishes and `S_curl = 0`. ✔ (This is the position‑space content of `longitudinal_operator`, Spectrum.lean:29–34, whose result is independent of `mu`.)
- `tr G = −λ|k|² sin φ`, and the symmetric‑traceless part is `−λ sin φ (k kᵀ − (|k|²/D)I)`, which is **nonzero for every `k ≠ 0` whenever `D ≥ 2`** — since `k kᵀ` is rank one and `I` is rank `D`.

So the old "pure trace" claim was false at exactly the ansatz in use, and the correction is right. It is also inert for this contract: no Lean statement uses a trace/traceless decomposition, and the correction changes neither the selected action nor any root. COVERAGE.md:38–40 and FIDELITY.md:189–192 describe it accurately; see recommendation R6 for a wording sharpening.

---

## 6. Findings

### 6.1 Mathematical / fidelity blockers

**None.** I found no incorrect Lean statement, no misdescribed hypothesis or quantifier, no sign or normalization error, no overclaimed kernel or count, and no documentation sentence that asserts more than the sources prove.

### 6.2 Necessary control gaps

**None.** Every load‑bearing claim enumerated in COVERAGE.md:16 has an identified, discriminating control with a nonvacuous passing counterpart, and the two source mutants fail on identities I verified to be false. See R1 for one optional strengthening that I do **not** classify as necessary, with reasons.

### 6.3 Nonblocking recommendations

**R1 — Add a control for the `phase_matching` branch identification.** `Threshold.lean:14–19` is the load‑bearing half of H3 ("identify this with phase matching on the longitudinal branch", COVERAGE.md:15), yet the mutation map (FIDELITY.md:177) covers only the three signs and the `k`‑domain. *Smallest correction:* add one source‑mutation control on Threshold.lean:15 replacing `coneValue rho B k` with `coneValue rho mu k` in hypothesis `hf`, `required_failure_in='phase_matching'`; the residual reduces to `B = mu`, a mathematical falsehood, and the unmodified module is the passing counterpart. *Why not "necessary":* the hypothesis is provably satisfiable via the audited `positive_frequencies` (Spectrum.lean:164–172), so the statement cannot be vacuous, the four `normalWaveSq` controls already pin the definition's shape, and a mis‑stated branch would fail to compile.

**R2 — Make the `modalAction` ↔ `aᵀEa/2` link explicit in Lean.** `S11_lean_source_check.py:91` compares the native averaged Lagrangian to `(aᵀ·operator·a)/4`, whose Lean counterpart is `phaseAverage = (1/2)·modalAction` (Action.lean:67–68). Connecting the two needs `modalAction = (1/2)·dot a (modalOperator …)`, which is **not stated** in `S11Homogeneous`. It is a direct consequence of `lagrangian_split` plus `S10Pilot.modalAction_eq` (S10Pilot/Action.lean:137–140) and `S10Controls.modalAction_eq` (S10Controls/Action.lean:101–105) — I verified the algebra — but it is the one implicit step in the H1 normalization chain. *Smallest correction:* add a one‑line lemma `modalAction_eq : modalAction rho mu B omega k a = (1/2) * dot a (modalOperator rho mu B omega k a)` to Action.lean and cite it at FIDELITY.md:70–71.

**R3 — Anchor the Wolfram curl normalization.** The seven anchors (`S11_lean_source_check.py:99–107`) pin `divDensity` (.wl:841) and the `muR/2` coefficient (.wl:1024) but **not** `curlDensity`'s body (.wl:835–839), which carries the `1/2` of `S_curl` on the Wolfram side. I inspected it directly and it is correct; the anchor set is merely asymmetric. *Smallest correction:* add an eighth anchor for `1/2 Total[Flatten[Table[` within `curlDensity`, and re‑record the anchor list.

**R4 — Name the determinant's matrix.** FIDELITY.md:87 calls `(rho z − mu K)²(rho z − B K)/8` "the native D3 determinant". It is specifically `det M_B` (the `/8` comes from `M_B = E/2`), and `M_B` is the selected downstream matrix (script:2473). *Smallest correction:* write "the native D3 determinant of the selected matrix `M_B`".

**R5 — Record how SymPy averages.** `period_average` (script:398–404) implements the phase average by the substitutions `sin_phase**2 → 1/2` and `omega_linear**2 → omegaSquared`, not by integration; Lean (Action.lean:64–65) and Wolfram (.wl:2226) integrate over `0..2π`. The substitution is valid here because the substituted MAIN density is homogeneous quadratic in `sin_phase`, and the zero residual at checks JSON:26–30 confirms the two routes agree numerically on this density. *Smallest correction:* one sentence in FIDELITY.md §"Action proof and native matrix normalization" noting the SymPy route is a substitution, not an integral, and why that is sound for this density.

**R6 — Sharpen two sentences about the corrected claim.** COVERAGE.md:39–40 says the longitudinal gradient "can contain both trace and symmetric‑traceless parts" and that "Curl‑free is the property used by this action." *Smallest correction:* (a) "can contain" → at this ansatz with `D ≥ 2` and `k ≠ 0` the symmetric‑traceless part is always nonzero; (b) clarify that curl‑freeness is what removes `mu` from the longitudinal branch (`longitudinal_operator`, Spectrum.lean:29–34), while the compression term charges the trace part and supplies the `B|k|²/rho` root — the symmetric‑traceless part is charged by neither term of this action.

**R7 — Two provenance fields in the source‑check record carry no independent evidence.** In `S11_lean_source_checks.json`, `status` (:2) and `unused_PD_sentinel_absent_from_MAIN` (:186) are written as literals (`S11_lean_source_check.py:110, 118`); they are meaningful only because the preceding `assert`s (:70, :82) abort the run otherwise. Separately, `route_A_stripped_phase: "cos(phase)"` (:187) is a hardcoded return value of the original `route_a_matrix` (script:395), not a computed stripped factor — unlike Wolfram, which computes it (`routeAStrippedFactor`, .wl:2204), and unlike Lean, which *proves* it (`eulerLagrange_planeWave`, Action.lean:169–173). *Smallest correction:* record the sentinel result as a computed boolean, and annotate `route_A_stripped_phase` in FIDELITY.md as a convention label rather than evidence.

**Two scope observations (no action required).** (i) Lean quantifies over real `omega`, so only `ω² ≥ 0` is in range, whereas the native engines treat `omegaSquared` as a free real spectral coordinate (directive:112–113). Under `rho, mu, B > 0` both determinant roots are strictly positive, so nothing is lost; the negative‑`omegaSquared` region is simply outside the Lean statement, consistent with FORMALIZATION_POLICY.md:178–181. (ii) `off_roots_excluded_positive` uses a purely transverse amplitude; the general theorem `off_roots_iff` covers all amplitudes, and a mixed amplitude such as `e1 + e3` would exercise both projections if a future revision wants a sharper concrete witness.

---

## 7. Work explicitly excluded from this clearance

Not reviewed and not certified here: the SO(D)/O(D) invariant census and the invariant‑count repairs; all other S11 packages (`XFORM_*`, `XCOEF_*`, `XKIN_ANISO`); S11b interface/passivity; S11c operators, calculations and files; bound‑state, radiation, leakage and coupling conclusions; systematic per‑output CAS expression bridging; the comparator/export pipeline, parser, registry and component‑census debts; variable coefficients, finite boundaries, nonlinear physics, confinement and general PDE/Fourier completeness; energy bases modulo total divergences. These are outside the bounded contract and I did **not** treat any of them as a prerequisite. The whole of `S10_exports.py` beyond the three named scalar records, and the audit engines beyond the selected constructor ranges, were likewise out of scope. Per FORMALIZATION_POLICY.md:107–115, H4 also requires a **second** independent non‑author review; this is one leg.

## 8. Verdict

**CLEAR** for the bounded H1–H4 homogeneous S11 contract at revision `S11 homogeneous bounded fidelity contract v1` (`4541cf82…`).

The supplied density, the coordinate and parameter maps, the zero‑inertia split, the integrated‑stationarity chain and the action‑derived operator are stated correctly and match the shared‑physics directive and both native constructors. The kernel classification is exhaustive, states full subspaces with genuine `finrank` dimensions rather than witness vectors, contains the coincidence locus rather than charting around it, and correctly separates `B = 0` and `k = 0`. `M_A = −E` and `M_B = E/2` are correct, and their signs and factors are retained rather than discarded. The threshold theorem is kinematic and is disclaimed as such. The mutation controls are meaningful, their passing counterparts are nonvacuous and admissible, and both source mutants fail on identities I independently confirmed to be false. The editorial correction about longitudinal gradients is mathematically right at the stated ansatz and is inert for the formal claims.

The seven recommendations above are documentation, anchor‑coverage and one‑line‑lemma improvements; none of them blocks fidelity clearance at this revision. This report is one non‑author review leg and does not by itself close H4.
