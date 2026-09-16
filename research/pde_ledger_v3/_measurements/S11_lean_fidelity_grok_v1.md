# Independent statement-fidelity review — S11 homogeneous bounded Lean contract

**Verdict: CLEAR** for the bounded H1–H4 contract at this revision.

**Packet:** `S11 homogeneous bounded fidelity contract v1`  
**Aggregate SHA256:** `4541cf8243b2058f813a09cf316b165a6d730ea53153f2284e27b39e5d4e51cf`  
**Reviewer:** Grok 4.6 (xAI), independent non-author. No other reviewer was consulted.  
**Computation:** none. No builds, no CAS, no hash recomputation, no script execution. Source inspection and algebraic reasoning only.

This is a statement-fidelity review of a frozen packet. Compilation, `PASS` labels, and hashes were not treated as proof of fidelity or of a meaningful mutation. Filed leftover goals were read as mathematics.

---

## Sources inspected and limits

**Read in full (or through the load-bearing end of the file):**  
`MANIFEST.json`; `lean/FORMALIZATION_POLICY.md`; `lean/AGENTS.md`; `lean/s11/COVERAGE.md`; `lean/s11/FIDELITY.md`; `lean/s11/README.md`; `lean/s11/VERIFICATION.txt`; `lean/s11/S11Homogeneous.lean`; `lean/s11/S11Homogeneous/Action.lean`; `lean/s11/S11Homogeneous/Spectrum.lean`; `lean/s11/S11Homogeneous/Threshold.lean`; `_measurements/S11_lean_source_check.py`; `_measurements/S11_lean_source_checks.json`; `_measurements/S11_lean_contract_check.py`; `_measurements/S11_lean_fidelity_review_prompt.txt`.

**S10/S9 actually used by the new modules, plus D3 curl support:**  
`S10Pilot/Action.lean`, `PlaneWave.lean`, `Spectrum.lean`, `PhaseAverage.lean`, `Variation.lean` (through the integrated stationarity theorem), `Analytic.lean` (SmoothField/TestField and compact calculus); `S10Controls/Action.lean`, `PlaneWave.lean`, `PhaseAverage.lean`, `Variation.lean`, `Analytic.lean`; `S10Pilot/Specialization.lean`; `S9Pilot/Action.lean` (`jetCurl` / density); `S9Pilot/Spectrum.lean` and `PlaneWave.lean` as D3 curl/phase support. Not reopened as S9 work.

**Native MAIN constructors and Q11, at the cited locations:**  
`scripts/S11_stray_longitudinal_sympy_audit.py` (selected constructor ranges from `S11_lean_source_checks.json`, plus `q11_objects` at 1672–1698); `mathematica/S11_stray_longitudinal_mathematica_audit.wl` (`curlDensity`/`divDensity` 835–841, `stiffnessBlueprint`/`kineticRecords` 1023–1038, `runCell` Lagrangian/EL/average/matrixB 2117–2241, `emitQ11` 2058–2114); `directives/S11_SHARED_PHYSICS.md` §§1–3, §7 MAIN pair, Q11 from line 880; `steps/S11_stray_longitudinal.md` Moves 2–3 and 6 (pure-trace correction and kinematic threshold).

**S10 export, only as specified:**  
`_RELATIONALS` / `_restore` at `scripts/S10_exports.py` 10–21, and the three scalar records at lines 27, 35, 1146 (`rho_br`, `mu_R`, `omegaSquared`). The rest of that export was not a review obligation.

**Filed verification evidence:**  
`_measurements/S11_lean_contract_checks.json` — axiom audit text, both source-mutant leftover goals, generated control sources, and `⊢ False` diagnostics for the concrete mutants that were opened. The JSON is large because it embeds full mutated modules; not every later control body was re-read once it matched the generator in `S11_lean_contract_check.py`.

**Not inspected, and not required for this contract:** remaining S10/S9 certificate/finite-action modules; the rest of `S10_exports.py`; other S11 packages; invariant-census bodies beyond confirming MAIN does not use `P_D`; S11c; interface/passivity; `docs/model_map.md` as a physics source (banner only).

Hashes in the manifest, `VERIFICATION.txt`, and the two measurement JSON files agree with each other for the named live paths. Hashes were not recomputed; they are provenance, not a certified Lean↔CAS translator.

---

## Contract identification (H1–H3)

### Action, coordinates, parameters

Lean density (`S11Homogeneous/Action.lean` 15–17):

\[
L=\frac{\rho}{2}|u_t|^2-\frac{\mu}{2}S_{\mathrm{curl}}(J)-\frac{B}{2}(\operatorname{div} J)^2,
\]

with S10 `stiffness` \(=\frac12\sum_{ij}(J_{i+1,j}-J_{j+1,i})^2\) and `divergence` \(=\sum_i J_{i+1,i}\). That is exactly MAIN in `S11_SHARED_PHYSICS.md` 954–967 and the native builders:

- SymPy `package_build` (`S11_stray_longitudinal_sympy_audit.py` 429–475): `Term(mu_R/2,"S_curl",S_curl)+Term(B_comp/2,"S_div",S_div)` subtracted from `rho_br/2\sum v^2`, with `S_curl=\frac12\sum(G_{ij}-G_{ji})^2` and `S_div=(\operatorname{tr} G)^2` (331–340).
- Wolfram `stiffnessBlueprint` / `kineticRecords` (1023–1038) and `lagrangianJet = Total[kineticTermsJet]-Total[stiffnessTermsJet]` (2177).

Parameter map is as claimed: `rho↔rho_br/rhoBr`, `mu↔mu_R/muR`, `B↔B_comp/bComp`, `cs↔c_s0/cs0`. Upstream scalars restored from `S10_exports.py` 27, 35, 1146 are `Symbol('rho_br', positive=True)`, `Symbol('mu_R', positive=True)`, `Symbol('omegaSquared', real=True)`. `B_comp` is declared locally, positive.

Coordinates match by named derivatives, not argument order:

| Object | Lean | Native |
|---|---|---|
| Point | time index 0, spatial `i.succ` | functions of `(x_1,…,x_D,t)` |
| Jet / `G` | `J j i = ∂_j u_i` | `G[i,j]=∂_{x_i} u_j` |
| Phase | `k·x-ωt` via `waveCovector = Fin.cases (-omega) k` (`S10Pilot/PlaneWave.lean` 22–24) | `Σ k_m x_m − ω t` (shared physics 81; Wolfram 2143–2144) |
| Plane wave | `a_i cos(phase)` | same real cosine ansatz |

`fieldJet(planeWave)=(-sin φ)•modeJet(-ω)k a` (`S10Pilot/PlaneWave.lean` 54–59) gives \(∂_t u= a\,ω\sinφ\) and \(∂_{x_i}u=-a_j k_i\sinφ\), which is the SymPy `route_b` substitution (411–416) and the Wolfram `planeGradientRules` / `planeVelocityRules` (2217–2222).

Kinetic energy is counted once: `lagrangian_split` (19–23) adds S10 curl at `(rho,mu)` to the divergence-only control at inertia `0` and stiffness `B`.

### Operator and native factors \(M_A=-E\), \(M_B=E/2\)

Lean operator (`Action.lean` 33–34):

\[
Ea=(\rho ω^2-μ|k|^2)a+(μ-B)(k·a)k,
\]

i.e. \(E=(ρz-μK)I+(μ-B)kk^T\). This is the sum of the two already derived operators (`modalOperator_split` 36–41). On a plane wave, the local EL (variational sign \(-\sum_j∂_j(∂L/∂J_j)\)) is \(\cosφ\cdot Ea\) (`eulerLagrange_planeWave` 169–173). Integrated compact-test stationarity is equivalent to pointwise EL (`actionStationary_iff_eulerLagrange` 145–167) and, for plane waves, to `ModalStationary` (`actionStationary_planeWave_iff` 175–182). Backgrounds need not decay; only tests are compact (`S10Pilot/Variation.lean` 3–6, 83–86, 162–163).

Native route A uses the opposite EL sign: SymPy `amp_eq` with \((-ω^2 a)\) and \((-k_i k_ℓ a)\) (381–395); Wolfram `eulerExpressions` is \(∂_t(∂L/∂v)+∂_x(∂L/∂G)-∂L/∂u\) (2183–2191). For this quadratic \(L\), that is \(-E\) after the cosine is stripped (`route_A_stripped_phase` recorded as `cos(phase)`). So \(M_A=-E\).

Native route B: Lean `phaseAverage` is the actual \((1/2π)\int_0^{2π}\) (`Action.lean` 63–73) and equals \(\frac12\) times `modalAction`. `modalAction=\frac12 a·Ea`, so the average is \(\frac14 a·Ea\). The Hessian of that quadratic is \(E/2\). Wolfram literally integrates then twice-differentiates (2225–2231) and selects `M_B` (2239–2241). SymPy `period_average` replaces \(\sin^2φ\mapsto 1/2\) (398–404); for a quadratic density that is the same average. Downstream matrix is \(M_B\).

Kernels of \(E\), \(-E\), and \(E/2\) coincide. Resolvent/residue normalization is correctly not claimed (`FIDELITY.md` 74–76).

Handwritten checker formulas (`S11_lean_source_check.py` 74–75, 86–92) are the Lean density with D3 `|curl|^2` and the Lean matrix \(E\). Algebraically they match the Lean definitions and the MAIN constructors. Filed residuals: five identities zero; reversed compression and omitted half-factor nonzero. Those residuals were not recomputed here.

D3 curl normalization, as supporting evidence only: `S10Pilot.stiffness_three` (`Specialization.lean` 16–18) equals `|S9Pilot.jetCurl|^2`, and `jetCurl` (`S9Pilot/Action.lean` 22) is the ordinary curl. The checker's `|curl|^2` handwritten density is therefore the same D3 object as S10 `stiffness`.

### Exhaustive kernels (H2)

On `rho≠0`, `k≠0`, `modeSpace_classification` (`Spectrum.lean` 122–140) is exhaustive, including coincidence:

| Hypothesis | Kernel | D3 `finrank` (`kernel_census_three` 142–162) |
|---|---|---|
| \(ω^2=μK/ρ\), \(B≠μ\) | `transverseSpace k` | 2 |
| \(ω^2=BK/ρ\), \(B≠μ\) | `longitudinalSpace k` | 1 |
| \(B=μ\) and on that common root | `⊤` | 3 |
| off both roots | `{0}` | 0 |

These are geometric dimensions of full subspaces, not sampled null vectors. `frequency_coincidence_iff` (36–47) is exactly \(B=μ\) on this domain. `transverse_operator` (25–27) holds for every \(B\); `longitudinal_operator` (29–34) is independent of \(μ\). `off_roots_iff` (78–94) is the whole complement, not a witness. `positive_frequencies` (164–172) supplies positive square roots under `rho,mu,B>0`, `k≠0`; both frequency signs sit in the \(ω^2\) classification. The zero vector is in every kernel and does not add a mode.

`B=0` recovers S10 on the nose (`zero_compression_action` / `_operator` / `_stationarity`, `Action.lean` 184–190 and `Spectrum.lean` 174–181). S10’s static longitudinal line is then a direct implication of `S10Pilot.zero_frequency_iff` (`S10Pilot/Spectrum.lean` 103–122) plus operator identity, with `mu≠0` inherited from S10. `k=0` is separate: `zero_wavevector_iff` (183–186) is \(ω=0∨a=0\), independent of \(μ,B\), and does not use the nonzero-\(k\) dimension formula.

Physical domain in the coverage contract is D=3, positive `rho,mu,B`, nonzero real \(k\). The classification lemmas are slightly more general (`rho≠0` rather than `>0`); that does not overclaim the D=3 completion.

### Kinematic threshold (H3)

Q11 supplies only \(ω^2=c_{s0}^2(K+k_w^2)\) and shared \(ω,k\), with no interface law (`S11_SHARED_PHYSICS.md` 880–894; SymPy 1678–1681; Wolfram 2064–2091). Lean `phase_matching` (`Threshold.lean` 14–19) substitutes the proved longitudinal root `coneValue rho B k` and gets

\[
k_w^2=K\Bigl(\frac{B}{ρ c_s^2}-1\Bigr)=ω^2/c_s^2-K.
\]

`threshold_classification` (21–37), under `rho>0`, `cs>0`, `k≠0`, is the three-way sign equivalence with \(B\lessgtr ρ c_s^2\). At \(k=0\) the square vanishes even off that locus, so the nonzero-\(k\) guard is necessary. The module header (3–4) and coverage text correctly refuse bound-state, leakage, and coupling claims.

### “Pure trace” correction at the homogeneous ansatz

`steps/S11_stray_longitudinal.md` 43–46 is right for this ansatz. Longitudinal \(a=αk\) gives \(G=-α\sinφ\,kk^T\), which is symmetric and curl-free. Its symmetric-traceless part is \(-α\sinφ\bigl(kk^T-(|k|^2/D)I\bigr)\), nonzero for \(D=3\), \(k≠0\). MAIN charges curl and trace only; curl-free is the property the selected action uses. The correction does not change MAIN and does not certify an invariant basis.

---

## Controls (H4)

Axiom audit (`S11Homogeneous.lean` 5–33; filed output in `S11_lean_contract_checks.json` 130): 29 audited declarations, only `propext`, `Classical.choice`, `Quot.sound`. That is a kernel-axiom report, not fidelity.

Source mutants fail on the intended identities, not on import/syntax/resource:

1. Compression sign (`Action.lean` 16–17 flipped to `+ B/2`; `lagrangian_split` retained). Leftover (`contract_checks.json` 154–160):  
   `B * divergenceOnlyStiffness J * (1/2) = B * divergenceOnlyStiffness J * (-1/2)`.  
   False at \(B=1\), \(\operatorname{div} J=1\).

2. Phase half (`Action.lean` 68, `1/2` replaced by `1`; integral retained). Leftover (190–196):  
   `modalAction = 2 * modalAction`.  
   False whenever the modal action is nonzero; the documented sample \(ρ=μ=1,B=4,ω=2,k=(0,0,1),a=(1,0,0)\) has modal action \(3/2\).

Concrete controls are generated from `S11_lean_contract_check.py` 46–62 and use `kernel_census_three` for 2/1/3 versus false 3/2/2 at an admissible \(k=(0,0,1)\), not witness-vector counts. Mode and threshold mutants that were opened leave `⊢ False` in `contract_control`. Passing counterparts are non-vacuous (dimensions 2, 1, 3; both stationary and non-stationary claims; all three threshold signs). `B=0` and `k=0` are separate samples. Transverse stationarity at \(B=9\) on the shear cone, and a longitudinal amplitude failing at \(ω=1\) when \(B=4\), match the stated independence claims.

The compact native checker’s unused `PD` sentinel (`S11_lean_source_check.py` 68–70) is absent from MAIN, so the generic census input is not silently feeding this action.

---

## Findings

### Mathematical / fidelity blockers

None for this bounded contract.

### Necessary control gaps

None. Every load-bearing claim in the `FIDELITY.md` mutation map has an identified control, and the source-mutant leftovers are the intended identities.

### Nonblocking recommendations

1. **`grazing_requires_nonzero_wavevector` mutant is a dummy.**  
   Generator: `S11_lean_contract_check.py` 59–62. The passing statement is  
   `normalWaveSq 1 1 2 0 = 0 ∧ (1:ℝ) ≠ 1*2^2`.  
   The second conjunct is unrelated; the mutant `1=4` fails without touching the \(k=0\) identity. The passing conjunct already exhibits \(k_w^2=0\) off the coefficient locus. Smallest fix: mutate to `normalWaveSq 1 1 2 0 ≠ 0`.

2. **`B=0` recovered line is identified by operator equality, and H4 only positively samples the longitudinal vector.**  
   `zero_compression_operator` (`Action.lean` 188–190) plus `S10Pilot.zero_frequency_iff` already give the full line under `mu≠0`. A static-transverse `B=0` exclusion control would make that visible in H4 without reading S10. Not required, given the identity theorems.

3. **Native `period_average` is a \(\sin^2\mapsto 1/2\) rewrite, not the integral Wolfram and Lean use.**  
   `S11_stray_longitudinal_sympy_audit.py` 398–404 versus Wolfram 2225–2227 and `Action.lean` 63–73. For quadratic MAIN they agree; the checker residual is the right compact test. Optional: one sentence in `FIDELITY.md` that the SymPy route is the quadratic rewrite, while Lean/Wolfram integrate.

### Excluded work (not blocking, not requested)

SO(D)/O(D) invariant completeness and `P_D`; other S11 packages (`XFORM_*`, `XCOEF_*`, `XKIN_ANISO`); S11b interface/passivity; S11c; bound-state / radiation / leakage; variable coefficients and boundaries; Lean polynomial-root multiplicity; systematic CAS bridging; whole-export execution; uniqueness of an invariant basis. Q11 closure tags `C1`–`C4` measure missing interface content and are out of scope. Historical engine debts are not waived and are not this contract.

---

## Disposition

The supplied homogeneous density, the combined variational operator, the native factors \(M_A=-E\) and \(M_B=E/2\), the exhaustive D3 full-amplitude kernels including coincidence and off-root exclusion, positive-frequency existence, \(B=0\) and \(k=0\) boundaries, and the three-way kinematic grazing threshold are the objects claimed in H1–H3, with the stated quantifiers and with mutation diagnostics that fail for mathematical reasons.

**CLEAR** for H1–H4 at packet revision v1, aggregate `4541cf8243b2058f813a09cf316b165a6d730ea53153f2284e27b39e5d4e51cf`. Independent review of this fixed packet is complete on this leg. Closing H4 still requires the second independent non-author review and disposition of any findings it raises.
