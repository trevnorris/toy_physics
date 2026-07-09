# II-G3a (ledger_stage016) source map — pathA_32 ℓ=2 SO(3) covariance theorem (ISOTROPY_CALIBRATED 1/2, the EARNED slice)

> Running-start prep captured 2026-07-08 (post stage015, in the /clear prep) so a fresh session can author the reshape
> directive without re-discovery. **All line refs below VERIFIED against the current sources 2026-07-08**
> (`pathA_32_grouped_p2_isotropy_sympy.py` = 1496 lines; `pathA_32_grouped_p2_isotropy.wl` = 696 lines;
> `reports/pathA_32_grouped_p2_isotropy.md`).
> Companion: `part2_gravity_atomic_split.md` (rows 016/017 + the pathA_32 trip-ups L84–86). **pathA_32 SPLITS 2-way:
> 016 (THIS: the SO(3) covariance theorem = the EARNED angular structure) + 017 (II-G3b grouped-P2 lane isotropy =
> the grouped {20,21,22} lanes, raw-D=0 PRIMARY, normalized-u cross-check, calibration partition; exports the ℓ=2 port
> kernel).** Build-order id **016**, Part II. Source top-line: **`ISOTROPY_CALIBRATED`** (the joint; 016 lands the 1/2
> EARNED component, 017 completes + lands the calibration).
> ⚠ **No prior source map existed for stage016** — this is the FIRST authored artifact of the stage (per the calibrated
> pipeline: source map → reshape directive → reviews → gate).

## ⭐ The THREE headline differences from 015 (READ FIRST)
1. **⚠ 016 is a 2-way split (016+017), NOT the 3-way pathA_31 split, and 016 is the EARNED-FIRST slice** (like 011 was
   the EARNED reduction and 012 the completing DtN ladder). 016 = the **angular** covariance theorem: real ℓ=2 harmonics,
   angular Gram = I₅, the COMPUTED `−Δ_S²` eigenvalues `λ_m=6` (genuine Laplacian + Rayleigh + eigenfunction residual,
   NOT typed), and the K₂ **angular** stiffness `K₂ = K̃ + λ_m·T̃_Ω` with the eigenvalue that ENTERS K₂ being the computed
   one. 017 = the **radial/response** lane gate (assemble lanes → raw-D=0 → normalized-u cross-check → calibration
   partition). The cut is angular-structure (016, earned) vs grouped-lane-response + calibration (017).
2. **⭐ The verdict-vs-PASS split: 016 EARNS the covariance theorem (adds NO calibration); the CALIBRATED label is 017's.**
   The source's single verdict `verdict_from_gates` (`.py` L280–299) lands `ISOTROPY_CALIBRATED` **because
   `calibration_inputs` is non-empty** (L297–299) — the radial/support scalars are FROZEN, not derived from the Gate-1
   `R0` support equation (report :5). ⚠ **This is the ISOTROPY_CALIBRATED essence: the ANGULAR structure is EARNED
   (016), the RADIAL profile is CALIBRATED/frozen (017).** So 016 prints `ISOTROPY_CALIBRATED (1/2) — SO(3) covariance
   theorem EARNED`, CITES the frozen packet it consumes, and DEFERS the calibration partition + the joint landing to 017
   (the 013-prints-partial-component pattern). Do NOT let 016 re-present 017's calibration accounting.
3. **⭐⭐ The pathA_32 trip-ups live HERE (016-owned): the aggregate probe battery + computed-not-typed eigenvalues +
   the T_Ω FAIL_DIMENSIONAL.** From the split's rig-history (`part2_gravity_atomic_split.md` L84–86): (a) **keep the
   aggregate probe battery intact** — neuter any one probe's computed flag → the aggregate `able_to_fail_ok` flips false
   (`.py` `able_to_fail_if_probe_neutered` L1313–1316, `neutering_any_probe_flips_false` L1445); (b) **the eigenvalues
   must be COMPUTED, not typed** — the `tautology_hash_collision` probe (identical harmonic inputs → non-distinct hashes →
   `FAIL_TAUTOLOGICAL`, `.py` L1133–1138/L1219–1227) + the eigenfunction residual `residuals_zero` (L817/L821); (c) **a
   dim-mutation on the sourced `T_Ω` must fire `FAIL_DIMENSIONAL`** (`.py` `build_dimensional_check` corrupt-probe
   L361–382, verdict L1245–1257). ⚠ These are 016's central able-to-fail obligations (per-tooth ablation WILL probe each).

## §1 The 016 slice (`.py` line ranges) — the CLEAN CUTS (all VERIFIED)
The whole computation lives in `symbolic_engine()` (`.py` L762–1458). 016 owns the **angular covariance** sub-block +
its dedicated probes; 017 owns the **grouped-lane response + calibration partition**. The clean cuts:

- **Real ℓ=2 harmonics construction — `harmonics` L777–783** (five real harmonics `20/21c/21s/22c/22s`), `order`/`ys`
  L784–785. The physical basis of the covariance theorem.
- **Angular helpers — `integrate_s2` L121–127** (the S² measure `sin θ dθ dφ`) + **`laplacian_s2` L130–134** (the real
  `−Δ_S²` operator). These are 016's engine primitives (017 reuses `integrate_s2` for self-overlaps — shared helper, not
  a cut boundary).
- **Angular Gram = I₅ — L787–788** (`gram = [[integrate_s2(ys[i]·ys[j])]]`, `gram_is_identity = bool(gram==eye(5))`). The
  orthonormality theorem.
- **⭐ The COMPUTED `−Δ_S²` eigenvalues λ_m=6 — L808–822** (the CRUX EARNED result). For each harmonic:
  `neg_lap = −laplacian_s2(y)` L813; `rayleigh = integrate_s2(y·neg_lap)/integrate_s2(y²)` L814 → `lambdas[name]` L815;
  `k_coeff_used[name] = rayleigh` L816 (⭐ **the coefficient that enters K₂ IS the computed eigenvalue**); the
  **eigenfunction residual** `residuals[name] = neg_lap − rayleigh·y` L817 (proves `y` is a genuine eigenfunction, not
  just a Rayleigh number); `k_coeff_equal[name]` L818. Aggregates: `lambda_all_six = all(λ−6==0)` L820,
  `residuals_zero = all(residual==0)` L821, `k_coefficients_equal_all` L822.
- **The K₂ ANGULAR stiffness form — `response_symbols` L1381–1387** (build-time payload): `K2_per_channel = "K̃ +
  λ_m·T̃_Ω"` L1383, `K2_on_isotropic_l2 = "K̃ + 6·T̃_Ω"` L1384. ⚠ **016 owns the bare angular stiffness `K̃ + λ_m·T̃_Ω`
  with λ_m the computed eigenvalue.** (017's `assemble_channel` L843 multiplies this by the per-lane angular self-overlap
  `pref`; for the "20" channel `pref = integrate_s2(Y20²) = 1` — a 016 Gram-diagonal fact — so `grouped_lanes["20"]["K2"]`
  L843 numerically equals the bare `K̃ + 6·T̃_Ω`. This is why the dimensional check can validate the core angular
  stiffness through the "20" lane; see below.)
- **The angular DIMENSIONAL check + the T_Ω FAIL_DIMENSIONAL probe — `build_dimensional_check` def L302–487, called
  L885–892.** Reconstructs the density-level integrand (`μ_η`, `T_w`, `K_η`, `T_Ω` densities × the S² measure
  `a²·dw·dΩ`, L316–321), checks `[M₂]=M`, `[K₂]=M T⁻²`, `[K₂/M₂]=T⁻²` and term-homogeneity (L350–359), and the
  **corrupt-`[T_Ω]` probe** L361–382 (`corrupt_dims[T_Ω] += (1,0,0)` and `corrupt_dims[T̃_Ω] += (1,0,0)` L362–363 →
  the `K₂` sum becomes inhomogeneous → `DimError` → `corrupt_ok=False` → `probe_verdict = FAIL_DIMENSIONAL` L382). ⚠
  **This is trip-up (c) — 016-owned.** The call passes `M2_expr=grouped_lanes["20"]["M2"]`, `K2_expr=grouped_lanes["20"]
  ["K2"]`, `lambda_m=6` L889–891 — see §1c for the clean-cut handling (the "20" lane = the bare angular stiffness).
- **`tautology_clear` (harmonics computed-not-typed) — L940** (`distinct_hashes and all(self_overlap==1)`); the input
  hashes + per-lane trace L929–939. ⚠ trip-up (b): the harmonics are genuinely distinct (`distinct_hashes` L930) and
  each self-overlap = 1 (the Gram diagonal = orthonormality). 016-owned.
- **The 016-owned able-to-fail probes (the covariance teeth):**
  - **`forced_eigenvalue_probe` def L1100–1114** → instantiated `wrong_eigen_probe = forced_eigenvalue_probe(2)` L1129
    + ablation `forced_eigenvalue_probe(6)` L1130. Forces the K₂ angular coefficient to 2 while the computed `λ_m`
    stays 6 → `coefficient_equals=False` → `covariant=False` → `FAIL_NOT_COVARIANT` (L1106/L1113). The `wrong_eigenvalue`
    probe payload L1192–1200. ⚠ trip-up (a)-member + the covariance headline tooth.
  - **`lane_hash_probe` def L1116–1125** → `tautology_probe = lane_hash_probe({all three lanes = harmonics["20"]})`
    L1133–1135 + ablation `{20,21c,22c}` L1136–1138. Identical inputs → non-distinct hashes → `FAIL_TAUTOLOGICAL`
    (L1119/L1124). The `tautology_hash_collision` probe payload L1219–1227. ⚠ trip-up (b).
  - **The `dimensional_corrupt_T_Omega` probe** L1245–1257 (wraps the `build_dimensional_check` corrupt-probe;
    `dimensional_verdict = case_verdict(dimensional_ok=with_mutation)` L1145 → `FAIL_DIMENSIONAL`; self-ablation restores
    L1250–1256). ⚠ trip-up (c).
- **The `covariant` gate — `baseline_gates["covariant"]` L1322** = `gram_is_identity and lambda_all_six and
  residuals_zero and k_coefficients_equal_all`. ⭐ **This one boolean IS the 016 EARNED headline** (Gram=I₅ ∧ λ=6 ∧ genuine
  eigenfunctions ∧ the K₂ coefficient = the computed λ). It gates `FAIL_NOT_COVARIANT` (L283–284).

- **⭐ CLEAN CUT — 016 owns the harmonics + Gram=I₅ + the computed eigenvalues + the K₂ angular-stiffness form + the
  angular dimensional check + the three covariance teeth (wrong-eigenvalue / tautology-hash / T_Ω-dim) + its own aggregate
  battery over THOSE probes. Do NOT pull in 017 territory:**
  - 017: `assemble_channel` L837–859, `ungrouped` L861, `grouped_lanes` L866–884, `cs_equal` L894–903, the raw-D gate
    L905–914, the normalized-u cross-check L916–927, the calibration/stability windows L942–953, the
    `derived_inputs`/`calibration_inputs` PARTITION L955–978, `d_common`/`c_group` L994–1003, and the response probes
    `pure_prefactor`/`sector_selective`/`m_dependent`/`degenerate_beta`/`singular_denominator`/`static_drop_inertia`
    (L1006–1098, L1127–1128, L1131–1132, L1139–1143, payloads L1148–1244). **017 CITES 016's harmonics + eigenvalues +
    K₂ form; do not re-derive them in 017 either.**
  - Shared machinery (both stages reference; 017 lands the full form): `verdict_from_gates` L280–299, `case_verdict`
    L980–992, the aggregate `expected_probe_verdicts`/`computed_probe_gate_flags`/`able_to_fail` L1260–1316, `baseline_gates`
    L1320–1329, the final `verdict` L1330. ⚠ **016 keeps its OWN aggregate over its 3 probes (neuter-one flips 016's
    aggregate); 017 keeps its OWN over its 6 probes.** The 9-probe union is reconstituted conceptually across 016∧017
    (the BREATHING_CALIBRATED=013∧014∧015 pattern) — do NOT force 016 to carry the full 9-probe battery.

## §1c The dimensional-check clean cut (READ — the one genuinely-shared object)
`build_dimensional_check` is called with 017's `grouped_lanes["20"]` M2/K2 (L889–891). But for the "20" channel the
angular self-overlap `pref = integrate_s2(Y20²) = 1` (a **016 Gram-diagonal fact**), so `grouped_lanes["20"]["M2"] = M̃`
and `["K2"] = K̃ + 6·T̃_Ω` — i.e. the **bare core angular stiffness**, no radial-lane content. ⭐ **In the reshape,
016 should build its dimensional check on the CORE angular stiffness `K₂ = K̃ + λ_m·T̃_Ω` (and `M₂ = M̃`) directly**
(citing `pref₂₀=1` from its own Gram diagonal), so 016 does not import 017's `grouped_lanes` object. The `.wl` already
does exactly this pattern (`actualK2Dim = dimOf[groupedLanes["20"]["K2"]]` L205, where `groupedLanes["20"]` is the "20"
channel with `cSelf=1`). This keeps the T_Ω FAIL_DIMENSIONAL tooth wholly inside 016.

## §1b The `.wl` 016 slice (VERIFIED) — the independent Wolfram route
⭐ **The pathA_32 `.wl` is ALREADY a genuinely independent engine** (native `Integrate`/`FullSimplify` S² integrals,
native `D[...]` Laplacian, native `Hash[...,"SHA256"]`, its OWN `verdictFromGates`/`caseVerdict`, its OWN dimensional
walker `dimOf`). It does **NOT** `Get`/`Import` the `.py`'s expressions. So the reshape KEEPS the native route and only
severs the scratch-YAML handoff (see §3). 016's `.wl` slice:
- harmonics L27–35; `intS2` L37–43 + `lapS2` L44–48 (native angular primitives).
- `gram` L93 + `gramIsIdentity` L94 (Gram=I₅).
- `negLaps` L105 + `lambdas` (Rayleigh) L106–112 + `residuals` L113–116 + `lambdaAllSix` L117 + `residualsZero` L118 +
  `kCoeffUsed`/`kCoeffEquals`/`kCoeffEqualsAll` L119–124 (the computed eigenvalues + residual + K₂-coefficient identity).
- The angular dimensional check L178–225 (`measureExpr`/`m2IntegralExpr`/`k2IntegralExpr` L178–183, `dimRules` L184–197,
  the native `dimOf` L57–81, `corruptDimRules` L213–216 → `corruptKOmegaDim`/`corruptDimensionalOk` L217–224,
  `dimensionProbeVerdict` L225) + **the standalone dim assert L227–229** (already print-only/`fail[]`-on-failure — a
  reshape-ready idiom).
- `covariantOk` L500 (= `gramIsIdentity && lambdaAllSix && residualsZero && kCoeffEqualsAll`).
- `tautologyClear` L494–497 (native `ToString[...,InputForm]` distinctness + `intS2[y²]==1`).
- The 016 probes: `forcedEigenvalueProbe` L404–414 (+ `wrongEigenProbe`/`Ablation` L426–427), `laneHashProbe` L416–422
  (+ `hashText` L415, `tautologyProbe`/`Ablation` L446–447), the dim probe verdicts `dimensionalProbeVerdict`/
  `dimensionalAblationVerdict` L453–454.
- **⚠ The bridge to sever (§3):** the `.wl` writes `yamlOut` (the scratch YAML `pathA_32_mathematica_results.yaml`) via
  `Export` L693 for the `.py` to read back → engine-agreement. SEVER: no scratch YAML, print-only + `fail[]`; the dual-
  engine agreement becomes **transcript-level** (both engines print the same computed λ=6 / Gram=I₅ / probe verdicts;
  the stage014 numeric-transcript pattern). KEEP the native route + the standalone dim assert + the arity discipline.
- ⚠ **17 territory in the `.wl` (EXCLUDE from 016):** `dCommon` L126–130, `assembleChannel` L132–155, `groupedLanes`
  L158–176, `csEqual` L231–240, raw/normalized defects L241–260, `calibrationInputs` L262–276, the response probes
  L321–390/L424–445/L448–452, `stabilityOk`/`denominatorGuardOk`/`laneCollapseOk`/`dynamicRetained` L492–499, and the
  full YAML writer L510–693.

## §2 The 016 claim-set (derive + assert; report quotes)
- **The SO(3) covariance theorem (EARNED — the headline).** The five real ℓ=2 harmonics are orthonormal on S²
  (`Gram = I₅`, report :9 `Gram matrix equals I5: True`) and are genuine `−Δ_S²` eigenfunctions with the SAME eigenvalue
  `λ_m = ℓ(ℓ+1) = 6` for every `m` (report :10 `Computed -Delta_S2 eigenvalues: {20:6, 21c:6, 21s:6, 22c:6, 22s:6}`),
  proven by BOTH the Rayleigh quotient AND the eigenfunction residual `−Δ_S²Y − λY = 0` (`.py` `residuals_zero` L821).
  ⭐ **This is the earned angular content: the ℓ=2 sector is a single 5-dimensional SO(3)-irrep, so any isotropic wall
  response is angular-degenerate (one `λ_m` across all five m) — the covariance that the grouped-lane isotropy (017)
  then rides.**
- **The K₂ angular stiffness uses the COMPUTED eigenvalue (EARNED).** `K₂ = K̃ + λ_m·T̃_Ω` with `λ_m` selected FROM the
  computed Laplacian spectrum (report :11 `K2 angular coefficient equals computed lambda_m: True`; `.py` `k_coeff_used =
  rayleigh` L816, `k_coefficients_equal_all` L822). ⚠ **NOT a typed 6** — the `forced_eigenvalue_probe(2)` tooth
  (report :22 `Wrong eigenvalue coefficient: FAIL_NOT_COVARIANT`) proves the gate rejects a typed coefficient ≠ the
  computed λ.
- **The angular dimensional consistency (EARNED able-to-fail).** `[K₂] = [K̃] = [λ_m·T̃_Ω] = M T⁻²`, `[M₂] = M`,
  `[K₂/M₂] = T⁻²` (report :53–55, the retrofit table); the sourced-`T_Ω` corrupt probe fires `FAIL_DIMENSIONAL`
  (report :26 + :56 the full probe-flip record). ⚠ trip-up (c).
- **Computed-not-typed genuineness (EARNED able-to-fail).** The harmonics are distinct (`distinct_hashes`) with
  self-overlap 1; the `tautology_hash_collision` probe (identical inputs) fires `FAIL_TAUTOLOGICAL` (report :24). ⚠
  trip-up (b).
- **The 016 aggregate battery over ITS three probes (EARNED able-to-fail).** `wrong_eigenvalue`, `tautology_hash_collision`,
  `dimensional_corrupt_T_Omega` each fire their own verdict with a self-ablation that suppresses the fail (report :31/:33/
  :35), and neutering any one flips the 016 aggregate (the `neutering_any_probe_flips_false` discipline, report :40 for
  the full battery — 016 keeps this over its subset). ⚠ trip-up (a).
- **The 016-scoped landing (partial component).** Land at the 016 component: `ISOTROPY_CALIBRATED (1/2) = the ℓ=2 SO(3)
  covariance theorem EARNED (real harmonics + Gram=I₅ + computed λ_m=6 + K₂ angular stiffness), with the grouped-lane
  isotropy + calibration partition = 017.` Do NOT print the joint as complete (that is 017's landing) and do NOT
  re-present 017's calibration accounting (cite the frozen packet as CONSUMED — §4).

## §3 Reshape cost (the bridge to sever) — cross-script scratch-YAML family, KEEP the native `.wl`
pathA_32 is the **cross-script runtime-YAML** reshape family (not the sympy-expr-import family): the `.py`'s
`symbolic_engine()` is pure/self-contained, but `main()` L1461–1495 writes `SYM_YAML` (L1464), reads the `.wl`'s
`MMA_YAML` (L1466), computes `compute_engine_agreement` (L618–759, L1467), and writes `RESULTS_YAML`/`REPORT_MD`/
`FEED_NOTE` (L1483–1487); the `.wl` writes its scratch YAML via `Export` L693. **Reshape = sever ALL file I/O both
directions:** drop `main`'s YAML/report/feed writers + the `yaml_read`/`yaml_write`/`compute_engine_agreement`
scratch-bridge (`.py`); drop the `Export` + the YAML-line assembly L510–693 (`.wl`). Each engine becomes standalone:
print-only, `expect_zero`-style asserts (`.py` local ledger idioms; the `.wl` already has `fail[]`), `raise`/`Exit[1]`
on failure. **KEEP the `.wl`'s already-independent native route** (§1b) — it is NOT a mirror; re-target it to assert its
OWN computed λ=6/Gram=I₅/probe verdicts. **Dual-engine agreement = transcript-level** (stage014 pattern). **Zero file
I/O.** Arity discipline (standing — def/call scan + unevaluated-leakage transcript scan).

## §4 Consumed / exported
- **Consumes (cite, dual-site integrity, don't re-derive):** the **frozen throat packet + Gate-1 D/N provenance** as
  CALIBRATION inputs (`.py` `calibration_inputs` L964–978): `R0(w)` linearized isotropic reference, `β₂(w)` ℓ=2 radial
  profile, the densities `μ_η(w)`, `T_w(w)`, `K_η(w)`, `T_Ω(w)`, the radial scalars `M̃`/`K̃`/`T̃_Ω`, the support scalars
  `B̃`/`Z̃`, the physical calibration window, and the **Gate-1 D/N boundary provenance** (= stage012's R28
  `BC_DEPENDENT`, and stage011's frozen-reduction record R26). ⚠ **`μ_η`, `T_w` are ALREADY counted (stage013 CALIB
  `{μ_η,T_w,β}`); `K_η = T_w β²` is DERIVED (stage013 R29).** Cite them dual-site (stage013 + stage011/012), do NOT
  re-count or re-derive. ⚠ **`c_S` NOT consumed** (matter-sector deferred, as 013/014/015). ⚠ **The ANGULAR structure
  016 earns is NOT consumed — it is derived here** (the covariance theorem); only the RADIAL profile/scalars are frozen
  (that is the ISOTROPY_CALIBRATED essence, §head #2).
- **Exports:** the ℓ=2 SO(3) covariance theorem (Gram=I₅ + `λ_m=6` + the K₂ angular-stiffness form `K̃ + λ_m·T̃_Ω`) →
  017 (which rides it for the grouped-lane isotropy) → then 017 exports the grouped-P2 ℓ=2 **port kernel** (`M₂`, the
  angular `K₂`, the `B̃/Z̃` support scalars, the `D`-lanes) → 018–021 (Λ₂ fingerprint / DtN Hankel) + 022/023 (ℓ=2
  cross-ℓ map) + 024 (wall mode). ⚠ **The port-kernel export is 017's, not 016's** — 016 exports only the angular
  covariance theorem.

## §5 Teeth candidates (016-specific)
1. **⭐ The covariance headline tooth (`wrong_eigenvalue`, FAIL_NOT_COVARIANT):** the K₂ coefficient must be the COMPUTED
   `λ_m` — ablation: force it to 2 → `coefficient_equals=False` → `FAIL_NOT_COVARIANT`; self-ablation `forced=6` suppresses.
   Per-tooth: mutate the forced coefficient and confirm the fail fires AT the covariant gate. ⚠ Also probe the eigenfunction
   `residuals_zero` genuineness (a wrong Rayleigh number without a true eigenfunction would still be caught by the residual).
2. **⭐ Computed-not-typed eigenvalues (`tautology_hash_collision`, FAIL_TAUTOLOGICAL):** identical harmonic inputs →
   non-distinct hashes → fail; self-ablation `{20,21c,22c}` distinct → suppressed. Per-tooth: the eigenvalues are read
   from the Laplacian spectrum, never string-typed.
3. **⭐ The sourced-`T_Ω` dimensional tooth (`dimensional_corrupt_T_Omega`, FAIL_DIMENSIONAL):** corrupt `[T_Ω]` (and its
   assembled `T̃_Ω` scalar) by one power of L → the K₂ sum goes inhomogeneous → `FAIL_DIMENSIONAL`; self-ablation restores.
   Per-tooth: the mutation participates in the verdict (`participates_in_verdict=True`, `.py` L470).
4. **The Gram=I₅ tooth:** the orthonormality is a genuine `integrate_s2` computation — a non-orthonormal basis breaks
   `gram_is_identity` → `covariant=False`. (Confirm this is reachable — it is part of the `covariant` gate L1322.)
5. **The 016 aggregate-battery discipline:** neuter any ONE of {wrong_eigenvalue, tautology_hash, T_Ω-dim} → 016's
   aggregate `able_to_fail_ok` flips false. ⚠ trip-up (a) — keep the battery intact (do NOT let a probe become a no-op).
6. **Consumed-packet dual-site** (frozen densities/scalars cited from 013 + 011/012; a corruption fires both engines) +
   `.wl` arity scan.
⚠ **NOT 016 (do not pull in as 016 teeth):** the lane-collapse / anisotropy-branch probes (pure_prefactor / sector_selective
/ m_dependent), the stability (degenerate_beta) / denominator (singular) / static-response probes — those are 017's battery.

## §6 Register expectation — ⭐ THE KEY 016 QUESTION (likely ZERO new counted knobs; CONFIRM)
Per the split pin (3): **the SO(3) covariance theorem is EARNED (angular structure); the calibration partition is 017's.**
So the honest pre-read (confirm at the register step + Codex-verify against the scripts):
- **016 likely adds ZERO new counted knobs** (like 011/012/014 — an EARNED/structural stage). The covariance theorem is
  pure angular math (Gram=I₅, `λ_m=6` from the Laplacian) — it introduces NO calibration. It CITES the frozen packet
  (μ_η/T_w already counted at 013; K_η=T_wβ² DERIVED R29) as CONSUMED calibration inputs.
- **⚠ The `T_Ω` (angular stiffness density) + `β₂(w)` (ℓ=2 radial profile) counting question is 017's, NOT 016's.**
  `T_Ω` FIRST APPEARS here (016 uses it dimensionally in the K₂ angular stiffness), and `β₂(w)` is the frozen ℓ=2 radial
  profile. Whether `T_Ω` is a NEW counted CALIB knob (the angular-stiffness scale — 017 would be the calibration-adding
  stage) or is absorbed/derived MUST be determined at 017's calibration partition + Codex-verified. ⭐ **016 should DEFER
  this** — cite `T_Ω` as consumed-frozen, resolve its count in 017. Do NOT pre-decide; do NOT let 016 silently count
  `T_Ω` (that would pre-empt 017's partition).
- **Likely a structural edge** recording the ℓ=2 SO(3) covariance provenance (the five real harmonics form one SO(3)
  irrep with `λ_m=6`; the K₂ angular stiffness uses the computed eigenvalue) — discharges NOTHING (it is earned angular
  structure, not a reduction of a debt). CONFIRM at registration.

## Verdict tokens + honest scope
016 carries the **SO(3) covariance-theorem component (1/2) of `ISOTROPY_CALIBRATED`** — the EARNED angular structure:
real ℓ=2 harmonics, angular Gram = I₅, the computed `−Δ_S²` eigenvalues `λ_m=6` (Rayleigh + eigenfunction residual), and
the K₂ angular stiffness `K̃ + λ_m·T̃_Ω` with the computed eigenvalue. EARNED = the angular covariance is DERIVED (not
calibrated); the CALIBRATED label + the grouped-lane isotropy + the calibration partition = 017. Caveats: the RADIAL
profile/scalars are FROZEN calibration inputs (the ISOTROPY_CALIBRATED essence — angular earned, radial calibrated;
report :5); `c_S` deferred (matter-sector, `kξ≪1`). ⭐ The pathA_32 trip-ups are 016's locus — the aggregate probe
battery intact, the eigenvalues COMPUTED (tautology-hash probe), the sourced `T_Ω` fires FAIL_DIMENSIONAL. ⚠ A
CONSUMED-packet citation guard: a det-identity is blind to an off-diagonal sign flip — but 016 consumes DENSITY SCALARS
(not a matrix), so the dual-site cite is value-level; still, guard each consumed datum individually (a coordinated
both-sites corruption is the escape to close, per the stage010/015 lesson).

## Process (unchanged, calibrated — the per-stage pipeline)
Author the II-G3a reshape directive (§1 the clean 016 slice + §1c the shared-dim-check handling + KEEP-016-only + §2
faithful cover + §3 bridge-strip incl. KEEP-native-`.wl` + transcript-level agreement + §5 the covariance/tautology/dim
teeth with per-tooth ablation + §6 the DEFER-T_Ω register question) → **Codex xhigh design-review** → fold to
`DIRECTIVE_CLEAN` → **⭐ final Grok-4.5 headless compute-verify pass** (assess + independently verify each catch — Grok
is noisy; it compute-verifies the Laplacian/Gram/eigenvalue math) → fold → **Codex confirm-pass on the folds** →
**pre-exec USER GATE** → Codex builds the two scripts (`--sandbox danger-full-access`, background, stdin, xhigh) →
dual-engine both exit 0 (repo root + foreign CWD) → arbiter re-run → full tri-review (fidelity + adversarial-with-
**per-tooth ablation**, ⭐ hunt the aggregate-battery integrity + any vacuous able-to-fail + a mirror-`.wl`) → remediate →
fresh-agent re-verify → bump counts 15→16 → parameter register (⭐ confirm ZERO new 016 knobs + DEFER T_Ω to 017) +
Codex-verify → note/card/`\input{stages/stage_016}` + registration → rebuild PDF → commit + docs/memory sync.
Orchestrator authors notes/cards/LaTeX/registration; Codex codes. Target stem: `ledger_stage016_l2_so3_covariance`
(confirm slug at authoring).
