# II-G3b (ledger_stage017) source map — pathA_32 grouped-P2 lane isotropy (ISOTROPY_CALIBRATED 2/2, the COMPLETING slice)

> Running-start prep captured 2026-07-09 (post stage016, before authoring the II-G3b reshape directive) so the directive
> can be written without re-discovery. **All line refs below VERIFIED against the current sources 2026-07-09**
> (`pathA_32_grouped_p2_isotropy_sympy.py` = 1495 lines; `pathA_32_grouped_p2_isotropy.wl` = 696 lines;
> `reports/pathA_32_grouped_p2_isotropy.md` = 64 lines). Same source pair as stage016 (the pathA_32 2-way split:
> 016 = the EARNED angular covariance theorem; **017 = THIS = the grouped-lane response + calibration partition, COMPLETES
> the joint `ISOTROPY_CALIBRATED` and EXPORTS the ℓ=2 port kernel**).
> Companions: `stage016_pathA32_l2_covariance_source_map.md` (the 016 slice — this map is its complement; the 016↔017 cut is
> documented there in §1 CLEAN CUT + §1b) and `part2_gravity_atomic_split.md` (rows 016/017 + the pathA_32 trip-ups L84–86).
> Build-order id **017**, Part II. Source top-line: **`ISOTROPY_CALIBRATED`** — 017 lands the JOINT (016 was the PARTIAL).

## ⭐ The FOUR headline differences from 016 (READ FIRST)

1. **017 is the COMPLETING slice — it lands the FULL `ISOTROPY_CALIBRATED`.** Where 016 printed `ISOTROPY_CALIBRATED (1/2)
   — SO(3) covariance theorem EARNED (PARTIAL)`, **017 runs the full `verdict_from_gates` (`.py` L280–299) and lands the
   joint** (all 8 gates ∧ non-empty `calibration_inputs`). 017 owns the **response** half of the cut: assemble the grouped
   {20,21,22} lanes → **raw-D defects = 0 (the PRIMARY isotropy test)** → **normalized u-defects = 0 (the CROSS-CHECK)** →
   the stability/denominator guards → the **calibration partition** → its **6-probe aggregate battery** → the joint landing.
2. **⭐ 017 CONSUMES 016's covariance theorem via a GENUINE cross-stage dual-site (contrast 016's provenance-only cite of
   013).** Key pin (1): **016↔017 are the SAME pathA_32 convention** (VOLUME densities on `a²dwdΩ`, dimensionless `β₂`) —
   *unlike* the 013→016 seam (line-vs-volume, `β=L⁻¹` vs dimensionless β₂) where the cross-stage relation was NON-transferable
   and 016 could only cite 013 as PROVENANCE. So here a genuine cross-stage integrity guard on 016's exported `λ_m=6` and the
   K₂ angular-stiffness form `K̃+λ_m·T̃_Ω` IS available and REQUIRED. **017 CITES (does NOT re-derive) 016's Gram=I₅ / λ_m=6 /
   K₂-form**; the `assemble_channel` lane `K2 = c_self·(K̃+λ·T̃_Ω)` rides the cited `λ_m=6` and `c_self=1` (the Gram diagonal, a
   016 fact). 017 does NOT recompute the Laplace–Beltrami / Rayleigh / eigenfunction-residual / full 5×5 Gram — those are
   016's earned derivation. (017 DOES build the harmonics + integrate the **anisotropy** coefficients `∫P₂·Y² dΩ` — its OWN
   response machinery, NOT part of 016's covariance theorem.)
3. **⭐⭐ THE calibration partition lives HERE (the KEY 017 register question).** The single verdict lands
   `ISOTROPY_CALIBRATED` (not `ISOTROPY_PASS`) **because `calibration_inputs` is non-empty** (`.py` L297–298, verdict
   `verdict_from_gates` L280–299; the partition `derived_inputs`/`calibration_inputs` L955–978; report :5, :58–63). ⚠ **017 is
   the stage that carries the calibration partition — so 017 is where the DEFERRED `T_Ω`/`β₂` counting (16 punted) gets
   RESOLVED.** Likely a **calibration-adding stage** (the 3rd Part-II one, after 013's `{μ_η,T_w,β}` + 015's `{Vp0/ℓ_c}`). Do
   NOT pre-decide the count — resolve at the register step + Codex-verify (§6).
4. **⭐ The pathA_32 trip-ups apply to 017's OWN 6-probe battery.** From the rig-history (`part2_gravity_atomic_split.md`
   L84–86, and the 016 trip-ups now on 017's side): (a) **keep the 6-probe aggregate battery intact** — neuter any one probe's
   computed flag → 017's aggregate `able_to_fail_ok` flips false (`.py` `able_to_fail_if_probe_neutered` L1313–1316,
   `neutering_any_probe_flips_false` behavior); (b) each probe = a COMPUTED conjunction (its mutation fires ITS verdict ∧ its
   self-ablation suppresses the fail), **no stamped literals** (per-tooth ablation WILL probe each — check the `computed_fail_gate`
   / `fail_fires` / `fail_suppressed` flags are computed from the real objects, cf. 016's stamped `participates_in_verdict`
   remediation); (c) a forced probe must NOT port foreign lane machinery (the process lesson banked this session).

## §1 The 017 slice (`.py` line ranges) — the CLEAN CUTS (all VERIFIED)

The whole computation lives in `symbolic_engine()` (`.py` L762–1458). **016 owns the angular covariance sub-block + its 3
covariance teeth; 017 owns the grouped-lane RESPONSE + calibration partition + the 6 response probes + the joint landing.**
The 017-owned cuts:

- **Support-scalar symbols — `symbols` L764–774** (adds `B0/Z0/B2/Z2/B4/Z4` tilde — the D-lane support/Maxwell scalars,
  017's D-lane content; 016 used only `M̃`/`K̃`/`T̃_Ω`). `subs_sample` L775.
- **The anisotropy machinery (017's response driver) — L790–806:** `p2_axis_z=(3cos²θ−1)/2` L790; `anisotropy_coefficients =
  ∫ P₂·Y² dΩ` per harmonic L791–794 (⭐ 017's OWN S² integrals — NOT the Gram/λ of 016); `anisotropy_self_overlaps` L795–798;
  `grouped_isotropic`/`grouped_linear_coeff`/`grouped_perturbed` L800–806 (the ε-perturbed grouped response).
- **Support scalars — L824–835** (`M`=`M̃`, `Kbase`=`K̃`, `Tomega`=`T̃_Ω` cited from 016; `B0..Z4`, `S0=B0+Z0`, `S2=B2+Z2`,
  `S4=B4+Z4` — 017's D-lane support).
- **⭐ `assemble_channel` L837–859** — the per-channel lane assembly: `c_self=∫(w·Y²)dΩ` L839 (isotropic weight → the Gram
  diagonal =1, a cited 016 fact), `lam=lambdas[name]` L840 (⭐ **the CONSUMED λ_m=6** — cite, don't re-derive), `pref=c_self·profile`
  L841, `M_lane=pref·M̃` L842, `K_lane=pref·(K̃+λ·T̃_Ω)` L843 (⭐ **rides the cited K₂ angular-stiffness form**), `B_lane`/`Z_lane`
  L844–845, and the **D-lanes** `D0=K−B0−Z0`, `D2=−(M+B2+Z2)`, `D4=−(B4+Z4)` L846–850 (the ℓ=2 port kernel's D content).
- **`ungrouped` L861** (assemble all 5 channels) + **`average_expr` L863–864** + **`grouped_lanes` L866–884** (the grouped
  {20 / 21=avg(21c,21s) / 22=avg(22c,22s)} M2/K2/D lanes — the exported port kernel).
- **The `build_dimensional_check` CALL L885–892** — ⚠ **016-owned** (the T_Ω FAIL_DIMENSIONAL probe; the def is L302–487). The
  call passes `M2_expr=grouped_lanes["20"]["M2"]`, `K2_expr=grouped_lanes["20"]["K2"]`, `lambda_m=6` — but for the "20" lane
  `pref=∫Y₂₀²=1` (Gram diagonal), so those are the **bare core** `M̃` / `K̃+6·T̃_Ω`. **016 built its dim check on the bare core
  directly (016 source map §1c) — so 017 CITES `dimensional_ok` from 016; do NOT rebuild the dim check in 017.**
- **`cs_equal` L894–903** — the c/s degeneracy within-group check (D21c=D21s, D22c=D22s at orders 0/2/4); part of
  `lane_collapse_ok`. **017-owned.**
- **⭐ The raw-D gate (PRIMARY isotropy) — L905–914:** `raw_triples` L905–908 (the {20,21,22} D-lane triple per order),
  `raw_defects` via `defect_pair` L909–913, `raw_defects_zero=all_zero(...)` L914. ⭐ **`defect_pair` (L137–141) = the
  A/B P₂-projection `(2·D20−D21−D22)/10`, `(D21−D22)/2`; raw-D defects = 0 IS the isotropy statement** (report :12). PRIMARY.
- **The normalized-u CROSS-CHECK — L916–927:** `u2_from_d=−D2/D0` L916–917, `u4_from_d=(D2²−D0·D4)/D0²` L919–920, `u2_lanes`/
  `u4_lanes` L922–923, `a2,b2,a4,b4=defect_pair(...)` L924–925, `normalized_defects` L926, `normalized_defects_zero` L927
  (report :13). ⚠ **CROSS-CHECK, not primary** — the honest caveat (report :5): the rung is CALIBRATED because the wall
  profile + radial/support scalars are FROZEN (not derived from the Gate-1 `R0` support equation).
- **`tautology_clear` L929–940** — ⚠ **016-owned** (harmonics distinctness + self-overlap==1); 017 cites it.
- **The stability/denominator guards (017-owned) — L942–953:** `K2_sample`/`M2_sample`/`omega_sample` L942–944 (numeric
  sample at `CALIBRATION_SAMPLE`), `k2_window_min`/`d0_window_min` L945–951 (over `CALIBRATION_WINDOW`), `stability_ok=(M̃_min>0
  ∧ K2_window_min>0)` L952, `denominator_guard_ok=(d0_window_min>0)` L953. ⚠ These are the FROZEN calibration window (report :14).
- **⭐⭐ THE calibration partition — L955–978:** `derived_inputs` L955–963 (the EARNED content: harmonics, 5×5 Gram, eigenvalues,
  K₂ coeff from computed λ_m [= 016's, cited], ungrouped self-overlaps, **raw-D + normalized defect algebra** [017's], probe
  verdicts) vs `calibration_inputs` L964–978 (the FROZEN inputs: `R0(w)`, `β₂(w)`, `μ_η`, `T_w`, `K_η`, `T_Ω`, `M̃`, `K̃`, `T̃_Ω`,
  `B̃0/2/4`, `Z̃0/2/4`, the physical window, the Gate-1 D/N provenance). ⭐ **This partition IS the register question (§6).**
- **`case_verdict` L980–992** (a scoped `verdict_from_gates` with overridable gates — the probe-verdict engine) — 017-owned
  (both engines' probes call it).
- **`d_common` L994–998 / `c_group` L999–1004** — the isotropic common D-lane `{D0=K̃+6T̃_Ω−S0, D2=−(M̃+S2), D4=−S4}` + the
  grouped anisotropy coefficients `{20,21,22}` (used by the response probes). 017-owned.
- **⭐ The 6 response probes (017's aggregate battery) — machinery L1006–1098, instantiations L1127–1146, payloads L1148–1258:**
  - **`pure_prefactor_anisotropy`** — machinery `pure_d_by_lane`/`pure_raw_defects`/`pure_samples`/`pure_u2/u4`/`pure_u_defects`/
    `pure_raw_moves`/`pure_linear_delta` L1006–1025; payload L1149–1161. ⭐ raw-D MOVES (linear in ε) ∧ normalized-u STAYS ZERO →
    `FAIL_ANISOTROPIC_BRANCH` (report :18). The decisive raw-vs-normalized discriminator.
  - **`sector_selective_anisotropy`** — machinery `sector_d_by_lane`/…/`sector_u_moves` L1027–1058; payload L1162–1172. raw-D AND
    normalized-u BOTH move → `FAIL_ANISOTROPIC_BRANCH` (report :19).
  - **`m_dependent_profile`** — machinery `profile_scales`/`profile_d_by_lane`/`profile_raw_defects`/`profile_sample`/
    `profile_raw_moves` L1060–1074; payload L1173–1182. m-dependent radial scale → raw-D moves → `FAIL_ANISOTROPIC_BRANCH`
    (report :20).
  - **`degenerate_beta_zero`** — `beta_stability_probe` def L1083–1098; `degenerate_beta_probe=beta_stability_probe(0)` L1127 +
    ablation `(1)` L1128; payload L1183–1191. `β₂=0 → M2=K2=0 → FAIL_STABILITY` (report :21).
  - **`singular_denominator`** — machinery `singular_subs`/`singular_d0_value`/guards L1076–1081; `singular_verdict` L1131 +
    ablation L1132; payload L1201–1218. `B0+Z0=K2 → D_A0=0 → FAIL_SINGULAR_RESPONSE` (report :23).
  - **`static_drop_inertia`** — `static_wrong_d2=−(B2+Z2)` L1139, `static_dynamic_retained=D2.has(M̃)` L1140, `static_verdict`
    L1141 + ablation L1142–1143; payload L1228–1244. drop `M̃` from D_A2 → `FAIL_STATIC_RESPONSE` (report :25).
  - ⚠ **NOT 017 (these 3 are 016-owned probes — do NOT rebuild in 017; cite 016's aggregate):** `wrong_eigenvalue`
    (`forced_eigenvalue_probe` def L1100–1114, inst L1129–1130, payload L1192–1200), `tautology_hash_collision` (`lane_hash_probe`
    def L1116–1125, inst L1133–1138, payload L1219–1227), `dimensional_corrupt_T_Omega` (inst L1144–1146, payload L1245–1257).
- **⭐ The joint landing (017-owned) — L1260–1330:** `expected_probe_verdicts` L1260–1270 (⚠ 017's has the 6; 016's has the 3),
  `expected_probe_verdicts_match` L1271–1273, `computed_probe_gate_flags` L1274–1308 (⚠ 017 aggregates ITS 6), `able_to_fail_from_flags`
  L1309–1310, `able_to_fail_ok` L1312, `able_to_fail_if_probe_neutered` L1313–1316 (the neuter-one-flips discipline),
  `dynamic_retained` L1318 (= `D2.has(M̃) ∧ not D0.has(M̃)` — 017-owned gate), `baseline_gates` L1320–1329, and the FULL
  **`verdict = verdict_from_gates(baseline_gates, calibration_inputs)` L1330 → `ISOTROPY_CALIBRATED`**.
  - ⚠ `baseline_gates` (L1320–1329) assembles all 8 gates. **In the split, `covariant` L1322 (= gram∧λ∧residuals∧k_coeff) +
    `dimensional_ok` L1321 + `tautology_clear` L1323 are CONSUMED from 016 (cited booleans, guarded by the §1c dual-site);**
    `dynamic_retained` L1324 + `stability_ok` L1325 + `denominator_guard_ok` L1326 + `lane_collapse_ok` L1327 (= `raw_defects_zero
    ∧ cs_equal`) + `able_to_fail_ok` L1328 (017's 6-probe aggregate ∧ 016's cited 3-probe aggregate for the JOINT battery) are
    017-COMPUTED. Then the joint lands. **The 9-probe battery = 016's 3 ∧ 017's 6, reconstituted across the two stages** (the
    BREATHING=013∧014∧015 pattern) — do NOT force 017 to re-run 016's 3 probes; cite 016's aggregate.

## §1b The `.wl` 017 slice (VERIFIED) — the independent Wolfram route (KEEP native, sever ONLY the YAML)

⭐ **The pathA_32 `.wl` is ALREADY a genuinely independent engine** (native `Integrate`/`FullSimplify`, native `D` Laplacian,
native `Hash[...,"SHA256"]`, its OWN `dimOf`/`verdictFromGates`/`caseVerdict`; NO `Get`/`Import` of the `.py`'s expressions).
The reshape KEEPS the native route and severs ONLY the scratch-YAML handoff (§3). 017's `.wl` slice:
- Shared helpers: `defectA`/`defectB` L49–50, `u2FromD`/`u4FromD` L51–52 (native); `p2AxisZ` L95, `anisotropyCoefficients=intS2[p2AxisZ·#²]`
  L96, `cGroup` L97–101, `groupedCoeff` L103 (017's anisotropy machinery).
- `dCommon` L126–130; `lanes` L131; **`assembleChannel` L132–155** (native `intS2`+`FullSimplify` lane assembly — cites `lambdas[name]`
  L136 = the consumed λ); `ungrouped` L156; `avgExpr` L157; **`groupedLanes` L158–176**.
- The dimensional check L178–229 (⚠ **016-owned**; the standalone assert L227–229) — 017 cites `dimensionalOk`.
- **`csEqual` L231–240**; **raw defects** `rawTriples`/`rawDefects`/`rawDefectsZero` L241–250; **normalized** `u2Lanes`/`u4Lanes`/
  `normalizedDefects`/`normalizedDefectsZero` L252–260; **`calibrationInputs` L262–276**.
- `verdictFromGates` L278–289 + `caseVerdict` L290–297 (native, its own); `sampleRules` L299–309 + `evalSample` L311 + sample
  helpers L312–319.
- The response-probe machinery: `pureDByLane`/`pureRawDefects`/`pureU2/U4`/`pureUDefects`/`pureRawMoves`/`pureUZero` L321–345;
  `sectorDByLane`/…/`sectorUMoves` L347–375; `profileScales`/…/`profileRawMoves` L377–390; `betaStabilityProbe` def L393–403;
  (`forcedEigenvalueProbe` L404–414 + `laneHashProbe` L416–422 = 016-owned defs); the singular machinery `k2Sample`/`singularRules`/
  `singularD0Sample`/guards + `singularVerdict`/ablation L428–445; `tautologyProbe`/`Ablation` L446–447 (016); `staticWrongD2`/
  `staticDynamicRetained`/`staticVerdict`/ablation L448–452 (017); `dimensionalProbeVerdict`/`Ablation` L453–454 (016).
- The aggregate: `probeVerdicts` L456–466, `expectedProbeVerdicts` L467–477, `expectedProbeVerdictsMatch` L478, `computedProbeGateFlags`
  L479–489, `ableToFailOk` L490.
- `stabilityOk`/`denominatorGuardOk` L492–493 (017), `tautologyClear` L494–497 (016), `laneCollapseOk` L498 (017), `dynamicRetained`
  L499 (017), `covariantOk` L500 (016), and the FULL `verdict = verdictFromGates[...]` **L501** → `ISOTROPY_CALIBRATED`.
- **The standalone calibrated-baseline assert L503–508** (`If[! ...verdict === "ISOTROPY_CALIBRATED", fail[...]]`) — a
  reshape-ready idiom (print-only / `fail[]` on failure). KEEP + re-target to 017's gates.
- **⚠ The bridge to sever (§3): the full YAML writer L510–696** (`matrixLines` L510+, `numText`/`quoteText` helpers, the
  `yamlOut` assembly, and the `Export` of the scratch YAML `pathA_32_mathematica_results.yaml`). SEVER — print-only + `fail[]`;
  the dual-engine agreement becomes **transcript-level** (both engines print the same raw-D=0 / normalized-u=0 / probe
  verdicts / joint `ISOTROPY_CALIBRATED`; the stage014/016 transcript pattern). KEEP the native route + the standalone assert +
  the arity discipline (def/call scan + unevaluated-leakage transcript scan).

## §1c The consumption of 016 (the GENUINE cross-stage dual-site) — READ

⭐ **Key pin (1): 017 CONSUMES 016's covariance theorem — cite, do NOT re-derive — via a genuine cross-stage dual-site
(same convention, so a checkable relation, unlike 016's provenance-only cite of 013).** 016 exported: `Gram=I₅`, the computed
`λ_m=6`, the K₂ angular-stiffness form `K̃+λ_m·T̃_Ω`, and (implicitly) the covariance gate booleans (`covariant`,
`dimensional_ok`, `tautology_clear`) + its 3-probe aggregate. In 017:
- **What 017 CITES (not re-derived):** `λ_m=6` (used in `assemble_channel` L840/L843), `c_self=1` for the isotropic
  self-overlaps (the Gram diagonal — 016's `Gram=I₅`), the K₂ form, and the three consumed gate booleans (`covariant` /
  `dimensional_ok` / `tautology_clear`) + 016's 3-probe aggregate. **017 does NOT recompute** the Laplace–Beltrami operator, the
  Rayleigh quotient, the eigenfunction residual, or the full 5×5 Gram — those are 016's earned derivation.
- **What 017 DOES compute (its own):** the harmonics (needed for the anisotropy integrals), the anisotropy coefficients
  `∫P₂·Y²dΩ` (its response driver — NOT part of 016's theorem), the lane assembly, the D-lanes, the raw-D + normalized
  defects, the stability/denominator guards, and its 6 response probes.
- **⭐ The REQUIRED genuine dual-site guard (directive to specify; Codex designs the mechanism):** because 016↔017 share
  the SAME convention, a corruption of the cited `λ_m` (or the K₂ form) MUST be CHECKABLE and fire in BOTH engines. Acceptance
  = a ONE-SITE corruption of the cited `λ_m`/K₂-form fires a guard (both engines); the coordinated-both-sites escape is closed
  by an anchor (the stage010/015 lesson). ⚠ **The trap to avoid: `λ` is dimensionless, so a `λ:6→4` corruption is
  dimensionally silent AND leaves stability/denominator windows positive — it may fire NO existing 017 gate.** So the dual-site
  cannot rely on an existing gate incidentally catching it; it needs an EXPLICIT integrity check (e.g. site A = the `λ` used in
  the lane K₂; site B = an independent tie to 016's export — e.g. reconstruct `c_self=∫Y₂₀²=1` freshly and assert the "20"-lane
  bare `K2 = K̃+λ_cited·T̃_Ω` matches 016's exported form, so a one-site `λ` drift is caught). Do NOT hand-design the route in
  the directive — state the REQUIREMENT + acceptance; Codex builds it; the per-tooth-ablation leg proves it able-to-fail.
- ⚠ **`c_S` NOT consumed** (matter-sector deferred, as 013/014/015/016).

## §2 The 017 claim-set (derive + assert; report quotes)

- **Grouped-lane isotropy — raw-D defects = 0 (the PRIMARY isotropy result, EARNED-on-cited-covariance).** With the grouped
  {20,21,22} lanes assembled from the CONSUMED λ_m=6 + K₂ form, the P₂-projection defects vanish: `{a_D0,b_D0,a_D2,b_D2,a_D4,b_D4}=0`
  (report :12). This is the isotropy statement — the ℓ=2 response is angularly degenerate across the three grouped lanes,
  riding 016's SO(3) covariance. The within-group c/s degeneracy (`cs_equal`) + `raw_defects_zero` = the `lane_collapse_ok` gate.
- **The normalized-u CROSS-CHECK — u2/u4 defects = 0.** `{a2,b2,a4,b4}=0` (report :13) confirms the isotropy at the
  normalized-response level (the DtN-observable u₂=−D2/D0, u₄=(D2²−D0D4)/D0²). ⚠ **Cross-check, NOT primary** — a pure-prefactor
  anisotropy MOVES raw-D but leaves normalized-u zero (the `pure_prefactor` probe, report :18), so normalized-u alone would
  MISS it. Keep BOTH; raw-D is decisive.
- **The stability + denominator guards (17-owned, FROZEN-window).** `stability_ok` (M̃>0 ∧ K2>0 over the window) +
  `denominator_guard_ok` (D_A0>0) hold (report :14) on the FROZEN `CALIBRATION_WINDOW`. ⚠ This frozen-ness is the
  `ISOTROPY_CALIBRATED`-not-`PASS` reason (report :5).
- **The 6-probe response battery (EARNED able-to-fail).** Each of `{pure_prefactor / sector_selective / m_dependent /
  degenerate_beta / singular_denominator / static_drop_inertia}` fires its own verdict (`FAIL_ANISOTROPIC_BRANCH` ×3 /
  `FAIL_STABILITY` / `FAIL_SINGULAR_RESPONSE` / `FAIL_STATIC_RESPONSE`; report :18–25) with a self-ablation that suppresses the
  fail (report :30/:32/:34 for the fixed self-ablations), and **neutering any one flips 017's aggregate** (report :40, the
  `neutering_any_probe_flips_false` discipline — 017 keeps this over ITS 6).
- **The calibration partition (THE 017 accounting).** `derived_inputs` (harmonics/Gram/eigenvalues/K₂-coeff [016, cited] +
  raw-D/normalized defect algebra + probe verdicts [017]) vs `calibration_inputs` (R0/β₂/μ_η/T_w/K_η/T_Ω/M̃/K̃/T̃_Ω/B̃/Z̃/window/
  D-N provenance) — report :58–63; the source's input-partition is report :60–61. ⭐ **The counting resolution is §6.**
- **The joint landing (COMPLETES `ISOTROPY_CALIBRATED`).** The full `verdict_from_gates` over the 8 gates (016's 3 cited ∧
  017's 4 computed ∧ the joint 9-probe `able_to_fail`) ∧ non-empty `calibration_inputs` → **`ISOTROPY_CALIBRATED`** (report :3).
  017 PRINTS the joint (not the partial) and states the 016∧017 composition explicitly.

## §3 Reshape cost (the bridge to sever) — cross-script scratch-YAML family, KEEP the native `.wl`

Same family as 016 (the cross-script runtime-YAML reshape, NOT the sympy-expr-import family). The `.py`'s `symbolic_engine()`
is pure/self-contained; `main()` L1461–1495 writes `SYM_YAML` (L1464), reads `MMA_YAML` (L1466), computes
`compute_engine_agreement` (L618–759, L1467), and writes `RESULTS_YAML`/`REPORT_MD`/`FEED_NOTE` (L1483–1487); the `.wl` writes
its scratch YAML via `Export` (~L693). **Reshape = sever ALL file I/O both directions:** drop `main`'s YAML/report/feed writers
+ `yaml_read`/`yaml_write`/`compute_engine_agreement` (`.py`); drop the `Export` + the YAML-line assembly L510–696 (`.wl`). Each
engine → standalone: print-only, `expect_zero`/`expect_bool`-style asserts (`.py` local ledger idioms), `fail[]`/`Exit[1]` on
failure (`.wl` already has `fail[]` + the standalone assert L503–508). **KEEP the `.wl`'s already-independent native route** —
re-target it to assert its OWN raw-D=0 / normalized-u=0 / 6-probe verdicts / joint `ISOTROPY_CALIBRATED`. **Dual-engine
agreement = transcript-level** (stage014/016 pattern). **Zero file I/O.** Arity discipline (standing — def/call scan +
unevaluated-leakage transcript scan; the stage007 lesson).

## §4 Consumed / exported

- **Consumes (cite, genuine cross-stage dual-site per §1c, don't re-derive):**
  - **016's covariance theorem** — `Gram=I₅`, `λ_m=6`, the K₂ angular-stiffness form `K̃+λ_m·T̃_Ω`, and the `covariant` /
    `dimensional_ok` / `tautology_clear` gate booleans + 016's 3-probe aggregate. ⭐ **This IS a checkable cross-stage
    relation (same convention) — a genuine dual-site, not provenance-only.** Cite `stage016`.
  - **The frozen throat packet + Gate-1 D/N provenance** (PROVENANCE, as 016): `μ_η`, `T_w` (ALREADY counted CALIB at 013 —
    cite as provenance, the density-vs-line convention makes the cross-stage relation NON-checkable, per the 016 lesson;
    `K_η=T_w β²` DERIVED R29 but convention-caveated); `R0(w)` linearized reference (stage011/012 Gate-1); the Gate-1 D/N
    boundary provenance (stage012 R28 `BC_DEPENDENT`, stage011 R26).
  - ⚠ **`c_S` NOT consumed** (matter-sector deferred).
- **Exports (⭐ the ℓ=2 PORT KERNEL — this is 017's, not 016's):** the grouped-P2 ℓ=2 port kernel = the grouped `M₂`, the
  angular `K₂` (`= K̃ + 6·T̃_Ω`), the support scalars `B̃/Z̃`, and the **D-lanes** (`D0/D2/D4`) → **018–021** (pathA_33 quadrupole:
  `QUAD_CALIBRATED` — the wall mode) + **022/023** (pathA_34 cross-ℓ ℓ=2 map) + **024** (pathA_43 density-port wall mode). Per
  the cross-stage flow (`part2_gravity_atomic_split.md` L106): "017 exports the grouped-P2 ℓ=2 port kernel + `K_η+2T_Ω` +
  support scalars → 018–021 + 022/023 + 024 (wall mode)." Cite the exact export contract at note-authoring.

## §5 Teeth candidates (017-specific, per-tooth ablation MANDATORY)

1. **⭐ The raw-D-defects=0 PRIMARY isotropy tooth.** `raw_defects_zero` is a genuine `all_zero` over the P₂-projection of the
   grouped D-triples — a non-isotropic lane (or a corrupted grouped average) makes a defect ≠ 0 → `lane_collapse_ok=False` →
   `FAIL_ANISOTROPIC_BRANCH`. Per-tooth: mutate one grouped D-lane, confirm the fail fires AT the lane-collapse gate.
2. **⭐ The `pure_prefactor` raw-vs-normalized discriminator (the decisive probe).** raw-D MOVES (linear in ε, `pure_linear_delta<1e-10`)
   ∧ normalized-u STAYS ZERO → `FAIL_ANISOTROPIC_BRANCH`. Per-tooth: this proves raw-D is the PRIMARY test (normalized-u alone
   misses a pure-prefactor anisotropy). Confirm both the `raw_D_moves` AND `normalized_u_defects_stay_zero` flags are COMPUTED
   from the real sampled defects, not stamped.
3. **The 5 other response probes (each a computed conjunction):** `sector_selective` (raw ∧ u both move), `m_dependent`
   (raw moves), `degenerate_beta_zero` (β₂=0 → FAIL_STABILITY), `singular_denominator` (D_A0=0 → FAIL_SINGULAR_RESPONSE),
   `static_drop_inertia` (drop M̃ → FAIL_STATIC_RESPONSE). Per-tooth: each `computed_fail_gate` reads the REAL object (not a
   literal) ∧ each `self_ablation` genuinely suppresses. ⚠ **Check for the 016-class stamped-literal defect** (016 had a
   stamped `participates_in_verdict` → remediated): confirm none of the `fail_fires`/`computed_fail_gate`/`fail_suppressed`
   flags are hardcoded True/False; each must be a computed comparison.
4. **⭐ The 6-probe aggregate-battery discipline.** Neuter any ONE of the 6 probes' computed flag → 017's aggregate
   `able_to_fail_ok` flips false (`able_to_fail_if_probe_neutered` L1313–1316). ⚠ trip-up (a) — keep the battery intact (a probe
   must NOT become a no-op).
5. **⭐ The genuine cross-stage dual-site on 016's λ_m/K₂-form (§1c).** A one-site corruption of the cited `λ_m` (or the K₂
   form) fires a guard in BOTH engines; the coordinated-both-sites escape is closed by an anchor. ⚠ **The λ-is-dimensionless
   trap: a `λ:6→4` corruption is dimensionally silent AND window-positive** — so the guard must be an EXPLICIT integrity check,
   not an incidental gate. Per-tooth: mutate the cited λ at one site, confirm the guard (not a downstream window) fires.
6. **The joint-landing tooth.** The full `verdict_from_gates` genuinely reaches `ISOTROPY_CALIBRATED` only when all 8 gates ∧
   the non-empty `calibration_inputs` hold; corrupt any 017-computed gate (dynamic/stability/denominator/lane/able) → the
   corresponding `FAIL_*` fires (the gate chain L281–298 is reachable). Confirm each rung is reachable.
7. **Consumed-provenance dual-site** (μ_η/T_w cited from 013 as provenance; R0/D-N from 011/012) + `.wl` arity scan.

⚠ **NOT 017 (do not rebuild as 017 teeth — 016 owns these):** `wrong_eigenvalue` / `tautology_hash_collision` /
`dimensional_corrupt_T_Omega`. 017 CITES 016's aggregate over those.

## §6 Register expectation — ⭐⭐ THE KEY 017 QUESTION (the calibration partition; likely calibration-ADDING; do NOT pre-decide)

Per key pin (2): **017 is the stage that carries the calibration partition — so the DEFERRED `T_Ω`/`β₂` counting (016 punted)
gets RESOLVED here.** The honest pre-read (⚠ CONFIRM at the register step + Codex-verify against the scripts; the source's
input-partition is report :60–61):

- **⭐ 017 is LIKELY a calibration-ADDING stage** (the 3rd Part-II one, after 013 `{μ_η,T_w,β}` + 015 `{Vp0/ℓ_c}`). The rung is
  `CALIBRATED` precisely BECAUSE `calibration_inputs` is non-empty and the wall/radial/support scalars are FROZEN (not derived
  from the Gate-1 `R0` support equation; report :5).
- **Candidate NEW counted CALIB knobs (RESOLVE — do not pre-decide):**
  - **`T_Ω` (the ℓ=2 angular-stiffness density/scale).** Physically a NEW ingredient — the wall's ROTATIONAL stiffness, absent
    from the ℓ=0 breathing sector (013 had only μ_η/T_w/K_η/β). First appeared dimensionally at 016; 016 DEFERRED its count to
    HERE. Likely a genuinely new counted CALIB knob. ⚠ Note the density-vs-line convention (016 lesson): `T_Ω` is a VOLUME
    density (M L⁻³ T⁻²) in pathA_32's convention.
  - **`β₂(w)` (the frozen ℓ=2 radial profile).** Distinct from 013's `β` (the ℓ=0 wall inverse-length, `β=L⁻¹`): `β₂` is a
    dimensionless radial PROFILE SHAPE (pathA_32 convention). Likely a new frozen CALIB input (a function, not a scalar) OR a
    manifestation of the wall geometry — resolve.
  - **The support scalars `B̃0/2/4`, `Z̃0/2/4` (the D-lane support/Maxwell content).** First become load-bearing at 017 (016
    did NOT consume them). Whether they are independent frozen CALIB knobs or manifestations of the D/N port + T_Ω must be
    resolved. ⚠ These are the exported port kernel's support content — their count feeds 018–024.
- **Likely DERIVED manifestations (not counted afresh):** `M̃`/`K̃`/`T̃_Ω` — the radial scalars = `∫(density × β₂²) dw` integrals
  of the frozen densities against the profile (like 013's `K_η=T_w β²` was a manifestation, R29). Record as manifestations with
  their reduction edge.
- **Cited provenance (NOT re-counted):** `μ_η`, `T_w` (counted at 013; density-convention cite, the cross-stage relation
  NON-checkable per the 016 lesson); `K_η=T_w β²` (DERIVED R29, convention-caveated); `R0(w)` (Gate-1, stage011/012); the
  Gate-1 D/N provenance (stage012 R28, stage011 R26).
- **Control/tracked-not-counted:** the `CALIBRATION_WINDOW`/`CALIBRATION_SAMPLE` numeric controls (like 014's numeric
  controls) + `EPS_VALUES`/`DELTA_PROFILE`/`PROBE_TOL` probe knobs.
- **Likely a new structural/reduction edge** recording (i) the calibration partition (frozen wall profile + radial/support
  scalars → `CALIBRATED` not `PASS`) and (ii) the ℓ=2 port-kernel export; plus the named-PENDING reduction debt (the frozen
  β₂/T_Ω/B̃/Z̃ derivation from the Gate-1 `R0` support equation — the sibling of R30/R33). CONFIRM at registration.

⚠ **Do NOT let 017 silently under-count OR over-count** — the DEFERRED `T_Ω`/`β₂` MUST be resolved here (016 explicitly
punted to 017). Resolve each of {T_Ω, β₂, B̃/Z̃} as counted-knob vs manifestation vs provenance, with dims + provenance class,
and Codex-verify (the register verify is the gate that catches an `IMPOSED`-dressed-as-`DERIVED` mislabel that would falsely
shrink the irreducible codimension count).

## Verdict tokens + honest scope

017 carries the **grouped-P2 lane-isotropy component (2/2) of `ISOTROPY_CALIBRATED` — the COMPLETING slice** — and lands the
JOINT: grouped {20,21,22} lanes, **raw-D defects = 0 (PRIMARY)**, normalized-u defects = 0 (CROSS-CHECK), the stability/
denominator guards, the 6-probe response battery, the calibration partition, and the full `verdict_from_gates` →
`ISOTROPY_CALIBRATED`. CONSUMES 016's covariance theorem via a genuine cross-stage dual-site (same convention). Caveats: the
rung is CALIBRATED (not PASS) because the wall profile + radial/support scalars are FROZEN calibration inputs (report :5), not
derived from the Gate-1 `R0` support equation; `c_S` deferred (matter-sector, `kξ≪1`); the 54/5 quadrupole normalization +
outgoing odd-N coefficients + solved nonlinear branch data remain Gate 4/5/6 (`G=GENUINE_BLOCKED`), report :63. ⭐ EXPORTS the
ℓ=2 port kernel (K̃ + angular-K₂ + support scalars + D-lanes) → 018–021 + 022/023 + 024. The pathA_32 trip-ups are 017's own
6-probe battery (aggregate intact, computed-not-stamped flags, no ported foreign lane machinery).

## Process (unchanged, calibrated — the per-stage pipeline)

Author the II-G3b reshape directive (§1 the clean 017 slice + §1b the native-`.wl` KEEP + §1c the 016 consumption / genuine
cross-stage dual-site + §2 faithful cover + §3 bridge-strip incl. sever-YAML + transcript-level agreement + §5 the raw-D-PRIMARY
/ 6-probe / dual-site teeth with per-tooth ablation + §6 the RESOLVE-`T_Ω`/`β₂` register question) → **Codex xhigh design-review**
→ fold to `DIRECTIVE_CLEAN` → **⭐ final Grok-4.5 headless compute-verify pass** (assess + independently verify each catch — Grok
compute-verifies the defect algebra / partition; it caught the 016 volume-vs-line convention mismatch, so watch the T_Ω/β₂
dimensional accounting here) → fold → **Codex confirm-pass on the folds** → **pre-exec USER GATE** → Codex builds the two scripts
(`--sandbox danger-full-access`, background, xhigh) → dual-engine both exit 0 (repo root + foreign CWD) → arbiter re-run → full
tri-review (fidelity + adversarial-with-**per-tooth ablation**; ⭐ hunt the 6-probe aggregate-battery integrity + any vacuous
able-to-fail/stamped-literal + a mirror-`.wl` + the genuine dual-site) → remediate → fresh-agent re-verify → bump counts 16→17 →
parameter register (⭐ RESOLVE the `T_Ω`/`β₂`/support-scalar calibration-partition question + dims + provenance class) +
Codex-verify → note/card/`\input{stages/stage_017}` + registration → rebuild PDF → commit + docs/memory sync. Orchestrator
authors notes/cards/LaTeX/registration; Codex codes. Target stem: `ledger_stage017_grouped_p2_lane_isotropy` (confirm slug at
authoring).
