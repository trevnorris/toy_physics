# II-G5a (ledger_stage022) source map — pathA_34 cross-ℓ fingerprints + Gate-4 non-regression (FAIL_UNDERDETERMINED_NOT_PREDICTIVE 1/2, the EARNED-first slice)

> Running-start prep captured 2026-07-09 (post stage021 = pathA_33 QUAD_CALIBRATED fold COMPLETE, before authoring
> the II-G5a reshape directive) so the directive can be written without re-discovery. **All line refs below VERIFIED
> against the current sources 2026-07-09** by a full structural read + targeted orchestrator re-reads of the 022-slice
> ranges (`pathA_34_cross_l_unification_sympy.py` = **1693 lines**; `pathA_34_cross_l_unification.wl` = **388 lines**;
> `reports/pathA_34_cross_l_unification.md` = 47 lines; `directives/pathA_34_cross_l_unification.md` = 337 lines).
> Companion: `part2_gravity_atomic_split.md` (rows **022/023** L42–43 = the pathA_34 2-way `FAIL_UNDERDETERMINED_NOT_PREDICTIVE`
> split + the pathA_34 trip-ups **L90–92** + the cross-stage flows L104–109 + the ▶ NEXT stage022 entry L463–470) and the
> stage016/018 source maps (the EARNED-FIRST PARTIAL-landing exemplars — 022 follows 016/018's EARNED-FIRST pattern) +
> the stage021 reshape directive (the freshest full-pipeline exemplar; the KEEP-native `.wl` + re-scope-to-local-verdict
> self-ablation discipline 022 reuses). Build-order id **022**, Part II. Source top-line (verbatim, report
> `reports/pathA_34_cross_l_unification.md:3`): **`FAIL_UNDERDETERMINED_NOT_PREDICTIVE`** — 022 lands the 1/2 EARNED-first
> component (the cross-ℓ fingerprints + the Gate-4 non-regression); 023 (the native nullspace underdetermination departure
> + the `Z0_ret/Z1_ret` selector need) DELIVERS the FAIL and completes the joint.
> ⚠ **No prior source map existed for stage022** — this is the FIRST authored artifact of the stage (per the calibrated
> pipeline: source map → reshape directive → Codex→Grok→Codex reviews → USER GATE).

## ⭐ The FIVE headline points (READ FIRST)

1. **⚠ pathA_34 splits 2-WAY (022/023), and 022 is the EARNED-FIRST slice — but the joint headline is a FAIL** (unlike
   016/018 whose joints were CALIBRATED). This is a NEW landing pattern: **022's OWN content genuinely PASSES** (the
   cross-ℓ fingerprints derive correctly + the Gate-4 non-regression holds — the report's `Earned:` list, `:47`), yet the
   token it lands is `FAIL_UNDERDETERMINED_NOT_PREDICTIVE (1/2)` because pathA_34's OVERALL gate fails. **The FAIL is
   DELIVERED by 023** (the native nullspace dim-8 / return-nullity-2 leaves `ε_eff` free → underdetermined; the missing
   Gate-6 `Z0_ret/Z1_ret` selector). 022 = the earned half of a gate that ultimately fails (the honest first-class
   representation — cf. pathA_36/stage003 = earned transverse photons + characterized FAIL_CAUCHY departure). ⚠ **Keep
   the two verdicts DISTINCT** (see §3): (a) 022's SCRIPT passes (exit 0) — the earned content is correct + every tooth
   fires; (b) the physics LANDING LABEL 022 prints is the joint FAIL-token as a **1/2 EARNED-first PARTIAL** (atomic split
   row 022 L42 assigns 022 exactly `FAIL_UNDERDETERMINED_NOT_PREDICTIVE (1/2)`).
2. **⭐ 022 is the CROSS-ℓ GENERALIZATION of 018's ℓ=2-only fingerprint + a NON-REGRESSION on the completed pathA_33
   fold.** Where 018 built the single ℓ=2 outgoing DtN response `Ŷ₂ᵒᵘᵗ=−3/Λ₂` and its `{u₂=1/9, u₄=4/81, v₅=1/27}` +
   χ_Q=+1, 022 builds the **cross-ℓ family** `Ŷ_ℓᵒᵘᵗ=−(ℓ+1)/Λ_ℓ` for ℓ=0,1,2 → the radiative fingerprint
   `{radiative_coeff: 1, 1/2, 1/27}` at leading orders `{ω¹, ω³, ω⁵}` (= `z^(2ℓ+1)`, via `z=aω/c_s`), COMPUTED from the
   spherical Hankel `h_ℓ=j_ℓ±i·y_ℓ` (`.py` L127–187, `.wl` L35–64) — **not typed**. The ℓ=2 leg (`−(ℓ+1)/Λ_ℓ|_{ℓ=2}=−3/Λ₂`)
   REPRODUCES pathA_33's earned quadrupole fingerprint — the **Gate-4 non-regression** (the cross-ℓ machinery does not
   regress the earned ℓ=2 result).
3. **⭐⭐ THE CLEANEST CUT: 022 = pure exterior-wave fingerprints + a FINGERPRINT-LEVEL non-regression — NO port kernel,
   NO 019 prefactor algebra, NO nullspace.** The source's `build_gate4_non_regression` (`.py` L225–274) and `.wl`
   `gate4Ok` (L107–113) BUNDLE the fingerprint non-regression (`fingerprint_matches` L254–259; `.wl` L108–111) WITH a
   re-derivation of **019's** squared-denominator prefactor `P(ω)=D₀N/D^cons²` (`.py` L231–252; `.wl` port kernel L88–91
   + prefactor block L93–105 + `resP0/resP2/resP4` in `gate4Ok` L112). ⭐ **STRIP the 019 prefactor echo** (019 DONE
   `f1c426f9` — do NOT re-litigate `P(ω)=D₀N/D^cons²`); keep 022's non-regression at the **FINGERPRINT level** (the ℓ=2
   leg reproduces pathA_33's `{u₂=1/9, u₄=4/81, v₅=1/27, χ_Q=+1}`). ⭐ **This removes 022's ONLY port-kernel dependency**
   (`build_port_kernel_for` L281–303 supplies `N0_from_port` ONLY to the stripped prefactor P0; the fingerprint_matches
   use only `out2["u2_z"]`…). Result: **022 is self-contained exterior spherical-Hankel algebra** (like 018), needing
   NEITHER the port kernel (→ 023's rank-audit row + 019/020/021 DONE) NOR any nullspace/return machinery (→ 023).
4. **⭐ The 022/023 seam (the shared machinery).** 022 owns `spherical_j`/`spherical_y` + `dtn_branch` +
   `build_fingerprints` + a FINGERPRINT-only non-regression + the `3e_break_gate4` tooth. 023 owns
   `build_transfers`/`build_residuals`/`build_rank_audit`/`build_provenance`/the selector + probes R1/3a/3b/3c/3d/3f/3g/3h.
   ⚠ **The verdict engine (`base_verdict`/`run_gate` L1008–1079), the dimensional checker (L702–896), the ablation
   framework (`ablation` L1082–1099), and the dual-engine compare span BOTH stages** — since 022 and 023 are SEPARATE
   standalone scripts, **022 DUPLICATES only the minimal machinery it needs** (an **022-LOCAL** verdict over just the
   fingerprints + non-regression, and its own DYNAMIC 022-local self-ablation), and does NOT import the joint
   `base_verdict` ladder (which computes 023's nullspace/FAIL_UNDERDETERMINED). Same re-scope discipline as stage021's
   021-local `dimensional_ok` verdict.
5. **⭐⭐ The pathA_34 trip-ups live at 023, but 022 must NOT reintroduce them (part2 L90–92).** (a) **The DEFAULT verdict is
   PASS — `UNDERDETERMINED` must be EARNED from a genuinely-computed native nullspace** (the v1 was tri-review-REJECTED for
   a rigged-to-UNDERDETERMINED zero-padded-constraint + flag-driven probes + a headline-only `.wl`). 022 owns NO nullspace,
   so it is structurally clean of the rigged-nullspace concern — but the directive must state 022 does NOT compute/print a
   FAIL from a rigged mechanism; it lands the EARNED partial and DEFERS the nullspace/underdetermination to 023.
   (b) **Back-solving `ε_eff`/`Z` from residuals is FORBIDDEN** (the `FAIL_TAUTOLOGICAL` firewall, `build_provenance`
   L899–970, is 023's) — 022's fingerprint non-regression must be a **derive-vs-typed** check (022's re-derived ℓ=2
   fingerprint vs pathA_33's earned TYPED values), NEVER a back-solve. (c) The `−(ℓ+1)/Λ_ℓ` shapes stay COMPUTED (mutate
   the Hankel derivation → the coeff changes → the tooth fires), never stamped — the same firewall 018 owned for its `27`.

## §1 The 022 slice (`.py` line ranges) — the CLEAN CUTS (all VERIFIED)

The whole computation is `build_*` helpers assembled by `run_gate` (L1027–1079) + `main` (L1670–1693).
**022 owns the fingerprint core + the fingerprint-level non-regression + the `3e_break_gate4` tooth; 023 owns the
transfers/residuals/rank-audit/provenance/selector + the other probes.** The 022-owned cuts:

- **`spherical_j(lval)` L127–134 + `spherical_y(lval)` L137–144** — the explicit ℓ=0,1,2 spherical Bessel functions
  (hand-built from `sin`/`cos`; NOT `SphericalHankelH1`). ℓ=2: `j₂=(3/z³−1/z)sin z−3cos z/z²`, `y₂=(1/z−3/z³)cos z−3sin z/z²`
  (L132–133, 142–143 — identical to 018's `spherical_j2`/`spherical_y2`). The physical basis of the exterior wave.
- **⭐ `dtn_branch(lval, kind)` L147–187** — the core cross-ℓ DtN log-derivative machinery (`@lru_cache`):
  - `h = j + i·y` (outgoing_hankel1 L151–153) / `h = j − i·y` (incoming_hankel2 L154–156) / `h = j` (standing_j L157–159);
  - **`lam = compact(z·diff(h,z)/h)` L163** — the DtN eigenvalue `Λ_ℓ=z·h_ℓ′/h_ℓ`;
  - **`yout = compact(−(lval+1)/lam)` L164** — `Ŷ_ℓ=−(ℓ+1)/Λ_ℓ` (⭐ the cross-ℓ generalization; the `−(ℓ+1)` is the
    outgoing static DtN slot, `normalization_factor=−(lval+1)` L174 → `−1/−2/−3` for ℓ=0/1/2);
  - `series_order=max(8, 2ℓ+4)` L165; `lam_series`/`y_series` L166–167;
  - **`radiative_power = 2·lval+1` L168** (= 1/3/5 for ℓ=0/1/2 — the ω¹/ω³/ω⁵ leading orders);
  - **`radiative_coeff = compact(y_series.coeff(z, radiative_power)/i)` L169** (the RADIATING part, `/i` extraction);
  - `static = compact(y_series.coeff(z,0))` L180 (=1); `raw_outgoing_order` L183;
  - the ℓ=2 sub-coefficients **`u2_z=coeff(z,2)` L184, `u4_z=coeff(z,4)` L185, `v5_z=coeff(z,5)/i` L186** (the pathA_33
    fingerprint slots — used by the non-regression).
- **⭐ `build_fingerprints()` L190–222** — assembles the cross-ℓ battery (`@lru_cache`):
  - **`expected_radiative = {0:1, 1:1/2, 2:1/27}` L192** (the derive-vs-typed targets; the ℓ=2 slot `1/27` IS pathA_33's
    earned `v₅ᶻ`), **`expected_order = {0:1, 1:3, 2:5}` L193** (the ω¹/ω³/ω⁵ orders);
  - `outgoing`/`incoming` = `dtn_branch(ell, …)` for ℓ∈{0,1,2} L194–195;
  - the per-ℓ matches L197–206: `ell{ℓ}_normalization` (=−(ℓ+1)), `ell{ℓ}_static` (=1), **`ell{ℓ}_radiative_coeff`
    (`bool_zero(radiative_coeff_z − expected_radiative[ℓ])` L200–202 — the CROSS-ℓ fingerprint, DERIVED-vs-typed)**,
    `ell{ℓ}_raw_order` (=expected_order[ℓ]), `ell{ℓ}_incoming_flips_radiative_sign` (`bool_zero(incoming+outgoing)` L204–206);
  - the ℓ=2 non-regression sub-checks **`ell2_u2=bool_zero(u2_z−1/9)` L207, `ell2_u4=…−4/81` L208, `ell2_v5=…−1/27` L209**,
    **`chi_q=compact(v5_z/(1/27))` L210, `chi_Q=bool_zero(chi_q−1)` L211** (χ_Q=+1 outgoing);
  - `ok = all(matches.values())` L220, `chi_Q=chi_q` L221.
- **⭐ `build_gate4_non_regression(fingerprints, port, mutation)` L225–274 — ⚠ FOLD ONLY THE FINGERPRINT HALF (L254–260 +
  L262, L271–273); STRIP THE 019 PREFACTOR ECHO (L229–252, L262–270 prefactor keys).** The 022-owned content:
  - `out2 = incoming[2] if mutation.break_gate4 else outgoing[2]` L228 (the `break_gate4` mutation switches the ℓ=2
    branch to incoming → the fingerprint sign flips — the non-regression tooth);
  - **`fingerprint_matches = {u2: bool_zero(out2["u2_z"]−1/9), u4: …−4/81, v5: …−1/27, chi_Q: bool_zero(v5_z/(1/27)−1)}`
    L254–259** — the Gate-4 NON-REGRESSION (the ℓ=2 leg reproduces pathA_33's earned `{1/9, 4/81, 1/27, χ_Q=+1}`);
  - `chi_Q=compact(out2["v5_z"]/(1/27))` L272, `branch_used=out2["kind"]` L262.
  - ⚠ **STRIP (019 DONE, `f1c426f9` — cite as PROVENANCE, do NOT re-derive):** `n0_eff=port["N0_from_port"]` L229, the
    `Nomega`/`Dcons`/`correct_obj=D0·Nomega/Dcons²`/`plain_obj` L230–235, the series `P0/P2/P4` L236–241, `expected`
    `(D0·N2−2·D2·n0)/D0²` etc L242–251, `residuals`/`prefactor_matches` L252–253, the prefactor-object keys L263–270.
    Stripping this ALSO removes `build_port_kernel`/`build_port_kernel_for` (L277–303) from 022 (their `N0_from_port`
    only fed the stripped P0). 022's non-regression `ok` = `all(fingerprint_matches.values())` ONLY.
- **The one 022 able-to-fail probe (inside `build_counterfactuals` L1102–1244):**
  - **`3e_break_gate4` L1179–1186** — `break_gate4=True` → incoming ℓ=2 → the fingerprint radiative sign flips
    (`v₅ᶻ→−1/27`) → the non-regression fails → `expected_fail="FAIL_QUAD_REGRESSION"`. Routes through the DYNAMIC
    `ablation` helper L1181 — ⚠ **RE-SCOPE to an 022-LOCAL verdict** (§3; the source `ablation` re-runs the joint
    `run_gate`/`base_verdict` = 023's ladder, which 022 does NOT build). (The prefactor half of the source's `break_gate4`
    — `plain_obj` instead of `D0·N/Dcons²` — is 019's; 022's `3e` is fingerprint-only.)
- **Shared helpers 022 uses (NOT cut boundaries):** `compact` L84, `series_no_o` L118–124, `bool_zero` L108–113. Symbol
  decls L44–63 (022 uses `z, omega, a, c_s`; NOT `K0c, K_eta, T_Omega, Z0ret, Z1ret, N2..D4, OmegaU..Sport, q_free` = 023's).

- **⭐ CLEAN CUT — 022 owns L127–222 (fingerprint core) + L254–260/262/271–273 (the fingerprint non-regression, the
  FINGERPRINT half of `build_gate4_non_regression`) + the `3e_break_gate4` probe L1179–1186 + the `ablation` helper
  L1082–1099 (re-scoped). It touches NONE of 023's territory. Do NOT pull in 023 (or the already-DONE 019 prefactor):**
  - **019 (DONE `f1c426f9`) — STRIP from `build_gate4_non_regression`:** the prefactor `P(ω)=D₀N/D^cons²` re-derivation
    L229–252 + prefactor keys L262–270 + `build_port_kernel`/`build_port_kernel_for` L277–303. Cite 019 DONE.
  - **023:** `selector_equations`/`selector_subs`/`selector_provenance` L306–358; `build_transfers` L361–461 (T_dc, ε_eff);
    `build_residuals` L464–508 (A0/A1 vs pathA_29); `GENERATOR_DOFS`/`build_rank_audit` L511–699 (the 11-dof nullity-8 /
    return-nullity-2 audit — the FAIL headline); the dimensional checker L702–896 (023's `3f`); `build_provenance` L899–970
    (the FAIL_TAUTOLOGICAL firewall); `detect_decoupling` L973–1005; `base_verdict`/`run_gate` L1008–1079 (the joint ladder);
    probes `R1` L1111–1123, `3a` L1124–1131, `3b` L1132–1162, `3c` L1163–1170, `3d` L1171–1178, `3f` L1187–1203,
    `3g` L1204–1215, `3h` L1216–1223.
  - **⚠ 022 CITES pathA_33's completed fold (018∧019∧020∧021) as the non-regression reference + 008's raw amplitudes +
    009/010's bulk mode as PROVENANCE (§4).** Do NOT let 022 re-present the prefactor / `54/5` magnitude / dim closure /
    nullspace accounting.

## §1b The `.wl` 022 slice (VERIFIED) — the independent Wolfram route (KEEP native, sever ONLY the YAML)

⭐ **The pathA_34 `.wl` is ALREADY a genuinely independent engine** (structural read + orchestrator re-read confirmed:
native `j0..y2` hand-built from `Sin`/`Cos` L35–40, native `lam=z D[h,z]/h` L44, native `-(ell+1)/lam` L45, native
`serZ`/`serW` `Series` L32–33, native `Coefficient` — it does **NOT** `Get`/`Import` the `.py`; the ONLY bridge is the
scratch-YAML `Export` at L385). So — LIKE 018/019/021 (KEEP-native + sever-YAML), UNLIKE 020 (which needed authoring) —
the reshape KEEPS the native route and severs ONLY the YAML. 022's `.wl` slice:
- `serZ`/`serW` L32–33 (native z/ω series); `j0..y2` L35–40 (native ℓ=0/1/2 spherical Bessel).
- **⭐ `branchData[ell, h, kind]` L42–64** — the native cross-ℓ DtN core: `lam=FullSimplify[z D[h,z]/h]` L44, `yout=-(ell+1)/lam`
  L45, `pow=2 ell+1` L48, `rad=Coefficient[yser,z,pow]/I` L49, `normalizationFactor=-(ell+1)` L53, `static=Coefficient[yser,z,0]`
  L56, `u2z=Coefficient[yser,z,2]` L60 / `u4z=…z,4]` L61 / `v5z=…z,5]/I` L62.
- `out0/out1/out2` L66–68, `in0/in1/in2` L69–71 (outgoing hankel1 / incoming hankel2).
- **`fingerprintOk` L73–86** — the cross-ℓ fingerprint asserts: `normalizationFactor==-1/-2/-3` L74–76, `static==1` L77–79,
  **`radiativeCoeff==1 / 1/2 / 1/27` L80–82**, `incoming==-outgoing` L83–85 (derive-vs-typed).
- **The Gate-4 non-regression (the FINGERPRINT half of `gate4Ok`):** `chiQ=out2["v5z"]/(1/27)` L106, **`out2["u2z"]==1/9`
  L108, `out2["u4z"]==4/81` L109, `out2["v5z"]==1/27` L110, `chiQ==1` L111** — the ℓ=2 reproduces pathA_33.
- **⚠ 019 territory in the `.wl` (STRIP — 019 DONE):** the port kernel `Pport`/`DeltaPort`/`N0port`/`P0port` L88–91 + the
  prefactor block `Nomega`/`Dcons`/`prefObj`/`p0`/`p2`/`p4`/`expectedP*`/`resP*` L93–105 + the `resP0==0 && resP2==0 &&
  resP4==0` conjuncts in `gate4Ok` L112. **022's reshaped non-regression guard = `out2` u2z/u4z/v5z/chiQ (L108–111) ONLY.**
- **⚠ 023 territory in the `.wl` (EXCLUDE):** `K1dc`/`T0dc`/`T1dc`/`eps*`/`A*lead`/residuals L115–126; `rankDofs`/`rankAuditFor`/
  `rankBaseline`/`rankSelector` L128–149; the verdict machinery `positiveBoundedQ`…`gateVerdictFor` L151–177; the dimensional
  checker L179–265; the verdict probes + `headlineOk` L267–298.
- **⚠ The bridge to sever (§3): the YAML writer L300–388** (`Export[yamlOut,…]` L385, the YAML-line assembly, the
  `scratchDir`/`yamlOut` setup L15–22). SEVER: no scratch YAML, print-only + `fail[]` on failure (the `.wl` has `fail[msg_]`
  L4). **Dual-engine agreement = transcript-level** (both engines print the SAME derived `radiativeCoeff = 1 / 1/2 / 1/27`,
  `normalizationFactor = -1/-2/-3`, `static = 1`, `u2z=1/9`/`u4z=4/81`/`v5z=1/27`/`chiQ=+1`, incoming-flips-sign; the
  stage018 transcript pattern). **Zero file I/O.** Arity discipline (standing — def/call scan + unevaluated-leakage
  transcript scan; the `.wl` has `Module`s in `branchData`).

## §1c The consumption resolution (READ — the non-regression = a checkable derive-vs-typed; the rest provenance)

⭐ **UNLIKE 018 (pure provenance, no cross-stage relation to check), 022 HAS a checkable cross-stage relation: the Gate-4
non-regression.** But the OTHER inputs are provenance-only (structural read §6-confirmed). Resolve each precisely:

- **pathA_33's completed fold (018∧019∧020∧021) — the NON-REGRESSION reference (a genuine derive-vs-typed check, NOT a
  provenance-only cite and NOT a theatrical dual-site).** 022 re-derives the ℓ=2 outgoing fingerprint from its OWN Hankel
  `h₂` (via `−(ℓ+1)/Λ_ℓ|_{ℓ=2}=−3/Λ₂`) and checks it reproduces pathA_33's earned TYPED values `{u₂=1/9, u₄=4/81, v₅=1/27,
  χ_Q=+1}` (`.py` L254–259, typed at L207–209/L255–258; `.wl` L108–111). ⭐ **This IS the emitted-vs-CHECKED test** (per
  018/019/021's discipline): the pathA_33 earned values GENUINELY enter 022's non-regression assert (CHECKED), so the
  `3e_break_gate4` tooth is genuine (mutate the ℓ=2 branch → the derived coeff changes → the match fails). ⚠ **The "second
  site" is a TYPED pathA_33-EARNED value, NOT a re-derivation from the SAME `−(ℓ+1)/Λ_ℓ` machinery** (that would be the 017
  subsumed-dual-site trap / X≡X). The earned VALUE of the non-regression: the cross-ℓ `−(ℓ+1)` ℓ-dependence, specialized to
  ℓ=2, gives the RIGHT `−3` slot → reproduces the quadrupole; a corruption `−(ℓ+1)→−ℓ` gives `−2/Λ₂` at ℓ=2 → wrong
  fingerprint → fires. Cite `stage018`/`stage019`/`stage020`/`stage021` (the COMPLETE QUAD_CALIBRATED fold, per the
  cross-stage flow `part2` L107 "018–021 export the Λ₂ fingerprint + χ_Q=1 + the 54/5 partition → 022 (non-regression)").
  ⚠ Cite the prefactor (019) / `54/5` magnitude (020) / μ̂₀-free dim closure (021) as DONE — 022's non-regression is
  **FINGERPRINT level ONLY** (the earned exterior signature), NOT the `54/5` magnitude re-derivation (020's CALIBRATED,
  not re-litigated) and NOT the prefactor algebra (019 DONE, STRIPPED per §1/§1b).
- **008's raw ℓ=0/1/2 outgoing amplitudes (PROVENANCE only).** 008 exports the raw DtN outgoing amplitudes ℓ=0/1/2 (orders
  p=1/3/5) + the `R0=−M0`/`R1=−D1` cancellation targets. 022's cross-ℓ fingerprints are the DtN responses `−(ℓ+1)/Λ_ℓ`,
  **self-derived** from the spherical Hankel `h_ℓ` (like 018 self-derived `h₂` from `j₂/y₂`) — 022 does NOT read 008's
  amplitude VALUES. In the source, `R0`/`R1` appear only as typed provenance tags (023's `build_provenance` L926–929), never
  ingested. Cite `stage008` as PROVENANCE (the ℓ=0/1/2 channel structure whose exterior DtN responses 022 computes). ⚠ Do
  NOT build a theatrical dual-site on 008's amplitudes (022's fingerprints are self-contained).
- **009/010's bulk Helmholtz mode (PROVENANCE only).** The exterior ℓ=0/1/2 outgoing solutions' bulk companion. 022
  reconstructs `h_ℓ` self-contained from `j_ℓ`/`y_ℓ`. Cite `stage009`/`stage010` as PROVENANCE.
- **`c_s` (the density sound speed) — PROVENANCE / units carrier (as at 018).** The cross-ℓ fingerprints are dimensionless
  z-space rationals; `z=aω/c_s` is the units-restoring dictionary that realizes the radiative order `z^(2ℓ+1)→ω^(2ℓ+1)`
  (`.py` L120–121 in 023's residual leg; 022's fingerprints stay abstract-z). Cite `stage005` R1 (`c_s²=5Kρ⁴/m`) as the
  PROVENANCE of what `c_s` IS (NOT a consumed value; the fingerprints are c_s-free). ⚠ Distinct from the frozen-wall `c_S`
  (011–017). **022 likely carries NO units-restored dim leg** (the fingerprints are dimensionless; the physical `z=aω/c_s`
  realization + the residual dim gate `3f` are 023's) — confirm at directive (018 restored `u₂=a²/9c_s²` with a dim leg;
  022 may keep the fingerprints in z-space only, deferring the ω-order units to 023's residuals).

## §2 The 022 claim-set (derive + assert; report/directive quotes)

- **The cross-ℓ outgoing DtN fingerprints (EARNED — the headline).** `Ŷ_ℓᵒᵘᵗ=−(ℓ+1)/Λ_ℓ`, `Λ_ℓ(z)=z·h_ℓ⁽¹⁾′(z)/h_ℓ⁽¹⁾(z)`,
  `z=aω/c_s`, series-expanded about z=0 → the radiative fingerprint `{ℓ=0: 1, ℓ=1: 1/2, ℓ=2: 1/27}` at leading orders
  `{ω¹, ω³, ω⁵}` (= `z^(2ℓ+1)`), with static slot `1` and normalization `−(ℓ+1)` (report `:47` "the `-(ell+1)/Lambda_l^out`
  fingerprints, raw radiative orders"; report `:5`). ⭐ **COMPUTED from the spherical Hankel `h_ℓ`, not typed** — the derived
  `radiative_coeff` is checked derive-vs-typed against `expected_radiative` (L200–202). The incoming `h_ℓ⁽²⁾` flips ONLY the
  radiative sign (`incoming_flips_radiative_sign` L204–206) → the χ_Q sign classification carries cross-ℓ.
- **The Gate-4 quadrupole non-regression (EARNED — the second headline; report `:5`, `:47`).** The ℓ=2 leg
  (`−(ℓ+1)/Λ_ℓ|_{ℓ=2}=−3/Λ₂`) reproduces pathA_33's earned outgoing fingerprint `{u₂=1/9, u₄=4/81, v₅=1/27, χ_Q=+1}`
  (`.py` L254–259; `.wl` L108–111) — the cross-ℓ machinery does NOT regress the earned quadrupole. ⭐ A derive-vs-typed
  check (022's re-derived ℓ=2 vs pathA_33's earned typed values, §1c), NEVER a back-solve (the FAIL_TAUTOLOGICAL firewall,
  023's, must not be reintroduced here). Breaking the ℓ=2 branch (incoming) ⇒ `FAIL_QUAD_REGRESSION` (probe `3e`).
- **The 022-scoped landing (EARNED-first PARTIAL component of a FAIL joint).** Land at the 022 component:
  `FAIL_UNDERDETERMINED_NOT_PREDICTIVE (1/2) = the cross-ℓ −(ℓ+1)/Λ_ℓ fingerprints for ℓ=0,1,2 (radiative {1, 1/2, 1/27} at
  orders {ω¹, ω³, ω⁵}) DERIVED + the Gate-4 quadrupole non-regression (ℓ=2 reproduces the COMPLETE pathA_33
  QUAD_CALIBRATED fold) — the EARNED-first slice; the native nullspace underdetermination departure that DELIVERS the FAIL
  (dim-8 / return-nullity-2, the missing Gate-6 Z0_ret/Z1_ret selector) = 023.` ⭐ **022's SCRIPT passes (exit 0)** — the
  earned content is correct + every tooth fires; the FAIL-token is the physics LANDING label of the joint gate, carried as a
  1/2 PARTIAL (the report's `Earned:` list `:47` is exactly 022; the `Deferred:` list is 023). Do NOT print the joint as
  complete / as an unqualified FAIL (that is 023's departure), and do NOT re-present 023's nullspace/residual/selector
  accounting or 019/020/021's prefactor/magnitude/dim accounting (cite them per §4). ⭐ Follows 016/018's EARNED-FIRST
  PARTIAL-landing pattern (016 printed `ISOTROPY_CALIBRATED (1/2) — SO(3) covariance theorem EARNED (PARTIAL)`; 018 printed
  `QUAD_CALIBRATED (1/4) — outgoing ℓ=2 DtN Hankel fingerprint EARNED (PARTIAL)`).

## §3 Reshape cost (the bridge to sever) + the shared-machinery duplication + the 022-LOCAL verdict

Same family as pathA_30–34 / 018–021 (the cross-script scratch-YAML reshape, NOT the sympy-expr-import family). No argparse
(the three-run protocol is driven by the presence of the MMA scratch file, `.py` L1674–1678). **Reshape = sever ALL file I/O
both directions:** drop `main`'s YAML/report writers (L1670–1693) + `yaml_read`/`yaml_write` (L100–107 region) +
`compare_engines`/`engine_summary`/`build_final_payload`/`build_report` (`.py`); drop the `Export` + the YAML-line assembly
L300–388 (`.wl`). Each engine → standalone: print-only, `expect_zero`/`bool_zero`/`expect_bool`/`expect_fail`-style asserts
(`.py` local ledger idioms), `fail[]`/`Exit[1]` on failure (`.wl` has `fail[msg_]` L4). **KEEP the `.wl`'s already-independent
native route** (§1b). **Dual-engine agreement = transcript-level** (stage018 pattern). **Zero file I/O.** Arity discipline.

**⭐ The shared-machinery duplication (the KEY reshape decision for 022).** The source's verdict engine
(`base_verdict`/`run_gate` L1008–1079), dimensional checker (L702–896), and `ablation` helper (L1082–1099) span BOTH 022 and
023. Since 022 and 023 are SEPARATE standalone ledger scripts, **022 duplicates ONLY the minimal machinery its earned
content needs**:
- **⭐ An 022-LOCAL verdict — NOT the joint `base_verdict` ladder.** 022's local gate = `all(cross-ℓ fingerprints correct)`
  ∧ `all(fingerprint non-regression matches)` → an EARNED token (e.g. `CROSS_L_FINGERPRINT_OK`); breaking the ℓ=2 branch →
  `FAIL_QUAD_REGRESSION`; corrupting a fingerprint derivation → `FAIL_FINGERPRINT`. ⚠ **022 does NOT import
  `base_verdict`/`run_gate`** (those compute 023's nullspace + emit `FAIL_UNDERDETERMINED_NOT_PREDICTIVE`, which is NOT
  built here) — same discipline as stage021's 021-local `dimensional_ok` verdict. The joint FAIL-token is the printed
  LANDING label (a string, the 1/2 partial), NOT the script's internal pass/fail gate.
- **⭐ A DYNAMIC 022-local self-ablation** (the `ablation` helper L1082–1099 re-targeted to the 022-LOCAL verdict): a
  two-verdict re-run (`with_mutation=FAIL_QUAD_REGRESSION` ≠ `without_mutation=CROSS_L_FINGERPRINT_OK`, `rerun_gate_logic=True`),
  NOT a constant (the pathA_34 v1 trip-up + the 018/021 self-ablation discipline). Neuter the mutation → `with==without` →
  the self-ablation flags NOT-able-to-fail.

**Acceptance (dual-engine, both exit 0, CWD-independent):**
- Run each engine from the **repo root** AND from a **foreign CWD** (e.g. `/tmp`), both print-only, both exit 0, no files
  written (verify with a `find` for new files; a `.wl` `Export` slip is the classic leak).
- Both engines emit the same fingerprint transcript: `radiativeCoeff = 1 / 1/2 / 1/27` at orders `ω¹/ω³/ω⁵`,
  `normalizationFactor = -1/-2/-3`, `static = 1`, incoming-flips-sign; the ℓ=2 non-regression `u2z=1/9`/`u4z=4/81`/`v5z=1/27`/
  `chiQ=+1`; `3e_break_gate4` `FAIL_QUAD_REGRESSION`.
- **All able-to-fail teeth fire at their own assert** (per-tooth ablation, §5) — the cross-ℓ fingerprint-derivation teeth
  (mutate the Hankel / the `−(ℓ+1)` ℓ-dependence → the coeff changes), the non-regression `3e` tooth, the χ_Q sign, the
  order/scaling, the DYNAMIC 022-local self-ablation; the `.wl` native-KEEP arity scan.

## §4 Consumed / exported

- **Consumes:**
  - **⭐ pathA_33's COMPLETE fold (018∧019∧020∧021) — a CHECKABLE derive-vs-typed NON-REGRESSION** (§1c): 022's re-derived
    ℓ=2 fingerprint reproduces pathA_33's earned `{u₂=1/9, u₄=4/81, v₅=1/27, χ_Q=+1}` (the typed target = pathA_33's earned
    values; a genuine `3e`-tooth-backed check, NOT provenance-only, NOT a theatrical dual-site). Cite `stage018` (the
    fingerprint) / `stage019` (prefactor, DONE) / `stage020` (`54/5`, DONE) / `stage021` (μ̂₀-free dim closure, DONE). ⚠ 022's
    non-regression is FINGERPRINT level ONLY — do NOT re-derive the prefactor / `54/5` / dim closure.
  - **008's raw ℓ=0/1/2 outgoing amplitudes** — cite `stage008` as PROVENANCE (the channel structure; 022 self-derives the
    DtN responses). No dual-site.
  - **009/010's bulk Helmholtz mode** — cite `stage009`/`stage010` as PROVENANCE. 022 reconstructs `h_ℓ` self-contained.
  - **`c_s`** — cite `stage005` R1 as the PROVENANCE of the units symbol (NOT a consumed value; fingerprints c_s-free). ⚠
    Distinct from the frozen-wall `c_S` (011–017).
- **Exports (→ 023 + Part VII):** the cross-ℓ `−(ℓ+1)/Λ_ℓ` fingerprints + the ℓ=0/1 radiative coefficients `{1, 1/2}` at
  orders `{ω¹, ω³}` (→ **023**'s residual amplitudes `A0/A1`, which realize `z^(2ℓ+1)→(aω/c_s)^(2ℓ+1)` and multiply by the
  return transmission `(1−T_ℓ)`) + the ℓ=2 non-regression result (the earned quadrupole survives the cross-ℓ generalization)
  → **023** (the nullspace departure builds ON the earned fingerprints) + the Part-VII cross-ℓ consistency record. Per the
  cross-stage flow (`part2` L104–108: "008 exports … raw amplitudes → 009/010/023/026"; "018–021 export the Λ₂ fingerprint
  … → 022 (non-regression)"). Cite the exact export contract at note-authoring.

## §5 Teeth candidates (022-specific, per-tooth ablation MANDATORY — mutate the named object, confirm exit-1 AT its own assert)

1. **⭐⭐ The cross-ℓ fingerprint-derivation teeth (`radiative_coeff = {1, 1/2, 1/27}` at orders `{ω¹, ω³, ω⁵}`, normalization
   `−(ℓ+1)`).** The derived `outgoing[ell]["radiative_coeff_z"]` (`.py` L169, L200–202; `.wl` L49, L80–82) is genuinely
   series-computed from the real `h_ℓ`; the check is derive-vs-typed (independent). ⚠ **Per-tooth: mutate the DERIVATION**
   (corrupt the `j_ℓ`/`y_ℓ` expression, or the `−(ℓ+1)` ℓ-dependence in `dtn_branch` L164 → e.g. `−ℓ/Λ_ℓ` → ℓ=2 gives the
   wrong coeff, or a wrong Hankel order) → the derived coefficient changes → the match fails AT its own assert. This is the
   central EARNED tooth; the derivation-from-`−(ℓ+1)/Λ_ℓ` must be EMITTED (the firewall: forbidden to hardcode `{1, 1/2,
   1/27}` and "check" against them), not just the compare.
2. **⭐ The Gate-4 non-regression tooth (`3e_break_gate4`, `FAIL_QUAD_REGRESSION`).** Break the ℓ=2 branch (incoming) → the
   ℓ=2 radiative sign flips (`v₅ᶻ→−1/27`) → the non-regression against pathA_33's `{1/9, 4/81, 1/27, χ_Q=+1}` fails →
   `FAIL_QUAD_REGRESSION` (`.py` L228, L254–259, probe L1179–1186; `.wl` L108–111). Per-tooth: mutate the ℓ=2 branch → the
   022-local verdict fires `FAIL_QUAD_REGRESSION` at its own assert; the correct outgoing ℓ=2 does NOT. ⚠ A derive-vs-typed
   check (022's re-derived ℓ=2 vs pathA_33's earned typed values), NEVER a back-solve of `ε_eff`/`Z` (the FAIL_TAUTOLOGICAL
   firewall, 023's, must not be reintroduced).
3. **⭐ The χ_Q=+1-outgoing / incoming-flips-sign tooth.** The incoming `h_ℓ⁽²⁾` flips ONLY the radiative sign
   (`incoming_flips_radiative_sign` L204–206; `.wl` L83–85) → for ℓ=2, `χ_Q=−1`. Per-tooth: mutate the branch
   (outgoing→incoming) → the sign flips AND a typed `χ_Q=1` cannot survive (the sign is COMPUTED from `j_ℓ±i·y_ℓ`). ⚠ Confirm
   `chi_Q` is COMPUTED (`.py` L210–211/L258; `.wl` L106/L111), never string-typed.
4. **The order/scaling tooth (`raw_outgoing_order = 2ℓ+1 = {1, 3, 5}`).** The radiative power is COMPUTED
   (`radiative_power=2·lval+1` L168; `.wl` `pow=2 ell+1` L48). Per-tooth: mutate the order (e.g. read `coeff(z, 2ℓ+2)`) → the
   `raw_order` match against `expected_order` fails. Reads the REAL series, not a typed tuple.
5. **The DYNAMIC 022-local self-ablation (`3e`) — NOT a constant.** RE-SCOPE `ablation` to the 022-LOCAL verdict over 022's
   fingerprint + non-regression gate only, NOT the joint `base_verdict` (which computes 023's nullspace, NOT built here).
   Two-verdict re-run (`with_mutation=FAIL_QUAD_REGRESSION` ≠ `without_mutation=CROSS_L_FINGERPRINT_OK`,
   `rerun_gate_logic=True`, `.py` L1092–1096). Per-tooth: neuter the mutation → `with==without` → the self-ablation flags
   NOT-able-to-fail.
6. **Provenance-cite integrity (light) + the `.wl` native-KEEP + arity.** pathA_33 fold cited as the CHECKABLE non-regression
   reference (§1c); 008 amplitudes / 009-010 bulk mode / `c_s` R1 cited as PROVENANCE; the `.wl` GENUINELY computes 022's
   fingerprint block (native `branchData`/`serZ`/`Coefficient`), a truly independent Wolfram route (KEPT native), NOT a `.py`
   mirror. Def/call arity scan + unevaluated-leakage transcript scan (stage007 lesson — the `.wl` has `Module`s in `branchData`).

⚠ **NOT 022 (do not rebuild as 022 teeth — 023 owns these, and 019 the prefactor):** `R1_port_kernel_dependency` (023 — the
port kernel), `3a_decouple_knobs`/`FAIL_DECOUPLED` (023), `3b_null_direction`/`inject_null` (023 — the nullspace detector +
the real selector flip), `3c_wrong_sign_antilocalizing`/`FAIL_EPSILON_MISMATCH` (023), `3d_perfect_return`/`FAIL_OVERCANCEL`
(023), `3f_corrupt_dimension`/`FAIL_DIMENSIONAL` + the free-carrier `q_free` control (023 — the residual dim gate),
`3g_assert_not_derive`/`FAIL_TAUTOLOGICAL` (023 — the firewall), `3h_no_consistent_return` (023). ⚠ **The squared-denominator
prefactor `P(ω)=D₀N/D^cons²` (`.py` L229–252; `.wl` L88–105) is 019's DONE algebra — STRIP it from 022's non-regression
(§1/§1b); the port kernel is 023's rank-audit row + 019/020/021's DONE provenance.**

## §6 Register expectation — ⭐ THE KEY 022 QUESTION (likely ZERO new counted knobs; CONFIRM)

Per headline #3 (the pure exterior-wave cut) + the split: **022 is the EARNED cross-ℓ fingerprint + non-regression slice;
the stiffnesses/return admittances that would introduce knobs (`K0c`, `K_eta+2·T_Omega`, `Z0_ret`, `Z1_ret`) all live in
023's residual/transfer/nullspace legs, NOT 022.** So the honest pre-read (⚠ CONFIRM at the register step + Codex-verify
against the scripts):

- **⭐ 022 likely adds ZERO new counted knobs** (like 016/018 / 011/012/014 — an EARNED/structural slice). The cross-ℓ DtN
  fingerprints are pure exterior-wave math (dimensionless z-space rationals from `−(ℓ+1)/Λ_ℓ` + the incoming-flips-sign);
  the non-regression cites pathA_33's already-earned/counted content. **After stripping the 019 prefactor echo, 022 does NOT
  even reference the port kernel or the ℓ=0/1 stiffnesses** — so it introduces NO calibration. `c_s` is R1-DERIVED (cited
  PROVENANCE, a units carrier — NOT a new knob, §1c); `a` is the `CONV` pin (R2-family).
- **⭐ Likely a new STRUCTURAL edge (call it R41 — confirm the next free number; R40 was 021's):** the cross-ℓ outgoing-DtN
  `−(ℓ+1)/Λ_ℓ` fingerprint provenance — the ℓ=0/1/2 radiative signature `{1, 1/2, 1/27}` at orders `{ω¹, ω³, ω⁵}` + the
  Gate-4 quadrupole non-regression (the ℓ=2 leg reproduces the COMPLETE pathA_33 fold) — **discharges NOTHING** (earned
  exterior-wave structure + a consistency non-regression, not a reduction of a debt; like R37 was for 018's fingerprint).
- **Cited provenance (NOT re-counted):** pathA_33's fold (`{T_Ω,β₂}` counted at 017; the `54/5`/`G`/dim already at 018–021);
  `c_s` (R1, stage005); 008's amplitudes; 009/010's bulk mode. ⚠ **Do NOT let 022 count the ℓ=0/1 stiffnesses `K0c`/
  `K_eta+2·T_Omega` or the return admittances `Z0_ret`/`Z1_ret`** — those are 023's (the transfers/nullspace), NOT 022's.
- **Control/tracked-not-counted:** the sample numeric-substitution dict (numeric cross-check controls, like 014's controls),
  if 022 carries one (likely NOT needed — the fingerprints are exact rationals; confirm at directive).

⚠ **Do NOT let 022 silently count `c_s`, the port scalars, or 023's stiffnesses/admittances.** Resolve `c_s` as
R1-DERIVED-provenance, pathA_33's fold as cited-DONE, 023's return machinery as NOT-022, and Codex-verify (the register
verify is the gate that catches an over-count that would falsely inflate — or a mislabel that would falsely shrink — the
irreducible codimension count).

## Verdict tokens + honest scope

022 carries the **cross-ℓ fingerprint + Gate-4 non-regression component (1/2) of `FAIL_UNDERDETERMINED_NOT_PREDICTIVE` — the
EARNED-FIRST slice**: the outgoing cross-ℓ DtN fingerprints `Ŷ_ℓᵒᵘᵗ=−(ℓ+1)/Λ_ℓ` for ℓ=0,1,2 → radiative `{1, 1/2, 1/27}` at
orders `{ω¹, ω³, ω⁵}` DERIVED from the spherical Hankel `h_ℓ` (not typed), + the Gate-4 quadrupole non-regression (the ℓ=2
leg `−3/Λ₂` reproduces the COMPLETE pathA_33 `QUAD_CALIBRATED` fold `{u₂=1/9, u₄=4/81, v₅=1/27, χ_Q=+1}`). EARNED = the
cross-ℓ exterior radiative signature + the non-regression are genuinely earned (the report's `Earned:` list `:47`); the
native nullspace underdetermination departure (dim-8 / return-nullity-2, the missing Gate-6 `Z0_ret/Z1_ret` selector, the
residuals vs pathA_29, the FAIL_TAUTOLOGICAL firewall) that DELIVERS the FAIL = 023 (the `Deferred:` list). ⭐ **022's SCRIPT
passes (exit 0) — the earned content is correct + every tooth fires; the FAIL-token is the physics LANDING label of the
joint gate, carried as a 1/2 PARTIAL.** Self-contained exterior spherical-Hankel algebra — the Gate-4 non-regression is a
CHECKABLE derive-vs-typed relation against pathA_33's earned values (018∧019∧020∧021, cited DONE); 008's amplitudes +
009/010's bulk mode + `c_s` (R1) cited as PROVENANCE. Caveats: `c_s` is a units carrier not a live value (the fingerprints
are c_s-free); the prefactor (019) / `54/5` (020) / μ̂₀-free dim closure (021) are DONE and NOT re-litigated (022's
non-regression is FINGERPRINT level); the scalar/dipole return MAGNITUDE + the nonzero prediction are DEFERRED (023 — the
native nullspace leaves `ε_eff` free at the linear Gate-5 level); Gate 6 (the `Z0_ret/Z1_ret` selector) stays sim-deferred.
⭐ The pathA_34 trip-ups (023's): the DEFAULT verdict is PASS (022 owns NO nullspace, so it lands the EARNED partial + DEFERS
the underdetermination to 023 — it must NOT compute/print a FAIL from a rigged mechanism); back-solving `ε_eff`/`Z` is
forbidden (022's non-regression is derive-vs-typed, NOT a back-solve); the `−(ℓ+1)/Λ_ℓ` shapes stay COMPUTED (not stamped).

## Process (unchanged, calibrated — the per-stage pipeline)

Author the II-G5a reshape directive (§1 the clean 022 slice / 2-way cut + the STRIP-019-prefactor + the exterior-wave-only
result + §1b the native-`.wl` KEEP + §1c the pathA_33-non-regression-is-checkable / 008-009-010-c_s-provenance framing + §2
faithful cover + §3 bridge-strip incl. sever-YAML + the shared-machinery-duplication + the 022-LOCAL verdict + transcript-
level agreement + §5 the fingerprint-derivation / non-regression `3e` / χ_Q-sign / order-scaling teeth with per-tooth
ablation + §6 the ZERO-new-knobs + R41-edge register question) → **Codex xhigh design-review** → fold to `DIRECTIVE_CLEAN` →
**⭐ final Grok-4.5 headless compute-verify pass** (assess + independently verify each catch — Grok compute-verifies the
cross-ℓ spherical-Hankel series `{1, 1/2, 1/27}` at `{ω¹, ω³, ω⁵}`, the `−(ℓ+1)` ℓ-dependence, the incoming sign flip, the
ℓ=2 non-regression; it caught the 016 volume-vs-line convention mismatch + the 020 rule-inversion, so watch the
non-regression genuineness [derive-vs-typed, not X≡X] + the 022-LOCAL-vs-joint-verdict distinction + the fingerprint-value
computed-ness) → fold → **Codex confirm-pass on the folds** → **pre-exec USER GATE** → Codex builds the two scripts
(`--sandbox danger-full-access`, background, xhigh) → dual-engine both exit 0 (repo root + foreign CWD) → arbiter re-run →
full tri-review (fidelity + adversarial-with-**per-tooth ablation**; ⭐ hunt the fingerprint-value stamped-vs-derived
genuineness + the non-regression derive-vs-typed [not a subsumed X≡X reconstruction] + the χ_Q computed-ness + a
mirror-`.wl` + any vacuous able-to-fail + confirm the 022-LOCAL verdict is NOT the joint `base_verdict`) → remediate →
fresh-agent re-verify → bump counts 21→22 → parameter register (⭐ confirm ZERO new 022 knobs + `c_s` R1-provenance + R41
edge + 023's stiffnesses/admittances NOT counted) + Codex-verify → note/card/`\input{stages/stage_022}` + registration →
rebuild PDF → commit + docs/memory sync (keep STATUS ▶ RESUME HERE THIN; append per-stage detail to part2 Progress).
Orchestrator authors notes/cards/LaTeX/registration; Codex codes. Target stem: `ledger_stage022_cross_l_fingerprints`
(confirm slug at authoring).
