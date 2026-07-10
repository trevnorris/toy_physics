# II-G4d (ledger_stage021) source map — pathA_33 μ̂₀-free `[P₀^phys]=1` dimensional closure (QUAD_CALIBRATED 4/4, COMPLETING leg)

> Running-start prep captured 2026-07-09 (post stage020, before authoring the II-G4d reshape directive) so the directive
> can be written without re-discovery. **All line refs below VERIFIED against the current sources 2026-07-09**
> (`pathA_33_quadrupole_normalization_sympy.py` = 1351 lines; `pathA_33_quadrupole_normalization.wl` = 256 lines;
> `reports/pathA_33_quadrupole_normalization.md` = 49 lines; `directives/pathA_33_quadrupole_normalization.md` = 534 lines).
> Companion: `part2_gravity_atomic_split.md` (row 021 = the pathA_33 4-way `QUAD_CALIBRATED` split, L41; the ▶ NEXT
> stage021 entry L427–439; the pathA_33 trip-ups L87–89; the cross-stage flows L106–108) and the **stage018 + stage019 +
> stage020 source maps + reshape directives** (the pathA_33 4-way-split exemplars — 021 is the FOURTH, COMPLETING leg, the
> one that LANDS the joint `QUAD_CALIBRATED` COMPLETE; ⭐ it reuses 018/019's **KEEP-native-`.wl`-sever-only-YAML** pattern —
> the source `.wl` HAS the 021 dim block — NOT 020's genuinely-authored-`.wl`). Build-order id **021**, Part II. Source
> top-line: **`QUAD_CALIBRATED`** — 021 lands the 4/4 component (the μ̂₀-free `[P₀^phys]=1` dimensional closure); 018
> (fingerprint, DONE `4872e8b7`) is 1/4; 019 (prefactor algebra, DONE `f1c426f9`) is 2/4; 020 (`54/5=2·27/5` provenance
> partition + CALIBRATED label, DONE `4b1247e4`) is 3/4. **021 COMPLETES the pathA_33 fold (018∧019∧020∧021).**

> **⚠ CORRECTION BANNER (2026-07-09, post Codex directive design-review — the directive `_scratch/ledger_stage021_reshape_directive.md`
> SUPERSEDES this map on these SIX points; folded there, pinned here per falsification-first):**
> 1. **The anti-v1 discriminator logic below was WRONG.** A μ̂₀-BACK-SOLVED gate re-solves `[μ̂₀]=(rhs−p0)/2` after ANY
>    corruption, so its `homogeneity_pass` stays True under {`[N₀]`,`[D₀]`,`[G]`,`[c_s]`} — it fires on NOTHING (a tautology).
>    So corrupt-`[G]` does NOT "fire a back-solved gate". The REAL discriminator: the μ̂₀-free gate `[P₀^phys]==ZERO` FIRES on
>    `[N₀]/[D₀]/[c_s]` (able-to-fail) while the back-solve fires on none; the `[G]` NO_FAIL row is only a dependency-scope
>    diagnostic (`G ∉ P₀^phys`). The decisive anti-v1 tooth = a read-set assertion (verdict excludes `mu_hat0`/`mu_dim`/
>    `homogeneity_pass`) + a wired back-solve mutant shown to stay NO_FAIL under all corruptions.
> 2. **Corrupt-`[N₀]` dim:** `[P₀^phys]=(1,−1,0)=L M⁻¹` (the `(c_s/a)²` factor REMAINS), NOT `−[D₀]=(1,−1,2)` (that's `[P₀_raw]`).
> 3. **`[D₀]` provenance:** the carried reduced static conservative denominator `D₀=K−B₀−Z₀` (pathA_43 `calibrated_anchor` /
>    pathA_34 carried), NOT a "continuity-balance denominator". `[N₀]` = the density/continuity-port numerator.
> 4. **Self-containment:** `build_dimensions(fingerprint)` + the `.wl` `YhatPhysical` reference 018's `u2/u4/v5` → supply them as
>    a provenance-frozen LOCAL FIXTURE; REMOVE the 018/019 blocks; "KEEP-native" = the dim ENGINE (L101–174), not the whole file.
> 5. **`.wl` headline:** do NOT include `homogeneityPass` in the local pass guard (it re-wires the μ̂₀ diagnostic into the
>    verdict); split out `yhatOk=(yhatDim==zeroDim)` as its own assert.
> 6. **`Yhat` tooth:** wrap `dim_of` in a structured `try_dim_of`/`Catch` so the named assert CONSUMES the caught mismatch (an
>    uncaught `DimError` mid-eval bypasses the named tooth).

## ⭐ The FIVE headline differences from 020 (READ FIRST)

1. **⚠ 021 is the FOURTH, COMPLETING leg of the SAME 4-way pathA_33 split (018/019/020/021) — the leg that LANDS the JOINT
   `QUAD_CALIBRATED` COMPLETE, NOT a PARTIAL.** Unlike 018/019/020 (each landing `QUAD_CALIBRATED` as a PARTIAL component),
   021 is the COMPLETING leg (the pattern of stage010/012/015/017): its landing prints the joint `QUAD_CALIBRATED`
   **COMPLETE** (all four legs earned), NOT a PARTIAL. ⭐ **BUT "COMPLETE" ≠ "PASS":** the joint verdict token STAYS
   `QUAD_CALIBRATED` (calibrated, not `QUAD_PASS`) — 020 already DETERMINED the CALIBRATED classification (the assembled
   `54/5` + `G` are `external_bridge_input`, `G=GENUINE_BLOCKED`); 021 does NOT flip it to PASS. 021 supplies the
   **dimensional closure** that makes the whole assembled magnitude dimensionally sound (`[P₀^phys]=1`) and completes the
   4-leg fold. So the landing = "the μ̂₀-free dimensional closure that makes the joint `QUAD_CALIBRATED` dimensionally
   consistent AND COMPLETE — all four Gate-4 legs earned." 021 owns ONE earned gate cluster: `build_dimensions`
   (the dim engine + the μ̂₀-free `[P₀^phys]=1` gate + the μ̂₀ diagnostic + the drop-normalization + corrupt-`[N₀]`
   mutations) + probes `3d`/`3d′`. **The 4-way cut: fingerprint (018, EARNED) vs prefactor algebra (019, EARNED) vs the
   `54/5` provenance partition + the CALIBRATED label (020) vs the μ̂₀-free dimensional closure (021).**

2. **⭐⭐ 021 IS UNITS-BEARING AND DIMENSION-CHECKING (the KEY contrast with 020).** 020 is units-bearing but does
   ALGEBRA+PROVENANCE, explicitly NOT the dimensional-homogeneity gate. **021 is the OTHER half of that operation-level
   cut: it does the DIMENSIONAL homogeneity gate `[P₀^phys]=1` and NOTHING of 020's algebra/provenance.** ⚠ The 020/021
   boundary is **operation-level, not expression-level**: both stages reference `P0_physical=(c_s/a)²·(N0/D0)`, `Gamma5`,
   `target_rhs` — 020 uses their VALUES in the algebraic bridge, **021 checks their DIMENSIONS via a real `[·]`
   dimensional-vector algebra** (`dim_of` over an `(L,M,T)` triple). Draw the cut by OPERATION (algebra+provenance = 020;
   dimension-homogeneity = 021), NOT by which symbols appear. 021's genuine tool is the `dim_of` dimensional-vector engine
   — the thing 020 was PROHIBITED from carrying (020 has a structural no-`dim_of` cut). **021 is the ONLY pathA_33 leg with
   a dimensional gate.**

3. **⭐⭐ THE CORE DELIVERABLE = the μ̂₀-FREE `[P₀^phys]=1` dimensional gate that CATCHES THE NATURAL-UNITS TRAP, and the
   corrupt-`[N₀]` vs corrupt-`[G]` DISCRIMINATOR that proves the gate reads the SOURCED port dims (NOT a back-solved μ̂₀).**
   This is the **v1 REJECTION locus** — the leg that was tri-review-REJECTED and remediated. Three joined earned things:
   - **(a) The μ̂₀-FREE gate.** The verdict-bearing check is `dimensional_ok = ([P₀^phys] == ZERO_DIM)` where
     `P₀^phys = (c_s/a)²·(N₀/D₀)` (`.py` L401, L413–414). Its dimension is read from the **SOURCED port dims**
     `[N₀]=L⁻¹M` (`SOURCED_N0_DIM=(-1,1,0)`, L379) and `[D₀]=L⁻¹T⁻²M` (`SOURCED_D0_DIM=(-1,1,-2)`, L380) plus
     `[(c_s/a)²]=T⁻²`. So `[P₀_raw]=[N₀/D₀]=T²`, `[(c_s/a)²]=T⁻²`, `[P₀^phys]=1` (report :22). ⚠ **`μ̂₀` (`mu_hat0`) does
     NOT enter the verdict** — the v1 tautology BACK-SOLVED the free carrier `μ̂₀` from `(rhs_dim − p0_dim)/2` (L428) to
     FORCE homogeneity, which can never fail. In the reshaped gate `μ̂₀` is confined to a labelled DIAGNOSTIC block
     (`mu_hat0_homogeneity_diagnostic`, L476–486) with `participates_in_verdict: False` (L478) and
     `label: "non-able-to-fail (mu_hat0 free carrier)"` (L477). **The gate must NOT become the v1 μ̂₀-back-solved
     homogeneity form** (`part2` L88, the pathA_33 trip-up 021 OWNS).
   - **(b) The natural-units-trap catch (the DIMENSIONAL MILESTONE).** The handoff's `P₀=N₀/D₀` silently DROPS a
     `(c_s/a)²` frequency-normalization factor (a natural-units artifact where `c_s/a` was set to 1). The gate CATCHES it:
     the `mutation_drop_frequency_normalization` (L487–495) checks `[P₀_raw]=[N₀/D₀]=T²≠1`, so the `3d` drop-normalization
     probe (L820–826) fires `FAIL_DIMENSIONAL`. This is the earned finding: `[P₀_raw]=T²`, `[(c_s/a)²]=T⁻²`,
     `[P₀^phys]=1` — the frequency normalization is REQUIRED for dimensional consistency (report :22, :24).
   - **(c) The corrupt-`[N₀]` vs corrupt-`[G]` DISCRIMINATOR (the natural-units-trap discriminator; the v1 remediation —
     KEEP genuine + able-to-fail).** The `mutation_corrupt_N0_dimension` (L496–506) sets `[N₀]→ZERO_DIM`, then
     `[P₀^phys]` becomes `−[D₀]=(1,−1,2)≠1`, so the `3d′` corrupt-port-dimension probe (L828–836) fires `FAIL_DIMENSIONAL`
     (report :25). ⭐ **The DISCRIMINATOR (part2 trip-up L88, STATUS "corrupting `[N₀]` now fires `FAIL_DIMENSIONAL`;
     `[G]` doesn't"):** corrupting `[N₀]` FIRES the gate (because `N₀ ∈ free_symbols(P₀^phys)`), but corrupting `[G]` does
     **NOT** fire it (because `G ∉ free_symbols(P₀^phys)` — `P₀^phys=(c_s/a)²·N₀/D₀` has no `G`). This PROVES the gate
     reads the SOURCED port dims `[N₀]/[D₀]`, NOT a μ̂₀-back-solve (which WOULD be `[G]`-sensitive, since the back-solve
     goes through `target_rhs=54Gc_s⁵/5a⁵c⁵`, L409, → `rhs_dim`, L427, → `mu_dim`, L428). ⚠ **The source has ONLY the
     corrupt-`[N₀]` mutation (which fires) — it does NOT ship an explicit corrupt-`[G]`-does-NOT-fire NEGATIVE control.
     021 MUST ADD it** (a computed runtime check: corrupt `[G]` → `dimensional_ok` UNCHANGED (True) → `NO_FAIL`, while the
     μ̂₀ DIAGNOSTIC's `homogeneity_pass` may change) — this is the second half of the natural-units-trap discriminator and
     the sharpest proof the gate is μ̂₀-free. See §5 tooth 3.

4. **⭐ The consumption of 018 (fingerprint `u₂/u₄/v₅`) + 019 (`P0=N0/D0`) + the SOURCED port dims — GENUINELY ENTERS
   021's SELF-CONTAINED dim check, but stays PROVENANCE-cited (no theatrical cross-stage dual-site).** Same discipline as
   018/019/020. 021's dim gate genuinely DEPENDS on three cited inputs, all used as PROVENANCE:
   - **The SOURCED port dims `[N₀]=L⁻¹M`, `[D₀]=L⁻¹T⁻²M`** (`.py` L379–380) — the port-kernel dimensional data the closure
     ASSUMES. These are the DENSITY-PORT dims: `[N₀]=L⁻¹M` is pathA_43's earned density-port numerator dimension
     (`DENSITY_PORT_HOSTED`, Cluster C stages 024/027, NOT yet built in build order), and `[D₀]=L⁻¹T⁻²M` is the
     continuity-balance denominator dim. ⚠ They are NOT literally 017's exported ℓ=2 D-lane dim (`[D0]=M T⁻²`, VOLUME-density
     pathA_32 convention) — a DIFFERENT convention. So cite them as **SOURCED port-dimension inputs** (pathA_43 density-port
     provenance), NOT a checkable cross-stage relation to 017's D-lanes. The genuine earned content is that GIVEN those
     sourced dims + `(c_s/a)²`, `[P₀^phys]=1`, and the natural-units trap is caught.
   - **018's fingerprint `u₂=a²/9c_s²`, `u₄=4a⁴/81c_s⁴`, `v₅=a⁵/27c_s⁵`** (via `fingerprint["outgoing"]["u2"/"u4"/"v5"]`,
     L402–406) — enter the `Yhat_out_physical_expansion` dimensionless check (L441): `1 + u₂ω² + u₄ω⁴ + i·v₅ω⁵` is
     dimensionless (each term: `[u₂ω²]=[a²/c_s²]·[ω²]=T²·T⁻²=1`). A genuine dim-consistency check tying 018's fingerprint to
     dimensionlessness. Cite `stage018` (the values are 018's DERIVED fingerprint; 021 dim-checks them, does NOT re-derive).
   - **019's `P0=N0/D0`** (the prefactor structure) — enters `P0_raw=N0/D0` (L399) → `P0_physical` (L401). Cite `stage019`
     (definitional/provenance; the abstract `N₀/D₀` structure).
   - **⭐ RESOLUTION (the honest framing — same discipline as 018/019/020):** cite 018's fingerprint / 019's `P0=N0/D0` /
     the SOURCED port dims as **PROVENANCE**; the `[P₀^phys]=1` gate + the natural-units-trap catch + the corrupt-`[N₀]`/`[G]`
     discriminator are 021's **SELF-CONTAINED earned checks**. ⚠ **Per the emitted-vs-checked test:** UNLIKE 019's
     `build_port_moments` (emitted-but-never-checked → a guard on them = a vacuous tooth), the SOURCED port dims
     `[N₀]/[D₀]` GENUINELY ENTER 021's `[P₀^phys]` gate (they ARE checked) — so the corrupt-`[N₀]` probe IS the genuine
     able-to-fail tooth (it reads the sourced dim), NOT a vacuous guard. There is still **NO theatrical cross-stage
     dual-site** (a "second site" re-reading 017's/pathA_43's port dim would be subsumed); the genuine content is 021's own
     dim gate, and 018/019/pathA_43 supply the provenance of WHAT the fingerprint / `P0` / port dims are.

5. **⭐⭐ The `.wl` must KEEP its NATIVE 021 dimensional block (L101–174) — the source `.wl` HAS the 021 content (the KEY
   contrast with 020).** ⚠ **VERIFIED (2026-07-09 full `.wl` read):** the source `pathA_33_quadrupole_normalization.wl`
   (256 lines) has a native, already-independent **dimensional block L101–174**: `zeroDim`/`dimAdd`/`dimScale`/`dimOf`
   (L101–124, a native Wolfram `Which`-based dimensional-vector engine — a genuinely different construction from the `.py`'s
   `dim_of`), `expText`/`dimText` (L125–133), `rawDims` (L135–145, native `<|…|>` Association with `N0->{-1,1,0}`,
   `D0->{-1,1,-2}`), `P0Raw`/`frequencyNormalization`/`P0Physical` (L146–148), `YhatPhysical` (L149), `Gamma5`/`targetRHS`
   (L150–151), `p0RawDim`/`frequencyNormDim`/`p0Dim`/`dimensionalOk` (L152–155, the μ̂₀-free gate), `dropNormDim`/`dropNormOk`/
   `dropNormVerdict` (L156–158, the drop-normalization mutation), `corruptN0Dims`/`corruptN0RawDim`/`corruptN0P0Dim`/
   `corruptN0Ok`/`corruptN0Verdict` (L159–163, the corrupt-`[N₀]` mutation), `rhsDim`/`muDim`/`requiredP0Dim` (L164–171, the
   μ̂₀ DIAGNOSTIC back-solve), `gamma5Dim`/`yhatDim`/`homogeneityPass` (L172–174). So — LIKE 018/019 (which KEPT the `.wl`'s
   already-native fingerprint/prefactor blocks + severed only the YAML) and UNLIKE 020 (which genuinely AUTHORED its `.wl`) —
   **021 KEEPS the native `.wl` dimensional block L101–174 as the independent Wolfram engine**; sever ONLY the scratch-YAML
   machinery (`scratchDir`/`yamlOut` L20–22 + `CreateDirectory` L22 + `Export[yamlOut,…]` L253 + the YAML-line assembly
   `lines={…}` L182–251 + the `headlineOk` guard L254). ⭐ **ADD to the `.wl`:** the corrupt-`[G]`-does-NOT-fire negative
   control (§5 tooth 3) native — the source `.wl` lacks it too (it has `corruptN0…` but no `corruptG…`).

## §1 The 021 slice (`.py` line ranges) — the CLEAN CUTS (all VERIFIED)

The whole computation is `build_*` helpers assembled by `build_context` (L972–1017) + `build_counterfactuals` (L703–902).
**021 owns the dim engine + `build_dimensions` + probes `3d`/`3d′` + the dim-scoped self-ablation; 018/019/020 own the
fingerprint / prefactor / scaling+equivalence+provenance blocks.** The 021-owned cuts:

- **The dimensional-vector engine `DimError` L297–298 + `dim_add`/`dim_sub`/`dim_scale` L301–310 + `dim_of` L313–340 +
  `dim_to_monomial` L343–344 + `exp_text` L347–351 + `dim_to_text` L354–365 + `dim_record` L368–376.** This is 021's genuine
  tool: a recursive `(L,M,T)`-triple dimension reader over SymPy exprs (Symbol/Mul/Pow/Add), raising `DimError` on a missing
  symbol dim, a non-numeric exponent, or a **sum-dimension mismatch** (L336–338 — the homogeneity teeth). ⚠ **Internal
  triple order is `(L, M, T)`** (L55–57, L302; `dim_to_text` L354–365 re-orders to the human `L,T,M` display). Keep this
  faithfully — the SOURCED dims below are in `(L,M,T)`.
- **`SOURCED_N0_DIM = (-1, 1, 0)` L379** (`[N₀]=L⁻¹M`, the density-port numerator) + **`SOURCED_D0_DIM = (-1, 1, -2)` L380**
  (`[D₀]=L⁻¹T⁻²M`, the continuity-balance denominator) + `dim_vector_text` L383–384.
- **⭐ `build_dimensions(fingerprint)` L387–507 — the 021 CORE:**
  - `raw_symbol_dims` L388–398: `a=(1,0,0)`, `c_s=(1,0,-1)`, `c=(1,0,-1)`, `G=(3,-1,-2)`, `omega=(0,0,-1)`, `D0=SOURCED_D0_DIM`,
    `N0=SOURCED_N0_DIM`, `chi_Q=ZERO_DIM`, `mtilde0=ZERO_DIM`. (⚠ `mu_hat0` is NOT in `raw_symbol_dims` — it is added ONLY to
    the diagnostic `symbol_dims` at L429–430, AFTER back-solving `mu_dim`, so the gate can NEVER read it.)
  - `P0_raw = N0/D0` L399 ; `frequency_normalization = (c_s/a)²` L400 ; `P0_physical = frequency_normalization·P0_raw` L401.
  - `yhat_physical = 1 + u2·ω² + u4·ω⁴ + i·v5·ω⁵` L402–407 (018's fingerprint values).
  - `Gamma5 = chi_Q·P0_physical·a⁵/(27·c_s⁵)` L408 ; `target_rhs = 54·G·c_s⁵/(5·a⁵·c⁵)` L409 (re-defined here for the
    DIMENSION check only — NOT 020's algebra).
  - `raw_p0_dim = dim_of(P0_raw)` L411 (`=T²`) ; `frequency_norm_dim = dim_of(frequency_normalization)` L412 (`=T⁻²`) ;
    **`p0_dim = dim_of(P0_physical)` L413 (`=ZERO_DIM`) ; `dimensional_ok = (p0_dim == ZERO_DIM)` L414 — THE μ̂₀-FREE GATE.**
  - **`drop_norm` mutation L416–418:** `drop_norm_dim = dim_of(P0_raw)` (`=T²≠ZERO`) → `drop_norm_ok=False` →
    `drop_norm_verdict="FAIL_DIMENSIONAL"`.
  - **`corrupt_n0` mutation L420–425:** `corrupt_n0_dims[N0]=ZERO_DIM` → `corrupt_n0_p0_dim = dim_of(P0_physical)`
    (`=−[D₀]=(1,−1,2)≠ZERO`) → `corrupt_n0_ok=False` → `corrupt_n0_verdict="FAIL_DIMENSIONAL"`.
  - **The μ̂₀ DIAGNOSTIC (non-verdict) L427–448:** `rhs_dim = dim_of(target_rhs)` L427 ;
    `mu_dim = (rhs_dim − p0_dim)/2` L428 (the BACK-SOLVE — confined to the diagnostic) ; `symbol_dims[mu_hat0]=mu_dim`
    L429–430 ; `lhs = (mu_hat0·mtilde0)²·P0_physical` L431 ; `lhs_raw_mutation = (mu_hat0·mtilde0)²·P0_raw` L432 ;
    `required_p0_dim = rhs_dim − 2·mu_dim` L433 ; the `table` L435–445 ; `lhs_dim`/`lhs_raw_dim` L446–447 ;
    `homogeneity_pass = (lhs_dim==rhs_dim ∧ p0_dim==required_p0_dim)` L448.
  - The returned dict L449–507: `dimensional_gate` / `dimensional_gate_depends_on_mu_hat0: False` L451–453 ; the
    `symbol_dimensions` table L454–463 ; `P0_raw_dimension`/`frequency_normalization_dimension`/`P0_dimension`/
    `P0_physical_dimension` L464–468 ; `dimensional_ok`/`dimensional_status` L469–470 ; `mu_hat0_dimension`/`Gamma5_dimension`/
    `lhs_dimension`/`rhs_dimension` L471–474 ; the `table` L475 ; **`mu_hat0_homogeneity_diagnostic` L476–486**
    (`label:"non-able-to-fail (mu_hat0 free carrier)"`, `participates_in_verdict:False`, `homogeneity_pass`) ;
    **`mutation_drop_frequency_normalization` L487–495** (`verdict`, `mutation_fires`) ; **`mutation_corrupt_N0_dimension`
    L496–506** (`verdict`, `mutation_fires`).
- **The 021 probes (in `build_counterfactuals` L703–902):**
  - **probe `3d_dimensional_break` L820–826** — reads `dimensions["mutation_drop_frequency_normalization"]`,
    `mutated_dimensional_ok`, `verdict` (`FAIL_DIMENSIONAL`); `self_ablation = dimensional_ablation` (L760, DYNAMIC).
  - **probe `3d_prime_corrupt_port_dimension` L828–836** — reads `dimensions["mutation_corrupt_N0_dimension"]`,
    `sourced_N0_dimension`/`corrupted_N0_dimension`, `mutated_P0_physical_dimension`, `verdict` (`FAIL_DIMENSIONAL`);
    `self_ablation = dimensional_ablation` (L760, DYNAMIC).
- **The dim-scoped self-ablation `dimensional_ablation` L760** — `ablation({"dimensional_ok": False}, expected_fail=
  "FAIL_DIMENSIONAL")` via the DYNAMIC `ablation` helper L728–755 (`base_verdict(mutated_gates)` vs `base_verdict(baseline)`,
  `with_mutation`≠`without_mutation`, `rerun_gate_logic:True`). ⚠ **SCOPE caveat (SAME as 018/019/020):** `base_verdict`
  (L676–694) re-runs the JOINT gate set incl. 018/019/020's gates (`fingerprint_ok`/`prefactor_ok`/`scaling_ok`/
  `equivalence_ok`/`outgoing_ok`/`provenance_ok`), NOT built here. **Re-scope the `3d`/`3d′` `with_mutation`/`without_mutation`
  to an 021-LOCAL verdict over 021's gate only** (`dimensional_ok` — the μ̂₀-free `[P₀^phys]=1` gate, plus the drop-norm +
  corrupt-`[N₀]` mutation-fires being able-to-fail), do NOT pull the joint `base_verdict`. It MUST be a real two-verdict
  re-run (`with_mutation`=`FAIL_DIMENSIONAL` ≠ `without_mutation`=`DIMENSIONAL_OK`), NOT a constant (the v1 rig trip-up,
  `part2` L88).
- **Shared helpers 021 uses (NOT cut boundaries):** `compact` L60, `hstr` L66, symbol decls L45–57 (021 uses `a, c_s, c, G,
  omega, D0, N0, chi_Q, mtilde0` + the dimension symbols `Ldim/Mdim/Tdim` L55 + `Dim`/`ZERO_DIM` L56–57 + `mu_hat0` L48 for
  the DIAGNOSTIC only). 021 does NOT need `z, j₂/y₂` (018's), `D2/D4/N2/N4` (019's prefactor), `lambda_G` (020's g-invariance),
  the fingerprint's `chi_Q_incoming` (018's) — it CITES 018's `u₂/u₄/v₅` (via the `fingerprint` arg) as provenance.

- **⭐ CLEAN CUT — 021 owns the dim engine `DimError`/`dim_*`/`dim_to_*`/`dim_record` L297–376 + `SOURCED_N0_DIM`/
  `SOURCED_D0_DIM`/`dim_vector_text` L379–384 + `build_dimensions` L387–507 + probes `3d` L820–826 / `3d′` L828–836 + the
  dim-scoped self-ablation L760. It touches NONE of `build_fingerprint`/`build_prefactor`/`build_scaling`/`build_equivalence`/
  the provenance machinery. Do NOT pull in 018/019/020 territory:**
  - **018 (DONE `4872e8b7`):** `spherical_j2`/`spherical_y2` L101–106 + `dtn_branch` L109–149 + `build_fingerprint` L152–187
    + `passivity_from_source` L666–673 + probes `3a` L774–797 / `3b` L799–812. (021 CITES 018's `u₂/u₄/v₅` fingerprint values
    via the `fingerprint` arg for the `Yhat` dimensionless check — §1c.)
  - **019 (DONE `f1c426f9`):** `build_port_moments` L190–209 + `build_prefactor` L212–273 + probe `3g` L866–875 +
    `prefactor_ablation` L762. (021 CITES 019's `P0=N0/D0` structure as provenance — it enters `P0_raw` — §1c.)
  - **020 (DONE `4b1247e4`):** `a_power`/`build_scaling` L276–294 + `build_equivalence` L510–535 + the provenance machinery
    L538–663 + probes `3c` L813–818 / `3e` L837–841 / `3f` L846–865 + `partition_ablation` L764–771 + the CALIBRATED verdict
    determination L690–694. ⚠ **`P0_physical`/`Gamma5`/`target_rhs` are RE-DEFINED in BOTH `build_equivalence` (020, algebra)
    AND `build_dimensions` (021, dimension) — 020 uses their VALUES in the bridge, 021 checks their DIMENSIONS. The cut is by
    OPERATION** (headline-2). Keep 021 free of ANY `compact(expr)==0` algebraic-residual bridge, any `classify_provenance`/
    partition, any `a_power`/scaling, any `/. G->lambdaG G` g-invariance — 021 does dimensions ONLY.

## §1b The `.wl` 021 slice (VERIFIED) — KEEP-native (sever only YAML), NOT genuine-authoring like 020

⭐ **VERIFIED (full `.wl` read 2026-07-09): the source `.wl` HAS the 021 dimensional block L101–174** (headline-5) — a
native, already-independent Wolfram dimensional-vector engine (`Which`-based `dimOf`, native `<|…|>` `rawDims` Association,
native `Series`-free dimension algebra). So — LIKE 018/019, UNLIKE 020 — **021 KEEPS the native `.wl` dim block as the
independent engine** and severs only the YAML machinery.

**Keep (native, already-independent):** the whole dimensional block L101–174 —
- `zeroDim`/`dimAdd`/`dimScale`/`dimOf` L101–124 (native `Which` over `Times`/`Power`/`Plus`, `fail[]` on missing
  dim/non-numeric exponent/sum mismatch — the homogeneity teeth); `expText`/`dimText` L125–133.
- `rawDims` L135–145 (`N0->{-1,1,0}`, `D0->{-1,1,-2}`, `a->{1,0,0}`, `cs->{1,0,-1}`, `c->{1,0,-1}`, `G->{3,-1,-2}`,
  `omega->{0,0,-1}`, `chiQsym->zeroDim`, `mtilde0->zeroDim`; ⚠ `muHat0` NOT in `rawDims` — added only at L166 for the
  diagnostic).
- `P0Raw`/`frequencyNormalization`/`P0Physical` L146–148 ; `YhatPhysical` L149 ; `Gamma5`/`targetRHS` L150–151.
- `p0RawDim`/`frequencyNormDim`/`p0Dim`/`dimensionalOk` L152–155 (the μ̂₀-free gate) ; `dropNormDim`/`dropNormOk`/
  `dropNormVerdict` L156–158 (drop-normalization mutation) ; `corruptN0Dims`/`corruptN0RawDim`/`corruptN0P0Dim`/
  `corruptN0Ok`/`corruptN0Verdict` L159–163 (corrupt-`[N₀]` mutation).
- `rhsDim`/`muDim`/`dims`/`lhs`/`lhsRawMutation`/`lhsDim`/`lhsRawDim`/`requiredP0Dim`/`gamma5Dim`/`yhatDim`/`homogeneityPass`
  L164–174 (the μ̂₀ back-solve DIAGNOSTIC + the `yhatDim` dimensionless check).

**Sever (the YAML bridge):** `scratchDir`/`yamlOut` L20–22 (+ the `CreateDirectory` L22) + the YAML-line assembly
`lines={…}` L182–251 + `Export[yamlOut,…]` L253 + the `headlineOk` guard L254 (L176–180 `headlineOk` currently AND-s
`u2Match…p4Match` (018/019's) with `dimensionalOk`, `!dropNormOk`, `!corruptN0Ok` — **re-scope the reshaped `.wl`'s pass
guard to 021's dim gates ONLY**: `dimensionalOk ∧ !dropNormOk ∧ !corruptN0Ok ∧ homogeneityPass ∧ (the ADDED corrupt-`[G]`
NO-fire control)`; drop the 018/019 `u*/p*Match` from 021's guard).

**Re-shape to ledger idioms:** print-only + `fail[]`/`Exit[1]` on failure (the `.wl` already has `fail[msg_]` L5). ⭐ **ADD
native:** (a) the corrupt-`[G]`-does-NOT-fire negative control (`corruptGDims = Join[KeyDrop[rawDims,G], <|G->zeroDim|>]`,
`dimOf[P0Physical, corruptGDims] == zeroDim` STILL True → `NO_FAIL`; because `G ∉ P0Physical`); (b) the DYNAMIC 021-local
self-ablation for `3d`/`3d′` (a two-verdict re-run over an 021-local verdict, mirroring the `.py`). **Dual-engine agreement =
transcript-level** (both engines print `[P₀_raw]=T²`, `[(c_s/a)²]=T⁻²`, `[P₀^phys]=1`, `dimensionalOk=True`, `3d`
`FAIL_DIMENSIONAL`, `3d′` `FAIL_DIMENSIONAL`, corrupt-`[G]` `NO_FAIL`, `homogeneityPass=True`). **Zero file I/O.** Arity
discipline (standing — def/call scan + unevaluated-leakage transcript scan; the stage007 lesson — the `.wl` has `Module`s
and the `dimOf` recursion, scan them).

## §1c The consumption resolution (READ — the honest PROVENANCE framing; the KEY 021 discipline)

⭐ **021's consumption of 018 (fingerprint `u₂/u₄/v₅`) + 019 (`P0=N0/D0`) + the SOURCED port dims — cite as PROVENANCE; the
`[P₀^phys]=1` gate + the natural-units-trap catch + the corrupt-`[N₀]`/`[G]` discriminator are 021's SELF-CONTAINED checks.
Do NOT manufacture a theatrical cross-stage dual-site.**

- **The SOURCED port dims `[N₀]=L⁻¹M` / `[D₀]=L⁻¹T⁻²M`** (`.py` L379–380): the port-kernel dimensional inputs the dim
  closure ASSUMES. Cite the DENSITY-PORT provenance — `[N₀]=L⁻¹M` is pathA_43's `DENSITY_PORT_HOSTED` numerator dim (Cluster
  C stages 024/027, not yet built), `[D₀]=L⁻¹T⁻²M` the continuity-balance denominator dim. ⚠ NOT 017's `[D0]=M T⁻²`
  (VOLUME-density convention — DIFFERENT). ⭐ **These genuinely ENTER 021's gate (they are CHECKED, not emitted-but-never-
  checked)** — so the corrupt-`[N₀]` probe IS a genuine able-to-fail tooth (per the emitted-vs-checked test, contrast 019's
  `build_port_moments`). No dual-site: a "second site" re-reading pathA_43's/017's port dim would be subsumed.
- **018's fingerprint `u₂/u₄/v₅`:** enter the `Yhat_out_physical_expansion` dimensionless check (each term dimensionless).
  Cite `stage018` (018's DERIVED fingerprint; 021 dim-checks, does NOT re-derive). Self-contained check.
- **019's `P0=N0/D0`:** enters `P0_raw=N0/D0` → `P0_physical`. Cite `stage019` (the prefactor structure); definitional/
  provenance in 021.
- ⚠ **NO theatrical cross-stage dual-site** on the fingerprint / `P0` / port dims (the 017 lesson `part2` L331–334; the
  018/019/020 §1c discipline). The genuine earned content is 021's self-contained `[P₀^phys]=1` gate + the natural-units-trap
  catch + the corrupt-`[N₀]`/`[G]` discriminator; 018/019/pathA_43 supply the provenance of WHAT the fingerprint / `P0` /
  port dims are. **The sourced values genuinely enter 021's dim algebra (so a corrupt cited `[N₀]` breaks `[P₀^phys]` — that
  IS the `3d′` tooth), but the check is 021's own dim algebra, not a cross-stage read.**

## §2 The 021 claim-set (derive + assert; report/directive quotes)

- **⭐ The μ̂₀-FREE `[P₀^phys]=1` dimensional gate (EARNED — the DIMENSIONAL MILESTONE; report :22, :3, :5).**
  `[P₀^phys]=(c_s/a)²·(N₀/D₀)` is dimensionless from the SOURCED port dims `[N₀]=L⁻¹M`, `[D₀]=L⁻¹T⁻²M`:
  `[P₀_raw]=[N₀/D₀]=T²`, `[(c_s/a)²]=T⁻²`, `[P₀^phys]=1`; `dimensional_ok=True`. ⚠ **μ̂₀ does NOT enter the verdict**
  (`dimensional_gate_depends_on_mu_hat0=False`, L453; `participates_in_verdict=False`, L478). **Do NOT let the gate become
  the v1 μ̂₀-back-solved homogeneity form.**
- **The natural-units-trap catch (EARNED — the milestone; report :22, :24, :37 `3d`).** The handoff's `P₀=N₀/D₀` silently
  drops `(c_s/a)²`; the gate catches it — `[P₀_raw]=T²≠1`, so `3d_dimensional_break` fires `FAIL_DIMENSIONAL` (the frequency
  normalization is REQUIRED for `[P₀^phys]=1`).
- **⭐ The corrupt-`[N₀]` vs corrupt-`[G]` discriminator (EARNED — the v1 remediation; report :25, :38 `3d′`).** Corrupting
  `[N₀]→(0,0,0)` makes `[P₀^phys]=−[D₀]≠1` → `3d_prime` fires `FAIL_DIMENSIONAL` (`N₀ ∈ P₀^phys`). ⭐ **ADD:** corrupting
  `[G]` does NOT fire (`G ∉ P₀^phys`) — the discriminator that the gate reads the SOURCED port dims, NOT a μ̂₀-back-solve
  (which WOULD depend on `[G]` via `target_rhs→rhs_dim→mu_dim`). §5 tooth 3.
- **The μ̂₀ homogeneity DIAGNOSTIC (NON-gate, EARNED-as-diagnostic; report :23).** The full closure `m̂₀²P₀=54Gc_s⁵/5a⁵c⁵`
  back-solves `[μ̂₀]=L⁻¹T⁻¹M⁻¹ᐟ²` from `(rhs_dim−p0_dim)/2`, giving `[lhs]=[rhs]=L⁻²T⁻²M⁻¹`, `homogeneity_pass=True` — a
  DIAGNOSTIC (`participates_in_verdict=False`), NOT the verdict-bearing gate. This is the v1-rejected tautology, correctly
  DEMOTED to a labelled non-gate; it must STAY non-verdict-driving.
- **The `Yhat` dimensionless check (EARNED).** `1 + u₂ω² + u₄ω⁴ + i·v₅ω⁵` is dimensionless (018's fingerprint × ω-powers),
  `[Yhat]=1` (the `.wl` `yhatDim` in `homogeneityPass`).
- **The able-to-fail probes (EARNED).** `3d_dimensional_break` (drop `(c_s/a)²` → `FAIL_DIMENSIONAL`, report :37),
  `3d_prime_corrupt_port_dimension` (corrupt `[N₀]` → `FAIL_DIMENSIONAL`, report :38), + the ADDED corrupt-`[G]`-NO-fire
  control, each with a DYNAMIC 021-local self-ablation.
- **⭐ The 021-scoped landing (the COMPLETING leg — LANDS THE JOINT COMPLETE).** `QUAD_CALIBRATED (4/4) = the μ̂₀-free
  `[P₀^phys]=1` dimensional closure ∧ the natural-units-trap catch ∧ the corrupt-`[N₀]`/`[G]` discriminator`, COMPLETING the
  joint `QUAD_CALIBRATED` (018 fingerprint ∧ 019 prefactor ∧ 020 provenance-partition+CALIBRATED ∧ 021 dim closure).
  ⭐ **021 PRINTS the joint COMPLETE** (all four legs earned) — but the joint TOKEN stays `QUAD_CALIBRATED` (calibrated, not
  PASS; 020's provenance + `G=GENUINE_BLOCKED`). Do NOT re-present 018/019/020's fingerprint/prefactor/provenance accounting;
  CITE them as DONE and declare the 4-leg fold COMPLETE.

## §3 Reshape cost (the bridge to sever) — cross-script scratch-YAML family; KEEP the native `.wl`

Same family as pathA_30–34 / 018 / 019 / 020 (the cross-script runtime-YAML reshape). No argparse. The `.py`'s `build_*`
helpers are pure/self-contained, but `main` L1321–1347 writes `SYM_YAML` (L1324), reads `MMA_YAML` (L1326), writes
`RESULTS_YAML`/`REPORT_MD`/`FEED_NOTE` (L1333–1335); `compare_engines` L1024–1099 cross-checks; the `.wl` writes its scratch
YAML via `Export` L253. **Reshape = sever ALL file I/O both directions:** drop `main`'s YAML/report/feed writers +
`yaml_read`/`yaml_write` (L72–86) + `compare_engines`/`engine_summary`/`build_final_payload`/`build_report`/`build_feed_note`
(`.py`); drop the `Export` + the YAML-line assembly L182–255 (`.wl`). Each engine → standalone: print-only,
`expect_zero`/`bool_from_residual`/`expect_bool`/`expect_fail`-style asserts (`.py` local ledger idioms —
`banner`/`subbanner`/`_record_pass`/`_record_fail`/`expect_zero`/`expect_bool`/`expect_fail`, a `Verdict labels:` block,
tallies, `OVERALL PASS`/nonzero exit), `fail[]`/`Exit[1]` on failure (`.wl`). ⚠ **UNLIKE 020 (genuine-authoring), the 021
`.wl` KEEPS its native dim block L101–174** (§1b): sever ONLY the YAML machinery, ADD the corrupt-`[G]` control + the
DYNAMIC self-ablation native. **Dual-engine agreement = transcript-level** (both engines print the same `[P₀_raw]=T²`,
`[(c_s/a)²]=T⁻²`, `[P₀^phys]=1`, the `3d`/`3d′` `FAIL_DIMENSIONAL`, the corrupt-`[G]` `NO_FAIL`, `homogeneity_pass=True`).
**Zero file I/O.** Arity discipline (standing).

## §4 Consumed / exported

- **Consumes (PROVENANCE only — cite, do NOT re-derive, do NOT build a theatrical dual-site; §1c):**
  - **The SOURCED port dims `[N₀]=L⁻¹M` / `[D₀]=L⁻¹T⁻²M`** — the density-port dimensional inputs (pathA_43 `DENSITY_PORT_HOSTED`
    numerator dim + the continuity-balance denominator dim; Cluster C stages 024/027 not yet built); cite the density-port
    provenance. They GENUINELY enter 021's `[P₀^phys]` gate (CHECKED, not emitted-but-never-checked) — the corrupt-`[N₀]`
    probe is the genuine tooth.
  - **018's fingerprint `u₂/u₄/v₅`** — enter the `Yhat` dimensionless check; cite `stage018`. Self-contained check.
  - **019's `P0=N0/D0`** — enters `P0_raw`; cite `stage019`. Definitional/provenance.
  - **020's assembled magnitude `54Gc_s⁵/5a⁵c⁵`** (`target_rhs`) — enters the μ̂₀ DIAGNOSTIC (`rhs_dim`) ONLY (not the gate);
    cite `stage020`. Provenance (the diagnostic's rhs).
- **Exports (→ 022/027 + Part VII):** the μ̂₀-free `[P₀^phys]=1` dimensional closure + the COMPLETED joint `QUAD_CALIBRATED`
  (all four Gate-4 legs) → **022** (Gate-4 non-regression: pathA_34 cross-ℓ cites the completed quadrupole normalization) +
  **027** (pathA_43 port checks + the `P0_phys=(c_s/a)²N0/D0` closure slot — the same dim closure, marked shared) + the
  Part-VII whole-system dimensional-firewall check (021's per-stage dim closure feeds it). Per the cross-stage flow
  (`part2` L107): "018–021 export the Λ₂ fingerprint + χ_Q=1 + the 54/5 partition → 022 + 027." Cite the exact export contract
  at note-authoring.

## §5 Teeth candidates (021-specific, per-tooth ablation MANDATORY)

1. **⭐⭐ The μ̂₀-free `[P₀^phys]=1` gate + the natural-units-trap catch (`3d_dimensional_break`, `FAIL_DIMENSIONAL`) — the
   CENTRAL 021 tooth.** Drop `(c_s/a)²` → `[P₀_raw]=T²≠1` → fires. Per-tooth ablation: mutate `P0_physical→P0_raw` (drop the
   norm) → the dim gate fires `FAIL_DIMENSIONAL` AT its own assert; the correct `P0_physical` does NOT. ⚠ Assert the gate
   reads the COMPUTED `dim_of(P0_physical)`, NOT a typed `[P₀^phys]="1"` string (a computed runtime guard, not a grep —
   [[feedback-grep-acceptance-dodgeable]]).
2. **⭐ The corrupt-`[N₀]` port-dimension tooth (`3d_prime_corrupt_port_dimension`, `FAIL_DIMENSIONAL`).** Corrupt
   `[N₀]→(0,0,0)` → `[P₀^phys]=−[D₀]≠1` → fires. Per-tooth: mutate the sourced `[N₀]` → fires; the correct sourced dim does
   NOT. This proves the gate reads the SOURCED port dim (not μ̂₀-back-solved).
3. **⭐⭐ ADD: the corrupt-`[G]`-does-NOT-fire NEGATIVE control (the natural-units-trap DISCRIMINATOR — the sharpest μ̂₀-free
   proof).** The source ships corrupt-`[N₀]` (fires) but NO corrupt-`[G]` negative control. ⭐ **021 must ADD it (computed):**
   corrupt `[G]→(0,0,0)` → `dimensional_ok` UNCHANGED (`[P₀^phys]=1` still, because `G ∉ free_symbols(P0_physical)`) →
   `NO_FAIL`; while the μ̂₀ DIAGNOSTIC's `homogeneity_pass` MAY change (the back-solve goes through `target_rhs`, which
   carries `G`). **This is a computed DISCRIMINATOR, not a grep:** assert (a) `G ∉ free_symbols(P0_physical)` (runtime), (b)
   corrupt-`[G]` → the gate's `dimensional_ok` stays True (`NO_FAIL`), (c) corrupt-`[N₀]` → the gate's `dimensional_ok` flips
   False (`FAIL_DIMENSIONAL`). The CONTRAST (N₀ fires / G doesn't) proves the gate is μ̂₀-free (a back-solved gate would be
   `[G]`-sensitive). Per-tooth ablation on the control: if the gate were the v1 back-solve, corrupt-`[G]` WOULD fire → the
   control would catch it. ⚠ **FULL truth-table** ([[feedback-grok-final-review-pass]]): tabulate {corrupt-`[N₀]`,
   corrupt-`[D₀]`, corrupt-`[G]`, corrupt-`[c_s]`, correct} → {FAIL, FAIL, NO_FAIL, FAIL, NO_FAIL} — pin the μ̂₀-free
   default with the MIDDLE cases (the ones that would flip under a back-solved gate), not just the endpoints.
4. **The μ̂₀ DIAGNOSTIC stays NON-verdict (the demotion tooth).** Assert `dimensional_gate_depends_on_mu_hat0==False` and
   `mu_hat0_homogeneity_diagnostic.participates_in_verdict==False` are COMPUTED (the gate's `dimensional_ok` is independent of
   `mu_dim`). Per-tooth: wire `mu_dim`/`homogeneity_pass` INTO the verdict → the demotion tooth fires (the verdict must NOT
   read the diagnostic). ⚠ This is the anti-v1-tautology tooth — make it a computed structural check (the verdict's inputs
   exclude `mu_hat0`/`mu_dim`), NOT a label assertion.
5. **The `dim_of` homogeneity teeth (the sum-mismatch `DimError`).** `dim_of` raises `DimError` on a sum-dimension mismatch
   (L336–338); the `Yhat` expansion + the `lhs`/`rhs` homogeneity exercise it. Per-tooth: corrupt one `Yhat` term's power
   (e.g. `u₂ω³`) → `[u₂ω³]≠1` → the sum-mismatch/dimensionless assert fires. Assert `Yhat` dimensionless is COMPUTED
   (`dim_of` reads the fingerprint × ω-powers), not typed.
6. **The DYNAMIC 021-local self-ablation (`3d`/`3d′`) — NOT a constant (the v1 rig trip-up).** RE-SCOPE the `ablation` to an
   021-LOCAL verdict over 021's gate only (`dimensional_ok` — the μ̂₀-free gate + drop-norm-fires + corrupt-`[N₀]`-fires),
   NOT the joint `base_verdict` (which re-runs 018/019/020's gates NOT built here). It MUST be a two-verdict re-run
   (`with_mutation`=`FAIL_DIMENSIONAL` ≠ `without_mutation`=`DIMENSIONAL_OK`, `rerun_gate_logic:True`), NOT a constant.
   Per-tooth: neuter the mutation → `with_mutation`==`without_mutation` → the self-ablation flags NOT-able-to-fail.
7. **The `.wl` native-KEEP + arity integrity (light).** Def/call arity scan + unevaluated-leakage transcript scan (stage007
   lesson — the `.wl` has `Module`s + the `dimOf` recursion); confirm the `.wl` GENUINELY computes 021's dim block (native
   `dimOf`/`rawDims`/`P0Physical`/`dropNorm`/`corruptN0` + the ADDED `corruptG` control + DYNAMIC self-ablation), a truly
   independent Wolfram route (KEPT native, §1b), NOT a `.py` mirror.

⚠ **The self-ablation RE-SCOPE (like 018's 3a/3b + 019's 3g + 020's 3c/3e/3f):** the `3d`/`3d′` self-ablations MUST be
DYNAMIC re-runs (`ablation` L728–755, `with_mutation`≠`without_mutation`, `rerun_gate_logic:True`) RE-SCOPED to an
**021-LOCAL verdict** over 021's dim gate only, NOT the joint `base_verdict`. ⚠ **NO constant `self_ablation`** (the v1 rig
trip-up, `part2` L88).

⚠ **NOT 021 (do not rebuild as 021 teeth — 018/019/020 own these):** `3a_wrong_bc`/`3b_imposed_dissipation` (018, DONE);
`3g_wrong_prefactor_object` (019, DONE); `3c_wrong_scaling`/`3e_equivalence_break`/`3f_partition_mislabel` (020, DONE). ⚠
**The `54/5=2·27/5` provenance partition + the Γ5/χ_Q equivalence + the g-invariance diagnostic are 020's — 021 owns NEITHER**
(021 does the dim-homogeneity gate, NOT the algebra/provenance). 021 references `target_rhs`/`Gamma5`/`P0_physical` for their
DIMENSIONS only.

## §6 Register expectation — ⭐ likely ZERO new counted knobs; a structural edge R40 (CONFIRM)

Per the split + key-pin 7: **021 is the dimensional-closure slice; it introduces NO new calibration knob.** The honest
pre-read (⚠ CONFIRM at the register step + Codex-verify against the scripts):

- **⭐ 021 likely adds ZERO new counted knobs** (like 018/019/020 / 016 / 011/012/014). The objects it first-makes-LIVE are
  NOT new counted MODEL knobs:
  - **`μ̂₀` (`mu_hat0`)** = a FREE CARRIER / non-verdict diagnostic (`participates_in_verdict=False`). NOT a counted knob — it
    is the closure normalization back-solved for the homogeneity DIAGNOSTIC only. Do NOT count it.
  - **The SOURCED port dims `[N₀]=L⁻¹M` / `[D₀]=L⁻¹T⁻²M`** = DIMENSIONAL DATA sourced from pathA_43's density port (numerator
    `[N₀]=L⁻¹M`) + the continuity balance (`[D₀]=L⁻¹T⁻²M`) — cited provenance, NOT new knobs; the abstract `N₀/D₀` structure is
    019's (already provenance-cited). ⚠ **CONFIRM at Codex-verify:** the sourced port dims are dimensional PROVENANCE (the
    Gate-6 port scalars' dims), NOT fresh counted parameters.
  - **`c`, `G`** = already registered (`c`=the `c_γ` GR-units bridge, cited benchmark; `G=GENUINE_BLOCKED`/`external_bridge_input`,
    the `force-magnitude norm` row). 021 dim-checks them, does NOT introduce them.
- **⭐ Likely a new STRUCTURAL edge R40 (confirm the next free number at registration; R39 was 020's):** the μ̂₀-free
  `[P₀^phys]=1` dimensional closure — `[P₀^phys]=(c_s/a)²·(N₀/D₀)` is dimensionless from the SOURCED port dims `[N₀]=L⁻¹M`,
  `[D₀]=L⁻¹T⁻²M` (the μ̂₀-free gate; μ̂₀ non-verdict); the natural-units trap (`P₀=N₀/D₀` drops `(c_s/a)²`) is CAUGHT
  (`3d` `FAIL_DIMENSIONAL`); the corrupt-`[N₀]` fires / corrupt-`[G]` does NOT (the gate reads the sourced port dims, not a
  μ̂₀-back-solve). **Discharges NOTHING** (a dimensional-consistency closure / the natural-units-trap catch, NOT a reduction;
  like R37/R38/R39 for 018/019/020's fingerprint/prefactor/provenance). ⭐ Also confirm R40 records: **021 COMPLETES the joint
  `QUAD_CALIBRATED`** (018∧019∧020∧021), which stays CALIBRATED (not PASS — 020's provenance + `G=GENUINE_BLOCKED`).
- **Cited provenance (NOT re-counted):** the SOURCED port dims (pathA_43 density-port), 018's `u₂/u₄/v₅`, 019's `P0=N0/D0`,
  020's `target_rhs`/`54/5`.
- **Part-II counted CALIB set stays 6** (`{μ_η,T_w,β}`(013)+`{Vp0/ℓ_c}`(015)+`{T_Ω,β₂}`(017)); 021 adds NONE.

⚠ **Codex-verify the register** (the gate that catches an over-count that would falsely inflate — or a mislabel that would
falsely shrink — the irreducible codimension count): confirm `μ̂₀`/`[N₀]`/`[D₀]`/`c`/`G` are NOT counted as new model knobs
(μ̂₀ = free-carrier diagnostic; the port dims = sourced dimensional data; `c`/`G` already registered), and R40 is `structural`
(discharges nothing). ⭐ **This is also the leg AFTER which the scheduled MIDWAY KNOB AUDIT runs** — the Part-II gravity
sector closes at the end of Cluster A (stage 023); 021 does not trigger it, but note the register's checkpoint (register
`### ⭐ MIDWAY KNOB AUDIT` — after the entire Part-II gravity sector closes).

## Verdict tokens + honest scope

021 carries the **μ̂₀-free `[P₀^phys]=1` dimensional-closure component (4/4) of `QUAD_CALIBRATED` — the COMPLETING leg**:
the μ̂₀-FREE dimensional gate `[P₀^phys]=(c_s/a)²·(N₀/D₀)` is dimensionless from the SOURCED port dims `[N₀]=L⁻¹M`,
`[D₀]=L⁻¹T⁻²M` (`[P₀_raw]=T²`, `[(c_s/a)²]=T⁻²`, `[P₀^phys]=1`; `dimensional_ok=True`); μ̂₀ does NOT enter the verdict
(`participates_in_verdict=False`, a labelled free-carrier DIAGNOSTIC that back-solves `[μ̂₀]=L⁻¹T⁻¹M⁻¹ᐟ²` for the homogeneity
check only). The natural-units trap (the handoff's `P₀=N₀/D₀` silently drops `(c_s/a)²`) is CAUGHT — `3d_dimensional_break`
fires `FAIL_DIMENSIONAL`. The corrupt-`[N₀]` probe fires `FAIL_DIMENSIONAL` (`N₀∈P₀^phys`) while corrupt-`[G]` does NOT
(`G∉P₀^phys`) — the discriminator that the gate reads the SOURCED port dims, NOT a μ̂₀-back-solve (which would be
`[G]`-sensitive). EARNED = the μ̂₀-free gate + the natural-units-trap catch + the corrupt-`[N₀]`/`[G]` discriminator + the
`Yhat` dimensionless check + the drop-normalization + corrupt-port probes (each with a DYNAMIC 021-local self-ablation).
**021 IS dimension-checking** (the OTHER half of the 020/021 operation-level cut — 020 did algebra+provenance, 021 does the
`[·]` dim-homogeneity gate). **Consumption is PROVENANCE** — the SOURCED port dims (pathA_43 density-port) + 018's `u₂/u₄/v₅`
+ 019's `P0=N0/D0` are cited (they genuinely enter 021's self-contained dim gate, but 021 does NOT re-derive them or read
018/019's output — no theatrical dual-site; the sourced dims ARE checked, so the corrupt-`[N₀]` tooth is genuine).
⭐ **021 COMPLETES the joint `QUAD_CALIBRATED`** (018 fingerprint ∧ 019 prefactor ∧ 020 provenance-partition+CALIBRATED ∧ 021
dim closure) — the landing PRINTS the joint COMPLETE (all four Gate-4 legs earned), but the joint token STAYS
`QUAD_CALIBRATED` (calibrated, not PASS; 020's provenance + `G=GENUINE_BLOCKED`). Caveats: the actual branch a-scaling + the
numerical `(D_n,N_n)` port scalars remain Gate-6 sim-deferred (report :49); the `54/5` magnitude + `G` are CALIBRATED
(`G=GENUINE_BLOCKED`; the dim closure makes the FORM dimensionally sound, it does NOT deliver `G`). ⭐ The pathA_33 trip-ups
021 OWNS (the v1 REJECTION locus): the μ̂₀-FREE dim gate must NOT back-solve μ̂₀ to force homogeneity (the v1 tautology —
demoted to a non-verdict diagnostic); the per-probe self_ablation must be a real two-verdict re-run (NOT a constant),
re-scoped to an 021-local verdict; the corrupt-`[N₀]` probe fires `FAIL_DIMENSIONAL` while corrupt-`[G]` does NOT (the
natural-units-trap discriminator — keep genuine + able-to-fail, ADD the corrupt-`[G]` negative control the source lacks); the
`.wl` KEEPS its native dim block (unlike 020). The joint `QUAD_CALIBRATED` is COMPLETE at 021; then Cluster A continues with
022/023 (pathA_34 cross-ℓ).

## Process (unchanged, calibrated — the per-stage pipeline)

Author the II-G4d reshape directive (§1 the clean 021 slice / 4-way cut + §1b the KEEP-native `.wl` dim block (sever only
YAML, ADD corrupt-`[G]` control + DYNAMIC self-ablation) + §1c the provenance-cite consumption framing + §2 faithful cover +
§3 bridge-strip incl. sever-YAML + transcript-level agreement + §5 the μ̂₀-free-gate/natural-units-trap/corrupt-`[N₀]`/`[G]`
teeth with per-tooth ablation + the FULL truth-table on the corrupt-dim discriminator + §6 the ZERO-new-knobs + R40 register
question) → **Codex xhigh design-review** → fold to `DIRECTIVE_CLEAN` → **⭐ final Grok-4.5 headless compute-verify pass**
(Grok SymPy-verifies `[P₀^phys]=1` from the sourced dims, `[P₀_raw]=T²`, `[(c_s/a)²]=T⁻²`, the drop-norm → `T²≠1`, the
corrupt-`[N₀]` → `−[D₀]≠1`, the corrupt-`[G]` → `[P₀^phys]=1` still; watch the μ̂₀-non-verdict demotion + the DYNAMIC
self-ablation + the full truth-table) → assess/verify each Grok catch → fold → **Codex confirm-pass on the folds** →
**pre-exec USER GATE** → Codex builds the two scripts (`--sandbox danger-full-access`, background, xhigh) → dual-engine both
exit 0 (repo root + foreign CWD) → arbiter re-run → full tri-review (fidelity + adversarial-with-**per-tooth ablation**;
⭐ hunt the μ̂₀-back-solve-not-returned, the corrupt-`[G]`-negative-control genuine, the DYNAMIC self-ablation, the KEEP-native
`.wl`, and any vacuous able-to-fail / grep-dodgeable acceptance) → remediate → fresh-agent re-verify → bump counts 20→21 →
parameter register (⭐ confirm ZERO new 021 knobs + μ̂₀/`[N₀]`/`[D₀]`/`c`/`G` not counted + R40 edge; note the scheduled
MIDWAY KNOB AUDIT is after Cluster A closes, not here) + Codex-verify → note/card/`\input{stages/stage_021}` + registration →
rebuild PDF → commit + docs/memory sync (keep STATUS ▶ RESUME HERE THIN; append per-stage detail to part2 Progress).
Orchestrator authors notes/cards/LaTeX/registration; Codex codes. Target stem: `ledger_stage021_dimensional_closure` (confirm
slug at authoring).
