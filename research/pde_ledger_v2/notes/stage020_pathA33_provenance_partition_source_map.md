# II-G4c (ledger_stage020) source map — pathA_33 `54/5=2·27/5` provenance partition + the CALIBRATED verdict label (QUAD_CALIBRATED 3/4)

> Running-start prep captured 2026-07-09 (post stage019, before authoring the II-G4c reshape directive) so the directive
> can be written without re-discovery. **All line refs below VERIFIED against the current sources 2026-07-09**
> (`pathA_33_quadrupole_normalization_sympy.py` = 1351 lines; `pathA_33_quadrupole_normalization.wl` = 256 lines;
> `reports/pathA_33_quadrupole_normalization.md` = 49 lines; `directives/pathA_33_quadrupole_normalization.md` = 534 lines).
> Companion: `part2_gravity_atomic_split.md` (row 020 = the pathA_33 4-way `QUAD_CALIBRATED` split, L40; the pathA_33
> trip-ups L87–89; the cross-stage flows L106–108) and the **stage018 + stage019 source maps + reshape directives** (the
> pathA_33 4-way-split exemplars — 020 continues the SAME split, is the THIRD leg, the one that LANDS the CALIBRATED
> verdict label, and reuses the provenance-cite discipline + KEEP/AUTHOR-native-`.wl`-sever-only-YAML pattern). Build-order
> id **020**, Part II. Source top-line: **`QUAD_CALIBRATED`** — 020 lands the 3/4 component (the `54/5=2·27/5` provenance
> partition + the CALIBRATED verdict label); 018 (fingerprint, DONE `4872e8b7`) is 1/4; 019 (prefactor algebra, DONE
> `f1c426f9`) is 2/4; 021 (μ̂₀-free `[P₀^phys]=1` dim closure) is 4/4 the COMPLETING leg.

## ⭐ The FIVE headline differences from 019 (READ FIRST)

1. **⚠ 020 is the THIRD leg of the SAME 4-way pathA_33 split (018/019/020/021) — the leg that LANDS THE CALIBRATED HEADLINE,
   still a PARTIAL.** Unlike 018/019 (each an EARNED slice landing `QUAD_CALIBRATED` as a PARTIAL), 020 owns the provenance
   machinery that DETERMINES the CALIBRATED *classification* — **why the verdict is `QUAD_CALIBRATED` and not `QUAD_PASS`**:
   the assembled magnitude `54/5` and Newton's `G` are classified `external_bridge_input` by a 4-way PROVENANCE partition
   (`derived_in_gate` / `external_bridge_input` / `deferred_branch_data` / `convention`), so the verdict lands CALIBRATED.
   It is still a **PARTIAL** of the joint `QUAD_CALIBRATED` — the μ̂₀-free `[P₀^phys]=1` dimensional closure is 021, the
   COMPLETING leg. 020 = THREE earned gates: the `a⁻⁵` **scaling** (`build_scaling`), the Γ5/χ_Q **equivalence** algebra
   (`build_equivalence`), and the **provenance partition** (`build_partition` + the g-invariance diagnostic), plus the
   probes `3c`/`3e`/`3f` and the CALIBRATED verdict determination (`base_verdict` L690–694). **The 4-way cut: fingerprint
   (018, EARNED) vs prefactor algebra (019, EARNED) vs the `54/5` provenance partition + the CALIBRATED label (020) vs the
   μ̂₀-free dim closure (021).**

2. **⭐⭐ 020 IS UNITS-BEARING (the KEY reversal from 019).** UNLIKE 019 (pure abstract algebra over `{ω, D0..N4}`, NO
   `c_s`/`a`/`G`/`μ̂₀`), 020 carries the units symbols `{c_s, a, c, G}` PLUS the abstract port scalars `{N0, D0}` PLUS `χ_Q`
   PLUS the earned integer `27` — because the `54/5` magnitude lives in `target_rhs = 54·G·c_s⁵/(5·a⁵·c⁵)` and the
   Burke–Thorne bridge target `2G/(5c⁵)`. **BUT 020 does ALGEBRAIC identities over these units-bearing expressions (the Γ5
   bridge residual, the `a⁻⁵` scaling, the `54/5=2·27/5` decomposition, the provenance classification) — it does NOT do the
   DIMENSIONAL homogeneity gate `[P₀^phys]=1` (that is 021's `build_dimensions`).** ⚠ The 020/021 boundary is
   **operation-level, not expression-level**: both stages reference `P0_physical=(c_s/a)²·(N0/D0)`, `Gamma5`, `target_rhs`
   — 020 uses their VALUES in the algebraic bridge, 021 checks their DIMENSIONS. Draw the cut by OPERATION (algebra+provenance
   = 020; dimension-homogeneity = 021), NOT by which symbols appear.

3. **⭐⭐ THE CORE DELIVERABLE = the `54/5=2·27/5` SymPy-VERIFIED identity (currently a TYPED STRING) + the PROVENANCE-DRIVEN
   verdict (NOT G→λG invariance).** Two joined earned things:
   - **(a) the `54/5=2·27/5` decomposition** — currently the source `build_partition` emits it as a **typed dict/STRING**
     (`.py` L622–627: `"identity": "54/5 = 2*27/5"`, `earned_factor={27, class}`, `calibrated_factor={"2/5", class}`).
     ⭐ 020 must make this a **SymPy-VERIFIED rational identity BOUND to 020's own assembled expressions** (⚠ NOT a
     bare-literal `Rational(54,5) − 2·Integer(27)/5` tautology — Grok catch): extract `mag = compact(target_rhs/(G·c_s⁵/
     (a⁵·c⁵)))` (→ 54/5, from the assembled magnitude) + `27_from_slot = compact(a⁵/(c_s⁵·v5_slot))` (→ 27, from 018's `v₅`
     slot), then `compact(mag − 2·27_from_slot/5)==0` (exact), with the `27` CITED as 018's `1/v₅ᶻ` (`derived_in_gate`; ⚠
     **the `27` stays 018's computed value — do NOT re-derive it here**) and the `2/5` = the GR Burke–Thorne `2G/(5c⁵)`
     coefficient (`external_bridge_input`,
     `G=GENUINE_BLOCKED`). The assembled `54/5` class = `external_bridge_input` (external dominates — see the classifier,
     headline-4).
   - **(b) the verdict is PROVENANCE-DRIVEN, NOT G-invariance-driven** — the source emits a **SEPARATE**
     `g_to_lambda_g_diagnostic` (`.py` L635–663) whose whole point is the **invariance-only TRAP**: `pure_54_over_5` IS
     G-invariant (`54/5` is a pure number, `.subs(G, λ_G·G)` unchanged) YET it is `external_bridge_input` by provenance
     (`invariance_only_trap_catches_54_over_5 = g_invariant ∧ class==external_bridge_input`, L661–662). So a G→λG
     invariance test would MISLABEL the `54/5` as earned. **The verdict is set by `classify_provenance` (the tag-based
     partition), NOT by `is_g_invariant`.** The `3f_partition_mislabel` probe (L846–865) fires `FAIL_PROVENANCE_PARTITION`
     in **BOTH directions** — reclassify external→derived (`{"G": "derived_in_gate"}`, L842) AND derived→external
     (`{"fingerprint_27": "external_bridge_input"}`, L843).

4. **⭐ The consumption of 018 (χ_Q, `27`) + 019 (`P0=N0/D0`) — GENUINELY ENTERS 020's SELF-CONTAINED algebra, but stays
   PROVENANCE-cited (no theatrical cross-stage dual-site).** UNLIKE 018/019 (whose consumed inputs were provenance-only
   symbolic carriers), 020's equivalence bridge genuinely DEPENDS on 018's derived values:
   - **018's χ_Q=+1** enters `forward = target_rhs·χ·a⁵/(27c_s⁵) − 2G/(5c⁵)` (`.py` L516). With `target_rhs=54Gc_s⁵/5a⁵c⁵`:
     `forward = 2G(χ−1)/(5c⁵)`, so **`forward=0 ⟺ χ=1`** — the equivalence closes to the +2G/5c⁵ (outgoing, positive
     quadrupole) target ONLY for 018's outgoing χ_Q=+1 (an incoming χ_Q=−1 gives the wrong sign).
   - **018's `27`** (= `1/v₅ᶻ`, the radiative slot; `a⁵/27c_s⁵` IS 018's `v₅`) enters the bridge `a⁵/(27c_s⁵)`. Because
     `target_rhs` carries `54` and the bridge carries `27`, `forward=0 ⟺ 54/27=2` — i.e. the equivalence residual is
     EXACTLY the `54/5=2·27/5` decomposition. (Corrupt the `27`→`26` and `forward = 2G/130c⁵ ≠ 0`.)
   - **019's `P0=N0/D0`** enters ONLY the `Gamma5 = χ_Q·P0_physical·a⁵/(27c_s⁵)` *definition* (`.py` L512–513); the
     equivalence `ok` residual (`forward`/`reverse`, L516–517, L520) uses `target_rhs`, NOT `P0_physical` — so `N0/D0`
     does **not** enter the checkable equivalence residual (it is definitional/provenance in 020).
   - **⭐ RESOLUTION (the honest framing — same discipline as 018/019, but the values genuinely enter 020's OWN algebra):**
     cite 018's χ_Q=+1 / `v₅ᶻ=1/27` and 019's `P0=N0/D0` as **PROVENANCE**; the equivalence bridge residual (`forward`/
     `reverse`) + the `54/5=2·27/5` SymPy identity are 020's **SELF-CONTAINED earned checks** (they verify `54/27=2` and
     `χ=1`-closure as 020's own algebra — NOT a cross-stage read of 018's output). ⚠ **Do NOT manufacture a theatrical
     cross-stage dual-site on χ_Q/`27`/`N0/D0`** (a "second site" that reads 018's computed χ_Q or re-derives it would be
     subsumed / re-derive 018's DONE work — the 017 lesson `part2` L331–334 + the 018/019 §1c discipline). The genuine
     content is 020's self-contained bridge + partition; 018/019 supply the provenance of WHY χ=1, WHY 27, WHAT `P0` is.

5. **⭐⭐ The `.wl` must be GENUINELY AUTHORED for 020's content — the source `.wl` has NONE of it.** ⚠ **VERIFIED (2026-07-09
   full `.wl` read):** the source `pathA_33_quadrupole_normalization.wl` (256 lines) has ONLY the fingerprint block (L30–71,
   018), the prefactor block (L73–92, 019), and the dimensional block (L101–174, 021). It has **NO scaling probe, NO
   equivalence RESIDUAL (`forward`/`reverse`), NO provenance partition, NO g-invariance diagnostic** — it defines `Gamma5`
   (L150) + `targetRHS` (L151) but ONLY dimension-checks them (`gamma5Dim` L172), never algebraically. So — UNLIKE 018/019
   (which KEPT the `.wl`'s already-native fingerprint/prefactor blocks + severed only the YAML) — **020's `.wl` is a GENUINE
   fresh independent-engine authoring**: native Wolfram `Exponent`/`Together`/`Cancel`/`FullSimplify` for the a⁻⁵ scaling + the equivalence bridge residual (⚠ NO series expansion in the 020 slice — `Series` is NOT an independence criterion), native `Rational`
   arithmetic for `54/5=2·27/5`, native `Association`/pattern classification for the provenance partition, native
   `expr /. G -> lambdaG G` for the g-invariance diagnostic. This is a truly independent route (Wolfram-native provenance
   classification, NOT a mirror of the `.py`'s tag dicts) — the strongest form of the dual-engine requirement. (021 will keep
   the existing `.wl` dimensional block.)

## §1 The 020 slice (`.py` line ranges) — the CLEAN CUTS (all VERIFIED)

The whole computation is `build_*` helpers assembled by `build_context` (L972–1017) + `build_counterfactuals` (L703–902).
**020 owns the scaling + equivalence-algebra + provenance-partition + probes `3c`/`3e`/`3f` + the CALIBRATED verdict
determination; 018/019/021 own the fingerprint / prefactor / dimensional-homogeneity blocks.** The 020-owned cuts:

- **`a_power()` L276–279** — helper: the `a`-exponent of a factored expr (`sp.factor(expr).as_powers_dict()`).
- **⭐ `build_scaling()` L282–294** — the `a⁻⁵` target scaling (UNITS-BEARING):
  - `target_rhs = 54·G·c_s⁵/(5·a⁵·c⁵)` L283 (the assembled quadrupole magnitude); `mutated_rhs = 54·G·c_s⁵/(5·a⁴·c⁵)` L284
    (the wrong-scaling probe object, `a⁵→a⁴`).
  - `power = a_power(target_rhs)` L285 (`=−5`); `P0_target_scaling = int(power)` L289; `target_scaling_ok = (power==−5)` L290.
  - `actual_branch_scaling_status = "DEFERRED_GATE_6_SYMBOLIC_Dn_Nn_NO_A_SCALINGS_SUPPLIED"` L291 (the ACTUAL branch
    a-scaling — from `N0/D0`'s a-dependence — is Gate-6 deferred; only the TARGET rhs a-power is checked here).
  - `wrong_scaling_probe_power = int(mutated_power)` L292 (`=−4`); `wrong_scaling_probe_ok = (mutated_power != −5)` L293.
  - ⚠ **THIS `3c` TOOTH IS WEAK BOOKKEEPING (key trip-up, `part2` L88 / sub-agent §8): STRENGTHEN IT.** The `−5` is
    TYPED into `target_rhs` (`a⁵` in the denominator), so `a_power(target_rhs)=−5` is near-tautological; the mutation just
    reads `a⁴→−4`. **Strengthening (see §5 tooth 1):** tie the `−5` to the EQUIVALENCE `a`-cancellation — the a-free-ness
    of the Burke–Thorne coefficient `gamma_target=2G/(5c⁵)` (`a_power=0`) FORCES `target_rhs = gamma_target·27c_s⁵/(χ·a⁵)`
    to carry `a⁻⁵`, so the `−5` is DERIVED from a-cancellation, not read off a typed target; the wrong-scaling probe then
    breaks a-cancellation (`gamma_target·27c_s⁵/(χ·a⁴)` is not a-free ⟺ the assembled residual ≠ 0).
- **⭐ `build_equivalence(fingerprint)` L510–535** — the Γ5/χ_Q equivalence algebra (UNITS-BEARING; the `54Gc_s⁵/5a⁵c⁵ ⟺
  2G/5c⁵` bridge):
  - `chi = fingerprint["chi_Q"]` L511 (018's DERIVED χ_Q=+1 — the ONE genuine cross-stage value consumed).
  - `P0_physical = (c_s/a)²·(N0/D0)` L512 ; `Gamma5 = χ_Q·P0_physical·a⁵/(27c_s⁵)` L513 (definitions; `N0/D0`=019's, the
    `27`=018's; these DEFINE Gamma5 but do NOT enter the `ok` residual).
  - `target_rhs = 54·G·c_s⁵/(5·a⁵·c⁵)` L514 ; `gamma_target = 2·G/(5·c⁵)` L515 (the Burke–Thorne PN coefficient).
  - **`forward = compact(target_rhs·chi·a⁵/(27c_s⁵) − gamma_target)` L516 (=0 ⟺ χ=1 ∧ 54/27=2)** ; `reverse =
    compact(gamma_target·27c_s⁵/(chi·a⁵) − target_rhs)` L517 (the inverse bridge). `ok = forward==0 ∧ reverse==0` L520.
  - `wrong_gamma = 3·G/(5·c⁵)` L518 (the `3e` probe: 2/5→3/5) ; `wrong_reverse = compact(wrong_gamma·27c_s⁵/(chi·a⁵) −
    target_rhs)` L519 ; `wrong_gamma_probe_fires = not bool_from_residual(wrong_reverse)` L534 (≠0 → fires).
- **⭐ The provenance machinery L538–663 — the 020 CORE:**
  - `DERIVED_TAGS`/`EXTERNAL_TAGS`/`DEFERRED_TAGS`/`CONVENTION_TAGS` L538–549 (the tag→class dictionaries;
    `EXTERNAL_TAGS={external_gr_constant, external_pn_bridge, einstein_bridge_identity}`).
  - `PROVENANCE_SOURCES` L552–572 — per-quantity tag assignment. ⭐ KEY rows: `fingerprint_27:[dtn_radiative_slot]`
    (→derived), `PN_2_over_5:[external_pn_bridge]` (→external), `G:[external_gr_constant]` (→external),
    `assembled_54_over_5_magnitude:[external_pn_bridge, dtn_radiative_slot]` (→external, because external dominates —
    L566), `D_n_N_n_numeric_values:[gate6_branch_solve]` (→deferred), `port_scalars:[gate6_branch_solve]` (→deferred).
  - **`classify_provenance(tags)` L575–585** — the DOMINANCE order: `deferred > external > derived > convention`. ⭐ **This
    is WHY the assembled `54/5` (mixed earned+external tags) is `external_bridge_input`: external DOMINATES derived.**
  - `group_partition` L588–598 ; **`build_partition(overrides)` L601–628** — for each quantity: `computed_class` vs
    `emitted_class` (overrides for the mislabel probes), `class_matches_computed`, `ok = all(class_matches_computed)`. ⭐
    **The `decomposition_54_over_5` L622–627 is currently a TYPED dict/STRING** (`"identity":"54/5 = 2*27/5"`,
    `earned_factor={27, items["fingerprint_27"]["computed_class"]}`, `calibrated_factor={"2/5",
    items["PN_2_over_5"]["computed_class"]}`, `assembled_magnitude_class=items["assembled_54_over_5_magnitude"]["computed_class"]`)
    — **020 must COMPUTE the identity in SymPy** (headline-3a).
  - `is_g_invariant(expr)` L631–632 (`expr.subs(G, λ_G·G) − expr == 0`) ; **`build_g_invariance_diagnostic(items)` L635–663**
    — the SEPARATE diagnostic: for `{G, 2G/5c⁵, pure 54/5, fingerprint 27}` compute `g_invariant` + `provenance_class`;
    `classifier="provenance_not_g_invariance"` L660 ; `invariance_only_trap_catches_54_over_5 = pure_54_over_5.g_invariant
    ∧ pure_54_over_5.class=="external_bridge_input"` L661–662.
- **The CALIBRATED verdict determination — `base_verdict(gates, partition)` L676–694** (SHARED, but 020 owns the CALIBRATED
  landing): after all ordered gate-failures pass, L690–694 reads `g_class = partition["items"]["G"]["computed_class"]` and
  `mag_class = partition["items"]["assembled_54_over_5_magnitude"]["computed_class"]`; if BOTH `derived_in_gate` →
  `QUAD_PASS`, ELSE → **`QUAD_CALIBRATED`**. ⭐ **Since `G` and the assembled `54/5` are `external_bridge_input`, the verdict
  is `QUAD_CALIBRATED` — THIS is the CALIBRATED landing 020 owns (provenance-driven, not G-invariance-driven).**
- **The 020 probes (in `build_counterfactuals` L703–902):**
  - `partition_ablation(mutated_partition)` L764–771 — the DYNAMIC self-ablation helper for `3f` (re-runs `ablation` with
    `{"provenance_ok": mutated_partition["ok"]}`).
  - **probe `3c_wrong_scaling` L813–818** — reads `scaling["wrong_scaling_probe_power"]`(=−4), `wrong_scaling_probe_ok` →
    `FAIL_SCALING`; `self_ablation = scaling_ablation` (L759, DYNAMIC). ⚠ WEAK — strengthen (§5 tooth 1).
  - **probe `3e_equivalence_break` L837–841** — reads `equivalence["wrong_gamma_probe_residual"]`,
    `wrong_gamma_probe_fires`(2/5→3/5) → `FAIL_EQUIVALENCE`; `self_ablation = equivalence_ablation` (L761, DYNAMIC).
  - **probe `3f_partition_mislabel` L846–865** — `mutated_external = build_partition({"G":"derived_in_gate"})` L842,
    `mutated_derived = build_partition({"fingerprint_27":"external_bridge_input"})` L843; BOTH → `partition_ok=False` →
    `FAIL_PROVENANCE_PARTITION`; `g_invariance_only_would_miss_54_over_5` L855–857 (the trap); `self_ablation` =
    `external_partition_ablation` with `with_mutation_cases` for both directions L858–864 (DYNAMIC).
- **The `ablation` helper L728–755** — the DYNAMIC self-ablation engine: `base_verdict(mutated_gates)` vs
  `base_verdict(baseline_gates)`, `fail_suppressed` computed (`rerun_gate_logic:True`, `with_mutation`≠`without_mutation`).
  ⚠ **Note the scope caveat (SAME as 018's 3a/3b + 019's 3g):** `base_verdict` (L676–694) re-runs the JOINT gate set incl.
  018/019/021's gates (fingerprint/prefactor/dimensional/outgoing), NOT built here. **Re-scope the `3c`/`3e`/`3f`
  `with_mutation`/`without_mutation` to an 020-LOCAL verdict over 020's gates only** (`scaling_ok ∧ equivalence_ok ∧
  provenance_ok`, plus the CALIBRATED-vs-PASS determination from the partition), do NOT pull the joint `base_verdict`.
- **Shared helpers 020 uses (NOT cut boundaries):** `compact` L60, `bool_from_residual` L89, `numeric` L93, `hstr`. Symbol
  decls L45–57 (020 uses `omega, a, c_s, c, G, N0, D0, χ_Q` + the earned integers `27`/`54`/`5`/`2`/`λ_G`; it does NOT need
  `z, j₂/y₂` (018's), `N2/N4/D2/D4` (019's prefactor — but note `PROVENANCE_SOURCES` names `prefactor_P0_P2_P4` as a
  provenance ITEM, symbolic, no port values), `mu_hat0/mtilde0` (021's)).

- **⭐ CLEAN CUT — 020 owns `a_power`/`build_scaling` L276–294 + `build_equivalence` L510–535 + the provenance machinery
  L538–663 + probes `3c` L813–818 / `3e` L837–841 / `3f` L846–865 + `partition_ablation` L764–771 + the CALIBRATED verdict
  determination L690–694. It touches NONE of `build_fingerprint`/`build_prefactor`/`build_dimensions`. Do NOT pull in
  018/019/021 territory:**
  - **018 (DONE `4872e8b7`):** `spherical_j2`/`spherical_y2` L101–106 + `dtn_branch` L109–149 + `build_fingerprint` L152–187
    + `passivity_from_source` L666–673 + probes `3a` L774–797 / `3b` L799–812. (020 CITES 018's χ_Q=+1 + `v₅ᶻ=1/27`/`27` as
    provenance — §1c.)
  - **019 (DONE `f1c426f9`):** `build_port_moments` L190–209 + `build_prefactor` L212–273 + probe `3g` L866–875 +
    `prefactor_ablation` L762. (020 CITES 019's `P0=N0/D0` as provenance — it enters ONLY the `Gamma5` def, not the
    equivalence `ok` residual.)
  - **021:** the dim engine `DimError`/`dim_*` L297–376 + `SOURCED_N0_DIM`/`SOURCED_D0_DIM` L379–380 + `build_dimensions`
    L387–507 (`P0_raw=N0/D0` L399, `P0_physical=(c_s/a)²·(N0/D0)` L401, the μ̂₀-free `[P0_phys]` gate + the μ̂₀ diagnostic +
    drop-norm / corrupt-N0 mutations) + probes `3d` L820–826, `3d′` L828–836. ⚠ **`P0_physical`/`Gamma5`/`target_rhs` are
    RE-DEFINED in BOTH `build_equivalence` (020, algebra) AND `build_dimensions` (021, dimension) — 020 uses their VALUES in
    the bridge, 021 checks their DIMENSIONS. The cut is by OPERATION** (headline-2). Keep 020 free of ANY `dim_of`/`[·]`
    homogeneity check and any `mu_hat0`.

## §1b The `.wl` 020 slice (VERIFIED) — GENUINE INDEPENDENT AUTHORING (not KEEP-native like 018/019)

⭐ **VERIFIED (full `.wl` read 2026-07-09): the source `.wl` has NO 020 content** (headline-5). The `.wl` has only:
fingerprint L30–71 (018), prefactor L73–92 (019), dimensional L101–174 (021). It DEFINES `Gamma5` L150 + `targetRHS` L151 +
`P0Physical` L148 but ONLY dimension-checks them (`gamma5Dim` L172; `rhsDim` L164) — there is **NO `forward`/`reverse`
equivalence residual, NO `a⁻⁵` scaling probe, NO `54/5=2·27/5` identity, NO provenance partition, NO g-invariance
diagnostic** anywhere in the `.wl`.

**So the 020 `.wl` is a GENUINE fresh authoring of an independent Wolfram route** (the strongest dual-engine form — a
Wolfram-native provenance classifier, NOT a mirror of the `.py` tag dicts):
- the **`a⁻⁵` scaling** — native `Exponent[targetRHS, a]` (or the a-cancellation derivation) on its OWN `targetRHS = 54 G
  cs^5/(5 a^5 c^5)`; the wrong-scaling probe `a⁵→a⁴`.
- the **Γ5/χ_Q equivalence** — native `FullSimplify[targetRHS·chiQ·a^5/(27 cs^5) − 2 G/(5 c^5)] == 0` (forward) + the
  reverse + the `wrong_gamma` (2/5→3/5) probe; using its OWN `targetRHS`/`gammaTarget` (units-bearing; NO dim gate).
- the **`54/5=2·27/5` identity** — BOUND to `targetRHS`/`v5Slot` (⚠ NOT a bare `Simplify[54/5 − 2·27/5]`): `mag=Simplify[
  targetRHS/(G cs^5/(a^5 c^5))]`→54/5, `27FromSlot=Simplify[a^5/(cs^5 v5Slot)]`→27, `Simplify[mag − 2·27FromSlot/5]==0`, with
  the `27` cited as 018's `1/v₅ᶻ` and the `2/5` as Burke–Thorne.
- the **provenance partition** — a native `Association`-based classifier (tag-sets → class by the same dominance order
  `deferred > external > derived > convention`) WITH the dominance truth-table + key-class + tag-mutation checks, the
  per-quantity `computedClass`/`emittedClass`/`matches`, the `ok`, the source-faithful verdict (`both derived→QUAD_PASS;
  else→QUAD_CALIBRATED`) with the QUAD_PASS + REQUIRED MIXED controls, and the two mislabel mutations (`G`→derived, `27`→external).
- the **g-invariance diagnostic** — native `expr /. G -> lambdaG G` equality test for `{G, 2G/5c⁵, 54/5, 27}`, the
  `invariance_only_trap_catches_54_over_5` boolean.
- ⚠ **The bridge to sever:** `scratchDir`/`yamlOut` L19–22 + the YAML-line assembly `lines={…}` L182–251 + `Export[yamlOut,…]`
  L253 + the `headlineOk` guard L254 — DROP all of it; print-only + `fail[]`/`Exit[1]` on failure. **Dual-engine agreement =
  transcript-level** (both engines print the same `a⁻⁵`, the `forward=0`/`reverse=0` equivalence, the `54/5=2·27/5` identity
  TRUE, the provenance classes {`G`:external, `27`:derived, assembled `54/5`:external}, and both `3f` directions firing).
  **Zero file I/O.** Arity discipline (standing — def/call scan + unevaluated-leakage transcript scan; the stage007 lesson).
- ⚠ **Sample-subs (like 019's Grok catch):** if a light numeric cross-check is included, its rules may include `{a,cs,c,G}`
  (020 IS units-bearing, unlike 019) but the EARNED content is exact symbolic/rational — do NOT float-ify the `54/5=2·27/5`
  identity or the provenance classes.

## §1c The consumption resolution (READ — the honest PROVENANCE framing; the KEY 020 discipline)

⭐ **020's consumption of 018 (χ_Q, `27`) + 019 (`P0=N0/D0`) — cite as PROVENANCE; the equivalence + decomposition are
020's SELF-CONTAINED checks (headline-4). Do NOT manufacture a theatrical cross-stage dual-site.**

- **018's χ_Q=+1:** genuinely enters `forward = 2G(χ−1)/5c⁵` (closes ⟺ χ=1). Cite `stage018` (the outgoing-branch χ_Q=+1;
  `+1` gives the +2G/5c⁵ outgoing target, `−1` incoming = wrong sign). The equivalence residual is 020's OWN self-contained
  check that the outgoing sign closes the GR bridge — NOT a re-derivation of χ_Q (018's) and NOT a read-back of 018's output.
- **018's `27`** (= `1/v₅ᶻ`, the radiative slot; `a⁵/27c_s⁵`=018's `v₅`): enters the bridge `a⁵/(27c_s⁵)`; `forward=0 ⟺
  54/27=2` (the decomposition). Cite `stage018` (`v₅ᶻ=a⁵/27c_s⁵`, DONE). ⚠ **The `27` stays 018's computed value — 020 does
  NOT re-derive it** (key trip-up); it CITES it and verifies `54/5=2·27/5` as a SymPy identity + the equivalence closure.
- **019's `P0=N0/D0`:** enters ONLY the `Gamma5` def (L512–513), NOT the equivalence `ok` residual. Cite `stage019`
  (the prefactor `P0=N0/D0`); definitional/provenance in 020.
- ⚠ **NO theatrical cross-stage dual-site** on χ_Q/`27`/`N0/D0`: a "second site" that reads 018's computed χ_Q or re-derives
  it is subsumed / re-does 018's DONE work (the 017 lesson `part2` L331–334; the 018/019 §1c discipline). The genuine earned
  content is 020's self-contained bridge + `54/5=2·27/5` identity + provenance partition. **The `27` and χ_Q genuinely enter
  020's algebra (so a corrupt cited value would break `forward` — that IS the equivalence tooth 3e/`forward`), but the check
  is 020's own algebra, not a cross-stage read.** This is the honest resolution: consumed values are CITED provenance; the
  checks are self-contained.

## §2 The 020 claim-set (derive + assert; report/directive quotes)

- **The `a⁻⁵` target scaling (EARNED — the assembled-magnitude a-structure; report :36 `3c`).** `target_rhs=54Gc_s⁵/5a⁵c⁵`
  has `a`-power `−5` (report :5 "the a^-5 target scaling"). ⭐ STRENGTHEN: derive `−5` from the equivalence a-cancellation
  (the a-free `gamma_target=2G/5c⁵` forces `a⁻⁵`), not a typed target power. The `3c` probe (`a⁵→a⁴`) breaks a-cancellation.
  ⚠ HONEST: the ACTUAL branch a-scaling (from `N0/D0`) is DEFERRED Gate-6 (`actual_branch_scaling_status`, L291; report :49).
- **The Γ5/χ_Q equivalence (EARNED — the GR bridge; report :5).** `m̂₀²P₀ = 54Gc_s⁵/(5a⁵c⁵) ⟺ γ_quad^eff = 2G/(5c⁵)`
  (the Burke–Thorne target). The forward/reverse residuals vanish (`forward=0 ⟺ χ=1 ∧ 54/27=2`); the `wrong_gamma` probe
  (2/5→3/5) fires `FAIL_EQUIVALENCE` (report :39 `3e`).
- **⭐ The `54/5=2·27/5` provenance partition (EARNED — the headline; report :29–30).** The decomposition `54/5=2·27/5`:
  the `27` is `derived_in_gate` (018's fingerprint, `1/v₅ᶻ`), the `2/5`+`G` are `external_bridge_input` (GR Burke–Thorne),
  the assembled `54/5` is `external_bridge_input` (external dominates). ⭐ COMPUTE the identity in SymPy (NOT the typed
  string L622–627). The 4-way partition (`derived_in_gate`/`external_bridge_input`/`deferred_branch_data`/`convention`)
  classifies every quantity from its provenance tags via `classify_provenance` (dominance order).
- **⭐ The CALIBRATED verdict label (the LANDING; report :3, :30).** `base_verdict` returns `QUAD_CALIBRATED` (not
  `QUAD_PASS`) BECAUSE `G` and the assembled `54/5` are `external_bridge_input` (L690–694). ⭐ **Provenance-driven, NOT
  G-invariance-driven** — the SEPARATE g-invariance diagnostic shows `54/5` is G-invariant yet calibrated (the invariance-only
  trap, report :30-adjacent / `.py` L661–662): an invariance-only test would MISLABEL it as earned.
- **The able-to-fail probes (EARNED).** `3c_wrong_scaling` (a⁵→a⁴ → `FAIL_SCALING`, report :36), `3e_equivalence_break`
  (2/5→3/5 → `FAIL_EQUIVALENCE`, report :39), `3f_partition_mislabel` (BOTH directions: G external→derived AND 27
  derived→external → `FAIL_PROVENANCE_PARTITION`, report :40), each with a DYNAMIC 020-local self-ablation.
- **The 020-scoped landing (PARTIAL component — LANDS THE CALIBRATED LABEL).** `QUAD_CALIBRATED (3/4) = the a⁻⁵ scaling
  ∧ the Γ5/χ_Q equivalence ∧ the 54/5=2·27/5 provenance partition + the CALIBRATED verdict label (G + assembled 54/5 =
  external_bridge_input → CALIBRATED not PASS; the invariance-only trap avoided)`, with the fingerprint=018 (DONE), the
  prefactor algebra=019 (DONE), and the μ̂₀-free `[P₀^phys]=1` dim closure=021 (the COMPLETING leg). Do NOT print the joint
  as complete (021 completes it); do NOT re-present 018/019/021's accounting.

## §3 Reshape cost (the bridge to sever) — cross-script scratch-YAML family; AUTHOR the native `.wl`

Same family as pathA_30–34 / 018 / 019 (the cross-script runtime-YAML reshape). No argparse. The `.py`'s `build_*` helpers
are pure/self-contained, but `main` L1321–1347 writes `SYM_YAML` (L1324), reads `MMA_YAML` (L1326), writes `RESULTS_YAML`/
`REPORT_MD`/`FEED_NOTE` (L1333–1335); `compare_engines` L1024–1099 cross-checks; the `.wl` writes its scratch YAML via
`Export` L253. **Reshape = sever ALL file I/O both directions:** drop `main`'s YAML/report/feed writers + `yaml_read`/
`yaml_write` (L72–86) + `compare_engines`/`engine_summary`/`build_final_payload`/`build_report`/`build_feed_note` (`.py`);
drop the `Export` + the YAML-line assembly L182–255 (`.wl`). Each engine → standalone: print-only, `expect_zero`/
`bool_from_residual`/`expect_bool`-style asserts (`.py` local ledger idioms — `banner`/`subbanner`/`_record_pass`/
`_record_fail`/`expect_zero`/`expect_bool`/`expect_fail`, a `Verdict labels:` block, tallies, `OVERALL PASS`/nonzero exit),
`fail[]`/`Exit[1]` on failure (`.wl`). ⚠ **UNLIKE 018/019 (KEEP-native), the 020 `.wl` must be GENUINELY AUTHORED** for the
scaling/equivalence-algebra/provenance-partition/g-invariance (the source `.wl` has NONE of it — §1b/§5): a native Wolfram
route (own `Exponent`/`Together`/`Cancel`/`FullSimplify`/`Rational`/`Association`-classifier/`/. G->lambdaG G`), NOT a transliteration of the `.py`. **Dual-engine
agreement = transcript-level** (both engines print `a⁻⁵`, `forward=0`/`reverse=0`, `54/5=2·27/5` TRUE, the provenance classes,
both `3f` directions). **Zero file I/O.** Arity discipline (standing).

## §4 Consumed / exported

- **Consumes (PROVENANCE only — cite, do NOT re-derive, do NOT build a theatrical dual-site; §1c):**
  - **018's χ_Q=+1** — enters the equivalence `forward` (closes ⟺ χ=1); cite `stage018`. Self-contained check (020's bridge),
    NOT a re-derivation or read-back.
  - **018's `27`** (= `1/v₅ᶻ`, the `a⁵/27c_s⁵` radiative slot) — enters the bridge + the `54/5=2·27/5` identity; cite
    `stage018`. ⚠ The `27` stays 018's computed value; 020 does NOT re-derive it.
  - **019's `P0=N0/D0`** — enters the `Gamma5` def (not the `ok` residual); cite `stage019`. Definitional/provenance.
  - **017's ℓ=2 port kernel** (the D-lanes/`{B̃,Z̃}`) — provenance of the abstract port scalars (via 019's chain); cite
    `stage017`. NO value consumed; NO dual-site.
- **Exports (→ 021/022/027):** the `54/5=2·27/5` provenance partition + the CALIBRATED verdict label + the Γ5/χ_Q
  equivalence chain + the `a⁻⁵` target scaling → **021** (the μ̂₀-free `[P₀^phys]=1` dim closure completes the joint) +
  **022** (Gate-4 non-regression) + **027** (pathA_43 closure slot). Per the cross-stage flow (`part2` L107): "018–021
  export the Λ₂ fingerprint + χ_Q=1 + the 54/5 partition → 022 + 027." Cite the exact export contract at note-authoring.

## §5 Teeth candidates (020-specific, per-tooth ablation MANDATORY)

1. **⭐ STRENGTHEN the `a⁻⁵` scaling tooth (`3c_wrong_scaling`, `FAIL_SCALING`) — the WEAK-bookkeeping trip-up.** The source
   `3c` only reads `a_power(target_rhs)=−5` (typed `a⁵`) and mutates to `a⁴`. ⭐ **Strengthen: the `−5` must be TIED to the
   equivalence a-cancellation** — the a-free `gamma_target=2G/5c⁵` (`a_power=0`) forces `target_rhs=gamma_target·27c_s⁵/(χ·a⁵)`
   to carry `a⁻⁵`; so the `−5` is DERIVED (not a typed target power), and the wrong-scaling mutation (`a⁴`) breaks the
   a-cancellation in the assembled equivalence residual (not just a bare `a_power≠−5` compare). Per-tooth ablation: mutate the
   scaling → the a-cancellation residual fires; the correct `a⁵` does NOT.
2. **⭐ The Γ5/χ_Q equivalence tooth (`3e_equivalence_break`, `FAIL_EQUIVALENCE`).** `forward`/`reverse` residuals vanish for
   the correct `2G/5c⁵`; the `wrong_gamma` (2/5→3/5) → `wrong_reverse≠0` → `FAIL_EQUIVALENCE`. Per-tooth: mutate the bridge
   coefficient → fires; the correct object does NOT. DYNAMIC 020-local self-ablation (re-scoped, §5-note).
3. **⭐⭐ The provenance-partition tooth (`3f_partition_mislabel`, `FAIL_PROVENANCE_PARTITION`) — the CENTRAL 020 tooth.**
   Fires in **BOTH directions**: `build_partition({"G":"derived_in_gate"})` (external→derived) AND
   `build_partition({"fingerprint_27":"external_bridge_input"})` (derived→external), each → `partition_ok=False`. Per-tooth
   ablation MUST confirm BOTH fire at the partition assert. ⭐ Also assert the SEPARATE g-invariance diagnostic
   (`invariance_only_trap_catches_54_over_5=True`) — the verdict is provenance-driven, and an invariance-only test would MISS
   the `54/5`. ⚠ Ablate: force the classifier to be G-invariance-based → the `54/5` would mislabel as earned → the
   partition/verdict tooth must catch it.
4. **⭐ The `54/5=2·27/5` SymPy-identity tooth (BOUND to 020's expressions).** COMPUTE `mag=compact(target_rhs/(G·c_s⁵/
   (a⁵·c⁵)))` (→54/5), `27_from_slot=compact(a⁵/(c_s⁵·v5_slot))` (→27), `compact(mag − 2·27_from_slot/5)==0` (exact); the `27`
   cited from 018 (`derived_in_gate`), the `2/5` the Burke–Thorne coeff (`external_bridge_input`), the assembled `54/5`
   `external_bridge_input`. ⚠ **NOT a typed string** (`.py` L622–627) and ⚠ **NOT a bare-literal `Rational(54,5)−2·27/5`
   tautology** (bind to `target_rhs`/`v5_slot` — Grok catch). Per-tooth: mutate the assembled magnitude or the identity (e.g.
   `2·26/5`) → `compact≠0` → fires; mutate the earned-factor class (27→external) → the partition catches it (subsumed by
   tooth 3, but assert the decomposition's factor-classes are READ from the partition, not typed).
5. **The CALIBRATED-verdict tooth.** `base_verdict` returns `QUAD_CALIBRATED` (not `QUAD_PASS`) because G + assembled `54/5`
   are `external_bridge_input`. Per-tooth: reclassify BOTH G→derived AND magnitude→derived (a double mislabel) → the verdict
   would flip to `QUAD_PASS`, but `3f` catches the mislabel (`partition_ok=False`), so the honest verdict stays CALIBRATED.
   Assert the verdict is READ from the partition classes (computed), not typed.
6. **The `.wl` native-authoring + arity integrity (light).** Def/call arity scan + unevaluated-leakage transcript scan
   (stage007 lesson); confirm the `.wl` GENUINELY computes 020's content (native `Exponent`/`Together`/`Cancel`/`FullSimplify`/`Rational`/`Association`-classifier/
   `/. G->lambdaG G`), NOT a mirror of the `.py` tag dicts (§1b).

⚠ **The self-ablation RE-SCOPE (like 018's 3a/3b + 019's 3g):** the `3c`/`3e`/`3f` self-ablations MUST be DYNAMIC re-runs
(`ablation` L728–755, `with_mutation`≠`without_mutation`, `rerun_gate_logic:True`) RE-SCOPED to an **020-LOCAL verdict** over
020's gates only (`scaling_ok ∧ equivalence_ok ∧ provenance_ok`, plus the CALIBRATED-vs-PASS determination from the
partition), NOT the joint `base_verdict` (which re-runs 018/019/021's gates NOT built here). ⚠ **NO constant `self_ablation`**
(the v1 rig trip-up, `part2` L88).

⚠ **NOT 020 (do not rebuild as 020 teeth — 018/019/021 own these):** `3a_wrong_bc`/`3b_imposed_dissipation` (018, DONE); `3g_wrong_prefactor_object` (019, DONE); `3d`/`3d′` dimensional-break (021). ⚠ **The μ̂₀-free `[P₀^phys]=1` gate + the
`c_s/a` dimensional restoration (`.py` `build_dimensions`) are 021's — 020 owns NEITHER** (020 is units-bearing but does
ALGEBRA, not the dim-homogeneity gate).

## §6 Register expectation — ⭐ likely ZERO new counted knobs; a structural edge R39 (CONFIRM)

Per the split + key-pin 5: **020 is the CALIBRATED-classification slice; it introduces NO new calibration knob.** The honest
pre-read (⚠ CONFIRM at the register step + Codex-verify against the scripts):

- **⭐ 020 likely adds ZERO new counted knobs** (like 018/019 / 016 / 011/012/014). The units symbols it first-makes-LIVE in
  Part II — `G` and `c` — are NOT new counted MODEL knobs:
  - **`G` (Newton's constant)** = `GENUINE_BLOCKED` / `external_bridge_input` — the thing the PDE does NOT deliver (already in
    the register as `G=GENUINE_BLOCKED`; the `force-magnitude norm` row, II-002). NOT a new knob; 020 CLASSIFIES it, does not
    introduce it.
  - **`c` (the GR/PN vacuum light speed in `2G/5c⁵` and `54Gc_s⁵/5a⁵c⁵`)** = the GR-matching / λγ units bridge (`P₀ ∝
    c_s⁵/c⁵ = 1/λγ⁵`; plan §1). In the model `c` = the light cone `c_γ` (already a register row, `FREE-UNREDUCED` with the
    Route-A/cone-lock reduction). ⚠ **CONFIRM at Codex-verify:** `c` here is the benchmark GR/PN light-speed CARRIER of the
    Burke–Thorne target (cited from the PN corpus, like 018's `c_s` was cited PROVENANCE), NOT a new counted knob — it is the
    `c_γ` in its GR-units-bridge role, or a benchmark units carrier. Do NOT let 020 count a fresh `c`.
  - **the `2/5`** = the GR Burke–Thorne coefficient (`external_bridge_input`, a benchmark bridge); **the `27`** = 018's
    derived fingerprint (`derived_in_gate`, cited); **`54/5`** = the assembled magnitude (`external_bridge_input`). NONE is a
    new knob.
- **⭐ Likely a new STRUCTURAL edge R39 (confirm the next free number at registration; R38 was 019's):** the `54/5=2·27/5`
  provenance-partition + the CALIBRATED-verdict provenance — the assembled quadrupole magnitude `54/5` and `G` are
  `external_bridge_input` (external dominates the mixed earned+external tags), the `27` is `derived_in_gate` (018's), so the
  verdict lands `QUAD_CALIBRATED` not `QUAD_PASS`; the SEPARATE g-invariance diagnostic shows the invariance-only trap (`54/5`
  is G-invariant yet calibrated). **Discharges NOTHING** (a classification landing / the honest CALIBRATED verdict, NOT a
  reduction; like R37/R38 for 018/019's fingerprint/prefactor). Confirm the edge captures: the provenance-dominance rule
  (external > derived), the invariance-only trap, and `G=GENUINE_BLOCKED`.
- **Cited provenance (NOT re-counted):** 018's χ_Q=+1 / `v₅ᶻ=1/27` / `27`, 019's `P0=N0/D0`, 017's port kernel D-lanes.
- **Part-II counted CALIB set stays 6** (`{μ_η,T_w,β}`(013)+`{Vp0/ℓ_c}`(015)+`{T_Ω,β₂}`(017)); 020 adds NONE.

⚠ **Codex-verify the register** (the gate that catches an over-count that would falsely inflate — or a mislabel that would
falsely shrink — the irreducible codimension count): confirm `G`/`c`/`2/5`/`27`/`54/5` are NOT counted as new model knobs
(they are external bridge inputs / cited benchmarks / 018's derived value), and R39 is `structural` (discharges nothing).

## Verdict tokens + honest scope

020 carries the **`54/5=2·27/5` provenance-partition + the CALIBRATED-verdict-label component (3/4) of `QUAD_CALIBRATED`**:
the assembled quadrupole magnitude decomposes as `54/5 = 2·27/5` (a SymPy-verified identity), where the `27` is
`derived_in_gate` (018's fingerprint `1/v₅ᶻ`, cited) and the `2/5`+`G` are `external_bridge_input` (the GR Burke–Thorne
`2G/5c⁵`, `G=GENUINE_BLOCKED`); a 4-way PROVENANCE partition classifies every quantity from its provenance tags (dominance:
deferred > external > derived > convention), so the assembled `54/5` is `external_bridge_input` (external dominates) and the
verdict lands **`QUAD_CALIBRATED` not `QUAD_PASS`** (`base_verdict` L690–694). ⭐ **The verdict is PROVENANCE-driven, NOT
G→λG-invariance-driven** — the SEPARATE g-invariance diagnostic shows `54/5` is G-invariant yet calibrated (the invariance-only
trap that an invariance test would fall into). EARNED = the `a⁻⁵` target scaling (strengthened via the equivalence
a-cancellation), the Γ5/χ_Q equivalence `54Gc_s⁵/5a⁵c⁵ ⟺ 2G/5c⁵` (closes ⟺ χ=1 ∧ 54/27=2), and the provenance partition +
its able-to-fail (`3c`/`3e`/`3f`, `3f` firing BOTH directions). **020 IS units-bearing** (`{c_s,a,c,G}`+abstract`{N0,D0}`+χ_Q+`27`)
but does ALGEBRA, not the dim-homogeneity gate (021's `[P₀^phys]=1`). **Consumption is PROVENANCE** — 018's χ_Q=+1/`27` +
019's `P0=N0/D0` are cited (they genuinely enter 020's self-contained equivalence bridge, but 020 does NOT re-derive them or
read 018/019's output — no theatrical dual-site). Caveats: the ACTUAL branch a-scaling + the numerical `(D_n,N_n)` port
scalars remain Gate-6 sim-deferred (report :49); the `54/5` magnitude + `G` are CALIBRATED (`G=GENUINE_BLOCKED`). ⭐ The
pathA_33 trip-ups 020 OWNS: make `54/5=2·27/5` a SymPy identity NOT the typed string (L622–627); the verdict from PROVENANCE
NOT g-invariance; STRENGTHEN the weak `3c` a⁻⁵ scaling tooth; the v1 rig (back-solved `μ̂₀` + constant `self_ablation`) must
NOT return — 020 owns NO `μ̂₀` (021's) and its `3c`/`3e`/`3f` self-ablations are DYNAMIC re-runs re-scoped to an 020-local
verdict. The joint `QUAD_CALIBRATED` COMPLETES at 021 (the μ̂₀-free dim closure).

## Process (unchanged, calibrated — the per-stage pipeline)

Author the II-G4c reshape directive (§1 the clean 020 slice / 4-way cut + §1b the GENUINELY-AUTHORED native `.wl` +
§1c the provenance-cite consumption framing + §2 faithful cover + §3 bridge-strip incl. sever-YAML + transcript-level
agreement + §5 the scaling/equivalence/provenance teeth with per-tooth ablation + STRENGTHEN 3c + §6 the ZERO-new-knobs +
R39 register question) → **Codex xhigh design-review** → fold to `DIRECTIVE_CLEAN` → **⭐ final Grok-4.5 headless
compute-verify pass** (Grok SymPy-verifies `54/5=2·27/5`, the Γ5 bridge `54/27=2` closure, the `classify_provenance`
dominance → assembled `54/5`=external, and the g-invariance trap `54/5` invariant-yet-calibrated; watch the identity
COMPUTED-not-typed + the 3f-both-directions + the strengthened 3c a-cancellation) → assess/verify each Grok catch → fold →
**Codex confirm-pass on the folds** → **pre-exec USER GATE** → Codex builds the two scripts (`--sandbox danger-full-access`,
background, xhigh) → dual-engine both exit 0 (repo root + foreign CWD) → arbiter re-run → full tri-review (fidelity +
adversarial-with-**per-tooth ablation**; ⭐ hunt the `54/5=2·27/5`-COMPUTED-not-typed, the provenance-vs-invariance verdict,
the 3f-both-directions, the strengthened 3c, a mirror-`.wl`, and any vacuous able-to-fail) → remediate → fresh-agent
re-verify → bump counts 19→20 → parameter register (⭐ confirm ZERO new 020 knobs + G/c/2-5/27/54-5 not counted + R39 edge) +
Codex-verify → note/card/`\input{stages/stage_020}` + registration → rebuild PDF → commit + docs/memory sync. Orchestrator
authors notes/cards/LaTeX/registration; Codex codes. Target stem: `ledger_stage020_provenance_partition` (confirm slug at
authoring).
