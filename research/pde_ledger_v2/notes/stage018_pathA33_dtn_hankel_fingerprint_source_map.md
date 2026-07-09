# II-G4a (ledger_stage018) source map — pathA_33 DtN Hankel fingerprint + χ_Q sign (QUAD_CALIBRATED 1/4, the EARNED slice)

> Running-start prep captured 2026-07-09 (post stage017, before authoring the II-G4a reshape directive) so the directive
> can be written without re-discovery. **All line refs below VERIFIED against the current sources 2026-07-09**
> (`pathA_33_quadrupole_normalization_sympy.py` = 1351 lines; `pathA_33_quadrupole_normalization.wl` = 256 lines;
> `reports/pathA_33_quadrupole_normalization.md` = 49 lines; `directives/pathA_33_quadrupole_normalization.md` = 457 lines).
> Companion: `part2_gravity_atomic_split.md` (rows 018–021 = the pathA_33 4-way `QUAD_CALIBRATED` split + the pathA_33
> trip-ups L87–89 + the cross-stage flows L106–108) and the stage016/017 source maps (the pathA_32 EARNED-first /
> completing exemplars — 018 follows 016's EARNED-FIRST PARTIAL-landing pattern). Build-order id **018**, Part II.
> Source top-line: **`QUAD_CALIBRATED`** — 018 lands the 1/4 EARNED component (the outgoing fingerprint shape + χ_Q sign);
> 019 (prefactor algebra) + 020 (54/5 partition + the CALIBRATED label) + 021 (μ̂₀-free dim closure) complete the joint.
> ⚠ **No prior source map existed for stage018** — this is the FIRST authored artifact of the stage (per the calibrated
> pipeline: source map → reshape directive → Codex→Grok→Codex reviews → gate).

## ⭐ The FOUR headline differences from 017 (READ FIRST)

1. **⚠ pathA_33 splits 4-WAY (018/019/020/021), and 018 is the EARNED-FIRST slice** — the same pattern as 016 (the
   EARNED angular covariance theorem, PARTIAL landing) was in the pathA_32 2-way split. 018 = the **outgoing ℓ=2 DtN
   spherical-Hankel fingerprint**: the dimensionless-z series coefficients `u₂ᶻ=1/9`, `u₄ᶻ=4/81`, `v₅ᶻ=1/27` (physical
   `u₂=a²/9c_s²`, `u₄=4a⁴/81c_s⁴`, `v₅=a⁵/27c_s⁵`) DERIVED from `Ŷ₂ᵒᵘᵗ=−3/Λ₂ᵒᵘᵗ`, `Λ₂ᵒᵘᵗ=z·h₂⁽¹⁾′/h₂⁽¹⁾`, **plus the
   outgoing-vs-incoming sign `χ_Q=+1` (outgoing) vs `−1` (incoming)** and the standing-branch contrast (no radiating
   term). 019 = the prefactor algebra `P(ω)=D₀N/D^cons²`; 020 = the `54/5=2·27/5` provenance partition **and the
   CALIBRATED verdict label** (the `27` COMPUTED / `2·G/5` `external_bridge_input`); 021 = the μ̂₀-free `[P₀^phys]=1`
   dim closure. **The cut is exterior-wave fingerprint (018, EARNED) vs prefactor (019) vs magnitude-provenance (020,
   the CALIBRATED landing) vs dim-closure (021).**
2. **⭐ 018 is SELF-CONTAINED exterior spherical-Hankel algebra — it does NOT literally consume 017's port kernel or
   009/010's bulk mode (a REFINEMENT of the running-start key pin 1).** The fingerprint core (`.py` L101–187) builds the
   ℓ=2 exterior outgoing solution `h₂⁽¹⁾=j₂+i·y₂` FROM explicit `j₂`/`y₂` (L101–106) and series-expands `−3/(z h′/h)` —
   it never references `build_port_moments` (L190–209, the port/wall-mode object, which carries the symbolic `N_n/D_n` →
   **019's** territory) nor any stored bulk-Helmholtz artifact. So **018 CITES 017's ℓ=2 port kernel + 009/010's bulk
   Helmholtz mode as PROVENANCE** (the physical two-port setup: 018's outgoing DtN `Ŷ₂ᵒᵘᵗ` is the EXTERIOR radiative
   response of the ℓ=2 channel whose INTERIOR wall mode 017 exports and whose bulk mode 009/010 export; the literal
   port-kernel D-lanes / N-moment symbols enter at **019**). ⚠ **Do NOT manufacture a theatrical cross-stage dual-site
   on objects the fingerprint does not use** — provenance cite only, like 016's provenance-only cite of 013 (the
   [[reference-conceptual-foundation-doc]] "cite as provenance across a boundary" lesson).
3. **⭐ c_s is a units-restoring free symbol, NOT a live consumed VALUE (a REFINEMENT of the running-start key pin 2).**
   The EARNED content — `u₂ᶻ=1/9`, `u₄ᶻ=4/81`, `v₅ᶻ=1/27` (`.py` L129–131) and `χ_Q=±1` (L168–169) — is **c_s-FREE**
   (dimensionless z-space rationals; `χ_Q`'s `c_s⁵` cancels). `c_s` is declared a bare positive free symbol (L46, NOT
   defined via `5Kρ⁴/m` in-file) and enters ONLY as the units-restoring carrier attached at L146–148 (`u₂=u₂ᶻ·a²/c_s²`,
   via `z=aω/c_s`). ⚠ So **018 carries `c_s` as an inert units symbol whose physical identity is the density sound speed
   (cite stage005 R1 `c_s²=5Kρ⁴/m` as PROVENANCE), but does NOT consume its value** — the fingerprint would be identical
   for any `c_s>0`. This is a softer resolution than "c_s becomes a live consumed symbol"; it keeps 018 an EARNED
   fingerprint slice with **ZERO new counted knobs** (`c_s` is R1-DERIVED, and here merely a units carrier).
4. **⭐⭐ The pathA_33 trip-ups live at 020/021, but 018 owns the "27-stays-COMPUTED / χ_Q-not-stamped" obligation.**
   From the rig-history (`part2_gravity_atomic_split.md` L87–89; directive L23–26): the v1 rig (a) **back-solved the FREE
   carrier `μ̂₀`** to force homogeneity, and (b) used a **constant `self_ablation`** not a re-run — BOTH are 021's / the
   probe-engine's concern. For 018: (i) **the `27` stays COMPUTED** — the earned `27` is `v₅ᶻ=1/27` series-derived at
   L131 (NOT the typed reference denominator `a⁵/(27c_s⁵)` at L168, which is only the canonical slot the derived `v₅`
   is compared against); (ii) **`χ_Q` must be COMPUTED, not stamped** — `χ_derived=out["v5"]/(a⁵/27c_s⁵)` L168 (a typed
   `χ_Q=+1` ⇒ tautology, directive §4/§2.4); (iii) **018's fingerprint slice must stay free of any `μ̂₀` back-solve**
   (it already is — `μ̂₀` enters only `build_dimensions`, 021's; the fingerprint slice has ZERO `μ̂₀` refs, sub-agent
   §5-confirmed). The μ̂₀-free dim gate + the real two-verdict self-ablations are 021's, not 018's.

## §1 The 018 slice (`.py` line ranges) — the CLEAN CUTS (all VERIFIED)

The whole computation is a set of `build_*` helpers assembled by `build_counterfactuals` (L703–902) + `main` (L1321–1347).
**018 owns the fingerprint core + its two able-to-fail probes; 019/020/021 own the prefactor/partition/dimension blocks.**
The 018-owned cuts:

- **`spherical_j2()` L101–103 + `spherical_y2()` L105–106** — the explicit ℓ=2 `j₂(z)=(3/z³−1/z)sin z−3cos z/z²`,
  `y₂(z)=(1/z−3/z³)cos z−3 sin z/z²`. The physical basis of the exterior wave.
- **⭐ `dtn_branch(kind)` L109–149** — the core DtN log-derivative machinery:
  - `h = j₂ + i·y₂` (outgoing L112–113) / `h = j₂ − i·y₂` (incoming L115–116) / `h = j₂` (standing L118–119);
  - **`lam = compact(z·diff(h,z)/h)` L124** — the DtN eigenvalue `Λ₂=z h′/h`;
  - **`yhat = compact(−3/lam)` L125** — `Ŷ₂=−3/Λ₂`;
  - the Taylor series `h_series`/`lam_series`/`yhat_series` L126–128 (via `series(...).removeO()`);
  - **`u2_z=yhat_series.coeff(z,2)` L129, `u4_z=...coeff(z,4)` L130, `v5_z=...coeff(z,5)/I` L131** (the radiative part
    extracted by `/I`); `static=...coeff(z,0)` L132;
  - the physical (units-restored) coefficients `u2=u2_z·a²/c_s²`, `u4=...·a⁴/c_s⁴`, `v5=...·a⁵/c_s⁵` L146–148.
- **⭐ `build_fingerprint()` L152–187** — assembles the three branches, the derive-then-check, and χ_Q:
  - `out=dtn_branch("outgoing")`, `incoming=...("incoming")`, `standing=...("standing")` L153–155;
  - the EXPECTED targets `expected={u2:a²/9c_s², u4:4a⁴/81c_s⁴, v5:a⁵/27c_s⁵}` L157–162 (typed reference — the
    derived `out[...]` and typed `expected[...]` are INDEPENDENT, so L164–167 `bool_from_residual(out[k]−expected[k])`
    is a genuine derive-then-check that could disagree);
  - **`chi_derived=compact(out["v5"]/(a⁵/(27c_s⁵)))` L168** (=+1 iff `v5_z=1/27`), **`chi_incoming=...incoming[...]` L169**
    (=−1, the incoming z⁵ coefficient is `−i/27`); `sample subs {a:3, c_s:2}` L170;
  - `ok=all(matches) and bool_from_residual(chi_derived−1)` L177; `chi_Q=chi_derived` L178; `chi_Q_incoming=chi_incoming`
    L179. ⭐ **The ±1 sign EMERGES from `j₂±i·y₂` (L113 vs L116) propagating through `Λ₂` + the `/I` extraction — not
    hardcoded** (cross-checked in the counterfactual at L719–720: `bool_from_residual(incoming["v5_z"]+out["v5_z"])`,
    `bool_from_residual(fingerprint["chi_Q_incoming"]+1)`).
- **`passivity_from_source()` L666–673** — the outgoing-vs-incoming "genuine_outgoing" gate (the §2.7 passivity claim:
  the radiating imaginary term arises from the outgoing BC, not an inserted sink). 018-owned.
- **The two 018 able-to-fail probes (inside `build_counterfactuals` L703–902):**
  - **`3a_wrong_bc` L774–797** — incoming (`χ_Q=−1`) + standing (no `v5`, `Λ_stand(0)=+2`) recompute; asserts the tuple
    changes in the predicted way (`FAIL_FINGERPRINT`). Routes through `fingerprint_ablation` L757 (a DYNAMIC re-run).
  - **`3b_imposed_dissipation` L799–812** — a phenomenological damping term vs the genuine `dtn_hankel1` outgoing BC;
    asserts the inserted-sink radiative part is DETECTED (`FAIL_NOT_OUTGOING`) and that removing the outgoing BC removes
    the imaginary `v5`. Routes through `outgoing_ablation` L758 (DYNAMIC).
- **Shared helpers 018 uses (NOT cut boundaries):** `compact` L60, `series_no_o` L97, `bool_from_residual` L89,
  `numeric` L93. Symbol decls L45–57 (018 uses `z, omega, a, c_s`; NOT `D0..N4, mu_hat0, mtilde0`).

- **⭐ CLEAN CUT — 018 owns L101–187 (fingerprint core) + L666–673 (passivity) + the 3a/3b probes (L774–812) + the
  two ablation helpers L757–758. It touches NONE of `build_dimensions`/`build_prefactor`/`build_scaling`/`build_partition`/
  `build_equivalence`. Do NOT pull in 019/020/021 territory:**
  - **019:** `build_port_moments` L190–209 (symbolic `N_n/D_n` port scalars, deferred/Gate-6) + `build_prefactor` L212–273
    (`P0/P2/P4` from `D₀N/D^cons²`, the squared-denominator algebra + N/D self-check) + probe `3g` L866–875.
  - **020:** `a_power`/`build_scaling` L276–294 (a⁻⁵ target scaling) + `build_equivalence` L510–535 (`Γ₅/χ_Q ↔ 2G/5c⁵`)
    + the provenance machinery L538–663 (`classify_provenance` L575–585, `build_partition` L601–628 incl. the STRING
    `decomposition_54_over_5` L622–627, the g-invariance diagnostic L635–663) + probes `3c` L813–818, `3e` L837–841,
    `3f` L846–865.
  - **021:** the dim engine `DimError`/`dim_*` L297–376 + `SOURCED_N0_DIM`/`SOURCED_D0_DIM` L379–380 + `build_dimensions`
    L387–507 (the μ̂₀-free `[P0_phys]` gate + the μ̂₀ diagnostic + drop-norm / corrupt-`N0` probes) + probes `3d` L820–826,
    `3d′` L828–836.
  - **⚠ 018 CITES 017's port kernel + 009/010's bulk mode as PROVENANCE (§4); the literal `N_n/D_n` consumption is
    019's.** Do NOT let 018 re-present the prefactor / magnitude / dimension accounting.

## §1b The `.wl` 018 slice (VERIFIED) — the independent Wolfram route (KEEP native, sever ONLY the YAML)

⭐ **The pathA_33 `.wl` is ALREADY a genuinely independent engine** (native `Series`/`Coefficient`/`FullSimplify` on
native `j2`/`y2` spherical-Bessel expressions, its OWN `dimOf`/`serZ`/`serW`, its OWN prefactor algebra; it does **NOT**
`Get`/`Import` the `.py`'s expressions). It is **NOT a payload-mirror** — the ONLY bridge is the scratch-YAML `Export` at
L253. The reshape KEEPS the native route and severs ONLY that YAML handoff (§3). 018's `.wl` slice:
- `j2` L30, `y2` L31 (native ℓ=2 spherical Bessel); `serZ` L27 (native z-series).
- **⭐ `branchData[name,h,source]` L33–60** — the native DtN core: `lam=z D[h,z]/h` L35, `yhat=−3/lam` L36, series
  L37–39, `u2z=Coefficient[yser,z,2]` L40 / `u4z=...z,4]` L41 / `v5z=Coefficient[yser,z,5]/I` L42 / `static=...z,0]` L43,
  and the physical `u2=u2z·a²/cs²` / `u4` / `v5` L56–58.
- `out=branchData["outgoing_hankel1", j2+I·y2, ...]` L62, `incoming=...j2−I·y2...]` L63, `standing=...j2...]` L64.
- **The fingerprint asserts L66–71:** `u2Match=TrueQ[out["u2z"]==1/9]` L66, `u4Match` L67, `v5Match` L68 (⚠ derived-vs-
  typed — the derived `out["u2z"]` is genuinely series-computed; see §5 the per-tooth-ablation note), `chiQ=out["v5"]/
  (a⁵/(27cs⁵))` L69, `chiQIncoming=incoming["v5"]/...` L70, `chiMatch=TrueQ[chiQ==1]` L71.
- **⚠ 019/021 territory in the `.wl` (EXCLUDE from 018):** the prefactor block `Nomega`/`Dcons`/`prefObj`/`p0`/`p2`/`p4`/
  `plainP2` L73–92 (→ 019); the dimensional block `zeroDim`/`dimOf`/`rawDims`/`P0Physical`/`Gamma5`/`targetRHS`/`muDim`/
  the drop-norm + corrupt-N0 probes L101–174 (→ 021); the combined `headlineOk` L176–180 (rebuild an 018-scoped headline
  over `u2Match∧u4Match∧v5Match∧chiMatch` only); the standing-branch static-slot lines are part of the YAML payload L208–210.
- **⚠ The bridge to sever (§3): the YAML writer L182–255** (`lines={...}` L182–251 assembling the scratch YAML, the
  `Export[yamlOut,...]` L253, the `headlineOk` guard L254, `scratchDir`/`yamlOut` setup L19–22). SEVER: no scratch YAML,
  print-only + `fail[]` on failure (the `.wl` already has `fail[]` L5). **Dual-engine agreement = transcript-level** (both
  engines print the same derived `u2z=1/9`/`u4z=4/81`/`v5z=1/27`/`χ_Q=+1`/incoming `χ_Q=−1`/standing static slot; the
  stage014/016/017 transcript pattern). **Zero file I/O.** Arity discipline (standing — def/call scan + unevaluated-leakage
  transcript scan; the stage007 lesson).

## §1c The consumption resolution (READ — the honest c_s + provenance framing)

⭐ **Both running-start pins (1) and (2) RESOLVE to provenance-only, not literal consumption** (sub-agent §3/§4-confirmed):

- **`c_s` (pin 2):** carried as a bare positive free symbol (units carrier via `z=aω/c_s`); the earned rationals + χ_Q are
  c_s-FREE. **Cite stage005 R1 (`c_s²=5Kρ⁴/m`) as PROVENANCE of what `c_s` IS** (the density sound speed — the first
  time the density sound speed appears as an object in the Part-II gravity radiative sector, a physical shift from the
  013–017 breathing/isotropy sector which deferred the matter mode `kξ≪1`), but 018 does NOT consume its VALUE. ⚠ Do NOT
  build a theatrical `c_s²=5Kρ⁴/m` dual-site — the fingerprint is value-independent of `c_s`; R1 is a one-line provenance
  cite, not a checkable in-stage relation. `[c_s]=LT⁻¹` (used in the units-restored physical coefficients — a genuine
  dim leg, able-to-fail if a coefficient's `c_s`-power is corrupted).
- **017's ℓ=2 port kernel + 009/010's bulk Helmholtz mode (pin 1):** cite as PROVENANCE (the physical two-port setup —
  018's outgoing DtN `Ŷ₂ᵒᵘᵗ` is the EXTERIOR radiative response of the ℓ=2 channel; 017's grouped-P2 wall mode is the
  INTERIOR, 009/010's is the bulk mode; the exterior/interior are MATCHED in **019's** prefactor `P=D₀N/D^cons²`, where
  the port kernel's `N_n/D_n` literally enter). 018 self-derives the exterior Hankel fingerprint from `j₂`/`y₂`; it does
  NOT read the D-lanes. ⚠ The `−3=−(ℓ+1)` outgoing static DtN slot and 016/017's `λ_m=6=ℓ(ℓ+1)` angular eigenvalue both
  encode ℓ=2 but are DIFFERENT objects (exterior DtN slot vs angular Laplacian eigenvalue) — 018 self-derives its `−3`;
  do NOT tie it to 017's `λ_m` (provenance only: both live in the ℓ=2 sector 017's port kernel defines).
- ⚠ **`c_S` (the frozen-wall Helmholtz speed, stages 011–017) is a DISTINCT object** from `c_s` (the bulk density sound
  speed here) — same R1 functional form but evaluated at different densities (`ρ*` wall vs `ρ0` bulk). 018's `c_s` is the
  bulk density sound speed. Keep the symbols distinct (an R27-style firewall note; do NOT substitute `c_S`→`c_s`).

## §2 The 018 claim-set (derive + assert; report/directive quotes)

- **The outgoing ℓ=2 DtN Hankel fingerprint (EARNED — the headline).** Series-expanding `Ŷ₂ᵒᵘᵗ=−3/Λ₂ᵒᵘᵗ`,
  `Λ₂ᵒᵘᵗ(z)=z·h₂⁽¹⁾′(z)/h₂⁽¹⁾(z)`, `z=aω/c_s` about `z=0` yields the dimensionless coefficients
  `u₂ᶻ=1/9, u₄ᶻ=4/81, v₅ᶻ=1/27` (report :9 SymPy `u2=a**2/(9*c_s**2)`, `u4=4*a**4/(81*c_s**4)`, `v5=a**5/(27*c_s**5)`;
  directive §2.1). ⭐ **The `27` (the radiative slot's denominator) is EARNED here as `v₅ᶻ=1/27` series-derived — this
  is the `derived_in_gate` "27" that 020's `54/5=2·27/5` partition rides.** The even `u₂,u₄` are the reactive static-
  response coefficients; the imaginary `v₅ω⁵` is the RADIATING part (the `/I` extraction, L131).
- **`χ_Q=+1` outgoing vs `−1` incoming (EARNED — the sign classification).** `χ_Q = derived v₅ / (a⁵/27c_s⁵)` comes out
  `+1` on the outgoing `h₂⁽¹⁾` branch and `−1` on the incoming `h₂⁽²⁾` branch (report :11 `Derived chi_Q: 1; incoming
  gives -1`; directive §2.4). ⭐ **COMPUTED, not stamped** — the sign propagates from `j₂±i·y₂`; a typed `χ_Q=1` is a
  tautology the firewall forbids (directive §4). The standing `j₂` branch has `Λ_stand(0)=+2` (not the outgoing `−3`),
  `Ŷ_stand(0)=−3/2≠1`, and NO imaginary term — proving `+1/27` is outgoing-BC-selected, not universal.
- **Passivity / not-static (EARNED able-to-fail).** The radiating imaginary term arises from the genuine outgoing-wave BC
  (`passivity_from_source` L666–673); `P₀=Ŷ₂ᵒᵘᵗ(ω→0)` is the static limit of the SAME outgoing response, not a static
  gain with damping appended (directive §2.7). Inserting a phenomenological dissipation instead ⇒ `FAIL_NOT_OUTGOING`
  (§3b).
- **The dim leg (EARNED able-to-fail, `c_s`-carrying).** `[u₂]=[a²/c_s²]=T²`, `[v₅]=[a⁵/c_s⁵]=T⁵` with `[a]=L`,
  `[c_s]=LT⁻¹` — a genuine units-restored dim check on the physical coefficients (corrupting a coefficient's `a`/`c_s`
  power fires). ⚠ **NOT** the μ̂₀-free `[P₀^phys]=1` gate (that is 021's — 018 does not touch `P₀`, `μ̂₀`, `N₀`, `D₀`).
- **The 018-scoped landing (PARTIAL component).** Land at the 018 component: `QUAD_CALIBRATED (1/4) = the outgoing ℓ=2
  DtN Hankel fingerprint EARNED (u₂=a²/9c_s², u₄=4a⁴/81c_s⁴, v₅=a⁵/27c_s⁵ + χ_Q=+1 outgoing / −1 incoming + standing
  contrast), with the prefactor algebra = 019, the 54/5 provenance partition + the CALIBRATED label = 020, and the
  μ̂₀-free dim closure = 021.` Do NOT print the joint as complete (that is 020's CALIBRATED landing) and do NOT re-present
  019/020/021's accounting (cite the port kernel as PROVENANCE — §4). ⭐ Follows 016's EARNED-FIRST PARTIAL-landing
  pattern (016 printed `ISOTROPY_CALIBRATED (1/2) — SO(3) covariance theorem EARNED (PARTIAL)`).

## §3 Reshape cost (the bridge to sever) — cross-script scratch-YAML family, KEEP the native `.wl`

Same family as pathA_30–34 (the cross-script runtime-YAML reshape, NOT the sympy-expr-import family). No argparse (the
two-phase behavior is driven purely by presence of the MMA scratch file, `.py` L1326–1330). The `.py`'s `build_*`
helpers are pure/self-contained, but `main` L1321–1347 writes `SYM_YAML` (L1324), reads `MMA_YAML` (L1326), and writes
`RESULTS_YAML`/`REPORT_MD`/`FEED_NOTE` (L1333–1335); `compare_engines` L1024–1099 cross-checks; the `.wl` writes its
scratch YAML via `Export` L253. **Reshape = sever ALL file I/O both directions:** drop `main`'s YAML/report/feed writers +
`yaml_read`/`yaml_write` (L72–86) + `compare_engines`/`engine_summary`/`build_final_payload`/`build_report`/
`build_feed_note` (`.py`); drop the `Export` + the YAML-line assembly L182–255 (`.wl`). Each engine → standalone:
print-only, `expect_zero`/`bool_from_residual`-style asserts (`.py` local ledger idioms), `fail[]`/`Exit[1]` on failure
(`.wl` already has `fail[]`). **KEEP the `.wl`'s already-independent native route** (§1b) — re-target it to assert its OWN
derived `u2z=1/9`/`u4z=4/81`/`v5z=1/27`/`χ_Q=+1`/incoming `−1`/standing static slot. **Dual-engine agreement =
transcript-level** (stage014/016/017 pattern). **Zero file I/O.** Arity discipline (standing).

## §4 Consumed / exported

- **Consumes (PROVENANCE only — cite, do NOT re-derive, do NOT build a theatrical dual-site; §1c):**
  - **017's ℓ=2 port kernel** (the grouped `M₂`, angular `K₂=K̃+6·T̃_Ω`, support scalars `B̃/Z̃`, D-lanes) — the INTERIOR
    wall mode whose EXTERIOR radiative response 018 computes. Cite `stage017`. ⚠ The literal `N_n/D_n` consumption is
    **019's** (`build_prefactor`), NOT 018's — 018's fingerprint is self-contained.
  - **009/010's bulk Helmholtz mode** — the exterior ℓ=2 outgoing solution's bulk companion. Cite `stage009`/`stage010`.
    018 reconstructs `h₂⁽¹⁾` self-contained from `j₂`/`y₂` (no stored-artifact read).
  - **`c_s`** (the density sound speed) — cite `stage005` R1 (`c_s²=5Kρ⁴/m`) as the PROVENANCE of the units symbol; NOT a
    consumed value (§1c). ⚠ Distinct from the frozen-wall `c_S` (011–017).
- **Exports (→ 019/020/022/027):** the outgoing ℓ=2 fingerprint = the radiative slot `v₅=a⁵/27c_s⁵` (→ the `27` that
  **020's** `54/5=2·27/5` partition rides), the reactive coefficients `u₂=a²/9c_s²`, `u₄=4a⁴/81c_s⁴` (→ **019**'s
  prefactor context), `χ_Q=+1` (→ **020**'s `Γ₅=2χ_Q·G/5c⁵` equivalence), and the outgoing/incoming/standing branch
  classification. Per the cross-stage flow (`part2_gravity_atomic_split.md` L107): "018–021 export the Λ₂ fingerprint +
  χ_Q=1 + the 54/5 partition → 022 (non-regression) + 027 (closure)." Cite the exact export contract at note-authoring.

## §5 Teeth candidates (018-specific, per-tooth ablation MANDATORY)

1. **⭐ The fingerprint-value teeth (`u₂ᶻ=1/9`, `u₄ᶻ=4/81`, `v₅ᶻ=1/27`).** The derived `out["u2z"]` etc. (L129–131,
   `.wl` L40–42) are genuinely series-computed from the real `h₂⁽¹⁾`; the L164–167 / `.wl` L66–68 check is derive-vs-typed
   (independent). ⚠ **Per-tooth ablation MUST confirm these are NOT stamped** (the §4 firewall: forbidden to hardcode
   `1/9,4/81,1/27` and "check" against them, or an unconstrained "outgoing solve") — mutate the DERIVATION (e.g. corrupt
   the `j₂`/`y₂` expression, or use the wrong Hankel order) → the derived coefficient changes → the match fails. This is
   the central EARNED tooth; the derivation-from-`−3/Λ₂` must be EMITTED (directive §4 firewall), not just the compare.
2. **⭐ The `χ_Q=+1`-outgoing / `−1`-incoming sign tooth (`3a_wrong_bc`, `FAIL_FINGERPRINT`).** The incoming `h₂⁽²⁾` flips
   ONLY the imaginary `v₅` sign (even coeffs unchanged → `χ_Q=−1`); standing `j₂` gives `Λ_stand(0)=+2`, `Ŷ_stand(0)=−3/2`,
   no imaginary term. Per-tooth: mutate the branch (outgoing→incoming) and confirm the tuple changes in the predicted way
   (χ_Q sign flip) AND that a typed `χ_Q=1` cannot survive (the sign is computed from `j₂±i·y₂`). ⚠ Confirm `χ_Q` is
   COMPUTED (L168–169), never string-typed.
3. **⭐ The passivity tooth (`3b_imposed_dissipation`, `FAIL_NOT_OUTGOING`).** Inserting a phenomenological damping term
   instead of the outgoing BC ⇒ `FAIL_NOT_OUTGOING`; removing the genuine outgoing BC removes the imaginary `v₅`.
   Per-tooth: the radiating term must come from the BC (the `dtn_hankel1` branch), not an inserted sink — the
   `outgoing_ablation` (L758) is a DYNAMIC re-run (confirm `with_mutation`≠`without_mutation`, not a constant
   `self_ablation` — the v1 trip-up).
4. **The dim leg tooth.** `[u₂]=[a²/c_s²]=T²`, `[v₅]=T⁵` — corrupt a coefficient's `a`/`c_s` power → the dim leg fails.
   Per-tooth: reads the REAL physical coefficient, not a typed tuple (the anti-vacuous mandate).
5. **Provenance-cite integrity (light).** 017 port kernel / 009-010 bulk mode / R1 `c_s` cited as PROVENANCE (§4) — a
   citation guard, NOT a theatrical dual-site (018 does not consume their values); + `.wl` arity scan.

⚠ **NOT 018 (do not rebuild as 018 teeth — 019/020/021 own these):** `3g_wrong_prefactor_object` (019); `3c_wrong_scaling`
(020 — the typed-target a⁻⁵ scaling, a weak bookkeeping tooth per sub-agent §8, to be STRENGTHENED at 020), `3e_equivalence_break`
(020 — the PN `2G/5` bridge), `3f_partition_mislabel` (020); `3d`/`3d′` dimensional-break (021 — the μ̂₀-free `[P₀^phys]`
gate + corrupt-`N₀`). ⚠ **The `54/5=2·27/5` STRING label (`.py` L623) is 020's** — it is currently a typed string, NOT a
SymPy-verified identity (sub-agent §8); 020 must make the decomposition COMPUTED. 018 only owns the `27=1/v₅ᶻ` half.

## §6 Register expectation — ⭐ THE KEY 018 QUESTION (likely ZERO new counted knobs; CONFIRM)

Per headline #3 + the split: **018 is the EARNED fingerprint slice (exterior-wave structure); the CALIBRATED label +
the `54/5` magnitude are 020's.** So the honest pre-read (⚠ CONFIRM at the register step + Codex-verify against the
scripts):

- **⭐ 018 likely adds ZERO new counted knobs** (like 016 / 011/012/014 — an EARNED/structural slice). The outgoing DtN
  Hankel fingerprint is pure exterior-wave math (dimensionless z-space rationals from `−3/Λ₂` + the χ_Q sign) — it
  introduces NO calibration. `c_s` is R1-DERIVED (cited PROVENANCE, a units carrier — NOT a new knob, §1c); `a` is the
  `CONV` pin (R2-family); the port scalars `N_n/D_n` are 019's deferred Gate-6 branch data (`deferred_branch_data`, not
  018's).
- **⭐ Likely a new STRUCTURAL edge (call it R37 — confirm the next free number at registration):** the outgoing-DtN
  ℓ=2 Hankel-fingerprint provenance — the five-slot exterior radiative signature `u₂=a²/9c_s², u₄=4a⁴/81c_s⁴, v₅=a⁵/27c_s⁵`
  (the `27` COMPUTED from `v₅ᶻ`) + `χ_Q=+1` outgoing / `−1` incoming + the standing-branch contrast — **discharges
  NOTHING** (earned exterior-wave structure, not a reduction of a debt; like R34 was for 016's covariance theorem).
- **⚠ Flag the χ_Q Part-VII reconciliation debt — do NOT reconcile here.** pathA_33 gives `χ_Q=+1` (the outgoing-DtN
  Hankel context); pathA_22b gave `χ_Q≈0.712` (the older minimal-combination context) — same name, DIFFERENT computations
  (blueprint §8 tracked item, L173/L307). 018 lands `χ_Q=+1` in ITS context and records the reconciliation as a Part-VII
  open item; it does NOT silently merge the two.
- **Cited provenance (NOT re-counted):** `c_s` (R1, stage005), 017's port kernel (`{T_Ω,β₂}` counted at 017, `{B̃,Z̃}`
  tracked/downstream-pinned), 009/010's bulk mode. ⚠ **Do NOT let 018 count the port-kernel `N_n/D_n` scalars** (019's
  deferred branch data) or re-count 017's `{T_Ω,β₂}`.
- **Control/tracked-not-counted:** the sample subs `{a:3, c_s:2}` (numeric evaluation controls, like 014's controls).

⚠ **Do NOT let 018 silently count `c_s` or the port scalars.** Resolve `c_s` as R1-DERIVED-provenance (units carrier),
the port scalars as 019's deferred branch data, and Codex-verify (the register verify is the gate that catches an
over-count that would falsely inflate — or a mislabel that would falsely shrink — the irreducible codimension count).

## Verdict tokens + honest scope

018 carries the **outgoing ℓ=2 DtN Hankel-fingerprint component (1/4) of `QUAD_CALIBRATED` — the EARNED-FIRST slice**:
the outgoing series `u₂=a²/9c_s², u₄=4a⁴/81c_s⁴, v₅=a⁵/27c_s⁵` DERIVED from `Ŷ₂ᵒᵘᵗ=−3/Λ₂ᵒᵘᵗ` (the `27` COMPUTED from
`v₅ᶻ=1/27`), the sign `χ_Q=+1` outgoing / `−1` incoming (COMPUTED from `j₂±i·y₂`), the standing-branch contrast, and the
passivity. EARNED = the exterior radiative signature is DERIVED (not calibrated); the prefactor algebra = 019, the
`54/5=2·27/5` provenance partition + the CALIBRATED label + `G=GENUINE_BLOCKED` = 020, the μ̂₀-free dim closure = 021.
SELF-CONTAINED exterior-wave algebra — cites 017's ℓ=2 port kernel + 009/010's bulk mode as PROVENANCE (the physical
two-port setup; literal `N_n/D_n` consumption is 019's), and `c_s` (R1, stage005) as the units-symbol provenance (NOT a
consumed value). Caveats: `c_s` is a units carrier not a live value (the earned rationals are c_s-free); the `54/5`
magnitude + `G` are CALIBRATED (020, `G=GENUINE_BLOCKED`); the actual-branch `a`-scaling + the numerical `N_n/D_n` port
scalars remain Gate-6 (report :49); the `χ_Q` reconciliation with pathA_22b's `≈0.712` is a tracked Part-VII item (not
merged here). ⭐ The pathA_33 trip-ups: the `27` stays COMPUTED (`v₅ᶻ=1/27` series-derived), `χ_Q` COMPUTED (not stamped),
the fingerprint slice free of any `μ̂₀` back-solve (the back-solved-`μ̂₀` + constant-`self_ablation` v1 rig is 021's/the
probe-engine's concern — must not return there).

## Process (unchanged, calibrated — the per-stage pipeline)

Author the II-G4a reshape directive (§1 the clean 018 slice / 4-way cut + §1b the native-`.wl` KEEP + §1c the
provenance-only c_s + port-kernel framing + §2 faithful cover + §3 bridge-strip incl. sever-YAML + transcript-level
agreement + §5 the fingerprint-value / χ_Q-sign / passivity teeth with per-tooth ablation + §6 the ZERO-new-knobs +
R37-edge + χ_Q-Part-VII-flag register question) → **Codex xhigh design-review** → fold to `DIRECTIVE_CLEAN` → **⭐ final
Grok-4.5 headless compute-verify pass** (assess + independently verify each catch — Grok compute-verifies the
spherical-Hankel series / the χ_Q sign / the c_s-free-ness of the rationals; it caught the 016 volume-vs-line convention
mismatch, so watch the `c_s`/`[c_s]` dim accounting + the fingerprint-value genuineness here) → fold → **Codex confirm-pass
on the folds** → **pre-exec USER GATE** → Codex builds the two scripts (`--sandbox danger-full-access`, background, xhigh)
→ dual-engine both exit 0 (repo root + foreign CWD) → arbiter re-run → full tri-review (fidelity + adversarial-with-
**per-tooth ablation**; ⭐ hunt the fingerprint-value stamped-vs-derived genuineness + the χ_Q computed-ness + a mirror-`.wl`
+ any vacuous able-to-fail) → remediate → fresh-agent re-verify → bump counts 17→18 → parameter register (⭐ confirm ZERO
new 018 knobs + `c_s` R1-provenance + R37 edge + the χ_Q Part-VII flag) + Codex-verify → note/card/`\input{stages/stage_018}`
+ registration → rebuild PDF → commit + docs/memory sync. Orchestrator authors notes/cards/LaTeX/registration; Codex codes.
Target stem: `ledger_stage018_dtn_hankel_fingerprint` (confirm slug at authoring).
