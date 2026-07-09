# II-G4b (ledger_stage019) source map — pathA_33 prefactor algebra `P(ω)=D₀N/D^cons²` (QUAD_CALIBRATED 2/4)

> Running-start prep captured 2026-07-09 (post stage018, before authoring the II-G4b reshape directive) so the directive
> can be written without re-discovery. **All line refs below VERIFIED against the current sources 2026-07-09**
> (`pathA_33_quadrupole_normalization_sympy.py` = 1351 lines; `pathA_33_quadrupole_normalization.wl` = 256 lines;
> `reports/pathA_33_quadrupole_normalization.md` = 49 lines; `directives/pathA_33_quadrupole_normalization.md` = 534 lines).
> Companion: `part2_gravity_atomic_split.md` (row 019 = the pathA_33 4-way `QUAD_CALIBRATED` split, L39; the pathA_33
> trip-ups L87–89; the cross-stage flows L106–108) and the **stage018 source map + reshape directive** (the pathA_33
> EARNED-FIRST exemplars — 019 continues the SAME 4-way split, is the SECOND leg, and reuses 018's provenance-only
> consumption discipline + KEEP-native-`.wl`-sever-only-YAML pattern). Build-order id **019**, Part II.
> Source top-line: **`QUAD_CALIBRATED`** — 019 lands the 2/4 EARNED component (the squared-denominator prefactor
> algebra); 018 (fingerprint, DONE `4872e8b7`) is 1/4; 020 (54/5 partition + the CALIBRATED label) is 3/4; 021 (μ̂₀-free
> dim closure) is 4/4.

## ⭐ The FIVE headline differences from 018 (READ FIRST)

1. **⚠ 019 is the SECOND leg of the SAME 4-way pathA_33 split (018/019/020/021) — a MIDDLE PARTIAL leg** (neither the
   EARNED-FIRST 018 nor the COMPLETING 021). Like 018 it lands `QUAD_CALIBRATED` as a **PARTIAL** (the CALIBRATED landing
   is 020's). 019 = the **squared-denominator prefactor algebra** `P(ω)=D₀·N(ω)/D^cons(ω)²`: the Taylor coefficients
   `P₀=N₀/D₀`, `P₂=(D₀N₂−2D₂N₀)/D₀²`, `P₄=(D₀²N₄−2D₀(D₂N₂+D₄N₀)+3D₂²N₀)/D₀³` SERIES-EXTRACTED from the actual series of
   `D₀·N/D^cons²`, **plus the N/D self-check** (plain `N/D` gives `P₂=(D₀N₂−D₂N₀)/D₀²` — missing the factor of 2 on
   `D₂N₀`; the squared denominator is the only object that produces the `−2`). **The cut is fingerprint (018, EARNED)
   vs prefactor algebra (019, EARNED) vs magnitude-provenance (020, the CALIBRATED landing) vs dim-closure (021).**

2. **⭐⭐ THE PORT-KERNEL CONSUMPTION RESOLVES TO PROVENANCE, NOT a checkable dual-site (the KEY 019 question — resolved).**
   This is where the literal `N_n/D_n` port-kernel consumption that 018 DEFERRED nominally "happens" — but the honest
   source reading is that **it happens at the PROVENANCE level, not as a value-consumption dual-site**:
   - `build_prefactor` (L212–273) proves the squared-denominator algebra over the **ABSTRACT, port-AGNOSTIC** symbols
     `D0, D2, D4, N0, N2, N4` (declared L47, generic `nonzero=True`). The algebra `P(ω)=D₀N/D^cons²` holds for ANY nonzero
     `D0..N4` — it never references the concrete port forms.
   - `build_port_moments` (L190–209) — the concrete ℓ=2 port N-moments `N_A0_r=P_port²/Δ²`, `N_A2_r=2P_port(P_port·S−Δg_W)/Δ³`,
     `N_A4_r=(…)/Δ⁴` (with `P_port=Ω_U²g_W+R·g_U`; symbols `Ω_U,g_W,R,g_U,Δ,S` L51–53) — is **EMITTED BUT NEVER CHECKED.**
     ⭐ **VERIFIED (2026-07-09 grep):** its output appears ONLY at `ctx["port_moments"]=build_port_moments()` (L1006) and
     `"outgoing_port_moments": stringify_expr_tree(ctx["port_moments"])` (L1150, stringified into the payload). It is
     referenced by NO residual, NO match, NO verdict gate; `build_prefactor` does NOT read it. So the port N-moments are
     **deferred Gate-6 branch data, carried symbolically, asserted against nothing** — exactly the running-start key-pin-1
     hypothesis, CONFIRMED. (The directive's reduction certificate says the same, L357–362: the symbolic `D_n` + the
     symbolic port scalars `(Ω_U,Ω_W,g_U,g_W,R,Δ,S)/N_n` are **FROZEN-INPUT**; the numerical `(D_n,N_n)` are **DEFERRED**
     Gate 6; the `P₀/P₂/P₄` algebra is **COMPUTED**.)
   - **⭐ RESOLUTION (like 018, like 016's convention-boundary lesson): PROVENANCE-only.** 017's exported ℓ=2 port kernel
     (the D-lanes `D0/D2/D4` + the support/Maxwell scalars `{B̃,Z̃}`) and `build_port_moments`' concrete N-moments are
     cited as the PROVENANCE of what the abstract `D0..N4` physically ARE (the ℓ=2 port), NOT consumed as values. There is
     **NO checkable numeric cross-stage relation** — the algebra is port-agnostic; the concrete port instance is deferred.
     ⚠ **Do NOT manufacture a theatrical dual-site on the port moments** (they are unused/deferred → a guard on them is a
     vacuous tooth, the exact class per-tooth-ablation catches; the 017 lesson banked at `part2` L331–334 + the 018 §1c
     discipline). This is the SAME provenance-only landing as 018, one leg later.
   - ⚠ **Latent structural observation (do NOT build a check on it): `build_port_moments`' `n0/n2/n4` have the SAME
     squared-denominator SHAPE as `P₀/P₂/P₄`** (`n0=P²/Δ²` ↔ `P₀=N₀/D₀`; the factor-of-2 in `n2` ↔ the `−2D₂N₀` in `P₂`;
     the factor-of-3 in `n4` ↔ the `+3D₂²N₀` in `P₄`) — i.e. `build_port_moments` IS the concrete port instance of the
     same abstract algebra `build_prefactor` proves. The SOURCE does NOT wire this into a check (deferred branch data), and
     the faithful reshape must NOT add an unrequested cross-check (new-physics/rig risk, [[feedback-claude-reviews-codex-codes]]).
     Note it as provenance narrative only; 019 asserts NOTHING against the concrete moments.

3. **⭐ 019 is c_s-FREE, a-FREE, G-FREE, μ̂₀-FREE — pure abstract algebra over `{ω, D0..N4}`.** UNLIKE 018 (which carried
   `c_s`/`a` as units symbols), 019's variables are ONLY `omega` (the expansion variable) + the six abstract port scalars
   `D0,D2,D4,N0,N2,N4`. **No `c_s`, no `a`, no `G`, no `μ̂₀`, no `z`, no `j₂`/`y₂`.** (The `c_s`/`a` units carriers are
   018's fingerprint; `G`/`54/5` are 020's; `μ̂₀`/`N₀`-dim/`c_s`-restore are 021's.) So 019 has **no new units symbol** and
   **no dimensional leg** (the dim closure `[P₀^phys]=1` via `(c_s/a)²·(N₀/D₀)` is 021's, source L512 `build_dimensions`).
   019's only dimensional content is the abstract algebraic structure of the coefficients — no units.

4. **⭐ The central EARNED tooth = the SQUARED-DENOMINATOR object `P(ω)=D₀·N/D^cons²` (NOT plain `N/D`).** The prefactor
   coefficients must be **series-EXTRACTED** — `coeff(ω,n)` off the actual `series_no_o(D₀·N/D^cons², ω, 6)` (source L217,
   L220–222; `.wl` `serW`/`Coefficient` L77–81), NOT typed. The `expected` P0/P2/P4 (L224–233) are the TYPED reference the
   derived `coeffs` are checked against (derive-then-check, genuinely independent — a corrupted object changes the
   series-extracted coeff and the residual fires). The **N/D self-check** (`plain_equals_correct_P2`, L263, `.wl` L92) is
   the sharp discriminator: `plain["P2"]=(D₀N₂−D₂N₀)/D₀²` ≠ `correct expected P2=(D₀N₂−2D₂N₀)/D₀²`; the difference
   `plain−correct = D₂N₀/D₀²` (nonzero, L264) → `plain_equals_correct_P2=False` → probe `3g` fires `FAIL_PREFACTOR_ALGEBRA`.
   This is the §2.2/§3g firewall discipline (cf. 018's fingerprint-derivation genuineness): the coefficients EMERGE from
   the series of the *right* object; the plain-`N/D` object is provably WRONG (the missing factor of 2).

5. **⭐⭐ The pathA_33 trip-ups live at 020/021, and 019 owns the "keep-it-clean" obligations** (`part2` L87–89;
   directive L23–26): (i) **the `27` stays 018's** (`v₅ᶻ=1/27` series-derived) — 019 does NOT touch the `27`/`54/5`/`G`
   magnitude (020's); (ii) **the v1 rig (back-solved free-carrier `μ̂₀` + constant `self_ablation`) must NOT return** —
   019's prefactor slice has ZERO `μ̂₀` refs (`μ̂₀` enters ONLY `build_dimensions`, 021's), and its `3g` self-ablation is a
   **DYNAMIC** re-run (source `prefactor_ablation = ablation({"prefactor_ok": False}, …)` L762 → `ablation` L728–755 which
   RE-RUNS the verdict gate logic, `rerun_gate_logic:True`, `with_mutation`≠`without_mutation` computed), NOT a constant.
   Keep 019 free of any `μ̂₀` back-solve and any `54/5`/`G` magnitude.

## §1 The 019 slice (`.py` line ranges) — the CLEAN CUTS (all VERIFIED)

The whole computation is `build_*` helpers assembled by `build_counterfactuals` (L703–902) + `main` (L1321–1347).
**019 owns the prefactor algebra + the port N-moments (deferred, provenance) + the `3g` prefactor probe; 018/020/021 own
the fingerprint/partition/dimension blocks.** The 019-owned cuts:

- **⭐ `build_prefactor()` L212–273** — THE EARNED core:
  - `Nomega = N0 + N2·ω² + N4·ω⁴` L213 (the port numerator); `Dcons = D0 + D2·ω² + D4·ω⁴` L214 (the D-lanes /
    "consumption" denominator).
  - **`correct_obj = D0·Nomega/Dcons²` L215** (the SQUARED-DENOMINATOR object) ; `plain_obj = Nomega/Dcons` L216 (the WRONG
    object, the N/D self-check).
  - `correct_series = series_no_o(correct_obj, ω, 6)` L217 ; `plain_series = series_no_o(plain_obj, ω, 6)` L218.
  - **`coeffs = {P0: series.coeff(ω,0), P2: series.coeff(ω,2), P4: series.coeff(ω,4)}` L219–223** (SERIES-EXTRACTED off
    `correct_series` — the genuine derive; `compact`-wrapped).
  - **`expected = {P0:N0/D0, P2:(D0·N2−2D2·N0)/D0², P4:(D0²N4−2D0(D2N2+D4N0)+3D2²N0)/D0³}` L224–233** (the TYPED reference).
  - `residuals = {coeffs[n] − expected[n]}` L234 ; `matches = {bool_from_residual(residual)}` L235 (derive-vs-typed).
  - `plain = {P0/P2/P4 from plain_series}` L236–240.
  - **`self_check` L260–265** — `correct_P2_D2N0_term="-2*D2*N0"`, `plain_P2_D2N0_term="-D2*N0"` (labels), the LIVE
    `plain_equals_correct_P2 = bool_from_residual(plain["P2"]−expected["P2"])` L263 (=False → the factor-of-2 gap), and
    `difference_plain_minus_correct_P2 = compact(plain["P2"]−expected["P2"])` L264 (=`D2·N0/D0²`, the computed gap).
  - `sample_subs = {D0:19, D2:23, D4:29, N0:11, N2:13, N4:17}` L241–248 (distinct primes) → `sample_values` L266–272 (a
    light NUMERIC cross-check of the symbolic coeffs; `numeric` L93).
  - `ok = all(matches.values())` L258.
- **`build_port_moments()` L190–209** — the concrete ℓ=2 port N-moments (`P_port`, `N_A0_r`/`N_A2_r`/`N_A4_r`,
  `isotropic_branch` string L208). ⭐ **DEFERRED Gate-6 branch data — EMITTED BUT NEVER CHECKED** (§headline-2). 019 carries
  it as a PROVENANCE/export narrative block (the concrete port instance the abstract `N0..N4` stand for), asserting
  NOTHING against it. **It is NOT a tooth** (a labeled non-check export, not a can't-fail check).
- **probe `3g_wrong_prefactor_object` L866–875** — reads `prefactor["plain_object"]`, `plain_N_over_D["P2"]`,
  `expected["P2"]`, and `plain_equals_correct_P2`; verdict `FAIL_PREFACTOR_ALGEBRA` iff `not plain_equals_correct_P2`
  (else `NO_FAIL`); `self_ablation = prefactor_ablation` (the DYNAMIC re-run, L762).
- **`prefactor_ablation` L762 + the `ablation` helper L728–755** — the DYNAMIC self-ablation: `ablation({"prefactor_ok":
  False}, expected_fail="FAIL_PREFACTOR_ALGEBRA")` RE-RUNS `base_verdict(mutated_gates)` vs `base_verdict(baseline_gates)`,
  computes `fail_suppressed` (`rerun_gate_logic:True`, `with_mutation`≠`without_mutation`). ⚠ **Note the scope caveat (like
  018's 3a/3b):** `base_verdict` (L676–694) re-runs the JOINT gate set incl. 018/020/021's gates, which are NOT built here.
  Re-scope the `3g` `with_mutation`/`without_mutation` to an **019-local verdict over 019's gate(s) only** (the prefactor
  match ∧ the N/D self-check), do NOT pull the joint `base_verdict`. (This is the same re-scoping the 018 directive applied
  to 3a/3b, directive §2e L305–309.)
- **Shared helpers 019 uses (NOT cut boundaries):** `compact` L60, `series_no_o` L97, `bool_from_residual` L89, `numeric`
  L93. Symbol decls L45–57 (019 uses `omega, D0, D2, D4, N0, N2, N4`; NOT `z, a, c_s, c, G, mu_hat0, mtilde0, chi_Q`).

- **⭐ CLEAN CUT — 019 owns L190–209 (port moments, provenance) + L212–273 (prefactor algebra, EARNED) + probe `3g`
  L866–875 + `prefactor_ablation` L762. It touches NONE of `build_fingerprint`/`build_scaling`/`build_equivalence`/
  `build_partition`/`build_dimensions`. Do NOT pull in 018/020/021 territory:**
  - **018 (DONE):** `spherical_j2`/`spherical_y2` L101–106 + `dtn_branch` L109–149 + `build_fingerprint` L152–187 +
    `passivity_from_source` L666–673 + probes `3a` L774–797 / `3b` L799–812. (019 CITES 018's fingerprint context as
    provenance only — §4.)
  - **020:** `a_power`/`build_scaling` L276–294 + `build_equivalence` L510–535 + the provenance machinery L538–663
    (`classify_provenance` L575–585, `build_partition` L601–628 incl. the STRING `decomposition_54_over_5` L622–627, the
    g-invariance diagnostic L635–663) + probes `3c` L813–818, `3e` L837–841, `3f` L846–865.
  - **021:** the dim engine `DimError`/`dim_*` L297–376 + `SOURCED_N0_DIM`/`SOURCED_D0_DIM` L379–380 + `build_dimensions`
    L387–507 (`P0_raw=N0/D0` L399, `P0_physical=(c_s/a)²·(N0/D0)` L512, the μ̂₀-free `[P0_phys]` gate + the μ̂₀ diagnostic +
    drop-norm / corrupt-`N0` probes) + probes `3d` L820–826, `3d′` L828–836. ⚠ **`N0`/`D0` are SHARED symbols** (declared
    L47) used by BOTH `build_prefactor` (019, abstractly) AND `build_dimensions` (021, with dimensions). 019 uses them
    algebraically (no units); the μ̂₀-free dim closure that assigns `[N0]/[D0]` is 021's. Keep 019 units-free.

## §1b The `.wl` 019 slice (VERIFIED) — the independent Wolfram route (KEEP native, sever ONLY the YAML)

⭐ **The pathA_33 `.wl` is ALREADY a genuinely independent engine** (native `serW`/`Coefficient`/`FullSimplify` on the
native prefactor object; its OWN algebra; it does **NOT** `Get`/`Import` the `.py`'s expressions). The ONLY bridge is the
scratch-YAML `Export` at L253. The reshape KEEPS the native route and severs ONLY that YAML handoff (§3). 019's `.wl`
slice = **the prefactor block L73–92**:
- `Nomega=N0+N2 ω²+N4 ω⁴` L73, `Dcons=D0+D2 ω²+D4 ω⁴` L74, `prefObj=D0 Nomega/Dcons²` L75, `plainObj=Nomega/Dcons` L76.
- `prefSeries=serW[prefObj,6]` L77, `plainSeries=serW[plainObj,6]` L78 (native `serW=Collect[Normal[Series[…]],ω,FullSimplify]`
  L28); `p0=Coefficient[prefSeries,ω,0]` L79 / `p2` L80 / `p4` L81 (SERIES-EXTRACTED, native).
- `expectedP0=N0/D0` L82, `expectedP2=(D0 N2−2 D2 N0)/D0²` L83, `expectedP4=(D0²N4−2 D0(D2 N2+D4 N0)+3 D2²N0)/D0³` L84
  (TYPED reference).
- `resP0/resP2/resP4 = FullSimplify[p_n − expected_n]` L85–87 ; `p0Match/p2Match/p4Match = TrueQ[res==0]` L88–90
  (derive-vs-typed).
- `plainP2=Coefficient[plainSeries,ω,2]` L91 ; `plainEqualsCorrectP2 = TrueQ[FullSimplify[plainP2==expectedP2]]` L92 (the
  N/D self-check — comes out False, the factor-of-2 gap).
- `sampleRules` L94–99 (`D0→19, D2→23, D4→29, N0→11, N2→13, N4→17`, + the 018/020 carriers `a,cs,c,G` which 019 ignores)
  → `evalSample` L100 (the light numeric cross-check).
- ⚠ **The `.wl` has NO port-moments block** (grep-confirmed: no `P_port`/`N_A0`/`OmegaU`/`Delta`-port in the `.wl`; the
  `Delta` at L122 is an unrelated dim-error message). So `build_port_moments` is `.py`-only deferred branch data; the `.wl`
  019 slice is PURELY the abstract prefactor algebra. Dual-engine agreement is on the prefactor coefficients + the N/D
  self-check.
- **⚠ EXCLUDE from 019 (018/021 territory in the `.wl`):** the fingerprint block L30–71 (→ 018, DONE); the dimensional
  block `zeroDim`/`dimOf`/`rawDims`/`P0Physical`/`Gamma5`/`targetRHS`/`muDim`/drop-norm+corrupt-N0 L101–174 (→ 021); the
  combined `headlineOk` L176–180 (rebuild an **019-scoped** headline over `p0Match∧p2Match∧p4Match∧(not plainEqualsCorrectP2)`
  only); the YAML payload L182–255.
- **⚠ The bridge to sever (§3):** `scratchDir`/`yamlOut` setup L19–22 + the YAML-line assembly `lines={…}` L182–251 + the
  `Export[yamlOut,…]` L253 + the `headlineOk` guard L254. SEVER: no scratch YAML, print-only + `fail[]`/`Exit[1]` on failure
  (the `.wl` already has `fail[]` L5). **Dual-engine agreement = transcript-level** (both engines print the same series-
  extracted `P0=N0/D0`, `P2=(D0N2−2D2N0)/D0²`, `P4=…`, the plain-`N/D` `P2=(D0N2−D2N0)/D0²`, and `plainEqualsCorrectP2=False`
  — the stage018 transcript pattern). **Zero file I/O.** Arity discipline (standing — def/call scan + unevaluated-leakage
  transcript scan; the stage007 lesson).

## §1c The consumption resolution (READ — the honest PROVENANCE framing; the KEY 019 discipline)

⭐ **The port-kernel consumption RESOLVES to provenance-only, exactly like 018 (grep-confirmed §headline-2):**

- **017's ℓ=2 port kernel (the D-lanes `D0/D2/D4` + support/Maxwell scalars `{B̃,Z̃}`) + `build_port_moments`' concrete
  N-moments:** cite as PROVENANCE — the abstract `D0..N4` in `build_prefactor` ARE 017's exported port kernel scalars,
  carried as generic symbols; the algebra is PORT-AGNOSTIC (holds for any nonzero `D0..N4`); the concrete port instance
  (`build_port_moments`) is DEFERRED Gate-6 branch data, asserted-against-nothing. **There is NO checkable numeric cross-
  stage relation** (no value is consumed). ⚠ Build NO guard/dual-site on the port kernel or the port moments (a guard on
  deferred/unused data is a vacuous tooth). This is the SAME provenance-only landing as 018 (§1c of the 018 source map),
  one leg later — and consistent with the 016 convention-boundary lesson ("cite as PROVENANCE across a boundary; do not
  manufacture a theatrical dual-site").
- **018's fingerprint context:** cite `stage018` as PROVENANCE — 019's prefactor `P(ω)` is the exterior/interior MATCH the
  018 fingerprint's `v₅=a⁵/27c_s⁵` slot feeds; but 019's algebra is over the abstract port scalars, NOT the fingerprint
  rationals. 019 does NOT re-derive or re-check 018's `1/9,4/81,1/27`/`χ_Q` (those are 018's, DONE). Provenance narrative
  only.
- ⚠ **NO `c_s`/`a`/`G`/`μ̂₀` in 019** (§headline-3): the abstract prefactor algebra is units-free. The `c_s`/`a` units
  restoration + the `[P₀^phys]=1` dim closure are 021's (`build_dimensions`); the `54/5`/`G` magnitude is 020's. 019 has
  NO dimensional leg and introduces NO units symbol.

## §2 The 019 claim-set (derive + assert; report/directive quotes)

- **The squared-denominator prefactor algebra (EARNED — the headline).** Series-expanding the SQUARED-DENOMINATOR object
  `P(ω)=D₀·N(ω)/D^cons(ω)²`, `N(ω)=N₀+N₂ω²+N₄ω⁴`, `D^cons(ω)=D₀+D₂ω²+D₄ω⁴` about `ω=0` yields the Taylor coefficients
  `P₀=N₀/D₀`, `P₂=(D₀N₂−2D₂N₀)/D₀²`, `P₄=(D₀²N₄−2D₀(D₂N₂+D₄N₀)+3D₂²N₀)/D₀³` (report :15–17; directive §2.2 L189–193). ⭐
  **The coefficients are SERIES-EXTRACTED** (`coeff(ω,n)` off the actual `series_no_o(D₀N/D^cons²)`, source L217–222; `.wl`
  `serW`/`Coefficient` L77–81), checked against the TYPED `expected` (derive-then-check; a corrupt object changes the
  extracted coeff and the residual fires). The `−2D₂N₀` term is the SIGNATURE of the squared denominator (a plain `N/D`
  never produces it).
- **The N/D self-check (EARNED — the sharp discriminator; the `3g` tooth).** Plain `N/D` gives `P₂=(D₀N₂−D₂N₀)/D₀²` —
  missing the factor of 2 on `D₂N₀`; the computed difference `plain−correct = D₂N₀/D₀² ≠ 0` (source L264; report :18). So
  `plain_equals_correct_P2=False` (L263) → probe `3g` fires `FAIL_PREFACTOR_ALGEBRA` (report :41), and the correct
  `D₀N/D^cons²` object does NOT fire (self-ablation). ⭐ This makes the squared-denominator object PROVABLY the right one
  (directive §3g L338–341, §4 firewall item (4) L354).
- **The port N-moments (PROVENANCE / deferred branch data — asserted-against-nothing).** `build_port_moments` (L190–209)
  emits the concrete ℓ=2 port N-moments (`N_A0_r=P_port²/Δ²` etc.; `isotropic_branch` = "N_20,n=N_21,n=N_22,n=N_n and
  D_20,n=D_21,n=D_22,n=D_n carried symbolically", L208) — the concrete port instance the abstract `N₀..N₄` stand for. ⭐
  Carried as a PROVENANCE/export narrative block, DEFERRED Gate-6 branch data (report :49 "the numerical branch data
  `(D_n,N_n)`, port scalars … remain outside Gate 4"); 019 asserts NOTHING against it. NOT a tooth.
- **The 019-scoped landing (PARTIAL component).** Land at the 019 component: `QUAD_CALIBRATED (2/4) = the squared-
  denominator prefactor algebra EARNED (P₀=N₀/D₀, P₂=(D₀N₂−2D₂N₀)/D₀², P₄=…, SERIES-EXTRACTED from D₀N/D^cons² + the N/D
  factor-of-2 self-check), with the fingerprint = 018 (DONE), the 54/5 provenance partition + the CALIBRATED label = 020,
  and the μ̂₀-free dim closure = 021.` Do NOT print the joint as complete (that is 020's CALIBRATED landing) and do NOT
  re-present 018/020/021's accounting. ⭐ Continues 018's EARNED PARTIAL-landing pattern.

## §3 Reshape cost (the bridge to sever) — cross-script scratch-YAML family, KEEP the native `.wl`

Same family as pathA_30–34 / 018 (the cross-script runtime-YAML reshape). No argparse. The `.py`'s `build_*` helpers are
pure/self-contained, but `main` L1321–1347 writes `SYM_YAML` (L1324), reads `MMA_YAML` (L1326), writes `RESULTS_YAML`/
`REPORT_MD`/`FEED_NOTE` (L1333–1335); `compare_engines` L1024–1099 cross-checks; the `.wl` writes its scratch YAML via
`Export` L253. **Reshape = sever ALL file I/O both directions:** drop `main`'s YAML/report/feed writers + `yaml_read`/
`yaml_write` (L72–86) + `compare_engines`/`engine_summary`/`build_final_payload`/`build_report`/`build_feed_note` (`.py`);
drop the `Export` + the YAML-line assembly L182–255 (`.wl`). Each engine → standalone: print-only, `expect_zero`/
`bool_from_residual`-style asserts (`.py` local ledger idioms — `banner`/`subbanner`/`_record_pass`/`_record_fail`/
`expect_zero`/`expect_bool`/`expect_fail`, a `Verdict labels:` block, tallies, `OVERALL PASS`/nonzero exit), `fail[]`/
`Exit[1]` on failure (`.wl`). **KEEP the `.wl`'s already-independent native route** (§1b) — re-target it to assert its OWN
series-extracted `P0/P2/P4` + the N/D self-check. **Dual-engine agreement = transcript-level** (stage018 pattern). **Zero
file I/O.** Arity discipline (standing).

## §4 Consumed / exported

- **Consumes (PROVENANCE only — cite, do NOT re-derive, do NOT build a theatrical dual-site; §1c):**
  - **017's ℓ=2 port kernel** (the D-lanes `D0/D2/D4` + `{B̃,Z̃}`) — the abstract `D0..N4` ARE 017's exported port
    scalars, carried as generic symbols. Cite `stage017`. ⚠ NO value consumed; NO dual-site (the algebra is port-agnostic).
  - **`build_port_moments`' concrete N-moments** — DEFERRED Gate-6 branch data (emitted, asserted-against-nothing).
  - **018's fingerprint context** — cite `stage018` (the exterior/interior two-port match; 019's `P(ω)` is the prefactor
    the 018 fingerprint feeds). Provenance narrative only.
- **Exports (→ 020/021/027):** the squared-denominator prefactor algebra `P₀=N₀/D₀`, `P₂=(D₀N₂−2D₂N₀)/D₀²`,
  `P₄=(D₀²N₄−2D₀(D₂N₂+D₄N₀)+3D₂²N₀)/D₀³` (SERIES-EXTRACTED) + the N/D self-check → **021** (the μ̂₀-free `[P₀^phys]=1` dim
  closure builds on `P₀=N₀/D₀`, source L399/L512) + **020** (the `54/5` partition's `P₀` context) + **027** (pathA_43
  closure slot). Per the cross-stage flow (`part2` L107): "018–021 export the Λ₂ fingerprint + χ_Q=1 + the 54/5 partition →
  022 + 027." Cite the exact export contract at note-authoring.

## §5 Teeth candidates (019-specific, per-tooth ablation MANDATORY)

1. **⭐ The prefactor-value teeth (`P₀=N₀/D₀`, `P₂=(D₀N₂−2D₂N₀)/D₀²`, `P₄=…`).** The derived `coeffs["P2"]` etc.
   (L219–223, `.wl` L79–81) are genuinely `coeff(ω,n)` off the actual series of `D₀N/D^cons²`; the L234–235 / `.wl` L85–90
   check is derive-vs-typed (independent). ⚠ **Per-tooth ablation MUST confirm these are NOT stamped** (the §4 firewall:
   forbidden to hardcode the formulas and "check" against them) — mutate the DERIVATION (e.g. corrupt the object to plain
   `N/D`, or perturb a `D_n`/`N_n` power) → the series-extracted coefficient changes → the match fires. This is the central
   EARNED tooth; the expansion of `D₀N/D^cons²` must be EMITTED (directive §4 firewall), not just the compare.
2. **⭐ The N/D squared-denominator self-check tooth (`3g_wrong_prefactor_object`, `FAIL_PREFACTOR_ALGEBRA`).** Replacing
   `D₀N/D^cons²` with plain `N/D` → `P₂=(D₀N₂−D₂N₀)/D₀²` (missing the factor of 2 on `D₂N₀`, difference `D₂N₀/D₀²≠0`) →
   `plain_equals_correct_P2=False` → `FAIL_PREFACTOR_ALGEBRA` fires; the correct object does NOT fire (self-ablation). ⭐
   Per-tooth: the `3g` self-ablation `prefactor_ablation` (L762) is a DYNAMIC re-run (`rerun_gate_logic:True`,
   `with_mutation`≠`without_mutation`), NOT a constant `self_ablation` (the v1 trip-up); ⚠ **re-scope the verdict to an
   019-local gate** (the prefactor match ∧ the N/D self-check), NOT the joint `base_verdict` (the 018 3a/3b re-scoping
   pattern). Confirm both directions fire: plain `N/D` → FAIL, correct object → NO_FAIL.
3. **The port-moments PROVENANCE integrity (light — NOT a tooth).** `build_port_moments`' N-moments are carried symbolic,
   DEFERRED Gate-6 branch data, asserted-against-nothing (§1c). ⚠ Confirm it is a **labeled non-check export**, NOT dressed
   as a checked relation (a "check" that reads back the deferred data would be a vacuous tooth). The provenance-cite
   integrity for 017's D-lanes / 018's fingerprint is a citation guard, NOT a dual-site.
4. **The `.wl` arity + native-route integrity (light).** Def/call arity scan + unevaluated-leakage transcript scan
   (stage007 lesson); confirm the `.wl` series-EXTRACTS its own `p0/p2/p4` (native `serW`/`Coefficient`), NOT a mirror of
   the `.py`.

⚠ **NOT 019 (do not rebuild as 019 teeth — 018/020/021 own these):** `3a_wrong_bc` / `3b_imposed_dissipation` (018, DONE);
`3c_wrong_scaling` (020), `3e_equivalence_break` (020), `3f_partition_mislabel` (020); `3d`/`3d′` dimensional-break (021).
⚠ **The `54/5=2·27/5` STRING label (`.py` L622–627) is 020's**; the μ̂₀-free `[P₀^phys]=1` gate + the `c_s/a` restoration
(`.py` L399/L512) are 021's — 019 owns NEITHER.

## §6 Register expectation — ⭐ THE KEY 019 QUESTION (likely ZERO new counted knobs; CONFIRM)

Per headline #2/#3 + the split: **019 is the EARNED prefactor-algebra slice (abstract squared-denominator structure over
017's exported port scalars); the CALIBRATED label + the `54/5` magnitude are 020's.** The honest pre-read (⚠ CONFIRM at
the register step + Codex-verify against the scripts):

- **⭐ 019 likely adds ZERO new counted knobs** (like 016/018 / 011/012/014 — an EARNED/structural slice). The squared-
  denominator prefactor algebra is pure abstract algebra over the ALREADY-EXPORTED port scalars `D0..N4` (017's ℓ=2 port
  kernel D-lanes + the deferred Gate-6 N-moments) — it introduces NO calibration and NO new units symbol (§headline-3:
  no `c_s`/`a`/`G`/`μ̂₀`). The port scalars `N_n/D_n` are 017's exported kernel (D-lanes) + `build_port_moments`' deferred
  branch data (`deferred_branch_data`) — NOT new knobs (they are NOT counted at 017 either: 017 counted `{T_Ω,β₂}`; the
  D-lanes/`{B̃,Z̃}` are "tracked/downstream-pinned").
- **⭐ Likely a new STRUCTURAL edge (call it R38 — confirm the next free number at registration; R37 was 018's):** the
  squared-denominator prefactor-algebra provenance — `P(ω)=D₀N/D^cons²` yields `P₀=N₀/D₀`, `P₂=(D₀N₂−2D₂N₀)/D₀²`, `P₄=…`
  (SERIES-EXTRACTED; the `−2D₂N₀` factor-of-2 is the squared-denominator signature, provably absent from plain `N/D`) —
  **discharges NOTHING** (earned abstract-algebra structure over the deferred port kernel, not a reduction of a debt; like
  R37 for 018's fingerprint, R34 for 016's covariance theorem). Confirm the edge captures: the algebra is port-agnostic
  (holds for any `D0..N4`), the concrete port instance (`build_port_moments`) is DEFERRED Gate-6 branch data, and the N/D
  self-check is the anti-mislabel tooth.
- **Cited provenance (NOT re-counted):** 017's ℓ=2 port kernel D-lanes + `{B̃,Z̃}` (tracked/downstream-pinned at 017;
  NOT counted), `build_port_moments`' N-moments (deferred Gate-6 branch data), 018's fingerprint context. ⚠ **Do NOT let
  019 count the port-kernel `N_n/D_n` scalars** (deferred branch data) or re-count 017's `{T_Ω,β₂}` / 013's `{μ_η,T_w,β}`.
- **Control/tracked-not-counted:** the sample subs `{D0:19, D2:23, D4:29, N0:11, N2:13, N4:17}` (numeric evaluation
  controls, like 018's `{a:3, c_s:2}` and 014's controls).

⚠ **Do NOT let 019 silently count the port scalars.** Resolve them as 017's exported kernel (D-lanes, tracked-not-counted)
+ deferred Gate-6 branch data, and Codex-verify (the register verify is the gate that catches an over-count that would
falsely inflate — or a mislabel that would falsely shrink — the irreducible codimension count). **Part-II counted CALIB
set stays 6** (`{μ_η,T_w,β}`(013)+`{Vp0/ℓ_c}`(015)+`{T_Ω,β₂}`(017)); 019 adds NONE.

## Verdict tokens + honest scope

019 carries the **squared-denominator prefactor-algebra component (2/4) of `QUAD_CALIBRATED`**: `P(ω)=D₀·N/D^cons²` yields
`P₀=N₀/D₀`, `P₂=(D₀N₂−2D₂N₀)/D₀²`, `P₄=(D₀²N₄−2D₀(D₂N₂+D₄N₀)+3D₂²N₀)/D₀³` (SERIES-EXTRACTED from the actual series of the
squared-denominator object — NOT typed), and the N/D self-check (plain `N/D` gives `−D₂N₀` vs the correct `−2D₂N₀` — the
factor of 2 only the squared denominator produces → `FAIL_PREFACTOR_ALGEBRA`). EARNED = the prefactor structure is DERIVED
(series-extracted) and the squared-denominator object is provably the right one (the plain-`N/D` object provably wrong).
The fingerprint = 018 (DONE); the `54/5=2·27/5` provenance partition + the CALIBRATED label + `G=GENUINE_BLOCKED` = 020;
the μ̂₀-free `[P₀^phys]=1` dim closure = 021. **Consumption is PROVENANCE only** — the abstract `D0..N4` ARE 017's exported
ℓ=2 port kernel scalars (D-lanes) carried as generic symbols, the algebra is PORT-AGNOSTIC, and `build_port_moments`'
concrete N-moments are DEFERRED Gate-6 branch data asserted-against-nothing (like 018 — no checkable cross-stage relation,
no dual-site). 019 is c_s-FREE / a-FREE / G-FREE / μ̂₀-FREE (no units symbol, no dimensional leg — those are 021's).
Caveats: the port scalars `N_n/D_n` remain symbolic/Gate-6 (the numerical branch data deferred, report :49); the `54/5`
magnitude + `G` are CALIBRATED (020, `G=GENUINE_BLOCKED`). ⭐ The pathA_33 trip-ups: the `27` stays 018's (019 does NOT
touch the `27`/`54/5`/`G`); the v1 rig (back-solved `μ̂₀` + constant `self_ablation`) must NOT return — 019's prefactor
slice is `μ̂₀`-free and its `3g` self-ablation is a DYNAMIC re-run.

## Process (unchanged, calibrated — the per-stage pipeline)

Author the II-G4b reshape directive (§1 the clean 019 slice / 4-way cut + §1b the native-`.wl` KEEP + §1c the provenance-
only port-kernel framing + §2 faithful cover + §3 bridge-strip incl. sever-YAML + transcript-level agreement + §5 the
prefactor-value / N/D-self-check teeth with per-tooth ablation + §6 the ZERO-new-knobs + R38-edge register question) →
**Codex xhigh design-review** → fold to `DIRECTIVE_CLEAN` → **⭐ final Grok-4.5 headless compute-verify pass** (Grok
compute-verifies the squared-denominator series → `P₀/P₂/P₄`, the `−2D₂N₀` vs plain `−D₂N₀` factor-of-2, and the port-
moments emitted-but-never-checked reading; it caught the 016 volume-vs-line convention mismatch + a kernel-preserving
residual-tooth on 013, so watch the series-extraction genuineness + the provenance-only port-kernel resolution here) →
fold → **Codex confirm-pass on the folds** → **pre-exec USER GATE** → Codex builds the two scripts (`--sandbox
danger-full-access`, background, xhigh) → dual-engine both exit 0 (repo root + foreign CWD) → arbiter re-run → full
tri-review (fidelity + adversarial-with-**per-tooth ablation**; ⭐ hunt the prefactor-value series-extracted-vs-typed
genuineness + the N/D factor-of-2 self-check + a mirror-`.wl` + any vacuous able-to-fail, esp. a port-moments "check" that
reads back deferred data) → remediate → fresh-agent re-verify → bump counts 18→19 → parameter register (⭐ confirm ZERO
new 019 knobs + port scalars = 017's exported kernel / deferred Gate-6 branch data + R38 edge) + Codex-verify → note/card/
`\input{stages/stage_019}` + registration → rebuild PDF → commit + docs/memory sync. Orchestrator authors notes/cards/
LaTeX/registration; Codex codes. Target stem: `ledger_stage019_prefactor_algebra` (confirm slug at authoring).
