# II-P4 (ledger_stage027) source map — pathA_43 density-port CHECKS + CLOSURE (`DENSITY_PORT_HOSTED` 4/4, the COMPLETING leg that OWNS the joint)

> Running-start prep captured 2026-07-10 (post stage026 = pathA_43 II-P3 continuity-lineage token-check DONE, commit
> `22ca42c7`; Cluster C OPEN at 3 of 4), before authoring the II-P4 reshape directive, so the directive can be written
> without re-discovery. **All line refs below VERIFIED against the current sources 2026-07-10** by a fresh-agent
> full read of BOTH engines (`pathA_43_density_quadrupole_port_sympy.py` = **860 lines**, `.wl` = **545 lines**) covering
> the 027 slice (the dim engine, `dtn_sign`, `closure_overlay`, `evaluate`/`assert_gate`/`controls`) PLUS a read of the
> three BUILT sibling scripts (024/025/026) to pin the exact host contract, the LOCAL/JOINT verdict strings, and the
> forward-ref exports 027 consumes — corroborated against the stage024/025/026 source maps' 027-slice pointers. No line
> drift (pathA_43 sources untouched since the 2026-07-06 commit).
> Companion: `part2_gravity_atomic_split.md` (Cluster C rows **024–027** L49–52 = the 4-way pathA_43 split; the **▶ NEXT
> stage027 entry L636–647**; the **per-gate trip-ups L93** [⚠ "pathA_43 (025/026): the computed-taint + host-set guard +
> lineage token-check are THE earned mechanisms — never collapse to name-checks or flags"; the pathA_33 (020/021) trip-up
> "the 27 stays COMPUTED; 54/5 asserted only as labeled `external_bridge_input`; the v1 rig back-solved free-carrier μ̂₀ +
> constant self_ablation must not return — μ̂₀-free dim gate + real two-verdict self-ablations"]; the reshape-cost map
> **L78–79** [pathA_43 = lightest reshape, contract-clean]; the cross-stage flows L104–109) and the **stage024/025/026
> pairs + source maps + directives** (the pathA_43 SIBLING legs 027 CONSUMES from: 024's `N0_den` + the `physical_relations`
> provenance + the density host-set; 025's vector-freedom verdict + `source_map_complete` certificate; 026's `moment_valid`
> + validated `I25` + lineage certificate) + **stage021** (the CLOSEST μ̂₀-free dim-check exemplar — 021's
> `[P0_phys]=(c_s/a)²·N₀/D₀` dimensional gate + `scale_power`/`BASE_DIMS` machinery + the corrupt-`[N₀]`/`[G]` truth-table +
> the μ̂₀-back-solve-is-a-tautology demotion is EXACTLY 027's dim/scaling slice, KEPT-NATIVE `.wl`) + **stage018** (the
> `dtn_sign`/`+i z⁵/27` exterior-Hankel fingerprint + `χ_Q=+1` outgoing / `−1` incoming — 027 CONSUMES the sign check;
> 018 KEPT its native Hankel `.wl`) + **stage020** (the `54/5=2·27/5` partition + `G=GENUINE_BLOCKED` + the CALIBRATED
> verdict — 027's closure `Γ̄₅=2G/(5c⁵)` rides the same calibrated magnitude). Build-order id **027**, Part II, Cluster C.
> Source top-line (verbatim, report `reports/pathA_43_density_quadrupole_port.md:1`): **`DENSITY_PORT_HOSTED`** — the JOINT
> verdict of the whole pathA_43 gate; **027 OWNS it** (the COMPLETING leg — unlike the PARTIAL legs 024/025/026 which
> printed the JOINT as a 1/4·2/4·3/4 partial). Proposed target stem: `ledger_stage027_port_checks_closure` (confirm slug at
> directive authoring).

## ⭐ The FIVE headline points (READ FIRST)

1. **⚠ 027 is the COMPLETING leg of the pathA_43 4-way split (024/025/026/027) — the ONE leg that OWNS the joint verdict.**
   024 DERIVED the density two-port `N0_den`; 025 PROVED it vector-free; 026 EARNED the continuity lineage of its ℓ=2
   moment `I25`. **027 runs the 6 able-to-fail PORT CHECKS + the K̄ CLOSURE overlay, ASSEMBLES the joint `origin_ok` /
   `DENSITY_PORT_HOSTED` from the 025 (vector-freedom) + 026 (continuity-lineage) certificates + its OWN
   dim/scaling/sign/nonzero/closure, and RECORDS the EM `A_w`/`U,W` scaffold RETIREMENT with the completed joint.** This is
   the DISTINCTIVE role — 024/025/026 each printed a JOINT PARTIAL (they do NOT own the joint); **027 emits the joint
   `DENSITY_PORT_HOSTED` as its OWN verdict** (the 023-pattern: the COMPLETING leg owns the joint; the 018/024 EARNED-first
   legs do not). ⚠ **COMPLETE ≠ PASS:** the joint token STAYS **CALIBRATED not PASS** (`G=GENUINE_BLOCKED`, 020's
   provenance) — 027 verifies the port is density-HOSTED (STRUCTURE/form earned, vector-free, continuity-lineaged,
   dimensionally + radiatively + scaling well-formed, closure-consistent), the MAGNITUDE stays CALIBRATED/SIM_DEFERRED. ⚠
   **The retirement is CONDITIONAL on 027 passing:** R43/R44 explicitly DEFERRED the EM-scaffold retirement to this joint —
   **if 027 FAILS dim/scaling/sign/closure the density port has NOT displaced the scaffold** (the joint would land a `FAIL_*`
   / `PORT_INCONCLUSIVE_SIM_DEFERRED` and the diagnostic sliver reopens). The gate must be able to emit that.
2. **⭐⭐ THE EARNED CONTENT = the 6 COMPUTED port checks + the K̄ closure overlay, each able-to-fail (never a typed/stamped
   pass).** The checks (`.py` `evaluate` L537–698; `assert_gate` L773–803):
   - **(i) the dimensional gate `dim_ok`** (L559–569) — `[N0_den] == (−1,1,0) = L⁻¹M` via the dim engine `dim_of`
     (L82–108) over `BASE_DIMS` (L200–233). ⚠ **Dimension tuple convention is `(L,M,T)`** (position 0=L, 1=M, 2=T;
     `dim_to_text` L116–123 / `.wl` L56–62) — NOT the register header's `[L,T,M]`; carry the convention explicitly (the
     stage016 volume-vs-line convention lesson: a cross-convention dim identity may NOT transfer — cite as PROVENANCE +
     rely on self-contained integrity).
   - **(ii) the a⁻⁵ scaling gate `scaling_ok`** (L571–576) — `P0_physical = (c_s/a)²·N0_den/D0` must carry a-power `−5`
     (`p0_power == −5`; `scale_power` L126–158 over `BASE_A_POWERS` L236–268). This is 021's `μ̂₀`-free machinery
     (`[P0_phys]=1` / a-power `−5`) — CITE 021, do NOT re-derive; the pathA_33 (020/021) trip-up bans the v1 back-solved-μ̂₀
     rig (μ̂₀-free gate + real two-verdict self-ablations).
   - **(iii) the outgoing DtN sign gate `sign_ok`** (`dtn_sign` L521–534, verdict L577–578) — the `+i z⁵/27` fingerprint:
     `coeff = series(ŷ₂,z⁷).coeff(z,5)/i == 1/27` for the OUTGOING branch (`h₂=j₂+i·y₂`); the INCOMING branch (`j₂−i·y₂`)
     flips the imaginary sign → `−1/27` → `sign_ok=False`. This is 018's `χ_Q=+1` outgoing / `−1` incoming — CONSUME 018's
     fingerprint (do NOT re-derive the whole series; 027's gate is the SIGN of the radiating coefficient).
   - **(iv) the nonzero-port gate `nonzero_ok`** (L615) — the port must be nonzero (`FAIL_PORT_VANISHES` when the coupling
     vanishes; the `zero_coupling` control). ⚠ Guarded by `vanished_continuity_coupling` (L593–599) so a legitimately
     vanishing coupling fails at `nonzero_ok` (`FAIL_PORT_VANISHES`), NOT at `origin_ok` (which would wrongly read
     `FAIL_NOT_DENSITY_DERIVED`) — 027 must reproduce that arm or the `zero_coupling` verdict flips (agent §E-3).
   - **(v) the `deferred_scalar` SIM-inconclusive branch** — the HONEST landing: if the branch scalar is uncertified
     (`deferred_uncertified=True` → `xi=Xi_deferred`, `powers[Xi_deferred]=None` → `scale_power → None` →
     `scaling_undecidable=True` L575), the verdict falls through EVERY FAIL/PASS branch to the final `else` →
     **`PORT_INCONCLUSIVE_SIM_DEFERRED`** (L633–634). The converse `proven_deferred=True` (`powers[Xi_deferred]=0` L557 →
     decidable) → `DENSITY_PORT_HOSTED`. This branch is neither a FAIL nor a spurious PASS — it is the "the sim-deferred
     literal is not yet certified" outcome; both controls must fire (the branch is honest, not baked toward hosted).
   - **(vi) the K̄ CLOSURE overlay `closure["ok"]`** (`closure_overlay` L701–722) — `K̄₄ − 4K̄₂²/K̄₀ == 0` ∧
     `Γ̄₅ − 2G/(5c⁵) == 0`, with `K̄₀=54Gc_s⁵/(5a⁵c⁵)`, `K̄₂=6Gc_s³/(5a³c⁵)`, `K̄₄=8Gc_s/(15ac⁵)`, `Γ̄₅=K̄₀·a⁵/(27c_s⁵)`
     (the `/27` ties to the DtN `z⁵/27`). ⚠ **This is the A3 2.5PN Burke–Thorne boundary — SHARED with stage028's
     match-back** (028 runs the full INV1–INV5 + the 11-mutation coherent-rescale matrix). Mark 027's closure slot as the
     PORT-closure CONSISTENCY (the two residuals) and 028 as the full match-back — **SHARED, NOT double-counted** (the
     directive must state the 027/028 cut so 027 does not rebuild 028's INV5 literal anchors and 028 does not rebuild the
     port derivation).
3. **⭐⭐ THE JOINT ASSEMBLY — `origin_ok` is a COMPOSITE 027 assembles from 025 + 026, NOT a pure-027 object (the #1
   design decision).** `origin_ok` (`.py` L600–609) = `("continuity_interface" ∈ tags)` [**026** tag] ∧
   `("vector_port" ∉ tags)` [**025**] ∧ `(not vector_host_symbols)` [**025**] ∧ `source_map_complete` [**025**] ∧
   `continuity_dependency_ok` [**026**] `OR vanished_continuity_coupling`. The verdict ALSO needs `vector_independence_ok`
   [**025**, L610–614] ∧ `nonzero_ok` ∧ `dim_ok` ∧ `scaling_ok` ∧ `sign_ok` ∧ `closure["ok"]` [027's own]. ⭐ **027 owns
   the verdict-ASSEMBLY machinery (`evaluate` L537–698, the verdict chain L618–634, `assert_gate` L773–803, `controls`
   L725–770) + the 5 NEW port checks + the closure — but it CONSUMES 025's vector-free certificate + 026's continuity
   certificate to form the joint.** The directive must state the cut: 027 does NOT rebuild the taint graph
   (`source_tag_map`/`taint_for_expr`/`vector_ablated_expr` = 025's) or the lineage token-check
   (`continuity_lineage_valid`/`continuity_moment_symbol` = 026's) — it CITES their earned verdicts (the exported
   `vector_free`/`source_map_complete` from 025; the `moment_valid`/validated `I25`/`continuity_dependency_ok` from 026) and
   assembles the joint. ⚠ How much of `origin_ok` 027 RE-COMPUTES vs CITES is the load-bearing directive decision (see §1d)
   — the cleanest cut: 027 CONSUMES the two certificates as booleans (with a checkable consumption-integrity tie to the
   real `N0_den`), and runs its OWN dim/scaling/sign/nonzero/closure teeth; a FAIL of a consumed certificate flows to the
   joint (the gate is able-to-fail if a sibling's certificate is corrupted).
4. **⭐ THE a-power `−7/2` SEAM (026↔027) — 027 tests the SCALING consequence, NOT the moment-earning (026's).** Config
   default `coupling_a_power = Rational(−7,2)` (L170). The `scaling` control (`coupling_a_power = −3`) trips BOTH concerns:
   it flips the earned moment to `I_wrong2` (026's `continuity_moment_symbol` L325–332) AND breaks `p0_power`. ⚠ **The
   subtle design 027 MUST preserve (agent §E-1): `I_wrong2` dim `(2,0,0)` + `a⁻³` DELIBERATELY COMPENSATE `I25` dim
   `(5/2,0,0)` + `a^(−7/2)` in the L-slot** (`2−3 = 5/2−7/2 = −1`), so `dim_ok` STAYS True and ONLY `scaling_wrong` fires →
   **`FAIL_PORT_MALFORMED(scaling)`** (single-tag, L618–634's comma-joined bad list). If 027 rebuilds the scaling control
   without honoring that compensation it will wrongly emit `FAIL_PORT_MALFORMED(dimensional,scaling)`. 027 CONSUMES 026's
   earned moment / `moment_valid` and tests only the scaling/dimension consequence — it does NOT re-litigate the
   `I25`-vs-`I_wrong2` earning (026 owns it; 027 uses the `−7/2` baseline as the port's moment).
5. **⭐ CONTRACT-CLEAN: there is NO bridge to sever (the lightest reshape class in Part II, shared with 024/025/026).** BOTH
   engines are already standalone print-only, zero file I/O (`.py` imports only `sympy`/`dataclasses`/`typing`; `.wl` grep
   for `Export|Put|Write|OpenWrite|>>|Save|Import` = ZERO hits). So the reshape "sever the bridge" step is a **no-op**. The
   027 work is: (a) the DECOMPOSITION (extract the dim engine + `scale_power` + `dtn_sign` + `closure_overlay` +
   `evaluate`/`assert_gate`/`controls` into a self-contained script that CONSUMES 024's `N0_den` + 025's vector-free
   certificate + 026's `moment_valid`/`I25`, runs the 6 checks + closure, ASSEMBLES the joint, and DROPS the 024
   derivation / 025 taint graph / 026 lineage token-check); (b) the `.wl` decision (§1b — the dim engine + `dtn_sign` are
   the SAME algebra 018/021 KEPT-NATIVE, so keep-native is DEFENSIBLE for those; `closure_overlay` + the `evaluate`
   assembly must pass the transliteration screen — route the per-function determination through the Codex→Grok→Codex
   bookend); (c) the JOINT-OWNING verdict framing (027 OWNS `DENSITY_PORT_HOSTED`, the 023-pattern — NOT a printed PARTIAL);
   (d) the per-tooth-ablation control battery over EACH of the 6 checks + the closure; (e) the LOCAL-ledger idioms
   (`banner`/`expect_zero`/`expect_bool`/`RigAssertion`/`routed_assert`/`exercise_rig`/tally/nonzero exit) replacing the
   pathA_43 `assert_gate` monolith. ⚠ **`.py`-only constructs to DROP:** the `non_dodge` table (L636–655, `.py`-only,
   documentary) and the `powers[Xi_Q]=0` no-op (L550–553); the `.wl` has neither.

## §1 The 027 slice (line ranges) — the CLEAN CUTS (all VERIFIED 2026-07-10)

pathA_43 is ONE script doing all four buckets; 027 owns the CHECKS + CLOSURE + the joint verdict assembly. **027 owns the
dim engine + `scale_power` + `dtn_sign` + `closure_overlay` + `evaluate`/`assert_gate`/`controls` (the verdict machinery)
+ the 6 port-check controls; it does NOT own the derivation (024), the taint graph (025), or the lineage token-check
(026) — it CONSUMES their exported certificates.**

### §1a The SymPy 027 slice (`.py`)
- **⭐⭐ The dimension engine (027's core, KEEP-NATIVE class per 021):**
  - `Dim` alias L21; `ZERO_DIM=(0,0,0)` L22; `dim_add` L70–71, `dim_sub` L74–75, `dim_scale` L78–79.
  - **`dim_of` L82–108** — recursive over `Mul`/`Pow`/`Add`; raises `DimError` on a symbol-not-in-`dims`, a non-numeric
    exponent, or an `Add` mismatch (the corrupt-`[N₀]` tooth rides this, à la 021).
  - `exp_text` L111–113, `dim_to_text` L116–123 (⚠ prints `L,T,M` from the `(L,M,T)` tuple — the convention seam).
  - **`scale_power` L126–158** — the a-power engine (`a→1` hardcoded L130–131; returns `None`=undecidable on a
    non-numeric exponent / `Add` mismatch → the `deferred_scalar` inconclusive branch).
  - **`BASE_DIMS` L200–233** — the port-symbol base dims: `a:(1,0,0)`, `c_s:(1,0,−1)`, `c:(1,0,−1)`, `G:(3,−1,−2)`,
    `D0:(−1,1,−2)`, `rho:(−3,1,0)`, `I25:(5/2,0,0)`, `I_wrong2:(2,0,0)`, `Xi_Q/Xi_deferred/eta_q/eta_phi:ZERO_DIM`
    (L209–212), `varpi_q2/varpi_Phi2/lambda_c:(0,0,−2)` (L213–215). ⭐ Expected `[N0_den] = (−1,1,0) = L⁻¹M` (asserted
    L565, L781).
  - **`BASE_A_POWERS` L236–268** — all density scalars `0`; `varpi_q2/varpi_Phi2/lambda_c:−2` (L248–250); the vector
    symbols documentary (off the density path).
- **⭐ `dtn_sign(z, outgoing) L521–534`** — `j2` L522, `y2` L523; `h=j2+I·y2` (outgoing) else `j2−I·y2` L524; `lam=z·h'/h`
  L525; `yhat=−3/lam` L526; `series(yhat,z,z⁷)` L527; `coeff=series.coeff(z,5)/I` L528; asserts
  `matches_outgoing = compact(coeff − 1/27) == 0` L533. (018's fingerprint; the sign gate.)
- **⭐ `closure_overlay(n0) L701–722`** — `p0_physical=(c_s/a)²·n0/D0` L704; `kbar0=54Gc_s⁵/(5a⁵c⁵)` L705;
  `kbar2=6Gc_s³/(5a³c⁵)` L706; `kbar4=8Gc_s/(15ac⁵)` L707; `kbar4_residual=kbar4−4·kbar2²/kbar0` L708;
  `gamma5=kbar0·a⁵/(27c_s⁵)` L709; `gamma5_residual=gamma5−2G/(5c⁵)` L710; returns `ok = (kbar4_residual==0 and
  gamma5_residual==0)` L721. (Both residuals are IDENTICALLY 0 by the 020 magnitude — VERIFIED: `4K̄₂²/K̄₀ = (8/15)Gc_s/(ac⁵)
  = K̄₄`; `Γ̄₅ = 54G/(5·27·c⁵) = 2G/(5c⁵)`. L713 records `"1 with chi_Q=1, mhat0=1, N_Q=1"`.)
- **⭐⭐ `evaluate(config) L537–698` — THE VERDICT MACHINERY (027's core assembly):**
  - dims + `corrupt_dimension` L545–547; powers + `deferred`/`proven` overrides L549–557 (⚠ the `powers[Xi_Q]=0` no-op
    L550–553 — DROP).
  - `dim_ok` L559–569; `n0_power`/`p0_power`/`scaling_wrong`/`scaling_ok`/`scaling_undecidable` L571–576; `sign_ok`
    L577–578.
  - `vector_host_symbols`/`source_map_complete` L580–584 (**025's** — consumed); `continuity_dependency_ok` L585–592
    (**026's** — consumed); `vanished_continuity_coupling` L593–599 (027's `coupling_zero` blended with 025+026 predicates);
    **`origin_ok` L600–609 (the COMPOSITE 027 assembles)**; `vector_independence_ok` L610–614 (**025's** — consumed).
  - `nonzero_ok` L615; `closure = closure_overlay(expr)` L616.
  - **the verdict chain L618–634:** `not origin_ok or not vector_independence_ok → FAIL_NOT_DENSITY_DERIVED`;
    `not nonzero_ok → FAIL_PORT_VANISHES`; `not dim_ok or scaling_wrong or not sign_ok → FAIL_PORT_MALFORMED(<comma-joined
    bad list: dimensional,scaling,sign>)`; all-true incl `closure["ok"] → DENSITY_PORT_HOSTED`; else →
    `PORT_INCONCLUSIVE_SIM_DEFERRED`.
  - `non_dodge` table L636–655 (`.py`-only — DROP); return dict L657–698.
- **⭐ `assert_gate` L773–803** — baseline `DENSITY_PORT_HOSTED` (L775); all checks true; taint ==
  `["continuity_interface","pathA_29_bulk","pathA_32_wall"]` (L778); `N0_den == "L^-1 M"` (L781);
  `P0_physical_a_power == "-5"` (L782); `Kbar4_residual == 0` (L783); then the `expected` control-verdict dict L785–800.
  ⚠ `Γ̄₅` residual is gated ONLY through `closure["ok"]` (no standalone assert) — 027 should add a standalone Γ̄₅ residual
  assert so BOTH closure residuals are independently ablatable (agent §E-4).
- **Symbol declarations 027 USES (of L24–45):** the density host-set `a, c_s, rho, I25, Xi_Q, eta_q, eta_phi, varpi_q,
  varpi_phi, lambda_c, q2, Phi2` (in the consumed `N0_den`), PLUS `c` (GR-units bridge, in the closure), `G` (`GENUINE_BLOCKED`,
  in the closure), `D0` (the reduced static conservative denominator, in `P0_physical`), `z` (the DtN variable — LIVE in
  `dtn_sign`, the ONE engine where `z` is not dead), `Xi_deferred` (the deferred-scalar branch), `I_wrong2` (the scaling
  control's wrong-symbol, dim `(2,0,0)`). ⚠ **NOT 027 (drop / leave to the siblings):** the vector scaffold `A_w,…,g_W`
  (025's control fixtures), the relabel-fixtures `omega_wall,…,g_qold` (025's), `sigma_hidden`/`free_carrier` (025's),
  `omega` (dead both engines).
- **⭐ The Config fields 027's control battery drives (L161–177):** `corrupt_dimension` (L166, the `dimensional` control),
  `incoming_sign` (L167, the `sign` control), `coupling_a_power` (L170, the `scaling` control at `−3`; default `−7/2`),
  `coupling_zero` (L165, the `zero_coupling`/nonzero control), `deferred_uncertified`/`proven_deferred` (the
  `deferred_scalar` + converse controls). These are the rig SOURCES; 027 re-expresses them as its able-to-fail teeth (§5).
- **The `controls()` fixtures + expected verdicts 027 OWNS (L725–803):** `zero_coupling` → `FAIL_PORT_VANISHES`;
  `dimensional` → `FAIL_PORT_MALFORMED(dimensional)`; `sign` → `FAIL_PORT_MALFORMED(sign)`; `scaling` →
  `FAIL_PORT_MALFORMED(scaling)`; `deferred_scalar` → `PORT_INCONCLUSIVE_SIM_DEFERRED`; `deferred_scalar_proven_converse`
  → `DENSITY_PORT_HOSTED`. (The other 8 controls — `vector_hosted`/`relabel_rig`/`ablation_isolation`/`attack5…`/
  `provenance_less_rider`/`free_carrier_dimension_corruption` = 025's; `fake_continuity`/`attack2_continuity_corruption` =
  026's — 027 does NOT re-run their machinery; a corrupted sibling CERTIFICATE flowing to `FAIL_NOT_DENSITY_DERIVED` is the
  joint-assembly tooth, §5.)
- **Shared helpers 027 uses (NOT cut boundaries):** `compact` L52–57, `hstr` L60–63, `bstr` L66–67. ⚠ **The taint machinery
  `VECTOR_SYMBOLS`/`BASE_SOURCE_TAGS`/`source_tag_map`/`taint_for_expr`/`source_graph_for_expr`/`vector_ablated_expr`
  (L179–190, 271–302, 335–389) are 025's; the lineage machinery `CONTINUITY_*`/`contains_all`/`continuity_lineage_valid`/
  `continuity_moment_symbol`/`lineage_for` (L192–197, 305–332, 392–417) are 026's; `schur_density_expression` L420–465 is
  024's** — 027 CONSUMES `N0_den` + the two certificates, runs NO taint/lineage/derivation.
- **⭐ CLEAN CUT (SymPy) — 027 owns L70–158 + L200–268 + L521–534 + L537–698 + L701–722 + L725–803 (the dim engine +
  `dtn_sign` + `closure_overlay` + `evaluate` + `assert_gate` + `controls`), plus the minimal driver + emit.** It touches
  NONE of `schur_density_expression` (024), the taint machinery (025), the lineage machinery (026). ⚠ **The mixed driver
  `derive` L482–518 is the SEAM** (it wires 026 lineage + 025 tag-map + 024 schur-expr + 025 source-graph → the trace/expr
  `evaluate` reads); 027 keeps the `evaluate`-facing wiring (it needs `N0_den` + the tags + the certificates) but the
  underlying taint/lineage/derivation FUNCTIONS are cited-not-rebuilt (consume their exported booleans; §1c).

### §1b The Mathematica 027 slice (`.wl`) — ⚠ MIXED determination (dim engine + `dtnSign` keep-native class per 018/021; the assembly + closure need the transliteration screen)
- **The source `.wl` 027 slice:** the dim engine `zeroDim` L51, `dimAdd/Sub/Scale` L52–54, `dimToText` L56–62, `dimOf`
  L64–86, `scalePower` L88–115, `baseDims` L117–129, `basePowers` L131–139 (same `(L,M,T)` convention + same values);
  `dtnSign` L341–351 (`j2`/`y2`/`h=j2±I y2`/`lam=z D[h,z]/h`/`yhat=−3/lam`/`Series[…,{z,0,6}]`/`Coefficient[…,z,5]/I`/
  `TrueQ[FullSimplify[coeff==1/27]]`); `closureOverlay` L353–364 (same K̄₀/K̄₂/K̄₄/Γ̄₅ + both residuals); `evaluate`
  L366–442 (the verdict `Which` L414–426, identical branch order/strings to `.py`); `assertGate` L465–494; `controls`
  L444–463. ⚠ **No `non_dodge` table, no `powers[xiQ]=0` no-op** (both `.py`-only — keep dropped).
- **⭐ The `.wl` determination (a directive-bookend decision, NOT pre-decided here):** the dim engine + `scalePower` +
  `dtnSign` are the SAME kind of exterior-Hankel / dim-tuple algebra that **stage018 KEPT-NATIVE** (native `j2`/`y2` +
  `Series`/`Coefficient` + own `dimOf`) and **stage021 KEPT-NATIVE** (native `Which`-based `dimOf` + `rawDims`
  Association + `Series`-free dim algebra) — so **keep-native is DEFENSIBLE for the dim engine + `dtnSign` + the dim gate**
  (they are not `.py` mirrors — Wolfram `Series`/`Coefficient`/`SphericalBesselJ` vs SymPy). **But `closureOverlay` +
  the `evaluate` verdict assembly may read as a transliteration** (identical branch order/strings) → the
  `MATHEMATICA_MIRROR_POLICY` screen must be applied PER-FUNCTION at the directive bookend. ⭐ Candidate independent
  routes (suggestions for the directive, NOT a pre-design — Codex owns the code): for `dtnSign`, a native
  `SphericalHankelH1`/`SphericalHankelH2` + `SeriesCoefficient` route (like stage022's re-author) instead of the
  hand-built `j2+I y2`; for the closure, a native `Cancel`/`Together`/`FullSimplify` residual on independently-typed
  K̄ moments; for the verdict, an `Association`/`Which` assembly that reads the consumed certificates. ⚠ **Directive routes
  the per-function keep-native-vs-re-author determination through the Codex→Grok→Codex bookend (the transliteration screen
  EXPLICITLY applied), function by function.**
- **`.wl` shared/emit helpers:** `fail`/`assertTrue` L14–15, `heading`/`subheading`, `clean`, the free-symbol extractor.
  ⚠ **Arity discipline** (the stage007 lesson): a def/call arity scan + an unevaluated-leakage transcript scan (the 027
  `.wl` has several `Module`s in `evaluate`/`closureOverlay`/`dtnSign`).

### §1c The consumption resolution (027 CONSUMES 024's `N0_den` + 025's + 026's certificates)
⭐ **027's SUBJECT is the assembled joint over 024's `N0_den`.** Since there is zero file I/O, 027 CITES 024's factored
export + the two sibling certificates (each a boolean/set, cited-not-rebuilt), with a checkable consumption-integrity
tie so the assembly runs over the REAL objects:
- **⭐ stage024 `N0_den`** (the density two-port numerator) — the object the dim/scaling/sign/nonzero/closure checks run
  over: `N0_den = I25²·Ξ_Q²·c_s⁴·rho_eff·(η_φ·ϖ_q2 + η_q·λ_c)²/(a⁷·(λ_c² − ϖ_Φ2·ϖ_q2)²)` (stage024 report `:9`; built
  stage024 `.py` `make_N0(response)`). 027 CITES the factored form + asserts its computed `free_symbols == the 10-symbol
  HOST_CONTRACT `{a, c_s, rho_eff, I25, Ξ_Q, η_q, η_φ, ϖ_q2, ϖ_Φ2, λ_c}`` (the consumption-integrity oracle, the same one
  025/026 use). ⚠ `rho = rho_eff` (reduced-3D MASS density `[M L⁻³]`, NOT stage005's `ρ0`; carry the rename). The dim/scaling
  checks run over the REAL `N0_den` so the corrupt-`[N₀]` and scaling teeth are genuine.
- **⭐ stage025 the vector-freedom certificate** — the exported `vector_free` / `source_map_complete` / `vector_port ∉
  taint` (025's `DENSITY_PORT_VECTOR_FREE` verdict). 027 CONSUMES these as booleans feeding `origin_ok` /
  `vector_independence_ok`. ⚠ 027 does NOT rebuild the taint graph (025's) — it CITES the certificate; a corrupted
  certificate → `FAIL_NOT_DENSITY_DERIVED` (the joint-assembly tooth). ⭐ 025's `DENSITY_PORT_VECTOR_FREE` was CONDITIONAL
  on `moment_valid=True` — 026 discharged that; 027 assembles the now-unconditional certificate.
- **⭐ stage026 the continuity-lineage certificate** — the exported `moment_valid=True` + validated `I25` +
  `continuity_dependency_ok` (026's `CONTINUITY_LINEAGE_EARNED`; the `lineage_certificate=PASS` L659 export). 027 CONSUMES
  these feeding `origin_ok` (`continuity_interface ∈ tags` + `continuity_dependency_ok`). ⚠ 027 does NOT rebuild the token
  check (026's) — it CITES; a corrupted certificate → `FAIL_NOT_DENSITY_DERIVED`.
- **⭐ stage018 the DtN fingerprint** (`χ_Q=+1` outgoing / `−1` incoming, `+i z⁵/27`) — 027's `dtn_sign` gate is the SIGN
  of the radiating coefficient (`coeff/i == 1/27` outgoing). CONSUME 018 (the fingerprint is 018's earned content; 027's
  gate is the outgoing-branch selection + the sign). PROVENANCE + the sign check.
- **⭐ stage021 the μ̂₀-free dim machinery** (`[P0_phys]=(c_s/a)²·N₀/D₀` = 1, a-power `−5`; `BASE_DIMS`/`scale_power`;
  the corrupt-`[N₀]`/`[G]` truth-table; the μ̂₀-back-solve-is-a-tautology demotion) — 027's dim/scaling gate IS 021's
  machinery applied to the density port. CITE 021's dim-check pattern (do NOT re-derive the `[P0_phys]=1` closure logic);
  the pathA_33 (020/021) trip-up L93 forbids the v1 back-solved-μ̂₀ rig + constant self_ablation.
- **⭐ stage020 the `54/5=2·27/5` partition + `G=GENUINE_BLOCKED`** — the closure `Γ̄₅=2G/(5c⁵)` rides the SAME calibrated
  magnitude (`54/5` = `external_bridge_input`, the `27` = `derived_in_gate` from 018). CITE 020 (the closure's `Γ̄₅`/`K̄₄`
  values are the calibrated moments; `G=GENUINE_BLOCKED`; the joint STAYS CALIBRATED not PASS).
- **stage028 the 2.5PN match-back** — the K̄ CLOSURE OVERLAY is the A3 boundary SHARED with 028. 027's closure slot = the
  port-closure CONSISTENCY (`K̄₄−4K̄₂²/K̄₀=0`, `Γ̄₅−2G/(5c⁵)=0`); 028 = the full INV1–INV5 match-back + the 11-mutation
  coherent-rescale matrix. ⚠ **SHARED, NOT double-counted** — the directive states the 027/028 cut (027 does NOT rebuild
  028's INV5 literal anchors; 028 does NOT rebuild the port derivation).
- **stage005 (`c_s`) + `a` (CONV) + `c` (GR-units bridge, benchmark) + `G` (`GENUINE_BLOCKED`) + `D0`** — the units/closure
  carriers.

## §1d The 027 cut (the checks-vs-siblings boundary — CONFIRM at directive)
The split (part2 L49–52): **024 = the derivation; 025 = the taint; 026 = the lineage; 027 = "the 6 able-to-fail checks
(`[N0]=L⁻¹M`, a⁻⁵, `+i z⁵/27`, χ_Q=1, …) + `P0_phys=(c_s/a)²N0/D0` + the K̄ closure slot (the A3 boundary — marked shared,
not double-counted)".** Concretely:
- **027 owns:** the dim engine (`dim_of`/`scale_power`/`BASE_DIMS`/`BASE_A_POWERS`), `dtn_sign`, `closure_overlay`,
  `evaluate` (the verdict chain), `assert_gate`, `controls`; the 5 NEW port checks (`dim_ok`, `scaling_ok`, `sign_ok`,
  `nonzero_ok`, `closure["ok"]`) + the `deferred_scalar` inconclusive branch; the joint `origin_ok`/`DENSITY_PORT_HOSTED`
  ASSEMBLY; the EM-scaffold RETIREMENT record. The controls: `zero_coupling`, `dimensional`, `sign`, `scaling`,
  `deferred_scalar`, `deferred_scalar_proven_converse`.
- **024 owns:** `N0_den` (027 CONSUMES it, checkable host-contract). **025 owns:** the taint graph + `vector_free`/
  `source_map_complete` (027 CONSUMES the certificate). **026 owns:** the lineage token-check + `moment_valid`/`I25`/
  `continuity_dependency_ok` (027 CONSUMES the certificate). ⚠ **027 assembles `origin_ok` FROM the 025+026 certificates**
  (the composite, §H3) — it does NOT re-litigate vector-freedom or lineage; a corrupted certificate flows to the joint FAIL.
- **⚠ The a-power `−7/2` seam (026↔027):** 026 owns the `I25`-vs-`I_wrong2` earning (the SYMBOL); 027 owns the SCALING
  VERDICT (`p0_power == −5`). The `scaling` control (a-power `−3`) is 027's — but it also flips 026's moment; 027 tests
  ONLY the scaling consequence (with the dim-compensation subtlety, headline 4), consuming 026's earned moment for the
  `−7/2` baseline. The directive must state this so 027 does not rebuild 026's moment gate and 026 does not rebuild 027's
  scaling verdict.
- **⚠ The A3 closure boundary (027↔028):** 027's `closure_overlay` = the two K̄ residuals (the port-closure consistency,
  the A3 SLOT); 028 = the full 2.5PN match-back (INV1–INV5, the 11-mutation matrix). SHARED, marked-not-double-counted.

## §2 The pathA_43 gate history 027 must carry (the split L93 anti-rig burden + the 020/021 trip-up)

**⭐ The CALIBRATED-not-PASS discipline (020's trip-up, L93; the closure's honest landing).** The joint `DENSITY_PORT_HOSTED`
lands the STRUCTURE (density-hosted, vector-free, continuity-lineaged, dim/scaling/sign/closure well-formed) — the MAGNITUDE
stays CALIBRATED (`G=GENUINE_BLOCKED`; `54/5`=`external_bridge_input`; the `27`=`derived_in_gate` from 018). The closure
`Γ̄₅=2G/(5c⁵)` is a CONSISTENCY over the calibrated moments, NOT a first-principles `Γ̄₅`/`G` derivation (the A3 scope,
SHARED with 028). ⚠ 027 must NOT dress the closure as an EARNED `Γ̄₅`/`G` — it is the reduced-closure consistency; the
directive carries the 020 `G=GENUINE_BLOCKED` + the `external_bridge_input` provenance.

**⭐ The joint-assembly anti-rig (the retirement is CONDITIONAL).** R43/R44 DEFERRED the EM `A_w`/`U,W` scaffold RETIREMENT
to this joint. 027's directive MUST make the retirement CONDITIONAL on the joint landing — a `FAIL_*` (dim/scaling/sign/
closure) or `PORT_INCONCLUSIVE_SIM_DEFERRED` must NOT record the retirement (the port has not displaced the scaffold; the
diagnostic sliver reopens). The gate must be able to emit those (the `dimensional`/`sign`/`scaling`/`zero_coupling`/
`deferred_scalar` controls prove it able-to-fail).

**⚠ The sibling-certificate consumption anti-rig (do NOT rebuild 024/025/026).** 027 CONSUMES the three certificates. The
consumption-integrity oracle (`N0_den.free_symbols == HOST_CONTRACT`) ties the assembly to the real port. A corrupted
certificate (a forged `vector_free=True` over a vector-hosted `N0_den`, a forged `moment_valid=True` over a mis-lineaged
port) must flow to `FAIL_NOT_DENSITY_DERIVED` — 027's joint-assembly tooth (§5). ⚠ Do NOT enforce the consumption with a
source grep ([[feedback-grep-acceptance-dodgeable]]) — the real property is the COMPUTED assembly over the real objects.

## §3 Reshape cost (NO bridge to sever) + the 027-JOINT verdict (027 OWNS it)

**⭐ CONTRACT-CLEAN (shared with 024/025/026; part2 L78–79).** BOTH engines already standalone print-only, zero file I/O. NO
`argparse`/JSON/scratch-YAML/`Export` to strip. The reshape = **DECOMPOSITION + the `.wl` per-function determination**:
extract the dim engine + `scale_power` + `dtn_sign` + `closure_overlay` + `evaluate`/`assert_gate`/`controls` slice into a
self-contained script that (1) CONSUMES 024's `N0_den` (checkable host-contract) + 025's vector-free certificate + 026's
`moment_valid`/`I25`/lineage certificate, (2) runs the 6 able-to-fail checks (`dim_ok`, `scaling_ok`, `sign_ok`,
`nonzero_ok`, the `deferred_scalar` inconclusive branch) + the K̄ closure (both residuals, each a standalone assert), (3)
ASSEMBLES the joint `origin_ok`/`DENSITY_PORT_HOSTED` + RECORDS the EM-scaffold RETIREMENT (conditional on the joint
landing), (4) runs the per-tooth control battery (§5), (5) DROPS the 024 derivation / 025 taint / 026 lineage machinery
(cite the certificates) + the `.py`-only `non_dodge` table + the `powers[Xi_Q]=0` no-op, (6) applies the `.wl` per-function
transliteration screen (headline 5 / §1b), (7) uses LOCAL-ledger idioms (`banner`/`subbanner`/`expect_zero`/`expect_bool`/
`RigAssertion`/`routed_assert`/`exercise_rig`/tally/`OVERALL PASS`/nonzero exit for `.py`; `fail[]`/`Exit[1]` for `.wl`)
instead of the pathA_43 `assert_gate` monolith. Arity discipline (def/call scan + unevaluated-leakage scan).

**⭐ THE 027-JOINT VERDICT (027 OWNS `DENSITY_PORT_HOSTED` — the 023-pattern, NOT a printed PARTIAL).** 027 is the COMPLETING
leg, so — UNLIKE 024/025/026 which printed a JOINT PARTIAL — **027 emits the joint `DENSITY_PORT_HOSTED` as its OWN
verdict** (the exit-0 gate; the 023-pattern where the completing leg owns the joint). It emits:
- the **JOINT verdict** `DENSITY_PORT_HOSTED` (the exit-0 gate; the assembled `origin_ok` ∧ `vector_independence_ok` ∧
  `nonzero_ok` ∧ `dim_ok` ∧ `scaling_ok` ∧ `sign_ok` ∧ `closure["ok"]`), with the honest scope printed: **CALIBRATED not
  PASS** (`G=GENUINE_BLOCKED`, magnitude SIM_DEFERRED; the STRUCTURE is what is hosted);
- a **completion statement**: the pathA_43 4-way split is COMPLETE (024 derivation ∧ 025 vector-freedom ∧ 026 lineage ∧ 027
  checks+closure), the EM `A_w`/`U,W` scaffold RETIRES (recorded with the completed joint), the pathA_43 diagnostic sliver
  CLOSES. ⚠ CONDITIONAL: only if the joint LANDS `DENSITY_PORT_HOSTED` (a FAIL/inconclusive → no retirement).

**Acceptance (dual-engine, both exit 0, CWD-independent):**
- Each engine runs from repo root AND a foreign CWD (`/tmp`), print-only, exit 0, no files written.
- Both transcripts print: the consumed `N0_den` (factored) + its `free_symbols == HOST_CONTRACT`; the consumed 025 + 026
  certificates; `[N0_den]=L⁻¹M` (the dim gate, units-restored); `P0_physical=(c_s/a)²N0_den/D0` + `p0_power=−5` (the scaling
  gate); the DtN sign `coeff/i=1/27` outgoing (`χ_Q=+1`); the nonzero-port; the closure residuals `K̄₄−4K̄₂²/K̄₀=0` +
  `Γ̄₅−2G/(5c⁵)=0`; the ASSEMBLED joint `DENSITY_PORT_HOSTED` (CALIBRATED); the retirement record; the per-control outcomes.
- Dual-engine agreement is transcript-level (both print the same checks + the same joint verdict + per-control verdicts);
  neither reads the other; the `.wl` passes the per-function transliteration screen (§1b).
- **Every check + control fires at its own named assert** (per-tooth ablation, §5); the `deferred_scalar` branch emits the
  honest `PORT_INCONCLUSIVE_SIM_DEFERRED` (and the proven converse → `DENSITY_PORT_HOSTED`); a corrupted sibling
  certificate → `FAIL_NOT_DENSITY_DERIVED`.

## §4 Consumed / exported

- **Consumes:**
  - **stage024 `N0_den`** (the density two-port numerator) — the SUBJECT of the 6 checks; cited factored + the checkable
    host-contract (`free_symbols == HOST_CONTRACT`; `rho=rho_eff`).
  - **stage025 the vector-freedom certificate** (`vector_free`/`source_map_complete`/`vector_port ∉ taint`) — CONSUMED as
    booleans → `origin_ok`/`vector_independence_ok`.
  - **stage026 the continuity-lineage certificate** (`moment_valid=True`/validated `I25`/`continuity_dependency_ok`) —
    CONSUMED → `origin_ok`.
  - **stage018 the DtN fingerprint** (`+i z⁵/27`, `χ_Q=+1` outgoing / `−1` incoming) — the `dtn_sign` gate. PROVENANCE +
    the sign check.
  - **stage021 the μ̂₀-free dim machinery** (`[P0_phys]=(c_s/a)²N₀/D₀=1`, a-power `−5`; `BASE_DIMS`/`scale_power`) — the
    dim/scaling gate. CITE (do NOT re-derive the closure logic).
  - **stage020 the `54/5=2·27/5` + `G=GENUINE_BLOCKED`** — the closure `Γ̄₅=2G/(5c⁵)` calibrated magnitude. CITE.
  - **stage005 (`c_s`) + `a` (CONV) + `c` (GR bridge benchmark) + `G` (`GENUINE_BLOCKED`) + `D0`** — units/closure carriers.
- **Exports (→ 028 + Part VII):** ⭐ **the completed joint `DENSITY_PORT_HOSTED`** (CALIBRATED) + the EM-scaffold RETIREMENT
  record (the diagnostic sliver CLOSED) + the K̄ moments `{K̄₀,K̄₂,K̄₄,Γ̄₅}` (→ 028's 2.5PN match-back, the A3 boundary —
  SHARED, not double-counted; per cross-stage flows part2 L109: "024–027 export N0_den + K̄ moments → 028"). ⭐ **This is
  where the pathA_43 gate LANDS — the density-native ℓ=2 radiative port is HOSTED (structure earned, vector-free,
  continuity-lineaged, well-formed), the EM vector scaffold RETIRES, the magnitude stays CALIBRATED (Gate-6/sim-deferred).**

## §5 Teeth candidates (027-specific, per-tooth ablation MANDATORY — mutate the named object, confirm exit-1 AT its own assert)

⭐⭐ **The split L93 + the 020/021 trip-ups: the 6 checks + closure are COMPUTED (never typed/stamped); the μ̂₀-free dim gate
+ real two-verdict self-ablations (the v1 back-solved-μ̂₀ rig banned); the joint assembly is able-to-fail (a corrupted
certificate flows to FAIL). Each check + control independently ablatable ([[feedback-per-tooth-ablation]]).**

1. **⭐⭐ The dimensional gate tooth (`dim_ok`, the `dimensional` control).** Corrupt `[N0_den]` (the `corrupt_dimension`
   hook) → `dim_of` fires `DimError` / `[N0_den] ≠ L⁻¹M` → `FAIL_PORT_MALFORMED(dimensional)`. ⚠ The μ̂₀-free discipline
   (021): the gate reads the SOURCED port dims (`[N₀]`, `[D₀]`), NOT a back-solved free carrier; the corrupt-`[N₀]`
   fires, corrupt-`[G]` is NO_FAIL (`G ∉ P0_physical.free_symbols`, a scope diagnostic). (`dimensional` control →
   `FAIL_PORT_MALFORMED(dimensional)`.)
2. **⭐⭐ The a⁻⁵ scaling gate tooth (`scaling_ok`, the `scaling` control).** The `coupling_a_power=−3` hook → `p0_power=−4
   ≠ −5` → `scaling_wrong` → `FAIL_PORT_MALFORMED(scaling)`. ⚠ **The dim-compensation subtlety (headline 4):** the `−3`
   case flips the moment to `I_wrong2` `(2,0,0)` + `a⁻³`, which COMPENSATES `I25` `(5/2,0,0)` + `a^(−7/2)` in the L-slot
   → `dim_ok` STAYS True → ONLY `scaling` fires (single-tag `FAIL_PORT_MALFORMED(scaling)`, NOT `(dimensional,scaling)`).
   The tooth must assert the SINGLE tag (proving the compensation is honored). (`scaling` control.)
3. **⭐⭐ The DtN sign gate tooth (`sign_ok`, the `sign` control).** The `incoming_sign` hook (`h₂=j₂−i·y₂`) → `coeff/i=
   −1/27 ≠ 1/27` → `FAIL_PORT_MALFORMED(sign)` (`χ_Q=−1` incoming). The outgoing baseline → `+1/27` (`χ_Q=+1`). ⚠ A typed
   `χ_Q` would be a tautology — the sign is COMPUTED from `j₂±i·y₂` (018's discipline). (`sign` control.)
4. **⭐ The nonzero-port tooth (`nonzero_ok`, the `zero_coupling` control).** `coupling_zero → g_base=0 → N0_den=0` →
   `nonzero_ok=False` → `FAIL_PORT_VANISHES`. ⚠ Guarded by `vanished_continuity_coupling` (L593–599) so it fails at
   `nonzero_ok` NOT at `origin_ok` (a legit vanishing coupling keeps `origin_ok` True) — the tooth asserts the verdict is
   `FAIL_PORT_VANISHES` (not `FAIL_NOT_DENSITY_DERIVED`). (`zero_coupling` control.)
5. **⭐⭐ The K̄ closure tooth (`closure["ok"]`, BOTH residuals standalone).** Corrupt `K̄₄` → `kbar4_residual ≠ 0` → the
   closure assert fires; corrupt `K̄₀`/`Γ̄₅` → `gamma5_residual ≠ 0` → fires. ⚠ **Add a STANDALONE `Γ̄₅` residual assert**
   (the source gates `Γ̄₅` only through `closure["ok"]`, L783 asserts only `K̄₄`; agent §E-4) — so BOTH residuals are
   independently ablatable. ⚠ The A3 boundary: this is the port-closure CONSISTENCY (SHARED with 028); do NOT rebuild 028's
   INV5 literal anchors.
6. **⭐⭐ The `deferred_scalar` inconclusive branch (BOTH directions — the honest-landing tooth).** `deferred_uncertified`
   → `Xi_deferred` power `None` → `scaling_undecidable` → `PORT_INCONCLUSIVE_SIM_DEFERRED` (neither FAIL nor PASS). The
   converse `proven_deferred` → power `0` → decidable → `DENSITY_PORT_HOSTED`. ⭐ Both must fire (the branch is honest, not
   baked toward hosted — the able-to-PASS-and-able-to-inconclusive reversibility, the stage020 MIXED-control discipline).
   (`deferred_scalar` → `PORT_INCONCLUSIVE_SIM_DEFERRED`; `deferred_scalar_proven_converse` → `DENSITY_PORT_HOSTED`.)
7. **⭐⭐ The joint-assembly tooth (a corrupted sibling certificate → FAIL).** Consume a forged 025 certificate
   (`vector_free=True` over a vector-hosted `N0_den`, or `source_map_complete=False`) or a forged 026 certificate
   (`moment_valid=False`) → `origin_ok=False` → `FAIL_NOT_DENSITY_DERIVED`. ⭐ Proves the joint is ASSEMBLED (able-to-fail
   if a sibling's certificate is corrupted), not a stamped pass. The consumption-integrity oracle (`free_symbols ==
   HOST_CONTRACT`) ties the assembly to the real port.
8. **⭐ The baseline-valid POSITIVE (the anti-over-rejection / reversibility control).** The genuine density port
   (`−7/2` a-power, outgoing sign, valid certificates) → all 6 checks PASS + closure ok → `DENSITY_PORT_HOSTED`
   (CALIBRATED). ⭐ Proves the gate is NOT rigged toward FAIL — a genuine density port is HOSTED (the stage020
   MIXED-control / stage025 able-to-PASS discipline).
9. **⭐ The retirement-is-conditional tooth.** Force a `FAIL_*` (any of the 6 checks) → the EM-scaffold retirement is NOT
   recorded (the joint did not land). ⭐ Proves the retirement is CONDITIONAL on the joint (R43/R44 deferred it here — a
   failed port does not displace the scaffold).
10. **The `.wl` per-function INDEPENDENCE + arity integrity (§1b).** Confirm the `.wl` dim engine + `dtnSign` are
    keep-native-defensible (like 018/021) and `closureOverlay`/`evaluate` pass the transliteration screen; def/call arity
    scan + unevaluated-leakage transcript scan (the `.wl` `Module`s).

⚠ **NOT 027 (do not rebuild — 024/025/026 own):** `schur_density_expression`/the 2×2 inverse (024 — 027 CONSUMES `N0_den`);
the taint graph `source_tag_map`/`taint_for_expr`/`source_graph_for_expr`/`vector_ablated_expr` + `VECTOR_SYMBOLS`/
`BASE_SOURCE_TAGS` (025 — 027 CONSUMES the certificate); `continuity_lineage_valid`/`continuity_moment_symbol`/the
`I25`-vs-`I_wrong2` earning gate + `fake_continuity`/`attack2` (026 — 027 CONSUMES `moment_valid`/`I25`). ⚠ **the moment
EARNING** (026 — 027 owns only the SCALING consequence of the a-power, not the `I25`-vs-`I_wrong2` symbol choice). ⚠
**Vestigial — DROP:** `omega` (dead both engines), the `.py`-only `non_dodge` table (L636–655) + `powers[Xi_Q]=0` no-op
(L550–553), the vector/relabel/free-carrier control fixtures (025's), the lineage controls (026's).

## §6 Register expectation (orchestrator authors; CONFIRM at register + Codex-verify)

- **⭐ ZERO new counted CALIB knobs** (027 is a CHECKS+CLOSURE/proof slice — it consumes 024's `N0_den` + 025/026
  certificates + 018/020/021 provenance, introduces NO new physical symbols). The closure K̄ moments `{K̄₀,K̄₂,K̄₄,Γ̄₅}` are
  functions of `{G, c_s, a, c}` — `G=GENUINE_BLOCKED` (registered), `c_s`=R1, `a`=`CONV`, `c`=the GR-units bridge (benchmark,
  registered at 020); `D0`=the reduced static conservative denominator (dimensional PROVENANCE, from 021). Part-II CALIB set
  stays **6** (`{μ_η,T_w,β}`+`{Vp0/ℓ_c}`+`{T_Ω,β₂}`). Like the EARNED/proof legs 018/020/021/024/025/026 — all ZERO new
  knobs.
- **⭐ New dims to VERIFY (027 runs the dim gate — the CHECK previously PROVENANCE-only at 024):** `[N0_den]=L⁻¹M` (=
  `(−1,1,0)` in the `(L,M,T)` convention) is now dual-engine-CHECKED here (upgrade the register `N0_den` row L187 from "the
  dim CHECK is 027's" → DERIVED/verified at 027); `[P0_physical]=1` (a-power `−5`, via `(c_s/a)²·N0_den/D0`, `[D0]=L⁻¹T⁻²M`);
  the DtN `ŷ₂` dimensionless. ⚠ **Carry the `(L,M,T)`-vs-register-`[L,T,M]` convention note** (the stage016 lesson — a
  convention mismatch is a cross-stage transfer trap; record the check self-contained + note the convention).
- **New structural edge (propose R46):** the density-port CHECKS + CLOSURE + the JOINT landing — `[N0_den]=L⁻¹M` ∧ a⁻⁵
  scaling (`P0_phys` a-power `−5`) ∧ `+i z⁵/27` outgoing DtN sign (`χ_Q=+1`) ∧ nonzero-port ∧ the K̄ closure
  (`K̄₄−4K̄₂²/K̄₀=0`, `Γ̄₅−2G/(5c⁵)=0`) → the joint `DENSITY_PORT_HOSTED` LANDS (CALIBRATED not PASS, `G=GENUINE_BLOCKED`),
  ASSEMBLED from 025's vector-freedom + 026's continuity-lineage certificates + 027's own checks. ⭐ **This is where the EM
  `A_w`/`U,W` vector-scaffold RETIREMENT is RECORDED with the completed joint** (R43/R44 explicitly DEFERRED it here — "the
  retirement is the JOINT 025/027 result; 027 lands it with the completed joint"). Structural (a proof/provenance + joint-
  landing edge, like R43/R44); discharges NOTHING at the KNOB level. ⚠ The K̄ closure is the A3 2.5PN boundary SHARED with
  028 (marked shared, not double-counted). ⚠ The magnitude (`Γ̄₅`/`G`/`54/5`) stays CALIBRATED/SIM_DEFERRED (Gate-6);
  record the retirement as the JOINT result, NOT a magnitude derivation.
- **⚠ CONFIRM at register:** `c`/`D0`/`Xi_deferred`/`I_wrong2` are NOT new knobs — `c`=GR bridge benchmark (registered),
  `D0`=reduced denominator (dimensional PROVENANCE, 021), `Xi_deferred`=the deferred-scalar control (tracked-not-counted,
  like `I_wrong2` the control-scar); the K̄ moments are the calibrated-magnitude functions (`G=GENUINE_BLOCKED`), NOT new
  knobs. No double-count of 018/020/021/024's registered symbols. ⭐ **Update the `N0_den` row (L187) — the dim CHECK now
  LANDS at 027** (was "the dim CHECK is 027's"); ⭐ **record the pathA_43 gate COMPLETE** (024∧025∧026∧027; the EM scaffold
  retired; the diagnostic sliver closed) as a Cluster-C milestone note.

## §7 Deliverables (per the calibrated pipeline)

`ledger_stage027_port_checks_closure_{sympy.py,.wl}` (Codex builds) + self-contained note
`notes/stages/ledger_stage027_port_checks_closure.md` (orchestrator; full inline — the 6 checks, the K̄ closure, the joint
assembly from the 025+026 certificates, the CALIBRATED-not-PASS landing, the EM-scaffold retirement, the a-power seam, the
deferred-scalar branch) + paper card `paper/stages/stage_027.tex` + `\input` (in
`paper/appendices/stage_appendix_part02.tex`) + registration (provenance index / coverage / manifest, count → 27) +
parameter register update (R46, ZERO new knobs, `[N0_den]=L⁻¹M` CHECK lands, the pathA_43 gate COMPLETE + retirement
recorded) + Codex-verify + PDF + commit + docs/memory sync. ⚠ **Cluster C continues 028 (2.5PN match-back — the A3 boundary
027's K̄ closure shares) → 029 (PN DOI-cite) → then the scheduled MIDWAY KNOB AUDIT (Part-II gravity sector CLOSES — the
pathA_40 `Δr=2` codimension dry-run over Parts I–II + the held-out vs irreducible-route-less tally).** ⭐ **027 CLOSES the
pathA_43 density-port fold (024∧025∧026∧027) — the last DERIVATION-class gate of the Part-II gravity sector; 028/029 are the
consistency + citation caps.**
