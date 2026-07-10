# II-P1 (ledger_stage024) source map — pathA_43 density-native ℓ=2 quadrupole two-port DERIVATION (`DENSITY_PORT_HOSTED` 1/4, the EARNED-FIRST derivation leg)

> Running-start prep captured 2026-07-10 (post stage023 = pathA_34 II-G5b nullspace underdetermination DONE, commit
> `a2e80aa8`; Gate-1–5 gravity ladder CLOSED), before authoring the II-P1 reshape directive, so the directive can be
> written without re-discovery. **All line refs below VERIFIED against the current sources 2026-07-10** by a full
> orchestrator read of the SymPy engine (`pathA_43_density_quadrupole_port_sympy.py` = **860 lines**) + a grep-level
> read of the Mathematica engine (`pathA_43_density_quadrupole_port.wl` = **545 lines**) + the directive
> (`directives/pathA_43_density_quadrupole_port.md` = 358 lines) + the report
> (`reports/pathA_43_density_quadrupole_port.md` = 115 lines), PLUS an independent fresh-agent structural distillation of
> BOTH engines that corroborated every mechanic and the 024/025/026/027 cut boundaries — **no line drift** (pathA_43
> sources untouched since the 2026-07-06 commit).
> Companion: `part2_gravity_atomic_split.md` (Cluster C rows **024–027** L49–52 = the 4-way pathA_43 split + the
> per-gate trip-ups **L93** [pathA_43 025/026 anti-rig] + the reshape-cost map **L78–79** [pathA_43 = lightest reshape]
> + the cross-stage flows L104–109 + the ▶ NEXT stage024 entry L537–544) and the **stage018 pair + source map +
> directive** (the CLOSEST exemplar: 018 is pathA_33's EARNED-FIRST 1/4 leg = an exterior DtN slice that KEPT its
> native `.wl`; 024 is pathA_43's EARNED-FIRST 1/4 leg = an interior two-port derivation that KEEPS its native
> Green-DtN `.wl`) + **stage022/023** (the freshest 2-way EARNED-first→FAIL split; 022's LOCAL-verdict-plus-JOINT-PARTIAL
> discipline is 024's template) + **stage009/010 + 016/017** (the CONSUMED bulk mode + wall mode — 024's operator
> entries `varpi_Phi2`/`varpi_q2` are the ℓ=2 continuations of exactly these). Build-order id **024**, Part II, Cluster C.
> Source top-line (verbatim, report `reports/pathA_43_density_quadrupole_port.md:1`): **`DENSITY_PORT_HOSTED`** — the
> JOINT verdict of the whole pathA_43 gate; **024 DELIVERS the DERIVATION half (the EARNED-FIRST 1/4 PARTIAL)** and
> exports `N0_den` to 025 (taint) / 026 (lineage) / 027 (checks+closure). Proposed target stem:
> `ledger_stage024_density_port_derivation` (confirm slug at directive authoring).

## ⭐ The FIVE headline points (READ FIRST)

1. **⚠ 024 is the EARNED-FIRST DERIVATION leg of the pathA_43 4-way split (024/025/026/027).** pathA_43 answered the
   Phase-A A1 question: is the ℓ=2 quadrupole radiative-port numerator `N₀` DENSITY-native, or was the old EM vector
   (`A_w`/`U,W`) channel load-bearing? Answer = **`DENSITY_PORT_HOSTED`** — a genuine density two-port on `(q₂` wall,
   `Φ₂` bulk-density`)`, the vector scaffold RETIRES. **024 builds the DERIVATION** (the 2×2 static coupling operator →
   invert → the two-port numerator `N0_den`). It is a PARTIAL landing (like 018/016/022): 024's script emits a LOCAL
   audit verdict (the derivation is genuinely computed + both engines agree) + the JOINT `DENSITY_PORT_HOSTED (1/4)` as
   a printed PARTIAL. The **anti-rig burden** — the computed vector-freedom taint (025), the continuity-moment lineage
   token-check (026), the 6 able-to-fail port checks + closure (027) — are the SIBLING legs; 024 owns the physics
   derivation + the two-independent-route agreement, and EXPORTS `N0_den` for them to guard.
2. **⭐⭐ THE EARNED CONTENT = a genuine 2×2 Schur / Green-DtN two-port, NOT a relabel.** The port collapses to the
   ℓ=2 projection of the wall motion coupled to the single bulk Helmholtz mode. `schur_density_expression` (`.py`
   L420–465) assembles the static Lagrangian operator over `(q₂, Φ₂)`
   `M = [[ϖ_q2, −λ_c], [−λ_c, ϖ_Φ2]]` (L437), a source vector `(g_q, g_φ) = g_base·(η_q, η_φ)` with
   `g_base = √ρ·c_s²·I25·Ξ_Q / a^(7/2)` (L423, L433–434), and the Φ₂ response `= (M⁻¹·source)[Φ₂] = P_den/Δ` (L439),
   giving the port numerator **`N0_den = P_den²/Δ²`** (L442) with `Δ = ϖ_q2·ϖ_Φ2 − λ_c²` (L440) and
   `P_den = ϖ_q2·g_φ + λ_c·g_q` (L441). Fully expanded (report `:9`):
   **`N0_den = I25²·Ξ_Q²·c_s⁴·ρ·(η_φ·ϖ_q2 + η_q·λ_c)² / (a⁷·(λ_c² − ϖ_Φ2·ϖ_q2)²)`**. The two-port is picture **(ii)**
   of the directive crux (a genuine `(q₂ wall, Φ₂ bulk)` two-port with a continuity-DERIVED coupling), with picture (i)
   single-channel as a Schur-compressed limit and picture (c) irreducibly-vector-hosted a LIVE first-class negative that
   the derivation, by construction from density variables ONLY, does not realize (025 PROVES it computationally).
3. **⭐⭐ THE KEY GENUINENESS UPGRADE the reshape MUST make (the #1 directive requirement).** In the CURRENT source the
   SymPy full-inverse `response` (`.py` L439) is **DECORATIVE** — `N0_den` is assembled directly from the hand-written
   `P_den` (L441) + `Δ` (L440) as `P_den²/Δ²` (L442), so the matrix inverse is NOT load-bearing (a corrupted inverse
   would not move `N0_den`). By contrast the Mathematica engine's `n0 = phiResponse²` (`.wl` L284) IS load-bearing
   (`phiResponse` is the DtN-eliminated solve). **The reshaped 024 SymPy must make the inverse load-bearing** — build
   `N0_den = response²` from `response = (M⁻¹·source)[Φ₂]` (so mutating any operator entry `ϖ_q2/ϖ_Φ2/λ_c` genuinely
   moves `N0_den`), and keep `P_den²/Δ²` as an INTERNAL factorization cross-check (`assert compact(response − P_den/Δ)
   == 0`). This is what makes the "SymPy inverse ≡ Mathematica Green-DtN elimination" agreement a REAL dual-derivation,
   not two typings of the same closed form. (Structural-distillation flag; confirmed by the orchestrator read of L439–442.)
4. **⭐ THE TWO INDEPENDENT ROUTES (genuinely different algorithms — the `.wl` KEEPS native).** SymPy: form the full 2×2
   matrix, take `static_operator.inv()`, read the Φ₂ component (`.py` L437–439). Mathematica: a genuinely different
   **elimination order** — solve the WALL (row-1) equation `ϖ_q2·q₂ − λ_c·Φ₂ = g_q` for `q₂` first
   (`qRule = q₂ → (λ_c·Φ₂ + g_q)/ϖ_q2`, `.wl` L279), substitute into the BULK (row-2) Φ₂ interface equation
   (`phiEq = ϖ_Φ2·Φ₂ − λ_c·q₂ == g_φ /. qRule`, L280), solve the reduced scalar equation
   (`phiResponse = Φ₂ /. First[Solve[phiEq, Φ₂]]`, L281), square (`n0 = phiResponse²`, L284). WL uses `Solve` +
   substitution (`qRule`/`phiEq`) that SymPy never invokes → **materially different route → keep-native is DEFENSIBLE**
   (unlike stage022/023 which re-authored; like stage018/019/021 which kept native). ⚠ The directive routes this
   keep-native decision through the Codex→Grok→Codex bookend (the transliteration screen EXPLICITLY applied), and once
   the SymPy side derives `N0_den` from `response` (headline 3) BOTH engines are genuine independent derivations.
5. **⭐ CONTRACT-CLEAN: there is NO bridge to sever (the lightest reshape in Part II).** BOTH engines are already
   standalone print-only (`.py` imports only `sympy`/`dataclasses`/`typing`; `.wl` grep for
   `Export|Put|Write|OpenWrite|>>|Save|Import` = ZERO hits — the `results.yaml` is written by an EXTERNAL orchestrator,
   not by either engine). So — UNLIKE the pathA_28/29 JSON-digest + cross-YAML bridges (008–010) and the pathA_30–34
   scratch-YAML payload-mirror (011–023) — the reshape "sever the bridge" step is a **no-op**. The 024 work is: (a) the
   DECOMPOSITION (extract the derivation slice into a self-contained script that emits `N0_den` + provenance for
   025/026/027, dropping the taint/lineage/checks machinery); (b) the headline-3 derivation-genuineness upgrade; (c) the
   LOCAL-verdict + JOINT-PARTIAL framing; (d) the derivation-local able-to-fail teeth; (e) the local-ledger idioms
   (`banner`/`expect_zero`/tally/nonzero-exit) replacing the pathA_43 `assert_gate` monolith.

## §1 The 024 slice (line ranges) — the CLEAN CUTS (all VERIFIED 2026-07-10)

pathA_43 is ONE script doing all four buckets; 024 owns ONLY the derivation. **024 owns the density two-port
construction + the two-route agreement; it does NOT own the taint graph (025), the continuity lineage token-check (026),
or the 6 checks + closure overlay (027).**

### §1a The SymPy 024 slice (`.py`)
- **⭐⭐ `schur_density_expression(config, lineage)` L420–465 — THE DERIVATION (the atomic step; math ends L442).**
  - `xi = Xi_Q` L421 (or `Xi_deferred` under the 027 deferred controls — NOT 024's baseline).
  - `q2_moment, moment_valid = continuity_moment_symbol(config, lineage)` L422 — the ONLY upstream call; returns `I25`
    when the lineage is valid (the `I25`-vs-`I_wrong2` gate is **026's**; 024 CITES `I25` as the ℓ=2 continuity moment
    symbol, §1c).
  - `g_base = √ρ·c_s²·q2_moment·xi / a^(7/2)` L423 — the continuity/interface source coupling amplitude (carries the
    `a^(−7/2)` from the structured ℓ=2 projection). The `coupling_a_power != −7/2` branch L424–427 (→ `/a³`) and the
    `vector_injected_density` branch L428–429 (`·Ω_U/Ω_W`) and `coupling_zero` L430–431 (`g_base=0`) are the mutation
    hooks (§5; verdicts decided in 025/027 but the mutations live HERE in the derivation body).
  - `g_q_den = g_base·η_q` L433, `g_phi_den = g_base·η_φ` L434 — the two-coordinate source vector.
  - `static_operator = Matrix([[ϖ_q2, −λ_c], [−λ_c, ϖ_Φ2]])` L437 — the 2×2 static Lagrangian operator.
  - `source_vector = Matrix([g_q_den, g_phi_den])` L438.
  - `response = (static_operator.inv()·source_vector)[1]` L439 — the Φ₂ component of the full-inverse response
    (**= `P_den/Δ`; currently DECORATIVE — headline 3: the reshape makes `N0_den = response²`**).
  - `delta = ϖ_q2·ϖ_Φ2 − λ_c²` L440 (the determinant); `p_den = ϖ_q2·g_phi_den + λ_c·g_q_den` L441 (the adjugate·source
    numerator); `n0 = p_den²/delta²` L442 (**= `response²`**).
  - `hidden_vector_intermediate` L444–445 (`n0·σ_hidden`) — the `ablation_isolation` control hook (025-verdict).
  - `trace` L447–464 — method/`static_operator`/`source_vector`/`Phi2_response`/`Delta_den`/`P_den`/`g_base`/`g_q_den`/
    `g_phi_den`/`continuity_moment_symbol`/`continuity_moment_valid` + ⭐ `physical_relations` L459–463 (the PROVENANCE:
    `varpi_q2 = K₂/M₂ = (c_s/a)²·κ_q` from **pathA_32 wall** ℓ=2 angular operator; `varpi_Phi2 = (c_s/a)²·(6 + (ma)²)`
    from **pathA_29 bulk** Helmholtz ℓ=2 mode at `c_s` [`6 = ℓ(ℓ+1)` for ℓ=2]; `lambda_c = (c_s/a)²·λ̂_Q` from
    **projected continuity/interface matching**). These three relations are the CONSUMPTION of 009/010 (bulk) + 016/017
    (wall) + the continuity operator (§1c).
- **Symbol declarations 024 USES (of L24–45):** `a, c_s, rho, I25, Xi_Q, eta_q, eta_phi` (L25–27), `varpi_q, varpi_phi,
  lambda_c` (L29–31), `q2, Phi2` (L32). ⚠ **NOT 024 (drop / leave to 025's guard):** the retired EM vector scaffold
  `A_w, F_muw, Jw, U, W, Omega_U, Omega_W, R_mix, g_U, g_W` (L34–37 — used ONLY by `vector_expression`/`VECTOR_SYMBOLS`,
  which are 025/control), the relabel-fixture symbols `omega_wall, omega_rho, r_mix, g_rho, g_qold` (L38–40, control-only),
  `sigma_hidden`/`Xi_deferred`/`free_carrier`/`I_wrong2` (L28,41–43, control/fallback), and `omega` (L45 — **dead in both
  engines**, only `z` is used, in the 027 `dtn_sign`).
- **Shared helpers 024 uses (NOT cut boundaries):** `compact` L52–57 (`factor(cancel(simplify()))`), `hstr` L60–63,
  `bstr` L66–67. **⚠ The dim engine `dim_of`/`dim_add/sub/scale` L70–108 + `scale_power` L126–158 + `BASE_DIMS` L200–233
  + `BASE_A_POWERS` L236–268 are 027's** (the dimension + a-scaling CHECKS) — 024 does NOT run a dim/scaling gate (that
  is 027's `[N0_den]=L⁻¹M` / `a⁻⁵` slice); 024 emits `N0_den` in a form 027 will dimension-check.
- **⭐ CLEAN CUT (SymPy) — 024 owns L420–465 (`schur_density_expression`) ONLY, plus the minimal driver wiring to call
  it and emit `N0_den` + trace.** It touches NONE of: `source_tag_map`/`taint_for_expr`/`source_graph_for_expr`/
  `vector_ablated_expr`/`VECTOR_SYMBOLS`/`BASE_SOURCE_TAGS` (L179–190, 271–302, 335–389 = **025**); `continuity_lineage_valid`/
  `lineage_for`/`CONTINUITY_*`/`contains_all` (L192–197, 305–322, 392–417 = **026**, EXCEPT it CITES the moment symbol
  `I25` via `continuity_moment_symbol` L325–332 as a typed input); `dtn_sign`/`dim_of`/`scale_power`/`closure_overlay`/
  `evaluate`/`controls`/`assert_gate` (L521–535, 82–158, 701–803 = **027**). ⚠ **The mixed driver `derive` L482–518 is
  the SEAM** (it wires 026 lineage + 025 tag-map + 024 schur-expr + 025 source-graph); for a self-contained 024 slice,
  keep ONLY the `schur_density_expression` call + its `I25` moment input, DROP the tag-map/source-graph wiring.

### §1b The Mathematica 024 slice (`.wl`) — ⚠ KEEP-NATIVE (a genuinely independent Green-DtN route)
- **⭐⭐ `densityExpression[c_]` L262–300 — THE DERIVATION (math ends L284).** `xi`/`q2Moment`/`momentValid` L264–266;
  `gBase = q2Moment·xi·cs²·√ρ / a^(7/2)` L267 (`/a³` mutation branch L268–271); `gqDen = gBase·etaq`, `gphiDen =
  gBase·etaphi` L274–275; then the **DtN/interface-matching route** (comment L277–278: "solve wall equation for q2
  first, then insert it into the bulk Phi2 interface equation"): `qRule = q2sym → FullSimplify[(lambdac·phi2sym +
  gqDen)/varpiq]` L279 (eliminate `q₂` from the wall row), `phiEq = FullSimplify[varpiphi·phi2sym − lambdac·q2sym ==
  gphiDen /. qRule]` L280 (substitute into the bulk row), `phiResponse = FullSimplify[phi2sym /. First[Solve[phiEq,
  phi2sym]]]` L281 (solve the reduced scalar interface equation), `delta`/`pden` L282–283 (trace only), **`n0 =
  FullSimplify[phiResponse², $Assumptions]` L284 (LOAD-BEARING — built from the DtN-solved response)**; trace L286–298.
- **⭐ INDEPENDENCE VERDICT = KEEP-NATIVE (defensible).** WL's `Solve`+substitution elimination order (`qRule`/`phiEq`)
  is a materially different algorithm than SymPy's symmetric `Matrix.inv()`. Once the SymPy side derives `N0_den` from
  `response` (headline 3), the two engines are genuine independent derivations of the SAME `N0_den` (SymPy full-inverse
  vs WL DtN eliminate-q2). ⚠ Directive MUST still route this through the transliteration screen (Codex→Grok→Codex) and
  confirm the WL `dtnSign`/`dimOf`/etc. that WOULD read as `.py` mirrors are 027's (out of 024's slice), so 024's WL is
  just `densityExpression` + the emit path.
- **Shared/emit helpers (WL):** `fail`/`assertTrue` L14, `exprText`/`pretty` L35–44, `$Assumptions` L46–49, `ClearAll`
  L12. ⚠ **Dead in WL — drop:** `listText` L18 (never called), the FIRST `exprText` L17 (shadowed by L44), `omega` L33.
- **⭐ CLEAN CUT (WL) — 024 owns L262–300 (`densityExpression`) + the minimal driver + emit.** NOT 024: `vectorSymbols`/
  `baseSourceTags`/`sourceTagMap`/`taintForExpr`/`sourceGraphForExpr`/`vectorAblatedExpr` L154–242 (**025**);
  `continuity*`/`lineageValidQ`/`continuityMomentSymbol`/`lineageFor` L155–196, 244–260 (**026**, cite the moment);
  `dtnSign`/`dimOf`/`scalePower`/`closureOverlay`/`evaluate`/`controls`/`assertGate` L64–115, 341–494 (**027**). The
  mixed driver `derive` L312–339 is the seam (same as SymPy).

## §1c The consumption resolution (the checkable derive-vs-typed vs the provenance-only cites)

⭐ **024 is a DERIVATION leg; its consumptions are mostly PROVENANCE cites (the operator entries' physical origins),
with the two-route agreement being the load-bearing internal check.** There is no cross-STAGE dual-site check that fires
(unlike 017/023) — the checkable content is INTRA-stage (the two independent routes + the inverse-is-load-bearing
factorization cross-check). The consumed objects:

- **⭐ pathA_29 bulk ℓ=2 Helmholtz mode → `varpi_Phi2 = (c_s/a)²·(6 + (ma)²)`** (`physical_relations` `.py` L461). This
  is the ℓ=2 continuation of stage009/010's bulk-density Helmholtz mode (`∂_w²Φ + (ω/c_s)²Φ = 0`, 3D radial
  `g''+(2/r)g'+κ²g`, `κ²=(ω/c_s)²−m²`, the `m=0` zero-mode → `p=2`; the `6 = ℓ(ℓ+1)` angular eigenvalue is stage016's
  `λ_m=6` covariance). Cite **stage009/010** (the pathA_29 bulk mode) + **stage016** (the `ℓ(ℓ+1)=6`). PROVENANCE
  (the operator entry `ϖ_Φ2` is an ABSTRACT symbol in 024's derivation, tagged `pathA_29_bulk`; its literal value is
  the cited provenance, not consumed numerically — the two-port structure holds for any `ϖ_Φ2`).
- **⭐ pathA_32 wall ℓ=2 angular operator → `varpi_q2 = K₂/M₂ = (c_s/a)²·κ_q`** (`.py` L460), with the wall angular
  stiffness `K₂ = ∫[T_w β₂'² + (K_η+6·T_Ω)β₂²]` = stage017's grouped-P₂ ℓ=2 port kernel + the `K_η+6·T_Ω` combination
  (the `6` = stage016's covariance eigenvalue). Cite **stage016/017** (the pathA_32 wall mode). PROVENANCE (`ϖ_q2`
  abstract, tagged `pathA_32_wall`). ⭐ **This is a downstream pin of 017's tracked-downstream port scalars** (register
  row 180: 017's `{B̃,Z̃}` port scalars were "tracked/downstream-pinned"; 024 is one of the pins).
- **⭐ The projected-continuity / interface coupling → `lambda_c = (c_s/a)²·λ̂_Q`** (`.py` L462), the wall→bulk mixing
  DERIVED from continuity/interface matching at `w=0` — the SAME projected-continuity operator that produced pathA_29's
  ℓ=0/1 moments (`M0`, `D1`), now carried to ℓ=2 via the `∫Y₂*·S_leak d³x` moment. ⚠ **The LINEAGE of this (the
  ℓ0→ℓ1→ℓ2 token-check that `I25` is a genuine continuity moment, not a vector relabel) is 026's** — 024 CITES the
  continuity/interface origin as PROVENANCE and uses `I25`/`λ_c` as symbols with their tags. The mixing symbol `λ_c`
  carries ALL THREE tags `{continuity_interface, pathA_29_bulk, pathA_32_wall}` (`.py` L284) — it is the interface
  object joining the two modes.
- **`I25` (the ℓ=2 continuity moment symbol, `[I25]=L^(5/2)`) + `Xi_Q` (dimensionless branch scalar).** 024 CITES
  `I25` as the ℓ=2 projected-continuity moment (lineage validated in **stage026**; a FORWARD reference — 026 is the
  sibling leg that PROVES `I25`'s `∫Y₂*S_leak` lineage; 024 uses it as a typed input, like 018 used the deferred port
  kernel `N_n/D_n`). `Xi_Q`, `eta_q`, `eta_phi`, `lambda_c`'s literal value = **SIM_DEFERRED** (report `:57` — "the
  literal value of the `Y2*` moment, `Xi_Q`, `eta_q`, `eta_phi`, `lambda_c` remains SIM_DEFERRED; what is hosted here is
  the reduced coupling STRUCTURE"). 024 derives the STRUCTURE, not the magnitude.
- **`c_s`/`a`** — R1 units carriers (stage005 `c_s²=5Kρ⁴/m`; `a`=`CONV` pin). ⚠ **`rho` = `rho_eff`, an effective
  reduced-3D MASS density `[M L⁻³]`, NOT stage005's `ρ0` (a 4D NUMBER density `[L⁻⁴]`, register L125)** — corrected per
  the directive bookend (Codex BLOCKING 2). It is the mass-normalization density of the reduced `(q₂,Φ₂)` coordinates
  (source docstring "mass-normalized"); provenance = pathA_29 bulk-mode mass-normalization (STRUCTURAL); the literal
  reduction `rho_eff ← {ρ0, m, geometry}` = SIM_DEFERRED/GAP (tracked, NOT a counted CALIB). `c_s` is the density/sound
  speed (the ripple speed = the analog "speed of light" for the GRAVITATIONAL sector — gravity quadrupole radiation
  rides the density mode, NOT an EM vector field; `docs/conceptual_foundation.md` §3/§5).
- **The EARNED exterior fingerprint (`χ_Q=1`, the `+i z⁵/27` outgoing DtN, `P0_physical=(c_s/a)²N₀/D₀`)** — these are
  **027's** consumption (the radiative-sign + dimensional + closure checks), stages 018/021. NOT 024's (024 is the
  INTERIOR two-port derivation; the exterior-wave match is the sibling check leg). ⚠ Do NOT pull 018's `dtn_sign` into
  024.

## §2 The 024 claim-set (derive + assert; report/directive quotes)

- **(a) ⭐ The density two-port derivation (EARNED — the headline; report `:5–14`, `:72`; directive §3A picture (ii)).**
  The 2×2 static Lagrangian operator `M=[[ϖ_q2,−λ_c],[−λ_c,ϖ_Φ2]]` over `(q₂` wall, `Φ₂` bulk-density`)` with source
  `g_base·(η_q,η_φ)`, `g_base=√ρ·c_s²·I25·Ξ_Q/a^(7/2)`, gives Φ₂-response `= (M⁻¹source)[Φ₂] = P_den/Δ` and port numerator
  `N0_den = P_den²/Δ² = I25²·Ξ_Q²·c_s⁴·ρ·(η_φ·ϖ_q2 + η_q·λ_c)²/(a⁷·(λ_c²−ϖ_Φ2·ϖ_q2)²)` (report `:9`). Picture (ii): a
  genuine `(q₂,Φ₂)` two-port with a continuity-DERIVED coupling, built from density variables ONLY (host-set `⊂
  {q₂,Φ₂,c_s,a,ρ,I25,Ξ_Q,η_q,η_φ,ϖ_q2,ϖ_Φ2,λ_c}`, none of `{A_w,F_μw,J^w,U,W,Ω_U,Ω_W,R_mix,g_U,g_W}`). ⭐ COMPUTED via
  `Matrix.inv()` (SymPy) / `Solve`+eliminate (WL), NOT typed. The `port_picture: ii two-port(q2,Phi2)` label (`.py` L811).
- **(b) ⭐ The two independent routes AGREE (EARNED — the dual-derivation).** SymPy full 2×2 inverse `response` and
  Mathematica DtN eliminate-q2 `phiResponse` yield the SAME `N0_den` (`= P_den/Δ` squared) — transcript-level agreement,
  neither engine consuming the other (report `:11–14`: "SymPy uses the full inverse of the real 2×2 static operator.
  Mathematica uses the independent Green-function/DtN route, eliminating q2 before substituting into the Phi2 equation.
  The agreement is symbolic; neither engine consumes the other's numerator or booleans").
- **(c) ⭐ The inverse/elimination is LOAD-BEARING (the reshape's genuineness UPGRADE — headline 3).** After the reshape:
  SymPy builds `N0_den = response²` from `response = (M⁻¹source)[Φ₂]` (so a corrupted operator entry moves `N0_den`), with
  `assert compact(response − P_den/Δ) == 0` as the factorization cross-check; WL's `n0 = phiResponse²` already
  load-bearing. **This is 024's core able-to-fail** (mutate `ϖ_q2/ϖ_Φ2/λ_c` → `N0_den` changes; mutate the inverse
  route → the cross-check fires).
- **(d) The port vanishes iff the coupling vanishes (EARNED — the derivation-genuineness boundary; the `zero_coupling`
  control).** `coupling_zero` → `g_base=0` → `N0_den=0` (`.py` L430–431 / `.wl` L273). 024 owns "the derived port is
  nonzero for nonzero coupling and vanishes iff `g_base→0`" as a derivation property (the standalone `nonzero`→
  `FAIL_PORT_VANISHES` VERDICT check is 027's; 024's tooth is that the vanishing is COMPUTED from the coupling, not
  asserted).
- **(e) The provenance is EXHIBITED (EARNED — the `physical_relations` trace).** `ϖ_q2`←pathA_32 wall, `ϖ_Φ2`←pathA_29
  bulk, `λ_c`←projected continuity/interface (`.py` L459–463; report `Computed Source Graph`). ⚠ The COMPUTED taint-set
  verification (that these tags are genuine, not self-asserted) is **025's**; 024 EXHIBITS the provenance, 025 PROVES it.
- **(f) Scope honesty (report `:57`, `:79–80`; directive §0).** 024 derives the reduced coupling STRUCTURE
  (vector-freedom, a-scaling structure, DtN sign structure) — the literal magnitudes (`I25` value, `Ξ_Q`, `η_q`, `η_φ`,
  `λ_c` throat value) are **SIM_DEFERRED**; `G`, `2/5`, `54/5` are **CALIBRATED** (`G=GENUINE_BLOCKED`). 024 is the
  form/structure derivation, NOT a magnitude derivation. The 054/5 partition + closure are **027's**.

## §3 Reshape cost (NO bridge to sever) + the shared-machinery + the 024-LOCAL verdict

**⭐ CONTRACT-CLEAN (the lightest reshape in Part II; part2 L78–79).** BOTH engines already standalone print-only, zero
file I/O (§ headline 5). So NO `argparse`/JSON/scratch-YAML/`Export` to strip (contrast 008–023). The reshape =
**DECOMPOSITION**: extract the derivation slice into a self-contained script that (1) builds `N0_den` via the
load-bearing inverse/elimination (headline 3), (2) asserts the two-route agreement + the factorization cross-check +
the operator-entry able-to-fail + the coupling-vanishes tooth, (3) emits `N0_den` (factored) + the static operator +
the `physical_relations` provenance + the host-set for 025/026/027 to consume, (4) DROPS the 025 taint-graph, 026
lineage token-check (cite `I25` as a typed continuity-moment input), and 027 dim/scaling/sign/closure machinery, (5)
uses the LOCAL-ledger idioms (`banner`/`subbanner`/`expect_zero`/`expect_bool`/tally/`OVERALL PASS`/nonzero exit for
`.py`; `fail[]`/`Exit[1]` for `.wl` — the stage018/022 template) instead of the pathA_43 `assert_gate` monolith.
Arity discipline (def/call scan + unevaluated-leakage transcript scan; the `.wl` has `Module`s in `densityExpression`).

**⭐ THE 024-LOCAL VERDICT + JOINT-PARTIAL (the 018/016/022 discipline).** 024 is the EARNED-FIRST leg of a 4-way split
that JOINTLY lands `DENSITY_PORT_HOSTED`. So 024 emits BOTH:
- a **LOCAL audit verdict** (the derivation-integrity token, exit-0 gate) — e.g. `DENSITY_TWO_PORT_DERIVED` (proposed;
  confirm at directive) — asserting `N0_den` is genuinely computed from the 2×2 inverse/elimination, both routes agree,
  the operator entries are load-bearing, the coupling-vanishes tooth fires, and the host-set is density-only
  (`⊂` the density set, no vector symbols — a manifest-construction check, NOT the computed taint graph which is 025's);
- the **JOINT LANDING (PARTIAL)** printed string: `DENSITY_PORT_HOSTED (1/4, DERIVATION — the density-native two-port
  N0_den; 025 = vector-freedom taint, 026 = continuity lineage, 027 = port checks + closure)`. ⚠ 024 does NOT emit the
  joint `DENSITY_PORT_HOSTED` as its OWN verdict (that lands at 027, the COMPLETING leg) — 024 prints it as the PARTIAL,
  the 018/016/022 pattern. (Contrast 023, which owned the joint FAIL because it was the FAIL-DELIVERING completing leg;
  024 is the EARNED-FIRST leg, so it's the 018/022 pattern.)

**Acceptance (dual-engine, both exit 0, CWD-independent):**
- Run each engine from the **repo root** AND from a **foreign CWD** (e.g. `/tmp`), both print-only, both exit 0, no
  files written (verify with `find` for new files — though there was never an `Export`, confirm the reshape introduced none).
- Both engines emit the same transcript: `N0_den` (the factored form), `port_picture: ii two-port(q2,Phi2)`, the Φ₂
  response = `P_den/Δ` (both routes), the two-route agreement, the operator/coupling provenance (`physical_relations`),
  the host-set (density-only), the LOCAL verdict `DENSITY_TWO_PORT_DERIVED`, the JOINT PARTIAL `DENSITY_PORT_HOSTED
  (1/4)`.
- **All able-to-fail teeth fire at their own assert** (per-tooth ablation, §5) — the inverse-load-bearing operator-entry
  mutations, the factorization cross-check, the two-route agreement, the coupling-vanishes tooth, the density-only
  host-set check; the `.wl` independent-route + arity.

## §4 Consumed / exported

- **Consumes:**
  - **stage009/010 (pathA_29 bulk mode) + stage016 (ℓ(ℓ+1)=6 covariance)** → `varpi_Phi2 = (c_s/a)²(6+(ma)²)`.
    PROVENANCE (abstract operator entry, tagged `pathA_29_bulk`).
  - **stage016/017 (pathA_32 wall mode)** → `varpi_q2 = (c_s/a)²·κ_q`, `K₂` grouped-P₂ ℓ=2 kernel. PROVENANCE (tagged
    `pathA_32_wall`; ⭐ downstream pin of 017's tracked port scalars, register row 180).
  - **The projected-continuity operator (pathA_29 ℓ=0/1 → ℓ=2)** → `lambda_c = (c_s/a)²·λ̂_Q`; the `I25` moment (lineage
    validated in **stage026** — forward ref, cite typed). PROVENANCE.
  - **stage005 (`c_s²=5Kρ⁴/m`) + `a` CONV** — units carriers. ⚠ `rho`=`rho_eff` (reduced-3D mass density `[M L⁻³]`,
    SIM_DEFERRED/GAP, NOT stage005's `ρ0` number density; §1c corrected).
  - ⚠ NOT consumed by 024: 018's `dtn_sign`/`χ_Q=1`, 021's dim gate, the closure `54/5` — those are **027's**.
- **Exports (→ 025/026/027 + Part VII):** `N0_den` (the density two-port numerator, factored) — the CENTRAL export
  consumed by **025** (its `free_symbols` feed the taint graph + the ablation), **026** (the `I25` moment must appear in
  `N0_den.free_symbols`), **027** (the dimension `[N0_den]=L⁻¹M`, the a-scaling `P0_physical` a-power −5, the radiative
  sign, and the `P0_physical=(c_s/a)²N0_den/D0` + K̄ closure). Also the static operator + `physical_relations`
  provenance + the density-only host-set. ⭐ **This is the "genuinely-new derivation" that RETIRES the EM `A_w`/`U,W`
  vector scaffold + closes the pathA_43 diagnostic sliver** (blueprint §3 Part II A1; the vector scaffold's retirement
  is EXHIBITED here as the density construction, PROVEN vector-free in 025). Per the cross-stage flows (`part2`
  L108–109: "009/010 export the bulk Helmholtz mode … → 024/026. 017 export … the ℓ=2 port kernel … → 024 (wall mode).
  024–027 export N0_den + K̄ moments → 028.").

## §5 Teeth candidates (024-specific, per-tooth ablation MANDATORY — mutate the named object, confirm exit-1 AT its own assert)

1. **⭐⭐ The inverse-is-load-bearing operator-entry teeth (headline 3 — 024's CORE).** `N0_den = response²`,
   `response = (M⁻¹source)[Φ₂]`. Per-tooth: mutate `ϖ_q2` / `ϖ_Φ2` / `λ_c` (an operator entry) → `response` and hence
   `N0_den` CHANGE → the "N0_den equals the derived-from-inverse value" assert / the downstream self-consistency fires.
   ⚠ **The firewall: `N0_den` MUST be built from `response` (the inverse), NOT independently typed as `P_den²/Δ²`** — else
   a corrupted inverse would not move `N0_den` and the derivation is decorative (the current-source weakness the reshape
   fixes). This is the analog of stage023's "genuine `sp.Matrix.rank()`, not zero-padded".
2. **⭐ The factorization cross-check tooth (`response ≡ P_den/Δ`).** `assert compact(response − P_den/Δ) == 0`
   (the full-inverse result equals the adjugate/determinant closed form). Per-tooth: corrupt `P_den` (drop a term) or
   `Δ` (wrong sign on `λ_c²`) → the cross-check fires. This keeps `P_den`/`Δ` honest (they are the trace/emit form) while
   `response` is the load-bearing object.
3. **⭐ The two-independent-route agreement tooth.** SymPy full-inverse `response` ≡ WL DtN-eliminate `phiResponse`
   (transcript-level; neither engine reads the other). Per-tooth (in-engine): each engine's `N0_den` is derived by ITS
   route; a corrupted operator entry moves BOTH; the arbiter confirms the two transcripts print the SAME `N0_den`. ⚠ The
   `.wl` route must stay the genuine `Solve`+eliminate-q2 (not a mirror of the SymPy inverse — the §1b transliteration
   screen).
4. **⭐ The coupling-vanishes tooth (`zero_coupling` → `N0_den=0`).** `g_base→0` → `N0_den→0` (COMPUTED, not asserted).
   Per-tooth: set `coupling_zero` → the derivation yields `N0_den=0` → the "nonzero-for-nonzero-coupling" derivation
   assert reflects it (024's derivation-boundary tooth; the standalone `FAIL_PORT_VANISHES` verdict is 027's). Neuter
   (keep `g_base` nonzero) → the port is nonzero.
5. **⭐ The density-only host-set tooth (manifest construction, NOT the 025 taint graph).** `N0_den.free_symbols ⊂
   {q₂,Φ₂ (coordinate hosts), c_s,a,ρ,I25,Ξ_Q,η_q,η_φ,ϖ_q2,ϖ_Φ2,λ_c}` and `∩ {A_w,F_μw,J^w,U,W,Ω_U,Ω_W,R_mix,g_U,g_W} =
   ∅`. Per-tooth: inject a vector symbol into the derivation (e.g. the `vector_injected_density` `·Ω_U/Ω_W` hook,
   `.py` L428–429) → the host-set intersection with the vector set is nonempty → the manifest check fires. ⚠ This is
   024's LIGHT construction-level check (024 builds from density variables); the COMPUTED taint-graph proof (that even a
   density-LOOKING relabel with a `vector_port` ancestor is caught) is **025's** — 024 does NOT own the taint graph, so
   do NOT rebuild `source_graph_for_expr`/`vector_ablated_expr` here.
6. **The a-scaling structure tooth (light — the STRUCTURE, the CHECK is 027's).** `g_base` carries `a^(−7/2)`; the
   `coupling_a_power=−3` mutation (`.py` L424–427) changes the a-structure. Per-tooth: 024 may assert the derived
   `N0_den` a-power is the structural `−3` (⟹ `P0_physical` a-power `−5` at 027) as a derivation-structure check —
   but the μ̂₀-free `[P0_phys]=1` + a-power −5 VERDICT gate is **027's** (do NOT rebuild `scale_power`/the dim engine in
   024 unless the directive scopes a light structural a-power check to 024; DEFAULT: leave a-scaling to 027, keep 024
   the pure two-port derivation). ⚠ Directive to DECIDE the 024/027 a-scaling boundary.
7. **The `.wl` INDEPENDENT-ROUTE + arity integrity (§1b).** Confirm the `.wl` `densityExpression` is the genuine
   `Solve`+eliminate-q2 route (NOT a mirror of the SymPy `Matrix.inv()`); def/call arity scan + unevaluated-leakage
   transcript scan (the `Module` in `densityExpression`).

⚠ **NOT 024 (do not rebuild — 025/026/027 own):** the computed taint graph `source_tag_map`/`taint_for_expr`/
`source_graph_for_expr`/`vector_ablated_expr` + the vector-independence ablation (025); the continuity lineage
`continuity_lineage_valid`/`lineage_for`/the `I25`-vs-`I_wrong2` gate + the `fake_continuity`/`attack2` controls (026);
`dtn_sign`/the dim engine/`scale_power`/`closure_overlay` + the dimension/scaling/sign checks + the 054/5 closure (027).
⚠ **Vestigial — DROP:** the EM vector scaffold symbols (025's guard, not 024's derivation), `omega` (dead both engines),
the relabel-fixture symbols, `sigma_hidden`/`Xi_deferred`/`free_carrier`/`I_wrong2` (control-only), WL `listText` + the
shadowed first `exprText`.

## §6 Register expectation (orchestrator authors; CONFIRM at register + Codex-verify)

- **⭐ Likely ZERO new counted CALIB knobs** (024 is a DERIVATION/structure slice, like 018/016/022 — the EARNED-first
  legs, all ZERO new knobs). The operator entries `ϖ_q2`/`ϖ_Φ2`/`λ_c` are ABSTRACT symbols whose VALUES are
  SIM_DEFERRED (report `:57`) — NOT calibrated knobs (their literal throat values are Gate-6/sim-deferred, tracked as a
  `GAP`/deferred, not a counted CALIB). `I25`/`Ξ_Q`/`η_q`/`η_φ` are dimensionless/moment structural symbols with
  SIM_DEFERRED magnitudes → tracked-not-counted (deferred branch data). Part-II CALIB set should stay **6**.
- **New structural edge (propose R43):** the density-native ℓ=2 two-port derivation provenance (the `N0_den` structure
  earned; the vector scaffold retired) — discharges NOTHING (a structure/provenance edge, like R37/R39/R41). ⚠ The
  literal port magnitude reduction (the `Ξ_Q`/`λ_c` throat values) is the Gate-6/sim-deferred obligation — record it as
  a `PENDING`/deferred note, NOT a discharge.
- **Dimensions to verify at register (dual-engine, units-restored):** `[N0_den]=L⁻¹M` (checked at 027 but the symbol
  first appears here — record its sourced-dim provenance: `[I25]=L^(5/2)`, `[c_s]=LT⁻¹`, `[ρ]=ML⁻³`, `[a]=L`,
  `[ϖ_q2]=[ϖ_Φ2]=[λ_c]=T⁻²`, `[η_q]=[η_φ]=[Ξ_Q]=1`). ⚠ Do NOT assert `[N0_den]` DERIVED unless the 024 script actually
  runs a dim leg — if 024 leaves the dim gate to 027, record the symbol dims as PROVENANCE and note the check lands at 027.
- **⚠ CONFIRM at register:** whether `varpi_q2`/`varpi_Phi2`/`lambda_c`/`I25`/`Xi_Q`/`eta_q`/`eta_phi` are the FIRST
  appearance of these symbols (they are — pathA_43 is new to the ledger) → add their rows (class = `GAP`/deferred for
  the SIM_DEFERRED magnitudes; the STRUCTURE is EARNED). Distinguish from the 013/017 wall packet `{μ_η,T_w,β,T_Ω,β₂}`
  (those are the ℓ=0/ℓ=2 wall CALIB knobs; `ϖ_q2 = K₂/M₂` is BUILT FROM them + the ℓ=2 geometry, a DERIVED
  manifestation at the value level — but 024 treats `ϖ_q2` as an abstract entry, so the reduction `ϖ_q2 ← {T_w,K_η,T_Ω,β₂}`
  is a provenance note, likely a `DERIVED` edge, CONFIRM).

## §7 Deliverables (per the calibrated pipeline)

`ledger_stage024_density_port_derivation_{sympy.py,.wl}` (Codex builds) + self-contained note
`notes/stages/ledger_stage024_density_port_derivation.md` (orchestrator; full derivation inline — the 2×2 operator,
the inverse/elimination, `N0_den`, the provenance) + paper card `paper/stages/stage_024.tex` + `\input` +
registration (provenance index / coverage / manifest, count → 24) + parameter register update + Codex-verify + PDF +
commit + docs/memory sync. ⚠ Cluster C continues 025 (taint) → 026 (lineage) → 027 (checks+closure, COMPLETES the joint
`DENSITY_PORT_HOSTED`) → 028 (2.5PN match-back) → 029 (PN DOI-cite) → then the scheduled MIDWAY KNOB AUDIT.
