# II-P2 (ledger_stage025) source map — pathA_43 vector-freedom taint / host-set guard (`DENSITY_PORT_HOSTED` 2/4, the anti-rig proof leg)

> Running-start prep captured 2026-07-10 (post stage024 = pathA_43 II-P1 density-port DERIVATION DONE, commit
> `3676bd24`; Cluster C OPEN), before authoring the II-P2 reshape directive, so the directive can be written without
> re-discovery. **All line refs below VERIFIED against the current sources 2026-07-10** by an orchestrator read of the
> SymPy engine (`pathA_43_density_quadrupole_port_sympy.py` = **860 lines**) covering the taint slice (L179–190,
> L271–302, L335–389) + the verdict machinery (L580–634) + the controls/assert_gate (L725–803), PLUS a fresh-agent
> structural distillation of BOTH engines (report + directive + `.wl`) that supplied the two-rig history, the `.wl`
> slice line ranges, and the ⭐ **keep-native-vs-re-author determination** (below). No line drift (pathA_43 sources
> untouched since the 2026-07-06 commit).
> Companion: `part2_gravity_atomic_split.md` (Cluster C rows **024–027** L49–52 = the 4-way pathA_43 split; the
> **per-gate trip-ups L93** [⚠ "pathA_43 (025/026): the computed-taint + host-set guard + lineage token-check are THE
> earned mechanisms — two caught rigs live in their history — never collapse to name-checks or flags"]; the reshape-cost
> map **L78–79** [pathA_43 = lightest reshape, contract-clean]; the cross-stage flows L104–109; the ▶ NEXT stage025
> entry L566–572) and the **stage024 pair + source map + directive** (the IMMEDIATE predecessor + SIBLING — 024 EXPORTS
> `N0_den`, which 025 consumes; 024's LOCAL-verdict + JOINT-PARTIAL discipline is 025's template) + **stage020/022/023**
> (the freshest RE-AUTHORED-`.wl` exemplars — 025's `.wl` re-authors like these, UNLIKE 024's keep-native) + **stage019**
> (the ⚠ pathA_41 string-concat anti-pattern lesson — do NOT hand Codex a grep-based acceptance; enforce the property
> with a COMPUTED runtime guard, [[feedback-grep-acceptance-dodgeable]]). Build-order id **025**, Part II, Cluster C.
> Source top-line (verbatim, report `reports/pathA_43_density_quadrupole_port.md:1`): **`DENSITY_PORT_HOSTED`** — the
> JOINT verdict of the whole pathA_43 gate; **025 DELIVERS the VECTOR-FREEDOM half (the 2/4 PARTIAL)** — it PROVES the
> density port `N0_den` (exported by 024) is computationally vector-free, retiring the EM `A_w`/`U,W` scaffold. Proposed
> target stem: `ledger_stage025_vector_freedom_taint` (confirm slug at directive authoring).

> **⚠⚠ ERRATA — SUPERSEDED IN PART BY THE DIRECTIVE v2 (the Codex→Grok→Codex bookend, 2026-07-10). The
> `_scratch/ledger_stage025_reshape_directive.md` v2 is AUTHORITATIVE where they differ.** Four framings below were
> corrected by the bookend (Codex 5 BLOCKING, Grok `DIRECTIVE_CLEAN`): **(1)** the expression-level ablation (§2 headline 2(iii),
> §5) is LOGICALLY SUBSUMED by the taint-set gate → it is a RETAINED redundant WITNESS, NOT an independent/decisive proof leg;
> the DECISIVE gate is the COMPUTED taint-set identity + `source_map_complete`; the relabel/hidden rigs are caught by the
> TAINT-SET gate (NOT by ablating the full `P²/Δ²` rational, whose ablation is singular → `nan`). **(2)** Split
> `baseline_ancestry_ok` (taint == EXACTLY the 3 density tags, BASELINE only) from the general `vector_free` predicate — the
> properly-tagged free-carrier control has a legitimate 4th tag, so it passes the general predicate NOT the exact-3-tag
> identity. **(3)** The consumption (§1c) CITES the factored `N0_den` + asserts `free_symbols == the exported host-set
> contract` — do NOT reconstruct via the 2×2 inverse (that dual-site is dropped — over-engineering, blurs the cut). **(4)**
> R44 records ONLY the vector-freedom CONJUNCT; the vector-scaffold RETIREMENT is recorded at 027 with the completed joint.
> Read the directive v2 §2/§5/§6 for the corrected mechanics; the line refs + slice boundaries below remain accurate.

## ⭐ The FIVE headline points (READ FIRST)

1. **⚠ 025 is the ANTI-RIG PROOF leg of the pathA_43 4-way split (024/025/026/027).** 024 DERIVED the density two-port
   numerator `N0_den` and EXHIBITED a density-only construction. 025 **PROVES** that construction is genuinely
   vector-free — that `N0_den` cannot be a disguised copy of the retired EM vector port. This is where the **EM
   `A_w`/`U,W` vector-scaffold RETIREMENT is COMPUTATIONALLY PROVEN** (R43's note: "the vector-scaffold retirement is
   the JOINT 025/027 result — 025 PROVES vector-freedom — NOT 024's"). 025 is a PARTIAL landing (like 018/016/022/024):
   its script emits a LOCAL audit verdict (the vector-freedom proof holds, computed) + the JOINT `DENSITY_PORT_HOSTED
   (2/4)` as a printed PARTIAL. It does NOT own the joint (lands at 027).
2. **⭐⭐ THE EARNED CONTENT = a COMPUTED ancestry/taint graph over `N0_den.free_symbols`, NOT a name-check or a flag.**
   The proof has three COMPUTED components (all over the consumed `N0_den`), each with an anti-rig burden the split L93
   forbids collapsing:
   - **(i) the taint set** — every free symbol of `N0_den` is mapped through a fixed provenance tag map
     (`BASE_SOURCE_TAGS`) and the tags are UNIONed; the result must be exactly `{continuity_interface, pathA_29_bulk,
     pathA_32_wall}` (the density lineage) with **`vector_port ∉ taint`** (SymPy `taint_for_expr` L351–359 →
     `evaluate` L618, `"vector_port" not in tags`);
   - **(ii) the host-set / `source_map_complete` guard** — every free symbol must be PRESENT in the tag map (no
     `missing_source_symbols`) AND carry a NON-EMPTY tag set (no `empty_source_symbols`); a symbol with no provenance
     cannot ride along (SymPy `source_graph_for_expr` L362–383 → `source_map_complete` L581–584); and the direct
     name-membership `vector_host_symbols = N0_den.free_symbols ∩ VECTOR_SYMBOLS = ∅` (L580);
   - **(iii) the expression-level vector-independence ablation** — substitute EVERY vector-tainted symbol → 0 and assert
     `compact(N0_den − ablated) == 0` (the port is UNCHANGED because it depends on none) (SymPy `vector_ablated_expr`
     L386–389 → `vector_independence_ok` L610–614).
3. **⭐⭐ THE #1 GENUINENESS REQUIREMENT — the taint must be COMPUTED OVER TAGS, not a name-check against `VECTOR_SYMBOLS`,
   and the ablation must have INDEPENDENT teeth (Rig 1 + the hidden-vector case).** The two rigs a name-check would MISS:
   - **the relabel rig** (`relabel_rig` control, `kind="relabel"`): a density-LOOKING two-port built from
     `{omega_wall, omega_rho, r_mix, g_rho, g_qold}` — symbols NOT in `VECTOR_SYMBOLS` but tagged `vector_port` in
     `BASE_SOURCE_TAGS` (L297–301). A `free_symbols ∩ VECTOR_SYMBOLS` name-check PASSES it (none of those names are in
     the set); the COMPUTED taint catches it (`vector_port ∈ taint`) AND the ablation catches it (`vector_ablated_expr`
     zeroes them via `vector_tainted_symbols`, L387 → expr → 0 ≠ expr). **This is the caught Rig 1 (below); the anti-fiat
     teeth.**
   - **the hidden-vector-intermediate rig** (`ablation_isolation` control, `hidden_vector_intermediate=True`):
     `N0_den · σ_hidden` where `sigma_hidden` is tagged `vector_port` (L286) but is NOT in `VECTOR_SYMBOLS` — a
     vector intermediate that survives into `N0_den` even though the host set otherwise "looks" density-only. The
     ablation catches it BECAUSE `vector_ablated_expr` ablates over `VECTOR_SYMBOLS ∪ {tag == vector_port}` (L387), not
     just `VECTOR_SYMBOLS`. **This proves the ablation has independent teeth — the split L93 anti-rig.**
   ⚠ **NEVER collapse (i)/(ii)/(iii) to a `set(names) & FORBIDDEN` check or a self-asserted boolean flag.** 024's tooth D
   (`live_names.isdisjoint(FORBIDDEN_VECTOR_NAMES)`, built stage024 `.py` L219–224) is exactly the LIGHT name-check that
   MISSES the relabel + hidden-vector rigs — 025 is the DEEPER computed proof that 024's tooth D deliberately deferred to
   (024 source map §5 tooth 5: "the COMPUTED taint-graph proof … is 025's").
4. **⭐ THE `.wl` — ⚠ RE-AUTHOR (the taint machinery is a TRANSLITERATION in the source; the OPPOSITE of 024's
   keep-native).** For the numerator DERIVATION lane (024's content), the source `.wl` `densityExpression` is a genuinely
   independent Green-DtN route (kept native at 024). **But the TAINT machinery is 025's content, and there the source
   `.wl` (`vectorSymbols` L154, `baseSourceTags` L162–177, `sourceTagMap` L198–206, `taintForExpr` L208–218,
   `sourceGraphForExpr` L220–235, `vectorAblatedExpr` L237–242, the verdict booleans L389–411) is a near line-for-line
   TRANSLITERATION of the `.py`** (identical symbol set, identical tag-map contents, identical taint algorithm, identical
   substitute-vector→0-and-compare ablation, identical `source_map_complete` predicate, identical verdict ordering —
   fresh-agent determination). So the `MATHEMATICA_MIRROR_POLICY` screen REQUIRES 025's `.wl` to be RE-AUTHORED as a
   genuinely independent taint computation (like stage020/022/023's re-authored `.wl`, UNLIKE 024/018/019/021's
   keep-native). ⭐ One candidate independent route (a suggestion for the directive, NOT a pre-design — Codex owns the
   code): compute vector-independence by DIFFERENTIATION rather than substitution — `∀ v ∈ vector-tainted-symbols:
   ∂(N0_den)/∂v ≡ 0` (an expression is independent of a symbol iff its partial derivative w.r.t. it vanishes
   identically) — a materially different algorithm proving the SAME property; plus a different tag-graph traversal (e.g.
   an explicit reachability walk over a `symbol → tags` graph vs a flat `Variables`+union). ⚠ Directive routes the
   RE-AUTHOR requirement through the Codex→Grok→Codex bookend (the transliteration screen EXPLICITLY applied).
5. **⭐ CONTRACT-CLEAN: there is NO bridge to sever (the lightest reshape class in Part II, shared with 024).** BOTH
   engines are already standalone print-only, zero file I/O (`.py` imports only `sympy`/`dataclasses`/`typing`; `.wl`
   grep for `Export|Put|Write|OpenWrite|>>|Save|Import` = ZERO hits — `results.yaml` is written by an EXTERNAL
   orchestrator, not either engine). So the reshape "sever the bridge" step is a **no-op**. The 025 work is: (a) the
   DECOMPOSITION (extract the taint/host-set/ablation slice into a self-contained script that CONSUMES 024's `N0_den`,
   runs the COMPUTED vector-freedom proof + the rig battery, and drops the 026 lineage-validation + 027
   dim/scaling/sign/closure machinery); (b) the `.wl` RE-AUTHORING (headline 4); (c) the LOCAL-verdict + JOINT-PARTIAL
   framing; (d) the per-tooth-ablation rig battery; (e) the LOCAL-ledger idioms (`banner`/`expect_zero`/tally/nonzero
   exit) replacing the pathA_43 `assert_gate` monolith.

## §1 The 025 slice (line ranges) — the CLEAN CUTS (all VERIFIED 2026-07-10)

pathA_43 is ONE script doing all four buckets; 025 owns ONLY the vector-freedom taint / host-set guard / ablation.
**025 owns the COMPUTED taint graph + host-set guard + expression-level ablation over 024's `N0_den`; it does NOT own the
derivation (024), the continuity-lineage token-check (026), or the dim/scaling/sign/closure checks (027).**

### §1a The SymPy 025 slice (`.py`)
- **⭐ `VECTOR_SYMBOLS` L179–190** — the 10-symbol forbidden-vector set `{A_w, F_muw, Jw, U, W, Omega_U, Omega_W, R_mix,
  g_U, g_W}` (the retired EM vector-port coordinates). Used by the host-set membership (`vector_host_symbols`, L580) and
  the ablation (`vector_ablated_expr`, L388).
- **⭐⭐ `BASE_SOURCE_TAGS` L271–302 — THE PROVENANCE LEDGER (025's core object).** The `symbol → {tags}` map. The density
  symbols carry the density lineage — `c_s`/`rho`/`varpi_phi`→`{pathA_29_bulk}`, `varpi_q`→`{pathA_32_wall}`,
  `Xi_Q`/`eta_q`/`eta_phi`→`{continuity_interface}`, `lambda_c`→`{continuity_interface, pathA_29_bulk, pathA_32_wall}`
  (the interface object joining the two modes), `a`→`{pathA_29_bulk, pathA_32_wall}`. ⚠ The vector + relabel symbols
  ALL carry `{vector_port}` (L286–301, incl. `sigma_hidden` L286 and the relabel-fixture `omega_wall`/…/`g_qold`
  L297–301 — the density-LOOKING-but-vector-tagged rig symbols); `free_carrier`→`set()` (empty, L285).
- **⭐ `source_tag_map(moment_symbol, moment_valid, *, free_carrier_tagged)` L335–348** — clones `BASE_SOURCE_TAGS`,
  optionally tags `free_carrier` as `pathA_34_dimensionless_free_carrier` (L342–343), and assigns the ℓ=2 moment symbol
  `{continuity_interface, pathA_32_wall}` if valid else `set()` (L344–347). ⚠ The `moment_valid` input comes from 026's
  `continuity_moment_symbol` — **025 CONSUMES the `moment_valid` boolean as a typed input (like 024 cited `I25`); the
  lineage VALIDATION is 026's.** 025 owns the tag-MAP construction + the taint over it.
- **⭐ `taint_for_expr(expr, tag_map) L351–359`** — the atomic taint: for each `sym ∈ expr.free_symbols`, if `sym ∉
  tag_map` add to `missing`, else union `tag_map[sym]`; return `(taint, missing)`.
- **⭐ `source_graph_for_expr(expr, tag_map, *, coordinate_hosts) L362–383`** — assembles `{symbol_tags, taint_set,
  missing_source_symbols, empty_source_symbols (in-map but empty tags, L374), vector_host_symbols (= symbols ∩
  VECTOR_SYMBOLS, L375), coordinate_hosts}`.
- **⭐ `vector_ablated_expr(expr, tag_map) L386–389`** — `vector_tainted_symbols = {sym : "vector_port" ∈ tag_map[sym]}`
  (L387, ⚠ over TAGS not just `VECTOR_SYMBOLS`); `ablate = (VECTOR_SYMBOLS ∪ vector_tainted_symbols) ∩ expr.free_symbols`;
  return `compact(expr.subs({sym: 0}))`.
- **⭐ The verdict machinery in `evaluate` L580–614 (025's booleans):**
  - `vector_host_symbols = source_graph["vector_host_symbols"]` (L580);
  - `source_map_complete = not missing_source_symbols and not empty_source_symbols` (L581–584);
  - `origin_ok` L600–609 — `("continuity_interface" ∈ tags) and ("vector_port" ∉ tags) and (not vector_host_symbols)
    and source_map_complete and continuity_dependency_ok` (OR the vanished-coupling branch). ⚠ `continuity_dependency_ok`
    (L585–592) is **026's** (lineage_valid + moment_valid + moment ∈ free_symbols); 025's slice of `origin_ok` is the
    vector-freedom + source-map-complete conjuncts. Directive to DECIDE how much of `origin_ok` 025 asserts vs 026 (see
    §"025/026 cut" below).
  - `vector_independence_ok = (not vector_host_symbols) and ("vector_port" ∉ tags) and (compact(expr − ablated) == 0)`
    (L610–614) — **025's headline gate.**
- **Symbol declarations 025 USES (of L24–45):** the density set `a, c_s, rho, I25, Xi_Q, eta_q, eta_phi, varpi_q,
  varpi_phi, lambda_c, q2, Phi2` (as the host-set of the consumed `N0_den`) PLUS — for the rig battery — the vector
  scaffold `Omega_U, Omega_W, R_mix, g_U, g_W, A_w, F_muw, Jw, U, W` (L34–37), the relabel-fixture `omega_wall,
  omega_rho, r_mix, g_rho, g_qold` (L38–40), `sigma_hidden` (L41), `free_carrier` (L43). ⚠ Unlike 024 (which DROPPED the
  vector/relabel symbols), **025 KEEPS them** — they are the rig fixtures 025's proof must catch. `Xi_deferred`/
  `I_wrong2`/`omega`/`z` are NOT 025 (027/026/dead).
- **⭐ The mutation hooks 025's rig battery drives (in `schur_density_expression` / `vector_expression` / `derive`):**
  `vector_injected_density` (`g_base·Ω_U/Ω_W`, L428–429), `hidden_vector_intermediate` (`n0·σ_hidden`, L444–445),
  `free_carrier_rider` (`expr·free_carrier`, L498–499), `free_carrier_tagged` (L488), `kind="vector"`/`"relabel"`
  (`vector_expression` L468–479). These are the rig SOURCES; 025 re-expresses them as its able-to-fail teeth (§5).
- **Shared helpers 025 uses (NOT cut boundaries):** `compact` L52–57, `hstr` L60–63, `bstr` L66–67. ⚠ **The dim engine
  `dim_of`/`scale_power`/`BASE_DIMS`/`BASE_A_POWERS`/`dtn_sign`/`closure_overlay` (L70–158, 200–268, 521–535, 701–722)
  are 027's** — 025 runs NO dim/scaling/sign/closure gate.
- **⭐ CLEAN CUT (SymPy) — 025 owns L179–190 + L271–302 + L335–389 + the L580–584/L610–614 booleans + the rig battery
  driven from `derive`/`vector_expression`, plus the minimal driver + emit.** It touches NONE of:
  `schur_density_expression` L420–465 (**024's** — but 025 CONSUMES its `N0_den`, §1c); `continuity_lineage_valid`/
  `lineage_for`/`continuity_moment_symbol`/`contains_all`/`CONTINUITY_*` L192–197, 305–332, 392–417 + `continuity_dependency_ok`
  L585–592 (**026's** — 025 CITES `moment_valid` as a typed input); `dtn_sign`/the dim engine/`scale_power`/`closure_overlay`
  (**027's**). ⚠ The mixed driver `derive` L482–518 is the SEAM (it wires 026 lineage + 025 tag-map + 024 schur-expr +
  025 source-graph); 025 keeps ONLY the tag-map + source-graph + ablation wiring, CITES 024's `N0_den` and 026's
  `moment_valid`.

### §1b The Mathematica 025 slice (`.wl`) — ⚠ RE-AUTHOR (the source `.wl` taint machinery is a TRANSLITERATION)
- **The source `.wl` taint slice (to be RE-AUTHORED, NOT mirrored):** `vectorSymbols` **L154**, `baseSourceTags`
  **L162–177**, `sourceTagMap` **L198–206**, `taintForExpr` **L208–218** (`Variables[{expr}]` + tag union + missing),
  `sourceGraphForExpr` **L220–235** (incl. `emptySourceSymbols` L224, `vectorHostSymbols = Intersection[symbols,
  vectorSymbols]`), `vectorAblatedExpr` **L237–242** (`Keys[Select[tagMap, MemberQ[#,"vector_port"]&]]` → subst → 0 →
  `FullSimplify`), and the verdict booleans in `evaluate` **L389–411** (`vectorHostSymbols` L389, `sourceMapComplete`
  L390–393, `originOk` L403–407, `vectorOk` L408–411). ⚠ **These are a near line-for-line transliteration of the `.py`**
  (identical algorithm) → the MIRROR_POLICY REQUIRES 025's `.wl` to re-author them as a genuinely independent taint
  computation (headline 4). Supporting lineage constructs (026's) L155–160, `lineageValidQ` L181–192,
  `continuityMomentSymbol` L194–196, `lineageFor` L244–260 — NOT 025.
- **⭐ RE-AUTHOR direction (candidate, not a mandate — Codex owns the code):** a materially different vector-independence
  test — `∀ v ∈ vector-tainted: D[N0_den, v] === 0` (partial-derivative independence) — plus a different tag-graph
  traversal for the taint set + `source_map_complete`. The dual-engine agreement is transcript-level (both print the
  same taint set, the same `vector_independence` verdict, the same per-rig outcomes); neither reads the other.
- **`.wl` shared/emit helpers:** `fail`/`assertTrue` L14–15, `heading`/`subheading`, `clean` (`Factor@Cancel@Together@
  FullSimplify`), the name-map `externalSymbolNames` + `globalSymbolNames` (built stage024 `.wl` L37–48, 77–80 — the free
  symbol extractor). ⚠ **Arity discipline** (the stage007 lesson): a def/call arity scan + an unevaluated-leakage
  transcript scan (the taint `Module`s + any `Solve`/`FullSimplify` leakage) — 025's `.wl` has several `Module`s.

### §1c The consumption resolution (025 CONSUMES 024's `N0_den` — the checkable derive-vs-typed)
⭐ **025's SUBJECT is 024's exported `N0_den` (the density two-port numerator).** Since there is zero file I/O, 025
reconstructs `N0_den` from the same two-port derivation 024 published (cite 024) so that the taint runs over the REAL
symbol content, NOT a hand-typed symbol list (a typed clean symbol list would be a rig — the free_symbols MUST be
COMPUTED from the actual expression).
- **The consumed object:** `N0_den = response² = I25²·Ξ_Q²·c_s⁴·rho_eff·(η_φ·ϖ_q2 + η_q·λ_c)²/(a⁷·(λ_c² −
  ϖ_Φ2·ϖ_q2)²)` (stage024 report `:9`; built stage024 `.py` `derive_density_two_port` → `make_N0(response)`). Its
  `free_symbols` = the density host-set `{a, c_s, rho_eff, I25, Xi_Q, eta_q, eta_phi, varpi_q2, varpi_Phi2, lambda_c}`
  (`q2, Phi2` are the coordinate hosts, projected out of `N0_den`).
- **⭐ The checkable derive-vs-typed (the consuming-stage dual-site discipline, cf. 022/023):** 025 obtains `N0_den` two
  ways and asserts they agree — (site A) the factored closed form cited from 024's export; (site B) reconstructed from
  the 2×2 static-operator inverse (`response = (M⁻¹·source)[Φ₂]`, `N0_den = response²`, the same derivation 024
  published). A one-site corruption (a dropped/renamed symbol) makes the two disagree → fires. ⚠ 025 does NOT re-VERIFY
  the derivation (that is 024's earned content — no double-count); it uses the reconstructed expression as the genuine
  SUBJECT of the taint proof so the `free_symbols` are real. ⭐ Directive to CONFIRM the exact consumption mechanism
  (Codex designs the code); the REQUIREMENT is: `N0_den.free_symbols` is COMPUTED from a genuine reconstruction, tied by
  a dual-site integrity check to the cited 024 form.
- **⚠ `rho` = `rho_eff`** (the reduced-3D MASS density `[M L⁻³]`, NOT stage005's number-density `ρ0` `[L⁻⁴]`) — carry
  the 024 rename consistently (the source symbol is `rho`; the export symbol is `rho_eff`).
- **The `moment_valid` boolean** — 025 CITES it as a typed input (from 026's `continuity_moment_symbol`), like 024 cited
  `I25`. 025 owns the tag-map + taint; 026 owns the lineage VALIDATION that earns `moment_valid=True`.

## §1d The 025/026 cut (the taint-vs-lineage boundary — CONFIRM at directive)
The split (part2 L50–51): **025 = "COMPUTED ancestry/taint over free_symbols + host-set/source_map_complete guard +
expression ablation"; 026 = "ℓ0(M0)→ℓ1(D1)→ℓ2(∫Y2*S_leak) moment ancestry via 29's operator; token-check computed."**
Concretely:
- **025 owns:** `VECTOR_SYMBOLS`, `BASE_SOURCE_TAGS`, `source_tag_map`, `taint_for_expr`, `source_graph_for_expr`,
  `vector_ablated_expr`; the booleans `vector_host_symbols`, `source_map_complete`, `vector_independence_ok`; and the
  taint-SET assertion (`taint == {continuity_interface, pathA_29_bulk, pathA_32_wall}`, `vector_port ∉ taint`). The rigs:
  `relabel_rig`, `ablation_isolation` (hidden_vector), `attack5_vector_injection`, `provenance_less_rider`
  (source_map_complete), `free_carrier_dimension_corruption` (the able-to-PASS reversibility control), and the raw
  `vector_hosted` (kind=vector) baseline-negative.
- **026 owns:** `continuity_lineage_valid`, `lineage_for`, `continuity_moment_symbol`, the `I25`-vs-`I_wrong2` gate,
  `continuity_dependency_ok`, the ℓ0→ℓ1→ℓ2 token-check. The rigs: `fake_continuity`, `attack2_continuity_corruption`.
- **⚠ The shared object `origin_ok` (L600–609)** mixes both (`continuity_interface ∈ tags` + vector-freedom +
  source_map_complete + `continuity_dependency_ok`). **025 asserts the vector-freedom + source-map-complete conjuncts +
  the taint-set identity; 026 asserts `continuity_dependency_ok` (the lineage).** 025 may PRINT the full `origin_ok`
  but its VERDICT-bearing asserts are the vector-freedom set; the lineage-earned half is a typed input (`moment_valid`),
  026 proves it. ⭐ Directive to state this cut explicitly so 025 does not rebuild 026's lineage validator and 026 does
  not rebuild 025's taint graph.
- **⚠ `fake_continuity` is caught by BOTH legs** (its underlying expr is the vector port → 025's ablation fires; its
  lineage is mis-tagged → 026's `lineageValidQ` fires). 025 MAY carry `fake_continuity` as a control its taint machinery
  ALSO catches (via the vector-port ablation), but the DEDICATED lineage-token-check that a `continuity_interface` tag
  is EARNED is 026's — do NOT rebuild `continuity_lineage_valid` in 025.

## §2 The two caught rigs in pathA_43's history (PIN THEM — the split L93 anti-rig burden)

**⭐ Rig 1 — the name-only / relabel ablation rig (025's; caught by Codex xhigh design-review, directive v0→v1
BLOCKER 2).** The ORIGINAL vector-independence check was by symbol NAME (or a post-hoc symbol substitution). A renamed
copy of the old vector 2-port `P_A = Ω_U²·g_W + R_mix·g_U` has density-LOOKING free symbols but is STILL the vector
formula — so a name/substitution check passes it "by fiat" (directive §0 red-flag: "Replacing `U/W/A_w` names with
density labels while preserving the old port formula by fiat"; §3A: "A symbol-name check is NOT sufficient: a renamed
copy of the old vector formula has density-looking free symbols but a `vector_port` ancestor"). **Remediation (the
earned mechanism 025 MUST preserve):** the ancestry/taint machinery — `N0_den` carries a COMPUTED taint set (over the
tag map), and the vector-independence check is a SOURCE-GRAPH ablation ("remove vector EQUATIONS/tainted symbols before
re-deriving, not a post-hoc symbol substitution"). The `relabel_rig` control routes to `FAIL_NOT_DENSITY_DERIVED`.
⭐ **025's directive MUST forbid the name-check regression** ([[feedback-grep-acceptance-dodgeable]] — enforce with a
COMPUTED taint, not a `set(names) & FORBIDDEN`; the relabel + hidden-vector rigs are exactly what a name-check misses).

**Rig 2 — the self-asserted `continuity_interface` tag / re-tag rig (026's; caught by GLM-5.2 tertiary, directive
v1→v2 NIT-1, "the material one").** Origin validated off a SELF-REPORTED tag: a devious executor writes a NEW equation
≡ the vector port `P_A`, tags it `continuity_interface`, wires it as `N0`'s ancestor → the origin check passes, the
ablation removes nothing (no bare vector symbols survive), the relabel control does not fire. **Remediation:** a BINDING
tag-verification lineage — a `continuity_interface` tag is valid ONLY if the equation is a structured ℓ=2 projection of
the SAME projected-continuity operator that produced pathA_29's ℓ=0/1 moments (`M0=∫S_leak d³x`, `D1_i=∫x_i S_leak +
∫j_i`), i.e. an `∫Y₂*S_leak d³x`-type moment; the validator checks structural TOKENS, not a shared `valid` flag.
**⚠ This is 026's remit** (the lineage token-check) — 025 CITES `moment_valid`; but 025's directive should NOTE Rig 2 so
the sibling 026 pins it, and so 025 does not accidentally accept a self-asserted `continuity_interface` tag as a
substitute for the vector-freedom proof.

## §3 Reshape cost (NO bridge to sever) + the 025-LOCAL verdict + JOINT-PARTIAL

**⭐ CONTRACT-CLEAN (shared with 024; part2 L78–79).** BOTH engines already standalone print-only, zero file I/O. NO
`argparse`/JSON/scratch-YAML/`Export` to strip. The reshape = **DECOMPOSITION + `.wl` RE-AUTHOR**: extract the
taint/host-set/ablation slice into a self-contained script that (1) CONSUMES 024's `N0_den` (the checkable
derive-vs-typed, §1c), (2) builds the tag map + COMPUTES the taint set / `source_map_complete` / vector-ablation over
`N0_den`, (3) runs the per-tooth rig battery (§5), (4) EMITS the taint set + the vector-freedom verdict for 026/027, (5)
DROPS the 024 derivation-verification (cite `N0_den`), the 026 lineage validation (cite `moment_valid`), the 027
dim/scaling/sign/closure machinery, (6) RE-AUTHORS the `.wl` taint machinery (headline 4), (7) uses LOCAL-ledger idioms
(`banner`/`subbanner`/`expect_zero`/`expect_bool`/tally/`OVERALL PASS`/nonzero exit for `.py`; `fail[]`/`Exit[1]` for
`.wl`) instead of the pathA_43 `assert_gate` monolith. Arity discipline (def/call scan + unevaluated-leakage scan).

**⭐ THE 025-LOCAL VERDICT + JOINT-PARTIAL (the 018/016/022/024 discipline).** 025 emits BOTH:
- a **LOCAL audit verdict** (the vector-freedom-proof token, exit-0 gate) — propose `DENSITY_PORT_VECTOR_FREE` (confirm
  slug at directive) — asserting the consumed `N0_den`'s COMPUTED taint set is exactly `{continuity_interface,
  pathA_29_bulk, pathA_32_wall}` (no `vector_port`), `source_map_complete` holds, `vector_host_symbols=∅`, the
  expression-level ablation is invariant (`vector_independence_ok`), AND the rig battery fires (every rig →
  `FAIL_NOT_DENSITY_DERIVED` at its own named assert; the properly-tagged free-carrier control PASSES);
- the **JOINT LANDING (PARTIAL)** printed string: `DENSITY_PORT_HOSTED (2/4, VECTOR-FREEDOM — N0_den is computationally
  vector-free; 024 = derivation, 026 = continuity lineage, 027 = port checks + closure)`. ⚠ 025 does NOT emit the joint
  `DENSITY_PORT_HOSTED` as its OWN verdict (that lands at 027) — it prints it as the PARTIAL (the 018/022/024 pattern).

**Acceptance (dual-engine, both exit 0, CWD-independent):**
- Each engine runs from repo root AND a foreign CWD (`/tmp`), print-only, exit 0, no files written.
- Both transcripts print: the consumed `N0_den` (factored), its computed taint set, the host-set + `source_map_complete`
  status, the vector-ablation delta (= 0 for the baseline), the per-rig outcomes, the LOCAL verdict, the JOINT PARTIAL.
- Dual-engine agreement is transcript-level (both print the same taint set + the same per-rig verdicts); neither reads
  the other; the `.wl` is a genuinely independent taint computation (headline 4).
- **Every rig tooth fires at its own named assert** (per-tooth ablation, §5); the able-to-PASS free-carrier-tagged
  control PASSES (reversibility).

## §4 Consumed / exported

- **Consumes:**
  - **stage024 `N0_den`** (the density two-port numerator) — the SUBJECT of the proof; reconstructed genuinely (§1c),
    free_symbols COMPUTED. The one checkable derive-vs-typed dual-site.
  - **stage026 `moment_valid`** (the `I25` continuity-moment validity boolean) — cited as a typed input (forward ref;
    026 earns it). PROVENANCE.
  - the density symbols' PROVENANCE tags (`BASE_SOURCE_TAGS`) — 025's own ledger, sourced from 024's `physical_relations`
    (ϖ_q2←pathA_32 wall, ϖ_Φ2←pathA_29 bulk, λ_c←continuity) + stage005 (`c_s`) + `a` CONV.
- **Exports (→ 026/027 + Part VII):** the COMPUTED taint set `{continuity_interface, pathA_29_bulk, pathA_32_wall}` +
  the `vector_independence` verdict (`N0_den` is vector-free) + the `source_map_complete` certificate. ⭐ **This is where
  the EM `A_w`/`U,W` vector-scaffold RETIREMENT is PROVEN** (024 EXHIBITED a density construction; 025 PROVES it cannot be
  a disguised vector port). 027 consumes the vector-freedom verdict as one conjunct of the joint `DENSITY_PORT_HOSTED`.

## §5 Teeth candidates (025-specific, per-tooth ablation MANDATORY — mutate the named object, confirm exit-1 AT its own assert)

⭐⭐ **The split L93 burden: the taint + host-set guard + ablation are COMPUTED — NEVER a name-check or a self-asserted
flag. Each rig below is a mutation a name-check would MISS but the computed taint catches.**

1. **⭐⭐ The relabel rig (Rig 1 — 025's CORE anti-fiat tooth).** A density-LOOKING two-port built from the relabel
   fixture `{omega_wall, omega_rho, r_mix, g_rho, g_qold}` (all tagged `vector_port`, L297–301; NOT in `VECTOR_SYMBOLS`).
   **Able-to-fail (fires at its own assert):** the COMPUTED taint contains `vector_port` → the taint-set-identity assert
   fires; AND `vector_ablated_expr` zeroes them → `expr − ablated ≠ 0` → the `vector_independence` assert fires. ⚠ A
   `free_symbols ∩ VECTOR_SYMBOLS` name-check PASSES it → this tooth PROVES the computed-taint-over-tags is load-bearing.
   (`relabel_rig` control → `FAIL_NOT_DENSITY_DERIVED`.)
2. **⭐⭐ The hidden-vector-intermediate rig.** `N0_den · σ_hidden` (`sigma_hidden` tagged `vector_port`, L286; NOT in
   `VECTOR_SYMBOLS`). **Able-to-fail:** `vector_ablated_expr` ablates `σ_hidden` (via `vector_tainted_symbols`, L387) →
   `expr − ablated ≠ 0` → the `vector_independence` assert fires; AND `vector_port ∈ taint`. ⚠ Proves the ablation has
   INDEPENDENT teeth (covers tag-vector symbols, not just `VECTOR_SYMBOLS` names). (`ablation_isolation` control.)
3. **⭐ The vector-injection rig.** `N0_den · Ω_U/Ω_W` (genuine `VECTOR_SYMBOLS`, the `vector_injected_density` hook
   L428–429). **Able-to-fail:** `vector_host_symbols = free_symbols ∩ VECTOR_SYMBOLS ≠ ∅` → the host-set-membership
   assert fires; AND the ablation fires. (`attack5_vector_injection` control.)
4. **⭐ The provenance-less rider rig (the `source_map_complete` / host-set guard).** `N0_den · free_carrier` where
   `free_carrier` has an EMPTY tag set (L285). **Able-to-fail:** `empty_source_symbols ≠ ∅` → `source_map_complete =
   False` → the `source_map_complete` assert fires (a symbol with no provenance cannot ride along). (`provenance_less_rider`
   control → `FAIL_NOT_DENSITY_DERIVED`.)
5. **⭐ The missing-symbol rig (the complementary `source_map_complete` guard).** Inject a symbol NOT present in the tag
   map at all → `missing_source_symbols ≠ ∅` → `source_map_complete = False` → fires. (Confirms both the missing- AND
   empty-tag halves of `source_map_complete` are live.)
6. **⭐ The able-to-PASS reversibility control (NOT a FAIL — the anti-over-rejection tooth).** `N0_den · free_carrier`
   where `free_carrier` is PROPERLY tagged `pathA_34_dimensionless_free_carrier` (the `free_carrier_tagged` hook, L342–343).
   **Able-to-PASS:** `source_map_complete` stays True (non-empty tag), not vector → the vector-freedom + source-map-complete
   asserts PASS. ⭐ This proves the guard is NOT rigged toward FAIL — a legitimately-provenanced dimensionless factor is
   ALLOWED; only provenance-less or vector-tainted factors are rejected (the stage020 MIXED-control / stage041 able-to-PASS
   discipline; the pathA_34 sourced-vs-free lesson). (`free_carrier_dimension_corruption` control → stays
   `DENSITY_PORT_HOSTED` in the full gate; 025 asserts its taint/host-set slice PASSES.)
7. **The baseline vector-port negative (`vector_hosted`, kind=vector).** The raw old vector port `P_A = Ω_U²·g_W +
   R_mix·g_U` (all `vector_port`). **Able-to-fail:** the whole taint is `{vector_port}`, `vector_host_symbols ≠ ∅`,
   ablation → 0 ≠ expr → all three 025 gates fire. The contrast case showing the genuine density `N0_den` PASSES where
   the vector port FAILS.
8. **The `.wl` INDEPENDENT-taint-route + arity integrity (§1b).** Confirm the `.wl` taint machinery is a genuinely
   independent route (NOT a transliteration of the `.py`); def/call arity scan + unevaluated-leakage transcript scan on
   the `.wl` `Module`s.

⚠ **NOT 025 (do not rebuild — 026/027 own):** the continuity lineage `continuity_lineage_valid`/`lineage_for`/the
`I25`-vs-`I_wrong2` gate + `fake_continuity`/`attack2_continuity_corruption` (026 — 025 CITES `moment_valid`; `fake_continuity`
is ALSO caught by 025's ablation but the dedicated lineage token-check is 026's); `dtn_sign`/the dim engine/`scale_power`/
`closure_overlay` + the dim/scaling/sign/054-5 machinery + `zero_coupling`/`dimensional`/`sign`/`scaling`/`deferred_scalar`
controls (027). ⚠ **Vestigial — DROP:** `omega`/`z` (dead in 025's slice), `Xi_deferred` (027), `I_wrong2` (026), the
built-024 name-check tooth D (025 SUPERSEDES it with the computed taint).

## §6 Register expectation (orchestrator authors; CONFIRM at register + Codex-verify)

- **⭐ ZERO new counted CALIB knobs** (025 is a PROOF/taint slice — it consumes 024's `N0_den`, introduces NO new physical
  symbols; the vector/relabel/free-carrier symbols are CONTROL fixtures, tracked-not-counted). Part-II CALIB set stays
  **6** (`{μ_η,T_w,β}`+`{Vp0/ℓ_c}`+`{T_Ω,β₂}`). Like the EARNED/proof legs 018/016/022/024 — all ZERO new knobs.
- **New structural edge (propose R44):** the density-port vector-freedom PROOF — `N0_den`'s COMPUTED taint set is exactly
  `{continuity_interface, pathA_29_bulk, pathA_32_wall}` (no `vector_port`), `source_map_complete`, and the
  expression-level vector-ablation is invariant → **the EM `A_w`/`U,W` vector scaffold RETIRES (PROVEN here, the JOINT
  025/027 result R43 deferred)**. Structural (a proof/provenance edge, like R43/R41/R39); discharges NOTHING. ⚠ **This is
  where the vector-scaffold retirement is recorded** (R43 explicitly deferred it — "the retirement is the JOINT 025/027
  result; 025 PROVES vector-freedom").
- **No new dims to verify** (025 runs no dim gate — the `[N0_den]=L⁻¹M` check is 027's). The `BASE_SOURCE_TAGS` tag map is
  a provenance ledger, not a dimensioned quantity.
- **⚠ CONFIRM at register:** the control-fixture symbols (`Omega_U`/…/`g_W`, `omega_wall`/…/`g_qold`, `sigma_hidden`,
  `free_carrier`) are NOT knobs — they are rig scaffolding used ONLY to exercise the anti-rig teeth (tracked-not-counted,
  like stage023's `q_free`/`eta_null`). No double-count of 024's `ϖ_q2`/`ϖ_Φ2`/`λ_c`/`I25`/`rho_eff` (already registered
  at 024, rows L182–187).

## §7 Deliverables (per the calibrated pipeline)

`ledger_stage025_vector_freedom_taint_{sympy.py,.wl}` (Codex builds) + self-contained note
`notes/stages/ledger_stage025_vector_freedom_taint.md` (orchestrator; full proof inline — the tag map, the computed
taint over `N0_den`, the source_map_complete guard, the expression-level ablation, the rig battery) + paper card
`paper/stages/stage_025.tex` + `\input` (in `paper/appendices/stage_appendix_part02.tex`) + registration (provenance
index / coverage / manifest, count → 25) + parameter register update (R44, ZERO new knobs) + Codex-verify + PDF +
commit + docs/memory sync. ⚠ Cluster C continues 026 (continuity lineage — the SIBLING that earns `moment_valid`) → 027
(checks + closure, COMPLETES the joint `DENSITY_PORT_HOSTED`) → 028 (2.5PN match-back) → 029 (PN DOI-cite) → then the
scheduled MIDWAY KNOB AUDIT (Part-II gravity sector closes).
