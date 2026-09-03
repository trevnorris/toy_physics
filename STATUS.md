# STATUS — where the Path-A program is (single front door)

> ⚠ **SHA NOTE (2026-09-01 DataLad rewrite):** the git history of `ledger-v3-rebuild` was rewritten — **every commit SHA cited in a clause dated before 2026-09-01 is PRE-rewrite and no longer resolves.** Find those commits by their message, not the SHA. Post-rewrite anchors: #89 directive `660a3055`, #88 `05cb1ea5`, #87 `bab2b828`, #89 engine checkpoint `f655ea65`, #89 clearance `9f40c18e`, **#89a WL basis `d4adbd99`**. `.out` files are now git-annex pointers (`datalad get` to read).

## ⭐⭐⭐ S11c-b CROSS-ENGINE RESIDUAL DEFERRED to a ≥64 GB box (STEP 0 OVERTURNED); S11c-b CLOSES on per-engine leg-verification; NEXT = honest step record (2026-09-03)
⛔ **STEP 0 was WRONG.** It measured only `evaluatedModel` + one `FINAL_KERNEL` (7.95 GB) and concluded the residual fit 30 GB. The full `S11CB_PRIMARIES_ONLY` production emit ALSO builds `mainKernelOrigins` (a FULL `extractCoupling` per origin, ~6/case, for `COUPLING_KERNEL_TERM_ORIGINS`) → ONE case hit **15.6 GB still growing** → OOM. ⇒ **the full cross-engine residual is DEFERRED to a ≥64 GB box** (user decision 2026-09-03: defer whole residual, do not take the lighter core-only path). S11c-b closes on the per-engine leg-verification (fold + #90 + #89/a/b, all 2-leg CLEAR) + the coarse cross-engine consistency established 2026-09-02. The P2a/P2b comparator decision lists are committed as the SPEC for the ≥64 GB run (`DEFERRED_HEAVY_RUNS.md`). ⚠ The two whole-row SIGN conventions (kinetic −K PY vs +K WL; face) + #90's two flags (closure-fold sign; uniform-Λ survivors) remain cross-engine-UNVALIDATED (deferred). NEXT = honest step record `steps/S11c_b_*.md` (2 legs) documenting the arc + the deferrals; exports handling; then S11c-c. Below is the (now superseded) in-flight plan, kept for the artifact map.

User chose the FULL reviewed cross-engine residual, THEN (after the OOM finding) chose to DEFER it. The pass was a P1-WL emit gate + two comparator builders (P2a slab-row/row_residual, P2b §3a bridge) → run the residual — the gate + decision lists are DONE; the builds + run are ≥64 GB work.
- ✅ **P1-WL gate committed `06048d15` + pushed** (origin+GIN synced): additive +17-line `S11CB_PRIMARIES_ONLY` in the WL engine (mirrors PY's), set-mode emits primaries for all 4 cases + skips controls/tower, unset byte-identical to HEAD. **2 build legs CLEAR** (⚠ TWO fresh Claude agents — Grok substituted for the gate's 2nd leg due to an xAI capacity outage; low-risk additive change, objective checks: 14 tags/0 controls, extractors pass, payloads byte-identical set-vs-HEAD, both FORM ablations bite). v1→v2 after 2 decision legs rejected a single-case selector (branch-scoped energy families; `row_residual` raises on non-aligned keys → the fix is the symmetric PRIMARIES_ONLY, both engines emit matching 4-case primaries).
- ✅ **Fresh PY `.out` done** (`~/.s11_build/S11c_b_step0_scope/py_primaries_fresh.out`, 183 MB, folded+#90, PRIMARIES_ONLY). 🔄 **Production WL `.out` running** (`~/.s11_build/S11c_b_production/wl_primaries.out`, task `bgxqpd4ok`, memory-watched, ~8 GB/case as STEP 0 predicted, ~2h). These are the residual inputs (scratch, not committed).
- ✅ **P2a v2 (slab-row join + row_residual) committed `06048d15`; P2b v2 (§3a bridge) committed `69bcc18d`.** Each folded from 2 decision legs (Codex+Grok, ~10 computation-backed defects each). KEY findings: (P2a) `row_residual.py` is in the blast radius — its #90-STALE one-sided face subtraction (`:427`) manufactures a false disagreement; μ_θ not consumed; `DOF` must be pinned; the closed disposition table is now IN the directive; origin-family name-mismatch adapters; the TWO whole-row SIGN conventions (kinetic −K PY vs +K WL; face) are SURFACED in the residual, NOT normalized (adjudicated in the step record — never design away the disagreement). (P2b) the bridge needs an EXPRESSION-valued scale substitution (`gammaWJET*→W_0·gamma_s11cb_*`; `I_PY=W_0·I_WL`), not a string rename; a factor-LOCKED energy-term bijection certificate; the applied→bare guard on the FULL parse (width/modulus/eW), pinned to S11c-b tables (not s11ca FIELD → S11c-a regression). ⛔ Sequencing: P2b lands BEFORE P2a's final validation.
- ⛔ **OOM LESSON (2026-09-03):** ran the WL production run (~8 GB) + the P2b Codex build (parsing the full 156 MB WL + 84 MB PY `.out` into SymPy = tens of GB) CONCURRENTLY → system OOM killed BOTH (watchdog never fired — min avail 8.6 GB; no dmesg access; memory recovered to 23 GB after). ⇒ NEVER run two memory-heavy CAS jobs concurrently on this 30 GB box; SERIALIZE. The incomplete P2b build changes were discarded (committed baseline intact); the P2b directive gained a memory-bound build note (test against a §3a-invariant-only / single-case extract, not the full pipeline).
- ⚠ **exports.py REGEN (uncommitted):** the PY PRIMARIES_ONLY run republished `S11c_b_exports.py` (a hash-chained input; a faithful deterministic regen reflecting the folded+#90 engine — surfaces a PRE-EXISTING staleness the fold/#90 commits deferred). LEFT UNCOMMITTED (regen-able); decide when/whether to commit (it propagates downstream) after the cross-engine residual validates fold+#90.
- **NEXT (serialized, memory-safe):** production `.out` completes → **P2b build** (Codex, bounded, alone) → 2 build legs (fresh Claude+Grok) → commit → **P2a build** (alone) → 2 build legs → commit → **P3: run `row_residual --py <fresh PY> --wl <fresh WL>`** → the cross-engine residual (fold's constraint-reduced rows + #90's closure-fold sign / uniform-Λ survivors; the whole-row sign conventions surface here for step-record adjudication) → #88 KINETIC+θ re-adjudication + 2 owed control-hardenings → honest step record `steps/S11c_b_*.md` + exports (2 legs) → close of S11c-b. ⚠ NO per-substep card (one S11c roll-up after S11c-e); S11c-c/d/e remain. Records: `_measurements/S11c_b_{p1_wl_residual_emit,p2_comparator_update,p2a_slab_row_join,p2b_gamma_bridge}_directive.md`. Post-06048d15 anchors: P1-WL `06048d15`, P2b `69bcc18d`.

## ⭐⭐⭐ S11c-b CROSS-ENGINE STEP 0 (memory scope) DONE — the single-case residual FITS this 30 GB box; integration pass is a 3-builder sequence; NEXT = P1 residual-mode single-case emit (2026-09-03)
The resume-prompt STEP 0 is measured and recorded (`directives/_measurements/S11c_b_step0_residual_scope.md`; `DEFERRED_HEAVY_RUNS.md` reconciled). A guarded single-case WL run of exactly the objects the residual needs — `evaluatedModel` (operator) + `extractCoupling…["FINAL_KERNEL"]` (primary kernel), NO tower-depth control variants — peaked **7.95 GB RSS / 0.99 GB in-kernel, min MemAvailable 14.94 GB, ~26 min** for one case. ⇒ the cross-engine residual is **doable on this box**; the ≥64 GB requirement bounds only the full 4-case in-band `.out` regen with the tower/heavy controls (built unconditionally at `…mathematica_audit.wl:2204-2231`), which the residual does not need.
- ⚠ **REFRAMING (verified):** `row_residual` and the comparator PARSE the two committed `.out` (`--py`/`--wl`, default the committed files) — they do NOT run the engines live. **Both committed `.out` are STALE:** WL at `d4adbd99` (#89a, pre-#89b → FROZEN operator) and PY at the migration checkpoint (PRIMARIES_ONLY, pre-fold, pre-#90). The stale PY `.out` even has the pre-fold `THETA_BALANCE = μ_θ` structure. ⇒ the integration pass MUST regenerate **fresh single-case `.out`** from both current engines first, then point the comparator at them.
- ⚠ **STEP 0 OPERATIONAL DISCOVERY:** the engine's normal emit path always builds the tower variants (~16 GB/case; only `S11CB_SKIP_HEAVY_CONTROLS` gates the 2 equivalence controls, not the variants). A comparator-PARSEABLE single-case `.out` at the measured ~8 GB footprint needs a **residual-mode single-case emit** — restrict to one case, call the engine's OWN `emitShared`/`modelRecord`/`kernelRecord` on the primary objects, skip the tower variants. Same objects, same emit functions.
- ⭐ **BRIDGE DELTA MAPPED (agent) + rule-13 VERIFIED (me):** (a) `extract_slab` has a now-FALSE `THETA_BALANCE → ("MU_THETA","THETA")` bridge — post-fold that row is `evolution_mass_balance − Σ closure_shape_deriv` (mass-evolution), so it would join PY's mass-evolution against WL's raw μ_θ = a MANUFACTURED disagreement; plus `MASS_EVOLUTION_ROW` under-bridged (PY maps only the ADVECTIVE summand) and WL `CENTER_FACE_GENERALIZED_ROW` unpartnered. (b) The §3a coefficient bridge is FULLY STALE: WL renamed its coefficients to `gammaWJET*`/`gammaMUJET*` at `d4adbd99`, AFTER the bridge tables were written (`70164909`), so the 12 `gammaWidth*` entries now match NOTHING ⇒ **0 of ~30 §3a coefficients are currently bridged**; 7 suffixes/source never had an entry. ⛔ PY coefficient names are POSITIONAL (`gamma_s11cb_*_NN`, runtime-quotient-selected) — the pairing MUST be read off the invariant each coefficient multiplies (PY emits `ENERGY_BASIS_NEW_INVARIANTS`), NEVER guessed by index. (c) Profile-jet + DOF bridges are COMPLETE — no gap.
- **NEXT = the integration pass, a 3-builder sequence (each leg-gated):** **P1** add a residual-mode single-case emit to BOTH engines (primary tags only, tower/heavy controls skipped) → **P2** update the comparator (`extract_slab` structural fix + §3a bridge by invariant-matching) — its OWN decision list (drafted, 13 testable properties) + 2 decision legs + 2 build legs → **P3** run `row_residual` on the fresh single-case `.out`, read the residual (OUR finding, no value gate) — validates the fold's constraint-reduced U/E_W rows AND #90's closure-fold sign/uniform-Λ-survivors. Then #88 KINETIC+θ re-adjudication + 2 owed control-hardenings → step record + exports (2 legs) → close of S11c-b. ⚠ NO per-substep card (one S11c roll-up after S11c-e); S11c-c/d/e still remain.

## ⭐⭐⭐ S11c-b #90 §3c COUPLING CONTENT BUILT + leg-gated + committed `7677aa18` (pushed origin+GIN) — the "first physics number" (face+response) is IN; NEXT = cross-engine integration pass (this box) (2026-09-02)
PY's coupling was bulk-only; #90 folds the settled §3c INCLUDE/INCLUDE content — the reversible tilted-face geometry + the irreversible face response (`Λ` symbolic). Record `directives/_measurements/S11c_b_90_coupling_content_directive.md`; directive `directives/S11c_b_90_coupling_content_directive.md`; leg evidence `~/.s11_build/S11c_b_90_coupling/`.
- ⛔ **Decision legs (Codex+Grok) REJECTED v1 architecturally:** `FACE_FLUX_BOUNDARY_OPERANDS` is the RAW T-substrate bundle, not operator rows — weak-restricting it is the §3c-forbidden PARALLEL ROUTE. Fix (both legs, convergent): compute the face GENERALIZED-FORCE rows from the consumed virtual work (coeffs of `δ_vu`/`δ_ve_W`, live `μ_θ` bound), ADD to the constraint-reduced operator rows (origin FACE, not through the θ-fold), then weak-restrict the full operator; fold the SKIPPED `closure_shape_deriv` (Λ_A/Λ_V) + `traction` (Λ_X). Plus the `A_T` token collision (PY's `A_T_s11cb` is the test potential, not the geometry), double-count/over-reach traps, §0 "every" overstatement, rule-5 leaks. v2 folded both (one pass).
- ✅ **Build legs (fresh Claude agent + Grok) BOTH VERDICT CLEAR (8/8 probes, form ablations BITE):** face computed INTO the operator (`U/E_EXPANDED = bulk + face`; no `FACE_FLUX_BLOCKS`; `build_kernel` never reads the raw bundle; Grok's independent face-row reconstruction matches the emit); form-ablating disjoint face sources moves disjoint kernel channels, non-face (projection) moves nothing; `Λ` symbolic (`Λ_I⁰/(1−iωτ_I)`, no `Z`/DtN/impedance); `μ_θ` bound (no reserved `mu_theta_L/M` in blocks); no `ζ_c`; adjointness over the enlarged blocks; exact-once; `operator_from_density`/`committed_strong_rows` byte-identical (#88 refs). ⚠ Codex's OWN self-review (Codex variants, "Claude/Grok unavailable") was DISCARDED as invalid.
- ⚠ **TWO cross-engine / step-record flags (NOT build defects — the construction is correct per §3c; a single engine can't settle the response CONTENT):** (1) the closure-fold sign/magnitude (Grok's T-i identity `TRUE_AREA − RELFLUX_SUM = 0` corroborates the recipe, but sign/magnitude is PY↔WL-only); (2) the §5c uniform-limit residual is NONZERO and Λ-bearing (`(2,4,0,0,0,0)`, γ-count 0) — whether the irreversible face-response legitimately survives at `(η,σ_W)=0` or would violate §1d's uniform decoupling needs the cross-engine residual + step-record (does WL carry the same Λ survivors?). §5c is a smoke test, no value gate.
- **NEXT:** the **cross-engine integration pass**: complete the S11c-b symbol bridge (~35 §3a basis coeffs + jet + DOF transliteration) + update `extract_slab` to the folded PY structure (`THETA_BALANCE`→mass-evolution-minus-closure, `MU_THETA` from `MU_THETA_FACE_BINDING`) → run `row_residual` (reviewed cross-engine JOIN) — this validates the fold's constraint-reduced rows AND #90's closure-fold sign/uniform-survivors → #88 KINETIC+θ re-adjudication + 2 owed control-hardenings → honest step record `steps/S11c_b_*.md` + **exports** (2 legs) → STATUS/memory close **of S11c-b**. ⚠ **The card:** S11c-b owes NO separate card (S11c-a precedent, `steps/S11c_a_*.md:268-269`); N1 (`directives/S11c_decisions.md:24-31`) specifies ONE S11c roll-up card, produced ONLY after S11c-e. ⚠ S11c-b is ONE sub-step — **S11c-c** (curved-bulk closure/DtN), then **S11c-d/e** still remain. ⚠ **RESOURCE SCOPE (Codex-corrected):** the measured **0.9 GB fits only the single-case WL OPERATOR/U-row probe** (`evaluatedModel`), NOT the #90 KERNEL or a full comparator pass — the #90 kernel builds took the legs ~45 min/case, and `DEFERRED_HEAVY_RUNS.md:58-84` still requires ≥64 GB for full WL operator/kernel regen. So the cross-engine pass's feasibility ON THIS BOX is established only for the U-row operator; the kernel/comparator scope must be MEASURED before claiming it, and `DEFERRED_HEAVY_RUNS.md` reconciled.

## ⭐⭐⭐ S11c-b PY CONSTRAINT-FOLD BUILT + leg-gated + committed `82f53828` (pushed origin+GIN) — pin (B) implemented; NEXT = single-case cross-engine `row_residual` (this box) → #90 (2026-09-02)
The PY engine now folds pin (B): the slab `U`/`e_W` rows are the CONSTRAINT-REDUCED equations (the imported non-uniform `virtual_constraint` eliminates virtual θ, the same held-fixed `μ_θ` feeds both reactions), the θ-row is the imported sourced `evolution_mass_balance` (⚠ AT `82f53828`; **#90 SUBSEQUENTLY changed the θ-row to `evolution_mass_balance − Σ closure_shape_deriv`** — folding the face-response closure into the mass row; current HEAD's θ-row is the mass balance MINUS the closure residual, and the U/e_W rows are constraint-reduced internal + face generalized forces), `μ_θ` stays a separate held-fixed operand, and the jet-depth cascade is raised (`STRONG=3, COUPLING=4, DEPTH_CONTROL=5`). `operator_from_density`/`committed_strong_rows` are byte-identical to HEAD (#88 raw reference preserved). Record `directives/_measurements/S11c_b_py_constraint_fold_directive.md`; directive `directives/S11c_b_py_constraint_fold_directive.md`; leg evidence `~/.s11_build/S11c_b_constraint_fold/`.
- ⛔ **The 2 decision legs REJECTED v1 for a rule-17 FREEZE** — I handed §1c:143's `(uniform linearisation)` as the constraint (froze `W_bg`); the binding object is the MATERIAL `δ_vΣ_mat=0` with live `W_bg` (imported as `virtual_constraint`). Plus a θ-row double-count, a #88-break interface, a missing depth cascade, and provenance/leak. v2 folded both legs (rule 7, one pass). The gate caught the cardinal error before one engine line changed.
- ✅ **2 build legs (fresh Claude agent + Grok): PHYSICS ROWS CLEAR on both** — the reaction is COMPUTED from the imported non-uniform constraint (a FORM ablation of the constraint source MOVES the rows; `live−ε∇μ_θ ≠ 0` but `uniform−ε∇μ_θ = 0`, proving the v1 freeze was NOT committed), θ = imported mass-evolution with no double-count, depths load-bearing (cap→2 loses all order-3; coupling cap→3 loses 174 order-4 terms), #88 raw refs intact, provenance `BULK_ENERGY`, no leak, smoke PASS (564s, 4 cases).
- ⚠ **ONE FINDING folded (Claude agent Probe 2, orchestrator CAS-proven rule 13):** the two-route (elimination vs Lagrange) residual is TAUTOLOGICAL — `λ=−ε·μ_θ/∂_θC` makes the routes algebraically identical for any linear constraint, so the residual is `0` by construction even for a WRONG constraint (operand theatre, not independence). Relabeled `ROUTE_RESIDUAL → CONSTRAINT_FOLD_TRANSCRIPTION_RESIDUAL` + the directive claim corrected; the rows are unchanged. The legs DISAGREED (Grok CLEAR via a weaker CODE-corruption test; the agent's WRONG-INPUT test was decisive). Real independence = the cross-engine comparator. **FLAG deferred to #88:** `HESSIAN_FREEZE_STRONG_ROW_RESIDUAL` now compares folded-live vs raw-committed ⇒ is the reaction (nonzero); nothing asserts it.
- ⭐ **CROSS-ENGINE ATTEMPT (2026-09-02) → KEY DISCOVERY: the ≥64 GB box is NOT needed for the residual.** A guarded single-case WL U-row run (EULERIAN/LAB_HELD/RHO4_CONSTANT, full basis) FINISHED at **0.9 GB peak** (not 16 GB) — the ~16 GB is only the full run's tower-depth control variants (`operatorActivated/Truncated/Extended` + kernel variants), which the residual doesn't need. ⇒ the whole cross-engine integration (single-case AND likely full-4-case primaries) is doable **on this 30 GB box**, gating only the tower controls. Coarse cross-engine consistency established (both U-rows order-3; both carry live-`W_bg` `eta_bg*w1` coupling — neither froze the constraint). ⚠ RULE 4: a prose shortcut (a clean shared `(η w1 −1)` factor) was REFUTED by the all-jets check (72/75 PY jet coeffs not cleanly divisible). A real COEFFICIENT-level residual needs the full S11c-b symbol bridge (~35 §3a basis coefficients + jet + DOF transliteration) + the `extract_slab` update to the folded PY structure — a reviewed physics-bearing cross-engine JOIN (⛔ never blanket-collapse), doable HERE, deferred to a dedicated integration pass. Evidence + record `~/.s11_build/S11c_b_constraint_fold/{cross_engine_single_case_attempt.md, wl_urow_labheld_rho4.txt, py_urow_labheld_rho4.txt}`.
- **NEXT:** **#90** — PY §3c coupling content fix (the under-extracted reversible tilted-face `A_T` + irreversible response `A_T·Λ(ω)`, per the #84 SETTLED verdict) + a §0 clarity pin (Λ = supplied flat-face closure, NOT the S11c-c bulk kernel), on the corrected basis + folded operator. Then the cross-engine integration pass HERE (transliteration + `extract_slab` update + `row_residual`, reviewed) + #88 KINETIC+θ re-adjudication + 2 owed control-hardenings → step record + S11c card + close. ⚠ The ≥64 GB box is NO LONGER on the S11c-b critical path (only the belt-and-suspenders full in-band control run remains big-box work).

## ⭐⭐⭐ S11c-b STRONG-ROW JET-DEPTH RECONCILED → SPEC PINNED (B): the slab U-row is CONSTRAINT-REDUCED; **PY is the engine that must change** (WL correct). NEXT = the PY constraint-fold BUILDER (own decision list + build legs) (2026-09-02)
The #89b PY-check flag ("WL emits order-3 in strong U-rows, PY caps at depth-2 — PY under-emits?") is RESOLVED, and NOT as first thought. Reconciliation record `directives/_measurements/S11c_b_strong_row_jet_depth_reconciliation.md`; three consistent computations (orchestrator PY probe + Codex + Grok consult, evidence under `directives/_measurements/s11c_b_jet_depth_consult_{codex,grok}/`) showed it is NOT a jet-depth freeze and NOT a physics error — a REPRESENTATION mismatch: the §3a energy has only FIRST background jets, so the RAW held-fixed `δU/δu` is order-2-complete on BOTH engines; WL's order-3 is ENTIRELY the θ-constraint reaction `∇μ_θ` (its `constrainedRowsWithLiveEnergyEL` eliminates virtual θ via the MATERIAL VIRTUAL CONSTRAINT `δ_vΣ_mat=0` [linearized] — distinct from the separate sourced mass-EVOLUTION equation), which PY keeps SEPARATE (raw EL + a separate `μ_θ` operand).
- ⛔ **SPEC-PINNED (B) — `directives/S11c_b_jet_depth_spec_pin_decision.md` (FOLDED VERDICT), 2 decision legs Codex+Grok BOTH B, raw verdict transcripts `~/.s11_build/S11c_b_jet_depth/{codex,grok}_leg_verdict.log` (`.log` gitignored).** My proposed (A)/PY-correct lean was OVERTURNED by the binding S11b rule (orchestrator-verified verbatim, rule 13): `S11b_SHARED_PHYSICS.md:341-342` "`δ_vθ,δ_v(δW),δ_vu` not independent … **Do NOT vary `U` with `θ` held fixed**" (constraint eq at `:337`) + `:426` "the in-plane equation must carry `−∇(δU/δθ)` … varying at fixed θ removes this contribution … selects the convention uniquely"; S11b's own engine `…coupling_law_mathematica_audit.wl:280` `constrainedUL = explicitUL + I k muTheta`. §1c's "θ may not be eliminated before this derivative" scopes the CONSTITUTIVE `μ_θ`, not the U-row; the separate `μ_θ` operand is a scalar face-affinity driver (not double-counting the vector `∇μ_θ`); θ's real EOM is MASS-EVOLUTION, not `μ_θ=0`; (A) would also break the §5c uniform-limit regression vs S11b's `INPLANE_EOM`.
- **BUILDER TARGET (next phase; its OWN decision list + build legs, rule 7/9 — no engine change yet):** PY `operator_from_density` (`scripts/S11c_b_brane_operator_sympy_audit.py:1968-2062`, currently the raw held-fixed EL) MUST fold S11b's material virtual constraint `δ_vθ+δ_ve_W+∇·δ_vu=0` into `U_MOMENTUM_ROWS` + thickness row; raise `STRONG_ROW_JET_DEPTH` 2→3 (NOW a genuine rule-17 freeze; PY already reproduces the 10 order-3 atoms once the constraint is applied); keep `μ_θ` a separate constitutive operand; θ-row = mass evolution. Add a §3b sentence (both legs drafted it). **WL unchanged** (`constrainedRowsWithLiveEnergyEL` is already this object; #89b's un-freeze was right). **#88 RE-ADJUDICATION:** #88 identified the strong U-rows with the raw held-fixed EL — that conflicts with the pin; its energy-basis-completion disturbance measurement STANDS, but the full-row/KINETIC adjudication is redone after PY carries the constraint reaction, and #88's θ result is a disturbance of `MU_THETA_OPERATOR`, not an independent θ equation. Then the deferred `row_residual` compares the two constraint-reduced (order-3) U-rows.

## ⭐⭐⭐ S11c-b #89b WL OPERATOR UN-FREEZE + REPAIR DONE — engine committed `a1be8d8f` + PUSHED (origin + GIN); full `.out` regen DEFERRED to a ≥64 GB box; PY sibling CLEAR (no freeze); NEXT = integration (reconcile strong-row jet depth) → #90 (2026-09-02)
The §3b variable-coefficient slab operator is un-frozen: coefficients stay LIVE through every order-raising step (EL + `activateSpatialDivergences` + `extractCoupling` + face EL + constraint), `Inactive[Div]` split preserved, jet-retaining reduction LAST (rule 17). Tractability (naively HUNG 2h+ → ~6 min): innermost-first Div activation (no `If`-on-Association) + per-summand `Series` linearity, `PossibleZeroQ`-verified. ⚠ **MEMORY WALL — INTRINSIC:** the correct un-frozen operator must hold its full jet tower un-reduced (~16 GB/case) until the final reduction, so BOTH the 2 heavy equivalence controls (behind `S11CB_SKIP_HEAVY_CONTROLS`) AND the whole 4-case in-band `.out` regen are DEFERRED to a **≥64 GB box** (`research/pde_ledger_v3/DEFERRED_HEAVY_RUNS.md`); the committed `.out` is UNCHANGED — **the deliverable is the ENGINE.**
- ⛔ **A RE-FREEZE BLOCKER was caught by the 2 build legs and repaired.** The EMITTED operator (what the cross-engine comparison reads) reduced `operatorLive` while its outer `Inactive[Div]` was still HELD → froze widthBase inside the Div → dropped the mixed-2nd/3rd U-row Hessian jets (the correct `activatedOperator`, computed at `:1377`, was DISCARDED). Plus 2 broken controls (§5.E dim walker non-functional [`Times`/`Plus` Flat/OneIdentity]; MATERIAL_ADVECTED independence → `base−base=0`). ⭐ The 2 legs DISAGREED — Grok caught it, the fresh Agent CLEARED it (its atom-PRESENCE test saw the jet inside the held Div but never the ACTIVATED coefficient); the orchestrator resolved it by its own computation on the real emit (rule 13, activate & reduce do NOT commute on the U row) → Grok right.
- **REPAIR** (`directives/S11c_b_89b_wl_operator_repair.md`): 2 decision legs (Grok 7 + Codex 4 = **11 directive-level gaps** folded — the fix ripples to the tower-depth operands, the frozen Hessian-witness across ALL slots, the uniform-limit SLAB **and** TRANSVERSE_DISPERSION, plus a rule-5 acceptance leak that could pass the frozen object); Codex repair (+28 lines, scope-verified = only the 3 fixes + ripples); **2 re-review legs (Grok + fresh Agent, rebuilt from `operatorLive` not the emit) both VERDICT CLEAR** (emit == activate-then-reduce reference, one-sided swap MOVES the residual [order-3 jets restored], all comparisons like-with-like, dim walker moves SELECTIVELY under a primitive-atom mutation, independence carries `VALIDATED_ON_REPRESENTATIVE_BRANCH` with LAB_HELD live, no regression). Record `directives/_measurements/S11c_b_89b_wl_operator_build_review.md`.
- **NEXT:** (1) ✅ **PY SIBLING FREEZE-CHECK DONE** (rule 16; record `directives/_measurements/S11c_b_89b_py_sibling_freeze_check.md`): PY is CORRECT (activate-then-reduce — emitted strong rows retain the Hessian, verified live-vs-frozen; un-frozen genuinely by #89's `total_derivative` tower), so the old agreement was NOT on a frozen PY object and **no PY freeze repair is needed.** ⚠ BUT a NEW reconciliation FLAG (spec question): PY caps strong rows at `STRONG_ROW_JET_DEPTH=2` (Hessian, no 3rd-order bg jets — 3rd lives in the coupling cascade) while WL #89b now emits **3rd-order** bg jets in the strong U-rows (the re-review's "order-3 restored") → `row_residual` WILL show a strong-row disagreement that is a **jet-DEPTH-CONVENTION question** (⛔ NOT pre-judged — it could be a genuine depth error in one engine OR a representational cap; adjudicate which strong-row jet depth the retained grade requires, don't assume WL over- or PY under-emits) BEFORE reading the residual as a physics finding. (2) **Integration** — rebuild both `.out` (⚠ WL side needs the ≥64 GB box) + re-run `row_residual` (~15 min) + **re-adjudicate KINETIC+θ** (#88 blast radius) + fold 2 owed #88 control-hardenings. (3) **#90** PY §3c content fix (face+response) + §0 pin on the corrected basis. (4) honest step record `steps/S11c_b_*.md` + S11c roll-up card + close. ⚠ Committed with EXPLICIT paths (16 files); the other session's repo `memory/` staged work was left untouched.

## ⭐⭐⭐ S11c-b #89a WL §3a BASIS REPAIR DONE — WL emits basis 40 (span verified = independent derivation), CLEARED by 2 build legs, COMMITTED `d4adbd99`; #89b WL OPERATOR unfreeze NEXT (2026-09-01)
The WL side splits: **#89a (DONE) = the §3a energy BASIS; #89b (NEXT) = the operator unfreeze.** #89a completed the WL engine's hand-coded 8/source spurion family to the full **O(3)-Kronecker** field-bilinear family (`directives/S11c_b_89_wl_3a_repair_directive.md`), imposing the exact thickness map on the new invariants; the count is computed by `MatrixRank` (`26` = the restrict-to-8 regression, the only public target). ✅ **CLEARED by 2 build legs** (fresh Claude agent + Grok, both VERDICT 1-finding/0-blockers): two INDEPENDENT SymPy derivations + the engine all agree the completed basis is **40 = 10 uniform + 15 + 15**, and — the decisive check — the engine's 15/source **span the SAME space** as the independent derivation (`rank(union)=15=each`; byte-identical scalars), ⛔ not a matching-integer coincidence. Form ablation load-bearing (restrict→26, +ε Levi-Civita→16 [parity exclusion real], drop a shear pairing→38/14); count computed-not-typed; thickness map imposed. Record `directives/_measurements/S11c_b_89_wl_basis_buildleg_clearance.md`; evidence `~/.s11_build/S11c_b_89a_wl_buildleg_{claude,grok}/`; decision-review record `_measurements/S11c_b_89_wl_decision_review.md` (2 decision legs, computation-backed, caught 8 defects incl. the operator-freeze miss + a tautological rule-17 control — folded once, rule 7). ⇒ **both engines (PY #89, WL #89a) now emit 40 on the CORRECT completed basis.**
- ⚠ **OPERATOR FREEZE — DIAGNOSED, NOT repaired (scope split; = #89b).** WL's slab operator `evaluatedModel` applies `applyProfile` (frozen 2nd/3rd jets) BEFORE the EL variation → frozen-EL rank 26 vs live-EL rank 40 (Δ14; the Hessian is the non-absorbable operator content, the WL analog of the #88 blast radius). #89a emits this as an HONEST diagnostic (`OPERATOR_BACKGROUND_FREEZE_DIAGNOSTIC`, not asserted equal, no repair claim) and DEFERS the fix to **#89b** (EL-before-freeze then live-path reduction; a different mechanism/surface — bundling was a rule-15 risk). PY #89 already unfroze its operator, so the two engines reach operator parity only after #89b.
- ⚠ **One converged control-quality nit (folded to #89b):** the §5.E new-invariant dimension residual (`…_audit.wl:2269–2290`) is a copy-check (both sides = `dimensionGradient + Σ factor dims` over the same `dofFactorSpecifications`) → structurally 0, cannot catch a wrong dimension ([[feedback_a_check_cannot_audit_its_own_input]]). Benign: the derived dimensions ARE correct (validated OUT-OF-BAND by both legs), the old hand-typed list IS removed. #89b replaces it with an independent route (expression-atom dimensional analysis) or deletes the tautological residual.
- **NEXT:** #89b WL operator unfreeze — **spec-confirm §3b WL-side FIRST** (a both-engine-shaped fix is a SPEC question first) → orchestrator-written directive → 2 decision legs [Codex+Grok, rule 7] → Codex build [WL, `--sandbox danger-full-access`, xhigh] → 2 build legs [fresh Claude + Grok; ⛔ SERIALIZE Mathematica — 2-seat licence]; the **live-operator rank is WITHHELD** from the builder (the basis 40 is ALREADY public from #89a), the FROZEN operator is the public regression. → **integration** (rebuild both `.out`, re-run `scripts/S11c_b_row_residual.py` ~15min, **re-adjudicate KINETIC+θ** [#88 blast radius] + fold 2 owed #88 control-hardenings) → **#90** PY §3c content fix (face+response) + §0 pin → honest step record + S11c card + close.
- ⚠ #89a is COMMITTED + **PUSHED** (HEAD `6f5c8c38` atop basis `d4adbd99`; `origin/ledger-v3-rebuild` + GIN both at `6f5c8c38`; the 156 MB `.out` has a `gin` annex copy). The other session's `memory/` work was untouched.

## ⭐⭐⭐ S11c-b #89 PY §3a REPAIR CLEARED by 2 build legs (both CLEAR, cross-engine); engine COMMITTED, emits basis 40 — [WL basis now DONE: see #89a clause above] (2026-09-01)
The #89 PY engine repair (un-freeze the background jet tower across the 4 frozen consumers: basis quotient, strong rows, MATERIAL pullback, coupling cascade) is BUILT (Codex, directive `directives/S11c_b_89_sympy_3a_repair_directive.md` committed `660a3055`, 2 decision legs folded) and EMITS the corrected basis **40 = 10+15+15** (confirmed against the #86 reference — established 4 ways — reduces to frozen 26). ✅ **CLEARED by 2 independent build legs** (fresh Claude agent + Grok, both VERDICT CLEAR with agreeing from-scratch derivations: 15/src live [nullity 0] → 40, 8/src → 26, the mandatory 40→26 form ablation bites, and lever-C value-exactness [0 value diffs; 17 srepr-form-only, comparator value-based]; record `directives/_measurements/S11c_b_89_sympy_buildleg_clearance.md`, evidence `~/.s11_build/S11c_b_89_buildleg_{claude,grok2}/` — ⚠ Grok/leg-B saved the full A1–C8 set; the Claude/leg-A dir saved A1–B7 and reported C8 soundness in its verdict). The engine is **COMMITTED** (pre-migration checkpoint `fce14c1a`, rewritten to `f655ea65` by the DataLad migration; `.out` now git-annex-backed, basis 40). ⚠ **Provenance caveat:** the committed `.out` was made with `S11CB_PRIMARIES_ONLY` so the HESSIAN_FREEZE / PROJECTION_EQUIVALENCE / FORM controls are SKIPPED in-band — the 40→26 reduction and lever-C identity are proven OUT-OF-BAND (both legs' ablations/numerics); the in-band control run is deferred pending the kernel-build optimization (full run OOMs).
- ⚠ **PERF SAGA (resolved; build-cost DEFERRED per user):** the full run HUNG then OOM'd (~3.5h still inside the FIRST control) — root cause = the MATERIAL coupling-kernel BUILDS (~525s each, `restrict_expression` hotspot) dominate the controls + PROJECTION_EQUIVALENCE (NOT the projection). Fixed the PROJECTION (`first_shape_series` `sp.series`→first-order Taylor "lever C", 8–25x, **PHYSICS-identical** but emits **REDUCED rational form** on ~30 scalars — build legs VETTED it (both CLEAR): value-exact, reduced-form only, and the PY-vs-WL comparator is value-based; Integral/Derivative scalars fall back to the exact reference). Added `S11CB_PRIMARIES_ONLY=1` (skips the build-heavy tasks; primaries + `exports.py` still emit) for the core output. **The kernel-build optimization is DEFERRED** (user: "as clean as we need for these physics; revisit later").
- **NEXT:** ✅ PY build legs DONE (both CLEAR) + review record committed. → **#89 WL side** (complete WL's enumeration 8→15, emit 40 [⛔ WITHHELD from the WL builder per rule 5; public 26 is the frozen-regression check]; spec-confirm §1d/§3a applies to WL → WL repair directive → 2 decision legs [Codex+Grok] → Codex build → 2 build legs [fresh Claude + Grok; ⛔ SERIALIZE the Mathematica/WL legs — 2-seat licence]) → rebuild both `.out` + re-run `row_residual` + **re-adjudicate KINETIC+θ** (the #88 blast radius) + fold 2 owed #88 control-hardenings → **#90** PY §3c content fix (face+response) + §0 pin → honest step record + S11c card + close.

## ⭐⭐⭐ S11c-b FOUNDATION FIXED THROUGH #88: #84 §3c VERDICT, #86 basis=40, #88 BLAST RADIUS (disturbs KINETIC+θ non-absorbably), #87 WL 8⊊15 — all SETTLED; #89 PY CLEARED by 2 legs (see clause above) (2026-08-31)
Worked #84 (the "first physics number" headline). Two results; **no engine had been changed as of this #84 clause** (⚠ historical — #89 has SINCE repaired the PY engine; see the top clause) — all findings, recorded in `directives/_measurements/S11c_b_coupling_84_{diagnosis,consult,basis_verification}.md`.

- **§3c CONTENT verdict — SETTLED: INCLUDE / INCLUDE ⇒ WL spec-correct, PY under-extracts.** The two engines' coupling kernels are structurally near-disjoint at the retained grades: PY = bulk stored-energy only (`gamma·profile-jet`); WL = that bulk PLUS reversible tilted-face geometry (`A_T`, ~1734 terms) PLUS irreversible frequency-dependent response coupling (`A_T·Λ(ω)`, ~660 terms). §3c mandates INCLUDING both: the "weak restriction under the stored-energy/kinetic pairing" is the EXTRACTION INSTRUMENT, not a content filter — §3c forbids the "parallel direct-variation route" (= PY's bulk-only recipe) and "filtering to a single channel", requires `TERM_ORIGINS` to classify face/flux, and its adjointness residual is expressly NOT ∂²U (rules out an energy-only object). Q2 (dissipation): §1c's "not by putting an irreversible response kernel in an ordinary action" is a construction-ROUTE rule not a deletion; §0's "memory kernel" ban targets the S11c-c curved-bulk solve (`δp=Z·v_bulk`), NOT the supplied flat-face `Λ` (verified via the S11c-a **T-i** seam, which states T-i is "not B0c's bulk-response assembly"). Λ stays SYMBOLIC (no bulk elimination). **Adversarially confirmed** by Codex+Grok across 2 rounds (each tried to break it; §3c prohibitions verified verbatim by me). User endorsed the intent. ⇒ repair PY to extract the full operator block + a one-sentence **§0 clarity pin** (Λ = supplied closure, not the S11c-c kernel).
- ⛔⛔ **FOUNDATIONAL DEFECT (bulk core = root-cause of deferred #85): PY's §3a basis independence test FREEZES the background spurion (CAS-VERIFIED); WL is undercomplete by a DIFFERENT mechanism (hand-codes 8 invariants + literal linear independence, NO quotient — NOT a frozen Hessian; see the #86 tail below).** For PY (`_measurements/..._basis_verification.md`): `basis_euler_signatures` builds its Euler–Lagrange total-derivative map ONLY from the DOF fields — the background jet `∂W_bg` is never differentiated (no `∂²W_bg` Hessian) ⇒ it uses the **uniform** divergence quotient that **§1d says does NOT lift to variable coefficients**. My reconstruction reproduces PY's exact 8-selection {01,04,05,06,07,09,10,13} (frozen rank 8), confirming the mechanism. ⇒ PY's §3a basis is **undercomplete**, and the ~118-term bulk-core coupling residual is **genuine §1d physics, NOT a coefficient-relabeling artifact** (PY kept reps 07/10 vs WL's differently hand-coded set — the reviewed comparator `research/pde_ledger_v3/scripts/S11c_b_handcoded_comparison.py:205` already localized this and refused to fold it). ✅ **#86 SETTLED (`_measurements/S11c_b_86_reference_result.md`): the corrected §3a basis = 40** (= 10 uniform + 15 `∂W_bg`-spurion + 15 `∂μ_R,bg`-spurion), reducing to the engine's frozen **26**. Per-source **15 is EXACT** (nullity 0; the earlier "over-counts / between 8 and 15" was WRONG — the would-be null-Lagrangians `g·h` are parity-ODD and excluded), and there is **NO "second facet"** (combined vs family-separated delta 0). The §1d fix = **+14** (spurion 8→15/source), carrier-independent. Established 4 ways (engine anchor 26 + Codex decision-leg + a Claude build-leg's from-scratch 1-carrier reconstruction + own crux CAS). ⚠ The rev-2 reference SCRIPT was defective (double-counted the uniform per-source → 50; removed); we accepted the triangulated 40 without rebuilding. **WL is ALSO undercomplete (8/source) but by a DIFFERENT mechanism** — code-verified it HAND-CODES 8 invariants (`newInvariantExpressions`) + takes linear independence (NO divergence quotient ⇒ NOT a frozen Hessian); 26=26 was a coincidence of two undercomplete mechanisms. #87's remaining CAS check (WL's 8 ⊂ the correct 15) is a quick confirm; #89 fixes DIFFER per engine.
- **PLAN — user chose FIX THE FOUNDATION FIRST** (tasks #86–#90, sequenced): **#86** ✅ **DONE** (corrected basis = 40, settled 4 ways) → **#88** ✅ **DONE** (blast radius SETTLED, commit `05cb1ea5`, record `_measurements/S11c_b_88_blast_radius_result.md`): correcting 26→40 disturbs **every strong stored-energy EL row** at retained grade **non-absorbably** — `RANK_GAIN 8/8/8/6/4 > 0` over the const-coeff span (u-momentum×3, θ, thickness), Hessian-zeroing ⇒ 0, so the non-absorbable content **IS the background curvature** (Hessian). Quadruple-confirmed (orchestrator anchor + Codex instrument + Grok + fresh-Claude legs, all matching; termcounts 98/98/98/108/186). ⇒ the **KINETIC (u-momentum + thickness) and θ-strong verdicts are INVALIDATED** — #89 must re-adjudicate across families, **not only coupling**. Scope: PY-side witness, LAB_HELD stored-energy rows only; mass/ADVECTIVE (constraint `ε·u_t·∇ρ_br`, not an energy-EL) likely spared; WL side + MATERIAL_ADVECTED + admissibility-ε⁰ operator are #89 (a zero here ≠ clearance). OWED (fold at #89 rebuild): harden CONTROL_ENGINE (source engine row independently of the shared density) + CONTROL_JACOBIAN (assert jet-identity). → **#87** ✅ **DONE** (commit `bab2b828`, record `_measurements/S11c_b_87_wl_subspace_result.md`): WL's 8 hand-coded invariants span a STRICT 8-dim subspace of the correct 15 (rank 8/15/15; all 8 in-span; same under the quotient) ⇒ undercomplete by exactly 7, validating the WL fix = complete to 15 → **#89 (PY ✅ CLEARED — see top clause; WL side NEXT)** both-engine §3a repair — fixes DIFFER per engine (PY: retain the Hessian in the quotient — DONE + CLEARED by 2 build legs, commit `9f40c18e`; WL: complete the enumeration to 15 — NOT started); spec-confirm intent → directive → 2 decision legs → build → 2 build legs each → rebuild operator, re-run row-residual + **re-adjudicate KINETIC + θ** (#88); repaired engines must emit 40) → **#90** PY §3c content fix (face+response) + §0 pin, on the CORRECTED basis. ⛔ COMPUTE don't assert (rule 4); ⛔ a (both-)engine fix is a SPEC question FIRST; ⛔ verify basis SPANS not counts.

## ⭐⭐⭐ S11c-b STEP 3+4 COMPLETE (2026-08-30) — ONE WL gap REPAIRED, the kinetic "gap" REVERSED
**The admissibility-θ repair is done and re-adjudicated to AGREEMENT; the proposed kinetic repair was WITHDRAWN (WL is correct).** Full rigor: directive → 2 decision legs → v2 → 2 legs → build → 2 build legs → commit → re-run → re-adjudicate.
- ⛔ **KINETIC REVERSED — WL CORRECT, NOT a gap.** Two independent decision-list CAS derivations (Codex+Grok) + §1a `e_W≡δW/W_0`: the e_W-row thickness inertia IS `μ_W W_0²` (∂/∂W_bg=0); the sibling's `W_bg²` conflates `e_W` with `e_W,bg=(W_0/W_bg)e_W`. Repairing WL would CORRUPT it. §8-item-2 ("genuine WL gap") was WRONG — a one-engine fix is a SPEC question FIRST. ⇒ ONE repair, not two.
- ⭐ **ADMISSIBILITY θ REPAIRED** (`b8443a48` decision-list, `b875cdde` build). Genuine WL under-retention confirmed, mechanism re-diagnosed: the full-field lift of the MIXED `∇θ·∇e_W` invariant (completing round-1 W1, which lifted only the pure-thickness invariant), NOT a ρ_4D density weighting. Codex 3-line change to `constructFullFieldBackgroundEnergy` (`gradLocalEw`→`gradFullEw=anchoredWidth^(-1) gradient[fullWidth]`); both build legs SOUND (real Mathematica FORM ablation zeros θ under revert; `/W_bg` normalization forced by consistency with the committed pure-thickness invariant [[7]]; N15 untouched; scope byte-identical; blindness intact).
- **STEP 4 re-run** (`scripts/S11c_b_row_residual.py` vs the repaired transcript; record `_measurements/S11c_b_step4_readjudication.md`; run `~/.s11_build/S11c_b_row_residual_rerun.out` exit 0, 26.76MB): **all 4 admissibility-θ `ROW_RESIDUAL → Integer(0)`** (WL now = PY, reconciled `kappa_theta_W·σ_W·(η w1−1)·∇²w1/(L_W W_0)`). Consistency: EXACTLY 4 residuals changed (the θ cases nonzero→0); every other family BYTE-IDENTICAL (transcript scope: only the 6 admissibility-family tags moved). No propagation.
- **PER-FAMILY STANDING:** (1) ADMISSIBILITY — REPAIRED (θ→0, agrees). (2) KINETIC — WL CORRECT; residual = representational `e_W`/`e_W,bg` normalization artifact, not a gap (⚠ owed: representational-vs-sibling-bug at family re-adjudication). (3) ADVECTIVE — representational (owed: confirm PY continuity constraint). (4) COUPLING — GENUINE in-scope disagreement, UNCHANGED (`IN_SCOPE_WEAK_REMAINDER` nonzero ×20; §3c which-engine verdict + ADJOINTNESS_RESIDUAL OWED = task #84, the "first physics number" headline). DEFERRED: energy-basis §1d quotient (task #85). **NEXT:** coupling §3c spec-adjudication (#84) + energy-basis quotient (#85); then honest step record + S11c card + full close. ⛔ COMPUTE don't assert; ⛔⛔ never blanket-collapse; ⛔ a one-engine fix is a SPEC question FIRST (this session REVERSED a "settled" gap on exactly that).

## ⭐⭐ [SUPERSEDED by the STEP 3+4 clause above — the KINETIC verdict below was REVERSED (WL correct); ADMISSIBILITY was REPAIRED (θ→0)] S11c-b STEP 2 ADJUDICATION COMPLETE (2026-08-30) — the four families are COMPUTED, not asserted
**Row-residual instrument BUILT + FIXED + COMMITTED `ef26084c`** (`scripts/S11c_b_row_residual.py`; run
`~/.s11_build/S11c_b_row_residual_fixrun.out` exit 0, 26.76MB; reporter `scripts/S11c_b_row_residual_report.py`).
**FOUR leg rounds, all folded** (directive→2 decision legs [Codex+Grok, rejected a uniform bucket-sum recipe,
9 findings; the requested truncation = retain iff c≤1∧a≤1∧b≤1 (coeffs Taylor-linearized in η), independently derived by both + me, NOT a spec ambiguity]; build→2 legs [fresh Claude SOUND + Grok found the coupling bridge-convention + a nondeterministic witness;
a mass-residual leg DISAGREEMENT resolved rule-13 = false alarm]; fix-directive→2 decision legs [rejected two
of my drafts — the Euler operator does NOT commute with the profile bridge — and converged on the construction];
fix-build→2 legs [both CLEAN, ablations bite]). The instrument reuses the layer's `extract_slab` semantic rows,
splits strong-exact vs weak-modulo-total-in-plane-divergence (coupling), computes the coupling in-scope verdict
via an instrument-side homotopy remainder (correct pre-bridge Euler). **PER-FAMILY VERDICT (record
`directives/_measurements/S11c_b_step2_adjudication.md` §8):** (1) **ADMISSIBILITY = genuine WL gap** — the ONLY
nonzero admissibility residual is the θ body force (WL≡0 vs PY's §3a ∇²w1); U/e_W/tractions all AGREE ⇒ REPAIR
WL. (2) **KINETIC = genuine WL gap** — WL froze thickness inertia W_0² vs PY W_bg² (§2a/§3b) ⇒ REPAIR WL. (3)
**ADVECTIVE = REPRESENTATIONAL** — mass-row residual = the continuity accumulation; PY imposes continuity as a
constraint, WL as an evolution row (same physics) ⇒ no repair (confirm PY constraint owed). (4) **COUPLING =
GENUINE in-scope disagreement** — IN_SCOPE_WEAK_REMAINDER nonzero ×20 (not a total divergence in-scope; genuine
bulk), consistent with v2's non-IBP-bulk certification; ⛔ a one-engine fix is a SPEC question FIRST ⇒ the §3c
which-engine verdict + the nonzero ADJOINTNESS_RESIDUAL are OWED (deep). **DEFERRED:** U_MOMENTUM/THICKNESS
complete strong rows differ genuinely (numeric: WL 3rd-derivs, PY not ∝ WL) = the §1d energy-basis
representatives (07/10 + gamma·DivGrad) that do NOT lift to variable coefficients — the owed energy-basis
quotient reconciliation (separate reps from any face-leftover first). **NEXT = STEP 3:** full-rigor WL repair of
the two confirmed gaps (admissibility θ body force ∇²w1; thickness kinetic W_bg²) → re-run comparator/adjudication;
then the coupling §3c spec-adjudication and the energy-basis quotient are the remaining sub-projects.


## ⭐⭐⭐ CURRENT FRONT — S11c-b (variable-coefficient brane operator + off-diagonal coupling kernel; S11c's FIRST physics number). **ADJUDICATION LAYER COMPLETE (3 phases reviewed-sound) + STEP 1 DONE; STEP 2 IN PROGRESS — only ADMISSIBILITY-θ is a SETTLED verdict, advective/kinetic/coupling are LIKELY/PENDING (2026-08-29). THREE commits: `bb0bfc02` comparator time-order fix (canon_jet_name recorded time differentiation as a Boolean `has_time` → WL `D[,{time,2}]`→`u_t` collapsed while PY `u_tt` stayed = asymmetric, flattening transverse-trial ∂²ₜ in COUPLING; count time tokens, emit `_tt`; both build legs SOUND, S11c-a NOT regressed [28/28 diffs canonically zero]; ⭐ discovered BY this adjudication's round-3 decision-list legs); `2e5c6755` adjudication layer v1 `scripts/S11c_b_adjudicated_comparison.py` (Bridge A `bRho↦B_rho_3/W_0` + sort-routing + total-bijection containers + exact-multiset accounting + jet-conservation; both build legs SOUND, 38 MATCH all genuine NO false MATCH via triple ablation, ⭐ Bridge A creates ZERO match — the 38 agreements are RENAME-LEVEL); `afc276b7` adjudication layer v2 (Bridge D = the engine's own committed `PROFILE_GRADE_SUBS` background expansion + IBP/total-in-plane-divergence classifier for weak-pairing DENSITIES ONLY + strong-operator EXACT + atom-gated protection; the v2 directive took 2 decision-list leg rounds catching a naive-chain-rule Bridge D [engines use `sigma_W`/`L_W` not `W_0·eta_bg`], an over-broad divergence classifier [valid only for weak densities not strong operators], and Bridge-D-∂ non-commutation; both build legs SOUND, one reimplementing the Euler/divergence test FROM SCRATCH). ⭐ MEASURED CROSS-ENGINE RESULT (rule 6 — a disagreement IS the measurement): route counts `MATCH=38, PROTECTED_UNREDUCED=32, FLAG=12, RESIDUAL_BULK=8, STRUCTURE_INCOMPLETE=57, COVERAGE=84` (231, multiset-equal, JET 231/0). The adjudication routed 38 rename-level MATCH (ADMISSIBILITY operator(16)/support(20)/COUNT(2)); 32 REPRESENTATIONAL [protected ENERGY_BASIS quotient reps 07/10 + gamma-DivGrad — non-unique by design]; and 20 INDEPENDENTLY-NONZERO operator/coupling differences (not rename/bridge-level — ⚠ step 2 refines which are genuine-physics vs representational): ADMISSIBILITY-THETA `∇²w1` (PY carries, WL 0), SLAB_TERM_ORIGINS ADVECTIVE, KINETIC `W_bg²−W_0²`, 8 COUPLING_KERNEL bulk certified non-IBP by the v2 legs. `REPRESENTATIONAL_DIVERGENCE=0` is REAL (both legs; a from-scratch Euler operator with exact rational zero-test certified all 8 coupling residuals as genuine bulk; the classifier recovers a verified `V` for known divergences up to 5218 ops). ⚠ the original "systematic higher-background-order WL truncation" hypothesis was ⛔ FALSIFIED by STEP 1. **⭐ STEP 1 DONE + COMMITTED `e5bf4122` (2026-08-29): the `(eta_bg,sigma_W)` MULTIGRADE INSTRUMENT `scripts/S11c_b_background_multigrade.py` (both build legs SOUND — fresh Claude agent + Grok each re-derived the coefficients by a DIFFERENT route with ZERO mismatch, form-ablated every guard live; one tautological guard RECONSTRUCTION folded out; extractor `scripts/S11c_b_grade_fingerprint.py`; run `~/.s11_build/S11c_b_multigrade_run.out`; records `_measurements/S11c_b_multigrade_instrument_{build_directive,directive_review,build_review}.md`). ⭐ COMPUTED FINGERPRINT (`(a,b)=(eta_bg^a,sigma_W^b)`; A=PY, B=blind WL) — the 20 differences are a per-family MIX, NOT uniform WL under-retention: ADMISSIBILITY-THETA A={(0,1),(1,1)} B=∅ (WL lacks ∇²w1 entirely); KINETIC THICKNESS_ROW A={(0,0),(1,0),(2,0)} B={(0,0)} (u-momentum rows AGREE); ADVECTIVE A={(0,1)} B={(0,0),(0,1),(1,0)} ⇒ WL carries ∇·u_t that PY LACKS (OPPOSITE direction); COUPLING BIDIRECTIONAL A={(0,1)..(4,1)}+rem(5,1) B={(0,1),(0,2),(1,1),(1,2)}. ⭐ STEP 2 ADJUDICATION (WIP, recorded `_measurements/S11c_b_step2_adjudication.md`; user chose FULL RIGOR all families 2026-08-29): (1) ADMISSIBILITY-θ = GENUINE WL under-retention [SETTLED] — §3a mandates retaining ∇²W_bg (2nd spatial deriv still first shape order), WL operand is 0. Mechanism (leading hypothesis, CONFIRM at repair by running WL): WL's `truncateBackground` Series-to-first-background-order truncation drops it (vs a missing energy term) while PY keeps it at σ_W¹. CLEAR repair candidate. (2) SLAB ADVECTIVE = LIKELY representational provenance re-bucketing (PY = ∇W_bg·u_t, WL = full ∇·(W_bg u_t); PY assigns W_bg∇·u_t to ACCUMULATION/FLUX) — CONFIRM via the deferred ENERGY_BASIS quotient reconciliation. (3) SLAB KINETIC = LIKELY genuine WL under-retention (W_bg²→W_0²; direction PY-has-more). (4) COUPLING = certified non-IBP bulk (v2 legs) with a MEASURED bidirectional grade pattern; the §2a per-cell "which engine is spec-correct" verdict is PENDING. ⛔ NEXT (full rigor): (A) SLAB advective/kinetic — the ENERGY_BASIS quotient reconciliation via VARIATIONAL-EQUIVALENCE (⛔⛔ NOT `classify_total_divergence` — v2 legs proved it invalid for strong operators); (B) COUPLING §2a/§3a per-cell adjudication; (C) then repair ONLY the confirmed genuine engine gaps (⇒ admissibility ∇²w1, likely kinetic; NOT advective if representational) → re-run WL engine + comparator + adjudication → re-adjudicate survivors (rule 13) → honest step record + S11c card + close. ⛔ COMPUTE don't assert (rule 4); ⛔⛔ NEVER blanket-collapse (the four families adjudicate DIFFERENTLY); ⛔ a one-engine fix is a SPEC question FIRST (if §3d/§2a ambiguous for a cell → fix spec). OWED (carry to card): v2 N1/N2/N3; 57 STRUCTURE_INCOMPLETE + 12 control NAMESPACE cross-engine-owed; S11c-a control-family keying; admissibility §5 control-coverage.** **[SUPERSEDED — historical; the reconcile phase below LED HERE. Read as prior state, not the result.] RECONCILE PHASE (2026-08-28): engines review-clean + transcripts committed; T7 comparator BUILT + BOTH BUILD LEGS SOUND (fresh Claude + Grok FORM-ablated every fold) + committed `5f01f9fa`; my run of the COMMITTED comparator is byte-identical to both legs (224727199 bytes; families=28 all-join, 7 unpaired, 0 verdict tokens). One comparator finding (Grok: FACE_FLUX "silent drop") verified + downgraded — visible py_only via COUPLING_KERNEL_TERM_ORIGINS, not integrity. Hand-coded reconcile `S11c_b_handcoded_comparison.py` BUILT + BOTH build legs SOUND (map COMPLETE, no false MATCH, `--drop-rename` surgical/load-bearing) + committed `82ec3b5f`; run (`~/.s11_build/S11c_b_reconcile_run.out`): **2 MATCH** (ADMISSIBILITY_SUPPORT 20/20, ENERGY_BASIS_COUNT bare-int 26), **14 FLAG** (all core objects), **12 NAMESPACE_INCOMPLETE** (control families, WL `<|…|>` unparsed → OWED). DIMENSIONS: the digest shows the dimension VALUES agree (a container-shape diff), but its reconcile verdict is still `FLAG` — confirm in adjudication, do not pre-call it. Both legs verified the 14 FLAGs are NOT naming noise (map complete); representational-vs-finding is the pending adjudication, not yet computed. **⛔ NEXT + FINAL = the ADJUDICATION (decides the number; the never-blanket-collapse danger zone):** apply the further REVIEWED bridges + re-check zero — (i) bRho→`B_rho*W_0` scale, (ii) gamma-index selected-representatives, (iii) WL `f(xOne,xTwo,xThree,time)`→PY bare GATED on a computed 0-spatial-jet check per field (S11c-a INERT_APPLIED gate), (iv) integral linearity; ⛔⛔ NEVER collapse `W_bg`/`mu_R_bg` (real background-gradient jets = the variable-coeff physics; blanket applied→bare hid the S11c-a in-plane current freeze). Reduces-to-0 ⇒ engines AGREE on that object (coupling-kernel number); survives ⇒ genuine finding (rule 1/6). ENERGY_BASIS = non-unique quotient, never fold a representative. Should be a REVIEWED artifact (reconcile-v2 with the gated bridges → 2 legs → commit) — COMPUTE, don't assert (rule 4). Likely mirrors S11c-a: agreement up to reviewed representational identities.** Journey: spec drafted → 2 legs (8 defects) `7023420d` → both blind engines built (Codex) → 4 build legs found a CONFIRMED rule-7 shared-blind-spot (both engines independently produced a vacuous ε→0 admissibility from an ambiguous §3d — a comparator alone would have read it as agreement) + a 4-leg tautological Clairaut adjoint + per-engine kernel-extraction bugs → my spec round-2 fix was STILL insufficient (2 re-legs) → **rule-15 author change: Codex rewrote §3a/§3c/§3d/§5a, 2 re-legs CLEAN** `0c0e9a8a` → engine repairs: WL W1-W3 clean both legs `49d9fad1`; SymPy B2/B3/B4 clean, B1 took a round-2 (the round-1 legs DISAGREED — Grok caught a scalar over-promotion the Claude-agent leg missed; the sibling WL engine's independent-correct §3d construction was the tiebreaker → SymPy-only fix) → **SymPy B1 round-2 clean both legs `103cdfdb`**. Both engines now review-clean. **OWED (non-blocking, carry to the card):** the admissibility operand is verified correct but has no EMITTED §5 discriminating control (independence omits on structural absence, rep-invariance is a background-order structural zero, §5b/§5c don't name it) — a control-coverage refinement, physics verified; plus S11c-a's owed control-family keying. **[SUPERSEDED — the live NEXT is the gated-bridge adjudication in the RECONCILE-PHASE front clause above; comparator + reconcile are built/legged/committed/run. The following is the historical plan record.]** Transcripts committed (SymPy `55abd09b` 167MB, WL `73b4e639` 97MB); T7 comparator (Codex `bqj3ftjqc`, brief `directives/S11c_b_comparator_build_directive.md`) was verified + 2 build legs + run + reconciled per the front. Finding records `_measurements/S11c_b_*`. [Historical spec-build detail below.] **SPEC DONE + 2 legs folded + committed `7023420d` (2026-08-27).** The SHARED PHYSICS spec `directives/S11c_b_SHARED_PHYSICS.md` was authored (mirrors the S11c-a template; inherits S11c-a §§1–3 by reference — SymPy imports `S11c_a_exports`, blind Wolfram re-derives). Two decision-list legs (orchestrator-written → Codex + Grok) ran BEFORE any engine (rule 7) and were **very productive**: they converged on 2 serious shared-blind-spots and surfaced 6 more — all 8 verified against source (rule 13) and folded ONCE (not iterated to green). The 8: **D1** ζ_s=ζ_c+sδW/2 not ζ_±=±δW/2 (was forcing ζ_c=0); **D2 (serious)** confine "u only through gradients" to the uniform quotient — the non-uniform background admits undifferentiated-u spurion couplings, the very N15 channels the step emits (my gloss would have made both engines agree the coupling is absent); **D3** declare which quantities vary (W_bg,μ_R,bg,density) vs constant moduli (N12); **D4** enumerate ONE thickness coordinate + impose e_W,bg map before the rank test; **D5 (deep)** the total-divergence quotient does NOT lift to variable coefficients (`c∇·F≡−(∇c)·F` ⇒ first-jet terms are PHYSICS) — dropped the leaky μ_⊥ citation + removed the pre-registered representative fold from the comparator (could have masked the coupling); **D6 (serious)** admissibility = background-order (ε⁰) balance in the 𝒮_hold⁰ pairing, not ε→0 of the wave operator (vacuously 0) — new §3d; **D7** sectors via local curl/div structure not a global projector (N5), adjoint w.r.t. the supplied pairing; **D8** separate one-source form ablations for ∂W_bg and ∂μ_R,bg (independent profiles). Finding record `directives/_measurements/S11c_b_spec_review.md`; leg prompt `_legs/S11c_b_spec_review.md`; raw transcripts `~/.s11_build/S11c_b_spec_review_{grok,codex}.txt`. Rule-5 leak scan clean (no coupling grade/sign, basis count, or μ_⊥ identity stated). **[historical, DONE — the plan at spec-close was]** the two blind engine builds against the spec (SymPy chained on `S11c_a_exports`; Wolfram blind), each with 2 build legs, then this sub-step's own T7 comparator + step record. OWED (carry in from S11c-a, not blocking): bridge the control-family keying (REP_INVARIANCE/CONTROL_INDEPENDENCE PY-missing-DENSITY; CONTROL_FORM FACE-granularity) + a CONTROL_FORM re-characterization.

## ✅ JUST CLOSED — S11c-a T7 STEP 5: FIVE engine fixes; the "engines AGREE / all residuals representational" claim was WRONG (a BLANKET-COLLAPSE classifier HID two findings, corrected 2026-08-27 — never blanket-collapse to judge representational-vs-finding; basic pass, HAND-CODE the rest, FLAG never massage). **FACE_SHIFT non-join RESOLVED + committed `cccb4f9e` (5th fix):** WL carried the T-e density background as an ungrounded free-premise function (`rhoBulkBackground`) — a §3c violation (same class as bg-current `6fae82b8`); grounded on the §2b representative (`rho4Profile`) across all 3 shifted-trace sites → FACE_SHIFT gained the DENSITY axis, now joins PY 16/16 at residual 0. Verified: 2 directive legs + 2 build legs (form ablation bites) + cross-engine check. Leg-2's MATERIAL_ADVECTED×RHOBR concern (both engines omit a material-transport term — rule-7 shared-blind-spot candidate) **adjudicated BENIGN** by 2 from-spec legs: §3c's normal-shift-only law correctly omits the in-plane transport (`∂_wρ⁰=0`); it lives in `δh_s` + T-f/T-g/T-h; T-e injection would double-count. Recorded non-physics nit: WL `EXACT_TRACE_SOURCE` shows lab `ρ(x)` vs §3a `ρ(χ)` — changes no computed object. **FIVE engine fixes:** `c36beac4` (shifted-trace), `6fae82b8` (WL bg-current), `49b5c525` (current-freeze #3), `8c1a5ed1` (in-plane current #4), `cccb4f9e` (WL density grounding #5). **S11c-a CLOSED — step record rewritten + committed `3b552426`** (2 review legs Codex+Grok caught 6 overclaims in the rewrite, all folded+verified). Final full comparator run (`~/.s11_build/comparator_final_cccb4f9e.out`, families=39, families_with_unpaired=11) confirms **every PHYSICS family joins clean or reduces under a reviewed identity**; the 11 unpaired are all **control/bookkeeping KEYING asymmetries** (not physics). Detail: step record + `directives/_measurements/S11c_a_faceshift_nonjoin.md`. **OWED (carry into S11c-b, NOT S11c-a-blocking):** bridge the control-family keying (REP_INVARIANCE/CONTROL_INDEPENDENCE PY-missing-DENSITY; CONTROL_FORM FACE-granularity) + a CONTROL_FORM re-characterization for the abstract-symbol objects. **NEXT = S11c-b** (variable-coefficient brane operator — the off-diagonal coupling kernel, the first physics number; `S11c_decisions.md` row S11c-b), at the SAME full rigor (uniform bar every substep — [[feedback-correctness-is-king-cost-irrelevant]] #6). The comparator + reconciliation schema (axis-typed keys, hand-code-not-blanket-collapse, ground-every-background, μ_θ/δρ/CONORMAL folds) are pinned for b–e.

⛔ **[SUPERSEDED / OVERCLAIMED — see the corrected front above] The block below reported a clean "THE TWO ENGINES AGREE, every residual representational." That conclusion depended on the blanket collapse that hid the #4 in-plane current freeze and the FACE_SHIFT non-join. Read it as the (partly wrong) prior state, not the result.**

⭐⭐⭐ **T7 CROSS-ENGINE COMPARATOR — BUILT + 2 BUILD LEGS PASS + COMMITTED `50f43123` (2026-08-27); RESULT: ENGINES AGREE.**
- **Instrument:** `research/pde_ledger_v3/scripts/S11c_a_cross_engine_comparator.py` (+ synthetic tests). Reads
  the two committed transcripts (PY `scripts/out/S11c_a_interface_geometry_sympy_audit.out` `afdc8158`, WL
  `mathematica/out/…`), and per case prints operand_A(PY)/operand_B(WL)/A−B + per-family accounting. rule 2:
  computes+prints, decides nothing. 6714 triples, exit 0.
- **THE PIVOT (rule 15):** the prose build directive FAILED THREE independent-leg rounds (rev-1 ~14, rev-2 ~8,
  rev-3 more) — the instrument is 39 bespoke per-family schemas that cannot be pre-enumerated in prose
  ("2 failed designs = wrong shape"). rev-4 = a BUILD BRIEF that locks the thrice-reviewed physics folds,
  DROPS the (B) held-context diagnostic (last #3-trap), reuses S11b's typed `residual`, and DELEGATES
  per-family extraction to the builder with mandatory accounting. Codex discovered the schemas from the
  payloads and surfaced the real axis mismatches honestly (FACE_SHIFT WL-missing-DENSITY, ADMISSIBILITY
  PY-missing-BRANCH, CONTROL_FORM 528/960) rather than manufacturing matches.
- **2 build legs (fresh Claude + Grok) CLEAR it with FORM-ablation teeth:** the current fold is a strict
  arity-preserving AST-head rename (arg-strip ablation collapses to the exact #3 false-agreement → teeth);
  the integral canon is a capture-safe BoundIntegral (bound-equality ablation collapses a real disagreement);
  object-nested controls extract; no asserts on measured payloads, no verdict emitted, MU per-branch, CONORMAL
  raw. Non-blocking: PROJECTION_TERM_ORIGINS labels DYNAMIC/STATIC positionally (symmetric → can only surface a
  false DISagreement, never false agreement).
- **RESULT (rule-13-classified; detail SCOUT §24):** 8 physics families reconcile to exactly 0; TRACTION /
  CLOSURE / VIRTUAL_WORK are notational-only (collapse μ_θ(x)→μ_θ ⇒ residual 0 — the whole TRACTION residual
  is `coeff·(mu_theta_L − mu_theta_L(x1,x2,x3,time))`); PROJECTION reduces under collapse + full integral linearity
  to the density-time term δρ_4D vs ρ_br·θ; CONORMAL = §3c verdict-A. ⭐ The OLD reference reconciler had HIDDEN
  the μ_θ difference by arg-STRIPPING (the FIELD `AppliedUndef→Symbol` smuggle) — the honest arg-preserving
  comparator SURFACES it. NO genuine physics disagreement (pending the adjudication of the representational
  identities). NEXT = adjudicate (are WL's face-eval args inert per §3a/§1b + is δρ=ρ_br·θ) → controls
  bite/hold → step record.

⭐⭐⭐ **T7 CURRENT-FREEZING FINDING (#3) — FIXED + COMMITTED `49b5c525` (2026-08-26).** The projection
integrand — which the prior roadmap filed as a "mechanical bridge, NOT a finding" — was a REAL PY defect,
caught ONLY because the user flagged building a "field-name map" as a known pitfall. PY's `projection_terms`
froze the perturbation current's normal component at its face value (w-constant `j_bulk[3]`), so
`WINDOW_NORMAL_CURRENT = -ε∫ j_w·∂_wΩ` was IDENTICALLY 0 (constant × ∫∂_wΩ = boundary = 0) — the entire
`∂_wδj_w` contribution ABSENT, violating §1b (`∇₄·δj` includes the normal divergence). §3c's zero-background/
trace language is scoped to TRACED face fields, NOT the bulk current under the `w`-integral.
- **CAS-verify (Codex+Grok, runnable SymPy+stdout):** Q1 (does the current's w-dependence matter?) = DIFFERENT
  (concrete nonzero witnesses); Q2 (window shape-derivative form) = SAME/benign ⇒ the disagreement is the
  current, not the window.
- **From-spec adjudication (Codex+Grok):** UNANIMOUS — §1b requires the full w-dependent bulk current (WL
  faithful, PY defect); an implementation defect, not a spec ambiguity. Orchestrator rule-13 self-verified.
- **Fix:** `δj_w = delta_j_bulk_4(w)` INSIDE the existing post-IBP normal term (⛔ NOT a second `∫Ω·∂_wδj_w`
  channel = double-count), both dynamic/static branches + `uniform_projection_reference`. 2 build legs (fresh
  Claude + Grok) ablation-PASS: `∂_wδj_w` enters (a w-varying probe moves the term, a constant collapses it to
  the old 0), single channel; pre-fix `WINDOW_NORMAL_CURRENT` was literally ≡ 0 (direct confirmation of #3).
- ⚠ **Directive took rev-2:** rev-1 over-specified a FACE-keyed recipe (`affine_bulk_perturbation`) for a
  FACE-LESS projection (rule 3 — name the object, not the recipe; the SAME over-specification I made on the
  comparator directive). Both caught by the 2-leg directive gate BEFORE any build. ⛔ A §5c uniform-limit
  "corroboration" was RETRACTED as a non-simplified-zero false alarm (a CAS-integral zero-test needs INTEGRAL
  LINEARITY, not expand/cancel).
- ⚠ [SUPERSEDED by the front block — the comparator IS built + run 50f43123] The ULTIMATE cross-engine projection AGREEMENT was NOT yet computed AT THAT TIME — needs the current name
  map (LEGITIMATE now: both engines carry the same w-dependent δj_w) + ITG/IBP canon = the DEFERRED comparator.

⭐⭐ **[SUPERSEDED — comparator BUILT + committed 50f43123 on 2026-08-27; see front block] COMPARATOR DEFERRED (2026-08-26 — user chose "lean measurement first, then batch fixes").** The rev-1 build
directive (uncommitted WIP `directives/S11c_a_comparator_build_directive.md` +twin +review leg) drew ~14 defects
from its own 2-leg gate: WRONG PY input (`exports.py` is a LEDGER; real PY = an uncommitted scratch transcript — now the FIXED `~/.s11_build/S11c_a_py_fixed_run2.out` (47 tags,
post-`49b5c525`), loaded by a `measure_reconcile.py` in the PREVIOUS session's /tmp),
false "S11b reassembles multi-line WL", dropped FACE_SHIFT, `_RESIDUAL_OPERAND` (engines emit `_RESIDUAL`),
**2 SMUGGLING folds** (5 `mu_theta_L/M→mu_theta` branch-collapse, 9 CONORMAL Taylor) → must go to a REGISTRY not
the map, `canon_key` CASE-LOSS (full-axis keying needed), DIMENSIONS/BACKGROUND_DENSITY_MAP mislabeled supplied.
Lean measurement instead cleared the CONTROLS (no new finding; REP_INVARIANCE/UNIFORM_LIMIT all-0 invariants,
CONTROL_INDEPENDENCE/CONTROL_FORM/HOMOGENEITY bite; 7 parser false-alarms caught by rule-13). FACE_SHIFT/origins/
bookkeeping = deferred to the comparator. Detail `~/.s11_build/S11c_a_T7_SCOUT_FINDINGS.md` §§17–22.

⭐⭐⭐ **WL BACKGROUND-CURRENT FINDING — RESOLVED + FIXED + COMMITTED `6fae82b8` (2026-08-26).** The T-f
projection's "possible background-current finding" was REAL: WL carried the rest-frame background bulk current
as FREE PREMISE symbols (`currentWBackground`/`currentXBackground{i}`) — violating §3c ("none may be introduced
as a free premise"), §1b (`j=ρ_4D v`), and its OWN `bulkVelocityZero→0`; it survived in the dynamic-window
projection (1660 hits, 8/8 cases; absent in flat STATIC). PY (post `c36beac4`) had `j⁰=ρ⁰·v⁰=0` by construction.
- 2 from-spec adjudication legs (Agent+Grok) confirmed WL diverges. ⚠ **USER PIVOT (rule-6):** symbolic-carry is
  arguably the MORE informative computation — PY hardcoded `v⁰=0` (assumed the answer). 2 physics consults
  (Codex+Grok, own CAS+stdout) then COMPUTED the survivor: after IBP in `w` it is `∫δΩ·(∇₄·j⁰)`, and static
  continuity `∇₄·j⁰=0` (`∂_tρ⁰=0`) ⇒ total `w`-derivative ⇒ boundary term ⇒ **0**. So NOT a physical
  drift-coupling — an artifact of unconstrained free symbols — but WL is internally inconsistent.
- **USER DECISION:** build `j⁰=ρ⁰·v⁰` in WL and RECORD the continuity-cancellation in the step record (not a
  silent hardcode). Fix: `bulkCurrentZero=rhoBulkZero·bulkVelocityZero`; `currentWZero=Last[…]`;
  `currentXZero[index_Integer]=…[[index]]` (the `_Integer` restriction fixed a `Part::pkspec1` symbolic-index
  bug that also ballooned the run to 14 GB). Directive rev 2 `1610b9e9` (2 legs: Codex found 3 real defects — a
  leaked downstream result, an acceptance that couldn't distinguish construction from `:=0`, weak §6 clauses;
  Grok clean). Fix `6fae82b8` (2 build legs Agent+Grok PASS — the velocity-probe FORM ablation shows the bg
  current TRACKS a nonzero probe ⇒ genuine `ρ·v`, not `:=0`). Regenerated `.out`: 40 tags, **0** bg-current
  (was 667), only 10 bg-current-consumer tags changed, 30 byte-identical.
- ⚠ OPS this round: 2 SPURIOUS mass-kills of Codex Mathematica jobs + 1 real 14 GB OOM (orphan kernel, killed
  by PID). Decoupled the heavy full engine run FROM Codex (orchestrator runs + monitors it; Codex does source
  edit + a tiny isolated probe only). The full run ≈14 GB / ≈14 min exceeds review-legs `timeout 600` ⇒ legs use
  the isolated velocity-probe + the provided `/tmp` `.out`, never the full engine.

⭐⭐⭐ **T7 STEP-5 COVERAGE (measure-first, 2026-08-26).** **11 tag-families reconcile EXACTLY to 0** (both blind
engines agree): FACE_NORMAL/MEASURE/VELOCITY, VIRTUAL_CONSTRAINT, RELATIVE_FLUX, TRACTION, CLOSURE, EVOLUTION,
KINEMATIC (OPERAND_A+B), VIRTUAL_WORK (16). Instrument `~/.s11_build/S11c_a_cov_all.py` (grows the reconcile_fixed
declared map; committed reconcile_fixed.py / reconcile_traction_check.py UNCHANGED). Map extensions: +4 params +
the **EVOLUTION face-detect** (WL preprocessor reads `{±1/2*W0}` → `XFACEX` tokens; EVOLUTION 8/8→0).
- **CONORMAL (T-a′): VERDICT A — same physics, NO finding** (2 from-spec CAS legs residual 0 + WL-source check:
  WL's `W_0²`/3rd-order = §3c flat-face Taylor of PY's background-face traces, η-correction retained;
  `conormalPerturbation` is the probe wave, not δn̂).
- ⛔ **PROJECTION "window drops a face" was a FALSE finding — rule-13 catch.** MY prompt called WL "single-arg"
  from a 70-char TRUNCATED fragment; the raw WL window is the identical 2-arg `𝒪(w−W0/2,−w−W0/2)`, both-face
  derivatives. Lesson: a leg prompt must carry the ACCURATE engine form.
- **[SUPERSEDED — now FOUR fixes; see corrected front] NET: ~~THREE~~ FOUR confirmed T7 findings, ALL FIXED** —
  shifted-trace (PY `c36beac4`) + free-premise bg-current (WL `6fae82b8`) + current-freezing/normal-jet
  (PY `49b5c525`) + **in-plane current freeze (PY `8c1a5ed1`, 2026-08-27 — the 4th, hidden by the blanket
  collapse)**. PY carried three §3c/§1b-class implementation errors, WL one; the two-engine method caught all
  four. FACE_SHIFT non-join + a full comparator re-run remain OPEN. Detail
  `~/.s11_build/S11c_a_T7_SCOUT_FINDINGS.md` §§10–25.
- **NEXT (SUPERSEDED — the projection integrand was NOT mechanical; it was finding #3, fixed `49b5c525`):**
  the DEFERRED comparator — rebuild the rev-1 directive grounded in reality (commit a PY `.out`, port the
  ephemeral loaders + full-axis extractors, folds 5 `mu_theta`/9 CONORMAL → a REGISTRY not the map), 2 legs,
  synthetic fixtures → freeze + run → the full T7 cross-engine result (which should now show the projection
  AGREEING post-fix, and sweeps FACE_SHIFT/origins/bookkeeping) → step record + family card (carry ALL THREE
  §3c fixes + the recorded bg-current continuity-cancellation).

⭐⭐⭐ **T7 SHIFTED-TRACE FINDING — ADJUDICATED + FIXED + COMMITTED `c36beac4` (2026-08-26).** The dual-engine
method caught **two real SymPy bugs; WL was correct.** Adjudicated by 2 independent from-spec CAS derivations
(fresh Agent + Grok, each residual 0) + a source check: **(D1)** PY fabricated background-normal-jet PREMISE
symbols the spec never supplies (`d_w_v_bulk_0`/`d_w_delta_p_0`/`d_w_j_0`/`d_w_rho_4D_0`) → a spurious §3c
term-2; **(D2)** PY froze the traced perturbation at the flat face → dropped the mandated ε·η term-1 correction.
- **Fixed spec-question-first:** §3c clarified in the shared spec (structural premises only — every background
  jet is `∂_w` of a supplied `𝔅⁰` member, genuinely zero for velocity/pressure/current + `w`-independent
  density; traced perturbation evaluated at `h_s⁰=sW_bg/2`), then the PY engine fix.
- **Full gate:** 2 directive-review legs (Codex+Grok) caught a leaky, under-scoped rev-1 directive (killer
  grep-dodge, corrected-monomial leak, missed FACE_SHIFT/CLOSURE/VIRTUAL_WORK + current/density jets) → folded
  once → rebuilt; then 2 ablation legs (fresh Agent+Grok) both PASS — their form ablations proved term-2 now
  DERIVES from the real background (fires when nonzero, vanishes when zero), not deleted blindly.
- **Cross-engine confirmation** (fixed PY vs committed WL, via `~/.s11_build/S11c_a_reconcile_fixed.py` =
  measure_reconcile + the mechanical `d_w_X`↔`X_dw` perturbation-jet rename): **RELATIVE_FLUX 8/8→0** and
  **TRACTION 16/16→0** (TRACTION also needs the mechanical `mu_theta_L→mu_theta` rename + `sp.cancel` for the
  λ_X complex-denominator factoring — naming/CAS-form, not physics); the 24 pure-geometry cases stay at 0.
  Measurement `directives/_measurements/S11c_a_py_shifted_trace_fix_directive.md`. Engines at THIS fix: WL
  `a7459cb8` UNCHANGED (was correct), PY corrected in `c36beac4`. ⚠ WL was LATER fixed to `6fae82b8` for the
  bg-current finding — see the CURRENT FRONT block at the top.
- ⚠ **EVOLUTION residual = a COMPARATOR harness gap, NOT a finding:** EVOLUTION's case key has no FACE axis,
  so the instrument applies one face to a two-face-SUMMED expression and drops the minus-face field evaluations;
  the raw WL expression is complete (both faces, brackets balanced). Fix = detect the face from the `{±1/2*W0}`
  evaluation-point argument (comparator hardening), not the engine.
- ⛔ **[SUPERSEDED by the CURRENT FRONT block above]** The shifted-trace-era "still to measure" plan is done:
  CONORMAL adjudicated VERDICT A; CLOSURE is ALGEBRAIC (16/16→0, NOT a window integral); the PROJECTION
  "window" was a FALSE finding (rule-13); the bg-current WAS the finding hiding here and is FIXED (`6fae82b8`).
  ⛔ CORRECTED (see top block): the projection integrand was NOT mechanical — it was finding #3
  (current-freezing, FIXED `49b5c525`); the controls are measured clean; only the DEFERRED comparator build
  remains (FACE_SHIFT/origins/bookkeeping sweep folds into it).

⭐ **T7 STEP 2 (engine patches, fixes on BOTH engines).** WL patch DONE: added the density-representative axis
to the 5 T-f projection objects (spec §1b/§3a/§4:398), dropped the spurious axis from KINEMATIC/RELATIVE_FLUX
(they use `ρ_m`, §3b:351), extended §5b/§5c control coverage to the 5 missing objects (§5b:490/§5c:497). A
repair fixed a patch-introduced blocker (the reset refactor cleared `materialShape` without redefining it →
§5a rep-invariance material route went inert for T-c/T-d/T-i). Two fresh legs (Agent + Grok) PASS; Grok's F1
(static-projection form control) resolved by computation to a working control (fires 6/24 where the RHOBR
material-advected density genuinely carries `∂W_bg`; honest invariant elsewhere) — a step-record note, not a
defect. **PY patch DONE** (full physical×virtual virtual-work grid — 16 cases, off-diagonals genuinely
computed; BACKGROUND_STATE boundary loads f_hold/t_hold; BACKGROUND_DENSITY_MAP branch drop). Two fresh legs
(Agent + Grok) PASS; both independently reached — and both, with the orchestrator, judged non-defect — the one
finding (BACKGROUND_STATE dimension tuple is declarative/emitted-only, a pre-patch convention that skips the
supplied zeros; step-record note + a comparator DIMENSIONS-bridge note). [SUPERSEDED — see the Step-3 front at top] NEXT = shallow reconciliation bridges
→ the trivial T7 cross-engine comparator → step record + family card (pin the schema forward for S11c-b…e).

⭐⭐ **S11c is a STAGED FAMILY a–e** (decision list `research/pde_ledger_v3/directives/S11c_decisions.md`
`24853f3a`, N1–N15, 2 legs): **a**=background & geometry, **b**=variable-coeff brane operator/kernel,
**c**=curved-interface bulk DtN closure, **d**=profile-conditioned spectrum, **e**=leakage/confinement/
bench-optics falsification. Building **a** first (user chose staged: checkpoint before the blind WL engine).
⭐ **S11c-a SPEC CLOSED `2926c71c`** (`directives/S11c_a_SHARED_PHYSICS.md`) after a long rule-15 arc
(draft→my fold BRED defects→**Codex re-author**→repair-1→repair-2→final; both final legs PASS + CAS-verified;
core: perturbation-`δp` ⇒ `J_s⁰=0` rest, drain off; N4 branch map; T-h own-law corruption). **SymPy build
directive `304fa46f`** (2 legs; fixed my `ρ_m` re-origination → bind `LEDGER['rho_m']`; comparator paraphrase
→ point at G8(a)+T7 whole). ⭐ μ_θ PIN: both engines use §3a's branched `μ_s^α`.
⭐⭐ **S11c-a SymPy ENGINE CLOSED — committed `9b6438fa`** (build `buo60i510` → reviewed baseline `488c2a65`
→ §5a repair `b41ueycyz`). `scripts/S11c_a_interface_geometry_sympy_audit.py` (~2096 lines) +
`scripts/S11c_a_exports.py` (2093 frozen rows). **Two review rounds, four fresh legs total.** Round 1 (Grok +
Agent) found a BLOCKING defect: §5a route-2 (material-coordinate) was a byte-identical alias of route-1 for
T-c/T-d/T-i/T-g — a control that could not fail. **Codex repair** (directive relayed both legs VERBATIM,
rule 15): reimplemented route-2 as the genuine `w′` face-flattening derivation (`build_material_face_source`,
covector mapped back via `(dx/dX)^-T`) + one-sided corruption on the direct route only + T-c′ two operands.
**Round 2 (fresh Grok + Agent): BOTH CLEAR — safe to commit.** Primary T-0…T-i correct (independent derivs
MATCH; orientation `s(n̂·ŵ)>0`; `ρ_m` bound; virtual `δ_vu`) + form-ablation clean; §5a route-2 now genuinely
independent for all five (one-sided AND two-sided corruption isolate); uniform-limit→S11b exact (19
primaries). Step-record notes (NOT defects): T-c′ residual is a definitional identity; T-g's §5a residual is
a structural determinant identity (block-triangular Jacobian ⇒ det=det(F)·thickness) — independent
constructions, weak only vs a shared density/anchor mistake ⇒ WL engine is the cross-check.
⭐⭐ **S11c-a blind WL ENGINE — COMMITTED `ddecdbc2` (repair-round-2 done + verified).**
`mathematica/S11c_a_interface_geometry_mathematica_audit.wl` (blind: imports nothing, 40 WL_S11CA_ tags) +
committed `mathematica/out/*.out`. Build `277f3fe7` → §5a/T-h repair `a15bc69c` → **repair-2 `ddecdbc2`**
(T-f divergent `∫1 dw` + T-0 σ_W grading). Repair-2 directive `bf3c3acc`→`30e96ee2` (Codex clarity leg1 +
Grok directive leg2 = 3 CAS-backed edits: held-factor T-0 recipe [the "retain W_bg" wording was a NO-OP],
stripped a rule-5 grade leak, broadened the T-f confirm). Codex build → **2 fresh WL legs (Grok + Agent,
serialized) BOTH "safe to commit, no blocking defect"** — T-0 (RHO4 σ_W first-jet `{0,0,1}`, RHOBR zero
`{{0,-∞,-∞}}`) + T-f (0 divergent ∫1, genuine per-term localization) independently re-derived; FORM ablations
move exactly the data-dependent tags; §5a independence genuine. My authoritative re-run (run4): exit 0, 40
tags, ZERO divergent ∫1, byte-identical to Codex definitive run (40/40); regression vs a15bc69c = 28 identical
/ 12 changed (exactly T-0+5×PROJECTION_*+6×CONTROL_FORM/UNIFORM_LIMIT). ⚠ LESSON: my "88 ∫1 in ORIGINS" alarm
was a STALE INTERMEDIATE SCRATCH file — verify against a FRESH run I control, not the builder's scratch stdout.
⭐⭐⭐ **T7 CROSS-ENGINE — the engines DISAGREE on CASE STRUCTURE; user chose FULL RECONCILE; fixes land on BOTH engines.**
Roadmap `directives/S11c_a_comparator_reemit_plan.md`. Arc: comparator-side key-alignment REJECTED by 2
directive legs → user chose re-emit → Codex found the emit-only assumption FALSE (genuine case-structure
divergences) → user chose **MATRIX-FIRST** then **FULL RECONCILE**.
**STEP 0 — adjudication MATRIX `3c7f9137`** (`directives/S11c_a_T7_adjudication_matrix.md` + twin + frozen
census `_measurements/s11ca_t7_census/`): a COMPUTED axis-set census over both hash-locked streams, 2 legs
(Codex+Grok) + fold. 39-tag join; classifies every divergence.
**STEP 1 — adjudication VERDICTS `3491a376`** (`directives/S11c_a_T7_adjudication_verdicts.md` + twin +
`_measurements/s11ca_t7_adjudication/`): every verdict = decidable computation + spec citation; 2 legs
reproduced every engine ASSIGNMENT (no reversal) + corrected severity/scope; folded, each re-verified (rule 13).
**VERDICTS (fixes on BOTH engines):** (B) density → **PY correct** (WL drop redundant axis on kinematic/flux
[use `ρ_m` §3b:351, rep-identical, SHALLOW]; WL add axis to rep-dependent projection SHAPE_DERIV/_DYNAMIC/
_RESIDUAL/_TERM_ORIGINS [COMPUTATION]); (C) virtual-work → **WL correct** (full grid §4 T-d:419; PY missing
pairing CASES, off-diagonals physical-DOF-REDUNDANT not new physics; PY emit full grid → cross-checks WL);
(H.1) coverage → **PY correct** (WL add form+uniform for 5 quantities; spec NOT under-specified); (BG) →
**WL correct** (PY add BOUNDARY_LOADS only, already has the zeros §2d:251); BACKGROUND_DENSITY_MAP → PY branch
axis redundant (2-per-rep, §2b).
**STEP 2 DONE — both engine patches committed (WL `a7459cb8`, PY `2d0f0055`; each build + 2 fresh legs PASS;
WL also needed a §5a `materialShape` repair).** [SUPERSEDED — see the Step-3 front at top] NEXT = shallow reconciliation bridges (leaf-rep reconstruction, decompositions,
encoding) → trivial join+residual comparator (frozen T7 = `scripts/S11b_cross_engine_comparator.py`) + 2 legs
→ step record → family card. ⭐ pin the adjudicated schema FORWARD for S11c-b…e. ⚠ PY emit feeds
`export_candidates` — any reformat mutates `S11c_a_exports.py` (decouple/preserve/review). See memory
`project_s11c_a_state`. ⚠ another session works untracked in `memory/` — commit EXPLICIT paths only. ⚠ codex/
grok jobs die SPURIOUSLY; a spurious kill can orphan a WolframKernel eating memory (kill by PID).

---

### S11b (CLOSED) — historical recap

⚠ **PATH BASE:** every S11b artifact pointer in this block (`scripts/…`, `mathematica/…`, `steps/…`,
`directives/…`, `paper/…`) is relative to **`research/pde_ledger_v3/`** — e.g. `scripts/S11b_cross_engine_comparator.py`
= `research/pde_ledger_v3/scripts/S11b_cross_engine_comparator.py`, and the card `paper/steps/S11b_interface_coupling_law.tex`
= `research/pde_ledger_v3/paper/steps/S11b_interface_coupling_law.tex`. Only `STATUS.md` and `CLAUDE.md` are repo-root.

⭐⭐⭐ **S11b IS CLOSED (`565b3fe8`, 2026-08-23).** It was rebuilt as ONE unified export-chain step ("the linear
brane–bulk interface coupling law"), subsuming the historical A (bulk face response) + B (homogeneous assembly)
execution stages; **C is now the separate later step `S11c`.** SymPy engine `864d6f41` + blind WL engine
`ec89f9df` cross-checked by the frozen T7 comparator (`17fe32c8`, re-run `fba6a34c` — physics agree; the
**X-1** energy-basis over-count 11→10 corrected `53fcd98d`) + WL repair `bd598ae7`; step record `8ddccb74`;
card `565b3fe8`. **All 13 step-run steps + export integrity done.** The live NEXT is the S11c block below;
the artifact history is retained for provenance. Old spec/directives cite the deleted `reduction/` + a rule-12
denylist.

⭐⭐ **Artifact #1 DONE — `directives/S11b_unified_decisions.md` (G1–G14), committed `ddd0ae4c`** — two
legs (Codex + Grok), 7 findings folded once, rule-2 `_measurements` twin (which caught 2 of my citation
errors). Both design snags resolved (Codex consult + my verification; the A→B difference is a
GENERALISATION, ⛔ not a rename):
- **`v₀`** — S11's (in-plane brane, pinned 0) ≠ S11b's (**normal bulk drain**): DIFFERENT quantities ⇒
  distinct names; ⛔ never reuse the bare `v_0` key (`F9` object-compares two `Symbol('v_0')` as **EQUAL**
  and silently merges them — Grok reproduced it). ⛔ NOT a premise override.
- **`Λ`** — `Λ_A⁰` (affinity `𝒜=μ_s−δp/ρ_m`) ≠ `Λ_p⁰` (raw pressure); keep distinct, **A recovered as a
  `μ_s=0` reduction slice**; three independent times `τ_A,τ_V,τ_X`; `Λ_X` a **SUPPLIED** channel.
- **Scope (G14, corrected after both legs)** — the **UNIFORM** longitudinal fate + grazing are THIS
  step's (B5); only the non-uniform variable-coefficient spectrum → C; the DC/harmonic/sideband radiation
  audit → the **nonlinear program** (⛔ not C).

⭐⭐ **Artifact #2 DONE — `directives/S11b_SHARED_PHYSICS.md`, committed `1a2395a3`** (the unified A+B spec;
⛔ replaces the superseded rev-2 whole-interface spec that sat at that path). Two legs (Codex 14 + Grok 11,
convergent on an **incomplete bulk-response provenance-move**); folded once (rule 7), each verified myself.
Folds: `Z`/regimes/dissipation **DERIVED in B0b, not stated** (§2 prohibition narrowed; `q_out`'s supplied
branch structure exempt); ⭐ **physics fix** — `μ_s=0` turns the flux pressure-driven but keeps the velocity
channel; `ρ_br⁰` = `rho_br` **settled** (frozen wall width); restored the anti-tautology clauses in the two
energy discriminators; **cut the `v₀/c_s0` convective-order leak** (derive-only, G12a); distinct drain glyph
`v_dr` (⛔ never `v_0` — the F9 false-merge hazard); added the k=0 **breathing-slice task**, the permeable
regime/parity audit, the flux-channel + basis-count tags, the exact one-parameter control cuts, a
non-passive-root power-source obligation; **locus protocol restricted to equation-system loci** (B2d region
& B5 grazing exempt); self-containment fixes (4D continuity law supplied; `μ_s=0` velocity-channel; B2c→B0c);
T7/exports summary → a pointer. rule-2 twin + review prompt committed.

⭐⭐ **Artifact #3 DONE — the two per-engine build directives, committed `9bd2f184`** (each wraps the spec
§§0–13, adding only per-engine wiring). **SymPy** (`directives/S11b_sympy_build_directive.md`): imports the
S11 LEDGER, binds `c_s0`/`μ_R`/`ρ_br⁰=rho_br` to the imported objects (⛔ never re-declared), `ρ_m`/`v_dr`
originate here, writes `S11b_exports.py` (F1 flat keys, F3 row evidence, three-valued F9 applied **whole**,
digest pin `{own source, S11_exports.py, spec}`, `MappingProxyType` freeze; the three **G8** deviations —
comparator SEPARATE to G8(a)+T7, D3 restored, `_RELATIONALS` included; F6 publish gate). **Wolfram**
(`directives/S11b_wl_build_directive.md`): BLIND (imports nothing, re-derives), blindness by **absence** ⛔
not a denylist, no VERDICT, one-kernel/two-seat run discipline. Two legs (Codex+Grok), folded **once**: both
killed the SymPy F9 **paraphrase** (it dropped F9's *"TOTAL over the imported row shapes"* clause — F9's own
measured failure) ⇒ point at F9 whole + only the S11b-specific wirings pointing can't supply: the `s11b_` F9c
prefix (⭐ `s11_` is **already taken** in the import ⇒ would collide) and the **§5→§8** locus-protocol remap;
Codex caught a **T7 mis-attribution** (join/residual/repoint are G8(a)'s, not T7's); self-caught a
`B0–B7,B9` range imprecision (→ "every §9 task except B8"). WL directive **CLEAR** both legs. rule-2 twins
script-generated; review prompt committed.

⭐⭐ **Artifact #4 (SymPy engine) DONE — committed `864d6f41`** (3 build rounds, hard stop honored).
`scripts/S11b_interface_coupling_law_sympy_audit.py` → `S11b_exports.py` (1916 rows: 1663 S11 carried
forward + S11b) + committed `.out`. Round-1 build → 2 script legs found the physics CORRECT (both **script legs'**
independent derivations match every load-bearing object incl. the G13 criterion `Λ_p⁰=−Λ_A⁰/ρ_m`; ⛔ NOT
cross-engine — the WL engine is not built yet) but two
emission defects: a TYPED load-bearing impedance (F4) and ~12 EXPORTED tautological check-rows. Round-1
repair (directive `0863645b`): impedance now COMPUTED (`BULK_ACOUSTIC`), the big §6 checks bite. Round-2
repair (directive `75560832`, user chose "clean the export" then a Codex consult corrected the disposition):
made genuine the 4 tail checks — kernel_orientation (separate retarded reference; global transpose now bites),
onsager crosscheck (independent all-flux/all-force solves + monic compare), both parity branches (live
integrand), the B8 control (per-residual). Both round-2 legs: "nothing survives the physics filter";
value-preservation clean (non-target objects byte-identical; carry-forward intact). ⛔ **HARD STOP — no
round-3.** ⚠ Lessons: I twice leaked an expected value into an acceptance criterion (rule 5); a review
agent reverted `exports.py` to HEAD in cleanup (regenerated); a spurious kill fired after the engine's work
completed (`[[feedback-background-tasks-can-die-spuriously]]`).

⭐⭐ **Artifact #5 (blind WL engine) DONE + REVIEWED — committed `ec89f9df`** (`mathematica/S11b_..._audit.wl`,
178 tags, imports nothing, no VERDICT; `.out` captured). 2 script legs: physics CORRECT, but **F-WL-1** = the
withheld G13 map `ZPERM_SLICE_MAP` emitted WRONG SIGN + WRONG LEVEL (`+Λ_A/ρ_m` dynamic; correct = static
`Λ_p⁰=−Λ_A⁰/ρ_m`, verified 4 ways), **F-WL-2/3** = a tautological/decoupled §6 energy + causality/kernel/
grazing check layer. Codex consult (`bee278a7`) found **X-1**: SymPy basis OVER-COUNTS **11 vs correct 10**
(spec §5 mandates the total-divergence quotient; `st_squared=(2/3)div²+(1/2)curl²` mod div) ⇒ X-1 REOPENS the
committed SymPy engine (count/coeff only; the physical EOM span is unchanged but the emitted EOM COORDINATES shift (mu_S absorbs into mu_R,B_div: mu_R_eff=mu_R+mu_S/2, B_div_eff=B_div+2mu_S/3 -- do NOT require EOM byte-preservation)). Disposition: `steps/S11b_wl_engine_review_disposition.md`.

⭐⭐ **Artifact #6 (T7 comparator) DONE + FROZEN + RUN — `17fe32c8`.** `scripts/S11b_cross_engine_comparator.py`
(reused S10; build `c95213b7` → fix1 `bd97a571` → fix2 `f0c796a9`, each 2 legs; precedence rule
`DISAGREE>UNCOMPARED>UNDECIDED>AGREE`). Frozen run `scripts/out/S11b_cross_engine_comparison.out`: **X-1 confirmed**
(ENERGY_BASIS_COUNT 11 vs 10, CONTENT), DIM_THICKNESS_RESPONSE off by `(-1,0,0)=W₀` (known convention), all other
DIM AGREE, DEGENERATE_LOCI a locus-form/qOut-vs-q artifact (not physics). ⚠ FINDING: the two engines physics-agree
but are NOT emission-parallel — 103 STRUCTURE (PY Tuple vs WL Association) + 47 UNPAIRED + naming (qOut/q); only
20 clean AGREE. Comparator confirms NO cross-engine physics disagreement beyond X-1 + the known F-WL-1 sign.
⚠ **That frozen run (`17fe32c8`) PREDATED both repairs; the comparator has since been RE-RUN — see Artifact #8.**

⭐⭐ **Artifact #7 (WL engine REPAIR) DONE — committed `bd598ae7`** (directive `e273398a` + fix-2 `a6b0b684`;
baseline `a5186dce`). All F-WL fixes GENUINE + physics CORRECT: **F-WL-1** now emits static `−Λ_A⁰/ρ_m` (the
withheld G13, representation-invariant, ω/τ-free); **F-WL-2** energy = two independent routes (pre-elim EOM vs
§2 bulk acoustic power); **F-WL-3a** kernel-orientation from the pre-elim closure; **F-WL-3b** causality
removed-record presence control bites (dead A-A dropped); **F-WL-3c** grazing classifies the ACTUAL q⁰
sound-cone block (rank 3 non-grazing; moves under a form flip; chi5-tuning removed, dead alias dropped). No
Scope regression (byte-identical). `.out` regenerated from the repaired tree (stderr dropped → no
config-warning FORMAT_ISSUE). ⚠ Review: Agent full ablation leg = sound; Grok's FULL ablation leg was blocked
by **3 spurious kills of long Mathematica jobs** (not OOM/watchdog; the Agent path completes, grok CLI does
not) — covered for this narrow fix by the Agent leg + Grok source audit + orchestrator diff-verification + the
build's scratch demo. ⚠ OPS: long Mathematica legs/builds get spuriously killed here — recover mechanically
(edits complete before the kill), and prefer the Agent path over grok for Mathematica ablation legs.

⭐⭐ **Artifact #8 (SymPy X-1 repair + comparator re-run) DONE — `53fcd98d` (build) + `fba6a34c` (re-run)**,
directive `b4c02381`. The energy-basis independence judgment now honors §5's total-divergence equivalence
(EL-signature rank at `sympy_audit.py:623`): the constructed 11 invariants reduce to a 10-dim quotient. The
redundant `(∇·u)²` is eliminated by REWRITING it into the retained curl²/strain invariants
(`(∇·u)² ≡ (3/2)st² − (3/4)curl²`), folding its `B_div` coefficient in — ⛔ NOT deleted; `ENERGY_BASIS_COUNT`
emits `Integer(10)`. Two build legs (fresh Agent + Grok, parallel): each wrote its OWN O(3) + Euler–Lagrange
enumeration and independently derived 11 pointwise / 10 mod-divergence (redundant `(∇·u)²`, matching the
orchestrator's own `x1_independent_basis_count.py`); FORM ablations (rule 14) all bit (forcing pointwise → 11;
divergence-equivalent perturbation stays 10; the two-route `ENERGY_REEXPRESSION_RESIDUAL` EL derivative is 0
and fails under a wrong fold; delete-instead-of-fold changes `U_LONG` by `B_div·η²/2`). Physics preserved
(transverse operator, in-plane EOM, breathing form, `Λ_p⁰=−Λ_A⁰/ρ_m` slice), export chain intact (1663
carried, new proof rows via D1, §11 imports `is`-bound, D3 `PROVED_EQUAL`, F6 published, no task skipped).
⭐ **Comparator RE-RUN (`fba6a34c`)** against repaired PY + repaired WL: **`ENERGY_BASIS_COUNT` now AGREE**
(10=10) and **`ZPERM_SLICE_MAP` (F-WL-1) physics-agrees** (both `−Λ_A⁰/ρ_m`; only tuple-vs-Association remains).
AGREE 20→21, DISAGREE 109→108, UNPAIRED 47→71 (the X-1 proof emissions + WL controls, one-engine-only).
Adjudicated (rule 13): **ALL 108 disagreements are emission-format/naming/representation, ZERO are physics** —
102 STRUCTURE (Tuple vs Association), 3 KEY (qOut/q), 1 DIM (`DIM_THICKNESS_RESPONSE`, W₀ convention), 2 CONTENT
(`DEGENERATE_LOCI` = same locus under an `i`-factor + `qOut==q`; `ENERGY_BASIS_OMISSIONS` = the two engines'
equally-valid representative split PY-`st²` vs WL-`(∇·u)²`, same span). `FINAL_OPERATIONAL_STATUS FAIL` = the
long-standing SymPy-Tuple vs Wolfram-Association non-parallelism, not a physics divergence.

⭐⭐ **Step record (step 12) DONE — `8ddccb74`** (`steps/S11b_interface_coupling_law.md`). Orchestrator-written;
2 legs Codex+Grok (they DISAGREED — Codex the more rigorous; I adjudicated each against the emitted objects,
rule 13). Folded once, corrections: transverse `Im ω=0` is CONDITIONAL on `μ⊥=μ_R+μ_S/2≥0` (μ_S free, no
positivity §0; only DISSIPATION is unconditionally 0); breathing slice needs all three cuts incl. `Λ_X⁰=0`;
the adjudication is "no physics CONTRADICTION" (not "all format") — there are COVERAGE gaps
(`DEGENERATE_LOCI_SOLUTION` PY-solves/WL-empty; `ONSAGER_DETERMINABLE` PY-undecided/WL-solved), `DEGENERATE_LOCI`
equal only on the regular domain q≠0 (differs at grazing q=0); coeff map WL-`(∇·u)²` ≡ SymPy `B_div+2μ_S/3`;
3 KEY = Association-key-name mismatches (not qOut/q); DIM_THICKNESS = drive-normalization; background-flow fail
`|q v₀/ω|≳1`. ⚠ OPS: the step-record legs were spuriously killed twice (memory healthy) — recovered by
SERIALIZING + a mechanical emitted-tag digest so a leg finishes inside a sweep window.

⭐⭐⭐ **S11b IS CLOSED — card re-point (step 13) DONE `565b3fe8`.** The unified export-chain step "the linear
brane–bulk interface coupling law" is complete: SymPy + blind WL engines cross-checked (physics agree, X-1
corrected), step record `8ddccb74`, and card `paper/steps/S11b_interface_coupling_law.tex` re-pointed
(directive `3fcf085d`, 2 legs, folded once) → Codex build → 2 card legs (fresh Agent + Grok, both PASS with CAS:
slice map residual 0, transverse `μ_⊥=μ_R` representative) → folded once (added the frozen-wall-width freeze +
`𝒜→𝒜_±`); card stub-compiles clean. All 13 step-run steps done; all export-integrity items ✅.

⛔ **NEXT (the S11b light sector is not finished until C, now `S11c`):** **S11c** — the **non-uniform**
variable-coefficient transverse coupling (is light's confinement unconditional?), a separate later step with
its own decision list + spec (⛔ do not reuse the `S11b` slug; G1 renamed C → S11c). ⭐ **Start at the scope
doc `steps/S11c_SCOPE.md`** — it consolidates the five re-validate-in-the-decision-list requirements (tilted
faces; Eulerian vs material density; plane waves are NOT eigenmodes — ⛔ no global dispersion relation;
uniform-limit control is known-vacuous; **falsification-first** — the coefficient is bounded by bench-top
optics), the carry-ins (`O(v₀|q_n|/ω)`, the frozen-wall reconciliation), the G14 scope boundary (nonlinear
radiation audit is NOT S11c), and the "split S11c finer than S11b" lesson. Then author the S11c decision list
(2 legs, rule 7 TRIGGER). ⭐ Broader frontier: remaining force sectors → comprehensive S1–S8 substrate
register → **throat interior**. Full state in memory `project-s11b-interface-law-result`. Run-checklist:
`steps/S11b_RUN_CHECKLIST.md`.

---

## ⭐⭐ WHERE WE ARE — 2026-08-19. ⭐⭐⭐ S11 CLOSED on its computed conclusion (generic mode spectrum, gated PASS, 8 legs). Census instruments CERTIFIED `cbc49029`; the exhaustive Q8a/Q8b strata audit's 917 under-decisions + 7 defects are DOCUMENTED CAS LIMITATIONS — the engine round to zero them was NOT run (rule 11, user's call 2026-08-19).

⭐⭐ **Why S11 closes without the engine round** (`steps/S11_stray_longitudinal.md` strata-audit closure; ⚠ **corrected 2026-08-19 after an independent Codex check** refuted the first draft's "100% strata" taxonomy): the conclusion — the generic `M_ij` eigenvalues (transverse `μ_R/ρ_br` ×`D−1`, longitudinal `B_comp/ρ_br` ×1, no cross-modulus, degeneracy exactly `B_comp=μ_R`, FORM-control unique) — is a generic-`k` result whose root VALUES are cross-engine agreed (`ROOT{1,2}_N7_RESIDUAL=0`). The char. poly is degree 3, leading coeff `−ρ_br³≠0` ⇒ **exactly 3 finite roots always**; a degeneracy merges known kernels, never creates a mode. The census findings/under-decisions are **two kinds**: **(a)** measure-zero eigenvalue-degeneracy strata (`RANK_DROP`/`STACKED_DROP`/`ROOT_COINCIDENCE`/`STRATUM`) = completeness bookkeeping; **(b)** the `KW_ZERO_LOCUS` phase-matching family (592 lines, 88 MAIN) = `c_L=c_s0`, the `k_w=0` grazing threshold (Move 6) — **physics, but the bound-vs-radiating question S11 EXPLICITLY DEFERS to S11b**, not an S11-conclusion gap and not something the engine round would touch. Neither (a) nor (b) changes the decoupled spectrum. ⚠ Honest caveats: this is the in-plane **frozen-wall-width** spectrum (not the full observable one); positivity is trivially true from premises but SymPy's sign probe is `UNDECIDED`; separation is a `D=3` fact; the `PASS` gates live in the acceptance harness, not the engine `.out`; the 7 register defects are engine **logic** bugs, not merely CAS-hard. Certified instruments + all runs preserved (`~/.s11_build/census_build4/`, `orch_verify4/`), resumable if a downstream step (or S11b's grazing case) needs a specific locus decided.

The defects round opened instruments-first (rule 4). Instrument-repair rounds, each = leg-verified defect list → two-leg directive review → fold once → Codex build → full census both records → two scoped script legs → orchestrator verification:

- **Round 2 (7 classes) — WIP `89ed80c9`**, directive `33babf8d`. Sheet mask, exact-match coverage, EmptySet-as-proof, zero-sampling-as-proof, unreachable premise check, 4 parser gaps, open reducer taxonomy. Census: 1,019→0 parse failures. Legs found 5 residual defects (incl. MY directive's premise wording and a NEW one only I caught: ⭐ **pre-substitution re/im expansion of `Element[Sqrt[...],Reals]` is branch-unsound** — both legs took the instrument's FALSE at face value; my principal-value arithmetic showed the atom ≡TRUE).
- **Round 3 (6 classes) — WIP `fd9a5835`**, directive folded `ef9085c6` (defined-union product after Codex REFUTED per-branch OR with a real piecewise candidate `k1=I·Abs(k2)`; grok had endorsed the OR — legs disagreed, computation decided). Census: WL omitted 13→3, witness failures 172→104 (⭐ ALL 74 premise-FALSE round-2 failures were expansion artifacts; the 104 survivors are all membership-driven), UNDECIDED_ZERO_SAMPLES 81→0, reducer 917/348/815 exact (both legs + me). Both round-3 legs converged on exactly ONE false verdict: D2 ROOT_COINCIDENCE false OMITTED — `simplify_residual` classifies a failure-to-simplify nested-radical EXACT ZERO (minpoly `_x`, ~1e-161) as NONZERO (`s11_census_math.py:550,:570`), and the sampled coverage fallback treats that as refutation.
- **Round 4 (1 class) — ⭐ CERTIFIED, WIP pin `3f229bbf`** (chain `90ab5e2d→89ed80c9→fd9a5835→3f229bbf`; user's ontology commit `254df530` interleaves harmlessly), review prompt `b6452518`, directive folded `fd38ebba`. Fix: constant zero-status by EXACT certificate (minpoly mandatory; `simplify_residual`→`_constant_zero_status` `s11_census_math.py:606`, 2 s guarded, degrades only to `UNDECIDED_EXACT_ROUTE_EXPIRED` never false-NONZERO); sampled coverage under per-branch AND (`COVERED_SAMPLED`/`COMPLETE_FACTOR_COVER_SAMPLED` RETIRED — a sample never proves coverage). Canonical census = **run #3** (optimized exact route, low-contention sequential): the builder rejected a parallel WL D2 **timeout** as false repair-evidence and re-ran clean (rule 2), preserving the discarded runs as `.parallel.*`/`.pre_exact_predicate.*`. Exactly **4 completeness verdicts changed, all justified**: WL D2 ROOT_COINCIDENCE `OMITTED→COMPLETENESS_UNDECIDED` (the target); 3× `COMPLETE_FACTOR_COVER_SAMPLED→COMPLETENESS_UNDECIDED` (WL D3_ROOT2/ROOT3, PY D4_ROOT1). `failures=917` & `ROUND_FAIL` unchanged; findings 348→347. **TWO LEGS** (fresh Agent + Grok, independent methods): NO finding survives the filter, **both Required-1/Required-2 instrument-copy ablations BITE**. **Orchestrator anchors** (`~/.s11_build/orch_verify4/`, rule 13): D2 candidates zero the coincidence discriminant exactly (honest undecided, not a masked omission); D3/D4 both families of both records on-variety = genuine omitted families (round-3 leg's cand[1] AND→OR **refuted**); reducer 347/816/917 object-level, sheet-excluded.

⭐ CERTIFIED (both legs + orchestrator, rounds 2–4): reducer object-level arithmetic; population reconciliation 0-gap; all 14 repair classes live by per-class ablation; every plant able-to-fail at the prior commit; every verdict transition traced to a repair class; both round-4 Required semantics load-bearing. ⭐⭐ GENUINE CENSUS CORES (now certified measured facts feeding the engine round): **917 decided-undecided (BOTH ENGINES UNDER-DECIDE)**; 136 WL + 35 PY = 171 spurious (XKIN STACKED_DROP/rank-drop; sampler-decorrelation catches, 4 new D3 tags); 104 membership-driven witness failures (= register `_REAL_WITNESS` defect); **omitted 72 = D3/D4 ROOT_COINCIDENCE (2 genuine families each) + 70 PY**; D2 ROOT_COINCIDENCE now certified honest `COMPLETENESS_UNDECIDED` (NOT an omission).

⛔ NEXT — the actual **physics frontier**, NOT more S11 machinery. DEFERRED-but-preserved (reopen only if a downstream step needs a specific degenerate locus decided — none does today): the engine-repair round (7 `DEFECT_REGISTER` entries + certified census facts), the frozen comparator contract T7 (`RANK_DROP_JOINT` WL 16-branch vs PY 2-branch complex-k), the obligation-4 census instrument. ⭐ **FRONTIER**: remaining force sectors → comprehensive S1–S8 substrate register → **throat interior** (`project_throat_interior_is_the_real_front` — the real front; lepton tower falsified). Archives: `~/.s11_build/census_build{,2,3,4}/`, `orch_verify4/`, `census_repair{,3,4}_leg_scripts_archive.tar.gz`, all leg reports `~/.s11_build/*_leg.txt`.


⭐⭐⭐ WL FIX ROUND 2 CLOSED `a4cf6539` — THE CANONICAL WL RECORD IS 21/21 COMPLETE, SPLICED, AND LEG-REVIEWED. Build self-reported FAILED per obligation 5 (final-code D4 memguard-killed at STRATUM10, kernel RSS 24.9 GB against the 1 GiB floor; the earlier "complete" run had only ~122 MB margin — RSS sidecar profiles measured the accumulation both times, ~17.6 GB entering the strata, +1.3 GB/stratum, peak ~25 GB). ⭐ USER RULING (2026-08-16): a guard-floor margin kill is an infrastructure limit, NOT a script failure — D4 re-run isolated with the floor at 256 MB (variant runner `~/.s11_build/fix2_final/run_guarded_cell_floor256.sh`; pinned original untouched): COMPLETE guard=NONE (5,412 em, 4,579 s, min available 755 MB — it dips ~290 MB BELOW the old floor, so the old floor was the only killer). D2 followed (967 s, 1,898 em, 7 GB margin). Final-engine records verified vs pre-patch runs: tag sequences identical; all payload diffs = the time-vs-memory expiry race (semantic content equal) or the final patch replacing arithmetic-over-`Failure` with clean propagation. Splice: canonical `mathematica/out/S11_...audit.out` = 16,587 lines, 21 RUN_PAIRS, untouched regions byte-identical (verified by diff). TWO SCRIPT LEGS run serialized on `d6316978` (fresh agent + grok-4.6, identical committed prompt `S11_wl_fix2_script_review_prompt.md`): round-2 diff SURVIVED — fallbacks uniform at every call site, measured-expiry-only selection, no identity gates, form ablations moved physics (721/455-line diffs). Obligation-6 HALT (44→16/12 arity on `ROOT{2,3}_RANK_DROP_JOINT`) DISCHARGED by both legs + orchestrator: same solution content, admissible sets identical, coverage improved (20 UNDECIDED→decided), the added `{sRho->0}` is a cleared-denominator artifact EXCLUDED downstream. ⭐ FOUR NEW DEFECTS REGISTERED (entries 4–7, every one verified by orchestrator re-computation; none is a round-2 regression): omitted `Solve` branch in `ROOT_COINCIDENCE_COEFF` (confirmed `{muR==bComp, sRho==1}` residuals 0); spurious branch `{muR->bComp}` in `RANK_DROP_JOINT` failing 1/16 minors both sheets (engine's own STRATUM2 rank=CONSTANT 3 corroborates; grok's `{sRho->1}` claim REFUTED 16/16; ⭐ SymPy emits a 2-branch complex k-only object for the SAME tag — the T7 comparator will flag this family); `_REAL_WITNESS` quantifies away non-solve parameters (witness fails emitted equations at generic k); latent route gaps (empty-factor-support degradation, 0 live occurrences; missing SOLUTION_ATTEMPTS on 3 secondary paths). ⚠ OBLIGATION 4 NEVER RAN: both instruments (WL probe census, SymPy containment) parse-fail on the real record format (0 probes ever executed; probe census exits 0 while dying) — calibrated only on planted synthetic records; the legs' independent containment substituted for the HALT loci + 5 cheap loci only. The instrument repair belongs to the defects round. Archives: `~/.s11_build/fix2_final/final_engine_complete_runs.tar.gz`, `fix2_d4_rerun_memguard_death.tar.gz`, `fix2_first_complete_runs_snapshot.tar.gz`, grok report `fix2_final/grok_script_leg.txt`, agent leg scripts `/tmp/claude-1000/s11_fix2_leg_agent/`. NEXT: (a) the registered defects' round (now SEVEN entries, incl. the obligation-4 instrument), (b) the FROZEN COMPARATOR CONTRACT T7 (semantic compare for expiry payloads + additive provenance + WL-only strata + entry-5's known cross-engine divergence), (c) S11 closes.

<details><summary>Round-2 build history (superseded detail)</summary>

Prior top block: BUILD MID-FLIGHT, ALL 21 WL CELLS COMPLETE FOR THE FIRST TIME (uncertified). Diagnosis (2 analysts, verified): D4 died in the unbounded `Solve` at wl:312 (one terminal call; the object decides in <1 s by exact routes); D2 died in `assumedRank[stacked]` (wl:1076, `FullSimplify` ZeroTest) after ~1,223 s of unbounded-call ACCUMULATION; six unbounded `Solve` sites + zero `MemoryConstrained` existed engine-wide. Brief folded ONCE `a8f26909` after two BLOCKING legs (Codex 7 + grok 4 findings, 7 defective-repair constructions, all verified): hard D2/D4 completion (no measurement exception), killer classes extended to the whole live path, route-uniformity generalized to fallbacks (+ out-of-cell exercise per branch), no return coercion (spec:245/271), armed probe (live operands on undecided count-class records, no-starvation budgets, planted-record calibration, residual recompute, SymPy completeness containment), guard sha256 pinned, RSS instrumented. THIRD defect REGISTERED (additive count-payload provenance extension — deliberate, round-1 fold). Codex builder launched `b5a9cbae`; measured so far: D3 269 s / 2,501 (byte-shape = baseline, rerun reproducible), ⭐ D4 COMPLETE FIRST TIME EVER (rc=0, guard=NONE, 5,412 emissions, 4,603 s), ⭐ D2 COMPLETE (rc=0, 1,898 emissions, 964 s — faster than its death run). ⚠ WIP SAFETY SNAPSHOT `14b2d56e` (user-requested) holds the repaired engine — NOT certified, no legs yet; first-complete records archived `~/.s11_build/fix2_first_complete_runs_snapshot.tar.gz` (guard runner OVERWRITES per stem). ⛔ REMAINING: builder's acceptance outputs (probe census after planted-record calibration, manifest census, partial-record compare, 19-cell byte regression, RSS profiles) → orchestrator verifies the DELIVERABLE (not the self-report) → archive builder scratches → TWO SCRIPT LEGS (fresh agent + grok, SERIALIZED, form ablation) → verify findings → scoped legs on any post-fold delta → splice canonical .out → real commit replaces the WIP label → STATUS/memory. Then: (b) the THREE registered defects' round, (c) the FROZEN COMPARATOR CONTRACT (semantic compare for expiry payloads + additive provenance fields + WL-only strata surface: SymPy `STRATUM_ORDERING` is EMPTY on all 21 cells — measured).

</details>

⭐⭐ **Full state: `research/pde_ledger_v3/REBUILD_HANDOFF.md`, top block.**

⭐ **`directives/S11_sympy_build_directive.md` — CLOSED `29e0e1d7`**, 57 lines, ⭐ 5 review rounds / 10 legs.
⭐ **`F9` in `directives/S11_export_chain_decisions_v2.md` is the naming rule** — the user's decision, and
⛔ it never changed across any round.
⭐ **Every step writes its ledger, S11 included.** ⛔ *"S11 writes no ledger"* is REVERSED.

⭐ **(done — the engine was built and has had three fix rounds; see the S11 PY ENGINE block below.)**

⛔⛔ **WHAT THIS COST, and it is the lesson: ~20 legs, 0 lines of engine.** ⭐ Rounds returned 5, 4, 4, 4
findings — ⛔ **not convergence.** ⭐ What ended it was **one census I had asserted from precedent since
round 1**: three findings collapsed into a single measured fact ⇒ **rule 2 now binds the orchestrator**
(`CLAUDE.md`), with a `_measurements/` file beside every decision doc and a **commit gate** at
`.claude/hooks/require_measurements.sh`.

⚠ **Still open, ⛔ none blocking the build:** `F4` (S10 regeneration); `F3`'s row shape; S10's four owed
items; `T7`'s comparator-side native-boolean rejection ⇒ **comparator contract, frozen before it sees
either output**.
⭐⭐ **S11 PY ENGINE — ROUND 3 CLOSED `9fb45365` (brief folded `1b2f8cf9`). ⭐ D1, D2 AND THE PUBLISH ARE
FIXED; ⛔ ONLY D4 (§Q8b) REMAINS, AND IT IS A BUILD.**

⭐ **What round 3 closed** — two legs on the folded brief (Codex + Grok: three brief defects, all
verified before folding), Codex build, two script legs (fresh agent + Grok), **no finding surviving the
physics filter from either**, orchestrator spot-probes agreeing:
- **The attempt-free refusal class is dead.** `compound_radical_present` is deleted; every locus solve
  is attempted (`emit_locus:824-839`) and an unavailability record carries the **raised failure**. The
  refused compound-radical solves return real branches — including all three `XKIN_ANISO` D3 sibling
  pairs, emitted **before** that cell's wall.
- **Status tokens decide.** Both routes, both directions (`evaluate_premise:694-719`); `MAIN` D2
  `_REAL_ADMISSIBLE` went `{UNDECIDED: 22}` → `{EXCLUDED: 12, UNDECIDED: 10}`, zero token↔test
  disagreements; surviving UNDECIDEDs carry genuinely undecidable premises. `_INCONSISTENT` is
  `Eq(locus_conditionset(equations, variables), EmptySet)` (`:858-861`) — one-sided corruption moved it
  with the equations **32/32** and left it **byte-identical** under solver-payload corruption.
- **`MAIN` publishes without waiting.** Publish fires right after `MAIN` with a publish-time cell-state
  record (`NOT_YET_ATTEMPTED_AT_PUBLISH_TIME` ≠ failed ≠ completed) that cannot be mistaken for the
  spec's post-sweep `RUN_PAIRS`/`SKIPPED_PAIRS`, which still emit after the sweep (spec `:1038`).
  Publish failures attribute to `PUBLISH`; a failed run still emits §10.
⭐ Form ablations moved **150–320 tags**; no typed payload found. Round-2 wins intact.

⭐⭐ **ROUNDS 4+5 (2026-08-14) — the XKIN_ANISO walls, each measured then removed.** Sweep-2 ground 8+ h
in D4 ROOT2's §Q4 radical null-space (undecided zero tests → no early exit → non-terminating `inv()`);
round 4 `ae105530` fixed it (conjugate-norm zero decision in the SHARED Q4/Q8 machinery, Cramer pivots,
swell-free det route, exact reducer at both Q4 sites) — the 8-h chain now runs ~98 s, D2 byte-stable
263/263, statuses only IMPROVE (28 UNDECIDED→PROVED_FALSE with rational witnesses). Sweep-3 then ground
5+ h in D4's five 4×4 radical STACKED minors (`det(bareiss)` intermediates square per level; the first
minor is ≡0 by construction); round 5 `9e392206` fixed it (guarded single-radical t-lift; dual oracles —
633 specialization checks + 138 Laplace cross-checks — each shown able-to-fail; byte-identical below the
cliff). ⭐ Process lessons that now BIND: `_measurements/` twins are GENERATED by script, never
transcribed (5 fabrications caught by legs in one hand-written twin); acceptance operands are DISCOVERED
at runtime with values verified by INDEPENDENT ORACLES, never against a route that cannot complete.
⭐⭐ **ROUND 6 (2026-08-15, overnight) — the deferred KW radical-locus wall, measured then removed
`a7e0d026`.** Sweep-4 went silent 7+ h inside D4 ROOT2's KW zero-locus `sp.solve` (5 unknowns,
`manual=True`): two analysts converged — the system loop solves its first symbol in <1 s, then spends
the hours in `checksol`/`expand` (~99% of profile) computing per-radicand-symbol solutions **its own
discard rule then eliminates**; the hang reproduces cold even at D3. Two legs on the brief found my
acceptance greenable a SIXTH time (empty-return vacuity on the only unsolved cell; soundness without
completeness; the analyst scratches physically holding the answer — archived out of /tmp before the
build). The guarded repair (single quadratic radical, radicand excluding exactly two solve variables,
canonical NOT_APPLICABLE) solves the stuck operand in ~2 s. Both script legs clean: guard census over
all 357 emitted loci moves exactly the 7 KW tags; 6/6 completed repaired cells byte-identical; both
unsolved D4 operands verified by dual corruption-validated oracles + two independent enumerations with
bidirectional set equivalence; **full fresh D3 cell 391/391 byte-identical (1117 s)** — a scoped cold
rerun is now measured-viable; D2 invariant and 2.4× faster (217 s vs 527 s).
⭐⭐ **RESOLVED 08-15 AM (user): sweep-4 killed** (8.5 h silent in the now-fixed stage; partial =
`~/.s11_build/sweep4_partial_through_D4R2_KWEQ.out`), **D5 dropped** (`90acafa0`, spec + both
engines), **sweep-5 rerun → COMPLETE `19591194`**: 21/21 cells, `SKIPPED_PAIRS` empty, zero cell
exceptions, 2h41m; the formerly 8-hour KW stages emitted in seconds; §10 + `RUN_PAIRS` verified.
⭐⭐ **THE WL SWEEP — single-kernel mode DOES NOT FIT the 30 GB box** (kernel killed by the
available-memory guard at 23 GB inside XKIN D2; partial `~/.s11_build/wl_sweep1_partial_singlekernel.out`).
⭐ **Per-cell mode is the architecture now**: the engine's own argv contract (`wolframscript -file <wl>
PACKAGE D`) runs one cell per fresh kernel, sequential (2-seat licence), driver
`~/.s11_build/wl_percell_driver.sh` concatenating into the canonical `.out` (engine output only;
driver bookkeeping in `wl_percell_driver.log`). Measured: all 18 non-XKIN cells in ~6 min; XKIN D3
complete FIRST-EVER (285 s); ⛔ **XKIN D2 is the hole — its own working set exceeds ~28 GB from a COLD
kernel** (killed at 21 min; so it is the CELL, not accumulation); XKIN D4 fits so far (~10 GB flat,
100% CPU, multi-hour radical-loci stages — rank-drop ROOT2 took ~3 h and completed with a JOINT
solution; stacked-drop, ROOT3's blocks, the KW loci, and the 8 STRATUM audits remain).
⛔ D2 NEXT: analysts decide **bounded working set vs runaway expansion** (the SymPy history warns all
three of its "needs more memory" walls were algorithmic) → bounded ⇒ rent a big-RAM box for ONE cell
(user open to it; ⚠ verify Mathematica licence activation there); runaway ⇒ a WL fix round via the
usual brief→legs→build→legs. ⚠ Until D2 exists in the WL record, SymPy's XKIN D2 is single-engine —
unverified by our own standard.
⭐ Killed partials preserved:
`~/.s11_build/sweep2_partial_through_D4R1.out`, `sweep3_partial_through_D4R1Q8.out` (first-ever D4 ROOT1
complete Q8 record). ⚠ Working-tree `S11_exports.py` = sweep-3's honest re-publish (round-4 normal
forms, values identical); commits with the successful sweep. ⛔ The `.out` is committed only from a
completed run.

⚠ **Recorded below the physics filter (both legs):** a latent **loud** `TypeError` if a live branch ever
makes a relational premise's comparison non-real (`premise.subs` at `:693`; never observed live); ⚠
`_INCONSISTENT` stays honestly UNDECIDED at regular loci even where the engine's own `PROVED_NONEMPTY`
witness entails `PROVED_FALSE` — an under-claim, never a mis-claim; pre-existing `PUBLISH:
unclassifiable free symbol G_1_1/G_1_2` diagnostics, untouched by this diff.

⚠⚠ **THE PREVIOUS ENGINE DID NOT "MISS" §Q8b — ⛔ there was no spec.** `S11_SHARED_PHYSICS.md` first exists
at `f49a1684` (*"the shared spec it never had"*), **this cycle**. ⭐ Measured on the as-built copy:
**955 lines, 41 emitted objects total, and 0 mentions of `STRATUM` / `ADMISSIBLE` / `RANK_DROP` /
`ROOT_COINCIDENCE`.** ⇒ ⛔ S11's earlier results were far less verified than they looked: the mode count on
exceptional strata was **never computed**, ⛔ not computed and found fine.

⭐⭐ **MEMORY WATCHDOG — ⛔ KILLED 2026-08-13. ⭐ RELAUNCH AT THE START OF AN S11 RUN, ⛔ never leave idle.**
`/tmp/s11_watchdog.sh` kills any S11 python over **6 GB RSS** (healthy runs sit at 120–460 MB; the blowup
hit **22.3 GB** and would take the user's other `codex` sessions with it).
⛔⛔ **Why it must NOT sit idle (user, 2026-08-13):** the hazard exists only *during* a run, and an idle
watchdog holds a **hardcoded threshold and process pattern**. ⚠ When it fires it produces a **`SIGKILL`
with no traceback** — ⛔ **the exact symptom behind four wrong diagnoses in one day**
⇒ [[feedback-measure-the-process-that-works]]. ⭐ An armed silent killer is a landmine, ⛔ not a guard.
⛔ It is NOT the memory cap the user rejected — that one sat inside the acceptance path.
⚠⚠ **To stop it, ⛔ NEVER `pkill -f <its name>` from a shell whose own command line contains that name** —
⚠ measured 2026-08-13: the pattern matched my own shell and killed it mid-command, so the edit and commit
that followed silently did not run ⇒ [[feedback-background-process-launch]]. ⭐ Kill the captured pid.

⭐⭐ **§Q8b IS BUILT AND LEG-CLEAN** — directive `20607fe6` (2 legs, 7 findings folded; the acceptance
instrument is the measured dormant defect: *an ADMISSIBLE entry did not promote*), build `94e14e42`
(658 lines; 2 legs: point-invariance byte-exact, equation corruption moves 162–800 tags, form ablation
moves ordering+dispositions, publish carries 1046 stratum rows row-for-row), fix 1 `36af7d95` (1 line:
`VALUE` reports what was **obtained**, not what constancy decided; brief hardened by 2 legs —
bidirectional equality + a non-vacuous driven witness — then 2 script legs clean, tag diff exactly the
six `VALUE` fields).
⚠ **Live cells still promote zero strata** (all candidates EXCLUDED/UNDECIDED — that is the measured
physics, not a defect); the machinery is verified through driven `/tmp` promotions covering all three
source families. ⚠ Below-filter observations recorded in the leg logs: a latent loud NaN `TypeError`
under artificial corruption; an ADMISSIBLE component whose exact point fails drops with only an ISSUES
line (unexercised); 4× exact determinant fallback recorded on driven component matrices.
⭐ Runtime ruling (user, 2026-08-13): **a cell may exceed 600 s while its output streams** — a long
SILENT stretch is the failure. `XKIN_ANISO` D3 (>600 s, streaming) is a cell to run through, ⛔ not a
wall to design around.

⭐⭐⭐ **THE LEDGER: `scripts/S11_exports.py`, committed `2f643ec3`** — 5.0 MB, **1663 rows**, every row
imports and restores; build digest aboard; `publish_time_cell_states` says exactly what was true (4×
MAIN `COMPLETED_AT_PUBLISH_TIME`, every control cell `NOT_YET_ATTEMPTED_AT_PUBLISH_TIME`); honest empty
`stratum_ordering_d2..d5` rows; ⛔ the spec's post-sweep objects are correctly ABSENT from the export.
⚠ The producing run (first real `run()` ever) was **stopped by the user after the publish**, mid
`XFORM_CURLONLY` D2, zero cell exceptions. ⇒ ⛔ **the working-tree `.out` is a partial real stream
(3.2 MB, MAIN complete) and stays UNCOMMITTED until a sweep runs to completion** — only the run record
and §10 depend on that, ⛔ never the ledger.

⭐⭐ **THE WL SIBLING ENGINE IS BUILT (BLIND) AND LEG-CLEAN** — directive `b8395704` (2 legs, 8 defects
folded pre-build), build `46ba77c2` (1530 lines from the spec alone; Codex; empty-dir byte-identity,
strace-verified zero run-time reads, form ablation 352/430 moved with all 78 non-movers legitimately
invariant, both Q2 routes independent under one-sided corruption, both legs' independent oracles match
spectrum + census), fix round 1 `45042d55`+`a4ee881c` (4 items: undefined `RationalQ` → §Q6 inert on
radicals; Q8b status typed by structure; silent point-drop; `SameQ[Reduce,False]` forging PROVED_FALSE
from a solver non-answer — the last PROMOTED from a below-filter footnote after a leg measured the
"unreachable" premise false, 14/79 radical `_INCONSISTENT` families). 2 legs clean; regression diff
430/430 tags, 25 moved, every one item-2 scope. ⚠ Below-filter, step record: `realStatus` demotes
proved-nonempty-without-witness to UNDECIDED (under-claim); strata emit last in a cell.
⛔ `mathematica/out/*.out` is still the OLD engine's — replaced only by the orchestrator's WL sweep.

⚠ **Pending, user's call on timing — the two full sweeps (both orchestrator-run, hours each):**
the SymPy rerun (post-sweep record + §10 + a committable `.out`; `XKIN_ANISO` D3–D5 scale 4–5×/dim,
streaming) and the first WL sweep (same shape; ⚠ WL `XKIN_ANISO`/`XFORM_EXTRA` cells hit 6 GB RSS at
D2 in ~4–5 min under the build's kill criteria — expect kills or a bigger budget decision at D≥3).
⭐ THEN the frozen comparator contract (`T7`) — written before it sees either output.
⚠ Logs `~/.s11_build/` (WL under `wl/`, evidence `wl/evidence_r1/`); Q8b probes `~/.s11_build/q8b/`;
the reproducer is `~/.s11_build/repro_d5.py` (~353 s) — ⛔ never run the full loop to test a fix.

⚠ **HYGIENE, ⛔ NOT A CONTROL:** `/tmp/s11*_leg_*/`, `/tmp/f9*_leg_*/`, `/tmp/s11_fold_leg/` hold review-leg
scratch that computed real S11 physics. ⛔ Do not commit it into the tree — it is scratch, and the tree is
the record. ⛔⛔ **It is NOT quarantine and must not be turned into any: do not move it, do not hide it from
a builder, do not claim blindness from it.** ⚠ **Measured 2026-08-12 — this exact line, in this exact file,
is what made me move 336 files and then reverse it**, the third occurrence of a mechanism `CLAUDE.md`
rule 12 cut. ⇒ the only blindness control is
`research/pde_ledger_v3/directives/S9_export_chain_rebuild_directive.md:17`.

---

## ⭐⭐ (superseded) WHERE WE WERE — 2026-08-10

⭐ **S11's shared spec is REPAIRED and gated** (`ab8cb50e`, 914 → 1005 lines): an inertia control, `C16`'s
stacked-matrix stratum source, `Q6r` repointed off the deleted `reduction/`, `Q3`'s multiplicity object
split, and a closed live-read exemption. ⚠ Two legs on its decision list, two more on the repaired file.

⚠ **Trying to launch the PY rewrite surfaced four defects in the export chain itself** (`C17`–`C20`), and
the first two routes proposed for them were **both wrong**.

⭐⭐ **THE ROUTE — `directives/S11_export_chain_decisions_v2.md` (`4d81e9de`).**
⛔⛔ **`S11_naming_and_chain_plan.md` is SUPERSEDED. ⛔ There is NO S10 retrofit.**
- ⭐ **Keys stay FLAT; `D5` is unchanged.** Before writing a key that exists in the imported `LEDGER`,
  **compare the OBJECT** — same object ⇒ re-derivation (both operands + residual, then guard); different
  object ⇒ ⛔ **fails loudly.** ⚠ `DEFECT_REGISTER.md:675` had already prescribed exactly this.
- ⭐ A re-derived row carries its **evidence in the row**; S10's export is **REGENERATED**, ⛔ not re-keyed,
  and its tag names do not move.
- ⭐ `C19` is a **real deviation** (`S10:197` orders *"Emit `M_A`, `M_B`"*) — the record **discloses** it and
  the rename is its **own gate**; ⛔ S11's build does not depend on it.

⛔⛔ **WHY THE FIRST TWO ROUTES FAILED, and it is the most useful thing on this page.** *"Retrofit S10 first,
S11 binds to its keys"* — ⛔ measured false: S11 spells **two** S9-origin knob rows and a pointer. *"The
overwrite violates rule 2"* — ⛔ measured false: the engine emits both operands and the residual at
`S10…sympy_audit.py:2089-2111`. *"Producer-scope the keys so collisions are impossible by construction"* —
⛔⛔ **that IS the defect**: under `D5` two steps deriving one object meet on **one key** so they can be
compared; scope them and **nothing compares them** ⇒ ⭐ rule 6, *don't make divergence impossible*.

⚠ **Next concrete step:** the **S11 PY decision list**, rewritten against `F1`–`F7` and its own five blocked
findings — ⛔ two legs before any builder moves.

## ⚠ [HISTORICAL SNAPSHOT 2026-07-31 — the LIVE front is the top block (S11b CLOSED; NEXT S11c). Read below only for v3-rebuild backstory; ⛔ its "S11 NEXT / A,B to rebuild / S11b-C never built" is stale.] v3 ledger opened

⛔⛔ **FRONT DOOR CHANGED 2026-08-04 — `research/pde_ledger_v3/REBUILD_HANDOFF.md`.** Read that first.
**ALL NEW LEDGER PHYSICS IS STOPPED**, S11b-C included, until the script rebuild closes. ⚠ The engines were
measured to emit physics conclusions as **typed sentences with no computation behind them**, at named lines
in **three independently-built steps**, and **eight review legs missed it.** ⛔ Do not build on a v3 script
result without checking it there first.

⇒ **Rebuild progress:** S9 ✅ · **S10 ▶ record `e167b07f` + PAPER CARD DONE** — chain
`e644876c`+`c84263ed`, comparator `82443c95`, record over 4 builds / 14 legs, card over 2 builds / 8 legs.
⚠ **S10 still owes: the D12 naming pass, the registers.** · ⭐⭐ **S11 ▶ NEXT — full rewrite, same pattern
as S9/S10** · S11b-A, S11b-B to rebuild · S11b-C never built.
⛔ Detail belongs in `REBUILD_HANDOFF.md`, not here.

⛔⛔ **S11 STARTS WITH A SPEC REPAIR, ⛔ NOT AN ENGINE.** ⚠ Its 914-line spec is **closed and incorrect** in
two ways S10's rebuild measured: it carries **no inertia control** (all 7 packages vary `W` only, so the
second structural premise S10 needed is unprobed and **uninstrumented**), and its stratum enumeration
reproduces **`C16`** — strata come from minors of `M_r`, ⛔ never the stacked `[M_r; kᵀ]` that governs the
transverse count. ⚠ **A closed spec is not a correct spec**, and a shared spec is physics-bearing.
⛔ Also: the PY engine **will not run** — `from registry_read` at `:21`, and `reduction/` is deleted.

⛔⛔ **THE LIGHT SECTOR'S REMAINING QUEUE IS PAPERWORK. ⭐ The physics that moves anything is `S11c`** (was
S11b-C; S11b's uniform half is now CLOSED) **and `S8`** — ⚠ `S11c` is the **MacCullagh differentiator** and
is still **unbuilt** (scope: `steps/S11c_SCOPE.md`); ⛔ `R-S8-04` decides
whether the curl-only functional is an **admissible continuum mechanics at all**, and **S8 has never been
run** ⇒ `SUBSTRATE_REQUIREMENTS.md`, `V3_STEP_PLAN.md` PHASE 4b.

⚠⚠ **OUR FAILURE HAS A NAME — the *consistent comparison problem*, 1986** ⇒ `docs/method_prior_art_findings.md`.
⛔ No amount of spec care closes it.
⛔⛔ **CORRECTED 2026-08-06 — ⛔ do NOT join on a fingerprint.** ⚠ An earlier version of this line said
comparison would move to joining on a **fingerprint**; that was **retracted**, and the counterexample is
in our own tree (`S11bB…:539` — every passing row emits an identical payload, so repeated zeros, dimension
vectors and booleans all collide **deterministically**). ⭐ **A fingerprint is an EQUALITY ORACLE for
objects ALREADY PAIRED; it cannot discover which object corresponds to which.**
⭐ **Correct form:** a stable shared `quantity_id` as the **join key**, with the exact evaluation vector as
a **typed equality oracle after pairing**. ⇒ what that deletes is the symbol-spelling negotiation layer,
⛔ never the tag inventory.

Then `research/pde_ledger_v3/NEXT_SESSION.md` for everything else.

- ✅ **Step ① DONE** (`407eed94`) — the `a`-pin is retired from everything that computes. It was never a
  wrong *number* (`ħ/(m c_s0)` **is** a healing length in a standard convention); it was a **category
  error** — a medium-wide constant used as a per-defect quantity.
- ✅ **v3 opened, sector-scoped: gravity · light · gravitomagnetism · charge · magnetism.** Same
  walkthrough method, a sector boundary instead of a phase ordering. ⛔ Not a third method change.
  ⚠ **Widened 2026-07-31 (user decision, `c13f9329`)** — charge and magnetism share the **same
  unsolved object** as gravity: one nonlinear throat solve (`docs/model_map.md#shared-r1-throat-solve`).
- ⭐⭐ **Steps are now walked SIDE BY SIDE with the user** — participation in the derivation, not
  approval at the end. v2's delegation is why three structural findings went unnoticed.
- ⚠ **The scope boundary is amended and may not hold.** The worldtube result is response-side,
  conditional on compactness; the Spin Problem says a compact defect cannot give correct frame dragging.
- ⛔ **Open and out of scope: what makes a muon a muon.**

⇒ Known error surface: `research/pde_ledger_v3/DEFECT_REGISTER.md`.


**A thin pointer, not a copy.** Per-stage detail belongs in the per-Part split docs and per-stage
notes. History is in git. If this file starts growing narrative, cut it.

---

## ⚠ ▶ SUPERSEDED — v2 position, kept for provenance (2026-07-30)

> ⛔ **HISTORY, NOT CURRENT.** The live position is the v3 block at the top of this file.
> ⇒ `research/pde_ledger_v3/NEXT_SESSION.md`.

**Current front: the DERIVATION WALKTHROUGH** (user decision, 2026-07-30) — one derivation step at a
time, forward from the medium's defining properties, recording at each *what it is · what it does ·
what's new*. The irreducible input count accumulates **by construction** rather than being inferred
backward from finished artifacts. ⛔ **Its method has ONE canonical home — do not restate it here:**
`docs/derivation_walkthrough_plan.md`. Step records: `research/pde_ledger_v2/walkthrough/`.
Read-first handoff: `research/pde_ledger_v2/_scratch/NEXT_SESSION.md`.

**▶ Position: phase 0, two steps done** — `walkthrough/00_medium_and_brane.md` (the medium, and the {#position}
brane as its ordered state) and `walkthrough/01_sound_speed.md` (the sound speed).

⭐ **Why the method changed — recorded so it is not re-proposed:** the backward census produced real
findings but made the physics unfollowable for its one reviewer, and after eleven commits no physics had
been verified. ⚠ This was the **second** occurrence of a failure `docs/development_pipeline.md` already
records (apparatus growing above the physics).

**Continuing BEHIND it: the DIMENSION REWRITE — 7 of 30 scripts converted** (stage004, 011, 012, 013,
016, 018, **023**) — and ⭐ **all seven are WAIVER-FREE**: `ARTIFACT_NAME_WAIVERS` is empty, so every
converted stage compares every name it emits, with no exemptions.
⚠ That is a **coverage** statement, not a strength one.

⭐⭐ **THE CROSS-CHECK IS EARNING ITS KEEP — measured on stage016.** Relabelling that stage's basis
leaves its **own** 82 assertions completely blind (exit 0, 82 PASS, printing `measure: 'M^3'`), while
the comparator catches **18 of 21**. The comparator is the *sole* instrument between a converted stage
and a relabelled basis — which is why the empty waiver registry matters.

⚠ **`NEEDS_ADJUDICATION` in the canonical table — 3 groups, and all three are correct.** `K_eta`
(three levels: 013 line `M L⁻¹T⁻²` / 016 volume `M L⁻³T⁻²` / 023 reduced scalar `M T⁻²`), `T_Omega`
(016 volume `M L⁻³T⁻²` vs 023 reduced scalar `M T⁻²`) and `mu_eta` (`M L⁻¹` vs `M L⁻³`). These are
**REDUCTION LEVELS, not drift**, surfaced rather than hidden because the variants
were not renamed apart (§7).
⛔⛔ **CORRECTED 2026-07-30 — `a` across 016/018/023 is a DIMENSIONAL COINCIDENCE, not agreement.** This
line used to read *"`a` now groups AGREE across 016/018/023 … the same throat-radius `L`"* and log it as
evidence of consistency. ⛔ It is not evidence of anything: the three rows agree on being **lengths**,
which two different lengths do trivially. The **pin `a`** (the nullity-1 residue of imposing four unit
pins) and the **throat radius** (a physical size) are different quantities sharing one symbol — a
same-name-different-quantity collision at the foundation of the model, filed as confirmation.
⚠ This entry is damage class **(c)** in `research/pde_ledger_v2/walkthrough/DECISIONS.md` **D-01a**;
⚠ **RESOLVED 2026-07-31 by removal** (`407eed94`) — the pin relation and quantity are deleted; `a` now
means the throat radius. ⛔ There is no open registry class to decide.
⚠ The route caveat below still governs every other group: a green comparator is **agreement
between implementations, not independence of route** — per `research/pde_ledger_v2/manifests/DIMENSION_REWRITE.md`
§1 (**THE CHARTER**, *"a green comparator shows that two implementations agree, not that they were
reached by two independent routes"*), and git establishes no derivation order for any engine pair.
⚠ `T_w` does **not** group (016 `T_w` → `TW` vs 013 `Tw` → `Tw`), so
that line-vs-volume debt stays invisible — a measured consequence of naming debt, not a theory.

⭐ **TRACK TWO COUNTERS, NOT ONE — they are different finish lines.**
- **converted** — on the shared module, cross-engine gate green: **7 of 30**.
- **physics-verification evidence** — ⚠ **quantity-level; there is NO defensible stage count yet.**
  Recorded: stage012 **14 CORRECT / 0 WRONG**, stage013 **9 / 0**, all six emitted stage018
  records, the three formerly-waived 011/012 records, and ⭐ **stage016's 21 / 0 verdict**
  (**12** declared literals in both engines plus **9** computed, **0** of the 9 from any physical
  input). ⚠ **stage016's verdict is PROSE, not a per-quantity route table** —
  `research/pde_ledger_v2/notes/stages/ledger_stage016_l2_so3_covariance.md:175` (**§1.6**, the step-(c1)
  physics leg); its dimension-object enumeration is **§1.5**, not §1.6. The 011/012/013/018 results
  predate §4-c1, and **none** of these were normalised into per-stage tracked verdict tables, so
  **do not call whole stages verified and do
  not infer a "25 remaining" complement.** §4-c1 exists so this becomes countable going forward.

⛔⛔ **DUAL-ENGINE AGREEMENT IS VACUOUS WHERE BOTH SIDES ARE HAND-DECLARED LITERALS — recorded
2026-07-30, and it is the most important thing on this page.** ⚠ **State it at the width the
measurement supports** (matching `_scratch/NEXT_SESSION.md`): **stage023 emits 29 records — 22 typed
`SOURCED_DIMS` declarations plus 7 live `dim_of` walks**, and those walks run over exactly those 22
literals and no other dimensional input. **stage016 emits 21 — 12 typed rule-table entries plus 9
computed** by `dim_of` over real expression inputs, but again sourced only from its own declarations
(and **3 of the 9 are self-referential**, walking a declaration back to the constant that defines it).
⛔ It is **not** "29/29 literals" and **not** "zero computed" — both overstate. What holds is that **no
dimensional input enters from outside the stage's own typed declarations**, so the comparator on these
stages catches a **transcription split between two typed copies of the same numbers, and nothing
else**. ⇒ **Nothing remaining independently RE-DERIVES the physics outside one fresh agent.** ⛔ Do not soften this and ⛔ do not read a fix into it: it is stated here
because it is what the **derivation walkthrough** exists to change — a forward step derives, where a
converted stage only re-declares (front decided 2026-07-30).
⚠ **Read it beside the cross-check measurement above, not instead of it** — the comparator is still the
sole instrument against a relabelled basis, so this is ⛔ not licence to cut it.

Cross-engine agreement is necessary and **not sufficient** — proven twice: it is blind to a
same-dimension different-quantity merge, and to a shared wrong declaration. The physics leg is
what establishes correctness, and **it does not depend on conversion** — it can run against any
stage's existing declarations. If the goal that matters is Part VII's firewall being trustworthy
rather than the manifests unblocking, that is the leg to sequence around. See `DIMENSION_REWRITE.md`
§4-c1.

Every SymPy audit script's dimension handling moves onto one shared module
`research/pde_ledger_v2/scripts/ledger_dimensions.py`, **one stage at a time**, each verified by an
independent `.py`-vs-`.wl` cross-check.

> ⭐ **READ `research/pde_ledger_v2/manifests/DIMENSION_REWRITE.md`** — the single canonical doc for
> this workstream. And read `docs/model_map.md` **before touching any script**.

✅ **The shared module is no longer self-attesting** (2026-07-27) — conversion is unblocked. ⚠ What that
digest earns is a **staleness** signal only (*the module changed and the stages were not re-run*); see the
digest block below.
✅ **stage023 IS CONVERTED — both halves.** Its `.py` is on the shared module, the comparator is green
(`py=29|wl=29|shared=29|py_only=0|wl_only=0|mismatches=0`, no waivers), the orchestrator regenerated
both engines' artifacts itself, and the sealed prediction is adjudicated (**5 fully confirmed · 1
falsified (P2) · 1 split** — P3's mechanism confirmed, its exclusivity falsified by U13).
Evidence: stage note §1.6 (§4-a enumeration), §1.7 (§4-c1 physics verdict), §5.1 (steps g/g2/h/h2).
⛔ The `.wl` never joins the module: the charter is **SymPy-only** and the Mathematica side is authored as
its own route (`DIMENSION_REWRITE.md` §1), which is why the comparison is a permanent standing cross-check.
⚠ Read that as the **authoring rule** it is — §1's own honest statement governs: a green comparator shows
**two implementations agree, not that they were reached independently**, and git establishes no order.
**▶ NEXT CONVERSION** — ⚠ next *within this workstream*, which now runs **behind** the walkthrough; the order
and hazards below stand unchanged — **per the recorded conversion order (`DIMENSION_REWRITE.md` §8):** (1) the
stage027-shape decision, (2) 027, (3) 021 (heaviest). Detail and the measured validator/harness hazards
are in `DIMENSION_REWRITE.md` §8/§9.
✅ **The ablation-fixture FREEZE AUTHORITY IS RETIRED** (user decision, 2026-07-29/30) — with it, the
coupling that made nine live dimension-rewrite paths untouchable. Convert freely; nothing here is frozen,
byte-perfect or under a custody rule. See `DIMENSION_REWRITE.md` §4.
**▶ QUEUED INSIDE THE DIMENSION-REWRITE WORKSTREAM — ⛔ this is NOT the project's next build.**
⛔ **Nor is the next build the walkthrough's next derivation step.** The order is fixed by user decision
**D-01a** (`research/pde_ledger_v2/walkthrough/DECISIONS.md`):
✅ ~~archive the old apparatus~~ (**DONE** 2026-07-30, `archive/`) → ① **repair the `a`-pin damage** → ② **resume the walkthrough**
(`docs/derivation_walkthrough_plan.md`). ⛔ Do not resume deriving before ② is done.
⭐ **Both items below are
KEPT** — not because the rewrite needs them, but because the walkthrough's own checks consume them
(plan §5): the ablation driver is the spec for its **check 6** (able-to-fail with a plausible physics
mutation), the `DIM|` emitter the spec for its **check 5** (dual-engine). A queue within that
workstream, TWO items, in this order:
**(1) the ablation driver**, **RE-SCOPED SMALL** (user decision, 2026-07-29/30, `DIMENSION_REWRITE.md`
§12b(b)): mutate a declaration, confirm the declared assert fires, record it — reviewed by one fresh
agent. ⛔ No contract, no frozen fixtures, no three-session shape. Requirements (trimmed) at
`research/pde_ledger_v2/notes/ablation_driver/REQUIREMENTS.md`; the old contract and its frozen fixture
suite were **deleted 2026-07-30** (in git history if ever wanted) — ⭐ what survived them is that file's
**Appendix**, incl. the stage023 legacy mapping as a **reference** for the retrofit, but ⛔ it is **not**
authoritative and A7 does not require agreement with it.
**(2) the shared Mathematica `DIM|` emitter** (`DIMENSION_REWRITE.md` §12b, closing block) — **just write
it**, against `research/pde_ledger_v2/notes/wl_emitter/REQUIREMENTS.md` (v2, reviewed CLEAN). ⛔ No
contract, no frozen fixtures, no three-session shape. Which stage converts next and which
tooling gets built next are separate sequences; in time **both** builds land **before** 027 begins.

**▶ WHAT STAGE023 LEFT OPEN.** Stage023's tracked physics leg records 34 quantity routes and the
two scoped tallies `24/0/10` on corpus identifications and `27/0/7` inside the stage-local closure
(`research/pde_ledger_v2/notes/stages/ledger_stage023_nullspace_underdetermination.md:379-440`,
**§1.7(1), per-quantity verdict, tally, and the 34-row route table**). Its seven unresolved derivation questions now have named
work routes in `research/pde_ledger_v2/manifests/DIMENSION_REWRITE.md`, all under **§12 OPEN / DEFERRED**
(`:931`). ⭐ **The `WORK-023-*` ID is the stable anchor; the line number is a convenience** — each is the
first line of its own bullet: `:939` (**WORK-023-MOMENT-CONVENTION**), `:999`
(**WORK-023-STAGE009-MOMENT0**), `:1028` (**WORK-023-D0-SEAM**), `:1062`
(**WORK-023-STIFFNESS-REDUCTION**), `:1101` (**WORK-023-L1-L2-PROFILE-IDENTITY**), `:1125`
(**WORK-023-CS-EVALUATION**), and `:1167` (**WORK-023-SOURCED-PROVENANCE**). W3 is outside that work list because its confirmed
arithmetic correction is already folded at
`research/pde_ledger_v2/manifests/DIMENSION_REWRITE.md:648-663`
(**§8, stage023 `gU/gW` correction**); `q_free` is outside it because §1.7 classifies that record as
an unread control rather than a competing identification
(`research/pde_ledger_v2/notes/stages/ledger_stage023_nullspace_underdetermination.md:428,442-452`,
**§1.7(1), `q_free` verdict and tally explanation**). ⚠ These §12 items are derivation work, not user
gates, and they do **not** hold back the counter above: stage023 is converted on **both** halves and
is counted in the **7 of 30**.

✅ **stage016 COMPLETE** — both engines, comparator `py=21|wl=21|shared=21|mismatches=0`.
✅ The four group-A stages have `.wl` reachability verdicts — that recorded GAP is **closed**;
⚠ but those counts are **PROVISIONAL, not completeness proofs**: 016's survey got its 21 emitted
quantities right, yet missed two **non-emitted** source/control-flow cases — so the method does not
close the broader inventory. Each stage's tracked §4-a enumeration + adversarial review is what does.
027 is the awkward one (**MIXED**: its **computed** vector never reaches top level — it dies in `runAll`'s
`Module`. ⚠ **Corrected 2026-07-30:** its **16 declared per-symbol `baseDims` vectors DO reach top level**,
so 027 is **not** confined to a 1-row stage; only reaching the computed vector needs new call sites).
⛔ **Measured corpus-wide: all 43 `.wl` files end in `Exit[]`,** so a print appended at end-of-file is
dead code. Emit before the terminal `Exit[]`. See `DIMENSION_REWRITE.md` §8/§9.
⛔ **Before converting 035/036/037, fix the canonical-table generator** — it *raises* on cross-engine
axis-order disagreement. ⚠ **What was actually measured, stated precisely:** the committed stage037 `.wl`
declares **no axis order at all** — it is a formal `massScale`/`lengthScale`/`timeScale` rescaling with no
ordered vector and no `axes=` header — so the crash was **not** measured on the committed pair. The
`M,L,T` order was observed only in the **gitignored** `_scratch/spike037/` prototype (⛔ since deleted and not
retained; the order is recorded in `research/pde_ledger_v2/notes/stage037_dimension_emission_spike.md`) (`dimensionAxes` keyed
`M`,`L`,`T`), against the committed `.py`'s `(L,T,M)` (`ledger_stage037_route_b_boost_structural_relation_sympy_audit.py:604`).
⇒ The disagreement is real and will appear the moment 037 emits, but it is a **spike-vs-`.py`**
measurement, not a committed-pair one.

Per stage: emit `DIM|` records into the `.wl` → ⭐ **run the PHYSICS leg FIRST** (it is
blocking, and a NAMING decision is a physics decision) → rewrite the `.py` onto the
module → **re-run the `.py`, then** compare axis-labelled → **regenerate the `.out` and byte-compare it**
(⭐ load-bearing: an uncaught `Throw` exits 0 with an empty or truncated transcript and only the byte-compare
catches it — `DIMENSION_REWRITE.md` §9. ⛔ **No longer an orchestrator-only duty** — whoever runs the stage
does it; the old *"agents cannot run Mathematica"* reason was measured false and the
independent-party replacement was never adjudicated, so it is cut as a custody step, `DIMENSION_REWRITE.md`
§4-(g2)) → **two mutually independent fresh review legs** (each fidelity + per-tooth ablation; a stage's
math is physics-bearing) → commit.
⛔ **Retired 2026-07-29/30:** commit-the-`.wl`-before-the-`.py` reference custody. ⭐ **What replaces it and
is NOT ceremony: the co-authorship guard — the party that wrote the `.py` must not be the party that adjusts
the `.wl` until the comparator agrees** (`DIMENSION_REWRITE.md` §4-(d)); tuning whichever side disagrees is
the LLM-shortcut-that-resembles-a-pass, not a fix.
⭐ **The prediction goes in the SESSION SCRATCHPAD and is folded into the stage note after the reviews land**
(restored 2026-07-30) — because *before launching any agent, ask what is reachable from the working tree that
states your expected answer; absence beats instruction* [`never-supply-the-expected-reason`]. ⛔ That is the
whole reason: no custody, no sealing, no ordering ritual. The physics leg is still told to **derive from the
model**, never to *check a claim*.

⭐ **The hard tail is bounded.** A spike **prototyped** an independent `.wl` route for **stage037**
(`ROUTE_EXISTS`, 21/21 quantities, a real comparator failure on a seeded error), so the old
"genuinely impossible" verdict is **false for 037**. **035 and 036 have identified routes that are
NOT yet prototyped** — source inspection only. ~0.5–1 engineer-day per stage is an **estimate**, not
a measurement. See `DIMENSION_REWRITE.md` §3b.

⭐ **Python sidecars are source-hash bound, and that is exactly a STALENESS check — which is all it is
asked to be.** It catches the real, common error: *the `.py` changed and the sidecar was not
regenerated.* ⚠ It does **not** prove the sidecar came from a run — a hand-written one carrying the
right digest reaches `PASS` (demonstrated 2026-07-27). ⛔ **Retired 2026-07-29/30: the orchestrator-
regenerates-the-sidecar control, and calling this an open hole.** Forging a sidecar is a *motivated-
adversary* move, and this project hardens against **drift and honest error**, not that
(`docs/development_pipeline.md`, *THE POSTURE*). What catches a wrong sidecar is the comparator plus the
blocking physics leg. See `DIMENSION_REWRITE.md` §9.

⭐ **THE SHARED-MODULE DIGEST — A STALENESS PING (downgraded from an "authority" by user decision,
2026-07-29/30).** `scripts/check_ledger_dimensions_pin.py` compares `scripts/ledger_dimensions.py`'s
current source bytes against `scripts/ledger_dimensions.accepted.sha256`, and fails on any difference.
⛔⛔ **State only that.** It **cannot** establish *"you edited the module and did not re-run the stages"*:
it inspects no producer and no run — it compares two hashes, and `--accept` rewrites the recorded hash from
the current module bytes without checking that anything was re-run. A red digest means **the module differs
from the last accepted baseline**, which is usually a module edit whose downstream has not been refreshed —
useful, cheap, and *not* evidence about the stages either way. ⭐ **The check that actually detects a stage
not re-run after a module edit is the SIDECAR binding** (the module digest is stamped into each stage's
sidecar header and the comparator recomputes it). Stubbing `dim_residual` does trip the pin in the standalone
control, the comparator and the generator (class `MODULE_PIN_MISMATCH`, distinct from sidecar staleness).
⛔ **RETIRED FRAMING — do not restore it:** "an authority **no producer writes**", the digest as a
**validator** whose green says anything about correctness, and **re-acceptance as a review event with a
recorded reason and a second witness**. A red digest after a legitimate module edit means **refresh and
re-run the producers** (`--accept`, then stage → comparator → generator); that is a reset, not a trust
decision. Procedure: `DIMENSION_REWRITE.md` §4.
⚠ **Its bound:** it covers no stage source, no `.wl`, no `.out` and no sidecar content, and it never executes
a stage. ⇒ **Read a green digest as "not stale", nothing more.** ⛔ **Clear `scripts/__pycache__/` after any
ablation edit/restore loop** — equal-size edits (sign flips, `sum`→`min`) let CPython reuse timestamp-valid
stale bytecode by accident; that one is live and practical. *(The bytecode / `sitecustomize` / trust-root
analysis that used to sit here is cut 2026-07-30: those are motivated-adversary routes around a staleness
ping, which the governing test does not buy.)* Detail: `DIMENSION_REWRITE.md` §9.

⭐ **A bare stage run is a PRODUCER, not a validator** (user decision, 2026-07-27) — its exit code and
`PASS` tally are **not** validation evidence. The validators are the comparator and the generator (the
module digest is a staleness check, not a validator). §11 had already measured why (a relabelled basis leaves stage016's own 82 assertions blind
at exit 0); this states it as a rule, and it is what lets the pin sit outside the module.
⚠ `run_all_audits.sh` tallies `Fail: N` but **exits 0** — it gates on the pin, not on audit failures,
and it never invokes the comparator or generator.

⚠ **Before resuming, read `DIMENSION_REWRITE.md` §1b (the D1–D5 decisions) and §3b (what those
decisions REOPENED).** Several recorded conclusions — three waivers, four "impossible" stages, a
coverage estimate — were correct only under constraints since lifted, and will read as settled.

⛔ **THE DERIVED-vs-DECLARED CENSUS IS RETIRED AS THE FRONT** (2026-07-30), superseded by the
walkthrough. Its artifacts were **archived 2026-07-30** to `archive/census/`
(agreed in principle, not yet executed). ⭐ **Its four surviving findings are carried forward in
`research/pde_ledger_v2/_scratch/NEXT_SESSION.md`** (⛔ OPEN CORPUS FINDINGS: the
`parameter_register`-vs-`stage023 source map` tier-1-vs-tier-3 contradiction · stage016's false
`CONSUMED-from-011/012/013` attribution · the wrong stage016 locus in four tracked files · zero
cross-artifact citations resolving to a locus) — ⛔ **they must not be archived with the apparatus.**
⚠ Still true independently of the census: **Part VII's stage046 row in
`research/pde_ledger_v2/notes/part7_integration_atomic_split.md` requires every constant
DERIVED/INPUT/gap/benchmark**, and that requirement is unmet.

**Why the DIMENSION REWRITE continues behind it — ⭐ ITS JUSTIFICATION CHANGED (user, 2026-07-30).** The
shared import's purpose is to be **the single place holding every input and what it derives from**, so
the **leftovers are readable** — ⛔ not merely "one representation, from one place".
⚠ **Therefore "7 of 30 converted" measures REPRESENTATION UNITY and does NOT measure progress toward
that grand check:** `scripts/ledger_dimensions.py` currently carries **dimensions**, not **defining
relations**. The relations live in `research/pde_ledger_v2/reduction/` and grow one walkthrough step at
a time (`docs/derivation_walkthrough_plan.md` §5a). ⛔ **Nor is it "consistent by construction" — that
overstates.** One
shared module buys **representation unity** (one basis type, one exponent type, one recovery path); it does
**not** buy correct dimensions, because **two stages can declare the same wrong exponents through one module**
just as easily as through thirteen idioms. Correctness comes from the blocking physics leg.
Thirteen dimension idioms across 43 scripts *is* drift; the decision
(`b5527062`, `aae5d389`) was to **fix the corpus, not weaken the check**. Part VII's whole-system
dimensional firewall is **claimed** to consume the module directly rather than the manifests — ⚠ **an
assertion, not an established fact:** stage046 is unbuilt, and `notes/part7_integration_atomic_split.md`
(the 046 row) names the firewall without naming its input source. Do not cite it as settled. The manifests' semantic
core continues in parallel, trimmed (§ "PAUSED" below).
⛔ **The old justification — "the 44-stage manifest fanout is blocked, dimension recovery covers only
~16 of 43 scripts" — is RETIRED as false.** Verified 2026-07-29/30: the composite checker recovers
dimensions from **10 of 43** scripts — exactly **7** carry a `class Dim` the recovery walks (005, 006,
007, 008, 009, 030, 031) plus **3** registered bare-tuple digests (032, 038, 042) — and **all seven
converted stages carry none**, because the shared module exports `Dimension`, not `Dim`. Conversion
therefore does not raise that count, and the fanout was never what the rewrite unblocks.
⚠ **Precisely:** conversion *lowers* recovery only for the **seven `class Dim` stages**; for every other
script it merely fails to raise it, because there was nothing recoverable there to lose.

## ⏸ PAUSED (⚠ now behind BOTH the walkthrough front and the rewrite; user confirms sequencing)

- **stage 044-v2** — redo stage044 with a DYNAMICAL-Σ sleeve (un-freeze `S_hold`, commit the
  `κ_bend / κ_anchor / collar-tension` bending knobs). User-decided, 044-LOCAL.
  Anchor: `research/pde_ledger_v2/notes/stage044_v2_unfreeze_prep.md`.
- **stage 045** (VII-2b) — the non-variational drain/return block + BCs + force partition, where the
  drain-placement crux and the USER mini-gate land. Drain = the dynamical `Γ_B`; frozen-wall ruled out.
  ⭐ **That mini-gate STANDS** — it is a modelling DECISION for the user, which the reduced process still stops
  for; it is not a per-chunk gate.
  Anchor: `research/pde_ledger_v2/notes/stage045_nonvariational_block_prep.md`.
- **Manifest / integration-test system** — built + committed (`e849e303`), 4 of 44 manifests extracted.
  ⭐ **CONTINUES, ON ITS SEMANTIC CORE** (user decision, 2026-07-29/30; corrected 2026-07-30): quantity
  identity, **citation integrity**, **export/lifecycle enforcement**, dimensional relations, **the
  lifecycle census**, the dependency graph and its cycles, mutation, **genesis**, consumption completeness.
  ⛔ **The trim that dropped citation-integrity / lifecycle / genesis is WITHDRAWN** — each catches a way the
  *physics* could be wrong (a changed consumed equation · consumption of retired physics · a wrong
  irreducible-count range · calibrated-or-target-matched genesis, i.e. fit-vs-derive), not bookkeeping.
  ⚠ It is **not** what the dimension rewrite unblocks (see the front's justification above) — the two are
  independent. Docs:
  ⛔ **ARCHIVED 2026-07-30** → `archive/manifests/` (`MANIFEST_README.md`, `EXTRACTION_PROTOCOL.md`,
  `composite_build.py`, the schemas, and the 4 extracted stage manifests). ⭐ `DIMENSION_REWRITE.md` and
  `DIM_ORDER_DECISION.md` **stay** at `research/pde_ledger_v2/manifests/` — paths deliberately unchanged.
  ⛔ **UNRESOLVED against the new front:** `docs/derivation_walkthrough_plan.md` §5 marks `manifests/`
  (**except** `DIMENSION_REWRITE.md`) for **archive** as a superseded route. That is not yet executed and
  needs a **per-file split, not a directory move** — the directory holds live code beside the active
  conversion doc. ⇒ Read "CONTINUES" as *not yet withdrawn*, not as *reaffirmed*.

## LEDGER BUILD STATUS

| Part | Sector | Status |
|---|---|---|
| 0 | Conceptual | scaffolding only |
| I | Medium | ✅ built (004–007) |
| II | Gravity | ✅ COMPLETE (001–002, 008–029) — sector closed |
| III | Light | ✅ DONE = stage003 (surviving-solution rule) |
| IV | Charge | ✅ COMPLETE (030–033) |
| V | Magnetism | ✅ COMPLETE (034–039) |
| VI | Knit | ✅ COMPLETE (040–042) |
| VII | Integration | 🔄 **2 of 7** — 043 ✅, 044 ✅; 045–049 remain |

⚠ **OPEN — stage043's note and its script disagree, and the SCRIPT is authoritative for counts.** The
script asserts **exactly 152** manifest IDs
(`research/pde_ledger_v2/scripts/ledger_stage043_irreducible_count_range_sympy_audit.py:436,442-443`) and
implements `reduction-debt-counted-once` / `extension-convention-open` as **disjoint peer categories**
(`:165-176` — `22 + 18 + 9 = 49`). The note says **"≈ 152"** and presents both as roll-up **sub-tags**
of the `[40,49]` continuous total
(`research/pde_ledger_v2/notes/stages/ledger_stage043_irreducible_count_range.md:265,272-273,281`).
Unresolved; recorded so the range is not quoted from the note alone.

## STANDING RULES

- ⭐⭐ **PHYSICS, NOT CEREMONY** (user decision, 2026-07-29/30). Two-person toy-physics self-consistency
  project; checks exist to catch a **wrong derivation**. ⛔ **Immutability is not a discipline here** —
  files are freely editable, and nothing is frozen, byte-perfect, or under a custody rule. **The test for
  any check or process rule: does it catch a way the PHYSICS could be wrong? → keep. Only a way the
  TOOLING could be wrong, or a motivated adversary? → cut.** Roles collapse to **one builder; one fresh
  reviewer for prose and process, TWO mutually independent fresh review legs for physics-bearing artifacts**
  (amended 2026-07-30), with the physics leg **blocking**. Owned by `docs/development_pipeline.md`.
- **Findings are the product; green is not the goal.** A result that breaks the concept is welcome and
  first-class. A clean "it all works" is suspicious.
- **The ledger shows the SURVIVING solution only.** Discarded approaches →
  `research/pde_ledger_v2/notes/ledger_exclusions_failures_paper_backlog.md`.
- **Never adjust the process because the corpus is inconvenient.**
- **Only held-out DIMENSIONLESS ratios test the model.** `G`, `c`, `ℏ`, `ℓ_P` are calibration.
- **AI prose never establishes a math fact.** The builder codes and runs dual-engine; verification happens
  on **fresh agents** — **two independent ones where the artifact is physics-bearing**, one for prose and
  process. ⚠ Stop for the user at a decision, a blocking finding or a no-go — **not** at
  every chunk boundary. See `docs/development_pipeline.md`.
- **One canonical doc per workstream.** Fold new findings in — do not write a new doc re-explaining
  ground an existing one covers. Delete dead docs (git preserves them), but extract unique content first.

## MAP — what you want → where it lives

| You want… | Look here |
|---|---|
| ⭐⭐ **The model** — throughline, per-sector derivation atlas, honest earned/calibrated/R1/departure ledger, glossary | `docs/model_map.md` |
| ⭐⭐ **The current front** — the derivation walkthrough: method, per-step record, checks, step order, the closing certification | `docs/derivation_walkthrough_plan.md` |
| ⭐⭐ **The read-first handoff** — where we are, open defects, measured traps | `research/pde_ledger_v2/_scratch/NEXT_SESSION.md` |
| ⭐ **The step records** — one file per derivation step | `research/pde_ledger_v2/walkthrough/` |
| ⭐⭐ **The walkthrough DECISIONS** — user calls that shape what gets counted; ⛔ **can overturn a ruling in the plan, and D-01 already has** | `research/pde_ledger_v2/walkthrough/DECISIONS.md` |
| ⭐ **The two scripts** — the shared import (every quantity, its defining equation, what it derives from) + the dimensions check | `research/pde_ledger_v2/reduction/` |
| ⭐ **The workstream continuing behind it** — the dimension rewrite | `research/pde_ledger_v2/manifests/DIMENSION_REWRITE.md` |
| ⭐ **How we work** — pipeline, roles, the review gauntlet | `docs/development_pipeline.md` |
| The ledger-build resume detail | `research/pde_ledger_v2/notes/RESUME_ROADMAP.md` |
| What's left across ALL sectors → "simulation-ready" | `docs/development_plan.md` |
| Per-Part build history + user-gate records | `research/pde_ledger_v2/notes/part{1..7}_*_atomic_split.md` |
| Per-stage notes | `research/pde_ledger_v2/notes/stages/` |
| Every knob: dimension, class, provenance, reduction debt | `research/pde_ledger_v2/notes/parameter_register.md` |
| The irreducible-count audit | `research/pde_ledger_v2/notes/midway_knob_audit.md` |
| Every value classified DERIVED / INPUT / gap / benchmark | `software/stage1_solver/decisions/14_value_provenance_and_calibration_map.md` |
| The simulation hand-off spec (equations + BC packet + R1→R4 ladder) | `docs/two_throat_simulation_handoff_spec.md` |
| Numbered decisions (esp. **16**, the `P`-retirement) | `software/stage1_solver/decisions/` |
| Retired approaches / the failures-paper backlog | `research/pde_ledger_v2/notes/ledger_exclusions_failures_paper_backlog.md` |
| The calibrate-predict methodology | `software/stage1_solver/decisions/09_calibrate_predict_methodology.md` |
| ⭐ **The EM-track record** — U1/U2 + Phase B/C, the `𝔅` boundary-operator verdict (144/144 UNRESOLVED). **Sole home of that narrative** | `docs/em_analog_next_phase_handoff.md` |
| The EM physical picture + MacCullagh template + leak findings | `software/stage1_solver/decisions/15_em_medium_native_physical_picture.md` |
| Full current state + resume-after-`/compact` pointer (**sync this file with it**) | `software/stage1_solver/decisions/13_emergent_constants_derivation.md` §0 |
| The defect-regime + held-out-surplus roadmap | `docs/defect_interaction_map.md` |
| The gravity moving-throat PDE gate checklist | `research/pde_ledger/notes/stages/moving_throat_pde_completion_ladder.md` |
| The pre-registration (what was committed to in advance) | `docs/pathA_preregistration.md` |

⛔ **Do NOT read `docs/conceptual_foundation.md`** — vision/history, superseded by `docs/model_map.md`. It
predates the EM reconsideration and the retired `P` field, and re-confuses.

## Reference — the `pathA_22b` verdict equation (the *earlier* verdict-count framing)

```
P0 · χ_Q · g_mhat² · λγ⁵ / g_G  =  54/5
 ✓     ✓     gap1     gap2  cal-on-G     (✓ = derived; gap1 g_mhat absorbs 54/5; gap2 λγ ← EM anchor)
G = (a·c_s²/m_GNLS)·g_G ,  m̂0 = (c_s/(a²·√m_GNLS))·g_mhat ,  c = λγ·c_s
```

Here `χ_Q ≈ 0.712` (pathA_22b Gate-3). The moving-throat ladder's Gate 4 (`pathA_33`) later derived
**`χ_Q = 1`** in the outgoing-DtN Hankel context — a different computation of the same-named factor.
**Reconciling the two is a live Part-VII debt**; the ladder reaches the same `54/5` via the
earned/calibrated split `54/5 = 2·27/5`, not via this product.