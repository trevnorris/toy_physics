# S11c-c2 — self-energy fold: closing the slab operator with the curved-bulk response (step record)

Slug `S11c_c2_self_energy_fold`. Step record for **S11c-c2**, the second half of the S11c-c curved-interface bulk
closure (the S11c-c decision-list row was split c1/c2 by user choice 2026-09-03; `directives/S11c_decisions.md`
N2). c2 folds c1's closed permeable face response `(δp_s,J_s,t_s)(V_s,μ_θ)` into the S11c-b variable-coefficient
slab operator `S11CB_SLAB_OPERATOR` — whose θ-row and mechanical rows still carry the face pressure `δp_s`
symbolically — and **re-extracts** the off-diagonal transverse↔`{θ,e_W,u_L}` coupling from the CLOSED full operator,
yielding the coupled **nonlocal self-energy operator**. Physics authority `directives/S11c_c2_SHARED_PHYSICS.md`
(spec v2, `16849fc6`; §5c corrected `30d4b72d`). The profile-conditioned spectrum/scattering + leakage are S11c-d
(`N5`), not this record.

⭐ This record is the **interpretation** layer. The engines PRINT objects and state no conclusions; every result
names the computed object and the commit/record behind it, and every quantitative claim about a run carries its
command in `_measurements/` (rule 2).

> ⭐⭐ **STATUS OF THE STEP (2026-09-09).**
> **PER-ENGINE (SymPy) — the self-energy fold is SOUND.** The fold **wiring** (substitute the closed `δp_s(V_s,μ_θ)`
> + its w-jets into the symbolic `δp_±` slots — ⛔ not a closed `J_s`; operator-inverse response; the `dtn_operator`→
> `dtn_kernel` bridge; `V_s→face_velocity`; computed w-jets; the ε-strip) and the A/C/D1–D6 constructions are
> **two-leg confirmed** (fresh Claude agent + Grok, review-until-clear, + my rule-13 verification; `8f3a017f`), and
> the emitted increment/kernel **VALUES are unaffected** by everything below. ⛔ **NOT "0 defects"** (the `8f3a017f`
> commit subject overstates — the adjudication record was itself corrected). ⛔ **F and G are NOT established
> physics:** their STEP-A adjudication instruments (`verify_F`/`verify_EG`) were **WITHDRAWN**
> (`_measurements/S11c_c2_FG_regrounding_deferred.md:28-30`), so the earlier "F = the genuine coupling decouples" and
> "G = directional/one-way" are **withdrawn interpretations**, ⛔ not adjudicated results (recorded as raw
> observations only). **F/G re-grounding is OWED** — a numeric-probe (Schwartz–Zippel) diagnostic
> (`..._FG_regrounding_deferred.md:31-36`) — but was **PAUSED INDEFINITELY** by user decision after the full-symbolic
> route hit a tractability wall (`cfb2494c`): a **standing, presently non-blocking debt** (⛔ do not resurface it as a
> BLOCKER; ⛔ pausing did NOT discharge it). A separate light §5e/§3c wording clarification is **OWED** (§5e still
> reads "must vanish"; the §3c increment retains the `−extract(open)` open-slot O(ε) piece by construction — a
> structural fact independent of the withdrawn `verify_F`).
> **CROSS-ENGINE (this box) = the N6 representation-invariance thread ONLY.** There is **no WL self-energy engine and
> no self-energy comparator**; the assembled self-energy operator's full cross-engine residual is **DEFERRED**
> (≥64 GB, with c1's four giants, `DEFERRED_HEAVY_RUNS.md`). What WAS done cross-engine is **N6**: a blind Wolfram N6
> engine + the N6 comparator + a leg-cleared Path-B **disposition**. **N6 is covariance (Reading B), dual-engine
> confirmed at retained order**; the surfaced cross-engine OPERAND agreement is carried as an explicit **DEBT** (⛔
> NOT "weak N6"; ⛔ NOT "c2 has everything it needs"). Governing disposition:
> `_measurements/S11c_c2_N6_reconcile_disposition.md` (`0bca95f3`).

## What the step computes
On the inherited S11c-a/S11c-b background (in-plane-varying thickness `W_bg(y)=W̄₀[1+η w₁(ξ)]`; anchorings
`α∈{LAB_HELD,MATERIAL_ADVECTED}`; density representatives `ρ∈{ρ_4D,ρ_br}`; the two faces `s∈{+,−}` the slab EOM
already sums; power counting `(ε,η,σ_W)`, `N12`): **close** the slab operator by substituting the c1 closed face
response into its symbolic `δp_s`/`d_w_δp_s` slots per `(α,ρ)`; **re-extract** the off-diagonal
transverse↔`{θ,e_W,u_L}` coupling from the CLOSED full operator by the S11c-b §3c weak variational restriction —
the **close-then-extract** ordering (extract/eliminate don't commute: close FIRST, then re-extract); and emit the
**nonlocal self-energy** as the **substitution increment** `S11CC2_SELF_ENERGY_INCREMENT = extract(close) −
extract(open-symbolic)` (both operands re-extracted from `SLAB_OPERATOR` with the same extract), the assembled
`S11CC2_CLOSED_SLAB_OPERATOR`, and the `S11CC2_CLOSED_COUPLING_KERNEL`. The `t_s` traction (carrying `Λ_X` + `δp_s`)
routes into the mechanical rows; the closed `δp_s` (+jets) closes the θ-row's already-#90-folded flux terms (⛔ there
is no `J_s` slot — adding one would double-count). c2 also **re-adjudicates the six c1 items the fold makes
load-bearing** (§3d — carried below, ⛔ surfaced not pre-adjudicated), and the **N6 representation invariance** of
the increment.

## The engines + the N6 comparator
- **SymPy self-energy engine** `scripts/S11c_c2_selfenergy_fold_sympy_audit.py` (`S11CC2_*`; astra/`gpt-6-astra`
  build). Reads the inherited model through `ledger_fold.load_model` over the atomic frozen base
  `scripts/S11c_b_exports.py` with the c1 delta `scripts/S11c_c1_exports.py` folded on top (§7). Physics 2-leg
  reviewed `8f3a017f`; publication-only export repair (drop the increment to EMIT-only, keep both closed operators,
  60 MB→21.4 MB) `aa76105a`. Export `scripts/S11c_c2_exports.py` committed `aa76105a`. ⚠ The full 499 MB audit
  `.out` is EPHEMERAL/reproducible (⛔ not committed); the export IS committed.
- **Blind Wolfram N6 engine** `mathematica/S11c_c2_N6_mathematica_audit.wl` (`WL_S11CC2_N6*`) — **blind**, imports
  nothing, re-derives the S11c-a substrate + S11c-b pressure-slot rows + c1 response + constitutive `μ_E`/`μ_M` from
  the sibling specs (`S9_export_chain_rebuild_directive.md:16-18` is the only cross-engine control), and reproduces
  the carrier reconcile `C_E−C_M` + the source-naturality `R_cov`. ⚠ It re-derives the N6 objects, ⛔ NOT the
  assembled self-energy increment. Rebuilt fresh from the scrubbed value-free directive `28f87dec` → certified
  engine `e11f2f82` (the `.wl` engine content sha256 `e5cea55b…`); the compared `.out` is committed at `ae73b884`
  (regenerated from that engine, annex/GIN — `datalad get` the `.out` after a fresh checkout), superseding the stale
  first-clearance `48a0b4e7`.
- **SymPy N6 engines** (the rep-invariance thread, all astra-built, all build-leg-cleared): the diagnostic
  `S11c_c2_N6_diagnostic_sympy.py` (`I_{M→E}`, `R_N6`), the reconcile `S11c_c2_N6_reconcile_sympy.py` (carrier
  bridge `C_E−C_M`, source bridge, 3-way split), and the covariance `S11c_c2_N6_covariance_sympy.py` (the decisive
  source-naturality `R_cov`). Committed `.out`s under `scripts/out/` (annex/GIN).
- **N6 comparator** `scripts/S11c_c2_N6_cross_engine_comparator.py` — the frozen-T7 join: joins by object name,
  pairs residual operands, three-valued, rejects a native boolean, PRINTS and decides nothing (rule 2). Built +
  re-reviewed SOUND `2d12f287`; RUN clean `a094b284` (0 deferred/parse-fail/zero-extract, peak ~330 MB; the ~1.6 GB
  output is ephemeral/reproducible; committed tally `_measurements/S11c_c2_N6_comparator_run_tally.txt`).
  ⛔ There is **no** cross-engine comparator for the self-energy operator itself.

## The arc (each result: commit/record + how verified)
- **Self-energy fold — the WIRING + A/C/D + increment VALUES are per-engine SOUND** (`8f3a017f`; adjudication
  `_measurements/S11c_c2_physics_review_adjudication.md`). Both legs (fresh Claude + Grok) agreed the fold wiring +
  A/C/D1–D6 with shown CAS. The legs SPLIT on B/E/F/G; I first resolved those with orchestrator-authored scripts —
  ⛔ **those F/G conclusions were later WITHDRAWN** (the corrected process: orchestrator never authors the CAS
  instrument, CLAUDE.md `6f8dbd34`): **F** (the raw observation was `.doit()` → integrand literally 0, "the
  closure-induced coupling decouples") and **G** (the raw observation was a directional/one-way increment, reverse
  block identically zero, both blocks emitted, no adjointness residual per §3b) are **withdrawn interpretations**,
  ⛔ NOT adjudicated physics; **B** was F's `−extract(open)` open-slot residue. ⚠ My first "0 defects / F,G,E all
  false-positive" verdict was an **over-reach** (Codex-sol compact-prep verify, rule 13). F/G were then re-grounded
  the corrected way (question-vet + a 4-round reviewed build directive) but hit a tractability wall → numeric
  re-grounding OWED, PAUSED (below). **E/N6** was itself superseded (§5c mis-spec, below). ⇒ what STANDS from this
  build is the wiring + A/C/D + the emitted increment/kernel VALUES.
- **§5c MIS-SPECIFIED N6, then CORRECTED** (`30d4b72d`). The self-energy engine's own `REP_INVARIANCE_RESIDUAL`
  compared the two ANCHORINGS (distinct physics per S11c-a §2c) — the WRONG object (a nonzero value is EXPECTED,
  ⛔ not a defect). The real N6 (parent S11c-a §5a / sibling c1) is **Eulerian-vs-material-coordinate within a FIXED
  anchoring**. §5c corrected (2-leg spec review-until-clear); c2's real N6 was then built in the separate N6 engines.
- **N6 per-engine (SymPy) RESOLVED — covariance (Reading B), user-adopted** (`d21c8ff5`;
  `_measurements/S11c_c2_N6_RESOLVED.md`). Evidence chain: (1) geometry reconciles — carrier `C_E=C_M` no-nonzero all
  4 cases (live control); (2) the residual is purely constitutive — `R_N6 = I_E − I_{M→E}` nonzero in ~18 columns,
  3 of 4 cases (**R_N6 = 18/288**), localizes ENTIRELY to the source channel (SPLIT_CHECK=0); (3) the decisive test —
  `R_cov = ms − source_terms(μ_E.subs(Φ), V_E)` shows **no nonzero found, all 4 cases** (conditional δ≈2.6e-22),
  knives bite ⇒ the material builder faithfully implements the declared map Φ; the nonzero `R_N6` is (to that bound)
  the Φ-image of the source. ⭐ **N6 PASSES as COMPLETE OPERATOR COVARIANCE** at retained order under the declared-map
  premise + PIT qualification — ⛔ NOT strict `R_N6=0` (the increment is a response-operator coupling block, expected
  to transform covariantly like a vector's components under a change of frame; equality in genuinely COMMON variables
  is mandatory and IS met). ⚠ My "localization ⇒ satisfied" was an OVER-CLEAR caught by the adj-review gate (Codex-sol
  RIGHT); closed via the `R_cov` sufficient test.
- **Blind Wolfram N6 — per-engine CLEAR** (`48a0b4e7` → rebuilt `e11f2f82`/`ae73b884`). All 5 ablation controls BITE
  one-sided; δ WL-derived (own primes `{1e9+9, 998244353, 1004535809}`); blindness = import-free + byte-identical
  isolated-vs-in-repo. ⚠ 2 residual expected-zero soft-leaks in the directive were scrubbed post-hoc, and the engine
  was REBUILT fresh from the scrubbed value-free directive (user, maximal rigor); genuineness is independently
  established by the biting FORM ablations (a forced zero cannot move under corruption; these do).
- **N6 comparator RUN + Path-B disposition** (`a094b284` → disposition `0bca95f3`). See *Established vs owed* and
  *Method notes*.

## Established (per-engine / cross-engine) vs owed (surfaced/deferred/paused)
- **ESTABLISHED — per-engine SymPy SOUND (2-leg):** the self-energy fold wiring + A/C/D1–D6 and the emitted
  closed-slab operator + closed coupling kernel + substitution-increment **VALUES**. ⛔ **NOT F or G** — those
  interpretations rest on the WITHDRAWN `verify_F`/`verify_EG` instruments (numeric re-grounding OWED/paused, below);
  they are recorded as raw observations, ⛔ not established physics. **PER-ENGINE N6 = operator covariance
  (Reading B)**, in BOTH engines (SymPy `R_cov` no-nonzero; blind WL reproduces carrier + `R_cov`).
- **ESTABLISHED — cross-engine (N6 thread, dual-engine, on this box):** every MATCHED comparison in the
  vanishing/control/premise subset agreed — `R_cov`/`R_cov_baseline`/`R_cov_control_delta` (160/0 each), the carrier
  bridge `CARRIER_BRIDGE_RESIDUAL` (320/0), `SOURCE_CONTROL_DELTA` (160/0), + genuine nonzero-operand/support
  agreements (4 `FROZEN_RELATIONS` premise leaves; 400 structural support). ⚠ The covariance-channel matched zeros
  are `(0)−(0)` — a **dual-engine confirmation of the VANISHING statement (Reading B), ⛔ NOT operand AGREE.**
- **OWED — the cross-engine operand DEBT (this step's central open item):** the surfaced cross-engine OPERAND
  residuals are UNADJUDICATED — the blind-WL-vs-imported CARRIER (40), the constitutive SOURCE (76) + Φ (18). The
  collapse-instrument that would test whether they reconcile representationally was **CLOSED at v3 NOT-SOUND** (the
  post-EL graded coefficient table is the WRONG OBJECT — `EL(T·L)≠T(EL·L)` for the position-dependent `R_W`, and the
  ungraded `R_W` can't be per-grade rewritten under the no-grade-mixing contract; 2-expert consult, both Path B).
  ⛔ **Do NOT let "representational-difference-UNADJUDICATED" become "known to be just thickness"** — the leftover
  SHAPE was NOT inspected; the alternatives include a genuine constitutive-convention mismatch OR an implementation
  error. `R_N6`, the channels, the reconcile-engine sources, and the guards were **never directly compared**
  (schema non-join, ⛔ not heaviness). Governing: `_measurements/S11c_c2_N6_reconcile_disposition.md`.
- **DEFERRED (≥64 GB, `DEFERRED_HEAVY_RUNS.md`):** the assembled self-energy operator's full cross-engine residual +
  c1's four giant families. c2 was constructible + N6-cross-engine-testable on this box for its own increment (§7);
  the full residual is the ≥64 GB work.
- **OWED but PAUSED INDEFINITELY (non-blocking):** the F/G numeric-probe re-grounding
  (`_measurements/S11c_c2_FG_regrounding_deferred.md:31-42`) — a standing debt the user paused 2026-09-06 (the
  full-symbolic route hit a tractability wall); ⛔ do not resurface it as a BLOCKER, ⛔ but pausing did NOT discharge
  it. The increment VALUES are unaffected; the withdrawn `verify_F`/`verify_EG` F/G conclusions do NOT stand.

## Method notes
- ⭐ **`I_{M→E}` is the NATIVE MATERIAL-COORDINATE-ROUTE increment at a FIXED anchoring `α`, ⛔ NOT an anchoring and
  ⛔ NOT a "mapped-to-Eulerian operand."** `I_{M→E}^{α,ρ} = extract(close(SLAB_M) − SLAB_M)`
  (`S11c_c2_N6_diagnostic_sympy.py:806-853`; SHARED_PHYSICS §5c:305-320) — N6's two routes are Eulerian vs
  material-COORDINATE at the SAME `(α,ρ)` (the routes are the **representation** axis, ⛔ **never** the anchoring axis
  `{LAB_HELD, MATERIAL_ADVECTED}` — conflating them is exactly the §5c mis-spec that `30d4b72d` corrected). The two
  increments are differenced directly with **no separate `T`/final pullback on the increment**, because the material
  builder ALREADY performs the internal covector conversion (inverse-transpose) into the **common Eulerian face
  basis**. The "`M→E` mapped-operand" label misleads by implying `I_{M→E}` is the full frame-transformed (Φ) image —
  it is not; the FULL frame-change faithfulness is checked SEPARATELY by `R_cov` (source-naturality). ⇒ preserve BOTH
  findings: `R_N6 = I_E − I_{M→E}` **nonzero (18/288)** AND `R_cov` **no-nonzero** — consistent under Reading B, not
  contradictory. [[feedback_reconcile_representational_bridge]]
- ⭐ **A nonzero cross-engine residual is not a disagreement** (the comparator prints raw, rule 2). The N6 disposition
  is Path B: the surfaced operand residuals are carried as a DEBT because the only sound reconcile instrument is
  UPSTREAM of the EL differentiation that produced the source operands (a new emit path in both engines, or a
  derivative-aware chain-rule bridge — either ≥ the doomed instrument's cost; a replay validates the replay, not the
  already-emitted `.out`). [[feedback_reconcile_bridge_must_precede_el]] [[feedback_handcode_comparison_never_blanket_collapse]]
- ⚠ **The 3 N6 premise caveats (Reading B does NOT close them — PREMISE checks, ⛔ not defects):** (1) is Φ itself
  physically correct — derive it from the actual material motion + density/measure transformation, not merely confirm
  the builder reproduces the declared Φ (`R_cov` cannot exclude an error SHARED by the declared premise AND both
  routes); (2) does the face velocity `V` transform correctly — `V_E≡V_M` (SHA-equal) is BUILDER agreement, and the
  prediction uses `V_E`, not `Φ(V_E)`; (3) extracted-block omitted-block leakage. [[feedback_never_freeze_a_varying_field]]
- ⚠ **A §5e/§3c wording clarification is OWED** (⛔ not applied): SHARED_PHYSICS §5e still says the uniform-limit
  increment "must vanish", but the §3c increment retains the `−extract(open)` open-slot O(ε) piece **by
  construction** (a structural fact about the increment's definition, independent of the withdrawn `verify_F`). ⇒
  the "must vanish" wording is imprecise and owes a clarification (review-until-clear); ⛔ the settled replacement is
  NOT yet established (whether "the genuine closure-induced coupling decouples" holds is exactly the F question the
  OWED numeric-probe re-grounding must decide — do not assert it here). After the spec edit, the export's
  `BUILD_INPUT_DIGESTS` (which pins the spec) goes stale ⇒ lawfully repin + reverify.
- ⛔ **Rule 17 (background density) — a SURFACED freeze re-adjudicated in c2, ⛔ not waved through** (c1 seal 5). c2's
  fold sums over the face where `ρ(x)`'s variation is load-bearing; §3d.1 binds `rho_br_bg_rho4_constant` to the
  `background_density_map` before the fold (the O(εη) channel `d(μ_s)/dη|₀ = −μ_θ w₁/ρ_br`), so the freeze cannot be
  emitted as a bare constant. [[feedback_basis_independence_must_not_freeze_spurion]]
- ⛔ **Serialize CAS jobs; watch RSS.** N6 measured LIGHT (comparator peak ~330 MB); the self-energy `.out` is 499 MB;
  the full cross-engine residual is the ≥64 GB work. `run_in_background` is REAPED on this box — detached launch
  (`setsid` + DONE-marker + Monitor). [[feedback_background_tasks_can_die_spuriously]]

## Carry-forward — surfaced for S11c-d, ⛔ NOT pre-adjudicated
These change what may be CLAIMED, ⛔ not any computed object. c2's engines/exports STAND.
- **The cross-engine operand DEBT** + the un-inspected leftover SHAPE (the N6 disposition).
- **The 3 N6 premise caveats** (Φ physical-correctness; `V` transform; extracted-block leakage).
- **The 2 S11c-b sign conventions** that multiply the substituted slots and do NOT cancel from c2's residual
  (`steps/S11c_b_variable_coefficient_operator.md:112-115`): the **face generalized-force** convention (PY `+diff`
  vs WL `−linearVirtualVariation`) and the **#90 closure-fold** sign. (The **kinetic** `−K`/`+K` convention is a
  bulk term independent of the response slots.) The §3d.4 mechanical-power pairing adjudicates the face-force sign;
  cross-engine surfacing is owed (no self-energy comparator exists).
- **The 6 §3d re-adjudications** (census below): (1) background density field-vs-field (rule 17); (2) `t_s` traction
  representation (native covector `t_s = −(δp_s+Λ_X𝒜_s)n̂_s`; WL 4-vec vs PY scalar UNDECIDED from c1); (3) DtN
  kernel-vs-whole-form (the AGREE'd two-momentum kernel is load-bearing; the raw `dtn_operator` whole-form is c1
  UNDECIDED); (4) traction-vs-slab mechanical-power pairing (adjudicates the face-force sign; a one-sided `t_s`-sign
  corruption must move the residual against the slab kinetic term); (5) flat-resolvent leg-labeling (PY output-leg vs
  WL input-leg, equal on `k=k′`); (6) `μ_R,bg` form control (c1-reserved for c2). Each is emitted per-engine SymPy;
  the c1-UNDECIDED imports (i)–(v) among them stay cross-engine-UNDECIDED.
- **c1 ENERGY** (PY closed-form vs WL far-field integral) — UNDECIDED, deferred with c1's giants.
- **F/G interpretations** — a standing OWED debt (the numeric-probe re-grounding), **PAUSED INDEFINITELY** and
  presently non-blocking (⛔ do not resurface it as a BLOCKER; ⛔ but pausing did not discharge it — the `verify_F`/
  `verify_EG` conclusions are withdrawn, the increment VALUES stand).

## Census (§3d control/premise leaves) — a FACT-LOOKUP, ⛔ no instrument
The 8 nonzero control/premise leaves the N6 comparator surfaced (`N6COV_ACTUAL_CONTROL_PARAMETERS` 4,
`N6COV_PHI_DOMAIN_CENSUS` 4 [+44 BOOL rejected]) are production-control equivalence + domain-coverage equivalence,
⛔ not "do they carry physics." The frozen 6-row crosswalk: `kappa_a↔ADVECTION`, `kappa_j↔JUNK`,
`imported_theta_e_jets↔DOMAIN`, `coverage↔COVERAGE`, `uncovered↔UNCOVERED`, `max_present_jet_rank↔MAX_RANK`; WL-only
`THICKNESS` / `MATERIAL_NORMAL` stay one-sided (PY `actual_amplitudes` retags only theta-advection + junk); the
`n6cov_J_mu`↔`junkMu` spelling. If a census leaf reflects differing μ-jet coverage it routes to the source channel,
⛔ not to bookkeeping. Full disposition: `_measurements/S11c_c2_N6_reconcile_disposition.md` §8.

## What's next
**S11c-d** (profile-conditioned spectrum/scattering + leakage, `N5`): consumes c2's `CLOSED_SLAB_OPERATOR` (the
closure-modified diagonal → spectrum/resolvent poles) + `CLOSED_COUPLING_KERNEL` (off-diagonal → Born/mixing +
gradient-driven leakage). ⚠ Because S11c-d uses gradient-driven mixing/leakage, the carried cross-engine operand
DEBT is **material to a consumer**, ⛔ not dismissible on the strength of covariance alone. ⚠ **NO per-substep card**
— N1 specifies ONE S11c roll-up card, produced only after S11c-e. The full cross-engine self-energy residual + c1's
four giants (≥64 GB) and the c1-UNDECIDED imports close there / at the S11c-e surviving-chain audit (trigger =
dependency, ⛔ not blanket S9→S11). Until then S11c-d must treat the self-energy operand agreement, the background
density, `t_s`, the `dtn_operator` whole-form, and ENERGY as NOT cross-engine-closed.
