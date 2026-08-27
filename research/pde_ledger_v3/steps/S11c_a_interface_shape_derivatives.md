# S11c-a — tilted-face interface shape derivatives (step record)

> ⛔⛔ **DRAFT — NOT FINAL, and it OVERCLAIMS as written (corrected 2026-08-27).** The body below concludes a
> clean "the two engines AGREE; every cross-engine residual is representational." That rested on an
> aggregation-side **blanket collapse** (`X(args)→X`) that DELETES real dependence — it **hid two findings**:
> **(1)** PY froze the **in-plane** bulk current in w — a real §1b defect, now FIXED + committed `8c1a5ed1`
> (a 4th fix; the current now matches WL's w-structure, residual = coordinate convention + δρ identity);
> **(2)** FACE_SHIFT (T-e) has `join=0` (WL missing the DENSITY axis) — the family is **UNCOMPARED**, not
> "bookkeeping" as §"Bookkeeping families" claims — STILL OPEN (task #69). Also to correct on rewrite: the
> control-battery §("every object bites, no dead control") is FALSE (the ∂W_bg form ablation bites only the
> geometric objects; FACE_SHIFT/PROJECTION/EVOLUTION/VIRTUAL_CONSTRAINT don't); the δρ mechanism is
> family-dependent (PY carries u explicitly in dynamic/shape, and the EVOLUTION "closure" is void for RHOBR);
> and the standing limit must read `v_bulk_normal_0`, not `v₀`. ⛔ Do not treat this as the DONE gate. Rewrite
> after FACE_SHIFT + a full comparator re-run vs the `8c1a5ed1` transcript. Method fix:
> hand-code the comparison, never blanket-collapse ([[feedback-handcode-comparison-never-blanket-collapse]]).

Slug `S11c_a_interface_shape_derivatives`. This is the step record for **S11c-a**, the first sub-step of the
S11c staged family (the non-uniform transverse coupling). Its scope, ratified in the S11c decision list
`directives/S11c_decisions.md` (N1/table row S11c-a), is **background & geometry**: the full tilted-face
shape-derivative of every interface law — outward normal, true-area measure, face velocity, relative flux,
shifted-face bulk trace, projected mass balance, traction, kinematic balance, and virtual work — on a
**non-uniform** background. It exports the **geometric boundary operators only**; the variable-coefficient
slab operator and its spectrum are S11c-b..e (spec §1b; decision list table). The physics authority is the
sub-step spec `directives/S11c_a_SHARED_PHYSICS.md`.

⭐ This record is the **interpretation** layer. The engines PRINT objects and state no conclusions; every
result below names the computed object and the commit that produced it, and every quantitative claim about
the run carries the command that produced it in `directives/_measurements/` (rule 2).

## What the step computes

Given the S11c background ansatz — in-plane-varying thickness `W_bg(y)≡W̄₀[1+η w₁(ξ)]` and response modulus
`μ_R,bg(y)`, two density representatives (`ρ_4D`-constant / `ρ_br`-constant), and two physical anchorings
(LAB_HELD `Q_bg(x)` and MATERIAL_ADVECTED `Q_bg(χ(x,t))`) — S11c-a computes the **first shape derivative**
of the tilted-face interface objects T-a..T-i (spec §4), multigraded by `(ε,η,σ_W)`: the wave amplitude `ε`,
the zero-jet contrast `η`, and the supplied first-jet bookkeeper `σ_W≡η W̄₀/L_W`. These are the geometric
substrate on which S11c-b builds the variable-coefficient coupling kernel. S11c-a computes no curved-bulk
response (spec §1b).

## The two engines and their agreement

- **SymPy engine** `scripts/S11c_a_interface_geometry_sympy_audit.py` — builds the shape derivatives from
  the spec's ansatz, face maps, and supplied laws; emits 47 `PY_S11CA_*` tag families. Committed through its
  two §3c-class fixes: shifted-trace `c36beac4`, projection current-freezing `49b5c525`. Its committed
  transcript `scripts/out/S11c_a_interface_geometry_sympy_audit.out` was regenerated from HEAD and committed
  `afdc8158` (the provenance note below).
- **Wolfram engine** `mathematica/S11c_a_interface_geometry_mathematica_audit.wl` — **blind**: imports
  nothing, re-derives every object from the spec. Committed through its fix: background-current `6fae82b8`.
  Its blindness is the only cross-engine control — it cannot transcribe the `.py` because it never reads it,
  so an agreement between the two is independent construction, not a copy. Transcript
  `mathematica/out/S11c_a_interface_geometry_mathematica_audit.out` (`6fae82b8`).
- **T7 cross-engine comparator** `scripts/S11c_a_cross_engine_comparator.py` (+ `test_S11c_a_cross_engine_comparator.py`,
  14 synthetic tests) — joins the two transcripts by emitted-object name with axis-typed keys, prints
  `operand_A` (PY) / `operand_B` (WL) / `A_minus_B` per case and per-family accounting `{join, py_only,
  wl_only, duplicate_key, parse_failed, axis_set_mismatch}`; it asserts nothing on measured payloads and
  emits no PASS/FAIL/AGREE verdict (rule 2). Built via a delegated build brief after three prose-directive
  rounds failed independent review (the rule-15 note below); cleared by two build legs (fresh Claude + Grok)
  with FORM-ablation teeth. Committed `50f43123`. Run over the two committed transcripts: 6714 operand
  triples, exit 0, `~/.s11_build/comparator_run.out`.

**The result: the two independently-built engines AGREE on T7.** Across all 39 tag families, after the three
§3c-class engine defects below were fixed, every nonzero cross-engine residual reduces to a **representational
identity** — a difference in how the two engines *name or expand* the same object — or the known CONORMAL
§3c form. No genuine T7 physics disagreement survives.

## The three §3c-class engine defects — each engine's own error, caught by the two-engine method

Each was a single-engine implementation error that the sibling engine did not share; the disagreement was the
measurement (rule 1). None was a spec ambiguity.

1. **Shifted-trace (PY)** `c36beac4`. PY mishandled the §3c trace-linearisation `δ[f(x,h_s)] = δf(x,h_s⁰) +
   δh_s ∂_w f⁰`; WL was right. Two SymPy bugs fixed.
2. **Background current (WL)** `6fae82b8`. WL carried a §3c-forbidden free-premise background current; the
   fix builds `j⁰=ρ⁰v⁰` from the supplied background rather than asserting it, and the survivor is zero under
   static continuity.
3. **Projection current-freezing (PY)** `49b5c525`. PY froze the perturbation current at the face, dropping
   the normal jet `∂_wδj_w` the §1b projection requires; WL retained the full w-dependent bulk current. A
   from-spec adjudication (Codex + Grok, unanimous) ruled WL faithful, PY defective; PY was fixed to carry
   `δj_w`'s w-dependence. Found only because a "field-name map" trap was flagged — the honest arg-preserving
   comparison surfaced what an argument-strip had hidden.

## Cross-engine comparison — the adjudication (rule 13)

The comparator prints raw residuals; classification happens on our side, and a naive string zero-test gives
false nonzeros (zero matrices/tuples, integral-linearity, applied-function-vs-symbol). Under the correct
canonicalizer — structural zero-test, collapse of an applied bulk quantity `X(args)→X`, and integral
linearity — the residuals fall into these classes (`directives/_measurements/S11c_a_adjudication_*`):

**Clean zero (both engines identical).** FACE_NORMAL, FACE_MEASURE_SHAPE_DERIV, FACE_VELOCITY, RELATIVE_FLUX,
KINEMATIC_BALANCE, EVOLUTION_MASS_BALANCE, EVOLUTION_TERM_ORIGINS, VIRTUAL_CONSTRAINT.

**Representational identity — the closure coefficient μ_θ** (TRACTION, CLOSURE_SHAPE_DERIV,
VIRTUAL_WORK_SHAPE_DERIV). WL writes `μ_θ^α` as an applied function `mu_theta_{L,M}(x1,x2,x3,time)`; PY as a
bare symbol. The whole residual is `Σ coeff·(mu_theta − mu_theta(x1,x2,x3,time))`, which collapses to 0 for
all 48 cases. Decisive probe: across the whole comparator run (`~/.s11_build/comparator_run.out`, 37 300
lines) there are **zero** jet-suffixed `mu_theta` and **zero** `Derivative(mu_theta)` — WL never
differentiates μ_θ. The spec makes the evaluation-point
arguments inert: §2a supplies in-plane derivative maps for `W_bg` and `μ_R,bg` but **not** for `μ_θ`; §3c
forbids introducing a background derivative as a free premise; §4 T-i keeps `μ_θ` a **named** operand; μ_θ
enters only algebraically through `μ_s^α=μ_θ^α/ρ_br,bg^{0,α}`. Bare ≡ applied. Benign.

**Representational identity — the projected density perturbation δρ_4D** (PROJECTION_*). After collapse and
integral linearity the in-plane current divergence `∂_iδj_i` cancels between the engines; the remainder is a
density-time term where PY carries a single named symbol `delta_rho_4D_bulk_t` and WL carries its explicit
§3a expansion. Binding PY's symbol to the spec-built `∂_tδρ_4D = ∂_t[ρ_4D,bg^{0,α}·θ]` sends the residual to
0 for all four (anchoring × representative) combinations, **including** the MATERIAL_ADVECTED × RHOBR_CONSTANT
combination, where the advected background `ρ_4D,bg⁰(χ(x,t))` contributes a first-order advection term
`∝ σ_W·u_{i,t}·∂_ξw₁`. `delta_rho_4D_bulk` is declared a **primitive** symbol in the PY engine (never
expanded) and carried unexpanded in the projection; WL, blind, expands it. Benign.

The whole μ_θ/δρ class was gated the #3 way: a from-spec CAS construction consult (Codex + Grok, each
building μ_θ and δρ from the spec in runnable SymPy and reporting the IFF condition) plus from-spec
adjudication legs (a fresh Claude agent + Grok, reading the spec on whether μ_θ varies and what δρ means).
**Five of five — the two consult engines, the two adjudication legs, and the orchestrator's own rule-13
verification — returned IDENTITY / benign for both classes, no engine blamed.** Two catches worth recording:

- Both Grok legs modelled the MATERIAL_ADVECTED advection as "carried separately by *both* engines and
  cancelling." Direct measurement of the committed operands disproves this — PY's projection operand carries
  **no** `u` terms; it folds the advection into the named symbol, which WL expands. Grok's verdict was right,
  its structural account measurably wrong; Codex's and the agent's account (PY symbol = full object) is the
  correct one.
- The adjudication agent flagged the one contingency it was structurally blind to: the MATERIAL_ADVECTED ×
  RHOBR_CONSTANT verdict flips to a *finding* if PY defined `delta_rho_4D_bulk` as the densification about the
  *local advected* background (`ρ_bg^{0,M}·θ`, dropping the advection) rather than the perturbation of the
  full field. **This is closed by measurement, not interpretation:** PY's EVOLUTION family builds the
  density-time *concretely* — `sigma_exact` plus an explicit `BACKGROUND_ADVECTION = ε·u_t·∇σ0` term (engine
  L1298) — and EVOLUTION_MASS_BALANCE + EVOLUTION_TERM_ORIGINS agree with WL **8/8 to exactly 0**, with PY's
  MATERIAL_ADVECTED operand carrying the `u_{i,t}` advection. PY's density object is therefore the full §3a
  perturbation; the projection's named symbol is a faithful unexpanded name for it.

**CONORMAL_DERIV** is the one remaining non-representational remainder, and it is the already-adjudicated
§3c verdict-A form: a probe-representation / Taylor-order difference (`W_0²`, third order) in the unmapped
`conormalBackground` jets, same physics, no finding (prior adjudication `directives/_measurements/S11c_a_T7_adjudication_verdicts.md`).

**Bookkeeping families** (DIMENSIONS, BACKGROUND_STATE, FACE_MAP, ADMISSIBILITY, FACE_SHIFT) surface as
`Mismatch(STRUCTURE_DISAGREE)` (PY tuple/scalar vs WL Association) or a documented `axis_set_mismatch`
(FACE_SHIFT WL-missing-DENSITY, ADMISSIBILITY PY-missing-BRANCH) — representation/coverage, surfaced honestly
by the accounting, not manufactured.

## The control battery holds and bites

The spec §5 controls were characterized by the same collapse + linearity pass
(`directives/_measurements/S11c_a_control_battery_result.md`):

- **HOLD controls — 0 genuine nonzero cross-engine.** REP_INVARIANCE_RESIDUAL (N4/N6 Eulerian≡material
  route-invariance): 48 zero, 16 structural, 24 unjoined, **0** nonzero. UNIFORM_LIMIT_RESIDUAL (the
  S11c-a→S11b regression): 148 zero, 24 structural, 240 unjoined, **0** nonzero. The route-invariance and
  the regression hold across engines wherever the two join. (The 240 unjoined uniform-limit cases are a
  cross-engine *pairing* asymmetry — a coverage gap the run accounting already flags with
  `families_with_unpaired=15` — not a physics residual.)
- **BITE controls — every object alive.** CONTROL_FORM_RESIDUAL (the C-1 one-direction source form ablation)
  bites for all 18 objects (≥12 each); CONTROL_INDEPENDENCE_RESIDUAL (one-sided route corruption) bites for
  all 6 objects. The `null` cases are legitimately null — a direction/DOF whose `W_bg` jet does not couple to
  that object — not a dead control; each object bites in at least one case.

## The reconciliation schema — pinned for S11c-b..e

S11c-a fixes the cross-engine reconciliation schema that later sub-steps inherit (decision list N8: each
sub-step gets its own comparator and exports, but the reconciliation folds are the family's):

1. **Axis-typed join keys.** Cases join on the object name plus an axis-typed key (BRANCH, DENSITY, FACE,
   DOF, DIRECTION, OBJECT), not a full tuple — PY's bare `DELTA_W` and WL's `DOF_DELTA_W` are the same axis.
2. **Name folds are head renames that PRESERVE arity.** An applied bulk quantity is mapped to a bare symbol
   only for the classification zero-test, never in the comparator's operands; ⛔ never map an `AppliedUndef`
   to a bare `Symbol` in construction — that argument-strip is exactly what hid the current-freezing defect.
3. **Integral linearity is the canonicalizer for projection residuals.** Pull binder-free factors out of the
   bound integral and combine same-`(lower,upper)` integrands before any zero-test.
4. **The closed representational folds:** the closure coefficient `μ_θ^α` (applied-at-face ↔ bare, inert
   args); the projected density perturbation `δρ_4D` (named symbol ↔ §3a expansion `ρ_4D,bg^{0,α}·θ`,
   including the MATERIAL_ADVECTED advection); the CONORMAL §3c verdict-A form. These identities are reviewed
   REGISTRY entries, not comparator name-map surgery.

## Standing limits, freezes, and owed items

- **The background-flow correction `O(v₀|q_n|/ω)`** is inherited as a standing limit from S11b, uncarried and
  unbounded where `|q v₀/ω|≳1` (S11c decision list §E carry-in). S11c-a computes geometry only and does not
  address it.
- **Falsification is deferred to S11c-e.** S11c-a is the geometric substrate; the flux-normalized
  dimensionless conversion FORM and its bench-optics bound are the S11c-e deliverable (decision list table).
- **The family roll-up card** closes S11c after all sub-steps (decision list N1); S11c-a owes no separate
  paper card, only its exports and this record.
- **Provenance (rule 2 bound the orchestrator here):** the comparator's PY input is a regeneration of the
  engine transcript from HEAD, committed `afdc8158` — the uncommitted scratch run it replaced was *not* the
  committed engine's output, and a fresh regen equals it only up to `srepr` term ordering (the engine math is
  deterministic, its serialization is not).

## Operational note

The comparator runs ~37 minutes; heavy PY/WL runs were detached (`setsid`) against spurious background-job
kills, with a watcher loop. Review and consult legs were run with one Grok instance at a time (the single
per-user lock) and Codex at `xhigh` with `--sandbox danger-full-access`. The comparator's build followed the
rule-15 pivot: three prose build-directive rounds each drew fresh independent-leg defects (39 bespoke
per-family schemas cannot be pre-enumerated in prose), so the build was delegated via a brief that locked the
reviewed physics folds and let the two build legs gate the working instrument.

## Provenance

- Spec `directives/S11c_a_SHARED_PHYSICS.md`; decision list `directives/S11c_decisions.md`; scope
  `steps/S11c_SCOPE.md`.
- PY engine `scripts/S11c_a_interface_geometry_sympy_audit.py` (`c36beac4`, `49b5c525`); transcript
  `scripts/out/S11c_a_interface_geometry_sympy_audit.out` (`afdc8158`).
- WL engine + transcript `mathematica/{,out/}S11c_a_interface_geometry_mathematica_audit.{wl,out}` (`6fae82b8`).
- Comparator `scripts/S11c_a_cross_engine_comparator.py` + tests (`50f43123`); build brief
  `directives/S11c_a_comparator_build_directive.md` (+ `_measurements/` twin); run `~/.s11_build/comparator_run.out`.
- Adjudication + control measurements `directives/_measurements/S11c_a_{adjudication_rule13_verify,
  adjudication_matadv_rhobr_verify,control_battery_result,control_hold_classify}.{py,stdout,md}`; leg prompts
  `directives/_legs/S11c_a_repr_identity_{cas_consult,adjudication}.md`.
