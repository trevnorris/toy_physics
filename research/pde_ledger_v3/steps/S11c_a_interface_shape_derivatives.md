# S11c-a — tilted-face interface shape derivatives (step record)

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

> ⚠ **Method correction that governs this whole record (2026-08-27).** An earlier draft concluded a clean
> "the two engines AGREE; every cross-engine residual is representational." That conclusion rested on an
> aggregation-side **blanket collapse** (map every applied function `X(args)→X`, then zero-test), which
> DELETES real dependence and can only ever HIDE a difference. It hid two genuine findings (the in-plane
> current freeze `8c1a5ed1` and the FACE_SHIFT free-premise density background `cccb4f9e`). The classifier
> was replaced by a **hand-coded** comparison — a basic mechanical pass for exact zeros and hand-verified
> renames, an explicit map for the one verified-inert collapse, everything else compared as-is and flagged.
> Every physics claim below is grounded in that hand-coded comparison and the honest per-family accounting,
> **not** the blanket collapse. ([[feedback-handcode-comparison-never-blanket-collapse]])

## What the step computes

Given the S11c background ansatz — in-plane-varying thickness `W_bg(y)≡W̄₀[1+η w₁(ξ)]` and response modulus
`μ_R,bg(y)`, two density representatives (`ρ_4D`-constant / `ρ_br`-constant), and two physical anchorings
(LAB_HELD `Q_bg(x)` and MATERIAL_ADVECTED `Q_bg(χ(x,t))`) — S11c-a computes the **first shape derivative**
of the tilted-face interface objects T-a..T-i (spec §4), multigraded by `(ε,η,σ_W)`: the wave amplitude `ε`,
the zero-jet contrast `η`, and the supplied first-jet bookkeeper `σ_W≡η W̄₀/L_W`. These are the geometric
substrate on which S11c-b builds the variable-coefficient coupling kernel. S11c-a computes no curved-bulk
response (spec §1b).

## The two blind engines and how they were reconciled

- **SymPy engine** `scripts/S11c_a_interface_geometry_sympy_audit.py` — builds the shape derivatives from
  the spec's ansatz, face maps, and supplied laws; emits the `PY_S11CA_*` tag families. Committed through
  three fixes: shifted-trace `c36beac4`, projection current-freezing #3 `49b5c525`, in-plane current #4
  `8c1a5ed1`. Its committed transcript `scripts/out/S11c_a_interface_geometry_sympy_audit.out` (`8c1a5ed1`).
- **Wolfram engine** `mathematica/S11c_a_interface_geometry_mathematica_audit.wl` — **blind**: imports
  nothing, re-derives every object from the spec. Committed through two fixes: background-current `6fae82b8`,
  density-grounding `cccb4f9e`. Blindness is the only cross-engine control — it cannot transcribe the `.py`
  because it never reads it, so an agreement is independent construction, not a copy. Transcript
  `mathematica/out/S11c_a_interface_geometry_mathematica_audit.out` (`cccb4f9e`).
- **T7 cross-engine comparator** `scripts/S11c_a_cross_engine_comparator.py` (+ synthetic tests,
  `50f43123`) — joins the two transcripts by emitted-object name with axis-typed keys, prints
  `operand_A`(PY)/`operand_B`(WL)/`A_minus_B` per case and per-family accounting `{join, py_only, wl_only,
  duplicate_key, parse_failed, axis_set_mismatch}`; it asserts nothing on measured payloads and emits no
  PASS/FAIL/AGREE verdict (rule 2). Built via a delegated build brief after three prose-directive rounds
  failed independent review (39 bespoke per-family schemas cannot be pre-enumerated in prose — the rule-15
  note below); cleared by two build legs (fresh Claude + Grok) with FORM-ablation teeth. The **hand-coded
  reconciliation layer** on top of it is `scripts/S11c_a_handcoded_comparison.py`.

**The measured result.** Independent construction let the two engines DISAGREE (rule 1), and over the course
of S11c-a the two-engine method surfaced **five single-engine defects** (§ below), each caught because the
sibling did not share it. After all five fixes, the final full comparator run vs the committed
both-engines-fixed state (`~/.s11_build/comparator_final_cccb4f9e.out`;
`RUN_ACCOUNTING families=39 families_with_join=38 families_with_unpaired=11 parse_failed=0 duplicate_key=0
runtime_seconds=2214`):

- **Every physics-bearing family cross-compares to a clean zero or a reviewed representational identity.**
  The geometry set (FACE_NORMAL, FACE_MEASURE_SHAPE_DERIV, FACE_VELOCITY, RELATIVE_FLUX, KINEMATIC_BALANCE),
  the projection set (PROJECTION_*), the evolution set (EVOLUTION_*), VIRTUAL_CONSTRAINT, and — after the
  fifth fix — **FACE_SHIFT (`join=160, py_only=0, wl_only=0, axis_set_mismatch=0`)** and
  **UNIFORM_LIMIT (`join=332, …=0`)** all join clean. TRACTION/CLOSURE/VIRTUAL_WORK reduce to the μ_θ
  identity; CONORMAL to the §3c verdict-A form.
- **The 11 families still showing unpaired cases are all pre-existing keying/schema asymmetries in the
  CONTROL and bookkeeping families — not physics disagreements**, and none was introduced by any fix:
  `ADMISSIBILITY_PREMISE` (the reserved S11c-b object, PY-missing-BRANCH); `REP_INVARIANCE_*` and
  `CONTROL_INDEPENDENCE_*` (PY-missing-DENSITY — the two engines key these controls on different axes);
  `CONTROL_FORM_*` (a FACE-granularity difference: PY keys `BOTH_FACES`, WL keys individual faces, plus a
  WL-missing-FACE); `DIMENSIONS` (one extra WL nested group). These are documented axis bridges, surfaced
  honestly by the accounting; reconciling those control-family keys is an owed item (below), not a finding.

So: **no genuine T7 physics disagreement survives** — established by the comparator's raw residuals plus the
documented per-family reconciliations (the μ_θ inert-collapse, the δρ identity, the FACE_SHIFT field-name
bridge in `scratchpad/faceshift_postfix_check.py`, the CONORMAL verdict-A, and the from-spec adjudications),
each carrying its own measurement — not a single classifier. (The committed hand-coded comparison script
`S11c_a_handcoded_comparison.py`, `bb2a050a`, packages the basic-pass-then-flag method but does not carry
every family's bridge, so it correctly FLAGs the bridge-reducible families rather than settling them; the
bridges are the reviewed REGISTRY folds below.) This is a narrower and more honest claim than the retracted
"39/39 clean."

## The five single-engine defects — each caught by the two-engine method (rule 1)

Each was a single-engine implementation error the sibling did not share; the disagreement was the
measurement. None was a spec ambiguity. Two (the current-freezes) were **substantive** — left in, they
would have carried a wrong current into S11c-b's coupling number. Three were **§3c free-premise /
faithfulness** defects that compute to a physically-forced value but made an emitted object diverge.

1. **Shifted-trace (PY)** `c36beac4`. PY mishandled the §3c trace-linearisation `δ[f(x,h_s)] = δf(x,h_s⁰) +
   δh_s ∂_w f⁰`; WL was right. Two SymPy bugs fixed.
2. **Background current (WL)** `6fae82b8` — *free premise*. WL carried the rest-frame background bulk current
   as §3c-forbidden free-premise symbols; the fix builds `j⁰=ρ⁰v⁰` from the supplied background, and the
   survivor is zero under static continuity (`∇₄·j⁰=0`).
3. **Projection current-freezing #3 (PY)** `49b5c525` — *substantive*. PY froze the perturbation current's
   **normal** component at the face, so `WINDOW_NORMAL_CURRENT` was identically 0 — the entire `∂_wδj_w`
   contribution absent, violating §1b. A from-spec adjudication (Codex + Grok, unanimous) ruled WL faithful,
   PY defective; PY was fixed to carry `δj_w`'s w-dependence.
4. **In-plane current #4 (PY)** `8c1a5ed1` — *substantive*. PY also froze the **in-plane** bulk-current
   divergence at the face (`current_divergence = Σ grad_j_bulk[i][i]`, bare); §1b makes the whole current
   `j=ρ_4D∇₄φ` a field of w and the projection integrates over w, so `∂_iδj_i` varies across the slab. The
   #3 fix had un-frozen only the normal jet. Fixed to `Σ Function(grad_j_bulk[i][i].name)(w)` at both spots.
   Cross-verified against WL: in all 5 PROJECTION families the current now matches WL's w-structure and
   cancels, leaving only the δρ density identity. This finding and #3 were the two the blanket collapse hid.
5. **Density-grounding (WL)** `cccb4f9e` — *free premise*. WL carried the T-e shifted-trace density
   background as an ungrounded, undefined free-premise function (`rhoBulkBackground`) — §3c forbids
   introducing a background normal derivative as a free premise, and §2b makes the density w-independent.
   PY (grounded) was faithful; WL was not, so FACE_SHIFT never cross-compared (`join=0`, WL missing the
   DENSITY axis). The fix grounds the density on the §2b representative `rho4Profile[density,coords]` via a
   real `D[·,w]`, threaded through all three shifted-trace sites (primary + form-control + uniform-limit),
   which gain the DENSITY axis. After it, FACE_SHIFT joins 16/16 (density) at residual 0 and 160/160 overall.
   The `MATERIAL_ADVECTED × RHOBR_CONSTANT` transport question this raised was adjudicated benign — see below.

Detail + literal operands: `directives/_measurements/S11c_a_faceshift_nonjoin.md` (#5) and the per-fix
measurements/commits.

## Cross-engine comparison — the hand-coded adjudication (rule 13)

The comparator prints raw residuals; classification happens on our side by hand. A naive string zero-test
gives false nonzeros (zero matrices/tuples, integral-linearity, applied-function-vs-symbol), and — the
lesson of this step — an **aggregation-side blanket collapse hides real dependence**. The correct instrument
is a basic mechanical pass (exact zeros + hand-verified pure renames) plus ONE explicitly verified inert
collapse (`μ_θ`), with everything else compared as-is and flagged. The surviving nonzero residuals fall into
these reviewed classes:

**Clean zero — comparator `A_minus_B` is literally a scalar/tuple zero.** FACE_NORMAL,
FACE_MEASURE_SHAPE_DERIV, FACE_VELOCITY, RELATIVE_FLUX, KINEMATIC_BALANCE, EVOLUTION_MASS_BALANCE,
EVOLUTION_TERM_ORIGINS, VIRTUAL_CONSTRAINT. (⚠ PROJECTION_*, FACE_SHIFT and UNIFORM_LIMIT **join clean but
their raw `A_minus_B` is a nonzero expression** that reduces to 0 only under a reviewed representational
bridge — they belong to the identity classes below, not here.)

**Representational identity — FACE_SHIFT (after #5).** All 160 cases join (`py_only=wl_only=0`); the raw
residual is nonzero because the two engines represent the traced perturbation differently. The **16 density
cases** reduce to exactly **0** under the field-name bridge (`rhoBulkPerturbation ↔ delta_rho_4D_face`,
`e_W ≡ δW/W_0 ↔ eta_bg·w1_profile`) — verified, `directives/_measurements/S11c_a_faceshift_postfix_check.py`
("nonzero density residuals: 0"). The **64 current cases** carry the identical **applied-at-the-face-point ↔
bare-symbol** difference: WL writes the current perturbation as `delta_j_bulk_i(x1,x2,x3,(sW_0/2,time))`
applied at the trace point (plus its `_dw` jet), PY as the bare face symbol `delta_j_bulk_i` (plus its `_dw`
jet) — `directives/_measurements/S11c_a_faceshift_current_check.py` shows the residual is exactly that
naming/application difference (`delta_j_bulk_i(x,face) − delta_j_bulk_i` + jet rename), the same μ_θ-class
applied↔bare pattern, since the T-e trace is a point evaluation at the shared in-plane `x`, not an
`x`-derivative. Physics identical, representation different. The committed hand-coded comparison
(`S11c_a_handcoded_comparison.py`, `bb2a050a`) correctly **FLAGs** all of these — per its own rule it flags
rather than massages; the bridges are the reviewed REGISTRY folds listed below.

**Representational identity — the closure coefficient μ_θ** (TRACTION, CLOSURE_SHAPE_DERIV,
VIRTUAL_WORK_SHAPE_DERIV). WL writes `μ_θ^α` as an applied function `mu_theta_{L,M}(x1,x2,x3,time)`; PY as a
bare symbol. The whole residual is `Σ coeff·(mu_theta − mu_theta(x1,x2,x3,time))`, collapsing to 0 for all 48
cases. Decisive probe: across the whole run there are **zero** jet-suffixed `mu_theta` and **zero**
`Derivative(mu_theta)` — WL never differentiates μ_θ. The spec makes the evaluation-point arguments inert:
§2a supplies in-plane derivative maps for `W_bg` and `μ_R,bg` but **not** `μ_θ`; §3c forbids a free-premise
background derivative; §4 T-i keeps `μ_θ` a **named** operand; μ_θ enters only algebraically through
`μ_s^α=μ_θ^α/ρ_br,bg^{0,α}`. Bare ≡ applied. This is the ONE verified-inert collapse the hand-coded
comparison is permitted to apply. Benign.

**Representational identity — the projected density perturbation δρ_4D** (PROJECTION_*). After #3+#4 and
integral linearity the in-plane current divergence `∂_iδj_i` cancels between the engines; the remainder is a
density-time term where PY carries a named symbol `delta_rho_4D_bulk_t` and WL carries its explicit §3a
expansion. Binding PY's symbol to the spec-built `∂_tδρ_4D = ∂_t[ρ_4D,bg^{0,α}·θ]` sends the residual to 0
for all four (anchoring × representative) combinations, including MATERIAL_ADVECTED × RHOBR_CONSTANT, where
the advected background `ρ_4D,bg⁰(χ(x,t))` contributes a first-order advection term. `delta_rho_4D_bulk` is a
**primitive** symbol in PY (never expanded); WL, blind, expands it. That this named symbol is the *full* §3a
perturbation (not the local-advected shortcut) is closed by measurement, not interpretation: PY's EVOLUTION
family builds the density-time concretely (`sigma_exact` + an explicit `BACKGROUND_ADVECTION = ε·u_t·∇σ0`,
engine L1298) and EVOLUTION_MASS_BALANCE + EVOLUTION_TERM_ORIGINS agree with WL **8/8 to 0**, PY's
MATERIAL_ADVECTED operand carrying the `u_{i,t}` advection. Benign.

**CONORMAL_DERIV** — the already-adjudicated §3c verdict-A form: a probe-representation / Taylor-order
difference (`W_0²`, third order) in the unmapped `conormalBackground` jets, same physics, no finding
(VERDICT A: `directives/_measurements/S11c_a_comparator_reemit_plan.md:82`, two from-spec CAS legs each
residual 0).

**Supplied / bookkeeping families — structural or naming differences, not derived physics.**
BACKGROUND_STATE, FACE_MAP_{LAB_HELD,MATERIAL_ADVECTED}, BACKGROUND_DENSITY_MAP, and ADMISSIBILITY_PREMISE
surface as `Mismatch(STRUCTURE_DISAGREE)` (PY tuple/scalar vs WL Association), unequal-length tuples, or a
naming difference (`w1Profile ↔ w1_profile`). These are objects the spec **supplies** (§2b/§2d/§3a) and each
engine emits as-is; ADMISSIBILITY is the reserved S11c-b object (§2d, not tested here). Their differences are
representation/coverage, surfaced honestly by the accounting — not a physics disagreement to reconcile in
S11c-a.

## FACE_SHIFT and the MATERIAL_ADVECTED transport adjudication

The fifth fix (`cccb4f9e`) grounded the WL shifted density trace, and FACE_SHIFT now joins 160/160 clean.
The build review of that fix (leg 2, Grok) raised a sharper question: §3a defines the MATERIAL_ADVECTED
background density as `ρ_4D,bg⁰(χ)`, yet both engines emit the MATERIAL×RHOBR density trace with **no
material-transport term** (identical in form to LAB_HELD). Because both engines agree, the cross-engine
residual (0) cannot adjudicate it — a **rule-7 shared-blind-spot candidate** (an error in a spec both
engines read makes them agree on the same wrong thing). Two from-spec adjudication legs (fresh Claude + Grok,
each with SymPy + literal stdout) settled it: **the omission is correct.** §3c's trace law is normal-shift
only (`δh_s·∂_w f⁰`) and `∂_wρ⁰=0`, so the in-plane transport `−u·∇ρ⁰` is not a §3c face-shift term; it
enters `δh_s` (the MATERIAL-branch face displacement, where it multiplies `∂_wρ⁰=0` and contributes nothing
to the density trace) and the density representative consumed by T-f/T-g/T-h; injecting it into T-e would
double-count. Not a finding — both engines correct.

One non-physics residual is **recorded, not fixed**: WL's `EXACT_TRACE_SOURCE` sub-field (the
un-differentiated trace) displays the MATERIAL background at lab `ρ(x)` where §3a says `ρ(χ)`. Both legs
confirmed this changes **no computed object** (at background `u=0 ⇒ χ=x`; `ρ(χ)` also has `∂_w=0`; the
shape-derivative operand is frame-independent). It is a display-fidelity point in an intermediate that is not
cross-compared. (`directives/_measurements/S11c_a_faceshift_nonjoin.md` §8–10.)

## The control battery — corrected characterization

The earlier draft's "every object bites, no dead control" was **wrong** — it came from a string liveness
test that counted sentinel operands (`TextAtom('UNDEFINED_UNJOINED')`, `Mismatch`) as bites. A correct
liveness test (a genuine bite = `operand_A` is a nonzero sympy expression carrying a `Symbol`; sentinels
and literal zeros excluded) run against the final committed-state run (`scratchpad/bite_liveness.py`,
stdout `~/.s11_build/bite_liveness_cccb4f9e.out`; the stale claim in
`directives/_measurements/S11c_a_control_battery_result.md` is superseded by its own CORRECTION section):

- **HOLD controls hold — 0 genuine nonzero cross-engine.** REP_INVARIANCE_RESIDUAL (N4/N6 route-invariance)
  has 0 nonzero residual wherever the two engines join; its unpaired cases are the keying asymmetry
  (PY-missing-DENSITY). UNIFORM_LIMIT_RESIDUAL (the S11c-a→S11b regression) is **fully joined after
  `cccb4f9e`** (`join=332`, no unpaired) at 0 nonzero — the fifth fix closed its former WL-missing-DENSITY
  gap. (These replace the pre-fix `S11c_a_control_hold_classify` counts, which predate `cccb4f9e`.)
- **BITE controls — CONTROL_INDEPENDENCE bites all 6 objects it covers** (CLOSURE, EVOLUTION_MASS_BALANCE,
  RELATIVE_FLUX, TRACTION, VIRTUAL_CONSTRAINT, VIRTUAL_WORK). **CONTROL_FORM (the `∂W_bg` jet ablation) bites
  16 of 18 objects** — the geometry set, FACE_VELOCITY, EVOLUTION_*, VIRTUAL_CONSTRAINT, and the PROJECTION
  operands — and is **DEAD (no bite) for exactly two: FACE_SHIFT and PROJECTION_STATIC_OPERAND**, whose
  operands carry the `W_bg` profile *value* (`w1_profile`) but not its in-plane *jet* (`w1_profile_d`), so a
  jet-direction ablation cannot move them. So neither "every object bites" nor "only geometric objects bite"
  is right. FACE_SHIFT is therefore not exercised by the §5 battery; it is validated by the cross-engine
  residual (16/16 → 0 post-`cccb4f9e`), the density-grounding build-leg form ablation, and the
  material-transport adjudication instead. A re-run CONTROL_FORM coverage table over the abstract-symbol
  objects is an owed item; the bite/dead split above is what the corrected liveness test measures.

## The reconciliation schema — pinned for S11c-b..e

S11c-a fixes the cross-engine reconciliation schema later sub-steps inherit (decision list N8):

1. **Axis-typed join keys.** Cases join on the object name plus an axis-typed key (BRANCH, DENSITY, FACE,
   DOF, DIRECTION, OBJECT), not a full tuple — PY's bare `DELTA_W` and WL's `DOF_DELTA_W` are the same axis.
2. **Never blanket-collapse to judge representational-vs-finding.** The classifier is a basic mechanical pass
   (exact zeros + hand-verified pure renames) plus a hand-coded layer; an applied bulk quantity is mapped to
   a bare symbol only for the classification zero-test, **only where verified inert** (μ_θ), and never in the
   comparator's operands. ⛔ Mapping an `AppliedUndef→Symbol` in construction is the argument-strip that hid
   the current-freezing and density-grounding defects.
3. **Integral linearity is the canonicalizer for projection residuals.** Pull binder-free factors out of the
   bound integral and combine same-`(lower,upper)` integrands before any zero-test.
4. **Every background face value / normal derivative must be grounded in the supplied state `𝔅⁰`, never a
   free premise** (§3c). The bg-current (`6fae82b8`) and density-grounding (`cccb4f9e`) fixes are both this
   rule; a free-premise background is a single-engine defect the residual only catches once the sibling
   grounds it.
5. **The closed representational folds** (reviewed REGISTRY entries, not comparator name-map surgery): the
   closure coefficient `μ_θ^α` (applied-at-face ↔ bare, inert args); the projected density perturbation
   `δρ_4D` (named symbol ↔ §3a expansion, including the MATERIAL_ADVECTED advection); the CONORMAL §3c
   verdict-A form. The material advection belongs to `δh_s` + T-f/T-g/T-h, **not** the T-e trace.

## Standing limits, freezes, and owed items

- **The background-flow correction `O(v_bulk_normal_0·|q_n|/ω)`** is inherited as a standing limit from S11b,
  uncarried and unbounded where `|q·v_bulk_normal_0/ω|≳1` (S11c decision list §E carry-in). S11c-a computes
  geometry only and does not address it. (`v_bulk_normal_0` is the inert rest-frame scope limit of §1, per
  §3c — not a distinct symbol `v₀`.)
- **Owed: control-family keying bridges.** The final run leaves 11 control/bookkeeping families with unpaired
  cases from keying-convention differences (REP_INVARIANCE/CONTROL_INDEPENDENCE PY-missing-DENSITY,
  CONTROL_FORM FACE-granularity, ADMISSIBILITY reserved, DIMENSIONS envelope). None is a physics
  disagreement, but the control cross-comparison is incomplete until those axes are bridged; and a post-fix
  control-form re-characterization is owed (above). These are the honest residual of S11c-a's cross-check.
- **Falsification is deferred to S11c-e.** S11c-a is the geometric substrate; the flux-normalized
  dimensionless conversion FORM and its bench-optics bound are the S11c-e deliverable (decision list table).
- **The family roll-up card** closes S11c after all sub-steps (decision list N1); S11c-a owes no separate
  paper card, only its exports and this record.

## Operational note

The comparator runs ~37 minutes; heavy PY/WL runs were detached (`setsid`) against spurious background-job
kills, with a watcher loop. Review and consult legs were run with one Grok instance at a time (the single
per-user lock) and Codex at `xhigh` with `--sandbox danger-full-access`. The comparator's build followed the
rule-15 pivot: three prose build-directive rounds each drew fresh independent-leg defects (39 bespoke
per-family schemas cannot be pre-enumerated in prose), so the build was delegated via a brief that locked the
reviewed physics folds and let the two build legs gate the working instrument. Each of the five fixes went
through its own directive→legs→build→legs→commit loop (WL fixes as Codex builds with a regenerated
transcript; the density-grounding fix directive itself drew a scope-gap catch from its Codex/Grok directive
legs — the two extra shifted-trace sites — folded before the build).

## Provenance

- Spec `directives/S11c_a_SHARED_PHYSICS.md`; decision list `directives/S11c_decisions.md`; scope
  `steps/S11c_SCOPE.md`.
- PY engine `scripts/S11c_a_interface_geometry_sympy_audit.py` (`c36beac4`, `49b5c525`, `8c1a5ed1`);
  transcript `scripts/out/S11c_a_interface_geometry_sympy_audit.out` (`8c1a5ed1`).
- WL engine + transcript `mathematica/{,out/}S11c_a_interface_geometry_mathematica_audit.{wl,out}`
  (`6fae82b8`, `cccb4f9e`).
- Comparator `scripts/S11c_a_cross_engine_comparator.py` + tests (`50f43123`); hand-coded layer
  `scripts/S11c_a_handcoded_comparison.py` (`bb2a050a`); final run `~/.s11_build/comparator_final_cccb4f9e.out`;
  per-family bridge/liveness scripts
  `directives/_measurements/S11c_a_{faceshift_postfix_check,faceshift_current_check,bite_liveness}.py` +
  runs `~/.s11_build/{handcoded,bite_liveness}_cccb4f9e.out`.
- Measurements `directives/_measurements/S11c_a_{faceshift_nonjoin,control_battery_result,
  T7_adjudication_verdicts}.*` and the per-fix `_measurements/` + `_legs/` finding records under
  `directives/` and `research/pde_ledger_v3/_measurements/`.
