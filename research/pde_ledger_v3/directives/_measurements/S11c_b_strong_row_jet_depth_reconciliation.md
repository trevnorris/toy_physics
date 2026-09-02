# S11c-b strong-row jet-depth reconciliation (the #89b PY-check FLAG) — IN PROGRESS

The #89b PY sibling freeze-check (`_measurements/S11c_b_89b_py_sibling_freeze_check.md`) surfaced a FLAG:
PY caps strong rows at `STRONG_ROW_JET_DEPTH = 2` (Hessian) while WL #89b "emits 3rd-order background jets in the
strong U-momentum rows," so `scripts/S11c_b_row_residual.py` would show a strong-row disagreement to adjudicate
as a spec question BEFORE reading it as physics. This record reconciles it. ⛔ Do NOT pre-judge which engine.

## Step 0 — the σ_W-grade cannot be the discriminator (source read, both settled)

- PY `background_jet_expression` (`sympy_audit.py:681-693`): a background jet of ANY order n carries exactly ONE
  power of `sigma_W` (`sigma_W * profile_atom / L_W**(n-1)`); orders differ only by powers of `1/L_W`.
- PY `retained_grade` → `first_shape_series` truncates ONLY `eta_bg <= 1` and `sigma_W <= 1` (`:883`, `:911`).
- ⇒ the retained grade does NOT separate an order-2 Hessian jet from an order-3 jet (both are σ_W¹). Whether
  depth-2 is complete is therefore a question of whether a strong row GENERATES nonzero order-3 jets, which is
  computable — not a grade convention. (WL order = `Total` of the index triple, `mathematica…:1857`, so e.g.
  `widthProfileJet[1,1,0]` is order 2, `widthProfileJet[2,1,0]` is order 3.)

## Step 1 — PY side: the strong rows are ORDER-2-COMPLETE (measured, self-verified; rule 13)

PY already has the instrument: `task_tower_depth_control` (`sympy_audit.py:3731`) compares retained-grade strong
rows built at background_depth 2 vs 3 vs 4. It is in `SKIPPED_TASKS` of the committed `.out` (PRIMARIES_ONLY), so
it was re-run in isolation on the committed engine for one representative case (LAB_HELD, RHO4_CONSTANT,
EULERIAN), depths 1/2/3/4, each through `retained_grade`.

**Command:** `python3 tower_depth_probe.py` (imports the committed engine; builds
`retained_grade(live_strong_rows(construct_energy(branch, background_depth=d).density, …, d))` for d∈{1,2,3,4};
prints pairwise residuals + the max σ_W power and background-jet orders present; ephemeral script
`…/scratchpad/tower_depth_probe.py`, engine `research/pde_ledger_v3/scripts/S11c_b_brane_operator_sympy_audit.py`).

**Literal stdout:**
```
depth1 leaves=5   depth2 leaves=5   depth3 leaves=5   depth4 leaves=5

=== RESIDUAL_21_HESSIAN_SANITY ===        (depth2 - depth1: is the instrument LIVE?)
RESIDUAL_21_HESSIAN_SANITY_NONZERO_LEAF_COUNT: 5 of 5
RESIDUAL_21_HESSIAN_SANITY_MAX_SIGMA_W_POWER: 1
RESIDUAL_21_HESSIAN_SANITY_BG_JET_ORDERS_PRESENT: [1, 2]

=== RESIDUAL_32_DECISIVE ===              (depth3 - depth2: the WL(3) vs PY-default(2) question)
RESIDUAL_32_DECISIVE_NONZERO_LEAF_COUNT: 0 of 5
RESIDUAL_32_DECISIVE_RESIDUAL: 0  (identical after retained grade)

=== RESIDUAL_43_TERMINATION ===           (depth4 - depth3)
RESIDUAL_43_TERMINATION_NONZERO_LEAF_COUNT: 0 of 5
RESIDUAL_43_TERMINATION_RESIDUAL: 0  (identical after retained grade)

=== UNREDUCED DEPTH-3 STRONG ROWS ===     (generated-then-cancel, or never generated?)
UNREDUCED_DEPTH3_BG_JET_ORDERS_PRESENT: [1, 2]
UNREDUCED_DEPTH3_ORDER>=3_ATOMS (0): []
```

**Reading (rule 2 — the record interprets, not the script):**
- The instrument is LIVE: RESIDUAL_21 is 5/5 nonzero, so depth-1 (freezing at the first jet) genuinely drops the
  Hessian — the zeros below are real, not a dead probe.
- RESIDUAL_32 = 0 and RESIDUAL_43 = 0: the retained-grade strong rows are IDENTICAL at depth 2, 3, and 4. PY's
  `STRONG_ROW_JET_DEPTH = 2` drops NO retained-grade term — it is NOT a freeze; it is complete.
- The un-reduced depth-3 strong rows contain NO order≥3 atom at all: order-3 is never GENERATED in the strong
  rows (not generated-then-cancelled). The strong (bulk-EL) U-rows structurally stop at the Hessian (order 2).
- ⇒ the FLAG's presumption ("PY under-emits at depth 2") is REVERSED: PY's bulk strong rows are order-2-complete.

## Step 2 — WL side: U_MOMENTUM_ROWS is a 3-slot SUM; order-3 enters via a non-bulk slot (OPEN)

WL emits order-3 in its `U_MOMENTUM_ROWS` (durable evidence: the #89b build-review Finding-1 orchestrator
verification found `widthProfileJet[2,1,0], [3,0,0], [1,2,0]` — all order 3 — in the correct U-row emit, dropped
by the re-freeze bug: `_measurements/S11c_b_89b_wl_operator_build_review.md`). But WL builds
`U_MOMENTUM_ROWS = kineticU + rows["U_INTERNAL"] + faceRows["U_FACE"]` (`mathematica…:1170`), separated in
`origins` as KINETIC / BULK_ENERGY / FACE (`:1182-1189`). PY's bulk strong rows correspond to `U_INTERNAL`.
The cross-engine comparator `scripts/S11c_b_row_residual.py` already attributes and EXCLUDES the FACE part from
the strong-row residual (`FACE_ATTRIBUTED_AND_EXCLUDED_FROM_RESIDUAL`) and marks the KINETIC adapter
`PROVENANCE_ONLY_ALREADY_IN_COMPLETE_ROW`.

**HYPOTHESIS (was): the order-3 comes from FACE/coupling, not bulk `U_INTERNAL` — an object-boundary artifact.
⇒ REFUTED by measurement.**

**WL provenance measurement (rule 16; self-run, reduced scale):** loaded the committed WL engine's definitions
(proven StringTake-between-markers loader, `basisRepresentativeIndices = {16}`), evaluated
`evaluatedModel["EULERIAN","LAB_HELD","RHO4_CONSTANT"]`, and reported the engine's OWN `backgroundJetOrder`
(`mathematica…:1857`, order = Total of the jet index) over each provenance slot of the reduced `ORIGINS`.

**Command:** `wolframscript -file wl_provenance_probe.wl` (ephemeral `…/scratchpad/wl_provenance_probe.wl`; engine
`research/pde_ledger_v3/mathematica/S11c_b_brane_operator_mathematica_audit.wl` @ HEAD).

**Literal stdout:**
```
ORIGINS_SLOT_KEYS: {"KINETIC", "BULK_ENERGY", "FACE", "FLUX", "ADVECTIVE", "ACCUMULATION"}
FULL_U_MOMENTUM_ROWS_MAX_JET_ORDER: 3
FULL_U_MOMENTUM_ROWS_N_ATOMS: 19
FULL_U_MOMENTUM_ROWS_ORDER3PLUS_ATOMS: {widthProfileJet[2,0,1], widthProfileJet[1,0,2],
  widthProfileJet[1,1,1], widthProfileJet[2,1,0], widthProfileJet[1,2,0], widthProfileJet[3,0,0],
  widthProfileJet[0,1,2], widthProfileJet[0,2,1], widthProfileJet[0,3,0], widthProfileJet[0,0,3]}
ORIGIN_KINETIC_U_MAX_JET_ORDER: 0        ORIGIN_KINETIC_U_N_ATOMS: 0
ORIGIN_BULK_ENERGY_U_MAX_JET_ORDER: 3    ORIGIN_BULK_ENERGY_U_N_ATOMS: 19   (all 10 order-3 atoms here)
ORIGIN_FACE_U_MAX_JET_ORDER: 0           ORIGIN_FACE_U_N_ATOMS: 0
FULL_MINUS_SLOTSUM_U (expect 0 vector): {0, 0, 0}
```

**Reading:** the order-3 jets sit ENTIRELY in the `BULK_ENERGY` slot (WL's `U_INTERNAL`) — the exact slot that
corresponds to PY's bulk strong rows. FACE and KINETIC carry NO background jets for this case; the slot sum
reconstructs the full U-row exactly ({0,0,0}). So this is NOT an object-boundary artifact.

## VERDICT (both sides measured, rule 13): a GENUINE cross-engine disagreement in the BULK momentum EL rows

- **WL** bulk U-momentum EL rows: max background-jet order **3** (10 order-3 `widthProfileJet` atoms).
- **PY** bulk strong U-momentum rows: max background-jet order **2**; order-3 is NEVER generated (depth-invariant
  2=3=4). PY's `STRONG_ROW_JET_DEPTH = 2` is not a truncation of a tower PY would otherwise fill — PY structurally
  applies one fewer background differentiation than WL to form the bulk momentum row.

This is a real finding (rule 1/6), and `scripts/S11c_b_row_residual.py` WILL show it as a bulk strong-row
disagreement. ⛔ Do NOT pre-judge which engine is correct. Note the #89b context both ways: (a) #89b's whole point
was that WL wrongly FROZE the background (reduce-before-differentiate) and the fix restored the live jets — so
WL's order-3 may be the correct consequence of keeping W_bg live, making PY's depth-2 the analogous freeze
(rule 17: a jet/derivative order held fixed IS the finding); OR (b) the un-freeze may OVER-restore an order-3 that
the physical retained grade does not carry, making PY's depth-2 correct. The measurement cannot yet say which.

## NEXT — adjudication (a SPEC question; rule 7 = 2 decision legs before any engine change)

Localize the ROOT before deciding: measure the max background-jet order of the ENERGY DENSITY (pre-EL) on BOTH
engines. If WL's energy density already carries order-2 (Hessian) where PY's carries only order-1 (grad), the gap
is upstream in the basis/energy construction (the EL then adds one → 3 vs 2); if the energy densities agree, the
gap is in the EL/divergence row-formation itself. Then frame the spec question — how many background
differentiations a bulk momentum EL row physically carries, and to what retained jet order — and put it through
two decision legs. Definitive in-band cross-engine confirmation is the deferred `S11c_b_row_residual.py` on both
regenerated `.out` (≥64 GB box, DEFERRED_HEAVY_RUNS.md).

## Step 3 — Codex + Grok consultation (orchestrator-written brief → the two-advisor pair): the ROOT

Brief: `directives/_legs/S11c_b_jet_depth_reconciliation_consult.md`. Both advisors derived independently (scripts +
literal stdout committed under `_measurements/s11c_b_jet_depth_consult_codex/` and `…_grok/`; raw verdict
transcripts at `~/.s11_build/S11c_b_jet_depth/consult_{codex,grok}_verdict.log`, `.log` being gitignored). They
AGREE with each other and with the orchestrator's own PY probe — three consistent computations.

**CONSENSUS (all three; the earlier "genuine disagreement" is re-diagnosed, NOT a jet-depth freeze):**
1. The energy basis carries only FIRST-jet spurions (`∇W_bg`, `∇μ_R,bg`); NO Hessian is an energy invariant
   (Grok PY: `DENSITY_BG_JET_ORDERS: [1]`, 40 terms). ⇒ the ∇u-flux coefficient is order 1 (p=1), and a single EL
   divergence raises it to order **2**. The UNCONSTRAINED bulk momentum EL row `δL/δu` is **order-2-complete on
   BOTH engines** (Codex+Grok measured `RAW_U_ENERGY_EL_MAX_JET_ORDER: 2` in WL too; PY depth 3 generates nothing).
2. PY's `STRONG_ROW_JET_DEPTH = 2` is **not a rule-17 freeze** for this object; raising it is a measured no-op.
3. WL's order-3 is ENTIRELY the θ-constraint reaction `∇μ_θ`: `constrainedRowsWithLiveEnergyEL` eliminates virtual
   θ via the linearized mass/interface constraint (`θ + (W0/W) e_W + u·∇W/W + Div u`), leaving an `Inactive[Div]`
   of a flux that IS `μ_θ`; #89b correctly un-froze that outer Div ⇒ one EXTRA derivative of the θ-EL ⇒ order 3
   (the 10 atoms). "Un-freezing a derivative of the wrong object still yields the wrong object" (Grok).
4. PY keeps θ SEPARATE (`operator_from_density` returns `(operator, mu_theta_amplitude)`; the U-row is the
   unconstrained EL). Cross-check: PY's OWN constraint-reduced row reproduces WL's 10 order-3 atoms EXACTLY
   (Codex PY `CONSTRAINED_U_RETAINED` = `w1_profile_d1d1d1…d3d3d3` ↔ WL `widthProfileJet[3,0,0]…[0,0,3]`). So the
   engines reproduce each other's objects — there is NO physics error, only a REPRESENTATION mismatch.
5. ⇒ `scripts/S11c_b_row_residual.py` on the bulk U-row WILL disagree, and reading it as "jet-depth convention" or
   "PY under-emits Hessians" is the WRONG diagnosis. The residual is the off-shell `∇μ_θ` constraint reaction.

## THE DECISION (orchestrator + Codex + Grok): spec-pin the OBJECT, then 2 decision legs, then align — NO engine
## change first

The remaining question is NOT undecidable and NOT a depth knob: it is **which object the slab U-momentum row is
specified to be** —
- **(a) unconstrained energy EL `δL/δu`, θ an independent DOF** (§3b lists `{u, θ, e_W}` independent): then PY's
  representation is the target; WL should NOT activate the constraint `Div` into `BULK_ENERGY`/`U_INTERNAL` (those
  `∇μ_θ` terms are already the θ-row). Order-3 then leaves the bulk U-row on both engines.
- **(b) constraint-reduced row, virtual θ eliminated**: then WL's representation is the target for that (different)
  object; PY would have to fold the constraint, and ONLY THEN would PY depth-2 become a genuine freeze (since
  `∇μ_θ` is order 3). That is a spec change of the object, not a depth repair.

Codex's concrete form of the same decision: emit `RAW_BULK_U_EL` and `CONSTRAINT_REACTION_U` as SEPARATE objects
so the comparator aligns like-with-like. Grok's lean (to be settled by the legs, not pre-decided): §3b's
independent-DOF listing points to (a) / PY, with WL mislocating `∇μ_θ` into the bulk slot. ⛔ Do NOT pre-decide;
this ties into the S11b interface/coupling law and #88 (KINETIC+θ) — the pin must be consistent with those.

## Status: RECONCILED, then PINNED. NOT a jet-depth physics disagreement — a representation question, now decided.

**SPEC-PIN RESULT (decision list `directives/S11c_b_jet_depth_spec_pin_decision.md`, 2 legs Codex+Grok, raw verdict
transcripts `~/.s11_build/S11c_b_jet_depth/{codex,grok}_leg_verdict.log`): PINNED (B).** The slab U-momentum row is the
constraint-reduced in-plane equation carrying `−∇μ_θ`; θ's virtual displacement is NOT independent. The earlier
lean toward (A) is OVERTURNED. Decisive governing text (orchestrator-verified verbatim): `S11b_SHARED_PHYSICS.md`
`:337` ("Do NOT vary U with θ held fixed") and `:426` ("the in-plane equation must carry `−∇(δU/δθ)` … varying at
fixed θ removes this contribution … selects the convention uniquely"); S11b's engine `…coupling_law_mathematica_
audit.wl:280` (`constrainedUL = explicitUL + I k muTheta`). ⇒ **WL is correct; PY is the engine that must change**
(fold S11b's virtual constraint into `U_MOMENTUM_ROWS` + thickness row, raise `STRONG_ROW_JET_DEPTH` 2→3, keep
`μ_θ` a separate constitutive operand, θ-row = mass evolution); §3b gets a clarifying sentence; #88's raw-EL
identification of the U-row conflicts with the pin (its basis-completion measurement stands, its full-row
adjudication is redone after the PY fold). Full folded verdict + engine target: the decision list's FOLDED VERDICT
section. NEXT = PY constraint-fold BUILDER (own decision list + build legs, rule 7/9); then deferred `row_residual`
on both `.out`. ⛔ No engine change yet.
