# ledger_stage039_b_t_time_reversal_even_departure

## Status

**Part V — Magnetism. V-6 (build-order 039; the SIXTH and FINAL stage of the 6-stage Part-V
split, user decision 2026-07-22) — ⭐⭐ the CHARACTERIZED DEPARTURE the magnetism sector carries.**
The self-contained parity-level characterization of the candidate magnetic field
`b_T = ∇×u_T`. Where stage 034 (V-1) paid for the moving-throat action row, stage 035 (V-2)
derived the native source law + the parity census, stages 036/037 (V-3/V-4) ran the two blind
far-field routes + the boost-consistency structural relation, and stage 038 (V-5) landed the
sector's honest terminal `R1_REQUIRED(electric_bc_selection)`, stage 039 asks the complementary
question — **how does the model's magnetism DEPART from exact Maxwell?** — and records the honest,
first-class answer:

- **`B_TIME_REVERSAL_EVEN`** — the candidate magnetic field `b_T = ∇×u_T` is time-reversal
  **EVEN** and **axial**, whereas a physical Maxwell `B` is time-reversal **ODD** (and axial). The
  parity is **DERIVED** from `u_T`'s: the curl is a spatial-derivative operator that maps a
  **polar** vector to an **axial** vector and **preserves** time-reversal parity, so applying the
  curl-parity map to `u_T`'s census tuple `(R_w, P_w, T, rotation) = (−1, −1, +1, polar_vector)`
  (035, `PARITY_ROTATION` / `PARITY_TIME_REVERSAL`) gives `b_T = (−1, −1, +1, axial_vector)` —
  axial + T-even, matching 035's cited census. Against the canonical Maxwell benchmark
  `{time_reversal: −1, rotation: axial_vector}` the two **DIFFER in the time-reversal axis ONLY**
  (`b_T` is correctly axial, like a real `B`), so the departure is a genuine T-parity effect, not a
  rotation artifact. The T-even parity is a **consequence** of the throat's **active-drain
  time-arrow `τ_d`** (T-odd; 034 / R68) through the source `J_T = q_T s η V` — a passive T-even
  throat would supply no `O(V)` moving row at all. Reads out as: **"magnetism is NOT exact
  Maxwell."** Verdict token **`B_TIME_REVERSAL_EVEN`** (both engines, both exit 0).

**⚠ CRITICAL TYPING (do not overstate).** This is a **CHARACTERIZED-DEPARTURE, recorded
first-class** — a load-bearing prediction, NEVER softened, NEVER "rescued." It is **not** a bug,
**not** a knob, **not** a reduction, and **not** an R1 deferral. It discharges **no** knob and
**does not shrink the irreducible count** (`docs/model_map.md` §4, the honest departure ledger). It
is the **magnetic twin** of the charge-sector departure **`NATIVE_P_NO_EMERGENT_GAUSS`** (033 /
R66, "exact U(1)/Maxwell Gauss is non-native") and the **sibling** of the light-sector departure
**`FAIL_CAUCHY_STRAY_LONGITUDINAL`** (Part III / stage 003 / pathA_36, "a second-class pair, NOT
Maxwell's first-class Gauss") — the three first-class "EM is NOT exact Maxwell" departures the
model carries.

**⚠ 039 CHARACTERIZES the departure; it does NOT re-derive the census.** Like 033 was the
**quadratic** Stage-1 constraint gate (not an all-orders proof), 039 is the **parity-census-level**
departure: the `b_T` parity FACTS (`b_T` axial via `PARITY_ROTATION`, `u_T` T-even via
`PARITY_TIME_REVERSAL`) are DECIDED at stage 035 and recorded RAW there. 039 **cites** those facts
and supplies **only** the magnetism-NEW content: the explicit Maxwell-`B` comparison, the T-axis
localization, the active-drain time-arrow requirement, and the internal self-consistency of the
T-even parity with the derived census + 038's terminal landing. It does **NOT** re-run the four
`PARITY_*` teeth, re-derive the source law, re-adjudicate the 038 landing, or re-count
`u_T`/`τ_d`/`q_T`/`b_T`/the census in the irreducible set.

**⚠ CRITICAL SCOPE — 039 records the PARITY departure; it is DISTINCT from 038's post-resolution
truth-table cell.** 038's SEALED §4 truth table *contains* a `MAGNETISM_DEPARTURE_CHARACTERIZED`
landing (**6 cells**), but that landing is a **post-anchor-resolution** cell reached only when
`anchor=closed ∧ comparison=routes_differ` — NOT 038's production verdict (which is
`R1_REQUIRED(electric_bc_selection)`). 039's `B_TIME_REVERSAL_EVEN` is the **actual
characterization** of the departure, decided at the parity level; the sector *magnitude/sign* stays
R1 (038 / R72) — the two live on **orthogonal axes**. A reader must never conflate 039's
parity departure with 038's post-resolution §4 cell, nor read either as resolving the sign.

**⭐⭐ This COMPLETES Part V.** After 039, magnetism is closed: 034 (action row, R67/R68) → 035
(native source + census, R69) → 036/037 (the two blind routes + boost-consistency, R70/R71) → 038
(sealed terminal R1 landing, R72) → **039 (the first-class characterized departure, R73)**.

## Purpose

Record, as a first-class characterized departure, that the model's candidate magnetic field
`b_T = ∇×u_T` **departs from exact Maxwell in its time-reversal parity**: `b_T` is T-EVEN whereas a
physical Maxwell `B` is T-ODD. The decisive objects are (i) the curl-parity map that **derives**
`b_T`'s parity from `u_T`'s, (ii) the explicit typed comparison against the canonical Maxwell
benchmark, (iii) the T-parity propagation chain `τ_d → q_T → J_T → u_T → b_T` grounding the T-even
result in the active-drain time-arrow, and (iv) an internal two-path cross-check that both
derivations of `b_T`'s T-parity agree and are orthogonal to 038's terminal R1 landing. Both engines
re-instantiate 035's census entries `{u_T, b_T, τ_d, q_T, s, V, J_T}`, apply the curl-parity map,
compare against the Maxwell benchmark, propagate the active-drain chain, and reach
`B_TIME_REVERSAL_EVEN` on their own — through the **identical** pipeline

```text
035 census entries  →  curl-parity map  (u_T polar,T-even → b_T axial,T-even)
                    →  Maxwell comparison  (b_T T-even ≠ Maxwell B T-odd; departure localized to T)
                    →  active-drain propagation chain  (τ_d T-odd → … → b_T T-even)
                    →  two-path self-consistency + orthogonality to the 038 R1 landing.
```

Consumes **nothing new numeric** — it **imports by citation** 035's parity census (`u_T` T-even,
`b_T` axial), 034's `τ_d` (R68), and 038's terminal landing (R72), and supplies only the departure
characterization. It is the charge-sector twin of `NATIVE_P_NO_EMERGENT_GAUSS` (033) and the
light-sector twin of `FAIL_CAUCHY_STRAY_LONGITUDINAL` (stage 003).

## 1. `b_T`'s parity is DERIVED from `u_T` via the curl (`B_T_AXIAL_T_EVEN`)

`b_T = ∇×u_T` is a spatial curl of the physical brane displacement `u_T`. The curl is a
spatial-derivative operator: it maps a **polar** vector to an **axial** vector and **preserves** the
time-reversal parity. From 035's census (`PARITY_ROTATION` / `PARITY_TIME_REVERSAL`), `u_T` is
**polar** and **time-reversal EVEN**, with census tuple

```text
u_T  = (R_w, P_w, T, rotation) = (−1, −1, +1, polar_vector).
```

Applying the curl-parity map (polar → axial; T-parity preserved) to `u_T`'s tuple **derives**

```text
b_T = ∇×u_T = (−1, −1, +1, axial_vector)      — axial + T-even,
```

matching 035's cited census value exactly. ⚠ This is a **DERIVED** object, not a stored-tuple
duplicate (the 030 `X≡X` lesson): the `b_T` side is obtained by transforming `u_T`'s tuple through
the curl map, then checked against 035's census. Corrupting the curl-parity map (e.g. curl
preserves polar, or curl flips T-parity) makes the derived `b_T` parity diverge from the cited
census value and the tooth fires. Cites 035 `PARITY_ROTATION` + `PARITY_TIME_REVERSAL` — does NOT
re-run them.

## 2. The explicit Maxwell comparison — the departure (`MAXWELL_B_T_ODD_DEPARTURE`)

A physical Maxwell magnetic field `B` is time-reversal **ODD** and rotation **axial** — the
canonical external benchmark, stated as the typed pair

```text
Maxwell B = { time_reversal: −1, rotation: axial_vector }.
```

Maxwell `B` has **no** model-specific `R_w`/`P_w` benchmark (`R_w`/`P_w` are throat-coordinate
reflections native to THIS model, not properties of an external Maxwell field), so the comparison is
computed over EXACTLY the comparable domain `{time_reversal, rotation}`. The model's `b_T` is
time-reversal **EVEN** (`+1`, derived in §1); the two **DIFFER** in the time-reversal parity:

```text
departure_holds = ( b_T_T_parity ≠ maxwell_B_T_parity ) = ( +1 ≠ −1 ) = True.
```

⚠ **Both sides of the inequality are corruptible** (no can't-fail conjunct, the 035 lesson): the
`b_T` side is the DERIVED parity of §1 (corruptible via the curl map / a `u_T` flip); the Maxwell
benchmark side is corruptible in this tooth's OWN ablation (set the benchmark to T-even ⇒ the
comparison finds spurious agreement ⇒ `departure_holds` recomputes False and fires). The stage ships
**no** bare `b_T_T_parity == +1` literal assert — the departure is enforced as the **inequality
against the named Maxwell benchmark**.

## 3. The departure is localized to the time-reversal axis (`DEPARTURE_LOCALIZED_TO_T`)

`b_T`'s rotation parity is **axial** — which **AGREES** with a real Maxwell `B` (also axial). So the
ONLY parity axis on which `b_T` and Maxwell `B` disagree is the time-reversal axis:

```text
b_T_rotation == maxwell_B_rotation == axial_vector          (rotation agrees)
disagreement_set over {time_reversal, rotation}  ==  { time_reversal }   (only T mismatches).
```

The departure is therefore a genuine **T-parity** effect, not a rotation-parity artifact. Mutation:
make `b_T` polar (or make the Maxwell benchmark polar) ⇒ the rotation-agreement conjunct fails / the
disagreement-set gains `rotation` ⇒ fires.

## 4. Magnetism requires the active-drain time-arrow `τ_d` (`ACTIVE_DRAIN_TIME_ARROW_REQUIRED`)

The T-even parity of `b_T` is a **consequence** of the throat being an **active drain** with a
time-arrow `τ_d` (T-odd; 034 / R68). A **passive** T-even throat would not supply the `O(V)` moving
row at all. The census fixes the T-parity propagation chain, computed by composing the `±1`
T-parities:

```text
τ_d (T-odd, −1)  →  q_T = λ_T τ_d (T-odd, −1)  →  J_T = q_T s η V  (odd·even·odd = EVEN, +1)
                 →  u_T (T-even, +1)  →  b_T = ∇×u_T (T-even, +1).
```

So `b_T` inherits its T-even parity from the active drain's T-odd arrow through the source (`s`
T-even, `V` T-odd). The tooth composes the chain from the census T-parities and asserts it closes on
`b_T_T_parity = +1` **given `τ_d` T-odd** — NOT a stored `b_T = +1` literal. Mutation: make `τ_d`
T-even (a passive throat) ⇒ the chain re-derives `q_T` T-even ⇒ `J_T` T-odd ⇒ `u_T` T-odd ⇒ `b_T`
T-odd (and the `O(V)` row would not be supplied) ⇒ the computed chain no longer yields the T-even
`b_T` ⇒ fires. Cites 034 / R68 + 035 `PARITY_TIME_REVERSAL`.

## 5. The departure is internally self-consistent (`DEPARTURE_SELF_CONSISTENT`)

`b_T`'s T-even parity is reached by **two independent derivations** that must AGREE, and the
parity-DECIDED departure is orthogonal to the sector's open landing. The tooth is the conjunction of
two independently-computed objects:

- **(i) two-path agreement.** The **curl-inheritance** path (from `u_T`, §1) and the
  **active-drain source-propagation** path (`τ_d → … → b_T`, §4) both yield **T-even** — an internal
  cross-check analogous to 038's `adjudicate == oracle`. Mutation: corrupt one of the two paths
  (break the curl-parity map on path A while leaving the active-drain path B intact) ⇒ path A ≠ path
  B ⇒ the two-path cross-check recomputes False and fires.
- **(ii) orthogonality to the 038 R1 landing.** The departure's decided-axis set `{time_reversal}`
  (§3's disagreement set) is DISJOINT from the sector's cited R72 unresolved-axis set
  `{sign, magnitude}`, so the T-even departure neither asserts nor contradicts
  `R1_REQUIRED(electric_bc_selection)`: the parity is DECIDED while the magnitude/sign stay R1.
  Mutation: inject `time_reversal` into the cited R72 unresolved-axis set ⇒ the two sets intersect
  ⇒ the disjointness object recomputes False and fires.

⚠ This is a genuine computed two-path agreement + a computed set-disjointness (NOT a self-set
framing flag — the paths/sets are independently computed objects, avoiding 033's redundant-framing
smell); BOTH sub-objects are independently able-to-fail — each evidenced by its OWN computed
diagnostic (`path_only` / `axis_only`) — though they share ONE `DEPARTURE_SELF_CONSISTENT`
env-switch (the compound-tooth granularity recorded as non-blocking note (2) in Verification).

## 6. Departure typing (recorded honestly)

- The verdict is a **CHARACTERIZED-DEPARTURE, first-class** (`docs/model_map.md` §4). It is a
  load-bearing prediction, kept clean and simple, NEVER softened or rescued. Reads out as
  **"magnetism is NOT exact Maxwell."**
- It is **not a bug**, **not a knob and discharges no knob**; it is **not a reduction / codimension
  edge**; it **does not shrink the irreducible count**. Registered as edge **R73** (a departure edge
  that discharges nothing) plus **departure-support** facts (the `b_T` axial + T-even parity, the
  Maxwell-`B` T-odd benchmark, the active-drain T-parity propagation chain — dimensionless
  structural facts, no `[L, T, M]` to reduce).
- It is the **parity-census-level** departure only (the analog of 033's quadratic Stage-1 gate). The
  departure is DECIDED at the parity level while the sector *magnitude/sign* remains R1 (038 / R72;
  `q_T` R67; `A_E` R63) — the two live on orthogonal axes. 039 does **not** execute an all-orders
  dynamical non-Maxwell proof, does **not** re-adjudicate the 038 landing, and does **not** resolve
  the electric/magnetic sign or magnitude.

## 7. Source-to-stage predicate manifest

Completeness certificate: **no silently-dropped source claim.** 039 owns **NO** source-build tooth —
it REUSES 035's parity results and authors its own NEW `b_T` T-even departure assertions,
cited-not-owned (as `part5_magnetism_atomic_split.md` records: "cites V-2's `PARITY_TIME_REVERSAL` +
`PARITY_ROTATION`; authors its own `b_T`-axial + T-even vs Maxwell-T-odd asserts"). Every source
tooth of the source build's 35-tooth `TOOTH_ORDER` (`magnetism_moving_throat_check.{py,wl}`, commit
`53cf049f`) lands as **CITED** (upstream provenance, re-instantiated but not re-run), **SCOPED_OUT**
(owned by another Part-V stage, each named LITERALLY), or **build-global re-instantiated**. Same
partition in **both** engines, computed at runtime from the same canonical `(id, disposition, owner)`
triples:

```text
partition = { CITED: 2, SCOPED_OUT: 30, build-global: 3 }   (35 total, the SOURCE-MANIFEST arithmetic).
```

- **CITED (2):** `PARITY_TIME_REVERSAL` (⇒ `u_T` T-even), `PARITY_ROTATION` (⇒ `b_T` axial) — owned
  UPSTREAM by 035 (`b55a7a65`), re-instantiated here as the departure's input facts, NOT re-run.
- **SCOPED_OUT — owned upstream (30):**
  - **V-1 (034, `109070da`):** `FIELD_IDENTITY_UNITS`, `ACTION_KINETIC`, `ACTION_COUPLING`,
    `ACTION_STABILITY`, `G0_DAMAGE`, `LEDGER_READY_ROW` (6)
  - **V-2 (035, `b55a7a65`):** `SOURCE_TRANSLATION_CONTINUITY`, `SOURCE_NOT_IMPORTED`,
    `SOURCE_BASIS`, `PARITY_RW`, `PARITY_PW` (5) — the two `PARITY_*` in scope are CITED above
  - **V-3 (036, `df045a74`):** `BOOST_PROJECTOR`, `BOOST_GENERAL_VELOCITIES`, `BOOST_NEXT_ORDER` (3)
  - **V-4 (037, `c8780f00`):** `DIRECT_SOURCE`, `DIRECT_PROJECTOR`, `DIRECT_EXCHANGE_SIGN`,
    `DIRECT_FALLOFF`, `DIRECT_VELOCITY_ORDER`, `ROUTE_INDEPENDENCE`, `BOOST_COMMON_VELOCITY`,
    `COMPARE_COMPUTED`, `DELTA_RATIO`, `CONE_RATIO`, `QMAG_R1` (11)
  - **V-5 (038):** `TRUTH_TOTALITY`, `TRUTH_PRECEDENCE`, `LANDING_OWNERSHIP`, `ACTIVE_FLUX_CAVEAT`,
    `HOOK_LORENTZ` (5)
- **Build-global re-instantiated (3):** `TARGET_BLINDNESS`, `DUAL_ENGINE_TERMS`, `UNITS_RESTORED`.

Accounting: **2 CITED + 30 scoped-out + 3 build-global = 35** ✓. ⚠ This `35` (the SOURCE-MANIFEST
arithmetic over the source build's 35 teeth) is DISTINCT from the **10 EXECUTABLE** stage teeth: the
FIVE authored departure teeth (`B_T_AXIAL_T_EVEN`, `MAXWELL_B_T_ODD_DEPARTURE`,
`DEPARTURE_LOCALIZED_TO_T`, `ACTIVE_DRAIN_TIME_ARROW_REQUIRED`, `DEPARTURE_SELF_CONSISTENT`) + the
THREE build-global (`TARGET_BLINDNESS`, `DUAL_ENGINE_TERMS`, `UNITS_RESTORED`) +
`VERDICT_REDERIVATION` + `SOURCE_TO_STAGE_MANIFEST`, each with its OWN env-switch mutation firing at
its OWN assert in both engines. The `SOURCE_TO_STAGE_MANIFEST` tooth's own ablation (remove a
scoped-out row / mis-scope a CITED parity tooth) recomputes the partition-completeness object False
and fires.

## 8. Consumes / cites

- **Consumes NOTHING new numeric.** 039 re-instantiates 035's already-built parity entries but does
  not re-derive them; it does not touch any far-field kernel, the source law, or the 1152-cell truth
  table.
- **Cites (provenance, NOT re-derived, NOT re-adjudicated, NOT re-counted; de-dup deferred to Part
  VII):**
  - **035's parity census (V-2, R69):** `u_T` T-even (`PARITY_TIME_REVERSAL`), `b_T` axial
    (`PARITY_ROTATION`), the full 24-cell census tuples — the departure's input facts (`b_T`'s
    axial + T-even parity was recorded RAW at 035; 039 supplies ONLY its departure interpretation).
  - **034's `τ_d` (V-1, R68) + `q_T` (R67) + the field identity `b_T = ∇×u_T`** — the active-drain
    time-arrow the T-even departure requires, and the identity for dims.
  - **038's sealed terminal landing (V-5, R72):** `R1_REQUIRED(electric_bc_selection)` — the
    sector's terminal verdict, of which the T-even departure is a HOOK (decided at the parity level
    while the sector landing stays R1 — orthogonal axes). 039 does NOT re-adjudicate the 1152-cell
    truth table and does NOT re-emit `MAGNETISM_DEPARTURE_CHARACTERIZED`.
  - **The canonical Maxwell-`B` benchmark** `{time_reversal: −1, rotation: axial_vector}` — a
    standard external reference (like a `benchmarks.yaml` literature value), stated explicitly, not
    derived from the model.
  - **The sibling departures:** the charge-sector `NATIVE_P_NO_EMERGENT_GAUSS` (033 / R66) and the
    light-sector `FAIL_CAUCHY_STRAY_LONGITUDINAL` (Part III / stage 003 / pathA_36) — the three
    first-class "EM is NOT exact Maxwell" departures in `docs/model_map.md` §4.

## Verification

- **Dual-engine, both exit 0, 10 executable teeth each, genuinely independent implementations.**
  `scripts/ledger_stage039_b_t_time_reversal_even_departure_sympy_audit.py` — **SymPy 10 teeth**.
  `mathematica/ledger_stage039_b_t_time_reversal_even_departure_mathematica_audit.wl` —
  **Mathematica 10 teeth** (`B_T_AXIAL_T_EVEN`, `MAXWELL_B_T_ODD_DEPARTURE`,
  `DEPARTURE_LOCALIZED_TO_T`, `ACTIVE_DRAIN_TIME_ARROW_REQUIRED`, `DEPARTURE_SELF_CONSISTENT`,
  `TARGET_BLINDNESS`, `DUAL_ENGINE_TERMS`, `UNITS_RESTORED`, `VERDICT_REDERIVATION`,
  `SOURCE_TO_STAGE_MANIFEST`). Standalone, print-only, assert-zero (`raise SystemExit(1)` /
  `Exit[1]`), no argparse harness, no JSON/YAML payload, **zero file-I/O between engines**;
  float-/machine-real-free integer-parity payload throughout. Each engine re-instantiates the census
  entries, derives `b_T`'s parity through the curl map, compares against the Maxwell benchmark,
  propagates the active-drain chain, and reaches `B_TIME_REVERSAL_EVEN` on its own — cross-engine
  agreement is that they independently produce the SAME derived `b_T` parity + `departure_holds` +
  verdict token (arbiter-confirmed by the orchestrator re-running BOTH engines, NOT an in-script
  compare pass). **No dual-engine disagreement.**
- **The `.wl` is a genuinely INDEPENDENT implementation of the curl-parity derivation.** ⚠ As
  flagged for the bookend, 039's core objects (`b_T` axial + T-even; Maxwell `B` T-odd + axial) are
  **fixed structural integers decided at 035** — there is limited freedom to "re-derive a different
  parity," so the two engines must agree on the SAME parity for every object (that agreement IS the
  cross-engine check), and independence is assessed on the derivation **organization**, not on
  inventing a different parity (fixed-parity-integer independence — the parity RESULT is correctly
  identical). The genuine independence lives in materially distinct organizations: the **SymPy**
  route applies the curl-parity map as an explicit per-axis **operator dictionary** and multiplies
  `±1` T-parities along the chain; the **Mathematica** route is materially distinct — **named-key
  `ReplaceAll` curl rules** + a **`Cross`/`Det` improper-rotation classifier** (using
  `(Qa)×(Qb) = det(Q) · Q(a×b)` to classify polar vs axial) + **`FoldList` parity propagation** for
  the active-drain chain. It is NOT a line-by-line port of the stage `.py` NOR of the source `.wl`'s
  index-extraction parity checks (`…check.wl:124` stores `PARITIES` as an `Association` of
  index-tuples and checks `PARITY_ROTATION`/`PARITY_TIME_REVERSAL` by literal index extraction — a
  correspondence to that would be a FAIL); the `.wl` carries its OWN native per-tooth ablations. This
  fixed-parity-logic independence-nuance was FLAGGED for + CONFIRMED at the Codex→Grok→Codex bookend.
- **Source-to-stage manifest.** Same partition both engines (`{CITED: 2, SCOPED_OUT: 30,
  build-global: 3}`, 35 total), the 30 scoped-out teeth named verbatim by owner stage (no wildcard
  families); the `SOURCE_TO_STAGE_MANIFEST` mutation removes a scoped-out row / mis-scopes a CITED
  parity tooth ⇒ the partition-completeness object recomputes False and fires.
- **Per-tooth ablation** (env switch `LEDGER_STAGE039_MUTATION`): all **10 teeth's mutations
  FIRED_AT_OWN_ASSERT per engine (20 mutation runs across the two engines).** ⛔ **NO can't-fail
  conjuncts** (the 035 lesson): `B_T_AXIAL_T_EVEN` derives `b_T` through the curl map (NOT a stored
  `b_T` literal); `MAXWELL_B_T_ODD_DEPARTURE` enforces the inequality against the named Maxwell
  benchmark with BOTH sides corruptible; `ACTIVE_DRAIN_TIME_ARROW_REQUIRED` composes the `±1`
  T-parity chain (NOT a stored `b_T` literal); `DEPARTURE_SELF_CONSISTENT` is a genuine
  two-independent-path cross-check + a computed set-disjointness (NOT a self-set framing flag). Each
  departure conjunct corrupts a derived object — the departure inequality is corruptible on BOTH
  sides (the `b_T` parity OR the Maxwell record).
- **The verdict tooth is non-tautological** — it mutates a **COMPUTED** parity object, NOT the final
  token (flipping the `B_TIME_REVERSAL_EVEN` literal tests nothing, the 030 `X≡X` lesson). Under
  `VERDICT_REDERIVATION` the pipeline CONSTRUCTS each mutated parity witness, re-derives `b_T`'s
  T-parity through the curl map / propagation chain, re-runs the Maxwell comparison, and asserts it
  **re-derives the named verdict**:

  | Computed input state (mutated parity) | Re-derived verdict |
  |---|---|
  | **production** (`u_T` T-even, `τ_d` T-odd) ⇒ `b_T` T-even, axial | **`B_TIME_REVERSAL_EVEN`** (the departure — `b_T` T-even ≠ Maxwell `B` T-odd; scope **CHARACTERIZED-DEPARTURE**, `\StatusOpen`) — the build-faithful headline token (named in `result.md` Appendix D + the ratified V-6 token), computed here from the census facts (NOT a source enum landing) |
  | flip the derived `u_T` T-parity (T-even → **T-odd**) ⇒ `b_T = curl(u_T)` re-derives **T-odd** | **`COUNTERFACTUAL_B_TIME_REVERSAL_ODD_MAXWELL_CONSISTENT`** — an AUTHORED stage-local counterfactual token (a hypothetical "if `b_T` were T-odd it would be Maxwell-consistent"; NOT a physical claim, NOT a source enum) |
  | flip **`τ_d`** to T-even (passive throat) ⇒ the propagation chain re-derives `b_T` **T-odd** (root-cause mutation, distinct computed input from the `u_T` flip) | **`COUNTERFACTUAL_B_TIME_REVERSAL_ODD_MAXWELL_CONSISTENT`** (same authored counterfactual via the active-drain root; PHYSICALLY a passive T-even throat supplies NO `O(V)` row, so this witness is a formal hypothetical, not a physical branch) |

  The `COUNTERFACTUAL_` prefix is load-bearing: it is emitted **ONLY by the verdict-tooth mutation
  witnesses**, NEVER as a physical model claim (the model only ever produces `b_T` T-even). Both the
  verdict-token computation and the authored counterfactual alternate were FLAGGED for + CONFIRMED at
  the Codex→Grok→Codex bookend (the bookend confirmed the census facts genuinely re-derive `b_T`
  T-even + axial, the Maxwell benchmark is correctly T-odd + axial so the departure is exactly in the
  T axis, and the counterfactual token name reads as a hypothetical). Same derived parity + verdict in
  BOTH engines.
- **`UNITS_RESTORED` — the dimensional firewall.** 039 derives no new dimensionful physics (parities
  are dimensionless structural facts); the whole-stage dimensional check operates on the **cited**
  field identity `b_T = ∇×u_T`: with `[u_T] = L` and `[∇] = L⁻¹` the relation forces `[b_T] = 1`
  (dimensionless), consistent with 034/035's `[b_T]=1`. Able-to-fail + free-carrier-independent:
  corrupting one cited dim (`[u_T]=L → 1`, or `[b_T]=1 → L⁻¹`) makes the `b_T = ∇×u_T` homogeneity
  object inhomogeneous and the units tooth fires. `TARGET_BLINDNESS` is the STRUCTURAL guard (no
  filesystem text-grep) that the cited magnetism symbols `{u_T, b_T, τ_d, q_T, s, V, J_T}` remain
  disjoint from the barred pathA_39 markers `{N_u, a_T, a′_T, a_L, q_A^T, q_L}` — the `.py` a
  structural symbol / live-object-graph inventory (a walk of the live SymPy-symbol / dataclass /
  mapping / collection objects — no `__code__`/`CodeType` inventory), the `.wl` a held-symbol /
  `DownValues` inventory (NOT an
  `Import[$InputFileName,"Text"]` scan); injecting a barred marker fires it. `DUAL_ENGINE_TERMS` is a
  LOCAL per-engine tooth (each engine validates its OWN computed-parity inventory — the derived `b_T`
  tuple, the Maxwell benchmark record, `departure_holds`, the T-axis disagreement set, the
  active-drain chain result, and the verdict token — and drops/corrupts a term to fire its OWN
  assert; it does NOT read the other engine).
- **Tri-review outcome (falsification-first — recorded transparently). ZERO remediation** (the
  hardened build applied the 035/036/037/038 lessons at authoring time, so there were no can't-fail
  conjuncts to remove).
  - **FIDELITY: CLEAN** — the curl-parity map genuinely derives `b_T = (−1, −1, +1, axial_vector)`
    from `u_T = (−1, −1, +1, polar_vector)`; the Maxwell benchmark is correctly
    `{time_reversal: −1, rotation: axial_vector}`; the disagreement set is exactly `{time_reversal}`;
    the active-drain chain `τ_d → q_T → J_T → u_T → b_T` composes to `b_T` T-even; the two-path
    cross-check agrees; the verdict-tooth re-derivation targets are correctly named (both `u_T`-flip
    and `τ_d`-flip re-derive `COUNTERFACTUAL_B_TIME_REVERSAL_ODD_MAXWELL_CONSISTENT`).
  - **ADVERSARIAL: CLEAN** — every one of the 10 teeth fires at its own assert; no vacuous / X≡X /
    tautological construct; the departure inequality is corruptible on both sides; the `.wl` confirmed
    genuinely independent of BOTH the stage `.py` and the source `.wl` (no line-by-line
    index-extraction correspondence).
  - **Documented non-blocking notes (verified-safe, NOT remediated):**
    1. `MAXWELL_B_T_ODD_DEPARTURE` bundles its `b_T`-side and Maxwell-benchmark-side sub-conjuncts
       under ONE env-switch (a minor granularity choice vs the directive's "ablate separately") — but
       EACH conjunct is independently able-to-fail, evidenced by the computed `b_side` / `maxwell_side`
       diagnostics.
    2. `DEPARTURE_SELF_CONSISTENT` likewise bundles its two-path-agreement and set-disjointness
       sub-objects under ONE env-switch — but EACH is independently able-to-fail, evidenced by the
       computed `path_only` / `axis_only` diagnostics.
    3. `DEPARTURE_SELF_CONSISTENT` conjunct-2 (`decided_axes = {time_reversal}` disjoint from R72's
       `{sign, magnitude}`) is always-true in production (the disagreement domain is only
       `{time_reversal, rotation}`, which never contains `sign`/`magnitude`) — but it is able-to-fail
       via corrupting the cited R72 benchmark set (injecting `time_reversal`), a legitimate
       benchmark-side corruption exactly like the Maxwell-benchmark side of tooth 2.

  Arbiter re-runs of both engines reproduce the derived parities, the departure, and the verdict token
  `B_TIME_REVERSAL_EVEN`; the tri-review leaves the departure decisive.

## Downstream consumers

- **⭐⭐ Part V is COMPLETE** after 039: 034 (action row, R67/R68) → 035 (native source + census,
  R69) → 036/037 (the two blind routes + boost-consistency, R70/R71) → 038 (sealed terminal R1
  landing, R72) → **039 (the first-class characterized departure `B_TIME_REVERSAL_EVEN`, R73)**. The
  remaining EM-sector work (the `c_E=c_γ` cone lock re-adjudication in Part VI `pathA_40`; the shared
  sim-deferred throat solve that would convert both the electric and magnetic R1 → resolved; the
  `F_e/F_g` hierarchy capstone; the G0-vs-Part-I de-dup) belongs to Parts VI / VII, NOT to any
  deferred-from-039 scope.
- **Parameter register:** edge **R73** (the `B_TIME_REVERSAL_EVEN` characterized departure —
  discharges NO knob, NOT a reduction, does NOT shrink the irreducible count) + the departure-support
  facts (the `b_T` axial + T-even parity, the Maxwell-`B` T-odd benchmark, the active-drain T-parity
  propagation chain — structural, dimensionless). Cross-referenced to the sibling departures R66
  (`NATIVE_P_NO_EMERGENT_GAUSS`, charge) and stage003 `FAIL_CAUCHY_STRAY_LONGITUDINAL` (light).
- **`docs/model_map.md` §4:** the `B_TIME_REVERSAL_EVEN` departure bullet (already present) — kept
  consistent: `b_T = ∇×u_T` is time-reversal EVEN (Maxwell `B` is T-odd), correctly axial; a concrete
  self-consistent prediction; magnetism requires the active-drain time-arrow `τ_d`.
- **Part VII:** enters the honest departure ledger alongside the charge-sector twin
  `NATIVE_P_NO_EMERGENT_GAUSS` and the light-sector twin `FAIL_CAUCHY_STRAY_LONGITUDINAL` — a
  first-class recorded departure, not a reduction credit.

## Provenance

- **Physics source:** `software/em_charge_attribute/magnetism_moving_throat_result.md` (the
  **Appendix D "Hooks"** `B_TIME_REVERSAL_EVEN` characterized-departure bullet + the **Q-CURRENT
  parity census** table — the `b_T`/`u_T`/`τ_d` rows + the `P_w`-reflection caveat) +
  `software/em_charge_attribute/magnetism_moving_throat_check.py` (the parity machinery: `PARITIES`
  `OrderedDict`, the tuple `b_T = (−1, −1, +1, "axial_vector")`; the two CITED parity teeth
  `PARITY_ROTATION` + `PARITY_TIME_REVERSAL`) +
  `software/em_charge_attribute/magnetism_moving_throat_check.wl` (the mirror `PARITIES` +
  `PARITY_ROTATION`/`PARITY_TIME_REVERSAL`), commit `53cf049f`. **`check.py`/`check.wl` are
  AUTHORITATIVE over the `result.md` prose for the parity integers; the departure *characterization*
  is the `result.md` Appendix-D bullet + the ratified V-6 token in
  `notes/part5_magnetism_atomic_split.md`.** Stage 039 owns NO source parity tooth — it CITES 035's
  already-built `PARITY_TIME_REVERSAL` (`u_T` T-even) + `PARITY_ROTATION` (`b_T` axial) and authors
  the NEW departure characterization (the Maxwell comparison, the T-axis localization, the
  active-drain requirement, the internal self-consistency). The source `argparse`/`--json-only`
  harness, the `MAGNETISM_PAYLOAD_MUTATION`/`MAGNETISM_JSON_ONLY` plumbing, the source-law/continuity
  machinery (035), the Route-A/Route-B far-field internals (036/037), the sealed §4 truth-table /
  adjudication cluster (038 — cited as the terminal landing, NOT re-adjudicated), and all JSON/log
  file writing are STRIPPED (print-only / zero-file-I/O / independent-tokens contract). The stage
  `.wl` is a materially distinct Wolfram implementation (named-key `ReplaceAll` curl rules +
  `Cross`/`Det` improper-rotation classifier + `FoldList` parity propagation), not a line-port of the
  `.py` and not the source `.wl`'s index-extraction parity checks.
- **Consumes:** nothing new numeric — imports by citation 035's parity census (`u_T` T-even, `b_T`
  axial), 034's `τ_d` (R68), 038's terminal landing (R72), and the canonical Maxwell benchmark.
- **Cites (provenance, NOT re-derived, NOT re-adjudicated, NOT re-counted; de-dup deferred to Part
  VII):** `u_T` T-even + `b_T` axial (035, R69); `τ_d` (034, R68) + `q_T` (R67) + `b_T = ∇×u_T`;
  `R1_REQUIRED(electric_bc_selection)` (038, R72); the Maxwell-`B` benchmark
  `{time_reversal: −1, rotation: axial_vector}`; the sibling departures
  `NATIVE_P_NO_EMERGENT_GAUSS` (033 / R66) + `FAIL_CAUCHY_STRAY_LONGITUDINAL` (Part III / stage003 /
  pathA_36).
- **Governing:** `notes/ledger_v2_blueprint.md` §5 (reshape spec) + §6 (per-tooth ablation);
  `notes/part5_magnetism_atomic_split.md` (V-6 = the `b_T` T-even departure; the tooth-allocation
  table = 5 authored departure teeth + build-global `TARGET_BLINDNESS`/`DUAL_ENGINE_TERMS`/
  `UNITS_RESTORED`); `docs/model_map.md` §3.5 + §4 (the honest departure ledger — the magnetic twin
  of `NATIVE_P_NO_EMERGENT_GAUSS` and the sibling of `FAIL_CAUCHY_STRAY_LONGITUDINAL`). Reshape
  directive + review trail:
  `research/pde_ledger_v2/_scratch/stage039_reshape_directive.md`.
