# Directive pathA_41 — NG5 `SECOND_MEDIUM_DRIFT`: does the four-sector model live on ONE medium, or several? (consistency-knit gate 2)

**Status:** DRAFT v2 (Codex design-review + GLM-5.2 tertiary both folded — see §10; `SOUND` no-BLOCKERs from GLM; pending Codex
confirm-to-green → user gate before execution). **⚠ USER ESCALATION (GLM F5):** the pathA_35 committed `SECOND_MEDIUM_DRIFT_AT_FREEZE(11)`
count is likely INFLATED (a derived functional `varrho_br[ρ]` was re-counted as the postulate `ρ_br`) — GLM recommends an ERRATUM NOTE on
the pathA_35 G0 report (not a re-freeze, not amending earned specs). Decision pending user. **Scope provenance:** `_scratch/consistency_knit_ng5_scope_v3.md` (`SOUND_AS_A_SCOPE`; Codex×2 + GLM-5.2 tertiary
folded — the `varrho_br[ρ]` lineage BLOCKER + 3 more). **Conceptual source:** `docs/conceptual_foundation.md`; `STATUS.md` ▶ RESUME HERE.
**Prior gate:** `pathA_40` cone-lock = `CONE_LOCK_CALIBRATED` (left `FREEDOM_UNCONSTRAINED{C_hu, ρ_br}` for THIS gate to resolve).

**This is consistency-knit gate 2** (after the cone-lock). It tests the program's central claim — "all four force-sectors come from ONE
brane+bulk medium" — by classifying every independent parameter as reducible-to-the-one-medium or irreducibly-separate. `pde_ledger`
assembly is gate 3 (separate).

---

## 0. What this gate is + honesty preregistration

Unlike the cone-lock (near-vacuous by construction), **NG5 is where the knit's real falsification power lives.** A raw new-input count is
useless (trivially many); the decisive test is **REDUCIBILITY**: is each independent parameter expressible in the one shared medium, or
irreducibly its own thing? An irreducibly-independent active parameter (or a sector-contradiction) is a **first-class falsification of the
one-medium claim** — NOT a retraction of any sector's internally-earned spec. **`SECOND_MEDIUM_DRIFT` is knit-only, non-retroactive:** each
sector (light `c_γ`, charge localization, magnetism force, gravity bundle) remains earned; DRIFT says only that they cannot be seated on
one shared parameter set.

**Honest prior (pre-registered):** likely `ONE_MEDIUM_CONDITIONAL` — conditional on two named registered deferred solves {B4 brane layer
profile, pathA_40 Route-A throat}. `ρ_br` reduces via `varrho_br[ρ]=∫_layer dn·m·ρ` (pathA_25, medium-derived), and `c_E` via Route-A.
BUT a genuine `SECOND_MEDIUM_DRIFT` is reachable (if the `varrho_br/ρ_br` same-object reduction is rejected, if the B4 profile itself
carries independent active inputs = "drift one level down," if a Family-L/R/C input is active-and-irreducible, or if a route leaves a
residual free factor). The pathA_25→pathA_35 lineage downgrade (a derived functional re-counted as a postulate) is recorded regardless of
the headline verdict (on the SAME-object branch, §3.1). **Do not launder an unearned pass through "sim-deferred"** — and do NOT hardcode
the reducible landing: the DIFFERENT-object and `B4_UNSOLVED` branches must stay genuinely reachable (Codex design-review BLOCKERs).

---

## 1. The one-medium BASE SUBSTRATE + two axes (scope §1)

**Base substrate (declared, target-blind):** GNLS `{ρ, K, m, a}` (`c_s=c_s(ρ)=√(5Kρ⁴/m)` is a FUNCTION of ρ, not a constant) + the free
brane-geometry width `ℓ_g` (state whether `ℓ_g` is accepted substrate, gets its own origin row, or reduces to `a` — not currently earned).
The TWO profile objects are NOT one "shared geometry": (i) codim-1 confinement `g_ℓ(w)` (pathA_35, width `ℓ_g`); (ii) B4 smectic layer
profile `Σ_n[ρ]`/measure `δ_Σ[ρ]` (pathA_25). Each gets a §2 row.

**Two axes:** `incidence(p)` (sectors using it — a reuse tag, shared counts once but does NOT auto-pass) SEPARATE from `origin(p)` ∈
{`REDUCIBLE_DERIVED`, `REDUCIBLE_IN_PRINCIPLE_SIM_DEFERRED`, `IRREDUCIBLY_INDEPENDENT`}, with `CONTRADICTORY_VALUES` a verdict-level no-go.

---

## 2. The origin ledger (build a row for EVERY parameter)

Row schema: `p | incidence | dims | canonical source(:line) | local status (input/import/derived) | earned equalities | candidate
registered deferred route | pathA_40 freedom relevance | origin verdict`. Required coverage (none omitted):
- **Brane-inertia lineage pair (the crux):** `varrho_br[ρ]` (pathA_25 `∫_layer dn m ρ`, medium-derived, 0 free) AND `ρ_br` (pathA_35
  postulated) — adjudicated in §3.1.
- **Continuum candidates:** `μ_R`(=`μ_br`), `ρ_B0, χ_c, Q_E, ℓ, b, C_hu` + derived composites `B_eff, K_h, q_h, M_h, c_E, c_γ` (each
  states derived-FROM-what; derived-inside-a-sector ≠ REDUCIBLE_DERIVED to the base substrate).
- **B4 layer objects:** `Σ_n[ρ]`, `δ_Σ[ρ]` — tagged as the registered deferred solve (§3.3), not silent substrate.
- **pathA_25 driver inputs:** `c_L1, c_L2, A_R, k_R, λ_Cdiv, χ_Cpin, J_Pu, κ_Pu` — for EACH, run the active-path test (§3.4): active in a
  four-sector path → candidate `IRREDUCIBLY_INDEPENDENT`; else explicit `OUT_OF_ACTIVE_NG5`.
- **Frozen/internal/source:** `λ_Pu, Ω_w` (`FROZEN_UNUSED_OUT_OF_ACTIVE_NG5`); `λ_N, λ_τ` (pathA_38 wall-internal/exclude); `Ν, a_T, a_Tp,
  a_L` (`SOURCE_NORMALIZATION`); `ℓ_g` vs `ℓ` (do not merge). **These status-tagged rows are outside the active-medium origin competition**
  (they are ledger-eligibility tags, per Codex) — but each MUST appear.

**No `ONE_MEDIUM_CONDITIONAL` may be earned by adjudicating only the hinges while passing the rest as shared.**

---

## 3. The decisive computations (the falsifiable core)

### 3.1 The `varrho_br[ρ]` vs `ρ_br` same-object adjudication (dual-engine — the crux)
COMPUTE whether pathA_25 `varrho_br[ρ]=∫_layer dn·m·ρ` and pathA_35 `ρ_br` denote the SAME functional from the frozen DEFINITIONS (the
symbolic functional forms — the profiles `Σ_n[ρ]` themselves are unbuilt, GLM F3): same
kinetic role (`½·(·)·(D_t u)²`), same dimension (M L⁻³ — verify the `∫_layer dn` layer-normal factor supplies the L that fixes bulk
`m·ρ` (M L⁻⁴) → surface M L⁻³), and whether `ρ_br` carries any residual free multiplier beyond `varrho_br`. Outcomes:
- **SAME object** → pathA_35 re-postulated as independent what pathA_25 derived → `ρ_br` origin inherits `varrho_br`'s (medium-derived,
  conditional on B4, §3.3); RECORD the lineage downgrade as `OVERCOUNT_OR_SMUGGLE` (an overcount of the pathA_35 G0 `SECOND_MEDIUM_DRIFT_AT_FREEZE(11)`
  ledger, OR an in-lineage smuggle) — a first-class finding regardless of headline.
- **DIFFERENT objects** → two brane inertia densities (one derived, one postulated); each keeps its own row/origin (pathA_35 `ρ_br` per
  the §3.2 DIFFERENT branch). RECORD `DIFFERENT_BRANE_INERTIA_OBJECTS` (NOT `OVERCOUNT_OR_SMUGGLE`, which is SAME-branch only).

The `OVERCOUNT_OR_SMUGGLE` finding is emitted ONLY on the SAME-object branch; it is NOT an unconditional rider (Codex BLOCKER).

### 3.2 The `ρ_br` origin — CONDITIONAL on §3.1's branch (Codex BLOCKER: do not hardcode)
The `ρ_br` origin is computed FROM the §3.1 outcome — it is NOT unconditionally "reducible":
- **§3.1 = SAME object** → `ρ_br` inherits `varrho_br`'s origin: `REDUCIBLE_IN_PRINCIPLE_SIM_DEFERRED` via B4 (§3.3), OR `REDUCIBLE_DERIVED`
  if NG5 closes the layer-mass integral on an accepted profile (→ **mandatory pathA_40 rerun**, §4).
- **§3.1 = DIFFERENT objects** → pathA_35 `ρ_br` keeps its OWN origin row and is `IRREDUCIBLY_INDEPENDENT` unless a registered route is
  found for it specifically → contributes to `SECOND_MEDIUM_DRIFT`. (This branch MUST stay reachable — it is the §3.1 falsifier.)

### 3.3 The B4 registered-deferred anchor + the "drift one level down" check — with an explicit UNKNOWN landing
Name **B4 (`Σ_n[ρ]`/`δ_Σ[ρ]`, pathA_25)** as the registered deferred solve `varrho_br` (and the layer measure) pass through. **The deepest
check (scope §9):** does the B4 layer profile ITSELF require independent active inputs (from the Family-L/R/C driver set, §3.4)? **This is
NOT fully computable from the pathA_25 freeze — the B4 solve is UNBUILT** (Codex BLOCKER). Therefore the PRODUCTION landing has THREE
outcomes:
- `B4_DETERMINED_BY_MEDIUM` (only if an earned B4 relation shows `Σ_n[ρ]` is fixed by medium/throat data alone) → `ρ_br` genuinely
  `SIM_DEFERRED`-reducible (or `DERIVED`).
- `B4_NEEDS_INDEPENDENT_INPUTS` (only if an earned B4 relation shows `Σ_n[ρ]` requires independent active coefficients) → drift **moves one
  level down** to those inputs (`IRREDUCIBLY_INDEPENDENT` → `SECOND_MEDIUM_DRIFT`).
- **`B4_UNSOLVED_UNKNOWN` (the EXPECTED production case — B4 is unbuilt) → `REDUCIBLE_IN_PRINCIPLE_SIM_DEFERRED`** (registered route named,
  outcome pending the solve). **Until B4 is solved, `ρ_br` (and via Route-A, `c_E`) stay `SIM_DEFERRED`, pathA_40 stays conditional, NO
  rerun.** Production does NOT emit a drift verdict from the unresolved-B4 case — a drift-one-level-down verdict is only reachable via the
  §6 source-level control, never from unbuilt-B4 uncertainty.

**Binding forward re-open trigger (GLM F2 — so `CONDITIONAL` cannot become a permanent laundered pass):** the `ONE_MEDIUM_CONDITIONAL`
landing carries an EXPIRY. **When B4 (or Route-A) is solved in ANY future gate, NG5 MUST re-open:** re-decide `B4_DETERMINED_BY_MEDIUM` vs
`B4_NEEDS_INDEPENDENT_INPUTS` (and the analogous Route-A outcome), and the affected `SIM_DEFERRED` rows RE-ENTER the origin competition
(§5) — which can then flip the verdict to `SECOND_MEDIUM_DRIFT`/`NO_GO` or confirm `ONE_MEDIUM_CONSISTENT`. Record this trigger in
`_results.yaml` alongside each sim-deferred row (row → registered solve → "re-adjudicate on solve"), mirroring the pathA_40 rerun-trigger.

### 3.4 Family-L/R/C active-path test — separate DECIDABLE branch-status from UNDECIDABLE knit-propagation (Codex BLOCKER)
For each of `c_L1, c_L2, A_R, k_R, λ_Cdiv, χ_Cpin, J_Pu, κ_Pu`, split the question:
- **Branch status (decidable from pathA_25):** is it an active-baseline input (`c_L1, c_L2, J_Pu, κ_Pu`) vs a branch/sensitivity input
  (`A_R, k_R, λ_Cdiv, χ_Cpin`)? Record from the reports.
- **Four-sector-knit propagation (may be UNDECIDABLE now):** does it feed an active four-sector path via `Σ_n[ρ]` (§3.3) or the gravity
  bundle (§8)? If provably active-and-irreducible → `SECOND_MEDIUM_DRIFT`. If the propagation runs through the unresolved B4/gravity solve
  → `SIM_DEFERRED` (do NOT turn unresolved-B4 propagation into a drift verdict); if provably not on any active path → `OUT_OF_ACTIVE_NG5`.

**Production-collapse transparency (GLM E1):** state plainly that in the EXPECTED production case — where the only active paths for these
drivers run through the unresolved B4/gravity solves — ALL 8 Family-L/R/C drivers land `SIM_DEFERRED` or `OUT_OF_ACTIVE_NG5`, and NONE
contributes a production drift signal. That is honest (the drift they could carry is gated on the same unbuilt B4/gravity solves, re-opened
by the §3.3 forward trigger), NOT a silent pass — the narrow production drift-surface is disclosed, not hidden.

### 3.5 `c_E` (unchanged from cone-lock) + `C_hu`
`c_E`: `NOT_DERIVED_NOW` (pathA_40); `SIM_DEFERRED` iff a target-blind registered route (Route-A `{R1..R5}`) determines both sides of
`c_E²=K_h/M_h=μ_R/ρ_br=c_γ²` from medium/throat data without inserting the cone equality (target-blind = declared outputs from
medium/throat data alone, no sector-specific measured/tuned lock); else `IRREDUCIBLY_INDEPENDENT`. `C_hu`: sim-deferred (throat overlap,
registered route) vs free knob; **may NOT be BOTH sim-deferred AND freedom-discharged** — the pathA_40 `freedom_tie` (C_hu = charge
residue) VIOLATES stability under `L_B` → `FREEDOM_CONTRADICTED`/NO_GO route.

---

## 4. pathA_40 freedom state-machine (the concrete hand-off + rerun triggers)
Resolve each of `{C_hu, ρ_br}` into exactly one, with the rerun consequence:
- `FREEDOM_CERTIFIED_CURRENT_LEDGER{p}` — no earned tie today (no rerun; future ties trigger one).
- `FREEDOM_REVOKED_DERIVED{p, relation}` — NG5 earns a relation → **immediate pathA_40 rerun** (may change the Δr count / flip to NO_GO;
  reuse the pathA_40 lock+stability machinery to check SAT under the new relation).
- `FREEDOM_SIM_DEFERRED{p, route}` — no current relation but a named registered solve may tie it → keep pathA_40 conditional, record the
  trigger. **Expected `ρ_br → FREEDOM_SIM_DEFERRED{B4}`** (or `REVOKED_DERIVED` if §3.2 closes the integral).
- `FREEDOM_CONTRADICTED{p}` — a derived tie conflicts with stability/lock → no-go/rerun (the `C_hu` charge-residue tie is the live case).

---

## 5. Verdict grammar (first-match; DRIFT is knit-only)
1. Any `NO_GO` (`CONTRADICTORY_VALUES`, or a pathA_40 rerun that flips) → `NO_GO(sector-contradiction)` / `NO_GO(cone-lock-feedback)` (+ the
   conflicting relation / unsat core).
2. Else any ACTIVE medium parameter `IRREDUCIBLY_INDEPENDENT` → `SECOND_MEDIUM_DRIFT(irreducible={…})` (knit-only; sectors stay earned).
3. Else any `REDUCIBLE_IN_PRINCIPLE_SIM_DEFERRED` → `ONE_MEDIUM_CONDITIONAL(sim-deferred={…})` — where each sim-deferred row LISTS its own
   attached registered solve (e.g. `ρ_br→B4`, `c_E→Route-A`); do NOT hardcode `{B4, Route-A}` for the whole verdict (Codex BLOCKER — a
   case with only one pending, or an unexpected route, must be reported per-row).
4. Else all `REDUCIBLE_DERIVED`/accepted substrate → `ONE_MEDIUM_CONSISTENT` (too-clean → extra scrutiny).
Every verdict carries: the §3.1 lineage finding (`OVERCOUNT_OR_SMUGGLE` on the SAME-object branch OR `DIFFERENT_BRANE_INERTIA_OBJECTS` on
the DIFFERENT branch — whichever §3.1 computed, NOT an unconditional overcount rider), the pathA_40 freedom states + rerun triggers (§4),
and the drift count (active `IRREDUCIBLY_INDEPENDENT` count).

---

## 6. Controls (able-to-fail; each states its SOURCE MUTATION + RECOMPUTED INVARIANT + EXPECTED TRANSITION; no direct verdict/flag inputs)
Not every control flips the HEADLINE verdict — some flip a per-parameter ORIGIN, a lineage FINDING, or a freedom STATE (labelled below);
that is correct, and the directive does not overclaim a headline flip (Codex BLOCKER). All controls mutate SOURCE facts and RECOMPUTE.
- **Same-object → DIFFERENT (headline-relevant):** mutate the source so the two inertia definitions have genuinely different dimension/role
  (e.g. a M L⁻⁴ bulk coefficient vs the M L⁻³ surface one) → §3.1 recomputes DIFFERENT → §3.2 gives pathA_35 `ρ_br` its own row → if no
  route, `IRREDUCIBLY_INDEPENDENT` → **headline flips to `SECOND_MEDIUM_DRIFT`.** (Proves §3.2 is not hardcoded reducible.)
- **Overcount FINDING (de-rigged — recompute, do NOT re-tag):** mutate source facts to carry {same role, same M L⁻³ dimension, no residual
  multiplier, tagged derived in the earlier gate, tagged postulated in the later gate} → the pipeline MUST RECOMPUTE `OVERCOUNT_OR_SMUGGLE`
  from those five invariants (NOT from a fed classification). Flips the lineage FINDING, not necessarily the headline.
- **Drift-one-level-down (headline; SOURCE-level, distinct from unbuilt-B4):** mutate the B4 profile source so `Σ_n[ρ]` explicitly depends
  on a synthetic independent coefficient → the pipeline COMPUTES that no registered route covers that coefficient (the "no registered
  route" is pipeline-DERIVED from the route inventory, NOT asserted by the control fixture — GLM F4) → §3.3 recomputes
  `B4_NEEDS_INDEPENDENT_INPUTS` → drift routes to that coefficient → **headline `SECOND_MEDIUM_DRIFT`.** (The production `B4_UNSOLVED_UNKNOWN`
  case is NOT this control and must NOT drift.)
- **Irreducible (headline):** a synthetic active parameter, no registered route, with a computed non-expressibility witness →
  `IRREDUCIBLY_INDEPENDENT` → `SECOND_MEDIUM_DRIFT`.
- **Reducible-derived (origin):** a synthetic parameter with an earned closed relation in `{ρ,K,m,a}` → recomputes `REDUCIBLE_DERIVED`.
- **Contradiction (headline):** two sectors pinned to clashing `ρ_br` values (or the `C_hu` charge-residue tie) → `NO_GO` + the conflicting
  relation / unsat core.
- **Freedom-revoked (STATE + rerun):** assert a `ρ_br=varrho_br` closed integral → `FREEDOM_REVOKED_DERIVED` → a pathA_40-rerun SAT check
  (the rerun may stay SAT — then the NG5 headline can remain CONDITIONAL because `c_E`/Route-A is still sim-deferred; the control proves
  the STATE transition + rerun fired, not necessarily a headline flip).
- **Baseline PRODUCTION run:** the real earned ledger → reports the origin table + per-row sim-deferred solves + drift count (NOT an
  expected-landing assertion; a non-`CONDITIONAL` computed result is preserved + triggers scrutiny).

---

## 7. Dual-engine (BINDING — split genuine DERIVATION from EXTRACTION+AUDIT; Codex BLOCKER)
**Genuinely dual-engine DERIVABLE (SymPy + INDEPENDENT Mathematica, second engine derives not echoes):** (a) the §3.1 DIMENSION check —
both engines independently confirm `∫_layer dn·m·ρ` reduces to M L⁻³ (the `∫dn` layer-normal L fixing bulk `m·ρ` M L⁻⁴ → surface M L⁻³),
and check whether `ρ_br` carries any residual free multiplier beyond `varrho_br`; (b) current-ledger non-entailment witnesses for
`c_E²ρ_br−μ_R=0` and any asserted tie; (c) the tie/rerun SAT algebra (reuses the pathA_40 lock+stability system); (d) the `ℓ_g`
substrate-status check (is `ℓ_g` reducible to `a`, or a free continuous parameter? — a SEPARATE check from (a)'s residual-multiplier).
**EXTRACTION + PROVENANCE AUDIT (both engines/agents parse the frozen actions + confirm, but this is NOT symbolic derivation):** the §3.1
ROLE/NAME lineage — that pathA_35 symbol `ρ_br` IS the same lineage object as pathA_25 `varrho_br[ρ]` (same kinetic-term coefficient
position, cited `:line`); the §3.4 branch-status; the registered-route-vs-imagined check. **Bookkeeping (single-engine, cited):** incidence,
drift count, table formatting. Frame the second engine as independent EXTRACTION+AUDIT wherever a stage is classification, NOT derivation.
Orchestrator arbiter re-run + transliteration-fidelity + adversarial-with-ablation on fresh agents.

---

## 8. Scope boundary + the "gravity 0 new" check
NG5 = the origin ledger + the `varrho_br/ρ_br` lineage adjudication + the pathA_40 freedom resolution + rerun triggers. NOT the `pde_ledger`
assembly (gate 3). **"Gravity 0 new" must be justified or dropped:** `c_L1, c_L2` are the active-baseline smectic driver coefficients and
the pathA_25 "B gravity bundle" gate is pending on B4's background — if the gravity bundle consumes Family-L inputs, gravity is not 0 new
(→ those inputs re-enter the §3.4 active-path test).

---

## 9. Deliverables + review plan
- `tools/pathA_41_ng5_second_medium_drift_{sympy.py,.wl}` (dual-engine: §3 core + §6 controls; exit 0; engine-agreement asserted).
- `reports/pathA_41_ng5_second_medium_drift.md` (verdict line 1) + `_results.yaml` (the full origin ledger, the same-object adjudication,
  the freedom states + rerun triggers, the drift count, the §3.1 lineage finding (`OVERCOUNT_OR_SMUGGLE` on SAME / `DIFFERENT_BRANE_INERTIA_OBJECTS`
  on DIFFERENT), controls, provenance).
- **Gauntlet:** this directive → Codex design-review (xhigh) → fold → GLM-5.2 tertiary → Codex confirm-to-green → **user gate** → execute
  (Codex codes, Claude reviews) → arbiter re-run + fidelity + adversarial-with-ablation → **user gate** → commit. If a pathA_40 rerun is
  triggered, run it before the roll-up. Never alter the calibrated process unilaterally.

---

## 10. Changelog
- **v1 → v2 (GLM-5.2 tertiary folded; verdict `SOUND`, no BLOCKERs).** F2 (significant, escalating): added the **binding forward re-open
  trigger** (§3.3) — `ONE_MEDIUM_CONDITIONAL` carries an expiry; when B4/Route-A is solved, NG5 MUST re-open and the sim-deferred rows
  re-enter the origin competition (so the conditional can't become a permanent laundered pass). E1: §3.4 **production-collapse
  transparency** (all 8 Family-L/R/C drivers land SIM_DEFERRED/OUT_OF_ACTIVE in the expected case — disclosed, not hidden). NITs: F3 (§3.1
  "frozen profiles" → "frozen DEFINITIONS", profiles unbuilt); F4 (§6 drift-one-level-down: "no registered route" is pipeline-DERIVED, not
  control-asserted). **F5 escalated to the user** (pathA_35 erratum — see Status). E2 (verify baseline reports per-row sim-deferred solves)
  deferred to the execution/tri-review stage. GLM confirmed the `C_hu`→`FREEDOM_CONTRADICTED`→`NO_GO` route is a genuine production-reachable
  failure path (E5) and no param is missing a row (E3).
- **v0 → v1 (Codex design-review xhigh folded).** BLOCKERs: (1) §3.2 `ρ_br` origin made CONDITIONAL on §3.1's branch — SAME inherits
  `varrho_br` (SIM_DEFERRED/DERIVED), DIFFERENT gives pathA_35 `ρ_br` its own row + possible `IRREDUCIBLY_INDEPENDENT` (the falsifier is no
  longer hardcoded away); (2) §3.3 adds the explicit `B4_UNSOLVED_UNKNOWN → SIM_DEFERRED` production landing (the B4 solve is unbuilt →
  not computable now → do NOT fabricate a drift verdict), and §3.4 splits DECIDABLE branch-status from UNDECIDABLE knit-propagation; (3)
  `OVERCOUNT_OR_SMUGGLE` made SAME-branch-only + `DIFFERENT_BRANE_INERTIA_OBJECTS` added, §5 rider de-hardcoded; (4) §6 controls rewritten
  with source-mutation + recomputed-invariant + expected-transition, distinguishing headline-flips from state/finding-flips, overcount
  control de-rigged (recompute from 5 invariants, not a re-tag); (5) §7 splits genuine dual-engine DERIVATION (dimension, residual-
  multiplier, SAT, ℓ_g-status) from EXTRACTION+AUDIT (the ρ_br=varrho_br role/name lineage — cited, not derived). NITs: §5 `registered-solves`
  listed per-row not hardcoded {B4,Route-A}; §0 typo pathA_35→pathA_25 fixed; §7 ℓ_g-status vs residual-multiplier separated; note pathA_35's
  committed `SECOND_MEDIUM_DRIFT_AT_FREEZE` may be overcounted (recording in NG5 + roll-up sufficient; no pathA_35 amendment prerequisite).
- **v0** — authored from `_scratch/consistency_knit_ng5_scope_v3.md` (`SOUND_AS_A_SCOPE`; Codex×2 + GLM-5.2 folded — the `varrho_br[ρ]`
  lineage BLOCKER B1 + substrate/ℓ_g B2 + B4-anchor B3 + Family-L/R/C B4 + the target-blind/non-retroactive/gravity-0-new/c_s-function NITs).
