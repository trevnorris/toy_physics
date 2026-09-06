# Build directive — S11c-c2 N6 representation-invariance control (the untested one)

## Role of this artifact
Orchestrator-written **build directive** governing a **new SymPy companion diagnostic**. It computes the
**corrected** N6 representation-invariance control that `directives/S11c_c2_SHARED_PHYSICS.md` §5c specifies but that
the c2 audit **never ran** — the audit's `REP_INVARIANCE_*` loop is the obsolete cross-anchoring proxy
(`scripts/S11c_c2_selfenergy_fold_sympy_audit.py:1087-1107`), the wrong object.

- **Deliverable:** `research/pde_ledger_v3/scripts/S11c_c2_N6_diagnostic_sympy.py` (a dedicated companion, sibling of
  `S11c_c2_FG_diagnostic_sympy.py`), plus its transcript OUTSIDE the repo.
- **Governing physics (SUPPLIED — the ONLY spec you read):** `directives/S11c_c2_SHARED_PHYSICS.md` §5c (the object),
  with §3a (the close operation), §3c (the weak increment + slot-linearity + `extract`), §1d (routing), §5d
  (density), §6 (dimensions/multigrade contract). ⛔ Do NOT read any `_measurements/` design or adjudication record —
  they are orchestrator notes that state a physics expectation you must not see.
- **This is a pure-SymPy build. No Mathematica.** Watch memory/time (§ "Budget").

⛔ **You are BUILDING the object §5c NAMES. You are NOT handed, and must NOT invent, any expected value for any
residual.** §5c: *"its computed value is the finding, ⛔ no target value is supplied, and the diff is adjudicated on
our side (⛔ never a builder exit condition)."* Your job ENDS at compute-emit. The residual test is a finite-field
PIT whose numeric output we adjudicate — ⛔ **residual-zero is NEVER a builder exit / assert / pass condition.**

---

## 1 · The object — §5c's two coordinate routes at a FIXED anchoring, MAP-THEN-EXTRACT

For **each fixed anchoring `α`** and density representative `ρ` (⛔ the two routes are the *representation* axis; the
anchoring and density are HELD FIXED across the two routes), construct the c2 self-energy increment two independent
ways and difference them:

```text
route 1 (Eulerian):    I_E^{α,ρ}     = extract( close(SLAB)   − SLAB   )
route 2 (material→E):   I_{M→E}^{α,ρ} = extract( T_{M→E}[ close(SLAB_M) − SLAB_M ] )
R_N6^{α,ρ} = I_E^{α,ρ} − I_{M→E}^{α,ρ}
```

⭐⭐ **MAP-THEN-EXTRACT is load-bearing — ⛔ NEVER extract-then-map.** `extract` (§3c,
`scripts/S11c_c2_selfenergy_fold_sympy_audit.py:325-342`) is a **lossy, intrinsically-Eulerian projection** (Helmholtz
`restricted(·,TRANSVERSE/LONGITUDINAL)`, div/curl, fixed Eulerian test functions `s11cc2Test*`), and `T_{M→E}` mixes
sectors (`θ → θ + u·∇ρ⁰/ρ`), so the two do **NOT commute**: extracting the material rows first and mapping afterward
silently drops the advective off-diagonal — exactly the channel N6 tests. Apply `T_{M→E}` to the material **carrier
increment** first, then the SINGLE Eulerian `extract` on both routes. **There is no `extract_M`.** This mirrors the
parent house pattern (`scripts/S11c_b_brane_operator_sympy_audit.py:2762-2783`: `material_pullback` on the density →
`operator_from_density`).

Emit, keyed by object + anchoring + density (§5c tags):
```text
S11CC2_REP_INVARIANCE_EULERIAN_OPERAND[α,ρ] = I_E^{α,ρ}
S11CC2_REP_INVARIANCE_MATERIAL_OPERAND[α,ρ] = I_{M→E}^{α,ρ}
S11CC2_REP_INVARIANCE_RESIDUAL[α,ρ]         = R_N6^{α,ρ}
```
⛔ `Δρ = δρ_E + u·∇ρ⁰` is the **intra-anchoring** representation map; it ⛔ NEVER bridges `LAB_HELD ↔
MATERIAL_ADVECTED` (category error, §5c/`c1 §5a`). Both routes are at the SAME `α`.

Definitions (all §5c, at the SAME fixed `α,ρ` — read them verbatim, ⛔ do not paraphrase into shortcuts):
- **`SLAB`** = the imported Eulerian `slab_operator` for this `α`.
- **`close(X)`** = the §3a substitution of the **imported same-`α` c1 closed face response** `s11c_c1_face_response`
  into `X`; **`extract`** = the single §3c Eulerian weak extraction. Route 1 is c2's current Eulerian increment.
- **`SLAB_M`** = the **native material-coordinate** δp-slot pressure carrier at this same `α` — built from the
  **S11c-a MATERIAL face sources** (traction / `closure_shape_deriv` / face velocity `V_s`) folded into the **same δp
  slots** the S11c-b face-fold uses (a NAMED existing construction, ⛔ not a new recipe). ⚠ The S11c-a material face
  builder **already maps its covector to Eulerian** (inverse-transpose) — account for that, ⛔ no double-map. ⛔⛔ There
  is **no ready-made callable** for `SLAB_M`: `build_operator(route="MATERIAL")` folds the bulk map in and keeps
  Eulerian faces — ⛔ it is NOT `SLAB_M`.
- **`T_{M→E}`** = the `N4` representation map: `θ → θ + u·∇ρ⁰/ρ`, the anchoring-branched `e_W` shift (`e_W → e_W +
  u·∇W/W` for `LAB_HELD`, unchanged for `MATERIAL`), × Jacobian `1+tr(∇u)`. ⛔⛔ **It is NOT `material_pullback`
  applied to the carrier/extracted rows** (see §3): that callable is a bulk-density quadratic map whose final 2nd
  scale-derivative ANNIHILATES the linear carrier rows, never touches the δp face slots, and (trial-only) misses the
  reverse u-row channel `δθ_M = δθ_E + δu·∇ρ/ρ` — the off-diagonal N6 tests. The map is realized at the **action /
  virtual-work SCALAR level with the Eulerian variation taken AFTER**, then the single Eulerian `extract`. ⛔ Do NOT
  re-transform the already-Eulerian imported `close` response.

---

## 2 · TRACTABILITY MANDATE — carrier-first + finite-field PIT, ⛔ never full-symbolic, ⛔ never `build_case` end-to-end

⚠ The naive build (materialize whole symbolic slab objects, difference them, symbolically zero-test the residual over
every retained grade) is the **exact wall the earlier F/G diagnostic hit** — ~11 min / ~150 MB per case,
un-adjudicable. Design tractable **from the start**:

**(a) CARRIER-FIRST, BOTH ROUTES.** By §3c linearity the increment is the pressure-slot carrier only
(`P = {δp_s, ∂_w δp_s}`; the pressure-independent bulk/kinetic base cancels in `close(·) − (·)`). Build only
`S_{P,s} = S − S|_{P=0}` per face, per route — reuse the F/G companion's `carrier()` / per-face split
(`scripts/S11c_c2_FG_diagnostic_sympy.py:99-134`). ⛔⛔ **Do NOT call `build_case` end-to-end for EITHER route** — it
closes+truncates+extracts full rows before returning, recreating the wall (the F/G companion's `build_case`-then-
carrier was part of its wall). Reuse only the **low-level** binding / face-response / `extract` helpers, and
**evaluate/truncate coefficients EARLY, before closure expansion** (a random/numeric substitution applied AFTER the
full rows expand does NOT solve the wall). Include every route-specific closure ingredient (`μ_θ`, face velocity,
normals, slot coefficients).

**(b) Slot-linearity is GUARDED, not assumed — via the SAME finite-field PIT.** Emit the guard (that the increment is
the pressure-slot carrier, base cancels) as a **computed** structural check evaluated by the PIT below. ⛔ Do NOT
materialize a second full-symbolic increment to "guard" the first (that IS the wall).

**(c) Residual test = EXACT FINITE-FIELD PIT, ⛔ not full-symbolic, ⛔ not bare float.** Test `R_N6` (and each control
residual, and the slot-linearity guard) by polynomial identity testing:
- **Split** by retained `(ε,η,σ_W)` grade / the six weak blocks / formal-integral signature / non-`Integral`
  remainder (normalize only the permitted compact-support in-plane IBP).
- **Evaluate** each block's arithmetic circuit over **several exact finite fields (primes)** at random generic sample
  points **shared across the two routes**; clear denominators, reject singular samples, cover every `Piecewise`
  branch cell away from transition/singular surfaces.
- **Record** seeds, primes, degree bounds, and the resulting **false-negative bound** δ.
- A **certified nonzero modular numerator is DECISIVE** — a real N6 finding. **All-zero** gives only *"no nonzero
  found; conditional FN prob ≤ δ"* — a bounded conditional pass, ⛔ NOT a deductive proof, ⛔ NOT a build success gate.
  Emit **compact numeric** operand/residual/provenance tables, ⛔ NOT the giant symbolic expressions.

**(d) What to randomize / hold** (so it cannot be an accidental zero):
- Free fields / jets / params = random **generic**, **shared across both routes**; reject singular denominators /
  special limits.
- Keep **BOTH** `δp_s` families (`plus`, `minus`) generic.
- **Alpha-align the integral bound dummies**, then sample kernel coords consistently — ⛔ **never `k_in = k_out`**, ⛔
  never evaluate the Fourier/DtN integral.
- `Z` / resolvents / Fourier = the **SAME** definitions on both routes, sampled generically-but-coherently — ⛔ never
  `Z→0`, ⛔ never scalar, ⛔ never independent-random per route.

---

## 3 · The two independent constructions — build discipline

**Route 1 (Eulerian) — carrier-first, reuse low-level audit helpers.** Build the Eulerian pressure carrier
`close(SLAB) − SLAB` from the imported Eulerian `slab_operator` closed with the imported same-`α` c1 response, then
the single §3c `extract`. Reuse `bind_inputs / load_model / extract / retained_shape / difference` (as the F/G
companion does), ⛔ but not `build_case` end-to-end (§2a).

**Route 2 (material→E) — you DESIGN the typed pipeline; this names the OBJECT + REFS + definition-of-done.**
⭐ The exact route-2 CAS pipeline is yours to design and build — ⛔ the orchestrator does NOT hand you a recipe (two
prior recipes were type-wrong). Build the object §5c names, satisfying every DoD clause below; the **build review legs
will gate the construction** (types, annihilation, reverse channel, double-map).

*Verified refs (reuse these; ⛔ do not re-derive):*
- native material δp-slot carrier: the **S11c-a MATERIAL face sources** (traction / `closure_shape_deriv` / face
  velocity `V_s`, `scripts/S11c_a_interface_geometry_sympy_audit.py:703-845`) folded into the **same δp slots** the
  **S11c-b face-fold** uses; ⚠ that S11c-a builder **already maps its covector to Eulerian** (inverse-transpose,
  `:742-755`) — account for it, ⛔ no double-map.
- the `N4` field redefinition / Jacobian: `material_pullback` (`scripts/S11c_b_brane_operator_sympy_audit.py:1916-1981`)
  — read it as the **definition of the field map**, ⛔ do NOT call it on carrier/extracted rows (its 2nd scale-
  derivative annihilates linear rows; it never touches δp slots).
- `close` = the same imported same-`α` c1 response (`s11c_c1_face_response`), with route-2 `V_s` = **material** face
  velocity; ⛔ do not reconstruct the DtN; ⛔ do not re-transform the already-Eulerian response.
- the trial-only `representation_pullback` (`scripts/S11c_c2_selfenergy_fold_sympy_audit.py:1136-1171`) — ⛔ NOT the
  map (it misses the reverse channel).

*Definition of done (the build legs will check each):*
1. **Native, not a re-pullback of route 1** — route 2 is built from material face sources + material `V_s`, ⛔ never by
   transforming route 1's operand.
2. **Map-then-extract, at the SCALAR level** — any bulk `N4` field redefinition is applied to the action / virtual-
   work scalar with the **Eulerian variation taken AFTER**, then the SINGLE Eulerian `extract`. ⛔ Never `extract` on
   unmapped material rows (Eulerian tests ⇒ ill-posed).
3. **No annihilation** — the material increment must be a nonzero linear carrier; ⛔ if a step (e.g. a 2nd scale-
   derivative on linear rows) sends `I_{M→E}→0`, that is a construction bug — emit the intermediate size/grade as a
   computed check.
4. **Reverse channel present** — the reverse u-row off-diagonal (`δθ_M = δθ_E + δu·∇ρ/ρ`) must appear, not just the
   forward θ-row; a trial-only map that drops it is wrong.
5. **No double-map** — do not map geometry the S11c-a material builder already mapped.
6. **Carrier-first, no `build_case` end-to-end** (§2a); truncate/evaluate coefficients before closure expansion.

**Independence is the whole point (build-skill corollary 3 + one-sided corruption).** Route 1 (native Eulerian) and
route 2 (native material) must be structurally different constructions; verify by §4's one-sided corruptions — if
corrupting one route moves the other's operand, they were never independent and `R_N6` is decoration.

---

## 4 · The independence controls (§5c "then …") — one representation route only

At the SAME fixed `α,ρ`, mutate **one representation route only at its source** and recompute it, leaving the other
unchanged. Emit `S11CC2_CONTROL_INDEPENDENCE_{BASE,CORRUPTED,RESIDUAL}[α,ρ,probe]` — **PRINT** operands and residual
(via the same PIT), ⛔ do not assert; the disposition (corrupted route moved, uncorrupted route did not) is
adjudicated on our side.

1. **Tilt probe** — reverse one face's first-jet slope term in `n̂_s` on the **Eulerian** route. ⛔ The imported
   operator is already materialized and the c2 slope override only alters the c1 DtN kernel jet — that does NOT move
   the slab carrier's normal-derived coefficients. Build BOTH the uncorrupted baseline and the slope-corrupted operand
   through the **same injectable Eulerian carrier factory** (at the S11c-a face-source level), and **PIT-check the
   uncorrupted factory output against the imported carrier** before using the control (so a nonzero residual is the
   mutation, ⛔ not reconstruction drift). The material operand keeps its provenance/sample values.
2. **N4 advection probe** — omit the map's advective term `u·∇ρ⁰/ρ` (the `N4` `θ`-shift) only inside `T_{M→E}`, AFTER
   constructing the native material carrier increment (⛔ NOT from the `δp_s`-independent slab base, which cancels per
   §3c; the Eulerian route is unchanged).

⛔ **NO `A−A` control.** That advective term is **structurally ABSENT for `RHO4_CONSTANT`** (`ρ₄=ρ_br/W₀` ⇒ `∇ρ₄=0`)
and PRESENT for `RHOBR_CONSTANT` (`ρ₄=ρ_br/W_bg`, live `W_bg`) — the **live probe is `RHOBR_CONSTANT`**; for
`RHO4_CONSTANT` **emit the computed absence as a structural fact**, ⛔ do not difference an operand against an
identical copy. ⛔ Corrupting one *anchoring* is NOT this test.

---

## 5 · The audit correction — REPLACE the obsolete cross-anchoring loop (an EDIT, not companion-only)

The c2 audit's `for density in DENSITIES:` block (`scripts/S11c_c2_selfenergy_fold_sympy_audit.py:1087-1107`) emits
`REP_INVARIANCE_{LAB,MATERIAL}_OPERAND` + `REP_INVARIANCE_RESIDUAL` as a **cross-anchoring** difference via the
trial-field-only `representation_pullback` — the wrong object. ⛔ A companion-only re-emission would leave those
mislabeled audit outputs alive. **Edit the audit:** delete/replace lines `1087-1107` and emit the raw, unmapped,
correctly-oriented object (§5c:366, mandatory, ⛔ not optional):
```text
S11CC2_ANCHORING_L_MINUS_M[ρ] = increment[LAB_HELD, ρ] − increment[MATERIAL_ADVECTED, ρ]
```
— raw per-anchoring increments, **no `representation_pullback`**, orientation **LAB − MATERIAL** (⛔ not the old
`MATERIAL − LAB` sign), on the same footing as `S11CC2_DENSITY_LIVE_MINUS_FROZEN` (§5d): a computed outcome, ⛔ **no
prescribed zero target**, ⛔ never labeled a representation-invariance residual. ⛔ Do NOT touch the density-difference
loop (`:1085-1086`, `DENSITY_LIVE_MINUS_FROZEN`) — §5d, stays. The old `REP_INDEPENDENCE_*` corruptions are
superseded by §4; do not carry the old `REP_INVARIANCE_*` names forward.

---

## 6 · The three non-negotiable script clauses — carry them verbatim

> **1. The script may PRINT computed objects. It may NOT state conclusions.** An emit payload must be a CAS object (an
> expression, a solved root, a boolean from a symbolic/modular test), ⛔ never prose describing a result.
> **2. PRINT the residual; do NOT assert it.** `assert residual == 0` is the builder writing down the expected
> output. Compute → emit → (never) assert. Here the residual is a PIT numeric table; residual-zero is ⛔ never an exit.
> **3. Interpretation belongs to the STEP RECORD.** The script does not editorialise.

**Four corollaries (all binding):**
1. ⛔ **A hand-typed CAS object is still hand-typed.** The ONLY place the physical symbols may be combined by hand is
   in constructing the ACTION and the ANSATZ. **Every other expression involving them must be REACHED BY
   COMPUTATION. Every control re-enters the chain at the ACTION/source, ⛔ never at a result.**
2. ⛔ **The tag name is output too.** A name may name **the object**, ⛔ never its value, sign, ratio, or the shape of
   the answer.
3. ⛔ **No tautological residual.** Emit a difference only where the two operands were produced by INDEPENDENT ROUTES
   (route 1 native-Eulerian, route 2 native-material-then-mapped — verified independent by §4's one-sided
   corruption). ⛔ If `q:=A/B` and you emit `A − q·B`, it is zero by construction for any input. Where no second route
   exists, emit the objects and say so.
4. ⛔ **Emission must never be conditional on a payload's value.** Whether a tag appears may depend only on which
   object/representative it belongs to — ⛔ never on whether its value is zero. Tag names unique; payloads may
   legitimately repeat, and that repetition is a finding.

---

## 7 · Output discipline + provenance + dimensions

- `emit(name, payload, **labels)` prints one JSON line per object (mirror the F/G companion's `emit`). Every object
  carries `anchoring`, `density`, and (for controls) `probe` labels.
- Payloads are **compact numeric** — modular evaluations, sample tables, structural booleans, provenance — ⛔ never
  giant symbolic slab dumps.
- **Provenance block:** seeds, primes, per-block degree bounds, the false-negative bound δ, rejected singular
  samples, branch cells covered, and the slot-linearity guard result.
- **Dimensions (§6 contract):** restore `[L,T,M]` on every emitted object with an **able-to-fail**
  dimensional-consistency result, and carry `(ε,η,σ_W)` multigrade (`N12`) on every object. A numeric-only build can
  hide a dimensionally invalid map-back — this check must be able to fail.
- Emit **operands before the residual**, then the residual — for the invariance test AND each control.
- Run `reduction/derived_or_declared.py` on the deliverable after the build; a script still mostly CONSTANT has not
  been built (triage, not a verdict).

## 8 · Budget / ops
- Single Python process; ⛔ no shell `timeout` wrapper (SIGKILL has cost 300k+ tokens). Print progress so a long run
  is observably alive (a long+silent run is the failure mode).
- The carrier-first + early-truncation design exists to keep memory bounded. If a single block still balloons (whole-
  object symbolic materialization), STOP and emit the block's size as a computed finding — ⛔ do NOT push to
  full-symbolic. ⛔ Never run a second memory-heavy CAS job alongside this.

## 9 · Supplied vs withheld (state it in the transcript)
- **SUPPLIED and therefore unfalsifiable within this build:** the §5c object definition, §3a close, §3c weak
  increment + slot-linearity + `extract`, the c1 imported face response, the S11c-a/b material builders. A passing
  build does NOT verify these — it computes the N6 residual GIVEN them.
- **WITHHELD:** any expected value / acceptance criterion for `R_N6` or the controls. There is none to supply; §5c
  fixes no target. The diff is adjudicated on our side. ⛔ Do not read the `_measurements/` records; ⛔ do not iterate
  toward a zero residual; a certified nonzero is the most valuable output available.
