# Build directive — S11c-c2 N6 representation-invariance control (the untested one)

## Role of this artifact
Orchestrator-written **build directive** governing a **new SymPy companion diagnostic**. It computes the
**corrected** N6 representation-invariance control that `directives/S11c_c2_SHARED_PHYSICS.md` §5c specifies but that
the c2 audit **never ran** — the audit's `REP_INVARIANCE_*` loop is the obsolete cross-anchoring proxy
(`scripts/S11c_c2_selfenergy_fold_sympy_audit.py:1087-1107`), the wrong object.

- **Deliverable:** `research/pde_ledger_v3/scripts/S11c_c2_N6_diagnostic_sympy.py` (a dedicated companion, sibling of
  `S11c_c2_FG_diagnostic_sympy.py`), plus its transcript OUTSIDE the repo.
- **Governing physics (SUPPLIED):** `directives/S11c_c2_SHARED_PHYSICS.md` §5c (the object), with §3a (close), §3c
  (weak increment + slot-linearity + `extract`), §1d (routing), §5d (density), §6 (dimensions/multigrade contract);
  and the **authoritative, reviewed-CLEAR route-2 construction spec**
  `research/pde_ledger_v3/_measurements/S11c_c2_N6_route2_spec_astra.md` (astra-authored, fresh-Claude + Grok CLEAR;
  target withheld) — **build route 2 exactly as it specifies.** ⛔ Do NOT read any OTHER `_measurements/` design or
  adjudication record — those state a physics expectation you must not see (the astra route-2 spec is the sole
  exception, and it withholds the target).
- **This is a pure-SymPy build. No Mathematica.** Watch memory/time (§ "Budget").

⛔ **You are BUILDING the object §5c NAMES. You are NOT handed, and must NOT invent, any expected value for any
residual.** §5c: *"its computed value is the finding, ⛔ no target value is supplied, and the diff is adjudicated on
our side (⛔ never a builder exit condition)."* Your job ENDS at compute-emit. The residual test is a finite-field
PIT whose numeric output we adjudicate — ⛔ **residual-zero is NEVER a builder exit / assert / pass condition.**

---

## 1 · The object — §5c's two coordinate routes at a FIXED anchoring (⛔ NO `T` on the increment)

For **each fixed anchoring `α`** and density representative `ρ` (⛔ the two routes are the *representation* axis; the
anchoring and density are HELD FIXED across the two routes), construct the c2 self-energy increment two independent
ways and difference them:

```text
route 1 (Eulerian):  I_E^{α,ρ}     = extract( close(SLAB)   − SLAB   )
route 2 (material):  I_{M→E}^{α,ρ} = extract( close(SLAB_M) − SLAB_M )
R_N6^{α,ρ} = I_E^{α,ρ} − I_{M→E}^{α,ρ}
```

⭐⭐ **There is NO separate `T` on the increment** — route 2 mirrors the parent `a.task_rep_invariance` (native material
vs Eulerian face quantities, differenced directly): build route 2's carrier from **native material-route face
ingredients** (the material→Eulerian map lives INSIDE those builders, via the covector inverse-transpose) closed with
**material μ_θ**, and difference directly against route 1. ⛔ Do NOT apply a field-redefinition to the carrier / closed
increment / extracted rows: the increment is δp-slot-only (`S_P` carries no `θ`), so a θ-shift is a no-op/annihilating,
and `material_pullback` is a bulk-density quadratic that ANNIHILATES linear rows. ⭐ **Build route 2 EXACTLY per the
cleared astra spec `_measurements/S11c_c2_N6_route2_spec_astra.md`** (the definitions below are the summary; that spec
is authoritative on the pipeline, the μ_θ binding, the controls, and the finite-field PIT contract).

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
  **S11c-a MATERIAL face sources** (traction / `closure_shape_deriv` / virtual work / face velocity `V_s` via the
  `route="MATERIAL"` builders, already Eulerian-mapped via the covector inverse-transpose) folded into the **same δp
  symbols** the S11c-b face-fold uses, closed with **material μ_θ**. ⚠ The material face builder ALREADY maps the
  covector — ⛔ no double-map. ⛔⛔ **No ready-made callable** for `SLAB_M`: `build_operator(route="MATERIAL")` folds the
  bulk map in but keeps Eulerian faces — a thin material face-fold adapter must be built (astra spec §3, §8).
- **μ_θ binding (the load-bearing decision):** route 2 binds **material μ_θ** (from `operator_from_density` on the
  material-route bulk density) at BOTH the face substrate AND the c1 response source; route 1 binds the imported
  **Eulerian** μ_θ. ⛔ NOT geometry-only (Eulerian μ_θ both routes) — that makes the N4 advection channel vacuous. N4
  enters through the **closed pressure source μ** (the open face rows' direct μ terms are δp-independent and cancel).
- **Reverse-u channel is GRADE-SUPPRESSED, ⛔ NOT a mandatory survival requirement:** the N4-induced reverse-U pressure
  carrier is 0 (`LAB_HELD`, U carrier ≡ 0) or `O(σ_W²)` (`MATERIAL_ADVECTED`, U carrier `O(σ_W)`) — outside the
  retained `(η^{≤1},σ_W^{≤1})` rectangle (verified on the imported operator). ⛔ Do NOT require a reverse-U N4
  contribution to survive; the live N4 signal is the **forward μ/pressure-source** couplings. Report reverse-block
  `t`-sensitivities + **permit computed absence**; ⛔ never credit `extract`'s existing U-row curl as an N4 witness.

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

**Route 2 (material) — build EXACTLY per the cleared astra spec `_measurements/S11c_c2_N6_route2_spec_astra.md`.**
That spec (astra-authored, fresh-Claude + Grok CLEAR) is authoritative on the pipeline; follow its §2–§8. This is the
summary + the definition-of-done the **build review legs will gate**:

*Verified refs (from the astra spec; ⛔ do not re-derive):*
- native material δp-slot carrier: **S11c-a MATERIAL face sources** (traction / `closure_shape_deriv` / virtual work /
  `V_s` via `route="MATERIAL"`, `scripts/S11c_a_interface_geometry_sympy_audit.py:703-845`) folded into the **same δp
  symbols** the **S11c-b face-fold** uses (`face_generalized_force_rows`, `bind_mu_theta_operand` — bind μ_θ BEFORE
  `diff(..., δv)`, `scripts/S11c_b_brane_operator_sympy_audit.py:2135-2179`); ⚠ the S11c-a builder already maps its
  covector to Eulerian (`:742-755`) — ⛔ no double-map.
- **material μ_θ** = `operator_from_density` amplitude on the material-pulled-back bulk *scalar*
  (`material_pullback` `:1916-1981` used ONLY on that scalar, ⛔ never on rows/slots); bound at BOTH the face and the
  c1-source. Route 1 uses the imported Eulerian μ_θ.
- `close` = the same imported same-`(α,ρ)` c1 response (`s11c_c1_face_response`), **opaque** (⛔ no θ-sub / Jacobian
  into `DELTA_P`/`Z`/resolvent); route-2 source slots use material `V_s` + material μ_θ.

*Definition of done (the build legs will check each):*
1. **Native, not a re-pullback of route 1** — material face sources + material `V_s` + material μ_θ, ⛔ never by
   transforming route 1's operand.
2. **No `T` on the increment; single Eulerian `extract`** — the material→Eulerian map lives inside the material face
   builders, ⛔ not a field redefinition applied to the carrier / closed increment / extracted rows.
3. **No annihilation** — `material_pullback` only on the bulk *scalar* to get μ_θ (⛔ never on linear rows/slots).
4. **Reverse-u is GRADE-SUPPRESSED, ⛔ NOT mandatory** — the N4 reverse-U pressure carrier is 0 (`LAB_HELD`) or
   `O(σ_W²)` (`MATERIAL_ADVECTED`), outside the retained rectangle. Report reverse-block `t`-sensitivities and **permit
   computed absence**; the live N4 signal is the forward μ/pressure-source path. ⛔ Never force a reverse channel; ⛔
   never credit `extract`'s U-row curl as an N4 witness.
5. **No double-map** — do not re-map geometry the S11c-a material builder already mapped.
6. **Carrier-first, no `build_case` end-to-end** (§2a); a thin material face-fold adapter is required (astra §8);
   truncate/evaluate coefficients before closure expansion.

**Independence is the whole point.** Route 1 (imported Eulerian) and route 2 (native material) must be structurally
different constructions; verify by §4's one-sided corruptions — if corrupting one route moves the other's operand,
they were never independent and `R_N6` is decoration.

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
2. **N4 advection probe** — a source tag `t` on the material bulk-density θ-advection `u·∇ρ₄/ρ₄` inside the **material
   μ_θ** construction (baseline `t=1`, corrupted `t=0`), the recomputed μ bound at BOTH the face and the c1-source
   (astra §6). ⛔ NOT omitting `BACKGROUND_ADVECTION` from the `δp`-independent mass base (cancels per §3c); the
   Eulerian route is unchanged.

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
   (route 1 native-Eulerian, route 2 native-material face-fold already Eulerian-mapped inside its source builders —
   verified independent by §4's one-sided corruption). ⛔ If `q:=A/B` and you emit `A − q·B`, it is zero by construction for any input. Where no second route
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
  fixes no target. The diff is adjudicated on our side. ⛔ Do not read the `_measurements/` records **other than the
  cleared route-2 spec `S11c_c2_N6_route2_spec_astra.md`** (which withholds the target); ⛔ do not iterate toward a
  zero residual; a certified nonzero is the most valuable output available.
