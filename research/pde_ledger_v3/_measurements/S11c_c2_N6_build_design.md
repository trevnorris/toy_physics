# N6 build correction — vetted TRACTABLE design (2026-09-06)

⚠ **Orchestrator design record — ⛔ NOT builder-facing** (it states the physics expectation; the builder reads only
the directive + §5c). The corrected N6 test is defined by **§5c** (axis correction `30d4b72d`; **route-2 ordering
correction** folded 2026-09-06 from the directive-review r1 adjudication,
`S11c_c2_N6_directive_review_r1_adjudication.md`): within each FIXED anchoring `α` and density `ρ`, compare the
self-energy increment built two ways — Eulerian `I_E = extract(close(SLAB) − SLAB)` vs material `I_{M→E} =
extract( T_{M→E}[ close(SLAB_M) − SLAB_M ] )` — residual `R_N6 = I_E − I_{M→E}` (a computed finding); plus a
one-sided independence corruption. ⭐⭐ **MAP-THEN-EXTRACT, ⛔ never extract-then-map:** `extract` is a lossy Eulerian
projection and `T_{M→E}` mixes sectors, so they do not commute (the parent does map-then-extract). The c2 audit never
ran this control (its "representation" loop compares the two anchorings — the obsolete cross-anchoring proxy).

**Codex-sol question-vet** (`scratchpad/codex_N6_build_vet.txt`, ~500 KB) established the tractable shape. The naive
full-symbolic build hits the exact F/G wall; the design below avoids it. Directive-review r1 (Codex-sol + Grok) then
found the extract-then-map ordering defect + the callable/leak/dimension fixes now folded here.

## The design — CARRIER-FIRST + per-grade exact finite-field PIT

1. **Carrier-first, ⛔ NOT a full symbolic slab, ⛔ NOT `build_case` end-to-end (BOTH routes).** By §3c linearity the
   increment is the pressure-slot carrier only (`P={δp_s,∂_w δp_s}`; the pressure-independent bulk/kinetic base
   cancels). Map-then-extract:
   ```
   I_E     = extract( Σ_s ( S^E_{P,s}[χ_s] − S^E_{P,s} ) )
   I_{M→E} = extract( T_{M→E}[ Σ_s ( S^M_{P,s}[χ_s] − S^M_{P,s} ) ] )
   ```
   Build only the carrier `S_{P,s} = S − S|_{P=0}` per face, per route (the F/G companion
   `S11c_c2_FG_diagnostic_sympy.py:99-134` implements this split). ⛔ Do NOT call `build_case` end-to-end for EITHER
   route (Codex-3: it closes+truncates+extracts full rows before returning — recreates the wall; the F/G companion's
   `build_case`-then-carrier was part of its wall). Reuse only the low-level binding/face-response/extract helpers;
   **evaluate/truncate coefficients EARLY, before closure expansion**. Compute the carrier's slot-linearity
   **guarded, not assumed — via the SAME finite-field PIT** (Grok-3: ⛔ never a second full-symbolic increment to
   "guard" the first). Include every route-specific closure ingredient (`μ_θ`, face velocity, normals, slot
   coefficients).
2. **Native material δp-slot carrier `SLAB_M` from MATERIAL FACE SOURCES, ⛔ not `build_operator(MATERIAL)`
   wholesale.** ⛔⛔ `build_operator(route="MATERIAL")` (`:2762-2783`) folds the bulk map in + keeps Eulerian faces — NOT
   `SLAB_M`. Build it from the **S11c-a MATERIAL face sources** (traction / `closure_shape_deriv` / face velocity `V_s`,
   `S11c_a_interface_geometry_sympy_audit.py:703-845`) folded into the **same δp slots** the S11c-b face-fold uses. ⚠
   That S11c-a builder **already maps the covector to Eulerian** (inverse-transpose `:742-755`) — account for it, ⛔ no
   double-map (Codex-r2-1 / Grok-r2-F1).
3. **Map-then-extract; ⛔ `T_{M→E}` is NOT `material_pullback` on the carrier rows.** ⛔⛔ `material_pullback`
   (`:1916-1981`) is a bulk-density quadratic map — its final 2nd scale-derivative **ANNIHILATES the linear carrier
   rows** (`I_{M→E}≈0`), it never touches δp face slots, and (trial-only) misses the **reverse u-row channel** (the
   N6 off-diagonal). Read it as the DEFINITION of the field map only. The `N4` representation difference is realized
   **at the action / virtual-work SCALAR level (native material faces + any bulk field redefinition on the scalar),
   Eulerian variation taken AFTER**, then the SINGLE Eulerian `extract` (⛔ never `extract` on unmapped material rows —
   ill-posed; ⛔ never the trial-only `representation_pullback` `selfenergy_fold:1136-1171`; ⛔ never re-transform the
   already-Eulerian `close` response). ⭐ **The exact route-2 pipeline is astra's to design (build legs gate);** open
   question — whether a bulk field map is needed at all, or native material faces + material `V_s` already give the
   material increment.
4. **Close both routes with the SAME imported same-`α` c1 response** (⛔ do not reconstruct the DtN).
5. **Residual test = exact finite-field PIT** (⛔ not full-symbolic, ⛔ not bare float): split by retained
   `(ε,η,σ_W)` grade / six weak blocks / formal-integral signature / non-`Integral` remainder (normalize only the
   permitted compact-support in-plane IBP); evaluate the residual arithmetic circuits over **several exact finite
   fields (primes)**, clearing denominators, rejecting singular samples, covering Piecewise branch cells, recording
   seeds/primes/degree-bounds + the resulting **false-negative bound**. A certified nonzero modular numerator is
   **decisive** (a real N6 finding); all-zero gives only "no nonzero found; conditional FN prob ≤ δ" (⛔ not a
   deductive proof). Emit **compact numeric** operand/residual tables + provenance, ⛔ NOT giant symbolic; residual-zero
   is ⛔ NEVER a builder exit condition (§5c "computed finding" discipline).
   - **What to randomize/hold:** free fields/jets/params = random generic, **shared across routes**, reject singular
     denominators + special limits; keep BOTH `δp_s` families generic; alpha-align integral dummies then sample kernel
     coords consistently (⛔ never `k_in=k_out`, ⛔ never evaluate the integral); `Z`/resolvents/Fourier = SAME
     definitions both routes, sampled generically-but-coherently (⛔ never `Z→0`/scalar/independent-random); test every
     admissible branch cell away from transition/singular surfaces.
6. **Controls:** **Tilt** — reverse one face's `n̂_s` slope on the **Eulerian** route, built through the **same
   injectable carrier factory** as that route's baseline, PIT-checked against the imported carrier first (⛔ NOT by
   editing the imported operator, ⛔ NOT via the c1-DtN-kernel-jet override, which does not move the slab carrier's
   normal-derived coefficients — Codex-r2-3). **Advection** — omit the map's advective term `u·∇ρ⁰/ρ` (the `N4`
   `θ`-shift) only inside `T_{M→E}` (Eulerian route unchanged). ⚠ **Structurally ABSENT for `RHO4_CONSTANT`**
   (`ρ₄=ρ_br/W₀` ⇒ `∇ρ₄=0`), PRESENT for `RHOBR_CONSTANT` (live `W_bg`) — the live probe is `RHOBR_CONSTANT`; for
   `RHO4_CONSTANT` emit the computed absence, ⛔ NOT an `A−A` (Codex-r2-2 / Grok-r2-F2; my directive had this reversed).
7. **Dimensions (Codex-5b, §6 contract):** restore `[L,T,M]` + an **able-to-fail** dimensional-consistency result +
   `(ε,η,σ_W)` multigrade (`N12`) on every emitted object.

## Where it lives
A **dedicated companion** script (like the F/G one), reusing binding/provenance, pressure-response/resolvent,
retained-grade algebra, weak `extract`, integral-channel alignment, and the carrier split. Plus a **required EDIT to
the c2 audit** (Codex-4, ⛔ NOT companion-only — else the old mislabeled outputs survive): **replace** its
cross-anchoring "representation" loop (`S11c_c2_selfenergy_fold_sympy_audit.py:1087-1107`) with the RAW unmapped
`ANCHORING_L_MINUS_M[ρ] = increment[LAB_HELD] − increment[MATERIAL_ADVECTED]` (§5c:366 orientation LAB−MATERIAL, ⛔ no
`representation_pullback`, ⛔ not the old MATERIAL−LAB sign); **keep** the density loop (`:1085-1086`,
`DENSITY_LIVE_MINUS_FROZEN`). ⛔ The builder reads only the directive + §5c — ⛔ NOT this design doc or the
adjudication records (they state the "must vanish" expectation; Codex-5a leak).

## Process (corrected)
Build directive (this design) → **2 decision legs** (Codex-sol + Grok, review-until-clear since physics-bearing) →
**astra build** → **2 build legs** (fresh Claude + Grok) → adjudicate. ⛔ Design tractable FROM THE START (carrier +
finite-field PIT) — do NOT let it regress to full-symbolic (the F/G wall).
