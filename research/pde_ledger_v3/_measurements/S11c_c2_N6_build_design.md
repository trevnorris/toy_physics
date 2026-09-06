# N6 build correction — vetted TRACTABLE design (2026-09-06)

The corrected N6 test is defined by **§5c** (committed `30d4b72d`): within each FIXED anchoring `α` and density `ρ`,
compare the self-energy increment built two ways — Eulerian `I_E` vs material-coordinate `I_{M→E}` (mapped back by
`Δρ=δρ_E+u·∇ρ⁰`), residual `R_N6 = I_E − I_{M→E}` must vanish; plus a one-sided independence corruption. The c2 audit
never ran this control (its "representation" loop compares the two anchorings — the obsolete cross-anchoring proxy).

**Codex-sol question-vet** (`scratchpad/codex_N6_build_vet.txt`, ~500 KB) established a sound, tractable design. The
naive full-symbolic build hits the exact F/G wall; the design below avoids it.

## The design — CARRIER-FIRST + per-grade exact finite-field PIT

1. **Carrier-first, ⛔ NOT a full symbolic `SLAB_M`.** By §3c linearity the increment is the pressure-slot carrier
   only (`P={δp_s,∂_w δp_s}`; the pressure-independent bulk/kinetic base cancels):
   ```
   I_E     = Σ_s extract( S^E_{P,s}[χ_s] − S^E_{P,s} )
   I_{M→E} = T_{M→E}[ Σ_s extract_M( S^M_{P,s}[χ_s] − S^M_{P,s} ) ]
   ```
   Build only `S_{P,s} = S − S|_{P=0}` per face, per route (the F/G companion `S11c_c2_FG_diagnostic_sympy.py:99-134`
   already implements this split). Compute the carrier's slot-linearity **guarded, not assumed**; include every
   route-specific closure ingredient (`μ_θ`, face velocity, normals, slot coefficients). **Sparse truncated-grade
   representation; evaluate coefficients EARLY** — random substitution AFTER `build_case` expands the rows does NOT
   solve the wall.
2. **Native material route, ⛔ not a re-pullback of the Eulerian.** Construct `S^M_{P,s}` from S11c-b's material
   builders (`S11c_b_brane_operator_sympy_audit.py`: material energy pullback + Jacobian ~1916-1981;
   `build_operator(route="MATERIAL")` ~2723-2777; S11c-a exact material-map/face-flattening
   `S11c_a_interface_geometry_sympy_audit.py:703-815`). ⛔ Do NOT take S11c-b's already-mapped material operand and
   apply the c2 pullback again — native `SLAB_M`, then exactly ONE map-back on the increment.
   ⚠ `build_case` CANNOT build the material route (it selects the imported Eulerian slab; no representation axis).
3. **Full weak-operator map-back `T_{M→E}`** transforms **trial fields AND test covectors + the Jacobian** — ⛔ the
   existing `representation_pullback` is **trial-field-only** (misses the test covector; recorded
   `_measurements/S11c_c2_N6_spec_adjudication_record.md:38-44`), so calling it on another `build_case` result is NOT
   the corrected N6.
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
6. **Controls:** **Tilt** — mutate one face's `n̂_s` slope at the **Eulerian** source, rebuild only that route (the
   material operand keeps its provenance/sample values). **Advection** — omit `u·∇Σ_E⁰` only in `T_{M→E}` AFTER
   constructing the native material increment (Eulerian route unchanged). ⚠ `RHOBR_CONSTANT` (frozen surface density)
   may have the advection source **structurally absent** — emit that absence, ⛔ NOT an `A−A` (§5c forbids A−A controls).

## Where it lives
A **dedicated companion** script (like the F/G one), reusing binding/provenance, pressure-response/resolvent,
retained-grade algebra, weak `extract`, integral-channel alignment, and the carrier split. Plus a **small semantic
correction to the c2 audit**: retire its cross-anchoring "representation" loop (`S11c_c2_selfenergy_fold_sympy_audit.py:1085-1107`)
or re-emit it as the separately-mandated `ANCHORING_L_MINUS_M` (§5c:362-370).

## Process (corrected)
Build directive (this design) → **2 decision legs** (Codex-sol + Grok, review-until-clear since physics-bearing) →
**astra build** → **2 build legs** (fresh Claude + Grok) → adjudicate. ⛔ Design tractable FROM THE START (carrier +
finite-field PIT) — do NOT let it regress to full-symbolic (the F/G wall).
