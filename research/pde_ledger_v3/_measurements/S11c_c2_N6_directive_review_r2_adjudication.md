# N6 build-directive review — round 2 adjudication (2026-09-06)

Corrected directive + §5c (`d0666953`) → 2 decision legs r2 (Codex-sol + Grok, review-until-clear). Reports:
`scratchpad/{codex,grok}_N6_dirreview_r2.txt`.

## Ordering fix VALIDATED (the r1 spec correction holds, ⛔ not re-litigated)
**Both legs independently confirmed map-then-extract** from the code: `extract` (`selfenergy_fold:325-342`) is a lossy
Eulerian projection (Helmholtz `restricted`, div/curl, fixed Eulerian tests); `material_pullback` (`brane:1941`) mixes
sectors; ⇒ `extract ∘ T ≠ T ∘ extract`; extract-then-map drops the advective off-diagonal; the parent does
map-then-extract; there is no pre-map `SLAB_M` callable. ✅ §5c route-2 ORDER is correct.

## Findings (converged; verified by orchestrator reading — rule 13)
### F1 (both legs) — identifying `T` with `material_pullback` applied to the carrier ROWS is type-wrong + N6-breaking
Three code-grounded failure modes (Grok): (1) **Annihilation** — `material_pullback` ends in a **2nd scale-derivative
over `WAVE_SYMBOLS`** (`brane:1970-1981`) that DROPS the linear carrier rows ⇒ `I_{M→E}≈0`, `R_N6≈I_E` (spurious).
(2) **Wrong domain** — `material_pullback` maps a **bulk scalar density**; the c2 increment is the **δp-slot FACE
carrier** (traction / `closure_shape_deriv`), which it never touches. (3) **Incomplete adjoint** — "T on trials,
extract supplies Eulerian tests" is exactly the trial-only `representation_pullback` shape the directive forbids; the
**reverse u-row channel** (`δθ_M = δθ_E + δu·∇ρ/ρ`) only appears if the field redefinition acts on the **action /
virtual-work SCALAR** and Eulerian variation is taken AFTER — that reverse off-diagonal is the channel N6 tests.
- **Converged fix (both legs): build route 2 as a native construction, ⛔ NOT `material_pullback` on rows** — native
  material δp-slot carrier from **S11c-a MATERIAL face sources** (traction / `closure_shape_deriv` / face velocity
  `V_s`) folded into the same δp slots via the S11c-b **face-fold** (a NAMED existing construction, ⛔ not a new
  recipe); `close` = the SAME imported same-`α` c1 response with route-2 `V_s` = material face velocity; any residual
  N4 field redefinition applied at the **action/scalar level → Eulerian variation → single `extract`**. Account for the
  S11c-a material face builder ALREADY mapping the covector to Eulerian (inverse-transpose, `interface:742-755`) — ⛔ no
  double-map.
- **DISPOSITION (user-approved 2026-09-06 — delegate the typed pipeline to astra):** the exact route-2 CAS pipeline is
  the INSTRUMENT — astra's to design in the build; the directive hands the **object + verified refs + a
  definition-of-done** (the converged constraints above + the open design question), and the **build legs gate the
  construction**. ⛔ The orchestrator stops hand-authoring the route-2 recipe (2 rounds wrong). Open design question for
  astra: whether a bulk N4 field map is still required, or native material faces + material `V_s` already yield the
  material increment.

### F2 (both legs) — the N4-probe structural-absence names the WRONG representative
`material_pullback`'s advection is `u·∇ρ₄/ρ₄`. `RHO4_CONSTANT` (`ρ₄=ρ_br/W₀`, `brane:1885-1890`) ⇒ `∇ρ₄=0` ⇒ term
**structurally ABSENT**. `RHOBR_CONSTANT` (`ρ₄=ρ_br/W_bg`) ⇒ `∇ρ₄≠0` via live `W_bg` ⇒ term PRESENT (the live probe).
My directive/§5c named `RHOBR_CONSTANT` as the absent one — **reversed**; also §5c's `u·∇Σ_E⁰` is a DIFFERENT toggle.
- **Fix:** the probe omits the map's advective term (`u·∇ρ⁰/ρ`); structural absence is **RHO4_CONSTANT**; delete the
  RHOBR_CONSTANT example; adjudicate the absence from emitted provenance (⛔ no A−A).

### F3 (Codex) — the Eulerian tilt control lacks a same-construction baseline
Route 1 uses the IMPORTED (already-materialized) Eulerian carrier; the c2 slope override alters the c1 DtN kernel jet
(`selfenergy_fold:384-387`), ⛔ not the slab carrier's normal-derived coefficients; the genuine source mutation lives
in the S11c-a face constructor (`interface:818-845`). Comparing an imported baseline vs an independently-reconstructed
corrupted operand conflates reconstruction drift with the mutation.
- **Fix:** commission an **injectable Eulerian carrier factory**; produce both uncorrupted + slope-corrupted through
  it; **PIT-check the uncorrupted factory output against the imported carrier** before using the control.

## Passed (not findings, both legs): audit edit + orientation, leak discipline, dimensions, tractability shape.

## Disposition
Fold: (§5c) drop the `T=material_pullback` identification, state the object + map-then-extract principle + native-face
construction constraints + reverse-channel requirement, fix the advection representative. (directive) restructure
route 2 to object + verified refs + DoD + explicit delegation to astra; fix F2 + F3. Then re-review §5c+directive r3
(Codex-sol + Grok, review-until-clear). Prior baseline `d0666953`.
