# N6 build-directive review — round 1 adjudication (2026-09-06)

Directive `directives/S11c_c2_N6_sympy_build_directive.md` (orchestrator-written) → 2 decision legs (Codex-sol xhigh
+ Grok, review-until-clear). Both returned **NOT CLEAR TO BUILD**. Reports:
`scratchpad/{codex,grok}_N6_dirreview.txt`. (Grok's first run read all sources then hit a transient xAI 500 at
capacity with no verdict; re-run exit 0.)

## The one genuine leg DISAGREEMENT (M1: the disagreement is the measurement) — RESOLVED on evidence
**Contested:** is §5c's route-2 ordering `I_{M→E} = T_{M→E}[ extract(close(SLAB_M) − SLAB_M) ]` (extract-then-map)
correct?
- **Codex:** NO — spec defect. The N4 map mixes sectors, so extract-then-map drops the advective off-diagonal;
  correct order is **map-then-extract** `extract( T[close(SLAB_M) − SLAB_M] )`.
- **Grok:** §5c is fine; the defect is purely the directive citing `build_operator("MATERIAL")` (which already maps).
  Fix at the directive level: `SLAB_M` = unmapped native carrier, `T` = full congruence applied once to the weak
  increment.

**Decisive fact (verified by orchestrator reading — rule 13, mechanical code-reading, no CAS authored):** the §3c
`extract` (`scripts/S11c_c2_selfenergy_fold_sympy_audit.py:325-342`) is a **lossy, intrinsically-Eulerian
projection** — it applies `restricted(·,'TRANSVERSE'/'LONGITUDINAL'/'THETA')` (Helmholtz split), takes div/curl, and
pairs against **fixed Eulerian test functions** (`s11cc2TestTheta/E/Phi/A` on the Eulerian coord `X`). It does NOT
retain the full operator. The N4 map mixes sectors (`material_pullback:1941`: `θ_material = θ + dot(u,∇ρ₄)/ρ₄`). ⇒
`extract ∘ T ≠ T ∘ extract`. Grok's own finding 2 concedes "Eulerian extract on unmapped material rows makes
`T*·O·T` ill-posed" — which is itself the argument for map-then-extract. The **parent house pattern is
map-then-extract**: `build_operator(route="MATERIAL")` (`S11c_b_brane_operator_sympy_audit.py:2762-2783`) calls
`material_pullback` on the density THEN `operator_from_density` extracts.

**Adjudication: Codex is correct — §5c's route-2 ordering is a spec defect.** Corrected object (map-then-extract),
which SUBSUMES Grok's directive-level fixes:
```text
route 1 (Eulerian):    I_E     = extract( close(SLAB) − SLAB )
route 2 (material→E):   I_{M→E} = extract( T_{M→E}[ close(SLAB_M) − SLAB_M ] )
R_N6 = I_E − I_{M→E}   (computed finding, ⛔ no target)
```
- `SLAB_M` = the **unmapped native** material-coordinate c2 pressure-slot carrier (S11c-a face-flattening geometry +
  material face sources), BEFORE the map. (Grok's fix.)
- `T_{M→E}` = `material_pullback`'s field redefinition (`θ→θ+u·∇ρ⁰/ρ`, anchoring-branched `e_W` shift) × Jacobian
  `1+tr(∇u)`, applied once to the carrier DENSITY, as an explicit equation. (Both legs' fix — cite `:1942-1981` as
  the definition of T, ⛔ never call `build_operator(route="MATERIAL")` wholesale.)
- a single Eulerian `extract` applied AFTER the map on BOTH routes ⇒ no ill-posed `extract_M`. (Codex's fix +
  dissolves Grok's finding-2 worry.)

This is orthogonal to, and stacks on, the already-cleared axis correction (`30d4b72d`: Eulerian-vs-material at a
FIXED anchoring, not cross-anchoring). It refines the route-2 FORMULA within that correct axis.

## The other verified findings (all fold into design + directive; none contest §5c except finding 1)
- **Codex-2 / Grok-1,2 (no pre-map callable):** `build_operator("MATERIAL")` already maps + supplies Eulerian face
  slots; there is no callable returning an unmapped `SLAB_M`. ⇒ build the native unmapped carrier from S11c-a
  `build_material_face_source` (normals/measure/velocity/traces) + open-slot coefficients; `T` = `material_pullback`
  as the explicit density map. CONFIRMED (code read).
- **Codex-3 (carrier-first vs `build_case`):** route 1 must NOT call `build_case` end-to-end (it closes+truncates+
  extracts full rows before returning — recreates the wall); reuse only low-level binding/face-response/extraction
  helpers, carrier-split + specialize/truncate BEFORE closure expansion, both routes. CONFIRMED.
- **Grok-3 (slot-linearity guard):** the guard must be the SAME finite-field PIT (shared samples, FN bound), ⛔ not a
  second full-symbolic increment (that IS the wall). CONFIRMED.
- **Codex-4 (audit correction):** require an audit EDIT (not companion-only): replace `:1087-1107` with the RAW
  unmapped `lab['increment'] − material['increment']` under `ANCHORING_L_MINUS_M` (§5c:366 orientation LAB−MATERIAL,
  ⛔ no `representation_pullback`, ⛔ old sign); keep the density loop `:1085-1086`. CONFIRMED (§5c:366 read).
- **Codex-5a (leak):** ⛔ do NOT point the builder at the design/adjudication docs (they state `R_N6 "must vanish"`);
  the directive is the sole builder-facing artifact; remove the directive's own "nonzero expected" anchoring-sign
  statement. CONFIRMED (leak).
- **Codex-5b (dimensions):** add the §6 contract — `[L,T,M]` restored + able-to-fail dimensional consistency +
  `(ε,η,σ_W)` multigrade on every emitted object (`N12`). CONFIRMED (§6:394-400).

## Disposition
Route as a **spec correction (§5c route-2) + design + directive rewrite**, then **re-review §5c+directive together**
(Codex-sol + Grok, review-until-clear). No commit of the directive until the re-review clears. The current §5c
baseline is committed at `30d4b72d`.
