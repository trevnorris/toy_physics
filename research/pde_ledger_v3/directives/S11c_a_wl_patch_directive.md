# S11c-a WL engine PATCH — T7 reconciliation, step 2 (Verdict B projection + H.1 coverage + kinematic/flux axis drop)

You are patching the **blind Wolfram (WL) engine** for S11c-a. This is a SCHEMA + CONTROL-COVERAGE patch
driven by the CLOSED spec `directives/S11c_a_SHARED_PHYSICS.md`. Three changes, each grounded in a spec
citation below. ⛔ You are changing WHICH CASES/AXES are computed and emitted — ⛔ you are NOT changing any
already-verified shape-derivative construction, and you must PRESERVE everything in the PRESERVE list.

## ⛔⛔ BLIND ENGINE — non-negotiable
This engine imports nothing and re-derives every object from the spec. This patch cites ONLY the spec. It
carries NO construction, symbol, expression, or result from any sibling engine. Compute every object from the
spec's supplied equations. ⛔ Do not seek out, import, or reconstruct any other engine's output.

## ⛔⛔ RULE 5 — no expected values
This directive states WHAT to compute (which cases, which axes, which controls), ⛔ never what any expression
equals, whether two representatives agree or differ, how many terms survive, or what any residual is. Your job
ENDS at compute-grade-emit. ⛔ Do not compare any output to an expected value or to another engine. ⛔ Do not
assert a residual is zero; PRINT it. The diff happens on the reviewer's side, where a mismatch is a finding,
not a build failure.

## ⭐⭐⭐ THE THREE SCRIPT CLAUSES — carry verbatim, they bind this patch
> **1. The script may PRINT computed objects. It may NOT state conclusions.** Every emitted payload is a CAS
> object (expression, graded value, boolean from a symbolic test), ⛔ never prose describing a result.
> **2. PRINT the residual; do NOT assert it.** Compute -> emit -> (only then) any guard. ⛔ No
> `assert residual == 0`; that writes the expected output into the check.
> **3. Interpretation belongs to the STEP RECORD, not the script.**

⭐ **The added representative cases MUST be REACHED BY COMPUTATION.** When you add the density-representative
axis to the projection cluster (Change 1), each representative's case is computed from THAT representative's
density source (the engine's existing `rho4Profile[density, ...]` / `densityNames` machinery). ⛔ Do NOT
synthesize the second representative by copying the first case under a second key — a hand-copied duplicate is
a fake axis. A value that comes out present-and-identical across the two representatives is a legitimate
COMPUTED invariant; a value that is identical because it was copied is the defect. Compute each independently.

## Baseline, deliverable, how to run
- **Baseline (committed):** `mathematica/S11c_a_interface_geometry_mathematica_audit.wl` at HEAD (WL engine
  `ddecdbc2`, repaired + verified: T-h multiline parse, T-c' two operands, §5a genuine w' route-2 for
  T-c/d/i/g, T-f window-inside-integral, T-0 sigma_W grading). All 40 WL_S11CA_ tags emit; stderr clean.
- **Deliverable:** the SAME file, patched in place.
- **Run it yourself to confirm it still evaluates cleanly** after the patch (Mathematica needs
  `--sandbox danger-full-access`; the full run is ~15-20 min and peaks ~9-10 GB; wrap the kernel run in
  `timeout` generously, run ONE kernel at a time — the licence has two seats). Report the literal stderr and
  the emitted tag inventory. ⛔ Do not read the numeric payloads as pass/fail; only confirm it evaluates and
  the schema/coverage below is present.

## ⛔ PRESERVE — do not touch these (verified physics)
- Every shape-derivative construction (`shapeDerivative`, `projectionShapeDerivative`,
  `derivedFaceCase`, `derivedProjectionCase`, the direct-face / evolution / virtual sources).
- §5a rep-invariance and independence controls, including the genuine material-coordinate w' route-2 for
  T-c/T-d/T-i/T-g and the T-g flattening Jacobian.
- T-h (`sigmaEulerianSource` = rho4*(1+theta)*graphThickness with explicit products, and the whole-file
  multiline-product fix), T-c' two-operand emission, T-0 sigma_W grading, T-f window-inside-Inactive[Integrate].
- The `Inactive[Integrate]` / opaque-integral machinery — projection integrals stay opaque to outer transforms.
- Blindness, the WL_S11CA_ tag grammar, `stripConditional`, and the flush discipline.

Only the case/axis structure and control coverage described below changes.

---

## CHANGE 1 — ADD the density-representative axis to the T-f projection cluster
**Spec.** §1b (lines 74-83): the projection object is S11b's conservation law `∂_tρ_4D + ∇₄·j = 0`,
`j = ρ_4D v_bulk`, integrated against the slab window `Ω`. The density that enters the projection is the bulk
4D density `ρ_4D`. §3a (lines 312-318): `ρ_4D^α ≡ ρ_4D,bg^{0,α}(1+θ)` is representative- and
anchoring-dependent, and §2b (lines 228-230) constructs the two representatives without simplifying one into
the other before T-f/T-g/T-h/T-i. §4 (lines 397-400): every item is computed "for ... both density
representatives wherever density enters." §4 T-f (lines 424-429) enumerates the five projection objects.

**Current state (verified).** `projectionTermsSource[branch_String, dynamic_]` (line 758) builds the
projection from a single abstract bulk density (`rhoBulkZero`/`rhoBulkWave`) with NO representative parameter.
All five T-f emitters key by `branch|DOF_dof` (or `branch|origin|DOF_dof` for origins) with NO density axis:
- `PROJECTION_SHAPE_DERIV` (emit ~line 1090)
- `PROJECTION_DYNAMIC_OPERAND` (~1100)
- `PROJECTION_STATIC_OPERAND` (~1110)
- `PROJECTION_RESIDUAL` (~1119)
- `PROJECTION_TERM_ORIGINS` (~1134)

**Patch.** The projection's bulk density must be the SELECTED representative's `ρ_4D`, drawn from the same
representative machinery the face/evolution objects already use (`rho4Profile[density, ...]`, `densityNames`).
Compute and emit EACH of the five T-f objects for BOTH density representatives `{RHO4_CONSTANT,
RHOBR_CONSTANT}`, and extend each object's key to carry the representative (matching the engine's existing
density-face key convention, e.g. `densityFaceCaseKey` / a `density` token in the key). Because the density
enters the projection (§1b) and §4 mandates both representatives wherever density enters, ALL FIVE T-f objects
carry the representative axis — including `PROJECTION_STATIC_OPERAND`. ⛔ Do not pre-decide that any of the
five is representative-independent; compute each representative's case from that representative's density
source and emit both (⭐ a present-and-identical result is a computed invariant, not a reason to drop the
axis — see the copied-duplicate warning above).

⚠ Contrast with Change 2: the projection uses `ρ_4D` (a representative), so it CARRIES the axis; the
flux/kinematic objects use `ρ_m` (a bound constant, Change 2), so they do NOT.

---

## CHANGE 2 — DROP the density-representative axis from KINEMATIC_BALANCE and RELATIVE_FLUX
**Spec.** §3b (lines 349-353): `J_s^α ≡ ρ_m (v_bulk,s − v_face,s^α)·n̂_s^α`, and
`n̂_s^α·v_bulk,s = V_s^α + J_s^α/ρ_m`. The density in the relative flux and the kinematic balance is `ρ_m`,
the bound brane-density constant (an S11b KNOB, `LEDGER['rho_m']`) — ⛔ NOT one of the density representatives
`ρ_4D`/`ρ_br`. These two objects therefore do not depend on the density-representative axis.

**Current state (verified).** `RELATIVE_FLUX` (emit ~line 1001) and `KINEMATIC_BALANCE` (~1024) both key by
`densityFaceCaseKey[branch, density, sign, dof]` and loop `{density, densityNames}` — a density-representative
axis that §3b shows is spurious for these objects.

**Patch.** Remove the `{density, densityNames}` axis from `RELATIVE_FLUX` and `KINEMATIC_BALANCE`. Key each by
`(branch, face, dof)` only (the engine's `faceCaseKey`). Where the source constructor requires a `density`
argument, pass a single canonical representative (the objects use `ρ_m`, so the choice does not enter the
result) — but the EMITTED object must carry no representative axis. Apply this drop everywhere these two
objects appear: the main emit AND their §5b form-control and §5c uniform-limit coverage (Change 3), so the
control keys stay consistent with the object keys.

---

## CHANGE 3 — EXTEND §5b form-ablation and §5c uniform-limit coverage to EVERY emitted T-object
**Spec.** §5b (lines 487-492): "For every T-object and each direction emit the baseline operand, the
independently recomputed ablated operand, and their residual." §5c (lines 494-499): "For each S11c-a object,
independently obtain the (η,σ_W)→(0,0) operand and the corresponding S11b object. Emit both and their
residual."

**Current state (verified).** Both controls cover 13 objects and handle the projection via
`PROJECTION_SHAPE_DERIV` only. MISSING from BOTH §5b (form, ~lines 1377-1509) and §5c (uniform, ~1537-1659):
- `EVOLUTION_TERM_ORIGINS`
- `PROJECTION_STATIC_OPERAND`
- `PROJECTION_DYNAMIC_OPERAND`
- `PROJECTION_RESIDUAL`
- `PROJECTION_TERM_ORIGINS`

**Patch.** Add §5b form-ablation (base / independently-recomputed ablated / residual, per direction) AND §5c
uniform-limit (S11c-a operand / S11b operand / residual) coverage for those five objects, computed
independently from their sources (⛔ not by reusing the main-emit payloads). Consistency with Changes 1-2:
- The projection objects' §5b/§5c coverage must carry the density-representative axis added in Change 1 (both
  representatives), matching the objects.
- `KINEMATIC_BALANCE` / `RELATIVE_FLUX` §5b/§5c coverage must drop the density-representative axis per Change 2.
- The form ablation for each added object is a genuine per-direction source-form change (as the existing form
  loops do it), independently recomputed — ⛔ never a coefficient rescale, ⛔ never a reuse of the emitted value.

---

## Acceptance — SCHEMA/COVERAGE presence only (⛔ no numeric pass condition, rule 5)
Confirm structurally, ⛔ not by any expected value:
1. The engine still evaluates: full run, literal stderr reported, all T-object tag families still emit.
2. Each of the five T-f projection tags (`PROJECTION_SHAPE_DERIV`, `_STATIC_OPERAND`, `_DYNAMIC_OPERAND`,
   `_RESIDUAL`, `_TERM_ORIGINS`) now has keys that carry a density-representative token for BOTH
   `RHO4_CONSTANT` and `RHOBR_CONSTANT`, and each representative's payload was computed from that
   representative's density source (show the two source paths differ, ⛔ not the values).
3. `KINEMATIC_BALANCE` and `RELATIVE_FLUX` keys carry NO density-representative token.
4. `CONTROL_FORM_*` and `UNIFORM_LIMIT_*` keys now include `EVOLUTION_TERM_ORIGINS`,
   `PROJECTION_STATIC_OPERAND`, `PROJECTION_DYNAMIC_OPERAND`, `PROJECTION_RESIDUAL`, `PROJECTION_TERM_ORIGINS`
   (projection control keys carry the representative axis; kinematic/flux control keys carry none).
5. The PRESERVE list is byte-unchanged in behaviour: the §5a routes, T-h, T-c', T-0, T-f constructions are
   untouched (report the diff is confined to the projection density axis, the flux/kinematic axis drop, and
   the control-coverage additions).
6. No new `assert`-before-emit; no hand-typed payloads; every added case reached by computation.

## Ops
- Codex at xhigh, background, `--sandbox danger-full-access` (Mathematica).
- Wrap every kernel run in a generous `timeout`; ONE kernel at a time (two-seat licence); watch RSS — a run
  over ~15 GB or a leftover kernel after the job ends is a build-quality finding (kill by PID, check `free -h`).
- Write your transcript OUTSIDE the repo. Patch only the one `.wl` file.
