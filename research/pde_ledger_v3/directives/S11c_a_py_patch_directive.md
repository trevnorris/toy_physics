# S11c-a SymPy engine PATCH — T7 reconciliation, step 2 (Verdict C full virtual-work grid + BG boundary loads + density-map branch drop)

You are patching the **SymPy engine** for S11c-a. Three changes, each grounded in the CLOSED spec
`directives/S11c_a_SHARED_PHYSICS.md`. ⛔ You are changing WHICH CASES are computed/emitted for three objects —
you are NOT changing any already-verified shape-derivative construction. PRESERVE everything else.

## ⛔⛔ RULE 5 — no expected values
This directive states WHAT to compute (which cases/axes), ⛔ never what any expression equals, which pairings
"should" survive, whether two cases agree, or what any residual is. Your job ENDS at compute-grade-emit.
⛔ Do not compare any output to an expected value or to another engine. ⛔ Do not assert a residual is zero;
PRINT it. The diff happens on the reviewer's side.

## ⭐⭐⭐ THE THREE SCRIPT CLAUSES — carry verbatim, they bind this patch
> **1. PRINT computed objects; do NOT state conclusions.** Every emitted payload is a CAS object, ⛔ never prose.
> **2. PRINT the residual; do NOT assert it.** Compute -> emit -> (only then) any guard. ⛔ no `assert ==0`.
> **3. Interpretation belongs to the STEP RECORD, not the script.**
⭐ Every added case MUST be REACHED BY COMPUTATION from the action/ansatz. ⛔ Do not hand-type a case or
synthesize one case by copying another. The ONLY place physical symbols are combined by hand is the action and
the ansatz; every other expression is reached by computation.

## Baseline, deliverable, how to run
- **Baseline (committed):** `scripts/S11c_a_interface_geometry_sympy_audit.py` at HEAD (SymPy engine `9b6438fa`,
  verified: T-0…T-i correct, §5a genuine material-coordinate route-2, T-c′ two operands, ρ_m bound).
- **Deliverable:** the SAME file, patched in place.
- **Run it yourself** (`python3 scripts/S11c_a_interface_geometry_sympy_audit.py > /tmp/py_out.txt`, ~3 min) to
  confirm it still evaluates and emits the schema below. ⛔ Do not read the numeric payloads as pass/fail.

## ⛔ PRESERVE — do not touch (verified physics)
- Every shape-derivative construction: `shape()`, `build_face_source`, `build_material_face_source`, the
  face/evolution/projection/closure/traction builders, `density_pair`, `branch_profile`.
- §5a rep-invariance / independence routes (EULERIAN vs genuine MATERIAL w′ route-2), T-h evolution, T-c′
  two-operand kinematic, T-0 grading, T-f window-inside-`sp.Integral`.
- `ρ_m` binding to `LEDGER['rho_m']` — ⛔ never re-originate.
Only the three objects below change their case structure.

---

## CHANGE 1 (Verdict C) — compute+emit the FULL physical×virtual virtual-work grid
**Spec.** §4 T-d (lines 417-419): "compute the traction object and the shape derivative of `δ_v𝒲_bulk^α`,
taking `δ_vx_s^α` only from the applicable `R_s^α` … **Which virtual-displacement pairings occur is part of
the computation.**" The engine must COMPUTE all physical×virtual pairings and let the result show which occur.

**Current state (verified).** `virtual_work_cases(branch, dof, representative)` (line 854) uses a SINGLE `dof`
that sets BOTH (a) the physical interface perturbation — through `build_face_source`: `dof`→`dh`→`normal_exact`
/`measure_exact` and the traction — AND (b) the virtual displacement (`source.virtual_displacement`, contracted
against the traction). The emit loop (lines 919-924) keys `(branch, dof, representative)` — the DIAGONAL only
(physical DOF = virtual DOF). Emitted via `task_td` (line 1269, `emit_primary("VIRTUAL_WORK_SHAPE_DERIV", …,
export=True)`).

**Patch.** Compute `δ_v𝒲_bulk^α` for EVERY (physical DOF, virtual DOF) pairing:
- the **physical DOF** drives the shape-derivative and the traction (the face source's geometry, `normal_exact`,
  `measure_exact`, pressure, and the resulting `exact_traction`);
- the **virtual DOF** drives the virtual displacement `δ_v` (the `virtual_displacement` that `build_face_source`
  builds for that DOF) contracted against the physical-DOF traction.
Emit every pairing, keyed `(branch, physical_dof, virtual_dof, representative)`. ⛔ Do not pre-select the
diagonal and ⛔ do not filter which pairings to emit — compute+emit them all; the emitted result shows which
occur (§4 T-d:419). The §5 controls that consume the virtual-work producer (form ablation, uniform limit,
rep-invariance/independence) must carry the SAME full grid — they call the same `build_geometry_quantity`
path, so update it consistently. ⚠ `TRACTION` (line 1268) is a separate object and is NOT gridded — leave it.

⚠ **Export mutation (deliberate).** `VIRTUAL_WORK_SHAPE_DERIV` is `export=True`; expanding its case map changes
the regenerated `scripts/S11c_a_exports.py` payload. That is intended. ⛔ Do not commit. Leave the regenerated
export in the working tree; separate the genuine payload change from cosmetic `Dummy`-index counter churn.

---

## CHANGE 2 (Verdict BG) — add the supplied boundary loads to BACKGROUND_STATE
**Spec.** §2d (lines 246-272): `𝔅⁰ ≡ {W_bg, μ_R,bg, ρ_4D,bg⁰, ρ_br,bg⁰, θ⁰, V_s⁰, J_s⁰, 𝒜_s⁰, boundary
loads}` — the background state INCLUDES the boundary loads. The supplied support bundle is
`𝒮_hold⁰ ≡ {f_hold⁰(x), t_hold,s⁰(x)}`. These are SUPPLIED objects, ⛔ not computed.

**Current state (verified).** `emit_supplied_objects` (lines ~1315-1327) builds `state` per (branch,
representative) carrying `profile_context()`, `θ⁰=0`, `V_0_*=0`, `J_0_*=0`, `A_0_*=0` — but NO boundary loads.
The supplied loads `support_body = f_hold_0`, `support_face_plus = t_hold_plus_0`, `support_face_minus =
t_hold_minus_0` (defined lines 174-176) appear ONLY in `ADMISSIBILITY_PREMISE` (line 1331). `BACKGROUND_STATE`
is emitted by plain `emit` (line 1335), NOT `emit_primary` — so it is NOT an export object.

**Patch.** Add the supplied boundary-load fields (`f_hold_0`, `t_hold_plus_0`, `t_hold_minus_0`) to
`BACKGROUND_STATE`, as supplied objects (§2d:251). ⛔ Do NOT duplicate the `θ⁰/V/J/A` zeros it already emits.
Leave `ADMISSIBILITY_PREMISE` as is (its support-stabilisation declaration is a separate role). Update the
`BACKGROUND_STATE` dimension entry (the dimensions map near line 1649) so its tuple matches the new field
count. ⚠ This object is not exported, so no export churn results.

---

## CHANGE 3 (density-map branch drop) — BACKGROUND_DENSITY_MAP keyed by representative only
**Spec.** §2b (lines 228-230): "For each representative, construct `Σ_E⁰(y) ≡ ρ_4D,bg⁰(y)W_bg(y)` … Emit the
**two** computed maps under `S11CA_BACKGROUND_DENSITY_MAP`." The two maps are per REPRESENTATIVE, on the
pre-anchoring `y`; the map does not carry an anchoring-branch axis.

**Current state (verified).** `build_background_density_raw` (line 897) loops `for branch in BRANCHES: for
representative in DENSITY_REPS` and keys `(branch, representative)` — four cases, where `sigma_e =
rho4_value * W_bg` is branch-independent (`density_pair(representative, W_bg)`). Emitted via `task_t0` (line
1232, `emit_primary("BACKGROUND_DENSITY_MAP", …, export=True)`).

**Patch.** Drop the branch axis: key `BACKGROUND_DENSITY_MAP` by `(representative,)` only (two cases). Update
the matching dimensions construction in `task_t0` (line ~1230) and the uniform reference in
`independent_uniform_reference` (lines 1611-1615) so their keys match `(representative,)`. ⚠ Export mutation
(deliberate): this object is `export=True`; the payload shrinks from four to two cases. ⛔ Do not commit; leave
it regenerated; separate from `Dummy` churn.

---

## Acceptance — SCHEMA/COVERAGE presence only (⛔ no numeric pass condition, rule 5)
Confirm structurally, ⛔ not by any expected value:
1. The engine still evaluates (~3 min); all tag families still emit.
2. `VIRTUAL_WORK_SHAPE_DERIV` keys now range over every (branch, physical_dof, virtual_dof, representative)
   pairing (both DOFs independent), and every case was reached by computation (physical DOF drives the
   traction; virtual DOF drives `δ_v`). The §5 controls over virtual work carry the same full grid.
3. `BACKGROUND_STATE` now carries the supplied boundary-load fields, and does NOT duplicate the existing
   `θ⁰/V/J/A` zeros; its dimension tuple matches.
4. `BACKGROUND_DENSITY_MAP` keys carry the representative only (no branch axis); its dimensions and uniform
   reference match.
5. The PRESERVE list is behaviourally unchanged (report that the diff is confined to those three objects, plus
   the deliberate export regeneration and cosmetic `Dummy` churn).
6. No new `assert`-before-emit; no hand-typed or copied case.

## Ops
- Codex at xhigh, background. Write your transcript OUTSIDE the repo. Patch only the one `.py` file (the
  regenerated `S11c_a_exports.py` is an expected side effect — ⛔ do not commit either file).
