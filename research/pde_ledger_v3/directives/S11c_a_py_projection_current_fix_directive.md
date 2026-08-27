# Fix directive (rev 2) — S11c-a PY dynamic-window projection: the perturbation current is frozen in w

## Object to correct
`research/pde_ledger_v3/scripts/S11c_a_interface_geometry_sympy_audit.py`:
- `projection_terms` (~lines 1114–1171), consumed by `build_projection_raw`, emitting the **five** tags
  `S11CA_PROJECTION_{SHAPE_DERIV, STATIC_OPERAND, DYNAMIC_OPERAND, RESIDUAL, TERM_ORIGINS}` — each keyed by
  `(branch, dof, representative)`; there is **no** face axis in the projection key;
- `uniform_projection_reference` (~lines 1659–1682), the independent reconstruction that supplies every
  `PROJECTION_*` §5c uniform-limit reference.

## The defect
In `projection_terms` the perturbation bulk current enters only through the `w`-independent symbols
`j_bulk[i]` (= `delta_j_bulk_{i}`) and the in-plane jets `grad_j_bulk[i][i]`. The normal term is
`WINDOW_NORMAL_CURRENT = -ε·∫ j_bulk[3]·window_normal dw`, where `window_normal` is the background window's
normal derivative `∂_wΩ`. Because `j_bulk[3]` is a SymPy `Symbol` constant in `w`, the projected conservation
law silently imposes `∂_w δj_w = 0`: the perturbation current's normal variation is absent. The normal-jet
symbols `dw_delta_j_bulk` are declared but used only in the face-trace construction, never in the projection.
`uniform_projection_reference` repeats the same `w`-constant current.

## Governing law (supplied, §1b — unfalsifiable here)
`∂_tρ_4D + ∇₄·j = 0`, `j = ρ_4D v_bulk`; the projection integrates this against `Ω` and integrates by parts
in `w`. `∇₄·δj` includes the normal-derivative term `∂_w δj_w`. §3c makes the **background** current
`ρ_4D⁰v_bulk⁰` vanish — and that zero stays zero under the integral — and governs **traced** fields at the
face; it sets **no** condition on the **perturbation** `δj` or its normal derivative `∂_w δj_w`.

## The correction — name the object, not the recipe
The perturbation current's **normal variation `∂_w δj_w`** must be present in the projected conservation law,
so the §1b integral of `∇₄·δj` is complete rather than truncated to its in-plane and time parts. How the
current's `w`-dependence is represented in the projection integrand is the engine's to **compute from the
supplied conservation law**. ⛔ Do **not** reuse the face-trace mechanism (`affine_bulk_perturbation` /
`dw_delta_j_bulk[face]`): those are keyed to a single face and expand about `±W₀/2`, while the projection is
face-less and integrates once over `w`; a single-face choice, a two-face sum, or an interpolation across the
slab is not authorized by §1b/§3c.

## Binding constraints
1. **No double count across the IBP.** The projection is already in post-IBP form (its normal term is
   `-∫ δj_w·∂_wΩ`). Correct this **either** by (a) rebuilding the normal contribution from the *pre-IBP*
   conservation law `∫ Ω·∂_w δj_w` and re-deriving the window-normal origin, **or** (b) letting `δj_w` carry
   its normal variation *inside the existing post-IBP normal term* — but ⛔ **never both**. Do not add a
   second `∫ Ω·∂_w δj_w` channel alongside `WINDOW_NORMAL_CURRENT`; the two are the same contribution.
2. **Only the normal divergence is missing.** The in-plane divergence `∂_i δj_i` (via `grad_j_bulk`) and the
   time term `∂_t δρ` are already present and correct; the correction concerns `∂_w δj_w`. ⛔ Do not introduce
   undeclared mixed jets such as `∂_i∂_w δj_i`; if a candidate representation needs symbols the engine does
   not supply, it is the wrong representation.
3. **Scope.** Apply the correction to **both** the `dynamic=True` and `dynamic=False` branches of
   `projection_terms`, and to `uniform_projection_reference` (§5c requires it to remain an *independently*
   computed reference — keep it independent, built on the same corrected current law). Form-ablation reuses
   `build_projection_raw`, so it propagates automatically; do not special-case it.

## Guards
- ⛔ Do not touch the trace construction (line 651 is correct), the (zero) background current, the window
  definition, or any other T-object.
- ⛔ Do not add, remove, or rename a `PROJECTION_*` tag, and ⛔ do not add a face axis to the projection key.
- The correction re-enters the computation at the supplied conservation law / the current's construction,
  never at a result (structural rule, §6).

## Discipline
The engine's emission rules are unchanged: it **prints** operands and the residual and asserts no conclusion;
a nonzero residual emits and the run exits 0. Whether, and which, normal-variation terms survive is
**computed** by the engine, not asserted here, and is not stated in this directive. §§1–3 are supplied and
unfalsifiable in this build.

## Builder report
Under 25 lines: the functions/lines changed, the five tags re-emitted, runtime, and any ambiguity in applying
§1b to the projection current. State that §§1–3 were supplied and unfalsifiable in this build.
