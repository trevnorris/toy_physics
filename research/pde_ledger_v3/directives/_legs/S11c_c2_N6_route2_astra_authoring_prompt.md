# Author the route-2 construction spec — S11c-c2 N6 representation-invariance control

You are astra. **Author** the correct route-2 construction for the c2 N6 representation-invariance control, anchored to
the parent pattern. Working dir `/var/projects/toy_physics`; paths repo-root-relative unless absolute. Read the code;
reason from it; you may run small SymPy checks to verify your own claims. **Write your spec to
`research/pde_ledger_v3/_measurements/S11c_c2_N6_route2_spec_astra.md`** and summarize it in your final message.

## Why you are authoring this (context — the orchestrator's 3 attempts were type-wrong)
The N6 control compares the c2 **self-energy increment** built two ways at a FIXED anchoring `α` and density `ρ`:
Eulerian vs material-coordinate, which must agree (representation invariance; residual is a computed finding, ⛔ NO
target value — never make residual-zero an exit/assert). The orchestrator wrote route 2 three times as
`extract( T_{M→E}[ close(SLAB_M) − SLAB_M ] )` with `T` a field-redefinition applied to the increment — **all three
were the wrong TYPE.** Two review legs established why (verified in code — treat as given, but re-check anything you
doubt):
- The increment is **δp-slot-only** (`extract(close − open)` is linear; bulk/kinetic cancels). `S_P` = coefficients of
  `{δp_s, ∂_w δp_s}`; it carries **no field `θ`**, so a bulk θ-shift on the carrier is a no-op/annihilating.
- The face δp-rows get their θ / N4 dependence through **μ_θ**, which is **bound into the virtual work BEFORE
  differentiation**: `face_generalized_force_rows` (`scripts/S11c_b_brane_operator_sympy_audit.py:2135-2179`) does
  `bind_mu_theta_operand(work, branch, mu_theta_amplitude)` then `diff(density, delta_v_u/delta_v_e_W)`.
- `material_pullback` (`scripts/S11c_b_brane_operator_sympy_audit.py:1916-1981`) is a **bulk-density quadratic** map
  (2nd scale-derivative annihilates linear rows; never writes `delta_p_*`). ⛔ It is NOT the route-2 map.

## The parent template you MUST anchor to (this is the correct face N6, no `T`)
`task_rep_invariance` (`scripts/S11c_a_interface_geometry_sympy_audit.py:1481-1496`) builds each face/geometry quantity
(`RELATIVE_FLUX, TRACTION, CLOSURE_SHAPE_DERIV, VIRTUAL_WORK_SHAPE_DERIV, VIRTUAL_CONSTRAINT, EVOLUTION_MASS_BALANCE`)
Eulerian (`RAW_PRIMARY`) vs `route="MATERIAL"` (`build_geometry_quantity(..., route="MATERIAL")` etc.) and **differences
them directly** (`object_difference`) — NO field-redefinition `T`, NO bulk Jacobian on the difference. The
material→Eulerian map lives INSIDE the material builders (covector inverse-transpose, `:742-755`; material `V_s`,
`:794`). Read `build_geometry_quantity`, `build_virtual_constraint_raw`, `build_evolution_raw`, and the `route`
plumbing.

## What to author (route 2 = the parent pattern applied to the c2 face-slot increment)
```text
route 2 (material):  I_{M→E}^{α,ρ} = extract( close(SLAB_M) − SLAB_M )        (no separate T on the increment)
```
Author, as equations grounded in named functions:
1. **`SLAB_M`** — the c2 δp-slot face carrier built with **native material-route face ingredients** (traction /
   `closure_shape_deriv` / virtual work / `V_s` via the S11c-a `route="MATERIAL"` builders, already Eulerian-mapped)
   folded into the **same δp symbols** the S11c-b face-fold uses. Say exactly which functions/quantities, and confirm
   the δp symbols coincide with route 1's so the difference `I_E − I_{M→E}` is meaningful.
2. **The μ_θ binding — DECIDE and JUSTIFY it (this is the load-bearing spec decision the reviews flagged).** Does route
   2 bind the imported **Eulerian** `mu_theta_operator`, or a **material** μ_θ (from `operator_from_density` on the
   material-route bulk)? Derive the answer from what "representation invariance of the increment" MEANS + the parent
   pattern. State the consequence: if material μ_θ, the N4 advection channel (`u·∇ρ₄/ρ₄`) is present and the advection
   probe is live on `RHOBR_CONSTANT`; if Eulerian μ_θ, the channel is structurally absent from the increment for both
   density reps and the advection probe emits that absence. ⛔ Do not leave this open.
3. **`close`** — the same imported **same-`(α,ρ)`** c1 response `s11c_c1_face_response`, opaque (⛔ not re-mapped), with
   route-2 `V_s` = material face velocity.
4. **The controls** — (a) **advection probe**: name the exact term omitted and inside which object (per your μ_θ
   decision), and which representative has it structurally absent (`RHO4_CONSTANT`: `ρ₄=ρ_br/W₀` ⇒ `∇ρ₄=0`;
   `RHOBR_CONSTANT`: `ρ₄=ρ_br/W_bg`, live) — ⛔ emit computed absence, never `A−A`. (b) **tilt probe**: an injectable
   Eulerian carrier factory, PIT-checked against the imported carrier (⛔ not the c1-DtN-jet override
   `selfenergy_fold:384-387`). (c) Address the r3 concern that the "reverse u-row channel" must NOT be satisfied
   vacuously by `extract`'s existing U-row curl (`selfenergy_fold:337-341`, present without N4) — say what genuinely
   distinguishes the routes.
5. **Independence** — the two routes must be structurally different so `I_E − I_{M→E}` is not zero by construction;
   state the one-sided corruption that proves it.

## Constraints (carry into the spec)
- **Tractability:** carrier-first (`S_P = S − S|_{P=0}`), ⛔ no `build_case` end-to-end, evaluate/truncate BEFORE closure
  expansion; residual + slot-linearity guard + controls tested by **exact finite-field PIT** (several primes, shared
  generic samples across routes, denominators cleared, singular/branch cells handled, FN bound). ⛔ never full-symbolic
  over all grades on these operands (the F/G wall); ⛔ never `k_in=k_out`; ⛔ never evaluate the Fourier/DtN integral;
  same `Z`/resolvent both routes.
- **Emit discipline:** PRINT operands then residual; ⛔ never assert; residual-zero ⛔ never a builder exit; compact
  numeric, ⛔ no giant symbolic; `[L,T,M]` + able-to-fail dimensional consistency + `(ε,η,σ_W)` multigrade on every
  object.
- ⛔ WITHHOLD: there is no expected value for the residual — do not invent one; the diff is adjudicated externally.

## Output
Write `research/pde_ledger_v3/_measurements/S11c_c2_N6_route2_spec_astra.md`: the route-2 construction (formula + the
μ_θ decision with justification + the fold recipe in named functions + the controls + independence + tractability),
precise enough to fold into §5c and a build directive. Flag any place the code does not yet support the construction
(a missing material builder, a δp-symbol mismatch) as a computed finding. This spec will be reviewed (fresh Claude +
Grok) before use — make it correct and evidence-grounded, not merely plausible.
