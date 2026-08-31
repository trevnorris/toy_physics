# S11c-b #88 — BLAST-RADIUS INSTRUMENT (build directive)

## Role
You are building a **diagnostic CAS instrument** (SymPy). It computes, per strong
DOF operator row, the retained-grade difference between (i) the physically-correct
variable-coefficient brane operator on the **completed** §3a energy basis and (ii)
the operator the current engine emits (which freezes the background spurion at first
order). It **computes and prints objects. It states no conclusion.** The point is to
expose whether completing the §3a basis changes each operator row; the *interpretation*
(which cross-engine family verdicts, if any, are affected) is done downstream by the
orchestrator, **not** by this script.

Deliverable: `research/pde_ledger_v3/scripts/S11c_b_88_blast_radius.py`, runnable as
`python3 research/pde_ledger_v3/scripts/S11c_b_88_blast_radius.py`, printing the tagged
objects below and exiting 0.

⛔ **Independent construction.** Build the corrected operator yourself from the physics
stated here and the **committed** engine patterns cited by file:line. Do **not** search
for, read, or import any scratch/probe script outside the repository. Do **not** copy
an existing "blast radius" script — there is none in the tree; if you find one under
`/tmp` or a home dir, it is not yours to read.

---

## The physics (supplied — this is the §1d/§3a spec content, already settled)

The §3a stored-energy basis on the **non-uniform** background contains, per background
source `S ∈ {W_BG, MU_R_BG}`, a family of "spurion" invariants of the form
`g_a · B(DOF, ∂DOF)`, where `g_a = ∂_a W_bg` (source `W_BG`) or `∂_a μ_R,bg` (source
`MU_R_BG`) is the **first background jet**, and `B` is a DOF bilinear. The engine
enumerates these as `FIRST_JET_CONTRACTION_01..15` (15 per source) in
`enumerate_new_candidates` (`scripts/S11c_b_brane_operator_sympy_audit.py:1073`).

Two spec facts you must honour (both already established; do not re-litigate):

1. **The corrected §3a basis retains all 15 candidates per source** (they are linearly
   independent under the variable-coefficient quotient — nullity 0). The current engine
   keeps only 8/source because its **basis-independence quotient**
   (`basis_euler_signatures` L936 + `quotient_independent_indices` L972) builds its
   total-divergence map from the DOF fields **only** (`basis_fields` L1025) and never
   differentiates the background jet `g_a` — so it treats `g_a` as constant and
   spuriously identifies 7 candidates as redundant. You do **not** need to re-derive the
   count; you are given: **completed basis = all 15/source**, each carrying its own
   independent free coefficient.

2. **The correct Euler–Lagrange operator differentiates the background spurion.** For an
   energy term `c · g_a · B`, the true EL operator's total-divergence step
   `−Σ_i ∂_i( ∂(term)/∂(∂_i φ) )` must differentiate `g_a` too, producing the **second
   background derivative** (Hessian) `∂_i g_a = ∂_i∂_a W_bg`. In the engine's grade
   bookkeeping this second derivative is
   `grad_W[j] ↦ sigma_W * w1_profile_d{·}d{·} / L_W` and
   `grad_mu[j] ↦ mu_R * sigma_W * m1_profile_d{·}d{·} / (W0 * L_W)` — **exactly the map
   entries already installed** in the committed routines `operator_dx`
   (`...sympy_audit.py:1850-1874`, specifically the two `derivative_map[grad_W[...]]` /
   `derivative_map[grad_mu[...]]` assignments) and `background_dx`
   (`...sympy_audit.py:2122-2137`). Reuse that pattern to build a **Hessian-retaining
   `dx`**.

The current strong-row operator does **not** do this: `operator_from_density`
(`...sympy_audit.py:1459-1511`) uses the global `dx` (`...sympy_audit.py:616`) whose
`DERIVATIVE_MAP` stops at the first background jet (`...sympy_audit.py:611-613`), so the
Hessian is frozen out of the emitted strong rows.

---

## Objects to reuse (committed engine machinery — import the audit module)

Import `S11c_b_brane_operator_sympy_audit` as the engine `M` (add its dir to
`sys.path`). Reuse, verbatim:
- `M.enumerate_new_candidates(g_vector)` — the 15 candidates/source (`:1073`).
- `M.basis_euler_signatures` + `M.quotient_independent_indices` — to obtain the engine's
  **selected(8)** and **omitted(7)** index split per source (`:936`, `:972`); you need
  the split only to label drivers, not to decide the basis.
- `M.construct_energy("LAB_HELD")` — the emitted (frozen) energy `EnergyBuild`; its
  `.density` is the frozen 8/source density (`:1228`).
- `M.live_basis_substitution(actual_vector)` — maps abstract candidate symbols to live
  fields (`:1178`); `actual_vector = M.grad_W` for `W_BG`, `M.grad_mu` for `MU_R_BG`.
- `M.first_shape_series` — the retained-grade truncation: substitute
  `PROFILE_GRADE_SUBS`, series in `eta_bg` to O(2), keep terms with
  `eta_bg power ≤ 1 AND sigma_W power ≤ 1` (`:713`). This is the "η,σ_W ≤ 1" rule.
- `M.DERIVATIVE_MAP`, `M.grad_W`, `M.grad_mu`, `M.sigma_W`, `M.eta_bg`, `M.L_W`,
  `M.W0`, `M.mu_R`, `M.sorted_index`, `M.INCOMING_LEDGER`, and the DOF symbols
  `M.u`, `M.grad_u`, `M.theta`, `M.grad_theta`, `M.e_W`, `M.grad_e`.
- Second-jet profile values: `M.INCOMING_LEDGER["w1_profile_d{i}d{j}"]["value"]` for the
  thickness Hessian; mint `m1_profile_d{i}d{j}` KNOB symbols for the modulus Hessian, as
  `operator_dx`/`background_dx` do.

**DOF→row map** (from the committed comparator `S11c_b_cross_engine_comparator.py`
`extract_slab` ~L760): the strong `SLAB_OPERATOR` rows are
`U_BODY_BALANCE→U_MOMENTUM_ROWS` (DOF `u`, all 3 components — there is **no**
transverse/longitudinal split at the strong-row level), `THETA_BALANCE→MU_THETA`
(DOF `theta`), `E_W_BALANCE→THICKNESS_ROW` (DOF `e_W`). Build the per-row EL exactly as
`operator_from_density` does: `algebraic = diff(density, field)`,
`row = expand(algebraic − Σ_i dx(diff(density, first_jet_i), i))` — but with the
**Hessian-retaining `dx`** for the corrected operator and the **global (frozen) `dx`**
for the frozen operator.

---

## Enumeration + completed density (pin the recipe — do NOT rebuild the selected terms)

Enumerate **once** on the abstract spurion: `candidates = M.enumerate_new_candidates(M.bg)`
(`M.bg` is the abstract first-jet symbol vector `bg_1..bg_3`, `:1014`). Get the split with
`selected, omitted = M.quotient_independent_indices(exprs, M.basis_euler_signatures(exprs, M.basis_fields))`.
Indices are **zero-based**; the label of index `k` is `FIRST_JET_CONTRACTION_{k+1:02d}`. Only
**after** selecting an index do you live-substitute: `M.live_basis_substitution(M.grad_W)` for
source `W_BG`, `M.live_basis_substitution(M.grad_mu)` for `MU_R_BG`. (Passing `M.grad_W` into
`enumerate_new_candidates` is wrong — the candidates would contain no `bg`, and the μ
substitution would then replace nothing. Enumerate on `M.bg`, substitute per source.)

⛔ **Keep the engine's frozen density verbatim.** `density_frozen =
M.construct_energy("LAB_HELD").density` (10 uniform + 8 selected spurion/source). Do **NOT**
rebuild the selected 8 or the uniform 10 with fresh coefficients — that would break coefficient
identity with the emitted engine and poison the residual with relabeling (the #86 carrier-
doubling failure class). Build only the **completion**:
`density_omitted = Σ_S Σ_{k∈omitted(S)} c^{S}_k · (candidate_k live-substituted for S)`, each
`c^S_k` a **fresh independent** free symbol; `density_correct = density_frozen + density_omitted`.
Also split the frozen density by provenance into `density_uniform` (the 10 uniform selected
terms) and `density_selected_spurion` (the 8 spurion selected terms/source), via the term
labels `construct_energy` already carries — you need this split for the driver decomposition.

## Row EL definition (stored-energy amplitude — ε-stripped, inertia-free)

For each strong row use the **stored-energy EL amplitude** (no `epsilon`, no inertia term):
`EL(density, field) = expand( diff(density, field) − Σ_i DX(diff(density, first_jet_i), i) )`,
with `DX = frozen_dx` (global `M.DERIVATIVE_MAP` only) for the frozen operator and
`DX = hessian_dx` (that map plus the committed second-jet entries) for the corrected operator.
Rows and their `(field, first_jets)`: `U_MOMENTUM_1/2/3 = (u[a], grad_u[a])`,
`MU_THETA = (theta, grad_theta)`, `THICKNESS_EW = (e_W, grad_e)`.
Note the engine multiplies each emitted row by `epsilon` and adds an inertia term
(`epsilon*rhobr*u_tt`, `epsilon*mu_W*W_bg**2*e_tt`; `mu_theta_amplitude` is itself ε-free,
`:1465-1491`). The inertia carries **no** background spurion, so it cancels from every residual;
work with the ε-free, inertia-free amplitude and reconcile with the engine only in CONTROL_ENGINE
below.

## Compute and print (per row ∈ {U_MOMENTUM_1, U_MOMENTUM_2, U_MOMENTUM_3, MU_THETA, THICKNESS_EW})

`EL_FROZEN = first_shape_series( EL(density_frozen, field) with frozen_dx )`,
`EL_CORRECT = first_shape_series( EL(density_correct, field) with hessian_dx )`.
Emit (tag `S11CB_88_ROW_<NAME>`) a record with:
- `EL_FROZEN`, `EL_CORRECT`, `RESIDUAL = expand(EL_CORRECT − EL_FROZEN)`, and
  `RESIDUAL_TERMCOUNT`.
- **Three disjoint drivers** (partition the residual by linearity of the EL in the density):
  `UNIFORM_HESSIAN = first_shape_series( EL(density_uniform, field, hessian_dx) − EL(density_uniform, field, frozen_dx) )`;
  `SELECTED_SPURION_HESSIAN = first_shape_series( EL(density_selected_spurion, field, hessian_dx) − EL(density_selected_spurion, field, frozen_dx) )`;
  `OMITTED_CORRECT = first_shape_series( EL(density_omitted, field, hessian_dx) )` (the omitted
  invariants' full corrected EL — it naturally carries their own Hessian). Emit each + its
  termcount. **Assert** the reconstruction identity
  `expand( RESIDUAL − (UNIFORM_HESSIAN + SELECTED_SPURION_HESSIAN + OMITTED_CORRECT) ) == 0`
  (this is exact by linearity — a legitimate harness invariant produced by an independent route,
  not a tautology; it is one of the allowed asserts).
- `NEW_SECOND_JET_ATOMS` = sorted second-background-derivative profile atoms
  (`w1_profile_*d{i}d{j}`, `m1_profile_*d{i}d{j}`) occurring in `RESIDUAL` but in **no** term of
  `EL_FROZEN`. This is a **sufficient** (not complete) non-absorbability witness: a second-jet
  atom the frozen row cannot emit is not a re-parametrisation of any frozen coefficient.
- `GRADE_SUPPORT` = sorted `(a,b)`, `a,b∈{0,1}`, at which `RESIDUAL` has a nonzero
  `(eta_bg^a, sigma_W^b)` coefficient; and `PRE_TRUNCATION_TERMCOUNT` (termcount of the
  un-truncated residual) so a truncation-order error is visible.

## Non-absorbability test (constant-coefficient span — tag `S11CB_88_ABSORB_<NAME>`)

⚠ A rank over DOF monomials alone with background factors left in the entries is **invalid**
(`g·q` vs `H·q` gives rank-gain 0 spuriously). Do it as a **constant-coefficient membership**
test over the **joint** monomial basis:
- Monomials = the joint set of DOF jets **and** background-profile jets/bookkeepers
  (`w1_profile*`, `m1_profile*`, `eta_bg`, `sigma_W`). The scalar field = background-**independent**
  constants only (rationals and the fixed moduli `W0, L_W, mu_R, kappa_*, mu_*, k_W, …`).
  ⛔ Never leave `bg`/`grad_W`/Hessian atoms/profile values/free coefficients in a matrix entry
  as if they were the scalar field.
- `F` = span over that field of `{ ∂EL_FROZEN/∂p : p a free coefficient of density_frozen }`
  (one template column per free frozen coefficient, via coefficient differentiation), expressed
  on the joint monomial basis. `R` = `{ ∂RESIDUAL/∂c^S_k }` (one column per fresh omitted
  coefficient). Print `DIM_F`, `DIM_F_PLUS_R = dim(F + span R)`, and `RANK_GAIN = DIM_F_PLUS_R − DIM_F`.
  Print the operands; **assert nothing** about the gain. (`RANK_GAIN` is the complete
  absorbability measure; `NEW_SECOND_JET_ATOMS` is the cheap sufficient witness.)

## Scope (tag `S11CB_88_SCOPE`)
Emit a scope note object: this instrument measures the **LAB_HELD** stored-energy EL rows only
and is a **PY-side disturbance witness**. A nonzero row residual ⇒ the emitted PY row is not the
§1d-correct operator there. ⛔ A **zero** row residual does **NOT** clear that family — the
WL-frozen row (a different hand-coded 8-subspace, per #86) may change where PY's does not, and
`MATERIAL_ADVECTED` rows (built with a frozen `material_pullback` `dx`, `:1348-1377`) are not
measured here. Both are #89 work.

---

## Controls (the allowed asserts — each a real invariant via an independent route)

1. `S11CB_88_CONTROL_ENGINE` — **the EL recipe is the engine's.** For each row, independently
   extract the engine's emitted frozen row from `M.build_operator("LAB_HELD")` (or
   `M.operator_from_density(density_frozen, include_kinetic=False)`), divide out the common
   `epsilon`, drop the inertia term, and `first_shape_series`; print it beside the instrument's
   own `EL_FROZEN`; **assert** their difference is `0`. (This is the non-tautological replacement
   for a freeze-vs-freeze check: it validates the instrument's EL against committed engine code.)
2. `S11CB_88_CONTROL_RECON` — the driver reconstruction identity per row (defined above),
   asserted `== 0`.
3. `S11CB_88_CONTROL_SENTINEL` — chain-rule sentinels: for `grad_W[0]*u[0]` and `grad_mu[0]*u[0]`,
   print `frozen_dx` and `hessian_dx` before and after `first_shape_series`, showing the Hessian
   term appears only under `hessian_dx` and survives at `sigma_W^1`. Print objects; no assert.
4. `S11CB_88_CONTROL_JACOBIAN` — the coefficient Jacobian of `density_omitted`: assert its
   termcount `> 0`, and print that each `(source, index)` yields a distinct nonzero template
   carrying the correct source jet (`grad_W` for `W_BG`, `grad_mu` for `MU_R_BG`). Guards against a
   swapped/dead live substitution.

⛔ Do **not** assert anything about `RESIDUAL`, the drivers' values, `RANK_GAIN`,
`NEW_SECOND_JET_ATOMS`, or `GRADE_SUPPORT`. Those are measurements. Print them; conclude nothing.
There is no expected sign, count, or nonzero-ness to hit — a row whose residual is **zero** is as
informative as one whose residual is large. The only asserts are the four controls above (each a
structural invariant, not a physics target).

---

## Print/assert discipline (the three clauses — verbatim, binding)
1. The script may **PRINT** computed CAS objects. It may **NOT** state conclusions. No
   verdict strings of any kind — e.g. no `print("<row> is <adjective>")`, no labelling a
   residual with an interpretation word. Emit the CAS object; let the reader judge it.
2. **PRINT** each residual/driver/absorbability object; the only `assert`s are the four
   structural controls (`CONTROL_ENGINE`, `CONTROL_RECON`, and the termcount guard in
   `CONTROL_JACOBIAN`), each produced by an **independent route** (engine-committed rows;
   linearity of the EL; a termcount) — never a physics target. Compute → emit → then assert.
3. Interpretation belongs to the step/adjudication record, not this script.
- The only place physical symbols are combined by hand is in constructing the completion
  density (adjoining the omitted candidates with fresh coefficients) and the Hessian `dx`
  map (the committed `operator_dx` pattern). Every row expression must be **reached by
  computation** (`diff`/`dx`/`first_shape_series`), never typed.
- Tag **names** encode the row and the object (`RESIDUAL`, `OMITTED_CORRECT`, `RANK_GAIN`),
  never a value, sign, or the shape of the answer.

---

## Deliverable acceptance (what "done" means — no expected value)
- The file exists, is non-empty, runs to exit 0, and prints all `S11CB_88_ROW_*`,
  `S11CB_88_ABSORB_*`, `S11CB_88_SCOPE`, and the four `S11CB_88_CONTROL_*` tags as CAS
  objects.
- The four structural asserts pass (`CONTROL_ENGINE` == 0 per row, `CONTROL_RECON` == 0 per
  row, `CONTROL_JACOBIAN` termcount > 0). If `CONTROL_ENGINE` fails, the instrument's EL
  recipe does not match the engine — fix the recipe, not the control.
- Each emitted row/driver object is reached by computation (`diff`/`dx`/`first_shape_series`),
  not a typed literal — verify by inspection that no `RESIDUAL`/driver is a hand-written
  expression.
- The build's job **ends at compute-and-print**. Do not compare against any expected table,
  do not iterate toward a target, do not state which rows changed.

## Withheld (do not seek, do not infer)
The per-family cross-engine verdicts (admissibility / kinetic / advective / coupling) and
any expectation of how many rows change are **withheld on purpose**. They are the
orchestrator's to apply against your printed measurements. Building toward them would
defeat the instrument.
