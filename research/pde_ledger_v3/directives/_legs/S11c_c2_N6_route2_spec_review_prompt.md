# Independent review — astra's S11c-c2 N6 route-2 construction spec

## Artifact (read LAST)
`/var/projects/toy_physics/research/pde_ledger_v3/_measurements/S11c_c2_N6_route2_spec_astra.md`

A **Codex(astra)-authored** construction spec for the c2 N6 representation-invariance control's route 2, to be folded
into §5c + a build directive. Working dir `/var/projects/toy_physics`; repo-root-relative paths unless absolute.
Physics-bearing document review — the companion script does NOT exist yet (defer executable script-control tests to
the build, but flag any construction that cannot be built as specified).

## Method — derive your OWN account FIRST, from the code, THEN read the spec
Form your own answer to "what is the correct route-2 construction for the c2 face-slot self-energy increment" from
the code below, then quote both sides for every finding. ⛔ You are NOT handed, and there is NOT, any expected value
for the N6 residual.

Sources of truth (read and reason from these):
- The parent face N6: `scripts/S11c_a_interface_geometry_sympy_audit.py` `task_rep_invariance` (1481-1496),
  `build_geometry_quantity` / `build_virtual_constraint_raw` / `build_evolution_raw` and the `route` plumbing,
  `build_material_face_source` (703-813, covector inverse-transpose 742-755, material V_s 794).
- The c2 increment + extract: `scripts/S11c_c2_selfenergy_fold_sympy_audit.py` `extract` (325-342), `build_face`
  (476-515), the current audit `REP_INVARIANCE_*` loop (1087-1107), `FLIP_FACE_SLOPE` (384-387).
- The face fold + μ_θ: `scripts/S11c_b_brane_operator_sympy_audit.py` `face_generalized_force_rows` (2135-2179),
  `bind_mu_theta_operand`, `build_operator` (2762-2851), `operator_from_density` (2290-2384), `material_pullback`
  (1916-1981), `density_pair` (1878-1892).

## The load-bearing questions this review MUST settle (reason from the code + parent)
1. **The μ_θ binding decision (§1 of the spec): route 2 binds MATERIAL μ_θ, route 1 Eulerian μ_θ.** Is this the
   correct N6 for the self-energy INCREMENT? The spec argues the parent tests only *geometry* representation (leaving
   μ_θ a placeholder in both routes), and that the "native material constitutive-and-face" control additionally binds
   the material constitutive μ. Is that the right control — or should both routes bind the same (Eulerian) μ_θ (a
   geometry-only N6)? What does "representation invariance of the increment" actually require? A wrong choice makes
   the advection probe either vacuous or testing the wrong thing. Substantiate from the parent + what μ_θ means.
2. **The reverse-u grade-suppression finding (§4, §6).** The spec claims: the OPEN pressure coefficients carry no θ/μ
   (N4 enters only via the CLOSED pressure source μ); and the N4-induced reverse-u carrier contribution is
   grade-suppressed OUT of the retained (ε,η,σ_W) rectangle (LAB_HELD U carrier = 0; MATERIAL_ADVECTED = O(σ_W); their
   N4-only product ≥ σ_W²), so a *mandatory* "reverse-u channel must survive" requirement is INCOMPATIBLE with this
   fold and must be replaced by provenance + blockwise sensitivity reporting. **Is this correct**, or is it an excuse
   that lets a vacuous route 2 through? Check the grade bookkeeping and whether the N4 signal genuinely lives in the
   forward μ/pressure-source couplings (not a reverse U carrier). This is the crux — a wrong call here either forces a
   nonexistent channel or hides a real one.
3. **Buildability / no double-map / no annihilation:** native material face builders folded into the SAME δp symbols;
   `material_pullback` used only on the bulk *scalar* to get μ (⛔ never on rows/slots); the Jacobian only in that
   quadratic scalar; the c1 response opaque; route 1 = imported expanded rows. Does the construction type-check
   against the actual function signatures? Are the code-support findings (§8: no material c2 carrier factory; μ-source
   binding must be added; normal-only tilt injection missing) accurate?

## Also verify
- **Controls:** advection = a tag t on the material bulk-density θ-substitution, recompute μ, bind at BOTH face and
  c1-source, t=1 baseline / t=0 corrupted; live on RHOBR_CONSTANT, computed absence on RHO4_CONSTANT (⛔ not A−A).
  Tilt = an injectable Eulerian carrier factory PIT-checked vs the imported carrier, normal-only mutation (⛔ not
  `FLIP_FACE_SLOPE`). Are these one-sided, at-source, and non-vacuous?
- **PIT contract (§7):** several primes, ≥8 draws/prime/cell, joint singular rejection, branch-cell selection before
  modular reduction, on-shell q/k not sampled independently, degree bounds, FN bound with the bad-prime condition
  handled (⛔ not "multiple primes ⇒ done"). Sound?
- **Emit/dimensions (§8):** operand-before-residual; residual a finding, never an exit/assert; `[L,T,M]` + able-to-fail
  + (ε,η,σ_W) on every object. Leak check: does anything encode the sign/shape of R_N6? Is the target genuinely
  withheld?

## Physics filter
Report a finding only if it catches a way the built N6 control could be **wrong, vacuous, intractable, or
answer-leaking**.

## Output
Findings each with file:line quote, the source it violates, why it matters, minimal fix. If nothing outstanding
changes what will be computed or claimed, say **SPEC CLEAR**. Evidence-first, brief.
