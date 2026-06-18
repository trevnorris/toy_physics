# Directive — Path-A build Chunk 1c: self-consistent-balance validation of the closed solve

**Program:** Path-A self-consistent closed solver, per `decisions/10_pathA_solver_build_design.md` (the 1c spec)
and the chunk-1b result (`reports/patha_chunk1b_closed_newton.md`). **You (Codex) CODE this; Claude reviews +
runs a transliteration-fidelity audit + adversarial review.** Sandbox: workspace-write; every script you RUN gets
`timeout 600` and you iterate it to exit 0. CPU-only.

**Target-blind ENGINEERING validation:** placeholder `S_Σ` (the `smooth_positive_placeholder_v1` family, the same
default 1b solves with), parameters AWAY from the `D0→0` near-pole. **NO calibration, NO `R_norm`/target
comparison, NO physical export, NO GATE-A forms, NO field→coefficient extraction map / `D0` construction.** NEVER
touch `research/pde_audit/simulation/` or the `physical_export_permitted` guard.

This chunk does NOT change the operator, the residual, the gates' thresholds, or the constitutive providers. It is
the VALIDATION layer over the converged closed solve: does the self-consistently-solved `R0(w)` actually satisfy
the derived static balance term-by-term, does the closed solve converge under grid refinement, and does promoting
`R0` leave the already-established matter/charge/Gauss conservation intact.

## Read first
- `decisions/10_pathA_solver_build_design.md` — the **1c** bullet (recompute wall LHS/source from the converged
  `{ψ,A,R0}`; reuse step-4 convergence + step-6 conservation style; R0-aware raw-field observables; no extraction
  map) and the stop conditions.
- `derivations/pathA_01_return_source_and_balance.md` — D3 (the static self-consistent balance whose terms 1c
  decomposes) and D4.
- `reports/patha_chunk1b_closed_newton.md` — the converged 1b result you validate (residual 1.07→1.7e-12, R0∈
  [0.988,0.998] emergent, JVP rel 4e-11) and the **four carried MINORs** (see Part B).
- The closed solve + terms you reuse: `src/stage1_solver/patha_closed_newton.py` (`run_patha_closed_newton`,
  `default_closed_branch_config`, `default_closed_s_sigma_spec`, `target_token_scan`, `return_source_fidelity_diagnostic`),
  `src/stage1_solver/coupled_branch.py` (`patha_closed_wall_terms` returns the full `StaticBalanceTerms`:
  `flux_divergence, gradient_square, potential_gradient, source, lhs, residual, face_fluxes, center_gradient`;
  `patha_radial_reduced_return_source`, `patha_return_source_density`,
  `confinement_wall_to_matter_coefficient_torch`, `unpack_closed_coupled_fields`, `initial_closed_branch_state`,
  `resample_closed_branch_state`, `patha_closed_branch_residual`), `src/stage1_solver/patha_static_balance.py`
  (`static_balance_terms`, `StaticBalanceTerms`).
- The validation MATH + style you REUSE (do NOT re-derive): `src/stage1_solver/convergence.py`
  (`observed_order_from_three`, `richardson_estimate`, the conservative volume-restriction + raw-field
  self-difference pattern, the null/floor verdict vocabulary) and `src/stage1_solver/conservation_diagnostics.py`
  (the step-6 stationary FV conservation diagnostics + the `_not_a_physics_gate` discipline).

## Part A — what to build (the 1c validation harness)

Create a new module `src/stage1_solver/patha_closed_validation.py` (a `run_patha_closed_validation(...)` +
`write_patha_closed_validation_report(...)` + a `main()`), a report
`reports/patha_chunk1c_self_consistent_validation.md`, and tests `tests/test_patha_closed_validation.py`. REUSE the
convergence/conservation math by import; do NOT duplicate the order/Richardson formulas or the FV divergence
operators.

### A1. Self-consistent-balance term decomposition (the core; must NOT be a tautology)
From a converged closed state `{ψ,A,R0,μ}` (solve once at the 1b default placeholder via the existing closed
machinery), call `patha_closed_wall_terms` and report, per the wall grid:
- the L∞ (and L2) magnitude of EACH physical term separately: `flux_divergence`, `gradient_square`,
  `potential_gradient`, and `source`;
- the dominant-term magnitude `D = max L∞ over {flux_divergence, gradient_square, potential_gradient, source}`;
- the balance residual `lhs − source` reported **relative to D** (`||residual||∞ / D`), NOT in isolation.

**Why this is not a tautology (state this explicitly in the report):** the Newton already drove the full residual
(which CONTAINS `wall.residual`) to ~1e-12, so re-printing the residual alone proves nothing (0≈0). The genuine
content is (i) that the individual terms are O(1)-comparable (a real cancellation, not all-near-zero), and (ii) an
**independent recomputation** of the LHS and source through a second code path that does NOT echo the residual:
  - **source:** recompute `S_j = Σ_i ΔV_i^r·(−k1·ρ0)` from the shared `confinement_wall_to_matter_coefficient_torch`
    + `grid.radial_shell_volumes` (as `return_source_fidelity_diagnostic` does) and compare to
    `patha_closed_wall_terms(...).source`;
  - **flux divergence:** recompute the conservative wall flux-divergence `−(F[j+1]−F[j])/Δw` with
    `F[j+½]=T_w(face)·ΔR/Δw` from `R0` + the provider `T_w` through a SEPARATE local reconstruction (e.g. an
    independent NumPy/torch stencil, including the one-sided mouth stencil and the zero-traction exit `F[-1]=0`),
    and compare to `StaticBalanceTerms.flux_divergence`;
  - report `max_abs_diff` for both independent recomputations (these guard against a faithful-but-wrong
    operator — the transliteration sin; MMS/Newton cannot catch that).
Genuine gates here: `balance_terms_nontrivial` (each of the three LHS terms AND the source has L∞ above a stated
floor, so the balance is non-degenerate), `independent_source_recompute_matches`, `independent_flux_recompute_matches`,
`balance_residual_relative_below_solver_floor` (`||residual||∞/D` ≤ the achieved Newton residual scale). Do NOT add
a gate that merely restates `||wall.residual||∞≈1e-12`; if you report that number at all, label it a diagnostic.

### A2. Grid-convergence study of the CLOSED solve (resolves 1b MINOR-1)
1b's convergence evidence rested on a single grid (the K=0.18 continuation step was a near no-op). 1c provides the
genuine multi-point evidence by **solving the closed self-consistent system at ≥3 grid levels on a fixed
refinement ratio** and measuring self-convergence:
- Reuse the closed-solve machinery (`solve_newton_jvp` + `patha_closed_branch_residual` + the closed colored-sparse
  preconditioner + `initial_closed_branch_state` + `resample_closed_branch_state` to warm-start each finer level).
  This is the CLOSED (R0-promoted) solve, NOT `run_branch_continuation` (which is the frozen-geometry path).
- **R0-aware restriction + raw-field self-difference:** extend the `convergence.py` conservative volume-restriction
  to the closed state — restrict `{ψ_R,ψ_I,a0,ar,aw}` as today AND restrict the 1D wall field `R0(w)` along `w`
  (volume-consistent with the 2D restriction's `w`-measure). Report the raw-field self-difference + its observed
  order (reuse `observed_order_from_three`, `richardson_estimate`).
- **R0-aware surrogate observables:** add `R0 L2 norm`, `R0 min/max/range`, and the **static-balance
  residual-relative-to-dominant-term** (from A1) as observables tracked across levels; the latter should stay at or
  below the solver floor (or decrease) under refinement — that, not a single-grid number, is the convergence claim.
- **Cost / the 600s rule (hard):** each runner script must finish under `timeout 600`. The 1b (6,6) closed solve
  was ~123 s, dominated by the colored-Jacobian rebuild `every_newton_step`; you may use cheaper, solution-
  PRESERVING levers to fit ≥3 levels — switch the preconditioner rebuild policy to `once_per_continuation` (it is
  only the GMRES `M`; the true JVP is preserved — the step-3b lesson) and/or drop the continuation to a single
  representative `K`. Choose a ladder (e.g. start coarse, fixed 2× ratio) that completes ≥3 levels in budget; if a
  finer level genuinely cannot fit, STOP at the laptop limit and REPORT it honestly (reuse the step-4
  laptop-limit/`max_level_seconds` reporting). If you must split per-level into separate scripts that each persist a
  machine-readable table plus an aggregator, that is acceptable. **Reaching ≥3 completed levels is required for the
  convergence verdict to PASS**; fewer than 3 is a finding, not something to silence.
Genuine gates: `minimum_three_levels`, `all_completed_levels_converged`, `fixed_ratio_ladder`,
`raw_field_self_difference_decreases_or_at_floor` (incl. the R0 channel),
`balance_residual_relative_stable_under_refinement`.

### A3. R0-aware conservation diagnostics
Confirm that promoting `R0` to a solved field leaves the matter/charge/localized-Gauss conservation already
established in step 6 INTACT on the closed background (the confinement geometry is now `confinement_radius=R0`
rather than the frozen profile). Reuse `conservation_diagnostics.py` machinery on the converged closed state's
matter/gauge sub-block (the closed→coupled field view) — the independent center-gradient Gauss reconstruction +
the FV budget closure, keeping the same `_not_a_physics_gate` labeling for the roundoff FV identity. Report whether
the closed (R0-solved) conservation residuals match the frozen-geometry baseline to the solver floor. Keep the
honest step-6 framing (isotropic stationary branch ⇒ spatial-current sectors are null/floor diagnostics; the
current-carrying test stays deferred). Do NOT invent a new "throat conservation law" — R0 is determined by the
static force balance (A1), not by a transport conservation law; say so.

## Part B — apply the three carried 1b report MINORs (honest-labeling only)
These are presentation/labeling fixes to the chunk-1b harness/report (`patha_closed_newton.py` report generator);
NO math, gate-threshold, or operator changes. Re-run `python -m stage1_solver.patha_closed_newton` to regenerate
`reports/patha_chunk1b_closed_newton.md` so file and generator agree, and confirm the counted gates / numbers are
unchanged.
- **MINOR-2 (`constitutive_positive_margin`).** Annotate it honestly as a constitutive-positivity *sanity guard*
  for the placeholder family (near true-by-construction for this family; retained as a stability precondition for
  the background solve), NOT independent physics evidence. Keep it counted; add the clarifying note in the report.
- **MINOR-3 (placeholder-derivative `max_rel == max_abs`).** Add a note that relative == absolute here because the
  relative denominator is floored at `max(1, ||analytic||∞)` and the analytic magnitudes are < 1 — i.e. it is the
  denominator-floor convention, not a coincidence.
- **MINOR-4 (relaxed closed-Newton atol).** Note in the report that the closed-Newton gate's `residual_atol`/`rtol`
  is `2e-9` while the achieved residual is `~1.7e-12` (far below), so the relaxed atol is not masking
  non-convergence.

## Acceptance criteria (target-blind; genuine gates only; restatements split as `_not_a_physics_gate`)
1. **Non-tautological balance validation (A1):** the term-decomposition + the TWO independent recomputations
   (source, flux-divergence) pass to a stated tol, the terms are demonstrably non-degenerate, and the relative
   residual is reported against the dominant term — NOT a re-print of `wall.residual`.
2. **≥3-level closed convergence (A2):** at least three completed, converged levels on a fixed ratio; raw-field
   self-difference (incl. R0) order/decrease reported; balance-residual-relative stable/decreasing under
   refinement; laptop limit reported honestly if a finer level is dropped.
3. **R0-aware conservation (A3):** closed (R0-solved) conservation residuals reported and shown consistent with the
   frozen-geometry baseline to the solver floor; FV roundoff identity kept as `_not_a_physics_gate`.
4. **1b MINORs applied (Part B):** the three labeling notes present in the regenerated 1b report; 1b gates/numbers
   unchanged.
5. **Additive / no regression:** full `pytest software/stage1_solver/tests -q` passes; effective-closure and 1a/1b
   paths unchanged.
6. **Target-token scan clean** over all new 1c files (extend the scan path list); no calibration; firewall + export
   guard untouched.

## Stop conditions (HALT and report, do NOT paper over)
- The independent source/flux recomputation does NOT match the operator term (a faithful-but-wrong operator — a
  finding); or the balance terms are degenerate (all near zero ⇒ the ~1e-12 residual was 0≈0, not a real balance);
  or fewer than 3 levels complete in budget; or `R0` goes nonpositive at any level; or conservation residuals on
  the closed background are materially worse than the frozen-geometry baseline (R0 promotion broke a conservation
  property). Any of these is a result to surface, not to silence.

## Recurring sins to avoid
No tautological/can't-fail gate (esp. the "residual≈1e-12" restatement — that is the headline trap of this chunk);
no circular self-check (the independent recomputations must build from the kernel/stencil, not echo
`StaticBalanceTerms.lhs/residual`); no hardcoded values; no silent fallback that degrades to pass; honest labels
(engineering validation at placeholder `S_Σ`, NOT physics; the FV budget identity is roundoff, not a conservation
law). Stay within `timeout 600` per script — a timeout is a failure to reformulate (coarser ladder / cheaper
preconditioner rebuild / split scripts), never a cap to raise.

## End your final message (this is your report to the reviewer)
with: (1) the balance term-decomposition numbers (per-term L∞, dominant term, relative residual) + both
independent-recompute `max_abs_diff` values and WHY this is not a tautology; (2) the convergence ladder actually
run (levels, ratio, per-level converged?/wall-clock, R0 + field self-difference orders, balance-residual-relative
across levels, laptop limit if any); (3) the R0-aware conservation result vs the frozen-geometry baseline; (4)
confirmation the three 1b MINOR notes are in the regenerated 1b report with gates/numbers unchanged; (5) the
genuine-gate list with PASS/FAIL; (6) any stop condition that tripped or anything in the 1c spec that was
under-specified for coding.
