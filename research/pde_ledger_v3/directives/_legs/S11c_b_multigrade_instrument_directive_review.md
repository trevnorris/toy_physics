# Decision-list review — S11c-b background-order MULTIGRADE INSTRUMENT directive

You are reviewing an **orchestrator-written build directive** (a decision list) BEFORE any builder runs. Your
job is to find defects in the directive that would cost a build round or produce a wrong/uninformative
instrument — ambiguities, missing constraints, a leaked expected value, an object named as a recipe, an
acceptance criterion that would pass with a defect still in place, or a "level-above" error (the directive
asks for the wrong object, or an object that cannot answer the question it exists to answer). Report a
findings list; if the directive is sound, say so and say what you checked. **Do not edit anything.**

## The directive under review
`/var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_b_multigrade_instrument_build_directive.md`

## What the instrument is for (the question it must serve, downstream)
The committed adjudication layer found twenty genuine cross-engine differences between the SymPy engine and
the blind Wolfram engine for the S11c-b variable-coefficient brane operator / coupling kernel. Downstream (in
a step record, by hand, against the spec) we must decide, per case, which background-bookkeeper orders each
difference lives at, so we can judge which engine's retention is correct. This instrument's ONLY job is to
**measure** the `(eta_bg, sigma_W)` multigrade of each of the twenty aligned operand pairs and their
residual. It must state no conclusion about which engine is right.

## What you are handed (read these; derive your own view before judging the directive)
- The adjudication layer the instrument must reuse (do NOT let the instrument re-derive physics):
  `/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_b_adjudicated_comparison.py`
  — in particular its `transform`, `_bridge_d`, `_family_cases`, `_arithmetic_residual`, `_kinetic_residual`,
  `PROFILE_GRADE_SUBS`, and how its `main` builds `transformed_values` (the aligned `left`/`right`).
- The comparator it imports: `.../scripts/S11c_b_cross_engine_comparator.py` (`load_py`, `load_wl`,
  `DEFAULT_PY`, `DEFAULT_WL`, `serialise`, `WL_TO_PY_RENAME`).
- The SymPy engine module: `.../scripts/S11c_b_brane_operator_sympy_audit.py`
  (`PROFILE_GRADE_SUBS`, the `eta_bg`/`sigma_W` symbols, the profile substitution).
- The spec (for background: what `eta_bg`, `sigma_W`, and the background-order grading mean):
  `.../directives/S11c_b_SHARED_PHYSICS.md` §2a, §3a, §3d.
- The complete, untruncated stdout of the fresh adjudication run that produced the twenty cases (route lines
  `FLAG`/`RESIDUAL_BULK`, and each case's `operand_A`/`operand_B`/`A_minus_B` in SymPy `srepr`):
  `/home/trevnorris/.s11_build/S11c_b_adjud_fullres_run.out`
  (if you cannot read that path, regenerate with `python3 scripts/S11c_b_adjudicated_comparison.py` from
  `research/pde_ledger_v3/`, ~10 min).

## Required method — verify against the actual operands, do not reason in prose
1. From the run output, confirm the directive's enumerated twenty cases (twelve `FLAG` core + eight
   `RESIDUAL_BULK`) are exactly the genuine set, and that the `DIVERGENCE_ROUTE_FIXTURE` `FLAG` line is
   correctly excluded. Flag any case the directive names that is not in the run, or any genuine case it omits.
2. Read the actual `operand_A`/`operand_B` `srepr` for a spread of the twenty (at least one ADVECTIVE, one
   KINETIC, one ADMISSIBILITY THETA, one COUPLING_KERNEL). Determine whether the directive's computation is
   **well-defined and faithful on the real operands**:
   - Is the `(eta_bg, sigma_W)` multigrade well-defined for every operand shape present (scalar Expr; the
     KINETIC tuple; any `Association` container)? Does the directive's recursion instruction cover them all?
   - Do any operands contain `eta_bg` or `sigma_W` in a denominator (so exact polynomial coefficient
     extraction fails and the Taylor-series-with-remainder path is required)? Is the directive's remainder
     handling adequate, or could it silently drop content? Does `N ≥ 4` suffice, or could a real operand
     carry a bookkeeper power the instrument would miss without flagging?
   - Will the prescribed reuse (`transform` with `bridge_a=True`, then `_bridge_d`, `collapse=None`,
     `active_names = dict(WL_TO_PY_RENAME)`) reproduce the SAME aligned `left`/`right` the layer routed for
     these cases? If the layer's main does something the directive omits (e.g. a rename/collapse/energy-gate
     branch) that changes the operand, that is a level-above defect — name it.
   - Is the residual the directive computes (via `_arithmetic_residual` / `_kinetic_residual`) the SAME object
     the layer routed as `A_minus_B`? If not, the multigrade would be of the wrong residual.
3. **Leak check (rule 5).** The instrument must not encode any expected grade population, any "first/higher
   order" classification, any coefficient value, any "which engine is missing/right", or any physics
   acceptance. Confirm the directive neither states these nor forces the builder toward them. Equally, confirm
   the guards it does specify (`RECONSTRUCTION`, `GRADE_DIFFERENCE`) are pure decomposition self-consistency
   (a correct extraction round-trips to the zero object regardless of the physics) and are **not** a physics
   acceptance a builder could tune toward. If a guard could be satisfied by a wrong instrument, say so.
4. **Object-vs-recipe (rule 3).** Is the object ("the `(eta_bg, sigma_W)` multigrade of each operand and
   residual") named cleanly, or does the directive over-specify a derivation path that could manufacture a
   particular answer? Conversely, is anything under-specified such that two builders would produce
   incomparable outputs?
5. **Acceptance sufficiency (rule 7).** Would the stated acceptance pass an instrument that computed the
   wrong thing (e.g. graded in the wrong variables, merged container components, or dropped a non-polynomial
   remainder)? If yes, name the missing acceptance clause.

## Physics filter
Report a finding only if it catches a way the instrument would give a wrong or uninformative measurement, or
would cost a build round. Do not report style. Do not propose the physics answer (which engine is right) —
that is not this instrument's job and not yours here.

## Output
A numbered findings list (each: what is wrong in the directive, why it costs a round or corrupts the
measurement, and the minimal fix), then a one-line verdict. Do not edit; do not build.
