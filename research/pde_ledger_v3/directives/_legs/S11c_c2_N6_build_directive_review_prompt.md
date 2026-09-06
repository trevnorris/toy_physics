# Independent review — S11c-c2 N6 build directive

## Artifact (read it LAST — see method)
`/var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_c2_N6_sympy_build_directive.md`

An **orchestrator-written build directive** governing a new SymPy companion diagnostic that will compute the
**corrected** N6 representation-invariance control. Working dir: `/var/projects/toy_physics`. Paths below are
repo-root-relative unless absolute.

## What you are checking
Whether this directive **faithfully and tractably** commissions the corrected N6 control, **without leaking an
answer or a recipe**, so it can go to the builder (astra) with no further change to what will be computed or claimed.
This is a physics-bearing directive; review it by **reading** (the script does not exist yet — defer executable
script-control tests to the build, but flag any control that cannot be built as specified).

## Sources of truth (read these FIRST, form your own view, THEN open the directive)
1. `directives/S11c_c2_SHARED_PHYSICS.md` §5c (the governing object; also skim §3a close, §3c weak increment +
   slot-linearity, §1d routing, §5d density). §5c was corrected + cleared this session (commit `30d4b72d`).
2. `_measurements/S11c_c2_N6_spec_adjudication_record.md` — WHY the old object was wrong: the two anchorings
   (`LAB_HELD` / `MATERIAL_ADVECTED`) are **distinct physical setups**, not Eulerian/material representations; the
   real N6 compares the two **coordinate constructions of the increment at a FIXED anchoring**, mapped back by
   `Δρ = δρ_E + u·∇ρ⁰`. It summarizes the parent pattern (S11c-a §2c/§5a, parent N4/N6, sibling c1 §5a).
3. `_measurements/S11c_c2_N6_build_design.md` — the vetted tractable shape (carrier-first + finite-field PIT).
4. The two code anchors the directive cites — verify they say what the directive claims:
   - the F/G companion carrier scaffold: `scripts/S11c_c2_FG_diagnostic_sympy.py:99-134`;
   - the audit's OBSOLETE cross-anchoring loop: `scripts/S11c_c2_selfenergy_fold_sympy_audit.py:1087-1107`
     (and the DENSITY loop at `:1085-1086` which must STAY);
   - the material builders it points astra at: `scripts/S11c_b_brane_operator_sympy_audit.py` (~1916-1981,
     ~2723-2777) and `scripts/S11c_a_interface_geometry_sympy_audit.py:703-815`.

Form your OWN account of what the corrected N6 object is before reading the directive; then quote both sides for
every finding. ⛔ You are NOT handed, and there is NOT, any expected value for the residual — §5c fixes no target.

## Substantiate these independently (a prose "looks right" is worth nothing — reason from §5c / the parents)
1. **Right object / right question.** Does the directive commission the corrected N6 = **Eulerian vs
   material-coordinate construction of the SAME self-energy increment at a FIXED anchoring `α` and density `ρ`**,
   `T_{M→E}` applied once to the increment — and NOT the cross-anchoring proxy? Is the anchoring/density genuinely
   held fixed across the two routes? Is `Δρ` used only as the intra-anchoring representation map (⛔ never to bridge
   `LAB_HELD ↔ MATERIAL_ADVECTED`)?
2. **Native second route.** Route 2 must be built **natively** from the S11c-b material builders + S11c-a
   face-flattening, ⛔ not a re-pullback of the Eulerian, and closed with the SAME imported same-`α` c1 response. Is
   `SLAB_M` the material operator BEFORE the map-back (⛔ not double-mapped)? Is the directive right that `build_case`
   cannot build the material route? Is the demand for a FULL map-back (trial fields AND test covectors AND Jacobian,
   ⛔ not the trial-field-only `representation_pullback`) correct and load-bearing?
3. **Independence / no tautology.** Are the two routes structurally independent, so `R_N6` is not zero by
   construction? Are the two one-sided controls (tilt on the Eulerian route; N4-advection omission in the map-back
   only) each mutating ONE route at its source, capable of moving that route's operand while leaving the other? Is
   the `A−A` prohibition (emit structural absence, e.g. `RHOBR_CONSTANT`, ⛔ not a self-difference) correct?
4. **Tractability.** Is carrier-first (`S_P = S − S|_{P=0}`, base cancels by §3c linearity) plus **exact finite-field
   PIT** (over several primes, shared generic samples across routes, denominators cleared, singular/branch cells
   handled, false-negative bound recorded) actually sufficient to avoid the F/G full-symbolic wall — or is there a
   step (e.g. building `SLAB_M`, or `T_{M→E}`) that still forces whole-object symbolic materialization? Is
   "evaluate/truncate coefficients EARLY, before expansion" load-bearing and correctly placed?
5. **The audit correction.** Is retiring the cross-anchoring `REP_INVARIANCE_*` loop (`:1087-1107`) and re-emitting
   it as `S11CC2_ANCHORING_L_MINUS_M` (§5c:362-370, no zero target) correct — and does the directive correctly leave
   the `DENSITY_LIVE_MINUS_FROZEN` loop (`:1085-1086`) intact?
6. **Leakage / iterate-to-target.** Does anything let the builder iterate toward a residual value it can see? Is the
   residual correctly framed as a computed finding with NO supplied target, NEVER a builder exit / assert / pass
   condition? Are all anti-examples placeholder-only (no real step content)? Is the object NAMED (not the recipe /
   not the answer)? Does any tag name, or any supplied text, encode the sign/shape of the answer?
7. **Missing premise / under-specification.** Anything the builder needs (a field, a map, a closure ingredient, a
   coordinate) that is prose where it should be an equation or a named source — the kind of gap that would make the
   built control silently vacuous or wrong.

## Physics filter
Report a finding only if it catches a way the **built N6 control could be wrong, vacuous, intractable, or
answer-leaking** — not a stylistic preference. If you believe §5c itself (not just the directive) is wrong, say so
explicitly and show why from the parents; that routes back to the spec.

## Output
List findings, each with: the directive quote (with line), the source-of-truth quote it violates (with file:line),
why it matters, and the minimal fix. If nothing outstanding changes what will be computed or claimed, say **CLEAR
TO BUILD**. Keep it evidence-first and brief.
