# Independent review ROUND 2 — S11c-b carrier ablation-harness knife-list (orchestrator-written)

You are one of two independent review legs. Working dir: `/var/projects/toy_physics`. Repo read access; you may
write and run small check scripts under `/tmp`. ⛔ Do not modify the working tree. ⛔ You are not handed the other
leg's output — assess independently. **Prose assertions are discounted unless grounded in a cited engine line or a
small check script whose path + literal stdout you report.**

## Context
Round 1 of this knife-list came back **NOT-SOUND** with blocking findings. The list has been re-authored and the
packet **split** into a blind builder directive plus an orchestrator-side adjudication key. Your job: verify each
round-1 finding is **correctly folded**, verify the **new design choices** are right and implementable, and hunt
**new** defects. Report a finding only if it catches a way the harness would certify the wrong thing or fail to run.

## Artifacts under review (the pair)
- Blind builder directive (what `gpt-6-astra` receives — must be implementable AND leak-free of expected outcomes):
  `/var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_b_carrier_ablation_harness_directive.md`
- Orchestrator adjudication key (NOT given to the builder; holds the cones/expected bites — must be correct):
  `/var/projects/toy_physics/research/pde_ledger_v3/_measurements/S11c_b_carrier_ablation_harness_adjudication_key.md`

## Engines (verify every site/claim against these; cite lines)
- SymPy: `/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_b_brane_operator_sympy_audit.py`
- WL:    `/var/projects/toy_physics/research/pde_ledger_v3/mathematica/S11c_b_brane_operator_mathematica_audit.wl`
- N6 reference (atoms/reduction; ⛔ not imported by the harness): `.../scripts/S11c_c2_N6_diagnostic_sympy.py`
  (`pressure_coefficients` `:410-412`, `pressure_slots` `:321-322`).

## Verify each round-1 fold is CORRECT (not just present)
1. **WL entrypoint / OOM.** The directive now loads WL **definitions only** (up to the `Main variable-coefficient
   objects` marker `:1878`) and calls `evaluatedModel["EULERIAN","MATERIAL_ADVECTED","RHO4_CONSTANT"]["OPERATOR"]`
   directly, never the emit `Do` / `extractCouplingData` / `kernelOriginsFromOrigins` / `frozenEvaluatedModel`.
   Confirm: (a) the emit `Do` `:2190` really runs the kernel extraction before the `PRIMARIES_ONLY` test `:2206`
   (so the old path OOMs); (b) `evaluatedModel` `:1317` is fully defined before `:1878` and is callable after a
   definitions-only load; (c) `evaluatedModel[…]["OPERATOR"]` yields `U_MOMENTUM_ROWS`/`MASS_EVOLUTION_ROW`/
   `THICKNESS_ROW` carrying `pressureUpper`/`pressureLower`, and is the ~0.9 GB path, not the OOM one.
2. **Case pin MATERIAL_ADVECTED.** The key changes the branch from LAB_HELD to MATERIAL_ADVECTED because U is
   pressure-free under LAB_HELD. Verify independently: is the U carrier zero in LAB_HELD and **nonzero** in
   MATERIAL_ADVECTED, in BOTH engines? Is EULERIAN×MATERIAL_ADVECTED×RHO4_CONSTANT a real emitted primary case? Does
   changing the branch keep C_E (route EULERIAN) intact? Are θ and e_W still live under MATERIAL_ADVECTED?
3. **K_T whole channel.** SymPy FORM = remove all four `face_u`/`face_e` additions at `:2998-3040`; WL FORM =
   `virtualWork = 0` at `:1082-1083`. Confirm this drops BOTH the bare-`p` and the `lambdaXResponse affinity`
   subchannels (closing the round-1 Λ_X coverage hole), and that the two engines' FORMs are now equivalent
   whole-channel removals. Confirm K_T no longer names one function while patching another.
4. **K_A Λ_A-only.** SymPy FORM = null `Lambda_A_0` `:408` in the closure so Λ_A is removed but `Lambda_V_0` is
   left; WL FORM = drop `lambdaAResponse affinity` `:1080` leaving `lambdaVResponse normalVelocity`. Confirm this
   moves the θ carrier and does NOT touch the Λ_V channel or U/e_W. Is nulling `Lambda_A_0` a genuine structural
   channel removal (FORM), or does it read as a coefficient rescale?
5. **K_W collapse (not exchange).** SymPy FORM = inject `delta_p_minus→delta_p_plus` + `d_w_delta_p_minus→
   d_w_delta_p_plus` into the substitution applied at `filtered_substrate` `:2045-2048` (reaches `closure_shape_deriv`
   / `virtual_work_shape_deriv`); WL FORM = `pressureField[-1] := pressureUpper[…]` `:1015`. Confirm a full +/−
   **exchange** is observationally dead (symmetric ± carriers) while the **collapse** is rank-changing and bites.
   Confirm the SymPy collapse actually reaches the delta_p-bearing substrates at `:2045-2048`.
6. **SymPy access / citations.** `build_operator` returns `casify(...)` (nested `sp.Tuple`, `:3261`/`:603-613`), so
   the observation uses `named_tuple_row(operator,"…")` `:2584`, not `operator["…"]`. The θ observation targets the
   post-fold rows `:2998`/`:3028`/`:3042`, NOT the energy templates `:2367`/`:2372`/`:2377`. Confirm both.
7. **Drift guard (SymPy).** No single-case `SLAB_OPERATOR` tag is emitted on import; the guard compares the retained
   single-case object from imported `build_operator` against an unablated temp copy. Confirm this is coherent.
8. **Self-tests pinned** (extractor-order P→0-before-∂; dead-path bulk site; live-rescale contrast). Are they
   concrete and sufficient to catch a self-reporting harness? Is the extractor-order test's logic right (∂ then →0
   vs →0 then ∂)?

## Verify the NEW design + the split
- **Directive blindness.** Does the builder directive leak any expected **value / count / ratio** or **pass-condition**
  (must-MOVE / must-NOT / zero-vs-nonzero interpretation / "should move")? It must not — those belong only in the
  key. Confirm the directive tells the harness to print the COMPLETE carrier (all rows × all slots) with no selection.
- **Adjudication key correctness.** Are the per-knife cones in the key right at the CARRIER level (K_A→θ move,
  U/e_W not; K_T→U+e_W move, θ not; K_W→collapsed-column rank change)? Is dropping the "Λ_V channel" and "bulk base"
  from the must-NOT set correct (they are pressure-free ⇒ carrier trivially zero ⇒ no information)? Is the
  disjointness claim (K_A vs K_T patch disjoint sites; neither the shared `affinity` `:1079`) still right?
- **Discipline.** Exactly ONE function + ONE FORM per knife? Three PRINT-not-assert clauses intact? Any site whose
  named function/line does not carry what is claimed? Any NEW defect the re-author introduced?

## Output
A per-item verdict with cited lines, then a final line: **SOUND** (ready to build) or **NOT-SOUND** (with the exact
site + fix for every blocking issue).
