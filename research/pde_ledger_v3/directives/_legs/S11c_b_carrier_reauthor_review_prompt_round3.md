# Independent review ROUND 3 — S11c-b carrier ablation-harness knife-list (orchestrator-written)

You are one of two independent review legs. Working dir: `/var/projects/toy_physics`. Repo read access; you may
write/run small check scripts under `/tmp`. ⛔ Do not modify the working tree. ⛔ You are not handed the other leg's
output — assess independently. **Prose assertions are discounted unless grounded in a cited engine line or a small
check script whose path + literal stdout you report.**

## Context
Round 2 came back NOT-SOUND on a small set of precision items (the whole design PASSED both legs). Those were
folded. Verify each fold is correct, verify the directive is now fully blind, and hunt any NEW defect. Report a
finding only if it catches a way the harness would certify the wrong thing or fail to run.

## Artifacts (the pair)
- Blind builder directive (what `gpt-6-astra` receives — must be implementable AND leak-free of expected outcomes):
  `/var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_b_carrier_ablation_harness_directive.md`
- Orchestrator adjudication key (NOT given to the builder; holds cones/expected bites — must be correct):
  `/var/projects/toy_physics/research/pde_ledger_v3/_measurements/S11c_b_carrier_ablation_harness_adjudication_key.md`

## Engines (verify every claim against these; cite lines)
- SymPy: `.../scripts/S11c_b_brane_operator_sympy_audit.py`   · WL: `.../mathematica/S11c_b_brane_operator_mathematica_audit.wl`
- N6 ref (atoms/reduction; ⛔ not imported): `.../scripts/S11c_c2_N6_diagnostic_sympy.py` (`pressure_coefficients` `:410-412`).

## Verify each round-2 fold is CORRECT
1. **K_A SymPy FORM is now a structural addend deletion** (directive `:110-115`): delete the `closure_shape_deriv`
   addends containing the inherited `Lambda_A_0` symbol at `:2815-2828`, leaving `Lambda_V_0` addends; ⛔ NOT
   `Lambda_A_0→0`, ⛔ NOT the `:408` bind. Verify: (a) the Λ_A and Λ_V addends are **cleanly separable** (no addend
   contains both), so the deletion is well-defined; (b) the deletion is a genuine **FORM** (leaves the family),
   distinct from the ×2 rescale used in self-test 3; (c) it moves the θ carrier and leaves Λ_V / U / e_W untouched.
2. **Directive blindness** (the whole builder packet): does it still leak any expected **value / count / ratio /
   pass-condition / target-row / must-MOVE / must-NOT** anywhere? The knife headings should name the SITE, not the
   row it feeds; the target-object and scope sections should not map channel→row; the "no builder-chosen knives"
   clause (`:147`) should read as implementation-fidelity, not a bite expectation. Confirm the harness is told to
   print the COMPLETE carrier (all rows × all slots) with no selection.
3. **Self-tests pinned** (directive `:150-159`): dead-path = delete a pressure-free kinetic addend (SymPy `:2970`
   `e_kinetic`; WL `:1347` `kineticEwLive`) — confirm both are **pressure-free** so the dead-path carrier diff is 0;
   live-rescale = ×2 at each knife's own site (K_A Λ_A addends, K_T face additions/`virtualWork`, K_W minus pressure
   slot) — confirm each names a real site and that ×2 differs from that knife's FORM.
4. **K_W single function** (directive `:133-143`): pinned to `substrate_substitutions` `:1995` (consumed only by
   `filtered_substrate` `:2032-2048`). Confirm one function, and that the collapse reaches the delta_p-bearing
   `closure_shape_deriv` / `virtual_work_shape_deriv`.
5. **Adjudication key K_W cone** (key `:56-71`): the collapse is `C_plus ← C_plus+C_minus`, `C_minus ← 0`; observed
   as minus columns → 0, plus VALUE column doubles (±-value symmetric), plus `d_w` column **cancels** (±-`d_w`
   antisymmetric), NO full-matrix rank change (base and collapsed column ranks both 1). Verify these symbolic
   relations independently (the ±-value symmetry `p_minus − p_plus ≡ 0` and the ±-`d_w` antisymmetry
   `d_w_minus + d_w_plus ≡ 0`). Confirm the key no longer says "plus absorbs both" or "rank-changing."

## Also confirm (round-2 PASSES that must remain true)
WL definitions-only load + `evaluatedModel["EULERIAN","MATERIAL_ADVECTED","RHO4_CONSTANT"]["OPERATOR"]` (never the
emit `Do`/`extractCouplingData`); MATERIAL_ADVECTED makes U a live carrier in both engines; K_T whole-channel
closes the Λ_X hole; SymPy `named_tuple_row` access (not `["…"]`); disjointness (neither K_A nor K_T patches the
shared `affinity` `:1079`); three PRINT-not-assert clauses intact.

## Output
Per-item verdict with cited lines, then a final line: **SOUND** (ready to build) or **NOT-SOUND** (exact site + fix
for every blocking issue).
