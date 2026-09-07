# Independent review ROUND 4 (confirmation) — S11c-b carrier ablation-harness knife-list

You are one of two independent review legs. Working dir: `/var/projects/toy_physics`. Repo read access; you may
write/run small check scripts under `/tmp`. ⛔ Do not modify the working tree. ⛔ You are not handed the other leg's
output. **Prose is discounted unless grounded in a cited engine line or a check script whose path + literal stdout
you report.**

## Context
Round 3 was one-leg-SOUND, one-leg-NOT-SOUND on three precision items. All three are now folded. This is a tight
**confirmation**: verify the three fixes are correctly applied, the directive is still implementable and blind, and
no new defect was introduced. Report a finding only if it catches a way the harness would certify the wrong thing
or fail to run.

## Artifacts (the pair)
- Blind builder directive: `/var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_b_carrier_ablation_harness_directive.md`
- Orchestrator adjudication key: `/var/projects/toy_physics/research/pde_ledger_v3/_measurements/S11c_b_carrier_ablation_harness_adjudication_key.md`
- Engines: `.../scripts/S11c_b_brane_operator_sympy_audit.py`, `.../mathematica/S11c_b_brane_operator_mathematica_audit.wl`; N6 ref `.../scripts/S11c_c2_N6_diagnostic_sympy.py`.

## Confirm the three fixes
1. **K_A `sp.expand` fix** (directive `:108-116`). The FORM now expands each face's `closure_shape_deriv` before
   deleting the `Lambda_A_0`-bearing terms:
   `patched = sp.Add(*(t for t in sp.Add.make_args(sp.expand(closure_shape_deriv)) if not t.has(Lambda_A_0)))`.
   Confirm against the engine: (a) each `closure_shape_deriv` is an unexpanded `Mul` so the expand is **required**
   (without it the filter drops the whole expression); (b) after expand, the Λ_A / Λ_V addends are cleanly
   separable (no shared addend) so the deletion keeps Λ_V and is a genuine **FORM**; (c) the K_A ×2 self-test
   (directive `:154-156`) now rescales the *same expanded* Λ_A addends, so it is distinct from the deletion FORM.
2. **Blindness tightening** (directive). The target-object/scope no longer maps channel→row or calls
   energy/constraint pressure-independent; K_A/K_T no longer name the terminal emitted row they feed; the K_T "still
   pressure-bearing / uncovered" pass-rationale is gone. Confirm no expected **value / count / ratio /
   pass-condition / must-MOVE / must-NOT / target-row-cone** remains in the builder packet — AND that every knife
   still names a **locatable site** (function + line) so it is implementable. (The observation section legitimately
   lists all rows to print — that is the complete-carrier target, not a cone.)
3. **Key consistency** (adjudication key). The stale "null `Lambda_A_0`" / "`Lambda_A_0→0`" descriptions of the K_A
   FORM are replaced by "delete the expanded Λ_A-bearing addends" (the only remaining `Lambda_A_0→0` mentions are in
   the round-3 change-note, which correctly describes the sharpening *away* from zeroing). Confirm the key's K_A FORM
   is now consistently the structural deletion.

## Also confirm (must remain true)
K_A/K_T disjointness (neither patches WL `affinity` `:1079`); K_T whole-channel (`virtualWork=0` / four SymPy
additions) covers both bare-p and Λ_X; MATERIAL_ADVECTED keeps U live; WL definitions-only + `evaluatedModel[...]
["OPERATOR"]` avoids the OOM path; SymPy `named_tuple_row` access; K_W collapse + its key cone (minus→0, plus-value
doubles, plus-`d_w` cancels, no full-matrix rank change); three PRINT-not-assert clauses intact; dead-path sites
`:2970`/`:1347` pressure-free.

## Output
Per-item verdict with cited lines, then a final line: **SOUND** (ready to build) or **NOT-SOUND** (exact site + fix).
