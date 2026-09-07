# S11c-b carrier ablation-harness knife-list — RE-AUTHOR review record (CLEARED, 2026-09-07)

**Artifact (accepted):** `directives/S11c_b_carrier_ablation_harness_directive.md` (blind builder packet) +
`_measurements/S11c_b_carrier_ablation_harness_adjudication_key.md` (orchestrator adjudication key — the cones /
expected bites, ⛔ NOT handed to the builder). Author: orchestrator. **Legs (G1 orchestrator-written → Codex +
Grok):** codex-sol (gpt-5.6-sol xhigh) + Grok (grok-4.6 high), identical rendered prompt each round
(`directives/_legs/S11c_b_carrier_reauthor_review_prompt{,_round2,_round3,_round4}.md`). **Method:**
review-until-clear (a wrong-object knife-list is a physics defect; ⛔ not one-pass). Raw leg logs (ephemeral
scratch): `scratchpad/s11cb_reauthor/{,r2/,r3/,r4/}{sol,grok}.log`. This supersedes the NOT-SOUND G2-gate baseline
`dc4d4977` (`_measurements/S11c_b_carrier_ablation_harness_knife_list_G2_gate.md`).

## Trajectory (converged: 13 blocking → 5 → 3 precision → 0)
- **Round 1 — both NOT-SOUND.** Central re-author CONFIRMED correct (carrier is face-law-sourced; energy +
  constraint pressure-free ⇒ old K1/K2 rightly dropped). ~5 blocking impl defects: WL entrypoint OOM, K_T two-FORM
  inconsistency, K_A zeroing kills Λ_V, K_W exchange dead, pinned-case U pressure-free, SymPy dict access, wrong
  observation citations, drift guard, pass-condition leaks into the builder packet.
- **Round 2 — both NOT-SOUND (narrow).** Whole design PASSED both legs; blocking: K_A `Lambda_A_0→0` reads as a
  coefficient (self-test can't distinguish FORM from arithmetic); K_W key cone "plus absorbs both" wrong; residual
  blindness leaks; self-tests unpinned; K_W named two functions.
- **Round 3 — grok SOUND, sol NOT-SOUND (3 items).** K_A needs `sp.expand` before the addend filter; two stale
  `Lambda_A_0→0` in the key; a few remaining directive-side cone disclosures.
- **Round 4 — BOTH SOUND.** All three round-3 folds confirmed against the live engines; every must-remain-true
  property re-confirmed. **CLEARED.**

## Grounding evidence (from the legs' committed check-script stdout — the design rests on measurement, not prose)
- **Carrier is face-law-sourced (drops of energy/constraint correct):** `ENERGY_PRESSURE_DEPS (F,F,F,F)`,
  `CONSTRAINT_PRESSURE_DEPS (F,F,F,F)`, `UNFOLDED_MASS_PRESSURE_DEPS (F,F,F,F)`, `CLOSURE_PRESSURE_DEPS
  ((T,F,T,F),(F,T,F,T))`, `VIRTUAL_WORK_PRESSURE_DEPS (T,T,T,T)`; WL `ENERGY_PRESSURE_FREE=True`,
  `CONSTRAINT_PRESSURE_FREE=True`.
- **Case pin MATERIAL_ADVECTED (U live) vs LAB_HELD (U dead):** SymPy `U_CARRIER_LAB_NZ [(F,F,F,F)×3]`,
  `U_CARRIER_MAT_NZ [(T,T,T,T)×3]`, `EW_CARRIER_MAT_NZ (T,T,T,T)`, `THETA_SUM_NZ (T,T,T,T)`; WL
  `WL_U_CARRIER_LIVE=True`. (LAB_HELD virtual normal displacement has no virtual-U, `WL:1040-1052`.)
- **K_A addend deletion (needs `sp.expand`; genuine FORM):** unexpanded `closure_shape_deriv` is one `Mul`
  (`NAIVE_PATCHED_IS_ZERO True` — naive filter drops everything); expanded → `LA_N 4, LV_N 4, SHARED_N 0`,
  `PATCHED_HAS_LA False, PATCHED_HAS_LV True`, and `RESCALE_EQ_PATCHED False` (×2 ≠ deletion).
- **K_T whole channel closes the Λ_X hole:** `EW_BARE_P_SUBCHANNEL_NZ (T,T,T,T)` AND `EW_LAMBDAX_SUBCHANNEL_NZ
  (T,T,T,T)`; `VW_HAS_LX True`. Removing all four `face_u`/`face_e` / `virtualWork=0` covers both.
- **K_W collapse (exchange is dead):** `P_MINUS_MINUS_P_PLUS_ALL_ZERO True` (±-value symmetric),
  `DW_MINUS_PLUS_DW_PLUS_ALL_ZERO True` (±-`d_w` antisymmetric) ⇒ collapse gives `COLLAPSE_MINUS_COLUMNS_ALL_ZERO
  True`, `COLLAPSE_PLUS_VALUE_DOUBLES True`, `COLLAPSE_PLUS_DW_ALL_ZERO True` (cancels, does NOT absorb),
  `BASE_COLLAPSED_COLUMN_RANKS 1 1` (⛔ no full-matrix rank change). Exchange `MATERIAL_SWAP_COLUMN_DIFF_NZ
  {False,False}` (WL) / `(F,T,F,T)` (SymPy d_w only).
- **WL entrypoint avoids OOM:** emit `Do` (`:2190`) calls `extractCouplingData` (`:2199`) +
  `kernelOriginsFromOrigins` (`:2203`) BEFORE the `PRIMARIES_ONLY` test (`:2206`) ⇒ definitions-only load + direct
  `evaluatedModel["EULERIAN","MATERIAL_ADVECTED","RHO4_CONSTANT"]["OPERATOR"]` (3-arg ⇒ `corrupted=False`).
- **SymPy access:** `build_operator` returns `casify(...)` nested `Tuple` (`CASIFY_IS_TUPLE True`,
  `DICT_ACCESS_WORKS False`) ⇒ `named_tuple_row`. **Dead-path sites pressure-free:** `E_KINETIC_PRESSURE_MASK
  (F,F,F,F)`, WL `kineticEwLive` no `pressureField`.

## Final design (accepted)
Two committed harnesses wrap the live engines (SymPy import + `build_operator`; WL definitions-only + `evaluatedModel`),
single case **EULERIAN × MATERIAL_ADVECTED × RHO4_CONSTANT**, and form the carrier `∂(emitted rows)/∂(native pressure
atoms)` then `→0`. Three knives — **K_A** (delete expanded Λ_A closure addends), **K_T** (remove the whole traction
channel), **K_W** (native-pressure face collapse) — plus three self-tests (extractor-order, dead-path, ×2 rescale).
The harness PRINTS the complete carrier (all rows × all slots) with no selection; interpretation (the cones in the
adjudication key) is the orchestrator's. ⛔ Builder gets ONLY the directive.

## Disposition
**CLEARED to build** (both round-4 legs SOUND). NEXT: a fresh `gpt-6-astra` builds the two harnesses from the
directive ALONE; then fresh-Claude + Grok review-until-clear (ablate the harness itself); then commit + run +
record. The N6 4-harness backfill and the c2 T7 comparator follow.
