# Decision-list review — S11c-b P2 cross-engine comparator update (extract_slab + §3a bridge)

## Artifact
`research/pde_ledger_v3/directives/S11c_b_p2_comparator_update_directive.md` — an orchestrator-written DECISION LIST
for a PHYSICS-BEARING comparator change (a mis-map manufactures false cross-engine agreement OR disagreement). You
are ONE of two independent decision legs (rule 7). The builder will TRUST this list. Find its defects now — a
productive review FINDS things.

## What the change is for
Since the comparator was written, the PY engine folded (constraint-reduced U/E_W rows; θ-row became the
mass-evolution balance `evolution_mass_balance − Σ closure_shape_deriv`, NOT μ_θ; μ_θ a separate operand; face
generalized-force rows) and #90 put face+response into the operator; the WL engine renamed its §3a coefficients to
`gammaWJET*`/`gammaMUJET*`. The comparator's `extract_slab` still maps `THETA_BALANCE → ("MU_THETA","THETA")` and its
§3a bridge tables still key on the old `gammaWidth*` names — so today it manufactures a spurious disagreement on
μ_θ/mass-evolution and 0 of ~30 §3a coefficients are actually bridged. The update must make the comparator join the
CURRENT structure correctly, WITHOUT manufacturing agreement (a wrong bridge / blanket-collapse) or disagreement (a
stale bridge joining unlike objects).

## Context you are handed (read the CODE, cite file:line)
- Comparator `research/pde_ledger_v3/scripts/S11c_b_cross_engine_comparator.py`: `extract_slab` (~760, the PY
  `row_map` at ~766-770 and WL side ~785-806), `extract_coupling` (~850), `extract_mu` (~826), `extract_energy`
  (~736), name tables `EXTRA_SYMBOL`/`_extra_basic` (~349-429).
- Bridge tables `research/pde_ledger_v3/scripts/S11c_b_handcoded_comparison.py` (`WL_TO_PY_RENAME` ~80-152, the
  `gammaWidth*→gamma_s11cb_*` entries ~107-118) and `S11c_b_adjudicated_comparison.py` (`PROTECTED_ATOM_NAMES` ~54-65).
- PY engine emit `research/pde_ledger_v3/scripts/S11c_b_brane_operator_sympy_audit.py`: `build_operator` (~2918-3074,
  the post-fold `operator` dict keys `U_BODY_BALANCE`/`E_W_BALANCE`/`THETA_BALANCE`/`FACE_GENERALIZED_FORCE_ROWS`/
  `MU_THETA_FACE_BINDING`), `NEW_COEFFICIENTS` (~357-367, `gamma_s11cb_{source}_{NN}`), `ENERGY_BASIS_NEW_INVARIANTS`
  emit (~4101), the runtime quotient selection (`quotient_independent_indices` ~1406-1439).
- WL engine emit `research/pde_ledger_v3/mathematica/S11c_b_brane_operator_mathematica_audit.wl`: `modelRecord`
  (~1570, the `ROWS`/`DIVERGENCE_FORM_SOURCE`/`MASS_EVOLUTION_ROW`/`CENTER_FACE_GENERALIZED_ROW` structure),
  `newCoefficientSymbol` (~627, `gamma <> "WJET"|"MUJET" <> suffix`), the new-invariant records (~2139-2145).

## What to check (find defects; ground every claim in code, cite file:line)
1. **extract_slab structural completeness (D1-D6).** Is the mapping between the post-fold PY `operator` keys and the
   WL `modelRecord` keys COMPLETE and CORRECT? Independently confirm: PY `THETA_BALANCE.EXPANDED` is now mass-evolution
   (grep `build_operator`/the θ-row assembly) and its WL partner is `MASS_EVOLUTION_ROW` (not `DIVERGENCE_FORM_SOURCE.
   MU_THETA`). Is any PY or WL slab sub-object left an unjoined orphan the directive does not name? Does D2's
   "μ_θ only through MU_THETA_OPERATOR" actually avoid a double-count?
2. **§3a bridge soundness (D7-D10) — the manufactured-agreement danger.** Is D8's "match the invariant, never guess
   the index" implementable from what the engines EMIT? Confirm PY emits the invariant per coefficient
   (`ENERGY_BASIS_NEW_INVARIANTS`) and WL emits its new-invariant records, and that the profile-jet/DOF renames
   needed to compare the two invariant expressions are already complete. Is there a residual/one-sided-corruption
   test that would CATCH a wrong pairing (D10)? Could a builder satisfy the directive with a positional-index guess
   that the legs would not catch?
3. **extract_coupling (D6b).** The directive flags this as OPEN. Independently inspect whether the #90 face-force /
   closure changes altered the kernel structure `extract_coupling` reads. Report whether it needs updating.
4. **Scope / case-set (D11).** Both engines now emit 4-case primaries (matching case-sets). Confirm `row_residual`'s
   union+align consumes matching case-sets without raising. Is anything the directive assumes (matching cases, no
   intersection needed) actually true of the committed loader?
5. Any way the update could MANUFACTURE agreement (a blanket applied→bare collapse, a rename that folds a real
   `W_bg`/`mu_R_bg` background-jet dependence, a quotient representative folded as physics) or MANUFACTURE
   disagreement. Any missing property, any ambiguity a builder could satisfy wrongly.

## Required method
DOCUMENT review — do NOT modify the tree. Ground every claim in code (file:line); where checkable by a command
(grep the PY θ-row assembly; grep the WL new-invariant emit; grep the bridge tables), RUN it and quote the LITERAL
output — a prose assertion is discarded. Report defects as a numbered list; for each name the directive item (D#) it
fixes and the file:line evidence. A clean pass is weak evidence.
