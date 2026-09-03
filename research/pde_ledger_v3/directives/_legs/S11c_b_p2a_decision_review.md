# Decision-list review — S11c-b P2a slab-row join + `row_residual` structural fixes

## Artifact
`research/pde_ledger_v3/directives/S11c_b_p2a_slab_row_join_directive.md` — orchestrator-written DECISION LIST for a
PHYSICS-BEARING comparator + residual change. You are ONE of two independent decision legs (rule 7). The builder
will TRUST this list. Find its defects now — a productive review FINDS things.

## Context
The PY engine folded (θ-row = mass-evolution `evolution_mass_balance − Σ closure_shape_deriv`; μ_θ separate; face
generalized-force rows) and #90 put face+response into the operator. The comparator's `extract_slab` and
`row_residual.py` were written before this. A prior 2-leg review already found: the stale `THETA_BALANCE→MU_THETA`
map, the #90-stale one-sided face subtraction in `row_residual.py:427`, μ_θ not consumed by `row_residual`, the
duplicate-axis raise risk, and orphan keys. This directive (A1-A11) is the fix. Your job: is it COMPLETE and
CORRECT, or can a builder still manufacture agreement/disagreement or crash the loader?

## Context you are handed (read the CODE, cite file:line)
- `scripts/S11c_b_cross_engine_comparator.py`: `extract_slab` (~760-806, PY `row_map` ~766, WL side ~785-806).
- `scripts/S11c_b_row_residual.py`: `SLAB_OBJECTS` (~43), `_load_aligned_cases`/`_family_cases` calls (~551-602),
  the face subtraction (~427, ~721), the hardcoded default paths (~1084/1095), cardinality guards (~1171).
- `scripts/S11c_b_adjudicated_comparison.py`: the union+align loader (~1098), the duplicate-axis raise (~1091-1096).
- PY engine `scripts/S11c_b_brane_operator_sympy_audit.py`: `build_operator` (~2918-3074), THETA_BALANCE assembly
  (~3042-3046), FACE_GENERALIZED_FORCE_ROWS (~3051), the origin family (~3195-3246).
- WL engine `research/pde_ledger_v3/mathematica/S11c_b_brane_operator_mathematica_audit.wl`: `modelRecord` (~1570),
  `MASS_EVOLUTION_ROW`/`CENTER_FACE_GENERALIZED_ROW` (~1173-1177), the face rows in the complete rows (~1345).

## What to check (find defects; ground every claim in code, cite file:line; run greps and quote literal output)
1. **A1-A3 (θ/mass/μ_θ):** does the directive fully prevent (a) double-counting the θ row (SOURCE_OPERAND ≡
   EXPANDED), (b) the duplicate-axis raise from A2's ADVECTIVE line, (c) the one-sided SLAB `MU_THETA` raise from
   A3? Is any residual path to a canonical `MU_THETA` on ONE engine still open?
2. **A4-A6 (face / center-face / disposition):** is the CENTER_FACE sign reconciliation specified enough that a
   builder joins the SAME virtual displacement (not opposite-signed)? Is the double-count prohibition (don't re-join
   FACE U/E_W already in the rows) airtight? Is the A6 disposition table actually exhaustive — enumerate every PY
   and WL slab sub-object yourself and check each is dispositioned; report any the directive misses.
3. **A7 (the #90 face subtraction):** is "face stays in both complete rows / attributed symmetrically" the correct
   fix? Independently inspect `row_residual.py:427` + `:721` — confirm both engines now carry face in the complete
   rows, so the one-sided subtraction is the bug. Could A7 as written under- or over-correct (e.g. drop face
   provenance the instrument needs)?
4. **A8-A9 (μ_θ load + guards):** does A9's expected-set guard actually prevent a silently-missing case (the loader
   only unions)? Is loading `MU_THETA_OPERATOR` into `row_residual` correctly specified (which comparison type)?
5. Any way the change MANUFACTURES agreement (a blanket collapse, a folded background-jet dependence) or
   DISAGREEMENT (a stale/one-sided operand, an unreconciled sign). Any missing property, any ambiguity.

## Required method
DOCUMENT review — do NOT modify the tree. Ground every claim in code (file:line); run greps and quote LITERAL
output. Report defects as a numbered list; name the directive item (A#) each fixes. A clean pass is weak evidence.
