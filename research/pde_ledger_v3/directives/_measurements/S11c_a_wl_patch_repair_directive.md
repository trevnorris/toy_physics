# Measurements twin — S11c_a_wl_patch_repair_directive.md

Commands establishing every current-state claim + the provenance verification. rule-5: structural only.
CONVENTIONS: repo root `/var/projects/toy_physics`;
`WL=research/pde_ledger_v3/mathematica/S11c_a_interface_geometry_mathematica_audit.wl` (patched, uncommitted);
`OUT=research/pde_ledger_v3/mathematica/out/S11c_a_interface_geometry_mathematica_audit.out` (committed baseline 82062bd).

## materialShape has no reset indirection; directShape does
CMD `sed -n '595,597p' "$WL"` -> `materialShape[…] := materialShape[…] = shapeDerivative[materialFaceObjectsSource[…]]`
(memoizing, single definition, no `materialShapeDefinition` indirection).
CMD `sed -n '582,592p' "$WL"` -> `resetDirectShapeCache[]`/`resetDirectShapeMemoCache[]` each `Clear[directShape]`
then REDEFINE it. No analogue exists for materialShape.

## materialShape is Cleared and never redefined (patch-introduced)
CMD `grep -nE 'Clear\[' "$WL"` + read each block -> two `Clear[…, materialShape]` sites (~1456 and ~1973), each
followed by `resetDirectShapeCache[]` (redefines directShape only), never a materialShape redefine.
CMD `sed -n '1332p' "$WL"` -> the §5a object loop consumes `materialShape[branch, sign, density, …]` (inert after the clear).
CMD `git diff "$WL" | grep -nE '^[+].*(materialShape|resetDirectShapeCache)'` -> the patch ADDED the reset block
incl. `Clear[…, materialShape]` + `resetDirectShapeCache[]`.

## Committed baseline computed the material route (patch broke it)
CMD `grep -c 'materialShape\[' "$OUT"` -> 0 (no inert token in the baseline output).
CMD `grep -o 'WL_S11CA_REP_INVARIANCE_MATERIAL_OPERAND[^Z]*' "$OUT" | head -c 400` -> shows a fully computed
expression (e.g. VIRTUAL_CONSTRAINT operand), not an inert symbol.

## Spec anchor
CMD `sed -n '452,474p' "$SPEC"` (SPEC=research/pde_ledger_v3/directives/S11c_a_SHARED_PHYSICS.md) -> §5a: two
independent routes, emit EULERIAN/MATERIAL operands + residual; one-sided independence corrupts one route.

## Discipline
Repair relays the single blocking finding verbatim (rule 15) + points at the closed spec §5a. Blind-safe (only
the spec cited); rule-5 (no expected residual). PRESERVE the three verified patch changes.
