# Measurements twin — S11c_a_wl_patch_directive.md

Every factual claim the directive makes about the WL engine's CURRENT state carries the command that produced
it. ⚠ rule-5: no output-expression value appears here; claims are structural (which axis a key carries, which
objects a control loop iterates), ⛔ never a computed value.

CONVENTIONS (all commands): run from repo root `/var/projects/toy_physics` with
`WL=research/pde_ledger_v3/mathematica/S11c_a_interface_geometry_mathematica_audit.wl`,
`SPEC=research/pde_ledger_v3/directives/S11c_a_SHARED_PHYSICS.md`. WL engine at HEAD = committed `ddecdbc2`;
spec `2926c71c`.

## Claim: projection cluster (T-f) carries NO density-representative axis
CMD: `sed -n '758,760p' "$WL"` -> `projectionTermsSource[branch_String, dynamic_:True]` takes branch +
dynamic only, no `density` parameter; its density is `rhoBulkZero`/`rhoBulkWave` (abstract), not
`rho4Profile[density,...]`.
CMD: `grep -n 'emitShared\["PROJECTION_' "$WL"` -> the five T-f emitters
(SHAPE_DERIV/DYNAMIC_OPERAND/STATIC_OPERAND/RESIDUAL/TERM_ORIGINS).
CMD (keys carry branch|DOF, no density token):
`sed -n '1081,1136p' "$WL" | grep -nE 'key = |branch <> "\|'` -> keys are `branch <> "|DOF_" <> dof`
(and `branch <> "|" <> origin <> "|DOF_" <> dof` for TERM_ORIGINS); no `density` in any projection key.

## Claim: KINEMATIC_BALANCE and RELATIVE_FLUX DO carry a density-representative axis
CMD: `sed -n '988,1001p' "$WL"` -> `RELATIVE_FLUX` payload keyed `densityFaceCaseKey[branch, density, sign,
dof]`, iterating `{density, densityNames}`.
CMD: `sed -n '1003,1024p' "$WL"` -> `KINEMATIC_BALANCE` payload keyed `densityFaceCaseKey[branch, density,
sign, dof]`, iterating `{density, densityNames}`.
CMD (flux/kinematic use rho_m): `sed -n '476,478p;512p' "$WL"` -> `kinematicOperandB = normalVelocity +
flux/rhoM`; flux enters via `rhoM`, not a representative.
CMD (representative names): `grep -n 'densityNames = ' "$WL"` -> `{"RHO4_CONSTANT", "RHOBR_CONSTANT"}` (line 108).

## Claim: §5b form control covers 13 objects; projection = SHAPE_DERIV only; misses the 5
CMD (the object-lists the form-ablation loops iterate):
`sed -n '1377,1509p' "$WL" | grep -nE '\{object, \{|"FACE_SHIFT|"PROJECTION_SHAPE_DERIV\|"|VIRTUAL_WORK_SHAPE_DERIV\|"|"VIRTUAL_CONSTRAINT|"EVOLUTION_MASS_BALANCE'`
-> loops cover: {FACE_NORMAL, CONORMAL_DERIV, FACE_MEASURE_SHAPE_DERIV, FACE_VELOCITY}; {RELATIVE_FLUX,
KINEMATIC_BALANCE, TRACTION, CLOSURE_SHAPE_DERIV}; FACE_SHIFT; VIRTUAL_WORK_SHAPE_DERIV; PROJECTION_SHAPE_DERIV
(only projection object); {VIRTUAL_CONSTRAINT, EVOLUTION_MASS_BALANCE}. = 13. ABSENT: EVOLUTION_TERM_ORIGINS,
PROJECTION_STATIC_OPERAND, PROJECTION_DYNAMIC_OPERAND, PROJECTION_RESIDUAL, PROJECTION_TERM_ORIGINS.

## Claim: §5c uniform-limit covers the same 13; projection = SHAPE_DERIV only; misses the same 5
CMD:
`sed -n '1537,1659p' "$WL" | grep -nE '\{object, \{|"FACE_SHIFT|"PROJECTION_SHAPE_DERIV\|"|VIRTUAL_WORK_SHAPE_DERIV\|"|"VIRTUAL_CONSTRAINT|"EVOLUTION_MASS_BALANCE|BACKGROUND_DENSITY_MAP'`
-> same object set as §5b (+ BACKGROUND_DENSITY_MAP per representative); projection = PROJECTION_SHAPE_DERIV
only, keyed `branch <> "|DOF_" <> dof`. Same 5 absent.

## Spec anchors cited by the directive
CMD: `sed -n '74,83p' "$SPEC"` (§1b: j = ρ_4D v_bulk, projection = integrate conservation law vs Ω).
CMD: `sed -n '228,230p' "$SPEC"` (§2b: construct the two representatives; don't simplify one into the other
before T-f/T-g/T-h/T-i).
CMD: `sed -n '312,318p' "$SPEC"` (§3a: ρ_4D^α ≡ ρ_4D,bg^{0,α}(1+θ), representative/anchoring dependent).
CMD: `sed -n '349,353p' "$SPEC"` (§3b: J_s^α ≡ ρ_m(...); n̂·v = V + J/ρ_m — flux/kinematic use ρ_m).
CMD: `sed -n '397,400p' "$SPEC"` (§4: "both density representatives wherever density enters").
CMD: `sed -n '424,429p' "$SPEC"` (§4 T-f: the five projection objects).
CMD: `sed -n '487,492p' "$SPEC"` (§5b: "For every T-object and each direction ...").
CMD: `sed -n '494,499p' "$SPEC"` (§5c: "For each S11c-a object ...").

## Discipline
Directive relays the step-1 verdicts (B projection + H.1 coverage + kinematic/flux drop) as spec-grounded
schema/coverage changes, ⛔ never as "match the sibling" and ⛔ never with an expected value. Blind-safe: only
the closed spec is cited; no sibling construction/result appears.
