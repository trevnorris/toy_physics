# Measurements — Grok leg, S11c-a comparator directive rev 3

Commands were ad-hoc Python against committed transcripts + `~/.s11_build/S11c_a_scratch_loader.py`.
Literal counts below are from those runs (session terminal logs).

## Tag inventory
- PY tags: 47 (8 `_LOCAL_*`); WL tags: 40 (1 `LOCAL_TAG_NAMES`) → 39 joinable families.

## VIRTUAL_WORK keys (literal)
PY: `(BRANCH, DELTA_W|ZETA_C, DELTA_W|ZETA_C, RHO*)` ×16 — positional DOF then VDOF.
WL: `(BRANCH, RHO*, DOF_*, VIRTUAL_DOF_*)` ×16.
Scratch `canon_key` frozenset-of-untyped-DOF join: **0/16** (dup_py=4).
Typed `(DOF,·)/(VDOF,·)` frozenset join: **16/16** (dup_py=0, dup_wl=0).
PY VALUE: labeled `UPPER`/`LOWER` pair (raw srepr).
WL leaf: single `SHAPE_DERIVATIVE.EXPRESSION` — **no UPPER/LOWER field**.

## CONTROL_FORM keys (literal)
PY cases (**no** `Str('VALUE')` wrapper; direct key→payload): **528** keys; lengths {5:144, 6:336, 7:48}.
Example len5: `('FACE_NORMAL','LAB_HELD','1','DELTA_W','1')` — last token ∈ `{1,2,3}` (direction), **not** ±1.
Face-like slot tokens: `{1,-1,BOTH_FACES}`; `BOTH_FACES` in **240** keys.
WL cases: **960** unique; OBJECT + optional ORIGIN/FIELD/VDOF + `DIRECTION_1/2/3`; **288** without `FACE_PLUS`/`FACE_MINUS`.
Naive typed join across all objects: **144** (py_only 384, wl_only 816).
`FACE_NORMAL` denseless 5-token PY vs WL: PY 8, WL denseless subset join incomplete (WL often denseless for this object).

## FACE_SHIFT
PY 16 aggregate `(BRANCH,±1,DOF,DENSITY)`; RHO4 vs RHOBR VALUE strings: **identical** (8/8 twins).
VALUE Tuple arity 4: [0] pressure scalar, [1] 4×1 velocity, [2] density scalar, [3] 4×1 currents
  `[3][0..2]=…delta_j_bulk_{1,2,3}`, `[3][3]=…delta_j_bulk_4` (bare Symbols inside).
WL 80 `(BRANCH,FACE,FIELD_*,DOF)` — no density; 10 FIELD_* names. 8 bases × 10 fields = 80 **iff density dropped**.

## Currents (payload tokens)
PROJECTION_SHAPE PY: `Function('delta_j_bulk_4')(Symbol('w'))` ×8; in-plane bare `Symbol('delta_j_bulk_{1,2,3}')` ×16 each.
FACE_SHIFT PY: bare `Symbol('delta_j_bulk_4')` ×16; bare in-plane ×16.
WL PROJECTION: `currentWPerturbation[x1, x2, x3, {normalCoordinate, time}]` (count 40 on that tag);
  `currentXPerturbation{1,2,3}[same]` (96 each) — **no underscore** before the index.
WL FACE_SHIFT: same 4-arg field forms (`currentW` ×24, `currentX*` ×72).
`normalCoordinate` count on WL PROJECTION_SHAPE_DERIV: **2140**.

## KINEMATIC WL leaves (verified)
`OPERAND_A_SHAPE_DERIVATIVE` / `OPERAND_B_SHAPE_DERIVATIVE` / `RESIDUAL_A_MINUS_B` each wrap `"EXPRESSION" -> …`.

## HOMOGENEITY leaves (verified)
PY objects: `RELATIVE_FLUX_LAW`, `KINEMATIC_BALANCE`, `TRACTION`, `VIRTUAL_WORK`, `PROJECTION`, `VIRTUAL_CONSTRAINT`, `EVOLUTION_BALANCE`, `CLOSURE` — not VALUE-case maps.
WL `HOMOGENEITY_BASE_OPERAND` / RELATIVE_FLUX fields include `SOURCE_TERM_DIMENSIONS`, `HOMOGENEITY_OBJECT` (not `BASE_OPERAND`).
WL `HOMOGENEITY_RESIDUAL` / RELATIVE_FLUX fields: `BASE_OPERAND`, `CONTROL_OPERAND`, `RESIDUAL`, `TEST_OBJECTS`.

## Axes checked sound (with token)
- FACE_VELOCITY / CLOSURE / PROJECTION_* case-key axis-typing with named BRANCH/FACE±1/DOF/DENSITY: scratch-canon join 8/8 or 16/16.
- FACE_SHIFT flatten index layout matches VALUE arity-4 as above.
- Dwin LHS→RHS in directive matches `S11c_a_run_projection.py` window-only rewrites.
- CORRECTION 2 rejects reference `ITG`/`W_INT` bound-stripping (`run_projection.py` `canon_ints`).
- MUMAP / CONORMAL Taylor excluded as required by §3a branch `μ_θ^α` / raw CONORMAL.
- No rule-5 expected residual target found in the directive text for MEASURED families.
- S11b reuse ban on `compare_records` / `FINAL_OPERATIONAL_STATUS` matches frozen script (`S11b_cross_engine_comparator.py` ~L999, ~L1376).
