# Measurements twin — S11c_a_T7_adjudication_matrix.md

CONVENTIONS (all commands below): run from repo root `/var/projects/toy_physics` with
`PY=/home/trevnorris/.s11_build/S11c_a_sympy_engine.out`,
`WL=research/pde_ledger_v3/mathematica/out/S11c_a_interface_geometry_mathematica_audit.out`,
`SPEC=research/pde_ledger_v3/directives/S11c_a_SHARED_PHYSICS.md`. Producer-citation lines that abbreviate
a path with `…` mean the full engine source: PY `scripts/S11c_a_interface_geometry_sympy_audit.py`,
WL `mathematica/S11c_a_interface_geometry_mathematica_audit.wl` (both under `research/pde_ledger_v3/`).

Inputs (hash-locked):
```
sha256  6386471555b1e99d0aeb0f716eea30f839d59be50c0cedd4677ea7b376b79129  ~/.s11_build/S11c_a_sympy_engine.out   (PY, fresh run of 9b6438fa)
sha256  82062bd36cfb07b1f18631077f0c63ac1cbce7834967686f680fa9f30019e4ec  research/pde_ledger_v3/mathematica/out/S11c_a_interface_geometry_mathematica_audit.out  (WL, ddecdbc2)
```
CMD (hashes): `sha256sum ~/.s11_build/S11c_a_sympy_engine.out research/pde_ledger_v3/mathematica/out/S11c_a_interface_geometry_mathematica_audit.out`

## Tag join (39 exact)
CMD: `comm -12 <(grep -oE '^PY_S11CA_[A-Z][A-Z0-9_]*' "$PY" | sed 's/^PY_//' | sort -u) <(grep -oE '^WL_S11CA_[A-Z][A-Z0-9_]*' "$WL" | sed 's/^WL_//' | sort -u) | wc -l`
⇒ 39. PY-only: 8 `*_LOCAL_*`; WL-only: 1 `*_LOCAL_TAG_NAMES`. One line per tag both sides
(`grep -cE '^PY_S11CA_' = 47 = non-blank total`; WL `40 = total`).

## The AXIS-SET census (backbone of every count/axis claim in the matrix)
Apparatus (frozen, committed): `_measurements/s11ca_t7_census/s11ca_axis_census.py`.
CMD: `python3 _measurements/s11ca_t7_census/s11ca_axis_census.py`
Literal stdout (frozen): `_measurements/s11ca_t7_census/axis_census.stdout`; leaf keys:
`_measurements/s11ca_t7_census/axis_census.keys.txt`. The census flattens PY quantity-nesting to leaf
cases, classifies every key token into an AXIS, and prints per tag: PYn, WLn, PY_axes, WL_axes, DIFF.
It PRINTS computed objects and asserts nothing (rule 2). ⛔ It does NOT decide which engine is correct.

Spot-checks the census was validated against (rule 13, hand-verified before trusting the sweep):
- CMD: `grep -m1 '^WL_S11CA_VIRTUAL_WORK_SHAPE_DERIV' "$WL" | grep -oE '"[^"]*VIRTUAL_DOF[^"]*"' | sort -u | wc -l` ⇒ **16** (the 16 branch×density×DOF×VIRTUAL_DOF keys, incl. 8 off-diagonal x≠y — listed literally in the turn log).
- CMD: `grep -m1 '^PY_S11CA_VIRTUAL_WORK_SHAPE_DERIV' "$PY" | grep -oE "Str\('(LAB_HELD|MATERIAL_ADVECTED)'\)" | wc -l` ⇒ **8**.
- CMD (KINEMATIC/RELATIVE_FLUX WL keys): `grep -m1 '^WL_S11CA_KINEMATIC_BALANCE' "$WL" | grep -oE '"[A-Z_]+\|[A-Z0-9_|]+"' | sort -u` ⇒ 16 keys `branch|density|FACE_±|DOF_*` (density present). PY key is `Tuple(Str(branch),Integer(±1),Str(dof))` (no density).
- ⚠ A FRAGILE-regex artifact caught in-session: `grep -oE '"[A-Z_|]*VIRTUAL_DOF_[A-Z_]*"'` returned 8 (the `4` in `RHO4` breaks the class) — the robust `"[^"]*VIRTUAL_DOF[^"]*"` returns 16. The census uses full pipe-split, not this pattern.

## Producer citations (grep into the two engine sources)
CMD: `grep -nE 'representatives = DENSITY_REPS' scripts/S11c_a_interface_geometry_sympy_audit.py`
⇒ `:916` and `:1466`: `... if quantity in {"TRACTION","VIRTUAL_WORK_SHAPE_DERIV","CLOSURE_SHAPE_DERIV"} else ("RHO4_CONSTANT",)`.
CMD: `grep -nE 'VIRTUAL_WORK|virtual' scripts/…sympy_audit.py` ⇒ `:854` def, `:919-924`
`cases[(branch, dof, representative)] = virtual_work_cases(...)` (diagonal). Projection density path
`:1037` `rho_shape = finalize(shape(rho4_bg_alpha + parameter*rho4_bulk_1, parameter))`, `rho0` `:1039`
(fold: corrected from :1044).
CMD: `grep -nE 'RHO4_CONSTANT|virtualDof|VIRTUAL_DOF' mathematica/…audit.wl` ⇒ density profiles
`:108,:846-874`; VIRTUAL_WORK loop `:1051` `{physicalDof, dofNames}, {virtualDof, dofNames}`, key
`|VIRTUAL_DOF_` `:1039`. Flux/kinematic emit (fold: separate grep) CMD
`grep -nE 'KINEMATIC_BALANCE|RELATIVE_FLUX' mathematica/…audit.wl` ⇒ RELATIVE_FLUX payload/emit `:991`,
`:1001`; KINEMATIC emit `:1024`.

## Bespoke-tag raw structures (Family G — from raw heads, flattener returned BESPOKE)
CMD: `grep -m1 '^PY_S11CA_ADMISSIBILITY_PREMISE:' "$PY" | head -c 500`
⇒ `Tuple(Str('SUPPORT_STABILISED_BACKGROUND'), Tuple(Symbol('f_hold_0'),Symbol('t_hold_plus_0'),Symbol('t_hold_minus_0')), Str('STATIONARITY_NOT_TESTED_IN_S11C_A'), Str('S11CB_ADMISSIBILITY_PACKAGE_RESERVED'))`
CMD: `grep -m1 '^PY_S11CA_REP_INVARIANCE_RESIDUAL:' "$PY" | head -c 700` ⇒ top-level keyed by bare
`Str('RELATIVE_FLUX')`,`Str('TRACTION')`,… each nesting a `(branch,face,dof)` grid (⇒ PY nests quantity;
WL flattens quantity into the pipe-string — the 6-vs-80 raw top-level count was a nesting artifact, and
the flattened leaf counts are 64 vs 80, DIFF = VIRTUAL_DOF).

## PY leaf-field vocabulary (Family F evidence)
CMD: `grep -oE "Str\('[A-Z][A-Z0-9_]{2,}'\)" "$PY" | sed "s/Str('//;s/')//" | sort | uniq -c | sort -rn | head`
⇒ VALUE 656, MULTIGRADE 656, DIMENSION_L_T_M 656 (the flat PY field triple), plus status/operand tokens.
WL field tokens (from `axis_census.keys.txt` + raw): EXPRESSION, EXACT_SOURCE, EXACT_TRUE_AREA_SOURCE,
MULTIGRADE_EPSILON_ETA_SIGMAW, DIMENSION_L_T_M, SHAPE_DERIVATIVE.

## Provenance / discipline
User chose "matrix first, then decide full-vs-scope-down" (2026-08-24, post-compact). This matrix
CLASSIFIES divergences and NAMES the step-1 tests; it does not run them and asserts no cross-engine
result (rule 5). Counts/axis-sets are engine facts (computed), not verdicts.

## ROUND-1 FOLD measurements (both legs, orchestrator-verified against the streams)
Leg census artifacts: Grok `~/.s11_build/s11ca_t7_independent_census.{py,stdout}`; Codex
`~/.s11_build/s11ca_t7_codex_independent_census.{py,stdout}` — both reproduce the primary counts + join.

- **Family H MISSING-COVERAGE (WL omits 5 quantities).** CMD:
  `grep -m1 '^PY_S11CA_CONTROL_FORM_BASE_OPERAND:' "$PY" | grep -oE "Str\('(<18-quantity vocab>)'\)" | sort -u`
  ⇒ 18; same for WL `"<q>|"` ⇒ 13; `comm` ⇒ WL-only = ∅, PY-only =
  {EVOLUTION_TERM_ORIGINS, PROJECTION_STATIC_OPERAND, PROJECTION_DYNAMIC_OPERAND, PROJECTION_RESIDUAL,
  PROJECTION_TERM_ORIGINS}. WL UNIFORM_LIMIT_S11CA_OPERAND quantity histogram lacks the same 5.
- **Family C parentage (REP 64→80 = B+C).** CMD:
  `grep -m1 '^WL_S11CA_REP_INVARIANCE_RESIDUAL:' "$WL" | grep -oE '"[A-Z_]+\|[A-Z0-9_|]+"' | sed 's/|.*//' | sort | uniq -c`
  ⇒ CLOSURE 16, EVOLUTION_MASS_BALANCE 8, RELATIVE_FLUX 16, TRACTION 16, VIRTUAL_CONSTRAINT 8,
  VIRTUAL_WORK 16 (=80). RELATIVE_FLUX keys carry RHO4+RHOBR (density); VIRTUAL_WORK keys carry
  VIRTUAL_DOF_*. PY side: RELATIVE_FLUX 8 (no density), VIRTUAL_WORK 8 (diagonal) ⇒ the two +8 gaps.
- **DIRECTION shared.** CMD: `grep -m1 '^WL_S11CA_CONTROL_FORM_BASE_OPERAND:' "$WL" | grep -oE 'DIRECTION[A-Z0-9_]*' | sort -u`
  ⇒ DIRECTION_1/2/3. PY encodes direction as terminal integer 1/2/3 (census mislabels 2/3 as OTHER).
- **EVOLUTION origins 4↔3.** CMD PY: `grep -m1 '^PY_S11CA_EVOLUTION_TERM_ORIGINS:' "$PY" | grep -oE "Str\('(DENSITY_TIME|VELOCITY_DILATATION|BACKGROUND_ADVECTION|TRUE_AREA_FACE_FLUX)'\)" | sort -u`
  ⇒ 4; WL `grep -oE 'ORIGIN_[A-Z_]+' | sort -u` ⇒ 3 (TIME_DERIVATIVE, INPLANE_TRANSPORT, TRUE_AREA_FACE_FLUX).
- **BACKGROUND_STATE content divergence.** CMD:
  `grep -m1 '^WL_S11CA_BACKGROUND_STATE:' "$WL" | grep -oE '"[A-Z_]+"' | sort -u`
  ⇒ incl. BOUNDARY_LOADS, AFFINITY_ZERO, FACE_FLUX_ZERO, FACE_VELOCITY_ZERO, THETA_ZERO. PY BACKGROUND_STATE
  `grep -oiE "(boundary|hold|load)[a-z_]*"` ⇒ (none).
- **FACE_MAP not bespoke.** CMD: `grep -m1 '^PY_S11CA_FACE_MAP_LAB_HELD:' "$PY" | head -c 300`
  ⇒ `Tuple(Tuple(Integer(1), <map>), Tuple(Integer(-1), <map>))` (2 leaves, integer face). WL
  `<|"PLUS" -> …, "MINUS" -> …|>` (2 leaves). Both 2/FACE.
- **DIMENSIONS=38 (rule-2 gap Codex flagged).** 38 = TOP-LEVEL entry count (census-v1 measure): split the
  outer `Tuple(...)` on depth-0 commas ⇒ 38 args (script `s11ca_t7_census/s11ca_case_census.py` reports
  `PYn=38` for DIMENSIONS). ⚠ NOT the distinct-token count (`grep -oE "Str('[A-Z0-9_]+')" | sort -u | wc -l`
  ⇒ 42; different measure). WL DIMENSIONS = 4 field-keyed entries. Whole-run object; census-v2 returns BESPOKE.
- Producer lines (fold): PY form-control `:1406`; WL form-control `:1379`; WL uniform `:1537`; WL
  BACKGROUND_STATE `:865`; PY BACKGROUND_STATE `:1319`.

## Fold discipline
Both legs found NO rule-5 leak. One prior-result line in the pre-leg draft was stripped before launch (a
step-1 residual outcome must not be pre-stated). All fold claims above are engine facts (counts/tokens),
not verdicts;
which engine is correct is deferred to step 1.
