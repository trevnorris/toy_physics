# Measurements twin — S11c_a_comparator_build_directive.md

Every schema claim in the directive was read off the two real tag streams:
- WL: `research/pde_ledger_v3/mathematica/out/S11c_a_interface_geometry_mathematica_audit.out` (committed `ddecdbc2`).
- PY: `~/.s11_build/S11c_a_sympy_engine.out` (fresh run of the SymPy engine at `9b6438fa`, exit 0; NOT committed — the comparator regenerates/uses it. ⚠ running the SymPy engine rewrites committed `S11c_a_exports.py` only in `Dummy dummy_index=` counters — restored via `git checkout`).

## Join namespace (recon agent af381949 + confirmed)
CMD: `grep -oE '^(PY|WL)_S11CA_[A-Z0-9_]+' <out> | sort -u` ⇒ 39 non-local stems EACH side, exact 1-to-1 match; WL 1 local (`LOCAL_TAG_NAMES`), PY 8 locals. No unpaired non-local rows at the name level.

## Inner field-name vocabulary
CMD: `for tok in VALUE EXPRESSION EXACT_SOURCE MULTIGRADE MULTIGRADE_EPSILON_ETA_SIGMAW DIMENSION_L_T_M …; do grep -c "'$tok'" PY; grep -c "\"$tok\"" WL; done` ⇒
  PY: VALUE=21, MULTIGRADE=21, DIMENSION_L_T_M=21, STATUS_TOKEN=1, OPERAND_A=9, OPERAND_B=9, RESIDUAL=9. EXPRESSION/EXACT_SOURCE/MULTIGRADE_EPSILON_ETA_SIGMAW=0.
  WL: EXPRESSION=37, EXACT_SOURCE=4, MULTIGRADE_EPSILON_ETA_SIGMAW=40, DIMENSION_L_T_M=40, RESIDUAL=1. VALUE/MULTIGRADE/OPERAND_A/B=0.
⇒ safe alias table = {VALUE→EXPRESSION, MULTIGRADE→MULTIGRADE_EPSILON_ETA_SIGMAW}; DIMENSION_L_T_M identical. WL-only EXACT_SOURCE (4 objs, exact pre-jet form) + PY-only OPERAND_A/B/RESIDUAL (KINEMATIC_BALANCE) surface as KEY differences (NOT aliased).

## Outer case-key encoding
CMD: `grep -m1 '^WL_S11CA_FACE_NORMAL' WL` ⇒ key `"LAB_HELD|FACE_PLUS|DOF_DELTA_W"` (pipe-string, FACE_PLUS/FACE_MINUS, DOF_ prefix).
CMD: `grep -m1 '^PY_S11CA_FACE_NORMAL' PY` ⇒ key `Tuple(Str('LAB_HELD'), Integer(1), Str('DELTA_W'))` (tuple, faces as Integer 1/-1, no DOF_ prefix).
Vocabularies (BRANCHES/FACES/DOFS/DENSITY_REPS, sympy src lines 44-47) disjoint ⇒ dimension-classification is unambiguous: FACE_PLUS⇔1, FACE_MINUS⇔-1, strip DOF_.

## Window + integral representation (the projection cluster — 11 non-local objects)
CMD: `for tok in 'Dummy(' 'Integral(' 'Subs(' 'O_window' 'Derivative(' windowFunction; do grep -c "$tok" PY; done` ⇒ PY each = 13 lines (11 non-local proj objects + 2 local export); PY windowFunction=0.
CMD: `for tok in 'Inactive[Integrate]' windowFunction 'Derivative[' O_window normalCoordinate; do grep -c -F "$tok" WL; done` ⇒ WL Inactive[Integrate]=10, windowFunction=10, Derivative[=29, O_window=0, normalCoordinate=10.
⇒ PY window = `O_window(G_+, G_-)` (2-arg, via Subs/Derivative of an abstract fn, Dummy bound args); WL window = `windowFunction[single composite arg]` (Inactive[Integrate] over normalCoordinate). Head rename O_window↔windowFunction is a name; ⛔ the 2-arg↔1-arg PARAMETERIZATION is NOT bridged — surfaces as a projection DISAGREE for orchestrator adjudication (directive delta #7). Objects carrying it: PROJECTION_* (5), CONTROL_FORM_* (3), UNIFORM_LIMIT_* (3).

## Dummy non-determinism
CMD: re-run SymPy engine ⇒ `git diff S11c_a_exports.py` = 5 lines, all `Dummy(..., dummy_index=NNNN)` counter changes ⇒ Dummy canonicalization by NAME (not index) is mandatory (delta #6).

## BLINDNESS / rule 5
Directive references SCHEMA (field names, key encoding, marker heads) — the join contract — NEVER a real computed VALUE and NEVER an expected cross-engine result. Fixtures are synthetic placeholders; the orchestrator runs the frozen comparator against the real pair after review.
