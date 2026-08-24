# Measurements twin — S11c_a_comparator_reemit_plan.md

The plan's divergence map was read off the two real tag streams and the two directive-review leg logs.
Inputs: WL `research/pde_ledger_v3/mathematica/out/S11c_a_interface_geometry_mathematica_audit.out`
(committed `ddecdbc2`); PY `~/.s11_build/S11c_a_sympy_engine.out` (fresh run of `9b6438fa`). Many rows here
duplicate `_measurements/S11c_a_comparator_build_directive.md` (same greps); the leg-derived rows carry the
leg log as their command.

## Name-join (39 exact)
CMD: `grep -oE '^(PY|WL)_S11CA_[A-Z0-9_]+' <out> | sort -u | comm` ⇒ 39 non-local stems each side, exact 1-1.

## Case-key encoding
CMD: `grep -m1 '^WL_S11CA_FACE_NORMAL' WL` ⇒ `<|"LAB_HELD|FACE_PLUS|DOF_DELTA_W" -> …` (pipe-string).
CMD: `grep -m1 '^PY_S11CA_FACE_NORMAL' PY` ⇒ `Tuple(Tuple(Tuple(Str('LAB_HELD'), Integer(1), Str('DELTA_W')), …` (positional tuple; face = Integer 1/-1).
CMD (leg): Grok directive-review `~/.s11_build/S11c_a_comparator_directive_grok.log` ⇒ "1190 unique WL outer keys tokenize to **49 unclassified tokens**" outside {branch,face,dof,density}: DIRECTION_1/2/3, FIELD_PRESSURE/BULK_DENSITY/NORMAL_CURRENT/CURRENT_X1…, VIRTUAL_DOF_DELTA_W/ZETA_C, ORIGIN_*, quantity names. FACE_MAP keys = `PLUS`/`MINUS` (not FACE_PLUS), sympy `faceLabel` WL src 547-548; PY FACE_MAP keys = Integer 1/-1. Form-control PY key `("FACE_NORMAL","LAB_HELD",1,"DELTA_W",1)` — integer 1 is face AND direction (Codex log `~/.s11_build/S11c_a_comparator_directive_codex.log`, PY src 1432-1444).

## Field-envelope + names
CMD: `for t in VALUE EXPRESSION EXACT_SOURCE MULTIGRADE MULTIGRADE_EPSILON_ETA_SIGMAW DIMENSION_L_T_M OPERAND_A OPERAND_B; do grep -c "'$t'" PY; grep -c "\"$t\"" WL; done` ⇒ PY {VALUE 21, MULTIGRADE 21, DIM 21, OPERAND_A 9, OPERAND_B 9}; WL {EXPRESSION 37, EXACT_SOURCE 4, MULTIGRADE_EPSILON_ETA_SIGMAW 40, DIM 40}. ⇒ VALUE↔EXPRESSION, MULTIGRADE↔MULTIGRADE_EPSILON_ETA_SIGMAW; DIM identical.
CMD (leg): WL nests record under `SHAPE_DERIVATIVE` beside `EXACT_SOURCE` (FACE_NORMAL WL out line 6; WL src 946); KINEMATIC WL siblings `OPERAND_A_SHAPE_DERIVATIVE`/`OPERAND_B_SHAPE_DERIVATIVE`/`RESIDUAL_A_MINUS_B` (WL out 11); PY OPERAND_A/B nested INSIDE VALUE (PY src 837-841). Grok+Codex logs.

## CAS representation
CMD: `for t in 'Dummy(' 'Integral(' 'Subs(' O_window windowFunction; do grep -c "$t" PY; done` ⇒ each 13 (11 non-local proj + 2 local); PY windowFunction=0.
CMD: `for t in 'Inactive[Integrate]' windowFunction 'Derivative[' O_window; do grep -c -F "$t" WL; done` ⇒ Inactive[Integrate]=10 lines(1120 apps per leg), windowFunction=10, Derivative[=29, O_window=0.
CMD (leg, Grok): windowFunction arity histogram = **{2: 660}** (WL definition `windowFunction[gplus,gminus]`, WL src 749-756) — window is 2-arg BOTH engines; my earlier "1-arg WL" claim was FALSE. Inactive[Equal]=1233, Inactive[Greater]=8; `parse_mathematica("Equal[0,0]")` → native True. Nested `Derivative[o…,{a,b}]` = 22656; S11b WL_DERIVATIVE regex misses them. Integer Assoc keys `<|1->,-1->|>` in BACKGROUND_STATE = 18 (S11b `_parse_association_key` rejects → UNCOMPARED).
CMD: re-run SymPy ⇒ `git diff S11c_a_exports.py` = only `Dummy(dummy_index=NNNN)` counters ⇒ Dummy canon by NAME, keep distinct.

## Decision provenance
User chose the shared-schema re-emit over the partial comparator (2026-08-24), after I gave the pay-later
(hand-adjudication of ~half the objects + b–e recurrence) vs re-emit-cost (§7 amendment + emit patches +
re-run + light re-review + trivial comparator, ~2-3 cycles; physics-preservation mechanically checkable).
