# Measurements twin — S11c_a_comparator_reemit_plan.md

Inputs (paths are literal, runnable): `WL=research/pde_ledger_v3/mathematica/out/S11c_a_interface_geometry_mathematica_audit.out`
(committed `ddecdbc2`); `PY=/home/trevnorris/.s11_build/S11c_a_sympy_engine.out` (fresh run of `9b6438fa`).

⚠ Fix (2026-08-24 audit): the two commands below were inaccurate as first written — corrected here.
Run from repo root; `PY`/`WL` as defined above.

## Name-join (39 exact)
CMD: `comm -3 <(grep -oE '^PY_S11CA_[A-Z][A-Z0-9_]*' "$PY" | grep -v _LOCAL_ | sed 's/^PY_//' | sort -u) <(grep -oE '^WL_S11CA_[A-Z][A-Z0-9_]*' "$WL" | grep -v _LOCAL_ | sed 's/^WL_//' | sort -u)`
⇒ empty (39 non-local stems each side, exact 1-to-1). (Without the `_LOCAL_` filter the raw `comm -3`
shows 9 local-only stems: PY 8 + WL 1 — excluded here.)

## CASE-STRUCTURE DISAGREEMENTS (Codex plan-review, orchestrator-verified — the headline)
CMD (robust regex — the earlier `"[A-Z_|]*VIRTUAL_DOF_[A-Z_]*"` returned 8 because `RHO4` contains a
digit): `grep -m1 '^WL_S11CA_VIRTUAL_WORK_SHAPE_DERIV' "$WL" | grep -oE '"[^"]*VIRTUAL_DOF[^"]*"' | sort -u | wc -l`
⇒ **16** WL cases (branch×density×DOF×VIRTUAL_DOF), incl. 8 off-diagonal `DOF_x|VIRTUAL_DOF_y` (x≠y).
CMD: `grep -m1 '^PY_S11CA_VIRTUAL_WORK_SHAPE_DERIV' "$PY" | grep -oE "Str\('(LAB_HELD|MATERIAL_ADVECTED)'\)" | wc -l`
⇒ **8** PY top-level cases (branch×dof×density; virtual DOF TIED to physical DOF, no off-diagonal). ⇒ WL has
16 the 8-case PY does not — a genuine case-set divergence, not serialization.
CMD: `grep -m1 '^PY_S11CA_PROJECTION_SHAPE_DERIV' "$PY" | grep -oE 'RHO4_CONSTANT|RHOBR_CONSTANT' | sort -u`
⇒ `RHO4_CONSTANT`,`RHOBR_CONSTANT` (PY projection IS density-decomposed).
CMD: `grep -m1 '^WL_S11CA_PROJECTION_SHAPE_DERIV' "$WL" | grep -oE '"LAB_HELD\|DOF_[A-Z_]*"' | head`
⇒ `"LAB_HELD|DOF_DELTA_W"` … (WL projection keys branch|DOF only — NO density axis). ⇒ Codex further ran the
two PY density cases and got residual `False` (they differ) — not a droppable duplicate. Spec must decide.
Further case-set/semantic rows (CONTROL_FORM/UNIFORM coverage; PY `(background,ε·deriv,total)` tuple vs WL
coefficient + waveOrder stripped into MULTIGRADE; FACE_SHIFT field explosion) cited in the plan with
sympy/WL src line numbers; verify each in the step-0 feasibility matrix.

## Field-name vocabulary
CMD: `for t in VALUE EXPRESSION EXACT_SOURCE MULTIGRADE MULTIGRADE_EPSILON_ETA_SIGMAW DIMENSION_L_T_M; do printf '%s PY=%s WL=%s\n' "$t" "$(grep -c "'$t'" "$PY")" "$(grep -c "\"$t\"" "$WL")"; done`
⇒ PY{VALUE 21,MULTIGRADE 21,DIM 21}; WL{EXPRESSION 37,EXACT_SOURCE 4,MULTIGRADE_EPSILON_ETA_SIGMAW 40,DIM 40}.

## CAS representation
CMD: `for t in 'Inactive[Integrate]' windowFunction 'Derivative[' O_window; do printf '%s WL=%s\n' "$t" "$(grep -c -F "$t" "$WL")"; done`
⇒ Inactive[Integrate]=10 lines, windowFunction=10, Derivative[=29, O_window=0.
CMD: `for t in 'Dummy(' 'Integral(' O_window windowFunction; do printf '%s PY=%s\n' "$t" "$(grep -c "$t" "$PY")"; done`
⇒ Dummy/Integral/O_window each 13; windowFunction 0.
CMD (leg, Grok `~/.s11_build/S11c_a_comparator_directive_grok.log`): WL window arity histogram `{2:660}` (2-arg
both engines); 49 unclassified case-key tokens; nested `Derivative[o,{a,b}]`=22656; integer Assoc keys
`<|1->,-1->|>` in BACKGROUND_STATE=18; `parse_mathematica("Equal[0,0]")`→native True.
CMD: re-run SymPy ⇒ `git diff S11c_a_exports.py` = only `Dummy(dummy_index=NNNN)` counters (Dummy = per-run
noise; canon by NAME). ⚠ the emitted payload also feeds `export_candidates` (sympy 333-343 → 1852-1878) ⇒ any
PY emit reformat mutates exports.

## Provenance
User chose re-emit (2026-08-24) over the partial comparator. Codex plan-review then found the emit-only
assumption false (case-structure + semantic divergences); plan revised to a feasibility/adjudication-matrix
first. No expected cross-engine RESULT is asserted; the case-structure divergences are engine facts (counts),
not comparator verdicts.

## Progress (2026-08-24)
Step 0 committed `3c7f9137` (matrix + 2 legs + fold); Step 1 committed `3491a376` (verdicts + 2 legs +
fold). Verifiable: `git log --oneline | grep -E '3c7f9137|3491a376'`. User chose FULL RECONCILE; NEXT =
step-2 engine patches (both engines). Provenance carried in the two committed docs + their twins.

## Step 4 (2026-08-26) — fixed-PY vs WL reconciliation confirmation (the shifted-trace fix)
Fixed PY committed `c36beac4`; WL `a7459cb8` unchanged. Reconciliation instrument =
`~/.s11_build/S11c_a_reconcile_fixed.py` (measure_reconcile over the fixed PY .out + the mechanical
`d_w_X`↔`X_dw` perturbation-jet rename). The RELATIVE_FLUX + TRACTION headline (TRACTION also folds the
mechanical `mu_theta_L→mu_theta` rename + `sp.cancel` for the λ_X complex denominators) reproduces via:
```
$ python3 /home/trevnorris/.s11_build/S11c_a_reconcile_traction_check.py
RELATIVE_FLUX            join=8 zero=8 nonzero=0
TRACTION                 join=16 zero=16 nonzero=0
```
Full per-tag record (incl. the EVOLUTION comparator single-face gap, NOT a finding) in
`_measurements/S11c_a_py_shifted_trace_fix_directive.md`.
