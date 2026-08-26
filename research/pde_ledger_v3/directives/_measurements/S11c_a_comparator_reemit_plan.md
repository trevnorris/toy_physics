# Measurements twin — S11c_a_comparator_reemit_plan.md

Inputs (paths are literal, runnable): `WL=research/pde_ledger_v3/mathematica/out/S11c_a_interface_geometry_mathematica_audit.out`
(committed `ddecdbc2` FOR THE STEP-3 measurements below; ⚠ that same path is now the WL fix `6fae82b8` — the Step 4/5/5b sections use the CURRENT committed tree); `PY=/home/trevnorris/.s11_build/S11c_a_sympy_engine.out` (fresh run of `9b6438fa`; Step 4/5 use the fixed PY `c36beac4`).

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

## Step 5 (2026-08-26) — measure-first coverage + the two structural adjudications
Instrument `~/.s11_build/S11c_a_cov_all.py` (grows the reconcile_fixed declared map; committed
reconcile_fixed.py/reconcile_traction_check.py UNCHANGED). rule 2: prints residuals, asserts nothing.
CMD (regression + EVOLUTION face-detect): `python3 ~/.s11_build/S11c_a_run_evolution.py`
⇒ `EVOLUTION_MASS_BALANCE  8  8  0  0` (face-detect: reads `{±1/2*W0}` eval-point arg → face-tagged symbol).
CMD (algebraic-bare primaries): `python3 ~/.s11_build/S11c_a_cov_primaries.py`
⇒ FACE_VELOCITY/NORMAL/MEASURE 8/8→0, RELATIVE_FLUX 8/8→0, TRACTION 16/16→0, VIRTUAL_CONSTRAINT 8/8→0.
⚠ CLOSURE is NOT reconciled by `cov_primaries.py` — that script lacks the +4 params (lambdaAZero/lambdaVZero/
tauA/tauV) and reports CLOSURE 16 NONZERO (`Lambda_A_0`↔`lambdaAZero` unmapped). CLOSURE 16/16→0 is reproduced
by the fuller map: CMD `python3 ~/.s11_build/S11c_a_cov_all.py` ⇒ `CLOSURE_SHAPE_DERIV  16  16  0  0`.
CMD (KINEMATIC operands): `python3 ~/.s11_build/S11c_a_run_kinematic.py` ⇒ OPERAND_A 8/8→0, OPERAND_B 8/8→0.
CMD (VIRTUAL_WORK, position-aware key): `python3 ~/.s11_build/S11c_a_run_vw.py` ⇒ 16/16→0.
CONORMAL (T-a′) VERDICT A: 2 from-spec CAS legs, each residual 0 — Agent `~/.s11_build/S11c_a_conormal_leg_agent.out`
(STEP 4b/4c/5/6 all 0), Grok `research/pde_ledger_v3/_measurements/s11c_a_conormal_adjudication_grok/*.stdout`.
Rule-13 check: WL `S11c_a_interface_geometry_mathematica_audit.wl:601-610` conormalSource = conormalBackground
+ waveOrder·conormalPerturbation ⇒ conormalPerturbation is the PROBE wave, not δn̂.
PROJECTION window = NOT a finding (rule-13): `python3 ~/.s11_build/S11c_a_verify_proj_window.py` ⇒ WL window
`windowFunction[normalCoordinate - W0/2, -normalCoordinate - W0/2]` = 2 args (identical to PY O_window plus/minus);
`python3 ~/.s11_build/S11c_a_verify_proj_derivs.py` ⇒ WL carries Derivative[{m,n}][windowFunction] with BOTH
{1,0} and {0,1} (both faces). PROJECTION integrand residual OPEN — only the mechanical PERTURBATION current map + IBP bridge remain
(background-current RESOLVED, see Step 5b + §16d); see `~/.s11_build/S11c_a_T7_SCOUT_FINDINGS.md` §§13-16d;
instrument `~/.s11_build/S11c_a_run_projection.py`.

(Doc-hygiene 2026-08-26: SCOUT_FINDINGS §§11-12 are historical/superseded by §§13-14; the roadmap Step-5 PLAN bullet de-staled — EVOLUTION face-detect + CONORMAL are DONE. Commands above unchanged.)

## Step 5b (2026-08-26) — BACKGROUND-CURRENT finding RESOLVED (WL fixed 6fae82b8)
The §15-issue-3 "possible background-current finding" was measured, adjudicated, and fixed. The WL `.out` under
`mathematica/out/` was REGENERATED at `6fae82b8` (the commands below run against the current committed tree).
SURVIVAL (historical, PRE-FIX): `S11c_a_bgcurrent_check.py` reads the committed WL `.out`; at the PRE-FIX .out it
reported PROJECTION_SHAPE_DERIV / DYNAMIC_OPERAND / RESIDUAL each 8/8 cases, 1660 bg-hits, STATIC 0. ⚠ NOT
reproducible against the CURRENT tree — the committed `.out` is now the FIX (`6fae82b8`), so the same command now
reports 0 bg-hits for all four tags (that zero IS the post-fix confirmation below).
CMD (both physics consults COMPUTED the survivor = 0 under continuity — the reason zeroing loses no physics):
Grok `~/.s11_build/S11c_a_bgcurrent_consult_grok.txt` (`FINAL_RESIDUAL_AFTER_CONTINUITY_IBP_DECAY: 0`);
Codex `~/.s11_build/S11c_a_bgcurrent_consult_codex.log` (`RESIDUAL_AFTER_BOUNDARY_DROP = 0`; total-w-derivative).
Consult scripts (continuity→0): Grok `/tmp/s11ca_tf_bgcurrent_q1.py` (`FINAL_RESIDUAL_AFTER_CONTINUITY_IBP_DECAY:
0`); Codex computed it in its consult (`~/.s11_build/S11c_a_bgcurrent_consult_codex.log`, `RESIDUAL_AFTER_BOUNDARY
_DROP = 0`, total-w-derivative). ⛔ NOT `/tmp/s11ca_bgcurrent_clean.py` — that was the earlier Agent ADJUDICATION
(showed the survivor EXISTS, not the continuity cancellation). This continuity result is the STEP-RECORD
justification; it was deliberately WITHHELD from the WL-fix build directive (a builder leak).
CMD (post-fix committed WL `.out` — the fix worked, no regression):
`grep -c '^WL_S11CA_' research/pde_ledger_v3/mathematica/out/S11c_a_interface_geometry_mathematica_audit.out` ⇒ 40;
`grep -c 'currentWBackground\|currentXBackground' …audit.out` ⇒ 0 (baseline 667);
`grep -c 'currentWPerturbation\|currentXPerturbation' …audit.out` ⇒ 656 (perturbation current intact).
CMD (2 build legs, velocity-probe FORM ablation — bg current TRACKS a nonzero probe ⇒ genuine ρ·v, not :=0):
Grok `research/pde_ledger_v3/mathematica/_measurements/s11c_a_wl_bgcurrent_review_grok/ablate_velocity_probe.stdout`
(`currentWZero = vpW*rhoBulkBackground[...]`, `TRACKS_probe_* = True`), reverse zero-velocity ⇒ 0, `Part_pkspec1
_message_count = 0`. Directive `S11c_a_wl_bgcurrent_fix_directive.md` (+twin), fix committed `6fae82b8`.
NET: two confirmed T7 findings, BOTH FIXED — shifted-trace (PY `c36beac4`) + free-premise bg-current (WL
`6fae82b8`). REMAINING coverage (mechanical): projection integrand perturbation-current map + IBP bridge →
sweep → Codex comparator. Detail `~/.s11_build/S11c_a_T7_SCOUT_FINDINGS.md` §§16–16d.
