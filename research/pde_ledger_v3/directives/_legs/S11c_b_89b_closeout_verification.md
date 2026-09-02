# Closeout verification — S11c-b #89b records vs the committed reality

You are verifying that the ORCHESTRATOR-WRITTEN closeout documents accurately and honestly describe what was
actually done and committed — no overclaim, no inaccuracy, no missing caveat. This is a correctness/faithfulness
check of the RECORDS against the code and git, not a re-review of the physics (that was leg-gated separately).

## What to verify (cite file:line / git evidence for every claim you make)
1. **The committed engine matches the records.** Commit `a1be8d8f` (`git show --stat a1be8d8f`; `git log --oneline -5`).
   Confirm the engine `research/pde_ledger_v3/mathematica/S11c_b_brane_operator_mathematica_audit.wl` at HEAD:
   - emits `finalBackgroundReduction[activatedOperator]` / `activatedOrigins` (the re-freeze fix), not the
     un-activated `operatorLive`/`originsLive`;
   - keeps `KERNEL_SOURCE_OPERATOR`/`KERNEL_SOURCE_ORIGINS` un-reduced-live and the `LIVE_DIVERGENCE_FORM_OPERAND`
     tower slot un-reduced;
   - uses `expression_Times`/`expression_Plus` in `primitiveExpressionDimension` (§5.E fix);
   - has the `VALIDATED_ON_REPRESENTATIVE_BRANCH` marker on the MATERIAL_ADVECTED independence package;
   - has the `S11CB_SKIP_HEAVY_CONTROLS` deferral gate.
   Flag any claim in the records that the code does not support.
2. **The records are accurate and not overclaimed.** Read:
   - `research/pde_ledger_v3/directives/_measurements/S11c_b_89b_wl_operator_build_review.md` (the findings +
     repair + re-review record),
   - `research/pde_ledger_v3/directives/_measurements/S11c_b_89b_py_sibling_freeze_check.md` (the PY check),
   - the `## S11c-b #89b …` section at the top of `STATUS.md`,
   - `research/pde_ledger_v3/DEFERRED_HEAVY_RUNS.md` (WL entry).
   Check: does each record's description match the engine + commit? Is the "deferred to a ≥64 GB box" framing
   honest (the committed `.out` is unchanged — confirm `git show a1be8d8f -- …out` shows no `.out` change)? Is the
   re-freeze mechanism described correctly? Is the PY-check's "PY is activate-then-reduce" claim consistent with
   the PY engine source (`research/pde_ledger_v3/scripts/S11c_b_brane_operator_sympy_audit.py` — the `total_derivative`
   / `retained_grade` order at the lines it cites)? Is the strong-row jet-DEPTH reconciliation FLAG (WL 3rd-order
   vs PY `STRONG_ROW_JET_DEPTH=2`) real and correctly characterized as a spec question, not a freeze?
3. **Nothing important is missing or misfiled.** Are the "NEXT" items (PY-check done, the depth-convention flag,
   integration, #90) consistent across STATUS and the measurements? Is any owed item dropped? Did the commit
   accidentally include unrelated work (e.g. the other session's `memory/*`)? (`git show --stat a1be8d8f`.)

## Output
A numbered list of any INACCURACY, OVERCLAIM, or MISSING CAVEAT you find — `severity — file:line — what the record
says — what the code/git actually shows`. Severity ∈ {MUST-FIX, SHOULD-FIX, NIT}. Then `VERDICT: N issues (M
must-fix)` or `VERDICT: RECORDS FAITHFUL`. Do not re-litigate the physics; only whether the records tell the truth
about the committed artifact.
