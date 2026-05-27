---
unit_id: 066
batch: III.3
verifier_model: claude-opus-4-7[1m]
verify_date: 2026-05-26T00:00:00Z
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 1
findings_total: 1
material_change: false
---

# Verification — unit 066

## Per-finding outcomes

### F1 — insufficient_verification

**Classification:** resolved

**What changed:**
The `MONOTONICITY` block in `scripts/moving_throat_pde_stage066_wall_figure_of_merit_sympy_audit.py` (lines 66-92) was replaced verbatim against the directive. Two missing derivative computations and prints were added (`dW_dJ1`, `dW_dTX`); all six derivatives are bound to named locals (`dW_dV0sq, dW_da, dW_dL, dW_dell, dW_dJ1, dW_dTX`); and six signed assertions of the form `assert sp.simplify(dW_dXxx > 0) is sp.true, "..."` (or `< 0` for `ell` and `T_X`) were appended. The comment was updated to reference notes section 3 / six signed directions. The diff (`exec_logs/stage_066_diff.patch`) shows changes are isolated to the F1 target range; no collateral edits.

**Assessment:**
The edit matches the directive exactly. The six new SymPy assertions exercise the same six signed monotonicity directions that the Mathematica engine already certifies at `mathematica/moving_throat_pde_stage066_wall_figure_of_merit_mathematica_audit.wl:73-78`, closing the asymmetric coverage gap. The asserts are non-tautological: each derivative simplifies to a distinct rational monomial in the (positive-real) base symbols, so a sign error in the underlying `W_wall` definition (e.g., a missing minus on `T_X` if it were moved to the numerator, or a sign flip on `ell`) would cause `sp.simplify(... > 0)` to no longer reduce to `sp.true` and the corresponding `assert ... is sp.true` would fire. The guard form `... is sp.true` correctly distinguishes "manifestly true under positivity assumptions" from "Relational left un-evaluated" or "false," so a future regression that yielded a non-monomial residual would also fail. Engine outputs for the six derivatives agree symbol-for-symbol with the Mathematica log (modulo naming: `J1 vs j1`, `T_X vs tx`, `L vs len`, `ell vs ell`, `V0 vs v0`).

## Exec log assessment

**SymPy:** exit=0. Notable lines:

- `dW/dJ1 = 4*pi*L**2*V0**2*a**2/(T_X*ell)` (new; matches directive's predicted output)
- `dW/dT_X = -4*pi*J1*L**2*V0**2*a**2/(T_X**2*ell)` (new; matches directive's predicted output)
- The four pre-existing derivative print lines (`dW/d(V0^2)`, `dW/da`, `dW/dL`, `dW/dell`) are unchanged in form and value.
- All six asserts pass silently (no `AssertionError` traceback); script reaches `STAGE 49 AUDIT PASSED` banner.

**Mathematica:** exit=0. Notable lines:

- `PASS: dW/d(V0^2) > 0`, `PASS: dW/da > 0`, `PASS: dW/dL > 0`, `PASS: dW/dell < 0`, `PASS: dW/dJ1 > 0`, `PASS: dW/dT_X < 0`. Mathematica engine was not modified for this finding and continues to certify the same six directions independently.

**Output freshness:** Mathematica saved output (`mathematica/output/.../...txt`, mtime 1778525747) is newer than its `.wl` (mtime 1778522213) — fresh. The SymPy script (`scripts/.../...py`, mtime 1779821523) is newer than its saved `.txt` (mtime 1778525056) — the saved `.txt` predates the F1 edit and still lacks the two new derivative prints. The orchestrator's exec log (`exec_logs/stage_066_sympy.log`, dated `2026-05-26T13:12:52-06:00`) is the post-fix authoritative run and contains the expected new lines. Per verify prompt, the exec log is the source of truth and shows exit 0 with the directive-predicted outputs. (See side observation 1.)

## Material-change assessment

`material_change`: false.

The edit adds positivity/negativity assertions on derivatives of `W_wall` but does not change `W_wall`, `W_fail`, `W_suff`, `V0_fail^2`, `V0_suff^2`, `W_H`, or any derived numeric/symbolic quantity. No constants, definitions, or substitutions changed. Downstream units that consume `W_wall` or the threshold inversions see identical results.

## Side observations (non-blocking)

1. The saved SymPy output file `scripts/output/moving_throat_pde_stage066_wall_figure_of_merit_sympy_audit.txt` is stale relative to the post-fix script (its header still reads `Date: 2026-05-11T12:44:16-06:00` and the file is missing the two new `dW/dJ1` and `dW/dT_X` print lines and the new monotonicity block layout). The post-fix exec log captured by the orchestrator (`exec_logs/stage_066_sympy.log`, 2026-05-26T13:12:52-06:00) is fresh and is what this verification relies on per the prompt's "use the exec logs the orchestrator captured" rule. The orchestrator likely refreshes saved outputs as part of a separate sync step.
2. The docstring header at line 3 still reads "Stage 49 SymPy audit" (and `banner("STAGE 49 ...")` at line 45) — naming-drift from before the renumbering to stage 066. The directive explicitly scopes this out ("Do not change the docstring header (the historical 'Stage 49' label is a separate naming-drift issue outside this finding's scope)"), so it is correctly left alone here.
3. The constant-compressibility section (lines 94-99) also lacks monotonicity assertions for `W_H` (`dW_H/dI_f`, `dW_H/dH_w`, etc.), but the original auditor did not raise this and notes §4 may not require it; not in scope for F1.

## Verdict justification

The single finding F1 was applied verbatim against the directive: the SymPy MONOTONICITY block now computes and asserts all six signed monotonicity directions matching the Mathematica engine's `expectTrue` coverage at `.wl:73-78`. The exec log shows exit 0 with the two new derivative prints in the predicted canonical form, the diff confirms changes are confined to the targeted range with no collateral edits, and both engines now independently certify the same six monotonicity signs, honoring the `\StatusExactClosure{}` policy. The new assertions are non-tautological: each tests the sign of a manifestly-signed rational monomial that would fail on any sign regression in `W_wall`. No regressions, no material change to derived quantities, no new findings raised.
