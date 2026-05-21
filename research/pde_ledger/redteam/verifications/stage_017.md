---
unit_id: 017
batch: I.2
verifier_model: claude-opus-4-7[1m]
verify_date: 2026-05-21T00:00:00Z
verdict: verified
sympy_exit: 0
mathematica_exit: n/a
findings_resolved: 3
findings_total: 3
material_change: false
---

# Verification — unit 017

## Per-finding outcomes

### F1 — tautological_check

**Classification:** resolved

**What changed:**

`scripts/moving_throat_pde_stage017_parent_throat_action_weak_axisym_sympy_audit.py:19-25` — the `if m == 0: return sp.Integer(1)` shortcut is removed; the `same_sign` cross-term guard is now gated by `if m != 0:`; the return expression is now the single `sp.simplify((sp.Integer(-1) ** m) * gaunt(2, 2, 2, 0, m, -m) / base)` for all m (diff lines 7-18).

**Assessment:**

Edit matches the directive's required change exactly. For m=0 the helper now computes `gaunt(2,2,2,0,0,0) / gaunt(2,2,2,0,0,0)` via SymPy rather than returning a literal `Integer(1)`, so the downstream `assert_zero('Y20 overlap lane 20', lam20 - 1)` at line 43 is no longer trivially satisfied — it now genuinely exercises SymPy's `gaunt` for the m=0 cell. The sympy exec log still shows `lambda_(20) = 1, lambda_(21) = 1/2, lambda_(22) = -1.` (line 13), confirming the m=0 branch still yields 1 through the Gaunt machinery. No collateral changes outside the named lines.

### F2 — tautological_check

**Classification:** resolved

**What changed:**

`scripts/moving_throat_pde_stage017_parent_throat_action_weak_axisym_sympy_audit.py` — the two tautological asserts at the old lines 79-80 (`wall-only K1 specialization`, `wall-only H_even specialization`) are deleted (diff lines 25-26). Six new lane-cross-checks are inserted after the `K1_gate_*` / `Hev_*` definitions at the script's current lines 87-97 (diff lines 34-44): `generic_K1 = -M1 + K1w/9`, `generic_Hev = 2 M1/3 - K1w/27`, and six `assert_zero("generic {K1,Hev} vs lane 2m", ... - eps*lam2m*generic_*)` lines.

**Assessment:**

Edit matches the directive exactly — both the deletions and the insertion at the requested location (after the `Hev_22` definition, before `# Solve the wall-only even gates.`). The new assertions are non-tautological in the right way: `K1_gate_2m` is built from `D21_2m + D01_2m/9` where `D21_2m = -eps*lam2m*M1` and `D01_2m = eps*lam2m*K1w` come from the lane-resolved Gaunt-derived `lam20, lam21, lam22`. The RHS `eps*lam2m*generic_K1` uses an independently hand-written `generic_K1 = -M1 + K1w/9`. If `generic_K1`'s coefficient (e.g., -M1 vs -M1/2) or any of the three `lam2m` values were ever wrong, the assertion would fail on at least one lane. The auditor's self-test (lines 173-177 of the original report) already pre-checked these reduce correctly under the verified (1, 1/2, -1) lane ratios. The `K1_wall`/`Hev_wall` variables remain in scope for the wall-matrix determinant and `sol_even` Solve, as the directive specified.

### F3 — missing_verification_script (missing_mathematica)

**Classification:** resolved

**What changed:**

New file `mathematica/moving_throat_pde_stage017_parent_throat_action_weak_axisym_mathematica_audit.wl` (118 lines) implements an independent Mathematica audit using direct `Integrate[Sin[theta] * SphericalHarmonicY[2,q1,theta,phi] * SphericalHarmonicY[2,q2,theta,phi] * SphericalHarmonicY[2,q3,theta,phi], ...]` over the sphere — the more-independent of the two routes the directive permitted. The file covers all 12 manifest claims M1-M12 with PASS/FAIL labels: lane ratios m=0,1,2 (M1-M3); same-sign cross-term vanishing (M4); wall Jacobian Det = 1/27 (M9); wall-only Solve unique-trivial (M10) checked via solution-count plus dK/dM values; K1/Hev lane gate cross-checks for each m (M7, M8); grouped trace + b=3a for mass, stiffness, Xi load, and prefactor (M5, M6, M11, M12). Ends with `Print["STATUS: PASS"]; Exit[0]`.

**Assessment:**

Substantive independence is correct — the script uses Mathematica's `SphericalHarmonicY`, `ComplexExpand`, and analytic `Integrate` over the sphere, not `ThreeJSymbol` and not SymPy's Wigner machinery. I grep'd the file for the forbidden strings `sympy`, `gaunt(`, `real_y20_square_ratio`, `grouped_trace_anomaly`, `assert_zero`, `assert_nonzero` — none present. Quick arithmetic on the Jacobian (line 63): rows are `{D[-dMsym + dKsym/9, dKsym], D[..., dMsym]} = {1/9, -1}` and `{D[2 dMsym/3 - dKsym/27, dKsym], D[..., dMsym]} = {-1/27, 2/3}`; Det = (1/9)(2/3) - (-1)(-1/27) = 2/27 - 1/27 = 1/27. Correct. The lane-gate cross-checks (lines 76-92) reproduce the same non-tautological structure as the SymPy F2 fix but built from `laneFactor[m]` computed via direct integration, so they constitute a genuine independent re-derivation. The output transcript at `mathematica/output/moving_throat_pde_stage017_parent_throat_action_weak_axisym_mathematica_audit.txt` shows residual = 0 and `PASS:` for all 23 labeled checks plus `STATUS: PASS`. No collateral edits outside the named files.

## Exec log assessment

**SymPy:** exit=0. Notable lines:

- `lambda_(20) = 1, lambda_(21) = 1/2, lambda_(22) = -1.` (line 13) — confirms F1's now-active Gaunt evaluation for m=0 still yields 1.
- `K1_(20,21,22)   = (eps*(K1w - 9*M1)/9, eps*(K1w - 9*M1)/18, eps*(-K1w + 9*M1)/9)` (line 21) — confirms the F2 lane gates still print correctly; the six new `generic {K1,Hev} vs lane *` asserts must have all passed (silent `assert_zero`).
- `Solving K1_wall = 0 and H_even,wall = 0 gives: [{dKsym: 0, dMsym: 0}]` (line 26) — trivial-solution branch intact.
- Final `STATUS: PASS` (line 40), `exit_code: 0` (line 42).

**Mathematica:** exit=n/a — the orchestrator did not capture `/var/projects/toy_physics/research/pde_ledger/redteam/exec_logs/stage_017_mathematica.log` for this iteration (file missing). However, the saved transcript `mathematica/output/moving_throat_pde_stage017_parent_throat_action_weak_axisym_mathematica_audit.txt` (mtime 2026-05-21 15:02, ~6 min newer than the `.wl` mtime 2026-05-21 13:26) shows all 23 named checks reporting `residual = 0` followed by `PASS:`, ending with `STATUS: PASS` on line 48. Substance is verifiable from the transcript even though the orchestrator's exec log envelope is absent. Flagging the missing exec log as a procedural gap (not a substantive blocker).

**Output freshness:** Confirmed. SymPy script mtime 1779391432, sympy output mtime 1779397233 (newer by ~97 min). Mathematica `.wl` mtime 1779391606, mathematica output mtime 1779397342 (newer by ~96 min). Both outputs were regenerated after the post-fix edits.

## Material-change assessment

`material_change`: false.

Rationale: F1 changed the m=0 computation path but the value `lambda_(20) = 1` is unchanged. F2 deleted two tautological asserts and added six non-tautological ones; no published derived quantity (K1_gate_*, Hev_*, lane ratios, det = 1/27, Mbar/Kbar/Xibar/Pbar, a/b coefficients) changed numerically — the script transcript (lines 13, 16-17, 21-22, 26, 30-32, 34-35) is identical in substance to what a downstream consumer would expect. F3 added a new Mathematica file with no downstream code dependencies. No prose tracker or paper-tex coupling within scripts-only scope.

## Side observations (non-blocking)

- The orchestrator did not write `redteam/exec_logs/stage_017_mathematica.log`, even though a Mathematica script now exists for stage 017. The saved `mathematica/output/.../audit.txt` transcript is fresh and complete, so verification is unblocked, but the run envelope (argv, date, exit_code line) that other stages capture is not present here. Worth fixing in the orchestrator so the standard "exec log + saved txt" pair is consistently produced when a `.wl` first appears.
- F2's six new asserts are inserted between the `Hev_22` definition and `# Solve the wall-only even gates.` block, exactly where the directive said to put them. Good placement — it keeps `K1_wall`, `Hev_wall` available for the subsequent Jacobian and Solve calls without reordering.
- The Mathematica script's `M10 wall-only solution count` check (`Length[wallSolve] - 1`) is a slightly stronger statement than the SymPy `if sol_even != [{dKsym: 0, dMsym: 0}]: raise` — it explicitly confirms uniqueness rather than relying on the dict-equality semantics. Good cross-engine spread.

## Verdict justification

All three findings are resolved with edits that match the directive's required-change text. F1 routes the m=0 Gaunt evaluation through SymPy's own `gaunt`, F2 replaces tautological construction-asserts with six lane-resolved cross-checks against the independently written generic formulas, and F3 adds a substantively independent Mathematica audit that re-derives all 12 manifest claims via direct `SphericalHarmonicY` integration and reports 23/23 PASS plus `STATUS: PASS`. SymPy exec log shows `exit_code: 0` and the canonical transcript content is intact. Mathematica exec log is missing as an orchestrator artifact, but the saved output transcript is fresh and complete; substance is verifiable. No regressions or collateral edits.
