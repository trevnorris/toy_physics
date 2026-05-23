---
unit_id: 067
batch: III.3
verifier_model: claude-opus-4-7-1m
verify_date: 2026-05-22T00:00:00Z
verdict: verified
sympy_exit: n/a
mathematica_exit: n/a
findings_resolved: 4
findings_total: 4
material_change: false
---

# Verification — unit 067

## Per-finding outcomes

### F1 — tautological_check

**Classification:** resolved

**What changed:**
- `scripts/moving_throat_pde_stage067_sech_gaussian_sympy_audit.py:91-93` adds the three-line comment immediately above the `expect_zero("C^2(r) - C^2(pi/r) under duality", C2_dual - C2_target)` call (now at line 94), with text matching the directive verbatim.
- `mathematica/moving_throat_pde_stage067_sech_gaussian_resonance_mathematica_audit.wl:64-66` adds the three-line `(* ... *)` comment immediately above the `expectZero["C^2(r) - C^2(pi/r) under duality", c2Dual - c2Target]` call (now at line 67), with text matching the directive verbatim.

**Assessment:**
Comment-only change in both files. Wording is exactly what the directive specified (down to capitalization of "ANY function I" and reference to "section 5" / "below"). The `expect_zero`/`expectZero` logic was not modified — assertion remains structurally the same algebraic-implication check, but is now honestly labeled. No collateral edits. Output transcripts contain the same `... under duality = 0` / `PASS` line as before the edit (sympy line 20; mathematica line 17), consistent with a comment-only change.

### F2 — tautological_check

**Classification:** resolved

**What changed:**
- `scripts/moving_throat_pde_stage067_sech_gaussian_sympy_audit.py:116-117` inserts the two-line "Tautological: ..." comment immediately above the `expect_zero("self-dual overlap-slope relation", ...)` call (now at lines 118-121).
- `scripts/moving_throat_pde_stage067_sech_gaussian_sympy_audit.py:127-130` inserts the four-line "Tautological: dC2_selfdual is a hand-written derivative formula..." comment immediately above the `expect_zero("stationary derivative of C^2 at the self-dual point", ...)` call (now at lines 131-134).
- `mathematica/moving_throat_pde_stage067_sech_gaussian_resonance_mathematica_audit.wl:89-92` inserts the four-line `(* Tautological: Solve[2*C2PrimeLeft == 0] returns C2PrimeLeft -> 0; ... *)` comment immediately above the `expectZero["self-dual C^2 stationary slope from symmetry solve", ...]` call (now at lines 93-96).

**Assessment:**
All three comment blocks match the directive text verbatim. Assertion logic is unchanged — no rewriting of `duality_tangent_at_rstar`, `dC2_selfdual`, or `c2StationarySolve`. The output transcripts still show `self-dual overlap-slope relation = 0`, `stationary derivative of C^2 at the self-dual point = 0`, and `self-dual C^2 stationary slope from symmetry solve = 0` / `PASS` — these remain tautological by construction, but are now correctly labeled as such. No collateral edits.

### F3 — hardcoded_result

**Classification:** resolved

**What changed:**
- `scripts/moving_throat_pde_stage067_sech_gaussian_sympy_audit.py:68-75` inserts the norm-derivation block immediately after `print("N_(phi phi)     =", Npp)` (line 66) and before the blank line preceding section 2. The block declares a real symbol `y`, computes `Nss_integral = sp.integrate((sp.sech(y/wf)**2).rewrite(sp.cosh), (y, -sp.oo, sp.oo))` and `Npp_integral = sp.integrate(sp.exp(-2*y**2/wg**2), (y, -sp.oo, sp.oo))`, prints both, and calls `expect_zero` against `Nss = 2*wf` and `Npp = wg*sqrt(pi/2)`.

**Assessment:**
Codex applied a documented deviation: it rewrote `sech(y/wf)**2` via `.rewrite(sp.cosh)` so SymPy's integrator resolves the closed form in this environment. The deviation is appropriate — `sp.sech` directly is not always integrated in current SymPy releases, while `1/cosh(...)**2` is. The deviation only affects the expression handed to `sp.integrate`; the final value compared against `Nss = 2*wf` is the same. The output transcript confirms substantive evaluation:
- Line 11: `integrate(sech(y/w_f)^2)        = 2*w_f`
- Line 12: `integrate(exp(-2 y^2/w_g^2))    = sqrt(2)*sqrt(pi)*w_g/2`
- Line 13: `N_(sigma sigma) integral - 2 w_f = 0`
- Line 14: `N_(phi phi) integral - w_g sqrt(pi/2) = 0`

These are non-tautological — the integrals are computed independently by SymPy from the actual sech^2 and exp(-2y^2/wg^2) profiles, then compared to the declared norms. If a future edit changed the profile or factors, these checks would catch it. The Mathematica side already had the equivalent `Integrate[...]` derivation (.wl lines 45-46), so both engines now anchor the norms.

### F4 — hardcoded_result

**Classification:** resolved

**What changed:**
- `mathematica/moving_throat_pde_stage067_sech_gaussian_resonance_mathematica_audit.wl:122-125` inserts the four-line `(* c2Target / presTarget are the sympy mpmath quad results from scripts/output/... *)` comment immediately above the `c2Target = ToExpression[...]` line (now at line 126).

**Assessment:**
Comment text matches the directive verbatim. The literal `c2Target` and `presTarget` values are unchanged, as required. The subsequent `expectApprox` calls are unchanged. Output transcript still shows `C_res^2 numeric check` PASS / `P_res numeric check` PASS, agreeing to ~35 digits. The labeling now correctly identifies these as cross-engine empirical targets rather than analytic benchmarks.

## Exec log assessment

**SymPy:** exec log file not present at `/var/projects/toy_physics/research/pde_ledger/redteam/exec_logs/stage_067_sympy.log`. The orchestrator did not capture a fresh sympy log alongside the diff for this iteration. However, the persisted output file `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage067_sech_gaussian_sympy_audit.txt` has mtime 2026-05-22 19:56, which post-dates the script mtime (19:54). The output text shows all assertions producing residual `= 0`:
- Line 13: `N_(sigma sigma) integral - 2 w_f = 0`
- Line 14: `N_(phi phi) integral - w_g sqrt(pi/2) = 0`
- Line 20: `C^2(r) - C^2(pi/r) under duality = 0`
- Line 27: `self-dual overlap-slope relation = 0`
- Line 28: `stationary derivative of C^2 at the self-dual point = 0`
- Line 47: numerical duality samples all `0.0` or `~1e-61`
The presence of these lines, plus completion of the `FINAL LEDGER` banner, indicates the script ran to completion (any `expect_zero` failure raises and aborts before the ledger prints).

**Mathematica:** exec log file not present at `/var/projects/toy_physics/research/pde_ledger/redteam/exec_logs/stage_067_mathematica.log`. Persisted output `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage067_sech_gaussian_resonance_mathematica_audit.txt` has mtime 19:56 (post-dating the .wl mtime 19:53). Output shows:
- Line 7-10: norm checks `PASS`
- Line 16-17: duality implication check `PASS`
- Line 24-25: self-dual stationary slope `PASS`
- Line 34-37: numeric cross-engine checks `PASS`
- Line 42-56: all 5 duality samples `PASS`
- Line 68-69, 75-76: monotonicity scans `PASS`
- Final line: `Stage 067 Mathematica audit passed.` (only printed if `Exit[0]` is reached without `Exit[1]` from `fail`).

**Output freshness:** confirmed. Sympy script mtime 19:54 < output mtime 19:56. Mathematica .wl mtime 19:53 < output mtime 19:56. Both outputs were regenerated post-edit.

Note: The verifier's expected exec-log files are missing, so the conventional sympy/mathematica exit codes are unavailable. The surrounding evidence (post-edit output files, all expected PASS lines present, FINAL LEDGER reached, "audit passed" footer reached) is consistent with exit=0 in both engines. Marked `sympy_exit` and `mathematica_exit` as `n/a` rather than guessing.

## Material-change assessment

`material_change`: false.

Three of four edits are comment-only and do not change any computed value. The F3 edit adds two new `expect_zero` assertions and two new print lines to the sympy output, but the underlying declared norms `Nss = 2*wf` and `Npp = wg*sqrt(pi/2)` are unchanged, all downstream uses of those symbols are unchanged, and the new assertions confirm rather than alter the existing values. No downstream unit can observe a different result from unit 067 after this iteration.

## Side observations (non-blocking)

- The orchestrator did not write `stage_067_sympy.log` or `stage_067_mathematica.log` alongside `stage_067_diff.patch`. If this is a pipeline regression (rather than intentional reliance on persisted `.txt` outputs), it may be worth investigating so future verifiers have direct exit codes rather than mtime-based inference.
- F3's documented deviation (`.rewrite(sp.cosh)`) is a reasonable workaround for SymPy's spotty native `sech` integration support; if a future SymPy upgrade resolves the direct case, the rewrite could be removed. Not a verification blocker.
- The Mathematica numeric output for `P_res` (line 32 of the .txt) ends in `...878586507986918943333334`60.`, which last differs from the hardcoded `presTarget` (`...605997534`60.`) starting around digit 36 (`...605998...` vs `...605878...`). The `expectApprox` tolerance is `10^-34`, and the printed diff (`0``49.696539346527814`) confirms agreement well inside tolerance, so this is fine — just noting the engines genuinely agree at ~36-digit level rather than full 60 digits, which is expected from NIntegrate vs mpmath.quad and the directive does not require closing this.

## Verdict justification

All four findings were applied verbatim per the directive. F1/F2/F4 are comment-only labeling fixes that correctly annotate three tautological/empirical checks without altering assertion logic, and the labels match the directive text exactly. F3 adds a substantive, non-tautological SymPy integration that independently derives the sech^2 and Gaussian norms and compares them against the declared values (residuals reported as `= 0` in the output), closing the asymmetry where only the Mathematica side anchored the norms. The single documented deviation (`.rewrite(sp.cosh)`) is a benign environment workaround that produces the correct closed forms, as the output transcript verifies. Exec log files are absent, but persisted output `.txt` files are fresher than their scripts, contain all expected PASS / `= 0` lines, and complete with the FINAL LEDGER / "audit passed" footers — strong indirect evidence of exit=0 in both engines. No regressions, no tautological replacements, no collateral edits. Verdict: `verified`.
