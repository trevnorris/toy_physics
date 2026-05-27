---
unit_id: 067
batch: III.3
verifier_model: claude-opus-4-7-1m
verify_date: 2026-05-26T00:00:00Z
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 1
findings_total: 1
material_change: false
---

# Verification — unit 067 (second pass)

## Per-finding outcomes

### F1 — paper_misalignment (banner relabel STAGE 50/050 → STAGE 067)

**Classification:** resolved

**What changed:**
The orchestrator applied this finding directly (codex declined because the directive carried a stale `## Resolve before fix_loop` block from the auto-generated user-resolution skeleton; the audit report itself stated the direction was unambiguous, so the orchestrator made the three string replacements without re-invoking codex).

Three string replacements, exactly matching the audit report's "Required change" block:

- `scripts/moving_throat_pde_stage067_sech_gaussian_sympy_audit.py:4`: docstring filename changed from `moving_throat_pde_stage50_sech_gaussian_sympy_audit.py` to `moving_throat_pde_stage067_sech_gaussian_sympy_audit.py` (confirmed in the file at line 4 and in the diff at `stage_067_diff.patch` lines 22-23).
- `scripts/moving_throat_pde_stage067_sech_gaussian_sympy_audit.py:53`: banner argument changed from `"STAGE 50 — EXACT SECH–GAUSSIAN COHERENCE BENCHMARK"` to `"STAGE 067 — EXACT SECH-GAUSSIAN COHERENCE BENCHMARK"` (confirmed in the file at line 53 and in the diff at lines 31-32).
- `mathematica/moving_throat_pde_stage067_sech_gaussian_resonance_mathematica_audit.wl:38`: banner argument changed from `"STAGE 050 — EXACT SECH-GAUSSIAN COHERENCE BENCHMARK"` to `"STAGE 067 — EXACT SECH-GAUSSIAN COHERENCE BENCHMARK"` (confirmed in the file at line 38 and in the diff at lines 9-10).

**Assessment:**
All three replacements are label-only. The diff (`redteam/exec_logs/stage_067_diff.patch`) shows exactly the three single-line changes described — no collateral edits, no whitespace churn beyond the targeted strings, and no modifications outside lines 4/53 of the .py and line 38 of the .wl. The audit's recommended target strings used an ASCII hyphen `STAGE 067 — EXACT SECH-GAUSSIAN ...` (note the en-dash in the prose `—` separator with ASCII `-` inside `SECH-GAUSSIAN`); the applied banners in both engines match that exact form, including in the case of the sympy banner where the prior text used `SECH–GAUSSIAN` (en-dash) and the new text uses `SECH-GAUSSIAN` (ASCII hyphen). The mathematica banner already used `SECH-GAUSSIAN`, so it remains stable in that respect.

Assertion content: no assertion was added, removed, or modified. The sympy script still defines `expect_zero` (line 46) and uses it identically at the lines cited in the audit (A1-A9). The mathematica script still defines `expectZero`/`expectApprox`/`expectTrue` (lines around 19-36) and uses them at B1-B9. The diff does not touch any line containing `expect_zero`, `expectZero`, `expectApprox`, `expectTrue`, `raise`, `Solve`, `NIntegrate`, `Integrate`, `mp.quad`, or any numeric tolerance.

Transcript banners post-fix (confirmed in fresh exec logs):

- `redteam/exec_logs/stage_067_sympy.log:8`: `STAGE 067 — EXACT SECH-GAUSSIAN COHERENCE BENCHMARK` (was `STAGE 50 — EXACT SECH–GAUSSIAN COHERENCE BENCHMARK`).
- `redteam/exec_logs/stage_067_mathematica.log:8`: `STAGE 067 — EXACT SECH-GAUSSIAN COHERENCE BENCHMARK` (was `STAGE 050 — EXACT SECH-GAUSSIAN COHERENCE BENCHMARK`).

Both engines still print all the expected interior `PASS` / `= 0` lines and complete with the `FINAL LEDGER` banner (sympy) and `Stage 067 Mathematica audit passed.` footer (mathematica). The mathematica .wl is no longer internally inconsistent (banner at line 38 and closing print at line 174 both now use "067").

Verdict: this is a textual-label fix only, with the expected scope (three lines, three files including the docstring), and all surrounding behavior is preserved.

## Exec log assessment

**SymPy:** exit=0. Notable lines from `redteam/exec_logs/stage_067_sympy.log`:

- Line 8: `STAGE 067 — EXACT SECH-GAUSSIAN COHERENCE BENCHMARK` (confirms the relabel).
- Line 18: `N_(sigma sigma) integral - 2 w_f = 0`
- Line 19: `N_(phi phi) integral - w_g sqrt(pi/2) = 0`
- Line 25: `C^2(r) - C^2(pi/r) under duality = 0`
- Line 32: `self-dual overlap-slope relation = 0`
- Line 33: `stationary derivative of C^2 at the self-dual point = 0`
- Lines 48-52: five duality samples at `r = 0.75, 1.0, 1.2, 1.5, 2.0` evaluate to `0.0` or `~1.6e-61`.
- Line 87: `# exit_code: 0`.

**Mathematica:** exit=0. Notable lines from `redteam/exec_logs/stage_067_mathematica.log`:

- Line 8: `STAGE 067 — EXACT SECH-GAUSSIAN COHERENCE BENCHMARK` (confirms the relabel).
- Lines 13-15: `PASS: N_(sigma sigma) - 2 w_f`, `PASS: N_(phi phi) - w_g sqrt(pi/2)`.
- Line 22: `PASS: C^2(r) - C^2(pi/r) under duality`.
- Line 30: `PASS: self-dual C^2 stationary slope from symmetry solve`.
- Lines 40, 42: `PASS: C_res^2 numeric check`, `PASS: P_res numeric check`.
- Lines 48-61: all five `PASS: duality sample ...`.
- Lines 74, 81: `PASS: constructive-branch increase up to r_*` and `PASS: constructive-branch decrease after r_*`.
- Line 83: `Stage 067 Mathematica audit passed.`
- Line 84: `# exit_code: 0`.

**Output freshness:** the orchestrator's exec logs at `redteam/exec_logs/stage_067_sympy.log` and `redteam/exec_logs/stage_067_mathematica.log` are dated `2026-05-26T13:13:08-06:00` and `2026-05-26T13:13:18-06:00` respectively, both after the script mtimes (sympy .py: 2026-05-26 13:11:12; mathematica .wl: 2026-05-26 13:11:29). So the exec logs reflect the relabel. The persisted `.txt` outputs under `scripts/output/` and `mathematica/output/` retain the prior-pass mtimes of 2026-05-22 19:56 — those snapshots were not re-saved by the orchestrator after this iteration's banner relabel. This is a cosmetic staleness in the saved transcript snapshots only; the fresh exec logs (used by this verifier) show the corrected banners, and the .txt snapshot content is otherwise identical except for the banner line. Flagged as a side observation, not a verification blocker.

## Material-change assessment

`material_change`: false.

The three edits are pure string replacements in a docstring and two banner-print arguments. No symbol, expression, integrator setting, tolerance, or assertion line was touched. Every numeric and symbolic output line in both transcripts is byte-identical to the prior pass except for the banner line that is the explicit subject of the fix. No downstream unit can observe a different result from unit 067 as a function of this edit.

## Side observations (non-blocking)

- The persisted `.txt` outputs under `scripts/output/moving_throat_pde_stage067_sech_gaussian_sympy_audit.txt` (mtime 2026-05-22 19:56) and `mathematica/output/moving_throat_pde_stage067_sech_gaussian_resonance_mathematica_audit.txt` (mtime 2026-05-22 19:56) were not refreshed after the banner relabel. They still display the prior `STAGE 50` / `STAGE 050` banner. The fresh exec logs in `redteam/exec_logs/` show the corrected banner, and the orchestrator likely treats the exec logs as the authoritative post-fix transcript, but if downstream review consumes the `scripts/output/` and `mathematica/output/` snapshots directly they will need to be regenerated (a simple re-run of the script, no math changes). Not a blocker for the verification because the audit explicitly authorized only the three label edits and the exec logs confirm correctness.
- The audit's recommended replacement text deliberately uses an ASCII hyphen in `SECH-GAUSSIAN` rather than the prior en-dash; both engines now use the ASCII form consistently. No issue.
- The notes file `notes/stages/moving_throat_pde_stage067_sech_gaussian_resonance.md:2` still carries the "Stage 50" header. The audit and the directive both flagged this as out of red-team scope (notes are not edited by the red-team), so it is not part of this finding and is not a blocker. The user can address it out-of-band.

## Verdict justification

The orchestrator applied the three banner/docstring relabels exactly as the audit report's "Required change" block specified. The diff confirms no collateral edits beyond those three lines, no assertion or numeric content changed, and both fresh exec logs show the corrected `STAGE 067 — EXACT SECH-GAUSSIAN COHERENCE BENCHMARK` banner with full `PASS` / `= 0` interior content and exit code 0. The codex-decline → orchestrator-apply path is appropriate here because the audit itself stated the direction was unambiguous (paper card and filenames are authoritative on "067") and the change is a label fix with no mathematical content. Verdict: `verified`, `material_change: false`.
