---
unit_id: 033
batch: II.1
verifier_model: claude-opus-4-7[1m]
verify_date: 2026-05-26T00:30:00Z
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 1
findings_total: 1
material_change: false
---

# Verification — unit 033 (batch II.1 v2)

## Per-finding outcomes

### F1 — tautological_check (Stage 16.6 decorative final identity + redundant mma numeric block)

**Classification:** resolved

**What changed:**

- `scripts/moving_throat_pde_stage033_microscopic_normalization_equation_sympy_audit.py:125-128`: the `expect_zero` label was rewritten from `"alpha_crit(mic) - alpha_0(mic) - gate_num_target/gate_den_claim"` to the two-line label `"gate_num_target/gate_den_claim - (alpha_crit_mic - alpha_0_mic) " "(tautological by reconstruction; substantive check is den_ratio.is_number above)"`. The residual expression (`sp.simplify(alpha_crit_mic - alpha0_mic - gate_num_target / gate_den_claim)`) is unchanged.
- `mathematica/moving_throat_pde_stage033_microscopic_normalization_equation_mathematica_audit.wl:128-131`: same label change applied to the `expectZero` call (now `"gate_num_target/gate_den_claim - (alpha_crit_mic - alpha_0_mic) (tautological by reconstruction; substantive check is NumericQ[denRatio] above)"`). Residual expression unchanged.
- `mathematica/moving_throat_pde_stage033_microscopic_normalization_equation_mathematica_audit.wl` `Do[...]` block at lines 146-154: the `gateNumeric = N[(alphaCritMic - alpha0Mic - gateNumTarget/gateDenClaim) /. rule, 30];` print + `If` sub-block (six lines) was removed. The `monotonicityNumeric` block and the iterator `{rule, {numericRule1, numericRule2}}` are intact.

**Assessment:**
The diff (`redteam/exec_logs/stage_033_diff.patch`) shows exactly two label-replacement hunks (one in `.py`, one in `.wl`) plus one deletion hunk (six lines removed from the `.wl` `Do[...]` body). Nothing else was touched. The relabelling matches the directive's "after" snippets character-for-character (modulo the Python implicit-string-concatenation line split, which preserves identical runtime string content). The Mathematica `gateNumeric` excision matches the directive's "after" snippet — the `Do[...]` body now contains only the `monotonicityNumeric` block, terminated by `,` followed by the iterator spec. No collateral changes appear.

The substantive verification load — the `den_ratio.is_number` guard (sympy line 120) and `NumericQ[denRatio]` guard (mma line 122) — remains untouched and still passes (`den_ratio = 9*pi**2` / `-9*Pi^2`). The structurally distinct `monotonicityNumeric` numerical check at two rational test points is preserved and still passes (`monotonicity numeric residual = 0``78.83...` and `0``78.99...`, both within the `10^-20` tolerance). The relabelled `expect_zero`/`expectZero` continues to pass (the residual is identically zero by construction, which is precisely what the new label discloses).

The directive's verification criteria (a)–(e) are all met:
- (a) sympy label at 125-128 now contains `"tautological by reconstruction"` — confirmed by file read and by sympy exec log line 138.
- (b) mma label at 128-131 contains the same substring — confirmed by file read and mma exec log lines 42-43.
- (c) the mma `Do[...]` block no longer contains a `gateNumeric` block — confirmed by file read and by absence of any `gate-identity numeric residual` line in the mma exec log.
- (d) both scripts exit 0 (sympy log `# exit_code: 0` at line 188; mma log `# exit_code: 0` at line 77).
- (e) sympy transcript still shows zero residual under the new label (line 138); mma transcript shows `PASS: gate_num_target/gate_den_claim - (alpha_crit_mic - alpha_0_mic) (tautological by reconstruction; substantive check is NumericQ[denRatio] above)` (line 43) and `PASS: monotonicity numeric residual zero at rule` twice (lines 57, 74) with no `gate-identity numeric residual` lines.

The fix follows the directive's option (a) (cosmetic relabel + drop redundant numeric block) — the minimum, safe fix. No new tautology is introduced; the existing tautology is honestly disclosed rather than disguised. The `Applied: F1` block declares `deviation: none`, which is accurate.

## Exec log assessment

**SymPy:** exit=0. Notable lines:
- `denominator ratio (must be parameter-free) = 9*pi**2` (line 137) — substantive guard still operative.
- `gate_num_target/gate_den_claim - (alpha_crit_mic - alpha_0_mic) (tautological by reconstruction; substantive check is den_ratio.is_number above) = 0` (line 138) — relabelled assertion present and passing.
- `All Stage 16 checks passed.` (line 187), `# exit_code: 0` (line 188).

**Mathematica:** exit=0. Notable lines:
- `denominator ratio (must be parameter-free) = -9*Pi^2`, `PASS: gate denominator matches claim up to parameter-free constant` (lines 39-40) — substantive guard intact.
- `PASS: gate_num_target/gate_den_claim - (alpha_crit_mic - alpha_0_mic) (tautological by reconstruction; substantive check is NumericQ[denRatio] above)` (line 43).
- `PASS: monotonicity numeric residual zero at rule` appears twice (lines 57, 74). No `gate-identity numeric residual` line appears anywhere — the redundant block is gone.
- `N::meprec` precision-limit warnings (lines 52, 68) are informational; the printed residuals (`0``78.83...`, `0``78.99...`) are far below the `10^-20` tolerance.
- `Stage 033 Mathematica audit passed.` (line 76), `# exit_code: 0` (line 77).

**Output freshness:** The orchestrator-captured exec logs (`redteam/exec_logs/stage_033_sympy.log` mtime 2026-05-26T00:19:43, `stage_033_mathematica.log` mtime 2026-05-26T00:19:59) are post-edit (script mtime 2026-05-25T23:49:52). The saved transcripts under `scripts/output/...txt` (mtime 2026-05-21T17:31:40) and `mathematica/output/...txt` (mtime 2026-05-21T17:31:50) are stale relative to the v2 edits — a grep against those files still shows the pre-fix label `alpha_crit(mic) - alpha_0(mic) - gate_num_target/gate_den_claim`. This is an orchestration housekeeping issue rather than a verification failure: the prompt directs the verifier to use the orchestrator-captured exec logs (which are fresh and post-fix), and the directive correctly instructed Codex not to run python/mathematica. Flagged under side observations.

## Material-change assessment

`material_change`: false.

The edits are label-text changes plus removal of a numerical residual block that was itself tautological by construction. No derived result, no symbolic identity, no closed-form value, and no substantive assertion changed. The `den_ratio = ±9*pi^2` ratio that the substantive guard verifies is identical before and after the patch. Downstream units depending on stage 033's symbolic outputs (the microscopic `alpha_crit_mic`, `K0_onset` closed form, weak-loading coefficients) see no change in computed values. The orchestrator can leave units > 033 unaffected by this verification — no targeted re-audit is required.

## Side observations (non-blocking)

1. The saved transcripts at `scripts/output/moving_throat_pde_stage033_microscopic_normalization_equation_sympy_audit.txt` and `mathematica/output/moving_throat_pde_stage033_microscopic_normalization_equation_mathematica_audit.txt` are stale (mtime 2026-05-21) and still contain pre-fix labels. The orchestrator's exec-log capture is fresh and authoritative; a follow-up `redteam exec-sympy 033` / `redteam exec-mathematica 033` run would refresh those `output/*.txt` artefacts. This does not block verification.

2. The mma `Do[...]` body now contains a single statement (the `monotonicityNumeric` block) followed by the iterator. The trailing `,` after the `If[...]` and before `{rule, ...}` is preserved exactly as in the directive's "after" snippet — Mathematica syntax remains valid.

3. The Python label is a two-string implicit concatenation (`"...load1..." "...load2..."`); its runtime value is the single string the directive specified. The sympy exec log line 138 prints the concatenated label exactly as expected.

## Verdict justification

Finding F1 was the sole finding in the v2 auditor's report and was applied exactly as the directive specified, with `deviation: none`. Both engine scripts exit 0; the substantive `den_ratio` guard and the structurally distinct `monotonicityNumeric` cross-check remain operative; the relabelled assertion now self-discloses its tautological status; the redundant Mathematica `gateNumeric` block is excised. The diff is minimal (two label hunks + one six-line deletion) with no collateral edits. No regressions, no new tautologies, no material change to any derived result. The patch correctly addresses the documentation-level risk the finding flagged.
