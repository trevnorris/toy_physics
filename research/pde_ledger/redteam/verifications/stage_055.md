---
unit_id: 055
batch: III.2
verifier_model: claude-opus-4-7[1m]
verify_date: 2026-05-26
verdict: verified
sympy_exit: n/a
mathematica_exit: n/a
findings_resolved: 2
findings_total: 2
material_change: false
---

# Verification — unit 055

## Per-finding outcomes

### F1 — tautological_check

**Classification:** resolved

**What changed:**
`mathematica/moving_throat_pde_stage055_explicit_reachability_mathematica_audit.wl:51` now reads
`xFloor = FullSimplify[x /. First[Solve[zetaMax == zetaReq, x]], Assumptions -> $Assumptions];`
(was the hardcoded `FullSimplify[4 - Pi^2/zetaReq, Assumptions -> zetaReq > 0]`). Lines 55, 58, 59, 60 are unchanged, matching the directive's "only line 51 changes" instruction. The diff at `redteam/exec_logs/stage_055_diff.patch` confirms exactly this single-line replacement.

**Assessment:**
The edit matches the directive verbatim, including the use of `$Assumptions` (which carries `0 < x < 4 && zetaReq > 0` from line 36, so the `4 zetaReq > Pi^2` condition from `Solve` is preserved as a `ConditionalExpression` wrapper). The saved Mathematica transcript at `mathematica/output/moving_throat_pde_stage055_explicit_reachability_mathematica_audit.txt:18` shows `x floor from zeta_max = zeta_req = ConditionalExpression[4 - Pi^2/zetaReq, 4*zetaReq > Pi^2]` — i.e. Solve produced the closed form independently of any hand-typed literal. Line 58's `expectZero` then strips the `ConditionalExpression` via the existing wrapper rule on script line 26 and reports `x floor = 4 - Pi^2/zeta_req = 0` / `PASS` (output lines 23-24). The assertion is no longer tautological: it now compares the Solve-derived expression against the paper-stated literal, so a slip in `zetaMax`, `aK`, or `omegaExp` would surface in the residual. Line 59's substantive cross-check still passes (output lines 25-26). No collateral edits.

### F2 — symbol_assumption_error

**Classification:** resolved

**What changed:**
`scripts/moving_throat_pde_stage055_explicit_reachability_sympy_audit.py:27-31` now carries a 4-line comment block above the unchanged `sp.symbols(...)` declaration documenting the paper-stated domain (`alpha > 0, 0 < x < 4, 0 < y < pi/2, zeta_req > 0`). The diff confirms a pure-insertion change with the `sp.symbols(...)` line itself byte-identical to before.

**Assessment:**
Matches the directive's option-A documentary form verbatim. The symbol declaration itself is unchanged, so the symbolic algebra and all printed expressions/residuals are unaffected. The SymPy transcript at `scripts/output/moving_throat_pde_stage055_explicit_reachability_sympy_audit.txt` shows all four `expect_zero` residuals printing `0` and the regime-split narrative identical to the pre-fix transcript, as the directive's verification clause required. The script's domain documentation is now aligned with the Mathematica `$Assumptions` (`0 < x < 4 && y > 0`) on `.wl:36` and with the paper's stated domain. No collateral edits.

## Exec log assessment

**SymPy:** exit=n/a. No `redteam/exec_logs/stage_055_sympy.log` was captured for this iteration; fell back to the persisted transcript `scripts/output/moving_throat_pde_stage055_explicit_reachability_sympy_audit.txt` (mtime 2026-05-26 03:02, script mtime 03:01). Notable lines:
- `symmetric twin point = 1` (D1)
- `closure maximum = -pi**2/(x - 4)` (= `pi^2/(4-x)`, D2)
- `twin value = 0`, `closure maximum = 0`, `x floor = 4 - pi^2/zeta_req = 0`, `KX/KW equivalence = 0` — all residuals zero, script ran to completion through the regime-split narrative.

**Mathematica:** exit=n/a. No `redteam/exec_logs/stage_055_mathematica.log` captured; fell back to `mathematica/output/moving_throat_pde_stage055_explicit_reachability_mathematica_audit.txt` (mtime 2026-05-26 03:02, script mtime 03:01). Notable lines:
- `x floor from zeta_max = zeta_req = ConditionalExpression[4 - Pi^2/zetaReq, 4*zetaReq > Pi^2]` (Solve produced the form, as required by the directive's verification clause)
- `PASS: twin value`, `PASS: closure maximum`, `PASS: x floor = 4 - Pi^2/zeta_req`, `PASS: zeta_max(x_floor) - zeta_req`, `PASS: KX/KW equivalence`
- `Stage 055 Mathematica audit passed.` (clean termination). The `Limit::alimv` warnings on output lines 9-13 are the same benign warnings flagged in the auditor report and unaffected by F1's edit.

**Output freshness:** Mathematica script mtime 2026-05-26 03:01, output mtime 03:02 (fresh by 1 min). SymPy script mtime 03:01, output mtime 03:02 (fresh by 1 min). Both outputs were regenerated after Codex's edits.

## Material-change assessment

`material_change`: false.

Neither edit alters any derived quantity. F1 rewrites the Mathematica `xFloor` from a hardcoded literal to a `Solve`-derived expression that simplifies to exactly the same closed form (`4 - Pi^2/zetaReq`), confirmed by the transcript at line 18. F2 is a comment-only insertion above an unchanged `sp.symbols(...)` call. No constants change, no symbolic algebra changes, no printed result changes. No downstream unit can observe these edits.

## Side observations (non-blocking)

- The orchestrator did not write `redteam/exec_logs/stage_055_sympy.log` or `stage_055_mathematica.log` for this iteration; only the persisted `scripts/output/` and `mathematica/output/` transcripts were available. Verification proceeded against those (mtimes confirm post-edit regeneration). The diff capture was present and matched.
- The SymPy transcript banner still reads `STAGE 38 — EXPLICIT LOWEST-LANE REACHABILITY` while the Mathematica banner reads `STAGE 038 —`. Both refer to the same unit; this is a pre-existing cosmetic inconsistency, not a finding.
- The Mathematica transcript line 18 surfaces a `ConditionalExpression[..., 4*zetaReq > Pi^2]` wrapper from `Solve`. This is benign — it is stripped by `expectZero`'s rule on `.wl:26` — and is exactly the case the project-wide Mathematica-script idiom anticipates.

## Verdict justification

Both findings were applied exactly as directed, with no collateral edits and no deviations. F1's Mathematica edit moves `xFloor` from a tautology-introducing literal to a `Solve`-derived expression while leaving lines 55, 58, 59, 60 untouched; the saved transcript confirms Solve produces the paper-stated closed form and all five `expectZero` checks still PASS. F2 is documentary and inserts the paper-domain comment block above an unchanged `sp.symbols(...)` declaration; all SymPy residuals remain zero. Outputs are fresh post-edit. No regressions in the diff. Verdict: `verified`.
