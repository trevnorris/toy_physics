---
unit_id: 019
batch: I.2
verifier_model: claude-opus-4-7[1m]
verify_date: 2026-05-25T22:15:00Z
verdict: verified
sympy_exit: 0
mathematica_exit: n/a
findings_resolved: 1
findings_total: 1
material_change: false
---

# Verification — unit 019

## Per-finding outcomes

### F1 — symbol_assumption_error

**Classification:** resolved

**What changed:**
`scripts/moving_throat_pde_stage019_parent_throat_action_isotropic_bundle_sympy_audit.py:25-29` — symbol declarations updated exactly per the directive. The diff at `redteam/exec_logs/stage_019_diff.patch` shows only those 5 lines changed:

- Lines 25-27: `KSigma, MSigma`, `B0..Z4`, `N0, N2, N4` now carry `real=True, nonzero=True` (previously `nonzero=True` only).
- Line 28: `mhat0, G, cs, a, c` now declared `positive=True` (was `nonzero=True`).
- Line 29: `eps` now `real=True, nonzero=True`.

The directive's prohibition against making KSigma/MSigma/B*/Z*/N* positive was respected — only the dimensionful prefactors got `positive=True`, preserving the response-sign criterion's ability to genuinely sweep `D_0` and `MSigma` signs. No collateral edits anywhere else in the file or repository.

**Assessment:**
The edit matches the directive's "After" block character-for-character. The Codex `## Applied: F1` block correctly reports `deviation: none`.

Saved-output checks against the directive's verification criteria (canonical `scripts/output/moving_throat_pde_stage019_parent_throat_action_isotropic_bundle_sympy_audit.txt`):

- Output line 26: `positive-root numeric u2 = 0.6324555320336759` — unchanged.
- Output line 27: `negative-root numeric u2 = -0.6324555320336759` — unchanged.
- Output line 22: `constant-prefactor mutation guards = PASS` — preserved.
- Output line 32: `one-pole numerical response-sign guard = PASS` — preserved.
- Output line 37: `STATUS: PASS` — preserved.

Output line 23 (`u2_on_positive_root = -sqrt(3)*sqrt(-(B4 + Z4)*(B0 - KSigma + Z0))/(3*B0 - 3*KSigma + 3*Z0)`) retains its pre-fix form; the directive explicitly permits this line to either collapse or remain — non-collapse is acceptable because SymPy's simplifier does not automatically rewrite the `sqrt(-...)` even with positive prefactors when the radicand involves real symbols of unconstrained sign (`B4+Z4`, `B0-KSigma+Z0`, which remain `real=True, nonzero=True` per the directive's deliberate choice). This is consistent with the finding's intent: assumption hygiene is now in place for future-proofing without changing any current assertion.

The edit strengthens the symbol domain only and does not introduce or weaken any check, so non-tautology is preserved trivially.

## Exec log assessment

**SymPy:** exit=0 (inferred from `STATUS: PASS` trailer on the canonical saved output at `scripts/output/moving_throat_pde_stage019_parent_throat_action_isotropic_bundle_sympy_audit.txt`, mtime 2026-05-25 22:08, which is newer than the script mtime 2026-05-25 22:07). The orchestrator's `redteam/exec_logs/stage_019_sympy.log` (mtime May 21 14:58) is stale relative to the fix and was not used, per the prompt's explicit instruction to fall back to canonical outputs. Notable lines from the canonical output:

```
constant-prefactor mutation guards = PASS
one-pole numerical response-sign guard = PASS
STATUS: PASS
```

All M-root, K-from-one-pole, K-from-normalization, P2/P4 closures, mutation-guard residuals, and three numeric response-sign samples report PASS.

**Mathematica:** exit=n/a. The directive scope was SymPy-only; the `.wl` file (mtime May 21 13:35) and its output (mtime May 21 15:02) are both untouched and were already passing per the original auditor report (`STATUS: PASS`, all M1-M12 OK). No regression possible.

**Output freshness:** SymPy script mtime is 2026-05-25 22:07; saved SymPy output mtime is 2026-05-25 22:08 — output is one minute newer than the script, confirming re-generation post-fix. Mathematica script/output pair predate the fix consistently (script May 21 13:35 < output May 21 15:02) and were correctly out of scope.

## Material-change assessment

`material_change`: false.

No derived result changed. The closed-form expressions printed in the saved output (K_from_one_pole, K_from_norm, P0_target, P2, P4, N2_const, N4_const, compatibility, M-root expressions, numeric u2 values, Gaussian wall result) are bit-identical to the pre-fix output the auditor reviewed. Only SymPy symbol-domain metadata changed inside the script. Each stage's audit script is self-contained — no downstream unit consumes this script's outputs — so no narrow upstream-stale concern beyond the orchestrator's default `upstream_stale: true` flag.

## Side observations (non-blocking)

- The orchestrator's `exec_logs/stage_019_sympy.log` is stale (May 21) relative to the post-fix run (May 25 22:08). The prompt directed me to the canonical `scripts/output/...txt` instead, which I used; flagging only so the orchestrator can decide whether to refresh exec_log capture for future audits.
- The directive's `applied_at` is `2026-05-26T04:07:07Z` while today is `2026-05-25` — UTC vs. local timezone artifact, not substantive.

## Verdict justification

The single finding (low-severity symbol-assumption hygiene) is fully resolved. Codex's edit matches the directive's "After" block character-for-character, touches only the 5 specified lines, makes no collateral changes, and the post-fix saved output preserves every PASS line and every numeric value the directive required to remain identical. Mathematica was correctly left untouched per directive scope. No regressions, no material change to derived results, no new findings.
