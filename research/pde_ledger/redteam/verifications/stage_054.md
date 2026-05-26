---
unit_id: 054
batch: III.2
verifier_model: claude-opus-4-7[1m]
verify_date: 2026-05-26
verdict: verified
sympy_exit: n/a
mathematica_exit: n/a
findings_resolved: 1
findings_total: 1
material_change: false
---

# Verification — unit 054

## Per-finding outcomes

### F1 — insufficient_verification

**Classification:** resolved

**What changed:**
- `scripts/moving_throat_pde_stage054_robin_softening_sympy_audit.py:72-79` — inserted a four-line comment block plus three statements:
  - `dAK_dy = sp.diff(AK_sym, y)`
  - `dAK_dy_expected = -2 * x * y / (pi**2 * (1 - x / 4 + x * y**2 / pi**2) ** 2)`
  - `expect_zero("dA_K/dy closed form", dAK_dy - dAK_dy_expected)`
  - a `print(...)` line stating the prefactor-positivity argument on the declared domain.
- `mathematica/moving_throat_pde_stage054_robin_softening_mathematica_audit.wl:81-88` — inserted the mirrored block:
  - `dAKdy = FullSimplify[D[aKSym, y], Assumptions -> $Assumptions];`
  - `dAKdyExpected = -2 x y / (Pi^2 (1 - x/4 + x y^2/Pi^2)^2);`
  - `expectZero["dA_K/dy closed form", dAKdy - dAKdyExpected];`
  - mirror `Print[...]` line.

Insertion points sit exactly between the existing `soft-mouth limit` assertion and the next banner / `ineqRhs` line, matching the directive's "between line 70 and line 72" (sympy) / "between line 79 and line 81" (mathematica) instruction. The `## Applied: F1` block in the directive lists both files with `deviation: none`, consistent with the captured diff (`redteam/exec_logs/stage_054_diff.patch`).

**Assessment:**
Codex's edit matches the directive verbatim — no algebraic substitutions, no helper changes, no reformatting of surrounding lines. The new assertion is non-tautological:

- `dAK_dy = sp.diff(AK_sym, y)` (or `D[aKSym, y]`) is whatever the engine's differentiation routine produces from the script's own simplified form of `AK_sym = 1/(1 - x/4 + x*y^2/pi^2)`.
- `dAK_dy_expected` is the independently-typed closed form `-u'/u^2` with `u = 1 - x/4 + x*y^2/pi^2` and `u' = 2*x*y/pi^2`.
- Equality of the two residuals holds iff the engine's chain-rule output on the *actual* denominator agrees with the chain-rule applied to the *paper's* denominator. If `AK_sym`'s denominator had a sign flip on `x/4` or on `x*y^2/pi^2`, the differentiation output would have the corresponding sign flip while `dAK_dy_expected` would not — the residual would be nonzero.

Together with the printed prefactor-positivity remark on `0 < x < 4, 0 < y < pi/2`, the new identity certifies `dA_K/dy < 0` on the declared domain, which is the load-bearing piece that closes the endpoint window `1 <= A_K <= 4/(4-x)` flagged by the auditor as exercised only at its two endpoints.

No collateral edits: the diff (`stage_054_diff.patch`) is exactly two hunks, both pure insertions, with no modification to any pre-existing assertion, print, or banner.

## Exec log assessment

**SymPy:** exit=n/a. No `redteam/exec_logs/stage_054_sympy.log` was captured by the orchestrator. Fallback per verifier prompt: `scripts/output/moving_throat_pde_stage054_robin_softening_sympy_audit.txt` (mtime 2026-05-26 03:00, newer than the script's 02:59). Notable lines from the transcript:

- `soft-mouth limit = 0`
- `dA_K/dy closed form = 0`
- `Prefactor 2*x*y/pi^2 > 0 on 0<x<4, 0<y<pi/2 => dA_K/dy < 0 (monotone decreasing).`
- `x floor = 4 - 4/zeta_req = 0`
- FINAL LEDGER block printed normally, indicating the script completed without an `expect_zero` failure (the helper raises on a nonzero residual, terminating the transcript early).

**Mathematica:** exit=n/a. No `redteam/exec_logs/stage_054_mathematica.log` was captured. Fallback: `mathematica/output/moving_throat_pde_stage054_robin_softening_mathematica_audit.txt` (mtime 2026-05-26 03:00). Notable lines:

- `PASS: soft-mouth limit`
- `dA_K/dy closed form = 0`
- `PASS: dA_K/dy closed form`
- `PASS: x floor = 4 - 4/zeta_req`
- `PASS: A_K,max(x_floor) - zeta_req`
- Final line: `Stage 054 Mathematica audit passed.`

All previously-passing assertions still pass; the new `PASS: dA_K/dy closed form` line is appended exactly where the directive predicted (immediately after `PASS: soft-mouth limit`). The pre-existing `Limit::alimv` warning at the `y -> 0+` step is unchanged and remains informational, as the auditor noted.

**Output freshness:** confirmed. Both scripts have mtime `2026-05-26 02:59`; both saved `.txt` outputs have mtime `2026-05-26 03:00`. Outputs are strictly newer than scripts post-edit.

## Material-change assessment

`material_change`: false.

The edit only adds a new symbolic identity assertion. It does not alter any expression, derived value, simplification path, or constant that downstream units could consume. `A_K(y)`, `A_K,max`, `x_floor`, the soft-mouth limit, and the Robin characteristic equation are all numerically and symbolically unchanged in both engine transcripts. Downstream units that consume the softening factor or window will see identical inputs from unit 054.

## Side observations (non-blocking)

- The orchestrator did not write `stage_054_sympy.log` or `stage_054_mathematica.log` under `redteam/exec_logs/`. The verifier prompt anticipates this case and allows falling back to `scripts/output/` and `mathematica/output/`, which I did. If the state machine requires those log files to be present, that is a tooling gap rather than a script defect.
- The Mathematica banner still prints `STAGE 037` and the SymPy banner still prints `STAGE 37`, while the filename and paper card use `stage_054`. The auditor noted this as a print-only renumbering artifact and explicitly chose not to flag it; I concur — it does not enter any assertion residual.
- The `expectZero["dA_K/dy closed form", dAKdy - dAKdyExpected]` Mathematica call leans on `FullSimplify` (with `$Assumptions` declaring `0 < x < 4, y > 0`) to reduce the residual to zero. The transcript records `dA_K/dy closed form = 0` and `PASS: dA_K/dy closed form`, confirming FullSimplify is sufficient here; no additional hardening is needed.

## Verdict justification

Codex applied F1 exactly as the directive specified — same file paths, same insertion line ranges, same code blocks, and the directive's `## Applied: F1` block records `deviation: none`. Both engines now log a non-tautological `dA_K/dy closed form = 0` residual immediately after the soft-mouth-limit check; together with the printed prefactor-positivity argument on the declared domain `0 < x < 4, 0 < y < pi/2`, this certifies monotonicity of `A_K(y)` and thereby closes the window claim `1 <= A_K <= 4/(4-x)` that the auditor identified as endpoint-only. The captured diff is clean (two pure-insertion hunks, no collateral edits), output `.txt` files are newer than the edited scripts, and every pre-existing PASS line is preserved. Verdict: `verified`.
