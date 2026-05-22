---
unit_id: 036
batch: II.1
verifier_model: claude-opus-4-7-1m
verify_date: 2026-05-22T00:00:00Z
verdict: verified
sympy_exit: n/a
mathematica_exit: n/a
findings_resolved: 6
findings_total: 6
material_change: false
---

# Verification — unit 036

## Per-finding outcomes

### F1 — tautological_check

**Classification:** resolved

**What changed:**
- SymPy: line 72 deleted (the `expect_zero("R_target - pi^2 A NQ/(8 beta0)", ...)` call). Diff shows the single-line removal; the `R_target` definition at line 57 and its `print` at line 67 remain intact.
- Mathematica: line 59 deleted (the `expectZero["R_target - Pi^2 A NQ/(8 beta0)", ...]` call). `rTarget` definition at line 50 and its `Print` at line 56 remain.

**Assessment:**
Edits match the directive exactly with no collateral changes. Saved outputs no longer contain the lines `R_target - pi^2 A NQ/(8 beta0) = 0` or `PASS: R_target - Pi^2 A NQ/(8 beta0)`. `R_target` still appears in the printed banner section ("R_target = pi**2*A*NQ/(8*beta0)" / "R_target = (A*NQ*Pi^2)/(8*beta0)"), so context is preserved.

### F2 — tautological_check

**Classification:** resolved

**What changed:**
- SymPy: new block inserted at lines 120-139 declaring `A_sym, beta0_sym`, building `N_sym` from the kappa expansion, then asserting `simplify(F - R_target_sym) == 0`.
- Mathematica: new block inserted at lines 124-138 with `Asym, beta0Sym`, `nSym`, `rTargetSym`, and an `expectZero["symbolic kappa derivation: F(xi,delta) - R_target_sym", ...]`.

**Assessment:**
The insertion text matches the directive's required block character-for-character in both engines. The new assertion is genuinely non-tautological: `F` is declared as a closed form in `(xi, delta)`, while `R_target_sym = N_sym * A_sym / (beta0_sym * KAPPA0_SQ)` is built from the kappa expansion `(KAPPA0_SQ*(x+deltaK) + KAPPA1_SQ*x)^4 / (KAPPA0_SQ*(A-x)*(...)^2)` and the `A_sym, beta0_sym` cancellation only works when `KAPPA0_SQ = 8/pi^2` and `KAPPA1_SQ = 16/(9*pi^2)` are both correct. Replacing either kappa coefficient with a perturbed value would break the simplification. Both saved outputs show the new PASS line at output line 27 (SymPy) / lines 37-38 (Mathematica). The pre-existing assertions at SymPy lines 68-71 and 87-90 (and Mathematica lines 59-62, 99-102) were kept as instructed; they are now bookkeeping confirmations of the factorization anchored by the new kappa identity.

### F3 — tautological_check

**Classification:** resolved

**What changed:**
- Mathematica: the `alphaCrit = FullSimplify[9*Pi^2*A*(1 + delta)/(8*(11 + 9*delta)), ...]` declaration and the `expectZero["(Pi^2 A / 8) G_max - alpha_crit", ...]` assertion are both deleted in the diff.

**Assessment:**
Diff matches directive verbatim. Saved Mathematica output no longer contains `(Pi^2 A / 8) G_max - alpha_crit = 0` or its PASS line. `alphaCrit` is not referenced elsewhere in the script, so the deletion is clean.

### F4 — tautological_check

**Classification:** resolved

**What changed:**
- SymPy: line 126 (`expect_true("inadmissible sample: R_target < 1 blocks the branch", bool(sp.Rational(9, 10) < 1), "R_target=9/10")`) deleted.
- Mathematica: line 120 (`expectTrue["inadmissible sample: R_target < 1 blocks the branch", 9/10 < 1, "R_target=9/10"]`) deleted.

**Assessment:**
Both lines absent from diff context after deletion. Outputs no longer contain the misleading `inadmissible sample: R_target < 1 blocks the branch` row in either engine. The remaining "support deficit blocks the branch" check (using real `Mmix_bad > G_sample` numeric comparison) retains genuine inadmissibility-branch coverage.

### F5 — mathematica_transliteration

**Classification:** resolved

**What changed:**
- (a) Lines 66-71 of the Mathematica script: `dGTarget` declaration replaced with derived `dGPolynomial = FullSimplify[Expand[dG*(9*delta + 11*xi)^2/9]]` and an assertion against the bare polynomial `11*xi^2 + 18*delta*xi + 9*delta^2`. A new discriminant check `discSimplified + 72*delta^2 == 0` was added.
- (b) Line 78: `gMaxTarget = FullSimplify[9*(1 + delta)/(9*delta + 11), ...]` replaced with `gMaxTarget = FullSimplify[gMax, ...]` (derived from the symbolic limit). The corresponding assertion was rewritten to `expectZero["G_max - 9(1+delta)/(9delta+11)", gMax - 9*(1 + delta)/(9*delta + 11)]`, placing the closed form on the test side rather than as a self-declared target.
- (c) Lines 128-131: `gSeriesTarget` hand-declaration removed and replaced with `Coefficient[gSeries, xi, k]` extraction for k=0,1,2, with three independent `expectZero` assertions against `0`, `1`, and `-2/(9*delta)` respectively.

**Assessment:**
All three sub-edits match the directive exactly, including the corrected discriminant target (`-72*delta^2`, not `0`). The new assertions are non-tautological:
- `dGPolynomial - (11*xi^2 + 18*delta*xi + 9*delta^2) == 0` forces Mathematica to differentiate `gTarget` and algebraically simplify the result; the coefficients 11, 18, 9 are no longer declared.
- The discriminant check confirms the numerator is sign-definite (`-72*delta^2 < 0` for `delta > 0`), corroborating positivity of `dG/dxi`. Saved output line 19 shows `discriminant (in xi) = -72*delta^2`, confirming Mathematica computed this independently.
- `gMax - 9*(1 + delta)/(9*delta + 11) == 0` puts the closed form on the test side, with `gMax` computed by `Limit[gTarget, xi -> 1, Direction -> "FromBelow"]` — a genuinely independent algebraic step.
- The series coefficients `c0, c1, c2` are read out of `gSeries` (Mathematica's `Series` expansion of `gTarget`) before being compared to the closed-form constants. Output line 48 shows `c0=0, c1=1, c2=-2/(9*delta)`, derived by Mathematica.

The Mathematica script still retains the literal `gTarget = 9*xi*(xi + delta)/(9*delta + 11*xi)` declaration at line 43, but the directive did not require deriving `gTarget` itself — only the four derived targets `dGTarget`, `gMaxTarget`, `gSeriesTarget`, and the now-removed `alphaCrit`. With `gTarget` as a shared starting point, the second engine's derivative, limit, and series machinery still must produce the same coefficients as SymPy's — a genuine independent check on the differentiation/limit/series operations.

### F6 — insufficient_verification

**Classification:** resolved

**What changed:**
Covered by the F2 insertion — the new SymPy block at lines 120-139 and Mathematica block at lines 124-138 introduce a symbolic `(xi, delta, A_sym, beta0_sym)` identity `F - R_target_sym == 0`, promoting the single-sample numeric witness at `(A=3, beta0=5, xi=1/2, delta=1)` to a parameter-wide symbolic statement.

**Assessment:**
The directive explicitly tied F6 to F2 ("If Codex has already applied F2, mark F6 as Applied: F6 with summary: covered by F2"). Codex did so. The PASS line for `symbolic kappa derivation: F(xi,delta) - R_target_sym` appears in both saved outputs, confirming the symbolic identity holds. The pre-existing single-sample host check at SymPy lines 102-118 / Mathematica lines 111-122 is retained as additional numeric corroboration.

## Exec log assessment

**SymPy:** exec_log file `redteam/exec_logs/stage_036_sympy.log` is absent. Per instructions ("the saved outputs are already fresh — read them"), I rely on the saved output at `scripts/output/moving_throat_pde_stage036_support_feasibility_frontier_sympy_audit.txt`. The output ends with `All Stage 19 checks passed.` (line 36), implying the script completed without raising any `AssertionError`. Notable assertion lines:
- `g_B,req^2/varpi^2 - (pi^2 A / 8) (G - M_mix) = 0` (line 9)
- `dG/dxi - manifestly positive form = 0` (line 15)
- `symbolic kappa derivation: F(xi,delta) - R_target_sym = 0` (line 27) — new, F2/F6 fix landed.
- `G near-onset series through O(xi^2) = 0` (line 35)
No FAIL strings and no Python tracebacks in the saved output.

**Mathematica:** exec_log file `redteam/exec_logs/stage_036_mathematica.log` is absent. Saved output at `mathematica/output/...` ends with `Stage 036 Mathematica audit passed.` (line 56), implying clean exit via `Exit[0]` (line 165 of the script). Every assertion is followed by a `PASS:` line. Notable new/changed PASS lines:
- `PASS: dG/dxi positivity polynomial: 9 dG/dxi (9d+11xi)^2 / 9 == 11 xi^2 + 18 delta xi + 9 delta^2` (line 18)
- `PASS: dG/dxi numerator discriminant equals -72 delta^2` (line 21)
- `PASS: G_max - 9(1+delta)/(9delta+11)` (line 26)
- `PASS: symbolic kappa derivation: F(xi,delta) - R_target_sym` (line 38)
- `PASS: near-onset c0 = 0`, `PASS: near-onset c1 = 1`, `PASS: near-onset c2 = -2/(9 delta)` (lines 50, 52, 54)
No FAIL strings. The removed assertions for `R_target - Pi^2 A NQ/(8 beta0)`, `(Pi^2 A / 8) G_max - alpha_crit`, and `inadmissible sample: R_target < 1 blocks the branch` are absent from the output, as required.

**Output freshness:**
- `scripts/moving_throat_pde_stage036_..._sympy_audit.py` mtime: 2026-05-21 17:37
- `scripts/output/moving_throat_pde_stage036_..._sympy_audit.txt` mtime: 2026-05-21 17:40
- `mathematica/moving_throat_pde_stage036_..._audit.wl` mtime: 2026-05-21 17:37
- `mathematica/output/moving_throat_pde_stage036_..._audit.txt` mtime: 2026-05-21 17:40

Both saved outputs are 3 minutes newer than the corresponding scripts, confirming the outputs were re-generated post-fix.

## Material-change assessment

`material_change`: false.

The edits remove tautological assertions and strengthen the verification (by deriving targets, adding discriminant/coefficient extractions, and adding a symbolic kappa cross-check). They do not change any derived quantity that downstream stages consume: `G(xi,delta)`, `F(xi,delta)`, `R_target`, `M_mix`, `dG/dxi`, `G_max`, the near-onset coefficients, and the host-sample numeric values are all unchanged. No exported constants are touched. Downstream stages depending on stage 036's claim manifest see the same numerical and symbolic content; only the proof of those claims is improved.

## Side observations (non-blocking)

- The SymPy script still hand-declares `dG_target` at line 75 and `G_series_target` at line 153 in the original closed-form style. The F5 directive limited the derivation-not-declaration restructuring to the Mathematica script (consistent with the cross-engine intent: forcing the second engine to derive independently). This scope is explicit in the directive, so no rework is needed. If future audit batches want to harden SymPy similarly, that would be a new finding.
- The Mathematica script's pre-existing `dG/dxi` displayed form `9*(1 + (18*delta^2)/(9*delta + 11*xi)^2)/11` (output line 15) is algebraically equivalent to but textually different from SymPy's `9*(-11*xi*(delta + xi) + (delta + 2*xi)*(9*delta + 11*xi))/(9*delta + 11*xi)**2` (sympy output line 14). Both reduce to the polynomial `(11 xi^2 + 18 delta xi + 9 delta^2) / (9 delta + 11 xi)^2 * (9/9)` form, so the engines agree on the underlying derivative. The textual divergence is a healthy sign of independent simplification routes.

## Verdict justification

All six findings were applied as directed. The deletions in F1, F3, F4 removed assertions that were `X - X = 0` by construction; the saved outputs no longer reference them. F2/F6 introduced a single symbolic kappa-based identity `F - R_target_sym == 0` in both engines whose truth depends on the specific kappa coefficients, providing genuine non-tautological coverage of the support-feasibility frontier across the full parameter range. F5 restructured the Mathematica derivative, endpoint, and series checks so the closed-form targets are now derived from Mathematica-native operations (differentiation, limit, series) rather than declared up front, plus an extra discriminant check anchors the positivity claim. The saved outputs are fresh (3 minutes newer than scripts), all assertions PASS, no tautological residue remains, and no downstream-visible quantities changed.
