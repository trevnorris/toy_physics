---
unit_id: 050
batch: III.2
verifier_model: claude-opus-4-7-1m
verify_date: 2026-05-22T00:00:00Z
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 4
findings_total: 4
material_change: false
---

# Verification — unit 050

## Per-finding outcomes

### F1 — tautological_check (suppression-factor cancellation)

**Classification:** resolved

**What changed:**
- `scripts/moving_throat_pde_stage050_zeta_threshold_comparison_sympy_audit.py:77-86` — the trivial `(2n+1)^2 zeta_n - 1/(1+x n(n+1))` cancellation is gone. It is replaced with an admissibility residual:
  `admissibility_num = sp.together(sp.simplify(-x_eq * n * (n + 1) - (((2 * n + 1) ** 2 * zeta_req - 1) / ((2 * n + 1) ** 2 * zeta_req))))`
  asserted to zero, followed by the print `"Therefore the necessary condition is exactly zeta_req <= 1/(2n+1)^2."`.
- `mathematica/moving_throat_pde_stage050_zeta_threshold_comparison_mathematica_audit.wl:72-80` — the corresponding mirror builds `admissibilityNum = FullSimplify[-xEq n (n + 1) - ((2 n + 1)^2 zetaReq - 1)/((2 n + 1)^2 zetaReq), ...]` and asserts zero under the same label.

**Assessment:**
The check is now substantive: it ties the directional inequality `(2n+1)^2 zeta_req <= 1` to the algebraic structure of `x_max` rather than self-cancelling `(2n+1)^2`. The residual is non-trivial — both engines' outputs show the intermediate residual being simplified to 0. The directive flagged the literal snippet wouldn't simplify cleanly and noted the corrected sign (`(2n+1)^2 zeta_req - 1` vs. `1 - (2n+1)^2 zeta_req`); the applied form is `-x_eq*n(n+1) - ((2n+1)^2 zeta_req - 1)/((2n+1)^2 zeta_req)`, which is the algebraically equivalent corrected residual that the directive's "Applied" deviation block explicitly authorizes. Non-tautological.

### F2 — tautological_check (Mathematica xEq self-definition)

**Classification:** resolved

**What changed:**
- `mathematica/moving_throat_pde_stage050_zeta_threshold_comparison_mathematica_audit.wl:52-56` — `xEq` is now derived via `Solve`:
  `xEqSolution = x /. First[Solve[zetaN == zetaReq, x]];`
  `xEqSolution = xEqSolution /. ConditionalExpression[e_, _] :> e;`
  `xEq = FullSimplify[xEqSolution, Assumptions -> $Assumptions];`
  `xEqClosedForm = (1/(((2 n + 1)^2) zetaReq) - 1)/(n (n + 1));`
- Line 67-70: the assertion is now `expectZero["x_max (from Solve) - [1/((2n+1)^2 zeta_req)-1]/[n(n+1)]", xEq - xEqClosedForm];` — comparing the Solve-derived value to the independent closed form.

**Assessment:**
The Mathematica `Solve[zetaN == zetaReq, x]` step now mirrors SymPy's `sp.solve(Eq(zeta_n, zeta_req), x)`. The new check is independent (Solve's output vs. the manually-written closed form), no longer a definitional identity. The orchestrator's two-line `ConditionalExpression` strip (inside `xEqSolution` and inside `expectZero`) is a mechanical idiom adjustment — it does not weaken the substance, because under `$Assumptions` (`sReq>0, 0<eps<1, x>0, n>=1 integer`) the `ConditionalExpression[0, cond]` form is identically zero on the declared domain. The output shows `PASS: x_max (from Solve) - [1/((2n+1)^2 zeta_req)-1]/[n(n+1)]`. The substantive Solve derivation is preserved. The pre-existing line `expectZero["zeta_n^(twin)(x_max) - zeta_req", (zetaN /. x -> xEq) - zetaReq]` is kept (now line 66) and still substantively confirms the solved x satisfies the threshold equation.

### F3 — insufficient_verification (ceiling claim)

**Classification:** resolved

**What changed:**
- `scripts/moving_throat_pde_stage050_zeta_threshold_comparison_sympy_audit.py:95-107` — after the `S_n^(twin)(x=0) - S_n^(max)` check, the script now computes `ceiling_diff = sp.simplify(S_n_max - S_n)` and asserts it equals the explicit factored target `(1-eps)(2n+1)^2 n(n+1) x / [((2n+1)^2 - eps)((2n+1)^2(1 + x n(n+1)) - eps)]`. Comment notes that under `0 < eps < 1, n >= 1, x >= 0` every factor is non-negative.
- `mathematica/moving_throat_pde_stage050_zeta_threshold_comparison_mathematica_audit.wl:87-95` — mirrored: `ceilingDiff = FullSimplify[sNMax - sN, ...]` asserted against the same closed form `ceilingDiffTarget`.

**Assessment:**
Both checks are substantive: they assert that the difference factors into a product of factors with manifest sign under the stated assumptions, rather than only matching at `x=0`. The expected closed form is independent of the input expressions (it's a separately-constructed expression compared against the simplified difference). Both outputs show 0 residual and PASS. Non-tautological.

### F4 — insufficient_verification (monotonicity / direction)

**Classification:** resolved

**What changed:**
- `scripts/moving_throat_pde_stage050_zeta_threshold_comparison_sympy_audit.py:61-68` — after `zeta_n = twin_support_ratio(n, x)`, the script computes `d_zeta_n_dx = sp.simplify(sp.diff(zeta_n, x))` and asserts it equals `-n(n+1)/((2n+1)^2 (1 + x n(n+1))^2)` (manifestly negative on the domain).
- `mathematica/moving_throat_pde_stage050_zeta_threshold_comparison_mathematica_audit.wl:59-64` — mirrored with `D[zetaN, x]` against the same target.

**Assessment:**
The derivative-sign check pins down monotonic decrease in `x`, which fixes the direction of the inequality `zeta_n^(twin) >= zeta_req iff x <= x_max`. The target expression is independently written; comparison is non-trivial. Both outputs show 0 residual and PASS.

## Exec log assessment

**SymPy:** exit=0 (inferred from the saved txt file ending "All Stage-33 symbolic checks passed." and the `expect_zero` helper that would raise on a non-zero residual). Notable lines:
- `S(1;eps) - 2 = 0`
- `d zeta_n^(twin) / dx + n(n+1) / [(2n+1)^2 (1 + x n(n+1))^2] = 0`
- `admissibility numerator residual = 0`; `x_max non-negativity reduces to (2n+1)^2 zeta_req <= 1 = 0`
- `S_n^(max) - S_n^(twin) factored form = 0`

**Mathematica:** exit=0 (saved txt ends "Stage 050 Mathematica audit passed." which is printed before `Exit[0]`, and every `expectZero` reports `PASS:`). Notable lines:
- `PASS: S(1;eps) - 2`
- `PASS: d zeta_n^(twin) / dx + n(n+1)/[(2n+1)^2 (1 + x n(n+1))^2]`
- `PASS: x_max (from Solve) - [1/((2n+1)^2 zeta_req)-1]/[n(n+1)]`
- `PASS: x_max non-negativity reduces to (2n+1)^2 zeta_req <= 1`
- `PASS: S_n^(max) - S_n^(twin) factored form`

**Output freshness:** Confirmed. SymPy script mtime = 2026-05-22 16:59:09, output txt mtime = 2026-05-22 17:01:11 (output newer). Mathematica script mtime = 2026-05-22 17:24:49 (the orchestrator's manual `ConditionalExpression` patch), output txt mtime = 2026-05-22 17:25:57 (output newer than the patched script). Both outputs are post-fix.

## Material-change assessment

`material_change`: false.

The stage 050 audit script encodes already-derived identities; none of the edits introduce new symbolic expressions that downstream stages import. The only import-bearing object is `twin_support_ratio` from stage 049, which is unchanged. The new assertions (admissibility residual, ceiling factorization, derivative-sign) close the directive's findings without altering any exported result. No downstream stale-flag concern beyond the orchestrator's blanket `upstream_stale` marker.

## Side observations (non-blocking)

- The SymPy docstring (lines 1-11) and the Mathematica banner still refer to "Stage 33"/"STAGE 033". The script and output filenames use `stage050`. This is cosmetic and unrelated to any finding; not a blocker.
- The Mathematica `expectZero` helper now applies `ConditionalExpression -> e` stripping followed by a second `FullSimplify`. Under the declared `$Assumptions` this is sound for the current script, but downstream scripts that reuse this idiom should be aware the helper now treats `ConditionalExpression[0, _]` as 0 unconditionally; this is correct only when assumptions cover the condition. Not relevant to F1-F4.
- The SymPy variable `x` is declared `positive=True, real=True` (line 39), which lets `sp.simplify` use that domain in the new derivative and admissibility checks. Consistent with the assumptions on the Mathematica side.

## Verdict justification

All four findings are `resolved`. The tautological `(2n+1)^2 / (2n+1)^2` cancellation and the Mathematica `xEq` self-definition tautology are both gone, replaced by substantive checks (admissibility residual via `-x_eq n(n+1)`, and `Solve`-derived `xEq` compared to the independent closed form). The ceiling and direction-of-inequality claims now have explicit factored / signed-derivative assertions in both engines. Both engines exit 0, outputs are fresh, and the orchestrator's manual `ConditionalExpression` idiom patch is mechanical and preserves the substance of Codex's Solve-based derivation. No regressions visible. Verdict: `verified`.
