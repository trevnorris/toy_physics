---
unit_id: 126
batch: IV.3
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-05-29T15:40:00-06:00
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 1
findings_total: 1
material_change: false
---

# Verification — unit 126

## Per-finding outcomes

### F1 — insufficient_verification (boundary-only positivity sampling missed interior negatives)

**Classification:** resolved

**What changed:**

- SymPy `scripts/moving_throat_pde_stage126_positive_source_families_sympy_audit.py:55-77`: the old `min_xi0` / `val_xi1` / `sigma_corner` boundary-and-corner block (which sampled only `xi=0` min, the `xi=1` value, and the `(z=L, xi=0)` corner) was replaced by a two-piece global-positivity argument:
  - (i) exact affine-in-xi test `d2_xi = sp.diff(sigma_xi, xi, 2)`, asserting `d2_xi == 0` (line 61-64);
  - (ii) both endpoint slices' minima over the FULL interval `Interval(0, L)`: `min_xi0` and `min_xi1`, each asserted `>= 0` is `sp.true` (lines 65-72), with trailing identity pins `min_xi0 == 0` and `sigma_xi(xi=1) == 1/L` retained (lines 73-77).
- Mathematica `mathematica/moving_throat_pde_stage126_positive_source_families_mathematica_audit.wl:64-83`: mirror replacement — `expectZero["sigma_xi affine in xi", d2Xi]` (affine test), the `sigma_match(lM)=0` boundary-min and `sigma_xi(xi=1)=1/lM` value checks, plus a whole-box `globalPos = Resolve[ForAll[{z,xi}, 0<=z<=1 && 0<=xi<=1, (sigmaXi /. lM->1) >= 0], Reals]` with a hard `fail` if not `TrueQ`. The redundant `sigma_match(0)=k` and `sigma_xi(z=lM,xi=0)=0` corner samples were dropped (subsumed by the affine + ForAll checks), as the directive permitted.

**Assessment:**

The edit matches the directive's "required change" for both engines, and the soundness of the argument holds:

- **Affineness (i) is genuine and the right structural fact.** `sigma_xi = (1-xi)*k*cos(k*z) + xi/L = k*cos(k*z) + xi*(1/L - k*cos(k*z))` is degree-1 in xi, so `d^2/dxi^2 = 0` exactly. SymPy reports `d^2 sigma_xi / d xi^2 = 0`; Mathematica reports `PASS: sigma_xi affine in xi`. A function affine in xi attains its min over `xi in [0,1]` at an endpoint — so checking both endpoint slices is sufficient over the full box.
- **Endpoint slices (ii) are genuinely nonneg on [0,L].** The `xi=0` slice is `k*cos(k*z)` with `k=Pi/(2L)`, so `k*z in [0, Pi/2]` and `cos >= 0`, min 0 at z=L (output: `min sigma_xi(xi=0) on [0,L] = 0`). The `xi=1` slice is the constant `1/L`, nonneg for `L>0` (output: `min sigma_xi(xi=1) on [0,L] = 1/L`). The minima are taken over `Interval(0, L)` / quantified over the full z-range, not at sampled corners. Together (i)+(ii) imply `sigma_xi >= 0` on the entire box.
- **Non-tautology / can-it-fail confirmed.** The affineness check is exact symbolic: the canonical interior-dip attack `sigma_xi + (-A)*xi*(1-xi)*sin(2*pi*z/L)/L` (which preserves both endpoint slices and the corner) has `d^2/dxi^2 = 2*A*sin(2*pi*z/L)/L != 0`, so it FAILS (i) — the exact gap the original corner-only block missed. A sign-changing endpoint slice would make its interval minimum negative, failing the `>= 0` / `ForAll` test (ii). Neither check is a forced identity.
- **lM->1 deviation is legitimate.** The recorded Mathematica deviation substitutes `lM -> 1` for the positivity quantifier because `Resolve` with free positive `lM` did not collapse to `True`. The family is genuinely scale-covariant in lM for positivity: with `z = lM*u`, `k*z = (Pi/2)*u` is lM-independent, and `sigmaXi = (1/lM)*[(1-xi)*(Pi/2)*cos((Pi/2)*u) + xi]`. The bracket equals `lM * (sigmaXi|_{lM=1})` at `z=u`, and the prefactor `1/lM > 0`, so positivity at `lM=1` on `z in [0,1]` implies positivity for all `lM>0` on `z in [0,lM]`. The deviation is sound; the `0<=z<=1` domain used at `lM=1` is the correct rescaled box. (The trailing `Simplify[..., lM>0]` is a harmless no-op since `globalPos` already resolved to literal `True`.)

No collateral edits beyond the F1 block: the diff touches only the named line ranges in each file.

## Exec log assessment

**SymPy:** exit=0. Notable lines:
- `d^2 sigma_xi / d xi^2 = 0`
- `min sigma_xi(xi=0) on [0,L] = 0`
- `min sigma_xi(xi=1) on [0,L] = 1/L`
- `Check 2/pi < g_- < pi/4 -> True`

**Mathematica:** exit=0. Notable lines:
- `PASS: sigma_xi affine in xi`
- `PASS: sigma_match(lM) = 0 (boundary min on [0,lM])`
- `PASS: sigma_xi(xi=1) = 1/lM`
- `ForAll sigma_xi >= 0 on box -> True`
- `Stage 126 Mathematica audit passed.`

**PASS-line tally (per the directive's expected lines):**
- SymPy — all three new lines present and as expected: `d^2 sigma_xi / d xi^2 = 0`, `min sigma_xi(xi=0) on [0,L] = 0`, `min sigma_xi(xi=1) on [0,L] = 1/L`; prior lines (`min sigma_match on [0,L] = 0`, `xi_* = ...`, `g_xi(xi_*) - g_- = 0`, `Check 2/pi < g_- < pi/4 -> True`) unchanged; exit 0.
- Mathematica — both new PASS lines present: `PASS: sigma_xi affine in xi` and `ForAll sigma_xi >= 0 on box -> True`; plus `PASS: sigma_xi(xi=1) = 1/lM` and `PASS: sigma_match(lM) = 0`; final banner `Stage 126 Mathematica audit passed.`; exit 0.

Passing logs are necessary but not sufficient; the substance (affineness exactness + full-interval/box endpoint nonnegativity + scale-covariance) is confirmed above. The checks are falsifiable, not tautological.

**Output freshness:** confirmed. Script mtimes are 2026-05-29 14:31:20 (SymPy `.py`) and 14:32:01 (Mathematica `.wl`); both committed `.txt` outputs are stamped 2026-05-29 15:27:13 — newer than the scripts, so the saved transcripts were regenerated post-fix.

## Material-change assessment

`material_change`: false.

This is a verification-surface strengthening only. The replaced block did not change any derived result: the source family `sigma_xi`, its normalization, the bias `g_xi`, `g_-`, `xi_*`, the traction ratio, and the `2/pi < g_- < pi/4` bracket are all untouched and still report identical values (output lines for `g_-^F1`, `xi_* = 0.18391840551153962831`, `g_xi(xi_*) - g_- = 0`, interval `True` are unchanged). The new checks only prove the positivity claim more rigorously. No downstream unit depends on the old corner-sampling internals, so no narrow re-audit is warranted.

## Reconcile-note rows (PASS rows confirmed untouched)

- **Row 1 (min_match):** `scripts/...sympy_audit.py:51-54` `min sigma_match on [0,L] = 0` with `min_match == 0` assertion — present, correct, untouched. Output line `min sigma_match on [0,L] = 0`. PASS.
- **Row 4 (interval bracket):** SymPy `:88-94` `bool(2/pi < g_- < pi/4)` with raise-on-fail, and Mathematica `:90-94` analogue with `fail` — present, correct, untouched. Both output `Check 2/pi < g_- < pi/4 -> True`. PASS.
- **Row 5 (banner):** `STAGE 126 — EXPLICIT POSITIVE SOURCE FAMILIES` at `:13` (SymPy) and `:26` (Mathematica) — correct Stage 126 label, untouched. PASS.

All three are outside the diff hunks and unchanged.

## Side observations (non-blocking)

- The SymPy nonnegativity assertions use the directive's primary form `(sp.simplify(min_xiK) >= 0) is sp.true` rather than the `bool(...)` fallback; given `L` is declared positive and the minima are `0` and `1/L`, this evaluates to `sp.true` cleanly, so no weakening to a tautology occurred. Consistent with the directive's note.
- The Mathematica box check uses the rescaled domain `0<=z<=1` at `lM->1` (not `0<=z<=lM`), which is the correct box after the scale substitution; the printed banner label still references `[0,lM]`, which is accurate via the scale-covariance argument. Cosmetic only.

## Verdict justification

The single finding F1 is fully resolved in both engines. The boundary-and-corner sampling that left an interior gap was replaced by a sound, falsifiable two-piece global-positivity argument (exact affine-in-xi structural test plus full-interval/whole-box endpoint nonnegativity), and the directive-approved `lM->1` scale-covariance deviation is mathematically valid. Both committed outputs were regenerated post-fix and exit 0 with all expected PASS lines; the reconcile-note PASS rows are untouched and correct; and no derived result changed (material_change=false). Verdict: verified.
