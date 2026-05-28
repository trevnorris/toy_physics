---
unit_id: 130
batch: IV.4
verifier_model: claude-opus-4-7[1m]
verify_date: 2026-05-27T00:00:00Z
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 2
findings_total: 2
material_change: false
---

# Verification — unit 130

## Per-finding outcomes

### F1 — insufficient_verification (strict monotonicity sweep)

**Classification:** resolved

**What changed:**
- `scripts/moving_throat_pde_stage130_mouth_bias_map_sympy_audit.py:38-44` — new 6-point sweep over `Pi ∈ {1/10, 1/2, 1, 15088/10000, 3, 10}`. Computes `sp.N(dgPi.subs(Pi, val), 30)`, prints each value, and `raise AssertionError` if `deriv_val <= 0`.
- `mathematica/moving_throat_pde_stage130_mouth_bias_map_mathematica_audit.wl:61-70` — new `Module` over the same 6 values; computes `N[dgPi /. piM -> v, 40]`, prints each, and routes to `pass[]`/`fail[]` via `If[TrueQ[dv > 0], …]`.

**Assessment:**
Insertion points match the directive exactly (after sympy line 33, after WL line 57). Both implementations use the same 6 positive Π values bracketing Π_* ≈ 1.5088. The assertions are non-tautological: each evaluates the closed-form derivative numerically and rejects on `≤ 0`. Exec logs confirm all six points yielded strictly positive values (monotone decreasing from 0.0865 at Π=0.1 down to 0.00465 at Π=10, all > 0). The sympy log shows six `dg/dPi at Pi=...` lines with no `AssertionError`; the mathematica log shows six `PASS: dg/dpiM > 0 at piM=...` lines. No collateral edits.

### F2 — insufficient_verification (closed-form match)

**Classification:** resolved

**What changed:**
- `scripts/moving_throat_pde_stage130_mouth_bias_map_sympy_audit.py:16-18` — adds `gPi_boxed = 2*Pi*(2*Pi*sp.exp(Pi) + sp.pi) / ((4*Pi**2 + sp.pi**2) * (sp.exp(Pi) - 1))` and `raise AssertionError` if `sp.simplify(gPi - gPi_boxed) != 0`.
- `mathematica/moving_throat_pde_stage130_mouth_bias_map_mathematica_audit.wl:43-44` — adds `gPiBoxed = …` and `expectZero["g_Pi matches paper boxed closed form", gPi - gPiBoxed]`.

**Assessment:**
Insertion points match the directive (after sympy line 15, after WL line 42). The boxed expressions are written in terms of the integral-evaluated `gPi` so the comparison is non-tautological — `gPi` is produced by `sp.integrate`/`Integrate` of `sigma*f` and the boxed form is the literal expression from the notes. SymPy passed silently (no AssertionError, residual implicitly = 0). Mathematica log line 5 prints `g_Pi matches paper boxed closed form = 0` followed by `PASS:` line 6. Both engines agree the symbolic difference simplifies to zero. No collateral edits beyond the two-line inserts.

## Exec log assessment

**SymPy:** exit=0. Notable lines:
- `g_Pi = 2*Pi*(2*Pi*exp(Pi) + pi)/((4*Pi**2 + pi**2)*(exp(Pi) - 1))` (closed form unchanged)
- `Covariance identity residual = 0`
- Six `dg/dPi at Pi=... = <positive value>` lines, all positive
- `Pi_* = 1.50882951349315861144664988336` (root-find unchanged)

**Mathematica:** exit=0. Notable lines:
- `PASS: g_Pi matches paper boxed closed form` (F2)
- `PASS: covariance identity`
- `PASS: uniform-source limit` and `PASS: point-source limit`
- Six `PASS: dg/dpiM > 0 at piM=...` lines (F1) — note that `Print[..., fmt[15088/10000], ...]` causes a rational like `943/625` to render on multiple lines, but each PASS line is intact
- `PASS: Family-1 compensation point` (residual ~8.9e-58)
- Trailer: `Stage 130 Mathematica audit passed.`

**Output freshness:** Confirmed. Script mtimes are 17:23; sympy output is 17:45 and mathematica output is 17:47 — outputs were re-generated post-fix.

## Material-change assessment

`material_change`: false.

Both new assertions are pure checks on already-derived expressions; they add coverage but do not change any constant, expression, or numerical result that downstream units consume. `gPi`, `Π_*`, and `g'(Π_*)` printed in the post-fix output are identical to the pre-fix values (the closed-form derivative sweep evaluates the same `gPi` that downstream stages depend on). No downstream propagation required.

## Side observations (non-blocking)

1. The Mathematica `Print["dg/dpiM at piM=", fmt[v], …]` statement on line 66 renders `1/2`, `15088/10000` and `1/10` as fractional rationals across multiple console lines, producing slightly noisy log layout (e.g., the `PASS` line for `piM=15088/10000` reads `PASS: dg/dpiM > 0 at piM=943` then `---` then `625` on separate lines because `ToString[15088/10000]` collapses to `943/625`). Each PASS still emits as a single logical `Print` so the assertion fires correctly, but the visual readout is split. Cosmetic only.
2. Cluster A renumbering from the IV.4 paper-alignment resolutions is in place: `.wl` line 32 banner now reads `STAGE 130` (not the prior `STAGE 113`); notes H1 was not in scope of this script audit and is tracked separately by the resolutions doc.
3. Sympy script does not print an explicit `PASS:` confirmation for either new check (silent success on no-`AssertionError`), while Mathematica does. Both behaviors are consistent with the directive — sympy uses `raise`, mathematica uses `expectZero`. Not a defect.

## Verdict justification

Both findings from the original report (F1 strict-sign sweep, F2 boxed-form equality) were applied verbatim to both engines at the directive-specified line ranges. Exec logs from the post-fix runs (newer mtimes than the scripts) show both engines exit 0 with all six monotonicity points strictly positive and the closed-form residual identically zero. The new assertions are non-tautological — they compare the integral-derived `gPi`/derivative against an independent literal expression of the boxed form. No regressions; no collateral edits. Verdict: `verified`.
