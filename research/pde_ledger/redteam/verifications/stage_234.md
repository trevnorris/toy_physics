---
unit_id: 234
batch: VII.2
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-06-02T22:30:00-06:00
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 2
findings_total: 2
material_change: false
---

# Verification — unit 234

## Per-finding outcomes

### F1 — missing_verification_script (missing_mathematica)

**Classification:** resolved

**What changed:**
New file `mathematica/moving_throat_pde_stage234_direct_branch_observable_static_gate_and_the_two_observable_kill_test_mathematica_audit.wl` (created; mtime 22:07:22, newer than nothing pre-existing). Implements M1–M7 as banner-grouped, independently-checked blocks with `expectZero`/`expectLogZero`/`expectTrue` helpers that call `fail[...] -> Exit[1]` on any nonzero residual or non-True condition.

**Assessment:**
Genuinely independent route, NOT a transliteration of the `.py`:
- M2 linearization: `linCoeff[expr] = FullSimplify[Coefficient[Normal[Series[expr,{eps,0,1}]], eps, 1]]` (wl:66-69, 116-118) — the directive's mandated `Series`/`Coefficient` route, NOT the SymPy `diff(expr,eps).subs(eps,0)` choreography.
- M7 balanced family: `Minimize[{R1^2+E1^2, -R1-cEta*E1==xi},{R1,E1},Reals]` (wl:187-190) — NOT the SymPy explicit Lagrange-multiplier `sp.solve` with `lam`. Yields the same closed form `(-xi/(1+cEta^2), -cEta*xi/(1+cEta^2))` after stripping the `Piecewise[..., epsEtaStar!=1]` wrapper.
- M1 roundtrip via `PowerExpand`/`FullSimplify` native log handling (wl:40-45, 98-100) vs SymPy `expand_log(...,force=True)`.
- M4 cancellation derived natively: `xi1Raw = qNt1 + (bStar/cStar)*qTr1` built from M2's Series coefficients, then BOTH `D[xi1Raw,r1]==0` and `FreeQ[FullSimplify[xi1Raw], r1]` (wl:146-148) — does NOT reuse a pre-simplified `-R1-c_eta E1`, so the cancellation is load-bearing (contingent on the exact bStar/cStar weight).
- camelCase symbols (`qTr/qNt/cStar/rtr`), banner-grouped M1–M7 ordering, distinct helper architecture — not a line-for-line port.
All M1–M7 covered and PASS in the exec log (exit 0). The log-chart roundtrip (M1), `D[Xi1,r1]` cancellation (M4), and balanced-family minimization (M7) are all derived natively. `ConditionalExpression`/`Piecewise` wrappers stripped before the zero test per the idiom requirement. No collateral edits (only the new `.wl`).

### F2 — tautological_check

**Classification:** resolved

**What changed:**
`scripts/...stage234..._sympy_audit.py:148-163`. Added three value-sensitive guards after the literal definitions (py:148-154): `0 < float(robust) < float(nonempty)`, `simplify(nonempty - robust) != 0`, and documented-gap `abs(float(nonempty)-float(robust)-0.369688735168111) < 1e-12`. Removed the four tautological endpoint `assert_zero(Xi1.subs(R1, edge) ∓ const)` identities and replaced them with two half-width identities `assert_zero(simplify(robust_upper-robust_lower) - 2*robust/Abs(vareps))` and the nonempty analogue (py:162-163). Prints retained.

**Assessment:**
The original four endpoint checks were tautological exactly as the auditor described: the carried ceiling literal cancels by construction for any value (`-(-c_eta*E1 - const) - c_eta*E1 = const`), so a swap or digit-drop passed green. The fix makes a wrong value genuinely trip:
- Swap robust<->nonempty: the ordering guard `0 < 0.7376 < 0.3679` is False -> AssertionError.
- Digit-drop/corruption of either literal: the gap no longer equals `0.369688735168111` -> AssertionError (1e-12 tolerance).
- Equal values: distinctness guard fails.
The two retained half-width identities remain value-tautological in isolation (the literal appears on both sides via the same symbol), but the directive explicitly acknowledged this and paired them with guards (1)+(3) to compensate; the combined check is non-vacuous. The `.wl` M6 mirrors the ordering guard (`expectTrue["M6 0 < robust < nonempty"]`, wl:177) and exercises all four edges. No new constant or paper claim introduced; the gap literal is just the arithmetic difference of two already-stated ceilings.

## Exec log assessment

**SymPy:** exit=0. Notable lines: "Robust strip: ... 0.367930328492646/Abs(vareps) ...", "Nonempty strip: ... 0.737619063660757/Abs(vareps) ...", "Balanced minimal norm: (R1, E1) = (-xi/(eps_eta_star**2/(1 - eps_eta_star)**2 + 1), ...)", "All Stage 234 symbolic checks passed." No failures/warnings.

**Mathematica:** exit=0. 25 explicit PASS lines (M1×3, M2×3, M3×4, M4×2, M5×3, M6×5, M7×5), zero FAIL lines, no parser-skip. Notable: "PASS: M4 D[Xi1, r1]", "PASS: M4 simplified Xi1 contains no r1", "PASS: M6 0 < robust < nonempty", "PASS: M7 balanced R1 minimizer", "All Stage 234 Mathematica checks passed." Balanced Minimize emitted the `Piecewise[..., epsEtaStar!=1]` form, correctly stripped before the equality test.

**Output freshness:** confirmed. sympy_output mtime 22:10:48 > .py 22:05:42; mathematica_output mtime 22:10:48 > .wl 22:07:22. Both regenerated post-fix. (MANIFEST still carries a stale pre-fix sympy_output mtime 1778525515 and null mathematica_output path — a manifest-bookkeeping lag, not a freshness problem; the on-disk outputs are fresh. Orchestrator should refresh the manifest file block on advance.)

## Material-change assessment

`material_change`: false. F2 only strengthened guard assertions on already-carried import constants (no derived value changed). F1 added a second engine that reproduces the same closed forms the SymPy engine already produced — no result changed, it is a corroborating cross-check. No downstream unit's inputs are altered.

## Side observations (non-blocking)

- MANIFEST 234 `files.mathematica_output.path` is `null`/`exists: false` and `sympy_output.mtime` is pre-fix; the orchestrator's normal advance should re-stamp these. Cosmetic, does not affect the verdict.
- Engine cross-check now genuine: both engines independently reproduce `Xi1 = -R1 - c_eta*E1`, `d Xi1/d(ln R_tr) = 0`, the rigid-mouth `q_nt = Xi1`, and the balanced minimizer — the F1 gap is fully closed.

## Verdict justification

Both findings are resolved. F1's new `.wl` is a genuinely independent route (Series/Coefficient linearization and Minimize balanced family, distinct symbol set and architecture, native log/cancellation derivation), covers M1–M7, and passes with a non-vacuous `Exit[1]` fail path. F2's tautology is closed by value/ordering/gap guards that demonstrably trip on swap or corruption, while retaining the strip rearrangement print. Both exec logs exit 0 with PASS counts matching the in-file check inventory and no parser-skips; outputs are fresh. No regressions in the diff, no material change to any derived result.
