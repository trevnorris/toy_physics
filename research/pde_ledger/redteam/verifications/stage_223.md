---
unit_id: 223
batch: VII.1
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-06-02T16:50:00-06:00
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 3
findings_total: 3
material_change: false
---

# Verification — unit 223

## Per-finding outcomes

### F1 — paper_misalignment (value_mismatch, user-resolved)

**Classification:** resolved

**What changed:**
Notes-only edit (per the user-authorized `## RESOLVED` block). In
`notes/stages/moving_throat_pde_stage223_..._sympy_audit.md:351` the λ_W=0.2
scan row now reads `| 0.2 | 0.000576970879843 | 29.3158464872314 | 138.814136942081 | 137.502546600713 |`.
Grep confirms `206.814136942081` and `205.502546600713` no longer appear anywhere
in the notes, and the two corrected values are present. The first two columns
(P0, K) were already correct and are unchanged; only the +68 transcription slip on
the two wall R_Q cells was repaired (fractional tails `...814136942081` /
`...502546600713` preserved). No script change was made for F1, as directed.

**Assessment:**
Correct and complete. The corrected values match the SymPy `.txt` output line 60
(`138.81413694208146`, `137.50254660071312`), which is the authoritative
re-derived source. The directive's acceptance criterion is met exactly. Per finding
#5 in the verify brief: the notes file is internally titled "Stage 240" (pre-renumber)
— this is incidental; the F1 fix targeted only the line-351 data cells and neither
depended on nor altered the title. The values, not the title, are what F1 concerned.
The new `.wl` (M9) independently recomputes the scan, but only for λ_W = 0.4..1.0
(the λ_W=0.2 row was correctly held out of the .wl pending F1 per the directive note);
the λ_W=0.4 row's wall R_Q (30.20 / 36.17) is reproduced by both engines, so the
scan formula that produced 138.8/137.5 is corroborated as the correct family.

### F2 — tautological_check

**Classification:** resolved

**What changed:**
SymPy `py:156` inserts `P0_target_compat_solved = sp.solve(sp.Eq(K_norm, K_pole), P0_target)[0]`,
and `py:160` replaces the old circular residual with
`compatibility_identity = sp.cancel(P0_target_compat_solved - P0_target_compat)`,
asserted `== 0` at py:161. `K_norm` (py:155, `N0/P0_target + B0 + Z0`) and `K_pole`
(py:154, `3(M+B2+Z2)²/(B4+Z4) + B0 + Z0`) are kept as the two independently-built
stiffnesses. The boxed closed form `P0_target_compat` (py:157) is retained as a
separate comparison symbol. Lines 162 and 164-167 left as-is per directive Steps 5-6.
Diff is exactly the three-line change described; no collateral edits.

**Assessment:**
De-tautology is GENUINE, not still-circular. The compatibility target is now an
INDEPENDENT result of `sp.solve` (inverting the linear-in-1/P0_target equation
`K_norm == K_pole`), compared against a separately-authored boxed form. If the boxed
form were wrong — wrong constant 3, B4↔B2 swap, missing (B4+Z4) factor — the solve
would still return the true inverse `N0(B4+Z4)/(3(M+B2+Z2)²)` and the `cancel`
residual would be nonzero, failing the assert. It can no longer reduce to
`a/(a*b/c)==c/b`. The output (`.txt` line 23) shows `P0_target_compat` as the full
nontrivial rational expression in all couplings, and the run exits 0 with the assert
passing. (The primitive-specialization assert at py:164-167 remains a benign
definition-inline, left untouched per directive Step 6 — it was not the load-bearing
circularity F2 targeted.)

### F3 — missing_verification_script (missing_mathematica)

**Classification:** resolved

**What changed:**
New independent `.wl` at the exact target path
`mathematica/moving_throat_pde_stage223_..._mathematica_audit.wl`, covering
claim manifest M1-M9. Mathematica exec log exits 0 with 23 PASS, 0 FAIL.

**Assessment:**
Independent re-derivation, NOT a transliteration. Concrete evidence of distinct route:
- M3 derives the surface via `Solve[kNormIso == kPoleIso, p0Target]` (the F2-corrected
  route) and compares to a separately-written `boxedTarget`, plus an equivalent
  N0/P0 form — non-tautological, same logic as the corrected SymPy.
- M5 builds the quartic from the determinant/product form and reads degree via
  `Exponent[fyExpanded, y]` and `CoefficientList` length (= 5), and additionally
  checks the leading (y^4) coefficient equals `mass` — a DIFFERENT decomposition than
  SymPy's `Poly(...).degree()`.
- M7 selects wall-like poles via `Sort[TakeLargest[roots, 2]]` (native ordering) vs.
  SymPy's `roots[-2:]` slice convention, and independently checks
  `Min[wallRoots] > Max[internalRoots]` — exactly the anti-transliteration guard the
  directive asked for. Root residuals are O(1e-15) against the polynomial.
- M2 reconstructs `u4 - 4 u2^2` from native `d0Iso/d2Iso/d4Iso` definitions; residual 0.
- M1 uses native `Integrate` for κ = 2√2/π; residual 0.
- M4 specializes via FullSimplify substitution; residual 0.
- M6/M8/M9 numeric values match the SymPy `.txt` to ~1e-13..1e-16 (M6 P0
  residual 2.9e-21, K 1.4e-16; M8 four R_Q within 2e-13; M9 thresholds and the four
  survival-window ceilings within ~6e-17). M9 scan correctly limited to λ_W=0.4..1.0.
All M1-M9 PASS; no fabricated literals — the sample-slice rules (λ_B=1/2, λ_U=3/10,
λ_W=2/5, λ_R=1/4, Ω_U=1, Ω_W=7/5, ϖ=2, M=1, a=c_s=1) match upstream Stage 222 and the
SymPy script; census/threshold/window literals are the engine-recomputed values, not
hand-pasted constants divorced from derivation.

## Exec log assessment

**SymPy:** exit=0. Notable lines:
- `kappa = 2*sqrt(2)/pi` (M1/A1 holds).
- `P0_target_compat = kap**2*varpi**2*(Omega_U**2*lam_W + lam_R*lam_U)**2*(...)` — the
  full nontrivial solved/boxed expression (confirms F2 target is not a trivial ratio).
- Scan row `(0.2, ..., 138.81413694208146, 137.50254660071312)` (matches corrected notes).
- `Stage 223 audit completed successfully.` exit_code 0.

**Mathematica:** exit=0. Notable lines:
- `PASS: M3 solved target minus boxed target` and `PASS: M3 equivalent N0/P0
  compatibility equation` (independent Solve-based surface).
- `PASS: M5 determinant/product polynomial is quartic` (degree 4, 5 coefficients) and
  `PASS: M5 leading quartic coefficient is mass`.
- `PASS: M7 wall-like roots are the largest two roots = True` (TakeLargest ordering).
- `All Stage 223 Mathematica checks M1-M9 completed successfully.` exit_code 0.
- 23 PASS lines total, 0 FAIL — matches the orchestrator's reported count.

**Output freshness:** Confirmed. SymPy `.txt` mtime 16:37:47 > script 16:23:02;
Mathematica `.txt` mtime 16:37 > `.wl` 16:26:18. Both outputs regenerated post-fix.

## Material-change assessment

`material_change`: **false**.

F1 is a notes-only typo correction that aligns the prose to the already-verified
script value (no derived quantity changed). F2 changed only HOW the compatibility
surface is verified (solve-and-compare vs. circular rewrite); `P0_target_compat`
keeps the identical symbolic value and every downstream `assert_close` literal is
unchanged and still passes. F3 added a corroborating second engine that reproduces
all existing results. No derived result that a downstream unit (>223) could depend on
was altered. Downstream units do not need re-audit on account of stage 223.

## Side observations (non-blocking)

- The notes file remains internally titled "Stage 240" (and references "Stage 240"
  in §7/§8 prose). This is the known pre-renumber artifact and is out of scripts-only
  verifier scope; flagging only — it is not a finding and does not affect verification.
- The primitive-specialization assert (SymPy py:164-167) is still a definition-inline
  tautology by construction, but it was explicitly left as-is per directive Step 6 and
  was never the load-bearing circularity; M4 in the `.wl` exercises the same claim via
  an independent FullSimplify, so the specialization is now cross-engine corroborated.

## Verdict justification

All three findings are resolved. F1's notes typo is corrected to the authoritative
script values with the old 206/205 figures fully removed. F2's de-tautology is genuine:
the compatibility target now comes from an independent `sp.solve` of the
stiffness-equality condition and is compared against the separately-written boxed
form, so the assert would fail if the surface formula were wrong. F3's new Mathematica
audit is a real independent re-derivation (Solve-based surface, determinant/product
quartic with Exponent/CoefficientList, TakeLargest wall-pole ordering) covering M1-M9,
exits 0 with 23 PASS / 0 FAIL, and its numerics agree with SymPy to ~1e-13 or better.
Both exec logs pass, outputs are fresh, the diff shows no regressions, and no derived
result changed (material_change false). Verdict: verified.
