---
unit_id: 233
batch: VII.2
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-06-02T22:30:00-06:00
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 4
findings_total: 4
material_change: false
---

# Verification — unit 233

## Per-finding outcomes

### F1 — tautological_check

**Classification:** resolved

**What changed:**
`scripts/...stage233...sympy_audit.py:118-121`. The two critical-pressure
definitions were changed from the circular round-trip
`Pcrit_robust = Pbar_num*(1+robust_budget)` / `Pcrit_nonempty = Pbar_num*(1+nonempty_budget)`
to independently-sourced literals:
```
Pcrit_robust   = sp.Float("0.0028313316855593175")
Pcrit_nonempty = sp.Float("0.0035965105896846573")
```
with an inline comment citing the source (Stage 224 `ceilings` dict). The
recovery then computes `budget = Pcrit/Pbar_num - 1` (lines 123-124) and asserts
it equals the carried budgets at tol 1e-15 (lines 126-127).

**Assessment:**
Non-circular, confirmed two independent ways:

1. **Source verification.** Both Pcrit literals appear verbatim in the cited
   independent source `scripts/...stage224...sympy_audit.py:137-138`
   (`both_10 = 0.0028313316855593175`, `one_10 = 0.0035965105896846573`), which
   carries them from Stage 223's dynamic-survival-window walls
   (`stage223...:415-416`, independently computed). In Stage 224 the budget is the
   *output* of Pcrit (`eps_xi_budget = Pcrit_val/barP0_compat - 1`, line 145), so
   Stage 233 now reuses Pcrit as a genuine upstream input rather than back-deriving
   it from the budget. The directive's preferred route (independent `P_crit` per
   branch, two distinct literals) was taken exactly.

2. **Residual signature.** A circular round-trip yields a residual of *exactly* 0
   (`Pbar*(1+b)/Pbar - 1 == b` identically). The refreshed outputs instead show
   tiny-but-nonzero residuals: robust `-4.35e-18`, nonempty `-3.24e-16` (the
   Mathematica log prints these as exact rationals
   `-900068804371/206979231806288500000000000000` etc.). A nonzero residual is
   only possible if Pcrit was sourced as a separately-truncated decimal — proving
   the check is no longer `b == b`. The two budgets recover from two distinct,
   independently-carried Pcrit values, not from one re-derived P_crit.

**Non-vacuous:** yes. `assert_close` raises on `abs(actual-expected) > 1e-15`. The
nonempty residual (3.2e-16) sits within but the same order as the tolerance, so a
single transcribed digit in either Pcrit literal shifts the residual well above
1e-15 and trips the assertion. The directive's perturbation self-test is satisfied.

### F2 — missing_verification_script (missing_mathematica)

**Classification:** resolved

**What changed:**
New file `mathematica/...stage233...mathematica_audit.wl` (208-line independent
audit), output regenerated, exit 0.

**Assessment:**
Genuinely independent, not a transliteration:
- **Different decomposition for M1/M2.** Where the `.py` proves the track lock by
  bare zero-substitution (`Theta1.subs(dln_Rtr:0)`), the `.wl` extracts
  `Coefficient[theta1, dlnRtr]` and `Coefficient[theta1, dlnRtr, 0]` and asserts
  the coefficient (`-1`/`+bstar`) and constant-term structure, then *derives* the
  lock from the coefficient data (`thetaConstant + thetaCoeff*0`). This is the
  coefficient-extraction route the directive (M1/M2) demanded.
- **Different M6 route.** The `.py` re-substitutes `Delta_norm:0` into `gate_rhs`
  and re-asserts the same hand form (the F3 weak check). The `.wl` routes through
  the Pbar form: `calibratedFromPbar = pbarGate /. deltaNorm->0` then FullSimplify
  against `pcrit*mhat0^2/tQuad - 1`.
- **Different infrastructure / vars.** Native primitives (`Coefficient`,
  `Together`, `FullSimplify`, `Cancel`), distinct var names
  (`dlnRtr`/`bstar`/`pcrit`/`tQuad`), exact-rational arithmetic via
  `exactDecimal`/`Rationalize`, and `expectSmall` tolerance testing for M8 rather
  than float `assert_close`. `cleanExpr` strips `ConditionalExpression[0,...]` per
  the idiom note.

**Non-vacuous:** yes. `fail[]` calls `Exit[1]`; `expectZero` fires fail unless
`TrueQ[res === 0]`; `expectSmall` fires fail unless `Abs[res] <= tol` holds under
assumptions. M8 residuals are nonzero and bounded-tested, not `True===True`.

### F3 — insufficient_verification (folded into F2)

**Classification:** resolved

**What changed:** No standalone `.py` edit (as directive specified). Addressed via
the F2 `.wl`'s M2 coefficient-extraction (replacing the trivial `Theta1` zero-sub)
and M6 route-through-Pbar (replacing the re-substitute-and-re-assert). Both
independence requirements are honored in the `.wl` as read.

**Assessment:** Correct. The directive folded F3 entirely into F2's M2/M6
independence; those routes are present and exercised (PASS: M2 ×4, PASS: M6).

### F4 — stale upstream stage labels

**Classification:** resolved

**What changed:** Four comment/print strings in
`scripts/...stage233...sympy_audit.py` updated: line 20/32 `Stage 188`→`Stage 239`
(compiler), line 113 `Stage 224`→`Stage 241` (budgets), line 129 (was 128)
`Stage 223 / Stage 224`→`Stage 240 / Stage 241`.

**Assessment:** Diff confirms these are pure string edits — no symbol, value, or
assertion changed. Confined to this one file; no blind project-wide renumber.
Mapping matches notes §9 (239 compiler / 240 compatibility point / 241
budgets-ceiling). The output reflects the new labels and still exits 0.

## Exec log assessment

**SymPy:** exit=0. Notable lines:
- `Implied Pcrit (robust) = 0.0028313316855593175` / `(nonempty) = 0.0035965105896846573` — the independent literals, no longer back-computed.
- `Recovered robust budget = 0.367930328492646003049615366365` vs carried `0.367930328492646` — matches to 1e-15.
- `All Stage 233 symbolic and numerical checks passed.`

**Mathematica:** exit=0. Notable lines:
- `M2 Theta1 track lock from coefficient data residual = 0` / `PASS` (coefficient route, not substitution).
- `M8 robust numeric budget residual = -900068804371/206979231806288500000000000000` / `PASS` (nonzero residual under tol — confirms genuine recovery, not tautology).
- `All Stage 233 Mathematica checks passed.`

**Pass counts:** 14 `PASS:` lines in the `.wl` output == 14 `expect{Zero,Small}` invocations in the `.wl` (16 matches minus the 2 `:=` definitions). No FAIL/Error anywhere. Manifest coverage: M1×2, M2×4, M3×1, M4×1, M5×1, M6×1, M7×1, M8×3 — all M1-M8 present, no silent parser-skip.

**Output freshness:** confirmed. `.py` mtime 22:08:26 < `.txt` 22:10:24; `.wl` mtime 22:06:28 < `.txt` 22:10:24. Both outputs regenerated post-fix.

## Material-change assessment

`material_change`: false.

No derived result changed. The F1 fix replaced a tautological check with a genuine
recovery of the *same* carried budgets (`0.367930...`, `0.737619...`) from the
*same* independently-carried Pcrit values that already live in Stage 224; no number
was altered, only the verification direction. F2 adds a second engine confirming
the existing identities. F4 is comment/label-only. Downstream units depending on
Stage 233's outputs see identical values; no upstream-stale propagation warranted on
substance.

## Side observations (non-blocking)

- Scope was clean: only the four expected paths changed (`.py`, its `.txt`, new
  `.wl`, new `.wl` `.txt`). No `paper.tex`/`notes/` touched, no out-of-scope source
  edits, no collateral renumbering beyond the four F4 strings.
- Ansatz/provenance note (informational): the Pcrit ceilings `0.0028313...` /
  `0.0035965...` trace to Stage 223's compatibility walls; they are computed, not
  postulated, so no new ansatz-catalog entry is implied by this stage.

## Verdict justification

All four findings are `resolved`. F1's circular round-trip is genuinely broken: the
two budgets now recover from two distinct, independently-sourced critical-pressure
literals that exist verbatim in the cited Stage 224 source (and Stage 223 walls),
and the tiny-but-nonzero recovery residuals (-4.3e-18, -3.2e-16) are the positive
proof that the check is no longer `b == b`. The new `.wl` is independent (coefficient
extraction for M1/M2, route-through-Pbar for M6, native primitives and distinct
choreography), with non-vacuous `Exit[1]` fail paths and a PASS count (14) matching
its manifest invocations (M1-M8). F3 is satisfied through F2's independence routes,
and F4 is a correct, scoped label fix. Both engines exit 0 with fresh outputs and no
regressions in the diff.
