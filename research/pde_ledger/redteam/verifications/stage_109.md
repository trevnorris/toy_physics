---
unit_id: 109
batch: IV.2
verifier_model: claude-opus-4-7[1m]
verify_date: 2026-05-27T00:00:00Z
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 3
findings_total: 3
material_change: false
---

# Verification — unit 109

## Per-finding outcomes

### F1 — tautological_check

**Classification:** resolved

**What changed:**
`scripts/moving_throat_pde_stage109_linearized_branch_selection_sympy_audit.py:53-54` adds the new anchored check between the existing `print('a5 preservation condition = ...')` (L52) and the existing tautological `expect_zero('preservation substitution', ...)` (L55):

```
expected_a5_sol = -sp.Rational(5, 9)*b - sp.Rational(1, 27)*a0
expect_zero('a5 preservation closed-form', sp.simplify(a5_sol - expected_a5_sol))
```

**Assessment:**
Matches the directive's "After" block exactly. The new assertion compares `a5_sol` (computed from `sp.solve(sp.Eq(coeff, 0), a5)`) against a closed-form hard-coded from the notes (`-5b/9 - a0/27`). This is non-tautological: a regression altering any of the three `coeff` coefficients (on `s`, `b`, `a0`, or `a5`) would propagate into `a5_sol` and break the anchor. The pre-existing tautological L55 check is retained as the directive instructed. The exec log line `a5 preservation closed-form = 0` is present at the expected position between `a5 preservation condition = -a0/27 - 5*b/9` and `preservation substitution = 0`.

### F2 — mathematica_transliteration

**Classification:** resolved

**What changed:**
`mathematica/moving_throat_pde_stage109_linearized_branch_selection_mathematica_audit.wl:40-63` replaces the previous series-of-the-ratio path with an independent numerator/denominator-separate path. Specifically:
- L42-47: builds `num`, `den`, `numLin`, `denLin`, `invDenLin`, then `chiSeries = numLin*invDenLin` (series-expanded).
- L55: `defectCoeff` is now extracted via `Coefficient[chiSeries, eps, 1]` rather than `(chiSeries - 1)/eps`.
- L60: `Solve[chiSeries - 1 == 0, a5, Reals]` operates directly on the series-minus-one instead of on the `coeff` intermediate.
- L63: substitution check now uses `(chiSeries - 1) /. a5 -> a5Pres` (no `coeff` reference).

**Assessment:**
Matches the directive's "After" block exactly. The `coeff = (chiSeries - 1)/eps` intermediate is fully gone — `grep` of the .wl shows only `defectCoeff` (a different quantity, derived via `Coefficient`). The algebraic route — separate-then-invert — is genuinely distinct from the SymPy `sp.series(chi, eps, 0, 2)` of the whole ratio. The four `PASS:` lines (`linearized chi law`, `overall scale cancels`, `a5 preservation condition + 5 b/9 + a0/27`, `preservation substitution`) are present in the exec log, exit code 0, and the printed `first-order defect coefficient = a0/3 + 9*a5 + 5*b` matches the SymPy side symbol-for-symbol. The `a5Pres = -1/27*a0 - (5*b)/9` value agrees with the closed-form. Scaffolding (banner, `$Assumptions`, ansatz L35-38, trailing Print/Exit) is intact.

### F3 — paper_misalignment / script_missing_paper_claim (Cluster A resolution)

**Classification:** resolved

**What changed:**
- `.py` L2-11: module docstring added stating that the card's secondary `Checks` items (Robin/mixed-pole limits, even-coefficient preservation under compensation) are exercised at stages 110, 111, and 112 respectively, with this stage establishing the linearized framework those downstream stages consume.
- `.wl` L28-30: mirror comment block in Mathematica with the same Cluster A carry-forward content (downstream stages 110/111/112).

**Assessment:**
The directive originally held F3 for user resolution ("Resolve before fix_loop") with three candidate directions; the user picked variant (c)-style routing — citing downstream stages where (b)/(c) are demonstrated. Both engines now carry mirror-image docstrings explaining the cross-stage referencing. This does not change any derived result; it documents the script-side scope so a future auditor can match `Checks` to verifying stages without alarm. No script logic added (consistent with not expanding this stage's `Inputs`).

### Cluster C — banner correction

**Classification:** resolved (collateral, expected)

**What changed:**
- `.py` L33: `STAGE 92` → `STAGE 109`
- `.wl` L26: `STAGE 092` → `STAGE 109`

**Assessment:**
Banner labels now correctly reflect the unit ID. Exec logs show the corrected banner. No mathematical content affected.

## Exec log assessment

**SymPy:** exit=0. Notable lines:
- `chi_Q series = a0*eps/3 + 9*a5*eps + 5*b*eps + 1`
- `linearized chi law = 0`
- `overall scale cancels = 0`
- `a5 preservation condition = -a0/27 - 5*b/9`
- `a5 preservation closed-form = 0`  (new F1 anchor)
- `preservation substitution = 0`

**Mathematica:** exit=0. Notable lines:
- `chi_Q series = 1 + (a0*eps)/3 + 9*a5*eps + 5*b*eps`
- `PASS: linearized chi law`
- `PASS: overall scale cancels`
- `first-order defect coefficient = a0/3 + 9*a5 + 5*b`
- `a5 preservation condition = -1/27*a0 - (5*b)/9`
- `PASS: a5 preservation condition + 5 b/9 + a0/27`
- `PASS: preservation substitution`
- `Stage 109 Mathematica audit passed.`

Four `PASS:` lines on the Mathematica side as expected.

**Output freshness:** Confirmed. Script mtimes: `.py` 2026-05-27 15:10:13, `.wl` 2026-05-27 15:10:26. Output mtimes: `.txt` (sympy) 2026-05-27 15:18:05, `.txt` (mathematica) 2026-05-27 15:24:26. Both outputs are newer than their corresponding scripts; the saved `.txt` outputs were re-generated post-fix.

## Material-change assessment

`material_change`: false.

The only derived quantities asserted by this stage — `Delta_Q = 5b + a_0/3 + 9 a_5` and `a_5 = -5b/9 - a_0/27` — are unchanged. F1 adds a new non-tautological anchor for an already-correct value; F2 reaches the same closed-form via an independent algebraic route; F3 is documentation-only; Cluster C corrects a label. No downstream unit's inputs change.

## Side observations (non-blocking)

- The carry-forward docstring asserts that stages 110/111/112 perform the Robin, mixed-pole, and even-coefficient checks respectively. This is a claim about downstream content that the verifier did not (and per scope cannot) cross-check against those stages' scripts here. If those stages do not in fact perform those checks, this docstring will become a future audit flag — but that is the auditor's job, not blocking F3 closure.
- The `.wl` retains the `expectZero` step for `(chiSeries - 1) /. a5 -> a5Pres` (still tautological by construction, since `a5Pres` solves that equation). The directive explicitly kept this — paired with the new anchored closed-form check, this is acceptable.

## Verdict justification

All three findings have correct, in-place edits matching the directive's "After" blocks exactly. F1's new SymPy anchor is non-tautological (the right-hand side is hard-coded from notes, not from `coeff`). F2's Mathematica path no longer shares the `coeff = (chiSeries - 1)/eps` intermediate and uses a genuinely distinct algebraic route (separate-then-invert vs. series-of-the-ratio). F3's Cluster A carry-forward documents the downstream verification routing without altering script logic. Cluster C banners are corrected. Both engines pass with exit 0, four `PASS:` lines on the Mathematica side, and outputs are fresh. No regressions in the diff. Verdict: `verified`.
