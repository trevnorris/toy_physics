---
unit_id: 173
batch: V.1
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-05-28T16:20:00-06:00
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 2
findings_total: 2
material_change: false
---

# Verification — unit 173

## Per-finding outcomes

### F1 — mathematica_transliteration

**Classification:** resolved

**What changed:**
`mathematica/moving_throat_pde_stage173_axisymmetric_loading_mismatch_mathematica_audit.wl`:
- Lines 37-39: the differentiate-then-set-eps-zero extraction was replaced with direct
  series-coefficient extraction. Old `.wl` built `u2A,u4A,p0A` via
  `Expand[Normal[Series[...,{eps,0,1}]]]` then took `(D[u2A,eps]/.eps->0)/lam`; new `.wl`
  inlines the expressions into `Coefficient[Series[-d2A/d0A,{eps,0,1}]//Normal, eps, 1]/lam`
  (and likewise for `u41`, `p1`). The three `u2A/u4A/p0A` intermediates were removed.
- Lines 67-69: the intermediate general-solve `d41Hidden = ... First[Solve[u41Can==8*u21Can/9, d41]]`
  (and its print) was removed; `d41Even` is now obtained by solving
  `(u41Can == 8 u21Can/9) /. d21 -> u21ZeroD21` for `d41` directly — i.e. `d21` is fixed to its
  even-preserving value *before* the solve, rather than substituting into a previously stored
  general solution.
- `u21ZeroD21` (line 64) was kept, as the directive explicitly permitted.

**Assessment:**
The change matches the directive's "required change" precisely. The git diff (8 insertions, 13
deletions) touches only the F1-named line ranges plus the F2 banner; there is no collateral edit.

The restructure achieves genuine route-independence:
- `.py` route: `series(...).removeO()` into named intermediates, then
  `diff(expr,eps).subs(eps,0)/lam`. `.wl` route: `Coefficient[Series[...],eps,1]/lam` — a
  distinct coefficient-extraction path, not a differentiate-then-evaluate path.
- `.py` even-preserving: solve hidden-even generally for `D41` (`D41_hidden`), solve
  `u21_can=0` for `D21`, then `D41_hidden.subs(D21, ...)`. `.wl` even-preserving: solve
  `u21Can==0` for `d21` first, pre-substitute, then a single solve for `d41`. Distinct order
  and mechanics; the mirrored sequential `Solve`-then-substitute choreography is gone.

The directive's three verification conditions all hold:
(a) `grep` confirms the `.wl` no longer contains `D[u2A`/`D[u4A`/`D[p0A` nor `d41Hidden`.
(b) All six `expectZero` lines (.wl:45,56,57,61,72,73) are byte-for-byte the paper closed-form
    targets they were before the edit — verified against the original report's cross-check table
    and the directive's "must stay byte-for-byte" requirement. None became tautological: each LHS
    (`u21`, `u41Can`, `p1Ratio`, `hiddenEvenResidual`, `u21ZeroD21`, `d41Even`) is independently
    derived via Series/Coefficient or Solve, then compared to the paper RHS.
(c) The captured transcript ends "Stage 173 Mathematica audit passed." with 6 PASS lines and no
    FAIL/`$Failed`/residual markers; `Exit[0]`.

Non-rubber-stamp check: the six engine final forms in the new `.wl` transcript still match the
`.py` transcript symbolically (u2^(1), u4^(1), P1/P0, Xi_load, lane forms, D21=-D01/9,
D41=-D01/27), confirming the new route reaches the same results — which is the point of an
independent second engine.

### F2 — stale_output (cosmetic banner mislabel)

**Classification:** resolved

**What changed:**
- `scripts/.../sympy_audit.py:30`: `banner("STAGE 156 — ...")` → `banner("STAGE 173 — ...")`.
- `mathematica/.../mathematica_audit.wl:26`: `banner["STAGE 156 — ..."]` → `banner["STAGE 173 — ..."]`.

**Assessment:**
Exactly the directive's required string change; only the stage number changed, the rest of the
banner text and the surrounding line are intact. Both regenerated transcripts now read
"STAGE 173 — WEAK-AXISYMMETRIC PHYSICAL-SLOPE TRANSPORT" on line 3 (the banner emits a blank
line + rule before the title, so the title lands at transcript line 3, not line 11 as the report
estimated — the content is correct regardless). Pure string edit; no math touched.

## Exec log assessment

**SymPy:** exit=0. Notable lines:
- `u2 slope identity = 0` (line 8)
- `u4 canonical formula = 0`, `P1/P0 formula = 0` (lines 16-17)
- `hidden-even residual = 0` (line 22)
- `D21 + D01/9 = 0`, `D41 + D01/27 = 0` (lines 30-31)
All six `expect_zero` residuals print `= 0`; the run completes through the carry-forward block
with no `AssertionError`/traceback, so the script exited 0.

**Mathematica:** exit=0. Notable lines:
- `PASS: u2 slope identity` (line 9)
- `PASS: u4 canonical formula`, `PASS: P1/P0 formula` (lines 18, 20)
- `PASS: hidden-even residual` (line 26)
- `PASS: D21 + D01/9`, `PASS: D41 + D01/27` (lines 34, 36)
- `Stage 173 Mathematica audit passed.` (line 54)
6 PASS lines, no FAIL/`$Failed`/`residual ->` markers; `Exit[0]`.

**Output freshness:** confirmed fresh post-fix. `.wl` script mtime 15:55:37, its transcript
16:13:43 (newer). `.py` script mtime 15:55:39, its transcript 16:10:20 (newer). Both outputs were
regenerated after the edits.

## Material-change assessment

`material_change`: false.

F1 changed only the Mathematica derivation *route*; the final `expectZero` targets (the paper
closed forms) are byte-for-byte identical, and the new transcript's symbolic results match the
old ones and the SymPy results. No derived quantity changed value. F2 is a pure transcript-label
edit. No downstream unit can observe any change in result. Nothing to flag for re-audit.

## Side observations (non-blocking)

- The report (F2 "Verification") predicted the banner title would sit on transcript line 11; in
  the actual transcripts the title is line 3 (the `banner` helper prints a leading blank line plus
  a rule before the title). This is an estimate slip in the report text, not a defect in the fix —
  the title string itself is correct in both transcripts. Non-blocking.

## Verdict justification

Both findings are `resolved`. F1's restructure matches the directive exactly: the mirrored
`D[...]/.eps->0` first-order extraction and the `d41Hidden` general-solve are gone, replaced by a
genuinely distinct `Coefficient[Series[...]]` route and a pre-fixed-`d21` single solve, while the
six `expectZero` paper-form targets remain byte-for-byte and non-tautological. F2 is the exact
banner string change in both files. Both engines exit 0 with all six checks passing, the
transcripts are fresh, the git diff shows no collateral edits, and no derived result changed
(`material_change: false`). Verdict: verified.
