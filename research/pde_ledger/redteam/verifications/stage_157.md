---
unit_id: 157
batch: IV.6
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-05-29T23:10:00-06:00
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 3
findings_total: 3
material_change: false
---

# Verification — unit 157

## Per-finding outcomes

### F1 — insufficient_verification (SymPy: canonical-even check re-solves the same pair)

**Classification:** resolved

**What changed:**
- `scripts/...stage157...sympy_audit.py:113-130` — the old duplicate re-solve
  (`sol_deltaC = sp.solve([...dE2,dE4...]); expect_zero("delta C from canonical-even Solve", ...)`)
  is REMOVED and replaced by `even_det = sp.Matrix([[1, -9*sigma_star], [5, -72*sigma_star]]).det()`
  with `expect_zero("canonical-even non-degeneracy: ... (det = -27*sigma_star)", even_det + 27*sigma_star)`.
- py:12-14 — docstring item 6 corrected from "Tangent motion kills delta C and forces delta kappa_W = 0"
  to the deferral wording ("Canonical-even preservation pins delta C = delta kappa_W = 0 (non-degenerate
  kernel); ... deferred to Stage 158").
- py:64 — banner corrected from "Stages 138-139" to "Stages 155-156" (housekeeping item).

**Assessment:**
Correct and matches the directive's consult-Q3 option (ii). The new assertion is GENUINELY
non-tautological: it operates on the coefficient-matrix determinant (an independent quantity),
NOT a re-solve of the homogeneous pair already asserted at py:108-111, and does NOT re-`expect_zero`
the already-asserted-zero `deltaC`. A mistyped coefficient (e.g. -72 → -71) gives det = -36*sigma_star,
so `even_det + 27*sigma_star = -9*sigma_star ≠ 0` and the check FAILS — a real fail mode. The
constraint-imposition solve+assert at py:108-111 is still present and passing. No collateral edits
beyond the three named changes.

### F2 — transliteration (Mathematica: mirror of the same literal 2×2 system)

**Classification:** resolved

**What changed:**
`mathematica/...stage157...audit.wl:102-114` — the mirrored
`Clear[dCsym,dKsym]; solDeltaC = Solve[{dCsym - 9 sigmaStar dKsym == 0, 5 dCsym - 72 sigmaStar dKsym == 0}, ...];
deltaCIndep = ...; expectZero["delta C from canonical-even Solve", deltaCIndep]` block is GONE,
replaced by the parallel `evenDet = Det[{{1, -9 sigmaStar}, {5, -72 sigmaStar}}]` with
`expectZero["canonical-even non-degeneracy: ... (det = -27 sigmaStar)", evenDet + 27 sigmaStar]`.

**Assessment:**
Correct and parallels F1. The `.wl` no longer presents the literal 9/72/5 numerator system as
an "independent" `expectZero`; it asserts the determinant with the same genuine loss-of-full-rank
fail mode. The even-preservation constraint solve+assert at wl:96-100 is still present and passing.

### F3 — symbol_assumption_error (Mathematica: missing physical branch domain)

**Classification:** resolved

**What changed:**
`wl:93` — `$Assumptions = Element[{sigmaStar, deltaC, dKappa}, Reals]` → `... && 0 < sigmaStar < 1`.
Scoped to the Section-3 even-preservation block only (Sections 1/2/4 `$Assumptions` untouched).

**Assessment:**
Correct (consult Q4). No bare `ConditionalExpression` appears in the printed residuals (the
Section-3 `expectZero` lines print "= 0" / "PASS" cleanly), and the script still exits 0.

## Exec log assessment

**SymPy:** exit=0. Notable lines:
`canonical-even preservation solutions = [{deltaC: 0, delta_kappa: 0}]`
`canonical-even non-degeneracy: trivial kernel forces delta C = 0 (det = -27*sigma_star) = 0`
`2. Carry-forward numerical basepoint from Stages 155-156`

**Mathematica:** exit=0. Notable lines:
`canonical-even non-degeneracy: trivial kernel forces delta C = 0 (det = -27 sigmaStar) = 0`
`PASS: canonical-even non-degeneracy: ...`  /  `Stage 157 Mathematica audit passed.`
No `ConditionalExpression` in any residual.

**Output freshness:** confirmed. Both saved `.txt` outputs (mtime 1780117330) are newer than the
edited scripts (mtime 1780116781). The committed outputs contain the new non-degeneracy line and
the "Stages 155-156" banner; the old "delta C from canonical-even Solve" string is absent from both.

## Material-change assessment

`material_change`: false. All three edits are how-it's-checked / labeling fixes: a tautological
re-solve and its Mathematica mirror were collapsed to an independent determinant non-degeneracy
assertion, a docstring/banner were aligned to the already-deferring published card, and a physical
branch domain was scoped onto an assumption. No derived numeric or symbolic result that downstream
units consume was changed. Per the directive's `## RESOLVED (consult batch 7)` block, the published
card (Open/Numerical, deviation-to-normalization map already deferred to Stage 158) is unchanged and
already correct — escalation was resolved against; this is not a conceptual change.

## Side observations (non-blocking)

- The directive notes `notes/stages/review/stage_157_review.md` is orphaned (contains a Stage 038
  review body). That is an orchestrator/notes repair, out of scripts-only scope — not blocking.

## Verdict justification

All three findings are resolved with edits matching the directive exactly and no collateral changes
(the diff shows only the four named edits). The replacement assertions are genuinely non-tautological
in both engines — they operate on the carried coefficient-matrix determinant with a real fail mode on
coefficient mistype/rank-loss, not a re-solve of the previously-asserted homogeneous kernel. The
docstring item 6, the SymPy banner, the Mathematica mirror removal, and the physical-branch
assumption are all confirmed; both engines exit 0 and the refreshed committed outputs reflect the
post-fix state. Verdict: verified.
