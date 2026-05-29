---
unit_id: 109
batch: IV.2
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-05-29T00:00:00Z
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 2
findings_total: 3
material_change: false
---

# Verification — unit 109

## Per-finding outcomes

### F1 — tautological_check (SymPy)

**Classification:** resolved

**What changed:**
`scripts/moving_throat_pde_stage109_linearized_branch_selection_sympy_audit.py:51-56`. The
self-solved root `a5_sol = sp.solve(sp.Eq(coeff, 0), a5)[0]` remains (L51) but the two follow-on
assertions are now de-tautologized:
- L53-54: `expected_a5_sol = -sp.Rational(5,9)*b - sp.Rational(1,27)*a0` and
  `expect_zero('a5 preservation closed-form', sp.simplify(a5_sol - expected_a5_sol))` — compares the
  solved root to the INDEPENDENT closed form from the notes.
- L56: `expect_zero('preservation substitution', sp.simplify(coeff.subs(a5, expected_a5_sol)))` — the
  substitution argument is now `expected_a5_sol`, NOT the self-solved `a5_sol`.

**Assessment:**
Correct; matches the directive's "After" block exactly. Both required checks present. The
"preservation substitution" is genuinely non-tautological now: it substitutes the hardcoded closed
form into `coeff` (derived from the ansatz at L46 = `(chi_series-1)/eps`). A structural error in
`coeff` would make this nonzero rather than silently passing. The closed-form check would fail if
the solved root disagreed with `a_5 = -5b/9 - a_0/27`. Diff confirms only L53-onward changed; no
collateral edits. Exec log: `a5 preservation closed-form = 0`, `preservation substitution = 0`, both
pass.

### F2 — mathematica_transliteration (.wl)

**Classification:** resolved

**What changed:**
`mathematica/moving_throat_pde_stage109_linearized_branch_selection_mathematica_audit.wl:40-64`. The
derivation path is algebraically distinct from the `.py`: `num`/`den` expanded separately (L42-43),
linearized separately (L44-45), denominator inverted via its own series `invDenLin =
Normal[Series[1/denLin, {eps,0,1}]]` (L46), then multiplied and re-seriesed (L47) — versus the
`.py`'s single direct ratio `sp.series` of the whole quotient. `defectCoeff` is read off via
`Coefficient[chiSeries, eps, 1]` (L55). The preservation solve is direct on `Solve[chiSeries - 1 ==
0, a5]` (L60). The "preservation substitution" (L64) now substitutes the literal closed form
`(chiSeries - 1) /. a5 -> (-5*b/9 - a0/27)`, NOT `a5Pres`.

**Assessment:**
Correct; matches the directive's "After" block. Case-sensitive grep confirms the `(chiSeries -
1)/eps` intermediate is gone (0 matches); the only `coeff`-bearing token is the distinct
`defectCoeff` plus prose comments. The separate-then-invert path is genuinely distinct from the
`.py` single-ratio series, satisfying the second-engine policy. The de-tautologized substitution
depends on `chiSeries` being correct — a structural error breaks it. Diff confirms only the L63→L64
substitution line changed (plus its comment); all surrounding scaffolding intact. Exec log: 4 PASS /
0 FAIL.

### F3 — paper_misalignment (resolved as cross-reference, NO script change)

**Classification:** resolved (correctly untouched)

**What changed:** Nothing in any script/paper/notes for F3. The directive has NO `## Applied: F3`
block (confirmed by grep). The diff touches only the F1 `.py` and F2 `.wl` files — no F3 edits.
Direction (c) routes the three secondary card Checks to sibling owners (scale/argument→108,
Robin→110, mixed-pole→111, compensated even+odd→112) with a paper-card cross-reference in the manual
paper pass (PAPER_CLEANUP_TRACKER).

**Assessment:** Correct. The script header docstrings already mention the downstream
cross-reference (110/111/112), but these predate this diff and document the agreed direction (c);
they are not part of the stage-109 diff under review and are not an F3 script "fix." Verified
untouched per directive.

## Load-bearing non-tautology check

Confirmed: the real anchor is `expectZero["linearized chi law", chiSeries - expected]` (.wl L51) and
`expect_zero('linearized chi law', chi_series - expected)` (.py L44), with `expected = 1 + eps*(5*b +
a0/3 + 9*a5)` hardcoded — the paper's first-order law — in both engines (.wl L48, .py L43).
`chiSeries`/`chi_series` is derived from the (S, beta, Sigma0, Sigma5) ansatz independently of
`expected`, so this check would FAIL if the linearization were wrong. The de-tautologized
preservation-substitution checks (.py L56, .wl L64) substitute the independent closed form into the
ansatz-derived `coeff`/`chiSeries`; a structural error breaks them. They are no longer X−X.

## Exec log assessment

**SymPy:** exit=0. Notable lines:
- `linearized chi law = 0` (load-bearing anchor passes)
- `a5 preservation condition = -a0/27 - 5*b/9` (matches notes' closed form)
- `a5 preservation closed-form = 0` and `preservation substitution = 0` (both de-tautologized checks pass)

**Mathematica:** exit=0, 4 PASS / 0 FAIL. Notable lines:
- `PASS: linearized chi law`
- `PASS: overall scale cancels`
- `PASS: a5 preservation condition + 5 b/9 + a0/27` (with `a5 preservation condition = -1/27*a0 - (5*b)/9`)
- `PASS: preservation substitution`
All four PASS lines are content-bearing (each named, preceded by its printed residual `= 0`).

**Output freshness:** confirmed. `scripts/output/...sympy_audit.txt` (mtime 2026-05-29 11:23:13) and
`mathematica/output/...mathematica_audit.txt` (mtime 2026-05-29 11:23:32) are both newer than their
source scripts (both 2026-05-29 10:57:28). Outputs regenerated post-fix.

## Material-change assessment

`material_change`: false. The edits strengthen assertions and reorganize the `.wl` derivation, but
the derived results are unchanged: the `chi_Q` first-order law (`5b + a0/3 + 9a5`) and the
preservation closed form (`a_5 = -5b/9 - a_0/27`) are identical pre- and post-fix. No downstream
unit consumes a changed value. (Orchestrator will still mark units > 109 upstream_stale per policy;
no specific concern here.)

## Side observations (non-blocking)

- `redteam/codex_reviews/stage_109.md` reflects an EARLIER iteration: its verdict table flags
  `.py:55`/`.wl:63` as still tautological (substituting `a5_sol`/`a5Pres`). The CURRENT source and
  diff show both were subsequently de-tautologized to the independent closed form. The current state
  supersedes that review (consistent with an iterate-to-clean loop).
- The prior `verifications/stage_109.md` (2026-05-27) likewise described the substitutions as
  retained-tautological; this verdict reflects the newer directive/source state where both are fixed.
- Both scripts carry a header docstring documenting the F3 downstream cross-reference (110/111/112).
  This documents direction (c), not a behavior change, and is not part of the stage-109 diff.

## Verdict justification

Both findings are genuinely resolved. F1 and F2 each replace the self-solved-root substitution with
the INDEPENDENT closed form `a_5 = -5b/9 - a_0/27`, and F1 adds a direct closed-form comparison of
the solved root; both now depend on the ansatz-derived `coeff`/`chiSeries` and would fail under a
structural error. The `.wl` path is algebraically distinct (separate num/den linearization plus
denominator-inverse series; the `(chiSeries-1)/eps` intermediate is gone — grep confirms 0). F3 was
correctly left untouched (no `## Applied: F3`, no F3 edits in the diff). Both exec logs exit 0 with
all content-bearing checks passing, and saved outputs are fresher than the sources. Verdict:
verified.
