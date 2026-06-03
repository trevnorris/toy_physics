---
unit_id: 230
batch: VII.1
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-06-02T00:00:00-06:00
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 2
findings_total: 2
material_change: false
---

# Verification — unit 230

## Per-finding outcomes

### F1 — missing_verification_script (subtype: missing_mathematica)

**Classification:** resolved

**What changed:**
New independent Mathematica audit at
`mathematica/moving_throat_pde_stage230_selected_branch_classifier_to_dynamic_window_compiler_and_static_first_theorem_mathematica_audit.wl`
(9203 bytes, mtime 17:35:04). Covers the full M1–M6 claim manifest. The committed
output `mathematica/output/...mathematica_audit.txt` (17:40:49) and the exec log both
show exit 0 with 44 PASS tokens, 0 FAIL.

**Assessment:**
Genuinely independent — a different decomposition than the `.py`, not a transliteration:
- **M1** monotonicity proved by `Reduce[0<=xi<1 && delta>0 && D[rND,xi]<0,{xi,delta}]`
  and `Resolve[ForAll[...]]` that the Reduce region EQUALS the stable strip — a real
  sign proof, whereas the `.py` only matches `dR_dxi` against a hand-written derivative form.
- **M2** monotonicity by `Reduce[D[S±,R]<0 && R>=0]` + `Resolve[ForAll[R, reduce ⟺ R>=0]]`,
  not by reading off a negative constant's sign (the `.py` route).
- **M3** R_* obtained by independent `Solve[sMinus==0,R]`, nonnegative root selected, then
  asserted equal to the formula `sMinusDen/(-sMinusNum)`; plus `Reduce`-proved sign on both
  sides of R_*. Non-tautological.
- **M4** onset δ_* by independent `Solve[rNDAtOnset==rStar,delta]`, positive root selected,
  asserted equal to `8/(9 rStar)` — the Mathematica analogue of the F2 fix, NOT a back-substitution.
- **M5** ceilings via `Log`, `Limit[...,R->Infinity]`, and a direct `Which[]` piecewise
  construct (the directive explicitly forbade mirroring the `.py` `dynamic_ceilings` helper;
  the `.wl` does not).
- **M6** the "everywhere weaker" static-first claim reduced to the infimum: proves the
  CEILINGS `B±` are monotone-decreasing via `Reduce[D[B±,R]<0]` (even stronger than the
  `.py`, which only proves `S±` monotone), then `N`-compares the `R->Infinity` infima against
  the static budgets. A genuine monotonicity-to-infimum reduction, not a sampled claim.

Carry-forward literals (4 per-unit slopes, 3 R_Q figures, 2 static budgets) match the `.py`
exactly and are frozen Stage-228 inputs per the auditor's cross-check; the `.wl` re-derives
only the downstream algebra. No fabricated literal.

### F2 — tautological_check

**Classification:** resolved

**What changed:**
`scripts/...stage230..._sympy_audit.py:119-121`. The diff removed the round-trip
`assert sp.simplify(onset.subs(delta, delta_dyn_star) - R_star) == 0` and replaced it with:
```python
delta_solutions = sp.solve(sp.Eq(onset, R_star), delta)
assert len(delta_solutions) == 1
assert sp.simplify(delta_solutions[0] - delta_dyn_star) == 0
```
Matches the directive's prescribed block verbatim. Line 118 (`S_-(R_*)=0`) untouched; the
diff shows ONLY this 3-line replacement — no collateral edits.

**Assessment:**
De-tautologized and now genuine. The onset threshold is solved INDEPENDENTLY from
`8/(9*delta) == R_star` (the solver inverts the relation) and the recovered root is compared
to `delta_dyn_star`. If `delta_dyn_star` were defined incorrectly the assertion would fail,
because the solve returns the true root `8/(9*R_star)` regardless of `delta_dyn_star`'s
definition. The `len(...) == 1` guard also pins uniqueness. Non-tautological. The exec log
(line 23) reports δ_* = 0.723111617875019116..., a nonzero value, confirming the fixed check
is non-trivially satisfied.

## Exec log assessment

**SymPy:** exit=0. Notable lines:
- `R_*          = 1.22925543846333598788938878199`
- `delta_dyn_*  = 0.723111617875019116300604645303`
- `inf B_dyn^(both)     > B_stat^(both)` / `inf B_dyn^(nonempty) > B_stat^(nonempty)`
- `All Stage 230 symbolic and numerical audits passed.`

**Mathematica:** exit=0, 44 PASS / 0 FAIL. Notable lines:
- `M3 N[R_star,30] ... diff = 1.21...*^-17` PASS (independent Solve root = formula)
- `M4 solved onset root equals delta_dyn_star = 0` PASS (independent Solve)
- `M6 B_plus decreases on R >= 0 by Reduce = True` / `M6 B_minus decreases on finite nonempty region by Reduce = True`
- `R_star = 1.2292554384633359878893887819926314141...` / `delta_dyn_star = 0.723111617875019116...`
Cross-engine: R_* and δ_* agree to all printed digits across SymPy and Mathematica; the
infimum ceilings (0.96728238936382175..., 0.99058181070523368...) also match the SymPy values.

Token non-vacuity: confirmed. `fail[]` (`.wl` lines 8–12) calls `Exit[1]`; every check routes
through `expectZero`/`expectTrue`/`expectClose`/`expectInfinity`, each of which calls `fail`
on a false/`TrueQ`-unresolved result. `Exit[0]` is reached only after all checks. The SymPy
side uses bare `assert`/`raise AssertionError` (non-zero exit on failure). A bad check forces
a non-zero exit in both engines.

**Output freshness:** confirmed. SymPy `.txt` (17:40:49) and Mathematica `.txt` (17:40:49)
are both newer than their scripts (`.py` 17:33:40, `.wl` 17:35:04). Committed `.txt` content
matches the exec logs.

## Material-change assessment

`material_change`: false.

No derived result moved. F2 swapped a definitional round-trip for an equivalent independent
solve that recovers the SAME δ_*; F1 added a NEW second engine that confirms all existing
values (R_*, δ_*, endpoint ceilings, infima, static-first inequalities) rather than changing
any; the authorized notes renumber is prose-only (stale Stage 245/246/247 → canonical
228/229/230, values unchanged) and lies outside the scripts-only verifier scope. No
downstream-consumable quantity changed, so no unit > 230 is substantively affected by this fix.

## Side observations (non-blocking)

- The authorized notes renumber is recorded in the directive's `## Applied: notes-renumber`
  block; per scripts-only scope I did not read the notes file, but the captured diff.patch
  shows ONLY the `.py` changed — confirming the script and paper card were untouched by the
  renumber. Correctness of the prose renumber itself is a notes-review item, not a script
  verification item.
- The empty Codex log was flagged as a logging anomaly; all artifacts (the new `.wl`, the
  fixed `.py`, both fresh `.txt` outputs, both exit-0 exec logs) exist and are internally
  consistent, so the anomaly did not affect application. The directive's `## Applied: F1/F2/
  notes-renumber` blocks are all present with `applied:true / 2 applied`.

## Verdict justification

Both findings are resolved. F2 is genuinely de-tautologized: py:119 now solves the onset
relation independently and would fail on a wrong `delta_dyn_star`. F1's new `.wl` is a genuine
second engine using native primitives (`D`, `Reduce`, `Resolve[ForAll]`, `Solve`, `Limit`,
`Log`, `N`) through a different decomposition than the `.py`; M1–M6 are non-tautological, the
M6 static-first claim is a real monotonicity-to-infimum reduction, and the success helper
genuinely `Exit[1]`s on failure. R_* and δ_* match across engines, carried literals trace to
Stage 228 with no fabrication, both engines exit 0 (44 Mathematica PASS, 0 FAIL), and outputs
are fresh. No regressions in the diff. Verdict: verified; material_change false.
