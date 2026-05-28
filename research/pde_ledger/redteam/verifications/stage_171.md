---
unit_id: 171
batch: V.1
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-05-28T16:40:00-06:00
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 2
findings_total: 2
material_change: false
---

# Verification — unit 171

Re-verification after the orchestrator reworked F1 (the first verification
returned `needs_rework` because the directive's prescribed F1 fix had rebuilt
`zCombFormula`/`nCombFormula` from the same `D[z2,...]`/`D[z0,...]` calls that
build `zCombExact`/`nCombExact`, making the bundle checks `X − X`). This pass
confirms the rework restores a load-bearing, non-tautological check and that F2
is still resolved.

## Per-finding outcomes

### F1 — mathematica_transliteration

**Classification:** resolved

**What changed:**
The two bundle-target literals in the `.wl` are no longer `D[...]`-based
reconstructions. The tautological iteration-1 helpers (`dZ2Indep`, `dZ0Indep`,
`dN0Indep`, `dZ0Indep2`) are deleted — none appear anywhere in the file. The
current state is:

- `mathematica/...audit.wl:117-122` — `zCombFormula` is a hand-collected
  closed-form literal:
  `(s/delta^2 + 1/(9*delta))*dQ + (q/delta^2)*dS - dG/delta + (gSym/delta^2 - q/(9*delta^2) - 2*q*s/delta^3)*dDelta`.
  It is built with NO `D[...]` call.
- `mathematica/...audit.wl:123` — `expectZero["Z obstruction bundle", zCombExact - zCombFormula]`,
  where `zCombExact = Expand[dZ2Exact + dZ0Exact/9]` (wl:113) and `dZ2Exact`
  (wl:100) / `dZ0Exact` (wl:99) are the engine's own `D[]`-derived total
  differentials of `z2`/`z0`.
- `mathematica/...audit.wl:127-130` — NEW second route: `zCombSeries` extracts
  the first-order `eps2` coefficient of `(z2 + z0/9)` after substituting
  `q→q+eps2*dQ`, `s→s+eps2*dS`, `gSym→gSym+eps2*dG`, `delta→delta+eps2*dDelta`
  via `Series`, then `expectZero["Z obstruction bundle (series route)", zCombSeries - zCombFormula]`.
- `mathematica/...audit.wl:134-137` — `nCombFormula` is likewise a hand-collected
  literal `(p0/delta)*dQ + (2*p/delta^2)*dP - (p0*q/delta^2 + 2*p^2/delta^3)*dDelta`,
  no `D[...]`.
- `mathematica/...audit.wl:138` — `expectZero["N obstruction bundle", nCombExact - nCombFormula]`,
  `nCombExact = Expand[dN0Exact + p0*dZ0Exact]` (wl:132) is engine-`D[]`-derived.
- `mathematica/...audit.wl:140-143` — NEW second route `nCombSeries` via `Series`
  on `(n0Expr + p0*z0)` with `q,p,delta` perturbed; compared to `nCombFormula`.

The SymPy file is untouched on F1 (correct — the asymmetry is supposed to live
in the `.wl`). Diff vs HEAD is exactly these additions plus the F2 banner; no
collateral edits.

**Assessment:**
Non-tautological and load-bearing now. Three independent points:

1. `zCombFormula`/`nCombFormula` are collected closed-form literals, NOT
   reassemblies of the `D[z2,...]`/`D[z0,...]`/`D[n0Expr,...]` outputs that feed
   `zCombExact`/`nCombExact`. The coefficients are pre-combined by hand — e.g.
   the `dQ` coefficient `s/delta^2 + 1/(9*delta)` merges `∂z2/∂q = s/delta^2`
   with `(1/9)∂z0/∂q = 1/(9 delta)`; the `dDelta` coefficient
   `gSym/delta^2 - q/(9*delta^2) - 2*q*s/delta^3` merges the `z2` and weighted
   `z0` delta-slopes. These literals are not produced by any differentiation
   call in the script, so `zCombExact - zCombFormula` is genuinely
   `(engine-derived) − (independent target)`, not `X − X`.

2. The check FAILS if a collected coefficient is wrong. Concrete probe: if the
   `dDelta` coefficient dropped the `- q/(9*delta^2)` piece (the `1/9`-weighted
   `z0` contribution), `zCombExact` (which includes it via `dZ0Exact/9`) would no
   longer match and `FullSimplify[Together[...]]` would leave a `q/(9 delta^2)`
   residual ⇒ `fail[...]` ⇒ `Exit[1]`. Same for the `dQ` `1/(9*delta)` term, the
   `dS` `q/delta^2`, the `dG` `-1/delta`, and every N-bundle coefficient. The
   check is sensitive to each collected coefficient.

3. The added `Series`-route check is a genuinely distinct mechanism: it
   Taylor-expands the substituted primitive expression in a perturbation
   parameter and reads off the linear coefficient, rather than summing
   per-variable partials. `zCombSeries`/`nCombSeries` are compared against the
   SAME collected target, so a wrong target coefficient also fails this second
   route — and it cross-checks the `D[]`-summed route against the Series route on
   the closed forms `z2`, `z0`, `n0Expr`. This is strictly stronger than the
   original auditor ask (which only wanted the second engine to re-derive the
   target). The series substitutions cover exactly the variables each expression
   depends on (`q,s,gSym,delta` for `z2+z0/9`; `q,p,delta` for `n0Expr+p0*z0`),
   with non-dependencies correctly omitted.

The exec log shows both `Z obstruction bundle`, `Z obstruction bundle (series
route)`, `N obstruction bundle`, and `N obstruction bundle (series route)`
printing `= 0` / `PASS` (mathematica `.txt:31-38`). The vacuous-tautology failure
mode flagged in the prior verification is gone. F1 resolved.

### F2 — insufficient_verification (mislabeled unit)

**Classification:** resolved

**What changed:**
- `scripts/...audit.py:27` — banner now `"STAGE 171 — MICROSCOPIC GROUPED OUTLET OBSTRUCTIONS"`.
- `scripts/...audit.py:5` — docstring now reads "obstructions for Stage 171."
  (the "from Stage 170" reference is gone).
- `mathematica/...audit.wl:26` — banner now `banner["STAGE 171 — MICROSCOPIC GROUPED OUTLET OBSTRUCTIONS"];`.

**Assessment:**
Correct and complete; unchanged from the prior pass where it was already
resolved. Both saved transcripts print `STAGE 171 — MICROSCOPIC GROUPED OUTLET
OBSTRUCTIONS` at `.txt:3` (sympy and mathematica). No assertion line was touched.
The Mathematica closing message at wl:172 already said "Stage 171" and is
consistent. Traceability defect fully closed.

## Exec log assessment

**SymPy:** exit=0. Notable lines:
- `.txt:3` — `STAGE 171 — MICROSCOPIC GROUPED OUTLET OBSTRUCTIONS` (F2 fix present)
- `.txt:18` — `Z obstruction bundle = 0`
- `.txt:19` — `N obstruction bundle = 0`
- All 21 checks print `= 0`; the carry-forward block (`.txt:27-30`) prints after
  the last check, so `expect_zero` never raised ⇒ exit 0. (The SymPy file was not
  edited for F1; its bundle checks remain genuine literal-vs-`sp.diff`
  comparisons, py:134-148.)

**Mathematica:** exit=0. Notable lines:
- `.txt:3` — `STAGE 171 — MICROSCOPIC GROUPED OUTLET OBSTRUCTIONS` (F2 fix present)
- `.txt:31-32` — `Z obstruction bundle = 0` / `PASS: Z obstruction bundle`
- `.txt:33-34` — `Z obstruction bundle (series route) = 0` / `PASS: ...` (new)
- `.txt:35-36` — `N obstruction bundle = 0` / `PASS: N obstruction bundle`
- `.txt:37-38` — `N obstruction bundle (series route) = 0` / `PASS: ...` (new)
- `.txt:57` — `Stage 171 Mathematica audit passed.`
- 23 PASS lines total (was 21 + the 2 new series routes), 0 FAIL. `expectZero`
  calls `fail[...]`→`Exit[1]` on nonzero; no FAIL appeared and the script reaches
  the final `Exit[0]` (wl:174) ⇒ exit 0.

**Output freshness:** confirmed. The `.wl` mtime is 2026-05-28 16:26:29; its
`.txt` mtime is 16:26:44 (regenerated ~15s after the edit). The SymPy `.py` mtime
is 15:55:18 with `.txt` at 16:10:18. Both outputs postdate their scripts, so the
logs reflect the post-rework state. MANIFEST records no per-run exit field for
this unit, but the logs are conclusive (no FAIL, clean completion in both).

## Material-change assessment

`material_change`: false.

No derived numeric or symbolic result changed. F2 is a banner/docstring label
fix (no assertion touched). The F1 rework changed only how the comparison
TARGET is expressed and added a second cross-check route inside the `.wl`; the
defining expressions (`z0`, `z2`, `n0Expr`), the engine-derived `zCombExact`/
`nCombExact`, the collected closed forms, and all carry-forward formulas are
unchanged, and the SymPy file is unchanged. The collected literal is identical
to the (unchanged) SymPy `Zcomb_formula`/`Ncomb_formula`, confirming the target
is the same paper bundle, not a new value. Downstream units consuming the Stage
171 carry-forward (`K_A`/`G_A` regrouping, BdG/Z/N bundles, weak-axisymmetric
collapse) see no changed value. No `upstream_stale` propagation warranted on
physics grounds.

## Side observations (non-blocking)

- The `.wl` `zCombFormula`/`nCombFormula` collected literals are still textually
  the same as the SymPy `Zcomb_formula`/`Ncomb_formula`. This is fine and
  intended: the original F1 (severity low) accepted that the genuine
  cross-engine signal lives in each engine running its own differentiation, and
  the rework's added `Series` route gives the `.wl` a second, mechanism-distinct
  derivation of the bundle. A copied-but-wrong target would now be caught by BOTH
  the `D[]`-derived `zCombExact` side and the `Series` side, so the original
  "copied literal" weakness is materially mitigated within the Mathematica engine.
- The regrouping checks (K_A split, G_A split, weak-axisymmetric lanes) remain
  the same linear-combination-vs-reparenthesization style the original report
  rated "weak but not tautological." Unchanged and out of scope; not a basis for
  rework.

## Verdict justification

Both findings are resolved. F2 (banner/docstring) was already correct and the
fresh transcripts both read Stage 171. The F1 rework is correct and
non-tautological: `zCombFormula`/`nCombFormula` are hand-collected closed-form
literals (no `D[...]`), so the `... bundle` checks compare the engine's own
`D[]`-derived total differential against an independent target and would fail on
any wrong collected coefficient; and the added `... bundle (series route)`
checks supply a genuinely distinct second derivation (Taylor `Series` in a
perturbation parameter) compared against the same target. The tautological
helpers from iteration 1 are deleted. The exec log shows all bundle PASS lines
including the two new series routes (23 PASS, 0 FAIL) and both engines exit 0,
with outputs freshly regenerated after the edits. No material change to any
derived result. Verdict: verified.
