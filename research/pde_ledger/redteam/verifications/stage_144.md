---
unit_id: 144
batch: IV.5
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-05-29T22:00:00Z
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 1
findings_total: 1
material_change: false
---

# Verification — unit 144

## Per-finding outcomes

### R1 — mathematica_transliteration (load-bearing transcendental root)

**Classification:** resolved

**What changed:**
Two edits to `mathematica/moving_throat_pde_stage144_unique_regular_canonical_branch_mathematica_audit.wl`, matching the directive verbatim (diff `redteam/exec_logs/stage_144_diff.patch`, 19 insertions / 2 deletions, only this `.wl` touched):

- **R1a (`.wl:52-54`):** The old transliterated `piStar = piM /. FindRoot[gPi == gMinus, {piM, 1.5}, ...]` (point seed, same rational equation as SymPy) is REPLACED by a cleared-denominator bracketed solve:
  ```
  gThresholdResidual[p_] := 2*p*(2*p*Exp[p] + Pi) - gMinus*(4*p^2 + Pi^2)*(Exp[p] - 1);
  piStar = piM /. FindRoot[gThresholdResidual[piM] == 0, {piM, 1.4, 1.6}, WorkingPrecision -> 80, ...];
  ```
  SIGN CHECK CONFIRMED: residual keeps `(Exp[p] - 1)` (the positive factor), NOT `(1 - Exp[p])` — the §131-F3 sign-error trap is avoided.
- **R1b (`.wl:79-84`):** The old single check against the truncated literal `1.50882951349316` is REPLACED by two asserts: (i) a residual-near-zero witness `piStarResidual = Chop[N[gThresholdResidual[piStar],30],10^-30]` with gate `=== 0 || Abs[...] < 10^-25`, and (ii) an anchor to the stage-131 OWNING value `piStarOwner = N[Rationalize[1.50882951349315558300555075595,0],30]` at tol `10^-12`.

**Assessment:**
Correct and non-tautological. Specifically:

- The OLD transliterated route is gone: `grep "FindRoot[gPi == gMinus"` returns no match; `grep "1.50882951349316"` returns no match anywhere in the `.wl` (old truncated anchor literal eliminated).
- Residual witness PASSES in the log (line 25: `PASS: Pi_* solves cleared-denominator residual (independent root)`). The witness gate can fail on a sign error or non-root convergence, so it is a genuine convergence/sign check, not a restatement. The residual is ~0 (Chop'd / `< 10^-25`), NOT ~6366 — the §131-F3 sign signature is absent.
- The owning-value anchor PASSES (log line 26: `PASS: Pi_* matches stage-131 owning value (independent route)`). The literal `1.50882951349315558300555075595` is exactly the stage-131 OWNING value sourced at `scripts/output/moving_throat_pde_stage131_parent_mouth_threshold_sympy_audit.txt:2` (`Pi_* = 1.50882951349315558300555075595`), verified directly. It is NOT 144's own SymPy nsolve output.
- Genuine independence (MIRROR POLICY): (a) different equation FORM — Mathematica solves the cleared-denominator polynomial-in-`(piM, Exp[piM])` residual; SymPy `nsolve`s the rational `gPi - gminus` with denominator left in; (b) different engine; (c) bracketed seed `{piM,1.4,1.6}` vs SymPy point seed `1.5`; (d) anchored to a DIFFERENT stage's (131) independent computation, not to 144's own solver. The numeric witness corroborates the independence: the Mathematica root lands on `...5830055507559...` (matching the stage-131 owner to >29 digits) while 144's own SymPy nsolve produces `...552747043511772` (committed sympy `.txt:17`), which diverges from the Mathematica/131 value at ~digit 16. A copied root would not show this — the Mathematica value tracks stage 131, NOT 144's SymPy, confirming the two routes are distinct and the anchor is non-self.
- The SymPy `.py` was correctly left UNTOUCHED (mtime 2026-05-27, last committed 3e2b5c0; `git diff --stat` shows only the `.wl` changed). Its branch/target asserts (`.py:53-67`) and PASS line (`.py:68`) still hold (sympy log line 23 / committed `.txt:18` PASS, exit 0). The SymPy `Pi_*` anchor remains the truncated notes display literal `1.50882951349316` at tol `10^-12`, which is correct for the SymPy side per the directive's Reconcile block (only the cross-engine Mathematica anchor is upgraded to the owning value).
- Collateral edits: none beyond the named R1 edit (diff confirms).

## Exec log assessment

**SymPy:** exit=0. Notable lines:
- `Pi_*          = 1.50882951349315552747043511772` (144's own nsolve — diverges from owner at ~digit 16, as expected; not used as the Mathematica anchor)
- `PASS: numerical-target assertions (upper/lower branches, Pi_*, Sigma0_*, That_*, Pi_match, That_match)`

**Mathematica:** exit=0. 8 PASS lines, two of which are the new Pi_* lines:
- `PASS: Pi_* solves cleared-denominator residual (independent root)` (residual ~0 witness — sign correct)
- `PASS: Pi_* matches stage-131 owning value (independent route)`
- plus `upper branch g_+^F1 > 1`, `lower branch bracket 2/pi < g_-^F1 < 1`, `That(Pi_*)`, `Sigma0(Pi_*)`, `Pi_match`, `That(Pi_match)`
- `Pi_* = 1.5088295134931555830055507559...` — matches stage-131 owner.

**Output freshness:** confirmed re-generated post-fix. `.wl` mtime 2026-05-29 21:30:01; committed mathematica `.txt` mtime 2026-05-29 21:34:09 (newer). SymPy `.py` mtime 2026-05-27 19:50 (untouched); committed sympy `.txt` mtime 2026-05-29 21:34:02 (refreshed by orchestrator re-run). Both exec logs dated 2026-05-29T21:34, exit_code 0, and the committed `.txt` PASS lines match the logs.

## Material-change assessment

`material_change`: false.

No derived numeric result changed. `Pi_*` is unchanged (the Mathematica root is the same `1.5088295134931555830055...`); only the DERIVATION ROUTE (cleared-denominator bracketed FindRoot vs transliterated point FindRoot) and the ANCHOR (stage-131 owning value vs 144's truncated notes literal) changed. No downstream unit can depend on a value that did not move.

## Side observations (non-blocking)

- The two new Pi_* asserts give the Mathematica engine 8 PASS lines vs the SymPy engine's single rolled-up `PASS: numerical-target assertions` line. Cosmetic granularity asymmetry only; the underlying claims are the same set. Not a finding.
- `piMatch` is still found by `FindRoot[gPi == Pi/4, {piM, 1.9}]` (a point-seed transliteration of the SymPy `Pi_match` route). This was explicitly OUT of scope for R1 (the directive scoped R1 to the load-bearing `Pi_*` transcendental root only and left the `Pi_match`/`That`/`Sigma0` checks as-is, corroborated by the SymPy F2 block). Flagging for awareness, not blocking.

## Verdict justification

The single finding R1 is fully resolved. Codex applied exactly the prescribed cleared-denominator bracketed `FindRoot` (sign correct, `(Exp[p]-1)`), added a residual-near-zero witness that passes (~0, not the §131-F3 ~6366 signature), and re-anchored the Pi_* check to stage-131's independently-computed owning value `1.50882951349315558300555075595` — eliminating both the old transliterated `FindRoot[gPi==gMinus]` and the old truncated `1.50882951349316` literal. The Mathematica root tracks stage 131 (not 144's own SymPy nsolve, which diverges at digit 16), giving genuine second-engine independence under the MIRROR POLICY. SymPy `.py` was correctly left untouched and still passes. Both engines exit 0; committed outputs are fresh. No regressions, no material change.
