---
unit_id: 243
batch: VIII.1
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-06-03T10:18:00-06:00
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 2
findings_total: 2
material_change: false
---

# Verification — unit 243

CHECKPOINT stage — higher bar applied: I confirmed the re-authored `.wl` is a
genuinely independent second-engine route (not a re-skinned transliteration),
that it no longer relies on the SymPy script's hardcoded `expected*` closed forms
as the sole reference, and that each load-bearing quantity is independently
derived/witnessed. Both engines exit 0 (SymPy + 43 Mathematica `expectZero`
checks, all PASS).

## Per-finding outcomes

### F1 — mathematica_transliteration

**Classification:** resolved

**What changed:**
The `.wl` (mathematica/...stage243..._mathematica_audit.wl) was re-authored block
by block, replacing the line-by-line port of the `.py`:

- Block I (wl:97-106): dropped `expectedSleak`/`expectedWwork` as the sole
  reference. Added an IBP-closure witness `Sleak + ibpInterior - boundary == 0`
  (wl:101), where `ibpInterior = Integrate[W D[jw,w], ...]` (wl:89-92) is computed
  independently of `Sleak = Integrate[D[W,w] jw, ...]`. Added a half-line
  even-integrand decomposition `Wwork - 2*Integrate[jw Ew,{w,0,Inf}]` (wl:103).
  The closed-form residuals remain only as *additional* checks (wl:102,104), which
  the directive explicitly permitted.
- Block II (wl:115-139): removed hardcoded `Uexpected`/`Vexpected`/`drainExpected`.
  Added back-substitution `stationarity /. uvSol -> {0,0}` (wl:131) and a native
  `LinearSolve[{{kU,-chiLam},{-chiLam,kV}},{fU,0}]` agreement check (wl:132). Drain
  reference rebuilt from the in-script `detH` (`chiLam^2 kV fU^2/detH^2`, wl:126-129,135).
- Block III (wl:148,156-161): removed the manual `{Cos[Pi z]->y, Cos[2 Pi z]->2y^2-1}`
  substitution; rewrite now confirmed functionally via `TrigExpand[varsigma -
  (varsigmaY /. y->Cos[Pi z])]` (wl:156-159), forcing Mathematica's own double-angle
  identity. Added `vertex is stationary`: `D[varsigmaY,y] /. y->yStar == 0` (wl:161).
- Block V (wl:213-251): added an asymptotic-series route
  (`Limit[Normal[Series[x*mode,{x,Inf,0}]], x->Inf]`) alongside the existing
  `Limit` checks.
- Banner (wl:59): `STAGE 226` -> `STAGE 243`.

**Assessment:**
Correct and genuinely independent. The IBP closure is the true product-rule
identity `∫W'jw + ∫W jw' = [W jw]`, cross-validating `Sleak` against an
independently-computed interior integral — not against a hardcoded constant.
LinearSolve is a different native primitive from the symbolic `Solve`. The
functional trig rewrite puts the double-angle burden on `TrigExpand` rather than a
hand-typed `.subs` mirror. The series route is a distinct asymptotic mechanism. The
hardcoded `expectedSleak`/`Uexpected`/`Vexpected`/`drainExpected` literals are
confirmed removed in the diff (patch lines 19,46-47,59 deleted). No tautology traps:
the back-substitution feeds the solution into the original gradient equations (the
load-bearing premises), not into the equation that produced it as `X-X`; the
`vertex is stationary` derivative is w.r.t. `y`, which `varsigmaY` actually contains
(no variable-independence trap). Codex's one documented deviation (completing the IBP
identity with the interior term) is mathematically correct and strengthens, not
weakens, the check.

### F2 — stale_output (banner label artifact)

**Classification:** resolved

**What changed:** wl:59 banner string changed `STAGE 226` -> `STAGE 243`
(patch lines 9-10). Regenerated output (mathematica output line 8) prints
`STAGE 243 — RELAXED-CONSTRAINT BRANCH DECLARATION AND SHORT-RANGE OPEN-SYSTEM COMPILER`.

**Assessment:** Correct and non-load-bearing; banner now matches the stage.

## Exec log assessment

**SymPy:** exit=0. Notable lines: `S_leak = -sqrt(2)*ell_w*j0/4` /
`Expected S_leak = -sqrt(2)*ell_w*j0/4`; `Integral over [0,1] = 1`;
`lim x*deltaV_stat = 0`. SymPy script unchanged (per directive: it was correct and
paper-aligned), values identical to pre-fix.

**Mathematica:** exit=0. All 43 `expectZero` checks PASS. Notable independent
witnesses: `IBP closure = 0 / PASS`; `work half-line symmetry = 0 / PASS`;
`U/V satisfy stationarity = MatrixForm[{{0,0}}] / PASS`; `U/V equal LinearSolve =
MatrixForm[{{0,0}}] / PASS`; `quadratic rewrite (functional) = 0 / PASS`;
`vertex is stationary = 0 / PASS`; `x QQ -> 0 (series) = 0 / PASS`. Count reconciles:
44 `expectZero[` tokens in the file minus 1 (the function definition at wl:48) = 43
invocations = 43 PASS lines; none failed (a `fail` would `Exit[1]`).

**Output freshness:** confirmed. `.wl` mtime 09:47:59, sympy `.py` mtime 05-11;
both regenerated `.txt` outputs mtime 09:54:22 — newer than their scripts. Fresh.

## Material-change assessment

`material_change`: false.

This was a transliteration -> independent-route rewrite of the SECOND engine only.
Every derived/carried numeric value (S_leak, W_w, U, V, drain, det H, V/U, vertex,
boundary values, the kernel span, and all limits) lands on exactly the same closed
forms as before — now reached by independent mechanisms rather than echoing the
`.py`. The SymPy script was not touched. No downstream-visible value changed; no
unit > 243 needs re-audit on account of this fix.

## Side observations (non-blocking)

- `interior stationary point` (wl:160) is `yStar + a/(4 b)` with
  `yStar = FullSimplify[-a/(4 b)]`, i.e. it confirms `yStar` equals its own
  definition (near-vacuous). This is unchanged pre-existing infrastructure and is
  no longer the load-bearing stationarity check — the new `vertex is stationary`
  (wl:161, derivative vanishes) carries that weight. Not a regression; noted only
  for completeness.
- The `QQ/QY/YY source product` checks (wl:239-241) confirm the engine's arithmetic
  on the kernel definitions `SQ=x^-3`, `SY=Exp[-2 kappa x]/x` (e.g. `SQ^2 == x^-6`),
  mild self-confirmation but matches the audited-accepted structure and the SymPy
  side; substantive M5 independence comes from the limit + series routes.

## Verdict justification

Both findings are resolved. F1's higher checkpoint bar is genuinely met: the `.wl`
now confirms each load-bearing quantity through a decomposition distinct from the
`.py` (IBP closure, half-line symmetry, stationarity back-substitution + native
LinearSolve, functional TrigExpand rewrite, stationary-vertex derivative, and
asymptotic-series limits), and the shared hardcoded `expected*` closed forms that
previously made the second engine a structural echo are removed — a wrong constant
would now be caught by the independent witnesses on the Mathematica side. No
tautology or variable-independence trap survives. F2 banner is fixed and reflected in
the fresh output. Exec logs pass on both engines with reconciled PASS counts and
fresh outputs. No regression in the diff. material_change is false (same values,
independent route). Verdict: verified.
