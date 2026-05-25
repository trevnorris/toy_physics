---
unit_id: 081
batch: III.4
verifier_model: claude-opus-4-7[1m]
verify_date: 2026-05-25T00:40:00Z
verdict: verified
sympy_exit: n/a
mathematica_exit: 0
findings_resolved: 4
findings_total: 4
material_change: false
---

# Verification — unit 081

## Per-finding outcomes

### F1 — hardcoded_result

**Classification:** resolved

**What changed:**
`mathematica/moving_throat_pde_stage081_family1_pi_thresholds_mathematica_audit.wl:52-55` — `qq` is now `FullSimplify[piOfZeta / cMix, Assumptions -> $Assumptions]` (then a `ConditionalExpression` strip on :53), and the new `expectZero["Q matches closed form", qq - (1 + zeta - 2*epsBlk*zeta)/(1 - epsBlk*zeta)]` ties the derived `qq` back to the closed form. The Solve-derived `piOfZeta` (line 46) is now consumed.

**Assessment:**
Matches the directive verbatim modulo the orchestrator's pre-disclosed cosmetic edits (the `ConditionalExpression[e_, _] :> e` strip applied to `piOfZeta`, to `qq`, and inside `expectZero`). The residual `qq - (1 + zeta - 2*epsBlk*zeta)/(1 - epsBlk*zeta)` evaluates to plain `0` in the output (line 5: `Q matches closed form = 0`), which means the Solve-derived inversion and the previously hand-encoded closed form are now provably the same rational function on the declared domain. This is genuinely non-tautological: if `Solve[zeta == zetaExpr, piTr]` returned a different root or the premise `zetaExpr` were retyped, this check would fire.

### F2 — tautological_check

**Classification:** resolved

**What changed:**
Both the `dqq = FullSimplify[D[qq, zeta], ...]` definition (originally line 42) and the `expectZero["dQ/dzeta exact formula", dqq - (1 - epsBlk)/(1 - epsBlk*zeta)^2]` assertion (originally line 48) are deleted. `grep` for `dqq` and `dQ/dzeta` in both the script and its output returns nothing.

**Assessment:**
Matches the directive exactly. No `dqq`/`dQ/dzeta` lines remain in the transcript. No side effects: nothing else referenced `dqq`.

### F3 — tautological_check

**Classification:** resolved

**What changed:**
Lines 80-84 of the post-edit script replace each `expectApprox[..., N[piXxxOverC /. epsBlk -> 0, 40], ToExpression["3.46…`25"], 10^-14]` with `expectApprox["… matches 1+zeta", N[(piXxxOverC - (1 + zetaXxx)) /. epsBlk -> 0, 40], 0, 10^-14]`.

**Assessment:**
Matches the directive. The substantive content of each check is now "the symbolic value `qq /. zeta -> zeta_* /. epsBlk -> 0` minus `1 + zeta_*` equals 0 numerically", which exercises `qq`'s functional form rather than comparing a constant to its typed expansion. Output (lines 18-27) shows five `PASS: …matches 1+zeta` lines with `diff = 0``…`. After F1, `qq` is derived from `piOfZeta`, so these checks now do propagate failures in the inversion. Non-tautological.

### F4 — tautological_check

**Classification:** resolved

**What changed:**
Lines 86-88 replace the literal `expectApprox["blocking ceiling numeric check", epsCeiling, ToExpression["0.40526368971137149977`25"], 10^-14]` with `expectApprox["blocking ceiling reciprocal", N[epsCeiling*zetaMaxF1 - 1, 40], 0, 10^-14]`.

**Assessment:**
Matches the directive exactly. The transcript (line 29) shows `blocking ceiling reciprocal diff = 0``19.69…` then PASS. The new check verifies the reciprocal identity `epsCeiling * zetaMaxF1 == 1`. While at the symbolic level this is the definition of `epsCeiling`, the surrounding numeric `N[..., 40]` round-trip turns this into a meaningful guard against `N[]` precision loss / wrong `zetaMaxF1` retyping — and is exactly what the directive prescribed.

## Exec log assessment

**SymPy:** exit=n/a. No SymPy script was modified for unit 081 (all findings were against the Mathematica mirror), and no `stage_081_sympy.log` is present in `redteam/exec_logs/`. The pre-existing SymPy output (`scripts/output/moving_throat_pde_stage081_family1_pi_thresholds_sympy_audit.txt`, mtime 2026-05-22) is still the canonical reference and matches both engines' closed form and numerics.

**Mathematica:** exit=0. The orchestrator-supplied log file `redteam/exec_logs/stage_081_mathematica.log` is missing (only `stage_081_diff.patch` is present in exec_logs), but the canonical refreshed transcript `mathematica/output/moving_throat_pde_stage081_family1_pi_thresholds_mathematica_audit.txt` ends with `Stage 081 Mathematica audit passed.` and shows every expected PASS line:

- `PASS: Q matches closed form` (new from F1)
- `PASS: Q(0)-1`, `PASS: Q(1)-2` (still pass against Solve-derived `qq`)
- Five `PASS: … matches 1+zeta` lines (new from F3, replacing magic-number checks)
- `PASS: blocking ceiling reciprocal` (new from F4, replacing magic-number check)
- No `dQ/dzeta exact formula` line (removed per F2)

The `Pi_of_zeta` symbolic print on line 7 of the output now shows `cMix*(1 + zeta - 2*epsBlk*zeta)/(1 - epsBlk*zeta)` — i.e. the Solve branch is the expected one, confirming the cross-engine inversion match. Q symbolic form on line 8 matches both engines.

**Output freshness:** `mathematica/output/...stage081...txt` mtime is 2026-05-25 00:21:34, post-dating the .wl mtime of 2026-05-23 10:36:58. Confirmed refreshed.

## Material-change assessment

`material_change`: false.

The Mathematica-side closed-form `Q(zeta;eps_blk) = (1 + zeta - 2*epsBlk*zeta)/(1 - epsBlk*zeta)` is unchanged; only the route by which it is obtained changed (derived from `piOfZeta/cMix` instead of hand-encoded). All printed numeric values (`Pi_xxx/cMix`, blocking ceiling `0.40526368971137148…`) are unchanged from the prior version of the script, and they still agree with the SymPy output. No downstream-visible derivation result moved.

## Side observations (non-blocking)

1. The `Pi_xxx^()/C_mix` symbolic prints (output lines 13-17) are now partial-fraction-decomposed rationals in `epsBlk` (e.g. `2 + 0.594…/(0.405… - 1.*epsBlk)`) rather than the SymPy-style `(a*epsBlk - b)/(c*epsBlk - 1)` ratio. This is `FullSimplify` choosing a different normal form on the Solve-derived `qq`. At `epsBlk -> 0` they reduce to the same `1 + zeta_*` numbers (verified by the F3 residual checks), so it is purely cosmetic — but a reader comparing SymPy/Mathematica outputs side-by-side will need to evaluate at `epsBlk = 0` to see the agreement. Not a finding.

2. Orchestrator's manual ConditionalExpression-strip in `expectZero` (the wrapper helper) is generic and will apply to all future stages using this helper. The strip is sound under `$Assumptions` because `ConditionalExpression[0, cond]` is identically zero on the declared domain — but if a future audit places an assertion outside the declared domain, the strip could mask a genuine failure. Out of scope for this verification.

## Verdict justification

All four findings are resolved exactly per directive: `qq` is now derived from the Solve-based `piOfZeta` and tied to the closed form by a new residual assertion (F1); the tautological derivative check is removed (F2); the five magic-number numeric comparisons are replaced with residuals against the `1 + zeta` functional form (F3); and the blocking-ceiling magic-number check is replaced with the reciprocal identity (F4). The refreshed Mathematica transcript shows every expected PASS line and exits 0. The orchestrator's pre-disclosed manual edit (the `ConditionalExpression` wrapper stripping) is cosmetic — it lets the new F1 residual print as bare `0` instead of `ConditionalExpression[0, cond]` and does not alter any numeric or symbolic result. No SymPy script was touched; the cross-engine agreement still holds at the closed-form and numeric levels. No regressions in the diff.

stage 081: verified
