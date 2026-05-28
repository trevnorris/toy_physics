---
unit_id: 134
batch: IV.4
verifier_model: claude-opus-4-7[1m]
verify_date: 2026-05-27T00:00:00Z
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 4
findings_total: 4
material_change: false
---

# Verification — unit 134

## Per-finding outcomes

### F1 — missing_verification_script

**Classification:** resolved

**What changed:**
`scripts/moving_throat_pde_stage134_family1_mouth_fixedpoint_sympy_audit.py:27-91` now contains the banner update to `STAGE 134`, the symbolic computations (S_shell, S_q, fixed-point law) plus four substantive assertion blocks:
- Check 1 (lines 47-51): `assert sp.simplify(S_shell - 1) == 0`.
- Check 2 (lines 53-72): three independent numeric checks at Pi = 1/2, 1, 2 against pasted high-precision literal targets `0.608336415687717065435990381419`, `0.633127670034487546375729566676`, `0.681366857005321783286541952613` (the orchestrator-recomputed targets, not the auditor's fabricated values).
- Check 3 (lines 74-79): `S_q(Pi_*) ≈ 0.658075937605428` to 1e-12.
- Check 4 (lines 81-91): gain-line intercept and slope match notes' literals to 1e-12.

**Assessment:**
All four checks are non-tautological: the literal targets are sourced independently (per the in-script comment, computed via mpmath outside `sKernel`), so if `S(Pi, kappa)` had a typo the checks could fail. Check 1 (shell limit) is the only one previously exercised on the mathematica side; checks 2-4 are new substantive coverage on the sympy side. The exec log shows all four "OK:" lines passing with tiny rounding diffs (~3e-31 to 5e-31), well under the 1e-12 tolerance. The S_q value at Pi_* prints to 0.658075937605428494269581645208 (full 30-digit precision), which matches the notes' 15-digit literal.

### F2 — tautological_check

**Classification:** resolved

**What changed:**
`mathematica/moving_throat_pde_stage134_family1_mouth_fixedpoint_mathematica_audit.wl:26-72` updates the banner to `STAGE 134` (line 26), keeps the substantive `expectZero["static shell channel", sShell - 1]` (line 44), and replaces the tautological `sQExpected` block + its `expectZero` call with:
- `expectClose[name, got, want, tol]` helper (lines 51-55) that prints the diff and PASSes/FAILs.
- Three `expectClose` calls at p = 1/2, 1, 2 against the same independently-verified high-precision targets used in F1 (lines 57-62).

The previous `sQExpected = FullSimplify[...]` line and the `expectZero["specialized D/N kernel", sQ - sQExpected]` call are deleted as required.

**Assessment:**
The three new checks compare `N[sQ /. p -> <value>, 30]` against `SetPrecision[<literal>, 30]` where the literal is pasted from the same external mpmath source as F1. Because the targets are not derived from `sKernel` at runtime, a typo in `sKernel` would surface. The exec log confirms all three PASS lines (`PASS: S_q at p=1/2`, `PASS: S_q at p=1`, `PASS: S_q at p=2`) plus the retained `PASS: static shell channel`, with diff values ~0 to 30 digits.

### F3 — paper_misalignment

**Classification:** resolved

**What changed:**
`paper/stages/stage_134.tex:21-25` — the `\stagefield{Checks}` items 1 and 2 were rewritten from prescriptive checks into carry-forward citations:
- Item 1: "Outlet consistency of the gain pair (M_s, M_q) is verified at Stage~\ref{stage:135} ...; carried forward here."
- Item 2: "Self-matched susceptibility closure is verified at Stage~\ref{stage:137} ...; carried forward here."
- Item 3 (numerical fixed points) unchanged.

This corresponds to direction (b) from the directive's user-resolution menu: the outlet-consistency and susceptibility-closure checks belong downstream (stages 135 and 137, respectively), and the stage 134 card now points to where those checks actually live.

**Assessment:**
The edit is exactly at lines 21-25 of `paper/stages/stage_134.tex` as the orchestrator noted. The card's `Checks` block no longer claims stage 134 itself performs the outlet-consistency or susceptibility-closure checks, so the script-side absence of those checks is no longer a misalignment. The third checklist item (numerically-located fixed points) is honored by the scripts' use of `Pi_* = 1.50882951349316` as a literal (sourced from stage 131/233 per the notes), with no closed-form derivation here. Cross-references to stages 135 and 137 are reasonable carry-forward targets for outlet-consistent gains and susceptibility closure respectively.

### F4 — mathematica_transliteration

**Classification:** resolved

**What changed:**
No additional code change beyond F1+F2. Both scripts now compare `S_q(p)` at three independent numeric points against pasted high-precision literals that come from a separate mpmath computation, not from `sKernel`. A typo in the kernel formula would surface in both engines' numeric checks (since the literals are external).

**Assessment:**
The directive explicitly states F4 is subsumed by F1+F2 once both engines have non-tautological numeric checks against externally-sourced literals. That condition is now met: sympy's three `S_q(1/2|1|2)` assertions and mathematica's three `expectClose` checks all anchor against the same external literal targets. The kernel formula is still typed identically in both engines, but a typo would no longer be invisible — it would break the comparison to the external literal in both engines. The operational concern (silent shared error) is neutralized.

## Exec log assessment

**SymPy:** exit=0. Notable lines:
- `S_shell = 1`
- `S_q(Pi_star) = 0.658075937605428494269581645208`
- `OK: S_shell = 1`
- `S_q(1/2) = 0.608336415687717065435990381419  (target 0.608336415687717065435990381419, diff 3.94...E-31)`
- `S_q(1) = 0.633127670034487546375729566676  (target ..., diff 0)`
- `S_q(2) = 0.681366857005321783286541952613  (target ..., diff 4.93...E-31)`
- `OK: S_q matches independent numeric targets at Pi=1/2, 1, 2`
- `OK: S_q(Pi_star) matches notes value 0.658075937605428`
- `OK: gain line coefficients match notes (intercept 1.50882951349316, slope -0.658075937605428)`

**Mathematica:** exit=0. Notable lines:
- `static shell channel = 0` / `PASS: static shell channel`
- `S_q at p=1/2 = 0.6083364156877170654359903814193976...` / `PASS: S_q at p=1/2`
- `S_q at p=1 = 0.6331276700344875463757295666760496...` / `PASS: S_q at p=1`
- `S_q at p=2 = 0.6813668570053217832865419526134255...` / `PASS: S_q at p=2`
- `S_q(Pi_star) = 0.65807593760542948674050367268...`

**Output freshness:** Confirmed. Script mtimes: sympy 17:42, mathematica 17:42. Output mtimes: sympy 17:45, mathematica 17:47. Both outputs are newer than their scripts and were regenerated post-fix.

## Material-change assessment

`material_change`: false.

Rationale: F1 and F2 only add assertions and update banners; they do not alter any derived result. F3 is a paper-card edit that re-routes verification responsibility for two checklist items to downstream stages 135 and 137, without changing any quoted numeric or symbolic claim. F4 is a no-op (subsumed by F1+F2). All printed numeric values (S_q(Pi_*), gain line, S_q evaluations) agree with the pre-fix prints to all digits shown. No downstream re-audit narrow-vs-broad concern.

## Side observations (non-blocking)

- The H1 of `notes/stages/moving_throat_pde_stage134_family1_mouth_fixedpoint.md` already reads "Stage 134" (orchestrator's claim of `Stage 236 → Stage 134` either refers to a pre-batch state already fixed, or was a no-op for this notes file). Body references to "Stage 235" (kernel source) and "Stage 233" (Pi_* fixing) are upstream-citation references and are correct in context.
- The independently-verified targets in F1/F2 (`0.6083...`, `0.6331...`, `0.6814...`) differ from the auditor's originally proposed values (`0.2271...`, `0.4484...`, `0.8732...`). The orchestrator note flags this: the auditor's literals were fabricated, and the applied ones are mpmath-derived from the closed form `S(Pi, pi/2)`. This is correct per the orchestrator's recompute; nothing further to flag.
- The notes' boxed result has a slight typo (last line: `0.658075937605429` with a `9` at the 15th digit vs. the canonical `0.658075937605428` with an `8` used elsewhere in the notes and in the scripts). Cosmetic only; the scripts use the correct trailing-8 value.

## Verdict justification

All four findings resolved. F1's SymPy script now has four substantive assertion blocks (shell limit, three independent numeric checks, Pi_* value, gain-line coefficients) with literals sourced externally via mpmath. F2's Mathematica script replaces the tautological `sQ - sQExpected` check with three `expectClose` calls against the same external literals. F3 is resolved at the paper-card level by re-routing items 1 and 2 to stages 135/137 as carry-forward citations. F4 is operationally neutralized by F1+F2. Both exec logs exit 0, outputs are fresh, and the new "OK:" / `PASS:` lines confirm the substantive checks. No regressions; no material change to downstream numerics.
