---
unit_id: 250
batch: VIII.1
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-06-03T13:20:00-06:00
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 2
findings_total: 2
material_change: false
---

# Verification — unit 250

## Per-finding outcomes

### F1 — missing_verification_script (missing_mathematica)

**Classification:** resolved

**What changed:**
A new independent-route audit was created at
`mathematica/moving_throat_pde_stage250_crossing_vs_collapse_goldilocks_window_compiler_from_the_stage248_event_chain_and_relaxed_wall_timescale_closure_mathematica_audit.wl`
(6725 bytes, mtime 12:59). MANIFEST `mathematica.path` is now populated (line 8950, previously null). It covers M1–M7 with project-standard helpers (`expectZero`/`expectTrue`/`expectApprox`, `stripConditional`).

**Assessment:**
Genuinely independent, NOT a transliteration of the `.py`. The decomposition is different and uses native universal-quantifier proofs where the SymPy script used sampling:
- M1 (.wl:74-93): `Resolve[ForAll[..., Implies[ms>0 && lambdaEff>0 && En>Vpeak, D[tcross,En] < 0]]]` over `Reals` — proves strict negativity over the WHOLE domain (both an E-form and an x-form), returning literal `True`.
- M3 (.wl:107-124): edge via both `Solve` AND `Reduce`, plus a `Resolve[ForAll]` uniqueness statement — different route than the `.py`'s `solve(Eq(S**2,1),E)[0]`.
- M4 (.wl:133-145, the headline): `Resolve[ForAll[..., Implies[positivity && En>Vpeak, Equivalent[Sratio<1, En>Eedge]]]]` returns `True`. This is a genuine GLOBAL strict biconditional `S(E)<1 ⟺ E>E_edge` over the full half-line — proven by `Resolve`, NOT a sample point, and both sides are strict. This is the strongest independent form of the one-sided/unique-edge claim and is exactly what the SymPy side lacked.
- M5 (.wl:149-155): `EedgeHeavy = Eedge /. muEta->alpha ms` genuinely cancels `ms` (Eedge's `ms/(2 muEta)` becomes `1/(2 alpha)`), so `D[EedgeHeavy, ms]==0` is non-vacuous; the post-cancellation value is independently confirmed.
- M6/M7: speed-space identity via `FullSimplify` to 0, and the eight Session-III benchmark `expectApprox`/`expectTrue` regressions, all matching the published literals.
The log shows all `expectTrue` checks reduced to literal `True` (otherwise `fail[]` would `Exit[1]`; exit is 0), and ~18 PASS lines cover M1–M7. No non-strict testing of the strict window inclusion; no sample-point substitute for the global statement.

### F2 — insufficient_verification

**Classification:** resolved

**What changed:**
Two new global strict-monotonicity assertions added to the `.py` (diff lines 9-16 and 25-44):
- §2 (.py:63-70): positive gap symbol `x` declared (`positive=True`); `dtcross_dE_gap = dtcross_dE.subs(E, Vpeak+x)`; an identity assert against the closed negative form, then `assert (dtcross_dE_gap < 0) == True`.
- §4 (.py:93, 104-112): `dS_dE_gap = dS_dE.subs(E, Vpeak+x)`; identity assert against the closed negative form, then `assert (dS_dE_gap < 0) == True`. A `dS/dE (x=E-Vpeak)` print was added.

**Assessment:**
Global strict monotonicity is now GENUINELY ASSERTED, not sampled. With `x` declared `positive=True` and all params (`chi_peak, gUV, lam_eff, m_s, mu_eta`) `positive=True`, SymPy's assumption engine evaluates `(dS_dE_gap < 0)` to `S.true` for ALL positive `x` (the expression is a negative constant times sign-definite factors). The `== True` guard is non-trivial: had SymPy left the relational unevaluated, `Relational == True` would be `False` and the assert would FAIL — so the passing exit-0 run confirms a real global reduction, not a vacuous pass. The check is falsifiable: a non-monotone `S(E)` (derivative changing sign) would not satisfy `(expr<0)==True`. This closes the gap where the old script relied only on `S_inf==0` (a limit) plus the single-point `Smax_num<1.0`. The `dt_cross/dE<0` companion (notes §1.2) is asserted identically in §2. No existing assert was removed or weakened and no benchmark literal was touched (diff is additive only).

## Exec log assessment

**SymPy:** exit=0. Notable lines:
- `dt_cross/dE (x=E-Vpeak) = -sqrt(2)*lambda_eff*sqrt(m_s)/(4*x**(3/2))` (gap form recorded)
- `dS/dE (x=E-Vpeak)     = -sqrt(2)*sqrt(chi_peak)*sqrt(g_UV)*lambda_eff*sqrt(m_s)/(4*sqrt(mu_eta)*x**(3/2))`
- `All symbolic and numerical checks passed.`

**Mathematica:** exit=0. Notable lines:
- `M1 global dt_cross/dE < 0 = True` / `M1 global dt_cross/dE < 0 for x>0 = True` (PASS)
- `M4 S(E)<1 iff E>E_edge globally = True` (PASS) — the headline global window proof
- `M5 dE_edge/dm_s after mu_eta=alpha m_s = 0` (PASS) — cancellation non-vacuous
- `All Stage 250 Mathematica checks passed.`

**Output freshness:** confirmed. `scripts/output/...sympy_audit.txt` (13:05:18) and `mathematica/output/...mathematica_audit.txt` (13:05:18) are both newer than their scripts (`.py` 12:54, `.wl` 12:59). Outputs were regenerated post-fix.

## Material-change assessment

`material_change`: false.

No derived result changed: the edge formula `E_edge = Vpeak + (m_s/mu_eta) g_UV chi_peak lambda_eff^2/2`, the speed-space identity, and every Session-III benchmark literal are identical to pre-fix. F1 added a second engine; F2 strengthened existing verification from sampled to global. Nothing downstream of 250 depends on a newly-changed constant.

## Side observations (non-blocking)

- MANIFEST `mathematica_output.path` for unit 250 is still `null` (line 8958) even though the output `.txt` exists on disk and is fresh. This is a MANIFEST bookkeeping nit, not a verification gap — the engine output was captured and the exec log exit is 0. Flagging for the orchestrator to backfill if desired; does not affect this verdict.

## Verdict justification

Both findings are resolved. F1's new `.wl` is a genuinely independent native route: it proves the headline one-sided window `S(E)<1 ⟺ E>E_edge` and the strict monotonicities via `Resolve[ForAll]` over the full domain (returning literal `True`), derives the edge via `Solve`+`Reduce` with a `Resolve` uniqueness proof, and confirms the heavy-throat cancellation non-vacuously — a different decomposition than the SymPy `subs`/`simplify` choreography, with no non-strict or sample-point substitute for the global statement. F2 now asserts `dS/dE<0` and `dt_cross/dE<0` globally and strictly (positive gap variable + positivity assumptions, `(expr<0)==True` falsifiable), replacing the prior limit-plus-single-sample. Both engines exit 0, outputs are fresh, the diff is purely additive (no literal or existing assert weakened), and `material_change` is false.
