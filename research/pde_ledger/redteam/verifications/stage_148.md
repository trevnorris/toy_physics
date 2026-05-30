---
unit_id: 148
batch: IV.5
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-05-29T23:30:00Z
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 2
findings_total: 2
material_change: true
---

# Verification — unit 148 (remediation pass, supersedes 2026-05-27 v2)

Authoritative finding source = directive `redteam/directives/stage_148.md` (F1, F2 with
`## Applied:` blocks), encoding `redteam/codex_reviews/stage_148.md` and consult batch 7
(Q1, Q2). The prior 2026-05-27 verification recorded 3 findings and flagged the
cross-engine `dT/eps` divergence as a NON-blocking side observation — that was the live
bug; it is the crux this pass confirms FIXED.

## Per-finding outcomes

### F1 — paper_misalignment / insufficient_assertion (ξ_* bridge: radical constant + tolerance)

**Classification:** resolved

**What changed:**
- SymPy py:90-109: the loose `assert abs(residual) < sp.Float("1e-15")` (on a ~7.82e-17
  residual) is replaced by the EXACT symbolic route. `xi_star_closed` uses
  `sqrt(4107 - 100*sp.pi**2)` (py:93). The bias-neutral point is rebuilt from the EXACT
  `gminus`: `rF1_exact` (py:102), `gminus_exact` (py:103),
  `one_minus_lam_exact = (pi/4 - gminus_exact)/(pi/4 - 2/pi)` (py:104), and
  `exact_resid = sp.simplify(one_minus_lam_exact - xi_star_closed)` (py:105). The
  load-bearing assert is `assert exact_resid == 0` (py:109). The numeric `residual`
  (py:107) is retained only for transcript display, NOT asserted.
- Mathematica wl:96-98: no-op confirmation — live line reads `Sqrt[4107 - 100*Pi^2]`
  (wl:96) and keeps the exact `expectZero`/`FullSimplify` zero (wl:98). wl:36
  `MaxIterations -> 100` (FindRoot cap) untouched.

**Assessment:** Correct and non-tautological. No `168` survives in either live script —
the sole `dSigmaOfDeltas` grep hit (wl:46) is a comment naming the OLD route that was
removed, not code. `100` is the rF1-forced constant (`12*(37/20)^2 = 4107/100`). SymPy
output line 30 shows `exact (1-lambda_(Pi,0)) - xi_* = 0` (true symbolic 0); the assert is
on `exact_resid`, not the loose numeric. Mathematica output shows
`(1-lambda_(Pi,0)) - xi_* = 0` / PASS. Same-source consistency check accepted per consult
Q2.

### F2 — insufficient_verification (LIVE cross-engine divergence)

**Classification:** resolved

**What changed:**
- Mathematica wl:43-60: the chain-term-dropping `dSigmaOfDeltas`/`dTOfDeltas` block was
  DELETED (diff -9..-14) and replaced with Mathematica's OWN symbolic differentiation:
  `Tm[pp_] := Sqrt[(9/20)*(pp/(1 - sFormula/4))] /. p -> pp` (wl:47),
  `aT = N[-(D[Tm[p], p] /. p -> pStar)/gPrimeStar, 30]` (wl:49) — a genuine `D[]` total
  p-derivative along the S=sFormula(p) curve, regenerating the `sFormula'` chain term.
  `bT` is the explicit fixed-Π S-sensitivity (wl:51). Paper-literal anchors at wl:57-60
  (`aT vs paper literal A_T` / `bT vs paper literal B_T`, vs verbatim part04:846/848,
  1e-11 tol).
- SymPy py:42-45: added EXTERNAL paper-literal asserts
  `abs(sp.N(AT,30) - sp.Float("-4.27263956256927")) < 1e-11` and the BT analogue. SymPy's
  `AT` (py:32-37) retains the chain term `Pi_star*Sp_star/(4*gp_star*(1-S_star/4)**2)`.

**Assessment:** The live bug is fixed.
- (3) Mathematica `aT`/`bT` come from its own `D[Tm[p], p]`, NOT a transliteration of
  SymPy's `AT` algebra — confirmed in the live `.wl`.
- (4) Both engines anchor to paper literals and PASS: Mathematica output lines 5-8
  (`aT vs paper literal A_T = 0`/PASS, `bT vs paper literal B_T = 0`/PASS); SymPy exit 0
  with the asserts present, no AssertionError. The literals match
  paper/appendices/stage_appendix_part04.tex:846,848 verbatim.
- (5) CROSS-ENGINE AGREEMENT FIXED. SymPy `dTU = 0.508756302215083911...` vs Mathematica
  `dTU = 0.5087563022150839246...` (agree ~16 digits); SymPy
  `dTD = -0.116943802151810766...` vs Mathematica `dTD = -0.1169438021518107487...`
  (~16 digits). The OLD divergent `0.4976692...`/`-0.1144451...` are GONE (grep finds no
  `0.4976`/`0.1144` in either script or output).
- (6) ANTI-FRAGILITY: NO baked SymPy literal in the `.wl`. grep for
  `sympyDTU`/`sympyDTD`/`0.5087`/`0.4976` in the `.wl` returns nothing; the only
  cross-engine literals are the EXTERNAL paper values, exactly as the directive mandated.
- (7) No double-counting (consult Q1b): the `D[]` derivation is independent, anchors
  external. A_T carries the Sp chain term, B_T multiplies the projected S-deviation in the
  same validated model. The fix ADDED the chain term to Mathematica — it did not remove
  it.

## Exec log assessment

**SymPy:** exit=0. Notable lines: `uniform: dT/eps = 0.508756302215083911...`;
`derivative: dT/eps = -0.116943802151810766...`; `exact (1-lambda_(Pi,0)) - xi_* = 0`;
`(1-lambda_(Pi,0)) - xi_* = 7.82e-17` (display only). No AssertionError → AT/BT
paper-literal asserts passed.

**Mathematica:** exit=0. Notable lines: `aT vs paper literal A_T = 0` / `PASS`;
`bT vs paper literal B_T = 0` / `PASS`; `uniform: dT/eps = 0.5087563022150839246...`;
`derivative: dT/eps = -0.1169438021518107487...`; `(1-lambda_(Pi,0)) - xi_* = 0` / `PASS`.

**Output freshness:** confirmed. Output txts regenerated 2026-05-29 23:02:08; scripts last
modified 22:52:44 (.py) / 22:53:25 (.wl) — outputs newer than scripts.

## Material-change assessment

`material_change`: true. The Mathematica `dTU`/`dTD`/`dTLam` numerics changed (corrected to
match SymPy and the paper). SymPy values were already correct, so downstream consumers that
trusted SymPy are unaffected; the change repairs the Mathematica mirror. Downstream units
> 148 may be marked upstream_stale, but no SymPy-derived numeric shifted, so a narrow
re-check (Mathematica-mirror only) suffices for any unit that depended on this stage.

## Side observations (non-blocking)

- SymPy transcript still prints the ~7.82e-17 numeric residual line; intentional
  display-only context, no longer load-bearing. Not a concern.
- wl:46 comment references `dSigmaOfDeltas` solely to document the removed route.

## Verdict justification

Both findings resolved. F1: the assert is now an exact symbolic `exact_resid == 0` built
from the exact `gminus`, with the rF1-forced `100*pi^2` radical and no surviving `168`.
F2: the live divergence is fixed — Mathematica derives `aT`/`bT` via its own `D[]` autodiff
(restoring the dropped chain term), both engines anchor to the external paper literals
A_T/B_T (all PASS), no baked SymPy literal was introduced, and the printed `dTU`/`dTD` now
agree across engines to ~16 digits with the old divergent values gone. Both scripts exit 0
with non-tautological assertions; outputs are fresh.
