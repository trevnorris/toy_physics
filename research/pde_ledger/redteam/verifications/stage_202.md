---
unit_id: 202
batch: VI.1
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-06-02T00:00:00-06:00
verdict: verified
sympy_exit: n/a
mathematica_exit: 0
findings_resolved: 1
findings_total: 1
material_change: false
---

# Verification — unit 202

## Per-finding outcomes

### F1 — missing_verification_script (subtype: missing_mathematica)

**Classification:** resolved

**What changed:**
Codex authored the new second-engine script
`mathematica/moving_throat_pde_stage202_free_quintuple_target_graph_mathematica_audit.wl`
(253 lines, git-untracked `??`, no prior version — confirmed it is a new file).
The orchestrator's independent re-run
(`redteam/exec_logs/stage_202_mathematica.log`) exits 0 with all in-file checks
PASS. The directive's `## Applied: F1` block records the iter-2 reformulation
(linear log-system solve), deviation "none". The git diff
(`stage_202_diff.patch`) is empty because the file is brand new (untracked), and
the SymPy `.py` was not touched — consistent with the directive's ban on editing
the existing engine.

**Assessment:**
Every manifest item M1–M7 appears as an explicit, non-tautological check tracing
to the claim manifest:

- **I / graph solve** — `depMatrix . {logTGraph,logKetaGraph,logMuGraph} - depRhs == 0`
  confirms the LinearSolve is self-consistent (log:18-20, .wl:142-145).
- **M1** — `logCtr[logTGraph] - Log[Ctrtgt] == 0` (.wl:149-152).
- **M2** — `logEta[logKetaGraph] - Log[epsEtatgt] == 0` (.wl:154-157).
- **M3** — `logCnt[logTGraph,logKetaGraph,logMuGraph] - Log[Cnttgt] == 0` (.wl:159-162).
- **M4** — `FreeQ[muGraphFlat,L] && FreeQ[muGraphFlat,Pi]` True, and
  `D[logMuGraph,L] == 0` (.wl:164-169) — two distinct L/Pi-independence checks.
- **M5** — three graph-error identities as a zero 3-vector
  `{ET - qtr/(1+chi0s), EK + qeta, Emu - (qnt - qeta + Fstar qtr/(1+chi0s))}`
  (.wl:181-188).
- **M6** — 8-component projection log residuals zero, plus threaded equation truth
  (.wl:214-218); `x_proj_can` is built independently from the monomial ratios
  `qtr,qnt,qeta`, NOT from the graph triple.
- **M7** — repair-vector rewrite residuals (8-vector) zero, and reduced-family
  packet `{chiQ-1,ET,EK,Emu}` on the graph + chiQ=1 zero (4-vector) (.wl:220-249).

These are genuine reconstructions: the graph triple is obtained from the
LinearSolve, then re-substituted into the *original* monomial logs (M1-M3), and
the error/repair identities are checked against the independently-formed quotient
chart. No `x==x` tautology. No assertion tests a claim the paper does not make.

**Independent-derivation check (load-bearing):**
CONFIRMED independent. The `.py` writes the closed-form graph triple by hand
(`deltaU_graph = (Ctr_tgt/(...))**(1/(1+chi0))`, `.py`:92-104) and substitutes.
The `.wl` instead builds the log-monomial-match system as an explicit
coefficient matrix and `LinearSolve`s it (`.wl`:112-129):

```
depMatrix = {{1 + chi0s, 0, 0}, {0, 1, 0}, {-Fstar, -1, 1}};
...
{logTGraph, logKetaGraph, logMuGraph} = normalize[LinearSolve[depMatrix, depRhs]];
```

The matrix entries are the monomial *exponents*, not pasted closed forms:
- Row 1 `{1+chi0s,0,0}` is the coefficient of `LT` in
  `logCtr[lt] = (1+deltaUs) logTrackingBase + (1+chi0s) logSplit[lt]` where
  `logSplit[lt] = 2Log[Pi] + lt - 2Log[L] - Log[KU]` → `LT` coeff `(1+chi0s)`;
  Keta/muW absent → 0,0.
- Row 2 `{0,1,0}` is the epsEta equation
  `logEta[lk] = 2Log[cetaU] - Log[KU] - lk` solved for `LK`.
- Row 3 `{-Fstar,-1,1}` is `logCnt`'s `{-Fstar LT, -LK, +LMu}` structure from
  `logCnt[lt,lk,lmu] = 2Log[lamW] + lmu - lk - 2Log[KW] + Estar logDriveBase - Fstar logSplit[lt]`.

The RHS vector (`.wl`:118-127) carries exactly the constant (non-dependent)
terms moved to the other side. So the "linear system" genuinely encodes the
log-monomial-match equations with exponents-as-coefficients — it is NOT vacuous
and NOT a transliteration of the `.py` hand substitution. The matrix is upper/
lower-triangular-solvable and the row-1 zero of the residual check (Section I)
confirms the solve is real. This is a *different decomposition* than the SymPy
route, satisfying the anti-transliteration guard.

**Timeout reformulation (iter-2):**
CONFIRMED addressed. No `Solve`/`Reduce`/`NSolve`/`Roots` anywhere in the `.wl`
(grep clean) — the iter-1 transcendental `Solve` that hit the 600s cap is gone,
replaced by the instant `LinearSolve`. No raw global `FullSimplify`; all
reductions route through `scalarReduce`, which wraps
`Simplify[Expand[PowerExpand[...]]]` in `TimeConstrained[..., 20, $Failed]`
(.wl:27-34) and `expectTrue` wraps its `Simplify` in `TimeConstrained[..., 20]`
(.wl:56-65), so no single reduction can stall the run. The exec log is complete
through the closing banner and exits 0; runtime is far under the cap.

## Exec log assessment

**SymPy:** exit=n/a. Directive forbade touching the `.py`; no sympy re-run for
this unit. The existing `.py` reference is unchanged.

**Mathematica:** exit=0. Notable lines:
```
PASS: graph log-system residual            (MatrixForm[{{0, 0, 0}}])
PASS: M1 tracking monomial log residual    (= 0)
PASS: M3 nontracking monomial log residual (= 0)
PASS: M4 mu_graph is free of L and Pi      (= True)
PASS: M6 projection component log residuals (MatrixForm[{{0,0,0,0,0,0,0,0}}])
PASS: M7 reduced-family packet on the target graph (MatrixForm[{{0,0,0,0}}])
STAGE 202 MATHEMATICA AUDIT PASSED
# exit_code: 0
```

**PASS count:** the run produced **11** PASS lines (not 12 as the orchestrator
note anticipated). The 11 lines fully cover every manifest item: graph
log-system residual (1) + M1/M2/M3 (3) + M4 ×2 (FreeQ and D/dL) + M5 (1) +
M6 ×2 (residuals and threaded equations) + M7 ×2 (repair vector and family
packet) = 11. All of M1–M7 are represented; the 11-vs-12 difference is a
counting expectation in the prompt, not a missing or failed check. Each PASS
prints an explicit zero residual or `True`, none are silent.

**Output freshness:** confirmed. `.wl` mtime 2026-06-02 09:52:12; output
`.txt` mtime 2026-06-02 09:57:09 — output is newer than the script and matches
the exec-log header date (09:56:54). Output `.txt` carries the same 11 PASS
lines and exit 0. Outputs were re-generated post-fix.

## Material-change assessment

`material_change`: false. F1 only ADDS a second-engine `.wl`; it does not modify
the existing SymPy reference, paper, or notes, and does not change any derived
result. The graph triple, error identities, repair vector, and family packet
recomputed in Mathematica agree with the SymPy values (all zero residuals), so
downstream Stages 203–218 see no changed input. No `upstream_stale` concern
beyond the routine bookkeeping.

## Side observations (non-blocking)

- The `.py` banners read "STAGE 185 / ... AUDIT COMPLETE" (internal-renumbering
  residue), already noted by the auditor as cosmetic; the new `.wl` correctly
  banners "STAGE 202". No load-bearing identity differs. Not a finding.
- The orchestrator's "expect 12" note overcounts by one; flagging only so the
  state machine isn't tripped by the 11 PASS lines. Coverage is complete.

## Verdict justification

The single finding F1 (missing Mathematica engine) is fully resolved: a new
independent `.wl` exists at the registered path, exits 0, and verifies all of
M1–M7 with explicit non-tautological zero-residual / truth checks. The
derivation is genuinely independent of the SymPy route — it builds the
log-monomial-match coefficient matrix and `LinearSolve`s it (exponents as
coefficients), rather than re-typing the `.py` closed forms. The iter-2
reformulation removed the transcendental `Solve` that caused the timeout and
caps every remaining reduction at 20s, so the run completes well under the 600s
cap. No regressions, no material change. Verdict: verified.
