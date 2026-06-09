---
unit_id: 200
batch: V.3
verifier_model: Opus 4.8 (1M context) [claude-opus-4-8[1m]]
verify_date: 2026-06-09T14:30:00-06:00
verdict: verified
sympy_exit: n/a
mathematica_exit: 0
findings_resolved: 1
findings_total: 1
material_change: false
---

# Verification — unit 200 (CHECKPOINT)

The directive carries one actionable finding (F1 = re-author the `.wl` to break the
`mathematica_transliteration`); F2 is informational (stale-output, resolved by the orchestrator
re-run, no Codex action). I verified the three independence requirements, the deliverable
preservation, and the clean run — quoting `.wl` vs `.py` line numbers — and hand-checked the
new §I weight-vector arithmetic.

## Per-finding outcomes

### F1 — mathematica_transliteration (load-bearing; re-author the `.wl`)

**Classification:** resolved

**What changed (per the diff `stage_200_diff.patch`):**

- **§I `M_*` (wl:160-169):** the autodiff Jacobian `Mderived = Table[D[qPair[[i]], Dvec[[j]]], …]`
  (old wl:157) is REPLACED by an assembly from five primitive exponent-weight vectors
  (`chiCoreWeights`/`thermalWeights`/`nontrackPrefactorWeights`/`epsWWeights`/`epsEtaWeights`,
  wl:160-164):
  `Mderived = {(1+deltaUs) chiCoreWeights + (1+chi0s) thermalWeights, nontrackPrefactorWeights + eStar epsWWeights - fStar thermalWeights, epsEtaWeights}` (wl:165-169).
  `qPair` (wl:153-157, still `Log` of the `ratioToLogs`-substituted ratios) now survives ONLY as
  the cross-check `qPair - Mderived . Dvec == 0` (wl:181), not as the derivation of `M_*`.
- **§III orbit (wl:235-256):** the hand-written closed forms `KetaOrbit/TOrbit/muOrbit`
  (old wl:223-231, mirroring `.py:286-298`) are REPLACED by a genuine SOLVE: build the three
  log-residuals `Log[monomial(Exp[orbitLog…])/Target]` (wl:236-245), read the system matrix via
  `Coefficient[Expand[…], depOrbitLogs]` (wl:246-249), take the constant offset (wl:250-252), and
  `orbitLogSolution = LinearSolve[orbitCoeffMatrix, -orbitOffset]` (wl:253), then exponentiate
  (wl:254-256). The coefficient matrix printed in the output `{{-1,0,0},{0,1+chi0s,0},{-1,-fStar,1}}`
  confirms a real linear solve (it would be undefined if the forms were posited).
- **§V Packet-A (wl:310-320):** the `Series[…,{eps,0,1}]` linearization (old wl:281-292, mirroring
  `.py:409-420`) is REPLACED by an analytic base-point derivative:
  `chiBase = chiPerturbed /. eps->0` (wl:316) and `chiSlopeAtBase = D[chiPerturbed, eps] /. eps->0`
  (wl:317), assembled as `chiBase - 1 + eps chiSlopeAtBase` (wl:318-320). No `Series` remains.

**Assessment (the three independence requirements):**

(A1) **§I independence — confirmed-independent.** `git`-diff confirms `Table[D[qPair…, Dvec…]]` is
gone; there is no `D[]` autodiff of the `ratioToLogs` logs feeding `M_*`. The `.py` derives `M_*`
by `Mderived = q_pair.jacobian(Dvec)` (`.py:154`) — a symbolic differentiation of the same logs.
The `.wl` now derives it from primitive exponent weights — a genuinely different operation (no
differentiation, a linear combination of integer/parameter weight vectors). **Self-test of the
weight arithmetic** (so this is a real reconstruction, not a hidden retyping of `Mexpected`):
- Row 1: `(1+δU)·{0,1,1,-1,0,0,0,0} + (1+χ0)·{0,0,0,-1,0,0,0,1}`
  `= {0, 1+δU, 1+δU, -(1+δU)-(1+χ0), 0,0,0, 1+χ0} = {0, 1+δU, 1+δU, -(2+χ0+δU), 0,0,0, 1+χ0}`
  = `Mexpected` row 1 (wl:171). ✓
- Row 2: `{2,0,0,0,-1,-2,1,0} + E·{2,0,2,-1,0,-1,0,0} - F·{0,0,0,-1,0,0,0,1}`
  `= {2+2E, 0, 2E, -E+F, -1, -2-E, 1, -F} = {2(1+E), 0, 2E, F-E, -1, -(2+E), 1, -F}`
  = `Mexpected` row 2 (wl:172). ✓
- Row 3: `{0,2,0,-1,-1,0,0,0}` = `Mexpected` row 3 (wl:173). ✓
The weight vectors are distinct objects from the `Mexpected` literal, and the algebraic
combination is required to land on each row — a wrong weight or coefficient would fail
`Mderived - Mexpected == 0` (wl:180). Genuine reconstruction, not retyping.

(A2) **§III orbit independence — confirmed-independent.** `KetaOrbit/TOrbit/muOrbit` are now
produced by `Coefficient`→`LinearSolve` of the log-linear residual system (wl:246-256), NOT
posited as the `.py`'s closed forms (`.py:286-298`). This is the same solve-vs-posit pattern the
directive cites for verified siblings 198/199. The printed log-solution (output lines 45-46) is the
solved expression, and the downstream chart checks (wl:271-276) consume it.

(A3) **§V Packet-A independence — confirmed-independent.** Uses base-point derivative
`D[chiPerturbed, eps] /. eps -> 0` (wl:317), NOT `Series`. The `.py` uses
`sp.series(…, eps, 0, 2).removeO()` (`.py:409-420`). Different operation.

§II (witness invariance) and §IV (cocycle) are untouched and inherit independence: §II is built
from the shared monomial premises and orbit law (no extracted method changed), and §IV consumes the
independently-assembled `Mderived` (wl:301-303), so both ride on §I's new route. The shared
monomial DEFINITIONS (`ctrMonomial`/`cntMonomial`/`epsEtaMonomial`, wl:57-66) are the premise the
directive explicitly permits to stay.

**Deliverables preserved — confirmed byte-identical in value:**
- `Mexpected` literal matrix (wl:170-174) — unchanged by the diff (the diff touches only the lines
  that *produce* `Mderived`, not `Mexpected`).
- §I assertion `M_* - Mexpected == 0` (wl:180) — preserved; output line 19-20 reads all-zero.
- §II witness-invariance zeros (wl:212-214) and `DeltaHPairWitness - DeltaHIntrinsic` (wl:231) —
  preserved; output lines 27-38 all `= 0`/`PASS`.
- §III mismatch chart `qMismatchExpected = {(1+chi0s)Log[mT], Log[mMu]-Log[mK]-fStar Log[mT], -Log[mK]}`
  (wl:280) and `M_* . Dmis - qMismatch == 0` (wl:285) — preserved; output lines 53-58 zero/`PASS`.
- §IV cocycle `q31 - q32 - q21 == 0` (wl:306) — preserved; output lines 63-66 zero.
- §V Packet-A coeff `DeltaQLinearExpected = eps(5 epsBeta + dSigma0/(3 S) + 9 dSigma5/S)` (wl:321)
  and assertion (wl:323-326) — preserved; output lines 71-72 `= 0`/`PASS`.
No checkpoint constant moved.

## Exec log assessment

**SymPy:** exit=n/a — the `.py` is the unchanged reference (not edited; the directive forbids
touching it). Its refreshed output (`scripts/output/…sympy_audit.txt`, mtime 14:09) is the engine
cross-check and shows every section zero/PASS with the `STAGE 200` banner.

**Mathematica:** exit=0 (implied: the script ends `Exit[0]`, every `expectZero` would `Exit[1]` on
failure, and the orchestrator's re-run produced a complete output through the `STAGE 200 LEDGER`
footer). Notable lines from `mathematica/output/…mathematica_audit.txt`:
- L19-22: `derived M_* - carried Stage 192 matrix = {{0,…},{0,…},{0,…}}`; `q^(2<-1) - M_* Δx = {0,0,0}`
- L43-44: `dependent orbit log coefficient matrix = {{-1,0,0},{0,1+chi0s,0},{-1,-fStar,1}}` (the
  real `LinearSolve` system — confirms §III is solved, not posited)
- L47-52: three §III chart `… = 0` / `PASS` lines
- L71-72: `Delta_Q^lin - eps(5 eps_beta + dSigma0/(3S) + 9 dSigma5/S) = 0` / `PASS`
Every `expectZero`/`PASS` line passes; no `FAIL`.

**Output freshness / banner:** both `.txt` mtimes are 2026-06-09 14:09, newer than the `.wl` (13:50)
and the `.py` (06-03). The F1-stale `STAGE 183` banner is GONE — both outputs read
`STAGE 200 — EXACT REFERENCE-FREE HOME-STRETCH THEOREM` (line 3) and `STAGE 200 LEDGER` (mathematica
line 77 / sympy line 259). F2 (stale_output) is therefore resolved by the re-run.

## Material-change assessment

`material_change`: false. Only the METHOD that produces `M_*`, the orbit triple, and the Packet-A
linear law changed; every emitted value is byte-identical to before (the `Mexpected` literal, the
three chart-conversion laws, the cocycle, and the Packet-A coefficient all still read `= 0`/`PASS`).
No downstream unit's carried value moves.

## Side observations (non-blocking)

- The `.py` ledger footer (sympy output L263/L274) still reads "logarithmic Jacobian", matching the
  unchanged `.py` route; the `.wl` footer (mathematica output L81/L92) correctly reads "primitive
  exponent ledger". This is the intended asymmetry (the two engines now use different routes), not a
  defect — the directive forbids touching the `.py`.
- §II's `PhiMuW`/`PhiTW` witness factors and §IV's cocycle remain structurally parallel to the
  `.py`, but the directive scopes them as downstream consequences of §I/§III; with §I/§III/§V now
  independent, this parallel structure is the intended second-engine shape, not a surviving port.

## Verdict justification

The re-author lands all three independence requirements with no retyping: §I assembles `M_*` from
primitive exponent-weight vectors (hand-verified to reproduce all three `Mexpected` rows, and the
old `Table[D[…]]` autodiff is confirmed deleted), §III solves the orbit triple via
`Coefficient`→`LinearSolve` (with a printed real coefficient matrix), and §V uses an analytic
base-point `D[…,eps]/.eps->0` instead of `Series`. Every checkpoint deliverable is preserved at its
original value (`Mexpected`, the three chart-conversion laws, the cocycle, the Packet-A coefficient),
the refreshed Mathematica output passes every `expectZero`/`PASS` with the correct `STAGE 200`
banner (stale `STAGE 183` gone), and both `.txt` mtimes postdate the scripts. Unlike the V.1-173
precedent (where a re-author remained insufficient), this re-author is genuinely independent of the
`.py`'s Jacobian/posit/Series route. Verdict: verified.
