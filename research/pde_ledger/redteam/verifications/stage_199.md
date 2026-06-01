---
unit_id: 199
batch: V.3
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-06-01T12:30:00-06:00
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 1
findings_total: 1
material_change: false
---

# Verification — unit 199

## Per-finding outcomes

### F1 — missing_mathematica

**Classification:** resolved

**What changed:**
Codex created a new independent-route Mathematica audit at
`mathematica/moving_throat_pde_stage199_pairwise_orbit_transport_law_mathematica_audit.wl`
(11,768 bytes, mtime 2026-06-01 12:14). The SymPy reference script
`scripts/moving_throat_pde_stage199_pairwise_orbit_transport_law_sympy_audit.py`
was NOT modified (mtime unchanged at 2026-05-11 11:58). The `.wl` carries 42 `PASS:` lines
across eight sections (I monomial compiler, II transport-factor Solve, III mismatch collapse,
IV q-chart + projector split, V restoration map, VI cocycles, VII orbit-lock, VIII fixed-free
reduction) and 0 `FAIL` lines; the committed output reports `# exit_code: 0`.

**Assessment:**

*Independence (genuine route, NOT transliteration).* The `.wl` re-derives the same conclusions
by a materially different decomposition than the `.py`:

- **M_* construction.** SymPy hardcodes the `Mstar` 3×8 matrix (py:211-217). The `.wl` instead
  assembles each row from additive exponent-weight vectors of the primitive coordinates
  (`chiCoreWeights`, `thermalWeights`, `nontrackPrefactorWeights`, `epsWWeights`, `epsEtaWeights`,
  wl:98-108), forms `compilerRows = (1+delU)*chiCoreWeights + (1+chi)*thermalWeights`, etc., and
  then *checks* `compilerRows - compilerTarget == 0` (wl:120). This is an exponent-ledger
  re-derivation of the map, not a re-typing of the matrix.
- **Transport factors.** SymPy asserts the closed forms `PhiT = rU*(rU/(rg*rc))^alpha`, etc.
  (py:133-139). The `.wl` *derives* them via a native `Solve[compilerRows.logVars == 0,
  depLogVars, Reals]` (the same-orbit condition, wl:144-153) and only then subtracts the SymPy
  closed form (wl:173-176). Solve produces the transport law independently rather than restating it.
- **Projectors.** SymPy inverts the 3×3 dependent sub-block (`Pdep.inv()`, `Edep*Pdep_inv`,
  py:244-251). The `.wl` uses `LinearSolve` on the augmented constrained system
  `Join[compilerRows, freeSelector]` solved column-by-column (wl:213-226) — a different
  linear-algebra route to the same quotient section.
- **House idioms only borrowed.** `expectZero`/`stripCE`/`$Assumptions` positivity and the
  `math -script` convention follow the established `.wl` house style, as the directive permitted.

*Substantiveness (no X−X, no vacuous checks).* Every `expectZero`/`expectTrue` compares two
independently-constructed objects: Solve/LinearSolve output vs. SymPy closed form, or
weight-vector assembly vs. target. Specifically on the items the directive flagged:

- *Transport / mismatch / projector split* (II–IV): all load-bearing. Section IV additionally adds
  idempotency and orthogonality checks (`Q^2-Q`, `O^2-O`, `Q.O`, `O.Q`, wl:242-245) that the SymPy
  script does NOT contain — extra genuine projector-property verification, all PASS as `{{0...}}`.
- *Cocycle laws* (VI): uses fully independent 21- and 32-ratio symbol sets with `r31 = r32*r21`
  (wl:259-261) and re-runs the Solve-derived `phiLogFor`/`mLogsFrom` on each, exercising the
  homomorphism over the real parameter space (not a special case).
- *Two-point orbit-lock equivalence* (VII): this is genuinely derived, NOT a vacuous `solve q==0`.
  It (1) solves both `compilerRows.logVars==0` and `quotientProjector.logVars==0` for the dependent
  logs and checks the two same-orbit laws **agree** (wl:289-299), cross-linking two independently
  built operators; and (2) constructs the explicit `mismatchToQ` chart matrix, computes
  `Det = 1+chi` (the nonzero invertibility fact that makes the lock an iff, wl:308-309), and only
  then solves `mismatchToQ.{tau,kappa,mu}==0` to show q=0 forces zero mismatch (wl:310-317). The
  invertible chart is the load-bearing content, so the equivalence is established, not assumed.

*Cross-engine agreement.* The `.wl` conclusions match the SymPy results in every section:
both engines produce `ln Phi_T = (-((1+delU)Log[rC*rGam]) + (2+chi+delU)Log[rU])/(1+chi)`,
`Phi_K = rC^2/rU`, the identical `m_T/m_K/m_mu`, identical `Q_quot`/`O_orb` supports, and the
identical cocycle/reduction outcomes (compare wl-output lines 24-26, 41-43, 60-61 against the
SymPy pretty-printed forms at py-output lines 33-95, 141-308). All comparison residuals are 0.

## Exec log assessment

**SymPy:** exit=0. The SymPy script was untouched; its exec log ends `# exit_code: 0` and the
committed `.txt` shows all Section I–VII residuals as `0` / zero matrices.

**Mathematica:** exit=0. `redteam/exec_logs/stage_199_mathematica.log` ends `# exit_code: 0`.
Notable lines from the committed output:
- `derived M_* rows - SymPy compiler rows = {{0,...},...}` → `PASS` (independent map re-derivation matches)
- `ln Phi_T from Solve = (-((1 + delU)*Log[rC*rGam]) + (2 + chi + delU)*Log[rU])/(1 + chi)` → `PASS` vs SymPy
- `det(mismatch-to-q chart) = 1 + chi` → `PASS` (invertible chart → orbit-lock is an iff)
- `q == 0 forces zero logarithmic mismatch = {0, 0, 0}` → `PASS`
- 42 `PASS:` lines, 0 `FAIL`.

**Output freshness:** confirmed. `.wl` mtime 2026-06-01 12:14:47; its committed `.txt` mtime
2026-06-01 12:17:54 (newer than the script). The SymPy `.txt` was also refreshed (mtime
2026-06-01 12:17:54) by the orchestrator re-run; the `.py` itself is unchanged at 2026-05-11.

## Material-change assessment

`material_change`: false. No derived result changed — the SymPy reference engine was not edited
and the new `.wl` only re-confirms the existing conclusions via a second engine. No downstream
unit's inputs are altered.

## Side observations (non-blocking)

- The SymPy script (untouched) still carries stale `STAGE 182` / `Stage 192` / `Stage 198` labels
  in banners and prints; the original auditor already recorded this as a cosmetic renumbering
  artifact touching no assertion. The new `.wl` correctly banners `STAGE 199`. Not a finding;
  not in scope to fix here.
- The `.wl` adds projector idempotency/orthogonality checks beyond the SymPy coverage — a strict
  superset of substantive checks, which strengthens the cross-engine verification.

## Verdict justification

The sole finding (missing dual-engine verification) is resolved: the new `.wl` is a genuinely
independent route — it re-derives `M_*` from additive exponent-weight vectors, obtains the
transport factors from a native `Solve` of the same-orbit condition, and builds the projectors via
`LinearSolve` on a constrained system, rather than re-typing the SymPy closed forms. Every check is
substantive (compares two independently constructed objects, with extra projector-property checks),
the orbit-lock equivalence is genuinely derived through the invertible `mismatchToQ` chart rather
than a vacuous `solve q==0`, and all conclusions agree with the SymPy engine at residual 0. Both
exec logs report exit 0 with all PASS, and the outputs are fresh. Verdict: verified.
