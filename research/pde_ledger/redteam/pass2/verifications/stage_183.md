---
unit_id: 183
batch: V.2
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-06-09T00:30:00Z
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 2
findings_total: 2
material_change: false
---

# Verification — unit 183

## Per-finding outcomes

### F1 — mathematica_transliteration (USER-AUTHORIZED full re-author)

**Classification:** resolved

**What changed:**
`mathematica/...stage183...wl` fully re-authored (diff @ exec_logs/stage_183_diff.patch; current file read in full). The old `.wl` hand-coded `theta1`/`xi1`/`r1` with `sigmaTr` as an opaque carried symbol (old wl:39-50) and ran nine `expectZero`/`expectNonzero` asserts in 1:1 order/labels with the `.py`. The re-author replaces this with a genuinely different mechanism:
- raw observables `thetaRaw`/`xiRaw`/`rhoRaw` are built from the MICROSCOPIC slippages `sigmaChi,sigmaDel,sigmaZ,sigmaEps,sigmaEta` (wl:100-117); `sigmaTr` no longer appears as an opaque symbol;
- the prefactors are EXTRACTED, not hand-defined: `cTr=cFromThetaChi`, `aTr=aFromXiChi`, `etaPref` via `Coefficient[...]` (wl:121-139), versus the `.py` which hand-writes `A_tr`/`C_tr` closed forms (py:85-86);
- `rowOf` builds coefficient-vector rows over the 5 raw slippages and the triangular map is checked by matrix factorization `rawObservableRows - triangularMatrix.branchCompiler` (`expectMatrixZero`, wl:170-173);
- the inverse is obtained by symbolic `Solve` of `{thetaObs,xiObs,rhoObs} == triangularMatrix.{coordTr,coordNt,coordEta}` (wl:187-193);
- triple-rigidity uses `Reduce[branchAssumptions && expr==0,...]` branch zero-set tests (`expectNoBranchZero`, wl:45-52, 239-241) plus a `Reduce`-based iff/biconditional (wl:243-262).

**Assessment:**
Independence is genuine and visibly present — the route (microscopic raw rows → `Coefficient` extraction → matrix factorization → `Solve`-based inverse → `Reduce` zero-set) is not the `.py`'s hand-coded forward-substitution sequence, and it can FAIL on a wrong coefficient (the extracted `cTr`/`aTr` flow into both the factorization check and the inverse asserts). Deliverables preserved and all PASS in the exec log:
- D1 Σ_nt: `Canonical branch-adapted Sigma_nt = 0` PASS (wl:152) — `xiRaw - aTr·trackingRaw` vs the canonical closed form.
- D2: `raw observables - triangular matrix.compiler` (matrix) PASS; `Xi raw - (A_tr Sigma_tr + Sigma_nt) row` PASS; `Rho raw + eps_eta/(1-eps_eta) Sigma_eta row` PASS.
- D3 (three inverses): `inverse Sigma_tr/Sigma_nt/Sigma_eta formula = 0` all PASS, AND `inverse matrix.normal-form matrix - identity` PASS — these compare the `Solve`d inverse against independently-stated closed forms, so a wrong `A_tr`/`C_tr`/`Σ_nt` coefficient fails. NOT a round-trip by construction (the `.py`'s `SigmaTr_inv` is `SigmaTr` by construction; the `.wl`'s is solved from the abstract matrix).
- D4: `A_tr/C_tr closed form = 0` PASS; printed `(2*(1 + chi0 + deltaU))/deltaU`.
- D5: `C_tr`/`A_tr`/`dressing` diagonal prefactor `zero-set = False` → all PASS on branch `chi0,deltaU>0, 0<epsEta<1`; bonus `zero observables imply zero normal-form defect` PASS (directly exercises the iff the original report flagged as only indirectly tested). No collateral edits beyond the re-author; `.py` untouched.

### F2 — stale_output

**Classification:** resolved

**What changed:**
No source edit (per directive). Both committed transcripts regenerated post-fix.

**Assessment:**
`scripts/output/...sympy_audit.txt:3` and `mathematica/output/...mathematica_audit.txt:3` both read `STAGE 183 — TRIANGULAR NORMAL FORM OF THE COHERENT DEFECT`; zero `FAIL` lines in either committed `.txt`; mtimes 2026-06-09 00:24 are newer than the `.wl` (06-08 22:36) and `.py` (06-03). Banner staleness closed.

## Exec log assessment

**SymPy:** exit=0. Banner `STAGE 183` (log:8); all `= 0` residuals (`Xi_1 - (A_tr Sigma_tr + Sigma_nt) = 0`, `A_tr/C_tr - 2(1+chi0+deltaU)/deltaU = 0`, `Sigma_tr/Sigma_nt/Sigma_eta inverse = 0`); nonzero prefactors printed. `.py` is unchanged (not in diff).

**Mathematica:** exit=0. Banner `STAGE 183` (log:8); every check PASS including `raw observables - triangular matrix.compiler` (log:29), the three `inverse ... formula` (log:43-49), `inverse matrix.normal-form matrix - identity` (log:51), the three branch zero-set `PASS` (log:60-63), and `zero observables imply zero normal-form defect` (log:66). Closes with `Stage 183 Mathematica audit passed.`

**Output freshness:** confirmed — both `.txt` mtimes (2026-06-09 00:24:03) are newer than the corresponding scripts; both banners read STAGE 183; no FAIL lines.

## Material-change assessment

`material_change`: false. Only the `.wl` was re-authored (independent route) and the two transcripts refreshed; no derived value, coefficient, or closed form changed. Every printed form matches the prior (correct) results and the `.py`. No downstream unit is affected.

## Side observations (non-blocking)

- The `.wl` independently RE-DERIVES the prefactors and the canonical `Σ_nt`, so it now meaningfully exceeds the `.py`'s coverage (it would catch a transcription error in the Stage-182 inputs that the old mirror could not). This addresses the original F2 rationale directly.
- The `inverseMatrix` is assembled via `Coefficient` over `Tuples` then `Partition[...,3]`; correctness is independently corroborated by the `inverse matrix . normal-form matrix - identity` PASS, so the assembly is not load-bearing on its own.

## Verdict justification

Both findings are `resolved`. The Mathematica `.wl` no longer mirrors the `.py`'s hand-coded `theta1`/`xi1`/`r1` + nine-assert sequence; an independent derivation (raw-slope row extraction → matrix factorization → symbolic `Solve` inverse → `Reduce` branch zero-sets) is visibly present and can fail on a wrong prefactor. All five deliverables (D1–D5) still verify, the inverse reconstructions genuinely test the prefactors rather than round-tripping by construction, both engines exit 0 with every in-file check PASS, the refreshed transcripts carry the STAGE 183 banner, and the diff confirms the `.py` was not touched. Verdict: verified.
