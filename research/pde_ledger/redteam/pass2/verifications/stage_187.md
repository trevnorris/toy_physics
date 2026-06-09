---
unit_id: 187
batch: V.2
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-06-09T00:00:00Z
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 2
findings_total: 2
material_change: false
---

# Verification — unit 187

## Per-finding outcomes

### F1 — insufficient_verification (SymPy engine parity)

**Classification:** resolved

**What changed:**
`scripts/...stage187..._sympy_audit.py:107` — one line added immediately after the
determinant `print` at py:106:
`expect_zero("selected minor determinant", minor.det() - (1 + chi))`.
The diff (`stage_187_diff.patch`, py hunk) shows this is the ONLY change to the `.py`:
a single `+` line, no collateral edits.

**Assessment:**
Correct and non-tautological. `minor` is the hand-built literal
`[[0,0,1+chi],[-1,1,-F],[-1,0,0]]` (py:101-105); its determinant expands to `1+chi`,
so `minor.det() - (1+chi)` is identically 0 only if the matrix is transcribed
correctly — a mis-transcription of any entry would make the residual nonzero and fail.
This brings SymPy to parity with the Mathematica assert (wl:95). The new line
`selected minor determinant = 0` appears in the committed SymPy `.txt` (line 16) and
the script exits 0. Matches the directive's required change and acceptance exactly.

### F2 — mathematica_transliteration (USER-AUTHORIZED full re-author)

**Classification:** resolved

**What changed:**
`mathematica/...stage187..._mathematica_audit.wl` fully re-authored (diff replaces the
whole body). The old hand-coded `rowTr/rowNt/rowEta` literals (old wl:34-36), the
hand-typed `m` matrix (old wl:44-48), the hand-typed `minor` matrix (old wl:56-60), and
the simultaneous `Solve[{rowTr==0,...}, {deta,dt,dm}]` (old wl:81) are all DELETED.

The new `.wl` derives everything from physical monomials:
- monomial invariants as functions of a state association: `cTrInvariant`/`cNtInvariant`/
  `epsEtaInvariant` (wl:51-64), built from the exponent structure (1+deltaStar), (1+chiStar),
  eStar, -fStar over `gamma c / KU`, `Pi^2 T / (ell^2 KU)`, etc.;
- `baseState` and `targetState = base * Exp[delta]` (wl:46-49);
- `finiteLogRatio` (wl:66-70) = `PowerExpand[Log[invariant[targetState]/invariant[baseState]]]`,
  producing `derivedRows` (wl:73) — the rows EMERGE from the monomial exponents, not hand-typed;
- `mStar` built by `Coefficient` extraction from the derived rows (wl:75-76), then
  `matrixRows = mStar.deltaVector` cross-checked against `derivedRows` (D1, wl:85-91);
- `minor = mStar[[All,{5,7,8}]]` (wl:93) sliced from the DERIVED matrix, not hand-typed;
- fibre solved by TRIANGULAR ELIMINATION via `solveLinearFor` (wl:97-104): eta from row3,
  dt from row1, dm from row2 with back-substitution — a genuinely distinct route from the
  `.py`'s `sp.solve(...)` simultaneous solve.

**Assessment:**
This is a genuine independent derivation, not a relabel or a port. A wrong exponent in any
invariant would propagate into `derivedRows`, break the D1 row==matrixRow check, the minor
determinant, AND the solved laws (which are compared against independently hand-written
`detaExpected/dtExpected/dmExpected`, wl:113-120). Output banner is `STAGE 187`, all 12
`PASS:` lines present in the committed `.txt`, `math -script` exits 0, `.py` untouched by F2.
Acceptance (a)-(d) all met. The 173-lesson independence check passes: the second engine no
longer solves the same hand-coded system — it reconstructs the rows from the physical
monomial ratios via a state-Exp[delta] map + Coefficient extraction, and solves by a
different algorithm.

## Exec log assessment

**SymPy:** exit=0. Notable lines:
`det selected minor (Delta_eta, Delta_mu, Delta_T) = chi0_star + 1` then the new
`selected minor determinant = 0`; all `... finite law = 0` and `row_* after solve = 0`.

**Mathematica:** exit=0. Notable lines:
`Exact finite log-ratio equations derived from physical monomials:`, `D1 C_tr/C_nt/
epsilon_eta log ratio - M_* Delta x = 0` each with `PASS:`, `selected minor determinant = 0`
+ `PASS:`, three `Delta_* finite law` PASS, three `D4 ... row after finite fibre solve` PASS.
Solved laws (Delta_eta = 2*dc - du, etc.) algebraically identical to the SymPy forms and to
the carry-forward block printed by both engines.

**Output freshness:** confirmed. `.py` mtime 1780985564 < `.txt` 1780986383; `.wl` mtime
1780985600 < `.txt` 1780986383. Both committed `.txt` outputs are newer than their scripts and
their content matches the exec logs (including the new assertion lines).

## Material-change assessment

`material_change`: false. F1 adds a guard assertion (no result changed). F2 is a method-only
re-author: every emitted closed form is preserved — the three rows, the matrix `M_*`, the
selected-minor determinant `= 1+chi0_star`, and all three orbit laws (Delta_eta, Delta_T,
Delta_mu) are byte-equivalent to the pre-fix deliverables and to the report's reconciliation
table (9 deliverable values, all MATCH). No downstream unit depends on a changed value.

## Side observations (non-blocking)

The `.py` retains its own independent row-derivation route (`Ctr_ratio`/`Cnt_ratio`/`eps_ratio`
+ `expand_log`, py:70-86), which is distinct in mechanism from the `.wl`'s new state-association
+ `Exp[delta]` + `Coefficient` route — so the two engines now derive the rows by genuinely
different constructions, strengthening the cross-check beyond the directive's minimum.

## Verdict justification

Both findings are resolved. F1 is a single non-tautological assert added to the `.py` (the only
`.py` change), present in the committed output, exit 0. F2's `.wl` re-author genuinely derives
the finite rows and matrix from physical monomial ratios (state association + Exp[delta] map +
Coefficient extraction) and solves by triangular elimination — the old hand-coded rows/matrix/
minor and the shared `Solve` are gone, so a wrong exponent now fails. All four deliverables
(D1-D4) still verify in both engines with banner STAGE 187 and exit 0; all 9 deliverable values
are preserved (method-only). Verdict: verified.
