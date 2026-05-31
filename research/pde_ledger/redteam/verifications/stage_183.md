---
unit_id: 183
batch: V.2
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-05-30T01:50:00-06:00
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 1
findings_total: 1
material_change: false
---

# Verification — unit 183

## Per-finding outcomes

### F1 — tautological_check

**Classification:** resolved

**What changed:**
Codex replaced the tautological forward-substitution rigidity block with a genuine
nonzero-prefactor check in both engines, and added matching nonzero helpers.

- SymPy `scripts/moving_throat_pde_stage183_triangular_normal_form_sympy_audit.py`:
  - Added `expect_nonzero` helper (lines 35–39): `sp.simplify(expr)`, prints it, raises
    `AssertionError` iff `expr == 0`.
  - Rigidity block (lines 113–122): the three old `expect_zero(... .subs(→0))` lines
    are gone; replaced by `expect_nonzero` on `C_tr`, `A_tr`, and
    `dressing_pref = sp.simplify(eps_eta/(1-eps_eta))`.
- Mathematica `mathematica/moving_throat_pde_stage183_triangular_normal_form_mathematica_audit.wl`:
  - Added `expectNonzero` helper (lines 26–30): `FullSimplify` under `$Assumptions`,
    prints, `fail` iff `res === 0` else `pass` (correctly inverted vs `expectZero`).
  - Rigidity block (lines 85–92): old three `expectZero` zero-substitutions replaced by
    `expectNonzero` on `cTr`, `aTr`, `dressingPref`.

The diff (`stage_183_diff.patch`) contains exactly two hunks per file — the helper
insertion and the rigidity-block replacement. No collateral edits; the imported
`Theta1/Xi1/R1`, `SigmaNT_def`, `A_tr`/`C_tr`, and all A1–A6 checks are untouched.
`deviation: none` is accurate.

**Assessment:**
The edit matches the directive verbatim (helper bodies and replacement blocks identical
to the prescribed text) and correctly addresses the finding. The new checks are
**non-tautological**:

- `C_tr = chi0*deltaU/((chi0+1)(deltaU+1)(chi0+deltaU+1))` and
  `A_tr = 2*chi0/((chi0+1)(deltaU+1))` are rational functions of the strictly-positive
  branch variables (`chi0, deltaU` declared `positive=True`). `sp.simplify` does NOT
  collapse them to 0, so `expr == 0` is `False` and the check passes — but it genuinely
  *can* fail: had a numerator been mistyped/degenerated to zero, `simplify` would return
  `0` and the `AssertionError` would fire. This is exactly the non-trivial `⟹` content
  (nonzero diagonal prefactors ⟹ map invertible ⟹ vanishing observables force vanishing
  slippages) that the auditor said was missing. It is the structural complement to the
  already-present inverse round-trips A3/A5/A6 that prove invertibility constructively.
- `dressing_pref = eps_eta/(1-eps_eta)`: the tested numerator `eps_eta` is `positive`,
  so nonzero on the branch; sympy renders `-eps_eta/(eps_eta-1)`, Mathematica
  `epsEta/(1-epsEta)` — both nonzero. The `(1-eps_eta)` pole at `eps_eta=1` is outside
  the branch and does not affect the numerator-nonzero test.
- The earlier "bare `SigmaNT` detached from `SigmaNT_def`" structural complaint is moot:
  the new block no longer references `SigmaNT`/`sigmaNT` at all (it tests prefactors), so
  the detachment is removed rather than papered over.

Mathematica `expectNonzero` pass/fail is correctly inverted relative to `expectZero`
(`fail` iff `=== 0`), confirmed against the log's three `PASS:` lines. Both engines print
the three prefactor *values* (not three `= 0` lines), satisfying the directive's
verification criterion (a) and (b).

## Exec log assessment

**SymPy:** exit=0. Notable lines (rigidity section now prints values, not zeros):
```
C_tr (Theta_1 <- Sigma_tr prefactor) nonzero on branch = chi0*deltaU/((chi0 + 1)*(deltaU + 1)*(chi0 + deltaU + 1))
A_tr (Xi_1 <- Sigma_tr feed-through) nonzero on branch = 2*chi0/((chi0 + 1)*(deltaU + 1))
eps_eta/(1-eps_eta) (R_1+Xi_1 <- Sigma_eta prefactor) nonzero on branch = -eps_eta/(eps_eta - 1)
```
Pre-existing load-bearing checks still pass: `Xi_1 - (A_tr Sigma_tr + Sigma_nt) = 0`,
`R_1 + Xi_1 + eps_eta/(1-eps_eta) Sigma_eta = 0`, all four inverse residuals `= 0`,
`A_tr/C_tr = 2*(chi0 + deltaU + 1)/deltaU`.

**Mathematica:** exit=0. Notable lines:
```
PASS: C_tr (Theta_1 <- Sigma_tr prefactor) nonzero on branch
PASS: A_tr (Xi_1 <- Sigma_tr feed-through) nonzero on branch
PASS: eps_eta/(1-eps_eta) (R_1+Xi_1 <- Sigma_eta prefactor) nonzero on branch
Stage 183 Mathematica audit passed.
```
All A1–A6 mirrors still `PASS`. Engines agree on the printed prefactor forms (modulo CAS
formatting: sympy `-eps_eta/(eps_eta-1)` ≡ wl `epsEta/(1-epsEta)`).

**Output freshness:** confirmed. Script mtimes 01:21:03 (py) / 01:21:16 (wl); output
mtimes 01:40:47 (py txt) / 01:40:57 (wl txt). Outputs are newer than scripts → fresh,
regenerated post-fix.

## Material-change assessment

`material_change`: false. The edit only strengthens a verification block (replacing a
tautology with a genuine nonzero test). No derived result, constant, or coordinate
definition changed — `Theta1`, `Xi1`, `R1`, `Sigma_nt`, `A_tr`, `C_tr`, and all inverse
formulas are byte-identical. The "Carry-forward formulas" print block is unchanged.
Downstream units depending on Stage 183 outputs see no change in any propagated quantity.

## Side observations (non-blocking)

- The banner still reads "STAGE 166 — TRIANGULAR NORMAL FORM …" in both engines (sympy:42,
  wl:32) while the unit is 183. This pre-existed the fix and is outside F1's scope; not a
  verification blocker, just a stale label.
- The new SymPy `expect_nonzero` helper omits the `expand` that `expect_zero` uses
  (`sp.simplify(expr)` vs `sp.simplify(sp.expand(expr))`). Irrelevant for the three
  rational-function prefactors tested here; noted only for completeness.

## Verdict justification

The single finding (F1) is fully resolved. Codex's edits match the directive verbatim in
both engines, the diff is scoped to exactly the helper additions and the rigidity-block
replacement with no collateral changes, both scripts exit 0, and outputs were regenerated
after the fix. Critically, the replacement is non-tautological: the three diagonal
prefactors are rational functions of the branch variables that the `expect_nonzero` /
`expectNonzero` assertions would genuinely flag if any degenerated to zero, supplying the
non-trivial `⟹` rigidity content that the original `⟸`-only substitutions lacked. No
regressions, no new findings that block. Verdict: `verified`.
