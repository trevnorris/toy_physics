---
unit_id: 035
batch: II.1
verifier_model: claude-opus-4-7[1m]
verify_date: 2026-05-22T00:00:00Z
verdict: verified
sympy_exit: n/a
mathematica_exit: n/a
findings_resolved: 2
findings_total: 2
material_change: false
---

# Verification — unit 035

## Per-finding outcomes

### F1 — tautological_check

**Classification:** resolved

**What changed:**
In `mathematica/moving_throat_pde_stage035_dimensionless_normalization_locus_mathematica_audit.wl`, the tautological assertion `expectZero["support split residual", gBReqSqOverVarpi2 - (alphaReqTarget - alphaMix)];` (formerly line 111) has been deleted. The diff (`stage_035_diff.patch` line 45) confirms the line was removed. The current file at lines 108-110 retains the surrounding context:

```
alphaMix = FullSimplify[Chi^2/(OmegaU^2*Delta0), Assumptions -> $Assumptions];
gBReqSqOverVarpi2 = FullSimplify[alphaReqTarget - alphaMix, Assumptions -> $Assumptions];
Print["g_B,req^2 / varpi^2 = ", fmt[gBReqSqOverVarpi2]];
```

**Assessment:**
Exactly the change requested. The Mathematica output (`moving_throat_pde_stage035_dimensionless_normalization_locus_mathematica_audit.txt` line 28) still emits `g_B,req^2 / varpi^2 = ...` and no longer contains a `support split residual = 0` / `PASS: support split residual` line pair (verified against the saved output). The value is still computed and printed; only the false-confidence assertion is gone. No collateral edits.

### F2 — mathematica_transliteration

**Classification:** resolved

**What changed:**
Seven literal-driven LHS computations in the `.wl` script were switched to engine-derived quantities, matching the directive's enumerated edits exactly. From the diff and current file state:

- Line 64: `dF = FullSimplify[D[fTarget, xi], …]` → `dF = FullSimplify[D[f, xi], …]` (F2.1)
- Line 76: `Limit[(1 - xi)*fTarget, xi -> 1, …]` → `Limit[(1 - xi)*f, xi -> 1, …]` (F2.2)
- Line 87: `Series[(epsSoft*fTarget /. xi -> 1 - epsSoft), …]` → `Series[(epsSoft*f /. xi -> 1 - epsSoft), …]` (F2.3)
- Line 99: `Limit[alphaReqTarget, xi -> 1, …]` → `Limit[alphaReq, xi -> 1, …]` (F2.5)
- Line 112: `Series[fTarget, {xi, 0, 2}]` → `Series[f, {xi, 0, 2}]` (F2.6)
- Line 117: `Series[alphaReqTarget, {xi, 0, 2}]` → `Series[alphaReq, {xi, 0, 2}]` (F2.7)
- F2.4 (`alphaReq = FullSimplify[alphaX, ...]`) explicitly required "no edit"; the file at line 92 is unchanged.

The literal closed-form targets (`fTarget`, `dFTarget`, `softConstTarget`, `alphaReqTarget`, `alphaCritTarget`, `fSeriesTarget`, `alphaSeriesTarget`) remain unchanged at their original lines — as the directive required, since those literals encode the claim under test.

**Assessment:**
The structural independence property is now established. The dependency chain from the physical premises (`kappa0Sq`, `kappa1Sq`, `alphaX`, `nX`) flows to `f` (line 50: `f = FullSimplify[nX/(beta0*kappa0Sq/A), …]`) and `alphaReq` (line 92: `alphaReq = FullSimplify[alphaX, …]`), and every LHS of an `expectZero` against a literal target is now derived from one of these physically-anchored quantities, not from a re-simplification of the same literal.

Non-tautology check, per literal target:
- `fTarget`: line 61 check is `f - fTarget` where `f` comes from `nX`. A wrong coefficient in `fTarget` (e.g., changing `11*xi` to `12*xi`) would surface here as nonzero. Non-tautological.
- `dFTarget`: line 71 check is `dF - dFTarget` where `dF = D[f, xi]`. A wrong literal coefficient (e.g., `121*xi^3 → 122*xi^3`) would not be matched by Mathematica's independent `D[f, xi]`. Non-tautological.
- `softConstTarget`: line 83 check is `softConst - softConstTarget` where `softConst = Limit[(1-xi)*f, xi -> 1, ...]`. Non-tautological.
- `softConstTarget` (via series): line 90 check is `softScaledSeries - softConstTarget` where `softScaledSeries` series-expands `epsSoft*f`. Non-tautological.
- `alphaReqTarget`: line 105 check is `alphaReq - alphaReqTarget` where `alphaReq = FullSimplify[alphaX, ...]` (alphaX is the physical loading definition on line 40). Non-tautological.
- `alphaCritTarget`: line 106 check is `alphaCrit - alphaCritTarget` where `alphaCrit = Limit[alphaReq, xi -> 1, ...]`. Non-tautological.
- `fSeriesTarget`: line 122 check is `fSeries - fSeriesTarget` where `fSeries = Series[f, {xi, 0, 2}]`. Non-tautological.
- `alphaSeriesTarget`: line 123 check is `alphaSeries - alphaSeriesTarget` where `alphaSeries = Series[alphaReq, {xi, 0, 2}]`. Non-tautological.

The Mathematica output transcript confirms all eight residuals still evaluate to `0` and `PASS`. The `dF/dxi` printed form (output line 12: `((9*delta + 11*xi)^3*(81*delta^3 + 297*delta*xi^2 + 121*xi^3 + 9*delta^2*(8 + 21*xi)))/(81*(-1 + xi)^2*(9*delta^2 + 18*delta*xi + 11*xi^2)^3)`) is the engine-canonicalized form of `D[f, xi]` rather than a re-simplification of the literal `dFTarget`, consistent with the directive. No collateral edits.

## Exec log assessment

**SymPy:** exit=n/a (no exec log captured by the orchestrator at `/var/projects/toy_physics/research/pde_ledger/redteam/exec_logs/stage_035_sympy.log`; the SymPy script was not modified and its saved output `.txt` is unchanged).

**Mathematica:** exit=n/a (no exec log captured at `/var/projects/toy_physics/research/pde_ledger/redteam/exec_logs/stage_035_mathematica.log`).

Notable lines from the saved Mathematica output `.txt` (which serves as the post-fix run record):
- Line 13-14: `dF/dxi - manifestly positive form = 0 / PASS: dF/dxi - manifestly positive form` (the LHS now from `D[f, xi]`).
- Line 18-19: `softening constant - closed form = 0 / PASS` (LHS from `Limit[(1-xi)*f, …]`).
- Lines 28-34: the `support split residual` PASS line pair is absent; `g_B,req^2 / varpi^2 = …` print line is preserved.
- Line 36: `Stage 035 Mathematica audit passed.` confirms exit 0.

**Output freshness:** confirmed. Mathematica script mtime = 17:34, Mathematica output mtime = 17:36 (newer). SymPy script mtime = Apr 3 12:05 (untouched), SymPy output mtime = 17:35 (regenerated this session). Both saved outputs were re-generated post-fix.

## Material-change assessment

`material_change`: false.

The edits are restricted to the Mathematica audit script's internal verification structure. No derived numerical or symbolic result that downstream units could consume has changed:
- The closed-form literals (`fTarget`, `dFTarget`, `softConstTarget`, `alphaReqTarget`, `alphaCritTarget`, `fSeriesTarget`, `alphaSeriesTarget`) are byte-identical.
- The printed values (`F(xi,delta)`, `R_target`, `dF/dxi`, `softening constant`, `alpha_req`, `alpha_crit`, `g_B,req^2 / varpi^2`, the series expansions) are all algebraically unchanged — the diff shows only LHS-source substitutions in `expectZero` arguments and the deletion of one assertion, none of which alter any quantity exported to other stages.
- The SymPy script and SymPy output are unchanged.

Downstream units > 035 do not need re-auditing.

## Side observations (non-blocking)

- The script's banner on line 26 reads `"STAGE 018 — DIMENSIONLESS NORMALIZATION LOCUS"` and the SymPy output groups the audit under `STAGE 18.1/18.2/18.3/18.4`. The unit is 035, not 018, so the banner appears stale relative to the unit numbering — but this was pre-existing (not introduced by Codex in this iteration) and is purely cosmetic. Not part of any finding; flagging only.
- The orchestrator did not capture `stage_035_sympy.log` or `stage_035_mathematica.log` exec logs at the expected paths. The saved `.txt` outputs are fresh and serve as the record, but a future verifier may want a real exec log for transcript fidelity.

## Verdict justification

Both findings are fully resolved. F1's tautological `expectZero` is gone from the script and from the saved output, while the surrounding print is retained. F2's seven enumerated LHS substitutions are all present in the diff and the current file, with no collateral edits; the structural independence the directive demanded — LHS quantities derived from the physical premises (`nX`, `alphaX`) rather than from re-simplified literal targets — now holds for every closed-form comparison. The Mathematica output transcript shows all eight `expectZero` checks still PASS and the script exits 0. The SymPy script was correctly left untouched. No material change to any downstream-visible result.
