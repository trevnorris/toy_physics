---
unit_id: 008
batch: I.1
verifier_model: claude-opus-4-7[1m]
verify_date: 2026-05-25T17:25:00-06:00
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 1
findings_total: 1
material_change: false
---

# Verification -- unit 008 (v2 iteration)

## Per-finding outcomes

### F1 -- tautological_check

**Classification:** resolved

**What changed:**

Codex applied six textual edits across the two scripts, matching the directive exactly. (Diff captured at `/var/projects/toy_physics/research/pde_ledger/redteam/exec_logs/stage_008_diff.patch`.)

SymPy (`scripts/moving_throat_pde_stage008_projected_maxwell_extension_sympy_audit.py`):
- Lines 186-189: matched-Gaussian block now defines a fresh `H_Z_match = sp.exp(-w**2 / lam**2)` (textually retyped, not the Python `Z` binding) and integrates `W_match * H_Z_match` instead of `W_match * Z`.
- Line 204: new pre-ratio guard `assert_zero("Gaussian H=Z integrals match (independent computation)", I_WH_HZ_match - I_WZ_match)` inserted before the existing H=Z gauge ratio assertion on line 205.
- Lines 242-244: independent-profile block likewise retypes `H_Z_indep = sp.exp(-w**2 / lam**2)` and integrates `W_indep * H_Z_indep`.
- Line 245: new pre-ratio guard `assert_zero("independent-profile H=Z integrals match (independent computation)", I_WH_HZ_indep - I_WZ_indep)` inserted before the existing H=Z gauge ratio assertion on line 247.

Mathematica (`mathematica/moving_throat_pde_stage008_projected_maxwell_extension_mathematica_audit.wl`):
- M2 Pair A: line 35 now reads `Integrate[gaussianWeight*Exp[-w^2/lambda^2], ...]` instead of `gaussianWeight*localizedGaussian`. Lines 38-41 add the equality guard `If[FullSimplify[gaussGaugeWeight - gaussOverlap] =!= 0, ... Exit[1]]; Print["PASS: M2 Pair A H=Z integrals independently match"]` before the existing gauge residual on line 42.
- M2 Pair B: same treatment at line 60 and lines 63-66.
- M4 matched block: lines 108-114 add a fresh `matchedGaugeOverlap = Integrate[matchedWeight*Exp[-w^2/lambda^2], ...]` plus equality guard before the residual list. Lines 119 and 134 replace `xi*matchedOverlap/matchedOverlap - xi` with `xi*matchedOverlap/matchedGaugeOverlap - xi` in both the residual list and the dedicated `If` check.
- M7 Lorentzian numeric block: lines 167-169 introduce a fresh `lorentzGaugeOverlap = Integrate[lorentzWeight*Exp[-w^2/lambda^2], ...]` and line 174 swaps `lorentzGaugeWeight` for `lorentzGaugeOverlap` in the numeric residual definition.

**Assessment:**

The edit matches the directive's "required change" exactly. The four Mathematica sites named in the directive (M2 Pair A, M2 Pair B, M4 matched, M7 numeric) and the two SymPy sites (matched-Gaussian, independent-Gaussian) are all converted: the H=Z integrand on each profile is now `Exp[-w^2/lambda^2]` / `sp.exp(-w**2 / lam**2)` written as an independent textual expression rather than a re-use of the existing `localizedGaussian` / `Z` binding. The new equality assertions `I_WH_HZ - I_WZ == 0` precede the ratio assertions in both engines so the integrator's closed form for the H integrand must agree with the closed form for the Z integrand on its own merits.

The new assertions are non-tautological: each one feeds the integrator two textually distinct expressions and demands their simplified difference be zero. If `Integrate` produced the wrong value for either side, the equality guard would now fail. The downstream ratio assertion `xi*I_WZ/I_WH_HZ - xi == 0` is preserved verbatim, so the H=Z = xi claim is still being asserted, only now backed by a substantive integral-equality check.

The H=1 side of the mutation guard (always genuine in v1) is preserved unchanged:
- SymPy line 199: `xi_eff_H1_match = sp.simplify(xi * I_WZ_match / I_WH_H1_match)` with `I_WH_H1_match = sp.simplify(sp.integrate(W_match, ...))` (line 185) — integrates W alone (H=1), structurally distinct from the W*Z integrals.
- SymPy line 206 still asserts `xi_eff_H1_match - xi/sqrt(2) == 0`.
- Mathematica M4 line 120 still asserts `xi*matchedOverlap/matchedNorm - xi/Sqrt[2]` (where `matchedNorm = Integrate[matchedWeight, ...]` is the H=1 case, line 100-101).
- Mathematica M4 dedicated If on line 137 unchanged.

No collateral edits: the diff is confined to the named line ranges. No paper/notes touched.

## Exec log assessment

**SymPy:** exit=0 (inferred from `STATUS: PASS` in the post-fix output txt; per-stage log `stage_008_sympy.log` is missing from `exec_logs/`). The output file `scripts/output/moving_throat_pde_stage008_projected_maxwell_extension_sympy_audit.txt` (mtime 2026-05-25 17:13, newer than the script's 17:11) ends with:
```
  4. With a delta-localized source instead, the coupling mismatch remains, even if H = Z.
STATUS: PASS
```
The `assert_zero` helper (lines 30-33 of the script) is silent on success and only raises on residue != 0, so the absence of explicit "H=Z integrals match" lines in the SymPy stdout is by design; the script reaching `STATUS: PASS` requires every `assert_zero`, including the two new equality guards on lines 204 and 245, to have evaluated to zero. If `I_WH_HZ_match - I_WZ_match` or `I_WH_HZ_indep - I_WZ_indep` had failed to simplify to zero, the script would have raised AssertionError before printing `STATUS: PASS`.

**Mathematica:** exit=0 (inferred from `STATUS: PASS` in the post-fix output txt; `stage_008_mathematica.log` in `exec_logs/` is stale from 2026-05-21 11:16 and predates Codex's edits, so it does not reflect the post-fix run). The fresh output file `mathematica/output/moving_throat_pde_stage008_projected_maxwell_extension_mathematica_audit.txt` (mtime 2026-05-25 17:13) shows the new PASS lines:
```
PASS: M2 Pair A H=Z integrals independently match
M2 Pair A normalization residual = 0
M2 Pair A H=Z residual = 0
PASS: M2 Pair A Gaussian-Gaussian
PASS: M2 Pair B H=Z integrals independently match
...
PASS: M4 H=Z matched integrals independently match
M4 residuals = {0, 0, 0, 0, 0, 0, 0}
...
M7 Lorentzian numeric H=Z residuals = {0, 0}
...
STATUS: PASS
```
M7's numeric `{0, 0}` residual now reflects the substantive `lorentzOverlap / lorentzGaugeOverlap` ratio at the two parameter pairs `{lambda->1, sigma->1/2}` and `{lambda->1, sigma->2}`, not the algebraic `X / X` cancellation that v1 had.

**Output freshness:** Confirmed. Both output txt files (mtimes 17:13) are newer than their corresponding script files (mtimes 17:11). The `stage_008_mathematica.log` and missing `stage_008_sympy.log` in `exec_logs/` are stale orchestrator artifacts; the canonical post-fix evidence is the output txt files, which were regenerated after Codex's edits at the same exec timestamp.

## Material-change assessment

`material_change`: false.

The verified mathematical relationship is exactly what the original v1 audit certified: `xi_eff_proj(H=Z) = xi` for matched-Gaussian, independent-Gaussian, and Lorentzian observer kernels, plus the Gaussian explicit constants `Z_int = sqrt(pi)*lambda`, `Z2_int = sqrt(pi/2)*lambda`, `I_WZ_match = sqrt(2)/2`, the `mu0/Z_int` matched-source coupling, the `sqrt(2)` delta-vs-matched ratio, the H=1 matched gauge `xi/sqrt(2)`, the H=1 independent gauge `lambda*xi/sqrt(lambda^2+sigma^2)`, and the regulator-limit-zero result. Codex changed only how the H=Z verification machinery is structured (integrand textually re-typed, pre-ratio equality guard added). No new constants, no new symbolic claims, no altered numeric values. Downstream stages that consume stage 008's exported matching conditions (eq:stage005-hz `H=Z, xi_eff_proj = xi` and eq:stage005-source-match `S=Z/Z_int, mu_eff_proj = mu0/Z_int`) see no change to those exports.

## Side observations (non-blocking)

- The orchestrator's `redteam/exec_logs/stage_008_sympy.log` was not regenerated for this iteration -- only the diff and the stale Mathematica log are present. The canonical post-fix evidence is the script's own output txt under `scripts/output/`. This is an orchestrator bookkeeping observation, not a verification blocker.
- The stale `stage_008_mathematica.log` (mtime May 21) still in `exec_logs/` predates the fix; it should be regenerated by the orchestrator at the next exec-mathematica call. Non-blocking for this verification because the actual post-fix run output is captured in the script's own output txt.
- M6 (line 160 `xi*zIntSym/zIntSym - xi`) was flagged in the v2 audit's assertion inventory as `partial (substitution identity)`. The directive deliberately excluded it because it is a symbolic-reduction algebra check (not a concrete-profile check), so its X/X structure is intentional. No action required; consistent with the directive's scope.

## Verdict justification

Every required edit in directive F1 is in place, the diff is confined to the named ranges, both scripts re-ran post-fix to `STATUS: PASS`, the new pre-ratio equality guards are non-tautological (they feed the integrator two textually distinct integrands and demand their closed forms agree), and the H=1 mutation-guard side that was already genuine is preserved. The mathematical content of stage 008 is unchanged; only the verification machinery is stronger. `material_change` is false.
