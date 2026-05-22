---
unit_id: 040
batch: III.1
verifier_model: claude-opus-4-7[1m]
verify_date: 2026-05-22T00:00:00-06:00
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 4
findings_total: 4
material_change: false
---

# Verification — unit 040

## Per-finding outcomes

### F1 — insufficient_verification

**Classification:** resolved

**What changed:**
- SymPy `scripts/moving_throat_pde_stage040_generalized_selected_branch_sympy_audit.py:68-79` inserts an explicit perturbed-matrix construction (`M_base`, `z_vec`, `M_perturbed = M_base - alpha_req_q * (z_vec * z_vec.T)`), builds the eigenvector `e_minus = (1, q*xi/(delta+xi))^T`, computes `eig_residual = simplify(M_perturbed * e_minus - lam_minus * e_minus)`, and asserts both components are zero via `expect_zero`.
- Mathematica `mathematica/moving_throat_pde_stage040_generalized_selected_branch_mathematica_audit.wl:67-77` adds the mirror block using `Outer[Times, zVec, zVec]` and `mPerturbed.eMinus - lamMinus eMinus`, with two `expectZero` calls.

**Assessment:**
The edits match the directive's required-change block essentially verbatim. The new assertions are non-tautological: they apply the perturbed matrix to a posited eigenvector and check the residual vector is zero — the residual depends on `alpha_req_q`, `lam_minus`, `r`, and the matrix construction, so a sign or factor error in any of those would yield a nonzero residual. Output transcripts show both `eigenvector residual row 0 = 0` and `row 1 = 0` (sympy) and the same with `PASS:` lines (mathematica). No collateral edits beyond the directive scope.

Note: in the Mathematica script, F1's hardcoded `eMinus = {1, q xi/(delta + xi)}` temporarily overrides the NullSpace-derived eMinus from the F3 block; the script restores `eMinus = eMinusFromNullSpace` on line 78 before the F3 overlap construction proceeds. This is the "deviation" Codex flagged in the Applied: F3 block — it is benign and preserves F3's derivation-based construction downstream.

### F2 — tautological_check

**Classification:** resolved

**What changed:**
- SymPy lines 124-149 replace the prior tautological block. Two independent paths now compute the first-order coefficients: Path 1 differentiates `F_U`/`G_U` w.r.t. `R_U`; Path 2 substitutes the eps-parametrized `(q_U_eps, eta_U_eps)` into the section-23.2 functions `F_general`/`G_general` and differentiates. The new assertions are `expect_zero("H_F cross-check (F_U vs F_general)", HF - HF_direct)` and the analogous `HG` line.
- Mathematica lines 131-150 mirror the same two-path construction with `D[...]`/`/.`, asserting `hF - hFDirect` and `hG - hGDirect`.

**Assessment:**
The cross-check is genuinely substantive: `HF` is derived from `F_U` (which was itself built by substituting `q_U(R_U), eta_U(R_U)` into `F_expected`), while `HF_direct` is derived from `F_general` (the overlap-construction route in section 23.2). The two routes share `F_expected`/`F_general` only insofar as the script established their equality in 23.2 — `F_U` then specializes via R_U-parametrized substitution while `HF_direct` re-derives from `F_general` with eps-parametrized substitution. The fact that both yield identical `4*xi*(27*delta^2 + 36*delta*xi + 11*xi^2)/((9*delta + 11*xi)*(9*delta^2 + 18*delta*xi + 11*xi^2))` (sympy output line 36-37) and `-4*xi/(9*delta + 11*xi)` confirms the two paths converge. Old `F_ratio`/`G_ratio` constructions are removed. Output transcripts confirm PASS in both engines.

### F3 — mathematica_transliteration

**Classification:** resolved

**What changed:**
- Mathematica lines 42-66 now derive `alphaReq` via `charEq = Det[mPert[alpha] - lamMinus IdentityMatrix[2]] == 0; alphaSol = Solve[charEq, alpha]`, and derive the eigenvector via `nsVec = NullSpace[mPert[alphaReq] - lamMinus IdentityMatrix[2]]`, then normalize by the first component. `r` is read off from the second component of the normalized eigenvector.
- Lines 80-106 then build `fGeneral = (a0/lamMinus) zOverlapSq sOverlapSq` and `gGeneral = (z0^2/a0) alphaReq` from the derived eigenvector (`zVec.eMinus`, etc.) rather than from a posited `r`.
- `Clear[...]` on line 35 now includes `alpha`.

**Assessment:**
The structural change is visible: `Solve[charEq, alpha]` and `NullSpace[...]` calls appear where they did not before (output line 9 prints "alpha_req (from Det = 0)" and line 10 prints "e1/e0 (from NullSpace)"), so the Mathematica path is no longer a transliteration. `fExpected`/`gExpected` are kept as a comparison target (which the directive permitted: "the comparison becomes between the two engines' derived F_general forms rather than against a shared hardcoded closed form" — the script reduces both to the same canonical form via FullSimplify and the `expectZero` confirms agreement). Output shows both `F_general - expected = 0` and `G_general - expected = 0`, demonstrating that the derivation-based construction matches the SymPy-route closed forms.

### F4 — hardcoded_result

**Classification:** resolved

**What changed:**
- SymPy lines 113-117 add a four-line comment block immediately above `F_stage18`/`G_stage19` citing `scripts/moving_throat_pde_stage035_dimensionless_normalization_locus_sympy_audit.py lines 46-58` and `scripts/moving_throat_pde_stage036_support_feasibility_frontier_sympy_audit.py lines 53-70`.
- Mathematica lines 115-119 add the analogous comment block citing `mathematica/moving_throat_pde_stage035_dimensionless_normalization_locus_mathematica_audit.wl lines 50-61` and `mathematica/moving_throat_pde_stage036_support_feasibility_frontier_mathematica_audit.wl lines 41-58`.

**Assessment:**
I cross-checked the cited upstream files:
- `scripts/moving_throat_pde_stage035_..._sympy_audit.py:54` defines `F_target = sp.simplify((9 * delta + 11 * xi) ** 4 / (81 * (1 - xi) * (9 * delta**2 + 18 * delta * xi + 11 * xi**2) ** 2))` — matches the stage-040 `F_stage18` literal exactly.
- `scripts/moving_throat_pde_stage036_..._sympy_audit.py:54` defines `G = sp.simplify(9 * xi * (xi + delta) / (9 * delta + 11 * xi))` — matches `G_stage19` exactly.
- Mathematica counterparts at the cited line ranges also contain matching `fTarget` (stage035 line 52) and `gTarget` (stage036 line 43).

The citations are concrete, non-vacuous, and land in the correct line ranges. No code change to the expressions themselves, as the directive specified. Provenance is now traceable.

## Exec log assessment

The orchestrator's `redteam/exec_logs/stage_040_*.log` files do not exist; per the verifier's path override the canonical post-fix transcripts at `scripts/output/...` and `mathematica/output/...` were used. The codex iteration loop would have halted before this verifier was dispatched if FAIL/Traceback/AssertionError/$Failed had appeared in those transcripts.

**SymPy:** exit=0 (no Traceback or AssertionError in the transcript). Notable lines:
- `eigenvector residual row 0 = 0`
- `eigenvector residual row 1 = 0`
- `F_general - expected = 0` and `G_general - expected = 0`
- `H_F cross-check (F_U vs F_general) = 0` and `H_G cross-check (G_U vs G_general) = 0`

**Mathematica:** exit=0 (transcript ends with `Stage 040 Mathematica audit passed.` and `Exit[0]` at the script level; no FAIL or $Failed lines). Notable lines:
- `alpha_req (from Det = 0) = (a0*xi*(delta + xi))/((delta + xi + q^2*xi)*z0^2)`
- `e1/e0 (from NullSpace) = (q*xi)/(delta + xi)` followed by `PASS: e1/e0 closed form`
- `PASS: eigenvector residual row 0`, `PASS: eigenvector residual row 1`
- `PASS: F_general - expected`, `PASS: G_general - expected`
- `PASS: H_F cross-check (F_U vs F_general)`, `PASS: H_G cross-check (G_U vs G_general)`

**Output freshness:** confirmed. mtimes:
- sympy script: 1779474522 < sympy output: 1779474740
- mathematica script: 1779474673 < mathematica output: 1779474749
Both `.txt` outputs were regenerated after their respective scripts were edited.

## Material-change assessment

`material_change`: false.

Rationale: every closed-form result the script computes is unchanged. `alpha_req`, `e1/e0`, `F_(q,eta)`, `G_q`, `F_U`, `G_U`, `H_F`, `H_G` all simplify to the same expressions they did pre-fix (the sympy output line 9-10 still shows `alpha_req = A0*xi*(delta + xi)/(z0**2*(delta + q**2*xi + xi))` and `e1/e0 = q*xi/(delta + xi)`, matching the audit report's table). The edits strengthen the verification surface (F1, F2, F3) and add provenance comments (F4); they do not change any derived quantity that downstream units would consume.

## Side observations (non-blocking)

- The Mathematica F1 residual block on lines 67-77 transiently rebinds `eMinus` to a posited hardcoded vector `{1, q xi/(delta + xi)}` for the residual check, then restores `eMinus = eMinusFromNullSpace` on line 78. This works correctly but is slightly clumsy — the residual check could equivalently have used `eMinusFromNullSpace` directly (the NullSpace result already equals `(1, q xi/(delta+xi))` up to FullSimplify, as the `expectZero["e1/e0 closed form", ...]` confirms). Not a finding; just a stylistic note.
- The `mBase = ...` and `zVec = ...` assignments on lines 70-71 are duplicates of lines 47-48 within the same scope. Harmless but redundant.
- Section 23.4 in SymPy defines `q_U_eps` and `eta_U_eps` then uses them in Path 2 — `F_general_eps` and `G_general_eps` correctly substitute these into the section-23.2 expressions, so the cross-check is genuinely independent at the section level.

## Verdict justification

All four findings are resolved cleanly. F1's new eigenvector residual check exercises the perturbed-matrix relationship that the previous assertion only implied algebraically. F2's tautological `series`-vs-`diff` check is replaced with a genuine two-path cross-check (`H_F` via `F_U` vs `H_F` via `F_general`) that yields identical factored expressions in both engines. F3's transliteration concern is addressed: the Mathematica script now uses `Solve[Det[...] == 0, alpha]` and `NullSpace[...]` to derive `alphaReq` and the eigenvector independently from the SymPy route, with `fGeneral`/`gGeneral` built from the derived vector. F4's provenance comments cite concrete upstream files (stage035 and stage036) at line ranges I independently verified contain the matching closed forms. Both engines exit cleanly with all `expect_zero`/`expectZero` checks passing. No regressions in the diff; no derived quantity changed, so `material_change = false`.
