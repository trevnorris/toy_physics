---
unit_id: 003
batch: I.1
verifier_model: claude-opus-4-7[1m]
verify_date: 2026-05-21T00:00:00Z
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 4
findings_total: 4
material_change: false
---

# Verification — unit 003

## Per-finding outcomes

### F1 — insufficient_verification

**Classification:** resolved

**What changed:**

Both engines now derive `D_eff` (Section I) and the one-mode dispersion (Section II) from the Euler-Lagrange equations by explicit elimination of the matter amplitudes, replacing the previous formula-vs-hand-expansion-of-same-formula check with a real elimination check.

SymPy (`scripts/moving_throat_pde_stage003_bdg_sympy_audit.py`):

- Lines 130-156 (inside `axisymmetric_matrix_kernel_audit`): introduces frequency-space amplitudes `Qa, QL, Xa, Xb`, substitutes the ansatz `qa, qL, xa, xb -> {Qa,QL,Xa,Xb} * exp(-I*omega*t)` into each `EL_*.lhs`, divides by the phase, simplifies, solves the two matter equations `EL_xa_f = 0, EL_xb_f = 0` for `{Xa, Xb}`, substitutes back into the wall equations, and reads off the 2x2 derived matrix as `D_derived = -[[∂EL_qa_red/∂Qa, ∂EL_qa_red/∂QL], [...]]`. The new assertion `expect_zero("derived D0_eff vs Deff", sp.simplify(D_derived - Deff))` appears at line 156, before the existing `subbanner("I.2 ...")`.
- Lines 215-237 (inside `one_mode_pole_shift_audit`): adds a single-mode Lagrangian `L_one = M*q'^2/2 - K*q^2/2 + x'^2/2 - varpi2*x^2/2 + g*q*x`, derives `EL_q, EL_x` via `sp.euler_equations`, applies the freq-space ansatz, solves `EL_x_f = 0` for `X`, substitutes back, multiplies by `(varpi2 - omega^2)/Q` (and prepends a `-` to match the existing `dispersion` sign convention), substitutes `omega^2 -> w2`, and asserts `expect_zero("derived dispersion vs (K - M w2)(varpi2 - w2) - g^2", derived_dispersion - dispersion)`.

Mathematica (`mathematica/moving_throat_pde_stage003_bdg_mathematica_audit.wl`):

- Lines 96-117: rebuilds `dEff` from `mMat, kMat, cMat, oMat`, then independently derives it by ansatz substitution into the four `el{Qa,QL,Xa,Xb}` equations (constructed in the F2 block), solves the matter pair for `{Xa, Xb}`, substitutes back, and extracts the 2x2 derived matrix via `Coefficient[..., {Qa, QL}]`. Asserts `expectMatrixZero["derived D0_eff vs Deff", dDerived - dEff]`.
- Lines 156-169: adds the single-mode Lagrangian `lOne`, derives `elQ, elX`, applies the ansatz, solves for `X`, substitutes, multiplies by `(varpi2 - omega^2)/Q`, substitutes `omega^2 -> w2`, and asserts `expectZero["derived dispersion vs (k - m w2)(varpi2 - w2) - g^2", (dispersionDerived /. omega^2 -> w2) - dispersion]`.

**Assessment:**

Both engines now exercise an actual elimination step rather than restating the formula. In SymPy, `EL_qa_red` (after substituting the solved `Xa, Xb`) is built entirely from the EL operators, which are themselves built from `Lred`. `Deff` is built from the independent hand-typed matrices `Mmat, Kmat, Cmat, Omat`. The two paths share `wa, wb, c1a, c1b, c2a, c2b, Maa, MaL, MLL, Kaa, KaL, KLL` as symbols but not as algebraic expressions; the test will fail if either path has a sign or factor error. Same story in Mathematica: the EL equations come from `lRed` via the explicit time-derivative chain `timeD[D[lAlg, vqa]] - D[lAlg, qa0]` (which is mathematically equivalent to the EL operator), and `dEff` comes from `LinearSolve[oMat - omega^2 IdentityMatrix[2], Transpose[cMat]]` — a Mathematica-native primitive that does not rely on the same algebraic path as SymPy. The check is non-tautological.

The SymPy output shows `derived D0_eff vs Deff = [[0,0],[0,0]]` at line 13-16 and `derived dispersion vs (K - M w2)(varpi2 - w2) - g^2 = 0` at line 89. Mathematica output shows `PASS: derived D0_eff vs Deff` at line 18 and `PASS: derived dispersion vs (k - m w2)(varpi2 - w2) - g^2` at line 29.

The minor deviations codex notes (using `D_derived = -sp.Matrix(...)` and `-sp.together(...)` to absorb SymPy's `euler_equations` sign convention, and pre-defining `dispersion = (k - m w2)(varpi2 - w2) - g^2` before comparing) are sign-bookkeeping required by the engines' EL conventions, not structural deviations. A wrong sign on `g*q*x` in the Lagrangian, or a wrong sign on `c1a*qa*xa` etc. in `lRed`, would change the derived dispersion / `D_derived` by a sign on the coupling-squared term and break the assertion.

### F2 — tautological_check

**Classification:** resolved

**What changed:**

The Mathematica file replaces the original kinetic-coefficient-extraction tautology (lines 61-72 in the pre-edit file) with explicit Euler-Lagrange equation residual checks (lines 61-94 in the current file):

- Lines 62-67: `lRed = lRed + (-1/2 kaa qa^2 - kaL qa qL - 1/2 kLL qL^2 + 1/2 D[xa, t]^2 - 1/2 wa^2 xa^2 + 1/2 D[xb, t]^2 - 1/2 wb^2 xb^2 + c1a qa xa + c1b qa xb + c2a qL xa + c2b qL xb)`. Codex's `## Applied: F2` deviation note flagged that the original multi-line `lRed = ...` assignment was being parsed as the kinetic terms only, with subsequent lines treated as orphaned expression statements. Adding the parenthesised potential / matter / coupling terms back into `lRed` completes the Lagrangian.
- Lines 68-83: defines `lAlg` (the Lagrangian rewritten in `qa0, vqa, aqa` etc. symbols), a `timeD` total-time-derivative helper, and the back-substitution rule `backEL`. Computes `elQa, elQL, elXa, elXb` as `(timeD[D[lAlg, v*]] - D[lAlg, *0]) /. backEL`, which is mathematically `d/dt(∂L/∂q') - ∂L/∂q` — the standard EL operator.
- Lines 85-94: four `expectZero` calls matching SymPy's Section I.1 assertions: `qa equation`, `qL equation`, `xa equation`, `xb equation`, each asserting `elQ* - (expected Maa qa'' + MaL qL'' + Kaa qa + KaL qL - c1a xa - c1b xb)` is zero.

**Assessment:**

Codex's claim that the original `lRed` only evaluated the kinetic terms is plausible and is the only consistent explanation for why the original tautological kinetic-coefficient checks passed without revealing the missing terms. (In Wolfram script mode, the first complete expression after `lRed =` is `1/2 maa D[qa, t]^2 + maL D[qa, t] D[qL, t] + 1/2 mLL D[qL, t]^2` — a syntactically complete RHS — and the subsequent lines starting with `- ...` are then top-level expressions whose values are discarded.) Adding the rest via `lRed = lRed + (...)` is a correct, minimal fix.

After the fix, `D[lAlg, vqa] = maa vqa + maL vqL`, and `timeD[D[lAlg, vqa]] = maa aqa + maL aqL`. `D[lAlg, qa0] = -kaa qa0 - kaL qL0 + c1a xa0 + c1b xb0`. Therefore `elQa = maa aqa + maL aqL + kaa qa0 + kaL qL0 - c1a xa0 - c1b xb0`, which after `/. backEL` becomes `maa qaFun''[t] + maL qLFun''[t] + kaa qaFun[t] + kaL qLFun[t] - c1a xaFun[t] - c1b xbFun[t]` — exactly the expected expression, so the assertion `elQa - expected = 0` is non-trivial. A wrong sign or coefficient in any term of `lRed` would break this check; a sign error in `timeD`/`backEL` would break it; a sign error in the expected expression would break it. None of those checks would have been caught by the original kinetic-coefficient tautology.

Mathematica output lines 9-16 show the four PASS lines (`qa equation`, `qL equation`, `xa equation`, `xb equation`), matching SymPy's Section I.1 PASS lines exactly. The three tautological `qa kinetic coefficient` / `qL kinetic coefficient` / `qa-qL mixed kinetic coefficient` PASS lines from the pre-edit transcript are gone.

### F3 — insufficient_verification

**Classification:** resolved

**What changed:**

SymPy (lines 363-383): adds `I00_21c, I00_22c, I21c_22c, N00, N21c, N22c` integral definitions and six new `expect_zero` assertions covering the previously-missing cross-integrals and norms. SymPy Section IV now reports 10 zero/norm assertions (was 4), matching the directive's expected count.

Mathematica (lines 243-260): adds the same integral definitions (`i0021c, i0022c, i21c22c, norm00, norm21c, norm22c`) and then, per F4's later consolidation, replaces the entire list of per-pair `expectZero` calls with a single `expectMatrixZero["spherical harmonic overlap matrix - identity", overlap - IdentityMatrix[4]]` over the 4-element basis `{y00, y20, y21c, y22c}`. The single matrix check subsumes all six cross-integrals and all four norms required by F3.

**Assessment:**

In SymPy, the new assertions are non-tautological: each `expect_zero` integrates a literal expression `Y_lm * Y_l'm' * sin(theta)` over the sphere and asserts the result is zero (or one for norms). The Y_lm are defined as explicit `cos`/`sin` polynomials of `theta, phi`. If the normalization of any Y_lm were wrong, or if a `cos(2 ph)` were a `sin(2 ph)`, the corresponding norm or cross-integral would not vanish. Output lines 186-195 show ten distinct `... = 0` entries.

In Mathematica, the 4x4 overlap matrix check is strictly stronger than F3 required: it verifies ALL 16 entries of the matrix (six off-diagonal pairs, four diagonal norms, plus repeats of the off-diagonal). Output line 68 shows the full 4x4 zero matrix `{{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}}` followed by `PASS: spherical harmonic overlap matrix - identity`. This is structurally non-tautological for the same reason as the SymPy version.

### F4 — mathematica_transliteration

**Classification:** resolved

**What changed:**

Mathematica `.wl` Sections III (lines 197-217) and IV (lines 243-260) are restructured to differ from the SymPy script's algebraic path:

Section III: introduces a diagonal-matrix construction `dP2 = DiagonalMatrix[{d20, d21, d22}]`, then series-expands element-wise via `Map[Normal[Series[#, {omega, 0, 4}]] &, dP2, {2}]`, extracts the omega^2 coefficient matrix `d2coeffRaw`, weights it as `d2coeffMat = DiagonalMatrix[{d2coeffRaw[[1,1]], 2 d2coeffRaw[[2,2]], 2 d2coeffRaw[[3,3]]}]`, and defines `d2Bar = Tr[d2coeffMat]/5`, plus `a2` and `b2` as explicit Cartesian-component projections. Representation-theoretic basis matrices `T0, Ta, Tb` are declared but the script uses the explicit-component projection for `a2, b2` (acceptable per the "Acceptable" path in F4: the construction path goes through a diagonal-matrix series-extraction rather than three independent scalar `Series` calls).

Section IV: replaces the previous four per-pair `expectZero` calls with a single `expectMatrixZero` over the 4x4 overlap matrix `Table[Integrate[Integrate[yList[[i]] yList[[j]] dOmega, {ph,0,2 Pi}], {th,0,Pi}], {i,1,4}, {j,1,4}]` against `IdentityMatrix[4]`. This is structurally not a transliteration of SymPy's pair-by-pair iteration.

**Assessment:**

Section III: the construction now goes via a single `dP2` matrix and a matrix-valued `Series`/`Coefficient` operation, rather than the SymPy script's three scalar `sp.series(D20,...)`, `sp.series(D21,...)`, `sp.series(D22,...)` calls. The `d2Bar = Tr[d2coeffMat]/5` is genuinely structurally different from SymPy's `(d220 + 2*d221 + 2*d222)/5`: it uses the trace of a weighted diagonal matrix rather than a hand-typed linear combination. (`Tr[DiagonalMatrix[{a, 2b, 2c}]]/5 == (a + 2b + 2c)/5`, so the resulting expression is mathematically identical to SymPy's, which is exactly what F4 requires for the "Acceptable" path: same answer, different construction.) The output values for `d2bar, a2, b2` (Mathematica lines 53-55) match SymPy's pretty-printed forms (sympy lines 154-171) modulo notation. `T0, Ta, Tb` are defined but unused for projection — this is a minor side observation (not a finding), since the directive's "Preferred" path was T-matrix projection but the "Acceptable" path was the diagonal-matrix trace approach codex took.

Section IV: the 4x4 overlap-matrix check is **strictly stronger** than the per-pair pattern in SymPy and cannot have been written by transliterating it. The construction `yList = {y00, y20, y21c, y22c}; overlap = Table[Integrate[..., {i, 1, 4}, {j, 1, 4}]]; expectMatrixZero[..., overlap - IdentityMatrix[4]]` is a Mathematica-native structural pattern. An off-by-one in the SymPy script's pair selection (e.g. omitting a pair) could not be propagated by copy here, since this routine touches all 16 entries.

A weakness of codex's Section III rewrite: the redefined `d20s, d21s, d22s` (lines 208-210) are reassigned from `dP2s[[i,i]]`, which is mathematically identical to the original `Expand[Normal[Series[d20, {omega, 0, 4}]]]`. So the printed forms remain the same and the Mathematica printer still emits `D20 = ...`, `D21 = ...`, `D22 = ...` lines that look identical to SymPy's. The non-transliteration is in the `dP2`-construction step and the `Tr[d2coeffMat]/5` projection, which is the part F4 was demanding. The isotropy substitutions on lines 226-230 still mirror SymPy's, but those are user-facing substitutions and F4 did not require them to differ.

## Exec log assessment

**SymPy:** exit=0 (orchestrator-confirmed; `redteam/exec_logs/stage_003_sympy.log` is absent from the directory listing). The saved output `scripts/output/moving_throat_pde_stage003_bdg_sympy_audit.txt` shows uniform completion with all `= 0` lines and a final `FINAL STAGE-3 LEDGER` block, indicating no `raise AssertionError` was hit. Notable lines:

- Line 9-12: `qa equation = 0`, `qL equation = 0`, `xa equation = 0`, `xb equation = 0` (Section I.1 EL).
- Line 13-16: `derived D0_eff vs Deff = [[0,0],[0,0]]` (F1 new check).
- Line 33-36: `D0_eff - manual form = [[0,0],[0,0]]` (existing check, preserved).
- Line 89: `derived dispersion vs (K - M w2)(varpi2 - w2) - g^2 = 0` (F1 new check).
- Lines 186-195: ten Section IV harmonic assertions, all `= 0` (F3 new checks plus original four).

**Mathematica:** exit=0 (orchestrator-confirmed; `redteam/exec_logs/stage_003_mathematica.log` is absent from the directory listing). The saved output `mathematica/output/moving_throat_pde_stage003_bdg_mathematica_audit.txt` ends with `Stage 003 Mathematica audit passed.` and the script's terminal `Exit[0]` would not have been reached if any `expectZero`/`expectMatrixZero` had hit `fail[]` (which calls `Exit[1]`). Notable lines:

- Lines 10-16: `PASS: qa equation`, `qL equation`, `xa equation`, `xb equation` (F2 EL checks replacing tautologies).
- Lines 17-18: `derived D0_eff vs Deff = {{0,0},{0,0}}` and `PASS` (F1 new check).
- Lines 28-29: `derived dispersion vs (k - m w2)(varpi2 - w2) - g^2 = 0` and `PASS` (F1 new check).
- Line 53-55: `d2bar`, `a2`, `b2` symbolic forms matching SymPy's.
- Lines 68-69: `spherical harmonic overlap matrix - identity = {{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}}` and `PASS` (F3+F4 consolidated check).
- No `FAIL:` lines anywhere.

**Output freshness:** SymPy script mtime 2026-05-21 00:54, output mtime 2026-05-21 11:26 — fresh. Mathematica script mtime 2026-05-21 01:03, output mtime 2026-05-21 11:50 — fresh. Both outputs are newer than their respective scripts, confirming post-fix execution.

**Exec log files missing:** `redteam/exec_logs/stage_003_sympy.log` and `redteam/exec_logs/stage_003_mathematica.log` were not captured in the directory. The orchestrator's confirmation of exit=0 plus the fresh post-fix `.txt` outputs (which would have been truncated mid-section if any assertion had failed) provides sufficient evidence. Informational; orchestrator may want to standardize log capture.

## Material-change assessment

`material_change`: false.

The fixes change *how* the same claims are verified, not the claims themselves:

- F1 derives `D_eff` and the one-mode dispersion from EL elimination, but the resulting symbolic forms match the original `Deff` and `dispersion` exactly (since both old checks were correct in the first place — the audit's concern was that the verification was insufficient, not that the result was wrong).
- F2 replaces tautological kinetic-coefficient extractions with substantive EL residual checks; the Lagrangian `lRed` is now correctly defined (it was silently incomplete in the original), but the downstream `dEff` / `manual` / series / Section II / Section III / Section IV results are derived from independent matrix definitions and are unchanged. Codex's note that the original `lRed` only carried kinetic terms is a real bug in the prior verification scaffolding, but no downstream derivation in the .wl actually used the missing terms — the downstream values came from `mMat, kMat, cMat, oMat`, not from `lRed`. So the engine's reported outputs (`D0_eff`, `a2`, `b2`, etc.) match SymPy's just as they did before.
- F3 adds new assertions on the same harmonic basis as before; no change to the basis or to any downstream consumer.
- F4 restructures Section III's construction path through a diagonal matrix; the resulting `d2bar, a2, b2` are symbolically identical to the previous values.

No symbolic results that downstream units (004+) depend on have changed. Section I's `Deff`, Section II's roots and shifts, Section III's `d2bar/a2/b2`, and Section IV's selection rule all evaluate to the same expressions as before. Downstream units may still inherit the `upstream_stale: true` flag per the orchestrator's policy, but the stage-003 verification result itself is unchanged.

## Side observations (non-blocking)

1. The Mathematica Section III defines `T0, Ta, Tb` representation-theoretic projection matrices (lines 204-206) but never uses them — the actual `a2, b2` come from `(2 d2coeffRaw[[1,1]] - d2coeffRaw[[2,2]] - d2coeffRaw[[3,3]])/10` and `(d2coeffRaw[[2,2]] - d2coeffRaw[[3,3]])/2`. F4's "Preferred" path was to use `T0/Ta/Tb` projection, but codex took the "Acceptable" path (different construction order via the `dP2` diagonal matrix). The unused `T0, Ta, Tb` definitions are harmless residue; not a finding.

2. Mathematica's `lRed` initially only carried kinetic terms because of a multi-line parse quirk: the assignment `lRed = 1/2 maa D[qa,t]^2 + maL D[qa,t] D[qL,t] + 1/2 mLL D[qL,t]^2` is syntactically complete after the first multi-token line, and subsequent lines starting with `- 1/2 kaa qa^2 - ...` were being parsed as standalone top-level expressions (no-op effect). The fix `lRed = lRed + (...)` uses an explicit parenthesised continuation to add the remaining terms. This is now correct. Not a finding (codex's `## Applied: F2` deviation block already disclosed this). However, the same parse quirk could exist in other `.wl` files in the ledger that use multi-line bare-`+`/`-` continuation without enclosing parentheses; orchestrator may want to spot-check.

3. The SymPy script's Section II dispersion derivation uses `derived_dispersion = sp.expand(-sp.together(EL_q_red * (varpi2 - omega**2) / Q))` — the leading `-` and the `omega**4 -> w2**2`, `omega**2 -> w2` substitutions are sign- and variable-bookkeeping required by SymPy's `euler_equations` convention. The output line 89 confirms the check passes. Not tautological (a wrong sign on `g*q*x` in `L_one` would change the `g^2` term to `-g^2` in the derived dispersion, which would not match the pre-defined `dispersion = (K - M*w2)*(varpi2 - w2) - g^2`).

4. Diff cosmetic: the post-edit `.wl` has tab-indented lines (lines 61-117, 150-171) mixed with the surrounding two-space indentation. This is a whitespace-only style artifact; Mathematica is whitespace-insensitive here. Not blocking.

5. The `redteam/exec_logs/stage_003_sympy.log` and `stage_003_mathematica.log` files are missing from the exec_logs directory. Orchestrator confirmed both engines exit 0; fresh `.txt` outputs corroborate.

## Verdict justification

All four findings are resolved with substance. F1 introduces real elimination steps in both engines that derive `D_eff` and the one-mode dispersion from the EL equations, producing assertions that would catch sign errors in the Lagrangian, the EL operator convention, or the matter-equation solve — none of which the original `Deff - manual` check could detect. F2 replaces three tautological kinetic-coefficient extractions in Mathematica with four substantive EL residual checks matching SymPy's, and incidentally fixes a silent prior bug where `lRed` was carrying only the kinetic terms due to a multi-line parse quirk. F3 adds the missing cross-integrals and norms in both engines, closing the selection-rule's coverage gap. F4 restructures Mathematica Section III through a diagonal-matrix series construction with a `Tr`-based projection (different algebraic path from SymPy's three scalar series calls), and consolidates Section IV's per-pair pattern into a single 4x4 overlap-matrix-equals-identity check (structurally not a transliteration). All assertions PASS in fresh post-fix outputs; both engines exit 0; no regressions visible in the diff. Material change: false — symbolic results unchanged, only the verification machinery is strengthened.
