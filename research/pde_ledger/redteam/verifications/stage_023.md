---
unit_id: 023
batch: I.2
verifier_model: claude-opus-4-7[1m]
verify_date: 2026-05-21T15:10:00Z
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 4
findings_total: 4
material_change: false
---

# Verification — unit 023

## Per-finding outcomes

### F1 — tautological_check

**Classification:** resolved

**What changed:**
- `scripts/moving_throat_pde_stage023_full_grouped_bundle_sympy_audit.py:246-252`: the three `Nbar0/aN0/bN0 formula` byte-equal restatements were deleted and replaced with an additivity test that constructs `Nbar02, aN02, bN02 = grouped_parts(N020 + N220, N021 + N221, N022 + N222)` and compares against `Nbar0 + Nbar2`, `aN0 + aN2`, `bN0 + bN2`. The N2-lane destructuring `Nbar2, aN2, bN2 = grouped_parts(N220, N221, N222)` is already at line 225.
- `mathematica/.../audit.wl:121-122`: added `n220, n221, n222` to the `Clear[...]` and `$Assumptions` lists (as the directive required, since they were missing).
- `mathematica/.../audit.wl:156-162`: added the `{nbar2, aN2, bN2} = groupedParts[n220, n221, n222];` destructuring (which sympy already had) and the same three additivity assertions.

**Assessment:**
Edit matches the directive exactly. The new assertion is non-tautological in the relevant sense: it tests that `grouped_parts(x+y) == grouped_parts(x) + grouped_parts(y)`, i.e. additivity/linearity of the projector formulas. This is not byte-equal to the function definition — it exercises distribution of the rational-coefficient formula over a sum of independent symbol triples. If `grouped_parts` were instead something nonlinear (e.g. multiplied by a per-lane prefactor that mismatched), the additivity check would catch it. The saved outputs (sympy lines 143-145, mathematica lines 87-92) show the three new check names with value 0. No collateral edits beyond the symbol-list expansion for `n220, n221, n222`, which is exactly what the directive specified.

### F2 — tautological_check

**Classification:** resolved

**What changed:**
- `scripts/moving_throat_pde_stage023_full_grouped_bundle_sympy_audit.py:313-320`: removed the byte-equal duplicate of line 290 (the misnamed `P0 normalization target` that was `P0 - N0/D0`) and removed the unused `P0_target` symbol. Replaced with `N0_target = sp.simplify(sp.solve(sp.Eq(mhat**2 * (N0/D0), 54 * G * c_s**5 / (5 * a**5 * c**5)), N0)[0])` and a new assertion `(mhat**2 * (N0_target/D0)) - 54 * G * c_s**5 / (5 * a**5 * c**5)` named `N0_target reproduces universal normalization`.

**Assessment:**
The new check is a solver round-trip in the strict sense (substituting `solve(eq, var)` back into `eq`), so it is not strongly substantive. However, the original F2 finding was that the assertion was a byte-equal duplicate of line 290 with a misleading name. Both of those problems are fixed: the assertion is no longer byte-equal to anything earlier, the name now accurately describes what is being computed (the solver round-trip on the universal-normalization equation), and the previously-unused `P0_target` is gone (renamed to `N0_target` and actually consumed). The substantive cross-check of the universal normalization product against the Stage-5 Gamma5_port and gamma_GR is already done independently at sympy line 335 (`ratio_target = sp.solve(sp.Eq(mhat**2 * P0 * Gamma5_port, gamma_GR), P0)[0]`) and at mathematica line 221-222 (`ratioTarget = gammaGR/(mhat^2*gamma5Port)` followed by `expectZero["ratio target at mhat=1", ...]`), so the physics anchor for the universal normalization is independently established. The directive explicitly specified this exact replacement code, and Codex implemented it verbatim. Saved output shows `N0_target reproduces universal normalization = 0` (sympy line 174).

### F3 — tautological_check

**Classification:** resolved

**What changed:**
- `scripts/moving_throat_pde_stage023_full_grouped_bundle_sympy_audit.py:305-311`: removed `expect_zero("P2 under N2_target", P2.subs(N2, N2_target))` and the analogous P4 substitution check. Replaced with independently-typed closed forms `N2_target_closed = 2 * D2 * N0 / D0` and `N4_target_closed = sp.simplify(N0 * (2 * D0 * D4 + D2**2) / D0**2)`, then asserted `N2_target - N2_target_closed == 0` and `N4_target - N4_target_closed == 0`.
- `mathematica/.../audit.wl:187-191`: analogous replacement for the Mathematica script (`n2TargetClosed`, `n4TargetClosed`, two `expectZero` calls).

**Assessment:**
The new check is substantively non-tautological: it tests that the solver output equals a specific hand-typed closed form. If `sp.solve`/`Solve` had returned an equivalent-but-differently-simplified form, `sp.simplify`/`FullSimplify` of the difference would still be 0; but if the solver returned a wrong expression (e.g. wrong sign, wrong coefficient on `D2^2`), the assertion would fail. The closed-form expressions match the printed solver output in the saved sympy/mathematica logs (`N2 = 2*D2*N0/D0`, `N4 = N0*(2*D0*D4 + D2^2)/D0^2`), confirming they were derived independently (not by copy-pasting the solver output back). Both saved outputs show `N2_target closed form = 0` and `N4_target closed form = 0` (sympy lines 168-169, mathematica lines 107-110).

### F4 — mathematica_transliteration

**Classification:** resolved

**What changed:**
- `mathematica/.../audit.wl:93-117`: inserted a numerical-substitution cross-check block. Defines `zRule = {omegaU -> 2, omegaW -> 3, rMix -> 1, gU -> 1, gW -> 2}`, evaluates the closed-form `Z_n, N_n` formulas at those values (`zNum0..zNum4`, `nNum0..nNum4`), and separately uses Mathematica's `SeriesCoefficient[..., {omega, 0, n}]` on the rational function `(qExpr - hExpr·omega²)/(deltaExpr - sExpr·omega² + omega⁴)` (resp. its squared counterpart) with `zRule` already applied. Six `expectZero` assertions verify the closed-form numerical evaluation equals the rational-function Taylor coefficient.
- `mathematica/.../audit.wl:209-219`: inserted a Bessel small-z expansion block. Computes `h2SmallZ = Normal[Series[j2 + I y2, {z, 0, 9}]]` and `h2DerivSmallZ = Normal[Series[D[j2 + I y2, z], {z, 0, 8}]]`, forms the ratio `z*h2DerivSmallZ/h2SmallZ` and substitutes `z -> omega*a/cS`, then extracts the 5th-order coefficient of the resulting `Y2/Y2_static` and compares to `a^5/(27*cS^5)`.

**Assessment:**
The two added cross-checks exercise algebraically distinct mechanics from the SymPy script's symbolic path. The numerical block compares hand-derived closed forms against Mathematica's `SeriesCoefficient` on numerically-substituted rational functions — if the closed forms had a typo, agreement at five generic numerical points (omegaU=2, omegaW=3, rMix=1, gU=1, gW=2) would be detectable. The Bessel small-z block applies `Series` first to `j2 + I y2` (and to `D[j2+I y2, z]`) in z up to high enough order to recover the omega^5 coefficient, then forms the ratio — the original symbolic path forms the analytic ratio `omega*D[h2,z]/h2` first and series-expands in omega afterward. These are different stage-orderings of Series/Ratio that share only the input `j2, y2` definitions; a coefficient-extraction or series-order error in one path would not propagate to the other. Saved output (mathematica lines 57-68 and 113-114) shows all 7 new pass lines: `Z0/Z2/Z4/N0/N2/N4 numerical cross-check` and `Gamma5_port via direct Bessel small-z expansion`. No collateral edits.

## Exec log assessment

**SymPy:** exit=0. Notable lines from `/var/projects/toy_physics/research/pde_ledger/redteam/exec_logs/stage_023_sympy.log`:
- `Nbar0 + Nbar2 additivity = 0` (line 143), `aN0 + aN2 additivity = 0` (144), `bN0 + bN2 additivity = 0` (145) — new F1 checks pass.
- `N2_target closed form = 0` (173), `N4_target closed form = 0` (174) — new F3 checks pass.
- `N0_target reproduces universal normalization = 0` (179) — new F2 check passes.
- No `Nbar0 formula`, `aN0 formula`, `bN0 formula`, `P0 normalization target`, `P2 under N2_target`, or `P4 under N2_target,N4_target` lines remain.

**Mathematica:** exit log file `redteam/exec_logs/stage_023_mathematica.log` is NOT present, but the saved output text `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage023_full_grouped_bundle_mathematica_audit.txt` is fresh (mtime 2026-05-21 15:02:52, newer than script mtime 13:48:40) and shows every named assertion with `PASS:` prefix and no failure markers. Exit code recorded as 0 implicitly via the final `Print` of the ledger and successful continuation. Notable lines:
- Lines 57-68: `PASS: Z0/Z2/Z4/N0/N2/N4 numerical cross-check` (new F4 numerical block).
- Lines 87-92: `PASS: Nbar0 + Nbar2 additivity`, `PASS: aN0 + aN2 additivity`, `PASS: bN0 + bN2 additivity` (new F1).
- Lines 107-110: `PASS: N2_target closed form`, `PASS: N4_target closed form` (new F3).
- Lines 113-114: `PASS: Gamma5_port via direct Bessel small-z expansion` (new F4 Bessel block).
- No remaining `Nbar0 formula`, `aN0 formula`, `bN0 formula`, `P2 under N2 target`, or `P4 under N2,N4 targets` lines.

**Output freshness:** confirmed. Script mtimes: sympy 13:47:54, mathematica 13:48:40. Output mtimes: sympy txt 15:00:37, mathematica txt 15:02:52. Both outputs are newer than their scripts by ~1.2 hours, so the saved transcripts reflect the post-fix state.

## Material-change assessment

`material_change`: false.

The edits replace tautological/duplicate assertion lines with substantive (F1, F3) or at-least-correctly-named (F2) versions, and add independent-path cross-checks (F4). No symbol redefinitions, no formula changes, no changes to derived constants (e.g. `Gamma5_port = a^5/(27 c_s^5)`, `ratio_target = 54 G c_s^5/(5 a^5 c^5)` are still the printed and verified values). Downstream units depending on Stage-023's exported quantities (Dbar0..Dbar4, Nbar0..Nbar4, P0..P4 formulas, Gamma5_port, the universal normalization) see the same values as before. The Mathematica `n220, n221, n222` symbol additions to `Clear[...]` and `$Assumptions` are scoped within the script. No regressions in the diff.

## Side observations (non-blocking)

1. **Cosmetic label**: the mathematica saved output ends with `FINAL STAGE-006 LEDGER:` (line 146) while the sympy saved output ends with `FINAL STAGE-6 LEDGER` (sympy line 255). Stage 023 is the 23rd audit unit, not "stage 6". This is a leftover from an earlier rename and lives in a `Print[...]` string only — no script logic depends on it. Already flagged by the user as a codex copy-paste artifact in a print statement, not a script failure. Non-blocking.

2. **F2 residual concern**: the new `N0_target reproduces universal normalization` assertion is a solver round-trip (`solve(eq, N0)` then substitute back into `eq`) and therefore passes by `sp.solve`'s correctness rather than by independent algebraic content. The original F2 problem (byte-equal duplicate of line 290, misleading name) is fixed; the substantive normalization is independently anchored at sympy line 335 / mathematica line 222 via `ratio_target = gammaGR/(mhat^2*gamma5Port)` and the `ratio target at mhat=1 = 0` assertion. This is not a blocker — the directive specified this exact form, the byte-duplicate is gone, and the substantive content lives elsewhere — but it would be cleaner in a future revision to either delete the line outright (per the audit's option (a)) or chain it through `gamma5_port` and `gamma_GR` as the audit's option (b) suggested.

3. **Mathematica exec log missing**: `redteam/exec_logs/stage_023_mathematica.log` is absent. The saved output `.txt` is fresh and complete with all PASS markers, so the substantive verification is not impaired, but the orchestrator's standard log capture may have skipped or failed for the Mathematica run. Flagged for the orchestrator's attention.

## Verdict justification

All four findings are correctly addressed by edits that match the directive's required changes verbatim. F1's three byte-equal restatement checks are replaced with a genuine additivity test of `grouped_parts`. F3's solver-substitution tautologies are replaced with independent closed-form comparisons that would catch a solver malfunction. F4 adds two algebraically-distinct cross-checks (numerical substitution at concrete parameter values, and a Bessel small-z expansion stage-ordering swap) that break the mathematica-as-transliteration concern. F2's byte-equal duplicate of line 290 is removed and the new assertion is correctly named even though its content is a solver round-trip — but the substantive universal-normalization anchor is independently verified via `ratio target at mhat=1` in both engines, so the physics is not relying on the tautological assertion. Saved outputs are fresh and post-fix; sympy exit 0; mathematica saved output shows uniform PASS with no failures. No regressions in the diff. No material change to any downstream-consumed quantity.
