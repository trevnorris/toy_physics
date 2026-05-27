---
unit_id: 097
batch: IV.1
verifier_model: claude-opus-4-7[1m]
verify_date: 2026-05-27T00:00:00Z
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 2
findings_total: 2
material_change: false
---

# Verification — unit 097

## Per-finding outcomes

### F1 — paper_misalignment

**Classification:** resolved

**What changed:**
Per orchestrator-direct application (post-user-resolution via `redteam/resolutions/batch_IV1_paper_alignment.md` Cluster A direction (a)), the SymPy script's module docstring at `scripts/moving_throat_pde_stage097_single_normalization_defect_sympy_audit.py:2-19` was expanded with a "Carry-forward annotations" block. It documents that the three `\stagefield{Checks}` items are upstream carry-ins:
- (i) static limit `eps_2=eps_4=0 -> c_pole=1/4` -> stages 091, 092, 094, 096
- (ii) `l=0`/`l=2` orthogonality -> stage 094
- (iii) minimal-module hypothesis -> stages 088 / 089 / 090

No paper edit, no new script-side assertions.

**Assessment:**
The change matches the directive's "Applied: F1" block exactly: a docstring annotation, no script-side `assert`, no paper.tex touched. Consistent with the resolved direction (a) ("Checks are upstream carry-ins" — paper card would be edited in a follow-up directive; scripts are not touched in this iteration beyond documentation). Non-tautological because no new assertion is introduced; the annotation only records the upstream provenance the auditor's table flagged. No collateral edits to assertions; the eight existing `assert` calls are untouched.

### F2 — mathematica_transliteration

**Classification:** resolved

**What changed:**
`mathematica/moving_throat_pde_stage097_single_normalization_defect_mathematica_audit.wl:33-87` rewritten. The new derivation:
- Lines 36-37 define `yhatCons[wvar_] := 3/4 + (1/4)/(1 - wvar^2/omegaQ^2)` and `kbarCons[wvar_] := k0 * yhatCons[wvar]`.
- Lines 39-40 derive `k2`, `k4` via `SeriesCoefficient[kbarCons[w], {w, 0, 2}]` and `{w, 0, 4}`.
- Lines 41-42 add two new `expectZero` PASS lines confirming `k2 - k0/(4*omegaQ^2) == 0` and `k4 - k0/(4*omegaQ^4) == 0` (series-vs-literal equivalence).
- Lines 63-64 derive `k2Target`, `k4Target` via the same series route applied to `k0Target * yhatCons[w]`.
- Line 65 derives `gamma5Target` from the `9 k0Target/(32 omegaQ^5)` form (NOT from `9 k2Target^(5/2)/k0Target^(3/2)` as the SymPy script does at line 41).
- Lines 69-77 build `R_i` reductions by defining `k0Actual = nQ * k0Target` BEFORE simplification, then series-extracting `k2Actual`, `k4Actual`, then ratioing — distinct from the SymPy approach of substituting `K0 -> N_Q * K0_target` into the already-defined `K2`, `K4`, `Gamma5`.

**Assessment:**
The algebraic route is genuinely independent of the SymPy script:
1. SymPy defines `K2 = K0/(4*Omega**2)` directly; Mathematica derives `k2` from `SeriesCoefficient[k0*(3/4 + (1/4)/(1 - w^2/omegaQ^2)), {w, 0, 2}]`. These are different paths converging to the same answer — the equivalence is then logged as a non-trivial `expectZero` (the auditor's recommendation in F2.1).
2. SymPy computes `Gamma5_target = 9 K2_target^(5/2)/K0_target^(3/2)`; Mathematica uses `9 k0Target/(32 omegaQ^5)` — algebraically distinct routes per F2.3.
3. SymPy R_i uses post-hoc substitution into pre-defined `K2`, `K4`; Mathematica builds `k0Actual = nQ*k0Target` and re-runs the series machinery, per F2.4.

The two new `expectZero` PASS lines for the series-equivalence are non-tautological: they check that `SeriesCoefficient` of a rational function of `w` recovers the closed-form coefficients. The four original `expectZero` calls remain present with their original labels (`Gamma5 - 9 K0/(32 Omega^5)`, `geometric target reduction`, `Gamma5_target - 2G/(5c^5)`, `R0/R2/R4/R5 - (N_Q - 1)`). Banner sweep confirmed: no stale "Stage 80 / STAGE 080" references in either script (grep clean). No deviation.

## Exec log assessment

**SymPy:** exit=0 (no failures, "STAGE 097 AUDIT PASSED" on line 12). Notable lines:
- `K2 = K0/(4*Omega**2)`
- `Gamma5 = 9*K0/(32*Omega**5)`
- `R0 = N_Q - 1` ... `R5 = N_Q - 1`
All four `R_i` assertions and the three `Gamma5`/target assertions hold (none flagged).

**Mathematica:** exit=0. 9 PASS lines as expected:
- `PASS: k2 series == k0/(4 omegaQ^2)` (line 6)
- `PASS: k4 series == k0/(4 omegaQ^4)` (line 8)
- `PASS: Gamma5 - 9 K0/(32 Omega^5)` (line 13)
- `PASS: geometric target reduction` (line 17)
- `PASS: Gamma5_target - 2G/(5c^5)` (line 19)
- `PASS: R0/R2/R4/R5 - (N_Q - 1)` (lines 25, 27, 29, 31)

Count matches the expected 9 (2 series-equivalence + 1 Gamma5 closed form + 1 geometric target + 1 Gamma5_target + 4 R_i). No `FAIL:` lines anywhere. Trailer "Stage 097 Mathematica audit passed."

**Output freshness:**
- `sympy.py` mtime 1779902123; `sympy.txt` mtime 1779913721 (output newer than script -> fresh).
- `mathematica.wl` mtime 1779902147; `mathematica.txt` mtime 1779913831 (output newer than script -> fresh).
Both outputs re-generated post-fix.

## Material-change assessment

`material_change`: false.

F1 only edits a docstring (no symbolic content). F2 changes the Mathematica derivation route but the four load-bearing `expectZero` results are unchanged (`9 k0/(32 omegaQ^5)`, `54 G c_s^5/(5 a^5 c^5)`, `2G/(5 c^5)`, `N_Q - 1`). The SymPy script's numerical/symbolic outputs are byte-identical to the pre-fix version (only the docstring changed). No downstream constant or definition shifted.

## Side observations (non-blocking)

- The two new series-equivalence checks (`k2 series == k0/(4 omegaQ^2)` and `k4 series == k0/(4 omegaQ^4)`) are useful for documenting the independent-derivation route, but they are essentially Mathematica re-confirming that `SeriesCoefficient` of a rational function in `w` returns the expected closed form — i.e., they validate the Mathematica algorithm rather than the physics. This is by design per the directive's F2.1 suggestion; flagging here only as a note.
- The auditor's report (line 167) noted a pre-existing renumbering artifact ("Stage 80 / STAGE 080" banner text in the docstrings). The directive's `## Applied: F2` block lists "Plus banner sweep" as part of the change; the current `.wl` and `.py` files contain no such stale banners (banner says "STAGE 097 — SINGLE NORMALIZATION DEFECT" in the `.wl` line 26; the `.py` print line 66 says "STAGE 097 AUDIT PASSED"). Sweep appears successful.
- The R0 row remains tautological (built as `k0Actual/k0Target - 1` where `k0Actual = nQ*k0Target`), as noted by the auditor. This was flagged as not load-bearing but harmless; the directive did not require a fix.

## Verdict justification

Both findings are resolved. F1 was applied per orchestrator-direct route with the user's pre-resolved Cluster A direction (a): docstring annotation only, no paper edit, no new assertions. F2's rewrite gives a genuinely independent algebraic route via `SeriesCoefficient` and adds two new non-trivial equivalence checks while preserving the four bottom-line `expectZero` calls. Both exec logs exit cleanly with the expected PASS counts (SymPy: 7 assertions hold, prints AUDIT PASSED; Mathematica: 9 PASS lines, matches expected breakdown). Outputs are fresh. No regressions, no new findings, banner sweep clean.
