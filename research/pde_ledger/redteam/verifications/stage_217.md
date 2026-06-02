---
unit_id: 217
batch: VI.1
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-06-02T00:00:00-06:00
verdict: verified
sympy_exit: n/a
mathematica_exit: 0
findings_resolved: 2
findings_total: 2
material_change: false
---

# Verification — unit 217

## Per-finding outcomes

### F1 — paper_misalignment (value_mismatch)

**Classification:** resolved

**What changed:**
F1 was USER-RESOLVED via direction (a): the script's `162` (= 3⁴·2) is correct
and arithmetically forced; the published card/appendix `179` and notes `230` were
typos. Codex corrected the prose targets (`paper/stages/stage_217.tex`,
`paper/appendices/stage_appendix_part06.tex`, and the stage-217 notes), which the
orchestrator has already reviewed and confirmed correct and isolated. The
scripts-only verifier scope here is to confirm the SymPy `.py` was NOT changed.

**Assessment:**
Confirmed. The SymPy script
`scripts/moving_throat_pde_stage217_..._sympy_audit.py` is unmodified in the
worktree: it does not appear in `git status --porcelain` and `git diff` against
HEAD is empty. Its assertion `if bezout_lift != 162:` (line 138) and the projected
`if bezout_proj != 750:` (line 167) are intact and unchanged. No script value was
touched for F1 — the paper now matches the already-correct script. Correct, no
collateral.

### F2 — missing_verification_script (missing_mathematica)

**Classification:** resolved

**What changed:**
Codex created the new
`mathematica/moving_throat_pde_stage217_..._mathematica_audit.wl` (untracked, `??`
in git status), independently asserting M1–M6. Output `.txt` mtime (12:32:41) is
newer than the `.wl` mtime (12:16:58) — freshly regenerated post-fix.

**Assessment — INDEPENDENCE (load-bearing, not a rubber-stamp):**

*M1 — distinct denominator-clearing route.* The `.wl` computes the cleared
derivative numerators SEPARATELY first
(`scaledDerivativeNumerators = FullSimplify[Cancel[Together[2 sqrtDelta S^(3/2) D[phi,#]]], …]`,
lines 111–114), prints them, then asserts `row[[2]] - row[[3]] == 0` where
`row[[3]] = 2 sqrtDelta Mv + Lv` (lines 110, 117–125). This is the
directive-prescribed independent `Together`/`Cancel` denominator-clearing route on
`D[phi,v]`, NOT the `.py`'s single fused
`simplify(expand(2*sqrt(Delta)*S**(3/2)*diff(Phi,r) - (2*Mr*sqrt(Delta)+Lr)))`
(`.py` lines 92–107). Different decomposition. Not a transliteration.

*Degrees — native primitive.* `totalDegree` uses
`Max[Total /@ CoefficientRules[Expand[poly], vars][[All,1]]]` (lines 49–54) — the
directive-named native primitive — versus the `.py`'s `sp.Poly(...).total_degree()`
(`.py` lines 121–125, 152–155). Different route.

*M3 computed = 162, not hardcoded.* `liftedProduct = Times @@ liftedDegrees`
(line 142) where `liftedDegrees` are the MEASURED total degrees of `{Fr,Fs,Ft,Fu,FΔ}`;
asserts both `liftedDegrees === {3,3,3,3,2}` (line 147) AND `liftedProduct == 3^4*2`
(line 148). Log: "M3 lifted Bezout product = 162" computed from `{3,3,3,3,2}`. No
`179`/`230` appears anywhere in the `.wl` (grep: NONE FOUND). 162/750 occur only as
the literal forms `3^4*2` / `5*5*5*6` compared against the computed products.

*M4 computed = 750.* `projectedProduct = Times @@ projectedDegrees` from measured
`{5,5,5,6}`; asserts tuple `=== {5,5,5,6}` and product `== 5*5*5*6` (lines 156–163).
Log confirms 750. Computed, agrees with `.py`.

*M2/M5/M6.* Faces via `Table[Delete[axes,i]]` with `Length==5` + per-face cardinality
4; M5 diagonal-isotropic substitution residuals (`L_v(diag)-2 Klin M_v` and M/L at
gradient ratios); M6 fivefold-symmetric equal-mix residuals. All independent
substitutions, all non-tautological.

Exec: exit 0, 30 PASS, 0 FAIL, covering M1(4) + M2(2) + M3(2) + M4(2) + M5(12) +
M6(8) = 30. No FAIL lines. Genuinely independent derivation confirmed.

## Exec log assessment

**SymPy:** n/a (no SymPy re-run for this unit; F1 is paper-only, `.py` unchanged).

**Mathematica:** exit=0. Notable lines:
- `M3 lifted Bezout product = 162` / `PASS: M3 lifted product equals 3^4*2`
- `M4 projected one-chart Bezout product = 750` / `PASS: M4 projected product equals 5*5*5*6`
- `M1 stationary numerator identity r = 0` … `PASS` (r,s,t,u)
- `All Stage 217 Mathematica audit checks passed.` / `# exit_code: 0`

**Output freshness:** confirmed — `.txt` (2026-06-02 12:32:41) is newer than `.wl`
(2026-06-02 12:16:58); saved output and exec log agree (both 30 PASS, 0 FAIL).

## Material-change assessment

`material_change`: false. No derived script result changed. The F1 paper/notes typo
fix makes the prose match the already-correct script `162`; the new `.wl`
independently confirms the same computed `162`/`750`. No downstream unit relies on a
changed value.

## Side observations (non-blocking)

The orchestrator's note anticipated 31 PASS lines; the actual count is 30 (M1:4,
M2:2, M3:2, M4:2, M5:12, M6:8), consistent between the exec log and the saved
output `.txt`. This is a minor off-by-one in the expectation note, not a coverage
gap — all six manifest items M1–M6 are independently asserted and pass. Not blocking.

## Verdict justification

Both findings are resolved. F1: the SymPy `.py` is confirmed unchanged (empty
`git diff`, `bezout_lift != 162` intact at line 138), so the paper-only typo fix is
correctly scoped. F2: the new `.wl` is a genuinely independent derivation — M1 via a
distinct `Cancel[Together[D[phi,v]]]` denominator-clearing route, degrees via native
`CoefficientRules`, with the `162`/`750` bounds COMPUTED as `Times @@` of measured
degree tuples (no hardcoded 179/230 present) — exits 0 with 30 non-tautological
PASSes and 0 FAILs covering all of M1–M6. Verdict: verified.
