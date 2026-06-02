---
unit_id: 213
batch: VI.1
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-06-02T00:00:00Z
verdict: verified
sympy_exit: n/a
mathematica_exit: 0
findings_resolved: 1
findings_total: 1
material_change: false
---

# Verification — unit 213

## Per-finding outcomes

### F1 — missing_verification_script (subtype: missing_mathematica)

**Classification:** resolved

**What changed:**
Codex created the previously-absent Mathematica audit at the registered path
`mathematica/moving_throat_pde_stage213_four_coordinate_mixed_simplex_and_support_cardinality_4_gate_mathematica_audit.wl`
(354 lines). It implements `expectZero`/`expectZeroUnder`/`expectTrue` helpers with a
`stripConditional` guard (lines 27, 38-49) per the directive, and verifies each claim-manifest
item M1–M9 with explicit labeled residuals. The orchestrator's independent re-run
(`exec_logs/stage_213_mathematica.log`) exits 0 with 62 `PASS:` lines and 0 `FAIL:`.

**Assessment:**
Edit is correct and addresses the finding. Every manifest claim is represented as an explicit
non-tautological check (PASS counts by item: M1=3, M2=8, M3=20, M4=4, M5=6, M6=1, M7=11, M8=1,
M9=8). No collateral edits: the new `.wl` is the only file touched (the captured
`stage_213_diff.patch` is empty because the file is untracked/new; the SymPy `.py` is unchanged).

The non-tautology requirement is met because the `.wl` derives — rather than re-asserts — the
load-bearing results, using native Mathematica primitives via a different decomposition than the
SymPy `.py`. See the independent-derivation check below.

## Independent-derivation check (Mathematica vs SymPy)

**Verdict: confirmed independent (not a transliteration).** Three corresponding sections:

1. **M3 gradient optimum.** SymPy *posits* the closed form and asserts properties of it:
   `avec_grad = sp.Matrix([ki/Kgrad, ...])` then `expect_zero("gradient-optimal slope value",
   (avec_grad.T*kvec)[0] - Kgrad)` (`.py` lines 133, 141). The `.wl` instead *derives* the optimum:
   it forms the Lagrangian `lagObjective - mu(lagNorm - 1)`, takes `D[...]` per variable, assembles
   `lagSystem`, `Solve`s over Reals, and `SelectFirst`s the positive-`mu` branch (`.wl` lines
   151-167); then checks `lagValue - kNorm == 0`. It also enumerates all 15 `Subsets[Range[4]]`
   active supports and proves full-support dominance (`.wl` 169-189) — content with no SymPy
   analogue. Genuinely different route.

2. **M5 cross-leverage bound.** SymPy never bounds `w_Σ`; it only substitutes three discrete points
   (`.py` 190-192). The `.wl` derives the bound by constrained optimization:
   `Maximize[{wSigma[lagVars], lagVars.lagVars == 1 && And @@ Thread[lagVars >= 0]}, lagVars, Reals]`
   (`.wl` 210-219), recovering value 3 and maximizer (1/2,1/2,1/2,1/2). Independent primitive
   (`Maximize`) and a strictly stronger claim than the `.py`.

3. **M7 discriminant coefficients.** SymPy types A…J by hand and subtracts a hand-assembled
   `Delta_sharp` (`.py` 230-252). The `.wl` extracts each coefficient from the expanded numerator
   via nested `Coefficient` (`coeff3`, `.wl` 51-60, 251-254) and validates them against an
   *independently constructed* matrix kernel `Outer[Times, kVec, kVec] - 2 H0 Hblock` (`.wl`
   255-267, 280). This is exactly the anti-transliteration extraction the directive mandated, not a
   re-typing of A…J.

Shared scaffolding (face parametrizations, the symmetric Hblock, the τ ratio form) is unavoidable
common problem structure, not ported choreography — the verification *routes* differ at every
load-bearing claim.

## Exec log assessment

**SymPy:** exit=n/a. SymPy side unchanged this iteration; no re-run required (the auditor already
confirmed the `.py` output fresh with EXIT_CODE 0).

**Mathematica:** exit=0. Notable lines:
- `PASS: M3 maximum value from positive branch minus ||k||` (residual `0`) — Lagrange-derived optimum equals ‖k‖.
- `M5 Maximize result = {3, {x1 -> 1/2, x2 -> 1/2, x3 -> 1/2, x4 -> 1/2}}` then `PASS: M5 constrained maximum value - 3`.
- `PASS: M7 reconstructed extracted Delta sharp equals numerator` (residual `0`) — coefficient-extraction route closes.
- `All Stage 213 Mathematica audit checks passed.` / `# exit_code: 0`.

62 `PASS:`, 0 `FAIL:`. (The orchestrator's "63" includes the trailing summary line, which is not a
`PASS:`-prefixed assertion.)

**Output freshness:** confirmed. Output mtime `2026-06-02 10:28:50` is newer than the `.wl` script
mtime `2026-06-02 10:19:14` — regenerated post-fix.

## Material-change assessment

`material_change`: false. The edit only *adds* a second-engine verification of identities already
established by the unchanged SymPy script; no derived result that a downstream unit depends on was
altered.

## Side observations (non-blocking)

- The SymPy `.py` banners still read "STAGE 196 …" (lines 42, 402), a stale stage label carried from
  a sibling. Cosmetic only — not in scope for this finding and not a verification blocker.
- M9 verifies the three (r,s,t)-patch face collapses plus the jkl face; the directive's optional
  Section-VI interval-gate checks were (correctly, per directive) not reproduced as they are trivial
  ordering relations.

## Verdict justification

The single finding F1 (`missing_mathematica`) is resolved: the `.wl` now exists at the registered
path, runs under `math -script`, exits 0 with 62 PASS / 0 FAIL, and represents every claim-manifest
item M1–M9 as an explicit, non-tautological check. Critically, the script is a genuinely independent
derivation — Lagrange `Solve` for the gradient optimum, `Maximize` for the cross-leverage bound, and
`Coefficient`-extraction plus an independent matrix-kernel cross-check for the discriminant — rather
than a line-by-line port of the SymPy algebra. No regressions (SymPy unchanged; new file is the only
addition). Verdict: verified.
