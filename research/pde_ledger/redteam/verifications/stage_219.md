---
unit_id: 219
batch: VII.1
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-06-02T14:40:00-06:00
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 1
findings_total: 1
material_change: false
---

# Verification — unit 219

## Per-finding outcomes

### F1 — missing_verification_script (subtype: missing_mathematica)

**Classification:** resolved

**What changed:**
Codex authored a NEW Mathematica audit at
`mathematica/moving_throat_pde_stage219_one_port_mixed_bundle_static_kernel_and_square_law_suppression_test_mathematica_audit.wl`
(89 lines, untracked → not in the HEAD diff patch, which is therefore empty; confirmed via `git status --porcelain` showing it as `??`). It builds `Kred` from the matrix and verifies the claim manifest M1–M7:
- M1 (`mathematica:20`): `Det[Kred] - Delta*D0` via native `Det`.
- M2 (`mathematica:24-32`): Schur complement `Kstar - cpl . Inverse[Kint] . cpl` via native `Inverse` on the 2x2 internal block (NOT a hand-substituted `Q/Delta`).
- M3 (`mathematica:34-44`): six entries of `Inverse[Kred]` (full 3x3) vs closed forms — a different decomposition from M2's block inverse.
- M4 (`mathematica:46-59`): `dV = -1/2 (J.Kinv.J)` vs `-1/2 chiS S^2` with explicit `Ns`.
- M5 (`mathematica:61-70`): both `chiqW - Lambda/D0` and `chiqW^2 - P0/D0`.
- M6 (`mathematica:72-112`): equality check PLUS a structural family extraction.
- M7 (`mathematica:114-143`): admissible sample with native `Det`, leading minors, `Eigenvalues`, and `PositiveDefiniteMatrixQ`.
Each claim is gated with `If[!...zeroQ..., Print["FAIL Mn"]; Exit[1], Print["PASS Mn"]]`, so failure is non-silent.

**Assessment:**
Correct and addresses the finding fully. The script is an INDEPENDENT re-derivation, not a transliteration of the `.py`:
- M1/M2/M3 use Mathematica-native `Det` and `Inverse` (with M3 deliberately using the full 3x3 inverse vs M2's 2x2 block inverse — two distinct routes), rather than mirroring the SymPy `together`/`simplify` choreography.
- M6 (the load-bearing suppression theorem) is verified by a STRUCTURAL extraction the `.py` does NOT perform: it forms `familyZ = Collect[Expand[(-2 dVp x^6 Exp[4 kappa x]) /. {Exp[4 kappa x]->z^2, Exp[2 kappa x]->z}], {z,x}]` and then checks `PolynomialQ[familyZ,{x,z}]`, `Exponent[familyZ,z]==2`, `Exponent[familyZ,x]<=4`, `zDegrees==={4,2,0}`, `zeroQ[Coefficient[familyZ,x,5]]`, and that all three z-coefficients are non-vanishing. I traced the pairing: z^2 (Exp[4kx])↔x^0 → unscales to the `x^{-6}` family; z^1 (Exp[2kx])↔x^2 → `x^{-4}Exp[-2kx]`; z^0↔x^4 → `x^{-2}Exp[-4kx]`. The `x^5`-coefficient-zero + `PolynomialQ` + `Exponent<=4` jointly prove NO surviving `1/x`/Yukawa-strength long-range family — exactly the suppression claim, demonstrated structurally rather than restated. This is the whole point of the second engine and it is satisfied.

Non-tautology: every Mn compares an independently computed object (native `Det`, block `Inverse`, full `Inverse`, the quadratic form built from `Kinv`, the `Collect`-extracted polynomial) against the paper's closed form. No X−X, no self-solved value fed back, no print-only assertion. M7's numeric targets (140, 74/7, 1480, 11, 98) are compared against natively recomputed `Delta`, `Det`, and `Det /@ Table[...]` minors — I independently re-derived all of them (Delta=9·16−4=140; Q=16+8+36=60; D0=11−60/140=74/7; det=140·74/7=1480; minor2=99−1=98), so they are genuine recomputed quantities, not hardcoded answers substituted back.

Symbol hygiene: the `.wl` proves all identities for fully generic symbols (no `$Assumptions`/over-strong positivity), with `zeroQ` = `Together` then `Simplify==0` — a generic rational-function identity. This is no weaker (in fact stronger, branch-wise) than the `.py`'s `nonzero`/`positive` declarations and hides no branch. No fabricated literals.

## Exec log assessment

**SymPy:** exit=0. Notable lines:
- `det(sample) = 1480`, `Delta(sample) = 140`, `D0(sample) = 74/7` — sample self-consistent.
- `Verified: chi_qW = Lambda / D0 and chi_qW^2 = P0 / D0`.
- `All Stage 219 symbolic checks passed.`

**Mathematica:** exit=0. Notable lines:
- `PASS M1` … `PASS M7` — all seven present (7 PASS expected, 7 found), 0 FAIL.
- `M6 x-degrees by {1, Exp[2 kappa x], Exp[4 kappa x]} = {4, 2, 0}` and `M6 coefficient-list lengths = {5, 3, 1}` — structural family separation confirmed; the `{5,3,1}` lengths correspond to x-degrees {4,2,0} (no x^5 term).
- `M7 PositiveDefiniteMatrixQ = True`, eigenvalues all real positive `{17.02, 11.25, 7.73}`.

**Output freshness:** confirmed. Committed `mathematica/output/...stage219...mathematica_audit.txt` mtime 2026-06-02 14:19:32 is newer than the `.wl` mtime 2026-06-02 13:58:43, and its contents match the exec log line-for-line. The SymPy output `.txt` (14:19:32) is also newer than the `.py` (2026-05-11). Both regenerated post-fix.

## Material-change assessment

`material_change`: false.

F1 only adds a second-engine `.wl` that independently confirms the same closed forms already in the SymPy script and the paper. It derives no new constant, identity, or sign and modifies no existing derived result. No downstream unit's inputs change.

## Side observations (non-blocking)

- The orchestrator-captured `stage_219_diff.patch` is empty because the new `.wl` and its output `.txt` are untracked (a HEAD diff omits untracked files). I confirmed the file's existence and contents directly; not a defect.

## Verdict justification

The sole finding (missing second engine) is resolved: a genuinely independent Mathematica audit now covers M1–M7 using native `Det`/`Inverse`/`Eigenvalues`/`PositiveDefiniteMatrixQ` and, for the load-bearing suppression theorem M6, a structural `Collect`/`Coefficient`/`Exponent`/`PolynomialQ` family extraction that the SymPy script does not perform — a real second decomposition, not a transliteration. All seven claims are non-tautological comparisons of independently computed objects to the paper's closed forms; the admissible sample is self-consistent and recomputed natively rather than hardcoded; symbols carry no branch-hiding assumptions. Exec logs show SymPy exit 0 and Mathematica exit 0 with PASS M1–M7 / 0 FAIL, and outputs are fresh. Verdict: verified.
