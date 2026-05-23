---
unit_id: 059
batch: III.2
auditor_model: claude-opus-4-7[1m]
audit_date: 2026-05-22T00:00:00Z
verdict: findings
stop_cold: null
findings_count: 4
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: true
---

# Audit unit 059 red-team report

## Files reviewed

- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage059_operator_branch_residual_bounds_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage059_operator_branch_residual_bounds_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage059_operator_branch_residual_bounds_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage059_operator_branch_residual_bounds_mathematica_audit.txt`

## What the script claims to verify

Per the SymPy docstring the script asserts four claims for Stage 42 / Stage 059 ("operator-branch residual bounds"):
(1) exact lower/upper support-ratio brackets `zeta_- = A_K * Omega(Xi*Delta0)^2`, `zeta_+ = A_K * Omega(Xi*DeltaInf)^2`;
(2) exact residual-bracket definitions `R_- = zeta_req - zeta_+`, `R_+ = zeta_req - zeta_-` together with the centering identity `R_+ - R_- = zeta_+ - zeta_-`;
(3) exact coupling thresholds `Xi_fail = Pe_req/DeltaInf`, `Xi_suff = Pe_req/Delta0`, including ordering `Xi_suff > Xi_fail`, `Pe_req > Xi_fail*Delta0`, `Pe_req < Xi_suff*DeltaInf`, plus numerical confirmation that these thresholds saturate the operator-branch equation `A_K*Omega(Xi*Delta)^2 = zeta_req`;
(4) the weak-coupling expansion `Omega_Pe^2 = 1 + ((4-pi)/pi)*Pe + O(Pe^2)`, asserted via the linear coefficient.

## Assertion inventory

| # | Script | Line | Form | Anchored to claim? |
|---|---|---|---|---|
| A1 | sympy | 63 | `expect_zero("residual bracket center identity", R_hi - R_lo - (zeta_hi - zeta_lo))` | no (tautological by definitions on lines 59-60) |
| A2 | sympy | 74 | `expect_positive("Xi_suff - Xi_fail on the ordered branch", Xi_suff_ordered - Xi_fail_ordered)` | partial (trivial given declared positivity of `delta_gap`, `Pe_req`, `Delta0`) |
| A3 | sympy | 75 | `expect_positive("Pe_req - Xi_fail Delta_0", Pe_req - Xi_fail_ordered*Delta0)` | partial (trivial sign by declared positivity) |
| A4 | sympy | 76-79 | `expect_positive("Xi_suff Delta_inf - Pe_req", ...)` | partial (trivial sign by declared positivity) |
| A5 | sympy | 89-92 | `expect_zero("Omega^2 linear coefficient", coeff - (4-pi)/pi)` | yes (substantive algebra check of the explicit Omega_Pe formula) |
| A6 | sympy | 117 | `expect_close("nsolve Xi_fail saturation", Xi_fail_solved, Xi_fail_expected)` | no (zeta_req_probe is constructed as `A_K*Omega(Pe_req)^2`, so the root is `Pe_req/DeltaInf` by construction and the nsolve seed *is* that answer) |
| A7 | sympy | 118 | `expect_close("nsolve Xi_suff saturation", Xi_suff_solved, Xi_suff_expected)` | no (same construction tautology as A6) |
| A8 | wl | 62 | `expectZero["residual bracket center identity", rHi - rLo - (zetaHi - zetaLo)]` | no (mirrors A1) |
| A9 | wl | 68 | `expectPositive["Xi_suff - Xi_fail on the ordered branch", ...]` | partial (mirrors A2) |
| A10 | wl | 69 | `expectPositive["Pe_req - Xi_fail Delta_0", ...]` | partial (mirrors A3) |
| A11 | wl | 70 | `expectPositive["Xi_suff Delta_inf - Pe_req", ...]` | partial (mirrors A4) |
| A12 | wl | 98 | `expectPositive["nsolve-style Xi_fail root stayed positive", xiFailSolved]` | no (only checks sign of a value seeded at its expected positive answer) |
| A13 | wl | 99 | `expectPositive["nsolve-style Xi_suff root stayed positive", xiSuffSolved]` | no (same as A12) |
| A14 | wl | 100 | `expectApprox["FindRoot Xi_fail saturation solver", xiFailSolved, xiFailExpected, 10^-40]` | no (mirrors A6 tautology) |
| A15 | wl | 101 | `expectApprox["FindRoot Xi_suff saturation solver", xiSuffSolved, xiSuffExpected, 10^-40]` | no (mirrors A7 tautology) |
| A16 | wl | 106 | `expectZero["Omega^2 linear coefficient", Coefficient[omegaSqSeries, pe, 1] - (4-Pi)/Pi]` | yes (mirrors A5) |

## Findings

### F1 — tautological_check

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage059_operator_branch_residual_bounds_sympy_audit.py:59-63`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage059_operator_branch_residual_bounds_mathematica_audit.wl:51-62`

**What's wrong:**
The "residual bracket center identity" check (sympy line 63, wl line 62) tests
`R_hi - R_lo - (zeta_hi - zeta_lo) == 0` where the *immediately preceding lines* define
`R_lo := zeta_req - zeta_hi` and `R_hi := zeta_req - zeta_lo` (sympy 59-60; wl 51-52).
Substituting: `R_hi - R_lo - (zeta_hi - zeta_lo) = (zeta_req - zeta_lo) - (zeta_req - zeta_hi) - (zeta_hi - zeta_lo) = 0` identically. This is pure algebra on the definitions just typed; no physical or model-specific content can ever cause this to fail.

**Why this matters:**
The check is labelled as the "centering identity" of the residual bracket, suggesting the bracket has been verified to have the claimed midpoint structure. In reality the centering is true by virtue of the definitions chosen one line earlier, so a downstream reader (or paper section citing this check) would believe the residual bracket structure was verified when only its definitions were echoed back.

**Required change:**
Either remove the tautological assertion entirely, or replace it with a non-tautological check that actually exercises the bracket's *intended* role — e.g., assert that under the ordered branch `DeltaInf > Delta0` (so `zeta_hi > zeta_lo`, given Omega is monotone increasing on the relevant range), `R_hi > R_lo` (i.e., `zeta_hi - zeta_lo > 0` from independent monotonicity, not from a definitional rewrite). Minimal fix: delete the offending `expect_zero(...)` line at sympy:63 and the `expectZero[...]` at wl:62. If a substantive check is desired in its place, replace with `expect_positive("R_+ - R_- positive on ordered branch", R_hi.subs(DeltaInf, Delta0 + delta_gap) - R_lo.subs(DeltaInf, Delta0 + delta_gap))` after substituting an explicit monotone Omega such as the Stage-39 `Omega_Pe` evaluated at concrete probe values, so the positivity follows from the function form rather than from cancellation.

**Verification:**
After Codex applies, the verifier re-runs the scripts. Either the line is gone (output no longer contains "residual bracket center identity") or it is replaced with the proposed positivity assertion whose passing requires monotonicity of `Omega_Pe`, not algebraic identity.

### F2 — tautological_check

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage059_operator_branch_residual_bounds_sympy_audit.py:69-79`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage059_operator_branch_residual_bounds_mathematica_audit.wl:63-70`

**What's wrong:**
The three `expect_positive` checks for the ordering relations on the ordered branch:
- `Xi_suff - Xi_fail = Pe_req*delta_gap/(Delta0*(Delta0+delta_gap))`
- `Pe_req - Xi_fail*Delta_0 = Pe_req*delta_gap/(Delta0+delta_gap)`
- `Xi_suff*DeltaInf - Pe_req = Pe_req*delta_gap/Delta0`

Every factor in each simplified form is one of `Pe_req`, `Delta0`, `delta_gap`, which are declared `positive=True` (sympy 47-49,69; wl 41-44,64). The positivity therefore follows by sign of the symbols alone; no inequality involving Omega, A_K, the support-ratio bracket, or any operator-branch property is exercised. The same holds in the Mathematica script via `$Assumptions` declaring all these symbols positive.

**Why this matters:**
These checks are advertised in the docstring as verifying "exact coupling thresholds Xi_fail and Xi_suff" and the ordered-branch ordering. In reality they are verifying only that ratios of declared-positive symbols are positive — they would pass for any pair of formulas with the same monomial sign structure, regardless of whether those formulas correctly encoded the saturation thresholds.

**Required change:**
The thresholds themselves (`Xi_fail = Pe_req/DeltaInf`, `Xi_suff = Pe_req/Delta0`) are introduced as definitions on sympy:65-66 / wl:53-54 and never tested against an independent inversion. Either (a) drop these three `expect_positive` lines (sympy 74-79 and wl 68-70), since they assert nothing physical; or (b) add a substantive check that the *thresholds themselves* are correct: pick a concrete probe Omega (`Omega_Pe`), an *independent* zeta_req value (not constructed via Omega(Pe_req)), use `nsolve`/`FindRoot` to invert `A_K*Omega(Pe)^2 = zeta_req` for `Pe`, and then verify `Xi_fail*DeltaInf` equals that root and `Xi_suff*Delta0` equals that root. That actually tests the threshold formulas. Minimal in-scope fix is (a): delete sympy:74-79 and wl:68-70.

**Verification:**
Output lines `Xi_suff - Xi_fail on the ordered branch = ...`, `Pe_req - Xi_fail Delta_0 = ...`, `Xi_suff Delta_inf - Pe_req = ...` are removed, or replaced with a check that uses an explicit Omega so the sign depends on the function shape rather than declared positivity.

### F3 — insufficient_verification

**Severity:** high
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage059_operator_branch_residual_bounds_sympy_audit.py:94-118`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage059_operator_branch_residual_bounds_mathematica_audit.wl:77-101`

**What's wrong:**
The nsolve/FindRoot "saturation" checks (sympy 105-118, wl 88-101) are constructed circularly:

1. Line sympy:102 / wl:85 defines `zeta_req_probe := A_K_probe * Omega_Pe(Pe_req_probe)^2`.
2. Lines sympy:105-110 / wl:88-92 then solve `A_K_probe*Omega_Pe(Xi*DeltaInf_probe)^2 - zeta_req_probe = 0` for `Xi`, seeded at `Xi_fail_expected = Pe_req_probe/DeltaInf_probe`.
3. The check `expect_close(Xi_fail_solved, Xi_fail_expected)` (sympy 117) is satisfied because `Xi = Pe_req_probe/DeltaInf_probe` is an *exact* root by construction: substituting gives `Omega_Pe(Pe_req_probe)^2 - Omega_Pe(Pe_req_probe)^2 = 0`.

The seed equals the answer, the answer is forced by the construction of `zeta_req_probe`, and the solver therefore returns the seed to within numerical precision regardless of whether the threshold formulas `Pe_req/DeltaInf` and `Pe_req/Delta0` have any physical meaning. The same circular pattern repeats in the Mathematica script (line wl:85 constructs `zetaReqProbe` from `Omega_Pe` at `Pe_req`, then wl:88-97 inverts the same equation seeded at the answer). The two `expectPositive` lines wl:98-99 are even weaker — they only check the sign of a quantity seeded as positive.

**Why this matters:**
The docstring claims this section verifies "exact coupling thresholds Xi_fail and Xi_suff". A naive reader sees `nsolve(... ) ≈ Pe_req/DeltaInf` and concludes the thresholds have been numerically validated against the operator-branch saturation equation. They have not — the test would pass even if `Xi_fail` and `Xi_suff` were defined as arbitrary positive numbers, because `zeta_req_probe` is constructed *via the same Omega evaluated at those numbers' product with the corresponding Delta*. There is no independent target zeta value.

**Required change:**
Break the circularity by giving `zeta_req_probe` an *independent* numerical value (not derived from `Omega_Pe(Pe_req_probe)`). For example, set `zeta_req_probe = sp.Rational(2, 5)` (any value in the operator-branch range), then independently solve `A_K_probe*Omega_Pe(Pe)^2 - zeta_req_probe = 0` for `Pe` via nsolve to obtain `Pe_star`, and then assert `expect_close("Xi_fail*DeltaInf = Pe_star", Xi_fail_solved*DeltaInf_probe, Pe_star)` and `expect_close("Xi_suff*Delta0 = Pe_star", Xi_suff_solved*Delta0_probe, Pe_star)`. This verifies the *threshold relations* (`Xi_fail = Pe_star/DeltaInf`, `Xi_suff = Pe_star/Delta0`) against an externally-chosen zeta. Apply the analogous restructuring to the Mathematica script. If restructuring is out of scope, delete the four offending `expect_close`/`expectApprox` lines (sympy 117-118, wl 100-101) and the two `expectPositive` lines wl:98-99, since they currently assert nothing physical.

**Verification:**
Either (a) the `nsolve Xi_fail saturation diff`/`nsolve Xi_suff saturation diff` lines no longer appear in the output (deletion), or (b) the script defines `zeta_req_probe` as a literal rational (e.g., 2/5) prior to constructing `Pe_star` via solver, and the new pair of `expect_close` assertions on `Xi_fail*DeltaInf` and `Xi_suff*Delta0` appear in the output.

### F4 — mathematica_transliteration

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage059_operator_branch_residual_bounds_mathematica_audit.wl` (whole file)
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage059_operator_branch_residual_bounds_sympy_audit.py` (whole file)

**What's wrong:**
The Mathematica script is a line-by-line port of the SymPy script rather than an independent derivation. Three corresponding sections:

(a) Definitions of A_K, zeta_lo, zeta_hi, R_lo, R_hi:
- sympy 52-60: `A_K = (kappa + pi**2/4)/(kappa + y**2)`, `zeta_lo = A_K*Omega(Xi*Delta0)**2`, `zeta_hi = A_K*Omega(Xi*DeltaInf)**2`, `R_lo = zeta_req - zeta_hi`, `R_hi = zeta_req - zeta_lo`.
- wl 47-52: `aK = (kappa + Pi^2/4)/(kappa + y^2)`, `zetaLo = aK*omegaFn[capitalXi*delta0]^2`, `zetaHi = aK*omegaFn[capitalXi*deltaInf]^2`, `rLo = zetaReq - zetaHi`, `rHi = zetaReq - zetaLo`. Identical algebraic forms, identical ordering, identical intermediate names (modulo CamelCase).

(b) Probe values for the nsolve/FindRoot block:
- sympy 95-100: `kappa_probe = 2`, `y_probe = 1`, `Delta0_probe = Rational(3,5)`, `delta_gap_probe = Rational(2,5)`, `DeltaInf_probe = Delta0_probe + delta_gap_probe`, `Pe_req_probe = Rational(7,10)`.
- wl 78-83: `kappaProbe = 2`, `yProbe = 1`, `delta0Probe = 3/5`, `deltaGapProbe = 2/5`, `deltaInfProbe = delta0Probe + deltaGapProbe`, `peReqProbe = 7/10`. Identical numerical choices.

(c) Omega_Pe explicit formula and series order:
- sympy 83-87: `Omega_Pe = pi*Pe*(2*Pe*exp(Pe) + pi)/((4*Pe**2 + pi**2)*(exp(Pe) - 1))`, then `series(Omega_Pe**2, Pe, 0, 2)`.
- wl 72-76: `omegaPe = Pi*pe*(2*pe*Exp[pe] + Pi)/((4*pe^2 + Pi^2)*(Exp[pe] - 1))`, then `Series[omegaPe^2, {pe, 0, 1}]`. Same closed-form Omega, same series degree.

Both scripts choose the same seeds for nsolve/FindRoot (`Xi_fail_expected`, `Xi_suff_expected`) and use the same precision (70-digit working precision, 80-digit `N`). The Mathematica script also reproduces the SymPy script's `(4 - pi)/pi` linear-coefficient target verbatim (wl:106 vs sympy:91-92). No independent re-derivation is performed — no alternative path through the Stage-39 series (e.g., expanding `Omega_Pe` differently then squaring, or factoring out the `(exp(Pe) - 1)^(-1)` and using a different limit technique) is taken. The Mathematica script's only structural difference from the SymPy script is the rename `Xi → capitalXi`, `Delta0 → delta0`, etc., and the use of `FindRoot` instead of `nsolve`.

**Why this matters:**
The two-engine policy requires both engines to derive the result independently from the physical premises. A line-by-line port can only catch typos in one engine, not algebraic errors common to both. If the Stage-39 `Omega_Pe` formula had a sign error in its construction or a wrong combinatorial factor in the exponential, both scripts would echo the same wrong formula and both would pass.

**Required change:**
The Mathematica script should derive at least one of the four numbered claims by a path structurally different from the SymPy approach. Minimal independent re-derivation that fits the existing scope: in the Mathematica script replace the explicit closed-form `omegaPe` (wl:72-75) with a derivation that uses `Series[]` on the *integrand* of the Stage-39 representation of `Omega(Pe)` (or, if that representation is not present, derive the small-Pe series of `omegaPe` from `Limit[D[omegaPe, pe], pe -> 0]` for the linear coefficient rather than via `Series[..., {pe, 0, 1}]`). The verifier check at wl:106 then becomes `Limit[D[omegaPe^2, pe], pe -> 0] - (4 - Pi)/Pi`, which exercises a different Mathematica code path. Apply this single independence break and leave the rest of the structural ordering intact (the policy requires *independent derivation of the result*, not full reorganization).

**Verification:**
Output line `Omega^2 linear coefficient = ...` continues to print 0, but the underlying Mathematica computation now uses `Limit[D[...], pe -> 0]` rather than `Coefficient[Series[...], pe, 1]`. Cross-engine agreement on the (4-Pi)/Pi value is preserved.

## Independent-derivation check (Mathematica)

See F4 above. The Mathematica script reproduces the SymPy script's algebraic choreography line-for-line (same `A_K`, same `zeta_lo`/`zeta_hi`, same `R_lo`/`R_hi`, same probe numerics 2,1,3/5,2/5,7/10, same closed-form Omega_Pe, same series order, same nsolve seed strategy). It is a transliteration, not an independent re-derivation.

## Engine cross-check

The two engines agree at the symbolic level shown in their saved outputs:

| Quantity | SymPy | Mathematica |
|---|---|---|
| zeta_- | `(kappa + pi**2/4)*Omega(Delta0*Xi)**2/(kappa + y**2)` | `((kappa + Pi^2/4)*Omega[capitalXi*delta0]^2)/(kappa + y^2)` |
| R_- | `(zeta_req*(kappa+y**2) - (4*kappa+pi**2)*Omega(DeltaInf*Xi)**2/4)/(kappa+y**2)` | `zetaReq - ((4*kappa + Pi^2)*Omega[capitalXi*deltaInf]^2)/(4*(kappa + y^2))` (algebraically equal) |
| Xi_fail | `Pe_req/DeltaInf` | `peReq/deltaInf` |
| Xi_suff - Xi_fail | `Pe_req*delta_gap/(Delta0*(Delta0 + delta_gap))` | `(deltaGap*peReq)/(delta0*(delta0 + deltaGap))` |
| Omega_Pe^2 series | `Pe*(-1 + 4/pi) + 1` | `1 + pe*(-1 + 4/Pi)` |
| Omega^2 linear coefficient | 0 | 0 |
| nsolve diff Xi_fail | 1.8e-72 | 0``49.85 |
| nsolve diff Xi_suff | 6.0e-72 | 0``49.63 |

No engine disagreement. (This agreement is unsurprising given F4 — the engines compute the same thing through the same intermediate steps.)

## Verdict justification

Both engines run cleanly and their outputs agree, but the agreement is largely on tautological or definition-rewriting assertions (F1, F2) and on a circularly-constructed numerical saturation check (F3). The single substantive check — the linear coefficient `(4-pi)/pi` of `Omega_Pe^2` — is the same algebra performed twice by the same `Series[]`/`series()` path (F4), so the second engine does not provide independent confirmation. The audit therefore identifies four findings: two tautology cases, one insufficient-verification (the nsolve/FindRoot saturation circularity), and one transliteration. No `UNFIXABLE` or `CRITICAL_DOWNSTREAM` flag — the fixes are local and do not change any quoted constant that downstream units would use; the only genuine result (the `(4-pi)/pi` linear coefficient) is already correct, and the proposed F4 fix preserves it.

## Self-test notes

Variable-independence: no proposed `diff`/`D` checks introduce variable-mismatch issues; F4's proposed `Limit[D[omegaPe^2, pe], pe -> 0]` is well-defined since `omegaPe` explicitly depends on `pe`. Parity: no symmetric-domain integrals are involved. Trivial-case pre-check: for F1 the residue `R_hi - R_lo - (zeta_hi - zeta_lo) = (zeta_req - zeta_lo) - (zeta_req - zeta_hi) - (zeta_hi - zeta_lo) = 0` collapses identically (confirmed tautological); for F3, substituting `zeta_req_probe = A_K*Omega(Pe_req)^2` into `A_K*Omega(Xi*Delta)^2 - zeta_req_probe = 0` yields `Omega(Xi*Delta)^2 = Omega(Pe_req)^2` whose monotone branch is exactly `Xi*Delta = Pe_req`, i.e., the seeded answer (confirmed circular). Path specifications: no `missing_verification_script` findings, so no `.py`/`.wl` path-guessing exposure.
