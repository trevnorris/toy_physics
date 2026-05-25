---
unit_id: 078
batch: III.4
created_at: 2026-05-22T00:00:00Z
findings_count: 4
stop_cold: null
applied: true
applied_at: 2026-05-23T05:28:30Z
findings_applied: 4
findings_blocked: 0
verification_status: pending
---

# Codex directive — unit 078

Apply each finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

If a finding's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named.

Do NOT run python or mathematica. Only edit files.

Do NOT touch paper.tex, notes/, or any prose documents. The red-team only modifies scripts.

## F1 — tautological_check

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage078_family1_branch_verdict_sympy_audit.py:41-44`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage078_family1_branch_verdict_mathematica_audit.wl:61-62`

**Issue:**
The two ordering inequalities (`Pe_suff_chi < Pe_fail_chi` and `Pe_suff_J < Pe_fail_J`) reduce algebraically to `Theta_fail_coeff < Theta_suff_coeff`; `Theta_chi`, `Theta_J`, `lambda_mu`, and `Pe_req` all cancel. The script cannot detect any corruption of the Stage-60 Theta values it claims to be verifying. Strengthen the assertion set so it depends non-trivially on `Theta_chi` and `Theta_J` (this is fully resolved together with F4 by adding the branch-verdict checks; F1's purpose here is to record that the *original* assertions remain in place and that the new ones supplement, not replace, them).

**Required change:**
Do *not* delete the existing assertions at SymPy lines 41-44 or Mathematica lines 61-62. Keep them as a structural sanity check, and rely on F4's branch-verdict assertions to actually exercise the Theta values. No direct edit is required for F1; it is satisfied once F4 is applied. Mark F1 as `Applied: F1` with `summary: no direct edit; supplanted by F4 branch-verdict assertions which depend on all four Theta coefficients` and `files_changed: none`.

**Verification:**
After F4 lands, the SymPy script contains at least three new assertions referencing both `Pe_suff_chi`/`Pe_suff_J` and `Pe_fail_chi`/`Pe_fail_J`, and the Mathematica script mirrors them. Verifier confirms via `redteam exec-sympy 078` / `redteam exec-mathematica 078` that the new lines run and exit 0.

## Applied: F1

- files_changed: none
- summary: no direct edit; supplanted by F4 branch-verdict assertions which depend on all four Theta coefficients
- deviation: none

## F2 — hardcoded_result

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage078_family1_branch_verdict_sympy_audit.py:26-29`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage078_family1_branch_verdict_mathematica_audit.wl:37-40,57-60`

**Issue:**
The four Theta coefficients enter as bare floats with no provenance. The Mathematica `expectApprox` targets at lines 57-60 are bit-identical copies of the SymPy printed output, so the cross-engine "check" is `sympy_value == sympy_value`. Add provenance comments on the SymPy side and (for the Mathematica side, see F3) replace the literal-decimal `expectApprox` targets with independently-computed numerics.

**Required change:**

(a) In `scripts/moving_throat_pde_stage078_family1_branch_verdict_sympy_audit.py`, replace lines 26-29

```python
Theta_chi = sp.Float("4.06863235008162") * lambda_mu**2
Theta_J = sp.Float("0.927552032539308") * lambda_mu**2
Theta_fail = sp.Float("3.62605617972939e-4") * Pe_req
Theta_suff = sp.Float("4.21495341569977e-2") * Pe_req
```

with

```python
# Stage-77 family-1 Theta extraction:
#   Theta_w^(chi) / lambda_mu^2 = 4.06863235008161550927... (chi^2 weighted datum)
#   Theta_w^(J)   / lambda_mu^2 = 0.92755203253930797184... (Jensen floor)
# Source: scripts/output/moving_throat_pde_stage077_family1_theta_extraction_sympy_audit.txt:22-23
Theta_chi = sp.Float("4.06863235008162") * lambda_mu**2
Theta_J = sp.Float("0.927552032539308") * lambda_mu**2
# Stage-75 family-1 threshold window:
#   Theta_fail / Pe_req = 0.00036260561797293886969 (= (37 cosh(111 sqrt(5)/5)
#       + 111 sqrt(5) sinh(111 sqrt(5)/5)/5) / (136900 (-1 + sqrt(5) sinh(111 sqrt(5)/5)/3
#       + cosh(111 sqrt(5)/5))))
#   Upsilon_suff / Pe_req = 4.2149534156997728721 ; Theta_suff = Upsilon_suff / 100
# Source: scripts/output/moving_throat_pde_stage075_family1_threshold_window_sympy_audit.txt:30,34
Theta_fail = sp.Float("3.62605617972939e-4") * Pe_req
Theta_suff = sp.Float("4.21495341569977e-2") * Pe_req
```

(b) The Mathematica fix for the `expectApprox` literal-target half of this finding is covered by F3 below (which replaces the copied decimals with independently-computed `Sinh`/`Cosh` values). Do not attempt a separate edit here.

**Verification:**
The four Theta assignment lines in the SymPy script now carry comments naming the upstream output file (`stage077_..._audit.txt:22-23`, `stage075_..._audit.txt:30,34`) and (for `Theta_fail`) the symbolic closed form. The four numeric values are unchanged.

## Applied: F2

- files_changed:
  - `scripts/moving_throat_pde_stage078_family1_branch_verdict_sympy_audit.py`
- summary: added provenance comments for the four SymPy Theta coefficients while preserving the numeric values
- deviation: none

## F3 — mathematica_transliteration

**Target:** `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage078_family1_branch_verdict_mathematica_audit.wl:32-60`

**Issue:**
The .wl file is a one-to-one port of the .py file: identical hardcoded decimals, identical variable choreography (`thetaChiCoeff` ↔ `Theta_chi`, etc.), identical arithmetic, and `expectApprox` targets that are literal SymPy output. Rewrite the constant-loading block so that Mathematica derives the four coefficients from their symbolic closed forms, and base the `expectApprox` checks on those independent computations.

**Required change:**

Replace lines 37-40 of the .wl file

```mathematica
thetaChiCoeff = SetPrecision[4.06863235008162, 20];
thetaJCoeff = SetPrecision[0.927552032539308, 20];
thetaFailCoeff = SetPrecision[3.62605617972939*^-4, 20];
thetaSuffCoeff = SetPrecision[4.21495341569977*^-2, 20];
```

with

```mathematica
(* Independent re-derivation of the four window coefficients from
   their Stage-75 symbolic closed forms (no SymPy input).             *)
thetaFailSym = (37 Cosh[111 Sqrt[5]/5] + (111 Sqrt[5]/5) Sinh[111 Sqrt[5]/5])
               / (136900 (-1 + (Sqrt[5]/3) Sinh[111 Sqrt[5]/5]
                              + Cosh[111 Sqrt[5]/5]));
thetaSuffSym = 100 thetaFailSym * (4.21495341569977*^-2 / 3.62605617972939*^-4);
(* The chi^2 and Jensen-floor Theta values are recorded numerically in
   stage077 output; we adopt them at high precision but verify their
   ratio chi:J matches the Stage-77 ratio.                              *)
thetaChiCoeffNum = ToExpression["4.0686323500816155092718546246574670820527`40"];
thetaJCoeffNum   = ToExpression["0.92755203253930797183993260663904217023`40"];
thetaChiCoeff  = thetaChiCoeffNum;
thetaJCoeff    = thetaJCoeffNum;
thetaFailCoeff = N[thetaFailSym, 30];
thetaSuffCoeff = N[thetaSuffSym, 30];
```

Then replace the four `expectApprox` lines at lines 57-60

```mathematica
expectApprox["Pe_suff^(chi) numeric check", peSuffChi, SetPrecision[96.528524726438575954, 25], 10^-12];
expectApprox["Pe_fail^(chi) numeric check", peFailChi, SetPrecision[11220.544162625905301, 25], 10^-9];
expectApprox["Pe_suff^(J) numeric check", peSuffJ, SetPrecision[22.006222633075413597, 25], 10^-12];
expectApprox["Pe_fail^(J) numeric check", peFailJ, SetPrecision[2558.0189234920526360, 25], 10^-10];
```

with

```mathematica
(* Independent targets computed in Mathematica from the symbolic
   closed form (thetaFailSym) and the chi/J Stage-77 numerics.   *)
peSuffChiTarget = N[thetaChiCoeffNum / thetaSuffSym, 30];
peFailChiTarget = N[thetaChiCoeffNum / thetaFailSym, 30];
peSuffJTarget   = N[thetaJCoeffNum   / thetaSuffSym, 30];
peFailJTarget   = N[thetaJCoeffNum   / thetaFailSym, 30];
expectApprox["Pe_suff^(chi) numeric check", peSuffChi, peSuffChiTarget, 10^-12];
expectApprox["Pe_fail^(chi) numeric check", peFailChi, peFailChiTarget, 10^-9];
expectApprox["Pe_suff^(J) numeric check",   peSuffJ,   peSuffJTarget,   10^-12];
expectApprox["Pe_fail^(J) numeric check",   peFailJ,   peFailJTarget,   10^-10];
```

The remaining lines (the `Element[...]` assumptions, the `thetaChi`/`thetaJ`/`thetaFail`/`thetaSuff` definitions at lines 42-45, and the `peSuffChi`/`peFailChi`/`peSuffJ`/`peFailJ` computations at lines 47-50) stay unchanged so that the same `peSuffChi` etc. quantities are still being tested.

**Verification:**
The .wl file now contains an explicit `Cosh[111 Sqrt[5]/5]`/`Sinh[111 Sqrt[5]/5]` expression for `thetaFailSym`, and the four `expectApprox` calls reference `peSuffChiTarget`/`peFailChiTarget`/`peSuffJTarget`/`peFailJTarget` rather than literal decimals. The verifier re-runs `redteam exec-mathematica 078`, the script exits 0, and the four `... numeric check diff = ...` output lines still print and PASS (with a small but nonzero diff reflecting the independent computation, not the bit-identical zero of a self-compare).

## Applied: F3

- files_changed:
  - `mathematica/moving_throat_pde_stage078_family1_branch_verdict_mathematica_audit.wl`
- summary: replaced the Mathematica literal coefficient and numeric-check target block with symbolic/independently computed coefficients and targets
- deviation: parenthesized thetaFailSym and removed the extra 100 factor from thetaSuffSym so the Mathematica coefficient stays on the Stage-75 Theta_suff scale

## F4 — insufficient_verification

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage078_family1_branch_verdict_sympy_audit.py:44` (append after)
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage078_family1_branch_verdict_mathematica_audit.wl:62` (append after)

**Issue:**
The unit is the "family-1 branch verdict" but never asserts a verdict comparing the chi datum and the Jensen floor. Add three substantive assertions (per engine) that depend non-trivially on `Theta_chi` and `Theta_J` so that the script genuinely tests the branch verdict it is named for.

**Required change:**

(a) Append to `scripts/moving_throat_pde_stage078_family1_branch_verdict_sympy_audit.py` after line 44:

```python

# --- Branch verdict (chi vs Jensen-floor) ---------------------------------
# Since Theta_J < Theta_chi and both Pe-thresholds share positive denominators,
# the Jensen-floor success threshold and failure ceiling both lie below the
# chi-datum's corresponding thresholds. The two windows nest with overlap.
if not (sp.N(Pe_suff_J) < sp.N(Pe_suff_chi)):
    raise AssertionError(
        "Expected Jensen-floor success threshold to lie below chi-datum's "
        "(requires Theta_J < Theta_chi)"
    )
if not (sp.N(Pe_fail_J) < sp.N(Pe_fail_chi)):
    raise AssertionError(
        "Expected Jensen-floor failure ceiling to lie below chi-datum's "
        "(requires Theta_J < Theta_chi)"
    )
if not (sp.N(Pe_suff_chi) < sp.N(Pe_fail_J)):
    raise AssertionError(
        "Expected chi-datum success threshold below Jensen-floor failure "
        "ceiling (window overlap)"
    )
print("Pe_suff_J < Pe_suff_chi  :", sp.N(Pe_suff_J) < sp.N(Pe_suff_chi))
print("Pe_fail_J < Pe_fail_chi  :", sp.N(Pe_fail_J) < sp.N(Pe_fail_chi))
print("Pe_suff_chi < Pe_fail_J  :", sp.N(Pe_suff_chi) < sp.N(Pe_fail_J))
```

(b) Append to `mathematica/moving_throat_pde_stage078_family1_branch_verdict_mathematica_audit.wl` after line 62 (i.e. before the existing `Print[""]` / `Print["Stage 078 Mathematica audit passed."]` / `Exit[0]` block, which currently spans lines 64-67 — preserve those at the end):

```mathematica

(* --- Branch verdict (chi vs Jensen-floor) ---------------------------- *)
expectTrue["Pe_suff^(J) < Pe_suff^(chi)", peSuffJ < peSuffChi];
expectTrue["Pe_fail^(J) < Pe_fail^(chi)", peFailJ < peFailChi];
expectTrue["Pe_suff^(chi) < Pe_fail^(J) (window overlap)", peSuffChi < peFailJ];
```

**Verification:**
After Codex applies, both scripts exit 0. The SymPy output text contains three new `Pe_suff_J < Pe_suff_chi  : True` / `Pe_fail_J < Pe_fail_chi  : True` / `Pe_suff_chi < Pe_fail_J  : True` lines, and the Mathematica output text contains three new `Pe_suff^(J) < Pe_suff^(chi) = True` / `Pe_fail^(J) < Pe_fail^(chi) = True` / `Pe_suff^(chi) < Pe_fail^(J) (window overlap) = True` lines with corresponding `PASS:` entries.

## Applied: F4

- files_changed:
  - `scripts/moving_throat_pde_stage078_family1_branch_verdict_sympy_audit.py`
  - `mathematica/moving_throat_pde_stage078_family1_branch_verdict_mathematica_audit.wl`
- summary: added chi-versus-Jensen branch-verdict assertions and output lines in both audit scripts
- deviation: none
