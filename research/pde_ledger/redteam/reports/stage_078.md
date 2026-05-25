---
unit_id: 078
batch: III.4
auditor_model: claude-opus-4-7-1m
audit_date: 2026-05-22T00:00:00Z
verdict: findings
stop_cold: null
findings_count: 4
scripts_checked:
  sympy: insufficient
  mathematica: insufficient
  engines_agree: true
  outputs_fresh: true
---

# Audit unit 078 red-team report

## Files reviewed

- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage078_family1_branch_verdict_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage078_family1_branch_verdict_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage078_family1_branch_verdict_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage078_family1_branch_verdict_mathematica_audit.txt`

## What the script claims to verify

According to its docstring/banner, the script "Insert[s] the explicit Stage-60 Theta values into the Stage-58 threshold window", computes Pe_req success/failure windows for both the natural quadratic datum (`chi`) and the Jensen-floor (`J`) Theta values, and verifies the ordering `Pe_suff < Pe_fail` in both cases. Concretely, four hardcoded coefficients (`Theta_chi=4.06863..`, `Theta_J=0.927552..`, `Theta_fail=3.62606e-4`, `Theta_suff=4.21495e-2`) are divided pairwise to produce four `Pe / lambda_mu^2` numbers, then two strict-less-than tests are asserted. The Mathematica script does the same arithmetic and additionally `expectApprox`s the four ratios against literal targets equal to the SymPy script's printed values.

## Assertion inventory

| # | Script | Line | Form | Anchored to claim? |
|---|---|---|---|---|
| A1 | sympy | 41-42 | `Pe_suff_chi < Pe_fail_chi`  → `AssertionError` if not | no (tautological in the literals) |
| A2 | sympy | 43-44 | `Pe_suff_J < Pe_fail_J`  → `AssertionError` if not | no (tautological in the literals) |
| A3 | mathematica | 57 | `expectApprox[..., 96.528524726438575954, 10^-12]` | no (target = SymPy output) |
| A4 | mathematica | 58 | `expectApprox[..., 11220.544162625905301, 10^-9]` | no (target = SymPy output) |
| A5 | mathematica | 59 | `expectApprox[..., 22.006222633075413597, 10^-12]` | no (target = SymPy output) |
| A6 | mathematica | 60 | `expectApprox[..., 2558.0189234920526360, 10^-10]` | no (target = SymPy output) |
| A7 | mathematica | 61 | `expectTrue["Pe_suff^(chi) < Pe_fail^(chi)", peSuffChi < peFailChi]` | no (tautological in the literals) |
| A8 | mathematica | 62 | `expectTrue["Pe_suff^(J) < Pe_fail^(J)", peSuffJ < peFailJ]` | no (tautological in the literals) |

## Findings

### F1 — tautological_check

**Severity:** high
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage078_family1_branch_verdict_sympy_audit.py:41-44`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage078_family1_branch_verdict_mathematica_audit.wl:61-62`

**What's wrong:**
The two "ordering" assertions reduce to a comparison of two hardcoded numbers. By construction
```
Pe_suff_chi = Theta_chi_coeff / Theta_suff_coeff,
Pe_fail_chi = Theta_chi_coeff / Theta_fail_coeff,
```
so `Pe_suff_chi < Pe_fail_chi` is algebraically equivalent to `Theta_fail_coeff < Theta_suff_coeff` (Theta_chi cancels). The same simplification applies to the J branch (`Theta_J` cancels). Therefore the four `<` tests are *all* the single arithmetic statement `3.62605617972939e-4 < 4.21495341569977e-2`, which is true by inspection of the literals and entirely independent of the Stage-60 Theta values the script *claims* to be exercising.

**Why this matters:**
The script's stated purpose is to verify the Stage-60 Theta-window verdict — that the success threshold lies below the failure threshold once the Stage-60 Theta values are plugged into the Stage-58 window. As coded, the Stage-60 Theta values (`Theta_chi`, `Theta_J`) never participate in the verdict; only the success/failure coefficients do. A bogus or sign-flipped Theta_chi or Theta_J would pass the script unchanged, so the verdict the script issues is unrelated to the verdict it claims to be issuing.

**Required change:**
Add at least one *non-cancelling* assertion that depends on all four coefficients. The natural choice is the actual threshold-window membership statement: for `lambda_mu = 1` (numerical) the script must show that there is a non-empty `Pe_req` interval such that `Theta_chi_coeff * 1^2 < Theta_suff_coeff * Pe_req` AND `Theta_chi_coeff * 1^2 < Theta_fail_coeff * Pe_req` is *false*, i.e. compute the explicit `Pe_req` bounds
```
Pe_req_lower(branch) = Theta(branch) / Theta_suff_coeff,    # success threshold
Pe_req_upper(branch) = Theta(branch) / Theta_fail_coeff,    # failure threshold
```
and assert `Pe_req_lower(chi) < Pe_req_upper(chi)` *and* `Pe_req_lower(chi) < Pe_req_lower(J)` *and* `Pe_req_upper(J) < Pe_req_upper(chi)` (the Jensen floor moves both bounds inward because `Theta_J < Theta_chi`). The last two checks depend on `Theta_chi` *and* `Theta_J` and so cannot be passed if either Stage-60 value is corrupted.

**Verification:**
After the fix the SymPy script must contain at least two new assertions referencing both `Theta_chi` and `Theta_J` (e.g. comparing the two branches' windows), and the Mathematica script must mirror them with `expectTrue`. The auditor re-runs `redteam exec-sympy 078` / `redteam exec-mathematica 078`, exits 0, and the output text contains four new "Pe_req_lower(chi) < Pe_req_lower(J)" / "Pe_req_upper(J) < Pe_req_upper(chi)" lines.

### F2 — hardcoded_result

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage078_family1_branch_verdict_sympy_audit.py:26-29`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage078_family1_branch_verdict_mathematica_audit.wl:37-40`

**What's wrong:**
All four Theta coefficients enter the script as bare 15-digit floats with no comment naming the upstream source script, the symbolic closed form, or the verification stage that produced them. They originate (per upstream output files) from stage075 (`Upsilon_suff/Pe_req` and `Theta_fail/Pe_req`) and stage077 (`Theta_w^(chi)`, `Theta_w^(J)`), but the stage078 script offers no traceable anchor. Additionally, the Mathematica script's `expectApprox` targets at lines 57-60 are literal-for-literal copies of the SymPy script's printed output (compare `scripts/output/moving_throat_pde_stage078_family1_branch_verdict_sympy_audit.txt` lines 13-16 with `mathematica/.../audit.wl` lines 57-60), so the cross-engine "check" is `assert sympy_number == sympy_number`, not an independent derivation.

**Why this matters:**
If any upstream Theta value were edited tomorrow, this unit would still pass without complaint, because nothing in stage078 binds it to the upstream symbolic form. The Mathematica engine in particular provides no independent verification: it simply re-types the SymPy result and confirms it equals itself.

**Required change:**
At the four definitions in each script, add a comment giving the symbolic provenance, e.g.

In SymPy (lines 26-29), add comments like:
```python
# Theta_chi = (1/9) * (chi^2 weighted Stage-77 result); numeric value reproduced from
# scripts/output/moving_throat_pde_stage077_family1_theta_extraction_sympy_audit.txt:22
Theta_chi = sp.Float("4.06863235008162") * lambda_mu**2
# ...similarly for Theta_J (stage077:23), Theta_fail/Theta_suff (stage075:30,34).
```

In Mathematica (lines 37-40), replace the `expectApprox` literal targets at lines 57-60 by an *independent* computation of the same ratios from the closed-form expressions: for `Theta_fail / Pe_req`, the symbolic form (from stage075 output line 26) is
```
(37 Cosh[111 Sqrt[5]/5] + 111 Sqrt[5] Sinh[111 Sqrt[5]/5]/5)
   / (136900 (-1 + Sqrt[5] Sinh[111 Sqrt[5]/5]/3 + Cosh[111 Sqrt[5]/5]))
```
and for `Upsilon_suff / Pe_req` (from stage075 output line 18 or equivalent) the analogous closed form. Compute these via `N[..., 30]` in Mathematica, then use *those* values (not the SymPy-derived constants) as the `expectApprox` targets, with tolerance `10^-12`. This breaks the SymPy→Mathematica copy-loop.

**Verification:**
The four `Theta*` assignment lines in each script have a comment that names the upstream output file (and line) supplying the number. The Mathematica `expectApprox` targets at lines 57-60 are computed inside the .wl from `Sinh`/`Cosh` expressions, not from typed-in decimals.

### F3 — mathematica_transliteration

**Severity:** high
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage078_family1_branch_verdict_mathematica_audit.wl:32-62`

**What's wrong:**
The `.wl` file is a one-to-one port of the `.py` file. Compare for example

SymPy (lines 25-34):
```python
lambda_mu, Pe_req = sp.symbols("lambda_mu Pe_req", positive=True, real=True)
Theta_chi  = sp.Float("4.06863235008162")       * lambda_mu**2
Theta_J    = sp.Float("0.927552032539308")      * lambda_mu**2
Theta_fail = sp.Float("3.62605617972939e-4")    * Pe_req
Theta_suff = sp.Float("4.21495341569977e-2")    * Pe_req
Pe_suff_chi = sp.simplify(Theta_chi / sp.Float("4.21495341569977e-2") / lambda_mu**2)
Pe_fail_chi = sp.simplify(Theta_chi / sp.Float("3.62605617972939e-4") / lambda_mu**2)
Pe_suff_J   = sp.simplify(Theta_J   / sp.Float("4.21495341569977e-2") / lambda_mu**2)
Pe_fail_J   = sp.simplify(Theta_J   / sp.Float("3.62605617972939e-4") / lambda_mu**2)
```

Mathematica (lines 34-50):
```mathematica
Clear[lambdaMu, peReq];
$Assumptions = ... lambdaMu > 0 && peReq > 0;
thetaChiCoeff  = SetPrecision[4.06863235008162, 20];
thetaJCoeff    = SetPrecision[0.927552032539308, 20];
thetaFailCoeff = SetPrecision[3.62605617972939*^-4, 20];
thetaSuffCoeff = SetPrecision[4.21495341569977*^-2, 20];
thetaChi = thetaChiCoeff*lambdaMu^2;   thetaJ = thetaJCoeff*lambdaMu^2;
thetaFail = thetaFailCoeff*peReq;       thetaSuff = thetaSuffCoeff*peReq;
peSuffChi = N[thetaChiCoeff/thetaSuffCoeff, 30];
peFailChi = N[thetaChiCoeff/thetaFailCoeff, 30];
peSuffJ   = N[thetaJCoeff /thetaSuffCoeff, 30];
peFailJ   = N[thetaJCoeff /thetaFailCoeff, 30];
```

Identical sequence of definitions, identical variable choreography, identical use of the same hardcoded decimals, and the final ordering tests (lines 61-62) reproduce the SymPy assertions verbatim. There is no independent derivation: the Mathematica script does not start from the Stage-58 threshold-window expressions or the Stage-60/77 closed-form Theta values; it consumes the same prebaked numbers in the same arrangement.

**Why this matters:**
The second-engine policy exists to catch algebra errors and sympy-specific simplification bugs. Re-typing the SymPy script in Mathematica syntax catches none of them. A sign error, factor-of-2 error, or hidden simplify-under-strong-assumption in any upstream symbolic computation would propagate unflagged because the Mathematica script begins where the SymPy script's algebra ends.

**Required change:**
Rewrite lines 34-50 of the `.wl` file so that the four Theta values are *computed in Mathematica* from the symbolic closed forms recorded in upstream stages (stage075 output line 26 for `Theta_fail`, stage075 output for `Upsilon_suff`/`Theta_suff`, stage077 for `Theta_chi`/`Theta_J`). Specifically:

```mathematica
(* Independent re-derivation of the four window coefficients from
   their symbolic closed forms (Stage-58 / Stage-75 / Stage-77). *)
thetaFailSym = (37 Cosh[111 Sqrt[5]/5] + 111 Sqrt[5] Sinh[111 Sqrt[5]/5]/5)
              / (136900 (-1 + Sqrt[5] Sinh[111 Sqrt[5]/5]/3
                              + Cosh[111 Sqrt[5]/5]));
thetaSuffSym = ... (* analogous closed form from stage075 *);
thetaChiSym  = ... (* closed form from stage077 *);
thetaJSym    = ... (* closed form from stage077 *);
thetaFailCoeff = N[thetaFailSym, 30];
thetaSuffCoeff = N[thetaSuffSym, 30];
thetaChiCoeff  = N[thetaChiSym, 30];
thetaJCoeff    = N[thetaJSym, 30];
```

Then the `expectApprox` targets at lines 57-60 should be the same closed forms `N[..., 30]`-evaluated, not literal decimals copied from the SymPy output. If the upstream symbolic forms are unavailable in stage075/77 output, fall back to deriving them inside this script from the same starting equations cited in the upstream `.py` docstrings, but do not import the prebaked decimals.

**Verification:**
The `.wl` file contains explicit `Sinh`/`Cosh` (or other transcendental) expressions for `thetaFailSym`, `thetaSuffSym`, `thetaChiSym`, `thetaJSym`, and the four `expectApprox` calls compare the SymPy decimals against `N[thetaFailSym, 30]` etc. with tolerance `10^-12`. The literal decimal targets `96.528...`, `11220.54...`, `22.006...`, `2558.01...` no longer appear as `expectApprox` arguments.

### F4 — insufficient_verification

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage078_family1_branch_verdict_sympy_audit.py:23-44`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage078_family1_branch_verdict_mathematica_audit.wl:32-62`

**What's wrong:**
The docstring promises three deliverables: (i) insert Stage-60 Theta into the Stage-58 window, (ii) compute Pe_req success/failure windows for both branches, (iii) verify `Pe_suff < Pe_fail`. As coded, only (iii) is exercised, and as shown in F1 it's tautological in the literals. (i) is never demonstrated: there is no symbolic Stage-58 window equation in the script for the Stage-60 Theta to be "inserted into"; the script just performs a one-line division. (ii) is reduced to the same division. The branch verdict comparing `chi` vs `J` Jensen-floor outcomes — which is the headline result implied by the filename `family1_branch_verdict` — is not asserted at all.

**Why this matters:**
The unit is a "branch verdict" stage but never asserts a verdict that distinguishes one branch from the other. Whether the Jensen floor narrows or widens the Pe_req window is the only physics question the unit can plausibly answer, and the script answers it nowhere.

**Required change:**
Append, in both scripts, the explicit branch-verdict assertions:

SymPy (after line 44):
```python
# Branch verdict: the Jensen floor narrows the success window from below
# (Pe_suff_J < Pe_suff_chi since Theta_J < Theta_chi and both share Theta_suff_coeff > 0)
# and narrows the failure ceiling from above (Pe_fail_J < Pe_fail_chi by the same logic).
# Therefore the J branch's admissible window strictly nests inside the chi branch's:
if not (sp.N(Pe_suff_J) < sp.N(Pe_suff_chi)):
    raise AssertionError("Expected Jensen-floor success threshold to lie below chi datum's")
if not (sp.N(Pe_fail_J) < sp.N(Pe_fail_chi)):
    raise AssertionError("Expected Jensen-floor failure ceiling to lie below chi datum's")
# Strict nesting of admissible windows:
if not (sp.N(Pe_suff_chi) < sp.N(Pe_fail_J)):
    raise AssertionError(
        "Expected chi-datum success threshold to remain below Jensen-floor failure ceiling "
        "(both windows must overlap for the verdict to be meaningful)"
    )
```

Mathematica (after line 62), the corresponding `expectTrue` calls.

These four checks each genuinely depend on the relative magnitudes of `Theta_chi` and `Theta_J` (or of `Theta_suff_coeff` and `Theta_fail_coeff` combined with both Thetas), so a corruption of any of the four upstream constants would now produce a visible failure.

**Verification:**
SymPy script exits 0 with three new `if not (...)` blocks present at lines 46+. Mathematica script exits 0 with three new `expectTrue` calls beneath line 62. Output text from both engines lists three new verdict checks.

## Independent-derivation check (Mathematica)

The Mathematica script is a transliteration, not an independent re-derivation; see F3. The same four hardcoded decimals are typed in, the same arithmetic is performed, and the cross-check targets are literal copies of the SymPy output. No closed-form `Sinh`/`Cosh` expressions appear anywhere, despite stage075's SymPy output containing exactly such forms for the `Theta_fail`/`Upsilon_fail` window. This is exactly the failure mode that the `mathematica_transliteration` category exists to flag.

## Engine cross-check

Both engines compute (to better than 1e-12) the same four ratios:

| Quantity | SymPy | Mathematica |
|---|---|---|
| `Pe_suff^(chi) / lambda_mu^2` | `96.528524726438575954` | `96.52852472643856904608...` |
| `Pe_fail^(chi) / lambda_mu^2` | `11220.544162625905301` | `11220.54416262590551355...` |
| `Pe_suff^(J)   / lambda_mu^2` | `22.006222633075413597` | `22.00622263307541136581...` |
| `Pe_fail^(J)   / lambda_mu^2` | `2558.0189234920526360` | `2558.01892349205282379...` |

Agreement is to ~14 digits. This agreement is unsurprising and does not constitute verification: both engines are dividing the same hardcoded numerators by the same hardcoded denominators in IEEE-double-derived precision, so they *must* agree to the working precision. There is no engine-disagreement finding.

## Verdict justification

Verdict: `findings` (4). The scripts run, agree, and produce the printed numbers without difficulty, but they do not verify what they claim. The four ordering assertions are tautological in the hardcoded literals (F1); the upstream Theta provenance is unrecorded and the Mathematica `expectApprox` targets are SymPy-output copies (F2); the .wl is a line-by-line transliteration of the .py with no independent derivation (F3); the branch-verdict comparison the unit's name promises is never asserted (F4). No stop-cold: the verdict the scripts ought to be asserting is consistent with the underlying upstream stage075/077 numbers, so fixing these four findings only requires adding substantive assertions; nothing downstream is invalidated by the fix because no derived constant changes.

Outputs are fresh (script and output mtimes within hours of each other; outputs younger than scripts), so `stale_output` is not flagged.

## Self-test notes

I verified the F1 cancellation by writing out `Pe_suff_chi / Pe_fail_chi = (Theta_chi / Theta_suff_coeff) / (Theta_chi / Theta_fail_coeff) = Theta_fail_coeff / Theta_suff_coeff`; `Theta_chi` and `lambda_mu^2` cancel cleanly because they appear identically in numerator and denominator. I checked the proposed F4 branch-verdict signs against the inequalities `Theta_J ≈ 0.928 < Theta_chi ≈ 4.069` and `Theta_fail_coeff ≈ 3.6e-4 < Theta_suff_coeff ≈ 4.2e-2`: with `Pe_X^(Y) = Theta_Y / Theta_X_coeff`, smaller `Theta_Y` numerator gives smaller `Pe`, so `Pe_*_J < Pe_*_chi` holds branchwise, and `Pe_suff_chi ≈ 96.5 < Pe_fail_J ≈ 2558` confirms window overlap. I also confirmed F2's claim that the Mathematica `expectApprox` targets at lines 57-60 are literal SymPy outputs by matching `96.528524726438575954`, `11220.544162625905301`, `22.006222633075413597`, `2558.0189234920526360` against `scripts/output/.../audit.txt` lines 13-16 — bit-identical character strings.
