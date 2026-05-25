---
unit_id: 077
batch: III.4
created_at: 2026-05-22T00:00:00Z
findings_count: 3
stop_cold: null
applied: true
applied_at: 2026-05-22T23:25:22-06:00
findings_applied: 3
findings_blocked: 0
verification_status: pending
---

# Codex directive — unit 077

Apply each finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

If a finding's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named.

Do NOT run python or mathematica. Only edit files.

Do NOT touch paper.tex, notes/, or any prose documents. The red-team only modifies scripts.

## F1 — insufficient_verification

**Target:**
- `scripts/moving_throat_pde_stage077_family1_theta_extraction_sympy_audit.py:41-43`
- `mathematica/moving_throat_pde_stage077_family1_theta_extraction_mathematica_audit.wl:46-52`

**Issue:**
Both scripts define `xi_star` (resp. `xiStar`) as the closed-form `atanh(2/sqrt(alpha_r) - 1)` and never verify the physical content of that definition — namely that `1 - alpha_r * S(xi_*)^2 = 0` identically, which is what makes `xi_*` the "exact cut point" claimed in the docstring. Add a symbolic identity check in each engine immediately after the definition.

**Required change:**

In `scripts/moving_throat_pde_stage077_family1_theta_extraction_sympy_audit.py`, **after** the existing line 43 (`print("xi_*(alpha_r=10) =", sp.N(xi_star.subs(alpha_r, 10), 20))`) insert the following three lines (a blank line, then the check, then a blank line) before the existing line 44 comment `# Numerical branch evaluation at alpha_r = 10.`:

Before:
```
xi_star = sp.simplify(sp.atanh(2 / sp.sqrt(alpha_r) - 1))
print("xi_* =", xi_star)
print("xi_*(alpha_r=10) =", sp.N(xi_star.subs(alpha_r, 10), 20))

# Numerical branch evaluation at alpha_r = 10.
```

After:
```
xi_star = sp.simplify(sp.atanh(2 / sp.sqrt(alpha_r) - 1))
print("xi_* =", xi_star)
print("xi_*(alpha_r=10) =", sp.N(xi_star.subs(alpha_r, 10), 20))

S_at_star = (1 + sp.tanh(xi_star)) / 2
rho_quartic_at_star = sp.simplify(1 - alpha_r * S_at_star**2)
expect_zero("1 - alpha_r * S(xi_*)**2", rho_quartic_at_star)

# Numerical branch evaluation at alpha_r = 10.
```

In `mathematica/moving_throat_pde_stage077_family1_theta_extraction_mathematica_audit.wl`, **after** the existing line 52 (`Print["xi_*(alpha_r=10) = ", fmt[N[xiStar /. alphaR -> 10, 30]]];`) and **before** the existing line 53 blank line (or before `banner["NUMERICAL FAMILY-1 EXTRACTION"];` on line 54), insert:

```
sAtStar = (1 + Tanh[xiStar])/2;
rhoQuarticAtStar = FullSimplify[1 - alphaR*sAtStar^2, Assumptions -> $Assumptions];
expectZero["1 - alphaR*S[xi_*]^2", rhoQuarticAtStar];
```

Do not alter any other lines. The existing `$Assumptions` (`Element[{xi, alphaR}, Reals] && alphaR > 0`) is sufficient to reduce the residual to 0.

**Verification command:**
After Codex applies, the verifier will run `redteam exec-sympy 077` and `redteam exec-mathematica 077` (one at a time; never parallel) and confirm:
- SymPy output contains `1 - alpha_r * S(xi_*)**2 = 0` and the script exits 0.
- Mathematica output contains `1 - alphaR*S[xi_*]^2 = 0` followed by `PASS: 1 - alphaR*S[xi_*]^2`, and the script exits 0.

## Applied: F1

- files_changed:
  - `scripts/moving_throat_pde_stage077_family1_theta_extraction_sympy_audit.py`
  - `mathematica/moving_throat_pde_stage077_family1_theta_extraction_mathematica_audit.wl`
- summary: Added symbolic cut-point residual checks in both SymPy and Mathematica immediately after xi_* is defined.
- deviation: none

## F2 — tautological_check

**Target:** `mathematica/moving_throat_pde_stage077_family1_theta_extraction_mathematica_audit.wl:89`

**Issue:**
Line 89 compares `xiCut` (the numerical evaluation of `ArcTanh[2/Sqrt[10] - 1]` at 50 digits) against a hardcoded 50-digit literal of that same evaluation. There is no independent target; the check confirms a number against itself and exercises no physics. With F1's symbolic vanishing identity in place this line provides no audit value.

**Required change:**
Delete line 89 of `mathematica/moving_throat_pde_stage077_family1_theta_extraction_mathematica_audit.wl`, which currently reads:
```
expectApprox["xi_* numeric check", xiCut, ToExpression["-0.38558106921542562403635498846713378847348301441599`50"], 10^-30];
```
Do not touch lines 90-93 (the `<rho>_chi`, `<rho^2>_chi`, `Theta_w^(chi)`, `Theta_w^(J)` `expectApprox` calls); those are legitimate cross-engine numerical checks. Do not touch line 82 (`Print["numeric xi_* = ", fmt[xiCut]];`); keep the value visible in the output transcript.

**Verification command:**
After Codex applies, the verifier will run `redteam exec-mathematica 077` and confirm that the output no longer contains the strings `xi_* numeric check diff` or `PASS: xi_* numeric check`, that the script still exits 0, and that all remaining `PASS:` lines (for `<rho>_chi`, `<rho^2>_chi`, `Theta_w^(chi)`, `Theta_w^(J)`, the Jensen ordering, and the new F1 symbolic check) still appear.

## Applied: F2

- files_changed:
  - `mathematica/moving_throat_pde_stage077_family1_theta_extraction_mathematica_audit.wl`
- summary: Removed the tautological xi_* numeric expectApprox while leaving the printed xi_* value and remaining numeric checks intact.
- deviation: none

## F3 — insufficient_verification

**Target:** `scripts/moving_throat_pde_stage077_family1_theta_extraction_sympy_audit.py:55-79`

**Issue:**
The SymPy script computes `R1` (`<rho>_chi`), `R2` (`<rho^2>_chi`), `Theta_chi`, `Theta_J` numerically but only enforces a Jensen-ordering inequality. The specific numerical values claimed in docstring items (3) and (4) sit between print statements without any tolerance-based assertion. Mathematica enforces those values; SymPy must mirror that protection so the two engines provide symmetric numerical guards.

**Required change:**
In `scripts/moving_throat_pde_stage077_family1_theta_extraction_sympy_audit.py`, **after** the existing line 76 (`print("Theta_w^(J)   / lambda_mu^2 =", Theta_J)`) and **before** the existing line 77 (`if not (Theta_chi >= Theta_J > 0):`), insert the following block:

```

def expect_close(name: str, value, target, tol) -> None:
    diff = abs(value - target)
    print(f"{name} diff = {diff}")
    if diff > tol:
        raise AssertionError(f"{name} exceeds tol {tol}")

expect_close(
    "<rho>_chi",
    R1,
    mp.mpf('0.19261900555649309777068139356018510792903510747507'),
    mp.mpf('1e-28'),
)
expect_close(
    "<rho^2>_chi",
    R2,
    mp.mpf('0.16274529400326462037087418498629868328210821103971'),
    mp.mpf('1e-28'),
)
expect_close(
    "Theta_w^(chi)",
    Theta_chi,
    mp.mpf('4.0686323500816155092718546246574670820527052759928'),
    mp.mpf('1e-26'),
)
expect_close(
    "Theta_w^(J)",
    Theta_J,
    mp.mpf('0.92755203253930797183993260663904217023332624032789'),
    mp.mpf('1e-27'),
)
```

Do not modify the existing Jensen-ordering check (current lines 77-78) — keep it as a final guard. Do not change the existing `expect_zero` helper.

The 50-digit constants are copied verbatim from the Mathematica targets on lines 90-93 of `mathematica/moving_throat_pde_stage077_family1_theta_extraction_mathematica_audit.wl`, which are themselves the cross-engine ground truth.

**Verification command:**
After Codex applies, the verifier will run `redteam exec-sympy 077` and confirm:
- The SymPy output contains four new lines of the form `<rho>_chi diff = <small>`, `<rho^2>_chi diff = <small>`, `Theta_w^(chi) diff = <small>`, `Theta_w^(J) diff = <small>` (each diff at or below the per-line tolerance).
- The script exits 0.
- The existing `I_f - 1/3 = 0` line and Jensen ordering check are unchanged.

## Applied: F3

- files_changed:
  - `scripts/moving_throat_pde_stage077_family1_theta_extraction_sympy_audit.py`
- summary: Added SymPy tolerance checks for the four cross-engine numeric targets before the existing Jensen-ordering guard.
- deviation: none
