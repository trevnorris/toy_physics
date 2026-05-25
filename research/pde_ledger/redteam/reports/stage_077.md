---
unit_id: 077
batch: III.4
auditor_model: claude-opus-4-7[1m]
audit_date: 2026-05-22T00:00:00Z
verdict: findings
stop_cold: null
findings_count: 3
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: true
---

# Audit unit 077 red-team report

## Files reviewed

- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage077_family1_theta_extraction_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage077_family1_theta_extraction_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage077_family1_theta_extraction_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage077_family1_theta_extraction_mathematica_audit.txt`

## What the script claims to verify

Per its docstring, the script verifies four items for the Family-1 radial wall profile `rho^4 = 1 - alpha_r * S(xi)^2` with `S(xi) = (1+tanh(xi))/2`:
(1) the exact cut point `xi_* = atanh(2/sqrt(alpha_r) - 1)` (where the profile vanishes);
(2) the canonical support normalization `I_f = ∫ chi^2 d xi = 1/3` with `chi = dS/dxi = sech^2(xi)/2`;
(3) the numerical shell-weighted moments `<rho>_chi` and `<rho^2>_chi` on the `alpha_r = 10` branch (integrated from `-∞` to `xi_*` weighted by `chi^2`, normalized by `∫_{-∞}^{∞} chi^2`);
(4) the numerical effective wall-depth datum `Theta_w^(chi) = 25 <rho^2>_chi` and conservative Jensen floor `Theta_w^(J) = 25 <rho>_chi^2`. SymPy enforces only (2) symbolically and a Jensen-like ordering `Theta_chi >= Theta_J > 0`; Mathematica reproduces (2) symbolically and additionally enforces the numerical values of (1), (3), (4) by `expectApprox` against hardcoded 50-digit constants.

## Assertion inventory

| #  | Script        | Line   | Form                                                                 | Anchored to claim? |
|----|---------------|--------|----------------------------------------------------------------------|--------------------|
| A1 | sympy         | 39     | `expect_zero("I_f - 1/3", If - 1/3)`                                 | yes (claim 2)      |
| A2 | sympy         | 77-78  | `if not (Theta_chi >= Theta_J > 0): raise`                           | partial (Jensen ordering only; numerical values themselves not asserted) |
| A3 | mathematica   | 50     | `expectZero["I_f - 1/3", ifMom - 1/3]`                               | yes (claim 2)      |
| A4 | mathematica   | 89     | `expectApprox["xi_* numeric check", xiCut, "-0.38558...`50", 1e-30]` | partial — see F1   |
| A5 | mathematica   | 90     | `expectApprox["<rho>_chi numeric check", r1, "0.192619...`50", 1e-28]` | yes (cross-engine numerical) |
| A6 | mathematica   | 91     | `expectApprox["<rho^2>_chi numeric check", r2, "0.162745...`50", 1e-28]` | yes (cross-engine numerical) |
| A7 | mathematica   | 92     | `expectApprox["Theta_w^(chi) numeric check", thetaChi, "4.06863...`50", 1e-26]` | yes (cross-engine numerical) |
| A8 | mathematica   | 93     | `expectApprox["Theta_w^(J) numeric check", thetaJ, "0.927552...`50", 1e-27]` | yes (cross-engine numerical) |
| A9 | mathematica   | 94     | `expectTrue["Theta_w^(chi) >= Theta_w^(J) > 0", thetaChi >= thetaJ && thetaJ > 0]` | partial (Jensen ordering only) |

## Findings

### F1 — insufficient_verification

**Severity:** medium
**Files:**
- `scripts/moving_throat_pde_stage077_family1_theta_extraction_sympy_audit.py:41-43`
- `mathematica/moving_throat_pde_stage077_family1_theta_extraction_mathematica_audit.wl:46-52,89`

**What's wrong:**
The docstring's claim (1) is "Exact cut point `xi_*` for the Family-1 radial wall profile." The physical content of "cut point" is: `xi_*` is the value of `xi` where `1 - alpha_r * S(xi)^2 = 0`, i.e., where the fourth-root profile `rho = (1 - alpha_r * S^2)^{1/4}` first vanishes. Neither script ever asserts this.

SymPy (line 41) simply defines
```
xi_star = sp.simplify(sp.atanh(2 / sp.sqrt(alpha_r) - 1))
```
and prints it (lines 42-43). There is no `expect_zero` confirming that `1 - alpha_r * S(xi_star)^2 == 0` as a symbolic identity over `alpha_r`. The closed form is asserted only by author fiat.

Mathematica is no better: it defines `xiStar = FullSimplify[ArcTanh[2/Sqrt[alphaR] - 1], ...]` (line 46) and then `expectApprox` on line 89 only verifies that `N[xiStar /. alphaR -> 10, 50]` matches the literal 50-digit string `-0.38558106921542562...`. That literal is just `N[ArcTanh[2/Sqrt[10]-1], 50]` reproduced — see F2 — and exercises no physics.

**Why this matters:**
The unit's title and item (1) of the docstring claim an **exact** cut-point identity, but the framework currently treats `xi_*` as a definition, not a derivation. If `alpha_r` had been chosen differently (or if a future edit perturbed the formula by a factor of two), no assertion in either engine would fail until the downstream numerical Theta values shifted, and even then the failure would be diagnosed in the wrong place.

**Required change:**
In `scripts/moving_throat_pde_stage077_family1_theta_extraction_sympy_audit.py`, immediately after line 41 (where `xi_star` is defined), add a symbolic check that the radial profile vanishes there. Concretely, with `xi, alpha_r` already declared:
```
S_at_star = ((1 + sp.tanh(xi_star)) / 2)
rho_quartic_at_star = sp.simplify(1 - alpha_r * S_at_star**2)
expect_zero("1 - alpha_r * S(xi_*)**2", rho_quartic_at_star)
```
Because `tanh(atanh(2/sqrt(alpha_r) - 1)) = 2/sqrt(alpha_r) - 1`, `S(xi_*) = 1/sqrt(alpha_r)`, and `1 - alpha_r * (1/sqrt(alpha_r))^2 = 0` identically. SymPy will reduce it to 0 under the existing `positive=True, real=True` assumptions on `alpha_r`.

In `mathematica/moving_throat_pde_stage077_family1_theta_extraction_mathematica_audit.wl`, add the analogous symbolic check immediately after line 46 (where `xiStar` is defined), independent of any numeric evaluation:
```
sAtStar = (1 + Tanh[xiStar])/2;
rhoQuarticAtStar = FullSimplify[1 - alphaR*sAtStar^2, Assumptions -> $Assumptions];
expectZero["1 - alphaR*S[xi_*]^2", rhoQuarticAtStar];
```
This must reduce to `0` symbolically under `alphaR > 0`.

**Verification:**
After Codex applies, the SymPy output should contain a new line `1 - alpha_r * S(xi_*)**2 = 0` (followed by no AssertionError), and the Mathematica output should contain `1 - alphaR*S[xi_*]^2 = 0` and `PASS: 1 - alphaR*S[xi_*]^2`. Both scripts must still exit 0.

### F2 — tautological_check

**Severity:** low
**Files:**
- `mathematica/moving_throat_pde_stage077_family1_theta_extraction_mathematica_audit.wl:89`

**What's wrong:**
Line 89 is
```
expectApprox["xi_* numeric check", xiCut,
  ToExpression["-0.38558106921542562403635498846713378847348301441599`50"],
  10^-30];
```
But `xiCut` is defined on line 57 as `N[xiStar /. alphaR -> alphaNum, 50]`, i.e. `N[ArcTanh[2/Sqrt[10] - 1], 50]`. The hardcoded right-hand side is simply a literal copy of that same evaluation to 50 digits. The check therefore confirms a number against itself: there is no independently-derived target, no physics. If Mathematica's `ArcTanh` returned a wrong value, the literal in the script would presumably have been generated by the same (wrong) routine at script-authoring time, so the assertion would still pass.

Unlike the other `expectApprox` calls on lines 90-93 — which compare an independent numerical quadrature (NIntegrate) against a value computed by a different engine (SymPy/mpmath quad), and therefore constitute a real cross-engine check — line 89 has no second-engine content.

**Why this matters:**
The check looks substantive but adds no verification leverage. It also masks the real missing check (F1): if F1 is applied (the symbolic vanishing identity), the residual purpose of line 89 is sanity-checking Mathematica's `ArcTanh` numerical evaluator, which is not within the unit's scope.

**Required change:**
Delete line 89 of `mathematica/moving_throat_pde_stage077_family1_theta_extraction_mathematica_audit.wl`. Keep the print of `numeric xi_*` (line 82) so the value is still visible in the output.

The symbolic identity added by F1 (`1 - alphaR*S[xi_*]^2 == 0`) replaces this with a real check.

**Verification:**
After Codex applies, line 89 is removed; the Mathematica output will no longer contain the line `xi_* numeric check diff = 0...` or `PASS: xi_* numeric check`. The script must still exit 0 and all remaining `expectApprox` calls must still pass.

### F3 — insufficient_verification

**Severity:** low
**Files:**
- `scripts/moving_throat_pde_stage077_family1_theta_extraction_sympy_audit.py:55-79`

**What's wrong:**
The SymPy script computes `<rho>_chi`, `<rho^2>_chi`, `Theta_w^(chi)`, and `Theta_w^(J)` numerically at `alpha_r = 10` (lines 47-76) but the only assertion exercising these numbers is
```
if not (Theta_chi >= Theta_J > 0):
    raise AssertionError("Expected Theta_w^(chi) >= Theta_w^(J) > 0")
```
(lines 77-78). That is a Jensen-inequality ordering, which is `<f^2> >= <f>^2` and holds for any reasonable quadrature result with `f >= 0`; it does not verify the specific numerical values claimed in items (3) and (4) of the docstring.

The Mathematica side does enforce the specific values via `expectApprox` (lines 90-93), so the unit as a whole has cross-engine coverage; but the SymPy script's `expect_zero`-style assertion infrastructure (lines 24-28) is never used on the actual numerical results that the script's `print` statements parade as the bottom line. A reader of the SymPy script alone cannot tell whether the printed `<rho>_chi = 0.192619...` reflects a real check or just noise from `mp.quad`.

**Why this matters:**
The unit's numerical results are the data inputs to whatever downstream Theta_w consumer exists. Within the SymPy script, those values currently sit between two `print` calls and outside any assertion. If `mp.quad` regressed or the breakpoint `-4` on lines 62-64 was changed inadvertently, the SymPy script would still exit 0 silently as long as the Jensen ordering survived.

**Required change:**
In `scripts/moving_throat_pde_stage077_family1_theta_extraction_sympy_audit.py`, after line 76 (where `Theta_chi` and `Theta_J` are printed) and before the Jensen check on lines 77-78, add numerical-assertion lines using mpmath comparisons against the same 50-digit constants the Mathematica script already uses (these constants are the cross-engine ground truth, not invented by SymPy):
```
def expect_close(name, value, target, tol):
    diff = abs(value - target)
    print(f"{name} diff = {diff}")
    if diff > tol:
        raise AssertionError(f"{name} exceeds tol {tol}")

expect_close("<rho>_chi",
             R1, mp.mpf('0.19261900555649309777068139356018510792903510747507'),
             mp.mpf('1e-28'))
expect_close("<rho^2>_chi",
             R2, mp.mpf('0.16274529400326462037087418498629868328210821103971'),
             mp.mpf('1e-28'))
expect_close("Theta_w^(chi)",
             Theta_chi, mp.mpf('4.0686323500816155092718546246574670820527052759928'),
             mp.mpf('1e-26'))
expect_close("Theta_w^(J)",
             Theta_J, mp.mpf('0.92755203253930797183993260663904217023332624032789'),
             mp.mpf('1e-27'))
```
Keep the Jensen-ordering check on lines 77-78 intact.

**Verification:**
After Codex applies, the SymPy output should contain four new `<...> diff = <small number>` lines and the script must still exit 0. If any numeric quadrature drift introduced by an unrelated edit pushes a result above its tolerance, the script will now fail rather than silently pass.

## Independent-derivation check (Mathematica)

Both engines target the same physics (Family-1 wall profile `rho = (1 - alpha_r * S^2)^{1/4}`, `S = (1+tanh(xi))/2`, `chi = dS/dxi`) but they do not transliterate. Three observations:

1. The integration choreography differs. SymPy uses `mp.quad` with explicit breakpoints `[-mp.inf, -4, xi_cut]` to help adaptive quadrature (lines 62-64); Mathematica uses a single `NIntegrate` from `-Infinity` to `xiCut` with `WorkingPrecision -> 60` (lines 65-75) and uses `Quiet[..., NIntegrate::precw]` rather than introducing manual breakpoints.
2. The handling of the moment integrand differs. SymPy builds `rho_num` as the fourth root and then squares it inside the integrand for `<rho^2>` (lines 56-57, 64). Mathematica avoids the square-then-root by defining `rhoSqNum[x] := Sqrt[1 - alphaNum*sNum[x]^2]` directly (line 62) — a numerically distinct evaluation path.
3. Assertion strategies differ. SymPy enforces only the symbolic `I_f = 1/3` identity and a Jensen ordering; Mathematica enforces the symbolic identity *and* the numerical values to 28+ digits.

No transliteration concern.

## Engine cross-check

The numerical results match across engines to the precision claimed:

```
SymPy mp.quad (dps=50)                        Mathematica NIntegrate (WP=60)
<rho>_chi   = 0.19261900555649309777068139356  0.19261900555649309777068139356
<rho^2>_chi = 0.16274529400326462037087418498  0.16274529400326462037087418498
denom        = 0.33333333333333333333333333333  0.33333333333333333333333333333
Theta_chi   = 4.0686323500816155092718546246    4.0686323500816155092718546246
Theta_J     = 0.92755203253930797183993260663   0.92755203253930797183993260663
xi_*(10)    = -0.38558106921542562403635498846  -0.38558106921542562403635498846
```

Symbolic results also agree (`I_f = 1/3`, `xi_* = -atanh(1 - 2/sqrt(alpha_r))` aka `-(1/2)*Log[-1+Sqrt[alphaR]]`). `engines_agree: true`.

## Verdict justification

Engines agree numerically and symbolically. The `I_f = 1/3` identity is non-tautologically verified in both engines. The numerical moments and Theta values are protected on the Mathematica side by genuine cross-engine `expectApprox` checks but not on the SymPy side (F3). The "exact cut point" claim (docstring item 1) is asserted as a definition rather than derived from `1 - alpha_r * S^2 = 0` — easy to repair with a one-line symbolic identity check in each engine (F1), at which point the redundant `xi_*` numeric-evaluator self-check on Mathematica line 89 (F2) loses its purpose and should be removed. None of the issues propagate downstream; the verdict is `findings` with three low-to-medium items.

## Self-test notes

For F1's proposed check `expect_zero("1 - alpha_r * S(xi_*)**2", 1 - alpha_r * S(xi_star)^2)`: substituting `tanh(atanh(z)) = z`, `S(xi_*) = (1 + (2/sqrt(alpha_r) - 1))/2 = 1/sqrt(alpha_r)`, so `1 - alpha_r * (1/sqrt(alpha_r))^2 = 1 - 1 = 0`. SymPy's `simplify` with `alpha_r` positive collapses `tanh(atanh(...))` correctly; Mathematica's `FullSimplify` does the same. For F3's `expect_close` calls, the literal targets are the existing 50-digit Mathematica targets verbatim; mpmath's `abs(R1 - target)` with `R1` already at dps=50 yields a residual at or below 1e-30, well within the 1e-28 tolerance band I set (matching the Mathematica tolerances). No variable-independence trap (every `sp.diff` already involves a real dependence); no parity trap (no new integrals introduced).
