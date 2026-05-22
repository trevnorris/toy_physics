---
unit_id: 031
batch: II.1
auditor_model: claude-opus-4-7[1m]
audit_date: 2026-05-21T00:00:00Z
verdict: findings
stop_cold: null
findings_count: 4
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: true
---

# Audit unit 031 red-team report

## Files reviewed

- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage031_selected_branch_reachability_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage031_selected_branch_reachability_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage031_selected_branch_reachability_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage031_selected_branch_reachability_mathematica_audit.txt`

## What the script claims to verify

The unit treats a 2x2 reduction whose two eigenvalues are `lambda_+/-(alpha) = (2A + DK - alpha*(x0+x1) +/- R)/2` with `R = sqrt((DK + alpha*(x0-x1))^2 + 4*alpha^2*x0*x1)`, and a selected `-` branch overlap `s_-(alpha) = -d lambda_-/d alpha` and prefactor `P0_-(alpha) = beta0 * s_-/lambda_-`. The scripts assert: (i) the exact derivative `ds_-/d alpha = 2*DK^2*x0*x1/R^3`; (ii) a quotient/Hellmann-Feynman identity for `dP0_-/d alpha`; (iii) the initial values `lambda_-(0)=A`, `s_-(0)=x0`, `P0_-(0)=beta0*x0/A`; (iv) the softening threshold `alpha_crit = A*(A+DK)/T0` with `T0 = (A+DK)*x0 + A*x1`, including the determinant factorization `lambda_- lambda_+ = A(A+DK) - alpha*T0` and a perfect-square identity for the threshold radicand; (v) the stable-side factorization `P0_- = beta0*s_-*lambda_+/(T0*(alpha_crit - alpha))`.

## Assertion inventory

| # | Script | Line | Form | Anchored to claim? |
|---|---|---|---|---|
| A1 | sympy | 50 | `simplify(ds_exact - ds_expected) == 0` | yes |
| A2 | sympy | 69 | `simplify(dP_generic - dP_expected) == 0` | no (generic quotient identity, not physical `D[P0_sel]`) |
| A3 | sympy | 76 | `simplify(lam_minus(0) - A) == 0` | yes |
| A4 | sympy | 77 | `simplify(s_minus(0) - x0) == 0` | yes |
| A5 | sympy | 78 | `simplify(P0_sel(0) - beta0*x0/A) == 0` | yes |
| A6 | sympy | 84 | `expand(lam_minus*lam_plus - (A*(A+DK) - alpha*T0)) == 0` | yes |
| A7 | sympy | 85 | `simplify((A*(A+DK) - alpha*T0).subs(alpha, alpha_crit)) == 0` | no (tautological by definition of `alpha_crit`) |
| A8 | sympy | 89-92 | `radcrit - (A^2 x1 + (A+DK)^2 x0)^2 == 0` | partial (hand-typed `radcrit`, no link to `T0^2 R(alpha_crit)^2`) |
| A9 | sympy | 101-104 | `expand(lam_minus*lam_plus - T0*(alpha_crit - alpha)) == 0` | yes |
| A10 | sympy | 108 | `simplify(P0_sel - P0_factored) == 0` | yes (follows from A6/A9) |
| M1 | mma | 42 | `FullSimplify[dsExact - dsExpected] == 0` | yes |
| M2 | mma | 49 | `FullSimplify[D[p0Sel,alpha] - dPExpected] == 0` | yes (direct physical derivative) |
| M3-M5 | mma | 53-55 | initial values | yes |
| M6 | mma | 60 | det factorization | yes |
| M7 | mma | 61 | `(a*(a+dK) - alpha*t0) /. alpha -> alphaCrit == 0` | no (same tautology as A7) |
| M8 | mma | 66 | `radCrit - ((a+dK)^2 x0 + a^2 x1)^2 == 0` | partial (same hand-typed polynomial as sympy) |
| M9 | mma | 71 | stable-side `lambda_- lambda_+ - t0(alphaCrit - alpha)` | yes |
| M10 | mma | 73 | `P0_- factorization` | yes |

## Findings

### F1 — insufficient_verification

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage031_selected_branch_reachability_sympy_audit.py:54-72`

**What's wrong:**
PART II of the SymPy script never directly verifies that the physical prefactor derivative `d/d alpha[beta0 * s_minus / lam_minus]` equals `beta0 * (ds_expected * lam_minus + s_minus**2) / lam_minus**2`. Instead, it constructs a generic identity using abstract `sp.Function("S")(alpha)` and `sp.Function("L")(alpha)` placeholders:

```
dP_generic = sp.diff(beta0 * sp.Function("S")(alpha) / sp.Function("L")(alpha), alpha)
dP_generic = dP_generic.subs({
    sp.diff(sp.Function("S")(alpha), alpha): DSsym,
    sp.Function("S")(alpha): Ssym,
    sp.diff(sp.Function("L")(alpha), alpha): -Ssym,
    sp.Function("L")(alpha): Lsym,
})
dP_expected = sp.simplify(beta0 * (DSsym * Lsym + Ssym**2) / Lsym**2)
expect_zero("generic quotient/HF identity", sp.simplify(dP_generic - dP_expected))
```

The check `dP_generic - dP_expected == 0` is a pure calculus identity (quotient rule with `L' = -S` substituted). It never invokes `s_minus`, `lam_minus`, or `ds_expected` from the script's actual eigenvalue setup. As a result it cannot detect, for example, a sign error in `s_minus`'s definition or a mis-derivation in `ds_expected`. The Mathematica script (line 47-49) does this correctly — `dPExact = D[p0Sel, alpha]` is the physical derivative and is compared to `dPExpected = beta0*(dsExpected*lamMinus + sMinus^2)/lamMinus^2`. The asymmetry is the smoking gun: the assertion the script needs is the one the Mathematica file already performs.

**Why this matters:**
The script's docstring lists "Exact derivative formula for the selected prefactor P0_-(alpha)" as a claim. PART II as written does not verify that formula for the actual `P0_-(alpha)`; it only verifies the calculus identity that would relate the formula to its derivative IF the substitutions were the true ones. Under the chain (A1 + HF-by-construction + calculus identity) the formula does hold, but the assertion is structurally too weak to catch any error in (a) the explicit form of `ds_expected`, (b) the sign convention on `s_minus`, or (c) an algebraic miswrite of `dP_expected`. A direct `sp.simplify(sp.diff(P0_sel, alpha) - beta0*(ds_expected*lam_minus + s_minus**2)/lam_minus**2)` is the substantive check.

**Required change:**
Add a direct assertion in `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage031_selected_branch_reachability_sympy_audit.py` immediately after the existing PART II block (insert before line 71, which prints `dP0_-/dalpha`). The new check must compute `sp.diff(P0_sel, alpha)` symbolically and compare to `beta0*(ds_expected*lam_minus + s_minus**2)/lam_minus**2`. See directive for exact insertion.

**Verification:**
After the patch, the saved sympy output must contain a new line `dP0_-/dalpha direct identity = 0` (between the existing `generic quotient/HF identity = 0` line and the `dP0_-/dalpha =` pretty-print). Script must continue to exit 0.

### F2 — tautological_check

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage031_selected_branch_reachability_sympy_audit.py:85`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage031_selected_branch_reachability_mathematica_audit.wl:61`

**What's wrong:**
The "alpha_crit condition" assertion is `(A*(A+DK) - alpha*T0).subs(alpha, alpha_crit) == 0` where `alpha_crit = A*(A+DK)/T0` was just defined one line earlier. Substituting `alpha_crit` into `A*(A+DK) - alpha*T0` gives `A*(A+DK) - (A*(A+DK)/T0)*T0 = 0` purely by algebra, regardless of physics. The Mathematica script has the identical tautology at line 61. Neither check exercises a physical claim — they only confirm `alpha_crit` is, by definition, the root of the expression it was defined from.

**Why this matters:**
The physical content of `alpha_crit` is that `lambda_-(alpha_crit) = 0` (the stable eigenvalue collapses at the threshold). The current assertion does not verify this; it only re-states a definitional rearrangement. Without the substantive check, a sign error in `T0` (e.g., writing `T0 = (A+DK)*x0 - A*x1` by typo) would propagate into a wrong `alpha_crit` and the existing assertion would still pass. PART V's `lambda_- * lambda_+ - T0*(alpha_crit - alpha) == 0` would then also pass (since both sides are equally wrong), so this catch is genuinely missed by the existing suite.

**Required change:**
Replace the tautological assertion with a substantive check that `lambda_-(alpha_crit) = 0`. See directive for exact replacement text and the corresponding fix in the Mathematica file.

**Verification:**
After the patch, the saved sympy output must show `lam_-(alpha_crit) = 0` in place of `alpha_crit condition = 0`, and the saved mathematica output must show `PASS: lam_-(alpha_crit)` in place of `PASS: alpha_crit condition`. Both scripts must continue to exit 0.

### F3 — hardcoded_result

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage031_selected_branch_reachability_sympy_audit.py:88-92`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage031_selected_branch_reachability_mathematica_audit.wl:64-66`

**What's wrong:**
The 9-term polynomial `radcrit = A**4*x0**2 + 2*A**4*x0*x1 + A**4*x1**2 + 4*A**3*DK*x0**2 + ...` is hand-typed in both scripts with no in-script derivation linking it to the physics. The assertion then verifies `radcrit == (A**2 x1 + (A+DK)**2 x0)**2` — a polynomial identity between two hand-supplied expressions, both pulled from the same external source. Neither expression is anchored to the radicand `R**2` evaluated at `alpha = alpha_crit`. A typo in `radcrit` that happened to preserve the perfect-square structure (or a typo in both hand-typed forms in the same direction) would slip through unnoticed. Both engines reuse the *same* hand-typed polynomial verbatim, so a transcription error in one would also corrupt the other.

**Why this matters:**
The "threshold radical square identity" is asserted to express `sqrt(R**2(alpha_crit))` in closed form (which is what makes `lambda_+^(crit)` printable as a clean rational function). The current check verifies the closed-form identity but does not connect `radcrit` to `T0**2 * R(alpha_crit)**2`. The connection itself — `T0**2 * R(alpha_crit)**2 == ((A+DK)**2 x0 + A**2 x1)**2` — is the substantive physical statement and is straightforward polynomial algebra well within `sp.expand` / `Expand` capability.

**Required change:**
Replace the hand-typed `radcrit` block with a derived check `T0**2 * R**2.subs(alpha, alpha_crit) == (A**2 x1 + (A+DK)**2 x0)**2` (or its Mathematica equivalent). See directive.

**Verification:**
After the patch, the new check appears in the saved output as e.g. `T0^2 * R(alpha_crit)^2 perfect square = 0` (sympy) and `PASS: T0^2 * R(alphaCrit)^2 perfect square` (mathematica). Scripts must exit 0.

### F4 — mathematica_transliteration

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage031_selected_branch_reachability_mathematica_audit.wl` (whole file)
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage031_selected_branch_reachability_sympy_audit.py` (whole file, for comparison)

**What's wrong:**
The Mathematica script is largely a line-by-line port of the SymPy script: identical PART I-V structure, identical assertion name strings ("ds_-/dalpha exact formula", "det factorization", "alpha_crit condition", "threshold radical square identity", "lambda_- * lambda_+ - T0*(alpha_crit-alpha)", "P0_- factorization"), the same variable naming choreography (`sigma`, `deltaKappa`, `kappaProd`, `r`, `lamMinus`, `lamPlus`, `sMinus`, `dsExact`, `dsExpected`, `p0Sel`, `t0`, `alphaCrit`, `lambdaPlusCrit`), the same identities chosen in the same order, and the identical hand-typed `radCrit` 9-term polynomial (lines 64-65 of `.wl` vs line 88 of `.py`). For example, compare:

SymPy lines 39-44:
```
sigma = x0 + x1
delta_kappa = x0 - x1
KappaProd = x0 * x1
R = sp.sqrt((DK + alpha * delta_kappa) ** 2 + 4 * alpha**2 * KappaProd)
lam_minus = sp.simplify((2 * A + DK - alpha * sigma - R) / 2)
lam_plus = sp.simplify((2 * A + DK - alpha * sigma + R) / 2)
```

Mathematica lines 32-38:
```
sigma = x0 + x1;
deltaKappa = x0 - x1;
kappaProd = x0*x1;
r = Sqrt[(dK + alpha*deltaKappa)^2 + 4*alpha^2*kappaProd];
lamMinus = FullSimplify[(2*a + dK - alpha*sigma - r)/2, ...];
lamPlus = FullSimplify[(2*a + dK - alpha*sigma + r)/2, ...];
```

The one substantive divergence is PART II, where Mathematica does the direct `D[p0Sel, alpha]` check that SymPy avoids (cf. F1). Outside PART II, the Mathematica file does not derive results independently; it transcribes the SymPy file's algebraic choreography.

**Why this matters:**
The two-engine policy requires each engine to derive results independently from the physical premises. Identical hand-typed polynomials in both files mean a single mistake propagates to both engines and goes undetected. The pattern also weakens the "engines agree" signal: of course they agree, they are computing the same expression strings.

**Required change:**
Replace the hand-typed `radCrit` polynomial in the Mathematica script with an in-script derivation from `t0^2 * r^2 /. alpha -> alphaCrit`. This is the same fix as F3 on the Mathematica side and is the minimal structural change that breaks the verbatim port without expanding scope. See directive for exact replacement.

**Verification:**
Mathematica file no longer contains the hand-typed 9-term polynomial; the saved Mathematica output continues to show `PASS: ds_-/dalpha exact formula` and the equivalent for all other assertions, plus the new derived check from F3.

## Independent-derivation check (Mathematica)

The Mathematica script is a structural transliteration of the SymPy script (see F4) except in PART II, where it does an honest direct check that SymPy ducks (`dPExact = D[p0Sel, alpha]` at `.wl:47` versus the abstract `Function`-and-substitute trick at `.py:60-69`). Both scripts share the identical hand-typed `radCrit` polynomial in PART IV, which is the clearest single piece of evidence that the Mathematica file did not derive results independently. The PART II divergence is the only place I can identify where the Mathematica file does original work.

## Engine cross-check

Both engines produce matching results at the level they assert.

PART I `ds_-/d alpha`:
- SymPy output: `(2*DK^2 * x0 * x1) / (4*alpha^2 * x0 * x1 + (DK + alpha*(x0 - x1))^2)^(3/2)`
- Mathematica output: `(2*dK^2*x0*x1)/((dK + alpha*(x0 - x1))^2 + 4*alpha^2*x0*x1)^(3/2)`
- Same modulo casing and term order.

PART IV `alpha_crit`:
- SymPy: `A*(A + DK) / (A*x1 + x0*(A + DK))`
- Mathematica: `(a*(a + dK)) / (dK*x0 + a*(x0 + x1))`
- These are equal: `A*x1 + x0*(A+DK) = A*x1 + A*x0 + DK*x0 = A*(x0+x1) + DK*x0`. ✓

PART IV `lambda_+^(crit)`:
- Mathematica: `((a + dK)^2*x0 + a^2*x1)/(dK*x0 + a*(x0 + x1))`
- SymPy: complicated form with `sqrt(radcrit)` not auto-denested but algebraically reduces to the same expression (verified by hand using `sqrt(radcrit) = (A+DK)^2 x0 + A^2 x1` and simplifying the numerator). ✓

Both saved outputs report `EXIT_CODE: 0` and PASS for every assertion. Engines agree.

## Verdict justification

Verdict is `findings`. Both engines run, agree on numerical/symbolic content, and assert PART I, III, V correctly. PART II of SymPy uses a generic calculus identity rather than the direct physical derivative (F1); the "alpha_crit condition" check in both engines is purely tautological (F2); the threshold radical square identity uses a hand-typed polynomial in both engines with no in-script connection to `R(alpha_crit)` (F3); the Mathematica script is a structural transliteration of the SymPy script outside PART II, with identical naming and identical hand-typed polynomial (F4). None of the findings invalidates the physical conclusions stated in the docstring — they all hold by the chain of identities present in the script — but the script as written does not robustly defend those conclusions against ordinary editing errors.

Attacks tried and failed: I verified by hand that (i) `ds_-/d alpha = 2 DK^2 x0 x1 / R^3` follows from differentiating `R^2 = (DK+alpha(x0-x1))^2 + 4 alpha^2 x0 x1` twice and using `R^2 (x0+x1)^2 - (R R')^2 = 4 DK^2 x0 x1`; (ii) `lambda_- lambda_+ = A(A+DK) - alpha T0` is exact (the `alpha^2` terms cancel via `(x0+x1)^2 - (x0-x1)^2 - 4 x0 x1 = 0`); (iii) `T0^2 R(alpha_crit)^2 = ((A+DK)^2 x0 + A^2 x1)^2` is exact (numerator of `R(alpha_crit)^2 * T0^2` factors as a perfect square); (iv) `lambda_-(alpha_crit) = 0` (substituting `alpha_crit` and `R(alpha_crit) = ((A+DK)^2 x0 + A^2 x1)/T0` into `(2A+DK - alpha sigma - R)/2` collapses to zero). All four substantive checks are real and the math holds — only the assertions need to be tightened.

No `UNFIXABLE` and no `CRITICAL_DOWNSTREAM`. The fixes do not change any quoted result; they tighten existing checks. Downstream units that consume `alpha_crit`, `ds_-/d alpha`, or `P0_-` will see identical numerical values.

## Self-test notes

I walked through each proposed directive change before writing it. For F1, I verified that the direct `sp.diff(P0_sel, alpha) - beta0*(ds_expected*lam_minus + s_minus**2)/lam_minus**2` evaluates to zero at `alpha = 0` (both sides give `beta0*(2*A*x0*x1/DK + x0**2)/A**2`), so the new assertion is genuinely satisfied by the script's actual quantities and is not vacuous. For F2, I verified `lam_minus.subs(alpha, alpha_crit) = 0` by direct hand substitution using `R(alpha_crit) = ((A+DK)^2 x0 + A^2 x1)/T0` and confirmed the numerator collapses to `A^2 (x0+x1) - A^2 x0 - A^2 x1 = 0` after cancellation. For F3, I verified `T0^2 * R(alpha_crit)^2 = ((A+DK)^2 x0 + A^2 x1)^2` symbolically and confirmed `sp.expand` will reduce the LHS-RHS to zero without needing aggressive simplification. For F4, the directive's change is identical to F3 on the Mathematica side, so the same self-test applies. No integrals, parity, or `sp.diff` against an absent variable are involved, so the variable-independence and symmetry traps do not apply here.
