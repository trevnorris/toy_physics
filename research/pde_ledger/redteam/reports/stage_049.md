---
unit_id: 049
batch: III.2
auditor_model: claude-opus-4-7[1m]
audit_date: 2026-05-22T00:00:00-06:00
verdict: findings
stop_cold: null
findings_count: 2
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: true
---

# Audit unit 049 red-team report

## Files reviewed

- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage049_dn_overlap_zeta_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage049_dn_overlap_zeta_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage049_dn_overlap_zeta_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage049_dn_overlap_zeta_mathematica_audit.txt`

## What the script claims to verify

The script claims (per docstring) seven results for a Dirichlet/Neumann (D/N) half-wave problem on [0, L]: (1) the half-wave momentum k_n = (n+1/2)π/L; (2) the closed-form uniform-source overlap I_n = √(2L)/((n+1/2)π) by direct integration of χ_n(s) = √(2/L) sin(k_n s) from 0 to L; (3) the ratio hierarchy I_n/I_0 = 1/(2n+1); (4) the microscopic coherent-support law ζ_n^(phys) = (K_W_eff/K_phi_eff)(I_n/I_0)²; (5) a same-operator twin stiffness identity K_phi_twin = K_W_twin + π²T_X n(n+1)/L²; (6) the twin-lane support ratio ζ_n^(twin) = 1/((2n+1)²(1 + x n(n+1))) with x = π²T_X/(L² K_W_twin); and (7) the n=0 limit ζ_0^(twin) = 1.

## Assertion inventory

| # | Script | Line | Form | Anchored to claim? |
|---|---|---|---|---|
| A1 | sympy | 65-68 | `simplify(k_n - (n+1/2)*pi/L) == 0` | no (tautological — `k_n` is defined as that very expression via the helper) |
| A2 | sympy | 70-74 | `simplify(integrate(chi_n,(s,0,L)) - sqrt(2L)/((n+1/2)pi)) == 0` | yes (independent symbolic integration vs closed form) |
| A3 | sympy | 76-80 | `simplify(overlap_formula/I0 - 1/(2n+1)) == 0` | yes (ratio reduction) |
| A4 | sympy | 82-87 | `simplify((λ*I_n)²K_W/((λ*I_0)²K_phi) - (K_W/K_phi)(1/(2n+1))²) == 0` | yes (cancels λ_star, recovers ratio²) |
| A5 | sympy | 89-95 | `simplify(K_phi_twin - K_W_twin - π²T_X n(n+1)/L²) == 0` | yes (algebraic identity (n+1/2)² = 1/4 + n(n+1)) |
| A6 | sympy | 97-102 | `simplify((K_W_twin/K_phi_twin)(1/(2n+1))² - twin_support_ratio(n,x)) == 0` | yes (links ratio to x-form) |
| A7 | sympy | 103 | `simplify(twin_support_ratio(n,x).subs(n,0) - 1) == 0` | yes (n=0 limit) |
| B1 | mathematica | 45 | `FullSimplify[kN - (n+1/2)Pi/l] == 0` | no (same tautology) |
| B2 | mathematica | 47-51 | `FullSimplify[Integrate[chiN,{s,0,l}] - ...] == 0` | yes |
| B3 | mathematica | 53-57 | `FullSimplify[ratio - 1/(2n+1)] == 0` | yes |
| B4 | mathematica | 59-64 | `FullSimplify[zetaPhys - zetaPhysExpected] == 0` | yes |
| B5 | mathematica | 66-69 | `FullSimplify[kPhiTwin - kWTwin - π² tX n(n+1)/l²] == 0` | yes |
| B6 | mathematica | 71-76 | `FullSimplify[zetaTwin - twinSupportRatio[n,xExpr]] == 0` | yes |
| B7 | mathematica | 77 | `FullSimplify[(twinSupportRatio[n,xExpr]/. n->0) - 1] == 0` | yes |

## Findings

### F1 — tautological_check

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage049_dn_overlap_zeta_sympy_audit.py:65-68`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage049_dn_overlap_zeta_mathematica_audit.wl:45`

**What's wrong:**
The first assertion in each engine reads (sympy):
```
def dn_halfwave_momentum(n, L):
    return sp.simplify((n + sp.Rational(1, 2)) * sp.pi / L)
...
k_n = dn_halfwave_momentum(n, L)
...
expect_zero("k_n - (n+1/2) pi / L", k_n - (n + sp.Rational(1, 2)) * sp.pi / L)
```
`k_n` is literally defined as `simplify((n+1/2)*pi/L)`. The assertion then computes `k_n - (n+1/2)*pi/L`, which is identically zero by construction — the assertion cannot fail regardless of physics. The Mathematica counterpart at line 45 has the identical structure: `kN = FullSimplify[(n+1/2) Pi/l]` followed by `expectZero["k_n - (n+1/2) pi / L", kN - (n+1/2) Pi/l]`. The docstring claim "(1) Exact D/N half-wave momentum and stiffness data" is therefore not actually exercised — the script tests `X − X = 0`.

**Why this matters:**
The docstring lists "exact D/N half-wave momentum" as the first claim. The current assertion does not exercise any physical content (it does not check that k_n satisfies the D/N boundary condition). A reader sees a green check that proves nothing about the wavenumber's origin.

**Required change:**
Replace the tautological identity with a substantive check that k_n in fact satisfies the D/N boundary condition. Since χ_n(s) = √(2/L) sin(k_n s) has Dirichlet boundary at s=0 (automatic) and Neumann boundary at s=L, the non-trivial requirement is `cos(k_n · L) = 0`. Substituting k_n*L = (n+1/2)π gives `cos((n+1/2)π) = 0` for integer n. In sympy: replace the existing `expect_zero(...)` block with `expect_zero("k_n satisfies D/N Neumann boundary", sp.cos(k_n * L))`. In Mathematica: replace `expectZero["k_n - (n+1/2) pi / L", kN - (n + 1/2) Pi/l]` with `expectZero["k_n satisfies D/N Neumann boundary", Cos[kN l]]`.

**Verification:**
After the patch, the sympy and mathematica scripts should each contain an assertion involving `cos(k_n*L)` (resp. `Cos[kN l]`) and no assertion of the form `k_n - (n+1/2)*pi/L`. The new check evaluates to 0 under `n` integer (since cos((n+1/2)π) = 0), so the script still exits 0; the saved output should show a printed line such as `k_n satisfies D/N Neumann boundary = 0`.

### F2 — mathematica_transliteration

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage049_dn_overlap_zeta_mathematica_audit.wl:26-77`

**What's wrong:**
The `.wl` script is a one-to-one structural port of the `.py` script: identical helper-function set, identical variable choreography, identical assertion sequence, identical printed banners, identical "Carry-forward formulas" block. Quoted corresponding sections:

Helpers (py:34-47 vs wl:26-29):
```
# sympy
def dn_halfwave_momentum(n, L): return sp.simplify((n + sp.Rational(1,2)) * sp.pi / L)
def uniform_dn_overlap(n, L):    return sp.simplify(sp.sqrt(2*L)/((n + sp.Rational(1,2))*sp.pi))
def overlap_ratio(n):            return sp.simplify(sp.Rational(1,1)/(2*n+1))
def twin_support_ratio(n,x):     return sp.simplify(1/((2*n+1)**2*(1 + x*n*(n+1))))
# mathematica
dnHalfwaveMomentum[n_,l_] := FullSimplify[(n+1/2) Pi/l];
uniformDnOverlap[n_,l_]   := FullSimplify[Sqrt[2 l]/((n+1/2) Pi)];
overlapRatio[n_]          := FullSimplify[1/(2 n + 1)];
twinSupportRatio[n_,x_]   := FullSimplify[1/((2 n + 1)^2 (1 + x n (n+1)))];
```

ζ_phys block (py:82-87 vs wl:59-64):
```
# sympy
lambda_W      = sp.simplify(lambda_star * I0)
lambda_phi_n  = sp.simplify(lambda_star * overlap_formula)
zeta_phys     = sp.simplify(lambda_phi_n**2 * KW_eff / (lambda_W**2 * Kphi_eff))
zeta_phys_exp = sp.simplify(KW_eff/Kphi_eff * overlap_ratio(n)**2)
# mathematica
lambdaW         = FullSimplify[lambdaStar i0];
lambdaPhi       = FullSimplify[lambdaStar overlapFormula];
zetaPhys        = FullSimplify[(lambdaPhi^2 kWeff)/(lambdaW^2 kPhiEff)];
zetaPhysExpected= FullSimplify[(kWeff/kPhiEff) overlapRatio[n]^2];
```

Twin-stiffness block (py:89-103 vs wl:66-77): identical sequence of `KW_twin`, `Kphi_twin`, `twin_gap`, `x_expr`, `zeta_twin`, `zeta_twin_expected`, with the n=0 substitution as the last assertion in both. The Mathematica script is not an independent re-derivation; it transcribes the SymPy workflow into Wolfram syntax. This violates the second-engine policy: both engines should derive the result from the physical premises by independent routes, not echo the same algebra.

**Why this matters:**
A transliteration cannot catch errors in the shared algebraic choreography. If the SymPy script chose a wrong intermediate identity or a wrong sign convention, the Mathematica run would reproduce the same wrong path and still report PASS. The "two engines agree" signal becomes only "the same algebra parses in two CAS dialects."

**Required change:**
In the Mathematica script, replace at least one shared intermediate step with an independent derivation that does not mirror the SymPy choreography. The minimal change with the least scope impact: derive `overlapFormula` not by re-stating the closed form via `uniformDnOverlap[n,l]` but by integrating χ_n symbolically and recording the result, then derive `i0` and the ratio from the integral result rather than from the helper. Concretely, in the `.wl` script: delete the helper `uniformDnOverlap[n_,l_] := FullSimplify[Sqrt[2 l]/((n+1/2) Pi)];` (line 27), set `overlapFormula = FullSimplify[Integrate[chiN, {s, 0, l}, Assumptions -> Element[n, Integers] && n >= 0 && l > 0]]` (replacing the line `overlapFormula = uniformDnOverlap[n, l];` at line 48), and set `i0 = FullSimplify[overlapFormula /. n -> 0]` (replacing the line `i0 = uniformDnOverlap[0, l];` at line 53). Leave the SymPy script unchanged. This forces the Mathematica overlap derivation to flow from the integrator rather than from a pre-stated formula echoing the SymPy helper.

**Verification:**
After the patch, the Mathematica script should (a) no longer contain the definition `uniformDnOverlap[n_, l_]`, (b) define `overlapFormula` directly as an `Integrate[...]` call with explicit integer assumption on n, and (c) define `i0` via substitution `n -> 0` on `overlapFormula`. The script must still exit 0 and the saved output should show `I_n closed form = (2*Sqrt[2]*Sqrt[l])/(Pi + 2*n*Pi)` (or an algebraically equivalent FullSimplify form) consistent with the SymPy result.

## Independent-derivation check (Mathematica)

The `.wl` script is a structural transliteration of the `.py` script: same helper-function names, same variable choreography, same assertion sequence, same banner and "Carry-forward formulas" text. See F2 for quoted corresponding sections. The only Mathematica step that is genuinely independent of SymPy's algebra is the symbolic integration of χ_n (line 47), since the integrator implementations differ; but the surrounding algebra is mirrored.

## Engine cross-check

Both engines produce identical zero residuals across all seven assertions:

| Assertion | sympy residual | mathematica residual |
|---|---|---|
| k_n − (n+1/2)π/L | 0 | 0 |
| uniform overlap integral | 0 | 0 |
| overlap ratio hierarchy | 0 | 0 |
| microscopic coherent-support law | 0 | 0 |
| same-operator twin stiffness relation | 0 | 0 |
| exact twin-lane support ratio | 0 | 0 |
| lowest twin half-wave | 0 | 0 |

Closed forms agree: SymPy `I_n = 2√2·√L/(π(2n+1))`, Mathematica `(2 Sqrt[2] Sqrt[l])/(Pi + 2 n Pi)` — algebraically identical. SymPy `x = 4π²T_X/(4K_X L² + π²T_X)`, Mathematica `(4 Pi² tX)/(4 kX l² + Pi² tX)` — identical. No engine disagreement.

## Verdict justification

Six of seven assertions in each engine substantively exercise their stated claims (closed-form vs. direct integration; ratio reduction; stiffness identity (n+1/2)² = 1/4 + n(n+1); n=0 limit). The first assertion in each engine is tautological by construction — k_n is defined as the very expression it is compared to — so the "exact D/N half-wave momentum" claim is not actually tested (F1). The `.wl` script is a structural transliteration of the `.py` script and is not an independent derivation (F2). Engines agree numerically and the saved outputs are fresh (output mtimes May 11, script mtimes Apr 21). Verdict: findings, not stop-cold; the corrections are local and do not propagate to downstream units (k_n and ζ_n forms themselves are unchanged — only the verification mechanism is strengthened).

## Self-test notes

For F1's replacement check, the trivial-case pre-check passes: with n=0, `cos(k_n·L) = cos(π/2) = 0`; with n=1, `cos(3π/2) = 0`. The check depends on `n` being integer (existing assumption in both scripts), and after substituting k_n·L = (n+1/2)π the expression is independent of L — so there is no spurious variable-independence issue. For F2's replacement of `overlapFormula` via `Integrate`, the integrand sin((n+1/2)π s/l) is odd in s about s = l/(2(n+1/2)) but the integration domain [0, l] is not symmetric, so the integral is not forced to zero by parity; the integer assumption on n is required (otherwise cos((n+1/2)π) ≠ 0 and the closed form would acquire an extra `(1 − cos((n+1/2)π))` factor). Path specs: F1 edits the existing `.py` and `.wl` files in `scripts/` and `mathematica/` respectively; F2 edits only the `.wl` in `mathematica/` — no new files.
