---
unit_id: 049
batch: III.2
auditor_model: claude-opus-4-7[1m]
audit_date: 2026-05-26T00:00:00-06:00
verdict: findings
stop_cold: null
findings_count: 1
paper_alignment: aligned
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files:
    - /var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage049_dn_overlap_zeta.md
  paper_appendix: present
---

# Audit unit 049 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_049.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage049_dn_overlap_zeta.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part03.tex` (row 76 plus `\input` at line 216)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage049_dn_overlap_zeta_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage049_dn_overlap_zeta_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage049_dn_overlap_zeta_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage049_dn_overlap_zeta_mathematica_audit.txt`

## What the paper claims

The paper card (`stage_049.tex`) states `\stagefield{Output}{Microscopic D/N support ratios \eqref{eq:app-stage049-zeta-phys}--\eqref{eq:app-stage049-zeta-twin}}` and boxes four numbered results: (i) `zeta_n^{phys} = (K_W^{eff}/K_{phi,n}^{eff})(I_n/I_0)^2` (eq:app-stage049-zeta-phys); (ii) for a uniform local source density, `I_n/I_0 = 1/(2n+1)` (eq:app-stage049-In-ratio); (iii) for a same-operator twin family, `zeta_n^{twin} = 1/((2n+1)^2 [1 + x n(n+1)])` with `x = pi^2 T_X / (L^2 K_W^{eff})` (eq:app-stage049-zeta-twin); and (iv) `zeta_0^{twin} = 1` (eq:app-stage049-lowest-twin). `\stagefield{Inputs}` declares `chi_n(s) = sqrt(2/L) sin[(n+1/2) pi s/L]` and the source-to-mode overlaps `I_n` as the carried-forward quantities. The notes add: `I_n = sqrt(2L)/[(n+1/2) pi]`; the same-operator stiffnesses `K_W^{eff} = K_X + pi^2 T_X/(4L^2)` and `K_{phi,n}^{eff} = K_X + (n+1/2)^2 pi^2 T_X/L^2`, with the gap `K_{phi,n}^{eff} - K_W^{eff} = pi^2 T_X n(n+1)/L^2`. The appendix row at line 76 echoes the same three deliverables (physical support ratio, uniform-source overlaps, same-operator twin family).

## What the script claims to verify

The SymPy and Mathematica scripts declare (per docstring / banner) seven items: D/N half-wave momentum `k_n = (n+1/2) pi/L`; the uniform-source overlap `I_n = sqrt(2L)/((n+1/2) pi)` from direct integration of `chi_n(s) = sqrt(2/L) sin(k_n s)` on `[0, L]`; the ratio hierarchy `I_n/I_0 = 1/(2n+1)`; the microscopic coherent-support law `zeta_n^{phys} = (K_W^{eff}/K_{phi,n}^{eff})(I_n/I_0)^2`; the same-operator twin stiffness identity `K_{phi,n}^{twin} = K_W^{twin} + pi^2 T_X n(n+1)/L^2`; the twin-lane support ratio `zeta_n^{twin} = 1/((2n+1)^2 (1 + x n(n+1)))` with `x = pi^2 T_X/(L^2 K_W^{twin})`; and the `n=0` limit `zeta_0^{twin} = 1`. Both scripts also verify a D/N Neumann boundary check `cos(k_n L) = 0` (the v1 F1 patch).

## Paper ↔ script cross-check

| Paper-side deliverable | Script-side check | Status |
|---|---|---|
| (i) `zeta_n^{phys} = (K_W^{eff}/K_{phi,n}^{eff})(I_n/I_0)^2` | sympy line 84-87 / wl line 64-69 | match |
| (ii) `I_n/I_0 = 1/(2n+1)` for uniform sigma | sympy line 76-80 / wl line 58-62 | match |
| (iii) `zeta_n^{twin} = 1/((2n+1)^2 (1 + x n(n+1)))`, `x = pi^2 T_X/(L^2 K_W^{eff})` | sympy line 89-102 / wl line 71-81 | match |
| (iv) `zeta_0^{twin} = 1` | sympy line 103 / wl line 82 | match |
| Notes carry-forward: `I_n = sqrt(2L)/((n+1/2) pi)` | sympy line 70-74 / wl line 52-56 | match (sympy substantive; mathematica regressed — see F1) |
| Notes carry-forward: D/N quantization `k_n L = (n+1/2) pi` (Neumann at s=L) | sympy line 65-68 / wl line 50 | match (v1 F1 patch in place) |
| Notes: `K_{phi,n}^{eff} = K_W^{eff} + pi^2 T_X n(n+1)/L^2` | sympy line 89-95 / wl line 71-74 | match |

`paper_alignment: aligned`. All four boxed paper results and the supporting notes carry-forwards have corresponding script-side checks. No paper claim is missing a script check and no script check verifies an identity absent from the paper.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 65-68 | `simplify(cos(k_n L)) == 0` | D/N quantization (notes carry-forward) | yes (depends on n integer; cos((n+1/2)pi) = 0 non-trivially) |
| A2 | sympy | 70-74 | `simplify(integrate(chi_n,(s,0,L)) - sqrt(2L)/((n+1/2)pi)) == 0` | (ii) and notes `I_n` | yes (integral computed by sympy.integrate vs independent closed form) |
| A3 | sympy | 76-80 | `simplify(overlap_formula/I0 - 1/(2n+1)) == 0` | (ii) `I_n/I_0 = 1/(2n+1)` | yes (algebraic reduction) |
| A4 | sympy | 82-87 | `simplify(lambda_phi^2 K_W/(lambda_W^2 K_phi) - (K_W/K_phi)(1/(2n+1))^2) == 0` | (i) `zeta_n^{phys}` | yes (lambda_star cancels; substitution of A3 is non-trivial) |
| A5 | sympy | 89-95 | `simplify(K_phi_twin - K_W_twin - pi^2 T_X n(n+1)/L^2) == 0` | notes twin stiffness | yes (identity `(n+1/2)^2 - 1/4 = n(n+1)`) |
| A6 | sympy | 97-102 | `simplify((K_W_twin/K_phi_twin)(1/(2n+1))^2 - 1/((2n+1)^2(1+x n(n+1)))) == 0` | (iii) `zeta_n^{twin}` | yes (substantive link via x-form) |
| A7 | sympy | 103 | `twin_support_ratio(n,x).subs(n,0) - 1 == 0` | (iv) `zeta_0^{twin}=1` | yes |
| B1 | mathematica | 50 | `FullSimplify[Cos[kN l]] == 0` | D/N quantization | yes |
| B2 | mathematica | 52-56 | `FullSimplify[overlapFromIntegral - overlapFormula] == 0` | (ii) and notes `I_n` | **no** — both sides now `FullSimplify[Integrate[chiN,{s,0,l}],...]` differing only by Assumptions (see F1) |
| B3 | mathematica | 58-62 | `FullSimplify[ratio - 1/(2n+1)] == 0` | (ii) | yes |
| B4 | mathematica | 64-69 | `FullSimplify[zetaPhys - zetaPhysExpected] == 0` | (i) | yes |
| B5 | mathematica | 71-74 | `FullSimplify[kPhiTwin - kWTwin - Pi^2 tX n(n+1)/l^2] == 0` | notes twin stiffness | yes |
| B6 | mathematica | 76-81 | `FullSimplify[zetaTwin - twinSupportRatio[n,xExpr]] == 0` | (iii) | yes |
| B7 | mathematica | 82 | `FullSimplify[(twinSupportRatio[n,xExpr]/. n->0) - 1] == 0` | (iv) | yes |

## Findings

### F1 — tautological_check

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage049_dn_overlap_zeta_mathematica_audit.wl:52-56`

**What's wrong:**
The v1 F2 fix removed the `uniformDnOverlap` helper to break a transliteration mirror and replaced `overlapFormula` with a direct integrator call. The current code now reads:

```wolfram
overlapFromIntegral = FullSimplify[Integrate[chiN, {s, 0, l}]];
overlapFormula = FullSimplify[Integrate[chiN, {s, 0, l}], Assumptions -> Element[n, Integers] && n >= 0 && l > 0];
...
expectZero["uniform overlap integral", overlapFromIntegral - overlapFormula];
```

Both `overlapFromIntegral` and `overlapFormula` are now defined by the same expression — `FullSimplify[Integrate[chiN, {s, 0, l}], ...]` — differing only in whether explicit local `Assumptions` are passed. The subsequent `expectZero` then takes their difference and calls `FullSimplify` under `$Assumptions` (which already declares `n` integer and `l > 0`, per line 39-43). The two FullSimplify calls therefore operate on the same input under effectively the same assumption set; the residual is identically zero by construction and the assertion cannot fail. The paper-side claim — `I_n = sqrt(2L)/((n+1/2) pi)` (notes carry-forward, line 24-25 of notes) — is no longer verified on the Mathematica side because the script no longer contains an independent closed form to test against.

By contrast, the SymPy script still verifies this substantively: line 70-74 compares the integrator output against the independently declared closed form `sqrt(2L)/((n+1/2) pi)` returned by `uniform_dn_overlap`. So the two engines are no longer redundantly checking this carry-forward — only SymPy is.

**Why this matters:**
The Mathematica `expectZero["uniform overlap integral", ...]` assertion is the engine's single bottom-line check that the closed-form overlap `I_n = sqrt(2L)/((n+1/2) pi)` matches what the integrator returns. Eliminating the independent closed form converts the test into an integrator self-comparison: it can no longer catch a wrong closed form, because there is no closed form on the Mathematica side. The script still prints `I_n closed form = (2*Sqrt[2]*Sqrt[l])/(Pi + 2*n*Pi)` in the saved output, but that printed form is itself the integrator's output, not an independent target.

**Required change:**
Re-introduce an independent closed-form target on the Mathematica side without re-introducing the helper name `uniformDnOverlap` (to preserve the v1 F2 fix's goal of breaking the SymPy/Mathematica structural mirror). Concretely, in `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage049_dn_overlap_zeta_mathematica_audit.wl`, replace lines 52-53

Before:
```wolfram
overlapFromIntegral = FullSimplify[Integrate[chiN, {s, 0, l}]];
overlapFormula = FullSimplify[Integrate[chiN, {s, 0, l}], Assumptions -> Element[n, Integers] && n >= 0 && l > 0];
```

After:
```wolfram
overlapFromIntegral = FullSimplify[Integrate[chiN, {s, 0, l}], Assumptions -> Element[n, Integers] && n >= 0 && l > 0];
overlapFormula = Sqrt[2 l]/((n + 1/2) Pi);
```

This keeps the integrator call on the LHS of the comparison (the engine-independent computation from `chiN`) and supplies an independent closed-form expression on the RHS (the notes' boxed `I_n`). The subsequent `expectZero["uniform overlap integral", overlapFromIntegral - overlapFormula]` then becomes a substantive check: integrator output minus an independently stated closed form must FullSimplify to zero under the global integer assumption on `n`.

Note: the closed form `Sqrt[2 l]/((n + 1/2) Pi)` is not algebraically literally `uniformDnOverlap[n, l]` was — the deleted helper applied `FullSimplify` to the same expression, which under `$Assumptions` would not alter it. The substantive content (the integrand-independent target) is what F2 removed. Re-introducing the bare closed form (not a helper that mirrors the SymPy helper's name) restores the test while keeping the structural mirror broken.

**Verification:**
After the patch, the Mathematica script should (a) have `overlapFromIntegral` defined as the integrator call with explicit `Assumptions` and (b) have `overlapFormula` defined as a bare closed-form expression `Sqrt[2 l]/((n + 1/2) Pi)`, not as a second integrator call. The saved output should still show `I_n from direct integral` and `I_n closed form` as algebraically equivalent forms (e.g. `(2*Sqrt[2]*Sqrt[l])/(Pi + 2*n*Pi)` and `Sqrt[2 l]/((1/2 + n)*Pi)` differ only by trivial rearrangement), and `uniform overlap integral = 0` should still print with `PASS`. The script must still exit 0.

## Independent-derivation check (Mathematica)

The v1 F2 fix removed the `uniformDnOverlap` helper, breaking one strand of the structural mirror between the SymPy and Mathematica scripts. The remaining helpers — `overlapRatio[n_]` and `twinSupportRatio[n_, x_]` — are mathematical target forms (`1/(2 n + 1)` and `1/((2 n + 1)^2 (1 + x n (n+1)))`) rather than algebraic choreography; both engines have to assert these forms regardless of language. The genuinely independent piece is the symbolic integration of `chiN`, which Mathematica's `Integrate[]` performs by its own internal route. Assertions B3-B7 each independently substitute into and simplify expressions inside the Mathematica engine (`FullSimplify` under `$Assumptions`), so the engines' final outputs do not echo each other's intermediate steps for the zeta-related results. F1 above identifies a regression in B2 introduced by the v1 fix: the integrator self-comparison is no longer an independent derivation cross-check.

## Engine cross-check

Both engines report zero residuals on every assertion (saved outputs confirm `PASS` / `= 0` for the seven assertions in each). Closed-form algebraic agreement: SymPy `I_n = 2*sqrt(2)*sqrt(L)/(pi*(2n+1))` and Mathematica `(2*Sqrt[2]*Sqrt[l])/(Pi + 2*n*Pi)` are identical. SymPy `x = 4*pi^2 T_X/(4 K_X L^2 + pi^2 T_X)` and Mathematica `(4*Pi^2*tX)/(4*kX*l^2 + Pi^2*tX)` are identical. SymPy `zeta_n^{twin} = (4 K_X L^2 + pi^2 T_X)/((2n+1)^2(4 K_X L^2 + 4 pi^2 T_X n(n+1) + pi^2 T_X))` and Mathematica `1/((1 + 2 n)^2 (1 + (4 n (1 + n) Pi^2 tX)/(4 kX l^2 + Pi^2 tX)))` differ by a factor multiplication in the denominator but are algebraically identical (`(1 + a/b) = (a + b)/b`). No engine disagreement.

## Verdict justification

The paper-side claims are all exercised, and every script-side check traces to a specific paper or notes deliverable, so `paper_alignment` is aligned. The seven assertions in each engine substantively verify (1) the D/N quantization condition `cos(k_n L) = 0`, (2)/(3) the integral and ratio reductions, (4)/(6) the load-bearing `zeta_n^{phys}` and `zeta_n^{twin}` forms after non-trivial substitution, (5) the twin stiffness algebraic identity, and (7) the n=0 boundary value. The Mathematica scripts use Mathematica's integrator independently of SymPy's. The one remaining defect is that the v1 F2 fix accidentally turned the Mathematica `uniform overlap integral` assertion into an integrator self-comparison (both sides now derived from `Integrate[chiN, {s, 0, l}]`), so on the Mathematica side the closed-form overlap is no longer independently anchored. SymPy still verifies it. Verdict: findings, not stop-cold; the correction is a local edit that does not touch any downstream-visible result (the final `zeta` and `x` values are unchanged and remain verified by the other assertions).

## Self-test notes

For F1's replacement check: the integrand `Sqrt[2/l] Sin[(n+1/2) Pi s/l]` integrated over `s in [0, l]` produces `Sqrt[2/l] · (-Cos[(n+1/2) Pi] + Cos[0]) · l/((n+1/2) Pi)`; under the integer-n assumption `Cos[(n+1/2) Pi] = 0`, leaving `Sqrt[2/l] · l/((n+1/2) Pi) = Sqrt[2 l]/((n+1/2) Pi)`, which matches the proposed closed form. The trivial-case pre-check at n=0 gives `Sqrt[2 l]/(Pi/2) = 2 Sqrt[2 l]/Pi`, matching the saved `I_0` line. Variable independence: `overlapFormula = Sqrt[2 l]/((n + 1/2) Pi)` does depend on both `n` and `l`, so no spurious zero-derivative trap. Paper round-trip: the closed form `Sqrt[2 l]/((n + 1/2) Pi)` matches the notes line 28 `I_n = sqrt(2L) / [ (n+1/2) pi ]` exactly — no new paper_misalignment introduced. Path spec: the edit is in `mathematica/moving_throat_pde_stage049_dn_overlap_zeta_mathematica_audit.wl` only; no new files.
