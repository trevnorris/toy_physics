---
unit_id: 151
batch: IV.6
auditor_model: claude-opus-4-7[1m]
audit_date: 2026-05-28T03:22:59Z
verdict: findings
stop_cold: null
findings_count: 3
paper_alignment: partial
scripts_checked:
  sympy: insufficient
  mathematica: insufficient
  engines_agree: true
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files:
    - /var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage151_first_order_selected_correction.md
  paper_appendix: present
---

# Audit unit 151 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_151.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage151_first_order_selected_correction.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part04.tex` (only `\input{stages/stage_151}` line at L1336; no extra summary text)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage151_first_order_selected_correction_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage151_first_order_selected_correction_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage151_first_order_selected_correction_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage151_first_order_selected_correction_mathematica_audit.txt`

## What the paper claims

The stage's quoted body claim is: *"Actual selected deformation is the centered residual \(R_*\); correction is fixed by two covariances."* From the notes, the explicit deliverables are: (a) the normalized first-order correction is `δΣ_act(x) = Σ_*(x)·[1 - tildeR_*(x)]` where `tildeR_*(x) := R_*(x) - <R_*>_*` (i.e., the centered residual); (b) the two covariance moment shifts `δg_act = ∫₀¹ c(x) δΣ_*(x) dx = -Cov_*(c, R_*)` and `δS_act = ∫₀¹ K_q(x) δΣ_*(x) dx = -Cov_*(K_q, R_*)` with `c(x)=cos(πx/2)`, `K_q(x)=cosh(π(1-x)/2)/cosh(π/2)`; (c) the bias retuning `δΠ_act = -δg_act/g_*' = Cov_*(c,R_*)/g_*'`; (d) the traction retuning `δT̂_m,act = A_T·δg_act + B_T·δS_act = -A_T·Cov_*(c,R_*) - B_T·Cov_*(K_q,R_*)`. The card's `Checks` list further demands verification that (i) profile deformations are centered before covariance formulas are applied, (ii) the rigidity kernel (not branch ambiguity) controls non-exponential corrections, (iii) the correction stays inside the reduced mouth-layer regime.

## What the script claims to verify

Both engines declare five abstract real symbols (`Rbar, cbar, Kbar, cRbar, KRbar`) standing in for the moments `<R>_*, <c>_*, <K>_*, <cR>_*, <KR>_*`, define `CovcR := cRbar - cbar*Rbar` and `CovKR := KRbar - Kbar*Rbar` to play the role of the two covariances, then *set* `delta_g := -(cRbar - cbar*Rbar)` and `delta_S := -(KRbar - Kbar*Rbar)` and assert `delta_g + CovcR = 0` and `delta_S + CovKR = 0`. After that they compute `deltaPi = -delta_g/gprime` and `deltaT = AT*delta_g + BT*delta_S` and merely `print` the results without assertion. The script's stated theorem in its closing print is the third paper deliverable. There is no integration step, no `c(x)`/`K_q(x)` symbols, no `Σ_*(x)`, and no check of the centering property `<δΣ>_* = 0`.

## Paper ↔ script cross-check

| Paper deliverable | Script-side check | Status |
|---|---|---|
| (a) Centered residual `tildeR_* = R_* - <R_*>_*` and `δΣ = -Σ_* tildeR_*` | None — no `Σ_*`, no `<R>_*`, no expansion verifying centering produces `<δΣ>_* = 0` | missing |
| (b1) `δg_act = ∫ c δΣ dx = -Cov(c,R)` (derived from the integral) | Script *defines* `delta_g = -CovcR` and asserts `delta_g + CovcR = 0` — by construction | mismatch (tautological substitute) |
| (b2) `δS_act = ∫ K_q δΣ dx = -Cov(K_q,R)` | Same — defines `delta_S = -CovKR` and asserts `delta_S + CovKR = 0` | mismatch (tautological substitute) |
| (c) `δΠ_act = -δg_act/g_*' = Cov(c,R)/g_*'` | Printed (`deltaPi`), not asserted | missing assertion |
| (d) `δT̂_m,act = A_T·δg + B_T·δS` | Printed (`deltaT`), not asserted | missing assertion |
| Check (i) centering before covariance | Not exercised | missing |
| Check (ii) rigidity kernel controls correction | Not exercised | missing |
| Check (iii) one-step correction stays within reduced regime | Not exercised | missing |

`paper_alignment` is `partial`: the script names the right covariance algebra but verifies none of the paper's integration / centering / rigidity-kernel content.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 39 | `expect_zero("delta g + Cov(c,R)", delta_g + CovcR)` where `delta_g = -(cRbar - cbar*Rbar)` and `CovcR = cRbar - cbar*Rbar` | Should exercise (b1) — `δg = -Cov(c,R)`; in practice none (identity by construction) | no |
| A2 | sympy | 40 | `expect_zero("delta S + Cov(K,R)", delta_S + CovKR)` similarly | Should exercise (b2); identity by construction | no |
| A3 | mathematica | 37 | `expectZero["delta g + Cov(c,R)", deltaG + covCR]` | mirror of A1 | no |
| A4 | mathematica | 38 | `expectZero["delta S + Cov(K,R)", deltaS + covKR]` | mirror of A2 | no |

There are no other assertions; `deltaPi` and `deltaT` are printed only.

## Findings

### F1 — tautological_check

**Severity:** high
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage151_first_order_selected_correction_sympy_audit.py:30-40`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage151_first_order_selected_correction_mathematica_audit.wl:31-38`

**What's wrong:**
The SymPy script defines

```
CovcR = cRbar - cbar*Rbar          # line 30
CovKR = KRbar - Kbar*Rbar          # line 31
delta_g = -(cRbar - cbar*Rbar)     # line 36
delta_S = -(KRbar - Kbar*Rbar)     # line 37
```

and then asserts

```
expect_zero("delta g + Cov(c,R)", delta_g + CovcR)   # line 39
expect_zero("delta S + Cov(K,R)", delta_S + CovKR)   # line 40
```

`delta_g + CovcR` is, by literal substitution of the defining expressions, `-(cRbar - cbar*Rbar) + (cRbar - cbar*Rbar) = 0`. The same is true for `delta_S + CovKR`. These checks cannot fail under any choice of the underlying physics — they are algebraic tautologies guaranteed by the way `delta_g`/`delta_S` are defined. The Mathematica script does the exact same thing at lines 31–38 with renamed symbols. The paper's actual deliverables (b1) and (b2) require *deriving* `δg = ∫₀¹ c(x) δΣ_*(x) dx = -Cov_*(c,R_*)` and similarly for `δS`, starting from the centered linearized correction `δΣ = -Σ_*·(R - <R>_*)`. The scripts skip the derivation step entirely and assert what they just wrote down.

**Why this matters:**
The two only assertions in this stage are guaranteed to pass irrespective of whether the paper's identities are correct. If the linearization sign were wrong, if the centering convention were inconsistent with `δg = ∫ c δΣ`, if the covariance formula were missing a normalization — none of those errors would be caught. The "PASS" transcript is meaningless evidence.

**Required change:**
Introduce a concrete integrable parametrization of the underlying functions and derive the two covariance identities by integration, then assert the derived expressions equal `-Cov(c,R)` and `-Cov(K_q,R)`. The minimal non-tautological version uses the paper's exact definitions:

- Define a symbolic `Pi_star` (positive real) and `Sigma_star(x) = Pi_star * exp(-Pi_star*x) / (1 - exp(-Pi_star))` on `x ∈ [0,1]`.
- Define `c(x) = cos(pi*x/2)` and `K_q(x) = cosh(pi*(1-x)/2)/cosh(pi/2)`.
- Define a test residual `R(x)` symbolic (e.g., `r1*x + r2*x**2` with `r1, r2` symbolic reals) so the covariances are nontrivial functions of `Pi_star, r1, r2`.
- Define `<f>_* := integrate(Sigma_star(x)*f(x), (x, 0, 1))`.
- Compute `Rbar = <R>_*`, `cbar = <c>_*`, `Kbar = <K_q>_*`, `cRbar = <c*R>_*`, `KRbar = <K_q*R>_*`.
- Compute `delta_Sigma = -Sigma_star * (R - Rbar)`.
- Assert `<delta_Sigma>_* == 0` (centering — paper Check i).
- Compute `delta_g_int = integrate(c(x)*delta_Sigma, (x,0,1))` and assert `simplify(delta_g_int - (-(cRbar - cbar*Rbar))) == 0`.
- Compute `delta_S_int = integrate(K_q(x)*delta_Sigma, (x,0,1))` and assert `simplify(delta_S_int - (-(KRbar - Kbar*Rbar))) == 0`.
- Then add the downstream-form asserts `simplify(deltaPi - CovcR/gprime) == 0` and `simplify(deltaT - (-AT*CovcR - BT*CovKR)) == 0` (these are themselves nontrivial because `deltaPi` and `deltaT` are now functions of the actually integrated `delta_g`, `delta_S`).

Mirror the same derivation independently in the Mathematica script (using `Integrate[..., {x, 0, 1}]` and `Cos[Pi*x/2]`, `Cosh[Pi*(1-x)/2]/Cosh[Pi/2]`), not a line-for-line transliteration of the SymPy.

**Verification:**
After the rewrite, `redteam exec-sympy 151` and `redteam exec-mathematica 151` must run to exit 0 and the saved transcripts must show non-zero intermediate residuals (specifically, `delta_g_int` should equal a nontrivial symbolic combination of `Pi_star, r1, r2` before simplification cancels it against `-CovcR`). Each of the new assertions `<delta_Sigma>_* == 0`, `delta_g_int + CovcR == 0`, `delta_S_int + CovKR == 0`, `deltaPi - CovcR/gprime == 0`, `deltaT - (-AT*CovcR - BT*CovKR) == 0` must be present in the script source.

### F2 — insufficient_verification

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage151_first_order_selected_correction_sympy_audit.py:42-48`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage151_first_order_selected_correction_mathematica_audit.wl:40-46`

**What's wrong:**
After computing `deltaPi = sp.simplify(-delta_g / gprime)` and `deltaT = sp.simplify(AT*delta_g + BT*delta_S)`, the script only `sp.pprint`s them. There is no assertion that `deltaPi` equals the paper's stated form `Cov(c,R)/g_*'` (notes Section 3, boxed equation) or that `deltaT` equals `-A_T·Cov(c,R) - B_T·Cov(K_q,R)` (notes Section 3, boxed equation). The Mathematica script is the same. So even setting aside F1, the script's own stated theorem ("the selected first-order source correction is completely determined by Cov_*(c,R_*) and Cov_*(K_q,R_*)") is only `print`-asserted, never checked.

**Why this matters:**
A sign flip or factor on `deltaT` would not be caught. The paper's deliverables (c) and (d) — the bias retuning and the traction retuning — are not exercised by any assertion at all.

**Required change:**
After the integrated derivation in F1, add the explicit assertions:
- SymPy: `expect_zero("deltaPi - Cov(c,R)/g'", deltaPi - CovcR/gprime)` and `expect_zero("deltaT + AT*CovcR + BT*CovKR", deltaT + AT*CovcR + BT*CovKR)`.
- Mathematica: `expectZero["deltaPi - Cov(c,R)/gPrime", deltaPi - covCR/gPrime]` and `expectZero["deltaT + aT*covCR + bT*covKR", deltaT + aT*covCR + bT*covKR]`.

**Verification:**
Each of the new `expect_zero`/`expectZero` lines must appear in the script and produce a `= 0` line in the saved transcript.

### F3 — mathematica_transliteration

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage151_first_order_selected_correction_mathematica_audit.wl:28-46`
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage151_first_order_selected_correction_sympy_audit.py:26-48`

**What's wrong:**
The Mathematica script is a line-by-line transliteration of the SymPy script. Compare:

SymPy (py:26–37):
```
Rbar, cbar, Kbar, cRbar, KRbar = sp.symbols("Rbar cbar Kbar cRbar KRbar", real=True)
gprime, AT, BT = sp.symbols("gprime AT BT", real=True, nonzero=True)
CovcR = cRbar - cbar*Rbar
CovKR = KRbar - Kbar*Rbar
delta_g = -(cRbar - cbar*Rbar)
delta_S = -(KRbar - Kbar*Rbar)
```

Mathematica (wl:28–35):
```
Clear[rBar, cBar, kBar, cRBar, kRBar, gPrime, aT, bT];
$Assumptions = Element[{rBar, cBar, kBar, cRBar, kRBar, gPrime, aT, bT}, Reals] && gPrime != 0;
covCR = cRBar - cBar*rBar;
covKR = kRBar - kBar*rBar;
deltaG = -(cRBar - cBar*rBar);
deltaS = -(kRBar - kBar*rBar);
```

The variable choreography, the order of definitions, the sign convention, and the assertion calls are identical (lines 39–40 of the SymPy and 37–38 of the WL are matched one-to-one). Neither engine derives anything independently from the paper's premises (`Σ_*(x) = Π_*e^{-Π_*x}/(1-e^{-Π_*})`, `c(x) = cos(πx/2)`, `K_q(x) = cosh(π(1-x)/2)/cosh(π/2)`); both engines just rename the same five expectation-value symbols and assert the same algebraic identity. This violates the second-engine policy: both engines must derive the result independently from the physical premises, not echo each other's algebra.

**Why this matters:**
A Mathematica transliteration of a SymPy script gives no independent corroboration — it confirms only that the same five-line algebraic identity holds in two CAS engines. Any error in the symbol setup propagates identically to both engines.

**Required change:**
When rewriting under F1 and F2, give the Mathematica script its own independent derivation path. Concretely: in Mathematica use `Integrate[Sigma[x]*f[x], {x, 0, 1}]` with `Sigma[x_] := Pi_* Exp[-Pi_* x]/(1 - Exp[-Pi_*])`, and define `c[x_] := Cos[Pi*x/2]`, `Kq[x_] := Cosh[Pi*(1-x)/2]/Cosh[Pi/2]` (the SymPy version uses `cos(pi*x/2)` and the equivalent `Lambdify`-friendly form). Use Mathematica's `Series[..., {epsilon, 0, 1}]` to expand `Exp[-(Pi_* x + epsilon R[x])]/(...normalization...)` to first order in `epsilon` as an independent re-derivation of the linearized correction `δΣ = -Σ_*(R - <R>_*)` (the SymPy version may instead expand by hand or via `sp.series`). The two engines should reach the same boxed identities by different code paths.

**Verification:**
Reviewer reads the two scripts side by side; corresponding sections should differ in approach (one expands by series, the other by direct substitution; one uses `Integrate`, the other `sp.integrate`), not just in syntax. Both transcripts must still report `delta_g_int + CovcR == 0` etc.

## Independent-derivation check (Mathematica)

The Mathematica script does not independently re-derive anything from the physical premises (`Σ_*(x)`, `c(x)`, `K_q(x)`). It defines the same five expectation-value symbols and the same two covariance combinations as the SymPy script, then asserts the identical algebraic identity. The corresponding SymPy lines 26–40 and Mathematica lines 28–38 are in one-to-one correspondence. This is captured by finding F3.

## Engine cross-check

Both engines print/assert the same final expressions and the saved transcripts agree:

SymPy:
```
delta g + Cov(c,R) = 0
delta S + Cov(K,R) = 0
deltaPi = (-Rbar*cbar + cRbar)/gprime
deltaT  = AT*(Rbar*cbar - cRbar) - BT*(KRbar - Kbar*Rbar)
```

Mathematica:
```
delta g + Cov(c,R) = 0
delta S + Cov(K,R) = 0
deltaPi = (cRBar - cBar*rBar)/gPrime
deltaT  = -(aT*cRBar) - bT*kRBar + aT*cBar*rBar + bT*kBar*rBar
```

Algebraically `-Rbar*cbar + cRbar = cRBar - cBar*rBar` and `AT*(Rbar*cbar - cRbar) - BT*(KRbar - Kbar*Rbar) = -aT*cRBar - bT*kRBar + aT*cBar*rBar + bT*kBar*rBar`. Engines agree; the issue is that they agree on a tautology, not that they agree about anything substantive (F1).

## Verdict justification

The two scripts assert exactly the two same tautological identities and both pass — but neither identity exercises any of the paper's stated deliverables. The paper card and notes require: a centering identity, two integral-derived covariance moment shifts, a bias retuning, and a traction retuning. The scripts cover none of these by assertion: they only state the covariance algebra as a definitional rearrangement and check the rearrangement is self-consistent. Additionally, the Mathematica script is a one-to-one transliteration of the SymPy script and provides no independent corroboration. Outputs are fresh (sympy script 11:56 → output 12:46; mathematica script 11:56 → output 13:16). No `stop_cold`: nothing about the paper is mathematically inconsistent, and no downstream stage will silently break from re-doing the audit honestly (the *form* of the boxed identities is already what the paper states; the issue is that nobody has checked them in-script).

## Self-test notes

I checked: (1) variable independence — the proposed new check uses `Integrate[..., {x, 0, 1}]`, so `x` is a bound variable and `Pi_*, r1, r2` are real free symbols; the integrand `Σ_*(x)·R(x)` is a proper function of `x` so the integral is well-defined and yields nontrivial expressions in `Pi_*, r1, r2`. (2) Trivial-case pre-check — with `R(x) ≡ 0`, `Rbar = 0`, `cRbar = 0`, `KRbar = 0`, all four asserted residuals are identically 0, but for `R(x) = r1*x + r2*x²` with generic `r1, r2`, `delta_g_int` evaluates to a nontrivial expression in `Pi_*, r1, r2` that cancels against `-CovcR` only after `simplify`, so the assertion is non-tautological. (3) Paper round-trip — the prescribed fix uses exactly the paper-stated `Σ_*(x)`, `c(x)`, `K_q(x)`, and the boxed identities from notes Sections 2–3, so it introduces no new constants nor any new symbols absent from the paper.
