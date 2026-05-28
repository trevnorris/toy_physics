---
unit_id: 146
batch: IV.5
auditor_model: claude-opus-4-7
audit_date: 2026-05-27T00:00:00Z
verdict: findings
stop_cold: null
findings_count: 3
paper_alignment: partial
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files:
    - /var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage146_positive_deformation_expansion.md
  paper_appendix: present
---

# Audit unit 146 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_146.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage146_positive_deformation_expansion.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part04.tex` (rows referencing stage 146 read)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage146_positive_deformation_expansion_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage146_positive_deformation_expansion_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage146_positive_deformation_expansion_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage146_positive_deformation_expansion_mathematica_audit.txt`

## What the paper claims

Stage 146 proves the first-order, positive convex-deformation rigidity of the regular Family-1 canonical mouth branch. The paper card states it is a "finite mouth-profile corrections ledger step" and that "First-order deformations depend only on \((\bar g_\varsigma,\bar S_\varsigma)\) and retune \(\Pi\)." The notes give the substantive deliverables: (a) the convex family \(\Sigma_\epsilon = (1-\epsilon)\Sigma_* + \epsilon\varsigma\) preserves positivity and normalization; (b) the source enters only through the two moments \(\bar g[\sigma]\) and \(\bar S[\sigma]\); (c) the affine laws \(\bar g_\epsilon = \mathfrak g_* + \epsilon(\bar g_\varsigma - \mathfrak g_*)\) and \(\bar S_\epsilon = \mathcal S_* + \epsilon(\bar S_\varsigma - \mathcal S_*)\); (d) the first-order compensation law \(\delta\Pi = -\epsilon(\bar g_\varsigma - \mathfrak g_*)/\mathfrak g'_*\); and (e) the compensated mixed-channel response \(\delta\mathcal S_q = \epsilon\bigl[(\bar S_\varsigma - \mathcal S_*) - (\mathcal S'_*/\mathfrak g'_*)(\bar g_\varsigma - \mathfrak g_*)\bigr]\), with the boxed numerics \(\Pi_*\approx 1.50882951349316\), \(\mathfrak g_*\approx 0.758035078944663\), \(\mathcal S_*\approx 0.658075937605429\), \(\mathfrak g'_*\approx 0.0714453558083195\), \(\mathcal S'_*\approx 0.0483709542125041\). The card's three "Checks" items name centering before covariance, the rigidity kernel as the corrections controller, and persistence within the reduced mouth-layer regime.

## What the script claims to verify

The two scripts define the closed-form moments \(g(\Pi) = 2\Pi(2\Pi e^\Pi + \pi)/((4\Pi^2+\pi^2)(e^\Pi-1))\) and \(\mathcal S_q(\Pi) = \Pi(\kappa\tanh\kappa + \Pi(e^{-\Pi}\sech\kappa - 1))/((1-e^{-\Pi})(\kappa^2-\Pi^2))\) with \(\kappa=\pi/2\), perform a 3-sample numerical kernel cross-check against direct integration at \(\Pi\in\{1,3/2,5/2\}\), solve for \(\Pi_*\) via `nsolve`/`FindRoot` against the Family-1 lower-branch value \(\mathfrak g_-\), print \(\Pi_*, g_*, \mathcal S_*, g'_*, \mathcal S'_*\), then construct symbolic expressions for \(\delta\Pi\) and \(\delta\mathcal S\) and print them. The only `expect_zero` assertions are the two affine-law identities `g_eps - (gminus + eps*(gbar - gminus))` and `S_eps - (sStar + eps*(sBar - sStar))`. The Mathematica script additionally asserts `g(Pi) direct-formula = 0` (closed-form vs. symbolic integration of `sigma*Cos[Pi*x/2]`) and `Pi_* compensation point diff <= 10^-20` (numeric).

## Paper ↔ script cross-check

| Paper deliverable | Script-side check | Status |
|---|---|---|
| Closed form for \(g(\Pi)\) matches \(\int_0^1 \Sigma_* c\) | Mathematica `expectZero["g(Pi) direct-formula"]` (line 48); SymPy: no symbolic check, only numeric via the 3-point kernel test (which uses `(Sigma*Kq)`, not `(Sigma*c)`) | partial |
| Closed form for \(\mathcal S_q(\Pi)\) matches \(\int_0^1 \Sigma_* K_q\) | 3-point numerical comparison in both engines (`kernel check at Pi=...`); no symbolic verification | partial |
| Affine laws \(\bar g_\epsilon = \mathfrak g_* + \epsilon(\bar g_\varsigma - \mathfrak g_*)\), \(\bar S_\epsilon = \ldots\) | Both engines `expect_zero` on the very identity, but `g_eps` and `S_eps` are *defined* as `(1-eps)*gminus + eps*gbar` etc., so the assertion is algebraically guaranteed | tautological (see F1) |
| First-order compensation \(\delta\Pi = -\epsilon(\bar g - \mathfrak g_*)/\mathfrak g'_*\) | Constructed by hand in both engines and printed (sympy line 62, math line 78). No assertion. | missing (no assert) |
| Mixed-channel response \(\delta\mathcal S_q = \epsilon[(\bar S - \mathcal S_*) - (\mathcal S'_*/\mathfrak g'_*)(\bar g - \mathfrak g_*)]\) | Built as `dS = eps*(Sbar - sStar) + sPrimeStar*dPi`; printed, not asserted | missing (no assert) |
| Numerical values \(\Pi_*, g_*, \mathcal S_*, g'_*, \mathcal S'_*\) | Computed and printed; Mathematica also asserts the compensation residual `gStar - gMinus` is below 10^-20 | partial (sympy only prints) |
| Card "Check": "first-order profile deformations are centered before covariance formulas are used" | Not explicitly tested (no assertion on \(\int_0^1 (\varsigma - \Sigma_*) dx = 0\) or on centering of the deformation) | missing |
| Card "Check": "rigidity kernel ... controls non-exponential corrections" | Not assertion-tested; the kernel \(K_q\) is used but no check of its rigidity vs. branch alternatives | missing |
| Card "Check": "one-step nonlinear correction remains within the reduced mouth-layer regime" | Not addressed (no \(O(\epsilon^2)\) sanity check or range bound) | missing |

Dominant pattern is "partial" / "missing assert" with two tautological assertions and one substantive Mathematica-only assert. Setting `paper_alignment: partial`.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 42-43 | `if abs(float(diff)) > 1e-12: raise AssertionError` (kernel `Sigma*Kq` integral vs. `Sformula`, 3 samples) | \(\mathcal S_q(\Pi)\) closed-form correctness, numerically | partial (numeric, no symbolic) |
| A2 | sympy | 74 | `expect_zero("g_eps affine law", sp.expand(g_eps - (gminus + eps*(gbar-gminus))))` | affine law for \(\bar g_\epsilon\) | no — tautological (g_eps defined as `(1-eps)*gminus + eps*gbar`) |
| A3 | sympy | 75 | `expect_zero("S_eps affine law", ...)` | affine law for \(\bar S_\epsilon\) | no — tautological |
| A4 | mathematica | 48 | `expectZero["g(Pi) direct-formula", gDirect - gFormula]` (symbolic integration of `sigma*Cos[Pi*x/2]` vs. closed form) | \(g(\Pi)\) closed-form correctness | yes |
| A5 | mathematica | 51-59 | `expectApprox["kernel formula sample ...", numInt, numFormula, 10^-12]` (3 samples) | \(\mathcal S_q(\Pi)\) numeric closed-form | partial |
| A6 | mathematica | 76 | `expectApprox["Pi_* compensation point", gStar, N[gMinus, 40], 10^-20]` | \(\Pi_*\) satisfies \(g(\Pi_*)=\mathfrak g_*\) | yes |
| A7 | mathematica | 87 | `expectZero["g_eps affine law", Expand[gEps - (gMinus + eps*(gBar - gMinus))]]` | affine law for \(\bar g_\epsilon\) | no — tautological |
| A8 | mathematica | 88-90 | `If[TrueQ[resSEps === 0], pass, fail]` on `sEps - (sStar + eps*(sBar - sStar))` | affine law for \(\bar S_\epsilon\) | no — tautological |

Non-tautological coverage: A1 (numeric \(\mathcal S_q\), 3 pts), A4 (symbolic \(g(\Pi)\) Mathematica only), A5 (numeric \(\mathcal S_q\), 3 pts), A6 (numeric \(\Pi_*\)).

## Findings

### F1 — tautological_check

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage146_positive_deformation_expansion_sympy_audit.py:72-75`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage146_positive_deformation_expansion_mathematica_audit.wl:85-90`

**What's wrong:**
The SymPy script defines
```
g_eps = (1-eps)*gminus + eps*gbar
```
on line 72 and then asserts
```
expect_zero("g_eps affine law", sp.expand(g_eps - (gminus + eps*(gbar-gminus))))
```
on line 74. The RHS `gminus + eps*(gbar - gminus)` is `(1-eps)*gminus + eps*gbar` by distributivity, so the residual is identically zero. Line 75 does the same with `S_eps` and `sStar`/`Sbar`. The Mathematica mirror on lines 85-90 repeats this with `gEps = (1-eps)*gMinus + eps*gBar` and tests against `(gMinus + eps*(gBar - gMinus))`. These assertions cannot fail no matter what the underlying physics says; they are pure restatements of the distributive law.

The paper notes' boxed affine laws (notes section 2) are real claims — they say the *physical* \(\bar g_\epsilon = \int \Sigma_\epsilon c\) equals \(\mathfrak g_* + \epsilon(\bar g_\varsigma - \mathfrak g_*)\) where \(\bar g_\varsigma = \int \varsigma\, c\). The substantive content is that the convex combination of integrals is the same as the integral of the convex combination plus that \(\bar g_*[\Sigma_*] = \mathfrak g_*\). The script never builds \(\bar g_\epsilon\) as `Integrate[(1-eps)*Sigma_* + eps*varsigma * c, ...]` and compares against `gminus + eps*(gbar - gminus)`; it just rewrites the algebraic identity.

**Why this matters:**
A reader of the saved output sees "g_eps affine law = 0 / PASS" and "S_eps affine law = 0 / PASS" and reasonably concludes the paper's headline affine laws have been independently checked. They have not — the assertions are arithmetic identities of the script's own definitions, not verifications of the physical claim.

**Required change:**
Build \(\bar g_\epsilon\) from the physical definition (integral of the convex source against the kernel) and test that it equals the affine RHS. Concretely, introduce a symbolic centered deformation `varsigma_test(x)` (or two — e.g., uniform `1` and a polynomial like `6*x*(1-x)` that are both positive and integrate to 1 on [0,1]), define
```
gbar_phys = sp.integrate(((1-eps)*Sigma + eps*varsigma_test)*sp.cos(sp.pi*x/2), (x, 0, 1))
gbar_v    = sp.integrate(varsigma_test*sp.cos(sp.pi*x/2), (x, 0, 1))
```
both at `Pi -> Pi_star`, then
```
expect_zero("g_eps affine identity", sp.simplify(gbar_phys - (gminus + eps*(gbar_v - gminus))))
```
and mirror for \(\bar S\). This puts the linearity of integration plus \(\bar g[\Sigma_*] = \mathfrak g_*\) on the test surface instead of pure distributivity.

**Verification:**
After the fix, the saved transcript should show the new identity expression (containing `Integrate[...]` evaluations or explicit numeric values from `Pi_star`) reducing to 0, rather than the current trivially-zero form.

### F2 — insufficient_verification

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage146_positive_deformation_expansion_sympy_audit.py:22-28` (and lines 37-43 numerical only)

**What's wrong:**
The notes' two headline closed-form moments are
\[g(\Pi) = \tfrac{2\Pi(2\Pi e^\Pi + \pi)}{(4\Pi^2+\pi^2)(e^\Pi-1)}, \qquad \mathcal S_q(\Pi) = \tfrac{\Pi(\kappa\tanh\kappa + \Pi(e^{-\Pi}\sech\kappa - 1))}{(1-e^{-\Pi})(\kappa^2-\Pi^2)}\]
with \(\kappa = \pi/2\). The Mathematica script verifies the first symbolically (`expectZero["g(Pi) direct-formula", gDirect - gFormula]` at line 48, where `gDirect = Integrate[sigma*Cos[Pi*x/2], {x, 0, 1}]`). The SymPy script never performs the analogous symbolic check for either \(g(\Pi)\) or \(\mathcal S_q(\Pi)\) — it only tests \(\mathcal S_q\) numerically at three samples (Pi = 1, 3/2, 5/2). Neither engine performs a symbolic check of \(\mathcal S_q\) (no `Integrate[sigma*kq, {x,0,1}]` compared to `sFormula`).

This means the bridge from "definition of the moment as an integral over the canonical exponential source" to "the closed-form expression printed by the script" is anchored only at 3 numeric points (for \(\mathcal S_q\)) in both engines, and only in Mathematica for \(g\).

**Why this matters:**
All subsequent results (the \(\Pi_*\) root, the \(\mathfrak g'_*, \mathcal S'_*\) derivatives, the boxed compensation laws) ride on those two closed forms. If either closed form had a bug (e.g., a stray sign or a missing factor that happens to vanish at the sampled Pi values), the SymPy audit could pass without detecting it.

**Required change:**
Add to the SymPy script, after line 28 and before the numerical kernel loop:
```
expect_zero("g(Pi) direct-formula", sp.simplify(sp.integrate(Sigma*sp.cos(sp.pi*x/2), (x, 0, 1)) - gPi))
expect_zero("S_q(Pi) direct-formula", sp.simplify(sp.integrate(Sigma*Kq, (x, 0, 1)) - Sformula))
```
and add to the Mathematica script the parallel `S_q(Pi) direct-formula` check after line 48:
```
sDirect = FullSimplify[Integrate[sigma*kq, {x, 0, 1}], Assumptions -> p > 0];
expectZero["S_q(Pi) direct-formula", sDirect - sFormula];
```

If SymPy `integrate` cannot resolve one of these symbolically in closed form, fall back to a high-precision numeric comparison at, say, 8 random \(\Pi\) values drawn from [0.5, 3.0] with `mpmath` precision 50 — but do at least one symbolic attempt first.

**Verification:**
The new transcript will contain `g(Pi) direct-formula = 0` (SymPy) and `S_q(Pi) direct-formula = 0` (both engines). If SymPy falls back to numeric, the transcript should show 8 sample points each with `|diff| < 1e-30`.

### F3 — insufficient_verification

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage146_positive_deformation_expansion_mathematica_audit.wl:32`

**What's wrong:**
The Mathematica banner reads `STAGE 129 — FINITE-CORRECTION EXPANSION FOR POSITIVE MOUTH-LAYER DEFORMATIONS` but this is the stage-146 audit script. The saved transcript at line 11 echoes the wrong stage number in the section header. Filename, paper card, output filename, and final `Print["Stage 146 Mathematica audit passed."]` (line 93) all confirm the intended unit is 146, so the `STAGE 129` literal is stale scaffolding from an earlier copy.

This is not a math defect, but in v2 the report flags it because (a) the audit transcript is read downstream and (b) the discrepancy suggests the .wl was duplicated from a stage-129 file without a full review, which raises the prior on `mathematica_transliteration` for other content (though, as documented below, the engine-cross-check section, the engine does perform an extra independent symbolic check (`gDirect`) that SymPy does not — so the file is not a pure transliteration).

**Why this matters:**
If a future reader greps the transcript for "Stage 129" they'll see a result that has nothing to do with stage 129's actual claim.

**Required change:**
Edit line 32 of `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage146_positive_deformation_expansion_mathematica_audit.wl` from
```
banner["STAGE 129 — FINITE-CORRECTION EXPANSION FOR POSITIVE MOUTH-LAYER DEFORMATIONS"];
```
to
```
banner["STAGE 146 — FIRST-ORDER EXPANSION FOR POSITIVE MOUTH-LAYER DEFORMATIONS"];
```

**Verification:**
The refreshed Mathematica transcript should show `STAGE 146 — FIRST-ORDER EXPANSION FOR POSITIVE MOUTH-LAYER DEFORMATIONS` in the second banner block.

## Independent-derivation check (Mathematica)

The .wl is not a pure transliteration of the .py. Evidence:
- The Mathematica script adds a symbolic identity check `expectZero["g(Pi) direct-formula", gDirect - gFormula]` (line 48) that the SymPy script does not have. The SymPy script never integrates `Sigma*Cos[Pi*x/2]` symbolically; the Mathematica one does. This is an independent derivation step beyond the SymPy choreography.
- The Mathematica script adds the explicit `expectApprox["Pi_* compensation point", gStar, N[gMinus, 40], 10^-20]` consistency check on the located root (line 76); SymPy just prints and moves on.
- Both engines use the same algebraic closed forms for `gFormula` and `sFormula`. That coincidence is unavoidable — both are encoding the same identity from the notes — and is not by itself a transliteration tell.

The shared scaffolding (`expect_zero` / `expectZero`, same variable names, same banner phrasing) is the project's standard idiom, not an engine-specific port. No `mathematica_transliteration` finding.

## Engine cross-check

Both engines print numerically agreeing values:

| Quantity | SymPy | Mathematica |
|---|---|---|
| \(\Pi_*\) | 1.50882951349315552747043511772 | 1.5088295134931555830055507559542749... |
| \(\mathfrak g_*\) | 0.758035078944662822951954793784 | 0.7580350789446628269196808904141104... |
| \(\mathcal S_*\) | 0.658075937605429271930315313437 | 0.6580759376054292746166018491599488... |
| \(\mathfrak g'_*\) | 0.0714453558083195219958997177670 | 0.0714453558083195211894603019881357... |
| \(\mathcal S'_*\) | 0.0483709542125040994477978749511 | 0.0483709542125040992653572761724701... |

Agreement to ~12 digits on \(\Pi_*\) and 9-11 digits on the derivatives (SymPy's `nsolve` was called with starting guess 1.5 and 30-digit target precision; Mathematica's `FindRoot` ran at WorkingPrecision -> 80 with target 30, then printed N[..., 30]). The mild trailing-digit drift is consistent with the different precision strategies and not a finding.

Both engines also agree on the boxed values from the notes (notes line 14-16: \(\Pi_*\approx 1.50882951349316\), \(\mathfrak g_*\approx 0.758035078944663\), \(\mathcal S_*\approx 0.658075937605429\); notes line 122-128: \(\mathfrak g'_*\approx 0.0714453558083195\), \(\mathcal S'_*\approx 0.0483709542125041\)). No `engine_disagreement` finding.

## Verdict justification

The unit is not clean. Three script-side findings stand:
1. The only `expect_zero` assertions in either engine (the affine laws) are algebraic identities of the script's own definitions and cannot fail. This is the headline assertion the saved transcripts proudly print "PASS" for, so it materially overstates how much the scripts test.
2. The substantive analytic step — that the closed-form moments match the integrals of the canonical exponential source against the kernels \(c(x)\) and \(K_q(x)\) — is symbolically tested in Mathematica only for \(g(\Pi)\). \(\mathcal S_q(\Pi)\) is tested only at three numeric samples in both engines, and \(g(\Pi)\) is not symbolically checked in SymPy at all. The closed forms are the load-bearing input for every downstream value in the stage.
3. The Mathematica banner names the wrong stage. Cosmetic, but it's a script-side defect.

Paper alignment is partial: the closed-form moments, the numeric \(\Pi_*\)/derivatives, and the affine-law form all match the notes' boxed expressions. The compensation laws \(\delta\Pi\) and \(\delta\mathcal S_q\) are constructed in the scripts in the form the notes give, but as derived expressions rather than tested identities — that is acceptable because the formulas are definitions given the linearization, not independent claims. No `paper_misalignment` finding.

Attacks that did not produce a finding: I checked that `Pi` (the symbol) is positive-real declared in both engines (sympy line 16, math line 35) and that the `Pi` symbol does not collide with `sp.pi` / `Pi` (Mathematica's built-in) in problematic ways — both scripts use `Pi`/`p` as the bias parameter and `sp.pi`/`Pi` as the mathematical constant; SymPy uses different names (`Pi` user symbol vs. `sp.pi`) so the collision is hidden but the symbol `Pi` shadows the `sympy.pi` constant only locally and `sp.pi` is used throughout the formulas, so no domain error. I checked the F1 lower-branch input `rF1 = sqrt(12*(37/20)^2/pi^2 - 1)`: with \((37/20)^2 = 1.8769\cdot 12 / \pi^2 \approx 2.2807\), \(rF1\approx 1.131\), \(\sqrt{1+rF1^2}/2\approx 0.7522\), so \(\mathfrak g_- = rF1 - \sqrt{1+rF1^2}/2 \approx 0.378\) — wait, the script outputs `g_* = 0.758035...` not 0.378. Let me re-check: actually rF1 - sqrt(1+rF1^2)/2 with rF1=1.131 gives 1.131 - 0.752 = 0.379, but the printed g_* is 0.758, twice that. This made me look more carefully: the F1 inputs are conventions from upstream and the published value 0.758 matches the notes (boxed \(\mathfrak g_*\approx 0.758035078944663\)). On rechecking the formula `gminus = rF1 - sqrt(1+rF1^2)/2`: with rF1^2 = 12*(37/20)^2/pi^2 - 1 = 12*1.8769/9.8696 - 1 = 1.2812, so rF1 = 1.1319, sqrt(1+rF1^2) = sqrt(2.2812) = 1.5103, /2 = 0.7551, so rF1 - sqrt(...)/2 = 1.1319 - 0.7551 = 0.3767. But Python evaluates `gminus = simplify(rF1 - sqrt(1+rF1^2)/2)` and `nsolve(gPi - gminus, 1.5)` — and the *script's printed* `g_* = gPi.subs(Pi, Pi_star)` is `0.758035...`. That is `gPi` evaluated at the root, not `gminus`. So either `gminus = 0.758...` symbolically simplifies to a different form, or my arithmetic is wrong. Re-doing: `rF1^2 = 12*1.8769/9.8696 - 1`. Actually 12*1.8769 = 22.5228, /9.8696 = 2.2820, -1 = 1.2820 — so rF1^2 = 1.2820, rF1 = 1.1322. sqrt(1+1.2820)/2 = sqrt(2.2820)/2 = 1.5106/2 = 0.7553. 1.1322 - 0.7553 = 0.3769. So gminus ≈ 0.3769. But g_* should equal gminus (because Pi_* is found by `nsolve(gPi - gminus, ...)`). And the script prints g_* = 0.758035. Discrepancy by factor of ~2! Let me check whether the formula in `gminus` matches the notes — I cannot find an explicit formula for `gminus` in the notes provided, only the numeric value `mathfrak g_* ≈ 0.758035`. So either (a) my arithmetic is off, or (b) there is a real bug. Re-doing more carefully: rF1 = sqrt(12*(37/20)**2 / pi**2 - 1). (37/20) = 1.85, squared = 3.4225. 12 * 3.4225 = 41.07. /pi^2 = 41.07/9.8696 = 4.1612. -1 = 3.1612. sqrt(3.1612) = 1.7780. So rF1 ≈ 1.778. sqrt(1+rF1^2) = sqrt(4.1612) = 2.0399, /2 = 1.0200. rF1 - sqrt(...)/2 = 1.778 - 1.020 = 0.758. That matches! I had (37/20)^2 wrong as 1.8769 — actually 37/20 = 1.85 so (37/20)^2 = 3.4225. So the formula is correct and gminus = 0.758035 matches the notes. No finding.

Other attacks tried: looked for sign errors in `sech(kap) - 1` vs `1 - sech(kap)` (notes formula matches script line 26-28: `Pi*(kap*tanh(kap) + Pi*(e^{-Pi}*sech(kap) - 1))/((1-e^{-Pi})*(kap^2-Pi^2))`); checked that Pi_star is on the lower branch by checking dPi sign convention (script uses `nsolve(..., 1.5)` to bias toward the Family-1 lower-branch root); confirmed the kernel `Kq = cosh(kap*(1-x))/cosh(kap)` matches the notes' \(K_q(x)\) definition. All consistent.

## Self-test notes

- Variable independence: my proposed F2 fix uses `sp.integrate(Sigma*Kq, (x,0,1))` where `Sigma` depends on `Pi` and `x` and `Kq` depends on `x` — the integration variable `x` is in both factors, so the integral is a function of `Pi` alone, comparable to `Sformula(Pi)`. No spurious independence.
- Symmetry/parity: my F1 fix uses `6*x*(1-x)` as a test positive normalized profile on [0,1]. \(\int_0^1 6x(1-x)dx = 1\). Together with `c(x) = cos(pi x/2)` (positive on [0,1]) the integrals are well-defined positive reals; no parity-zero traps.
- Trivial-case pre-check: for F1's proposed `gbar_phys - (gminus + eps*(gbar_v - gminus))` at \(\epsilon=0\), `gbar_phys -> integral over Sigma_* * c = gminus` by construction, RHS -> `gminus`, residual -> 0. At \(\epsilon=1\), `gbar_phys -> integral over varsigma_test * c = gbar_v`, RHS -> `gbar_v`, residual -> 0. So the assertion can fail only if the convex-combination structure breaks, which is exactly the substantive content we want to test.
- Path specifications: no `missing_verification_script` finding; both engines are present.
- Paper round-trip: F1 and F2 add new substantive checks; F3 only relabels a banner; none introduces a new constant or a constant that disagrees with the paper.
