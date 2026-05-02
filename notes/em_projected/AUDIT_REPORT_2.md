# Critical Referee Audit — Second Round

**Audited:** all 18 `.md` ↔ `_sympy.py` pairs in `/home/trevnorris/Downloads/em_projected/` after the author's revision in response to `AUDIT_REPORT.md`.
**Method:** six independent referee agents, each given (a) the prior audit's findings for their bundle and (b) the relevant priority-fix list, then asked to verify whether each fix was *substantive* or *cosmetic* and whether new issues were introduced.
**Scripts:** all 18 execute and print `STATUS: PASS`. The author added Wigner/Gaunt machinery, several genuine `sp.integrate` calls, and many new `assert_zero` calls.
**Bottom line:** the bundle moved meaningfully — from "11 of 18 vacuous" to roughly "0 of 18 vacuous, 13 weak, 5 partially substantive." Two of the prior audit's ten priority-list items are clearly addressed; four are partially addressed; **four are not addressed** — most importantly Priority #1 (derive the boxed nonlinear Euler–Lagrange equation in step 14 §2 in sympy) and Priority #2 (step 19 ship an actual branch profile, not a tuned-constants packet).

---

## 1. Updated verdict matrix

| Step | Prior | Updated | Direction | Headline reason |
|---|---|---|---|---|
| 01 | WEAK | WEAK | flat | Genuine √2 ratio still the only substantive content; Faraday and IBP "checks" remain Leibniz-rule tautologies. |
| 02 | VACUOUS | WEAK | up | New §5 has 4 real integration checks; §1–§4 still print-and-delete with implicit decay. |
| 03 | WEAK / VACUOUS | WEAK | up (cosmetic) | Faraday/Ampere asserts still `x − x = 0` after substitution; **`Zg` declared and abandoned** — Maxwell content still uninstantiated. |
| 04 | SUBSTANTIVE / WEAK | SUBSTANTIVE / WEAK | up (real) | **Case B δ(0) regularization fixed substantively**; new §1 "factorization" asserts are anti-pattern §2.9. |
| 05 | WEAK + VACUOUS core | WEAK | up (cosmetic) | The flagged `simplify(ξ·I_WZ/I_WZ)` line **is still in the script verbatim**; new "generic-kernel" block re-encodes `x − x = 0` with `I_WH_when_HZ = sp.Integral(W*Z, …)` typed identically to `I_WZ_generic`. |
| 07 | WEAK + VACUOUS | WEAK→SUBSTANTIVE | up (real) | **Priority #8 fixed**: `sp.series` on `erfc` form cross-checked against hand-typed series and against direct Taylor integration. §3 IBP now real. §4 `H=Z` still typed-twice. |
| 08 | WEAK | WEAK + 1 substantive theorem | up | Y₂₀ lane signature now derived from `gaunt(2,2,2,…)`; δP₀ check real; one z₀-independence assertion remains vacuous. |
| 09 | VACUOUS | WEAK | up (real) | 0 → 15 assertions; `Ξ₁ = n/N + z/D` now derived; **2 new `solve→subst→0` tautologies introduced** at lines 112–113. |
| 10 | VACUOUS | WEAK | up (mostly cosmetic) | 0 → 11 assertions, but **8 are §2.9 anti-pattern** (independence-from-absent-symbols). |
| 11 | WEAK | WEAK | up (one real test) | One new universality test with `μ₁ = 2` (sibling gamma kernel); rest of script unchanged, anti-patterns intact. |
| 12 | VACUOUS | WEAK | up (cosmetic) | 0 → 13 assertions, but **10 of 13 are §2.1, §2.2, or §2.9 anti-patterns**; 3 are mild coefficient checks. Compliance without comprehension. |
| 13 | VACUOUS (1 line) | WEAK | up | Y₂₀ lane derived via `gaunt`; **but a hand-applied `−` sign on the m=1 lane** with no symbolic provenance. IBP still typed. Wall-only no-go still `det = 1/27 ≠ 0`. |
| 14 | VACUOUS | WEAK | up (barely) | 0 → 5 assertions. **Priority #1 NOT addressed**: the boxed nonlinear EL in §2 of the md still has zero computational support. |
| 15 | VACUOUS | WEAK | up | 8 new assertions; Y₂₀ lane via `gaunt`. Same hand-flipped sign as step 13. md addition oversells what the script does. |
| 16 | VACUOUS | SUBSTANTIVE-LITE | up (real) | New algebraic-identity check + Jacobian-determinant assertion + real Gaussian wall integration (`M_Σ = √π`, `K_Σ = 3√π/2`). |
| 17 | WEAK / borderline VACUOUS | WEAK | up | 0 → 11 assertions; **N₄ form mismatch resolved** in both md and py. Branch-sign "test" still tautological. |
| 18 | WEAK | SUBSTANTIVE-LITE | up (real) | Wigner/Gaunt triple overlap genuinely computed; sign on m=±1 lane still smuggled as `−gaunt(...)`. |
| 19 | VACUOUS | WEAK | up (partial) | One real `sp.integrate` of a Gaussian β₂(w); but `R_0(w)` absent, μ_η/T_w/T_Ω/K_η reduced to *constants*, B/Z hand-tuned to satisfy the one-pole target. **Priority #2 not met.** |

**Substantive : Substantive-lite : Weak : Vacuous = 1 : 4 : 13 : 0** (vs. prior 1 : 0 : 6 : 11).

The "no scripts have any assertions" zone is gone. But the WEAK majority is dominated by the same anti-patterns wrapped in new `assert_zero` calls.

---

## 2. Priority-list scorecard

The prior audit listed 10 prioritized fixes. Here is what happened to each:

| # | Prior priority | Status | Notes |
|---|---|---|---|
| 1 | **Step 14 — derive the boxed nonlinear EL in §2 by `sp.diff` + IBP in sympy** | **NOT addressed** | The full nonlinear Euler–Lagrange equation remains hand-typed in the md with zero computational support. Only the eps-quadratic check exists. |
| 2 | **Step 19 — ship an actual branch profile (`R_0`, `β_2`, `μ_η`, `T_w`, `T_Ω`, `K_η`)** | **Partially addressed** | `β_2(w) = e^{−w²/2}` is a real `w`-profile and `M_Σ, K_Σ` are integrated. But `R_0(w)` is absent, the other four wall coefficients are reduced to *constants* (`1, 1, 1/6, 0`), and `B_n, Z_n` are hand-tuned to satisfy the one-pole target. The "branch" is a tuned-constants packet built *from* the targets and then asserted against them. |
| 3 | **Step 17 — resolve N₄ form mismatch and test branch sign** | Half-addressed | N₄ mismatch fully resolved (md + py both make the equivalence explicit). Branch sign "test" added but is itself a tautology (asserts that both ± roots solve the one-pole equation). |
| 4 | **Step 12 — add `diff`-based dependency assertions parallel to step 11 lines 67–72** | Followed literally — and that is the problem | Step 11 lines 67–72 were already the §2.9 "independence from absent symbols" anti-pattern. The author obediently mimicked the *form* and reproduced the anti-pattern at scale. 10 of step 12's 13 new assertions are §2.1, §2.2, or §2.9 anti-patterns. |
| 5 | **Steps 02, 03 — convert prints to asserts; collapse IBP boundary terms; exhibit one concrete leakage term** | Partially addressed | Asserts added (§2.4 fixed). Boundary terms still print-and-deleted. **`Zg` declared in step 03 §5 and abandoned** — the `Leak_ν = −Proj_{W'}[Z F^{wν}]` instantiation requested by the prior audit is *not delivered*. |
| 6 | **Steps 13, 15 — compute Y₂₀·Y₂ₘ² triple overlap symbolically** | Half-addressed | `gaunt(2,2,2,…)` is now invoked. **But the m=1 sign is hand-flipped** (`−gaunt(2,2,2,0,1,−1)`), with no symbolic derivation of the `(−1)^m` factor that converts complex-Y triple overlaps to real-Y squared overlaps. Same in step 18. |
| 7 | **Steps 09, 10 — add asserts to bridge identities** | Addressed in form, mixed in substance | Step 09: 0 → 15 asserts; `Ξ₁` genuinely derived; but two new `solve→subst→0` tautologies introduced at lines 112–113. Step 10: 0 → 11 asserts, but 8 are §2.9 (independence from absent symbols). |
| 8 | **Step 07 — replace hand-typed mouth ℓ-series with `sp.series` against `erfc` closed form** | **Substantively addressed** | Lines 215–217 call `sp.series` on the typed `erfc` form; line 222 asserts equality with the hand-typed series; line 223 cross-checks against direct Taylor integration. The strongest single fix in the bundle. |
| 9 | **Step 11 — extend Taylor universality test to a kernel with μ₁ ≠ 1** | Addressed minimally | One new test against `W₂ = u·e^{−u}` (μ₁ = 2). Sibling kernel from the same gamma family. A symmetric kernel (μ₁ = 0) is still untested. |
| 10 | **Step 04 — derive `I_WZ ∂f = μ₀ I_WS j` symbolically; replace δ(0) Case B with smoothed source** | Half-addressed | **Case B δ(0)-by-fiat is fixed** with a Gaussian regulator and proper `sp.limit`. The reduced equation itself is still consumed as input — new factorization asserts at lines 65–72 verify only that a `w`-independent factor pulls out of `∫ W(w)·factor·... dw` (a §2.9 instance dressed as a derivation). |

**Score:** 2/10 substantive (Priorities #8, #10b). 4/10 half-addressed (#2, #3, #5, #6, #7, #9 — six total, but several with significant carve-outs). **2/10 not addressed (#1, #4 — though #4 was technically followed).** The four non-addressed or anti-fixed items are exactly the ones with the highest physics content.

---

## 3. The pattern: compliance without comprehension

Across the six referee reports, the same diagnosis recurs: **the author moved the bar from "Python finished" to "sympy assertions held," but a large fraction of the new assertions verify only that typed-in expressions are internally consistent with their own typed-in IBP / substitution / derivation steps.** The most concrete examples:

- **Step 12 line 80:** `assert_zero("z0 derivative map", z0d − (Δ·Q_x − Q·D_x)/Δ²)`. Both sides are the same expression after a `q1 → μ₁·Q_x, d1 → μ₁·D_x, /μ₁` rename. Sympy is verifying a rename.
- **Step 14 line 67:** `assert_zero("l=2 angular stiffness", (K_eta + ℓ(ℓ+1)·T_Ω₀).subs(ℓ, 2) − (K_eta + 6·T_Ω₀))`. This asserts that 2·3 = 6.
- **Step 10 lines 89–92:** 8 assertions that `z_0` and `n_0` don't depend on `s_1, h_1, p_1, g_1`, where the typed primitive formulas `Z_0 = Q/Δ` and `N_0 = P²/Δ²` make this true *by definition*.
- **Step 12 line 99:** `assert_zero("Z2 algebraic alias", (Δ·H − Q·S₂) + Δ²·((Q·S₂ − H·Δ)/Δ²))`. Reduces to `A + (−A) = 0`.
- **Step 17 lines 64–67:** "branch sign test" asserts that both ± roots of the one-pole numerator equation satisfy that equation. By definition of `solve`.
- **Step 18 line 39:** `assert_zero("Y20 overlap normalization", overlap0/overlap0 − 1)`. The `x/x = 1` antipattern in pure form.

Wrapping a tautology in `assert_zero(label, ...)` does not make it a check. A regression that flips a sign in any of the typed-in formulas would, in many of these cases, still pass `STATUS: PASS`.

---

## 4. The Wigner/Gaunt sign-flip — pattern across steps 13, 15, 18

The single most-used new piece of physics machinery in this revision is `gaunt(2,2,2,…)`. Three scripts now invoke it to derive the lane signature `λ = (1, 1/2, −1)`. **All three apply a hand sign-flip on the m=1 lane:**

```python
# step_13 line 58, step_15 line 26, step_18 lines ~37
overlap1 = sp.simplify(-gaunt(2, 2, 2, 0, 1, -1))   # leading minus by hand
overlap2 = sp.simplify( gaunt(2, 2, 2, 0, 2, -2))   # no minus
```

Raw values:
- `gaunt(2,2,2,0,0,0) = √5/(7√π)`
- `gaunt(2,2,2,0,1,−1) = −√5/(14√π)`  (negative)
- `gaunt(2,2,2,0,2,−2) = −√5/(7√π)`  (negative)

So raw ratios are `(1, −1/2, −1)`. The author flips the m=1 sign by hand to produce `(1, +1/2, −1)`. The physical justification — that converting complex-Y triple overlaps to real-Y squared overlaps introduces a `(−1)^m` factor via `Y*_{2,m} = (−1)^m Y_{2,−m}` — is **never coded**. A reader cannot tell from any of the three scripts whether the leading `−` is a defended convention or a fudge.

The fix is one expression away (`sp.Rational(-1)**m * gaunt(2,2,2,0,m,-m)/gaunt(2,2,2,0,0,0)` in a list comprehension); not making it leaves the most prominent new physics content with one editorial step in the middle of it.

---

## 5. Cross-cutting issues that persist

1. **The boxed "z₀ drops out of δC" theorem still has no real test** (steps 8, 9, 10). In each, the compatibility surface `compat_p` is constructed *without* `z_0` in its free symbols (because `Ptarget` is treated as a free symbol, not as `N_0/D_0`). So `sp.diff(d_compat, z_0) = 0` is true by construction, not by cancellation. The cancellation surface that gives this result physical content is still not exercised.
2. **The wall-only no-go theorem** (steps 13, 15) still reduces to "the 2×2 system `[(−1, 1/9), (2/3, −1/27)]` has determinant 1/27 ≠ 0, so only the trivial solution." The structural concern stands: this is a triviality, not a no-go. Neither md nor py acknowledges this.
3. **IBP boundary terms remain silently dropped** (steps 02, 03, 12, 13, 14). No `boundary` token in the IBP'd quantities; no `assumptions.assume(W → 0 at ±∞)`; no `Q.zero(...)` for asymptotic decay. The IBP is performed as a typed substitution.
4. **The full nonlinear Euler–Lagrange equation in step 14 §2** (the most central physics object in the parent-action bundle) **still has zero computational support**. No `sp.diff(L, R)`, no `sympy.calculus.euler.euler_equations`. The eps-quadratic expansion is the only sympy work; the rest of the boxed EL is hand-typed and printed.
5. **Maxwell content still absent from code** (steps 01–05, 07, 12). No 4+1 bulk `F^{MN}(t,x,y,z,w)` is ever projected to derive the projected law symbolically. `Z(w)` is now declared in step 03 §5 (as `Zg`) — and **abandoned**: the `Leak_ν = −Proj_{W'}[Z F^{wν}]` instantiation called for in Priority #5 is not delivered. The polynomial example `−Pgp(w) = 1` in step 03 line 226 is the same polynomial computation as step 02 §5, relabelled.
6. **`H = Z ⇒ ξ_eff = ξ` for a generic kernel still verified `x/x = 1`-style** (step 05 lines 102, 121–127; step 07 lines 183–195). The hypothesis `H = Z` is applied at the keyboard by typing identical strings on both sides, not by `sp.subs(H_function, Z_function)` on a symbolic `H`. The line `simplify(ξ · I_WZ / I_WZ)` flagged in the prior audit *is still in the script verbatim*.
7. **Step 19's `R_0(w)` is missing.** The step that promised "actual branch export" supplies one `w`-profile (β₂) and reduces every other wall coefficient to a *constant*. The geometric "parent action → wall integrals" derivation that is supposed to motivate the entire steps 13–19 chain is bypassed.

---

## 6. New issues introduced by the revision

The revision did not strictly improve every script. Several "fixes" introduced fresh anti-patterns:

1. **Step 04 lines 65–72** — new "factorization" asserts verify that a `w`-independent factor pulls out of a `w`-integral. Labelled "zero-mode projection factors brane operator" but performs no projection step. This is anti-pattern §2.9 freshly introduced.
2. **Step 09 lines 112–113** — two new `solve→subst→0` tautologies on the constant-prefactor branch. The original step 09 had zero assertions and therefore couldn't have this anti-pattern; the revision introduced it.
3. **Step 10 lines 89–92** — 8 new component-wise `assert sp.diff(z_0, sym) == 0` for `sym ∈ {s_1, h_1, p_1, g_1}` and `n_0` similarly. Anti-pattern §2.9 in its purest form. These didn't exist before.
4. **Step 12 lines 80–99** — 13 new assertions of which 10 are §2.1, §2.2, or §2.9 anti-patterns. The instruction "add diff-based dependency assertions parallel to step 11 lines 67–72" was followed literally; step 11 lines 67–72 were the WEAK pattern in the prior audit, so the revision propagated it.
5. **Step 14 line 67** — new "l=2 angular stiffness" assertion checks that 2·(2+1) = 6 by `subs`. Pure §2.9.
6. **Step 17 lines 64–67** — "branch sign test" asserts that both ± roots solve the one-pole equation. By definition of `solve`. Tautology dressed as a §2.10 fix.
7. **Step 18 line 39** — `assert_zero("Y20 overlap normalization", overlap0/overlap0 − 1)` — the cleanest possible §2.1 instance.

---

## 7. What was genuinely fixed (the wins)

To balance the criticism: there are five clean wins worth flagging.

1. **Step 04 Case B regularization** (lines 102–116): replaced distributionally-undefined `δ(0)` with a Gaussian regulator and a real `sp.limit`. Substantive.
2. **Step 07 mouth ℓ-series cross-check** (lines 215–223): `sp.series` on the `erfc` closed form, asserted against the hand-typed series, AND cross-checked via direct Taylor integration. Defense in depth. The strongest single fix.
3. **Step 17 N₄ form mismatch** (md lines 281–288, py line 59): the factor-of-9 discrepancy from the prior audit is explicitly reconciled in both md and py.
4. **Step 18 Wigner/Gaunt computation** (lines 36–41): `gaunt(2,2,2,0,m,−m)` is genuinely computed for `m = 0, 1, 2`. The "lane signature is hardcoded" critique is materially weakened. (The m=1 sign-flip is the residual concern.)
5. **Step 19 toy export** (lines 65–119): wall integrals `M_Σ = √π`, `K_Σ = 3√π/2` are now real `sp.integrate` outputs, propagated end-to-end through the algebraic packet with consistency assertions. Falls short of the priority-2 ask, but is a real step from VACUOUS toward WEAK.

These five are not cosmetic. They constitute genuine evidence that a regression in one of those specific places would now be caught.

---

## 8. What still needs to happen for a third pass

In rough order of importance:

1. **Derive the boxed nonlinear EL in step 14 §2 in sympy** — Priority #1 from the prior audit, skipped. This is the most central derivation of the parent-action bundle. `sympy.calculus.euler.euler_equations` or hand `∂L/∂R − d/dt(∂L/∂R_t) − …` via `sp.diff` would do it.
2. **Step 19 — supply a real `R_0(w)` and four real `w`-profile wall coefficients.** Currently four of six profiles are constants. Without a non-trivial `μ_η(w), T_w(w), T_Ω(w), K_η(w)`, the wall integrals are sterile.
3. **Step 19 — derive `B_n, Z_n` from a BdG/port computation independent of the targets**, instead of hand-tuning `B_4 = Z_4 = √π` to satisfy the one-pole equation. As written, the script *uses targets to build inputs and then asserts inputs against targets*.
4. **The `(−1)^m` Wigner sign convention** — derive it from real-Y harmonics and apply it as a symbolic `sp.Rational(−1)**m` rather than typing `−gaunt(...)` (steps 13, 15, 18).
5. **The `H = Z ⇒ ξ_eff = ξ` proof** — declare `H = sp.Function("H")(w)`, build `I_WH = sp.Integral(W·H, w)` symbolically, and assert `I_WH.subs(H, Z) − I_WZ = 0`. The substitution should live in sympy, not at the keyboard (steps 05, 07).
6. **Track IBP boundary terms** — either declare an asymptotic-decay assumption symbolically (`sp.Q.zero(W.subs(w, oo))`) and propagate it, or carry the boundary term in the assertion (steps 02, 03, 13, 14).
7. **Exhibit one concrete leakage term** — multiply the abandoned `Zg` into `ZFw1 = Zg · w` in step 03 §5 and confirm `−Proj_{W'}[Zg · w] ≠ 0` for a specific case. The prior Priority #5 sentence "exhibit one concrete leakage term" still has not been delivered.
8. **Resolve the z₀-drops-out theorem properly** — construct `compat_p` with `Ptarget = N_0/D_0` (so `z_0` enters via `D_0 = K − B_0 − Z_0`) and verify that the `z_0` terms cancel through the algebra rather than by absence.
9. **Address the wall-only no-go theorem honestly** — either acknowledge in the md that it reduces to "the 2×2 coefficient matrix has nonzero determinant," or replace it with a genuine no-go statement that captures something physical.
10. **Negative-failure tests** — for each new substantive assertion, mutate one coefficient and confirm the assertion fires. Currently the suite is silent on this; a buggy revision would not reveal itself.

---

## 9. Net evidentiary delta vs. round 1

The bundle gained roughly **40+ new assertions**, of which:
- ~10 are genuinely substantive new tests (Wigner/Gaunt, mouth ℓ-series cross-check, regulator limit, `M_Σ`/`K_Σ` Gaussian integrals, N₄ form equivalence on the one-pole surface, Ξ₁ derivation via series expansion).
- ~10 are mildly substantive (derivative-coefficient extractions, Frechet cross-checks, primitive Xi static algebraic identity).
- ~20+ are anti-pattern §2.1, §2.2, §2.3, or §2.9 — most of them freshly added in this revision.

The bundle's ratio of evidentiary content to assertion count went *up* on a relative scale (every assertion in step 04 §1 factorization, step 10's component-wise independence, step 12's z0/n0 derivative-map, step 14's l=2 angular stiffness, and step 17's branch-sign test could be deleted with no loss of information). On an absolute scale, the bundle now contains 5 genuinely substantive new physics computations that did not exist before. That is the meaningful delta.

---

## 10. Headline restatement

**The author moved the bundle from VACUOUS-dominated to WEAK-dominated**, and made five clean substantive fixes (step 04 Case B, step 07 ℓ-series, step 17 N₄ mismatch, step 18 Wigner/Gaunt, step 19 toy export). **But the four highest-physics-content items on the priority list — the full nonlinear EL in step 14, the actual branch in step 19, the leakage instantiation in step 03, and the proper z₀-cancellation test in steps 8–9 — were skipped, half-skipped, or replaced with tautologies.** The pattern across the new assertions suggests the revision was carried out by mechanically translating prior-audit prose into `assert_zero` calls without re-deriving the underlying physics — exactly the "compliance without comprehension" failure mode that the prior audit's `solve→subst→0` and `x/x = 1` flags were warning against.

A reader who sees `STATUS: PASS` across 18 scripts in this revision should treat the pass as evidence for: (a) Wigner/Gaunt-derived lane signatures (with one editorial sign), (b) one real `sp.series` cross-check on the mouth Gaussian, (c) one regulated δ(0) limit, (d) one Gaussian wall-integral packet, and (e) one resolved N₄ form discrepancy. They should *not* take it as evidence for the boxed nonlinear EL, the projected Maxwell derivation, the actual physical branch, the open-system leakage thesis, or the z₀-drops-out theorem — the bundle's headline claims. Those remain typed assertions verified against typed counterparts.

---

*Audit complete. 6 referee agents, 18 step pairs, 100% script execution. Net evidentiary delta: roughly 5–10 genuinely substantive new computations against ~20+ freshly-introduced anti-pattern assertions. Two of ten priority-list items addressed; four half-addressed; **four either skipped or replaced with tautologies, including the highest-physics-content item (step 14 nonlinear EL).** Next pass should focus there.*
