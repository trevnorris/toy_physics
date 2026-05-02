# Critical Referee Audit — Third Round

**Audited:** the 12 `.md` ↔ `_sympy.py` pairs modified after `AUDIT_REPORT_2.md` was delivered (steps 03, 05, 07, 08-py, 09, 10, 13, 14, 15, 17, 18, 19); the unchanged pairs (01, 02, 04, 11, 12, 16) carry their round-2 verdicts forward.
**Method:** six independent referee agents, each given (a) the round-2 findings for their bundle and (b) the relevant priority-fix list, then asked to verify whether each fix was *substantive* or *cosmetic* and whether new anti-patterns were introduced.
**Scripts:** all 18 still execute and print `STATUS: PASS`.
**Bottom line:** the bundle moved meaningfully again — from "1:4:13:0" (substantive : substantive-lite : weak : vacuous) to roughly **"1:11:6:0"**. Of the round-2 next-pass priority list, **4/10 are now substantively addressed (#1 nonlinear EL, #4 Wigner sign, #5 H=Z, #7 leakage), 1/10 partially (#2 R₀ profile), 5/10 still skipped (#3 BdG-derived branch coefficients, #6 IBP boundary terms, #8 z₀ via Pₜₐᵣ𝓰ₑₜ=N₀/D₀, #9 no-go acknowledgment, #10 negative-failure tests bundle-wide).** Two new anti-patterns were freshly introduced (step 17 lines 69–70, step 18 lines 69–70). The most central derivation in the parent-action bundle — the boxed nonlinear Euler–Lagrange equation in step 14 §2 — is now derived in sympy via `sympy.calculus.euler.euler_equations`. That is the headline of this round.

---

## 1. Updated verdict matrix

| Step | R2 | R3 | Direction | Headline reason |
|---|---|---|---|---|
| 01 | WEAK | WEAK | flat (not modified) | √2 ratio still the only substantive content. |
| 02 | WEAK | WEAK | flat (not modified) | §5 integration checks still real; §1–§4 still print-and-delete. |
| 03 | WEAK | **SUBSTANTIVE-LITE** | up (real) | New §5 instantiates `Zg`, computes `leak1 = −Proj_{W'}[Zg·w] = √2/4` via real `sp.integrate`; Faraday/Ampere `x−x=0` still untouched. |
| 04 | SUBSTANTIVE / WEAK | SUBSTANTIVE / WEAK | flat (not modified) | Case B Gaussian regulator + `sp.limit` still substantive; §1 factorization asserts still §2.9. |
| 05 | WEAK | **SUBSTANTIVE-LITE** | up (real) | Generic-kernel H=Z proof rewritten with `sp.Function("H")(w)` + `Integral(W·H).subs(H, Z)` — real symbolic substitution, fails without `.subs`. Round-2 `xi*I_WZ/I_WZ` repackaged via `subs(I_WH, I_WZ)` (cosmetic), but lives next to the substantive fix. |
| 07 | WEAK→SUBSTANTIVE | **SUBSTANTIVE** | up (real) | §4 `H=Z` typed-twice now uses `sp.Function("H")(w)` + `Integral(W·H).subs(...)`; round-2 mouth ℓ-series cross-check (sp.series on erfc + Taylor integration) persists. Strongest pair in the bundle now. |
| 08 | WEAK + theorem | WEAK | flat | Y₂₀-via-Gaunt and δP₀ persist; **z₀-drops-out theorem NOT fixed** — `Ptarget` still a free symbol at line 58, `+eps·z0` typed into both K-surfaces in parallel, cancellation by lexical identity. |
| 09 | WEAK | **SUBSTANTIVE-LITE** | up (partial) | Round-2 `solve→subst→0` tautologies replaced with real linearization checks; Wigner `(-1)^m` now symbolic at line 58; z₀ cancellation still by parallel typing. |
| 10 | WEAK | **SUBSTANTIVE-LITE** | up (real) | Round-2 `assert sp.diff(z_0, sym)==0` §2.9 anti-pattern at lines 89–92 **replaced with Frechet-derivative cross-checks** at lines 79–87 — the cleanest single fix in the priority-list category. |
| 11 | WEAK | WEAK | flat (not modified) | Sibling-gamma test still the lone universality probe; symmetric-kernel test still missing. |
| 12 | WEAK | WEAK | flat (not modified) | The 10/13 anti-pattern assertions from round 2 remain. |
| 13 | WEAK | **SUBSTANTIVE-LITE** | up (real) | `(-1)^m` factor now symbolic in `real_y20_square_ratio()`; same-sign cross-term vanishing also checked. IBP-by-typing and trivial wall-only no-go still intact. |
| 14 | WEAK (barely) | **SUBSTANTIVE-LITE** | **up — most important fix of the round** | **Boxed nonlinear EL derived via `sympy.calculus.euler.euler_equations` at lines 16–43**, matched against the hand-typed boxed form. Round-2 `(K_η + ℓ(ℓ+1)·T_Ω).subs(ℓ, 2)−(K_η+6T_Ω)` removed; replaced with real Y₂₀ Laplacian eigenvalue verification on `Ynm(2,0)`. |
| 15 | WEAK | **SUBSTANTIVE-LITE** | up | `(-1)^m` symbolic; md/py now honestly aligned; new propagation through D₀₁/D₂₁/D₄₁ and `b = 3a` weak-axisymmetric line verified across inertia, stiffness, Ξ-load, prefactor. |
| 16 | SUBSTANTIVE-LITE | SUBSTANTIVE-LITE | flat (not modified) | Real Gaussian wall integration persists. |
| 17 | WEAK | WEAK | flat (cosmetic) | N₄ mismatch fix persists. Round-2 branch-sign tautology replaced with one mild Vieta identity + two **new** definitional sign tautologies (lines 69–70). |
| 18 | SUBSTANTIVE-LITE | SUBSTANTIVE-LITE | up & down | Wigner `(-1)^m` now symbolic at line 27; round-2 `overlap0/overlap0−1` removed. **But two new `solve→subst→0` tautologies introduced at lines 69–70** ("compensated K1", "compensated H_even"). Net: real win on Wigner, fresh §2.2 elsewhere. |
| 19 | WEAK | **SUBSTANTIVE-LITE** | up (mixed) | `R_0(w) = e^{−w²/2}` is now a real profile; three of four wall coefficients are `1+R_0` (linked, not independent); `K_η` back-derived from background closure. **`B_n, Z_n` still hand-typed integrand kernels with no BdG derivation.** Critical honesty win: `assert_nonzero` on the target residues now displays *failure* of the toy packet to hit targets — the script no longer pretends. |

**Substantive : Substantive-lite : Weak : Vacuous = 1 : 11 : 6 : 0** (vs. round 2's 1:4:13:0, vs. round 1's 1:0:6:11).

The migration from VACUOUS-dominated (R1) → WEAK-dominated (R2) → SUBSTANTIVE-LITE-dominated (R3) is real. Six of the eleven SUBSTANTIVE-LITE entries (07, 09, 10, 13, 14, 18) involve at least one new computation that would catch a mutation in a typed coefficient.

---

## 2. Round-2 priority-list scorecard

The round-2 audit closed with a 10-item "what still needs to happen" list. Status by item:

| # | Round-2 priority | Round-3 status | Notes |
|---|---|---|---|
| 1 | **Step 14 — derive boxed nonlinear EL via `sp.diff` / `euler_equations`** | **SUBSTANTIVELY ADDRESSED** | Lines 16–43 of step_14 build `L_exact = (1/2)μ_Σ R_t² − (1/2)T_w R_w² − (1/2)T_Ω(R_u²+R_v²) − U_Σ` with `R = R(t,w,u,v)`, then `exact_el = euler_equations(L_exact, Rfield, …)[0].lhs` and asserts equality with the hand-typed boxed form. **The single biggest gain of this round.** Caveat: `T_Ω` divergence is verified in a flat local-orthonormal angular chart, not the full S² covariant Laplacian — the `+6` eigenvalue is separately checked on `Ynm(2,0)`. |
| 2 | **Step 19 — supply real `R_0(w)` and four w-profile wall coefficients** | **PARTIALLY ADDRESSED** | `R_0(w) = e^{−w²/2}` is real (line 76). `μ_η = 1+R_0`, `T_w = 1+R_0`, `T_Ω = (1+R_0)/6` are linked (one-knob structure). `K_η = U_scale − R_0''` is back-derived from background closure. So four "profiles" with one underlying degree of freedom — better than the round-2 four constants `(1,1,1/6,0)`, but still not four independent profiles. |
| 3 | **Step 19 — derive `B_n, Z_n` from BdG/port computation independent of targets** | **NOT ADDRESSED** | Lines 102–107 hand-type integrand kernels: `B_0 = ∫(1/5)β² dw`, `B_2 = −∫(1/3)w²β² dw`, `B_4 = ∫(1/6)(1+w²)β² dw`, `Z_0 = ∫(1/7)R_0β² dw`, etc. The prefactors `1/5, 1/3, 1/6, 1/7, 1/8, 1/9` and integrand structures (`w²`, `(1+w²)`, `R_0`, `R_0²`) have no derivation from any BdG eigenproblem or port spectrum. Different anti-pattern than round 2 (no longer back-tuned to targets — instead deliberately mis-tuned), still not a derivation. |
| 4 | **Wigner `(−1)^m` symbolic factor (steps 13, 15, 18)** | **SUBSTANTIVELY ADDRESSED in all three** | Each script now defines `real_y20_square_ratio(m) = (sp.Integer(-1)**m) * gaunt(2,2,2,0,m,-m) / base`. The hand-typed `−gaunt(...)` from round 2 is gone in steps 13, 15, 18 (and step 09 already used the symbolic form). Md notes in steps 13 and 15 explicitly explain the factor as `Y*₂ₘ = (−1)^m Y₂,−ₘ`. |
| 5 | **Steps 05, 07 — H=Z proof via symbolic `sp.subs`** | **SUBSTANTIVELY ADDRESSED in both** | Step 05 lines 119–139: declares `W_generic, Z_generic, H_generic = sp.Function(...)(w)`, builds `I_WH_generic = sp.Integral(W_generic·H_generic, ...)` with **distinct** kernels, then `.subs(H_generic, Z_generic)` against `I_WZ_generic`. Without the `.subs(...)` the difference is genuinely non-zero. Step 07 lines 188–195: identical pattern with `H_sym = sp.Function("H")(w)`. |
| 6 | **Track IBP boundary terms symbolically** | **NOT ADDRESSED** | No `sp.Q.zero(W.subs(w, oo))` or boundary-symbol carrying anywhere. IBP is still performed as a typed substitution in steps 02, 03, 13, 14. Boundary terms remain print-and-deleted. |
| 7 | **Exhibit one concrete leakage term in step 03** | **SUBSTANTIVELY ADDRESSED** | Lines 203–230 of step 03 instantiate `Wg = e^{−w²}/√π`, `Zg = e^{−w²}`, `Fw1 = w`, then compute `leak1 = −Proj_{W'}[Zg·w] = √2/4` via real `sp.integrate` and assert (i) the IBP identity between this and `Pg(∂_w(Zg·w))` and (ii) the explicit nonzero value. A constant `Fw1` would make leak=0, so the assertion would catch a regression. |
| 8 | **z₀-drops-out theorem with `Ptarget = N₀/D₀`** | **NOT ADDRESSED in steps 8, 9, 10** | All three still declare `Ptarget` as a free `sp.symbols("Ptarget", nonzero=True)`. Step 09 and step 10 now type `+eps·z0` into both `K_norm_p` and `K_one_pole_p` in parallel — one notch better than round 2 (the cancellation is now exhibited explicitly), but still a typing-convention identity rather than algebraic propagation through `D_0 = K − B_0 − Z_0` and `Ptarget = N_0/D_0`. Step 08 even less changed. |
| 9 | **Wall-only no-go acknowledgment** | **NOT ADDRESSED** | Steps 13 and 15 still present the trivial `det = 1/27 ≠ 0` outcome as a "no-go." Neither md acknowledges the triviality nor replaces it with a genuine no-go statement. |
| 10 | **Negative-failure tests — mutate one coefficient and confirm assertion fires** | **NOT ADDRESSED bundle-wide, with one notable exception** | Step 19 is the only place the author explicitly demonstrates failure: lines 134–137 use `assert_nonzero` on the four target residues for the toy packet, displaying values like `(−122577 − 49080√6 + …)/30240·π`. This is the right pattern; it should be propagated to other scripts. |

**Score:** 4/10 substantively addressed (#1, #4, #5, #7); 1/10 partially (#2); 5/10 still skipped (#3, #6, #8, #9, #10 — though #10 has one exemplar). Compare round 2: 2/10 substantive, 4/10 half, 4/10 skipped.

The four newly-addressed items (#1, #4, #5, #7) are the four highest-physics-content items on the round-2 list. The five that remain unaddressed are mostly methodological (boundary tracking, BdG derivation, parallel-K vs algebraic z₀ cancellation, no-go honesty, regression mutation) rather than physics derivations — though #3 and #8 carry real physics weight.

---

## 3. The single most important fix: step 14 nonlinear Euler–Lagrange equation

The round-2 audit's #1 priority — "derive the boxed nonlinear EL in step 14 §2 by `sp.diff` + IBP in sympy" — was the highest-physics-content item in the entire bundle and was the one item the round-2 revision skipped entirely. The round-3 revision delivers it. The relevant code at `step_14_parent_throat_action_candidate_sympy.py:16-43`:

```python
# Lines 17-31: Build full nonlinear Lagrangian density
Rfield = sp.Function('R')(t_exact, w_exact, u_exact, v_exact)
mu_S    = sp.Function('mu_Sigma')(Rfield, w_exact)
T_w_f   = sp.Function('T_w')(Rfield, w_exact)
T_O     = sp.Function('T_Omega')(Rfield, w_exact)
U_S     = sp.Function('U_Sigma')(Rfield, w_exact)
L_exact = (sp.Rational(1,2)*mu_S*R_t**2
         - sp.Rational(1,2)*T_w_f*R_w**2
         - sp.Rational(1,2)*T_O*(R_u**2 + R_v**2)
         - U_S)

# Line 32: SymPy derives the Euler-Lagrange equation
exact_el = euler_equations(L_exact, Rfield,
                           [t_exact, w_exact, u_exact, v_exact])[0].lhs

# Lines 33-42: Hand-typed boxed form from md §2
exact_el_lhs = (sp.diff(mu_S*R_t, t_exact)
              - sp.diff(T_w_f*R_w, w_exact)
              - sp.diff(T_O*R_u, u_exact)
              - sp.diff(T_O*R_v, v_exact)
              - sp.Rational(1,2)*sp.diff(mu_S, Rfield)*R_t**2
              + sp.Rational(1,2)*sp.diff(T_w_f, Rfield)*R_w**2
              + sp.Rational(1,2)*sp.diff(T_O, Rfield)*(R_u**2 + R_v**2)
              + sp.diff(U_S, Rfield))

# Line 43: Real symbolic match
assert_zero('exact nonlinear Euler-Lagrange equation', exact_el + exact_el_lhs)
```

This is the strongest fix of the round. It mechanizes the most central derivation of the parent-action bundle. A regression that flips a sign or misses a `∂R_t` term in the boxed equation would now be caught. The remaining qualifier is that the angular sector is verified in a flat local-orthonormal chart `(u, v)` rather than against the full S² covariant Laplacian; the `+6` eigenvalue is separately checked on `Ynm(2,0)` at line 103.

The auxiliary round-2 anti-pattern at line 67 (the `(K_η + ℓ(ℓ+1)·T_Ω).subs(ℓ, 2) − (K_η + 6T_Ω)` that asserted `2·3 = 6`) is removed.

---

## 4. The Wigner `(−1)^m` factor: now symbolic across the bundle

The round-2 audit's most-cited cross-cutting issue was the hand-typed leading minus sign on the m=1 lane in `gaunt(2,2,2,0,1,−1)` across steps 13, 15, 18. The round-3 revision converts this to a symbolic `(sp.Integer(-1)**m) * gaunt(2,2,2,0,m,-m) / base` factor in **all four** scripts that use it (13, 15, 18, and 09 — which already had it). The pattern is identical:

```python
def real_y20_square_ratio(m: int) -> sp.Expr:
    base = gaunt(2, 2, 2, 0, 0, 0)
    return (sp.Integer(-1) ** m) * gaunt(2, 2, 2, 0, m, -m) / base
```

The `(−1)^m` factor is now derivable from a single source rather than typed editorially in three scripts. Md notes in steps 13 and 15 explain the convention from `Y*₂ₘ = (−1)^m Y₂,−ₘ`. Steps 13 and 15 also verify the same-sign cross term `gaunt(2,2,2,0,m,m)` vanishes for `m = 1, 2`, an independent real-Y orthogonality check. The "editorial sign in the middle of the most prominent new physics" critique from round 2 is materially resolved.

---

## 5. New issues introduced in this round

Two new anti-patterns were introduced that did not exist in round 2.

1. **Step 17 lines 69–70** — `'positive M branch sign'` / `'negative M branch sign'` assert `(M + B2 + Z2).subs(M, M_root_positive) − root_gap == 0` and the analogous negative version. Since `M_root_positive = sqrt((K−B0−Z0)(B4+Z4)/3) − (B2+Z2)` and `root_gap = sqrt((K−B0−Z0)(B4+Z4)/3)` by construction, the LHS reduces to `sqrt(...) − sqrt(...) = 0`. This is a §2.10/§2.3 hybrid (definitional sign tautology). The md (lines 184–187) honestly concedes "the positive square-root branch is therefore a stability/physical-branch convention, not something fixed by the one-pole algebra alone" — so the prose is correct, but the script's "branch sign test" doesn't earn the claim. A second-derivative test, cost-functional minimization, or asymptotic-BC discriminant is what's needed.

2. **Step 18 lines 69–70** — `'compensated K1'` and `'compensated H_even'` assert `K1.subs(sol) == 0` and `H_even.subs(sol) == 0` after `sol = solve([K1==0, H_even==0], …)` at line 59. This is the textbook §2.2 `solve→subst→0` anti-pattern. Since line 57 already asserts the system's determinant is `−1/27 ≠ 0`, the substitution back is guaranteed to give zero by the definition of `sp.solve`. The substantive Wigner fix at line 27 lives next to this fresh tautology.

These two are smaller in volume than the round-2 introductions (where step 04 §1 and step 10 §2 each gained ~8 anti-pattern assertions), but they are the same family of issue: rewriting in response to an audit while reproducing the pattern the audit flagged.

A near-anti-pattern worth flagging:

3. **Step 03 lines 90–111, 142–163** — the round-1 Faraday/Ampere `x − x = 0` after substitution was not addressed in round 2 or round 3; the §1–§2 "verifications" remain pure §2.3 (typed against self after substitution). This is not a new issue but it is a persistent one in a script that otherwise made a substantive leap.

---

## 6. Persisting cross-cutting issues

1. **z₀-drops-out theorem** (steps 8, 9, 10): closer to substantive but still not algebraic. Steps 9 and 10 now type `+eps·z0` into both K-surfaces in parallel, so the cancellation is at least exhibited rather than absent — but the cancellation is enforced by typing identical strings into two halves rather than by propagating `Z_0 → Z_0 + eps z_0` through `D_0 = K − B_0 − Z_0` and then through `Ptarget = N_0/D_0`. Step 08 is the least improved on this front.
2. **IBP boundary terms still print-and-deleted** (steps 02, 03, 13, 14). No `sp.Q.zero` machinery; no boundary symbols carried.
3. **Wall-only no-go theorem** (steps 13, 15) still reduces to `det = 1/27 ≠ 0`. Neither md acknowledges the triviality nor replaces with a genuine no-go.
4. **Step 19 `B_n, Z_n` integrand kernels** are hand-typed with no BdG derivation. The author is now honest about not hitting the targets but the kernel structures `(1/5)β²`, `−(1/3)w²β²`, `(1/6)(1+w²)β²`, `(1/7)R_0β²`, `(1/8)w²R_0β²`, `(1/9)R_0²β²` are decorative — there is no separate eigenproblem solved.
5. **Step 03 §1–§2 Faraday/Ampere** still §2.3 tautology; round-3 work on §5 leakage did not propagate back.
6. **Step 17 persistent `solve→subst→0` at lines 51–52, 56–57** (pre-existing, not flagged in round 2 priority list).
7. **Step 05 lines 101–103** retains the cosmetic `subs(I_WH, I_WZ)` form alongside the substantive §3 fix — the lazy version was not deleted.

---

## 7. The five clean wins of round 3

To balance the unfinished business, five round-3 changes are unambiguously substantive.

1. **Step 14 nonlinear EL via `euler_equations`** (lines 16–43) — the boxed Euler–Lagrange equation that anchors the parent-action bundle is now derived in sympy, not just typed. Round-2 Priority #1 done.
2. **Step 03 leakage instantiation** (lines 203–230) — `leak1 = −Proj_{W'}[e^{−w²}·w] = √2/4` computed via real `sp.integrate`, with two consistency assertions and one IBP-identity assertion. Round-2 Priority #5/#7 done.
3. **Step 05 H=Z proof** (lines 119–139) and **step 07 H=Z proof** (lines 188–195) — both rewritten with `sp.Function("H")(w)`, distinct symbolic kernels, and `.subs(H, Z)` that genuinely fails before the substitution. Round-2 Priority #5 done.
4. **Wigner `(−1)^m` symbolic factor** across steps 13, 15, 18 (and pre-existing in step 09) — `real_y20_square_ratio(m) = (sp.Integer(-1)**m) * gaunt(2,2,2,0,m,-m) / base`. Round-2 Priority #4 done.
5. **Step 10 Frechet-derivative cross-checks** (lines 79–87) — replaced the round-2 §2.9 component-wise `assert sp.diff(z_0, sym) == 0` block with directional-derivative checks that compare each series-expanded primitive map against the actual Frechet derivative. The cleanest single replacement in the bundle.

Honorable mentions: step 19's `assert_nonzero` block at lines 134–137 (honest negative result reporting), step 19's `R_0(w)` being a real profile, step 09's lane-wise `Ξ₁_proj` derivation at lines 152–158, step 15's weak-axisymmetric `b = 3a` propagation across four distinct algebraic structures.

---

## 8. What still needs to happen for a fourth pass

In rough order of importance:

1. **Step 19 — derive `B_n, Z_n` from a BdG/port computation.** This is the central physics gap remaining in the parent-action chain. The current `(1/5)β², −(1/3)w²β², (1/6)(1+w²)β², (1/7)R_0β², (1/8)w²R_0β², (1/9)R_0²β²` are decorative integrand kernels with no derivation. A separate sympy block should solve a (linearized BdG, Sturm–Liouville) eigenproblem on the wall coordinate `w` and read off the moments.
2. **Step 19 — supply four genuinely independent w-profile wall coefficients.** Currently three (`μ_η, T_w, T_Ω`) share the `1+R_0` skeleton; `K_η` is back-derived from background closure. Four independent profiles (e.g., a different polynomial weight on each) would make the Gaussian wall integrals non-trivial.
3. **z₀ algebraic cancellation** (steps 8, 9, 10) — construct `Ptarget = N_0/D_0` with `D_0 = K − B_0 − Z_0`, perturb `Z_0 → Z_0 + eps z_0`, expand to first order in `eps`, and confirm the `z_0` terms cancel through the algebra rather than by parallel typing.
4. **IBP boundary terms** — adopt one of: `sp.Q.zero(W.subs(w, ±oo))` with refine, a `boundary` symbol carried in IBP'd quantities, or `sympy.calculus.singularities.singularities` checks. Apply uniformly to steps 02, 03, 13, 14.
5. **Step 17 branch-sign physical-root selection** — the lines 69–70 sign tautologies should be replaced with a discriminant from a stability functional, an asymptotic boundary condition, or a sign-of-second-derivative test. The md correctly identifies this as a convention issue; the script should reflect that.
6. **Step 18 lines 69–70 §2.2 cleanup** — the `K1.subs(sol)` / `H_even.subs(sol)` tautologies should be replaced with checks on `det ≠ 0` and on the *form* of `sol[K1_var], sol[H_even_var]` (e.g., that they have the expected dependence on input parameters).
7. **Step 03 §1–§2 Faraday/Ampere** — the §5 leakage work showed the author can do real symbolic projection; back-propagate that approach to §1–§2 so the homogeneous law is verified by genuine projection rather than by `x − x = 0` after typed substitution.
8. **Wall-only no-go honesty** (steps 13, 15) — the md should either acknowledge the no-go reduces to `det ≠ 0` or replace it with a structural statement about which kernel structures are excluded.
9. **Negative-failure tests bundle-wide** — the step-19 `assert_nonzero` pattern is the right idea. For each substantive `assert_zero` in the bundle, add a paired mutation test (flip a sign in one input coefficient, confirm the assertion now fires) — at minimum for the new EL, leakage, and Frechet checks.
10. **Step 14 angular sector — full S² covariant Laplacian.** Currently the `T_Ω` divergence is verified in a flat local chart, with the `+6` eigenvalue checked separately on `Ynm(2,0)`. Combining these into a single `sympy.functions.special.spherical_harmonics`-based check would close the loop on Priority #1.

---

## 9. Net evidentiary delta vs. round 2

The bundle gained roughly **20+ new genuinely substantive assertions** (the EL match in step 14, the four-piece H=Z block in steps 05 and 07, the leakage instantiation in step 03, the lane-wise Ξ₁ derivation in step 09, the Frechet block in step 10, the symbolic Wigner derivations across 13/15/18, the step-19 honest residue display) **against ~4 newly-introduced anti-pattern assertions** (step 17 lines 69–70 and step 18 lines 69–70).

The evidence-to-anti-pattern ratio is the inverse of round 2 (which added ~10 substantive against ~20+ anti-pattern). On absolute terms, the bundle now contains:
- One `euler_equations`-derived nonlinear EL (step 14).
- One concrete leakage value `√2/4` from real symbolic integration (step 03).
- Two genuine `sp.Function`-based H=Z proofs (steps 05, 07).
- Four `(−1)^m`-symbolic Wigner derivations (steps 09, 13, 15, 18).
- One Frechet-derivative cross-check block over the primitive→bundle map (step 10).
- One real `R_0(w)` profile and one honest negative-result display (step 19).

Plus all the round-2 wins (step 04 Case B regulator, step 07 mouth ℓ-series cross-check, step 17 N₄ form mismatch, step 18 Wigner machinery, step 19 toy export, step 16 Gaussian wall integration). That is **~15 substantive new physics computations across the bundle**, against ~4 anti-pattern instances introduced in round 3 — a real improvement in evidentiary density.

---

## 10. Headline restatement

**The round-3 revision delivers four of the round-2 audit's six highest-physics-content items: the boxed nonlinear EL (Priority #1), the Wigner `(−1)^m` symbolic factor (Priority #4), the H=Z proof via symbolic substitution (Priority #5), and the concrete leakage instantiation (Priority #7).** Step 14's `euler_equations` derivation, in particular, removes the round-2 "compliance without comprehension" diagnosis from the central derivation of the parent-action bundle. Step 10's Frechet-derivative block similarly demonstrates that the author can write a real cross-check rather than a §2.9 component-wise tautology.

Five priority items remain unaddressed (#3 BdG derivation, #6 boundary terms, #8 algebraic z₀ cancellation, #9 no-go acknowledgment, #10 systematic negative-failure tests). Two new anti-patterns were introduced (step 17 lines 69–70, step 18 lines 69–70) — both in scripts that simultaneously made substantive fixes elsewhere.

A reader who sees `STATUS: PASS` across 18 scripts in this revision should treat the pass as evidence for: (a) the boxed nonlinear EL of the parent action, (b) Wigner-derived lane signatures with the `(−1)^m` factor symbolic, (c) the H=Z effective-coupling theorem, (d) a concrete leakage value from the projected Maxwell law, (e) the Frechet equivalence of the primitive→bundle map, (f) the previously-confirmed Case B regularization, mouth ℓ-series cross-check, N₄ form equivalence, Y₂₀-via-Gaunt lane signature, and Gaussian wall-integral packet. They should *not* yet take it as evidence for the BdG-derived `B_n, Z_n` coefficients (still hand-typed kernels), the algebraic z₀-drops-out theorem (still parallel-typed cancellation), the IBP boundary-term tracking (still print-and-deleted), or the wall-only no-go theorem (still trivial determinant). Those remain the items for a fourth pass.

The bundle's trajectory is favourable. Round 1 → round 2 was "VACUOUS-dominated → WEAK-dominated." Round 2 → round 3 is "WEAK-dominated → SUBSTANTIVE-LITE-dominated, with the single highest-priority item now done." The remaining gap is methodological (boundary terms, algebraic propagation, mutation testing) rather than physics derivation, with one major exception: the BdG/port computation that should give `B_n, Z_n` substance in step 19.

---

*Audit complete. 6 referee agents, 12 modified step pairs (+ 6 unchanged carried forward), 100% script execution. Net evidentiary delta: ~15 cumulative substantive new physics computations against ~4 newly-introduced anti-pattern assertions in this round. Four of ten round-2 priority items now substantively addressed (#1, #4, #5, #7); one partially (#2); five remain (#3, #6, #8, #9, #10).*
