# Critical Referee Audit — Fourth Round

**Audited:** the 11 `.md` ↔ `_sympy.py` pairs modified after `AUDIT_REPORT_3.md` was delivered (steps 02, 03, 08, 09, 10, 13, 14, 15-md, 17, 18-py, 19); the unchanged pairs (01, 04, 05, 07, 11, 12, 16) carry their round-3 verdicts forward.
**Method:** five independent referee agents, each given (a) the round-3 findings for their bundle and (b) the relevant priority-fix list from round 3, then asked to verify whether each fix is *substantive* or *cosmetic* and whether new anti-patterns were introduced.
**Scripts:** all 18 still execute and print `STATUS: PASS`.
**Bottom line:** another genuine evidentiary gain — from "1:11:6:0" (round 3) to roughly **"6:8:3:0"** (substantive : substantive-lite : weak : vacuous, with step 04 still mixed). **The single largest fix of the round is step 19's BdG/Sturm–Liouville port construction**, which closes the round-3 priority #1 (BdG-derived `B_n, Z_n`) — the central physics gap remaining after round 3. Steps 18 (`det == −1/27` Jacobian replacement of the §2.2 tautology) and 02 (real symbolic boundary tracking with `sp.limit`) and 03 (full end-to-end concrete projection of Faraday/Gauss/Ampere on an explicit four-component bulk potential) are also substantive upgrades. Of the round-3 ten-item next-pass list, **6/10 are now substantively or partially addressed (#1 BdG done, #2 four independent profiles, #6 boundary terms partial, #8 z₀ partial, #9 no-go honesty, #10 negative tests)**, **3/10 remain (#3 step-17 branch-sign physical-root, #4 IBP boundary in steps 13/14, #5 step-17 §2.2 cleanup at lines 73–74)**, and **1/10 was sidestepped via retreat-and-honesty rather than upgrade (#5/step-17 branch-sign — Vieta identities replaced the tautology but no physical-root selection test was added)**.

---

## 1. Updated verdict matrix

| Step | R3 | R4 | Direction | Headline reason |
|---|---|---|---|---|
| 01 | WEAK | WEAK | flat (not modified) | √2 ratio still the only substantive content. |
| 02 | WEAK | **SUBSTANTIVE-LITE** | up (real) | New `boundary_value(expr) = sp.limit(·, w, ±oo)` helper at line 80; `boundary(W·Q)` now carried in §1 line 96; concrete IBP-with-boundary identity verified at line 231; `assert_nonzero` wrong-sign sentinel at line 233; pinned leakage value `leak − 1` at line 278. Round-3 priority #6 substantively addressed. |
| 03 | SUBSTANTIVE-LITE | **SUBSTANTIVE** | up | §5 expanded from single-point leakage to full end-to-end concrete projection on explicit four-component bulk potential `A_M(t,x,y,z,w)`: real `sp.integrate(Wg·F, …)` projection of div B (line 243), three Faraday components (lines 244–255), Gauss + leakage (lines 295–298), three Ampere + leakage components (lines 299–310), plus per-component IBP-with-boundary checks (lines 274–282) and two wrong-sign sentinels (lines 266, 311–314). Pinned `leak1 = √2/4` at line 265. §1–§2 typed-against-self anti-patterns persist but are now redundant decoration. |
| 04 | SUBSTANTIVE / WEAK | SUBSTANTIVE / WEAK | flat (not modified) | Case B Gaussian regulator + `sp.limit` still substantive; §1 factorization asserts still §2.9. |
| 05 | SUBSTANTIVE-LITE | SUBSTANTIVE-LITE | flat (not modified) | Generic-kernel H=Z proof persists. |
| 07 | SUBSTANTIVE | SUBSTANTIVE | flat (not modified) | H=Z proof + mouth ℓ-series cross-check persist. Strongest pair pre-19. |
| 08 | WEAK | **SUBSTANTIVE-LITE** | up (real) | Now at parity with 09/10 on z₀ front: K-elimination identity `(K_norm_p − K_one_pole_p) ≡ compat_direct_p` against z₀-free target (lines 64–93), plus mutation guard. Caveat: `Ptarget` still free symbol — round-3 priority #8 not fully resolved (no `Ptarget = N₀/D₀` propagation). |
| 09 | SUBSTANTIVE-LITE | **SUBSTANTIVE-LITE** | flat-to-up | K-elimination identity added (lines 105–141); round-3 lane-wise `Ξ₁_proj` (lines 173–186) and real linearization checks (lines 86–103) preserved. Same Ptarget caveat. |
| 10 | SUBSTANTIVE-LITE | **SUBSTANTIVE** | up | Round-3 Frechet cross-checks **expanded** from partial to full set (z₀, z₂, z₄, n₀, n₂, n₄ — all six, lines 85–93) via clean `frechet()` helper. K-elimination identity at primitive level (lines 116–125). Strongest of the 8/9/10 trio. |
| 11 | WEAK | WEAK | flat (not modified) | Sibling-gamma test still the lone universality probe. |
| 12 | WEAK | WEAK | flat (not modified) | Round-2 anti-pattern assertions remain. |
| 13 | SUBSTANTIVE-LITE | **SUBSTANTIVE-LITE** | flat-to-up | md (lines 90–98) now **honestly downgrades** the wall-only no-go to "an obstruction in this truncated system rather than a deeper standalone no-go theorem"; mutation `assert_nonzero` added at lines 48–51 (IBP sign flip) and 77–80 (`K_η` sign flip). Script no-go itself still trivial 2×2 with hand-coded `(1/9, 2/3, 1/27)` decoupled from kernels. IBP boundary terms still print-and-deleted. |
| 14 | SUBSTANTIVE-LITE | **SUBSTANTIVE** | up | Round-3 EL fix intact and **strengthened with sign-mutation negative** (line 51). Y₂₀ Laplacian now via `sympy.functions.special.spherical_harmonics.Ynm(2,0)` (lines 133–140) with `−6` eigenvalue check and `+5` mutation negative. Angular sector partially closed: real `Ynm` Laplacian computation now exists, but EL itself still uses flat `(u,v)` chart — the two computations remain stapled together rather than fused. IBP boundary terms still print-and-deleted. |
| 15 | SUBSTANTIVE-LITE | **SUBSTANTIVE-LITE** | up (md honesty) | md (lines 150–154) now **honestly qualifies** wall-only no-go as "linear obstruction inside the reduced wall-only gate system… not being claimed here as a deeper universal no-go theorem." md and py now aligned. Round-3 priority #9 satisfied for step 15. |
| 16 | SUBSTANTIVE-LITE | SUBSTANTIVE-LITE | flat (not modified) | Real Gaussian wall integration persists. |
| 17 | WEAK | **SUBSTANTIVE-LITE** | up (partial) | Round-3 §2.10 lines 69–70 branch-sign tautologies **removed**; replaced with Vieta sum/product identities and quadratic factorisation (lines 78–90). md (lines 182–188) honestly concedes physical-root selection is a convention. New concrete wall-integral example with `M_Σ = √π, K_Σ = 3√π/2` (lines 94–103). **However**: round-3 persistent §2.2 at old lines 51–52, 56–57 has migrated forward verbatim to current lines 73–74 — *not addressed*. Direction is up but Priority #5 is satisfied by retreat-and-honesty rather than physical-root selection upgrade. |
| 18 | SUBSTANTIVE-LITE | **SUBSTANTIVE** | up | Round-3 lines 69–70 §2.2 tautologies **substantively replaced** with: (a) Jacobian determinant pinned to `−1/27` (line 57), and (b) closed-form structural form-checks on `sol[dKΣ]`, `sol[dMΣ]`, compensated D₀₁/D₂₁/D₄₁, and `Ξ₁_comp` (lines 66–73). Exactly what round-3 priority #6 prescribed. Wigner `(−1)^m` symbolic factor preserved (line 27) and reinforced by explicit lane assertions `(1, 1/2, −1)` at lines 49–51. Residual soft §2.3 at lines 71–72 (closed-form gate residuals after closed form already shown equal to `sol`) is the only new blemish. |
| 19 | SUBSTANTIVE-LITE | **SUBSTANTIVE** | up — biggest fix of the round | Round-3 hand-typed integrand kernels `(1/5)β², −(1/3)w²β², (1/6)(1+w²)β², (1/7)R_0β², …` **eliminated**. Replaced with explicit Sturm–Liouville port construction (lines 125–211): `L_ho = −∂_w² + w²`, normalized Hermite eigenfunctions `ψ_n`, eigenvalues `2n+1` (line 166), orthonormality (lines 167–170), resolvent `R(ω) = Σ |n><n|/(eigenvalue+ω²)` (line 153), resolvent equation `(L_ho + ω²) R φ = Π φ` per port (lines 178–180), spectral-sum vs `∫ Π φ · R φ` cross-check (line 158), `B_n, Z_n` extracted via `sp.series(..., omega, 0, 6)` (lines 197–211). Three genuinely independent w-profiles `mu_shape = 1+w²/4, tw_shape = 1+w²/6, to_shape = 1+w²/8` (lines 98–100); `K_η` legitimately fixed by background closure (lines 111–113). Honesty win preserved: `assert_nonzero` on four target residues (lines 235–238) — packet still misses targets, script still doesn't pretend. **Caveat**: spectral basis is the *bare* harmonic oscillator, not a BdG operator linearized around a true parent background — script and md both honest about this. |

**Substantive : Substantive-lite : Weak : Vacuous = 6 : 8 : 3 : 0** (vs. round 3's 1:11:6:0, vs. round 2's 1:4:13:0, vs. round 1's 1:0:6:11). Step 04 remains in mixed (substantive Case B / weak §1) territory and is not counted in either column.

The migration is now striking: VACUOUS-dominated (R1) → WEAK-dominated (R2) → SUBSTANTIVE-LITE-dominated (R3) → **balanced SUBSTANTIVE / SUBSTANTIVE-LITE (R4)** with 6 of 17 unmixed pairs in the SUBSTANTIVE column.

---

## 2. Round-3 priority-list scorecard

The round-3 audit closed with a 10-item next-pass list. Status by item after round 4:

| # | Round-3 priority | Round-4 status | Notes |
|---|---|---|---|
| 1 | **Step 19 — derive `B_n, Z_n` from BdG/port computation** | **SUBSTANTIVELY ADDRESSED** | Sturm–Liouville eigenproblem on `L_ho = −∂_w² + w²` with normalized Hermite basis. Resolvent equation verified per port. `B_n, Z_n` are now extracted from `sp.series` of spectral matrix elements rather than typed. **Caveat:** bare HO, not parent-linearized BdG. Notes acknowledge this. The decorative-kernel anti-pattern is eliminated. |
| 2 | **Step 19 — four independent w-profile wall coefficients** | **SUBSTANTIVELY ADDRESSED** | `mu_shape = 1 + w²/4`, `tw_shape = 1 + w²/6`, `to_shape = 1 + w²/8` are three genuinely distinct polynomial weights. `K_η` is closure-derived from background equation — legitimate dependency, not free fourth knob. Round-3 critique that three of four coefficients shared the `1+R₀` skeleton is resolved. |
| 3 | **z₀ algebraic cancellation via `Ptarget = N₀/D₀`** | **PARTIALLY ADDRESSED** | All three (steps 8, 9, 10) now exhibit z₀ cancellation as an *elimination identity* `(K_norm_p − K_one_pole_p) ≡ compat_direct_p` against a z₀-free target, plus mutation guards. This is genuinely better than round-3 parallel typing — a mistyped sign on B₀ inside one of the `solve(...)` calls breaks the assertion. **However**, `Ptarget` is still a free symbol in all three; the deeper claim that z₀ propagates through `Ptarget = N₀/D₀ = N₀/(K − B₀ − Z₀)` is not exercised. |
| 4 | **IBP boundary terms tracked symbolically** | **PARTIALLY ADDRESSED** | Step 02 now does this properly: `boundary_value(expr) = sp.limit(·, w, ±oo)` helper (line 80), `boundary(W·Q)` carried in §1 (line 96), concrete IBP-with-boundary identity at line 231, sign-mutation sentinel at line 233. Step 03 §5 follows the same pattern. **Steps 13 and 14 still print-and-delete** — the boundary density appears as text in printed output but no symbolic carrier propagates through downstream assertions. |
| 5 | **Step 17 branch-sign physical-root selection** | **NOT ADDRESSED (sidestepped)** | Round-3 §2.10 lines 69–70 tautologies are gone, but they were not replaced with a second-derivative test, cost-functional minimization, or asymptotic-BC discriminant. Instead they were replaced by Vieta sum/product identities (lines 86, 89) and a quadratic factorisation check (line 84) — which are mathematically correct algebraic identities but do not select a physical branch. md (lines 182–188) and script (line 128) now honestly *concede* that branch selection is not fixed by the one-pole algebra alone. The script no longer claims to do branch selection; the round-3 critique was that it shouldn't, but the upgraded test isn't there either. |
| 6 | **Step 18 lines 69–70 §2.2 cleanup** | **SUBSTANTIVELY ADDRESSED** | Exactly per round-3 prescription: `det == −1/27` (line 57) plus closed-form structural form-checks on `sol[dKΣ], sol[dMΣ]`, the compensated D₀₁/D₂₁/D₄₁, and `Ξ₁_comp` (lines 66–73). The blunt §2.2 is gone. Residual soft §2.3 at lines 71–72 (closed-form gate residuals after the closed form is already shown equal to `sol`) is a minor blemish, not the same anti-pattern. |
| 7 | **Step 03 §1–§2 Faraday/Ampere via real projection** | **ADDRESSED (in §5, not §1–§2)** | The author did not rewrite §1–§2; instead §5 was massively expanded to include full end-to-end concrete projection of div B, three Faraday, Gauss + leakage, three Ampere on an explicit four-component bulk potential `A_M`. The §1–§2 typed-against-self anti-patterns persist at lines 98–111 and 142–163 but are now **redundant decoration**, since §5 verifies the same physics by genuine projection that would catch real coefficient mutations. Net: priority satisfied by addition rather than rewrite. |
| 8 | **Wall-only no-go honesty (steps 13, 15)** | **SUBSTANTIVELY ADDRESSED** | Step 13 md (lines 90–98): "an obstruction in this truncated system rather than a deeper standalone no-go theorem." Step 15 md (lines 150–154): "a linear obstruction inside the reduced wall-only gate system… not being claimed here as a deeper universal no-go theorem." Both prose clarifications are exactly what was asked. Scripts themselves remain trivial 2×2 systems, but the prose no longer over-claims. |
| 9 | **Negative-failure tests bundle-wide** | **SUBSTANTIVELY ADDRESSED** | The `assert_nonzero` mutation pattern from round-3 step 19 has propagated. Round-4 mutation tests now appear in: step 02 line 233 (wrong-sign IBP), step 03 lines 266, 311–314 (wrong-sign IBP, mutated Faraday), step 13 lines 48–51, 77–80 (IBP sign flip, K_η sign flip), step 14 lines 50–51 and 140 (mutated `U_Σ` derivative, mutated Y₂₀ eigenvalue), steps 08/09/10 (z₄-sign flip in K-elimination compat). Pattern is now standard rather than exceptional. |
| 10 | **Step 14 angular sector — full S² covariant Laplacian** | **PARTIALLY ADDRESSED** | The Y₂₀ Laplacian is now computed via `sympy.functions.special.spherical_harmonics.Ynm(2, 0, θ, φ)` with the genuine `(θ,φ)` Laplacian operator and `−6` eigenvalue check at lines 133–140, plus `+5` mutation negative. **However**, the EL derivation itself (lines 22–49) still uses a flat local-orthonormal `(u,v)` chart for the angular kinetic term — the `Ynm` check is a separate sanity verification, not yet fused into a single EL derivation. The two computations are now stapled together (both real, both substantive) rather than unified. |

**Score:** 6/10 substantively addressed (#1, #2, #6, #8, #9, plus #4 partial / #10 partial); 1/10 partial without strong upgrade (#3); 1/10 sidestepped (#5); 1/10 satisfied by addition rather than rewrite (#7); 0/10 entirely untouched. Compare round 3: 4/10 substantive, 1/10 partial, 5/10 skipped.

The four newly-addressed items (#1, #2, #6, #8) plus the partially-addressed #4 and #10 cover the highest-physics-content remaining gaps. The two persisting weaknesses are: the deeper z₀-via-Ptarget propagation (#3) and the step-17 physical-root selection (#5).

---

## 3. The single most important fix: step 19 BdG/Sturm–Liouville port construction

The round-3 audit's #1 priority — "derive the `B_n, Z_n` from a BdG/port computation independent of targets" — was the central remaining physics gap in the parent-action chain. Round 4 closes it. The relevant code at `step_19_parent_throat_action_actual_branch_export_sympy.py:125-211`:

```python
# Lines 133-136: Normalized Hermite eigenfunctions
psi_n = lambda n: sp.hermite(n, w) * sp.exp(-w**2/2) / (sp.pi**Rational(1,4) * sp.sqrt(2**n * factorial(n)))

# Line 142: Self-adjoint Sturm-Liouville operator
def ho_apply(expr): return -sp.diff(expr, w, 2) + w**2 * expr

# Line 139: Eigenvalues 2n+1
# Line 166: Verify L_ho psi_n = (2n+1) psi_n  -- substantive
assert_zero('Hermite eigenequation', ho_apply(psi_n(n)) - (2*n+1)*psi_n(n))

# Lines 167-170: Verify orthonormality
# Line 153: Resolvent
def resolvent(phi, omega):
    return sum(psi_n(n) * sp.integrate(psi_n(n)*phi, (w, -oo, oo)) / ((2*n+1) + omega**2)
               for n in modes)

# Lines 178-180: Verify (L_ho + omega^2) R phi = Pi phi per port
assert_zero('resolvent equation', ho_apply(R_phi) + omega**2 * R_phi - Pi_phi)

# Line 158: Spectral sum vs explicit matrix element cross-check
assert_zero('spectral/matrix-element port agreement',
            spectral_sum - sp.integrate(Pi_phi * R_phi, (w, -oo, oo)))

# Lines 197-211: Read off B_n, Z_n from low-frequency series of B(omega), Z(omega), N(omega)
B_series = sp.series(B_omega, omega, 0, 6).removeO()
B_0, B_2, B_4 = B_series.coeff(omega, 0), B_series.coeff(omega, 2), B_series.coeff(omega, 4)
# ... same for Z, N
```

This is the strongest single fix of the round and arguably of the entire audit chain. The hand-typed prefactors `1/5, 1/3, 1/6, 1/7, 1/8, 1/9` from round 3 are gone, replaced by natural denominators `(ω²+1)(ω²+5)(ω²+9)` arising from the `2n+1` eigenvalues for `n=0,2,4`. The `B_n, Z_n` are now objects whose values follow from real `sp.integrate` matrix elements against an orthonormal basis, with the resolvent equation explicitly verified per port. A regression that mistypes one of the polynomial profile weights `mu_shape, tw_shape, to_shape` would alter the `B_n, Z_n` series and fire the downstream wall-coefficient assertions.

The honest qualifier: the spectral basis is the *bare* harmonic oscillator `−∂_w² + w²`, not a BdG operator linearized around the actual parent background. The notes do not claim otherwise — they describe this as a "concrete operator-resolvent export" rather than the BdG of *this* parent action. The script delivers a working port construction with the right algebraic structure; it does not yet deliver "the BdG of *this* potential gives these coefficients." That is the natural round-5 ask.

---

## 4. Other clean wins of round 4

1. **Step 02 — symbolic boundary tracking with `sp.limit`.** Lines 80–82 add a `boundary_value(expr) = sp.limit(·, w, ±oo)` helper. Line 96 carries `boundary(W·Q)` as a first-class term in the IBP identity rather than dropping it. Lines 229–231 verify `boundary_Q = boundary_value(Wg·Q_concrete) == 0` for the Gaussian profile. Line 231 verifies the IBP-with-boundary identity. Line 233 is a sign-mutation sentinel that fires under wrong-sign IBP. Round-3 priority #6 done at the script level — the boundary terms are no longer print-and-deleted in step 02.

2. **Step 03 §5 end-to-end concrete projection.** Lines 229–242 instantiate an explicit four-component bulk potential `A_M(t,x,y,z,w)`; lines 243–255 verify projected div B and three Faraday components by real `sp.integrate(Wg·F, (w, -oo, oo))`; lines 257–266 pin the leakage value at `√2/4`; lines 268–310 verify projected Gauss + three Ampere components with leakage on the same potential; lines 311–314 are a wrong-sign mutation sentinel for Faraday. The §1–§2 typed-against-self anti-patterns persist but are now redundant decoration. Round-3 priority #7 satisfied by addition.

3. **Step 10 — full Frechet derivative cross-check sweep.** Lines 85–93 introduce a clean `frechet(expr) = Σ ∂f/∂var · slip` helper and assert `dlin(F) == frechet(F)` for **all six** induced bundle quantities `z₀, z₂, z₄, n₀, n₂, n₄`. A mistyped formula in any primitive `Z_n` or `N_n` definition would fire the assertion. Round-3 had partial coverage; round-4 has the complete set with a uniform helper.

4. **Step 14 — Y₂₀ Laplacian via `sympy.functions.special.spherical_harmonics`.** Lines 133–140 compute `Ynm(2, 0, θ, φ)` and apply the genuine spherical Laplacian `(1/sinθ) ∂_θ(sinθ ∂_θ Y) + (1/sin²θ) ∂_φ² Y`, then assert the eigenvalue is `−6` and the `+5` mutation negative fails. This upgrades the round-3 separately-verified `+6` eigenvalue into a real `Ynm`-based check.

5. **Step 18 — `det == −1/27` Jacobian replacement of §2.2 tautology.** Line 57 pins the coefficient-matrix determinant to a specific nonzero rational. Lines 66–67 are closed-form structural checks on `sol[dKΣ], sol[dMΣ]` that would fire on any mistyped coefficient in the gate definitions. Lines 68–70 add compensated denominator-moment closed forms. Line 73 adds the `Ξ₁_comp` closed-form check. This is the textbook right-way replacement of the round-3 §2.2 anti-pattern.

6. **Steps 13, 15, 17 — wall-only no-go honesty.** All three md files now qualify the no-go as a narrow obstruction in a reduced wall-only system rather than a universal theorem. The prose no longer over-claims.

7. **Step 17 — concrete wall-integral example.** Lines 94–103 add a closed-form `M_Σ = √π, K_Σ = 3√π/2` evaluation against an explicit Gaussian background, cross-checked symbolically — non-trivial new substantive content.

8. **Bundle-wide propagation of `assert_nonzero` mutation tests.** Round 3 had this pattern only in step 19. Round 4 has it in steps 02, 03, 08, 09, 10, 13, 14, 17, and 19 — at least one mutation guard per script in nearly every load-bearing assertion family. This is the round-3 priority #10 satisfied.

---

## 5. New issues introduced or persisting in round 4

### Persisting from round 3:

1. **Step 17 lines 73–74** — the round-3 persistent `solve→subst→0` at old lines 51–52, 56–57 (`P2.subs(N2, N2_const_closed) == 0` and `P4.subs({N2:..., N4:...}) == 0`) **migrated forward verbatim**. This is the textbook §2.2 anti-pattern, untouched across two audit rounds despite being flagged. The §2.10 branch-sign tautology was removed; this older one was not.

2. **Step 13 lines 82–89 wall-only no-go** — the script still solves a 2×2 linear system with hand-coded `(1/9, 2/3, 1/27)` constants and asserts the trivial-solution literal. The constants are decoupled from any actual kernel structure. md is honest now; script is unchanged. Round-4 mutation guards at lines 48–51 and 77–80 cover the IBP and `K_η` checks but not the wall-only no-go itself.

3. **Steps 13 and 14 IBP boundary terms** — still print-and-deleted. Round-4 step 02 demonstrated the `sp.limit` + `boundary(W·Q)` carrier pattern works; it has not been propagated to steps 13 and 14.

4. **`Ptarget` as free symbol** in steps 08, 09, 10 — the K-elimination identity is genuine, but the deeper round-3 priority of propagating z₀ through `Ptarget = N₀/D₀` remains unsatisfied.

5. **Step 14 angular-sector unification** — the EL derivation uses a flat `(u,v)` chart; the `Ynm(2,0)` eigenvalue check is a separate computation. They are not yet fused into a single covariant derivation.

### New in round 4 (small):

1. **Step 18 lines 71–72** — soft §2.3 redundancy. `K1.subs({dKΣ: closed_dK, dMΣ: closed_dM}) == 0` and the `H_even` analogue, where `closed_dK, closed_dM` are by lines 66–67 already shown equal to `sol[dKΣ], sol[dMΣ]`. The substitution back into the gates is therefore redundant with the form-check at lines 66–67. Not the same blunt §2.2 as round 3 — it is a §2.3 typed-against-equivalent rather than a `solve→subst→0` — but it adds nothing the form-check did not already give. The round-3 anti-pattern was honestly upgraded; this is a minor leftover.

2. **Step 13 wall-only no-go script unchanged** — flagged above. Strictly speaking not "new" since round 3, but worth noting that the md honesty edit was not accompanied by a script edit linking the constants to actual kernels. A genuine no-go would have those `(1/9, 2/3, 1/27)` arise from the kernel-derived β-overlap structure rather than be hand-typed.

These are far smaller than the round-3 introductions (step 17 §2.10 lines 69–70 and step 18 §2.2 lines 69–70). Round 4 introduced essentially no fresh anti-patterns of comparable severity.

---

## 6. Net evidentiary delta vs. round 3

The bundle gained roughly **15–20 new genuinely substantive assertions** in round 4:

- Step 19 BdG/Sturm–Liouville block: ~10 new substantive assertions (Hermite eigenequation, orthonormality 6-pack, resolvent equation per port, spectral-vs-matrix-element cross-check, Parseval norm split, mutation test for dropped `n=4` mode, closure derivation of `U_scale`, `M_Σ, K_Σ` against three independent profiles).
- Step 03 §5 end-to-end: ~12 new substantive assertions (concrete projected div B, three Faraday, Gauss+leakage, three Ampere with leakage components, per-component IBP-with-boundary, two wrong-sign sentinels).
- Step 02 boundary tracking: 4–5 substantive assertions (boundary discharge, IBP-with-boundary, pinned leakage, sign-mutation sentinel, continuity-IBP-with-boundary).
- Step 10 Frechet sweep: 6 substantive assertions (z₀, z₂, z₄, n₀, n₂, n₄ all via `frechet()` helper).
- Step 14 Y₂₀ Laplacian: 2 substantive assertions (`Ynm(2,0)` eigenvalue, `+5` mutation).
- Step 18 Jacobian + form-check: 5 substantive assertions (`det == −1/27`, four closed-form structural checks).
- Steps 8/9/10 K-elimination: 3 substantive assertions and 3 mutation guards.
- Bundle-wide `assert_nonzero` mutation guards: ~6–8 new ones across steps 02, 03, 13, 14, 17.

Against ~2 newly-introduced soft §2.3 assertions (step 18 lines 71–72) and one persistence (step 17 lines 73–74, carried forward unchanged).

**Evidence-to-anti-pattern ratio is dramatically inverted again** — round 4 added roughly 40+ substantive assertions for ~3 new/persistent anti-patterns. This is the strongest round in the chain on evidentiary density.

---

## 7. What still needs to happen for a fifth pass

In rough order of importance:

1. **Step 17 lines 73–74 §2.2 cleanup** — the persistent `P2.subs(N2, N2_const_closed) == 0` and `P4.subs({N2:..., N4:...}) == 0` should be replaced with the same pattern step 18 used in round 4: assert the linear-system determinant or matrix non-degeneracy, then form-check `N2_const_closed` and `N4_const_closed` against an independently derived structural target.
2. **Step 17 branch-sign physical-root selection** — the Vieta identities are correct but do not select a branch. A second-derivative test on the action `S[M]`, an asymptotic-BC discriminant, or a stability-functional minimization should be added. Alternative: explicitly remove the "branch sign" framing from the script and md, accepting that the script computes the algebra and the convention is exogenous.
3. **Steps 13, 14 IBP boundary terms via `sp.limit`** — propagate the step-02 pattern. A `boundary_value` helper plus a carried `boundary(quad_density)` symbol would lift these from print-and-delete to the same level as step 02.
4. **Step 14 EL ∘ Ynm fusion** — define `R_proj(t, w, θ, φ) = a(t, w) Y₂₀(θ, φ)` and substitute into the angular kinetic term `T_Ω · ∇_S² R · ∇_S² R` directly, deriving the EL with the `+6` factor in a single computation. This closes the round-3 priority #10 fully.
5. **z₀ propagation through `Ptarget = N₀/D₀`** — define `Ptarget = (N₀ + ε n₀)/(K − B₀ − Z₀ − ε z₀)` and recompute the K-surfaces with this substitution. Confirm z₀ cancellation algebraically through this stronger structure. Steps 08, 09, 10 all share the same one-line modification.
6. **Step 19 BdG-around-parent-background** — replace the bare HO `−∂_w² + w²` with the linearization of the actual parent action around `R_0(w) = e^{−w²/2}`, i.e. `L_BdG = −(d/dw)(T_w d/dw) + (V'' evaluated on background)`. The current script's structure should accommodate this swap with minimal rewriting; the spectral cross-check pattern would carry over directly.
7. **Step 13 wall-only no-go via real kernel structures** — replace the hand-coded `(1/9, 2/3, 1/27)` with overlap integrals of actual β-profiles. Alternatively, delete the script block and let the md statement stand on its own.
8. **Step 03 §1–§2 rewrite (cosmetic)** — the typed-against-self lines 98–111 and 142–163 are now redundant given §5; they should either be deleted or rewritten to use the projection pattern from §5 for consistency.
9. **Negative-failure pattern propagation to steps 04, 05, 11, 12, 15, 16** — the `assert_nonzero` mutation guard pattern reached most of the modified scripts in round 4 but did not back-propagate to the scripts that were not modified this round. At minimum, each substantive `assert_zero` in those scripts should gain a paired mutation guard.
10. **Step 04 §1 §2.9 cleanup** (carried over from round 2) — still unaddressed. Sibling-gamma test in step 11, anti-pattern bundle in step 12, also still unaddressed.

The first three items are the highest priority. After fifth-pass items 1–4, the scoreboard would likely move to roughly 9:7:1:0.

---

## 8. Headline restatement

**The round-4 revision delivers six of the round-3 priority list's ten items, including the central physics gap of the bundle (step-19 `B_n, Z_n` derivation via Sturm–Liouville port construction).** The hand-typed integrand kernels are gone, replaced by a real spectral construction with verified eigenequation, orthonormality, resolvent equation, and matrix-element cross-check. The four wall coefficients are now derived from three genuinely independent polynomial profiles. The honesty win — `assert_nonzero` on target residues — is preserved.

Step 18's `det == −1/27` Jacobian replacement of the round-3 §2.2 tautology is the cleanest example of the right-way correction in the chain so far: the audit's specific prescription was implemented exactly. Step 02's `sp.limit`-based boundary tracking does the same for round-3 priority #6, at least within step 02 itself.

The persisting weaknesses are narrower than after any prior round: step 17's persistent §2.2 at lines 73–74 (carried forward from round 3 unchanged), step 17's branch-sign physical-root selection (sidestepped by retreat), `Ptarget = N₀/D₀` deep propagation in steps 8/9/10, and IBP boundary tracking in steps 13/14. None of these are central physics gaps; all are methodological tightenings.

A reader who sees `STATUS: PASS` across 18 scripts in this revision should treat the pass as evidence for: (a) the boxed nonlinear EL of the parent action, with sign-mutation guard; (b) Wigner-derived lane signatures with the `(−1)^m` factor symbolic; (c) the H=Z effective-coupling theorem; (d) end-to-end concrete projection of all four projected Maxwell laws (div B, Faraday, Gauss + leakage, Ampere + leakage) on an explicit four-component bulk potential, with a pinned leakage value; (e) the Frechet equivalence of the primitive→bundle map across all six induced quantities; (f) z₀ cancellation as an elimination identity against a z₀-free target, with sign-mutation guard; (g) **a concrete Sturm–Liouville port construction yielding the `B_n, Z_n` low-frequency coefficients from a verified spectral resolvent** (with the caveat that the basis is the bare harmonic oscillator, not a parent-linearized BdG); (h) symbolic IBP-with-boundary in step 02; (i) the `det == −1/27` even-gate non-degeneracy with form-checks on the closed-form solution. They should *not* yet take it as evidence for: branch-selection physics in step 17 (still convention-dependent), the `Ptarget = N₀/D₀` z₀-propagation theorem (still parallel-typed at the K-level), IBP boundary tracking in steps 13/14 (still print-and-deleted), or a full BdG-around-parent-background spectral derivation in step 19 (still bare HO).

The bundle's trajectory remains favourable. R1 → R2 → R3 → R4 is "VACUOUS-dominated → WEAK-dominated → SUBSTANTIVE-LITE-dominated → balanced SUBSTANTIVE/SUBSTANTIVE-LITE." Round 4 closed the central physics gap (BdG-style port for `B_n, Z_n`) and propagated the negative-failure-test pattern to most of the bundle. The remaining ask is a handful of methodological tightenings plus the parent-background BdG swap.

---

*Audit complete. 5 referee agents, 11 modified step pairs (+ 7 unchanged carried forward), 100% script execution. Net evidentiary delta: ~40+ new substantive assertions against ~3 new or persistent anti-pattern instances. Six of ten round-3 priority items now substantively addressed (#1, #2, #6, #8, #9, plus partial #4, #10); one partial without strong upgrade (#3); one sidestepped (#5); one satisfied by addition (#7); zero entirely untouched.*
