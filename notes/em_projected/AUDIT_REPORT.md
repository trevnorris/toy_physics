# Critical Referee Audit — Projected Maxwell / Parent-Throat Bundle

**Audited:** all 18 `.md` ↔ `_sympy.py` pairs in `/home/trevnorris/Downloads/em_projected/` (steps 01–19, no step 06).
**Method:** six independent referee agents, each adversarially auditing a sub-bundle for hidden assumptions, tautological tests, and md-vs-py mismatches. All scripts were run; all execute without error.
**Bottom line:** `STATUS: PASS` across the bundle is, in most files, evidentially empty. The code is dominated by tautologies (`x/x = 1`, `solve(eq) → subst → eq = 0`, `simplify(A − A)`), hand-typed answers verified against themselves, and **prints with no assertions at all** in 8 of 18 scripts. The genuine physics — projection of Maxwell from a 4+1 bulk, mouth/throat geometry, Y₂ triple overlaps, wall integrals `M_Σ, K_Σ` — is **not exercised in code**.

---

## 1. Verdict matrix

| Step | File pair | Verdict | One-line reason |
|---|---|---|---|
| 01 | `_readme*` | **WEAK** | Real Gaussian arithmetic for one ratio; Faraday block is `x − x = 0` after substitution; IBP "check" is Leibniz-rule tautology. Misses claims C1.1–C1.4 entirely. |
| 02 | `_covariant*` | **VACUOUS** | **No `assert` anywhere.** Pure pretty-printer. Boundary terms typed but never collapsed. |
| 03 | `_vector*` | **WEAK / VACUOUS** | Faraday/Ampere "checks" are tautologies after substitution; results printed but never asserted. Leak/gauge terms inserted by hand as fresh `Function` symbols. Dead/duplicate code lines 132–138. |
| 04 | `_comparison*` | **SUBSTANTIVE on Gaussians, WEAK on the claim** | Real symbolic integrals; but the bridge `I_WZ ∂f = μ₀ I_WS j` is consumed as input. Case B (sharp slice) regularizes `δ(0)` by fiat. Works in `F^{wν}=0, Γ=0` limit, undermining the "open-system" framing. |
| 05 | `_extension*` | **WEAK with VACUOUS core** | Boxed gauge result `H=Z ⇒ ξ_eff=ξ` "checked" via `simplify(ξ·I_WZ/I_WZ)` — sympy cancelling identical symbols typed by the human. Real content limited to Gaussian numerics and one regulator limit. |
| 07 | `_near_throat*` | **WEAK with VACUOUS portions** | Symmetric-kernel moment expansion is dead code; answer typed in. Mouth-Gaussian ℓ-series `1 − 2ℓ²/λ² + 12ℓ⁴/λ⁴ − 120ℓ⁶/λ⁶` is hand-typed and never compared to the `erfc` closed form sympy did compute. `H=Z` and `S=CZ` "exact" identities reduce to `x/x`. |
| 08 | `_push_bundle*` | **WEAK** | One real check (δP₀ formula) and one grouped-anomaly arithmetic. Several "indep of z₀" assertions check sensitivity w.r.t. symbols absent-by-construction. |
| 09 | `_p2_bridge*` | **VACUOUS** | **Zero assertions.** Boxed claims (z₀ cancellation, b=3a, Ξ₁ formula) derived "for the eye" only. |
| 10 | `_primitive_bridge*` | **VACUOUS** as a check | Zero assertions. Consumes step-09 bridge formulas and primitive formulas as input; outputs primitive linearizations matching the md by inspection only. |
| 11 | `_mouth_taylor_master*` | **WEAK** | Taylor "test" hardcodes `W = e^{−u}` (`μ₁ = 1`); the universal first-moment claim is not tested. Mechanism sieve is a hand-input 2×2 with nonzero determinant. |
| 12 | `_gate_bridge*` | **VACUOUS** | **Zero assertions.** Every md formula typed-in and pretty-printed back. |
| 13 | `_action_master*` | **VACUOUS (with one non-trivial line)** | The `K_η` formula is typed and "verified" against itself. Grouped signature `(1, 1/2, −1)` is the input, not the output. The triple-overlap is never computed. Boundary term silently dropped in IBP. |
| 14 | `_action_candidate*` | **VACUOUS** | **No assertions.** Boxed nonlinear EL, background equation, and `K_η` formula all hand-typed and printed. Sympy never performs the variation that the file claims to certify. |
| 15 | `_weak_axisym*` | **VACUOUS** | **No assertions.** Lane weights `λ = (1, 1/2, −1)` are hardcoded literals; `b = 3a` then follows from grade-school arithmetic. Wall-block obstruction duplicates step 13's `det ≠ 0`. |
| 16 | `_bundle_master*` | **VACUOUS** | All assertions are linear-system round-trips (`solve(eq) → subst → eq = 0`). Adds no original content beyond restating steps 17–18. |
| 17 | `_isotropic_bundle*` | **WEAK / borderline VACUOUS** | **Zero assertions.** Branch sign choice (positive square root for `M_Σ`) never tested. Printed `N₄` differs from md `N₄` by a factor of 9 — not flagged. |
| 18 | `_packet*` | **WEAK** | Linear `solve` round-trip is correct, but the `Y₂₀` angular content and `λ = (1, 1/2, −1)` signature — the only things that distinguish "axisymmetric packet" from "two free symbols" — are absent from the code. No expectation values; the "packet" framing is rhetorical. |
| 19 | `_branch_export*` | **VACUOUS as a branch export** | All assertions are root-substitution tautologies. **No actual branch is exported** — no concrete `R₀(w)`, `β₂(w)`, `μ_η(w)`, `T_w(w)`, `T_Ω(w)`, `K_η(w)` is ever instantiated. |

**Substantive : Weak : Vacuous = 1 : 6 : 11** (step 04 is the one substantive script, in its Gaussian arithmetic).

---

## 2. Anti-patterns observed across the bundle

### 2.1 The `x/x = 1` "structural identity" (steps 05, 07, 17)

The bundle's headline gauge-sector results are:
- step 05: `H = Z ⇒ ξ_eff = ξ` (md §C5.5)
- step 07: `H = Z ⇒ ξ_eff = ξ` and `S = CZ ⇒ μ_eff = Cμ₀` (md §C7.7)

The scripts "verify" these via:

```python
xi_eff_HZ = sp.simplify(xi * I_WZ / I_WZ)            # step 05 line 95
mu_eff_proportional = sp.simplify(mu0 * (C * (z0+ℓz1+ℓ²z2)) / (z0+ℓz1+ℓ²z2))  # step 07 line 166
```

Both are `x / x = 1`. The substitution that replaces `I_WH` (or `⟨H·B⟩`) with `I_WZ` (or `⟨Z·B⟩`) under the hypothesis `H = Z` is performed *by the author at the keyboard* and then handed to sympy as the same symbol on both sides. Sympy is verifying its own ability to simplify a fraction.

A real test would compute both sides as integrals of a generic projection kernel `W` against the *same* function (under the hypothesis `H = Z`) and confirm the symbolic equality of the two integrals.

### 2.2 The `solve(eq) → subst → eq = 0` round-trip (steps 16, 18, 19)

Standard pattern:

```python
sol = sp.solve(eq, var)
assert_zero("...", eq.subs(sol))
```

This asserts only that `sp.solve` did its job; it carries no physics content. Every assertion in step 16 and step 19 has this form, and the two `Check_*` lines in step 18 do as well.

### 2.3 Hand-typed answer verified against itself (steps 01, 03, 13, 14)

```python
# step_13 line 44–46
K_eta = URR0 - d_TwR_R0p + TwRR0*R0p**2/2          # md formula typed in
canonical_L2 = ... K_eta * eta**2 / 2 ...          # uses same K_eta
L2_after_ibp = ... URR0 - d_TwR_R0p + TwRR0*R0p**2/2 ... eta**2 / 2  # also typed in
assert_zero("quadratic wall action after IBP", L2_after_ibp - canonical_L2)
```

The two expressions are written by the same human and are identical by inspection; the IBP that should turn `L₂_raw → L₂_ibp` lives in the human's head. Sympy never performs the IBP.

### 2.4 No assertions at all (steps 02, 04, 09, 10, 12, 14, 15, 17)

**8 of 18 scripts contain no `assert`, no `==` check, no `if … raise`.** They are pretty-printers. "PASS" means "Python finished," not "any claim was checked." A regression that flips a sign or drops a term in any boxed md formula would not be detected by these scripts.

### 2.5 Headline physics objects are pure free symbols

Across all four parent-throat scripts (16, 17, 18, 19) the central new physics objects — the wall integrals
- `M_Σ = ∫ μ_η β₂² dw`
- `K_Σ = ∫ [T_w (β₂′)² + (K_η + 6 T_Ω) β₂²] dw`

— are declared as free `sp.symbols("KSigma MSigma", ...)`. Nowhere in the bundle is `μ_η(w)`, `β₂(w)`, `T_w(w)`, `T_Ω(w)`, `K_η(w)` instantiated as a concrete profile and integrated. The "promotion of K and M from abstract knobs to explicit branch integrals" — advertised as the executive result of the parent-action bundle (step 16 md lines 56, 154–158) — is **a renaming of free symbols**, not a derivation. Step 19's "actual branch export" exports a schema, not data.

### 2.6 Maxwell content is absent from code (steps 01–04, 08–10)

The bundle is titled "projected Maxwell." But:
- No script ever takes a 4+1 bulk `F^{MN}(t,x,y,z,w)`, applies the projection `Proj_W`, and *derives* the projected law symbolically.
- Leakage and gauge-driver terms `Leak_i`, `Gauge_i` are introduced as fresh `sp.Function` symbols (step 03 lines 50–57) and **slotted into the equations by hand**. The structural identity `Leak_i = −Proj_{W'}[Z F^{wi}]` is never instantiated for any concrete `W, Z, F^{wν}`.
- The "homogeneous laws project cleanly" claim is supported by Faraday/Ampere `x − x` checks on a 3+1 antisymmetric tensor that has nothing to do with projection.

### 2.7 Spherical-harmonic / angular content is absent from code (steps 11, 13, 15, 18)

The grouped P₂ lane signature `λ = (1, 1/2, −1)` is the empirical target the parent action is meant to explain. In every script that mentions it, it is **hard-coded as a triple of literal numbers**:
```python
lam20, lam21, lam22 = 1, sp.Rational(1, 2), -1     # step 15 lines 18–20
```
The `Y₂₀ · Y₂ₘ · Y₂ₘ` triple-overlap (which sympy can compute via `sympy.physics.wigner` or `sympy.functions.special.spherical_harmonics`) is never evaluated. The "b = 3a" theorem of the weak-axisymmetric reduction reduces to elementary arithmetic on those three numbers.

### 2.8 Hand-typed series substituted for genuine series expansion (step 07)

```python
# step_07 line 185
IWZ_mouth_gauss_series = sp.simplify(1 + ell**2*(-2/lam**2) + ell**4*(12/lam**4) + ell**6*(-120/lam**6))
```

The line above (line 184) does compute the closed form `(√π λ /(2ℓ)) e^{λ²/(4ℓ²)} erfc(λ/(2ℓ))` symbolically. The script then *hand-types* the asymptotic series and prints them side by side — but **sympy is never asked to verify they match** (e.g., via `sp.series(IWZ_mouth_gauss, ell, 0, 8)` or numeric evaluation at a few `ℓ`). The series happens to be correct (I verified by hand asymptotic of `erfc`), but the script cannot detect a typo.

### 2.9 Sensitivity assertions on absent symbols

```python
# step_08 line 73
assert_zero("Xi_static independent of S2", sp.diff(Xi_static, S2))
```

`S2` was never typed into `Xi_static`, so the derivative is zero by construction. This is asserting "I didn't insert a symbol I didn't insert" rather than "the symbol drops out via cancellation." It is true but vacuous as a referee artifact.

### 2.10 Branch / sign smuggling (step 17)

Step 17 md line 178 picks "the positive square-root branch appropriate to a stable positive-mass wall lane" for `M_Σ`. The script never tests that the negative root is excluded by stability or any other algebraic constraint. The branch choice is editorial.

---

## 3. Md ↔ Py mismatches worth flagging

1. **Step 03 dead code** (lines 132–138): `curlH1_std` is assigned twice with `*0` factors; the first assignment is overwritten and never used. Indicates incomplete editing.

2. **Step 17 N₄ form mismatch** (md line 274 vs. py output): md gives `N_4 = -5(M_Σ+B_2+Z_2)²/(K_Σ-B_0-Z_0)² · N_0`; the script prints `-5 N_0 (B_4+Z_4)²/[9(B_2+M_Σ+Z_2)²]`. After the one-pole substitution `D_0(B_4+Z_4) = 3(M_Σ+B_2+Z_2)²` the two are algebraically equivalent, but the script never executes that substitution and never compares. A casual reader will notice the missing `1/9` and may suspect an error.

3. **Step 02 unfolded boundary term** (printed output): the "exact continuity law" prints `boundary(W·J_w) - Pwprime(J_w) + Pwprime(J_w)` artifacts because the integrals are never collapsed. With no assertion, this misleading display goes unflagged.

4. **Step 12 undocumented Compat construction** (line 59): introduces free symbols `S, T, P_target` not present in the md. No claim is attached.

5. **Steps 02, 03, 12** — the boundary terms in IBP are dropped in the md's reduction step but not tracked in the py. The decay hypothesis is never stated as a sympy assumption (e.g., `Z(±∞) = 0`).

6. **Step 04** uses the closed-system limit (`F^{wν} = 0`, `Γ^ν = 0`) to derive its √2 mismatch; this is the *opposite* regime from the open-system thesis advertised in step 02. The bundle's main takeaway "projection ≠ controlled reduction" is supported by a kernel-overlap demonstration in the closed-system limit, not by a leakage demonstration.

---

## 4. What a substantive bundle would do

The following changes would convert "PASS" from a green-text artifact into evidence:

1. **Convert prints to assertions.** Every `print(sp.Eq(...))` in steps 02, 09, 10, 12, 14, 15, 17 should become `assert_zero(...)` against an independent re-derivation. If both sides are typed by the same human, the assertion adds nothing — see #2.

2. **Derive once, check once.** For every boxed md identity, perform the derivation in sympy from prior inputs (bulk equation, raw Lagrangian, raw moments) and assert the result equals a closed form *independently typed*. The bundle currently types both sides of every asserted equality from the same source.

3. **Compute the angular content.** `sympy.physics.wigner.wigner_3j` or `sympy.functions.special.spherical_harmonics` will produce `λ = (1, 1/2, −1)` from `Y_{20} · Y_{2m} · Y_{2m}`. Hardcoding these as literals begs the central physics question.

4. **Instantiate one branch.** Pick a concrete `R₀(w)`, `β₂(w)`, `μ_η(w)`, `T_w(w)`, `T_Ω(w)`, `K_η(w)` (e.g., Gaussian profiles with named widths). Compute `M_Σ` and `K_Σ` by `sp.integrate`. Verify they satisfy the one-pole and normalization constraints numerically. Step 19 promises this and does not deliver.

5. **Exhibit the leakage.** Pick one bulk `Z(w) F^{wν}` and one `W(w)`. Compute `Leak_ν = −Proj_{W'}[Z F^{wν}]` by integration. Confirm it is nonzero. The "leakage controls the failure of brane theory to close" framing is currently rhetorical because no concrete instance is computed.

6. **Track boundary terms.** Either add an explicit assumption `W → 0, Z → 0 as |w| → ∞` and verify integrands decay, or carry the boundary terms in the assertions. The current treatment drops them silently.

7. **Negative-failure tests.** For each assertion, mutate one coefficient and confirm the assertion would fire. This catches the "x/x = 1" anti-pattern: if changing a sign in the input does not change the assertion's outcome, the assertion isn't testing what its label claims.

8. **Series consistency.** Replace hand-typed asymptotic series (step 07 line 185) with `sp.series(...)` against the closed form sympy already computed.

---

## 5. Per-step priority list for fixes (ordered by severity)

1. **Step 14** — entire boxed nonlinear EL has no computational support. Most central derivation in the parent-action notes; should be derived by `sp.diff` + IBP in sympy, not typed.
2. **Step 19** — promises "actual branch export"; ships only a schema. Add a concrete profile and integrate.
3. **Step 17** — boxed `N_4` printed in inequivalent form; resolve and assert the equivalence. Test branch sign choice.
4. **Step 12** — zero assertions, every md formula consumed as input. Add `diff`-based dependency assertions parallel to step 11's lines 67–72.
5. **Steps 02, 03** — convert prints to asserts; collapse IBP boundary terms; exhibit one concrete leakage term.
6. **Steps 13, 15** — compute `Y_{20} · Y_{2m}²` triple overlap symbolically rather than hardcoding `λ = (1, 1/2, −1)`.
7. **Step 09, 10** — add asserts to the bridge identities that step 08 already gestures at.
8. **Step 07** — replace hand-typed mouth ℓ-series with `sp.series` against the `erfc` closed form.
9. **Step 11** — extend Taylor universality test to a kernel with `μ₁ ≠ 1`.
10. **Step 04** — derive the reduced equation `I_WZ ∂f = μ₀ I_WS j` symbolically; replace the distribution-undefined Case B with a smoothed source.

---

## 6. What does "PASS" mean across this bundle?

To summarize the strongest claim that can be made about each script's `STATUS: PASS`:

- Step 04, 05 §5–§6, 07 §3 and §5 closed forms: **sympy verified the indicated Gaussian / `erfc` integrals.** Real evidence.
- Step 08 (lines 41, 59–62), 11 (lines 67–72, 90–92, 94–107): **sympy verified small linear-algebra corollaries of hand-typed input.** Evidence about the consistency of the input, not about the input itself.
- Step 16, 18, 19: **sympy verified that `solve` solves the equations it was given.** No physics content.
- Steps 02, 09, 10, 12, 14, 15, 17: **Python did not crash.** No verification content whatsoever.

A reader who treats the green-text "PASS" lines as endorsing the boxed md identities is being misled. The boxed identities are typed; sympy is then asked to confirm `x − x = 0`, `x / x = 1`, or `solve(eq).subs(eq) = 0`.

---

*Audit complete. 6 referee agents, 18 step pairs, 100% script execution. Net evidentiary content: roughly two genuinely substantive Gaussian-integration computations and a handful of small linear-algebra corollaries; the remainder is bookkeeping.*
