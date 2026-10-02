**Verdict: CLEAR FOR THIS BOUNDED FULL-WEAK METHOD**

This is a method clearance only. It is not worker clearance, a computed result, or physical acceptance.

## Substantive blockers
None.

## Why the old `local` object is ambiguous and the new partition is not
- `build_case` in `native/c2.py` (lines 643–648) closes the rows by substituting the four slot atoms with the pressure responses. This happens before `extract` and `ReducedActionAssembly` run.
- `ReducedActionAssembly.construct` then collects every term that is linear in a probe jet into `LOCAL`, with no tag for where the term came from.
- The old `LOCAL_MATRICES` is therefore a post-closure collection. It can absorb any reduced pressure piece that is polynomial in the probe jets. Its name cannot prove it is the pre-pressure local part.
- The proposed method partitions the open native rows instead. Slot atoms are still free symbols there, and the partition files show `unknownPressureAtoms: []` and complete name coverage.
- Child counts add up to the full constructors: THETA is 376 local plus 4 pressure, which equals 380 (indices 0–379). U0, U1 and U2 have no pressure children. This matches the saved `native U pressure absence` requirement at `source/composition-worker.py:380-381`.
- Pressure appears only in THETA and E_W, so the method's "four native pressure slots" means four slots per row. The 12 children (4 in THETA, 8 in E_W) are the pressure-bearing terms for those slots. The wording is slightly loose but the structure is right.
- The new method can therefore avoid the ambiguity. The pre-pressure local part is exactly the children without slot atoms, which is the same as the row with its slots set to zero.

## Necessity and sufficiency
- **Necessary:** `pressure/analytic-conclusion.json` states `localSlabPartIncluded: false`, and `pressure/method.md` §1 excludes the local slab part. Without this join there is no complete retained operator.
- **Sufficient for weak continuity:** the five-row operator is the sum of two parts, so continuity follows from these three pieces:
  - the local children, which are finite-order differential operators with bounded smooth coefficients;
  - the inherited pressure certificate;
  - the exactly-once assembly of both.
- **Sufficient only for that:** it does not give coercivity, invertibility, plane-wave extension, constant-height reduction, current, loss or leakage. The method states these exclusions correctly.

## Checks against the actual sources
- **Local symbols:** the census contains no speed or `c_s` names (`speedSymbolNames: []`) and no `rho_m`. The only local constructors are Add, Mul, Pow, Symbol, Integer and Rational.
- **Denominators:** the only local Adds I found are the constant denominators `ωρ_brτ + iρ_br`, which are nonzero at ω=3 and τ=1/10. U rows have no Adds. So the rational-quotient step reduces to constant-denominator checks, and the profile-dependent-denominator branch is a correct refusal case.
- **Bindings:** every local symbol is in `physical-input.json` except `sigma_W`, which is deliberately symbolic. The census's `additionalPhysicalBindingNames` lists exactly the extensions the method describes.
- **Wave jets:** `wave_jet` accepts the `_t`, `_tt` and `_t_dN` suffixes, so the method's `(-3i)^n_t (i/5)^n2 (i/10)^n3 ∂ₓ^n1` rule is consistent with it. The `grad_theta_i → theta_d_i` alias is also in `wave_jet`.
- **Profile rule:** the L-power rule at `c2.py:264-265` matches the saved `profile-scale-join.json`. The `dx` rule at `c2.py:285-296` is consistent with `w1_profile_d1 = L ∂ₓw1`.
- **Grade truncation:** `retained_shape` truncates to a≤1, b≤1 in eta and sigma, which is the same ideal (η², σ²) the method uses.
- **Polynomial and seminorm argument:** the derivative recurrence P_{j+1} = (1−T²)P_j′/10 is correct, so all profile jets are polynomials in T and all their derivatives are bounded. Bound (2) is valid, since ∫(1+|x|)⁻² = 2 and |a| ≤ C_a. The endpoint claim also holds: every jet of order ≥1 carries a factor (1−T²), so it vanishes at T=±1.
- **Local wave derivatives:** local children go up to order 3 (for example `u_1_d1d1d1`), against the pressure sources' maximum of 2. Bound (2) covers any finite n.

## Nonblocking limitations
1. **Pressure scope is inherited.** Continuity of the complete operator holds only within the pressure certificate's scope: cs∈[1,2], real internal momenta, distributional δ→0 limit, ω=3, and the rectangle G. The local part is cs-independent, so nothing about continuity in cs is added here.
2. **The grades are a formal rectangle.** η and σ are independent, and η_bg=1/100 in the physical input is never bound. The native identity σ = W₀η/L (in the saved density context) is not applied. Any later pilot specialization needs its own reviewed grade collapse.
3. **Physical values zero out parts of the rows.**
   - `Lambda_X_0=0` and `Lambda_V_0=0`, so most E_W Λ_X children and the THETA Λ_V child vanish, including four E_W pressure children.
   - Retain their ancestry, as the method requires.
   - Choose the mixed-grade control addresses among surviving nonzero children. This applies to control 1 in particular.
4. **Two density symbols are used.** Local children use `rho_br`=1, while E_W pressure children use `rho_m`=1/10. The method requires agreement only on the intersection, so this inherited difference is not reconciled here.
5. **Time-sign consistency has no dedicated control.**
   - The rule `D_t = −3i` should agree with the explicit `1/(ωτ+i)` memory factors in the local children and with the pressure `1−iΩτ` and Im Ω>0 convention.
   - A wrong sign would not be caught by the three listed controls.
   - I recommend one structural join check in the worker spec.
6. **Worker spec should reuse the existing identity checks.** The pressure children must sum to the saved row-raw pressure part (`composition-worker.py:358-385`). Each slot coefficient must equal the inventory's consumer field. The new worker should add that the local children sum to the row minus the slot terms. This makes "consumed exactly once" an identity rather than an index-disjointness check.
7. **Integrity checks are not physical controls.** The method already says this about the duplicate, missing and wrong-mapping checks. The worker spec should keep it.

## Examined
- **Read in full:** `evidence-guide.md`, `method.md`, `review-prompt.md`, `packet-index.json`, `physical-input.json`, `native/binding-context.json`, `native/profile-scale-join.json`, `legacy/reduced-action-excerpt.json`, `pressure/method.md`, `pressure/analytic-conclusion.json`, `pressure/new-weak-duality-and-order.json`, `inventory/typed-direct-objects.json`.
- **`native/local-source-census.json`:**
  - Read the U0, THETA and E_W symbol tables, constructor-call counts and the tail summary.
  - U1 and U2 counts were seen only through grep, and not their symbol tables.
- **Partitions:**
  - Read the THETA header, one child, the tail index lists and the pressure-index lists.
  - Grepped the E_W pressure children. Looked at only the header of U0.
- **`native/c2.py`:** lines 120–310, 555–700 and 839–880, and the section around `traction_pairing` (940–980).
- **`source/composition-worker.py`:** lines 95–125, 180–260, 322–392 and 960–1010.
- **Only grepped or heads read:** `legacy/numerical-binding.py` (grep), `inventory/address-metadata.json` (first 150 lines), `pressure/address-coverage-metadata.json` (first 120 lines).

## Material coverage gaps
I did not open the following, so I did not check them:
- the other pressure certificate JSONs;
- `pressure/worker.py`;
- the rest of `inventory/`, including `grade-coverage.json`;
- the 4.5 MB address-representatives file;
- the U1 and U2 partition constructors;
- the full 13 MB native constructors.

The census was taken as given and was not independently recomputed. Pressure-side claims rest on the saved certificates and `pressure/method.md`.