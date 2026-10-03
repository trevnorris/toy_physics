**Verdict: NEEDS REVISION**

This is a source-only review. I executed nothing, read no peer report or reviewer state, and ran no science. Source clearance here is not result acceptance.

## Blockers

**B1. Route B cannot converge on the inverse-square-root endpoints, so it will refuse deterministically.**
- Every point has `qh=q(k+t)=0` at the cuts `t=±κ−k`. The cuts come from `exact_cuts` (`numerical-library.py:36`) and become panel endpoints in `panels()`.
- J, D_height and D_quadratic carry `1/qh` (`numerical-library.py:46,49,50`), so they are `~c/√(t−t0)` at a panel endpoint.
- For such a panel the G7–K15 difference scales like `c·h^{1/2}`. The tolerance halves on every bisection (`numerical-library.py:166`, `tol/2`), so it scales like `h`. The acceptance test at line 162 therefore fails at every depth.
- The recursion only ends at `require(a<mid<b,'physical adaptive precision stagnation')` (line 151), at roughly 165 levels deep. The first point refuses in its B route, after A24, A48 and A48T124 are done. That consumes the single authorized run, with no 38-point result.
- D_reflected is `c0+c1√x`, so its error scales like `h^{3/2}` and it would converge. The failure affects three of the four primitives.
- `tooling-tests.py` has no quadrature test, so nothing synthetic exposes this.
- **Required correction:** give Route B a predeclared endpoint treatment that terminates on `1/√` endpoints, with its own tolerance allocation. Keep it independent of A's `z²` map.
  - One option is a different regularizing substitution.
  - Another is an explicit physical-variable subtraction. `qh²` is exactly quadratic in `t`, so the local coefficient is available.
  - The sum-of-leaf-error accounting in `adaptive` must still hold.

**B2. The panel-coverage `require` is likely to refuse on an ulp mismatch.**
- `panels()` builds `points=[a+(b-a)*i/n …]` (`numerical-library.py:120`). The last point `a+(b-a)*n/n` is not guaranteed to equal `b`.
- This is most likely on wide intervals such as `[−122, …]`, where `ulp(b−a) ≫ ulp(b)`.
- The next interval starts at the exact `b`, so `panels[i][1]==panels[i+1][0]` (line 122) can fail. That refuses on the first point's A24 route. It is a spurious failure of the "open coverage" check, not a real gap.
- **Required correction:** set `points[0]=a` and `points[-1]=b` explicitly. Check coverage against the exact cut list.

**B3. The slope scale join rests on the author's chosen substitution.**
- `worker.py:149-153` substitutes `w1_profile_d1 → (1−tanh²(x/10))/2`, the ξ-derivative, and checks it against `jp`. The residual therefore proves only algebra, not the scale.
- The saved metadata fact is `nativeRule: L_W**number_of_native_spatial_indices`, with length 10, in `inventory/native-profile-scale-join.json`. The worker copies that file but never reads it.
- My reading is that this rule is consistent with the substituted value, so the risk of a wrong value is low. It is still not a runtime join. A 10× error in j would pass silently, because the saved J/H already assume this j.
- **Required correction:** `require` the saved rule text and length 10. Derive `w1_profile_d1 = L_W·∂_x w(x/L)` from them instead of hand-substituting.

## Optional improvements

- **Numeric constants are not tied to the symbolic ones.**
  - The library's literals (`a=1/(1−3j/10)`, `mu=3/10`, `beta=(30+9j)/109`, `595/100`, the profile formula, `W=1`, `L=10`) are only joined through the symbolic `J.zero`.
  - `runtime-mu` (`worker.py:138`) hardcodes `3` instead of `physical['frequency']`.
  - I checked these by hand and found no mismatch. A numeric evaluation at one point of the saved `oldJ`/D adapters against `middle()` would close the gap at runtime.
- **The H-contact and lower-normal controls are tautological.** Their movement is computed from the same closed form being "mutated" (`numerical-library.py:243-247`). They are not run through the pipeline. Persist the pressure/normal movement instead. Only the reflected-root control actually re-evaluates the kernel.
- **The A48T124 check is insensitive.** The tail bounds are about 1e-29, far below the 1e-9 comparison tolerance, so it can only detect cut-handling bugs.
- **Rule restoration checks are structural only.** They cover length, openness, positivity and tuple identity. Weight-sum and moment checks at load would catch a swapped or corrupted payload. The B-rule centre node is `−3.3e-51`, not zero, which is harmless.
- **Journal volume is unestimated.** Full per-node records likely reach several GB per point. Page-cache growth could repeat the historical cgroup `max` events. The build correctly does not assume they were harmless.

## What I checked and found sound

- **H/J/D joins**
  - The saved J, D and H adapters match `kernel_components` term by term: prefactors, `(qi+β)(qo+β)`, the three separate D terms, and distinct `qh=q(k+t)` versus `qs=q(l−t)`.
  - I confirmed `β=(30+9i)/109` and `a·μ²·WL/4` against the saved coefficients.
- **Physical H product**
  - `hj = j/4 − (10/8)j'` holds.
  - `ĵ=5A`, `ĥ = δ/4 + A/(2it)` and product-to-convolution factor 1 hold under the 1/(2π) convention.
  - The `H(0)=5/(8π)` check passes.
- **Contact and cancellation:** the contact `5A(Q)/4`, the paired density and the even `χ`/`ĵ(Q)` cancellation are correct.
- **Route A mapping**
  - The `z²` map, positive Jacobian `2Lz`, half weights and open nodes are correct.
  - They resolve the `1/√` endpoints, including at coincident roots, where the integrand is analytic in the √ variable.
- **Tails**
  - `|qh+qi| ≥ |qh|, |qi|` and `|q+β| ≥ Re β = 30/109 > b` hold.
  - The polynomial-exponential moments reproduce the J, D_ref, D_height and D_quad constants. The numeric size is about 1e-29.
  - `121e^{-|t|}` follows from the saved `11e^{-|t|}` bound.
- **Collision coverage**
  - Exact `(a,b)` coefficient keys merge exactly coincident labels, and κ is irrational, so no distinct cuts can alias.
  - The 38 points exercise `k+l=0, ±2κ` and `l=k` at the line and at ±1/64, ±1/128 on both sides.
  - Grazing is approached from both sides for `k→±κ` and `l→±κ`. The `(0,0)` point also hits `k+l=0` and `l=k` together.
- **Comparison and evidence**
  - Per-primitive and sum comparisons against A48 use the unchanged tolerance.
  - The journal uses FULL synchronous transactions with failed-panel prefixes preserved and no retry.
  - Nothing claims a packet value, current, loss or global quadrature bound.

## Coverage limitations (not blockers)

- There is no exact external grazing, and the other coordinate is fixed at 1/5.
- Collision offsets shift `l` only, with `k` fixed.
- No corner overlaps are sampled, such as `k=l=±κ` or grazing on `k+l=±2κ`.
- Route B shares the analytic kernel with A, so agreement cannot validate the kernel. This is acknowledged in `build.md`.
- Controls are kernel-level, not addressed packet controls.
- The Fourier, summand-unit, outer-action and outer-resolution work remains open.
- The empirical B error sum is not a quadrature proof.

I could not evaluate whether the symbolic `J.zero` joins succeed in sympy. That is a runtime obligation, and I found no structural reason for them to fail.

**Literal verdict: NEEDS REVISION**