I'll review the bounded full-weak method from the packet itself: first the evidence guide and method, then the source and operand files they name.The packet index names the method and the source operands. I'll read those joins next: partition, wave and profile rules, bindings, and the inherited pressure scope.**CLEAR FOR THIS BOUNDED FULL-WEAK METHOD**

The missing piece of the retained weak operator is the native pre-pressure local slab. Joining those children, once, to the completed pressure form is necessary for that operator and sufficient as a method for it. The argument stays a method: a future instrument still has to emit the polynomial and partition certificates, and nothing here accepts scattering, an inverse, leakage, or a near-unity defect run.

Substantive blockers: none.

## What the join is

The target in `method.md` is the complex-bilinear form on \(S(\mathbb{R})^5\times S(\mathbb{R})^5\),

\[
B_g(v,u)=\sum_{r,A,n}\int v_r(x)\,a_{r,A,n,g}(x)\,\partial_x^n u_A(x)\,dx+B_{\mathrm{pressure},g}(v,u),
\]

with real \(\omega=3\), edge momenta \(1/5\) and \(1/10\), strict rest bulk, \(W=1\), \(L=10\), and retained grades \(G=\{00,10,01,11\}\). Eta and sigma stay independent. The integrals are the definition of the form; this method does not evaluate them.

That local sum is necessary. `pressure/analytic-conclusion.json` sets `localSlabPartIncluded` to false, and `pressure/checks.json` limits the existing claim to a Schwartz pressure contribution, with `scienceAcceptance`, `planeWaveScattering`, `finiteInverse`, and `loss` all false. The native rows still contain the local children counted below.

It is sufficient for this operator target because the local children are the complement of the four native slots, the pressure form already contains native iteration once and the closed direct convolution once, and the method adds those local children under the same jet, profile, grade, and epsilon conventions. Continuity of a finite sum of continuous bilinear forms is then the method argument. It becomes a runtime certificate only after the polynomial reconstructions, the constant nonzero grade-zero denominators, and the once-each address join exist. It does not become physical acceptance.

## Old `local` versus the fresh partition

The old object can contain reduced pressure effects. `legacy/reduced-action-excerpt.json` builds `LOCAL_MATRICES` by splitting probe-column expressions into probe jets and integrals. `legacy/numerical-binding.py` then stores `local` from addresses `('local', n, i, j)` on an already indexed finite row. Either stage sits after reduction, so a coefficient can be a reduced pressure remnant and still be called local.

The fresh reading uses a different object. `source/composition-worker.py` (`broad_row_census`) splits the literal `slab_operator` summands for `LAB_HELD` / `RHO4_CONSTANT`. A child is pressure when a symbol name contains `delta_p` or `d_w_`; every other child is local. The saved flags are `localPartUnchanged: true` and `localPartRestoredAsScience: false`.

Checked against the partition files:

| Row | Local indices | Pressure indices |
|---|---|---|
| U0, U1, U2 | 517 each, indices \(0..516\) | empty |
| THETA_BALANCE | 376 | \(51,52,101,102\) |
| E_W_BALANCE | 981 | \(67,68,117,118,395,396,399,400\) |

Totals are 2,908 local and 12 pressure, matching `native/local-source-census.json` (`totalLocalChildren`, `totalPressureChildren`). In THETA the local index list jumps \(50,53\) and \(100,103\), so the two sets are complements. Pressure child 52 is affine in `d_w_delta_p_plus` and uses `rho_m`. The next child, 53, has empty hits and is a local term in `e_W`, `eta_bg`, `sigma_W`, `w1_profile`, and `w1_profile_d1d1`, with `rho_br` in the denominator. E_W child 68 is the same kind of slot term; child 69 is local. Selecting the hit-free children before the four slots are replaced keeps the pressure slots out of \(a_{r,A,n,g}\).

## Conventions that the formulas actually join

The position-space jet in the method matches `wave()` in `source/composition-worker.py`:

```python
out = (-sp.I*3)**spec['timeOrder']
for n, mom in zip(spec['spatialOrders'], (momentum, Rational(1,5), Rational(1,10))):
    out *= (sp.I*mom)**n
```

Under the inherited transform \(\widehat{\partial_x f}(p)=(ip)\hat f(p)\), the direction-1 factor \((ip)^{n_1}\) is \(\partial_x^{n_1}\) in the \(x\)-integral. Time and the two edge directions remain the constant factors \((-i\cdot 3)^{n_t}\), \((i/5)^{n_2}\), and \((i/10)^{n_3}\). `native/c2.py` `wave_jet` and `jet_spec` both rewrite `grad_theta_k` to `theta_dk` before counting orders. Those names are present in the local census (`grad_theta_1`, `grad_theta_2`, `grad_theta_3`). The unrestricted sector is the right one: the transverse and longitudinal branches in `wave_jet` are a different decomposition, and the method does not use them.

The one-dimensional profile map matches `profile()` in the same composition worker, with \(L=10\):

- \(w(x)=(1+\tanh(x/10))/2\), \(m(x)=(1-\tanh(x/10)^2)/3\)
- a direction-1 jet of order \(j\) becomes \(L^j\partial_x^j\) of that profile
- any index 2 or 3, including mixed jets, becomes 0

`inventory/fourier-and-unit-provenance.json` records `profileDependsOnlyOnCoordinate: 1`. Transverse profile names are really in the local source (`w1_profile_d2` and `m1_profile_d3d3` occur, including U0 child 329), so the zero map has to be applied and the zero cells kept. Eta and sigma stay symbols through that substitution. `Inputs` skips `sigma` when it builds `profiles`, and `bind()` applies `densityMap` and `profileEqualities` only. The stored density tuple does contain the equality \(\sigma_W=W_0\eta_{\mathrm{bg}}/L_W\), and that equality is not one of the applied maps.

`Inputs.at_source` is only the factor \(L^{\mathrm{len(indices)}}\) on a three-argument coefficient function. It does not zero transverse derivatives. The written one-dimensional map is the rule to execute.

## Bindings, quotient, and the seminorm

`speedSymbolNames` is empty, and the global symbol union has no `c_s0`. Effective \(c_s\in[1,2]\) stays inside the inherited pressure response. The physical-input value \(c_s0=10\) and historical \(\omega=1\) stay out; the binding-context numeric value \(\omega=3\) is the override. Every local child census count for `epsilon_shape` equals the local child count (2,908 altogether). Degree one is still a runtime check, because a count is not a power.

Local-only names that are absent from the pressure numeric subset and present in `physical-input.json` are `G_W_u`, `k_W`, `kappa_W`, `mu_S`, `mu_W`, and the gamma coefficients listed in `additionalPhysicalBindingNames`, including both `*_15` families. Their values belong to that same file, with equality on the intersection. `eta_bg` is on that census list and remains a grade symbol; the physical value \(1/100\) is not a substitution. Shared names already in the binding context, including `Lambda_V_0=0`, `tau_A`, `tau_X`, `rho_br`, and `Lambda_A_0`, bind exactly. A term killed by the exact zero `Lambda_V_0` stays in the ancestry as an exact zero.

The quotient is the existing grade recurrence modulo \((\eta^2,\sigma^2)\), retaining only \(G\), with the excluded remainder stored. A grade-zero denominator that still depends on a profile is unresolved. The samples that do bind cleanly, such as THETA child 7, have denominators in \(\omega,\rho_{\mathrm{br}},\tau_A\) and \(I\) only, hence constant complex numbers after the numeric map. Bare `I` appears in those constructor strings. The call census itself is only `Add`, `Mul`, `Pow`, `Symbol`, `Integer`, and `Rational` (727 `Add`, 9 `Rational`).

Once a coefficient is an exact polynomial \(P(T)\), \(T=\tanh(x/10)\), with constant coefficients and a finite nonzero constant denominator, \(P_{j+1}=(1-T^2)P_j'(T)/10\) stays in that polynomial ring and is bounded for \(|T|\le 1\). With \(C_a\) the sum of the absolute coefficient values,

\[
\Bigl|\int v\,a\,\partial^n u\Bigr|\le 2\,C_a\sup_x(1+|x|)^2|v(x)|\,\sup_x|\partial^n u(x)|,
\]

because \(\int_{\mathbb{R}}(1+|x|)^{-2}\,dx=2\). The monochromatic and edge factors are constants inside \(a\). The finite sum over cells is continuous on \(S(\mathbb{R})^5\times S(\mathbb{R})^5\). Adding the inherited continuous \(B_{\mathrm{pressure}}\), under the once-each join, preserves continuity. The positive imaginary shift remains response-only. No coercivity, self-adjointness, invertibility, or scattering uniqueness follows.

Endpoint values at \(T=\pm 1\) are local bookkeeping. Every positive-order profile jet vanishes there by the same recurrence. That limit is not a constant-height reduction of the nonlocal pressure response and does not extend the form from Schwartz tests to incoming plane waves. `source/composition-worker.py` already records `nonzeroConstantHeightResponseReduction: NOT_COMPUTED`.

## Assembly and controls

The combined record is the 2,908 hit-free children once, the 12 pressure children once through the inherited form, and `certified_closed_direct_whole_convolution` once. `pressure/whole-definition-certificate.json` has `wholeDirectOnce` and `nativeIterationOnce`. `inventory/typed-direct-objects.json` keeps `Rprod`, the factored denominator, and `Dwhole` distinct, with `multiplyWholeTagByResolvents` false. U-row pressure coverage is `EXACT_ZERO_CONSUMER`, so those rows do not already carry a pressure block. Reading a status flag is not the join; row, field, grade, epsilon, and the saved law have to match.

The three controls are aimed at rows that exist:

1. Mixed grade is present before binding. THETA child 53 contains both `eta_bg` and `sigma_W`. The witness has to be a cell whose movement is exactly nonzero after binding.
2. Direction-1 profile jets are present (`w1_profile_d1d1` on local children 53 and E_W 69). Dropping that \(L^j\) is a real map change. Explicit powers such as `L_W**(-1)` already sitting in a product are part of the coefficient, separate from the jet map.
3. Positive \(x\)-order wave jets occur together with nonconstant profile factors (`u_1_d1`, `e_W_d1`, and profile symbols in the same rows). The Leibniz difference \(a\partial^n u\) against \(\partial^n(au)\) tests derivative placement. A transverse jet that has already become a constant factor is not that witness.

If binding cancels the selected movement, the instrument reports that the witness is unavailable. It does not invent one, and it does not reuse the pressure control that multiplies by a flat-response sample.

## Nonblocking limitations

- Execute the written one-dimensional profile map. `at_source` supplies the \(L\) power and still differentiates a function of all three coordinates.
- `additionalPhysicalBindingNames` includes `eta_bg`. Eta and sigma remain independent grade symbols. Intersection checks keep \(\omega=3\) and keep \(c_s0=10\) out of both the local rows and the pressure response.
- Global continuity inherits the pressure certificate’s stated conditions: Schwartz tests, response-only regularization, and no plane-wave or constant-height extension. The local polynomial identities and the disjoint once-each join are still future certificates.
- Controls 1 and 2 are fail-closed on an exactly nonzero post-binding movement. Control 3 already says to report an absent witness. A cancelled candidate is reported the same way.
- Endpoint matrices are local limits only.

## Files examined, and what was not

Examined: `method.md`, `evidence-guide.md`, `review-prompt.md`, `packet-index.json`, `physical-input.json`, `native/binding-context.json`, `native/profile-scale-join.json`, `native/local-source-census.json` (row headers, call census, symbol union, totals, and sampled children), partition records for U0, THETA_BALANCE, and E_W_BALANCE (index lists and local/pressure constructor samples), `native/c2.py` (`WAVE_NAMES`, `wave_jet`, `Inputs.at_source`, `Inputs.dx`, `infer_dimensions`), `source/composition-worker.py` (`broad_row_census`, `bind`, `profile`, `wave`, grade quotient, direct-whole assembly), `legacy/reduced-action-excerpt.json`, `legacy/numerical-binding.py` through the `local` matrix assembly, and the pressure scope records `analytic-conclusion.json`, `checks.json`, `global-parameter-domain.json`, `whole-definition-certificate.json`, `grade-route-certificate.json`, `new-weak-duality-and-order.json`, `method.md`, and the relevant `worker.py` sections, plus `inventory/typed-direct-objects.json`, `fourier-and-unit-provenance.json`, and the opening of `grade-coverage.json`.

Not read expression-by-expression: all 2,908 local constructor strings, the 4.5 MB `inventory/address-representatives.json` payload, and every pressure bound file beyond the scope and join records above. No coefficient was restored or recomputed. Those omissions do not leave the partition, jet, profile, binding, or continuity argument unspecified.