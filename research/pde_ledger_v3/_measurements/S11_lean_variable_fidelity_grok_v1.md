**CLEAR**

Packet identifier (as supplied, not recomputed): **S11 variable-coefficient and flat-interface VC1–VC4 fidelity contract v1**, `aggregate_sha256` `8152df86a31ecaabe89713a3fb178c7c469885f905b9e486207bbc695c02cf45`.

This review is of statement fidelity, hypotheses, index/normalization conventions, and proof logic against the declared local D3 family and D4 odd density. It is not a kernel rerun and not an S11c-wide operator review.

## Assumptions and conventions checked

- `G_ij = ∂_i u_j` with spatial derivative rows. `Point D = Fin (D+1) → ℝ` is time then D spatial coordinates; spatial rows use `i.succ`. Native `coordinate_substitution` uses the same row convention.
- Profiles are prescribed and frozen in the jet: `D3.lagrangian` / `D4.lagrangian` evaluate coefficients at the spacetime point and do not vary them. Hypotheses are `ContDiff ℝ ∞` profiles and `SmoothField u` only. No positivity, transversality, stationarity, or background field equation.
- Time momentum is the `Fin.cases 0` row of the derivative-defined momenta. That row is identically zero even if profiles or fields depend on time. Native checks use spatial coefficient profiles; Lean identities allow spacetime profiles.
- Actions: D3 `L = -[a(div u)^2 + b tr(G^2) + c ∑_{ij} G_ij^2]/2`; D4 `L = -β P/2` with `P = F01 F23 - F02 F13 + F03 F12`, `F = G-G^T`, native `P_D = P`. Local EL is `-div p` from actual jet derivatives. No new integrated first-variation theorem is claimed for variable profiles.
- Interface theorem is a 1D oriented two-interval identity with separate real-line flux representatives, two-sided `HasDerivAt` on each `uIcc`, and a common test trace. Normal is minus to plus. Traction is `n_i ∂L/∂G_ij`; stiffness traction has the opposite sign.

## VC1 — D3 variable-profile EL

The reused momenta in `S11D3Bulk.momentum_eq` are

`p_{r i} = -(a δ_{ri} div u + b ∂_i u_r + c ∂_r u_i)`

for spatial derivative index `r` and field component `i`. That is `p_ij = -(a δ_ij div + b ∂_j u_i + c ∂_i u_j)` in the coverage notation. Time row is zero. `D3.momentum_identity` and `D3.pointwise_variation` freeze the local coefficient value.

`D3.eulerLagrange` is minus the full coordinate divergence of these momenta. `D3.eulerLagrange_eq` plus `correction` is exactly

```
EL_j = (a+b) ∂_j(div u) + c Δ u_j
     + (∂_j a) div u
     + ∑_i (∂_i b) ∂_j u_i
     + ∑_i (∂_i c) ∂_i u_j.
```

The `b`-gradient contracts `∂_i b` with `fieldJet … j.succ i = ∂_j u_i`, which is the transposed index required by `b ∂_j u_i` in the momentum, not `∂_i u_j`. The native wrong-index witness `a=c=0`, `b=x1`, `u=(x2,0,0)` gives actual EL `(0,1,0)` and transposed claim `0`, matching that algebra.

`D3.constant_profile` is definitional recovery of the reviewed constant-coefficient operator. `D3.null_profile_residual` keeps

`(∂_j a) div u − ∑_i (∂_i a) ∂_j u_i`

for `a=-b`, `c=0`. The smooth witness `u=(0,x2,0)`, `a=x1` has `div u = 1` and residual component `EL_0 = 1`; native checks the full vector `(1,0,0)`.

No required correction.

## VC2 — weighted currents and D4 gradient response

`weighted_divergence` is the spatial product rule `div(aJ) = ∇a·J + a div J` (`i.succ` only). `weighted_density` rearranges it without dropping the gradient pairing.

D3 current is the reviewed `J_i = ∑_j[u_i ∂_j u_j − u_j ∂_j u_i]`. `D3.weighted_null_density` is

`L = −(1/2) div(aJ) + (1/2) ∇a·J`.

D4 current is `K_i = (1/2) ∑_j u_j M_ij`. `D4.weighted_odd_density` is the same pattern with `β` and `K`. The `+1/2` gradient contraction is the action factor; constant-coefficient bulk cancellation cannot remove it. Native identities are written at the `Q`/`P` level (`a div J`, `β P`) and are consistent with `L = −Q/2`.

D4 momenta are `−β M/2` on spatial rows, zero time row (`S11D4Odd.momentum_eq`, reused by `D4.momentum_identity`). `D4.eulerLagrange_eq` uses `dualCurl_divergence_zero` (`∑_i ∂_i M_ij = 0`) and leaves

`EL_j = (1/2) ∑_i (∂_i β) M_ij`.

`D4.constant_profile` recovers zero EL. `M` matches `∂P/∂G` via the reviewed dual matrix and `∑_{ij} G_ij M_ij = 2P`. The fully summed Levi-Civita contraction on `G` equaling `2P` is the previously reviewed D4 normalization, restated here, not a new classification of vanishing variable-profile responses.

No required correction.

## VC3 — normal slice and traction

`split_integration_by_parts` instantiates Mathlib `integral_mul_deriv_eq_deriv_mul` on each side. Premises are `HasDerivAt` of `pm,h` on `uIcc a c` and of `pp,h` on `uIcc c b`, plus `IntervalIntegrable` of `dm, dh` on `a..c` and `dp, dh` on `c..b`. The identity is

```
∫_a^c pm h' + ∫_c^b pp h'
  = pp(b)h(b) − pm(a)h(a) + (pm(c)−pp(c)) h(c)
    − ∫_a^c dm h − ∫_c^b dp h.
```

That retains endpoints and the exact interface sign. Oriented integrals impose no `a≤c≤b`; the adjacent-interval reading does. Fluxes are separate `ℝ → ℝ` representatives, so `HasDerivAt` at closed endpoints is a two-sided extension hypothesis, stronger than one-sided weak traces, as documented. Values at `c` need not agree. The same `h` is the common trace.

`compact_endpoint_split` assumes only `h(a)=h(b)=0`, not compact support. The name is slightly broader than the hypothesis; the theorem and fidelity record are accurate. `interfaceWitness_eq` is an actual integral of `p_-=2`, `p_+=5`, `h=1-x^2` on `[-1,0]∪[0,1]`, equal to `-3`, obtained from that IBP theorem.

`jumpPair_zero_iff` quantifies over every finite-dimensional test `h` and is equivalent to flux equality (proved by `Pi.single`). `jumpPair_swap` is the normal reversal. D3/D4 `traction` is `n_i` times spatial momentum. Witnesses: D3 `a=1,b=-1,c=0`, `G_{22}=1`, `n=e_1` gives `-1`; D4 `β=2`, reviewed `witnessJet`, `n=e_1` gives `-1`. These are momentum fluxes; stiffness traction is the opposite sign.

The scalar IBP is not a vector/multidimensional transmission theorem, weak-trace existence result, field-continuity statement, surface-action law, or solution-trace surjectivity theorem. Those limits match the theorems. They are not missing completion hypotheses for this bounded slice identity.

No required correction.

## VC4 — controls, native link, adequacy

The twelve paired mutants in `S11_lean_variable_contract_check.py` are concrete `contract_control` statement mutations, not edits of the general proofs. Author-recorded run2: each rejected mutant has exactly one diagnostic, unsolved `False` in `contract_control`, exit 1. Sixteen positives include the twelve partners plus constant D3 null profile, constant D4 `β=-7`, equal traces, and the quantified trace criterion. Spot-checked records (`d3_omit_gradient_mutant`, `d3_gradient_sign_mutant`, `d4_omit_gradient_mutant`, `interface_omission_mutant`, `interface_sign_mutant`) match that pattern.

Lean witnesses prove the specified nonzero components (`d3_nonzero_response`: `EL_0=1`; `d4_nonzero_response`: `EL_1=1/2`). Native checks report full vectors `(1,0,0)` and `(0,1/2,0,0)`. Derivative-index sensitivity is the tenth native control, not one of the twelve Lean pairs; that split is stated and is adequate for the load-bearing index claim.

The compact native instrument AST-extracts original Q9/coordinate helpers, identifies the D3 trace basis and `P_D=P`, and itself forms `p=∂L/∂G` then `-div p` for prescribed profiles. It does not call `euler_lagrange_from_placeholders` for this extension. Seven identities and ten controls are recorded PASS. Correspondence is symbolic translation checking, not kernel-certified CAS or an S11c-wide link.

The bounded suite is adequate for the stated local EL, weighted-current, slice-IBP, and traction claims, and is accurately limited.

## Required vs optional

No required fidelity corrections.

Optional stronger results, not needed for this contract: a Lean full-vector witness; a Lean statement-mutant of the `b`-index; a variable-profile integrated first variation; one-sided/weak-trace IBP; profile-valued traction as a separate definition; a classification of all profiles with vanishing residual; multidimensional transmission; D5, curved/Sobolev interface, general null Lagrangians, spectrum/scattering, or systematic CAS bridging.

## Verification limits

I did not run Lean, the native instrument, hash recomputation, or a Mathlib rebuild. Compilation, axiom audit (31 standard-axiom declarations), mutation outcomes, reuse of forty recorded objects, package-pin cleanliness, and aggregate hashing remain author evidence. Mathlib is not in the packet; recorded pins and direct import source/object hashes are the stated trust baseline. Unchanged D3/D4 objects may differ from earlier contracts after dependency rebuilds; that does not affect this statement-fidelity review.

Independent check performed here: definitions, quantifiers, index conventions, EL and product-rule algebra, IBP boundary algebra, traction signs/factors, control design, native identification path, and the documented application boundary.
