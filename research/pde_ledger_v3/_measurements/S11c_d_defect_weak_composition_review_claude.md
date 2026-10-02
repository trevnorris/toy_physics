CLEAR FOR THIS GLOBAL WEAK-COMPOSITION METHOD

I found no substantive math, domain or claim blocker. I recomputed the constants by hand from the saved formulas. I did not execute anything.

## Packet coverage

- **Read in full:** `method.md`, `evidence-guide.md`, `review-prompt.md`, `saved/direct/closed-density.json`, `saved/reference/retained-response-census.json`, `native-middle-sheet.json`, `first-shape-native-transport.json`, `left-height-subtracted-PV.json`, `continuation/new-bounds.json`, `inventory/native-profile-scale-join.json`, and `inventory/applicability-obligations.json`.
- **Read in part:**
  - The first 651 of 2667 lines of `packet-index.json` (hashes and paths only).
  - The first ~120 lines of `inventory/fields.json`, plus a grep over the whole file. Every non-constant field's text contains only `tanh(composition_x/10)` and rational/complex constants.
  - The first 120 lines of `responsive-formal-controls.json`.
  - A grep of `source/composition-worker.py`. It states the multiplication order as consumer, then response, then source followed by the wave derivative. It uses the `1/(2π)` forward normalization.
  - A grep of `address-representatives.json` and the first 80 lines of `address-metadata.json`. These show the `pressure` slot and the `NATIVE_FLAT` component, with 608 formal-nonzero addresses.
- **Not inspected:** the 320-entry grade census, the per-address `factors/*` files, the typed/whole-tag/argument-map files, and most of the 122 `saved/` files. My checks rest on the formulas actually read. I did not rely on any report's claim.

## Checks against the actual saved formulas

- **Direct density.** The saved `Bc` and `density` match method §4's `Bc` and `GD` term by term, with `qs=q(l-t)` and `qh=q(k+t)` distinct. The prefactor `-25i·t·10(l-t-k)/(16 sinh sinh)` equals `(WL/4i)·A(t)A(l-k-t)` for `L=10`.
- **J.** The saved Jdensity has `qm=q(k+t)`, `(k+t)(2k+t)`, the five resolvent factors and `(qi+qm)`. Its prefactor reduces to `aμ²WL/4·A·A` (`1/640` on both sides).
- **Census coefficients.** The flat, height, slope and mixed coefficients match F0, Fh, Fj and C (mixed iteration `=C·H+∫J`). The direct whole convolution enters once, as `certified_closed_direct_whole_convolution`, and is not multiplied by another resolvent or middle integral. The normal jets are `±i·qo·(...)`, so the face sign is a pure prefactor.
- **Root and profile constants.**
  - `|q_δ|≤|p|+4`: the radicand is at most 9.06, so the bound holds.
  - Same-quadrant `|q1-q2|²≤|q1²-q2²|` and `|q/(q+β)|≤1`: both hold, including on the evanescent imaginary axis at δ=0.
  - `Re β=30/D≥3000/11101` at δ=0.1: this matches `b`.
  - κ_δ lies in `[√879/20, 3]`.
  - Equation (2) is valid even though Ω is complex, because `|q²|≥|Re q²|=|κ_δ²−p²|`.
- **Equations (3) to (6).**
  - Equation (3) holds.
  - In (4), `sup w≤4/e<2`, `∫w=10` and `4·2+10=18` per root give `Cq=36/√a_*`.
  - `KJ=(4/5)·121·Cq/b³` holds, given `0.4·2` from `|a|·|μ|²·WL/4` and `|m(m+k)|`.
  - `KD`: the three numerator terms contribute 2, 2 and 16, so they fit under 18 per route, and `|μ|WL/4≤1`.
  - Both direct routes are covered, and `1/|qs+qo|≤1/|qs|` holds.
- **Height PV (7).** `|B|≤8Pk/(5b)` holds. The difference is `μ·qi·(qi−qo)/((qo+β)(qi+β))` and is bounded by `(8√2/5b²)Pk^{3/2}√|Q|`, so the stated √3 version holds. The normal kernel is `fμ·qi·qo/(qo+β)` with `|qo/(qo+β)−qi/(qi+β)|≤|β||qo−qi|/b²`, which gives the stated `16√3/(25b²)` bound. The proof uses the coefficient's own Hölder-½ bound, so q(l)Y is never treated as C1.

## Answers to the specific questions

1. **Polynomial bound without an exponential prefactor.** Yes. The envelope is a sum of two inverse-square-root terms, and the weight `w` absorbs `(1+|t|)²`. The bound is shift-uniform, so there is no `exp(30π)`-type factor like the one in the old `tailEnvelope`. That old factor came from the `|l−k|≤6` compact tail.
2. **Continuity.** Continuity in `(k,l,cs,δ)` is needed only pointwise, for a.e. `(k,l)`. The global weak continuity then follows from the uniform polynomial bound and dominated convergence against Schwartz `X,Y`. Uniform integrability of `|t−z|^{-1/2}` with moving z (Vitali, p<2) and a.e. convergence supply that pointwise continuity.
3. **Height PV.** Yes. Equation (7) is the inherited `h=(W/2)[δ/2+PV A/(is)]` rewritten with a symmetric subtraction (the `Aχ/Q` integral vanishes by oddness).
4. **Auxiliary regularization.** It is legitimate but not essential. The kernel bounds and quadrant facts hold at δ=0 directly. Real-frequency q is `√` or `i√`, as in the saved Piecewise, and it is the closed-quadrant limit. Source and consumer fields stay real because only the response kernels depend on δ.
5. **Coefficient absorption.** `X=\hat{bD_ju}` and `Y=2π\hat{cv}(−l)=∫e^{ilx}cv` give `∫cv·bD_ju dx` for the delta kernel. There is no stray 2π. The normal jet `i f q(l)` acts on the response `R̂(l)`, and the consumer multiplies afterward, as in the worker's stated order.
6. **tanh recurrence.** Yes. `T'=(1−T²)/10` keeps every `P_n(T)` bounded, so multiplication preserves S.
7. **Grade restrictions and signs.** These are algebraic and prior to the Fourier step, so a linear weak form leaves them untouched.

## Substantive blockers

None.

## Implementation obligations

These are not blockers. A certificate must record each of them with operands.

- Certify all 34 fields as polynomials in `T=tanh(x/10)` with exact constants. Check that no x, unbound symbol or other function remains, and persist the literal residual.
- Re-derive the `H` bound `|H|≤100`, the dropped `|Q|≤6` restriction and the `A≤11e^{−|s|}` step with the actual `L=10` join.
- Assert that the consumer slot set is exactly pressure and normal. A consumer that differentiates the response tangentially would need a different weak treatment, such as an `il` factor or integration by parts onto `cv`.
- Verify the inherited 2π and Fourier normalization of Fh, Fj and H against the constant-coefficient delta limit.
- Verify `J` and `GD` against the saved argument maps at every address, including the minus face. A separate join (not only the hand draft) must check the reflected `qs=q(l−t)`.
- Responsive controls have to be applicable at actual addresses:
  - Swapping `bD_ju` and `D_j(bu)` needs a nonconstant source field with a nonzero jet.
  - The `q(r)` versus `q(l)` control needs a nonconstant consumer on a normal address.
  - The control that removes the H contact must be shown to move the formal coefficient.
  - The reflected-root corruption must be applied at nonzero momenta.
  - Prior controls may be reused only where they are identical.
- Label any formal-coefficient movement as formal, not as a loss signal.

## Optional wording

- "Compact-external tail" is stronger than needed. Pointwise continuity plus the global bound suffices.
- State plainly that the local-slab part, incoming plane waves and the cutoffs 4/6 are outside the claim.

## What a successful bounded certificate could establish

- A bounded, `S×S`-continuous weak form for each retained pressure and normal block. It would carry polynomial weights of order at most 3 and the Hölder-½ height-PV action. It would be continuous in cs and δ at both matches.
- Source-joined envelope constants and a coefficient-field class, plus a verification that the saved addresses are consistent with the weak formulas.
- It would not establish a numerical value, a scattering inverse, loss, coercivity or an operator norm. It would not establish differentiability at matching, and the constants are existence bounds.

## Remaining barriers to any finite near-unity pilot

- A justified trial and test space and asymptotic excitation. Incoming plane waves are not Schwartz.
- Numerical action and endpoint controls (inverse-sqrt endpoints, PV subtraction).
- Finite-inverse and current sensitivity.
- Regulator policy, and the unresolved composed grazing limit and local-slab part.
- Dimensional and unit audit of the source coefficient (inherited, not reconstructed).
- Applicable independent reviews of any finite solver.
- Calibration. cs∈[1,2] is an effective family, not a primitive one.