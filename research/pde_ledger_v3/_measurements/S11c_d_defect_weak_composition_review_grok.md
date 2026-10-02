I'll read the method packet, evidence guide, and index first, then check the proposed estimates against the saved response formulas and join data.The method is a proposed global weak form. I'll check it against the saved response census, closed density, and the inventory maps rather than the draft alone.The census and closed density match the draft kernels at the symbol level. Next I’ll check the profile tails, height prescription, and whether the old compact-domain constants still cover all real momenta.The saved J density keeps the quadrant ratio, and the direct density keeps a separate reflected root. I’ll now check the Fourier address convention and whether the global height and profile bounds stay valid off the old compact set.The analytic bounds already match the saved formulas. I still need the address records and control applicability before any verdict.CLEAR FOR THIS GLOBAL WEAK-COMPOSITION METHOD

The shifted inverse-square-root estimate supplies a polynomial global bound with no exponential factor in `|k|+|l|`. Local uniform integrability on each fixed compact external set, plus that compact's tail, is the continuity argument the method actually uses for the integrated kernels. The global Schwartz pairing is a separate dominated-convergence step from the polynomial envelopes. The first-height principal-value prescription defines the height action and the normal kernel without treating `q(l)Y(l)` as `C^1`. Absorbing the coefficient fields into `b D_j u` and `c v` matches the saved Fourier addresses and the `2π` duality of the stated convention. The retained grade restriction, both face signs, one whole direct factor, and one native iteration survive. The genuinely new controls have saved addresses they can hit.

No substantive blocker.

## Bounds against the saved formulas

Real frequency is `ω=3`. The historical physical-input value `ω=1` is marked historical. Edges stay `1/5` and `1/10`. `Ω=3+iδ` with `0<δ≤1/10` enters only the response. Address maps record `frequency: 3` and `positiveRegulatorContinuation: false`.

`β=Ω/(10-iΩ)` agrees with both saved writings, `Ω/(1000(-iΩ/1000+1/100))` and `Ω/(100(-iΩ/100+1/10))`. The saved real part is `30/((10+δ)^2+9)`, minimized at `δ=1/10` by `b=3000/11101`. On this rectangle `|μ|≤2/5`, `|a|<1`, `|β|≤3/√109≤2/5`, and `Re β≥b`. First-quadrant roots give `|q+β|≥b`, `|q/(q+β)|≤1`, and `|qi/(qm+qi)|≤1`. The same-quadrant difference bound `|q(p)-q(k)|≤√(|p-k||p+k|)` follows from `|q(p)+q(k)|≥|q(p)-q(k)|`.

`|Ω^2|=9+δ^2≤9.01` and `|Ω^2/cs^2-1/20|≤9.06`, so `|q|≤|p|+√9.06<|p|+4≤4(1+|p|)`. That replaces the compact certificate's `|q|≤5` and endpoint enclosure `[-6,6]`. With `κ_δ=√((9-δ^2)/cs^2-1/20)` and `a_*=√879/20`, `|radicand|≥|κ_δ^2-p^2|` and `max(|p-κ_δ|,|p+κ_δ|)≥κ_δ`, so

`1/|q| ≤ a_*^{-1/2} Σ_{s=±1} |p-s κ_δ|^{-1/2}`

is a valid almost-everywhere envelope. For `δ>0` the radicand has positive imaginary part, so the real comparison points are not complex branch points.

`A(s)=10 s/(4 sinh(5π s))`, `A(0)=1/(2π)`, and `|A'|≤L/4=5/2` hold for every real `s` because `0≤x cosh x-sinh x≤sinh^2 x` for `x=5π s≥0`. The new bound `|A|≤11 e^{-|s|}` splits at `|s|=1` and uses the inherited `A≤10|s|e^{-5π|s|}` only for `|s|≥1`. Then `|A(t)A(Q-t)|≤121 e^{-|t|}` with no factor `e^{|Q|}`. The old `e^{30π}` tail at `|Q|≤6` is not reused.

With `w(t)=(1+|t|)^2 e^{-|t|}`, `sup w=4/e≤2` and `∫w=10`. Splitting at distance `1` gives `∫ w(t)|t-z|^{-1/2} dt≤18` for every `z`, hence `C_q=36/√a_*` for both `q(k+t)` and `q(l-t)`.

The saved iteration density in `saved/reference/right-height-PV-operands.json` is

`qi t Ω^2 (k+t)(2k+t) 10(l-k-t) / [6400 qm (qi+qm) (qi+β)(qm+β)(qo+β) (-iΩ/1000+1/100) sinh(5π t) sinh(5π(l-k-t))]`.

The factor `(-iΩ/1000+1/100)=(10-iΩ)/1000` converts this into the method's `(a μ^2 W L/4) A(t)A(l-k-t)` times `m(m+k) qi/[qm(qi+qm)(qi+β)(qm+β)(qo+β)]`, `m=k+t`. The quadrant ratio stays together. There is no reflected root in `J`. The saved difference identity is `qm-qi=-t(2k+t)/(qi+qm)`. `|a μ^2 W L/4|≤2/5` and `|m(m+k)|≤2 P^2(1+|t|)^2` give `K_J=(4/5)·121·C_q/b^3`.

The saved direct density keeps three roots: `qs=q(l-t)` in `k(2l-t)/(qo+qs)`, and `qh=q(k+t)` in `qi^2/qh` and `k(t+2k) qi/[qh(qh+qi)]`. The prefactor `W L/(4i)` times the two `A` factors reproduces the `2π`-free `5/32` and `sinh` normalization at `W=1`, `L=10`. The `qs` numerator is at most `2P^2+P|t|`. The two `qh` numerators together are at most `18P^2+P|t|`. Each route is therefore bounded by `18 P^2(1+|t|)`, and `(1+|t|)≤(1+|t|)^2`, so `K_D=18·121·(2 C_q)/b^2`. These constants are existence envelopes.

`|i f q(l)|≤4P`, so the normal images of `J` and the direct density grow at most like `P^3`. After `|qo/(qo+β)|≤1`, the iteration coefficient `C` is linear in `|k|` and its normal image is quadratic. Flat pressure and flat normal stay bounded. `inventory/factors/address-full-factor-10-operands.json`, address `1716`, is the minus mixed iteration `C Hwhole(l-k,1,10)+Jwhole_minus` with pressure normal `1`. It does not multiply `Dwhole` by another resolvent or a second middle integral.

`|H(Q)|≤100` for every real `Q` follows from the contact `(W/4)j(Q)`, the bound `|j(Q-s)-j(Q)|≤L/(2π)=5/π`, and `|A(s)/s|≤10 e^{-5π|s|}` outside `|s|≤1`. The saved contact `5(-10k+10l)/(16 sinh(5π(l-k)))` is that contact. The old restriction `|Q|≤6` is not required.

## The five method questions

The shifted estimate does give the polynomial bound without an exponential prefactor in `|k|+|l|`. Continuity in `(k,l,cs,δ)` on each fixed compact uses uniform absolute continuity on bounded `t` intervals, a uniform tail for that compact, and almost-everywhere convergence. A single global pointwise dominant is not assumed. Equations `(5)` and `(6)`, uniform in `δ` and `cs` through `b`, `a_*`, `|μ|`, and `|β|`, then dominate the pairing against Schwartz `X` and `Y`. The resulting object is a continuous bilinear form on `S(R)×S(R)`, with weak continuity in `cs` through both selected matches. It is not an operator-norm bound, a smoothness bound, or a numerical error bar.

The inner height action is `(W/4)f(0)+W/(2i)∫ A(Q)[f(Q)-f(0)χ(Q)]/Q dQ`, with `χ=1` on `|Q|≤1`. Because `A` and `χ` are even, the principal value of `A χ/Q` vanishes, and the subtraction is exactly the saved pairing against `h`. On `|Q|≤1`,

`|B|≤(8/(5b))P_k`, `|B(k+Q,k)-B(k,k)|≤(8√3/(5b^2)) P_k^{3/2} √|Q|`,

using `|qi|≤4P_k` and `|2k+Q|≤3P_k`. The product with `Y` adds `(8/(5b))P_k ‖Y'‖_∞ |Q|`. After division by `Q` the integrand is `O(|Q|^{-1/2})`, including at a real square-root point. For the normal coefficient `Bn=f μ qi qo/(qo+β)`, `|Bn|≤(8/5)P_k` and the variation is `(16√3/(25b^2)) P_k^{3/2} √|Q|`, from `|β|≤2/5` and the quotient `qo/(qo+β)`. The same subtracted integral applies on the lower face. Census normal jets are `+I qo` and `-I qo`. Address `9163` carries normal multiplier `I q(composition_l)` on a slope slot.

Holding the source and consumer fields at real frequency `3` while regularizing only the response is a limiting-absorption definition of that response. It is not analytic continuation of the composed operator.

The forward convention `hat f(p)=(1/(2π))∫ e^{-ipx} f dx` has `hat(fg)(k)=∫ hat f(k-p) hat g(p) dp` and `∫ f g dx=2π ∫ hat f(-l) hat g(l) dl`. Thus `Y(l)=2π hat(c v)(-l)` is the bilinear duality factor, and `X(k)=hat(b D_j u)(k)` is the source product. Address `7997` stores that order explicitly: jet `e_W_d1d1` with spatial orders `[2,0,0]`, nonconstant field `(30-100 I)(I tanh(x/10)+I)/10900`, transfer `l-p` with null constant value, and wave multiplier `-composition_p^2`. Differentiation stays at the original momentum `p`. The transverse and time factors `(i/5)`, `(i/10)`, and `(-i·3)` are mode multipliers. All `34` fields in `inventory/fields.json` are constants or polynomials in `tanh(x/10)` of degree at most `5`, with constant rational coefficients. The recurrence `P_{n+1}=(1-T^2)P_n'(T)/10` is the chain rule, and `|T|≤1` bounds every derivative, so each multiplication is continuous on `S(R)`.

Grades remain `{(0,0),(1,0),(0,1),(1,1)}`, with `(2,0)` and `(0,2)` discarded. `eta` and `sigma` stay independent. Direct multiplicity is `1`, `middleIntegrationOfWholeDirect` is false, and `multiplyWholeTagByResolvents` is false. `Dwhole` and the native iteration enter grade `(1,1)` as separate components. The normal factor is the saved jet `i f q(l)`, not another resolvent. Both faces are present in the census, in the `26` formal controls, and in sampled minus addresses `1716` and `1989`.

The `26` saved controls are algebraic tag movements with nonzero constant rational values: `source-p-to-k`, `omit-source10`, `normal-q-l-to-r`, `omit-consumer10`, `omit-direct`, `double-direct`, and `lower-jet-sign` on the `THETA` and `E_W` faces. `source-p-to-k` moves the wave momentum on an already nonconstant `x`-jet. It is not the Leibniz interchange. `normal-q-l-to-r` already replaces `q(l)` by `q(r)` on a normal slope address whose consumer is nonconstant; address `9163` is that shape, with consumer `(12-40 I)(-I tanh(x/10)-I)/1744`, transfer `r-l`, and null constant value. The method's reuse sentence covers it. `wrong-middle-root` corrupts the iteration middle root `qm` at `(k,l,m)=(3/2,2,30/13)`. It does not touch `qs`. Removing the explicit `H` contact, and corrupting `qs=q(l-t)` inside the saved direct density at nonzero momenta, are new and have saved operands. The direct argument map sends `grazing_qs` to `q(l-t)` and `grazing_qh` to `q(k+t)`.

## Issue classes

Substantive math, domain, or claim issues: none.

Implementation obligations already required by the method, before the conclusion is established: certify all `34` fields and the derivative recurrence; record the global `q`, `β`, profile, and shift constants and the separate pressure and normal degrees; join contact, principal value, and whole tags on both faces with `eta` and `sigma` uncollapsed; run the new addressed controls on address `7997`'s jet class, on address `9163`'s normal consumer, on the saved `H` contact, and on `qs` in the direct density. Restore saved operands only. Label algebraic movements as coefficient movements. Persist operands before guards.

Pure tooling: the later instrument stays under the existing no-deadline pooled guard. This review executed nothing and restored no symbolic payload.

Optional wording: the `q(r)` check can name the existing `normal-q-l-to-r` control so it is reused rather than copied. The constants `100` and `18` are loose existence bounds.

## What a bounded certificate can establish

A successful instrument can establish that every saved coefficient is a constant or a polynomial in `tanh(x/10)`, that the recurrence puts every derivative in a bounded class, and that multiplication preserves `S(R)`. It can establish the envelopes `(2)` through `(6)`, the height and normal Hölder constants, and polynomial growth of degree at most three for the ordinary whole kernels. It can establish the contact and subtracted principal value, the direct and iteration argument maps, both face signs, and one direct factor with one native iteration. Together with the dominated-convergence argument, that gives a finite Schwartz-seminorm bound for a continuous bilinear retained-pressure form on `S(R)×S(R)`, and continuity in `cs∈[1,2]` through both selected matches, including the outgoing limit `δ→0+`.

It does not establish an integral value, a scattering inverse, a loss, a plane-wave solution, a drain, a primitive calibration, an operator norm, or smooth leakage.

## Barriers before any finite near-unity pilot

The pilot still needs a justified trial and test space, an asymptotic excitation, numerical action and endpoint controls, current and finite-inverse sensitivity, a regulator policy, and the applicable reviews. Incoming plane waves are outside Schwartz space. The local slab part and the drain are outside this pressure-response form. Method clearance is not result acceptance, and the proposal does not authorize a finite solve or a defect sweep.

## Packet coverage

The formulas above were read from `method.md`, the response census, `right-height-PV-operands.json`, `left-height-subtracted-PV.json`, `domain-bound-certificate.json`, `whole-tag-definitions.json`, the direct argument map, factor `10`, and `composition-worker.py` through the address constructor and the `26` control definitions. `fields.json` was read in full: `8` constants and `26` tanh polynomials. `responsive-formal-controls.json` was read through `E_W_BALANCE-lower-jet-sign`. `address-representatives.json` is the exact selection of the `608` formal nonzero addresses from `THETA_BALANCE` and `E_W_BALANCE` (`U0`, `U1`, and `U2` contribute no selected indices); source hashes and the index lists are in `packet-index.json` and `address-metadata.json`. Records actually opened include addresses `7956`, `7967`, `7995`, `7997`, `7999`, `8001`, `8020`, `8071`, `8276`, `8736`, `9163`, and `10928`. The `320`-entry grade census and the `13260` original addresses were not read line by line. Their counts, zero statuses (`11232` exact-zero consumer, `1420` exact-zero source jet, `608` formal nonzero), and the constructor that emits them were.