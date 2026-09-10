# Grok (grok-4.6 high) review — S11c-d SHARED PHYSICS spec v2 (round 2, DOCUMENT branch)

VERDICT: NOT-SOUND (7 must-fix). Verbatim (a 1-line reasoning-trace preamble precedes the report).

---

I'll form an independent view from the governing sources first, then read the S11c-d spec and report anything that would change what an engine computes or what the spec may claim.Sources of truth next: c2 house format, N6 disposition, and the profile-class consult — still without opening the spec.I’ll keep building the independent view from S11c-a/c2 house format, N6, and the power-counting/confinement pieces.I have the independent picture. Next I’ll read the consults for context, then the spec itself, and run the order/bridge checks.I’ll verify the real c2 export keys and the S11c-b mixing channels, then run the order and bridge checks.# Independent physics review — S11c-d SHARED PHYSICS spec (v2)

**Artifact:** `research/pde_ledger_v3/directives/S11c_d_SHARED_PHYSICS.md`  
**Role:** orchestrator-written physics spec (review-until-clear).  
**Verdict: NOT-SOUND.**

Independent view was formed from `S11c_decisions.md` (N5, N6, N7, N10–N15, S11c-d row), the c2 step record + N6 disposition, `S11c_c2_SHARED_PHYSICS.md`, `S11c_a_SHARED_PHYSICS.md` §2a/§5a, and the real `scripts/S11c_c2_exports.py` keys, then the spec was read against that view. Profile-class consults were used as context only.

---

## Must-fix

### F1. Consume-set still names a nonexistent c2 key, and lists the increment as a consumed operand

**Spec §1a (ll. 92–100):**

> the **field carriers** `s11cc2FieldTheta`, `s11cc2FieldeW`, `s11cc2Fieldu{1,2,3}` …  
> ⛔ **There are NO c2 term-origin, parity, self-energy-increment, or six-§3d-re-adjudication export rows**

**Spec §1b (ll. 108–110):**

> the emitted **`s11cc2ClosedSlabOperator` + `s11cc2ClosedCouplingKernel` + increment VALUES**. These are the operands S11c-d consumes.

**Why it is wrong.** Mechanical lookup on the real delta (`scripts/S11c_c2_exports.py`, 70 keys):

| named in spec | in export |
|---|---|
| `s11cc2ClosedSlabOperator` | PRESENT |
| `s11cc2ClosedCouplingKernel` | PRESENT |
| `s11cc2FieldTheta` | **ABSENT** |
| `s11cc2Fieldtheta` | PRESENT |
| `s11cc2SelfEnergyIncrement` | **ABSENT** (emit-only; dropped from the delta) |

This is the same defect class v1 was folded for (“supplied c2 export rows that do not exist”). `s11cc2FieldTheta` is not a camelCase variant an engine can guess through F9: the live atom is `s11cc2Fieldtheta`. And §1b’s “these are the operands S11c-d consumes” puts the **increment** in the consume-set after §1a just forbade binding it. A builder that binds the spec’s names fails; a builder that binds the increment binds a row that is not there.

**Minimal fix.** Write the live keys (`s11cc2Fieldtheta`, …). Strike “increment VALUES” from the consumed-operand sentence; c2 SOUND status of the increment is a parent fact, not an S11c-d import.

---

### F2. The homotopy parameter `λ` is not a defined object, so the weak coefficients are not unique

**Spec §3c (ll. 290–294):**

> On a **one-parameter physical contrast homotopy** `(η,σ_W) → λ(η,σ_W)` at **fixed** `L_W` (fixed shapes `w₁`, `m₁`), the amplitude is `O(ελ)` … `C = J_conv/J_in = O(λ²)`

**Spec §3d (ll. 307–308):**

> the leading **amplitude** coefficient `∂_λ(A_H/ε)|₀` … **and** the leading **fraction** coefficient `lim_{λ→0} C/λ² = ½C''(0)`

**Why it is wrong.** S11c-a §2a supplies the identity `σ_W ≡ η W̄₀/L_W`. At **fixed** `L_W` there is only one remaining bookkeeper: `η`, with `σ_W` tracking. A map `λ(η,σ_W)` at fixed `L_W` is overconstrained (two arguments, one degree of freedom) and never names `λ`.

That is not cosmetic. Along the physical ray `σ_W = κ η` with `κ = W̄₀/L_W`, a first-shape amplitude `A/ε = a_η η + a_σ σ_W` gives

```text
∂_η(A/ε)|₀ = a_η + κ a_σ
∂_{σ_W}(A/ε)|₀ = a_η/κ + a_σ
```

The two “leading amplitude coefficients” differ by `κ`. The same factor squared hits `C/λ²`. Two engines that pick `λ=η` vs `λ=σ_W` emit different matching coefficients for the same physics. N12 asked for an `(ε,η)` order-count, not an unnamed reparametrization.

**Minimal fix.** Name the path: `λ ≡ η` at fixed `L_W` and fixed shapes (hence `σ_W = η W̄₀/L_W` tracks; that is a path, not a freeze of the formal `(ε,η,σ_W)` multigrade). Keep the independent-grade amplitude as a separate emitted object.

---

### F3. Asserting the Born coefficient is “nonzero” leaks a target and is false on the named class

**Spec §3d (ll. 306–308):**

> Emit **separately** … the leading **amplitude** coefficient `∂_λ(A_H/ε)|₀` **(nonzero, the Born vertex content)**

**Spec §3 header (ll. 234–236):** “states no component value, sign, order, parity, or grade in this document.”

**Why it is wrong.** M2 withholds expected values. “Nonzero” is a builder-visible acceptance criterion (N7’s exact failure mode).

It is also not a theorem of the class. The class discriminant is the **zero-transfer** moment `m̂₁′(0) = m₁(+∞)−m₁(−∞) ≠ 0` (§1c, correct). The on-shell vertex is the form factor at the physical momentum transfer. For a genuine interface (jump `= 1`) whose gradient is a compact ramp,

```text
FT(m₁′)(q) = sinc(q) ,   FT(0) = 1 ,   FT(nπ) = 0  (n = ±1, ±2, …).
```

So `∂_λ(A_H/ε)|₀` can vanish at finite `kL_W` with `∫ m₁′ ≠ 0`. An engine that treats a computed zero as a defect will reject a legal kinematic node, or iterate until the coefficient is nonzero.

**Minimal fix.** Emit the coefficient as a computed object with no value/sign claim. The class condition stays `m̂₁′(0) ≠ 0` on the §5c discriminant, which is not the on-shell Born coefficient.

---

### F4. The N6 control has no legal mutation site, and `Φ` is applied to objects it is not defined on

**Spec §5a (ll. 366–378):**

> **Both routes construct the S11c-d MIXING RESPONSE from the imported CLOSED operator/kernel** (⛔ NOT a fresh face level-set/graph linearization …)  
> route 1 … imported Eulerian closed operator/kernel  
> route 2 … **native material-coordinate closed construction**  
> `S11CD_REP_COVARIANCE_RESIDUAL = route1 − Φ(route2)` after transforming the SOURCE, the incoming/outgoing **modes**, and the **FLUX pairing** under `Φ`

Then the tilt probe (ll. 390–392): reverse one face’s slope “through that route’s own **carrier factory**, ⛔ **not by editing the imported operator**.”

**Why it is wrong.**

1. **No site for the N3/N4 probes.** The consume-set is one Eulerian closed operator/kernel per `(α,ρ)` (F1). S11c-d has no face-normal carrier factory; that factory is c2/S11c-a. Mutating a face slope without editing the imported operator and without a fresh level-set rebuild has nowhere to act. The vacuous alternative is `A−A` on a kernel that was never rebuilt — exactly the control N6 forbids (`decisions:98–104`; c2 §5c: “⛔ Do not emit an `A−A` control where … the mutated source is structurally absent”).

2. **Route 2 is not in the imports.** c2 exported one closed operator in the Eulerian face basis. The material closed operand lived in the N6 companions, not in `S11c_c2_exports.py`. “Native material-coordinate closed construction” is an unnamed c2-scale rebuild.

3. **`Φ(modes, flux)` is the c2 type error, relocated.** c2 §5c / step record: `I_{M→E}` is already in the common Eulerian face basis; **no separate `T` on the increment**; `R_cov` is source-naturality of `μ`, `ms − source_terms(μ_E.subs(Φ), V_E)`. `Φ` is the constitutive field map. It is not defined on scattering modes or on the flux bilinear form. Asking for `route1 − Φ(route2)` on the mixing response lets each engine invent a `Φ`-action (or apply a θ-shift that annihilates, the measured c2 failure).

N6 as required by the decision list — two independent constructions of **this** step’s object, then one-sided corruption of tilt vs advection, `∇W_bg→0` rejected, `RHO4_CONSTANT` as computed absence — is the right *control*. This write-up does not specify a constructible pair of routes for the mixing.

**Minimal fix.** Fix the object: two constructions of the **mixing insertion** (distorted waves + flux pairing) on the *imported* closed kernel, in Eulerian vs material coordinates at fixed anchoring, differenced in the common Eulerian basis with **no** `Φ` on the amplitude. Define the covariance residual, if any, as naturality of the **named maps that actually act** (incoming data / measure), not `Φ(A_H)`. Put the N3/N4 one-sided probes on a source that exists at this step (e.g. the imported kernel’s explicit `∂w₁` / `u·∇ρ₄` factors, rebuilt one-sided from those atoms), or say the kernel-N6 was c2’s and S11c-d only tests insertion independence. Keep `RHO4_CONSTANT` as computed absence, not `A−A`. Do not corrupt an anchoring.

---

### F5. `J_conv` dumps bound capture into a flux, contradicting N13 and §3b

**Spec §3c (ll. 285–288):**

> `J_conv` = that channel’s **lost** flux = continuum (thickness/bulk) **plus** bound capture.

**Spec §3b (ll. 252–267) / N13 (`decisions:130–133`):**

> the two are **distinct emitted objects** … a true bound state carries **no** asymptotic outgoing flux, so a **capture rate** requires a named protocol … ⛔ a pole residue + coupling alone is **not** a conversion probability — either define the capture protocol or emit only the spectral overlap.

**Why it is wrong.** A stationary scattering flux is a continuum object. A true bound pole has no outgoing flux. Adding “bound capture” into `J_conv` either (i) puts a non-flux overlap into a flux (dimensionally illegal under the overlap-only option §3b explicitly allows), or (ii) invents a capture rate the spec refused to define. Either way the conversion fraction `C = J_conv/J_in` is no longer the N7 flux-normalized FORM, and the two N13 channels are no longer distinct emitted objects. “Energy stays in the slab” is exactly the reading N13 forbids; folding capture into the same `J` as radiation is that reading with a plus sign.

**Minimal fix.** `J_conv` = continuum / bulk outgoing flux only. Bound capture stays `S11CD_BOUND_MODE_SPECTRAL_TEST` (+ protocol or overlap). If a total photon-loss fraction is wanted, name it separately as a sum of two **same-dimension** probabilities after the protocol exists — not as a term in `J`.

---

### F6. “Supported only where `∇μ_R ≠ 0`” drops channels the same spec requires in the vertex

**Spec §1a (ll. 89–91):**

> the **full** mixing vertex (tilt `∇w₁`, modulus-gradient `∇m₁`, and N4 advection channels; ⛔ not `∇w₁` alone)

**Spec §2 (ll. 202–203) and §5d (ll. 431–432):**

> supported only where `∇μ_R,bg ≠ 0`  
> FORM … `∝ k·a`, supported where `∇μ_R,bg ≠ 0`

**Why it is wrong.** N3: an in-plane displacement across `∇W_bg` is normal face motion at the coupling’s own gradient order and changes every interface law. S11c-b extracts the **whole** off-diagonal block, not a `μ_R`-only projection. V3_STEP_PLAN:1179 says `∇μ_R ≠ 0` *mixes* the families; it does not say tilt/advection are absent. The word **only** lets an engine project `K` onto `∇m₁` and drop the tilt and N4 pieces that §1a just required. That changes the mixing and the FORM.

**Minimal fix.** Keep “constitutive mixing is driven by `∇μ_R,bg`” as the N7/V3 support statement for that channel. State that the emitted vertex/FORM is the **full** imported off-diagonal (tilt + modulus-gradient + N4), and that a `μ_R`-only projection is an ablation, not the object.

---

### F7. `m₁(±∞)` and `kL_W` presuppose an unnamed interface axis

**Spec §1c (ll. 159–164):**

> asymptotically constant, `∫∂_ξ · dξ ≠ 0` … `m₁(−∞) ≠ m₁(+∞)`  
> `m̂₁′(0) = m₁(+∞) − m₁(−∞)`

**Source (S11c-a §2a):** `ξ ≡ y/L_W` with **three** in-plane coordinates `y`. Limits `±∞` and a scalar `kL_W` are not defined on a 3-vector argument.

**Why it is wrong.** A localized interface is a 1D transition along an edge-normal, uniform along the edge. Without naming that axis, `m₁(+∞)`, the zero-transfer moment, and the form factor `m̂₁′(kL_W)` are different objects for two engines that Fourier along different in-plane directions. That is a missing premise, not “wrong on a different input.”

**Minimal fix.** State the class as: profiles that depend on a **named** in-plane normal `n̂` (e.g. `y¹`), with `ξ = n̂·y/L_W`, limits `ξ→±∞`, and `kL_W` built from the momentum transfer along `n̂`. Transverse-to-the-edge momenta stay live parameters.

---

## Physics that does hold (not findings)

**Profile class / grades.** A localized *interface* (`∫ m₁′ ≠ 0`) rather than a bump (`= 0`) is the right non-global object for a slit edge (N5). `w₁`/`m₁` independent, interface condition on the mixer `m₁`, no invented constitutive tie — correct (N14/N15). `η`, `σ_W`, `kL_W` are distinct (`σ_W = η W̄₀/L_W` independent of `η` via `L_W`; `kL_W` via `k`). An extra `σ_W→0` / `L_W→∞` expansion *does* change the class (WKB / Riemann-Lebesgue on `m̂₁′(kL_W)`). Setting `η→O(1)` inside first-shape-order operators *is* an incomplete operator where omitted `O(η²)` pieces compete. Bound pole as **conditional**, 1D-well theorem refused, empty set permitted: correct for an interface (different asymptotes, nonlocal/non-Hermitian thickness operator). Confinement as a computed object, not a claim: correct (N10/N13).

**Order bookkeeping (item 2).** Linear mixing:

```text
A = ε (a_η η + a_σ σ_W) + O(ε · 2nd shape)
J_in = O(ε²) ,   J_conv = |A|² = O(ε² · (η,σ_W)²) ,   C = J_conv/J_in = O((η,σ_W)²)
```

`ε²` cancels; `C` is independent of `ε`. Emitting **both** absolute `O(ε²η²)` (or `O(ε²λ²)` once `λ` is named) and fractional `O(η²)` is the right N12/N10 split: an `O(ε²)` *vertex* would be the excluded nonlinear program; an `O(ε²η²)` *observable of a linear amplitude* is not. Do not collapse to a single `O(εη)` / `O(ε²η²)` label — the spec is right to refuse that, once F2 is fixed.

**Strong-edge bridge (item 3).** Verified:

```text
A_H = −i ε sin(η G) ,   C = sin²(η G) = η² G² + O(η⁴) ,   C(ηG = nπ) = 0
```

Born slope `∂_η(A_H/ε)|₀ = −i G ≠ 0`, yet `C=0` at finite coupling. A nonzero weak coefficient is not a lower bound on `C_strong`. `C'(0)=0` when `C=O(λ²)` along a smooth path through 0, so an amplitude slope is not a fraction slope. Order-unity edge as a **new** construction (piecewise-uniform / matched interior-exterior), not a reduction of the first-jet kernel, and `O(1)` reductio withheld: honest (N7). Do not call the lab bound “`F(1)` of an amplitude.”

**c2 honesty (item 4).** Per-engine SOUND vs N6-thread-only cross-engine; matched covariance zeros as `(0)−(0)` not operand AGREE; operand DEBT 40/76/18 UNADJUDICATED and material to this consumer; leftover SHAPE not “just thickness”; F/G WITHDRAWN; both `R_N6 = 18/288` and `R_cov` no-nonzero preserved; uniform limit computed, not the withdrawn F. That block is the right supplied premise (modulo F1’s increment slip).

**N11, chain, T7, blind WL, F9, no denylist, N15-at-the-right-layer, FORM-at-d as an N2 boundary refinement** — present and not weakened.

---

## Nits (do not by themselves hold SOUND)

1. **§2 “uniform Born is not equivalent to two-asymptote DWBA.”** At retained first-shape order, away from thresholds, they coincide: `G = G₀ + O(η)`, `K = O(η)`, `ψ = ψ₀ + O(η)` ⇒ mixing `= η G₀ K₁ ψ₀ + O(η²)`. After the spec’s own re-expansion, the continuum object *is* the uniform-mode matrix element of `K`. Keep two-asymptote as threshold/regular-domain organization and for `L±` kinematics; do not claim a retained-order difference. Thresholds remain the named exception; bound-pole resummation remains exempt.

2. **§5a `∇W₀→0`.** `W_0` is the reserved constant (`N14`); the live field is `W_bg`. Write `∇W_bg→0` / `η→0`.

3. **§3b “Not a Bloch band holds iff a pole is found.”** Bloch is out of class whether or not a pole exists. Say: a found pole is a bound/resonance object, not a Bloch band; an empty set is not a band either.

4. **§5c “require the mixing to move.”** House form-ablation language; N6 probes correctly refuse a nonzero target. Same discipline here: emit the residual, do not assert motion.

---

## Verdict

**NOT-SOUND.** Must-fix F1–F7 all change what an engine binds or computes, or what the spec may claim:

| ID | What changes |
|---|---|
| F1 | bind failure / wrong consume-set (v1 class) |
| F2 | non-unique weak matching coefficients |
| F3 | leaked nonzero target; false at form-factor nodes |
| F4 | N6 not constructible; `Φ` type error on modes/flux |
| F5 | conversion fraction mixes two N13 observables |
| F6 | tilt/N4 channels droppable from vertex/FORM |
| F7 | `±∞` / `kL_W` undefined in 3-in-plane |

Nits 1–4 do not block clearance by themselves. Nothing in F1–F7 is “it would be wrong on a different input.”
