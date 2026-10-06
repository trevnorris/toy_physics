# S9b — what brane light needs in order to bend and be delayed like GR (question spec, v8)

**Authors:** Codex (v5–v6), with v7–v8 edits by Claude (orchestrator). **Status:** v8, cleared for the build
on 2026-10-06.

## Why

S9 states what light needs from the medium in order to exist and to stay confined. It does not state what light
needs in order to be deflected and delayed by a mass. This step states that requirement as conditions on the
local wave equation that brane light sees near a mass. It does not derive that equation from the substrate:
the nonuniform brane operator is unfinished S11c work, and the substrate comes at the knit.

## Governing object (supplied; the build cannot test it)

Light is the transverse branch of the brane's in-plane displacement `u`, and `u` is the material displacement
of the stuff whose density is `ρ_br` (`steps/S11_stray_longitudinal.md:32–35`; S10). Near the mass, at
leading eikonal order, that branch obeys the supplied local dispersion relation

```
(ω − V^i k_i)² = c_γ(x)² · g^{ij}(x) k_i k_j ,      g_ij = δ_ij + ∂_iξ_w ∂_jξ_w .
```

Each piece is a supplied identification. Flag any result that depends on one.
- **Local light speed `c_γ(x)`.**

  ```
  c_γ(x)² ≡ μ_⊥(x)/ρ_br(x) .
  ```

  Here `μ_⊥` is the basis-invariant transverse stiffness (`steps/S11b_interface_coupling_law.md:74–87`);
  S9's basis writes it as `μ_R` (`steps/S9_light_requires_shear.md:78–79`). Both `μ_⊥(x)` and `ρ_br(x)`
  remain live. The same `ρ_br(x)` appears in the mass balance below. There is one isotropic speed, the same
  for both polarizations, measured relative to the shear-carrying material.
- **Advection.** `V(x)` is the in-plane velocity of the material that `u` displaces. Its kinetic symbol is
  `(ω − V^i k_i)²`: the coefficient of `ω²` is `1`, the coefficient of `ω` is `−2V^i k_i`, and the `kk`
  slot contains `V^iV^j k_i k_j`.
- **Steady brane mass balance.** The supplied balance is

  ```
  ∇·(ρ_br(x)V(x)) = −j_n(x) .
  ```

  `ρ_br(x)` and `j_n(x)` are live radial profiles. `j_n` is the brane's normal exchange with the bulk and is
  owned by the gravity sector or S12. S11b's uniform background normal drain `v_dr` is a different object
  (`directives/S11b_SHARED_PHYSICS.md:99–111`); the relation between `j_n` and `v_dr` is open.
- **Embedding.** `ξ_w(x)` is the brane's displacement into the bulk direction, as a length. The ledger's `h`
  is dimensionless, with `ξ_w = ℓh` (`directives/S11b_SHARED_PHYSICS.md:95–96`). `g^{ij}` is the inverse of
  the induced spatial metric `g_ij` above.
- **Anchor.** With `V = 0`, `ξ_w = 0` and constant coefficients, the relation reduces to S11b's uniform
  transverse branch, `ρ_br⁰ω² = μ_⊥k²` (`steps/S11b_interface_coupling_law.md:74–87`).
- **Background anchoring of `c_γ`.** The steady object in this step is the LAB_HELD profile `c_γ(x)`
  (`directives/S11c_a_SHARED_PHYSICS.md:232–242`). The S11c-a MATERIAL_ADVECTED anchoring is not used:
  for a steady radial profile, `D_t c_γ = 0` and `∂_t c_γ = 0` give `V_r ∂_r c_γ = 0`, inconsistent with
  keeping both a radial drain and a radially varying `c_γ` live.

**Outside this step.** These were narrowed out with the user's approval, and are recorded as open:
- direction-dependent (radial versus tangential) stiffness;
- coupling to thickness or bulk fields;
- polarization-dependent propagation;
- mixed `ωk` content other than advection.

**Reference speed.** `c₀` is the asymptotic value of `c_γ`, identified with the measured light speed. The bulk
sound speed `c_s` is separate and is not identified with `c₀`. Local light speed is not an observable here,
because rulers and clocks are made of the same medium. Only the far-field observables below are compared.

**Order.**
- **Eikonal.** The retained object is the dispersion relation above, with its position-dependent
  coefficients, and its rays. Excluded from this step's claim and not computed: the explicit subprincipal
  terms of the underlying operator, which affect amplitude and polarization transport.
- **Smallness.** Define `δ(x) ≡ c_γ(x)/c₀ − 1`. The retained set is every monomial

  ```
  δ^a (V/c₀)^b ((∂ξ_w)²)^c,      0 ≤ a ≤ 1,  0 ≤ b ≤ 2,  0 ≤ c ≤ 1,
  ```

  with nonnegative integer `a`, `b` and `c`. Each retained monomial is printed separately, and nothing outside
  this set is computed.
- **Order counting** (supplied; the build cannot test it; owned by the gravity sector or S12). In the far zone,
  let `ε(r) ≡ GM/(c₀² r)`, with `GM` as supplied below. The supplied counting is

  ```
  δ = O(ε),      (∂ξ_w)² = O(ε),      V/c₀ = O(ε^{1/2})   ⇒   (V/c₀)² = O(ε) .
  ```

  Under it, every other monomial in the retained box is `o(ε)`. Within these bounds, `δ`, `V` and `ξ_w` stay
  free radial profiles. Part B compares at first order in `ε`. For each Part B observable, it sums the
  contributions of the grades `δ`, `V/c₀`, `(V/c₀)²` and `(∂ξ_w)²` and subtracts the reference. The `V/c₀`
  contribution is computed and printed at its own order, `ε^{1/2}`, not assumed. Every other retained monomial
  is printed as its own higher-order object and is not compared with the first-order references.
- **Profiles.** Prefer general radial functions. If an engine restricts itself to a family, it says so and keeps
  every exponent and coefficient symbolic.
- **Inherited limits and freezes.** S9 took the sharp zero-width-sheet limit, no dissipation,
  frequency-independent moduli, the continuum limit and amplitude `→ 0`
  (`steps/S9_light_requires_shear.md:349–350`). It also took two distinct background limits
  (`directives/S9_wl_rebuild_directive.md:359`; `steps/S10_two_transverse_photons.md:816–819`):
  - S9's in-plane background flow `v₀ → 0` removed all convective terms. This step lifts that freeze by
    keeping the in-plane material velocity `V` live.
  - S9's background-strain freeze `strain → 0` is named separately. This step lifts its isotropic part by
    keeping it live in `c_γ(x)`; its direction-dependent part is the anisotropic stiffness excluded above.

  S11b separately froze its normal background drain `v_dr`: the correction recorded there as
  `O(v₀|q_n|/ω)` is uncarried (`steps/S11b_interface_coupling_law.md:158–164`). That correction remains
  outside this step, the `v_dr` freeze is open, and `v_dr` is distinct from both the live in-plane `V` and
  the live normal-exchange profile `j_n`. The effect of lifting every other inherited limit remains outside
  this step's claim.

## The observables

**Quantifier.** Each Part B comparison is required to hold for **every** far-zone `b`, not only at selected
`b`.

**Setting.** One isolated, spherically symmetric mass at rest, with its drain flowing: steady, not frozen,
with `V` live. Far field; linear waves. Time is the lab time of the brane's far-field rest frame.
1. **Δθ(b):** the total turning angle of a full flyby with impact parameter `b`.
2. **Round-trip (radar) excess time.** An emitter at distance `Z_E` on one side of the mass and a reflector
   at `Z_R` on the other, along the line. Subtract the flat round-trip time. Here
   `r_E = √(b² + Z_E²)` and `r_R = √(b² + Z_R²)`. Part A prints the full excess time; the Part B comparison
   uses only the coefficient of `ln(1/b²)` in the regime `Z_E, Z_R ≫ b`, for every profile, including tails
   other than `1/r`. The radar claim is limited to this logarithmic component.
3. **The two one-way excess times** between the same endpoints, and their **nonreciprocal part** (half the
   difference). Print whether that part depends on the path or only on the endpoints.

**Branch existence.** First print both the condition under which the branch is real and propagating everywhere
along the ray and the condition under which a ray can traverse the full flyby path and both legs of the round
trip in the required directions. Outside either condition, report the branch type (growing, decaying, absent,
or unable to traverse in a required direction), print `NOT_ESTABLISHED` for the observables, and do not compute
them.

## Supplied references (inputs the build cannot test)

- **`GM`** (owned by the gravity sector). The far-field Newtonian mass parameter, as measured by the orbits of
  slowly moving test matter. It is an independent symbol, not identified with any profile amplitude. Every
  Part B and Part C condition is stated relative to it.
- **GR reference, in PPN form with `c₀`:**
  - full-flyby deflection: `Δθ = (1+γ)·2GM/(b c₀²)`;
  - for `Z_E, Z_R ≫ b`, the round-trip logarithmic term is
    `2(1+γ)(GM/c₀³)·ln(4 r_E r_R/b²)`. Part B uses only its coefficient of `ln(1/b²)`;
  - `γ = 1` in GR.
- **Part C bulk:**
  - `P = Kρ^n`, with `n` symbolic, `ρ` the bulk number density and `m` the particle mass (v3 convention). So
    `c_s² = nKρ^(n−1)/m`.
  - Define the fractional bulk-density change at the brane by
    `f(x) ≡ δρ(x)/ρ₀ = ρ(x)/ρ₀ − 1`. It is a symbolic radial profile owned by the gravity sector or S12.
    `ρ₀` is the asymptotic bulk number density, `c_s0` is the corresponding asymptotic sound speed, and Part C
    keeps first order in `f`.

## Deliverables (printed by each engine; no conclusions in the scripts)

- **Part A.** The branch-existence and path-traversal conditions. Then, for the LAB_HELD profile, as expressions
  or functionals of `δ`, `V` and `ξ_w`, with each retained monomial printed separately:
  - `Δθ`;
  - the round-trip excess time;
  - the two one-way excess times and their nonreciprocal part.
- **Part B.**
  - An effective `γ` from `Δθ`, and one from the coefficient of `ln(1/b²)` in the round-trip excess time for
    `Z_E, Z_R ≫ b`.
  - Their residuals against the references.
  - For each observable, the condition on the profiles, relative to `GM`, that sets the first-order sum (see
    "Order counting") minus its reference to zero. Higher-order monomials are printed but not set against the
    references.
  - For each condition, the `j_n` it implies through `∇·(ρ_br V) = −j_n`, with `ρ_br` symbolic.
  - The difference between the two `γ`s, printed as an object.
- **Part C.** The Part B conditions rewritten, to first order in `f`, under three supplied responses of `c_γ`
  to the local bulk density, with `V` and `ξ_w` left live in each:
  - `c_γ(x) ≡ c₀` (`δ(x) ≡ 0`);
  - pointwise fixed ratio: `c_γ(x)/c_s(x) = c₀/c_s0`;
  - symbolic power response: `c_γ(x)/c₀ = (ρ(x)/ρ₀)^s`, with `s` symbolic.

  Here `ρ(x)` is the bulk number density at the brane and `ρ₀` is its asymptotic value. For each condition,
  print the `j_n` it implies through the supplied mass balance, with `ρ_br` symbolic. Print where `n` enters,
  if it enters at all.

## Engines, review, scope

- **Engines.** SymPy, plus a blind Wolfram engine that imports nothing. No Lean (CLAUDE.md L5).
- **Build review.** Codex-written, so a fresh Claude agent and Grok, each with a mandatory FORM ablation.
- **Model point.** As in "Setting": leading eikonal with the retained multigraded set above. Results do not
  transfer to:
  - the strong field;
  - the throat mouth or interior;
  - a moving or rotating mass;
  - a drain that changes while light crosses;
  - polarization transport;
  - anisotropic or coupled branches.
- **Deferred to the build** (implementation, not new physics; the build directive owns it):
  - the symbolic handling of the every-`b` requirement;
  - the mass balance on the induced metric. Under the supplied counting, its difference from the flat form is a
    relative `O(ε)` correction to the implied `j_n`.
- **Stop and report**, without choosing, when any of these happens:
  - a second method failure;
  - a sub-problem this spec does not name;
  - a premise this spec does not supply.
- **The step record** interprets the results.
