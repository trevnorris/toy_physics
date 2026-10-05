I've finished reading the staged sources and am now writing the verdict.

**REQUIRES A PHYSICAL CHOICE**

The proposal is coherent as a conditional framework. The finite elastic comparator it centres on is a free, relaxing mouth, and on paper it lands in the report's FIXED_SOURCE class. The sign that class gives is the opposite of the one the user wants. Getting the fixed-value sign needs a named reservoir or a different core mechanism, and that choice belongs to the user.

## Q1. Variable distinctions

Mostly correct:
- **Mouth value and thickness.** x = h_A = ξ_w/ℓ is kept distinct from the thickness field e_W.
- **Centre coordinate.** ζ_c and x are not equated by name.
- **Trapped mode vs. incident light.** The trapped mode (Q,P) is kept apart from the incident transverse waves, and a one-mode trial does not remove either polarization.
- **Orientation.** s is a label, not protection or quantization.

Four gaps:
1. **In-plane position.** The brane translation coordinate X (so R = X₁−X₂) is never named. The force test needs it, and it differs from both ζ_c and x.
2. **Name clash.** The parent source amplitude Q_χ[r_Σ,s] in the report's mouth functional (`electric-sign-result.md`, App. A) shares a letter with the proposal's trapped-mode Q. They must not be identified, and the proposal does not say so.
3. **Open ports.** The oscillating mouth value x radiates exterior h-branon (and u_L) waves at ω3. This port is not named separately from "outgoing transverse and other material waves".
4. **Parity.** If the trapped mode is R_w-even, Ω(a,x;s) is even in x. Trapped-mode energy then renormalizes an even stiffness and cannot supply an orientation-odd force. Only an explicit s·x term can, and that term is a new input.

## Q2. Is the comparator a meaningful, fully specified test?

**Meaningful, but it tests a different thing than the proposal says.** Write U_bias = ½κx² − κφ* s x + const.
- The quadratic part adds to the existing mouth self-stiffness: ½k_m x² is already in the parent functional (`electric-sign-result.md` App. A, η(k_m h − g s)). Only K = k_m + κ matters, so κ_m is not independent.
- The linear part is an odd source term identical in form to the committed −g s x. It acts as G = g + κφ* (this is the proposal's single new odd "holding row").
- So the comparator is the report's free-mouth control (FIXED_SOURCE in the App. C table) with a larger committed source. It is not a fixed-value mouth.

**Missing data.** A force is not yet derivable from the text:
- The proposal says to differentiate "the appropriate complete potential" while using x = ∂E₀/∂y. That mixes the y-representation (E₀) with the x-representation (U_bias(x)). The potential must be written in one set of variables.
- The remaining datum is the ownership of the ADD source work. Is Q_χ a frozen sleeve datum, independent of x and R? The report's S_hold freezes Σ, which suggests yes. If so, the source term −g s x is an external fixed force with no explicit R dependence, and the declaration is cheap. If Q_χ is dynamical, it needs its own energy.
- The report's m_gg already folds in k_m (through z_g) and relaxes u_L (through z_b). Using it also requires declaring that the mouth u_L condition is natural/free.

No microscopic theory is needed beyond those declarations.

## Q3. The variational obstruction

**Definitions.**
- W(x;R) = ½ xᵀ S⁻¹(R) x. This is the exterior field energy at given mouth values, the Legendre transform of E₀(y) = ½ yᵀSy.
- S has off-diagonal element ε = m_gg s₁s₂/(4πR) with m_gg > 0, which the report has for D > 0. So the pair part of S⁻¹ is −ε/S_gg², negative.

**Comparator energy** (x-representation, everything internal):

Π(x₁,x₂;R) = W(x;R) + Σᵢ[½K xᵢ² − G sᵢ xᵢ].

**Hypotheses.**
- Quadratic exterior with the report's positive kernel.
- R-independent sources.
- Far field R ≫ a.
- Reflection covariance, so xᵢ* = sᵢ x* with x* ≠ 0.
- A conservative core whose energy depends on x only.
- No reservoir.

**Result.** By the envelope theorem, F_out = −∂_R Π|_{x*} = −∂_R W|_{x*}. This gives

A_comp = −m_gg (y*)² + O(1/R) corrections,

where y* = S⁻¹x* is the isolated net source.
- Same orientations attract and opposite orientations repel, for any κ > 0, φ* > 0 and g > 0.
- That is the J-class sign, −m_gg(j+g)² in the report.
- The result holds for any convex or double-well U(x), as long as x₁*x₂* has the sign of s₁s₂.
- Π is convex here, so the isolated mouth trivially has a unique, nonzero, locally stable state. This is not a throat-stability statement.

**Stiff limit.** As κ → ∞, x → sφ*. This matches the V class's geometric value, but the force is −∂_R W|_x. That is the report's "wrong-functional control: bare E₀ gives the negative" (App. B, V row). The spring energy y²/2κ vanishes while the reaction y stays finite.
- The term the comparator lacks is the reservoir work −y·x in Ω_V = E₀ − y h = −W.
- The test can display exactly this difference.

**What I did not verify.** I did not find in the report a derivation of why the clamp's reservoir energy is −y h rather than zero. The packet states Ω_V as the conjugate functional, and I take it as stated, without adjudicating or replacing it. A mechanically fixed x does no work under R-motion, which is why a conservative spring cannot borrow V's sign.

**Scope.**
- This is not a no-go for all vortex mechanisms.
- A core energy that depends on R or on the exterior gradient escapes it.
- So does a driven or flux-maintained mouth with a named supply.
- So does a different exterior kernel.
- Only the conservative, value-only, quadratic-exterior class is excluded.

## Q4. Trapped-mode storage and bulk accounting

Yes, it can be developed without hidden assumptions, but only if these points are handled.

**Fixed-state derivative.**
- ∂_R H at fixed (Q,P) is an oscillating force.
- Its fast-orbit average at fixed action is −I ∂_RΩ, with I = E/Ω. M_tr drops out.
- At fixed energy, the derivative is zero by construction. That is not the support force.
- "Fixed amplitude" and "fixed energy" give different pictures.

**Adiabatic caveat.**
- This needs a free or weakly exchanging mode.
- It needs Ω̇/Ω² ≪ 1 and no resonance.
- A probe at ω3 drives the mode and exchanges energy with it, so I is not conserved.
- Near ω3 ≈ Ω the averaging fails.
- The probe response therefore needs the driven problem, not I = const.

**Double counting.**
- U_bias is the x-dependence of U_material. Adding both counts it twice.
- Ω(a,x;s) also depends on x, so a populated mode adds ∂ₓ(IΩ) to the stationarity condition and shifts x*. Its sign is open.
- Imposed bulk pressure work must not be counted in both U_material and the port flux. The energy flux needs to be enthalpy-type, counted once.

**Hidden supply.** A frozen Q_χ is a source of infinite impedance. It does −g s δx work when the probe moves x, and that is a reservoir.
- It has to be declared, with its response.
- For the conservative comparator plus radiation to infinity, linear passivity is plausible. It does not extend to any driven or flux-maintained variant.

**Instant response.** The spring has no inertia, so a response at ω3 needs mouth inertia and exterior radiation damping. A static force test does not substitute for it.

## Q5. Smallest next step

**Paper step with declared symbolic inputs, no worker.** I derived this reasoning on paper and have not run it.
1. Fix the representation: write Π(x;R) as above.
2. Declare the inputs: S_gg, m_gg, K = k_m + κ, G = g + κφ*, frozen Q_χ, and a free u_L.
3. Obtain x*, y* and A_comp = −m_gg y*².
4. State the stiff limit and display Ω_V − Π_stiff = −y·x (summed over mouths).
5. Record that the comparator has no free parameter that can flip the sign. This is a legitimate wrong-sign outcome and should be preserved, not repaired.

**Genuine user choice, a new physical mechanism.** Neither option may be tuned afterwards:
- A named reservoir or drive maintaining the mouth, with its energy law.
- A core energy that depends on R or on the exterior gradient.

## Scope and what I inspected

**Read in full.**
- `proposal.md`
- `guide.txt`
- `electric-sign-result.md`

**Read in part.**
- `native-interpretation.md`: lines 1–1135, 1339–1380, 1717–1826 and 3569–3628. That covers §5.9, §7.5–7.6 and §14.5. I did not read §0.2 in detail, nor §3.2, §3.6 or §7.1–7.4 beyond the roadmap.
- `electric-sign-check.py`: grep hits plus lines 300–340 (the ADD/REPLACE amendment definitions).

**Not inspected.**
- `electric-sign-check.wl`
- `interface-geometry.py`
- `review-prompt.md`
- `packet-index.json`

I did not verify any numerical result. The sign claims above are paper derivations from the report's stated functionals.