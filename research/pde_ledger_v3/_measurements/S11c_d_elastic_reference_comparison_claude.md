**Verdict: COHERENT CONDITIONAL COMPARISON.** It needs the five corrections below. None of them requires a physical choice for the comparison itself.

## 1. Well-posedness as a limited comparator

The paper is a limited comparator. It says the following openly:
- B_* is a new datum.
- The curl norm is symbolic.
- There is no relaxation law.
- Nothing is minimized over B_*.
- The reduction from four dimensions to the brane (P6) is still pending.

Neither case changes positions. Both keep u_d, v, χ and the geometry fixed, and only B_* differs. That matches s11ca 1a, where x=X+u and the inertia ½ρ_br⁰|∂ₜu|² are untouched.

Neither case adds a reservoir or a sink. In the retained case, χ: 0→1 at fixed B costs μ|r_b|²/2 per volume. That cost is a conversion term inside μ_χ, which P8 defines as δF/δχ. It is funded by the existing bulk re-ordering, and A2 (conversion moves order between phases and never creates medium) is not violated.

Two hidden inputs remain, and the paper should name them:
- **Retained case.** B_old is a spatial field on the parcel, carried through the disordered interval. Conversion in A2 moves order between phases, so there is no obvious parcel identity for B_old to ride on in a flowing, disordered phase. The paper's "open carrier" caveat is correct. In practice, choosing the profile of r_b is choosing a pre-stress. It is a support datum in disguise unless it is derived.
- **Relaxed case.** "B_*=B_b at formation" needs a freeze rule. χ rises continuously, so "t_b" is a threshold choice. If B_* tracks B until some χ_c and then freezes, the stored energy at χ>χ_c depends on χ_c and on how fast B changes during the ramp. The relaxed limit is therefore the limit of infinitely fast relaxation until freeze, and the freeze point is a datum. It is not a relaxation *rate*, so the paper's "no relaxation rate" statement survives, but "formation datum" should be defined this way.

## 2. Formation energy and the common Hessian

The algebra is correct:
- e_R and e_N at formation are correct.
- ∂_B e_R = χμ r_b and ∂_B e_N = 0 are correct.
- The quadratic expansions in δB are correct.
- The B–B block is χμI, which becomes χμCᵀC in displacement gradients.

The restriction "fixed χ, reference and geometry" is stated correctly, and the paper does not claim equal spectra or stability.

**Correction 1.** With χ perturbed, the incremental Hessian is not common. Since e is linear in χ:
- ∂²e/∂χ² = 0.
- ∂²e/∂χ∂B = μ r. This is μ r_b in the retained case and 0 in the relaxed case.
- ∂²e/∂χ∂B_* = −μ r.

So the cross-coupling between order and shear perturbations differs between the two limits. The paper's closing remark that coupled order perturbations "remain relevant" should be made explicit. This coupling is the one place the two cases already differ at quadratic order.

## 3. First variation, boundary terms and work

**Checked and correct:**
- δχ coefficient: f_B′−κΔχ+μ|r|²/2, with boundary term κ∂_Nχ δχ. The shear part matches P8's μ_χ=δF/δχ.
- δB_* coefficient: −χμ r.
- P_{iα}=χμ C_{Aiα} r_A.
- Interior term −∂_iP_{iα}δu_α and cut load N_iP_{iα}.
- The partial rate D_t e = (μ/2)|r|²D_tχ + χμ r·D_tB − χμ r·D_tB_*, and the ∫e ∇·v term.
- Reset of ordered material: Δe = −χμ|r|²/2.
- Cost at χ=0: zero in this term only.

The sign and normalization are consistent. The paper correctly separates the birth datum at χ=0 from resetting ordered material. The first is free in this term, and the second releases or absorbs real energy.

**Correction 2: ∇χ in the divergence.**
- ∂_iP_{iα} = μ ∂_iχ C_{Aiα} r_A + χμ C_{Aiα}∂_i r_A.
- The first term is an interface body force proportional to ∇χ and r.
- It is nonzero at the front in the retained case, because r_b≠0. In the relaxed case it vanishes at formation.
- The paper should state it, because it is where an r_b-dependent force on the wall would first appear.
- Also, "fixed B_*" in the first variation means Eulerian-fixed in the linear patch. In a nonlinear map, B_* carried with the parcel moves under δu, so the variation formulas hold only in the linear patch.

**Correction 3: double counting.**
- The first term of D_t e_sh, (μ/2)|r|²D_tχ, is the shear part of P_order=∫μ_χ D_tχ in P8.
- The second term is the mechanical stress power. It belongs to the kinetic and stress balance.
- The third term, the reference-work term, is a new state change.
- A total ledger must count each once. The paper states that the adjunct is not added on top, but it should say outright that the first term *is* the shear piece of P_order. Otherwise a reader could add them.
- A nonzero D_tB_* must have a named destination, or it is a hidden sink. The paper says as much for resetting.

## 4. Separation from existing structure

The separation is correct:
- B_* is not the bulk modulus B, a microrotation, the retired polar-P construction, or a renamed u_d. This matches native 3.2 (u_a is physical displacement, and a micropolar field would be new parent-action content) and 5.7 (the pathA_35 exclusion).
- The paper does not use a total-mass sink or a fixed-wall solve.
- Inertia stays homogeneous per 5.3.

**Correction 4: torque and frame (R-S8-04 and R-S8-05).**
- For a curl-type C, the stress is antisymmetric, and the torque density is ε_{jiα}P_{iα} = 2χμ r_j.
- The retained case therefore carries a static formation torque 2χμ r_b that needs a couple-stress or spin carrier (R-S8-04). The relaxed case has none at formation, only 2χμ δB later.
- This is a real discriminator, but it is the original obligation made static. It is not a new one.
- For a rigid rotation, B=2ω. The relaxed rule absorbs formation-time rotation into B_*, which ties the reference frame to the formation history. The retained rule ties it to the pre-disordering history.
- Neither answers R-S8-05. Under later rigid rotation, both leave r nonzero.
- The paper should say that the relaxed case relocates the preferred frame to the formation epoch and does not remove it.

**Correction 5: orientation and energy parity.**
- Under R_w (Appendix C), the tangent displacement is even and its w-profile is reflected. If B_* is reflected together with B, then e is a scalar, even under both R_w and C_s.
- A positive scalar stored-energy comparison therefore cannot carry attraction or repulsion between like and unlike orientations.
- The paper makes the same point in prose. It should add that any orientation-odd effect must come from the core coupling (the odd s-weighted couplings), and it must not be pushed into B_*.
- If B_* is not reflected with the state, it silently becomes a preferred-orientation datum. That is the R-S8-05 failure again.

Beyond those, the "reference field is not necessarily a curl" caveat is correct. If B_* has no compatible displacement, defects or extra state appear, and none is inserted.

## 5. Next discriminator

No physical choice is needed for the comparison as a paper object. The two limits can be compared as stated, once Corrections 1–5 are made.

Adopting either limit as a law needs one datum, and it is the same for both: **the formation rule for B_***. The smallest next item is a conditional paper test and is not a solver:
1. Fix the linearization patch.
2. Specify B_old (retained case) and B_b together with the freeze criterion in χ (relaxed case) as given data. Do not tune them to a force sign.
3. Evaluate, symbolically, three quantities for each case: the static torque 2χμ r, the interface force μ∇χ·C r, and the χ–B cross-coupling μ r.

Then check which of these a core with no couple carrier can accept. If R-S8-04 is answered "no carrier", the retained case needs r_b=0 or a divergence-free torque, and it collapses to the relaxed case. That would decide the question *for the existing substrate* without any choice about electric sign.

The corresponding question for the electric mechanism is separate. It asks whether the core coupling, odd in s, converts the surviving loads into a normal mouth stress that is orientation-odd, with B_* reflected consistently. Paper consistency is established without that.

## Source coverage and limits

**Inspected:**
- proposal.md and guide.txt in full.
- native-interpretation: the heading outline, plus §3.2, §5.3, §5.7 and Appendix C (C.1–C.3) in full.
- stage006-material-record: P1–P13 in full, and A2 (A1 and the start of A3 partly).
- s11ca-physical-definitions §1a and §1c.
- substrate-requirements: R-S8-04 and R-S8-05 in full. R-S1-02 only through its first lines, which show the on-failure statement. The rest of that entry is not read.

**Not inspected:**
- v3-step-plan (the banked re-ordering and S12), so the no-sink statement rests on the paper's and guide's wording and on A2.
- stage006-material-source.py and s11ca-geometry-source.py.
- stage006 corrections and errata.
- s11ca §3 (the face maps).
- native-interpretation §14.3.
- review-prompt.md and packet-index.json.

No code was run, and I did not recompute the algebra beyond checking it by hand on paper. The torque identity assumes C is the three-dimensional curl; the four-dimensional norm stays symbolic, so that result is conditional on the component convention.