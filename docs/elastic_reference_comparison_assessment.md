# Elastic reference comparison findings

The user-authorized comparison is complete at the level of a conditional paper calculation. Claude's [literal report](../research/pde_ledger_v3/_measurements/S11c_d_elastic_reference_comparison_claude.md) says **COHERENT CONDITIONAL COMPARISON** and requests five qualifications. The [reviewed paper](elastic_reference_limit_comparison.md) remains frozen. This assessment accepts the checked algebra, adds the relevant mixed-order terms explicitly, and separates those findings from several unsupported physical conclusions in the report.

No material law has been selected. No scientific worker, restored expression, old function, numerical mode, power integral or leakage calculation was run. The expressions below are paper consequences of the declared comparator, not native runtime certificates.

## What differs physically

For the conditional energy `e = chi*mu*|B-B_star|^2/2`, let `r = B-B_star`. At formation, the retained-reference case has `r=r_b`; the newly relaxed case has `r=0`. The comparison holds current geometry, physical displacement, velocity, density, order and coefficient fixed. It holds reference data fixed during the independent shear variation, rather than minimizing over them.

| Quantity from this partial energy | Retained reference at formation | Newly relaxed reference at formation |
|---|---|---|
| Stored shear energy | `chi*mu*|r_b|^2/2` | `0` |
| Conjugate shear load `partial_B e` | `chi*mu*r_b` | `0` |
| Incremental shear block `partial_B partial_B e` | `chi*mu*I` | `chi*mu*I` |
| Mixed order–shear block `partial_chi partial_B e` | `mu*r_b` | `0` |
| Shear contribution to the order derivative | `mu*|r_b|^2/2` | `0` |

Thus forming unstressed material does not erase its later incremental stiffness. Retaining a deformed reference also supplies pre-existing loading and a possible order–shear coupling. This is not a finding that either material supports a stable throat or has equal full wave spectra.

The mixed term is worth making explicit. For independent small perturbations `eta=delta chi` and `b=delta B`, with reference fixed, the quadratic part of this shear contribution is

    e^(2) = (chi*mu/2) |b|^2 + mu*eta*(r dot b).

There is no pure `eta^2` term from this contribution because it is linear in `chi`. The original order-well and gradient energy supply other terms; their full stability and dynamics have not been evaluated here. At formation the second term is present in the retained case only when `r_b dot b` is nonzero. It is zero in the newly relaxed case, but may develop after subsequent deformation. This is a conditional coupling coefficient, not a leakage probability, dissipative rate or guaranteed overlap with a particular light mode.

The original paper's first variation and material shear-energy rate are accepted in their stated independent-field, fixed-chart scope. The term `(mu/2)|r|^2 D_t chi` is already the shear contribution to `mu_chi D_t chi`: it must be counted once in any total work balance, not added again as a separate conversion supply.

## Boundary loading has a conditional direction

For the same symbolic curl map `C` as the paper, the conjugate gradient load is

    P_iα = chi*mu*sum_A C_Aiα*r_A.
    partial_i P_iα = mu*sum_Ai (partial_i chi)*C_Aiα*r_A
                     + chi*mu*sum_Ai C_Aiα*partial_i r_A.

The first term exposes how an order gradient can enter this partial displacement load. The report calls it nonzero whenever `r_b` is nonzero. That conclusion is too strong: its contraction with the normal gradient and the relevant polarization can vanish. Both terms, the cut contribution and the actual field/geometry map matter. The sign of the negative variational derivative must also be distinguished from `partial_i P_iα` itself.

These are contributions to a specified displacement variation. They have not been joined to the complete physical Cauchy stress, interface force, native normal mouth displacement or electric work. A nonzero `r_b` is therefore neither proof of outward throat support nor evidence of repulsion.

## Which review claims are not accepted

**Re-ordering is not a demonstrated energy source.** The report says that the retained case's additional energy is “funded by the existing bulk re-ordering.” The conserved mass/order balances establish the meaning and rate of conversion, not its available energy. The bare order well has equal pure-phase minima at fixed density. Other energies or work could supply the storage change, but their actual balance remains to be shown. The original paper explicitly kept this gap open; preserve it.

**Material parcels and memory are different questions.** Total continuity and a prescribed smooth material velocity can define material paths even when order changes. Conversion does not by itself remove parcel identity. What is unsupplied is a physical state capable of retaining the proposed elastic reference through disordering, together with its transport and work. Neither same-medium ontology nor a material label proves such memory exists.

**A prescribed birth datum is not automatically an infinite-relaxation law.** The comparator specifies an initialization event and subsequent reference history; it does not require `B_star` to track `B` throughout a conversion ramp. A continuous constitutive implementation would need an event/formation rule and ramp behavior. There is no derived order threshold `chi_c`, and none should be inserted as if the original model supplied it. The formation event may be treated as external test data for a conditional comparison, with its physical selection left open. If initialization is performed at nonzero order, its reference work must be kept.

**Fixed material reference and fixed spatial field are not interchangeable.** The first variation treats `B_star(x)` as an independent field held fixed at spatial coordinates. An actual virtual material motion of an inhomogeneous carried reference generally induces a spatial change of that field. This issue can already matter in linear perturbations about an inhomogeneous background; it is not restricted to a future nonlinear correction. The partial variation is valid as stated, but physical mouth-force use needs the induced variation derived from the selected transport and tensor convention.

**The three-dimensional torque formula is not a native four-dimensional result.** Claude discloses that its `2 chi mu r` axial-vector formula assumes the three-dimensional curl. The paper intentionally retained the original four-dimensional curl convention symbolically. That formula therefore cannot be inserted into this comparison as a native scalar/vector certificate. Moreover, the conjugate cut load has not yet been proved to be the complete physical stress. The actual angular-momentum and rotational-reference obligations in the [native interpretation, §5.7](native_light_em_and_vortex_throat_interpretation.md) and [R-S8-04/05](../research/pde_ledger_v3/SUBSTRATE_REQUIREMENTS.md) remain open.

The report's suggestion that “no couple carrier” would select the relaxed case is also not established. Being unstrained at formation does not remove the rotational-stress obligation when later waves deform it. A locally balanced full angular-momentum law cannot be replaced by calling an unjoined torque “divergence-free.” No general rejection or selection of either reference limit follows.

**Scalar energy does not forbid orientation-dependent forces.** The report overstates what joint reflection implies. A scalar interaction can depend on relative orientation: the schematic invariant `s_1*s_2*V(R)` is unchanged when both signs flip, but distinguishes like and unlike pairs. This is a symmetry counterexample, not a proposed new interaction law. Also, the native Appendix C explicitly says the complete orientation map `C_s` is not assumed identical to spatial reflection `R_w`. The comparison has no derived source/core map and therefore predicts no force sign; it does not establish that reference-mediated orientation dependence is forbidden. Reference data must transform with the physical state, rather than being tuned afterward to produce the desired sign.

## The smallest useful next test

The comparison identifies a more precise physical target than choosing a memory law by intuition: **whether the retained residual shear has an overlap with an actual order-changing deformation of the throat**. For a specified patch, independent order variation `eta` and displacement variation `xi`, the conditional mixed work is

    M[eta,xi] = Integral(mu*eta*r_b dot (C xi) d^4X).

Its vanishing or sign depends on profiles and admissible variations; it is not fixed by `r_b != 0`. The newly relaxed formation snapshot sets this partial contribution to zero, but neither all later coupling nor the other physical work channels vanish as a consequence.

The next bounded preparation should obtain those variations from the existing native throat/brane maps, or state exactly which core profile and reference history are still missing. It should distinguish independent field variation from a material shape variation and include the latter's induced reference transport. Keep the formation history supplied symbolically; do not invent a threshold, preload, relaxation rate or torque carrier merely to produce a number. If actual native map/profile data remain open, explain that concrete physical choice to the user before adopting a law.

This completes the requested comparison without selecting a winner. It needs no wording-only resubmission or general verification framework. A later substantive physical map or worker still requires its applicable assessment and readiness. Leakage workers, selected centre-drive implementation and mixed numerical recovery remain parked.

## Evidence preservation

The literal report and reviewed paper are unchanged. Review duration was 92.586 seconds, exit 0, with empty reviewer/coordinator stderr and no permission denials. All 11 packet/private/archive files, eight original sources, authority and helper hashes matched. Claude read the proposal/guide and the specified material-definition/reference passages; it did not inspect the original executable sources, S12 or native curved-face definitions. This was paper assessment, not runtime validation.

Preserved review metadata are `research/pde_ledger_v3/_measurements/S11c_d_elastic_reference_comparison_review_{sources,authority,readiness,launch}.json`, together with the literal report linked above. A separate completion-record JSON was not finalized before the scope pause; do not treat the formerly planned filename as an existing certificate. No scientific computation or replay occurred, and no new review or scientific job is started by this assessment.
