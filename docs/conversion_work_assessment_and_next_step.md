# Conversion work: assessment and the next physical choice

Status: source-only assessment of the [frozen paper](moving_interface_conversion_work_proposal.md) and its [literal Claude report](../research/pde_ledger_v3/_measurements/S11c_d_conversion_face_work_claude.md). The report says **COHERENT CONDITIONAL BRIDGE**, with partial source coverage. No scientific worker, symbolic restoration, old calculation replay or numerical result is involved. This note does not alter the reviewed paper.

## What this review settles

The moving-volume balances and the conditional travelling-front sign are coherent consequences of the supplied number/order balances. They preserve total material, boundary motion and accumulation. In the stated ordered-to-disordered normal coordinate, the steady-front relation is `j = Delta(K_n) - Integral(rho Gamma dz)`. Zero endpoint order-current makes positive re-ordering correspond to negative outward relative mass flux. Neither a steady solution nor its native curved-face mapping has been established. A DC conversion rate is not the harmonic S11b flux.

The corrected partial order-work formula and its integration by parts are also coherent. They do not supply a complete moving-volume energy law. The native face power has been used with its original pressure/chemical normalization; it has not been identified with order work. Finally, the `chi_B f_shear` gate supplies order-dependent storage, but it does not determine energy transfer to a trapped mode or an electric force.

These are assessed paper identities. They do not require another wording-only review, a generic verifier, or replay of an accepted calculation. The next uncertainty concerns the material model.

## Qualifications retained from the actual sources

- **Rate convention:** the paper uses the rate in the original conservation balance. It does not resolve the separately labelled P8 kinetic adjunct. If `Gamma_bal` denotes the balance rate and `Gamma_adj` the adjunct source, compatibility of the two displayed source equations would require `Gamma_bal - div(K)/rho = -M_chi mu_chi + Gamma_adj`. One cannot casually move relaxation into a transport current without a constitutive and boundary prescription. No such prescription has been adopted.
- **Front hypotheses:** steady-in-the-moving-frame includes density, not just order. Zero endpoint order-current is an explicit extra hypothesis. The review also asks for zero endpoint gradients; these are not necessary for the integrated identity and do not by themselves prove zero current without a current law. Do not add them as an unexplained compulsory condition.
- **Gradient boundary work:** for constant kappa, the fixed-domain first variation of the order-gradient term has boundary contribution `Integral_boundary(kappa * partial_N chi * delta chi dA)` as well as the volume derivative. Replacing a variation with a physical rate on a moving, advected domain needs the corresponding transport terms. The report's suggested boundary-power expression is not by itself a complete moving-domain energy law.
- **Native pressure is defined:** S11b already specifies bulk pressure and bulk mass density. The missing ingredient is a thermodynamic/kinematic map between different state variables and their conjugates, not a need to invent which density or pressure S11b meant. `mu_chi/rho` and `mu_theta/rho_br` have compatible units but different variational definitions.
- **Scope of reviewer coverage:** Claude inspected portions of stage006, S12, S11b shared physics and task_b2d. It did not independently inspect the full native interpretation or curved geometry. Full supplied bytes were available, but availability is not coverage. Its verdict is not a clearance of those uninspected constructions.

## The source follow-up narrows the missing physics

The September [native interpretation, sections 3.2 and 5.3](native_light_em_and_vortex_throat_interpretation.md) already declares in-plane `u` to be physical material displacement and supplies its homogeneous quadratic inertia. The original [S11c-a definition, sections 1a and 1c](../research/pde_ledger_v3/directives/S11c_a_SHARED_PHYSICS.md) is even more explicit:

    x(X,t) = X + u(X,t)
    T_u = (1/2) rho_br^0 |partial_t u|^2.

We therefore should not ask the user whether the established brane displacement is a material motion, or pretend that the homogeneous sector lacks an inertia. These are supplied effective inputs. What remains unsupplied by these records is how the material reference and elastic response extend through **new ordering**, and how the four-dimensional order model reduces to that effective brane kinetic law near a throat.

The same native interpretation, sections 5.7 and 14.3, explicitly leaves the reference structure of the rotational stiffness under investigation. The [substrate requirements](../research/pde_ledger_v3/SUBSTRATE_REQUIREMENTS.md), R-S1-02 and R-S8-04/05, record the order-dependent shear response, angular-momentum balance and physical rotational reference as open obligations. These are not closed merely by transporting coordinates. The older [MacCullagh differentiation essay](s11_maccullagh_differentiation.md) labels itself framing, not a physics authority; its broad survival prose does not override these later qualifications. The failed polar-P route remains failed and is not being reopened.

The limited source search has found an existing displacement map and an explicitly open reference obligation, but no supplied front-formation law in the inspected sources. This is a scoped finding, not a claim to have exhaustively searched every historical file.

## Physically meaningful alternatives to compare

The question is **what elastic reference newly ordered material acquires**, not whether return occurs. These are candidate classes, not completed constitutive laws or an exhaustive list:

| Candidate | Physical meaning | Consequence that a model must account for |
|---|---|---|
| Retained elastic reference | The same material carries a reference through its disordered interval, and re-ordering restores stiffness relative to that reference. | Old deformation may become stored elastic energy again. The state that carries the reference, its transport and the work needed when stiffness returns must be supplied. A shear-free bulk need not thereby carry a propagating shear wave, but neither does being shear-free establish such memory. |
| Newly relaxed reference | Newly ordered material adopts its current local configuration as its initial unstressed reference for the selected elastic measure. | New material initially adds no elastic loading from past deformation, while subsequent motion can strain it. Compatibility with adjacent ordered material, rotational reference, and any removed energy or momentum still require accounting. |
| Finite reference relaxation | Some memory survives but evolves during conversion over a stated time or length scale. | This adds a reference-evolution law and parameters. If elastic energy relaxes, its destination and possible feedback must appear in the work balance. Neither damping nor driving may be assigned without that law. |

Claude's shorthand option “reset u_d to zero” must not be adopted literally. The existing `u` enters the physical map `x=X+u`; changing it arbitrarily changes material positions. A relaxed **reference** must be distinguished from a physical displacement reset. Likewise, “conserve shear energy across the front” is a proposed outcome or constraint, not a sufficient transport/reference law; it leaves momentum, phase, other energy channels and interface work unresolved.

The recommended next preparation is a comparison of the **retained-reference and newly-relaxed limits**, leaving finite relaxation as a later model if needed. This avoids choosing an extra rate merely to obtain a preferred force. It is a recommendation to compare hypotheses, not adoption of either. The comparison should use the same prescribed geometry and physical displacement, state the reference data separately, and expose which boundary work and new fields each case needs. It must preserve the current MacCullagh/reference qualifications; neither case automatically supplies a physically admissible rotational substrate.

Before a numerical mode or core calculation, the small paper deliverable is the first variation of the existing order/shear energy for specified profiles, with its conjugate boundary tractions and explicit missing core/reference contributions. A full normal energy-momentum jump or physical mouth force cannot simply be inferred from that partial variation. Native geometry, inertia, compression, reference dynamics and throat terms must be joined before that stronger claim. No new force-sign test should silently add a reservoir after seeing its answer.

## Preservation and review record

The review completed in 82.984 seconds, exit 0, with empty reviewer and coordinator stderr. All 13 packet files, all 13 private copies, archive members, 10 original source files, authority and existing helper hashes matched. Literal and raw reports are preserved; packet SHA256 is `35746e5c4acba8b154af7efd736f2c809f32954d8e3d03552c4241e4ae3fa185`. Full receipts and the assessment record live under `research/pde_ledger_v3/_measurements/S11c_d_conversion_face_work_*`.

No new review or scientific execution is started by this assessment. Leakage workers, selected centre-drive implementation and mixed numerical recovery remain parked. Reflection still counts as survival. There is no trapped-mode solution, force sign, expansion rate or leakage factor from this paper. All previous physical inputs and results remain unchanged.
