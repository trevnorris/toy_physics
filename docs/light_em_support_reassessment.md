# Shared light and electric support reassessment

Status: documentary reassessment and proposed research direction, following the user's October 4, 2026 request. No new physical law is selected, no scientific calculation is performed, and no independent review or runtime clearance is claimed.

The next useful step is to specify how the same material system responds to light and to an oriented throat. Completing electromagnetism is not a prerequisite for finishing the linear light analysis. A physical leakage factor does, however, require a declared support response and its energy exchange. Repeated checks of the existing operator cannot provide that missing law.

New leakage workers and implementation of the latest selected centre-drive method are parked during this reassessment. The earlier mixed numerical recovery remains parked. Existing calculations, reports, failures, inputs, physical parameters and source files retain their original status.

## Interpretation and evidence

The September accepted [native interpretation](native_light_em_and_vortex_throat_interpretation.md) governs vocabulary and cross-sector interpretation; current calculation records govern computed results. Its sections 3.2, 3.6, 7 and 14 are central here. The July [charge-attribute requirements](em_charge_attribute_requirements.md) describe an earlier additional-sector hypothesis. Their own sections 8 and 9 distinguish conditional spin-ice consistency from a derivation in the existing medium. That route is not silently added to the native material model.

The [v3 charter](../research/pde_ledger_v3/CHARTER.md#two-halves) allows an explicit postulated effective medium and asks whether the sectors' requirements can coexist. A postulate is an input to test, not a derived microscopic mechanism. The newer interpretation also requires the throat's effective mouth condition eventually to be realized by its core dynamics.

Useful saved results remain:

- Native transverse end flux (`archive/pre-cleanup-2026-10-04:research/pde_ledger_v3/_measurements/S11c_d_first_order_transverse_flux_continue_result_record.json`): the first-order flux correction K1 vanishes for the saved incident doublet, with the native end-current correction included. This removes one obstruction to the proposed expansion; it is not a leakage value.
- Receiving regularity (`archive/pre-cleanup-2026-10-04:research/pde_ledger_v3/_measurements/S11c_d_first_order_receiving_regular_continue2_result_record.json`): the three-field receiving block has the recorded real-axis domain and pole exclusions at the fixed inputs, and the opposite-leg cross-current vanishes. Full five-field transverse poles remain; neither a full field nor far-bulk power was computed.
- Background force comparison (`archive/pre-cleanup-2026-10-04:research/pde_ledger_v3/_measurements/S11c_d_background_force_pairing_record.json`): finite formal comparisons passed, but the physical background traction measure, complete force pairing, total centre load and support response remain unresolved.
- Latest centre-drive method assessment (`archive/pre-cleanup-2026-10-04:research/pde_ledger_v3/_measurements/S11c_d_centre_drive_method_review_record.json`): a limited selected-source compatibility method was assessed. No worker or runtime result exists for that method. A zero selected drive would not supply a centre equation or a support law.

These facts justify preserving the work while changing the immediate priority. They do not establish an incompatibility between light and electricity.

## What can move

| Motion or property | Physical meaning | Consequence for this reassessment |
|---|---|---|
| Tangential shear, u_T | Material moves within the brane in the transverse light sector | Keep the supplied light modes and their possible conversion channels. |
| Tangential compression, u_L | Material compresses within the brane | It is a distinct physical mode, not the electric normal displacement. |
| Thickness motion, e_W | The two faces move oppositely in laboratory w; thickness changes | It is retained in the selected five-field light calculation. |
| Centre motion, zeta_c | Both faces move together in laboratory w | The native geometry retains this independent coordinate; the selected five-field reduction does not supply its complete dynamics. |
| Electric mediator, h | The interpretation's normal brane displacement, with geometrical displacement ell times h | Its precise connection to the finite-thickness centre coordinate needs a source and normalization map. No equality with e_W is assumed. |
| Throat orientation, s | A defect opens toward +w or -w | This labels candidate charge sign; it does not itself determine force, magnitude or conservation. |
| Trapped support mode | Internal wave energy may exert stress that helps maintain the throat | Its energy and response must be distinguished from the incident probe wave. |

The thickness and centre distinction is explicit in the existing [native geometry source](../research/pde_ledger_v3/scripts/S11c_a_interface_geometry_sympy_audit.py): `dof_fields` assigns opposite face displacements to DELTA_W and a common displacement to ZETA_C. This was source inspection only; no constructor was called.

A selected light input failing to drive the centre would not show that the centre is absent, nor that an oriented throat cannot source normal displacement. Conversely, retaining a centre coordinate does not prove a healthy, localized electric mediator.

## Three different support questions

**Background support:** what holds the chosen nonuniform profile stationary, and how does that agency respond when waves arrive? LAB_HELD specifies how the background profile is sampled. It does not by its name specify a conservative support, a clamp, a controller or zero work.

**Throat support:** the existing physical picture uses outward trapped-wave stress against inward tension and bulk backpressure. The accepted interpretation treats a stable realization as unresolved. If an oscillating support mode maintains a stationary average shape, the complete branch may be periodic. It is not automatically a static energy minimum.

**Electric mouth response:** what does the throat maintain when another throat perturbs it: a displacement, a source, a flux, or a compliant combination? This response must be connected to the support mechanism. It cannot be inferred from the name of the light-test background.

The throat and brane are open subsystems. Internal support therefore does not imply isolation or passivity. Drain, return, phase conversion and reservoir exchange remain explicit where present; none is switched off by this reassessment.

## Joint physical requirements

The user's electric target is like-orientation repulsion and opposite-orientation attraction. Section 7.5 of the accepted interpretation records conditional exterior calculations with different signs for different mouth conditions: fixed value gives the desired sign, fixed source the opposite sign, and fixed-monopole or mixed cases depend on their data. Core and reservoir work are part of that distinction. The force sign cannot be selected from positive exterior stiffness alone.

For a proposed shared response, assess these requirements together:

1. The existing light sector retains positive-energy propagating modes in its claimed regime, with conversion and loss accounted for.
2. A localized normal response remains available over the intended electric range. Maintaining a local throat mouth is not the same operation as fixing normal displacement everywhere.
3. The same throat model produces both orientations and the target electric interaction sign, including its source and reservoir work. Its coupling to moving throats must remain a later check against the same transverse sector.
4. Every support force has an identified displacement or rate with which it exchanges work. Stable shape, correct force sign and small light leakage remain separate tests.

There is no new proof here that these requirements coexist. There is also no present basis for declaring that they conflict.

## Where wave energy can go

The survival observable counts both transmitted and reflected transverse flux. Its deficit can represent conversion into other brane modes, outward bulk radiation, dissipation, or energy retained in the throat. Support and bulk reservoirs can also supply or remove energy. These contributions must be separated before calling the deficit a positive intrinsic leakage factor.

The required accounting therefore includes transverse end flux, other material and bulk fluxes, interface memory storage and dissipation, chemical and mass-transfer work, changes in trapped-mode energy, and any work needed to maintain the background. A stationary mean profile does not establish that these exchanges vanish. This list specifies a balance to derive; it is not an asserted identity for the current retained equations.

The independent background grades retained so far also do not supply every physical term at second order along the chosen defect path. The missing 20/02 terms, corresponding end-current contributions and threshold remainder must be derived, supplied as model content or shown unnecessary for the specific observable. They cannot be dropped because K1 is zero.

## Candidate interpretations and recommendation

| Candidate | Physical picture | What it would let us claim |
|---|---|---|
| Finite stored support energy | A trapped wave contributes pressure; its energy can change when radiation or conversion occurs | A conditional lifetime or slowly evolving supported state, if that regime and its balance are established. Exact permanent support is not presumed. |
| Autonomous support with throughput | The same medium maintains a steady or periodic throat through internal dynamics and explicit bulk exchange | A response around an open steady state, including any energy drawn from the background. It is not automatically passive. |
| Externally maintained profile | A specified holding apparatus maintains the chosen test geometry | A useful supported light experiment, with apparatus work included; it does not establish an autonomous particle. |

Recommendation: formulate the next candidate as an autonomous material response with explicit energy storage and reservoir ports. Keep the finite-store case as a limiting question, and use an externally held experiment only with its conditional interpretation stated. This is a proposed study direction, not adoption of an unprovided support law.

The smallest useful next proposal is a short effective-core specification: the mouth and shape coordinates, the trapped energy that supplies stress, the response to a slow signed displacement and to a wave, and all energy/flux ports. Any effective compliance or coupling introduced must be labelled an input. Static sign and wave conversion must be tested with those same inputs. A desired fixed-value sign alone is insufficient grounds to impose fixed value.

Before another worker, this specification should answer one concrete question: can the proposed core maintain an oriented mouth while allowing the required exterior normal and transverse responses, with explicit work accounting? The source and normalization map connecting h to the finite-thickness coordinates belongs in that proposal. A full nonlinear throat solution remains a later realization test.

No new numerical computation, review submission, scheduler, validator or framework is part of this reassessment. Necessary scientific proposals still receive the standing Claude-only assessment and the existing guarded readiness process. Earlier next-action fields remain historical; they do not trigger centre or leakage work while this user-directed reassessment is active.
