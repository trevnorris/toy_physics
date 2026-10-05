# Oriented throat response and energy proposal

Status: new physics proposal for Claude-only assessment. No core model has been adopted or executed. This extends the user-approved [shared light and electric reassessment](light_em_support_reassessment.md); the existing S11c equations, parameters and saved results remain unchanged. Centre-drive implementation, leakage workers and mixed numerical recovery remain parked.

The proposed next question is: can an internally maintained signed mouth produce the desired electric interaction while giving an explicit account of energy exchanged with waves and bulk flow? We should test the force and work interpretation before building a throat solver.

## Physical content of the proposal

Use one finite throat with a mouth, deformable core, trapped shear mode and named bulk exchange ports. Its opposite orientation is the reflected version of the same model, not a separately adjusted parameter set. Following the accepted [native interpretation](native_light_em_and_vortex_throat_interpretation.md), trapped-wave stress is a support candidate; a stable supported branch and its energy budget are still to be established.

The candidate state description is deliberately small:

| Variable | Meaning | Status |
|---|---|---|
| a | Mouth size or one retained shape coordinate | New effective coordinate; a solved geometry would have to justify the reduction. |
| x | Dimensionless signed mouth value h_A; physical displacement is ell times x | Uses the electric report's convention. It is not the thickness field e_W. |
| s = plus or minus one | Throat orientation | Supplied branch label; no topological protection or quantized magnitude follows. |
| Q and P | Amplitude coordinate and conjugate momentum of one candidate trapped support mode | New reduced description, not a derived native mode or either of the two incident light amplitudes. |
| Bulk and interface state | Pressure, chemical/mass exchange, conversion and return data | Explicit ports; no replenishment or zero work is assumed. |

The finite-thickness centre coordinate zeta_c and x need an actual geometric projection and normalization map. They must not be equated by their names. The support mode also differs from the incident test wave. A one-mode support trial does not remove either incident transverse polarization.

## A storage model to specify before a response law

As a candidate conservative storage sector, write

    H_storage = T_shape + U_material(a,x;s)
                + P^2/(2 M_tr) + M_tr Omega(a,x;s)^2 Q^2/2.

M_tr is positive, Q has length units, P has momentum units, and Omega is a positive inverse time on a putative trapped branch. These functions and the reduction to one oscillator are proposed inputs, not outputs of the existing model. T_shape and its inertia must be specified if shape dynamics are retained. U_material includes whichever local elastic and backpressure contributions are actually internal; imposed bulk pressure work must not also be counted there without its reservoir convention.

This is a template, not a closed Hamiltonian: the exterior fields, existing mouth coupling and explicit open ports must be joined once before forces are assigned. A throat that radiates or converts material needs their dynamics in addition to storage. No damping term is introduced solely to make it stable, and a mean stationary shape is not evidence of equilibrium.

Retain Q and P for the wave response. A slow, weakly exchanging harmonic mode can instead be described approximately by its action I = E_tr/Omega, not by assuming its energy and amplitude both stay fixed. The standard adiabatic result has hypotheses that have not been demonstrated for this throat ([MIT classical adiabatic invariant notes](https://ocw.mit.edu/courses/8-06-quantum-physics-iii-spring-2018/3d0d5959338277143edefca3103eeff9_qxBhW2DRnPg.pdf)). Such averaging would be a separately declared limit; it does not apply automatically to the existing omega3 probe. No quantization is inferred.

In physical terms, the wave can support the core only if changing the shape changes its energy in the appropriate direction under the declared state constraint. The direction and magnitude of that stress remain to be calculated. They are not obtained by calling the energy trapped.

## The first small test is about force and work

The existing [electric-sign calculation](../software/em_charge_attribute/puncture_deflection_electric_sign_result.md), Appendix B, is unusually useful here. Its positive fixed-value coefficient uses a conjugate functional that includes work maintaining the mouth value. Its explicit wrong-functional control gives the opposite sign from bare stored energy. An internal elastic preference for a signed mouth is therefore not automatically equivalent to that fixed-value calculation, even in a very stiff limit.

Proposed first comparator: at fixed a and fixed internal state, give each mouth a finite elastic preference

    U_bias(x;s) = kappa_m (x - s phi_*)^2 / 2,
    kappa_m > 0, phi_* > 0.

Both x and phi_* are dimensionless; kappa_m has energy units. These are two explicit additional inputs. The signed preference and the label s are assumptions of this comparator, not a derivation of orientation protection. This bias is only the local x dependence of a trial core storage energy, not a claim to have derived U_material or the trapped stress. It is defined relative to the same asymptotic brane reference as h_A; it does not postulate an external clamp on normal displacement everywhere.

For this comparator choose the report's ADD organization: retain the existing orientation source g and its original mouth self term, then add U_bias exactly once. This is a scoped trial completion, not selection of ADD for the native throat. Preserve the full original source convention. The existing source work must occur once in the total candidate energy or named reservoir work. If its origin cannot be assigned consistently, the result is unresolved rather than a force prediction.

Use the saved two-mouth response relation x_i = partial E_0/partial y_i, where E_0 is the exterior self-plus-interaction storage and y_i is the total conjugate mouth source. Restore the original normalization, self coefficients, coupling and source terms before any new comparison; no old exterior solve or matrix inversion is proposed for replay. The work of this new test would be to solve the NEW finite-compliance stationarity problem and differentiate the appropriate complete potential with respect to separation, at declared fixed core state. No new calculation has been made here.

The assessment must first establish whether the source assignment, variables held fixed and energy functional are sufficient for that question. Then a future bounded calculation could ask:

1. Does an isolated mouth have a nonzero locally stable signed state in this comparator? This is local mouth stability, not existence or stability of a throat.
2. What is the conservative orientation-dependent force for both same and opposite signs, after the mouth response and core energy are included? Keep orientation-even forces separate. No parameters are chosen after seeing the sign.
3. In the stiff limit, does this internal bias approach the report's fixed-value thermodynamic problem, or only its geometric boundary value? Display the missing or additional work term if they differ.
4. Can any claimed desired sign be attributed to a named source of work, rather than to an omitted energy term?

This test can reject the simple elastic implementation without rejecting every vortex-core mechanism. A wrong sign is evidence, not authority to insert a compensating reservoir or change the potential. If an active or flux-maintained mouth is needed, its specific law becomes a new physical proposal. Restoring a boundary label alone cannot demonstrate it.

## Bulk exchange and light response

The proposed full model has a finite control volume containing the core, trapped mode, interface storage and the near exterior field. Its energy accounting must identify input from probe waves and bulk/chemical channels, outgoing transverse and other material waves, far-bulk radiation, exported heat, and any change in stored energy. Internal conversion cancels between components when both are inside the volume. Chemical and mechanical work must not be added a second time if already carried by the chosen bulk energy flux. These are required conventions for a future balance, not an identity derived for the existing equations.

Specify the baseline and perturbation of each exchange. A steady input balancing steady radiation does not determine how that input responds to a passing wave. Conversely, zero net background power does not establish zero incremental work. Replenishment, if present, must come from a named material process and its response; no freely adjustable controller may be added after a loss result.

The slow two-throat force test and the finite-frequency light response must use the same core couplings and declared reservoir law. Causal memory or additional state variables may be required; no instantaneous spring is presumed adequate at omega3. The existing interface memory and permeability are preserved, with no identification between their coefficients and the new trial mouth stiffness.

Transmitted and reflected transverse flux both count as survival. The support-mode energy loss, bulk radiation fraction and transverse-survival deficit are distinct observables until a complete balance connects them. The stored K1 zero and scoped receiving regularity do not supply that connection. Physical 20/02 terms, end currents, far-bulk asymptotics and the threshold remainder remain open in the existing leakage track.

## Proposed order and decisions

The immediate deliverable is assessment of the finite elastic comparator and its energy convention. It is the smallest test of the tempting inference that an internally stiff signed mouth automatically gives the known desired force sign. It does not need a new PDE solver, Fourier integral, particle simulation or general verification framework.

If that comparator is well posed, a separately reviewed and guarded finite test can classify it. If its work convention is insufficient, specify the exact missing port or core law before computing. If its force sign fails, preserve the failure and assess a different physically motivated mechanism explicitly; do not select one by tuning the sign. No worker launch follows automatically from a method verdict.

Only a supported realization should proceed to the shape/trapped-wave balance and then the wave response. These require actual U_material, Omega and reservoir laws, with a finite list of parameters and validity assumptions. Their present placeholders are not an executable model, and arbitrary function freedom cannot count as a successful fit. A postulated effective realization would remain conditional pending its connection to a native throat branch.

This proposal introduces no localized Maxwell or spin-ice sector, no identification of circulation with charge, and no moving-current law. Magnetism remains a later test of translating the same supported branch through the transverse sector. No scope has been opened for a defect sweep, drain calibration or full nonlinear simulation.
