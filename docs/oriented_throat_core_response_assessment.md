# Oriented throat response: concept review assessment

Status: conceptual assessment only. Claude's literal verdict is **REQUIRES A PHYSICAL CHOICE**. The [complete report](../research/pde_ledger_v3/_measurements/S11c_d_core_response_concept_claude.md) and the [reviewed proposal](oriented_throat_core_response_proposal.md) remain unchanged. No native core model, force, equilibrium, wave response or leakage has been calculated or accepted here. There is no scientific worker or READY gate.

The useful finding is that **a stiff signed mouth does not establish the desired electric-force sign**. A simple conservative mouth-displacement model has a conditional sign obstruction. That is a reason to examine its physical coupling and work, not a universal rejection of the vortex throat, and not proof that an external power supply is necessary.

The September [native interpretation](native_light_em_and_vortex_throat_interpretation.md), especially sections 3.2–3.6, 7 and Appendix C, remains the interpretation baseline. The goal remains like-orientation repulsion and opposite-orientation attraction from the same throat mechanism, with the same mechanism's light response and energy exchanges accounted for.

## What the report supports

The candidate needs a complete potential or balance law in one representation before differentiating with respect to separation. A core energy expressed in mouth displacement cannot simply be added to an exterior source functional without joining their variables and ownership. The existing mouth self term and source work must each enter once.

The physical brane positions of the two throats should be named **X₁, X₂**, with separation **R = |X₁ − X₂|**. They differ from normal mouth values **xᵢ = h_A,i**, thickness and centre displacement. Use **q_tr, p_tr** for a support oscillator in future drafts, reserving the original **Q_chi** for its original source attribute. These are notation clarifications, not changes to the frozen proposal.

The proposed bias expands into an even stiffness and a signed linear term. It therefore does not, by itself, demonstrate a dynamically realized fixed-value ensemble. Its stiff limit constrains geometry; it does not automatically realize the work convention of the old fixed-value calculation. Static response also does not supply a finite-frequency response at the existing omega3 probe.

## A conditional paper diagnosis, with signs kept separate

The following is an illustrative paper argument for the restricted comparator, not a restored result or an independently cleared native force calculation.

Suppose a conservative exterior has an **unsigned** two-mouth compliance

    S(R) = [[chi, c/R], [c/R, chi]],   chi > 0, c > 0,
    W(x,R) = (1/2) x^T S(R)^(-1) x.

Assume sufficiently large R that S is positive, and add local, R-independent core energies U_i(x_i;s_i). All internal source and self contributions must be included exactly once in this representation. Assume a smooth, stable isolated branch x_i = s_i x0 with x0 nonzero, and a differentiable large-R expansion. There are no additional separation-dependent core couplings or external work terms in this comparator.

To first order in 1/R, the exterior cross term is

    W_pair = -c x1 x2 / (chi^2 R).

Stationarity of each isolated local energy means the induced first-order change of x does not add another first-order energy term. Thus the leading relaxed pair energy is

    E_pair = -c (x0/chi)^2 s1 s2 / R.

With outward force defined as minus the separation derivative, this restricted model attracts like orientations and repels opposite ones. Making the local bias arbitrarily stiff does not reverse that leading sign under these assumptions.

Here chi and c are abstract comparator coefficients. This assessment does **not** identify them with a bare or dressed native coefficient, add k_m to an already dressed response, or assign any S11c parameter to them. That source join remains necessary. The argument also does not cover branch changes, general geometric or gradient couplings, dynamical internal-state changes, or flow-maintained states. An arbitrary double-well potential is not a global no-go theorem.

## Corrections and limits of the literal review

1. **Orientation is counted once.** The original [Python source](../software/em_charge_attribute/puncture_deflection_electric_sign_check.py), `fixed_value_coefficient`, uses off-diagonal `eps*mgg` and signed values `s1*phi`, `s2*phi`. Claude's stated off-diagonal also includes s1*s2 while retaining signed values. That double placement is inconsistent with the source. The paper diagnosis above uses an unsigned matrix.
2. **Bare and dressed self terms cannot be mixed.** The [electric report](../software/em_charge_attribute/puncture_deflection_electric_sign_result.md), Appendix A, already places k_m in z_g and the dressed kernel. Claude's proposal to use that response and add k_m again is not accepted without an explicit split. The new bias and old self term can be distinct physical inputs even if only their sum appears in a particular static parametrization.
3. **Joint reflection is not fixed-orientation evenness.** Appendix C maps both s and h under normal reflection, and may also map internal data. It does not establish that Omega(a,x;s) is even in x at fixed s. An s*x dependence can respect joint reflection. No such new coupling has been adopted merely because symmetry permits it.
4. **A frozen parameter does not prove an external reservoir.** A conservative material coupling may contain a fixed label and a term proportional to -g*s*x. The ownership of that energy must be stated. Freezing Q_chi alone does not establish an active supply, infinite impedance, or its response law.
5. **Possible wave ports are not calculated couplings.** Coupling of a moving mouth to h, u_L and other channels, and coupling of the probe to the trapped mode, must be derived. They cannot be declared nonzero or resonant from their names. A derivative taken at artificially fixed total energy is also not a derivation of zero mechanical force.
6. **The old positive fixed-value coefficient remains conditional.** Appendix B computes its displayed conjugate functional. Claude explicitly did not derive the physical holder or reservoir law justifying that functional. The computed coefficient is preserved; a native mechanical force does not follow solely from choosing its sign. A mechanical displacement constraint and a reservoir with a conjugate thermodynamic control are not interchangeable without the actual work pairing.
7. **A similar sign does not prove a literal boundary classification.** The old classifier separately admits mixed value/conormal relations. A finite compliant core is not promoted to a specific old boundary class solely by expanding its local spring energy.

Claude read the proposal and electric report fully, but inspected only parts of the interpretation and Python source, and did not inspect the supplied Wolfram source or interface geometry. This assessment checks the specific source discrepancies above; it does not recertify the old calculations. The review's final two-option list is not an exhaustive physical theorem.

## The actual modeling choice

The preceding user-approved direction was an autonomous material response with explicit storage and exchange ports. It did not specify whether a continuing energy supply maintains the particle. The source documents admit an open throat but do not supply its core storage functions, endpoint/return law or their incremental response.

| Hypothesis to test | Physical picture | Main obligation |
|---|---|---|
| Autonomous flow plus stored energy | Trapped waves and core structure coexist with an explicit flow and return through the same medium | Identify the energy and momentum source and return, then derive force and wave response from that same law. Steady flow alone does not prove power input or the correct force sign. |
| Finite stored energy | The core's wave/vortex energy supports a state without a continuing power supply | Establish its support and coupling; any outgoing energy must reduce storage or be balanced by incident energy. The simple displacement-only comparator has the conditional sign problem above. |
| Externally held benchmark | A specified apparatus maintains the chosen throat or profile | Include its work. This can define a supported wave experiment but does not establish an autonomous particle. |

The first is the recommended research hypothesis because it fits the native open-throat direction already chosen, **not** because its force sign is known. The finite-store case remains a useful comparison. The question to the user is which physical hypothesis to investigate, not which desired force sign to impose or whether another script may run.

The smallest subsequent physics deliverable is one explicit energy/momentum and mouth-work specification for the chosen hypothesis. It must name what is held during a slow separation change, identify each internal and external work term, retain both orientations with the same parameters, and connect to the actual source convention. It must not import the positive fixed-value sign, invent a reservoir to repair a result, or choose new parameters after seeing the force. If source ownership remains unresolved, record that limit before a force calculation.

No new review is launched solely to obtain different wording or a CLEAR verdict. No general verifier, numerical solver, leakage worker, centre-drive implementation or mixed numerical recovery is resumed. All previous results and failures remain unchanged and conditional; no leakage factor follows from this assessment.
