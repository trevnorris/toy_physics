# Order conversion and throat energy bridge

Status: source inventory and bounded preparation following the user's `S_leak` correction. No symbolic restoration, old constructor, numerical calculation or new scientific worker was run. This note identifies what the existing model already supplies before proposing further physics. It does not select a new support law or certify a native force or power.

The main finding is that **the material model already connects order to shear storage and names order-conversion work**. We need to connect those terms to the actual throat and trapped mode. Re-inventing a generic return loop or assuming a pump would skip this existing content.

## Material return is already specified in meaning

The [V3 step plan](../research/pde_ledger_v3/V3_STEP_PLAN.md), S11b banked user postulate and S12, specifies re-ordering of the same bulk medium into the brane state and the dynamical order balance:

    ∂t(χ_B n) + ∇₄·(χ_B n u + J_χ) = n Γ_B,
    Γ_B = Γ_return − Γ_drain.

Total constituent density has its separate conserved balance. Positive Γ_B increases brane order; negative Γ_B de-structures it. The projected ordered-density source contains transverse transport and conversion, with `S_convert = Integral(W n Γ_B dw)`. The original `S_leak` is recovered under its recorded reduction, not relabelled as Γ_B itself.

The earlier [stage-045 preparation](../research/pde_ledger_v2/notes/stage045_nonvariational_block_prep.md) also explicitly records the user's selection of dynamical order conversion and rejection of the frozen-wall mass-sink plus remote-return construction. It is a historical preparation note, not a completed drain solver or current execution directive. The V3 plan keeps Γ_return, Γ_drain, J_χ and return-controller forms open; their physical meaning is not open in the same way as their constitutive functions.

## Existing energy ingredients

The operative amended [stage-006 record](../research/pde_ledger_v2/notes/stages/ledger_stage006_two_phase_chiB_ontology.md), P1–P12 and its corrections, supplies the following classifications:

| Existing term | What it provides | Limit to preserve |
|---|---|---|
| Constituent kinetic energy and U(n) | Material motion and compressional energy | n is a number density; the kinetic term requires m_GNLS. |
| f_B(χ_B) = a_B χ_B²(1−χ_B)² | Postulated two-phase order energy | The two pure-phase minima have equal bare energy at the same density. This term alone does not quantify power released by ongoing re-ordering. |
| κ_B (∇χ_B)²/2 | Interface energy | A recorded single-kink check is not a solved stable finite slab or throat. |
| χ_B f_shear | Order-dependent shear storage | A candidate connection to trapped shear support; the actual mode, geometry and overlap remain necessary. |
| f_throat | Place for the core/wall coupling | Explicitly a deferred placeholder in the inspected material specification. |
| f_mix | Other coupled-sector energy in that historical specification | The old handoff includes Maxwell/gauge content here. It cannot be imported as a derived native h/u_T electric coupling. |

The two-phase action and order field are postulated inputs. Recorded identities and wall checks are earned relative to that closure; they do not derive the medium or its autonomous throat.

The stored shear term is relevant for a simple physical reason: a region's ability to carry shear and its shear energy depend on its order state. This makes the existing term a candidate place to calculate exchange as the interface changes. It does **not** yet show that re-ordering pumps a particular support mode, fixes its amplitude, or chooses a signed mouth boundary condition.

## Order work and interface work are different objects

The corrected material-state convention is

    μ_χ = δF/δχ_B,
    P_order = Integral(μ_χ D_tχ_B d⁴X).

The [handoff erratum](../notes/brane_bulk_handoff.md) and stage-006 P8 explicitly reject the extra factor n in the older `μ_χ n Γ_B` expression. P_order is the order-change contribution in the recorded convention. Advection, changing geometry, material energy transport and boundary work must still be included when forming a complete control-volume balance. It is not automatically input power to a trapped wave.

The stage-006 relaxation equation is labelled an adjunct, not a completed native kinetic law. If a later model includes both a relaxation term and Γ_B, it must state which processes Γ_B already includes and remain consistent with the original order balance and J_χ. One must not count the same conversion a second time by combining equations with different definitions of the total rate.

The existing [S11b shared physics](../research/pde_ledger_v3/directives/S11b_SHARED_PHYSICS.md), §§3–6, instead describes finite-frequency **face** response:

    J_s = Λ_A(ω) A_s + Λ_V(ω) V_s,
    A_s = μ_s − δp_s/ρ_m,
    μ_s = μ_θ/ρ_br⁰.

J_s is outward relative **mass** flux per face area, V_s is outward face velocity, and μ_s is specific chemical energy. They must not be identified by name with the constituent-number source nΓ_B or the order derivative μ_χ. A finite-interface projection, mass conversion and geometric normalization are needed to relate them.

The saved S11b two-port power expression is implemented in [the original source, task_b2d](../research/pde_ledger_v3/scripts/S11b_interface_coupling_law_sympy_audit.py). It separates mechanical traction power and chemical mass-flux power, and rewrites them into bulk acoustic power, affinity work and the reciprocal traction contribution. The source was read, not called. This provides an existing accounting structure to reuse with its real-frequency, linear-face and sign conventions; it is not a new nonlinear core balance.

The [S11b result record](../research/pde_ledger_v3/steps/S11b_interface_coupling_law.md) treats passivity as a conditional region, not a ban on every active response. A response outside that region requires a named reservoir and stated budget. Its normal background drain remains an **uncarried** convective scope limit in rest-bulk wave equations. Consequently, the previous strict-rest-bulk S11c results cannot be promoted into a calculation of power supplied by a finite steady return flow.

## The missing connection is narrower than a new support ontology

The source inventory yields four specific joins to investigate:

1. **Conversion to face transfer:** how the finite-interface χ_B balance, relative transport and conversion reduce to the actual outward face mass fluxes. In-place re-ordering and material crossing a moving face must remain distinguishable.
2. **Order work to core work:** which parts of the existing energy and stress act on throat shape, normal mouth displacement and a support mode. The f_throat placeholder cannot be replaced by a fitted stiffness or a reservoir just to obtain the desired electric sign.
3. **The support mode itself:** its native profile and energy normalization, and the work it receives or gives back under interface motion. Calling energy trapped does not establish the sign or magnitude of its support stress.
4. **Background versus perturbation:** how the committed slow re-ordering background and the wave response share one set of laws. A wave-period approximation does not set background conversion or incremental work to zero.

These are preparation questions, not claims that the desired answers are absent everywhere in the repository. The inspected material record explicitly defers its core functions; the native body-dynamics [report](../software/em_charge_attribute/reports/u1_body_dynamics.md) also carries core and return inputs as open. Further source tracing should follow these entries before adding new constitutive content. Older G0 or Maxwell-based calculations retain their conditional status and cannot override the later dynamical-order choice.

The smallest next deliverable is a source-bound specification of these joins and their energy/momentum partition, with any still-missing function named explicitly. It should answer where the work comes from and where it goes before requesting a numerical throat calculation. A new physical law, if actually required, receives substantive Claude-only assessment. No existing scientific calculation needs to be rerun for this inventory, and no new verifier, integration bank, leakage worker or force-sign fit is warranted by it.
