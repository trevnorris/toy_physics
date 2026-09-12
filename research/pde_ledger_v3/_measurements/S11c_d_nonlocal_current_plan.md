# S11c-d closed energy current and channel normalization

2026-09-12. Start from checkpoint `18f3236a`, which saved the bulk exceptional
geometry and generalized threshold chains. The PROGRAM_BRIEF, both cleared
authorities, S11b energy-accounting rules, and retained eight-point solver/export
contract govern this construction.

1. Derive the slab energy balance from the inherited energy operand after the
   engine's own tangential/end reduction. Keep the instantaneous, zero-transfer
   virtual constraint distinct from the actual mass evolution. Retain the
   existing conservative boundary current. Compute the unrestricted energy
   variation, chemical functional derivative, actual time-translation boundary
   term, and the correction involving the material mass-rate defect. Emit the
   variation and integration-by-parts residuals before guarding them. Response
   kernels enter only after variation.
2. Construct the outgoing acoustic field ansatz from S11b section 2. Derive
   pressure, velocity, wave equation, energy balance and normal current. Solve
   the supplied face closure with the chemical derivative from step 1 and the
   reduced face geometry. Compare the reconstructed closure contributions with
   the actual reduced mass/thickness rows, retaining independent relaxation
   times. Integrate the bulk current in outward depth where it converges;
   preserve finite-depth operands and explicit exceptional/convergence domains.
   Do not replace a divergent integral by an analytic finite flux without a
   supplied prescription.
3. Evaluate the combined current on every resolved physical mode/subspace.
   Derive the left-field/row maps and the identity involving the actual pencil's
   frequency and normal-momentum derivatives. Compute full normalization and
   current matrices within each degenerate space before selecting a basis.
   Flux normalization is conditional on computed real signed currents and
   admissible sheet/channel domains. Do not replace physical current by a bare
   group velocity, or normalize an unresolved/closed mode as an open channel.
4. Develop against the source-pinned LAB_HELD/RHO4_CONSTANT cache. Use separate
   controls for mass transfer off, bulk loading, gradient-energy boundary work,
   and finite/infinite depth. After the construction is runnable and verified,
   integrate its emissions and run all four cases once with frozen sources.
   Inventory metadata, old-object preservation and domains, publish outputs
   atomically, and update the concise report. An intermediate preflight is a
   runnable checkpoint, not completion of the current/flux TODO.

The finite spectrum/threshold/sheet records are inputs to this work, not global
coverage certificates. Constant-end normal-momentum poles remain distinct from
section 3b's profile-dependent frequency poles. The supplied section 1 premises,
c2 operand/sign debt and existing shear-normalization debt remain explicit.
If a new premise or upstream repair becomes necessary, stop and tell the user
before changing the work. No S10/Lean edits, review legs, comparator, Wolfram,
downstream run or placeholder export. The initial checkpoint commit was
explicitly requested; subsequent work stays in build -> run -> report -> STOP.
