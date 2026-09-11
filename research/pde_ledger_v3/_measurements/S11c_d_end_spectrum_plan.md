# S11c-d end-spectrum construction plan

2026-09-10. Start from preservation checkpoint `c5af3181`. The user requested
a plan and implementation of the next construction. This remains the staged
builder lane: build, run, report, stop. No commit or review is part of this run.

Authority: `directives/S11c_d_sympy_build_PROGRAM_BRIEF.md`, the cleared
`directives/S11c_d_SHARED_PHYSICS.md` sections 2 and 3a, and
`directives/S11c_d_sympy_build_directive.md`. Retain the eight-point approved
solver/export contract in `S11c_d_sympy_builder_report.md`.

## Objective and limits

Extend the existing engine beyond regular reference-mode jets by computing
end-spectrum coverage from the full reduced pencils at explicit parameter
inputs. Account for candidate roots, multiplicities, denominator exclusions,
coordinate-chart failures, and left/right nullspaces. Preserve candidates on
both radical sheets and the existing fixed-frequency continuation record.
An unresolved domain is not an absent mode. Mode direction by the full energy
current, generic frequency/sheet continuation, and flux normalization remain
their separately named constructions.

Use the existing independent thickness-step and modulus-bump input in
`S11c_d_channel_preflight_input.json`, including its declared L/T/M frame and
physical eta/sigma homotopy, as the first physical instance. Algebraic PIT
records remain separate. A complete spectrum at one bound input does not
establish coverage of all profiles or parameter domains. Do not remove the
full-spectrum TODO until the implemented scope supports that change.

## Sequence

1. **Construction preflight, one case.** Read the pinned reference/left/right
   symbol caches from the completed inverse-Fourier run. Compare the physical
   five-field symbol with the potential-coordinate quotient used for regular
   jets. Compute determinant denominators, elimination support and coordinate
   singularities; measure their effect on candidate roots. Verify the cache,
   producer transcript and source bindings before use. This identifies the
   regular coverage domain and any roots that need separate treatment.
2. **Root coverage and modes.** Extend the engine with a carrier-first end
   spectrum census that consumes those computed reduced symbols. Compute the
   finite characteristic roots with multiplicities and independent root-count
   or isolation evidence, evaluate the original pencil at the candidates, and
   retain exceptional loci explicitly. Compute left/right mode spaces and
   classifier applicability; never identify a coordinate artifact as a
   physical channel. Print operands and residuals, including failed cases.
3. **Focused verification.** Develop on LAB_HELD/RHO4_CONSTANT. Compare with
   the existing regular-jet roots on their common domain, refine numerical
   precision, and exercise an exceptional coordinate or root domain through
   changes to the bound input. Emit all new objects with restored dimensions,
   grades, source/input hashes and compact fingerprints for heavy tensors.
4. **Integrated checkpoint.** After the focused path is stable, run the full
   four-case engine once, inventory the new records and check preservation of
   the reduction, input wiring and earlier computed objects. Publish completed
   output by atomic replacement under `scripts/out/`, preserving annex payloads.
   Update the concise builder report and construction/run records. Leave any
   unfinished construction explicitly in the live TODO and stop.

## Stop condition

If the preflight or construction exposes an incorrect earlier operand or
requires a different physical premise or a change to the approved approach,
preserve the evidence and stop to explain it to the user before repairing it.
Routine implementation and explicit treatment of the spectrum's already-open
domain limits proceed under this plan. Section 1 supplied inputs and the
existing upstream cross-engine/shear debts remain explicit.

## Completed checkpoint, 2026-09-11

Steps 1–4 completed for regular finite algebraic coverage at bound inputs.
Five focused checks and one fresh four-case run completed, with 24 polynomial
root-isolation certificates and 432 native candidate records. All 12 previous
pencil symbols, all 528 legacy candidate records and the inverse-Fourier
checks are preserved. The six successful transcripts are under `scripts/out/`;
the [construction report](/var/projects/toy_physics/research/pde_ledger_v3/_measurements/S11c_d_end_spectrum_report.md)
and run record document their hashes and scope. No upstream repair was needed.

The full-spectrum TODO remains open: generic sheet continuation, cut-bank and
continuum treatment, exceptional threshold/denominator domains, mixed
degeneracies and defective-root generalized modes still precede a complete
physical channel space. Current/flux normalization and scattering remain later
constructions. No export, review or commit was performed in this checkpoint.
