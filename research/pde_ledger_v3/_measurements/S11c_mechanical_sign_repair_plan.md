# S11c mechanical-load repair scope

The user authorized this repair on 2026-09-12 after the bounded sign audit and
checkpoint `03491aca`. This is the execution plan, not a record of completion.
The source/export/transcript baseline is pinned in
`_measurements/S11c_mechanical_repair_baseline.json`.

1. Preserve a source/export/transcript baseline and the retained S11c-d solver
   contract. Repair the placement of prescribed external virtual work in
   S11c-b's mechanical balance. Keep the extracted physical generalized force
   distinct from its contribution to the stored left-hand-side row. Derive the
   relative normalization from the supplied kinetic and constrained stored-energy
   variation; do not select it by a dispersion root or a desired stability sign.
   Preserve the separate mass-evolution row, independent relaxation kernels,
   material virtual constraint, and the completed inertia repair.

2. Carry that computed mechanical load through local/expanded rows and the
   `FACE_VIRTUAL_WORK` provenance used by weak restriction. Rebuild the operator
   and complete coupling blocks from those operands. Expected affected b roots
   are `slab_operator`, `slab_operator_term_origins`, `coupling_kernel`, and
   `coupling_kernel_term_origins`; determine the actual changed-root set by
   serialization and symbolic comparisons. Preserve the physical-force source
   operand's meaning. A narrow sign repair should not change the energy basis,
   chemical functional derivative, geometry, or background-order support;
   verify these before publishing.

3. Strengthen c2's traction/slab power control. Compute kinetic-and-stored power
   independently from the supplied energy variation and material virtual
   constraint. Do not define it by subtracting the assembled face force from
   assembled slab power: that subtraction cancels the row it is supposed to
   validate. Compare the independently computed balance with the closed
   mechanical rows, holding prescribed responses fixed under virtual variation.
   Exercise separate source-level changes to traction and to its mechanical-row
   routing; report both operands and residuals. Keep this upstream diagnostic
   separate from S11c-d's own reduced energy current.

4. Develop on LAB_HELD/RHO4_CONSTANT, then verify all four native cases with
   nonuniform profiles. Check both faces before summing, all four mechanical
   components, energy/inertia orientation, and the full independent memory
   channels. Include impermeable and zero-face-source controls and the uniform
   S11b face-work comparison. Compute unchanged mass-row, kinetic, constitutive,
   and non-face mechanical differences against the baseline. Preserve reciprocal
   denominators and their domains; zero rational residuals do not extend the
   expressions to an excluded denominator or certify spectral strata.

5. Regenerate serially: b, c1, c2, then the focused d reduction/current checks.
   c1 consumes b's `mu_theta_operator` and inherited geometry/closure inputs,
   not `slab_operator` or either coupling/provenance root. Thus a narrow face-load
   repair predicts a provenance refresh with unchanged c1 physical values;
   verify the direct-input census and all 44 c1 values after the rebuild.
   c2 substitutes pressure and pressure jets into b's rows, so its closed slab
   and coupling outputs require physical regeneration. Preserve existing
   deferred-heavy-control boundaries explicitly rather than reporting a full
   upstream clearance.

6. Rebuild d's reduced rows from the new exports. Repeat the current preflight
   first: compare the mass row, full five-coefficient mechanical load, and
   independent stiffness/inertia orientation. Recompute every spectrum, sheet,
   threshold, and current record that depends on changed rows. Historical root
   counts and availability classifications must not become new expected answers.
   Follow the existing per-root/subspace and exceptional-locus coverage plan.
   Only after the focused checks finish should the final full four-case d run
   publish a new transcript. Resume current/flux normalization and the remaining
   scattering program afterward; do not create an incomplete export.

No current evidence calls for a physical change to S11b or S11c-a. Verify their
source pins and preserve their outputs. This is not cross-engine adjudication:
the separate shear normalization and the carried c1/c2 representation, response,
and sign debts remain open. If the repair exposes a new upstream discrepancy or
requires changing an authority, stop and identify that change of scope.

Use one memory-heavy CAS job at a time, source-pinned temporary outputs, measured
resources, exact metadata/residual inventories, and atomic destination-local
publication. Keep `.out` files under `scripts/out/` for DataLad/git-annex at the
next requested checkpoint; use Git for code and reports. No review leg,
comparator, Wolfram run, S10/Lean edit, downstream task, or automatic commit.
