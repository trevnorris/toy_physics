# S11c-d adjoint-row / physical-field current map

Checkpoint `c7f2d879` saves the exact normal-reality repair and completed source
controls. Continue the approved current-map step on the same source-pinned
LAB_HELD/RHO4_CONSTANT reference packet; do not regenerate its root isolation
or replace any earlier field/current tensors.

1. Use the plus-row power map B obtained from the actual source-work ansatz.
   On each entire root subspace, compute its singular values and rank. Where
   B is invertible, solve B^dagger A = L_omega for the physical test-field
   representation of the frequency-normalized row-dual modes. Retain rank
   failures explicitly; do not substitute a pseudoinverse as an invertible map.
2. Form the power-weighted row representation Q = B P. Differentiate this
   computed product along the existing radical curve, including B_omega P
   and B_k P. Check the right and adjoint kernels, the row-dual/field map, and
   the full nonlinear frequency and normal pairings on every basis vector.
   Q is an invertible row representation on this domain, not a new action or
   a replacement for the physical pencil. All matrices use the already
   declared epsilon-squared coefficient convention.
3. Contract the derived finite-depth energy/current on A and R. Reconstruct
   those bilinears from the full differentiated balance, retaining the other
   source-work, interface, depth-boundary and phase-derivative terms. Do not
   identify a mixed adjoint-field bilinear with physical right-field flux by
   deleting a phase, bridge, defect or energy weight.
4. Lift the existing power-covector projection defect back to fields. Verify
   R = A G_omega^dagger + F and reconstruct R^dagger J R and R^dagger E R
   from the mixed forms and the defect term. For disk-certified bulk decay,
   retain the infinite-depth forms. Carry the already computed signed-current
   maps through the same identity; preserve all original normalization gates.
5. Check full matrix covariance under a nonunitary change of basis within
   each degenerate space: R -> R T and L_omega -> L_omega T^{-dagger}.
   This is a coordinate regression, not an independent physical derivation.
6. Emit bounded fingerprints, literal residuals, grade/unit metadata and a
   final source-line index. Source/input/cache hashes guard reuse. Publish
   atomically under scripts/out, preserving the committed reference transcript.
   Record the scoped result and remaining both-end/profile-matching work.

The initial pinned numerical preflight finds rank five for B at all 18
reference candidates. Its direct inverse and current reconstruction diagnostics
are small; it supplies no generic invertibility or exceptional-domain theorem.
The algebraic field representation has its own units and phase convention.
No upstream physics premise or authority change is indicated by this preflight.
If the construction requires one, stop and report it. Global/exceptional
coverage, profile-frequency poles, matching/scattering and the complete export
remain outside this focused checkpoint. No S10/Lean edits or review/downstream
run; one heavy CAS process at a time. The new work after c7f2d879 is uncommitted.

## Execution checkpoint

The six steps above are implemented for the pinned reference packet. All 18
field maps and 22 basis directions are retained; the 50 symbolic product-rule
residual scalars are zero, and the largest numerical residual norm is
3.679e-13. The validated transcript is published under scripts/out; see the
[report](S11c_d_adjoint_current_report.md) and
[inventory](S11c_d_adjoint_current_checkpoint.json). Next is the same current
construction at both ends, followed by one-case profile matching. No additional
physical premise or upstream repair was required.
