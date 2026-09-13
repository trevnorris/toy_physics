# S11c-d combined modal current continuation

## Completed focused checkpoint

The user approved the acoustic/face-lift split after the previous pause. It is
implemented and committed as `5b007ddb`: 802 literal zero scalar residuals in 34
families, with complete metadata. The local acoustic balance is divided on
both wave equations before the closed face lifts. Full five-field row joins,
bulk/port matrix reconstructions, the slab balance and finite-depth composed
balances are checked, including the separate equal-depth branch and the
original denominator domain.

Regular-sheet normal/frequency derivative operands, local wave tangency,
source product-rule reconstructions and the equal-depth integral derivative
are also computed and checked. The original report and inventory specify
which checks are symbolic in material parameters and which use the declared
reference material binding. This does not establish global spectral, sheet,
threshold or defective-mode coverage. No new upstream repair is indicated by
this checkpoint.

The new `ModalCurrentSubspaces` continuation covers all 18 isolated candidates
in the LAB_HELD/RHO4_CONSTANT reference packet, with both momentum lifts and
the complete 14 scalar and four two-dimensional nullspaces. All 18 frequency
pairings have full computed rank. It contracts the finite-depth current and
energy and reconstructs them from the actual pencil derivatives and row-power
maps, preserving their source, interface, upper-boundary and projection-defect
terms. Eight candidates have disk-certified bulk decay. Two physical-sheet,
exactly real-normal two-dimensional spaces admit signed physical right-current
normalization: four basis directions. See the
[subspace report](S11c_d_modal_subspace_report.md) and its inventory for domains
and numerical tolerances. This is a focused reference result.

The canonical row-dual frequency pairing and the physical right-field current
remain separate typed objects. Their computed power-map bridge and defect
are retained. The canonical adjoint-field current export still needs that
explicit map; a bare row-dual vector cannot be inserted into a field-current
slot without it. No upstream discrepancy or additional premise has been
identified by this calculation.

## Completed reference implementation sequence

1. Reuse the repaired reference's isolated-root records with their source and
   input pins. Reconstruct every full right/left subspace for each isolated
   root and both momentum lifts, preserving basis-rank and projector checks.
   The first focused packet has 18 candidates; an aggregate witness is not a
   substitute for these records. Preserve all exceptional/domain outcomes.
2. Contract the actual two-frequency physical current and the full nonlinear
   frequency derivative on those subspaces. Finish their relation using the
   computed field-to-row power maps, retaining interface exchange, the upper
   boundary and non-Hermitian terms. Handle conjugate physical legs and the
   equal-depth limit explicitly; use the computed integral derivative when
   differentiating before coincidence. Do not substitute bare group velocity
   or assume the physical-current and left-eigenvector pairings coincide.
3. Compute complete subspace current matrices before selecting a basis. Assign
   signed flux normalization only where the actual sheet, propagation,
   reality, convergence and rank conditions allow it. Record evanescent,
   zero-current, threshold, defective and nonconvergent cases individually.
   The regular-sheet radical derivative does not apply at its zero denominator.

## Next implementation sequence

1. Exercise the planned mass-transfer, bulk-load and boundary-work source
   controls in the focused case, mutating the action/imported operands and
   reconstructing the affected current and pencil. Preserve each control's
   changed spectrum and domain instead of reusing unaltered mode witnesses.
2. Finish the canonical adjoint-row/physical-field current join using the
   computed power maps, retaining projection defects and rank conditions.
   Extend the verified full-subspace calculation to both ends and all cases;
   preserve per-root exceptional outcomes and degenerate current matrices.
3. Continue the variable-profile current and matching work, controls and
   bookkeeping under the retained solver contract. Integrate the completed
   new emissions and run all four cases once when that stage is ready, then
   write the complete export. Do not replace unfinished objects with
   placeholders or promote a focused checkpoint to engine completion.

Preserve the native import wiring and Fourier reduction. Use the actual
repaired closed pencil and independently retained response kernels; do not
vary a response or tune inputs to manufacture channels. Heavy objects use
carrier/PIT fingerprints with dimensions and multigrades. Work with one heavy
CAS process and one development case, source-check caches, and publish
transcripts atomically under `scripts/out/`.

No S10/Lean edits, authority changes, review legs, comparator, Wolfram or
downstream runs. All ten broad TODOs remain. Constant-end momentum poles are
separate from §3b's profile-dependent frequency bound poles, whose existence
requires a targeted spectral-locus calculation. If a new upstream discrepancy
or additional physical premise requires changing this approach, stop and
report it before making that change. No new change of approach is requested
at this checkpoint.
